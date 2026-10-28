using System.Diagnostics;
using Sdcb.Arithmetic.Core;

namespace Sdcb.Arithmetic;

public readonly partial struct BigInt
{
    // ---------------------------------------------------------------
    // Bitwise operators (infinite two's-complement semantics, like mpz)
    // ---------------------------------------------------------------

    public static BigInt operator &(BigInt a, BigInt b) => BitwiseOp(a, b, 0);
    public static BigInt operator |(BigInt a, BigInt b) => BitwiseOp(a, b, 1);
    public static BigInt operator ^(BigInt a, BigInt b) => BitwiseOp(a, b, 2);

    /// <summary>One's complement: ~x == -x - 1.</summary>
    public static BigInt operator ~(BigInt value) => -value - One;

    private static BigInt BitwiseOp(BigInt a, BigInt b, int op)
    {
        if (a._size >= 0 && b._size >= 0)
        {
            // Fast path: plain magnitude logic.
            int an = a.LimbCount, bn = b.LimbCount;
            if (an < bn)
            {
                (a, b) = (b, a);
                (an, bn) = (bn, an);
            }
            ulong[] d = new ulong[an];
            switch (op)
            {
                case 0:
                    Nat.And(d, a.Limbs, b.Limbs, bn);
                    break;
                case 1:
                    Nat.Or(d, a.Limbs, b.Limbs, bn);
                    a.Limbs.Slice(bn).CopyTo(d.AsSpan(bn));
                    break;
                default:
                    Nat.Xor(d, a.Limbs, b.Limbs, bn);
                    a.Limbs.Slice(bn).CopyTo(d.AsSpan(bn));
                    break;
            }
            return FromMagnitude(d, an, 1);
        }

        // General path via two's complement.
        int len = Math.Max(a.LimbCount, b.LimbCount) + 1;
        ulong[] ta = a.ToTwosComplement(len);
        ulong[] tb = b.ToTwosComplement(len);
        switch (op)
        {
            case 0:
                Nat.And(ta, ta, tb, len);
                break;
            case 1:
                Nat.Or(ta, ta, tb, len);
                break;
            default:
                Nat.Xor(ta, ta, tb, len);
                break;
        }
        return FromTwosComplement(ta);
    }

    private ulong[] ToTwosComplement(int len)
    {
        Debug.Assert(len > LimbCount);
        ulong[] t = new ulong[len];
        Limbs.CopyTo(t);
        if (_size < 0)
        {
            NegateTwosComplementInPlace(t);
        }
        return t;
    }

    private static void NegateTwosComplementInPlace(Span<ulong> t)
    {
        ulong carry = 1;
        for (int i = 0; i < t.Length; i++)
        {
            ulong s = ~t[i] + carry;
            t[i] = s;
            carry = carry != 0 && s == 0 ? 1UL : 0UL;
        }
    }

    private static BigInt FromTwosComplement(ulong[] t)
    {
        int len = t.Length;
        if ((t[len - 1] >> 63) == 0)
        {
            return FromMagnitude(t, len, 1);
        }
        NegateTwosComplementInPlace(t);
        return FromMagnitude(t, len, -1);
    }

    // ---------------------------------------------------------------
    // Individual bit access (mpz_tstbit / setbit / clrbit / combit)
    // ---------------------------------------------------------------

    /// <summary>Tests bit <paramref name="index"/> of the infinite two's-complement representation.</summary>
    public bool TestBit(long index)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(index);
        long limbIndex = index >> 6;
        int bit = (int)(index & 63);
        ReadOnlySpan<ulong> m = Limbs;
        int n = m.Length;
        if (limbIndex >= n)
        {
            return _size < 0; // zero / sign extension
        }
        int limb = (int)limbIndex;
        if (_size >= 0)
        {
            return ((m[limb] >> bit) & 1) != 0;
        }
        // Bit of ~(|v| - 1): borrow from subtracting 1 reaches limb only if all lower limbs are zero.
        bool lowerZero = true;
        for (int i = 0; i < limb; i++)
        {
            if (m[i] != 0)
            {
                lowerZero = false;
                break;
            }
        }
        ulong t = lowerZero ? m[limb] - 1 : m[limb];
        return ((~t >> bit) & 1) != 0;
    }

    /// <summary>Returns a copy with bit <paramref name="index"/> set (mpz_setbit).</summary>
    public BigInt SetBit(long index) => WithBit(index, true);

    /// <summary>Returns a copy with bit <paramref name="index"/> cleared (mpz_clrbit).</summary>
    public BigInt ClearBit(long index) => WithBit(index, false);

    /// <summary>Returns a copy with bit <paramref name="index"/> flipped (mpz_combit).</summary>
    public BigInt FlipBit(long index) => WithBit(index, !TestBit(index));

    private BigInt WithBit(long index, bool value)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(index);
        if (TestBit(index) == value)
        {
            return this;
        }
        // Changing a bit beyond the current length needs limb + 2 limbs; overflow means the result cannot be represented.
        int limb = checked((int)(index >> 6));
        int bit = (int)(index & 63);
        if (_size >= 0)
        {
            int len = Math.Max(LimbCount, checked(limb + 1));
            ulong[] d = new ulong[len];
            Limbs.CopyTo(d);
            d[limb] ^= 1UL << bit;
            return FromMagnitude(d, len, 1);
        }
        // Negative: operate on two's complement.
        int tlen = checked(Math.Max(LimbCount, limb + 1) + 1);
        ulong[] t = ToTwosComplement(tlen);
        t[limb] ^= 1UL << bit;
        return FromTwosComplement(t);
    }

    // ---------------------------------------------------------------
    // Population count / scanning
    // ---------------------------------------------------------------

    /// <summary>Number of set bits; defined for non-negative values only (mpz_popcount).</summary>
    public long PopCount()
    {
        if (_size < 0)
        {
            throw new InvalidOperationException("PopCount is defined for non-negative values only.");
        }
        return Nat.PopCount(Limbs, LimbCount);
    }

    /// <summary>Hamming distance between two values of the same sign (mpz_hamdist).</summary>
    public static long HammingDistance(BigInt a, BigInt b)
    {
        if (a._size >= 0 && b._size >= 0)
        {
            int an = a.LimbCount, bn = b.LimbCount;
            if (an < bn)
            {
                (a, b) = (b, a);
                (an, bn) = (bn, an);
            }
            long count = Nat.HammingDistance(a.Limbs, b.Limbs, bn);
            count += Nat.PopCount(a.Limbs.Slice(bn), an - bn);
            return count;
        }
        if (a._size < 0 && b._size < 0)
        {
            int len = Math.Max(a.LimbCount, b.LimbCount) + 1;
            ulong[] ta = a.ToTwosComplement(len);
            ulong[] tb = b.ToTwosComplement(len);
            return Nat.HammingDistance(ta, tb, len);
        }
        throw new InvalidOperationException("HammingDistance requires operands of the same sign.");
    }

    /// <summary>Index of the first 1 bit at or above <paramref name="startBit"/>; -1 if none (mpz_scan1 for non-negative values).</summary>
    public long IndexOfSetBit(long startBit)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(startBit);
        if (_size >= 0)
        {
            return Nat.IndexOfSetBit(Limbs, LimbCount, startBit);
        }
        // Negative numbers have infinitely many 1 bits; scan two's complement.
        long idx = startBit;
        while (!TestBit(idx))
        {
            idx++;
        }
        return idx;
    }

    /// <summary>Index of the first 0 bit at or above <paramref name="startBit"/> (mpz_scan0); -1 if none (negative values end in infinite 1s).</summary>
    public long IndexOfClearBit(long startBit)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(startBit);
        if (_size >= 0)
        {
            return Nat.IndexOfClearBit(Limbs, LimbCount, startBit);
        }
        long limit = GetBitLength() + 1;
        for (long idx = startBit; idx <= limit; idx++)
        {
            if (!TestBit(idx))
            {
                return idx;
            }
        }
        return -1;
    }

    /// <summary>Number of trailing zero bits (mpz_scan1 from 0); throws for zero.</summary>
    public long TrailingZeroCount()
    {
        if (_size == 0)
        {
            throw new InvalidOperationException("TrailingZeroCount is undefined for zero.");
        }
        return Nat.TrailingZeroCount(Limbs, LimbCount);
    }
}
