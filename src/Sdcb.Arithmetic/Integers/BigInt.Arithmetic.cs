using System.Diagnostics;
using Sdcb.Arithmetic.Core;

namespace Sdcb.Arithmetic;

public readonly partial struct BigInt
{
    /// <summary>Allocates a limb buffer without zeroing; every limb must be written before use.</summary>
    [System.Runtime.CompilerServices.MethodImpl(System.Runtime.CompilerServices.MethodImplOptions.AggressiveInlining)]
    internal static ulong[] NewLimbs(int n) => GC.AllocateUninitializedArray<ulong>(n);

    // ---------------------------------------------------------------
    // Negation / absolute value
    // ---------------------------------------------------------------

    public static BigInt operator -(BigInt value) => new(value._d, -value._size);
    public static BigInt operator +(BigInt value) => value;

    public BigInt Negate() => new(_d, -_size);
    public BigInt Abs() => new(_d, Math.Abs(_size));

    public static BigInt Abs(BigInt value) => value.Abs();

    // ---------------------------------------------------------------
    // Addition / subtraction
    // ---------------------------------------------------------------

    public static BigInt operator +(BigInt a, BigInt b) => AddSigned(a, b, negateB: false);
    public static BigInt operator -(BigInt a, BigInt b) => AddSigned(a, b, negateB: true);

    public static BigInt operator ++(BigInt value) => value + One;
    public static BigInt operator --(BigInt value) => value - One;

    private static BigInt AddSigned(BigInt a, BigInt b, bool negateB)
    {
        int bSize = negateB ? -b._size : b._size;
        if (a._size == 0)
        {
            return new BigInt(b._d, bSize);
        }
        if (bSize == 0)
        {
            return a;
        }

        int an = a.LimbCount, bn = Math.Abs(bSize);
        ReadOnlySpan<ulong> am = a.Limbs, bm = b.Limbs;
        bool sameSign = (a._size > 0) == (bSize > 0);
        if (sameSign)
        {
            // Magnitude addition.
            if (an < bn)
            {
                ReadOnlySpan<ulong> tm = am; am = bm; bm = tm;
                (an, bn) = (bn, an);
            }
            ulong[] d = NewLimbs(an + 1);
            ulong carry = Nat.Add(d, am, an, bm, bn);
            d[an] = carry;
            return FromMagnitude(d, an + 1, a._size > 0 ? 1 : -1);
        }
        else
        {
            // Magnitude subtraction: bigger - smaller.
            int cmp = Nat.CompareNormalized(am, an, bm, bn);
            if (cmp == 0)
            {
                return default;
            }
            int sign;
            if (cmp < 0)
            {
                ReadOnlySpan<ulong> tm = am; am = bm; bm = tm;
                (an, bn) = (bn, an);
                sign = bSize > 0 ? 1 : -1;
            }
            else
            {
                sign = a._size > 0 ? 1 : -1;
            }
            ulong[] d = NewLimbs(an);
            ulong borrow = Nat.Sub(d, am, an, bm, bn);
            Debug.Assert(borrow == 0);
            return FromMagnitude(d, an, sign);
        }
    }

    // ---------------------------------------------------------------
    // Multiplication
    // ---------------------------------------------------------------

    public static BigInt operator *(BigInt a, BigInt b)
    {
        if (a._size == 0 || b._size == 0)
        {
            return default;
        }
        int an = a.LimbCount, bn = b.LimbCount;
        int sign = (a._size ^ b._size) < 0 ? -1 : 1;
        ulong[] d = NewLimbs(an + bn);
        if (a._d == b._d && an == bn)
        {
            Nat.Square(d, a.Limbs, an);
        }
        else if (an >= bn)
        {
            Nat.Mul(d, a.Limbs, an, b.Limbs, bn);
        }
        else
        {
            Nat.Mul(d, b.Limbs, bn, a.Limbs, an);
        }
        return FromMagnitude(d, an + bn, sign);
    }

    /// <summary>Returns value squared (uses the dedicated squaring path).</summary>
    public BigInt Square()
    {
        if (_size == 0)
        {
            return default;
        }
        int n = LimbCount;
        ulong[] d = NewLimbs(2 * n);
        Nat.Square(d, Limbs, n);
        return FromMagnitude(d, 2 * n, 1);
    }

    // ---------------------------------------------------------------
    // Division (truncating operators, plus floor/ceiling variants)
    // ---------------------------------------------------------------

    /// <summary>Truncating division (rounds toward zero), like C# integers and mpz_tdiv_q.</summary>
    public static BigInt operator /(BigInt a, BigInt b)
    {
        DivRem(a, b, out BigInt q, out _);
        return q;
    }

    /// <summary>Truncating remainder; result has the sign of the dividend (mpz_tdiv_r).</summary>
    public static BigInt operator %(BigInt a, BigInt b)
    {
        DivRem(a, b, out _, out BigInt r);
        return r;
    }

    /// <summary>Truncating division with remainder: a = q*b + r, |r| &lt; |b|, sign(r) == sign(a).</summary>
    public static void DivRem(BigInt a, BigInt b, out BigInt quotient, out BigInt remainder)
    {
        if (b._size == 0)
        {
            throw new DivideByZeroException();
        }
        int an = a.LimbCount, bn = b.LimbCount;
        if (Nat.CompareNormalized(a.Limbs, an, b.Limbs, bn) < 0)
        {
            quotient = default;
            remainder = a;
            return;
        }

        int qSign = (a._size ^ b._size) < 0 ? -1 : 1;
        int rSign = a._size < 0 ? -1 : 1;

        if (bn == 1)
        {
            ulong[] q = NewLimbs(an);
            ulong rem = Nat.DivRemScalar(q, a.Limbs, an, b._d![0]);
            quotient = FromMagnitude(q, an, qSign);
            remainder = rem == 0 ? default : new BigInt([rem], rSign);
            return;
        }

        ulong[] qd = NewLimbs(an - bn + 1);
        ulong[] rd = NewLimbs(bn);
        Nat.DivRem(qd, rd, a.Limbs, an, b.Limbs, bn);
        quotient = FromMagnitude(qd, an - bn + 1, qSign);
        remainder = FromMagnitude(rd, bn, rSign);
    }

    /// <summary>Division with the quotient rounded as requested (mpz_{t,f,c}div_q).</summary>
    public static BigInt Divide(BigInt a, BigInt b, IntegerRounding rounding)
    {
        DivRem(a, b, out BigInt q, out _, rounding);
        return q;
    }

    /// <summary>Ieee754Remainder of <see cref="Divide(BigInt, BigInt, IntegerRounding)"/> (mpz_{t,f,c}div_r).</summary>
    public static BigInt Remainder(BigInt a, BigInt b, IntegerRounding rounding)
    {
        DivRem(a, b, out _, out BigInt r, rounding);
        return r;
    }

    /// <summary>Division with remainder, a = q*b + r, with the quotient rounded as requested (mpz_{t,f,c}div_qr).</summary>
    public static void DivRem(BigInt a, BigInt b, out BigInt quotient, out BigInt remainder, IntegerRounding rounding)
    {
        switch (rounding)
        {
            case IntegerRounding.Truncate: DivRem(a, b, out quotient, out remainder); break;
            case IntegerRounding.Floor: FloorDivRem(a, b, out quotient, out remainder); break;
            case IntegerRounding.Ceiling: CeilingDivRem(a, b, out quotient, out remainder); break;
            default: throw new ArgumentOutOfRangeException(nameof(rounding));
        }
    }

    /// <summary>Floor division (rounds toward negative infinity), like mpz_fdiv_q.</summary>
    internal static BigInt FloorDivide(BigInt a, BigInt b)
    {
        FloorDivRem(a, b, out BigInt q, out _);
        return q;
    }

    /// <summary>Non-negative-toward-divisor remainder: result has the sign of the divisor (mpz_fdiv_r; Python %).</summary>
    public static BigInt Mod(BigInt a, BigInt b)
    {
        FloorDivRem(a, b, out _, out BigInt r);
        return r;
    }

    /// <summary>Floor division with remainder (mpz_fdiv_qr): a = q*b + r with 0 &lt;= |r| &lt; |b| and sign(r) == sign(b).</summary>
    internal static void FloorDivRem(BigInt a, BigInt b, out BigInt quotient, out BigInt remainder)
    {
        DivRem(a, b, out quotient, out remainder);
        if (!remainder.IsZero && (remainder._size < 0) != (b._size < 0))
        {
            quotient -= One;
            remainder += b;
        }
    }

    /// <summary>Ceiling division with remainder (mpz_cdiv_qr): a = q*b + r with sign(r) == -sign(b).</summary>
    internal static void CeilingDivRem(BigInt a, BigInt b, out BigInt quotient, out BigInt remainder)
    {
        DivRem(a, b, out quotient, out remainder);
        if (!remainder.IsZero && (remainder._size < 0) == (b._size < 0))
        {
            quotient += One;
            remainder -= b;
        }
    }

    /// <summary>Ceiling division (rounds toward positive infinity), like mpz_cdiv_q.</summary>
    internal static BigInt CeilingDivide(BigInt a, BigInt b)
    {
        CeilingDivRem(a, b, out BigInt q, out _);
        return q;
    }

    /// <summary>Exact division (mpz_divexact): b must divide a; the result is unspecified otherwise. About twice as fast as general division.</summary>
    public static BigInt DivExact(BigInt a, BigInt b)
    {
        if (b._size == 0)
        {
            throw new DivideByZeroException();
        }
        int an = a.LimbCount, bn = b.LimbCount;
        if (an < bn)
        {
            Debug.Assert(an == 0, "DivExact called with a non-exact quotient");
            return default;
        }
        if (bn == 1)
        {
            return DivRem(a, b._d![0], out _) is var q1 && b._size < 0 ? -q1 : q1;
        }
        int qn = an - bn + 1;
        ulong[] q = new ulong[qn];
        Nat.DivExact(q, a.Limbs, an, b.Limbs, bn);
        return FromMagnitude(q, qn, (a._size ^ b._size) < 0 ? -1 : 1);
    }

    /// <summary>Returns true and the quotient if b divides a (mpz_divisible_p + division).</summary>
    public static bool IsDivisible(BigInt a, BigInt b)
    {
        if (b._size == 0)
        {
            return a._size == 0;
        }
        DivRem(a, b, out _, out BigInt r);
        return r.IsZero;
    }

    /// <summary>Truncating division by a ulong, returning the remainder (mpz_tdiv_qr_ui-style fast path).</summary>
    public static BigInt DivRem(BigInt a, ulong b, out ulong remainder)
    {
        if (b == 0)
        {
            throw new DivideByZeroException();
        }
        if (a._size == 0)
        {
            remainder = 0;
            return default;
        }
        int an = a.LimbCount;
        ulong[] q = NewLimbs(an);
        remainder = Nat.DivRemScalar(q, a.Limbs, an, b);
        return FromMagnitude(q, an, a._size < 0 ? -1 : 1);
    }

    /// <summary>Ieee754Remainder of |a| divided by d (always non-negative).</summary>
    public static ulong RemainderAbs(BigInt a, ulong d)
    {
        if (d == 0)
        {
            throw new DivideByZeroException();
        }
        if (a._size == 0)
        {
            return 0;
        }
        return Nat.RemScalar(a.Limbs, a.LimbCount, d);
    }

    // ---------------------------------------------------------------
    // Shifts
    // ---------------------------------------------------------------

    /// <summary>Multiplies by 2^shift (mpz_mul_2exp).</summary>
    public static BigInt operator <<(BigInt value, int shift)
    {
        if (shift < 0)
        {
            return value >> (-shift);
        }
        return value.ShiftLeft(shift);
    }

    /// <summary>Arithmetic right shift: floor division by 2^shift (like System.Numerics.BigInteger).</summary>
    public static BigInt operator >>(BigInt value, int shift)
    {
        if (shift < 0)
        {
            return value << (-shift);
        }
        return value.ShiftRightFloor(shift);
    }

    /// <summary>Returns value * 2^shift.</summary>
    public BigInt ShiftLeft(long shift)
    {
        if (_size == 0 || shift == 0)
        {
            return this;
        }
        ArgumentOutOfRangeException.ThrowIfNegative(shift);
        int n = LimbCount;
        int limbShift = checked((int)(shift >> 6));
        int bitShift = (int)(shift & 63);
        int len = checked(n + limbShift + 1);
        ulong[] d = NewLimbs(len);
        d.AsSpan(0, limbShift).Clear();
        if (bitShift == 0)
        {
            Limbs.CopyTo(d.AsSpan(limbShift));
            d[len - 1] = 0;
        }
        else
        {
            d[n + limbShift] = Nat.ShiftLeft(d.AsSpan(limbShift, n), Limbs, n, bitShift);
        }
        return FromMagnitude(d, len, _size < 0 ? -1 : 1);
    }

    /// <summary>Divides by 2^shift with the quotient rounded as requested; the default (floor) matches <c>&gt;&gt;</c> (mpz_{t,f,c}div_q_2exp).</summary>
    public BigInt ShiftRight(long shift, IntegerRounding rounding = IntegerRounding.Floor) => rounding switch
    {
        IntegerRounding.Truncate => ShiftRightTruncate(shift),
        IntegerRounding.Floor => ShiftRightFloor(shift),
        IntegerRounding.Ceiling => ShiftRightCeiling(shift),
        _ => throw new ArgumentOutOfRangeException(nameof(rounding)),
    };

    /// <summary>Truncating right shift: |value| / 2^shift with the original sign (mpz_tdiv_q_2exp).</summary>
    internal BigInt ShiftRightTruncate(long shift)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(shift);
        if (_size == 0 || shift == 0)
        {
            return this;
        }
        int n = LimbCount;
        int limbShift = (int)Math.Min(shift >> 6, n);
        int bitShift = (int)(shift & 63);
        int len = n - limbShift;
        if (len <= 0)
        {
            return default;
        }
        ulong[] d = NewLimbs(len);
        if (bitShift == 0)
        {
            Limbs.Slice(limbShift).CopyTo(d);
        }
        else
        {
            Nat.ShiftRight(d, Limbs.Slice(limbShift), len, bitShift);
        }
        return FromMagnitude(d, len, _size < 0 ? -1 : 1);
    }

    /// <summary>Floor right shift: floor(value / 2^shift) (mpz_fdiv_q_2exp).</summary>
    internal BigInt ShiftRightFloor(long shift)
    {
        BigInt t = ShiftRightTruncate(shift);
        if (_size >= 0 || shift == 0)
        {
            return t;
        }
        // Negative values round toward -inf when any shifted-out bit was set.
        if (AnyBitsBelow(shift))
        {
            return t - One;
        }
        return t;
    }

    private bool AnyBitsBelow(long shift)
    {
        ReadOnlySpan<ulong> m = Limbs;
        int n = m.Length;
        long fullLimbs = Math.Min(shift >> 6, n);
        for (int i = 0; i < fullLimbs; i++)
        {
            if (m[i] != 0)
            {
                return true;
            }
        }
        int bitShift = (int)(shift & 63);
        long limbIndex = shift >> 6;
        if (bitShift != 0 && limbIndex < n && (m[(int)limbIndex] & ((1UL << bitShift) - 1)) != 0)
        {
            return true;
        }
        return false;
    }

    // ---------------------------------------------------------------
    // Powers
    // ---------------------------------------------------------------

    /// <summary>Returns value^exponent (mpz_pow_ui).</summary>
    public static BigInt Pow(BigInt value, ulong exponent)
    {
        if (exponent == 0)
        {
            return One;
        }
        if (exponent == 1 || value._size == 0 || value.IsOne)
        {
            return value;
        }
        if (value.IsMinusOne)
        {
            return (exponent & 1) == 0 ? One : MinusOne;
        }

        BigInt result = One;
        BigInt basePow = value;
        while (true)
        {
            if ((exponent & 1) != 0)
            {
                result = result.IsOne ? basePow : result * basePow;
            }
            exponent >>= 1;
            if (exponent == 0)
            {
                return result;
            }
            basePow = basePow.Square();
        }
    }

    /// <summary>Returns base^exponent for ulong base (mpz_ui_pow_ui).</summary>
    public static BigInt Pow(ulong value, ulong exponent) => Pow(new BigInt(value), exponent);
}
