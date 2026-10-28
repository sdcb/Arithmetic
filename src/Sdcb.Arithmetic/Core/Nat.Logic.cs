using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.Arm;
using System.Runtime.Intrinsics.X86;

namespace Sdcb.Arithmetic.Core;

internal static partial class Nat
{
    // ---------------------------------------------------------------
    // Bitwise logic over limb vectors: Vector512 → Vector256 → Vector128 → scalar,
    // selected through Isa (JIT-time constants). r may alias a or b exactly.
    // ---------------------------------------------------------------

    private interface ILimbOp
    {
        static abstract ulong Op(ulong a, ulong b);
        static abstract Vector128<ulong> Op(Vector128<ulong> a, Vector128<ulong> b);
        static abstract Vector256<ulong> Op(Vector256<ulong> a, Vector256<ulong> b);
        static abstract Vector512<ulong> Op(Vector512<ulong> a, Vector512<ulong> b);
    }

    private readonly struct AndOp : ILimbOp
    {
        public static ulong Op(ulong a, ulong b) => a & b;
        public static Vector128<ulong> Op(Vector128<ulong> a, Vector128<ulong> b) => a & b;
        public static Vector256<ulong> Op(Vector256<ulong> a, Vector256<ulong> b) => a & b;
        public static Vector512<ulong> Op(Vector512<ulong> a, Vector512<ulong> b) => a & b;
    }

    private readonly struct OrOp : ILimbOp
    {
        public static ulong Op(ulong a, ulong b) => a | b;
        public static Vector128<ulong> Op(Vector128<ulong> a, Vector128<ulong> b) => a | b;
        public static Vector256<ulong> Op(Vector256<ulong> a, Vector256<ulong> b) => a | b;
        public static Vector512<ulong> Op(Vector512<ulong> a, Vector512<ulong> b) => a | b;
    }

    private readonly struct XorOp : ILimbOp
    {
        public static ulong Op(ulong a, ulong b) => a ^ b;
        public static Vector128<ulong> Op(Vector128<ulong> a, Vector128<ulong> b) => a ^ b;
        public static Vector256<ulong> Op(Vector256<ulong> a, Vector256<ulong> b) => a ^ b;
        public static Vector512<ulong> Op(Vector512<ulong> a, Vector512<ulong> b) => a ^ b;
    }

    private readonly struct AndNotOp : ILimbOp
    {
        public static ulong Op(ulong a, ulong b) => a & ~b;
        public static Vector128<ulong> Op(Vector128<ulong> a, Vector128<ulong> b) => Vector128.AndNot(a, b);
        public static Vector256<ulong> Op(Vector256<ulong> a, Vector256<ulong> b) => Vector256.AndNot(a, b);
        public static Vector512<ulong> Op(Vector512<ulong> a, Vector512<ulong> b) => Vector512.AndNot(a, b);
    }

    private readonly struct OrNotOp : ILimbOp
    {
        public static ulong Op(ulong a, ulong b) => a | ~b;
        public static Vector128<ulong> Op(Vector128<ulong> a, Vector128<ulong> b) => a | ~b;
        public static Vector256<ulong> Op(Vector256<ulong> a, Vector256<ulong> b) => a | ~b;
        public static Vector512<ulong> Op(Vector512<ulong> a, Vector512<ulong> b) => a | ~b;
    }

    private readonly struct NandOp : ILimbOp
    {
        public static ulong Op(ulong a, ulong b) => ~(a & b);
        public static Vector128<ulong> Op(Vector128<ulong> a, Vector128<ulong> b) => ~(a & b);
        public static Vector256<ulong> Op(Vector256<ulong> a, Vector256<ulong> b) => ~(a & b);
        public static Vector512<ulong> Op(Vector512<ulong> a, Vector512<ulong> b) => ~(a & b);
    }

    private readonly struct NorOp : ILimbOp
    {
        public static ulong Op(ulong a, ulong b) => ~(a | b);
        public static Vector128<ulong> Op(Vector128<ulong> a, Vector128<ulong> b) => ~(a | b);
        public static Vector256<ulong> Op(Vector256<ulong> a, Vector256<ulong> b) => ~(a | b);
        public static Vector512<ulong> Op(Vector512<ulong> a, Vector512<ulong> b) => ~(a | b);
    }

    private readonly struct XnorOp : ILimbOp
    {
        public static ulong Op(ulong a, ulong b) => ~(a ^ b);
        public static Vector128<ulong> Op(Vector128<ulong> a, Vector128<ulong> b) => ~(a ^ b);
        public static Vector256<ulong> Op(Vector256<ulong> a, Vector256<ulong> b) => ~(a ^ b);
        public static Vector512<ulong> Op(Vector512<ulong> a, Vector512<ulong> b) => ~(a ^ b);
    }

    private static void Bitwise<TOp>(Span<ulong> r, ReadOnlySpan<ulong> a, ReadOnlySpan<ulong> b, int n) where TOp : struct, ILimbOp
    {
        if ((uint)n > (uint)r.Length || (uint)n > (uint)a.Length || (uint)n > (uint)b.Length)
        {
            throw new ArgumentOutOfRangeException(nameof(n));
        }
        ref ulong rp = ref MemoryMarshal.GetReference(r);
        ref ulong ap = ref MemoryMarshal.GetReference(a);
        ref ulong bp = ref MemoryMarshal.GetReference(b);
        nuint i = 0, len = (nuint)n;
        // Each tier consumes what the wider one left, so tails stay vectorized.
        if (Isa.Hw512)
        {
            for (; i + 8 <= len; i += 8)
            {
                TOp.Op(Vector512.LoadUnsafe(ref ap, i), Vector512.LoadUnsafe(ref bp, i)).StoreUnsafe(ref rp, i);
            }
        }
        if (Isa.Hw256)
        {
            for (; i + 4 <= len; i += 4)
            {
                TOp.Op(Vector256.LoadUnsafe(ref ap, i), Vector256.LoadUnsafe(ref bp, i)).StoreUnsafe(ref rp, i);
            }
        }
        if (Isa.Hw128)
        {
            for (; i + 2 <= len; i += 2)
            {
                TOp.Op(Vector128.LoadUnsafe(ref ap, i), Vector128.LoadUnsafe(ref bp, i)).StoreUnsafe(ref rp, i);
            }
        }
        for (; i < len; i++)
        {
            Unsafe.Add(ref rp, i) = TOp.Op(Unsafe.Add(ref ap, i), Unsafe.Add(ref bp, i));
        }
    }

    /// <summary>r[0..n) = a[0..n) &amp; b[0..n).</summary>
    public static void And(Span<ulong> r, ReadOnlySpan<ulong> a, ReadOnlySpan<ulong> b, int n) => Bitwise<AndOp>(r, a, b, n);

    /// <summary>r[0..n) = a[0..n) | b[0..n).</summary>
    public static void Or(Span<ulong> r, ReadOnlySpan<ulong> a, ReadOnlySpan<ulong> b, int n) => Bitwise<OrOp>(r, a, b, n);

    /// <summary>r[0..n) = a[0..n) ^ b[0..n).</summary>
    public static void Xor(Span<ulong> r, ReadOnlySpan<ulong> a, ReadOnlySpan<ulong> b, int n) => Bitwise<XorOp>(r, a, b, n);

    /// <summary>r[0..n) = a[0..n) &amp; ~b[0..n).</summary>
    public static void AndNot(Span<ulong> r, ReadOnlySpan<ulong> a, ReadOnlySpan<ulong> b, int n) => Bitwise<AndNotOp>(r, a, b, n);

    /// <summary>r[0..n) = a[0..n) | ~b[0..n).</summary>
    public static void OrNot(Span<ulong> r, ReadOnlySpan<ulong> a, ReadOnlySpan<ulong> b, int n) => Bitwise<OrNotOp>(r, a, b, n);

    /// <summary>r[0..n) = ~(a[0..n) &amp; b[0..n)).</summary>
    public static void Nand(Span<ulong> r, ReadOnlySpan<ulong> a, ReadOnlySpan<ulong> b, int n) => Bitwise<NandOp>(r, a, b, n);

    /// <summary>r[0..n) = ~(a[0..n) | b[0..n)).</summary>
    public static void Nor(Span<ulong> r, ReadOnlySpan<ulong> a, ReadOnlySpan<ulong> b, int n) => Bitwise<NorOp>(r, a, b, n);

    /// <summary>r[0..n) = ~(a[0..n) ^ b[0..n)).</summary>
    public static void Xnor(Span<ulong> r, ReadOnlySpan<ulong> a, ReadOnlySpan<ulong> b, int n) => Bitwise<XnorOp>(r, a, b, n);

    /// <summary>r[0..n) = ~a[0..n). r may alias a.</summary>
    public static void Not(Span<ulong> r, ReadOnlySpan<ulong> a, int n) => Bitwise<NandOp>(r, a, a, n);

    // ---------------------------------------------------------------
    // Population count: AVX-512BW/AVX2/SSSE3 nibble-LUT + PSADBW, AdvSimd CNT,
    // scalar POPCNT fallback.
    // ---------------------------------------------------------------

    /// <summary>Total number of set bits across a[0..n).</summary>
    public static long PopCount(ReadOnlySpan<ulong> a, int n) => PopCountCore(a, a, n, xor: false);

    /// <summary>Number of bit positions where a and b differ (equal-length inputs).</summary>
    public static long HammingDistance(ReadOnlySpan<ulong> a, ReadOnlySpan<ulong> b, int n) => PopCountCore(a, b, n, xor: true);

    private static long PopCountCore(ReadOnlySpan<ulong> a, ReadOnlySpan<ulong> b, int n, bool xor)
    {
        if ((uint)n > (uint)a.Length || (uint)n > (uint)b.Length)
        {
            throw new ArgumentOutOfRangeException(nameof(n));
        }
        ref ulong ap = ref MemoryMarshal.GetReference(a);
        ref ulong bp = ref MemoryMarshal.GetReference(b);
        nuint i = 0, len = (nuint)n;
        long count = 0;
        // Tiers cascade: each consumes what the wider one left.
        if (Isa.Hw512 && len - i >= 16)
        {
            Vector512<byte> lut = Vector512.Create(PopLut64);
            Vector512<byte> low = Vector512.Create((byte)0x0F);
            Vector512<ulong> acc = Vector512<ulong>.Zero;
            for (; i + 8 <= len; i += 8)
            {
                Vector512<ulong> v = Vector512.LoadUnsafe(ref ap, i);
                if (xor)
                {
                    v ^= Vector512.LoadUnsafe(ref bp, i);
                }
                Vector512<byte> bytes = v.AsByte();
                Vector512<byte> c = Avx512BW.Shuffle(lut, bytes & low) + Avx512BW.Shuffle(lut, Avx512BW.ShiftRightLogical(bytes.AsUInt16(), 4).AsByte() & low);
                acc += Avx512BW.SumAbsoluteDifferences(c, Vector512<byte>.Zero).AsUInt64();
            }
            count += (long)Vector512.Sum(acc);
        }
        if (Isa.Hw256 && len - i >= 4)
        {
            Vector256<byte> lut = Vector256.Create(PopLut64.AsSpan(0, 32));
            Vector256<byte> low = Vector256.Create((byte)0x0F);
            Vector256<ulong> acc = Vector256<ulong>.Zero;
            for (; i + 4 <= len; i += 4)
            {
                Vector256<ulong> v = Vector256.LoadUnsafe(ref ap, i);
                if (xor)
                {
                    v ^= Vector256.LoadUnsafe(ref bp, i);
                }
                Vector256<byte> bytes = v.AsByte();
                Vector256<byte> c = Avx2.Shuffle(lut, bytes & low) + Avx2.Shuffle(lut, Avx2.ShiftRightLogical(bytes.AsUInt16(), 4).AsByte() & low);
                acc += Avx2.SumAbsoluteDifferences(c, Vector256<byte>.Zero).AsUInt64();
            }
            count += (long)Vector256.Sum(acc);
        }
        if (Isa.Hw128 && AdvSimd.Arm64.IsSupported)
        {
            for (; i + 2 <= len; i += 2)
            {
                Vector128<ulong> v = Vector128.LoadUnsafe(ref ap, i);
                if (xor)
                {
                    v ^= Vector128.LoadUnsafe(ref bp, i);
                }
                count += AdvSimd.Arm64.AddAcross(AdvSimd.PopCount(v.AsByte())).ToScalar();
            }
        }
        else if (Isa.Hw128 && Ssse3.IsSupported && len - i >= 4)
        {
            Vector128<byte> lut = Vector128.Create(PopLut64.AsSpan(0, 16));
            Vector128<byte> low = Vector128.Create((byte)0x0F);
            Vector128<ulong> acc = Vector128<ulong>.Zero;
            for (; i + 2 <= len; i += 2)
            {
                Vector128<ulong> v = Vector128.LoadUnsafe(ref ap, i);
                if (xor)
                {
                    v ^= Vector128.LoadUnsafe(ref bp, i);
                }
                Vector128<byte> bytes = v.AsByte();
                Vector128<byte> c = Ssse3.Shuffle(lut, bytes & low) + Ssse3.Shuffle(lut, Sse2.ShiftRightLogical(bytes.AsUInt16(), 4).AsByte() & low);
                acc += Sse2.SumAbsoluteDifferences(c, Vector128<byte>.Zero).AsUInt64();
            }
            count += (long)Vector128.Sum(acc);
        }
        for (; i < len; i++)
        {
            ulong v = Unsafe.Add(ref ap, i);
            if (xor)
            {
                v ^= Unsafe.Add(ref bp, i);
            }
            count += BitOperations.PopCount(v);
        }
        return count;
    }

    private static readonly byte[] PopLut64 =
    [
        0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
        0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
        0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
        0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
    ];

    /// <summary>Index of the lowest set bit at or above position bit; -1 if none below n limbs.</summary>
    public static long IndexOfSetBit(ReadOnlySpan<ulong> a, int n, long bit)
    {
        if (bit >> 6 >= n)
        {
            return -1;
        }
        int limb = (int)(bit >> 6);
        ulong w = a[limb] & (ulong.MaxValue << (int)(bit & 63));
        if (w != 0)
        {
            return ((long)limb << 6) + BitOperations.TrailingZeroCount(w);
        }
        int idx = a.Slice(limb + 1, n - limb - 1).IndexOfAnyExcept(0UL);
        if (idx < 0)
        {
            return -1;
        }
        limb += 1 + idx;
        return ((long)limb << 6) + BitOperations.TrailingZeroCount(a[limb]);
    }

    /// <summary>Index of the lowest clear bit at or above position bit (treating the value as infinitely zero-extended).</summary>
    public static long IndexOfClearBit(ReadOnlySpan<ulong> a, int n, long bit)
    {
        if (bit >> 6 >= n)
        {
            return bit;
        }
        int limb = (int)(bit >> 6);
        ulong w = ~a[limb] & (ulong.MaxValue << (int)(bit & 63));
        if (w != 0)
        {
            return ((long)limb << 6) + BitOperations.TrailingZeroCount(w);
        }
        int idx = a.Slice(limb + 1, n - limb - 1).IndexOfAnyExcept(ulong.MaxValue);
        if (idx < 0)
        {
            return (long)n << 6;
        }
        limb += 1 + idx;
        return ((long)limb << 6) + BitOperations.TrailingZeroCount(~a[limb]);
    }

    /// <summary>Total bit length of a normalized magnitude (0 for zero).</summary>
    public static long BitLength(ReadOnlySpan<ulong> a, int n)
    {
        if (n == 0)
        {
            return 0;
        }
        return ((long)(n - 1) << 6) + (64 - BitOperations.LeadingZeroCount(a[n - 1]));
    }

    /// <summary>Number of trailing zero bits (input must be nonzero).</summary>
    public static long TrailingZeroCount(ReadOnlySpan<ulong> a, int n)
    {
        int i = a.Slice(0, n).IndexOfAnyExcept(0UL);
        return ((long)i << 6) + BitOperations.TrailingZeroCount(a[i]);
    }
}
