using System.Buffers;
using System.Diagnostics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

namespace Sdcb.Arithmetic.Core;

internal static partial class Nat
{
    /// <summary>Operand size (in limbs) at which Karatsuba beats the basecase.</summary>
    public static readonly int MulKaratsubaThreshold = EnvInt("SDCB_ARITHMETIC_KARA_MUL") ?? (Isa.Hw256 ? 80 : 24);

    /// <summary>Balanced operand size (in limbs) at which Toom-3 beats Karatsuba.</summary>
    public static readonly int MulToom3Threshold = EnvInt("SDCB_ARITHMETIC_TOOM3") ?? 128;

    /// <summary>Smaller-operand size (limbs) from which the NTT multiplication is used; depends on the SIMD tier.</summary>
    internal static int MulNttThreshold => s_nttMulEnv ?? NttThresholdForTier(3300);

    /// <summary>Operand size (limbs) from which squaring uses the NTT.</summary>
    internal static int SqrNttThreshold => s_nttSqrEnv ?? NttThresholdForTier(3500);

    private static readonly int? s_nttMulEnv = EnvInt("SDCB_ARITHMETIC_NTT_MUL");
    private static readonly int? s_nttSqrEnv = EnvInt("SDCB_ARITHMETIC_NTT_SQR");

    // The 32-bit-lane transform only beats Toom-3 once it is vectorized; with AVX2 the vector
    // basecase also speeds up Toom-3, which pushes the crossover up (measured on Zen 3).
    private static int NttThresholdForTier(int hw256) => Isa.Hw512 ? hw256 * 3 / 4 : Isa.Hw256 ? hw256 : Isa.Hw128 ? 1500 : 40000;

    private static int? EnvInt(string name)
        => int.TryParse(Environment.GetEnvironmentVariable(name), out int v) ? v : null;

    // ---------------------------------------------------------------
    // ---------------------------------------------------------------
    // Full multiplication
    // ---------------------------------------------------------------

    /// <summary>
    /// r[0..an+bn) = a[0..an) * b[0..bn); requires an &gt;= bn &gt;= 1 and r not
    /// overlapping a or b. Dispatches to basecase / Karatsuba / Toom-3.
    /// </summary>
    public static void Mul(Span<ulong> r, ReadOnlySpan<ulong> a, int an, ReadOnlySpan<ulong> b, int bn)
    {
        Debug.Assert(an >= bn && bn >= 1);
        if (bn == 1)
        {
            r[an] = MulScalar(r, a, an, b[0]);
            return;
        }
        if (bn <= MulKaratsubaThreshold)
        {
            MulSchoolbook(r, a, an, b, bn);
            return;
        }
        int need = MulScratchSize(an, bn);
        ulong[] pool = ArrayPool<ulong>.Shared.Rent(need);
        MulRecursive(r, a, an, b, bn, pool.AsSpan(0, need));
        ArrayPool<ulong>.Shared.Return(pool);
    }

    /// <summary>Conservative scratch bound for <see cref="MulRecursive"/>.</summary>
    public static int MulScratchSize(int an, int bn) => 6 * (an + bn) + 128;


    internal static void MulRecursive(Span<ulong> r, ReadOnlySpan<ulong> a, int an, ReadOnlySpan<ulong> b, int bn, Span<ulong> scratch)
    {
        Debug.Assert(an >= bn && bn >= 1);
        if (bn == 1)
        {
            r[an] = MulScalar(r, a, an, b[0]);
            return;
        }
        if (bn <= MulKaratsubaThreshold)
        {
            MulSchoolbook(r, a, an, b, bn);
            return;
        }

        int k = (an + 1) >> 1;
        if (bn >= MulNttThreshold && Ntt.Fits(an, bn))
        {
            Ntt.Mul(r, a, an, b, bn);
        }
        else if (bn <= k)
        {
            MulUnbalanced(r, a, an, b, bn, scratch);
        }
        else if (an == bn && an >= MulToom3Threshold)
        {
            Toom3Mul(r, a, b, an, scratch);
        }
        else
        {
            KaratsubaMul(r, a, an, b, bn, scratch);
        }
    }

    /// <summary>Handles an &gt;&gt; bn by multiplying bn-sized chunks of a with b.</summary>
    private static void MulUnbalanced(Span<ulong> r, ReadOnlySpan<ulong> a, int an, ReadOnlySpan<ulong> b, int bn, Span<ulong> scratch)
    {
        Span<ulong> tmp = scratch.Slice(0, 2 * bn);
        Span<ulong> rest = scratch.Slice(2 * bn);

        // First chunk writes directly into r.
        MulRecursive(r.Slice(0, 2 * bn), a.Slice(0, bn), bn, b, bn, rest);
        int ofs = bn;
        while (an - ofs >= bn)
        {
            MulRecursive(tmp, a.Slice(ofs, bn), bn, b, bn, rest);
            ulong cy = Add(r.Slice(ofs), r.Slice(ofs), tmp, bn);
            tmp.Slice(bn, bn).CopyTo(r.Slice(ofs + bn));
            AddScalarInPlace(r.Slice(ofs + bn), bn, cy);
            ofs += bn;
        }
        int rem = an - ofs;
        if (rem > 0)
        {
            // b (bn limbs) x tail (rem < bn limbs); swap so first operand is longer.
            MulRecursive(tmp.Slice(0, bn + rem), b, bn, a.Slice(ofs, rem), rem, rest);
            ulong cy = Add(r.Slice(ofs), r.Slice(ofs), tmp, bn);
            tmp.Slice(bn, rem).CopyTo(r.Slice(ofs + bn));
            AddScalarInPlace(r.Slice(ofs + bn), rem, cy);
        }
    }

    /// <summary>Karatsuba for an &gt;= bn &gt; ceil(an/2).</summary>
    private static void KaratsubaMul(Span<ulong> r, ReadOnlySpan<ulong> a, int an, ReadOnlySpan<ulong> b, int bn, Span<ulong> scratch)
    {
        int k = (an + 1) >> 1;
        int h = an - k; // 1 <= h <= k
        int l = bn - k; // 1 <= l <= h
        Debug.Assert(bn > k);

        Span<ulong> da = scratch.Slice(0, k);
        Span<ulong> db = scratch.Slice(k, k);
        Span<ulong> z1 = scratch.Slice(2 * k, 2 * k);
        Span<ulong> t2 = scratch.Slice(4 * k, 2 * k);
        Span<ulong> rest = scratch.Slice(6 * k);

        int sa = AbsSub(da, a.Slice(0, k), k, a.Slice(k), h);
        int sb = AbsSub(db, b.Slice(0, k), k, b.Slice(k), l);

        // z1 = |a0-a1| * |b0-b1|
        MulRecursive(z1, da, k, db, k, rest);
        // z0 = a0*b0 into r[0..2k)
        MulRecursive(r.Slice(0, 2 * k), a.Slice(0, k), k, b.Slice(0, k), k, rest);
        // z2 = a1*b1 into r[2k..an+bn)
        MulRecursive(r.Slice(2 * k, h + l), a.Slice(k), h, b.Slice(k), l, rest);

        // t2 = z0 + z2 (2k limbs), carry cc
        ulong cc = Add(t2, r.Slice(0, 2 * k), 2 * k, r.Slice(2 * k, h + l), h + l);
        if (sa == sb)
        {
            // middle = z0 + z2 - z1 >= 0
            cc -= SubInPlace(t2, 2 * k, z1, 2 * k);
        }
        else
        {
            cc += Add(t2, t2, z1, 2 * k);
        }

        // r[k..] += (t2, cc)
        int rn = an + bn;
        ulong cy = Add(r.Slice(k), r.Slice(k), t2, 2 * k);
        cy += cc;
        if (3 * k < rn)
        {
            cy = AddScalarInPlace(r.Slice(3 * k), rn - 3 * k, cy);
        }
        Debug.Assert(cy == 0);
    }

    /// <summary>Balanced Toom-3 multiplication for an == bn == n. b may alias a (squaring).</summary>
    private static void Toom3Mul(Span<ulong> r, ReadOnlySpan<ulong> a, ReadOnlySpan<ulong> b, int n, Span<ulong> scratch)
    {
        int k = (n + 2) / 3;
        int s = n - 2 * k; // 0 < s <= k
        Debug.Assert(s > 0 && s <= k);
        int e = k + 1;   // evaluation size
        int L = 2 * e;   // product size

        ReadOnlySpan<ulong> a0 = a.Slice(0, k), a1 = a.Slice(k, k), a2 = a.Slice(2 * k, s);
        ReadOnlySpan<ulong> b0 = b.Slice(0, k), b1 = b.Slice(k, k), b2 = b.Slice(2 * k, s);

        Span<ulong> ta = scratch.Slice(0, e);
        Span<ulong> tb = scratch.Slice(e, e);
        Span<ulong> e1a = scratch.Slice(2 * e, e);
        Span<ulong> e1b = scratch.Slice(3 * e, e);
        Span<ulong> em1a = scratch.Slice(4 * e, e);
        Span<ulong> em1b = scratch.Slice(5 * e, e);
        Span<ulong> e2a = scratch.Slice(6 * e, e);
        Span<ulong> e2b = scratch.Slice(7 * e, e);
        Span<ulong> v1 = scratch.Slice(8 * e, L);
        Span<ulong> vm1 = scratch.Slice(8 * e + L, L);
        Span<ulong> v2 = scratch.Slice(8 * e + 2 * L, L);
        Span<ulong> rest = scratch.Slice(8 * e + 3 * L);

        // ta = a0 + a2 ; e1a = ta + a1 ; em1a = |ta - a1| (sign sa)
        bool square = Unsafe.AreSame(ref MemoryMarshal.GetReference(a), ref MemoryMarshal.GetReference(b));
        ta[k] = Add(ta, a0, k, a2, s);
        e1a[k] = ta[k] + Add(e1a, ta.Slice(0, k), k, a1, k);
        int sa = AbsSub(em1a, ta, e, a1, k);
        // e2a = ((a2*2 + a1)*2 + a0)
        Eval2(e2a, a0, a1, a2, k, s);

        int svm1;
        if (square)
        {
            svm1 = 1;
            SquareRecursive(v1, e1a, e, rest);
            SquareRecursive(vm1, em1a, e, rest);
            SquareRecursive(v2, e2a, e, rest);
            SquareRecursive(r.Slice(0, 2 * k), a0, k, rest);
            SquareRecursive(r.Slice(4 * k, 2 * s), a2, s, rest);
        }
        else
        {
            tb[k] = Add(tb, b0, k, b2, s);
            e1b[k] = tb[k] + Add(e1b, tb.Slice(0, k), k, b1, k);
            int sb = AbsSub(em1b, tb, e, b1, k);
            Eval2(e2b, b0, b1, b2, k, s);
            svm1 = sa * sb;

            MulRecursive(v1, e1a, e, e1b, e, rest);
            MulRecursive(vm1, em1a, e, em1b, e, rest);
            MulRecursive(v2, e2a, e, e2b, e, rest);
            // v0 = a0*b0 into r[0..2k); vinf = a2*b2 into r[4k..2n)
            MulRecursive(r.Slice(0, 2 * k), a0, k, b0, k, rest);
            MulRecursive(r.Slice(4 * k, 2 * s), a2, s, b2, s, rest);
        }

        ReadOnlySpan<ulong> v0 = r.Slice(0, 2 * k);
        ReadOnlySpan<ulong> vinf = r.Slice(4 * k, 2 * s);

        // Interpolation (Bodrato's sequence for 5 points):
        // v2 <- (v2 - vm1) / 3
        if (svm1 >= 0)
        {
            SubInPlace(v2, L, vm1, L);
        }
        else
        {
            AddInPlace(v2, L, vm1, L);
        }
        DivExactBy3InPlace(v2, L);
        // vm1 <- (v1 - vm1) / 2   (== c1 + c3 >= 0)
        if (svm1 >= 0)
        {
            Sub(vm1, v1, vm1, L);
        }
        else
        {
            Add(vm1, v1, vm1, L);
        }
        ShiftRightOneInPlace(vm1, L);
        // v1 <- v1 - v0
        SubInPlace(v1, L, v0, 2 * k);
        // v2 <- (v2 - v1)/2 - 2*vinf   (== c3)
        SubInPlace(v2, L, v1, L);
        ShiftRightOneInPlace(v2, L);
        SubInPlace(v2, L, vinf, 2 * s);
        SubInPlace(v2, L, vinf, 2 * s);
        // v1 <- v1 - vm1 - vinf        (== c2)
        SubInPlace(v1, L, vm1, L);
        SubInPlace(v1, L, vinf, 2 * s);
        // vm1 <- vm1 - v2              (== c1)
        SubInPlace(vm1, L, v2, L);

        // Recomposition: r = v0 + vm1*B^k + v1*B^2k + v2*B^3k + vinf*B^4k
        int rn = 2 * n;
        r.Slice(2 * k, 2 * k).Clear();
        ulong cy;
        int len = TrimmedLength(vm1, L);
        Debug.Assert(len <= rn - k);
        cy = AddInPlace(r.Slice(k), rn - k, vm1, len);
        Debug.Assert(cy == 0);
        len = TrimmedLength(v1, L);
        Debug.Assert(len <= rn - 2 * k);
        cy = AddInPlace(r.Slice(2 * k), rn - 2 * k, v1, len);
        Debug.Assert(cy == 0);
        len = TrimmedLength(v2, L);
        Debug.Assert(len <= rn - 3 * k);
        cy = AddInPlace(r.Slice(3 * k), rn - 3 * k, v2, len);
        Debug.Assert(cy == 0);
    }

    /// <summary>e2 = a0 + 2*a1 + 4*a2 over e = k+1 limbs.</summary>
    private static void Eval2(Span<ulong> e2, ReadOnlySpan<ulong> a0, ReadOnlySpan<ulong> a1, ReadOnlySpan<ulong> a2, int k, int s)
    {
        int e = k + 1;
        a2.CopyTo(e2);
        e2.Slice(s, e - s).Clear();
        // e2 = 2*e2 + a1
        e2[k] = (e2[k] << 1) | ShiftLeft(e2.Slice(0, k), e2.Slice(0, k), k, 1);
        AddInPlace(e2, e, a1, k);
        // e2 = 2*e2 + a0
        ulong outBit = ShiftLeft(e2, e2, e, 1);
        Debug.Assert(outBit == 0);
        AddInPlace(e2, e, a0, k);
    }

    // ---------------------------------------------------------------
    // Squaring
    // ---------------------------------------------------------------

    /// <summary>Squaring size threshold below which the basecase is used.</summary>
    public static readonly int SqrKaratsubaThreshold = EnvInt("SDCB_ARITHMETIC_KARA_SQR") ?? (Isa.Hw256 ? 128 : 32);

    /// <summary>Size (limbs) from which squaring uses Toom-3.</summary>
    public static readonly int SqrToom3Threshold = EnvInt("SDCB_ARITHMETIC_TOOM3_SQR") ?? (Isa.Hw256 ? 384 : 128);

    /// <summary>r[0..2n) = a[0..n)^2; r must not overlap a.</summary>
    public static void Square(Span<ulong> r, ReadOnlySpan<ulong> a, int n)
    {
        Debug.Assert(n >= 1);
        if (n <= SqrKaratsubaThreshold)
        {
            SquareSchoolbook(r, a, n);
            return;
        }
        int need = MulScratchSize(n, n);
        ulong[] pool = ArrayPool<ulong>.Shared.Rent(need);
        SquareRecursive(r, a, n, pool.AsSpan(0, need));
        ArrayPool<ulong>.Shared.Return(pool);
    }

    internal static void SquareRecursive(Span<ulong> r, ReadOnlySpan<ulong> a, int n, Span<ulong> scratch)
    {
        if (n <= SqrKaratsubaThreshold)
        {
            SquareSchoolbook(r, a, n);
            return;
        }
        if (n >= SqrNttThreshold && Ntt.Fits(n, n))
        {
            Ntt.Mul(r, a, n, a, n);
            return;
        }
        if (n >= SqrToom3Threshold)
        {
            Toom3Mul(r, a, a, n, scratch);
            return;
        }

        // Karatsuba squaring: middle = z0 + z2 - (a0-a1)^2 (always subtract).
        int k = (n + 1) >> 1;
        int h = n - k;
        Span<ulong> da = scratch.Slice(0, k);
        Span<ulong> z1 = scratch.Slice(k, 2 * k);
        Span<ulong> t2 = scratch.Slice(3 * k, 2 * k);
        Span<ulong> rest = scratch.Slice(5 * k);

        AbsSub(da, a.Slice(0, k), k, a.Slice(k), h);
        SquareRecursive(z1, da, k, rest);
        SquareRecursive(r.Slice(0, 2 * k), a.Slice(0, k), k, rest);
        SquareRecursive(r.Slice(2 * k, 2 * h), a.Slice(k), h, rest);

        ulong cc = Add(t2, r.Slice(0, 2 * k), 2 * k, r.Slice(2 * k, 2 * h), 2 * h);
        cc -= SubInPlace(t2, 2 * k, z1, 2 * k);

        int rn = 2 * n;
        ulong cy = Add(r.Slice(k), r.Slice(k), t2, 2 * k);
        cy += cc;
        if (3 * k < rn)
        {
            cy = AddScalarInPlace(r.Slice(3 * k), rn - 3 * k, cy);
        }
        Debug.Assert(cy == 0);
    }

    // ---------------------------------------------------------------
    // Helpers
    // ---------------------------------------------------------------

    /// <summary>d[0..an) = |A - B| where A = a[0..an), B = b[0..bn), an &gt;= bn. Returns sign of A-B.</summary>
    public static int AbsSub(Span<ulong> d, ReadOnlySpan<ulong> a, int an, ReadOnlySpan<ulong> b, int bn)
    {
        int c = CompareLonger(a, an, b, bn);
        if (c >= 0)
        {
            Sub(d, a, an, b, bn);
            return c == 0 ? 0 : 1;
        }
        // B > A implies a's limbs above bn are all zero.
        Sub(d, b, a, bn);
        d.Slice(bn, an - bn).Clear();
        return -1;
    }

    private const ulong Inv3 = 0xAAAAAAAAAAAAAAABUL;       // 3^-1 mod 2^64
    private const ulong OneThirdFloor = 0x5555555555555555UL; // floor((2^64-1)/3)

    /// <summary>In-place exact division by 3 (value must be divisible by 3).</summary>
    public static void DivExactBy3InPlace(Span<ulong> r, int n)
    {
        ulong c = 0;
        for (int i = 0; i < n; i++)
        {
            ulong l = r[i];
            ulong s = l - c;
            c = l < c ? 1UL : 0UL;
            l = s * Inv3;
            r[i] = l;
            if (l > OneThirdFloor) c++;
            if (l > OneThirdFloor * 2) c++;
        }
        Debug.Assert(c == 0);
    }

    /// <summary>In-place right shift by exactly one bit; the shifted-out bit must be 0 for exact halving.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static void ShiftRightOneInPlace(Span<ulong> r, int n)
    {
        ShiftRight(r, r, n, 1);
    }
}
