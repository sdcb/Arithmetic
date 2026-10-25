using System.Diagnostics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics.X86;

namespace Sdcb.Arithmetic.Core;

internal static partial class Nat
{
    /// <summary>
    /// 64x64→128 multiply. On x64 uses IMUL for the low half and single-output MULX for
    /// the high half: Math.BigMul's two-output MULX forces the low word through memory.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal static ulong BigMul(ulong a, ulong b, out ulong lo)
    {
        if (Bmi2.X64.IsSupported)
        {
            lo = a * b;
            return Bmi2.X64.MultiplyNoFlags(a, b);
        }
        return Math.BigMul(a, b, out lo);
    }

    // ---------------------------------------------------------------
    // Single-limb multiply primitives. Carry detection in managed code costs
    // cmp+setb per add, so every kernel is ordered to keep only one add+compare
    // on the loop-carried carry chain; the rest runs in parallel.
    // ---------------------------------------------------------------

    /// <summary>r[0..n) = a[0..n) * m. Returns the high limb of the product.</summary>
    public static ulong MulScalar(Span<ulong> r, ReadOnlySpan<ulong> a, int n, ulong m)
    {
        Debug.Assert(r.Length >= n && a.Length >= n);
        ref ulong rp = ref MemoryMarshal.GetReference(r);
        ref ulong ap = ref MemoryMarshal.GetReference(a);
        ulong carry = 0;
        nint i = 0, len = n;
        for (; i + 2 <= len; i += 2)
        {
            ulong hi0 = BigMul(Unsafe.Add(ref ap, i), m, out ulong lo0);
            ulong hi1 = BigMul(Unsafe.Add(ref ap, i + 1), m, out ulong lo1);
            lo0 += carry;
            hi0 += lo0 < carry ? 1UL : 0UL;
            lo1 += hi0;
            hi1 += lo1 < hi0 ? 1UL : 0UL;
            Unsafe.Add(ref rp, i) = lo0;
            Unsafe.Add(ref rp, i + 1) = lo1;
            carry = hi1;
        }
        if (i < len)
        {
            ulong hi = BigMul(Unsafe.Add(ref ap, i), m, out ulong lo);
            lo += carry;
            hi += lo < carry ? 1UL : 0UL;
            Unsafe.Add(ref rp, i) = lo;
            carry = hi;
        }
        return carry;
    }

    /// <summary>r[0..n) += a[0..n) * m. Returns the carry limb.</summary>
    public static ulong AddMulScalar(Span<ulong> r, ReadOnlySpan<ulong> a, int n, ulong m)
    {
        Debug.Assert(r.Length >= n && a.Length >= n);
        return AddMulScalar(ref MemoryMarshal.GetReference(r), ref MemoryMarshal.GetReference(a), n, m);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static ulong AddMulScalar(ref ulong rp, ref ulong ap, nint len, ulong m)
    {
        ulong carry = 0;
        nint i = 0;
        for (; i + 2 <= len; i += 2)
        {
            ulong hi0 = BigMul(Unsafe.Add(ref ap, i), m, out ulong lo0);
            ulong hi1 = BigMul(Unsafe.Add(ref ap, i + 1), m, out ulong lo1);
            ref ulong r0 = ref Unsafe.Add(ref rp, i);
            ref ulong r1 = ref Unsafe.Add(ref rp, i + 1);
            // Off the carry chain: fold r[i] into the product first.
            ulong t0 = r0, t1 = r1;
            lo0 += t0;
            hi0 += lo0 < t0 ? 1UL : 0UL;
            lo1 += t1;
            hi1 += lo1 < t1 ? 1UL : 0UL;
            // On the chain: a single add + compare per limb.
            lo0 += carry;
            hi0 += lo0 < carry ? 1UL : 0UL;
            lo1 += hi0;
            hi1 += lo1 < hi0 ? 1UL : 0UL;
            r0 = lo0;
            r1 = lo1;
            carry = hi1;
        }
        if (i < len)
        {
            ulong hi = BigMul(Unsafe.Add(ref ap, i), m, out ulong lo);
            ref ulong r0 = ref Unsafe.Add(ref rp, i);
            ulong t = r0;
            lo += t;
            hi += lo < t ? 1UL : 0UL;
            lo += carry;
            hi += lo < carry ? 1UL : 0UL;
            r0 = lo;
            carry = hi;
        }
        return carry;
    }

    /// <summary>r[0..n) -= a[0..n) * m. Returns the borrow limb.</summary>
    public static ulong SubMulScalar(Span<ulong> r, ReadOnlySpan<ulong> a, int n, ulong m)
    {
        Debug.Assert(r.Length >= n && a.Length >= n);
        ref ulong rp = ref MemoryMarshal.GetReference(r);
        ref ulong ap = ref MemoryMarshal.GetReference(a);
        ulong carry = 0;
        nint i = 0, len = n;
        for (; i + 2 <= len; i += 2)
        {
            ulong hi0 = BigMul(Unsafe.Add(ref ap, i), m, out ulong lo0);
            ulong hi1 = BigMul(Unsafe.Add(ref ap, i + 1), m, out ulong lo1);
            ref ulong r0 = ref Unsafe.Add(ref rp, i);
            ref ulong r1 = ref Unsafe.Add(ref rp, i + 1);
            ulong t0 = r0, t1 = r1;
            ulong d0 = t0 - lo0;
            hi0 += t0 < lo0 ? 1UL : 0UL;
            ulong d1 = t1 - lo1;
            hi1 += t1 < lo1 ? 1UL : 0UL;
            ulong e0 = d0 - carry;
            hi0 += d0 < carry ? 1UL : 0UL;
            ulong e1 = d1 - hi0;
            hi1 += d1 < hi0 ? 1UL : 0UL;
            r0 = e0;
            r1 = e1;
            carry = hi1;
        }
        if (i < len)
        {
            ulong hi = BigMul(Unsafe.Add(ref ap, i), m, out ulong lo);
            ref ulong r0 = ref Unsafe.Add(ref rp, i);
            ulong t = r0;
            ulong d = t - lo;
            hi += t < lo ? 1UL : 0UL;
            ulong e = d - carry;
            hi += d < carry ? 1UL : 0UL;
            r0 = e;
            carry = hi;
        }
        return carry;
    }

    // ---------------------------------------------------------------
    // Basecase multiplication / squaring: product scanning (Comba). Each
    // output column is summed into a 3-limb accumulator whose three short
    // dependency chains overlap across products, instead of one long
    // row-wise carry chain per AddMulScalar pass.
    // ---------------------------------------------------------------

    /// <summary>Schoolbook multiplication; r must not overlap operands. Requires an &gt;= bn &gt;= 1.</summary>
    public static void MulSchoolbook(Span<ulong> r, ReadOnlySpan<ulong> a, int an, ReadOnlySpan<ulong> b, int bn)
    {
        Debug.Assert(an >= bn && bn >= 1 && r.Length >= an + bn && a.Length >= an && b.Length >= bn);
        if (Isa.Hw256 && bn >= MulVecMin && bn < 1 << 14)
        {
            MulSchoolbookVec256(r, a, an, b, bn);
            return;
        }
        MulSchoolbookScalar(r, a, an, b, bn);
    }

    /// <summary>Scalar product-scanning schoolbook multiplication (same contract as <see cref="MulSchoolbook"/>).</summary>
    internal static void MulSchoolbookScalar(Span<ulong> r, ReadOnlySpan<ulong> a, int an, ReadOnlySpan<ulong> b, int bn)
    {
        ref ulong rp = ref MemoryMarshal.GetReference(r);
        ref ulong ap = ref MemoryMarshal.GetReference(a);
        ref ulong bp = ref MemoryMarshal.GetReference(b);
        ulong c0 = 0, c1 = 0, c2 = 0;
        nint rn = an + bn;
        for (nint k = 0; k < rn - 1; k++)
        {
            nint jlo = k - an + 1;
            if (jlo < 0)
            {
                jlo = 0;
            }
            nint jhi = k < bn - 1 ? k : bn - 1;
            ref ulong ak = ref Unsafe.Add(ref ap, k); // a[k - j] == Unsafe.Subtract(ref ak, j)
            nint j = jlo;
            for (; j + 1 <= jhi; j += 2)
            {
                ulong h0 = BigMul(Unsafe.Subtract(ref ak, j), Unsafe.Add(ref bp, j), out ulong l0);
                ulong h1 = BigMul(Unsafe.Subtract(ref ak, j + 1), Unsafe.Add(ref bp, j + 1), out ulong l1);
                c0 += l0;
                h0 += c0 < l0 ? 1UL : 0UL;
                c1 += h0;
                c2 += c1 < h0 ? 1UL : 0UL;
                c0 += l1;
                h1 += c0 < l1 ? 1UL : 0UL;
                c1 += h1;
                c2 += c1 < h1 ? 1UL : 0UL;
            }
            if (j <= jhi)
            {
                ulong h0 = BigMul(Unsafe.Subtract(ref ak, j), Unsafe.Add(ref bp, j), out ulong l0);
                c0 += l0;
                h0 += c0 < l0 ? 1UL : 0UL;
                c1 += h0;
                c2 += c1 < h0 ? 1UL : 0UL;
            }
            Unsafe.Add(ref rp, k) = c0;
            c0 = c1;
            c1 = c2;
            c2 = 0;
        }
        Unsafe.Add(ref rp, rn - 1) = c0;
        Debug.Assert(c1 == 0);
    }

    /// <summary>Schoolbook squaring: per column, the off-diagonal sum is doubled once and the diagonal square added.</summary>
    public static void SquareSchoolbook(Span<ulong> r, ReadOnlySpan<ulong> a, int n)
    {
        Debug.Assert(n >= 1 && r.Length >= 2 * n && a.Length >= n);
        if (Isa.Hw256 && n >= SqrVecMin && n < 1 << 13)
        {
            SquareSchoolbookVec256(r, a, n);
            return;
        }
        ref ulong rp = ref MemoryMarshal.GetReference(r);
        ref ulong ap = ref MemoryMarshal.GetReference(a);
        if (n == 1)
        {
            rp = BigMul(ap, ap, out ulong lo0);
            Unsafe.Add(ref rp, 1) = rp;
            rp = lo0;
            return;
        }
        ulong c0 = 0, c1 = 0, c2 = 0;
        nint len = n, rn = 2 * len;
        for (nint k = 0; k < rn - 1; k++)
        {
            // Off-diagonal pairs (i, k-i) with i < k-i.
            nint ilo = k - len + 1;
            if (ilo < 0)
            {
                ilo = 0;
            }
            nint ihi = (k - 1) >> 1;
            ulong d0 = 0, d1 = 0, d2 = 0;
            ref ulong ak = ref Unsafe.Add(ref ap, k);
            for (nint i = ilo; i <= ihi; i++)
            {
                ulong h = BigMul(Unsafe.Add(ref ap, i), Unsafe.Subtract(ref ak, i), out ulong l);
                d0 += l;
                h += d0 < l ? 1UL : 0UL;
                d1 += h;
                d2 += d1 < h ? 1UL : 0UL;
            }
            // Double the off-diagonal sum.
            d2 = (d2 << 1) | (d1 >> 63);
            d1 = (d1 << 1) | (d0 >> 63);
            d0 <<= 1;
            if ((k & 1) == 0)
            {
                ulong x = Unsafe.Add(ref ap, k >> 1);
                ulong h = BigMul(x, x, out ulong l);
                d0 += l;
                h += d0 < l ? 1UL : 0UL;
                d1 += h;
                d2 += d1 < h ? 1UL : 0UL;
            }
            c0 += d0;
            ulong cy = c0 < d0 ? 1UL : 0UL;
            d1 += cy;
            d2 += d1 < cy ? 1UL : 0UL;
            c1 += d1;
            c2 += d2 + (c1 < d1 ? 1UL : 0UL);
            Unsafe.Add(ref rp, k) = c0;
            c0 = c1;
            c1 = c2;
            c2 = 0;
        }
        Unsafe.Add(ref rp, rn - 1) = c0;
    }
}
