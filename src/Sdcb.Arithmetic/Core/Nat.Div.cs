using System.Buffers;
using System.Diagnostics;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.Intrinsics.X86;

namespace Sdcb.Arithmetic.Core;

internal static partial class Nat
{
    // ---------------------------------------------------------------
    // 128/64 division primitives
    // ---------------------------------------------------------------

    /// <summary>
    /// Divides the 128-bit value (hi:lo) by d, requiring hi &lt; d.
    /// Uses the hardware DIV instruction when available, else a portable routine.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static ulong DivRem128(ulong hi, ulong lo, ulong d, out ulong rem)
    {
        Debug.Assert(hi < d);
        if (X86Base.X64.IsSupported)
        {
            (ulong q, ulong r) = X86Base.X64.DivRem(lo, hi, d);
            rem = r;
            return q;
        }
        return DivRem128Portable(hi, lo, d, out rem);
    }

    /// <summary>Portable 128/64 division (Knuth-style using 32-bit halves).</summary>
    internal static ulong DivRem128Portable(ulong n1, ulong n0, ulong d, out ulong rem)
    {
        int shift = BitOperations.LeadingZeroCount(d);
        if (shift > 0)
        {
            d <<= shift;
            n1 = (n1 << shift) | (n0 >> (64 - shift));
            n0 <<= shift;
        }

        ulong d1 = d >> 32, d0 = d & 0xFFFFFFFFUL;

        ulong q1 = n1 / d1;
        ulong r1 = n1 - q1 * d1;
        ulong m = q1 * d0;
        r1 = (r1 << 32) | (n0 >> 32);
        if (r1 < m)
        {
            q1--; r1 += d;
            if (r1 >= d && r1 < m)
            {
                q1--; r1 += d;
            }
        }
        r1 -= m;

        ulong q0 = r1 / d1;
        ulong r0 = r1 - q0 * d1;
        m = q0 * d0;
        r0 = (r0 << 32) | (n0 & 0xFFFFFFFFUL);
        if (r0 < m)
        {
            q0--; r0 += d;
            if (r0 >= d && r0 < m)
            {
                q0--; r0 += d;
            }
        }
        r0 -= m;

        rem = r0 >> shift;
        return (q1 << 32) | q0;
    }

    /// <summary>
    /// Computes floor((B^2 - 1) / d) - B for a normalized divisor d (bit 63 set).
    /// The reciprocal from Möller and Granlund, "Improved division by invariant integers" (2011).
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static ulong Reciprocal(ulong d)
    {
        Debug.Assert((d >> 63) != 0);
        return DivRem128(~d, ulong.MaxValue, d, out _);
    }

    /// <summary>
    /// Divides (n1:n0) by normalized d using the precomputed reciprocal dinv.
    /// Requires n1 &lt; d. Returns the quotient; the remainder goes to r.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static ulong DivRemByReciprocal(ulong n1, ulong n0, ulong d, ulong dinv, out ulong r)
    {
        Debug.Assert(n1 < d);
        ulong qh = Nat.BigMul(n1, dinv, out ulong ql);
        // (qh:ql) += (n1 + 1 : n0)
        ql += n0;
        qh += n1 + 1 + (ql < n0 ? 1UL : 0UL);
        ulong rr = n0 - qh * d;
        if (rr > ql)
        {
            qh--;
            rr += d;
        }
        if (rr >= d)
        {
            qh++;
            rr -= d;
        }
        r = rr;
        return qh;
    }

    // ---------------------------------------------------------------
    // Division by a single limb
    // ---------------------------------------------------------------

    /// <summary>q[0..n) = a[0..n) / d; returns the remainder. q may alias a.</summary>
    public static ulong DivRemScalar(Span<ulong> q, ReadOnlySpan<ulong> a, int n, ulong d)
    {
        Debug.Assert(n >= 1 && d != 0);
        if ((d & (d - 1)) == 0)
        {
            // Power of two.
            int bits = BitOperations.TrailingZeroCount(d);
            ulong rem = a[0] & (d - 1);
            if (bits == 0)
            {
                a.Slice(0, n).CopyTo(q);
            }
            else
            {
                ShiftRight(q, a, n, bits);
            }
            return rem;
        }

        int shift = BitOperations.LeadingZeroCount(d);
        ulong r;
        if (shift == 0)
        {
            ulong dinv = Reciprocal(d);
            r = 0;
            for (int i = n - 1; i >= 0; i--)
            {
                q[i] = DivRemByReciprocal(r, a[i], d, dinv, out r);
            }
            return r;
        }
        else
        {
            ulong dn = d << shift;
            ulong dinv = Reciprocal(dn);
            int tnc = 64 - shift;
            r = a[n - 1] >> tnc;
            for (int i = n - 1; i >= 1; i--)
            {
                ulong nl = (a[i] << shift) | (a[i - 1] >> tnc);
                q[i] = DivRemByReciprocal(r, nl, dn, dinv, out r);
            }
            q[0] = DivRemByReciprocal(r, a[0] << shift, dn, dinv, out r);
            return r >> shift;
        }
    }

    /// <summary>Returns a[0..n) mod d without storing the quotient.</summary>
    public static ulong RemScalar(ReadOnlySpan<ulong> a, int n, ulong d)
    {
        Debug.Assert(n >= 1 && d != 0);
        if ((d & (d - 1)) == 0)
        {
            return a[0] & (d - 1);
        }
        int shift = BitOperations.LeadingZeroCount(d);
        ulong r;
        if (shift == 0)
        {
            ulong dinv = Reciprocal(d);
            r = 0;
            for (int i = n - 1; i >= 0; i--)
            {
                DivRemByReciprocal(r, a[i], d, dinv, out r);
            }
            return r;
        }
        else
        {
            ulong dn = d << shift;
            ulong dinv = Reciprocal(dn);
            int tnc = 64 - shift;
            r = a[n - 1] >> tnc;
            for (int i = n - 1; i >= 1; i--)
            {
                ulong nl = (a[i] << shift) | (a[i - 1] >> tnc);
                DivRemByReciprocal(r, nl, dn, dinv, out r);
            }
            DivRemByReciprocal(r, a[0] << shift, dn, dinv, out r);
            return r >> shift;
        }
    }

    // ---------------------------------------------------------------
    // General division (Knuth algorithm D with reciprocal estimation)
    // ---------------------------------------------------------------

    /// <summary>
    /// Computes q = a / d and r = a mod d.
    /// q needs an-dn+1 limbs, r needs dn limbs; requires an &gt;= dn &gt;= 1 and d
    /// normalized (top limb nonzero). q/r must not overlap a/d.
    /// </summary>
    public static void DivRem(Span<ulong> q, Span<ulong> r, ReadOnlySpan<ulong> a, int an, ReadOnlySpan<ulong> d, int dn)
    {
        Debug.Assert(an >= dn && dn >= 1 && d[dn - 1] != 0);
        if (dn == 1)
        {
            r[0] = DivRemScalar(q, a, an, d[0]);
            return;
        }
        if (dn >= BZThreshold)
        {
            int qn = an - dn + 1;
            if (2 * qn + 4 < dn)
            {
                DivRemShortQuotient(q, r, a, an, d, dn);
            }
            else if (dn >= DivNewtonThreshold)
            {
                DivRemNewton(q, r, a, an, d, dn);
            }
            else
            {
                DivRemBurnikelZiegler(q, r, a, an, d, dn);
            }
            return;
        }

        int shift = BitOperations.LeadingZeroCount(d[dn - 1]);
        int un = an + 1;
        ulong[] pool = ArrayPool<ulong>.Shared.Rent(un + dn);
        Span<ulong> u = pool.AsSpan(0, un);
        Span<ulong> dv = pool.AsSpan(un, dn);

        if (shift == 0)
        {
            a.Slice(0, an).CopyTo(u);
            u[an] = 0;
            d.Slice(0, dn).CopyTo(dv);
        }
        else
        {
            u[an] = ShiftLeft(u.Slice(0, an), a.Slice(0, an), an, shift);
            ShiftLeft(dv, d.Slice(0, dn), dn, shift);
        }

        DivRemNormalized(q, u, un, dv, dn);

        if (shift == 0)
        {
            u.Slice(0, dn).CopyTo(r);
        }
        else
        {
            ShiftRight(r, u.Slice(0, dn), dn, shift);
        }
        ArrayPool<ulong>.Shared.Return(pool);
    }

    /// <summary>
    /// Quotient-only division: q (an-dn+1 limbs) gets q̂ with floor(a/d) &lt;= q̂ &lt;= floor(a/d) + 2
    /// (exact outside the Burnikel-Ziegler range). Same requirements as <see cref="DivRem"/>.
    /// </summary>
    public static void DivApprox(Span<ulong> q, ReadOnlySpan<ulong> a, int an, ReadOnlySpan<ulong> d, int dn)
    {
        int qn = an - dn + 1;
        if (dn >= BZThreshold && 2 * qn + 4 >= dn && dn < DivNewtonThreshold)
        {
            DivApproxBurnikelZiegler(q, a, an, d, dn);
            return;
        }
        ulong[] pool = ArrayPool<ulong>.Shared.Rent(dn);
        DivRem(q, pool.AsSpan(0, dn), a, an, d, dn);
        ArrayPool<ulong>.Shared.Return(pool);
    }

    /// <summary>
    /// Schoolbook division with normalized divisor (bit 63 of dv[dn-1] set).
    /// u has un = an+1 limbs and is replaced by the remainder in u[0..dn).
    /// q receives un-dn = an-dn+1 limbs.
    /// </summary>
    internal static void DivRemNormalized(Span<ulong> q, Span<ulong> u, int un, ReadOnlySpan<ulong> dv, int dn)
    {
        Debug.Assert((dv[dn - 1] >> 63) != 0 && un > dn);
        ulong dtop = dv[dn - 1];
        ulong dinv = Reciprocal(dtop);

        for (int j = un - 1 - dn; j >= 0; j--)
        {
            ulong n2 = u[j + dn];
            ulong qhat;
            if (n2 >= dtop)
            {
                // Can only be equal (invariant: remainder < divisor).
                qhat = ulong.MaxValue;
            }
            else
            {
                qhat = DivRemByReciprocal(n2, u[j + dn - 1], dtop, dinv, out _);
            }

            ulong borrow = SubMulScalar(u.Slice(j, dn), dv, dn, qhat);
            ulong t = n2 - borrow;
            bool neg = n2 < borrow;
            while (neg)
            {
                qhat--;
                ulong cy = Add(u.Slice(j, dn), u.Slice(j, dn), dv, dn);
                if (t == ulong.MaxValue && cy == 1)
                {
                    t = 0;
                    neg = false;
                }
                else
                {
                    t += cy;
                }
            }
            Debug.Assert(t == 0);
            u[j + dn] = 0;
            q[j] = qhat;
        }
    }
}
