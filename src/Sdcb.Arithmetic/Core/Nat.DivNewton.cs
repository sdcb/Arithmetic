using System.Buffers;
using System.Diagnostics;
using System.Numerics;

namespace Sdcb.Arithmetic.Core;

internal static partial class Nat
{
    /// <summary>Divisor size (limbs) from which division uses a Newton reciprocal instead of Burnikel-Ziegler.</summary>
    internal static int DivNewtonThreshold => s_divNewtonEnv ?? (Isa.Hw256 ? 12000 : Isa.Hw128 ? 16000 : int.MaxValue);

    private static readonly int? s_divNewtonEnv = EnvInt("SDCB_ARITHMETIC_DIV_NEWTON");

    /// <summary>Below this size the reciprocal is computed by a direct division.</summary>
    private const int RecipBaseCase = 48;

    /// <summary>
    /// Approximate reciprocal (Brent-Zimmermann, Modern Computer Arithmetic, Alg. 3.5):
    /// for a normalized a[0..n) (top bit set) computes x[0..n] with
    /// a*x &lt; B^(2n) &lt;= a*(x+2).
    /// </summary>
    internal static void Recip(Span<ulong> x, ReadOnlySpan<ulong> a, int n)
    {
        Debug.Assert(n >= 1 && (a[n - 1] >> 63) != 0 && x.Length >= n + 1);
        if (n <= RecipBaseCase)
        {
            ulong[] pool = ArrayPool<ulong>.Shared.Rent(4 * n + 1);
            Span<ulong> num = pool.AsSpan(0, 2 * n);
            Span<ulong> rem = pool.AsSpan(2 * n, n);
            Span<ulong> quo = pool.AsSpan(3 * n, n + 1);
            num.Fill(ulong.MaxValue);
            if (n >= BZThreshold)
            {
                DivRemBurnikelZiegler(quo, rem, num, 2 * n, a, n);
            }
            else
            {
                DivRem(quo, rem, num, 2 * n, a, n);
            }
            quo.CopyTo(x);
            ArrayPool<ulong>.Shared.Return(pool);
            return;
        }

        int l = (n - 1) >> 1, h = n - l;
        Span<ulong> xh = x.Slice(l, h + 1);
        Recip(xh, a.Slice(l, h), h);

        int tn = n + h + 1;
        ulong[] tp = ArrayPool<ulong>.Shared.Rent(tn + 2 * h + 2);
        Span<ulong> t = tp.AsSpan(0, tn);
        Span<ulong> u = tp.AsSpan(tn, 2 * h + 2);

        // T = A * X_h, pulled below B^(n+h).
        Mul(t, a, n, xh, h + 1);
        while (t[n + h] != 0)
        {
            SubScalarInPlace(xh, h + 1, 1);
            SubInPlace(t, tn, a, n);
        }
        // T = B^(n+h) - T  (<= 2 B^n by the recursive guarantee).
        Not(t, t, n + h);
        AddScalarInPlace(t, n + h, 1);
        int tmn = TrimmedLength(t.Slice(l, n + h - l), n + h - l);
        Debug.Assert(tmn <= h + 1);

        // X = X_h * B^l + floor(T_m * X_h / B^(2h - l)).
        x.Slice(0, l).Clear();
        if (tmn > 0)
        {
            Span<ulong> tm = t.Slice(l, tmn);
            if (tmn >= h + 1)
            {
                Mul(u, tm, tmn, xh, h + 1);
            }
            else
            {
                Mul(u, xh, h + 1, tm, tmn);
            }
            int un = tmn + h + 1;
            int shift = 2 * h - l;
            if (un > shift)
            {
                // xh aliases x[l..]; the correction is added into the whole of x.
                ulong cy = AddInPlace(x, n + 1, u.Slice(shift), un - shift);
                Debug.Assert(cy == 0);
            }
        }
        ArrayPool<ulong>.Shared.Return(tp);
    }

    /// <summary>
    /// Division with a Newton reciprocal: q[0..an-dn+1) = a / d, r[0..dn) = a mod d.
    /// Each dn-limb block of quotient costs two dn x dn products.
    /// </summary>
    internal static void DivRemNewton(Span<ulong> q, Span<ulong> r, ReadOnlySpan<ulong> a, int an, ReadOnlySpan<ulong> d, int dn)
    {
        Debug.Assert(an >= dn && dn >= 2 && d[dn - 1] != 0);
        int n = dn;
        int shift = BitOperations.LeadingZeroCount(d[n - 1]);
        long uBits = BitLength(a, an) + shift;
        int blocks = (int)((uBits + 64L * n - 1) / (64L * n));
        if (blocks < 1)
        {
            blocks = 1;
        }
        int uLen = Math.Max(blocks * n, an + 1);

        ulong[] pool = ArrayPool<ulong>.Shared.Rent(n + (n + 1) + uLen + uLen + (2 * n + 1) + 2 * n);
        int o = 0;
        Span<ulong> dv = pool.AsSpan(o, n); o += n;
        Span<ulong> xr = pool.AsSpan(o, n + 1); o += n + 1;
        Span<ulong> u = pool.AsSpan(o, uLen); o += uLen;
        Span<ulong> qBig = pool.AsSpan(o, uLen); o += uLen;
        Span<ulong> p = pool.AsSpan(o, 2 * n + 1); o += 2 * n + 1;
        Span<ulong> w = pool.AsSpan(o, 2 * n);

        if (shift == 0)
        {
            d.Slice(0, n).CopyTo(dv);
            a.Slice(0, an).CopyTo(u);
            u.Slice(an).Clear();
        }
        else
        {
            ShiftLeft(dv, d, n, shift);
            ulong top = ShiftLeft(u, a, an, shift);
            u.Slice(an).Clear();
            u[an] = top;
        }
        Recip(xr, dv, n);

        // w holds the running 2n-limb partial dividend (R, next block); R starts as the
        // top block, reduced below d by at most one subtraction.
        qBig.Clear();
        Span<ulong> rem = w.Slice(n, n);
        u.Slice((blocks - 1) * n, n).CopyTo(rem);
        if (Compare(rem, dv, n) >= 0)
        {
            Sub(rem, rem, dv, n);
            qBig[(blocks - 1) * n] = 1;
        }
        for (int c = blocks - 2; c >= 0; c--)
        {
            Span<ulong> qc = qBig.Slice(c * n, n);
            u.Slice(c * n, n).CopyTo(w);
            int rn = TrimmedLength(rem, n);
            if (rn == 0)
            {
                qc.Clear();
            }
            else
            {
                // qhat = floor(R * X / B^n) <= true quotient, short by at most 4.
                Mul(p, xr, n + 1, rem, rn);
                p.Slice(n + 1 + rn).Clear();
                Debug.Assert(p[2 * n] == 0);
                p.Slice(n, n).CopyTo(qc);
                int qn = TrimmedLength(qc, n);
                if (qn > 0)
                {
                    Span<ulong> prod = p; // R*X is no longer needed
                    MulInto(prod, qc, qn, dv, n);
                    ulong borrow = SubInPlace(w, 2 * n, prod, qn + n);
                    Debug.Assert(borrow == 0);
                }
            }
            // Fix up: while W >= d, W -= d, q++.
            while (!IsZero(w.Slice(n, n), n) || Compare(w, dv, n) >= 0)
            {
                SubInPlace(w, 2 * n, dv, n);
                AddScalarInPlace(qc, n, 1);
            }
            w.Slice(0, n).CopyTo(rem);
        }

        int qnOut = an - dn + 1;
        Debug.Assert(TrimmedLength(qBig, uLen) <= qnOut);
        qBig.Slice(0, qnOut).CopyTo(q);
        if (shift == 0)
        {
            rem.CopyTo(r);
        }
        else
        {
            ShiftRight(r, rem, n, shift);
        }
        ArrayPool<ulong>.Shared.Return(pool);
    }

    /// <summary>prod[0..an+bn) = a * b choosing operand order.</summary>
    private static void MulInto(Span<ulong> prod, ReadOnlySpan<ulong> a, int an, ReadOnlySpan<ulong> b, int bn)
    {
        if (an >= bn)
        {
            Mul(prod, a, an, b, bn);
        }
        else
        {
            Mul(prod, b, bn, a, an);
        }
    }

    /// <summary>
    /// Division with a short quotient (qn much smaller than dn): the quotient is estimated
    /// from the top qn+2 limbs of both operands, then corrected against the full divisor.
    /// </summary>
    internal static void DivRemShortQuotient(Span<ulong> q, Span<ulong> r, ReadOnlySpan<ulong> a, int an, ReadOnlySpan<ulong> d, int dn)
    {
        int qn = an - dn + 1;
        int k = qn + 2;              // divisor limbs kept
        int cut = dn - k;            // limbs dropped from both
        Debug.Assert(cut > 0);
        ulong[] pool = ArrayPool<ulong>.Shared.Rent((qn + 1) + k + (an + 1) + (an + 1));
        Span<ulong> qe = pool.AsSpan(0, qn + 1);
        Span<ulong> re = pool.AsSpan(qn + 1, k);
        Span<ulong> prod = pool.AsSpan(qn + 1 + k, an + 1);
        Span<ulong> w = pool.AsSpan(qn + 1 + k + an + 1, an + 1);

        // Estimate from the truncated operands: off by at most a couple of units either way.
        DivRem(qe, re, a.Slice(cut), an - cut, d.Slice(cut), k);
        qe[qn] = 0;
        int qen = TrimmedLength(qe, qn);

        a.Slice(0, an).CopyTo(w);
        w[an] = 0;
        if (qen > 0)
        {
            MulInto(prod, d, dn, qe, qen);
            int pn = TrimmedLength(prod, dn + qen);
            while (CompareNormalized(w, TrimmedLength(w, an + 1), prod, pn) < 0)
            {
                SubScalarInPlace(qe, qn, 1);
                SubInPlace(prod, pn, d, dn);
                pn = TrimmedLength(prod, pn);
            }
            SubInPlace(w, an + 1, prod, pn);
        }
        int wn = TrimmedLength(w, an + 1);
        while (CompareNormalized(w, wn, d, dn) >= 0)
        {
            SubInPlace(w, wn, d, dn);
            wn = TrimmedLength(w, wn);
            AddScalarInPlace(qe, qn + 1, 1);
        }
        Debug.Assert(qe[qn] == 0);
        qe.Slice(0, qn).CopyTo(q);
        r.Slice(0, dn).Clear();
        w.Slice(0, Math.Min(wn, dn)).CopyTo(r);
        ArrayPool<ulong>.Shared.Return(pool);
    }
}
