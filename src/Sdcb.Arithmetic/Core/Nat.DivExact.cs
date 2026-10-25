using System.Buffers;
using System.Diagnostics;
using System.Numerics;

namespace Sdcb.Arithmetic.Core;

internal static partial class Nat
{
    /// <summary>Quotient size (limbs) from which exact division uses a Newton-Hensel inverse.</summary>
    private static readonly int s_divexactHensel = EnvInt("SDCB_ARITHMETIC_DIVEXACT_H") ?? 20000;
    private static int DivExactHenselThreshold => s_divexactHensel;

    /// <summary>
    /// Exact division: q[0..an-dn+1) = a / d, valid only when d divides a.
    /// Works from the low end (Hensel/Jebelean), so no quotient estimation or remainder
    /// is needed and only the limbs below the quotient length are ever updated.
    /// </summary>
    public static void DivExact(Span<ulong> q, ReadOnlySpan<ulong> a, int an, ReadOnlySpan<ulong> d, int dn)
    {
        Debug.Assert(an >= dn && dn >= 1 && d[dn - 1] != 0);
        int qn = an - dn + 1;

        // Strip common low zero limbs, then the low zero bits of d.
        int z = 0;
        while (d[z] == 0)
        {
            z++;
        }
        a = a.Slice(z, an - z);
        d = d.Slice(z, dn - z);
        an -= z;
        dn -= z;
        int tz = BitOperations.TrailingZeroCount(d[0]);

        ulong[] pool = ArrayPool<ulong>.Shared.Rent(an + dn + 2);
        Span<ulong> w = pool.AsSpan(0, an + 1);
        Span<ulong> dd = pool.AsSpan(an + 1, dn);
        if (tz == 0)
        {
            a.Slice(0, an).CopyTo(w);
            d.Slice(0, dn).CopyTo(dd);
        }
        else
        {
            ShiftRight(w, a, an, tz);
            ShiftRight(dd, d, dn, tz);
        }
        int ddn = TrimmedLength(dd, dn);
        Span<ulong> qo = q.Slice(0, qn);
        qo.Clear();
        // After the shift the quotient still has at most qn limbs (and at least an - dn).
        if (qn >= DivExactHenselThreshold && ddn >= 32)
        {
            DivExactHensel(qo, w, qn, dd, ddn);
        }
        else if (qn > s_divexactDc && ddn > 4)
        {
            int need = 14 * qn + 256;
            ulong[] sc = ArrayPool<ulong>.Shared.Rent(need);
            DivExactDc(qo, w, qn, dd, ddn, sc.AsSpan(0, need));
            ArrayPool<ulong>.Shared.Return(sc);
        }
        else
        {
            DivExactBasecase(qo, w, qn, dd, ddn);
        }
        ArrayPool<ulong>.Shared.Return(pool);
    }

    /// <summary>Quotient size (limbs) below which divide-and-conquer exact division uses the basecase.</summary>
    private static readonly int s_divexactDc = EnvInt("SDCB_ARITHMETIC_DIVEXACT_DC") ?? 24;

    /// <summary>
    /// Divide-and-conquer Hensel division: the low half of the quotient depends only on the low
    /// limbs, so it is found recursively, its contribution is removed from the upper limbs with
    /// one fast multiplication, and the high half is found recursively from the rest.
    /// w (&gt;= qn limbs) is destroyed; only w[0..qn) and d[0..min(dn, qn)) are used.
    /// </summary>
    private static void DivExactDc(Span<ulong> q, Span<ulong> w, int qn, ReadOnlySpan<ulong> d, int dn, Span<ulong> scratch)
    {
        dn = Math.Min(dn, qn);
        if (qn <= s_divexactDc || dn <= 4)
        {
            DivExactBasecase(q, w, qn, d, dn);
            return;
        }
        int lo = qn >> 1, hi = qn - lo;
        DivExactDc(q, w, lo, d, dn, scratch);
        int qln = TrimmedLength(q, lo);
        if (qln > 0)
        {
            // w[lo..qn) -= (q_lo * d)[lo..qn); borrows past qn are irrelevant.
            int dl = Math.Min(dn, qn);
            int pn = qln + dl;
            Span<ulong> p = scratch.Slice(0, pn);
            Span<ulong> rest = scratch.Slice(pn);
            if (qln >= dl)
            {
                MulRecChecked(p, q, qln, d, dl, rest);
            }
            else
            {
                MulRecChecked(p, d, dl, q, qln, rest);
            }
            int top = Math.Min(pn, qn);
            if (top > lo)
            {
                SubInPlace(w.Slice(lo, hi), hi, p.Slice(lo, top - lo), top - lo);
            }
        }
        DivExactDc(q.Slice(lo), w.Slice(lo), hi, d, dn, scratch);
    }

    private static void MulRecChecked(Span<ulong> r, ReadOnlySpan<ulong> a, int an, ReadOnlySpan<ulong> b, int bn, Span<ulong> scratch)
    {
        Debug.Assert(scratch.Length >= MulScratchSize(an, bn));
        MulRecursive(r, a, an, b, bn, scratch);
    }

    /// <summary>Jebelean's right-to-left exact division; w (&gt;= qn limbs) is destroyed.</summary>
    private static void DivExactBasecase(Span<ulong> q, Span<ulong> w, int qn, ReadOnlySpan<ulong> d, int dn)
    {
        ulong dinv = MontgomeryContext.Inv64(d[0]);
        for (int i = 0; i < qn; i++)
        {
            ulong qi = w[i] * dinv;
            q[i] = qi;
            int len = Math.Min(dn, qn - i);
            ulong borrow = SubMulScalar(w.Slice(i, len), d, len, qi);
            int next = i + len;
            if (borrow != 0 && next < qn)
            {
                SubScalarInPlace(w.Slice(next), qn - next, borrow);
            }
        }
    }

    /// <summary>q = (w mod B^qn) * d^-1 mod B^qn with the 2-adic inverse from Newton iteration.</summary>
    private static void DivExactHensel(Span<ulong> q, ReadOnlySpan<ulong> w, int qn, ReadOnlySpan<ulong> d, int dn)
    {
        ulong[] pool = ArrayPool<ulong>.Shared.Rent(qn + 4 * qn + 4);
        Span<ulong> v = pool.AsSpan(0, qn);
        Span<ulong> tmp = pool.AsSpan(qn, 4 * qn + 4);
        HenselInverse(v, d, dn, qn, tmp);
        MulLow(q, w, v, qn, tmp);
        ArrayPool<ulong>.Shared.Return(pool);
    }

    /// <summary>v[0..k) = d^-1 mod B^k for odd d (Newton: v' = v + v(1 - d v)).</summary>
    private static void HenselInverse(Span<ulong> v, ReadOnlySpan<ulong> d, int dn, int k, Span<ulong> tmp)
    {
        v.Slice(0, k).Clear();
        v[0] = MontgomeryContext.Inv64(d[0]);
        int cur = 1;
        while (cur < k)
        {
            int nxt = Math.Min(2 * cur, k);
            // e = d*v mod B^nxt; its low cur limbs are 1, 0, ..., 0.
            Span<ulong> e = tmp.Slice(0, nxt);
            MulLow(e, d.Slice(0, Math.Min(dn, nxt)), v.Slice(0, cur), nxt, tmp.Slice(nxt));
            int h = nxt - cur;
            // v[cur..nxt) = -(v * e_hi) mod B^h
            Span<ulong> s = tmp.Slice(nxt, h);
            MulLow(s, v.Slice(0, Math.Min(cur, h)), e.Slice(cur, h), h, tmp.Slice(nxt + h));
            Span<ulong> hi = v.Slice(cur, h);
            Not(hi, s, h);
            AddScalarInPlace(hi, h, 1);
            cur = nxt;
        }
    }

    /// <summary>r[0..n) = (a * b) mod B^n for a, b of up to n limbs (full product, truncated).</summary>
    private static void MulLow(Span<ulong> r, ReadOnlySpan<ulong> a, ReadOnlySpan<ulong> b, int n, Span<ulong> scratch)
    {
        int an = TrimmedLength(a, Math.Min(a.Length, n));
        int bn = TrimmedLength(b, Math.Min(b.Length, n));
        if (an == 0 || bn == 0)
        {
            r.Slice(0, n).Clear();
            return;
        }
        ulong[]? big = null;
        Span<ulong> p = an + bn <= scratch.Length ? scratch.Slice(0, an + bn) : (big = ArrayPool<ulong>.Shared.Rent(an + bn)).AsSpan(0, an + bn);
        MulInto(p, a, an, b, bn);
        int c = Math.Min(n, an + bn);
        p.Slice(0, c).CopyTo(r);
        r.Slice(c, n - c).Clear();
        if (big is not null)
        {
            ArrayPool<ulong>.Shared.Return(big);
        }
    }
}
