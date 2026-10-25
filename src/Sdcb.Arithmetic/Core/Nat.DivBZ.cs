using System.Buffers;
using System.Diagnostics;

namespace Sdcb.Arithmetic.Core;

internal static partial class Nat
{
    /// <summary>Divisor size (limbs) at which Burnikel-Ziegler divide-and-conquer division kicks in.</summary>
    public static readonly int BZThreshold = EnvInt("SDCB_ARITHMETIC_BZ") ?? (Isa.Hw256 ? 20 : 40);

    /// <summary>
    /// Burnikel-Ziegler divide-and-conquer division.
    /// q gets an-dn+1 limbs, r gets dn limbs; requires an &gt;= dn &gt;= BZThreshold.
    /// </summary>
    public static void DivRemBurnikelZiegler(Span<ulong> q, Span<ulong> r, ReadOnlySpan<ulong> a, int an, ReadOnlySpan<ulong> d, int dn)
        => DivBurnikelZieglerCore(q, r, a, an, d, dn, approx: false);

    /// <summary>
    /// Quotient-only Burnikel-Ziegler division: q (an-dn+1 limbs) gets q̂ with
    /// floor(a/d) &lt;= q̂ &lt;= floor(a/d) + 2. The last half-block skips its q*d0 product and
    /// correction, which only serve the remainder.
    /// </summary>
    public static void DivApproxBurnikelZiegler(Span<ulong> q, ReadOnlySpan<ulong> a, int an, ReadOnlySpan<ulong> d, int dn)
        => DivBurnikelZieglerCore(q, default, a, an, d, dn, approx: true);

    private static void DivBurnikelZieglerCore(Span<ulong> q, Span<ulong> r, ReadOnlySpan<ulong> a, int an, ReadOnlySpan<ulong> d, int dn, bool approx)
    {
        Debug.Assert(an >= dn && dn >= 2 && d[dn - 1] != 0);

        // Choose the block unit: m = 2^k such that j = ceil(dn/m) lands near the threshold.
        int m = 1;
        while ((dn + m - 1) / m > 2 * BZThreshold)
        {
            m <<= 1;
        }
        int j = (dn + m - 1) / m;
        int nPrime = j * m;

        long dBits = BitLength(d, dn);
        long sigma = (long)nPrime * 64 - dBits;
        int limbShift = (int)(sigma >> 6);
        int bitShift = (int)(sigma & 63);

        long aBits = BitLength(a, an);
        // Significant limbs of the shifted dividend. Quotient limbs above the last whole
        // block come from one short schoolbook step on the top nPrime+extra limbs when
        // extra is small; otherwise the top is zero-padded to one more whole block.
        // (Whole blocks only would waste an entire block division when the dividend is
        // just over a block boundary, e.g. 2n+1 limbs as in correctly rounded division.)
        int un = (int)((aBits + sigma + 63) / 64);
        int k = (un - nPrime) / nPrime;
        int extra = un - nPrime - k * nPrime;
        bool shortTop = extra <= nPrime / 8 + 2;
        int uLen = shortTop ? un + 1 : (k + 2) * nPrime;
        int qLen = uLen;
        int scratchLen = 2 * nPrime + 64;
        ulong[] pool = ArrayPool<ulong>.Shared.Rent(nPrime + uLen + qLen + scratchLen);
        Span<ulong> dv = pool.AsSpan(0, nPrime);
        Span<ulong> u = pool.AsSpan(nPrime, uLen);
        Span<ulong> qBig = pool.AsSpan(nPrime + uLen, qLen);
        Span<ulong> scratch = pool.AsSpan(nPrime + uLen + qLen, scratchLen);

        // dv = d << sigma (exactly nPrime limbs, top bit set).
        dv.Clear();
        if (bitShift == 0)
        {
            d.Slice(0, dn).CopyTo(dv.Slice(limbShift));
        }
        else
        {
            ulong outBits = ShiftLeft(dv.Slice(limbShift, dn), d.Slice(0, dn), dn, bitShift);
            if (limbShift + dn < nPrime)
            {
                dv[limbShift + dn] = outBits;
            }
            else
            {
                Debug.Assert(outBits == 0);
            }
        }
        Debug.Assert((dv[nPrime - 1] >> 63) != 0);

        // u = a << sigma (zero-padded to t blocks).
        u.Clear();
        if (bitShift == 0)
        {
            a.Slice(0, an).CopyTo(u.Slice(limbShift));
        }
        else
        {
            ulong outBits = ShiftLeft(u.Slice(limbShift, an), a.Slice(0, an), an, bitShift);
            if (limbShift + an < uLen)
            {
                u[limbShift + an] = outBits;
            }
            else
            {
                Debug.Assert(outBits == 0);
            }
        }

        qBig.Clear();
        if (shortTop)
        {
            DivRemNormalized(qBig.Slice(k * nPrime), u.Slice(k * nPrime, nPrime + extra + 1), nPrime + extra + 1, dv, nPrime);
        }
        else
        {
            // Top block holds extra < nPrime significant limbs, so it is below dv.
            k++;
        }
        for (int i = k - 1; i >= 0; i--)
        {
            if (approx && i == 0)
            {
                DivTwoByOneApprox(qBig.Slice(0, nPrime), u.Slice(0, 2 * nPrime), dv, nPrime, scratch);
            }
            else
            {
                DivTwoByOne(qBig.Slice(i * nPrime, nPrime), u.Slice(i * nPrime, 2 * nPrime), dv, nPrime, scratch);
            }
        }

        // Quotient is unaffected by the common shift. qn may exceed the block
        // area by one all-zero high limb (when a's top limb is small).
        int qn = an - dn + 1;
        int copyLen = Math.Min(qn, qLen);
        Debug.Assert(TrimmedLength(qBig, qLen) <= qn);
        qBig.Slice(0, copyLen).CopyTo(q);
        if (copyLen < qn)
        {
            q.Slice(copyLen, qn - copyLen).Clear();
        }

        if (approx)
        {
            ArrayPool<ulong>.Shared.Return(pool);
            return;
        }

        // Ieee754Remainder = (u[0..nPrime) >> sigma).
        if (bitShift != 0)
        {
            ShiftRight(u.Slice(0, nPrime), u.Slice(0, nPrime), nPrime, bitShift);
        }
        u.Slice(limbShift, dn).CopyTo(r);

        ArrayPool<ulong>.Shared.Return(pool);
    }

    /// <summary>
    /// Divides 2n limbs (in u, destroyed) by n normalized limbs; quotient (n limbs)
    /// into q, remainder left in u[0..n) with u[n..2n) cleared.
    /// Requires the high n limbs of u to be &lt; d.
    /// </summary>
    private static void DivTwoByOne(Span<ulong> q, Span<ulong> u, ReadOnlySpan<ulong> d, int n, Span<ulong> scratch)
    {
        if ((n & 1) != 0 || n <= BZThreshold)
        {
            DivRemNormalized(q, u, 2 * n, d.Slice(0, n), n);
            return;
        }
        int h = n >> 1;
        DivThreeByTwo(q.Slice(h, h), u.Slice(h, 3 * h), d, h, scratch);
        DivThreeByTwo(q.Slice(0, h), u.Slice(0, 3 * h), d, h, scratch);
    }

    /// <summary>
    /// <see cref="DivTwoByOne"/> without the remainder: the high half-quotient is exact, the low one is
    /// the (a2,a1)/d1 estimate, which exceeds the true digit by at most 2. u is destroyed.
    /// </summary>
    private static void DivTwoByOneApprox(Span<ulong> q, Span<ulong> u, ReadOnlySpan<ulong> d, int n, Span<ulong> scratch)
    {
        if ((n & 1) != 0 || n <= BZThreshold)
        {
            DivRemNormalized(q, u, 2 * n, d.Slice(0, n), n);
            return;
        }
        int h = n >> 1;
        DivThreeByTwo(q.Slice(h, h), u.Slice(h, 3 * h), d, h, scratch);
        Span<ulong> ul = u.Slice(0, 3 * h);
        ReadOnlySpan<ulong> d1 = d.Slice(h, h);
        if (Compare(ul.Slice(2 * h, h), d1, h) < 0)
        {
            DivTwoByOne(q.Slice(0, h), ul.Slice(h, 2 * h), d1, h, scratch);
        }
        else
        {
            q.Slice(0, h).Fill(ulong.MaxValue);
        }
    }

    /// <summary>
    /// Divides 3n limbs (a2,a1,a0) in u by the 2n-limb normalized divisor d = (d1,d0);
    /// quotient (n limbs) into q, remainder (2n limbs) left in u[0..2n) with u[2n..3n) cleared.
    /// Requires the high 2n limbs of u to be &lt; d.
    /// </summary>
    private static void DivThreeByTwo(Span<ulong> q, Span<ulong> u, ReadOnlySpan<ulong> d, int n, Span<ulong> scratch)
    {
        ReadOnlySpan<ulong> d1 = d.Slice(n, n);
        ReadOnlySpan<ulong> d0 = d.Slice(0, n);
        ulong c = 0; // virtual limb of the partial remainder above u[2n)

        if (Compare(u.Slice(2 * n, n), d1, n) < 0)
        {
            // (a2,a1) / d1 -> q, remainder into u[n..2n), u[2n..3n) cleared.
            DivTwoByOne(q, u.Slice(n, 2 * n), d1, n, scratch);
        }
        else
        {
            // a2 == d1 (since (a2,a1) < (d1,d0)); q = B^n - 1,
            // R1 = (a2,a1) - q*d1 = a1 + d1, carried in (u[n..2n), c).
            q.Fill(ulong.MaxValue);
            c = Add(u.Slice(n, n), u.Slice(n, n), d1, n);
            u.Slice(2 * n, n).Clear();
        }

        // u[0..2n) = (R1, a0) - q*d0
        Span<ulong> prod = scratch.Slice(0, 2 * n);
        Mul(prod, q, n, d0, n);
        ulong borrow = Sub(u.Slice(0, 2 * n), u.Slice(0, 2 * n), prod, 2 * n);

        long top = (long)c - (long)borrow;
        while (top < 0)
        {
            // Quotient was too large: q -= 1, remainder += d. (At most 2 iterations.)
            SubScalarInPlace(q, n, 1);
            top += (long)Add(u.Slice(0, 2 * n), u.Slice(0, 2 * n), d, 2 * n);
        }
        Debug.Assert(top == 0);
    }
}
