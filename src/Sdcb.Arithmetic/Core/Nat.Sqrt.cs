using System.Buffers;
using System.Diagnostics;
using System.Numerics;

namespace Sdcb.Arithmetic.Core;

internal static partial class Nat
{
    /// <summary>
    /// Square root of a normalized 2-limb value (np[1] &gt;= 2^62): sp[0] = floor(sqrt),
    /// np[0] = low limb of the remainder; returns the remainder's high limb (0 or 1).
    /// </summary>
    internal static ulong SqrtRem128(Span<ulong> sp, Span<ulong> np)
    {
        ulong hi = np[1], lo = np[0];
        Debug.Assert(hi >= 1UL << 62);
        UInt128 n = new(hi, lo);
        double d = Math.Sqrt(hi * 18446744073709551616.0 + lo);
        ulong s = d >= 18446744073709551615.0 ? ulong.MaxValue : (ulong)d;
        if (hi >= s)
        {
            s = ulong.MaxValue;
        }
        else
        {
            // One Newton step from a 52-bit-accurate estimate lands within 1 of the root.
            ulong q = DivRem128(hi, lo, s, out _);
            s = (ulong)(((UInt128)s + q) >> 1);
        }
        while ((UInt128)s * s > n)
        {
            s--;
        }
        while (s != ulong.MaxValue && (UInt128)(s + 1) * (s + 1) <= n)
        {
            s++;
        }
        UInt128 r = n - (UInt128)s * s;
        sp[0] = s;
        np[0] = (ulong)r;
        return (ulong)(r >> 64);
    }

    /// <summary>
    /// Zimmermann's Karatsuba square root (INRIA RR-3805, 1999): sp[0..n) = floor(sqrt({np, 2n})),
    /// np[0..n) = low n limbs of the remainder; returns its high limb. {np, 2n} must be
    /// normalized (np[2n-1] &gt;= 2^62). scratch needs n/2+1 limbs, tmp n limbs.
    /// </summary>
    private static ulong SqrtRemRecursive(Span<ulong> sp, Span<ulong> np, int n, Span<ulong> scratch, Span<ulong> tmp, bool approx = false)
    {
        Debug.Assert(n > 1 && np[2 * n - 1] >= 1UL << 62);
        int l = n / 2, h = n - l;
        ulong q = h == 1
            ? SqrtRem128(sp.Slice(l), np.Slice(2 * l))
            : SqrtRemRecursive(sp.Slice(l), np.Slice(2 * l), h, scratch, tmp);
        if (q != 0)
        {
            ulong bw = Sub(np.Slice(2 * l), np.Slice(2 * l), sp.Slice(l), h);
            Debug.Assert(bw == 1);
        }
        // (q', r) = divrem(np[l..l+n), s'); the true quotient is q'/2.
        np.Slice(l, n).CopyTo(tmp);
        DivRem(scratch, np.Slice(l, h), tmp, n, sp.Slice(l), h);
        q += scratch[l];
        long c = (long)(scratch[0] & 1);
        ShiftRight(sp, scratch, l, 1);
        sp[l - 1] |= q << 63;
        q >>= 1;
        if (approx)
        {
            // The root is (sp + q B^l) or one less; the remainder pass only decides which.
            if (q != 0 && AddScalar(sp.Slice(l), sp.Slice(l), h, q) != 0)
            {
                sp.Slice(0, n).Fill(ulong.MaxValue);
            }
            return 0;
        }
        if (c != 0)
        {
            c = (long)Add(np.Slice(l), np.Slice(l), sp.Slice(l), h);
        }
        Square(np.Slice(n), sp, l);
        long b = (long)q + (long)Sub(np, np, np.Slice(n), 2 * l);
        c -= l == h ? b : (long)SubScalar(np.Slice(2 * l), np.Slice(2 * l), 1, (ulong)b);
        if (c < 0)
        {
            q = AddScalar(sp.Slice(l), sp.Slice(l), h, q);
            c += (long)AddMulScalar(np, sp, n, 2) + 2 * (long)q;
            c -= (long)SubScalar(np, np, n, 1);
            q -= SubScalar(sp, sp, n, 1);
        }
        return (ulong)c;
    }

    /// <summary>
    /// sp[0..ceil(nn/2)) = floor(sqrt(np[0..nn))) and, when rp is non-empty, rp = remainder
    /// (needs ceil(nn/2)+1 limbs). np[nn-1] must be nonzero. Returns the remainder's size.
    /// </summary>
    internal static int SqrtRem(Span<ulong> sp, Span<ulong> rp, ReadOnlySpan<ulong> np, int nn)
    {
        Debug.Assert(nn >= 1 && np[nn - 1] != 0);
        ulong high = np[nn - 1];
        int c = (high & 0xC000000000000000UL) != 0 ? 0 : BitOperations.LeadingZeroCount(high) / 2;
        int tn = (nn + 1) / 2;
        int poolSize = 2 * tn + 2 + (tn / 2 + 1) + tn;
        ulong[] pool = ArrayPool<ulong>.Shared.Rent(poolSize);
        Span<ulong> tp = pool.AsSpan(0, 2 * tn + 2);
        Span<ulong> scratch = pool.AsSpan(2 * tn + 2, tn / 2 + 1);
        Span<ulong> tmp = pool.AsSpan(2 * tn + 2 + tn / 2 + 1, tn);
        tp.Clear();
        int rn;
        if ((nn & 1) != 0 || c > 0)
        {
            // Shift left by 2k bits so that {tp, 2tn} is normalized; k = c (+32 for odd nn).
            if (c != 0)
            {
                ShiftLeft(tp.Slice(nn & 1), np, nn, 2 * c);
            }
            else
            {
                np.Slice(0, nn).CopyTo(tp.Slice(nn & 1));
            }
            int k = c + ((nn & 1) != 0 ? 32 : 0);
            ulong mask = (1UL << k) - 1;
            ulong rl = tn == 1 ? SqrtRem128(sp, tp) : SqrtRemRecursive(sp, tp, tn, scratch, tmp);
            // 2^(2k) N = S^2 + R = (S - s0)^2 + 2 S s0 - s0^2 + R with s0 = S mod 2^k.
            ulong s0 = sp[0] & mask;
            rl += AddMulScalar(tp, sp, tn, 2 * s0);
            ulong hi = Math.BigMul(s0, s0, out ulong lo);
            ulong bw = SubWithBorrow(tp[0], lo, 0, out tp[0]);
            if (tn > 1)
            {
                bw = SubWithBorrow(tp[1], hi, bw, out tp[1]);
                rl -= tn > 2 ? SubScalarInPlace(tp.Slice(2), tn - 2, bw) : bw;
            }
            else
            {
                rl -= hi + bw;
            }
            ShiftRight(sp, sp, tn, k);
            tp[tn] = rl;
            if (!rp.IsEmpty)
            {
                int shift = 2 * k, start = 0, len = tn + 1;
                if (shift >= 64)
                {
                    start = 1;
                    len = tn;
                    shift -= 64;
                }
                if (shift != 0)
                {
                    ShiftRight(rp, tp.Slice(start), len, shift);
                }
                else
                {
                    tp.Slice(start, len).CopyTo(rp);
                }
                rn = TrimmedLength(rp, len);
            }
            else
            {
                rn = 0;
            }
        }
        else
        {
            np.Slice(0, nn).CopyTo(tp);
            ulong rh = tn == 1 ? SqrtRem128(sp, tp) : SqrtRemRecursive(sp, tp, tn, scratch, tmp);
            if (!rp.IsEmpty)
            {
                tp.Slice(0, tn).CopyTo(rp);
                rp[tn] = rh;
                rn = TrimmedLength(rp, tn + 1);
            }
            else
            {
                rn = 0;
            }
        }
        ArrayPool<ulong>.Shared.Return(pool);
        return rn;
    }

    /// <summary>
    /// sp[0..ceil(nn/2)) = S with floor(sqrt(N)) &lt;= S &lt;= floor(sqrt(N)) + 1 for N = np[0..nn)
    /// (np[nn-1] nonzero, nn &gt;= 2): the top Karatsuba step skips the remainder squaring.
    /// </summary>
    internal static void SqrtApprox(Span<ulong> sp, ReadOnlySpan<ulong> np, int nn)
    {
        Debug.Assert(nn >= 3 && np[nn - 1] != 0);
        ulong high = np[nn - 1];
        int c = (high & 0xC000000000000000UL) != 0 ? 0 : BitOperations.LeadingZeroCount(high) / 2;
        int tn = (nn + 1) / 2;
        int poolSize = 2 * tn + 2 + (tn / 2 + 1) + tn;
        ulong[] pool = ArrayPool<ulong>.Shared.Rent(poolSize);
        Span<ulong> tp = pool.AsSpan(0, 2 * tn + 2);
        Span<ulong> scratch = pool.AsSpan(2 * tn + 2, tn / 2 + 1);
        Span<ulong> tmp = pool.AsSpan(2 * tn + 2 + tn / 2 + 1, tn);
        tp.Clear();
        // TrimmedLength by 2^(2k) as in Sqrt; (S' >> k) stays within one of floor(sqrt(N)).
        int k = c + ((nn & 1) != 0 ? 32 : 0);
        if (c != 0)
        {
            ShiftLeft(tp.Slice(nn & 1), np, nn, 2 * c);
        }
        else
        {
            np.Slice(0, nn).CopyTo(tp.Slice(nn & 1));
        }
        SqrtRemRecursive(sp, tp, tn, scratch, tmp, approx: true);
        if (k != 0)
        {
            ShiftRight(sp, sp, tn, k);
        }
        ArrayPool<ulong>.Shared.Return(pool);
    }
}