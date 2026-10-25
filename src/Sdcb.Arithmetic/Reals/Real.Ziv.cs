using System.Diagnostics;

namespace Sdcb.Arithmetic;

public readonly partial struct Real
{
    /// <summary>
    /// An approximation step: computes the value at working precision wp, returning
    /// a signed fixed-point/scaled mantissa m with value = m * 2^e and an error bound
    /// of errUlps units in the last place of m. m == 0 with errUlps &lt; 0 signals "exact zero".
    /// </summary>
    internal delegate (BigInt m, long e, long errUlps) ZivStep(int wp);

    /// <summary>
    /// Ziv's onion-peeling loop: retries with growing working precision until the
    /// approximation rounds unambiguously to prec bits.
    /// </summary>
    internal static Real Ziv(int prec, RealRounding rnd, ZivStep step)
    {
        for (int wp = prec + 40; ; wp += Math.Max(64, wp >> 1))
        {
            (BigInt m, long e, long errUlps) = step(wp);
            if (m.IsZero)
            {
                if (errUlps < 0)
                {
                    return MakeZero(negative: false, prec);
                }
                continue; // too small to see at this precision; grow
            }
            bool neg = m.Sign < 0;
            BigInt mag = m.Abs();
            if (CanRound(mag, errUlps, prec, rnd))
            {
                return Round(mag, e, sticky: true, neg, prec, rnd, out _);
            }
            if (wp > 32L * prec + (1 << 16))
            {
                // Defensive cap: reachable only for exactly-representable results
                // that slipped past the algebraic-exactness checks.
                return Round(mag, e, sticky: true, neg, prec, rnd, out _);
            }
        }
    }

    /// <summary>
    /// True when a magnitude m known to within ±errUlps ulps rounds unambiguously
    /// to prec bits in all rounding modes (conservative version of mpfr_can_round).
    /// </summary>
    internal static bool CanRound(BigInt m, long errUlps, int prec, RealRounding rnd)
    {
        long bits = m.GetBitLength();
        long drop = bits - prec - 1; // keep 1 spare bit against msb changes
        if (drop < 2)
        {
            return false;
        }
        if (errUlps is >= 0 and < 1L << 62)
        {
            return CanRoundLimbs(m.Limbs, drop, (ulong)errUlps + 1, rnd);
        }
        BigInt err = new(errUlps + 1);
        BigInt window = BigInt.One << (int)Math.Min(drop, int.MaxValue);
        BigInt d = m & (window - BigInt.One);
        // Ambiguous when the dropped field is within err of 0, half, or the top.
        if (d.CompareTo(err) < 0 || (window - d).CompareTo(err) <= 0)
        {
            return false;
        }
        if (rnd == RealRounding.NearestEven)
        {
            BigInt half = window >> 1;
            BigInt diff = (d - half).Abs();
            if (diff.CompareTo(err) <= 0)
            {
                return false;
            }
        }
        return true;
    }

    // ---------------------------------------------------------------
    // Fixed-point helpers: F represents value * 2^wp as a signed BigInt.
    // ---------------------------------------------------------------

    /// <summary>Fixed-point multiply: (a*b) >> wp, truncated toward zero.</summary>
    internal static BigInt FixMul(BigInt a, BigInt b, int wp)
        => FixMulCore(a, b, wp);

    /// <summary>Fixed-point square.</summary>
    internal static BigInt FixSqr(BigInt a, int wp)
        => FixMulCore(a, a, wp);

    /// <summary>Fixed-point divide: (a &lt;&lt; wp) / b, truncated.</summary>
    internal static BigInt FixDiv(BigInt a, BigInt b, int wp)
        => FixDivCore(a, b, wp);

    /// <summary>Fixed-point inverse square root input helper: isqrt(a &lt;&lt; wp).</summary>
    internal static BigInt FixSqrt(BigInt a, int wp)
        => FixSqrtCore(a, wp);

    /// <summary>
    /// Evaluates S = Σ_{i&gt;=0} x^i / Π_{t=1..i} d(t) in fixed point (wp fractional bits) by
    /// rectangular splitting (Smith's method): with the powers x^1..x^k precomputed, each
    /// block of k terms costs one full multiplication plus k single-limb divisions, so the
    /// whole series takes about 2√N multiplications instead of N. |x| must be below 1 and
    /// the terms must decrease. The truncation and rounding error is a few ulps per block.
    /// </summary>
    internal static BigInt SeriesRect(BigInt x, int wp, Func<long, ulong> d, out long errUlps)
        => SeriesRect(x, wp, null, d, out errUlps);

    /// <summary>
    /// As above for S = Σ x^i Π_{t=1..i} num(t)/d(t) (num(t) &lt;= d(t), both word-sized);
    /// the Horner step multiplies by num(t) before the single-limb division.
    /// </summary>
    internal static BigInt SeriesRect(BigInt x, int wp, Func<long, ulong>? num, Func<long, ulong> d, out long errUlps)
    {
        BigInt one = BigInt.One << wp;
        if (x.IsZero)
        {
            errUlps = 1;
            return one;
        }
        // Number of terms: smallest N with |x|^N Π num/d < 2^-(wp+2), estimated in log2.
        double lx = (double)x.Abs().GetBitLength() - wp; // log2|x| upper bound
        double acc = 0;
        long n = 0;
        while (acc > -(wp + 2) && n < 1L << 30)
        {
            n++;
            acc += lx - Math.Log2(d(n)) + (num is null ? 0 : Math.Log2(num(n)));
        }
        int k = Math.Max(2, (int)Math.Ceiling(Math.Sqrt(n)));
        long m = (n + k - 1) / k;
        var pw = new BigInt[k + 1];
        pw[0] = one;
        pw[1] = x;
        for (int l = 2; l <= k; l++)
        {
            pw[l] = (l & 1) == 0 ? FixSqr(pw[l / 2], wp) : FixMul(pw[l - 1], x, wp);
        }
        BigInt s = BigInt.Zero;
        using var u = new FixAccumulator((wp >> 6) + 4);
        for (long j = m - 1; j >= 0; j--)
        {
            u.Set(j == m - 1 ? BigInt.Zero : FixMul(pw[k], s, wp));
            for (int l = k; l >= 1; l--)
            {
                long t = j * k + l;
                if (num is not null && !u.IsZero)
                {
                    u.MulSmall(num(t));
                }
                u.DivSmall(d(t));
                u.Add(pw[l - 1]);
            }
            s = u.ToBigInt();
        }
        errUlps = 4 * (m + k) + 8;
        return s;
    }

    /// <summary>Converts this Real to fixed point at wp fractional bits (truncated).</summary>
    internal BigInt ToFixed(int wp)
    {
        Debug.Assert(_kind == Kind.Normal);
        long shift = _exp + wp;
        BigInt mag = shift >= 0 ? _man.ShiftLeft(shift) : _man.ShiftRightTruncate(-shift);
        return IsNeg ? -mag : mag;
    }

    /// <summary>Builds a Real from a fixed-point value (round at prec bits).</summary>
    internal static Real FromFixed(BigInt f, int wp, int prec, RealRounding rnd)
    {
        if (f.IsZero)
        {
            return MakeZero(negative: false, prec);
        }
        return Round(f.Abs(), -wp, sticky: false, f.Sign < 0, prec, rnd, out _);
    }

    /// <summary>
    /// Approximation "f(x) = x + delta" for tiny x where |delta| &lt;&lt; ulp(x) at wp bits:
    /// returns x nudged by one sub-ulp in the direction of deltaSign (a value-space sign).
    /// </summary>
    private static (BigInt, long, long) NearX(Real x, int wp, int deltaSign)
    {
        Debug.Assert(x._kind == Kind.Normal);
        int pad = (int)Math.Clamp(wp - x.Precision + 16, 16, 1 << 24);
        BigInt m = x._man.ShiftLeft(pad);
        if (x.IsNeg)
        {
            m = -m;
        }
        m += deltaSign > 0 ? BigInt.One : BigInt.MinusOne;
        return (m, x._exp - pad, 4);
    }
}
