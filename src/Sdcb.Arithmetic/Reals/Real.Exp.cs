using System.Diagnostics;

namespace Sdcb.Arithmetic;

public readonly partial struct Real
{
    // ---------------------------------------------------------------
    // Exponential family (mpfr_exp, exp2, exp10, expm1)
    // ---------------------------------------------------------------

    /// <summary>e^x correctly rounded (mpfr_exp).</summary>
    public static Real Exp(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = precision != 0 ? precision : x.Precision;
        switch (x._kind)
        {
            case Kind.NaN:
                return MakeNaN(prec);
            case Kind.Inf:
                return x.IsNeg ? MakeZero(negative: false, prec) : MakeInf(negative: false, prec);
            case Kind.Zero:
                return Round(BigInt.One, 0, sticky: false, negative: false, prec, rounding, out _);
        }
        if (x.BinaryExponent >= 62)
        {
            // |x| is astronomically large: the result exponent would not fit.
            return x.IsNeg ? MakeZero(negative: false, prec) : MakeInf(negative: false, prec);
        }
        return Ziv(prec, rounding, wp => ExpStep(x, wp));
    }

    /// <summary>Core exp approximation: returns exp(x) as (m, e, err).</summary>
    private static (BigInt, long, long) ExpStep(Real x, int wp)
    {
        // With rectangular splitting the series costs ~2√N multiplications, so fewer
        // halvings (each undone by a squaring) than the classic √wp pay off.
        int j = (int)(Math.Sqrt(wp) * ExpReductionFactor) + 1;
        int wpi = wp + 2 * j + 32;

        // Argument reduction: x = n*ln2 + r, |r| <= ~0.7
        long n = (long)Math.Round(x.ToDouble() / 0.6931471805599453);
        int nBits = 64 - System.Numerics.BitOperations.LeadingZeroCount((ulong)Math.Abs(n) | 1);
        BigInt ln2 = Ln2Fixed(wpi + nBits + 4);
        BigInt rf = x.ToFixed(wpi + nBits + 4) - new BigInt(n) * ln2;
        rf = rf.ShiftRightTruncate(nBits + 4);

        // Second reduction: r' = r / 2^j
        rf = rf.ShiftRightTruncate(j);

        // Taylor sum of e^(r') = Σ r'^i / i!
        BigInt sum = SeriesRect(rf, wpi, static t => (ulong)t, out long errS);
        // Undo the second reduction by squaring j times.
        for (int k = 0; k < j; k++)
        {
            sum = FixSqr(sum, wpi);
        }
        // Error: series + reduction ulps, doubled per squaring.
        long err = (errS + 16) << Math.Min(j, 40);
        return (sum, n - wpi, err);
    }

    internal static double ExpReductionFactor = double.TryParse(Environment.GetEnvironmentVariable("SDCB_ARITHMETIC_EXPRED"), out double f) ? f : 0.3;

    /// <summary>Signed truncating division by a small positive integer.</summary>
    private static BigInt DivRoundZ(BigInt v, ulong d) => BigInt.DivRem(v, d, out _);

    /// <summary>2^x (mpfr_exp2).</summary>
    public static Real Exp2(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
        => ExpScaled(x, 2, precision, rounding);

    /// <summary>10^x (mpfr_exp10).</summary>
    public static Real Exp10(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
        => ExpScaled(x, 10, precision, rounding);

    private static Real ExpScaled(Real x, int b, int precision, RealRounding rounding)
    {
        int prec = precision != 0 ? precision : x.Precision;
        switch (x._kind)
        {
            case Kind.NaN:
                return MakeNaN(prec);
            case Kind.Inf:
                return x.IsNeg ? MakeZero(negative: false, prec) : MakeInf(negative: false, prec);
            case Kind.Zero:
                return Round(BigInt.One, 0, sticky: false, negative: false, prec, rounding, out _);
        }
        // Exact power for integer x with base 2.
        if (b == 2 && x.IsInteger && x.BinaryExponent < 62)
        {
            long k = (long)x.ToBigInt();
            return Round(BigInt.One, k, sticky: false, negative: false, prec, rounding, out _);
        }
        if (b == 10 && x.IsInteger && x.BinaryExponent < 21)
        {
            // 10^k = 5^k * 2^k exactly (or a single correctly-rounded division for k < 0).
            long k = (long)x.ToBigInt();
            if (k >= 0)
            {
                return Round(BigInt.Pow(5, (ulong)k), k, sticky: false, negative: false, prec, rounding, out _);
            }
            BigInt den = BigInt.Pow(5, (ulong)(-k));
            long s = prec + 3 + den.GetBitLength();
            BigInt.DivRem(BigInt.One << (int)Math.Min(s, int.MaxValue), den, out BigInt q, out BigInt r);
            return Round(q, k - s, sticky: !r.IsZero, negative: false, prec, rounding, out _);
        }
        if (x.BinaryExponent >= 56)
        {
            return x.IsNeg ? MakeZero(negative: false, prec) : MakeInf(negative: false, prec);
        }
        return Ziv(prec, rounding, wp =>
        {
            // b^x = exp(x * ln b)
            int wpi = wp + 64;
            BigInt lnb = b == 2 ? Ln2Fixed(wpi) : LogOfSmallInt(b, wpi);
            BigInt xf = x.ToFixed(wpi);
            BigInt prodF = FixMul(xf, lnb, wpi); // x*ln b at wpi
            Real y = FromFixed(prodF, wpi, wpi, RealRounding.NearestEven);
            if (y._kind != Kind.Normal)
            {
                return (BigInt.One << wp, -wp, 8); // x*lnb ~ 0 → result ~ 1
            }
            (BigInt m, long e, long err) = ExpStep(y, wp);
            return (m, e, err + 256);
        });
    }

    /// <summary>ln(k) for a small integer k at wp fractional bits.</summary>
    private static BigInt LogOfSmallInt(int k, int wp)
    {
        Real r = FromBigInt(new BigInt(k), wp + 8);
        (BigInt m, long e, _) = LogStep(r, wp + 8);
        long shift = e + wp;
        return shift >= 0 ? m.ShiftLeft(shift) : m.ShiftRightTruncate(-shift);
    }

    /// <summary>e^x - 1 with full accuracy near zero (mpfr_expm1).</summary>
    public static Real ExpM1(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = precision != 0 ? precision : x.Precision;
        switch (x._kind)
        {
            case Kind.NaN:
                return MakeNaN(prec);
            case Kind.Inf:
                return x.IsNeg ? Round(BigInt.One, 0, sticky: false, negative: true, prec, rounding, out _) : MakeInf(negative: false, prec);
            case Kind.Zero:
                return MakeZero(x.IsNeg, prec);
        }
        if (x.BinaryExponent >= 62)
        {
            return x.IsNeg
                ? Round(BigInt.One, 0, sticky: true, negative: true, prec, rounding, out _)
                : MakeInf(negative: false, prec);
        }
        return Ziv(prec, rounding, wp =>
        {
            // For the cancellation case, work at wp + |exponent of x| extra bits.
            long extra = Math.Max(0, -x.BinaryExponent + 4);
            int wpi = wp + (int)Math.Min(extra, 1 << 24) + 32;
            if (x.IsNeg && x.ToDouble() < -(wpi + 16) * 0.6931)
            {
                // exp(x) underflows the working window: result = -1 + tiny.
                return (BigInt.One - (BigInt.One << wpi), -wpi, 2);
            }
            (BigInt m, long e, long err) = ExpStep(x, wpi);
            // m*2^e - 1; when e > 0 the -1 is below one ulp of m (folded into err).
            long sh = -e;
            BigInt one = sh >= 0 ? BigInt.One << (int)Math.Min(sh, int.MaxValue) : BigInt.Zero;
            BigInt diff = m - one;
            return (diff, e, err + 2);
        });
    }

    // ---------------------------------------------------------------
    // Logarithm family (mpfr_log, log2, log10, log1p)
    // ---------------------------------------------------------------

    /// <summary>Natural logarithm, correctly rounded (mpfr_log).</summary>
    public static Real Log(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = precision != 0 ? precision : x.Precision;
        Real? special = LogSpecial(x, prec);
        if (special.HasValue)
        {
            return special.Value;
        }
        // ln(1) == +0 exactly.
        if (IsPowerOfTwoValue(x, out long pk) && pk == 0)
        {
            return MakeZero(negative: false, prec);
        }
        return Ziv(prec, rounding, wp => LogStep(x, wp));
    }

    private static Real? LogSpecial(Real x, int prec)
    {
        switch (x._kind)
        {
            case Kind.NaN:
                return MakeNaN(prec);
            case Kind.Zero:
                return MakeInf(negative: true, prec); // log(±0) = -inf
            case Kind.Inf:
                return x.IsNeg ? MakeNaN(prec) : MakeInf(negative: false, prec);
        }
        if (x.IsNeg)
        {
            return MakeNaN(prec);
        }
        return null;
    }

    private static bool IsPowerOfTwoValue(Real x, out long k)
    {
        k = 0;
        if (x._kind != Kind.Normal || !x._man.IsPowerOfTwo)
        {
            return false;
        }
        k = x.BinaryExponent - 1;
        return true;
    }

    /// <summary>Working precision (bits) from which log uses the AGM instead of square roots + series.</summary>
    internal static int LogAgmThreshold = int.TryParse(Environment.GetEnvironmentVariable("SDCB_ARITHMETIC_LOGAGM"), out int t0) ? t0 : 700;

    /// <summary>Core log approximation.</summary>
    private static (BigInt, long, long) LogStep(Real x, int wp)
    {
        Debug.Assert(x._kind == Kind.Normal && !x.IsNeg);
        if (wp >= LogAgmThreshold)
        {
            // Near 1 the AGM formula cancels; the series path below handles that cheaply.
            long ex = x.BinaryExponent;
            long cancel = 0;
            if (ex is 0 or 1)
            {
                Real d = Subtract(x, One, x.Precision + 2);
                cancel = d._kind == Kind.Zero ? long.MaxValue : -d.BinaryExponent;
            }
            if (cancel < wp / 8)
            {
                return LogAgmStep(x, wp, (int)Math.Max(0, cancel));
            }
        }
        return LogSeriesStep(x, wp);
    }

    /// <summary>
    /// ln x = π / (2 AGM(1, 4/s)) - m ln 2 with s = x 2^m &gt;= 2^(wpi/2) (Brent-Salamin,
    /// as in mpfr_log): O(log p) multiplications and square roots.
    /// </summary>
    private static (BigInt, long, long) LogAgmStep(Real x, int wp, int cancel)
    {
        int lg = 64 - System.Numerics.BitOperations.LeadingZeroCount((ulong)wp);
        int p = wp + cancel + 2 * lg + 32;
        long m = (p + 1) / 2 + 4 - x.BinaryExponent;
        // Floating-point AGM (b0 = 4/s is tiny, so fixed point would lose its relative precision).
        Real a = FromBigInt(BigInt.One, p);
        Real b = Divide(FromBigInt(new BigInt(4)), x.ShiftLeft(m), p);
        int iters = 0;
        while (true)
        {
            Real diff = Subtract(a, b, p);
            // Converged when (a-b)^2/a is below one ulp: then AGM = (a+b)/2 to full precision.
            if (diff._kind == Kind.Zero || diff.BinaryExponent < a.BinaryExponent - p / 2 + 1)
            {
                break;
            }
            Real a1 = Add(a, b, p).ShiftRight(1);
            b = Sqrt(Multiply(a, b, p), p);
            a = a1;
            iters++;
        }
        Real agm2 = Add(a, b, p); // 2 * AGM
        Real pi = FromFixed(PiFixed(p + 8), p + 8, p + 8, RealRounding.NearestEven);
        BigInt lns = Divide(pi, agm2, p).ToFixed(p);
        BigInt v = lns - new BigInt(m) * Ln2Fixed(p);
        // Relative error ~ (iters+4) ulps on ln s (~ p/2 * 0.7 in magnitude), plus m ulps from m ln2.
        long err = 4L * ((iters + 6L) * (p + 64) + Math.Abs(m) + 16);
        return (v, -p, err);
    }

    /// <summary>Log via square roots toward 1 and the atanh series (best at low precision or near 1).</summary>
    private static (BigInt, long, long) LogSeriesStep(Real x, int wp)
    {
        // Square roots cost ~3x a multiply, so use fewer of them than the classic sqrt(wp).
        int t = (int)(Math.Sqrt(wp) * 0.55) + 2;
        int wpi = wp + 2 * t + 48;

        long k = x.BinaryExponent - 1; // x = m * 2^k with m in [1, 2)
        Real m = x.ShiftLeft(-k);
        BigInt mf = m.ToFixed(wpi);
        BigInt one = BigInt.One << wpi;

        if (mf.Equals(one))
        {
            // x is exactly 2^k: ln x = k ln2
            if (k == 0)
            {
                return (BigInt.Zero, 0, -1); // exact zero
            }
            BigInt v0 = new BigInt(k) * Ln2Fixed(wpi);
            return (v0.ShiftRightTruncate(wpi - wp), -wp, 64);
        }

        // t square roots pull m toward 1.
        for (int i = 0; i < t; i++)
        {
            mf = FixSqrt(mf, wpi);
        }
        // atanh series: ln u = 2 atanh((u-1)/(u+1))
        BigInt z = FixDiv(mf - one, mf + one, wpi);
        BigInt z2 = FixSqr(z, wpi);
        BigInt sum = z;
        BigInt term = z;
        for (ulong i = 1; ; i++)
        {
            term = FixMul(term, z2, wpi);
            if (term.IsZero)
            {
                break;
            }
            sum += DivRoundZ(term, 2 * i + 1);
        }
        BigInt lnm = sum.ShiftLeft(t + 1); // * 2^t * 2
        BigInt v = lnm;
        if (k != 0)
        {
            v += new BigInt(k) * Ln2Fixed(wpi);
        }
        long err = 1L << Math.Min(t + 10, 56);
        return (v.ShiftRightTruncate(2 * t + 48), -wp, err);
    }

    /// <summary>Base-2 logarithm (mpfr_log2).</summary>
    public static Real Log2(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
        => LogScaled(x, 2, precision, rounding);

    /// <summary>Base-10 logarithm (mpfr_log10).</summary>
    public static Real Log10(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
        => LogScaled(x, 10, precision, rounding);

    private static Real LogScaled(Real x, int b, int precision, RealRounding rounding)
    {
        int prec = precision != 0 ? precision : x.Precision;
        Real? special = LogSpecial(x, prec);
        if (special.HasValue)
        {
            return special.Value;
        }
        if (b == 2 && IsPowerOfTwoValue(x, out long k2))
        {
            // Exact integer result.
            return FromBigInt(new BigInt(k2), prec, rounding);
        }
        if (b == 10 && x.IsInteger)
        {
            // log10(10^k) = k exactly.
            BigInt v = x.ToBigInt();
            BigInt stripped = BigInt.RemoveFactor(v, BigInt.Ten, out ulong k10);
            if (stripped.IsOne && k10 >= 1)
            {
                return FromBigInt(new BigInt(k10), prec, rounding);
            }
        }
        if (IsPowerOfTwoValue(x, out long k0) && k0 == 0)
        {
            return MakeZero(negative: false, prec);
        }
        return Ziv(prec, rounding, wp =>
        {
            int wpi = wp + 64;
            (BigInt lm, long le, long lerr) = LogStep(x, wpi);
            BigInt den = b == 2 ? Ln2Fixed(wpi) : LogOfSmallInt(b, wpi);
            // (lm * 2^le) / (den * 2^-wpi) — align to fixed point wpi.
            long sh = le + wpi;
            BigInt num = sh >= 0 ? lm.ShiftLeft(sh) : lm.ShiftRightTruncate(-sh);
            BigInt q = FixDiv(num, den, wpi);
            return (q.ShiftRightTruncate(64), -wp, lerr + 16);
        });
    }

    /// <summary>ln(1+x) accurate near zero (mpfr_log1p): 2*atanh(x/(2+x)).</summary>
    public static Real LogP1(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = precision != 0 ? precision : x.Precision;
        switch (x._kind)
        {
            case Kind.NaN:
                return MakeNaN(prec);
            case Kind.Zero:
                return MakeZero(x.IsNeg, prec);
            case Kind.Inf:
                return x.IsNeg ? MakeNaN(prec) : MakeInf(negative: false, prec);
        }
        int cmpm1 = CompareValues(x, MinusOne);
        if (cmpm1 < 0)
        {
            return MakeNaN(prec);
        }
        if (cmpm1 == 0)
        {
            return MakeInf(negative: true, prec);
        }
        return Ziv(prec, rounding, wp =>
        {
            int wpi = wp + 64;
            if (x.BinaryExponent >= 0)
            {
                // |x| >= 1/2: ln(1 + x) directly (the atanh series below would need z = x/(2+x)
                // near ±1). Rounding 1 + x to wpi+160 bits perturbs the log by at most
                // 2^-(wpi+160) absolute, below one ulp of LogStep's internal precision; the
                // result is at least ln(3/2) in magnitude here.
                Real y = Add(One, x, wpi + 160, RealRounding.NearestEven);
                (BigInt lm, long le, long lerr) = LogStep(y, wpi);
                return (lm, le, lerr + 2);
            }
            BigInt xf = x.ToFixed(wpi);
            BigInt two = BigInt.One << (wpi + 1);
            BigInt z = FixDiv(xf, two + xf, wpi);
            if (z.IsZero)
            {
                // |x| far below 2^-wpi: ln(1+x) ≈ x
                return (xf, -wpi, 4);
            }
            BigInt z2 = FixSqr(z, wpi);
            BigInt sum = z;
            BigInt term = z;
            for (ulong i = 1; ; i++)
            {
                term = FixMul(term, z2, wpi);
                if (term.IsZero)
                {
                    break;
                }
                sum += DivRoundZ(term, 2 * i + 1);
            }
            return (sum << 1, -wpi, 64);
        });
    }
}
