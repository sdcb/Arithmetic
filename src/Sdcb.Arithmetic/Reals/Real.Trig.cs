using System.Diagnostics;

namespace Sdcb.Arithmetic;

public readonly partial struct Real
{
    // ---------------------------------------------------------------
    // Trigonometric functions (mpfr_sin/cos/tan/asin/acos/atan/atan2)
    // ---------------------------------------------------------------

    /// <summary>Correctly-rounded sine (mpfr_sin).</summary>
    public static Real Sin(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, x.Precision);
        if (x._kind is Kind.NaN or Kind.Inf)
        {
            return MakeNaN(prec);
        }
        if (x._kind == Kind.Zero)
        {
            return MakeZero(x.IsNeg, prec); // sin(±0) = ±0
        }
        return Ziv(prec, rounding, wp =>
        {
            if (x.BinaryExponent < -(wp / 2) - 16)
            {
                return NearX(x, wp, -x.Sign); // sin x = x - x³/6…
            }
            (BigInt s, BigInt c, int wpi, long err) = SinCosCore(x, wp);
            _ = c;
            return (s, (long)-wpi, err);
        });
    }

    /// <summary>Correctly-rounded cosine (mpfr_cos).</summary>
    public static Real Cos(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, x.Precision);
        if (x._kind is Kind.NaN or Kind.Inf)
        {
            return MakeNaN(prec);
        }
        if (x._kind == Kind.Zero)
        {
            return Round(BigInt.One, 0, sticky: false, negative: false, prec, rounding, out _);
        }
        return Ziv(prec, rounding, wp =>
        {
            (BigInt s, BigInt c, int wpi, long err) = SinCosCore(x, wp);
            _ = s;
            return (c, (long)-wpi, err);
        });
    }

    /// <summary>Correctly-rounded tangent (mpfr_tan).</summary>
    public static Real Tan(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, x.Precision);
        if (x._kind is Kind.NaN or Kind.Inf)
        {
            return MakeNaN(prec);
        }
        if (x._kind == Kind.Zero)
        {
            return MakeZero(x.IsNeg, prec);
        }
        return Ziv(prec, rounding, wp =>
        {
            if (x.BinaryExponent < -(wp / 2) - 16)
            {
                return NearX(x, wp, x.Sign); // tan x = x + x³/3…
            }
            (BigInt s, BigInt c, int wpi, long err) = SinCosCore(x, wp + 16);
            if (c.IsZero)
            {
                return (BigInt.Zero, 0, long.MaxValue >> 1); // force retry
            }
            BigInt t = FixDiv(s, c, wpi);
            return (t, (long)-wpi, 4 * err + 8);
        });
    }

    /// <summary>Computes sin and cos with one reduction (mpfr_sin_cos).</summary>
    public static (Real Sin, Real Cos) SinCos(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
        => (Sin(x, precision, rounding), Cos(x, precision, rounding));

    /// <summary>
    /// Core sin/cos evaluation at working precision: reduces mod π/2, evaluates a
    /// short Taylor series after j halvings, then applies double-angle steps.
    /// </summary>
    private static (BigInt s, BigInt c, int wpi, long err) SinCosCore(Real x, int wp)
    {
        int j = (int)(Math.Sqrt(wp) * ExpReductionFactor) + 2;
        int wpi = wp + 3 * j + 48;

        // Reduction: x = k*(π/2) + r, |r| <= π/4.
        long xe = Math.Max(0, x.BinaryExponent);
        int g = (int)Math.Min(xe, 1 << 24) + 32;
        BigInt pi2 = PiFixed(wpi + g).ShiftRightTruncate(1);
        BigInt xf = x.ToFixed(wpi + g);
        BigInt.DivRem(xf, pi2, out BigInt k, out BigInt rr);
        // Pull remainder into [-π/4, π/4].
        if ((rr.Abs() << 1).CompareTo(pi2) > 0)
        {
            if (rr.Sign > 0)
            {
                rr -= pi2;
                k += BigInt.One;
            }
            else
            {
                rr += pi2;
                k -= BigInt.One;
            }
        }
        int quadrant = (int)BigInt.RemainderAbs(k, 4);
        if (k.Sign < 0 && quadrant != 0)
        {
            quadrant = 4 - quadrant;
        }
        BigInt rf = rr.ShiftRightTruncate(g);
        BigInt one = BigInt.One << wpi;
        BigInt s, c;
        long err;
        if (rf.GetBitLength() - wpi < -j)
        {
            // |r| < 2^-j: both series converge fast enough without halving.
            BigInt mr2 = -FixSqr(rf, wpi);
            s = FixMul(rf, SeriesRect(mr2, wpi, static t => (ulong)(2 * t) * (ulong)(2 * t + 1), out long errSin), wpi);
            c = SeriesRect(mr2, wpi, static t => (ulong)(2 * t - 1) * (ulong)(2 * t), out long errCos);
            err = Math.Max(errSin, errCos) + 64;
        }
        else
        {
            // Halve j times and evaluate only cos; double back on v = 1 - cos with
            // v(2θ) = 4v - 2v² (one squaring per step, relative error preserved), then
            // sin r = sign(r) sqrt(v (2 - v)) where |r| >= 2^-j bounds the cancellation.
            BigInt h = rf.ShiftRightTruncate(j);
            BigInt v = one - SeriesRect(-FixSqr(h, wpi), wpi, static t => (ulong)(2 * t - 1) * (ulong)(2 * t), out long errCos);
            for (int i = 0; i < j; i++)
            {
                v = (v << 2) - (FixSqr(v, wpi) << 1);
            }
            c = one - v;
            s = BigInt.Sqrt(v * ((one << 1) - v));
            if (rf.Sign < 0)
            {
                s = -s;
            }
            // v's absolute error grows 4x per doubling; the square root divides it by ~sin r.
            long smallR = Math.Max(0, wpi - rf.GetBitLength());
            err = (errCos + 64) << (int)Math.Min(2 * j + smallR + 3, 50);
        }

        // Apply the quadrant.
        (BigInt rs, BigInt rc) = quadrant switch
        {
            0 => (s, c),
            1 => (c, -s),
            2 => (-s, -c),
            _ => (-c, s),
        };
        return (rs, rc, wpi, err);
    }

    // ---------------------------------------------------------------
    // Inverse trigonometric functions
    // ---------------------------------------------------------------

    /// <summary>Correctly-rounded arctangent (mpfr_atan).</summary>
    public static Real Atan(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, x.Precision);
        switch (x._kind)
        {
            case Kind.NaN:
                return MakeNaN(prec);
            case Kind.Zero:
                return MakeZero(x.IsNeg, prec);
            case Kind.Inf:
                // ±π/2
                return Ziv(prec, rounding, wp =>
                {
                    BigInt v = PiFixed(wp).ShiftRightTruncate(1);
                    return (x.IsNeg ? -v : v, (long)-wp, 64);
                });
        }
        return Ziv(prec, rounding, wp =>
        {
            if (x.BinaryExponent < -(wp / 2) - 16)
            {
                return NearX(x, wp, -x.Sign); // atan x = x - x³/3…
            }
            int wpi = wp + 64;
            BigInt v = AtanFixed(x.Abs().ToFixed(wpi), wpi);
            return (x.IsNeg ? -v : v, (long)-wpi, 1L << 40);
        });
    }

    internal static double AtanReductionFactor = double.TryParse(Environment.GetEnvironmentVariable("SDCB_ARITHMETIC_ATANRED"), out double f) ? f : 0.5;

    /// <summary>atan of a non-negative fixed-point value (any magnitude).</summary>
    private static BigInt AtanFixed(BigInt xf, int wpi)
    {
        Debug.Assert(xf.Sign >= 0);
        BigInt one = BigInt.One << wpi;
        bool invert = xf.CompareTo(one) > 0;
        if (invert)
        {
            xf = FixDiv(one, xf, wpi); // atan x = π/2 - atan(1/x)
        }

        // Halve the argument t times: z <- z / (1 + sqrt(1+z²)). Each halving costs a sqrt
        // and a division (~5 multiplications) while the rectangular-split series costs only
        // ~2√N, so a few halvings suffice.
        int t = (int)(Math.Cbrt(wpi) * AtanReductionFactor) + 1;
        for (int i = 0; i < t; i++)
        {
            BigInt sq = FixSqrt(one + FixSqr(xf, wpi), wpi);
            xf = FixDiv(xf, one + sq, wpi);
        }

        // atan z = z Σ (-z²)^i / (2i+1), term ratio (2t-1)/(2t+1).
        BigInt series = SeriesRect(-FixSqr(xf, wpi), wpi, static k => (ulong)(2 * k - 1), static k => (ulong)(2 * k + 1), out _);
        BigInt sum = FixMul(xf, series, wpi);
        BigInt result = sum.ShiftLeft(t);
        if (invert)
        {
            result = PiFixed(wpi).ShiftRightTruncate(1) - result;
        }
        return result;
    }

    /// <summary>Correctly-rounded arcsine (mpfr_asin); NaN outside [-1, 1].</summary>
    public static Real Asin(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, x.Precision);
        switch (x._kind)
        {
            case Kind.NaN:
            case Kind.Inf:
                return MakeNaN(prec);
            case Kind.Zero:
                return MakeZero(x.IsNeg, prec);
        }
        int cmp = CompareValues(x.Abs(), One);
        if (cmp > 0)
        {
            return MakeNaN(prec);
        }
        if (cmp == 0)
        {
            // ±π/2
            return Ziv(prec, rounding, wp =>
            {
                BigInt v = PiFixed(wp).ShiftRightTruncate(1);
                return (x.IsNeg ? -v : v, (long)-wp, 64);
            });
        }
        return Ziv(prec, rounding, wp =>
        {
            if (x.BinaryExponent < -(wp / 2) - 16)
            {
                return NearX(x, wp, x.Sign); // asin x = x + x³/6…
            }
            // asin x = atan(x / sqrt(1-x²)); extra bits when |x| near 1.
            long extra = Math.Max(0, -(Subtract(One, Square(x.Abs(), 80), 80).BinaryExponent)) / 2 + 8;
            int wpi = wp + (int)Math.Min(extra, 1 << 24) + 64;
            BigInt one = BigInt.One << wpi;
            BigInt xf = x.Abs().ToFixed(wpi);
            BigInt den = FixSqrt(one - FixSqr(xf, wpi), wpi);
            if (den.IsZero)
            {
                return (BigInt.Zero, 0, long.MaxValue >> 1);
            }
            BigInt q = FixDiv(xf, den, wpi);
            BigInt v = AtanFixed(q, wpi);
            return (x.IsNeg ? -v : v, (long)-wpi, 1L << 40);
        });
    }

    /// <summary>Correctly-rounded arccosine (mpfr_acos); NaN outside [-1, 1].</summary>
    public static Real Acos(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, x.Precision);
        switch (x._kind)
        {
            case Kind.NaN:
            case Kind.Inf:
                return MakeNaN(prec);
            case Kind.Zero:
                return Ziv(prec, rounding, wp => (PiFixed(wp).ShiftRightTruncate(1), (long)-wp, 64));
        }
        int cmp = CompareValues(x.Abs(), One);
        if (cmp > 0)
        {
            return MakeNaN(prec);
        }
        if (cmp == 0)
        {
            return x.IsNeg
                ? Ziv(prec, rounding, wp => (PiFixed(wp), (long)-wp, 64))
                : MakeZero(negative: false, prec); // acos(1) = +0
        }
        return Ziv(prec, rounding, wp =>
        {
            // acos x = π/2 - asin-core(x); cancellation near x→1 handled by wp growth.
            long extra = Math.Max(0, -(Subtract(One, x.Abs(), 80).BinaryExponent)) + 16;
            int wpi = wp + (int)Math.Min(extra, 1 << 24) + 64;
            BigInt one = BigInt.One << wpi;
            BigInt xf = x.Abs().ToFixed(wpi);
            BigInt den = FixSqrt(one - FixSqr(xf, wpi), wpi);
            if (den.IsZero)
            {
                return (BigInt.Zero, 0, long.MaxValue >> 1);
            }
            BigInt asinAbs = AtanFixed(FixDiv(xf, den, wpi), wpi);
            BigInt half = PiFixed(wpi).ShiftRightTruncate(1);
            BigInt v = x.IsNeg ? half + asinAbs : half - asinAbs;
            return (v, (long)-wpi, 1L << 40);
        });
    }

    /// <summary>Two-argument arctangent (mpfr_atan2) following the usual quadrant conventions.</summary>
    public static Real Atan2(Real y, Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, Math.Max(y.Precision, x.Precision));
        if (y._kind == Kind.NaN || x._kind == Kind.NaN)
        {
            return MakeNaN(prec);
        }
        // Zero / infinity special cases (mirrors IEEE atan2).
        if (y._kind == Kind.Zero)
        {
            bool xNegSide = x.IsNeg && x._kind != Kind.NaN;
            if (!xNegSide)
            {
                return MakeZero(y.IsNeg, prec);
            }
            return Ziv(prec, rounding, wp => ((y.IsNeg ? -PiFixed(wp) : PiFixed(wp)), (long)-wp, 64));
        }
        if (x._kind == Kind.Zero)
        {
            return Ziv(prec, rounding, wp =>
            {
                BigInt v = PiFixed(wp).ShiftRightTruncate(1);
                return (y.IsNeg ? -v : v, (long)-wp, 64);
            });
        }
        if (y._kind == Kind.Inf || x._kind == Kind.Inf)
        {
            if (x._kind == Kind.Inf && y._kind != Kind.Inf && !x.IsNeg)
            {
                return MakeZero(y.IsNeg, prec); // atan2(finite, +inf) = ±0
            }
            return Ziv(prec, rounding, wp =>
            {
                BigInt pi = PiFixed(wp);
                BigInt v;
                if (y._kind == Kind.Inf && x._kind == Kind.Inf)
                {
                    v = x.IsNeg ? pi - pi.ShiftRightTruncate(2) : pi.ShiftRightTruncate(2); // 3π/4 or π/4
                }
                else if (y._kind == Kind.Inf)
                {
                    v = pi.ShiftRightTruncate(1);
                }
                else
                {
                    v = pi; // x = -inf, finite y: ±π
                }
                return (y.IsNeg ? -v : v, (long)-wp, 64);
            });
        }
        // Both finite nonzero: work on q = |y/x| as a Real to handle extreme ratios.
        return Ziv(prec, rounding, wp =>
        {
            int wpi = wp + 64;
            Real q = Divide(y.Abs(), x.Abs(), wpi + 8, RealRounding.TowardZero);
            (BigInt m, long e, long err) r;
            if (q.BinaryExponent < -(wp / 2) - 16)
            {
                if (!x.IsNeg)
                {
                    r = NearX(q, wp, -1); // atan2 ≈ q - O(q³)
                }
                else
                {
                    r = (PiFixed(wpi) - q.ToFixed(wpi), -wpi, 16); // ≈ π - q
                }
            }
            else if (q.BinaryExponent > wp + 16)
            {
                BigInt half = PiFixed(wpi).ShiftRightTruncate(1);
                Real inv = Divide(One, q, wpi + 8, RealRounding.TowardZero);
                BigInt v = half - (inv._kind == Kind.Normal ? inv.ToFixed(wpi) : BigInt.Zero);
                if (x.IsNeg)
                {
                    v = PiFixed(wpi) - v;
                }
                r = (v, -wpi, 64);
            }
            else
            {
                BigInt a = AtanFixed(q.ToFixed(wpi), wpi);
                if (x.IsNeg)
                {
                    a = PiFixed(wpi) - a;
                }
                r = (a, -wpi, 1L << 41);
            }
            return (y.IsNeg ? -r.m : r.m, r.e, r.err);
        });
    }
}
