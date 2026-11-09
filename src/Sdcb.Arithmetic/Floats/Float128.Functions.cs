namespace Sdcb.Arithmetic;

// Elementary functions through Real. Approach (no double rounding, subnormals included):
// the argument converts exactly to a 113-bit Real; the function is evaluated by Real, which is
// correctly rounded, at WorkPrecision = 128 bits toward zero. Rounding that truncated value to
// binary128 (at the normal or subnormal ulp) is correct for every mode once we know whether the
// truncation was inexact; that only matters when every bit below the target guard bit is zero, in
// which case the same function is evaluated away from zero and the two results are compared
// (equal means exact). Arguments whose results overflow or underflow far beyond the binary128
// range are answered directly so Real never builds astronomically large exponents, and so are
// arguments below 2^-200 for functions whose result there is x(1 ± tiny) or 1 ± tiny.
public readonly partial struct Float128
{
    private const int WorkPrecision = 128;

    private static Float128 Correct<TState>(TState state, Func<TState, int, RealRounding, Real> f, RealRounding rounding)
    {
        if ((uint)rounding > (uint)RealRounding.AwayFromZero)
        {
            throw new ArgumentOutOfRangeException(nameof(rounding));
        }
        Real z = f(state, WorkPrecision, RealRounding.TowardZero);
        if (!z.IsNumber)
        {
            return new(z.ToBinary128Bits(rounding)); // NaN, ±∞ and zeros are exact
        }
        BigInt mant = z.Mantissa;
        long e = z.Exponent;
        int p = z.Precision;
        long lsb = Math.Max(e + p - 1 - FractionBits, -16494);
        long drop = lsb - e; // significand bits below the binary128 ulp (>= p - 113)
        bool inexact = true;
        if (drop - 1 < p && mant.TrailingZeroCount() >= drop - 1)
        {
            Real za = f(state, WorkPrecision, RealRounding.AwayFromZero);
            inexact = za.CompareTo(z) != 0;
        }
        return RoundBig(z.SignBit, mant, e, inexact, rounding);
    }

    private static Float128 Unary(Float128 x, Func<Real, int, RealRounding, Real> f, RealRounding rounding) =>
        Correct(x.ToReal(), f, rounding);

    private static Float128 Binary(Float128 x, Float128 y, Func<(Real X, Real Y), int, RealRounding, Real> f, RealRounding rounding) =>
        Correct((x.ToReal(), y.ToReal()), f, rounding);

    /// <summary>A nonzero value far below half the smallest subnormal, rounded per mode (0 or ±Epsilon).</summary>
    private static Float128 Tiny(bool negative, RealRounding rounding) =>
        RoundWide(negative, UInt128.One, -20000, sticky: true, rounding);

    /// <summary>±(1 - δ) for a tiny positive δ (well below 2^-114), rounded per mode.</summary>
    private static Float128 OneMinusTiny(bool negative, RealRounding rounding) =>
        RoundWide(negative, UInt128.MaxValue, -128, sticky: true, rounding);

    /// <summary>1 ± δ (δ far below 2^-127), rounded per mode.</summary>
    private static Float128 NearOne(bool above, RealRounding rounding) =>
        RoundWide(false, above ? UInt128.One << 127 : (UInt128.One << 127) - 1, -127, sticky: true, rounding);

    /// <summary>x·(1 ± δ) for a tiny relative δ (far below 2^-123), rounded per mode; up = larger magnitude.</summary>
    private static Float128 NearX(Float128 x, bool up, RealRounding rounding)
    {
        int exp = x.BiasedExponent;
        UInt128 sig = exp == 0 ? NormSubnormal(x.Fraction, out exp) : x.Fraction | Hidden;
        UInt128 m = sig << 10;
        if (!up)
        {
            m--;
        }
        return RoundWide(x.SignBit, m, exp - 16495L - 10, sticky: true, rounding);
    }

    /// <summary>
    /// |x| &lt; 2^-200 (nonzero, finite): f(x) = x(1 + O(x²)) or 1 + O(x), so the result is decided by the
    /// sign of the first correction term without evaluating f (Real's Ziv loop is slow for such inputs).
    /// </summary>
    private static bool IsTinyArg(Float128 x) => IsFinite(x) && !IsZero(x) && x.BiasedExponent < ExponentBias - 200;

    // ---------------------------------------------------------------
    // Exponentials
    // ---------------------------------------------------------------

    /// <summary>e^x, correctly rounded (nearest-even).</summary>
    public static Float128 Exp(Float128 x) => Exp(x, RealRounding.NearestEven);

    /// <summary>e^x, correctly rounded.</summary>
    public static Float128 Exp(Float128 x, RealRounding rounding)
    {
        if (IsTinyArg(x))
        {
            return NearOne(!x.SignBit, rounding);
        }
        if (IsFinite(x))
        {
            if (x > (Float128)11357)
            {
                return Huge(false, rounding);
            }
            if (x < (Float128)(-11440))
            {
                return Tiny(false, rounding);
            }
        }
        return Unary(x, static (r, p, m) => Real.Exp(r, p, m), rounding);
    }

    /// <summary>2^x, correctly rounded (nearest-even).</summary>
    public static Float128 Exp2(Float128 x) => Exp2(x, RealRounding.NearestEven);

    /// <summary>2^x, correctly rounded.</summary>
    public static Float128 Exp2(Float128 x, RealRounding rounding)
    {
        if (IsTinyArg(x))
        {
            return NearOne(!x.SignBit, rounding);
        }
        if (IsFinite(x))
        {
            if (x > (Float128)16384)
            {
                return Huge(false, rounding);
            }
            if (x < (Float128)(-16600))
            {
                return Tiny(false, rounding);
            }
        }
        return Unary(x, static (r, p, m) => Real.Exp2(r, p, m), rounding);
    }

    /// <summary>10^x, correctly rounded (nearest-even).</summary>
    public static Float128 Exp10(Float128 x) => Exp10(x, RealRounding.NearestEven);

    /// <summary>10^x, correctly rounded.</summary>
    public static Float128 Exp10(Float128 x, RealRounding rounding)
    {
        if (IsTinyArg(x))
        {
            return NearOne(!x.SignBit, rounding);
        }
        if (IsFinite(x))
        {
            if (x > (Float128)4933)
            {
                return Huge(false, rounding);
            }
            if (x < (Float128)(-4970))
            {
                return Tiny(false, rounding);
            }
        }
        return Unary(x, static (r, p, m) => Real.Exp10(r, p, m), rounding);
    }

    /// <summary>e^x - 1, correctly rounded (nearest-even).</summary>
    public static Float128 ExpM1(Float128 x) => ExpM1(x, RealRounding.NearestEven);

    /// <summary>e^x - 1, correctly rounded.</summary>
    public static Float128 ExpM1(Float128 x, RealRounding rounding)
    {
        if (IsTinyArg(x))
        {
            return NearX(x, !x.SignBit, rounding);
        }
        if (IsFinite(x))
        {
            if (x > (Float128)11357)
            {
                return Huge(false, rounding);
            }
            if (x < (Float128)(-100))
            {
                return OneMinusTiny(true, rounding);
            }
        }
        return Unary(x, static (r, p, m) => Real.ExpM1(r, p, m), rounding);
    }

    /// <summary>2^x - 1, correctly rounded (nearest-even).</summary>
    public static Float128 Exp2M1(Float128 x) => Exp2M1(x, RealRounding.NearestEven);

    /// <summary>2^x - 1, correctly rounded.</summary>
    public static Float128 Exp2M1(Float128 x, RealRounding rounding)
    {
        if (IsFinite(x))
        {
            if (x > (Float128)16384)
            {
                return Huge(false, rounding);
            }
            if (x < (Float128)(-150))
            {
                return OneMinusTiny(true, rounding);
            }
        }
        return Unary(x, static (r, p, m) => Real.Exp2M1(r, p, m), rounding);
    }

    /// <summary>10^x - 1, correctly rounded (nearest-even).</summary>
    public static Float128 Exp10M1(Float128 x) => Exp10M1(x, RealRounding.NearestEven);

    /// <summary>10^x - 1, correctly rounded.</summary>
    public static Float128 Exp10M1(Float128 x, RealRounding rounding)
    {
        if (IsFinite(x))
        {
            if (x > (Float128)4933)
            {
                return Huge(false, rounding);
            }
            if (x < (Float128)(-50))
            {
                return OneMinusTiny(true, rounding);
            }
        }
        return Unary(x, static (r, p, m) => Real.Exp10M1(r, p, m), rounding);
    }

    // ---------------------------------------------------------------
    // Logarithms
    // ---------------------------------------------------------------

    /// <summary>ln x, correctly rounded (nearest-even).</summary>
    public static Float128 Log(Float128 x) => Log(x, RealRounding.NearestEven);

    /// <summary>ln x, correctly rounded.</summary>
    public static Float128 Log(Float128 x, RealRounding rounding) => Unary(x, static (r, p, m) => Real.Log(r, p, m), rounding);

    /// <summary>Logarithm of x in base newBase, correctly rounded (nearest-even).</summary>
    public static Float128 Log(Float128 x, Float128 newBase) => Log(x, newBase, RealRounding.NearestEven);

    /// <summary>Logarithm of x in base newBase, correctly rounded.</summary>
    public static Float128 Log(Float128 x, Float128 newBase, RealRounding rounding) =>
        Binary(x, newBase, static (s, p, m) => Real.Log(s.X, s.Y, p, m), rounding);

    /// <summary>log2 x, correctly rounded (nearest-even).</summary>
    public static Float128 Log2(Float128 x) => Log2(x, RealRounding.NearestEven);

    /// <summary>log2 x, correctly rounded.</summary>
    public static Float128 Log2(Float128 x, RealRounding rounding) => Unary(x, static (r, p, m) => Real.Log2(r, p, m), rounding);

    /// <summary>log10 x, correctly rounded (nearest-even).</summary>
    public static Float128 Log10(Float128 x) => Log10(x, RealRounding.NearestEven);

    /// <summary>log10 x, correctly rounded.</summary>
    public static Float128 Log10(Float128 x, RealRounding rounding) => Unary(x, static (r, p, m) => Real.Log10(r, p, m), rounding);

    /// <summary>ln(1 + x), correctly rounded (nearest-even).</summary>
    public static Float128 LogP1(Float128 x) => LogP1(x, RealRounding.NearestEven);

    /// <summary>ln(1 + x), correctly rounded.</summary>
    public static Float128 LogP1(Float128 x, RealRounding rounding) =>
        IsTinyArg(x) ? NearX(x, x.SignBit, rounding) : Unary(x, static (r, p, m) => Real.LogP1(r, p, m), rounding);

    /// <summary>log2(1 + x), correctly rounded (nearest-even).</summary>
    public static Float128 Log2P1(Float128 x) => Log2P1(x, RealRounding.NearestEven);

    /// <summary>log2(1 + x), correctly rounded.</summary>
    public static Float128 Log2P1(Float128 x, RealRounding rounding) => Unary(x, static (r, p, m) => Real.Log2P1(r, p, m), rounding);

    /// <summary>log10(1 + x), correctly rounded (nearest-even).</summary>
    public static Float128 Log10P1(Float128 x) => Log10P1(x, RealRounding.NearestEven);

    /// <summary>log10(1 + x), correctly rounded.</summary>
    public static Float128 Log10P1(Float128 x, RealRounding rounding) => Unary(x, static (r, p, m) => Real.Log10P1(r, p, m), rounding);

    // ---------------------------------------------------------------
    // Powers and roots
    // ---------------------------------------------------------------

    /// <summary>x^y with IEEE 754 pow special cases, correctly rounded (nearest-even).</summary>
    public static Float128 Pow(Float128 x, Float128 y) => Pow(x, y, RealRounding.NearestEven);

    /// <summary>x^y with IEEE 754 pow special cases, correctly rounded.</summary>
    public static Float128 Pow(Float128 x, Float128 y, RealRounding rounding)
    {
        if (IsFinite(x) && IsFinite(y) && !IsZero(x) && !IsZero(y) && Abs(x) != One && (!x.SignBit || IsInteger(y)))
        {
            // Results far outside the binary128 range are decided from log2|x^y| ≈ y log2|x|.
            double est = y.ToDouble() * Log2Abs(x);
            if (est > 16500 || est < -16700)
            {
                bool neg = x.SignBit && IsOddInteger(y);
                return est > 0 ? Huge(neg, rounding) : Tiny(neg, rounding);
            }
        }
        return Binary(x, y, static (s, p, m) => Real.Pow(s.X, s.Y, p, m), rounding);
    }

    /// <summary>A double estimate of log2|x| (finite nonzero x), accurate also for x near ±1.</summary>
    private static double Log2Abs(Float128 x)
    {
        Float128 a = Abs(x);
        int e = ILogB(a);
        if (e is 0 or -1)
        {
            double d = (a - One).ToDouble(); // exact subtraction (Sterbenz)
            double l = Math.Abs(d) < 1e-8 ? d - d * d / 2 : Math.Log(1 + d);
            return l / Math.Log(2);
        }
        return e + Math.Log2(ScaleB(a, -e).ToDouble());
    }

    /// <summary>The n-th root of x (IEEE 754 rootn), correctly rounded (nearest-even).</summary>
    public static Float128 RootN(Float128 x, int n) => RootN(x, n, RealRounding.NearestEven);

    /// <summary>The n-th root of x (IEEE 754 rootn), correctly rounded.</summary>
    public static Float128 RootN(Float128 x, int n, RealRounding rounding) =>
        Correct((x.ToReal(), (long)n), static (s, p, m) => Real.RootN(s.Item1, s.Item2, p, m), rounding);

    /// <summary>Cube root, correctly rounded (nearest-even).</summary>
    public static Float128 Cbrt(Float128 x) => Cbrt(x, RealRounding.NearestEven);

    /// <summary>Cube root, correctly rounded.</summary>
    public static Float128 Cbrt(Float128 x, RealRounding rounding) => Unary(x, static (r, p, m) => Real.Cbrt(r, p, m), rounding);

    /// <summary>sqrt(x² + y²) without intermediate overflow, correctly rounded (nearest-even).</summary>
    public static Float128 Hypot(Float128 x, Float128 y) => Hypot(x, y, RealRounding.NearestEven);

    /// <summary>sqrt(x² + y²) without intermediate overflow, correctly rounded.</summary>
    public static Float128 Hypot(Float128 x, Float128 y, RealRounding rounding) =>
        Binary(x, y, static (s, p, m) => Real.Hypot(s.X, s.Y, p, m), rounding);

    // ---------------------------------------------------------------
    // Trigonometric functions
    // ---------------------------------------------------------------

    /// <summary>sin x, correctly rounded (nearest-even).</summary>
    public static Float128 Sin(Float128 x) => Sin(x, RealRounding.NearestEven);

    /// <summary>sin x, correctly rounded.</summary>
    public static Float128 Sin(Float128 x, RealRounding rounding) =>
        IsTinyArg(x) ? NearX(x, false, rounding) : Unary(x, static (r, p, m) => Real.Sin(r, p, m), rounding);

    /// <summary>cos x, correctly rounded (nearest-even).</summary>
    public static Float128 Cos(Float128 x) => Cos(x, RealRounding.NearestEven);

    /// <summary>cos x, correctly rounded.</summary>
    public static Float128 Cos(Float128 x, RealRounding rounding) =>
        IsTinyArg(x) ? NearOne(false, rounding) : Unary(x, static (r, p, m) => Real.Cos(r, p, m), rounding);

    /// <summary>tan x, correctly rounded (nearest-even).</summary>
    public static Float128 Tan(Float128 x) => Tan(x, RealRounding.NearestEven);

    /// <summary>tan x, correctly rounded.</summary>
    public static Float128 Tan(Float128 x, RealRounding rounding) =>
        IsTinyArg(x) ? NearX(x, true, rounding) : Unary(x, static (r, p, m) => Real.Tan(r, p, m), rounding);

    /// <summary>(sin x, cos x), each correctly rounded (nearest-even).</summary>
    public static (Float128 Sin, Float128 Cos) SinCos(Float128 x) => (Sin(x), Cos(x));

    /// <summary>asin x, correctly rounded (nearest-even).</summary>
    public static Float128 Asin(Float128 x) => Asin(x, RealRounding.NearestEven);

    /// <summary>asin x, correctly rounded.</summary>
    public static Float128 Asin(Float128 x, RealRounding rounding) =>
        IsTinyArg(x) ? NearX(x, true, rounding) : Unary(x, static (r, p, m) => Real.Asin(r, p, m), rounding);

    /// <summary>acos x, correctly rounded (nearest-even).</summary>
    public static Float128 Acos(Float128 x) => Acos(x, RealRounding.NearestEven);

    /// <summary>acos x, correctly rounded.</summary>
    public static Float128 Acos(Float128 x, RealRounding rounding) => Unary(x, static (r, p, m) => Real.Acos(r, p, m), rounding);

    /// <summary>atan x, correctly rounded (nearest-even).</summary>
    public static Float128 Atan(Float128 x) => Atan(x, RealRounding.NearestEven);

    /// <summary>atan x, correctly rounded.</summary>
    public static Float128 Atan(Float128 x, RealRounding rounding) =>
        IsTinyArg(x) ? NearX(x, false, rounding) : Unary(x, static (r, p, m) => Real.Atan(r, p, m), rounding);

    /// <summary>atan2(y, x), correctly rounded (nearest-even).</summary>
    public static Float128 Atan2(Float128 y, Float128 x) => Atan2(y, x, RealRounding.NearestEven);

    /// <summary>atan2(y, x), correctly rounded.</summary>
    public static Float128 Atan2(Float128 y, Float128 x, RealRounding rounding) =>
        Binary(y, x, static (s, p, m) => Real.Atan2(s.X, s.Y, p, m), rounding);

    /// <summary>sin(πx), correctly rounded (nearest-even).</summary>
    public static Float128 SinPi(Float128 x) => SinPi(x, RealRounding.NearestEven);

    /// <summary>sin(πx), correctly rounded.</summary>
    public static Float128 SinPi(Float128 x, RealRounding rounding) => Unary(x, static (r, p, m) => Real.SinPi(r, p, m), rounding);

    /// <summary>cos(πx), correctly rounded (nearest-even).</summary>
    public static Float128 CosPi(Float128 x) => CosPi(x, RealRounding.NearestEven);

    /// <summary>cos(πx), correctly rounded.</summary>
    public static Float128 CosPi(Float128 x, RealRounding rounding) => Unary(x, static (r, p, m) => Real.CosPi(r, p, m), rounding);

    /// <summary>tan(πx), correctly rounded (nearest-even).</summary>
    public static Float128 TanPi(Float128 x) => TanPi(x, RealRounding.NearestEven);

    /// <summary>tan(πx), correctly rounded.</summary>
    public static Float128 TanPi(Float128 x, RealRounding rounding) => Unary(x, static (r, p, m) => Real.TanPi(r, p, m), rounding);

    /// <summary>(sin πx, cos πx), each correctly rounded (nearest-even).</summary>
    public static (Float128 SinPi, Float128 CosPi) SinCosPi(Float128 x) => (SinPi(x), CosPi(x));

    /// <summary>asin(x)/π, correctly rounded (nearest-even).</summary>
    public static Float128 AsinPi(Float128 x) => AsinPi(x, RealRounding.NearestEven);

    /// <summary>asin(x)/π, correctly rounded.</summary>
    public static Float128 AsinPi(Float128 x, RealRounding rounding) => Unary(x, static (r, p, m) => Real.AsinPi(r, p, m), rounding);

    /// <summary>acos(x)/π, correctly rounded (nearest-even).</summary>
    public static Float128 AcosPi(Float128 x) => AcosPi(x, RealRounding.NearestEven);

    /// <summary>acos(x)/π, correctly rounded.</summary>
    public static Float128 AcosPi(Float128 x, RealRounding rounding) => Unary(x, static (r, p, m) => Real.AcosPi(r, p, m), rounding);

    /// <summary>atan(x)/π, correctly rounded (nearest-even).</summary>
    public static Float128 AtanPi(Float128 x) => AtanPi(x, RealRounding.NearestEven);

    /// <summary>atan(x)/π, correctly rounded.</summary>
    public static Float128 AtanPi(Float128 x, RealRounding rounding) => Unary(x, static (r, p, m) => Real.AtanPi(r, p, m), rounding);

    /// <summary>atan2(y, x)/π, correctly rounded (nearest-even).</summary>
    public static Float128 Atan2Pi(Float128 y, Float128 x) => Atan2Pi(y, x, RealRounding.NearestEven);

    /// <summary>atan2(y, x)/π, correctly rounded.</summary>
    public static Float128 Atan2Pi(Float128 y, Float128 x, RealRounding rounding) =>
        Binary(y, x, static (s, p, m) => Real.Atan2Pi(s.X, s.Y, p, m), rounding);

    // ---------------------------------------------------------------
    // Hyperbolic functions
    // ---------------------------------------------------------------

    /// <summary>sinh x, correctly rounded (nearest-even).</summary>
    public static Float128 Sinh(Float128 x) => Sinh(x, RealRounding.NearestEven);

    /// <summary>sinh x, correctly rounded.</summary>
    public static Float128 Sinh(Float128 x, RealRounding rounding)
    {
        if (IsTinyArg(x))
        {
            return NearX(x, true, rounding);
        }
        if (IsFinite(x) && Abs(x) > (Float128)11358)
        {
            return Huge(x.SignBit, rounding);
        }
        return Unary(x, static (r, p, m) => Real.Sinh(r, p, m), rounding);
    }

    /// <summary>cosh x, correctly rounded (nearest-even).</summary>
    public static Float128 Cosh(Float128 x) => Cosh(x, RealRounding.NearestEven);

    /// <summary>cosh x, correctly rounded.</summary>
    public static Float128 Cosh(Float128 x, RealRounding rounding)
    {
        if (IsTinyArg(x))
        {
            return NearOne(true, rounding);
        }
        if (IsFinite(x) && Abs(x) > (Float128)11358)
        {
            return Huge(false, rounding);
        }
        return Unary(x, static (r, p, m) => Real.Cosh(r, p, m), rounding);
    }

    /// <summary>tanh x, correctly rounded (nearest-even).</summary>
    public static Float128 Tanh(Float128 x) => Tanh(x, RealRounding.NearestEven);

    /// <summary>tanh x, correctly rounded.</summary>
    public static Float128 Tanh(Float128 x, RealRounding rounding)
    {
        if (IsTinyArg(x))
        {
            return NearX(x, false, rounding);
        }
        if (IsFinite(x) && Abs(x) > (Float128)100)
        {
            return OneMinusTiny(x.SignBit, rounding);
        }
        return Unary(x, static (r, p, m) => Real.Tanh(r, p, m), rounding);
    }

    /// <summary>asinh x, correctly rounded (nearest-even).</summary>
    public static Float128 Asinh(Float128 x) => Asinh(x, RealRounding.NearestEven);

    /// <summary>asinh x, correctly rounded.</summary>
    public static Float128 Asinh(Float128 x, RealRounding rounding) =>
        IsTinyArg(x) ? NearX(x, false, rounding) : Unary(x, static (r, p, m) => Real.Asinh(r, p, m), rounding);

    /// <summary>acosh x, correctly rounded (nearest-even).</summary>
    public static Float128 Acosh(Float128 x) => Acosh(x, RealRounding.NearestEven);

    /// <summary>acosh x, correctly rounded.</summary>
    public static Float128 Acosh(Float128 x, RealRounding rounding) => Unary(x, static (r, p, m) => Real.Acosh(r, p, m), rounding);

    /// <summary>atanh x, correctly rounded (nearest-even).</summary>
    public static Float128 Atanh(Float128 x) => Atanh(x, RealRounding.NearestEven);

    /// <summary>atanh x, correctly rounded.</summary>
    public static Float128 Atanh(Float128 x, RealRounding rounding) =>
        IsTinyArg(x) ? NearX(x, true, rounding) : Unary(x, static (r, p, m) => Real.Atanh(r, p, m), rounding);
}
