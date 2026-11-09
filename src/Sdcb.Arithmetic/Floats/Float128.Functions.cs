namespace Sdcb.Arithmetic;

// Elementary functions through BigFloat. Approach (no double rounding, subnormals included):
// the argument converts exactly to a 113-bit BigFloat; the function is evaluated by BigFloat, which is
// correctly rounded, at WorkPrecision = 128 bits toward zero. Rounding that truncated value to
// binary128 (at the normal or subnormal ulp) is correct for every mode once we know whether the
// truncation was inexact; that only matters when every bit below the target guard bit is zero, in
// which case the same function is evaluated away from zero and the two results are compared
// (equal means exact). Arguments whose results overflow or underflow far beyond the binary128
// range are answered directly so BigFloat never builds astronomically large exponents, and so are
// arguments below 2^-200 for functions whose result there is x(1 ± tiny) or 1 ± tiny (below 2^-116 or
// 2^-60 for the functions whose correction term is that small, see IsTinyArg).
//
// Exp, Exp2, Log, Log2, Log10, Pow, Sin, Cos, Tan, Asin, Acos, Atan, Atan2, Sinh, Cosh and Tanh first try
// a triple-double fast path (Float128.FastPath*.cs) that returns only results it can prove correctly
// rounded; everything else (about 0.03% of typical arguments, exact results in directed modes, special
// values, subnormal or overflowing results, huge trig arguments) continues here.
public readonly partial struct Float128
{
    private const int WorkPrecision = 128;

    private static Float128 Correct<TState>(TState state, Func<TState, int, RoundingMode, BigFloat> f, RoundingMode rounding)
    {
        if ((uint)rounding > (uint)RoundingMode.AwayFromZero)
        {
            throw new ArgumentOutOfRangeException(nameof(rounding));
        }
        BigFloat z = f(state, WorkPrecision, RoundingMode.TowardZero);
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
            BigFloat za = f(state, WorkPrecision, RoundingMode.AwayFromZero);
            inexact = za.CompareTo(z) != 0;
        }
        return RoundBig(z.SignBit, mant, e, inexact, rounding);
    }

    private static Float128 Unary(Float128 x, Func<BigFloat, int, RoundingMode, BigFloat> f, RoundingMode rounding) =>
        Correct(x.ToBigFloat(), f, rounding);

    private static Float128 Binary(Float128 x, Float128 y, Func<(BigFloat X, BigFloat Y), int, RoundingMode, BigFloat> f, RoundingMode rounding) =>
        Correct((x.ToBigFloat(), y.ToBigFloat()), f, rounding);

    /// <summary>A nonzero value far below half the smallest subnormal, rounded per mode (0 or ±Epsilon).</summary>
    private static Float128 Tiny(bool negative, RoundingMode rounding) =>
        RoundWide(negative, UInt128.One, -20000, sticky: true, rounding);

    /// <summary>±(1 - δ) for a tiny positive δ (well below 2^-114), rounded per mode.</summary>
    private static Float128 OneMinusTiny(bool negative, RoundingMode rounding) =>
        RoundWide(negative, UInt128.MaxValue, -128, sticky: true, rounding);

    /// <summary>1 ± δ for 0 &lt; δ &lt; 2^-114 (every such δ rounds alike: the nearest boundaries are 1 + 2^-113 and 1 - 2^-114), rounded per mode.</summary>
    private static Float128 NearOne(bool above, RoundingMode rounding) =>
        RoundWide(false, above ? UInt128.One << 127 : (UInt128.One << 127) - 1, -127, sticky: true, rounding);

    /// <summary>
    /// x·(1 ± δ) for a relative 0 &lt; δ &lt; 2^-116 (no rounding boundary lies that close to a binary128 x, the nearest
    /// being the midpoint 2^-115·|x| below a power of two), rounded per mode; up = larger magnitude.
    /// </summary>
    private static Float128 NearX(Float128 x, bool up, RoundingMode rounding)
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
    /// sign of the first correction term without evaluating f (BigFloat's Ziv loop is slow for such inputs).
    /// </summary>
    private static bool IsTinyArg(Float128 x) => IsTinyArg(x, -200);

    /// <summary>
    /// |x| &lt; 2^maxExp (nonzero, finite). The callers use -116 for f = 1 + c·x (|c| ≤ ln 10, so δ &lt; 2^-114.8)
    /// and x(1 + x/2), and -60 for f = x(1 + c·x²) or 1 + c·x² with |c| ≤ 1/2 (δ &lt; 2^-121), all inside the
    /// ranges where <see cref="NearOne"/> and <see cref="NearX"/> round exactly like f.
    /// </summary>
    private static bool IsTinyArg(Float128 x, int maxExp) => IsFinite(x) && !IsZero(x) && x.BiasedExponent < ExponentBias + maxExp;

    // ---------------------------------------------------------------
    // Exponentials
    // ---------------------------------------------------------------

    /// <summary>e^x, correctly rounded (nearest-even).</summary>
    public static Float128 Exp(Float128 x) => Exp(x, RoundingMode.NearestEven);

    /// <summary>e^x, correctly rounded.</summary>
    public static Float128 Exp(Float128 x, RoundingMode rounding)
    {
        if (IsTinyArg(x, -116))
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
        return TryExpFast(x, rounding, out Float128 fast) ? fast : Unary(x, static (r, p, m) => BigFloat.Exp(r, p, m), rounding);
    }

    /// <summary>2^x, correctly rounded (nearest-even).</summary>
    public static Float128 Exp2(Float128 x) => Exp2(x, RoundingMode.NearestEven);

    /// <summary>2^x, correctly rounded.</summary>
    public static Float128 Exp2(Float128 x, RoundingMode rounding)
    {
        if (IsTinyArg(x, -116))
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
        return TryExp2Fast(x, rounding, out Float128 fast) ? fast : Unary(x, static (r, p, m) => BigFloat.Exp2(r, p, m), rounding);
    }

    /// <summary>10^x, correctly rounded (nearest-even).</summary>
    public static Float128 Exp10(Float128 x) => Exp10(x, RoundingMode.NearestEven);

    /// <summary>10^x, correctly rounded.</summary>
    public static Float128 Exp10(Float128 x, RoundingMode rounding)
    {
        if (IsTinyArg(x, -116))
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
        return Unary(x, static (r, p, m) => BigFloat.Exp10(r, p, m), rounding);
    }

    /// <summary>e^x - 1, correctly rounded (nearest-even).</summary>
    public static Float128 ExpM1(Float128 x) => ExpM1(x, RoundingMode.NearestEven);

    /// <summary>e^x - 1, correctly rounded.</summary>
    public static Float128 ExpM1(Float128 x, RoundingMode rounding)
    {
        if (IsTinyArg(x, -116))
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
        return Unary(x, static (r, p, m) => BigFloat.ExpM1(r, p, m), rounding);
    }

    /// <summary>2^x - 1, correctly rounded (nearest-even).</summary>
    public static Float128 Exp2M1(Float128 x) => Exp2M1(x, RoundingMode.NearestEven);

    /// <summary>2^x - 1, correctly rounded.</summary>
    public static Float128 Exp2M1(Float128 x, RoundingMode rounding)
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
        return Unary(x, static (r, p, m) => BigFloat.Exp2M1(r, p, m), rounding);
    }

    /// <summary>10^x - 1, correctly rounded (nearest-even).</summary>
    public static Float128 Exp10M1(Float128 x) => Exp10M1(x, RoundingMode.NearestEven);

    /// <summary>10^x - 1, correctly rounded.</summary>
    public static Float128 Exp10M1(Float128 x, RoundingMode rounding)
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
        return Unary(x, static (r, p, m) => BigFloat.Exp10M1(r, p, m), rounding);
    }

    // ---------------------------------------------------------------
    // Logarithms
    // ---------------------------------------------------------------

    /// <summary>ln x, correctly rounded (nearest-even).</summary>
    public static Float128 Log(Float128 x) => Log(x, RoundingMode.NearestEven);

    /// <summary>ln x, correctly rounded.</summary>
    public static Float128 Log(Float128 x, RoundingMode rounding) =>
        TryLogFast(x, rounding, out Float128 fast) ? fast : Unary(x, static (r, p, m) => BigFloat.Log(r, p, m), rounding);

    /// <summary>Logarithm of x in base newBase, correctly rounded (nearest-even).</summary>
    public static Float128 Log(Float128 x, Float128 newBase) => Log(x, newBase, RoundingMode.NearestEven);

    /// <summary>Logarithm of x in base newBase, correctly rounded.</summary>
    public static Float128 Log(Float128 x, Float128 newBase, RoundingMode rounding) =>
        Binary(x, newBase, static (s, p, m) => BigFloat.Log(s.X, s.Y, p, m), rounding);

    /// <summary>log2 x, correctly rounded (nearest-even).</summary>
    public static Float128 Log2(Float128 x) => Log2(x, RoundingMode.NearestEven);

    /// <summary>log2 x, correctly rounded.</summary>
    public static Float128 Log2(Float128 x, RoundingMode rounding) =>
        TryLog2Fast(x, rounding, out Float128 fast) ? fast : Unary(x, static (r, p, m) => BigFloat.Log2(r, p, m), rounding);

    /// <summary>log10 x, correctly rounded (nearest-even).</summary>
    public static Float128 Log10(Float128 x) => Log10(x, RoundingMode.NearestEven);

    /// <summary>log10 x, correctly rounded.</summary>
    public static Float128 Log10(Float128 x, RoundingMode rounding) =>
        TryLog10Fast(x, rounding, out Float128 fast) ? fast : Unary(x, static (r, p, m) => BigFloat.Log10(r, p, m), rounding);

    /// <summary>ln(1 + x), correctly rounded (nearest-even).</summary>
    public static Float128 LogP1(Float128 x) => LogP1(x, RoundingMode.NearestEven);

    /// <summary>ln(1 + x), correctly rounded.</summary>
    public static Float128 LogP1(Float128 x, RoundingMode rounding) =>
        IsTinyArg(x, -116) ? NearX(x, x.SignBit, rounding) : Unary(x, static (r, p, m) => BigFloat.LogP1(r, p, m), rounding);

    /// <summary>log2(1 + x), correctly rounded (nearest-even).</summary>
    public static Float128 Log2P1(Float128 x) => Log2P1(x, RoundingMode.NearestEven);

    /// <summary>log2(1 + x), correctly rounded.</summary>
    public static Float128 Log2P1(Float128 x, RoundingMode rounding) => Unary(x, static (r, p, m) => BigFloat.Log2P1(r, p, m), rounding);

    /// <summary>log10(1 + x), correctly rounded (nearest-even).</summary>
    public static Float128 Log10P1(Float128 x) => Log10P1(x, RoundingMode.NearestEven);

    /// <summary>log10(1 + x), correctly rounded.</summary>
    public static Float128 Log10P1(Float128 x, RoundingMode rounding) => Unary(x, static (r, p, m) => BigFloat.Log10P1(r, p, m), rounding);

    // ---------------------------------------------------------------
    // Powers and roots
    // ---------------------------------------------------------------

    /// <summary>x^y with IEEE 754 pow special cases, correctly rounded (nearest-even).</summary>
    public static Float128 Pow(Float128 x, Float128 y) => Pow(x, y, RoundingMode.NearestEven);

    /// <summary>x^y with IEEE 754 pow special cases, correctly rounded.</summary>
    public static Float128 Pow(Float128 x, Float128 y, RoundingMode rounding)
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
        return TryPowFast(x, y, rounding, out Float128 fast) ? fast : Binary(x, y, static (s, p, m) => BigFloat.Pow(s.X, s.Y, p, m), rounding);
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
    public static Float128 RootN(Float128 x, int n) => RootN(x, n, RoundingMode.NearestEven);

    /// <summary>The n-th root of x (IEEE 754 rootn), correctly rounded.</summary>
    public static Float128 RootN(Float128 x, int n, RoundingMode rounding) =>
        Correct((x.ToBigFloat(), (long)n), static (s, p, m) => BigFloat.RootN(s.Item1, s.Item2, p, m), rounding);

    /// <summary>Cube root, correctly rounded (nearest-even).</summary>
    public static Float128 Cbrt(Float128 x) => Cbrt(x, RoundingMode.NearestEven);

    /// <summary>Cube root, correctly rounded.</summary>
    public static Float128 Cbrt(Float128 x, RoundingMode rounding) => Unary(x, static (r, p, m) => BigFloat.Cbrt(r, p, m), rounding);

    /// <summary>sqrt(x² + y²) without intermediate overflow, correctly rounded (nearest-even).</summary>
    public static Float128 Hypot(Float128 x, Float128 y) => Hypot(x, y, RoundingMode.NearestEven);

    /// <summary>sqrt(x² + y²) without intermediate overflow, correctly rounded.</summary>
    public static Float128 Hypot(Float128 x, Float128 y, RoundingMode rounding) =>
        Binary(x, y, static (s, p, m) => BigFloat.Hypot(s.X, s.Y, p, m), rounding);

    // ---------------------------------------------------------------
    // Trigonometric functions
    // ---------------------------------------------------------------

    /// <summary>sin x, correctly rounded (nearest-even).</summary>
    public static Float128 Sin(Float128 x) => Sin(x, RoundingMode.NearestEven);

    /// <summary>sin x, correctly rounded.</summary>
    public static Float128 Sin(Float128 x, RoundingMode rounding) =>
        IsTinyArg(x, -60) ? NearX(x, false, rounding) : TrySinFast(x, rounding, out Float128 fast) ? fast : Unary(x, static (r, p, m) => BigFloat.Sin(r, p, m), rounding);

    /// <summary>cos x, correctly rounded (nearest-even).</summary>
    public static Float128 Cos(Float128 x) => Cos(x, RoundingMode.NearestEven);

    /// <summary>cos x, correctly rounded.</summary>
    public static Float128 Cos(Float128 x, RoundingMode rounding) =>
        IsTinyArg(x, -60) ? NearOne(false, rounding) : TryCosFast(x, rounding, out Float128 fast) ? fast : Unary(x, static (r, p, m) => BigFloat.Cos(r, p, m), rounding);

    /// <summary>tan x, correctly rounded (nearest-even).</summary>
    public static Float128 Tan(Float128 x) => Tan(x, RoundingMode.NearestEven);

    /// <summary>tan x, correctly rounded.</summary>
    public static Float128 Tan(Float128 x, RoundingMode rounding) =>
        IsTinyArg(x, -60) ? NearX(x, true, rounding) : TryTanFast(x, rounding, out Float128 fast) ? fast : Unary(x, static (r, p, m) => BigFloat.Tan(r, p, m), rounding);

    /// <summary>(sin x, cos x), each correctly rounded (nearest-even).</summary>
    public static (Float128 Sin, Float128 Cos) SinCos(Float128 x) => (Sin(x), Cos(x));

    /// <summary>asin x, correctly rounded (nearest-even).</summary>
    public static Float128 Asin(Float128 x) => Asin(x, RoundingMode.NearestEven);

    /// <summary>asin x, correctly rounded.</summary>
    public static Float128 Asin(Float128 x, RoundingMode rounding) =>
        IsTinyArg(x, -60) ? NearX(x, true, rounding) : TryAsinFast(x, rounding, out Float128 fast) ? fast : Unary(x, static (r, p, m) => BigFloat.Asin(r, p, m), rounding);

    /// <summary>acos x, correctly rounded (nearest-even).</summary>
    public static Float128 Acos(Float128 x) => Acos(x, RoundingMode.NearestEven);

    /// <summary>acos x, correctly rounded.</summary>
    public static Float128 Acos(Float128 x, RoundingMode rounding) =>
        TryAcosFast(x, rounding, out Float128 fast) ? fast : Unary(x, static (r, p, m) => BigFloat.Acos(r, p, m), rounding);

    /// <summary>atan x, correctly rounded (nearest-even).</summary>
    public static Float128 Atan(Float128 x) => Atan(x, RoundingMode.NearestEven);

    /// <summary>atan x, correctly rounded.</summary>
    public static Float128 Atan(Float128 x, RoundingMode rounding) =>
        IsTinyArg(x, -60) ? NearX(x, false, rounding) : TryAtanFast(x, rounding, out Float128 fast) ? fast : Unary(x, static (r, p, m) => BigFloat.Atan(r, p, m), rounding);

    /// <summary>atan2(y, x), correctly rounded (nearest-even).</summary>
    public static Float128 Atan2(Float128 y, Float128 x) => Atan2(y, x, RoundingMode.NearestEven);

    /// <summary>atan2(y, x), correctly rounded.</summary>
    public static Float128 Atan2(Float128 y, Float128 x, RoundingMode rounding) =>
        TryAtan2Fast(y, x, rounding, out Float128 fast) ? fast : Binary(y, x, static (s, p, m) => BigFloat.Atan2(s.X, s.Y, p, m), rounding);

    /// <summary>sin(πx), correctly rounded (nearest-even).</summary>
    public static Float128 SinPi(Float128 x) => SinPi(x, RoundingMode.NearestEven);

    /// <summary>sin(πx), correctly rounded.</summary>
    public static Float128 SinPi(Float128 x, RoundingMode rounding) => Unary(x, static (r, p, m) => BigFloat.SinPi(r, p, m), rounding);

    /// <summary>cos(πx), correctly rounded (nearest-even).</summary>
    public static Float128 CosPi(Float128 x) => CosPi(x, RoundingMode.NearestEven);

    /// <summary>cos(πx), correctly rounded.</summary>
    public static Float128 CosPi(Float128 x, RoundingMode rounding) => Unary(x, static (r, p, m) => BigFloat.CosPi(r, p, m), rounding);

    /// <summary>tan(πx), correctly rounded (nearest-even).</summary>
    public static Float128 TanPi(Float128 x) => TanPi(x, RoundingMode.NearestEven);

    /// <summary>tan(πx), correctly rounded.</summary>
    public static Float128 TanPi(Float128 x, RoundingMode rounding) => Unary(x, static (r, p, m) => BigFloat.TanPi(r, p, m), rounding);

    /// <summary>(sin πx, cos πx), each correctly rounded (nearest-even).</summary>
    public static (Float128 SinPi, Float128 CosPi) SinCosPi(Float128 x) => (SinPi(x), CosPi(x));

    /// <summary>asin(x)/π, correctly rounded (nearest-even).</summary>
    public static Float128 AsinPi(Float128 x) => AsinPi(x, RoundingMode.NearestEven);

    /// <summary>asin(x)/π, correctly rounded.</summary>
    public static Float128 AsinPi(Float128 x, RoundingMode rounding) => Unary(x, static (r, p, m) => BigFloat.AsinPi(r, p, m), rounding);

    /// <summary>acos(x)/π, correctly rounded (nearest-even).</summary>
    public static Float128 AcosPi(Float128 x) => AcosPi(x, RoundingMode.NearestEven);

    /// <summary>acos(x)/π, correctly rounded.</summary>
    public static Float128 AcosPi(Float128 x, RoundingMode rounding) => Unary(x, static (r, p, m) => BigFloat.AcosPi(r, p, m), rounding);

    /// <summary>atan(x)/π, correctly rounded (nearest-even).</summary>
    public static Float128 AtanPi(Float128 x) => AtanPi(x, RoundingMode.NearestEven);

    /// <summary>atan(x)/π, correctly rounded.</summary>
    public static Float128 AtanPi(Float128 x, RoundingMode rounding) => Unary(x, static (r, p, m) => BigFloat.AtanPi(r, p, m), rounding);

    /// <summary>atan2(y, x)/π, correctly rounded (nearest-even).</summary>
    public static Float128 Atan2Pi(Float128 y, Float128 x) => Atan2Pi(y, x, RoundingMode.NearestEven);

    /// <summary>atan2(y, x)/π, correctly rounded.</summary>
    public static Float128 Atan2Pi(Float128 y, Float128 x, RoundingMode rounding) =>
        Binary(y, x, static (s, p, m) => BigFloat.Atan2Pi(s.X, s.Y, p, m), rounding);

    // ---------------------------------------------------------------
    // Hyperbolic functions
    // ---------------------------------------------------------------

    /// <summary>sinh x, correctly rounded (nearest-even).</summary>
    public static Float128 Sinh(Float128 x) => Sinh(x, RoundingMode.NearestEven);

    /// <summary>sinh x, correctly rounded.</summary>
    public static Float128 Sinh(Float128 x, RoundingMode rounding)
    {
        if (IsTinyArg(x, -60))
        {
            return NearX(x, true, rounding);
        }
        if (IsFinite(x) && Abs(x) > (Float128)11358)
        {
            return Huge(x.SignBit, rounding);
        }
        return TrySinhFast(x, rounding, out Float128 fast) ? fast : Unary(x, static (r, p, m) => BigFloat.Sinh(r, p, m), rounding);
    }

    /// <summary>cosh x, correctly rounded (nearest-even).</summary>
    public static Float128 Cosh(Float128 x) => Cosh(x, RoundingMode.NearestEven);

    /// <summary>cosh x, correctly rounded.</summary>
    public static Float128 Cosh(Float128 x, RoundingMode rounding)
    {
        if (IsTinyArg(x, -60))
        {
            return NearOne(true, rounding);
        }
        if (IsFinite(x) && Abs(x) > (Float128)11358)
        {
            return Huge(false, rounding);
        }
        return TryCoshFast(x, rounding, out Float128 fast) ? fast : Unary(x, static (r, p, m) => BigFloat.Cosh(r, p, m), rounding);
    }

    /// <summary>tanh x, correctly rounded (nearest-even).</summary>
    public static Float128 Tanh(Float128 x) => Tanh(x, RoundingMode.NearestEven);

    /// <summary>tanh x, correctly rounded.</summary>
    public static Float128 Tanh(Float128 x, RoundingMode rounding)
    {
        if (IsTinyArg(x, -60))
        {
            return NearX(x, false, rounding);
        }
        if (IsFinite(x) && Abs(x) >= (Float128)41)
        {
            // tanh x = 1 - δ with 0 < δ < 2e^-82 < 2^-117: every mode rounds it like 1 - tiny.
            return OneMinusTiny(x.SignBit, rounding);
        }
        return TryTanhFast(x, rounding, out Float128 fast) ? fast : Unary(x, static (r, p, m) => BigFloat.Tanh(r, p, m), rounding);
    }

    /// <summary>asinh x, correctly rounded (nearest-even).</summary>
    public static Float128 Asinh(Float128 x) => Asinh(x, RoundingMode.NearestEven);

    /// <summary>asinh x, correctly rounded.</summary>
    public static Float128 Asinh(Float128 x, RoundingMode rounding) =>
        IsTinyArg(x, -60) ? NearX(x, false, rounding) : Unary(x, static (r, p, m) => BigFloat.Asinh(r, p, m), rounding);

    /// <summary>acosh x, correctly rounded (nearest-even).</summary>
    public static Float128 Acosh(Float128 x) => Acosh(x, RoundingMode.NearestEven);

    /// <summary>acosh x, correctly rounded.</summary>
    public static Float128 Acosh(Float128 x, RoundingMode rounding) => Unary(x, static (r, p, m) => BigFloat.Acosh(r, p, m), rounding);

    /// <summary>atanh x, correctly rounded (nearest-even).</summary>
    public static Float128 Atanh(Float128 x) => Atanh(x, RoundingMode.NearestEven);

    /// <summary>atanh x, correctly rounded.</summary>
    public static Float128 Atanh(Float128 x, RoundingMode rounding) =>
        IsTinyArg(x, -60) ? NearX(x, true, rounding) : Unary(x, static (r, p, m) => BigFloat.Atanh(r, p, m), rounding);
}
