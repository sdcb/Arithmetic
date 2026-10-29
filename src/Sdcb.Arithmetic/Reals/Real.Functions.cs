using System.Diagnostics;

namespace Sdcb.Arithmetic;

public readonly partial struct Real
{
    // ---------------------------------------------------------------
    // Helpers for functions built from correctly rounded sub-results
    // ---------------------------------------------------------------

    /// <summary>Turns a Real computed with errUlps ulps of error into a Ziv approximation triple.</summary>
    private static (BigInt, long, long) Approx(Real v, long errUlps)
    {
        if (v._kind != Kind.Normal)
        {
            return (BigInt.Zero, 0, 1); // zero/overflow at this precision: retry higher
        }
        return (v.IsNeg ? -v._man : v._man, v._exp, errUlps);
    }

    private static Real ExactOne(int prec, RealRounding rnd, bool negative = false)
        => Round(BigInt.One, 0, sticky: false, negative, prec, rnd, out _);

    // ---------------------------------------------------------------
    // Reciprocal trigonometric / hyperbolic functions
    // ---------------------------------------------------------------

    /// <summary>Secant 1/cos(x) (mpfr_sec).</summary>
    public static Real Sec(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, x.Precision);
        return x._kind switch
        {
            Kind.NaN or Kind.Inf => MakeNaN(prec),
            Kind.Zero => ExactOne(prec, rounding),
            _ => Ziv(prec, rounding, wp => Approx(Divide(One, Cos(x, wp + 8), wp + 8), 4)),
        };
    }

    /// <summary>Cosecant 1/sin(x) (mpfr_csc).</summary>
    public static Real Csc(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, x.Precision);
        return x._kind switch
        {
            Kind.NaN or Kind.Inf => MakeNaN(prec),
            Kind.Zero => MakeInf(x.IsNeg, prec),
            _ => Ziv(prec, rounding, wp => Approx(Divide(One, Sin(x, wp + 8), wp + 8), 4)),
        };
    }

    /// <summary>Cotangent cos(x)/sin(x) (mpfr_cot).</summary>
    public static Real Cot(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, x.Precision);
        return x._kind switch
        {
            Kind.NaN or Kind.Inf => MakeNaN(prec),
            Kind.Zero => MakeInf(x.IsNeg, prec),
            _ => Ziv(prec, rounding, wp => Approx(Divide(One, Tan(x, wp + 8), wp + 8), 4)),
        };
    }

    /// <summary>Hyperbolic secant 1/cosh(x) (mpfr_sech).</summary>
    public static Real Sech(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, x.Precision);
        return x._kind switch
        {
            Kind.NaN => MakeNaN(prec),
            Kind.Inf => MakeZero(negative: false, prec),
            Kind.Zero => ExactOne(prec, rounding),
            _ => Ziv(prec, rounding, wp => Approx(Divide(One, Cosh(x, wp + 8), wp + 8), 4)),
        };
    }

    /// <summary>Hyperbolic cosecant 1/sinh(x) (mpfr_csch).</summary>
    public static Real Csch(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, x.Precision);
        return x._kind switch
        {
            Kind.NaN => MakeNaN(prec),
            Kind.Inf => MakeZero(x.IsNeg, prec),
            Kind.Zero => MakeInf(x.IsNeg, prec),
            _ => Ziv(prec, rounding, wp => Approx(Divide(One, Sinh(x, wp + 8), wp + 8), 4)),
        };
    }

    /// <summary>Hyperbolic cotangent 1/tanh(x) (mpfr_coth).</summary>
    public static Real Coth(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, x.Precision);
        return x._kind switch
        {
            Kind.NaN => MakeNaN(prec),
            Kind.Inf => ExactOne(prec, rounding, x.IsNeg),
            Kind.Zero => MakeInf(x.IsNeg, prec),
            _ => Ziv(prec, rounding, wp => Approx(Divide(One, Tanh(x, wp + 8), wp + 8), 4)),
        };
    }

    /// <summary>Hyperbolic sine and cosine of x with one call (mpfr_sinh_cosh).</summary>
    public static (Real Sinh, Real Cosh) SinhCosh(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
        => (Sinh(x, precision, rounding), Cosh(x, precision, rounding));

    // ---------------------------------------------------------------
    // sinpi / cospi / tanpi (IEEE 754-2019 sinPi family; mpfr 4.2)
    // ---------------------------------------------------------------

    /// <summary>
    /// Splits a finite x into n + r with n the nearest integer (ties up) and r exact in
    /// [-1/2, 1/2); returns whether n is odd.
    /// </summary>
    private static bool SplitNearestInteger(Real x, out Real r)
    {
        Debug.Assert(x._kind == Kind.Normal);
        if (x.BinaryExponent < 0)
        {
            r = x; // |x| < 1/2
            return false;
        }
        if (x._exp >= 0)
        {
            r = MakeZero(x.IsNeg, x.Precision);
            return x._exp == 0 && x._man.IsOdd;
        }
        long f = -x._exp;
        BigInt m = x.IsNeg ? -x._man : x._man;
        BigInt n = (m + (BigInt.One << (int)Math.Min(f - 1, int.MaxValue))).ShiftRightFloor(f);
        BigInt rm = m - n.ShiftLeft(f);
        r = rm.IsZero ? MakeZero(x.IsNeg, x.Precision) : Round(rm.Abs(), x._exp, sticky: false, rm.Sign < 0, x.Precision + 1, RealRounding.NearestEven, out _);
        return n.IsOdd;
    }

    private static int CompareAbsQuarter(Real r) // |r| vs 1/4
        => r._kind == Kind.Zero ? -1 : CompareValues(r.Abs(), FromBigInt(BigInt.One, MinPrecision).ShiftRight(2));

    /// <summary>sin(πx) (mpfr_sinpi).</summary>
    public static Real SinPi(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
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
        bool odd = SplitNearestInteger(x, out Real r);
        if (r._kind == Kind.Zero)
        {
            return MakeZero(x.IsNeg, prec); // sinPi(n) = ±0 with the sign of x
        }
        bool neg = r.IsNeg ^ odd;
        Real ra = r.Abs();
        if (CompareValues(ra, Half) == 0)
        {
            return ExactOne(prec, rounding, neg);
        }
        bool useCos = CompareAbsQuarter(ra) > 0;
        Real arg = useCos ? Subtract(Half, ra, ra.Precision + 2) : ra; // exact
        return Ziv(prec, rounding, wp =>
        {
            Real t = Multiply(Pi(wp + 16), arg, wp + 16);
            Real v = useCos ? Cos(t, wp + 16) : Sin(t, wp + 16);
            return Approx(neg ? -v : v, 8);
        });
    }

    /// <summary>cos(πx) (mpfr_cospi).</summary>
    public static Real CosPi(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, x.Precision);
        switch (x._kind)
        {
            case Kind.NaN:
            case Kind.Inf:
                return MakeNaN(prec);
            case Kind.Zero:
                return ExactOne(prec, rounding);
        }
        bool odd = SplitNearestInteger(x, out Real r);
        if (r._kind == Kind.Zero)
        {
            return ExactOne(prec, rounding, odd);
        }
        Real ra = r.Abs();
        if (CompareValues(ra, Half) == 0)
        {
            return MakeZero(negative: false, prec); // cosPi(n + 1/2) = +0
        }
        bool useSin = CompareAbsQuarter(ra) > 0;
        Real arg = useSin ? Subtract(Half, ra, ra.Precision + 2) : ra;
        return Ziv(prec, rounding, wp =>
        {
            Real t = Multiply(Pi(wp + 16), arg, wp + 16);
            Real v = useSin ? Sin(t, wp + 16) : Cos(t, wp + 16);
            return Approx(odd ? -v : v, 8);
        });
    }

    /// <summary>tan(πx) (mpfr_tanpi).</summary>
    public static Real TanPi(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
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
        bool odd = SplitNearestInteger(x, out Real r);
        if (r._kind == Kind.Zero)
        {
            // tanPi(n): +0 for positive even / negative odd n, -0 otherwise.
            return MakeZero(x.IsNeg ^ odd, prec);
        }
        Real ra = r.Abs();
        if (CompareValues(ra, Half) == 0)
        {
            // r = -1/2: x = (n-1) + 1/2, and tanPi(k + 1/2) is +inf for even k, -inf for odd k.
            return MakeInf(negative: !odd, prec);
        }
        int q = CompareAbsQuarter(ra);
        if (q == 0)
        {
            return ExactOne(prec, rounding, r.IsNeg);
        }
        bool useCot = q > 0;
        Real arg = useCot ? Subtract(Half, ra, ra.Precision + 2) : ra;
        return Ziv(prec, rounding, wp =>
        {
            Real t = Multiply(Pi(wp + 16), arg, wp + 16);
            Real v = Tan(t, wp + 16);
            if (useCot)
            {
                v = Divide(One, v, wp + 16);
            }
            return Approx(r.IsNeg ? -v : v, 8);
        });
    }

    // ---------------------------------------------------------------
    // 2^x - 1, 10^x - 1, log2(1+x), log10(1+x) (mpfr 4.2)
    // ---------------------------------------------------------------

    /// <summary>2^x - 1 accurate near 0 (mpfr_exp2m1).</summary>
    public static Real Exp2M1(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
        => ExpBaseM1(x, 2, precision, rounding);

    /// <summary>10^x - 1 accurate near 0 (mpfr_exp10m1).</summary>
    public static Real Exp10M1(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
        => ExpBaseM1(x, 10, precision, rounding);

    private static Real ExpBaseM1(Real x, int b, int precision, RealRounding rounding)
    {
        int prec = ResolvePrecision(precision, x.Precision);
        switch (x._kind)
        {
            case Kind.NaN:
                return MakeNaN(prec);
            case Kind.Inf:
                return x.IsNeg ? ExactOne(prec, rounding, negative: true) : MakeInf(negative: false, prec);
            case Kind.Zero:
                return MakeZero(x.IsNeg, prec);
        }
        if (x.IsInteger && x.BinaryExponent <= 24)
        {
            long k = (long)x.ToBigInt();
            if (b == 2)
            {
                // 2^k - 1 is exact (dyadic) for every integer k.
                return Subtract(RoundPow2(k, negative: false, MinPrecision, RealRounding.NearestEven), One, prec, rounding);
            }
            if (k >= 0 && k <= 100000)
            {
                return FromBigInt(BigInt.Pow(10, (ulong)k) - BigInt.One, prec, rounding);
            }
        }
        if (x.BinaryExponent >= 62)
        {
            return x.IsNeg ? Round(BigInt.One, 0, sticky: true, negative: true, prec, rounding, out _) : MakeInf(negative: false, prec);
        }
        return Ziv(prec, rounding, wp =>
        {
            // Relative condition number of expm1 at y = x ln b is about max(1, |y|).
            int wpi = wp + 16 + (int)Math.Max(0, x.BinaryExponent + 4);
            Real lnb = b == 2 ? Ln2(wpi) : Log(FromBigInt(new BigInt(b)), wpi);
            Real y = Multiply(x, lnb, wpi);
            return Approx(ExpM1(y, wpi), 8);
        });
    }

    /// <summary>log2(1 + x) accurate near 0 (mpfr_log2p1).</summary>
    public static Real Log2P1(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
        => LogBaseP1(x, 2, precision, rounding);

    /// <summary>log10(1 + x) accurate near 0 (mpfr_log10p1).</summary>
    public static Real Log10P1(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
        => LogBaseP1(x, 10, precision, rounding);

    private static Real LogBaseP1(Real x, int b, int precision, RealRounding rounding)
    {
        int prec = ResolvePrecision(precision, x.Precision);
        switch (x._kind)
        {
            case Kind.NaN:
                return MakeNaN(prec);
            case Kind.Inf:
                return x.IsNeg ? MakeNaN(prec) : MakeInf(negative: false, prec);
            case Kind.Zero:
                return MakeZero(x.IsNeg, prec);
        }
        int c = CompareValues(x, MinusOne);
        if (c <= 0)
        {
            return c == 0 ? MakeInf(negative: true, prec) : MakeNaN(prec);
        }
        // Exact when 1 + x is a power of the base.
        long bits = Math.Max(x.BinaryExponent, 1) - Math.Min(x._exp, 0) + 2;
        if (bits <= 1 << 20)
        {
            Real y = Add(One, x, (int)bits);
            if (b == 2 && IsPowerOfTwoValue(y, out long k2))
            {
                return FromBigInt(new BigInt(k2), prec, rounding);
            }
            if (b == 10 && y.IsInteger)
            {
                BigInt v = y.ToBigInt();
                BigInt stripped = BigInt.RemoveFactor(v, BigInt.Ten, out ulong k10);
                if (stripped.IsOne && k10 >= 1)
                {
                    return FromBigInt(new BigInt(k10), prec, rounding);
                }
            }
        }
        return Ziv(prec, rounding, wp =>
        {
            int wpi = wp + 16;
            Real l = LogP1(x, wpi);
            Real lnb = b == 2 ? Ln2(wpi) : Log(FromBigInt(new BigInt(b)), wpi);
            return Approx(Divide(l, lnb, wpi), 4);
        });
    }

    /// <summary>Natural logarithm of an unsigned integer (mpfr_log_ui).</summary>
    public static Real Log(ulong n, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, DefaultPrecision);
        if (n == 0)
        {
            return MakeInf(negative: true, prec);
        }
        return Log(FromBigInt(new BigInt(n), 66), prec, rounding);
    }

    // ---------------------------------------------------------------
    // Arithmetic-geometric mean and constants
    // ---------------------------------------------------------------

    /// <summary>Arithmetic-geometric mean of a and b (mpfr_agm); NaN for negative arguments.</summary>
    public static Real ArithmeticGeometricMean(Real a, Real b, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrec(a, b, precision);
        if (a._kind == Kind.NaN || b._kind == Kind.NaN)
        {
            return MakeNaN(prec);
        }
        if ((a.IsNeg && a._kind != Kind.Zero) || (b.IsNeg && b._kind != Kind.Zero))
        {
            return MakeNaN(prec);
        }
        if (a._kind == Kind.Inf || b._kind == Kind.Inf)
        {
            return a._kind == Kind.Zero || b._kind == Kind.Zero ? MakeNaN(prec) : MakeInf(negative: false, prec);
        }
        if (a._kind == Kind.Zero || b._kind == Kind.Zero)
        {
            return MakeZero(negative: false, prec);
        }
        if (CompareValues(a, b) == 0)
        {
            return a.WithPrecision(prec, rounding);
        }
        return Ziv(prec, rounding, wp =>
        {
            int lg = 64 - System.Numerics.BitOperations.LeadingZeroCount((ulong)wp);
            // The iteration count grows with log |log(a/b)|; extra bits cover it.
            long spread = Math.Abs(a.BinaryExponent - b.BinaryExponent);
            int lgs = 64 - System.Numerics.BitOperations.LeadingZeroCount((ulong)spread + 1);
            int p = wp + 2 * (lg + lgs) + 16;
            Real u = a.WithPrecision(p), v = b.WithPrecision(p);
            int iters = 0;
            while (true)
            {
                Real diff = Subtract(u, v, p);
                if (diff._kind == Kind.Zero || diff.BinaryExponent < Math.Max(u.BinaryExponent, v.BinaryExponent) - p / 2 + 1)
                {
                    break;
                }
                Real u1 = Add(u, v, p).ShiftRight(1);
                v = Sqrt(Multiply(u, v, p), p);
                u = u1;
                iters++;
            }
            return Approx(Add(u, v, p).ShiftRight(1), 4L * (iters + 4));
        });
    }

    private static BigInt s_eulerCache;
    private static int s_eulerBits;
    private static BigInt s_catalanCache;
    private static int s_catalanBits;

    /// <summary>EulerGamma-Mascheroni constant γ (mpfr_const_euler).</summary>
    public static Real EulerGamma(int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, DefaultPrecision);
        return Ziv(prec, rounding, wp => (EulerFixed(wp), -wp, 64));
    }

    /// <summary>CatalanConstant's constant G (mpfr_const_catalan).</summary>
    public static Real CatalanConstant(int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, DefaultPrecision);
        return Ziv(prec, rounding, wp => (CatalanFixed(wp), -wp, 64));
    }

    /// <summary>
    /// γ at wp fractional bits by the Brent-McMillan algorithm: with n ~ wp ln2 / 4,
    /// B_k = B_{k-1} n²/k², A_k = (A_{k-1} n²/k + B_k)/k, A_0 = -ln n, B_0 = 1, and
    /// γ = ΣA/ΣB + O(e^{-4n}).
    /// </summary>
    internal static BigInt EulerFixed(int wp)
    {
        lock (ConstLock)
        {
            if (s_eulerBits >= wp)
            {
                return s_eulerCache.ShiftRightTruncate(s_eulerBits - wp);
            }
        }
        int cp = wp + 64;
        ulong n = (ulong)Math.Ceiling(cp * 0.6931471805599453 / 4) + 2;
        // Terms grow to about e^{2n}/(2πn); carry that many integer bits on top of cp.
        int extra = (int)(2.8854 * n) + 16;
        int w = cp + extra;
        Real lnN = Log(FromBigInt(new BigInt(n)), w + 64);
        BigInt a = -lnN.ToFixed(w);
        BigInt bk = BigInt.One << w;
        BigInt sa = a, sb = bk;
        ulong n2 = n * n;
        for (ulong k = 1; ; k++)
        {
            bk = BigInt.DivRem(bk * n2, k * k, out _);
            a = BigInt.DivRem(BigInt.DivRem(a * n2, k, out _) + bk, k, out _);
            if (bk.IsZero && a.IsZero)
            {
                break;
            }
            sa += a;
            sb += bk;
        }
        BigInt g = (sa << cp) / sb;
        lock (ConstLock)
        {
            if (cp > s_eulerBits)
            {
                s_eulerCache = g;
                s_eulerBits = cp;
            }
        }
        return g.ShiftRightTruncate(cp - wp);
    }

    /// <summary>
    /// CatalanConstant's constant at wp bits: G = (π/8) ln(2+√3) + (3/8) Σ_{k&gt;=0} (k!)²/((2k)!(2k+1)²).
    /// </summary>
    internal static BigInt CatalanFixed(int wp)
    {
        lock (ConstLock)
        {
            if (s_catalanBits >= wp)
            {
                return s_catalanCache.ShiftRightTruncate(s_catalanBits - wp);
            }
        }
        int cp = wp + 64;
        // a_k = (k!)²/(2k)!: a_k = a_{k-1} * k / (2(2k-1)).
        BigInt term = BigInt.One << cp;
        BigInt sum = term;
        for (ulong k = 1; ; k++)
        {
            term = BigInt.DivRem(term * k, 2 * (2 * k - 1), out _);
            if (term.IsZero)
            {
                break;
            }
            sum += BigInt.DivRem(term, (2 * k + 1) * (2 * k + 1), out _);
        }
        Real two3 = Add(FromBigInt(BigInt.Two), Sqrt(FromBigInt(new BigInt(3)), cp + 16), cp + 16);
        BigInt l = Log(two3, cp + 16).ToFixed(cp);
        BigInt g = (FixMul(PiFixed(cp), l, cp) >> 3) + ((sum * 3) >> 3);
        lock (ConstLock)
        {
            if (cp > s_catalanBits)
            {
                s_catalanCache = g;
                s_catalanBits = cp;
            }
        }
        return g.ShiftRightTruncate(cp - wp);
    }

    // ---------------------------------------------------------------
    // Factorial and power variants
    // ---------------------------------------------------------------

    /// <summary>n! correctly rounded (mpfr_fac_ui).</summary>
    public static Real Factorial(ulong n, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, DefaultPrecision);
        if (n <= 100000)
        {
            return FromBigInt(BigInt.Factorial(n), prec, rounding);
        }
        // Very large n: exp(lngamma(n+1)).
        return Exp(LogGamma(FromBigInt(new BigInt(n + 1)), prec + 64), prec, rounding);
    }

    /// <summary>x^n for an unsigned integer exponent (mpfr_pow_ui).</summary>
    internal static Real PowUi(Real x, ulong n, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
        => n <= long.MaxValue ? Pow(x, (long)n, precision, rounding) : Pow(x, FromBigInt(new BigInt(n)), precision, rounding);


    /// <summary>x^z for a BigInt exponent (mpfr_pow_z).</summary>
    public static Real Pow(Real x, BigInt z, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        if (z.GetBitLength() <= 62)
        {
            return Pow(x, (long)z, precision, rounding);
        }
        return Pow(x, FromBigInt(z, (int)Math.Max(MinPrecision, z.GetBitLength())), precision, rounding);
    }

    /// <summary>n^m for unsigned integers (mpfr_ui_pow_ui).</summary>
    internal static Real UiPowUi(ulong n, ulong m, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
        => PowUi(FromBigInt(new BigInt(n), 66), m, ResolvePrecision(precision, DefaultPrecision), rounding);

    /// <summary>n^x for an unsigned integer base (mpfr_ui_pow).</summary>
    internal static Real UiPow(ulong n, Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
        => Pow(FromBigInt(new BigInt(n), 66), x, ResolvePrecision(precision, x.Precision), rounding);

    /// <summary>
    /// x^y defined as exp(y log x) (mpfr_powr, IEEE powr): NaN for x &lt; 0, and the
    /// indeterminate forms 0^0, inf^0, 1^inf are NaN.
    /// </summary>
    public static Real Powr(Real x, Real y, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, Math.Max(x.Precision, y.Precision));
        if (x._kind == Kind.NaN || y._kind == Kind.NaN || (x.IsNeg && x._kind != Kind.Zero))
        {
            return MakeNaN(prec);
        }
        bool xZero = x._kind == Kind.Zero, xInf = x._kind == Kind.Inf;
        bool xOne = x._kind == Kind.Normal && CompareValues(x, One) == 0;
        if (y._kind == Kind.Zero)
        {
            return xZero || xInf ? MakeNaN(prec) : ExactOne(prec, rounding);
        }
        if (y._kind == Kind.Inf && xOne)
        {
            return MakeNaN(prec);
        }
        if (xZero)
        {
            // powr(±0, y) = +0 for y > 0, +inf for y < 0.
            return y.IsNeg ? MakeInf(negative: false, prec) : MakeZero(negative: false, prec);
        }
        return Pow(x, y, prec, rounding);
    }

    /// <summary>(1 + x)^n (mpfr_compound_si); NaN for x &lt; -1.</summary>
    public static Real Compound(Real x, long n, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, x.Precision);
        if (n == 0)
        {
            return ExactOne(prec, rounding);
        }
        if (x._kind == Kind.NaN)
        {
            return MakeNaN(prec);
        }
        if (x._kind == Kind.Zero)
        {
            return ExactOne(prec, rounding);
        }
        if (x._kind == Kind.Inf)
        {
            return x.IsNeg ? MakeNaN(prec) : n > 0 ? MakeInf(negative: false, prec) : MakeZero(negative: false, prec);
        }
        int c = CompareValues(x, MinusOne);
        if (c < 0)
        {
            return MakeNaN(prec);
        }
        if (c == 0)
        {
            return n > 0 ? MakeZero(negative: false, prec) : MakeInf(negative: false, prec);
        }
        // Exact base 1 + x when it fits a reasonable precision: then an exact integer power.
        long bits = Math.Max(x.BinaryExponent, 1) - Math.Min(x._exp, 0) + 2;
        if (bits <= 1 << 16)
        {
            return Pow(Add(One, x, (int)bits), n, prec, rounding);
        }
        return Ziv(prec, rounding, wp =>
        {
            int lgn = 64 - System.Numerics.BitOperations.LeadingZeroCount((ulong)Math.Abs(n));
            int wpi = wp + lgn + 16;
            Real l = Multiply(LogP1(x, wpi), FromBigInt(new BigInt(n)), wpi);
            return Approx(Exp(l, wpi), 8);
        });
    }

    /// <summary>k-th root for a signed k (mpfr_rootn_si): x^(1/k).</summary>
    public static Real RootN(Real x, long k, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        if (k > 0 && k <= uint.MaxValue)
        {
            return RootN(x, (uint)k, precision, rounding);
        }
        int prec = ResolvePrecision(precision, x.Precision);
        if (k == 0)
        {
            return MakeNaN(prec);
        }
        bool odd = (k & 1) != 0, recip = k < 0;
        switch (x._kind)
        {
            case Kind.NaN:
                return MakeNaN(prec);
            case Kind.Zero:
                return recip ? MakeInf(x.IsNeg && odd, prec) : MakeZero(x.IsNeg && odd, prec);
            case Kind.Inf:
                if (x.IsNeg && !odd)
                {
                    return MakeNaN(prec);
                }
                return recip ? MakeZero(x.IsNeg && odd, prec) : MakeInf(x.IsNeg && odd, prec);
        }
        if (x.IsNeg && !odd)
        {
            return MakeNaN(prec);
        }
        if (k == -1)
        {
            return Divide(One, x, prec, rounding);
        }
        if (k == -2)
        {
            return ReciprocalSqrt(x, prec, rounding);
        }
        ulong kk = recip ? 0UL - (ulong)k : (ulong)k; // |k|, exact for long.MinValue
        // x^(1/k) is dyadic only when x is a power of two 2^e2 with |k| | e2: return it exactly
        // (Ziv cannot terminate on an exact result).
        if (x._man.IsPowerOfTwo)
        {
            long e2 = x.BinaryExponent - 1;
            ulong ae2 = e2 < 0 ? 0UL - (ulong)e2 : (ulong)e2;
            if (ae2 % kk == 0)
            {
                long q = (long)(ae2 / kk); // kk >= 2, so |q| <= 2^62
                Real p2 = FromBigInt(BigInt.One, prec).ShiftLeft((e2 < 0) != recip ? -q : q);
                return x.IsNeg ? -p2 : p2;
            }
        }
        if (kk <= uint.MaxValue)
        {
            return Ziv(prec, rounding, wp =>
            {
                Real root = RootN(x.Abs(), (uint)kk, wp + 8);
                Real v = Divide(One, root, wp + 8);
                return Approx(x.IsNeg ? -v : v, 4);
            });
        }
        // |k| > 2^32: exp(log|x| / k). |log|x|| < 2^63, so |t| < 2^31 and the absolute error of t
        // (which is the relative error of the result) stays below 2^(31 - wpi).
        return Ziv(prec, rounding, wp =>
        {
            int wpi = wp + 48;
            Real t = Divide(Log(x.Abs(), wpi), FromBigInt(new BigInt(k), 64), wpi);
            Real v = Exp(t, wpi);
            return Approx(x.IsNeg ? -v : v, 8);
        });
    }

    /// <summary>Square root of an unsigned integer (mpfr_sqrt_ui).</summary>
    public static Real Sqrt(ulong n, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
        => Sqrt(FromBigInt(new BigInt(n), 66), ResolvePrecision(precision, DefaultPrecision), rounding);

    // ---------------------------------------------------------------
    // Decomposition and remainders
    // ---------------------------------------------------------------

    /// <summary>Returns the integral part of x and its fractional part in <paramref name="fraction"/>, both with the sign of x (mpfr_modf).</summary>
    public static Real Truncate(Real x, out Real fraction, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, x.Precision);
        switch (x._kind)
        {
            case Kind.NaN:
                fraction = MakeNaN(prec);
                return MakeNaN(prec);
            case Kind.Inf:
                fraction = MakeZero(x.IsNeg, prec);
                return MakeInf(x.IsNeg, prec);
            case Kind.Zero:
                fraction = MakeZero(x.IsNeg, prec);
                return MakeZero(x.IsNeg, prec);
        }
        fraction = x.FractionalPart(prec, rounding);
        return x.Truncate().WithPrecision(prec, rounding);
    }

    /// <summary>
    /// IEEE remainder x - n y (n nearest x/y, ties even) together with the low 63 bits of
    /// n carrying the sign of x/y (mpfr_remquo).
    /// </summary>
    public static Real Ieee754Remainder(Real x, Real y, out long quotient, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        Real r = Ieee754Remainder(x, y, precision, rounding);
        quotient = QuotientLowBits(x, y, nearest: true);
        return r;
    }

    /// <summary>fmod(x, y) together with the low 63 bits of trunc(x/y) carrying its sign (mpfr_fmodquo).</summary>
    public static Real TruncatedRemainder(Real x, Real y, out long quotient, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        Real r = TruncatedRemainder(x, y, precision, rounding);
        quotient = QuotientLowBits(x, y, nearest: false);
        return r;
    }

    /// <summary>Low 63 bits of the integer quotient (nearest-even or truncated) of x/y with its sign.</summary>
    private static long QuotientLowBits(Real x, Real y, bool nearest)
    {
        if (x._kind != Kind.Normal || y._kind != Kind.Normal)
        {
            return 0;
        }
        // At the common scale e: X = |x|/2^e, Y = |y|/2^e; work modulo 2^64 * Y so the
        // quotient's low bits survive even when X is astronomically large.
        long e = Math.Min(x._exp, y._exp);
        BigInt ya = y._man.ShiftLeft(y._exp - e);
        BigInt mod = ya << 64;
        long shift = x._exp - e;
        BigInt xm = shift == 0 ? BigInt.Mod(x._man, mod) : BigInt.Mod(x._man * BigInt.ModPow(BigInt.Two, new BigInt(shift), mod), mod);
        BigInt.DivRem(xm, ya, out BigInt q, out BigInt r);
        if (nearest)
        {
            int c = (r << 1).CompareTo(ya);
            if (c > 0 || (c == 0 && q.IsOdd))
            {
                q += BigInt.One;
            }
        }
        ulong low = (ulong)(q & ((BigInt.One << 63) - 1));
        return x.IsNeg ^ y.IsNeg ? -(long)low : (long)low;
    }

    /// <summary>
    /// Returns m in [1/2, 1) and sets exponent so that x = m 2^exponent (mpfr_frexp);
    /// m is rounded to the requested precision (which may carry it to 1, then renormalized).
    /// </summary>
    public static Real Decompose(Real x, out long exponent, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, x.Precision);
        exponent = 0;
        if (x._kind != Kind.Normal)
        {
            return x.WithPrecisionSpecial(prec);
        }
        Real m = Round(x._man, x._exp, sticky: false, x.IsNeg, prec, rounding, out _);
        exponent = m.BinaryExponent;
        return m.ShiftLeft(-exponent);
    }

    private Real WithPrecisionSpecial(int prec) => _kind switch
    {
        Kind.NaN => MakeNaN(prec),
        Kind.Inf => MakeInf(IsNeg, prec),
        _ => MakeZero(IsNeg, prec),
    };

    // ---------------------------------------------------------------
    // Correctly rounded sums and dot products
    // ---------------------------------------------------------------

    /// <summary>Correctly rounded sum of all values with a single rounding (mpfr_sum).</summary>
    public static Real Sum(ReadOnlySpan<Real> values, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, DefaultPrecision);
        return SumExact(values, prec, rounding);
    }

    /// <summary>Correctly rounded dot product Σ a[i] b[i] with a single rounding (mpfr_dot).</summary>
    public static Real Dot(ReadOnlySpan<Real> a, ReadOnlySpan<Real> b, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        if (a.Length != b.Length)
        {
            throw new ArgumentException("Vectors must have the same length.");
        }
        int prec = ResolvePrecision(precision, DefaultPrecision);
        var products = new Real[a.Length];
        for (int i = 0; i < a.Length; i++)
        {
            Real x = a[i], y = b[i];
            products[i] = x._kind == Kind.Normal && y._kind == Kind.Normal
                ? ExactProduct(x, y)
                : Multiply(x, y, MinPrecision); // special values only
        }
        return SumExact(products, prec, rounding);
    }

    private static Real ExactProduct(Real x, Real y)
    {
        BigInt m = x._man * y._man;
        return new Real(m, x._exp + y._exp, (int)m.GetBitLength(), Kind.Normal, x.IsNeg ^ y.IsNeg);
    }

    private static Real SumExact(ReadOnlySpan<Real> values, int prec, RealRounding rnd)
    {
        int n = values.Length;
        bool posInf = false, negInf = false, nan = false, allZero = true, allNegZero = true, anyNegZero = false;
        long minExp = long.MaxValue;
        for (int i = 0; i < n; i++)
        {
            Real v = values[i];
            switch (v._kind)
            {
                case Kind.NaN:
                    nan = true;
                    break;
                case Kind.Inf:
                    if (v.IsNeg)
                    {
                        negInf = true;
                    }
                    else
                    {
                        posInf = true;
                    }
                    break;
                case Kind.Zero:
                    allNegZero &= v.IsNeg;
                    anyNegZero |= v.IsNeg;
                    break;
                default:
                    allZero = false;
                    allNegZero = false;
                    minExp = Math.Min(minExp, v._exp);
                    break;
            }
        }
        if (nan || (posInf && negInf))
        {
            return MakeNaN(prec);
        }
        if (posInf || negInf)
        {
            return MakeInf(negInf, prec);
        }
        if (allZero)
        {
            // IEEE: -0 when all are -0; mixed signs give +0 (-0 toward -inf); empty sum is +0.
            return MakeZero(n > 0 && (allNegZero || (anyNegZero && rnd == RealRounding.TowardNegative)), prec);
        }
        BigInt acc = BigInt.Zero;
        for (int i = 0; i < n; i++)
        {
            Real v = values[i];
            if (v._kind == Kind.Normal)
            {
                BigInt m = v._man.ShiftLeft(v._exp - minExp);
                acc += v.IsNeg ? -m : m;
            }
        }
        if (acc.IsZero)
        {
            return MakeZero(rnd == RealRounding.TowardNegative, prec);
        }
        return Round(acc.Abs(), minExp, sticky: false, acc.Sign < 0, prec, rnd, out _);
    }

    // ---------------------------------------------------------------
    // Exact construction / deconstruction
    // ---------------------------------------------------------------

    /// <summary>m * 2^exponent rounded to the given precision (mpfr_set_z_2exp).</summary>
    public static Real FromBigInt(BigInt m, long exponent, int precision, RealRounding rounding = RealRounding.NearestEven)
    {
        ValidatePrecision(precision);
        return m.IsZero ? MakeZero(negative: false, precision) : Round(m.Abs(), exponent, sticky: false, m.Sign < 0, precision, rounding, out _);
    }

    /// <summary>A rational number correctly rounded (mpfr_set_q).</summary>
    public static Real FromRational(BigRational q, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, DefaultPrecision);
        BigInt num = q.Numerator, den = q.Denominator;
        if (num.IsZero)
        {
            return MakeZero(negative: false, prec);
        }
        long s = prec + 3 + den.GetBitLength() - num.GetBitLength();
        BigInt.DivRem(s >= 0 ? num.Abs().ShiftLeft(s) : num.Abs(), den, out BigInt qq, out BigInt r);
        if (s < 0)
        {
            s = 0;
        }
        return Round(qq, -s, sticky: !r.IsZero, num.Sign < 0, prec, rounding, out _);
    }

    /// <summary>The exact value as a rational (mpfr_get_q); throws for NaN/Inf.</summary>
    public BigRational ToRational()
    {
        if (_kind is Kind.NaN or Kind.Inf)
        {
            throw new OverflowException("Cannot convert NaN or Infinity to a rational.");
        }
        if (_kind == Kind.Zero)
        {
            return BigRational.Zero;
        }
        BigInt m = IsNeg ? -_man : _man;
        return _exp >= 0 ? new BigRational(m.ShiftLeft(_exp)) : new BigRational(m, BigInt.One.ShiftLeft(-_exp));
    }
}
