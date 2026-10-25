using System.Diagnostics;

namespace Sdcb.Arithmetic;

public readonly partial struct Real
{
    // ---------------------------------------------------------------
    // Infrastructure for special functions
    //
    // Each function computes at internal precision p = wp + guard with round-to-nearest
    // operations, and reports the largest magnitude (log2) among the terms it summed.
    // When the result is smaller than that by more than the guard allows (cancellation),
    // the evaluation is repeated with more guard bits. Ziv then sees only the trusted top
    // wp+20 bits with a generous error, so the rounding test is robust.
    // ---------------------------------------------------------------

    /// <summary>One evaluation at precision p: the value and log2 of the largest summand magnitude.</summary>
    private delegate (Real Value, long MaxLog2) Evaluator(int p);

    private static (BigInt, long, long) Guarded(int wp, int initialGuard, Evaluator eval)
    {
        int guard = Math.Max(initialGuard, 64);
        while (true)
        {
            int p = (int)Math.Min((long)wp + guard, MaxPrecision);
            (Real v, long maxLog2) = eval(p);
            if (v._kind != Kind.Normal)
            {
                if (guard > wp + (1 << 20))
                {
                    return (BigInt.Zero, 0, -1); // treat as exact zero after exhausting precision
                }
                guard = guard * 2 + 64;
                continue;
            }
            long cancel = Math.Max(0, maxLog2 - v.BinaryExponent);
            if (cancel + 44 <= guard || p == MaxPrecision)
            {
                return Trusted(v, wp + 20);
            }
            guard = (int)Math.Min(cancel + 96, (long)MaxPrecision);
        }
    }

    /// <summary>Keeps the top bits of v as a Ziv approximation with a 64-ulp error.</summary>
    private static (BigInt, long, long) Trusted(Real v, int bits)
    {
        long drop = Math.Max(0, v.Precision - bits);
        BigInt m = drop == 0 ? v._man : v._man.ShiftRightTruncate(drop);
        return (v.IsNeg ? -m : m, v._exp + drop, 64);
    }

    private static long Log2Of(Real v) => v._kind == Kind.Normal ? v.BinaryExponent : long.MinValue / 4;

    /// <summary>a + b computed exactly (precision spans both operands; capped at 2^24 bits).</summary>
    private static Real ExactAdd(Real a, Real b)
    {
        if (a._kind != Kind.Normal || b._kind != Kind.Normal)
        {
            return Add(a, b, Math.Max(a.Precision, b.Precision));
        }
        long hi = Math.Max(a.BinaryExponent, b.BinaryExponent), lo = Math.Min(a._exp, b._exp);
        return Add(a, b, (int)Math.Clamp(hi - lo + 2, MinPrecision, 1 << 24));
    }

    private static Real Int(long n) => FromBigInt(new BigInt(n), 66);

    private static Real Int(BigInt n) => FromBigInt(n, (int)Math.Max(MinPrecision, n.GetBitLength()));

    // ---------------------------------------------------------------
    // Bernoulli numbers via tangent numbers (Brent-Harvey)
    // ---------------------------------------------------------------

    private static BigInt[] s_tangent = [];
    private static readonly object TangentLock = new();

    /// <summary>Tangent numbers T_1..T_n (index k holds T_k): B_{2k} = (-1)^(k-1) 2k T_k / (4^k (4^k - 1)).</summary>
    private static BigInt[] TangentNumbers(int n)
    {
        lock (TangentLock)
        {
            if (s_tangent.Length > n)
            {
                return s_tangent;
            }
            int m = Math.Max(n, 2 * (s_tangent.Length - 1));
            var t = new BigInt[m + 1];
            t[1] = BigInt.One;
            for (int k = 2; k <= m; k++)
            {
                t[k] = t[k - 1] * (k - 1);
            }
            for (int k = 2; k <= m; k++)
            {
                for (int j = k; j <= m; j++)
                {
                    t[j] = t[j - 1] * (j - k) + t[j] * (j - k + 2);
                }
            }
            s_tangent = t;
            return t;
        }
    }

    /// <summary>B_{2k} / (2k) as a Real at precision p.</summary>
    private static Real BernoulliOver2k(int k, int p)
    {
        BigInt t = TangentNumbers(k)[k];
        BigInt d = (BigInt.One << (2 * k)) * ((BigInt.One << (2 * k)) - 1);
        Real v = Divide(Int(t), Int(d), p);
        return (k & 1) == 0 ? -v : v;
    }

    // ---------------------------------------------------------------
    // Gamma family
    // ---------------------------------------------------------------

    /// <summary>
    /// ln Γ(x) for x &gt; 0 at precision p by shifting to X = x + N and Stirling's series:
    /// lnΓ(X) = (X - 1/2) ln X - X + ln(2π)/2 + Σ B_{2k}/(2k(2k-1) X^{2k-1}).
    /// </summary>
    private static (Real Value, long MaxLog2) LnGammaPositive(Real x, int p)
    {
        Debug.Assert(x._kind == Kind.Normal && !x.IsNeg);
        double xd = Math.Min(x.ToDouble(), 1e300);
        double x0 = Math.Max(10, p * 0.3);
        long n = xd >= x0 ? 0 : (long)Math.Ceiling(x0 - xd);
        int pp = p + 2 * BitLen(n + 1) + 8;
        Real bigX = n == 0 ? x : Add(x, Int(n), pp);
        // Product P = x (x+1) ... (x+N-1), taken in log.
        Real lnP = MakeZero(negative: false, pp);
        if (n > 0)
        {
            Real prod = x.WithPrecision(pp);
            for (long i = 1; i < n; i++)
            {
                prod = Multiply(prod, Add(x, Int(i), pp), pp);
            }
            lnP = Log(prod, pp);
        }
        Real lnX = Log(bigX, pp);
        Real main = Subtract(Multiply(Subtract(bigX, Half, pp), lnX, pp), bigX, pp);
        Real half2pi = Log(Multiply(Pi(pp), FromBigInt(BigInt.Two), pp), pp).ShiftRight(1);
        Real sum = Add(main, half2pi, pp);
        long maxLog = Math.Max(Log2Of(main), Log2Of(lnP));
        // Stirling tail: Σ_k B_{2k}/(2k(2k-1)) X^{1-2k}; enveloping, stop below 2^-(pp+8) of sum.
        Real invX = Divide(One, bigX, pp), invX2 = Square(invX, pp);
        Real pw = invX;
        long tol = Math.Max(Log2Of(sum), 0) - pp - 8;
        Real prevTerm = default;
        for (int k = 1; k < 1 << 20; k++)
        {
            Real c = Divide(BernoulliOver2k(k, pp), Int(2 * k - 1), pp);
            Real term = Multiply(c, pw, pp);
            if (term._kind != Kind.Normal || term.BinaryExponent < tol)
            {
                break;
            }
            if (k > 2 && prevTerm._kind == Kind.Normal && term.BinaryExponent > prevTerm.BinaryExponent)
            {
                maxLog = Log2Of(sum) + pp; // not converged: force a retry at higher precision
                break;
            }
            sum = Add(sum, term, pp);
            prevTerm = term;
            pw = Multiply(pw, invX2, pp);
        }
        return (Subtract(sum, lnP, pp), maxLog);
    }

    private static int BitLen(long v) => 64 - System.Numerics.BitOperations.LeadingZeroCount((ulong)Math.Abs(v) | 1);

    /// <summary>log |Γ(x)| at precision p for finite non-pole x (reflection for x &lt; 0).</summary>
    private static (Real Value, long MaxLog2) LnAbsGamma(Real x, int p)
    {
        if (!x.IsNeg)
        {
            return LnGammaPositive(x, p);
        }
        // |Γ(x)| = π / (|sin(πx)| Γ(1-x))
        Real oneMinus = Subtract(One, x, p + 8);
        (Real lg, long ml) = LnGammaPositive(oneMinus, p + 8);
        Real s = SinPi(x, p + 8).Abs();
        Real v = Subtract(Subtract(Log(Pi(p + 8), p + 8), Log(s, p + 8), p + 8), lg, p + 8);
        return (v, Math.Max(ml, Math.Max(Log2Of(lg), 2)));
    }

    private static bool IsNonPositiveInteger(Real x) => x._kind == Kind.Zero || (x._kind == Kind.Normal && x.IsNeg && x.IsInteger);

    /// <summary>Sign of Γ(x) for finite non-pole x: positive for x &gt; 0, (-1)^(floor(-x)+1) for x &lt; 0.</summary>
    private static bool GammaIsNegative(Real x)
    {
        if (!x.IsNeg)
        {
            return false;
        }
        BigInt fl = x.Abs().Floor().ToBigInt();
        return fl.IsEven;
    }

    /// <summary>Gamma function Γ(x) (mpfr_gamma).</summary>
    public static Real Gamma(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = precision != 0 ? precision : x.Precision;
        switch (x._kind)
        {
            case Kind.NaN:
                return MakeNaN(prec);
            case Kind.Inf:
                return x.IsNeg ? MakeNaN(prec) : MakeInf(negative: false, prec);
            case Kind.Zero:
                return MakeInf(x.IsNeg, prec);
        }
        if (x.IsInteger)
        {
            if (x.IsNeg)
            {
                return MakeNaN(prec);
            }
            BigInt n = x.ToBigInt();
            if (n.GetBitLength() <= 17)
            {
                return FromBigInt(BigInt.Factorial((ulong)n - 1), prec, rounding);
            }
        }
        if (x.BinaryExponent >= 60)
        {
            return MakeInf(negative: false, prec); // astronomically large
        }
        bool neg = GammaIsNegative(x);
        return Ziv(prec, rounding, wp =>
        {
            // Γ = exp(ln|Γ|): the absolute error of ln|Γ| becomes the relative error, so
            // carry log2 |ln Γ| extra bits (|ln Γ| ~ x ln x for large x, -ln x for tiny x).
            double ax = Math.Abs(x.ToDouble());
            double lgMag = Math.Max(ax * Math.Log(ax + 2), Math.Abs(x.BinaryExponent) * 0.7) + 2;
            int extra = (int)Math.Log2(lgMag) + 4;
            return Guarded(wp, 64 + extra, p =>
            {
                (Real lg, long ml) = LnAbsGamma(x, p);
                Real g = Exp(lg, p);
                return (neg ? -g : g, Log2Of(g) + Math.Max(0, Math.Max(ml, Log2Of(lg))));
            });
        });
    }

    /// <summary>ln Γ(x) (mpfr_lngamma): NaN where Γ(x) &lt; 0, +inf at poles.</summary>
    public static Real LogGamma(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = precision != 0 ? precision : x.Precision;
        if (x._kind == Kind.NaN)
        {
            return MakeNaN(prec);
        }
        if (x._kind == Kind.Inf || IsNonPositiveInteger(x))
        {
            return MakeInf(negative: false, prec);
        }
        if (GammaIsNegative(x))
        {
            return MakeNaN(prec);
        }
        return LgammaCore(x, prec, rounding);
    }

    /// <summary>log |Γ(x)| and the sign of Γ(x) (mpfr_lgamma).</summary>
    public static Real LogGamma(Real x, out int sign, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = precision != 0 ? precision : x.Precision;
        sign = 1;
        if (x._kind == Kind.NaN)
        {
            return MakeNaN(prec);
        }
        if (x._kind == Kind.Inf)
        {
            sign = 1;
            return MakeInf(negative: false, prec);
        }
        if (IsNonPositiveInteger(x))
        {
            // mpfr: sign -1 only for -0 (Γ(-0) = -inf); +1 at the negative-integer poles.
            sign = x._kind == Kind.Zero && x.IsNeg ? -1 : 1;
            return MakeInf(negative: false, prec);
        }
        sign = GammaIsNegative(x) ? -1 : 1;
        return LgammaCore(x, prec, rounding);
    }

    private static Real LgammaCore(Real x, int prec, RealRounding rounding)
    {
        // lnΓ(1) = lnΓ(2) = +0 exactly.
        if (CompareValues(x, One) == 0 || CompareValues(x, Two) == 0)
        {
            return MakeZero(negative: false, prec);
        }
        if (x.BinaryExponent >= 60)
        {
            return MakeInf(negative: false, prec);
        }
        return Ziv(prec, rounding, wp => Guarded(wp, 64, p => LnAbsGamma(x, p)));
    }

    /// <summary>Digamma ψ(x) = Γ'(x)/Γ(x) (mpfr_digamma).</summary>
    public static Real Digamma(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = precision != 0 ? precision : x.Precision;
        switch (x._kind)
        {
            case Kind.NaN:
                return MakeNaN(prec);
            case Kind.Inf:
                return x.IsNeg ? MakeNaN(prec) : MakeInf(negative: false, prec);
            case Kind.Zero:
                return MakeInf(!x.IsNeg, prec); // ψ(+0) = -inf, ψ(-0) = +inf
        }
        if (x.IsNeg && x.IsInteger)
        {
            return MakeNaN(prec);
        }
        return Ziv(prec, rounding, wp => Guarded(wp, 64, p => DigammaEval(x, p)));
    }

    private static (Real, long) DigammaEval(Real x, int p)
    {
        if (x.IsNeg)
        {
            // ψ(x) = ψ(1 - x) - π cot(πx)
            (Real v1, long m1) = DigammaPositive(Subtract(One, x, p + 8), p);
            Real c = Divide(Multiply(Pi(p), CosPi(x, p), p), SinPi(x, p), p);
            return (Subtract(v1, c, p), Math.Max(m1, Math.Max(Log2Of(v1), Log2Of(c))));
        }
        return DigammaPositive(x, p);
    }

    /// <summary>ψ(x), x &gt; 0: ψ(x) = ψ(X) - Σ_{i&lt;N} 1/(x+i), ψ(X) ~ ln X - 1/(2X) - Σ B_{2k}/(2k X^{2k}).</summary>
    private static (Real, long) DigammaPositive(Real x, int p)
    {
        double xd = Math.Min(x.ToDouble(), 1e300);
        double x0 = Math.Max(10, p * 0.3);
        long n = xd >= x0 ? 0 : (long)Math.Ceiling(x0 - xd);
        int pp = p + BitLen(n + 1) + 8;
        Real bigX = n == 0 ? x : Add(x, Int(n), pp);
        Real recip = MakeZero(negative: false, pp);
        long maxLog = 0;
        for (long i = 0; i < n; i++)
        {
            Real t = Divide(One, Add(x, Int(i), pp), pp);
            maxLog = Math.Max(maxLog, Log2Of(t));
            recip = Add(recip, t, pp);
        }
        maxLog = Math.Max(maxLog, Log2Of(recip));
        Real lnX = Log(bigX, pp);
        maxLog = Math.Max(maxLog, Log2Of(lnX));
        Real invX = Divide(One, bigX, pp), invX2 = Square(invX, pp);
        Real sum = Subtract(lnX, invX.ShiftRight(1), pp);
        Real pw = invX2;
        long tol = Math.Max(Log2Of(sum), -pp) - pp - 8;
        for (int k = 1; k < 1 << 20; k++)
        {
            Real term = Multiply(BernoulliOver2k(k, pp), pw, pp);
            if (term._kind != Kind.Normal || term.BinaryExponent < tol)
            {
                break;
            }
            sum = Subtract(sum, term, pp);
            pw = Multiply(pw, invX2, pp);
        }
        return (Subtract(sum, recip, pp), maxLog);
    }

    /// <summary>Beta function B(a, b) = Γ(a)Γ(b)/Γ(a+b) (mpfr_beta).</summary>
    public static Real Beta(Real a, Real b, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrec(a, b, precision);
        if (a._kind == Kind.NaN || b._kind == Kind.NaN)
        {
            return MakeNaN(prec);
        }
        if (a._kind != Kind.Normal || b._kind != Kind.Normal)
        {
            // Zeros and infinities: follow the limit when it exists.
            if (a._kind == Kind.Inf || b._kind == Kind.Inf)
            {
                return MakeNaN(prec);
            }
            return MakeInf(a._kind == Kind.Zero ? a.IsNeg : b.IsNeg, prec);
        }
        // B(n, b) = (n-1)! / (b (b+1) ... (b+n-1)) exactly for a small positive integer n.
        if (TryBetaInteger(a, b, prec, rounding, out Real r) || TryBetaInteger(b, a, prec, rounding, out r))
        {
            return r;
        }
        Real s = ExactAdd(a, b);
        if (IsNonPositiveInteger(a) || IsNonPositiveInteger(b))
        {
            return MakeNaN(prec); // pole of Γ(a) or Γ(b) with no cancelling pole (mpfr)
        }
        if (IsNonPositiveInteger(s))
        {
            return MakeZero(negative: false, prec);
        }
        bool neg = GammaIsNegative(a) ^ GammaIsNegative(b) ^ GammaIsNegative(s);
        return Ziv(prec, rounding, wp =>
        {
            double mag = Math.Abs(a.ToDouble()) + Math.Abs(b.ToDouble()) + 2;
            int extra = (int)Math.Max(0, Math.Log2(mag * Math.Log(mag))) + 4;
            return Guarded(wp, 64 + extra, p =>
            {
                (Real la, long ma) = LnAbsGamma(a, p);
                (Real lb, long mb) = LnAbsGamma(b, p);
                (Real ls, long ms) = LnAbsGamma(s, p);
                Real l = Subtract(Add(la, lb, p), ls, p);
                Real v = Exp(l, p);
                // Absolute error of l (relative error of v) scales with its largest component.
                long comp = Math.Max(Math.Max(ma, mb), Math.Max(ms, Math.Max(Log2Of(la), Math.Max(Log2Of(lb), Log2Of(ls)))));
                return (neg ? -v : v, Log2Of(v) + Math.Max(0, comp));
            });
        });
    }

    private static bool TryBetaInteger(Real a, Real b, int prec, RealRounding rounding, out Real result)
    {
        result = default;
        if (!a.IsInteger || a.IsNeg || a.BinaryExponent > 12)
        {
            return false;
        }
        long n = (long)a.ToBigInt();
        Real den = b;
        for (long i = 1; i < n; i++)
        {
            Real t = ExactAdd(b, Int(i));
            if (t._kind == Kind.Zero)
            {
                // A zero factor is an uncompensated pole of Γ(b+i) (mpfr → NaN).
                result = MakeNaN(prec);
                return true;
            }
            int need = (int)Math.Min(MaxPrecision, (long)den.Precision + t.Precision + 2);
            den = Multiply(den, t, need); // exact
        }
        result = Divide(FromBigInt(BigInt.Factorial((ulong)(n - 1))), den, prec, rounding);
        return true;
    }

    // ---------------------------------------------------------------
    // Error function
    // ---------------------------------------------------------------

    /// <summary>Error function erf(x) (mpfr_erf).</summary>
    public static Real Erf(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = precision != 0 ? precision : x.Precision;
        switch (x._kind)
        {
            case Kind.NaN:
                return MakeNaN(prec);
            case Kind.Inf:
                return ExactOne(prec, rounding, x.IsNeg);
            case Kind.Zero:
                return MakeZero(x.IsNeg, prec);
        }
        if (x.BinaryExponent > 40)
        {
            return Round((BigInt.One << (prec + 2)) - BigInt.One, -(prec + 2), sticky: true, x.IsNeg, prec, rounding, out _);
        }
        return Ziv(prec, rounding, wp =>
        {
            if (x.BinaryExponent < -(wp / 2) - 16)
            {
                // erf x = (2/√π)(x - x³/3 ...): 2x/√π with a tiny correction toward 0.
                Real v = Divide(x.ShiftLeft(1), Sqrt(Pi(wp + 32), wp + 32), wp + 32);
                return Approx(v, 8);
            }
            return Guarded(wp, 64, p =>
            {
                Real xa = x.Abs();
                double xd = xa.ToDouble();
                Real v;
                if (xd * xd < (p + 16) * 0.6931471805599453)
                {
                    v = ErfSeries(xa, p);
                }
                else
                {
                    v = Subtract(One, ErfcAsymptotic(xa, p, out _), p);
                }
                return (x.IsNeg ? -v : v, Log2Of(v));
            });
        });
    }

    /// <summary>erf(x) for x &gt; 0: (2x/√π) e^{-x²} Σ (2x²)^k / (1·3···(2k+1)) (all terms positive).</summary>
    private static Real ErfSeries(Real x, int p)
    {
        int lg = BitLen(p);
        // Fixed point with the series' integer growth (~ x² log2 e bits) on top of p.
        double xd = x.ToDouble();
        int wpf = p + lg + 16;
        BigInt y = Multiply(Square(x, wpf + 8), Two, wpf + 8).ToFixed(wpf);
        BigInt s = SeriesRect(y, wpf, static t => (ulong)(2 * t + 1), out _);
        Real sr = FromFixed(s, wpf, p + 8, RealRounding.NearestEven);
        Real e = Exp(-Square(x, p + 8), p + 8);
        Real c = Divide(x.ShiftLeft(1), Sqrt(Pi(p + 8), p + 8), p + 8);
        _ = xd;
        return Multiply(Multiply(c, e, p + 8), sr, p);
    }

    /// <summary>
    /// erfc(x) for large x: e^{-x²}/(x√π) Σ (-1)^k (2k-1)!!/(2x²)^k, enveloping; ok reports
    /// whether the terms got below 2^-(p+8) before their minimum.
    /// </summary>
    private static Real ErfcAsymptotic(Real x, int p, out bool ok)
    {
        int pp = p + 16;
        Real inv2x2 = Divide(One, Multiply(Square(x, pp), Two, pp), pp);
        Real term = ExactOne(pp, RealRounding.NearestEven);
        Real sum = term;
        ok = false;
        for (long k = 1; k < 1 << 24; k++)
        {
            Real next = Multiply(Multiply(term, Int(2 * k - 1), pp), inv2x2, pp);
            if (next.BinaryExponent < -pp - 8)
            {
                ok = true;
                break;
            }
            if (next.BinaryExponent >= term.BinaryExponent)
            {
                break; // reached the smallest term
            }
            term = next;
            sum = (k & 1) != 0 ? Subtract(sum, term, pp) : Add(sum, term, pp);
        }
        Real pre = Divide(Exp(-Square(x, pp), pp), Multiply(x, Sqrt(Pi(pp), pp), pp), pp);
        return Multiply(pre, sum, p);
    }

    /// <summary>Complementary error function erfc(x) = 1 - erf(x) (mpfr_erfc).</summary>
    public static Real Erfc(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = precision != 0 ? precision : x.Precision;
        switch (x._kind)
        {
            case Kind.NaN:
                return MakeNaN(prec);
            case Kind.Inf:
                return x.IsNeg ? Round(BigInt.Two, 0, sticky: false, negative: false, prec, rounding, out _) : MakeZero(negative: false, prec);
            case Kind.Zero:
                return ExactOne(prec, rounding);
        }
        if (x.IsNeg && x.BinaryExponent > 40)
        {
            return Round((BigInt.One << (prec + 2)) - BigInt.One, -(prec + 1), sticky: true, negative: false, prec, rounding, out _);
        }
        if (!x.IsNeg && x.BinaryExponent > 30)
        {
            // e^{-x²} far below anything the exponent range can show meaningfully.
            return Round(BigInt.One, long.MinValue / 8, sticky: true, negative: false, prec, rounding, out _);
        }
        return Ziv(prec, rounding, wp =>
        {
            double xd = x.ToDouble();
            int extra = x.IsNeg ? 0 : (int)Math.Min(xd * xd * 1.4427 + 8, 1 << 24);
            return Guarded(wp, 64 + extra, p =>
            {
                if (x.IsNeg)
                {
                    Real v = Add(One, ErfSeriesOrOne(x.Abs(), p), p);
                    return (v, 1);
                }
                if (xd * xd >= (p + 24) * 0.6931471805599453)
                {
                    Real a = ErfcAsymptotic(x, p, out bool ok);
                    if (ok)
                    {
                        return (a, Log2Of(a));
                    }
                }
                Real e = ErfSeries(x, p);
                return (Subtract(One, e, p), 1);
            });
        });
    }

    private static Real ErfSeriesOrOne(Real xa, int p)
    {
        double xd = xa.ToDouble();
        return xd * xd < (p + 16) * 0.6931471805599453 ? ErfSeries(xa, p) : Subtract(One, ErfcAsymptotic(xa, p, out _), p);
    }

    // ---------------------------------------------------------------
    // Riemann zeta function
    // ---------------------------------------------------------------

    /// <summary>Riemann zeta function ζ(s) for real s (mpfr_zeta).</summary>
    public static Real Zeta(Real s, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = precision != 0 ? precision : s.Precision;
        switch (s._kind)
        {
            case Kind.NaN:
                return MakeNaN(prec);
            case Kind.Inf:
                return s.IsNeg ? MakeNaN(prec) : ExactOne(prec, rounding);
            case Kind.Zero:
                return Round(BigInt.One, -1, sticky: false, negative: true, prec, rounding, out _); // -1/2
        }
        int c1 = CompareValues(s, One);
        if (c1 == 0)
        {
            return MakeInf(negative: false, prec);
        }
        if (s.IsNeg && s.IsInteger && s.ToBigInt().IsEven)
        {
            return MakeZero(negative: false, prec); // trivial zeros
        }
        if (!s.IsNeg && s.BinaryExponent > 1 && CompareValues(s, Int(prec + 8)) > 0)
        {
            // ζ(s) = 1 + 2^-s + ... with 2^-s below half an ulp of 1.
            return Round(BigInt.One << (prec + 2), -(prec + 2), sticky: true, negative: false, prec, rounding, out _);
        }
        if (s.IsNeg && s.BinaryExponent > 40)
        {
            return MakeNaN(prec); // beyond the supported range
        }
        return Ziv(prec, rounding, wp => Guarded(wp, 64, p => s.IsNeg ? ZetaReflect(s, p) : ZetaEulerMaclaurin(s, p)));
    }

    /// <summary>ζ(s) for s &lt; 0: 2^s π^(s-1) sin(πs/2) Γ(1-s) ζ(1-s).</summary>
    private static (Real, long) ZetaReflect(Real s, int p)
    {
        int pp = p + 16;
        Real oneMinus = ExactAdd(One, -s);
        (Real z1, _) = ZetaEulerMaclaurin(oneMinus, pp);
        (Real lg, _) = LnGammaPositive(oneMinus, pp);
        Real l = Add(Multiply(s, Ln2(pp), pp), Add(Multiply(Subtract(s, One, pp), Log(Pi(pp), pp), pp), lg, pp), pp);
        Real mag = Exp(l, pp);
        Real v = Multiply(Multiply(mag, SinPi(s.ShiftRight(1), pp), pp), z1, p);
        // exp(l): absolute error of l is relative error of v.
        return (v, Log2Of(v) + Math.Max(0, Log2Of(l)) + Math.Max(0, Log2Of(lg)));
    }

    /// <summary>
    /// EulerGamma-Maclaurin: ζ(s) = Σ_{n&lt;N} n^-s + N^(1-s)/(s-1) + N^-s/2
    /// + Σ_k B_{2k}/(2k)! s(s+1)...(s+2k-2) N^(1-s-2k).
    /// </summary>
    private static (Real, long) ZetaEulerMaclaurin(Real s, int p)
    {
        double sd = s.ToDouble();
        int pp = p + 16;
        long n = Math.Max(4, (long)Math.Ceiling(pp / 6.0 + Math.Abs(sd) / 4.0));
        bool intS = s.IsInteger && s.BinaryExponent <= 30;
        long si = intS ? (long)s.ToBigInt() : 0;
        // n^-s for n < N via the smallest prime factor (exp/log only at primes).
        var pw = new Real[n + 1];
        Real sum = MakeZero(negative: false, pp);
        long maxLog = 0;
        for (long k = 1; k <= n; k++)
        {
            Real v;
            if (k == 1)
            {
                v = One;
            }
            else
            {
                long q = SmallestFactor(k);
                if (q == k)
                {
                    v = intS ? Pow(Int(k), -si, pp) : Exp(Multiply(-s, Log(Int(k), pp), pp), pp);
                }
                else
                {
                    v = Multiply(pw[q], pw[k / q], pp);
                }
            }
            pw[k] = v;
            if (k < n)
            {
                sum = Add(sum, v, pp);
                maxLog = Math.Max(maxLog, Log2Of(v));
            }
        }
        Real nR = Int(n);
        Real nms = pw[n]; // N^-s
        Real t1 = Divide(Multiply(nms, nR, pp), Subtract(s, One, pp), pp);
        Real t2 = nms.ShiftRight(1);
        sum = Add(Add(sum, t1, pp), t2, pp);
        maxLog = Math.Max(maxLog, Math.Max(Log2Of(t1), Log2Of(sum)));
        // Tail: term_k = B_{2k}/(2k)! (s)_{2k-1} N^{1-s-2k}; (s)_{2k-1} N^{1-2k} kept in f.
        Real f = Divide(s, nR, pp); // (s)_1 N^-1
        Real invN2 = Divide(One, Square(nR, pp), pp);
        long tol = Log2Of(sum) - pp - 8;
        BigInt fact = BigInt.Two; // (2k)!
        Real prev = default;
        for (int k = 1; k < 1 << 20; k++)
        {
            if (k > 1)
            {
                f = Multiply(Multiply(Multiply(f, Add(s, Int(2 * k - 3), pp), pp), Add(s, Int(2 * k - 2), pp), pp), invN2, pp);
                fact *= (2 * k - 1) * (2L * k);
            }
            Real b2k = Multiply(BernoulliOver2k(k, pp), Int(2 * k), pp);
            Real term = Multiply(Multiply(Divide(b2k, Int(fact), pp), f, pp), nms, pp);
            if (term._kind != Kind.Normal || term.BinaryExponent < tol)
            {
                break;
            }
            if (k > 2 && prev._kind == Kind.Normal && term.BinaryExponent > prev.BinaryExponent)
            {
                maxLog = Log2Of(sum) + pp; // asymptotic minimum reached: force a retry
                break;
            }
            sum = Add(sum, term, pp);
            prev = term;
        }
        return (sum, maxLog);
    }

    private static long SmallestFactor(long k)
    {
        if ((k & 1) == 0)
        {
            return 2;
        }
        for (long d = 3; d * d <= k; d += 2)
        {
            if (k % d == 0)
            {
                return d;
            }
        }
        return k;
    }

    /// <summary>ζ(n) for an unsigned integer n (mpfr_zeta_ui).</summary>
    public static Real Zeta(ulong n, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
        => Zeta(FromBigInt(new BigInt(n), 66), precision != 0 ? precision : DefaultPrecision, rounding);

    // ---------------------------------------------------------------
    // Exponential integral and dilogarithm
    // ---------------------------------------------------------------

    /// <summary>Exponential integral Ei(x) (mpfr_eint); for x &lt; 0 this is -E1(-x).</summary>
    public static Real ExponentialIntegral(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = precision != 0 ? precision : x.Precision;
        switch (x._kind)
        {
            case Kind.NaN:
                return MakeNaN(prec);
            case Kind.Inf:
                return x.IsNeg ? MakeZero(negative: true, prec) : MakeInf(negative: false, prec);
            case Kind.Zero:
                return MakeInf(negative: true, prec);
        }
        if (x.BinaryExponent > 40)
        {
            return x.IsNeg ? MakeZero(negative: true, prec) : MakeInf(negative: false, prec);
        }
        return Ziv(prec, rounding, wp =>
        {
            double xd = x.ToDouble();
            int extra = x.IsNeg ? (int)Math.Min(-xd * 2.886 + 8, 1 << 24) : 0;
            return Guarded(wp, 64 + extra, p => EintEval(x, p));
        });
    }

    private static (Real, long) EintEval(Real x, int p)
    {
        double xd = x.ToDouble();
        int pp = p + 16;
        if (Math.Abs(xd) >= (pp + 32) * 0.6931471805599453)
        {
            // Asymptotic: Ei(x) ~ e^x/x Σ k!/x^k (enveloping for x < 0, stopped well before the minimum).
            Real invX = Divide(One, x, pp);
            Real term = One, sum = One;
            for (long k = 1; k < 1 << 24; k++)
            {
                Real next = Multiply(Multiply(term, Int(k), pp), invX, pp);
                if (next.BinaryExponent < -pp - 8 || next.BinaryExponent >= term.BinaryExponent)
                {
                    break;
                }
                term = next;
                sum = Add(sum, term, pp);
            }
            Real v = Multiply(Divide(Exp(x, pp), x, pp), sum, p);
            return (v, Log2Of(v));
        }
        // Ei(x) = γ + ln|x| + x Σ_{j>=0} x^j / ((j+1)(j+1)!)
        int wpf = pp + 16;
        BigInt xf = x.ToFixed(wpf);
        BigInt s = SeriesRect(xf, wpf, static t => (ulong)t, static t => (ulong)(t + 1) * (ulong)(t + 1), out _);
        Real ser = Multiply(x, FromFixed(s, wpf, pp, RealRounding.NearestEven), pp);
        Real lg = Add(EulerGamma(pp), Log(x.Abs(), pp), pp);
        Real r = Add(lg, ser, p);
        long maxLog = Math.Max(Log2Of(lg), Math.Max(Log2Of(ser), Math.Abs(xd) > 1 ? (long)(Math.Abs(xd) * 1.4427) + 2 : 1));
        return (r, maxLog);
    }

    /// <summary>Dilogarithm Dilogarithm(x) = -∫_0^x log(1-t)/t dt (real part for x &gt; 1; mpfr_li2).</summary>
    public static Real Dilogarithm(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = precision != 0 ? precision : x.Precision;
        switch (x._kind)
        {
            case Kind.NaN:
                return MakeNaN(prec);
            case Kind.Inf:
                return MakeInf(negative: true, prec);
            case Kind.Zero:
                return MakeZero(x.IsNeg, prec);
        }
        if (x.BinaryExponent > 40)
        {
            // Dilogarithm(x) ~ -ln²|x|/2 (+ π²/3 for x > 1): let the general path handle via Ziv.
        }
        return Ziv(prec, rounding, wp =>
        {
            if (x.BinaryExponent < -(wp / 2) - 16)
            {
                return NearX(x, wp, x.Sign); // Dilogarithm(x) = x + x²/4 + ...
            }
            return Guarded(wp, 64, p => Li2Eval(x, p));
        });
    }

    private static (Real, long) Li2Eval(Real x, int p)
    {
        int pp = p + 16;
        Real pi2over6 = Divide(Square(Pi(pp), pp), Int(6), pp);
        if (x.IsNeg)
        {
            // Landen: Dilogarithm(x) = -Dilogarithm(x/(x-1)) - ln²(1-x)/2, with x/(x-1) in (0, 1).
            Real y = Divide(x, Subtract(x, One, pp), pp);
            (Real ly, long m1) = Li2Unit(y, pp);
            Real l = LogP1(-x, pp);
            Real t = Square(l, pp).ShiftRight(1);
            Real v = Subtract(-ly, t, p);
            return (v, Math.Max(m1, Math.Max(Log2Of(ly), Log2Of(t))));
        }
        int c1 = CompareValues(x, One);
        if (c1 == 0)
        {
            return (pi2over6.WithPrecision(p), Log2Of(pi2over6));
        }
        if (c1 < 0)
        {
            return Li2Unit(x, pp);
        }
        // x > 1: Re Dilogarithm(x) = π²/3 - ln²(x)/2 - Dilogarithm(1/x)
        Real inv = Divide(One, x, pp);
        (Real li, long mi) = Li2Unit(inv, pp);
        Real lx = Log(x, pp);
        Real a = pi2over6.ShiftLeft(1), b = Square(lx, pp).ShiftRight(1);
        Real r = Subtract(Subtract(a, b, pp), li, p);
        return (r, Math.Max(mi, Math.Max(Log2Of(a), Math.Max(Log2Of(b), Log2Of(li)))));
    }

    /// <summary>Dilogarithm(x) for 0 &lt; x &lt; 1 (reflection above 1/2, Bernoulli series in u = -ln(1-x) below).</summary>
    private static (Real, long) Li2Unit(Real x, int pp)
    {
        if (CompareValues(x, Half) > 0)
        {
            // Dilogarithm(x) = π²/6 - ln(x) ln(1-x) - Dilogarithm(1-x)
            Real omx = Subtract(One, x, pp);
            (Real l2, long m2) = Li2Series(omx, pp);
            Real prod = Multiply(Log(x, pp), Log(omx, pp), pp);
            Real pi2over6 = Divide(Square(Pi(pp), pp), Int(6), pp);
            Real v = Subtract(Subtract(pi2over6, prod, pp), l2, pp);
            return (v, Math.Max(m2, Math.Max(Log2Of(pi2over6), Log2Of(prod))));
        }
        return Li2Series(x, pp);
    }

    /// <summary>Dilogarithm(x) = Σ B_n u^(n+1)/(n+1)! with u = -ln(1-x), 0 &lt; x &lt;= 1/2 (u &lt;= ln 2).</summary>
    private static (Real, long) Li2Series(Real x, int pp)
    {
        Real u = -LogP1(-x, pp);
        Real u2 = Square(u, pp);
        // u - u²/4 + Σ_{k>=1} B_{2k} u^{2k+1}/(2k+1)!
        Real sum = Subtract(u, u2.ShiftRight(2), pp);
        Real pw = Multiply(u, u2, pp); // u^3
        BigInt fact = new BigInt(6); // (2k+1)!
        long tol = Log2Of(sum) - pp - 8;
        for (int k = 1; k < 1 << 20; k++)
        {
            if (k > 1)
            {
                fact *= (2 * k) * (2L * k + 1);
            }
            Real b2k = Multiply(BernoulliOver2k(k, pp), Int(2 * k), pp);
            Real term = Divide(Multiply(b2k, pw, pp), Int(fact), pp);
            if (term._kind != Kind.Normal || term.BinaryExponent < tol)
            {
                break;
            }
            sum = Add(sum, term, pp);
            pw = Multiply(pw, u2, pp);
        }
        return (sum, Math.Max(Log2Of(u), Log2Of(sum)));
    }

    // ---------------------------------------------------------------
    // Bessel functions of the first and second kind
    // ---------------------------------------------------------------

    /// <summary>Bessel function of the first kind of order 0 (mpfr_j0).</summary>
    internal static Real J0(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven) => BesselJ(0, x, precision, rounding);

    /// <summary>Bessel function of the first kind of order 1 (mpfr_j1).</summary>
    internal static Real J1(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven) => BesselJ(1, x, precision, rounding);

    /// <summary>Bessel function of the first kind of order n (mpfr_jn).</summary>
    public static Real BesselJ(long n, Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = precision != 0 ? precision : x.Precision;
        switch (x._kind)
        {
            case Kind.NaN:
                return MakeNaN(prec);
            case Kind.Inf:
                return MakeZero(negative: false, prec);
            case Kind.Zero:
                if (n == 0)
                {
                    return ExactOne(prec, rounding);
                }
                // J_n(±0) = ±0 for odd n (sign of x^n), +0 for even n (with the (-1)^n of negative n).
                return MakeZero(((n & 1) != 0 && x.IsNeg) ^ (n < 0 && (n & 1) != 0), prec);
        }
        // J_{-n}(x) = (-1)^n J_n(x); J_n(-x) = (-1)^n J_n(x).
        bool neg = (n < 0 && (n & 1) != 0) ^ (x.IsNeg && (n & 1) != 0);
        long an = Math.Abs(n);
        Real xa = x.Abs();
        if (x.BinaryExponent > 40)
        {
            return MakeNaN(prec); // beyond the supported range
        }
        return Ziv(prec, rounding, wp =>
        {
            if (an == 0 && x.BinaryExponent < -(wp / 2) - 16)
            {
                // J0(x) = 1 - x²/4 + ...: just below 1.
                return ((BigInt.One << (wp + 16)) - BigInt.One, -(wp + 16), 2);
            }
            double xd = xa.ToDouble();
            int extra = (int)Math.Min(xd * 1.45 + 8, 1 << 24);
            return Guarded(wp, 64 + extra, p =>
            {
                (Real v, long m) = BesselJ(an, xa, p);
                return (neg ? -v : v, m);
            });
        });
    }

    /// <summary>J_n(x), x &gt; 0: series (x/2)^n Σ (-x²/4)^k/(k!(n+k)!), or Hankel's expansion for large x.</summary>
    private static (Real, long) BesselJ(long n, Real x, int p)
    {
        int pp = p + 16;
        double xd = x.ToDouble();
        if (xd > (pp + 32) * 0.5 && xd > 4.0 * n * n + 32)
        {
            (Real pa, Real qa) = HankelPQ(n, x, pp);
            Real chi = Subtract(x, Multiply(Pi(pp), Divide(Int(2 * n + 1), Int(4), pp), pp), pp);
            Real amp = Sqrt(Divide(Two, Multiply(Pi(pp), x, pp), pp), pp);
            Real v = Multiply(amp, Subtract(Multiply(pa, Cos(chi, pp), pp), Multiply(qa, Sin(chi, pp), pp), pp), p);
            return (v, Log2Of(amp) + 1);
        }
        int wpf = pp + 16;
        BigInt y = -Square(x.ShiftRight(1), wpf + 8).ToFixed(wpf);
        // Σ y^k / (k! (n+k)!) * n!  : term ratio 1/(k (n+k)).
        BigInt s = SeriesRect(y, wpf, t => (ulong)t * (ulong)(n + t), out _);
        Real ser = FromFixed(s, wpf, pp, RealRounding.NearestEven);
        Real pre = Divide(Pow(x.ShiftRight(1), n, pp), FromBigInt(BigInt.Factorial((ulong)n)), pp);
        Real r = Multiply(pre, ser, p);
        long maxLog = Log2Of(pre) + (long)(xd * 1.4427) + 2;
        return (r, maxLog);
    }

    /// <summary>
    /// Hankel's asymptotic P and Q for order n: P = Σ (-1)^k a_{2k}/x^{2k}, Q = Σ (-1)^k a_{2k+1}/x^{2k+1}
    /// with a_k = Π_{j=1..k} (4n² - (2j-1)²) / (k! 8^k), stopped at the smallest term.
    /// </summary>
    private static (Real P, Real Q) HankelPQ(long n, Real x, int pp)
    {
        Real mu = Int(4 * n * n);
        Real inv8x = Divide(One, Multiply(x, Int(8), pp), pp);
        Real pSum = One, qSum = MakeZero(negative: false, pp);
        Real a = One; // a_k / x^k with signs folded later
        for (long k = 1; k < 1 << 20; k++)
        {
            Real next = Divide(Multiply(Multiply(a, Subtract(mu, Int((2 * k - 1) * (2 * k - 1)), pp), pp), inv8x, pp), Int(k), pp);
            if (next._kind != Kind.Normal || next.BinaryExponent < -pp - 8 || (k > 2 * n + 2 && next.BinaryExponent >= a.BinaryExponent))
            {
                break;
            }
            a = next;
            // k even -> P with sign (-1)^(k/2); k odd -> Q with sign (-1)^((k-1)/2)
            bool negTerm = ((k / 2) & 1) != 0;
            if ((k & 1) == 0)
            {
                pSum = negTerm ? Subtract(pSum, a, pp) : Add(pSum, a, pp);
            }
            else
            {
                qSum = negTerm ? Subtract(qSum, a, pp) : Add(qSum, a, pp);
            }
        }
        return (pSum, qSum);
    }

    /// <summary>Bessel function of the second kind of order 0 (mpfr_y0).</summary>
    internal static Real Y0(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven) => BesselY(0, x, precision, rounding);

    /// <summary>Bessel function of the second kind of order 1 (mpfr_y1).</summary>
    internal static Real Y1(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven) => BesselY(1, x, precision, rounding);

    /// <summary>Bessel function of the second kind of order n (mpfr_yn); defined for x &gt; 0.</summary>
    public static Real BesselY(long n, Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = precision != 0 ? precision : x.Precision;
        switch (x._kind)
        {
            case Kind.NaN:
                return MakeNaN(prec);
            case Kind.Inf:
                return x.IsNeg ? MakeNaN(prec) : MakeZero(negative: false, prec);
            case Kind.Zero:
                // Y_n(±0) = -inf, with (-1)^n for negative n.
                return MakeInf(!(n < 0 && (n & 1) != 0), prec);
        }
        if (x.IsNeg || x.BinaryExponent > 40)
        {
            return MakeNaN(prec);
        }
        bool neg = n < 0 && (n & 1) != 0; // Y_{-n} = (-1)^n Y_n
        long an = Math.Abs(n);
        return Ziv(prec, rounding, wp =>
        {
            double xd = x.ToDouble();
            int extra = (int)Math.Min(xd * 1.45 + 8, 1 << 24);
            return Guarded(wp, 64 + extra, p =>
            {
                (Real v, long m) = BesselY(an, x, p);
                return (neg ? -v : v, m);
            });
        });
    }

    /// <summary>
    /// Y_n(x) = (2/π) J_n(x) ln(x/2) - (1/π) Σ_{k&lt;n} (n-k-1)!/k! (x/2)^(2k-n)
    /// - (1/π) (x/2)^n Σ_k (ψ(k+1) + ψ(n+k+1)) (-x²/4)^k / (k!(n+k)!), or Hankel for large x.
    /// </summary>
    private static (Real, long) BesselY(long n, Real x, int p)
    {
        int pp = p + 16;
        double xd = x.ToDouble();
        if (xd > (pp + 32) * 0.5 && xd > 4.0 * n * n + 32)
        {
            (Real pa, Real qa) = HankelPQ(n, x, pp);
            Real chi = Subtract(x, Multiply(Pi(pp), Divide(Int(2 * n + 1), Int(4), pp), pp), pp);
            Real amp = Sqrt(Divide(Two, Multiply(Pi(pp), x, pp), pp), pp);
            Real v = Multiply(amp, Add(Multiply(pa, Sin(chi, pp), pp), Multiply(qa, Cos(chi, pp), pp), pp), p);
            return (v, Log2Of(amp) + 1);
        }
        int wpi = pp + (int)(xd * 1.45) + BitLen(n) * 2 + 16;
        Real pi = Pi(wpi);
        Real half = x.ShiftRight(1);
        (Real jn, long mj) = BesselJ(n, x, wpi);
        Real t1 = Divide(Multiply(Multiply(jn, Log(half, wpi), wpi), Two, wpi), pi, wpi);
        long maxLog = Math.Max(Log2Of(t1), mj + 2);
        // Finite sum for k < n.
        Real fin = MakeZero(negative: false, wpi);
        if (n > 0)
        {
            Real h2 = Square(half, wpi);
            Real pw = Pow(half, -n, wpi); // (x/2)^(-n)
            for (long k = 0; k < n; k++)
            {
                Real c = Divide(FromBigInt(BigInt.Factorial((ulong)(n - k - 1))), FromBigInt(BigInt.Factorial((ulong)k)), wpi);
                Real term = Multiply(c, pw, wpi);
                maxLog = Math.Max(maxLog, Log2Of(term));
                fin = Add(fin, term, wpi);
                pw = Multiply(pw, h2, wpi);
            }
            fin = Divide(fin, pi, wpi);
        }
        // Infinite sum with digamma weights: ψ(m+1) = -γ + H_m.
        Real gamma = EulerGamma(wpi);
        Real hk = MakeZero(negative: false, wpi);        // H_k
        Real hnk = MakeZero(negative: false, wpi);       // H_{n+k}
        for (long j = 1; j <= n; j++)
        {
            hnk = Add(hnk, Divide(One, Int(j), wpi), wpi);
        }
        Real y = -Square(half, wpi);
        Real coef = Divide(One, FromBigInt(BigInt.Factorial((ulong)n)), wpi); // y^k/(k!(n+k)!)
        Real sum = MakeZero(negative: false, wpi);
        long tol = long.MinValue;
        for (long k = 0; k < 1 << 24; k++)
        {
            if (k > 0)
            {
                coef = Divide(Multiply(coef, y, wpi), Multiply(Int(k), Int(n + k), wpi), wpi);
                hk = Add(hk, Divide(One, Int(k), wpi), wpi);
                hnk = Add(hnk, Divide(One, Int(n + k), wpi), wpi);
            }
            Real w = Subtract(Add(hk, hnk, wpi), gamma.ShiftLeft(1), wpi);
            Real term = Multiply(coef, w, wpi);
            if (term._kind == Kind.Normal)
            {
                maxLog = Math.Max(maxLog, Log2Of(term) + Log2Of(Pow(half, n, 64)));
                if (k > xd + 4 && term.BinaryExponent < tol)
                {
                    break;
                }
            }
            sum = Add(sum, term, wpi);
            tol = Log2Of(sum) - wpi - 8;
        }
        Real t3 = Divide(Multiply(sum, Pow(half, n, wpi), wpi), pi, wpi);
        Real r = Subtract(Subtract(t1, fin, wpi), t3, p);
        return (r, Math.Max(maxLog, Math.Max(Log2Of(fin), Log2Of(t3))));
    }

    // ---------------------------------------------------------------
    // Airy function AiryAi
    // ---------------------------------------------------------------

    /// <summary>Airy function Ai(x) (mpfr_ai).</summary>
    public static Real AiryAi(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = precision != 0 ? precision : x.Precision;
        switch (x._kind)
        {
            case Kind.NaN:
                return MakeNaN(prec);
            case Kind.Inf:
                return MakeZero(negative: false, prec); // AiryAi(±inf) = +0 (mpfr)
        }
        if (x._kind == Kind.Normal && x.BinaryExponent > 20)
        {
            return x.IsNeg ? MakeNaN(prec) : MakeZero(negative: false, prec);
        }
        return Ziv(prec, rounding, wp =>
        {
            double xd = x._kind == Kind.Zero ? 0 : Math.Abs(x.ToDouble());
            int extra = (int)Math.Min(Math.Pow(xd, 1.5) * 1.93 + 8, 1 << 24);
            return Guarded(wp, 64 + extra, p => AiSeries(x, p));
        });
    }

    /// <summary>
    /// AiryAi(x) = c1 f(x) - c2 g(x), f = Σ 3^k (1/3)_k x^{3k}/(3k)!, g = Σ 3^k (2/3)_k x^{3k+1}/(3k+1)!,
    /// c1 = 1/(3^{2/3} Γ(2/3)), c2 = 1/(3^{1/3} Γ(1/3)).
    /// </summary>
    private static (Real, long) AiSeries(Real x, int p)
    {
        int pp = p + 16;
        Real three = Int(3);
        Real c1 = Divide(One, Multiply(Exp(Divide(Multiply(Two, Log(three, pp), pp), three, pp), pp), Gamma(Divide(Two, three, pp), pp), pp), pp);
        Real c2 = Divide(One, Multiply(Exp(Divide(Log(three, pp), three, pp), pp), Gamma(Divide(One, three, pp), pp), pp), pp);
        if (x._kind == Kind.Zero)
        {
            return (c1.WithPrecision(p), Log2Of(c1));
        }
        // f: term ratio x³ / ((3k-1)(3k)) ; g: x³ / ((3k)(3k+1)) (the 3^k (a)_k factors cancel to these).
        Real x3 = Multiply(Square(x, pp), x, pp);
        Real f = One, g = One, tf = One, tg = One;
        long maxLog = 0;
        for (long k = 1; k < 1 << 24; k++)
        {
            tf = Divide(Multiply(tf, x3, pp), Int((3 * k - 1) * (3 * k)), pp);
            tg = Divide(Multiply(tg, x3, pp), Int((3 * k) * (3 * k + 1)), pp);
            if (tf._kind != Kind.Normal || (tf.BinaryExponent < -pp - 8 && tg.BinaryExponent < -pp - 8 && k > 2))
            {
                break;
            }
            maxLog = Math.Max(maxLog, Math.Max(Log2Of(tf), Log2Of(tg)));
            f = Add(f, tf, pp);
            g = Add(g, tg, pp);
        }
        Real a = Multiply(c1, f, pp), b = Multiply(Multiply(c2, x, pp), g, pp);
        Real r = Subtract(a, b, p);
        return (r, Math.Max(maxLog, Math.Max(Log2Of(a), Log2Of(b))));
    }

    // ---------------------------------------------------------------
    // Random numbers (mpfr_urandom, mpfr_nrandom, mpfr_erandom)
    // ---------------------------------------------------------------

    /// <summary>
    /// A uniformly distributed real in [0, 1), correctly rounded from an infinite random bit
    /// stream (mpfr_urandom). May round up to 1 under round-to-nearest or upward rounding.
    /// </summary>
    public static Real Random(System.Random random, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        ArgumentNullException.ThrowIfNull(random);
        int prec = precision != 0 ? precision : DefaultPrecision;
        long exp = 0;
        BigInt m;
        // Skip leading zero words (each with probability 2^-64), then draw the significand.
        while ((m = BigInt.Random(random, 64)).IsZero)
        {
            exp -= 64;
        }
        int lead = (int)(64 - m.GetBitLength());
        exp -= lead;
        BigInt sig = BigInt.Random(random, prec + 64) | (BigInt.One << (prec + 63)); // leading 1 + random continuation
        return Round(sig, exp - (prec + 64), sticky: true, negative: false, prec, rounding, out _);
    }

    /// <summary>A standard normal deviate (Box-Muller), rounded (mpfr_nrandom).</summary>
    public static Real RandomNormal(System.Random random, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = precision != 0 ? precision : DefaultPrecision;
        int p = prec + 32;
        Real u1 = Random(random, p, RealRounding.TowardZero), u2 = Random(random, p, RealRounding.TowardZero);
        if (u1._kind == Kind.Zero)
        {
            u1 = Round(BigInt.One, -p, sticky: false, negative: false, p, RealRounding.NearestEven, out _);
        }
        Real r = Sqrt(Multiply(Log(u1, p), Int(-2), p), p);
        return Multiply(r, CosPi(u2.ShiftLeft(1), p), prec, rounding);
    }

    /// <summary>An exponential deviate with rate 1, rounded (mpfr_erandom).</summary>
    public static Real RandomExponential(System.Random random, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = precision != 0 ? precision : DefaultPrecision;
        Real u = Random(random, prec + 32, RealRounding.TowardZero);
        if (u._kind == Kind.Zero)
        {
            u = Round(BigInt.One, -(prec + 32), sticky: false, negative: false, prec + 32, RealRounding.NearestEven, out _);
        }
        RealRounding flipped = rounding switch
        {
            RealRounding.TowardPositive => RealRounding.TowardNegative,
            RealRounding.TowardNegative => RealRounding.TowardPositive,
            _ => rounding,
        };
        return -Log(u, prec, flipped);
    }
}
