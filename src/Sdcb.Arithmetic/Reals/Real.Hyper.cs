namespace Sdcb.Arithmetic;

public readonly partial struct Real
{
    // ---------------------------------------------------------------
    // Hyperbolic functions (mpfr_sinh/cosh/tanh and inverses)
    // ---------------------------------------------------------------

    /// <summary>Hyperbolic sine (mpfr_sinh).</summary>
    public static Real Sinh(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, x.Precision);
        switch (x._kind)
        {
            case Kind.NaN:
                return MakeNaN(prec);
            case Kind.Inf:
                return MakeInf(x.IsNeg, prec);
            case Kind.Zero:
                return MakeZero(x.IsNeg, prec);
        }
        if (x.BinaryExponent >= 62)
        {
            return MakeInf(x.IsNeg, prec);
        }
        return Ziv(prec, rounding, wp =>
        {
            if (x.BinaryExponent < -(wp / 2) - 16)
            {
                return NearX(x, wp, x.Sign); // sinh x = x + x³/6…
            }
            // sinh = sign(x) (e^|x| - e^-|x|)/2 (e^|x| >= 1 keeps full relative precision in
            // fixed point); cancellation for small x adds |Exp2| bits.
            long extra = Math.Max(0, -x.BinaryExponent) + 8;
            int wpi = wp + (int)Math.Min(extra, 1 << 24) + 32;
            BigInt ef = ExpFixed(x.Abs(), wpi);
            BigInt inv = FixDiv(BigInt.One << wpi, ef, wpi);
            BigInt v = (ef - inv).ShiftRightTruncate(1);
            return (x.IsNeg ? -v : v, (long)-wpi, 1L << 20);
        });
    }

    /// <summary>Hyperbolic cosine (mpfr_cosh).</summary>
    public static Real Cosh(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, x.Precision);
        switch (x._kind)
        {
            case Kind.NaN:
                return MakeNaN(prec);
            case Kind.Inf:
                return MakeInf(negative: false, prec);
            case Kind.Zero:
                return Round(BigInt.One, 0, sticky: false, negative: false, prec, rounding, out _);
        }
        if (x.BinaryExponent >= 62)
        {
            return MakeInf(negative: false, prec);
        }
        return Ziv(prec, rounding, wp =>
        {
            int wpi = wp + 48;
            BigInt ef = ExpFixed(x.Abs(), wpi);
            BigInt inv = FixDiv(BigInt.One << wpi, ef, wpi);
            return ((ef + inv).ShiftRightTruncate(1), (long)-wpi, 1L << 20);
        });
    }

    /// <summary>Hyperbolic tangent (mpfr_tanh).</summary>
    public static Real Tanh(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, x.Precision);
        switch (x._kind)
        {
            case Kind.NaN:
                return MakeNaN(prec);
            case Kind.Inf:
                return Round(BigInt.One, 0, sticky: false, x.IsNeg, prec, rounding, out _);
            case Kind.Zero:
                return MakeZero(x.IsNeg, prec);
        }
        // Saturation: |tanh| -> 1 for large |x|.
        if (x.BinaryExponent >= 40)
        {
            return Round((BigInt.One << 8) - BigInt.One, -8, sticky: true, x.IsNeg, prec, rounding, out _);
        }
        return Ziv(prec, rounding, wp =>
        {
            if (x.BinaryExponent < -(wp / 2) - 16)
            {
                return NearX(x, wp, -x.Sign); // tanh x = x - x³/3…
            }
            // tanh = (e^2x - 1)/(e^2x + 1); cancellation for small x.
            long extra = Math.Max(0, -x.BinaryExponent) + 8;
            int wpi = wp + (int)Math.Min(extra, 1 << 24) + 32;
            if (x.BinaryExponent > 0 && x.ToDouble() * 2 > (wpi + 8) * 0.6931)
            {
                // e^2x overflows the fixed-point window: tanh = ±(1 - tiny)
                BigInt one = BigInt.One << wpi;
                return ((x.IsNeg ? -(one - BigInt.One) : one - BigInt.One), (long)-wpi, 2);
            }
            BigInt e2 = ExpFixed(x.ShiftLeft(1), wpi);
            BigInt one2 = BigInt.One << wpi;
            BigInt v = FixDiv(e2 - one2, e2 + one2, wpi);
            return (v, (long)-wpi, 1L << 20);
        });
    }

    /// <summary>exp(x) in fixed point at wp bits (positive result).</summary>
    private static BigInt ExpFixed(Real x, int wpi)
    {
        (BigInt m, long e, _) = ExpStep(x, wpi);
        long sh = e + wpi;
        return sh >= 0 ? m.ShiftLeft(sh) : m.ShiftRightTruncate(-sh);
    }

    /// <summary>Inverse hyperbolic sine: ln(x + sqrt(x²+1)) (mpfr_asinh).</summary>
    public static Real Asinh(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, x.Precision);
        switch (x._kind)
        {
            case Kind.NaN:
                return MakeNaN(prec);
            case Kind.Inf:
                return MakeInf(x.IsNeg, prec);
            case Kind.Zero:
                return MakeZero(x.IsNeg, prec);
        }
        return Ziv(prec, rounding, wp =>
        {
            if (x.BinaryExponent < -(wp / 2) - 16)
            {
                return NearX(x, wp, -x.Sign); // asinh x = x - x³/6…
            }
            int wpi = wp + 64;
            BigInt one = BigInt.One << wpi;
            BigInt xf = x.Abs().ToFixed(wpi);
            BigInt s = FixSqrt(FixSqr(xf, wpi) + one, wpi);
            Real u = FromFixed(xf + s, wpi, wpi + 8, RealRounding.NearestEven);
            (BigInt lm, long le, long lerr) = LogStep(u, wp + 8);
            return (x.IsNeg ? -lm : lm, le, lerr + 64);
        });
    }

    /// <summary>Inverse hyperbolic cosine: ln(x + sqrt(x²-1)) for x >= 1 (mpfr_acosh).</summary>
    public static Real Acosh(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, x.Precision);
        if (x._kind == Kind.NaN || x._kind == Kind.Zero || x.IsNeg)
        {
            return MakeNaN(prec);
        }
        if (x._kind == Kind.Inf)
        {
            return MakeInf(negative: false, prec);
        }
        int cmp = CompareValues(x, One);
        if (cmp < 0)
        {
            return MakeNaN(prec);
        }
        if (cmp == 0)
        {
            return MakeZero(negative: false, prec); // acosh(1) = +0
        }
        return Ziv(prec, rounding, wp =>
        {
            // Cancellation when x → 1: x²-1 loses ~|Exp2(x-1)| bits.
            long extra = Math.Max(0, -Subtract(x, One, 80).BinaryExponent) + 16;
            int wpi = wp + (int)Math.Min(extra, 1 << 24) + 64;
            BigInt one = BigInt.One << wpi;
            BigInt xf = x.ToFixed(wpi);
            BigInt s = FixSqrt(FixSqr(xf, wpi) - one, wpi);
            Real u = FromFixed(xf + s, wpi, wpi + 8, RealRounding.NearestEven);
            (BigInt lm, long le, long lerr) = LogStep(u, wp + 8);
            return (lm, le, lerr + 64);
        });
    }

    /// <summary>Inverse hyperbolic tangent for |x| &lt; 1 (mpfr_atanh).</summary>
    public static Real Atanh(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
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
            return MakeInf(x.IsNeg, prec);
        }
        return Ziv(prec, rounding, wp =>
        {
            if (x.BinaryExponent < -(wp / 2) - 16)
            {
                return NearX(x, wp, x.Sign); // atanh x = x + x³/3…
            }
            // atanh = ln((1+x)/(1-x))/2; near ±1 needs extra bits.
            long extra = Math.Max(0, -Subtract(One, x.Abs(), 80).BinaryExponent) + 16;
            int wpi = wp + (int)Math.Min(extra, 1 << 24) + 64;
            BigInt one = BigInt.One << wpi;
            BigInt xf = x.Abs().ToFixed(wpi);
            Real u = FromFixed(FixDiv(one + xf, one - xf, wpi), wpi, wpi + 8, RealRounding.NearestEven);
            (BigInt lm, long le, long lerr) = LogStep(u, wp + 8);
            return (x.IsNeg ? -lm : lm, le - 1, lerr + 64);
        });
    }
}
