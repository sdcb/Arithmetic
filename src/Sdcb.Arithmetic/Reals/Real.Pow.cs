using System.Diagnostics;

namespace Sdcb.Arithmetic;

public readonly partial struct Real
{
    // ---------------------------------------------------------------
    // Power functions (mpfr_pow, hypot, cbrt, rootn, rec_sqrt)
    // ---------------------------------------------------------------

    /// <summary>x^y correctly rounded (mpfr_pow), with the full IEEE 754 special-case table.</summary>
    public static Real Pow(Real x, Real y, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, Math.Max(x.Precision, y.Precision));

        // pow(x, ±0) = 1 and pow(+1, y) = 1 even for NaN operands (IEEE 754 / MPFR).
        if (y._kind == Kind.Zero)
        {
            return Round(BigInt.One, 0, sticky: false, negative: false, prec, rounding, out _);
        }
        if (x._kind == Kind.Normal && !x.IsNeg && CompareValues(x, One) == 0)
        {
            return Round(BigInt.One, 0, sticky: false, negative: false, prec, rounding, out _);
        }
        if (x._kind == Kind.NaN || y._kind == Kind.NaN)
        {
            return MakeNaN(prec);
        }

        bool yInt = y.IsInteger;
        bool yOdd = yInt && IsOddInteger(y);

        if (x._kind == Kind.Zero)
        {
            // pow(±0, y): 0 for y > 0 (sign = x's for odd y), inf for y < 0.
            bool neg = yOdd && x.IsNeg;
            return y.IsNeg ? MakeInf(neg, prec) : MakeZero(neg, prec);
        }
        if (x._kind == Kind.Inf)
        {
            bool neg = yOdd && x.IsNeg;
            return y.IsNeg ? MakeZero(neg, prec) : MakeInf(neg, prec);
        }
        if (y._kind == Kind.Inf)
        {
            int cmpAbs = CompareValues(x.Abs(), One);
            if (cmpAbs == 0)
            {
                return Round(BigInt.One, 0, sticky: false, negative: false, prec, rounding, out _); // (±1)^±inf = 1
            }
            bool big = cmpAbs > 0;
            return big == !y.IsNeg ? MakeInf(negative: false, prec) : MakeZero(negative: false, prec);
        }

        if (x.IsNeg && !yInt)
        {
            return MakeNaN(prec);
        }

        // Integer exponent within range: powering with sign handling.
        if (yInt && y.BinaryExponent <= 31)
        {
            long n = (long)y.ToBigInt();
            return Pow(x, n, prec, rounding);
        }
        if (x.IsNeg)
        {
            // Huge odd/even integer exponents: |x|^y with sign from parity.
            Real mag = PowExpLog(x.Abs(), y, prec, rounding, negResult: yOdd);
            return mag;
        }
        // Algebraic exactness: y = n/2^d (dyadic); x^y is representable iff x is a
        // perfect 2^d-th power — detect and dispatch to the exact integer-power path.
        Real? exact = TryPowExact(x, y, prec, rounding);
        if (exact.HasValue)
        {
            return exact.Value;
        }
        return PowExpLog(x, y, prec, rounding, negResult: false);
    }

    /// <summary>Detects exactly-representable x^y for positive x and dyadic non-integer y.</summary>
    private static Real? TryPowExact(Real x, Real y, int prec, RealRounding rounding)
    {
        Debug.Assert(x._kind == Kind.Normal && !x.IsNeg && y._kind == Kind.Normal);
        long tzY = y._man.IsEven ? y._man.TrailingZeroCount() : 0;
        long d = -(y._exp + tzY);
        if (d <= 0)
        {
            return null; // integer y handled elsewhere
        }
        BigInt nOdd = y._man.ShiftRightTruncate(tzY);
        long tzX = x._man.IsEven ? x._man.TrailingZeroCount() : 0;
        BigInt odd = x._man.ShiftRightTruncate(tzX);
        long g = x._exp + tzX;

        if (odd.IsOne)
        {
            // x = 2^g: exact iff 2^d divides g*n (n odd ⇒ iff 2^d | g).
            if (d < 63 && (g & ((1L << (int)d) - 1)) == 0)
            {
                Real b2 = RoundPow2(g >> (int)d, negative: false, MinPrecision + 62, RealRounding.NearestEven);
                return PowIntOfDyadicRoot(b2, nOdd, y.IsNeg, prec, rounding);
            }
            return null;
        }
        if (d > 40 || (g & ((1L << (int)d) - 1)) != 0)
        {
            return null; // perfect 2^d-th powers of odd > 1 need ≥ 2^d+1 bits; g must divide too
        }
        BigInt root = odd;
        for (long i = 0; i < d; i++)
        {
            root = BigInt.Sqrt(root, out BigInt rem);
            if (!rem.IsZero)
            {
                return null;
            }
        }
        int rp = (int)Math.Max(MinPrecision, root.GetBitLength());
        Real rootReal = Round(root, g >> (int)d, sticky: false, negative: false, rp, RealRounding.NearestEven, out _);
        return PowIntOfDyadicRoot(rootReal, nOdd, y.IsNeg, prec, rounding);
    }

    private static Real? PowIntOfDyadicRoot(Real root, BigInt n, bool negativeExp, int prec, RealRounding rounding)
    {
        if (n.GetBitLength() > 40)
        {
            return null; // let Ziv converge; the value is not representable at sane precision
        }
        long nv = (long)n;
        return Pow(root, negativeExp ? -nv : nv, prec, rounding);
    }

    /// <summary>x^n for a machine integer n via squaring under Ziv (mpfr_pow_si).</summary>
    public static Real Pow(Real x, long n, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, x.Precision);
        if (n == 0)
        {
            return Round(BigInt.One, 0, sticky: false, negative: false, prec, rounding, out _);
        }
        switch (x._kind)
        {
            case Kind.NaN:
                return MakeNaN(prec);
            case Kind.Zero:
                {
                    bool neg = (n & 1) != 0 && x.IsNeg;
                    return n < 0 ? MakeInf(neg, prec) : MakeZero(neg, prec);
                }
            case Kind.Inf:
                {
                    bool neg = (n & 1) != 0 && x.IsNeg;
                    return n < 0 ? MakeZero(neg, prec) : MakeInf(neg, prec);
                }
        }
        if (n == 1)
        {
            return Round(x._man, x._exp, sticky: false, x.IsNeg, prec, rounding, out _);
        }
        bool negResult = (n & 1) != 0 && x.IsNeg;

        // Exact powers: ±2^t, or any base whose exact power stays manageable.
        if (x._man.IsPowerOfTwo)
        {
            long t = x.BinaryExponent - 1;
            return RoundPow2(t * n, negResult, prec, rounding);
        }
        if (n > 0 && (long)x.Precision * n <= (1 << 22))
        {
            BigInt em = BigInt.Pow(x._man, (ulong)n);
            return Round(em, x._exp * n, sticky: false, negResult, prec, rounding, out _);
        }
        if (n < 0 && (long)x.Precision * (-n) <= (1 << 22))
        {
            // 1 / exact positive power: a single correctly-rounded division.
            BigInt em = BigInt.Pow(x._man, (ulong)(-n));
            long s = prec + 3 + em.GetBitLength();
            BigInt.DivRem(BigInt.One << (int)Math.Min(s, int.MaxValue), em, out BigInt q, out BigInt r);
            return Round(q, -s - x._exp * (-n), sticky: !r.IsZero, negResult, prec, rounding, out _);
        }
        ulong e = (ulong)Math.Abs(n);

        return Ziv(prec, rounding, wp =>
        {
            int wpi = wp + 64 + 2 * (64 - System.Numerics.BitOperations.LeadingZeroCount(e));
            // Fixed-point square-and-multiply on the mantissa, tracking the binary exponent.
            BigInt one = BigInt.One << wpi;
            BigInt acc = one;
            long accExp = 0;
            long sh0 = wpi - x.Precision;
            BigInt baseF = sh0 >= 0 ? x._man.ShiftLeft(sh0) : x._man.ShiftRightTruncate(-sh0);
            long baseExp = x.BinaryExponent; // value = baseF * 2^(baseExp - wpi)
            ulong k = e;
            while (true)
            {
                if ((k & 1) != 0)
                {
                    acc = FixMul(acc, baseF, wpi);
                    accExp += baseExp;
                    // renormalize acc into [2^(wpi-1), 2^wpi)
                    long bl = acc.GetBitLength();
                    if (bl < wpi)
                    {
                        acc = acc.ShiftLeft(wpi - bl);
                        accExp -= wpi - bl;
                    }
                }
                k >>= 1;
                if (k == 0)
                {
                    break;
                }
                baseF = FixSqr(baseF, wpi);
                baseExp *= 2;
                long bbl = baseF.GetBitLength();
                if (bbl < wpi)
                {
                    baseF = baseF.ShiftLeft(wpi - bbl);
                    baseExp -= wpi - bbl;
                }
            }
            if (n < 0)
            {
                acc = FixDiv(one, acc, wpi);
                accExp = -accExp;
                long bl = acc.GetBitLength();
                if (bl < wpi)
                {
                    acc = acc.ShiftLeft(wpi - bl);
                    accExp -= wpi - bl;
                }
            }
            if (negResult)
            {
                acc = -acc;
            }
            return (acc, accExp - wpi, 128 * (long)(64 - System.Numerics.BitOperations.LeadingZeroCount(e | 1)) + 16);
        });
    }

    /// <summary>Rounds ±2^k to the requested precision.</summary>
    private static Real RoundPow2(long k, bool negative, int prec, RealRounding rounding)
        => Round(BigInt.One, k, sticky: false, negative, prec, rounding, out _);

    /// <summary>General positive-base power via exp(y ln x).</summary>
    private static Real PowExpLog(Real x, Real y, int prec, RealRounding rounding, bool negResult)
    {
        Debug.Assert(x._kind == Kind.Normal && !x.IsNeg);
        return Ziv(prec, rounding, wp =>
        {
            int wpi = wp + 96;
            (BigInt lm, long le, long lerr) = LogStep(x, wpi);
            // z = y * ln x as a Real at high precision.
            Real lnx = lm.IsZero ? MakeZero(negative: false, wpi) : Round(lm.Abs(), le, sticky: true, lm.Sign < 0, wpi, RealRounding.NearestEven, out _);
            Real z = Multiply(y, lnx, wpi, RealRounding.NearestEven);
            if (z._kind != Kind.Normal)
            {
                return (BigInt.One << wp, -wp, 8); // y*lnx ≈ 0 → 1
            }
            if (z.BinaryExponent >= 62)
            {
                // exp overflow/underflow beyond representable exponents.
                return z.IsNeg ? (BigInt.One, EminUlp, 2) : (BigInt.One << 62, EmaxUlp, 2);
            }
            (BigInt em, long ee, long eerr) = ExpStep(z, wp + 16);
            if (negResult)
            {
                em = -em;
            }
            return (em, ee, eerr + lerr + 1024);
        });
    }

    /// <summary>sqrt(x² + y²) without intermediate overflow (mpfr_hypot).</summary>
    public static Real Hypot(Real x, Real y, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, Math.Max(x.Precision, y.Precision));
        if (x._kind == Kind.Inf || y._kind == Kind.Inf)
        {
            return MakeInf(negative: false, prec); // hypot(±inf, NaN) = +inf per IEEE
        }
        if (x._kind == Kind.NaN || y._kind == Kind.NaN)
        {
            return MakeNaN(prec);
        }
        if (x._kind == Kind.Zero)
        {
            return y._kind == Kind.Zero
                ? MakeZero(negative: false, prec)
                : Round(y._man, y._exp, sticky: false, negative: false, prec, rounding, out _);
        }
        if (y._kind == Kind.Zero)
        {
            return Round(x._man, x._exp, sticky: false, negative: false, prec, rounding, out _);
        }
        long diffE = Math.Abs(x.BinaryExponent - y.BinaryExponent);
        if (diffE <= 2L * (x.Precision + y.Precision) + 64)
        {
            return HypotExact(x, y, prec, rounding);
        }
        return Ziv(prec, rounding, wp =>
        {
            // Scale by the larger exponent to keep fixed point well-ranged.
            long scale = Math.Max(x.BinaryExponent, y.BinaryExponent);
            int wpi = wp + 64;
            BigInt xf = x.Abs().ShiftLeft(-scale).ToFixed(wpi);
            BigInt yf = y.Abs().ShiftLeft(-scale).ToFixed(wpi);
            BigInt s = FixSqrt(FixSqr(xf, wpi) + FixSqr(yf, wpi), wpi);
            return (s, -wpi + scale, 16);
        });
    }

    /// <summary>Exact-path hypot: correctly rounded from the exact sum of squares.</summary>
    private static Real HypotExact(Real x, Real y, int prec, RealRounding rounding)
    {
        long e0 = 2 * Math.Min(x._exp, y._exp);
        BigInt sum = x._man.Square().ShiftLeft(2 * x._exp - e0) + y._man.Square().ShiftLeft(2 * y._exp - e0);
        long need = 2L * (prec + 2) - sum.GetBitLength();
        long s = Math.Max(0, need);
        if (((e0 - s) & 1) != 0)
        {
            s++;
        }
        BigInt r = BigInt.Sqrt(sum.ShiftLeft(s), out BigInt rem);
        return Round(r, (e0 - s) / 2, sticky: !rem.IsZero, negative: false, prec, rounding, out _);
    }

    /// <summary>Cube root, defined for negative arguments (mpfr_cbrt).</summary>
    public static Real Cbrt(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
        => RootN(x, 3, precision, rounding);

    /// <summary>k-th root (mpfr_rootn_ui): negative x allowed for odd k.</summary>
    public static Real RootN(Real x, uint k, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        ArgumentOutOfRangeException.ThrowIfZero(k);
        int prec = ResolvePrecision(precision, x.Precision);
        if (k == 1)
        {
            return x.WithPrecision(prec, rounding);
        }
        bool odd = (k & 1) != 0;
        switch (x._kind)
        {
            case Kind.NaN:
                return MakeNaN(prec);
            case Kind.Zero:
                return MakeZero(x.IsNeg && odd, prec);
            case Kind.Inf:
                if (x.IsNeg && !odd)
                {
                    return MakeNaN(prec);
                }
                return MakeInf(x.IsNeg && odd, prec);
        }
        if (x.IsNeg && !odd)
        {
            return MakeNaN(prec);
        }
        if (k == 2)
        {
            return Sqrt(x, prec, rounding); // after the zero case: rootn(-0, 2) = +0, unlike sqrt(-0) = -0
        }
        // Direct correctly-rounded root: r = floor((m << s)^(1/k)), sticky from the remainder.
        long need = (long)k * (prec + 2);
        long s = Math.Max(0, need - x._man.GetBitLength());
        long e = x._exp - s;
        long fix = ((e % k) + k) % k;
        s += fix; // make (x.exp - s) divisible by k
        e = x._exp - s;
        Debug.Assert(e % k == 0);
        BigInt r = BigInt.Root(x._man.ShiftLeft(s), k, out BigInt rem);
        return Round(r, e / k, sticky: !rem.IsZero, x.IsNeg, prec, rounding, out _);
    }

    /// <summary>1/sqrt(x) correctly rounded (mpfr_rec_sqrt).</summary>
    public static Real ReciprocalSqrt(Real x, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, x.Precision);
        switch (x._kind)
        {
            case Kind.NaN:
                return MakeNaN(prec);
            case Kind.Zero:
                return MakeInf(negative: false, prec); // mpfr: rec_sqrt(±0) = +inf
            case Kind.Inf:
                return x.IsNeg ? MakeNaN(prec) : MakeZero(negative: false, prec);
        }
        if (x.IsNeg)
        {
            return MakeNaN(prec);
        }
        // 1/sqrt(x) is dyadic only for x = 4^k, where it is exactly 2^-k.
        if (IsPowerOfTwoValue(x, out long p2) && (p2 & 1) == 0)
        {
            return Round(BigInt.One, -(p2 / 2), sticky: false, negative: false, prec, rounding, out _);
        }
        return Ziv(prec, rounding, wp =>
        {
            int wpi = wp + 64;
            // Scale x into [1, 4) by an even exponent shift.
            long k = x.BinaryExponent - 1;
            if ((k & 1) != 0)
            {
                k--;
            }
            BigInt xf = x.ShiftLeft(-k).ToFixed(2 * wpi);
            BigInt s = BigInt.Sqrt(xf); // sqrt of value*2^(2wpi) = sqrt(value)*2^wpi
            BigInt inv = FixDiv(BigInt.One << wpi, s, wpi);
            return (inv, -wpi - k / 2, 16);
        });
    }
}
