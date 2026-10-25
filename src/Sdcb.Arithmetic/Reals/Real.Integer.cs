namespace Sdcb.Arithmetic;

public readonly partial struct Real
{
    // ---------------------------------------------------------------
    // Rounding to integer values (mpfr_rint family)
    // ---------------------------------------------------------------

    /// <summary>Rounds to an integer value using the given mode (mpfr_rint).</summary>
    public Real RoundToInteger(RealRounding rounding)
    {
        if (_kind != Kind.Normal)
        {
            return this;
        }
        if (_exp >= 0)
        {
            return this; // already an integer
        }
        long frac = -_exp;
        long intBits = Precision - frac; // may be <= 0 for |x| < 1
        BigInt ip = intBits > 0 ? _man.ShiftRightTruncate(frac) : BigInt.Zero;
        bool guard = _man.TestBit(frac - 1);
        bool sticky = frac > 1 && _man.TrailingZeroCount() < frac - 1;
        if (!guard && !sticky)
        {
            // Exact integer already.
            return FromIntegerMagnitude(ip, IsNeg);
        }
        bool up = rounding switch
        {
            RealRounding.NearestEven => guard && (sticky || ip.IsOdd),
            RealRounding.TowardZero => false,
            RealRounding.AwayFromZero => true,
            RealRounding.TowardPositive => !IsNeg,
            RealRounding.TowardNegative => IsNeg,
            _ => throw new ArgumentOutOfRangeException(nameof(rounding)),
        };
        if (up)
        {
            ip += BigInt.One;
        }
        return FromIntegerMagnitude(ip, IsNeg);
    }

    private Real FromIntegerMagnitude(BigInt mag, bool neg)
    {
        if (mag.IsZero)
        {
            return MakeZero(neg, Precision);
        }
        // Keep this value's precision; integers within precision stay exact.
        return Round(mag, 0, sticky: false, neg, Precision, RealRounding.NearestEven, out _);
    }

    /// <summary>Largest integer value not above this (mpfr_floor).</summary>
    public Real Floor() => RoundToInteger(RealRounding.TowardNegative);

    /// <summary>Smallest integer value not below this (mpfr_ceil).</summary>
    public Real Ceiling() => RoundToInteger(RealRounding.TowardPositive);

    /// <summary>Integer value truncated toward zero (mpfr_trunc).</summary>
    public Real Truncate() => RoundToInteger(RealRounding.TowardZero);

    /// <summary>Nearest integer, ties away from zero (mpfr_round).</summary>
    internal Real RoundHalfAway()
    {
        if (_kind != Kind.Normal || _exp >= 0)
        {
            return this;
        }
        long frac = -_exp;
        BigInt ip = Precision - frac > 0 ? _man.ShiftRightTruncate(frac) : BigInt.Zero;
        if (_man.TestBit(frac - 1))
        {
            ip += BigInt.One; // frac >= 0.5: away from zero
        }
        return FromIntegerMagnitude(ip, IsNeg);
    }

    /// <summary>Nearest integer, ties to even (mpfr_roundeven).</summary>
    internal Real RoundHalfEven() => RoundToInteger(RealRounding.NearestEven);

    /// <summary>Fractional part with the sign of the value (mpfr_frac): x - trunc(x).</summary>
    public Real FractionalPart(int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        if (_kind != Kind.Normal)
        {
            return _kind == Kind.Inf ? MakeNaN(Precision) : this;
        }
        if (_exp >= 0)
        {
            return MakeZero(IsNeg, precision != 0 ? precision : Precision);
        }
        long frac = Math.Min(-_exp, Precision); // fractional bits present in the mantissa
        BigInt mask = (BigInt.One << (int)frac) - BigInt.One;
        BigInt fm = _man & mask;
        if (fm.IsZero)
        {
            return MakeZero(IsNeg, precision != 0 ? precision : Precision);
        }
        return Round(fm, _exp, sticky: false, IsNeg, precision != 0 ? precision : Precision, rounding, out _);
    }

    /// <summary>True when the value is an exact integer (mpfr_integer_p).</summary>
    public bool IsInteger
    {
        get
        {
            if (_kind == Kind.Zero)
            {
                return true;
            }
            if (_kind != Kind.Normal)
            {
                return false;
            }
            return _exp >= 0 || _man.TrailingZeroCount() >= -_exp;
        }
    }

    // ---------------------------------------------------------------
    // fmod / remainder (mpfr_fmod, mpfr_remainder)
    // ---------------------------------------------------------------

    /// <summary>Exact x - trunc(x/y)*y with the sign of x (mpfr_fmod), then rounded.</summary>
    public static Real TruncatedRemainder(Real x, Real y, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrec(x, y, precision);
        if (x._kind == Kind.NaN || y._kind == Kind.NaN || x._kind == Kind.Inf || y._kind == Kind.Zero)
        {
            return MakeNaN(prec);
        }
        if (x._kind == Kind.Zero || y._kind == Kind.Inf)
        {
            return x._kind == Kind.Zero ? MakeZero(x.IsNeg, prec) : Round(x._man, x._exp, sticky: false, x.IsNeg, prec, rounding, out _);
        }
        BigInt r = ModMagnitude(x, y, out _);
        if (r.IsZero)
        {
            return MakeZero(x.IsNeg, prec);
        }
        long e = Math.Min(x._exp, y._exp);
        return Round(r, e, sticky: false, x.IsNeg, prec, rounding, out _);
    }

    /// <summary>IEEE remainder: x - n*y with n = nearest integer of x/y, ties to even (mpfr_remainder).</summary>
    public static Real Ieee754Remainder(Real x, Real y, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrec(x, y, precision);
        if (x._kind == Kind.NaN || y._kind == Kind.NaN || x._kind == Kind.Inf || y._kind == Kind.Zero)
        {
            return MakeNaN(prec);
        }
        if (x._kind == Kind.Zero || y._kind == Kind.Inf)
        {
            return x._kind == Kind.Zero ? MakeZero(x.IsNeg, prec) : Round(x._man, x._exp, sticky: false, x.IsNeg, prec, rounding, out _);
        }
        BigInt r = ModMagnitude(x, y, out BigInt yAligned);
        long e = Math.Min(x._exp, y._exp);
        // r in [0, |y|); adjust to (-|y|/2, |y|/2] with ties toward even quotient.
        BigInt twice = r << 1;
        int cmp = twice.CompareTo(yAligned);
        bool subtract = cmp > 0 || (cmp == 0 && IsQuotientOdd(x, y, e, r, yAligned));
        bool neg = x.IsNeg;
        if (subtract)
        {
            r = yAligned - r;
            neg = !neg;
        }
        if (r.IsZero)
        {
            return MakeZero(x.IsNeg, prec);
        }
        return Round(r, e, sticky: false, neg, prec, rounding, out _);
    }

    private static bool IsQuotientOdd(Real x, Real y, long e, BigInt r, BigInt yAligned)
    {
        // q = (|x| - r) / |y| at scale e; oddness of the integer quotient.
        BigInt xa = x._man.ShiftLeft(x._exp - e);
        BigInt q = (xa - r) / yAligned;
        return q.IsOdd;
    }

    /// <summary>|x| mod |y| computed exactly on aligned mantissas; also returns |y| at the common scale.</summary>
    private static BigInt ModMagnitude(Real x, Real y, out BigInt yAligned)
    {
        long e = Math.Min(x._exp, y._exp);
        // xa = |x| / 2^e; ya = |y| / 2^e. xa can be huge when x.exp >> e; use modular reduction.
        yAligned = y._man.ShiftLeft(y._exp - e);
        long shift = x._exp - e;
        // (man * 2^shift) mod ya  — use ModPow for the power of two to keep it cheap.
        if (shift == 0)
        {
            return BigInt.Mod(x._man, yAligned);
        }
        BigInt p2 = BigInt.ModPow(BigInt.Two, new BigInt(shift), yAligned);
        return BigInt.Mod(x._man * p2, yAligned);
    }

    // ---------------------------------------------------------------
    // nextafter (mpfr_nextabove / mpfr_nextbelow)
    // ---------------------------------------------------------------

    /// <summary>The next representable value above this at its precision (mpfr_nextabove).</summary>
    internal Real NextAbove()
    {
        switch (_kind)
        {
            case Kind.NaN:
                return this;
            case Kind.Inf:
                return IsNeg ? LargestNegative() : this;
            case Kind.Zero:
                return SmallestPositive(false);
        }
        if (!IsNeg)
        {
            BigInt m = _man + BigInt.One;
            long e = _exp;
            if (m.GetBitLength() > Precision)
            {
                m = m.ShiftRightTruncate(1);
                e++;
            }
            return new Real(m, e, _prec, Kind.Normal, negative: false);
        }
        else
        {
            BigInt m = _man - BigInt.One;
            long e = _exp;
            if (m.IsZero)
            {
                return MakeZero(negative: true, Precision);
            }
            if (m.GetBitLength() < Precision)
            {
                m = (m << 1) + BigInt.One;
                e--;
            }
            return new Real(m, e, _prec, Kind.Normal, negative: true);
        }
    }

    /// <summary>The next representable value below this (mpfr_nextbelow).</summary>
    internal Real NextBelow() => -(-this).NextAbove();

    private Real SmallestPositive(bool neg)
    {
        int p = Precision;
        BigInt m = BigInt.One << (p - 1);
        return new Real(m, EminUlp, p, Kind.Normal, neg);
    }

    private Real LargestNegative()
    {
        int p = Precision;
        BigInt m = (BigInt.One << p) - BigInt.One;
        return new Real(m, EmaxUlp - p, p, Kind.Normal, negative: true);
    }

    /// <summary>Conventional exponent bounds used by nextafter at the edges (exponents are otherwise unbounded).</summary>
    internal const long EminUlp = long.MinValue / 4;
    internal const long EmaxUlp = long.MaxValue / 4;
}
