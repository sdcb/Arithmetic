using System.Diagnostics;

namespace Sdcb.Arithmetic;

/// <summary>Rounding modes for <see cref="Real"/> operations (MPFR's mpfr_rnd_t).</summary>
public enum RealRounding
{
    /// <summary>Round to nearest, ties to even (RNDN, IEEE default).</summary>
    NearestEven = 0,
    /// <summary>Round toward zero (RNDZ).</summary>
    TowardZero = 1,
    /// <summary>Round toward +infinity (RNDU).</summary>
    TowardPositive = 2,
    /// <summary>Round toward -infinity (RNDD).</summary>
    TowardNegative = 3,
    /// <summary>Round away from zero (RNDA).</summary>
    AwayFromZero = 4,
}

/// <summary>
/// An immutable correctly-rounded arbitrary-precision binary floating-point number,
/// the equivalent of MPFR's mpfr_t: every operation returns the representable value
/// nearest to the exact result under the requested rounding mode, with full support
/// for NaN, ±infinity and signed zero.
/// The significand of a finite nonzero value always has exactly Precision bits
/// (top bit set): value = Sign * Mantissa * 2^Exponent.
/// </summary>
public readonly partial struct Real :
    IEquatable<Real>,
    IComparable<Real>,
    IComparable,
    ISpanFormattable,
    ISpanParsable<Real>
{
    internal const int MinPrecision = 2;
    internal const int MaxPrecision = int.MaxValue / 64;

    private static int s_defaultPrecision = 128;

    /// <summary>Default precision in bits for values created without an explicit precision (mpfr_set_default_prec).</summary>
    public static int DefaultPrecision
    {
        get => s_defaultPrecision;
        set
        {
            ArgumentOutOfRangeException.ThrowIfLessThan(value, MinPrecision);
            ArgumentOutOfRangeException.ThrowIfGreaterThan(value, MaxPrecision);
            s_defaultPrecision = value;
        }
    }

    internal enum Kind : sbyte
    {
        Zero = 0,
        Normal = 1,
        Inf = 2,
        NaN = 3,
    }

    private readonly BigInt _man;  // positive significand with exactly _prec bits (Kind.Normal only)
    private readonly long _exp;    // value = sign * _man * 2^_exp
    private readonly int _prec;    // 0 means "DefaultPrecision at the time of use"
    private readonly Kind _kind;
    private readonly sbyte _neg;   // 1 when the sign bit is set

    internal Real(BigInt man, long exp, int prec, Kind kind, bool negative)
    {
        Debug.Assert(kind != Kind.Normal || (man.Sign > 0 && man.GetBitLength() == prec));
        _man = man;
        _exp = exp;
        _prec = prec;
        _kind = kind;
        _neg = negative ? (sbyte)1 : (sbyte)0;
    }

    // ---------------------------------------------------------------
    // Basic properties
    // ---------------------------------------------------------------

    /// <summary>The working precision of this value in bits.</summary>
    public int Precision => _prec == 0 ? DefaultPrecision : _prec;

    internal Kind KindOf => _kind;
    internal bool IsNeg => _neg != 0;

    /// <summary>The positive significand (exactly Precision bits); undefined for non-normal values.</summary>
    public BigInt Mantissa => _man;

    /// <summary>The exponent of the least significant significand bit: value = Sign * Mantissa * 2^Exponent.</summary>
    public long Exponent => _exp;

    public bool IsNaN => _kind == Kind.NaN;
    public bool IsInfinity => _kind == Kind.Inf;
    public bool IsPositiveInfinity => _kind == Kind.Inf && !IsNeg;
    public bool IsNegativeInfinity => _kind == Kind.Inf && IsNeg;
    public bool IsZero => _kind == Kind.Zero;
    public bool IsFinite => _kind is Kind.Zero or Kind.Normal;

    /// <summary>True for finite nonzero values.</summary>
    public bool IsNumber => _kind == Kind.Normal;

    /// <summary>-1, 0, +1 by value sign (0 for both zeros; throws for NaN like mpfr_sgn semantics).</summary>
    public int Sign
    {
        get
        {
            if (_kind == Kind.NaN)
            {
                throw new ArithmeticException("Sign of NaN is undefined.");
            }
            return _kind == Kind.Zero ? 0 : IsNeg ? -1 : 1;
        }
    }

    /// <summary>The IEEE sign bit (true for -0, negative values and -inf; mpfr_signbit).</summary>
    public bool SignBit => IsNeg;

    /// <summary>Exponent of the most significant bit + 1, i.e. value in [2^(Exp2-1), 2^Exp2) (mpfr_get_exp).</summary>
    public long BinaryExponent => _exp + Precision;

    // ---------------------------------------------------------------
    // Special values
    // ---------------------------------------------------------------

    public static Real NaN => new(default, 0, 0, Kind.NaN, negative: false);
    public static Real PositiveInfinity => new(default, 0, 0, Kind.Inf, negative: false);
    public static Real NegativeInfinity => new(default, 0, 0, Kind.Inf, negative: true);
    public static Real Zero => default;
    public static Real NegativeZero => new(default, 0, 0, Kind.Zero, negative: true);
    public static Real One => FromBigInt(BigInt.One, 0);
    public static Real Two => FromBigInt(BigInt.Two, 0);
    public static Real MinusOne => FromBigInt(BigInt.MinusOne, 0);
    public static Real Half => FromBigInt(BigInt.One, 0).ShiftRight(1);

    /// <summary>An explicit precision, validated to [MinPrecision, MaxPrecision]; 0 selects <paramref name="fallback"/>.</summary>
    internal static int ResolvePrecision(int precision, int fallback)
    {
        if (precision == 0)
        {
            return fallback;
        }
        ValidatePrecision(precision);
        return precision;
    }

    private static void ValidatePrecision(int precision)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(precision, MinPrecision);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(precision, MaxPrecision);
    }

    internal static Real MakeZero(bool negative, int prec) => new(default, 0, prec, Kind.Zero, negative);
    internal static Real MakeInf(bool negative, int prec) => new(default, 0, prec, Kind.Inf, negative);
    internal static Real MakeNaN(int prec) => new(default, 0, prec, Kind.NaN, negative: false);

    // ---------------------------------------------------------------
    // The rounding core
    // ---------------------------------------------------------------

    /// <summary>
    /// Builds a Real from sign * (m + θ) * 2^e where m &gt; 0 and θ ∈ (0,1) iff
    /// sticky, correctly rounded to prec bits. ternary is the sign of
    /// (rounded - exact), 0 when exact (the MPFR ternary value).
    /// </summary>
    internal static Real Round(BigInt m, long e, bool sticky, bool negative, int prec, RealRounding rnd, out int ternary)
    {
        Debug.Assert(m.Sign > 0);
        if (prec == 0)
        {
            prec = DefaultPrecision;
        }
        return RoundLimbs(m.Limbs, e, sticky, negative, prec, rnd, out ternary);
    }

    internal static Real Round(BigInt m, long e, bool sticky, bool negative, int prec, RealRounding rnd)
        => Round(m, e, sticky, negative, prec, rnd, out _);

    /// <summary>The rounding core on a 128-bit significand (prec &lt;= 126): one allocation for the result.</summary>
    internal static Real RoundSmall(UInt128 m, long e, bool sticky, bool negative, int prec, RealRounding rnd, out int ternary)
    {
        Debug.Assert(m != 0 && prec >= MinPrecision && prec <= 126);
        int bits = 128 - (int)UInt128.LeadingZeroCount(m);
        int drop = bits - prec;
        bool guard = false, stickyLow = sticky;
        if (drop <= 0)
        {
            m <<= -drop;
        }
        else
        {
            guard = ((m >> (drop - 1)) & 1) != 0;
            stickyLow |= drop > 1 && (m & ((UInt128.One << (drop - 1)) - 1)) != 0;
            m >>= drop;
        }
        e = checked(e + drop); // an exponent beyond long range must not wrap
        ternary = 0;
        if (guard || stickyLow)
        {
            bool roundUp = rnd switch
            {
                RealRounding.NearestEven => guard && (stickyLow || (m & 1) != 0),
                RealRounding.TowardZero => false,
                RealRounding.AwayFromZero => true,
                RealRounding.TowardPositive => !negative,
                RealRounding.TowardNegative => negative,
                _ => throw new ArgumentOutOfRangeException(nameof(rnd)),
            };
            if (roundUp)
            {
                m++;
                if (m >> prec != 0)
                {
                    m >>= 1; // was 2^prec: exact halving
                    e = checked(e + 1);
                }
                ternary = negative ? -1 : 1;
            }
            else
            {
                ternary = negative ? 1 : -1;
            }
        }
        return new Real(new BigInt(m), e, prec, Kind.Normal, negative);
    }

    /// <summary>Rounds a signed exact value man * 2^e to prec bits.</summary>
    internal static Real RoundSigned(BigInt man, long e, int prec, RealRounding rnd, out int ternary)
    {
        if (man.IsZero)
        {
            ternary = 0;
            return MakeZero(negative: false, prec);
        }
        return Round(man.Abs(), e, sticky: false, man.Sign < 0, prec, rnd, out ternary);
    }

    // ---------------------------------------------------------------
    // Construction
    // ---------------------------------------------------------------

    /// <summary>Creates an exact Real from an integer (rounded only if it exceeds the precision).</summary>
    public static Real FromBigInt(BigInt value, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        // Exact by default: use enough precision for the whole integer.
        int prec = ResolvePrecision(precision, (int)Math.Min(Math.Max(DefaultPrecision, value.GetBitLength()), MaxPrecision));
        if (value.IsZero)
        {
            return MakeZero(negative: false, precision);
        }
        return Round(value.Abs(), 0, sticky: false, value.Sign < 0, prec, rounding, out _);
    }

    /// <summary>Creates an exact Real from an integer (precision is max(DefaultPrecision, bit length)).</summary>
    public Real(long value) : this()
    {
        this = FromBigInt(new BigInt(value));
    }

    public Real(double value, int precision = 0) : this()
    {
        this = FromDouble(value, precision);
    }

    /// <summary>Converts a double exactly (a double always fits in 53 bits).</summary>
    public static Real FromDouble(double value, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = ResolvePrecision(precision, DefaultPrecision);
        if (double.IsNaN(value))
        {
            return MakeNaN(prec);
        }
        if (double.IsPositiveInfinity(value))
        {
            return MakeInf(negative: false, prec);
        }
        if (double.IsNegativeInfinity(value))
        {
            return MakeInf(negative: true, prec);
        }
        if (value == 0)
        {
            return MakeZero(double.IsNegative(value), prec);
        }
        ulong bits = BitConverter.DoubleToUInt64Bits(value);
        bool neg = (bits >> 63) != 0;
        int biased = (int)((bits >> 52) & 0x7FF);
        ulong mant = bits & 0xFFFFFFFFFFFFFUL;
        long e;
        if (biased == 0)
        {
            e = -1074;
        }
        else
        {
            mant |= 1UL << 52;
            e = biased - 1075;
        }
        return Round(new BigInt(mant), e, sticky: false, neg, prec, rounding, out _);
    }

    /// <summary>Returns this value rounded to a new precision (mpfr_prec_round).</summary>
    public Real WithPrecision(int precision, RealRounding rounding = RealRounding.NearestEven)
    {
        ValidatePrecision(precision);
        return _kind switch
        {
            Kind.Zero => MakeZero(IsNeg, precision),
            Kind.Inf => MakeInf(IsNeg, precision),
            Kind.NaN => MakeNaN(precision),
            _ => Round(_man, _exp, sticky: false, IsNeg, precision, rounding, out _),
        };
    }

    /// <summary>Exact multiplication by 2^shift (mpfr_mul_2si).</summary>
    public Real ShiftLeft(long shift) => _kind == Kind.Normal ? new Real(_man, _exp + shift, _prec, _kind, IsNeg) : this;

    /// <summary>Exact division by 2^shift (mpfr_div_2si).</summary>
    public Real ShiftRight(long shift) => ShiftLeft(-shift);

    // ---------------------------------------------------------------
    // Conversion to double
    // ---------------------------------------------------------------

    /// <summary>Converts to double with the given rounding (mpfr_get_d).</summary>
    public double ToDouble(RealRounding rounding = RealRounding.NearestEven) => ToIeee(IeeeFormat.Binary64, rounding);

    /// <summary>Converts to float with the given rounding (correctly rounded, including subnormals).</summary>
    public float ToSingle(RealRounding rounding = RealRounding.NearestEven) => (float)ToIeee(IeeeFormat.Binary32, rounding);

    /// <summary>Converts to System.Half with the given rounding (correctly rounded, including subnormals).</summary>
    public System.Half ToHalf(RealRounding rounding = RealRounding.NearestEven) => (System.Half)ToIeee(IeeeFormat.Binary16, rounding);

    private double ToIeee(IeeeFormat format, RealRounding rounding) => _kind switch
    {
        Kind.NaN => double.NaN,
        Kind.Inf => IsNeg ? double.NegativeInfinity : double.PositiveInfinity,
        Kind.Zero => IsNeg ? -0.0 : 0.0,
        _ => RoundToIeee(_man, _exp, sticky: false, IsNeg, format, rounding),
    };

    /// <summary>Truncates to a BigInt (mpfr_get_z with RNDZ); throws on NaN/Inf.</summary>
    public BigInt ToBigInt(RealRounding rounding = RealRounding.TowardZero)
    {
        if (_kind is Kind.NaN or Kind.Inf)
        {
            throw new OverflowException("Cannot convert NaN or Infinity to BigInt.");
        }
        if (_kind == Kind.Zero)
        {
            return default;
        }
        Real r = RoundToInteger(rounding);
        BigInt mag = r._exp >= 0 ? r._man.ShiftLeft(r._exp) : r._man.ShiftRightTruncate(-r._exp);
        return r.IsNeg ? -mag : mag;
    }
}
