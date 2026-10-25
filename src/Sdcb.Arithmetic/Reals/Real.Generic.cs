using System.Buffers.Binary;
using System.Globalization;
using System.Numerics;
using Sdcb.Arithmetic.Core;

namespace Sdcb.Arithmetic;

/// <summary>
/// .NET generic math: <see cref="IBinaryFloatingPointIeee754{TSelf}"/>. In generic code unary
/// functions round to the argument's precision, binary operations to the larger operand precision,
/// and constants (Zero, One, Pi, E, Tau, Epsilon, ...) use <see cref="DefaultPrecision"/>; all
/// results are correctly rounded to nearest-even. The bitwise operators required by
/// <see cref="IBinaryNumber{TSelf}"/> are available only through the interface.
/// </summary>
public readonly partial struct Real :
    IBinaryFloatingPointIeee754<Real>,
    IUtf8SpanFormattable,
    IUtf8SpanParsable<Real>
{
    private const RealRounding Rne = RealRounding.NearestEven;

    static Real IAdditiveIdentity<Real, Real>.AdditiveIdentity => Zero;
    static Real IMultiplicativeIdentity<Real, Real>.MultiplicativeIdentity => One;
    static int INumberBase<Real>.Radix => 2;

    /// <summary>-1 (the <see cref="ISignedNumber{TSelf}"/> name for <see cref="MinusOne"/>).</summary>
    public static Real NegativeOne => MinusOne;

    /// <summary>
    /// The smallest positive value at <see cref="DefaultPrecision"/> (BitIncrement(Zero)). Exponents
    /// are otherwise unbounded, so this is a conventional floor, not a hard underflow threshold.
    /// </summary>
    public static Real Epsilon
    {
        get
        {
            int p = DefaultPrecision;
            return new Real(BigInt.One << (p - 1), EminUlp, p, Kind.Normal, negative: false);
        }
    }

    static Real IFloatingPointConstants<Real>.E => E();
    static Real IFloatingPointConstants<Real>.Pi => Pi();
    static Real IFloatingPointConstants<Real>.Tau => Tau();
    static Real IBinaryNumber<Real>.AllBitsSet => NaN;

    /// <summary>2π, correctly rounded.</summary>
    public static Real Tau(int precision = 0, RealRounding rounding = Rne) => Pi(precision, rounding).ShiftLeft(1);

    // ---------------------------------------------------------------
    // INumberBase predicates
    // ---------------------------------------------------------------

    static bool INumberBase<Real>.IsCanonical(Real value) => true;
    static bool INumberBase<Real>.IsComplexNumber(Real value) => false;
    static bool INumberBase<Real>.IsFinite(Real value) => value.IsFinite;
    static bool INumberBase<Real>.IsImaginaryNumber(Real value) => false;
    static bool INumberBase<Real>.IsInfinity(Real value) => value.IsInfinity;
    static bool INumberBase<Real>.IsInteger(Real value) => value.IsInteger;
    static bool INumberBase<Real>.IsNaN(Real value) => value.IsNaN;
    static bool INumberBase<Real>.IsNegativeInfinity(Real value) => value.IsNegativeInfinity;
    static bool INumberBase<Real>.IsPositiveInfinity(Real value) => value.IsPositiveInfinity;
    static bool INumberBase<Real>.IsRealNumber(Real value) => !value.IsNaN;
    static bool INumberBase<Real>.IsSubnormal(Real value) => false;
    static bool INumberBase<Real>.IsZero(Real value) => value.IsZero;

    /// <summary>True when the sign bit is set (including -0, -∞).</summary>
    public static bool IsNegative(Real value) => value.IsNeg;

    /// <summary>True when the sign bit is clear (including +0, +∞).</summary>
    public static bool IsPositive(Real value) => !value.IsNeg;

    /// <summary>True for finite nonzero values (Real has no subnormals).</summary>
    public static bool IsNormal(Real value) => value._kind == Kind.Normal;

    public static bool IsEvenInteger(Real value) =>
        value._kind == Kind.Zero || (value.IsInteger && value._kind == Kind.Normal && !UnitBit(value));

    public static bool IsOddInteger(Real value) =>
        value._kind == Kind.Normal && value.IsInteger && UnitBit(value);

    private static bool UnitBit(Real v) => v._exp == 0 ? v._man.IsOdd : v._exp < 0 && v._man.TestBit(-v._exp);

    /// <summary>True for positive powers of two.</summary>
    public static bool IsPow2(Real value) => value._kind == Kind.Normal && !value.IsNeg && value._man.IsPowerOfTwo;

    // ---------------------------------------------------------------
    // INumber
    // ---------------------------------------------------------------

    static int INumber<Real>.Sign(Real value) => value.Sign;

    public static Real Abs(Real value) => value.Abs();

    /// <summary>The magnitude of value with the sign bit of sign.</summary>
    public static Real CopySign(Real value, Real sign) => value.CopySign(sign);

    public static Real Clamp(Real value, Real min, Real max)
    {
        if (min > max || min.IsNaN || max.IsNaN)
        {
            throw new ArgumentException("min must be less than or equal to max.", nameof(min));
        }
        if (value.IsNaN)
        {
            return value;
        }
        return value < min ? min : value > max ? max : value;
    }

    /// <summary>IEEE 754:2019 maximum: NaN if either argument is NaN; +0 &gt; -0.</summary>
    static Real INumber<Real>.Max(Real x, Real y)
    {
        if (x._kind == Kind.NaN)
        {
            return x;
        }
        if (y._kind == Kind.NaN)
        {
            return y;
        }
        return Max(x, y);
    }

    /// <summary>IEEE 754:2019 minimum: NaN if either argument is NaN; -0 &lt; +0.</summary>
    static Real INumber<Real>.Min(Real x, Real y)
    {
        if (x._kind == Kind.NaN)
        {
            return x;
        }
        if (y._kind == Kind.NaN)
        {
            return y;
        }
        return Min(x, y);
    }

    /// <summary>IEEE 754:2019 maximumNumber (same as <see cref="Max(Real, Real)"/>): a number wins over NaN.</summary>
    public static Real MaxNumber(Real x, Real y) => Max(x, y);

    /// <summary>IEEE 754:2019 minimumNumber (same as <see cref="Min(Real, Real)"/>): a number wins over NaN.</summary>
    public static Real MinNumber(Real x, Real y) => Min(x, y);

    private static int CompareMagnitude(Real x, Real y) => CompareValues(x.Abs(), y.Abs());

    public static Real MaxMagnitude(Real x, Real y)
    {
        if (x.IsNaN)
        {
            return x;
        }
        if (y.IsNaN)
        {
            return y;
        }
        int c = CompareMagnitude(x, y);
        return c > 0 ? x : c < 0 ? y : x.IsNeg ? y : x;
    }

    public static Real MinMagnitude(Real x, Real y)
    {
        if (x.IsNaN)
        {
            return x;
        }
        if (y.IsNaN)
        {
            return y;
        }
        int c = CompareMagnitude(x, y);
        return c < 0 ? x : c > 0 ? y : x.IsNeg ? x : y;
    }

    public static Real MaxMagnitudeNumber(Real x, Real y)
    {
        if (y.IsNaN)
        {
            return x;
        }
        if (x.IsNaN)
        {
            return y;
        }
        int c = CompareMagnitude(x, y);
        return c > 0 ? x : c < 0 ? y : x.IsNeg ? y : x;
    }

    public static Real MinMagnitudeNumber(Real x, Real y)
    {
        if (y.IsNaN)
        {
            return x;
        }
        if (x.IsNaN)
        {
            return y;
        }
        int c = CompareMagnitude(x, y);
        return c < 0 ? x : c > 0 ? y : x.IsNeg ? x : y;
    }

    // ---------------------------------------------------------------
    // Remaining operators
    // ---------------------------------------------------------------

    /// <summary>Truncated remainder with the sign of the dividend (like C# double %, mpfr_fmod).</summary>
    public static Real operator %(Real a, Real b) => TruncatedRemainder(a, b);

    /// <summary>value + 1 rounded to value's precision.</summary>
    public static Real operator ++(Real value) => Add(value, One, value.Precision, Rne);

    /// <summary>value - 1 rounded to value's precision.</summary>
    public static Real operator --(Real value) => Subtract(value, One, value.Precision, Rne);

    // ---------------------------------------------------------------
    // Rounding (IFloatingPoint)
    // ---------------------------------------------------------------

    public static Real Floor(Real x) => x.Floor();
    public static Real Ceiling(Real x) => x.Ceiling();
    public static Real Truncate(Real x) => x.Truncate();

    /// <summary>Nearest integer, ties to even.</summary>
    public static Real Round(Real x) => x.RoundHalfEven();

    public static Real Round(Real x, MidpointRounding mode) => mode switch
    {
        MidpointRounding.ToEven => x.RoundHalfEven(),
        MidpointRounding.AwayFromZero => x.RoundHalfAway(),
        MidpointRounding.ToZero => x.Truncate(),
        MidpointRounding.ToNegativeInfinity => x.Floor(),
        MidpointRounding.ToPositiveInfinity => x.Ceiling(),
        _ => throw new ArgumentException("Invalid MidpointRounding value.", nameof(mode)),
    };

    public static Real Round(Real x, int digits) => Round(x, digits, MidpointRounding.ToEven);

    /// <summary>
    /// Rounds to <paramref name="digits"/> decimal places: the exact decimal rounding of x, then
    /// correctly rounded back to x's precision.
    /// </summary>
    public static Real Round(Real x, int digits, MidpointRounding mode)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(digits);
        if (digits == 0)
        {
            return Round(x, mode);
        }
        if (x._kind != Kind.Normal || x._exp >= 0)
        {
            return x;
        }
        int prec = x.Precision;
        long floor = -(4L * digits + 16);
        // Far below 10^-digits only the sign matters: substitute a small stand-in to keep the rational tiny.
        BigRational q = x.BinaryExponent < floor
            ? new BigRational(x.IsNeg ? BigInt.MinusOne : BigInt.One, BigInt.One.ShiftLeft(-floor + 1))
            : x.ToRational();
        BigRational r = BigRational.Round(q, digits, mode);
        return r.IsZero ? MakeZero(x.IsNeg, prec) : FromRational(r, prec, Rne);
    }

    // ---------------------------------------------------------------
    // IFloatingPointIeee754
    // ---------------------------------------------------------------

    public static Real BitIncrement(Real x) => x.NextAbove();
    public static Real BitDecrement(Real x) => x.NextBelow();
    public static Real FusedMultiplyAdd(Real left, Real right, Real addend) => FusedMultiplyAdd(left, right, addend, 0, Rne);
    public static Real Ieee754Remainder(Real left, Real right) => Ieee754Remainder(left, right, 0, Rne);
    public static Real ScaleB(Real x, int n) => x.ShiftLeft(n);

    /// <summary>floor(log2|x|) clamped to int; int.MinValue for zero, int.MaxValue for NaN and infinity.</summary>
    public static int ILogB(Real x) => x._kind switch
    {
        Kind.Zero => int.MinValue,
        Kind.NaN or Kind.Inf => int.MaxValue,
        _ => (int)Math.Clamp(x.BinaryExponent - 1, int.MinValue + 1, int.MaxValue - 1),
    };

    // ---------------------------------------------------------------
    // IFloatingPoint layout members (exponent E with value = significand * 2^(E - (Precision - 1)))
    // ---------------------------------------------------------------

    private long IeeeExponent => _kind == Kind.Normal ? BinaryExponent - 1 : 0;

    int IFloatingPoint<Real>.GetExponentByteCount() => sizeof(long);

    int IFloatingPoint<Real>.GetExponentShortestBitLength()
    {
        long e = IeeeExponent;
        return e >= 0 ? 64 - BitOperations.LeadingZeroCount((ulong)e) : 65 - BitOperations.LeadingZeroCount((ulong)~e);
    }

    int IFloatingPoint<Real>.GetSignificandBitLength() => Precision;

    int IFloatingPoint<Real>.GetSignificandByteCount() => (Precision + 7) >> 3;

    bool IFloatingPoint<Real>.TryWriteExponentBigEndian(Span<byte> destination, out int bytesWritten) =>
        TryWriteLong(destination, IeeeExponent, bigEndian: true, out bytesWritten);

    bool IFloatingPoint<Real>.TryWriteExponentLittleEndian(Span<byte> destination, out int bytesWritten) =>
        TryWriteLong(destination, IeeeExponent, bigEndian: false, out bytesWritten);

    bool IFloatingPoint<Real>.TryWriteSignificandBigEndian(Span<byte> destination, out int bytesWritten) =>
        TryWriteSignificand(destination, bigEndian: true, out bytesWritten);

    bool IFloatingPoint<Real>.TryWriteSignificandLittleEndian(Span<byte> destination, out int bytesWritten) =>
        TryWriteSignificand(destination, bigEndian: false, out bytesWritten);

    private static bool TryWriteLong(Span<byte> destination, long value, bool bigEndian, out int bytesWritten)
    {
        bool ok = bigEndian
            ? BinaryPrimitives.TryWriteInt64BigEndian(destination, value)
            : BinaryPrimitives.TryWriteInt64LittleEndian(destination, value);
        bytesWritten = ok ? sizeof(long) : 0;
        return ok;
    }

    private bool TryWriteSignificand(Span<byte> destination, bool bigEndian, out int bytesWritten)
    {
        int count = (Precision + 7) >> 3;
        if (destination.Length < count)
        {
            bytesWritten = 0;
            return false;
        }
        Span<byte> d = destination.Slice(0, count);
        d.Clear();
        if (_kind == Kind.Normal)
        {
            byte[] bytes = _man.ToMagnitudeBytes(bigEndian: false);
            bytes.AsSpan(0, Math.Min(bytes.Length, count)).CopyTo(d);
            if (bigEndian)
            {
                d.Reverse();
            }
        }
        bytesWritten = count;
        return true;
    }

    // ---------------------------------------------------------------
    // Bitwise operators (IBinaryNumber; interface-only)
    // ---------------------------------------------------------------

    static Real IBitwiseOperators<Real, Real, Real>.operator &(Real left, Real right) => Bitwise(left, right, 0);
    static Real IBitwiseOperators<Real, Real, Real>.operator |(Real left, Real right) => Bitwise(left, right, 1);
    static Real IBitwiseOperators<Real, Real, Real>.operator ^(Real left, Real right) => Bitwise(left, right, 2);

    /// <summary>Complements the significand bits within the value's precision window and flips the sign; NaN for NaN/∞.</summary>
    static Real IBitwiseOperators<Real, Real, Real>.operator ~(Real value)
    {
        switch (value._kind)
        {
            case Kind.NaN:
            case Kind.Inf:
                return MakeNaN(value.Precision);
            case Kind.Zero:
                return MakeZero(!value.IsNeg, value.Precision);
        }
        int p = value.Precision;
        BigInt m = ((BigInt.One << p) - BigInt.One) - value._man;
        return m.IsZero ? MakeZero(!value.IsNeg, p) : Round(m, value._exp, sticky: false, !value.IsNeg, p, Rne, out _);
    }

    /// <summary>
    /// Applies AND/OR/XOR to the exact binary expansions of |left| and |right| (aligned at the
    /// lowest significand bit) and to the sign bits, rounding to the larger precision.
    /// </summary>
    private static Real Bitwise(Real a, Real b, int op)
    {
        int prec = ResolvePrec(a, b, 0);
        if (a._kind is Kind.NaN or Kind.Inf || b._kind is Kind.NaN or Kind.Inf)
        {
            return MakeNaN(prec);
        }
        bool neg = op switch
        {
            0 => a.IsNeg & b.IsNeg,
            1 => a.IsNeg | b.IsNeg,
            _ => a.IsNeg ^ b.IsNeg,
        };
        bool az = a._kind == Kind.Zero, bz = b._kind == Kind.Zero;
        if (az || bz)
        {
            if (op == 0 || (az && bz))
            {
                return MakeZero(neg, prec);
            }
            Real nz = az ? b : a;
            return Round(nz._man, nz._exp, sticky: false, neg, prec, Rne, out _);
        }
        // Disjoint bit ranges: AND is empty, OR/XOR are the exact sum.
        bool disjoint = a._exp >= b._exp + b.Precision || b._exp >= a._exp + a.Precision;
        if (disjoint)
        {
            if (op == 0)
            {
                return MakeZero(neg, prec);
            }
            Real s = Add(a.Abs(), b.Abs(), prec, Rne);
            return neg ? -s : s;
        }
        long e = Math.Min(a._exp, b._exp);
        BigInt ma = a._man.ShiftLeft(a._exp - e), mb = b._man.ShiftLeft(b._exp - e);
        BigInt m = op switch
        {
            0 => ma & mb,
            1 => ma | mb,
            _ => ma ^ mb,
        };
        return m.IsZero ? MakeZero(neg, prec) : Round(m, e, sticky: false, neg, prec, Rne, out _);
    }

    // ---------------------------------------------------------------
    // Single-argument function overloads (generic-math signatures)
    // ---------------------------------------------------------------

    public static Real Sqrt(Real x) => Sqrt(x, 0, Rne);
    public static Real Cbrt(Real x) => Cbrt(x, 0, Rne);
    public static Real Hypot(Real x, Real y) => Hypot(x, y, 0, Rne);
    public static Real RootN(Real x, int n) => RootN(x, (long)n, 0, Rne);
    public static Real Pow(Real x, Real y) => Pow(x, y, 0, Rne);

    public static Real Exp(Real x) => Exp(x, 0, Rne);
    public static Real ExpM1(Real x) => ExpM1(x, 0, Rne);
    public static Real Exp2(Real x) => Exp2(x, 0, Rne);
    public static Real Exp2M1(Real x) => Exp2M1(x, 0, Rne);
    public static Real Exp10(Real x) => Exp10(x, 0, Rne);
    public static Real Exp10M1(Real x) => Exp10M1(x, 0, Rne);

    public static Real Log(Real x) => Log(x, 0, Rne);
    public static Real LogP1(Real x) => LogP1(x, 0, Rne);
    public static Real Log2(Real x) => Log2(x, 0, Rne);
    public static Real Log2P1(Real x) => Log2P1(x, 0, Rne);
    public static Real Log10(Real x) => Log10(x, 0, Rne);
    public static Real Log10P1(Real x) => Log10P1(x, 0, Rne);

    public static Real Sinh(Real x) => Sinh(x, 0, Rne);
    public static Real Cosh(Real x) => Cosh(x, 0, Rne);
    public static Real Tanh(Real x) => Tanh(x, 0, Rne);
    public static Real Asinh(Real x) => Asinh(x, 0, Rne);
    public static Real Acosh(Real x) => Acosh(x, 0, Rne);
    public static Real Atanh(Real x) => Atanh(x, 0, Rne);

    public static Real Sin(Real x) => Sin(x, 0, Rne);
    public static Real Cos(Real x) => Cos(x, 0, Rne);
    public static Real Tan(Real x) => Tan(x, 0, Rne);
    public static Real Asin(Real x) => Asin(x, 0, Rne);
    public static Real Acos(Real x) => Acos(x, 0, Rne);
    public static Real Atan(Real x) => Atan(x, 0, Rne);
    public static Real Atan2(Real y, Real x) => Atan2(y, x, 0, Rne);
    public static Real SinPi(Real x) => SinPi(x, 0, Rne);
    public static Real CosPi(Real x) => CosPi(x, 0, Rne);
    public static Real TanPi(Real x) => TanPi(x, 0, Rne);
    public static (Real Sin, Real Cos) SinCos(Real x) => SinCos(x, 0, Rne);
    public static Real AsinPi(Real x) => AsinPi(x, 0, Rne);
    public static Real AcosPi(Real x) => AcosPi(x, 0, Rne);
    public static Real AtanPi(Real x) => AtanPi(x, 0, Rne);
    public static Real Atan2Pi(Real y, Real x) => Atan2Pi(y, x, 0, Rne);

    /// <summary>(sin(πx), cos(πx)), each correctly rounded.</summary>
    public static (Real SinPi, Real CosPi) SinCosPi(Real x) => (SinPi(x), CosPi(x));

    // ---------------------------------------------------------------
    // New correctly rounded functions: asinpi/acospi/atanpi/atan2pi, degree conversion, log base b
    // ---------------------------------------------------------------

    private static Real ExactDyadic(long numerator, int log2Den, bool negative, int prec, RealRounding rnd)
    {
        Real r = Round(new BigInt(numerator), -log2Den, sticky: false, negative, prec, rnd, out _);
        return r;
    }

    /// <summary>asin(x)/π (mpfr_asinpi).</summary>
    public static Real AsinPi(Real x, int precision = 0, RealRounding rounding = Rne)
    {
        int prec = precision != 0 ? precision : x.Precision;
        if (x._kind == Kind.NaN || x._kind == Kind.Inf)
        {
            return MakeNaN(prec);
        }
        if (x._kind == Kind.Zero)
        {
            return MakeZero(x.IsNeg, prec);
        }
        int c = CompareMagnitude(x, One);
        if (c > 0)
        {
            return MakeNaN(prec);
        }
        if (c == 0)
        {
            return ExactDyadic(1, 1, x.IsNeg, prec, rounding);
        }
        return Ziv(prec, rounding, wp => Approx(Divide(Asin(x, wp + 8), Pi(wp + 8), wp + 8), 4));
    }

    /// <summary>acos(x)/π (mpfr_acospi).</summary>
    public static Real AcosPi(Real x, int precision = 0, RealRounding rounding = Rne)
    {
        int prec = precision != 0 ? precision : x.Precision;
        if (x._kind == Kind.NaN || x._kind == Kind.Inf)
        {
            return MakeNaN(prec);
        }
        if (x._kind == Kind.Zero)
        {
            return ExactDyadic(1, 1, negative: false, prec, rounding);
        }
        int c = CompareMagnitude(x, One);
        if (c > 0)
        {
            return MakeNaN(prec);
        }
        if (c == 0)
        {
            return x.IsNeg ? ExactDyadic(1, 0, negative: false, prec, rounding) : MakeZero(negative: false, prec);
        }
        return Ziv(prec, rounding, wp => Approx(Divide(Acos(x, wp + 8), Pi(wp + 8), wp + 8), 4));
    }

    /// <summary>atan(x)/π (mpfr_atanpi).</summary>
    public static Real AtanPi(Real x, int precision = 0, RealRounding rounding = Rne)
    {
        int prec = precision != 0 ? precision : x.Precision;
        switch (x._kind)
        {
            case Kind.NaN:
                return MakeNaN(prec);
            case Kind.Zero:
                return MakeZero(x.IsNeg, prec);
            case Kind.Inf:
                return ExactDyadic(1, 1, x.IsNeg, prec, rounding);
        }
        if (CompareMagnitude(x, One) == 0)
        {
            return ExactDyadic(1, 2, x.IsNeg, prec, rounding);
        }
        return Ziv(prec, rounding, wp => Approx(Divide(Atan(x, wp + 8), Pi(wp + 8), wp + 8), 4));
    }

    /// <summary>atan2(y, x)/π (mpfr_atan2pi).</summary>
    public static Real Atan2Pi(Real y, Real x, int precision = 0, RealRounding rounding = Rne)
    {
        int prec = ResolvePrec(y, x, precision);
        if (y._kind == Kind.NaN || x._kind == Kind.NaN)
        {
            return MakeNaN(prec);
        }
        bool yn = y.IsNeg;
        // Exact results: multiples of 1/4 (y or x zero/infinite, or |y| == |x|).
        if (y._kind == Kind.Zero)
        {
            return x.IsNeg ? ExactDyadic(1, 0, yn, prec, rounding) : MakeZero(yn, prec);
        }
        if (y._kind == Kind.Inf)
        {
            return x._kind == Kind.Inf
                ? ExactDyadic(x.IsNeg ? 3 : 1, 2, yn, prec, rounding)
                : ExactDyadic(1, 1, yn, prec, rounding);
        }
        if (x._kind == Kind.Zero)
        {
            return ExactDyadic(1, 1, yn, prec, rounding);
        }
        if (x._kind == Kind.Inf)
        {
            return x.IsNeg ? ExactDyadic(1, 0, yn, prec, rounding) : MakeZero(yn, prec);
        }
        if (CompareMagnitude(y, x) == 0)
        {
            return ExactDyadic(x.IsNeg ? 3 : 1, 2, yn, prec, rounding);
        }
        return Ziv(prec, rounding, wp => Approx(Divide(Atan2(y, x, wp + 8), Pi(wp + 8), wp + 8), 4));
    }

    /// <summary>x·π/180, correctly rounded to x's precision.</summary>
    public static Real DegreesToRadians(Real degrees)
    {
        int prec = degrees.Precision;
        if (degrees._kind != Kind.Normal)
        {
            return degrees;
        }
        Real d180 = FromBigInt(180);
        return Ziv(prec, Rne, wp => Approx(Divide(Multiply(degrees, Pi(wp + 8), wp + 8), d180, wp + 8), 4));
    }

    /// <summary>x·180/π, correctly rounded to x's precision.</summary>
    public static Real RadiansToDegrees(Real radians)
    {
        int prec = radians.Precision;
        if (radians._kind != Kind.Normal)
        {
            return radians;
        }
        Real d180 = FromBigInt(180);
        return Ziv(prec, Rne, wp => Approx(Divide(Multiply(radians, d180, wp + 8), Pi(wp + 8), wp + 8), 4));
    }

    /// <summary>Logarithm of x in base newBase, correctly rounded to the larger precision.</summary>
    public static Real Log(Real x, Real newBase) => Log(x, newBase, 0, Rne);

    /// <summary>Logarithm of x in base newBase (log(x)/log(newBase), correctly rounded).</summary>
    public static Real Log(Real x, Real newBase, int precision, RealRounding rounding = Rne)
    {
        int prec = ResolvePrec(x, newBase, precision);
        Real b = newBase;
        if (x._kind != Kind.Normal || b._kind != Kind.Normal || x.IsNeg || b.IsNeg)
        {
            // Special values follow log(x)/log(b) with IEEE semantics.
            if (x._kind == Kind.NaN || b._kind == Kind.NaN || (x.IsNeg && x._kind != Kind.Zero) || (b.IsNeg && b._kind != Kind.Zero))
            {
                return MakeNaN(prec);
            }
            return Divide(Log(x, prec + 8), Log(b, prec + 8), prec, rounding);
        }
        int cb = CompareValues(b, One);
        if (cb == 0)
        {
            return MakeNaN(prec);
        }
        if (CompareValues(x, One) == 0)
        {
            return MakeZero(negative: false, prec);
        }
        if (b._man.IsPowerOfTwo && x._man.IsPowerOfTwo)
        {
            // log_{2^j}(2^i) = i / j exactly.
            return Divide(FromBigInt(x.BinaryExponent - 1), FromBigInt(b.BinaryExponent - 1), prec, rounding);
        }
        if (TryExactIntegerLog(x, b, prec, rounding, out Real exact))
        {
            return exact;
        }
        return Ziv(prec, rounding, wp => Approx(Divide(Log(x, wp + 8), Log(b, wp + 8), wp + 8), 4));
    }

    /// <summary>Detects x == b^k or b == x^k for small positive k, where Ziv's loop could not terminate.</summary>
    private static bool TryExactIntegerLog(Real x, Real b, int prec, RealRounding rounding, out Real result)
    {
        result = default;
        double lx = Math.Log2(Math.Abs(x.ToDouble(RealRounding.TowardZero))) , lb = Math.Log2(Math.Abs(b.ToDouble(RealRounding.TowardZero)));
        if (!double.IsFinite(lx) || !double.IsFinite(lb) || lx == 0 || lb == 0)
        {
            return false;
        }
        double ratio = lx / lb;
        if (ratio >= 1)
        {
            double k = Math.Round(ratio);
            if (k <= 64 && Math.Abs(ratio - k) < 1e-6 && IsExactPower(b, (long)k, x))
            {
                result = FromBigInt((long)k, prec, rounding);
                return true;
            }
        }
        else if (ratio > 0)
        {
            double k = Math.Round(1 / ratio);
            if (k <= 64 && Math.Abs(1 / ratio - k) < 1e-6 && IsExactPower(x, (long)k, b))
            {
                result = Divide(One, FromBigInt((long)k), prec, rounding);
                return true;
            }
        }
        return false;
    }

    private static bool IsExactPower(Real b, long k, Real target)
    {
        long p = (long)b.Precision * k + 64;
        if (p > MaxPrecision)
        {
            return false;
        }
        Real pw = Pow(b, k, (int)p, RealRounding.TowardZero);
        return pw._kind == Kind.Normal && CompareValues(pw, target) == 0;
    }

    // ---------------------------------------------------------------
    // NumberStyles parsing / UTF-8
    // ---------------------------------------------------------------

    /// <summary>Parses with .NET <see cref="NumberStyles"/> (hex/binary specifiers are not supported), rounding to DefaultPrecision.</summary>
    public static Real Parse(ReadOnlySpan<char> s, NumberStyles style, IFormatProvider? provider) =>
        TryParse(s, style, provider, out Real r) ? r : throw new FormatException("The value could not be parsed as a Real.");

    public static Real Parse(string s, NumberStyles style, IFormatProvider? provider)
    {
        ArgumentNullException.ThrowIfNull(s);
        return Parse(s.AsSpan(), style, provider);
    }

    public static bool TryParse(string? s, NumberStyles style, IFormatProvider? provider, out Real result) =>
        TryParse(s.AsSpan(), style, provider, out result);

    public static bool TryParse(ReadOnlySpan<char> s, NumberStyles style, IFormatProvider? provider, out Real result)
    {
        NumberText.ValidateStyle(style);
        if ((style & (NumberStyles.AllowHexSpecifier | NumberStyles.AllowBinarySpecifier)) != 0)
        {
            throw new ArgumentException("Hex and binary specifiers are not supported for Real.", nameof(style));
        }
        result = default;
        ReadOnlySpan<char> t = NumberText.TrimForSymbols(s, style);
        NumberFormatInfo nfi = NumberFormatInfo.GetInstance(provider);
        int prec = DefaultPrecision;
        if (t.Equals(nfi.NaNSymbol, StringComparison.OrdinalIgnoreCase))
        {
            result = MakeNaN(prec);
            return true;
        }
        if (t.Equals(nfi.PositiveInfinitySymbol, StringComparison.OrdinalIgnoreCase))
        {
            result = MakeInf(negative: false, prec);
            return true;
        }
        if (t.Equals(nfi.NegativeInfinitySymbol, StringComparison.OrdinalIgnoreCase))
        {
            result = MakeInf(negative: true, prec);
            return true;
        }
        if (NumberText.TryNormalize(s, style, provider, out string text, out _))
        {
            return TryParse(text, 0, Rne, out result);
        }
        // Also accept the invariant nan / inf / infinity spellings with an optional sign.
        ReadOnlySpan<char> word = t.Length > 0 && t[0] is '+' or '-' ? t.Slice(1) : t;
        if (word.Equals("nan", StringComparison.OrdinalIgnoreCase) || word.Equals("inf", StringComparison.OrdinalIgnoreCase) ||
            word.Equals("infinity", StringComparison.OrdinalIgnoreCase))
        {
            return TryParse(t, 0, Rne, out result);
        }
        return false;
    }

    public bool TryFormat(Span<byte> utf8Destination, out int bytesWritten, ReadOnlySpan<char> format, IFormatProvider? provider) =>
        NumberText.TryWriteUtf8(ToString(), utf8Destination, out bytesWritten);

    public static Real Parse(ReadOnlySpan<byte> utf8Text, IFormatProvider? provider) =>
        TryParse(utf8Text, provider, out Real r) ? r : throw new FormatException("The value could not be parsed as a Real.");

    public static bool TryParse(ReadOnlySpan<byte> utf8Text, IFormatProvider? provider, out Real result) =>
        NumberText.TryParseUtf8(utf8Text, (ReadOnlySpan<char> s, out Real r) => TryParse(s, provider, out r), out result);

    public static Real Parse(ReadOnlySpan<byte> utf8Text, NumberStyles style, IFormatProvider? provider) =>
        TryParse(utf8Text, style, provider, out Real r) ? r : throw new FormatException("The value could not be parsed as a Real.");

    public static bool TryParse(ReadOnlySpan<byte> utf8Text, NumberStyles style, IFormatProvider? provider, out Real result) =>
        NumberText.TryParseUtf8(utf8Text, (ReadOnlySpan<char> s, out Real r) => TryParse(s, style, provider, out r), out result);
}
