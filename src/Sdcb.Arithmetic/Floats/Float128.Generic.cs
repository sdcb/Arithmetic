using System.Buffers.Binary;
using System.Numerics;

namespace Sdcb.Arithmetic;

/// <summary>.NET generic math: <see cref="IBinaryFloatingPointIeee754{TSelf}"/> with <see cref="double"/>-like semantics.</summary>
public readonly partial struct Float128 :
    IBinaryFloatingPointIeee754<Float128>,
    IMinMaxValue<Float128>,
    ISpanFormattable,
    ISpanParsable<Float128>,
    IUtf8SpanFormattable,
    IUtf8SpanParsable<Float128>
{
    static Float128 IAdditiveIdentity<Float128, Float128>.AdditiveIdentity => Zero;
    static Float128 IMultiplicativeIdentity<Float128, Float128>.MultiplicativeIdentity => One;
    static int INumberBase<Float128>.Radix => 2;
    static Float128 IBinaryNumber<Float128>.AllBitsSet => new(UInt128.MaxValue);

    // ---------------------------------------------------------------
    // INumberBase predicates
    // ---------------------------------------------------------------

    static bool INumberBase<Float128>.IsCanonical(Float128 value) => true;
    static bool INumberBase<Float128>.IsComplexNumber(Float128 value) => false;
    static bool INumberBase<Float128>.IsImaginaryNumber(Float128 value) => false;

    /// <summary>True for every value except NaN.</summary>
    public static bool IsRealNumber(Float128 value) => !IsNaN(value);

    /// <summary>True for finite integral values (±0 included).</summary>
    public static bool IsInteger(Float128 value)
    {
        int exp = value.BiasedExponent;
        if (exp == MaxBiasedExponent)
        {
            return false;
        }
        if (exp >= ExponentBias + FractionBits)
        {
            return true;
        }
        if (exp < ExponentBias)
        {
            return IsZero(value);
        }
        int fracBits = ExponentBias + FractionBits - exp;
        return (value._bits & ((UInt128.One << fracBits) - 1)) == 0;
    }

    /// <summary>True for even integers (±0 included).</summary>
    public static bool IsEvenInteger(Float128 value) => IsInteger(value) && !UnitBit(value);

    /// <summary>True for odd integers.</summary>
    public static bool IsOddInteger(Float128 value) => IsInteger(value) && UnitBit(value);

    private static bool UnitBit(Float128 value)
    {
        int exp = value.BiasedExponent;
        if (exp < ExponentBias || exp > ExponentBias + FractionBits)
        {
            return false;
        }
        int fracBits = ExponentBias + FractionBits - exp;
        UInt128 sig = value.Fraction | Hidden;
        return ((sig >> fracBits) & 1) != 0;
    }

    /// <summary>True for positive powers of two (subnormal ones included).</summary>
    public static bool IsPow2(Float128 value)
    {
        if (value.SignBit)
        {
            return false;
        }
        int exp = value.BiasedExponent;
        UInt128 f = value.Fraction;
        if (exp == MaxBiasedExponent)
        {
            return false;
        }
        return exp == 0 ? f != 0 && (f & (f - 1)) == 0 : f == 0;
    }

    // ---------------------------------------------------------------
    // INumber / INumberBase (double semantics)
    // ---------------------------------------------------------------

    /// <summary>-1, 0 or +1; throws <see cref="ArithmeticException"/> for NaN (like <see cref="Math.Sign(double)"/>).</summary>
    public static int Sign(Float128 value)
    {
        if (IsNaN(value))
        {
            throw new ArithmeticException("Function does not accept floating point Not-a-Number values.");
        }
        return IsZero(value) ? 0 : value.SignBit ? -1 : 1;
    }

    /// <summary>|value| (clears the sign bit).</summary>
    public static Float128 Abs(Float128 value) => new(value._bits & ~SignMask);

    /// <summary>The magnitude of value with the sign bit of sign.</summary>
    public static Float128 CopySign(Float128 value, Float128 sign) =>
        new((value._bits & ~SignMask) | (sign._bits & SignMask));

    /// <summary>Clamps like <see cref="Math.Clamp(double, double, double)"/>; throws when min &gt; max.</summary>
    public static Float128 Clamp(Float128 value, Float128 min, Float128 max)
    {
        if (min > max)
        {
            throw new ArgumentException($"'{min}' cannot be greater than {max}.");
        }
        if (value < min)
        {
            return min;
        }
        if (value > max)
        {
            return max;
        }
        return value;
    }

    /// <summary>IEEE 754:2019 maximum: NaN if either is NaN; +0 &gt; -0.</summary>
    public static Float128 Max(Float128 x, Float128 y)
    {
        if (x != y)
        {
            return !IsNaN(x) ? (y < x ? x : y) : x;
        }
        return y.SignBit ? x : y;
    }

    /// <summary>IEEE 754:2019 minimum: NaN if either is NaN; -0 &lt; +0.</summary>
    public static Float128 Min(Float128 x, Float128 y)
    {
        if (x != y)
        {
            return !IsNaN(x) ? (x < y ? x : y) : x;
        }
        return x.SignBit ? x : y;
    }

    /// <summary>IEEE 754:2019 maximumNumber: a number wins over NaN.</summary>
    public static Float128 MaxNumber(Float128 x, Float128 y)
    {
        if (x != y)
        {
            return !IsNaN(y) ? (y < x ? x : y) : x;
        }
        return y.SignBit ? x : y;
    }

    /// <summary>IEEE 754:2019 minimumNumber: a number wins over NaN.</summary>
    public static Float128 MinNumber(Float128 x, Float128 y)
    {
        if (x != y)
        {
            return !IsNaN(y) ? (x < y ? x : y) : x;
        }
        return x.SignBit ? x : y;
    }

    public static Float128 MaxMagnitude(Float128 x, Float128 y)
    {
        Float128 ax = Abs(x), ay = Abs(y);
        if (ax > ay || IsNaN(ax))
        {
            return x;
        }
        if (ax == ay)
        {
            return x.SignBit ? y : x;
        }
        return y;
    }

    public static Float128 MaxMagnitudeNumber(Float128 x, Float128 y)
    {
        Float128 ax = Abs(x), ay = Abs(y);
        if (ax > ay || IsNaN(ay))
        {
            return x;
        }
        if (ax == ay)
        {
            return x.SignBit ? y : x;
        }
        return y;
    }

    public static Float128 MinMagnitude(Float128 x, Float128 y)
    {
        Float128 ax = Abs(x), ay = Abs(y);
        if (ax < ay || IsNaN(ax))
        {
            return x;
        }
        if (ax == ay)
        {
            return x.SignBit ? x : y;
        }
        return y;
    }

    public static Float128 MinMagnitudeNumber(Float128 x, Float128 y)
    {
        Float128 ax = Abs(x), ay = Abs(y);
        if (ax < ay || IsNaN(ay))
        {
            return x;
        }
        if (ax == ay)
        {
            return x.SignBit ? x : y;
        }
        return y;
    }

    // ---------------------------------------------------------------
    // IFloatingPoint rounding
    // ---------------------------------------------------------------

    /// <summary>Largest integer &lt;= x (exact).</summary>
    public static Float128 Floor(Float128 x) => RoundToInt(x, IntRounding.Floor);

    /// <summary>Smallest integer &gt;= x (exact).</summary>
    public static Float128 Ceiling(Float128 x) => RoundToInt(x, IntRounding.Ceiling);

    /// <summary>Integer part of x, rounding toward zero (exact).</summary>
    public static Float128 Truncate(Float128 x) => RoundToInt(x, IntRounding.Truncate);

    /// <summary>Nearest integer, ties to even (exact).</summary>
    public static Float128 Round(Float128 x) => RoundToInt(x, IntRounding.HalfEven);

    /// <summary>Nearest integer with the given midpoint rule (directed modes round toward that direction).</summary>
    public static Float128 Round(Float128 x, MidpointRounding mode) => mode switch
    {
        MidpointRounding.ToEven => RoundToInt(x, IntRounding.HalfEven),
        MidpointRounding.AwayFromZero => RoundToInt(x, IntRounding.HalfAway),
        MidpointRounding.ToZero => RoundToInt(x, IntRounding.Truncate),
        MidpointRounding.ToNegativeInfinity => RoundToInt(x, IntRounding.Floor),
        MidpointRounding.ToPositiveInfinity => RoundToInt(x, IntRounding.Ceiling),
        _ => throw new ArgumentException("Invalid MidpointRounding value.", nameof(mode)),
    };

    /// <summary>Rounds to the given number of decimal places, ties to even (exact decimal rounding, then correctly rounded).</summary>
    public static Float128 Round(Float128 x, int digits) => Round(x, digits, MidpointRounding.ToEven);

    /// <summary>
    /// Rounds to <paramref name="digits"/> decimal places: the exact decimal rounding of x under
    /// <paramref name="mode"/>, then correctly rounded (nearest-even) back to binary128.
    /// </summary>
    public static Float128 Round(Float128 x, int digits, MidpointRounding mode)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(digits);
        if (digits == 0)
        {
            return Round(x, mode);
        }
        if (!IsFinite(x) || IsZero(x) || x.BiasedExponent >= ExponentBias + FractionBits)
        {
            return IsNaN(x) ? NaN : x;
        }
        // Far below 10^-digits only the sign matters: a small stand-in keeps the rational tiny.
        long floor = -(4L * digits + 16);
        BigRational q = ILogB(x) < floor
            ? new BigRational(x.SignBit ? BigInt.MinusOne : BigInt.One, BigInt.One.ShiftLeft(-floor + 1))
            : (BigRational)x;
        BigRational r = BigRational.Round(q, digits, mode);
        return r.IsZero ? ZeroOf(x.SignBit) : FromRational(r);
    }

    // ---------------------------------------------------------------
    // IFloatingPoint layout (double semantics: Exponent = biased - bias, Significand includes the hidden bit)
    // ---------------------------------------------------------------

    private short IeeeExponent => (short)(BiasedExponent - ExponentBias);

    private UInt128 IeeeSignificand => BiasedExponent != 0 ? Fraction | Hidden : Fraction;

    int IFloatingPoint<Float128>.GetExponentByteCount() => sizeof(short);

    int IFloatingPoint<Float128>.GetExponentShortestBitLength()
    {
        short e = IeeeExponent;
        return e >= 0 ? 16 - BitOperations.LeadingZeroCount((uint)e) + 16 : 17 - BitOperations.LeadingZeroCount((uint)(ushort)~e) + 16;
    }

    int IFloatingPoint<Float128>.GetSignificandBitLength() => SignificandBits;

    int IFloatingPoint<Float128>.GetSignificandByteCount() => 16;

    bool IFloatingPoint<Float128>.TryWriteExponentBigEndian(Span<byte> destination, out int bytesWritten)
    {
        bool ok = BinaryPrimitives.TryWriteInt16BigEndian(destination, IeeeExponent);
        bytesWritten = ok ? sizeof(short) : 0;
        return ok;
    }

    bool IFloatingPoint<Float128>.TryWriteExponentLittleEndian(Span<byte> destination, out int bytesWritten)
    {
        bool ok = BinaryPrimitives.TryWriteInt16LittleEndian(destination, IeeeExponent);
        bytesWritten = ok ? sizeof(short) : 0;
        return ok;
    }

    bool IFloatingPoint<Float128>.TryWriteSignificandBigEndian(Span<byte> destination, out int bytesWritten)
    {
        bool ok = BinaryPrimitives.TryWriteUInt128BigEndian(destination, IeeeSignificand);
        bytesWritten = ok ? 16 : 0;
        return ok;
    }

    bool IFloatingPoint<Float128>.TryWriteSignificandLittleEndian(Span<byte> destination, out int bytesWritten)
    {
        bool ok = BinaryPrimitives.TryWriteUInt128LittleEndian(destination, IeeeSignificand);
        bytesWritten = ok ? 16 : 0;
        return ok;
    }

    // ---------------------------------------------------------------
    // Bitwise operators on the bit pattern (interface-only, like double)
    // ---------------------------------------------------------------

    static Float128 IBitwiseOperators<Float128, Float128, Float128>.operator &(Float128 left, Float128 right) => new(left._bits & right._bits);
    static Float128 IBitwiseOperators<Float128, Float128, Float128>.operator |(Float128 left, Float128 right) => new(left._bits | right._bits);
    static Float128 IBitwiseOperators<Float128, Float128, Float128>.operator ^(Float128 left, Float128 right) => new(left._bits ^ right._bits);
    static Float128 IBitwiseOperators<Float128, Float128, Float128>.operator ~(Float128 value) => new(~value._bits);
}
