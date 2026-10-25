using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.Numerics;
using Sdcb.Arithmetic.Core;

namespace Sdcb.Arithmetic;

/// <summary>.NET generic math (<see cref="INumber{TSelf}"/>, <see cref="ISignedNumber{TSelf}"/>) and conversions.</summary>
public readonly partial struct BigRational :
    INumber<BigRational>,
    ISignedNumber<BigRational>,
    IUtf8SpanFormattable,
    IUtf8SpanParsable<BigRational>
{
    static BigRational IAdditiveIdentity<BigRational, BigRational>.AdditiveIdentity => default;
    static BigRational IMultiplicativeIdentity<BigRational, BigRational>.MultiplicativeIdentity => One;

    public static BigRational NegativeOne => new(BigInt.MinusOne);

    static int INumberBase<BigRational>.Radix => 2;

    // ---------------------------------------------------------------
    // Predicates
    // ---------------------------------------------------------------

    static bool INumberBase<BigRational>.IsCanonical(BigRational value) => true;
    static bool INumberBase<BigRational>.IsComplexNumber(BigRational value) => false;
    public static bool IsEvenInteger(BigRational value) => value.IsInteger && value._num.IsEven;
    static bool INumberBase<BigRational>.IsFinite(BigRational value) => true;
    static bool INumberBase<BigRational>.IsImaginaryNumber(BigRational value) => false;
    static bool INumberBase<BigRational>.IsInfinity(BigRational value) => false;
    static bool INumberBase<BigRational>.IsInteger(BigRational value) => value.IsInteger;
    static bool INumberBase<BigRational>.IsNaN(BigRational value) => false;
    public static bool IsNegative(BigRational value) => value._num.Sign < 0;
    static bool INumberBase<BigRational>.IsNegativeInfinity(BigRational value) => false;
    static bool INumberBase<BigRational>.IsNormal(BigRational value) => !value.IsZero;
    public static bool IsOddInteger(BigRational value) => value.IsInteger && value._num.IsOdd;
    public static bool IsPositive(BigRational value) => value._num.Sign >= 0;
    static bool INumberBase<BigRational>.IsPositiveInfinity(BigRational value) => false;
    static bool INumberBase<BigRational>.IsRealNumber(BigRational value) => true;
    static bool INumberBase<BigRational>.IsSubnormal(BigRational value) => false;
    static bool INumberBase<BigRational>.IsZero(BigRational value) => value.IsZero;

    static int INumber<BigRational>.Sign(BigRational value) => value.Sign;

    // ---------------------------------------------------------------
    // INumber helpers
    // ---------------------------------------------------------------

    public static BigRational Abs(BigRational value) => value.Abs();

    public static BigRational Clamp(BigRational value, BigRational min, BigRational max)
    {
        if (min > max)
        {
            throw new ArgumentException("min must be less than or equal to max.", nameof(min));
        }
        return value < min ? min : value > max ? max : value;
    }

    public static BigRational CopySign(BigRational value, BigRational sign) => sign.Sign < 0 ? -value.Abs() : value.Abs();

    public static BigRational Max(BigRational x, BigRational y) => x >= y ? x : y;
    public static BigRational Min(BigRational x, BigRational y) => x <= y ? x : y;
    public static BigRational MaxNumber(BigRational x, BigRational y) => Max(x, y);
    public static BigRational MinNumber(BigRational x, BigRational y) => Min(x, y);

    public static BigRational MaxMagnitude(BigRational x, BigRational y)
    {
        int c = x.Abs().CompareTo(y.Abs());
        return c > 0 ? x : c < 0 ? y : x.Sign < 0 ? y : x;
    }

    public static BigRational MinMagnitude(BigRational x, BigRational y)
    {
        int c = x.Abs().CompareTo(y.Abs());
        return c < 0 ? x : c > 0 ? y : x.Sign < 0 ? x : y;
    }

    public static BigRational MaxMagnitudeNumber(BigRational x, BigRational y) => MaxMagnitude(x, y);
    public static BigRational MinMagnitudeNumber(BigRational x, BigRational y) => MinMagnitude(x, y);

    // ---------------------------------------------------------------
    // Remaining operators
    // ---------------------------------------------------------------

    /// <summary>Truncated remainder a - trunc(a/b)*b (sign of the dividend, like C# %).</summary>
    public static BigRational operator %(BigRational a, BigRational b)
    {
        if (b.IsZero)
        {
            throw new DivideByZeroException();
        }
        BigRational t = Truncate(a / b);
        return a - t * b;
    }

    public static BigRational operator ++(BigRational value) => value + One;
    public static BigRational operator --(BigRational value) => value - One;

    // ---------------------------------------------------------------
    // Rounding to integers
    // ---------------------------------------------------------------

    /// <summary>Largest integer not above the value.</summary>
    public static BigRational Floor(BigRational value) => new(RoundToInteger(value, MidpointRounding.ToNegativeInfinity));

    /// <summary>Smallest integer not below the value.</summary>
    public static BigRational Ceiling(BigRational value) => new(RoundToInteger(value, MidpointRounding.ToPositiveInfinity));

    /// <summary>The integer part, rounded toward zero.</summary>
    public static BigRational Truncate(BigRational value) => new(RoundToInteger(value, MidpointRounding.ToZero));

    /// <summary>Nearest integer, ties to even.</summary>
    public static BigRational Round(BigRational value) => new(RoundToInteger(value, MidpointRounding.ToEven));

    /// <summary>Rounds to an integer with the given midpoint / directed mode.</summary>
    public static BigRational Round(BigRational value, MidpointRounding mode) => new(RoundToInteger(value, mode));

    /// <summary>Rounds to <paramref name="digits"/> decimal places (exactly; the result is a rational with denominator dividing 10^digits).</summary>
    public static BigRational Round(BigRational value, int digits, MidpointRounding mode = MidpointRounding.ToEven)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(digits);
        if (digits == 0)
        {
            return Round(value, mode);
        }
        BigInt scale = BigInt.Pow(BigInt.Ten, (ulong)digits);
        return new BigRational(RoundToInteger(value * new BigRational(scale), mode), scale);
    }

    /// <summary>Rounds to a BigInt with the given midpoint / directed mode.</summary>
    public static BigInt RoundToInteger(BigRational value, MidpointRounding mode)
    {
        BigInt num = value._num, den = value.Denominator;
        if (den.IsOne)
        {
            return num;
        }
        switch (mode)
        {
            case MidpointRounding.ToZero:
                return num / den;
            case MidpointRounding.ToNegativeInfinity:
                return BigInt.FloorDivide(num, den);
            case MidpointRounding.ToPositiveInfinity:
                return BigInt.CeilingDivide(num, den);
            case MidpointRounding.ToEven:
            case MidpointRounding.AwayFromZero:
                BigInt.FloorDivRem(num, den, out BigInt q, out BigInt r); // value = q + r/den, 0 < r < den
                int c = (r << 1).CompareTo(den);
                if (c > 0 || (c == 0 && (mode == MidpointRounding.ToEven ? q.IsOdd : num.Sign > 0)))
                {
                    q += BigInt.One;
                }
                return q;
            default:
                throw new ArgumentException("Invalid MidpointRounding value.", nameof(mode));
        }
    }

    // ---------------------------------------------------------------
    // Conversions
    // ---------------------------------------------------------------

    public static implicit operator BigRational(sbyte value) => new((BigInt)value);
    public static implicit operator BigRational(byte value) => new((BigInt)value);
    public static implicit operator BigRational(short value) => new((BigInt)value);
    public static implicit operator BigRational(ushort value) => new((BigInt)value);
    public static implicit operator BigRational(char value) => new((BigInt)value);
    public static implicit operator BigRational(uint value) => new((BigInt)value);
    public static implicit operator BigRational(ulong value) => new((BigInt)value);
    public static implicit operator BigRational(nint value) => new((BigInt)value);
    public static implicit operator BigRational(nuint value) => new((BigInt)value);
    public static implicit operator BigRational(Int128 value) => new((BigInt)value);
    public static implicit operator BigRational(UInt128 value) => new((BigInt)value);
    public static implicit operator BigRational(BigInteger value) => new((BigInt)value);

    /// <summary>Exact conversion from decimal.</summary>
    public static implicit operator BigRational(decimal value) => Real.DecimalToRational(value);

    /// <summary>Exact; throws <see cref="OverflowException"/> for NaN or infinity.</summary>
    public static explicit operator BigRational(float value) => new((double)value);

    /// <summary>Exact; throws <see cref="OverflowException"/> for NaN or infinity.</summary>
    public static explicit operator BigRational(Half value) => new((double)value);

    public static explicit operator sbyte(BigRational value) => (sbyte)(BigInt)value;
    public static explicit operator byte(BigRational value) => (byte)(BigInt)value;
    public static explicit operator short(BigRational value) => (short)(BigInt)value;
    public static explicit operator ushort(BigRational value) => (ushort)(BigInt)value;
    public static explicit operator char(BigRational value) => (char)(BigInt)value;
    public static explicit operator int(BigRational value) => (int)(BigInt)value;
    public static explicit operator uint(BigRational value) => (uint)(BigInt)value;
    public static explicit operator long(BigRational value) => (long)(BigInt)value;
    public static explicit operator ulong(BigRational value) => (ulong)(BigInt)value;
    public static explicit operator nint(BigRational value) => (nint)(BigInt)value;
    public static explicit operator nuint(BigRational value) => (nuint)(BigInt)value;
    public static explicit operator Int128(BigRational value) => (Int128)(BigInt)value;
    public static explicit operator UInt128(BigRational value) => (UInt128)(BigInt)value;
    public static explicit operator BigInteger(BigRational value) => (BigInteger)(BigInt)value;
    public static explicit operator float(BigRational value) => value.ToSingle();
    public static explicit operator Half(BigRational value) => value.ToHalf();

    /// <summary>Nearest decimal (ties to even); throws <see cref="OverflowException"/> when out of range.</summary>
    public static explicit operator decimal(BigRational value) =>
        Real.TryRationalToDecimal(value._num, value.Denominator, out decimal d) ? d : throw new OverflowException("Value does not fit in Decimal.");

    static bool INumberBase<BigRational>.TryConvertFromChecked<TOther>(TOther value, out BigRational result) =>
        TryConvertFrom(value, ConvertMode.Checked, out result);

    static bool INumberBase<BigRational>.TryConvertFromSaturating<TOther>(TOther value, out BigRational result) =>
        TryConvertFrom(value, ConvertMode.Saturating, out result);

    static bool INumberBase<BigRational>.TryConvertFromTruncating<TOther>(TOther value, out BigRational result) =>
        TryConvertFrom(value, ConvertMode.Truncating, out result);

    static bool INumberBase<BigRational>.TryConvertToChecked<TOther>(BigRational value, [MaybeNullWhen(false)] out TOther result) =>
        TryConvertTo(value, ConvertMode.Checked, out result);

    static bool INumberBase<BigRational>.TryConvertToSaturating<TOther>(BigRational value, [MaybeNullWhen(false)] out TOther result) =>
        TryConvertTo(value, ConvertMode.Saturating, out result);

    static bool INumberBase<BigRational>.TryConvertToTruncating<TOther>(BigRational value, [MaybeNullWhen(false)] out TOther result) =>
        TryConvertTo(value, ConvertMode.Truncating, out result);

    private static bool TryConvertFrom<TOther>(TOther value, ConvertMode mode, out BigRational result)
        where TOther : INumberBase<TOther>
    {
        if (BigInt.TryFromIntegerPrimitive(value, out BigInt i))
        {
            result = new BigRational(i);
            return true;
        }
        if (BigInt.IsBinaryFloatPrimitive<TOther>())
        {
            double d = double.CreateTruncating(value);
            if (double.IsNaN(d) && mode != ConvertMode.Checked)
            {
                result = default;
                return true;
            }
            result = new BigRational(d); // throws OverflowException for NaN / infinity
            return true;
        }
        if (typeof(TOther) == typeof(decimal))
        {
            result = (decimal)(object)value;
            return true;
        }
        if (typeof(TOther) == typeof(BigInteger))
        {
            result = (BigInteger)(object)value;
            return true;
        }
        if (typeof(TOther) == typeof(BigInt))
        {
            result = (BigInt)(object)value;
            return true;
        }
        result = default;
        return false;
    }

    private static bool TryConvertTo<TOther>(BigRational value, ConvertMode mode, [MaybeNullWhen(false)] out TOther result)
        where TOther : INumberBase<TOther>
    {
        if (BigInt.IsIntegerPrimitive<TOther>())
        {
            result = BigInt.ConvertToIntegerPrimitive<TOther>((BigInt)value, mode);
            return true;
        }
        if (BigInt.IsBinaryFloatPrimitive<TOther>())
        {
            IeeeFormat f = typeof(TOther) == typeof(double) ? IeeeFormat.Binary64
                : typeof(TOther) == typeof(float) ? IeeeFormat.Binary32 : IeeeFormat.Binary16;
            result = TOther.CreateTruncating(value.ToIeee(f));
            return true;
        }
        if (typeof(TOther) == typeof(decimal))
        {
            if (!Real.TryRationalToDecimal(value._num, value.Denominator, out decimal d))
            {
                if (mode == ConvertMode.Checked)
                {
                    throw new OverflowException("Value does not fit in Decimal.");
                }
                d = value.Sign < 0 ? decimal.MinValue : decimal.MaxValue;
            }
            result = (TOther)(object)d;
            return true;
        }
        if (typeof(TOther) == typeof(BigInteger))
        {
            result = (TOther)(object)(BigInteger)(BigInt)value;
            return true;
        }
        if (typeof(TOther) == typeof(BigInt))
        {
            result = (TOther)(object)(BigInt)value;
            return true;
        }
        result = default;
        return false;
    }

    // ---------------------------------------------------------------
    // NumberStyles parsing / UTF-8
    // ---------------------------------------------------------------

    /// <summary>Parses "num/den" (each side an integer in the given style) or a decimal number allowed by the style.</summary>
    public static BigRational Parse(ReadOnlySpan<char> s, NumberStyles style, IFormatProvider? provider) =>
        TryParse(s, style, provider, out BigRational r) ? r : throw new FormatException("The value could not be parsed as a BigRational.");

    public static BigRational Parse(string s, NumberStyles style, IFormatProvider? provider)
    {
        ArgumentNullException.ThrowIfNull(s);
        return Parse(s.AsSpan(), style, provider);
    }

    public static bool TryParse(string? s, NumberStyles style, IFormatProvider? provider, out BigRational result) =>
        TryParse(s.AsSpan(), style, provider, out result);

    public static bool TryParse(ReadOnlySpan<char> s, NumberStyles style, IFormatProvider? provider, out BigRational result)
    {
        result = default;
        int slash = s.IndexOf('/');
        if (slash >= 0)
        {
            NumberStyles left = style & ~NumberStyles.AllowTrailingWhite, right = style & ~NumberStyles.AllowLeadingWhite;
            if (!BigInt.TryParse(s.Slice(0, slash), left, provider, out BigInt num) ||
                !BigInt.TryParse(s.Slice(slash + 1), right, provider, out BigInt den) ||
                den.IsZero)
            {
                return false;
            }
            result = new BigRational(num, den);
            return true;
        }
        if ((style & (NumberStyles.AllowHexSpecifier | NumberStyles.AllowBinarySpecifier)) != 0)
        {
            if (!BigInt.TryParse(s, style, provider, out BigInt i))
            {
                return false;
            }
            result = new BigRational(i);
            return true;
        }
        return NumberText.TryNormalize(s, style, provider, out string text, out _) &&
            TryParse(text, CultureInfo.InvariantCulture, out result);
    }

    public bool TryFormat(Span<byte> utf8Destination, out int bytesWritten, ReadOnlySpan<char> format, IFormatProvider? provider) =>
        NumberText.TryWriteUtf8(ToString(), utf8Destination, out bytesWritten);

    public static BigRational Parse(ReadOnlySpan<byte> utf8Text, IFormatProvider? provider) =>
        TryParse(utf8Text, provider, out BigRational r) ? r : throw new FormatException("The value could not be parsed as a BigRational.");

    public static bool TryParse(ReadOnlySpan<byte> utf8Text, IFormatProvider? provider, out BigRational result) =>
        NumberText.TryParseUtf8(utf8Text, (ReadOnlySpan<char> s, out BigRational r) => TryParse(s, provider, out r), out result);

    public static BigRational Parse(ReadOnlySpan<byte> utf8Text, NumberStyles style, IFormatProvider? provider) =>
        TryParse(utf8Text, style, provider, out BigRational r) ? r : throw new FormatException("The value could not be parsed as a BigRational.");

    public static bool TryParse(ReadOnlySpan<byte> utf8Text, NumberStyles style, IFormatProvider? provider, out BigRational result) =>
        NumberText.TryParseUtf8(utf8Text, (ReadOnlySpan<char> s, out BigRational r) => TryParse(s, style, provider, out r), out result);

    // ---------------------------------------------------------------
    // Public CreateChecked / CreateSaturating / CreateTruncating (like BigInteger)
    // ---------------------------------------------------------------

    public static BigRational CreateChecked<TOther>(TOther value) where TOther : INumberBase<TOther>
    {
        if (typeof(TOther) == typeof(BigRational))
        {
            return (BigRational)(object)value;
        }
        if (!TryConvertFrom(value, ConvertMode.Checked, out BigRational result) && !TOther.TryConvertToChecked(value, out result))
        {
            throw new NotSupportedException($"Conversion from {typeof(TOther).Name} to BigRational is not supported.");
        }
        return result;
    }

    public static BigRational CreateSaturating<TOther>(TOther value) where TOther : INumberBase<TOther>
    {
        if (typeof(TOther) == typeof(BigRational))
        {
            return (BigRational)(object)value;
        }
        if (!TryConvertFrom(value, ConvertMode.Saturating, out BigRational result) && !TOther.TryConvertToSaturating(value, out result))
        {
            throw new NotSupportedException($"Conversion from {typeof(TOther).Name} to BigRational is not supported.");
        }
        return result;
    }

    public static BigRational CreateTruncating<TOther>(TOther value) where TOther : INumberBase<TOther>
    {
        if (typeof(TOther) == typeof(BigRational))
        {
            return (BigRational)(object)value;
        }
        if (!TryConvertFrom(value, ConvertMode.Truncating, out BigRational result) && !TOther.TryConvertToTruncating(value, out result))
        {
            throw new NotSupportedException($"Conversion from {typeof(TOther).Name} to BigRational is not supported.");
        }
        return result;
    }
}