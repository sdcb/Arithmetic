using System.Diagnostics.CodeAnalysis;
using System.Numerics;

namespace Sdcb.Arithmetic;

public readonly partial struct BigFloat
{
    // ---------------------------------------------------------------
    // Implicit (exact) conversions. Integers get max(DefaultPrecision, bit length) bits;
    // binary floating-point values get max(DefaultPrecision, their significand width).
    // ---------------------------------------------------------------

    public static implicit operator BigFloat(sbyte value) => FromBigInt(value);
    public static implicit operator BigFloat(byte value) => FromBigInt(value);
    public static implicit operator BigFloat(short value) => FromBigInt(value);
    public static implicit operator BigFloat(ushort value) => FromBigInt(value);
    public static implicit operator BigFloat(char value) => FromBigInt(value);
    public static implicit operator BigFloat(int value) => FromBigInt(value);
    public static implicit operator BigFloat(uint value) => FromBigInt(value);
    public static implicit operator BigFloat(long value) => FromBigInt(value);
    public static implicit operator BigFloat(ulong value) => FromBigInt(value);
    public static implicit operator BigFloat(nint value) => FromBigInt(value);
    public static implicit operator BigFloat(nuint value) => FromBigInt(value);
    public static implicit operator BigFloat(Int128 value) => FromBigInt(value);
    public static implicit operator BigFloat(UInt128 value) => FromBigInt(value);
    public static implicit operator BigFloat(BigInt value) => FromBigInt(value);
    public static implicit operator BigFloat(BigInteger value) => FromBigInt(value);
    public static implicit operator BigFloat(System.Half value) => FromDouble((double)value, Math.Max(DefaultPrecision, 11));
    public static implicit operator BigFloat(float value) => FromDouble(value, Math.Max(DefaultPrecision, 24));
    public static implicit operator BigFloat(double value) => FromDouble(value, Math.Max(DefaultPrecision, 53));

    // ---------------------------------------------------------------
    // Explicit (rounding) conversions
    // ---------------------------------------------------------------

    /// <summary>Rounds the decimal to DefaultPrecision bits (nearest-even).</summary>
    public static explicit operator BigFloat(decimal value) => FromDecimal(value);

    /// <summary>Rounds the rational to DefaultPrecision bits (nearest-even).</summary>
    public static explicit operator BigFloat(BigRational value) => FromRational(value);

    /// <summary>A decimal correctly rounded to the given precision.</summary>
    public static BigFloat FromDecimal(decimal value, int precision = 0, RoundingMode rounding = RoundingMode.NearestEven)
    {
        BigFloat r = FromRational(DecimalToRational(value), precision, rounding);
        return r.IsZero && decimal.IsNegative(value) ? -r : r;
    }

    public static explicit operator sbyte(BigFloat value) => (sbyte)value.TruncateForInteger();
    public static explicit operator byte(BigFloat value) => (byte)value.TruncateForInteger();
    public static explicit operator short(BigFloat value) => (short)value.TruncateForInteger();
    public static explicit operator ushort(BigFloat value) => (ushort)value.TruncateForInteger();
    public static explicit operator char(BigFloat value) => (char)value.TruncateForInteger();
    public static explicit operator int(BigFloat value) => (int)value.TruncateForInteger();
    public static explicit operator uint(BigFloat value) => (uint)value.TruncateForInteger();
    public static explicit operator long(BigFloat value) => (long)value.TruncateForInteger();
    public static explicit operator ulong(BigFloat value) => (ulong)value.TruncateForInteger();
    public static explicit operator nint(BigFloat value) => (nint)value.TruncateForInteger();
    public static explicit operator nuint(BigFloat value) => (nuint)value.TruncateForInteger();
    public static explicit operator Int128(BigFloat value) => (Int128)value.TruncateForInteger();
    public static explicit operator UInt128(BigFloat value) => (UInt128)value.TruncateForInteger();

    /// <summary>Truncates toward zero; throws <see cref="OverflowException"/> for NaN or infinity.</summary>
    public static explicit operator BigInt(BigFloat value) => value.ToBigInt();

    /// <summary>Truncates toward zero; throws <see cref="OverflowException"/> for NaN or infinity.</summary>
    public static explicit operator BigInteger(BigFloat value) => (BigInteger)value.ToBigInt();

    /// <summary>The exact value; throws <see cref="OverflowException"/> for NaN or infinity.</summary>
    public static explicit operator BigRational(BigFloat value) => value.ToRational();

    public static explicit operator double(BigFloat value) => value.ToDouble();
    public static explicit operator float(BigFloat value) => value.ToSingle();
    public static explicit operator System.Half(BigFloat value) => value.ToHalf();

    /// <summary>Nearest decimal (ties to even); throws <see cref="OverflowException"/> for NaN, infinity or out-of-range values.</summary>
    public static explicit operator decimal(BigFloat value) => value.ToDecimal();

    /// <summary>Nearest decimal (ties to even); throws <see cref="OverflowException"/> for NaN, infinity or out-of-range values.</summary>
    public decimal ToDecimal() =>
        TryToDecimal(out decimal d) ? d : throw new OverflowException("Value does not fit in Decimal.");

    private bool TryToDecimal(out decimal result)
    {
        result = 0m;
        switch (_kind)
        {
            case Kind.NaN:
            case Kind.Inf:
                return false;
            case Kind.Zero:
                return true;
        }
        long msb = BinaryExponent;
        if (msb > 97)
        {
            return false;
        }
        if (msb < -100)
        {
            return true; // far below 10^-28 / 2
        }
        BigRational q = ToRational();
        return TryRationalToDecimal(q.Numerator, q.Denominator, out result);
    }

    /// <summary>ToBigInt for primitive-integer targets: avoids materializing astronomically large integers.</summary>
    private BigInt TruncateForInteger()
    {
        if (_kind == Kind.Normal && BinaryExponent > 130)
        {
            throw new OverflowException("Value is too large for the target integer type.");
        }
        return ToBigInt();
    }

    // ---------------------------------------------------------------
    // Generic-math conversions
    // ---------------------------------------------------------------

    static bool INumberBase<BigFloat>.TryConvertFromChecked<TOther>(TOther value, out BigFloat result) => TryConvertFrom(value, out result);
    static bool INumberBase<BigFloat>.TryConvertFromSaturating<TOther>(TOther value, out BigFloat result) => TryConvertFrom(value, out result);
    static bool INumberBase<BigFloat>.TryConvertFromTruncating<TOther>(TOther value, out BigFloat result) => TryConvertFrom(value, out result);

    static bool INumberBase<BigFloat>.TryConvertToChecked<TOther>(BigFloat value, [MaybeNullWhen(false)] out TOther result) =>
        TryConvertTo(value, ConvertMode.Checked, out result);

    static bool INumberBase<BigFloat>.TryConvertToSaturating<TOther>(BigFloat value, [MaybeNullWhen(false)] out TOther result) =>
        TryConvertTo(value, ConvertMode.Saturating, out result);

    static bool INumberBase<BigFloat>.TryConvertToTruncating<TOther>(BigFloat value, [MaybeNullWhen(false)] out TOther result) =>
        TryConvertTo(value, ConvertMode.Truncating, out result);

    private static bool TryConvertFrom<TOther>(TOther value, out BigFloat result)
        where TOther : INumberBase<TOther>
    {
        if (BigInt.TryFromIntegerPrimitive(value, out BigInt i))
        {
            result = FromBigInt(i);
            return true;
        }
        if (typeof(TOther) == typeof(double))
        {
            result = (double)(object)value;
            return true;
        }
        if (typeof(TOther) == typeof(float))
        {
            result = (float)(object)value;
            return true;
        }
        if (typeof(TOther) == typeof(System.Half))
        {
            result = (System.Half)(object)value;
            return true;
        }
        if (typeof(TOther) == typeof(decimal))
        {
            result = FromDecimal((decimal)(object)value);
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
        if (typeof(TOther) == typeof(BigRational))
        {
            result = FromRational((BigRational)(object)value);
            return true;
        }
        result = default;
        return false;
    }

    private static bool TryConvertTo<TOther>(BigFloat value, ConvertMode mode, [MaybeNullWhen(false)] out TOther result)
        where TOther : INumberBase<TOther>
    {
        bool special = value._kind is Kind.NaN or Kind.Inf;
        if (BigInt.IsIntegerPrimitive<TOther>())
        {
            bool huge = value._kind == Kind.Normal && value.BinaryExponent > 130;
            if (special || huge)
            {
                if (mode == ConvertMode.Checked)
                {
                    throw new OverflowException($"Value does not fit in {typeof(TOther).Name}.");
                }
                if (mode == ConvertMode.Truncating && huge)
                {
                    // Two's-complement truncation keeps the low 128 bits, all zero once the lsb is >= 2^128.
                    result = value._exp >= 128 ? TOther.Zero : BigInt.ConvertToIntegerPrimitive<TOther>(value.ToBigInt(), mode);
                    return true;
                }
                // NaN -> 0, ±Inf / huge -> saturate (the BCL rule for double sources).
                double d = special ? value.ToDouble() : value.IsNeg ? double.NegativeInfinity : double.PositiveInfinity;
                result = mode == ConvertMode.Saturating ? TOther.CreateSaturating(d) : TOther.CreateTruncating(d);
                return true;
            }
            result = BigInt.ConvertToIntegerPrimitive<TOther>(value.ToBigInt(), mode);
            return true;
        }
        if (BigInt.IsBinaryFloatPrimitive<TOther>())
        {
            IeeeFormat f = typeof(TOther) == typeof(double) ? IeeeFormat.Binary64
                : typeof(TOther) == typeof(float) ? IeeeFormat.Binary32 : IeeeFormat.Binary16;
            result = TOther.CreateTruncating(value.ToIeee(f, RoundingMode.NearestEven));
            return true;
        }
        if (typeof(TOther) == typeof(decimal))
        {
            if (!value.TryToDecimal(out decimal m))
            {
                if (mode == ConvertMode.Checked)
                {
                    throw new OverflowException("Value does not fit in Decimal.");
                }
                m = value._kind == Kind.NaN ? 0m : value.IsNeg ? decimal.MinValue : decimal.MaxValue;
            }
            result = (TOther)(object)m;
            return true;
        }
        if (typeof(TOther) == typeof(BigInteger) || typeof(TOther) == typeof(BigInt) || typeof(TOther) == typeof(BigRational))
        {
            if (value._kind == Kind.NaN && mode != ConvertMode.Checked)
            {
                result = TOther.Zero;
                return true;
            }
            if (typeof(TOther) == typeof(BigRational))
            {
                result = (TOther)(object)value.ToRational();
            }
            else
            {
                BigInt b = value.ToBigInt();
                result = typeof(TOther) == typeof(BigInt) ? (TOther)(object)b : (TOther)(object)(BigInteger)b;
            }
            return true;
        }
        result = default;
        return false;
    }

    // ---------------------------------------------------------------
    // Public CreateChecked / CreateSaturating / CreateTruncating (like BigInteger)
    // ---------------------------------------------------------------

    public static BigFloat CreateChecked<TOther>(TOther value) where TOther : INumberBase<TOther>
    {
        if (typeof(TOther) == typeof(BigFloat))
        {
            return (BigFloat)(object)value;
        }
        if (!TryConvertFrom(value, out BigFloat result) && !TOther.TryConvertToChecked(value, out result))
        {
            throw new NotSupportedException($"Conversion from {typeof(TOther).Name} to BigFloat is not supported.");
        }
        return result;
    }

    public static BigFloat CreateSaturating<TOther>(TOther value) where TOther : INumberBase<TOther>
    {
        if (typeof(TOther) == typeof(BigFloat))
        {
            return (BigFloat)(object)value;
        }
        if (!TryConvertFrom(value, out BigFloat result) && !TOther.TryConvertToSaturating(value, out result))
        {
            throw new NotSupportedException($"Conversion from {typeof(TOther).Name} to BigFloat is not supported.");
        }
        return result;
    }

    public static BigFloat CreateTruncating<TOther>(TOther value) where TOther : INumberBase<TOther>
    {
        if (typeof(TOther) == typeof(BigFloat))
        {
            return (BigFloat)(object)value;
        }
        if (!TryConvertFrom(value, out BigFloat result) && !TOther.TryConvertToTruncating(value, out result))
        {
            throw new NotSupportedException($"Conversion from {typeof(TOther).Name} to BigFloat is not supported.");
        }
        return result;
    }
}