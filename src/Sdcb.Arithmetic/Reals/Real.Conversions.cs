using System.Diagnostics.CodeAnalysis;
using System.Numerics;

namespace Sdcb.Arithmetic;

public readonly partial struct Real
{
    // ---------------------------------------------------------------
    // Implicit (exact) conversions. Integers get max(DefaultPrecision, bit length) bits;
    // binary floating-point values get max(DefaultPrecision, their significand width).
    // ---------------------------------------------------------------

    public static implicit operator Real(sbyte value) => FromBigInt(value);
    public static implicit operator Real(byte value) => FromBigInt(value);
    public static implicit operator Real(short value) => FromBigInt(value);
    public static implicit operator Real(ushort value) => FromBigInt(value);
    public static implicit operator Real(char value) => FromBigInt(value);
    public static implicit operator Real(int value) => FromBigInt(value);
    public static implicit operator Real(uint value) => FromBigInt(value);
    public static implicit operator Real(long value) => FromBigInt(value);
    public static implicit operator Real(ulong value) => FromBigInt(value);
    public static implicit operator Real(nint value) => FromBigInt(value);
    public static implicit operator Real(nuint value) => FromBigInt(value);
    public static implicit operator Real(Int128 value) => FromBigInt(value);
    public static implicit operator Real(UInt128 value) => FromBigInt(value);
    public static implicit operator Real(BigInt value) => FromBigInt(value);
    public static implicit operator Real(BigInteger value) => FromBigInt(value);
    public static implicit operator Real(System.Half value) => FromDouble((double)value, Math.Max(DefaultPrecision, 11));
    public static implicit operator Real(float value) => FromDouble(value, Math.Max(DefaultPrecision, 24));
    public static implicit operator Real(double value) => FromDouble(value, Math.Max(DefaultPrecision, 53));

    // ---------------------------------------------------------------
    // Explicit (rounding) conversions
    // ---------------------------------------------------------------

    /// <summary>Rounds the decimal to DefaultPrecision bits (nearest-even).</summary>
    public static explicit operator Real(decimal value) => FromDecimal(value);

    /// <summary>Rounds the rational to DefaultPrecision bits (nearest-even).</summary>
    public static explicit operator Real(BigRational value) => FromRational(value);

    /// <summary>A decimal correctly rounded to the given precision.</summary>
    public static Real FromDecimal(decimal value, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        Real r = FromRational(DecimalToRational(value), precision, rounding);
        return r.IsZero && decimal.IsNegative(value) ? -r : r;
    }

    public static explicit operator sbyte(Real value) => (sbyte)value.TruncateForInteger();
    public static explicit operator byte(Real value) => (byte)value.TruncateForInteger();
    public static explicit operator short(Real value) => (short)value.TruncateForInteger();
    public static explicit operator ushort(Real value) => (ushort)value.TruncateForInteger();
    public static explicit operator char(Real value) => (char)value.TruncateForInteger();
    public static explicit operator int(Real value) => (int)value.TruncateForInteger();
    public static explicit operator uint(Real value) => (uint)value.TruncateForInteger();
    public static explicit operator long(Real value) => (long)value.TruncateForInteger();
    public static explicit operator ulong(Real value) => (ulong)value.TruncateForInteger();
    public static explicit operator nint(Real value) => (nint)value.TruncateForInteger();
    public static explicit operator nuint(Real value) => (nuint)value.TruncateForInteger();
    public static explicit operator Int128(Real value) => (Int128)value.TruncateForInteger();
    public static explicit operator UInt128(Real value) => (UInt128)value.TruncateForInteger();

    /// <summary>Truncates toward zero; throws <see cref="OverflowException"/> for NaN or infinity.</summary>
    public static explicit operator BigInt(Real value) => value.ToBigInt();

    /// <summary>Truncates toward zero; throws <see cref="OverflowException"/> for NaN or infinity.</summary>
    public static explicit operator BigInteger(Real value) => (BigInteger)value.ToBigInt();

    /// <summary>The exact value; throws <see cref="OverflowException"/> for NaN or infinity.</summary>
    public static explicit operator BigRational(Real value) => value.ToRational();

    public static explicit operator double(Real value) => value.ToDouble();
    public static explicit operator float(Real value) => value.ToSingle();
    public static explicit operator System.Half(Real value) => value.ToHalf();

    /// <summary>Nearest decimal (ties to even); throws <see cref="OverflowException"/> for NaN, infinity or out-of-range values.</summary>
    public static explicit operator decimal(Real value) => value.ToDecimal();

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

    static bool INumberBase<Real>.TryConvertFromChecked<TOther>(TOther value, out Real result) => TryConvertFrom(value, out result);
    static bool INumberBase<Real>.TryConvertFromSaturating<TOther>(TOther value, out Real result) => TryConvertFrom(value, out result);
    static bool INumberBase<Real>.TryConvertFromTruncating<TOther>(TOther value, out Real result) => TryConvertFrom(value, out result);

    static bool INumberBase<Real>.TryConvertToChecked<TOther>(Real value, [MaybeNullWhen(false)] out TOther result) =>
        TryConvertTo(value, ConvertMode.Checked, out result);

    static bool INumberBase<Real>.TryConvertToSaturating<TOther>(Real value, [MaybeNullWhen(false)] out TOther result) =>
        TryConvertTo(value, ConvertMode.Saturating, out result);

    static bool INumberBase<Real>.TryConvertToTruncating<TOther>(Real value, [MaybeNullWhen(false)] out TOther result) =>
        TryConvertTo(value, ConvertMode.Truncating, out result);

    private static bool TryConvertFrom<TOther>(TOther value, out Real result)
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

    private static bool TryConvertTo<TOther>(Real value, ConvertMode mode, [MaybeNullWhen(false)] out TOther result)
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
            result = TOther.CreateTruncating(value.ToIeee(f, RealRounding.NearestEven));
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

    public static Real CreateChecked<TOther>(TOther value) where TOther : INumberBase<TOther>
    {
        if (typeof(TOther) == typeof(Real))
        {
            return (Real)(object)value;
        }
        if (!TryConvertFrom(value, out Real result) && !TOther.TryConvertToChecked(value, out result))
        {
            throw new NotSupportedException($"Conversion from {typeof(TOther).Name} to Real is not supported.");
        }
        return result;
    }

    public static Real CreateSaturating<TOther>(TOther value) where TOther : INumberBase<TOther>
    {
        if (typeof(TOther) == typeof(Real))
        {
            return (Real)(object)value;
        }
        if (!TryConvertFrom(value, out Real result) && !TOther.TryConvertToSaturating(value, out result))
        {
            throw new NotSupportedException($"Conversion from {typeof(TOther).Name} to Real is not supported.");
        }
        return result;
    }

    public static Real CreateTruncating<TOther>(TOther value) where TOther : INumberBase<TOther>
    {
        if (typeof(TOther) == typeof(Real))
        {
            return (Real)(object)value;
        }
        if (!TryConvertFrom(value, out Real result) && !TOther.TryConvertToTruncating(value, out result))
        {
            throw new NotSupportedException($"Conversion from {typeof(TOther).Name} to Real is not supported.");
        }
        return result;
    }
}