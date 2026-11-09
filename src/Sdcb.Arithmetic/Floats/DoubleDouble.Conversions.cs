using System.Diagnostics.CodeAnalysis;
using System.Numerics;
using System.Runtime.CompilerServices;

namespace Sdcb.Arithmetic;

public readonly partial struct DoubleDouble
{
    /// <summary>Enough bits for any exact double-double (2^1024 down to 2^-1074) plus guard bits.</summary>
    private const int BigFloatExactPrecision = 1100;

    // ---------------------------------------------------------------
    // Implicit (exact) conversions
    // ---------------------------------------------------------------

    public static implicit operator DoubleDouble(double value) => new(value, 0.0, normalized: true);
    public static implicit operator DoubleDouble(float value) => new(value, 0.0, normalized: true);
    public static implicit operator DoubleDouble(System.Half value) => new((double)value, 0.0, normalized: true);
    public static implicit operator DoubleDouble(sbyte value) => new(value, 0.0, normalized: true);
    public static implicit operator DoubleDouble(byte value) => new(value, 0.0, normalized: true);
    public static implicit operator DoubleDouble(short value) => new(value, 0.0, normalized: true);
    public static implicit operator DoubleDouble(ushort value) => new(value, 0.0, normalized: true);
    public static implicit operator DoubleDouble(char value) => new(value, 0.0, normalized: true);
    public static implicit operator DoubleDouble(int value) => new(value, 0.0, normalized: true);
    public static implicit operator DoubleDouble(uint value) => new(value, 0.0, normalized: true);

    /// <summary>Exact: a 64-bit integer always fits in 106 bits.</summary>
    public static implicit operator DoubleDouble(long value)
    {
        double hi = value;
        // hi may round up to 2^63, which has no long representation.
        double lo = hi >= 9.2233720368547758E18 ? (double)(value - long.MaxValue) - 1.0 : (double)(value - (long)hi);
        return new DoubleDouble(hi, lo, normalized: true);
    }

    /// <summary>Exact: a 64-bit integer always fits in 106 bits.</summary>
    public static implicit operator DoubleDouble(ulong value)
    {
        double hi = value;
        // hi may round up to 2^64; otherwise the wrapped difference is a small signed integer.
        double lo = hi >= 1.8446744073709552E19 ? -(double)(ulong.MaxValue - value) - 1.0 : (double)unchecked((long)(value - (ulong)hi));
        return new DoubleDouble(hi, lo, normalized: true);
    }

    public static implicit operator DoubleDouble(nint value) => (long)value;
    public static implicit operator DoubleDouble(nuint value) => (ulong)value;

    // ---------------------------------------------------------------
    // Explicit (rounding) conversions into DoubleDouble
    // ---------------------------------------------------------------

    /// <summary>Nearest double-double (the 128-bit value may need more than 106 bits).</summary>
    public static explicit operator DoubleDouble(Int128 value) => FromBigFloat(BigFloat.FromBigInt(new BigInt(value), 0, BigFloatExactPrecision));

    /// <summary>Nearest double-double (the 128-bit value may need more than 106 bits).</summary>
    public static explicit operator DoubleDouble(UInt128 value) => FromBigFloat(BigFloat.FromBigInt(new BigInt(value), 0, BigFloatExactPrecision));

    /// <summary>Nearest double-double of the decimal value.</summary>
    public static explicit operator DoubleDouble(decimal value)
    {
        DoubleDouble r = FromBigFloat(BigFloat.FromDecimal(value, 256));
        return r._hi == 0.0 && decimal.IsNegative(value) ? CopySign(r, NegativeOne) : r;
    }

    /// <summary>Nearest double-double; ±∞ beyond the double range.</summary>
    public static explicit operator DoubleDouble(BigInt value) => FromBigFloat(BigFloat.FromBigInt(value));

    /// <summary>Nearest double-double; ±∞ beyond the double range.</summary>
    public static explicit operator DoubleDouble(BigInteger value) => FromBigFloat(BigFloat.FromBigInt((BigInt)value));

    /// <summary>Nearest double-double of the rational.</summary>
    public static explicit operator DoubleDouble(BigRational value) => FromBigFloat(BigFloat.FromRational(value, 256));

    /// <summary>Rounds to the nearest double-double (see <see cref="FromBigFloat(BigFloat)"/>).</summary>
    public static explicit operator DoubleDouble(BigFloat value) => FromBigFloat(value);

    /// <summary>
    /// Rounds a <see cref="BigFloat"/> to double-double: Hi = RN(value), Lo = RN(value - Hi).
    /// Values beyond <see cref="MaxValue"/> become ±∞ (or ±MaxValue when they would round into it).
    /// </summary>
    public static DoubleDouble FromBigFloat(BigFloat value)
    {
        double hi = value.ToDouble();
        if (!double.IsFinite(hi) || hi == 0.0)
        {
            return new DoubleDouble(hi, 0.0, normalized: true);
        }
        double lo = BigFloat.Subtract(value, BigFloat.FromDouble(hi, 53), 53).ToDouble();
        double s = hi + lo;
        if (!double.IsFinite(s))
        {
            // lo reached exactly half an ulp of MaxValue: step back into range.
            return new DoubleDouble(hi, Math.BitDecrement(Math.Abs(lo)) * Math.Sign(lo), normalized: true);
        }
        return new DoubleDouble(hi, lo);
    }

    // ---------------------------------------------------------------
    // Exact conversion to BigFloat
    // ---------------------------------------------------------------

    /// <summary>
    /// Exact: Hi + Lo as a <see cref="BigFloat"/> with precision max(<see cref="BigFloat.DefaultPrecision"/>,
    /// the bits the sum needs).
    /// </summary>
    public static implicit operator BigFloat(DoubleDouble value) => value.ToBigFloatExact(BigFloat.DefaultPrecision);

    /// <summary>The exact value as a <see cref="BigFloat"/> (same as the implicit conversion).</summary>
    public BigFloat ToBigFloat() => ToBigFloatExact(BigFloat.DefaultPrecision);

    /// <summary>Hi + Lo correctly rounded to float (nearest-even by default).</summary>
    public float ToSingle(RoundingMode rounding = RoundingMode.NearestEven) =>
        rounding == RoundingMode.NearestEven ? (float)this : ToBigFloatExact().ToSingle(rounding);

    /// <summary>Hi + Lo correctly rounded to <see cref="System.Half"/> (nearest-even by default).</summary>
    public System.Half ToHalf(RoundingMode rounding = RoundingMode.NearestEven) =>
        rounding == RoundingMode.NearestEven ? (System.Half)this : ToBigFloatExact().ToHalf(rounding);

    /// <summary>The exact value as a BigFloat with at least <paramref name="minPrecision"/> bits.</summary>
    internal BigFloat ToBigFloatExact(int minPrecision = 2)
    {
        if (!double.IsFinite(_hi) || _hi == 0.0 || _lo == 0.0)
        {
            return BigFloat.FromDouble(_hi, Math.Max(minPrecision, 53));
        }
        BigInt mh = Significand(_hi, out long eh);
        BigInt ml = Significand(_lo, out long el);
        long e = Math.Min(eh, el);
        BigInt m = mh.ShiftLeft(eh - e) + ml.ShiftLeft(el - e);
        int bits = (int)m.Abs().GetBitLength();
        return BigFloat.FromBigInt(m, e, Math.Max(minPrecision, bits));
    }

    /// <summary>Signed integer significand and lsb exponent of a finite nonzero double.</summary>
    private static BigInt Significand(double d, out long exponent)
    {
        ulong bits = BitConverter.DoubleToUInt64Bits(d);
        int biased = (int)((bits >> 52) & 0x7FF);
        ulong mant = bits & 0xFFFFFFFFFFFFFUL;
        if (biased == 0)
        {
            exponent = -1074;
        }
        else
        {
            mant |= 1UL << 52;
            exponent = biased - 1075;
        }
        BigInt m = new(mant);
        return (long)bits < 0 ? -m : m;
    }

    // ---------------------------------------------------------------
    // Explicit conversions out of DoubleDouble
    // ---------------------------------------------------------------

    /// <summary>Hi + Lo rounded to double (that is <see cref="Hi"/>).</summary>
    public static explicit operator double(DoubleDouble value) => value._hi;

    /// <summary>
    /// Hi + Lo correctly rounded to double. Nearest-even is <see cref="Hi"/> itself; the directed
    /// modes round the exact value.
    /// </summary>
    public double ToDouble(RoundingMode rounding = RoundingMode.NearestEven) =>
        rounding == RoundingMode.NearestEven ? _hi : ToBigFloatExact().ToDouble(rounding);

    /// <summary>Correctly rounded to float (no double rounding through <see cref="Hi"/>).</summary>
    public static explicit operator float(DoubleDouble value)
    {
        float f = (float)value._hi;
        if (value._lo == 0.0 || !IsFloatMidpoint(value._hi, f))
        {
            return f;
        }
        return value.ToBigFloatExact().ToSingle();
    }

    /// <summary>Correctly rounded to Half.</summary>
    public static explicit operator System.Half(DoubleDouble value) =>
        value._lo == 0.0 ? (System.Half)value._hi : value.ToBigFloatExact().ToHalf();

    private static bool IsFloatMidpoint(double d, float f)
    {
        float other = d > f ? MathF.BitIncrement(f) : MathF.BitDecrement(f);
        return ((double)f + other) * 0.5 == d;
    }

    /// <summary>Nearest decimal; throws <see cref="OverflowException"/> for NaN, ±∞ or out-of-range values.</summary>
    public static explicit operator decimal(DoubleDouble value)
    {
        if (!double.IsFinite(value._hi))
        {
            throw new OverflowException("Value does not fit in Decimal.");
        }
        return value.ToBigFloatExact().ToDecimal();
    }

    /// <summary>Truncates toward zero; throws <see cref="OverflowException"/> for NaN or ±∞.</summary>
    public static explicit operator BigInt(DoubleDouble value) => value.ToBigIntTruncated();

    /// <summary>Truncates toward zero; throws <see cref="OverflowException"/> for NaN or ±∞.</summary>
    public static explicit operator BigInteger(DoubleDouble value) => (BigInteger)value.ToBigIntTruncated();

    /// <summary>The exact value; throws <see cref="OverflowException"/> for NaN or ±∞.</summary>
    public static explicit operator BigRational(DoubleDouble value)
    {
        if (!double.IsFinite(value._hi))
        {
            throw new OverflowException("Cannot convert NaN or Infinity to BigRational.");
        }
        return value.ToBigFloatExact().ToRational();
    }

    private BigInt ToBigIntTruncated()
    {
        if (!double.IsFinite(_hi))
        {
            throw new OverflowException("Cannot convert NaN or Infinity to an integer.");
        }
        DoubleDouble t = Truncate(this);
        return new BigInt(t._hi) + new BigInt(t._lo);
    }

    // Integer conversions truncate toward zero and throw OverflowException for NaN, ±∞ and
    // out-of-range values (like BigFloat).
    public static explicit operator sbyte(DoubleDouble value) => ToInteger<sbyte>(value, ConvertMode.Checked);
    public static explicit operator byte(DoubleDouble value) => ToInteger<byte>(value, ConvertMode.Checked);
    public static explicit operator short(DoubleDouble value) => ToInteger<short>(value, ConvertMode.Checked);
    public static explicit operator ushort(DoubleDouble value) => ToInteger<ushort>(value, ConvertMode.Checked);
    public static explicit operator char(DoubleDouble value) => ToInteger<char>(value, ConvertMode.Checked);
    public static explicit operator int(DoubleDouble value) => ToInteger<int>(value, ConvertMode.Checked);
    public static explicit operator uint(DoubleDouble value) => ToInteger<uint>(value, ConvertMode.Checked);
    public static explicit operator long(DoubleDouble value) => ToInteger<long>(value, ConvertMode.Checked);
    public static explicit operator ulong(DoubleDouble value) => ToInteger<ulong>(value, ConvertMode.Checked);
    public static explicit operator nint(DoubleDouble value) => ToInteger<nint>(value, ConvertMode.Checked);
    public static explicit operator nuint(DoubleDouble value) => ToInteger<nuint>(value, ConvertMode.Checked);
    public static explicit operator Int128(DoubleDouble value) => ToInteger<Int128>(value, ConvertMode.Checked);
    public static explicit operator UInt128(DoubleDouble value) => ToInteger<UInt128>(value, ConvertMode.Checked);

    /// <summary>
    /// Truncating conversion to an integer primitive. Checked throws on NaN/±∞/overflow; saturating
    /// and truncating both clamp to the target range with NaN -> 0 (the BCL rule for double).
    /// </summary>
    private static T ToInteger<T>(DoubleDouble value, ConvertMode mode) where T : INumberBase<T>
    {
        double hi = value._hi;
        if (Math.Abs(hi) < 4.611686018427388E18) // 2^62: hi + lo cannot leave long
        {
            DoubleDouble t = Truncate(value);
            long l = (long)t._hi + (long)t._lo;
            return mode == ConvertMode.Checked ? T.CreateChecked(l) : T.CreateSaturating(l);
        }
        if (!double.IsFinite(hi))
        {
            if (mode == ConvertMode.Checked)
            {
                throw new OverflowException($"Value does not fit in {typeof(T).Name}.");
            }
            return T.CreateSaturating(hi);
        }
        if (Math.Abs(hi) > 1e40)
        {
            if (mode == ConvertMode.Checked)
            {
                throw new OverflowException($"Value does not fit in {typeof(T).Name}.");
            }
            return T.CreateSaturating(hi);
        }
        return BigInt.ConvertToIntegerPrimitive<T>(value.ToBigIntTruncated(), mode == ConvertMode.Checked ? ConvertMode.Checked : ConvertMode.Saturating);
    }

    // ---------------------------------------------------------------
    // Mixed 64-bit integer operators: without these, C# would pick the double overloads and round
    // the integer to 53 bits first.
    // ---------------------------------------------------------------

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static bool FitsDouble(long v) => (ulong)(v + (1L << 53)) <= (1UL << 54);

    public static DoubleDouble operator +(DoubleDouble a, long b) => FitsDouble(b) ? a + (double)b : a + (DoubleDouble)b;
    public static DoubleDouble operator +(long a, DoubleDouble b) => b + a;
    public static DoubleDouble operator -(DoubleDouble a, long b) => FitsDouble(b) ? a - (double)b : a - (DoubleDouble)b;
    public static DoubleDouble operator -(long a, DoubleDouble b) => FitsDouble(a) ? (double)a - b : (DoubleDouble)a - b;
    public static DoubleDouble operator *(DoubleDouble a, long b) => FitsDouble(b) ? a * (double)b : a * (DoubleDouble)b;
    public static DoubleDouble operator *(long a, DoubleDouble b) => b * a;
    public static DoubleDouble operator /(DoubleDouble a, long b) => FitsDouble(b) ? a / (double)b : a / (DoubleDouble)b;
    public static DoubleDouble operator /(long a, DoubleDouble b) => (DoubleDouble)a / b;

    public static DoubleDouble operator +(DoubleDouble a, ulong b) => a + (DoubleDouble)b;
    public static DoubleDouble operator +(ulong a, DoubleDouble b) => (DoubleDouble)a + b;
    public static DoubleDouble operator -(DoubleDouble a, ulong b) => a - (DoubleDouble)b;
    public static DoubleDouble operator -(ulong a, DoubleDouble b) => (DoubleDouble)a - b;
    public static DoubleDouble operator *(DoubleDouble a, ulong b) => a * (DoubleDouble)b;
    public static DoubleDouble operator *(ulong a, DoubleDouble b) => (DoubleDouble)a * b;
    public static DoubleDouble operator /(DoubleDouble a, ulong b) => a / (DoubleDouble)b;
    public static DoubleDouble operator /(ulong a, DoubleDouble b) => (DoubleDouble)a / b;

    // ---------------------------------------------------------------
    // Generic-math conversions
    // ---------------------------------------------------------------

    static bool INumberBase<DoubleDouble>.TryConvertFromChecked<TOther>(TOther value, out DoubleDouble result) => TryConvertFrom(value, out result);
    static bool INumberBase<DoubleDouble>.TryConvertFromSaturating<TOther>(TOther value, out DoubleDouble result) => TryConvertFrom(value, out result);
    static bool INumberBase<DoubleDouble>.TryConvertFromTruncating<TOther>(TOther value, out DoubleDouble result) => TryConvertFrom(value, out result);

    static bool INumberBase<DoubleDouble>.TryConvertToChecked<TOther>(DoubleDouble value, [MaybeNullWhen(false)] out TOther result) =>
        TryConvertTo(value, ConvertMode.Checked, out result);

    static bool INumberBase<DoubleDouble>.TryConvertToSaturating<TOther>(DoubleDouble value, [MaybeNullWhen(false)] out TOther result) =>
        TryConvertTo(value, ConvertMode.Saturating, out result);

    static bool INumberBase<DoubleDouble>.TryConvertToTruncating<TOther>(DoubleDouble value, [MaybeNullWhen(false)] out TOther result) =>
        TryConvertTo(value, ConvertMode.Truncating, out result);

    private static bool TryConvertFrom<TOther>(TOther value, out DoubleDouble result) where TOther : INumberBase<TOther>
    {
        if (BigInt.IsSignedIntegerPrimitive<TOther>())
        {
            result = long.CreateTruncating(value);
            return true;
        }
        if (BigInt.IsUnsignedIntegerPrimitive<TOther>())
        {
            result = ulong.CreateTruncating(value);
            return true;
        }
        if (typeof(TOther) == typeof(Int128))
        {
            result = (DoubleDouble)(Int128)(object)value;
            return true;
        }
        if (typeof(TOther) == typeof(UInt128))
        {
            result = (DoubleDouble)(UInt128)(object)value;
            return true;
        }
        if (BigInt.IsBinaryFloatPrimitive<TOther>())
        {
            result = double.CreateTruncating(value); // exact for float and Half
            return true;
        }
        if (typeof(TOther) == typeof(decimal))
        {
            result = (DoubleDouble)(decimal)(object)value;
            return true;
        }
        if (typeof(TOther) == typeof(BigInteger))
        {
            result = (DoubleDouble)(BigInteger)(object)value;
            return true;
        }
        if (typeof(TOther) == typeof(BigInt))
        {
            result = (DoubleDouble)(BigInt)(object)value;
            return true;
        }
        if (typeof(TOther) == typeof(BigRational))
        {
            result = (DoubleDouble)(BigRational)(object)value;
            return true;
        }
        if (typeof(TOther) == typeof(BigFloat))
        {
            result = FromBigFloat((BigFloat)(object)value);
            return true;
        }
        result = default;
        return false;
    }

    private static bool TryConvertTo<TOther>(DoubleDouble value, ConvertMode mode, [MaybeNullWhen(false)] out TOther result)
        where TOther : INumberBase<TOther>
    {
        if (BigInt.IsIntegerPrimitive<TOther>())
        {
            result = ToInteger<TOther>(value, mode);
            return true;
        }
        if (typeof(TOther) == typeof(double))
        {
            result = (TOther)(object)value._hi;
            return true;
        }
        if (typeof(TOther) == typeof(float))
        {
            result = (TOther)(object)(float)value;
            return true;
        }
        if (typeof(TOther) == typeof(System.Half))
        {
            result = (TOther)(object)(System.Half)value;
            return true;
        }
        if (typeof(TOther) == typeof(decimal))
        {
            decimal m;
            if (!double.IsFinite(value._hi) || Math.Abs(value._hi) >= 7.9228162514264338E28)
            {
                if (mode == ConvertMode.Checked)
                {
                    throw new OverflowException("Value does not fit in Decimal.");
                }
                m = double.IsNaN(value._hi) ? 0m : value._hi < 0 ? decimal.MinValue : decimal.MaxValue;
            }
            else
            {
                m = (decimal)value;
            }
            result = (TOther)(object)m;
            return true;
        }
        if (typeof(TOther) == typeof(BigFloat))
        {
            result = (TOther)(object)(BigFloat)value;
            return true;
        }
        if (typeof(TOther) == typeof(BigInteger) || typeof(TOther) == typeof(BigInt) || typeof(TOther) == typeof(BigRational))
        {
            if (double.IsNaN(value._hi) && mode != ConvertMode.Checked)
            {
                result = TOther.Zero;
                return true;
            }
            if (typeof(TOther) == typeof(BigRational))
            {
                result = (TOther)(object)(BigRational)value;
            }
            else
            {
                BigInt b = value.ToBigIntTruncated();
                result = typeof(TOther) == typeof(BigInt) ? (TOther)(object)b : (TOther)(object)(BigInteger)b;
            }
            return true;
        }
        result = default;
        return false;
    }

    public static DoubleDouble CreateChecked<TOther>(TOther value) where TOther : INumberBase<TOther>
    {
        if (typeof(TOther) == typeof(DoubleDouble))
        {
            return (DoubleDouble)(object)value;
        }
        if (!TryConvertFrom(value, out DoubleDouble result) && !TOther.TryConvertToChecked(value, out result))
        {
            throw new NotSupportedException($"Conversion from {typeof(TOther).Name} to DoubleDouble is not supported.");
        }
        return result;
    }

    public static DoubleDouble CreateSaturating<TOther>(TOther value) where TOther : INumberBase<TOther>
    {
        if (typeof(TOther) == typeof(DoubleDouble))
        {
            return (DoubleDouble)(object)value;
        }
        if (!TryConvertFrom(value, out DoubleDouble result) && !TOther.TryConvertToSaturating(value, out result))
        {
            throw new NotSupportedException($"Conversion from {typeof(TOther).Name} to DoubleDouble is not supported.");
        }
        return result;
    }

    public static DoubleDouble CreateTruncating<TOther>(TOther value) where TOther : INumberBase<TOther>
    {
        if (typeof(TOther) == typeof(DoubleDouble))
        {
            return (DoubleDouble)(object)value;
        }
        if (!TryConvertFrom(value, out DoubleDouble result) && !TOther.TryConvertToTruncating(value, out result))
        {
            throw new NotSupportedException($"Conversion from {typeof(TOther).Name} to DoubleDouble is not supported.");
        }
        return result;
    }
}
