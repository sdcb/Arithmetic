using System.Diagnostics.CodeAnalysis;
using System.Numerics;

namespace Sdcb.Arithmetic;

public readonly partial struct Float128
{
    // ---------------------------------------------------------------
    // From binary floating point (exact)
    // ---------------------------------------------------------------

    /// <summary>Converts a double exactly (NaN becomes the canonical NaN).</summary>
    public static Float128 FromDouble(double value)
    {
        ulong b = BitConverter.DoubleToUInt64Bits(value);
        bool sign = (long)b < 0;
        int exp = (int)(b >> 52) & 0x7FF;
        ulong frac = b & 0x000F_FFFF_FFFF_FFFFUL;
        if (exp == 0x7FF)
        {
            return frac != 0 ? NaN : InfOf(sign);
        }
        if (exp == 0)
        {
            if (frac == 0)
            {
                return ZeroOf(sign);
            }
            int m = 63 - BitOperations.LeadingZeroCount(frac); // value = frac * 2^-1074
            UInt128 f = ((UInt128)frac << (112 - m)) & FracMask;
            return new(((UInt128)(uint)(m - 1074 + ExponentBias) << 112) | f | SignBits(sign));
        }
        return new(((UInt128)(uint)(exp - 1023 + ExponentBias) << 112) | ((UInt128)frac << 60) | SignBits(sign));
    }

    private static Float128 FromUInt64(bool negative, ulong mag)
    {
        if (mag == 0)
        {
            return Zero;
        }
        int m = 63 - BitOperations.LeadingZeroCount(mag);
        UInt128 f = ((UInt128)mag << (112 - m)) & FracMask;
        return new(((UInt128)(uint)(m + ExponentBias) << 112) | f | SignBits(negative));
    }

    private static Float128 FromUInt128(bool negative, UInt128 mag, RoundingMode rounding)
    {
        if (mag == 0)
        {
            return Zero;
        }
        return RoundWide(negative, mag, 0, sticky: false, rounding);
    }

    /// <summary>Converts an <see cref="Int128"/>, correctly rounded.</summary>
    public static Float128 FromInt128(Int128 value, RoundingMode rounding = RoundingMode.NearestEven) =>
        Int128.IsNegative(value)
            ? FromUInt128(true, (UInt128)(-(value + 1)) + 1, rounding)
            : FromUInt128(false, (UInt128)value, rounding);

    /// <summary>Converts a <see cref="UInt128"/>, correctly rounded.</summary>
    public static Float128 FromUInt128(UInt128 value, RoundingMode rounding = RoundingMode.NearestEven) =>
        FromUInt128(false, value, rounding);

    /// <summary>Converts a <see cref="BigInt"/>, correctly rounded (±∞ or ±MaxValue on overflow per mode).</summary>
    public static Float128 FromBigInt(BigInt value, RoundingMode rounding = RoundingMode.NearestEven)
    {
        if (value.IsZero)
        {
            return Zero;
        }
        bool neg = value.Sign < 0;
        BigInt mag = value.Abs();
        if (mag.GetBitLength() > 16400)
        {
            return Huge(neg, rounding);
        }
        return RoundBig(neg, mag, 0, sticky: false, rounding);
    }

    /// <summary>Converts a <see cref="BigRational"/>, correctly rounded.</summary>
    public static Float128 FromRational(BigRational value, RoundingMode rounding = RoundingMode.NearestEven)
    {
        BigInt num = value.Numerator;
        if (num.IsZero)
        {
            return Zero;
        }
        return RoundRational(num.Sign < 0, num.Abs(), value.Denominator, 0, rounding);
    }

    /// <summary>Converts a decimal, correctly rounded (a negative zero stays negative).</summary>
    public static Float128 FromDecimal(decimal value, RoundingMode rounding = RoundingMode.NearestEven)
    {
        if (value == 0m)
        {
            return ZeroOf(decimal.IsNegative(value));
        }
        return FromRational(BigFloat.DecimalToRational(value), rounding);
    }

    /// <summary>Converts a <see cref="BigFloat"/>, correctly rounded to binary128 (subnormals and overflow included).</summary>
    public static Float128 FromBigFloat(BigFloat value, RoundingMode rounding = RoundingMode.NearestEven) =>
        new(value.ToBinary128Bits(rounding));

    /// <summary>The exact value as a <see cref="BigFloat"/> with 113 bits of precision (NaN and ±∞ map to their BigFloat counterparts).</summary>
    public BigFloat ToBigFloat() => BigFloat.FromBinary128Bits(_bits, SignificandBits);

    // ---------------------------------------------------------------
    // Implicit (exact) conversions
    // ---------------------------------------------------------------

    public static implicit operator Float128(double value) => FromDouble(value);
    public static implicit operator Float128(float value) => FromDouble(value);
    public static implicit operator Float128(Half value) => FromDouble((double)value);
    public static implicit operator Float128(byte value) => FromUInt64(false, value);
    public static implicit operator Float128(sbyte value) => FromInt64(value);
    public static implicit operator Float128(short value) => FromInt64(value);
    public static implicit operator Float128(ushort value) => FromUInt64(false, value);
    public static implicit operator Float128(char value) => FromUInt64(false, value);
    public static implicit operator Float128(int value) => FromInt64(value);
    public static implicit operator Float128(uint value) => FromUInt64(false, value);
    public static implicit operator Float128(long value) => FromInt64(value);
    public static implicit operator Float128(ulong value) => FromUInt64(false, value);
    public static implicit operator Float128(nint value) => FromInt64(value);
    public static implicit operator Float128(nuint value) => FromUInt64(false, value);

    private static Float128 FromInt64(long value) =>
        value < 0 ? FromUInt64(true, (ulong)(-(value + 1)) + 1) : FromUInt64(false, (ulong)value);

    // ---------------------------------------------------------------
    // Explicit (rounding) conversions into Float128
    // ---------------------------------------------------------------

    /// <summary>Correctly rounded (nearest-even).</summary>
    public static explicit operator Float128(Int128 value) => FromInt128(value);

    /// <summary>Correctly rounded (nearest-even).</summary>
    public static explicit operator Float128(UInt128 value) => FromUInt128(value);

    /// <summary>Correctly rounded (nearest-even).</summary>
    public static explicit operator Float128(decimal value) => FromDecimal(value);

    /// <summary>Correctly rounded (nearest-even); ±∞ on overflow.</summary>
    public static explicit operator Float128(BigInt value) => FromBigInt(value);

    /// <summary>Correctly rounded (nearest-even); ±∞ on overflow.</summary>
    public static explicit operator Float128(BigRational value) => FromRational(value);

    /// <summary>Correctly rounded (nearest-even), including subnormals and overflow.</summary>
    public static explicit operator Float128(BigFloat value) => FromBigFloat(value);

    // ---------------------------------------------------------------
    // Conversions out of Float128
    // ---------------------------------------------------------------

    /// <summary>Exact (113-bit <see cref="BigFloat"/>).</summary>
    public static implicit operator BigFloat(Float128 value) => value.ToBigFloat();

    /// <summary>Correctly rounded (nearest-even), including double subnormals and overflow to ±∞.</summary>
    public static explicit operator double(Float128 value) => value.ToDouble();

    /// <summary>Correctly rounded (nearest-even), single rounding.</summary>
    public static explicit operator float(Float128 value) => value.ToSingle();

    /// <summary>Correctly rounded (nearest-even), single rounding.</summary>
    public static explicit operator Half(Float128 value) => value.ToHalf();

    /// <summary>Converts to double with the given rounding.</summary>
    public double ToDouble(RoundingMode rounding = RoundingMode.NearestEven) =>
        BitConverter.UInt64BitsToDouble(Narrow(_bits, 52, 11, rounding));

    /// <summary>Converts to float with the given rounding (one rounding step, no double rounding).</summary>
    public float ToSingle(RoundingMode rounding = RoundingMode.NearestEven) =>
        BitConverter.UInt32BitsToSingle((uint)Narrow(_bits, 23, 8, rounding));

    /// <summary>Converts to <see cref="Half"/> with the given rounding (one rounding step, no double rounding).</summary>
    public Half ToHalf(RoundingMode rounding = RoundingMode.NearestEven) =>
        BitConverter.UInt16BitsToHalf((ushort)Narrow(_bits, 10, 5, rounding));

    /// <summary>Rounds a binary128 value to a narrower IEEE binary format and returns its bit pattern.</summary>
    private static ulong Narrow(UInt128 bits, int fracBits, int expBits, RoundingMode rounding)
    {
        bool sign = SignOf(bits);
        int exp = ExpOf(bits);
        UInt128 sig = FracOf(bits);
        ulong signOut = sign ? 1UL << (fracBits + expBits) : 0;
        int maxExpOut = (1 << expBits) - 1;
        ulong infOut = (ulong)maxExpOut << fracBits;
        if (exp == MaxBiasedExponent)
        {
            return sig != 0 ? infOut | (1UL << (fracBits - 1)) : signOut | infOut;
        }
        if (exp == 0)
        {
            if (sig == 0)
            {
                return signOut;
            }
            sig = NormSubnormal(sig, out exp);
        }
        sig |= Hidden;
        int biasOut = (1 << (expBits - 1)) - 1;
        long e = (long)exp - ExponentBias + biasOut; // target biased exponent of the leading bit
        if (e >= maxExpOut)
        {
            return signOut | NarrowHuge(sign, rounding, infOut);
        }
        int shift = FractionBits - fracBits;
        if (e <= 0)
        {
            shift += (int)Math.Min(1 - e, 200);
            e = 1;
        }
        UInt128 q;
        bool half, above, any;
        if (shift >= 114)
        {
            q = UInt128.Zero;
            half = false;
            above = false;
            any = true;
        }
        else
        {
            q = sig >> shift;
            UInt128 rem = sig & ((UInt128.One << shift) - 1);
            UInt128 h = UInt128.One << (shift - 1);
            half = rem == h;
            above = rem > h;
            any = rem != 0;
        }
        bool inc = rounding switch
        {
            RoundingMode.NearestEven => above || (half && (q & 1) != 0),
            RoundingMode.TowardZero => false,
            RoundingMode.TowardPositive => any && !sign,
            RoundingMode.TowardNegative => any && sign,
            RoundingMode.AwayFromZero => any,
            _ => throw new ArgumentOutOfRangeException(nameof(rounding)),
        };
        if (inc)
        {
            q++;
        }
        ulong result = ((ulong)(e - 1) << fracBits) + (ulong)q;
        if (result >= infOut)
        {
            return signOut | NarrowHuge(sign, rounding, infOut);
        }
        return signOut | result;
    }

    private static ulong NarrowHuge(bool negative, RoundingMode rounding, ulong infOut)
    {
        bool toInf = rounding switch
        {
            RoundingMode.TowardZero => false,
            RoundingMode.TowardPositive => !negative,
            RoundingMode.TowardNegative => negative,
            _ => true,
        };
        return toInf ? infOut : infOut - 1;
    }

    // ---------------------------------------------------------------
    // To integers (truncation toward zero; OverflowException like BigFloat for NaN, ±∞ and out-of-range values)
    // ---------------------------------------------------------------

    /// <summary>|trunc(x)| when |x| &lt; 2^128 and x is finite; false otherwise.</summary>
    private bool TryTruncateMagnitude(out UInt128 magnitude)
    {
        int exp = BiasedExponent;
        magnitude = UInt128.Zero;
        if (exp == MaxBiasedExponent)
        {
            return false;
        }
        int e = exp - ExponentBias;
        if (e < 0)
        {
            return true;
        }
        if (e >= 128)
        {
            return false;
        }
        UInt128 sig = Fraction | Hidden;
        magnitude = e >= 112 ? sig << (e - 112) : sig >> (112 - e);
        return true;
    }

    private enum IntMode
    {
        Checked,
        Saturating,
    }

    /// <summary>Truncates to an integer primitive: checked (throws) or saturating (NaN → 0).</summary>
    private static TInt ToInteger<TInt>(Float128 value, IntMode mode) where TInt : IBinaryInteger<TInt>, IMinMaxValue<TInt>
    {
        if (IsNaN(value))
        {
            return mode == IntMode.Checked ? throw new OverflowException("Cannot convert NaN to an integer.") : TInt.Zero;
        }
        bool neg = value.SignBit;
        if (!value.TryTruncateMagnitude(out UInt128 mag))
        {
            return mode == IntMode.Checked
                ? throw new OverflowException($"Value is outside the range of {typeof(TInt).Name}.")
                : neg ? TInt.MinValue : TInt.MaxValue;
        }
        if (!neg || mag == 0)
        {
            return mode == IntMode.Checked ? TInt.CreateChecked(mag) : TInt.CreateSaturating(mag);
        }
        if (mag > (UInt128)Int128.MaxValue + 1)
        {
            return mode == IntMode.Checked
                ? throw new OverflowException($"Value is outside the range of {typeof(TInt).Name}.")
                : TInt.MinValue;
        }
        Int128 s = unchecked(-(Int128)mag);
        return mode == IntMode.Checked ? TInt.CreateChecked(s) : TInt.CreateSaturating(s);
    }

    public static explicit operator sbyte(Float128 value) => ToInteger<sbyte>(value, IntMode.Checked);
    public static explicit operator byte(Float128 value) => ToInteger<byte>(value, IntMode.Checked);
    public static explicit operator short(Float128 value) => ToInteger<short>(value, IntMode.Checked);
    public static explicit operator ushort(Float128 value) => ToInteger<ushort>(value, IntMode.Checked);
    public static explicit operator char(Float128 value) => ToInteger<char>(value, IntMode.Checked);
    public static explicit operator int(Float128 value) => ToInteger<int>(value, IntMode.Checked);
    public static explicit operator uint(Float128 value) => ToInteger<uint>(value, IntMode.Checked);
    public static explicit operator long(Float128 value) => ToInteger<long>(value, IntMode.Checked);
    public static explicit operator ulong(Float128 value) => ToInteger<ulong>(value, IntMode.Checked);
    public static explicit operator nint(Float128 value) => ToInteger<nint>(value, IntMode.Checked);
    public static explicit operator nuint(Float128 value) => ToInteger<nuint>(value, IntMode.Checked);
    public static explicit operator Int128(Float128 value) => ToInteger<Int128>(value, IntMode.Checked);
    public static explicit operator UInt128(Float128 value) => ToInteger<UInt128>(value, IntMode.Checked);

    /// <summary>Truncates toward zero; throws <see cref="OverflowException"/> for NaN or infinity.</summary>
    public static explicit operator BigInt(Float128 value)
    {
        if (!IsFinite(value))
        {
            throw new OverflowException("Cannot convert NaN or Infinity to BigInt.");
        }
        if (value.TryTruncateMagnitude(out UInt128 mag))
        {
            BigInt m = new(mag);
            return value.SignBit ? -m : m;
        }
        return value.ToBigFloat().ToBigInt();
    }

    /// <summary>The exact value; throws <see cref="OverflowException"/> for NaN or infinity.</summary>
    public static explicit operator BigRational(Float128 value) => value.ToBigFloat().ToRational();

    /// <summary>Nearest decimal (ties to even); throws <see cref="OverflowException"/> for NaN, infinity or out-of-range values.</summary>
    public static explicit operator decimal(Float128 value) => value.ToBigFloat().ToDecimal();

    // ---------------------------------------------------------------
    // Generic-math conversions
    // ---------------------------------------------------------------

    static bool INumberBase<Float128>.TryConvertFromChecked<TOther>(TOther value, out Float128 result) => TryConvertFrom(value, out result);
    static bool INumberBase<Float128>.TryConvertFromSaturating<TOther>(TOther value, out Float128 result) => TryConvertFrom(value, out result);
    static bool INumberBase<Float128>.TryConvertFromTruncating<TOther>(TOther value, out Float128 result) => TryConvertFrom(value, out result);

    static bool INumberBase<Float128>.TryConvertToChecked<TOther>(Float128 value, [MaybeNullWhen(false)] out TOther result) =>
        TryConvertTo(value, IntMode.Checked, out result);

    static bool INumberBase<Float128>.TryConvertToSaturating<TOther>(Float128 value, [MaybeNullWhen(false)] out TOther result) =>
        TryConvertTo(value, IntMode.Saturating, out result);

    /// <summary>Like <see cref="double"/>, truncating conversions to integers saturate (NaN → 0).</summary>
    static bool INumberBase<Float128>.TryConvertToTruncating<TOther>(Float128 value, [MaybeNullWhen(false)] out TOther result) =>
        TryConvertTo(value, IntMode.Saturating, out result);

    private static bool TryConvertFrom<TOther>(TOther value, out Float128 result) where TOther : INumberBase<TOther>
    {
        if (typeof(TOther) == typeof(double))
        {
            result = (double)(object)value;
        }
        else if (typeof(TOther) == typeof(float))
        {
            result = (float)(object)value;
        }
        else if (typeof(TOther) == typeof(Half))
        {
            result = (Half)(object)value;
        }
        else if (typeof(TOther) == typeof(int))
        {
            result = (int)(object)value;
        }
        else if (typeof(TOther) == typeof(long))
        {
            result = (long)(object)value;
        }
        else if (typeof(TOther) == typeof(uint))
        {
            result = (uint)(object)value;
        }
        else if (typeof(TOther) == typeof(ulong))
        {
            result = (ulong)(object)value;
        }
        else if (typeof(TOther) == typeof(short))
        {
            result = (short)(object)value;
        }
        else if (typeof(TOther) == typeof(ushort))
        {
            result = (ushort)(object)value;
        }
        else if (typeof(TOther) == typeof(sbyte))
        {
            result = (sbyte)(object)value;
        }
        else if (typeof(TOther) == typeof(byte))
        {
            result = (byte)(object)value;
        }
        else if (typeof(TOther) == typeof(char))
        {
            result = (char)(object)value;
        }
        else if (typeof(TOther) == typeof(nint))
        {
            result = (nint)(object)value;
        }
        else if (typeof(TOther) == typeof(nuint))
        {
            result = (nuint)(object)value;
        }
        else if (typeof(TOther) == typeof(Int128))
        {
            result = FromInt128((Int128)(object)value);
        }
        else if (typeof(TOther) == typeof(UInt128))
        {
            result = FromUInt128((UInt128)(object)value);
        }
        else if (typeof(TOther) == typeof(decimal))
        {
            result = FromDecimal((decimal)(object)value);
        }
        else if (typeof(TOther) == typeof(BigInt))
        {
            result = FromBigInt((BigInt)(object)value);
        }
        else if (typeof(TOther) == typeof(BigInteger))
        {
            result = FromBigInt((BigInteger)(object)value);
        }
        else if (typeof(TOther) == typeof(BigRational))
        {
            result = FromRational((BigRational)(object)value);
        }
        else if (typeof(TOther) == typeof(BigFloat))
        {
            result = FromBigFloat((BigFloat)(object)value);
        }
        else
        {
            result = default;
            return false;
        }
        return true;
    }

    private static bool TryConvertTo<TOther>(Float128 value, IntMode mode, [MaybeNullWhen(false)] out TOther result)
        where TOther : INumberBase<TOther>
    {
        if (typeof(TOther) == typeof(double))
        {
            result = (TOther)(object)value.ToDouble();
        }
        else if (typeof(TOther) == typeof(float))
        {
            result = (TOther)(object)value.ToSingle();
        }
        else if (typeof(TOther) == typeof(Half))
        {
            result = (TOther)(object)value.ToHalf();
        }
        else if (typeof(TOther) == typeof(int))
        {
            result = (TOther)(object)ToInteger<int>(value, mode);
        }
        else if (typeof(TOther) == typeof(long))
        {
            result = (TOther)(object)ToInteger<long>(value, mode);
        }
        else if (typeof(TOther) == typeof(uint))
        {
            result = (TOther)(object)ToInteger<uint>(value, mode);
        }
        else if (typeof(TOther) == typeof(ulong))
        {
            result = (TOther)(object)ToInteger<ulong>(value, mode);
        }
        else if (typeof(TOther) == typeof(short))
        {
            result = (TOther)(object)ToInteger<short>(value, mode);
        }
        else if (typeof(TOther) == typeof(ushort))
        {
            result = (TOther)(object)ToInteger<ushort>(value, mode);
        }
        else if (typeof(TOther) == typeof(sbyte))
        {
            result = (TOther)(object)ToInteger<sbyte>(value, mode);
        }
        else if (typeof(TOther) == typeof(byte))
        {
            result = (TOther)(object)ToInteger<byte>(value, mode);
        }
        else if (typeof(TOther) == typeof(char))
        {
            result = (TOther)(object)ToInteger<char>(value, mode);
        }
        else if (typeof(TOther) == typeof(nint))
        {
            result = (TOther)(object)ToInteger<nint>(value, mode);
        }
        else if (typeof(TOther) == typeof(nuint))
        {
            result = (TOther)(object)ToInteger<nuint>(value, mode);
        }
        else if (typeof(TOther) == typeof(Int128))
        {
            result = (TOther)(object)ToInteger<Int128>(value, mode);
        }
        else if (typeof(TOther) == typeof(UInt128))
        {
            result = (TOther)(object)ToInteger<UInt128>(value, mode);
        }
        else if (typeof(TOther) == typeof(decimal))
        {
            if (!TryToDecimal(value, out decimal m))
            {
                if (mode == IntMode.Checked)
                {
                    throw new OverflowException("Value does not fit in Decimal.");
                }
                m = IsNaN(value) ? 0m : value.SignBit ? decimal.MinValue : decimal.MaxValue;
            }
            result = (TOther)(object)m;
        }
        else if (typeof(TOther) == typeof(BigInt) || typeof(TOther) == typeof(BigInteger) || typeof(TOther) == typeof(BigRational))
        {
            if (IsNaN(value) && mode != IntMode.Checked)
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
                BigInt b = (BigInt)value;
                result = typeof(TOther) == typeof(BigInt) ? (TOther)(object)b : (TOther)(object)(BigInteger)b;
            }
        }
        else if (typeof(TOther) == typeof(BigFloat))
        {
            result = (TOther)(object)value.ToBigFloat();
        }
        else
        {
            result = default;
            return false;
        }
        return true;
    }

    private static bool TryToDecimal(Float128 value, out decimal result)
    {
        if (!IsFinite(value) || ILogB(value) >= 96)
        {
            result = 0m;
            return false;
        }
        try
        {
            // Rare: values just below 2^96 that round up to 2^96 overflow inside BigFloat.ToDecimal.
            result = value.ToBigFloat().ToDecimal();
            return true;
        }
        catch (OverflowException)
        {
            result = 0m;
            return false;
        }
    }

    // ---------------------------------------------------------------
    // Public CreateChecked / CreateSaturating / CreateTruncating
    // ---------------------------------------------------------------

    public static Float128 CreateChecked<TOther>(TOther value) where TOther : INumberBase<TOther>
    {
        if (typeof(TOther) == typeof(Float128))
        {
            return (Float128)(object)value;
        }
        if (!TryConvertFrom(value, out Float128 result) && !TOther.TryConvertToChecked(value, out result))
        {
            throw new NotSupportedException($"Conversion from {typeof(TOther).Name} to Float128 is not supported.");
        }
        return result;
    }

    public static Float128 CreateSaturating<TOther>(TOther value) where TOther : INumberBase<TOther>
    {
        if (typeof(TOther) == typeof(Float128))
        {
            return (Float128)(object)value;
        }
        if (!TryConvertFrom(value, out Float128 result) && !TOther.TryConvertToSaturating(value, out result))
        {
            throw new NotSupportedException($"Conversion from {typeof(TOther).Name} to Float128 is not supported.");
        }
        return result;
    }

    public static Float128 CreateTruncating<TOther>(TOther value) where TOther : INumberBase<TOther>
    {
        if (typeof(TOther) == typeof(Float128))
        {
            return (Float128)(object)value;
        }
        if (!TryConvertFrom(value, out Float128 result) && !TOther.TryConvertToTruncating(value, out result))
        {
            throw new NotSupportedException($"Conversion from {typeof(TOther).Name} to Float128 is not supported.");
        }
        return result;
    }
}
