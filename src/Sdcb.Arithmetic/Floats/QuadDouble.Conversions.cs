using System.Diagnostics.CodeAnalysis;
using System.Numerics;
using System.Runtime.CompilerServices;

namespace Sdcb.Arithmetic;

public readonly partial struct QuadDouble
{
    // ---------------------------------------------------------------
    // Implicit conversions (exact)
    // ---------------------------------------------------------------

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static implicit operator QuadDouble(double value) => Single(value);

    public static implicit operator QuadDouble(float value) => Single(value);
    public static implicit operator QuadDouble(System.Half value) => Single((double)value);
    public static implicit operator QuadDouble(sbyte value) => Single(value);
    public static implicit operator QuadDouble(byte value) => Single(value);
    public static implicit operator QuadDouble(short value) => Single(value);
    public static implicit operator QuadDouble(ushort value) => Single(value);
    public static implicit operator QuadDouble(char value) => Single(value);
    public static implicit operator QuadDouble(int value) => Single(value);
    public static implicit operator QuadDouble(uint value) => Single(value);

    /// <summary>Exact (a long needs at most two components).</summary>
    public static implicit operator QuadDouble(long value)
    {
        double hi = value;
        if ((long)hi == value && hi != 9.2233720368547758E+18)
        {
            return Single(hi);
        }
        return FromInt128(value);
    }

    /// <summary>Exact (a ulong needs at most two components).</summary>
    public static implicit operator QuadDouble(ulong value) => value <= (1UL << 53) ? Single(value) : FromUInt128(value, negative: false);

    public static implicit operator QuadDouble(nint value) => (long)value;
    public static implicit operator QuadDouble(nuint value) => (ulong)value;

    /// <summary>Exact (128 bits fit in three components).</summary>
    public static implicit operator QuadDouble(Int128 value) => FromInt128(value);

    /// <summary>Exact (128 bits fit in three components).</summary>
    public static implicit operator QuadDouble(UInt128 value) => FromUInt128(value, negative: false);

    private static QuadDouble FromInt128(Int128 value)
    {
        bool neg = value < 0;
        UInt128 mag = neg ? (UInt128)(-(value + 1)) + 1 : (UInt128)value;
        return FromUInt128(mag, neg);
    }

    /// <summary>Splits the magnitude into 53-bit chunks (each an exact double) and sums them in QD.</summary>
    private static QuadDouble FromUInt128(UInt128 mag, bool negative)
    {
        const ulong mask = (1UL << 53) - 1;
        double c = (double)(ulong)(mag & mask);
        double b = Math.ScaleB((double)(ulong)((mag >> 53) & mask), 53);
        double a = Math.ScaleB((double)(ulong)(mag >> 106), 106);
        QuadDouble r = Add(Add(Single(a), b), c);
        return negative ? -r : r;
    }

    // ---------------------------------------------------------------
    // Explicit conversions into QuadDouble
    // ---------------------------------------------------------------

    /// <summary>The decimal rounded to QD (via a 320-bit intermediate).</summary>
    public static explicit operator QuadDouble(decimal value)
    {
        QuadDouble r = FromBigFloat(BigFloat.FromDecimal(value, 320));
        return r._x0 == 0 && decimal.IsNegative(value) ? NegativeZero : r;
    }

    /// <summary>The integer rounded to QD; ±∞ beyond the double range.</summary>
    public static explicit operator QuadDouble(BigInt value)
    {
        long bits = value.GetBitLength();
        if (bits > 1100)
        {
            return value.Sign < 0 ? NegativeInfinity : PositiveInfinity;
        }
        return FromBigFloat(BigFloat.FromBigInt(value, (int)Math.Max(bits, 64)));
    }

    /// <summary>The integer rounded to QD; ±∞ beyond the double range.</summary>
    public static explicit operator QuadDouble(BigInteger value) => (QuadDouble)(BigInt)value;

    /// <summary>The rational rounded to QD (via a 320-bit intermediate).</summary>
    public static explicit operator QuadDouble(BigRational value) => FromBigFloat(BigFloat.FromRational(value, 320));

    /// <summary>
    /// The BigFloat rounded to QD by successive nearest-double extraction of the exact residuals
    /// (x0 = round(v), x1 = round(v - x0), ...).
    /// </summary>
    public static explicit operator QuadDouble(BigFloat value) => FromBigFloat(value);

    /// <summary>
    /// Rounds a <see cref="BigFloat"/> to quad-double by successive nearest-double extraction of the
    /// exact residuals (same as the explicit conversion).
    /// </summary>
    public static QuadDouble FromBigFloat(BigFloat value)
    {
        if (value.IsNaN)
        {
            return NaN;
        }
        if (value.IsInfinity)
        {
            return value.SignBit ? NegativeInfinity : PositiveInfinity;
        }
        if (value.IsZero)
        {
            return value.SignBit ? NegativeZero : Zero;
        }
        int prec = Math.Max(value.Precision, 64);
        double c0 = value.ToDouble();
        if (!double.IsFinite(c0))
        {
            return Single(c0);
        }
        BigFloat r = BigFloat.Subtract(value, BigFloat.FromDouble(c0, 64), prec);
        double c1 = r.ToDouble();
        r = BigFloat.Subtract(r, BigFloat.FromDouble(c1, 64), prec);
        double c2 = r.ToDouble();
        r = BigFloat.Subtract(r, BigFloat.FromDouble(c2, 64), prec);
        double c3 = r.ToDouble();
        return new QuadDouble(c0, c1, c2, c3, Tag.Raw);
    }

    // ---------------------------------------------------------------
    // Conversions out of QuadDouble
    // ---------------------------------------------------------------

    /// <summary>The exact value (the precision is chosen from the components' exponent span).</summary>
    public static implicit operator BigFloat(QuadDouble value) => ToBigFloatExact(value);

    /// <summary>The exact value as a <see cref="BigFloat"/> (same as the implicit conversion).</summary>
    public BigFloat ToBigFloat() => ToBigFloatExact(this);

    /// <summary>The value correctly rounded to float (nearest-even by default).</summary>
    public float ToSingle(RoundingMode rounding = RoundingMode.NearestEven) =>
        rounding == RoundingMode.NearestEven ? (float)this : ToBigFloatExact(this).ToSingle(rounding);

    /// <summary>The value correctly rounded to <see cref="System.Half"/> (nearest-even by default).</summary>
    public System.Half ToHalf(RoundingMode rounding = RoundingMode.NearestEven) =>
        rounding == RoundingMode.NearestEven ? (System.Half)this : ToBigFloatExact(this).ToHalf(rounding);

    internal static BigFloat ToBigFloatExact(QuadDouble a)
    {
        double x0 = a._x0;
        if (double.IsNaN(x0))
        {
            return BigFloat.NaN;
        }
        if (double.IsInfinity(x0))
        {
            return x0 > 0 ? BigFloat.PositiveInfinity : BigFloat.NegativeInfinity;
        }
        if (x0 == 0)
        {
            return double.IsNegative(x0) ? BigFloat.NegativeZero.WithPrecision(64) : BigFloat.Zero.WithPrecision(64);
        }
        double last = a._x3 != 0 ? a._x3 : a._x2 != 0 ? a._x2 : a._x1 != 0 ? a._x1 : x0;
        int lsb = Math.Max(Math.ILogB(last) - 52, -1074);
        int msb = Math.ILogB(x0) + 1;
        int prec = Math.Max(msb - lsb + 2, 64);
        BigFloat s = BigFloat.FromDouble(x0, prec);
        if (a._x1 != 0)
        {
            s = BigFloat.Add(s, BigFloat.FromDouble(a._x1, 64), prec);
        }
        if (a._x2 != 0)
        {
            s = BigFloat.Add(s, BigFloat.FromDouble(a._x2, 64), prec);
        }
        if (a._x3 != 0)
        {
            s = BigFloat.Add(s, BigFloat.FromDouble(a._x3, 64), prec);
        }
        return s;
    }

    /// <summary>The nearest double (ties to even).</summary>
    public static explicit operator double(QuadDouble value) => ToDoubleNearest(value);

    /// <summary>
    /// The value correctly rounded to double. Nearest-even uses the component sum directly; the
    /// directed modes round the exact value.
    /// </summary>
    public double ToDouble(RoundingMode rounding = RoundingMode.NearestEven) =>
        rounding == RoundingMode.NearestEven ? ToDoubleNearest(this) : ToBigFloatExact(this).ToDouble(rounding);

    /// <summary>The nearest float (ties to even; no double rounding).</summary>
    public static explicit operator float(QuadDouble value)
    {
        double d = ToDoubleNearest(value);
        if (IsSingleMidpoint(d))
        {
            return ToBigFloatExact(value).ToSingle();
        }
        return (float)d;
    }

    /// <summary>The nearest Half (ties to even; no double rounding).</summary>
    public static explicit operator System.Half(QuadDouble value)
    {
        double d = ToDoubleNearest(value);
        System.Half h = (System.Half)d;
        double back = (double)h;
        if (back != d && double.IsFinite(back))
        {
            System.Half other = d > back ? System.Half.BitIncrement(h) : System.Half.BitDecrement(h);
            if (double.IsFinite((double)other) && (double)other - d == d - back)
            {
                return ToBigFloatExact(value).ToHalf();
            }
        }
        return h;
    }

    /// <summary>True when d lies exactly halfway between two adjacent floats (where rounding twice could be wrong).</summary>
    private static bool IsSingleMidpoint(double d)
    {
        float f = (float)d;
        double back = f;
        if (back == d || !double.IsFinite(back))
        {
            return false;
        }
        float other = d > back ? MathF.BitIncrement(f) : MathF.BitDecrement(f);
        return double.IsFinite(other) && (double)other - d == d - back;
    }

    // Integer conversions truncate toward zero and saturate like the .NET 9+ double conversions:
    // NaN gives 0, out-of-range values give the type's MinValue / MaxValue.

    public static explicit operator sbyte(QuadDouble value) => ToInteger<sbyte>(value, ConvertMode.Saturating);
    public static explicit operator byte(QuadDouble value) => ToInteger<byte>(value, ConvertMode.Saturating);
    public static explicit operator short(QuadDouble value) => ToInteger<short>(value, ConvertMode.Saturating);
    public static explicit operator ushort(QuadDouble value) => ToInteger<ushort>(value, ConvertMode.Saturating);
    public static explicit operator char(QuadDouble value) => ToInteger<char>(value, ConvertMode.Saturating);
    public static explicit operator int(QuadDouble value) => ToInteger<int>(value, ConvertMode.Saturating);
    public static explicit operator uint(QuadDouble value) => ToInteger<uint>(value, ConvertMode.Saturating);
    public static explicit operator long(QuadDouble value) => ToInteger<long>(value, ConvertMode.Saturating);
    public static explicit operator ulong(QuadDouble value) => ToInteger<ulong>(value, ConvertMode.Saturating);
    public static explicit operator nint(QuadDouble value) => ToInteger<nint>(value, ConvertMode.Saturating);
    public static explicit operator nuint(QuadDouble value) => ToInteger<nuint>(value, ConvertMode.Saturating);
    public static explicit operator Int128(QuadDouble value) => ToInteger<Int128>(value, ConvertMode.Saturating);
    public static explicit operator UInt128(QuadDouble value) => ToInteger<UInt128>(value, ConvertMode.Saturating);

    /// <summary>The nearest decimal; throws <see cref="OverflowException"/> for NaN, ∞ or out-of-range values (like double).</summary>
    public static explicit operator decimal(QuadDouble value)
    {
        if (!double.IsFinite(value._x0))
        {
            throw new OverflowException("Value was either too large or too small for a Decimal.");
        }
        return ToBigFloatExact(value).ToDecimal();
    }

    /// <summary>Truncates toward zero; throws <see cref="OverflowException"/> for NaN or ∞.</summary>
    public static explicit operator BigInt(QuadDouble value) => ToBigFloatExact(value).ToBigInt();

    /// <summary>Truncates toward zero; throws <see cref="OverflowException"/> for NaN or ∞.</summary>
    public static explicit operator BigInteger(QuadDouble value) => (BigInteger)ToBigFloatExact(value).ToBigInt();

    /// <summary>The exact value; throws <see cref="OverflowException"/> for NaN or ∞.</summary>
    public static explicit operator BigRational(QuadDouble value) => ToBigFloatExact(value).ToRational();

    /// <summary>
    /// Truncates toward zero into an integer primitive. Checked throws on NaN or overflow;
    /// saturating and truncating clamp (NaN → 0), matching the BCL rule for double sources.
    /// </summary>
    private static T ToInteger<T>(QuadDouble value, ConvertMode mode) where T : INumberBase<T>
    {
        if (double.IsNaN(value._x0))
        {
            return mode == ConvertMode.Checked ? throw new OverflowException($"NaN cannot be converted to {typeof(T).Name}.") : T.Zero;
        }
        QuadDouble t = Truncate(value);
        const double two127 = 1.7014118346046923E+38, two128 = 3.4028236692093846E+38;
        double x0 = t._x0;
        if (Math.Abs(x0) <= two128)
        {
            // Every component is an integer; the lower ones are below 2^76 in magnitude.
            Int128 rest = (Int128)t._x1 + (Int128)t._x2 + (Int128)t._x3;
            if (x0 > -two127 && x0 < two127)
            {
                Int128 v = (Int128)x0 + rest;
                return mode == ConvertMode.Checked ? T.CreateChecked(v) : T.CreateSaturating(v);
            }
            if (x0 == -two127 && rest >= 0)
            {
                Int128 v = Int128.MinValue + rest;
                return mode == ConvertMode.Checked ? T.CreateChecked(v) : T.CreateSaturating(v);
            }
            if (x0 > 0 && (x0 < two128 || rest < 0))
            {
                UInt128 u = x0 < two128
                    ? (rest < 0 ? (UInt128)x0 - (UInt128)(-rest) : (UInt128)x0 + (UInt128)rest)
                    : UInt128.MaxValue - (UInt128)(-(rest + 1));
                return mode == ConvertMode.Checked ? T.CreateChecked(u) : T.CreateSaturating(u);
            }
        }        if (mode == ConvertMode.Checked)
        {
            throw new OverflowException($"Value does not fit in {typeof(T).Name}.");
        }
        return T.CreateSaturating(t._x0 > 0 ? double.PositiveInfinity : double.NegativeInfinity);
    }

    // ---------------------------------------------------------------
    // Generic-math conversions
    // ---------------------------------------------------------------

    static bool INumberBase<QuadDouble>.TryConvertFromChecked<TOther>(TOther value, out QuadDouble result) => TryConvertFrom(value, out result);
    static bool INumberBase<QuadDouble>.TryConvertFromSaturating<TOther>(TOther value, out QuadDouble result) => TryConvertFrom(value, out result);
    static bool INumberBase<QuadDouble>.TryConvertFromTruncating<TOther>(TOther value, out QuadDouble result) => TryConvertFrom(value, out result);

    static bool INumberBase<QuadDouble>.TryConvertToChecked<TOther>(QuadDouble value, [MaybeNullWhen(false)] out TOther result) =>
        TryConvertTo(value, ConvertMode.Checked, out result);

    static bool INumberBase<QuadDouble>.TryConvertToSaturating<TOther>(QuadDouble value, [MaybeNullWhen(false)] out TOther result) =>
        TryConvertTo(value, ConvertMode.Saturating, out result);

    static bool INumberBase<QuadDouble>.TryConvertToTruncating<TOther>(QuadDouble value, [MaybeNullWhen(false)] out TOther result) =>
        TryConvertTo(value, ConvertMode.Truncating, out result);

    private static bool TryConvertFrom<TOther>(TOther value, out QuadDouble result) where TOther : INumberBase<TOther>
    {
        if (typeof(TOther) == typeof(double)) { result = (double)(object)value; return true; }
        if (typeof(TOther) == typeof(float)) { result = (float)(object)value; return true; }
        if (typeof(TOther) == typeof(System.Half)) { result = (System.Half)(object)value; return true; }
        if (typeof(TOther) == typeof(sbyte)) { result = (sbyte)(object)value; return true; }
        if (typeof(TOther) == typeof(byte)) { result = (byte)(object)value; return true; }
        if (typeof(TOther) == typeof(short)) { result = (short)(object)value; return true; }
        if (typeof(TOther) == typeof(ushort)) { result = (ushort)(object)value; return true; }
        if (typeof(TOther) == typeof(char)) { result = (char)(object)value; return true; }
        if (typeof(TOther) == typeof(int)) { result = (int)(object)value; return true; }
        if (typeof(TOther) == typeof(uint)) { result = (uint)(object)value; return true; }
        if (typeof(TOther) == typeof(long)) { result = (long)(object)value; return true; }
        if (typeof(TOther) == typeof(ulong)) { result = (ulong)(object)value; return true; }
        if (typeof(TOther) == typeof(nint)) { result = (nint)(object)value; return true; }
        if (typeof(TOther) == typeof(nuint)) { result = (nuint)(object)value; return true; }
        if (typeof(TOther) == typeof(Int128)) { result = (Int128)(object)value; return true; }
        if (typeof(TOther) == typeof(UInt128)) { result = (UInt128)(object)value; return true; }
        if (typeof(TOther) == typeof(decimal)) { result = (QuadDouble)(decimal)(object)value; return true; }
        if (typeof(TOther) == typeof(BigInt)) { result = (QuadDouble)(BigInt)(object)value; return true; }
        if (typeof(TOther) == typeof(BigInteger)) { result = (QuadDouble)(BigInteger)(object)value; return true; }
        if (typeof(TOther) == typeof(BigRational)) { result = (QuadDouble)(BigRational)(object)value; return true; }
        if (typeof(TOther) == typeof(BigFloat)) { result = FromBigFloat((BigFloat)(object)value); return true; }
        result = default;
        return false;
    }

    private static bool TryConvertTo<TOther>(QuadDouble value, ConvertMode mode, [MaybeNullWhen(false)] out TOther result)
        where TOther : INumberBase<TOther>
    {
        if (BigInt.IsIntegerPrimitive<TOther>())
        {
            result = ToInteger<TOther>(value, mode);
            return true;
        }
        if (typeof(TOther) == typeof(double)) { result = (TOther)(object)(double)value; return true; }
        if (typeof(TOther) == typeof(float)) { result = (TOther)(object)(float)value; return true; }
        if (typeof(TOther) == typeof(System.Half)) { result = (TOther)(object)(System.Half)value; return true; }
        if (typeof(TOther) == typeof(decimal))
        {
            decimal m;
            if (double.IsNaN(value._x0))
            {
                m = mode == ConvertMode.Checked ? throw new OverflowException("Value does not fit in Decimal.") : 0m;
            }
            else
            {
                try
                {
                    m = (decimal)value;
                }
                catch (OverflowException) when (mode != ConvertMode.Checked)
                {
                    m = value._x0 < 0 ? decimal.MinValue : decimal.MaxValue;
                }
            }
            result = (TOther)(object)m;
            return true;
        }
        if (typeof(TOther) == typeof(BigFloat))
        {
            result = (TOther)(object)ToBigFloatExact(value);
            return true;
        }
        if (typeof(TOther) == typeof(BigInteger) || typeof(TOther) == typeof(BigInt) || typeof(TOther) == typeof(BigRational))
        {
            if (double.IsNaN(value._x0) && mode != ConvertMode.Checked)
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
            return true;
        }
        result = default;
        return false;
    }

    public static QuadDouble CreateChecked<TOther>(TOther value) where TOther : INumberBase<TOther>
    {
        if (typeof(TOther) == typeof(QuadDouble))
        {
            return (QuadDouble)(object)value;
        }
        if (!TryConvertFrom(value, out QuadDouble result) && !TOther.TryConvertToChecked(value, out result))
        {
            throw new NotSupportedException($"Conversion from {typeof(TOther).Name} to QuadDouble is not supported.");
        }
        return result;
    }

    public static QuadDouble CreateSaturating<TOther>(TOther value) where TOther : INumberBase<TOther>
    {
        if (typeof(TOther) == typeof(QuadDouble))
        {
            return (QuadDouble)(object)value;
        }
        if (!TryConvertFrom(value, out QuadDouble result) && !TOther.TryConvertToSaturating(value, out result))
        {
            throw new NotSupportedException($"Conversion from {typeof(TOther).Name} to QuadDouble is not supported.");
        }
        return result;
    }

    public static QuadDouble CreateTruncating<TOther>(TOther value) where TOther : INumberBase<TOther>
    {
        if (typeof(TOther) == typeof(QuadDouble))
        {
            return (QuadDouble)(object)value;
        }
        if (!TryConvertFrom(value, out QuadDouble result) && !TOther.TryConvertToTruncating(value, out result))
        {
            throw new NotSupportedException($"Conversion from {typeof(TOther).Name} to QuadDouble is not supported.");
        }
        return result;
    }
}
