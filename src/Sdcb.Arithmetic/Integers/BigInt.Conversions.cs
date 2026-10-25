using System.Buffers.Binary;
using System.Diagnostics.CodeAnalysis;
using System.Numerics;
using System.Runtime.InteropServices;

namespace Sdcb.Arithmetic;

/// <summary>Overflow behaviour of a generic-math conversion (CreateChecked / CreateSaturating / CreateTruncating).</summary>
internal enum ConvertMode
{
    Checked,
    Saturating,
    Truncating,
}

public readonly partial struct BigInt
{
    // ---------------------------------------------------------------
    // Primitive conversions (mirroring System.Numerics.BigInteger)
    // ---------------------------------------------------------------

    public static implicit operator BigInt(sbyte value) => new((long)value);
    public static implicit operator BigInt(byte value) => new((ulong)value);
    public static implicit operator BigInt(short value) => new((long)value);
    public static implicit operator BigInt(ushort value) => new((ulong)value);
    public static implicit operator BigInt(char value) => new((ulong)value);
    public static implicit operator BigInt(nint value) => new((long)value);
    public static implicit operator BigInt(nuint value) => new((ulong)value);

    /// <summary>Truncates toward zero; throws <see cref="OverflowException"/> for NaN or infinity.</summary>
    public static explicit operator BigInt(float value) => new((double)value);

    /// <summary>Truncates toward zero; throws <see cref="OverflowException"/> for NaN or infinity.</summary>
    public static explicit operator BigInt(Half value) => new((double)value);

    /// <summary>Truncates toward zero.</summary>
    public static explicit operator BigInt(decimal value) => FromDecimal(value);

    public static explicit operator sbyte(BigInt value) => checked((sbyte)(long)value);
    public static explicit operator byte(BigInt value) => checked((byte)(ulong)value);
    public static explicit operator short(BigInt value) => checked((short)(long)value);
    public static explicit operator ushort(BigInt value) => checked((ushort)(ulong)value);
    public static explicit operator char(BigInt value) => checked((char)(ulong)value);
    public static explicit operator nint(BigInt value) => checked((nint)(long)value);
    public static explicit operator nuint(BigInt value) => checked((nuint)(ulong)value);

    /// <summary>Correctly rounded (ties to even); ±Infinity on overflow.</summary>
    public static explicit operator float(BigInt value) => value.ToSingle();

    /// <summary>Correctly rounded (ties to even); ±Infinity on overflow.</summary>
    public static explicit operator Half(BigInt value) => value.ToHalf();

    /// <summary>Exact; throws <see cref="OverflowException"/> when |value| &gt;= 2^96.</summary>
    public static explicit operator decimal(BigInt value) =>
        value.TryToDecimal(out decimal d) ? d : throw new OverflowException("Value does not fit in Decimal.");

    /// <summary>Converts to float with correct round-to-nearest-even; returns ±Infinity on overflow.</summary>
    public float ToSingle() => (float)ToIeee(IeeeFormat.Binary32);

    /// <summary>Converts to Half with correct round-to-nearest-even; returns ±Infinity on overflow.</summary>
    public Half ToHalf() => (Half)ToIeee(IeeeFormat.Binary16);

    private double ToIeee(IeeeFormat format) =>
        _size == 0 ? 0.0 : Real.RoundToIeee(Abs(), 0, sticky: false, _size < 0, format, RealRounding.NearestEven);

    private bool TryToDecimal(out decimal result)
    {
        int n = LimbCount;
        if (n > 2 || (n == 2 && _d![1] > uint.MaxValue))
        {
            result = 0m;
            return false;
        }
        ulong lo = n > 0 ? _d![0] : 0;
        ulong hi = n > 1 ? _d![1] : 0;
        result = new decimal((int)(uint)lo, (int)(uint)(lo >> 32), (int)(uint)hi, _size < 0, 0);
        return true;
    }

    private static BigInt FromDecimal(decimal value)
    {
        Span<int> bits = stackalloc int[4];
        decimal.GetBits(decimal.Truncate(value), bits);
        UInt128 mant = ((UInt128)(uint)bits[2] << 64) | ((UInt128)(uint)bits[1] << 32) | (uint)bits[0];
        int scale = (bits[3] >> 16) & 0xFF;
        while (scale-- > 0)
        {
            mant /= 10;
        }
        BigInt m = new(mant);
        return bits[3] < 0 ? -m : m;
    }

    // ---------------------------------------------------------------
    // System.Numerics.BigInteger interop
    // ---------------------------------------------------------------

    /// <summary>Lossless conversion from <see cref="BigInteger"/>.</summary>
    public static implicit operator BigInt(BigInteger value)
    {
        if (value.IsZero)
        {
            return default;
        }
        int sign = value.Sign;
        BigInteger mag = BigInteger.Abs(value);
        int byteCount = mag.GetByteCount(isUnsigned: true);
        int limbs = (byteCount + 7) >> 3;
        ulong[] d = new ulong[limbs];
        mag.TryWriteBytes(MemoryMarshal.AsBytes(d.AsSpan()), out _, isUnsigned: true, isBigEndian: false);
        if (!BitConverter.IsLittleEndian)
        {
            for (int i = 0; i < d.Length; i++)
            {
                d[i] = BinaryPrimitives.ReverseEndianness(d[i]);
            }
        }
        return FromMagnitude(d, limbs, sign);
    }

    /// <summary>Lossless conversion to <see cref="BigInteger"/> (explicit to keep mixed BigInt/BigInteger expressions unambiguous).</summary>
    public static explicit operator BigInteger(BigInt value)
    {
        if (value._size == 0)
        {
            return BigInteger.Zero;
        }
        BigInteger m = BitConverter.IsLittleEndian
            ? new BigInteger(MemoryMarshal.AsBytes(value.Limbs), isUnsigned: true, isBigEndian: false)
            : new BigInteger(value.ToMagnitudeBytes(), isUnsigned: true, isBigEndian: false);
        return value._size < 0 ? -m : m;
    }

    // ---------------------------------------------------------------
    // Generic-math conversion plumbing
    // ---------------------------------------------------------------

    internal static bool IsSignedIntegerPrimitive<T>() =>
        typeof(T) == typeof(sbyte) || typeof(T) == typeof(short) || typeof(T) == typeof(int) ||
        typeof(T) == typeof(long) || typeof(T) == typeof(nint);

    internal static bool IsUnsignedIntegerPrimitive<T>() =>
        typeof(T) == typeof(byte) || typeof(T) == typeof(ushort) || typeof(T) == typeof(char) ||
        typeof(T) == typeof(uint) || typeof(T) == typeof(ulong) || typeof(T) == typeof(nuint);

    internal static bool IsIntegerPrimitive<T>() =>
        IsSignedIntegerPrimitive<T>() || IsUnsignedIntegerPrimitive<T>() ||
        typeof(T) == typeof(Int128) || typeof(T) == typeof(UInt128);

    internal static bool IsBinaryFloatPrimitive<T>() =>
        typeof(T) == typeof(double) || typeof(T) == typeof(float) || typeof(T) == typeof(Half);

    /// <summary>Converts an integer primitive exactly; false for any other type.</summary>
    internal static bool TryFromIntegerPrimitive<TOther>(TOther value, out BigInt result)
        where TOther : INumberBase<TOther>
    {
        if (IsSignedIntegerPrimitive<TOther>())
        {
            result = new BigInt(long.CreateTruncating(value));
            return true;
        }
        if (IsUnsignedIntegerPrimitive<TOther>())
        {
            result = new BigInt(ulong.CreateTruncating(value));
            return true;
        }
        if (typeof(TOther) == typeof(Int128))
        {
            result = new BigInt((Int128)(object)value);
            return true;
        }
        if (typeof(TOther) == typeof(UInt128))
        {
            result = new BigInt((UInt128)(object)value);
            return true;
        }
        result = default;
        return false;
    }

    private static bool TryConvertFrom<TOther>(TOther value, ConvertMode mode, out BigInt result)
        where TOther : INumberBase<TOther>
    {
        if (TryFromIntegerPrimitive(value, out result))
        {
            return true;
        }
        if (IsBinaryFloatPrimitive<TOther>())
        {
            double d = double.CreateTruncating(value); // exact for float and Half
            if (double.IsNaN(d) && mode != ConvertMode.Checked)
            {
                result = default;
                return true;
            }
            result = new BigInt(d); // throws OverflowException for NaN / infinity
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
        result = default;
        return false;
    }

    private static bool TryConvertTo<TOther>(BigInt value, ConvertMode mode, [MaybeNullWhen(false)] out TOther result)
        where TOther : INumberBase<TOther>
    {
        if (IsIntegerPrimitive<TOther>())
        {
            result = ConvertToIntegerPrimitive<TOther>(value, mode);
            return true;
        }
        if (IsBinaryFloatPrimitive<TOther>())
        {
            IeeeFormat f = typeof(TOther) == typeof(double) ? IeeeFormat.Binary64
                : typeof(TOther) == typeof(float) ? IeeeFormat.Binary32 : IeeeFormat.Binary16;
            result = TOther.CreateTruncating(value.ToIeee(f)); // exact: already rounded to the target format
            return true;
        }
        if (typeof(TOther) == typeof(decimal))
        {
            if (!value.TryToDecimal(out decimal d))
            {
                if (mode == ConvertMode.Checked)
                {
                    throw new OverflowException("Value does not fit in Decimal.");
                }
                d = value._size < 0 ? decimal.MinValue : decimal.MaxValue;
            }
            result = (TOther)(object)d;
            return true;
        }
        if (typeof(TOther) == typeof(BigInteger))
        {
            result = (TOther)(object)(BigInteger)value;
            return true;
        }
        result = default;
        return false;
    }

    /// <summary>Converts to any integer primitive with checked / saturating / two's-complement-truncating semantics.</summary>
    internal static TOther ConvertToIntegerPrimitive<TOther>(BigInt value, ConvertMode mode)
        where TOther : INumberBase<TOther>
    {
        if (mode == ConvertMode.Truncating)
        {
            return TOther.CreateTruncating(value.LowBits128());
        }
        if (value.LimbCount <= 2)
        {
            UInt128 mag = value.Abs().LowBits128();
            if (value._size >= 0)
            {
                return mode == ConvertMode.Checked ? TOther.CreateChecked(mag) : TOther.CreateSaturating(mag);
            }
            if (mag <= (UInt128)Int128.MaxValue + 1)
            {
                Int128 s = unchecked(-(Int128)mag);
                return mode == ConvertMode.Checked ? TOther.CreateChecked(s) : TOther.CreateSaturating(s);
            }
        }
        if (mode == ConvertMode.Checked)
        {
            throw new OverflowException($"Value does not fit in {typeof(TOther).Name}.");
        }
        return value._size > 0 ? TOther.CreateSaturating(UInt128.MaxValue) : TOther.CreateSaturating(Int128.MinValue);
    }

    /// <summary>The low 128 bits of the infinite two's-complement representation.</summary>
    internal UInt128 LowBits128()
    {
        int n = LimbCount;
        ulong lo = n > 0 ? _d![0] : 0;
        ulong hi = n > 1 ? _d![1] : 0;
        UInt128 mag = ((UInt128)hi << 64) | lo;
        return _size < 0 ? unchecked(UInt128.Zero - mag) : mag;
    }
}
