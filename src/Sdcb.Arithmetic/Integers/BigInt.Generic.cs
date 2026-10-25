using System.Diagnostics.CodeAnalysis;
using System.Globalization;
using System.Numerics;
using System.Runtime.InteropServices;
using Sdcb.Arithmetic.Core;

namespace Sdcb.Arithmetic;

/// <summary>
/// .NET generic math (<see cref="IBinaryInteger{TSelf}"/>, <see cref="ISignedNumber{TSelf}"/>).
/// Fixed-width notions (leading zeros, rotation, unsigned shift, popcount of negatives) use
/// the minimal two's-complement representation rounded up to whole 64-bit limbs, the same
/// convention <see cref="BigInteger"/> uses with its 32-bit words.
/// </summary>
public readonly partial struct BigInt :
    IBinaryInteger<BigInt>,
    ISignedNumber<BigInt>,
    IUtf8SpanFormattable,
    IUtf8SpanParsable<BigInt>
{
    static BigInt IAdditiveIdentity<BigInt, BigInt>.AdditiveIdentity => default;
    static BigInt IMultiplicativeIdentity<BigInt, BigInt>.MultiplicativeIdentity => One;

    /// <summary>-1 (the <see cref="ISignedNumber{TSelf}"/> name for <see cref="MinusOne"/>).</summary>
    public static BigInt NegativeOne => MinusOne;

    static int INumberBase<BigInt>.Radix => 2;

    // ---------------------------------------------------------------
    // INumberBase predicates
    // ---------------------------------------------------------------

    static bool INumberBase<BigInt>.IsCanonical(BigInt value) => true;
    static bool INumberBase<BigInt>.IsComplexNumber(BigInt value) => false;
    public static bool IsEvenInteger(BigInt value) => value.IsEven;
    static bool INumberBase<BigInt>.IsFinite(BigInt value) => true;
    static bool INumberBase<BigInt>.IsImaginaryNumber(BigInt value) => false;
    static bool INumberBase<BigInt>.IsInfinity(BigInt value) => false;
    static bool INumberBase<BigInt>.IsInteger(BigInt value) => true;
    static bool INumberBase<BigInt>.IsNaN(BigInt value) => false;
    static bool INumberBase<BigInt>.IsNegative(BigInt value) => value._size < 0;
    static bool INumberBase<BigInt>.IsNegativeInfinity(BigInt value) => false;
    static bool INumberBase<BigInt>.IsNormal(BigInt value) => value._size != 0;
    public static bool IsOddInteger(BigInt value) => value.IsOdd;
    static bool INumberBase<BigInt>.IsPositive(BigInt value) => value._size >= 0;
    static bool INumberBase<BigInt>.IsPositiveInfinity(BigInt value) => false;
    static bool INumberBase<BigInt>.IsRealNumber(BigInt value) => true;
    static bool INumberBase<BigInt>.IsSubnormal(BigInt value) => false;
    static bool INumberBase<BigInt>.IsZero(BigInt value) => value._size == 0;

    /// <summary>True for positive powers of two.</summary>
    public static bool IsPow2(BigInt value) => value.IsPowerOfTwo;

    // ---------------------------------------------------------------
    // INumber / magnitude helpers
    // ---------------------------------------------------------------

    static int INumber<BigInt>.Sign(BigInt value) => value.Sign;

    public static BigInt Clamp(BigInt value, BigInt min, BigInt max)
    {
        if (min > max)
        {
            throw new ArgumentException("min must be less than or equal to max.", nameof(min));
        }
        return value < min ? min : value > max ? max : value;
    }

    /// <summary>|value| with the sign of <paramref name="sign"/> (non-negative when sign is zero).</summary>
    public static BigInt CopySign(BigInt value, BigInt sign) => sign._size < 0 ? -value.Abs() : value.Abs();

    public static BigInt MaxNumber(BigInt x, BigInt y) => Max(x, y);
    public static BigInt MinNumber(BigInt x, BigInt y) => Min(x, y);

    public static BigInt MaxMagnitude(BigInt x, BigInt y)
    {
        int c = CompareAbs(x, y);
        return c > 0 ? x : c < 0 ? y : x._size < 0 ? y : x;
    }

    public static BigInt MinMagnitude(BigInt x, BigInt y)
    {
        int c = CompareAbs(x, y);
        return c < 0 ? x : c > 0 ? y : x._size < 0 ? x : y;
    }

    public static BigInt MaxMagnitudeNumber(BigInt x, BigInt y) => MaxMagnitude(x, y);
    public static BigInt MinMagnitudeNumber(BigInt x, BigInt y) => MinMagnitude(x, y);

    // ---------------------------------------------------------------
    // IBinaryInteger
    // ---------------------------------------------------------------

    /// <summary>Truncating division returning (quotient, remainder) like <see cref="Math.DivRem(int, int)"/>.</summary>
    public static (BigInt Quotient, BigInt Remainder) DivRem(BigInt left, BigInt right)
    {
        DivRem(left, right, out BigInt q, out BigInt r);
        return (q, r);
    }

    /// <summary>floor(log2(value)); 0 for zero. Throws for negative values.</summary>
    public static BigInt Log2(BigInt value)
    {
        if (value._size < 0)
        {
            throw new ArgumentOutOfRangeException(nameof(value), "Log2 is undefined for negative values.");
        }
        return value._size == 0 ? default : new BigInt(value.GetBitLength() - 1);
    }

    /// <summary>Leading zeros within the top 64-bit limb; 0 for negative values, 64 for zero.</summary>
    public static BigInt LeadingZeroCount(BigInt value)
    {
        if (value._size < 0)
        {
            return default;
        }
        if (value._size == 0)
        {
            return new BigInt(64);
        }
        return new BigInt(BitOperations.LeadingZeroCount(value._d![value._size - 1]));
    }

    /// <summary>Set bits of the two's-complement representation (over whole limbs for negative values).</summary>
    public static BigInt PopCount(BigInt value)
    {
        if (value._size >= 0)
        {
            return new BigInt(Nat.PopCount(value.Limbs, value.LimbCount));
        }
        int n = value.TwosComplementLimbCount();
        ulong[] t = value.TwosComplementLimbs(n);
        return new BigInt(Nat.PopCount(t, n));
    }

    /// <summary>Trailing zero bits; 64 for zero.</summary>
    public static BigInt TrailingZeroCount(BigInt value) =>
        new(value._size == 0 ? 64 : Nat.TrailingZeroCount(value.Limbs, value.LimbCount));

    /// <summary>
    /// Rotates within the minimal whole-limb two's-complement width. Non-negative inputs are
    /// treated as unsigned (the result stays non-negative); negative inputs as signed.
    /// </summary>
    public static BigInt RotateLeft(BigInt value, int rotateAmount)
    {
        if (value._size == 0)
        {
            return value;
        }
        int n = value.TwosComplementLimbCount();
        long width = 64L * n;
        long r = ((rotateAmount % width) + width) % width;
        if (r == 0)
        {
            return value;
        }
        BigInt u = FromMagnitude(value.TwosComplementLimbs(n), n, 1);
        BigInt mask = (One << (int)width) - One;
        BigInt rotated = (u.ShiftLeft(r) | u.ShiftRightTruncate(width - r)) & mask;
        if (value._size < 0 && rotated.TestBit(width - 1))
        {
            rotated -= One << (int)width;
        }
        return rotated;
    }

    public static BigInt RotateRight(BigInt value, int rotateAmount)
    {
        if (value._size == 0)
        {
            return value;
        }
        long width = 64L * value.TwosComplementLimbCount();
        long r = ((rotateAmount % width) + width) % width;
        return RotateLeft(value, (int)((width - r) % width));
    }

    /// <summary>Logical right shift over the minimal whole-limb two's-complement width.</summary>
    public static BigInt operator >>>(BigInt value, int shiftAmount)
    {
        if (value._size >= 0 || shiftAmount == 0)
        {
            return value >> shiftAmount;
        }
        if (shiftAmount < 0)
        {
            return value << -shiftAmount;
        }
        int n = value.TwosComplementLimbCount();
        if ((long)shiftAmount >= 64L * n)
        {
            return default;
        }
        return FromMagnitude(value.TwosComplementLimbs(n), n, 1).ShiftRightTruncate(shiftAmount);
    }

    /// <summary>Bits in the shortest two's-complement form, excluding the sign bit for non-negative values (like BigInteger).</summary>
    public long GetShortestBitLength()
    {
        if (_size >= 0)
        {
            return GetBitLength();
        }
        long bits = GetBitLength();
        return IsMagnitudePowerOfTwo() ? bits : bits + 1;
    }

    /// <summary>Bytes written by <see cref="TryWriteLittleEndian"/>: the shortest signed two's-complement form (at least one byte).</summary>
    public int GetByteCount()
    {
        long bits = _size >= 0 ? GetBitLength() + 1 : GetShortestBitLength();
        return (int)((bits + 7) >> 3);
    }

    int IBinaryInteger<BigInt>.GetShortestBitLength() => checked((int)GetShortestBitLength());

    /// <summary>Writes the shortest signed two's-complement little-endian form (<see cref="GetByteCount"/> bytes).</summary>
    public bool TryWriteLittleEndian(Span<byte> destination, out int bytesWritten)
    {
        int count = GetByteCount();
        if (destination.Length < count)
        {
            bytesWritten = 0;
            return false;
        }
        int n = (count + 7) >> 3;
        ulong[] t = TwosComplementLimbs(n);
        if (BitConverter.IsLittleEndian)
        {
            MemoryMarshal.AsBytes(t.AsSpan()).Slice(0, count).CopyTo(destination);
        }
        else
        {
            for (int i = 0; i < count; i++)
            {
                destination[i] = (byte)(t[i >> 3] >> ((i & 7) * 8));
            }
        }
        bytesWritten = count;
        return true;
    }

    /// <summary>Writes the shortest signed two's-complement big-endian form (<see cref="GetByteCount"/> bytes).</summary>
    public bool TryWriteBigEndian(Span<byte> destination, out int bytesWritten)
    {
        if (!TryWriteLittleEndian(destination, out bytesWritten))
        {
            return false;
        }
        destination.Slice(0, bytesWritten).Reverse();
        return true;
    }

    /// <summary>Reads a little-endian integer; signed input uses two's complement.</summary>
    public static BigInt ReadLittleEndian(ReadOnlySpan<byte> source, bool isUnsigned)
    {
        if (source.IsEmpty)
        {
            return default;
        }
        int n = (source.Length + 7) >> 3;
        ulong[] d = new ulong[n];
        if (BitConverter.IsLittleEndian)
        {
            source.CopyTo(MemoryMarshal.AsBytes(d.AsSpan()));
        }
        else
        {
            for (int i = 0; i < source.Length; i++)
            {
                d[i >> 3] |= (ulong)source[i] << ((i & 7) * 8);
            }
        }
        BigInt v = FromMagnitude(d, n, 1);
        if (!isUnsigned && (source[^1] & 0x80) != 0)
        {
            v -= One << (8 * source.Length);
        }
        return v;
    }

    /// <summary>Reads a big-endian integer; signed input uses two's complement.</summary>
    public static BigInt ReadBigEndian(ReadOnlySpan<byte> source, bool isUnsigned)
    {
        byte[] tmp = source.ToArray();
        Array.Reverse(tmp);
        return ReadLittleEndian(tmp, isUnsigned);
    }

    static bool IBinaryInteger<BigInt>.TryReadLittleEndian(ReadOnlySpan<byte> source, bool isUnsigned, out BigInt value)
    {
        value = ReadLittleEndian(source, isUnsigned);
        return true;
    }

    static bool IBinaryInteger<BigInt>.TryReadBigEndian(ReadOnlySpan<byte> source, bool isUnsigned, out BigInt value)
    {
        value = ReadBigEndian(source, isUnsigned);
        return true;
    }

    private bool IsMagnitudePowerOfTwo() => _size != 0 && Nat.TrailingZeroCount(Limbs, LimbCount) == GetBitLength() - 1;

    /// <summary>Limbs needed for the minimal two's-complement form (magnitude limbs for non-negative values).</summary>
    private int TwosComplementLimbCount() => _size >= 0 ? LimbCount : (int)((GetShortestBitLength() + 63) >> 6);

    /// <summary>The two's-complement form truncated / sign-extended to len limbs.</summary>
    private ulong[] TwosComplementLimbs(int len)
    {
        ulong[] t = new ulong[len];
        ReadOnlySpan<ulong> m = Limbs;
        m.Slice(0, Math.Min(m.Length, len)).CopyTo(t);
        if (_size < 0)
        {
            NegateTwosComplementInPlace(t);
        }
        return t;
    }

    // ---------------------------------------------------------------
    // Fused operations
    // ---------------------------------------------------------------

    /// <summary>a * b + c computed into a single result buffer.</summary>
    public static BigInt MultiplyAdd(BigInt a, BigInt b, BigInt c) => MulAddCore(a, b, c, negateC: false);

    /// <summary>a * b - c computed into a single result buffer.</summary>
    public static BigInt MultiplySubtract(BigInt a, BigInt b, BigInt c) => MulAddCore(a, b, c, negateC: true);

    /// <summary>(a * b) mod |m| in [0, |m|) (mpz_mul + mpz_mod).</summary>
    public static BigInt ModMultiply(BigInt a, BigInt b, BigInt m)
    {
        if (m._size == 0)
        {
            throw new DivideByZeroException();
        }
        BigInt r = (a * b) % m;
        return r._size < 0 ? r + m.Abs() : r;
    }

    private static BigInt MulAddCore(BigInt a, BigInt b, BigInt c, bool negateC)
    {
        if (a._size == 0 || b._size == 0)
        {
            return negateC ? -c : c;
        }
        if (c._size == 0)
        {
            return a * b;
        }
        int an = a.LimbCount, bn = b.LimbCount, cn = c.LimbCount;
        int pn = an + bn;
        int len = Math.Max(pn, cn) + 1;
        ulong[] d = NewLimbs(len);
        if (an >= bn)
        {
            Nat.Mul(d, a.Limbs, an, b.Limbs, bn);
        }
        else
        {
            Nat.Mul(d, b.Limbs, bn, a.Limbs, an);
        }
        d.AsSpan(pn).Clear();
        int pSign = (a._size ^ b._size) < 0 ? -1 : 1;
        int cSign = (c._size < 0) != negateC ? -1 : 1;
        if (pSign == cSign)
        {
            Nat.AddInPlace(d, len, c.Limbs, cn);
            return FromMagnitude(d, len, pSign);
        }
        pn = Nat.TrimmedLength(d, pn);
        int cmp = Nat.CompareNormalized(d, pn, c.Limbs, cn);
        if (cmp == 0)
        {
            return default;
        }
        if (cmp > 0)
        {
            Nat.SubInPlace(d, pn, c.Limbs, cn);
            return FromMagnitude(d, pn, pSign);
        }
        Nat.Sub(d, c.Limbs, cn, d, pn);
        return FromMagnitude(d, cn, cSign);
    }

    // ---------------------------------------------------------------
    // Generic conversions
    // ---------------------------------------------------------------

    static bool INumberBase<BigInt>.TryConvertFromChecked<TOther>(TOther value, out BigInt result) =>
        TryConvertFrom(value, ConvertMode.Checked, out result);

    static bool INumberBase<BigInt>.TryConvertFromSaturating<TOther>(TOther value, out BigInt result) =>
        TryConvertFrom(value, ConvertMode.Saturating, out result);

    static bool INumberBase<BigInt>.TryConvertFromTruncating<TOther>(TOther value, out BigInt result) =>
        TryConvertFrom(value, ConvertMode.Truncating, out result);

    static bool INumberBase<BigInt>.TryConvertToChecked<TOther>(BigInt value, [MaybeNullWhen(false)] out TOther result) =>
        TryConvertTo(value, ConvertMode.Checked, out result);

    static bool INumberBase<BigInt>.TryConvertToSaturating<TOther>(BigInt value, [MaybeNullWhen(false)] out TOther result) =>
        TryConvertTo(value, ConvertMode.Saturating, out result);

    static bool INumberBase<BigInt>.TryConvertToTruncating<TOther>(BigInt value, [MaybeNullWhen(false)] out TOther result) =>
        TryConvertTo(value, ConvertMode.Truncating, out result);

    // ---------------------------------------------------------------
    // NumberStyles parsing
    // ---------------------------------------------------------------

    /// <summary>
    /// Parses with .NET <see cref="NumberStyles"/> semantics. Hex and binary digits are read as
    /// two's complement (a leading digit with the top bit set is negative), like BigInteger.
    /// Decimal points and exponents are accepted when the style allows them and the value is integral.
    /// </summary>
    public static BigInt Parse(ReadOnlySpan<char> s, NumberStyles style, IFormatProvider? provider) =>
        TryParse(s, style, provider, out BigInt r) ? r : throw new FormatException("The value could not be parsed as a BigInt.");

    public static BigInt Parse(string s, NumberStyles style, IFormatProvider? provider)
    {
        ArgumentNullException.ThrowIfNull(s);
        return Parse(s.AsSpan(), style, provider);
    }

    public static bool TryParse(string? s, NumberStyles style, IFormatProvider? provider, out BigInt result) =>
        TryParse(s.AsSpan(), style, provider, out result);

    public static bool TryParse(ReadOnlySpan<char> s, NumberStyles style, IFormatProvider? provider, out BigInt result)
    {
        result = default;
        if (!NumberText.TryNormalize(s, style, provider, out string text, out NumberText.Radix radix))
        {
            return false;
        }
        switch (radix)
        {
            case NumberText.Radix.Hex:
                if (!TryParse(text, 16, out result))
                {
                    return false;
                }
                if (HexDigitValue(text[0]) >= 8)
                {
                    result -= One << (4 * text.Length);
                }
                return true;
            case NumberText.Radix.Binary:
                if (!TryParse(text, 2, out result))
                {
                    return false;
                }
                if (text[0] == '1')
                {
                    result -= One << text.Length;
                }
                return true;
        }
        if (text.AsSpan().IndexOfAny('.', 'e') < 0)
        {
            return TryParse(text, 10, out result);
        }
        if (!BigRational.TryParse(text, CultureInfo.InvariantCulture, out BigRational q) || !q.IsInteger)
        {
            return false;
        }
        result = q.Numerator;
        return true;
    }

    private static int HexDigitValue(char c) => c <= '9' ? c - '0' : (c | 0x20) - 'a' + 10;

    // ---------------------------------------------------------------
    // UTF-8
    // ---------------------------------------------------------------

    public bool TryFormat(Span<byte> utf8Destination, out int bytesWritten, ReadOnlySpan<char> format, IFormatProvider? provider) =>
        NumberText.TryWriteUtf8(ToString(format.IsEmpty ? null : format.ToString(), provider), utf8Destination, out bytesWritten);

    public static BigInt Parse(ReadOnlySpan<byte> utf8Text, IFormatProvider? provider) =>
        TryParse(utf8Text, provider, out BigInt r) ? r : throw new FormatException("The value could not be parsed as a BigInt.");

    public static bool TryParse(ReadOnlySpan<byte> utf8Text, IFormatProvider? provider, out BigInt result) =>
        TryParse(utf8Text, NumberStyles.Integer, provider, out result);

    public static BigInt Parse(ReadOnlySpan<byte> utf8Text, NumberStyles style, IFormatProvider? provider) =>
        TryParse(utf8Text, style, provider, out BigInt r) ? r : throw new FormatException("The value could not be parsed as a BigInt.");

    public static bool TryParse(ReadOnlySpan<byte> utf8Text, NumberStyles style, IFormatProvider? provider, out BigInt result) =>
        NumberText.TryParseUtf8(utf8Text, (ReadOnlySpan<char> s, out BigInt r) => TryParse(s, style, provider, out r), out result);

    // ---------------------------------------------------------------
    // Public CreateChecked / CreateSaturating / CreateTruncating (like BigInteger)
    // ---------------------------------------------------------------

    public static BigInt CreateChecked<TOther>(TOther value) where TOther : INumberBase<TOther>
    {
        if (typeof(TOther) == typeof(BigInt))
        {
            return (BigInt)(object)value;
        }
        if (!TryConvertFrom(value, ConvertMode.Checked, out BigInt result) && !TOther.TryConvertToChecked(value, out result))
        {
            throw new NotSupportedException($"Conversion from {typeof(TOther).Name} to BigInt is not supported.");
        }
        return result;
    }

    public static BigInt CreateSaturating<TOther>(TOther value) where TOther : INumberBase<TOther>
    {
        if (typeof(TOther) == typeof(BigInt))
        {
            return (BigInt)(object)value;
        }
        if (!TryConvertFrom(value, ConvertMode.Saturating, out BigInt result) && !TOther.TryConvertToSaturating(value, out result))
        {
            throw new NotSupportedException($"Conversion from {typeof(TOther).Name} to BigInt is not supported.");
        }
        return result;
    }

    public static BigInt CreateTruncating<TOther>(TOther value) where TOther : INumberBase<TOther>
    {
        if (typeof(TOther) == typeof(BigInt))
        {
            return (BigInt)(object)value;
        }
        if (!TryConvertFrom(value, ConvertMode.Truncating, out BigInt result) && !TOther.TryConvertToTruncating(value, out result))
        {
            throw new NotSupportedException($"Conversion from {typeof(TOther).Name} to BigInt is not supported.");
        }
        return result;
    }
}