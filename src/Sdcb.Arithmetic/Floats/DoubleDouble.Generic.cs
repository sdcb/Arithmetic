using System.Buffers.Binary;
using System.Numerics;

namespace Sdcb.Arithmetic;

/// <summary>
/// .NET generic math. <see cref="IFloatingPointIeee754{TSelf}"/> and <see cref="IBinaryNumber{TSelf}"/>
/// are deliberately not implemented: a double-double is not an IEEE 754 format (no fixed significand
/// width, not correctly rounded).
/// </summary>
public readonly partial struct DoubleDouble :
    IFloatingPoint<DoubleDouble>,
    ISignedNumber<DoubleDouble>,
    IMinMaxValue<DoubleDouble>,
    IRootFunctions<DoubleDouble>,
    IExponentialFunctions<DoubleDouble>,
    ILogarithmicFunctions<DoubleDouble>,
    IPowerFunctions<DoubleDouble>,
    ITrigonometricFunctions<DoubleDouble>,
    IHyperbolicFunctions<DoubleDouble>,
    ISpanFormattable,
    ISpanParsable<DoubleDouble>,
    IUtf8SpanFormattable,
    IUtf8SpanParsable<DoubleDouble>
{
    static DoubleDouble IAdditiveIdentity<DoubleDouble, DoubleDouble>.AdditiveIdentity => Zero;
    static DoubleDouble IMultiplicativeIdentity<DoubleDouble, DoubleDouble>.MultiplicativeIdentity => One;
    static int INumberBase<DoubleDouble>.Radix => 2;

    static bool INumberBase<DoubleDouble>.IsCanonical(DoubleDouble value) => true;
    static bool INumberBase<DoubleDouble>.IsComplexNumber(DoubleDouble value) => false;
    static bool INumberBase<DoubleDouble>.IsImaginaryNumber(DoubleDouble value) => false;
    static bool INumberBase<DoubleDouble>.IsRealNumber(DoubleDouble value) => !IsNaN(value);

    // ---------------------------------------------------------------
    // IFloatingPoint layout: the exact value |x| = significand · 2^(exponent - (bits - 1)), where
    // significand is the exact integer Hi + Lo needs (its bit length varies with the gap between
    // the words) and exponent = floor(log2 |x|).
    // ---------------------------------------------------------------

    private int LayoutExponent => _hi == 0.0 || !double.IsFinite(_hi) ? 0 : ILogB(this);

    private BigInt LayoutSignificand => _hi == 0.0 || !double.IsFinite(_hi) ? BigInt.Zero : ToRealExact().Mantissa;

    int IFloatingPoint<DoubleDouble>.GetExponentByteCount() => sizeof(int);

    int IFloatingPoint<DoubleDouble>.GetExponentShortestBitLength()
    {
        int e = LayoutExponent;
        return e >= 0 ? 32 - BitOperations.LeadingZeroCount((uint)e) : 33 - BitOperations.LeadingZeroCount((uint)~e);
    }

    int IFloatingPoint<DoubleDouble>.GetSignificandBitLength() => (int)LayoutSignificand.GetBitLength();

    int IFloatingPoint<DoubleDouble>.GetSignificandByteCount() => (int)((LayoutSignificand.GetBitLength() + 7) >> 3);

    bool IFloatingPoint<DoubleDouble>.TryWriteExponentBigEndian(Span<byte> destination, out int bytesWritten) =>
        TryWriteInt(destination, LayoutExponent, bigEndian: true, out bytesWritten);

    bool IFloatingPoint<DoubleDouble>.TryWriteExponentLittleEndian(Span<byte> destination, out int bytesWritten) =>
        TryWriteInt(destination, LayoutExponent, bigEndian: false, out bytesWritten);

    bool IFloatingPoint<DoubleDouble>.TryWriteSignificandBigEndian(Span<byte> destination, out int bytesWritten) =>
        TryWriteSignificand(destination, bigEndian: true, out bytesWritten);

    bool IFloatingPoint<DoubleDouble>.TryWriteSignificandLittleEndian(Span<byte> destination, out int bytesWritten) =>
        TryWriteSignificand(destination, bigEndian: false, out bytesWritten);

    private static bool TryWriteInt(Span<byte> destination, int value, bool bigEndian, out int bytesWritten)
    {
        bool ok = bigEndian
            ? BinaryPrimitives.TryWriteInt32BigEndian(destination, value)
            : BinaryPrimitives.TryWriteInt32LittleEndian(destination, value);
        bytesWritten = ok ? sizeof(int) : 0;
        return ok;
    }

    private bool TryWriteSignificand(Span<byte> destination, bool bigEndian, out int bytesWritten)
    {
        BigInt m = LayoutSignificand;
        int count = (int)((m.GetBitLength() + 7) >> 3);
        if (destination.Length < count)
        {
            bytesWritten = 0;
            return false;
        }
        Span<byte> d = destination.Slice(0, count);
        d.Clear();
        if (count > 0)
        {
            byte[] bytes = m.ToMagnitudeBytes(bigEndian: false);
            bytes.AsSpan(0, Math.Min(bytes.Length, count)).CopyTo(d);
            if (bigEndian)
            {
                d.Reverse();
            }
        }
        bytesWritten = count;
        return true;
    }
}
