using System.Numerics;
using System.Runtime.CompilerServices;
using static Sdcb.Arithmetic.Core.ErrorFree;

namespace Sdcb.Arithmetic;

/// <summary>
/// .NET generic math for <see cref="QuadDouble"/>: <see cref="IFloatingPoint{TSelf}"/> plus the root,
/// exponential, logarithmic, power, trigonometric and hyperbolic function interfaces. QuadDouble is
/// not an IEEE binary format, so <see cref="IBinaryNumber{TSelf}"/> and
/// <see cref="IFloatingPointIeee754{TSelf}"/> are deliberately not implemented.
/// </summary>
public readonly partial struct QuadDouble :
    IFloatingPoint<QuadDouble>,
    ISignedNumber<QuadDouble>,
    IMinMaxValue<QuadDouble>,
    IRootFunctions<QuadDouble>,
    IExponentialFunctions<QuadDouble>,
    ILogarithmicFunctions<QuadDouble>,
    IPowerFunctions<QuadDouble>,
    ITrigonometricFunctions<QuadDouble>,
    IHyperbolicFunctions<QuadDouble>,
    IUtf8SpanFormattable,
    IUtf8SpanParsable<QuadDouble>,
    ISpanFormattable,
    ISpanParsable<QuadDouble>,
    IEquatable<QuadDouble>,
    IComparable<QuadDouble>,
    IComparable
{
    static QuadDouble IAdditiveIdentity<QuadDouble, QuadDouble>.AdditiveIdentity => Zero;
    static QuadDouble IMultiplicativeIdentity<QuadDouble, QuadDouble>.MultiplicativeIdentity => One;
    static int INumberBase<QuadDouble>.Radix => 2;

    static bool INumberBase<QuadDouble>.IsCanonical(QuadDouble value) => true;
    static bool INumberBase<QuadDouble>.IsComplexNumber(QuadDouble value) => false;
    static bool INumberBase<QuadDouble>.IsImaginaryNumber(QuadDouble value) => false;

    /// <summary>True unless NaN.</summary>
    public static bool IsRealNumber(QuadDouble value) => !double.IsNaN(value._x0);

    // ---------------------------------------------------------------
    // Sign, magnitude, min / max
    // ---------------------------------------------------------------

    /// <summary>|value| (exact).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static QuadDouble Abs(QuadDouble value) => double.IsNegative(value._x0) ? -value : value;

    /// <summary>The magnitude of value with the sign bit of sign.</summary>
    public static QuadDouble CopySign(QuadDouble value, QuadDouble sign) =>
        double.IsNegative(value._x0) == double.IsNegative(sign._x0) ? value : -value;

    /// <summary>-1, 0 or +1; throws <see cref="ArithmeticException"/> for NaN (like <see cref="Math.Sign(double)"/>).</summary>
    public static int Sign(QuadDouble value) => Math.Sign(value._x0);

    /// <summary>Clamps value to [min, max]; NaN value stays NaN.</summary>
    public static QuadDouble Clamp(QuadDouble value, QuadDouble min, QuadDouble max)
    {
        if (!(min <= max))
        {
            throw new ArgumentException("min must be less than or equal to max.", nameof(min));
        }
        if (IsNaN(value))
        {
            return value;
        }
        return value < min ? min : value > max ? max : value;
    }

    /// <summary>The larger value (IEEE 754:2019 maximum: NaN propagates, +0 &gt; -0), like <see cref="Math.Max(double, double)"/>.</summary>
    public static QuadDouble Max(QuadDouble x, QuadDouble y)
    {
        if (IsNaN(x))
        {
            return x;
        }
        if (IsNaN(y))
        {
            return y;
        }
        int c = CompareNonNaN(x, y);
        return c > 0 ? x : c < 0 ? y : double.IsNegative(x._x0) ? y : x;
    }

    /// <summary>The smaller value (IEEE 754:2019 minimum: NaN propagates, -0 &lt; +0), like <see cref="Math.Min(double, double)"/>.</summary>
    public static QuadDouble Min(QuadDouble x, QuadDouble y)
    {
        if (IsNaN(x))
        {
            return x;
        }
        if (IsNaN(y))
        {
            return y;
        }
        int c = CompareNonNaN(x, y);
        return c < 0 ? x : c > 0 ? y : double.IsNegative(x._x0) ? x : y;
    }

    /// <summary>IEEE 754:2019 maximumNumber: a number wins over NaN.</summary>
    public static QuadDouble MaxNumber(QuadDouble x, QuadDouble y) => IsNaN(y) ? x : IsNaN(x) ? y : Max(x, y);

    /// <summary>IEEE 754:2019 minimumNumber: a number wins over NaN.</summary>
    public static QuadDouble MinNumber(QuadDouble x, QuadDouble y) => IsNaN(y) ? x : IsNaN(x) ? y : Min(x, y);

    private static int CompareMagnitude(QuadDouble x, QuadDouble y) => CompareNonNaN(Abs(x), Abs(y));

    /// <summary>The value with the larger magnitude (NaN propagates; on ties the positive one).</summary>
    public static QuadDouble MaxMagnitude(QuadDouble x, QuadDouble y)
    {
        if (IsNaN(x))
        {
            return x;
        }
        if (IsNaN(y))
        {
            return y;
        }
        int c = CompareMagnitude(x, y);
        return c > 0 ? x : c < 0 ? y : double.IsNegative(x._x0) ? y : x;
    }

    /// <summary>The value with the smaller magnitude (NaN propagates; on ties the negative one).</summary>
    public static QuadDouble MinMagnitude(QuadDouble x, QuadDouble y)
    {
        if (IsNaN(x))
        {
            return x;
        }
        if (IsNaN(y))
        {
            return y;
        }
        int c = CompareMagnitude(x, y);
        return c < 0 ? x : c > 0 ? y : double.IsNegative(x._x0) ? x : y;
    }

    /// <summary>Like <see cref="MaxMagnitude"/> but a number wins over NaN.</summary>
    public static QuadDouble MaxMagnitudeNumber(QuadDouble x, QuadDouble y)
    {
        if (IsNaN(y))
        {
            return x;
        }
        if (IsNaN(x))
        {
            return y;
        }
        return MaxMagnitude(x, y);
    }

    /// <summary>Like <see cref="MinMagnitude"/> but a number wins over NaN.</summary>
    public static QuadDouble MinMagnitudeNumber(QuadDouble x, QuadDouble y)
    {
        if (IsNaN(y))
        {
            return x;
        }
        if (IsNaN(x))
        {
            return y;
        }
        return MinMagnitude(x, y);
    }

    // ---------------------------------------------------------------
    // Scaling and exponent
    // ---------------------------------------------------------------

    /// <summary>x · 2^n, scaling each component (exact unless the result over- or underflows).</summary>
    public static QuadDouble ScaleB(QuadDouble x, int n) =>
        Make(Math.ScaleB(x._x0, n), Math.ScaleB(x._x1, n), Math.ScaleB(x._x2, n), Math.ScaleB(x._x3, n));

    /// <summary>floor(log2 |x|); <see cref="int.MinValue"/> for zero, <see cref="int.MaxValue"/> for NaN and ∞ (like <see cref="Math.ILogB"/>).</summary>
    public static int ILogB(QuadDouble x)
    {
        int e = Math.ILogB(x._x0);
        if (double.IsFinite(x._x0) && x._x0 != 0 && double.IsPow2(Math.Abs(x._x0)))
        {
            double t = x._x1 != 0 ? x._x1 : x._x2 != 0 ? x._x2 : x._x3;
            if (t != 0 && (t < 0) != (x._x0 < 0))
            {
                e--; // a power of two minus a positive tail
            }
        }
        return e;
    }

    // ---------------------------------------------------------------
    // Rounding to integers
    // ---------------------------------------------------------------

    /// <summary>The largest integer ≤ x (exact).</summary>
    public static QuadDouble Floor(QuadDouble a)
    {
        double x0 = Math.Floor(a._x0), x1 = 0, x2 = 0, x3 = 0;
        if (x0 == a._x0 && double.IsFinite(x0))
        {
            x1 = Math.Floor(a._x1);
            if (x1 == a._x1)
            {
                x2 = Math.Floor(a._x2);
                if (x2 == a._x2)
                {
                    x3 = Math.Floor(a._x3);
                }
            }
            Renormalize(ref x0, ref x1, ref x2, ref x3);
            if (x0 == 0)
            {
                return CopySign(Zero, a);
            }
        }
        return Make(x0, x1, x2, x3);
    }

    /// <summary>The smallest integer ≥ x (exact).</summary>
    public static QuadDouble Ceiling(QuadDouble a)
    {
        double x0 = Math.Ceiling(a._x0), x1 = 0, x2 = 0, x3 = 0;
        if (x0 == a._x0 && double.IsFinite(x0))
        {
            x1 = Math.Ceiling(a._x1);
            if (x1 == a._x1)
            {
                x2 = Math.Ceiling(a._x2);
                if (x2 == a._x2)
                {
                    x3 = Math.Ceiling(a._x3);
                }
            }
            Renormalize(ref x0, ref x1, ref x2, ref x3);
            if (x0 == 0)
            {
                return CopySign(Zero, a);
            }
        }
        return Make(x0, x1, x2, x3);
    }

    /// <summary>x rounded toward zero (exact).</summary>
    public static QuadDouble Truncate(QuadDouble x) => x._x0 >= 0 ? Floor(x) : Ceiling(x);

    /// <summary>x rounded to the nearest integer, ties to even.</summary>
    public static QuadDouble Round(QuadDouble x) => Round(x, MidpointRounding.ToEven);

    /// <summary>x rounded to an integer with the given midpoint rule.</summary>
    public static QuadDouble Round(QuadDouble x, MidpointRounding mode)
    {
        switch (mode)
        {
            case MidpointRounding.ToZero:
                return Truncate(x);
            case MidpointRounding.ToNegativeInfinity:
                return Floor(x);
            case MidpointRounding.ToPositiveInfinity:
                return Ceiling(x);
            case MidpointRounding.ToEven:
            case MidpointRounding.AwayFromZero:
                break;
            default:
                throw new ArgumentException("Invalid MidpointRounding value.", nameof(mode));
        }
        if (!double.IsFinite(x._x0) || IsInteger(x))
        {
            return x;
        }
        bool neg = x._x0 < 0;
        QuadDouble ax = Abs(x);
        QuadDouble f = Floor(ax);
        int c = CompareNonNaN(Subtract(ax, f), Single(0.5));
        bool up = c > 0 || (c == 0 && (mode == MidpointRounding.AwayFromZero || OddParity(f)));
        QuadDouble r = up ? Add(f, 1.0) : f;
        if (r._x0 == 0)
        {
            r = Zero;
        }
        return neg ? -r : r;
    }

    /// <summary>x rounded to <paramref name="digits"/> decimal places, ties to even.</summary>
    public static QuadDouble Round(QuadDouble x, int digits) => Round(x, digits, MidpointRounding.ToEven);

    /// <summary>x rounded to <paramref name="digits"/> decimal places (exact decimal rounding via <see cref="Real"/>, then to QD).</summary>
    public static QuadDouble Round(QuadDouble x, int digits, MidpointRounding mode)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(digits);
        if (digits == 0)
        {
            return Round(x, mode);
        }
        if (!double.IsFinite(x._x0) || IsInteger(x))
        {
            return x;
        }
        Real r = Real.Round(ToRealExact(x).WithPrecision(Math.Max(ToRealExact(x).Precision, 400)), digits, mode);
        QuadDouble q = FromReal(r);
        return q._x0 == 0 ? CopySign(Zero, x) : q;
    }

    // ---------------------------------------------------------------
    // IFloatingPoint layout members (delegated to the exact Real value)
    // ---------------------------------------------------------------

    private IFloatingPoint<Real> Layout => ToRealExact(this);

    int IFloatingPoint<QuadDouble>.GetExponentByteCount() => Layout.GetExponentByteCount();
    int IFloatingPoint<QuadDouble>.GetExponentShortestBitLength() => Layout.GetExponentShortestBitLength();
    int IFloatingPoint<QuadDouble>.GetSignificandBitLength() => Layout.GetSignificandBitLength();
    int IFloatingPoint<QuadDouble>.GetSignificandByteCount() => Layout.GetSignificandByteCount();

    bool IFloatingPoint<QuadDouble>.TryWriteExponentBigEndian(Span<byte> destination, out int bytesWritten) =>
        Layout.TryWriteExponentBigEndian(destination, out bytesWritten);

    bool IFloatingPoint<QuadDouble>.TryWriteExponentLittleEndian(Span<byte> destination, out int bytesWritten) =>
        Layout.TryWriteExponentLittleEndian(destination, out bytesWritten);

    bool IFloatingPoint<QuadDouble>.TryWriteSignificandBigEndian(Span<byte> destination, out int bytesWritten) =>
        Layout.TryWriteSignificandBigEndian(destination, out bytesWritten);

    bool IFloatingPoint<QuadDouble>.TryWriteSignificandLittleEndian(Span<byte> destination, out int bytesWritten) =>
        Layout.TryWriteSignificandLittleEndian(destination, out bytesWritten);
}
