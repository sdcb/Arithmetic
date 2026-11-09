using System.Runtime.CompilerServices;
using Sdcb.Arithmetic.Core;

namespace Sdcb.Arithmetic;

/// <summary>
/// A double-double floating-point number: the unevaluated sum <c>Hi + Lo</c> of two
/// <see cref="double"/> values, normalized so that <c>Hi == Hi + Lo</c> in double arithmetic,
/// giving about 106 significand bits (~32 decimal digits) with the exponent range of
/// <see cref="double"/>. Algorithms follow the QD library's dd_real (Hida, Li, Bailey) and the
/// FMA-based double-word algorithms of Joldes, Muller and Popescu.
/// <para>
/// This is a fast approximation, NOT a correctly rounded format: basic operations have a small
/// proven relative error bound (addition 3·2^-106, multiplication 2^-104, division about
/// 2^-102, square root about 2^-104) and transcendental functions are accurate to about 2^-100
/// over their usual ranges. There is no subnormal double-double range: below ~2^-969 the low
/// word underflows and precision degrades towards that of <see cref="double"/>. Use
/// <see cref="Real"/> when correct rounding or a wider exponent range is required.
/// </para>
/// </summary>
public readonly partial struct DoubleDouble :
    IEquatable<DoubleDouble>,
    IComparable<DoubleDouble>,
    IComparable
{
    private readonly double _hi;
    private readonly double _lo;

    /// <summary>Creates hi + lo, renormalizing (any two doubles are accepted).</summary>
    public DoubleDouble(double hi, double lo)
    {
        double s = ErrorFree.TwoSum(hi, lo, out double e);
        if (!double.IsFinite(s))
        {
            e = 0.0;
        }
        _hi = s;
        _lo = e;
    }

    /// <summary>Raw constructor: the caller guarantees (hi, lo) is already normalized.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private DoubleDouble(double hi, double lo, bool normalized)
    {
        _ = normalized;
        _hi = hi;
        _lo = lo;
    }

    /// <summary>The leading (high) word: the value rounded to the nearest double.</summary>
    public double Hi => _hi;

    /// <summary>The trailing (low) word: the rounding error of <see cref="Hi"/>, |Lo| &lt;= ulp(Hi)/2.</summary>
    public double Lo => _lo;

    /// <summary>Builds (s, e) from a Fast2Sum result, flushing a non-finite s to (s, 0).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static DoubleDouble Finish(double s, double e) =>
        double.IsFinite(s) ? new DoubleDouble(s, e, normalized: true) : new DoubleDouble(s, 0.0, normalized: true);

    // ---------------------------------------------------------------
    // Classification
    // ---------------------------------------------------------------

    public static bool IsNaN(DoubleDouble value) => double.IsNaN(value._hi);
    public static bool IsInfinity(DoubleDouble value) => double.IsInfinity(value._hi);
    public static bool IsPositiveInfinity(DoubleDouble value) => double.IsPositiveInfinity(value._hi);
    public static bool IsNegativeInfinity(DoubleDouble value) => double.IsNegativeInfinity(value._hi);
    public static bool IsFinite(DoubleDouble value) => double.IsFinite(value._hi);
    public static bool IsZero(DoubleDouble value) => value._hi == 0.0;

    /// <summary>True when the sign bit of <see cref="Hi"/> is set (including -0 and -∞).</summary>
    public static bool IsNegative(DoubleDouble value) => double.IsNegative(value._hi);

    /// <summary>True when the sign bit of <see cref="Hi"/> is clear (including +0 and +∞).</summary>
    public static bool IsPositive(DoubleDouble value) => double.IsPositive(value._hi);

    /// <summary>True when <see cref="Hi"/> is a normal double.</summary>
    public static bool IsNormal(DoubleDouble value) => double.IsNormal(value._hi);

    /// <summary>True when <see cref="Hi"/> is a subnormal double.</summary>
    public static bool IsSubnormal(DoubleDouble value) => double.IsSubnormal(value._hi);

    /// <summary>True for finite integral values.</summary>
    public static bool IsInteger(DoubleDouble value) => double.IsInteger(value._hi) && double.IsInteger(value._lo);

    /// <summary>True for even integers.</summary>
    public static bool IsEvenInteger(DoubleDouble value) =>
        IsInteger(value) && (value._lo == 0.0 ? double.IsEvenInteger(value._hi) : double.IsEvenInteger(value._lo));

    /// <summary>True for odd integers.</summary>
    public static bool IsOddInteger(DoubleDouble value) =>
        IsInteger(value) && (value._lo == 0.0 ? double.IsOddInteger(value._hi) : double.IsOddInteger(value._lo));

    /// <summary>True for positive powers of two.</summary>
    public static bool IsPow2(DoubleDouble value) => value._lo == 0.0 && double.IsPow2(value._hi);

    // ---------------------------------------------------------------
    // Equality and ordering (lexicographic on the normalized pair)
    // ---------------------------------------------------------------

    public static bool operator ==(DoubleDouble left, DoubleDouble right) => left._hi == right._hi && left._lo == right._lo;
    public static bool operator !=(DoubleDouble left, DoubleDouble right) => !(left == right);

    public static bool operator <(DoubleDouble left, DoubleDouble right) =>
        left._hi < right._hi || (left._hi == right._hi && left._lo < right._lo);

    public static bool operator >(DoubleDouble left, DoubleDouble right) =>
        left._hi > right._hi || (left._hi == right._hi && left._lo > right._lo);

    public static bool operator <=(DoubleDouble left, DoubleDouble right) =>
        left._hi < right._hi || (left._hi == right._hi && left._lo <= right._lo);

    public static bool operator >=(DoubleDouble left, DoubleDouble right) =>
        left._hi > right._hi || (left._hi == right._hi && left._lo >= right._lo);

    /// <summary>Value equality where NaN equals NaN and +0 equals -0 (like <see cref="double.Equals(double)"/>).</summary>
    public bool Equals(DoubleDouble other) => _hi.Equals(other._hi) && _lo.Equals(other._lo);

    public override bool Equals(object? obj) => obj is DoubleDouble other && Equals(other);

    public override int GetHashCode() => HashCode.Combine(_hi, _lo);

    /// <summary>Total order like <see cref="double.CompareTo(double)"/>: NaN sorts first, -0 equals +0.</summary>
    public int CompareTo(DoubleDouble other)
    {
        int c = _hi.CompareTo(other._hi);
        return c != 0 ? c : _lo.CompareTo(other._lo);
    }

    public int CompareTo(object? obj) => obj switch
    {
        null => 1,
        DoubleDouble other => CompareTo(other),
        _ => throw new ArgumentException("Object must be of type DoubleDouble.", nameof(obj)),
    };
}
