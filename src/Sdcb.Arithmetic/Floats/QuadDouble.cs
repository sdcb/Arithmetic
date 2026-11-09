using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using static Sdcb.Arithmetic.Core.ErrorFree;

namespace Sdcb.Arithmetic;

/// <summary>
/// A quad-double floating-point number: the unevaluated sum of four <see cref="double"/>
/// components <c>X0 + X1 + X2 + X3</c> that do not overlap and decrease in magnitude, giving about
/// 212 significand bits (≈ 64 decimal digits) with the exponent range of <see cref="double"/>.
/// </summary>
/// <remarks>
/// <para>
/// The algorithms follow the QD library of Hida, Li and Bailey (<c>qd_real</c>). Results are fast
/// approximations, <b>not</b> correctly rounded: addition and subtraction use QD's accurate
/// ("IEEE") addition when the operands' signs differ, so cancellation is handled well, and QD's
/// sloppy addition when they agree (no cancellation is possible then); multiplication uses QD's
/// accurate multiplication; division is long division with five quotient digits, and the square
/// root the analogous digit-by-digit recurrence, both with remainders updated by error-free
/// transforms. Measured worst relative errors on random operands (the tests compare against
/// <see cref="Real"/> at 500 bits): about 2^-213 for + and -, about 2^-214 for ×, about 2^-215 for ÷
/// and √, and between 2^-208.5 and 2^-215 for
/// the transcendental functions in normal ranges (Pow's error grows like |y·ln x|·2^-212).
/// </para>
/// <para>
/// The exponent range is that of <see cref="double"/>. Values whose leading component is close to
/// <see cref="double.MaxValue"/> may overflow in intermediate steps, and values below about 2^-860
/// lose precision because their lower components become subnormal or underflow to zero.
/// A non-finite value is stored as a single component (the others are zero).
/// </para>
/// <para>
/// Components are exposed through <see cref="X0"/>, <see cref="X1"/>, <see cref="X2"/> and <see cref="X3"/>;
/// there is no indexer.
/// </para>
/// </remarks>
[StructLayout(LayoutKind.Sequential)]
public readonly partial struct QuadDouble
{
    private readonly double _x0;
    private readonly double _x1;
    private readonly double _x2;
    private readonly double _x3;

    /// <summary>
    /// Creates the quad-double nearest (in QD arithmetic) to the exact sum
    /// <paramref name="x0"/> + <paramref name="x1"/> + <paramref name="x2"/> + <paramref name="x3"/>.
    /// The components may overlap or come in any order; they are summed and renormalized.
    /// </summary>
    public QuadDouble(double x0, double x1 = 0, double x2 = 0, double x3 = 0)
    {
        if (x1 == 0 && x2 == 0 && x3 == 0)
        {
            _x0 = x0;
            _x1 = 0;
            _x2 = 0;
            _x3 = 0;
            return;
        }
        QuadDouble s = AddDouble(AddDouble(AddDouble(new QuadDouble(x0, 0, 0, 0, Tag.Raw), x1), x2), x3);
        this = s;
    }

    private enum Tag : byte
    {
        Raw,
    }

    /// <summary>Stores the components as given (callers guarantee they are normalized).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private QuadDouble(double x0, double x1, double x2, double x3, Tag _)
    {
        _x0 = x0;
        _x1 = x1;
        _x2 = x2;
        _x3 = x3;
    }

    /// <summary>Builds a value from normalized components; a non-finite leading component clears the others.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static QuadDouble Make(double x0, double x1, double x2, double x3) =>
        double.IsFinite(x0) ? new QuadDouble(x0, x1, x2, x3, Tag.Raw) : new QuadDouble(x0, 0, 0, 0, Tag.Raw);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static QuadDouble Single(double x0) => new(x0, 0, 0, 0, Tag.Raw);

    /// <summary>The leading (largest) component, which is also an approximation of the value as a double.</summary>
    public double X0 => _x0;

    /// <summary>The second component, |X1| ≤ ulp(X0).</summary>
    public double X1 => _x1;

    /// <summary>The third component, |X2| ≤ ulp(X1).</summary>
    public double X2 => _x2;

    /// <summary>The fourth (smallest) component, |X3| ≤ ulp(X2).</summary>
    public double X3 => _x3;

    // ---------------------------------------------------------------
    // Constants (correctly rounded values split into four doubles)
    // ---------------------------------------------------------------

    /// <summary>0.</summary>
    public static QuadDouble Zero => default;

    /// <summary>-0.</summary>
    public static QuadDouble NegativeZero => Single(-0.0);

    /// <summary>1.</summary>
    public static QuadDouble One => Single(1.0);

    /// <summary>-1.</summary>
    public static QuadDouble NegativeOne => Single(-1.0);

    /// <summary>Not a number.</summary>
    public static QuadDouble NaN => Single(double.NaN);

    /// <summary>+∞.</summary>
    public static QuadDouble PositiveInfinity => Single(double.PositiveInfinity);

    /// <summary>-∞.</summary>
    public static QuadDouble NegativeInfinity => Single(double.NegativeInfinity);

    /// <summary>π, the correctly rounded 212+ bit value split into four doubles.</summary>
    public static QuadDouble Pi => new(3.141592653589793, 1.2246467991473532E-16, -2.9947698097183397E-33, 1.1124542208633653E-49, Tag.Raw);

    /// <summary>2π.</summary>
    public static QuadDouble Tau => new(6.283185307179586, 2.4492935982947064E-16, -5.989539619436679E-33, 2.2249084417267306E-49, Tag.Raw);

    /// <summary>Euler's number e.</summary>
    public static QuadDouble E => new(2.718281828459045, 1.4456468917292502E-16, -2.1277171080381768E-33, 1.5156301598412191E-49, Tag.Raw);

    /// <summary>ln 2.</summary>
    public static QuadDouble Ln2 => new(0.6931471805599453, 2.3190468138462996E-17, 5.707708438416212E-34, -3.5824322106018114E-50, Tag.Raw);

    /// <summary>ln 10.</summary>
    public static QuadDouble Ln10 => new(2.302585092994046, -2.1707562233822494E-16, -9.984262454465777E-33, -4.023357454450206E-49, Tag.Raw);

    /// <summary>
    /// The smallest positive value, 2^-1074 (same as <see cref="double.Epsilon"/>; this follows the .NET
    /// convention, not QD's machine epsilon 2^-209).
    /// </summary>
    public static QuadDouble Epsilon => Single(double.Epsilon);

    /// <summary>The relative precision 2^-209, the same convention as the QD library's qd_real::_eps.</summary>
    public static QuadDouble MachineEpsilon => Single(Math.ScaleB(1.0, -209));

    /// <summary>
    /// The largest finite value (QD's <c>_max</c>): <see cref="double.MaxValue"/> followed by the largest
    /// non-overlapping tails, (2^53 - 1)·(2^971 + 2^917 + 2^863 + 2^809).
    /// </summary>
    public static QuadDouble MaxValue => new(1.7976931348623157E+308, 9.979201547673598E+291, 5.5395696628011126E+275, 3.075078893078405E+259, Tag.Raw);

    /// <summary>The most negative finite value, -<see cref="MaxValue"/>.</summary>
    public static QuadDouble MinValue => new(-1.7976931348623157E+308, -9.979201547673598E+291, -5.5395696628011126E+275, -3.075078893078405E+259, Tag.Raw);

    internal static QuadDouble PiOver2 => new(1.5707963267948966, 6.123233995736766E-17, -1.4973849048591698E-33, 5.562271104316826E-50, Tag.Raw);
    internal static QuadDouble PiOver4 => new(0.7853981633974483, 3.061616997868383E-17, -7.486924524295849E-34, 2.781135552158413E-50, Tag.Raw);
    internal static QuadDouble InvPi => new(0.3183098861837907, -1.9678676675182486E-17, -1.0721436282893004E-33, 8.053563926594112E-50, Tag.Raw);
    internal static QuadDouble Log2E => new(1.4426950408889634, 2.0355273740931033E-17, -1.0614659956117258E-33, -1.3836716780181402E-50, Tag.Raw);
    internal static QuadDouble Log10E => new(0.4342944819032518, 1.098319650216765E-17, 3.717181233110959E-34, 7.734484346504299E-51, Tag.Raw);
    internal static QuadDouble PiOver180 => new(0.017453292519943295, 2.9486522708701687E-19, -1.3427726813345382E-35, 1.4287195201881262E-52, Tag.Raw);
    internal static QuadDouble D180OverPi => new(57.29577951308232, -1.9878495670576283E-15, -1.6833394980391744E-31, -5.5659577936984464E-49, Tag.Raw);

    // ---------------------------------------------------------------
    // Classification
    // ---------------------------------------------------------------

    /// <summary>True for NaN.</summary>
    public static bool IsNaN(QuadDouble value) => double.IsNaN(value._x0);

    /// <summary>True for ±∞.</summary>
    public static bool IsInfinity(QuadDouble value) => double.IsInfinity(value._x0);

    /// <summary>True for +∞.</summary>
    public static bool IsPositiveInfinity(QuadDouble value) => double.IsPositiveInfinity(value._x0);

    /// <summary>True for -∞.</summary>
    public static bool IsNegativeInfinity(QuadDouble value) => double.IsNegativeInfinity(value._x0);

    /// <summary>True for zero and finite nonzero values.</summary>
    public static bool IsFinite(QuadDouble value) => double.IsFinite(value._x0);

    /// <summary>True for ±0.</summary>
    public static bool IsZero(QuadDouble value) => value._x0 == 0;

    /// <summary>True when the sign bit is set (including -0 and -∞; NaN follows its sign bit like <see cref="double.IsNegative"/>).</summary>
    public static bool IsNegative(QuadDouble value) => double.IsNegative(value._x0);

    /// <summary>True when the sign bit is clear.</summary>
    public static bool IsPositive(QuadDouble value) => double.IsPositive(value._x0);

    /// <summary>True when the leading component is a normal double.</summary>
    public static bool IsNormal(QuadDouble value) => double.IsNormal(value._x0);

    /// <summary>True when the leading component is a subnormal double.</summary>
    public static bool IsSubnormal(QuadDouble value) => double.IsSubnormal(value._x0);

    /// <summary>True for finite integral values (every component is an integer).</summary>
    public static bool IsInteger(QuadDouble value) =>
        double.IsInteger(value._x0) && double.IsInteger(value._x1) && double.IsInteger(value._x2) && double.IsInteger(value._x3);

    /// <summary>True for even integers (including zero).</summary>
    public static bool IsEvenInteger(QuadDouble value) => IsInteger(value) && !OddParity(value);

    /// <summary>True for odd integers.</summary>
    public static bool IsOddInteger(QuadDouble value) => IsInteger(value) && OddParity(value);

    /// <summary>Parity of an integral value: the XOR of its components' parities (components ≥ 2^53 are even).</summary>
    private static bool OddParity(QuadDouble v) => IsOdd(v._x0) ^ IsOdd(v._x1) ^ IsOdd(v._x2) ^ IsOdd(v._x3);

    private static bool IsOdd(double c) => Math.Abs(c) < 9007199254740992.0 && Math.Abs(c % 2.0) == 1.0;

    /// <summary>True for positive powers of two.</summary>
    public static bool IsPow2(QuadDouble value) =>
        value._x0 > 0 && double.IsFinite(value._x0) && value._x1 == 0 && double.IsPow2(value._x0);

    // ---------------------------------------------------------------
    // Equality, ordering and hashing
    // ---------------------------------------------------------------

    /// <summary>
    /// Three-way comparison of two non-NaN values. Normalized representations are almost always
    /// unique, so equal components decide quickly; otherwise the sign of the exact-sign difference decides.
    /// </summary>
    private static int CompareNonNaN(QuadDouble a, QuadDouble b)
    {
        if (a._x0 != b._x0)
        {
            double d = a._x0 - b._x0;
            // Tails are below ulp(x0), so a gap wider than a few ulps of the larger leading component decides.
            double m = Math.Max(Math.Abs(a._x0), Math.Abs(b._x0));
            if (!double.IsFinite(m) || Math.Abs(d) > m * 8.881784197001252E-16 /* 2^-50 */)
            {
                return d > 0 ? 1 : -1;
            }
        }
        else if (a._x1 == b._x1 && a._x2 == b._x2 && a._x3 == b._x3)
        {
            return 0;
        }
        double s = Subtract(a, b)._x0;
        return s > 0 ? 1 : s < 0 ? -1 : 0;
    }

    /// <summary>Value equality (IEEE: NaN != NaN, +0 == -0).</summary>
    public static bool operator ==(QuadDouble left, QuadDouble right) =>
        !double.IsNaN(left._x0) && !double.IsNaN(right._x0) && CompareNonNaN(left, right) == 0;

    /// <summary>Value inequality (IEEE: NaN != anything).</summary>
    public static bool operator !=(QuadDouble left, QuadDouble right) => !(left == right);

    public static bool operator <(QuadDouble left, QuadDouble right) =>
        !double.IsNaN(left._x0) && !double.IsNaN(right._x0) && CompareNonNaN(left, right) < 0;

    public static bool operator >(QuadDouble left, QuadDouble right) =>
        !double.IsNaN(left._x0) && !double.IsNaN(right._x0) && CompareNonNaN(left, right) > 0;

    public static bool operator <=(QuadDouble left, QuadDouble right) =>
        !double.IsNaN(left._x0) && !double.IsNaN(right._x0) && CompareNonNaN(left, right) <= 0;

    public static bool operator >=(QuadDouble left, QuadDouble right) =>
        !double.IsNaN(left._x0) && !double.IsNaN(right._x0) && CompareNonNaN(left, right) >= 0;

    /// <summary>Equality like <see cref="double.Equals(double)"/>: NaN equals NaN and +0 equals -0.</summary>
    public bool Equals(QuadDouble other)
    {
        bool n1 = double.IsNaN(_x0), n2 = double.IsNaN(other._x0);
        return n1 || n2 ? n1 && n2 : CompareNonNaN(this, other) == 0;
    }

    public override bool Equals(object? obj) => obj is QuadDouble q && Equals(q);

    /// <summary>Hash of the value rounded to the nearest double, so equal values (including ±0) hash alike.</summary>
    public override int GetHashCode()
    {
        if (double.IsNaN(_x0))
        {
            return double.NaN.GetHashCode();
        }
        double d = ToDoubleNearest(this);
        return d == 0 ? 0 : d.GetHashCode();
    }

    /// <summary>Total order like <see cref="double.CompareTo(double)"/>: NaN is less than every other value and equal to itself.</summary>
    public int CompareTo(QuadDouble other)
    {
        bool n1 = double.IsNaN(_x0), n2 = double.IsNaN(other._x0);
        if (n1 || n2)
        {
            return n1 && n2 ? 0 : n1 ? -1 : 1;
        }
        return CompareNonNaN(this, other);
    }

    public int CompareTo(object? obj) => obj switch
    {
        null => 1,
        QuadDouble q => CompareTo(q),
        _ => throw new ArgumentException("Object must be of type QuadDouble.", nameof(obj)),
    };

    /// <summary>
    /// The double nearest to the exact value (ties to even), using that the tail X1 + X2 + X3 is
    /// below ulp(X0) so only the neighbours of X0 are candidates.
    /// </summary>
    internal static double ToDoubleNearest(QuadDouble a)
    {
        double s = a._x0, r = a._x1;
        if (r == 0 || !double.IsFinite(s))
        {
            return s;
        }
        // The tail's sign is that of X1; the candidate is the neighbour of s in that direction.
        double n = r > 0 ? Math.BitIncrement(s) : Math.BitDecrement(s);
        double half = double.IsInfinity(n) ? 9.9792015476736E+291 /* 2^970 */ : Math.Abs(n - s) * 0.5;
        if (half == 0)
        {
            return s;
        }
        double d = Math.Abs(r) - half; // exact when |r| is within a factor two of half; otherwise its sign is right
        int c;
        if (d != 0)
        {
            c = d > 0 ? 1 : -1;
        }
        else
        {
            double t = a._x2 != 0 ? a._x2 : a._x3;
            c = t == 0 ? 0 : (t > 0) == (r > 0) ? 1 : -1;
        }
        if (c < 0)
        {
            return s;
        }
        if (c > 0)
        {
            return n;
        }
        // Exact tie: the neighbour with an even significand.
        return (BitConverter.DoubleToInt64Bits(s) & 1) == 0 ? s : n;
    }

    // ---------------------------------------------------------------
    // Kernels shared across partials
    // ---------------------------------------------------------------

    /// <summary>The components as a read-only span (no copy).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static ReadOnlySpan<double> Components(in QuadDouble a) =>
        MemoryMarshal.CreateReadOnlySpan(ref Unsafe.AsRef(in a._x0), 4);

    /// <summary>
    /// QD's five-input renorm, identical to <see cref="Core.ErrorFree.Renormalize(ref double, ref double, ref double, ref double, ref double)"/>
    /// but force-inlined so the components stay in registers on the hot paths. Callers handle
    /// non-finite leading components themselves.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void Renorm(ref double c0, ref double c1, ref double c2, ref double c3, double c4)
    {
        double s0 = FastTwoSum(c3, c4, out c4);
        s0 = FastTwoSum(c2, s0, out c3);
        s0 = FastTwoSum(c1, s0, out c2);
        c0 = FastTwoSum(c0, s0, out c1);

        double s1, s2 = 0.0, s3 = 0.0;
        s0 = FastTwoSum(c0, c1, out s1);
        if (s1 != 0.0)
        {
            s1 = FastTwoSum(s1, c2, out s2);
            if (s2 != 0.0)
            {
                s2 = FastTwoSum(s2, c3, out s3);
                if (s3 != 0.0)
                {
                    s3 += c4;
                }
                else
                {
                    s2 = FastTwoSum(s2, c4, out s3);
                }
            }
            else
            {
                s1 = FastTwoSum(s1, c3, out s2);
                if (s2 != 0.0)
                {
                    s2 = FastTwoSum(s2, c4, out s3);
                }
                else
                {
                    s1 = FastTwoSum(s1, c4, out s2);
                }
            }
        }
        else
        {
            s0 = FastTwoSum(s0, c2, out s1);
            if (s1 != 0.0)
            {
                s1 = FastTwoSum(s1, c3, out s2);
                if (s2 != 0.0)
                {
                    s2 = FastTwoSum(s2, c4, out s3);
                }
                else
                {
                    s1 = FastTwoSum(s1, c4, out s2);
                }
            }
            else
            {
                s0 = FastTwoSum(s0, c3, out s1);
                if (s1 != 0.0)
                {
                    s1 = FastTwoSum(s1, c4, out s2);
                }
                else
                {
                    s0 = FastTwoSum(s0, c4, out s1);
                }
            }
        }

        c0 = s0;
        c1 = s1;
        c2 = s2;
        c3 = s3;
    }

    /// <summary>QD three_sum: (a, b, c) → (a + b + c, error, error of error), exact.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void ThreeSum(ref double a, ref double b, ref double c)
    {
        double t1 = TwoSum(a, b, out double t2);
        a = TwoSum(c, t1, out double t3);
        b = TwoSum(t2, t3, out c);
    }

    /// <summary>QD three_sum2: like <see cref="ThreeSum"/> but folds the two error terms into b.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void ThreeSum2(ref double a, ref double b, double c)
    {
        double t1 = TwoSum(a, b, out double t2);
        a = TwoSum(c, t1, out double t3);
        b = t2 + t3;
    }
}
