using System.Runtime.CompilerServices;
using Sdcb.Arithmetic.Core;
using static Sdcb.Arithmetic.Core.ErrorFree;

namespace Sdcb.Arithmetic;

public readonly partial struct QuadDouble
{
    // ---------------------------------------------------------------
    // Addition
    // ---------------------------------------------------------------

    /// <summary>
    /// QD's quick_three_accum: adds c into the double-length accumulator (a, b); returns the
    /// completed leading component when both accumulator words are nonzero, else 0.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double QuickThreeAccum(ref double a, ref double b, double c)
    {
        double s = TwoSum(b, c, out b);
        s = TwoSum(a, s, out a);
        bool za = a != 0.0, zb = b != 0.0;
        if (za && zb)
        {
            return s;
        }
        if (!zb)
        {
            b = a;
            a = s;
        }
        else
        {
            a = s;
        }
        return 0.0;
    }

    /// <summary>
    /// Accurate ("IEEE") addition from the QD library: merges the eight components by decreasing
    /// magnitude into a double-length accumulator, so the result has a small relative error even
    /// under massive cancellation.
    /// </summary>
    public static QuadDouble Add(QuadDouble a, QuadDouble b)
    {
        double naive = a._x0 + b._x0;
        if (!double.IsFinite(naive))
        {
            return Single(naive);
        }
        if (b._x0 == 0)
        {
            return a._x0 == 0 ? Single(naive) : a;
        }
        if (a._x0 == 0)
        {
            return b;
        }

        ref double pa = ref Unsafe.AsRef(in a._x0);
        ref double pb = ref Unsafe.AsRef(in b._x0);
        int i = 0, j = 0, k = 0;
        double u, v, t, s;
        double x0 = 0, x1 = 0, x2 = 0, x3 = 0;

        if (Math.Abs(pa) > Math.Abs(pb))
        {
            u = pa;
            i = 1;
        }
        else
        {
            u = pb;
            j = 1;
        }
        if (Math.Abs(Unsafe.Add(ref pa, i)) > Math.Abs(Unsafe.Add(ref pb, j)))
        {
            v = Unsafe.Add(ref pa, i++);
        }
        else
        {
            v = Unsafe.Add(ref pb, j++);
        }
        u = FastTwoSum(u, v, out v);

        while (k < 4)
        {
            if (i >= 4 && j >= 4)
            {
                SetX(k, u, ref x0, ref x1, ref x2, ref x3);
                if (k < 3)
                {
                    SetX(++k, v, ref x0, ref x1, ref x2, ref x3);
                }
                break;
            }
            if (i >= 4)
            {
                t = Unsafe.Add(ref pb, j++);
            }
            else if (j >= 4)
            {
                t = Unsafe.Add(ref pa, i++);
            }
            else if (Math.Abs(Unsafe.Add(ref pa, i)) > Math.Abs(Unsafe.Add(ref pb, j)))
            {
                t = Unsafe.Add(ref pa, i++);
            }
            else
            {
                t = Unsafe.Add(ref pb, j++);
            }

            s = QuickThreeAccum(ref u, ref v, t);
            if (s != 0.0)
            {
                SetX(k++, s, ref x0, ref x1, ref x2, ref x3);
            }
        }

        for (int m = i; m < 4; m++)
        {
            x3 += Unsafe.Add(ref pa, m);
        }
        for (int m = j; m < 4; m++)
        {
            x3 += Unsafe.Add(ref pb, m);
        }

        Renorm(ref x0, ref x1, ref x2, ref x3, 0.0);
        if (x0 == 0)
        {
            return Single(naive == 0 ? naive : 0.0);
        }
        return double.IsFinite(x0) ? new QuadDouble(x0, x1, x2, x3, Tag.Raw) : Single(naive);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void SetX(int k, double value, ref double x0, ref double x1, ref double x2, ref double x3)
    {
        switch (k)
        {
            case 0: x0 = value; break;
            case 1: x1 = value; break;
            case 2: x2 = value; break;
            default: x3 = value; break;
        }
    }

    /// <summary>Accurate subtraction a - b (see <see cref="Add(QuadDouble, QuadDouble)"/>).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static QuadDouble Subtract(QuadDouble a, QuadDouble b) => Add(a, -b);

    /// <summary>
    /// QD's sloppy addition: component-wise TwoSums then renormalization. Its error is bounded
    /// relative to |a| + |b| rather than to the result, so it is used internally only where no
    /// catastrophic cancellation can occur.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal static QuadDouble AddSloppy(QuadDouble a, QuadDouble b)
    {
        double s0 = TwoSum(a._x0, b._x0, out double t0);
        double s1 = TwoSum(a._x1, b._x1, out double t1);
        double s2 = TwoSum(a._x2, b._x2, out double t2);
        double s3 = TwoSum(a._x3, b._x3, out double t3);

        s1 = TwoSum(s1, t0, out t0);
        ThreeSum(ref s2, ref t0, ref t1);
        ThreeSum2(ref s3, ref t0, t2);
        t0 = t0 + t1 + t3;

        Renorm(ref s0, ref s1, ref s2, ref s3, t0);
        return double.IsFinite(s0) ? new QuadDouble(s0, s1, s2, s3, Tag.Raw) : Single(a._x0 + b._x0);
    }

    /// <summary>a + b for a double b (QD's accurate qd + double).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static QuadDouble Add(QuadDouble a, double b)
    {
        double c0 = TwoSum(a._x0, b, out double e);
        double c1 = TwoSum(a._x1, e, out e);
        double c2 = TwoSum(a._x2, e, out e);
        double c3 = TwoSum(a._x3, e, out e);
        Renorm(ref c0, ref c1, ref c2, ref c3, e);
        if (c0 == 0 || !double.IsFinite(c0))
        {
            return Single(a._x0 + b); // signed zero or the IEEE special value
        }
        return new QuadDouble(c0, c1, c2, c3, Tag.Raw);
    }

    /// <summary>Internal alias of <see cref="Add(QuadDouble, double)"/> used by the constructor.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static QuadDouble AddDouble(QuadDouble a, double b) => Add(a, b);

    // ---------------------------------------------------------------
    // Multiplication
    // ---------------------------------------------------------------

    /// <summary>a * b for a double b (QD's qd * double, accurate to about 2^-211).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static QuadDouble Multiply(QuadDouble a, double b)
    {
        double p0 = TwoProd(a._x0, b, out double q0);
        if (p0 == 0 || !double.IsFinite(p0))
        {
            return Single(p0);
        }
        double p1 = TwoProd(a._x1, b, out double q1);
        double p2 = TwoProd(a._x2, b, out double q2);
        double p3 = a._x3 * b;

        double s0 = p0;
        double s1 = TwoSum(q0, p1, out double s2);
        ThreeSum(ref s2, ref q1, ref p2);
        ThreeSum2(ref q1, ref q2, p3);
        double s3 = q1;
        double s4 = q2 + p2;

        Renorm(ref s0, ref s1, ref s2, ref s3, s4);
        return double.IsFinite(s0) ? new QuadDouble(s0, s1, s2, s3, Tag.Raw) : Single(p0);
    }

    /// <summary>
    /// QD's accurate multiplication: all partial products up to order ε³ are formed exactly with
    /// TwoProd and accumulated with error-free sums; order-ε⁴ terms are added in plain arithmetic.
    /// </summary>
    public static QuadDouble Multiply(QuadDouble a, QuadDouble b)
    {
        double p0 = TwoProd(a._x0, b._x0, out double q0);
        if (p0 == 0 || !double.IsFinite(p0))
        {
            return Single(p0);
        }
        double p1 = TwoProd(a._x0, b._x1, out double q1);
        double p2 = TwoProd(a._x1, b._x0, out double q2);

        ThreeSum(ref p1, ref p2, ref q0);

        // O(eps^2) terms.
        double p3 = TwoProd(a._x0, b._x2, out double q3);
        double p4 = TwoProd(a._x1, b._x1, out double q4);
        double p5 = TwoProd(a._x2, b._x0, out double q5);

        ThreeSum(ref p2, ref q1, ref q2);
        ThreeSum(ref p3, ref p4, ref p5);

        double s0 = TwoSum(p2, p3, out double t0);
        double s1 = TwoSum(q1, p4, out double t1);
        double s2 = q2 + p5;
        s1 = TwoSum(s1, t0, out t0);
        s2 += t0 + t1;

        // O(eps^3) terms.
        double p6 = TwoProd(a._x0, b._x3, out double q6);
        double p7 = TwoProd(a._x1, b._x2, out double q7);
        double p8 = TwoProd(a._x2, b._x1, out double q8);
        double p9 = TwoProd(a._x3, b._x0, out double q9);

        q0 = TwoSum(q0, q3, out q3);
        q4 = TwoSum(q4, q5, out q5);
        p6 = TwoSum(p6, p7, out p7);
        p8 = TwoSum(p8, p9, out p9);

        t0 = TwoSum(q0, q4, out t1);
        t1 += q3 + q5;

        double r0 = TwoSum(p6, p8, out double r1);
        r1 += p7 + p9;

        q3 = TwoSum(t0, r0, out q4);
        q4 += t1 + r1;

        t0 = TwoSum(q3, s1, out t1);
        t1 += q4;

        // O(eps^4) terms.
        t1 += a._x1 * b._x3 + a._x2 * b._x2 + a._x3 * b._x1 + q6 + q7 + q8 + q9 + s2;

        Renorm(ref p0, ref p1, ref s0, ref t0, t1);
        return double.IsFinite(p0) ? new QuadDouble(p0, p1, s0, t0, Tag.Raw) : Single(a._x0 * b._x0);
    }

    /// <summary>QD's sloppy multiplication (drops the order-ε³ TwoProd errors); used internally by series evaluation.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal static QuadDouble MulSloppy(QuadDouble a, QuadDouble b)
    {
        double p0 = TwoProd(a._x0, b._x0, out double q0);
        if (p0 == 0 || !double.IsFinite(p0))
        {
            return Single(p0);
        }
        double p1 = TwoProd(a._x0, b._x1, out double q1);
        double p2 = TwoProd(a._x1, b._x0, out double q2);
        double p3 = TwoProd(a._x0, b._x2, out double q3);
        double p4 = TwoProd(a._x1, b._x1, out double q4);
        double p5 = TwoProd(a._x2, b._x0, out double q5);

        ThreeSum(ref p1, ref p2, ref q0);
        ThreeSum(ref p2, ref q1, ref q2);
        ThreeSum(ref p3, ref p4, ref p5);

        double s0 = TwoSum(p2, p3, out double t0);
        double s1 = TwoSum(q1, p4, out double t1);
        double s2 = q2 + p5;
        s1 = TwoSum(s1, t0, out t0);
        s2 += t0 + t1;

        s1 += a._x0 * b._x3 + a._x1 * b._x2 + a._x2 * b._x1 + a._x3 * b._x0 + q0 + q3 + q4 + q5;
        Renorm(ref p0, ref p1, ref s0, ref s1, s2);
        return double.IsFinite(p0) ? new QuadDouble(p0, p1, s0, s1, Tag.Raw) : Single(a._x0 * b._x0);
    }

    /// <summary>a² (QD's sqr; cheaper than <see cref="Multiply(QuadDouble, QuadDouble)"/>).</summary>
    public static QuadDouble Square(QuadDouble a)
    {
        double p0 = TwoSquare(a._x0, out double q0);
        if (p0 == 0 || !double.IsFinite(p0))
        {
            return Single(p0);
        }
        double p1 = TwoProd(2.0 * a._x0, a._x1, out double q1);
        double p2 = TwoProd(2.0 * a._x0, a._x2, out double q2);
        double p3 = TwoSquare(a._x1, out double q3);

        p1 = TwoSum(q0, p1, out q0);
        q0 = TwoSum(q0, q1, out q1);
        p2 = TwoSum(p2, p3, out p3);

        double s0 = TwoSum(q0, p2, out double t0);
        double s1 = TwoSum(q1, p3, out double t1);
        s1 = TwoSum(s1, t0, out t0);
        t0 += t1;

        s1 = FastTwoSum(s1, t0, out t0);
        p2 = FastTwoSum(s0, s1, out t1);
        p3 = FastTwoSum(t1, t0, out q0);

        double p4 = 2.0 * a._x0 * a._x3;
        double p5 = 2.0 * a._x1 * a._x2;

        p4 = TwoSum(p4, p5, out p5);
        q2 = TwoSum(q2, q3, out q3);

        t0 = TwoSum(p4, q2, out t1);
        t1 = t1 + p5 + q3;

        p3 = TwoSum(p3, t0, out p4);
        p4 = p4 + q0 + t1;

        Renorm(ref p0, ref p1, ref p2, ref p3, p4);
        return double.IsFinite(p0) ? new QuadDouble(p0, p1, p2, p3, Tag.Raw) : Single(a._x0 * a._x0);
    }

    /// <summary>Exact multiplication by a power of two (each component scaled; exact barring under/overflow).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static QuadDouble MulPwr2(QuadDouble a, double p2) => Make(a._x0 * p2, a._x1 * p2, a._x2 * p2, a._x3 * p2);

    // ---------------------------------------------------------------
    // Division
    // ---------------------------------------------------------------

    /// <summary>
    /// Long division (QD's accurate_div): five quotient digits q_k = r_0 / b_0 with the remainder
    /// r ← r - b·q_k updated in QD arithmetic, then renormalized.
    /// </summary>
    public static QuadDouble Divide(QuadDouble a, QuadDouble b)
    {
        double q0 = a._x0 / b._x0;
        if (q0 == 0 || !double.IsFinite(q0) || !double.IsFinite(b._x0))
        {
            return Single(q0);
        }
        double inv = 1.0 / b._x0;
        QuadDouble r = AddSloppy(a, -Multiply(b, q0));
        double q1 = r._x0 * inv;
        r = AddSloppy(r, -Multiply(b, q1));
        double q2 = r._x0 * inv;
        r = AddSloppy(r, -Multiply(b, q2));
        double q3 = r._x0 * inv;
        // The fifth digit only needs a few correct bits: estimate the leading component of r - b·q3
        // from r0 - b0·q3 (exact by Sterbenz) plus the next-order terms.
        double p = TwoProd(b._x0, q3, out double e);
        double q4 = (((r._x0 - p) - e) + (r._x1 - b._x1 * q3)) * inv;

        Renorm(ref q0, ref q1, ref q2, ref q3, q4);
        return Make(q0, q1, q2, q3);
    }

    /// <summary>a / b for a double b (long division with the exact products q_k·b).</summary>
    public static QuadDouble Divide(QuadDouble a, double b)
    {
        double q0 = a._x0 / b;
        if (q0 == 0 || !double.IsFinite(q0) || !double.IsFinite(b))
        {
            return Single(q0);
        }
        double inv = 1.0 / b;
        QuadDouble r = SubProduct(a, q0, b);
        double q1 = r._x0 * inv;
        r = SubProduct(r, q1, b);
        double q2 = r._x0 * inv;
        r = SubProduct(r, q2, b);
        double q3 = r._x0 * inv;
        double p = TwoProd(q3, b, out double e);
        double q4 = (((r._x0 - p) - e) + r._x1) * inv;

        Renorm(ref q0, ref q1, ref q2, ref q3, q4);
        return Make(q0, q1, q2, q3);
    }

    /// <summary>r - q·b with the product formed exactly (TwoProd), accumulated in sloppy QD.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static QuadDouble SubProduct(QuadDouble r, double q, double b)
    {
        double p = TwoProd(q, b, out double e);
        return AddSloppy(r, new QuadDouble(-p, -e, 0, 0, Tag.Raw));
    }

    /// <summary>a / b for a double numerator.</summary>
    public static QuadDouble Divide(double a, QuadDouble b) => Divide(Single(a), b);

    /// <summary>1 / a.</summary>
    public static QuadDouble Reciprocal(QuadDouble a) => Divide(One, a);

    // ---------------------------------------------------------------
    // Operators
    // ---------------------------------------------------------------

    /// <summary>Negation (exact).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static QuadDouble operator -(QuadDouble value) => new(-value._x0, -value._x1, -value._x2, -value._x3, Tag.Raw);

    public static QuadDouble operator +(QuadDouble value) => value;

    public static QuadDouble operator +(QuadDouble left, QuadDouble right) => Add(left, right);
    public static QuadDouble operator +(QuadDouble left, double right) => Add(left, right);
    public static QuadDouble operator +(double left, QuadDouble right) => Add(right, left);

    public static QuadDouble operator -(QuadDouble left, QuadDouble right) => Add(left, -right);
    public static QuadDouble operator -(QuadDouble left, double right) => Add(left, -right);
    public static QuadDouble operator -(double left, QuadDouble right) => Add(-right, left);

    public static QuadDouble operator *(QuadDouble left, QuadDouble right) => Multiply(left, right);
    public static QuadDouble operator *(QuadDouble left, double right) => Multiply(left, right);
    public static QuadDouble operator *(double left, QuadDouble right) => Multiply(right, left);

    public static QuadDouble operator /(QuadDouble left, QuadDouble right) => Divide(left, right);
    public static QuadDouble operator /(QuadDouble left, double right) => Divide(left, right);
    public static QuadDouble operator /(double left, QuadDouble right) => Divide(left, right);

    /// <summary>Truncated remainder with the sign of the dividend (like C# <see cref="double"/> %).</summary>
    public static QuadDouble operator %(QuadDouble left, QuadDouble right) => TruncatedRemainder(left, right);

    public static QuadDouble operator ++(QuadDouble value) => Add(value, 1.0);
    public static QuadDouble operator --(QuadDouble value) => Add(value, -1.0);

    // ---------------------------------------------------------------
    // Remainder, FMA
    // ---------------------------------------------------------------

    /// <summary>
    /// x - trunc(x / y)·y with the sign of x (like <see cref="Math.IEEERemainder"/>'s truncating sibling,
    /// C# <see cref="double"/> %). The exact remainder is computed through <see cref="Real"/> and then
    /// rounded to QD, so the result is accurate even for huge quotients.
    /// </summary>
    public static QuadDouble TruncatedRemainder(QuadDouble x, QuadDouble y)
    {
        if (double.IsNaN(x._x0) || double.IsNaN(y._x0) || double.IsInfinity(x._x0) || y._x0 == 0)
        {
            return NaN;
        }
        if (double.IsInfinity(y._x0) || x._x0 == 0)
        {
            return x;
        }
        if (CompareNonNaN(Abs(x), Abs(y)) < 0)
        {
            return x;
        }
        Real rx = ToRealExact(x), ry = ToRealExact(y);
        int p = Math.Max(rx.Precision, ry.Precision);
        QuadDouble r = FromReal(Real.TruncatedRemainder(rx, ry, p));
        return r._x0 == 0 ? CopySign(Zero, x) : r;
    }

    /// <summary>
    /// (left · right) + addend. The product is formed in QD arithmetic; when the sum cancels by
    /// more than 8 bits it is recomputed through <see cref="Real"/> (correctly rounded at 320 bits),
    /// so the result is accurate, though not correctly rounded, in every case.
    /// </summary>
    public static QuadDouble FusedMultiplyAdd(QuadDouble left, QuadDouble right, QuadDouble addend)
    {
        QuadDouble p = Multiply(left, right);
        if (!double.IsFinite(p._x0) || !double.IsFinite(addend._x0) || p._x0 == 0 || addend._x0 == 0)
        {
            return double.IsFinite(p._x0) ? Add(p, addend) : Single(Math.FusedMultiplyAdd(left._x0, right._x0, addend._x0));
        }
        QuadDouble s = Add(p, addend);
        if (Math.Abs(s._x0) >= Math.Max(Math.Abs(p._x0), Math.Abs(addend._x0)) * 0.00390625 /* 2^-8 */)
        {
            return s;
        }
        Real exact = Real.FusedMultiplyAdd(ToRealExact(left), ToRealExact(right), ToRealExact(addend), 320);
        return FromReal(exact);
    }

    // ---------------------------------------------------------------
    // Square root
    // ---------------------------------------------------------------

    /// <summary>
    /// √a: a double-double estimate (one Newton step from <see cref="Math.Sqrt"/>) refined by one
    /// Newton step whose residual a - y² is formed in QD arithmetic. -0 gives -0, negative values NaN.
    /// </summary>
    public static QuadDouble Sqrt(QuadDouble a)
    {
        double a0 = a._x0;
        if (!(a0 > 0) || double.IsPositiveInfinity(a0))
        {
            return Single(a0 == 0 ? a0 : Math.Sqrt(a0));
        }

        // Keep y² and the residual comfortably inside the normal range.
        int scale = 0;
        if (a0 > 8.452712498170644E+270 /* 2^900 */ || a0 < 1.499696813895631E-241 /* 2^-800 */)
        {
            scale = Math.ILogB(a0) / 2;
            a = ScaleB(a, -2 * scale);
            a0 = a._x0;
        }

        double s = Math.Sqrt(a0);
        double e = Math.FusedMultiplyAdd(-s, s, a0);
        double inv2s = 0.5 / s;
        double t = (e + a._x1) * inv2s;
        double yh = FastTwoSum(s, t, out double yl);

        QuadDouble y2 = Square(new QuadDouble(yh, yl, 0, 0, Tag.Raw));
        QuadDouble r = Subtract(a, y2);

        // c = r / (2 y) in double-double precision.
        double d0 = 2.0 * yh, d1 = 2.0 * yl;
        double ch = r._x0 / d0;
        double ph = TwoProd(ch, d0, out double pe);
        double rem = ((r._x0 - ph) - pe) + r._x1 - ch * d1;
        double cl = rem / d0;

        QuadDouble y = new(yh, yl, 0, 0, Tag.Raw);
        // √a = y + δ - δ²/(2y) + O(δ³) with δ = r/(2y) = ch + cl.
        double c2 = -ch * ch * (0.5 / yh);
        y = Add(Add(Add(y, ch), cl), c2);
        return scale == 0 ? y : ScaleB(y, scale);
    }
}
