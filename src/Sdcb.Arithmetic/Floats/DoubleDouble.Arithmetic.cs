using System.Runtime.CompilerServices;
using Sdcb.Arithmetic.Core;

namespace Sdcb.Arithmetic;

public readonly partial struct DoubleDouble
{
    // Error bounds below are relative to the exact result, u = 2^-53, from Joldes, Muller and
    // Popescu, "Tight and rigorous error bounds for basic building blocks of double-word
    // arithmetic" (ACM TOMS 2017).

    // ---------------------------------------------------------------
    // Sign
    // ---------------------------------------------------------------

    public static DoubleDouble operator -(DoubleDouble value) => new(-value._hi, -value._lo, normalized: true);
    public static DoubleDouble operator +(DoubleDouble value) => value;

    public static DoubleDouble Abs(DoubleDouble value) => double.IsNegative(value._hi) ? -value : value;

    /// <summary>The magnitude of value with the sign of sign.</summary>
    public static DoubleDouble CopySign(DoubleDouble value, DoubleDouble sign) =>
        double.IsNegative(value._hi) != double.IsNegative(sign._hi) ? -value : value;

    /// <summary>-1, 0 or +1; throws <see cref="ArithmeticException"/> for NaN (like <see cref="Math.Sign(double)"/>).</summary>
    public static int Sign(DoubleDouble value) => Math.Sign(value._hi);

    // ---------------------------------------------------------------
    // Addition (AccurateDWPlusDW, 3u^2; DWPlusFP, 2u^2)
    // ---------------------------------------------------------------

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static DoubleDouble operator +(DoubleDouble a, DoubleDouble b)
    {
        double s = ErrorFree.TwoSum(a._hi, b._hi, out double se);
        double t = ErrorFree.TwoSum(a._lo, b._lo, out double te);
        se += t;
        s = ErrorFree.FastTwoSum(s, se, out se);
        se += te;
        s = ErrorFree.FastTwoSum(s, se, out se);
        return double.IsFinite(s) ? new DoubleDouble(s, se, normalized: true) : Overflowed(a._hi + b._hi);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static DoubleDouble operator -(DoubleDouble a, DoubleDouble b) => a + (-b);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static DoubleDouble operator +(DoubleDouble a, double b)
    {
        double s = ErrorFree.TwoSum(a._hi, b, out double e);
        e += a._lo;
        s = ErrorFree.FastTwoSum(s, e, out e);
        return double.IsFinite(s) ? new DoubleDouble(s, e, normalized: true) : Overflowed(a._hi + b);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static DoubleDouble operator +(double a, DoubleDouble b) => b + a;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static DoubleDouble operator -(DoubleDouble a, double b) => a + (-b);

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static DoubleDouble operator -(double a, DoubleDouble b) => (-b) + a;

    /// <summary>
    /// The result of an operation whose rounded high word is not finite: <paramref name="h"/>
    /// (the plain double result) when that is ±∞ or NaN, otherwise ±∞ (an overflow that only
    /// appeared once the low words were added); the error term is never propagated.
    /// </summary>
    private static DoubleDouble Overflowed(double h) =>
        new(double.IsFinite(h) ? double.CopySign(double.PositiveInfinity, h) : h, 0.0, normalized: true);

    // ---------------------------------------------------------------
    // Multiplication (DWTimesDW3, 4u^2; DWTimesFP3, 2u^2)
    // ---------------------------------------------------------------

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static DoubleDouble operator *(DoubleDouble a, DoubleDouble b)
    {
        double p = ErrorFree.TwoProd(a._hi, b._hi, out double e);
        double t = Math.FusedMultiplyAdd(a._lo, b._hi, Math.FusedMultiplyAdd(a._hi, b._lo, a._lo * b._lo));
        e += t;
        p = ErrorFree.FastTwoSum(p, e, out e);
        return double.IsFinite(p) ? new DoubleDouble(p, e, normalized: true) : Overflowed(a._hi * b._hi);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static DoubleDouble operator *(DoubleDouble a, double b)
    {
        double p = ErrorFree.TwoProd(a._hi, b, out double e);
        e = Math.FusedMultiplyAdd(a._lo, b, e);
        p = ErrorFree.FastTwoSum(p, e, out e);
        return double.IsFinite(p) ? new DoubleDouble(p, e, normalized: true) : Overflowed(a._hi * b);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static DoubleDouble operator *(double a, DoubleDouble b) => b * a;

    /// <summary>x² (cheaper than x * x, same 4u² bound).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static DoubleDouble Square(DoubleDouble x)
    {
        double p = ErrorFree.TwoSquare(x._hi, out double e);
        e = Math.FusedMultiplyAdd(x._hi + x._hi, x._lo, e);
        p = ErrorFree.FastTwoSum(p, e, out e);
        return double.IsFinite(p) ? new DoubleDouble(p, e, normalized: true) : Overflowed(x._hi * x._hi);
    }

    /// <summary>Exact multiplication by 2^n for a small |n| where no over/underflow can occur.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static DoubleDouble MulPow2(DoubleDouble x, double pow2) => new(x._hi * pow2, x._lo * pow2, normalized: true);

    // ---------------------------------------------------------------
    // Division (DWDivDW2, 15u^2 + 56u^3; DWDivFP3, 3u^2)
    // ---------------------------------------------------------------

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static DoubleDouble operator /(DoubleDouble a, DoubleDouble b)
    {
        double th = a._hi / b._hi;
        if (!double.IsFinite(th) || th == 0.0)
        {
            return DivideSpecial(a, b, th);
        }
        // r = b * th (DWTimesFP3)
        double rh = ErrorFree.TwoProd(b._hi, th, out double rl);
        rl = Math.FusedMultiplyAdd(b._lo, th, rl);
        rh = ErrorFree.FastTwoSum(rh, rl, out rl);
        double ph = a._hi - rh;          // exact (Sterbenz)
        double dl = a._lo - rl;
        double d = ph + dl;
        double tl = d / b._hi;
        double s = ErrorFree.FastTwoSum(th, tl, out double e);
        return double.IsFinite(s) ? new DoubleDouble(s, e, normalized: true) : Overflowed(th);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static DoubleDouble operator /(DoubleDouble a, double b)
    {
        double th = a._hi / b;
        if (!double.IsFinite(th) || th == 0.0)
        {
            return DivideSpecial(a, new DoubleDouble(b, 0.0, normalized: true), th);
        }
        double ph = ErrorFree.TwoProd(th, b, out double pl);
        double dh = a._hi - ph;
        double dt = dh - pl;
        double d = dt + a._lo;
        double tl = d / b;
        double s = ErrorFree.FastTwoSum(th, tl, out double e);
        return double.IsFinite(s) ? new DoubleDouble(s, e, normalized: true) : Overflowed(th);
    }

    public static DoubleDouble operator /(double a, DoubleDouble b) => new DoubleDouble(a, 0.0, normalized: true) / b;

    /// <summary>
    /// Quotients whose leading estimate is 0, ±∞ or NaN: exact zeros/specials follow double
    /// semantics; a finite nonzero quotient that under- or overflowed only in the estimate is
    /// recomputed with both operands scaled by powers of two.
    /// </summary>
    private static DoubleDouble DivideSpecial(DoubleDouble a, DoubleDouble b, double th)
    {
        if (double.IsNaN(th) || a._hi == 0.0 || b._hi == 0.0 || !double.IsFinite(a._hi) || !double.IsFinite(b._hi))
        {
            return new DoubleDouble(th, 0.0, normalized: true);
        }
        // Finite nonzero operands: rescale so the leading quotient is near 1, then scale back.
        int ea = Math.ILogB(a._hi), eb = Math.ILogB(b._hi);
        DoubleDouble q = ScaleB(a, -ea) / ScaleB(b, -eb);
        return ScaleB(q, ea - eb);
    }

    /// <summary>1 / x.</summary>
    public static DoubleDouble Reciprocal(DoubleDouble x) => One / x;

    // ---------------------------------------------------------------
    // Increment / decrement / remainder
    // ---------------------------------------------------------------

    public static DoubleDouble operator ++(DoubleDouble value) => value + 1.0;
    public static DoubleDouble operator --(DoubleDouble value) => value - 1.0;

    /// <summary>
    /// Truncated remainder with the sign of the dividend (like <see cref="double"/> %): computed
    /// exactly through <see cref="BigFloat"/> and then rounded to double-double.
    /// </summary>
    public static DoubleDouble operator %(DoubleDouble a, DoubleDouble b)
    {
        if (!double.IsFinite(a._hi) || double.IsNaN(b._hi) || b._hi == 0.0)
        {
            return NaN;
        }
        if (double.IsInfinity(b._hi) || Abs(a) < Abs(b))
        {
            return a;
        }
        BigFloat r = BigFloat.TruncatedRemainder(a.ToBigFloatExact(), b.ToBigFloatExact(), 2 * BigFloatExactPrecision);
        DoubleDouble result = FromBigFloat(r);
        return result._hi == 0.0 ? CopySign(Zero, a) : result;
    }

    // ---------------------------------------------------------------
    // Square root, cube root, FMA
    // ---------------------------------------------------------------

    /// <summary>√x (Karp/Markstein style correction of the double estimate; about 2^-104).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static DoubleDouble Sqrt(DoubleDouble x)
    {
        double xh = x._hi;
        if (!(xh > 1e-280) || !double.IsFinite(xh))
        {
            return SqrtSpecial(x);
        }
        double sh = Math.Sqrt(xh);
        double r = Math.FusedMultiplyAdd(-sh, sh, xh) + x._lo;
        double sl = r / (sh + sh);
        double s = ErrorFree.FastTwoSum(sh, sl, out double e);
        return new DoubleDouble(s, e, normalized: true);
    }

    private static DoubleDouble SqrtSpecial(DoubleDouble x)
    {
        double xh = x._hi;
        if (xh == 0.0 || double.IsPositiveInfinity(xh))
        {
            return x;
        }
        if (xh > 0.0)
        {
            // Tiny inputs: the residual would underflow, so scale by an even power of two first.
            return ScaleB(Sqrt(ScaleB(x, 600)), -300);
        }
        return NaN;
    }

    /// <summary>
    /// ∛x: two Newton steps from <see cref="Math.Cbrt"/>. The seed's accuracy is platform-dependent
    /// (glibc's cbrt is off by a few ulps), and one quadratic step would only square that error.
    /// </summary>
    public static DoubleDouble Cbrt(DoubleDouble x)
    {
        double xh = x._hi;
        if (xh == 0.0 || !double.IsFinite(xh))
        {
            return x;
        }
        int e = Math.ILogB(xh);
        int k = 0;
        if (e > 600 || e < -600)
        {
            k = e / 3;
            x = ScaleB(x, -3 * k);
        }
        double y = Math.Cbrt(x._hi);
        // y += (x - y^3) / (3 y^2)
        DoubleDouble y2 = new(ErrorFree.TwoSquare(y, out double y2l), y2l, normalized: true);
        DoubleDouble r = x - y2 * y;
        DoubleDouble result = new DoubleDouble(y, 0.0, normalized: true) + r._hi / (3.0 * y2._hi);
        // Second step at double-double precision: the correction only needs ~53 bits.
        DoubleDouble sq = Square(result);
        r = x - sq * result;
        result += r._hi / (3.0 * sq._hi);
        return k == 0 ? result : ScaleB(result, k);
    }

    /// <summary>
    /// (left · right) + addend with the product kept unrounded, accurate to about 2^-106 relative
    /// even under cancellation (e.g. residuals a·b − fl(a·b)). The product is split exactly into
    /// eight doubles by four TwoProd. When the signs or magnitudes rule out cancellation they are
    /// added directly (~13 ns); otherwise the ten terms (with the addend's two words) are distilled
    /// by VecSum passes until a certified bound shows the leading pair is accurate (~40 ns for a
    /// residual). Cancellation beyond ~190 bits falls back to an exact <see cref="BigFloat"/>
    /// evaluation. In generic code that only needs speed, use <see cref="MultiplyAddEstimate"/>
    /// (the unfused form, about twice as fast).
    /// </summary>
    public static DoubleDouble FusedMultiplyAdd(DoubleDouble left, DoubleDouble right, DoubleDouble addend)
    {
        double ah = left._hi, al = left._lo, bh = right._hi, bl = right._lo, ch = addend._hi;
        double p0 = ErrorFree.TwoProd(ah, bh, out double e0);
        if (!double.IsFinite(p0) || !double.IsFinite(ch))
        {
            return Math.FusedMultiplyAdd(ah, bh, ch); // NaN, ±∞ and ∞·0 follow double semantics
        }
        if (p0 == 0.0 || ch == 0.0)
        {
            return left * right + addend; // a zero term: nothing to cancel, keeps signed-zero rules
        }
        double ap = Math.Abs(p0), ac = Math.Abs(ch);
        if (ap < 1e-241 /* ~2^-800 */ || ap > 1e301 /* ~2^1000 */ || ac > 1e301)
        {
            // TwoProd/TwoSum lose exactness near underflow and may overflow near the top of the range.
            return FusedMultiplyAddExact(left, right, addend);
        }

        double cl = addend._lo;
        double p1 = ErrorFree.TwoProd(ah, bl, out double e1);
        double p2 = ErrorFree.TwoProd(al, bh, out double e2);

        // No cancellation possible (same signs, or magnitudes ≥ 4× apart so |result| ≥ max/2):
        // the product to ~2^-159·|ab| as p0 + m + low, then one accurate double-double addition.
        // Absolute errors are ~2^-106·max(|ab|, |c|), so ≤ ~2^-105 relative to the result.
        if ((p0 > 0) == (ch > 0) || ap >= 4 * ac || ac >= 4 * ap)
        {
            double m = ErrorFree.TwoSum(p1, p2, out double r1);
            m = ErrorFree.TwoSum(m, e0, out double r2);
            double low = r1 + r2 + e1 + e2 + al * bl;
            double s = ErrorFree.TwoSum(p0, ch, out double se);
            double v = ErrorFree.TwoSum(m, cl, out double ve);
            v = ErrorFree.TwoSum(v, se, out double ve2);
            double hi0 = ErrorFree.FastTwoSum(s, v, out double lo0);
            return new DoubleDouble(hi0, lo0 + (low + ve + ve2));
        }

        // Possible cancellation: the exact sum of these ten terms is a·b + c. Distil it with VecSum
        // passes (x9 becomes fl(Σ), every other term the exact error of one partial sum) until the
        // terms below (x9, x8) are negligible.
        double x0 = ErrorFree.TwoProd(al, bl, out double e3);
        double x9 = p0, x8 = ch, x7 = e0, x6 = cl, x5 = e2, x4 = p2, x3 = e1, x2 = p1, x1 = x0;
        x0 = e3;
        for (int pass = 0; pass < 5; pass++)
        {
            x1 = ErrorFree.TwoSum(x0, x1, out x0);
            x2 = ErrorFree.TwoSum(x1, x2, out x1);
            x3 = ErrorFree.TwoSum(x2, x3, out x2);
            x4 = ErrorFree.TwoSum(x3, x4, out x3);
            x5 = ErrorFree.TwoSum(x4, x5, out x4);
            x6 = ErrorFree.TwoSum(x5, x6, out x5);
            x7 = ErrorFree.TwoSum(x6, x7, out x6);
            x8 = ErrorFree.TwoSum(x7, x8, out x7);
            x9 = ErrorFree.TwoSum(x8, x9, out x8);
            // After a pass |x8| ≤ ulp(x9)/2, so (x9, x8) is a double-double; it is the answer once
            // the other eight terms are negligible: tail ≤ 2^-108·|x9| keeps their total, plus the
            // final rounding of x8 + Σtail (≤ u·|x8| ≈ 2^-106·|x9|), within ~2^-105.7 relative.
            // The 2^-900 floor keeps the terms clear of subnormal rounding.
            double tail = Math.Abs(x0) + Math.Abs(x1) + Math.Abs(x2) + Math.Abs(x3)
                        + Math.Abs(x4) + Math.Abs(x5) + Math.Abs(x6) + Math.Abs(x7);
            double ah9 = Math.Abs(x9);
            if (tail <= ah9 * (1.0 / (1L << 54) / (1L << 54)) /* 2^-108 */ && (ah9 >= 1e-271 /* ~2^-900 */ || x9 == 0.0))
            {
                return new DoubleDouble(x9, x8 + (((x0 + x1) + (x2 + x3)) + ((x4 + x5) + (x6 + x7))));
            }
        }
        return FusedMultiplyAddExact(left, right, addend);
    }

    private static DoubleDouble FusedMultiplyAddExact(DoubleDouble left, DoubleDouble right, DoubleDouble addend) =>
        FromBigFloat(BigFloat.FusedMultiplyAdd(left.ToBigFloatExact(), right.ToBigFloatExact(), addend.ToBigFloatExact(), 320));

    /// <summary>
    /// left * right + addend, unfused: the product is first rounded to double-double, so the result
    /// carries the multiplication and the addition error bounds (fast; INumberBase).
    /// </summary>
    public static DoubleDouble MultiplyAddEstimate(DoubleDouble left, DoubleDouble right, DoubleDouble addend) => left * right + addend;

    // ---------------------------------------------------------------
    // Scaling
    // ---------------------------------------------------------------

    /// <summary>x · 2^n, exact unless the result over- or underflows.</summary>
    public static DoubleDouble ScaleB(DoubleDouble x, int n)
    {
        double h = Math.ScaleB(x._hi, n);
        if (!double.IsFinite(h) || h == 0.0)
        {
            return new DoubleDouble(h, 0.0, normalized: true);
        }
        double l = Math.ScaleB(x._lo, n);
        h = ErrorFree.FastTwoSum(h, l, out l);
        return new DoubleDouble(h, l, normalized: true);
    }

    /// <summary>floor(log2 |x|); int.MinValue for 0, int.MaxValue for NaN and ±∞ (like <see cref="Math.ILogB(double)"/>).</summary>
    public static int ILogB(DoubleDouble x)
    {
        int e = Math.ILogB(x._hi);
        if (x._lo != 0.0 && double.IsFinite(x._hi) && double.IsPow2(Math.Abs(x._hi)) && (x._lo < 0.0) != (x._hi < 0.0))
        {
            e--; // |hi| is a power of two and lo pulls the magnitude below it
        }
        return e;
    }

    // ---------------------------------------------------------------
    // Rounding to integers
    // ---------------------------------------------------------------

    public static DoubleDouble Floor(DoubleDouble x)
    {
        double hi = Math.Floor(x._hi);
        if (hi != x._hi)
        {
            return new DoubleDouble(hi, 0.0, normalized: true);
        }
        if (!double.IsFinite(hi))
        {
            return x;
        }
        double lo = Math.Floor(x._lo);
        hi = ErrorFree.FastTwoSum(hi, lo, out lo);
        return new DoubleDouble(hi, lo, normalized: true);
    }

    public static DoubleDouble Ceiling(DoubleDouble x)
    {
        double hi = Math.Ceiling(x._hi);
        if (hi != x._hi)
        {
            return new DoubleDouble(hi, 0.0, normalized: true);
        }
        if (!double.IsFinite(hi))
        {
            return x;
        }
        double lo = Math.Ceiling(x._lo);
        hi = ErrorFree.FastTwoSum(hi, lo, out lo);
        return new DoubleDouble(hi, lo, normalized: true);
    }

    public static DoubleDouble Truncate(DoubleDouble x)
    {
        DoubleDouble t = double.IsNegative(x._hi) ? Ceiling(x) : Floor(x);
        return t._hi == 0.0 ? new DoubleDouble(double.CopySign(0.0, x._hi), 0.0, normalized: true) : t;
    }

    /// <summary>Nearest integer, ties to even.</summary>
    public static DoubleDouble Round(DoubleDouble x) => Round(x, MidpointRounding.ToEven);

    public static DoubleDouble Round(DoubleDouble x, MidpointRounding mode)
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
        if (!double.IsFinite(x._hi) || IsInteger(x))
        {
            return x;
        }
        DoubleDouble f = Floor(x);
        DoubleDouble frac = x - f; // exact, in [0, 1)
        DoubleDouble r;
        if (frac._hi > 0.5 || (frac._hi == 0.5 && frac._lo > 0.0))
        {
            r = f + 1.0;
        }
        else if (frac._hi < 0.5 || (frac._hi == 0.5 && frac._lo < 0.0))
        {
            r = f;
        }
        else if (mode == MidpointRounding.AwayFromZero)
        {
            r = double.IsNegative(x._hi) ? f : f + 1.0;
        }
        else
        {
            r = IsEvenInteger(f) ? f : f + 1.0;
        }
        return r._hi == 0.0 ? new DoubleDouble(double.CopySign(0.0, x._hi), 0.0, normalized: true) : r;
    }

    public static DoubleDouble Round(DoubleDouble x, int digits) => Round(x, digits, MidpointRounding.ToEven);

    /// <summary>
    /// Rounds to <paramref name="digits"/> decimal places: the exact decimal rounding (computed
    /// with <see cref="BigFloat"/>), then rounded to the nearest double-double.
    /// </summary>
    public static DoubleDouble Round(DoubleDouble x, int digits, MidpointRounding mode)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(digits);
        if (digits == 0)
        {
            return Round(x, mode);
        }
        if (!double.IsFinite(x._hi) || IsInteger(x))
        {
            return x;
        }
        DoubleDouble r = FromBigFloat(BigFloat.Round(x.ToBigFloatExact(), digits, mode));
        return r._hi == 0.0 ? new DoubleDouble(double.CopySign(0.0, x._hi), 0.0, normalized: true) : r;
    }

    // ---------------------------------------------------------------
    // Min / Max (IEEE 754:2019 semantics, like Math.Min/Math.Max)
    // ---------------------------------------------------------------

    public static DoubleDouble Clamp(DoubleDouble value, DoubleDouble min, DoubleDouble max)
    {
        if (!(min <= max))
        {
            throw new ArgumentException("min must be less than or equal to max.", nameof(min));
        }
        return value < min ? min : value > max ? max : value;
    }

    /// <summary>The larger value; NaN if either is NaN; +0 is larger than -0.</summary>
    public static DoubleDouble Max(DoubleDouble x, DoubleDouble y)
    {
        if (IsNaN(x))
        {
            return x;
        }
        if (IsNaN(y))
        {
            return y;
        }
        if (x == y)
        {
            return double.IsNegative(x._hi) ? y : x;
        }
        return x > y ? x : y;
    }

    /// <summary>The smaller value; NaN if either is NaN; -0 is smaller than +0.</summary>
    public static DoubleDouble Min(DoubleDouble x, DoubleDouble y)
    {
        if (IsNaN(x))
        {
            return x;
        }
        if (IsNaN(y))
        {
            return y;
        }
        if (x == y)
        {
            return double.IsNegative(x._hi) ? x : y;
        }
        return x < y ? x : y;
    }

    /// <summary>Like <see cref="Max(DoubleDouble, DoubleDouble)"/> but a number wins over NaN.</summary>
    public static DoubleDouble MaxNumber(DoubleDouble x, DoubleDouble y) => IsNaN(y) ? x : IsNaN(x) ? y : Max(x, y);

    /// <summary>Like <see cref="Min(DoubleDouble, DoubleDouble)"/> but a number wins over NaN.</summary>
    public static DoubleDouble MinNumber(DoubleDouble x, DoubleDouble y) => IsNaN(y) ? x : IsNaN(x) ? y : Min(x, y);

    public static DoubleDouble MaxMagnitude(DoubleDouble x, DoubleDouble y)
    {
        if (IsNaN(x))
        {
            return x;
        }
        if (IsNaN(y))
        {
            return y;
        }
        DoubleDouble ax = Abs(x), ay = Abs(y);
        return ax > ay ? x : ax < ay ? y : double.IsNegative(x._hi) ? y : x;
    }

    public static DoubleDouble MinMagnitude(DoubleDouble x, DoubleDouble y)
    {
        if (IsNaN(x))
        {
            return x;
        }
        if (IsNaN(y))
        {
            return y;
        }
        DoubleDouble ax = Abs(x), ay = Abs(y);
        return ax < ay ? x : ax > ay ? y : double.IsNegative(x._hi) ? x : y;
    }

    public static DoubleDouble MaxMagnitudeNumber(DoubleDouble x, DoubleDouble y) => IsNaN(y) ? x : IsNaN(x) ? y : MaxMagnitude(x, y);

    public static DoubleDouble MinMagnitudeNumber(DoubleDouble x, DoubleDouble y) => IsNaN(y) ? x : IsNaN(x) ? y : MinMagnitude(x, y);
}
