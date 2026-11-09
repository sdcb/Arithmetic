using System.Runtime.CompilerServices;
using static Sdcb.Arithmetic.Core.ErrorFree;

namespace Sdcb.Arithmetic;

/// <summary>
/// Elementary functions evaluated natively in QD arithmetic (QD-library style): exp by reduction with
/// k·ln 2/64 and a table of 2^(j/64) followed by a Taylor series; log by one Newton correction of a double
/// estimate through expm1 plus a short log1p series; sin/cos by reduction modulo π/256 with tables of
/// sin(kπ/256) and cos(kπ/256) and short Taylor series; atan2 by one rotation correcting a double estimate.
/// Arguments of trigonometric functions beyond about 2^20 are reduced exactly through <see cref="Real"/>.
/// </summary>
public readonly partial struct QuadDouble
{
    /// <summary>1/n! for n = 0..24, computed once in QD arithmetic.</summary>
    private static class Coefficients
    {
        internal static readonly QuadDouble[] InvFact = Build();

        private static QuadDouble[] Build()
        {
            var t = new QuadDouble[25];
            t[0] = One;
            for (int n = 1; n < t.Length; n++)
            {
                t[n] = Divide(t[n - 1], n);
            }
            return t;
        }
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static QuadDouble Entry(ReadOnlySpan<double> table, int k) =>
        new(table[4 * k], table[4 * k + 1], table[4 * k + 2], table[4 * k + 3], Tag.Raw);

    // ---------------------------------------------------------------
    // Exponential
    // ---------------------------------------------------------------

    private const double SixtyFourOverLn2 = 92.33248261689366;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double DdMul(double ah, double al, double bh, double bl, out double lo)
    {
        double p = TwoProd(ah, bh, out double e);
        e += ah * bl + al * bh;
        return FastTwoSum(p, e, out lo);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double DdAdd(double ah, double al, double bh, double bl, out double lo)
    {
        double s = TwoSum(ah, bh, out double e);
        e += al + bl;
        return FastTwoSum(s, e, out lo);
    }

    /// <summary>
    /// expm1(t) for |t| ≲ ln2/128 by its Taylor series to t^21 in Horner form. The inner levels
    /// (t^12 and beyond) only need about 100 bits and run in double-double; levels 2..11 run in QD.
    /// </summary>
    private static QuadDouble ExpM1Taylor(QuadDouble t)
    {
        QuadDouble[] f = Coefficients.InvFact;
        double th = t._x0, tl = t._x1;
        double ph = f[21]._x0, pl = f[21]._x1;
        for (int n = 20; n >= 12; n--)
        {
            double mh = DdMul(th, tl, ph, pl, out double ml);
            ph = DdAdd(f[n]._x0, f[n]._x1, mh, ml, out pl);
        }
        QuadDouble p = new(ph, pl, 0, 0, Tag.Raw);
        for (int n = 11; n >= 2; n--)
        {
            p = AddSloppy(f[n], MulSloppy(t, p));
        }
        return AddSloppy(t, Multiply(Square(t), p));
    }

    /// <summary>x - k·(ln2/64) with k·(ln2/64) formed exactly from five parts (k·part is exact as TwoProd for |k| ≤ 2^17).</summary>
    private static QuadDouble SubtractLn2Over64Multiple(QuadDouble x, double k)
    {
        if (k == 0)
        {
            return x;
        }
        ReadOnlySpan<double> l = Ln2Parts;
        QuadDouble r = x;
        for (int i = 0; i < 5; i++)
        {
            double p = TwoProd(k, l[i] * (1.0 / 64), out double e);
            r = Add(r, -p);
            if (i < 4)
            {
                r = Add(r, -e);
            }
        }
        return r;
    }

    /// <summary>
    /// Given x = k·ln2/64 + t, splits k = 64m + j and returns 2^(j/64)·(1 + expm1(t)) (or that minus
    /// 2^-m when <paramref name="minusOne"/>, before scaling) as w with the result w·2^m.
    /// </summary>
    private static QuadDouble ExpAssemble(QuadDouble t, double k, out int m, bool minusOne)
    {
        QuadDouble e1 = ExpM1Taylor(t);
        double md = Math.Floor(k * (1.0 / 64));
        int j = (int)(k - 64 * md);
        m = (int)md;
        QuadDouble tj = j == 0 ? One : Entry(Exp2Table, j);
        QuadDouble w = j == 0 ? e1 : Multiply(tj, e1);
        if (minusOne)
        {
            // 2^m·tj - 1 is formed first so that the small result keeps its relative accuracy.
            QuadDouble head = Add(ScaleB(tj, m), -1.0);
            m = 0;
            return Add(head, ScaleB(w, (int)md));
        }
        return Add(w, tj);
    }

    /// <summary>exp(x) = 2^m · w with w ∈ [0.99, 2.02]; x must be finite and in range.</summary>
    private static QuadDouble ExpScaled(QuadDouble x, out int m)
    {
        double k = Math.Round(x._x0 * SixtyFourOverLn2);
        return ExpAssemble(SubtractLn2Over64Multiple(x, k), k, out m, minusOne: false);
    }

    /// <summary>e^x (reduction by k·ln2/64 with a table of 2^(j/64), then a Taylor series).</summary>
    public static QuadDouble Exp(QuadDouble x)
    {
        double x0 = x._x0;
        if (double.IsNaN(x0))
        {
            return x;
        }
        if (x0 > 709.8)
        {
            return PositiveInfinity;
        }
        if (x0 < -745.2)
        {
            return Zero;
        }
        if (x0 == 0)
        {
            return One;
        }
        QuadDouble w = ExpScaled(x, out int m);
        return ScaleB(w, m);
    }

    /// <summary>e^x - 1, accurate near zero.</summary>
    public static QuadDouble ExpM1(QuadDouble x)
    {
        double x0 = x._x0;
        if (double.IsNaN(x0) || x0 == 0)
        {
            return x;
        }
        if (x0 > 709.8)
        {
            return PositiveInfinity;
        }
        if (x0 < -160)
        {
            return Add(Exp(x), -1.0);
        }
        if (Math.Abs(x0) < 7.703719777548943E-34 /* 2^-110 */)
        {
            return Add(x, MulPwr2(Square(x), 0.5));
        }
        double k = Math.Round(x0 * SixtyFourOverLn2);
        if (k == 0)
        {
            return ExpM1Taylor(x);
        }
        return ExpAssemble(SubtractLn2Over64Multiple(x, k), k, out _, minusOne: true);
    }

    /// <summary>2^x (the integer part is applied exactly).</summary>
    public static QuadDouble Exp2(QuadDouble x)
    {
        double x0 = x._x0;
        if (double.IsNaN(x0))
        {
            return x;
        }
        if (x0 > 1024.5)
        {
            return PositiveInfinity;
        }
        if (x0 < -1080)
        {
            return Zero;
        }
        double n = Math.Round(x0);
        QuadDouble f = Add(x, -n);
        if (f._x0 == 0)
        {
            return ScaleB(One, (int)n);
        }
        // 2^f = e^(f·ln2) with |f·ln2| ≤ ln2/2; f·64 rounds to the table index directly.
        double k = Math.Round(f._x0 * 64);
        QuadDouble t = Multiply(Add(f, -k * (1.0 / 64)), Ln2);
        QuadDouble w = ExpAssemble(t, k, out int m, minusOne: false);
        return ScaleB(w, m + (int)n);
    }

    /// <summary>2^x - 1, accurate near zero.</summary>
    public static QuadDouble Exp2M1(QuadDouble x)
    {
        if (double.IsNaN(x._x0) || x._x0 == 0)
        {
            return x;
        }
        return Math.Abs(x._x0) < 0.5 ? ExpM1(Multiply(x, Ln2)) : Add(Exp2(x), -1.0);
    }

    /// <summary>10^x, with x·ln10 - k·ln2/64 formed from exact partial products.</summary>
    public static QuadDouble Exp10(QuadDouble x)
    {
        double x0 = x._x0;
        if (double.IsNaN(x0))
        {
            return x;
        }
        if (x0 > 308.3)
        {
            return PositiveInfinity;
        }
        if (x0 < -323.7)
        {
            return Zero;
        }
        if (x0 == 0)
        {
            return One;
        }
        double kd = Math.Round(x0 * 212.6033980727912); // 64·log2(10)
        ReadOnlySpan<double> L = Ln10Parts, l2 = Ln2Parts;
        ReadOnlySpan<double> xs = Components(x);
        const double s64 = 1.0 / 64;

        // Leading terms first: x0·L0 - k·ln2_0/64 cancels to |r| ≲ ln2/128.
        double p = TwoProd(xs[0], L[0], out double e);
        double q = TwoProd(kd, l2[0] * s64, out double eq);
        QuadDouble r = Add(Single(p), -q);
        r = Add(r, e);
        r = Add(r, -eq);
        for (int k = 1; k < 5; k++)
        {
            q = TwoProd(kd, l2[k] * s64, out eq);
            r = Add(r, -q);
            if (k < 4)
            {
                r = Add(r, -eq);
            }
            // Products x_i · L_j with i + j = k.
            for (int i = 0; i <= Math.Min(k, 3); i++)
            {
                int j = k - i;
                if (k < 4)
                {
                    p = TwoProd(xs[i], L[j], out e);
                    r = Add(r, p);
                    r = Add(r, e);
                }
                else
                {
                    r = Add(r, xs[i] * L[j]);
                }
            }
        }
        QuadDouble w = ExpAssemble(r, kd, out int m, minusOne: false);
        return ScaleB(w, m);
    }

    /// <summary>10^x - 1, accurate near zero.</summary>
    public static QuadDouble Exp10M1(QuadDouble x)
    {
        if (double.IsNaN(x._x0) || x._x0 == 0)
        {
            return x;
        }
        return Math.Abs(x._x0) < 0.15 ? ExpM1(Multiply(x, Ln10)) : Add(Exp10(x), -1.0);
    }
    // ---------------------------------------------------------------
    // Logarithms
    // ---------------------------------------------------------------

    /// <summary>
    /// log1p(t) for |t| ≲ 0.42: z0 = log1p(t) in double, y = (1 + t)·e^(-z0) - 1 = t + em + t·em with
    /// em = expm1(-z0) (|y| ≲ 2^-50), then z0 + y - y²/2 + y³/3 - y⁴/4.
    /// </summary>
    private static QuadDouble Log1pCore(QuadDouble t)
    {
        double td = t._x0 + t._x1;
        if (td == 0)
        {
            return t;
        }
        double u = 1.0 + td;
        double z0 = u == 1.0 ? td : Math.Log(u) * td / (u - 1.0);
        QuadDouble em = ExpM1(Single(-z0));
        QuadDouble y = Add(Add(t, em), Multiply(t, em));
        if (y._x0 == 0)
        {
            return Single(z0);
        }
        QuadDouble y2 = Square(y);
        // y³/3 in double-double precision.
        double y3h = TwoProd(y2._x0, y._x0, out double y3l);
        y3l += y2._x1 * y._x0 + y2._x0 * y._x1;
        double c3h = y3h / 3.0;
        double c3l = (Math.FusedMultiplyAdd(-c3h, 3.0, y3h) + y3l) / 3.0;
        double y4 = y2._x0 * y2._x0 * 0.25;
        QuadDouble s = AddSloppy(y, MulPwr2(y2, -0.5));
        s = Add(Add(Add(s, c3h), c3l), -y4);
        return Add(s, z0);
    }

    /// <summary>Splits a finite positive a as 2^e · (1 + t) with 1 + t ∈ [√½, √2).</summary>
    private static QuadDouble LogDecompose(QuadDouble a, out int e)
    {
        e = Math.ILogB(a._x0);
        if (Math.ScaleB(a._x0, -e) > 1.4142135623730951)
        {
            e++;
        }
        return Add(ScaleB(a, -e), -1.0);
    }

    /// <summary>e·ln 2 to full accuracy for an integer |e| ≤ 2^11 (exact partial products).</summary>
    private static QuadDouble Ln2Multiple(double e)
    {
        ReadOnlySpan<double> l = Ln2Parts;
        double p = TwoProd(e, l[0], out double err);
        QuadDouble r = Add(Single(p), err);
        for (int i = 1; i < 5; i++)
        {
            p = TwoProd(e, l[i], out err);
            r = Add(r, p);
            if (i < 4)
            {
                r = Add(r, err);
            }
        }
        return r;
    }

    /// <summary>Handles NaN, negative, zero and +∞ arguments of the logarithms; false for finite positive values.</summary>
    private static bool LogSpecial(QuadDouble x, out QuadDouble result)
    {
        double x0 = x._x0;
        if (double.IsNaN(x0) || x0 < 0)
        {
            result = NaN;
            return true;
        }
        if (x0 == 0)
        {
            result = NegativeInfinity;
            return true;
        }
        if (double.IsPositiveInfinity(x0))
        {
            result = x;
            return true;
        }
        result = default;
        return false;
    }

    /// <summary>Natural logarithm.</summary>
    public static QuadDouble Log(QuadDouble x)
    {
        if (LogSpecial(x, out QuadDouble special))
        {
            return special;
        }
        QuadDouble t = LogDecompose(x, out int e);
        QuadDouble l = Log1pCore(t);
        return e == 0 ? l : Add(Ln2Multiple(e), l);
    }

    /// <summary>log(1 + x), accurate near zero.</summary>
    public static QuadDouble LogP1(QuadDouble x)
    {
        double x0 = x._x0;
        if (double.IsNaN(x0) || x0 < -1)
        {
            return NaN;
        }
        if (x0 == 0)
        {
            return x;
        }
        if (x0 >= -0.29 && x0 <= 0.41)
        {
            return Log1pCore(x);
        }
        return Log(Add(x, 1.0));
    }

    /// <summary>Base-2 logarithm (exact for powers of two).</summary>
    public static QuadDouble Log2(QuadDouble x)
    {
        if (LogSpecial(x, out QuadDouble special))
        {
            return special;
        }
        QuadDouble t = LogDecompose(x, out int e);
        QuadDouble l = t._x0 == 0 ? Zero : Multiply(Log1pCore(t), Log2E);
        return e == 0 ? l : Add(l, e);
    }

    /// <summary>Base-10 logarithm.</summary>
    public static QuadDouble Log10(QuadDouble x)
    {
        if (LogSpecial(x, out QuadDouble special))
        {
            return special;
        }
        return Multiply(Log(x), Log10E);
    }

    /// <summary>log2(1 + x).</summary>
    public static QuadDouble Log2P1(QuadDouble x) => Multiply(LogP1(x), Log2E);

    /// <summary>log10(1 + x).</summary>
    public static QuadDouble Log10P1(QuadDouble x) => Multiply(LogP1(x), Log10E);

    /// <summary>Logarithm of x in base newBase (log x / log newBase).</summary>
    public static QuadDouble Log(QuadDouble x, QuadDouble newBase) => Divide(Log(x), Log(newBase));

    // ---------------------------------------------------------------
    // Powers and roots
    // ---------------------------------------------------------------

    /// <summary>
    /// x^n. Small |n| (≤ 64), or bases with |ln x| ≥ 1, use binary powering (exact for small exact
    /// results, relative error up to about |n|·2^-212); other bases use exp(n·log x).
    /// </summary>
    public static QuadDouble Pow(QuadDouble x, int n)
    {
        if (n == 0)
        {
            return One;
        }
        double x0 = x._x0;
        if (double.IsNaN(x0))
        {
            return x;
        }
        bool binary = n is >= -64 and <= 64 || x0 == 0 || !double.IsFinite(x0) || Math.Abs(Math.Log(Math.Abs(x0))) >= 1;
        if (!binary)
        {
            QuadDouble r = Exp(Multiply(Log(Abs(x)), n));
            return x0 < 0 && (n & 1) != 0 ? -r : r;
        }
        ulong k = n < 0 ? (ulong)(-(long)n) : (ulong)n;
        QuadDouble b = x, acc = One;
        bool first = true;
        while (true)
        {
            if ((k & 1) != 0)
            {
                acc = first ? b : Multiply(acc, b);
                first = false;
            }
            k >>= 1;
            if (k == 0)
            {
                break;
            }
            b = Square(b);
        }
        return n < 0 ? Divide(One, acc) : acc;
    }

    /// <summary>
    /// x^y with the special cases of <see cref="Math.Pow"/>. Integral y uses <see cref="Pow(QuadDouble, int)"/>;
    /// otherwise exp(y·log |x|), whose relative error grows like |y·ln x|·2^-212.
    /// </summary>
    public static QuadDouble Pow(QuadDouble x, QuadDouble y)
    {
        double x0 = x._x0, y0 = y._x0;
        if (y0 == 0 || (x0 == 1 && x._x1 == 0))
        {
            return One;
        }
        if (double.IsNaN(x0) || double.IsNaN(y0))
        {
            return NaN;
        }
        if (double.IsInfinity(y0))
        {
            int c = CompareNonNaN(Abs(x), One);
            return c == 0 ? One : (c > 0) == (y0 > 0) ? PositiveInfinity : Zero;
        }
        bool yInt = IsInteger(y);
        bool yOdd = yInt && OddParity(y);
        if (double.IsInfinity(x0))
        {
            if (x0 > 0)
            {
                return y0 > 0 ? PositiveInfinity : Zero;
            }
            return yOdd ? (y0 > 0 ? NegativeInfinity : NegativeZero) : (y0 > 0 ? PositiveInfinity : Zero);
        }
        if (x0 == 0)
        {
            if (y0 > 0)
            {
                return yOdd ? x : Zero;
            }
            return yOdd && double.IsNegative(x0) ? NegativeInfinity : PositiveInfinity;
        }
        if (x0 < 0 && !yInt)
        {
            return NaN;
        }
        if (yInt && Math.Abs(y0) <= int.MaxValue)
        {
            return Pow(x, (int)(y._x0 + y._x1));
        }
        QuadDouble r = Exp(Multiply(y, Log(Abs(x))));
        return x0 < 0 && yOdd ? -r : r;
    }

    /// <summary>Cube root (three Newton steps from <see cref="Math.Cbrt"/>); odd, so negative inputs give negative roots.</summary>
    public static QuadDouble Cbrt(QuadDouble x)
    {
        double x0 = x._x0;
        if (x0 == 0 || !double.IsFinite(x0))
        {
            return x;
        }
        bool neg = x0 < 0;
        QuadDouble a = neg ? -x : x;
        int scale = 0;
        if (a._x0 > 8.452712498170644E+270 /* 2^900 */ || a._x0 < 1.499696813895631E-241 /* 2^-800 */)
        {
            scale = Math.ILogB(a._x0) / 3;
            a = ScaleB(a, -3 * scale);
        }
        QuadDouble y = Math.Cbrt(a._x0);
        for (int i = 0; i < 3; i++)
        {
            y = Add(y, Divide(Subtract(Divide(a, Square(y)), y), 3.0));
        }
        if (scale != 0)
        {
            y = ScaleB(y, scale);
        }
        return neg ? -y : y;
    }

    /// <summary>The n-th root (IEEE rootn): exp(log|x|/n) refined by one Newton step; odd n accept negative x.</summary>
    public static QuadDouble RootN(QuadDouble x, int n)
    {
        double x0 = x._x0;
        if (n == 0 || double.IsNaN(x0))
        {
            return NaN;
        }
        if (n == 1)
        {
            return x;
        }
        if (n == -1)
        {
            return Divide(One, x);
        }
        bool odd = (n & 1) != 0;
        if (x0 < 0 && !odd)
        {
            return NaN;
        }
        if (x0 == 0)
        {
            QuadDouble z = odd ? x : Zero;
            return n > 0 ? z : Divide(One, z);
        }
        if (double.IsInfinity(x0))
        {
            return n > 0 ? x : odd ? CopySign(Zero, x) : Zero;
        }
        if (n == 2)
        {
            return Sqrt(x);
        }
        if (n == 3)
        {
            return Cbrt(x);
        }
        bool neg = x0 < 0;
        QuadDouble a = neg ? -x : x;
        long m = Math.Abs((long)n);
        QuadDouble y = Exp(Divide(Log(a), m));
        if (double.IsFinite(y._x0) && y._x0 != 0)
        {
            // y ← y + y·(a / y^m - 1) / m
            QuadDouble ym = Pow(y, (int)Math.Min(m, int.MaxValue));
            if (double.IsFinite(ym._x0) && ym._x0 != 0)
            {
                QuadDouble c = Divide(Multiply(y, Add(Divide(a, ym), -1.0)), m);
                y = Add(y, c);
            }
        }
        if (n < 0)
        {
            y = Divide(One, y);
        }
        return neg ? -y : y;
    }

    /// <summary>√(x² + y²) without intermediate overflow or underflow.</summary>
    public static QuadDouble Hypot(QuadDouble x, QuadDouble y)
    {
        if (double.IsInfinity(x._x0) || double.IsInfinity(y._x0))
        {
            return PositiveInfinity;
        }
        if (double.IsNaN(x._x0) || double.IsNaN(y._x0))
        {
            return NaN;
        }
        QuadDouble ax = Abs(x), ay = Abs(y);
        if (ax._x0 == 0)
        {
            return ay;
        }
        if (ay._x0 == 0)
        {
            return ax;
        }
        int e = Math.Max(Math.ILogB(ax._x0), Math.ILogB(ay._x0));
        QuadDouble xs = ScaleB(ax, -e), ys = ScaleB(ay, -e);
        return ScaleB(Sqrt(Add(Square(xs), Square(ys))), e);
    }

    // ---------------------------------------------------------------
    // Trigonometric functions
    // ---------------------------------------------------------------

    /// <summary>
    /// x = j·π/256 + r with |r| ≲ π/512; returns r and j mod 512. Uses six-part π with exact
    /// products for |j| ≤ 2^28 when r is not tiny, otherwise reduces exactly through <see cref="Real"/>.
    /// </summary>
    private static QuadDouble ReduceTrig(QuadDouble x, out int octant)
    {
        double jd = Math.Round(x._x0 * 81.48733086305042); // 256/π
        if (jd == 0)
        {
            octant = 0;
            return x;
        }
        if (Math.Abs(jd) <= 268435456.0 /* 2^28 */)
        {
            ReadOnlySpan<double> pi = PiParts;
            QuadDouble r = x;
            for (int i = 0; i < 6; i++)
            {
                double p = TwoProd(jd, pi[i] * (1.0 / 256), out double e);
                r = Add(r, -p);
                if (i < 5)
                {
                    r = Add(r, -e);
                }
            }
            if (Math.Abs(r._x0) > 8.470329472543003E-22 /* 2^-70 */)
            {
                octant = (int)((long)jd & 511);
                return r;
            }
        }
        // Exact reduction through Real for huge arguments or near-multiples of π/256.
        Real rx = ToRealExact(x);
        int prec = Math.Max(Math.ILogB(x._x0), 0) + rx.Precision + 600;
        Real step = Real.Pi(prec).ShiftRight(8);
        Real q = Real.Round(Real.Divide(rx, step, prec));
        Real rr = Real.Subtract(rx, Real.Multiply(q, step, prec), prec);
        BigInt qi = q.ToBigInt();
        octant = (int)(long)(qi & new BigInt(511));
        return FromReal(rr);
    }

    /// <summary>
    /// sin(r) and cos(r) for |r| ≲ π/512 by Taylor series in r² (eleven terms each). Levels 6..11 need
    /// under 100 bits and run in double-double; levels 1..5 run in QD.
    /// </summary>
    private static void SinCosSeries(QuadDouble r, out QuadDouble sin, out QuadDouble cos)
    {
        QuadDouble[] f = Coefficients.InvFact;
        QuadDouble r2 = Square(r);
        double rh = r2._x0, rl = r2._x1;
        // Coefficients (-1)^k/(2k+1)! and (-1)^k/(2k)!, starting at k = 11 (negative).
        double sh = -f[23]._x0, sl = -f[23]._x1, ch = -f[22]._x0, cl = -f[22]._x1;
        for (int k = 10; k >= 6; k--)
        {
            double sign = (k & 1) != 0 ? -1.0 : 1.0;
            double mh = DdMul(rh, rl, sh, sl, out double ml);
            sh = DdAdd(sign * f[2 * k + 1]._x0, sign * f[2 * k + 1]._x1, mh, ml, out sl);
            mh = DdMul(rh, rl, ch, cl, out ml);
            ch = DdAdd(sign * f[2 * k]._x0, sign * f[2 * k]._x1, mh, ml, out cl);
        }
        QuadDouble ps = new(sh, sl, 0, 0, Tag.Raw), pc = new(ch, cl, 0, 0, Tag.Raw);
        for (int k = 5; k >= 1; k--)
        {
            QuadDouble cs = f[2 * k + 1], cc = f[2 * k];
            if ((k & 1) != 0)
            {
                cs = -cs;
                cc = -cc;
            }
            ps = AddSloppy(cs, MulSloppy(r2, ps));
            pc = AddSloppy(cc, MulSloppy(r2, pc));
        }
        sin = AddSloppy(r, Multiply(r, Multiply(r2, ps)));
        cos = AddSloppy(Multiply(r2, pc), One);
    }
    /// <summary>(sin x, cos x) for finite nonzero x.</summary>
    private static void SinCosCore(QuadDouble x, out QuadDouble sin, out QuadDouble cos)
    {
        QuadDouble r = ReduceTrig(x, out int j);
        SinCosSeries(r, out QuadDouble sr, out QuadDouble cr);
        int k = j & 127;
        QuadDouble s, c;
        if (k == 0)
        {
            s = sr;
            c = cr;
        }
        else
        {
            QuadDouble sk, ck;
            if (k <= 64)
            {
                sk = Entry(SinTable, k);
                ck = Entry(CosTable, k);
            }
            else
            {
                sk = Entry(CosTable, 128 - k);
                ck = Entry(SinTable, 128 - k);
            }
            s = Add(Multiply(sk, cr), Multiply(ck, sr));
            c = Subtract(Multiply(ck, cr), Multiply(sk, sr));
        }
        switch (j >> 7)
        {
            case 0:
                sin = s;
                cos = c;
                break;
            case 1:
                sin = c;
                cos = -s;
                break;
            case 2:
                sin = -s;
                cos = -c;
                break;
            default:
                sin = -c;
                cos = s;
                break;
        }
    }

    /// <summary>Simultaneous sine and cosine.</summary>
    public static (QuadDouble Sin, QuadDouble Cos) SinCos(QuadDouble x)
    {
        if (!double.IsFinite(x._x0))
        {
            return (NaN, NaN);
        }
        if (x._x0 == 0)
        {
            return (x, One);
        }
        SinCosCore(x, out QuadDouble s, out QuadDouble c);
        return (s, c);
    }

    /// <summary>Sine.</summary>
    public static QuadDouble Sin(QuadDouble x) => SinCos(x).Sin;

    /// <summary>Cosine.</summary>
    public static QuadDouble Cos(QuadDouble x) => SinCos(x).Cos;

    /// <summary>Tangent (sin/cos).</summary>
    public static QuadDouble Tan(QuadDouble x)
    {
        (QuadDouble s, QuadDouble c) = SinCos(x);
        return x._x0 == 0 ? x : Divide(s, c);
    }

    /// <summary>(sin πx, cos πx) with exact reduction modulo 2 (exact zeros and ±1 at multiples of 1/2).</summary>
    public static (QuadDouble SinPi, QuadDouble CosPi) SinCosPi(QuadDouble x)
    {
        double x0 = x._x0;
        if (!double.IsFinite(x0))
        {
            return (NaN, NaN);
        }
        if (x0 == 0)
        {
            return (x, One);
        }
        // n = round(2x), f = x - n/2 ∈ [-1/4, 1/4] (exact), q = n mod 4.
        QuadDouble n = Round(MulPwr2(x, 2.0));
        QuadDouble f = Subtract(x, MulPwr2(n, 0.5));
        double qd = 0;
        foreach (double c in Components(n))
        {
            qd += c % 4.0; // each component is an integer; the remainders are exact
        }
        int q = (int)(((long)qd % 4 + 4) % 4);
        QuadDouble s, co;
        if (f._x0 == 0)
        {
            QuadDouble zero = CopySign(Zero, x);
            return q switch
            {
                0 => (zero, One),
                1 => (One, Zero),
                2 => (zero, NegativeOne),
                _ => (NegativeOne, Zero),
            };
        }
        SinCosCore(Multiply(f, Pi), out s, out co);
        return q switch
        {
            0 => (s, co),
            1 => (co, -s),
            2 => (-s, -co),
            _ => (-co, s),
        };
    }

    /// <summary>sin(πx).</summary>
    public static QuadDouble SinPi(QuadDouble x) => SinCosPi(x).SinPi;

    /// <summary>cos(πx).</summary>
    public static QuadDouble CosPi(QuadDouble x) => SinCosPi(x).CosPi;

    /// <summary>tan(πx); ±∞ at odd multiples of 1/2.</summary>
    public static QuadDouble TanPi(QuadDouble x)
    {
        (QuadDouble s, QuadDouble c) = SinCosPi(x);
        return x._x0 == 0 ? x : Divide(s, c);
    }

    /// <summary>
    /// atan2(y, x): z0 = Math.Atan2 on the leading components, then the residual angle of (x, y)
    /// rotated by -z0 (tiny, ≲ 2^-50) is added via atan(t) ≈ t - t³/3.
    /// </summary>
    public static QuadDouble Atan2(QuadDouble y, QuadDouble x)
    {
        double y0 = y._x0, x0 = x._x0;
        if (double.IsNaN(y0) || double.IsNaN(x0))
        {
            return NaN;
        }
        if (y0 == 0)
        {
            return double.IsNegative(x0) ? CopySign(Pi, y) : y;
        }
        if (double.IsInfinity(y0))
        {
            QuadDouble a = double.IsInfinity(x0) ? (x0 > 0 ? PiOver4 : Multiply(PiOver4, 3.0)) : PiOver2;
            return CopySign(a, y);
        }
        if (x0 == 0)
        {
            return CopySign(PiOver2, y);
        }
        if (double.IsInfinity(x0))
        {
            return x0 > 0 ? CopySign(Zero, y) : CopySign(Pi, y);
        }
        int ey = Math.ILogB(y0), ex = Math.ILogB(x0);
        if (ey - ex > 900)
        {
            return CopySign(PiOver2, y); // |x/y| < 2^-899
        }
        if (ex - ey > 900)
        {
            return x0 > 0 ? Divide(y, x) : CopySign(Pi, y);
        }
        int e = Math.Max(ex, ey);
        QuadDouble xs = ScaleB(x, -e), ys = ScaleB(y, -e);
        double z0 = Math.Atan2(ys._x0, xs._x0);
        SinCosCore(Single(z0), out QuadDouble s, out QuadDouble c);
        QuadDouble num = Subtract(Multiply(ys, c), Multiply(xs, s));
        QuadDouble den = Add(Multiply(xs, c), Multiply(ys, s));
        QuadDouble t = Divide(num, den);
        // atan(t) = t - t³/3 + O(t⁵), t³/3 in double-double.
        double t2h = TwoSquare(t._x0, out double t2l);
        double t3h = TwoProd(t2h, t._x0, out double t3l);
        t3l += t2l * t._x0 + 3.0 * t2h * t._x1;
        double c3h = t3h / 3.0;
        double c3l = (Math.FusedMultiplyAdd(-c3h, 3.0, t3h) + t3l) / 3.0;
        QuadDouble d = Add(Add(t, -c3h), -c3l);
        return Add(d, z0);
    }

    /// <summary>Arctangent.</summary>
    public static QuadDouble Atan(QuadDouble x)
    {
        if (x._x0 == 0 || double.IsNaN(x._x0))
        {
            return x;
        }
        return Atan2(x, One);
    }

    /// <summary>Arcsine (atan2(x, √((1-x)(1+x)))); NaN outside [-1, 1].</summary>
    public static QuadDouble Asin(QuadDouble x)
    {
        if (x._x0 == 0 || double.IsNaN(x._x0))
        {
            return x;
        }
        int c = CompareNonNaN(Abs(x), One);
        if (c > 0)
        {
            return NaN;
        }
        if (c == 0)
        {
            return CopySign(PiOver2, x);
        }
        return Atan2(x, Sqrt(Multiply(Add(-x, 1.0), Add(x, 1.0))));
    }

    /// <summary>Arccosine (atan2(√((1-x)(1+x)), x)); NaN outside [-1, 1].</summary>
    public static QuadDouble Acos(QuadDouble x)
    {
        if (double.IsNaN(x._x0))
        {
            return x;
        }
        int c = CompareNonNaN(Abs(x), One);
        if (c > 0)
        {
            return NaN;
        }
        if (c == 0)
        {
            return x._x0 > 0 ? Zero : Pi;
        }
        return Atan2(Sqrt(Multiply(Add(-x, 1.0), Add(x, 1.0))), x);
    }

    /// <summary>asin(x)/π (exact ±1/2 at ±1).</summary>
    public static QuadDouble AsinPi(QuadDouble x)
    {
        if (x._x0 == 0 || double.IsNaN(x._x0))
        {
            return x;
        }
        int c = CompareNonNaN(Abs(x), One);
        return c > 0 ? NaN : c == 0 ? CopySign(Single(0.5), x) : Divide(Asin(x), Pi);
    }

    /// <summary>acos(x)/π (exact 0 and 1 at ±1, 1/2 at 0).</summary>
    public static QuadDouble AcosPi(QuadDouble x)
    {
        if (double.IsNaN(x._x0))
        {
            return x;
        }
        if (x._x0 == 0)
        {
            return Single(0.5);
        }
        int c = CompareNonNaN(Abs(x), One);
        return c > 0 ? NaN : c == 0 ? (x._x0 > 0 ? Zero : One) : Divide(Acos(x), Pi);
    }

    /// <summary>atan(x)/π (exact ±1/4 at ±1, ±1/2 at ±∞).</summary>
    public static QuadDouble AtanPi(QuadDouble x)
    {
        if (x._x0 == 0 || double.IsNaN(x._x0))
        {
            return x;
        }
        if (double.IsInfinity(x._x0))
        {
            return CopySign(Single(0.5), x);
        }
        return CompareNonNaN(Abs(x), One) == 0 ? CopySign(Single(0.25), x) : Divide(Atan(x), Pi);
    }

    /// <summary>atan2(y, x)/π (exact at the axes, the diagonals and infinities).</summary>
    public static QuadDouble Atan2Pi(QuadDouble y, QuadDouble x)
    {
        double y0 = y._x0, x0 = x._x0;
        if (double.IsNaN(y0) || double.IsNaN(x0))
        {
            return NaN;
        }
        if (y0 == 0)
        {
            return double.IsNegative(x0) ? CopySign(One, y) : y;
        }
        if (double.IsInfinity(y0))
        {
            return CopySign(Single(double.IsInfinity(x0) ? (x0 > 0 ? 0.25 : 0.75) : 0.5), y);
        }
        if (x0 == 0)
        {
            return CopySign(Single(0.5), y);
        }
        if (double.IsInfinity(x0))
        {
            return x0 > 0 ? CopySign(Zero, y) : CopySign(One, y);
        }
        if (CompareMagnitude(y, x) == 0)
        {
            return CopySign(Single(x0 > 0 ? 0.25 : 0.75), y);
        }
        return Divide(Atan2(y, x), Pi);
    }

    /// <summary>degrees · π/180.</summary>
    public static QuadDouble DegreesToRadians(QuadDouble degrees) => Multiply(degrees, PiOver180);

    /// <summary>radians · 180/π.</summary>
    public static QuadDouble RadiansToDegrees(QuadDouble radians) => Multiply(radians, D180OverPi);

    // ---------------------------------------------------------------
    // Hyperbolic functions
    // ---------------------------------------------------------------

    /// <summary>Hyperbolic sine (expm1-based near zero, overflow-safe scaling for large |x|).</summary>
    public static QuadDouble Sinh(QuadDouble x)
    {
        double x0 = x._x0;
        if (!double.IsFinite(x0) || x0 == 0)
        {
            return x;
        }
        QuadDouble a = Abs(x);
        QuadDouble r;
        if (a._x0 < 0.5)
        {
            // sinh = em (em + 2) / (2 (em + 1)), em = expm1(|x|)
            QuadDouble em = ExpM1(a);
            r = Divide(Multiply(em, Add(em, 2.0)), MulPwr2(Add(em, 1.0), 2.0));
        }
        else if (a._x0 > 75)
        {
            r = a._x0 > 711 ? PositiveInfinity : ScaleB(ExpScaled(a, out int m), m - 1);
        }
        else
        {
            QuadDouble e = Exp(a);
            r = MulPwr2(Subtract(e, Divide(One, e)), 0.5);
        }
        return x0 < 0 ? -r : r;
    }

    /// <summary>Hyperbolic cosine.</summary>
    public static QuadDouble Cosh(QuadDouble x)
    {
        double x0 = x._x0;
        if (double.IsNaN(x0))
        {
            return x;
        }
        if (double.IsInfinity(x0))
        {
            return PositiveInfinity;
        }
        QuadDouble a = Abs(x);
        if (a._x0 > 75)
        {
            return a._x0 > 711 ? PositiveInfinity : ScaleB(ExpScaled(a, out int m), m - 1);
        }
        QuadDouble e = Exp(a);
        return MulPwr2(Add(e, Divide(One, e)), 0.5);
    }

    /// <summary>Hyperbolic tangent (expm1(2x) / (expm1(2x) + 2)).</summary>
    public static QuadDouble Tanh(QuadDouble x)
    {
        double x0 = x._x0;
        if (double.IsNaN(x0) || x0 == 0)
        {
            return x;
        }
        if (Math.Abs(x0) > 75)
        {
            return x0 > 0 ? One : NegativeOne;
        }
        QuadDouble em = ExpM1(MulPwr2(x, 2.0));
        return Divide(em, Add(em, 2.0));
    }

    /// <summary>Inverse hyperbolic sine (log1p form near zero).</summary>
    public static QuadDouble Asinh(QuadDouble x)
    {
        double x0 = x._x0;
        if (!double.IsFinite(x0) || x0 == 0)
        {
            return x;
        }
        QuadDouble a = Abs(x);
        QuadDouble r;
        if (a._x0 > 8.112963841460668E+31 /* 2^106 */)
        {
            r = Add(Log(a), Ln2);
        }
        else
        {
            // asinh(a) = log1p(a + a² / (1 + √(1 + a²)))
            QuadDouble a2 = Square(a);
            r = LogP1(Add(a, Divide(a2, Add(Sqrt(Add(a2, 1.0)), 1.0))));
        }
        return x0 < 0 ? -r : r;
    }

    /// <summary>Inverse hyperbolic cosine; NaN below 1.</summary>
    public static QuadDouble Acosh(QuadDouble x)
    {
        double x0 = x._x0;
        if (double.IsNaN(x0) || x0 < 1)
        {
            return NaN;
        }
        if (double.IsPositiveInfinity(x0))
        {
            return x;
        }
        if (x0 > 8.112963841460668E+31 /* 2^106 */)
        {
            return Add(Log(x), Ln2);
        }
        QuadDouble t = Add(x, -1.0);
        if (t._x0 < 0)
        {
            return NaN;
        }
        // acosh(1 + t) = log1p(t + √(t (t + 2)))
        return LogP1(Add(t, Sqrt(Multiply(t, Add(t, 2.0)))));
    }

    /// <summary>Inverse hyperbolic tangent; NaN outside [-1, 1], ±∞ at ±1.</summary>
    public static QuadDouble Atanh(QuadDouble x)
    {
        double x0 = x._x0;
        if (double.IsNaN(x0) || x0 == 0)
        {
            return x;
        }
        int c = CompareNonNaN(Abs(x), One);
        if (c > 0)
        {
            return NaN;
        }
        if (c == 0)
        {
            return x0 > 0 ? PositiveInfinity : NegativeInfinity;
        }
        // atanh(x) = log1p(2x / (1 - x)) / 2
        return MulPwr2(LogP1(Divide(MulPwr2(x, 2.0), Add(-x, 1.0))), 0.5);
    }
}
