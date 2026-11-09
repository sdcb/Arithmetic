using System.Runtime.CompilerServices;
using Sdcb.Arithmetic.Core;

namespace Sdcb.Arithmetic;

public readonly partial struct DoubleDouble
{
    // ---------------------------------------------------------------
    // Exponential kernel
    // ---------------------------------------------------------------

    /// <summary>
    /// expm1(r) for |r| &lt;= ~0.36: r is scaled by 2^-k (|r'| &lt; 2^-6), a degree-14 Taylor
    /// polynomial is evaluated (DD Horner for the 1/2!..1/7! terms, double Horner for the tail,
    /// which contributes below 2^-57 relative), then k doublings s = 2s + s² undo the scaling while
    /// keeping the expm1 form (so the result stays accurate relative to expm1(r) itself).
    /// </summary>
    private static DoubleDouble ExpM1Kernel(DoubleDouble r)
    {
        int k = Math.Max(0, Math.ILogB(r._hi) + 7);
        if (k > 0)
        {
            r = MulPow2(r, Math.ScaleB(1.0, -k));
        }
        ReadOnlySpan<double> fh = InvFactHi, fl = InvFactLo;
        double rh = r._hi;
        double q = fh[14];
        q = Math.FusedMultiplyAdd(q, rh, fh[13]);
        q = Math.FusedMultiplyAdd(q, rh, fh[12]);
        q = Math.FusedMultiplyAdd(q, rh, fh[11]);
        q = Math.FusedMultiplyAdd(q, rh, fh[10]);
        q = Math.FusedMultiplyAdd(q, rh, fh[9]);
        q = Math.FusedMultiplyAdd(q, rh, fh[8]);
        DoubleDouble p = new DoubleDouble(fh[7], fl[7], normalized: true) + rh * q;
        p = new DoubleDouble(fh[6], fl[6], normalized: true) + r * p;
        p = new DoubleDouble(fh[5], fl[5], normalized: true) + r * p;
        p = new DoubleDouble(fh[4], fl[4], normalized: true) + r * p;
        p = new DoubleDouble(fh[3], fl[3], normalized: true) + r * p;
        p = new DoubleDouble(0.5, 0.0, normalized: true) + r * p;
        DoubleDouble s = r + Square(r) * p;
        for (int i = 0; i < k; i++)
        {
            s = MulPow2(s, 2.0) + Square(s);
        }
        return s;
    }

    /// <summary>x - n·ln2 with a three-word ln 2 (each n·Ln2_i product is exact for |n| &lt; 2^11 or made exact by TwoProd).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static DoubleDouble ReduceLn2(DoubleDouble x, double n)
    {
        double p0 = ErrorFree.TwoProd(n, Ln2_0, out double e0);
        DoubleDouble r = x + new DoubleDouble(-p0, -e0, normalized: true);
        double p1 = ErrorFree.TwoProd(n, Ln2_1, out double e1);
        r += new DoubleDouble(-p1, -e1, normalized: true);
        return r - n * Ln2_2;
    }

    /// <summary>(1 + s)·2^n.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static DoubleDouble OnePlusScaled(DoubleDouble s, double n)
    {
        DoubleDouble v = s + 1.0;
        return n == 0.0 ? v : ScaleB(v, (int)n);
    }

    /// <summary>Result of exp-type functions for NaN, ±∞ and arguments beyond the overflow/underflow thresholds.</summary>
    private static bool ExpSpecial(DoubleDouble x, double overflow, double underflow, out DoubleDouble result)
    {
        double h = x._hi;
        if (double.IsNaN(h))
        {
            result = x;
            return true;
        }
        if (h > overflow)
        {
            result = PositiveInfinity;
            return true;
        }
        if (h < underflow)
        {
            result = Zero;
            return true;
        }
        result = default;
        return false;
    }

    // ---------------------------------------------------------------
    // Exponentials
    // ---------------------------------------------------------------

    /// <summary>e^x.</summary>
    public static DoubleDouble Exp(DoubleDouble x)
    {
        if (ExpSpecial(x, 709.79, -746.0, out DoubleDouble special))
        {
            return special;
        }
        return ExpCore(x, 0);
    }

    /// <summary>e^x · 2^shift for finite |x| &lt;= ~750 (no range checks; the scaling is exact).</summary>
    private static DoubleDouble ExpCore(DoubleDouble x, int shift)
    {
        double n = Math.Round(x._hi * Log2E0);
        DoubleDouble r = n == 0.0 ? x : ReduceLn2(x, n);
        return OnePlusScaled(ExpM1Kernel(r), n + shift);
    }

    /// <summary>e^x - 1, accurate relative to the result also for tiny x.</summary>
    public static DoubleDouble ExpM1(DoubleDouble x)
    {
        if (Math.Abs(x._hi) < 0.34)
        {
            return ExpM1Kernel(x);
        }
        return Exp(x) - 1.0;
    }

    /// <summary>2^x.</summary>
    public static DoubleDouble Exp2(DoubleDouble x)
    {
        if (ExpSpecial(x, 1024.5, -1076.0, out DoubleDouble special))
        {
            return special;
        }
        double n = Math.Round(x._hi);
        DoubleDouble f = x - n; // exact, |f| <= 1/2
        return OnePlusScaled(ExpM1Kernel(f * Ln2), n);
    }

    /// <summary>2^x - 1.</summary>
    public static DoubleDouble Exp2M1(DoubleDouble x) =>
        Math.Abs(x._hi) < 0.5 ? ExpM1Kernel(x * Ln2) : Exp2(x) - 1.0;

    /// <summary>10^x (the argument reduction x·ln10 - n·ln2 is carried out with three-word constants).</summary>
    public static DoubleDouble Exp10(DoubleDouble x)
    {
        if (ExpSpecial(x, 308.26, -324.0, out DoubleDouble special))
        {
            return special;
        }
        return OnePlusScaled(ExpM1Kernel(ReduceLn10(x, out double n)), n);
    }

    /// <summary>10^x - 1.</summary>
    public static DoubleDouble Exp10M1(DoubleDouble x) =>
        Math.Abs(x._hi) < 0.15 ? ExpM1Kernel(x * Ln10) : Exp10(x) - 1.0;

    /// <summary>r = x·ln10 - n·ln2 with |r| &lt;= ~ln2/2, to about 2^-106 absolute.</summary>
    private static DoubleDouble ReduceLn10(DoubleDouble x, out double n)
    {
        n = Math.Round(x._hi * Log2_10_0);
        double ph = ErrorFree.TwoProd(x._hi, Ln10_0, out double pl);
        double nh = ErrorFree.TwoProd(n, Ln2_0, out double nl);
        DoubleDouble r = new DoubleDouble(ph, pl, normalized: true) + new DoubleDouble(-nh, -nl, normalized: true);
        double ah = ErrorFree.TwoProd(x._hi, Ln10_1, out double al);
        r += new DoubleDouble(ah, al, normalized: true);
        double bh = ErrorFree.TwoProd(x._lo, Ln10_0, out double bl);
        r += new DoubleDouble(bh, bl, normalized: true);
        double ch = ErrorFree.TwoProd(n, Ln2_1, out double cl);
        r += new DoubleDouble(-ch, -cl, normalized: true);
        return r + (Math.FusedMultiplyAdd(x._lo, Ln10_1, x._hi * Ln10_2) - n * Ln2_2);
    }

    // ---------------------------------------------------------------
    // Logarithms
    // ---------------------------------------------------------------

    /// <summary>log1p(u) of a double, accurate to a few ulps (Goldberg's trick).</summary>
    private static double Log1pApprox(double u)
    {
        double w = 1.0 + u;
        return w == 1.0 ? u : Math.Log(w) * (u / (w - 1.0));
    }

    /// <summary>
    /// log(1 + u) for |u| &lt;= ~0.42: one Newton step from the double estimate x0, written as
    /// x0 + (u·e^-x0 + expm1(-x0)) so that no cancellation against 1 occurs.
    /// </summary>
    private static DoubleDouble LogP1Kernel(DoubleDouble u)
    {
        double x0 = Log1pApprox(u._hi);
        DoubleDouble m = ExpM1Kernel(new DoubleDouble(-x0, 0.0, normalized: true));
        DoubleDouble res = u * (m + 1.0) + m;
        return res + x0;
    }

    /// <summary>Splits positive finite x into 2^e · m with m in [√½, √2) and returns log(m).</summary>
    private static DoubleDouble LogReduce(DoubleDouble x, out int e)
    {
        e = Math.ILogB(x._hi);
        DoubleDouble m = e == 0 ? x : ScaleB(x, -e);
        if (m._hi > 1.4142135623730951)
        {
            m = MulPow2(m, 0.5);
            e++;
        }
        return LogP1Kernel(m - 1.0);
    }

    /// <summary>Results of log-type functions for x &lt;= 0, NaN and +∞; false for positive finite x.</summary>
    private static bool LogSpecial(DoubleDouble x, out DoubleDouble result)
    {
        double h = x._hi;
        if (h > 0.0 && double.IsFinite(h))
        {
            result = default;
            return false;
        }
        result = h == 0.0 ? NegativeInfinity : double.IsPositiveInfinity(h) ? x : NaN;
        return true;
    }

    /// <summary>Natural logarithm.</summary>
    public static DoubleDouble Log(DoubleDouble x)
    {
        if (LogSpecial(x, out DoubleDouble special))
        {
            return special;
        }
        DoubleDouble lm = LogReduce(x, out int e);
        return e == 0 ? lm : Ln2 * (double)e + lm;
    }

    /// <summary>Base-2 logarithm (exact for powers of two).</summary>
    public static DoubleDouble Log2(DoubleDouble x)
    {
        if (LogSpecial(x, out DoubleDouble special))
        {
            return special;
        }
        DoubleDouble lm = LogReduce(x, out int e);
        return lm * new DoubleDouble(Log2E0, Log2E1, normalized: true) + (double)e;
    }

    /// <summary>Base-10 logarithm.</summary>
    public static DoubleDouble Log10(DoubleDouble x)
    {
        if (LogSpecial(x, out DoubleDouble special))
        {
            return special;
        }
        DoubleDouble lm = LogReduce(x, out int e);
        DoubleDouble r = lm * new DoubleDouble(Log10E0, Log10E1, normalized: true);
        return e == 0 ? r : new DoubleDouble(Log10_2_0, Log10_2_1, normalized: true) * (double)e + r;
    }

    /// <summary>Logarithm of x in base newBase: log(x) / log(newBase).</summary>
    public static DoubleDouble Log(DoubleDouble x, DoubleDouble newBase) => Log(x) / Log(newBase);

    /// <summary>log(1 + x), accurate relative to the result also for tiny x.</summary>
    public static DoubleDouble LogP1(DoubleDouble x)
    {
        double h = x._hi;
        if (h >= -0.29 && h <= 0.41)
        {
            return LogP1Kernel(x);
        }
        return Log(x + 1.0);
    }

    /// <summary>log2(1 + x).</summary>
    public static DoubleDouble Log2P1(DoubleDouble x) => LogP1(x) * new DoubleDouble(Log2E0, Log2E1, normalized: true);

    /// <summary>log10(1 + x).</summary>
    public static DoubleDouble Log10P1(DoubleDouble x) => LogP1(x) * new DoubleDouble(Log10E0, Log10E1, normalized: true);

    // ---------------------------------------------------------------
    // Powers and roots
    // ---------------------------------------------------------------

    /// <summary>
    /// x^n by binary powering (relative error grows roughly like |n|·2^-104, the conditioning of
    /// x^n itself); Pow(x, 0) is 1 for every x, including NaN.
    /// </summary>
    public static DoubleDouble Pow(DoubleDouble x, int n)
    {
        if (n == 0)
        {
            return One;
        }
        long m = Math.Abs((long)n);
        DoubleDouble result = One, b = x;
        bool first = true;
        while (true)
        {
            if ((m & 1) != 0)
            {
                result = first ? b : result * b;
                first = false;
            }
            m >>= 1;
            if (m == 0)
            {
                break;
            }
            b = Square(b);
        }
        return n < 0 ? One / result : result;
    }

    /// <summary>
    /// x^y with IEEE pow special cases. Small integer exponents use binary powering; otherwise
    /// exp(y·log|x|), whose relative error is about |y·log x|·2^-104 (the error of the logarithm
    /// is amplified by the exponential).
    /// </summary>
    public static DoubleDouble Pow(DoubleDouble x, DoubleDouble y)
    {
        double xh = x._hi, yh = y._hi;
        if (yh == 0.0 || x == One)
        {
            return One;
        }
        if (double.IsNaN(xh) || double.IsNaN(yh))
        {
            return NaN;
        }
        bool yInt = IsInteger(y);
        if (double.IsInfinity(yh))
        {
            DoubleDouble ax = Abs(x);
            if (ax == One)
            {
                return One;
            }
            return (ax > One) == (yh > 0) ? PositiveInfinity : Zero;
        }
        bool negate = double.IsNegative(xh) && yInt && IsOddInteger(y);
        if (xh == 0.0 || double.IsInfinity(xh))
        {
            double r = Math.Pow(Math.Abs(xh), yh); // 0 or +∞ by the sign of y
            return negate ? -r : r;
        }
        if (xh < 0.0 && !yInt)
        {
            return NaN;
        }
        DoubleDouble result;
        if (yInt && Math.Abs(yh) <= 64)
        {
            result = Pow(Abs(x), (int)yh);
        }
        else
        {
            result = Exp(y * Log(Abs(x)));
        }
        return negate ? -result : result;
    }

    /// <summary>
    /// The n-th root: Newton iteration on y^n = x after scaling x by 2^(n·k) towards 1 (falls back
    /// to exp(log(x)/n) when the scaled value would be subnormal, i.e. for very large n).
    /// </summary>
    public static DoubleDouble RootN(DoubleDouble x, int n)
    {
        if (n == 0 || IsNaN(x))
        {
            return NaN;
        }
        if (n < 0)
        {
            DoubleDouble pos = n == int.MinValue ? RootN(Sqrt(x), -(n / 2)) : RootN(x, -n);
            return One / pos;
        }
        if (n == 1)
        {
            return x;
        }
        bool odd = (n & 1) != 0;
        double xh = x._hi;
        if (xh < 0.0)
        {
            return odd ? -RootN(-x, n) : NaN;
        }
        if (n == 2)
        {
            return Sqrt(x);
        }
        if (xh == 0.0 || double.IsInfinity(xh))
        {
            return odd ? x : Abs(x);
        }
        int e = Math.ILogB(xh);
        int k = (int)Math.Round((double)e / n);
        DoubleDouble xs = k == 0 ? x : ScaleB(x, -k * n);
        if (double.IsSubnormal(xs._hi))
        {
            return Exp(Log(x) / (double)n);
        }
        DoubleDouble y = Math.Pow(xs._hi, 1.0 / n);
        int steps = n <= 5 ? 1 : 2;
        for (int i = 0; i < steps; i++)
        {
            DoubleDouble t = Pow(y, n);
            y += y * ((xs - t) / (t * (double)n));
        }
        return k == 0 ? y : ScaleB(y, k);
    }

    /// <summary>√(x² + y²) without intermediate overflow or underflow.</summary>
    public static DoubleDouble Hypot(DoubleDouble x, DoubleDouble y)
    {
        if (IsInfinity(x) || IsInfinity(y))
        {
            return PositiveInfinity;
        }
        if (IsNaN(x) || IsNaN(y))
        {
            return NaN;
        }
        DoubleDouble a = Abs(x), b = Abs(y);
        if (a < b)
        {
            (a, b) = (b, a);
        }
        if (b._hi == 0.0)
        {
            return a;
        }
        int e = Math.ILogB(a._hi);
        DoubleDouble sa = ScaleB(a, -e), sb = ScaleB(b, -e);
        return ScaleB(Sqrt(Square(sa) + Square(sb)), e);
    }
}
