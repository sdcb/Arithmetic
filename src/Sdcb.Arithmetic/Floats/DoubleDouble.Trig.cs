using Sdcb.Arithmetic.Core;

namespace Sdcb.Arithmetic;

public readonly partial struct DoubleDouble
{
    // ---------------------------------------------------------------
    // Kernels on |r| <= ~π/4 (z = r², |z| <= 0.62)
    // ---------------------------------------------------------------

    /// <summary>sin(r): r + r·z·S(z), Taylor to r^29; the terms from r^17 on are summed in double (below 2^-53 relative).</summary>
    private static DoubleDouble SinKernel(DoubleDouble r, DoubleDouble z)
    {
        if (r._hi == 0.0)
        {
            return r; // keeps the sign of zero
        }
        ReadOnlySpan<double> fh = InvFactHi, fl = InvFactLo;
        double zh = z._hi;
        double q = fh[29];
        q = Math.FusedMultiplyAdd(q, zh, -fh[27]);
        q = Math.FusedMultiplyAdd(q, zh, fh[25]);
        q = Math.FusedMultiplyAdd(q, zh, -fh[23]);
        q = Math.FusedMultiplyAdd(q, zh, fh[21]);
        q = Math.FusedMultiplyAdd(q, zh, -fh[19]);
        q = Math.FusedMultiplyAdd(q, zh, fh[17]);
        DoubleDouble p = new DoubleDouble(-fh[15], -fl[15], normalized: true) + z * q;
        p = new DoubleDouble(fh[13], fl[13], normalized: true) + z * p;
        p = new DoubleDouble(-fh[11], -fl[11], normalized: true) + z * p;
        p = new DoubleDouble(fh[9], fl[9], normalized: true) + z * p;
        p = new DoubleDouble(-fh[7], -fl[7], normalized: true) + z * p;
        p = new DoubleDouble(fh[5], fl[5], normalized: true) + z * p;
        p = new DoubleDouble(-fh[3], -fl[3], normalized: true) + z * p;
        return r + r * (z * p);
    }

    /// <summary>cos(r): 1 + z·C(z), Taylor to r^28; the terms from r^18 on are summed in double.</summary>
    private static DoubleDouble CosKernel(DoubleDouble z)
    {
        ReadOnlySpan<double> fh = InvFactHi, fl = InvFactLo;
        double zh = z._hi;
        double q = fh[28];
        q = Math.FusedMultiplyAdd(q, zh, -fh[26]);
        q = Math.FusedMultiplyAdd(q, zh, fh[24]);
        q = Math.FusedMultiplyAdd(q, zh, -fh[22]);
        q = Math.FusedMultiplyAdd(q, zh, fh[20]);
        q = Math.FusedMultiplyAdd(q, zh, -fh[18]);
        DoubleDouble p = new DoubleDouble(fh[16], fl[16], normalized: true) + z * q;
        p = new DoubleDouble(-fh[14], -fl[14], normalized: true) + z * p;
        p = new DoubleDouble(fh[12], fl[12], normalized: true) + z * p;
        p = new DoubleDouble(-fh[10], -fl[10], normalized: true) + z * p;
        p = new DoubleDouble(fh[8], fl[8], normalized: true) + z * p;
        p = new DoubleDouble(-fh[6], -fl[6], normalized: true) + z * p;
        p = new DoubleDouble(fh[4], fl[4], normalized: true) + z * p;
        p = new DoubleDouble(-0.5, 0.0, normalized: true) + z * p;
        return z * p + 1.0;
    }

    // ---------------------------------------------------------------
    // Argument reduction modulo π/2
    // ---------------------------------------------------------------

    /// <summary>
    /// r = x - q·π/2 with |r| &lt;= ~π/4 and returns q mod 4. Up to |q| = 2^20 a Cody-Waite
    /// reduction with a three-word π/2 is used; for larger arguments, or when r cancels below
    /// 2^-28 (x very close to a multiple of π/2), the reduction is done exactly with <see cref="BigFloat"/>.
    /// </summary>
    private static int ReduceHalfPi(DoubleDouble x, out DoubleDouble r)
    {
        double xh = x._hi;
        if (Math.Abs(xh) <= PiOver4_0)
        {
            r = x;
            return 0;
        }
        double q = Math.Round(xh * TwoOverPi0);
        if (Math.Abs(q) <= 1048576.0)
        {
            double p0 = ErrorFree.TwoProd(q, PiOver2_0, out double e0);
            double s = xh - p0; // exact (Sterbenz)
            r = new DoubleDouble(s, x._lo) - e0;
            double p1 = ErrorFree.TwoProd(q, PiOver2_1, out double e1);
            r += new DoubleDouble(-p1, -e1, normalized: true);
            r -= q * PiOver2_2;
            if (Math.Abs(r._hi) >= 3.725290298461914E-09)
            {
                return (int)((long)q & 3);
            }
        }
        return ReduceHalfPiExact(x, out r);
    }

    private static int ReduceHalfPiExact(DoubleDouble x, out DoubleDouble r)
    {
        int prec = 320 + Math.Max(0, Math.ILogB(x._hi));
        BigFloat xr = x.ToBigFloatExact();
        BigFloat halfPi = BigFloat.Pi(prec).ShiftRight(1);
        BigFloat q = BigFloat.Round(BigFloat.Divide(xr, halfPi, prec));
        BigFloat rr = BigFloat.Subtract(xr, BigFloat.Multiply(q, halfPi, prec + 64), prec);
        r = FromBigFloat(rr);
        BigFloat q4 = BigFloat.TruncatedRemainder(q, BigFloat.FromDouble(4.0, 8), 8);
        return ((int)q4.ToDouble() + 4) & 3;
    }

    // ---------------------------------------------------------------
    // Trigonometric functions
    // ---------------------------------------------------------------

    /// <summary>sin x.</summary>
    public static DoubleDouble Sin(DoubleDouble x)
    {
        if (!double.IsFinite(x._hi))
        {
            return NaN;
        }
        int q = ReduceHalfPi(x, out DoubleDouble r);
        DoubleDouble z = Square(r);
        return q switch
        {
            0 => SinKernel(r, z),
            1 => CosKernel(z),
            2 => -SinKernel(r, z),
            _ => -CosKernel(z),
        };
    }

    /// <summary>cos x.</summary>
    public static DoubleDouble Cos(DoubleDouble x)
    {
        if (!double.IsFinite(x._hi))
        {
            return NaN;
        }
        int q = ReduceHalfPi(x, out DoubleDouble r);
        DoubleDouble z = Square(r);
        return q switch
        {
            0 => CosKernel(z),
            1 => -SinKernel(r, z),
            2 => -CosKernel(z),
            _ => SinKernel(r, z),
        };
    }

    /// <summary>(sin x, cos x) with a single argument reduction.</summary>
    public static (DoubleDouble Sin, DoubleDouble Cos) SinCos(DoubleDouble x)
    {
        if (!double.IsFinite(x._hi))
        {
            return (NaN, NaN);
        }
        int q = ReduceHalfPi(x, out DoubleDouble r);
        DoubleDouble z = Square(r);
        DoubleDouble s = SinKernel(r, z), c = CosKernel(z);
        return q switch
        {
            0 => (s, c),
            1 => (c, -s),
            2 => (-s, -c),
            _ => (-c, s),
        };
    }

    /// <summary>tan x.</summary>
    public static DoubleDouble Tan(DoubleDouble x)
    {
        if (!double.IsFinite(x._hi))
        {
            return NaN;
        }
        int q = ReduceHalfPi(x, out DoubleDouble r);
        DoubleDouble z = Square(r);
        DoubleDouble s = SinKernel(r, z), c = CosKernel(z);
        return (q & 1) == 0 ? s / c : -c / s;
    }

    // ---------------------------------------------------------------
    // π-scaled trigonometric functions (exact reduction)
    // ---------------------------------------------------------------

    /// <summary>x = n/2 + r exactly with |r| &lt;= 1/4; returns n mod 4.</summary>
    private static int ReduceHalfInteger(DoubleDouble x, out DoubleDouble r)
    {
        // sin(πx) has period 2: reduce both words modulo 2 exactly first.
        double a = Math.Abs(x._hi) < 4.0 ? x._hi : x._hi % 2.0;
        double b = Math.Abs(x._lo) < 4.0 ? x._lo : x._lo % 2.0;
        double yh = ErrorFree.TwoSum(a, b, out double yl);
        double n = Math.Round(yh * 2.0);
        double d = yh - n * 0.5; // exact
        double rh = ErrorFree.TwoSum(d, yl, out double rl);
        r = new DoubleDouble(rh, rl, normalized: true);
        return (int)n & 3;
    }

    /// <summary>sin(πx), with exact argument reduction (exactly 0 at integers).</summary>
    public static DoubleDouble SinPi(DoubleDouble x)
    {
        if (!double.IsFinite(x._hi))
        {
            return NaN;
        }
        int q = ReduceHalfInteger(x, out DoubleDouble r);
        DoubleDouble t = r * Pi;
        DoubleDouble z = Square(t);
        DoubleDouble v = q switch
        {
            0 => SinKernel(t, z),
            1 => CosKernel(z),
            2 => -SinKernel(t, z),
            _ => -CosKernel(z),
        };
        return v._hi == 0.0 ? CopySign(Zero, x) : v;
    }

    /// <summary>cos(πx), with exact argument reduction (exactly 0 at half-integers).</summary>
    public static DoubleDouble CosPi(DoubleDouble x)
    {
        if (!double.IsFinite(x._hi))
        {
            return NaN;
        }
        int q = ReduceHalfInteger(x, out DoubleDouble r);
        DoubleDouble t = r * Pi;
        DoubleDouble z = Square(t);
        DoubleDouble v = q switch
        {
            0 => CosKernel(z),
            1 => -SinKernel(t, z),
            2 => -CosKernel(z),
            _ => SinKernel(t, z),
        };
        return v._hi == 0.0 ? Zero : v;
    }

    /// <summary>(sin(πx), cos(πx)).</summary>
    public static (DoubleDouble SinPi, DoubleDouble CosPi) SinCosPi(DoubleDouble x) => (SinPi(x), CosPi(x));

    /// <summary>tan(πx); +∞ at 1/2 + 2k and -∞ at 3/2 + 2k.</summary>
    public static DoubleDouble TanPi(DoubleDouble x)
    {
        if (!double.IsFinite(x._hi))
        {
            return NaN;
        }
        int q = ReduceHalfInteger(x, out DoubleDouble r);
        if (r._hi == 0.0)
        {
            return q switch
            {
                1 => PositiveInfinity,
                3 => NegativeInfinity,
                _ => CopySign(Zero, x),
            };
        }
        DoubleDouble t = r * Pi;
        DoubleDouble z = Square(t);
        DoubleDouble s = SinKernel(t, z), c = CosKernel(z);
        return (q & 1) == 0 ? s / c : -c / s;
    }

    // ---------------------------------------------------------------
    // Inverse trigonometric functions (Newton on sin/cos from the double estimate)
    // ---------------------------------------------------------------

    private static DoubleDouble PiOver2 => new(PiOver2_0, PiOver2_1, normalized: true);
    private static DoubleDouble PiOver4 => new(PiOver4_0, PiOver4_1, normalized: true);
    private static DoubleDouble ThreePiOver4 => new(ThreePiOver4_0, ThreePiOver4_1, normalized: true);
    private static DoubleDouble InvPi => new(InvPi0, InvPi1, normalized: true);

    /// <summary>atan2(y, x) in [-π, π], with the IEEE special cases of <see cref="Math.Atan2(double, double)"/>.</summary>
    public static DoubleDouble Atan2(DoubleDouble y, DoubleDouble x)
    {
        double yh = y._hi, xh = x._hi;
        if (double.IsNaN(yh) || double.IsNaN(xh))
        {
            return NaN;
        }
        if (yh == 0.0 || double.IsInfinity(xh) || xh == 0.0 || double.IsInfinity(yh))
        {
            return Atan2Special(y, x);
        }
        if (x == y)
        {
            return yh > 0 ? PiOver4 : -ThreePiOver4;
        }
        if (x == -y)
        {
            return yh > 0 ? ThreePiOver4 : -PiOver4;
        }
        // Common power-of-two scaling keeps the squares in range; atan2 is invariant under it.
        int e = Math.Max(Math.ILogB(xh), Math.ILogB(yh));
        DoubleDouble xs = ScaleB(x, -e), ys = ScaleB(y, -e);
        DoubleDouble rad = Sqrt(Square(xs) + Square(ys));
        DoubleDouble xx = xs / rad, yy = ys / rad;
        DoubleDouble z = Math.Atan2(yh, xh);
        (DoubleDouble s, DoubleDouble c) = SinCos(z);
        if (Math.Abs(xx._hi) > Math.Abs(yy._hi))
        {
            z += (yy - s) / c;
        }
        else
        {
            z -= (xx - c) / s;
        }
        return z;
    }

    private static DoubleDouble Atan2Special(DoubleDouble y, DoubleDouble x)
    {
        double yh = y._hi, xh = x._hi;
        bool yNeg = double.IsNegative(yh);
        DoubleDouble Signed(DoubleDouble v) => yNeg ? -v : v;
        if (double.IsInfinity(yh))
        {
            if (double.IsInfinity(xh))
            {
                return Signed(xh > 0 ? PiOver4 : ThreePiOver4);
            }
            return Signed(PiOver2);
        }
        if (yh == 0.0)
        {
            return double.IsNegative(xh) ? Signed(Pi) : y;
        }
        if (xh == 0.0)
        {
            return Signed(PiOver2);
        }
        // x = ±∞, y finite nonzero
        return xh > 0 ? Signed(Zero) : Signed(Pi);
    }

    /// <summary>atan x.</summary>
    public static DoubleDouble Atan(DoubleDouble x)
    {
        if (double.IsInfinity(x._hi))
        {
            return x._hi > 0 ? PiOver2 : -PiOver2;
        }
        return x._hi == 0.0 ? x : Atan2(x, One);
    }

    /// <summary>asin x (NaN outside [-1, 1]).</summary>
    public static DoubleDouble Asin(DoubleDouble x)
    {
        DoubleDouble a = Abs(x);
        if (!(a <= One))
        {
            return NaN;
        }
        if (a == One)
        {
            return x._hi > 0 ? PiOver2 : -PiOver2;
        }
        if (x._hi == 0.0)
        {
            return x;
        }
        return Atan2(x, Sqrt((One - x) * (x + 1.0)));
    }

    /// <summary>acos x (NaN outside [-1, 1]).</summary>
    public static DoubleDouble Acos(DoubleDouble x)
    {
        DoubleDouble a = Abs(x);
        if (!(a <= One))
        {
            return NaN;
        }
        if (a == One)
        {
            return x._hi > 0 ? Zero : Pi;
        }
        if (x._hi == 0.0)
        {
            return PiOver2;
        }
        return Atan2(Sqrt((One - x) * (x + 1.0)), x);
    }

    /// <summary>asin(x)/π.</summary>
    public static DoubleDouble AsinPi(DoubleDouble x)
    {
        if (Abs(x) == One)
        {
            return x._hi > 0 ? 0.5 : -0.5;
        }
        return Asin(x) * InvPi;
    }

    /// <summary>acos(x)/π.</summary>
    public static DoubleDouble AcosPi(DoubleDouble x)
    {
        if (x == NegativeOne)
        {
            return One;
        }
        if (x._hi == 0.0)
        {
            return 0.5;
        }
        return Acos(x) * InvPi;
    }

    /// <summary>atan(x)/π.</summary>
    public static DoubleDouble AtanPi(DoubleDouble x)
    {
        if (double.IsInfinity(x._hi))
        {
            return x._hi > 0 ? 0.5 : -0.5;
        }
        if (Abs(x) == One)
        {
            return x._hi > 0 ? 0.25 : -0.25;
        }
        return Atan(x) * InvPi;
    }

    /// <summary>atan2(y, x)/π.</summary>
    public static DoubleDouble Atan2Pi(DoubleDouble y, DoubleDouble x)
    {
        DoubleDouble a = Atan2(y, x);
        if (a == Pi || a == -Pi)
        {
            return a._hi > 0 ? One : NegativeOne;
        }
        if (a == PiOver2 || a == -PiOver2)
        {
            return a._hi > 0 ? 0.5 : -0.5;
        }
        if (a == PiOver4 || a == -PiOver4)
        {
            return a._hi > 0 ? 0.25 : -0.25;
        }
        if (a == ThreePiOver4 || a == -ThreePiOver4)
        {
            return a._hi > 0 ? 0.75 : -0.75;
        }
        return a * InvPi;
    }

    /// <summary>degrees · π/180.</summary>
    public static DoubleDouble DegreesToRadians(DoubleDouble degrees) => degrees * new DoubleDouble(DegToRad0, DegToRad1, normalized: true);

    /// <summary>radians · 180/π.</summary>
    public static DoubleDouble RadiansToDegrees(DoubleDouble radians) => radians * new DoubleDouble(RadToDeg0, RadToDeg1, normalized: true);

    // ---------------------------------------------------------------
    // Hyperbolic functions
    // ---------------------------------------------------------------

    /// <summary>sinh x (via expm1 for accuracy near 0).</summary>
    public static DoubleDouble Sinh(DoubleDouble x)
    {
        if (!double.IsFinite(x._hi) || x._hi == 0.0)
        {
            return x;
        }
        DoubleDouble a = Abs(x);
        DoubleDouble r;
        if (a._hi > 40.0)
        {
            r = HalfExp(a); // e^-a is below 2^-115 relative
        }
        else
        {
            DoubleDouble em = ExpM1(a);
            r = MulPow2(em + em / (em + 1.0), 0.5);
        }
        return double.IsNegative(x._hi) ? -r : r;
    }

    /// <summary>e^a / 2 for a &gt; 0, finite up to the sinh/cosh overflow threshold (~710.48).</summary>
    private static DoubleDouble HalfExp(DoubleDouble a) => a._hi > 710.5 ? PositiveInfinity : ExpCore(a, -1);

    /// <summary>cosh x.</summary>
    public static DoubleDouble Cosh(DoubleDouble x)
    {
        if (double.IsNaN(x._hi))
        {
            return x;
        }
        DoubleDouble a = Abs(x);
        if (a._hi > 40.0)
        {
            return HalfExp(a);
        }
        DoubleDouble e = Exp(a);
        return MulPow2(e + One / e, 0.5);
    }

    /// <summary>tanh x.</summary>
    public static DoubleDouble Tanh(DoubleDouble x)
    {
        if (double.IsNaN(x._hi) || x._hi == 0.0)
        {
            return x;
        }
        DoubleDouble a = Abs(x);
        DoubleDouble r;
        if (a._hi < 0.5)
        {
            DoubleDouble em = ExpM1(MulPow2(a, 2.0));
            r = em / (em + 2.0);
        }
        else if (a._hi > 40.0)
        {
            r = One - MulPow2(Exp(-MulPow2(a, 2.0)), 2.0);
        }
        else
        {
            r = One - 2.0 / (Exp(MulPow2(a, 2.0)) + 1.0);
        }
        return double.IsNegative(x._hi) ? -r : r;
    }

    /// <summary>asinh x.</summary>
    public static DoubleDouble Asinh(DoubleDouble x)
    {
        if (!double.IsFinite(x._hi) || x._hi == 0.0)
        {
            return x;
        }
        DoubleDouble a = Abs(x);
        DoubleDouble r;
        if (a._hi > 1.152921504606847E18) // 2^60: sqrt(a² + 1) = a to beyond 2^-120
        {
            r = Log(a) + Ln2;
        }
        else if (a._hi < 0.5)
        {
            DoubleDouble a2 = Square(a);
            r = LogP1(a + a2 / (Sqrt(a2 + 1.0) + 1.0));
        }
        else
        {
            r = Log(a + Sqrt(Square(a) + 1.0));
        }
        return double.IsNegative(x._hi) ? -r : r;
    }

    /// <summary>acosh x (NaN below 1).</summary>
    public static DoubleDouble Acosh(DoubleDouble x)
    {
        if (!(x >= One))
        {
            return NaN;
        }
        if (double.IsInfinity(x._hi))
        {
            return x;
        }
        if (x._hi > 1.152921504606847E18)
        {
            return Log(x) + Ln2;
        }
        DoubleDouble t = x - 1.0; // exact
        return LogP1(t + Sqrt(t * (t + 2.0)));
    }

    /// <summary>atanh x (±∞ at ±1, NaN beyond).</summary>
    public static DoubleDouble Atanh(DoubleDouble x)
    {
        if (x._hi == 0.0 || double.IsNaN(x._hi))
        {
            return x;
        }
        DoubleDouble a = Abs(x);
        if (!(a <= One))
        {
            return NaN;
        }
        DoubleDouble r = a == One ? PositiveInfinity : MulPow2(LogP1(MulPow2(a, 2.0) / (One - a)), 0.5);
        return double.IsNegative(x._hi) ? -r : r;
    }
}
