using System.Runtime.CompilerServices;
using static Sdcb.Arithmetic.Core.ErrorFree;

namespace Sdcb.Arithmetic;

// Ziv-style fast path for the elementary functions (see Float128.FastPath.Kernels.cs).
//
// A kernel evaluates f in triple-double arithmetic (three non-overlapping doubles, ~159 bits) to an
// approximation y with a proven relative error bound ε, and TryRound returns the binary128 rounding of y
// only when every value of [y(1-ε), y(1+ε)] rounds to the same binary128 value under the requested mode
// (rounding is monotonic, so f's correct rounding is then that value). Otherwise, and whenever an argument
// or result leaves a kernel's proven range (subnormal or overflowing results, huge trig arguments, special
// values), the caller falls back to the correctly rounded BigFloat path. Exact results sit on a rounding
// boundary in the directed modes and RNE midpoints sit on one in RNE, so both always fall back.
//
// Error bounds of the primitives below are relative to the exact result of the operation on the given
// operands. They assume normalized operands: |M| ≤ 2^-53|H| and |L| ≤ 2^-104|H| (every primitive
// renormalizes its result with error-free TwoSums, which keeps this invariant).
public readonly partial struct Float128
{
    /// <summary>A triple-double H + M + L (non-overlapping, decreasing magnitudes).</summary>
    private readonly struct Td
    {
        public readonly double H, M, L;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Td(double h, double m, double l)
        {
            H = h;
            M = m;
            L = l;
        }

        public static Td operator -(Td a) => new(-a.H, -a.M, -a.L);
    }

    // ---------------------------------------------------------------
    // Triple-double primitives
    // ---------------------------------------------------------------

    /// <summary>Error-free: returns (h, m, l) with h + m + l = a + b + c exactly and |m| ≤ ulp(h)/2.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static Td TdRenorm(double a, double b, double c)
    {
        double s = TwoSum(b, c, out double t);
        double h = TwoSum(a, s, out double u);
        double m = TwoSum(u, t, out double l);
        h = TwoSum(h, m, out m);
        return new Td(h, m, l);
    }

    /// <summary>
    /// a + b + c + d (decreasing magnitude levels, arbitrary cancellation in a + b): the four values are
    /// summed with TwoSums; the only rounding is in the last component, so the error is at most
    /// 2^-106·|c| + 2^-158·|result|.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static Td TdRenorm4(double a, double b, double c, double d)
    {
        double s = TwoSum(c, d, out double t3);
        s = TwoSum(b, s, out double t2);
        double h = TwoSum(a, s, out double t1);
        double m = TwoSum(t1, t2, out double u);
        double l = u + t3;
        h = TwoSum(h, m, out m);
        return new Td(h, m, l);
    }

    /// <summary>a + b; error ≤ 2^-153·(|a| + |b|): only the third-level sum a.L + b.L + e1 + e2 is rounded.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static Td TdAdd(in Td a, in Td b)
    {
        double s0 = TwoSum(a.H, b.H, out double e0);
        double s1 = TwoSum(a.M, b.M, out double e1);
        s1 = TwoSum(s1, e0, out double e2);
        double s2 = a.L + b.L + e1 + e2;
        return TdRenorm(s0, s1, s2);
    }

    /// <summary>a - b; error ≤ 2^-153·(|a| + |b|).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static Td TdSub(in Td a, in Td b) => TdAdd(a, -b);

    /// <summary>a + (bh + bl) for a double-double; error ≤ 2^-153·(|a| + |b|).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static Td TdAdd2(in Td a, double bh, double bl)
    {
        double s0 = TwoSum(a.H, bh, out double e0);
        double s1 = TwoSum(a.M, bl, out double e1);
        s1 = TwoSum(s1, e0, out double e2);
        double s2 = a.L + e1 + e2;
        return TdRenorm(s0, s1, s2);
    }

    /// <summary>a + b for a double b; error ≤ 2^-155·(|a| + |b|) (exact when a + b fits in the three components).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static Td TdAddD(in Td a, double b)
    {
        double s0 = TwoSum(a.H, b, out double e0);
        double s1 = TwoSum(a.M, e0, out double e1);
        double s2 = a.L + e1;
        return TdRenorm(s0, s1, s2);
    }

    /// <summary>a·d for a double d; error ≤ 2^-155·|a·d| (the dropped products are all formed; only a.L·d + q1 + e is rounded).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static Td TdMulD(in Td a, double d)
    {
        double p0 = TwoProd(a.H, d, out double q0);
        double p1 = TwoProd(a.M, d, out double q1);
        double s1 = TwoSum(q0, p1, out double e);
        double s2 = Math.FusedMultiplyAdd(a.L, d, q1 + e);
        return TdRenorm(p0, s1, s2);
    }

    /// <summary>
    /// a·b; error ≤ 2^-150·|a·b|. HH, HM and MH are formed exactly (TwoProd); the order-2^-104 terms
    /// (their errors, HL, MM, LH) are summed in double (≤ 2^-152 after six roundings) and ML, LM, LL
    /// (≤ 2^-156) are dropped.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static Td TdMul(in Td a, in Td b)
    {
        double p0 = TwoProd(a.H, b.H, out double q0);
        double p1 = TwoProd(a.H, b.M, out double q1);
        double p2 = TwoProd(a.M, b.H, out double q2);
        double s1 = TwoSum(p1, p2, out double e1);
        s1 = TwoSum(s1, q0, out double e2);
        double s2 = Math.FusedMultiplyAdd(a.H, b.L, Math.FusedMultiplyAdd(a.M, b.M, Math.FusedMultiplyAdd(a.L, b.H, q1 + q2 + e1 + e2)));
        return TdRenorm(p0, s1, s2);
    }

    /// <summary>a²; error ≤ 2^-150·a² (as <see cref="TdMul"/>).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static Td TdSqr(in Td a)
    {
        double p0 = TwoSquare(a.H, out double q0);
        double h2 = 2.0 * a.H;
        double p1 = TwoProd(h2, a.M, out double q1);
        double s1 = TwoSum(p1, q0, out double e1);
        double s2 = Math.FusedMultiplyAdd(h2, a.L, Math.FusedMultiplyAdd(a.M, a.M, q1 + e1));
        return TdRenorm(p0, s1, s2);
    }

    /// <summary>
    /// a/b by long division with three quotient digits; error ≤ 2^-150·|a/b|. Each digit q = r.H·(1/b.H)
    /// has relative error ≤ 2^-51.4, so the remainders shrink to 2^-51.4|a| and 2^-102.8|a|; the first
    /// remainder carries the TdMulD/TdAdd rounding (≤ 2^-152.5|a|), the last digit errs by ≤ 2^-154|a/b|.
    /// </summary>
    private static Td TdDiv(in Td a, in Td b)
    {
        double inv = 1.0 / b.H;
        double q0 = a.H * inv;
        Td r = TdSub(a, TdMulD(b, q0));
        double q1 = r.H * inv;
        r = TdSub(r, TdMulD(b, q1));
        double q2 = r.H * inv;
        return TdRenorm(q0, q1, q2);
    }

    /// <summary>
    /// √a for a &gt; 0: s = √a.H, d = (a - s²)/(2s) in double-double, √a = s + d - d²/(2s) + O(2^-156·s);
    /// error ≤ 2^-150·√a (the residual's TdAdd rounding contributes 2^-153.9·s, the rest ≤ 2^-155·s).
    /// </summary>
    private static Td TdSqrt(in Td a)
    {
        double s = Math.Sqrt(a.H);
        double sq = TwoSquare(s, out double sqe);
        Td r = TdAdd(a, new Td(-sq, -sqe, 0));
        double inv2s = 0.5 / s;
        double d0 = r.H * inv2s;
        double pp = TwoProd(d0, 2.0 * s, out double pe);
        double d1 = (((r.H - pp) - pe) + r.M) * inv2s;
        double c = -d0 * d0 * inv2s;
        return TdRenorm(s, d0, d1 + c);
    }

    /// <summary>a·2^n (exact for the ranges used here).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static Td TdScale(in Td a, int n)
    {
        double f = Pow2(n);
        return new Td(a.H * f, a.M * f, a.L * f);
    }

    // ---------------------------------------------------------------
    // Double-double helpers (error ≤ 2^-104 relative for DdMul; DdAdd ≤ 2^-104·(|a| + |b|))
    // ---------------------------------------------------------------

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double DdMul(double ah, double al, double bh, double bl, out double lo)
    {
        double p = TwoProd(ah, bh, out double e);
        e = Math.FusedMultiplyAdd(ah, bl, Math.FusedMultiplyAdd(al, bh, e));
        return FastTwoSum(p, e, out lo);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double DdAdd(double ah, double al, double bh, double bl, out double lo)
    {
        double s = TwoSum(ah, bh, out double e);
        e += al + bl;
        return FastTwoSum(s, e, out lo);
    }

    /// <summary>One double-double Horner step: c + x·v.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double DdHorner(double ch, double cl, double xh, double xl, double vh, double vl, out double lo)
    {
        double mh = DdMul(xh, xl, vh, vl, out double ml);
        return DdAdd(ch, cl, mh, ml, out lo);
    }

    /// <summary>2^n as a double for n in [-1022, 1023].</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double Pow2(int n) => BitConverter.Int64BitsToDouble((long)(n + 1023) << 52);

    // ---------------------------------------------------------------
    // Conversions
    // ---------------------------------------------------------------

    /// <summary>The 113-bit significand of a normal value (hidden bit included).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static UInt128 Significand(Float128 x) => (x._bits & FracMask) | Hidden;

    /// <summary>
    /// A 113-bit significand sig times 2^(e - 112) as an exact triple-double: the top 53 bits, the next 53 and
    /// the last 7 (requires e - 112 ≥ -1022 and e ≤ 1023).
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static Td SigToTd(UInt128 sig, int e, bool negative)
    {
        ulong hi = (ulong)(sig >> 60);
        ulong mid = (ulong)(sig >> 7) & ((1UL << 53) - 1);
        ulong lo = (ulong)sig & 127;
        double h = hi * Pow2(e - 52), m = mid * Pow2(e - 105), l = lo * Pow2(e - 112);
        return negative ? new Td(-h, -m, -l) : new Td(h, m, l);
    }

    /// <summary>x·2^shift as an exact triple-double; x must be normal with -910 ≤ ILogB(x) + shift ≤ 1023.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static Td ToTd(Float128 x, int shift = 0) => SigToTd(Significand(x), x.BiasedExponent - ExponentBias + shift, x.SignBit);

    // ---------------------------------------------------------------
    // Rounding test
    // ---------------------------------------------------------------

    /// <summary>
    /// Rounds y·2^scale to binary128 when the result is decided by y's error bound: f = y(1 + θ) with
    /// |θ| ≤ eps. Returns false (caller falls back) when the interval contains a rounding boundary, when y is
    /// zero or not finite, or when the result would be subnormal or within one ulp of overflow.
    /// </summary>
    /// <remarks>
    /// With F = 2^(115 - e) (e the exponent of y.H after normalization) y·F = A + t exactly up to the
    /// rounding of t: A = y.H·F + floor(y.M·F) is an integer in [2^115 - 2^62, 2^116 + 2^62] and
    /// t = frac(y.M·F) + y.L·F a double (|t| ≲ 2^12, rounding ≤ 2^-40). The binary128 significand is the top
    /// 113 bits of A; the tail T (bits below them plus t, in units of 2^E) must stay at distance &gt; δ from the
    /// rounding boundaries, δ = eps·2^118 + 2^-36 ≥ eps·y·F plus every rounding of t and T.
    /// </remarks>
    private static bool TryRound(in Td y, int scale, double eps, RoundingMode rounding, out Float128 result)
    {
        result = default;
        double h = y.H, m = y.M, l = y.L;
        if (!(Math.Abs(h) < double.PositiveInfinity) || h == 0 || !(eps < 2.0E-30))
        {
            return false;
        }
        bool neg = h < 0;
        if (neg)
        {
            h = -h;
            m = -m;
            l = -l;
        }
        h = FastTwoSum(h, m, out m); // |m| ≤ ulp(h)/2
        ulong hb = BitConverter.DoubleToUInt64Bits(h);
        int eh = (int)(hb >> 52) - 1023;
        if (eh < -900 || eh > 900)
        {
            return false;
        }
        double f = Pow2(115 - eh);
        double mf = m * f; // exact, |mf| ≤ 2^62
        double mi = Math.Floor(mf);
        double t = (mf - mi) + l * f;
        Int128 a = ((Int128)((hb & 0xF_FFFF_FFFF_FFFFUL) | (1UL << 52)) << 63) + (long)mi;
        int bl = 128 - (int)Int128.LeadingZeroCount(a);
        int s = bl - 113; // 2..4
        Int128 q = a >> s;
        double w = s switch { 2 => 4.0, 3 => 8.0, _ => 16.0 };
        double tail = (double)(long)(a & ((Int128.One << s) - 1)) + t;
        if (!(tail >= 0 && tail < w))
        {
            double k = Math.Floor(tail / w);
            if (!(Math.Abs(k) < 1e6))
            {
                return false;
            }
            q += (long)k;
            tail -= k * w;
        }
        if (q < (Int128.One << 112) || q >= (Int128.One << 113))
        {
            return false; // the tail moved the value across a binade
        }
        double delta = eps * 3.3230699894622794E+35 /* 2^118 */ + 1.4551915228366852E-11 /* 2^-36 */;
        bool up;
        if (rounding == RoundingMode.NearestEven)
        {
            double half = 0.5 * w;
            if (tail + delta < half)
            {
                up = false;
            }
            else if (tail - delta > half)
            {
                up = true;
            }
            else
            {
                return false;
            }
        }
        else
        {
            if (!(tail - delta > 0) || !(tail + delta < w))
            {
                return false;
            }
            up = DirectedUp(neg, rounding);
        }
        // value = q·2^(s - 115 + eh + scale); binary128 with a 113-bit significand q: biased exponent below.
        long be = (long)s - 115 + eh + scale + 16495;
        if (be < 1 || be > 0x7FFD)
        {
            return false;
        }
        UInt128 bits = ((UInt128)(ulong)(be - 1) << 112) + (UInt128)q + (up ? UInt128.One : UInt128.Zero);
        result = new(bits | SignBits(neg));
        return true;
    }

    // ---------------------------------------------------------------
    // Constants and tables (computed once from BigFloat at 320 bits)
    // ---------------------------------------------------------------

    private const int TablePrecision = 320;

    /// <summary>The nearest triple-double to v (error ≤ 2^-160·|v|).</summary>
    private static Td TdOf(BigFloat v)
    {
        double h = v.ToDouble();
        BigFloat r = BigFloat.Subtract(v, BigFloat.FromDouble(h), 1000);
        double m = r.ToDouble();
        r = BigFloat.Subtract(r, BigFloat.FromDouble(m), 1000);
        return new Td(h, m, r.ToDouble());
    }

    /// <summary>v as n doubles of decreasing magnitude (each the nearest double to the remainder).</summary>
    private static double[] PartsOf(BigFloat v, int n)
    {
        var p = new double[n];
        for (int i = 0; i < n; i++)
        {
            p[i] = v.ToDouble();
            v = BigFloat.Subtract(v, BigFloat.FromDouble(p[i]), 1000);
        }
        return p;
    }

    private static BigFloat Rat(double num, double den) => BigFloat.Divide(BigFloat.FromDouble(num), BigFloat.FromDouble(den), TablePrecision);

    /// <summary>Shared constants.</summary>
    private static class FpConst
    {
        internal static readonly double Ln2P0, Ln2P1, Ln2P2;
        internal static readonly Td Ln2, InvLn2, InvLn10, Pi, PiHalf, OneThird, OneSixth;

        static FpConst()
        {
            BigFloat ln2 = BigFloat.Ln2(TablePrecision);
            double[] p = PartsOf(ln2, 3);
            (Ln2P0, Ln2P1, Ln2P2) = (p[0], p[1], p[2]);
            Ln2 = TdOf(ln2);
            InvLn2 = TdOf(BigFloat.Divide(BigFloat.FromDouble(1), ln2, TablePrecision));
            InvLn10 = TdOf(BigFloat.Divide(BigFloat.FromDouble(1), BigFloat.Log(BigFloat.FromDouble(10), TablePrecision), TablePrecision));
            BigFloat pi = BigFloat.Pi(TablePrecision);
            Pi = TdOf(pi);
            PiHalf = TdOf(pi.ShiftRight(1));
            OneThird = TdOf(Rat(1, 3));
            OneSixth = TdOf(Rat(1, 6));
        }
    }

    /// <summary>Exponential: ln2/4096 in four parts, 2^(j/64) and 2^(j/4096) for j = 0..63, 1/k! coefficients.</summary>
    private static class FpExp
    {
        internal const double K = 4096 / 0.6931471805599453; // ≈ 4096/ln2 (only selects N)
        internal static readonly double L0, L1, L2, L3;
        internal static readonly Td[] T1 = new Td[64], T2 = new Td[64];
        // 1/k! for k = 3..6 as double-doubles, k = 7..9 as doubles.
        internal static readonly double F3h, F3l, F4h, F4l, F5h, F5l, F6h, F6l, F7, F8, F9;

        static FpExp()
        {
            double[] l = PartsOf(BigFloat.Ln2(TablePrecision).ShiftRight(12), 4);
            (L0, L1, L2, L3) = (l[0], l[1], l[2], l[3]);
            for (int j = 0; j < 64; j++)
            {
                T1[j] = TdOf(BigFloat.Exp2(BigFloat.FromDouble(j / 64.0), TablePrecision));
                T2[j] = TdOf(BigFloat.Exp2(BigFloat.FromDouble(j / 4096.0), TablePrecision));
            }
            double[] f3 = PartsOf(Rat(1, 6), 2), f4 = PartsOf(Rat(1, 24), 2), f5 = PartsOf(Rat(1, 120), 2), f6 = PartsOf(Rat(1, 720), 2);
            (F3h, F3l, F4h, F4l, F5h, F5l, F6h, F6l) = (f3[0], f3[1], f4[0], f4[1], f5[0], f5[1], f6[0], f6[1]);
            F7 = 1.0 / 5040;
            F8 = 1.0 / 40320;
            F9 = 1.0 / 362880;
        }
    }

    /// <summary>
    /// Logarithm: step 1 uses c1 = RN(128/v) for v = round(128m) ∈ [90, 182] with -log c1; step 2 uses
    /// c2 = RN(1/(1 + v2·2^-14)) for v2 ∈ [-96, 96] with -log c2; log1p coefficients.
    /// </summary>
    private static class FpLog
    {
        internal static readonly double[] C1 = new double[93], C2 = new double[193];
        internal static readonly Td[] NL1 = new Td[93], NL2 = new Td[193];
        /// <summary>floor(√2·2^112): significands above it are halved so m ∈ [√½, √2).</summary>
        internal static readonly UInt128 Sqrt2Sig;
        // Coefficients of P(z) = 1/3 - z/4 + z²/5 - z³/6 (double-double) + z⁴(1/7 - z/8 + z²/9 - z³/10) (double).
        internal static readonly double P3h, P3l, P4h, P4l, P5h, P5l, P6h, P6l;
        internal const double P7 = 1.0 / 7, P8 = -1.0 / 8, P9 = 1.0 / 9, P10 = -1.0 / 10;

        static FpLog()
        {
            for (int v = 90; v <= 182; v++)
            {
                double c = 128.0 / v;
                C1[v - 90] = c;
                NL1[v - 90] = TdOf(-BigFloat.Log(BigFloat.FromDouble(c), TablePrecision));
            }
            for (int v = -96; v <= 96; v++)
            {
                double c = 1.0 / (1.0 + v / 16384.0);
                C2[v + 96] = c;
                NL2[v + 96] = TdOf(-BigFloat.Log(BigFloat.FromDouble(c), TablePrecision));
            }
            BigFloat s2 = BigFloat.Sqrt(BigFloat.FromDouble(2), 200).ShiftLeft(112);
            Sqrt2Sig = (UInt128)BigFloat.Floor(s2).ToBigInt();
            double[] a = PartsOf(Rat(1, 3), 2), b = PartsOf(Rat(-1, 4), 2), c5 = PartsOf(Rat(1, 5), 2), d = PartsOf(Rat(-1, 6), 2);
            (P3h, P3l, P4h, P4l, P5h, P5l, P6h, P6l) = (a[0], a[1], b[0], b[1], c5[0], c5[1], d[0], d[1]);
        }
    }

    /// <summary>Trigonometric: π/256 in four parts, sin(jπ/256) and cos(jπ/256) for j = 0..127, series coefficients.</summary>
    private static class FpTrig
    {
        internal const double K = 81.48733086305042; // ≈ 256/π (only selects N)
        internal static readonly double P0, P1, P2, P3;
        internal static readonly Td[] S = new Td[128], C = new Td[128];
        // sin: S5 = 1/5! - r²/7! + r⁴/9! (double-double) + r⁶(-1/11! + r²/13! - r⁴/15!) (double)
        internal static readonly double S5h, S5l, S7h, S7l, S9h, S9l, S11, S13, S15;
        // cos: C4 = 1/4! - r²/6! + r⁴/8! (double-double) + r⁶(-1/10! + r²/12! - r⁴/14!) (double)
        internal static readonly double C4h, C4l, C6h, C6l, C8h, C8l, C10, C12, C14;

        static FpTrig()
        {
            BigFloat pi = BigFloat.Pi(TablePrecision + 20);
            double[] p = PartsOf(pi.ShiftRight(8), 4);
            (P0, P1, P2, P3) = (p[0], p[1], p[2], p[3]);
            for (int j = 0; j < 128; j++)
            {
                BigFloat a = BigFloat.Multiply(pi, BigFloat.FromDouble(j / 256.0), TablePrecision + 20);
                S[j] = TdOf(BigFloat.Sin(a, TablePrecision));
                C[j] = TdOf(BigFloat.Cos(a, TablePrecision));
            }
            double[] s5 = PartsOf(Rat(1, 120), 2), s7 = PartsOf(Rat(-1, 5040), 2), s9 = PartsOf(Rat(1, 362880), 2);
            (S5h, S5l, S7h, S7l, S9h, S9l) = (s5[0], s5[1], s7[0], s7[1], s9[0], s9[1]);
            S11 = -1.0 / 39916800;
            S13 = 1.0 / 6227020800;
            S15 = -1.0 / 1307674368000;
            double[] c4 = PartsOf(Rat(1, 24), 2), c6 = PartsOf(Rat(-1, 720), 2), c8 = PartsOf(Rat(1, 40320), 2);
            (C4h, C4l, C6h, C6l, C8h, C8l) = (c4[0], c4[1], c6[0], c6[1], c8[0], c8[1]);
            C10 = -1.0 / 3628800;
            C12 = 1.0 / 479001600;
            C14 = -1.0 / 87178291200;
        }
    }

    /// <summary>Arctangent: atan(j/128) for j = 0..128 and series coefficients.</summary>
    private static class FpAtan
    {
        internal static readonly Td[] A = new Td[129];
        // A5 = 1/5 - u²/7 + u⁴/9 - u⁶/11 (double-double) + u⁸(1/13 - u²/15 + u⁴/17 - u⁶/19) (double)
        internal static readonly double A5h, A5l, A7h, A7l, A9h, A9l, A11h, A11l;
        internal const double A13 = 1.0 / 13, A15 = -1.0 / 15, A17 = 1.0 / 17, A19 = -1.0 / 19;

        static FpAtan()
        {
            for (int j = 0; j <= 128; j++)
            {
                A[j] = TdOf(BigFloat.Atan(BigFloat.FromDouble(j / 128.0), TablePrecision));
            }
            double[] a5 = PartsOf(Rat(1, 5), 2), a7 = PartsOf(Rat(-1, 7), 2), a9 = PartsOf(Rat(1, 9), 2), a11 = PartsOf(Rat(-1, 11), 2);
            (A5h, A5l, A7h, A7l, A9h, A9l, A11h, A11l) = (a5[0], a5[1], a7[0], a7[1], a9[0], a9[1], a11[0], a11[1]);
        }
    }
}
