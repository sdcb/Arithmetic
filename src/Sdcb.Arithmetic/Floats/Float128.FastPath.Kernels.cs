using System.Runtime.CompilerServices;
using static Sdcb.Arithmetic.Core.ErrorFree;

namespace Sdcb.Arithmetic;

// Triple-double kernels of the fast path (see Float128.FastPath.cs for the primitives, their error bounds
// and the rounding test). Each K* method returns y and a relative bound eps with f = y·2^scale·(1 + θ),
// |θ| ≤ eps, or false when the argument is outside the range its bound was derived for. The Eps* constants
// are the bounds used in the rounding test: each is the derived bound (stated next to it) times a safety
// factor of at least 2^6, and at least 2^8 times the largest error observed against Real at 300 bits
// (Float128FastPathTests.KernelErrorsStayFarBelowTheirBounds).
public readonly partial struct Float128
{
    /// <summary>exp kernel (derived 2^-143.5, see <see cref="ExpAssemble"/>).</summary>
    private const double EpsExp = 1.1479437019748901E-41; // 2^-136

    /// <summary>expm1 kernel for |x| &lt; 1 (derived 2^-141.9, see <see cref="ExpM1Core"/>).</summary>
    private const double EpsExpM1 = 4.591774807899561E-41; // 2^-134

    /// <summary>sin/cos kernel excluding the reduction term (derived 2^-134, see <see cref="SinCosOf"/>).</summary>
    private const double EpsTrig = 2.938735877055719E-39; // 2^-128

    /// <summary>atan(t), 0 ≤ t ≤ 1, and its π/2 complement (derived 2^-134.9, see <see cref="AtanSmall"/>).</summary>
    private const double EpsAtan = 2.938735877055719E-39; // 2^-128

    private const double P2m149 = 1.4012984643248171E-45, P2m150 = 7.006492321624085E-46, P2m151 = 3.503246160812043E-46;
    private const double P2m153 = 8.758115402030107E-47, P2m155 = 2.1895288505075267E-47, P2m156 = 1.0947644252537633E-47;
    private const double P2m157 = 5.473822126268817E-48, P2m200 = 6.223015277861142E-61;
    private const double P2m101 = 3.944304526105059E-31, P2m143 = 8.968310171678829E-44;

    /// <summary>Safety factor applied to the run-time (absolute-error based) bounds of the log family and pow.</summary>
    private const double DynamicSafety = 128;

    private static Td TdOne => new(1.0, 0, 0);

    // ---------------------------------------------------------------
    // Cody-Waite reduction shared by exp and the trigonometric functions
    // ---------------------------------------------------------------

    /// <summary>
    /// x - nd·C for C = c0 + c1 + c2 + c3 (≈ 212 bits): the products nd·c_i are formed exactly (TwoProd) and
    /// summed by magnitude level (|x|, 2^-53|x|, 2^-106|x|, 2^-159|x|) with TwoSums, so the only roundings are in
    /// the fourth level and in TdRenorm4: error ≤ 2^-208·|x| + 2^-157·|result| (including the 2^-212 truncation of C).
    /// </summary>
    private static Td ReduceCw(in Td x, double nd, double c0, double c1, double c2, double c3)
    {
        double p0 = TwoProd(nd, c0, out double q0);
        double p1 = TwoProd(nd, c1, out double q1);
        double p2 = TwoProd(nd, c2, out double q2);
        double p3 = nd * c3;
        double a = TwoSum(x.H, -p0, out double ea);
        double b = TwoSum(x.M, -q0, out double eb1);
        b = TwoSum(b, -p1, out double eb2);
        b = TwoSum(b, ea, out double eb3);
        double c = TwoSum(x.L, -q1, out double ec1);
        c = TwoSum(c, -p2, out double ec2);
        c = TwoSum(c, eb1, out double ec3);
        c = TwoSum(c, eb2, out double ec4);
        c = TwoSum(c, eb3, out double ec5);
        double d = (ec1 + ec2) + (ec3 + ec4) + (ec5 - q2 - p3);
        return TdRenorm4(a, b, c, d);
    }

    // ---------------------------------------------------------------
    // Exponential
    // ---------------------------------------------------------------

    /// <summary>e^x = W·2^k for |x| ≤ 11400 (N = round(4096x/ln2), x = N·ln2/4096 + r, |r| ≤ 2^-13.5).</summary>
    private static Td ExpCore(in Td x, out int k)
    {
        double nd = Math.Round(x.H * FpExp.K);
        Td r = nd == 0 ? x : ReduceCw(x, nd, FpExp.L0, FpExp.L1, FpExp.L2, FpExp.L3);
        long n = (long)nd;
        k = (int)(n >> 12);
        return ExpAssemble(r, (int)(n & 4095));
    }

    /// <summary>
    /// 2^(idx/4096)·e^r = T + T·p with T = 2^(j1/64)·2^(j2/4096) and p = e^r - 1 (<see cref="ExpPoly"/>).
    /// Error relative to W: p (2^-144.4) + TdMul(T1, T2) (2^-150) + tables (2^-159) + TdMul(T, p) (2^-163)
    /// + TdAdd (2^-153) + the reduction (2^-158 + 2^-171) ≤ 2^-143.5.
    /// </summary>
    private static Td ExpAssemble(in Td r, int idx)
    {
        Td p = ExpPoly(r);
        if (idx == 0)
        {
            return TdAddD(p, 1.0);
        }
        Td t = ExpTable(idx);
        return TdAdd(t, TdMul(t, p));
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static Td ExpTable(int idx)
    {
        int j1 = idx >> 6, j2 = idx & 63;
        return j2 == 0 ? FpExp.T1[j1] : j1 == 0 ? FpExp.T2[j2] : TdMul(FpExp.T1[j1], FpExp.T2[j2]);
    }

    /// <summary>
    /// e^r - 1 for |r| ≤ 2^-13.5 as r + r²/2 + r³·S(r), S = 1/3! + r/4! + ... + r⁶/9! (the r^10/10! tail is
    /// below 2^-157). r² is a triple-double; r³ and S are double-double (S's levels 1/7!..1/9! in double:
    /// they multiply r^7 ≤ 2^-94.7, so a 2^-53 relative error stays below 2^-160). |r³S| ≤ 2^-43.1 with a
    /// relative error ≤ 2^-101.3, so the absolute error is ≤ 2^-144.4.
    /// </summary>
    private static Td ExpPoly(in Td r)
    {
        Td r2 = TdSqr(r);
        double rh = r.H, rm = r.M;
        double r3h = DdMul(r2.H, r2.M, rh, rm, out double r3l);
        double v = FpExp.F7 + rh * (FpExp.F8 + rh * FpExp.F9);
        double sh = DdHorner(FpExp.F6h, FpExp.F6l, rh, rm, v, 0, out double sl);
        sh = DdHorner(FpExp.F5h, FpExp.F5l, rh, rm, sh, sl, out sl);
        sh = DdHorner(FpExp.F4h, FpExp.F4l, rh, rm, sh, sl, out sl);
        sh = DdHorner(FpExp.F3h, FpExp.F3l, rh, rm, sh, sl, out sl);
        double qh = DdMul(r3h, r3l, sh, sl, out double ql);
        Td p = TdAdd(r, new Td(0.5 * r2.H, 0.5 * r2.M, 0.5 * r2.L));
        return TdAdd2(p, qh, ql);
    }

    /// <summary>
    /// e^r - 1 with a small relative error for |r| ≤ 2^-13.5: r + r²/2 + r³/6 + r⁴·S4(r) with r³/6 a
    /// triple-double and S4 = 1/4! + ... + r⁵/9! double-double (|r⁴S4| ≤ 2^-58.6, error ≤ 2^-159.6);
    /// relative error ≤ 2^-151.
    /// </summary>
    private static Td ExpM1Poly(in Td r)
    {
        Td r2 = TdSqr(r);
        Td r3 = TdMul(r2, r);
        double rh = r.H, rm = r.M;
        double r4h = DdMul(r2.H, r2.M, r2.H, r2.M, out double r4l);
        double v = FpExp.F7 + rh * (FpExp.F8 + rh * FpExp.F9);
        double sh = DdHorner(FpExp.F6h, FpExp.F6l, rh, rm, v, 0, out double sl);
        sh = DdHorner(FpExp.F5h, FpExp.F5l, rh, rm, sh, sl, out sl);
        sh = DdHorner(FpExp.F4h, FpExp.F4l, rh, rm, sh, sl, out sl);
        double qh = DdMul(r4h, r4l, sh, sl, out double ql);
        Td p = TdAdd(r, new Td(0.5 * r2.H, 0.5 * r2.M, 0.5 * r2.L));
        p = TdAdd(p, TdMul(r3, FpConst.OneSixth));
        return TdAdd2(p, qh, ql);
    }

    /// <summary>
    /// e^x - 1 for |x| &lt; 1 as (t - 1) + t·p, t = 2^k·T. Absolute error ≤ 2^-148.5 (T's TdMul/table error
    /// 2^-148.9 dominates) when j1 ≠ 0, i.e. |x| ≥ 0.0108 and |e^x - 1| ≥ 2^-6.6; otherwise T is a single table
    /// entry (2^-158.6). With |e^x - 1| ≥ 2^-13.6 for N ≠ 0 the relative error is ≤ 2^-141.9 (2^-149 for N = 0).
    /// </summary>
    private static Td ExpM1Core(in Td x)
    {
        double nd = Math.Round(x.H * FpExp.K);
        if (nd == 0)
        {
            return ExpM1Poly(x);
        }
        Td r = ReduceCw(x, nd, FpExp.L0, FpExp.L1, FpExp.L2, FpExp.L3);
        Td p = ExpM1Poly(r);
        long n = (long)nd;
        int idx = (int)(n & 4095);
        Td t = TdScale(idx == 0 ? TdOne : ExpTable(idx), (int)(n >> 12));
        return TdAdd(TdAddD(t, -1.0), TdMul(t, p));
    }

    /// <summary>|x| in [2^-200, 11400], normal.</summary>
    private static bool ExpArgOk(Float128 x)
    {
        int be = x.BiasedExponent;
        return be >= ExponentBias - 200 && (be < ExponentBias + 13 || (be == ExponentBias + 13 && Abs(x) <= (Float128)11400));
    }

    private static bool KExp(Float128 x, out Td y, out int scale, out double eps)
    {
        eps = EpsExp;
        if (!ExpArgOk(x))
        {
            y = default;
            scale = 0;
            return false;
        }
        y = ExpCore(ToTd(x), out scale);
        return true;
    }

    /// <summary>
    /// 2^x: N = round(4096x), f = x - N/4096 is exact (both are multiples of ulp(x.H) when |x| ≥ 2^-12, and
    /// f = x otherwise), r = f·ln2 (TdMul, absolute error ≤ 2^-162); then as exp (bound EpsExp).
    /// </summary>
    private static bool KExp2(Float128 x, out Td y, out int scale, out double eps)
    {
        eps = EpsExp;
        int be = x.BiasedExponent;
        if (be < ExponentBias - 200 || be > ExponentBias + 14 || !IsFinite(x))
        {
            y = default;
            scale = 0;
            return false;
        }
        Td t = ToTd(x);
        double nd = Math.Round(t.H * 4096);
        Td r = TdMul(nd == 0 ? t : TdAddD(t, -nd * (1.0 / 4096)), FpConst.Ln2);
        long n = (long)nd;
        scale = (int)(n >> 12);
        y = ExpAssemble(r, (int)(n & 4095));
        return true;
    }

    // ---------------------------------------------------------------
    // Logarithms
    // ---------------------------------------------------------------

    /// <summary>
    /// log m for x = 2^e·m (x positive normal), m ∈ [√½, √2), with an absolute error bound err. Two table
    /// steps: z1 = m·c1 - 1 (|z1| ≤ 2^-7.5; exact when c1 = 1), z = z1·c2 + (c2 - 1) (|z| ≤ 2^-14.9; z = z1
    /// when c2 = 1), then log m = -log c1 - log c2 + log1p(z) with log1p(z) = z - z²/2 + z³·P(z) (z² triple-
    /// double, z³ and P double-double, P's levels from z⁴ on in double; truncation after z^10).
    /// Error: log1p ≤ |z|·2^-151.4 + |z|³·2^-101 (the double-double z³P term); the first table step adds
    /// ≤ 2^-152.5 (z1's rounding 2^-153, the table add 2^-154.4), the second ≤ 2^-159.
    /// </summary>
    private static Td LogCore(Float128 x, out int e, out double err)
    {
        UInt128 sig = Significand(x);
        e = x.BiasedExponent - ExponentBias;
        int me = 0;
        if (sig > FpLog.Sqrt2Sig)
        {
            e++;
            me = -1;
        }
        Td m = SigToTd(sig, me, false);
        int v = (int)Math.Round(m.H * 128);
        bool step1 = v != 128;
        Td z = step1 ? TdAddD(TdMulD(m, FpLog.C1[v - 90]), -1.0) : TdAddD(m, -1.0);
        int v2 = (int)Math.Round(z.H * 16384);
        bool step2 = v2 != 0;
        if (step2)
        {
            double c2 = FpLog.C2[v2 + 96];
            z = TdAddD(TdMulD(z, c2), c2 - 1.0);
        }
        Td zsq = TdSqr(z);
        double zh = z.H, zm = z.M;
        double z3h = DdMul(zsq.H, zsq.M, zh, zm, out double z3l);
        double w = FpLog.P7 + zh * (FpLog.P8 + zh * (FpLog.P9 + zh * FpLog.P10));
        double ph = DdHorner(FpLog.P6h, FpLog.P6l, zh, zm, w, 0, out double pl);
        ph = DdHorner(FpLog.P5h, FpLog.P5l, zh, zm, ph, pl, out pl);
        ph = DdHorner(FpLog.P4h, FpLog.P4l, zh, zm, ph, pl, out pl);
        ph = DdHorner(FpLog.P3h, FpLog.P3l, zh, zm, ph, pl, out pl);
        double qh = DdMul(z3h, z3l, ph, pl, out double ql);
        Td lp = TdAdd2(TdAdd(z, new Td(-0.5 * zsq.H, -0.5 * zsq.M, -0.5 * zsq.L)), qh, ql);
        double az = Math.Abs(zh);
        err = az * P2m150 + az * az * az * P2m101;
        if (step2)
        {
            lp = TdAdd(FpLog.NL2[v2 + 96], lp);
            err += P2m157;
        }
        if (step1)
        {
            lp = TdAdd(FpLog.NL1[v - 90], lp);
            err += P2m151;
        }
        return lp;
    }

    /// <summary>e·ln2 as a triple-double; error ≤ |e|·2^-156 (ln2 in three parts: 2^-160 truncation).</summary>
    private static Td Ln2Multiple(int e)
    {
        double p0 = TwoProd(e, FpConst.Ln2P0, out double q0);
        double p1 = TwoProd(e, FpConst.Ln2P1, out double q1);
        double s1 = TwoSum(q0, p1, out double e1);
        return TdRenorm(p0, s1, Math.FusedMultiplyAdd(e, FpConst.Ln2P2, q1 + e1));
    }

    /// <summary>ln x = e·ln2 + log m with its absolute error bound (x positive normal).</summary>
    private static Td LnCore(Float128 x, out double err)
    {
        Td lm = LogCore(x, out int e, out err);
        if (e == 0)
        {
            return lm;
        }
        Td l = TdAdd(Ln2Multiple(e), lm);
        err += Math.Abs(e) * (P2m156 + P2m153) + P2m153;
        return l;
    }

    private static bool LogArgOk(Float128 x) => !x.SignBit && x.BiasedExponent is > 0 and < MaxBiasedExponent;

    private static bool KLog(Float128 x, out Td y, out int scale, out double eps)
    {
        scale = 0;
        if (!LogArgOk(x))
        {
            y = default;
            eps = 0;
            return false;
        }
        y = LnCore(x, out double err);
        eps = RelBound(err, y);
        return true;
    }

    /// <summary>log2 x = e + log m·(1/ln2): TdMul adds 2^-149·|t| (with the constant's 2^-160), the add 2^-155(|t| + |e|).</summary>
    private static bool KLog2(Float128 x, out Td y, out int scale, out double eps)
    {
        scale = 0;
        if (!LogArgOk(x))
        {
            y = default;
            eps = 0;
            return false;
        }
        Td lm = LogCore(x, out int e, out double err);
        Td t = TdMul(lm, FpConst.InvLn2);
        double at = Math.Abs(t.H);
        err = err * 1.4426950408889634 * (1 + 1e-10) + at * P2m149;
        if (e != 0)
        {
            t = TdAddD(t, e);
            err += (at + Math.Abs(e)) * P2m155;
        }
        y = t;
        eps = RelBound(err, y);
        return true;
    }

    /// <summary>log10 x = ln x·(1/ln10) (TdMul and constant: 2^-149·|result|).</summary>
    private static bool KLog10(Float128 x, out Td y, out int scale, out double eps)
    {
        scale = 0;
        if (!LogArgOk(x))
        {
            y = default;
            eps = 0;
            return false;
        }
        Td l = LnCore(x, out double err);
        y = TdMul(l, FpConst.InvLn10);
        eps = RelBound(err * 0.43429448190325182 * (1 + 1e-10) + Math.Abs(y.H) * P2m149, y);
        return true;
    }

    /// <summary>A relative bound from an absolute one (|y| ≥ |y.H|(1 - 2^-52)); +∞ for y = 0.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double RelBound(double err, in Td y) => err / Math.Abs(y.H) * DynamicSafety;

    // ---------------------------------------------------------------
    // Power
    // ---------------------------------------------------------------

    /// <summary>
    /// x^y = ±e^z, z = y·ln|x| (TdMul): z's absolute error is |y|·err(ln|x|) + 2^-149·|z|, which is the
    /// relative error it causes in e^z (up to a factor 1 + 2^-100); the exp kernel adds EpsExp.
    /// </summary>
    private static bool KPow(Float128 x, Float128 y, out Td res, out int scale, out double eps)
    {
        res = default;
        scale = 0;
        eps = 0;
        int bx = x.BiasedExponent, by = y.BiasedExponent;
        if (bx is 0 or MaxBiasedExponent || by < ExponentBias - 800 || by > ExponentBias + 900 || Abs(x) == One)
        {
            return false;
        }
        bool neg = false;
        if (x.SignBit)
        {
            if (!IsInteger(y))
            {
                return false;
            }
            neg = IsOddInteger(y);
        }
        Td l = LnCore(Abs(x), out double err);
        Td yt = ToTd(y);
        Td z = TdMul(yt, l);
        if (!(Math.Abs(z.H) < 11400))
        {
            return false;
        }
        res = ExpCore(z, out scale);
        if (neg)
        {
            res = -res;
        }
        eps = EpsExp + (Math.Abs(yt.H) * err + Math.Abs(z.H) * P2m149) * DynamicSafety;
        return true;
    }

    // ---------------------------------------------------------------
    // Trigonometric functions
    // ---------------------------------------------------------------

    /// <summary>
    /// sin r and cos r - 1 for |r| ≤ π/512·(1 + 2^-40) ≤ 2^-7.35. sin r = r - r³/6 + r⁵·S5(r²) with r³/6 a
    /// triple-double and S5 double-double (levels from r^11 in double, truncation after r^15): |r⁵S5| ≤ 2^-43.6
    /// with relative error 2^-101.3, so sr errs by ≤ |r|·2^-137.5. cos r - 1 = -r²/2 + r⁴·C4(r²) (truncation
    /// after r^14): |r⁴C4| ≤ 2^-34, so cm1 errs by ≤ 2^-135.3 absolutely.
    /// </summary>
    private static void SinCosPoly(in Td r, bool needSin, bool needCos, out Td sr, out Td cm1)
    {
        Td r2 = TdSqr(r);
        double zh = r2.H, zl = r2.M;
        sr = default;
        cm1 = default;
        if (needSin)
        {
            Td r3 = TdMul(r2, r);
            double ds = FpTrig.S11 + zh * (FpTrig.S13 + zh * FpTrig.S15);
            double sh = DdHorner(FpTrig.S9h, FpTrig.S9l, zh, zl, ds, 0, out double sl);
            sh = DdHorner(FpTrig.S7h, FpTrig.S7l, zh, zl, sh, sl, out sl);
            sh = DdHorner(FpTrig.S5h, FpTrig.S5l, zh, zl, sh, sl, out sl);
            double r5h = DdMul(r3.H, r3.M, zh, zl, out double r5l);
            double qh = DdMul(r5h, r5l, sh, sl, out double ql);
            sr = TdAdd2(TdSub(r, TdMul(r3, FpConst.OneSixth)), qh, ql);
        }
        if (needCos)
        {
            double dc = FpTrig.C10 + zh * (FpTrig.C12 + zh * FpTrig.C14);
            double ch = DdHorner(FpTrig.C8h, FpTrig.C8l, zh, zl, dc, 0, out double cl);
            ch = DdHorner(FpTrig.C6h, FpTrig.C6l, zh, zl, ch, cl, out cl);
            ch = DdHorner(FpTrig.C4h, FpTrig.C4l, zh, zl, ch, cl, out cl);
            double r4h = DdMul(zh, zl, zh, zl, out double r4l);
            double qh = DdMul(r4h, r4l, ch, cl, out double ql);
            cm1 = TdAdd2(new Td(-0.5 * r2.H, -0.5 * r2.M, -0.5 * r2.L), qh, ql);
        }
    }

    /// <summary>
    /// x = N·π/256 + r (|x| &lt; 2^30, so the reduction errs by ≤ 2^-178 + 2^-157·|r|), a = (N mod 128)·π/256.
    /// Returns sin(a + r) (wantCos false) and/or cos(a + r) and the quadrant q = (N div 128) mod 4.
    /// For a ≠ 0: sin(a + r) = S + (S·cm1 + C·sr), cos(a + r) = C + (C·cm1 - S·sr); since
    /// |sin(a + r)| ≥ max(S/2, sin(π/512)) (and likewise for cos with a ≤ 127π/256), the errors of cm1 and sr
    /// give ≤ 2·2^-135.3 + 2^-137.5·2^-7.35/2^-7.35 relative, plus 2^-150 for the products and adds: ≤ 2^-134.
    /// For a = 0 the result is sr (≤ 2^-137.5 relative) or 1 + cm1 (≤ 2^-135.2).
    /// </summary>
    private static void SinCosOf(in Td x, bool wantSin, bool wantCos, out Td sinA, out Td cosA, out int q)
    {
        double nd = Math.Round(x.H * FpTrig.K);
        Td r = nd == 0 ? x : ReduceCw(x, nd, FpTrig.P0, FpTrig.P1, FpTrig.P2, FpTrig.P3);
        int j = (int)((long)nd & 511);
        q = j >> 7;
        int jj = j & 127;
        sinA = default;
        cosA = default;
        if (jj == 0)
        {
            SinCosPoly(r, wantSin, wantCos, out sinA, out Td cm1);
            if (wantCos)
            {
                cosA = TdAddD(cm1, 1.0);
            }
            return;
        }
        SinCosPoly(r, true, true, out Td sr, out Td cm);
        Td s = FpTrig.S[jj], c = FpTrig.C[jj];
        if (wantSin)
        {
            sinA = TdAdd(s, TdAdd(TdMul(s, cm), TdMul(c, sr)));
        }
        if (wantCos)
        {
            cosA = TdAdd(c, TdSub(TdMul(c, cm), TdMul(s, sr)));
        }
    }

    private static bool TrigArgOk(Float128 x)
    {
        int be = x.BiasedExponent;
        return be >= ExponentBias - 200 && be < ExponentBias + 30;
    }

    /// <summary>The reduction's absolute error (2^-178 + 2^-157|r|) relative to a result y, with margin.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double TrigReductionBound(in Td x, in Td y) => P2m200 * Math.Abs(x.H) / Math.Abs(y.H) + P2m153;

    private static bool KSin(Float128 x, out Td y, out int scale, out double eps)
    {
        scale = 0;
        if (!TrigArgOk(x))
        {
            y = default;
            eps = 0;
            return false;
        }
        Td tx = ToTd(x);
        int q;
        Td s, c;
        // sin x = sin(qπ/2 + a + r): q = 0 → sinA, 1 → cosA, 2 → -sinA, 3 → -cosA. Only the needed one is formed,
        // so the quadrant must be known first: compute it from N directly.
        double nd = Math.Round(tx.H * FpTrig.K);
        bool odd = (((long)nd >> 7) & 1) != 0;
        SinCosOf(tx, !odd, odd, out s, out c, out q);
        y = odd ? c : s;
        if (q >= 2)
        {
            y = -y;
        }
        eps = EpsTrig + TrigReductionBound(tx, y);
        return true;
    }

    private static bool KCos(Float128 x, out Td y, out int scale, out double eps)
    {
        scale = 0;
        if (!TrigArgOk(x))
        {
            y = default;
            eps = 0;
            return false;
        }
        Td tx = ToTd(x);
        // cos x = cos(qπ/2 + a + r): q = 0 → cosA, 1 → -sinA, 2 → -cosA, 3 → sinA.
        double nd = Math.Round(tx.H * FpTrig.K);
        bool odd = (((long)nd >> 7) & 1) != 0;
        SinCosOf(tx, odd, !odd, out Td s, out Td c, out int q);
        y = odd ? s : c;
        if (q is 1 or 2)
        {
            y = -y;
        }
        eps = EpsTrig + TrigReductionBound(tx, y);
        return true;
    }

    /// <summary>tan x = sinA/cosA (q even) or -cosA/sinA (q odd): two trig bounds plus TdDiv's 2^-150.</summary>
    private static bool KTan(Float128 x, out Td y, out int scale, out double eps)
    {
        scale = 0;
        if (!TrigArgOk(x))
        {
            y = default;
            eps = 0;
            return false;
        }
        Td tx = ToTd(x);
        SinCosOf(tx, true, true, out Td s, out Td c, out int q);
        y = (q & 1) == 0 ? TdDiv(s, c) : -TdDiv(c, s);
        eps = 2 * EpsTrig + TrigReductionBound(tx, s) + TrigReductionBound(tx, c) + P2m149;
        return true;
    }

    // ---------------------------------------------------------------
    // Inverse trigonometric functions
    // ---------------------------------------------------------------

    /// <summary>
    /// atan(u) for |u| ≤ 2^-7.9: u - u³/3 + u⁵·A5(u²) with u³/3 a triple-double and A5 double-double (levels
    /// from u^13 in double, truncation after u^19): |u⁵A5| ≤ 2^-42.3 with relative error 2^-101.3, i.e.
    /// ≤ 2^-143.6 absolutely and ≤ 2^-135.6 relative to u.
    /// </summary>
    private static Td AtanPoly(in Td u)
    {
        Td u2 = TdSqr(u);
        Td u3 = TdMul(u2, u);
        double zh = u2.H, zl = u2.M;
        double d = FpAtan.A13 + zh * (FpAtan.A15 + zh * (FpAtan.A17 + zh * FpAtan.A19));
        double ah = DdHorner(FpAtan.A11h, FpAtan.A11l, zh, zl, d, 0, out double al);
        ah = DdHorner(FpAtan.A9h, FpAtan.A9l, zh, zl, ah, al, out al);
        ah = DdHorner(FpAtan.A7h, FpAtan.A7l, zh, zl, ah, al, out al);
        ah = DdHorner(FpAtan.A5h, FpAtan.A5l, zh, zl, ah, al, out al);
        double u5h = DdMul(u3.H, u3.M, zh, zl, out double u5l);
        double qh = DdMul(u5h, u5l, ah, al, out double ql);
        return TdAdd2(TdSub(u, TdMul(u3, FpConst.OneThird)), qh, ql);
    }

    /// <summary>
    /// atan(t) for 0 ≤ t ≤ 1 + 2^-50: c = j/128 (j = round(128t)), u = (t - c)/(1 + tc) (|u| ≤ 2^-8; TdDiv
    /// 2^-150, the numerator and denominator ≤ 2^-154 absolutely), atan t = atan(c) + atan(u). The result is
    /// ≥ max(2^-8, |atan u|) with at most a factor 2 cancellation, so the relative error is
    /// ≤ 2·2^-135.6·... ≤ 2^-134.9 (j = 0: atan u alone, ≤ 2^-135.6).
    /// </summary>
    private static Td AtanSmall(in Td t)
    {
        int j = (int)Math.Round(t.H * 128);
        if (j == 0)
        {
            return AtanPoly(t);
        }
        double c = j * (1.0 / 128);
        Td u = TdDiv(TdAddD(t, -c), TdAddD(TdMulD(t, c), 1.0));
        return TdAdd(FpAtan.A[j], AtanPoly(u));
    }

    /// <summary>atan(n/d) ∈ [0, π/2] for n, d &gt; 0 (π/2 - atan(d/n) when n &gt; d: no cancellation, result ≥ π/4).</summary>
    private static Td AtanRatio(in Td n, in Td d) =>
        n.H <= d.H ? AtanSmall(TdDiv(n, d)) : TdSub(FpConst.PiHalf, AtanSmall(TdDiv(d, n)));

    private static bool KAtan(Float128 x, out Td y, out int scale, out double eps)
    {
        scale = 0;
        eps = EpsAtan;
        int be = x.BiasedExponent;
        if (be < ExponentBias - 200 || be > ExponentBias + 900)
        {
            y = default;
            return false;
        }
        Td a = ToTd(Abs(x));
        y = a.H <= 1 ? AtanSmall(a) : TdSub(FpConst.PiHalf, AtanSmall(TdDiv(TdOne, a)));
        if (x.SignBit)
        {
            y = -y;
        }
        return true;
    }

    /// <summary>
    /// atan2(y, x) = ±atan(|y|/|x|) or ±(π - atan(|y|/|x|)) for x &lt; 0 (≥ π/2, no cancellation), both scaled
    /// by the same power of two; the ratio's TdDiv adds 2^-150.
    /// </summary>
    private static bool KAtan2(Float128 yv, Float128 xv, out Td res, out int scale, out double eps)
    {
        res = default;
        scale = 0;
        eps = EpsAtan + P2m149;
        int by = yv.BiasedExponent, bx = xv.BiasedExponent;
        if (by is 0 or MaxBiasedExponent || bx is 0 or MaxBiasedExponent || Math.Abs(by - bx) > 600)
        {
            return false;
        }
        int shift = ExponentBias - Math.Max(by, bx);
        Td n = ToTd(Abs(yv), shift), d = ToTd(Abs(xv), shift);
        Td a = AtanRatio(n, d);
        if (xv.SignBit)
        {
            a = TdSub(FpConst.Pi, a);
        }
        res = yv.SignBit ? -a : a;
        return true;
    }

    /// <summary>
    /// √(1 - a²) for 0 ≤ a &lt; 1 as √((1 - a)(1 + a)) (1 ± a: TdAddD, ≤ 2^-154 relative; product and root
    /// 2^-150 each): relative error ≤ 2^-149.
    /// </summary>
    private static Td SqrtOneMinusSquare(in Td a) => TdSqrt(TdMul(TdAddD(-a, 1.0), TdAddD(a, 1.0)));

    /// <summary>asin x = atan(|x| / √(1 - x²)) (the ratio errs by ≤ 2^-148, atan's condition number is ≤ 1).</summary>
    private static bool KAsin(Float128 x, out Td y, out int scale, out double eps)
    {
        scale = 0;
        eps = EpsAtan + P2m149 * 4;
        int be = x.BiasedExponent;
        if (be < ExponentBias - 200 || be >= ExponentBias)
        {
            y = default;
            return false;
        }
        Td a = ToTd(Abs(x));
        y = AtanRatio(a, SqrtOneMinusSquare(a));
        if (x.SignBit)
        {
            y = -y;
        }
        return true;
    }

    /// <summary>acos x = atan(√(1 - x²)/|x|), or π minus that for x &lt; 0 (≥ π/2).</summary>
    private static bool KAcos(Float128 x, out Td y, out int scale, out double eps)
    {
        scale = 0;
        eps = EpsAtan + P2m149 * 4;
        int be = x.BiasedExponent;
        if (be < ExponentBias - 200 || be >= ExponentBias)
        {
            y = default;
            return false;
        }
        Td a = ToTd(Abs(x));
        y = AtanRatio(SqrtOneMinusSquare(a), a);
        if (x.SignBit)
        {
            y = TdSub(FpConst.Pi, y);
        }
        return true;
    }

    // ---------------------------------------------------------------
    // Hyperbolic functions
    // ---------------------------------------------------------------

    /// <summary>
    /// sinh: |x| &lt; 0.5: (em + em/(em + 1))/2 with em = expm1|x| (EpsExpM1 + 2^-149); |x| &lt; 50:
    /// (E - 1/E)/2 with cancellation factor ≤ coth(0.5) &lt; 2.2 (2.2·(2^-143.5 + 2^-143.5 + 2^-150)
    /// ≤ 2^-141.4, within EpsExp·2.2); |x| ≥ 50: E/2, neglecting e^-2|x| ≤ 2^-144.
    /// </summary>
    private static bool KSinh(Float128 x, out Td y, out int scale, out double eps)
    {
        y = default;
        scale = 0;
        eps = 0;
        if (!ExpArgOk(x))
        {
            return false;
        }
        Td a = ToTd(Abs(x));
        if (a.H < 0.5)
        {
            Td em = ExpM1Core(a);
            y = TdScale(TdAdd(em, TdDiv(em, TdAddD(em, 1.0))), -1);
            eps = EpsExpM1 + P2m149;
        }
        else if (a.H < 50)
        {
            Td w = ExpCore(a, out int k);
            Td e = TdScale(w, k);
            y = TdScale(TdSub(e, TdDiv(TdOne, e)), -1);
            eps = 2.2 * EpsExp;
        }
        else
        {
            y = ExpCore(a, out int k);
            scale = k - 1;
            eps = EpsExp + P2m143;
        }
        if (x.SignBit)
        {
            y = -y;
        }
        return true;
    }

    /// <summary>cosh: |x| &lt; 50: (E + 1/E)/2 (no cancellation: EpsExp + 2^-149); |x| ≥ 50: E/2 (+2^-143).</summary>
    private static bool KCosh(Float128 x, out Td y, out int scale, out double eps)
    {
        y = default;
        scale = 0;
        eps = 0;
        if (!ExpArgOk(x))
        {
            return false;
        }
        Td a = ToTd(Abs(x));
        Td w = ExpCore(a, out int k);
        if (a.H < 50)
        {
            Td e = TdScale(w, k);
            y = TdScale(TdAdd(e, TdDiv(TdOne, e)), -1);
            eps = EpsExp + P2m149;
        }
        else
        {
            y = w;
            scale = k - 1;
            eps = EpsExp + P2m143;
        }
        return true;
    }

    /// <summary>
    /// tanh: |x| &lt; 0.5: em/(em + 2) with em = expm1(2|x|) (2·EpsExpM1 + 2^-149); |x| &lt; 41:
    /// (E - 1)/(E + 1) with E = e^(2|x|) ≥ e (E/(E - 1) ≤ 1.6: 2.6·EpsExp + 2^-149); larger |x| are answered
    /// by the caller (tanh x = 1 - δ, δ &lt; 2^-117).
    /// </summary>
    private static bool KTanh(Float128 x, out Td y, out int scale, out double eps)
    {
        y = default;
        scale = 0;
        eps = 0;
        int be = x.BiasedExponent;
        if (be < ExponentBias - 200 || be > ExponentBias + 5 || !(Abs(x) < (Float128)41))
        {
            return false;
        }
        Td a2 = ToTd(Abs(x), 1);
        if (a2.H < 1)
        {
            Td em = ExpM1Core(a2);
            y = TdDiv(em, TdAddD(em, 2.0));
            eps = 2 * EpsExpM1 + P2m149;
        }
        else
        {
            Td w = ExpCore(a2, out int k);
            Td e = TdScale(w, k);
            y = TdDiv(TdAddD(e, -1.0), TdAddD(e, 1.0));
            eps = 2.6 * EpsExp + P2m149;
        }
        if (x.SignBit)
        {
            y = -y;
        }
        return true;
    }

    // ---------------------------------------------------------------
    // Entry points
    // ---------------------------------------------------------------

    /// <summary>Functions with a fast path (test hooks <see cref="FastPathKernel"/> and <see cref="FastPathTry"/> take their index).</summary>
    internal enum FastFunction
    {
        Exp, Exp2, Log, Log2, Log10, Pow, Sin, Cos, Tan, Asin, Acos, Atan, Atan2, Sinh, Cosh, Tanh,
    }

    private static bool Kernel(FastFunction f, Float128 x, Float128 y, out Td r, out int scale, out double eps) => f switch
    {
        FastFunction.Exp => KExp(x, out r, out scale, out eps),
        FastFunction.Exp2 => KExp2(x, out r, out scale, out eps),
        FastFunction.Log => KLog(x, out r, out scale, out eps),
        FastFunction.Log2 => KLog2(x, out r, out scale, out eps),
        FastFunction.Log10 => KLog10(x, out r, out scale, out eps),
        FastFunction.Pow => KPow(x, y, out r, out scale, out eps),
        FastFunction.Sin => KSin(x, out r, out scale, out eps),
        FastFunction.Cos => KCos(x, out r, out scale, out eps),
        FastFunction.Tan => KTan(x, out r, out scale, out eps),
        FastFunction.Asin => KAsin(x, out r, out scale, out eps),
        FastFunction.Acos => KAcos(x, out r, out scale, out eps),
        FastFunction.Atan => KAtan(x, out r, out scale, out eps),
        FastFunction.Atan2 => KAtan2(x, y, out r, out scale, out eps),
        FastFunction.Sinh => KSinh(x, out r, out scale, out eps),
        FastFunction.Cosh => KCosh(x, out r, out scale, out eps),
        FastFunction.Tanh => KTanh(x, out r, out scale, out eps),
        _ => throw new ArgumentOutOfRangeException(nameof(f)),
    };

    /// <summary>Test hook: the kernel's y = h + m + l, scale and relative bound, or false outside the fast path.</summary>
    internal static bool FastPathKernel(int function, Float128 x, Float128 y, out double h, out double m, out double l, out int scale, out double eps)
    {
        bool ok = Kernel((FastFunction)function, x, y, out Td r, out scale, out eps);
        (h, m, l) = (r.H, r.M, r.L);
        return ok;
    }

    /// <summary>Test hook: the fast path alone (false when it would fall back).</summary>
    internal static bool FastPathTry(int function, Float128 x, Float128 y, RealRounding rounding, out Float128 result)
    {
        result = default;
        return (uint)rounding <= (uint)RealRounding.AwayFromZero
            && Kernel((FastFunction)function, x, y, out Td r, out int scale, out double eps)
            && TryRound(r, scale, eps, rounding, out result);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static bool ValidMode(RealRounding rounding) => (uint)rounding <= (uint)RealRounding.AwayFromZero;

    private static bool TryExpFast(Float128 x, RealRounding rounding, out Float128 result)
    {
        result = default;
        return ValidMode(rounding) && KExp(x, out Td y, out int s, out double e) && TryRound(y, s, e, rounding, out result);
    }

    private static bool TryExp2Fast(Float128 x, RealRounding rounding, out Float128 result)
    {
        result = default;
        return ValidMode(rounding) && KExp2(x, out Td y, out int s, out double e) && TryRound(y, s, e, rounding, out result);
    }

    private static bool TryLogFast(Float128 x, RealRounding rounding, out Float128 result)
    {
        result = default;
        return ValidMode(rounding) && KLog(x, out Td y, out int s, out double e) && TryRound(y, s, e, rounding, out result);
    }

    private static bool TryLog2Fast(Float128 x, RealRounding rounding, out Float128 result)
    {
        result = default;
        return ValidMode(rounding) && KLog2(x, out Td y, out int s, out double e) && TryRound(y, s, e, rounding, out result);
    }

    private static bool TryLog10Fast(Float128 x, RealRounding rounding, out Float128 result)
    {
        result = default;
        return ValidMode(rounding) && KLog10(x, out Td y, out int s, out double e) && TryRound(y, s, e, rounding, out result);
    }

    private static bool TryPowFast(Float128 x, Float128 y, RealRounding rounding, out Float128 result)
    {
        result = default;
        return ValidMode(rounding) && KPow(x, y, out Td r, out int s, out double e) && TryRound(r, s, e, rounding, out result);
    }

    private static bool TrySinFast(Float128 x, RealRounding rounding, out Float128 result)
    {
        result = default;
        return ValidMode(rounding) && KSin(x, out Td y, out int s, out double e) && TryRound(y, s, e, rounding, out result);
    }

    private static bool TryCosFast(Float128 x, RealRounding rounding, out Float128 result)
    {
        result = default;
        return ValidMode(rounding) && KCos(x, out Td y, out int s, out double e) && TryRound(y, s, e, rounding, out result);
    }

    private static bool TryTanFast(Float128 x, RealRounding rounding, out Float128 result)
    {
        result = default;
        return ValidMode(rounding) && KTan(x, out Td y, out int s, out double e) && TryRound(y, s, e, rounding, out result);
    }

    private static bool TryAsinFast(Float128 x, RealRounding rounding, out Float128 result)
    {
        result = default;
        return ValidMode(rounding) && KAsin(x, out Td y, out int s, out double e) && TryRound(y, s, e, rounding, out result);
    }

    private static bool TryAcosFast(Float128 x, RealRounding rounding, out Float128 result)
    {
        result = default;
        return ValidMode(rounding) && KAcos(x, out Td y, out int s, out double e) && TryRound(y, s, e, rounding, out result);
    }

    private static bool TryAtanFast(Float128 x, RealRounding rounding, out Float128 result)
    {
        result = default;
        return ValidMode(rounding) && KAtan(x, out Td y, out int s, out double e) && TryRound(y, s, e, rounding, out result);
    }

    private static bool TryAtan2Fast(Float128 y, Float128 x, RealRounding rounding, out Float128 result)
    {
        result = default;
        return ValidMode(rounding) && KAtan2(y, x, out Td r, out int s, out double e) && TryRound(r, s, e, rounding, out result);
    }

    private static bool TrySinhFast(Float128 x, RealRounding rounding, out Float128 result)
    {
        result = default;
        return ValidMode(rounding) && KSinh(x, out Td y, out int s, out double e) && TryRound(y, s, e, rounding, out result);
    }

    private static bool TryCoshFast(Float128 x, RealRounding rounding, out Float128 result)
    {
        result = default;
        return ValidMode(rounding) && KCosh(x, out Td y, out int s, out double e) && TryRound(y, s, e, rounding, out result);
    }

    private static bool TryTanhFast(Float128 x, RealRounding rounding, out Float128 result)
    {
        result = default;
        return ValidMode(rounding) && KTanh(x, out Td y, out int s, out double e) && TryRound(y, s, e, rounding, out result);
    }
}
