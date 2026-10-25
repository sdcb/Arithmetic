using Xunit;

namespace Sdcb.Arithmetic.Tests;

/// <summary>
/// Special functions against published digits, exact special values and identities.
/// (Bit-exact agreement with native MPFR is checked by `tests/Sdcb.Arithmetic.Benchmarks -- exact`.)
/// </summary>
public class SpecialFunctionTests
{
    private const int P = 256;

    // Published constants (OEIS A001620, A006752, A002193-style expansions).
    private const string EulerDigits = "0.57721566490153286060651209008240243104215933593992359880576723488486772677766467093694706329174674951463144724980708248096050401448654283622417399764492353625350033374293733773767394279259525824709491600873520394816567085323315177661152862119950150798479374508570574002992135478614669402960432542151905877553526733139925401296742051375413954911";
    private const string CatalanDigits = "0.91596559417721901505460351493238411077414937428167213426649811962176301977625476947935651292611510624857442261919619957903589880332585905943159473748115840699533202877331946051903872747816408786590902";
    private const string Zeta3Digits = "1.20205690315959428539973816151144999076498629234049888179227155534183820578631309018645587360933525814619915779526071941849199599867328321377639683720790016145394178294936006671919157552224249424396157";
    private const string Ai0Digits = "0.35502805388781723926006318600418317639797917419917724058332651030081004245";
    private const string J0of1Digits = "0.76519768655796655145";

    private static void AssertAgrees(Real actual, string expected, int digits)
    {
        Real e = Real.Parse(expected, P + 64);
        Real diff = Real.Subtract(actual, e, 64).Abs();
        if (diff.IsZero)
        {
            return;
        }
        Real bound = Real.Multiply(e.Abs(), Real.Exp10(Real.FromDouble(-digits, 64), 64), 64);
        Assert.True(diff.CompareTo(bound) <= 0, $"{actual} vs {expected}");
    }

    private static Real R(double v) => Real.FromDouble(v, P);

    [Fact]
    public void Constants()
    {
        AssertAgrees(Real.EulerGamma(P), EulerDigits, 70);
        AssertAgrees(Real.CatalanConstant(P), CatalanDigits, 70);
        AssertAgrees(Real.Zeta(R(3), P), Zeta3Digits, 70);
        AssertAgrees(Real.AiryAi(Real.Zero, P), Ai0Digits, 60);
        AssertAgrees(Real.BesselJ(0, R(1), P), J0of1Digits, 18);
    }

    [Fact]
    public void Gamma_Identities()
    {
        // Γ(n) = (n-1)! exactly.
        Assert.Equal(Real.FromBigInt(BigInt.Factorial(19), P), Real.Gamma(R(20), P));
        // Γ(1/2)² = π, Γ(x+1) = x Γ(x), lnΓ(1) = lnΓ(2) = +0.
        AssertAgrees(Real.Square(Real.Gamma(R(0.5), P), P), Real.Pi(P).ToString(), 70);
        Real x = R(3.25);
        AssertAgrees(Real.Gamma(R(4.25), P), Real.Multiply(x, Real.Gamma(x, P), P).ToString(), 70);
        Assert.True(Real.LogGamma(Real.One, P).IsZero && Real.LogGamma(R(2), P).IsZero);
        // Reflection: Γ(x) Γ(1-x) = π / sin(πx).
        Real y = R(-2.75);
        Real lhs = Real.Multiply(Real.Gamma(y, P), Real.Gamma(Real.Subtract(Real.One, y, P), P), P);
        AssertAgrees(lhs, Real.Divide(Real.Pi(P), Real.SinPi(y, P), P).ToString(), 65);
        // Signs and poles.
        Assert.True(Real.Gamma(R(-0.5), P).IsNegative());
        Assert.True(Real.Gamma(R(-3), P).IsNaN);
        Assert.True(Real.LogGamma(R(-0.5), P).IsNaN); // Γ(-0.5) < 0
        Real lg = Real.LogGamma(R(-0.5), out int sign, P);
        Assert.Equal(-1, sign);
        AssertAgrees(lg, Real.Log(Real.Gamma(R(-0.5), P).Abs(), P).ToString(), 70);
        // ψ(1) = -γ, ψ(x+1) = ψ(x) + 1/x.
        AssertAgrees(Real.Digamma(Real.One, P), (-Real.EulerGamma(P)).ToString(), 70);
        AssertAgrees(Real.Digamma(R(4.5), P), Real.Add(Real.Digamma(R(3.5), P), Real.Divide(Real.One, R(3.5), P), P).ToString(), 70);
        // B(2, 3) = 1/12 and symmetry.
        AssertAgrees(Real.Beta(R(2), R(3), P), Real.Divide(Real.One, R(12), P).ToString(), 70);
        Assert.Equal(Real.Beta(R(1.5), R(2.25), P), Real.Beta(R(2.25), R(1.5), P));
    }

    [Fact]
    public void ErfZetaEintLi2_Identities()
    {
        Real x = R(0.8125);
        AssertAgrees(Real.Add(Real.Erf(x, P), Real.Erfc(x, P), P), "1", 70);
        Assert.Equal(-Real.Erf(x, P), Real.Erf(-x, P));
        Assert.True(Real.Erfc(R(30), 64).IsPositive() && Real.Erfc(R(30), 64).BinaryExponent < -1290);
        // ζ(2) = π²/6, ζ(4) = π⁴/90, ζ(0) = -1/2, ζ(-1) = -1/12, ζ(-2) = 0.
        Real pi = Real.Pi(P);
        AssertAgrees(Real.Zeta(R(2), P), Real.Divide(Real.Square(pi, P), R(6), P).ToString(), 70);
        AssertAgrees(Real.Zeta(R(4), P), Real.Divide(Real.Pow(pi, 4, P), R(90), P).ToString(), 70);
        Assert.Equal(R(-0.5), Real.Zeta(Real.Zero, P));
        AssertAgrees(Real.Zeta(R(-1), P), Real.Divide(Real.MinusOne, R(12), P).ToString(), 70);
        Assert.True(Real.Zeta(R(-2), P).IsZero);
        // Dilogarithm(1) = π²/6, Dilogarithm(-1) = -π²/12, Dilogarithm(1/2) = π²/12 - ln²2/2.
        Real pi2 = Real.Square(pi, P);
        AssertAgrees(Real.Dilogarithm(Real.One, P), Real.Divide(pi2, R(6), P).ToString(), 70);
        AssertAgrees(Real.Dilogarithm(Real.MinusOne, P), Real.Divide(-pi2, R(12), P).ToString(), 70);
        Real half = Real.Subtract(Real.Divide(pi2, R(12), P), Real.Square(Real.Ln2(P), P).ShiftRight(1), P);
        AssertAgrees(Real.Dilogarithm(R(0.5), P), half.ToString(), 70);
        // Ei'(x) = e^x / x: central difference agrees.
        Real h = Real.FromDouble(1.0 / (1 << 20), P);
        Real d = Real.Divide(Real.Subtract(Real.ExponentialIntegral(Real.Add(x, h, P), P), Real.ExponentialIntegral(Real.Subtract(x, h, P), P), P), h.ShiftLeft(1), P);
        AssertAgrees(d, Real.Divide(Real.Exp(x, P), x, P).ToString(), 10);
    }

    [Fact]
    public void Bessel_Wronskian()
    {
        // J_{n+1}(x) Y_n(x) - J_n(x) Y_{n+1}(x) = 2 / (π x)
        foreach (double xv in (double[])[0.75, 3.5, 17.25])
        {
            Real x = R(xv);
            for (int n = 0; n < 3; n++)
            {
                Real w = Real.Subtract(Real.Multiply(Real.BesselJ(n + 1, x, P), Real.BesselY(n, x, P), P), Real.Multiply(Real.BesselJ(n, x, P), Real.BesselY(n + 1, x, P), P), P);
                AssertAgrees(w, Real.Divide(R(2), Real.Multiply(Real.Pi(P), x, P), P).ToString(), 60);
            }
        }
        Assert.Equal(-Real.BesselJ(1, R(2.5), P), Real.BesselJ(1, R(-2.5), P));
    }

    [Fact]
    public void PiFunctionsAndM1Variants_ExactCases()
    {
        Assert.Equal(Real.One.WithPrecision(53), Real.SinPi(R(0.5), 53));
        Assert.Equal(Real.MinusOne.WithPrecision(53), Real.CosPi(R(3), 53));
        Assert.Equal(Real.One.WithPrecision(53), Real.TanPi(R(0.25), 53));
        Assert.True(Real.SinPi(R(-4), 53) is { IsZero: true, SignBit: true });
        Assert.True(Real.CosPi(R(2.5), 53) is { IsZero: true, SignBit: false });
        Assert.True(Real.TanPi(R(0.5), 53).IsPositiveInfinity);
        Assert.Equal(R(7), Real.Exp2M1(R(3), P));
        Assert.Equal(R(999), Real.Exp10M1(R(3), P));
        Assert.Equal(R(3), Real.Log2P1(R(7), P));
        Assert.Equal(R(2), Real.Log10P1(R(99), P));
        AssertAgrees(Real.SinPi(R(1.0 / 6), P), Real.Sin(Real.Divide(Real.Pi(P), R(6), P), P).ToString(), 15);
    }

    [Fact]
    public void ReciprocalTrigAgmSumDot()
    {
        Real x = R(1.375);
        AssertAgrees(Real.Multiply(Real.Sec(x, P), Real.Cos(x, P), P), "1", 70);
        AssertAgrees(Real.Multiply(Real.Cot(x, P), Real.Tan(x, P), P), "1", 70);
        AssertAgrees(Real.Multiply(Real.Coth(x, P), Real.Tanh(x, P), P), "1", 70);
        // Gauss: 1/agm(1, √2) = 0.8346268416740731862814297...
        AssertAgrees(Real.Divide(Real.One, Real.ArithmeticGeometricMean(Real.One, Real.Sqrt(R(2), P), P), P), "0.83462684167407318628142973279904680899399301349034700244982737010368199270952641186969116035127532412906785", 70);
        // Exact sum with cancellation: 2^100 + 1 - 2^100 = 1.
        Real big = R(1).ShiftLeft(100);
        Assert.Equal(R(1).WithPrecision(8), Real.Sum([big, Real.One, -big], 8));
        Assert.Equal(R(14).WithPrecision(8), Real.Dot([R(1), R(2), R(3)], [R(1), R(2), R(3)], 8));
        Real r = Real.Ieee754Remainder(R(17), R(5), out long q, 53);
        Assert.True(r == R(2));
        Assert.Equal(3, q);
        Real ip = Real.Truncate(R(-3.75), out Real fp, 53);
        Assert.True(ip == R(-3));
        Assert.True(fp == R(-0.75));
        Real m = Real.Decompose(R(12), out long e, 53);
        Assert.True(m == R(0.75));
        Assert.Equal(4, e);
    }

    [Fact]
    public void Rationals_RoundTrip()
    {
        var q = new BigRational(new BigInt(1), new BigInt(3));
        Real third = Real.FromRational(q, 64, RealRounding.TowardZero);
        Assert.Equal(Real.Divide(Real.One, R(3), 64, RealRounding.TowardZero), third);
        Assert.Equal(new BigRational(new BigInt(-3), new BigInt(4)), R(-0.75).ToRational());
    }
}

internal static class RealTestExtensions
{
    public static bool IsNegative(this Real r) => r.IsNumber && r.SignBit;

    public static bool IsPositive(this Real r) => r.IsNumber && !r.SignBit;
}
