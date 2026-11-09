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

    private static void AssertAgrees(BigFloat actual, string expected, int digits)
    {
        BigFloat e = BigFloat.Parse(expected, P + 64);
        BigFloat diff = BigFloat.Subtract(actual, e, 64).Abs();
        if (diff.IsZero)
        {
            return;
        }
        BigFloat bound = BigFloat.Multiply(e.Abs(), BigFloat.Exp10(BigFloat.FromDouble(-digits, 64), 64), 64);
        Assert.True(diff.CompareTo(bound) <= 0, $"{actual} vs {expected}");
    }

    private static BigFloat R(double v) => BigFloat.FromDouble(v, P);

    [Fact]
    public void Constants()
    {
        AssertAgrees(BigFloat.EulerGamma(P), EulerDigits, 70);
        AssertAgrees(BigFloat.CatalanConstant(P), CatalanDigits, 70);
        AssertAgrees(BigFloat.Zeta(R(3), P), Zeta3Digits, 70);
        AssertAgrees(BigFloat.AiryAi(BigFloat.Zero, P), Ai0Digits, 60);
        AssertAgrees(BigFloat.BesselJ(0, R(1), P), J0of1Digits, 18);
    }

    [Fact]
    public void Gamma_Identities()
    {
        // Γ(n) = (n-1)! exactly.
        Assert.Equal(BigFloat.FromBigInt(BigInt.Factorial(19), P), BigFloat.Gamma(R(20), P));
        // Γ(1/2)² = π, Γ(x+1) = x Γ(x), lnΓ(1) = lnΓ(2) = +0.
        AssertAgrees(BigFloat.Square(BigFloat.Gamma(R(0.5), P), P), BigFloat.Pi(P).ToString(), 70);
        BigFloat x = R(3.25);
        AssertAgrees(BigFloat.Gamma(R(4.25), P), BigFloat.Multiply(x, BigFloat.Gamma(x, P), P).ToString(), 70);
        Assert.True(BigFloat.LogGamma(BigFloat.One, P).IsZero && BigFloat.LogGamma(R(2), P).IsZero);
        // Reflection: Γ(x) Γ(1-x) = π / sin(πx).
        BigFloat y = R(-2.75);
        BigFloat lhs = BigFloat.Multiply(BigFloat.Gamma(y, P), BigFloat.Gamma(BigFloat.Subtract(BigFloat.One, y, P), P), P);
        AssertAgrees(lhs, BigFloat.Divide(BigFloat.Pi(P), BigFloat.SinPi(y, P), P).ToString(), 65);
        // Signs and poles.
        Assert.True(BigFloat.Gamma(R(-0.5), P).IsNegative());
        Assert.True(BigFloat.Gamma(R(-3), P).IsNaN);
        Assert.True(BigFloat.LogGamma(R(-0.5), P).IsNaN); // Γ(-0.5) < 0
        BigFloat lg = BigFloat.LogGamma(R(-0.5), out int sign, P);
        Assert.Equal(-1, sign);
        AssertAgrees(lg, BigFloat.Log(BigFloat.Gamma(R(-0.5), P).Abs(), P).ToString(), 70);
        // ψ(1) = -γ, ψ(x+1) = ψ(x) + 1/x.
        AssertAgrees(BigFloat.Digamma(BigFloat.One, P), (-BigFloat.EulerGamma(P)).ToString(), 70);
        AssertAgrees(BigFloat.Digamma(R(4.5), P), BigFloat.Add(BigFloat.Digamma(R(3.5), P), BigFloat.Divide(BigFloat.One, R(3.5), P), P).ToString(), 70);
        // B(2, 3) = 1/12 and symmetry.
        AssertAgrees(BigFloat.Beta(R(2), R(3), P), BigFloat.Divide(BigFloat.One, R(12), P).ToString(), 70);
        Assert.Equal(BigFloat.Beta(R(1.5), R(2.25), P), BigFloat.Beta(R(2.25), R(1.5), P));
    }

    [Fact]
    public void ErfZetaEintLi2_Identities()
    {
        BigFloat x = R(0.8125);
        AssertAgrees(BigFloat.Add(BigFloat.Erf(x, P), BigFloat.Erfc(x, P), P), "1", 70);
        Assert.Equal(-BigFloat.Erf(x, P), BigFloat.Erf(-x, P));
        Assert.True(BigFloat.Erfc(R(30), 64).IsPositive() && BigFloat.Erfc(R(30), 64).BinaryExponent < -1290);
        // ζ(2) = π²/6, ζ(4) = π⁴/90, ζ(0) = -1/2, ζ(-1) = -1/12, ζ(-2) = 0.
        BigFloat pi = BigFloat.Pi(P);
        AssertAgrees(BigFloat.Zeta(R(2), P), BigFloat.Divide(BigFloat.Square(pi, P), R(6), P).ToString(), 70);
        AssertAgrees(BigFloat.Zeta(R(4), P), BigFloat.Divide(BigFloat.Pow(pi, 4, P), R(90), P).ToString(), 70);
        Assert.Equal(R(-0.5), BigFloat.Zeta(BigFloat.Zero, P));
        AssertAgrees(BigFloat.Zeta(R(-1), P), BigFloat.Divide(BigFloat.MinusOne, R(12), P).ToString(), 70);
        Assert.True(BigFloat.Zeta(R(-2), P).IsZero);
        // Dilogarithm(1) = π²/6, Dilogarithm(-1) = -π²/12, Dilogarithm(1/2) = π²/12 - ln²2/2.
        BigFloat pi2 = BigFloat.Square(pi, P);
        AssertAgrees(BigFloat.Dilogarithm(BigFloat.One, P), BigFloat.Divide(pi2, R(6), P).ToString(), 70);
        AssertAgrees(BigFloat.Dilogarithm(BigFloat.MinusOne, P), BigFloat.Divide(-pi2, R(12), P).ToString(), 70);
        BigFloat half = BigFloat.Subtract(BigFloat.Divide(pi2, R(12), P), BigFloat.Square(BigFloat.Ln2(P), P).ShiftRight(1), P);
        AssertAgrees(BigFloat.Dilogarithm(R(0.5), P), half.ToString(), 70);
        // Ei'(x) = e^x / x: central difference agrees.
        BigFloat h = BigFloat.FromDouble(1.0 / (1 << 20), P);
        BigFloat d = BigFloat.Divide(BigFloat.Subtract(BigFloat.ExponentialIntegral(BigFloat.Add(x, h, P), P), BigFloat.ExponentialIntegral(BigFloat.Subtract(x, h, P), P), P), h.ShiftLeft(1), P);
        AssertAgrees(d, BigFloat.Divide(BigFloat.Exp(x, P), x, P).ToString(), 10);
    }

    [Fact]
    public void Bessel_Wronskian()
    {
        // J_{n+1}(x) Y_n(x) - J_n(x) Y_{n+1}(x) = 2 / (π x)
        foreach (double xv in (double[])[0.75, 3.5, 17.25])
        {
            BigFloat x = R(xv);
            for (int n = 0; n < 3; n++)
            {
                BigFloat w = BigFloat.Subtract(BigFloat.Multiply(BigFloat.BesselJ(n + 1, x, P), BigFloat.BesselY(n, x, P), P), BigFloat.Multiply(BigFloat.BesselJ(n, x, P), BigFloat.BesselY(n + 1, x, P), P), P);
                AssertAgrees(w, BigFloat.Divide(R(2), BigFloat.Multiply(BigFloat.Pi(P), x, P), P).ToString(), 60);
            }
        }
        Assert.Equal(-BigFloat.BesselJ(1, R(2.5), P), BigFloat.BesselJ(1, R(-2.5), P));
    }

    [Fact]
    public void PiFunctionsAndM1Variants_ExactCases()
    {
        Assert.Equal(BigFloat.One.WithPrecision(53), BigFloat.SinPi(R(0.5), 53));
        Assert.Equal(BigFloat.MinusOne.WithPrecision(53), BigFloat.CosPi(R(3), 53));
        Assert.Equal(BigFloat.One.WithPrecision(53), BigFloat.TanPi(R(0.25), 53));
        Assert.True(BigFloat.SinPi(R(-4), 53) is { IsZero: true, SignBit: true });
        Assert.True(BigFloat.CosPi(R(2.5), 53) is { IsZero: true, SignBit: false });
        Assert.True(BigFloat.TanPi(R(0.5), 53).IsPositiveInfinity);
        Assert.Equal(R(7), BigFloat.Exp2M1(R(3), P));
        Assert.Equal(R(999), BigFloat.Exp10M1(R(3), P));
        Assert.Equal(R(3), BigFloat.Log2P1(R(7), P));
        Assert.Equal(R(2), BigFloat.Log10P1(R(99), P));
        AssertAgrees(BigFloat.SinPi(R(1.0 / 6), P), BigFloat.Sin(BigFloat.Divide(BigFloat.Pi(P), R(6), P), P).ToString(), 15);
    }

    [Fact]
    public void ReciprocalTrigAgmSumDot()
    {
        BigFloat x = R(1.375);
        AssertAgrees(BigFloat.Multiply(BigFloat.Sec(x, P), BigFloat.Cos(x, P), P), "1", 70);
        AssertAgrees(BigFloat.Multiply(BigFloat.Cot(x, P), BigFloat.Tan(x, P), P), "1", 70);
        AssertAgrees(BigFloat.Multiply(BigFloat.Coth(x, P), BigFloat.Tanh(x, P), P), "1", 70);
        // Gauss: 1/agm(1, √2) = 0.8346268416740731862814297...
        AssertAgrees(BigFloat.Divide(BigFloat.One, BigFloat.ArithmeticGeometricMean(BigFloat.One, BigFloat.Sqrt(R(2), P), P), P), "0.83462684167407318628142973279904680899399301349034700244982737010368199270952641186969116035127532412906785", 70);
        // Exact sum with cancellation: 2^100 + 1 - 2^100 = 1.
        BigFloat big = R(1).ShiftLeft(100);
        Assert.Equal(R(1).WithPrecision(8), BigFloat.Sum([big, BigFloat.One, -big], 8));
        Assert.Equal(R(14).WithPrecision(8), BigFloat.Dot([R(1), R(2), R(3)], [R(1), R(2), R(3)], 8));
        BigFloat r = BigFloat.Ieee754Remainder(R(17), R(5), out long q, 53);
        Assert.True(r == R(2));
        Assert.Equal(3, q);
        BigFloat ip = BigFloat.Truncate(R(-3.75), out BigFloat fp, 53);
        Assert.True(ip == R(-3));
        Assert.True(fp == R(-0.75));
        BigFloat m = BigFloat.Decompose(R(12), out long e, 53);
        Assert.True(m == R(0.75));
        Assert.Equal(4, e);
    }

    [Fact]
    public void Rationals_RoundTrip()
    {
        var q = new BigRational(new BigInt(1), new BigInt(3));
        BigFloat third = BigFloat.FromRational(q, 64, RoundingMode.TowardZero);
        Assert.Equal(BigFloat.Divide(BigFloat.One, R(3), 64, RoundingMode.TowardZero), third);
        Assert.Equal(new BigRational(new BigInt(-3), new BigInt(4)), R(-0.75).ToRational());
    }
}

internal static class BigFloatTestExtensions
{
    public static bool IsNegative(this BigFloat r) => r.IsNumber && r.SignBit;

    public static bool IsPositive(this BigFloat r) => r.IsNumber && !r.SignBit;
}
