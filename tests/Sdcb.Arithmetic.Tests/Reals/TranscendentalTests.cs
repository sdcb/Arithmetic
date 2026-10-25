using System.Numerics;
using Xunit;

namespace Sdcb.Arithmetic.Tests;

/// <summary>High-precision checks against published constant/function digits and identities.</summary>
public class TranscendentalTests
{
    private const string PiDigits =
        "3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303820";
    private const string EDigits =
        "2.7182818284590452353602874713526624977572470936999595749669676277240766303535475945713821785251664274274663919320030599218174135966290435729003342952605956307381323286279434907632338298807531952510190";
    private const string Ln2Digits =
        "0.6931471805599453094172321214581765680755001343602552541206800094933936219696947156058633269964186875420014810205706857336855202357581305570326707516350759619307275708283714351903070386238916734711233";
    private const string Sqrt2Digits =
        "1.4142135623730950488016887242096980785696718753769480731766797379907324784621070388503875343276415727350138462309122970249248360558507372126441214970999358314132226659275055927557999505011527820605715";

    private const int P = 512;

    [Theory]
    [InlineData(256.0, -8L, 0.5)]
    [InlineData(256.0, -4L, 0.25)]
    [InlineData(0.0625, -4L, 2.0)]
    [InlineData(-8.0, -3L, -0.5)]
    public void RootN_NegativeK_ExactPowerOfTwo(double x, long k, double expected)
    {
        foreach (RealRounding rnd in Enum.GetValues<RealRounding>())
        {
            foreach (int prec in new[] { 2, 53, 300 })
            {
                Real r = Real.RootN(Real.FromDouble(x, prec), k, prec, rnd);
                Assert.Equal(expected, r.ToDouble());
            }
        }
    }

    private static void AssertAgrees(Real actual, string expectedDecimal, int matchDigits)
    {
        Real expected = Real.Parse(expectedDecimal, P + 64);
        Real diff = Real.Subtract(actual, expected, 64).Abs();
        if (diff.IsZero)
        {
            return;
        }
        Real scale = expected.Abs();
        // |diff| <= |expected| * 10^-matchDigits
        Real bound = Real.Multiply(scale, Real.Exp10(Real.FromDouble(-matchDigits, 64), 64), 64);
        Assert.True(diff.CompareTo(bound) <= 0,
            $"values differ beyond 1e-{matchDigits}: {actual} vs {expectedDecimal[..40]}...");
    }

    [Fact]
    public void Constants_MatchPublishedDigits()
    {
        AssertAgrees(Real.Pi(P), PiDigits, 150);
        AssertAgrees(Real.E(P), EDigits, 150);
        AssertAgrees(Real.Ln2(P), Ln2Digits, 150);
        AssertAgrees(Real.Sqrt(Real.Two.WithPrecision(P), P), Sqrt2Digits, 150);
    }

    [Fact]
    public void ExpLog_LowPrecision_MatchesMath()
    {
        var rng = new Random(10);
        for (int i = 0; i < 400; i++)
        {
            double x = (rng.NextDouble() - 0.5) * 40;
            double expected = Math.Exp(x);
            double actual = Real.Exp(Real.FromDouble(x, 53), 53).ToDouble();
            AssertUlps(expected, actual, 1);
            if (expected is > 0 and < double.MaxValue)
            {
                double lx = Math.Abs(x) + 0.001;
                AssertUlps(Math.Log(lx), Real.Log(Real.FromDouble(lx, 53), 53).ToDouble(), 1);
            }
        }
    }

    [Fact]
    public void ExpLog_Identities_HighPrecision()
    {
        var rng = new Random(11);
        for (int i = 0; i < 25; i++)
        {
            double xd = (rng.NextDouble() - 0.5) * 20;
            Real x = Real.FromDouble(xd, P);
            Real y = Real.Log(Real.Exp(x, P + 32), P + 32);
            // |y - x| should be ~2^-(P) relative.
            Real diff = Real.Subtract(y, x, 64).Abs();
            if (!diff.IsZero && !x.IsZero)
            {
                Assert.True(diff.BinaryExponent <= x.BinaryExponent - P + 8,
                    $"log(exp(x)) drifted: x={xd}, diffExp={diff.BinaryExponent}");
            }
        }
    }

    [Fact]
    public void Exp_SpecialCases()
    {
        Assert.Equal(1.0, Real.Exp(Real.Zero).ToDouble());
        Assert.True(Real.Exp(Real.PositiveInfinity).IsPositiveInfinity);
        Assert.True(Real.Exp(Real.NegativeInfinity).IsZero);
        Assert.True(Real.Exp(Real.NaN).IsNaN);
        Assert.True(Real.Log(Real.Zero).IsNegativeInfinity);
        Assert.True(Real.Log(Real.MinusOne).IsNaN);
        Assert.True(Real.Log(Real.One).IsZero);
        Assert.False(Real.Log(Real.One).SignBit);
        Assert.Equal(10.0, Real.Log2(Real.FromDouble(1024.0)).ToDouble());
        Assert.Equal(3.0, Real.Log10(Real.FromDouble(1000.0, 96), 96).ToDouble());
    }

    [Fact]
    public void TrigFunctions_MatchMath()
    {
        var rng = new Random(12);
        for (int i = 0; i < 300; i++)
        {
            double x = (rng.NextDouble() - 0.5) * 20;
            AssertUlps(Math.Sin(x), Real.Sin(Real.FromDouble(x, 53), 53).ToDouble(), 2);
            AssertUlps(Math.Cos(x), Real.Cos(Real.FromDouble(x, 53), 53).ToDouble(), 2);
            AssertUlps(Math.Tan(x), Real.Tan(Real.FromDouble(x, 53), 53).ToDouble(), 2);
            AssertUlps(Math.Atan(x), Real.Atan(Real.FromDouble(x, 53), 53).ToDouble(), 2);
        }
        for (int i = 0; i < 200; i++)
        {
            double x = rng.NextDouble() * 1.99 - 0.995;
            AssertUlps(Math.Asin(x), Real.Asin(Real.FromDouble(x, 53), 53).ToDouble(), 2);
            AssertUlps(Math.Acos(x), Real.Acos(Real.FromDouble(x, 53), 53).ToDouble(), 2);
        }
    }

    [Fact]
    public void Trig_Identity_HighPrecision()
    {
        var rng = new Random(13);
        for (int i = 0; i < 15; i++)
        {
            double xd = (rng.NextDouble() - 0.5) * 6;
            Real x = Real.FromDouble(xd, P);
            Real s = Real.Sin(x, P + 16);
            Real c = Real.Cos(x, P + 16);
            Real sum = Real.Add(Real.Square(s, P + 16), Real.Square(c, P + 16), 64);
            Real diff = Real.Subtract(sum, Real.One, 64).Abs();
            Assert.True(diff.IsZero || diff.BinaryExponent <= -P + 8, $"sin²+cos² != 1 at x={xd}");
        }
    }

    [Fact]
    public void Trig_LargeArgumentReduction()
    {
        // sin(10^30) needs π to ~100 extra digits; verify vs high-precision identity.
        Real big = Real.Parse("1e30", 256);
        Real s = Real.Sin(big, 128);
        Real c = Real.Cos(big, 128);
        Real sum = Real.Add(Real.Square(s, 160), Real.Square(c, 160), 64);
        Real diff = Real.Subtract(sum, Real.One, 64).Abs();
        Assert.True(diff.IsZero || diff.BinaryExponent <= -120);
        // Cross-check the double value.
        AssertUlps(Math.Sin(1e30), Real.Sin(Real.FromDouble(1e30, 53), 53).ToDouble(), 4);
    }

    [Fact]
    public void Atan2_Quadrants()
    {
        AssertUlps(Math.Atan2(1, 1), Real.Atan2(Real.One, Real.One, 53).ToDouble(), 2);
        AssertUlps(Math.Atan2(1, -1), Real.Atan2(Real.One, Real.MinusOne, 53).ToDouble(), 2);
        AssertUlps(Math.Atan2(-1, 1), Real.Atan2(Real.MinusOne, Real.One, 53).ToDouble(), 2);
        AssertUlps(Math.Atan2(-1, -1), Real.Atan2(Real.MinusOne, Real.MinusOne, 53).ToDouble(), 2);
        Assert.True(Real.Atan2(Real.Zero, Real.One).IsZero);
        AssertUlps(Math.PI, Real.Atan2(Real.Zero, Real.MinusOne, 53).ToDouble(), 2);
        AssertUlps(Math.PI / 2, Real.Atan2(Real.One, Real.Zero, 53).ToDouble(), 2);
        var rng = new Random(14);
        for (int i = 0; i < 200; i++)
        {
            double y = (rng.NextDouble() - 0.5) * 100, x = (rng.NextDouble() - 0.5) * 100;
            AssertUlps(Math.Atan2(y, x), Real.Atan2(Real.FromDouble(y, 53), Real.FromDouble(x, 53), 53).ToDouble(), 2);
        }
    }

    [Fact]
    public void Hyperbolic_MatchMath()
    {
        var rng = new Random(15);
        for (int i = 0; i < 200; i++)
        {
            double x = (rng.NextDouble() - 0.5) * 12;
            AssertUlps(Math.Sinh(x), Real.Sinh(Real.FromDouble(x, 53), 53).ToDouble(), 2);
            AssertUlps(Math.Cosh(x), Real.Cosh(Real.FromDouble(x, 53), 53).ToDouble(), 2);
            AssertUlps(Math.Tanh(x), Real.Tanh(Real.FromDouble(x, 53), 53).ToDouble(), 2);
            AssertUlps(Math.Asinh(x), Real.Asinh(Real.FromDouble(x, 53), 53).ToDouble(), 2);
        }
        for (int i = 0; i < 100; i++)
        {
            double x = 1 + rng.NextDouble() * 20;
            AssertUlps(Math.Acosh(x), Real.Acosh(Real.FromDouble(x, 53), 53).ToDouble(), 2);
            double t = rng.NextDouble() * 1.98 - 0.99;
            AssertUlps(Math.Atanh(t), Real.Atanh(Real.FromDouble(t, 53), 53).ToDouble(), 2);
        }
    }

    [Fact]
    public void Pow_MatchesMath()
    {
        var rng = new Random(16);
        for (int i = 0; i < 300; i++)
        {
            double x = rng.NextDouble() * 10 + 0.01;
            double y = (rng.NextDouble() - 0.5) * 20;
            AssertUlps(Math.Pow(x, y), Real.Pow(Real.FromDouble(x, 53), Real.FromDouble(y, 53), 53).ToDouble(), 2);
        }
        // Special cases.
        Assert.Equal(1.0, Real.Pow(Real.NaN, Real.Zero).ToDouble());
        Assert.Equal(1.0, Real.Pow(Real.One, Real.NaN).ToDouble());
        Assert.True(Real.Pow(Real.MinusOne, Real.FromDouble(0.5)).IsNaN);
        Assert.Equal(-8.0, Real.Pow(Real.FromDouble(-2.0), Real.FromDouble(3.0), 53).ToDouble());
        Assert.Equal(0.25, Real.Pow(Real.FromDouble(-2.0), Real.FromDouble(-2.0), 53).ToDouble());
        Assert.True(Real.Pow(Real.Zero, Real.FromDouble(-1.0)).IsPositiveInfinity);
        Assert.True(Real.Pow(Real.NegativeZero, Real.FromDouble(-3.0)).IsNegativeInfinity);
        Assert.Equal(1.0, Real.Pow(Real.PositiveInfinity, Real.Zero).ToDouble());
    }

    [Fact]
    public void PowInt_ExactCases()
    {
        Assert.Equal(1024.0, Real.Pow(Real.Two, 10, 53).ToDouble());
        Assert.Equal(Math.Pow(3, 20), Real.Pow(Real.FromDouble(3, 53), 20, 53).ToDouble());
        Assert.Equal(Math.Pow(2, -10), Real.Pow(Real.Two, -10, 53).ToDouble());
        Assert.Equal(-27.0, Real.Pow(Real.FromDouble(-3, 53), 3, 53).ToDouble());
    }

    [Fact]
    public void HypotCbrtRootN()
    {
        var rng = new Random(17);
        for (int i = 0; i < 200; i++)
        {
            double x = (rng.NextDouble() - 0.5) * 1e10;
            double y = (rng.NextDouble() - 0.5) * 1e10;
            // Exact oracles: libm cbrt (e.g. glibc) is not correctly rounded, and x*x + y*y rounds.
            var (mx, ex) = Decompose(x);
            var (my, ey) = Decompose(y);
            int e2 = Math.Min(2 * ex, 2 * ey);
            BigInteger sumSq = (mx * mx << (2 * ex - e2)) + (my * my << (2 * ey - e2));
            AssertNearestRoot(Real.Hypot(Real.FromDouble(x, 53), Real.FromDouble(y, 53), 53).ToDouble(), 2, sumSq, e2, $"hypot({x:R}, {y:R})");
            double cbrt = Real.Cbrt(Real.FromDouble(x, 53), 53).ToDouble();
            Assert.Equal(Math.Sign(x), Math.Sign(cbrt));
            AssertNearestRoot(Math.Abs(cbrt), 3, mx, ex, $"cbrt({x:R})");
        }
        Assert.Equal(2.0, Real.RootN(Real.FromDouble(32.0, 53), 5, 53).ToDouble());
        Assert.Equal(-2.0, Real.Cbrt(Real.FromDouble(-8.0, 53), 53).ToDouble());
        Assert.True(Real.RootN(Real.FromDouble(-8.0, 53), 4, 53).IsNaN);
        AssertUlps(1 / Math.Sqrt(2), Real.ReciprocalSqrt(Real.FromDouble(2.0, 53), 53).ToDouble(), 1);
    }

    [Fact]
    public void Expm1Log1p_AccurateNearZero()
    {
        double x = 1e-30;
        Real r = Real.ExpM1(Real.FromDouble(x, 53), 53);
        // expm1(1e-30) ≈ 1e-30 exactly at double precision.
        Assert.Equal(x, r.ToDouble());
        Real l = Real.LogP1(Real.FromDouble(x, 53), 53);
        Assert.Equal(x, l.ToDouble());
        var rng = new Random(18);
        for (int i = 0; i < 100; i++)
        {
            double t = (rng.NextDouble() - 0.5) * 2;
            Real xt = Real.FromDouble(t, 53);
            // Reference via high-precision arithmetic (double Math.Exp(t)-1 suffers cancellation).
            double expm1Ref = Real.Subtract(Real.Exp(xt, 120), Real.One, 53).ToDouble();
            AssertUlps(expm1Ref, Real.ExpM1(xt, 53).ToDouble(), 1);
            if (t > -1)
            {
                double log1pRef = Real.Log(Real.Add(Real.One, xt, 120), 53).ToDouble();
                AssertUlps(log1pRef, Real.LogP1(xt, 53).ToDouble(), 1);
            }
        }
    }

    [Fact]
    public void StringRoundTrip()
    {
        var rng = new Random(19);
        for (int i = 0; i < 300; i++)
        {
            double d = Math.ScaleB(rng.NextDouble() * 2 - 1, rng.Next(-80, 80));
            Real r = Real.FromDouble(d, 53);
            Real back = Real.Parse(r.ToString(), 53);
            Assert.Equal(r.ToDouble(), back.ToDouble());
        }
        // High-precision round trip.
        Real pi = Real.Pi(300);
        Real back2 = Real.Parse(pi.ToString(), 300);
        Assert.True(Real.Subtract(pi, back2, 64).IsZero);
        Assert.Equal("NaN", Real.NaN.ToString());
        Assert.Equal("Infinity", Real.PositiveInfinity.ToString());
        Assert.Equal("-0", Real.NegativeZero.ToString());
        Assert.True(Real.Parse("-inf").IsNegativeInfinity);
        Assert.True(Real.Parse("nan").IsNaN);
    }

    [Fact]
    public void TinyArguments_FastPath()
    {
        // Functions ≈ x for very small x must terminate quickly; under round-to-nearest
        // the result equals x, and directed modes must land on the correct side.
        Real tiny = Real.FromDouble(1.0, 200).ShiftRight(100000); // 2^-100000
        Assert.Equal(0, Real.Sin(tiny, 64).CompareTo(tiny));
        Assert.True(Real.Sin(tiny, 64, RealRounding.TowardZero).CompareTo(tiny) < 0);  // sin x < x
        Assert.True(Real.Tan(tiny, 64, RealRounding.AwayFromZero).CompareTo(tiny) > 0); // tan x > x
        Assert.True(Real.Atan(tiny, 64, RealRounding.TowardZero).CompareTo(tiny) < 0);
        Assert.True(Real.Sinh(tiny, 64, RealRounding.AwayFromZero).CompareTo(tiny) > 0);
        Assert.True(Real.Tanh(tiny, 64, RealRounding.TowardZero).CompareTo(tiny) < 0);
        Assert.True(Real.Asinh(tiny, 64, RealRounding.TowardZero).CompareTo(tiny) < 0);
        Assert.True(Real.Atanh(tiny, 64, RealRounding.AwayFromZero).CompareTo(tiny) > 0);
    }

    /// <summary>|d| = m * 2^e exactly, for finite d.</summary>
    private static (BigInteger m, int e) Decompose(double d)
    {
        long bits = BitConverter.DoubleToInt64Bits(Math.Abs(d));
        int be = (int)(bits >> 52);
        long m = bits & 0xFFFFFFFFFFFFFL;
        return be == 0 ? (m, -1074) : (m | 1L << 52, be - 1075);
    }

    /// <summary>
    /// Asserts r is the round-to-nearest n-th root of t * 2^te, i.e. t * 2^te lies between
    /// (r - ulp⁻/2)^n and (r + ulp⁺/2)^n, using exact integer arithmetic.
    /// </summary>
    private static void AssertNearestRoot(double r, int n, BigInteger t, int te, string what)
    {
        Assert.True(r > 0 && double.IsNormal(r), $"{what}: unexpected result {r:R}");
        var (m, e) = Decompose(r);
        // In units of 2^(e-2): r = 4m, half the upper gap = 2, half the lower gap = 2 (or 1 at a power of two).
        BigInteger lo = 4 * m - (m == BigInteger.One << 52 ? 1 : 2), hi = 4 * m + 2;
        BigInteger loN = BigInteger.Pow(lo, n), hiN = BigInteger.Pow(hi, n);
        int shift = n * (e - 2) - te;
        if (shift >= 0)
        {
            loN <<= shift;
            hiN <<= shift;
        }
        else
        {
            t <<= -shift;
        }
        Assert.True(loN <= t && t <= hiN, $"{what}: {r:R} is not the correctly rounded root");
    }

    private static void AssertUlps(double expected, double actual, int maxUlps)
    {
        if (double.IsNaN(expected))
        {
            Assert.True(double.IsNaN(actual));
            return;
        }
        if (double.IsInfinity(expected) || expected == actual)
        {
            Assert.Equal(expected, actual);
            return;
        }
        long ulps = Math.Abs(
            BitConverter.DoubleToInt64Bits(expected) - BitConverter.DoubleToInt64Bits(actual));
        Assert.True(ulps <= maxUlps, $"expected {expected:R}, actual {actual:R}, ulps={ulps}");
    }
}
