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
        foreach (RoundingMode rnd in Enum.GetValues<RoundingMode>())
        {
            foreach (int prec in new[] { 2, 53, 300 })
            {
                BigFloat r = BigFloat.RootN(BigFloat.FromDouble(x, prec), k, prec, rnd);
                Assert.Equal(expected, r.ToDouble());
            }
        }
    }

    private static void AssertAgrees(BigFloat actual, string expectedDecimal, int matchDigits)
    {
        BigFloat expected = BigFloat.Parse(expectedDecimal, P + 64);
        BigFloat diff = BigFloat.Subtract(actual, expected, 64).Abs();
        if (diff.IsZero)
        {
            return;
        }
        BigFloat scale = expected.Abs();
        // |diff| <= |expected| * 10^-matchDigits
        BigFloat bound = BigFloat.Multiply(scale, BigFloat.Exp10(BigFloat.FromDouble(-matchDigits, 64), 64), 64);
        Assert.True(diff.CompareTo(bound) <= 0,
            $"values differ beyond 1e-{matchDigits}: {actual} vs {expectedDecimal[..40]}...");
    }

    [Fact]
    public void Constants_MatchPublishedDigits()
    {
        AssertAgrees(BigFloat.Pi(P), PiDigits, 150);
        AssertAgrees(BigFloat.E(P), EDigits, 150);
        AssertAgrees(BigFloat.Ln2(P), Ln2Digits, 150);
        AssertAgrees(BigFloat.Sqrt(BigFloat.Two.WithPrecision(P), P), Sqrt2Digits, 150);
    }

    [Fact]
    public void ExpLog_LowPrecision_MatchesMath()
    {
        var rng = new Random(10);
        for (int i = 0; i < 400; i++)
        {
            double x = (rng.NextDouble() - 0.5) * 40;
            double expected = Math.Exp(x);
            double actual = BigFloat.Exp(BigFloat.FromDouble(x, 53), 53).ToDouble();
            AssertUlps(expected, actual, 1);
            if (expected is > 0 and < double.MaxValue)
            {
                double lx = Math.Abs(x) + 0.001;
                AssertUlps(Math.Log(lx), BigFloat.Log(BigFloat.FromDouble(lx, 53), 53).ToDouble(), 1);
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
            BigFloat x = BigFloat.FromDouble(xd, P);
            BigFloat y = BigFloat.Log(BigFloat.Exp(x, P + 32), P + 32);
            // |y - x| should be ~2^-(P) relative.
            BigFloat diff = BigFloat.Subtract(y, x, 64).Abs();
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
        Assert.Equal(1.0, BigFloat.Exp(BigFloat.Zero).ToDouble());
        Assert.True(BigFloat.Exp(BigFloat.PositiveInfinity).IsPositiveInfinity);
        Assert.True(BigFloat.Exp(BigFloat.NegativeInfinity).IsZero);
        Assert.True(BigFloat.Exp(BigFloat.NaN).IsNaN);
        Assert.True(BigFloat.Log(BigFloat.Zero).IsNegativeInfinity);
        Assert.True(BigFloat.Log(BigFloat.MinusOne).IsNaN);
        Assert.True(BigFloat.Log(BigFloat.One).IsZero);
        Assert.False(BigFloat.Log(BigFloat.One).SignBit);
        Assert.Equal(10.0, BigFloat.Log2(BigFloat.FromDouble(1024.0)).ToDouble());
        Assert.Equal(3.0, BigFloat.Log10(BigFloat.FromDouble(1000.0, 96), 96).ToDouble());
    }

    [Fact]
    public void TrigFunctions_MatchMath()
    {
        var rng = new Random(12);
        for (int i = 0; i < 300; i++)
        {
            double x = (rng.NextDouble() - 0.5) * 20;
            AssertUlps(Math.Sin(x), BigFloat.Sin(BigFloat.FromDouble(x, 53), 53).ToDouble(), 2);
            AssertUlps(Math.Cos(x), BigFloat.Cos(BigFloat.FromDouble(x, 53), 53).ToDouble(), 2);
            AssertUlps(Math.Tan(x), BigFloat.Tan(BigFloat.FromDouble(x, 53), 53).ToDouble(), 2);
            AssertUlps(Math.Atan(x), BigFloat.Atan(BigFloat.FromDouble(x, 53), 53).ToDouble(), 2);
        }
        for (int i = 0; i < 200; i++)
        {
            double x = rng.NextDouble() * 1.99 - 0.995;
            AssertUlps(Math.Asin(x), BigFloat.Asin(BigFloat.FromDouble(x, 53), 53).ToDouble(), 2);
            AssertUlps(Math.Acos(x), BigFloat.Acos(BigFloat.FromDouble(x, 53), 53).ToDouble(), 2);
        }
    }

    [Fact]
    public void Trig_Identity_HighPrecision()
    {
        var rng = new Random(13);
        for (int i = 0; i < 15; i++)
        {
            double xd = (rng.NextDouble() - 0.5) * 6;
            BigFloat x = BigFloat.FromDouble(xd, P);
            BigFloat s = BigFloat.Sin(x, P + 16);
            BigFloat c = BigFloat.Cos(x, P + 16);
            BigFloat sum = BigFloat.Add(BigFloat.Square(s, P + 16), BigFloat.Square(c, P + 16), 64);
            BigFloat diff = BigFloat.Subtract(sum, BigFloat.One, 64).Abs();
            Assert.True(diff.IsZero || diff.BinaryExponent <= -P + 8, $"sin²+cos² != 1 at x={xd}");
        }
    }

    [Fact]
    public void Trig_LargeArgumentReduction()
    {
        // sin(10^30) needs π to ~100 extra digits; verify vs high-precision identity.
        BigFloat big = BigFloat.Parse("1e30", 256);
        BigFloat s = BigFloat.Sin(big, 128);
        BigFloat c = BigFloat.Cos(big, 128);
        BigFloat sum = BigFloat.Add(BigFloat.Square(s, 160), BigFloat.Square(c, 160), 64);
        BigFloat diff = BigFloat.Subtract(sum, BigFloat.One, 64).Abs();
        Assert.True(diff.IsZero || diff.BinaryExponent <= -120);
        // Cross-check the double value.
        AssertUlps(Math.Sin(1e30), BigFloat.Sin(BigFloat.FromDouble(1e30, 53), 53).ToDouble(), 4);
    }

    [Fact]
    public void Atan2_Quadrants()
    {
        AssertUlps(Math.Atan2(1, 1), BigFloat.Atan2(BigFloat.One, BigFloat.One, 53).ToDouble(), 2);
        AssertUlps(Math.Atan2(1, -1), BigFloat.Atan2(BigFloat.One, BigFloat.MinusOne, 53).ToDouble(), 2);
        AssertUlps(Math.Atan2(-1, 1), BigFloat.Atan2(BigFloat.MinusOne, BigFloat.One, 53).ToDouble(), 2);
        AssertUlps(Math.Atan2(-1, -1), BigFloat.Atan2(BigFloat.MinusOne, BigFloat.MinusOne, 53).ToDouble(), 2);
        Assert.True(BigFloat.Atan2(BigFloat.Zero, BigFloat.One).IsZero);
        AssertUlps(Math.PI, BigFloat.Atan2(BigFloat.Zero, BigFloat.MinusOne, 53).ToDouble(), 2);
        AssertUlps(Math.PI / 2, BigFloat.Atan2(BigFloat.One, BigFloat.Zero, 53).ToDouble(), 2);
        var rng = new Random(14);
        for (int i = 0; i < 200; i++)
        {
            double y = (rng.NextDouble() - 0.5) * 100, x = (rng.NextDouble() - 0.5) * 100;
            AssertUlps(Math.Atan2(y, x), BigFloat.Atan2(BigFloat.FromDouble(y, 53), BigFloat.FromDouble(x, 53), 53).ToDouble(), 2);
        }
    }

    [Fact]
    public void Hyperbolic_MatchMath()
    {
        var rng = new Random(15);
        for (int i = 0; i < 200; i++)
        {
            double x = (rng.NextDouble() - 0.5) * 12;
            AssertUlps(Math.Sinh(x), BigFloat.Sinh(BigFloat.FromDouble(x, 53), 53).ToDouble(), 2);
            AssertUlps(Math.Cosh(x), BigFloat.Cosh(BigFloat.FromDouble(x, 53), 53).ToDouble(), 2);
            AssertUlps(Math.Tanh(x), BigFloat.Tanh(BigFloat.FromDouble(x, 53), 53).ToDouble(), 2);
            AssertUlps(Math.Asinh(x), BigFloat.Asinh(BigFloat.FromDouble(x, 53), 53).ToDouble(), 2);
        }
        for (int i = 0; i < 100; i++)
        {
            double x = 1 + rng.NextDouble() * 20;
            AssertUlps(Math.Acosh(x), BigFloat.Acosh(BigFloat.FromDouble(x, 53), 53).ToDouble(), 2);
            double t = rng.NextDouble() * 1.98 - 0.99;
            AssertUlps(Math.Atanh(t), BigFloat.Atanh(BigFloat.FromDouble(t, 53), 53).ToDouble(), 2);
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
            AssertUlps(Math.Pow(x, y), BigFloat.Pow(BigFloat.FromDouble(x, 53), BigFloat.FromDouble(y, 53), 53).ToDouble(), 2);
        }
        // Special cases.
        Assert.Equal(1.0, BigFloat.Pow(BigFloat.NaN, BigFloat.Zero).ToDouble());
        Assert.Equal(1.0, BigFloat.Pow(BigFloat.One, BigFloat.NaN).ToDouble());
        Assert.True(BigFloat.Pow(BigFloat.MinusOne, BigFloat.FromDouble(0.5)).IsNaN);
        Assert.Equal(-8.0, BigFloat.Pow(BigFloat.FromDouble(-2.0), BigFloat.FromDouble(3.0), 53).ToDouble());
        Assert.Equal(0.25, BigFloat.Pow(BigFloat.FromDouble(-2.0), BigFloat.FromDouble(-2.0), 53).ToDouble());
        Assert.True(BigFloat.Pow(BigFloat.Zero, BigFloat.FromDouble(-1.0)).IsPositiveInfinity);
        Assert.True(BigFloat.Pow(BigFloat.NegativeZero, BigFloat.FromDouble(-3.0)).IsNegativeInfinity);
        Assert.Equal(1.0, BigFloat.Pow(BigFloat.PositiveInfinity, BigFloat.Zero).ToDouble());
    }

    [Fact]
    public void PowInt_ExactCases()
    {
        Assert.Equal(1024.0, BigFloat.Pow(BigFloat.Two, 10, 53).ToDouble());
        Assert.Equal(Math.Pow(3, 20), BigFloat.Pow(BigFloat.FromDouble(3, 53), 20, 53).ToDouble());
        Assert.Equal(Math.Pow(2, -10), BigFloat.Pow(BigFloat.Two, -10, 53).ToDouble());
        Assert.Equal(-27.0, BigFloat.Pow(BigFloat.FromDouble(-3, 53), 3, 53).ToDouble());
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
            AssertNearestRoot(BigFloat.Hypot(BigFloat.FromDouble(x, 53), BigFloat.FromDouble(y, 53), 53).ToDouble(), 2, sumSq, e2, $"hypot({x:R}, {y:R})");
            double cbrt = BigFloat.Cbrt(BigFloat.FromDouble(x, 53), 53).ToDouble();
            Assert.Equal(Math.Sign(x), Math.Sign(cbrt));
            AssertNearestRoot(Math.Abs(cbrt), 3, mx, ex, $"cbrt({x:R})");
        }
        Assert.Equal(2.0, BigFloat.RootN(BigFloat.FromDouble(32.0, 53), 5, 53).ToDouble());
        Assert.Equal(-2.0, BigFloat.Cbrt(BigFloat.FromDouble(-8.0, 53), 53).ToDouble());
        Assert.True(BigFloat.RootN(BigFloat.FromDouble(-8.0, 53), 4, 53).IsNaN);
        AssertUlps(1 / Math.Sqrt(2), BigFloat.ReciprocalSqrt(BigFloat.FromDouble(2.0, 53), 53).ToDouble(), 1);
    }

    [Fact]
    public void Expm1Log1p_AccurateNearZero()
    {
        double x = 1e-30;
        BigFloat r = BigFloat.ExpM1(BigFloat.FromDouble(x, 53), 53);
        // expm1(1e-30) ≈ 1e-30 exactly at double precision.
        Assert.Equal(x, r.ToDouble());
        BigFloat l = BigFloat.LogP1(BigFloat.FromDouble(x, 53), 53);
        Assert.Equal(x, l.ToDouble());
        var rng = new Random(18);
        for (int i = 0; i < 100; i++)
        {
            double t = (rng.NextDouble() - 0.5) * 2;
            BigFloat xt = BigFloat.FromDouble(t, 53);
            // Reference via high-precision arithmetic (double Math.Exp(t)-1 suffers cancellation).
            double expm1Ref = BigFloat.Subtract(BigFloat.Exp(xt, 120), BigFloat.One, 53).ToDouble();
            AssertUlps(expm1Ref, BigFloat.ExpM1(xt, 53).ToDouble(), 1);
            if (t > -1)
            {
                double log1pRef = BigFloat.Log(BigFloat.Add(BigFloat.One, xt, 120), 53).ToDouble();
                AssertUlps(log1pRef, BigFloat.LogP1(xt, 53).ToDouble(), 1);
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
            BigFloat r = BigFloat.FromDouble(d, 53);
            BigFloat back = BigFloat.Parse(r.ToString(), 53);
            Assert.Equal(r.ToDouble(), back.ToDouble());
        }
        // High-precision round trip.
        BigFloat pi = BigFloat.Pi(300);
        BigFloat back2 = BigFloat.Parse(pi.ToString(), 300);
        Assert.True(BigFloat.Subtract(pi, back2, 64).IsZero);
        Assert.Equal("NaN", BigFloat.NaN.ToString());
        Assert.Equal("Infinity", BigFloat.PositiveInfinity.ToString());
        Assert.Equal("-0", BigFloat.NegativeZero.ToString());
        Assert.True(BigFloat.Parse("-inf").IsNegativeInfinity);
        Assert.True(BigFloat.Parse("nan").IsNaN);
    }

    [Fact]
    public void TinyArguments_FastPath()
    {
        // Functions ≈ x for very small x must terminate quickly; under round-to-nearest
        // the result equals x, and directed modes must land on the correct side.
        BigFloat tiny = BigFloat.FromDouble(1.0, 200).ShiftRight(100000); // 2^-100000
        Assert.Equal(0, BigFloat.Sin(tiny, 64).CompareTo(tiny));
        Assert.True(BigFloat.Sin(tiny, 64, RoundingMode.TowardZero).CompareTo(tiny) < 0);  // sin x < x
        Assert.True(BigFloat.Tan(tiny, 64, RoundingMode.AwayFromZero).CompareTo(tiny) > 0); // tan x > x
        Assert.True(BigFloat.Atan(tiny, 64, RoundingMode.TowardZero).CompareTo(tiny) < 0);
        Assert.True(BigFloat.Sinh(tiny, 64, RoundingMode.AwayFromZero).CompareTo(tiny) > 0);
        Assert.True(BigFloat.Tanh(tiny, 64, RoundingMode.TowardZero).CompareTo(tiny) < 0);
        Assert.True(BigFloat.Asinh(tiny, 64, RoundingMode.TowardZero).CompareTo(tiny) < 0);
        Assert.True(BigFloat.Atanh(tiny, 64, RoundingMode.AwayFromZero).CompareTo(tiny) > 0);
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
