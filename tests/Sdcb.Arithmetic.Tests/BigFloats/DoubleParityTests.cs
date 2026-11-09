using Xunit;

namespace Sdcb.Arithmetic.Tests;

/// <summary>
/// At 53-bit precision every correctly-rounded BigFloat operation must agree
/// bit-for-bit with IEEE double arithmetic.
/// </summary>
public class DoubleParityTests
{
    private static double RandomDouble(Random rng)
    {
        // Finite values across a wide but overflow-safe exponent range.
        double m = rng.NextDouble() * 2 - 1;
        int e = rng.Next(-300, 300);
        return Math.ScaleB(m, e);
    }

    [Fact]
    public void AddSubMulDiv_MatchDouble_At53Bits()
    {
        var rng = new Random(1);
        for (int i = 0; i < 5000; i++)
        {
            double a = RandomDouble(rng), b = RandomDouble(rng);
            BigFloat ra = BigFloat.FromDouble(a, 53), rb = BigFloat.FromDouble(b, 53);
            Assert.Equal(a + b, BigFloat.Add(ra, rb, 53).ToDouble());
            Assert.Equal(a - b, BigFloat.Subtract(ra, rb, 53).ToDouble());
            Assert.Equal(a * b, BigFloat.Multiply(ra, rb, 53).ToDouble());
            if (b != 0)
            {
                Assert.Equal(a / b, BigFloat.Divide(ra, rb, 53).ToDouble());
            }
        }
    }

    [Fact]
    public void Sqrt_MatchesDouble_At53Bits()
    {
        var rng = new Random(2);
        for (int i = 0; i < 2000; i++)
        {
            double a = Math.Abs(RandomDouble(rng));
            Assert.Equal(Math.Sqrt(a), BigFloat.Sqrt(BigFloat.FromDouble(a, 53), 53).ToDouble());
        }
    }

    [Fact]
    public void Fma_MatchesDouble_At53Bits()
    {
        var rng = new Random(3);
        for (int i = 0; i < 2000; i++)
        {
            double a = RandomDouble(rng), b = RandomDouble(rng), c = RandomDouble(rng);
            Assert.Equal(Math.FusedMultiplyAdd(a, b, c),
                BigFloat.FusedMultiplyAdd(BigFloat.FromDouble(a, 53), BigFloat.FromDouble(b, 53), BigFloat.FromDouble(c, 53), 53).ToDouble());
        }
    }

    [Fact]
    public void DoubleRoundTrip_Exact()
    {
        var rng = new Random(4);
        for (int i = 0; i < 3000; i++)
        {
            double d = RandomDouble(rng);
            Assert.Equal(d, BigFloat.FromDouble(d).ToDouble());
        }
        Assert.True(double.IsNaN(BigFloat.NaN.ToDouble()));
        Assert.Equal(double.PositiveInfinity, BigFloat.PositiveInfinity.ToDouble());
        Assert.Equal(double.NegativeInfinity, BigFloat.NegativeInfinity.ToDouble());
        Assert.Equal(-0.0, BigFloat.NegativeZero.ToDouble());
        Assert.True(double.IsNegative(BigFloat.NegativeZero.ToDouble()));
    }

    [Fact]
    public void RoundingModes_DirectedDivision()
    {
        // 1/3 at 8 bits: 0.0101010101...b → check all modes.
        BigFloat oneThird(RoundingMode r) => BigFloat.Divide(BigFloat.One, BigFloat.FromBigInt(3), 8, r);
        BigFloat down = oneThird(RoundingMode.TowardZero);
        BigFloat up = oneThird(RoundingMode.TowardPositive);
        Assert.True(down.CompareTo(up) < 0);
        Assert.Equal(down, oneThird(RoundingMode.TowardNegative));
        Assert.Equal(up, oneThird(RoundingMode.AwayFromZero));
        // Negative: directions flip.
        BigFloat m = BigFloat.Divide(BigFloat.MinusOne, BigFloat.FromBigInt(3), 8, RoundingMode.TowardPositive);
        Assert.Equal((-down).ToDouble(), m.ToDouble());
    }

    [Fact]
    public void SpecialValues_Arithmetic()
    {
        Assert.True((BigFloat.PositiveInfinity + BigFloat.NegativeInfinity).IsNaN);
        Assert.True((BigFloat.Zero * BigFloat.PositiveInfinity).IsNaN);
        Assert.True((BigFloat.PositiveInfinity / BigFloat.PositiveInfinity).IsNaN);
        Assert.True((BigFloat.One / BigFloat.Zero).IsPositiveInfinity);
        Assert.True((BigFloat.MinusOne / BigFloat.Zero).IsNegativeInfinity);
        Assert.True((BigFloat.One / BigFloat.NegativeZero).IsNegativeInfinity);
        Assert.True(BigFloat.Sqrt(BigFloat.MinusOne).IsNaN);
        Assert.True(BigFloat.Sqrt(BigFloat.NegativeZero).IsZero);
        Assert.True(BigFloat.Sqrt(BigFloat.NegativeZero).SignBit);
        // 0 + -0 = +0 except toward negative.
        Assert.False((BigFloat.Zero + BigFloat.NegativeZero).SignBit);
        Assert.True(BigFloat.Add(BigFloat.Zero, BigFloat.NegativeZero, 53, RoundingMode.TowardNegative).SignBit);
        // NaN comparisons.
        Assert.False(BigFloat.NaN == BigFloat.NaN);
        Assert.True(BigFloat.NaN != BigFloat.NaN);
        Assert.False(BigFloat.NaN < BigFloat.One);
        Assert.False(BigFloat.NaN >= BigFloat.One);
        Assert.True(BigFloat.NaN.Equals(BigFloat.NaN)); // bitwise Equals
    }

    [Fact]
    public void Comparison_TotalOrder()
    {
        Assert.True(BigFloat.NegativeZero.CompareTo(BigFloat.Zero) < 0);
        Assert.True(BigFloat.NaN.CompareTo(BigFloat.NegativeInfinity) < 0);
        Assert.True(BigFloat.NegativeInfinity < BigFloat.MinusOne);
        Assert.True(BigFloat.PositiveInfinity > BigFloat.FromBigInt(BigInt.Parse("1" + new string('0', 100))));
        Assert.True(BigFloat.NegativeZero == BigFloat.Zero); // IEEE equality
    }

    [Fact]
    public void NextAboveBelow()
    {
        BigFloat one = BigFloat.FromDouble(1.0, 53);
        BigFloat up = BigFloat.BitIncrement(one);
        Assert.Equal(Math.BitIncrement(1.0), up.ToDouble());
        Assert.Equal(1.0, BigFloat.BitDecrement(up).ToDouble());
        BigFloat down = BigFloat.BitDecrement(one);
        Assert.Equal(Math.BitDecrement(1.0), down.ToDouble());
    }

    [Fact]
    public void IntegerRounding_MatchesDouble()
    {
        var rng = new Random(5);
        for (int i = 0; i < 1000; i++)
        {
            double d = (rng.NextDouble() - 0.5) * 1000;
            BigFloat r = BigFloat.FromDouble(d, 53);
            Assert.Equal(Math.Floor(d), r.Floor().ToDouble());
            Assert.Equal(Math.Ceiling(d), r.Ceiling().ToDouble());
            Assert.Equal(Math.Truncate(d), r.Truncate().ToDouble());
            Assert.Equal(Math.Round(d, MidpointRounding.AwayFromZero), BigFloat.Round(r, MidpointRounding.AwayFromZero).ToDouble());
            Assert.Equal(Math.Round(d), BigFloat.Round(r).ToDouble());
        }
        Assert.Equal(2.0, BigFloat.Round(BigFloat.FromDouble(2.5, 53)).ToDouble());
        Assert.Equal(3.0, BigFloat.Round(BigFloat.FromDouble(2.5, 53), MidpointRounding.AwayFromZero).ToDouble());
    }

    [Fact]
    public void FmodRemainder_MatchDouble()
    {
        var rng = new Random(6);
        for (int i = 0; i < 1000; i++)
        {
            double a = (rng.NextDouble() - 0.5) * 1e6;
            double b = (rng.NextDouble() - 0.5) * 100;
            if (b == 0)
            {
                continue;
            }
            BigFloat ra = BigFloat.FromDouble(a, 53), rb = BigFloat.FromDouble(b, 53);
            Assert.Equal(Math.IEEERemainder(a, b), BigFloat.Ieee754Remainder(ra, rb, 53).ToDouble());
            Assert.Equal(a % b, BigFloat.TruncatedRemainder(ra, rb, 53).ToDouble());
        }
    }
}
