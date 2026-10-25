using Xunit;

namespace Sdcb.Arithmetic.Tests;

/// <summary>
/// At 53-bit precision every correctly-rounded Real operation must agree
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
            Real ra = Real.FromDouble(a, 53), rb = Real.FromDouble(b, 53);
            Assert.Equal(a + b, Real.Add(ra, rb, 53).ToDouble());
            Assert.Equal(a - b, Real.Subtract(ra, rb, 53).ToDouble());
            Assert.Equal(a * b, Real.Multiply(ra, rb, 53).ToDouble());
            if (b != 0)
            {
                Assert.Equal(a / b, Real.Divide(ra, rb, 53).ToDouble());
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
            Assert.Equal(Math.Sqrt(a), Real.Sqrt(Real.FromDouble(a, 53), 53).ToDouble());
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
                Real.FusedMultiplyAdd(Real.FromDouble(a, 53), Real.FromDouble(b, 53), Real.FromDouble(c, 53), 53).ToDouble());
        }
    }

    [Fact]
    public void DoubleRoundTrip_Exact()
    {
        var rng = new Random(4);
        for (int i = 0; i < 3000; i++)
        {
            double d = RandomDouble(rng);
            Assert.Equal(d, Real.FromDouble(d).ToDouble());
        }
        Assert.True(double.IsNaN(Real.NaN.ToDouble()));
        Assert.Equal(double.PositiveInfinity, Real.PositiveInfinity.ToDouble());
        Assert.Equal(double.NegativeInfinity, Real.NegativeInfinity.ToDouble());
        Assert.Equal(-0.0, Real.NegativeZero.ToDouble());
        Assert.True(double.IsNegative(Real.NegativeZero.ToDouble()));
    }

    [Fact]
    public void RoundingModes_DirectedDivision()
    {
        // 1/3 at 8 bits: 0.0101010101...b → check all modes.
        Real oneThird(RealRounding r) => Real.Divide(Real.One, Real.FromBigInt(3), 8, r);
        Real down = oneThird(RealRounding.TowardZero);
        Real up = oneThird(RealRounding.TowardPositive);
        Assert.True(down.CompareTo(up) < 0);
        Assert.Equal(down, oneThird(RealRounding.TowardNegative));
        Assert.Equal(up, oneThird(RealRounding.AwayFromZero));
        // Negative: directions flip.
        Real m = Real.Divide(Real.MinusOne, Real.FromBigInt(3), 8, RealRounding.TowardPositive);
        Assert.Equal((-down).ToDouble(), m.ToDouble());
    }

    [Fact]
    public void SpecialValues_Arithmetic()
    {
        Assert.True((Real.PositiveInfinity + Real.NegativeInfinity).IsNaN);
        Assert.True((Real.Zero * Real.PositiveInfinity).IsNaN);
        Assert.True((Real.PositiveInfinity / Real.PositiveInfinity).IsNaN);
        Assert.True((Real.One / Real.Zero).IsPositiveInfinity);
        Assert.True((Real.MinusOne / Real.Zero).IsNegativeInfinity);
        Assert.True((Real.One / Real.NegativeZero).IsNegativeInfinity);
        Assert.True(Real.Sqrt(Real.MinusOne).IsNaN);
        Assert.True(Real.Sqrt(Real.NegativeZero).IsZero);
        Assert.True(Real.Sqrt(Real.NegativeZero).SignBit);
        // 0 + -0 = +0 except toward negative.
        Assert.False((Real.Zero + Real.NegativeZero).SignBit);
        Assert.True(Real.Add(Real.Zero, Real.NegativeZero, 53, RealRounding.TowardNegative).SignBit);
        // NaN comparisons.
        Assert.False(Real.NaN == Real.NaN);
        Assert.True(Real.NaN != Real.NaN);
        Assert.False(Real.NaN < Real.One);
        Assert.False(Real.NaN >= Real.One);
        Assert.True(Real.NaN.Equals(Real.NaN)); // bitwise Equals
    }

    [Fact]
    public void Comparison_TotalOrder()
    {
        Assert.True(Real.NegativeZero.CompareTo(Real.Zero) < 0);
        Assert.True(Real.NaN.CompareTo(Real.NegativeInfinity) < 0);
        Assert.True(Real.NegativeInfinity < Real.MinusOne);
        Assert.True(Real.PositiveInfinity > Real.FromBigInt(BigInt.Parse("1" + new string('0', 100))));
        Assert.True(Real.NegativeZero == Real.Zero); // IEEE equality
    }

    [Fact]
    public void NextAboveBelow()
    {
        Real one = Real.FromDouble(1.0, 53);
        Real up = Real.BitIncrement(one);
        Assert.Equal(Math.BitIncrement(1.0), up.ToDouble());
        Assert.Equal(1.0, Real.BitDecrement(up).ToDouble());
        Real down = Real.BitDecrement(one);
        Assert.Equal(Math.BitDecrement(1.0), down.ToDouble());
    }

    [Fact]
    public void IntegerRounding_MatchesDouble()
    {
        var rng = new Random(5);
        for (int i = 0; i < 1000; i++)
        {
            double d = (rng.NextDouble() - 0.5) * 1000;
            Real r = Real.FromDouble(d, 53);
            Assert.Equal(Math.Floor(d), r.Floor().ToDouble());
            Assert.Equal(Math.Ceiling(d), r.Ceiling().ToDouble());
            Assert.Equal(Math.Truncate(d), r.Truncate().ToDouble());
            Assert.Equal(Math.Round(d, MidpointRounding.AwayFromZero), Real.Round(r, MidpointRounding.AwayFromZero).ToDouble());
            Assert.Equal(Math.Round(d), Real.Round(r).ToDouble());
        }
        Assert.Equal(2.0, Real.Round(Real.FromDouble(2.5, 53)).ToDouble());
        Assert.Equal(3.0, Real.Round(Real.FromDouble(2.5, 53), MidpointRounding.AwayFromZero).ToDouble());
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
            Real ra = Real.FromDouble(a, 53), rb = Real.FromDouble(b, 53);
            Assert.Equal(Math.IEEERemainder(a, b), Real.Ieee754Remainder(ra, rb, 53).ToDouble());
            Assert.Equal(a % b, Real.TruncatedRemainder(ra, rb, 53).ToDouble());
        }
    }
}
