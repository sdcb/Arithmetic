using Xunit;

namespace Sdcb.Arithmetic.Tests;

/// <summary>ToDouble(RoundingMode) on DoubleDouble/QuadDouble and the fused DoubleDouble.FusedMultiplyAdd.</summary>
public class FixedFloatApiTests
{
    private static readonly RoundingMode[] Modes =
        [RoundingMode.NearestEven, RoundingMode.TowardZero, RoundingMode.TowardPositive, RoundingMode.TowardNegative, RoundingMode.AwayFromZero];

    private static double RandomDouble(Random rng, int minExp, int maxExp) =>
        (rng.Next(2) == 0 ? 1 : -1) * (1 + rng.NextDouble()) * Math.ScaleB(1.0, rng.Next(minExp, maxExp));

    [Fact]
    public void DoubleDoubleToDoubleMatchesExactRounding()
    {
        var rng = new Random(21);
        for (int i = 0; i < 4000; i++)
        {
            double hi = RandomDouble(rng, -1000, 1000);
            double lo = i % 7 == 0 ? 0.0 : hi * Math.ScaleB(rng.NextDouble() - 0.5, -53);
            DoubleDouble x = new(hi, lo);
            foreach (RoundingMode m in Modes)
            {
                Assert.Equal(((BigFloat)x).ToDouble(m), x.ToDouble(m));
            }
            Assert.Equal((double)x, x.ToDouble());
        }
    }

    [Fact]
    public void QuadDoubleToDoubleMatchesExactRounding()
    {
        var rng = new Random(22);
        for (int i = 0; i < 4000; i++)
        {
            double x0 = RandomDouble(rng, -800, 800);
            QuadDouble q = i % 7 == 0 ? new QuadDouble(x0) : QuadDouble.Sqrt(new QuadDouble(Math.Abs(x0))) * Math.Sign(x0);
            foreach (RoundingMode m in Modes)
            {
                Assert.Equal(((BigFloat)q).ToDouble(m), q.ToDouble(m));
            }
            Assert.Equal((double)q, q.ToDouble());
        }
    }

    [Fact]
    public void DoubleDoubleToDoubleSpecialValues()
    {
        Assert.True(double.IsNaN(DoubleDouble.NaN.ToDouble(RoundingMode.TowardZero)));
        Assert.Equal(double.PositiveInfinity, DoubleDouble.PositiveInfinity.ToDouble(RoundingMode.TowardNegative));
        Assert.True(double.IsNegative(((DoubleDouble)(-0.0)).ToDouble(RoundingMode.TowardPositive)));
        Assert.True(double.IsNegative(((QuadDouble)(-0.0)).ToDouble(RoundingMode.TowardPositive)));
    }

    private static double Log2RelErr(DoubleDouble got, BigFloat exact)
    {
        if (exact.IsZero)
        {
            return ((BigFloat)got).IsZero ? double.NegativeInfinity : 0;
        }
        BigFloat err = BigFloat.Abs(BigFloat.Divide(BigFloat.Subtract(got, exact, 600), exact, 64));
        return err.IsZero ? double.NegativeInfinity : BigFloat.Log2(err, 64).ToDouble();
    }

    [Fact]
    public void FusedMultiplyAddIsAccurateUnderCancellation()
    {
        var rng = new Random(23);
        double worst = double.NegativeInfinity;
        for (int i = 0; i < 4000; i++)
        {
            DoubleDouble a = DoubleDouble.Sqrt(new DoubleDouble(1 + rng.NextDouble() * 1e6, 0)) * Math.ScaleB(1.0, rng.Next(-300, 300));
            DoubleDouble b = DoubleDouble.Sqrt(new DoubleDouble(1 + rng.NextDouble() * 1e6, 0)) * Math.ScaleB(1.0, rng.Next(-300, 300));
            // Cases: residual (total cancellation of the rounded product), partial cancellation, random addend.
            DoubleDouble c = (i % 3) switch
            {
                0 => -(a * b),
                1 => -(a * b) * (1 + Math.ScaleB(rng.NextDouble(), -40)),
                _ => new DoubleDouble(RandomDouble(rng, -10, 10), 0) * (a * b),
            };
            BigFloat exact = BigFloat.FusedMultiplyAdd(a, b, c, 600);
            double e = Log2RelErr(DoubleDouble.FusedMultiplyAdd(a, b, c), exact);
            worst = Math.Max(worst, e);
            Assert.True(e < -103, $"FMA error 2^{e} at ({a}, {b}, {c})");
        }
        Assert.True(worst < -103);
    }

    [Fact]
    public void MultiplyAddEstimateIsTheUnfusedForm()
    {
        DoubleDouble a = DoubleDouble.Sqrt(2), b = DoubleDouble.Sqrt(3);
        DoubleDouble c = -(a * b);
        Assert.Equal(DoubleDouble.Zero, DoubleDouble.MultiplyAddEstimate(a, b, c));
        Assert.NotEqual(DoubleDouble.Zero, DoubleDouble.FusedMultiplyAdd(a, b, c));
    }

    [Fact]
    public void FusedMultiplyAddSpecialValues()
    {
        Assert.True(DoubleDouble.IsNaN(DoubleDouble.FusedMultiplyAdd(DoubleDouble.PositiveInfinity, DoubleDouble.Zero, DoubleDouble.One)));
        Assert.Equal(DoubleDouble.PositiveInfinity, DoubleDouble.FusedMultiplyAdd(DoubleDouble.MaxValue, 2.0, DoubleDouble.One));
        Assert.Equal(new DoubleDouble(3, 0), DoubleDouble.FusedMultiplyAdd(DoubleDouble.Zero, 5.0, 3.0));
        Assert.Equal(new DoubleDouble(7, 0), DoubleDouble.FusedMultiplyAdd(2.0, 3.0, DoubleDouble.One));
    }
}
