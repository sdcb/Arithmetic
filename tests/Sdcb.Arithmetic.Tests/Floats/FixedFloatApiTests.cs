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
    public void ToSingleToHalfAndBigFloatRoundTrip()
    {
        var rng = new Random(26);
        for (int i = 0; i < 3000; i++)
        {
            DoubleDouble dd = new(RandomDouble(rng, -140, 120), RandomDouble(rng, -200, -150));
            QuadDouble qd = QuadDouble.Sqrt(new QuadDouble(Math.Abs(dd.Hi))) * Math.Sign(dd.Hi);
            BigFloat edd = dd, eqd = qd;
            Assert.Equal(0, dd.ToBigFloat().CompareTo(edd));
            Assert.Equal(0, qd.ToBigFloat().CompareTo(eqd));
            Assert.Equal(dd, DoubleDouble.FromBigFloat(dd.ToBigFloat()));
            Assert.Equal(qd, QuadDouble.FromBigFloat(qd.ToBigFloat()));
            foreach (RoundingMode m in Modes)
            {
                Assert.Equal(edd.ToSingle(m), dd.ToSingle(m));
                Assert.Equal(eqd.ToSingle(m), qd.ToSingle(m));
                Assert.Equal(edd.ToHalf(m), dd.ToHalf(m));
                Assert.Equal(eqd.ToHalf(m), qd.ToHalf(m));
            }
        }
    }

    [Fact]
    public void QuadDoubleMultiplyAddEstimateIsUnfused()
    {
        QuadDouble a = QuadDouble.Sqrt(2), b = QuadDouble.Sqrt(3);
        Assert.Equal(a * b + QuadDouble.One, QuadDouble.MultiplyAddEstimate(a, b, QuadDouble.One));
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
            Assert.True(e < -105, $"FMA error 2^{e} at ({a}, {b}, {c})");
        }
        Assert.True(worst < -105);
    }

    private static DoubleDouble RandomDd(Random rng, int minExp, int maxExp)
    {
        double hi = RandomDouble(rng, minExp, maxExp);
        // Low words of every size, including gaps far below ulp(hi) and zero.
        double lo = rng.Next(5) == 0 ? 0.0 : hi * Math.ScaleB(rng.NextDouble() - 0.5, -53 - rng.Next(0, 60));
        return new DoubleDouble(hi, lo);
    }

    [Fact]
    public void FusedMultiplyAddMatchesExactAcrossInputClasses()
    {
        var rng = new Random(24);
        for (int i = 0; i < 60000; i++)
        {
            DoubleDouble a = RandomDd(rng, -300, 300), b = RandomDd(rng, -300, 300);
            DoubleDouble c = (i % 6) switch
            {
                0 => RandomDd(rng, -600, 600),                                            // unrelated magnitude
                1 => -(a * b),                                                            // residual of the rounded product
                2 => -(DoubleDouble)((QuadDouble)a * b),                                  // deeper cancellation
                3 => -(a * b) + (a * b) * Math.ScaleB(rng.NextDouble() - 0.5, -rng.Next(1, 100)), // partial
                4 => (a * b) * RandomDouble(rng, -2, 2),                                  // same order, random sign
                _ => (a * b) * Math.ScaleB(1.0, rng.Next(-200, 200)),                     // large gaps
            };
            BigFloat exact = BigFloat.FusedMultiplyAdd(a, b, c, 800);
            DoubleDouble got = DoubleDouble.FusedMultiplyAdd(a, b, c);
            double e = Log2RelErr(got, exact);
            Assert.True(e < -105, $"case {i % 6}: FMA error 2^{e} at ({a.Hi}, {a.Lo}) ({b.Hi}, {b.Lo}) ({c.Hi}, {c.Lo})");
        }
    }

    [Fact]
    public void FusedMultiplyAddEdgeRanges()
    {
        // Exact zero result is +0.
        DoubleDouble a = new(3.0, 0), b = new(5.0, 0);
        DoubleDouble r = DoubleDouble.FusedMultiplyAdd(a, b, new DoubleDouble(-15.0, 0));
        Assert.Equal(0.0, r.Hi);
        Assert.False(double.IsNegative(r.Hi));

        var rng = new Random(25);
        foreach (int scale in new[] { -1000, -900, -850, -800, 950, 1000, 1010 })
        {
            for (int i = 0; i < 300; i++)
            {
                DoubleDouble x = RandomDd(rng, scale / 2 - 2, scale / 2 + 2), y = RandomDd(rng, scale / 2 - 2, scale / 2 + 2);
                DoubleDouble z = i % 2 == 0 ? -(x * y) : RandomDd(rng, scale - 3, scale + 3);
                if (!double.IsFinite((x * y).Hi) || !double.IsFinite(z.Hi))
                {
                    continue;
                }
                BigFloat exact = BigFloat.FusedMultiplyAdd(x, y, z, 2400);
                DoubleDouble got = DoubleDouble.FusedMultiplyAdd(x, y, z);
                BigFloat expected = DoubleDouble.FromBigFloat(exact);
                // Near the subnormal range double-double itself only has ~53 bits; compare to the
                // nearest double-double there, elsewhere require the full bound.
                if (Math.Abs(((double)expected)) < Math.ScaleB(1.0, -969))
                {
                    Assert.Equal(expected.ToDouble(), got.Hi);
                }
                else
                {
                    Assert.True(Log2RelErr(got, exact) < -105, $"scale {scale}: ({x}, {y}, {z})");
                }
            }
        }
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
