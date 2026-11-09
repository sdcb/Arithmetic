using Xunit;

namespace Sdcb.Arithmetic.Tests;

public class FloatInteropTests
{
    private const int P = 400;

    [Fact]
    public void DoubleDoubleToQuadDoubleIsExact()
    {
        var rng = new Random(11);
        for (int i = 0; i < 2000; i++)
        {
            DoubleDouble dd = new(rng.NextDouble() * Math.ScaleB(1.0, rng.Next(-900, 900)), rng.NextDouble() * Math.ScaleB(1.0, rng.Next(-1000, -60)));
            QuadDouble qd = dd;
            Assert.Equal(0, ((Real)dd).CompareTo((Real)qd));
            Assert.Equal(dd, (DoubleDouble)qd);
        }
    }

    [Fact]
    public void QuadDoubleToDoubleDoubleIsNearlyRounded()
    {
        var rng = new Random(12);
        for (int i = 0; i < 2000; i++)
        {
            QuadDouble q = QuadDouble.Sqrt(new QuadDouble(1.0 + rng.NextDouble() * 1e6));
            DoubleDouble d = (DoubleDouble)q;
            Real err = Real.Abs(Real.Subtract(d, q, P) / (Real)q);
            Assert.True(err < Real.FromBigInt(BigInt.One, -105, 2), $"{q}");
        }
    }

    [Fact]
    public void DoubleDoubleToFloat128IsCorrectlyRounded()
    {
        var rng = new Random(13);
        for (int i = 0; i < 2000; i++)
        {
            DoubleDouble dd = new(rng.NextDouble() * Math.ScaleB(1.0, rng.Next(-1000, 1000)), rng.NextDouble() * Math.ScaleB(1.0, rng.Next(-1074, -900)));
            Float128 f = (Float128)dd;
            Assert.Equal(Float128.FromReal(dd).ToBits(), f.ToBits());
        }
    }

    [Fact]
    public void Float128ToQuadDoubleIsExactAndBack()
    {
        var rng = new Random(14);
        for (int i = 0; i < 2000; i++)
        {
            Float128 f = Float128.Sqrt((Float128)(1.0 + rng.NextDouble())) * Math.ScaleB(1.0, rng.Next(-900, 900));
            QuadDouble q = (QuadDouble)f;
            Assert.Equal(0, ((Real)f).CompareTo((Real)q));
            Assert.Equal(f.ToBits(), ((Float128)q).ToBits());
        }
    }

    [Fact]
    public void Float128ToDoubleDoubleIsNearlyRounded()
    {
        Float128 f = Float128.Sqrt(2);
        DoubleDouble d = (DoubleDouble)f;
        Assert.Equal((double)f, d.Hi);
        Real err = Real.Abs(Real.Subtract(d, f, P));
        Assert.True(err < Real.FromBigInt(BigInt.One, -105, 2));
    }

    [Fact]
    public void SpecialValues()
    {
        Assert.True(double.IsPositiveInfinity(((QuadDouble)Float128.MaxValue).X0));
        Assert.True(QuadDouble.IsNaN((QuadDouble)Float128.NaN));
        Assert.True(Float128.IsNaN((Float128)DoubleDouble.NaN));
        Assert.True(Float128.IsPositiveInfinity((Float128)DoubleDouble.PositiveInfinity));
        Assert.True(DoubleDouble.IsNegative((DoubleDouble)(Float128)(-0.0)));
        Assert.Equal(0.0, ((DoubleDouble)Float128.Epsilon).Hi);
        Assert.True(double.IsNegative(((DoubleDouble)(QuadDouble)(-0.0)).Hi));
        Assert.True(Float128.IsNegative((Float128)(DoubleDouble)(-0.0)));
        Assert.True(Float128.IsNegative((Float128)(QuadDouble)(-0.0)));
    }
}
