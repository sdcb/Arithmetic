using Xunit;

namespace Sdcb.Arithmetic.Tests;

public class RationalTests
{
    [Fact]
    public void Canonicalization()
    {
        var r = new BigRational(6, -4);
        Assert.Equal(new BigInt(-3), r.Numerator);
        Assert.Equal(new BigInt(2), r.Denominator);
        Assert.Equal(BigRational.Zero, new BigRational(0, 7));
        Assert.Equal(BigInt.One, BigRational.Zero.Denominator);
        Assert.Throws<DivideByZeroException>(() => new BigRational(1, 0));
    }

    [Fact]
    public void Arithmetic_MatchesRationalIdentities()
    {
        var rng = new Random(90);
        for (int i = 0; i < 300; i++)
        {
            var a = new BigRational(rng.NextInt64(-100000, 100000), rng.NextInt64(1, 100000));
            var b = new BigRational(rng.NextInt64(-100000, 100000), rng.NextInt64(1, 100000));
            // (a+b) - b == a
            Assert.Equal(a, a + b - b);
            if (!b.IsZero)
            {
                Assert.Equal(a, a * b / b);
            }
            // Cross-check against double arithmetic loosely.
            double da = a.ToDouble(), db = b.ToDouble();
            Assert.Equal(da + db, (a + b).ToDouble(), 6);
            Assert.Equal(da * db, (a * b).ToDouble(), 6);
        }
    }

    [Fact]
    public void CompareTo_Works()
    {
        Assert.True(new BigRational(1, 3) < new BigRational(1, 2));
        Assert.True(new BigRational(-1, 3) > new BigRational(-1, 2));
        Assert.True(new BigRational(2, 4) == new BigRational(1, 2));
        Assert.Equal(0, new BigRational(7, 7).CompareTo(BigRational.One));
    }

    [Fact]
    public void DoubleRoundTrip_Exact()
    {
        var rng = new Random(91);
        for (int i = 0; i < 500; i++)
        {
            double d = (rng.NextDouble() - 0.5) * Math.Pow(10, rng.Next(-20, 20));
            var r = new BigRational(d);
            Assert.Equal(d, r.ToDouble());
        }
        Assert.Equal(0.5, new BigRational(1, 2).ToDouble());
        Assert.Equal(1.0 / 3.0, new BigRational(1, 3).ToDouble());
    }

    [Fact]
    public void ParseAndFormat()
    {
        Assert.Equal(new BigRational(1, 2), BigRational.Parse("1/2"));
        Assert.Equal(new BigRational(-3, 4), BigRational.Parse("-3/4"));
        Assert.Equal(new BigRational(5, 4), BigRational.Parse("1.25"));
        Assert.Equal(new BigRational(125, 1), BigRational.Parse("1.25e2"));
        Assert.Equal(new BigRational(1, 80), BigRational.Parse("1.25e-2"));
        Assert.Equal("1/2", new BigRational(2, 4).ToString());
        Assert.Equal("7", new BigRational(7, 1).ToString());
        Assert.Equal("-a/b", new BigRational(-10, 11).ToString(16));
    }

    [Fact]
    public void Inverse_And_Pow()
    {
        Assert.Equal(new BigRational(3, 2), new BigRational(2, 3).Inverse());
        Assert.Equal(new BigRational(-3, 2), new BigRational(-2, 3).Inverse());
        Assert.Equal(new BigRational(8, 27), BigRational.Pow(new BigRational(2, 3), 3));
        Assert.Equal(new BigRational(27, 8), BigRational.Pow(new BigRational(2, 3), -3));
        Assert.Equal(BigRational.One, BigRational.Pow(new BigRational(5, 7), 0));
    }

    [Fact]
    public void ShiftLeft_Works()
    {
        Assert.Equal(new BigRational(4), new BigRational(1).ShiftLeft(2));
        Assert.Equal(new BigRational(1, 4), new BigRational(1).ShiftLeft(-2));
        Assert.Equal(new BigRational(3, 2), new BigRational(3, 8).ShiftLeft(2));
    }
}

public class RandomTests
{
    [Fact]
    public void Deterministic_WithSeed()
    {
        var r1 = new Random(12345);
        var r2 = new Random(12345);
        for (int i = 0; i < 100; i++)
        {
            Assert.Equal(BigInt.Random(r1, 300), BigInt.Random(r2, 300));
            Assert.Equal(Real.Random(r1, 200), Real.Random(r2, 200));
        }
    }

    [Fact]
    public void BigIntRandom_InRange()
    {
        var r = new Random(1);
        for (int i = 0; i < 200; i++)
        {
            long bits = i % 700;
            BigInt v = BigInt.Random(r, bits);
            Assert.True(v.Sign >= 0);
            Assert.True(v.GetBitLength() <= bits);
        }
    }

    [Fact]
    public void BigIntRandomBelow_InRange()
    {
        var r = new Random(2);
        BigInt bound = BigInt.Parse("123456789012345678901234567890");
        for (int i = 0; i < 500; i++)
        {
            BigInt v = BigInt.RandomBelow(r, bound);
            Assert.True(v.Sign >= 0 && v < bound);
        }
        Assert.Throws<ArgumentOutOfRangeException>(() => BigInt.RandomBelow(r, BigInt.Zero));
    }

    [Fact]
    public void RealRandom_Distributions()
    {
        var r = new Random(3);
        double sumU = 0, sumN = 0, sumE = 0;
        const int N = 2000;
        for (int i = 0; i < N; i++)
        {
            Real u = Real.Random(r, 100, RealRounding.TowardZero);
            Assert.True(u >= Real.Zero && u < Real.One);
            Assert.Equal(100, u.Precision);
            Real e = Real.RandomExponential(r, 80);
            Assert.True(e >= Real.Zero);
            sumU += u.ToDouble();
            sumN += Real.RandomNormal(r, 80).ToDouble();
            sumE += e.ToDouble();
        }
        Assert.InRange(sumU / N, 0.45, 0.55);
        Assert.InRange(sumN / N, -0.1, 0.1);
        Assert.InRange(sumE / N, 0.9, 1.1);
    }
}
