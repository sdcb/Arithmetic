using System.Globalization;
using System.Numerics;
using System.Text;
using Xunit;

namespace Sdcb.Arithmetic.Tests.Generic;

public class BigRationalGenericTests
{
    private static readonly MidpointRounding[] Modes =
    [
        MidpointRounding.ToEven, MidpointRounding.AwayFromZero, MidpointRounding.ToZero,
        MidpointRounding.ToNegativeInfinity, MidpointRounding.ToPositiveInfinity,
    ];

    /// <summary>Random rationals with terminating decimal expansions, mirrored exactly as decimal.</summary>
    private static IEnumerable<(BigRational Q, decimal D)> DecimalSamples(int seed, int count = 400)
    {
        var rng = new Random(seed);
        int[] dens = [1, 2, 4, 5, 8, 10, 16, 20, 25, 40, 100, 125, 1000, 1024];
        foreach (decimal d in new[] { 0m, 0.5m, -0.5m, 1.5m, -1.5m, 2.5m, -2.5m, 0.125m, -7.25m })
        {
            yield return (d, d);
        }
        for (int i = 0; i < count; i++)
        {
            long num = rng.NextInt64(-1_000_000, 1_000_000);
            int den = dens[rng.Next(dens.Length)];
            yield return (new BigRational(num, den), (decimal)num / den);
        }
    }

    [Fact]
    public void Rounding_MatchesDecimal()
    {
        foreach ((BigRational q, decimal d) in DecimalSamples(20))
        {
            Assert.Equal((BigRational)decimal.Floor(d), BigRational.Floor(q));
            Assert.Equal((BigRational)decimal.Ceiling(d), BigRational.Ceiling(q));
            Assert.Equal((BigRational)decimal.Truncate(d), BigRational.Truncate(q));
            Assert.Equal((BigRational)decimal.Round(d), BigRational.Round(q));
            foreach (MidpointRounding m in Modes)
            {
                Assert.Equal((BigRational)decimal.Round(d, m), BigRational.Round(q, m));
                Assert.Equal((BigRational)decimal.Round(d, 1, m), BigRational.Round(q, 1, m));
                Assert.Equal((BigRational)decimal.Round(d, 2, m), BigRational.Round(q, 2, m));
                Assert.Equal((BigRational)decimal.Round(d, m), Generic(q, m));
            }
        }

        static T Generic<T>(T v, MidpointRounding m) where T : INumber<T>
            => T.CreateChecked(BigRational.Round(BigRational.CreateChecked(v), m));
    }

    [Fact]
    public void Remainder_MatchesDecimal()
    {
        var samples = DecimalSamples(21, 200).ToArray();
        for (int i = 0; i + 1 < samples.Length; i++)
        {
            (BigRational a, decimal da) = samples[i];
            (BigRational b, decimal db) = samples[i + 1];
            if (db == 0)
            {
                Assert.Throws<DivideByZeroException>(() => a % b);
                continue;
            }
            Assert.Equal((BigRational)(da % db), a % b);
        }
        BigRational x = new(7, 3);
        Assert.Equal(new BigRational(10, 3), ++x);
        Assert.Equal(new BigRational(7, 3), --x);
    }

    [Fact]
    public void DecimalConversions()
    {
        foreach ((BigRational q, decimal d) in DecimalSamples(22))
        {
            Assert.Equal(d, (decimal)q);
            Assert.Equal(d, decimal.CreateChecked(q));
            Assert.Equal(q, BigRational.CreateChecked(d));
        }
        Assert.Equal(1m / 3m, (decimal)new BigRational(1, 3));
        Assert.Equal(2m / 3m, (decimal)new BigRational(2, 3));
        Assert.Equal(-1m / 7m, (decimal)new BigRational(-1, 7));
        Assert.Equal(decimal.MaxValue, (decimal)(BigRational)decimal.MaxValue);
        Assert.Equal(decimal.MinValue, (decimal)(BigRational)decimal.MinValue);
        BigRational tooBig = (BigRational)decimal.MaxValue + 1;
        Assert.Throws<OverflowException>(() => (decimal)tooBig);
        Assert.Equal(decimal.MaxValue, decimal.CreateSaturating(tooBig));
        Assert.Equal(0m, (decimal)new BigRational(BigInt.One, BigInt.One << 200));
    }

    /// <summary>d is the nearest T to q (ties to even) iff neither neighbour is strictly closer and ties pick the even significand.</summary>
    private static void AssertNearest<T>(BigRational q, T d) where T : IBinaryFloatingPointIeee754<T>, IMinMaxValue<T>
    {
        Assert.True(T.IsFinite(d) || BigRational.Abs(q) > new BigRational(double.CreateChecked(T.MaxValue)), $"{q} -> {d}");
        if (!T.IsFinite(d))
        {
            return;
        }
        BigRational err = BigRational.Abs(q - BigRational.CreateChecked(d));
        foreach (T n in new[] { T.BitIncrement(d), T.BitDecrement(d) })
        {
            if (!T.IsFinite(n))
            {
                continue;
            }
            BigRational e2 = BigRational.Abs(q - BigRational.CreateChecked(n));
            Assert.True(err <= e2, $"{q}: {d} not nearest (neighbour {n})");
        }
    }

    [Fact]
    public void BinaryFloatConversions_AreCorrectlyRounded()
    {
        var rng = new Random(23);
        for (int i = 0; i < 400; i++)
        {
            BigInt num = TestHelpers.RandomPair(rng, 2).A;
            BigInt den = TestHelpers.RandomPair(rng, 2, allowNegative: false).A;
            if (den.IsZero)
            {
                continue;
            }
            int shift = rng.Next(-1200, 1200);
            BigRational q = new BigRational(num, den).ShiftLeft(shift);
            AssertNearest(q, (double)q);
            AssertNearest(q, double.CreateChecked(q));
            BigRational small = new BigRational(num, den).ShiftLeft(rng.Next(-170, 140));
            AssertNearest(small, (float)small);
            AssertNearest(small, float.CreateChecked(small));
            BigRational tiny = new BigRational(num, den).ShiftLeft(rng.Next(-40, 20));
            AssertNearest(tiny, (Half)tiny);
            AssertNearest(tiny, Half.CreateChecked(tiny));
        }
        // Subnormal ties and overflow.
        Assert.Equal(double.Epsilon, (double)new BigRational(BigInt.One, BigInt.One << 1074));
        Assert.Equal(0.0, (double)new BigRational(BigInt.One, BigInt.One << 1075)); // tie -> even (zero)
        Assert.Equal(double.Epsilon * 2, (double)new BigRational(new BigInt(3), BigInt.One << 1075)); // 1.5 ulp -> 2
        Assert.Equal(double.PositiveInfinity, (double)new BigRational(BigInt.One << 1024));
        Assert.Equal(float.Epsilon, (float)new BigRational(BigInt.One, BigInt.One << 149));
    }

    [Fact]
    public void IntegerConversions_TruncateTowardZero()
    {
        BigRational q = new(-7, 2);
        Assert.Equal(-3, (int)q);
        Assert.Equal(-3, int.CreateChecked(q));
        Assert.Equal((byte)0, byte.CreateSaturating(q));
        Assert.Equal(unchecked((byte)-3), byte.CreateTruncating(q));
        Assert.Throws<OverflowException>(() => (uint)q);
        Assert.Equal(new BigInteger(-3), (BigInteger)q);
        Assert.Equal(new BigInt(-3), BigInt.CreateChecked(q));
        Assert.Equal(new BigRational(5), BigRational.CreateChecked((byte)5));
        Assert.Equal(new BigRational(-1, 4), BigRational.CreateChecked(-0.25f));
        Assert.Equal(new BigRational(1, 1024), (BigRational)(Half)(1.0 / 1024));
        Assert.Equal(BigRational.Zero, BigRational.CreateSaturating(double.NaN));
        Assert.Throws<OverflowException>(() => BigRational.CreateChecked(double.NaN));
        BigRational fromBigInteger = BigInteger.Pow(10, 40);
        Assert.Equal(new BigRational(BigInt.Pow(10, 40)), fromBigInteger);
    }

    [Fact]
    public void PredicatesAndMagnitude()
    {
        Assert.True(BigRational.IsEvenInteger(4));
        Assert.False(BigRational.IsEvenInteger(new BigRational(4, 3)));
        Assert.True(BigRational.IsOddInteger(-3));
        Assert.True(BigRational.IsNegative(new BigRational(-1, 2)));
        Assert.True(BigRational.IsPositive(BigRational.Zero));
        Assert.Equal(new BigRational(-3, 2), BigRational.MaxMagnitude(new BigRational(-3, 2), new BigRational(1, 1)));
        Assert.Equal(new BigRational(3, 2), BigRational.MaxMagnitude(new BigRational(-3, 2), new BigRational(3, 2)));
        Assert.Equal(new BigRational(-3, 2), BigRational.MinMagnitude(new BigRational(-3, 2), new BigRational(3, 2)));
        Assert.Equal(new BigRational(-5, 2), BigRational.CopySign(new BigRational(5, 2), -1));
        Assert.Equal(BigRational.One, BigRational.Clamp(7, -1, 1));
        Assert.Equal(-1, Sign(new BigRational(-1, 9)));

        static int Sign<T>(T v) where T : INumber<T> => T.Sign(v);
    }

    [Theory]
    [InlineData("1,234.5", NumberStyles.Number, "2469/2")]
    [InlineData("(3/4)", NumberStyles.Integer, null)]
    [InlineData("-3/4", NumberStyles.Integer, "-3/4")]
    [InlineData(" 1.25e-2 ", NumberStyles.Float, "1/80")]
    [InlineData("FF", NumberStyles.HexNumber, "-1")]
    [InlineData("1/0", NumberStyles.Integer, null)]
    [InlineData("1.5", NumberStyles.Integer, null)]
    public void ParseStyles(string s, NumberStyles style, string? expected)
    {
        bool ok = BigRational.TryParse(s, style, CultureInfo.InvariantCulture, out BigRational r);
        Assert.Equal(expected is not null, ok);
        if (ok)
        {
            Assert.Equal(BigRational.Parse(expected!), r);
            Assert.Equal(r, BigRational.Parse(Encoding.UTF8.GetBytes(s), style, CultureInfo.InvariantCulture));
        }
    }

    [Fact]
    public void Utf8()
    {
        BigRational q = new(-22, 7);
        byte[] buf = new byte[16];
        Assert.True(q.TryFormat(buf, out int n, default, null));
        Assert.Equal("-22/7", Encoding.UTF8.GetString(buf, 0, n));
        Assert.Equal(q, BigRational.Parse("-22/7"u8, null));
    }
}
