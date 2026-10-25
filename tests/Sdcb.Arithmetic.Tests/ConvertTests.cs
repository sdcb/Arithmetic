using System.Numerics;
using Xunit;

namespace Sdcb.Arithmetic.Tests;

public class ConvertTests
{
    [Fact]
    public void ToString_Base10_MatchesBigInteger()
    {
        var rng = new Random(60);
        for (int i = 0; i < 300; i++)
        {
            var (a, ba) = TestHelpers.RandomPair(rng, 20);
            Assert.Equal(ba.ToString(), a.ToString());
        }
    }

    [Fact]
    public void ToString_LargeNumber_DivideAndConquer()
    {
        var rng = new Random(61);
        var (a, ba) = TestHelpers.RandomPair(rng, 800, allowNegative: false);
        Assert.Equal(ba.ToString(), a.ToString());
    }

    [Fact]
    public void ParseAndFormat_RoundTripsAllBases()
    {
        var rng = new Random(62);
        for (int b = 2; b <= 62; b++)
        {
            for (int i = 0; i < 20; i++)
            {
                var (a, _) = TestHelpers.RandomPair(rng, 8);
                string s = a.ToString(b);
                BigInt back = BigInt.Parse(s, b);
                Assert.Equal(a, back);
            }
        }
    }

    [Fact]
    public void Parse_MatchesBigInteger()
    {
        var rng = new Random(63);
        for (int i = 0; i < 200; i++)
        {
            var (_, ba) = TestHelpers.RandomPair(rng, 15);
            string s = ba.ToString();
            TestHelpers.AssertEqual(ba, BigInt.Parse(s), $"parse #{i}: {s}");
        }
    }

    [Fact]
    public void Parse_LargeDecimal_DivideAndConquer()
    {
        var rng = new Random(64);
        var (_, ba) = TestHelpers.RandomPair(rng, 700, allowNegative: false);
        string s = ba.ToString();
        TestHelpers.AssertEqual(ba, BigInt.Parse(s), "large parse");
    }

    [Fact]
    public void Parse_AutoBasePrefixes()
    {
        Assert.Equal(new BigInt(255), BigInt.Parse("0xff", 0));
        Assert.Equal(new BigInt(255), BigInt.Parse("0XFF", 0));
        Assert.Equal(new BigInt(5), BigInt.Parse("0b101", 0));
        Assert.Equal(new BigInt(8), BigInt.Parse("0o10", 0));
        Assert.Equal(new BigInt(8), BigInt.Parse("010", 0));
        Assert.Equal(new BigInt(-42), BigInt.Parse("-42", 0));
        Assert.Equal(new BigInt(1_000_000), BigInt.Parse("1_000_000"));
    }

    [Fact]
    public void Parse_Invalid_ReturnsFalse()
    {
        Assert.False(BigInt.TryParse("", out _));
        Assert.False(BigInt.TryParse("abc", 10, out _));
        Assert.False(BigInt.TryParse("-", out _));
        Assert.False(BigInt.TryParse("12z", 16, out _));
        Assert.True(BigInt.TryParse("-0", out BigInt zero));
        Assert.True(zero.IsZero);
    }

    [Fact]
    public void HexFormat_Works()
    {
        var v = BigInt.Parse("deadbeefcafebabe0123456789", 16);
        Assert.Equal("deadbeefcafebabe0123456789", v.ToString(16));
        Assert.Equal("DEADBEEFCAFEBABE0123456789", v.ToString("X", null));
        Assert.Equal("-ff", (-new BigInt(255)).ToString(16));
    }

    [Fact]
    public void PowerOfTwoBases_MatchGeneralPath()
    {
        var rng = new Random(65);
        for (int i = 0; i < 50; i++)
        {
            var (a, ba) = TestHelpers.RandomPair(rng, 12, allowNegative: false);
            if (a.IsZero)
            {
                continue;
            }
            // Cross-check base 16 against BigInteger's hex formatting.
            string expected = ba.ToString("x").TrimStart('0');
            if (expected.Length == 0)
            {
                expected = "0";
            }
            Assert.Equal(expected, a.ToString(16));
        }
    }

    [Fact]
    public void MagnitudeBytes_RoundTrip()
    {
        var rng = new Random(66);
        for (int i = 0; i < 200; i++)
        {
            var (a, _) = TestHelpers.RandomPair(rng, 6);
            Assert.Equal(a, BigInt.FromMagnitudeBytes(a.ToMagnitudeBytes(), a.Sign));
            Assert.Equal(a, BigInt.FromMagnitudeBytes(a.ToMagnitudeBytes(bigEndian: true), a.Sign, bigEndian: true));
        }
    }
}
