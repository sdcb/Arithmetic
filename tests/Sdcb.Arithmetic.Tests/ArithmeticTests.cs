using System.Numerics;
using Xunit;

namespace Sdcb.Arithmetic.Tests;

public class ArithmeticTests
{
    [Fact]
    public void AddSubMul_MatchesBigInteger()
    {
        var rng = new Random(42);
        for (int i = 0; i < 3000; i++)
        {
            int maxLimbs = i < 2500 ? 12 : 400; // exercise Karatsuba/Toom paths too
            var (a, ba) = TestHelpers.RandomPair(rng, maxLimbs);
            var (b, bb) = TestHelpers.RandomPair(rng, maxLimbs);
            TestHelpers.AssertEqual(ba + bb, a + b, $"add #{i}");
            TestHelpers.AssertEqual(ba - bb, a - b, $"sub #{i}");
            TestHelpers.AssertEqual(ba * bb, a * b, $"mul #{i}");
        }
    }

    [Fact]
    public void LargeBalancedMul_HitsToom3()
    {
        var rng = new Random(7);
        for (int i = 0; i < 8; i++)
        {
            var (a, ba) = TestHelpers.RandomPair(rng, 600, allowNegative: false);
            var (b, bb) = TestHelpers.RandomPair(rng, 600, allowNegative: false);
            TestHelpers.AssertEqual(ba * bb, a * b, $"toom mul #{i}");
            TestHelpers.AssertEqual(ba * ba, a.Square(), $"toom sqr #{i}");
        }
    }

    [Fact]
    public void UnbalancedMul_Chunking()
    {
        var rng = new Random(11);
        for (int i = 0; i < 20; i++)
        {
            var (a, ba) = TestHelpers.RandomPair(rng, 500, allowNegative: false);
            var (b, bb) = TestHelpers.RandomPair(rng, 40, allowNegative: false);
            TestHelpers.AssertEqual(ba * bb, a * b, $"unbalanced #{i}");
        }
    }

    [Fact]
    public void DivRem_MatchesBigInteger()
    {
        var rng = new Random(43);
        for (int i = 0; i < 2000; i++)
        {
            var (a, ba) = TestHelpers.RandomPair(rng, 30);
            var (b, bb) = TestHelpers.RandomPair(rng, 10);
            if (b.IsZero)
            {
                continue;
            }
            BigInt.DivRem(a, b, out BigInt q, out BigInt r);
            TestHelpers.AssertEqual(BigInteger.Divide(ba, bb), q, $"div #{i}");
            TestHelpers.AssertEqual(BigInteger.Remainder(ba, bb), r, $"rem #{i}");
            Assert.Equal(a, q * b + r);
        }
    }

    [Fact]
    public void DivRem_Large_BurnikelZiegler()
    {
        var rng = new Random(53);
        for (int i = 0; i < 60; i++)
        {
            // Divisors well above the BZ threshold, various dividend/divisor shapes.
            var (a, ba) = TestHelpers.RandomPair(rng, 700, allowNegative: false);
            var (b, bb) = TestHelpers.RandomPair(rng, 40 + rng.Next(300), allowNegative: false);
            if (b.IsZero || a.IsZero)
            {
                continue;
            }
            if (BigInt.CompareAbs(a, b) < 0)
            {
                (a, b) = (b, a);
                (ba, bb) = (bb, ba);
            }
            BigInt.DivRem(a, b, out BigInt q, out BigInt r);
            TestHelpers.AssertEqual(BigInteger.Divide(ba, bb), q, $"bz div #{i}");
            TestHelpers.AssertEqual(BigInteger.Remainder(ba, bb), r, $"bz rem #{i}");
        }
        // Adversarial: dividend = divisor * x + small, and all-ones patterns.
        for (int i = 0; i < 30; i++)
        {
            var (b, _) = TestHelpers.RandomPair(rng, 64 + rng.Next(200), allowNegative: false);
            var (x, _) = TestHelpers.RandomPair(rng, 64 + rng.Next(200), allowNegative: false);
            if (b.IsZero)
            {
                continue;
            }
            BigInt a = b * x + (x.IsZero ? BigInt.Zero : BigInt.Mod(x, b));
            BigInt.DivRem(a, b, out BigInt q, out BigInt r);
            Assert.Equal(a, q * b + r);
            Assert.True(BigInt.CompareAbs(r, b) < 0);
        }
        // Divisor with top limb 0xFFFF... to stress the qhat = B-1 paths.
        BigInt ones = (BigInt.One << (64 * 100)) - BigInt.One;
        BigInt big = (ones << (64 * 100)) + ones - BigInt.One;
        BigInt.DivRem(big, ones, out BigInt q2, out BigInt r2);
        Assert.Equal(big, q2 * ones + r2);
        Assert.True(BigInt.CompareAbs(r2, ones) < 0);
        // Dividend with a tiny top limb but full divisor (quotient one limb shorter
        // than an-dn+1) — regression for the BZ quotient copy bounds.
        BigInt dFull = (BigInt.One << (64 * 80)) - BigInt.One;
        BigInt aSmallTop = (BigInt.One << (64 * 159 + 1)) + 12345;
        BigInt.DivRem(aSmallTop, dFull, out BigInt q3, out BigInt r3);
        Assert.Equal(aSmallTop, q3 * dFull + r3);
        Assert.True(BigInt.CompareAbs(r3, dFull) < 0);
    }

    [Fact]
    public void FloorAndCeilingDiv_Identities()
    {
        var rng = new Random(44);
        for (int i = 0; i < 1000; i++)
        {
            var (a, _) = TestHelpers.RandomPair(rng, 8);
            var (b, _) = TestHelpers.RandomPair(rng, 4);
            if (b.IsZero)
            {
                continue;
            }
            BigInt.DivRem(a, b, out BigInt fq, out BigInt fr, IntegerRounding.Floor);
            Assert.Equal(a, fq * b + fr);
            Assert.True(fr.IsZero || fr.Sign == b.Sign, $"floor rem sign #{i}");
            Assert.True(BigInt.CompareAbs(fr, b) < 0);

            BigInt.DivRem(a, b, out BigInt cq, out BigInt cr, IntegerRounding.Ceiling);
            Assert.Equal(a, cq * b + cr);
            Assert.True(cr.IsZero || cr.Sign == -b.Sign, $"ceil rem sign #{i}");
        }
    }

    [Fact]
    public void Shifts_MatchBigInteger()
    {
        var rng = new Random(45);
        for (int i = 0; i < 2000; i++)
        {
            var (a, ba) = TestHelpers.RandomPair(rng, 10);
            int sh = rng.Next(0, 300);
            TestHelpers.AssertEqual(ba << sh, a << sh, $"lshift #{i} by {sh}");
            TestHelpers.AssertEqual(ba >> sh, a >> sh, $"rshift #{i} by {sh}");
        }
    }

    [Fact]
    public void Bitwise_MatchesBigInteger()
    {
        var rng = new Random(46);
        for (int i = 0; i < 2000; i++)
        {
            var (a, ba) = TestHelpers.RandomPair(rng, 8);
            var (b, bb) = TestHelpers.RandomPair(rng, 8);
            TestHelpers.AssertEqual(ba & bb, a & b, $"and #{i}");
            TestHelpers.AssertEqual(ba | bb, a | b, $"or #{i}");
            TestHelpers.AssertEqual(ba ^ bb, a ^ b, $"xor #{i}");
            TestHelpers.AssertEqual(~ba, ~a, $"not #{i}");
        }
    }

    [Fact]
    public void TestBit_SetBit_MatchBigInteger()
    {
        var rng = new Random(47);
        for (int i = 0; i < 500; i++)
        {
            var (a, ba) = TestHelpers.RandomPair(rng, 5);
            int bit = rng.Next(0, 400);
            bool expected = ((ba >> bit) & 1) == 1;
            Assert.Equal(expected, a.TestBit(bit));
            TestHelpers.AssertEqual(ba | (BigInteger.One << bit), a.SetBit(bit), $"setbit #{i}");
            TestHelpers.AssertEqual(ba & ~(BigInteger.One << bit), a.ClearBit(bit), $"clrbit #{i}");
            TestHelpers.AssertEqual(ba ^ (BigInteger.One << bit), a.FlipBit(bit), $"flipbit #{i}");
        }
    }

    [Fact]
    public void Compare_MatchesBigInteger()
    {
        var rng = new Random(48);
        for (int i = 0; i < 2000; i++)
        {
            var (a, ba) = TestHelpers.RandomPair(rng, 6);
            var (b, bb) = TestHelpers.RandomPair(rng, 6);
            Assert.Equal(ba.CompareTo(bb), a.CompareTo(b));
            Assert.Equal(ba == bb, a == b);
        }
    }

    [Fact]
    public void Pow_MatchesBigInteger()
    {
        var rng = new Random(49);
        for (int i = 0; i < 200; i++)
        {
            var (a, ba) = TestHelpers.RandomPair(rng, 3);
            int e = rng.Next(0, 30);
            TestHelpers.AssertEqual(BigInteger.Pow(ba, e), BigInt.Pow(a, (ulong)e), $"pow #{i}");
        }
    }

    [Fact]
    public void Conversions_RoundTrip()
    {
        var rng = new Random(50);
        for (int i = 0; i < 1000; i++)
        {
            long v = (long)(rng.NextDouble() * long.MaxValue) * (rng.Next(2) == 0 ? 1 : -1);
            Assert.Equal(v, (long)new BigInt(v));
            ulong u = (ulong)(rng.NextDouble() * ulong.MaxValue);
            Assert.Equal(u, (ulong)new BigInt(u));
        }
        Assert.Equal(long.MinValue, (long)new BigInt(long.MinValue));
        Assert.Equal(long.MaxValue, (long)new BigInt(long.MaxValue));
        Assert.Equal(ulong.MaxValue, (ulong)new BigInt(ulong.MaxValue));
        Assert.Throws<OverflowException>(() => (long)(new BigInt(long.MaxValue) + BigInt.One));
        Assert.Equal(Int128.MinValue, (Int128)new BigInt(Int128.MinValue));
        Assert.Equal(UInt128.MaxValue, (UInt128)new BigInt(UInt128.MaxValue));
    }

    [Fact]
    public void DoubleConversion_TruncatesTowardZero()
    {
        Assert.Equal(new BigInt(3), new BigInt(3.99));
        Assert.Equal(new BigInt(-3), new BigInt(-3.99));
        Assert.Equal(BigInt.Zero, new BigInt(0.5));
        Assert.Equal(new BigInt(1L << 52), new BigInt(Math.Pow(2, 52)));
        var big = new BigInt(1.5e30);
        Assert.Equal(new BigInteger(1.5e30), big.ToBigInteger());
        Assert.Equal(1.5e30, big.ToDouble());
        Assert.Throws<OverflowException>(() => new BigInt(double.NaN));

        var rng = new Random(51);
        for (int i = 0; i < 500; i++)
        {
            var (a, ba) = TestHelpers.RandomPair(rng, 6);
            // BigInteger's own conversion is double-rounded via a 64-bit
            // intermediate, so allow a 1-ulp difference.
            double expected = (double)ba;
            double actual = a.ToDouble();
            long ulps = Math.Abs(BitConverter.DoubleToInt64Bits(Math.Abs(expected)) - BitConverter.DoubleToInt64Bits(Math.Abs(actual)));
            Assert.True(ulps <= 1, $"ToDouble #{i}: expected {expected:R}, actual {actual:R}");
        }
    }
}
