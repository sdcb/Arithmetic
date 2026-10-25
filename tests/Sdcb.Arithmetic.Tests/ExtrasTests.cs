using System.Numerics;
using Xunit;

namespace Sdcb.Arithmetic.Tests;

/// <summary>mpz_*_2exp variants, congruences, double comparisons and prime stepping.</summary>
public class ExtrasTests
{
    [Fact]
    public void PowerOfTwoDivisions_MatchBigIntegerSemantics()
    {
        var rng = new Random(500);
        for (int i = 0; i < 400; i++)
        {
            var (a, ba) = TestHelpers.RandomPair(rng, 6);
            int k = rng.Next(0, 420);
            BigInteger p = BigInteger.One << k;
            BigInteger tq = BigInteger.Divide(ba, p), tr = ba - tq * p;
            BigInteger fq = tr.Sign < 0 ? tq - 1 : tq, fr = ba - fq * p;
            BigInteger cq = tr.Sign > 0 ? tq + 1 : tq, cr = ba - cq * p;
            TestHelpers.AssertEqual(tr, a.RemainderPow2(k, IntegerRounding.Truncate), $"tdiv_r_2exp {i}");
            TestHelpers.AssertEqual(fr, a.ModPow2(k), $"fdiv_r_2exp {i}");
            TestHelpers.AssertEqual(cq, a.ShiftRight(k, IntegerRounding.Ceiling), $"cdiv_q_2exp {i}");
            TestHelpers.AssertEqual(cr, a.RemainderPow2(k, IntegerRounding.Ceiling), $"cdiv_r_2exp {i}");
            Assert.Equal(tr.IsZero, a.IsDivisibleByPow2(k));
            TestHelpers.AssertEqual(fr, a.RemainderPow2(k, IntegerRounding.Floor), $"floor r_2exp {i}");
            TestHelpers.AssertEqual(fq, a.ShiftRight(k), $"default ShiftRight is floor {i}");
            TestHelpers.AssertEqual(tq, a.ShiftRight(k, IntegerRounding.Truncate), $"tdiv_q_2exp {i}");
            Assert.Equal(a >> k, a.ShiftRight(k));
        }
    }

    [Fact]
    public void DivideAndRemainder_AllRoundings()
    {
        var rng = new Random(501);
        for (int i = 0; i < 300; i++)
        {
            var (a, ba) = TestHelpers.RandomPair(rng, 5);
            var (b, bb) = TestHelpers.RandomPair(rng, 3);
            if (bb.IsZero) continue;
            BigInteger tq = BigInteger.Divide(ba, bb), tr = ba - tq * bb;
            BigInteger fq = tr.Sign != 0 && (tr.Sign < 0) != (bb.Sign < 0) ? tq - 1 : tq;
            BigInteger cq = tr.Sign != 0 && (tr.Sign < 0) == (bb.Sign < 0) ? tq + 1 : tq;
            foreach ((IntegerRounding mode, BigInteger q) in (ReadOnlySpan<(IntegerRounding, BigInteger)>)[(IntegerRounding.Truncate, tq), (IntegerRounding.Floor, fq), (IntegerRounding.Ceiling, cq)])
            {
                TestHelpers.AssertEqual(q, BigInt.Divide(a, b, mode), $"{mode} q {i}");
                TestHelpers.AssertEqual(ba - q * bb, BigInt.Remainder(a, b, mode), $"{mode} r {i}");
            }
        }
    }

    [Fact]
    public void Congruences()
    {
        Assert.True(BigInt.IsCongruent(17, 5, 12));
        Assert.True(BigInt.IsCongruent(-7, 5, 12));
        Assert.False(BigInt.IsCongruent(17, 6, 12));
        Assert.True(BigInt.IsCongruent(9, 9, 0));
        Assert.False(BigInt.IsCongruent(9, 8, 0));
        Assert.True(BigInt.IsCongruentModPow2(BigInt.One << 70, BigInt.Zero, 70));
        Assert.False(BigInt.IsCongruentModPow2((BigInt.One << 70) + 1, BigInt.Zero, 70));
    }

    [Fact]
    public void DoubleComparisonsAndDecomposition()
    {
        BigInt big = BigInt.One << 80;
        Assert.Equal(0, big.CompareTo(Math.Pow(2, 80)));
        Assert.Equal(1, (big + 1).CompareTo(Math.Pow(2, 80)));
        Assert.Equal(-1, new BigInt(5).CompareTo(5.5));
        Assert.Equal(1, new BigInt(-5).CompareTo(-5.5));
        Assert.Equal(-1, new BigInt(-5).CompareTo(double.PositiveInfinity));
        Assert.Equal(1, BigInt.One.CompareTo(double.NegativeInfinity));
        double m = new BigInt(12).ToDouble(out long e);
        Assert.Equal(0.75, m);
        Assert.Equal(4, e);
        m = (-(BigInt.One << 100) - 1).ToDouble(out e);
        Assert.Equal(-0.5, m); // truncated toward zero
        Assert.Equal(101, e);
        Assert.True(new BigInt(long.MinValue).FitsInt64);
        Assert.False((new BigInt(long.MaxValue) + 1).FitsInt64);
        Assert.True(new BigInt(uint.MaxValue).FitsUInt32);
        Assert.False(new BigInt(-1).FitsUInt64);
    }

    [Fact]
    public void NextAndPrevPrime()
    {
        Assert.Equal(new BigInt(2), BigInt.Zero.NextPrime());
        Assert.Equal(new BigInt(3), new BigInt(2).NextPrime());
        Assert.Equal(new BigInt(101), new BigInt(97).NextPrime());
        Assert.Equal(new BigInt(97), new BigInt(101).PrevPrime());
        Assert.False(new BigInt(2).TryPrevPrime(out _));
        Assert.Equal(new BigInt(2), new BigInt(3).PrevPrime());
        // Around 2^64: 2^64 + 13 is the next prime, 2^64 - 59 the previous one.
        BigInt p64 = BigInt.One << 64;
        Assert.Equal(p64 + 13, p64.NextPrime());
        Assert.Equal(p64 - 59, p64.PrevPrime());
        // Large: next prime after 10^100 is 10^100 + 267.
        BigInt g = BigInt.Pow(10, 100);
        Assert.Equal(g + 267, g.NextPrime());
    }
}
