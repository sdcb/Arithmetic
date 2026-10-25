using System.Numerics;
using Xunit;

namespace Sdcb.Arithmetic.Tests;

/// <summary>
/// Kernel cross-checks against BigInteger at the sizes where SIMD and NTT paths engage. Every
/// vector tier must give bit-identical results: CI reruns the suite with
/// DOTNET_PreferredVectorBitWidth=128 and DOTNET_EnableHWIntrinsic=0 to cover the narrower tiers.
/// </summary>
public class KernelConsistencyTests
{
    [Fact]
    public void LinearKernels_MatchBigInteger()
    {
        var rng = new Random(900);
        for (int i = 0; i < 200; i++)
        {
            var (a, ba) = TestHelpers.RandomPair(rng, 90);
            var (b, bb) = TestHelpers.RandomPair(rng, 90);
            TestHelpers.AssertEqual(ba + bb, a + b, $"add #{i}");
            TestHelpers.AssertEqual(ba - bb, a - b, $"sub #{i}");
            TestHelpers.AssertEqual(ba & bb, a & b, $"and #{i}");
            TestHelpers.AssertEqual(ba | bb, a | b, $"or #{i}");
            TestHelpers.AssertEqual(ba ^ bb, a ^ b, $"xor #{i}");
            int s = rng.Next(1, 300);
            TestHelpers.AssertEqual(ba << s, a << s, $"shl #{i}");
            TestHelpers.AssertEqual(ba >> s, a >> s, $"shr #{i}");
            BigInt pa = a.Abs(), pb = b.Abs();
            long expected = 0;
            foreach (byte x in BigInteger.Abs(ba).ToByteArray())
            {
                expected += BitOperations.PopCount(x);
            }
            Assert.Equal(expected, pa.PopCount());
            Assert.Equal(expected + pb.PopCount() - 2 * (pa & pb).PopCount(), BigInt.HammingDistance(pa, pb));
        }
    }

    [Fact]
    public void BasecaseMulSqr_AllShapes_MatchBigInteger()
    {
        // Covers the scalar and AVX2 basecases (and the Karatsuba/Toom splits above them),
        // including all-ones operands, which maximize every column sum.
        var rng = new Random(905);
        static byte[] Limbs(Random rng, int n, int kind)
        {
            byte[] x = new byte[n * 8];
            if (kind == 0)
            {
                rng.NextBytes(x);
            }
            else
            {
                x.AsSpan().Fill(0xFF);
            }
            x[^1] |= 0x80;
            return x;
        }
        for (int an = 1; an <= 140; an += an < 40 ? 1 : 7)
        {
            for (int bn = 1; bn <= an; bn += bn < 20 ? 1 : 5)
            {
                for (int kind = 0; kind < 2; kind++)
                {
                    byte[] x = Limbs(rng, an, kind), y = Limbs(rng, bn, kind);
                    BigInt a = BigInt.FromMagnitudeBytes(x), b = BigInt.FromMagnitudeBytes(y);
                    BigInteger ba = new(x, isUnsigned: true), bb = new(y, isUnsigned: true);
                    TestHelpers.AssertEqual(ba * bb, a * b, $"mul {an}x{bn} kind {kind}");
                }
            }
        }
        for (int n = 1; n <= 400; n += n < 64 ? 1 : 13)
        {
            for (int kind = 0; kind < 2; kind++)
            {
                byte[] x = Limbs(rng, n, kind);
                BigInt a = BigInt.FromMagnitudeBytes(x);
                BigInteger ba = new(x, isUnsigned: true);
                TestHelpers.AssertEqual(ba * ba, a.Square(), $"sqr {n} kind {kind}");
                TestHelpers.AssertEqual(ba * ba, a * a, $"self-mul {n} kind {kind}");
            }
        }
    }

    [Fact]
    public void ToomSizes_MatchBigInteger()
    {
        // Toom-3 split shapes (n mod 3) up to the NTT crossover, all-ones operands and zero parts.
        var rng = new Random(925);
        foreach (int n in (int[])[384, 385, 600, 601, 602, 603, 700, 701, 1000, 1203, 2500, 3299])
        {
            for (int kind = 0; kind < 3; kind++)
            {
                byte[] x = new byte[n * 8], y = new byte[n * 8];
                rng.NextBytes(x);
                rng.NextBytes(y);
                if (kind == 1)
                {
                    x.AsSpan().Fill(0xFF);
                    y.AsSpan().Fill(0xFF);
                }
                else if (kind == 2)
                {
                    x.AsSpan(x.Length / 2, x.Length / 2 - 8).Clear();
                    y.AsSpan(8, y.Length / 3).Clear();
                }
                x[^1] |= 1;
                y[^1] |= 1;
                BigInt a = BigInt.FromMagnitudeBytes(x), b = BigInt.FromMagnitudeBytes(y);
                BigInteger ba = new(x, isUnsigned: true), bb = new(y, isUnsigned: true);
                TestHelpers.AssertEqual(ba * bb, a * b, $"mul {n} kind {kind}");
                TestHelpers.AssertEqual(ba * ba, a.Square(), $"sqr {n} kind {kind}");
            }
        }
    }

    [Fact]
    public void DivExact_AllShapes()
    {
        // Basecase, divide-and-conquer and Hensel paths; even divisors and zero low limbs.
        var rng = new Random(915);
        (int qn, int dn)[] shapes = [(1, 1), (5, 3), (24, 24), (25, 25), (26, 3), (49, 30), (100, 100), (100, 7),
            (257, 256), (300, 1000), (1025, 1024), (1500, 40), (2049, 2048), (5000, 3000), (21000, 20000)];
        foreach (var (qn, dn) in shapes)
        {
            for (int kind = 0; kind < (qn > 5000 ? 1 : 3); kind++)
            {
                byte[] x = new byte[qn * 8], y = new byte[dn * 8];
                rng.NextBytes(x);
                rng.NextBytes(y);
                if (kind == 1)
                {
                    y[0] &= 0xF0; // even divisor
                }
                else if (kind == 2)
                {
                    x.AsSpan().Fill(0xFF);
                    y.AsSpan(0, Math.Min(16, y.Length - 8)).Clear(); // zero low limbs
                }
                y[^1] |= 1;
                BigInt q = BigInt.FromMagnitudeBytes(x), d = BigInt.FromMagnitudeBytes(y);
                BigInt a = q * d;
                Assert.Equal(q, BigInt.DivExact(a, d));
                Assert.Equal(-q, BigInt.DivExact(-a, d));
            }
        }
    }

    [Fact]
    public void LargeDivision_Consistent()
    {
        var rng = new Random(970);
        (int an, int dn)[] shapes = [(3000, 1500), (2600, 1300), (5000, 700), (1800, 1700), (4100, 2000), (1400, 90), (700, 650), (1300, 640)];
        foreach (var (an, dn) in shapes)
        {
            byte[] x = new byte[an * 8], y = new byte[dn * 8];
            rng.NextBytes(x);
            rng.NextBytes(y);
            if (rng.Next(2) == 0)
            {
                y.AsSpan(8, y.Length / 2).Fill(0xFF); // divisor near a power of two
            }
            if (rng.Next(3) == 0)
            {
                y[^1] = 1; // badly normalized divisor
            }
            BigInt a = BigInt.FromMagnitudeBytes(x), d = BigInt.FromMagnitudeBytes(y);
            BigInt.DivRem(a, d, out BigInt q, out BigInt r);
            Assert.True(r.Sign >= 0 && r < d, $"remainder range {an}/{dn}");
            Assert.Equal(a, q * d + r);
            BigInteger bq = BigInteger.DivRem(new BigInteger(x, isUnsigned: true), new BigInteger(y, isUnsigned: true), out _);
            TestHelpers.AssertEqual(bq, q, $"quotient {an}/{dn}");
            // Exact multiples.
            BigInt m = q * d;
            Assert.Equal(q, BigInt.DivExact(m, d));
        }
    }

    [Fact]
    public void Division_BlockBoundaryShapes()
    {
        var rng = new Random(980);
        foreach (int dn in (int[])[40, 41, 64, 65, 100, 129, 256, 300])
        {
            foreach (int extra in (int[])[0, 1, 2, 5, dn / 8, dn / 8 + 3, dn / 2, dn - 1, dn + 1])
            {
                int an = 2 * dn + extra;
                byte[] x = new byte[an * 8 - rng.Next(8)], y = new byte[dn * 8 - rng.Next(8)];
                rng.NextBytes(x);
                rng.NextBytes(y);
                BigInt a = BigInt.FromMagnitudeBytes(x), d = BigInt.FromMagnitudeBytes(y);
                BigInt.DivRem(a, d, out BigInt q, out BigInt r);
                BigInteger bq = BigInteger.DivRem(new BigInteger(x, isUnsigned: true), new BigInteger(y, isUnsigned: true), out BigInteger br);
                TestHelpers.AssertEqual(bq, q, $"q {an}/{dn}");
                TestHelpers.AssertEqual(br, r, $"r {an}/{dn}");
            }
        }
    }

    [Fact]
    public void LargeMultiplication_MatchesBigInteger()
    {
        var rng = new Random(950);
        foreach (int limbs in (int[])[1600, 2100, 4000, 9000])
        {
            byte[] x = new byte[limbs * 8], y = new byte[(limbs - rng.Next(limbs / 2)) * 8];
            rng.NextBytes(x);
            rng.NextBytes(y);
            if (limbs == 2100)
            {
                x.AsSpan(100, 4000).Fill(0xFF); // long carry runs
                y.AsSpan(0, 3000).Fill(0xFF);
            }
            BigInt a = BigInt.FromMagnitudeBytes(x), b = BigInt.FromMagnitudeBytes(y);
            BigInteger ba = new(x, isUnsigned: true), bb = new(y, isUnsigned: true);
            TestHelpers.AssertEqual(ba * bb, a * b, $"mul {limbs}");
            TestHelpers.AssertEqual(ba * ba, a.Square(), $"sqr {limbs}");
            BigInt allOnes = (BigInt.One << (64 * limbs)) - 1;
            BigInteger bOnes = (BigInteger.One << (64 * limbs)) - 1;
            TestHelpers.AssertEqual(bOnes * bOnes, allOnes * allOnes, $"max {limbs}");
        }
    }
}
