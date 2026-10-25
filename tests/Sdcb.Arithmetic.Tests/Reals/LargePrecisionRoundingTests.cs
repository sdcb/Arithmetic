using Xunit;

namespace Sdcb.Arithmetic.Tests;

/// <summary>
/// Multi-limb division and square root against exact rational / integer oracles, in every
/// rounding mode, including exactly representable results and results at or next to a tie.
/// These sizes take the remainder-free fast paths.
/// </summary>
public class LargePrecisionRoundingTests
{
    private static readonly RealRounding[] Modes =
        [RealRounding.NearestEven, RealRounding.TowardZero, RealRounding.TowardPositive, RealRounding.TowardNegative, RealRounding.AwayFromZero];

    private static readonly int[] Precs = [300, 640, 1000, 1300, 2000, 4097, 9000];

    private static BigInt RandomBits(Random rng, int bits)
    {
        byte[] x = new byte[(bits + 7) / 8 + 1];
        rng.NextBytes(x);
        BigInt v = BigInt.FromMagnitudeBytes(x) & ((BigInt.One << bits) - BigInt.One);
        return v | (BigInt.One << (bits - 1));
    }

    private static Real Exact(BigInt m, long exp = 0) => Real.FromBigInt(m, exp, (int)Math.Max(2, m.Abs().GetBitLength()));

    private static void AssertSame(Real expected, Real actual, string msg)
    {
        Assert.True(expected.Equals(actual) && expected.Precision == actual.Precision, $"{msg}: expected {expected}, got {actual}");
    }

    [Fact]
    public void Div_MatchesRationalOracle()
    {
        var rng = new Random(1207);
        foreach (int p in Precs)
        {
            for (int i = 0; i < 6; i++)
            {
                BigInt ym = RandomBits(rng, p + rng.Next(-200, 200));
                BigInt qm = RandomBits(rng, p + (i < 3 ? 0 : 1));
                // i % 3: exact quotient, just above, just below (q*y ± 1).
                BigInt xm = qm * ym + (i % 3 == 0 ? BigInt.Zero : i % 3 == 1 ? BigInt.One : BigInt.MinusOne);
                Real x = Exact(rng.Next(2) == 0 ? xm : -xm, rng.Next(-100, 100)), y = Exact(ym, rng.Next(-100, 100));
                foreach (RealRounding rnd in Modes)
                {
                    Real expected = Real.FromRational(x.ToRational() / y.ToRational(), p, rnd);
                    AssertSame(expected, Real.Divide(x, y, p, rnd), $"div p={p} case={i} {rnd}");
                }
            }
        }
    }

    [Fact]
    public void Sqrt_MatchesIntegerOracle()
    {
        var rng = new Random(1208);
        foreach (int p in Precs)
        {
            for (int i = 0; i < 6; i++)
            {
                // i % 3: perfect square, just above, just below.
                BigInt r = RandomBits(rng, p + (i < 3 ? 0 : 1));
                BigInt nm = r * r + (i % 3 == 0 ? BigInt.Zero : i % 3 == 1 ? BigInt.One : BigInt.MinusOne);
                long e = 2 * rng.Next(-100, 100);
                Real x = Exact(nm, e);
                foreach (RealRounding rnd in Modes)
                {
                    // sqrt(nm * 2^e) at precision p: floor(sqrt(nm * 4^k)) with k large plus sticky.
                    long k = p + 8;
                    BigInt scaled = nm << (int)(2 * k);
                    BigInt root = BigInt.Sqrt(scaled, out BigInt rem);
                    Real expected = Real.FromRational(new BigRational(root * 2 + (rem.IsZero ? BigInt.Zero : BigInt.One), BigInt.One << (int)(k + 1)) * (e >= 0 ? new BigRational(BigInt.One << (int)(e / 2)) : new BigRational(BigInt.One, BigInt.One << (int)(-e / 2))), p, rnd);
                    AssertSame(expected, Real.Sqrt(x, p, rnd), $"sqrt p={p} case={i} {rnd}");
                }
            }
        }
    }
}
