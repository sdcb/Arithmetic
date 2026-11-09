using Xunit;

namespace Sdcb.Arithmetic.Tests.BigFloats;

/// <summary>long exponents / root degrees must not be truncated or wrap (issue #6).</summary>
public class BigFloatLongArgumentRegressionTests
{
    [Theory]
    [InlineData(1, long.MinValue)] // normalizing needs long.MinValue - 1
    [InlineData(4, long.MaxValue)] // normalizing needs long.MaxValue + 1
    public void FromBigInt_ExponentMustNotWrap(int mantissa, long exponent)
    {
        Assert.Throws<OverflowException>(() => BigFloat.FromBigInt(new BigInt(mantissa), exponent, 2));
    }

    [Fact]
    public void FromBigInt_ExtremeExponentWithinRange()
    {
        BigFloat x = BigFloat.FromBigInt(BigInt.One, -long.MaxValue, 2);
        Assert.Equal(-long.MaxValue + 1, x.BinaryExponent);
    }

    // A large even degree must not be clamped to the odd uint.MaxValue.
    [Theory]
    [InlineData(1L << 32)]
    [InlineData((1L << 32) + 2)]
    public void RootN_EvenLongDegree_NegativeInfinityIsNaN(long degree)
    {
        Assert.True(BigFloat.RootN(BigFloat.NegativeInfinity, degree, 64).IsNaN);
    }

    [Theory]
    [InlineData(1L << 32)]
    [InlineData((1L << 32) + 2)]
    public void RootN_EvenLongDegree_NegativeZeroBecomesPositive(long degree)
    {
        BigFloat result = BigFloat.RootN(BigFloat.NegativeZero, degree, 64);
        Assert.True(result.IsZero);
        Assert.False(result.SignBit);
    }

    [Theory]
    [InlineData(1L << 32)]
    [InlineData((1L << 32) + 1)]
    [InlineData(-(1L << 33) - 1)]
    public void RootN_HugeDegree_MatchesExpLog(long degree)
    {
        BigFloat x = BigFloat.FromDouble(12345.678, 128);
        BigFloat expected = BigFloat.Exp(BigFloat.Divide(BigFloat.Log(x, 256), BigFloat.FromBigInt(new BigInt(degree), 64), 256), 128);
        Assert.Equal(expected, BigFloat.RootN(x, degree, 128));
    }

    [Fact]
    public void RootN_HugeOddDegree_NegativeArgument()
    {
        long degree = (1L << 32) + 1;
        Assert.Equal(-BigFloat.RootN(BigFloat.FromDouble(3, 64), degree, 64), BigFloat.RootN(BigFloat.FromDouble(-3, 64), degree, 64));
        Assert.True(BigFloat.RootN(BigFloat.FromDouble(-3, 64), degree + 1, 64).IsNaN);
    }

    [Fact]
    public void RootN_HugeDegree_ExactPowerOfTwo()
    {
        long degree = 1L << 40;
        BigFloat x = BigFloat.FromBigInt(BigInt.One, 3 * degree, 64); // 2^(3*2^40)
        Assert.Equal(BigFloat.FromBigInt(new BigInt(8), 64), BigFloat.RootN(x, degree, 64));
        Assert.Equal(BigFloat.FromBigInt(BigInt.One, -3, 64), BigFloat.RootN(x, -degree, 64));
    }

    [Fact]
    public void RootN_MinimumLongDegree_DoesNotUseFalseExactnessShortcut()
    {
        // |k| = 2^63 used to be clamped and cast back to a signed value, taking the exact branch.
        // x^(1/k) = 2^(1 - 2^-63), which rounds to 2 at 2 bits (nearest-even).
        BigFloat x = BigFloat.FromBigInt(BigInt.One, -long.MaxValue, 2);
        Assert.Equal(BigFloat.FromBigInt(new BigInt(2), 2), BigFloat.RootN(x, long.MinValue, 2));
    }

    [Fact]
    public void RootN_MinimumLongDegree_ExactPowerOfTwo()
    {
        BigFloat x = BigFloat.FromBigInt(BigInt.One, long.MinValue + 1, 2); // 2^(-2^63)
        Assert.Equal(BigFloat.FromBigInt(new BigInt(2), 2), BigFloat.RootN(x, long.MinValue, 2));
    }

    // k = 2 delegated to Sqrt and kept sqrt's -0 rule.
    [Fact]
    public void RootN_DegreeTwo_NegativeZeroBecomesPositive()
    {
        BigFloat result = BigFloat.RootN(BigFloat.NegativeZero, 2L, 64);
        Assert.True(result.IsZero);
        Assert.False(result.SignBit);
        Assert.True(BigFloat.Sqrt(BigFloat.NegativeZero, 64).SignBit); // sqrt(-0) stays -0
    }
}
