using Xunit;

namespace Sdcb.Arithmetic.Tests.Reals;

/// <summary>long exponents / root degrees must not be truncated or wrap (issue #6).</summary>
public class RealLongArgumentRegressionTests
{
    [Theory]
    [InlineData(1, long.MinValue)] // normalizing needs long.MinValue - 1
    [InlineData(4, long.MaxValue)] // normalizing needs long.MaxValue + 1
    public void FromBigInt_ExponentMustNotWrap(int mantissa, long exponent)
    {
        Assert.Throws<OverflowException>(() => Real.FromBigInt(new BigInt(mantissa), exponent, 2));
    }

    [Fact]
    public void FromBigInt_ExtremeExponentWithinRange()
    {
        Real x = Real.FromBigInt(BigInt.One, -long.MaxValue, 2);
        Assert.Equal(-long.MaxValue + 1, x.BinaryExponent);
    }

    // A large even degree must not be clamped to the odd uint.MaxValue.
    [Theory]
    [InlineData(1L << 32)]
    [InlineData((1L << 32) + 2)]
    public void RootN_EvenLongDegree_NegativeInfinityIsNaN(long degree)
    {
        Assert.True(Real.RootN(Real.NegativeInfinity, degree, 64).IsNaN);
    }

    [Theory]
    [InlineData(1L << 32)]
    [InlineData((1L << 32) + 2)]
    public void RootN_EvenLongDegree_NegativeZeroBecomesPositive(long degree)
    {
        Real result = Real.RootN(Real.NegativeZero, degree, 64);
        Assert.True(result.IsZero);
        Assert.False(result.SignBit);
    }

    [Theory]
    [InlineData(1L << 32)]
    [InlineData((1L << 32) + 1)]
    [InlineData(-(1L << 33) - 1)]
    public void RootN_HugeDegree_MatchesExpLog(long degree)
    {
        Real x = Real.FromDouble(12345.678, 128);
        Real expected = Real.Exp(Real.Divide(Real.Log(x, 256), Real.FromBigInt(new BigInt(degree), 64), 256), 128);
        Assert.Equal(expected, Real.RootN(x, degree, 128));
    }

    [Fact]
    public void RootN_HugeOddDegree_NegativeArgument()
    {
        long degree = (1L << 32) + 1;
        Assert.Equal(-Real.RootN(Real.FromDouble(3, 64), degree, 64), Real.RootN(Real.FromDouble(-3, 64), degree, 64));
        Assert.True(Real.RootN(Real.FromDouble(-3, 64), degree + 1, 64).IsNaN);
    }

    [Fact]
    public void RootN_HugeDegree_ExactPowerOfTwo()
    {
        long degree = 1L << 40;
        Real x = Real.FromBigInt(BigInt.One, 3 * degree, 64); // 2^(3*2^40)
        Assert.Equal(Real.FromBigInt(new BigInt(8), 64), Real.RootN(x, degree, 64));
        Assert.Equal(Real.FromBigInt(BigInt.One, -3, 64), Real.RootN(x, -degree, 64));
    }

    [Fact]
    public void RootN_MinimumLongDegree_DoesNotUseFalseExactnessShortcut()
    {
        // |k| = 2^63 used to be clamped and cast back to a signed value, taking the exact branch.
        // x^(1/k) = 2^(1 - 2^-63), which rounds to 2 at 2 bits (nearest-even).
        Real x = Real.FromBigInt(BigInt.One, -long.MaxValue, 2);
        Assert.Equal(Real.FromBigInt(new BigInt(2), 2), Real.RootN(x, long.MinValue, 2));
    }

    [Fact]
    public void RootN_MinimumLongDegree_ExactPowerOfTwo()
    {
        Real x = Real.FromBigInt(BigInt.One, long.MinValue + 1, 2); // 2^(-2^63)
        Assert.Equal(Real.FromBigInt(new BigInt(2), 2), Real.RootN(x, long.MinValue, 2));
    }

    // k = 2 delegated to Sqrt and kept sqrt's -0 rule.
    [Fact]
    public void RootN_DegreeTwo_NegativeZeroBecomesPositive()
    {
        Real result = Real.RootN(Real.NegativeZero, 2L, 64);
        Assert.True(result.IsZero);
        Assert.False(result.SignBit);
        Assert.True(Real.Sqrt(Real.NegativeZero, 64).SignBit); // sqrt(-0) stays -0
    }
}
