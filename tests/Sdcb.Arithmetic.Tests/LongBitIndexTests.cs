using Xunit;

namespace Sdcb.Arithmetic.Tests;

/// <summary>Bit indices / shift counts whose limb index exceeds int range must not wrap (issue #5).</summary>
public class LongBitIndexTests
{
    [Theory]
    [InlineData(1, 1L << 38, false)] // limb index used to wrap to 0
    [InlineData(-2, 1L << 38, true)] // sign extension
    [InlineData(0, 1L << 37, false)] // limb index used to wrap negative
    [InlineData(-1, 1L << 37, true)]
    [InlineData(5, long.MaxValue, false)]
    [InlineData(-5, long.MaxValue, true)]
    public void TestBit_HighIndex(int value, long index, bool expected)
    {
        Assert.Equal(expected, new BigInt(value).TestBit(index));
    }

    [Fact]
    public void SetClearBit_NoOpOnHighIndex_PreservesValue()
    {
        Assert.Equal(BigInt.One, BigInt.One.ClearBit(1L << 38));
        Assert.Equal(new BigInt(-2), new BigInt(-2).SetBit(1L << 38));
        Assert.Equal(new BigInt(-2), new BigInt(-2).SetBit(long.MaxValue));
    }

    [Fact]
    public void ChangingUnrepresentableBit_Throws()
    {
        Assert.Throws<OverflowException>(() => BigInt.One.SetBit(1L << 38));
        Assert.Throws<OverflowException>(() => BigInt.One.FlipBit(1L << 38));
        Assert.Throws<OverflowException>(() => new BigInt(-2).ClearBit(1L << 38));
    }

    [Fact]
    public void IndexOfBit_HighStart()
    {
        Assert.Equal(-1, BigInt.One.IndexOfSetBit(1L << 38));
        Assert.Equal(1L << 38, BigInt.One.IndexOfClearBit(1L << 38));
        Assert.Equal(1L << 38, new BigInt(-2).IndexOfSetBit(1L << 38));
        Assert.Equal(-1, new BigInt(-2).IndexOfClearBit(1L << 38));
    }

    [Fact]
    public void ShiftLeft_UnrepresentableCount_Throws()
    {
        Assert.Throws<OverflowException>(() => BigInt.One.ShiftLeft(1L << 38));
        Assert.Throws<OverflowException>(() => BigInt.MinusOne.ShiftLeft(long.MaxValue));
    }

    [Fact]
    public void ShiftRight_HighCount()
    {
        Assert.Equal(BigInt.Zero, new BigInt(3).ShiftRight((1L << 37) + 1));
        Assert.Equal(BigInt.MinusOne, new BigInt(-3).ShiftRight((1L << 37) + 1));
        Assert.Equal(BigInt.One, new BigInt(3).ShiftRight((1L << 37) + 1, IntegerRounding.Ceiling));
    }

    [Fact]
    public void RemainderPow2_HighBits()
    {
        // |value| < 2^bits: the remainder is the value itself (or zero) and needs no 2^bits allocation.
        Assert.Equal(BigInt.One, BigInt.One.ModPow2(1L << 40));
        Assert.Equal(BigInt.MinusOne, BigInt.MinusOne.RemainderPow2(1L << 40, IntegerRounding.Ceiling));
        Assert.Equal(BigInt.MinusOne, BigInt.MinusOne.RemainderPow2(1L << 40, IntegerRounding.Truncate));
        // Otherwise the result is 2^bits - |value|, which cannot be represented: throw instead of truncating bits.
        Assert.Throws<OverflowException>(() => BigInt.MinusOne.ModPow2(1L << 40));
        Assert.Throws<OverflowException>(() => BigInt.MinusOne.RemainderPow2(1L << 40, IntegerRounding.Floor));
        Assert.Throws<OverflowException>(() => BigInt.One.RemainderPow2(1L << 40, IntegerRounding.Ceiling));
    }
}
