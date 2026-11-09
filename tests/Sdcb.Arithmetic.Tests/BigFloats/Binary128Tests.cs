using Xunit;

namespace Sdcb.Arithmetic.Tests;

public class Binary128Tests
{
    private static readonly RoundingMode[] Modes =
        [RoundingMode.NearestEven, RoundingMode.TowardZero, RoundingMode.TowardPositive, RoundingMode.TowardNegative, RoundingMode.AwayFromZero];

    private static UInt128 Hex(ulong hi, ulong lo) => new(hi, lo);

    private static readonly UInt128 PosInf = Hex(0x7FFF_0000_0000_0000, 0);
    private static readonly UInt128 MaxFinite = Hex(0x7FFE_FFFF_FFFF_FFFF, ulong.MaxValue);
    private static readonly BigFloat MinSubnormal = BigFloat.FromBigInt(BigInt.One, -16494, 2);

    [Theory]
    [InlineData(1.0, 0x3FFF_0000_0000_0000UL, 0UL)]
    [InlineData(-2.0, 0xC000_0000_0000_0000UL, 0UL)]
    [InlineData(0.5, 0x3FFE_0000_0000_0000UL, 0UL)]
    [InlineData(-0.0, 0x8000_0000_0000_0000UL, 0UL)]
    [InlineData(0.0, 0UL, 0UL)]
    [InlineData(double.PositiveInfinity, 0x7FFF_0000_0000_0000UL, 0UL)]
    [InlineData(double.NegativeInfinity, 0xFFFF_0000_0000_0000UL, 0UL)]
    [InlineData(0.1, 0x3FFB_9999_9999_9999UL, 0xA000_0000_0000_0000UL)]
    public void KnownEncodings(double value, ulong hi, ulong lo)
    {
        UInt128 bits = Hex(hi, lo);
        Assert.Equal(bits, ((BigFloat)value).ToBinary128Bits());
        BigFloat back = BigFloat.FromBinary128Bits(bits);
        Assert.Equal(BitConverter.DoubleToUInt64Bits(value), BitConverter.DoubleToUInt64Bits(back.ToDouble()));
    }

    [Fact]
    public void KnownConstants()
    {
        Assert.Equal(Hex(0x4000_921F_B544_42D1, 0x8469_898C_C517_01B8), BigFloat.Pi(300).ToBinary128Bits());
        Assert.Equal(Hex(0x3FFD_5555_5555_5555, 0x5555_5555_5555_5555), (BigFloat.One / BigFloat.FromBigInt(3, 300)).ToBinary128Bits());
        Assert.Equal(MaxFinite, BigFloat.FromBinary128Bits(MaxFinite).ToBinary128Bits());
        Assert.Equal(UInt128.One, MinSubnormal.ToBinary128Bits());
        Assert.Equal(Hex(0x0001_0000_0000_0000, 0), BigFloat.FromBigInt(BigInt.One, -16382, 2).ToBinary128Bits());
    }

    [Fact]
    public void SpecialValues()
    {
        Assert.True(BigFloat.FromBinary128Bits(Hex(0x7FFF_8000_0000_0000, 0)).IsNaN);
        Assert.True(BigFloat.FromBinary128Bits(Hex(0xFFFF_0000_0000_0000, 1)).IsNaN);
        Assert.Equal(Hex(0x7FFF_8000_0000_0000, 0), BigFloat.NaN.ToBinary128Bits());
        Assert.True(BigFloat.FromBinary128Bits(PosInf).IsPositiveInfinity);
        Assert.True(BigFloat.FromBinary128Bits(PosInf | (UInt128.One << 127)).IsNegativeInfinity);
        BigFloat negZero = BigFloat.FromBinary128Bits(UInt128.One << 127);
        Assert.True(negZero.IsZero && negZero.SignBit);
    }

    [Fact]
    public void DefaultPrecisionIsExact()
    {
        int saved = BigFloat.DefaultPrecision;
        try
        {
            BigFloat.DefaultPrecision = 53;
            BigFloat x = BigFloat.FromBinary128Bits(Hex(0x3FFF_0000_0000_0000, 1));
            Assert.Equal(113, x.Precision);
            Assert.Equal(Hex(0x3FFF_0000_0000_0000, 1), x.ToBinary128Bits());
            BigFloat r = BigFloat.FromBinary128Bits(Hex(0x3FFF_0000_0000_0000, 1), 53, RoundingMode.TowardPositive);
            Assert.Equal(53, r.Precision);
            Assert.Equal(0, BigFloat.FromBigInt((BigInt.One << 52) + 1, -52, 53).CompareTo(r));
        }
        finally
        {
            BigFloat.DefaultPrecision = saved;
        }
    }

    [Fact]
    public void RandomBitPatternsRoundTrip()
    {
        Random rng = new(128);
        Span<byte> buf = stackalloc byte[16];
        for (int i = 0; i < 20000; i++)
        {
            rng.NextBytes(buf);
            UInt128 bits = System.Buffers.Binary.BinaryPrimitives.ReadUInt128LittleEndian(buf);
            if (i % 4 == 1)
            {
                bits &=  (UInt128.One << 127) | ((UInt128.One << 112) - 1);
            }
            int biased = (int)(bits >> 112) & 0x7FFF;
            if (biased == 0x7FFF)
            {
                continue;
            }
            BigFloat x = BigFloat.FromBinary128Bits(bits);
            Assert.Equal(bits, x.ToBinary128Bits());
            Assert.Equal(bits, BigFloat.FromBinary128Bits(bits, 500).ToBinary128Bits());
        }
    }

    /// <summary>Oracle: round |x| to a multiple of 2^max(BinaryExponent - 113, -16494), then check the range.</summary>
    private static BigFloat Expected(BigFloat x, RoundingMode rnd, out bool overflow)
    {
        long ulp = Math.Max(x.BinaryExponent - 113, -16494);
        BigInt k = x.ShiftLeft(-ulp).ToBigInt(rnd);
        overflow = false;
        if (k.IsZero)
        {
            return x.SignBit ? BigFloat.NegativeZero : BigFloat.Zero;
        }
        BigFloat r = BigFloat.FromBigInt(k, ulp, (int)Math.Max(2, k.Abs().GetBitLength()));
        if (r.Abs().BinaryExponent > 16384)
        {
            overflow = true;
        }
        return r;
    }

    [Fact]
    public void RandomRoundingMatchesOracle()
    {
        Random rng = new(113);
        long[] exps = [-16700, -16550, -16500, -16494, -16490, -16400, -16382, -16300, -200, 0, 300, 16200, 16270, 16271, 16272, 16300];
        foreach (long baseExp in exps)
        {
            for (int i = 0; i < 300; i++)
            {
                int prec = rng.Next(2, 260);
                byte[] mb = new byte[(prec + 7) / 8 + 1];
                rng.NextBytes(mb);
                BigInt m = (BigInt.FromMagnitudeBytes(mb) & ((BigInt.One << prec) - BigInt.One)) | (BigInt.One << (prec - 1));
                if (i % 5 == 0)
                {
                    m = (m >> (prec - 1) << (prec - 1)) | (prec > 113 ? BigInt.One << (prec - 114) : BigInt.Zero); // exact ties at the normal boundary
                }
                bool neg = rng.Next(2) == 0;
                BigFloat x = BigFloat.FromBigInt(neg ? -m : m, baseExp + rng.Next(-120, 120) - prec + 113, prec);
                foreach (RoundingMode rnd in Modes)
                {
                    UInt128 bits = x.ToBinary128Bits(rnd);
                    BigFloat exp = Expected(x, rnd, out bool overflow);
                    if (overflow)
                    {
                        Assert.True(bits == (PosInf | (neg ? UInt128.One << 127 : 0)) || bits == (MaxFinite | (neg ? UInt128.One << 127 : 0)), $"{x} {rnd}");
                        bool toInf = rnd switch
                        {
                            RoundingMode.TowardZero => false,
                            RoundingMode.TowardPositive => !neg,
                            RoundingMode.TowardNegative => neg,
                            _ => true,
                        };
                        Assert.Equal(toInf, (bits & ~(UInt128.One << 127)) == PosInf);
                        continue;
                    }
                    BigFloat got = BigFloat.FromBinary128Bits(bits);
                    Assert.True(exp.CompareTo(got) == 0 && exp.SignBit == got.SignBit, $"x={x} (prec {prec}) {rnd}: expected {exp}, got {got}");
                }
            }
        }
    }

    [Fact]
    public void TiesAtSubnormalBoundary()
    {
        BigFloat half = MinSubnormal.ShiftRight(1);
        Assert.Equal(UInt128.Zero, half.ToBinary128Bits());
        Assert.Equal(UInt128.One, half.ToBinary128Bits(RoundingMode.AwayFromZero));
        Assert.Equal(UInt128.One << 127, (-half).ToBinary128Bits(RoundingMode.TowardPositive));
        BigFloat threeHalves = MinSubnormal * BigFloat.FromDouble(1.5);
        Assert.Equal((UInt128)2, threeHalves.ToBinary128Bits());
        BigFloat big = BigFloat.FromBigInt(BigInt.One, 16384, 2);
        Assert.Equal(PosInf, big.ToBinary128Bits());
        Assert.Equal(MaxFinite, big.ToBinary128Bits(RoundingMode.TowardZero));
    }
}
