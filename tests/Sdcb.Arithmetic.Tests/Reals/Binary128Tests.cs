using Xunit;

namespace Sdcb.Arithmetic.Tests;

public class Binary128Tests
{
    private static readonly RealRounding[] Modes =
        [RealRounding.NearestEven, RealRounding.TowardZero, RealRounding.TowardPositive, RealRounding.TowardNegative, RealRounding.AwayFromZero];

    private static UInt128 Hex(ulong hi, ulong lo) => new(hi, lo);

    private static readonly UInt128 PosInf = Hex(0x7FFF_0000_0000_0000, 0);
    private static readonly UInt128 MaxFinite = Hex(0x7FFE_FFFF_FFFF_FFFF, ulong.MaxValue);
    private static readonly Real MinSubnormal = Real.FromBigInt(BigInt.One, -16494, 2);

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
        Assert.Equal(bits, ((Real)value).ToBinary128Bits());
        Real back = Real.FromBinary128Bits(bits);
        Assert.Equal(BitConverter.DoubleToUInt64Bits(value), BitConverter.DoubleToUInt64Bits(back.ToDouble()));
    }

    [Fact]
    public void KnownConstants()
    {
        Assert.Equal(Hex(0x4000_921F_B544_42D1, 0x8469_898C_C517_01B8), Real.Pi(300).ToBinary128Bits());
        Assert.Equal(Hex(0x3FFD_5555_5555_5555, 0x5555_5555_5555_5555), (Real.One / Real.FromBigInt(3, 300)).ToBinary128Bits());
        Assert.Equal(MaxFinite, Real.FromBinary128Bits(MaxFinite).ToBinary128Bits());
        Assert.Equal(UInt128.One, MinSubnormal.ToBinary128Bits());
        Assert.Equal(Hex(0x0001_0000_0000_0000, 0), Real.FromBigInt(BigInt.One, -16382, 2).ToBinary128Bits());
    }

    [Fact]
    public void SpecialValues()
    {
        Assert.True(Real.FromBinary128Bits(Hex(0x7FFF_8000_0000_0000, 0)).IsNaN);
        Assert.True(Real.FromBinary128Bits(Hex(0xFFFF_0000_0000_0000, 1)).IsNaN);
        Assert.Equal(Hex(0x7FFF_8000_0000_0000, 0), Real.NaN.ToBinary128Bits());
        Assert.True(Real.FromBinary128Bits(PosInf).IsPositiveInfinity);
        Assert.True(Real.FromBinary128Bits(PosInf | (UInt128.One << 127)).IsNegativeInfinity);
        Real negZero = Real.FromBinary128Bits(UInt128.One << 127);
        Assert.True(negZero.IsZero && negZero.SignBit);
    }

    [Fact]
    public void DefaultPrecisionIsExact()
    {
        int saved = Real.DefaultPrecision;
        try
        {
            Real.DefaultPrecision = 53;
            Real x = Real.FromBinary128Bits(Hex(0x3FFF_0000_0000_0000, 1));
            Assert.Equal(113, x.Precision);
            Assert.Equal(Hex(0x3FFF_0000_0000_0000, 1), x.ToBinary128Bits());
            Real r = Real.FromBinary128Bits(Hex(0x3FFF_0000_0000_0000, 1), 53, RealRounding.TowardPositive);
            Assert.Equal(53, r.Precision);
            Assert.Equal(0, Real.FromBigInt((BigInt.One << 52) + 1, -52, 53).CompareTo(r));
        }
        finally
        {
            Real.DefaultPrecision = saved;
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
            Real x = Real.FromBinary128Bits(bits);
            Assert.Equal(bits, x.ToBinary128Bits());
            Assert.Equal(bits, Real.FromBinary128Bits(bits, 500).ToBinary128Bits());
        }
    }

    /// <summary>Oracle: round |x| to a multiple of 2^max(BinaryExponent - 113, -16494), then check the range.</summary>
    private static Real Expected(Real x, RealRounding rnd, out bool overflow)
    {
        long ulp = Math.Max(x.BinaryExponent - 113, -16494);
        BigInt k = x.ShiftLeft(-ulp).ToBigInt(rnd);
        overflow = false;
        if (k.IsZero)
        {
            return x.SignBit ? Real.NegativeZero : Real.Zero;
        }
        Real r = Real.FromBigInt(k, ulp, (int)Math.Max(2, k.Abs().GetBitLength()));
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
                Real x = Real.FromBigInt(neg ? -m : m, baseExp + rng.Next(-120, 120) - prec + 113, prec);
                foreach (RealRounding rnd in Modes)
                {
                    UInt128 bits = x.ToBinary128Bits(rnd);
                    Real exp = Expected(x, rnd, out bool overflow);
                    if (overflow)
                    {
                        Assert.True(bits == (PosInf | (neg ? UInt128.One << 127 : 0)) || bits == (MaxFinite | (neg ? UInt128.One << 127 : 0)), $"{x} {rnd}");
                        bool toInf = rnd switch
                        {
                            RealRounding.TowardZero => false,
                            RealRounding.TowardPositive => !neg,
                            RealRounding.TowardNegative => neg,
                            _ => true,
                        };
                        Assert.Equal(toInf, (bits & ~(UInt128.One << 127)) == PosInf);
                        continue;
                    }
                    Real got = Real.FromBinary128Bits(bits);
                    Assert.True(exp.CompareTo(got) == 0 && exp.SignBit == got.SignBit, $"x={x} (prec {prec}) {rnd}: expected {exp}, got {got}");
                }
            }
        }
    }

    [Fact]
    public void TiesAtSubnormalBoundary()
    {
        Real half = MinSubnormal.ShiftRight(1);
        Assert.Equal(UInt128.Zero, half.ToBinary128Bits());
        Assert.Equal(UInt128.One, half.ToBinary128Bits(RealRounding.AwayFromZero));
        Assert.Equal(UInt128.One << 127, (-half).ToBinary128Bits(RealRounding.TowardPositive));
        Real threeHalves = MinSubnormal * Real.FromDouble(1.5);
        Assert.Equal((UInt128)2, threeHalves.ToBinary128Bits());
        Real big = Real.FromBigInt(BigInt.One, 16384, 2);
        Assert.Equal(PosInf, big.ToBinary128Bits());
        Assert.Equal(MaxFinite, big.ToBinary128Bits(RealRounding.TowardZero));
    }
}
