using System.Globalization;
using System.Numerics;
using System.Text;
using Xunit;

namespace Sdcb.Arithmetic.Tests.Generic;

public class BigIntGenericTests
{
    private static IEnumerable<(BigInt A, BigInteger B)> Samples(int seed, int count = 400, int maxLimbs = 5)
    {
        var rng = new Random(seed);
        long[] edges = [0, 1, -1, 2, -2, 255, 256, -255, -256, long.MaxValue, long.MinValue, int.MinValue, uint.MaxValue];
        foreach (long e in edges)
        {
            yield return (e, e);
        }
        BigInteger p64 = BigInteger.One << 64, p63 = BigInteger.One << 63, p128 = BigInteger.One << 128;
        foreach (BigInteger e in new[] { p64, -p64, p64 - 1, -(p64 - 1), p63, -p63, p128, -p128, p128 - 1, -(p128 + 1) })
        {
            yield return ((BigInt)e, e);
        }
        for (int i = 0; i < count; i++)
        {
            yield return TestHelpers.RandomPair(rng, maxLimbs);
        }
    }

    private static int WidthLimbs(BigInteger b)
    {
        long bits = b.Sign >= 0 ? (long)b.GetBitLength() : ((IBinaryInteger<BigInteger>)b).GetShortestBitLength();
        return (int)((bits + 63) / 64);
    }

    private static BigInteger Unsigned(BigInteger b, int width) => ((b % (BigInteger.One << width)) + (BigInteger.One << width)) % (BigInteger.One << width);

    [Fact]
    public void BitCounts_MatchLimbWidthDefinition()
    {
        foreach ((BigInt a, BigInteger b) in Samples(10))
        {
            string ctx = b.ToString();
            // trailing zeros: same as BigInteger except for zero (one limb width).
            Assert.Equal(b.IsZero ? 64 : (BigInteger)BigInteger.TrailingZeroCount(b), (BigInteger)BigInt.TrailingZeroCount(a));

            long bits = (long)b.GetBitLength();
            BigInteger lzc = b.Sign < 0 ? 0 : b.IsZero ? 64 : 64 * ((bits + 63) / 64) - bits;
            Assert.True(lzc == (BigInteger)BigInt.LeadingZeroCount(a), ctx);

            int w = 64 * Math.Max(1, WidthLimbs(b));
            BigInteger u = b.Sign >= 0 ? b : Unsigned(b, w);
            long pop = 0;
            foreach (byte by in u.ToByteArray(isUnsigned: true))
            {
                pop += BitOperations.PopCount(by);
            }
            Assert.True(pop == (BigInteger)BigInt.PopCount(a), ctx);

            IBinaryInteger<BigInteger> bi = b;
            Assert.Equal(bi.GetShortestBitLength(), a.GetShortestBitLength());
            if (b.Sign > 0)
            {
                Assert.Equal(BigInteger.Log2(b), (BigInteger)BigInt.Log2(a));
            }
            Assert.Equal(BigInteger.IsPow2(b), BigInt.IsPow2(a));
            Assert.Equal(b.IsEven, BigInt.IsEvenInteger(a));
            Assert.Equal(!b.IsEven, BigInt.IsOddInteger(a));
        }
        Assert.Throws<ArgumentOutOfRangeException>(() => BigInt.Log2(-5));
        Assert.Equal(BigInt.Zero, BigInt.Log2(0));
    }

    [Fact]
    public void RotateAndUnsignedShift_MatchTwosComplementDefinition()
    {
        var rng = new Random(11);
        foreach ((BigInt a, BigInteger b) in Samples(12, 200))
        {
            if (b.IsZero)
            {
                Assert.Equal(BigInt.Zero, BigInt.RotateLeft(a, 5));
                continue;
            }
            int width = 64 * WidthLimbs(b);
            BigInteger u = Unsigned(b, width);
            foreach (int k in new[] { 1, 7, 63, 64, 65, width - 1, rng.Next(0, 3 * width) })
            {
                int r = k % width;
                BigInteger rot = ((u << r) | (u >> (width - r))) % (BigInteger.One << width);
                if (b.Sign < 0 && rot >= BigInteger.One << (width - 1))
                {
                    rot -= BigInteger.One << width;
                }
                string ctx = $"{b} rotl {k}";
                Assert.True(rot == (BigInteger)BigInt.RotateLeft(a, k), ctx);
                Assert.True(rot == (BigInteger)BigInt.RotateRight(a, -k), ctx);

                BigInteger urs = b.Sign >= 0 ? b >> k : k >= width ? 0 : u >> k;
                Assert.True(urs == (BigInteger)(a >>> k), $"{b} >>> {k}");
            }
        }
    }

    [Fact]
    public void Bytes_MatchBigIntegerToByteArray()
    {
        foreach ((BigInt a, BigInteger b) in Samples(13))
        {
            byte[] le = b.ToByteArray();
            Assert.Equal(le.Length, a.GetByteCount());
            byte[] buf = new byte[le.Length + 3];
            Assert.True(a.TryWriteLittleEndian(buf, out int n));
            Assert.Equal(le, buf.AsSpan(0, n).ToArray());
            Assert.False(a.TryWriteLittleEndian(new byte[le.Length - 1], out _));

            byte[] be = b.ToByteArray(isBigEndian: true);
            Assert.True(a.TryWriteBigEndian(buf, out n));
            Assert.Equal(be, buf.AsSpan(0, n).ToArray());

            Assert.Equal(a, BigInt.ReadLittleEndian(le, isUnsigned: false));
            Assert.Equal(a, BigInt.ReadBigEndian(be, isUnsigned: false));
            Assert.Equal(new BigInteger(le, isUnsigned: true), (BigInteger)BigInt.ReadLittleEndian(le, isUnsigned: true));

            // Through the generic interface too.
            Assert.Equal(a, ReadGeneric<BigInt>(le));
        }
    }

    private static T ReadGeneric<T>(byte[] le) where T : IBinaryInteger<T> => T.ReadLittleEndian(le, isUnsigned: false);

    [Fact]
    public void DivRemAndMagnitudeHelpers()
    {
        var rng = new Random(14);
        for (int i = 0; i < 300; i++)
        {
            (BigInt a, BigInteger ab) = TestHelpers.RandomPair(rng, 6);
            (BigInt d, BigInteger db) = TestHelpers.RandomPair(rng, 3);
            if (db.IsZero)
            {
                continue;
            }
            (BigInt q, BigInt r) = BigInt.DivRem(a, d);
            (BigInteger qb, BigInteger rb) = BigInteger.DivRem(ab, db);
            Assert.Equal(qb, (BigInteger)q);
            Assert.Equal(rb, (BigInteger)r);
            Assert.Equal(BigInteger.MaxMagnitude(ab, db), (BigInteger)BigInt.MaxMagnitude(a, d));
            Assert.Equal(BigInteger.MinMagnitude(ab, db), (BigInteger)BigInt.MinMagnitude(a, d));
            Assert.Equal(BigInteger.CopySign(ab, db), (BigInteger)BigInt.CopySign(a, d));
        }
        Assert.Equal(new BigInt(3), BigInt.MaxMagnitude(-3, 3));
        Assert.Equal(new BigInt(-3), BigInt.MinMagnitude(-3, 3));
        Assert.Equal(new BigInt(3), BigInt.Clamp(5, 1, 3));
        Assert.Throws<ArgumentException>(() => BigInt.Clamp(5, 3, 1));
    }

    [Fact]
    public void Fused_MulAddMulSubMulMod()
    {
        var rng = new Random(15);
        for (int i = 0; i < 400; i++)
        {
            BigInt a = TestHelpers.RandomPair(rng, 6).A;
            BigInt b = TestHelpers.RandomPair(rng, 6).A;
            BigInt c = TestHelpers.RandomPair(rng, 14).A;
            Assert.Equal(a * b + c, BigInt.MultiplyAdd(a, b, c));
            Assert.Equal(a * b - c, BigInt.MultiplySubtract(a, b, c));
            Assert.Equal(a * b - a * b, BigInt.MultiplySubtract(a, b, a * b));
            if (!c.IsZero)
            {
                BigInt m = BigInt.ModMultiply(a, b, c);
                Assert.True(m.Sign >= 0 && m < c.Abs());
                Assert.Equal(BigInt.Mod(a * b, c.Abs()), m);
            }
        }
        Assert.Throws<DivideByZeroException>(() => BigInt.ModMultiply(3, 4, 0));
    }

    [Fact]
    public void BigIntegerInterop_RoundTrips()
    {
        foreach ((BigInt a, BigInteger b) in Samples(16))
        {
            Assert.Equal(b, (BigInteger)a);
            BigInt back = b;
            Assert.Equal(a, back);
            Assert.Equal(a, BigInt.CreateChecked(b));
            Assert.Equal(b, BigInteger.CreateChecked(a));
        }
    }

    [Fact]
    public void PrimitiveConversions_MatchBigInteger()
    {
        foreach ((BigInt a, BigInteger b) in Samples(17, 300, 3))
        {
            Check<sbyte>(a, b); Check<byte>(a, b); Check<short>(a, b); Check<ushort>(a, b); Check<char>(a, b);
            Check<int>(a, b); Check<uint>(a, b); Check<long>(a, b); Check<ulong>(a, b);
            Check<nint>(a, b); Check<nuint>(a, b); Check<Int128>(a, b); Check<UInt128>(a, b);

            string s = b.ToString(CultureInfo.InvariantCulture);
            Assert.Equal(double.Parse(s, CultureInfo.InvariantCulture), double.CreateChecked(a));
            Assert.Equal(double.Parse(s, CultureInfo.InvariantCulture), (double)a);
            Assert.Equal(float.Parse(s, CultureInfo.InvariantCulture), float.CreateChecked(a));
            Assert.Equal(float.Parse(s, CultureInfo.InvariantCulture), (float)a);
            Assert.Equal(Half.Parse(s, CultureInfo.InvariantCulture), Half.CreateChecked(a));
            if (BigInteger.Abs(b) < BigInteger.One << 96)
            {
                Assert.Equal((decimal)b, decimal.CreateChecked(a));
                Assert.Equal((decimal)b, (decimal)a);
                Assert.Equal(a, BigInt.CreateChecked((decimal)b));
            }
            else
            {
                Assert.Throws<OverflowException>(() => decimal.CreateChecked(a));
                Assert.Equal(decimal.CreateSaturating(b), decimal.CreateSaturating(a));
            }
        }

        static void Check<T>(BigInt a, BigInteger b) where T : INumberBase<T>
        {
            Assert.Equal(T.CreateTruncating(b), T.CreateTruncating(a));
            Assert.Equal(T.CreateSaturating(b), T.CreateSaturating(a));
            bool fits;
            T expected = default!;
            try
            {
                expected = T.CreateChecked(b);
                fits = true;
            }
            catch (OverflowException)
            {
                fits = false;
            }
            if (fits)
            {
                Assert.Equal(expected, T.CreateChecked(a));
                Assert.Equal(a, BigInt.CreateChecked(expected));
            }
            else
            {
                Assert.Throws<OverflowException>(() => T.CreateChecked(a));
            }
        }
    }

    [Fact]
    public void FloatConversions_FromFloatingPoint()
    {
        Assert.Equal(new BigInt(-1), BigInt.CreateTruncating(-1.9));
        Assert.Equal(new BigInt(2), BigInt.CreateChecked(2.9m));
        Assert.Equal(new BigInt(65504), BigInt.CreateChecked(Half.MaxValue));
        Assert.Equal((BigInt)(BigInteger)float.MaxValue, BigInt.CreateChecked(float.MaxValue));
        Assert.Equal(BigInt.Zero, BigInt.CreateSaturating(double.NaN));
        Assert.Throws<OverflowException>(() => BigInt.CreateChecked(double.NaN));
        Assert.Throws<OverflowException>(() => BigInt.CreateSaturating(double.PositiveInfinity));
        Assert.Equal(new BigInt(-7), (BigInt)(-7.99f));
        Assert.Equal(new BigInt(123), (BigInt)123.999m);
    }

    [Theory]
    [InlineData("FF", NumberStyles.HexNumber)]
    [InlineData("0FF", NumberStyles.HexNumber)]
    [InlineData("  7fffffffffffffffffff  ", NumberStyles.HexNumber)]
    [InlineData("8000000000000000000000", NumberStyles.AllowHexSpecifier)]
    [InlineData("101", NumberStyles.BinaryNumber)]
    [InlineData("0101", NumberStyles.BinaryNumber)]
    [InlineData("1,234,567,890,123,456,789", NumberStyles.Number)]
    [InlineData("  -12  ", NumberStyles.Integer)]
    [InlineData("(12)", NumberStyles.AllowParentheses)]
    [InlineData("1.00", NumberStyles.Number)]
    [InlineData("1.50", NumberStyles.Number)]
    [InlineData("1e3", NumberStyles.Float)]
    [InlineData("-12.5e1", NumberStyles.Float)]
    [InlineData("12-", NumberStyles.AllowTrailingSign)]
    [InlineData("12", NumberStyles.None)]
    [InlineData(" 12", NumberStyles.None)]
    [InlineData("+99999999999999999999999999", NumberStyles.Integer)]
    [InlineData("1,2", NumberStyles.Integer)]
    [InlineData("", NumberStyles.Integer)]
    [InlineData("12x", NumberStyles.Integer)]
    public void ParseStyles_MatchBigInteger(string s, NumberStyles style)
    {
        bool ok = BigInteger.TryParse(s, style, CultureInfo.InvariantCulture, out BigInteger expected);
        Assert.Equal(ok, BigInt.TryParse(s, style, CultureInfo.InvariantCulture, out BigInt actual));
        if (ok)
        {
            Assert.Equal(expected, (BigInteger)actual);
            Assert.Equal(actual, BigInt.Parse(Encoding.UTF8.GetBytes(s), style, CultureInfo.InvariantCulture));
        }
    }

    [Fact]
    public void Parse_CultureSymbols()
    {
        var de = CultureInfo.GetCultureInfo("de-DE");
        Assert.Equal(new BigInt(1234567), BigInt.Parse("1.234.567", NumberStyles.Number, de));
        Assert.Equal(new BigInt(12), BigInt.Parse("12,000", NumberStyles.Number, de));
        Assert.Throws<ArgumentException>(() => BigInt.Parse("12", NumberStyles.HexNumber | NumberStyles.AllowLeadingSign, null));
    }

    [Fact]
    public void Utf8_FormatAndParse()
    {
        BigInt v = BigInt.Parse("-123456789012345678901234567890");
        byte[] buf = new byte[64];
        Assert.True(v.TryFormat(buf, out int n, default, null));
        Assert.Equal("-123456789012345678901234567890", Encoding.UTF8.GetString(buf, 0, n));
        Assert.False(v.TryFormat(new byte[5], out _, default, null));
        Assert.True(v.TryFormat(buf, out n, "x", null));
        Assert.Equal(v.ToString("x", null), Encoding.UTF8.GetString(buf, 0, n));
        Assert.Equal(v, BigInt.Parse(Encoding.UTF8.GetBytes("-123456789012345678901234567890"), null));
        Assert.False(BigInt.TryParse("12a"u8, null, out _));
        Assert.Equal($"[{v}]", $"[{v}]"); // interpolation via ISpanFormattable
    }

    [Fact]
    public void GenericIdentitiesAndPredicates()
    {
        Assert.Equal(BigInt.Zero, Zero<BigInt>());
        Assert.Equal(BigInt.One, One<BigInt>());
        Assert.Equal(BigInt.MinusOne, BigInt.NegativeOne);
        Assert.True(IsNegative<BigInt>(-3));
        Assert.True(IsPositive<BigInt>(0));
        Assert.Equal(-1, Sign<BigInt>(-8));
        Assert.Equal(new BigInt(-1), AllBits<BigInt>());

        static T Zero<T>() where T : INumberBase<T> => T.AdditiveIdentity;
        static T One<T>() where T : INumberBase<T> => T.MultiplicativeIdentity;
        static bool IsNegative<T>(T v) where T : INumberBase<T> => T.IsNegative(v);
        static bool IsPositive<T>(T v) where T : INumberBase<T> => T.IsPositive(v);
        static int Sign<T>(T v) where T : INumber<T> => T.Sign(v);
        static T AllBits<T>() where T : IBinaryNumber<T> => T.AllBitsSet;
    }
}
