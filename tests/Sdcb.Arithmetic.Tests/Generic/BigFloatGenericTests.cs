using System.Globalization;
using System.Numerics;
using System.Text;
using Xunit;

namespace Sdcb.Arithmetic.Tests.Generic;

public class BigFloatGenericTests
{
    private static BigFloat R(double v, int prec = 53) => BigFloat.FromDouble(v, prec);

    /// <summary>The exact decimal expansion of man * 2^exp (always terminating).</summary>
    private static string ExactDecimal(BigInteger man, long exp)
    {
        if (exp >= 0)
        {
            return (man << (int)exp).ToString(CultureInfo.InvariantCulture);
        }
        int k = (int)-exp;
        BigInteger m = man * BigInteger.Pow(5, k);
        string s = BigInteger.Abs(m).ToString(CultureInfo.InvariantCulture).PadLeft(k + 1, '0');
        int dot = s.Length - k;
        return (m.Sign < 0 ? "-" : "") + s[..dot] + "." + s[dot..];
    }

    [Fact]
    public void BinaryFloatConversions_MatchCorrectlyRoundedParse()
    {
        var rng = new Random(30);
        (int lo, int hi, int fmt)[] ranges = [(-1140, 1030, 64), (-160, 130, 32), (-30, 20, 16)];
        for (int i = 0; i < 1500; i++)
        {
            (int lo, int hi, int fmt) = ranges[i % 3];
            int prec = rng.Next(2, 90);
            BigInteger man = new BigInteger(rng.NextInt64(1, long.MaxValue)) << 40 | rng.NextInt64(0, 1L << 40);
            if (rng.Next(2) == 0)
            {
                man = -man;
            }
            BigFloat x = BigFloat.FromBigInt((BigInt)man, rng.Next(lo, hi) - 104, prec);
            BigRational q = x.ToRational();
            string exact = ExactDecimal((BigInteger)q.Numerator, -(long)BigInt.Log2(q.Denominator));
            switch (fmt)
            {
                case 64:
                    double d = double.Parse(exact, CultureInfo.InvariantCulture);
                    Assert.Equal(d, x.ToDouble());
                    Assert.Equal(d, double.CreateChecked(x));
                    Assert.Equal(d, (double)q);
                    break;
                case 32:
                    float f = float.Parse(exact, CultureInfo.InvariantCulture);
                    Assert.Equal(f, x.ToSingle());
                    Assert.Equal(f, (float)x);
                    Assert.Equal(f, float.CreateChecked(x));
                    Assert.Equal(f, (float)q);
                    break;
                default:
                    Half h = Half.Parse(exact, CultureInfo.InvariantCulture);
                    Assert.Equal(h, x.ToHalf());
                    Assert.Equal(h, (Half)x);
                    Assert.Equal(h, Half.CreateChecked(x));
                    Assert.Equal(h, (Half)q);
                    break;
            }
        }
    }

    [Fact]
    public void ToDouble_DirectedRoundingAtTheEdges()
    {
        BigFloat big = BigFloat.FromBigInt(BigInt.One, 1024, 60); // 2^1024 overflows
        Assert.Equal(double.PositiveInfinity, big.ToDouble());
        Assert.Equal(double.MaxValue, big.ToDouble(RoundingMode.TowardZero));
        BigFloat justAbove = BigFloat.BitIncrement(BigFloat.FromDouble(double.MaxValue, 80));
        Assert.Equal(double.MaxValue, justAbove.ToDouble(RoundingMode.TowardZero));
        Assert.Equal(double.PositiveInfinity, justAbove.ToDouble(RoundingMode.AwayFromZero));
        BigFloat halfMin = BigFloat.FromBigInt(BigInt.One, -1075, 10);
        Assert.Equal(0.0, halfMin.ToDouble());
        Assert.Equal(double.Epsilon, halfMin.ToDouble(RoundingMode.TowardPositive));
        Assert.Equal(double.Epsilon, BigFloat.FromBigInt(new BigInt(3), -1076, 10).ToDouble()); // 0.75 ulp -> 1
        Assert.Equal(-0.0, (-halfMin).ToDouble(RoundingMode.TowardPositive));
        Assert.True(double.IsNegative((-halfMin).ToDouble()));
    }

    [Fact]
    public void ImplicitConversions_AreExact()
    {
        // DefaultPrecision is global and tests run in parallel, so only check the precision floor.
        BigFloat d = 0.1;
        Assert.Equal(0.1, d.ToDouble());
        Assert.True(d.Precision >= 53);
        BigFloat l = long.MaxValue;
        Assert.Equal(long.MaxValue, (long)l);
        BigFloat u = UInt128.MaxValue;
        Assert.Equal(UInt128.MaxValue, (UInt128)u);
        Assert.True(u.Precision >= 128);
        BigFloat f = 1.1f;
        Assert.Equal(1.1f, (float)f);
        BigFloat h = (Half)0.3;
        Assert.Equal((Half)0.3, (Half)h);
        BigFloat bi = BigInteger.Pow(3, 300);
        Assert.Equal(BigInteger.Pow(3, 300), (BigInteger)bi);
    }

    [Fact]
    public void IntegerConversions_FollowBclRules()
    {
        Assert.Equal(-3, (int)R(-3.99));
        Assert.Equal(-3, int.CreateChecked(R(-3.99)));
        Assert.Equal(0, int.CreateSaturating(BigFloat.NaN));
        Assert.Equal(int.MaxValue, int.CreateSaturating(BigFloat.PositiveInfinity));
        Assert.Equal(int.MinValue, int.CreateSaturating(BigFloat.NegativeInfinity));
        Assert.Equal(int.MaxValue, int.CreateTruncating(BigFloat.PositiveInfinity));
        Assert.Throws<OverflowException>(() => int.CreateChecked(BigFloat.NaN));
        Assert.Throws<OverflowException>(() => (int)BigFloat.PositiveInfinity);
        Assert.Throws<OverflowException>(() => (byte)R(256));
        Assert.Equal((byte)255, byte.CreateSaturating(R(1e300)));
        Assert.Equal((byte)0, byte.CreateSaturating(R(-1e300)));
        Assert.Equal(long.MaxValue, long.CreateSaturating(BigFloat.FromBigInt(BigInt.One, 1L << 40, 10)));
        Assert.Equal(0L, long.CreateTruncating(BigFloat.FromBigInt(BigInt.One, 1L << 40, 10)));
        Assert.Throws<OverflowException>(() => long.CreateChecked(BigFloat.FromBigInt(BigInt.One, 1L << 40, 10)));
        // Truncating keeps the low two's-complement bits.
        BigFloat x = BigFloat.FromBigInt((BigInt.One << 200) + 0x1234, 0, 220);
        Assert.Equal((ushort)0x1234, ushort.CreateTruncating(x));
        Assert.Equal(unchecked((int)0xFFFFEDCC), int.CreateTruncating(-x));

        foreach (long v in new long[] { 0, 1, -1, 1234567, long.MinValue, long.MaxValue })
        {
            BigFloat r = v;
            Assert.Equal(v, long.CreateChecked(r));
            Assert.Equal((BigInt)v, BigInt.CreateChecked(r));
            Assert.Equal((BigInteger)v, BigInteger.CreateChecked(r));
            Assert.True(r == BigFloat.CreateChecked(v));
        }
        Assert.Equal(BigInt.Zero, BigInt.CreateSaturating(BigFloat.NaN));
        Assert.Throws<OverflowException>(() => BigInt.CreateChecked(BigFloat.PositiveInfinity));
    }

    [Fact]
    public void DecimalConversions()
    {
        Assert.Equal(0.1, BigFloat.FromDecimal(0.1m, 53).ToDouble());
        Assert.Equal(0.5m, (decimal)R(0.5));
        Assert.Equal(-1234.25m, decimal.CreateChecked(R(-1234.25)));
        Assert.Equal(0.1000000000000000055511151231m, (decimal)R(0.1));
        Assert.Equal(decimal.MaxValue, (decimal)(BigFloat)decimal.MaxValue);
        Assert.Throws<OverflowException>(() => (decimal)BigFloat.NaN);
        Assert.Throws<OverflowException>(() => (decimal)R(1e30));
        Assert.Equal(decimal.MaxValue, decimal.CreateSaturating(R(1e30)));
        Assert.Equal(0m, decimal.CreateSaturating(BigFloat.NaN));
        Assert.Equal(0m, (decimal)BigFloat.FromBigInt(BigInt.One, -500, 10));
        Assert.True(R(2.5) == BigFloat.CreateChecked(2.5m));
    }

    [Fact]
    public void CrossTypeConversions()
    {
        BigRational third = new(1, 3);
        BigFloat r = (BigFloat)third;
        Assert.Equal(BigFloat.DefaultPrecision, r.Precision);
        Assert.True(BigFloat.Divide(BigFloat.One, 3) == r);
        Assert.Equal(new BigRational(5, 4), (BigRational)R(1.25));
        Assert.Equal(new BigRational(5, 4), BigRational.CreateChecked(R(1.25)));
        Assert.True(R(1.25) == BigFloat.CreateChecked(new BigRational(5, 4)));
        Assert.Equal(new BigInt(-2), (BigInt)R(-2.75));
        Assert.Throws<OverflowException>(() => (BigRational)BigFloat.NaN);
    }

    [Fact]
    public void BitIncrementDecrement_ScaleB_ILogB()
    {
        BigFloat x = R(1.0);
        Assert.Equal(1.0 + Math.Pow(2, -52), BigFloat.BitIncrement(x).ToDouble());
        Assert.Equal(1.0 - Math.Pow(2, -53), BigFloat.BitDecrement(x).ToDouble());
        Assert.True(BigFloat.BitDecrement(BigFloat.BitIncrement(x)) == x);
        Assert.True(BigFloat.BitIncrement(BigFloat.Zero) == BigFloat.Epsilon);
        Assert.True(BigFloat.Epsilon > BigFloat.Zero);
        Assert.True(BigFloat.BitDecrement(BigFloat.Epsilon) < BigFloat.Epsilon);
        Assert.True(BigFloat.BitIncrement(BigFloat.PositiveInfinity).IsPositiveInfinity);
        Assert.True(BigFloat.BitIncrement(BigFloat.NaN).IsNaN);
        Assert.Equal(12.0, BigFloat.ScaleB(R(3), 2).ToDouble());
        foreach (double v in new[] { 1.0, 0.75, 3e100, -5e-300, double.Epsilon, 0.0, double.NaN, double.PositiveInfinity })
        {
            Assert.Equal(double.ILogB(v), BigFloat.ILogB(R(v)));
        }
    }

    [Fact]
    public void Rounding_ModesAndDigits()
    {
        Assert.Equal(2.0, BigFloat.Round(R(2.5)).ToDouble());
        Assert.Equal(3.0, BigFloat.Round(R(2.5), MidpointRounding.AwayFromZero).ToDouble());
        Assert.Equal(-3.0, BigFloat.Round(R(-2.5), MidpointRounding.AwayFromZero).ToDouble());
        Assert.Equal(-2.0, BigFloat.Round(R(-2.5), MidpointRounding.ToZero).ToDouble());
        Assert.Equal(-3.0, BigFloat.Round(R(-2.5), MidpointRounding.ToNegativeInfinity).ToDouble());
        Assert.Equal(-2.0, BigFloat.Round(R(-2.5), MidpointRounding.ToPositiveInfinity).ToDouble());
        Assert.Equal(1.2, BigFloat.Round(R(1.25), 1).ToDouble());
        Assert.Equal(1.3, BigFloat.Round(R(1.25), 1, MidpointRounding.AwayFromZero).ToDouble());
        Assert.Equal(-1.3, BigFloat.Round(R(-1.21), 1, MidpointRounding.ToNegativeInfinity).ToDouble());
        Assert.Equal(2.67, BigFloat.Round(R(2.675), 2).ToDouble()); // 2.675 is really 2.67499999...
        Assert.Equal(0.001, BigFloat.Round(BigFloat.FromBigInt(BigInt.One, -1000, 53), 3, MidpointRounding.ToPositiveInfinity).ToDouble());
        Assert.True(BigFloat.Round(-BigFloat.FromBigInt(BigInt.One, -1000, 10), 3).SignBit);
        Assert.Equal(-2.0, BigFloat.Floor(R(-1.5)).ToDouble());
        Assert.Equal(-1.0, BigFloat.Ceiling(R(-1.5)).ToDouble());
        Assert.Equal(-1.0, BigFloat.Truncate(R(-1.5)).ToDouble());

        foreach (double v in new[] { 0.5, 1.5, 2.5, -0.5, -1.5, 3.7, -3.7, 1e20 + 0.0 })
        {
            foreach (MidpointRounding m in Enum.GetValues<MidpointRounding>())
            {
                Assert.Equal(Math.Round(v, m), BigFloat.Round(R(v), m).ToDouble());
            }
        }
    }

    [Fact]
    public void Predicates()
    {
        Assert.True(BigFloat.IsEvenInteger(R(4)));
        Assert.True(BigFloat.IsEvenInteger(R(0)));
        Assert.False(BigFloat.IsEvenInteger(R(4.5)));
        Assert.True(BigFloat.IsOddInteger(R(-3)));
        Assert.True(BigFloat.IsOddInteger(R(1)));
        Assert.False(BigFloat.IsOddInteger(R(1e20)));
        Assert.True(BigFloat.IsPow2(R(0.25)));
        Assert.False(BigFloat.IsPow2(R(-0.25)));
        Assert.True(BigFloat.IsNegative(BigFloat.NegativeZero));
        Assert.True(BigFloat.IsNormal(R(1e-300)));
        Assert.False(BigFloat.IsNormal(BigFloat.Zero));
        Assert.True(Check<BigFloat>(R(3)));
        Assert.True(Check<double>(3));

        static bool Check<T>(T three) where T : IFloatingPointIeee754<T> =>
            T.IsFinite(three) && T.IsInteger(three) && !T.IsZero(three) && T.IsPositive(three) &&
            T.IsNaN(T.NaN) && T.IsInfinity(T.NegativeInfinity) && T.IsNegative(T.NegativeZero) &&
            T.Sign(-three) == -1 && T.IsRealNumber(three) && !T.IsSubnormal(three) && T.IsOddInteger(three);
    }

    [Fact]
    public void MinMax_GenericSemanticsMatchDouble()
    {
        double[] vals = [1, -1, 0.0, -0.0, double.NaN, 2, -2, double.PositiveInfinity];
        foreach (double a in vals)
        {
            foreach (double b in vals)
            {
                Same(double.Max(a, b), Max(R(a), R(b)));
                Same(double.Min(a, b), Min(R(a), R(b)));
                Same(double.MaxNumber(a, b), BigFloat.MaxNumber(R(a), R(b)));
                Same(double.MinNumber(a, b), BigFloat.MinNumber(R(a), R(b)));
                Same(double.MaxMagnitude(a, b), BigFloat.MaxMagnitude(R(a), R(b)));
                Same(double.MinMagnitude(a, b), BigFloat.MinMagnitude(R(a), R(b)));
                Same(double.MaxMagnitudeNumber(a, b), BigFloat.MaxMagnitudeNumber(R(a), R(b)));
                Same(double.MinMagnitudeNumber(a, b), BigFloat.MinMagnitudeNumber(R(a), R(b)));
            }
        }

        static T Max<T>(T a, T b) where T : INumber<T> => T.Max(a, b);
        static T Min<T>(T a, T b) where T : INumber<T> => T.Min(a, b);
        static void Same(double expected, BigFloat actual) =>
            Assert.Equal(BitConverter.DoubleToInt64Bits(expected), BitConverter.DoubleToInt64Bits(actual.ToDouble()));
    }

    [Fact]
    public void Bitwise_OnExactBinaryExpansions()
    {
        Assert.Equal(0.5, And(R(0.75), R(0.5)).ToDouble());
        Assert.Equal(3.5, Or(R(3), R(0.5)).ToDouble());
        Assert.Equal(5.0, Xor(R(6), R(3)).ToDouble());
        Assert.Equal(0.0, And(R(4), R(0.5)).ToDouble());
        Assert.Equal(4.5, Xor(R(4), R(0.5)).ToDouble());
        Assert.Equal(-6.0, Or(R(-6), R(2)).ToDouble());
        Assert.Equal(6.0, Xor(R(-4), R(-2)).ToDouble());
        Assert.True(And(BigFloat.NaN, R(1)).IsNaN);
        // ~ complements the 2-bit significand 11 -> 00 (zero) and 10 -> 01.
        Assert.True(Not(R(3, 2)).IsZero);
        Assert.Equal(-1.0, Not(R(2, 2)).ToDouble());
        Assert.True(AllBits<BigFloat>().IsNaN);

        static T And<T>(T a, T b) where T : IBitwiseOperators<T, T, T> => a & b;
        static T Or<T>(T a, T b) where T : IBitwiseOperators<T, T, T> => a | b;
        static T Xor<T>(T a, T b) where T : IBitwiseOperators<T, T, T> => a ^ b;
        static T Not<T>(T a) where T : IBitwiseOperators<T, T, T> => ~a;
        static T AllBits<T>() where T : IBinaryNumber<T> => T.AllBitsSet;
    }

    private static void AssertUlps(double expected, BigFloat actual, int ulps = 1)
    {
        double a = actual.ToDouble();
        if (expected == a)
        {
            return;
        }
        long d = Math.Abs(BitConverter.DoubleToInt64Bits(expected) - BitConverter.DoubleToInt64Bits(a));
        Assert.True(d <= ulps, $"expected {expected:R}, got {a:R}");
    }

    [Fact]
    public void PiFunctions_AndDegrees()
    {
        var rng = new Random(31);
        for (int i = 0; i < 200; i++)
        {
            double x = rng.NextDouble() * 2 - 1;
            AssertUlps(double.AsinPi(x), BigFloat.AsinPi(R(x)), 2);
            AssertUlps(double.AcosPi(x), BigFloat.AcosPi(R(x)), 2);
            double t = (rng.NextDouble() - 0.5) * 100;
            AssertUlps(double.AtanPi(t), BigFloat.AtanPi(R(t)), 2);
            AssertUlps(double.Atan2Pi(x, t), BigFloat.Atan2Pi(R(x), R(t)), 2);
            AssertUlps(double.DegreesToRadians(t), BigFloat.DegreesToRadians(R(t)), 3);
            AssertUlps(double.RadiansToDegrees(x), BigFloat.RadiansToDegrees(R(x)), 3);
        }
        Assert.Equal(0.5, BigFloat.AsinPi(R(1)).ToDouble());
        Assert.Equal(-0.5, BigFloat.AsinPi(R(-1)).ToDouble());
        Assert.Equal(1.0, BigFloat.AcosPi(R(-1)).ToDouble());
        Assert.Equal(0.5, BigFloat.AcosPi(R(0)).ToDouble());
        Assert.True(BigFloat.AcosPi(R(1)).IsZero);
        Assert.Equal(0.25, BigFloat.AtanPi(R(1)).ToDouble());
        Assert.Equal(-0.5, BigFloat.AtanPi(BigFloat.NegativeInfinity).ToDouble());
        Assert.Equal(-0.75, BigFloat.Atan2Pi(R(-2), R(-2)).ToDouble());
        Assert.Equal(1.0, BigFloat.Atan2Pi(R(0), R(-3)).ToDouble());
        Assert.Equal(-1.0, BigFloat.Atan2Pi(R(-0.0), R(-0.0)).ToDouble());
        Assert.Equal(0.75, BigFloat.Atan2Pi(BigFloat.PositiveInfinity, BigFloat.NegativeInfinity).ToDouble());
        Assert.True(BigFloat.AsinPi(R(2)).IsNaN);
        Assert.True(BigFloat.DegreesToRadians(BigFloat.FromBigInt(180, 0, 128)) == BigFloat.Pi(128));
        (BigFloat s, BigFloat c) = BigFloat.SinCosPi(R(0.25));
        Assert.True(s == c);
    }

    [Fact]
    public void LogWithBase()
    {
        Assert.Equal(3.0, BigFloat.Log(R(8), R(2)).ToDouble());
        Assert.Equal(-3.0, BigFloat.Log(R(0.125), R(2)).ToDouble());
        Assert.Equal(1.5, BigFloat.Log(R(8), R(4)).ToDouble());
        Assert.Equal(0.5, BigFloat.Log(R(3), R(9)).ToDouble());
        Assert.Equal(2.0, BigFloat.Log(R(100), R(10)).ToDouble());
        Assert.Equal(5.0, BigFloat.Log(R(243), R(3)).ToDouble());
        Assert.True(BigFloat.Log(R(1), R(7)).IsZero);
        Assert.True(BigFloat.Log(R(5), R(1)).IsNaN);
        Assert.True(BigFloat.Log(R(-5), R(2)).IsNaN);
        Assert.True(BigFloat.Log(R(0), R(2)).IsNegativeInfinity);
        var rng = new Random(32);
        for (int i = 0; i < 100; i++)
        {
            double x = rng.NextDouble() * 1000, b = rng.NextDouble() * 20 + 0.01;
            AssertUlps(Math.Log(x) / Math.Log(b), BigFloat.Log(R(x), R(b)), 4);
        }
    }

    [Fact]
    public void FunctionOverloads_MatchPrecisionVersions()
    {
        BigFloat x = BigFloat.FromDouble(0.7, 100);
        Assert.Equal(BigFloat.Exp(x, 0), BigFloat.Exp(x));
        Assert.Equal(BigFloat.ExpM1(x, 0), Exp<BigFloat>(x));
        Assert.Equal(BigFloat.Exp2M1(x, 0), BigFloat.Exp2M1(x));
        Assert.Equal(BigFloat.Exp10M1(x, 0), BigFloat.Exp10M1(x));
        Assert.Equal(BigFloat.LogP1(x, 0), BigFloat.LogP1(x));
        Assert.Equal(BigFloat.Log2P1(x, 0), BigFloat.Log2P1(x));
        Assert.Equal(BigFloat.Log10P1(x, 0), BigFloat.Log10P1(x));
        Assert.Equal(BigFloat.Sin(x, 0), BigFloat.Sin(x));
        Assert.Equal(BigFloat.Atanh(x, 0), BigFloat.Atanh(x));
        Assert.Equal(BigFloat.RootN(x, 5L, 0), BigFloat.RootN(x, 5));
        Assert.Equal(BigFloat.Pow(x, x, 0), BigFloat.Pow(x, x));
        Assert.Equal(BigFloat.Hypot(x, x, 0), BigFloat.Hypot(x, x));
        Assert.Equal(100, BigFloat.Exp(x).Precision);
        Assert.Equal(BigFloat.Pi(0).ShiftLeft(1), Tau<BigFloat>());
        Assert.Equal(BigFloat.E(0), E<BigFloat>());
        Assert.True(BigFloat.FusedMultiplyAdd(x, x, x) == BigFloat.FusedMultiplyAdd(x, x, x));
        Assert.True(BigFloat.Ieee754Remainder(R(7), R(2)) == R(-1));
        Assert.True(R(-7) % R(2) == R(-1));
        BigFloat y = R(1.5);
        y++;
        Assert.Equal(2.5, y.ToDouble());
        Assert.Equal(53, y.Precision);

        static T Exp<T>(T v) where T : IExponentialFunctions<T> => T.ExpM1(v);
        static T Tau<T>() where T : IFloatingPointConstants<T> => T.Tau;
        static T E<T>() where T : IFloatingPointConstants<T> => T.E;
    }

    [Fact]
    public void Layout_ExponentAndSignificand()
    {
        IFloatingPoint<BigFloat> one = R(1);
        Assert.Equal(0, one.GetExponentShortestBitLength());
        Assert.Equal(53, one.GetSignificandBitLength());
        Assert.Equal(7, one.GetSignificandByteCount());
        byte[] buf = new byte[8];
        Assert.True(((IFloatingPoint<BigFloat>)R(0.75)).TryWriteExponentLittleEndian(buf, out int n));
        Assert.Equal(8, n);
        Assert.Equal(-1L, BitConverter.ToInt64(buf));
        Assert.True(((IFloatingPoint<BigFloat>)R(3, 8)).TryWriteSignificandLittleEndian(buf, out n));
        Assert.Equal(1, n);
        Assert.Equal(0xC0, buf[0]);
    }

    [Theory]
    [InlineData("1,234.5", NumberStyles.Number, 1234.5)]
    [InlineData("(1.5)", NumberStyles.Float | NumberStyles.AllowParentheses, -1.5)]
    [InlineData(" -2.5e-3 ", NumberStyles.Float, -0.0025)]
    [InlineData("Infinity", NumberStyles.Float, double.PositiveInfinity)]
    [InlineData("-Infinity", NumberStyles.Float, double.NegativeInfinity)]
    [InlineData("NaN", NumberStyles.Float, double.NaN)]
    [InlineData("-inf", NumberStyles.Float, double.NegativeInfinity)]
    public void ParseStyles(string s, NumberStyles style, double expected)
    {
        BigFloat r = BigFloat.Parse(s, style, CultureInfo.InvariantCulture);
        Assert.Equal(expected, r.ToDouble());
        Assert.Equal(expected, BigFloat.Parse(Encoding.UTF8.GetBytes(s), style, CultureInfo.InvariantCulture).ToDouble());
        Assert.Equal(expected, double.Parse(s.Replace("inf", "Infinity"), style, CultureInfo.InvariantCulture));
    }

    [Fact]
    public void ParseStyles_Rejects()
    {
        Assert.False(BigFloat.TryParse("1,5", NumberStyles.Float, CultureInfo.InvariantCulture, out _));
        Assert.False(BigFloat.TryParse("1.5", NumberStyles.Integer, CultureInfo.InvariantCulture, out _));
        Assert.Throws<ArgumentException>(() => BigFloat.TryParse("ff", NumberStyles.HexNumber, null, out _));
        Assert.Equal(1.5, BigFloat.Parse("1,5", NumberStyles.Float, CultureInfo.GetCultureInfo("de-DE")).ToDouble());
    }

    [Fact]
    public void Utf8()
    {
        BigFloat x = R(0.375);
        byte[] buf = new byte[64];
        Assert.True(x.TryFormat(buf, out int n, default, null));
        Assert.Equal(x.ToString(), Encoding.UTF8.GetString(buf, 0, n));
        Assert.True(BigFloat.Parse("0.375"u8, null) == x);
    }
}
