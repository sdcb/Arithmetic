using System.Globalization;
using System.Numerics;
using System.Text;
using Xunit;

namespace Sdcb.Arithmetic.Tests.Generic;

public class RealGenericTests
{
    private static Real R(double v, int prec = 53) => Real.FromDouble(v, prec);

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
            Real x = Real.FromBigInt((BigInt)man, rng.Next(lo, hi) - 104, prec);
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
        Real big = Real.FromBigInt(BigInt.One, 1024, 60); // 2^1024 overflows
        Assert.Equal(double.PositiveInfinity, big.ToDouble());
        Assert.Equal(double.MaxValue, big.ToDouble(RealRounding.TowardZero));
        Real justAbove = Real.BitIncrement(Real.FromDouble(double.MaxValue, 80));
        Assert.Equal(double.MaxValue, justAbove.ToDouble(RealRounding.TowardZero));
        Assert.Equal(double.PositiveInfinity, justAbove.ToDouble(RealRounding.AwayFromZero));
        Real halfMin = Real.FromBigInt(BigInt.One, -1075, 10);
        Assert.Equal(0.0, halfMin.ToDouble());
        Assert.Equal(double.Epsilon, halfMin.ToDouble(RealRounding.TowardPositive));
        Assert.Equal(double.Epsilon, Real.FromBigInt(new BigInt(3), -1076, 10).ToDouble()); // 0.75 ulp -> 1
        Assert.Equal(-0.0, (-halfMin).ToDouble(RealRounding.TowardPositive));
        Assert.True(double.IsNegative((-halfMin).ToDouble()));
    }

    [Fact]
    public void ImplicitConversions_AreExact()
    {
        // DefaultPrecision is global and tests run in parallel, so only check the precision floor.
        Real d = 0.1;
        Assert.Equal(0.1, d.ToDouble());
        Assert.True(d.Precision >= 53);
        Real l = long.MaxValue;
        Assert.Equal(long.MaxValue, (long)l);
        Real u = UInt128.MaxValue;
        Assert.Equal(UInt128.MaxValue, (UInt128)u);
        Assert.True(u.Precision >= 128);
        Real f = 1.1f;
        Assert.Equal(1.1f, (float)f);
        Real h = (Half)0.3;
        Assert.Equal((Half)0.3, (Half)h);
        Real bi = BigInteger.Pow(3, 300);
        Assert.Equal(BigInteger.Pow(3, 300), (BigInteger)bi);
    }

    [Fact]
    public void IntegerConversions_FollowBclRules()
    {
        Assert.Equal(-3, (int)R(-3.99));
        Assert.Equal(-3, int.CreateChecked(R(-3.99)));
        Assert.Equal(0, int.CreateSaturating(Real.NaN));
        Assert.Equal(int.MaxValue, int.CreateSaturating(Real.PositiveInfinity));
        Assert.Equal(int.MinValue, int.CreateSaturating(Real.NegativeInfinity));
        Assert.Equal(int.MaxValue, int.CreateTruncating(Real.PositiveInfinity));
        Assert.Throws<OverflowException>(() => int.CreateChecked(Real.NaN));
        Assert.Throws<OverflowException>(() => (int)Real.PositiveInfinity);
        Assert.Throws<OverflowException>(() => (byte)R(256));
        Assert.Equal((byte)255, byte.CreateSaturating(R(1e300)));
        Assert.Equal((byte)0, byte.CreateSaturating(R(-1e300)));
        Assert.Equal(long.MaxValue, long.CreateSaturating(Real.FromBigInt(BigInt.One, 1L << 40, 10)));
        Assert.Equal(0L, long.CreateTruncating(Real.FromBigInt(BigInt.One, 1L << 40, 10)));
        Assert.Throws<OverflowException>(() => long.CreateChecked(Real.FromBigInt(BigInt.One, 1L << 40, 10)));
        // Truncating keeps the low two's-complement bits.
        Real x = Real.FromBigInt((BigInt.One << 200) + 0x1234, 0, 220);
        Assert.Equal((ushort)0x1234, ushort.CreateTruncating(x));
        Assert.Equal(unchecked((int)0xFFFFEDCC), int.CreateTruncating(-x));

        foreach (long v in new long[] { 0, 1, -1, 1234567, long.MinValue, long.MaxValue })
        {
            Real r = v;
            Assert.Equal(v, long.CreateChecked(r));
            Assert.Equal((BigInt)v, BigInt.CreateChecked(r));
            Assert.Equal((BigInteger)v, BigInteger.CreateChecked(r));
            Assert.True(r == Real.CreateChecked(v));
        }
        Assert.Equal(BigInt.Zero, BigInt.CreateSaturating(Real.NaN));
        Assert.Throws<OverflowException>(() => BigInt.CreateChecked(Real.PositiveInfinity));
    }

    [Fact]
    public void DecimalConversions()
    {
        Assert.Equal(0.1, Real.FromDecimal(0.1m, 53).ToDouble());
        Assert.Equal(0.5m, (decimal)R(0.5));
        Assert.Equal(-1234.25m, decimal.CreateChecked(R(-1234.25)));
        Assert.Equal(0.1000000000000000055511151231m, (decimal)R(0.1));
        Assert.Equal(decimal.MaxValue, (decimal)(Real)decimal.MaxValue);
        Assert.Throws<OverflowException>(() => (decimal)Real.NaN);
        Assert.Throws<OverflowException>(() => (decimal)R(1e30));
        Assert.Equal(decimal.MaxValue, decimal.CreateSaturating(R(1e30)));
        Assert.Equal(0m, decimal.CreateSaturating(Real.NaN));
        Assert.Equal(0m, (decimal)Real.FromBigInt(BigInt.One, -500, 10));
        Assert.True(R(2.5) == Real.CreateChecked(2.5m));
    }

    [Fact]
    public void CrossTypeConversions()
    {
        BigRational third = new(1, 3);
        Real r = (Real)third;
        Assert.Equal(Real.DefaultPrecision, r.Precision);
        Assert.True(Real.Divide(Real.One, 3) == r);
        Assert.Equal(new BigRational(5, 4), (BigRational)R(1.25));
        Assert.Equal(new BigRational(5, 4), BigRational.CreateChecked(R(1.25)));
        Assert.True(R(1.25) == Real.CreateChecked(new BigRational(5, 4)));
        Assert.Equal(new BigInt(-2), (BigInt)R(-2.75));
        Assert.Throws<OverflowException>(() => (BigRational)Real.NaN);
    }

    [Fact]
    public void BitIncrementDecrement_ScaleB_ILogB()
    {
        Real x = R(1.0);
        Assert.Equal(1.0 + Math.Pow(2, -52), Real.BitIncrement(x).ToDouble());
        Assert.Equal(1.0 - Math.Pow(2, -53), Real.BitDecrement(x).ToDouble());
        Assert.True(Real.BitDecrement(Real.BitIncrement(x)) == x);
        Assert.True(Real.BitIncrement(Real.Zero) == Real.Epsilon);
        Assert.True(Real.Epsilon > Real.Zero);
        Assert.True(Real.BitDecrement(Real.Epsilon) < Real.Epsilon);
        Assert.True(Real.BitIncrement(Real.PositiveInfinity).IsPositiveInfinity);
        Assert.True(Real.BitIncrement(Real.NaN).IsNaN);
        Assert.Equal(12.0, Real.ScaleB(R(3), 2).ToDouble());
        foreach (double v in new[] { 1.0, 0.75, 3e100, -5e-300, double.Epsilon, 0.0, double.NaN, double.PositiveInfinity })
        {
            Assert.Equal(double.ILogB(v), Real.ILogB(R(v)));
        }
    }

    [Fact]
    public void Rounding_ModesAndDigits()
    {
        Assert.Equal(2.0, Real.Round(R(2.5)).ToDouble());
        Assert.Equal(3.0, Real.Round(R(2.5), MidpointRounding.AwayFromZero).ToDouble());
        Assert.Equal(-3.0, Real.Round(R(-2.5), MidpointRounding.AwayFromZero).ToDouble());
        Assert.Equal(-2.0, Real.Round(R(-2.5), MidpointRounding.ToZero).ToDouble());
        Assert.Equal(-3.0, Real.Round(R(-2.5), MidpointRounding.ToNegativeInfinity).ToDouble());
        Assert.Equal(-2.0, Real.Round(R(-2.5), MidpointRounding.ToPositiveInfinity).ToDouble());
        Assert.Equal(1.2, Real.Round(R(1.25), 1).ToDouble());
        Assert.Equal(1.3, Real.Round(R(1.25), 1, MidpointRounding.AwayFromZero).ToDouble());
        Assert.Equal(-1.3, Real.Round(R(-1.21), 1, MidpointRounding.ToNegativeInfinity).ToDouble());
        Assert.Equal(2.67, Real.Round(R(2.675), 2).ToDouble()); // 2.675 is really 2.67499999...
        Assert.Equal(0.001, Real.Round(Real.FromBigInt(BigInt.One, -1000, 53), 3, MidpointRounding.ToPositiveInfinity).ToDouble());
        Assert.True(Real.Round(-Real.FromBigInt(BigInt.One, -1000, 10), 3).SignBit);
        Assert.Equal(-2.0, Real.Floor(R(-1.5)).ToDouble());
        Assert.Equal(-1.0, Real.Ceiling(R(-1.5)).ToDouble());
        Assert.Equal(-1.0, Real.Truncate(R(-1.5)).ToDouble());

        foreach (double v in new[] { 0.5, 1.5, 2.5, -0.5, -1.5, 3.7, -3.7, 1e20 + 0.0 })
        {
            foreach (MidpointRounding m in Enum.GetValues<MidpointRounding>())
            {
                Assert.Equal(Math.Round(v, m), Real.Round(R(v), m).ToDouble());
            }
        }
    }

    [Fact]
    public void Predicates()
    {
        Assert.True(Real.IsEvenInteger(R(4)));
        Assert.True(Real.IsEvenInteger(R(0)));
        Assert.False(Real.IsEvenInteger(R(4.5)));
        Assert.True(Real.IsOddInteger(R(-3)));
        Assert.True(Real.IsOddInteger(R(1)));
        Assert.False(Real.IsOddInteger(R(1e20)));
        Assert.True(Real.IsPow2(R(0.25)));
        Assert.False(Real.IsPow2(R(-0.25)));
        Assert.True(Real.IsNegative(Real.NegativeZero));
        Assert.True(Real.IsNormal(R(1e-300)));
        Assert.False(Real.IsNormal(Real.Zero));
        Assert.True(Check<Real>(R(3)));
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
                Same(double.MaxNumber(a, b), Real.MaxNumber(R(a), R(b)));
                Same(double.MinNumber(a, b), Real.MinNumber(R(a), R(b)));
                Same(double.MaxMagnitude(a, b), Real.MaxMagnitude(R(a), R(b)));
                Same(double.MinMagnitude(a, b), Real.MinMagnitude(R(a), R(b)));
                Same(double.MaxMagnitudeNumber(a, b), Real.MaxMagnitudeNumber(R(a), R(b)));
                Same(double.MinMagnitudeNumber(a, b), Real.MinMagnitudeNumber(R(a), R(b)));
            }
        }

        static T Max<T>(T a, T b) where T : INumber<T> => T.Max(a, b);
        static T Min<T>(T a, T b) where T : INumber<T> => T.Min(a, b);
        static void Same(double expected, Real actual) =>
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
        Assert.True(And(Real.NaN, R(1)).IsNaN);
        // ~ complements the 2-bit significand 11 -> 00 (zero) and 10 -> 01.
        Assert.True(Not(R(3, 2)).IsZero);
        Assert.Equal(-1.0, Not(R(2, 2)).ToDouble());
        Assert.True(AllBits<Real>().IsNaN);

        static T And<T>(T a, T b) where T : IBitwiseOperators<T, T, T> => a & b;
        static T Or<T>(T a, T b) where T : IBitwiseOperators<T, T, T> => a | b;
        static T Xor<T>(T a, T b) where T : IBitwiseOperators<T, T, T> => a ^ b;
        static T Not<T>(T a) where T : IBitwiseOperators<T, T, T> => ~a;
        static T AllBits<T>() where T : IBinaryNumber<T> => T.AllBitsSet;
    }

    private static void AssertUlps(double expected, Real actual, int ulps = 1)
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
            AssertUlps(double.AsinPi(x), Real.AsinPi(R(x)), 2);
            AssertUlps(double.AcosPi(x), Real.AcosPi(R(x)), 2);
            double t = (rng.NextDouble() - 0.5) * 100;
            AssertUlps(double.AtanPi(t), Real.AtanPi(R(t)), 2);
            AssertUlps(double.Atan2Pi(x, t), Real.Atan2Pi(R(x), R(t)), 2);
            AssertUlps(double.DegreesToRadians(t), Real.DegreesToRadians(R(t)), 3);
            AssertUlps(double.RadiansToDegrees(x), Real.RadiansToDegrees(R(x)), 3);
        }
        Assert.Equal(0.5, Real.AsinPi(R(1)).ToDouble());
        Assert.Equal(-0.5, Real.AsinPi(R(-1)).ToDouble());
        Assert.Equal(1.0, Real.AcosPi(R(-1)).ToDouble());
        Assert.Equal(0.5, Real.AcosPi(R(0)).ToDouble());
        Assert.True(Real.AcosPi(R(1)).IsZero);
        Assert.Equal(0.25, Real.AtanPi(R(1)).ToDouble());
        Assert.Equal(-0.5, Real.AtanPi(Real.NegativeInfinity).ToDouble());
        Assert.Equal(-0.75, Real.Atan2Pi(R(-2), R(-2)).ToDouble());
        Assert.Equal(1.0, Real.Atan2Pi(R(0), R(-3)).ToDouble());
        Assert.Equal(-1.0, Real.Atan2Pi(R(-0.0), R(-0.0)).ToDouble());
        Assert.Equal(0.75, Real.Atan2Pi(Real.PositiveInfinity, Real.NegativeInfinity).ToDouble());
        Assert.True(Real.AsinPi(R(2)).IsNaN);
        Assert.True(Real.DegreesToRadians(Real.FromBigInt(180, 0, 128)) == Real.Pi(128));
        (Real s, Real c) = Real.SinCosPi(R(0.25));
        Assert.True(s == c);
    }

    [Fact]
    public void LogWithBase()
    {
        Assert.Equal(3.0, Real.Log(R(8), R(2)).ToDouble());
        Assert.Equal(-3.0, Real.Log(R(0.125), R(2)).ToDouble());
        Assert.Equal(1.5, Real.Log(R(8), R(4)).ToDouble());
        Assert.Equal(0.5, Real.Log(R(3), R(9)).ToDouble());
        Assert.Equal(2.0, Real.Log(R(100), R(10)).ToDouble());
        Assert.Equal(5.0, Real.Log(R(243), R(3)).ToDouble());
        Assert.True(Real.Log(R(1), R(7)).IsZero);
        Assert.True(Real.Log(R(5), R(1)).IsNaN);
        Assert.True(Real.Log(R(-5), R(2)).IsNaN);
        Assert.True(Real.Log(R(0), R(2)).IsNegativeInfinity);
        var rng = new Random(32);
        for (int i = 0; i < 100; i++)
        {
            double x = rng.NextDouble() * 1000, b = rng.NextDouble() * 20 + 0.01;
            AssertUlps(Math.Log(x) / Math.Log(b), Real.Log(R(x), R(b)), 4);
        }
    }

    [Fact]
    public void FunctionOverloads_MatchPrecisionVersions()
    {
        Real x = Real.FromDouble(0.7, 100);
        Assert.Equal(Real.Exp(x, 0), Real.Exp(x));
        Assert.Equal(Real.ExpM1(x, 0), Exp<Real>(x));
        Assert.Equal(Real.Exp2M1(x, 0), Real.Exp2M1(x));
        Assert.Equal(Real.Exp10M1(x, 0), Real.Exp10M1(x));
        Assert.Equal(Real.LogP1(x, 0), Real.LogP1(x));
        Assert.Equal(Real.Log2P1(x, 0), Real.Log2P1(x));
        Assert.Equal(Real.Log10P1(x, 0), Real.Log10P1(x));
        Assert.Equal(Real.Sin(x, 0), Real.Sin(x));
        Assert.Equal(Real.Atanh(x, 0), Real.Atanh(x));
        Assert.Equal(Real.RootN(x, 5L, 0), Real.RootN(x, 5));
        Assert.Equal(Real.Pow(x, x, 0), Real.Pow(x, x));
        Assert.Equal(Real.Hypot(x, x, 0), Real.Hypot(x, x));
        Assert.Equal(100, Real.Exp(x).Precision);
        Assert.Equal(Real.Pi(0).ShiftLeft(1), Tau<Real>());
        Assert.Equal(Real.E(0), E<Real>());
        Assert.True(Real.FusedMultiplyAdd(x, x, x) == Real.FusedMultiplyAdd(x, x, x));
        Assert.True(Real.Ieee754Remainder(R(7), R(2)) == R(-1));
        Assert.True(R(-7) % R(2) == R(-1));
        Real y = R(1.5);
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
        IFloatingPoint<Real> one = R(1);
        Assert.Equal(0, one.GetExponentShortestBitLength());
        Assert.Equal(53, one.GetSignificandBitLength());
        Assert.Equal(7, one.GetSignificandByteCount());
        byte[] buf = new byte[8];
        Assert.True(((IFloatingPoint<Real>)R(0.75)).TryWriteExponentLittleEndian(buf, out int n));
        Assert.Equal(8, n);
        Assert.Equal(-1L, BitConverter.ToInt64(buf));
        Assert.True(((IFloatingPoint<Real>)R(3, 8)).TryWriteSignificandLittleEndian(buf, out n));
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
        Real r = Real.Parse(s, style, CultureInfo.InvariantCulture);
        Assert.Equal(expected, r.ToDouble());
        Assert.Equal(expected, Real.Parse(Encoding.UTF8.GetBytes(s), style, CultureInfo.InvariantCulture).ToDouble());
        Assert.Equal(expected, double.Parse(s.Replace("inf", "Infinity"), style, CultureInfo.InvariantCulture));
    }

    [Fact]
    public void ParseStyles_Rejects()
    {
        Assert.False(Real.TryParse("1,5", NumberStyles.Float, CultureInfo.InvariantCulture, out _));
        Assert.False(Real.TryParse("1.5", NumberStyles.Integer, CultureInfo.InvariantCulture, out _));
        Assert.Throws<ArgumentException>(() => Real.TryParse("ff", NumberStyles.HexNumber, null, out _));
        Assert.Equal(1.5, Real.Parse("1,5", NumberStyles.Float, CultureInfo.GetCultureInfo("de-DE")).ToDouble());
    }

    [Fact]
    public void Utf8()
    {
        Real x = R(0.375);
        byte[] buf = new byte[64];
        Assert.True(x.TryFormat(buf, out int n, default, null));
        Assert.Equal(x.ToString(), Encoding.UTF8.GetString(buf, 0, n));
        Assert.True(Real.Parse("0.375"u8, null) == x);
    }
}
