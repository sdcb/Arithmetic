using System.Globalization;
using System.Numerics;
using System.Text;
using Xunit;
using Xunit.Abstractions;

namespace Sdcb.Arithmetic.Tests.Floats;

public class QuadDoubleTests(ITestOutputHelper output)
{
    private const int P = 500;

    // ---------------------------------------------------------------
    // Helpers
    // ---------------------------------------------------------------

    /// <summary>log2 of the relative error |actual - expected| / |expected| (-1000 when exact).</summary>
    private static double ErrBits(QuadDouble actual, BigFloat expected)
    {
        BigFloat a = actual;
        BigFloat diff = BigFloat.Subtract(a, expected, P);
        if (diff.IsZero)
        {
            return -1000;
        }
        if (expected.IsZero)
        {
            return 1000;
        }
        BigFloat rel = BigFloat.Divide(diff, expected, 64).Abs();
        return Math.Log2(rel.ToDouble());
    }

    private static double ErrBits(QuadDouble actual, BigFloat expected, BigFloat scale)
    {
        BigFloat diff = BigFloat.Subtract(actual, expected, P);
        if (diff.IsZero)
        {
            return -1000;
        }
        BigFloat rel = BigFloat.Divide(diff, scale, 64).Abs();
        return Math.Log2(rel.ToDouble());
    }

    /// <summary>A random normalized QD with full 212-bit significand and exponent in [emin, emax].</summary>
    private static QuadDouble RandomQd(Random rng, int emin, int emax, bool allowNegative = true)
    {
        int e = rng.Next(emin, emax + 1);
        double m0 = 1 + rng.NextDouble();
        double m1 = (rng.NextDouble() - 0.5) * Math.ScaleB(1, -52);
        double m2 = (rng.NextDouble() - 0.5) * Math.ScaleB(1, -105);
        double m3 = (rng.NextDouble() - 0.5) * Math.ScaleB(1, -158);
        QuadDouble q = new(m0, m1, m2, m3);
        q = QuadDouble.ScaleB(q, e);
        return allowNegative && rng.Next(2) == 0 ? -q : q;
    }

    private void Check(string name, double maxBits, double bound)
    {
        output.WriteLine($"{name}: max error 2^{maxBits:F2}");
        Assert.True(maxBits <= bound, $"{name}: max relative error 2^{maxBits:F2} exceeds 2^{bound}");
    }

    private static BigFloat R(QuadDouble q) => q;

    // ---------------------------------------------------------------
    // Constants
    // ---------------------------------------------------------------

    [Fact]
    public void ConstantsMatchCorrectlyRoundedBigFloat()
    {
        AssertSplit(BigFloat.Pi(800), QuadDouble.Pi);
        AssertSplit(BigFloat.Pi(800).ShiftLeft(1), QuadDouble.Tau);
        AssertSplit(BigFloat.E(800), QuadDouble.E);
        AssertSplit(BigFloat.Ln2(800), QuadDouble.Ln2);
        AssertSplit(BigFloat.Log(BigFloat.FromDouble(10, 800), 800), QuadDouble.Ln10);
    }

    private static void AssertSplit(BigFloat exact, QuadDouble q)
    {
        QuadDouble expected = (QuadDouble)exact;
        Assert.Equal(expected.X0, q.X0);
        Assert.Equal(expected.X1, q.X1);
        Assert.Equal(expected.X2, q.X2);
        Assert.Equal(expected.X3, q.X3);
    }

    [Fact]
    public void SpecialConstants()
    {
        Assert.True(QuadDouble.IsNaN(QuadDouble.NaN));
        Assert.True(QuadDouble.IsPositiveInfinity(QuadDouble.PositiveInfinity));
        Assert.True(QuadDouble.IsNegativeInfinity(QuadDouble.NegativeInfinity));
        Assert.Equal(double.Epsilon, QuadDouble.Epsilon.X0);
        Assert.Equal(double.MaxValue, (double)QuadDouble.MaxValue);
        Assert.True(QuadDouble.MaxValue > double.MaxValue);
        Assert.Equal(-QuadDouble.MaxValue, QuadDouble.MinValue);
        BigFloat max = QuadDouble.MaxValue;
        BigFloat expected = BigFloat.Multiply(BigFloat.FromDouble(9007199254740991, 64), BigFloat.Add(BigFloat.Add(BigFloat.One.ShiftLeft(971), BigFloat.One.ShiftLeft(917), 400), BigFloat.Add(BigFloat.One.ShiftLeft(863), BigFloat.One.ShiftLeft(809), 400), 400), 400);
        Assert.True(max == expected);
        Assert.Equal(QuadDouble.One, -QuadDouble.NegativeOne);
    }

    // ---------------------------------------------------------------
    // Basic arithmetic vs BigFloat
    // ---------------------------------------------------------------

    [Fact]
    public void AddSubtractRandom()
    {
        var rng = new Random(1234);
        double worst = -1000;
        for (int i = 0; i < 4000; i++)
        {
            QuadDouble a = RandomQd(rng, -60, 60);
            QuadDouble b;
            switch (i % 4)
            {
                case 0:
                    b = RandomQd(rng, -60, 60);
                    break;
                case 1:
                    // Heavy cancellation: b ≈ -a with a perturbation at a random depth.
                    b = -a + RandomQd(rng, -60, 60) * Math.ScaleB(1, -rng.Next(1, 200));
                    break;
                case 2:
                    b = -QuadDouble.ScaleB(a, rng.Next(-3, 4)) + QuadDouble.ScaleB(RandomQd(rng, 0, 0), QuadDouble.ILogB(a) - rng.Next(40, 230));
                    break;
                default:
                    b = RandomQd(rng, -300, 300);
                    break;
            }
            BigFloat ea = BigFloat.Add(a, b, P);
            if (!ea.IsZero)
            {
                worst = Math.Max(worst, ErrBits(a + b, ea));
            }
            BigFloat es = BigFloat.Subtract(a, b, P);
            if (!es.IsZero)
            {
                worst = Math.Max(worst, ErrBits(a - b, es));
            }
            double d = RandomQd(rng, -60, 60).X0;
            worst = Math.Max(worst, ErrBits(a + d, BigFloat.Add(a, BigFloat.FromDouble(d, 64), P)));
        }
        Check("add/sub", worst, -208);
    }

    [Fact]
    public void AddExactCancellation()
    {
        QuadDouble a = new QuadDouble(1) / 3;
        Assert.Equal(QuadDouble.Zero, a - a);
        Assert.False(QuadDouble.IsNegative(a - a));
        QuadDouble x = new(1, Math.ScaleB(1, -80));
        Assert.Equal(Math.ScaleB(1, -80), (double)(x - 1));
    }

    [Fact]
    public void MultiplyRandom()
    {
        var rng = new Random(42);
        double worst = -1000, worstD = -1000, worstSq = -1000;
        for (int i = 0; i < 4000; i++)
        {
            QuadDouble a = RandomQd(rng, -200, 200), b = RandomQd(rng, -200, 200);
            worst = Math.Max(worst, ErrBits(a * b, BigFloat.Multiply(a, b, P)));
            double d = RandomQd(rng, -50, 50).X0;
            worstD = Math.Max(worstD, ErrBits(a * d, BigFloat.Multiply(a, BigFloat.FromDouble(d, 64), P)));
            worstSq = Math.Max(worstSq, ErrBits(QuadDouble.Square(a), BigFloat.Multiply(a, a, P)));
        }
        Check("mul", worst, -208);
        Check("mul double", worstD, -208);
        Check("square", worstSq, -208);
    }

    [Fact]
    public void DivideRandom()
    {
        var rng = new Random(7);
        double worst = -1000, worstD = -1000;
        for (int i = 0; i < 4000; i++)
        {
            QuadDouble a = RandomQd(rng, -200, 200), b = RandomQd(rng, -200, 200);
            worst = Math.Max(worst, ErrBits(a / b, BigFloat.Divide(a, b, P)));
            double d = RandomQd(rng, -50, 50).X0;
            worstD = Math.Max(worstD, ErrBits(a / d, BigFloat.Divide(a, BigFloat.FromDouble(d, 64), P)));
        }
        Check("div", worst, -208);
        Check("div double", worstD, -208);
    }

    [Fact]
    public void SqrtCbrtRandom()
    {
        var rng = new Random(99);
        double worst = -1000, worstC = -1000;
        for (int i = 0; i < 3000; i++)
        {
            QuadDouble a = RandomQd(rng, -1000, 1000, allowNegative: false);
            worst = Math.Max(worst, ErrBits(QuadDouble.Sqrt(a), BigFloat.Sqrt(a, P)));
            QuadDouble c = i % 2 == 0 ? a : -a;
            worstC = Math.Max(worstC, ErrBits(QuadDouble.Cbrt(c), BigFloat.Cbrt(c, P)));
        }
        Check("sqrt", worst, -208);
        Check("cbrt", worstC, -206);
        Assert.Equal(QuadDouble.Zero, QuadDouble.Sqrt(0));
        Assert.True(QuadDouble.IsNegative(QuadDouble.Sqrt(QuadDouble.NegativeZero)));
        Assert.True(QuadDouble.IsNaN(QuadDouble.Sqrt(-1)));
        Assert.True(QuadDouble.IsPositiveInfinity(QuadDouble.Sqrt(double.PositiveInfinity)));
        Assert.Equal(new QuadDouble(3), QuadDouble.Sqrt(9));
    }

    /// <summary>A QD whose significand has only the leading <paramref name="bits"/> bits (an exact multiple of 2^(e - bits + 1)).</summary>
    private static QuadDouble ShortQd(Random rng, int bits, int e)
    {
        BigFloat m = BigFloat.Zero;
        for (int i = 0; i < bits; i += 50)
        {
            int n = Math.Min(50, bits - i);
            BigFloat chunk = BigFloat.FromDouble(rng.NextInt64(1L << (n - 1), 1L << n), 64);
            m = BigFloat.Add(m.ShiftLeft(n), chunk, P);
        }
        return QuadDouble.ScaleB((QuadDouble)m, e - (bits - 1));
    }

    [Fact]
    public void DivideStress()
    {
        var rng = new Random(2024);
        double wRandom = -1000, wPow2 = -1000, wNearOne = -1000, wShort = -1000, wExtreme = -1000, wExact = -1000;
        for (int i = 0; i < 3000; i++)
        {
            // Random operands across a wide exponent range.
            // (Quotients stay above 2^-800, inside the documented full-precision range.)
            QuadDouble a = RandomQd(rng, -400, 400), b = RandomQd(rng, -400, 400);
            wRandom = Math.Max(wRandom, ErrBits(a / b, BigFloat.Divide(a, b, P)));

            // Divisors at or next to a power of two (the leading component may be 2^k exactly).
            int k = rng.Next(-300, 300);
            QuadDouble p2 = QuadDouble.ScaleB(QuadDouble.One, k);
            QuadDouble tiny = QuadDouble.ScaleB(RandomQd(rng, 0, 0), k - rng.Next(1, 215));
            QuadDouble bp = (i % 3) switch { 0 => p2, 1 => p2 + tiny, _ => p2 - tiny };
            wPow2 = Math.Max(wPow2, ErrBits(a / bp, BigFloat.Divide(a, bp, P)));
            wPow2 = Math.Max(wPow2, ErrBits(bp / a, BigFloat.Divide(bp, a, P)));

            // Quotients close to one: b = a (1 + δ) with δ at a random depth (also δ = 0).
            QuadDouble near = i % 50 == 0 ? a : a + a * QuadDouble.ScaleB(RandomQd(rng, 0, 0), -rng.Next(1, 230));
            wNearOne = Math.Max(wNearOne, ErrBits(a / near, BigFloat.Divide(a, near, P)));
            wNearOne = Math.Max(wNearOne, ErrBits(near / a, BigFloat.Divide(near, a, P)));

            // Operands with short significands (zero lower components).
            QuadDouble sa = ShortQd(rng, rng.Next(1, 160), rng.Next(-50, 50));
            QuadDouble sb = ShortQd(rng, rng.Next(1, 160), rng.Next(-50, 50));
            wShort = Math.Max(wShort, ErrBits(sa / sb, BigFloat.Divide(sa, sb, P)));
            wShort = Math.Max(wShort, ErrBits(a / sb, BigFloat.Divide(a, sb, P)));
            wShort = Math.Max(wShort, ErrBits(sa / b, BigFloat.Divide(sa, b, P)));

            // Very large and very small magnitudes with an in-range quotient.
            QuadDouble big = RandomQd(rng, 700, 1000), small = RandomQd(rng, -780, -700);
            QuadDouble mid = RandomQd(rng, -20, 20), hb = RandomQd(rng, 600, 700);
            wExtreme = Math.Max(wExtreme, ErrBits(big / mid, BigFloat.Divide(big, mid, P)));
            wExtreme = Math.Max(wExtreme, ErrBits(small / mid, BigFloat.Divide(small, mid, P)));
            wExtreme = Math.Max(wExtreme, ErrBits(mid / small, BigFloat.Divide(mid, small, P)));
            wExtreme = Math.Max(wExtreme, ErrBits(big / hb, BigFloat.Divide(big, hb, P)));
            wExtreme = Math.Max(wExtreme, ErrBits(hb / big, BigFloat.Divide(hb, big, P)));

            // Exact quotients: (c·m) / m with a short c must give c exactly.
            QuadDouble c = ShortQd(rng, rng.Next(1, 50), rng.Next(-30, 30));
            QuadDouble m = ShortQd(rng, rng.Next(1, 100), rng.Next(-30, 30));
            QuadDouble cm = c * m; // exact: at most 150 significant bits
            Assert.True(R(cm) == BigFloat.Multiply(c, m, P));
            QuadDouble qq = cm / m;
            wExact = Math.Max(wExact, ErrBits(qq, c));
        }
        Check("div random", wRandom, -208);
        Check("div near 2^k", wPow2, -208);
        Check("div quotient near 1", wNearOne, -208);
        Check("div short operands", wShort, -208);
        Check("div extreme magnitudes", wExtreme, -208);
        Check("div exact quotients", wExact, -208);

        // Special values.
        QuadDouble one = 1, two = 2, inf = QuadDouble.PositiveInfinity;
        Assert.Equal(new QuadDouble(0.5), one / two);
        Assert.Equal(new QuadDouble(3), new QuadDouble(6) / two);
        Assert.True(QuadDouble.IsNaN(inf / inf));
        Assert.True(QuadDouble.IsNaN(QuadDouble.NaN / one));
        Assert.True(QuadDouble.IsNaN(one / QuadDouble.NaN));
        Assert.True(QuadDouble.IsPositiveInfinity(inf / two));
        Assert.True(QuadDouble.IsNegativeInfinity(inf / -two));
        Assert.True(QuadDouble.IsNegative(-one / inf));
        Assert.Equal(QuadDouble.Zero, QuadDouble.Zero / two);
        Assert.True(QuadDouble.IsNegative(QuadDouble.NegativeZero / two));
        Assert.True(QuadDouble.IsNegative(QuadDouble.Zero / -two));
        Assert.True(QuadDouble.IsPositiveInfinity(QuadDouble.MaxValue / 0.5));
        Assert.Equal(0.0, (QuadDouble.MaxValue / 0.5).X1);
        Assert.True(QuadDouble.IsFinite(QuadDouble.MaxValue / 2));
        Assert.True(R(QuadDouble.MaxValue / 1) == R(QuadDouble.MaxValue));
    }

    [Fact]
    public void SqrtStress()
    {
        var rng = new Random(4096);
        double wRandom = -1000, wSquares = -1000, wNear = -1000, wShort = -1000;
        for (int i = 0; i < 3000; i++)
        {
            // Random operands over the whole exponent range, odd and even exponents.
            QuadDouble a = RandomQd(rng, -1020, 1020, allowNegative: false);
            wRandom = Math.Max(wRandom, ErrBits(QuadDouble.Sqrt(a), BigFloat.Sqrt(a, P)));
            QuadDouble a1 = RandomQd(rng, -3, 3, allowNegative: false);
            wRandom = Math.Max(wRandom, ErrBits(QuadDouble.Sqrt(a1), BigFloat.Sqrt(a1, P)));

            // Perfect squares n² (n with up to 106 bits, so n² is exact in QD): the root is n.
            QuadDouble n = ShortQd(rng, rng.Next(1, 107), rng.Next(-400, 400));
            BigFloat n2 = BigFloat.Multiply(n, n, P);
            QuadDouble sq = (QuadDouble)n2;
            Assert.True(R(sq) == n2);
            wSquares = Math.Max(wSquares, ErrBits(QuadDouble.Sqrt(sq), n));

            // Arguments next to powers of four and of two (the root's leading component is near 2^k).
            int k = rng.Next(-500, 500);
            QuadDouble p = QuadDouble.ScaleB(QuadDouble.One, k);
            QuadDouble d = QuadDouble.ScaleB(RandomQd(rng, 0, 0, allowNegative: false), k - rng.Next(1, 215));
            QuadDouble x = (i % 3) switch { 0 => p + d, 1 => p - d, _ => p };
            wNear = Math.Max(wNear, ErrBits(QuadDouble.Sqrt(x), BigFloat.Sqrt(x, P)));

            // Short significands.
            QuadDouble s = ShortQd(rng, rng.Next(1, 212), rng.Next(-900, 900));
            wShort = Math.Max(wShort, ErrBits(QuadDouble.Sqrt(s), BigFloat.Sqrt(s, P)));
        }
        Check("sqrt random", wRandom, -208);
        Check("sqrt perfect squares", wSquares, -208);
        Check("sqrt near 2^k", wNear, -208);
        Check("sqrt short operands", wShort, -208);

        for (int v = 1; v < 3000; v++)
        {
            Assert.Equal(new QuadDouble(v), QuadDouble.Sqrt(new QuadDouble((double)v * v)));
        }
        Assert.Equal(new QuadDouble(Math.ScaleB(1, -537)), QuadDouble.Sqrt(QuadDouble.Epsilon));
        Assert.Equal(new QuadDouble(Math.ScaleB(1, 511)), QuadDouble.Sqrt(new QuadDouble(Math.ScaleB(1, 1022))));
        Assert.True(ErrBits(QuadDouble.Sqrt(QuadDouble.MaxValue), BigFloat.Sqrt(QuadDouble.MaxValue, P)) <= -208);
        Assert.True(ErrBits(QuadDouble.Sqrt(double.MaxValue), BigFloat.Sqrt(BigFloat.FromDouble(double.MaxValue, 64), P)) <= -208);
        Assert.True(QuadDouble.IsNaN(QuadDouble.Sqrt(QuadDouble.NaN)));
        Assert.True(QuadDouble.IsNaN(QuadDouble.Sqrt(QuadDouble.NegativeInfinity)));
        Assert.True(QuadDouble.IsNaN(QuadDouble.Sqrt(-QuadDouble.Epsilon)));
        Assert.True(QuadDouble.IsPositiveInfinity(QuadDouble.Sqrt(QuadDouble.PositiveInfinity)));
        Assert.Equal(0.0, QuadDouble.Sqrt(QuadDouble.PositiveInfinity).X1);
    }

    [Fact]
    public void AddStress()
    {
        var rng = new Random(77);
        double wSame = -1000, wMixed = -1000, wShort = -1000;
        for (int i = 0; i < 4000; i++)
        {
            // Same sign, any exponent gap (including gaps beyond 212 bits).
            QuadDouble a = RandomQd(rng, -100, 100, allowNegative: false);
            QuadDouble b = QuadDouble.ScaleB(RandomQd(rng, 0, 0, allowNegative: false), QuadDouble.ILogB(a) - rng.Next(-260, 260));
            if (i % 2 == 0)
            {
                (a, b) = (-a, -b);
            }
            wSame = Math.Max(wSame, ErrBits(a + b, BigFloat.Add(a, b, P)));

            // Opposite signs with partial cancellation at the leading components.
            QuadDouble c = -QuadDouble.ScaleB(a, rng.Next(-2, 3)) + QuadDouble.ScaleB(RandomQd(rng, 0, 0), QuadDouble.ILogB(a) - rng.Next(0, 230));
            BigFloat ec = BigFloat.Add(a, c, P);
            if (!ec.IsZero)
            {
                wMixed = Math.Max(wMixed, ErrBits(a + c, ec));
            }

            // Short significands: zero lower components on one or both sides.
            QuadDouble s = ShortQd(rng, rng.Next(1, 212), rng.Next(-120, 120));
            QuadDouble t = rng.Next(2) == 0 ? -s : s;
            QuadDouble u = ShortQd(rng, rng.Next(1, 212), rng.Next(-120, 120));
            BigFloat es = BigFloat.Add(t, u, P);
            if (!es.IsZero)
            {
                wShort = Math.Max(wShort, ErrBits(t + u, es));
            }
            BigFloat ea = BigFloat.Add(a, u, P);
            if (!ea.IsZero)
            {
                wShort = Math.Max(wShort, ErrBits(a + u, ea));
            }
        }
        Check("add same sign", wSame, -208);
        Check("add mixed signs", wMixed, -208);
        Check("add short operands", wShort, -208);
        Assert.Equal(QuadDouble.Zero, QuadDouble.Pi + -QuadDouble.Pi);
        Assert.True(QuadDouble.IsNegative(QuadDouble.NegativeZero + QuadDouble.NegativeZero));
        Assert.True(QuadDouble.IsPositiveInfinity(QuadDouble.MaxValue + QuadDouble.MaxValue));
        Assert.Equal(0.0, (QuadDouble.MaxValue + QuadDouble.MaxValue).X1);
        Assert.True(QuadDouble.IsNegativeInfinity(QuadDouble.MinValue + QuadDouble.MinValue));
        Assert.True(QuadDouble.IsNaN(QuadDouble.PositiveInfinity + QuadDouble.NegativeInfinity));
        Assert.Equal(QuadDouble.MaxValue, QuadDouble.MaxValue + QuadDouble.Zero);
    }

    [Fact]
    public void NearOverflowAndUnderflow()
    {
        var rng = new Random(5);
        double worst = -1000;
        for (int i = 0; i < 500; i++)
        {
            QuadDouble a = RandomQd(rng, 1000, 1018), b = RandomQd(rng, -5, 2);
            worst = Math.Max(worst, ErrBits(a * b, BigFloat.Multiply(a, b, P)));
            worst = Math.Max(worst, ErrBits(a / b, BigFloat.Divide(a, b, P)));
            QuadDouble c = RandomQd(rng, 1000, 1022);
            BigFloat s = BigFloat.Add(a, c, P);
            if (s.BinaryExponent < 1024)
            {
                worst = Math.Max(worst, ErrBits(a + c, s));
            }
            QuadDouble t = RandomQd(rng, -840, -800);
            worst = Math.Max(worst, ErrBits(t * b, BigFloat.Multiply(t, b, P)));
        }
        Check("near overflow/underflow", worst, -206);
        Assert.True(QuadDouble.IsPositiveInfinity(QuadDouble.MaxValue * 2));
        Assert.True(QuadDouble.IsPositiveInfinity(QuadDouble.MaxValue + QuadDouble.MaxValue));
        Assert.Equal(QuadDouble.Zero, QuadDouble.Epsilon / 4);
    }

    [Fact]
    public void SpecialValueArithmetic()
    {
        QuadDouble inf = QuadDouble.PositiveInfinity, nan = QuadDouble.NaN, one = 1;
        Assert.True(QuadDouble.IsNaN(inf - inf));
        Assert.True(QuadDouble.IsNaN(inf * 0));
        Assert.True(QuadDouble.IsNaN(nan + one));
        Assert.True(QuadDouble.IsNaN(QuadDouble.Zero / QuadDouble.Zero));
        Assert.True(QuadDouble.IsPositiveInfinity(one / QuadDouble.Zero));
        Assert.True(QuadDouble.IsNegativeInfinity(one / QuadDouble.NegativeZero));
        Assert.Equal(QuadDouble.Zero, one / inf);
        Assert.True(QuadDouble.IsPositiveInfinity(inf + one));
        Assert.True(QuadDouble.IsNegative(QuadDouble.NegativeZero * one));
        Assert.True(QuadDouble.IsNegative(QuadDouble.NegativeZero + QuadDouble.NegativeZero));
        Assert.False(QuadDouble.IsNegative(QuadDouble.NegativeZero + QuadDouble.Zero));
        Assert.Equal(0.0, (inf * 2).X1);
    }

    [Fact]
    public void RemainderAndFma()
    {
        var rng = new Random(11);
        for (int i = 0; i < 300; i++)
        {
            QuadDouble a = RandomQd(rng, -10, 60), b = RandomQd(rng, -10, 10);
            QuadDouble r = a % b;
            BigFloat exact = BigFloat.TruncatedRemainder(a, b, 2000);
            Assert.True(ErrBits(r, exact) < -210, $"{a} % {b}");
        }
        Assert.Equal(new QuadDouble(1), new QuadDouble(7) % 3);
        Assert.Equal(new QuadDouble(-1), new QuadDouble(-7) % 3);
        Assert.True(QuadDouble.IsNaN(new QuadDouble(1) % 0));

        QuadDouble x = QuadDouble.One / 3, y = 3;
        QuadDouble f = QuadDouble.FusedMultiplyAdd(x, y, -1);
        BigFloat fe = BigFloat.FusedMultiplyAdd(x, y, BigFloat.MinusOne, 600);
        Assert.True(ErrBits(f, fe) < -200);
    }

    // ---------------------------------------------------------------
    // Rounding, comparison, classification
    // ---------------------------------------------------------------

    [Fact]
    public void RoundingFunctions()
    {
        QuadDouble big = new QuadDouble(Math.ScaleB(1, 80)) + 0.5;
        Assert.Equal(new QuadDouble(Math.ScaleB(1, 80)), QuadDouble.Floor(big));
        Assert.Equal(new QuadDouble(Math.ScaleB(1, 80)) + 1, QuadDouble.Ceiling(big));
        Assert.Equal(new QuadDouble(Math.ScaleB(1, 80)), QuadDouble.Round(big)); // tie to even
        Assert.Equal(new QuadDouble(Math.ScaleB(1, 80)) + 1, QuadDouble.Round(big, MidpointRounding.AwayFromZero));
        Assert.Equal(new QuadDouble(2), QuadDouble.Round(new QuadDouble(2.5)));
        Assert.Equal(new QuadDouble(4), QuadDouble.Round(new QuadDouble(3.5)));
        Assert.Equal(new QuadDouble(-4), QuadDouble.Round(new QuadDouble(-3.5)));
        Assert.Equal(new QuadDouble(-3), QuadDouble.Truncate(new QuadDouble(-3.7)));
        Assert.True(QuadDouble.IsNegative(QuadDouble.Round(new QuadDouble(-0.3))));
        Assert.True(QuadDouble.IsNegative(QuadDouble.Floor(QuadDouble.NegativeZero)));
        Assert.True(QuadDouble.IsNegative(QuadDouble.Ceiling(new QuadDouble(-0.5))));
        // 2.5 + tiny rounds up even though it looks like a tie in the leading component.
        Assert.Equal(new QuadDouble(3), QuadDouble.Round(new QuadDouble(2.5, Math.ScaleB(1, -100))));
        Assert.Equal(new QuadDouble(2), QuadDouble.Round(new QuadDouble(2.5, -Math.ScaleB(1, -100))));
        QuadDouble third = QuadDouble.One / 3;
        Assert.True(ErrBits(QuadDouble.Round(third, 3), BigFloat.Divide(BigFloat.FromDouble(333, 64), BigFloat.FromDouble(1000, 64), P)) < -212);
        Assert.True(ErrBits(QuadDouble.Round(third, 2, MidpointRounding.ToPositiveInfinity), BigFloat.Divide(BigFloat.FromDouble(34, 64), BigFloat.FromDouble(100, 64), P)) < -212);
        Assert.True(QuadDouble.IsInteger(big - 0.5));
        Assert.False(QuadDouble.IsInteger(big));
        Assert.True(QuadDouble.IsOddInteger(new QuadDouble(Math.ScaleB(1, 80)) + 1));
        Assert.True(QuadDouble.IsEvenInteger(new QuadDouble(Math.ScaleB(1, 80))));
    }

    [Fact]
    public void ComparisonAndEquality()
    {
        QuadDouble a = new(1, Math.ScaleB(1, -60)), b = new(1, Math.ScaleB(1, -61));
        Assert.True(a > b);
        Assert.True(b < a);
        Assert.True(a != b);
        QuadDouble a2 = a;
        Assert.True(a >= a2 && a <= a2);
        Assert.True(QuadDouble.Zero == QuadDouble.NegativeZero);
        Assert.Equal(QuadDouble.Zero.GetHashCode(), QuadDouble.NegativeZero.GetHashCode());
        Assert.False(QuadDouble.NaN == QuadDouble.NaN);
        Assert.True(QuadDouble.NaN.Equals(QuadDouble.NaN));
        Assert.True(QuadDouble.NaN.CompareTo(QuadDouble.Zero) < 0);
        Assert.False(QuadDouble.NaN < 1 || QuadDouble.NaN > 1 || QuadDouble.NaN <= 1);
        // Non-canonical representations of one value compare and hash equal.
        double u = Math.ScaleB(1, -53);
        QuadDouble c1 = new(1 + 2 * u, -u / 2, 0, 0);
        QuadDouble c2 = new(1 + u * 2, -u / 2);
        Assert.Equal(c1, c2);
        Assert.Equal(c1.GetHashCode(), c2.GetHashCode());
        Assert.Equal(QuadDouble.One, QuadDouble.Max(QuadDouble.One, QuadDouble.Zero));
        Assert.True(QuadDouble.IsNaN(QuadDouble.Max(QuadDouble.NaN, 1)));
        Assert.Equal(QuadDouble.One, QuadDouble.MaxNumber(QuadDouble.NaN, 1));
        Assert.False(QuadDouble.IsNegative(QuadDouble.Max(QuadDouble.NegativeZero, QuadDouble.Zero)));
        Assert.True(QuadDouble.IsNegative(QuadDouble.Min(QuadDouble.NegativeZero, QuadDouble.Zero)));
        Assert.Equal(new QuadDouble(-3), QuadDouble.MaxMagnitude(2, -3));
        Assert.Equal(new QuadDouble(2), QuadDouble.MinMagnitude(2, -3));
        Assert.Equal(-1, QuadDouble.Sign(-a));
        Assert.Throws<ArithmeticException>(() => QuadDouble.Sign(QuadDouble.NaN));
        Assert.Equal(-a, QuadDouble.CopySign(a, -1));
        Assert.Equal(0, QuadDouble.ILogB(a));
        Assert.Equal(-1, QuadDouble.ILogB(new QuadDouble(1, -Math.ScaleB(1, -70))));
        Assert.Equal(int.MinValue, QuadDouble.ILogB(0));
    }

    // ---------------------------------------------------------------
    // Conversions
    // ---------------------------------------------------------------

    [Fact]
    public void IntegerConversionsAreExact()
    {
        Int128[] values = [Int128.MaxValue, Int128.MinValue, Int128.MaxValue - 12345, (Int128)long.MaxValue, (Int128)long.MinValue, 0, -1, (Int128)1 << 100];
        foreach (Int128 v in values)
        {
            QuadDouble q = v;
            Assert.True(R(q) == (BigFloat)v, v.ToString());
            Assert.Equal(v, (Int128)q);
        }
        UInt128[] uvalues = [UInt128.MaxValue, UInt128.MaxValue - 1, (UInt128)Int128.MaxValue + 1, ulong.MaxValue];
        foreach (UInt128 v in uvalues)
        {
            QuadDouble q = v;
            Assert.True(R(q) == (BigFloat)v, v.ToString());
            Assert.Equal(v, (UInt128)q);
        }
        Assert.Equal(long.MaxValue, (long)(QuadDouble)long.MaxValue);
        Assert.Equal(long.MinValue, (long)(QuadDouble)long.MinValue);
        Assert.Equal(ulong.MaxValue, (ulong)(QuadDouble)ulong.MaxValue);
        Assert.Equal(long.MaxValue - 1, (long)((QuadDouble)long.MaxValue - 1));
        // Saturation like double.
        Assert.Equal(int.MaxValue, (int)new QuadDouble(1e20));
        Assert.Equal(0, (int)QuadDouble.NaN);
        Assert.Equal(long.MinValue, (long)QuadDouble.NegativeInfinity);
        Assert.Equal(Int128.MaxValue, (Int128)new QuadDouble(1e300));
        Assert.Equal(UInt128.MinValue, (UInt128)new QuadDouble(-5));
        Assert.Equal(-7, (int)new QuadDouble(-7.9));
        Assert.Throws<OverflowException>(() => int.CreateChecked(new QuadDouble(1e20)));
        Assert.Equal(int.MaxValue, int.CreateSaturating(new QuadDouble(1e20)));
        Assert.Equal(12345678901234567L, long.CreateChecked((QuadDouble)12345678901234567L));
        Assert.Equal((QuadDouble)12345678901234567L, QuadDouble.CreateChecked(12345678901234567L));
        Assert.Equal((QuadDouble)Int128.MaxValue, QuadDouble.CreateTruncating(Int128.MaxValue));
    }

    [Fact]
    public void FloatingConversions()
    {
        var rng = new Random(3);
        for (int i = 0; i < 500; i++)
        {
            QuadDouble q = RandomQd(rng, -1000, 1000);
            BigFloat exact = q;
            Assert.Equal(exact.ToDouble(), (double)q);
            Assert.Equal(exact.ToSingle(), (float)q);
            Assert.Equal(exact.ToHalf(), (Half)q);
        }
        // Tie handling: 1 + 2^-53 + tiny rounds up, exact tie rounds to even.
        double u = Math.ScaleB(1, -53);
        Assert.Equal(1 + 2 * u, (double)new QuadDouble(1, u, Math.ScaleB(1, -120)));
        Assert.Equal(1.0, (double)new QuadDouble(1, u));
        Assert.Equal(1.0, (double)new QuadDouble(1, u, -Math.ScaleB(1, -120)));
        // float double-rounding trap: just above a float midpoint.
        float f = 1.0f;
        double mid = 1.0 + Math.ScaleB(1, -24);
        Assert.Equal(MathF.BitIncrement(f), (float)new QuadDouble(mid, Math.ScaleB(1, -100)));
        Assert.Equal(f, (float)new QuadDouble(mid, -Math.ScaleB(1, -100)));
        Assert.Equal((Half)1.5, (Half)new QuadDouble(1.5));
        QuadDouble fromF = 0.1f, fromH = (Half)0.1;
        Assert.Equal((double)0.1f, fromF.X0);
        Assert.Equal((double)(Half)0.1, fromH.X0);
    }

    [Fact]
    public void DecimalBigIntRationalBigFloatConversions()
    {
        decimal d = 1234567890.123456789012345678m;
        QuadDouble q = (QuadDouble)d;
        BigFloat exactD = BigFloat.FromDecimal(d, 400);
        Assert.True(ErrBits(q, exactD) < -210);
        Assert.Equal(d, (decimal)q);
        Assert.Equal(decimal.MaxValue, (decimal)(QuadDouble)decimal.MaxValue);
        Assert.Equal(decimal.MinValue, (decimal)(QuadDouble)decimal.MinValue);
        Assert.Throws<OverflowException>(() => (decimal)new QuadDouble(1e30));
        Assert.True(QuadDouble.IsNegative((QuadDouble)(-0.0m)));

        BigInt big = BigInt.Parse("123456789012345678901234567890123456789012345678901234567890");
        QuadDouble qb = (QuadDouble)big;
        Assert.True(ErrBits(qb, BigFloat.FromBigInt(big, 400)) < -210);
        Assert.True(QuadDouble.IsPositiveInfinity((QuadDouble)(BigInt.One << 2000)));
        Assert.Equal(new BigInt(-12), (BigInt)new QuadDouble(-12.9));
        Assert.Equal((BigInteger)12, (BigInteger)new QuadDouble(12.9));

        BigRational third = new(BigInt.One, new BigInt(3));
        QuadDouble qt = (QuadDouble)third;
        Assert.Equal(QuadDouble.One / 3, qt);
        BigRational back = (BigRational)new QuadDouble(0.75);
        Assert.Equal(new BigRational(new BigInt(3), new BigInt(4)), back);

        // Exact round trip through BigFloat, including gap values.
        QuadDouble gap = new(1, Math.ScaleB(1, -400));
        BigFloat rg = gap;
        Assert.True(rg == BigFloat.Add(BigFloat.One, BigFloat.One.ShiftRight(400), 500));
        Assert.Equal(gap, (QuadDouble)rg);
        BigFloat pi = BigFloat.Pi(1000);
        Assert.True(ErrBits((QuadDouble)pi, pi) < -212);
    }

    // ---------------------------------------------------------------
    // Formatting and parsing
    // ---------------------------------------------------------------

    [Fact]
    public void ToStringParseRoundTrip()
    {
        var rng = new Random(77);
        for (int i = 0; i < 1000; i++)
        {
            QuadDouble q = RandomQd(rng, -1000, 1000);
            if (i % 3 == 0)
            {
                q = q / RandomQd(rng, -10, 10);
            }
            else if (i % 3 == 1 && i % 2 == 0)
            {
                // Components separated by a gap (spanning far more than 212 bits).
                q = new QuadDouble(q.X0, q.X0 * Math.ScaleB(1 + rng.NextDouble(), -rng.Next(60, 400)));
            }
            string s = q.ToString();
            QuadDouble back = QuadDouble.Parse(s, CultureInfo.InvariantCulture);
            Assert.True(q == back, $"{s}: {q.X0} {q.X1} {q.X2} {q.X3} vs {back.X0} {back.X1} {back.X2} {back.X3}");
            Assert.True(i % 6 == 4 || s.Split('e')[0].TrimStart('-').Replace(".", "").TrimStart('0').Length <= 75, s);
        }
        string pi = QuadDouble.Pi.ToString();
        Assert.StartsWith("3.14159265358979323846264338327950288419716939937510582097494459", pi);
        Assert.Equal("NaN", QuadDouble.NaN.ToString());
        Assert.Equal("-Infinity", QuadDouble.NegativeInfinity.ToString());
        Assert.Equal("-0", QuadDouble.NegativeZero.ToString());
        Assert.Equal("0.1", new QuadDouble(1).Divide10().ToString("G30", null)[..3]);
        QuadDouble gap = new(1, Math.ScaleB(1, -400));
        Assert.Equal(gap, QuadDouble.Parse(gap.ToString()));
    }

    [Fact]
    public void ParseVariants()
    {
        Assert.Equal(new QuadDouble(1.5), QuadDouble.Parse("1.5"));
        Assert.Equal(new QuadDouble(-1250), QuadDouble.Parse("-1.25e3"));
        Assert.True(QuadDouble.IsNaN(QuadDouble.Parse("NaN")));
        Assert.True(QuadDouble.IsNegativeInfinity(QuadDouble.Parse("-Infinity")));
        Assert.True(QuadDouble.IsNegative(QuadDouble.Parse("-0")));
        Assert.False(QuadDouble.TryParse("abc", null, out _));
        QuadDouble tenth = QuadDouble.Parse("0.1");
        Assert.True(ErrBits(tenth, BigFloat.Divide(BigFloat.One, BigFloat.FromDouble(10, P), P)) < -211);
        Assert.Equal(tenth, QuadDouble.Parse("0,1", NumberStyles.Float, new CultureInfo("de-DE")));
        Assert.Equal(new QuadDouble(1234567.5), QuadDouble.Parse("1,234,567.5", NumberStyles.Float | NumberStyles.AllowThousands, CultureInfo.InvariantCulture));
        Assert.Equal(tenth, QuadDouble.Parse("0.1"u8, null));
        Assert.True(QuadDouble.TryParse(" 0.1 "u8, NumberStyles.Float, CultureInfo.InvariantCulture, out QuadDouble u8));
        Assert.Equal(tenth, u8);
        Assert.Equal(tenth, QuadDouble.Parse("0.1".AsSpan(), null));

        Span<char> buf = stackalloc char[128];
        Assert.True(tenth.TryFormat(buf, out int n, default, null));
        Assert.Equal(tenth.ToString(), buf[..n].ToString());
        Assert.False(tenth.TryFormat(buf[..5], out _, default, null));
        Span<byte> ubuf = stackalloc byte[128];
        Assert.True(tenth.TryFormat(ubuf, out int nb, default, null));
        Assert.Equal(tenth.ToString(), Encoding.UTF8.GetString(ubuf[..nb]));
        Assert.Equal("0.1", tenth.ToString("G1", CultureInfo.InvariantCulture));
    }

    // ---------------------------------------------------------------
    // Transcendental functions vs BigFloat
    // ---------------------------------------------------------------

    private delegate QuadDouble QdFunc(QuadDouble x);
    private delegate BigFloat BigFloatFunc(BigFloat x, int precision);

    private double MaxError(int seed, int count, Func<Random, QuadDouble> gen, QdFunc f, BigFloatFunc g)
    {
        var rng = new Random(seed);
        double worst = -1000;
        for (int i = 0; i < count; i++)
        {
            QuadDouble x = gen(rng);
            BigFloat exact = g(x, P);
            if (exact.IsZero || !exact.IsFinite || exact.BinaryExponent > 1024 || exact.BinaryExponent < -850)
            {
                continue;
            }
            double e = ErrBits(f(x), exact);
            if (e > worst)
            {
                worst = e;
            }
        }
        return worst;
    }

    private const int N = 400;

    [Fact]
    public void ExpFamily()
    {
        Check("exp", MaxError(1, N, r => RandomQd(r, -12, 8), QuadDouble.Exp, (x, p) => BigFloat.Exp(x, p)), -205);
        Check("exp tiny", MaxError(2, N, r => RandomQd(r, -200, -20), QuadDouble.Exp, (x, p) => BigFloat.Exp(x, p)), -205);
        Check("expm1", MaxError(3, N, r => RandomQd(r, -120, 6), QuadDouble.ExpM1, (x, p) => BigFloat.ExpM1(x, p)), -205);
        Check("exp2", MaxError(4, N, r => RandomQd(r, -12, 10), QuadDouble.Exp2, (x, p) => BigFloat.Exp2(x, p)), -205);
        Check("exp10", MaxError(5, N, r => RandomQd(r, -12, 8), QuadDouble.Exp10, (x, p) => BigFloat.Exp10(x, p)), -205);
        Check("exp2m1", MaxError(6, N, r => RandomQd(r, -60, 3), QuadDouble.Exp2M1, (x, p) => BigFloat.Exp2M1(x, p)), -205);
        Check("exp10m1", MaxError(7, N, r => RandomQd(r, -60, 3), QuadDouble.Exp10M1, (x, p) => BigFloat.Exp10M1(x, p)), -205);
        Assert.Equal(QuadDouble.One, QuadDouble.Exp(0));
        Assert.True(QuadDouble.IsPositiveInfinity(QuadDouble.Exp(710)));
        Assert.Equal(QuadDouble.Zero, QuadDouble.Exp(-746));
        Assert.Equal(new QuadDouble(1024), QuadDouble.Exp2(10));
        Assert.Equal(new QuadDouble(0.125), QuadDouble.Exp2(-3));
        Assert.True(ErrBits(QuadDouble.Exp(709.7), BigFloat.Exp(BigFloat.FromDouble(709.7, P), P)) < -205);
        Assert.Equal(QuadDouble.NegativeOne, QuadDouble.ExpM1(-1000));
    }

    [Fact]
    public void LogFamily()
    {
        Check("log", MaxError(11, N, r => RandomQd(r, -1000, 1000, false), QuadDouble.Log, (x, p) => BigFloat.Log(x, p)), -205);
        Check("log near 1", MaxError(12, N, r => 1 + RandomQd(r, -150, -2), QuadDouble.Log, (x, p) => BigFloat.Log(x, p)), -205);
        Check("log1p", MaxError(13, N, r => RandomQd(r, -150, 4) is var v && v < -0.9 ? -v : v, QuadDouble.LogP1, (x, p) => BigFloat.LogP1(x, p)), -205);
        Check("log2", MaxError(14, N, r => RandomQd(r, -1000, 1000, false), QuadDouble.Log2, (x, p) => BigFloat.Log2(x, p)), -205);
        Check("log10", MaxError(15, N, r => RandomQd(r, -1000, 1000, false), QuadDouble.Log10, (x, p) => BigFloat.Log10(x, p)), -205);
        Assert.Equal(QuadDouble.Zero, QuadDouble.Log(1));
        Assert.Equal(new QuadDouble(-30), QuadDouble.Log2(Math.ScaleB(1, -30)));
        Assert.True(QuadDouble.IsNaN(QuadDouble.Log(-1)));
        Assert.True(QuadDouble.IsNegativeInfinity(QuadDouble.Log(0)));
        Assert.True(QuadDouble.IsNegativeInfinity(QuadDouble.LogP1(-1)));
        Assert.True(ErrBits(QuadDouble.Log(QuadDouble.E), BigFloat.One) < -207);
    }

    [Fact]
    public void PowAndRoots()
    {
        var rng = new Random(21);
        double worst = -1000, worstI = -1000, worstR = -1000, worstH = -1000;
        for (int i = 0; i < 300; i++)
        {
            QuadDouble x = RandomQd(rng, -6, 6, false), y = RandomQd(rng, -4, 3);
            worst = Math.Max(worst, ErrBits(QuadDouble.Pow(x, y), BigFloat.Pow(x, y, P)));
            int n = rng.Next(-40, 41);
            QuadDouble b = RandomQd(rng, -3, 3);
            worstI = Math.Max(worstI, ErrBits(QuadDouble.Pow(b, n), BigFloat.Pow(b, (long)n, P)));
            int k = rng.Next(2, 30) * (rng.Next(4) == 0 ? -1 : 1);
            QuadDouble c = RandomQd(rng, -300, 300, false);
            worstR = Math.Max(worstR, ErrBits(QuadDouble.RootN(c, k), BigFloat.RootN(c, (long)k, P)));
            QuadDouble h1 = RandomQd(rng, -500, 500), h2 = RandomQd(rng, -500, 500);
            worstH = Math.Max(worstH, ErrBits(QuadDouble.Hypot(h1, h2), BigFloat.Hypot(h1, h2, P)));
        }
        Check("pow", worst, -205);
        Check("pow int", worstI, -205);
        Check("rootn", worstR, -205);
        Check("hypot", worstH, -207);
        Assert.Equal(new QuadDouble(1024), QuadDouble.Pow(2, 10));
        Assert.Equal(new QuadDouble(-8), QuadDouble.Pow(-2, new QuadDouble(3)));
        Assert.True(QuadDouble.IsNaN(QuadDouble.Pow(-2, new QuadDouble(0.5))));
        Assert.Equal(QuadDouble.One, QuadDouble.Pow(QuadDouble.NaN, QuadDouble.Zero));
        Assert.True(QuadDouble.IsPositiveInfinity(QuadDouble.Pow(0, new QuadDouble(-2))));
        Assert.True(QuadDouble.IsNegativeInfinity(QuadDouble.Pow(QuadDouble.NegativeZero, new QuadDouble(-3))));
        Assert.True(ErrBits(QuadDouble.RootN(-32, 5), BigFloat.FromDouble(-2, 64)) < -210);
        Assert.True(QuadDouble.IsNaN(QuadDouble.RootN(-16, 4)));
        Assert.Equal(new QuadDouble(5), QuadDouble.Hypot(3, -4));
    }

    [Fact]
    public void TrigFamily()
    {
        Check("sin", MaxError(31, N, r => RandomQd(r, -30, 6), QuadDouble.Sin, (x, p) => BigFloat.Sin(x, p)), -205);
        Check("cos", MaxError(32, N, r => RandomQd(r, -30, 6), QuadDouble.Cos, (x, p) => BigFloat.Cos(x, p)), -205);
        Check("tan", MaxError(33, N, r => RandomQd(r, -30, 6), QuadDouble.Tan, (x, p) => BigFloat.Tan(x, p)), -204);
        Check("sin large", MaxError(34, 100, r => RandomQd(r, 20, 900), QuadDouble.Sin, (x, p) => BigFloat.Sin(x, p)), -205);
        Check("cos large", MaxError(35, 100, r => RandomQd(r, 20, 900), QuadDouble.Cos, (x, p) => BigFloat.Cos(x, p)), -205);
        Check("asin", MaxError(36, N, r => RandomQd(r, -40, -1), QuadDouble.Asin, (x, p) => BigFloat.Asin(x, p)), -205);
        Check("acos", MaxError(37, N, r => RandomQd(r, -40, -1), QuadDouble.Acos, (x, p) => BigFloat.Acos(x, p)), -205);
        Check("atan", MaxError(38, N, r => RandomQd(r, -40, 40), QuadDouble.Atan, (x, p) => BigFloat.Atan(x, p)), -205);
        Check("sinpi", MaxError(39, N, r => RandomQd(r, -30, 10), QuadDouble.SinPi, (x, p) => BigFloat.SinPi(x, p)), -205);
        Check("cospi", MaxError(40, N, r => RandomQd(r, -30, 10), QuadDouble.CosPi, (x, p) => BigFloat.CosPi(x, p)), -205);
        Check("tanpi", MaxError(41, N, r => RandomQd(r, -30, 10), QuadDouble.TanPi, (x, p) => BigFloat.TanPi(x, p)), -204);
        Check("asinpi", MaxError(42, N, r => RandomQd(r, -40, -1), QuadDouble.AsinPi, (x, p) => BigFloat.AsinPi(x, p)), -205);
        Check("acospi", MaxError(43, N, r => RandomQd(r, -40, -1), QuadDouble.AcosPi, (x, p) => BigFloat.AcosPi(x, p)), -205);
        Check("atanpi", MaxError(44, N, r => RandomQd(r, -40, 40), QuadDouble.AtanPi, (x, p) => BigFloat.AtanPi(x, p)), -205);

        // Every table entry: sin/cos at j·π/256 + small offsets.
        double worst = -1000;
        var rng = new Random(45);
        for (int j = -520; j <= 520; j++)
        {
            QuadDouble x = QuadDouble.Pi * j / 256 + RandomQd(rng, -12, -9);
            worst = Math.Max(worst, ErrBits(QuadDouble.Sin(x), BigFloat.Sin(x, P)));
            worst = Math.Max(worst, ErrBits(QuadDouble.Cos(x), BigFloat.Cos(x, P)));
        }
        Check("sin/cos table sweep", worst, -205);

        // Near multiples of π the reduction must stay relatively accurate.
        Check("sin(k·pi)", MaxError(46, 50, r => QuadDouble.Pi * r.Next(1, 100000), QuadDouble.Sin, (x, p) => BigFloat.Sin(x, p)), -205);

        double worst2 = -1000;
        for (int i = 0; i < N; i++)
        {
            QuadDouble y = RandomQd(rng, -50, 50), x = RandomQd(rng, -50, 50);
            worst2 = Math.Max(worst2, ErrBits(QuadDouble.Atan2(y, x), BigFloat.Atan2(y, x, P)));
            worst2 = Math.Max(worst2, ErrBits(QuadDouble.Atan2Pi(y, x), BigFloat.Atan2Pi(y, x, P)));
        }
        Check("atan2/atan2pi", worst2, -205);

        (QuadDouble s, QuadDouble c) = QuadDouble.SinCos(1);
        Assert.Equal(QuadDouble.Sin(1), s);
        Assert.Equal(QuadDouble.Cos(1), c);
        Assert.Equal(QuadDouble.One, QuadDouble.SinPi(0.5));
        Assert.Equal(QuadDouble.Zero, QuadDouble.CosPi(0.5));
        Assert.Equal(QuadDouble.NegativeOne, QuadDouble.CosPi(1));
        Assert.True(QuadDouble.IsNegative(QuadDouble.Sin(QuadDouble.NegativeZero)));
        Assert.True(QuadDouble.IsNaN(QuadDouble.Sin(QuadDouble.PositiveInfinity)));
        Assert.Equal(new QuadDouble(0.5), QuadDouble.AsinPi(1));
        Assert.Equal(new QuadDouble(0.25), QuadDouble.AtanPi(1));
        Assert.Equal(QuadDouble.One, QuadDouble.AcosPi(-1));
        Assert.Equal(new QuadDouble(-0.75), QuadDouble.Atan2Pi(-1, -1));
        Assert.Equal(QuadDouble.Pi, QuadDouble.Atan2(0, -1));
        Assert.True(QuadDouble.IsNaN(QuadDouble.Asin(1.5)));
        Assert.True(ErrBits(QuadDouble.DegreesToRadians(180), BigFloat.Pi(P)) < -210);
        Assert.True(ErrBits(QuadDouble.RadiansToDegrees(QuadDouble.Pi), BigFloat.FromDouble(180, P)) < -209);
        (QuadDouble sp, QuadDouble cp) = QuadDouble.SinCosPi(0.25);
        Assert.True(ErrBits(sp, BigFloat.Sqrt(BigFloat.Half, P)) < -207);
        Assert.Equal(sp, cp);
    }

    [Fact]
    public void HyperbolicFamily()
    {
        Check("sinh", MaxError(51, N, r => RandomQd(r, -100, 8), QuadDouble.Sinh, (x, p) => BigFloat.Sinh(x, p)), -205);
        Check("cosh", MaxError(52, N, r => RandomQd(r, -100, 8), QuadDouble.Cosh, (x, p) => BigFloat.Cosh(x, p)), -205);
        Check("tanh", MaxError(53, N, r => RandomQd(r, -100, 6), QuadDouble.Tanh, (x, p) => BigFloat.Tanh(x, p)), -205);
        Check("asinh", MaxError(54, N, r => RandomQd(r, -100, 200), QuadDouble.Asinh, (x, p) => BigFloat.Asinh(x, p)), -205);
        Check("acosh", MaxError(55, N, r => 1 + RandomQd(r, -100, 200, false), QuadDouble.Acosh, (x, p) => BigFloat.Acosh(x, p)), -205);
        Check("atanh", MaxError(56, N, r => RandomQd(r, -100, -1), QuadDouble.Atanh, (x, p) => BigFloat.Atanh(x, p)), -205);
        Assert.True(ErrBits(QuadDouble.Sinh(710), BigFloat.Sinh(BigFloat.FromDouble(710, P), P)) < -205);
        Assert.True(QuadDouble.IsPositiveInfinity(QuadDouble.Cosh(712)));
        Assert.Equal(QuadDouble.One, QuadDouble.Tanh(100));
        Assert.True(QuadDouble.IsPositiveInfinity(QuadDouble.Atanh(1)));
        Assert.True(QuadDouble.IsNaN(QuadDouble.Acosh(0.5)));
    }

    // ---------------------------------------------------------------
    // Generic math
    // ---------------------------------------------------------------

    private static T Horner<T>(T x, params T[] c) where T : IFloatingPoint<T>
    {
        T r = T.Zero;
        foreach (T k in c)
        {
            r = r * x + k;
        }
        return r;
    }

    private static T GenericSqrtExp<T>(T x) where T : IRootFunctions<T>, IExponentialFunctions<T>, ITrigonometricFunctions<T> =>
        T.Sqrt(x) + T.Exp(x) + T.Sin(x);

    [Fact]
    public void GenericMath()
    {
        QuadDouble x = QuadDouble.One / 7;
        Assert.Equal(x * x + 2 * x + 3, Horner(x, (QuadDouble)1, 2, 3));
        Assert.Equal(QuadDouble.Sqrt(x) + QuadDouble.Exp(x) + QuadDouble.Sin(x), GenericSqrtExp(x));
        Assert.Equal(QuadDouble.Pi, GetPi<QuadDouble>());
        Assert.Equal(QuadDouble.MaxValue, GetMax<QuadDouble>());
        Assert.Equal(new QuadDouble(2.5), QuadDouble.CreateChecked(2.5m));
        Assert.Equal(2.5, double.CreateChecked(new QuadDouble(2.5)));
        Assert.Equal(new QuadDouble(3), QuadDouble.CreateChecked((BigInt)3));
        Assert.Equal(new QuadDouble(3), QuadDouble.CreateChecked(BigFloat.FromDouble(3, 64)));
        Assert.Equal(new BigInt(3), BigInt.CreateChecked(new QuadDouble(3.9)));
        Assert.Equal(QuadDouble.NegativeOne, GetNegOne<QuadDouble>());
        Assert.Equal(new QuadDouble(4), QuadDouble.Clamp(10, 1, 4));
        IFloatingPoint<QuadDouble> fp = new QuadDouble(6);
        Assert.Equal(2, fp.GetExponentShortestBitLength());
    }

    private static T GetPi<T>() where T : IFloatingPointConstants<T> => T.Pi;
    private static T GetMax<T>() where T : IMinMaxValue<T> => T.MaxValue;
    private static T GetNegOne<T>() where T : ISignedNumber<T> => T.NegativeOne;
}

internal static class QuadDoubleTestExtensions
{
    public static QuadDouble Divide10(this QuadDouble q) => q / 10;
}
