using System.Globalization;
using System.Numerics;
using System.Text;
using Xunit;
using Xunit.Abstractions;

namespace Sdcb.Arithmetic.Tests;

public class DoubleDoubleTests(ITestOutputHelper output)
{
    private const int P = 300; // reference precision for Real

    // ---------------------------------------------------------------
    // Helpers
    // ---------------------------------------------------------------

    /// <summary>The exact value (independent of Real.DefaultPrecision, which other tests mutate).</summary>
    private static Real R(DoubleDouble x) => (Real)x;

    /// <summary>log2 of the relative error of got against exact (-inf when exact).</summary>
    private static double Log2RelErr(DoubleDouble got, Real exact)
    {
        Real diff = Real.Subtract(R(got), exact, P);
        if (diff.IsZero)
        {
            return double.NegativeInfinity;
        }
        return Math.Log2(Math.Abs(Real.Divide(diff, exact, 64).ToDouble()));
    }

    private static Real Rd(double d) => Real.FromDouble(d, P);

    /// <summary>A random double-double with a full random low word and binary exponent in [minExp, maxExp].</summary>
    private static DoubleDouble RandomDD(Random rng, int minExp, int maxExp, bool allowNegative = true)
    {
        double m = 1.0 + rng.NextDouble();
        double hi = Math.ScaleB(m, rng.Next(minExp, maxExp + 1));
        if (allowNegative && rng.Next(2) == 0)
        {
            hi = -hi;
        }
        double lo = hi * (rng.NextDouble() - 0.5) * Math.ScaleB(1.0, -52);
        return new DoubleDouble(hi, lo);
    }

    private void AssertBound(string name, IEnumerable<(DoubleDouble Got, Real Exact, string Input)> cases, double log2Bound)
    {
        double worst = double.NegativeInfinity;
        string worstInput = "";
        int n = 0;
        foreach ((DoubleDouble got, Real exact, string input) in cases)
        {
            n++;
            double e = Log2RelErr(got, exact);
            if (e > worst)
            {
                worst = e;
                worstInput = input;
            }
        }
        output.WriteLine($"{name}: {n} cases, max relative error 2^{worst:F2} (bound 2^{log2Bound})");
        Assert.True(worst <= log2Bound, $"{name}: relative error 2^{worst:F2} exceeds 2^{log2Bound} at {worstInput}");
    }

    private static string S(DoubleDouble x) => $"({x.Hi:R}, {x.Lo:R})";

    // ---------------------------------------------------------------
    // Representation and constants
    // ---------------------------------------------------------------

    [Fact]
    public void ConstructorNormalizes()
    {
        DoubleDouble x = new(1.0, 1.0);
        Assert.Equal(2.0, x.Hi);
        Assert.Equal(0.0, x.Lo);
        DoubleDouble y = new(1e-20, 1.0);
        Assert.Equal(1.0, y.Hi);
        Assert.Equal(1e-20, y.Lo);
        DoubleDouble z = new(double.MaxValue, double.MaxValue);
        Assert.True(DoubleDouble.IsPositiveInfinity(z));
        Assert.Equal(0.0, z.Lo);
        Assert.True(DoubleDouble.IsNaN(new DoubleDouble(double.PositiveInfinity, double.NegativeInfinity)));
    }

    public static TheoryData<string> ConstantNames => ["Pi", "Tau", "E", "Ln2", "Ln10"];

    [Theory]
    [MemberData(nameof(ConstantNames))]
    public void ConstantsMatchReal(string name)
    {
        (DoubleDouble dd, Real exact) = name switch
        {
            "Pi" => (DoubleDouble.Pi, Real.Pi(P)),
            "Tau" => (DoubleDouble.Tau, Real.Pi(P).ShiftLeft(1)),
            "E" => (DoubleDouble.E, Real.E(P)),
            "Ln2" => (DoubleDouble.Ln2, Real.Ln2(P)),
            "Ln10" => (DoubleDouble.Ln10, Real.Log(Rd(10), P)),
            _ => throw new ArgumentException(name),
        };
        double hi = exact.ToDouble();
        double lo = Real.Subtract(exact, Rd(hi), 53).ToDouble();
        Assert.Equal(hi, dd.Hi);
        Assert.Equal(lo, dd.Lo);
        Assert.Equal(DoubleDouble.FromReal(exact), dd);
        Assert.True(Log2RelErr(dd, exact) < -106);
    }

    [Fact]
    public void LimitConstants()
    {
        DoubleDouble max = DoubleDouble.MaxValue;
        Assert.Equal(double.MaxValue, max.Hi);
        Assert.Equal(double.MaxValue, max.Hi + max.Lo); // normalized
        Assert.True(double.IsInfinity(max.Hi + Math.BitIncrement(max.Lo))); // and the largest such low word
        Assert.Equal(-max, DoubleDouble.MinValue);
        Assert.Equal(Math.ScaleB(1.0, -104), DoubleDouble.Epsilon.Hi);
        Assert.True(DoubleDouble.One + DoubleDouble.Epsilon > DoubleDouble.One);
        Assert.Equal(1.0, (double)DoubleDouble.One);
        Assert.Equal(-1.0, (double)DoubleDouble.NegativeOne);
        Assert.True(DoubleDouble.IsZero(DoubleDouble.Zero));
        Assert.True(DoubleDouble.IsNaN(DoubleDouble.NaN));
        Assert.True(DoubleDouble.IsNegativeInfinity(DoubleDouble.NegativeInfinity));
    }

    // ---------------------------------------------------------------
    // Basic arithmetic vs Real
    // ---------------------------------------------------------------

    private static IEnumerable<(DoubleDouble A, DoubleDouble B)> RandomPairs(int seed, int count, int minExp = -60, int maxExp = 60)
    {
        Random rng = new(seed);
        for (int i = 0; i < count; i++)
        {
            DoubleDouble a = RandomDD(rng, minExp, maxExp);
            DoubleDouble b = rng.Next(8) switch
            {
                // Cancellation: b close to -a
                0 => -a + RandomDD(rng, -40, -20) * a,
                1 => -a + new DoubleDouble(0.0, a.Lo * rng.NextDouble()),
                _ => RandomDD(rng, minExp, maxExp),
            };
            yield return (a, b);
        }
    }

    [Fact]
    public void AdditionAndSubtraction()
    {
        AssertBound("add", RandomPairs(1, 4000).Select(p => (p.A + p.B, Real.Add(R(p.A), R(p.B), P), S(p.A) + S(p.B))), -104);
        AssertBound("sub", RandomPairs(2, 4000).Select(p => (p.A - p.B, Real.Subtract(R(p.A), R(p.B), P), S(p.A) + S(p.B))), -104);
        AssertBound("add-double", RandomPairs(3, 4000).Select(p => (p.A + p.B.Hi, Real.Add(R(p.A), Rd(p.B.Hi), P), S(p.A) + p.B.Hi)), -104);
        AssertBound("double-sub", RandomPairs(4, 4000).Select(p => (p.B.Hi - p.A, Real.Subtract(Rd(p.B.Hi), R(p.A), P), S(p.A) + p.B.Hi)), -104);
    }

    [Fact]
    public void AdditionCancellationIsExact()
    {
        DoubleDouble a = new(1.0, 1e-20);
        Assert.Equal(new DoubleDouble(1e-20, 0), a - 1.0);
        Assert.Equal(DoubleDouble.Zero, a - a);
        DoubleDouble b = DoubleDouble.Pi - 3.0;
        Assert.True(Log2RelErr(b, Real.Subtract(R(DoubleDouble.Pi), Rd(3), P)) == double.NegativeInfinity);
    }

    [Fact]
    public void Multiplication()
    {
        AssertBound("mul", RandomPairs(5, 4000, -300, 300).Select(p => (p.A * p.B, Real.Multiply(R(p.A), R(p.B), P), S(p.A) + S(p.B))), -104);
        AssertBound("mul-double", RandomPairs(6, 4000, -300, 300).Select(p => (p.A * p.B.Hi, Real.Multiply(R(p.A), Rd(p.B.Hi), P), S(p.A) + p.B.Hi)), -104);
        AssertBound("square", RandomPairs(7, 4000, -300, 300).Select(p => (DoubleDouble.Square(p.A), Real.Multiply(R(p.A), R(p.A), P), S(p.A))), -104);
    }

    [Fact]
    public void Division()
    {
        // DWDivDW2 is proven to 15u² + 56u³ (< 2^-102.09); DWDivFP3 to 3u².
        AssertBound("div", RandomPairs(8, 4000, -300, 300).Select(p => (p.A / p.B, Real.Divide(R(p.A), R(p.B), P), S(p.A) + S(p.B))), -102);
        AssertBound("div-double", RandomPairs(9, 4000, -300, 300).Select(p => (p.A / p.B.Hi, Real.Divide(R(p.A), Rd(p.B.Hi), P), S(p.A) + p.B.Hi)), -104);
        AssertBound("double-div", RandomPairs(10, 4000, -300, 300).Select(p => (p.B.Hi / p.A, Real.Divide(Rd(p.B.Hi), R(p.A), P), S(p.A) + p.B.Hi)), -102);
    }

    [Fact]
    public void SquareRootAndCubeRoot()
    {
        Random rng = new(11);
        DoubleDouble[] xs = Enumerable.Range(0, 3000).Select(_ => DoubleDouble.Abs(RandomDD(rng, -1000, 1000))).ToArray();
        AssertBound("sqrt", xs.Select(x => (DoubleDouble.Sqrt(x), Real.Sqrt(R(x), P), S(x))), -104);
        AssertBound("cbrt", xs.Select(x => (DoubleDouble.Cbrt(x), Real.Cbrt(R(x), P), S(x))), -103);
        AssertBound("cbrt-neg", xs.Select(x => (DoubleDouble.Cbrt(-x), Real.Cbrt(-R(x), P), S(x))), -103);
        AssertBound("rootn", xs.Take(500).SelectMany(x => new[] { 3, 4, 5, 7, 10, 33 }.Select(n =>
            (DoubleDouble.RootN(x, n), Real.RootN(R(x), (long)n, P), S(x) + " n=" + n))), -102);
        Assert.Equal(DoubleDouble.Sqrt(new DoubleDouble(4.0, 0.0)), (DoubleDouble)2.0);
        Assert.True(DoubleDouble.IsNaN(DoubleDouble.Sqrt(-1.0)));
        Assert.True(double.IsNegative(DoubleDouble.Sqrt(-0.0).Hi));
        Assert.True(DoubleDouble.IsPositiveInfinity(DoubleDouble.Sqrt(double.PositiveInfinity)));
        Assert.True(Log2RelErr(DoubleDouble.Sqrt(double.Epsilon), Real.Sqrt(Rd(double.Epsilon), P)) < -50);
    }

    [Fact]
    public void SpecialValueArithmetic()
    {
        DoubleDouble inf = DoubleDouble.PositiveInfinity, max = DoubleDouble.MaxValue;
        Assert.Equal(inf, inf + 1.0);
        Assert.Equal(inf, inf + max);
        Assert.Equal(0.0, (inf + max).Lo);
        Assert.True(DoubleDouble.IsNaN(inf - inf));
        Assert.Equal(inf, max + max);
        Assert.Equal(inf, max * 2.0);
        Assert.Equal(DoubleDouble.NegativeInfinity, max * -2.0);
        Assert.Equal(inf, inf * inf);
        Assert.True(DoubleDouble.IsNaN(inf * 0.0));
        Assert.Equal(inf, DoubleDouble.One / 0.0);
        Assert.Equal(DoubleDouble.NegativeInfinity, DoubleDouble.One / -0.0);
        Assert.True(DoubleDouble.IsNaN(DoubleDouble.Zero / DoubleDouble.Zero));
        Assert.Equal(DoubleDouble.Zero, DoubleDouble.One / inf);
        Assert.Equal(inf, max / 0.5);
        Assert.Equal(0.0, (max / 0.5).Lo);
        Assert.Equal(max, max * 1.0);
        Assert.Equal(max, max / 1.0);
        // Quotient that overflows only in the leading estimate is still found when it is finite.
        DoubleDouble small = Math.ScaleB(1.0, -1000), big = Math.ScaleB(1.0, 1000);
        Assert.Equal((DoubleDouble)Math.ScaleB(1.0, -1000), small / DoubleDouble.One);
        Assert.Equal(0.0, (small / big).Hi);
        DoubleDouble nan = DoubleDouble.NaN;
        Assert.True(DoubleDouble.IsNaN(nan + 1.0) && DoubleDouble.IsNaN(nan * 2.0) && DoubleDouble.IsNaN(nan / 2.0));
    }

    [Fact]
    public void Remainder()
    {
        Random rng = new(12);
        for (int i = 0; i < 300; i++)
        {
            DoubleDouble a = RandomDD(rng, -20, 60), b = RandomDD(rng, -20, 20);
            Real expected = Real.TruncatedRemainder(R(a), R(b), P);
            Assert.Equal(DoubleDouble.FromReal(expected), a % b);
        }
        Assert.Equal((DoubleDouble)1.5, (DoubleDouble)7.5 % 3.0);
        Assert.Equal((DoubleDouble)(-1.5), (DoubleDouble)(-7.5) % 3.0);
        Assert.True(DoubleDouble.IsNaN(DoubleDouble.One % 0.0));
        Assert.Equal((DoubleDouble)2.0, (DoubleDouble)2.0 % double.PositiveInfinity);
    }

    // ---------------------------------------------------------------
    // Comparison, equality, classification
    // ---------------------------------------------------------------

    [Fact]
    public void ComparisonAndEquality()
    {
        DoubleDouble a = new(1.0, 1e-20), b = new(1.0, -1e-20);
        Assert.True(a > b && b < a && a >= b && b <= a && a != b);
        Assert.True(a > 1.0 && b < 1.0);
        Assert.Equal(1, a.CompareTo(b));
        Assert.Equal(0, a.CompareTo((object)a));
        Assert.Equal(1, a.CompareTo(null));
        DoubleDouble nan = DoubleDouble.NaN, nan2 = DoubleDouble.NaN;
        Assert.False(nan == nan2);
        Assert.True(nan != nan2);
        Assert.True(nan.Equals(nan));
        Assert.False(nan < a || nan > a || nan <= a || nan >= a);
        Assert.Equal(-1, nan.CompareTo(a));
        DoubleDouble pz = 0.0, nz = -0.0;
        Assert.True(pz == nz);
        Assert.True(pz.Equals(nz));
        Assert.Equal(pz.GetHashCode(), nz.GetHashCode());
        Assert.Equal(new DoubleDouble(3.0, 1e-17).GetHashCode(), new DoubleDouble(3.0, 1e-17).GetHashCode());
        Assert.True(DoubleDouble.IsNegative(nz) && !DoubleDouble.IsNegative(pz));
    }

    [Fact]
    public void Classification()
    {
        DoubleDouble big = (DoubleDouble)long.MaxValue;
        Assert.True(DoubleDouble.IsInteger(big));
        Assert.True(DoubleDouble.IsOddInteger(big));
        Assert.True(DoubleDouble.IsEvenInteger(big + 1L));
        Assert.True(DoubleDouble.IsInteger(new DoubleDouble(Math.ScaleB(1.0, 60), 3.0)));
        Assert.False(DoubleDouble.IsInteger(new DoubleDouble(Math.ScaleB(1.0, 60), 0.5)));
        Assert.False(DoubleDouble.IsInteger(DoubleDouble.PositiveInfinity));
        Assert.True(DoubleDouble.IsEvenInteger(DoubleDouble.Zero));
        Assert.True(DoubleDouble.IsPow2(8.0));
        Assert.False(DoubleDouble.IsPow2(new DoubleDouble(8.0, 1e-20)));
        Assert.True(DoubleDouble.IsNormal(1.0) && !DoubleDouble.IsNormal(0.0));
        Assert.True(DoubleDouble.IsSubnormal(double.Epsilon));
        Assert.True(DoubleDouble.IsFinite(DoubleDouble.MaxValue) && !DoubleDouble.IsFinite(DoubleDouble.NaN));
        Assert.Equal(3, DoubleDouble.ILogB(8.5));
        Assert.Equal(2, DoubleDouble.ILogB(new DoubleDouble(8.0, -1e-20)));
        Assert.Equal(-4, DoubleDouble.ILogB(new DoubleDouble(-0.125, 1e-30))); // magnitude just below 1/8
        Assert.Equal(int.MinValue, DoubleDouble.ILogB(0.0));
        Assert.Equal(-1, DoubleDouble.Sign(-3.0));
        Assert.Equal(0, DoubleDouble.Sign(-0.0));
        Assert.Throws<ArithmeticException>(() => DoubleDouble.Sign(DoubleDouble.NaN));
    }

    [Fact]
    public void MinMaxAndSign()
    {
        DoubleDouble a = new(1.0, 1e-20), b = new(1.0, -1e-20), nan = DoubleDouble.NaN;
        Assert.Equal(a, DoubleDouble.Max(a, b));
        Assert.Equal(b, DoubleDouble.Min(a, b));
        Assert.True(DoubleDouble.IsNaN(DoubleDouble.Max(a, nan)));
        Assert.Equal(a, DoubleDouble.MaxNumber(a, nan));
        Assert.Equal(a, DoubleDouble.MinNumber(nan, a));
        Assert.False(double.IsNegative(DoubleDouble.Max(-0.0, 0.0).Hi));
        Assert.True(double.IsNegative(DoubleDouble.Min(0.0, -0.0).Hi));
        Assert.Equal((DoubleDouble)(-3.0), DoubleDouble.MaxMagnitude(2.0, -3.0));
        Assert.Equal((DoubleDouble)2.0, DoubleDouble.MinMagnitude(2.0, -3.0));
        Assert.Equal((DoubleDouble)3.0, DoubleDouble.MaxMagnitude(-3.0, 3.0));
        Assert.Equal((DoubleDouble)(-3.0), DoubleDouble.MinMagnitude(-3.0, 3.0));
        Assert.Equal(-a, DoubleDouble.CopySign(a, -2.0));
        Assert.Equal(a, DoubleDouble.Abs(-a));
        Assert.Equal((DoubleDouble)2.0, DoubleDouble.Clamp(5.0, 1.0, 2.0));
        Assert.Throws<ArgumentException>(() => DoubleDouble.Clamp(5.0, 2.0, 1.0));
        DoubleDouble x = 1.5;
        x++;
        Assert.Equal((DoubleDouble)2.5, x);
        x--;
        x--;
        Assert.Equal((DoubleDouble)0.5, x);
    }

    // ---------------------------------------------------------------
    // Rounding
    // ---------------------------------------------------------------

    [Fact]
    public void Rounding()
    {
        DoubleDouble half = new(2.5, 0.0), above = new(2.5, 1e-20), below = new(2.5, -1e-20);
        Assert.Equal((DoubleDouble)2.0, DoubleDouble.Round(half));
        Assert.Equal((DoubleDouble)3.0, DoubleDouble.Round(half, MidpointRounding.AwayFromZero));
        Assert.Equal((DoubleDouble)3.0, DoubleDouble.Round(above));
        Assert.Equal((DoubleDouble)2.0, DoubleDouble.Round(below, MidpointRounding.AwayFromZero));
        Assert.Equal((DoubleDouble)(-2.0), DoubleDouble.Round(-half));
        Assert.Equal((DoubleDouble)(-3.0), DoubleDouble.Round(-half, MidpointRounding.AwayFromZero));
        Assert.Equal((DoubleDouble)2.0, DoubleDouble.Floor(above));
        Assert.Equal((DoubleDouble)3.0, DoubleDouble.Ceiling(above));
        Assert.Equal((DoubleDouble)(-3.0), DoubleDouble.Floor(-above));
        Assert.Equal((DoubleDouble)(-2.0), DoubleDouble.Truncate(-above));
        Assert.True(double.IsNegative(DoubleDouble.Truncate(-0.5).Hi));
        Assert.True(double.IsNegative(DoubleDouble.Round(-0.25).Hi));
        // Integer high word, fractional low word.
        DoubleDouble big = new(Math.ScaleB(1.0, 60), -0.5);
        Assert.Equal(new DoubleDouble(Math.ScaleB(1.0, 60), -1.0), DoubleDouble.Floor(big));
        Assert.Equal(new DoubleDouble(Math.ScaleB(1.0, 60), 0.0), DoubleDouble.Ceiling(big));
        Assert.Equal(new DoubleDouble(Math.ScaleB(1.0, 60), 0.0), DoubleDouble.Round(big)); // tie to even
        Assert.Equal(new DoubleDouble(Math.ScaleB(1.0, 60), -1.0), DoubleDouble.Round(big, MidpointRounding.ToZero));
        Assert.Equal(DoubleDouble.Parse("3.14159"), DoubleDouble.Round(DoubleDouble.Pi, 5));
        Assert.Equal(DoubleDouble.Parse("2.72"), DoubleDouble.Round(DoubleDouble.E, 2, MidpointRounding.AwayFromZero));
        Assert.True(DoubleDouble.IsPositiveInfinity(DoubleDouble.Round(DoubleDouble.PositiveInfinity)));
        Random rng = new(13);
        for (int i = 0; i < 500; i++)
        {
            DoubleDouble x = RandomDD(rng, -3, 80);
            Assert.Equal(DoubleDouble.FromReal(Real.Floor(R(x))), DoubleDouble.Floor(x));
            Assert.Equal(DoubleDouble.FromReal(Real.Ceiling(R(x))), DoubleDouble.Ceiling(x));
            Assert.Equal(DoubleDouble.FromReal(Real.Round(R(x))), DoubleDouble.Round(x));
        }
    }

    [Fact]
    public void ScaleBIsExact()
    {
        DoubleDouble x = DoubleDouble.Pi;
        Assert.Equal(new DoubleDouble(Math.ScaleB(x.Hi, 100), Math.ScaleB(x.Lo, 100)), DoubleDouble.ScaleB(x, 100));
        Assert.True(DoubleDouble.IsPositiveInfinity(DoubleDouble.ScaleB(x, 1100)));
        Assert.Equal(0.0, DoubleDouble.ScaleB(x, 1100).Lo);
        Assert.Equal(0.0, DoubleDouble.ScaleB(x, -1200).Hi);
    }

    // ---------------------------------------------------------------
    // Transcendental functions vs Real
    // ---------------------------------------------------------------

    private static DoubleDouble[] Sample(int seed, int count, int minExp, int maxExp, bool positive = false)
    {
        Random rng = new(seed);
        return Enumerable.Range(0, count).Select(_ => RandomDD(rng, minExp, maxExp, !positive)).ToArray();
    }

    private static DoubleDouble[] Uniform(int seed, int count, double lo, double hi)
    {
        Random rng = new(seed);
        return Enumerable.Range(0, count).Select(_ =>
        {
            double h = lo + (hi - lo) * rng.NextDouble();
            return new DoubleDouble(h, h * (rng.NextDouble() - 0.5) * Math.ScaleB(1.0, -52));
        }).ToArray();
    }

    private void Check(string name, IEnumerable<DoubleDouble> xs, Func<DoubleDouble, DoubleDouble> f, Func<Real, Real> exact, double log2Bound = -100) =>
        AssertBound(name, xs.Select(x => (f(x), exact(R(x)), S(x))), log2Bound);

    [Fact]
    public void ExponentialFunctions()
    {
        DoubleDouble[] wide = [.. Uniform(20, 1500, -670, 709), .. Sample(21, 500, -60, 2)];
        Check("exp", wide, DoubleDouble.Exp, x => Real.Exp(x, P));
        Check("exp2", [.. Uniform(22, 1000, -965, 1023), .. Sample(23, 300, -60, 0)], DoubleDouble.Exp2, x => Real.Exp2(x, P));
        Check("exp10", [.. Uniform(24, 1000, -290, 308), .. Sample(25, 300, -60, 0)], DoubleDouble.Exp10, x => Real.Exp10(x, P));
        Check("expm1", [.. Sample(26, 1500, -300, 3), .. Uniform(27, 500, -50, 50)], DoubleDouble.ExpM1, x => Real.ExpM1(x, P));
        Check("exp2m1", Sample(28, 800, -300, 3), DoubleDouble.Exp2M1, x => Real.Exp2M1(x, P));
        Check("exp10m1", Sample(29, 800, -300, 2), DoubleDouble.Exp10M1, x => Real.Exp10M1(x, P));

        Assert.Equal(DoubleDouble.One, DoubleDouble.Exp(0.0));
        Assert.Equal(DoubleDouble.One, DoubleDouble.Exp(-0.0));
        Assert.Equal(DoubleDouble.E, DoubleDouble.Exp(1.0));
        Assert.True(DoubleDouble.IsPositiveInfinity(DoubleDouble.Exp(710.0)));
        Assert.True(DoubleDouble.IsPositiveInfinity(DoubleDouble.Exp(double.PositiveInfinity)));
        Assert.Equal(DoubleDouble.Zero, DoubleDouble.Exp(-750.0));
        Assert.Equal(DoubleDouble.Zero, DoubleDouble.Exp(double.NegativeInfinity));
        Assert.True(DoubleDouble.IsNaN(DoubleDouble.Exp(double.NaN)));
        Assert.Equal((DoubleDouble)1024.0, DoubleDouble.Exp2(10.0));
        Assert.True(Log2RelErr(DoubleDouble.Exp10(22.0), Rd(1e22)) < -104);
        Assert.Equal(DoubleDouble.NegativeOne, DoubleDouble.ExpM1(double.NegativeInfinity));
        Assert.True(Log2RelErr(DoubleDouble.Exp(709.78), Real.Exp(Rd(709.78), P)) < -100);
    }

    [Fact]
    public void LogarithmicFunctions()
    {
        DoubleDouble[] pos = [.. Sample(30, 1500, -1020, 1020, positive: true), .. Uniform(31, 500, 0.5, 2.0)];
        Check("log", pos, DoubleDouble.Log, x => Real.Log(x, P));
        Check("log2", pos, DoubleDouble.Log2, x => Real.Log2(x, P));
        Check("log10", pos, DoubleDouble.Log10, x => Real.Log10(x, P));
        // Values just around 1 (the cancellation-prone region).
        DoubleDouble[] near1 = Sample(32, 1000, -100, -2).Select(d => d + 1.0).ToArray();
        Check("log-near-1", near1, DoubleDouble.Log, x => Real.Log(x, P));
        DoubleDouble[] p1 = [.. Sample(33, 1500, -300, 6), .. Uniform(34, 300, -0.999, 0)];
        Check("logp1", p1.Where(x => x > -1.0), DoubleDouble.LogP1, x => Real.LogP1(x, P));
        Check("log2p1", p1.Where(x => x > -1.0), DoubleDouble.Log2P1, x => Real.Log2P1(x, P));
        Check("log10p1", p1.Where(x => x > -1.0), DoubleDouble.Log10P1, x => Real.Log10P1(x, P));
        Check("log-subnormal", [double.Epsilon, 3 * double.Epsilon, Math.ScaleB(1.5, -1050)], DoubleDouble.Log, x => Real.Log(x, P));

        Assert.Equal(DoubleDouble.Zero, DoubleDouble.Log(1.0));
        Assert.Equal((DoubleDouble)(-3.0), DoubleDouble.Log2(0.125));
        Assert.Equal((DoubleDouble)1000.0, DoubleDouble.Log2(Math.ScaleB(1.0, 1000)));
        Assert.True(DoubleDouble.IsNegativeInfinity(DoubleDouble.Log(0.0)));
        Assert.True(DoubleDouble.IsNaN(DoubleDouble.Log(-1.0)));
        Assert.True(DoubleDouble.IsPositiveInfinity(DoubleDouble.Log(double.PositiveInfinity)));
        Assert.True(DoubleDouble.IsNegativeInfinity(DoubleDouble.LogP1(-1.0)));
        Assert.True(DoubleDouble.IsNaN(DoubleDouble.LogP1(-2.0)));
        Assert.True(Log2RelErr(DoubleDouble.Log(8.0, 2.0), Rd(3)) < -102);
    }

    [Fact]
    public void PowerFunctions()
    {
        Random rng = new(40);
        List<(DoubleDouble X, DoubleDouble Y)> pairs = [];
        for (int i = 0; i < 1500; i++)
        {
            DoubleDouble x = DoubleDouble.Abs(RandomDD(rng, -20, 20));
            // keep |y·log x| <= ~16 so that the exp(y·log x) error stays below 2^-100
            double ly = Math.Abs(Math.Log(x.Hi));
            DoubleDouble y = Uniform(rng.Next(), 1, -16 / Math.Max(ly, 1), 16 / Math.Max(ly, 1))[0];
            pairs.Add((x, y));
        }
        AssertBound("pow", pairs.Select(p => (DoubleDouble.Pow(p.X, p.Y), Real.Pow(R(p.X), R(p.Y), P), S(p.X) + S(p.Y))), -100);
        AssertBound("pow-int", Sample(41, 1000, -20, 20).SelectMany(x => new[] { 2, 3, 7, -5, 30 }.Select(n =>
            (DoubleDouble.Pow(x, n), Real.Pow(R(x), (long)n, P), S(x) + " n=" + n))), -100);

        Assert.Equal((DoubleDouble)243.0, DoubleDouble.Pow(3.0, 5.0));
        Assert.Equal((DoubleDouble)(-243.0), DoubleDouble.Pow(-3.0, 5.0));
        Assert.Equal(DoubleDouble.One, DoubleDouble.Pow(DoubleDouble.NaN, 0.0));
        Assert.Equal(DoubleDouble.One, DoubleDouble.Pow(1.0, DoubleDouble.NaN));
        Assert.True(DoubleDouble.IsNaN(DoubleDouble.Pow(-2.0, 0.5)));
        Assert.True(DoubleDouble.IsPositiveInfinity(DoubleDouble.Pow(0.0, -1.0)));
        Assert.True(DoubleDouble.IsNegativeInfinity(DoubleDouble.Pow(-0.0, -3.0)));
        Assert.Equal(DoubleDouble.Zero, DoubleDouble.Pow(2.0, double.NegativeInfinity));
        Assert.True(DoubleDouble.IsPositiveInfinity(DoubleDouble.Pow(new DoubleDouble(1.0, 1e-20), double.PositiveInfinity)));
        Assert.Equal(DoubleDouble.One, DoubleDouble.Pow(-1.0, double.PositiveInfinity));
        Assert.Equal((DoubleDouble)0.25, DoubleDouble.Pow(2.0, -2));
        Assert.Equal((DoubleDouble)2.0, DoubleDouble.RootN(8.0, 3));
        Assert.Equal((DoubleDouble)(-2.0), DoubleDouble.RootN(-8.0, 3));
        Assert.Equal((DoubleDouble)0.5, DoubleDouble.RootN(16.0, -4));
        Assert.True(DoubleDouble.IsNaN(DoubleDouble.RootN(-16.0, 4)));
        Assert.True(Log2RelErr(DoubleDouble.RootN(double.Epsilon, 5000), Real.RootN(Rd(double.Epsilon), 5000L, P)) < -100);
        Assert.True(Log2RelErr(DoubleDouble.RootN(double.MaxValue, 7), Real.RootN(Rd(double.MaxValue), 7L, P)) < -100);
    }

    [Fact]
    public void Hypot()
    {
        AssertBound("hypot", RandomPairs(42, 2000, -500, 500).Select(p => (DoubleDouble.Hypot(p.A, p.B), Real.Hypot(R(p.A), R(p.B), P), S(p.A) + S(p.B))), -103);
        Assert.Equal((DoubleDouble)5.0, DoubleDouble.Hypot(3.0, -4.0));
        Assert.Equal(DoubleDouble.Hypot(Math.ScaleB(3.0, 1000), Math.ScaleB(4.0, 1000)), (DoubleDouble)Math.ScaleB(5.0, 1000));
        Assert.True(DoubleDouble.IsPositiveInfinity(DoubleDouble.Hypot(DoubleDouble.NaN, double.NegativeInfinity)));
    }

    [Fact]
    public void TrigonometricFunctions()
    {
        DoubleDouble[] xs = [.. Uniform(50, 1500, -10, 10), .. Sample(51, 500, -300, 0), .. Uniform(52, 300, -1e6, 1e6)];
        Check("sin", xs, DoubleDouble.Sin, x => Real.Sin(x, P));
        Check("cos", xs, DoubleDouble.Cos, x => Real.Cos(x, P));
        Check("tan", xs, DoubleDouble.Tan, x => Real.Tan(x, P));
        Check("sincos.sin", xs, x => DoubleDouble.SinCos(x).Sin, x => Real.Sin(x, P));
        Check("sincos.cos", xs, x => DoubleDouble.SinCos(x).Cos, x => Real.Cos(x, P));
        // Huge arguments and near-multiples of π/2 go through the exact Real reduction.
        DoubleDouble[] hard =
        [
            .. Sample(53, 100, 20, 1020),
            new DoubleDouble(Math.PI / 2, 6.123233995736766E-17),
            DoubleDouble.Pi * 1000.0,
            DoubleDouble.Pi * 1e300,
            1e22, 6381956970095103.0 * Math.ScaleB(1.0, 797), // classic hard double for reduction
        ];
        Check("sin-hard", hard, DoubleDouble.Sin, x => Real.Sin(x, P));
        Check("cos-hard", hard, DoubleDouble.Cos, x => Real.Cos(x, P));

        DoubleDouble[] unit = [.. Uniform(54, 1500, -1, 1), .. Sample(55, 300, -300, -1)];
        Check("asin", unit, DoubleDouble.Asin, x => Real.Asin(x, P));
        Check("acos", unit, DoubleDouble.Acos, x => Real.Acos(x, P));
        DoubleDouble[] any = [.. Uniform(56, 1000, -5, 5), .. Sample(57, 800, -300, 300)];
        Check("atan", any, DoubleDouble.Atan, x => Real.Atan(x, P));
        AssertBound("atan2", RandomPairs(58, 2000, -300, 300).Select(p => (DoubleDouble.Atan2(p.A, p.B), Real.Atan2(R(p.A), R(p.B), P), S(p.A) + S(p.B))), -100);

        Assert.True(DoubleDouble.IsNaN(DoubleDouble.Sin(double.PositiveInfinity)));
        Assert.True(double.IsNegative(DoubleDouble.Sin(-0.0).Hi));
        Assert.Equal(DoubleDouble.One, DoubleDouble.Cos(0.0));
        Assert.Equal(DoubleDouble.Pi, DoubleDouble.Atan2(0.0, -1.0));
        Assert.Equal(-DoubleDouble.Pi, DoubleDouble.Atan2(-0.0, -1.0));
        Assert.True(double.IsNegative(DoubleDouble.Atan2(-0.0, 1.0).Hi));
        Assert.Equal(DoubleDouble.Pi / 2.0, DoubleDouble.Atan2(1.0, 0.0));
        Assert.Equal(DoubleDouble.Pi / 4.0, DoubleDouble.Atan(1.0));
        Assert.Equal(DoubleDouble.Pi / 2.0, DoubleDouble.Asin(1.0));
        Assert.Equal(DoubleDouble.Pi, DoubleDouble.Acos(-1.0));
        Assert.True(DoubleDouble.IsNaN(DoubleDouble.Asin(1.5)));
        Assert.Equal(DoubleDouble.Pi * 0.75, DoubleDouble.Atan2(1.0, -1.0));
    }

    [Fact]
    public void PiScaledTrigonometricFunctions()
    {
        DoubleDouble[] xs = [.. Uniform(60, 1500, -10, 10), .. Sample(61, 300, -300, 0), .. Uniform(62, 200, -1e15, 1e15)];
        Check("sinpi", xs, DoubleDouble.SinPi, x => Real.SinPi(x, P));
        Check("cospi", xs, DoubleDouble.CosPi, x => Real.CosPi(x, P));
        Check("tanpi", xs, DoubleDouble.TanPi, x => Real.TanPi(x, P));
        DoubleDouble[] unit = Uniform(63, 800, -1, 1);
        Check("asinpi", unit, DoubleDouble.AsinPi, x => Real.AsinPi(x, P));
        Check("acospi", unit, DoubleDouble.AcosPi, x => Real.AcosPi(x, P));
        Check("atanpi", Uniform(64, 800, -50, 50), DoubleDouble.AtanPi, x => Real.AtanPi(x, P));
        AssertBound("atan2pi", RandomPairs(65, 800, -30, 30).Select(p => (DoubleDouble.Atan2Pi(p.A, p.B), Real.Atan2Pi(R(p.A), R(p.B), P), S(p.A) + S(p.B))), -100);

        Assert.Equal(DoubleDouble.Zero, DoubleDouble.SinPi(3.0));
        Assert.Equal(DoubleDouble.NegativeOne, DoubleDouble.CosPi(3.0));
        Assert.Equal(DoubleDouble.Zero, DoubleDouble.CosPi(2.5));
        Assert.Equal(DoubleDouble.One, DoubleDouble.SinPi(0.5));
        Assert.True(DoubleDouble.IsPositiveInfinity(DoubleDouble.TanPi(0.5)));
        Assert.True(DoubleDouble.IsNegativeInfinity(DoubleDouble.TanPi(1.5)));
        Assert.Equal(DoubleDouble.One, DoubleDouble.TanPi(0.25));
        Assert.Equal(DoubleDouble.Zero, DoubleDouble.SinPi(Math.ScaleB(1.0, 80)));
        Assert.Equal((DoubleDouble)0.5, DoubleDouble.AsinPi(1.0));
        Assert.Equal((DoubleDouble)0.25, DoubleDouble.AtanPi(1.0));
        Assert.Equal((DoubleDouble)0.75, DoubleDouble.Atan2Pi(1.0, -1.0));
        Assert.Equal(DoubleDouble.Pi, DoubleDouble.DegreesToRadians(180.0));
        Assert.True(Log2RelErr(DoubleDouble.RadiansToDegrees(DoubleDouble.Pi), Rd(180)) < -102);
    }

    [Fact]
    public void HyperbolicFunctions()
    {
        DoubleDouble[] xs = [.. Uniform(70, 1500, -30, 30), .. Sample(71, 500, -300, 1), .. Uniform(72, 200, -709, 709)];
        Check("sinh", xs, DoubleDouble.Sinh, x => Real.Sinh(x, P));
        Check("cosh", xs, DoubleDouble.Cosh, x => Real.Cosh(x, P));
        Check("tanh", xs, DoubleDouble.Tanh, x => Real.Tanh(x, P));
        Check("asinh", [.. xs, .. Sample(73, 300, 0, 1000)], DoubleDouble.Asinh, x => Real.Asinh(x, P));
        DoubleDouble[] ge1 = [.. Sample(74, 800, -300, -1, positive: true).Select(d => d + 1.0), .. Sample(75, 800, 0, 1000, positive: true)];
        Check("acosh", ge1, DoubleDouble.Acosh, x => Real.Acosh(x, P));
        DoubleDouble[] unit = [.. Uniform(76, 1000, -0.999999, 0.999999), .. Sample(77, 500, -300, -1)];
        Check("atanh", unit, DoubleDouble.Atanh, x => Real.Atanh(x, P));

        Assert.True(DoubleDouble.IsPositiveInfinity(DoubleDouble.Cosh(-800.0)));
        Assert.True(DoubleDouble.IsNegativeInfinity(DoubleDouble.Sinh(-800.0)));
        Assert.Equal(DoubleDouble.One, DoubleDouble.Tanh(1000.0));
        Assert.True(DoubleDouble.IsPositiveInfinity(DoubleDouble.Atanh(1.0)));
        Assert.True(DoubleDouble.IsNaN(DoubleDouble.Atanh(1.5)));
        Assert.True(DoubleDouble.IsNaN(DoubleDouble.Acosh(0.5)));
        Assert.Equal(DoubleDouble.Zero, DoubleDouble.Acosh(1.0));
    }

    // ---------------------------------------------------------------
    // Formatting and parsing
    // ---------------------------------------------------------------

    [Fact]
    public void ToStringRoundTrips()
    {
        Random rng = new(80);
        DoubleDouble[] xs =
        [
            .. Enumerable.Range(0, 2000).Select(_ => RandomDD(rng, -1000, 1000)),
            new DoubleDouble(1.0, Math.ScaleB(1.0, -200)), // widely separated words
            DoubleDouble.MaxValue, DoubleDouble.MinValue, double.Epsilon, -0.0, 0.0,
            DoubleDouble.Pi, DoubleDouble.E,
        ];
        foreach (DoubleDouble x in xs)
        {
            string s = x.ToString();
            DoubleDouble back = DoubleDouble.Parse(s);
            Assert.True(back.Equals(x), $"{S(x)} -> {s} -> {S(back)}");
            Assert.Equal(double.IsNegative(x.Hi), double.IsNegative(back.Hi));
        }
        Assert.Equal("3.1415926535897932384626433832795", DoubleDouble.Pi.ToString());
        Assert.Equal("0.1", DoubleDouble.Parse("0.1").ToString());
        Assert.Equal("NaN", DoubleDouble.NaN.ToString());
        Assert.Equal("-Infinity", DoubleDouble.NegativeInfinity.ToString());
        Assert.Equal("3.1416", DoubleDouble.Pi.ToString("G5", null));
        Assert.Throws<FormatException>(() => DoubleDouble.Pi.ToString("X", null));
    }

    [Fact]
    public void ParseIsNearest()
    {
        DoubleDouble tenth = DoubleDouble.Parse("0.1");
        Real exact = Real.Divide(Rd(1), Rd(10), P);
        Assert.Equal(DoubleDouble.FromReal(exact), tenth);
        Assert.True(Log2RelErr(tenth, exact) < -106);
        Assert.Equal(DoubleDouble.One / 3.0, DoubleDouble.Parse("0.33333333333333333333333333333333333333333"));
        Assert.Equal(DoubleDouble.Pi, DoubleDouble.Parse("3.14159265358979323846264338327950288419716939937510"));
        Assert.True(DoubleDouble.IsPositiveInfinity(DoubleDouble.Parse("1e400")));
        Assert.Equal(DoubleDouble.Zero, DoubleDouble.Parse("1e-400"));
        Assert.True(DoubleDouble.IsNaN(DoubleDouble.Parse("NaN")));
        Assert.True(DoubleDouble.IsNegativeInfinity(DoubleDouble.Parse("-Infinity")));
        Assert.False(DoubleDouble.TryParse("abc", null, out _));
        Assert.Throws<FormatException>(() => DoubleDouble.Parse("1.2.3"));
    }

    [Fact]
    public void ParseWithStylesAndUtf8()
    {
        CultureInfo de = CultureInfo.GetCultureInfo("de-DE");
        Assert.Equal((DoubleDouble)1234.5, DoubleDouble.Parse("1.234,5", NumberStyles.Float | NumberStyles.AllowThousands, de));
        Assert.Equal((DoubleDouble)(-12.0), DoubleDouble.Parse("(12)", NumberStyles.AllowParentheses, CultureInfo.InvariantCulture));
        Assert.False(DoubleDouble.TryParse("1e5", NumberStyles.Integer, CultureInfo.InvariantCulture, out _));
        Assert.Throws<ArgumentException>(() => DoubleDouble.Parse("ff", NumberStyles.HexNumber, null));
        Assert.Equal(DoubleDouble.Pi, DoubleDouble.Parse(Encoding.UTF8.GetBytes(DoubleDouble.Pi.ToString()), null));
        Assert.Equal((DoubleDouble)2.5, DoubleDouble.Parse("2.5"u8, NumberStyles.Float, CultureInfo.InvariantCulture));

        Span<char> chars = stackalloc char[64];
        Assert.True(DoubleDouble.E.TryFormat(chars, out int n, default, null));
        Assert.Equal(DoubleDouble.E.ToString(), chars[..n].ToString());
        Assert.False(DoubleDouble.E.TryFormat(chars[..5], out _, default, null));
        Span<byte> bytes = stackalloc byte[64];
        Assert.True(DoubleDouble.E.TryFormat(bytes, out int nb, "G10", null));
        Assert.Equal("2.718281828", Encoding.UTF8.GetString(bytes[..nb]));
    }

    // ---------------------------------------------------------------
    // Conversions
    // ---------------------------------------------------------------

    [Fact]
    public void IntegerConversionsAreExact()
    {
        long[] longs = [long.MaxValue, long.MinValue, long.MaxValue - 1, (1L << 53) + 1, -(1L << 60) - 7, 0, -1];
        foreach (long v in longs)
        {
            DoubleDouble d = v;
            Assert.Equal(R(d), Real.FromBigInt(new BigInt(v), 128));
            Assert.Equal(v, (long)d);
            Assert.Equal(d.Hi, d.Hi + d.Lo);
        }
        ulong[] ulongs = [ulong.MaxValue, ulong.MaxValue - 1, (1UL << 63) + 1, 12345];
        foreach (ulong v in ulongs)
        {
            DoubleDouble d = v;
            Assert.Equal(v, (ulong)d);
            Assert.Equal(d.Hi, d.Hi + d.Lo);
        }
        Assert.Equal(long.MaxValue, (long)((DoubleDouble)long.MaxValue));
        // Mixed operators keep 64-bit integers exact.
        Assert.Equal((DoubleDouble)long.MaxValue, DoubleDouble.Zero + long.MaxValue);
        Assert.Equal(long.MaxValue - 1, (long)((DoubleDouble)long.MaxValue - 1L));
        Assert.Equal((DoubleDouble)ulong.MaxValue, DoubleDouble.One * ulong.MaxValue);
        Assert.Equal(UInt128.MaxValue >> 30, (UInt128)(DoubleDouble)(UInt128.MaxValue >> 30));
        Assert.Equal(Int128.MinValue, (Int128)(DoubleDouble)Int128.MinValue);
        Assert.Equal((Int128)1 << 100, (Int128)(DoubleDouble)((Int128)1 << 100));
        Assert.Equal(7, (int)(DoubleDouble)7.9);
        Assert.Equal(-7, (int)(DoubleDouble)(-7.9));
        Assert.Throws<OverflowException>(() => (int)(DoubleDouble)3e9);
        Assert.Throws<OverflowException>(() => (long)DoubleDouble.NaN);
        Assert.Throws<OverflowException>(() => (ulong)(DoubleDouble)(-1.0));
        Assert.Equal(new BigInt(1) << 1000, (BigInt)(DoubleDouble)Math.ScaleB(1.0, 1000));
        Assert.Equal(BigInt.Parse("12345678901234567890123"), (BigInt)DoubleDouble.Parse("12345678901234567890123.75"));
    }

    [Fact]
    public void FloatingAndDecimalConversions()
    {
        Assert.Equal(DoubleDouble.FromReal(Real.FromDecimal(0.1m, P)), (DoubleDouble)0.1m);
        Assert.Equal(79228162514264337593543950335m, (decimal)(DoubleDouble)79228162514264337593543950335m);
        Assert.Equal(1.2345678901234567890123456789m, (decimal)(DoubleDouble)1.2345678901234567890123456789m);
        Assert.Throws<OverflowException>(() => (decimal)(DoubleDouble)1e30);
        Assert.True(double.IsNegative(((DoubleDouble)new decimal(0, 0, 0, true, 1)).Hi));
        // float conversion rounds the exact value, not Hi (which is a float midpoint here).
        double mid = 1.0 + Math.ScaleB(1.0, -24);
        Assert.Equal(MathF.BitIncrement(1.0f), (float)new DoubleDouble(mid, 1e-30));
        Assert.Equal(1.0f, (float)new DoubleDouble(mid, -1e-30));
        Assert.Equal((System.Half)1.5, (System.Half)(DoubleDouble)1.5);
        Assert.Equal(DoubleDouble.Pi.Hi, (double)DoubleDouble.Pi);
        Assert.Equal(DoubleDouble.Pi, (DoubleDouble)Real.Pi(P));
        Assert.Equal(new DoubleDouble(1.0, Math.ScaleB(1.0, -200)), (DoubleDouble)(BigRational)new DoubleDouble(1.0, Math.ScaleB(1.0, -200)));
        Assert.Equal((DoubleDouble)0.5, (DoubleDouble)new BigRational(new BigInt(1), new BigInt(2)));
        Assert.True(DoubleDouble.IsPositiveInfinity((DoubleDouble)Real.FromBigInt(BigInt.One << 2000)));
        Assert.Equal(DoubleDouble.MaxValue, DoubleDouble.FromReal(R(DoubleDouble.MaxValue)));
        Real exactSum = R(new DoubleDouble(1.0, Math.ScaleB(1.0, -1000)));
        Assert.True(Real.Subtract(Real.Add(Rd(1), Rd(Math.ScaleB(1.0, -1000)), 2000), exactSum, 2000).IsZero);
    }

    [Fact]
    public void GenericMathConversions()
    {
        Assert.Equal((DoubleDouble)42.0, DoubleDouble.CreateChecked(42));
        Assert.Equal((DoubleDouble)long.MinValue, DoubleDouble.CreateChecked(long.MinValue));
        Assert.Equal((DoubleDouble)1.5, DoubleDouble.CreateChecked(1.5f));
        Assert.Equal((DoubleDouble)0.25m, DoubleDouble.CreateChecked(0.25m));
        Assert.Equal((DoubleDouble)12.0, DoubleDouble.CreateChecked(new BigInteger(12)));
        Assert.Equal(int.MaxValue, int.CreateSaturating((DoubleDouble)1e20));
        Assert.Equal(0, int.CreateSaturating(DoubleDouble.NaN));
        Assert.Equal(int.MinValue, int.CreateTruncating(DoubleDouble.NegativeInfinity));
        Assert.Throws<OverflowException>(() => int.CreateChecked((DoubleDouble)1e20));
        Assert.Equal(byte.MaxValue, byte.CreateSaturating((DoubleDouble)300.0));
        Assert.Equal(UInt128.MaxValue, UInt128.CreateSaturating((DoubleDouble)1e60));
        Assert.Equal(long.MaxValue, long.CreateChecked((DoubleDouble)long.MaxValue));
        Assert.Equal(DoubleDouble.Pi.Hi, double.CreateChecked(DoubleDouble.Pi));
        Assert.Equal(decimal.MaxValue, decimal.CreateSaturating((DoubleDouble)1e40));
        Assert.Equal(new BigInteger(3), BigInteger.CreateChecked((DoubleDouble)3.7));
        Assert.Equal(DoubleDouble.Pi, DoubleDouble.CreateChecked(Real.Pi(P)));
    }

    private static T Horner<T>(T x, params T[] coefficients) where T : IFloatingPoint<T>
    {
        T acc = T.Zero;
        foreach (T c in coefficients)
        {
            acc = acc * x + c;
        }
        return acc;
    }

    private static T Generic<T>(T x) where T : IFloatingPoint<T>, IRootFunctions<T>, IExponentialFunctions<T>, ILogarithmicFunctions<T>, ITrigonometricFunctions<T>
        => T.Sqrt(T.Exp(T.Log(x))) + T.Sin(T.Pi / T.CreateChecked(6)) + T.Round(x * T.One, 3) - T.Floor(x);

    [Fact]
    public void WorksInGenericCode()
    {
        DoubleDouble x = DoubleDouble.Parse("1.25");
        Assert.Equal(Horner(1.25, 1.0, 2.0, 3.0), (double)Horner(x, 1.0, 2.0, 3.0));
        DoubleDouble g = Generic(x);
        double gd = Generic(1.25);
        Assert.True(Math.Abs((double)g - gd) < 1e-14);
        Assert.Equal(DoubleDouble.Pi, DoubleDouble.CreateChecked(Math.PI) + DoubleDouble.Pi.Lo);
        IFloatingPoint<DoubleDouble> fp = DoubleDouble.Pi;
        Assert.Equal(1, fp.GetExponentShortestBitLength() >= 1 ? 1 : 0);
        Assert.InRange(fp.GetSignificandBitLength(), 53, 108);
        Span<byte> buf = stackalloc byte[32];
        Assert.True(fp.TryWriteExponentLittleEndian(buf, out int w) && w == 4 && buf[0] == 1);
    }
}
