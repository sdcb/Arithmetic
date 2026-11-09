using System.Reflection;
using Xunit;
using Xunit.Abstractions;

namespace Sdcb.Arithmetic.Tests;

/// <summary>
/// The triple-double fast path of the Float128 elementary functions: bit-exact agreement with the Real oracle in
/// all five rounding modes on random and hard inputs, fallback rates, and kernel errors against their bounds.
/// </summary>
public class Float128FastPathTests(ITestOutputHelper output)
{
    private static readonly RealRounding[] Modes =
        [RealRounding.NearestEven, RealRounding.TowardZero, RealRounding.TowardPositive, RealRounding.TowardNegative, RealRounding.AwayFromZero];

    private static readonly UInt128 Sign = UInt128.One << 127;
    private static readonly UInt128 QNaN = new(0x7FFF_8000_0000_0000, 0);

    private static Float128 F(UInt128 bits) => Float128.FromBits(bits);
    private static string H(UInt128 b) => $"0x{(ulong)(b >> 64):X16}_{(ulong)b:X16}";

    // ---------------------------------------------------------------
    // Test hooks (internal, reached by reflection)
    // ---------------------------------------------------------------

    private delegate bool KernelHook(int function, Float128 x, Float128 y, out double h, out double m, out double l, out int scale, out double eps);
    private delegate bool TryHook(int function, Float128 x, Float128 y, RealRounding rounding, out Float128 result);

    private static readonly KernelHook Kernel = typeof(Float128).GetMethod("FastPathKernel", BindingFlags.NonPublic | BindingFlags.Static)!
        .CreateDelegate<KernelHook>();
    private static readonly TryHook TryFast = typeof(Float128).GetMethod("FastPathTry", BindingFlags.NonPublic | BindingFlags.Static)!
        .CreateDelegate<TryHook>();

    // Indices of Float128.FastFunction.
    private enum Fn { Exp, Exp2, Log, Log2, Log10, Pow, Sin, Cos, Tan, Asin, Acos, Atan, Atan2, Sinh, Cosh, Tanh }

    private sealed record Case(Fn Fn, Func<Float128, Float128, RealRounding, Float128> Api, Func<Real, Real, int, RealRounding, Real> Oracle,
        Func<Random, (Float128 X, Float128 Y)> Gen);

    // ---------------------------------------------------------------
    // Oracle: Real at 160 bits toward zero and away from zero; equal results are exact, otherwise the
    // toward-zero value with its last bit forced to 1 is the round-to-odd value, whose rounding to binary128
    // is the correctly rounded result in every mode (as in Float128Tests).
    // ---------------------------------------------------------------

    private const int OracleBits = 160;

    private static UInt128[] OracleAllModes(Func<int, RealRounding, Real> f)
    {
        var r = new UInt128[Modes.Length];
        Real z = f(OracleBits, RealRounding.TowardZero);
        if (z.IsZero || !z.IsNumber)
        {
            for (int i = 0; i < Modes.Length; i++)
            {
                r[i] = (z.IsZero ? f(OracleBits, Modes[i]) : z).ToBinary128Bits(Modes[i]);
            }
            return r;
        }
        Real a = f(OracleBits, RealRounding.AwayFromZero);
        Real v = z;
        if (a.CompareTo(z) != 0)
        {
            BigInt m = z.Mantissa | BigInt.One;
            v = Real.FromBigInt(z.SignBit ? -m : m, z.Exponent, OracleBits);
        }
        for (int i = 0; i < Modes.Length; i++)
        {
            r[i] = v.ToBinary128Bits(Modes[i]);
        }
        return r;
    }

    private static UInt128 Canon(UInt128 b) => Float128.IsNaN(F(b)) ? QNaN : b;

    // ---------------------------------------------------------------
    // Inputs
    // ---------------------------------------------------------------

    private static UInt128 RandomFrac(Random rng)
    {
        UInt128 f = new UInt128((ulong)rng.NextInt64(), (ulong)rng.NextInt64() ^ ((ulong)rng.Next() << 63)) & ((UInt128.One << 112) - 1);
        switch (rng.Next(10))
        {
            case 0: // sparse
                f = UInt128.Zero;
                for (int i = rng.Next(4); i >= 0; i--)
                {
                    f |= UInt128.One << rng.Next(112);
                }
                break;
            case 1: // runs of ones
                f = ((UInt128.One << 112) - 1) >> rng.Next(112);
                break;
        }
        return f;
    }

    /// <summary>±m·2^e with e uniform in [emin, emax] and a random 113-bit significand.</summary>
    private static Float128 RandomExp(Random rng, int emin, int emax, bool allowNegative = true)
    {
        int be = 16383 + rng.Next(emin, emax + 1);
        UInt128 b = ((UInt128)(uint)be << 112) | RandomFrac(rng);
        if (allowNegative && rng.Next(2) == 0)
        {
            b |= Sign;
        }
        return F(b);
    }

    /// <summary>Uniform in [lo, hi) with a full random significand.</summary>
    private static Float128 Uniform(Random rng, double lo, double hi)
    {
        Float128 u = F(((UInt128)0x3FFF << 112) | RandomFrac(rng)) - 1;
        return (Float128)lo + u * (Float128)(hi - lo);
    }

    /// <summary>A few ulps around v.</summary>
    private static Float128 Near(Random rng, Float128 v)
    {
        UInt128 b = v.ToBits();
        int d = rng.Next(-1000, 1001);
        return F(d >= 0 ? b + (UInt128)d : b - (UInt128)(-d));
    }

    private static Float128 Pick(Random rng, params Func<Float128>[] gens) => gens[rng.Next(gens.Length)]();

    private static readonly Float128 PiHalf = Float128.Pi / 2;

    private static Float128 NearMultipleOfPiHalf(Random rng, int maxK)
    {
        int k = rng.Next(-maxK, maxK + 1);
        Float128 v = F(Real.Multiply(Real.Pi(400), Real.FromDouble(k / 2.0), 400).ToBinary128Bits());
        return k == 0 ? RandomExp(rng, -200, -10) : Near(rng, v);
    }

    private static Real R(Float128 x) => x.ToReal();

    private static readonly Case[] Cases =
    [
        new(Fn.Exp, (x, _, m) => Float128.Exp(x, m), (x, _, p, m) => Real.Exp(x, p, m),
            rng => (Pick(rng, () => Uniform(rng, -60, 60), () => Uniform(rng, -11350, 11350), () => RandomExp(rng, -200, 13)), 0)),
        new(Fn.Exp2, (x, _, m) => Float128.Exp2(x, m), (x, _, p, m) => Real.Exp2(x, p, m),
            rng => (Pick(rng, () => Uniform(rng, -100, 100), () => Uniform(rng, -16380, 16380), () => RandomExp(rng, -200, 13)), 0)),
        new(Fn.Log, (x, _, m) => Float128.Log(x, m), (x, _, p, m) => Real.Log(x, p, m),
            rng => (Pick(rng, () => RandomExp(rng, -16382, 16383, false), () => Uniform(rng, 0.5, 2), () => 1 + RandomExp(rng, -112, -8)), 0)),
        new(Fn.Log2, (x, _, m) => Float128.Log2(x, m), (x, _, p, m) => Real.Log2(x, p, m),
            rng => (Pick(rng, () => RandomExp(rng, -16382, 16383, false), () => Uniform(rng, 0.5, 2), () => 1 + RandomExp(rng, -112, -8)), 0)),
        new(Fn.Log10, (x, _, m) => Float128.Log10(x, m), (x, _, p, m) => Real.Log10(x, p, m),
            rng => (Pick(rng, () => RandomExp(rng, -16382, 16383, false), () => Uniform(rng, 0.5, 2), () => 1 + RandomExp(rng, -112, -8)), 0)),
        new(Fn.Pow, Float128.Pow, Real.Pow,
            rng =>
            {
                Float128 x = Pick(rng, () => Uniform(rng, 0, 100), () => RandomExp(rng, -100, 100, false), () => 1 + RandomExp(rng, -60, -2));
                return (x, Pick(rng, () => Uniform(rng, -40, 40), () => RandomExp(rng, -30, 10), () => (Float128)rng.Next(-60, 61)));
            }),
        new(Fn.Sin, (x, _, m) => Float128.Sin(x, m), (x, _, p, m) => Real.Sin(x, p, m),
            rng => (Pick(rng, () => Uniform(rng, -10, 10), () => RandomExp(rng, -200, 29), () => NearMultipleOfPiHalf(rng, 1000)), 0)),
        new(Fn.Cos, (x, _, m) => Float128.Cos(x, m), (x, _, p, m) => Real.Cos(x, p, m),
            rng => (Pick(rng, () => Uniform(rng, -10, 10), () => RandomExp(rng, -200, 29), () => NearMultipleOfPiHalf(rng, 1000)), 0)),
        new(Fn.Tan, (x, _, m) => Float128.Tan(x, m), (x, _, p, m) => Real.Tan(x, p, m),
            rng => (Pick(rng, () => Uniform(rng, -10, 10), () => RandomExp(rng, -200, 29), () => NearMultipleOfPiHalf(rng, 1000)), 0)),
        new(Fn.Asin, (x, _, m) => Float128.Asin(x, m), (x, _, p, m) => Real.Asin(x, p, m),
            rng => (Pick(rng, () => Uniform(rng, -1, 1), () => RandomExp(rng, -200, -1), () => (rng.Next(2) == 0 ? 1 : -1) * (1 - RandomExp(rng, -113, -2, false))), 0)),
        new(Fn.Acos, (x, _, m) => Float128.Acos(x, m), (x, _, p, m) => Real.Acos(x, p, m),
            rng => (Pick(rng, () => Uniform(rng, -1, 1), () => RandomExp(rng, -200, -1), () => (rng.Next(2) == 0 ? 1 : -1) * (1 - RandomExp(rng, -113, -2, false))), 0)),
        new(Fn.Atan, (x, _, m) => Float128.Atan(x, m), (x, _, p, m) => Real.Atan(x, p, m),
            rng => (Pick(rng, () => Uniform(rng, -10, 10), () => RandomExp(rng, -200, 200), () => Uniform(rng, -1.01, 1.01)), 0)),
        new(Fn.Atan2, Float128.Atan2, Real.Atan2,
            rng => (Pick(rng, () => Uniform(rng, -100, 100), () => RandomExp(rng, -300, 300)), Pick(rng, () => Uniform(rng, -100, 100), () => RandomExp(rng, -300, 300)))),
        new(Fn.Sinh, (x, _, m) => Float128.Sinh(x, m), (x, _, p, m) => Real.Sinh(x, p, m),
            rng => (Pick(rng, () => Uniform(rng, -20, 20), () => Uniform(rng, -1, 1), () => RandomExp(rng, -200, 13)), 0)),
        new(Fn.Cosh, (x, _, m) => Float128.Cosh(x, m), (x, _, p, m) => Real.Cosh(x, p, m),
            rng => (Pick(rng, () => Uniform(rng, -20, 20), () => Uniform(rng, -1, 1), () => RandomExp(rng, -200, 13)), 0)),
        new(Fn.Tanh, (x, _, m) => Float128.Tanh(x, m), (x, _, p, m) => Real.Tanh(x, p, m),
            rng => (Pick(rng, () => Uniform(rng, -20, 20), () => Uniform(rng, -1, 1), () => RandomExp(rng, -200, 6)), 0)),
    ];

    public static TheoryData<int> AllFunctions()
    {
        var d = new TheoryData<int>();
        foreach (Case c in Cases)
        {
            d.Add((int)c.Fn);
        }
        return d;
    }

    // ---------------------------------------------------------------
    // Random inputs: every mode against the oracle
    // ---------------------------------------------------------------

    private const int RandomCount = 20000;

    [Theory]
    [MemberData(nameof(AllFunctions))]
    public void RandomInputsMatchOracleInAllModes(int fn)
    {
        Case c = Cases[fn];
        Random rng = new(1000 + fn);
        int attempts = 0, fallbacks = 0;
        for (int i = 0; i < RandomCount; i++)
        {
            (Float128 x, Float128 y) = c.Gen(rng);
            Real rx = R(x), ry = R(y);
            UInt128[] expected = OracleAllModes((p, m) => c.Oracle(rx, ry, p, m));
            for (int k = 0; k < Modes.Length; k++)
            {
                RealRounding mode = Modes[k];
                string ctx = $"{c.Fn}({H(x.ToBits())}, {H(y.ToBits())}) {mode}";
                if (TryFast(fn, x, y, mode, out Float128 fast))
                {
                    Assert.True(Canon(expected[k]) == fast.ToBits(), $"fast path {ctx}: expected {H(expected[k])}, got {H(fast.ToBits())}");
                }
                else if (Kernel(fn, x, y, out _, out _, out _, out _, out _))
                {
                    fallbacks++;
                }
                if (Kernel(fn, x, y, out _, out _, out _, out _, out _))
                {
                    attempts++;
                }
                Float128 api = c.Api(x, y, mode);
                Assert.True(Canon(expected[k]) == api.ToBits(), $"{ctx}: expected {H(expected[k])}, got {H(api.ToBits())}");
            }
        }
        // These generators deliberately include hard inputs (tiny arguments handled by the callers' shortcuts, near
        // multiples of π/2, exact powers), so the rate is only reported; FallbackRateOnTypicalInputs bounds it.
        output.WriteLine($"{c.Fn}: {attempts} fast-path attempts, {fallbacks} fallbacks ({100.0 * fallbacks / Math.Max(attempts, 1):F4}%)");
        Assert.True(attempts > RandomCount * Modes.Length / 2, $"{c.Fn}: only {attempts} inputs reached the fast path");
    }

    // ---------------------------------------------------------------
    // Kernel error against Real at 300 bits, relative to the bound used by the rounding test
    // ---------------------------------------------------------------

    private const int ErrorSamples = 4000;

    [Theory]
    [MemberData(nameof(AllFunctions))]
    public void KernelErrorsStayFarBelowTheirBounds(int fn)
    {
        Case c = Cases[fn];
        Random rng = new(2000 + fn);
        double worst = 0, worstRatio = 0, minEps = double.MaxValue;
        for (int i = 0; i < ErrorSamples; i++)
        {
            (Float128 x, Float128 y) = c.Gen(rng);
            if (!Kernel(fn, x, y, out double h, out double m, out double l, out int scale, out double eps))
            {
                continue;
            }
            Real exact = c.Oracle(R(x), R(y), 300, RealRounding.NearestEven);
            Real approx = Real.Add(Real.Add(Real.FromDouble(h), Real.FromDouble(m), 400), Real.FromDouble(l), 400).ShiftLeft(scale);
            double rel = Math.Abs(Real.Divide(Real.Subtract(approx, exact, 400), exact, 64).ToDouble());
            worst = Math.Max(worst, rel);
            worstRatio = Math.Max(worstRatio, rel / eps);
            minEps = Math.Min(minEps, eps);
            Assert.True(rel * 256 <= eps, $"{c.Fn}({H(x.ToBits())}, {H(y.ToBits())}): error 2^{Math.Log2(rel):F1} vs bound 2^{Math.Log2(eps):F1}");
        }
        output.WriteLine($"{c.Fn}: max error 2^{Math.Log2(worst):F1}, min bound 2^{Math.Log2(minEps):F1}, max error/bound 2^{Math.Log2(worstRatio):F1}");
    }

    // ---------------------------------------------------------------
    // Fallback rate on typical inputs (uniform over the usual domain)
    // ---------------------------------------------------------------

    /// <summary>Uniform in [lo, hi) with every significand bit random (no sparse patterns, which make results exact).</summary>
    private static Float128 Dense(Random rng, double lo, double hi)
    {
        UInt128 frac = new UInt128((ulong)rng.NextInt64() ^ ((ulong)rng.Next() << 40), (ulong)rng.NextInt64() ^ ((ulong)rng.Next() << 63)) & ((UInt128.One << 112) - 1);
        Float128 u = F(((UInt128)0x3FFF << 112) | frac) - 1;
        return (Float128)lo + u * (Float128)(hi - lo);
    }

    private static (Float128, Float128) Typical(Fn fn, Random rng) => fn switch
    {
        Fn.Exp => (Dense(rng, -700, 700), 0),
        Fn.Exp2 => (Dense(rng, -1000, 1000), 0),
        Fn.Log or Fn.Log2 or Fn.Log10 => (rng.Next(2) == 0 ? Dense(rng, 0, 1000) : Float128.ScaleB(Dense(rng, 1, 2), rng.Next(-1000, 1000)), 0),
        Fn.Pow => (Dense(rng, 0, 100), Dense(rng, -40, 40)),
        Fn.Sin or Fn.Cos or Fn.Tan or Fn.Atan => (Dense(rng, -100, 100), 0),
        Fn.Asin or Fn.Acos => (Dense(rng, -1, 1), 0),
        Fn.Atan2 => (Dense(rng, -100, 100), Dense(rng, -100, 100)),
        Fn.Sinh or Fn.Cosh => (Dense(rng, -50, 50), 0),
        _ => (Dense(rng, -20, 20), 0),
    };

    [Theory]
    [MemberData(nameof(AllFunctions))]
    public void FallbackRateOnTypicalInputs(int fn)
    {
        Random rng = new(3000 + fn);
        const int n = 200_000;
        int fallbacks = 0;
        for (int i = 0; i < n; i++)
        {
            (Float128 x, Float128 y) = Typical((Fn)fn, rng);
            if (!TryFast(fn, x, y, Modes[i % Modes.Length], out _))
            {
                fallbacks++;
            }
        }
        output.WriteLine($"{(Fn)fn}: {fallbacks} fallbacks in {n} ({100.0 * fallbacks / n:F4}%)");
        Assert.True(fallbacks * 10_000 < n * 10, $"{(Fn)fn}: fallback rate {fallbacks}/{n}"); // < 0.1%
    }

    // ---------------------------------------------------------------
    // Hard cases
    // ---------------------------------------------------------------
    private static void CheckAllModes(Case c, Float128 x, Float128 y)
    {
        Real rx = R(x), ry = R(y);
        UInt128[] expected = OracleAllModes((p, m) => c.Oracle(rx, ry, p, m));
        for (int k = 0; k < Modes.Length; k++)
        {
            Float128 api = c.Api(x, y, Modes[k]);
            Assert.True(Canon(expected[k]) == api.ToBits(), $"{c.Fn}({x}, {y}) {Modes[k]}: expected {H(expected[k])}, got {H(api.ToBits())}");
            if (TryFast((int)c.Fn, x, y, Modes[k], out Float128 fast))
            {
                Assert.True(Canon(expected[k]) == fast.ToBits(), $"fast {c.Fn}({x}, {y}) {Modes[k]}");
            }
        }
    }

    [Fact]
    public void ExactResultsFallBackInDirectedModes()
    {
        (Fn, Float128, Float128)[] exact =
        [
            (Fn.Exp2, 10, 0), (Fn.Exp2, -300, 0), (Fn.Exp2, 16383, 0), (Fn.Exp2, -16382, 0),
            (Fn.Log2, 1024, 0), (Fn.Log2, Float128.ScaleB(1, -16000), 0), (Fn.Log2, Float128.ScaleB(1, 16383), 0),
            (Fn.Log10, 1000, 0), (Fn.Log10, (Float128)1e30, 0), (Fn.Log10, 1, 0), (Fn.Log, 1, 0), (Fn.Log2, 1, 0),
            (Fn.Pow, 2, 10), (Fn.Pow, 4, 0.5), (Fn.Pow, 9, 1.5), (Fn.Pow, -2, 3), (Fn.Pow, 0.25, -0.5), (Fn.Pow, 10, 20),
            (Fn.Pow, (Float128)1.5, 2), (Fn.Pow, 7, 1), (Fn.Pow, 2, -1074), (Fn.Pow, 8, (Float128)1 / 3),
        ];
        foreach ((Fn fn, Float128 x, Float128 y) in exact)
        {
            Case c = Cases[(int)fn];
            CheckAllModes(c, x, y);
            for (int k = 1; k < Modes.Length; k++)
            {
                // Exactly representable results lie on a boundary in the directed modes: the fast path must decline.
                Real rx = R(x), ry = R(y);
                bool isExact = c.Oracle(rx, ry, OracleBits, RealRounding.TowardZero).CompareTo(c.Oracle(rx, ry, OracleBits, RealRounding.AwayFromZero)) == 0;
                if (isExact)
                {
                    Assert.False(TryFast((int)fn, x, y, Modes[k], out _), $"{fn}({x}, {y}) {Modes[k]} took the fast path");
                }
            }
        }
    }

    [Fact]
    public void RangeBoundaries()
    {
        Random rng = new(77);
        Float128[] exp = [11355, (Float128)11355.5, 11356, (Float128)11356.5, -11355, (Float128)(-11355.1), (Float128)(-11355.2), -11356, 11400, (Float128)709.78, (Float128)(-745.13)];
        Float128[] exp2 = [16383, (Float128)16383.5, (Float128)16383.999, -16382, (Float128)(-16381.5), (Float128)(-16382.5), 16384];
        Float128 tiny = Float128.ScaleB(1, -200);
        Float128[] small = [tiny, -tiny, Float128.BitDecrement(tiny), Float128.BitIncrement(tiny), Float128.ScaleB(1, -199), Float128.ScaleB((Float128)1.5, -150)];
        Float128[] trigBig = [Float128.ScaleB(1, 30), Float128.BitDecrement(Float128.ScaleB(1, 30)), Float128.ScaleB(Float128.Pi, 28), (Float128)1e9, (Float128)(-1e9)];
        Float128[] hyp = [(Float128)0.5, Float128.BitDecrement((Float128)0.5), 50, Float128.BitDecrement((Float128)50), 41, Float128.BitDecrement((Float128)41), (Float128)0.25, Float128.BitIncrement((Float128)0.25), 11356, 11357, (Float128)11358.5];
        foreach (Float128 v in exp)
        {
            CheckAllModes(Cases[(int)Fn.Exp], v, 0);
            CheckAllModes(Cases[(int)Fn.Exp], Near(rng, v), 0);
        }
        foreach (Float128 v in exp2)
        {
            CheckAllModes(Cases[(int)Fn.Exp2], v, 0);
            CheckAllModes(Cases[(int)Fn.Exp2], Near(rng, v), 0);
        }
        foreach (Float128 v in small)
        {
            foreach (Fn fn in new[] { Fn.Exp, Fn.Exp2, Fn.Sin, Fn.Cos, Fn.Tan, Fn.Asin, Fn.Acos, Fn.Atan, Fn.Sinh, Fn.Cosh, Fn.Tanh })
            {
                CheckAllModes(Cases[(int)fn], v, 0);
            }
            CheckAllModes(Cases[(int)Fn.Atan2], v, 1);
            CheckAllModes(Cases[(int)Fn.Atan2], 1, v);
            CheckAllModes(Cases[(int)Fn.Log], 1 + v * Float128.ScaleB(1, 120), 0);
        }
        foreach (Float128 v in trigBig)
        {
            foreach (Fn fn in new[] { Fn.Sin, Fn.Cos, Fn.Tan })
            {
                CheckAllModes(Cases[(int)fn], v, 0);
                CheckAllModes(Cases[(int)fn], -v, 0);
            }
        }
        foreach (Float128 v in hyp)
        {
            foreach (Fn fn in new[] { Fn.Sinh, Fn.Cosh, Fn.Tanh })
            {
                CheckAllModes(Cases[(int)fn], v, 0);
                CheckAllModes(Cases[(int)fn], -v, 0);
            }
        }
        // Atan2 exponent difference around the 600 limit, asin/acos near ±1, log at the ends of the range.
        foreach (int d in new[] { 590, 600, 601, 610 })
        {
            CheckAllModes(Cases[(int)Fn.Atan2], Float128.ScaleB((Float128)1.25, -d), 3);
            CheckAllModes(Cases[(int)Fn.Atan2], 3, Float128.ScaleB(-(Float128)1.25, -d));
            CheckAllModes(Cases[(int)Fn.Atan2], -Float128.ScaleB((Float128)1.25, -d), -3);
        }
        foreach (Float128 v in new[] { Float128.BitDecrement(1), Float128.BitDecrement(Float128.BitDecrement(1)), 1 - Float128.ScaleB(1, -60), (Float128)0.5, Float128.BitDecrement((Float128)0.5) })
        {
            foreach (Fn fn in new[] { Fn.Asin, Fn.Acos })
            {
                CheckAllModes(Cases[(int)fn], v, 0);
                CheckAllModes(Cases[(int)fn], -v, 0);
            }
        }
        foreach (Float128 v in new[] { Float128.MaxValue, Float128.MinNormal, Float128.BitIncrement(1), Float128.BitDecrement(1), Float128.Sqrt(2), Float128.BitIncrement(Float128.Sqrt(2)) })
        {
            foreach (Fn fn in new[] { Fn.Log, Fn.Log2, Fn.Log10 })
            {
                CheckAllModes(Cases[(int)fn], v, 0);
            }
        }
    }

    [Fact]
    public void TrigNearMultiplesOfPiOverTwo()
    {
        // The closest binary128 neighbours of kπ/2 (and kπ/256 table points) for k up to large values.
        Random rng = new(78);
        Real pi = Real.Pi(600);
        foreach (long k in new long[] { 1, 2, 3, 4, 5, 7, 8, 100, 255, 256, 257, 1000, 12345, 1_000_000, 123_456_789, 600_000_000 })
        {
            foreach (int denom in new[] { 2, 256 })
            {
                Real v = Real.Divide(Real.Multiply(pi, Real.FromBigInt(new BigInt(k)), 600), Real.FromDouble(denom), 600);
                Float128 c = Float128.FromBits(v.ToBinary128Bits());
                foreach (Float128 x in new[] { c, Float128.BitIncrement(c), Float128.BitDecrement(c), Near(rng, c) })
                {
                    foreach (Fn fn in new[] { Fn.Sin, Fn.Cos, Fn.Tan })
                    {
                        CheckAllModes(Cases[(int)fn], x, 0);
                        CheckAllModes(Cases[(int)fn], -x, 0);
                    }
                }
            }
        }
    }

    [Fact]
    public void SpecialValuesKeepTheirResults()
    {
        Float128[] specials = [0, Float128.NegativeZero, Float128.PositiveInfinity, Float128.NegativeInfinity, Float128.NaN, Float128.Epsilon, -Float128.Epsilon, 1, -1];
        foreach (Case c in Cases)
        {
            foreach (Float128 x in specials)
            {
                foreach (Float128 y in c.Fn is Fn.Pow or Fn.Atan2 ? specials : [0])
                {
                    CheckAllModes(c, x, y);
                }
            }
        }
    }

    [Fact]
    public void InvalidRoundingModeStillThrows()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => Float128.Exp(1, (RealRounding)17));
        Assert.Throws<ArgumentOutOfRangeException>(() => Float128.Sin(1, (RealRounding)17));
        Assert.Throws<ArgumentOutOfRangeException>(() => Float128.Pow(3, (Float128)0.7, (RealRounding)17));
    }
}
