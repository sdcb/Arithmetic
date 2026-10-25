using System.Numerics;
using System.Runtime.InteropServices;
using BenchmarkDotNet.Attributes;
using BenchmarkDotNet.Configs;
using BenchmarkDotNet.Jobs;
using BenchmarkDotNet.Running;

namespace Sdcb.Arithmetic.Benchmarks;

public static class Program
{
    public static void Main(string[] args)
    {
        if (args.Length > 0 && args[0] == "smoke")
        {
            Smoke();
            return;
        }
        if (args.Length > 0 && args[0] == "smoke-mpfr")
        {
            SmokeMpfr();
            return;
        }
        if (args.Length > 0 && args[0] == "gmp-exact")
        {
            GmpExact.Run(args);
            return;
        }
        if (args.Length > 0 && args[0] == "exact")
        {
            MpfrExact.Run(args);
            return;
        }
        if (args.Length > 0 && args[0] == "perf")
        {
            Perf.Run(args);
            return;
        }
        if (args.Length > 0 && args[0] == "quick")
        {
            Quick();
            return;
        }
        BenchmarkSwitcher
            .FromTypes([typeof(MulBench), typeof(AddBench), typeof(DivBench), typeof(ToStringBench), typeof(PowModBench), typeof(GcdBench), typeof(RealArithBench), typeof(RealTransBench)])
            .Run(args, DefaultConfig.Instance.AddJob(Job.ShortRun.WithWarmupCount(2).WithIterationCount(5)));
    }

    /// <summary>Fast manual timing loop for inner-loop tuning (not a rigorous benchmark).</summary>
    private static void Quick()
    {
        var rng = new Random(1);
        foreach (int n in (int[])[16, 64, 256, 1024, 4096])
        {
            ulong[] la = RandomLimbs(rng, n), lb = RandomLimbs(rng, n);
            BigInt a = FromLimbs(la), b = FromLimbs(lb);
            var za = Gmp.FromLimbs(la);
            var zb = Gmp.FromLimbs(lb);
            Gmp.MpzT zr = default;
            Gmp.Init2(ref zr, (ulong)(2 * n + 2) * 64);

            int iters = Math.Max(8, 2_000_000 / (n * n / 8 + 1));
            // warmup
            for (int i = 0; i < iters / 4 + 1; i++) { _ = a * b; Gmp.Mul(ref zr, ref za, ref zb); }

            var sw = System.Diagnostics.Stopwatch.StartNew();
            for (int i = 0; i < iters; i++) { _ = a * b; }
            double ours = sw.Elapsed.TotalNanoseconds / iters;

            sw.Restart();
            for (int i = 0; i < iters; i++) { Gmp.Mul(ref zr, ref za, ref zb); }
            double gmp = sw.Elapsed.TotalNanoseconds / iters;

            Console.WriteLine($"mul {n,5} limbs: ours {ours,12:F1} ns | gmp {gmp,12:F1} ns | ratio {ours / gmp:F2}");
            Gmp.Clear(ref za); Gmp.Clear(ref zb); Gmp.Clear(ref zr);
        }
    }

    /// <summary>Cross-checks Real against native MPFR: bit-exact at 53 bits, value-exact at 200 bits.</summary>
    private static void SmokeMpfr()
    {
        var rng = new Random(7);
        int checks = 0;

        // 53-bit lane: both sides are correctly rounded, so doubles must match bit for bit.
        for (int i = 0; i < 500; i++)
        {
            double x = Math.ScaleB(rng.NextDouble() * 2 - 1, rng.Next(-40, 40));
            double y = Math.ScaleB(rng.NextDouble() * 2 - 1, rng.Next(-20, 20));
            Real rx = Real.FromDouble(x, 53), ry = Real.FromDouble(y, 53);
            var mx = Mpfr.FromDouble(x, 53);
            var my = Mpfr.FromDouble(y, 53);
            Mpfr.MpfrT mr = default;
            Mpfr.Init2(ref mr, 53);

            void Check(string what, Real ours)
            {
                double a = ours.ToDouble(), b = Mpfr.GetD(ref mr, 0);
                if (!a.Equals(b) && !(double.IsNaN(a) && double.IsNaN(b)))
                {
                    throw new InvalidOperationException($"{what}(x={x:R}, y={y:R}): ours={a:R} mpfr={b:R}");
                }
                checks++;
            }

            Mpfr.Add(ref mr, ref mx, ref my, 0); Check("add", Real.Add(rx, ry, 53));
            Mpfr.Mul(ref mr, ref mx, ref my, 0); Check("mul", Real.Multiply(rx, ry, 53));
            Mpfr.Div(ref mr, ref mx, ref my, 0); Check("div", Real.Divide(rx, ry, 53));
            Mpfr.Exp(ref mr, ref mx, 0); Check("exp", Real.Exp(rx, 53));
            Mpfr.Sin(ref mr, ref mx, 0); Check("sin", Real.Sin(rx, 53));
            Mpfr.Cos(ref mr, ref mx, 0); Check("cos", Real.Cos(rx, 53));
            Mpfr.Atan(ref mr, ref mx, 0); Check("atan", Real.Atan(rx, 53));
            if (x > 0)
            {
                Mpfr.Sqrt(ref mr, ref mx, 0); Check("sqrt", Real.Sqrt(rx, 53));
                Mpfr.Log(ref mr, ref mx, 0); Check("log", Real.Log(rx, 53));
                Mpfr.Pow(ref mr, ref mx, ref my, 0); Check("pow", Real.Pow(rx, ry, 53));
            }
            Mpfr.Clear(ref mx); Mpfr.Clear(ref my); Mpfr.Clear(ref mr);
        }

        // 200-bit lane: round-trip our decimal output into mpfr and require exact equality.
        const int P = 200;
        for (int i = 0; i < 40; i++)
        {
            double x = Math.ScaleB(rng.NextDouble() * 2 - 1, rng.Next(-8, 8));
            var mx = Mpfr.FromDouble(x, P);
            Mpfr.MpfrT mr = default, mo = default;
            Mpfr.Init2(ref mr, P);
            Mpfr.Init2(ref mo, P);

            void CheckHi(string what, Real ours)
            {
                Mpfr.SetStr(ref mo, ours.ToString(), 10, 0);
                double da = Mpfr.GetD(ref mo, 0), db = Mpfr.GetD(ref mr, 0);
                Mpfr.Sub(ref mo, ref mo, ref mr, 0);
                double diff = Math.Abs(Mpfr.GetD(ref mo, 0));
                double scale = Math.Max(Math.Abs(db), double.Epsilon);
                if (diff > scale * Math.Pow(2, -190))
                {
                    throw new InvalidOperationException($"{what}@200(x={x:R}): ours={da:R} mpfr={db:R} diff={diff}");
                }
                checks++;
            }

            Mpfr.Exp(ref mr, ref mx, 0); CheckHi("exp", Real.Exp(Real.FromDouble(x, P), P));
            Mpfr.Sin(ref mr, ref mx, 0); CheckHi("sin", Real.Sin(Real.FromDouble(x, P), P));
            Mpfr.Atan(ref mr, ref mx, 0); CheckHi("atan", Real.Atan(Real.FromDouble(x, P), P));
            if (x > 0)
            {
                Mpfr.Log(ref mr, ref mx, 0); CheckHi("log", Real.Log(Real.FromDouble(x, P), P));
            }
            Mpfr.Clear(ref mx); Mpfr.Clear(ref mr); Mpfr.Clear(ref mo);
        }
        {
            Mpfr.MpfrT mpi = default;
            Mpfr.Init2(ref mpi, P);
            Mpfr.ConstPi(ref mpi, 0);
            Mpfr.MpfrT mo = default;
            Mpfr.Init2(ref mo, P);
            Mpfr.SetStr(ref mo, Real.Pi(P).ToString(), 10, 0);
            Mpfr.Sub(ref mo, ref mo, ref mpi, 0);
            if (Math.Abs(Mpfr.GetD(ref mo, 0)) > Math.Pow(2, -190))
            {
                throw new InvalidOperationException("pi mismatch");
            }
            Mpfr.Clear(ref mpi); Mpfr.Clear(ref mo);
            checks++;
        }
        Console.WriteLine($"smoke-mpfr OK: {checks} checks, Real matches native MPFR");
    }

    /// <summary>Cross-checks Sdcb.Arithmetic against native GMP on random inputs.</summary>
    private static unsafe void Smoke()
    {
        var rng = new Random(1);
        for (int size = 1; size <= 512; size *= 2)
        {
            ulong[] la = RandomLimbs(rng, size);
            ulong[] lb = RandomLimbs(rng, size);
            BigInt a = FromLimbs(la), b = FromLimbs(lb);
            var za = Gmp.FromLimbs(la);
            var zb = Gmp.FromLimbs(lb);
            Gmp.MpzT zr = default;
            Gmp.Init2(ref zr, (ulong)(2 * size + 2) * 64);

            Gmp.Mul(ref zr, ref za, ref zb);
            Check(a * b, ref zr, $"mul {size}");
            Gmp.Add(ref zr, ref za, ref zb);
            Check(a + b, ref zr, $"add {size}");
            Gmp.Gcd(ref zr, ref za, ref zb);
            Check(BigInt.Gcd(a, b), ref zr, $"gcd {size}");

            Gmp.Clear(ref za);
            Gmp.Clear(ref zb);
            Gmp.Clear(ref zr);
        }
        Console.WriteLine("smoke OK: Sdcb.Arithmetic matches native GMP");
    }

    private static unsafe void Check(BigInt expected, ref Gmp.MpzT actual, string what)
    {
        IntPtr s = Gmp.GetStr(IntPtr.Zero, 16, ref actual);
        string gmpStr = Marshal.PtrToStringAnsi(s)!;
        Gmp.FreeStr(s, (nuint)(gmpStr.Length + 1));
        string ours = expected.ToString(16);
        if (gmpStr != ours)
        {
            throw new InvalidOperationException($"{what}: mismatch\ngmp:  {gmpStr}\nours: {ours}");
        }
        Console.WriteLine($"{what}: OK");
    }

    internal static ulong[] RandomLimbs(Random rng, int n)
    {
        ulong[] d = new ulong[n];
        foreach (ref ulong x in d.AsSpan())
        {
            x = (ulong)rng.NextInt64() | ((ulong)rng.NextInt64() << 62);
        }
        d[n - 1] |= 1UL << 63; // full length
        return d;
    }

    internal static BigInt FromLimbs(ulong[] limbs)
    {
        byte[] bytes = MemoryMarshal.AsBytes(limbs.AsSpan()).ToArray();
        return BigInt.FromMagnitudeBytes(bytes);
    }

    internal static BigInteger BigIntegerFromLimbs(ulong[] limbs)
    {
        byte[] bytes = MemoryMarshal.AsBytes(limbs.AsSpan()).ToArray();
        return new BigInteger(bytes, isUnsigned: true, isBigEndian: false);
    }
}

public class BenchBase
{
    protected static ulong[] Limbs(int n, int seed) => Program.RandomLimbs(new Random(seed), n);
}

[MemoryDiagnoser]
public class MulBench : BenchBase
{
    [Params(16, 64, 256, 1024, 4096)]
    public int N;

    private BigInt _a, _b;
    private BigInteger _na, _nb;
    private Gmp.MpzT _za, _zb, _zr;

    [GlobalSetup]
    public void Setup()
    {
        ulong[] la = Limbs(N, 100), lb = Limbs(N, 200);
        _a = Program.FromLimbs(la);
        _b = Program.FromLimbs(lb);
        _na = Program.BigIntegerFromLimbs(la);
        _nb = Program.BigIntegerFromLimbs(lb);
        _za = Gmp.FromLimbs(la);
        _zb = Gmp.FromLimbs(lb);
        Gmp.Init2(ref _zr, (ulong)(2 * N + 2) * 64);
    }

    [Benchmark(Baseline = true)]
    public void GmpNative() => Gmp.Mul(ref _zr, ref _za, ref _zb);

    [Benchmark]
    public BigInt SdcbArithmetic() => _a * _b;

    [Benchmark]
    public BigInteger SystemNumerics() => _na * _nb;
}

[MemoryDiagnoser]
public class AddBench : BenchBase
{
    [Params(256, 4096)]
    public int N;

    private BigInt _a, _b;
    private BigInteger _na, _nb;
    private Gmp.MpzT _za, _zb, _zr;

    [GlobalSetup]
    public void Setup()
    {
        ulong[] la = Limbs(N, 101), lb = Limbs(N, 201);
        _a = Program.FromLimbs(la);
        _b = Program.FromLimbs(lb);
        _na = Program.BigIntegerFromLimbs(la);
        _nb = Program.BigIntegerFromLimbs(lb);
        _za = Gmp.FromLimbs(la);
        _zb = Gmp.FromLimbs(lb);
        Gmp.Init2(ref _zr, (ulong)(N + 2) * 64);
    }

    [Benchmark(Baseline = true)]
    public void GmpNative() => Gmp.Add(ref _zr, ref _za, ref _zb);

    [Benchmark]
    public BigInt SdcbArithmetic() => _a + _b;

    [Benchmark]
    public BigInteger SystemNumerics() => _na + _nb;
}

[MemoryDiagnoser]
public class DivBench : BenchBase
{
    [Params(256, 1024, 4096, 16384)]
    public int N;

    private BigInt _a, _b;
    private BigInteger _na, _nb;
    private Gmp.MpzT _za, _zb, _zq, _zr;

    [GlobalSetup]
    public void Setup()
    {
        ulong[] la = Limbs(2 * N, 102), lb = Limbs(N, 202);
        _a = Program.FromLimbs(la);
        _b = Program.FromLimbs(lb);
        _na = Program.BigIntegerFromLimbs(la);
        _nb = Program.BigIntegerFromLimbs(lb);
        _za = Gmp.FromLimbs(la);
        _zb = Gmp.FromLimbs(lb);
        Gmp.Init2(ref _zq, (ulong)(N + 2) * 64);
        Gmp.Init2(ref _zr, (ulong)(N + 2) * 64);
    }

    [Benchmark(Baseline = true)]
    public void GmpNative() => Gmp.TdivQr(ref _zq, ref _zr, ref _za, ref _zb);

    [Benchmark]
    public BigInt SdcbArithmetic()
    {
        BigInt.DivRem(_a, _b, out BigInt q, out _);
        return q;
    }

    [Benchmark]
    public BigInteger SystemNumerics() => BigInteger.DivRem(_na, _nb, out _);
}

[MemoryDiagnoser]
public class ToStringBench : BenchBase
{
    [Params(64, 1024)]
    public int N;

    private BigInt _a;
    private BigInteger _na;
    private Gmp.MpzT _za;

    [GlobalSetup]
    public void Setup()
    {
        ulong[] la = Limbs(N, 103);
        _a = Program.FromLimbs(la);
        _na = Program.BigIntegerFromLimbs(la);
        _za = Gmp.FromLimbs(la);
    }

    [Benchmark(Baseline = true)]
    public unsafe string GmpNative()
    {
        IntPtr s = Gmp.GetStr(IntPtr.Zero, 10, ref _za);
        string result = Marshal.PtrToStringAnsi(s)!;
        Gmp.FreeStr(s, (nuint)(result.Length + 1));
        return result;
    }

    [Benchmark]
    public string SdcbArithmetic() => _a.ToString();

    [Benchmark]
    public string SystemNumerics() => _na.ToString();
}

[MemoryDiagnoser]
public class PowModBench : BenchBase
{
    [Params(16, 32)] // 1024-bit / 2048-bit modular exponentiation
    public int N;

    private BigInt _b, _e, _m;
    private BigInteger _nb, _ne, _nm;
    private Gmp.MpzT _zb, _ze, _zm, _zr;

    [GlobalSetup]
    public void Setup()
    {
        ulong[] lb = Limbs(N, 104), le = Limbs(N, 204), lm = Limbs(N, 304);
        lm[0] |= 1; // odd modulus
        _b = Program.FromLimbs(lb);
        _e = Program.FromLimbs(le);
        _m = Program.FromLimbs(lm);
        _nb = Program.BigIntegerFromLimbs(lb);
        _ne = Program.BigIntegerFromLimbs(le);
        _nm = Program.BigIntegerFromLimbs(lm);
        _zb = Gmp.FromLimbs(lb);
        _ze = Gmp.FromLimbs(le);
        _zm = Gmp.FromLimbs(lm);
        Gmp.Init2(ref _zr, (ulong)(N + 2) * 64);
    }

    [Benchmark(Baseline = true)]
    public void GmpNative() => Gmp.Powm(ref _zr, ref _zb, ref _ze, ref _zm);

    [Benchmark]
    public BigInt SdcbArithmetic() => BigInt.ModPow(_b, _e, _m);

    [Benchmark]
    public BigInteger SystemNumerics() => BigInteger.ModPow(_nb, _ne, _nm);
}

[MemoryDiagnoser]
public class GcdBench : BenchBase
{
    [Params(64, 1024)]
    public int N;

    private BigInt _a, _b;
    private BigInteger _na, _nb;
    private Gmp.MpzT _za, _zb, _zr;

    [GlobalSetup]
    public void Setup()
    {
        ulong[] la = Limbs(N, 105), lb = Limbs(N, 205);
        _a = Program.FromLimbs(la);
        _b = Program.FromLimbs(lb);
        _na = Program.BigIntegerFromLimbs(la);
        _nb = Program.BigIntegerFromLimbs(lb);
        _za = Gmp.FromLimbs(la);
        _zb = Gmp.FromLimbs(lb);
        Gmp.Init2(ref _zr, (ulong)(N + 2) * 64);
    }

    [Benchmark(Baseline = true)]
    public void GmpNative() => Gmp.Gcd(ref _zr, ref _za, ref _zb);

    [Benchmark]
    public BigInt SdcbArithmetic() => BigInt.Gcd(_a, _b);

    [Benchmark]
    public BigInteger SystemNumerics() => BigInteger.GreatestCommonDivisor(_na, _nb);
}
