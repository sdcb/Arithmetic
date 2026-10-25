using System.Diagnostics;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

namespace Sdcb.Arithmetic.Benchmarks;

/// <summary>
/// Fast side-by-side timing table (ours vs native GMP/MPFR). Not as rigorous as
/// BenchmarkDotNet but quick enough to drive tuning: `-- perf [filter]`.
/// </summary>
internal static unsafe class Perf
{
    private static string? s_filter;
    // PERF_SIZES=16,24,32 restricts the rows to those sizes (limbs / bits / k).
    private static readonly int[]? s_sizes = Environment.GetEnvironmentVariable("PERF_SIZES")?.Split(',').Select(int.Parse).ToArray();
    // Integer rows also time System.Numerics.BigInteger (speedup = BigInteger ns / ours ns); PERF_BCL=0 skips it.
    private static readonly bool s_noBcl = Environment.GetEnvironmentVariable("PERF_BCL") == "0";

    public static void Run(string[] args)
    {
        s_filter = args.Length > 1 ? args[1] : null;
        Console.WriteLine($"ISA: Vector512={System.Runtime.Intrinsics.Vector512.IsHardwareAccelerated} Vector256={System.Runtime.Intrinsics.Vector256.IsHardwareAccelerated} Vector128={System.Runtime.Intrinsics.Vector128.IsHardwareAccelerated}");
        Console.WriteLine($"{"op",-22}{"size",8}{"ours ns",14}{"native ns",14}{"ratio",8}{"BigInteger ns",14}{"speedup",9}");
        GmpOps();
        MpfrOps();
    }

    private static double Time(Action a)
    {
        // Warm up long enough for tiered compilation (call counting + 100 ms delay +
        // background rejit) to promote every method on the path to optimized code.
        // The delay restarts whenever new tier-0 code is jitted, so keep going until the
        // per-call time stops improving.
        var sw = Stopwatch.StartNew();
        var total = Stopwatch.StartNew();
        double prev = double.MaxValue;
        int stable = 0;
        while (total.ElapsedMilliseconds < 3000)
        {
            long k = 0;
            sw.Restart();
            while (sw.ElapsedMilliseconds < 40 || k < 16)
            {
                a();
                k++;
            }
            double per = sw.Elapsed.TotalNanoseconds / k;
            stable = per > prev * 0.97 && total.ElapsedMilliseconds > 200 ? stable + 1 : 0;
            prev = Math.Min(prev, per);
            if (stable >= 2)
            {
                break;
            }
        }
        long iters = 1;
        sw.Restart();
        while (true)
        {
            sw.Restart();
            for (long i = 0; i < iters; i++) a();
            if (sw.Elapsed.TotalMilliseconds >= 20) break;
            iters *= 2;
        }
        double best = double.MaxValue;
        for (int round = 0; round < 5; round++)
        {
            sw.Restart();
            for (long i = 0; i < iters; i++) a();
            best = Math.Min(best, sw.Elapsed.TotalNanoseconds / iters);
        }
        return best;
    }

    private static void Row(string op, int size, Action ours, Action native, Action? bcl = null)
    {
        if (s_filter is not null && !s_filter.Split(',').Any(f => f.StartsWith('=') ? op.Equals(f[1..], StringComparison.OrdinalIgnoreCase) : op.Contains(f, StringComparison.OrdinalIgnoreCase)))
        {
            return;
        }
        if (s_sizes is not null && !s_sizes.Contains(size))
        {
            return;
        }
        double o = Time(ours), n = Time(native);
        string bclCols = $"{"-",14}{"-",9}";
        if (bcl is not null && !s_noBcl)
        {
            double b = Time(bcl);
            bclCols = $"{b,14:F1}{b / o,9:F2}";
        }
        Console.WriteLine($"{op,-22}{size,8}{o,14:F1}{n,14:F1}{o / n,8:F2}{bclCols}");
    }

    private static void GmpOps()
    {
        var rng = new Random(5);
        foreach (int n in s_sizes ?? [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 16384, 65536])
        {
            ulong[] la = Program.RandomLimbs(rng, n), lb = Program.RandomLimbs(rng, n), l2 = Program.RandomLimbs(rng, 2 * n);
            BigInt a = Program.FromLimbs(la), b = Program.FromLimbs(lb), a2 = Program.FromLimbs(l2);
            Gmp.MpzT za = Gmp.FromLimbs(la), zb = Gmp.FromLimbs(lb), z2 = Gmp.FromLimbs(l2), zr = default, zq = default;
            Gmp.Init2(ref zr, (ulong)(4 * n + 4) * 64);
            Gmp.Init2(ref zq, (ulong)(4 * n + 4) * 64);
            BigInteger ba = (BigInteger)a, bb = (BigInteger)b, ba2 = (BigInteger)a2, bsink;
            BigInt sink;
            if (n >= 64)
            {
                Row("add", n, () => sink = a + b, () => Gmp.Add(ref zr, ref za, ref zb), () => bsink = ba + bb);
                Row("sub", n, () => sink = a - b, () => Gmp.Sub(ref zr, ref za, ref zb), () => bsink = ba - bb);
                Row("xor", n, () => sink = a ^ b, () => Gmp.Xor(ref zr, ref za, ref zb), () => bsink = ba ^ bb);
                Row("shl", n, () => sink = a << 13, () => Gmp.Mul2Exp(ref zr, ref za, 13), () => bsink = ba << 13);
                Row("shr", n, () => sink = a >> 13, () => Gmp.TdivQ2Exp(ref zr, ref za, 13), () => bsink = ba >> 13);
                Row("popcount", n, () => a.PopCount(), () => Gmp.Popcount(ref za), () => bsink = BigInteger.PopCount(ba));
            }
            Row("mul", n, () => sink = a * b, () => Gmp.Mul(ref zr, ref za, ref zb), () => bsink = ba * bb);
            Row("sqr", n, () => sink = a.Square(), () => Gmp.Mul(ref zr, ref za, ref za), () => bsink = ba * ba);
            if (n <= 16384)
            {
                Row("divrem 2n/n", n, () => BigInt.DivRem(a2, b, out sink, out _), () => Gmp.TdivQr(ref zq, ref zr, ref z2, ref zb), () => bsink = BigInteger.DivRem(ba2, bb, out _));
                BigInt prod = a * b;
                BigInteger bprod = (BigInteger)prod;
                Gmp.MpzT zp = default;
                Gmp.Init2(ref zp, (ulong)(2 * n + 2) * 64);
                Gmp.Mul(ref zp, ref za, ref zb);
                Row("divexact", n, () => sink = BigInt.DivExact(prod, b), () => Gmp.DivExact(ref zq, ref zp, ref zb), () => bsink = bprod / bb);
                Gmp.Clear(ref zp);
                Row("sqrt 2n", n, () => sink = BigInt.Sqrt(a2), () => Gmp.Sqrt(ref zr, ref z2));
            }
            if (n <= 4096)
            {
                Row("gcd", n, () => sink = BigInt.Gcd(a, b), () => Gmp.Gcd(ref zr, ref za, ref zb), () => bsink = BigInteger.GreatestCommonDivisor(ba, bb));
                Gmp.MpzT zs = default, zt = default;
                Gmp.Init2(ref zs, (ulong)(n + 2) * 64);
                Gmp.Init2(ref zt, (ulong)(n + 2) * 64);
                Row("gcdext", n, () => BigInt.ExtendedGcd(a, b, out _, out _), () => Gmp.GcdExt(ref zr, ref zs, ref zt, ref za, ref zb));
                Gmp.Clear(ref zs); Gmp.Clear(ref zt);
            }
            if (n <= 4096)
            {
                Row("tostring10", n, () => a.ToString(), () =>
                {
                    IntPtr s = Gmp.GetStr(IntPtr.Zero, 10, ref za);
                    Gmp.FreeStr(s, (nuint)(Marshal.PtrToStringAnsi(s)!.Length + 1));
                }, () => ba.ToString());
                string dec = a.ToString();
                Row("parse10", n, () => sink = BigInt.Parse(dec), () => Gmp.SetStr(ref zr, dec, 10), () => bsink = BigInteger.Parse(dec));
            }
            if (n is >= 4 and <= 64)
            {
                ulong[] lm = Program.RandomLimbs(rng, n);
                lm[0] |= 1;
                BigInt m = Program.FromLimbs(lm);
                Gmp.MpzT zm = Gmp.FromLimbs(lm);
                BigInteger bm = (BigInteger)m;
                Row("powmod", n, () => sink = BigInt.ModPow(a, b, m), () => Gmp.Powm(ref zr, ref za, ref zb, ref zm), () => bsink = BigInteger.ModPow(ba, bb, bm));
                Row("invert", n, () => BigInt.TryModInverse(a, m, out sink), () => Gmp.Invert(ref zr, ref za, ref zm));
                Row("jacobi", n, () => BigInt.Jacobi(a, m), () => Gmp.Jacobi(ref za, ref zm));
                Row("isprime", n, () => m.GetPrimality(25), () => Gmp.ProbabPrimeP(ref zm, 25));
                Gmp.Clear(ref zm);
            }
            Gmp.Clear(ref za); Gmp.Clear(ref zb); Gmp.Clear(ref z2); Gmp.Clear(ref zr); Gmp.Clear(ref zq);
        }
        {
            Gmp.MpzT zr = default;
            Gmp.Init2(ref zr, 64);
            BigInt sink;
            foreach (ulong k in (ulong[])[1000, 100000])
            {
                Row("factorial", (int)k, () => sink = BigInt.Factorial(k), () => Gmp.FacUi(ref zr, k));
                Row("fibonacci", (int)k, () => sink = BigInt.Fibonacci(k), () => Gmp.FibUi(ref zr, k));
                Row("binomial 2k,k", (int)k, () => sink = BigInt.Binomial(2 * k, k), () => Gmp.BinUiUi(ref zr, 2 * k, k));
            }
            Gmp.Clear(ref zr);
        }
    }

    private static void MpfrOps()
    {
        foreach (int p in (int[])[53, 128, 1024, 4096, 16384])
        {
            Real x = Real.Sqrt(Real.FromDouble(2.0, p), p), y = Real.Sqrt(Real.FromDouble(3.0, p), p);
            // Native structs live in unmanaged memory: lambdas capture only the pointers, so
            // the GC can never move what mpfr is writing to.
            Mpfr.MpfrT* mx = Alloc(p), my = Alloc(p), mr = Alloc(p);
            Mpfr.SetStr(ref *mx, x.ToString(), 10, 0);
            Mpfr.SetStr(ref *my, y.ToString(), 10, 0);
            Real sink;
            Row("r.add", p, () => sink = Real.Add(x, y, p), () => Mpfr.Add(ref *mr, ref *mx, ref *my, 0));
            Row("r.mul", p, () => sink = Real.Multiply(x, y, p), () => Mpfr.Mul(ref *mr, ref *mx, ref *my, 0));
            Row("r.div", p, () => sink = Real.Divide(x, y, p), () => Mpfr.Div(ref *mr, ref *mx, ref *my, 0));
            Row("r.sqrt", p, () => sink = Real.Sqrt(x, p), () => Mpfr.Sqrt(ref *mr, ref *mx, 0));
            bool big = p > 4096;
            if (!big || s_filter is not null)
            {
                Unary("r.exp", p, x, Real.Exp, "mpfr_exp", mx, mr);
                Unary("r.log", p, x, Real.Log, "mpfr_log", mx, mr);
                Unary("r.sin", p, x, Real.Sin, "mpfr_sin", mx, mr);
                Unary("r.atan", p, x, Real.Atan, "mpfr_atan", mx, mr);
                Row("r.pow", p, () => sink = Real.Pow(x, y, p), () => Mpfr.Pow(ref *mr, ref *mx, ref *my, 0));
                Row("r.pi", p, () => sink = Real.Pi(p + 1), () => { Mpfr.FreeCache(); Mpfr.ConstPi(ref *mr, 0); });
                Row("r.tostring", p, () => x.ToString(), () => Mpfr.GetDigits(ref *mx, 0));
            }
            Free(mx); Free(my); Free(mr);
        }
    }

    private static Mpfr.MpfrT* Alloc(int prec)
    {
        var p = (Mpfr.MpfrT*)NativeMemory.AllocZeroed((nuint)sizeof(Mpfr.MpfrT));
        Mpfr.Init2(ref *p, prec);
        return p;
    }

    private static void Free(Mpfr.MpfrT* p)
    {
        Mpfr.Clear(ref *p);
        NativeMemory.Free(p);
    }

    private static void Unary(string op, int p, Real x, Func<Real, int, RealRounding, Real> ours, string native, Mpfr.MpfrT* px, Mpfr.MpfrT* pr)
    {
        var f = Mpfr.Unary(native);
        Row(op, p, () => ours(x, p, RealRounding.NearestEven), () => f(ref *pr, ref *px, 0));
    }
}
