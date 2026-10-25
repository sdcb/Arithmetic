using BenchmarkDotNet.Attributes;

namespace Sdcb.Arithmetic.Benchmarks;

[MemoryDiagnoser]
public class RealArithBench
{
    [Params(128, 1024, 8192)]
    public int P;

    private Real _a, _b;
    private Mpfr.MpfrT _ma, _mb, _mr;

    [GlobalSetup]
    public void Setup()
    {
        // Full-precision irrational-ish operands via sqrt.
        _a = Real.Sqrt(Real.FromDouble(2.0, P), P);
        _b = Real.Sqrt(Real.FromDouble(3.0, P), P);
        _ma = default; Mpfr.Init2(ref _ma, P);
        _mb = default; Mpfr.Init2(ref _mb, P);
        _mr = default; Mpfr.Init2(ref _mr, P);
        Mpfr.SetStr(ref _ma, _a.ToString(), 10, 0);
        Mpfr.SetStr(ref _mb, _b.ToString(), 10, 0);
    }

    [Benchmark(Baseline = true)]
    public void MpfrMul() => Mpfr.Mul(ref _mr, ref _ma, ref _mb, 0);

    [Benchmark]
    public Real OursMul() => Real.Multiply(_a, _b, P);

    [Benchmark]
    public void MpfrDiv() => Mpfr.Div(ref _mr, ref _ma, ref _mb, 0);

    [Benchmark]
    public Real OursDiv() => Real.Divide(_a, _b, P);

    [Benchmark]
    public void MpfrSqrt() => Mpfr.Sqrt(ref _mr, ref _ma, 0);

    [Benchmark]
    public Real OursSqrt() => Real.Sqrt(_a, P);

    [Benchmark]
    public void MpfrAdd() => Mpfr.Add(ref _mr, ref _ma, ref _mb, 0);

    [Benchmark]
    public Real OursAdd() => Real.Add(_a, _b, P);
}

[MemoryDiagnoser]
public class RealTransBench
{
    [Params(128, 1024, 4096)]
    public int P;

    private Real _x;
    private Mpfr.MpfrT _mx, _mr;

    [GlobalSetup]
    public void Setup()
    {
        _x = Real.Sqrt(Real.FromDouble(2.0, P), P); // ~1.414
        _mx = default; Mpfr.Init2(ref _mx, P);
        _mr = default; Mpfr.Init2(ref _mr, P);
        Mpfr.SetStr(ref _mx, _x.ToString(), 10, 0);
    }

    [Benchmark(Baseline = true)]
    public void MpfrExp() => Mpfr.Exp(ref _mr, ref _mx, 0);

    [Benchmark]
    public Real OursExp() => Real.Exp(_x, P);

    [Benchmark]
    public void MpfrLog() => Mpfr.Log(ref _mr, ref _mx, 0);

    [Benchmark]
    public Real OursLog() => Real.Log(_x, P);

    [Benchmark]
    public void MpfrSin() => Mpfr.Sin(ref _mr, ref _mx, 0);

    [Benchmark]
    public Real OursSin() => Real.Sin(_x, P);

    [Benchmark]
    public void MpfrAtan() => Mpfr.Atan(ref _mr, ref _mx, 0);

    [Benchmark]
    public Real OursAtan() => Real.Atan(_x, P);

    [Benchmark]
    public void MpfrPow() => Mpfr.Pow(ref _mr, ref _mx, ref _mx, 0);

    [Benchmark]
    public Real OursPow() => Real.Pow(_x, _x, P);
}
