using BenchmarkDotNet.Attributes;
using BenchmarkDotNet.Configs;

namespace Sdcb.Arithmetic.Benchmarks;

/// <summary>
/// DoubleDouble (~106 bits) against double, Real at 106 bits and native MPFR at 106 bits.
/// Each category's baseline is the DoubleDouble operation.
/// </summary>
[MemoryDiagnoser]
[GroupBenchmarksBy(BenchmarkLogicalGroupRule.ByCategory)]
[CategoriesColumn]
public class DoubleDoubleBench
{
    private const int P = 106;

    private DoubleDouble _a, _b;
    private double _da, _db;
    private Real _ra, _rb;
    private Mpfr.MpfrT _ma, _mb, _mr;

    [GlobalSetup]
    public void Setup()
    {
        _a = DoubleDouble.Sqrt(2.0);  // ~1.414, full 106-bit operands
        _b = DoubleDouble.Sqrt(3.0);
        _da = _a.Hi;
        _db = _b.Hi;
        _ra = Real.Sqrt(Real.FromDouble(2.0, P), P);
        _rb = Real.Sqrt(Real.FromDouble(3.0, P), P);
        _ma = default; Mpfr.Init2(ref _ma, P);
        _mb = default; Mpfr.Init2(ref _mb, P);
        _mr = default; Mpfr.Init2(ref _mr, P);
        Mpfr.SetStr(ref _ma, _ra.ToString(), 10, 0);
        Mpfr.SetStr(ref _mb, _rb.ToString(), 10, 0);
    }

    [GlobalCleanup]
    public void Cleanup()
    {
        Mpfr.Clear(ref _ma);
        Mpfr.Clear(ref _mb);
        Mpfr.Clear(ref _mr);
    }

    [Benchmark(Baseline = true), BenchmarkCategory("Add")] public DoubleDouble DD_Add() => _a + _b;
    [Benchmark, BenchmarkCategory("Add")] public double Double_Add() => _da + _db;
    [Benchmark, BenchmarkCategory("Add")] public Real Real_Add() => Real.Add(_ra, _rb, P);
    [Benchmark, BenchmarkCategory("Add")] public int Mpfr_Add() => Mpfr.Add(ref _mr, ref _ma, ref _mb, 0);

    [Benchmark(Baseline = true), BenchmarkCategory("Mul")] public DoubleDouble DD_Mul() => _a * _b;
    [Benchmark, BenchmarkCategory("Mul")] public double Double_Mul() => _da * _db;
    [Benchmark, BenchmarkCategory("Mul")] public Real Real_Mul() => Real.Multiply(_ra, _rb, P);
    [Benchmark, BenchmarkCategory("Mul")] public int Mpfr_Mul() => Mpfr.Mul(ref _mr, ref _ma, ref _mb, 0);

    [Benchmark(Baseline = true), BenchmarkCategory("Div")] public DoubleDouble DD_Div() => _a / _b;
    [Benchmark, BenchmarkCategory("Div")] public double Double_Div() => _da / _db;
    [Benchmark, BenchmarkCategory("Div")] public Real Real_Div() => Real.Divide(_ra, _rb, P);
    [Benchmark, BenchmarkCategory("Div")] public int Mpfr_Div() => Mpfr.Div(ref _mr, ref _ma, ref _mb, 0);

    [Benchmark(Baseline = true), BenchmarkCategory("Sqrt")] public DoubleDouble DD_Sqrt() => DoubleDouble.Sqrt(_a);
    [Benchmark, BenchmarkCategory("Sqrt")] public double Double_Sqrt() => Math.Sqrt(_da);
    [Benchmark, BenchmarkCategory("Sqrt")] public Real Real_Sqrt() => Real.Sqrt(_ra, P);
    [Benchmark, BenchmarkCategory("Sqrt")] public int Mpfr_Sqrt() => Mpfr.Sqrt(ref _mr, ref _ma, 0);

    [Benchmark(Baseline = true), BenchmarkCategory("Exp")] public DoubleDouble DD_Exp() => DoubleDouble.Exp(_a);
    [Benchmark, BenchmarkCategory("Exp")] public double Double_Exp() => Math.Exp(_da);
    [Benchmark, BenchmarkCategory("Exp")] public Real Real_Exp() => Real.Exp(_ra, P);
    [Benchmark, BenchmarkCategory("Exp")] public int Mpfr_Exp() => Mpfr.Exp(ref _mr, ref _ma, 0);

    [Benchmark(Baseline = true), BenchmarkCategory("Log")] public DoubleDouble DD_Log() => DoubleDouble.Log(_a);
    [Benchmark, BenchmarkCategory("Log")] public double Double_Log() => Math.Log(_da);
    [Benchmark, BenchmarkCategory("Log")] public Real Real_Log() => Real.Log(_ra, P);
    [Benchmark, BenchmarkCategory("Log")] public int Mpfr_Log() => Mpfr.Log(ref _mr, ref _ma, 0);

    [Benchmark(Baseline = true), BenchmarkCategory("Sin")] public DoubleDouble DD_Sin() => DoubleDouble.Sin(_a);
    [Benchmark, BenchmarkCategory("Sin")] public double Double_Sin() => Math.Sin(_da);
    [Benchmark, BenchmarkCategory("Sin")] public Real Real_Sin() => Real.Sin(_ra, P);
    [Benchmark, BenchmarkCategory("Sin")] public int Mpfr_Sin() => Mpfr.Sin(ref _mr, ref _ma, 0);
}
