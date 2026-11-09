using BenchmarkDotNet.Attributes;

namespace Sdcb.Arithmetic.Benchmarks;

/// <summary>
/// QuadDouble (≈212 bits) against Real and MPFR at 212 bits, with double for reference.
/// Operands: a = √2, b = √3; the unary functions take a.
/// </summary>
[MemoryDiagnoser]
public class QuadDoubleBench
{
    private const int P = 212;

    private QuadDouble _qa, _qb;
    private Real _ra, _rb;
    private double _da, _db;
    private Mpfr.MpfrT _ma, _mb, _mr;

    [GlobalSetup]
    public void Setup()
    {
        _qa = QuadDouble.Sqrt(2);
        _qb = QuadDouble.Sqrt(3);
        _ra = Real.Sqrt(Real.FromDouble(2.0, P), P);
        _rb = Real.Sqrt(Real.FromDouble(3.0, P), P);
        _da = Math.Sqrt(2);
        _db = Math.Sqrt(3);
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

    [Benchmark] public QuadDouble QdAdd() => _qa + _qb;
    [Benchmark] public Real RealAdd() => Real.Add(_ra, _rb, P);
    [Benchmark] public void MpfrAdd() => Mpfr.Add(ref _mr, ref _ma, ref _mb, 0);
    [Benchmark] public double DoubleAdd() => _da + _db;

    // Opposite signs: QuadDouble takes the cancellation-safe merge path here.
    [Benchmark] public QuadDouble QdSub() => _qa - _qb;
    [Benchmark] public void MpfrSub() => Mpfr.Sub(ref _mr, ref _ma, ref _mb, 0);

    [Benchmark] public QuadDouble QdMul() => _qa * _qb;
    [Benchmark] public Real RealMul() => Real.Multiply(_ra, _rb, P);
    [Benchmark] public void MpfrMul() => Mpfr.Mul(ref _mr, ref _ma, ref _mb, 0);
    [Benchmark] public double DoubleMul() => _da * _db;

    [Benchmark] public QuadDouble QdDiv() => _qa / _qb;
    [Benchmark] public Real RealDiv() => Real.Divide(_ra, _rb, P);
    [Benchmark] public void MpfrDiv() => Mpfr.Div(ref _mr, ref _ma, ref _mb, 0);
    [Benchmark] public double DoubleDiv() => _da / _db;

    [Benchmark] public QuadDouble QdSqrt() => QuadDouble.Sqrt(_qa);
    [Benchmark] public Real RealSqrt() => Real.Sqrt(_ra, P);
    [Benchmark] public void MpfrSqrt() => Mpfr.Sqrt(ref _mr, ref _ma, 0);
    [Benchmark] public double DoubleSqrt() => Math.Sqrt(_da);

    [Benchmark] public QuadDouble QdExp() => QuadDouble.Exp(_qa);
    [Benchmark] public Real RealExp() => Real.Exp(_ra, P);
    [Benchmark] public void MpfrExp() => Mpfr.Exp(ref _mr, ref _ma, 0);
    [Benchmark] public double DoubleExp() => Math.Exp(_da);

    [Benchmark] public QuadDouble QdLog() => QuadDouble.Log(_qa);
    [Benchmark] public Real RealLog() => Real.Log(_ra, P);
    [Benchmark] public void MpfrLog() => Mpfr.Log(ref _mr, ref _ma, 0);
    [Benchmark] public double DoubleLog() => Math.Log(_da);

    [Benchmark] public QuadDouble QdSin() => QuadDouble.Sin(_qa);
    [Benchmark] public Real RealSin() => Real.Sin(_ra, P);
    [Benchmark] public void MpfrSin() => Mpfr.Sin(ref _mr, ref _ma, 0);
    [Benchmark] public double DoubleSin() => Math.Sin(_da);
}
