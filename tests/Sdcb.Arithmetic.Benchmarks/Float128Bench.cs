using System.Globalization;
using BenchmarkDotNet.Attributes;

namespace Sdcb.Arithmetic.Benchmarks;

/// <summary>Float128 (IEEE binary128) against double, Real at 113 bits and MPFR at 113 bits.</summary>
[MemoryDiagnoser]
public class Float128Bench
{
    private const int P = 113;

    private Float128 _a, _b, _c;
    private double _da, _db, _dc;
    private Real _ra, _rb, _rc;
    private Mpfr.MpfrT _ma, _mb, _mc, _mr;
    private string _text = "";

    [GlobalSetup]
    public void Setup()
    {
        _a = Float128.Sqrt(2);
        _b = Float128.Sqrt(3);
        _c = Float128.One / 7;
        _da = Math.Sqrt(2);
        _db = Math.Sqrt(3);
        _dc = 1.0 / 7;
        _ra = _a;
        _rb = _b;
        _rc = _c;
        _ma = default; Mpfr.Init2(ref _ma, P);
        _mb = default; Mpfr.Init2(ref _mb, P);
        _mc = default; Mpfr.Init2(ref _mc, P);
        _mr = default; Mpfr.Init2(ref _mr, P);
        Mpfr.SetStr(ref _ma, _a.ToString(CultureInfo.InvariantCulture), 10, 0);
        Mpfr.SetStr(ref _mb, _b.ToString(CultureInfo.InvariantCulture), 10, 0);
        Mpfr.SetStr(ref _mc, _c.ToString(CultureInfo.InvariantCulture), 10, 0);
        _text = _a.ToString(CultureInfo.InvariantCulture);
    }

    // Add
    [Benchmark] public double DoubleAdd() => _da + _db;
    [Benchmark] public Float128 F128Add() => _a + _b;
    [Benchmark] public Real RealAdd() => Real.Add(_ra, _rb, P);
    [Benchmark] public int MpfrAdd() => Mpfr.Add(ref _mr, ref _ma, ref _mb, 0);

    // Sub (same exponent: exercises the cancellation/renormalization path)
    [Benchmark] public Float128 F128Sub() => _a - _b;
    [Benchmark] public int MpfrSub() => Mpfr.Sub(ref _mr, ref _ma, ref _mb, 0);

    // Mul
    [Benchmark] public double DoubleMul() => _da * _db;
    [Benchmark] public Float128 F128Mul() => _a * _b;
    [Benchmark] public Real RealMul() => Real.Multiply(_ra, _rb, P);
    [Benchmark] public int MpfrMul() => Mpfr.Mul(ref _mr, ref _ma, ref _mb, 0);

    // Div
    [Benchmark] public double DoubleDiv() => _da / _db;
    [Benchmark] public Float128 F128Div() => _a / _b;
    [Benchmark] public Real RealDiv() => Real.Divide(_ra, _rb, P);
    [Benchmark] public int MpfrDiv() => Mpfr.Div(ref _mr, ref _ma, ref _mb, 0);

    // Sqrt
    [Benchmark] public double DoubleSqrt() => Math.Sqrt(_db);
    [Benchmark] public Float128 F128Sqrt() => Float128.Sqrt(_b);
    [Benchmark] public Real RealSqrt() => Real.Sqrt(_rb, P);
    [Benchmark] public int MpfrSqrt() => Mpfr.Sqrt(ref _mr, ref _mb, 0);

    // FMA
    [Benchmark] public double DoubleFma() => Math.FusedMultiplyAdd(_da, _db, _dc);
    [Benchmark] public Float128 F128Fma() => Float128.FusedMultiplyAdd(_a, _b, _c);
    [Benchmark] public Real RealFma() => Real.FusedMultiplyAdd(_ra, _rb, _rc, P);

    // Formatting / parsing
    [Benchmark] public string DoubleToString() => _da.ToString(CultureInfo.InvariantCulture);
    [Benchmark] public string F128ToString() => _a.ToString(CultureInfo.InvariantCulture);
    [Benchmark] public string RealToString() => _ra.ToString();
    [Benchmark] public double DoubleParse() => double.Parse("1.4142135623730951", CultureInfo.InvariantCulture);
    [Benchmark] public Float128 F128Parse() => Float128.Parse(_text, CultureInfo.InvariantCulture);
    [Benchmark] public Real RealParse() => Real.Parse(_text, P);

    // Transcendentals (Float128 goes through Real in this version)
    [Benchmark] public double DoubleExp() => Math.Exp(_da);
    [Benchmark] public Float128 F128Exp() => Float128.Exp(_a);
    [Benchmark] public Real RealExp() => Real.Exp(_ra, P);
    [Benchmark] public int MpfrExp() => Mpfr.Exp(ref _mr, ref _ma, 0);

    [Benchmark] public double DoubleSin() => Math.Sin(_da);
    [Benchmark] public Float128 F128Sin() => Float128.Sin(_a);
    [Benchmark] public Real RealSin() => Real.Sin(_ra, P);
    [Benchmark] public int MpfrSin() => Mpfr.Sin(ref _mr, ref _ma, 0);
}
