# Sdcb.Arithmetic [![CI](https://github.com/sdcb/Arithmetic/actions/workflows/ci.yml/badge.svg)](https://github.com/sdcb/Arithmetic/actions/workflows/ci.yml) [![NuGet](https://img.shields.io/nuget/v/Sdcb.Arithmetic.svg)](https://www.nuget.org/packages/Sdcb.Arithmetic) [![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://github.com/sdcb/Arithmetic/blob/main/LICENSE) [![QQ](https://img.shields.io/badge/QQ_Group-495782587-52B6EF?style=social&logo=tencent-qq&logoColor=000&logoWidth=20)](https://qm.qq.com/cgi-bin/qm/qr?_wv=1027&k=&authKey=&noverify=0&group_code=495782587)

**English** | [中文](https://github.com/sdcb/Arithmetic/blob/main/README_CN.md)

> Looking for the GMP/MPFR P/Invoke bindings (`Sdcb.Arithmetic.Gmp` / `Sdcb.Arithmetic.Mpfr`) that used to
> live in this repository? They moved to [sdcb/Sdcb.Arithmetic](https://github.com/sdcb/Sdcb.Arithmetic).

High-performance, **pure managed C#** arbitrary-precision arithmetic with GMP/MPFR-class
functionality and no native dependencies. It runs anywhere .NET 10 runs (x64, ARM64, AOT),
uses AVX-512/AVX2/SSE/NEON where available, and falls back to portable code elsewhere.

| Type | Models | Description |
|---|---|---|
| `BigInt` | GMP `mpz_t` | Immutable arbitrary-precision signed integer |
| `BigRational` | GMP `mpq_t` | Immutable rational, always canonical |
| `Real` | MPFR `mpfr_t` | Correctly-rounded binary floating point: five rounding modes, NaN/±∞/±0, 100+ functions |

All three implement .NET generic math (`IBinaryInteger<T>`, `INumber<T>`,
`IBinaryFloatingPointIeee754<T>`) and interoperate with `System.Numerics.BigInteger`.

## Install

```
dotnet add package Sdcb.Arithmetic
```

## Quick start

```csharp
using Sdcb.Arithmetic;

// Integers (mpz)
BigInt a = BigInt.Parse("123456789012345678901234567890");
BigInt p = BigInt.Pow(2, 521) - 1;                  // Mersenne prime M521
bool prime = p.IsProbablePrime();                   // BPSW + Miller-Rabin
BigInt g = BigInt.Gcd(a, p);
BigInt m = BigInt.ModPow(a, p, BigInt.Pow(10, 30));
Console.WriteLine(p.ToString(16));                  // any base 2..62

// Rationals (mpq), always canonical
BigRational r = new BigRational(1, 3) + new BigRational(1, 6); // 1/2
double d = r.ToDouble();                            // correctly rounded

// Correctly-rounded floating point (mpfr)
Real.DefaultPrecision = 256;                        // bits
Real pi = Real.Pi();
Real y = Real.Sin(Real.Parse("1.5"));
Real z = Real.Exp(y, precision: 512, RealRounding.TowardPositive);
Console.WriteLine(z);                               // round-trip decimal output
```

## Performance

AMD Ryzen 7 5800X (AVX2), .NET 10, operand sizes in 64-bit limbs.
Full tables and methodology: [docs/performance.md](https://github.com/sdcb/Arithmetic/blob/main/docs/performance.md).

| Operation | Speedup over `BigInteger` | Time vs native GMP 6.3 (lower is better) |
|---|---|---|
| Multiply | 1.9–2.8× (≤ 16 limbs), 3.7–5× (32–2048), 7–28× (4096–65536) | 1.3–1.6× (≥ 16 limbs) |
| Divide 2N/N | 2.2–4.3× (16–4096 limbs), 5.9× (16384) | 1.5–2.6× |
| Modular power | 4.2–7× | 1.3–2.8× |
| GCD | 1.5–5.8× | 0.9–1.4× |
| Add / sub / shift | 2.5–11× | 1.1–3.1× |
| Decimal ToString / Parse | up to 16× / 4× (large); ToString 0.9× at ≤ 4 limbs | 0.5–1.9× / 0.6–1.4× |

`Real` against native MPFR 4.2 (53–16384 bits, both sides correctly rounded):
add/mul 1.3–2.5×, div/sqrt 1.3–3.3×, exp/sin/atan 0.9–1.8×, log/pow 1.6–2.7×.

Every `BigInt` operation is cross-checked against `BigInteger` and native GMP, and every
`Real` function bit-for-bit against MPFR in all rounding modes.

## Documentation

- [Examples](https://github.com/sdcb/Arithmetic/blob/main/examples/README.md): 1,000,000 digits of π, RSA, rounding modes, generic math (also shipped as LINQPad samples in the NuGet package)
- [Real: features and functions](https://github.com/sdcb/Arithmetic/blob/main/docs/real.md)
- [Generic math and conversion semantics](https://github.com/sdcb/Arithmetic/blob/main/docs/generic-math.md)
- [GMP / MPFR → `BigInt` / `Real` API mapping](https://github.com/sdcb/Arithmetic/blob/main/docs/gmp-mapping.md)
- [Performance](https://github.com/sdcb/Arithmetic/blob/main/docs/performance.md)
- [Design, algorithms and verification](https://github.com/sdcb/Arithmetic/blob/main/docs/design.md)
- [References](https://github.com/sdcb/Arithmetic/blob/main/docs/references.md)

## WeChat group

![](https://io.starworks.cc:88/cv-public/2026/ocr-wxg-qr.png?1012)

If the WeChat QR code has expired, join the QQ group [.NET骚操作 495782587](https://qm.qq.com/cgi-bin/qm/qr?_wv=1027&k=&authKey=&noverify=0&group_code=495782587).

## License

MIT. `Sdcb.Arithmetic` is an independent C# implementation of published algorithms and
contains no GMP or MPFR source code. Native GMP/MPFR are used only as black-box references
for testing and benchmarking. See [docs/references.md](https://github.com/sdcb/Arithmetic/blob/main/docs/references.md) for the literature.
