# Sdcb.Arithmetic [![CI](https://github.com/sdcb/Arithmetic/actions/workflows/ci.yml/badge.svg)](https://github.com/sdcb/Arithmetic/actions/workflows/ci.yml) [![NuGet](https://img.shields.io/nuget/v/Sdcb.Arithmetic.svg)](https://www.nuget.org/packages/Sdcb.Arithmetic) [![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE) [![QQ](https://img.shields.io/badge/QQ_Group-495782587-52B6EF?style=social&logo=tencent-qq&logoColor=000&logoWidth=20)](https://qm.qq.com/cgi-bin/qm/qr?_wv=1027&k=&authKey=&noverify=0&group_code=495782587)

[English](README.md) | **中文**

> 以前放在本仓库的 GMP/MPFR P/Invoke 绑定（`Sdcb.Arithmetic.Gmp` / `Sdcb.Arithmetic.Mpfr`）已迁移至
> [sdcb/Sdcb.Arithmetic](https://github.com/sdcb/Sdcb.Arithmetic)。

高性能、**纯托管 C#** 的任意精度算术库，功能对标 GMP/MPFR，不依赖任何原生库。凡是 .NET 10
能运行的地方（x64、ARM64、AOT）都能用；有 AVX-512/AVX2/SSE/NEON 时自动使用，其余情况回退到可移植实现。

| 类型 | 对应 | 说明 |
|---|---|---|
| `BigInt` | GMP `mpz_t` | 不可变的任意精度有符号整数 |
| `BigRational` | GMP `mpq_t` | 不可变有理数，始终保持最简形式 |
| `Real` | MPFR `mpfr_t` | 正确舍入的二进制浮点数：五种舍入模式、NaN/±∞/±0、100 多个函数 |

三种类型都实现了 .NET 泛型数学接口（`IBinaryInteger<T>`、`INumber<T>`、
`IBinaryFloatingPointIeee754<T>`），并可与 `System.Numerics.BigInteger` 互转。

## 安装

```
dotnet add package Sdcb.Arithmetic
```

## 快速上手

```csharp
using Sdcb.Arithmetic;

// 整数（mpz）
BigInt a = BigInt.Parse("123456789012345678901234567890");
BigInt p = BigInt.Pow(2, 521) - 1;                  // 梅森素数 M521
bool prime = p.IsProbablePrime();                   // BPSW + Miller-Rabin
BigInt g = BigInt.Gcd(a, p);
BigInt m = BigInt.ModPow(a, p, BigInt.Pow(10, 30));
Console.WriteLine(p.ToString(16));                  // 支持 2..62 进制

// 有理数（mpq），始终为最简形式
BigRational r = new BigRational(1, 3) + new BigRational(1, 6); // 1/2
double d = r.ToDouble();                            // 正确舍入

// 正确舍入的浮点数（mpfr）
Real.DefaultPrecision = 256;                        // 单位：比特
Real pi = Real.Pi();
Real y = Real.Sin(Real.Parse("1.5"));
Real z = Real.Exp(y, precision: 512, RealRounding.TowardPositive);
Console.WriteLine(z);                               // 可往返的十进制输出
```

## 性能

测试环境：AMD Ryzen 7 5800X（AVX2），.NET 10，操作数大小以 64 位 limb 计。
完整数据和测试方法见 [docs/performance.md](docs/performance.md)（英文）。

| 操作 | 相比 `BigInteger` 的加速比 | 相对原生 GMP 6.3 的耗时（越低越好） |
|---|---|---|
| 乘法 | 1.9–2.8×（≤ 16 limbs），3.7–5×（32–2048），7–28×（4096–65536） | 1.3–1.6×（≥ 16 limbs） |
| 除法 2N/N | 2.2–4.3×（16–4096 limbs），5.9×（16384） | 1.5–2.6× |
| 模幂 | 4.2–7× | 1.3–2.8× |
| GCD | 1.5–5.8× | 0.9–1.4× |
| 加 / 减 / 移位 | 2.5–11× | 1.1–3.1× |
| 十进制 ToString / Parse | 大数最高 16× / 4×；≤ 4 limbs 时 ToString 为 0.9× | 0.5–1.9× / 0.6–1.4× |

`Real` 相对原生 MPFR 4.2（53–16384 位，两边都是正确舍入）：
add/mul 1.3–2.5×，div/sqrt 1.3–3.3×，exp/sin/atan 0.9–1.8×，log/pow 1.6–2.7×。

每个 `BigInt` 操作都与 `BigInteger` 和原生 GMP 交叉验证；每个 `Real` 函数都在全部舍入模式下与
MPFR 逐位比对。

## 文档（英文）

- [示例](examples/README.md)：计算 100 万位 π、RSA、舍入模式、泛型数学（NuGet 包内也附带同样的 LINQPad 示例）
- [Real：特性与函数](docs/real.md)
- [泛型数学与转换语义](docs/generic-math.md)
- [GMP / MPFR → `BigInt` / `Real` API 对照](docs/gmp-mapping.md)
- [性能](docs/performance.md)
- [设计、算法与验证](docs/design.md)
- [参考文献](docs/references.md)

## 微信群

![](https://io.starworks.cc:88/cv-public/2026/ocr-wxg-qr.png?1012)

如果微信群二维码过期了，请加入 QQ 群 [.NET骚操作 495782587](https://qm.qq.com/cgi-bin/qm/qr?_wv=1027&k=&authKey=&noverify=0&group_code=495782587)。

## 许可证

MIT。`Sdcb.Arithmetic` 是基于公开算法文献的独立 C# 实现，不包含任何 GMP 或 MPFR 源代码。
原生 GMP/MPFR 仅作为黑盒参照，用于测试比对和性能对比。所参考的文献见 [docs/references.md](docs/references.md)。
