# Design, algorithms and verification

## Runtime and hardware

- **Pure managed code**: targets .NET 10, AOT-friendly, no P/Invoke.
- **Automatic ISA fallback**: every vectorized kernel cascades
  `Vector512 → Vector256 → Vector128 → scalar` with explicit capability checks
  (Avx512F/Avx2/Sse41/AdvSimd). Wider tiers process the bulk and narrower tiers
  take the tail. Tier checks are JIT-time constants, so unused tiers are
  compiled out. Lower the tier with the standard .NET knobs, e.g.
  `DOTNET_PreferredVectorBitWidth=128`, `DOTNET_EnableAVX2=0` or
  `DOTNET_EnableHWIntrinsic=0`.
- **Hardware-accelerated on x64**: `MULX` via `Math.BigMul`, hardware 128/64
  division via `X86Base.X64.DivRem`, vectorized add/sub/shift/bitwise/popcount,
  with portable fallbacks everywhere (ARM64, WASM and x86 all work).

## Algorithms

Citations for each algorithm are in [references.md](references.md).

- **Multiplication**: basecase (AVX2 `vpmuludq` kernel over 32×16-bit digit
  products, scalar product scanning elsewhere) → Karatsuba → Toom-3 → SIMD NTT
  (three-prime convolution + Garner CRT) for large operands, with dedicated
  squaring paths.
- **Division**: reciprocal (`invert_limb`) preinversion + Knuth algorithm D,
  Burnikel-Ziegler divide-and-conquer, and a remainder-free approximate quotient
  for `Real`. Exact division uses Jebelean's right-to-left scheme, a
  divide-and-conquer variant and Newton-Hensel.
- **GCD**: subquadratic HGCD (half-GCD) over Lehmer's 62-bit inner loop, with
  canonical cofactor normalization.
- **Square root**: Karatsuba recursive `mpn_sqrtrem`, plus a remainder-free
  approximate root for `Real`.
- **Modular exponentiation**: Montgomery REDC with sliding windows; large moduli
  reduce with two multiplications instead of a row-by-row REDC.
- **Primality**: trial division + BPSW + extra Miller-Rabin rounds; next/previous
  prime via an incremental small-prime sieve.
- **Base conversion**: divide-and-conquer in both directions, bases 2–62, with a
  cached power table.
- **Combinatorics**: factorial via prime-swing, binomial via Kummer's method,
  primorial via product trees, Fibonacci/Lucas via fast doubling.

Algorithm thresholds can be overridden with `SDCB_ARITHMETIC_*` environment
variables for tuning (see `Core/Nat.Mul.cs` and friends).

## Verification

- 90+ property-based test suites cross-check every `BigInt` operation against
  `System.Numerics.BigInteger` on randomized inputs (including carry-chain
  stress patterns and Karatsuba/Toom-3/NTT-sized operands).
- `gmp-exact` cross-verifies every `mpz` operation against native GMP
  (7,637 checks: exact agreement across randomized sizes and signs).
- `exact` cross-verifies every `Real` function against native MPFR
  (56,370 checks: mantissa + exponent + kind, all five rounding modes,
  precisions 2–640).
- At 53-bit precision, `Real`'s +, −, ×, ÷, sqrt and fma agree bit-for-bit with
  IEEE double arithmetic on thousands of random inputs.
- Every SIMD tier computes identical results. CI runs the suite at the default
  tier, with `DOTNET_PreferredVectorBitWidth=512`/`256`/`128` and with
  `DOTNET_EnableHWIntrinsic=0`, on x64 and ARM64, plus the native GMP/MPFR
  cross-checks on Windows.

## Project layout

```
src/Sdcb.Arithmetic/
  Core/Nat.*.cs                          low-level natural-number limb kernels
  Core/MontgomeryContext.cs              REDC engine for ModPow
  Core/Isa.cs                            vector tier checks (JIT-time constants)
  Integers/BigInt.*.cs                   mpz equivalent (partial struct)
  Rationals/BigRational.cs               mpq equivalent
  Reals/Real.cs                          mpfr equivalent: representation + rounding core
  Reals/Real.Arithmetic|Integer|Convert  ring ops, rint family, decimal I/O
  Reals/Real.Ziv.cs                      Ziv loops + fixed-point helpers
  Reals/Real.Kernels.cs                  limb-level rounding/ring/fixed-point kernels on Core/Nat
  Reals/Real.Constants|Exp|Trig|Hyper|Pow  transcendental functions
tests/Sdcb.Arithmetic.Tests/             xunit suite
tests/Sdcb.Arithmetic.Benchmarks/        BenchmarkDotNet + perf/exact harnesses vs native GMP/MPFR
examples/                                runnable console demos (π, RSA, rounding, generic math)
```
