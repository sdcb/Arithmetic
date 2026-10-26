# Performance

## Setup

Measured on win-x64, AMD Ryzen 7 5800X (Zen 3, AVX2 / `Vector256` tier), .NET 10,
against **native GMP 6.3** (`Sdcb.Arithmetic.Gmp.runtime.win-x64`, fat build with
Zen-tuned assembly) and `System.Numerics.BigInteger`. Ratios are from the `perf` mode
in `tests/Sdcb.Arithmetic.Benchmarks` (best of 5 rounds, single run, so expect roughly ±10–30% noise between runs;
higher = slower than GMP; sizes from 1 to 65536 limbs, limbs = 64-bit words, N = limb count).

## BigInt vs native GMP 6.3 and BigInteger

The `BigInteger` column is the speedup of `BigInt` over `System.Numerics.BigInteger`
(BigInteger time / ours, **higher = faster than BigInteger**), from the same `perf` run;
"—" means `BigInteger` has no equivalent operation.

| operation | vs native GMP (ours / GMP) | vs `BigInteger` (speedup) |
|---|---|---|
| add / sub | ~1.1–2.3× | ~2.5–3× (64–256 limbs); ~3.4–5× (≥ 512) |
| shl / shr | ~1.2–3.1× (worst around 256 limbs) | ~3.9–7.6× / ~5.7–11× |
| xor | ~1.0–1.6× | ~4.7–7.3× |
| popcount | ~0.2–0.43× (faster) | ~1.5–3.1× |
| mul | ~1.3–1.6× from 16 limbs up (~1.2–1.4× at 65536); ~1.2–2.3× below 16 limbs | ~1.9–2.8× (≤ 16 limbs); ~3.7–5× (32–2048); ~7× / ~15× / ~28× at 4096 / 16384 / 65536 |
| sqr | ~1.6–2.1× from 32 limbs up; ~2.6× at 8–16 limbs | ~1.8–2.3× (≤ 16); ~3.3–5.3× (32–2048); ~8× / ~18× / ~36× at 4096 / 16384 / 65536 |
| divrem 2N/N | ~1.5–2.6× | ~1.1–2.1× (≤ 8 limbs); ~2.2–2.7× (16–64); ~3.1–4.3× (128–4096); ~5.9× at 16384 |
| divexact | ~1.7–2.4×; ~2.9× at 16384 limbs | ~1.1–2.5× (≤ 8, vs `/`); ~4–5.3× (16–4096); ~6× at 16384 |
| gcd | ~0.9–1.4× | ~1.5–2.4× (≤ 512 limbs); ~2.8–5.8× (1024–4096) |
| gcdext | ~1.1–2.1× | — |
| sqrt 2N | ~1.5–2.1× (≤ 16 limbs); ~2.1–3.4× (32–16384) | — |
| powmod (4–16 / 32–64 limbs) | ~2.5–2.8× / ~1.3–1.6× | ~4.2–7× / ~5.3–6.2× |
| jacobi / invert | ~2.3–3.1× / ~1.3–1.7× | — |
| factorial / binomial / fibonacci | ~1.9–2× / ~1.9–3.4× / ~1.5–2.8× | — |
| ToString base 10 | ~0.5–1× (≤ 32 limbs); ~1.35–1.9× (≥ 64) | ~0.9× (≤ 4 limbs, slower); ~1.1–2× (8–64); ~2.9–16× (128–4096) |
| parse base 10 | ~0.6–0.95× (≤ 128 limbs); ~1.1–1.4× (≥ 256) | ~1.4–1.8× (≤ 64 limbs); ~2.2–4× (128–4096) |

Note that the native side writes into preallocated `mpz_t` destinations, while
`BigInt` is immutable and allocates every result, so small-operand ratios include
allocation/GC cost.

Native GMP's inner loops use `MULX`/`ADCX`/`ADOX` dual carry chains that cannot
be expressed in managed code, so the scalar row kernels (`addmul_1`, schoolbook
division, REDC rows) stay ~2–2.5× behind. Multiplication avoids that gap with an
AVX2 basecase: 32×16-bit `vpmuludq` products accumulate carry-free in 64-bit
lanes. Most operations land at ~1.2–2.5× of GMP, and a few small or mid-size
kernels reach ~3×.

## Real vs native MPFR 4.2

All results are correctly rounded on both sides, verified bit-for-bit by the
`exact` harness (mantissa + exponent + NaN/±inf/±0, all five rounding modes,
precisions 2–640, 56,370 checks, 0 failures). Times vs `libmpfr-6` (same machine
and `perf` mode as above):

| function | 53 bits | 1024 bits | 4096 bits | 16384 bits |
|---|---|---|---|---|
| add | ~1.3× | ~2× | ~2.2× | ~2.5× |
| mul | ~1.8× | ~1.5× | ~1.6× | ~1.6× |
| div | ~2× | ~3× | ~2.6–3.3× | ~1.8–2.1× |
| sqrt | ~1.3× | ~1.9× | ~2.7–3.1× | ~2.4–2.9× |
| exp | ~0.9× | ~1.5× | ~1.4× | — |
| sin | ~1.35× | ~1.5× | ~1.3× | — |
| log | ~1.8× | ~2.3× | ~2.7× (AGM) | — |
| atan | ~0.9× | ~1× | ~1.8× | — |
| pow | ~1.6× | ~2.2× | ~2.1× | — |
| decimal ToString | ~0.2× | ~1.1× | ~2× | — |

(— : `perf` only times the transcendental rows up to 4096 bits.)

`Real` runs its rounding core, ring operations, fixed-point helpers and series
accumulators directly on the internal limb kernels with pooled/stack scratch, so each
operation allocates only its result. Above ~4096 bits the remaining gap is the
integer kernels themselves (multiplication/division), not `Real` overhead.

π via Chudnovsky binary splitting: 262144 bits in ~0.1 s and 1M decimal digits (3.3M bits) in ~0.8 s (first call in a fresh process,
including JIT). The `r.pi` row printed by `perf` is not a like-for-like ratio: `Real.Pi` reuses its
internal constant cache while the MPFR side frees its cache on every call.

## Running the benchmarks

```
dotnet run --project tests/Sdcb.Arithmetic.Benchmarks -c Release -- perf           # timing table vs GMP + MPFR (+ BigInteger)
dotnet run --project tests/Sdcb.Arithmetic.Benchmarks -c Release -- gmp-exact      # bit-exact vs native GMP (7637 checks)
dotnet run --project tests/Sdcb.Arithmetic.Benchmarks -c Release -- exact          # bit-exact vs native MPFR (56370 checks)
dotnet run --project tests/Sdcb.Arithmetic.Benchmarks -c Release -- exact big      # + 4097/9000-bit precisions (60540 checks, slow)
```

`perf` accepts an op filter (`perf "=mul,=sqr,r.div"`, where `=` means an exact match)
and `PERF_SIZES=64,256` limits the sizes. Integer rows also time `BigInteger`;
set `PERF_BCL=0` to skip those columns. The `Perf` GitHub workflow
(`.github/workflows/perf.yml`, run manually or weekly) prints the same table, plus the
runner's CPU and vector tier, to the job summary. It only reports and gates nothing.
