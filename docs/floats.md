# Fixed-width floats: Float128, DoubleDouble, QuadDouble

Three allocation-free `readonly struct` types for when `double` is not precise enough but
arbitrary precision (`BigFloat`) is more than you need.

| | `Float128` | `DoubleDouble` | `QuadDouble` |
|---|---|---|---|
| Representation | IEEE 754 binary128 bit pattern (`UInt128`) | hi + lo, two doubles | x0 + x1 + x2 + x3, four doubles |
| Precision | 113 bits | ~106 bits | ~212 bits |
| Exponent range | 2^-16494 … 2^16384 | double's (precision degrades below ~2^-969) | double's (precision degrades below ~2^-850) |
| Rounding | Correctly rounded, five `RoundingMode` modes | Approximate (error bounds below) | Approximate (error bounds below) |
| Generic math | `IBinaryFloatingPointIeee754<T>` | `IFloatingPoint<T>` + function interfaces | `IFloatingPoint<T>` + function interfaces |
| Models | IEEE 754-2019, C `_Float128` | QD library `dd_real` | QD library `qd_real` |

Pick `Float128` when results must be reproducible bit for bit or exchanged with other IEEE
binary128 implementations. Pick `DoubleDouble`/`QuadDouble` for speed: their transcendental
functions are native and several times faster than `BigFloat` or MPFR at the same precision.

## Float128

- `+ - * /`, `Sqrt`, `FusedMultiplyAdd`, `%` (exact fmod), `Ieee754Remainder`, conversions and
  all rounding to integers are correctly rounded with subnormals, ±0, ±∞ and NaN. Operators use
  nearest-even; `Add`, `Subtract`, `Multiply`, `Divide`, `Sqrt`, `FusedMultiplyAdd` and the
  `FromXxx`/`ToXxx` conversions take a `RoundingMode`.
- Core arithmetic is pure integer code (unpack, compute with guard/sticky bits, one shared
  round-and-pack), structured after Berkeley SoftFloat's f128.
- `ToString()` prints the shortest string that round-trips; `Parse` is correctly rounded from
  the exact decimal value.
- Transcendental functions (34, each with a `RoundingMode` overload) are correctly rounded.
  Exp, Exp2, Log, Log2, Log10, Pow, Sin, Cos, Tan, Asin, Acos, Atan, Atan2, Sinh, Cosh and Tanh
  first evaluate a triple-double (~159-bit) kernel with a proven error bound and round directly
  when the error interval cannot straddle a rounding boundary (Ziv's strategy). Otherwise —
  exact results, midpoints, subnormal or near-overflow results, trig arguments ≥ 2^30, about
  0.02–0.04% of random trig inputs — and for the remaining functions, they are evaluated through
  `BigFloat` and rounded once.
- NaN payloads are not preserved: arithmetic returns the canonical quiet NaN.
- Verified bit for bit against `BigFloat` (round-to-odd oracle, so no double rounding) on ~180,000
  random add/sub/mul/div/sqrt/FMA cases across all five rounding modes, biased toward
  subnormals, overflow, large exponent gaps, cancellation and ties.

## DoubleDouble and QuadDouble

Algorithms from Hida, Li and Bailey's QD library: accurate ("IEEE") addition, FMA-based
multiplication, division and square root (`QuadDouble` uses digit-by-digit division and square
root with an error-free remainder expansion and a single final renormalization). Exp, log, trig and hyperbolic
functions are native (argument reduction + tables + series); arguments that need very precise
reduction (huge trig arguments, `%`, `Round(x, digits)`) fall back to `BigFloat` for correctness.

Worst relative error measured against `BigFloat` at 300/500 bits on seeded random inputs:

| | `DoubleDouble` | `QuadDouble` |
|---|---|---|
| add / sub / mul | 2^-104.5 | 2^-214 |
| div / sqrt | 2^-103.8 / 2^-104.8 | 2^-215 / 2^-215 |
| exp / log families | 2^-104.3 or better | 2^-211.4 or better |
| trig / inverse trig / hyperbolic | 2^-103.4 or better | 2^-210.4 or better |
| pow | 2^-101.2 (\|y·log x\| ≤ 16) | 2^-209.7 |

`Pow` is exp(y·log x), so its error grows with \|y·log x\|. `DoubleDouble.FusedMultiplyAdd` is
not fused (the product is rounded to double-double first).

`Epsilon` is the smallest positive value (2^-1074, as for `double`); `MachineEpsilon` is the
relative precision (2^-104 for `DoubleDouble`, 2^-209 for `QuadDouble`).

## Conversions between the three

| From → To | Behavior |
|---|---|
| `DoubleDouble` → `QuadDouble` | implicit, exact |
| `QuadDouble` → `DoubleDouble` | explicit, leading ~106 bits |
| `DoubleDouble` → `Float128` | explicit, correctly rounded (hi + lo rounds once) |
| `QuadDouble` → `Float128` | explicit, correctly rounded (via exact `BigFloat`) |
| `Float128` → `QuadDouble` | explicit, exact within double's exponent range |
| `Float128` → `DoubleDouble` | explicit, leading ~106 bits |

All three convert exactly to `BigFloat` (implicit) and from `BigFloat`, `BigInt`, `BigRational`,
`decimal` and every primitive numeric type.

## Performance

See [performance.md](performance.md#fixed-width-floats).
