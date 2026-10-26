# Real: correctly-rounded floating point

`Real` is the managed equivalent of MPFR's `mpfr_t`: an immutable binary floating-point
number with a per-value precision in bits.

## Features

- **Correct rounding, always**: every operation returns the representable value
  nearest the exact mathematical result under the requested mode, verified
  bit-for-bit against native MPFR (the `exact` harness, see [performance.md](performance.md)).
- **Five rounding modes**: NearestEven, TowardZero, TowardPositive,
  TowardNegative, AwayFromZero.
- **Full IEEE special values**: NaN, ±infinity, signed zero with the MPFR
  special-case tables for every function.
- **Function coverage** (all correctly rounded; see [gmp-mapping.md](gmp-mapping.md) for the MPFR names):
  - arithmetic: `Add`, `Subtract`, `Multiply`, `Divide`, `Square`, `FusedMultiplyAdd`/`FusedMultiplySubtract`,
    single-rounding `Sum`/`Dot`, `Hypot`, `PositiveDifference`, `Min`/`Max`, `CopySign`;
  - roots and powers: `Sqrt`, `ReciprocalSqrt`, `Cbrt`, `RootN`, `Pow` (`Real`, `long` or `BigInt` exponent),
    `Powr`, `Compound`, with algebraic exactness detection;
  - exponentials and logarithms: `Exp`, `Exp2`, `Exp10`, `ExpM1`, `Exp2M1`, `Exp10M1`, `Log`, `Log2`, `Log10`,
    `LogP1`, `Log2P1`, `Log10P1`, logarithm to an arbitrary base;
  - trigonometric: `Sin`/`Cos`/`Tan`/`SinCos`, `Asin`/`Acos`/`Atan`/`Atan2`, `Sec`/`Csc`/`Cot`, the `…Pi` variants,
    degree/radian conversion; hyperbolic counterparts and their inverses;
  - special functions: `Gamma`, `LogGamma`, `Digamma`, `Beta`, `Erf`/`Erfc`, `Zeta`, `ExponentialIntegral`,
    `Dilogarithm`, `BesselJ`/`BesselY`, `AiryAi`, `ArithmeticGeometricMean`, `Factorial`;
  - integer and remainder operations: `Floor`/`Ceiling`/`Truncate`/`Round`/`RoundToInteger`, `FractionalPart`,
    `TruncatedRemainder` (`%`), `Ieee754Remainder`, `Decompose`, `BitIncrement`/`BitDecrement`;
  - constants `Pi`, `E`, `Ln2`, `EulerGamma`, `CatalanConstant`; random `Random`/`RandomNormal`/`RandomExponential`;
    decimal formatting and parsing.
- **Ziv's onion-peeling** with per-function error bounds, a guarded
  cancellation detector (intermediate magnitudes are tracked so heavy
  cancellation automatically escalates working precision), and algebraic
  exactness detection (`2^10 = 1024`, `log10(1000) = 3`, `pow(4, 0.5) = 2`
  return instantly and exactly).
- **Asymptotically fast transcendentals**: rectangular-splitting series for
  exp/sin/atan, AGM logarithm above `SDCB_ARITHMETIC_LOGAGM` bits (default 700),
  Chudnovsky π and ln 2 / e by binary splitting above
  `SDCB_ARITHMETIC_CONSTAGM` bits (default 2000) — 262144-bit π in ~0.1 s, 1M decimal digits in ~0.8 s.

## Example

```csharp
using Sdcb.Arithmetic;

Real.DefaultPrecision = 256;                       // bits
Real pi = Real.Pi();                               // correctly rounded π
Real x = Real.Parse("1.5");
Real y = Real.Sin(x);                              // correctly rounded
Real z = Real.Exp(y, precision: 512, RealRounding.TowardPositive);
Console.WriteLine(z);                              // round-trip decimal output
double d = z.ToDouble(RealRounding.NearestEven);   // correctly rounded projection
```
