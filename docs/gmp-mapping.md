# GMP / MPFR → Sdcb.Arithmetic API mapping

Method names follow .NET conventions (`System.Numerics.BigInteger`, `Math`, generic math).
These tables help when porting code written against GMP or MPFR.

## `mpz_*` → `BigInt`

| GMP | Sdcb.Arithmetic |
|---|---|
| `mpz_add/sub/mul/tdiv_qr` | `+ - * / %`, `BigInt.DivRem(a, b, out q, out r)` |
| `mpz_{t,f,c}div_{q,r,qr}` | `BigInt.Divide/Remainder/DivRem(a, b, …, IntegerRounding)`; `BigInt.Mod` (floor remainder) |
| `mpz_mul_2exp`, `mpz_fdiv_q_2exp` | `<<`, `>>`, `ShiftLeft`, `ShiftRight(k)` (floor, same as `>>`) |
| `mpz_{t,f,c}div_{q,r}_2exp` | `ShiftRight(k, IntegerRounding)`, `RemainderPow2(k, IntegerRounding)`, `ModPow2(k)` |
| `mpz_divisible_2exp_p`, `mpz_congruent(_2exp)_p` | `IsDivisibleByPow2`, `IsCongruent` / `IsCongruentModPow2` |
| `mpz_divexact` | `BigInt.DivExact` |
| `mpz_cmp_d`, `mpz_get_d_2exp`, `mpz_fits_*` | `CompareTo(double)`, `ToDouble(out long exponent)`, `Fits*` |
| `mpz_and/ior/xor/com` | `& \| ^ ~` (infinite two's complement) |
| `mpz_tstbit/setbit/clrbit/combit` | `TestBit/SetBit/ClearBit/FlipBit` |
| `mpz_popcount/hamdist/scan0/scan1` | `PopCount/HammingDistance/IndexOfClearBit/IndexOfSetBit` |
| `mpz_pow_ui`, `mpz_powm` | `Pow`, `ModPow` |
| `mpz_sqrt(rem)`, `mpz_root(rem)` | `Sqrt(x)` / `Sqrt(x, out remainder)`, `Root(x, k)` / `Root(x, k, out remainder)` |
| `mpz_perfect_square_p/power_p` | `IsPerfectSquare/IsPerfectPower` |
| `mpz_gcd(ext)`, `mpz_lcm`, `mpz_invert` | `Gcd`, `ExtendedGcd`, `Lcm`, `ModInverse` / `TryModInverse` |
| `mpz_jacobi/legendre/kronecker` | `Jacobi/Legendre/Kronecker` |
| `mpz_probab_prime_p`, `mpz_nextprime` | `IsProbablePrime()` (bool), `GetPrimality()` (`Primality` enum), `NextPrime`/`PrevPrime` |
| `mpz_fac_ui`, `mpz_2fac_ui`, `mpz_mfac_uiui` | `Factorial`, `DoubleFactorial`, `MultiFactorial` |
| `mpz_primorial_ui`, `mpz_bin_ui(ui)` | `Primorial`, `Binomial` |
| `mpz_fib(2)_ui`, `mpz_lucnum(2)_ui` | `Fibonacci` / `FibonacciPair`, `Lucas` / `LucasPair` (tuples) |
| `mpz_get_str/set_str` (base 2..62) | `ToString(int)`, `Parse(s, fromBase)` (base 0 = auto prefix) |
| `mpz_import/export` | `FromMagnitudeBytes`, `ToMagnitudeBytes`; two's complement via `TryWriteLittleEndian/BigEndian`, `ReadLittleEndian/BigEndian` |
| `mpz_addmul`, `mpz_submul` | `BigInt.MultiplyAdd`, `BigInt.MultiplySubtract`; `BigInt.ModMultiply` for (a·b) mod m |
| `mpz_urandomb/urandomm` | `BigInt.Random/RandomBelow` (any `System.Random`) |

## `mpfr_*` → `Real`

Every function takes optional `precision` (bits) and `rounding` (`RealRounding`) arguments.

| MPFR | Sdcb.Arithmetic |
|---|---|
| `mpfr_add/sub/mul/div/sqr` | `Real.Add/Subtract/Multiply/Divide/Square` (and the operators) |
| `mpfr_fma/fms` | `Real.FusedMultiplyAdd/FusedMultiplySubtract` |
| `mpfr_sqrt/sqrt_ui/rec_sqrt/cbrt/rootn_si` | `Real.Sqrt` (also `ulong`), `ReciprocalSqrt`, `Cbrt`, `RootN` |
| `mpfr_pow/pow_si/pow_z/powr/compound_si` | `Real.Pow` (`Real`, `long` or `BigInt` exponent), `Powr`, `Compound` |
| `mpfr_log/log_ui/log2/log10/log1p` | `Real.Log` (also `ulong`), `Log2`, `Log10`, `LogP1` (+ `Log2P1`, `Log10P1`) |
| `mpfr_exp/exp2/exp10/expm1` | `Real.Exp`, `Exp2`, `Exp10`, `ExpM1` (+ `Exp2M1`, `Exp10M1`) |
| `mpfr_fmod/fmodquo` | `Real.TruncatedRemainder` (same as `%`), with `out long quotient` |
| `mpfr_remainder/remquo` | `Real.Ieee754Remainder`, with `out long quotient` |
| `mpfr_modf`, `mpfr_frexp`, `mpfr_frac` | `Real.Truncate(x, out fraction)`, `Real.Decompose(x, out exponent)`, `x.FractionalPart()` |
| `mpfr_rint/round/roundeven` | `x.RoundToInteger(RealRounding)`, `Real.Round(x, MidpointRounding)` |
| `mpfr_nextabove/nextbelow` | `Real.BitIncrement/BitDecrement` |
| `mpfr_dim` | `Real.PositiveDifference` |
| `mpfr_agm` | `Real.ArithmeticGeometricMean` |
| `mpfr_lngamma/lgamma` | `Real.LogGamma(x)` / `Real.LogGamma(x, out int sign)` |
| `mpfr_zeta/zeta_ui`, `mpfr_eint`, `mpfr_li2`, `mpfr_ai` | `Real.Zeta` (also `ulong`), `ExponentialIntegral`, `Dilogarithm`, `AiryAi` |
| `mpfr_j0/j1/jn`, `mpfr_y0/y1/yn` | `Real.BesselJ(n, x)`, `Real.BesselY(n, x)` |
| `mpfr_const_pi/log2/euler/catalan` | `Real.Pi`, `Ln2`, `EulerGamma`, `CatalanConstant` |
| `mpfr_urandom/nrandom/erandom` | `Real.Random`, `RandomNormal`, `RandomExponential` |
| `mpfr_set_d/get_d`, `mpfr_prec_round` | `Real.FromDouble` / `ToDouble`, `WithPrecision` |
