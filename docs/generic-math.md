# Generic math

All three types plug into .NET generic math:

| Type | Interfaces |
|---|---|
| `BigInt` | `IBinaryInteger<BigInt>`, `ISignedNumber<BigInt>` (like `BigInteger`) |
| `BigRational` | `INumber<BigRational>`, `ISignedNumber<BigRational>` |
| `Real` | `IBinaryFloatingPointIeee754<Real>` (all exp/log/trig/hyperbolic/root/power function interfaces) |

All three are also `IUtf8SpanFormattable` / `IUtf8SpanParsable` and parse with `NumberStyles`.

```csharp
static T Horner<T>(T x, params int[] c) where T : INumberBase<T>
{
    T r = T.Zero;
    foreach (int k in c) r = r * x + T.CreateChecked(k);
    return r;
}

BigInt     i = Horner((BigInt)10, 3, -7, 11);             // exact
BigRational q = Horner(new BigRational(2, 3), 3, -7, 11); // exact: 23/3
Real       r = Horner(Real.FromDouble(0.1, 200), 3, -7, 11);
double     d = double.CreateChecked(r);                   // correctly rounded
int        n = int.CreateSaturating(Real.PositiveInfinity); // int.MaxValue
```

## Semantics

- **Real precision in generic code**:
  - Unary functions (`T.Sin(x)`, `T.Sqrt(x)`, …) round to `x`'s precision.
  - Binary operators round to the larger operand precision.
  - Constants (`T.Zero`, `T.One`, `T.Pi`, `T.E`, `T.Tau`, `T.Epsilon`) use `Real.DefaultPrecision`.
  - Every result is correctly rounded to nearest-even.
  - The overloads with precision/rounding parameters (`Real.Sin(x, 512, RealRounding.TowardZero)`) are still available.
  - `Real.Log(x, 10)` binds to the *precision* overload. Write `Real.Log(x, (Real)10)` (or call `T.Log(x, b)` in generic code) for an arbitrary base.
- **Real ↔ primitives**:
  - Implicit conversions are exact. Integers, `BigInt` and `BigInteger` keep their full bit length; `Half`/`float`/`double` keep their significand width.
  - `decimal` and `BigRational` convert explicitly, rounded to `DefaultPrecision`.
  - Conversions to `double`/`float`/`Half`/`decimal` are correctly rounded, including subnormals.
  - IEEE binary128 (quad precision, e.g. C `__float128`) is exchanged as raw bits: `Real.FromBinary128Bits(UInt128)` is exact by default (precision `max(DefaultPrecision, 113)`), and `ToBinary128Bits(rounding)` is correctly rounded, including subnormals and overflow. NaN payloads are not preserved.
  - `CreateSaturating`/`CreateTruncating` follow the BCL rules for `double` (NaN → 0, ±∞ → Min/Max).
- **`Real.Min`/`Max`** keep the MPFR semantics: a number wins over NaN, the same as `MinNumber`/`MaxNumber`. `INumber<Real>.Min`/`Max` propagate NaN like `double.Min`/`Max`.
- **`Real.Epsilon`** is `BitIncrement(Zero)` at `DefaultPrecision`. Real exponents are otherwise unbounded, so this is a convention, not an underflow threshold.
- **Real bitwise operators** (`& | ^ ~`, required by `IBinaryNumber`) exist only through the interface:
  - `&`, `|` and `^` act on the exact binary expansions of the magnitudes and on the sign bits.
  - `~` complements the significand within its precision.
- **BigInt fixed-width operations** (`LeadingZeroCount`, `PopCount` of negatives, `RotateLeft/Right`, `>>>`) use the minimal two's-complement form rounded up to whole 64-bit limbs. This mirrors `BigInteger`, which uses 32-bit words.
- **BigInt bytes**: `GetByteCount`/`TryWriteLittleEndian` write the shortest signed two's-complement form, the same bytes as `BigInteger.ToByteArray()`.
- **BigInt hex/binary `NumberStyles` parsing** is two's complement (`"FF"` → -1), like `BigInteger`.
- **`BigInteger` interop**:
  - Implicit `BigInteger → BigInt`; explicit `BigInt → BigInteger`, so mixed expressions stay unambiguous.
  - `BigRational` and `Real` also accept `BigInteger` implicitly.
- **Extras**:
  - `BigRational.Floor/Ceiling/Truncate/Round` (with `MidpointRounding` and decimal digits) and `%`.
  - `Real.Round(x, digits, mode)` rounds the exact value to decimal digits.
  - Fused `BigInt.MultiplyAdd`/`MultiplySubtract` (one result buffer) and `BigInt.ModMultiply`.
