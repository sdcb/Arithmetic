# Migrating from 1.x to 2.0

## Renames

| 1.x | 2.0 |
|---|---|
| `Real` | `BigFloat` |
| `RealRounding` | `RoundingMode` |

Members, behavior and the binary layout are otherwise unchanged. `BigFloat` joins `BigInt` and
`BigRational` as the arbitrary-precision family; the new fixed-width types are `Float128`,
`DoubleDouble` and `QuadDouble` (see [floats.md](floats.md)).

## Quick migration

Either rename the uses (a whole-word replace of `RealRounding` → `RoundingMode`, then
`Real` → `BigFloat`), or keep existing code compiling with two aliases in any one file of
your project (C# 12+):

```csharp
global using Real = Sdcb.Arithmetic.BigFloat;
global using RealRounding = Sdcb.Arithmetic.RoundingMode;
```

## New in 2.0

- `Float128`: IEEE 754 binary128, correctly rounded in all five rounding modes.
- `DoubleDouble` / `QuadDouble`: fast ~106 / ~212-bit floating point (QD library algorithms).
