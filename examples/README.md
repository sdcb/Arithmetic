# Examples

Each folder is a small console app that references the library project directly.
The timing examples print a warning when started from a Debug build or with a debugger attached;
use `-c Release` for representative numbers.

| Example | Shows | Run |
|---|---|---|
| [Pi](Pi/Program.cs) | 1,000,000 digits of π two ways: `Real.Pi` (AGM) and a Chudnovsky binary-splitting series written with `BigInt`, checked against each other | `dotnet run --project examples/Pi -c Release -- [digits] [output-file]` |
| [Rsa](Rsa/Program.cs) | Textbook RSA from `BigInt` primitives (`NextPrime`, `ModInverse`, `ModPow`, CRT) and PKCS#1 v1.5 signatures verified by .NET's `RSA` in both directions (for learning only) | `dotnet run --project examples/Rsa -c Release -- [bits]` |
| [Rounding](Rounding/Program.cs) | 0.1 + 0.2, the five rounding modes, interval bounds via directed rounding, and Rump's example where `double` is silently wrong | `dotnet run --project examples/Rounding` |
| [GenericMath](GenericMath/Program.cs) | The same generic algorithms on `long`/`double` and `BigInt`/`BigRational`/`Real` | `dotnet run --project examples/GenericMath` |

## LINQPad samples

The same four examples ship inside the NuGet package as LINQPad queries
([src/Sdcb.Arithmetic/linqpad-samples](../src/Sdcb.Arithmetic/linqpad-samples)). In LINQPad 9+ (the
package targets .NET 10), add the `Sdcb.Arithmetic` NuGet package to a query and open the
**Samples** tab: they appear under *NuGet → Sdcb.Arithmetic*.
