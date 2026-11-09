// One piece of generic code, many number types: BigInt, BigRational and BigFloat implement the
// .NET generic math interfaces, so they drop into the same algorithms as int/long/double.
// Usage: dotnet run --project examples/GenericMath
using System.Numerics;
using Sdcb.Arithmetic;
using Sdcb.Arithmetic.Examples;

PerfWarning.Print();

BigFloat.DefaultPrecision = 256;

Console.WriteLine("1. Factorial(30) with IBinaryInteger<T>");
Console.WriteLine($"   long   : {TryFactorial<long>(30)}");
Console.WriteLine($"   BigInt : {TryFactorial<BigInt>(30)}");

Console.WriteLine();
Console.WriteLine("2. Newton's method for √2 with INumber<T>: x ← (x + 2/x) / 2");
for (int i = 1; i <= 5; i++)
{
    Console.WriteLine($"   step {i}: BigRational {NewtonSqrt2<BigRational>(i),-28} double {NewtonSqrt2<double>(i):R}");
}
BigFloat r = NewtonSqrt2<BigFloat>(9);
Console.WriteLine($"   BigFloat after 9 steps : {r.ToString(70)}");
Console.WriteLine($"   error vs BigFloat.Sqrt : {BigFloat.Abs(r - BigFloat.Sqrt(2)).ToString(3)}");

Console.WriteLine();
Console.WriteLine("3. (1 + 1/n)^n → e for n = 10^15 with IPowerFunctions<T>");
Console.WriteLine($"   double : {CompoundE<double>(1e15):R}   (1 + 1e-15 already lost most of its digits)");
Console.WriteLine($"   BigFloat   : {CompoundE<BigFloat>(BigInt.Pow(10, 15)).ToString(20)}");
Console.WriteLine($"   e      : {BigFloat.Exp(1).ToString(20)}  ((1 + 1/n)^n itself is below e by ≈ e/(2n) ≈ 1.4e-15)");
Console.WriteLine();
Console.WriteLine("4. Conversions between types go through INumberBase<T>.CreateXxx");
BigFloat third = BigFloat.FromRational(new BigRational(1, 3));
Console.WriteLine($"   double.CreateChecked(BigFloat 1/3)       : {double.CreateChecked(third):R}");
Console.WriteLine($"   int.CreateSaturating(BigFloat +∞)        : {int.CreateSaturating(BigFloat.PositiveInfinity)}");
Console.WriteLine($"   BigInt.CreateTruncating(BigFloat 1e30+.5): {BigInt.CreateTruncating(BigFloat.Parse("1e30") + BigFloat.Half)}");

// n! computed with checked arithmetic: overflows for fixed-width integers, never for BigInt.
static string TryFactorial<T>(int n) where T : IBinaryInteger<T>
{
    try
    {
        T result = T.One;
        for (int k = 2; k <= n; k++) result = checked(result * T.CreateChecked(k));
        return result.ToString()!;
    }
    catch (OverflowException)
    {
        return "OverflowException";
    }
}

static T NewtonSqrt2<T>(int steps) where T : INumber<T>
{
    T two = T.CreateChecked(2), x = T.One;
    for (int i = 0; i < steps; i++) x = (x + two / x) / two;
    return x;
}

static T CompoundE<T>(T n) where T : IPowerFunctions<T>, INumber<T> => T.Pow(T.One + T.One / n, n);
