// Correct rounding with BigFloat: what double cannot tell you.
//   1. 0.1 + 0.2 in double vs. BigFloat.
//   2. The five rounding modes.
//   3. Interval arithmetic: directed rounding gives guaranteed bounds.
//   4. Rump's example: double returns a wrong result and nothing warns you.
// Usage: dotnet run --project examples/Rounding
using System.Numerics;
using Sdcb.Arithmetic;

Console.WriteLine("1. 0.1 + 0.2");
Console.WriteLine($"   double           : {0.1 + 0.2:R}");
Console.WriteLine($"   BigFloat, 53 bits    : {BigFloat.Add(0.1, 0.2, 53)}   (same bits as double: BigFloat at 53 bits is IEEE arithmetic)");
BigFloat tenth = BigFloat.Parse("0.1", 200), fifth = BigFloat.Parse("0.2", 200);
Console.WriteLine($"   BigFloat, 200 bits   : {BigFloat.Add(tenth, fifth, 200).ToString(40)}   (printed to 40 digits; the error is below 1e-60)");

Console.WriteLine();
Console.WriteLine("2. -2/3 rounded to 8 bits in each mode");
foreach (RoundingMode mode in Enum.GetValues<RoundingMode>())
{
    BigFloat q = BigFloat.Divide(-2, 3, 8, mode);
    Console.WriteLine($"   {mode,-15}: {q,-12} = -{q.Mantissa}·2^{q.Exponent}");
}

Console.WriteLine();
Console.WriteLine("3. Harmonic number H(1000) at only 24 bits, with guaranteed bounds");
const int prec = 24;
BigFloat lo = 0, hi = 0, nearest = 0;
BigRational exact = 0;
for (int k = 1; k <= 1000; k++)
{
    lo = BigFloat.Add(lo, BigFloat.Divide(1, k, prec, RoundingMode.TowardNegative), prec, RoundingMode.TowardNegative);
    hi = BigFloat.Add(hi, BigFloat.Divide(1, k, prec, RoundingMode.TowardPositive), prec, RoundingMode.TowardPositive);
    nearest = BigFloat.Add(nearest, BigFloat.Divide(1, k, prec), prec);
    exact += new BigRational(1, k);
}
BigFloat exactBigFloat = BigFloat.FromRational(exact, 64);
Console.WriteLine($"   lower bound      : {lo.ToString(12)}");
Console.WriteLine($"   exact (rational) : {exactBigFloat.ToString(12)}");
Console.WriteLine($"   upper bound      : {hi.ToString(12)}");
Console.WriteLine($"   round-to-nearest : {nearest.ToString(12)}  (accumulated error, but which way?)");
Console.WriteLine($"   lo <= exact <= hi: {lo <= exactBigFloat && exactBigFloat <= hi}");

Console.WriteLine();
Console.WriteLine("4. Rump's example f(77617, 33096); true value -0.827396059946821368...");
Console.WriteLine($"   double           : {Rump(77617.0, 33096.0):R}");
foreach (int bits in (int[])[64, 100, 122, 128, 200])
{
    BigFloat.DefaultPrecision = bits;
    Console.WriteLine($"   BigFloat, {bits,3} bits   : {Rump<BigFloat>(77617, 33096).ToString(18)}");
}

// f(a, b) = 333.75·b⁶ + a²·(11·a²·b² − b⁶ − 121·b⁴ − 2) + 5.5·b⁸ + a / (2b), written once for any T.
static T Rump<T>(T a, T b) where T : IFloatingPoint<T>
{
    T b2 = b * b, b4 = b2 * b2, b6 = b4 * b2, b8 = b4 * b4, a2 = a * a;
    T c(double v) => T.CreateChecked(v);
    return c(333.75) * b6 + a2 * (c(11) * a2 * b2 - b6 - c(121) * b4 - c(2)) + c(5.5) * b8 + a / (c(2) * b);
}
