<Query Kind="Program">
  <NuGetReference Prerelease="true">Sdcb.Arithmetic</NuGetReference>
  <Namespace>Sdcb.Arithmetic</Namespace>
  <Namespace>System.Numerics</Namespace>
</Query>

// Correct rounding with Real: what double cannot tell you.
//   1. 0.1 + 0.2 in double vs. Real.
//   2. The five rounding modes.
//   3. Interval arithmetic: directed rounding gives guaranteed bounds.
//   4. Rump's example: double returns a wrong result and nothing warns you.
void Main()
{
	Console.WriteLine("1. 0.1 + 0.2");
	Console.WriteLine($"   double           : {0.1 + 0.2:R}");
	Console.WriteLine($"   Real, 53 bits    : {Real.Add(0.1, 0.2, 53)}   (same bits as double: Real at 53 bits is IEEE arithmetic)");
	Real tenth = Real.Parse("0.1", 200), fifth = Real.Parse("0.2", 200);
	Console.WriteLine($"   Real, 200 bits   : {Real.Add(tenth, fifth, 200).ToString(40)}   (printed to 40 digits; the error is below 1e-60)");

	Console.WriteLine();
	Console.WriteLine("2. -2/3 rounded to 8 bits in each mode");
	foreach (RealRounding mode in Enum.GetValues<RealRounding>())
	{
		Real q = Real.Divide(-2, 3, 8, mode);
		Console.WriteLine($"   {mode,-15}: {q,-12} = -{q.Mantissa}·2^{q.Exponent}");
	}

	Console.WriteLine();
	Console.WriteLine("3. Harmonic number H(1000) at only 24 bits, with guaranteed bounds");
	const int prec = 24;
	Real lo = 0, hi = 0, nearest = 0;
	BigRational exact = 0;
	for (int k = 1; k <= 1000; k++)
	{
		lo = Real.Add(lo, Real.Divide(1, k, prec, RealRounding.TowardNegative), prec, RealRounding.TowardNegative);
		hi = Real.Add(hi, Real.Divide(1, k, prec, RealRounding.TowardPositive), prec, RealRounding.TowardPositive);
		nearest = Real.Add(nearest, Real.Divide(1, k, prec), prec);
		exact += new BigRational(1, k);
	}
	Real exactReal = Real.FromRational(exact, 64);
	Console.WriteLine($"   lower bound      : {lo.ToString(12)}");
	Console.WriteLine($"   exact (rational) : {exactReal.ToString(12)}");
	Console.WriteLine($"   upper bound      : {hi.ToString(12)}");
	Console.WriteLine($"   round-to-nearest : {nearest.ToString(12)}  (accumulated error, but which way?)");
	Console.WriteLine($"   lo <= exact <= hi: {lo <= exactReal && exactReal <= hi}");

	Console.WriteLine();
	Console.WriteLine("4. Rump's example f(77617, 33096); true value -0.827396059946821368...");
	Console.WriteLine($"   double           : {Rump(77617.0, 33096.0):R}");
	int saved = Real.DefaultPrecision;
	foreach (int bits in (int[])[64, 100, 122, 128, 200])
	{
		Real.DefaultPrecision = bits;
		Console.WriteLine($"   Real, {bits,3} bits   : {Rump<Real>(77617, 33096).ToString(18)}");
	}
	Real.DefaultPrecision = saved;
}

// f(a, b) = 333.75·b⁶ + a²·(11·a²·b² − b⁶ − 121·b⁴ − 2) + 5.5·b⁸ + a / (2b), written once for any T.
static T Rump<T>(T a, T b) where T : IFloatingPoint<T>
{
	T b2 = b * b, b4 = b2 * b2, b6 = b4 * b2, b8 = b4 * b4, a2 = a * a;
	T c(double v) => T.CreateChecked(v);
	return c(333.75) * b6 + a2 * (c(11) * a2 * b2 - b6 - c(121) * b4 - c(2)) + c(5.5) * b8 + a / (c(2) * b);
}
