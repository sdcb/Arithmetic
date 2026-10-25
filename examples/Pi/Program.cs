// Computes π to N decimal digits (default 1,000,000) two ways and checks that they agree:
//   1. Real.Pi: the built-in correctly-rounded constant (Brent–Salamin AGM).
//   2. Chudnovsky series with binary splitting, written here using only BigInt.
// Usage: dotnet run --project examples/Pi -c Release -- [digits] [output-file]
using System.Diagnostics;
using Sdcb.Arithmetic;
using Sdcb.Arithmetic.Examples;

PerfWarning.Print();

int digits = args.Length > 0 ? int.Parse(args[0]) : 1_000_000;
string? outFile = args.Length > 1 ? args[1] : null;

// 1. Real.Pi: log2(10) ≈ 3.3219 bits per decimal digit, plus a few guard bits.
var sw = Stopwatch.StartNew();
int bits = (int)Math.Ceiling(digits * Math.Log2(10)) + 64;
Real pi = Real.Pi(bits);
string agm = pi.ToString(digits + 1, RealRounding.TowardZero).Replace(".", "")[..(digits + 1)];
Console.WriteLine($"{$"Real.Pi ({bits} bits)",-30}: {sw.ElapsedMilliseconds,6} ms");

// 2. Chudnovsky: π = 426880·√10005·Q / T, ~14.18 digits per term.
sw.Restart();
long terms = digits / 14 + 2;
int guard = 10; // extra digits so truncation errors cannot reach the last printed digit
(BigInt _, BigInt q, BigInt t) = Split(0, terms);
BigInt scale = BigInt.Pow(10, (ulong)(digits + guard));
BigInt sqrt10005 = BigInt.Sqrt(10005 * scale * scale);        // ⌊√10005 · 10^(digits+guard)⌋
BigInt piScaled = 426880 * sqrt10005 * q / t; // ≈ π · 10^(digits+guard)
string chud = piScaled.ToString()[..(digits + 1)];
Console.WriteLine($"{$"Chudnovsky ({terms} terms)",-30}: {sw.ElapsedMilliseconds,6} ms");

bool same = agm == chud;
Console.WriteLine($"{$"agree on {digits} digits",-30}: {same}");
Console.WriteLine($"π = {agm[0]}.{agm[1..51]}...");
Console.WriteLine($"digits {digits - 9}..{digits}: {agm[(digits - 9)..(digits + 1)]}");

if (outFile is not null)
{
    File.WriteAllText(outFile, $"{agm[0]}.{agm[1..]}");
    Console.WriteLine($"written to {outFile}");
}
return same ? 0 : 1;

// Binary splitting of the Chudnovsky series over terms [a, b): returns (P, Q, T).
static (BigInt P, BigInt Q, BigInt T) Split(long a, long b)
{
    if (b - a == 1)
    {
        if (a == 0)
        {
            return (BigInt.One, BigInt.One, 13591409);
        }
        BigInt p = (BigInt)(6 * a - 5) * (2 * a - 1) * (6 * a - 1);
        BigInt q = (BigInt)a * a * a * 10939058860032000L;
        BigInt t = p * (13591409 + 545140134 * (BigInt)a);
        return (p, q, (a & 1) == 1 ? -t : t);
    }
    long m = (a + b) / 2;
    (BigInt p1, BigInt q1, BigInt t1) = Split(a, m);
    (BigInt p2, BigInt q2, BigInt t2) = Split(m, b);
    return (p1 * p2, q1 * q2, t1 * q2 + p1 * t2);
}
