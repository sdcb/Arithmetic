using System.Diagnostics;
using System.Numerics;
using Sdcb.Arithmetic.Core;

namespace Sdcb.Arithmetic;

public readonly partial struct BigInt
{
    // ---------------------------------------------------------------
    // Factorials & binomials via prime factorization (no big divisions)
    // ---------------------------------------------------------------

    /// <summary>Largest n handled by the prime-factorization paths (bounds the sieve to ~n/16 bytes).</summary>
    private const ulong PrimeMethodLimit = 1UL << 31;

    /// <summary>n! (mpz_fac_ui): Luschny's prime-swing recursion on the odd part, then one shift.</summary>
    public static BigInt Factorial(ulong n)
    {
        if (n <= 20)
        {
            return new BigInt(SmallFactorial(n));
        }
        if (n > PrimeMethodLimit)
        {
            throw new ArgumentOutOfRangeException(nameof(n), "Factorial argument too large.");
        }
        uint[] primes = SmallPrimes.PrimesUpTo((uint)n);
        // v2(n!) = n - popcount(n).
        return OddFactorial(n, primes).ShiftLeft((long)(n - (ulong)BitOperations.PopCount(n)));
    }

    private static ulong SmallFactorial(ulong n)
    {
        ulong f = 1;
        for (ulong i = 2; i <= n; i++)
        {
            f *= i;
        }
        return f;
    }

    /// <summary>The odd part of n!: oddfac(n) = oddfac(n/2)^2 * oddswing(n).</summary>
    private static BigInt OddFactorial(ulong n, uint[] primes)
    {
        if (n <= 20)
        {
            ulong f = SmallFactorial(n);
            return new BigInt(f >> BitOperations.TrailingZeroCount(f));
        }
        BigInt half = OddFactorial(n / 2, primes);
        return half.Square() * OddSwing(n, primes);
    }

    /// <summary>n! / (floor(n/2)!)^2 without its powers of two: prod p^(Σ_k floor(n/p^k) mod 2) over odd primes p &lt;= n.</summary>
    private static BigInt OddSwing(ulong n, uint[] primes)
    {
        var words = new WordProduct();
        ulong sqrtN = SmallPrimes.Sqrt64(n);
        for (int i = 1; i < primes.Length && primes[i] <= n; i++)
        {
            ulong p = primes[i];
            if (p <= sqrtN)
            {
                ulong pe = 1, q = n;
                while (q >= p)
                {
                    q /= p;
                    if ((q & 1) != 0)
                    {
                        pe *= p;
                    }
                }
                words.Add(pe);
            }
            else if (((n / p) & 1) != 0)
            {
                words.Add(p);
            }
        }
        return words.Result();
    }

    /// <summary>Accumulates a product of word-sized factors, packing them into full limbs before a balanced product tree.</summary>
    private struct WordProduct
    {
        private List<ulong>? _list;
        private ulong _acc;

        public void Add(ulong f)
        {
            if (f <= 1)
            {
                return;
            }
            if (_acc == 0)
            {
                _acc = f;
                return;
            }
            ulong hi = Math.BigMul(_acc, f, out ulong lo);
            if (hi == 0)
            {
                _acc = lo;
                return;
            }
            (_list ??= []).Add(_acc);
            _acc = f;
        }

        public BigInt Result()
        {
            if (_acc != 0)
            {
                (_list ??= []).Add(_acc);
                _acc = 0;
            }
            if (_list is null || _list.Count == 0)
            {
                return One;
            }
            ulong[] w = [.. _list];
            return ProductOfWords(w, 0, w.Length);
        }
    }

    /// <summary>Balanced product of words[start..start+count).</summary>
    private static BigInt ProductOfWords(ulong[] words, int start, int count)
    {
        if (count <= 16)
        {
            // Schoolbook accumulation with in-place MulScalar into one buffer.
            ulong[] d = new ulong[count + 1];
            d[0] = words[start];
            int n = 1;
            for (int i = 1; i < count; i++)
            {
                ulong c = Nat.MulScalar(d, d, n, words[start + i]);
                if (c != 0)
                {
                    d[n++] = c;
                }
            }
            return FromMagnitude(d, n, 1);
        }
        int half = count / 2;
        return ProductOfWords(words, start, half) * ProductOfWords(words, start + half, count - half);
    }

    /// <summary>Binomial coefficient C(n, k) (mpz_bin_uiui).</summary>
    public static BigInt Binomial(ulong n, ulong k)
    {
        if (k > n)
        {
            return default;
        }
        if (k > n - k)
        {
            k = n - k;
        }
        if (k == 0)
        {
            return One;
        }
        if (k == 1)
        {
            return new BigInt(n);
        }
        if (k <= 24 || n > PrimeMethodLimit)
        {
            // Few factors (or a sieve too large): exact product / k!.
            BigInt num = ProductRange(n - k + 1, n);
            return k <= 20 ? DivRem(num, SmallFactorial(k), out _) : DivExact(num, Factorial(k));
        }

        // Kummer: the exponent of p in C(n, k) is the number of carries when adding k and
        // n-k in base p, i.e. Σ_i floor(n/p^i) - floor(k/p^i) - floor((n-k)/p^i).
        uint[] primes = SmallPrimes.PrimesUpTo((uint)n);
        ulong nk = n - k;
        var words = new WordProduct();
        ulong sqrtN = SmallPrimes.Sqrt64(n);
        for (int i = 1; i < primes.Length; i++)
        {
            ulong p = primes[i];
            if (p > nk)
            {
                // Primes in (n-k, n] divide the numerator exactly once and not k!.
                words.Add(p);
            }
            else if (p <= sqrtN)
            {
                ulong pe = 1, a = n, b = k, c = nk;
                while (a >= p)
                {
                    a /= p; b /= p; c /= p;
                    if (a - b - c != 0)
                    {
                        pe *= p;
                    }
                }
                words.Add(pe);
            }
            else if (n / p - k / p - nk / p != 0)
            {
                words.Add(p);
            }
        }
        long twos = BitOperations.PopCount(k) + BitOperations.PopCount(nk) - BitOperations.PopCount(n);
        return words.Result().ShiftLeft(twos);
    }

    /// <summary>Product of lo..hi inclusive.</summary>
    private static BigInt ProductRange(ulong lo, ulong hi)
    {
        Debug.Assert(lo <= hi);
        ulong count = hi - lo + 1;
        if (count <= 32)
        {
            var words = new WordProduct();
            for (ulong v = lo; v <= hi; v++)
            {
                words.Add(v);
            }
            return words.Result();
        }
        ulong mid = lo + count / 2;
        return ProductRange(lo, mid - 1) * ProductRange(mid, hi);
    }

    // ---------------------------------------------------------------
    // Fibonacci / Lucas numbers (mpz_fib_ui etc.)
    // ---------------------------------------------------------------

    /// <summary>F(n) (mpz_fib_ui).</summary>
    public static BigInt Fibonacci(ulong n)
    {
        if (n < (ulong)FibTable.Length)
        {
            return new BigInt(FibTable[n]);
        }
        // Last doubling step with one multiplication instead of two squarings:
        // F(2k) = F(k) (F(k) + 2F(k-1)), F(2k+1) = (2F(k) + F(k-1)) (2F(k) - F(k-1)) + 2(-1)^k.
        ulong k = n >> 1;
        FibonacciPairCore(k, out BigInt fk, out BigInt fk1);
        if ((n & 1) == 0)
        {
            return fk * (fk + (fk1 << 1));
        }
        BigInt twoFk = fk << 1;
        return (twoFk + fk1) * (twoFk - fk1) + ((k & 1) == 0 ? Two : -Two);
    }

    /// <summary>F(n) and F(n-1) (mpz_fib2_ui).</summary>
    public static (BigInt Current, BigInt Previous) FibonacciPair(ulong n)
    {
        if (n == 0)
        {
            return (Zero, One); // F(-1) = 1
        }
        FibonacciPairCore(n, out BigInt fn, out BigInt fnm1);
        return (fn, fnm1);
    }

    /// <summary>F(0..93) (F(93) is the largest that fits in 64 bits).</summary>
    private static readonly ulong[] FibTable = BuildFibTable();

    private static ulong[] BuildFibTable()
    {
        ulong[] t = new ulong[94];
        t[1] = 1;
        for (int i = 2; i < t.Length; i++)
        {
            t[i] = t[i - 1] + t[i - 2];
        }
        return t;
    }

    /// <summary>
    /// (F(n), F(n-1)) for n &gt;= 1 by doubling with two squarings per step:
    /// F(2k-1) = F(k)^2 + F(k-1)^2, F(2k+1) = 4F(k)^2 - F(k-1)^2 + 2(-1)^k, F(2k) = F(2k+1) - F(2k-1).
    /// </summary>
    private static void FibonacciPairCore(ulong n, out BigInt fn, out BigInt fnm1)
    {
        Debug.Assert(n >= 1);
        if (n < (ulong)FibTable.Length)
        {
            fn = new BigInt(FibTable[n]);
            fnm1 = new BigInt(FibTable[n - 1]);
            return;
        }
        // Iterative doubling over pooled limb buffers: no temporaries per step.
        int sh = 63 - BitOperations.LeadingZeroCount(n) - 5; // top 6 bits: k < 64 is in the table
        ulong k = n >> sh;
        // F(n) < phi^n with log2(phi) < 0.6943; +2 limbs of slack for the 4*F(k)^2 intermediate.
        int cap = (int)(n * 0.6943 / 64) + 3;
        ulong[] pool = System.Buffers.ArrayPool<ulong>.Shared.Rent(8 * cap + 8);
        Span<ulong> f = pool.AsSpan(0, cap), g = pool.AsSpan(cap, cap);
        Span<ulong> s1 = pool.AsSpan(2 * cap, 2 * cap + 2), s2 = pool.AsSpan(4 * cap + 2, 2 * cap + 2);
        Span<ulong> t = pool.AsSpan(6 * cap + 4, cap + 2);
        f[0] = FibTable[k];
        g[0] = FibTable[k - 1];
        int fnl = 1, gnl = g[0] == 0 ? 0 : 1;
        for (int i = sh - 1; i >= 0; i--)
        {
            // s1 = F(k)^2, s2 = F(k-1)^2
            Nat.Square(s1, f, fnl);
            int s1n = Nat.TrimmedLength(s1, 2 * fnl);
            int s2n = 0;
            if (gnl > 0)
            {
                Nat.Square(s2, g, gnl);
                s2n = Nat.TrimmedLength(s2, 2 * gnl);
            }
            // g = F(2k-1) = s1 + s2
            ulong c = s2n > 0 ? Nat.Add(g, s1, s1n, s2, s2n) : CopyLimbs(g, s1, s1n);
            g[s1n] = c;
            int an = Nat.TrimmedLength(g, s1n + 1);
            // f = F(2k+1) = 4 s1 - s2 + 2(-1)^k
            f[s1n] = Nat.ShiftLeft(f, s1, s1n, 2);
            int bn = s1n + 1;
            if (s2n > 0)
            {
                Nat.SubInPlace(f, bn, s2, s2n);
            }
            if ((k & 1) == 0)
            {
                Nat.AddScalarInPlace(f, bn, 2);
            }
            else
            {
                Nat.SubScalarInPlace(f, bn, 2);
            }
            bn = Nat.TrimmedLength(f, bn);
            // t = F(2k) = F(2k+1) - F(2k-1)
            Nat.Sub(t, f, bn, g, an);
            int cn = Nat.TrimmedLength(t, bn);
            if (((n >> i) & 1) != 0)
            {
                // (F(2k+1), F(2k)): f already holds F(2k+1).
                t.Slice(0, cn).CopyTo(g);
                fnl = bn;
                gnl = cn;
                k = 2 * k + 1;
            }
            else
            {
                // (F(2k), F(2k-1)): g already holds F(2k-1).
                t.Slice(0, cn).CopyTo(f);
                fnl = cn;
                gnl = an;
                k = 2 * k;
            }
        }
        ulong[] fd = NewLimbs(fnl), gd = NewLimbs(Math.Max(gnl, 1));
        f.Slice(0, fnl).CopyTo(fd);
        g.Slice(0, gnl).CopyTo(gd);
        System.Buffers.ArrayPool<ulong>.Shared.Return(pool);
        fn = FromMagnitude(fd, fnl, 1);
        fnm1 = FromMagnitude(gd, gnl, 1);

        static ulong CopyLimbs(Span<ulong> dst, ReadOnlySpan<ulong> src, int len)
        {
            src.Slice(0, len).CopyTo(dst);
            return 0;
        }
    }

    /// <summary>Lucas number L(n) (mpz_lucnum_ui).</summary>
    public static BigInt Lucas(ulong n)
    {
        if (n == 0)
        {
            return Two;
        }
        // L(n) = F(n) + 2 F(n-1)
        FibonacciPairCore(n, out BigInt fn, out BigInt fnm1);
        return fn + (fnm1 << 1);
    }

    /// <summary>L(n) and L(n-1) (mpz_lucnum2_ui).</summary>
    public static (BigInt Current, BigInt Previous) LucasPair(ulong n)
    {
        if (n == 0)
        {
            return (Two, MinusOne); // L(-1) = -1
        }
        FibonacciPairCore(n, out BigInt fn, out BigInt fnm1);
        // L(n) = F(n) + 2 F(n-1); L(n-1) = 2 F(n) - F(n-1)
        return (fn + (fnm1 << 1), (fn << 1) - fnm1);
    }
}
