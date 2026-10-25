namespace Sdcb.Arithmetic.Core;

/// <summary>Small prime utilities shared by primality testing, NextPrime and Primorial.</summary>
internal static class SmallPrimes
{
    /// <summary>Primes below 1024, used for trial division.</summary>
    public static readonly uint[] TrialPrimes = Sieve(1024);

    /// <summary>Simple sieve of Eratosthenes returning all primes below limit.</summary>
    public static uint[] Sieve(int limit)
    {
        if (limit < 3)
        {
            return [];
        }
        bool[] composite = new bool[limit];
        var primes = new List<uint>();
        for (int i = 2; i < limit; i++)
        {
            if (!composite[i])
            {
                primes.Add((uint)i);
                for (long j = (long)i * i; j < limit; j += i)
                {
                    composite[j] = true;
                }
            }
        }
        return [.. primes];
    }

    /// <summary>All primes &lt;= n via an odd-only bit sieve (memory n/16 bytes).</summary>
    public static uint[] PrimesUpTo(uint n)
    {
        if (n < 2)
        {
            return [];
        }
        // Bit i represents the odd number 2i+1.
        uint half = n / 2 + 1;
        ulong[] composite = new ulong[(half + 63) / 64];
        for (uint i = 1; ; i++)
        {
            ulong p = 2UL * i + 1;
            if (p * p > n)
            {
                break;
            }
            if ((composite[i >> 6] & (1UL << (int)(i & 63))) != 0)
            {
                continue;
            }
            for (ulong j = p * p / 2; j < half; j += p)
            {
                composite[j >> 6] |= 1UL << (int)(j & 63);
            }
        }
        int count = 1;
        for (uint i = 1; i < half && 2 * i + 1 <= n; i++)
        {
            if ((composite[i >> 6] & (1UL << (int)(i & 63))) == 0)
            {
                count++;
            }
        }
        uint[] primes = new uint[count];
        primes[0] = 2;
        int k = 1;
        for (uint i = 1; i < half && 2 * i + 1 <= n; i++)
        {
            if ((composite[i >> 6] & (1UL << (int)(i & 63))) == 0)
            {
                primes[k++] = 2 * i + 1;
            }
        }
        return primes;
    }

    /// <summary>Binary GCD on 64-bit values.</summary>
    public static ulong Gcd64(ulong a, ulong b)
    {
        if (a == 0)
        {
            return b;
        }
        if (b == 0)
        {
            return a;
        }
        int shift = System.Numerics.BitOperations.TrailingZeroCount(a | b);
        a >>= System.Numerics.BitOperations.TrailingZeroCount(a);
        do
        {
            // Branch-free step (cmov): (a, b) <- (min, |a - b|), then strip b's twos.
            b >>= System.Numerics.BitOperations.TrailingZeroCount(b);
            ulong d = b - a;
            ulong m = a < b ? a : b;
            b = a < b ? d : a - b;
            a = m;
        }
        while (b != 0);
        return a << shift;
    }

    /// <summary>Floor integer square root of a 64-bit value.</summary>
    public static ulong Sqrt64(ulong n)
    {
        if (n == 0)
        {
            return 0;
        }
        ulong r = (ulong)Math.Sqrt(n);
        if (r > 0xFFFFFFFF)
        {
            r = 0xFFFFFFFF;
        }
        while (r * r > n)
        {
            r--;
        }
        while (r < 0xFFFFFFFF && (r + 1) * (r + 1) <= n)
        {
            r++;
        }
        return r;
    }
}
