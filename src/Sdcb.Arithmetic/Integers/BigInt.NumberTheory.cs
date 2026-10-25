using System.Diagnostics;
using System.Numerics;
using Sdcb.Arithmetic.Core;

namespace Sdcb.Arithmetic;

public readonly partial struct BigInt
{
    // ---------------------------------------------------------------
    // Modular exponentiation (mpz_powm)
    // ---------------------------------------------------------------

    /// <summary>
    /// Computes value^exponent mod |modulus| (result in [0, |modulus|)).
    /// Negative exponents require the base to be invertible (like mpz_powm).
    /// </summary>
    public static BigInt ModPow(BigInt value, BigInt exponent, BigInt modulus)
    {
        BigInt m = modulus.Abs();
        if (m._size == 0)
        {
            throw new DivideByZeroException("Modulus must be nonzero.");
        }
        if (m.IsOne)
        {
            return Zero;
        }
        if (exponent._size < 0)
        {
            value = ModInverse(value, m);
            exponent = -exponent;
        }
        if (exponent._size == 0)
        {
            return One;
        }
        BigInt b = Mod(value, m);
        if (b._size == 0)
        {
            return Zero;
        }
        if (m.IsOdd)
        {
            return PowModMontgomery(b, exponent, m);
        }
        return PowModGeneric(b, exponent, m);
    }

    private static BigInt PowModMontgomery(BigInt b, BigInt e, BigInt m)
    {
        var ctx = new MontgomeryContext(m.Limbs);
        int n = ctx.N;
        ulong[] t = new ulong[2 * n + 1];
        ulong[] xm = new ulong[n];
        ulong[] result = new ulong[n];

        Span<ulong> bPad = n <= 128 ? stackalloc ulong[n] : new ulong[n];
        bPad.Clear();
        b.Limbs.CopyTo(bPad);
        ctx.ToMont(xm, bPad, t);

        long bits = e.GetBitLength();
        int window = bits switch
        {
            <= 8 => 1,
            <= 64 => 3,
            <= 512 => 4,
            <= 2048 => 5,
            _ => 6,
        };

        // Precompute odd powers xm^1, xm^3, ..., xm^(2^window - 1).
        int tableSize = 1 << (window - 1);
        ulong[][] table = new ulong[tableSize][];
        table[0] = (ulong[])xm.Clone();
        if (tableSize > 1)
        {
            ulong[] sq = new ulong[n];
            ctx.MontSqr(sq, xm, t);
            for (int i = 1; i < tableSize; i++)
            {
                table[i] = new ulong[n];
                ctx.MontMul(table[i], table[i - 1], sq, t);
            }
        }

        ctx.RMod.CopyTo(result); // Montgomery form of 1
        long i2 = bits - 1;
        bool started = false;
        while (i2 >= 0)
        {
            if (!e.TestBit(i2))
            {
                if (started)
                {
                    ctx.MontSqr(result, result, t);
                }
                i2--;
                continue;
            }
            // Gather a window ending on a set bit.
            int len = (int)Math.Min(window, i2 + 1);
            long j = i2 - len + 1;
            while (!e.TestBit(j))
            {
                j++;
                len--;
            }
            int val = 0;
            for (long k = i2; k >= j; k--)
            {
                val = (val << 1) | (e.TestBit(k) ? 1 : 0);
            }
            if (started)
            {
                for (int k = 0; k < len; k++)
                {
                    ctx.MontSqr(result, result, t);
                }
                ctx.MontMul(result, result, table[val >> 1], t);
            }
            else
            {
                table[val >> 1].CopyTo(result, 0);
                started = true;
            }
            i2 = j - 1;
        }

        ctx.FromMont(result, result, t);
        return FromMagnitude(result, n, 1);
    }

    private static BigInt PowModGeneric(BigInt b, BigInt e, BigInt m)
    {
        BigInt result = One;
        BigInt acc = b;
        long bits = e.GetBitLength();
        for (long i = 0; i < bits; i++)
        {
            if (e.TestBit(i))
            {
                result = Mod(result * acc, m);
            }
            if (i + 1 < bits)
            {
                acc = Mod(acc.Square(), m);
            }
        }
        return result;
    }

    /// <summary>value^exponent mod modulus for a ulong exponent (mpz_powm_ui).</summary>
    public static BigInt ModPow(BigInt value, ulong exponent, BigInt modulus) => ModPow(value, new BigInt(exponent), modulus);

    // ---------------------------------------------------------------
    // Primality (BPSW + extra Miller-Rabin, mpz_probab_prime_p)
    // ---------------------------------------------------------------

    /// <summary>
    /// Probabilistic primality test on |value|. Runs trial division, a base-2
    /// Miller-Rabin test, a strong Lucas test (BPSW), and additional
    /// Miller-Rabin rounds (<paramref name="reps"/>). Values below 2^64 are decided
    /// deterministically (<see cref="Primality.Prime"/>).
    /// </summary>
    public Primality GetPrimality(int reps = 25)
    {
        BigInt n = Abs();
        if (n.LimbCount <= 1)
        {
            ulong v = n._size == 0 ? 0 : n._d![0];
            return IsPrime64(v) ? Primality.Prime : Primality.Composite;
        }

        if (n.IsEven)
        {
            return Primality.Composite;
        }
        // Trial division by small primes.
        foreach (uint p in SmallPrimes.TrialPrimes)
        {
            if (Nat.RemScalar(n.Limbs, n.LimbCount, p) == 0)
            {
                return Primality.Composite;
            }
        }
        // Miller-Rabin base 2.
        BigInt nm1 = n - One;
        long s = nm1.TrailingZeroCount();
        BigInt d = nm1.ShiftRightTruncate(s);
        if (!MillerRabinRound(n, nm1, d, s, Two))
        {
            return Primality.Composite;
        }
        // Strong Lucas test (completes the BPSW test).
        if (!StrongLucasProbablePrime(n))
        {
            return Primality.Composite;
        }
        // Extra Miller-Rabin rounds with pseudo-random bases.
        ulong seed = 0x9E3779B97F4A7C15UL ^ n._d![0];
        for (int i = 0; i < reps; i++)
        {
            seed = seed * 6364136223846793005UL + 1442695040888963407UL;
            BigInt a = Mod(new BigInt(seed | 2), n - Two) + Two;
            if (!MillerRabinRound(n, nm1, d, s, a))
            {
                return Primality.Composite;
            }
        }
        return Primality.ProbablePrime;
    }

    /// <summary>True when |value| is prime or probably prime; see <see cref="GetPrimality"/>.</summary>
    public bool IsProbablePrime(int reps = 25) => GetPrimality(reps) != Primality.Composite;

    private static bool MillerRabinRound(BigInt n, BigInt nm1, BigInt d, long s, BigInt a)
    {
        BigInt x = ModPow(a, d, n);
        if (x.IsOne || x == nm1)
        {
            return true;
        }
        for (long i = 1; i < s; i++)
        {
            x = Mod(x.Square(), n);
            if (x == nm1)
            {
                return true;
            }
            if (x.IsOne)
            {
                return false;
            }
        }
        return false;
    }

    /// <summary>Strong Lucas probable prime test with Selfridge's parameters.</summary>
    private static bool StrongLucasProbablePrime(BigInt n)
    {
        // Find D in 5, -7, 9, -11, ... with Jacobi(D, n) == -1.
        long dCandidate = 5;
        while (true)
        {
            BigInt dBig = new(dCandidate);
            int j = Jacobi(Mod(dBig, n), n);
            if (j == -1)
            {
                break;
            }
            if (j == 0 && CompareAbs(dBig, n) < 0)
            {
                return false; // shares a factor
            }
            if (dCandidate == 13 && IsPerfectSquare(n))
            {
                return false;
            }
            dCandidate = dCandidate > 0 ? -(dCandidate + 2) : -(dCandidate - 2);
        }
        BigInt D = Mod(new BigInt(dCandidate), n);
        // P = 1, Q = (1 - D) / 4
        BigInt q = Mod(FloorDivide(One - new BigInt(dCandidate), new BigInt(4)), n);

        BigInt nPlus1 = n + One;
        long s = nPlus1.TrailingZeroCount();
        BigInt dOdd = nPlus1.ShiftRightTruncate(s);

        // Compute U_dOdd, V_dOdd via binary chain (Joye-Quisquater style doubling).
        BigInt u = One, v = One, qk = q; // U_1 = 1, V_1 = P = 1
        long bits = dOdd.GetBitLength();
        for (long i = bits - 2; i >= 0; i--)
        {
            // Double: U_2k = U_k * V_k ; V_2k = V_k^2 - 2*Q^k
            u = Mod(u * v, n);
            v = Mod(v.Square() - (qk << 1), n);
            qk = Mod(qk.Square(), n);
            if (dOdd.TestBit(i))
            {
                // Increment: U_{2k+1} = (P*U + V)/2 ; V_{2k+1} = (D*U + P*V)/2
                BigInt u2 = u + v;
                BigInt v2 = D * u + v;
                if (u2.IsOdd)
                {
                    u2 += n;
                }
                if (v2.IsOdd)
                {
                    v2 += n;
                }
                u = Mod(u2.ShiftRightTruncate(1), n);
                v = Mod(v2.ShiftRightTruncate(1), n);
                qk = Mod(qk * q, n);
            }
        }

        if (u._size == 0 || v._size == 0)
        {
            return true;
        }
        for (long r = 1; r < s; r++)
        {
            v = Mod(v.Square() - (qk << 1), n);
            if (v._size == 0)
            {
                return true;
            }
            qk = Mod(qk.Square(), n);
        }
        return false;
    }

    /// <summary>Deterministic Miller-Rabin for 64-bit values.</summary>
    internal static bool IsPrime64(ulong n)
    {
        if (n < 2)
        {
            return false;
        }
        foreach (ulong p in (ReadOnlySpan<ulong>)[2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37])
        {
            if (n % p == 0)
            {
                return n == p;
            }
        }
        ulong d = n - 1;
        int s = BitOperations.TrailingZeroCount(d);
        d >>= s;
        foreach (ulong a in (ReadOnlySpan<ulong>)[2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37])
        {
            ulong x = PowMod64(a % n, d, n);
            if (x == 1 || x == n - 1)
            {
                continue;
            }
            bool composite = true;
            for (int i = 1; i < s; i++)
            {
                x = MulMod64(x, x, n);
                if (x == n - 1)
                {
                    composite = false;
                    break;
                }
            }
            if (composite)
            {
                return false;
            }
        }
        return true;
    }

    internal static ulong MulMod64(ulong a, ulong b, ulong m)
    {
        ulong hi = Nat.BigMul(a, b, out ulong lo);
        if (hi == 0)
        {
            return lo % m;
        }
        Nat.DivRem128(hi % m, lo, m, out ulong r);
        return r;
    }

    internal static ulong PowMod64(ulong b, ulong e, ulong m)
    {
        ulong result = 1;
        b %= m;
        while (e != 0)
        {
            if ((e & 1) != 0)
            {
                result = MulMod64(result, b, m);
            }
            b = MulMod64(b, b, m);
            e >>= 1;
        }
        return result;
    }

    /// <summary>Smallest probable prime strictly greater than this value (mpz_nextprime).</summary>
    public BigInt NextPrime()
    {
        if (CompareTo(Two) < 0)
        {
            return Two;
        }
        BigInt candidate = this + One;
        if (candidate.IsEven)
        {
            candidate += One;
        }
        return candidate.CompareTo(3L) <= 0 ? candidate : SieveSearch(candidate, step: 2);
    }

    /// <summary>Removes all occurrences of a factor (mpz_remove); returns the stripped value.</summary>
    public static BigInt RemoveFactor(BigInt value, BigInt factor, out ulong count)
    {
        if (factor.Abs().CompareTo(One) <= 0)
        {
            throw new ArgumentOutOfRangeException(nameof(factor), "Factor must have |factor| > 1.");
        }
        count = 0;
        if (value._size == 0)
        {
            return value;
        }
        // Fast path for factor == 2.
        if (factor.LimbCount == 1 && factor._d![0] == 2)
        {
            long tz = value.TrailingZeroCount();
            count = (ulong)tz;
            return value.ShiftRightTruncate(tz);
        }
        while (true)
        {
            DivRem(value, factor, out BigInt q, out BigInt r);
            if (!r.IsZero)
            {
                return value;
            }
            value = q;
            count++;
        }
    }

    /// <summary>
    /// Number of digits needed in the given base, like mpz_sizeinbase: exact for
    /// powers of two, otherwise exact or one too big.
    /// </summary>
    public long CountDigits(int inBase)
    {
        if (inBase < 2 || inBase > 62)
        {
            throw new ArgumentOutOfRangeException(nameof(inBase));
        }
        if (_size == 0)
        {
            return 1;
        }
        long bits = GetBitLength();
        if ((inBase & (inBase - 1)) == 0)
        {
            int bpd = System.Numerics.BitOperations.TrailingZeroCount((uint)inBase);
            return (bits + bpd - 1) / bpd;
        }
        return (long)(bits / Math.Log2(inBase)) + 1;
    }

    // ---------------------------------------------------------------
    // Roots (mpz_sqrt / mpz_sqrtrem / mpz_root / mpz_rootrem)
    // ---------------------------------------------------------------

    /// <summary>Floor square root (mpz_sqrt); value must be non-negative.</summary>
    public static BigInt Sqrt(BigInt value) => SqrtCore(value, wantRemainder: false, out _);

    /// <summary>Floor square root with remainder: value = root^2 + rem (mpz_sqrtrem).</summary>
    public static BigInt Sqrt(BigInt value, out BigInt remainder) => SqrtCore(value, wantRemainder: true, out remainder);

    private static BigInt SqrtCore(BigInt value, bool wantRemainder, out BigInt remainder)
    {
        if (value._size < 0)
        {
            throw new ArithmeticException("Square root of a negative value.");
        }
        remainder = default;
        if (value._size == 0)
        {
            return default;
        }
        int nn = value.LimbCount, tn = (nn + 1) / 2;
        ulong[] s = NewLimbs(tn);
        if (!wantRemainder)
        {
            Nat.SqrtRem(s, default, value.Limbs, nn);
            return new BigInt(s, tn);
        }
        ulong[] r = new ulong[tn + 1];
        int rn = Nat.SqrtRem(s, r, value.Limbs, nn);
        remainder = rn == 0 ? default : new BigInt(r, rn);
        return new BigInt(s, tn);
    }

    /// <summary>True when the value is a perfect square (mpz_perfect_square_p).</summary>
    public static bool IsPerfectSquare(BigInt value)
    {
        if (value._size < 0)
        {
            return false;
        }
        if (value._size == 0)
        {
            return true;
        }
        // Quick residue filter mod 64.
        ulong low = value._d![0];
        // Squares mod 64 ∈ {0,1,4,9,16,17,25,33,36,41,49,57}
        const ulong mask64 = (1UL << 0) | (1UL << 1) | (1UL << 4) | (1UL << 9) | (1UL << 16) | (1UL << 17)
            | (1UL << 25) | (1UL << 33) | (1UL << 36) | (1UL << 41) | (1UL << 49) | (1UL << 57);
        if (((mask64 >> (int)(low & 63)) & 1) == 0)
        {
            return false;
        }
        Sqrt(value, out BigInt rem);
        return rem._size == 0;
    }

    /// <summary>Floor k-th root (mpz_root). Negative values allowed for odd k.</summary>
    public static BigInt Root(BigInt value, uint k) => Root(value, k, out _);

    /// <summary>Floor k-th root with remainder: value = root^k + rem (mpz_rootrem).</summary>
    public static BigInt Root(BigInt value, uint k, out BigInt remainder)
    {
        if (k == 0)
        {
            throw new ArgumentOutOfRangeException(nameof(k));
        }
        if (value._size < 0)
        {
            if ((k & 1) == 0)
            {
                throw new ArithmeticException("Even root of a negative value.");
            }
            BigInt rootNeg = Root(-value, k, out BigInt remAbs);
            remainder = -remAbs;
            return -rootNeg;
        }
        if (k == 1 || value._size == 0 || value.IsOne)
        {
            remainder = default;
            return value;
        }
        if (k == 2)
        {
            return Sqrt(value, out remainder);
        }

        long bits = value.GetBitLength();
        if (bits <= k)
        {
            // root is 1
            remainder = value - One;
            return One;
        }
        // Initial upper bound: 2^ceil(bits/k)
        BigInt x = One.ShiftLeft((bits + k - 1) / k);
        BigInt km1 = new(k - 1);
        while (true)
        {
            BigInt xkm1 = Pow(x, k - 1);
            BigInt next = (km1 * x + value / xkm1) / new BigInt(k);
            if (next.CompareTo(x) >= 0)
            {
                break;
            }
            x = next;
        }
        // x is now floor(root) or slightly above; fix up.
        while (Pow(x, k).CompareTo(value) > 0)
        {
            x -= One;
        }
        remainder = value - Pow(x, k);
        Debug.Assert(remainder._size >= 0);
        return x;
    }

    /// <summary>True when value == a^k for some a and k &gt; 1 (mpz_perfect_power_p).</summary>
    public static bool IsPerfectPower(BigInt value)
    {
        if (value._size == 0 || value.IsOne || value.IsMinusOne)
        {
            return true;
        }
        BigInt abs = value.Abs();
        long bits = abs.GetBitLength();
        foreach (uint p in SmallPrimes.TrialPrimes)
        {
            if (p > bits)
            {
                break;
            }
            if (value._size < 0 && p == 2)
            {
                continue; // negative values need odd exponents
            }
            Root(abs, p, out BigInt rem);
            if (rem._size == 0)
            {
                return true;
            }
        }
        return false;
    }

    // ---------------------------------------------------------------
    // Factorials & combinatorics
    // ---------------------------------------------------------------

    /// <summary>Double factorial n!! (mpz_2fac_ui).</summary>
    public static BigInt DoubleFactorial(ulong n)
    {
        if (n < 2)
        {
            return One;
        }
        return ProductStride(n < 2 ? 2 : (2 - (n & 1)), n, 2);
    }

    /// <summary>Multifactorial n!^(m) (mpz_mfac_uiui).</summary>
    public static BigInt MultiFactorial(ulong n, ulong m)
    {
        ArgumentOutOfRangeException.ThrowIfZero(m);
        if (n < 2)
        {
            return One;
        }
        ulong first = n % m;
        if (first == 0)
        {
            first = m;
        }
        while (first < 2)
        {
            first += m;
        }
        return first > n ? One : ProductStride(first, n, m);
    }

    /// <summary>Product of all primes &lt;= n (mpz_primorial_ui).</summary>
    public static BigInt Primorial(ulong n)
    {
        if (n < 2)
        {
            return One;
        }
        if (n > int.MaxValue - 2)
        {
            throw new ArgumentOutOfRangeException(nameof(n));
        }
        uint[] primes = SmallPrimes.Sieve((int)n + 1);
        return ProductOf(primes, 0, primes.Length);
    }

    private static BigInt ProductOf(uint[] values, int start, int count)
    {
        if (count == 0)
        {
            return One;
        }
        if (count <= 8)
        {
            BigInt acc = new((ulong)values[start]);
            for (int i = 1; i < count; i++)
            {
                acc = MulSigned(acc, values[start + i]);
            }
            return acc;
        }
        int half = count / 2;
        return ProductOf(values, start, half) * ProductOf(values, start + half, count - half);
    }

    /// <summary>Product of first, first+stride, ..., up to hi inclusive.</summary>
    private static BigInt ProductStride(ulong first, ulong hi, ulong stride)
    {
        ulong count = (hi - first) / stride + 1;
        if (count <= 16)
        {
            BigInt acc = One;
            ulong chunk = 1;
            for (ulong i = 0; i < count; i++)
            {
                ulong v = first + i * stride;
                if (chunk > ulong.MaxValue / v)
                {
                    acc = MulSignedU(acc, chunk);
                    chunk = 1;
                }
                chunk *= v;
            }
            return MulSignedU(acc, chunk);
        }
        ulong halfCount = count / 2;
        return ProductStride(first, first + (halfCount - 1) * stride, stride)
             * ProductStride(first + halfCount * stride, hi, stride);
    }

    private static BigInt MulSignedU(BigInt v, ulong f)
    {
        if (f == 0 || v._size == 0)
        {
            return default;
        }
        if (f == 1)
        {
            return v;
        }
        int n = v.LimbCount;
        ulong[] d = new ulong[n + 1];
        d[n] = Nat.MulScalar(d, v.Limbs, n, f);
        return FromMagnitude(d, n + 1, v._size < 0 ? -1 : 1);
    }

    /// <summary>Binomial coefficient C(n, k) for BigInt n (mpz_bin_ui). Supports negative n via C(n,k) = (-1)^k C(k-n-1, k).</summary>
    public static BigInt Binomial(BigInt n, ulong k)
    {
        if (n._size < 0)
        {
            BigInt r = Binomial(k - 1 - n, k);
            return (k & 1) == 0 ? r : -r;
        }
        if (n.LimbCount <= 1)
        {
            ulong nv = n._size == 0 ? 0 : n._d![0];
            return Binomial(nv, k);
        }
        if (k == 0)
        {
            return One;
        }
        return DivExact(ProductConsecutive(n - new BigInt(k) + One, k), Factorial(k));
    }

    /// <summary>start * (start+1) * ... * (start+count-1) by a balanced product tree.</summary>
    private static BigInt ProductConsecutive(BigInt start, ulong count)
    {
        if (count <= 8)
        {
            BigInt acc = start;
            for (ulong i = 1; i < count; i++)
            {
                acc *= start + new BigInt(i);
            }
            return acc;
        }
        ulong half = count / 2;
        return ProductConsecutive(start, half) * ProductConsecutive(start + new BigInt(half), count - half);
    }
}
