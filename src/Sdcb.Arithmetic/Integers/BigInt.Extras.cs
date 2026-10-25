using System.Numerics;
using Sdcb.Arithmetic.Core;

namespace Sdcb.Arithmetic;

public readonly partial struct BigInt
{
    // ---------------------------------------------------------------
    // Powers-of-two quotients and remainders (mpz_{t,f,c}div_{q,r}_2exp)
    // ---------------------------------------------------------------

    /// <summary>Ceiling division by 2^shift (mpz_cdiv_q_2exp).</summary>
    internal BigInt ShiftRightCeiling(long shift)
    {
        BigInt t = ShiftRightTruncate(shift);
        return _size > 0 && shift > 0 && AnyBitsBelow(shift) ? t + One : t;
    }

    /// <summary>Ieee754Remainder of dividing by 2^bits with the quotient rounded as in <see cref="ShiftRight(long, IntegerRounding)"/> (mpz_{t,f,c}div_r_2exp).</summary>
    public BigInt RemainderPow2(long bits, IntegerRounding rounding = IntegerRounding.Truncate) => rounding switch
    {
        IntegerRounding.Truncate => RemainderPow2Truncate(bits),
        IntegerRounding.Floor => ModPow2(bits),
        IntegerRounding.Ceiling => RemainderPow2Ceiling(bits),
        _ => throw new ArgumentOutOfRangeException(nameof(rounding)),
    };

    /// <summary>Ieee754Remainder of truncating division by 2^bits, with the sign of the value (mpz_tdiv_r_2exp).</summary>
    internal BigInt RemainderPow2Truncate(long bits)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(bits);
        BigInt low = LowBitsMagnitude(bits);
        return _size < 0 ? -low : low;
    }

    /// <summary>Non-negative remainder of floor division by 2^bits, i.e. value mod 2^bits (mpz_fdiv_r_2exp).</summary>
    public BigInt ModPow2(long bits)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(bits);
        BigInt low = LowBitsMagnitude(bits);
        return _size < 0 && !low.IsZero ? (One << (int)bits) - low : low;
    }

    /// <summary>Non-positive remainder of ceiling division by 2^bits (mpz_cdiv_r_2exp).</summary>
    internal BigInt RemainderPow2Ceiling(long bits)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(bits);
        BigInt low = LowBitsMagnitude(bits);
        return _size > 0 && !low.IsZero ? low - (One << (int)bits) : -low;
    }

    /// <summary>|value| mod 2^bits.</summary>
    private BigInt LowBitsMagnitude(long bits)
    {
        int n = LimbCount;
        if (n == 0 || bits == 0)
        {
            return default;
        }
        if (bits >= (long)n * 64)
        {
            return Abs();
        }
        int full = (int)(bits >> 6), rem = (int)(bits & 63);
        int len = full + (rem != 0 ? 1 : 0);
        ulong[] d = NewLimbs(len);
        Limbs.Slice(0, len).CopyTo(d);
        if (rem != 0)
        {
            d[len - 1] &= (1UL << rem) - 1;
        }
        return FromMagnitude(d, len, 1);
    }

    /// <summary>True when 2^bits divides the value (mpz_divisible_2exp_p).</summary>
    public bool IsDivisibleByPow2(long bits) => _size == 0 || TrailingZeroCount() >= bits;

    /// <summary>True when a ≡ c (mod d); d = 0 means a == c (mpz_congruent_p).</summary>
    public static bool IsCongruent(BigInt a, BigInt c, BigInt d)
        => d.IsZero ? a == c : IsDivisible(a - c, d);

    /// <summary>True when a ≡ c (mod 2^bits) (mpz_congruent_2exp_p).</summary>
    public static bool IsCongruentModPow2(BigInt a, BigInt c, long bits) => (a - c).IsDivisibleByPow2(bits);

    // ---------------------------------------------------------------
    // Comparisons and conversions with double / machine integers
    // ---------------------------------------------------------------

    /// <summary>Compares with a double exactly (mpz_cmp_d); NaN throws, ±infinity compare as expected.</summary>
    public int CompareTo(double other)
    {
        if (double.IsNaN(other))
        {
            throw new ArgumentException("Cannot compare with NaN.", nameof(other));
        }
        if (double.IsInfinity(other))
        {
            return other > 0 ? -1 : 1;
        }
        BigInt t = new(other); // truncated toward zero
        int c = CompareTo(t);
        if (c != 0)
        {
            return c;
        }
        double frac = other - Math.Truncate(other);
        return frac > 0 ? -1 : frac < 0 ? 1 : 0;
    }

    /// <summary>Compares |value| with |d| (mpz_cmpabs_d).</summary>
    public int CompareAbs(double d) => Abs().CompareTo(Math.Abs(d));

    /// <summary>
    /// Returns m with 0.5 &lt;= |m| &lt; 1 truncated toward zero and sets exponent so that
    /// value ≈ m 2^exponent (mpz_get_d_2exp); returns 0 for zero.
    /// </summary>
    public double ToDouble(out long exponent)
    {
        if (_size == 0)
        {
            exponent = 0;
            return 0;
        }
        long bits = GetBitLength();
        exponent = bits;
        ulong top = bits <= 53 ? Limbs[0] << (int)(53 - bits) : ExtractBitsAt(Limbs, bits - 53);
        double m = Math.ScaleB((double)(top & ((1UL << 53) - 1)), -53);
        return _size < 0 ? -m : m;
    }

    public bool FitsInt64 => LimbCount <= 1 && (_size == 0 || (_size > 0 ? _d![0] <= long.MaxValue : _d![0] <= (ulong)long.MaxValue + 1));
    public bool FitsUInt64 => _size >= 0 && _size <= 1;
    public bool FitsInt32 => FitsInt64 && (long)this is >= int.MinValue and <= int.MaxValue;
    public bool FitsUInt32 => _size >= 0 && (_size == 0 || (_size == 1 && _d![0] <= uint.MaxValue));

    // ---------------------------------------------------------------
    // Prime search with a small-prime sieve (mpz_nextprime / mpz_prevprime)
    // ---------------------------------------------------------------

    /// <summary>Largest probable prime strictly below this value (mpz_prevprime); false when there is none (value &lt;= 2).</summary>
    public bool TryPrevPrime(out BigInt prime)
    {
        prime = default;
        if (CompareTo(3L) <= 0)
        {
            if (CompareTo(3L) == 0)
            {
                prime = Two;
                return true;
            }
            return false;
        }
        BigInt c = this - One;
        if (c.IsEven)
        {
            c -= One;
        }
        prime = SieveSearch(c, step: -2);
        return true;
    }

    /// <summary>Largest probable prime strictly below this value (mpz_prevprime); throws when the value is &lt;= 2.</summary>
    public BigInt PrevPrime() => TryPrevPrime(out BigInt p) ? p : throw new ArithmeticException("No prime below the value.");

    private const int SievePrimeCount = 600;
    private static readonly uint[] SievePrimes = SmallPrimes.PrimesUpTo(4500)[1..SievePrimeCount];

    /// <summary>Walks odd candidates from c (odd) by step ±2, skipping those with a small factor via incremental residues.</summary>
    private static BigInt SieveSearch(BigInt c, int step)
    {
        if (c.CompareTo(new BigInt(SievePrimes[^1]) * SievePrimes[^1]) <= 0)
        {
            while (!c.IsProbablePrime())
            {
                c += step;
            }
            return c;
        }
        uint[] primes = SievePrimes;
        uint[] res = new uint[primes.Length];
        for (int i = 0; i < primes.Length; i++)
        {
            res[i] = (uint)RemainderAbs(c, primes[i]);
        }
        long offset = 0;
        while (true)
        {
            bool composite = false;
            for (int i = 0; i < primes.Length; i++)
            {
                if (res[i] == 0)
                {
                    composite = true;
                    break;
                }
            }
            if (!composite)
            {
                BigInt cand = c + new BigInt(offset);
                if (cand.IsProbablePrime())
                {
                    return cand;
                }
            }
            offset += step;
            for (int i = 0; i < primes.Length; i++)
            {
                uint p = primes[i];
                uint r = step > 0 ? res[i] + 2 : res[i] + p - 2;
                res[i] = r >= p ? r - p : r;
            }
        }
    }
}
