using System.Runtime.InteropServices;

namespace Sdcb.Arithmetic;

public readonly partial struct BigInt
{
    /// <summary>Uniformly random BigInt in [0, 2^bits) (mpz_urandomb).</summary>
    public static BigInt Random(System.Random random, long bits)
    {
        ArgumentNullException.ThrowIfNull(random);
        ArgumentOutOfRangeException.ThrowIfNegative(bits);
        if (bits == 0)
        {
            return default;
        }
        int limbs = checked((int)((bits + 63) >> 6));
        ulong[] d = new ulong[limbs];
        random.NextBytes(MemoryMarshal.AsBytes(d.AsSpan()));
        int topBits = (int)(bits & 63);
        if (topBits != 0)
        {
            d[limbs - 1] &= (1UL << topBits) - 1;
        }
        return FromMagnitude(d, limbs, 1);
    }

    /// <summary>Uniformly random BigInt in [0, bound) (mpz_urandomm).</summary>
    public static BigInt RandomBelow(System.Random random, BigInt bound)
    {
        ArgumentNullException.ThrowIfNull(random);
        if (bound.Sign <= 0)
        {
            throw new ArgumentOutOfRangeException(nameof(bound), "Bound must be positive.");
        }
        long bits = bound.GetBitLength();
        // Rejection sampling; expected < 2 iterations.
        while (true)
        {
            BigInt candidate = Random(random, bits);
            if (candidate.CompareTo(bound) < 0)
            {
                return candidate;
            }
        }
    }
}
