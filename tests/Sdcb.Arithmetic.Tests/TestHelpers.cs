using System.Numerics;

namespace Sdcb.Arithmetic.Tests;

/// <summary>Shared helpers: deterministic random values mirrored as System.Numerics.BigInteger.</summary>
public static class TestHelpers
{
    public static (BigInt A, BigInteger B) RandomPair(Random rng, int maxLimbs, bool allowNegative = true)
    {
        int limbs = rng.Next(0, maxLimbs + 1);
        byte[] bytes = new byte[limbs * 8 + rng.Next(0, 8)];
        rng.NextBytes(bytes);
        // Occasionally force runs of 0x00/0xFF to exercise carry chains.
        if (bytes.Length > 16 && rng.Next(4) == 0)
        {
            int start = rng.Next(bytes.Length / 2);
            int len = rng.Next(bytes.Length - start);
            bytes.AsSpan(start, len).Fill(rng.Next(2) == 0 ? (byte)0x00 : (byte)0xFF);
        }
        int sign = allowNegative && rng.Next(2) == 0 ? -1 : 1;
        BigInt a = BigInt.FromMagnitudeBytes(bytes, sign);
        BigInteger b = new(bytes, isUnsigned: true, isBigEndian: false);
        if (sign < 0)
        {
            b = -b;
        }
        return (a, b);
    }

    public static BigInteger ToBigInteger(this BigInt v)
    {
        byte[] mag = v.ToMagnitudeBytes();
        BigInteger b = new(mag, isUnsigned: true, isBigEndian: false);
        return v.Sign < 0 ? -b : b;
    }

    public static BigInt ToBigInt(this BigInteger v)
    {
        byte[] mag = BigInteger.Abs(v).ToByteArray(isUnsigned: true, isBigEndian: false);
        return BigInt.FromMagnitudeBytes(mag, v.Sign < 0 ? -1 : 1);
    }

    public static void AssertEqual(BigInteger expected, BigInt actual, string context)
    {
        BigInteger actualBI = actual.ToBigInteger();
        if (expected != actualBI)
        {
            throw new Xunit.Sdk.XunitException(
                $"{context}\nexpected: {expected}\nactual:   {actualBI}");
        }
    }
}
