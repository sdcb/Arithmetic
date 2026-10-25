using System.Diagnostics;

namespace Sdcb.Arithmetic;

/// <summary>An IEEE 754 binary interchange format: significand bits, exponent of the smallest subnormal ulp, and the overflow threshold 2^MaxExponent.</summary>
internal readonly record struct IeeeFormat(int Precision, int MinUlpExponent, int MaxExponent)
{
    public static readonly IeeeFormat Binary64 = new(53, -1074, 1024);
    public static readonly IeeeFormat Binary32 = new(24, -149, 128);
    public static readonly IeeeFormat Binary16 = new(11, -24, 16);
}

public readonly partial struct Real
{
    /// <summary>
    /// Correctly rounds sign * (mag + θ) * 2^exp (θ ∈ (0,1) iff sticky) to an IEEE binary format,
    /// including subnormals and overflow, and returns it as a double (exact for every narrower
    /// format, so a final cast to float/Half is exact). When sticky is set the caller must
    /// supply at least two more bits than the format keeps.
    /// </summary>
    internal static double RoundToIeee(BigInt mag, long exp, bool sticky, bool negative, IeeeFormat f, RealRounding rnd)
    {
        Debug.Assert(mag.Sign > 0);
        long bits = mag.GetBitLength();
        if (exp + bits > f.MaxExponent + 1L)
        {
            return IeeeHuge(negative, rnd, f);
        }
        long drop = Math.Max(bits - f.Precision, f.MinUlpExponent - exp);
        ulong q;
        long e;
        if (drop <= 0)
        {
            Debug.Assert(!sticky);
            q = (ulong)mag;
            e = exp;
        }
        else
        {
            q = drop >= bits ? 0 : (ulong)mag.ShiftRightTruncate(drop);
            bool guard = drop - 1 < bits && mag.TestBit(drop - 1);
            bool rest = sticky || mag.TrailingZeroCount() < drop - 1;
            e = exp + drop;
            if (guard || rest)
            {
                bool up = rnd switch
                {
                    RealRounding.NearestEven => guard && (rest || (q & 1) != 0),
                    RealRounding.TowardZero => false,
                    RealRounding.AwayFromZero => true,
                    RealRounding.TowardPositive => !negative,
                    RealRounding.TowardNegative => negative,
                    _ => throw new ArgumentOutOfRangeException(nameof(rnd)),
                };
                if (up)
                {
                    q++;
                }
            }
        }
        if (q == 0)
        {
            return negative ? -0.0 : 0.0;
        }
        if (e + (64 - System.Numerics.BitOperations.LeadingZeroCount(q)) > f.MaxExponent)
        {
            return IeeeHuge(negative, rnd, f);
        }
        double d = Math.ScaleB((double)q, (int)e);
        return negative ? -d : d;
    }

    private static double IeeeHuge(bool negative, RealRounding rnd, IeeeFormat f)
    {
        bool toInf = rnd switch
        {
            RealRounding.TowardZero => false,
            RealRounding.TowardPositive => !negative,
            RealRounding.TowardNegative => negative,
            _ => true,
        };
        if (toInf)
        {
            return negative ? double.NegativeInfinity : double.PositiveInfinity;
        }
        double max = Math.ScaleB((double)((1UL << f.Precision) - 1), f.MaxExponent - f.Precision);
        return negative ? -max : max;
    }

    // ---------------------------------------------------------------
    // IEEE 754 binary128 (quadruple precision) interchange
    // ---------------------------------------------------------------

    private const int Binary128Precision = 113;
    private const int Binary128ExponentBias = 16383;
    private const long Binary128MinUlpExponent = -16494;   // 2^-16494 is the smallest subnormal
    private const long Binary128MaxExponent = 16384;       // values >= 2^16384 overflow
    private static readonly UInt128 s_binary128SignBit = UInt128.One << 127;
    private static readonly UInt128 s_binary128FractionMask = (UInt128.One << 112) - 1;
    private static readonly UInt128 s_binary128Infinity = (UInt128)0x7FFF << 112;
    private static readonly UInt128 s_binary128QuietNaN = ((UInt128)0x7FFF << 112) | (UInt128.One << 111);

    /// <summary>
    /// Decodes an IEEE 754 binary128 (quadruple precision) value from its bit pattern, e.g. the
    /// bits of a C <c>__float128</c>/<c>_Float128</c> or Fortran <c>real(16)</c>.
    /// The conversion is exact when <paramref name="precision"/> is at least 113 bits (the default
    /// is <c>max(DefaultPrecision, 113)</c>); otherwise the value is rounded with <paramref name="rounding"/>.
    /// NaN payloads and the NaN sign are not preserved.
    /// </summary>
    public static Real FromBinary128Bits(UInt128 bits, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = precision != 0 ? precision : Math.Max(DefaultPrecision, Binary128Precision);
        ArgumentOutOfRangeException.ThrowIfLessThan(prec, MinPrecision);
        bool neg = (bits & s_binary128SignBit) != 0;
        int biased = (int)(bits >> 112) & 0x7FFF;
        UInt128 frac = bits & s_binary128FractionMask;
        if (biased == 0x7FFF)
        {
            return frac != 0 ? MakeNaN(prec) : MakeInf(neg, prec);
        }
        if (biased == 0)
        {
            return frac == 0
                ? MakeZero(neg, prec)
                : Round(new BigInt(frac), Binary128MinUlpExponent, sticky: false, neg, prec, rounding, out _);
        }
        UInt128 mant = frac | (UInt128.One << 112);
        long e = biased - Binary128ExponentBias - 112L;
        return Round(new BigInt(mant), e, sticky: false, neg, prec, rounding, out _);
    }

    /// <summary>
    /// Encodes this value as an IEEE 754 binary128 (quadruple precision) bit pattern, correctly rounded
    /// with <paramref name="rounding"/>, including subnormals and overflow to ±∞ (or the largest finite
    /// value for directed roundings toward zero). NaN encodes as the positive quiet NaN.
    /// </summary>
    public UInt128 ToBinary128Bits(RealRounding rounding = RealRounding.NearestEven)
    {
        UInt128 sign = IsNeg ? s_binary128SignBit : UInt128.Zero;
        switch (_kind)
        {
            case Kind.NaN:
                return s_binary128QuietNaN;
            case Kind.Inf:
                return sign | s_binary128Infinity;
            case Kind.Zero:
                return sign;
        }

        bool negative = IsNeg;
        BigInt mag = _man;
        long exp = _exp;
        long bits = mag.GetBitLength();
        if (exp + bits > Binary128MaxExponent + 1L)
        {
            return Binary128Huge(negative, rounding);
        }
        long drop = Math.Max(bits - Binary128Precision, Binary128MinUlpExponent - exp);
        UInt128 q;
        long e;
        if (drop <= 0)
        {
            q = (UInt128)mag;
            e = exp;
        }
        else
        {
            q = drop >= bits ? UInt128.Zero : (UInt128)mag.ShiftRightTruncate(drop);
            bool guard = drop - 1 < bits && mag.TestBit(drop - 1);
            bool rest = mag.TrailingZeroCount() < drop - 1;
            e = exp + drop;
            if (guard || rest)
            {
                bool up = rounding switch
                {
                    RealRounding.NearestEven => guard && (rest || (q & 1) != 0),
                    RealRounding.TowardZero => false,
                    RealRounding.AwayFromZero => true,
                    RealRounding.TowardPositive => !negative,
                    RealRounding.TowardNegative => negative,
                    _ => throw new ArgumentOutOfRangeException(nameof(rounding)),
                };
                if (up)
                {
                    q++;
                }
            }
        }
        if (q == 0)
        {
            return sign;
        }

        // TrimmedLength q to 113 bits (or leave it subnormal at the minimum exponent).
        int qBits = 128 - (int)UInt128.LeadingZeroCount(q);
        if (qBits > Binary128Precision)
        {
            Debug.Assert(qBits == Binary128Precision + 1 && (q & 1) == 0);
            q >>= 1;
            e++;
        }
        else if (qBits < Binary128Precision && e > Binary128MinUlpExponent)
        {
            int shift = (int)Math.Min(Binary128Precision - qBits, e - Binary128MinUlpExponent);
            q <<= shift;
            e -= shift;
        }
        if (e + Binary128Precision > Binary128MaxExponent)
        {
            return Binary128Huge(negative, rounding);
        }
        if (q >> 112 == 0)
        {
            Debug.Assert(e == Binary128MinUlpExponent);
            return sign | q; // subnormal
        }
        UInt128 biased = (UInt128)(ulong)(e + 112 + Binary128ExponentBias);
        return sign | (biased << 112) | (q & s_binary128FractionMask);
    }

    private static UInt128 Binary128Huge(bool negative, RealRounding rnd)
    {
        UInt128 sign = negative ? s_binary128SignBit : UInt128.Zero;
        bool toInf = rnd switch
        {
            RealRounding.TowardZero => false,
            RealRounding.TowardPositive => !negative,
            RealRounding.TowardNegative => negative,
            _ => true,
        };
        return toInf ? sign | s_binary128Infinity : sign | (s_binary128Infinity - 1);
    }

    /// <summary>Correctly rounds the positive rational num/den (round to nearest, ties to even) to an IEEE binary format.</summary>
    internal static double RationalToIeee(BigInt num, BigInt den, bool negative, IeeeFormat f)
    {
        Debug.Assert(num.Sign > 0 && den.Sign > 0);
        long e0 = num.GetBitLength() - den.GetBitLength(); // value in (2^(e0-1), 2^(e0+1))
        if (e0 > f.MaxExponent + 2L)
        {
            return IeeeHuge(negative, RealRounding.NearestEven, f);
        }
        if (e0 + 1 < f.MinUlpExponent - 1L)
        {
            return negative ? -0.0 : 0.0; // below half the smallest subnormal
        }
        long u = e0 - f.Precision - 3; // quotient carries >= Precision + 3 bits
        BigInt q, r;
        if (u <= 0)
        {
            BigInt.DivRem(num.ShiftLeft(-u), den, out q, out r);
        }
        else
        {
            BigInt.DivRem(num, den.ShiftLeft(u), out q, out r);
        }
        return RoundToIeee(q, u, !r.IsZero, negative, f, RealRounding.NearestEven);
    }

    /// <summary>Exact decimal → rational (mantissa / 10^scale).</summary>
    internal static BigRational DecimalToRational(decimal value)
    {
        Span<int> bits = stackalloc int[4];
        decimal.GetBits(value, bits);
        UInt128 mant = ((UInt128)(uint)bits[2] << 64) | ((UInt128)(uint)bits[1] << 32) | (uint)bits[0];
        int scale = (bits[3] >> 16) & 0xFF;
        BigInt m = new(mant);
        if (bits[3] < 0)
        {
            m = -m;
        }
        return scale == 0 ? new BigRational(m) : new BigRational(m, BigInt.Pow(BigInt.Ten, (ulong)scale));
    }

    private static readonly BigInt s_two96 = BigInt.One << 96;

    /// <summary>Nearest decimal (ties to even) to num/den; false when the magnitude does not fit.</summary>
    internal static bool TryRationalToDecimal(BigInt num, BigInt den, out decimal result)
    {
        result = 0m;
        if (num.IsZero)
        {
            return true;
        }
        bool neg = num.Sign < 0;
        BigInt a = num.Abs();
        for (int scale = 28; scale >= 0; scale--)
        {
            BigInt scaled = scale == 0 ? a : a * BigInt.Pow(BigInt.Ten, (ulong)scale);
            BigInt.DivRem(scaled, den, out BigInt q, out BigInt r);
            int c = (r << 1).CompareTo(den);
            if (c > 0 || (c == 0 && q.IsOdd))
            {
                q += BigInt.One;
            }
            if (q >= s_two96)
            {
                continue;
            }
            int s = scale;
            while (s > 0 && !q.IsZero && BigInt.RemainderAbs(q, 10) == 0)
            {
                q = BigInt.DivRem(q, 10UL, out _);
                s--;
            }
            if (q.IsZero)
            {
                return true;
            }
            UInt128 m = (UInt128)q;
            result = new decimal((int)(uint)m, (int)(uint)(m >> 32), (int)(uint)(m >> 64), neg, (byte)s);
            return true;
        }
        return false;
    }
}
