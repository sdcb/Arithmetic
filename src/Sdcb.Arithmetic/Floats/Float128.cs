using System.Diagnostics;
using System.Runtime.CompilerServices;

namespace Sdcb.Arithmetic;

/// <summary>
/// An IEEE 754-2019 binary128 (quadruple precision) floating-point number: 1 sign bit, 15 exponent
/// bits (bias 16383) and a 112-bit fraction, i.e. a 113-bit significand. Bit-compatible with C
/// <c>__float128</c>/<c>_Float128</c> and Fortran <c>real(16)</c>.
/// </summary>
/// <remarks>
/// <para>
/// The basic operations (<c>+ - * /</c>, <see cref="Sqrt(Float128, RealRounding)"/>,
/// <see cref="FusedMultiplyAdd(Float128, Float128, Float128, RealRounding)"/>, remainders, rounding to
/// integers and all conversions) are correctly rounded in pure integer arithmetic with no allocations,
/// including subnormals, signed zeros, infinities and overflow per rounding mode. Operators round to
/// nearest-even; the static methods take an optional <see cref="RealRounding"/>.
/// </para>
/// <para>
/// NaN: every operation with a NaN operand (and every invalid operation) returns the canonical
/// positive quiet NaN (<see cref="NaN"/>); payloads and NaN signs are not propagated. Sign-bit
/// operations (negation, <see cref="Abs"/>, <see cref="CopySign"/>) and the bitwise operators act on
/// the raw bits, NaN included.
/// </para>
/// <para>
/// Elementary functions (exp, log, trigonometric, ...) are evaluated through <see cref="Real"/> and are
/// also correctly rounded (see <see cref="Exp(Float128, RealRounding)"/>).
/// </para>
/// </remarks>
public readonly partial struct Float128 : IEquatable<Float128>, IComparable<Float128>, IComparable
{
    private readonly UInt128 _bits;

    private Float128(UInt128 bits) => _bits = bits;

    // ---------------------------------------------------------------
    // Layout constants
    // ---------------------------------------------------------------

    internal const int SignificandBits = 113;
    internal const int FractionBits = 112;
    internal const int ExponentBias = 16383;
    internal const int MaxBiasedExponent = 0x7FFF;

    private const ulong SignMaskHi = 0x8000_0000_0000_0000;
    private const ulong ExpMaskHi = 0x7FFF_0000_0000_0000;
    private const ulong QuietBitHi = 0x0000_8000_0000_0000;
    private const ulong HiddenHi = 0x0001_0000_0000_0000;
    private const ulong FracMaskHi = 0x0000_FFFF_FFFF_FFFF;

    private static UInt128 SignMask { [MethodImpl(MethodImplOptions.AggressiveInlining)] get => new(SignMaskHi, 0); }
    private static UInt128 FracMask { [MethodImpl(MethodImplOptions.AggressiveInlining)] get => new(FracMaskHi, ulong.MaxValue); }
    private static UInt128 Hidden { [MethodImpl(MethodImplOptions.AggressiveInlining)] get => new(HiddenHi, 0); }
    private static UInt128 InfBits { [MethodImpl(MethodImplOptions.AggressiveInlining)] get => new(ExpMaskHi, 0); }
    private static UInt128 NaNBits { [MethodImpl(MethodImplOptions.AggressiveInlining)] get => new(ExpMaskHi | QuietBitHi, 0); }
    private static UInt128 MaxBits { [MethodImpl(MethodImplOptions.AggressiveInlining)] get => new(0x7FFE_FFFF_FFFF_FFFF, ulong.MaxValue); }
    private static UInt128 Max113 { [MethodImpl(MethodImplOptions.AggressiveInlining)] get => new(0x0001_FFFF_FFFF_FFFF, ulong.MaxValue); }

    // ---------------------------------------------------------------
    // Raw bits
    // ---------------------------------------------------------------

    /// <summary>Creates a value from its IEEE 754 binary128 bit pattern.</summary>
    public static Float128 FromBits(UInt128 bits) => new(bits);

    /// <summary>The IEEE 754 binary128 bit pattern of this value.</summary>
    public UInt128 ToBits() => _bits;

    // ---------------------------------------------------------------
    // Constants
    // ---------------------------------------------------------------

    /// <summary>+0.</summary>
    public static Float128 Zero => default;

    /// <summary>-0.</summary>
    public static Float128 NegativeZero => new(SignMask);

    /// <summary>1.</summary>
    public static Float128 One => new(new UInt128(0x3FFF_0000_0000_0000, 0));

    /// <summary>-1.</summary>
    public static Float128 NegativeOne => new(new UInt128(0xBFFF_0000_0000_0000, 0));

    /// <summary>π correctly rounded to nearest.</summary>
    public static Float128 Pi => new(new UInt128(0x4000_921F_B544_42D1, 0x8469_898C_C517_01B8));

    /// <summary>2π correctly rounded to nearest.</summary>
    public static Float128 Tau => new(new UInt128(0x4001_921F_B544_42D1, 0x8469_898C_C517_01B8));

    /// <summary>Euler's number e correctly rounded to nearest.</summary>
    public static Float128 E => new(new UInt128(0x4000_5BF0_A8B1_4576, 0x9535_5FB8_AC40_4E7A));

    /// <summary>The smallest positive (subnormal) value, 2^-16494 (the <see cref="double.Epsilon"/> analogue).</summary>
    public static Float128 Epsilon => new(UInt128.One);

    /// <summary>The smallest positive normal value, 2^-16382.</summary>
    public static Float128 MinNormal => new(new UInt128(HiddenHi, 0));

    /// <summary>The largest finite value, (2 - 2^-112) * 2^16383 ≈ 1.19e4932.</summary>
    public static Float128 MaxValue => new(MaxBits);

    /// <summary>The most negative finite value, -<see cref="MaxValue"/>.</summary>
    public static Float128 MinValue => new(MaxBits | SignMask);

    /// <summary>The canonical quiet NaN (positive sign, quiet bit set, zero payload).</summary>
    public static Float128 NaN { [MethodImpl(MethodImplOptions.AggressiveInlining)] get =>  new(NaNBits); }

    /// <summary>+∞.</summary>
    public static Float128 PositiveInfinity => new(InfBits);

    /// <summary>-∞.</summary>
    public static Float128 NegativeInfinity => new(InfBits | SignMask);

    // ---------------------------------------------------------------
    // Field access
    // ---------------------------------------------------------------

    private ulong Hi => (ulong)(_bits >> 64);

    private bool SignBit => (long)Hi < 0;

    private int BiasedExponent => (int)(Hi >> 48) & 0x7FFF;

    private UInt128 Fraction => _bits & FracMask;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int ExpOf(UInt128 b) => (int)((ulong)(b >> 64) >> 48) & 0x7FFF;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static UInt128 FracOf(UInt128 b) => b & FracMask;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static bool SignOf(UInt128 b) => (long)(ulong)(b >> 64) < 0;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static UInt128 SignBits(bool negative) => negative ? SignMask : UInt128.Zero;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static Float128 ZeroOf(bool negative) => new(SignBits(negative));

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static Float128 InfOf(bool negative) => new(InfBits | SignBits(negative));

    /// <summary>Normalizes a nonzero subnormal fraction: returns the significand with bit 112 set and its biased exponent (1 - shift).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static UInt128 NormSubnormal(UInt128 frac, out int exp)
    {
        Debug.Assert(frac != 0 && frac >> 112 == 0);
        int shift = (int)UInt128.LeadingZeroCount(frac) - 15;
        exp = 1 - shift;
        return frac << shift;
    }

    // ---------------------------------------------------------------
    // Rounding core
    // ---------------------------------------------------------------

    /// <summary>
    /// The largest finite value (directed rounding toward zero) or ±∞, for an overflowing result.
    /// </summary>
    private static Float128 Huge(bool negative, RealRounding rounding)
    {
        bool toInf = rounding switch
        {
            RealRounding.NearestEven or RealRounding.AwayFromZero => true,
            RealRounding.TowardZero => false,
            RealRounding.TowardPositive => !negative,
            RealRounding.TowardNegative => negative,
            _ => throw new ArgumentOutOfRangeException(nameof(rounding)),
        };
        return new((toInf ? InfBits : MaxBits) | SignBits(negative));
    }

    /// <summary>Shifts the 192-bit value (sig:extra) right by n, OR-ing every dropped bit into bit 0 of extra.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void ShiftRightJam(ref UInt128 sig, ref ulong extra, int n)
    {
        Debug.Assert(n > 0);
        if (n < 64)
        {
            ulong e = ((ulong)sig << (64 - n)) | (extra >> n) | ((extra << (64 - n)) != 0 ? 1UL : 0UL);
            sig >>= n;
            extra = e;
            return;
        }
        ShiftRightJamLong(ref sig, ref extra, n);
    }

    private static void ShiftRightJamLong(ref UInt128 sig, ref ulong extra, int n)
    {
        if (n >= 192)
        {
            extra = (sig != 0 || extra != 0) ? 1UL : 0UL;
            sig = UInt128.Zero;
            return;
        }
        while (n >= 64)
        {
            extra = (ulong)sig | (extra != 0 ? 1UL : 0UL);
            sig >>= 64;
            n -= 64;
        }
        if (n > 0)
        {
            ShiftRightJam(ref sig, ref extra, n);
        }
    }

    /// <summary>Shifts right by n (n &gt; 0), OR-ing the dropped bits into bit 0.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static UInt128 ShiftRightJam(UInt128 a, int n)
    {
        Debug.Assert(n > 0);
        if (n < 128)
        {
            return (a >> n) | ((a << (128 - n)) != 0 ? UInt128.One : UInt128.Zero);
        }
        return a != 0 ? UInt128.One : UInt128.Zero;
    }

    /// <summary>
    /// Rounds sign * (sig + extra / 2^64) * 2^(exp - 16494) to binary128 and packs it.
    /// Requires sig &lt; 2^113 with bit 112 set whenever exp &gt; 0 (exp is the biased exponent minus one, as
    /// in SoftFloat's roundPackToF128). Handles subnormal results (exp &lt; 0), overflow and every rounding mode.
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static Float128 RoundPack(bool sign, int exp, UInt128 sig, ulong extra, RealRounding rounding)
    {
        if (rounding == RealRounding.NearestEven && (uint)exp < 0x7FFD)
        {
            // Normal result, nearest-even: a carry out of bit 112 bumps the exponent field and cannot reach ∞.
            ulong lo = (ulong)sig, hi = (ulong)(sig >> 64) + ((ulong)(uint)exp << 48);
            if ((long)extra < 0)
            {
                lo++;
                if (lo == 0)
                {
                    hi++;
                }
                if (extra == 0x8000_0000_0000_0000)
                {
                    lo &= ~1UL;
                }
            }
            return new(new UInt128(hi | (sign ? SignMaskHi : 0), lo));
        }
        return RoundPackSlow(sign, exp, sig, extra, rounding);
    }

    [MethodImpl(MethodImplOptions.NoInlining)]
    private static Float128 RoundPackSlow(bool sign, int exp, UInt128 sig, ulong extra, RealRounding rounding)
    {
        if ((uint)exp >= 0x7FFD)
        {
            if (exp < 0)
            {
                ShiftRightJam(ref sig, ref extra, -exp);
                exp = 0;
            }
            else if (exp > 0x7FFD)
            {
                return Huge(sign, rounding);
            }
        }
        if (rounding == RealRounding.NearestEven)
        {
            if ((long)extra < 0)
            {
                sig++;
                if (extra == 0x8000_0000_0000_0000)
                {
                    sig &= ~UInt128.One; // tie: round to even
                }
            }
        }
        else if (extra != 0 && DirectedUp(sign, rounding))
        {
            sig++;
        }
        UInt128 bits = ((UInt128)(uint)exp << 112) + sig;
        if (exp >= 0x7FFD && bits >= InfBits)
        {
            return Huge(sign, rounding);
        }
        return new(bits | SignBits(sign));
    }

    /// <summary>Whether an inexact result rounds away from zero under a directed mode.</summary>
    private static bool DirectedUp(bool sign, RealRounding rounding) => rounding switch
    {
        RealRounding.TowardZero => false,
        RealRounding.TowardPositive => !sign,
        RealRounding.TowardNegative => sign,
        RealRounding.AwayFromZero => true,
        _ => throw new ArgumentOutOfRangeException(nameof(rounding)),
    };

    /// <summary>
    /// Normalizes a nonzero sig (any width) with value sign * sig * 2^(exp - 16494), then rounds and packs it.
    /// Exact when sig has at most 113 significant bits and the result is in range.
    /// </summary>
    private static Float128 NormRoundPack(bool sign, int exp, UInt128 sig, RealRounding rounding)
    {
        Debug.Assert(sig != 0);
        int shift = (int)UInt128.LeadingZeroCount(sig) - 15;
        exp -= shift;
        if (shift >= 0)
        {
            sig <<= shift;
            if ((uint)exp < 0x7FFD)
            {
                return new((((UInt128)(uint)exp << 112) + sig) | SignBits(sign));
            }
            return RoundPack(sign, exp, sig, 0, rounding);
        }
        int r = -shift;
        ulong extra = (ulong)sig << (64 - r);
        sig >>= r;
        return RoundPack(sign, exp, sig, extra, rounding);
    }

    /// <summary>
    /// Rounds sign * (m + θ) * 2^e2 (θ ∈ (0,1) iff sticky; m &gt; 0). When sticky is set m must carry at least
    /// 115 significant bits (or the result must be far below the subnormal range).
    /// </summary>
    internal static Float128 RoundWide(bool negative, UInt128 m, long e2, bool sticky, RealRounding rounding)
    {
        Debug.Assert(m != 0);
        int p = 127 - (int)UInt128.LeadingZeroCount(m);
        UInt128 sig;
        ulong extra;
        long exp;
        if (p > 112)
        {
            int r = p - 112;
            extra = (ulong)(m & ((UInt128.One << r) - 1)) << (64 - r);
            sig = m >> r;
            exp = e2 + r + 16494;
        }
        else
        {
            sig = m << (112 - p);
            extra = 0;
            exp = e2 - (112 - p) + 16494;
        }
        if (sticky)
        {
            extra |= 1;
        }
        exp = Math.Clamp(exp, -100_000, 100_000);
        return RoundPack(negative, (int)exp, sig, extra, rounding);
    }

    /// <summary>Rounds sign * (mag + θ) * 2^e2 for an arbitrary BigInt magnitude (θ ∈ (0,1) iff sticky).</summary>
    internal static Float128 RoundBig(bool negative, BigInt mag, long e2, bool sticky, RealRounding rounding)
    {
        Debug.Assert(mag.Sign > 0);
        long bits = mag.GetBitLength();
        if (bits <= 128)
        {
            return RoundWide(negative, (UInt128)mag, e2, sticky, rounding);
        }
        long shift = bits - 128;
        UInt128 top = (UInt128)mag.ShiftRightTruncate(shift);
        bool st = sticky || mag.TrailingZeroCount() < shift;
        return RoundWide(negative, top, e2 + shift, st, rounding);
    }

    /// <summary>Correctly rounds sign * num / den * 2^e2 (num, den &gt; 0).</summary>
    internal static Float128 RoundRational(bool negative, BigInt num, BigInt den, long e2, RealRounding rounding)
    {
        Debug.Assert(num.Sign > 0 && den.Sign > 0);
        long e0 = num.GetBitLength() - den.GetBitLength() + e2; // value in (2^(e0-1), 2^(e0+1))
        if (e0 > 16390)
        {
            return Huge(negative, rounding);
        }
        if (e0 < -16600)
        {
            return RoundWide(negative, UInt128.One, -20000, sticky: true, rounding);
        }
        long u = num.GetBitLength() - den.GetBitLength() - 120;
        BigInt q, r;
        if (u <= 0)
        {
            BigInt.DivRem(num.ShiftLeft(-u), den, out q, out r);
        }
        else
        {
            BigInt.DivRem(num, den.ShiftLeft(u), out q, out r);
        }
        return RoundBig(negative, q, u + e2, !r.IsZero, rounding);
    }

    // ---------------------------------------------------------------
    // Classification
    // ---------------------------------------------------------------

    /// <summary>True for NaN.</summary>
    public static bool IsNaN(Float128 value) => (value._bits & ~SignMask) > InfBits;

    /// <summary>True for ±∞.</summary>
    public static bool IsInfinity(Float128 value) => (value._bits & ~SignMask) == InfBits;

    /// <summary>True for +∞.</summary>
    public static bool IsPositiveInfinity(Float128 value) => value._bits == InfBits;

    /// <summary>True for -∞.</summary>
    public static bool IsNegativeInfinity(Float128 value) => value._bits == (InfBits | SignMask);

    /// <summary>True for finite values (zero, subnormal or normal).</summary>
    public static bool IsFinite(Float128 value) => value.BiasedExponent != MaxBiasedExponent;

    /// <summary>True for ±0.</summary>
    public static bool IsZero(Float128 value) => (value._bits & ~SignMask) == 0;

    /// <summary>True for normal (finite, nonzero, not subnormal) values.</summary>
    public static bool IsNormal(Float128 value)
    {
        int e = value.BiasedExponent;
        return e != 0 && e != MaxBiasedExponent;
    }

    /// <summary>True for subnormal values.</summary>
    public static bool IsSubnormal(Float128 value) => value.BiasedExponent == 0 && value.Fraction != 0;

    /// <summary>True when the sign bit is set (including -0, -∞ and negative NaNs).</summary>
    public static bool IsNegative(Float128 value) => value.SignBit;

    /// <summary>True when the sign bit is clear (including +0 and +∞).</summary>
    public static bool IsPositive(Float128 value) => !value.SignBit;

    // ---------------------------------------------------------------
    // Comparison and equality (IEEE 754 semantics: NaN is unordered, -0 == +0)
    // ---------------------------------------------------------------

    /// <summary>IEEE equality: false when either operand is NaN; -0 == +0.</summary>
    public static bool operator ==(Float128 left, Float128 right)
    {
        if (IsNaN(left) || IsNaN(right))
        {
            return false;
        }
        return left._bits == right._bits || ((left._bits | right._bits) & ~SignMask) == 0;
    }

    /// <summary>IEEE inequality: true when either operand is NaN.</summary>
    public static bool operator !=(Float128 left, Float128 right) => !(left == right);

    /// <summary>IEEE less-than (false when either operand is NaN).</summary>
    public static bool operator <(Float128 left, Float128 right)
    {
        if (IsNaN(left) || IsNaN(right))
        {
            return false;
        }
        return LessThanOrdered(left._bits, right._bits);
    }

    /// <summary>IEEE greater-than (false when either operand is NaN).</summary>
    public static bool operator >(Float128 left, Float128 right) => right < left;

    /// <summary>IEEE less-than-or-equal (false when either operand is NaN).</summary>
    public static bool operator <=(Float128 left, Float128 right)
    {
        if (IsNaN(left) || IsNaN(right))
        {
            return false;
        }
        return !LessThanOrdered(right._bits, left._bits);
    }

    /// <summary>IEEE greater-than-or-equal (false when either operand is NaN).</summary>
    public static bool operator >=(Float128 left, Float128 right) => right <= left;

    private static bool LessThanOrdered(UInt128 a, UInt128 b)
    {
        bool sa = SignOf(a), sb = SignOf(b);
        if (sa != sb)
        {
            return sa && ((a | b) & ~SignMask) != 0;
        }
        return a != b && (sa ^ (a < b));
    }

    /// <summary>
    /// The IEEE 754 totalOrder relation as a comparison: -NaN &lt; -∞ &lt; ... &lt; -0 &lt; +0 &lt; ... &lt; +∞ &lt; +NaN
    /// (NaNs ordered by payload). Returns a negative, zero or positive value.
    /// </summary>
    public static int CompareTotalOrder(Float128 x, Float128 y)
    {
        static UInt128 Key(UInt128 b) => SignOf(b) ? ~b : b | SignMask;
        return Key(x._bits).CompareTo(Key(y._bits));
    }

    /// <summary>Compares like <see cref="double.CompareTo(double)"/>: NaN sorts below everything and equals NaN; -0 equals +0.</summary>
    public int CompareTo(Float128 other)
    {
        if (this < other)
        {
            return -1;
        }
        if (this > other)
        {
            return 1;
        }
        if (this == other)
        {
            return 0;
        }
        return IsNaN(this) ? (IsNaN(other) ? 0 : -1) : 1;
    }

    /// <inheritdoc/>
    public int CompareTo(object? obj)
    {
        if (obj is null)
        {
            return 1;
        }
        if (obj is Float128 f)
        {
            return CompareTo(f);
        }
        throw new ArgumentException("Object must be of type Float128.", nameof(obj));
    }

    /// <summary>Equality like <see cref="double.Equals(double)"/>: NaN equals NaN, -0 equals +0.</summary>
    public bool Equals(Float128 other) => this == other || (IsNaN(this) && IsNaN(other));

    /// <inheritdoc/>
    public override bool Equals(object? obj) => obj is Float128 f && Equals(f);

    /// <summary>Hash consistent with <see cref="Equals(Float128)"/> (all NaNs and both zeros hash alike).</summary>
    public override int GetHashCode()
    {
        UInt128 b = _bits;
        if (IsNaN(this) || IsZero(this))
        {
            b &= InfBits;
        }
        return b.GetHashCode();
    }
}
