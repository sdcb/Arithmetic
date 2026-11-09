using System.Diagnostics;
using System.Runtime.CompilerServices;
using Sdcb.Arithmetic.Core;

namespace Sdcb.Arithmetic;

// Core arithmetic in the style of Berkeley SoftFloat 3e (f128_add/sub/mul/div/sqrt/mulAdd):
// unpack to (sign, biased exponent, 113-bit significand), compute exactly with a 64-bit
// guard/sticky word, then one shared RoundPack. Exponent convention inside RoundPack:
// value = sig * 2^(exp - 16494), i.e. exp is the biased exponent minus one for a normalized sig.
public readonly partial struct Float128
{
    // ---------------------------------------------------------------
    // Operators (round to nearest, ties to even)
    // ---------------------------------------------------------------

    /// <summary>Correctly rounded sum (nearest-even).</summary>
    public static Float128 operator +(Float128 left, Float128 right) => Add(left, right);

    /// <summary>Correctly rounded difference (nearest-even).</summary>
    public static Float128 operator -(Float128 left, Float128 right) => Subtract(left, right);

    /// <summary>Correctly rounded product (nearest-even).</summary>
    public static Float128 operator *(Float128 left, Float128 right) => Multiply(left, right);

    /// <summary>Correctly rounded quotient (nearest-even).</summary>
    public static Float128 operator /(Float128 left, Float128 right) => Divide(left, right);

    /// <summary>Truncated remainder with the sign of the dividend (C fmod, like <see cref="double"/> %); always exact.</summary>
    public static Float128 operator %(Float128 left, Float128 right) => Remainder(left, right, ieee: false);

    /// <summary>Flips the sign bit.</summary>
    public static Float128 operator -(Float128 value) => new(value._bits ^ SignMask);

    /// <summary>Returns the value unchanged.</summary>
    public static Float128 operator +(Float128 value) => value;

    /// <summary>value + 1 (nearest-even).</summary>
    public static Float128 operator ++(Float128 value) => Add(value, One);

    /// <summary>value - 1 (nearest-even).</summary>
    public static Float128 operator --(Float128 value) => Subtract(value, One);

    // ---------------------------------------------------------------
    // Addition / subtraction
    // ---------------------------------------------------------------

    /// <summary>Correctly rounded a + b.</summary>
    public static Float128 Add(Float128 a, Float128 b, RealRounding rounding = RealRounding.NearestEven)
    {
        bool signA = a.SignBit;
        return signA == b.SignBit
            ? AddMags(a._bits, b._bits, signA, rounding)
            : SubMags(a._bits, b._bits, signA, rounding);
    }

    /// <summary>Correctly rounded a - b.</summary>
    public static Float128 Subtract(Float128 a, Float128 b, RealRounding rounding = RealRounding.NearestEven)
    {
        bool signA = a.SignBit;
        return signA != b.SignBit
            ? AddMags(a._bits, b._bits, signA, rounding)
            : SubMags(a._bits, b._bits, signA, rounding);
    }

    private static Float128 AddMags(UInt128 uiA, UInt128 uiB, bool signZ, RealRounding rounding)
    {
        int expA = ExpOf(uiA), expB = ExpOf(uiB);
        UInt128 sigA = FracOf(uiA), sigB = FracOf(uiB);
        int expDiff = expA - expB;
        ulong extra = 0;
        int expZ;
        UInt128 sigZ;
        if (expDiff == 0)
        {
            if (expA == MaxBiasedExponent)
            {
                return (sigA | sigB) != 0 ? NaN : new(uiA);
            }
            sigZ = sigA + sigB;
            if (expA == 0)
            {
                return new(sigZ | SignBits(signZ)); // subnormal sum is exact; a carry becomes the normal exponent
            }
            expZ = expA;
            sigZ |= Hidden << 1;
            ShiftRightJam(ref sigZ, ref extra, 1);
            return RoundPack(signZ, expZ, sigZ, extra, rounding);
        }
        if (expDiff < 0)
        {
            if (expB == MaxBiasedExponent)
            {
                return sigB != 0 ? NaN : InfOf(signZ);
            }
            expZ = expB;
            if (expA != 0)
            {
                sigA |= Hidden;
            }
            else
            {
                ++expDiff;
            }
            if (expDiff != 0)
            {
                ShiftRightJam(ref sigA, ref extra, -expDiff);
            }
        }
        else
        {
            if (expA == MaxBiasedExponent)
            {
                return sigA != 0 ? NaN : InfOf(signZ);
            }
            expZ = expA;
            if (expB != 0)
            {
                sigB |= Hidden;
            }
            else
            {
                --expDiff;
            }
            if (expDiff != 0)
            {
                ShiftRightJam(ref sigB, ref extra, expDiff);
            }
        }
        // The larger operand's hidden bit is added through sigA (the smaller one is already below 2^112).
        sigZ = (sigA | Hidden) + sigB;
        --expZ;
        if (sigZ >> 113 != 0)
        {
            ++expZ;
            ShiftRightJam(ref sigZ, ref extra, 1);
        }
        return RoundPack(signZ, expZ, sigZ, extra, rounding);
    }

    private static Float128 SubMags(UInt128 uiA, UInt128 uiB, bool signZ, RealRounding rounding)
    {
        int expA = ExpOf(uiA), expB = ExpOf(uiB);
        UInt128 sigA = FracOf(uiA) << 4, sigB = FracOf(uiB) << 4;
        UInt128 hidden4 = Hidden << 4;
        int expDiff = expA - expB;
        int expZ;
        UInt128 sigZ;
        if (expDiff > 0)
        {
            if (expA == MaxBiasedExponent)
            {
                return sigA != 0 ? NaN : new(uiA);
            }
            if (expB != 0)
            {
                sigB |= hidden4;
            }
            else
            {
                --expDiff;
            }
            if (expDiff != 0)
            {
                sigB = ShiftRightJam(sigB, expDiff);
            }
            expZ = expA;
            sigA |= hidden4;
            sigZ = sigA - sigB;
        }
        else if (expDiff < 0)
        {
            if (expB == MaxBiasedExponent)
            {
                return sigB != 0 ? NaN : InfOf(!signZ);
            }
            if (expA != 0)
            {
                sigA |= hidden4;
            }
            else
            {
                ++expDiff;
            }
            if (expDiff != 0)
            {
                sigA = ShiftRightJam(sigA, -expDiff);
            }
            expZ = expB;
            sigB |= hidden4;
            signZ = !signZ;
            sigZ = sigB - sigA;
        }
        else
        {
            if (expA == MaxBiasedExponent)
            {
                return NaN; // ∞ - ∞ or a NaN operand
            }
            expZ = expA == 0 ? 1 : expA;
            if (sigA > sigB)
            {
                sigZ = sigA - sigB;
            }
            else if (sigA < sigB)
            {
                signZ = !signZ;
                sigZ = sigB - sigA;
            }
            else
            {
                return ZeroOf(rounding == RealRounding.TowardNegative);
            }
        }
        return NormRoundPack(signZ, expZ - 5, sigZ, rounding);
    }

    // ---------------------------------------------------------------
    // Multiplication
    // ---------------------------------------------------------------

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static UInt128 Mul64(ulong a, ulong b)
    {
        ulong hi = Math.BigMul(a, b, out ulong lo);
        return new UInt128(hi, lo);
    }

    /// <summary>The exact 256-bit product of two 113-bit significands (result &lt; 2^226).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void Mul113(UInt128 a, UInt128 b, out UInt128 hi, out UInt128 lo)
    {
        ulong aL = (ulong)a, aH = (ulong)(a >> 64), bL = (ulong)b, bH = (ulong)(b >> 64);
        // aH, bH < 2^49, so the high words of the cross products are < 2^49 and their sum cannot overflow.
        ulong p0Hi = Math.BigMul(aL, bL, out ulong p0Lo);
        ulong m1Hi = Math.BigMul(aL, bH, out ulong m1Lo);
        ulong m2Hi = Math.BigMul(aH, bL, out ulong m2Lo);
        ulong p3Hi = Math.BigMul(aH, bH, out ulong p3Lo);
        ulong w1 = p0Hi + m1Lo;
        ulong c = w1 < p0Hi ? 1UL : 0UL;
        ulong w1b = w1 + m2Lo;
        c += w1b < w1 ? 1UL : 0UL;
        ulong w2 = m1Hi + m2Hi + c;
        ulong w2b = w2 + p3Lo;
        ulong w3 = p3Hi + (w2b < w2 ? 1UL : 0UL);
        lo = new UInt128(w1b, p0Lo);
        hi = new UInt128(w3, w2b);
    }

    /// <summary>Correctly rounded a * b.</summary>
    public static Float128 Multiply(Float128 a, Float128 b, RealRounding rounding = RealRounding.NearestEven)
    {
        UInt128 uiA = a._bits, uiB = b._bits;
        bool signZ = SignOf(uiA) ^ SignOf(uiB);
        int expA = ExpOf(uiA), expB = ExpOf(uiB);
        UInt128 sigA = FracOf(uiA), sigB = FracOf(uiB);
        if (expA == MaxBiasedExponent || expB == MaxBiasedExponent)
        {
            if ((expA == MaxBiasedExponent && sigA != 0) || (expB == MaxBiasedExponent && sigB != 0))
            {
                return NaN;
            }
            bool zeroA = expA == 0 && sigA == 0, zeroB = expB == 0 && sigB == 0;
            return zeroA || zeroB ? NaN : InfOf(signZ);
        }
        if (expA == 0)
        {
            if (sigA == 0)
            {
                return ZeroOf(signZ);
            }
            sigA = NormSubnormal(sigA, out expA);
        }
        if (expB == 0)
        {
            if (sigB == 0)
            {
                return ZeroOf(signZ);
            }
            sigB = NormSubnormal(sigB, out expB);
        }
        int expZ = expA + expB - 0x4000;
        Mul113(sigA | Hidden, sigB | Hidden, out UInt128 hi, out UInt128 lo);
        // product in [2^224, 2^226): keep 113 bits, then 64 guard bits, then sticky.
        UInt128 sigZ;
        ulong extra;
        if (hi >> 97 != 0)
        {
            sigZ = (hi << 15) | (lo >> 113);
            extra = (ulong)(lo >> 49) | (((ulong)lo & ((1UL << 49) - 1)) != 0 ? 1UL : 0UL);
            expZ++;
        }
        else
        {
            sigZ = (hi << 16) | (lo >> 112);
            extra = (ulong)(lo >> 48) | (((ulong)lo & ((1UL << 48) - 1)) != 0 ? 1UL : 0UL);
        }
        return RoundPack(signZ, expZ, sigZ, extra, rounding);
    }

    // ---------------------------------------------------------------
    // Division
    // ---------------------------------------------------------------

    /// <summary>
    /// Divides the 192-bit (u2:u1:u0) by the normalized (d1:d0) (bit 127 set), requiring (u2:u1) &lt; (d1:d0):
    /// Knuth's algorithm D step with a 128/64 hardware division (portable fallback in <see cref="Nat"/>).
    /// </summary>
    private static ulong Div3By2(ulong u2, ulong u1, ulong u0, ulong d1, ulong d0, out UInt128 rem)
    {
        Debug.Assert((long)d1 < 0 && new UInt128(u2, u1) < new UInt128(d1, d0));
        ulong q;
        if (u2 < d1)
        {
            q = Nat.DivRem128(u2, u1, d1, out ulong r);
            // D3: qhat is at most two too large; the test leaves it at most one too large.
            UInt128 p = Mul64(q, d0);
            if (p > new UInt128(r, u0))
            {
                q--;
                ulong r2 = r + d1;
                if (r2 >= r && p - d0 > new UInt128(r2, u0))
                {
                    q--;
                }
            }
        }
        else
        {
            q = ulong.MaxValue; // u2 == d1: rhat = u1 + d1 >= 2^63, the D3 test is skipped only if it overflows
            UInt128 rr = (UInt128)u1 + d1;
            if ((ulong)(rr >> 64) == 0)
            {
                ulong r = (ulong)rr;
                UInt128 p = Mul64(q, d0);
                if (p > new UInt128(r, u0))
                {
                    q--;
                    ulong r2 = r + d1;
                    if (r2 >= r && p - d0 > new UInt128(r2, u0))
                    {
                        q--;
                    }
                }
            }
        }
        // D4/D6: multiply-subtract, add back once if q was one too large.
        UInt128 plo = Mul64(q, d0);
        UInt128 sub = Mul64(q, d1) + (plo >> 64);
        ulong lo64 = (ulong)plo;
        ulong r0 = u0 - lo64;
        if (u0 < lo64)
        {
            sub++;
        }
        UInt128 top = new(u2, u1);
        UInt128 rTop = top - sub;
        if (top < sub)
        {
            q--;
            rem = new UInt128((ulong)rTop, r0) + new UInt128(d1, d0);
        }
        else
        {
            Debug.Assert(rTop >> 64 == 0);
            rem = new UInt128((ulong)rTop, r0);
        }
        return q;
    }

    /// <summary>Correctly rounded a / b.</summary>
    public static Float128 Divide(Float128 a, Float128 b, RealRounding rounding = RealRounding.NearestEven)
    {
        UInt128 uiA = a._bits, uiB = b._bits;
        bool signZ = SignOf(uiA) ^ SignOf(uiB);
        int expA = ExpOf(uiA), expB = ExpOf(uiB);
        UInt128 sigA = FracOf(uiA), sigB = FracOf(uiB);
        if (expA == MaxBiasedExponent)
        {
            if (sigA != 0 || expB == MaxBiasedExponent)
            {
                return NaN; // NaN or ∞/∞ (or ∞/NaN)
            }
            return InfOf(signZ);
        }
        if (expB == MaxBiasedExponent)
        {
            return sigB != 0 ? NaN : ZeroOf(signZ);
        }
        if (expB == 0)
        {
            if (sigB == 0)
            {
                return expA == 0 && sigA == 0 ? NaN : InfOf(signZ);
            }
            sigB = NormSubnormal(sigB, out expB);
        }
        if (expA == 0)
        {
            if (sigA == 0)
            {
                return ZeroOf(signZ);
            }
            sigA = NormSubnormal(sigA, out expA);
        }
        int expZ = expA - expB + 0x3FFE;
        sigA |= Hidden;
        sigB |= Hidden;
        // Q = floor(sigA * 2^127 / sigB) in (2^126, 2^128): two 3-by-2 steps on the normalized divisor.
        UInt128 d = sigB << 15;
        ulong d1 = (ulong)(d >> 64), d0 = (ulong)d;
        UInt128 u = sigA << 14;
        ulong q1 = Div3By2((ulong)(u >> 64), (ulong)u, 0, d1, d0, out UInt128 r);
        ulong q0 = Div3By2((ulong)(r >> 64), (ulong)r, 0, d1, d0, out r);
        UInt128 q = new(q1, q0);
        UInt128 sigZ;
        ulong extra;
        if ((long)q1 < 0)
        {
            sigZ = q >> 15;
            extra = q0 << 49;
        }
        else
        {
            sigZ = q >> 14;
            extra = q0 << 50;
            expZ--;
        }
        if (r != 0)
        {
            extra |= 1;
        }
        return RoundPack(signZ, expZ, sigZ, extra, rounding);
    }

    // ---------------------------------------------------------------
    // Square root
    // ---------------------------------------------------------------

    /// <summary>floor(sqrt(m)) and m - root² for m in [2^114, 2^116).</summary>
    private static ulong ISqrt116(UInt128 m, out ulong rem)
    {
        // Estimate from the top 64 bits (error ~2^6), one integer Newton step (error < 2), exact fix-up.
        ulong s = (ulong)(Math.Sqrt((double)(ulong)(m >> 52)) * 67108864.0);
        ulong q = Nat.DivRem128((ulong)(m >> 64), (ulong)m, s, out _);
        s = (ulong)(((UInt128)s + q) >> 1);
        UInt128 sq = Mul64(s, s);
        while (sq > m)
        {
            s--;
            sq = Mul64(s, s);
        }
        while (m - sq > 2 * (UInt128)s)
        {
            s++;
            sq = Mul64(s, s);
        }
        rem = (ulong)(m - sq);
        return s;
    }

    /// <summary>Correctly rounded square root (nearest-even); sqrt(-0) = -0, NaN for negative operands.</summary>
    public static Float128 Sqrt(Float128 x) => Sqrt(x, RealRounding.NearestEven);

    /// <summary>Correctly rounded square root; sqrt(-0) = -0, NaN for negative operands.</summary>
    public static Float128 Sqrt(Float128 x, RealRounding rounding)
    {
        UInt128 ui = x._bits;
        bool sign = SignOf(ui);
        int exp = ExpOf(ui);
        UInt128 sig = FracOf(ui);
        if (exp == MaxBiasedExponent)
        {
            return sig != 0 || sign ? NaN : x;
        }
        if (sign)
        {
            return (exp == 0 && sig == 0) ? x : NaN;
        }
        if (exp == 0)
        {
            if (sig == 0)
            {
                return x;
            }
            sig = NormSubnormal(sig, out exp);
        }
        sig |= Hidden;
        // x = sig * 2^t; N = sig * 2^k with t - k even and N in [2^230, 2^232).
        int t = exp - 16495;
        int k = ((t - 118) & 1) == 0 ? 118 : 119;
        // Zimmermann's SqrtRem with β = 2^58 on N = (sig << (k - 116)) * β^2:
        UInt128 m = sig << (k - 116);
        ulong s1 = ISqrt116(m, out ulong r1);
        UInt128 num = (UInt128)r1 << 58;
        ulong q = Nat.DivRem128((ulong)(num >> 64), (ulong)num, s1 << 1, out ulong u);
        UInt128 s = ((UInt128)s1 << 58) + q;
        Int128 rem = (Int128)((UInt128)u << 58) - (Int128)Mul64(q, q);
        if (rem < 0)
        {
            rem += (Int128)(2 * s - 1);
            s--;
        }
        // s = floor(sqrt(N)) in [2^115, 2^116)
        UInt128 sigZ = s >> 3;
        ulong extra = ((ulong)s << 61) | (rem != 0 ? 1UL : 0UL);
        int expZ = (t - k) / 2 + 3 + 16494;
        return RoundPack(false, expZ, sigZ, extra, rounding);
    }

    // ---------------------------------------------------------------
    // Fused multiply-add
    // ---------------------------------------------------------------

    /// <summary>(left * right) + addend with a single rounding (nearest-even).</summary>
    public static Float128 FusedMultiplyAdd(Float128 left, Float128 right, Float128 addend) =>
        FusedMultiplyAdd(left, right, addend, RealRounding.NearestEven);

    /// <summary>(left * right) + addend with a single rounding.</summary>
    public static Float128 FusedMultiplyAdd(Float128 left, Float128 right, Float128 addend, RealRounding rounding)
    {
        UInt128 uiA = left._bits, uiB = right._bits, uiC = addend._bits;
        bool signP = SignOf(uiA) ^ SignOf(uiB), signC = SignOf(uiC);
        int expA = ExpOf(uiA), expB = ExpOf(uiB), expC = ExpOf(uiC);
        UInt128 sigA = FracOf(uiA), sigB = FracOf(uiB), sigC = FracOf(uiC);
        bool nanA = expA == MaxBiasedExponent && sigA != 0;
        bool nanB = expB == MaxBiasedExponent && sigB != 0;
        bool nanC = expC == MaxBiasedExponent && sigC != 0;
        if (nanA || nanB || nanC)
        {
            return NaN;
        }
        bool zeroA = expA == 0 && sigA == 0, zeroB = expB == 0 && sigB == 0, zeroC = expC == 0 && sigC == 0;
        if (expA == MaxBiasedExponent || expB == MaxBiasedExponent)
        {
            if (zeroA || zeroB)
            {
                return NaN; // ∞ * 0
            }
            if (expC == MaxBiasedExponent && signC != signP)
            {
                return NaN; // ∞ - ∞
            }
            return InfOf(signP);
        }
        if (expC == MaxBiasedExponent)
        {
            return addend;
        }
        if (zeroA || zeroB)
        {
            if (zeroC)
            {
                return ZeroOf(signP == signC ? signC : rounding == RealRounding.TowardNegative);
            }
            return addend;
        }
        if (expA == 0)
        {
            sigA = NormSubnormal(sigA, out expA);
        }
        if (expB == 0)
        {
            sigB = NormSubnormal(sigB, out expB);
        }
        Mul113(sigA | Hidden, sigB | Hidden, out UInt128 pHi, out UInt128 pLo);
        int sP = expA + expB - 32990; // product = (pHi:pLo) * 2^sP
        if (zeroC)
        {
            return RoundPack256(signP, pHi, pLo, sP, rounding);
        }
        if (expC == 0)
        {
            sigC = NormSubnormal(sigC, out expC);
        }
        sigC |= Hidden;
        // Both terms get two zero low bits so that a jammed (odd) smaller term never lands the
        // difference exactly on a rounding boundary: the true value stays strictly inside (R-1, R+1).
        ShiftLeft256(ref pHi, ref pLo, 2);
        sP -= 2;
        UInt128 cHi = sigC >> 13, cLo = sigC << 115; // sigC * 2^115, aligned with the product's top
        int sC = expC - 16610;
        int d = sC - sP;
        int s;
        if (d >= 0)
        {
            if (d <= 27)
            {
                ShiftLeft256(ref cHi, ref cLo, d);
                s = sP;
            }
            else
            {
                ShiftRightJam256(ref pHi, ref pLo, d);
                s = sC;
            }
        }
        else
        {
            if (d >= -27)
            {
                ShiftLeft256(ref pHi, ref pLo, -d);
                s = sC;
            }
            else
            {
                ShiftRightJam256(ref cHi, ref cLo, -d);
                s = sP;
            }
        }
        UInt128 rHi, rLo;
        bool signZ;
        if (signP == signC)
        {
            Add256(pHi, pLo, cHi, cLo, out rHi, out rLo);
            signZ = signP;
        }
        else
        {
            int c = pHi != cHi ? pHi.CompareTo(cHi) : pLo.CompareTo(cLo);
            if (c == 0)
            {
                return ZeroOf(rounding == RealRounding.TowardNegative);
            }
            if (c > 0)
            {
                Sub256(pHi, pLo, cHi, cLo, out rHi, out rLo);
                signZ = signP;
            }
            else
            {
                Sub256(cHi, cLo, pHi, pLo, out rHi, out rLo);
                signZ = signC;
            }
        }
        return RoundPack256(signZ, rHi, rLo, s, rounding);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void Add256(UInt128 aHi, UInt128 aLo, UInt128 bHi, UInt128 bLo, out UInt128 hi, out UInt128 lo)
    {
        ulong a0 = (ulong)aLo, a1 = (ulong)(aLo >> 64), a2 = (ulong)aHi, a3 = (ulong)(aHi >> 64);
        ulong b0 = (ulong)bLo, b1 = (ulong)(bLo >> 64), b2 = (ulong)bHi, b3 = (ulong)(bHi >> 64);
        ulong r0 = a0 + b0, c = r0 < a0 ? 1UL : 0UL;
        ulong t1 = a1 + b1, r1 = t1 + c;
        c = (t1 < a1 ? 1UL : 0UL) | (r1 < t1 ? 1UL : 0UL);
        ulong t2 = a2 + b2, r2 = t2 + c;
        c = (t2 < a2 ? 1UL : 0UL) | (r2 < t2 ? 1UL : 0UL);
        ulong r3 = a3 + b3 + c;
        lo = new UInt128(r1, r0);
        hi = new UInt128(r3, r2);
    }

    /// <summary>(aHi:aLo) - (bHi:bLo), requiring a &gt;= b.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void Sub256(UInt128 aHi, UInt128 aLo, UInt128 bHi, UInt128 bLo, out UInt128 hi, out UInt128 lo)
    {
        ulong a0 = (ulong)aLo, a1 = (ulong)(aLo >> 64), a2 = (ulong)aHi, a3 = (ulong)(aHi >> 64);
        ulong b0 = (ulong)bLo, b1 = (ulong)(bLo >> 64), b2 = (ulong)bHi, b3 = (ulong)(bHi >> 64);
        ulong r0 = a0 - b0, c = a0 < b0 ? 1UL : 0UL;
        ulong t1 = a1 - b1, r1 = t1 - c;
        c = (a1 < b1 ? 1UL : 0UL) | (t1 < c ? 1UL : 0UL);
        ulong t2 = a2 - b2, r2 = t2 - c;
        c = (a2 < b2 ? 1UL : 0UL) | (t2 < c ? 1UL : 0UL);
        ulong r3 = a3 - b3 - c;
        lo = new UInt128(r1, r0);
        hi = new UInt128(r3, r2);
    }

    private static void ShiftLeft256(ref UInt128 hi, ref UInt128 lo, int n)
    {
        Debug.Assert(n >= 0 && n < 128);
        if (n == 0)
        {
            return;
        }
        hi = (hi << n) | (lo >> (128 - n));
        lo <<= n;
    }

    private static void ShiftRightJam256(ref UInt128 hi, ref UInt128 lo, int n)
    {
        Debug.Assert(n > 0);
        if (n >= 256)
        {
            lo = (hi | lo) != 0 ? UInt128.One : UInt128.Zero;
            hi = UInt128.Zero;
            return;
        }
        if (n >= 128)
        {
            int m = n - 128;
            bool sticky = lo != 0 || (m > 0 && (hi << (128 - m)) != 0);
            lo = m == 0 ? hi : hi >> m;
            hi = UInt128.Zero;
            if (sticky)
            {
                lo |= 1;
            }
            return;
        }
        bool st = (lo << (128 - n)) != 0;
        lo = (lo >> n) | (hi << (128 - n));
        hi >>= n;
        if (st)
        {
            lo |= 1;
        }
    }

    /// <summary>Rounds sign * (hi:lo) * 2^s (nonzero 256-bit magnitude).</summary>
    private static Float128 RoundPack256(bool sign, UInt128 hi, UInt128 lo, int s, RealRounding rounding)
    {
        int lz = hi != 0 ? (int)UInt128.LeadingZeroCount(hi) : 128 + (int)UInt128.LeadingZeroCount(lo);
        Debug.Assert(lz < 256);
        // Normalize so bit 255 is set: top 113 bits = sig, next 64 = extra, the rest is sticky.
        if (lz >= 128)
        {
            hi = lo << (lz - 128);
            lo = UInt128.Zero;
        }
        else if (lz > 0)
        {
            hi = (hi << lz) | (lo >> (128 - lz));
            lo <<= lz;
        }
        UInt128 sig = hi >> 15;
        ulong extra = ((ulong)hi << 49) | (ulong)(lo >> 79);
        if ((lo & ((UInt128.One << 79) - 1)) != 0)
        {
            extra |= 1;
        }
        long exp = 143L - lz + s + 16494;
        exp = Math.Clamp(exp, -100_000, 100_000);
        return RoundPack(sign, (int)exp, sig, extra, rounding);
    }

    // ---------------------------------------------------------------
    // Remainders (always exact)
    // ---------------------------------------------------------------

    /// <summary>IEEE 754 remainder: x - n*y with n = x/y rounded to nearest-even integer (exact).</summary>
    public static Float128 Ieee754Remainder(Float128 left, Float128 right) => Remainder(left, right, ieee: true);

    /// <summary>The IEEE 754 remainder, same as <see cref="Ieee754Remainder(Float128, Float128)"/> (the <see cref="Math.IEEERemainder(double, double)"/> name).</summary>
    public static Float128 IEEERemainder(Float128 x, Float128 y) => Remainder(x, y, ieee: true);

    private static Float128 Remainder(Float128 a, Float128 b, bool ieee)
    {
        UInt128 uiA = a._bits, uiB = b._bits;
        bool signA = SignOf(uiA);
        int expA = ExpOf(uiA), expB = ExpOf(uiB);
        UInt128 sigA = FracOf(uiA), sigB = FracOf(uiB);
        if ((expA == MaxBiasedExponent && sigA != 0) || (expB == MaxBiasedExponent && sigB != 0))
        {
            return NaN;
        }
        if (expA == MaxBiasedExponent || (expB == 0 && sigB == 0))
        {
            return NaN; // ∞ rem y, x rem 0
        }
        if (expB == MaxBiasedExponent || (expA == 0 && sigA == 0))
        {
            return a; // x rem ∞ = x; ±0 rem y = ±0
        }
        if (expA == 0)
        {
            sigA = NormSubnormal(sigA, out expA);
        }
        if (expB == 0)
        {
            sigB = NormSubnormal(sigB, out expB);
        }
        sigA |= Hidden;
        sigB |= Hidden;
        // |a| = sigA * 2^(expA - 16495), |b| = sigB * 2^(expB - 16495)
        if (expA < expB)
        {
            if (!ieee || expA < expB - 1)
            {
                return a; // |a| < |b| (and |a| < |b|/2 for the IEEE remainder)
            }
            // expA == expB - 1: compare |a| with |b|/2 = sigB * 2^(expA - 16495).
            if (sigA <= sigB)
            {
                return a;
            }
            return NormRoundPack(!signA, expA - 1, (sigB << 1) - sigA, RealRounding.NearestEven);
        }
        UInt128 r = sigA;
        ulong lastQ = 0;
        if (r >= sigB)
        {
            r -= sigB;
            lastQ = 1;
        }
        int diff = expA - expB;
        if (diff > 0)
        {
            UInt128 dn = sigB << 15;
            ulong d1 = (ulong)(dn >> 64), d0 = (ulong)dn;
            while (diff > 0)
            {
                int k = Math.Min(diff, 64);
                UInt128 rr = r << 15; // < dn
                ulong u2, u1, u0;
                if (k == 64)
                {
                    u2 = (ulong)(rr >> 64);
                    u1 = (ulong)rr;
                    u0 = 0;
                }
                else
                {
                    u2 = (ulong)(rr >> (128 - k));
                    UInt128 low = rr << k;
                    u1 = (ulong)(low >> 64);
                    u0 = (ulong)low;
                }
                lastQ = Div3By2(u2, u1, u0, d1, d0, out UInt128 rem);
                r = rem >> 15;
                diff -= k;
            }
        }
        bool sign = signA;
        if (ieee)
        {
            UInt128 twice = r << 1;
            if (twice > sigB || (twice == sigB && (lastQ & 1) != 0))
            {
                r = sigB - r;
                sign = !sign;
            }
        }
        if (r == 0)
        {
            return ZeroOf(signA);
        }
        return NormRoundPack(sign, expB - 1, r, RealRounding.NearestEven);
    }

    // ---------------------------------------------------------------
    // Rounding to integers
    // ---------------------------------------------------------------

    private enum IntRounding
    {
        HalfEven,
        HalfAway,
        Truncate,
        Floor,
        Ceiling,
    }

    private static Float128 RoundToInt(Float128 x, IntRounding mode)
    {
        UInt128 ui = x._bits;
        int exp = ExpOf(ui);
        if (exp >= 0x3FFF + 112)
        {
            return exp == MaxBiasedExponent && FracOf(ui) != 0 ? NaN : x;
        }
        bool sign = SignOf(ui);
        if (exp < 0x3FFF)
        {
            if ((ui & ~SignMask) == 0)
            {
                return x;
            }
            bool one = mode switch
            {
                IntRounding.HalfEven => exp == 0x3FFE && FracOf(ui) != 0,
                IntRounding.HalfAway => exp == 0x3FFE,
                IntRounding.Truncate => false,
                IntRounding.Floor => sign,
                _ => !sign,
            };
            return one ? new(One._bits | SignBits(sign)) : ZeroOf(sign);
        }
        int fracBits = 0x3FFF + 112 - exp; // 1..112
        UInt128 lastBit = UInt128.One << fracBits;
        UInt128 roundMask = lastBit - 1;
        if ((ui & roundMask) == 0)
        {
            return x;
        }
        switch (mode)
        {
            case IntRounding.HalfEven:
                ui += lastBit >> 1;
                if ((ui & roundMask) == 0)
                {
                    ui &= ~lastBit; // exact tie
                }
                break;
            case IntRounding.HalfAway:
                ui += lastBit >> 1;
                break;
            case IntRounding.Floor:
                if (sign)
                {
                    ui += roundMask;
                }
                break;
            case IntRounding.Ceiling:
                if (!sign)
                {
                    ui += roundMask;
                }
                break;
        }
        return new(ui & ~roundMask);
    }

    // ---------------------------------------------------------------
    // Exponent manipulation and neighbors
    // ---------------------------------------------------------------

    /// <summary>x * 2^n, correctly rounded (nearest-even) on underflow/overflow.</summary>
    public static Float128 ScaleB(Float128 x, int n) => ScaleB(x, n, RealRounding.NearestEven);

    /// <summary>x * 2^n, correctly rounded with the given mode on underflow/overflow.</summary>
    public static Float128 ScaleB(Float128 x, int n, RealRounding rounding)
    {
        UInt128 ui = x._bits;
        int exp = ExpOf(ui);
        UInt128 sig = FracOf(ui);
        if (exp == MaxBiasedExponent)
        {
            return sig != 0 ? NaN : x;
        }
        if (exp == 0)
        {
            if (sig == 0)
            {
                return x;
            }
            sig = NormSubnormal(sig, out exp);
        }
        int e = (int)Math.Clamp((long)exp + n, -100_000, 100_000);
        return RoundPack(SignOf(ui), e - 1, sig | Hidden, 0, rounding);
    }

    /// <summary>floor(log2|x|); <see cref="int.MinValue"/> for zero, <see cref="int.MaxValue"/> for NaN and infinity.</summary>
    public static int ILogB(Float128 x)
    {
        int exp = x.BiasedExponent;
        if (exp == MaxBiasedExponent)
        {
            return int.MaxValue;
        }
        if (exp == 0)
        {
            UInt128 f = x.Fraction;
            if (f == 0)
            {
                return int.MinValue;
            }
            NormSubnormal(f, out exp);
        }
        return exp - ExponentBias;
    }

    /// <summary>The next value toward +∞ (like <see cref="Math.BitIncrement(double)"/>).</summary>
    public static Float128 BitIncrement(Float128 x)
    {
        UInt128 b = x._bits;
        if (x.BiasedExponent == MaxBiasedExponent)
        {
            return b == (InfBits | SignMask) ? MinValue : x;
        }
        if (b == SignMask)
        {
            return Epsilon;
        }
        return new(SignOf(b) ? b - 1 : b + 1);
    }

    /// <summary>The next value toward -∞ (like <see cref="Math.BitDecrement(double)"/>).</summary>
    public static Float128 BitDecrement(Float128 x)
    {
        UInt128 b = x._bits;
        if (x.BiasedExponent == MaxBiasedExponent)
        {
            return b == InfBits ? MaxValue : x;
        }
        if (b == 0)
        {
            return -Epsilon;
        }
        return new(SignOf(b) ? b + 1 : b - 1);
    }
}
