using System.Globalization;
using System.Numerics;
using Xunit;

namespace Sdcb.Arithmetic.Tests;

public class Float128Tests
{
    private static readonly RealRounding[] Modes =
        [RealRounding.NearestEven, RealRounding.TowardZero, RealRounding.TowardPositive, RealRounding.TowardNegative, RealRounding.AwayFromZero];

    private static UInt128 Hex(ulong hi, ulong lo) => new(hi, lo);
    private static readonly UInt128 Sign = UInt128.One << 127;
    private static readonly UInt128 PosInf = Hex(0x7FFF_0000_0000_0000, 0);
    private static readonly UInt128 QNaN = Hex(0x7FFF_8000_0000_0000, 0);
    private static readonly UInt128 MaxFinite = Hex(0x7FFE_FFFF_FFFF_FFFF, ulong.MaxValue);

    private static Float128 F(UInt128 bits) => Float128.FromBits(bits);
    private static Real R(Float128 x) => x.ToReal();
    private static string H(UInt128 b) => $"0x{(ulong)(b >> 64):X16}_{(ulong)b:X16}";

    // ---------------------------------------------------------------
    // Oracle: Real evaluated toward zero and away from zero at 160 bits. Equal results are exact;
    // otherwise the toward-zero result with its last bit forced to 1 is the round-to-odd value, whose
    // rounding to binary128 (<= 113 bits, subnormals included) is the correctly rounded result.
    // ---------------------------------------------------------------

    private const int OracleBits = 160;

    private static UInt128 Oracle(Func<int, RealRounding, Real> f, RealRounding mode)
    {
        Real z = f(OracleBits, RealRounding.TowardZero);
        if (z.IsZero)
        {
            return f(OracleBits, mode).ToBinary128Bits(mode); // exact zero: its sign depends on the mode
        }
        if (!z.IsNumber)
        {
            return z.ToBinary128Bits(mode);
        }
        Real a = f(OracleBits, RealRounding.AwayFromZero);
        if (a.CompareTo(z) == 0)
        {
            return z.ToBinary128Bits(mode);
        }
        BigInt m = z.Mantissa | BigInt.One;
        Real odd = Real.FromBigInt(z.SignBit ? -m : m, z.Exponent, OracleBits);
        return odd.ToBinary128Bits(mode);
    }

    private static UInt128 Canon(UInt128 b) => Float128.IsNaN(F(b)) ? QNaN : b;

    private static void AssertBits(UInt128 expected, Float128 actual, string context)
    {
        if (Canon(expected) != actual.ToBits())
        {
            Assert.Fail($"{context}: expected {H(expected)}, got {H(actual.ToBits())}");
        }
    }

    // ---------------------------------------------------------------
    // Random operands biased toward the hard cases
    // ---------------------------------------------------------------

    private static UInt128 RandomFrac(Random rng)
    {
        ulong hi = (ulong)rng.NextInt64() ^ ((ulong)rng.Next() << 33);
        ulong lo = (ulong)rng.NextInt64() ^ ((ulong)rng.Next() << 63);
        UInt128 f = new UInt128(hi, lo) & ((UInt128.One << 112) - 1);
        switch (rng.Next(8))
        {
            case 0: // sparse: a few bits
                f = UInt128.Zero;
                for (int i = rng.Next(4); i >= 0; i--)
                {
                    f |= UInt128.One << rng.Next(112);
                }
                break;
            case 1: // runs of ones from the top or bottom
                f = rng.Next(2) == 0 ? ((UInt128.One << 112) - 1) >> rng.Next(112) : ((UInt128.One << 112) - 1) ^ ((UInt128.One << rng.Next(112)) - 1);
                break;
            case 2:
                f = UInt128.Zero;
                break;
        }
        return f;
    }

    private static UInt128 RandomBits(Random rng)
    {
        int exp = rng.Next(12) switch
        {
            0 => 0,
            1 => 1 + rng.Next(3),
            2 => 0x7FFE - rng.Next(3),
            3 => rng.Next(1, 0x7FFF),
            4 => rng.Next(1, 120),
            5 => 0x7FFE - rng.Next(120),
            _ => 0x3FFF + rng.Next(-60, 60),
        };
        UInt128 b = ((UInt128)(uint)exp << 112) | RandomFrac(rng);
        if (rng.Next(2) == 0)
        {
            b |= Sign;
        }
        return b;
    }

    /// <summary>A second operand related to the first: close exponent, close bits, or independent.</summary>
    private static UInt128 RelatedBits(Random rng, UInt128 a)
    {
        switch (rng.Next(6))
        {
            case 0: // close bit pattern: cancellation / near ties
                UInt128 delta = (UInt128)(ulong)rng.Next(0, 1 << rng.Next(1, 30));
                UInt128 b = rng.Next(2) == 0 ? a + delta : a - delta;
                return rng.Next(2) == 0 ? b ^ Sign : b;
            case 1: // same exponent, other fraction
            {
                UInt128 e = a & ~((UInt128.One << 112) - 1);
                return e | RandomFrac(rng);
            }
            case 2: // exponent difference 100..130 (alignment edge)
            {
                int exp = (int)((a >> 112) & 0x7FFF) - rng.Next(100, 130);
                if (exp < 0)
                {
                    exp = 0;
                }
                return ((UInt128)(uint)exp << 112) | RandomFrac(rng) | (rng.Next(2) == 0 ? Sign : UInt128.Zero);
            }
            default:
                return RandomBits(rng);
        }
    }

    private static UInt128 SpecialOrRandom(Random rng, UInt128[] specials) =>
        rng.Next(20) == 0 ? specials[rng.Next(specials.Length)] : RandomBits(rng);

    private static readonly UInt128[] Specials =
    [
        UInt128.Zero, Sign, UInt128.One, Sign | 1, Hex(0x0001_0000_0000_0000, 0), Hex(0x0000_FFFF_FFFF_FFFF, ulong.MaxValue),
        MaxFinite, MaxFinite | Sign, PosInf, PosInf | Sign, QNaN, Hex(0x3FFF_0000_0000_0000, 0), Hex(0xBFFF_0000_0000_0000, 0),
        Hex(0x7FFF_0000_0000_0000, 1), Hex(0x4000_0000_0000_0000, 0), Hex(0x3FFE_FFFF_FFFF_FFFF, ulong.MaxValue),
    ];

    // ---------------------------------------------------------------
    // Bit-exact basic operations against the Real oracle
    // ---------------------------------------------------------------

    private static void CheckBinary(int seed, int count, Func<Float128, Float128, RealRounding, Float128> op,
        Func<Real, Real, int, RealRounding, Real> oracle, string name)
    {
        Random rng = new(seed);
        for (int i = 0; i < count; i++)
        {
            UInt128 a = SpecialOrRandom(rng, Specials);
            UInt128 b = rng.Next(3) == 0 ? SpecialOrRandom(rng, Specials) : RelatedBits(rng, a);
            Real ra = R(F(a)), rb = R(F(b));
            foreach (RealRounding mode in Modes)
            {
                UInt128 expected = Oracle((p, m) => oracle(ra, rb, p, m), mode);
                AssertBits(expected, op(F(a), F(b), mode), $"{name}({H(a)}, {H(b)}) {mode}");
            }
        }
    }

    [Fact]
    public void AddMatchesOracle() => CheckBinary(1, 6000, Float128.Add, Real.Add, "Add");

    [Fact]
    public void SubtractMatchesOracle() => CheckBinary(2, 6000, Float128.Subtract, Real.Subtract, "Sub");

    [Fact]
    public void MultiplyMatchesOracle() => CheckBinary(3, 6000, Float128.Multiply, Real.Multiply, "Mul");

    [Fact]
    public void DivideMatchesOracle() => CheckBinary(4, 6000, Float128.Divide, Real.Divide, "Div");

    [Fact]
    public void SqrtMatchesOracle()
    {
        Random rng = new(5);
        for (int i = 0; i < 6000; i++)
        {
            UInt128 a = SpecialOrRandom(rng, Specials);
            if (i % 3 == 0)
            {
                a &= ~Sign;
            }
            Real ra = R(F(a));
            foreach (RealRounding mode in Modes)
            {
                UInt128 expected = Oracle((p, m) => Real.Sqrt(ra, p, m), mode);
                AssertBits(expected, Float128.Sqrt(F(a), mode), $"Sqrt({H(a)}) {mode}");
            }
        }
        // Exact squares.
        for (int i = 0; i < 2000; i++)
        {
            Float128 x = Float128.FromBits(RandomBits(rng) & ~Sign);
            Float128 t = Float128.Truncate(Float128.ScaleB(Float128.Abs(x), -Float128.ILogB(x) + 50));
            Float128 sq = t * t;
            Assert.Equal(t.ToBits(), Float128.Sqrt(sq).ToBits());
        }
    }

    [Fact]
    public void FusedMultiplyAddMatchesOracle()
    {
        Random rng = new(6);
        for (int i = 0; i < 6000; i++)
        {
            UInt128 a = SpecialOrRandom(rng, Specials);
            UInt128 b = SpecialOrRandom(rng, Specials);
            UInt128 c;
            int kind = rng.Next(4);
            if (kind == 0)
            {
                // c close to -(a*b): massive cancellation.
                Float128 p = Float128.Multiply(F(a), F(b), RealRounding.TowardZero);
                c = (p.ToBits() ^ Sign) + (UInt128)(ulong)rng.Next(-3, 4);
            }
            else
            {
                c = kind == 1 ? RelatedBits(rng, a) : SpecialOrRandom(rng, Specials);
            }
            Real ra = R(F(a)), rb = R(F(b)), rc = R(F(c));
            foreach (RealRounding mode in Modes)
            {
                UInt128 expected = Oracle((p, m) => Real.FusedMultiplyAdd(ra, rb, rc, p, m), mode);
                AssertBits(expected, Float128.FusedMultiplyAdd(F(a), F(b), F(c), mode), $"Fma({H(a)}, {H(b)}, {H(c)}) {mode}");
            }
        }
    }

    // ---------------------------------------------------------------
    // Boundaries of the normal-operand fast paths (alignment shifts, normal/subnormal threshold,
    // overflow threshold, exact ties)
    // ---------------------------------------------------------------

    private static UInt128 Make(bool negative, int exp, UInt128 frac) =>
        ((UInt128)(uint)exp << 112) | frac | (negative ? Sign : UInt128.Zero);

    private static void CheckArith(UInt128 a, UInt128 b, bool add, bool mul, bool div)
    {
        Real ra = R(F(a)), rb = R(F(b));
        foreach (RealRounding mode in Modes)
        {
            if (add)
            {
                AssertBits(Oracle((p, m) => Real.Add(ra, rb, p, m), mode), Float128.Add(F(a), F(b), mode), $"Add({H(a)}, {H(b)}) {mode}");
                AssertBits(Oracle((p, m) => Real.Subtract(ra, rb, p, m), mode), Float128.Subtract(F(a), F(b), mode), $"Sub({H(a)}, {H(b)}) {mode}");
            }
            if (mul)
            {
                AssertBits(Oracle((p, m) => Real.Multiply(ra, rb, p, m), mode), Float128.Multiply(F(a), F(b), mode), $"Mul({H(a)}, {H(b)}) {mode}");
            }
            if (div)
            {
                AssertBits(Oracle((p, m) => Real.Divide(ra, rb, p, m), mode), Float128.Divide(F(a), F(b), mode), $"Div({H(a)}, {H(b)}) {mode}");
            }
        }
    }

    [Fact]
    public void AddAlignmentBoundaries()
    {
        Random rng = new(11);
        int[] shifts = [0, 1, 2, 3, 62, 63, 64, 65, 66, 111, 112, 113, 114, 115, 116, 126, 127, 128, 129, 130, 191, 192, 193, 400];
        foreach (int d in shifts)
        {
            for (int i = 0; i < 60; i++)
            {
                int expA = rng.Next(4) switch
                {
                    0 => d + 1 + rng.Next(3), // difference lands at or just above the subnormal threshold
                    1 => 0x7FFE - rng.Next(3), // overflow threshold
                    _ => 0x3FFF + rng.Next(-40, 40),
                };
                if (expA - d < 1 || expA > 0x7FFE)
                {
                    continue;
                }
                UInt128 a = Make(rng.Next(2) == 0, expA, RandomFrac(rng));
                UInt128 b = Make(rng.Next(2) == 0, expA - d, RandomFrac(rng));
                CheckArith(a, b, add: true, mul: false, div: false);
                CheckArith(b, a, add: true, mul: false, div: false);
            }
        }
    }

    [Fact]
    public void ProductAndQuotientRangeBoundaries()
    {
        Random rng = new(12);
        // Target biased result exponents around the subnormal threshold, deep underflow and overflow.
        int[] targets = [-116, -115, -114, -113, -112, -2, -1, 0, 1, 2, 3, 0x7FFC, 0x7FFD, 0x7FFE, 0x7FFF, 0x8000];
        foreach (int e in targets)
        {
            for (int i = 0; i < 80; i++)
            {
                int expA = rng.Next(2) == 0 ? rng.Next(1, 0x7FFF) : (e < 0x4000 ? rng.Next(1, 200) : rng.Next(0x7F00, 0x7FFF));
                int expMul = e - expA + 0x3FFF + rng.Next(-1, 2);
                int expDiv = expA - e + 0x3FFF + rng.Next(-1, 2);
                UInt128 fa = RandomFrac(rng), fb = RandomFrac(rng);
                bool sa = rng.Next(2) == 0, sb = rng.Next(2) == 0;
                if (expMul >= 1 && expMul <= 0x7FFE)
                {
                    CheckArith(Make(sa, expA, fa), Make(sb, expMul, fb), add: false, mul: true, div: false);
                }
                if (expDiv >= 1 && expDiv <= 0x7FFE)
                {
                    CheckArith(Make(sa, expA, fa), Make(sb, expDiv, fb), add: false, mul: false, div: true);
                }
            }
        }
    }

    [Fact]
    public void ExactTiesRoundToEven()
    {
        Random rng = new(13);
        UInt128 fracMask = (UInt128.One << 112) - 1;
        for (int i = 0; i < 400; i++)
        {
            // sigA odd and < 2^114 / 3: sigA * 3 has 114 bits ending in 1, an exact tie at 113 bits.
            UInt128 frac = (new UInt128((ulong)rng.NextInt64(), (ulong)rng.NextInt64()) & fracMask) >> 2 | 1;
            int expA = rng.Next(3) switch
            {
                0 => rng.Next(1, 120), // the product may land in the subnormal range
                1 => 0x7FFE - rng.Next(3),
                _ => 0x3FFF + rng.Next(-40, 40),
            };
            UInt128 a = Make(rng.Next(2) == 0, expA, frac);
            UInt128 b = Make(rng.Next(2) == 0, 0x3FFF + rng.Next(-3, 2), UInt128.One << 111); // ±1.5 * 2^k
            CheckArith(a, b, add: false, mul: true, div: false);
            // a + b with b just below a's last place: exact half-ulp sums and differences.
            UInt128 c = Make(rng.Next(2) == 0, expA > 113 ? expA - 113 : 1, rng.Next(2) == 0 ? UInt128.Zero : RandomFrac(rng));
            CheckArith(a, c, add: true, mul: false, div: false);
        }
    }

    [Fact]
    public void FusedMultiplyAddAlignmentBoundaries()
    {
        Random rng = new(14);
        int[] offsets = [-300, -257, -256, -255, -200, -129, -128, -127, -65, -64, -63, -29, -28, -27, -26, -2, -1, 0, 1, 2, 26, 27, 28, 29, 63, 64, 65, 127, 128, 129, 200, 255, 256, 257, 300];
        foreach (int d in offsets)
        {
            for (int i = 0; i < 40; i++)
            {
                int expA = rng.Next(3) switch
                {
                    0 => rng.Next(1, 200),
                    1 => rng.Next(0x7F00, 0x7FFF),
                    _ => 0x3FFF + rng.Next(-30, 30),
                };
                int expB = rng.Next(2) == 0 ? 0x3FFF + rng.Next(-30, 30) : rng.Next(1, 0x7FFF);
                // C's alignment offset relative to the product in the fast path: expC - expA - expB + 16382.
                int expC = d + expA + expB - 16382;
                if (expC < 1 || expC > 0x7FFE)
                {
                    continue;
                }
                UInt128 a = Make(rng.Next(2) == 0, expA, RandomFrac(rng));
                UInt128 b = Make(rng.Next(2) == 0, expB, RandomFrac(rng));
                UInt128 c = Make(rng.Next(2) == 0, expC, RandomFrac(rng));
                Real ra = R(F(a)), rb = R(F(b)), rc = R(F(c));
                foreach (RealRounding mode in Modes)
                {
                    UInt128 expected = Oracle((p, m) => Real.FusedMultiplyAdd(ra, rb, rc, p, m), mode);
                    AssertBits(expected, Float128.FusedMultiplyAdd(F(a), F(b), F(c), mode), $"Fma({H(a)}, {H(b)}, {H(c)}) {mode}");
                }
            }
        }
    }

    [Fact]
    public void SqrtRangeBoundaries()
    {
        Random rng = new(15);
        int[] exps = [1, 2, 3, 4, 0x3FFE, 0x3FFF, 0x4000, 0x7FFC, 0x7FFD, 0x7FFE];
        foreach (int e in exps)
        {
            for (int i = 0; i < 60; i++)
            {
                UInt128 a = Make(false, e, i == 0 ? UInt128.Zero : i == 1 ? (UInt128.One << 112) - 1 : RandomFrac(rng));
                Real ra = R(F(a));
                foreach (RealRounding mode in Modes)
                {
                    AssertBits(Oracle((p, m) => Real.Sqrt(ra, p, m), mode), Float128.Sqrt(F(a), mode), $"Sqrt({H(a)}) {mode}");
                }
            }
        }
    }

    [Fact]
    public void SpecialOperandPairs()
    {
        foreach (UInt128 a in Specials)
        {
            foreach (UInt128 b in Specials)
            {
                Real ra = R(F(a)), rb = R(F(b));
                foreach (RealRounding mode in Modes)
                {
                    AssertBits(Oracle((p, m) => Real.Add(ra, rb, p, m), mode), Float128.Add(F(a), F(b), mode), $"Add {H(a)} {H(b)} {mode}");
                    AssertBits(Oracle((p, m) => Real.Subtract(ra, rb, p, m), mode), Float128.Subtract(F(a), F(b), mode), $"Sub {H(a)} {H(b)} {mode}");
                    AssertBits(Oracle((p, m) => Real.Multiply(ra, rb, p, m), mode), Float128.Multiply(F(a), F(b), mode), $"Mul {H(a)} {H(b)} {mode}");
                    AssertBits(Oracle((p, m) => Real.Divide(ra, rb, p, m), mode), Float128.Divide(F(a), F(b), mode), $"Div {H(a)} {H(b)} {mode}");
                    foreach (UInt128 c in Specials)
                    {
                        Real rc = R(F(c));
                        AssertBits(Oracle((p, m) => Real.FusedMultiplyAdd(ra, rb, rc, p, m), mode),
                            Float128.FusedMultiplyAdd(F(a), F(b), F(c), mode), $"Fma {H(a)} {H(b)} {H(c)} {mode}");
                    }
                }
            }
        }
    }

    [Fact]
    public void SignedZerosAndInfinities()
    {
        Float128 pz = Float128.Zero, nz = Float128.NegativeZero, inf = Float128.PositiveInfinity;
        Float128 x = (Float128)1.25;
        Assert.Equal(Sign, (nz + nz).ToBits());
        Assert.Equal(UInt128.Zero, (pz + nz).ToBits());
        Assert.Equal(UInt128.Zero, (x - x).ToBits());
        Assert.Equal(Sign, Float128.Subtract(x, x, RealRounding.TowardNegative).ToBits());
        Assert.Equal(Sign, Float128.Add(x, -x, RealRounding.TowardNegative).ToBits());
        Assert.Equal(Sign, (nz * x).ToBits());
        Assert.Equal(Sign, Float128.Sqrt(nz).ToBits());
        Assert.True(Float128.IsNaN(inf - inf));
        Assert.True(Float128.IsNaN(inf * pz));
        Assert.True(Float128.IsNaN(pz / pz));
        Assert.True(Float128.IsNaN(Float128.Sqrt(-x)));
        Assert.Equal(QNaN, (Float128.FromBits(Hex(0xFFFF_0000_0000_0000, 7)) + x).ToBits()); // canonical NaN
        Assert.Equal(PosInf | Sign, (x / nz).ToBits());
        Assert.Equal(PosInf, (inf + x).ToBits());
        Assert.Equal(UInt128.Zero, (x / inf).ToBits());
    }

    [Fact]
    public void OverflowAndUnderflowPerMode()
    {
        Float128 max = Float128.MaxValue;
        foreach (bool neg in new[] { false, true })
        {
            Float128 m = neg ? -max : max;
            foreach (RealRounding mode in Modes)
            {
                bool toInf = mode switch
                {
                    RealRounding.TowardZero => false,
                    RealRounding.TowardPositive => !neg,
                    RealRounding.TowardNegative => neg,
                    _ => true,
                };
                UInt128 expected = (toInf ? PosInf : MaxFinite) | (neg ? Sign : UInt128.Zero);
                Assert.Equal(expected, Float128.Multiply(m, (Float128)2, mode).ToBits());
                Assert.Equal(expected, Float128.Add(m, m, mode).ToBits());
            }
        }
        Float128 eps = Float128.Epsilon;
        Assert.Equal(UInt128.Zero, (eps * (Float128)0.5).ToBits()); // tie to even (0)
        Assert.Equal((UInt128)2, (eps * (Float128)1.5).ToBits()); // tie to even (2)
        Assert.Equal(UInt128.One, Float128.Multiply(eps, (Float128)0.25, RealRounding.TowardPositive).ToBits());
        Assert.Equal(Sign, Float128.Multiply(-eps, (Float128)0.25, RealRounding.TowardPositive).ToBits());
        Assert.Equal(Sign | 1, Float128.Multiply(-eps, (Float128)0.25, RealRounding.AwayFromZero).ToBits());
        Assert.Equal(UInt128.Zero, Float128.Multiply(eps, eps, RealRounding.TowardZero).ToBits());
        Assert.Equal(UInt128.One, Float128.Multiply(eps, eps, RealRounding.TowardPositive).ToBits());
        // min normal / 2 is exactly the largest-ish subnormal power of two
        Assert.Equal(Hex(0x0000_8000_0000_0000, 0), (Float128.MinNormal / (Float128)2).ToBits());
    }

    [Fact]
    public void SqrtSpotCheck()
    {
        Assert.Equal(Real.Sqrt(2, 113).ToBinary128Bits(), Float128.Sqrt(2).ToBits());
        Assert.Equal(Real.Sqrt(3, 113).ToBinary128Bits(), Float128.Sqrt(3).ToBits());
        Assert.Equal((Float128)3, Float128.Sqrt(9));
    }

    // ---------------------------------------------------------------
    // Remainders and rounding to integers (exact)
    // ---------------------------------------------------------------

    [Fact]
    public void RemaindersAreExact()
    {
        Random rng = new(7);
        for (int i = 0; i < 4000; i++)
        {
            UInt128 a = SpecialOrRandom(rng, Specials);
            UInt128 b = rng.Next(2) == 0 ? RelatedBits(rng, a) : SpecialOrRandom(rng, Specials);
            if (i % 4 == 0)
            {
                // moderate exponent difference keeps the oracle cheap but covers multi-chunk loops
                b = (b & ~(UInt128.MaxValue << 112 & ~Sign)) | ((UInt128)(uint)Math.Max(1, (int)((a >> 112) & 0x7FFF) - rng.Next(0, 400)) << 112);
            }
            Float128 fa = F(a), fb = F(b);
            Real ra = R(fa), rb = R(fb);
            int ea = Float128.ILogB(fa), eb = Float128.ILogB(fb);
            if (Float128.IsFinite(fa) && Float128.IsFinite(fb) && !Float128.IsZero(fa) && !Float128.IsZero(fb) && ea - eb > 2000)
            {
                continue; // keep Real's exact quotient small
            }
            UInt128 fmod = Oracle((p, m) => Real.TruncatedRemainder(ra, rb, p, m), RealRounding.NearestEven);
            UInt128 ieee = Oracle((p, m) => Real.Ieee754Remainder(ra, rb, p, m), RealRounding.NearestEven);
            AssertBits(fmod, fa % fb, $"fmod({H(a)}, {H(b)})");
            AssertBits(ieee, Float128.Ieee754Remainder(fa, fb), $"remainder({H(a)}, {H(b)})");
        }
        // Huge exponent difference: 2^16000 mod 3 = (2^2)^8000 mod 3 = 1.
        Assert.Equal((Float128)1, Float128.ScaleB(1, 16000) % 3);
        Assert.Equal((Float128)1, Float128.Ieee754Remainder(Float128.ScaleB(1, 16000), 3));
        Assert.Equal((Float128)(-1), Float128.Ieee754Remainder(Float128.ScaleB(1, 16001), 3)); // 2 - 3
        Assert.Equal((Float128)(-1), Float128.Ieee754Remainder(5, 3));
        Assert.Equal((Float128)1, Float128.Ieee754Remainder(5, 4));
        Assert.Equal((Float128)(-2), Float128.Ieee754Remainder(6, 4)); // quotient 1.5 rounds to 2
    }

    [Fact]
    public void RoundingToIntegers()
    {
        Random rng = new(8);
        for (int i = 0; i < 4000; i++)
        {
            UInt128 a = RandomBits(rng);
            if (i % 2 == 0)
            {
                a = (a & ~(UInt128.MaxValue << 112 & ~Sign)) | ((UInt128)(uint)(0x3FFF + rng.Next(-3, 116)) << 112);
            }
            Float128 x = F(a);
            Real r = R(x);
            CheckInt(r.Floor(), Float128.Floor(x), x, "Floor");
            CheckInt(r.Ceiling(), Float128.Ceiling(x), x, "Ceiling");
            CheckInt(r.Truncate(), Float128.Truncate(x), x, "Truncate");
            CheckInt(Real.Round(r), Float128.Round(x), x, "Round");
            CheckInt(Real.Round(r, MidpointRounding.AwayFromZero), Float128.Round(x, MidpointRounding.AwayFromZero), x, "RoundAway");
        }
        Assert.Equal((Float128)2, Float128.Round((Float128)2.5));
        Assert.Equal((Float128)3, Float128.Round((Float128)2.5, MidpointRounding.AwayFromZero));
        Assert.Equal(Sign, Float128.Round((Float128)(-0.5)).ToBits());
        Assert.Equal(Sign, Float128.Ceiling((Float128)(-0.5)).ToBits());
        Assert.Equal((Float128)(-1), Float128.Floor((Float128)(-0.5)));
        Assert.Equal(Float128.Parse("1.23"), Float128.Round(Float128.Parse("1.2345"), 2));
        Assert.Equal(Float128.Parse("1.24"), Float128.Round(Float128.Parse("1.235"), 2, MidpointRounding.AwayFromZero));
    }

    private static void CheckInt(Real expected, Float128 actual, Float128 x, string name)
    {
        UInt128 e = expected.ToBinary128Bits();
        if (expected.IsZero)
        {
            e = Float128.IsNegative(x) ? Sign : UInt128.Zero;
        }
        AssertBits(e, actual, $"{name}({H(x.ToBits())})");
    }

    // ---------------------------------------------------------------
    // Conversions
    // ---------------------------------------------------------------

    [Theory]
    [InlineData(1.0, 0x3FFF_0000_0000_0000UL, 0UL)]
    [InlineData(-2.0, 0xC000_0000_0000_0000UL, 0UL)]
    [InlineData(0.5, 0x3FFE_0000_0000_0000UL, 0UL)]
    [InlineData(-0.0, 0x8000_0000_0000_0000UL, 0UL)]
    [InlineData(0.0, 0UL, 0UL)]
    [InlineData(double.PositiveInfinity, 0x7FFF_0000_0000_0000UL, 0UL)]
    [InlineData(double.NegativeInfinity, 0xFFFF_0000_0000_0000UL, 0UL)]
    [InlineData(0.1, 0x3FFB_9999_9999_9999UL, 0xA000_0000_0000_0000UL)]
    public void KnownEncodings(double value, ulong hi, ulong lo)
    {
        Float128 f = value;
        Assert.Equal(Hex(hi, lo), f.ToBits());
        Assert.Equal(BitConverter.DoubleToUInt64Bits(value), BitConverter.DoubleToUInt64Bits((double)f));
    }

    [Fact]
    public void Constants()
    {
        Assert.Equal(Real.Pi(300).ToBinary128Bits(), Float128.Pi.ToBits());
        Assert.Equal(Real.E(300).ToBinary128Bits(), Float128.E.ToBits());
        Assert.Equal(Real.Pi(300).ShiftLeft(1).ToBinary128Bits(), Float128.Tau.ToBits());
        Assert.Equal(Hex(0x3FFD_5555_5555_5555, 0x5555_5555_5555_5555), (Float128.One / 3).ToBits());
        Assert.Equal(MaxFinite, Float128.MaxValue.ToBits());
        Assert.Equal(MaxFinite | Sign, Float128.MinValue.ToBits());
        Assert.Equal(UInt128.One, Float128.Epsilon.ToBits());
        Assert.Equal(Hex(0x0001_0000_0000_0000, 0), Float128.MinNormal.ToBits());
        Assert.Equal(QNaN, Float128.NaN.ToBits());
        Assert.True(Float128.IsSubnormal(Float128.Epsilon));
        Assert.True(Float128.IsNormal(Float128.MinNormal));
        Assert.Equal(-16494, Float128.ILogB(Float128.Epsilon));
        Assert.Equal(16383, Float128.ILogB(Float128.MaxValue));
        Assert.Equal(Float128.PositiveInfinity, Float128.BitIncrement(Float128.MaxValue));
        Assert.Equal(Float128.Epsilon, Float128.BitIncrement(Float128.NegativeZero));
        Assert.Equal(-Float128.Epsilon, Float128.BitDecrement(Float128.Zero));
        Assert.Equal(Float128.MaxValue, Float128.BitDecrement(Float128.PositiveInfinity));
    }

    [Fact]
    public void DoubleConversionsRoundTripAndRound()
    {
        Random rng = new(9);
        for (int i = 0; i < 20000; i++)
        {
            double d = BitConverter.Int64BitsToDouble(rng.NextInt64() ^ ((long)rng.Next() << 62));
            if (i % 5 == 0)
            {
                d = BitConverter.Int64BitsToDouble((long)((ulong)rng.NextInt64() & 0x800F_FFFF_FFFF_FFFF)); // subnormals
            }
            Float128 f = d;
            if (double.IsNaN(d))
            {
                Assert.True(Float128.IsNaN(f));
                continue;
            }
            Assert.Equal(Real.FromDouble(d, 53).ToBinary128Bits(), f.ToBits());
            Assert.Equal(BitConverter.DoubleToUInt64Bits(d), BitConverter.DoubleToUInt64Bits((double)f));
        }
        for (int i = 0; i < 20000; i++)
        {
            UInt128 b = RandomBits(rng);
            if (i % 3 == 0)
            {
                // exponents around the double/float/half ranges
                int[] centers = [0x3FFF - 1022, 0x3FFF - 1074, 0x3FFF + 1023, 0x3FFF - 126, 0x3FFF - 149, 0x3FFF + 127, 0x3FFF - 14, 0x3FFF - 24, 0x3FFF + 15];
                int e = centers[rng.Next(centers.Length)] + rng.Next(-3, 4);
                b = (b & ~(UInt128.MaxValue << 112 & ~Sign)) | ((UInt128)(uint)e << 112);
            }
            Float128 f = F(b);
            Real r = R(f);
            foreach (RealRounding mode in Modes)
            {
                double ed = r.ToDouble(mode), ad = f.ToDouble(mode);
                Assert.True(BitConverter.DoubleToUInt64Bits(ed) == BitConverter.DoubleToUInt64Bits(ad) || (double.IsNaN(ed) && double.IsNaN(ad)),
                    $"double {H(b)} {mode}: {ed:R} vs {ad:R}");
                float es = r.ToSingle(mode), @as = f.ToSingle(mode);
                Assert.True(BitConverter.SingleToUInt32Bits(es) == BitConverter.SingleToUInt32Bits(@as) || (float.IsNaN(es) && float.IsNaN(@as)),
                    $"float {H(b)} {mode}");
                Half eh = r.ToHalf(mode), ah = f.ToHalf(mode);
                Assert.True(BitConverter.HalfToUInt16Bits(eh) == BitConverter.HalfToUInt16Bits(ah) || (Half.IsNaN(eh) && Half.IsNaN(ah)),
                    $"half {H(b)} {mode}");
            }
        }
        Assert.Equal(1f, (float)(Float128)1f);
        Assert.Equal((Half)1.5, (Half)(Float128)(Half)1.5);
    }

    [Fact]
    public void IntegerConversions()
    {
        Random rng = new(10);
        for (int i = 0; i < 5000; i++)
        {
            Int128 v = new((ulong)rng.NextInt64() >> rng.Next(64), (ulong)rng.NextInt64());
            if (rng.Next(2) == 0)
            {
                v = -v;
            }
            UInt128 u = (UInt128)rng.NextInt64() << rng.Next(65) | (ulong)rng.NextInt64();
            foreach (RealRounding mode in Modes)
            {
                Assert.Equal(Real.FromBigInt(new BigInt(v)).ToBinary128Bits(mode), Float128.FromInt128(v, mode).ToBits());
                Assert.Equal(Real.FromBigInt(new BigInt(u)).ToBinary128Bits(mode), Float128.FromUInt128(u, mode).ToBits());
            }
            long l = rng.NextInt64() >> rng.Next(64);
            Float128 fl = l;
            Assert.Equal(l, (long)fl);
            Assert.Equal(Real.FromBigInt(new BigInt(l)).ToBinary128Bits(), fl.ToBits());
            ulong ul = (ulong)rng.NextInt64() * 2 + 1;
            Assert.Equal(ul, (ulong)(Float128)ul);
            // truncation of arbitrary finite values
            Float128 x = F(RandomBits(rng));
            if (Float128.IsFinite(x))
            {
                BigInt t = R(x).ToBigInt();
                Assert.Equal(t, (BigInt)x);
                if (t >= new BigInt(long.MinValue) && t <= new BigInt(long.MaxValue))
                {
                    Assert.Equal((long)t, (long)x);
                }
                else
                {
                    Assert.Throws<OverflowException>(() => (long)x);
                    Assert.Equal(x < 0 ? long.MinValue : long.MaxValue, long.CreateSaturating(x));
                }
            }
        }
        Assert.Equal(long.MinValue, (long)(Float128)long.MinValue);
        Assert.Equal(Int128.MinValue, (Int128)(Float128)(-Float128.ScaleB(1, 127)));
        Assert.Throws<OverflowException>(() => (Int128)Float128.ScaleB(1, 127));
        Assert.Equal(UInt128.MaxValue - (UInt128.MaxValue >> 113), (UInt128)Float128.BitDecrement(Float128.ScaleB(1, 128)));
        Assert.Throws<OverflowException>(() => (int)Float128.NaN);
        Assert.Throws<OverflowException>(() => (int)Float128.PositiveInfinity);
        Assert.Throws<OverflowException>(() => (byte)(Float128)(-1));
        Assert.Equal((byte)0, (byte)(Float128)(-0.75));
        Assert.Equal(-3, (int)(Float128)(-3.99));
        Assert.Equal(0, int.CreateSaturating(Float128.NaN));
        Assert.Equal(int.MaxValue, int.CreateTruncating(Float128.PositiveInfinity));
        Assert.Throws<OverflowException>(() => int.CreateChecked(Float128.NaN));
        Assert.Equal((Float128)42, Float128.CreateChecked(42));
        Assert.Equal((Float128)(-7), Float128.CreateTruncating((short)-7));
    }

    [Fact]
    public void OtherConversions()
    {
        Assert.Equal(Float128.Parse("0.1"), (Float128)0.1m);
        Assert.Equal(Float128.Parse("79228162514264337593543950335"), (Float128)decimal.MaxValue);
        Assert.Equal(1.5m, (decimal)(Float128)1.5);
        Assert.Equal(0.1m, (decimal)Float128.Parse("0.1"));
        Assert.Throws<OverflowException>(() => (decimal)Float128.ScaleB(1, 100));
        Assert.Equal(Sign, ((Float128)(-0.0m)).ToBits());
        Assert.Equal((Float128)1 / 3, (Float128)new BigRational(1, 3));
        Assert.Equal(new BigRational(3, 4), (BigRational)(Float128)0.75);
        BigInt big = BigInt.One << 20000;
        Assert.Equal(Float128.PositiveInfinity, (Float128)big);
        Assert.Equal(Float128.MaxValue, Float128.FromBigInt(big, RealRounding.TowardZero));
        BigInt n = (BigInt.One << 200) + 1;
        Assert.Equal(Real.FromBigInt(n).ToBinary128Bits(RealRounding.TowardPositive), Float128.FromBigInt(n, RealRounding.TowardPositive).ToBits());
        Real r = Real.Sqrt(2, 300);
        foreach (RealRounding mode in Modes)
        {
            Assert.Equal(r.ToBinary128Bits(mode), Float128.FromReal(r, mode).ToBits());
        }
        Float128 x = Float128.Pi;
        Real back = x;
        Assert.Equal(113, back.Precision);
        Assert.Equal(x, (Float128)back);
        Assert.Equal(1.0, double.CreateChecked(Float128.One));
        Assert.Equal((Half)0.5, Half.CreateChecked((Float128)0.5));
    }

    // ---------------------------------------------------------------
    // Comparisons and generic surface
    // ---------------------------------------------------------------

    [Fact]
    public void ComparisonsFollowIeee()
    {
        Float128 nan = Float128.NaN, one = 1, two = 2;
        Float128 nan2 = Float128.NaN;
        Assert.False(nan == nan2);
        Assert.True(nan != nan2);
        Assert.False(nan < one || nan > one || nan <= one || nan >= one);
        Assert.True(Float128.Zero == Float128.NegativeZero);
        Assert.True(Float128.NegativeZero <= Float128.Zero && !(Float128.NegativeZero < Float128.Zero));
        Assert.True(-two < -one && -one < Float128.NegativeZero && one < two && Float128.MaxValue < Float128.PositiveInfinity);
        Assert.True(nan.Equals(nan));
        Assert.Equal(Float128.Zero.GetHashCode(), Float128.NegativeZero.GetHashCode());
        Assert.Equal(-1, nan.CompareTo(one));
        Assert.Equal(0, nan.CompareTo(nan));
        Assert.Equal(1, one.CompareTo(nan));
        Assert.Equal(0, Float128.Zero.CompareTo(Float128.NegativeZero));
        Assert.True(Float128.CompareTotalOrder(Float128.NegativeZero, Float128.Zero) < 0);
        Assert.True(Float128.CompareTotalOrder(Float128.PositiveInfinity, nan) < 0);
        Assert.True(Float128.CompareTotalOrder(-nan, Float128.NegativeInfinity) < 0);

        Random rng = new(11);
        for (int i = 0; i < 5000; i++)
        {
            Float128 a = F(RandomBits(rng)), b = F(rng.Next(3) == 0 ? RandomBits(rng) : RelatedBits(rng, a.ToBits()));
            if (Float128.IsNaN(a) || Float128.IsNaN(b))
            {
                continue;
            }
            int c = Float128.IsZero(a) && Float128.IsZero(b) ? 0 : R(a).CompareTo(R(b)); // Real orders -0 below +0
            Assert.True(c < 0 == a < b && (c == 0) == (a == b), $"{H(a.ToBits())} {H(b.ToBits())} c={c}");
            Assert.Equal(c > 0, a > b);
            Assert.Equal(c <= 0, a <= b);
            Assert.Equal(Math.Sign(c), Math.Sign(a.CompareTo(b)));
        }
    }

    [Fact]
    public void MinMaxAndPredicates()
    {
        Float128 nan = Float128.NaN, one = 1, m1 = -1;
        Assert.True(Float128.IsNaN(Float128.Max(nan, one)));
        Assert.Equal(one, Float128.MaxNumber(nan, one));
        Assert.Equal(m1, Float128.MinNumber(m1, nan));
        Assert.Equal(UInt128.Zero, Float128.Max(Float128.NegativeZero, Float128.Zero).ToBits());
        Assert.Equal(Sign, Float128.Min(Float128.Zero, Float128.NegativeZero).ToBits());
        Assert.Equal(m1, Float128.MinMagnitude(m1, (Float128)2));
        Assert.Equal(one, Float128.MaxMagnitude(m1, one));
        Assert.Equal(one, Float128.CopySign(m1, Float128.Zero));
        Assert.True(Float128.IsInteger(4) && Float128.IsEvenInteger(4) && Float128.IsOddInteger(3) && !Float128.IsInteger((Float128)0.5));
        Assert.True(Float128.IsOddInteger(Float128.ScaleB(1, 112) + 1));
        Assert.True(Float128.IsEvenInteger(Float128.MaxValue));
        Assert.True(Float128.IsPow2(Float128.Epsilon) && Float128.IsPow2(8) && !Float128.IsPow2(6));
        Assert.Equal(-1, Float128.Sign(m1));
        Assert.Throws<ArithmeticException>(() => Float128.Sign(nan));
        Assert.Equal((Float128)1.5, Float128.Clamp((Float128)1.5, 1, 2));
        Assert.Equal(Float128.ScaleB(1, -16494), Float128.Epsilon);
        Assert.Equal((Float128)3, (Float128)7 % 4);
        Assert.Equal((Float128)(-3), (Float128)(-7) % 4);
        Float128 inc = 2;
        inc++;
        Assert.Equal((Float128)3, inc);
    }

    private static T GenericSum<T>(T a, T b) where T : IBinaryFloatingPointIeee754<T> => T.FusedMultiplyAdd(a, T.One, b);

    [Fact]
    public void GenericMath()
    {
        Assert.Equal((Float128)3, GenericSum<Float128>(1, 2));
        Float128 x = 1.5;
        IFloatingPoint<Float128> fp = x;
        Assert.Equal(113, fp.GetSignificandBitLength());
        Span<byte> buf = stackalloc byte[16];
        Assert.True(fp.TryWriteSignificandLittleEndian(buf, out int n));
        Assert.Equal(16, n);
        Assert.Equal(Hex(0x0001_8000_0000_0000, 0), System.Buffers.Binary.BinaryPrimitives.ReadUInt128LittleEndian(buf));
        Assert.True(fp.TryWriteExponentBigEndian(buf, out n));
        Assert.Equal(2, n);
        Assert.Equal(0, System.Buffers.Binary.BinaryPrimitives.ReadInt16BigEndian(buf));
        Assert.Equal(Hex(0x0000_8000_0000_0000, 0), BitAnd(x, Float128.FromBits(Hex(0x0000_8000_0000_0000, 0))).ToBits());
        Assert.Equal(2, ConvInt<Float128>((Float128)2.7));
    }

    private static int ConvInt<T>(T x) where T : IFloatingPoint<T> => T.ConvertToInteger<int>(x);

    private static T BitAnd<T>(T a, T b) where T : IBitwiseOperators<T, T, T> => a & b;

    // ---------------------------------------------------------------
    // Formatting and parsing
    // ---------------------------------------------------------------

    [Theory]
    [InlineData("0.1", "0.1")]
    [InlineData("1", "1")]
    [InlineData("-2.5", "-2.5")]
    [InlineData("100", "100")]
    [InlineData("1e33", "1E+33")]
    [InlineData("123456789012345678901234567890123", "123456789012345678901234567890123")]
    [InlineData("0.0001", "0.0001")]
    [InlineData("0.00001", "1E-05")]
    [InlineData("3.14159265358979323846264338327950288", "3.1415926535897932384626433832795028")]
    [InlineData("-0", "-0")]
    [InlineData("1.18973149535723176508575932662800702e4932", "1.189731495357231765085759326628007E+4932")]
    [InlineData("6.47517511943802511092443895822764655e-4966", "6E-4966")]
    public void ShortestToString(string input, string expected)
    {
        Float128 x = Float128.Parse(input, CultureInfo.InvariantCulture);
        Assert.Equal(expected, x.ToString(CultureInfo.InvariantCulture));
        Assert.Equal(x.ToBits(), Float128.Parse(expected, CultureInfo.InvariantCulture).ToBits());
    }

    [Fact]
    public void ToStringRoundTripsAndIsShortest()
    {
        Random rng = new(12);
        for (int i = 0; i < 3000; i++)
        {
            UInt128 b = RandomBits(rng);
            Float128 x = F(b);
            if (!Float128.IsFinite(x))
            {
                continue;
            }
            string s = x.ToString(CultureInfo.InvariantCulture);
            Assert.True(Float128.Parse(s, CultureInfo.InvariantCulture).ToBits() == b, $"{H(b)} -> {s}");
            if (Float128.IsZero(x) || i % 3 != 0)
            {
                continue;
            }
            // No decimal with one digit fewer lies in the rounding interval: neither neighbor round-trips.
            string mant = s.TrimStart('-').Split('E')[0].Replace(".", "").Trim('0');
            int n = mant.Length;
            if (n <= 1)
            {
                continue;
            }
            Real r = R(x);
            foreach (RealRounding mode in new[] { RealRounding.TowardZero, RealRounding.AwayFromZero })
            {
                string shorter = r.ToString(n - 1, mode);
                Assert.True(Float128.Parse(shorter, CultureInfo.InvariantCulture).ToBits() != b, $"{H(b)}: {s} but {shorter} also round-trips");
            }
        }
    }

    [Fact]
    public void ParseIsCorrectlyRounded()
    {
        Random rng = new(13);
        for (int i = 0; i < 3000; i++)
        {
            int digits = rng.Next(1, 60);
            char[] cs = new char[digits];
            for (int k = 0; k < digits; k++)
            {
                cs[k] = (char)('0' + rng.Next(10));
            }
            if (i % 5 == 0)
            {
                // exact binary128 midpoints
                Float128 v = F(RandomBits(rng) & ~Sign);
                if (!Float128.IsFinite(Float128.BitIncrement(v)) || Float128.IsZero(v))
                {
                    continue;
                }
                Real mid = Real.Add(R(v), R(Float128.BitIncrement(v)), 300).ShiftRight(1);
                string ms = ExactDecimal(mid.ToRational());
                CheckParse(ms);
                continue;
            }
            int exp = rng.Next(5) switch
            {
                0 => rng.Next(-4990, -4900),
                1 => rng.Next(4900, 4940),
                _ => rng.Next(-60, 60),
            };
            string s = (rng.Next(2) == 0 ? "-" : "") + new string(cs) + "e" + exp.ToString(CultureInfo.InvariantCulture);
            CheckParse(s);
        }

        static void CheckParse(string s)
        {
            foreach (RealRounding mode in Modes)
            {
                UInt128 expected = Oracle((p, m) => Real.Parse(s, p, m), mode);
                AssertBits(expected, Float128.Parse(s, mode), $"Parse({s}) {mode}");
            }
        }
    }

    /// <summary>The exact decimal expansion of a dyadic rational.</summary>
    private static string ExactDecimal(BigRational q)
    {
        BigInt num = q.Numerator, den = q.Denominator;
        int k = (int)den.TrailingZeroCount(); // den = 2^k
        BigInt scaled = num * BigInt.Pow(5, (ulong)k); // value = scaled / 10^k
        return scaled.ToString() + "e-" + k.ToString(CultureInfo.InvariantCulture);
    }

    [Fact]
    public void ParseSpecialsAndStyles()
    {
        Assert.True(Float128.IsNaN(Float128.Parse("NaN")));
        Assert.Equal(Float128.PositiveInfinity, Float128.Parse("Infinity", CultureInfo.InvariantCulture));
        Assert.Equal(Float128.NegativeInfinity, Float128.Parse("-Infinity", CultureInfo.InvariantCulture));
        Assert.Equal(Float128.NegativeInfinity, Float128.Parse("-inf", CultureInfo.InvariantCulture));
        Assert.Equal(Float128.PositiveInfinity, Float128.Parse("1e5000", CultureInfo.InvariantCulture));
        Assert.Equal(Float128.MaxValue, Float128.Parse("1e5000", RealRounding.TowardZero));
        Assert.Equal(UInt128.Zero, Float128.Parse("1e-5000", CultureInfo.InvariantCulture).ToBits());
        Assert.Equal(UInt128.One, Float128.Parse("1e-5000", RealRounding.TowardPositive).ToBits());
        Assert.Equal(Sign, Float128.Parse("-0.000", CultureInfo.InvariantCulture).ToBits());
        Assert.Equal((Float128)1234.5, Float128.Parse("1,234.5", CultureInfo.InvariantCulture));
        Assert.Equal((Float128)1234.5, Float128.Parse("1.234,5", CultureInfo.GetCultureInfo("de-DE")));
        Assert.Equal((Float128)0.5, Float128.Parse("  .5  ", CultureInfo.InvariantCulture));
        Assert.False(Float128.TryParse("1.2.3", CultureInfo.InvariantCulture, out _));
        Assert.False(Float128.TryParse("abc", CultureInfo.InvariantCulture, out _));
        Assert.False(Float128.TryParse("", CultureInfo.InvariantCulture, out _));
        Assert.Equal((Float128)2.5, Float128.Parse("2.5"u8, CultureInfo.InvariantCulture));
        Assert.True(Float128.TryParse("1e1000000000000000000000", CultureInfo.InvariantCulture, out Float128 big) && Float128.IsPositiveInfinity(big));
    }

    [Fact]
    public void FormatStrings()
    {
        CultureInfo inv = CultureInfo.InvariantCulture;
        Float128 x = Float128.Parse("1234.5678", inv);
        Assert.Equal("1234.57", x.ToString("F2", inv));
        Assert.Equal("1,234.568", x.ToString("N3", inv));
        Assert.Equal("1.234568E+003", x.ToString("E", inv));
        Assert.Equal("1.23e+003", x.ToString("e2", inv));
        Assert.Equal("1234.6", x.ToString("G5", inv));
        Assert.Equal("1.2E+03", x.ToString("G2", inv));
        Assert.Equal("1234.5678", x.ToString("R", inv));
        Assert.Equal("0.12", ((Float128)0.125).ToString("F2", inv)); // exact ties to even, like double on .NET 10
        Assert.Equal(0.125.ToString("F2", inv), ((Float128)0.125).ToString("F2", inv));
        Assert.Equal("0.1000000000000000055511151231257827", ((Float128)0.1).ToString("G34", inv));
        Assert.Equal("-0.00", ((Float128)(-0.001)).ToString("F2", inv));
        Assert.Equal("1234,57", x.ToString("F2", CultureInfo.GetCultureInfo("de-DE")));
        Assert.Equal("NaN", Float128.NaN.ToString(inv));
        Assert.Equal("-Infinity", Float128.NegativeInfinity.ToString(inv));
        Assert.Throws<FormatException>(() => x.ToString("X", inv));
        Span<char> buf = stackalloc char[64];
        Assert.True(x.TryFormat(buf, out int n, "F1", inv));
        Assert.Equal("1234.6", buf[..n].ToString());
        Span<byte> u8 = stackalloc byte[64];
        Assert.True(x.TryFormat(u8, out n, default, inv));
        Assert.Equal("1234.5678", System.Text.Encoding.UTF8.GetString(u8[..n]));
        Random rng = new(14);
        for (int i = 0; i < 300; i++)
        {
            double d = rng.NextDouble() * Math.Pow(10, rng.Next(-8, 20)) * (rng.Next(2) == 0 ? -1 : 1);
            Float128 f = d;
            foreach (string fmt in new[] { "E10", "F4", "G7", "E0", "N2" })
            {
                Assert.Equal(d.ToString(fmt, inv), f.ToString(fmt, inv));
            }
        }
    }

    // ---------------------------------------------------------------
    // Elementary functions (correct rounding through Real)
    // ---------------------------------------------------------------

    private static Float128 RandomArg(Random rng, double lo, double hi)
    {
        if (rng.Next(6) == 0)
        {
            // small arguments around the 2^-200 shortcut threshold
            return F(((UInt128)(uint)(0x3FFF - rng.Next(150, 260)) << 112) | RandomFrac(rng) | (rng.Next(2) == 0 ? Sign : UInt128.Zero));
        }
        Float128 u = Float128.FromBits(((UInt128)0x3FFF << 112) | RandomFrac(rng)) - 1; // [0, 1)
        return (Float128)lo + u * (Float128)(hi - lo);
    }

    private static void CheckUnary(int seed, int count, double lo, double hi, Func<Float128, RealRounding, Float128> f,
        Func<Real, int, RealRounding, Real> g, string name)
    {
        Random rng = new(seed);
        for (int i = 0; i < count; i++)
        {
            Float128 x = RandomArg(rng, lo, hi);
            Real rx = R(x);
            RealRounding mode = i % 2 == 0 ? RealRounding.NearestEven : Modes[rng.Next(Modes.Length)];
            AssertBits(Oracle((p, m) => g(rx, p, m), mode), f(x, mode), $"{name}({H(x.ToBits())}) {mode}");
        }
    }

    [Fact]
    public void ExpLogFamily()
    {
        CheckUnary(20, 150, -60, 60, Float128.Exp, Real.Exp, "Exp");
        CheckUnary(21, 60, -11450, 11360, Float128.Exp, Real.Exp, "ExpWide");
        CheckUnary(22, 120, -100, 100, Float128.Exp2, Real.Exp2, "Exp2");
        CheckUnary(23, 120, -40, 40, Float128.Exp10, Real.Exp10, "Exp10");
        CheckUnary(24, 120, -50, 50, Float128.ExpM1, Real.ExpM1, "ExpM1");
        CheckUnary(25, 150, 0, 1e6, Float128.Log, Real.Log, "Log");
        CheckUnary(26, 120, 0, 1e6, Float128.Log2, Real.Log2, "Log2");
        CheckUnary(27, 120, 0, 1e6, Float128.Log10, Real.Log10, "Log10");
        CheckUnary(28, 120, -0.9, 100, Float128.LogP1, Real.LogP1, "LogP1");
        CheckUnary(29, 100, 0, 1e30, Float128.Cbrt, Real.Cbrt, "Cbrt");
        CheckUnary(30, 100, -20, 20, Float128.Sinh, Real.Sinh, "Sinh");
        CheckUnary(31, 100, -20, 20, Float128.Cosh, Real.Cosh, "Cosh");
        CheckUnary(32, 100, -20, 20, Float128.Tanh, Real.Tanh, "Tanh");
        CheckUnary(33, 100, -1e6, 1e6, Float128.Asinh, Real.Asinh, "Asinh");
        CheckUnary(34, 100, 1, 1e6, Float128.Acosh, Real.Acosh, "Acosh");
        CheckUnary(35, 100, -1, 1, Float128.Atanh, Real.Atanh, "Atanh");
        // Subnormal / overflow edges
        Assert.Equal(PosInf, Float128.Exp(11400).ToBits());
        Assert.Equal(MaxFinite, Float128.Exp(11400, RealRounding.TowardZero).ToBits());
        Assert.Equal(UInt128.Zero, Float128.Exp(-12000).ToBits());
        Assert.Equal(UInt128.One, Float128.Exp(-12000, RealRounding.TowardPositive).ToBits());
        Assert.Equal(UInt128.One, Float128.Exp2(-16494).ToBits());
        Assert.Equal(Float128.NegativeOne, Float128.ExpM1(-1000));
        Assert.Equal(Float128.BitIncrement(Float128.NegativeOne), Float128.ExpM1(-1000, RealRounding.TowardZero));
        Assert.Equal(Float128.One, Float128.Exp(0));
        Assert.Equal(Float128.Zero, Float128.Log(1));
        Assert.Equal(Float128.NegativeInfinity, Float128.Log(0));
        Assert.True(Float128.IsNaN(Float128.Log(-1)));
        Assert.Equal((Float128)10, Float128.Log2(1024));
        Assert.Equal((Float128)2, Float128.Cbrt(8));
    }

    [Fact]
    public void TinyArgumentShortcuts()
    {
        // Deep-tiny and subnormal arguments take the Taylor-sign shortcut; check it against the oracle.
        Float128[] xs = [Float128.ScaleB((Float128)1.375, -201), -Float128.ScaleB(Float128.One / 3, -260), Float128.ScaleB(Float128.Pi, -700)];
        var fs = new (string, Func<Float128, RealRounding, Float128>, Func<Real, int, RealRounding, Real>)[]
        {
            ("Exp", Float128.Exp, Real.Exp), ("ExpM1", Float128.ExpM1, Real.ExpM1), ("LogP1", Float128.LogP1, Real.LogP1),
            ("Sin", Float128.Sin, Real.Sin), ("Cos", Float128.Cos, Real.Cos), ("Tan", Float128.Tan, Real.Tan),
            ("Atan", Float128.Atan, Real.Atan), ("Sinh", Float128.Sinh, Real.Sinh), ("Tanh", Float128.Tanh, Real.Tanh),
        };
        foreach (var (name, f, g) in fs)
        {
            foreach (Float128 x in xs)
            {
                Real rx = R(x);
                foreach (RealRounding mode in Modes)
                {
                    AssertBits(Oracle((p, m) => g(rx, p, m), mode), f(x, mode), $"{name}({H(x.ToBits())}) {mode}");
                }
            }
        }
        Float128 e3 = Float128.Epsilon * 3;
        Assert.Equal(e3, Float128.Sin(e3));
        Assert.Equal(Float128.Epsilon * 2, Float128.Sin(e3, RealRounding.TowardZero));
        Assert.Equal(Float128.Epsilon * 4, Float128.Atanh(e3, RealRounding.TowardPositive));
        Assert.Equal(-Float128.Epsilon * 2, Float128.ExpM1(-e3, RealRounding.TowardZero));
        Assert.Equal(Float128.BitDecrement(Float128.One), Float128.Cos(e3, RealRounding.TowardZero));
        Assert.Equal(Float128.BitIncrement(Float128.One), Float128.Exp(e3, RealRounding.TowardPositive));
    }

    [Fact]
    public void TrigFamily()
    {
        CheckUnary(40, 150, -10, 10, Float128.Sin, Real.Sin, "Sin");
        CheckUnary(41, 150, -10, 10, Float128.Cos, Real.Cos, "Cos");
        CheckUnary(42, 100, -10, 10, Float128.Tan, Real.Tan, "Tan");
        CheckUnary(43, 40, -1e30, 1e30, Float128.Sin, Real.Sin, "SinLarge");
        CheckUnary(44, 100, -1, 1, Float128.Asin, Real.Asin, "Asin");
        CheckUnary(45, 100, -1, 1, Float128.Acos, Real.Acos, "Acos");
        CheckUnary(46, 100, -1e6, 1e6, Float128.Atan, Real.Atan, "Atan");
        CheckUnary(47, 80, -10, 10, Float128.SinPi, Real.SinPi, "SinPi");
        CheckUnary(48, 80, -10, 10, Float128.CosPi, Real.CosPi, "CosPi");
        CheckUnary(49, 60, -10, 10, Float128.TanPi, Real.TanPi, "TanPi");
        CheckUnary(50, 60, -1, 1, Float128.AsinPi, Real.AsinPi, "AsinPi");
        CheckUnary(51, 60, -1, 1, Float128.AcosPi, Real.AcosPi, "AcosPi");
        CheckUnary(52, 60, -100, 100, Float128.AtanPi, Real.AtanPi, "AtanPi");
        Assert.Equal(Float128.Zero, Float128.Sin(0));
        Assert.Equal(Float128.One, Float128.SinPi((Float128)0.5));
        Assert.Equal(Float128.Sin(Float128.Pi), Float128.SinCos(Float128.Pi).Sin);
    }

    [Fact]
    public void BinaryFunctions()
    {
        Random rng = new(60);
        for (int i = 0; i < 150; i++)
        {
            Float128 x = RandomArg(rng, 0, 100), y = RandomArg(rng, -40, 40);
            Real rx = R(x), ry = R(y);
            RealRounding mode = Modes[i % Modes.Length];
            AssertBits(Oracle((p, m) => Real.Pow(rx, ry, p, m), mode), Float128.Pow(x, y, mode), $"Pow({x}, {y}) {mode}");
            AssertBits(Oracle((p, m) => Real.Atan2(ry, rx, p, m), mode), Float128.Atan2(y, x, mode), $"Atan2 {mode}");
            AssertBits(Oracle((p, m) => Real.Atan2Pi(ry, rx, p, m), mode), Float128.Atan2Pi(y, x, mode), $"Atan2Pi {mode}");
            AssertBits(Oracle((p, m) => Real.Hypot(rx, ry, p, m), mode), Float128.Hypot(x, y, mode), $"Hypot {mode}");
            int n = rng.Next(-7, 8);
            if (n != 0)
            {
                AssertBits(Oracle((p, m) => Real.RootN(rx, (long)n, p, m), mode), Float128.RootN(x, n, mode), $"RootN({x}, {n}) {mode}");
            }
        }
        Assert.Equal((Float128)8, Float128.Pow(2, 3));
        Assert.Equal((Float128)(-8), Float128.Pow(-2, 3));
        Assert.Equal((Float128)5, Float128.Hypot(3, 4));
        Assert.Equal(Float128.PositiveInfinity, Float128.Pow(10, 5000));
        Assert.Equal(Float128.NegativeInfinity, Float128.Pow(-10, 5001));
        Assert.Equal(Float128.Zero, Float128.Pow(10, -5000));
        Assert.Equal(Float128.PositiveInfinity, Float128.Pow(Float128.BitIncrement(1), Float128.ScaleB(1, 130)));
        Assert.Equal(Float128.Pow((Float128)1.5, 7), Float128.Pow((Float128)1.5, 7));
        AssertBits(Oracle((p, m) => Real.Pow(R(1.5), R(-28100), p, m), RealRounding.NearestEven), Float128.Pow((Float128)1.5, -28100), "Pow subnormal");
    }
}
