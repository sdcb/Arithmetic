using System.Numerics;
using Xunit;

namespace Sdcb.Arithmetic.Tests;

public class NumberTheoryTests
{
    [Fact]
    public void Gcd_MatchesBigInteger()
    {
        var rng = new Random(70);
        for (int i = 0; i < 300; i++)
        {
            int maxLimbs = i < 250 ? 8 : 100;
            var (a, ba) = TestHelpers.RandomPair(rng, maxLimbs);
            var (b, bb) = TestHelpers.RandomPair(rng, maxLimbs);
            TestHelpers.AssertEqual(BigInteger.GreatestCommonDivisor(ba, bb), BigInt.Gcd(a, b), $"gcd #{i}");
        }
    }

    [Fact]
    public void Gcd_WithCommonFactor()
    {
        var rng = new Random(71);
        for (int i = 0; i < 50; i++)
        {
            var (f, bf) = TestHelpers.RandomPair(rng, 20, allowNegative: false);
            var (a, ba) = TestHelpers.RandomPair(rng, 20, allowNegative: false);
            var (b, bb) = TestHelpers.RandomPair(rng, 20, allowNegative: false);
            TestHelpers.AssertEqual(
                BigInteger.GreatestCommonDivisor(ba * bf, bb * bf),
                BigInt.Gcd(a * f, b * f), $"gcd common #{i}");
        }
    }

    [Fact]
    public void GcdExt_BezoutIdentity()
    {
        var rng = new Random(72);
        for (int i = 0; i < 300; i++)
        {
            var (a, _) = TestHelpers.RandomPair(rng, 6);
            var (b, _) = TestHelpers.RandomPair(rng, 6);
            BigInt g = BigInt.ExtendedGcd(a, b, out BigInt x, out BigInt y);
            Assert.Equal(g, a * x + b * y);
            Assert.Equal(BigInt.Gcd(a, b), g);
        }
    }

    [Fact]
    public void GcdExt_LargeOperands()
    {
        var rng = new Random(172);
        for (int i = 0; i < 120; i++)
        {
            var (a, ba) = TestHelpers.RandomPair(rng, 60);
            var (b, bb) = TestHelpers.RandomPair(rng, 60);
            if (i % 5 == 0)
            {
                var (f, _) = TestHelpers.RandomPair(rng, 10, allowNegative: false);
                a *= f;
                b *= f;
                ba = a.ToBigInteger();
                bb = b.ToBigInteger();
            }
            BigInt g = BigInt.ExtendedGcd(a, b, out BigInt x, out BigInt y);
            Assert.Equal(g, a * x + b * y);
            TestHelpers.AssertEqual(BigInteger.GreatestCommonDivisor(ba, bb), g, $"gcdext #{i}");
            // Euclid-minimal cofactors: |x| <= |b|/(2g) + 1.
            if (!g.IsZero && !b.IsZero)
            {
                Assert.True(BigInt.CompareAbs(x, b.Abs() / g) <= 0, $"cofactor size #{i}");
            }
        }
    }

    [Fact]
    public void Gcd_SubquadraticSizes_MatchBigInteger()
    {
        var rng = new Random(174);
        foreach (int limbs in (int[])[3, 40, 120, 250, 400, 700, 1500, 3000])
        {
            for (int rep = 0; rep < 3; rep++)
            {
                byte[] x = new byte[limbs * 8], y = new byte[(limbs - rng.Next(Math.Min(limbs, 30))) * 8], f = new byte[rng.Next(1, 40 * 8)];
                rng.NextBytes(x);
                rng.NextBytes(y);
                rng.NextBytes(f);
                BigInt cf = BigInt.FromMagnitudeBytes(f) | 1;
                BigInt a = BigInt.FromMagnitudeBytes(x) * cf, b = BigInt.FromMagnitudeBytes(y) * cf;
                if (rep == 1)
                {
                    b = -b;
                }
                BigInteger ba = a.ToBigInteger(), bb = b.ToBigInteger();
                BigInteger expected = BigInteger.GreatestCommonDivisor(ba, bb);
                TestHelpers.AssertEqual(expected, BigInt.Gcd(a, b), $"gcd {limbs}");
                BigInt g = BigInt.ExtendedGcd(a, b, out BigInt s, out BigInt t);
                TestHelpers.AssertEqual(expected, g, $"gcdext {limbs}");
                Assert.Equal(g, a * s + b * t);
                // Canonical (mpz_gcdext) cofactor bounds.
                Assert.True(BigInt.CompareAbs(s << 1, b.Abs() / g) <= 0, $"|s| bound {limbs}");
                Assert.True(BigInt.CompareAbs(t << 1, a.Abs() / g) <= 0, $"|t| bound {limbs}");
            }
        }
    }

    [Fact]
    public void GcdExt_CanonicalCofactorsSmall()
    {
        // Values from mpz_gcdext's documented normalization.
        (long a, long b, long g, long s, long t)[] cases =
        [
            (4, 6, 2, -1, 1), (6, 4, 2, 1, -1), (5, 5, 5, 0, 1), (-5, 5, 5, 0, 1), (5, -5, 5, 0, -1),
            (7, 0, 7, 1, 0), (-7, 0, 7, -1, 0), (0, -7, 7, 0, -1), (3, 6, 3, 1, 0), (6, 3, 3, 0, 1),
            (240, 46, 2, -9, 47), (2, 4, 2, 1, 0), (4, 2, 2, 0, 1), (1, 2, 1, 1, 0), (2, 1, 1, 0, 1),
        ];
        foreach (var (a, b, g, s, t) in cases)
        {
            BigInt gg = BigInt.ExtendedGcd(a, b, out BigInt ss, out BigInt tt);
            Assert.Equal((BigInt)g, gg);
            Assert.True((BigInt)s == ss && (BigInt)t == tt, $"gcdext({a},{b}) = ({gg},{ss},{tt}), expected ({g},{s},{t})");
        }
    }

    [Fact]
    public void Jacobi_LargeMatchesEulerAndMultiplicativity()
    {
        var rng = new Random(173);
        for (int i = 0; i < 60; i++)
        {
            var (p0, _) = TestHelpers.RandomPair(rng, 1 + i / 4, allowNegative: false);
            BigInt p = (p0 + 3).NextPrime();
            var (q0, _) = TestHelpers.RandomPair(rng, 1 + i / 6, allowNegative: false);
            BigInt q = (q0 + 3).NextPrime();
            for (int k = 0; k < 5; k++)
            {
                var (a, _) = TestHelpers.RandomPair(rng, 2 + i / 3);
                BigInt e = BigInt.ModPow(BigInt.Mod(a, p), (p - 1) >> 1, p);
                int euler = e.IsZero ? 0 : e.IsOne ? 1 : -1;
                Assert.Equal(euler, BigInt.Jacobi(a, p));
                Assert.Equal(BigInt.Jacobi(a, p) * BigInt.Jacobi(a, q), BigInt.Jacobi(a, p * q));
            }
        }
    }

    [Fact]
    public void ModInverse_Works()
    {
        var rng = new Random(73);
        for (int i = 0; i < 200; i++)
        {
            var (a, _) = TestHelpers.RandomPair(rng, 5);
            var (m, _) = TestHelpers.RandomPair(rng, 5, allowNegative: false);
            if (m.IsZero || m.IsOne)
            {
                continue;
            }
            if (BigInt.TryModInverse(a, m, out BigInt inv))
            {
                Assert.Equal(BigInt.One, BigInt.Mod(a * inv, m));
            }
            else
            {
                Assert.NotEqual(BigInt.One, BigInt.Gcd(a, m));
            }
        }
    }

    [Fact]
    public void PowMod_MatchesBigInteger()
    {
        var rng = new Random(74);
        for (int i = 0; i < 300; i++)
        {
            var (b, bb) = TestHelpers.RandomPair(rng, 6);
            var (e, be) = TestHelpers.RandomPair(rng, 3, allowNegative: false);
            var (m, bm) = TestHelpers.RandomPair(rng, 6, allowNegative: false);
            if (m.IsZero || m.IsOne)
            {
                continue;
            }
            BigInteger expected = BigInteger.ModPow(bb, be, bm);
            if (expected < 0)
            {
                expected += bm;
            }
            TestHelpers.AssertEqual(expected, BigInt.ModPow(b, e, m), $"powmod #{i} (odd={m.IsOdd})");
        }
    }

    [Fact]
    public void PowMod_OddModulusSizes_MatchBigInteger()
    {
        // Montgomery with row REDC (small n) and REDC by multiplication (larger n), including
        // moduli just below a power of two and with zero high/low limbs in the operands.
        var rng = new Random(76);
        foreach (int n in (int[])[1, 2, 3, 8, 16, 23, 24, 25, 31, 40, 64, 79, 80, 81, 100, 130, 170])
        {
            for (int kind = 0; kind < 3; kind++)
            {
                byte[] mb = new byte[n * 8];
                rng.NextBytes(mb);
                if (kind == 1)
                {
                    mb.AsSpan().Fill(0xFF);
                }
                mb[^1] |= 0x80;
                mb[0] |= 1;
                byte[] bbytes = new byte[n * 8];
                rng.NextBytes(bbytes);
                if (kind == 2)
                {
                    bbytes.AsSpan(0, bbytes.Length / 2).Clear();
                }
                bbytes[^1] &= 0x3F;
                BigInteger bm = new(mb, isUnsigned: true), bb = new(bbytes, isUnsigned: true);
                BigInteger be = new BigInteger(rng.NextInt64(1, long.MaxValue)) * rng.Next(1, 1000);
                BigInt m = BigInt.FromMagnitudeBytes(mb), b = BigInt.FromMagnitudeBytes(bbytes);
                TestHelpers.AssertEqual(BigInteger.ModPow(bb, be, bm), BigInt.ModPow(b, (BigInt)be, m), $"powmod n={n} kind={kind}");
            }
        }
    }

    [Fact]
    public void PowMod_LargeOddModulus()
    {
        var rng = new Random(75);
        var (b, bb) = TestHelpers.RandomPair(rng, 32, allowNegative: false);
        var (e, be) = TestHelpers.RandomPair(rng, 32, allowNegative: false);
        var (m, bm) = TestHelpers.RandomPair(rng, 32, allowNegative: false);
        m = m.IsEven ? m + BigInt.One : m;
        bm = m.ToBigInteger();
        TestHelpers.AssertEqual(BigInteger.ModPow(bb, be, bm), BigInt.ModPow(b, e, m), "large powmod");
    }

    [Fact]
    public void PowMod_NegativeExponent()
    {
        BigInt m = 101;
        BigInt r = BigInt.ModPow(3, new BigInt(-5), m);
        Assert.Equal(BigInt.One, BigInt.Mod(r * BigInt.ModPow(3, 5, m), m));
    }

    [Fact]
    public void Jacobi_MatchesNaive()
    {
        // Naive Jacobi via definition for small values.
        static int NaiveJacobi(long a, long n)
        {
            a %= n;
            if (a < 0)
            {
                a += n;
            }
            int result = 1;
            while (a != 0)
            {
                while ((a & 1) == 0)
                {
                    a >>= 1;
                    if ((n & 7) == 3 || (n & 7) == 5)
                    {
                        result = -result;
                    }
                }
                (a, n) = (n, a);
                if ((a & 3) == 3 && (n & 3) == 3)
                {
                    result = -result;
                }
                a %= n;
            }
            return n == 1 ? result : 0;
        }

        var rng = new Random(76);
        for (int i = 0; i < 500; i++)
        {
            long a = rng.NextInt64(-100000, 100000);
            long n = rng.NextInt64(1, 100000) | 1;
            Assert.Equal(NaiveJacobi(a, n), BigInt.Jacobi(a, n));
        }
    }

    [Fact]
    public void Kronecker_BasicProperties()
    {
        Assert.Equal(1, BigInt.Kronecker(1, 0));
        Assert.Equal(0, BigInt.Kronecker(2, 0));
        Assert.Equal(0, BigInt.Kronecker(4, 2));
        Assert.Equal(1, BigInt.Kronecker(7, 2));  // 7 ≡ -1 mod 8
        Assert.Equal(-1, BigInt.Kronecker(3, 2)); // 3 ≡ 3 mod 8
        // Multiplicativity spot check: (a/mn) = (a/m)(a/n) for odd m, n > 0.
        var rng = new Random(77);
        for (int i = 0; i < 200; i++)
        {
            long a = rng.NextInt64(-1000, 1000);
            long m = rng.NextInt64(1, 1000) | 1;
            long n = rng.NextInt64(1, 1000) | 1;
            Assert.Equal(BigInt.Kronecker(a, m) * BigInt.Kronecker(a, n), BigInt.Kronecker(a, m * n));
        }
    }

    [Fact]
    public void IsProbablePrime_KnownValues()
    {
        Assert.Equal(Primality.Composite, BigInt.Zero.GetPrimality());
        Assert.Equal(Primality.Composite, BigInt.One.GetPrimality());
        Assert.Equal(Primality.Prime, BigInt.Two.GetPrimality());
        Assert.Equal(Primality.Prime, new BigInt(3).GetPrimality());
        Assert.Equal(Primality.Composite, new BigInt(4).GetPrimality());
        Assert.Equal(Primality.Prime, new BigInt(2147483647).GetPrimality()); // 2^31-1
        Assert.Equal(Primality.Composite, new BigInt(2147483649).GetPrimality());
        // Carmichael numbers must be rejected.
        Assert.Equal(Primality.Composite, new BigInt(561).GetPrimality());
        Assert.Equal(Primality.Composite, new BigInt(41041).GetPrimality());
        // M127 = 2^127 - 1 is prime.
        BigInt m127 = (BigInt.One << 127) - BigInt.One;
        Assert.True(m127.IsProbablePrime());
        // 2^127 + 1 = 3 * ... composite
        Assert.Equal(Primality.Composite, ((BigInt.One << 127) + BigInt.One).GetPrimality());
        // A 256-bit prime (secp256k1 field prime).
        BigInt p256 = BigInt.Parse("fffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f", 16);
        Assert.True(p256.IsProbablePrime());
        Assert.Equal(Primality.Composite, (p256 * p256).GetPrimality());
    }

    [Fact]
    public void NextPrime_Works()
    {
        Assert.Equal(BigInt.Two, BigInt.Zero.NextPrime());
        Assert.Equal(new BigInt(3), BigInt.Two.NextPrime());
        Assert.Equal(new BigInt(11), new BigInt(8).NextPrime());
        Assert.Equal(new BigInt(101), new BigInt(97).NextPrime());
        BigInt p = (BigInt.One << 128).NextPrime();
        Assert.True(p.IsProbablePrime());
        Assert.True(p > BigInt.One << 128);
    }

    [Fact]
    public void SqrtRem_Identity()
    {
        var rng = new Random(78);
        for (int i = 0; i < 300; i++)
        {
            int maxLimbs = i < 250 ? 8 : 80;
            var (a, _) = TestHelpers.RandomPair(rng, maxLimbs, allowNegative: false);
            BigInt s = BigInt.Sqrt(a, out BigInt rem);
            Assert.True(rem.Sign >= 0);
            Assert.Equal(a, s.Square() + rem);
            // (s+1)^2 > a
            Assert.True((s + BigInt.One).Square() > a);
        }
        Assert.Equal(new BigInt(3), BigInt.Sqrt(new BigInt(15)));
        Assert.Equal(new BigInt(4), BigInt.Sqrt(new BigInt(16)));
    }

    [Fact]
    public void RootRem_Identity()
    {
        var rng = new Random(79);
        for (int i = 0; i < 200; i++)
        {
            var (a, _) = TestHelpers.RandomPair(rng, 10, allowNegative: false);
            uint k = (uint)rng.Next(2, 10);
            BigInt root = BigInt.Root(a, k, out BigInt rem);
            Assert.True(rem.Sign >= 0);
            Assert.Equal(a, BigInt.Pow(root, k) + rem);
            Assert.True(BigInt.Pow(root + BigInt.One, k) > a);
        }
        // Negative odd roots.
        Assert.Equal(new BigInt(-3), BigInt.Root(new BigInt(-27), 3));
        Assert.Equal(new BigInt(-2), BigInt.Root(new BigInt(-10), 3, out BigInt r2));
        Assert.Equal(new BigInt(-2), r2); // -10 = (-2)^3 - 2
        Assert.Throws<ArithmeticException>(() => BigInt.Root(new BigInt(-4), 2));
    }

    [Fact]
    public void PerfectSquareAndPower()
    {
        var rng = new Random(80);
        for (int i = 0; i < 100; i++)
        {
            var (a, _) = TestHelpers.RandomPair(rng, 5, allowNegative: false);
            Assert.True(BigInt.IsPerfectSquare(a.Square()));
            Assert.True(BigInt.IsPerfectPower(a.IsZero ? BigInt.Zero : BigInt.Pow(a, 3)));
        }
        Assert.False(BigInt.IsPerfectSquare(new BigInt(17)));
        Assert.True(BigInt.IsPerfectPower(new BigInt(-27)));
        Assert.False(BigInt.IsPerfectPower(new BigInt(-4)));
    }

    [Fact]
    public void Factorial_KnownValues()
    {
        Assert.Equal(BigInt.One, BigInt.Factorial(0));
        Assert.Equal(new BigInt(120), BigInt.Factorial(5));
        Assert.Equal(new BigInt(3628800), BigInt.Factorial(10));
        Assert.Equal(BigInt.Parse("265252859812191058636308480000000"), BigInt.Factorial(30));
        // 1000! has 2568 digits; verify against BigInteger.
        BigInteger expected = BigInteger.One;
        for (int i = 2; i <= 1000; i++)
        {
            expected *= i;
        }
        TestHelpers.AssertEqual(expected, BigInt.Factorial(1000), "1000!");
    }

    [Fact]
    public void DoubleFactorial_KnownValues()
    {
        Assert.Equal(new BigInt(105), BigInt.DoubleFactorial(7));  // 1*3*5*7
        Assert.Equal(new BigInt(384), BigInt.DoubleFactorial(8));  // 2*4*6*8
        Assert.Equal(new BigInt(280), BigInt.MultiFactorial(10, 3)); // 10*7*4*1... = 280
    }

    [Fact]
    public void Primorial_KnownValues()
    {
        Assert.Equal(new BigInt(2), BigInt.Primorial(2));
        Assert.Equal(new BigInt(30), BigInt.Primorial(6));   // 2*3*5
        Assert.Equal(new BigInt(210), BigInt.Primorial(7));  // 2*3*5*7
        Assert.Equal(new BigInt(210), BigInt.Primorial(10));
    }

    [Fact]
    public void Binomial_KnownValues()
    {
        Assert.Equal(new BigInt(252), BigInt.Binomial(10, 5));
        Assert.Equal(BigInt.One, BigInt.Binomial(10, 0));
        Assert.Equal(BigInt.Zero, BigInt.Binomial(5, 10));
        Assert.Equal(BigInt.Parse("40225345056"), BigInt.Binomial(40, 15));
        // Pascal identity.
        var rng = new Random(81);
        for (int i = 0; i < 50; i++)
        {
            ulong n = (ulong)rng.Next(2, 200);
            ulong k = (ulong)rng.Next(1, (int)n);
            Assert.Equal(BigInt.Binomial(n, k), BigInt.Binomial(n - 1, k - 1) + BigInt.Binomial(n - 1, k));
        }
        // Negative n: C(-3, 2) = 6
        Assert.Equal(new BigInt(6), BigInt.Binomial(new BigInt(-3), 2));
        Assert.Equal(new BigInt(-10), BigInt.Binomial(new BigInt(-3), 3));
    }

    [Fact]
    public void SqrtRem_AllSizesAndEdgeCases()
    {
        var rng = new Random(401);
        void Check(BigInt v, string what)
        {
            BigInt s = BigInt.Sqrt(v, out BigInt r);
            Assert.True(r.Sign >= 0 && s * s + r == v && (s + 1) * (s + 1) > v, $"sqrtrem {what}");
            Assert.Equal(s, BigInt.Sqrt(v));
        }
        for (int limbs = 1; limbs <= 300; limbs += limbs < 40 ? 1 : 17)
        {
            for (int rep = 0; rep < 4; rep++)
            {
                byte[] x = new byte[limbs * 8 - rng.Next(8)];
                rng.NextBytes(x);
                BigInt v = BigInt.FromMagnitudeBytes(x);
                if (v.IsZero)
                {
                    continue;
                }
                Check(v, $"random {limbs}");
                BigInt sq = BigInt.Sqrt(v);
                Check(sq * sq, $"square {limbs}");
                Check(sq * sq - 1, $"square-1 {limbs}");
                Check((sq + 1) * (sq + 1) - 1, $"next square-1 {limbs}");
            }
            Check((BigInt.One << (64 * limbs)) - 1, $"all ones {limbs}");
            Check(BigInt.One << (64 * limbs - 1), $"power of two {limbs}");
            Check(BigInt.One << (64 * limbs - 2), $"power of four {limbs}");
        }
        byte[] big = new byte[2000 * 8];
        rng.NextBytes(big);
        Check(BigInt.FromMagnitudeBytes(big), "2000 limbs");
        foreach (ulong u in (ulong[])[1, 2, 3, 4, 15, 16, 17, ulong.MaxValue, ulong.MaxValue - 1, 1UL << 62, (1UL << 32) * (1UL << 32 - 1)])
        {
            Check(new BigInt(u), $"word {u}");
        }
    }

    [Fact]
    public void FactorialBinomial_MatchNaiveProducts()
    {
        BigInteger f = BigInteger.One;
        for (int n = 1; n <= 700; n++)
        {
            f *= n;
            if (n < 60 || n % 37 == 0)
            {
                TestHelpers.AssertEqual(f, BigInt.Factorial((ulong)n), $"{n}!");
            }
        }
        var rng = new Random(400);
        for (int i = 0; i < 120; i++)
        {
            ulong n = (ulong)rng.Next(0, 3000);
            ulong k = (ulong)rng.Next(0, (int)n + 2);
            BigInteger expected = BigInteger.Zero;
            if (k <= n)
            {
                ulong kk = Math.Min(k, n - k);
                BigInteger num = BigInteger.One, den = BigInteger.One;
                for (ulong j = 0; j < kk; j++)
                {
                    num *= n - j;
                    den *= j + 1;
                }
                expected = num / den;
            }
            TestHelpers.AssertEqual(expected, BigInt.Binomial(n, k), $"C({n},{k})");
        }
        // Large n, small k (product path) and huge n via BigInt overload.
        TestHelpers.AssertEqual(BigInteger.Parse("4950"), BigInt.Binomial(100UL, 2UL), "C(100,2)");
        BigInt big = BigInt.Pow(10, 30);
        BigInt c3 = BigInt.Binomial(big, 3);
        Assert.Equal(big * (big - 1) * (big - 2) / 6, c3);
    }

    [Fact]
    public void Fibonacci_LargeIdentities()
    {
        // Cassini: F(n-1) F(n+1) - F(n)^2 = (-1)^n, across the table/doubling boundary.
        foreach (ulong n in (ulong[])[90, 92, 93, 94, 95, 187, 188, 1000, 4097, 100001])
        {
            (BigInt fn1, BigInt fn) = BigInt.FibonacciPair(n + 1);
            BigInt fm1 = fn1 - fn;
            Assert.Equal((n & 1) == 0 ? BigInt.One : BigInt.MinusOne, fm1 * fn1 - fn.Square());
            (BigInt ln, BigInt lnm1) = BigInt.LucasPair(n);
            Assert.Equal(fm1 + fn1, ln);
            Assert.Equal(fn + (fn - fm1), lnm1); // L(n-1) = F(n) + F(n-2)
        }
    }

    [Fact]
    public void Fibonacci_MatchesRecurrence()
    {
        BigInteger prev = 0, cur = 1; // F(0), F(1)
        for (ulong n = 1; n <= 3000; n++)
        {
            if (n < 600 || n % 97 == 0 || n % 128 <= 1)
            {
                TestHelpers.AssertEqual(cur, BigInt.Fibonacci(n), $"F({n})");
                (BigInt fn, BigInt fnm1) = BigInt.FibonacciPair(n);
                TestHelpers.AssertEqual(cur, fn, $"F2({n})");
                TestHelpers.AssertEqual(prev, fnm1, $"F2-1({n})");
                TestHelpers.AssertEqual(cur + 2 * prev, BigInt.Lucas(n), $"L({n})");
            }
            (prev, cur) = (cur, prev + cur);
        }
    }

    [Fact]
    public void Fibonacci_KnownValues()
    {
        Assert.Equal(BigInt.Zero, BigInt.Fibonacci(0));
        Assert.Equal(BigInt.One, BigInt.Fibonacci(1));
        Assert.Equal(new BigInt(55), BigInt.Fibonacci(10));
        Assert.Equal(BigInt.Parse("354224848179261915075"), BigInt.Fibonacci(100));
        (BigInt f10, BigInt f9) = BigInt.FibonacciPair(10);
        Assert.Equal(new BigInt(55), f10);
        Assert.Equal(new BigInt(34), f9);
        // Recurrence check at larger index.
        Assert.Equal(BigInt.Fibonacci(1000), BigInt.Fibonacci(999) + BigInt.Fibonacci(998));
    }

    [Fact]
    public void Lucas_KnownValues()
    {
        Assert.Equal(BigInt.Two, BigInt.Lucas(0));
        Assert.Equal(BigInt.One, BigInt.Lucas(1));
        Assert.Equal(new BigInt(123), BigInt.Lucas(10));
        (BigInt l10, BigInt l9) = BigInt.LucasPair(10);
        Assert.Equal(new BigInt(123), l10);
        Assert.Equal(new BigInt(76), l9);
        // L(n) = F(n-1) + F(n+1)
        Assert.Equal(BigInt.Lucas(500), BigInt.Fibonacci(499) + BigInt.Fibonacci(501));
    }

    [Fact]
    public void RemoveFactor_And_CountDigits()
    {
        BigInt v = BigInt.Pow(3, 7) * 250;
        Assert.Equal(new BigInt(250), BigInt.RemoveFactor(v, 3, out ulong c3));
        Assert.Equal(7UL, c3);
        Assert.Equal(new BigInt(125) * BigInt.Pow(3, 7), BigInt.RemoveFactor(v, 2, out ulong c2));
        Assert.Equal(1UL, c2);
        BigInt.RemoveFactor(new BigInt(7), 5, out ulong c5);
        Assert.Equal(0UL, c5);

        Assert.Equal(3, new BigInt(255).CountDigits(10));
        Assert.Equal(2, new BigInt(255).CountDigits(16));
        Assert.Equal(8, new BigInt(255).CountDigits(2));
        Assert.Equal(1, BigInt.Zero.CountDigits(10));
    }

    [Fact]
    public void Lcm_Works()
    {
        Assert.Equal(new BigInt(12), BigInt.Lcm(4, 6));
        Assert.Equal(BigInt.Zero, BigInt.Lcm(0, 5));
        var rng = new Random(82);
        for (int i = 0; i < 100; i++)
        {
            var (a, _) = TestHelpers.RandomPair(rng, 4);
            var (b, _) = TestHelpers.RandomPair(rng, 4);
            if (a.IsZero || b.IsZero)
            {
                continue;
            }
            BigInt l = BigInt.Lcm(a, b);
            Assert.Equal(a.Abs() * b.Abs(), l * BigInt.Gcd(a, b));
        }
    }
}
