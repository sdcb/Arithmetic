using System.Numerics;
using Xunit;

namespace Sdcb.Arithmetic.Tests.Generic;

/// <summary>Generic algorithms written only against the .NET generic-math interfaces.</summary>
public class GenericAlgorithmTests
{
    private static T Sum<T>(IEnumerable<T> values) where T : INumberBase<T>
    {
        T s = T.Zero;
        foreach (T v in values)
        {
            s += v;
        }
        return s;
    }

    private static T Horner<T>(T x, params int[] coefficients) where T : INumberBase<T>
    {
        T r = T.Zero;
        foreach (int c in coefficients)
        {
            r = r * x + T.CreateChecked(c);
        }
        return r;
    }

    private static T Gcd<T>(T a, T b) where T : IBinaryInteger<T>
    {
        a = T.Abs(a);
        b = T.Abs(b);
        while (!T.IsZero(b))
        {
            (a, b) = (b, a % b);
        }
        return a;
    }

    private static T NewtonSqrt<T>(T x) where T : IFloatingPointIeee754<T>
    {
        T g = x;
        T two = T.One + T.One;
        for (int i = 0; i < 400; i++)
        {
            T n = (g + x / g) / two;
            if (n == g)
            {
                break;
            }
            g = n;
        }
        return g;
    }

    private static T Factorial<T>(int n) where T : INumber<T>
    {
        T r = T.One;
        for (T i = T.One; i <= T.CreateChecked(n); i++)
        {
            r *= i;
        }
        return r;
    }

    [Fact]
    public void Sum_AllTypes()
    {
        var rng = new Random(1);
        long[] xs = Enumerable.Range(0, 200).Select(_ => rng.NextInt64(long.MinValue / 4, long.MaxValue / 4)).ToArray();
        BigInteger expected = Sum(xs.Select(x => (BigInteger)x));
        Assert.Equal(expected, (BigInteger)Sum(xs.Select(x => (BigInt)x)));
        Assert.Equal((BigRational)expected, Sum(xs.Select(x => (BigRational)x)));
        Assert.True((BigFloat)expected == Sum(xs.Select(x => BigFloat.FromBigInt(x, 256))));

        // Harmonic number H_20 exactly.
        BigRational h = Sum(Enumerable.Range(1, 20).Select(k => new BigRational(1, k)));
        Assert.Equal(new BigRational(new BigInt(55835135), new BigInt(15519504)), h);
    }

    [Fact]
    public void Horner_AllTypes()
    {
        int[] coeffs = [3, -7, 0, 11, -2, 5];
        for (int x = -20; x <= 20; x++)
        {
            BigInteger expected = Horner((BigInteger)x, coeffs);
            Assert.Equal(expected, (BigInteger)Horner((BigInt)x, coeffs));
            Assert.Equal((BigRational)expected, Horner((BigRational)x, coeffs));
            Assert.True((BigFloat)expected == Horner((BigFloat)x, coeffs));
        }
        Assert.Equal(new BigRational(205, 32), Horner(new BigRational(1, 2), coeffs));
        Assert.Equal(205.0 / 32, Horner((BigFloat)0.5, coeffs).ToDouble());
    }

    [Fact]
    public void Gcd_GenericMatchesBuiltIn()
    {
        var rng = new Random(2);
        for (int i = 0; i < 200; i++)
        {
            BigInt a = TestHelpers.RandomPair(rng, 10).A;
            BigInt b = TestHelpers.RandomPair(rng, 10).A;
            Assert.Equal(BigInt.Gcd(a, b), Gcd(a, b));
        }
    }

    [Fact]
    public void Factorial_GenericMatchesBuiltIn()
    {
        Assert.Equal(BigInt.Factorial(40), Factorial<BigInt>(40));
        Assert.Equal(new BigRational(BigInt.Factorial(25)), Factorial<BigRational>(25));
        Assert.True(BigFloat.FromBigInt(BigInt.Factorial(30)) == Factorial<BigFloat>(30));
    }

    [Fact]
    public void NewtonSqrt_BigFloatConvergesToSqrt()
    {
        // Tests run in parallel, so stay at the (unchanged) DefaultPrecision.
        int p = BigFloat.DefaultPrecision;
        foreach (int v in new[] { 2, 3, 10, 12345, 1 << 20 })
        {
            BigFloat x = BigFloat.FromBigInt(v, p);
            BigFloat g = NewtonSqrt(x);
            BigFloat exact = BigFloat.Sqrt(x);
            Assert.True((g - exact).Abs() <= BigFloat.ScaleB(exact, 3 - p), $"sqrt({v})");
        }
        Assert.Equal(Math.Sqrt(2), NewtonSqrt(2.0), 15);
    }

    [Fact]
    public void Interfaces_AreImplemented()
    {
        Assert.True(typeof(IBinaryInteger<BigInt>).IsAssignableFrom(typeof(BigInt)));
        Assert.True(typeof(ISignedNumber<BigInt>).IsAssignableFrom(typeof(BigInt)));
        Assert.True(typeof(INumber<BigRational>).IsAssignableFrom(typeof(BigRational)));
        Assert.True(typeof(ISignedNumber<BigRational>).IsAssignableFrom(typeof(BigRational)));
        Assert.True(typeof(IBinaryFloatingPointIeee754<BigFloat>).IsAssignableFrom(typeof(BigFloat)));
        Assert.True(typeof(IUtf8SpanFormattable).IsAssignableFrom(typeof(BigFloat)));
        Assert.True(typeof(IUtf8SpanParsable<BigInt>).IsAssignableFrom(typeof(BigInt)));
    }
}
