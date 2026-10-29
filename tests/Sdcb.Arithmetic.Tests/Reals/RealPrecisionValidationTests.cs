using System.Reflection;
using Xunit;

namespace Sdcb.Arithmetic.Tests.Reals;

/// <summary>An explicit precision outside [2, int.MaxValue / 64] must be rejected; 0 means the default (issue #7).</summary>
public class RealPrecisionValidationTests
{
    private const int MaxPrecision = int.MaxValue / 64;

    [Theory]
    [InlineData(0, -1)]
    [InlineData(1, -1)]
    [InlineData(0, 1)]
    [InlineData(1, 1)]
    [InlineData(0, MaxPrecision + 1)]
    [InlineData(1, MaxPrecision + 1)]
    public void FromBigInt_InvalidPrecision_Throws(int m, int p)
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => Real.FromBigInt(new BigInt(m), p));
    }

    [Fact]
    public void FromBigInt_PrecisionZeroIsDefault()
    {
        Assert.True(Real.FromBigInt(BigInt.Zero, 0).IsZero);
        Assert.Equal(Real.DefaultPrecision, Real.FromBigInt(BigInt.One, 0).Precision);
        Assert.Equal(2, Real.FromBigInt(BigInt.One, 2).Precision);
    }

    // RoundSmall / RoundLimbs: padding, truncation and rounding carry must not wrap the exponent.
    public static IEnumerable<object[]> OverflowCases()
    {
        yield return [BigInt.One, long.MinValue, 2];
        yield return [new BigInt(4), long.MaxValue, 2];
        yield return [new BigInt(7), long.MaxValue - 1, 2];
        yield return [BigInt.One, long.MinValue, 127];
        yield return [BigInt.One << 128, long.MaxValue, 127];
        yield return [(BigInt.One << 128) - BigInt.One, long.MaxValue - 1, 127];
    }

    [Theory]
    [MemberData(nameof(OverflowCases))]
    public void FromBigInt_NormalizedExponentOverflow_Throws(BigInt m, long e, int p)
    {
        Assert.Throws<OverflowException>(() => Real.FromBigInt(m, e, p));
    }

    public static IEnumerable<object[]> PrecisionMethods() =>
        typeof(Real).GetMethods(BindingFlags.Public | BindingFlags.Static | BindingFlags.Instance | BindingFlags.DeclaredOnly)
            .Where(m => !m.IsGenericMethodDefinition && m.GetParameters().Any(p => p.Name == "precision" && p.ParameterType == typeof(int)))
            .Where(m => m.GetParameters().All(p => !p.ParameterType.IsByRefLike && !(p.ParameterType.IsByRef && p.ParameterType.GetElementType()!.IsByRefLike)))
            .Select(m => new object[] { m.ToString()! });

    [Theory]
    [MemberData(nameof(PrecisionMethods))]
    public void EveryPublicPrecisionParameter_IsValidated(string signature)
    {
        MethodInfo method = typeof(Real).GetMethods(BindingFlags.Public | BindingFlags.Static | BindingFlags.Instance | BindingFlags.DeclaredOnly)
            .Single(m => m.ToString() == signature);
        foreach (int bad in new[] { 1, -1, MaxPrecision + 1 })
        {
            object?[] args = method.GetParameters().Select(p => p.Name == "precision" ? bad : Sample(p.ParameterType)).ToArray();
            object? target = method.IsStatic ? null : Real.FromDouble(1.5, 64);
            TargetInvocationException ex = Assert.Throws<TargetInvocationException>(() => method.Invoke(target, args));
            Assert.IsType<ArgumentOutOfRangeException>(ex.InnerException);
        }
    }

    private static object? Sample(Type t)
    {
        if (t.IsByRef)
        {
            return null; // out parameter
        }
        if (t == typeof(Real)) return Real.FromDouble(1.5, 64);
        if (t == typeof(BigInt)) return new BigInt(3);
        if (t == typeof(BigRational)) return new BigRational(3, 2);
        if (t == typeof(double)) return 1.5;
        if (t == typeof(float)) return 1.5f;
        if (t == typeof(decimal)) return 1.5m;
        if (t == typeof(Half)) return (Half)1.5;
        if (t == typeof(string)) return "1.5";
        if (t == typeof(bool)) return false;
        if (t == typeof(Random)) return new Random(1);
        if (t == typeof(RealRounding)) return RealRounding.NearestEven;
        if (t == typeof(IFormatProvider)) return null;
        if (t == typeof(int)) return 3;
        if (t == typeof(long)) return 3L;
        if (t == typeof(uint)) return 3u;
        if (t == typeof(ulong)) return 3UL;
        if (t == typeof(UInt128)) return (UInt128)3;
        if (t == typeof(Int128)) return (Int128)3;
        if (t.IsEnum) return Enum.GetValues(t).GetValue(0);
        if (t.IsArray) return Array.CreateInstance(t.GetElementType()!, 0);
        throw new NotSupportedException(t.ToString());
    }
}
