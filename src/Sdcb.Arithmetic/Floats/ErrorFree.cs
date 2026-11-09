using System.Runtime.CompilerServices;

namespace Sdcb.Arithmetic.Core;

/// <summary>
/// Error-free transformations on <see cref="double"/> shared by DoubleDouble and
/// QuadDouble: each returns s and the exact rounding error e with a + b = s + e
/// (resp. a * b = p + e), barring overflow.
/// </summary>
internal static class ErrorFree
{
    /// <summary>Knuth's TwoSum: s = fl(a + b), e = (a + b) - s exactly, for any a, b.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double TwoSum(double a, double b, out double e)
    {
        double s = a + b;
        double bb = s - a;
        e = (a - (s - bb)) + (b - bb);
        return s;
    }

    /// <summary>Dekker's FastTwoSum: requires |a| &gt;= |b| (or a == 0).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double FastTwoSum(double a, double b, out double e)
    {
        double s = a + b;
        e = b - (s - a);
        return s;
    }

    /// <summary>TwoProd via FMA: p = fl(a * b), e = a * b - p exactly.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double TwoProd(double a, double b, out double e)
    {
        double p = a * b;
        e = Math.FusedMultiplyAdd(a, b, -p);
        return p;
    }

    /// <summary>Exact square: p = fl(a * a), e = a * a - p.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double TwoSquare(double a, out double e)
    {
        double p = a * a;
        e = Math.FusedMultiplyAdd(a, a, -p);
        return p;
    }

    /// <summary>
    /// Renormalizes four overlapping components (decreasing magnitude, roughly) plus a fifth
    /// carry into a non-overlapping expansion (QD library's renorm with five inputs).
    /// </summary>
    public static void Renormalize(ref double c0, ref double c1, ref double c2, ref double c3, ref double c4)
    {
        if (double.IsInfinity(c0))
        {
            return;
        }

        double s0, s1, s2 = 0.0, s3 = 0.0;

        s0 = FastTwoSum(c3, c4, out c4);
        s0 = FastTwoSum(c2, s0, out c3);
        s0 = FastTwoSum(c1, s0, out c2);
        c0 = FastTwoSum(c0, s0, out c1);

        s0 = c0;
        s1 = c1;

        s0 = FastTwoSum(c0, c1, out s1);
        if (s1 != 0.0)
        {
            s1 = FastTwoSum(s1, c2, out s2);
            if (s2 != 0.0)
            {
                s2 = FastTwoSum(s2, c3, out s3);
                if (s3 != 0.0)
                {
                    s3 += c4;
                }
                else
                {
                    s2 = FastTwoSum(s2, c4, out s3);
                }
            }
            else
            {
                s1 = FastTwoSum(s1, c3, out s2);
                if (s2 != 0.0)
                {
                    s2 = FastTwoSum(s2, c4, out s3);
                }
                else
                {
                    s1 = FastTwoSum(s1, c4, out s2);
                }
            }
        }
        else
        {
            s0 = FastTwoSum(s0, c2, out s1);
            if (s1 != 0.0)
            {
                s1 = FastTwoSum(s1, c3, out s2);
                if (s2 != 0.0)
                {
                    s2 = FastTwoSum(s2, c4, out s3);
                }
                else
                {
                    s1 = FastTwoSum(s1, c4, out s2);
                }
            }
            else
            {
                s0 = FastTwoSum(s0, c3, out s1);
                if (s1 != 0.0)
                {
                    s1 = FastTwoSum(s1, c4, out s2);
                }
                else
                {
                    s0 = FastTwoSum(s0, c4, out s1);
                }
            }
        }

        c0 = s0;
        c1 = s1;
        c2 = s2;
        c3 = s3;
    }

    /// <summary>Four-input renormalization (Renormalize with a zero carry).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static void Renormalize(ref double c0, ref double c1, ref double c2, ref double c3)
    {
        double c4 = 0.0;
        Renormalize(ref c0, ref c1, ref c2, ref c3, ref c4);
    }
}
