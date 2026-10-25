using System.Diagnostics;
using Sdcb.Arithmetic.Core;

namespace Sdcb.Arithmetic;

public readonly partial struct Real
{
    // ---------------------------------------------------------------
    // Sign operations (exact)
    // ---------------------------------------------------------------

    public static Real operator -(Real v) => new(v._man, v._exp, v._prec, v._kind, !v.IsNeg);
    public static Real operator +(Real v) => v;

    /// <summary>Absolute value (mpfr_abs; exact here since precision is kept).</summary>
    public Real Abs() => new(_man, _exp, _prec, _kind, negative: false);

    /// <summary>Returns a value with the magnitude of this and the sign of source (mpfr_copysign).</summary>
    public Real CopySign(Real source) => new(_man, _exp, _prec, _kind, source.IsNeg);

    // ---------------------------------------------------------------
    // Addition / subtraction
    // ---------------------------------------------------------------

    public static Real operator +(Real a, Real b) => Add(a, b, 0, RealRounding.NearestEven);
    public static Real operator -(Real a, Real b) => Subtract(a, b, 0, RealRounding.NearestEven);

    public static Real Add(Real a, Real b, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
        => Add(a, b, precision, rounding, out _);

    public static Real Add(Real a, Real b, int precision, RealRounding rounding, out int ternary)
        => AddCore(a, b, negateB: false, ResolvePrec(a, b, precision), rounding, out ternary);

    public static Real Subtract(Real a, Real b, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
        => Subtract(a, b, precision, rounding, out _);

    public static Real Subtract(Real a, Real b, int precision, RealRounding rounding, out int ternary)
        => AddCore(a, b, negateB: true, ResolvePrec(a, b, precision), rounding, out ternary);

    private static int ResolvePrec(Real a, Real b, int precision)
        => precision != 0 ? precision : Math.Max(a.Precision, b.Precision);

    private static Real AddCore(Real a, Real b, bool negateB, int prec, RealRounding rnd, out int ternary)
    {
        ternary = 0;
        bool bNeg = b.IsNeg ^ negateB;

        if (a._kind == Kind.NaN || b._kind == Kind.NaN)
        {
            return MakeNaN(prec);
        }
        if (a._kind == Kind.Inf || b._kind == Kind.Inf)
        {
            if (a._kind == Kind.Inf && b._kind == Kind.Inf)
            {
                return a.IsNeg == bNeg ? MakeInf(a.IsNeg, prec) : MakeNaN(prec);
            }
            return a._kind == Kind.Inf ? MakeInf(a.IsNeg, prec) : MakeInf(bNeg, prec);
        }
        if (a._kind == Kind.Zero && b._kind == Kind.Zero)
        {
            // (+0) + (-0) = +0 except toward -inf; (-0) + (-0) = -0.
            bool neg = a.IsNeg && bNeg || (a.IsNeg != bNeg && rnd == RealRounding.TowardNegative);
            return MakeZero(neg, prec);
        }
        if (a._kind == Kind.Zero)
        {
            return Round(b._man, b._exp, sticky: false, bNeg, prec, rnd, out ternary);
        }
        if (b._kind == Kind.Zero)
        {
            return Round(a._man, a._exp, sticky: false, a.IsNeg, prec, rnd, out ternary);
        }

        // Order so a has the larger magnitude msb.
        bool aNeg = a.IsNeg;
        if (a.BinaryExponent < b.BinaryExponent)
        {
            (a, b) = (b, a);
            (aNeg, bNeg) = (bNeg, aNeg);
        }

        long msbB = b.BinaryExponent;
        // Extend a so it has >= prec+3 bits; if b lies entirely below that lsb, it only
        // decides the direction of a sticky unit (a ± θ ulp with 0 < θ < 1).
        long sh = Math.Max(3, prec + 3 - a._man.GetBitLength());
        if (msbB < a._exp - sh)
        {
            return AddStickyRound(a._man, sh, a._exp - sh, subtract: aNeg != bNeg, aNeg, prec, rnd, out ternary);
        }

        // Small significands: exact 128-bit alignment and sum.
        if (prec <= 126 && a._man.GetBitLength() <= 64 && b._man.GetBitLength() <= 64 && Math.Abs(a._exp - b._exp) <= 63)
        {
            long f0 = Math.Min(a._exp, b._exp);
            UInt128 xa = (UInt128)(ulong)a._man << (int)(a._exp - f0);
            UInt128 xb = (UInt128)(ulong)b._man << (int)(b._exp - f0);
            if (aNeg == bNeg)
            {
                return RoundSmall(xa + xb, f0, sticky: false, aNeg, prec, rnd, out ternary);
            }
            if (xa == xb)
            {
                return MakeZero(rnd == RealRounding.TowardNegative, prec);
            }
            return xa > xb
                ? RoundSmall(xa - xb, f0, sticky: false, aNeg, prec, rnd, out ternary)
                : RoundSmall(xb - xa, f0, sticky: false, bNeg, prec, rnd, out ternary);
        }

        // Exact alignment in scratch.
        return AddRound(a._man, a._exp, aNeg, b._man, b._exp, bNeg, prec, rnd, out ternary);
    }

    // ---------------------------------------------------------------
    // Multiplication
    // ---------------------------------------------------------------

    public static Real operator *(Real a, Real b) => Multiply(a, b, 0, RealRounding.NearestEven);

    public static Real Multiply(Real a, Real b, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
        => Multiply(a, b, precision, rounding, out _);

    public static Real Multiply(Real a, Real b, int precision, RealRounding rounding, out int ternary)
    {
        int prec = ResolvePrec(a, b, precision);
        ternary = 0;
        if (a._kind == Kind.NaN || b._kind == Kind.NaN)
        {
            return MakeNaN(prec);
        }
        bool neg = a.IsNeg ^ b.IsNeg;
        if (a._kind == Kind.Inf || b._kind == Kind.Inf)
        {
            if (a._kind == Kind.Zero || b._kind == Kind.Zero)
            {
                return MakeNaN(prec); // 0 * inf
            }
            return MakeInf(neg, prec);
        }
        if (a._kind == Kind.Zero || b._kind == Kind.Zero)
        {
            return MakeZero(neg, prec);
        }
        if (prec <= 126 && a._man.GetBitLength() <= 64 && b._man.GetBitLength() <= 64)
        {
            return RoundSmall((UInt128)(ulong)a._man * (ulong)b._man, a._exp + b._exp, sticky: false, neg, prec, rounding, out ternary);
        }
        return MulRound(a._man, b._man, a._exp + b._exp, neg, prec, rounding, out ternary);
    }

    /// <summary>Squares the value (mpfr_sqr).</summary>
    public static Real Square(Real a, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        if (a._kind != Kind.Normal)
        {
            return Multiply(a, a, precision, rounding, out _);
        }
        int prec = precision != 0 ? precision : a.Precision;
        return MulRound(a._man, a._man, 2 * a._exp, negative: false, prec, rounding, out _);
    }

    // ---------------------------------------------------------------
    // Division
    // ---------------------------------------------------------------

    public static Real operator /(Real a, Real b) => Divide(a, b, 0, RealRounding.NearestEven);

    public static Real Divide(Real a, Real b, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
        => Divide(a, b, precision, rounding, out _);

    public static Real Divide(Real a, Real b, int precision, RealRounding rounding, out int ternary)
    {
        int prec = ResolvePrec(a, b, precision);
        ternary = 0;
        if (a._kind == Kind.NaN || b._kind == Kind.NaN)
        {
            return MakeNaN(prec);
        }
        bool neg = a.IsNeg ^ b.IsNeg;
        if (a._kind == Kind.Inf)
        {
            return b._kind == Kind.Inf ? MakeNaN(prec) : MakeInf(neg, prec);
        }
        if (b._kind == Kind.Inf)
        {
            return MakeZero(neg, prec);
        }
        if (b._kind == Kind.Zero)
        {
            return a._kind == Kind.Zero ? MakeNaN(prec) : MakeInf(neg, prec);
        }
        if (a._kind == Kind.Zero)
        {
            return MakeZero(neg, prec);
        }

        long bitsA = a._man.GetBitLength(), bitsB = b._man.GetBitLength();
        long s = prec + 3 + bitsB - bitsA;
        if (s < 0)
        {
            s = 0;
        }
        if (prec <= 60 && bitsA <= 64 && bitsB <= 64)
        {
            // Quotient has prec+3 or prec+4 <= 64 bits: one 128/64 division.
            // The quotient fits a limb, so hi < d and one hardware 128/64 division suffices.
            UInt128 n = (UInt128)(ulong)a._man << (int)s;
            ulong d = (ulong)b._man;
            ulong q = Nat.DivRem128((ulong)(n >> 64), (ulong)n, d, out ulong rem);
            return RoundSmall(q, a._exp - b._exp - s, sticky: rem != 0, neg, prec, rounding, out ternary);
        }
        return DivRound(a._man, s, b._man, a._exp - b._exp - s, neg, prec, rounding, out ternary);
    }

    // ---------------------------------------------------------------
    // Square root / reciprocal square root
    // ---------------------------------------------------------------

    /// <summary>Correctly-rounded square root (mpfr_sqrt).</summary>
    public static Real Sqrt(Real a, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = precision != 0 ? precision : a.Precision;
        switch (a._kind)
        {
            case Kind.NaN:
                return MakeNaN(prec);
            case Kind.Zero:
                return MakeZero(a.IsNeg, prec); // sqrt(±0) = ±0
            case Kind.Inf:
                return a.IsNeg ? MakeNaN(prec) : MakeInf(negative: false, prec);
        }
        if (a.IsNeg)
        {
            return MakeNaN(prec);
        }
        long need = 2L * (prec + 2);
        long s = Math.Max(0, need - a._man.GetBitLength());
        if (((a._exp - s) & 1) != 0)
        {
            s++;
        }
        if (prec <= 61 && a._man.GetBitLength() <= 64 && a._man.GetBitLength() + s <= 127)
        {
            UInt128 n = (UInt128)(ulong)a._man << (int)s;
            ulong root = Isqrt128(n);
            return RoundSmall(root, (a._exp - s) / 2, sticky: (UInt128)root * root != n, negative: false, prec, rounding, out _);
        }
        return SqrtRound(a._man, s, (a._exp - s) / 2, prec, rounding, out _);
    }

    /// <summary>floor(sqrt(n)) for n &lt; 2^128 (double estimate + exact integer correction).</summary>
    private static ulong Isqrt128(UInt128 n)
    {
        double d = Math.Sqrt((double)n);
        ulong s = d >= 18446744073709551615.0 ? ulong.MaxValue : (ulong)d;
        if (s > 0 && n >> 64 < s)
        {
            // One Newton step from the 52-bit estimate lands within one of the root.
            s = (ulong)(((UInt128)s + Nat.DivRem128((ulong)(n >> 64), (ulong)n, s, out _)) >> 1);
        }
        while ((UInt128)s * s > n)
        {
            s--;
        }
        while (s != ulong.MaxValue && (UInt128)(s + 1) * (s + 1) <= n)
        {
            s++;
        }
        return s;
    }

    /// <summary>Fused multiply-add a*b + c with a single rounding (mpfr_fma).</summary>
    public static Real FusedMultiplyAdd(Real a, Real b, Real c, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = precision != 0 ? precision : Math.Max(Math.Max(a.Precision, b.Precision), c.Precision);
        if (a._kind == Kind.Normal && b._kind == Kind.Normal)
        {
            BigInt pm = a._man * b._man;
            var product = new Real(pm, a._exp + b._exp, (int)pm.GetBitLength(), Kind.Normal, a.IsNeg ^ b.IsNeg);
            return AddCore(product, c, negateB: false, prec, rounding, out _);
        }
        // Fall back through the (special-value-aware) primitive ops.
        Real p = Multiply(a, b, prec + 8, RealRounding.TowardZero);
        return AddCore(p, c, negateB: false, prec, rounding, out _);
    }

    /// <summary>Fused multiply-subtract a*b - c (mpfr_fms).</summary>
    public static Real FusedMultiplySubtract(Real a, Real b, Real c, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
        => FusedMultiplyAdd(a, b, -c, precision, rounding);

    // ---------------------------------------------------------------
    // Comparison
    // ---------------------------------------------------------------

    /// <summary>IEEE-style equality (NaN != NaN, -0 == +0).</summary>
    public static bool operator ==(Real a, Real b) => a.CompareIeee(b) == 0;
    public static bool operator !=(Real a, Real b) => a.CompareIeee(b) != 0;
    public static bool operator <(Real a, Real b) => a.CompareIeee(b) is < 0 and > int.MinValue;
    public static bool operator >(Real a, Real b) => a.CompareIeee(b) is > 0 and < int.MaxValue;
    public static bool operator <=(Real a, Real b) => a.CompareIeee(b) is <= 0 and > int.MinValue;
    public static bool operator >=(Real a, Real b) => a.CompareIeee(b) is >= 0 and < int.MaxValue;

    /// <summary>Returns int.MinValue/int.MaxValue markers for unordered (NaN) comparisons.</summary>
    private int CompareIeee(Real other)
    {
        if (_kind == Kind.NaN || other._kind == Kind.NaN)
        {
            return int.MinValue; // unordered: makes <, >, ==, <= and >= all false
        }
        return CompareValues(this, other);
    }

    /// <summary>Numeric comparison; both operands must not be NaN.</summary>
    internal static int CompareValues(Real a, Real b)
    {
        int sa = a._kind == Kind.Zero ? 0 : a.IsNeg ? -1 : 1;
        int sb = b._kind == Kind.Zero ? 0 : b.IsNeg ? -1 : 1;
        if (sa != sb)
        {
            return sa < sb ? -1 : 1;
        }
        if (sa == 0)
        {
            return 0;
        }
        if (a._kind == Kind.Inf || b._kind == Kind.Inf)
        {
            int c = a._kind == Kind.Inf ? (b._kind == Kind.Inf ? 0 : 1) : -1;
            return sa > 0 ? c : -c;
        }
        // Both normal, same nonzero sign.
        long pa = a.BinaryExponent, pb = b.BinaryExponent;
        if (pa != pb)
        {
            int c = pa < pb ? -1 : 1;
            return sa > 0 ? c : -c;
        }
        long d = a._exp - b._exp;
        BigInt ma = a._man, mb = b._man;
        if (d > 0)
        {
            ma = ma.ShiftLeft(d);
        }
        else if (d < 0)
        {
            mb = mb.ShiftLeft(-d);
        }
        int cm = ma.CompareTo(mb);
        return sa > 0 ? cm : -cm;
    }

    /// <summary>Total-order comparison like double.CompareTo: NaN first, then -inf .. +inf; -0 &lt; +0.</summary>
    public int CompareTo(Real other)
    {
        bool an = _kind == Kind.NaN, bn = other._kind == Kind.NaN;
        if (an || bn)
        {
            return an == bn ? 0 : an ? -1 : 1;
        }
        int c = CompareValues(this, other);
        if (c != 0)
        {
            return c;
        }
        if (_kind == Kind.Zero && other._kind == Kind.Zero && IsNeg != other.IsNeg)
        {
            return IsNeg ? -1 : 1;
        }
        return 0;
    }

    int IComparable.CompareTo(object? obj) => obj switch
    {
        null => 1,
        Real r => CompareTo(r),
        _ => throw new ArgumentException("Object must be a Real.", nameof(obj)),
    };

    /// <summary>Bit-exact equality including NaN == NaN and -0 != +0 (unlike operator ==).</summary>
    public bool Equals(Real other)
    {
        if (_kind != other._kind || IsNeg != other.IsNeg)
        {
            return false;
        }
        if (_kind != Kind.Normal)
        {
            return true;
        }
        return _exp == other._exp && _man.Equals(other._man);
    }

    public override bool Equals(object? obj) => obj is Real r && Equals(r);

    public override int GetHashCode() => _kind switch
    {
        Kind.NaN => int.MinValue,
        Kind.Inf => IsNeg ? -2 : 2,
        Kind.Zero => 0,
        _ => HashCode.Combine(_man, _exp, IsNeg),
    };

    /// <summary>Smaller of two values; NaN if either is NaN unless the other is a number (mpfr_min semantics).</summary>
    public static Real Min(Real a, Real b)
    {
        if (a._kind == Kind.NaN)
        {
            return b._kind == Kind.NaN ? a : b;
        }
        if (b._kind == Kind.NaN)
        {
            return a;
        }
        int c = CompareValues(a, b);
        if (c != 0)
        {
            return c < 0 ? a : b;
        }
        // -0 preferred over +0.
        return a.IsNeg ? a : b;
    }

    /// <summary>Larger of two values (mpfr_max semantics).</summary>
    public static Real Max(Real a, Real b)
    {
        if (a._kind == Kind.NaN)
        {
            return b._kind == Kind.NaN ? a : b;
        }
        if (b._kind == Kind.NaN)
        {
            return a;
        }
        int c = CompareValues(a, b);
        if (c != 0)
        {
            return c > 0 ? a : b;
        }
        return a.IsNeg ? b : a;
    }

    /// <summary>Positive difference: max(a-b, +0) (mpfr_dim).</summary>
    public static Real PositiveDifference(Real a, Real b, int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        if (a._kind == Kind.NaN || b._kind == Kind.NaN)
        {
            return MakeNaN(ResolvePrec(a, b, precision));
        }
        Real d = Subtract(a, b, precision, rounding);
        return d._kind != Kind.NaN && (d._kind == Kind.Zero || d.IsNeg) ? MakeZero(negative: false, d.Precision) : d;
    }
}
