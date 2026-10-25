using System.Diagnostics;

namespace Sdcb.Arithmetic;

/// <summary>
/// An immutable arbitrary-precision rational number (the equivalent of GMP's mpq_t),
/// always stored in canonical form: gcd(numerator, denominator) == 1, denominator &gt; 0.
/// </summary>
public readonly partial struct BigRational :
    IEquatable<BigRational>,
    IComparable<BigRational>,
    IComparable,
    ISpanFormattable,
    ISpanParsable<BigRational>
{
    private readonly BigInt _num;
    private readonly BigInt _den; // canonical: > 0; default(BigRational) stores 0 here and means denominator 1

    /// <summary>The numerator (carries the sign).</summary>
    public BigInt Numerator => _num;

    /// <summary>The denominator (always positive).</summary>
    public BigInt Denominator => _den.IsZero ? BigInt.One : _den;

    private BigRational(BigInt num, BigInt den, bool canonical)
    {
        Debug.Assert(canonical && (den.Sign > 0 || (num.IsZero && den.IsZero)));
        _num = num;
        _den = den;
    }

    public BigRational(BigInt numerator) : this(numerator, BigInt.One, canonical: true) { }

    /// <summary>Creates numerator/denominator, canonicalizing (mpq_canonicalize).</summary>
    public BigRational(BigInt numerator, BigInt denominator)
    {
        if (denominator.IsZero)
        {
            throw new DivideByZeroException("Rational denominator cannot be zero.");
        }
        if (numerator.IsZero)
        {
            _num = default;
            _den = BigInt.One;
            return;
        }
        if (denominator.Sign < 0)
        {
            numerator = -numerator;
            denominator = -denominator;
        }
        BigInt g = BigInt.Gcd(numerator, denominator);
        if (!g.IsOne)
        {
            numerator = BigInt.DivExact(numerator, g);
            denominator = BigInt.DivExact(denominator, g);
        }
        _num = numerator;
        _den = denominator;
    }

    public BigRational(long numerator, long denominator) : this(new BigInt(numerator), new BigInt(denominator)) { }

    public BigRational(long value) : this(new BigInt(value)) { }

    /// <summary>Exact conversion from a finite double (mpq_set_d).</summary>
    public BigRational(double value)
    {
        if (double.IsNaN(value) || double.IsInfinity(value))
        {
            throw new OverflowException("Cannot convert NaN or Infinity to BigRational.");
        }
        if (value == 0)
        {
            _num = default;
            _den = BigInt.One;
            return;
        }
        ulong bits = BitConverter.DoubleToUInt64Bits(value);
        int sign = (bits >> 63) != 0 ? -1 : 1;
        int exp = (int)((bits >> 52) & 0x7FF);
        ulong mant = bits & 0xFFFFFFFFFFFFFUL;
        int e;
        if (exp == 0)
        {
            e = -1074; // subnormal
        }
        else
        {
            mant |= 1UL << 52;
            e = exp - 1075;
        }
        BigInt m = new(mant);
        if (sign < 0)
        {
            m = -m;
        }
        BigRational r = e >= 0
            ? new BigRational(m.ShiftLeft(e))
            : new BigRational(m, BigInt.One.ShiftLeft(-e));
        _num = r._num;
        _den = r._den;
    }

    internal static BigRational FromCanonical(BigInt num, BigInt den) => new(num, den, canonical: true);

    // ---------------------------------------------------------------
    // Well-known values / predicates
    // ---------------------------------------------------------------

    public static BigRational Zero => default;
    public static BigRational One => new(BigInt.One);

    public int Sign => _num.Sign;
    public bool IsZero => _num.IsZero;
    public bool IsInteger => Denominator.IsOne;

    // ---------------------------------------------------------------
    // Conversions
    // ---------------------------------------------------------------

    public static implicit operator BigRational(BigInt value) => new(value);
    public static implicit operator BigRational(long value) => new(new BigInt(value));
    public static implicit operator BigRational(int value) => new(new BigInt(value));
    public static explicit operator BigRational(double value) => new(value);

    /// <summary>Truncates toward zero (mpq → mpz).</summary>
    public static explicit operator BigInt(BigRational value) => value._num / value.Denominator;

    public static explicit operator double(BigRational value) => value.ToDouble();

    /// <summary>Converts to the nearest double (mpq_get_d), correctly rounded.</summary>
    public double ToDouble() => ToIeee(IeeeFormat.Binary64);

    /// <summary>Converts to the nearest float, correctly rounded (including subnormals).</summary>
    public float ToSingle() => (float)ToIeee(IeeeFormat.Binary32);

    /// <summary>Converts to the nearest Half, correctly rounded (including subnormals).</summary>
    public Half ToHalf() => (Half)ToIeee(IeeeFormat.Binary16);

    private double ToIeee(IeeeFormat format) =>
        _num.IsZero ? 0.0 : Real.RationalToIeee(_num.Abs(), Denominator, _num.Sign < 0, format);

    // ---------------------------------------------------------------
    // Arithmetic
    // ---------------------------------------------------------------

    public static BigRational operator -(BigRational v) => FromCanonical(-v._num, v._den);
    public static BigRational operator +(BigRational v) => v;

    public BigRational Abs() => FromCanonical(_num.Abs(), _den);

    /// <summary>Reciprocal (mpq_inv).</summary>
    public BigRational Inverse()
    {
        if (_num.IsZero)
        {
            throw new DivideByZeroException("Cannot invert zero.");
        }
        BigInt n = Denominator, d = _num;
        if (d.Sign < 0)
        {
            n = -n;
            d = -d;
        }
        return FromCanonical(n, d);
    }

    public static BigRational operator +(BigRational a, BigRational b)
    {
        BigInt ad = a.Denominator, bd = b.Denominator;
        BigInt g = BigInt.Gcd(ad, bd);
        if (g.IsOne)
        {
            return new BigRational(a._num * bd + b._num * ad, ad * bd);
        }
        BigInt adg = BigInt.DivExact(ad, g);
        BigInt bdg = BigInt.DivExact(bd, g);
        BigInt num = a._num * bdg + b._num * adg;
        BigInt g2 = BigInt.Gcd(num, g);
        return g2.IsOne
            ? FromCanonicalReduced(num, adg * bd)
            : FromCanonicalReduced(BigInt.DivExact(num, g2), adg * BigInt.DivExact(bd, g2));
    }

    private static BigRational FromCanonicalReduced(BigInt num, BigInt den)
    {
        if (num.IsZero)
        {
            return default;
        }
        Debug.Assert(BigInt.Gcd(num, den).IsOne);
        return FromCanonical(num, den);
    }

    public static BigRational operator -(BigRational a, BigRational b) => a + (-b);

    public static BigRational operator *(BigRational a, BigRational b)
    {
        if (a.IsZero || b.IsZero)
        {
            return default;
        }
        BigInt g1 = BigInt.Gcd(a._num, b.Denominator);
        BigInt g2 = BigInt.Gcd(b._num, a.Denominator);
        BigInt num = BigInt.DivExact(a._num, g1) * BigInt.DivExact(b._num, g2);
        BigInt den = BigInt.DivExact(a.Denominator, g2) * BigInt.DivExact(b.Denominator, g1);
        return FromCanonical(num, den);
    }

    public static BigRational operator /(BigRational a, BigRational b) => a * b.Inverse();

    /// <summary>Raises to an integer power (negative exponents invert).</summary>
    public static BigRational Pow(BigRational value, int exponent)
    {
        if (exponent == 0)
        {
            return One;
        }
        if (exponent < 0)
        {
            value = value.Inverse();
        }
        ulong e = (ulong)Math.Abs((long)exponent);
        return FromCanonical(BigInt.Pow(value._num, e), BigInt.Pow(value.Denominator, e));
    }

    public static BigRational operator *(BigRational a, BigInt b) => a * new BigRational(b);
    public static BigRational operator /(BigRational a, BigInt b) => a / new BigRational(b);

    /// <summary>Multiplies by 2^shift (mpq_mul_2exp / mpq_div_2exp for negative shift).</summary>
    public BigRational ShiftLeft(long shift)
    {
        if (IsZero || shift == 0)
        {
            return this;
        }
        if (shift > 0)
        {
            long denTwos = Denominator.IsEven ? Denominator.TrailingZeroCount() : 0;
            long fromDen = Math.Min(denTwos, shift);
            return FromCanonical(_num.ShiftLeft(shift - fromDen), Denominator.ShiftRightTruncate(fromDen));
        }
        long s = -shift;
        long numTwos = _num.IsEven ? _num.TrailingZeroCount() : 0;
        long fromNum = Math.Min(numTwos, s);
        return FromCanonical(_num.ShiftRightTruncate(fromNum), Denominator.ShiftLeft(s - fromNum));
    }

    // ---------------------------------------------------------------
    // Comparison
    // ---------------------------------------------------------------

    public int CompareTo(BigRational other)
    {
        int sa = Sign, sb = other.Sign;
        if (sa != sb)
        {
            return sa < sb ? -1 : 1;
        }
        if (sa == 0)
        {
            return 0;
        }
        return (_num * other.Denominator).CompareTo(other._num * Denominator);
    }

    int IComparable.CompareTo(object? obj) => obj switch
    {
        null => 1,
        BigRational r => CompareTo(r),
        _ => throw new ArgumentException("Object must be a BigRational.", nameof(obj)),
    };

    public bool Equals(BigRational other) => _num.Equals(other._num) && Denominator.Equals(other.Denominator);
    public override bool Equals(object? obj) => obj is BigRational r && Equals(r);
    public override int GetHashCode() => HashCode.Combine(_num, Denominator);

    public static bool operator ==(BigRational a, BigRational b) => a.Equals(b);
    public static bool operator !=(BigRational a, BigRational b) => !a.Equals(b);
    public static bool operator <(BigRational a, BigRational b) => a.CompareTo(b) < 0;
    public static bool operator >(BigRational a, BigRational b) => a.CompareTo(b) > 0;
    public static bool operator <=(BigRational a, BigRational b) => a.CompareTo(b) <= 0;
    public static bool operator >=(BigRational a, BigRational b) => a.CompareTo(b) >= 0;

    // ---------------------------------------------------------------
    // Formatting / parsing ("num/den", integers, or decimal like "1.25e3")
    // ---------------------------------------------------------------

    public override string ToString() => IsInteger ? _num.ToString() : $"{_num}/{Denominator}";

    /// <summary>Formats in the given base as "num/den" (mpq_get_str).</summary>
    public string ToString(int toBase) => IsInteger ? _num.ToString(toBase) : $"{_num.ToString(toBase)}/{Denominator.ToString(toBase)}";

    public string ToString(string? format, IFormatProvider? formatProvider) => ToString();

    public bool TryFormat(Span<char> destination, out int charsWritten, ReadOnlySpan<char> format, IFormatProvider? provider)
    {
        string s = ToString();
        if (s.Length > destination.Length)
        {
            charsWritten = 0;
            return false;
        }
        s.CopyTo(destination);
        charsWritten = s.Length;
        return true;
    }

    public static BigRational Parse(string s, IFormatProvider? provider = null) => Parse(s.AsSpan(), provider);

    public static BigRational Parse(ReadOnlySpan<char> s, IFormatProvider? provider = null)
    {
        if (!TryParse(s, provider, out BigRational result))
        {
            throw new FormatException("The value could not be parsed as a BigRational.");
        }
        return result;
    }

    public static bool TryParse(string? s, IFormatProvider? provider, out BigRational result) =>
        TryParse(s.AsSpan(), provider, out result);

    public static bool TryParse(ReadOnlySpan<char> s, IFormatProvider? provider, out BigRational result)
    {
        result = default;
        s = s.Trim();
        int slash = s.IndexOf('/');
        if (slash >= 0)
        {
            if (!BigInt.TryParse(s.Slice(0, slash), 10, out BigInt num) ||
                !BigInt.TryParse(s.Slice(slash + 1), 10, out BigInt den) ||
                den.IsZero)
            {
                return false;
            }
            result = new BigRational(num, den);
            return true;
        }
        // Decimal form: [sign]digits[.digits][e±exp]
        int ePos = s.IndexOfAny('e', 'E');
        long exp10 = 0;
        if (ePos >= 0)
        {
            if (!long.TryParse(s.Slice(ePos + 1), out exp10))
            {
                return false;
            }
            s = s.Slice(0, ePos);
        }
        int dot = s.IndexOf('.');
        if (dot >= 0)
        {
            int fracDigits = s.Length - dot - 1;
            Span<char> merged = s.Length <= 256 ? stackalloc char[s.Length - 1] : new char[s.Length - 1];
            s.Slice(0, dot).CopyTo(merged);
            s.Slice(dot + 1).CopyTo(merged.Slice(dot));
            exp10 -= fracDigits;
            s = new string(merged); // spans can't re-point at stackalloc portably here
        }
        if (!BigInt.TryParse(s, 10, out BigInt mant))
        {
            return false;
        }
        BigInt pow = BigInt.Pow(BigInt.Ten, (ulong)Math.Abs(exp10));
        result = exp10 >= 0 ? new BigRational(mant * pow) : new BigRational(mant, pow);
        return true;
    }
}
