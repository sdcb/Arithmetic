using System.Diagnostics;
using System.Numerics;
using System.Runtime.CompilerServices;
using Sdcb.Arithmetic.Core;

namespace Sdcb.Arithmetic;

/// <summary>
/// An immutable arbitrary-precision signed integer (the equivalent of GMP's mpz_t),
/// stored as a sign + little-endian array of 64-bit limbs.
/// </summary>
public readonly partial struct BigInt :
    IEquatable<BigInt>,
    IComparable<BigInt>,
    IComparable,
    ISpanFormattable,
    ISpanParsable<BigInt>
{
    /// <summary>Magnitude limbs, little-endian; exactly |_size| significant limbs. Null when zero.</summary>
    private readonly ulong[]? _d;

    /// <summary>Signed limb count: sign(value) * number of significant limbs.</summary>
    private readonly int _size;

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    internal BigInt(ulong[]? d, int size)
    {
        Debug.Assert(size == 0 || (d is not null && d.Length >= Math.Abs(size) && d[Math.Abs(size) - 1] != 0));
        _d = d;
        _size = size;
    }

    /// <summary>Creates a BigInt from a raw magnitude buffer, trimming high zero limbs.</summary>
    internal static BigInt FromMagnitude(ulong[] d, int len, int sign)
    {
        len = Nat.TrimmedLength(d, len);
        if (len == 0)
        {
            return default;
        }
        return new BigInt(d, sign < 0 ? -len : len);
    }

    internal int Size => _size;
    internal int LimbCount => Math.Abs(_size);
    internal ReadOnlySpan<ulong> Limbs => _d is null ? default : _d.AsSpan(0, Math.Abs(_size));

    // ---------------------------------------------------------------
    // Well-known values
    // ---------------------------------------------------------------

    public static BigInt Zero => default;
    public static BigInt One => new(OneLimb, 1);
    public static BigInt MinusOne => new(OneLimb, -1);
    public static BigInt Two => new(TwoLimb, 1);
    public static BigInt Ten => new(TenLimb, 1);

    private static readonly ulong[] OneLimb = [1];
    private static readonly ulong[] TwoLimb = [2];
    private static readonly ulong[] TenLimb = [10];

    // ---------------------------------------------------------------
    // Predicates / simple properties
    // ---------------------------------------------------------------

    /// <summary>-1, 0 or +1 according to the sign of the value.</summary>
    public int Sign => Math.Sign(_size);

    public bool IsZero => _size == 0;
    public bool IsOne => _size == 1 && _d![0] == 1;
    public bool IsMinusOne => _size == -1 && _d![0] == 1;
    public bool IsEven => _size == 0 || (_d![0] & 1) == 0;
    public bool IsOdd => !IsEven;
    public bool IsNegative => _size < 0;
    public bool IsPositive => _size > 0;

    /// <summary>True when the value is a power of two (and positive).</summary>
    public bool IsPowerOfTwo
    {
        get
        {
            if (_size <= 0)
            {
                return false;
            }
            ReadOnlySpan<ulong> m = Limbs;
            int n = m.Length;
            ulong top = m[n - 1];
            if ((top & (top - 1)) != 0)
            {
                return false;
            }
            for (int i = 0; i < n - 1; i++)
            {
                if (m[i] != 0)
                {
                    return false;
                }
            }
            return true;
        }
    }

    /// <summary>Number of bits in the absolute value (0 for zero). Same as mpz_sizeinbase(x, 2) for nonzero x.</summary>
    public long GetBitLength() => Nat.BitLength(Limbs, LimbCount);

    // ---------------------------------------------------------------
    // Construction from primitives
    // ---------------------------------------------------------------

    public BigInt(int value) : this((long)value) { }

    public BigInt(uint value) : this((ulong)value) { }

    public BigInt(long value)
    {
        if (value == 0)
        {
            _d = null;
            _size = 0;
        }
        else if (value > 0)
        {
            _d = [(ulong)value];
            _size = 1;
        }
        else
        {
            _d = [unchecked((ulong)(-value))];
            _size = -1;
        }
    }

    public BigInt(ulong value)
    {
        if (value == 0)
        {
            _d = null;
            _size = 0;
        }
        else
        {
            _d = [value];
            _size = 1;
        }
    }

    public BigInt(Int128 value)
    {
        if (value == 0)
        {
            _d = null;
            _size = 0;
            return;
        }
        int sign = 1;
        UInt128 mag;
        if (value < 0)
        {
            sign = -1;
            mag = unchecked((UInt128)(-value));
        }
        else
        {
            mag = (UInt128)value;
        }
        ulong hi = (ulong)(mag >> 64), lo = (ulong)mag;
        if (hi == 0)
        {
            _d = [lo];
            _size = sign;
        }
        else
        {
            _d = [lo, hi];
            _size = 2 * sign;
        }
    }

    public BigInt(UInt128 value)
    {
        ulong hi = (ulong)(value >> 64), lo = (ulong)value;
        if (hi != 0)
        {
            _d = [lo, hi];
            _size = 2;
        }
        else if (lo != 0)
        {
            _d = [lo];
            _size = 1;
        }
        else
        {
            _d = null;
            _size = 0;
        }
    }

    /// <summary>Creates a BigInt by truncating a finite double toward zero (like mpz_set_d).</summary>
    public BigInt(double value)
    {
        if (double.IsNaN(value) || double.IsInfinity(value))
        {
            throw new OverflowException("Cannot convert NaN or Infinity to BigInt.");
        }
        if (Math.Abs(value) < 1.0)
        {
            _d = null;
            _size = 0;
            return;
        }
        ulong bits = BitConverter.DoubleToUInt64Bits(value);
        int sign = (bits >> 63) != 0 ? -1 : 1;
        int exp = (int)((bits >> 52) & 0x7FF);
        ulong mant = bits & 0xFFFFFFFFFFFFFUL;
        Debug.Assert(exp != 0); // |value| >= 1 excludes subnormals
        mant |= 1UL << 52;
        // value = mant * 2^(exp - 1075)
        int e = exp - 1075;
        if (e <= 0)
        {
            ulong m = mant >> (-e);
            _d = [m];
            _size = sign;
            return;
        }
        int limbShift = e >> 6;
        int bitShift = e & 63;
        // mant has 53 bits; after bit shift may occupy 2 limbs.
        ulong lo = mant << bitShift;
        ulong hi = bitShift == 0 ? 0 : mant >> (64 - bitShift);
        int len = limbShift + (hi != 0 ? 2 : 1);
        ulong[] d = new ulong[len];
        d[limbShift] = lo;
        if (hi != 0)
        {
            d[limbShift + 1] = hi;
        }
        _d = d;
        _size = sign * len;
    }

    // ---------------------------------------------------------------
    // Implicit / explicit conversions
    // ---------------------------------------------------------------

    public static implicit operator BigInt(int value) => new(value);
    public static implicit operator BigInt(uint value) => new(value);
    public static implicit operator BigInt(long value) => new(value);
    public static implicit operator BigInt(ulong value) => new(value);
    public static implicit operator BigInt(Int128 value) => new(value);
    public static implicit operator BigInt(UInt128 value) => new(value);
    public static explicit operator BigInt(double value) => new(value);

    public static explicit operator long(BigInt value)
    {
        int n = value.LimbCount;
        if (n == 0)
        {
            return 0;
        }
        if (n > 1)
        {
            throw new OverflowException("Value does not fit in Int64.");
        }
        ulong m = value._d![0];
        if (value._size > 0)
        {
            if (m > long.MaxValue)
            {
                throw new OverflowException("Value does not fit in Int64.");
            }
            return (long)m;
        }
        if (m > (ulong)long.MaxValue + 1)
        {
            throw new OverflowException("Value does not fit in Int64.");
        }
        return unchecked(-(long)m);
    }

    public static explicit operator ulong(BigInt value)
    {
        int n = value.LimbCount;
        if (n == 0)
        {
            return 0;
        }
        if (n > 1 || value._size < 0)
        {
            throw new OverflowException("Value does not fit in UInt64.");
        }
        return value._d![0];
    }

    public static explicit operator int(BigInt value) => checked((int)(long)value);
    public static explicit operator uint(BigInt value) => checked((uint)(ulong)value);

    public static explicit operator Int128(BigInt value)
    {
        int n = value.LimbCount;
        if (n == 0)
        {
            return 0;
        }
        if (n > 2)
        {
            throw new OverflowException("Value does not fit in Int128.");
        }
        UInt128 mag = n == 1 ? value._d![0] : ((UInt128)value._d![1] << 64) | value._d![0];
        if (value._size > 0)
        {
            if (mag > (UInt128)Int128.MaxValue)
            {
                throw new OverflowException("Value does not fit in Int128.");
            }
            return (Int128)mag;
        }
        if (mag > (UInt128)Int128.MaxValue + 1)
        {
            throw new OverflowException("Value does not fit in Int128.");
        }
        return unchecked(-(Int128)mag);
    }

    public static explicit operator UInt128(BigInt value)
    {
        int n = value.LimbCount;
        if (n == 0)
        {
            return 0;
        }
        if (n > 2 || value._size < 0)
        {
            throw new OverflowException("Value does not fit in UInt128.");
        }
        return n == 1 ? value._d![0] : ((UInt128)value._d![1] << 64) | value._d![0];
    }

    public static explicit operator double(BigInt value) => value.ToDouble();

    /// <summary>Converts to double with correct round-to-nearest-even; returns ±Infinity on overflow.</summary>
    public double ToDouble()
    {
        int n = LimbCount;
        if (n == 0)
        {
            return 0.0;
        }
        ReadOnlySpan<ulong> m = Limbs;
        double result;
        long bits = Nat.BitLength(m, n);
        if (bits <= 53)
        {
            result = m[0]; // exact
        }
        else
        {
            // Extract the top 53 bits, then apply round-to-nearest-even using
            // the guard bit and a sticky bit computed from everything below.
            long shift = bits - 53;
            ulong m53 = ExtractBitsAt(m, shift);
            bool guard = ((ExtractBitsAt(m, shift - 1)) & 1) != 0;
            bool sticky = guard && Nat.TrailingZeroCount(m, n) < shift - 1;
            if (guard && (sticky || (m53 & 1) != 0))
            {
                m53++;
            }
            result = Math.ScaleB(m53, (int)Math.Min(shift, int.MaxValue));
        }
        return _size < 0 ? -result : result;
    }

    /// <summary>Extracts 64 bits of the magnitude starting at bit index shift (shift &gt;= 0).</summary>
    private static ulong ExtractBitsAt(ReadOnlySpan<ulong> m, long shift)
    {
        int limb = (int)(shift >> 6);
        int ofs = (int)(shift & 63);
        ulong lo = limb < m.Length ? m[limb] : 0;
        ulong hi = limb + 1 < m.Length ? m[limb + 1] : 0;
        return ofs == 0 ? lo : (lo >> ofs) | (hi << (64 - ofs));
    }

    /// <summary>Returns the low 64 bits of the magnitude with the sign applied (like mpz_get_si semantics, wrapping).</summary>
    public long ToInt64Unchecked()
    {
        if (_size == 0)
        {
            return 0;
        }
        ulong low = _d![0];
        return _size > 0 ? unchecked((long)low) : unchecked(-(long)low);
    }

    // ---------------------------------------------------------------
    // Comparison / equality
    // ---------------------------------------------------------------

    public int CompareTo(BigInt other)
    {
        if (_size != other._size)
        {
            return _size < other._size ? -1 : 1;
        }
        int c = Nat.Compare(Limbs, other.Limbs, LimbCount);
        return _size < 0 ? -c : c;
    }

    public int CompareTo(long other)
    {
        return CompareTo(new BigInt(other));
    }

    int IComparable.CompareTo(object? obj)
    {
        if (obj is null)
        {
            return 1;
        }
        if (obj is BigInt b)
        {
            return CompareTo(b);
        }
        throw new ArgumentException("Object must be a BigInt.", nameof(obj));
    }

    public bool Equals(BigInt other)
    {
        if (_size != other._size)
        {
            return false;
        }
        return Nat.Compare(Limbs, other.Limbs, LimbCount) == 0;
    }

    public override bool Equals(object? obj) => obj is BigInt other && Equals(other);

    public override int GetHashCode()
    {
        HashCode hash = default;
        hash.Add(_size);
        ReadOnlySpan<ulong> m = Limbs;
        for (int i = 0; i < m.Length; i++)
        {
            hash.Add(m[i]);
        }
        return hash.ToHashCode();
    }

    public static bool operator ==(BigInt left, BigInt right) => left.Equals(right);
    public static bool operator !=(BigInt left, BigInt right) => !left.Equals(right);
    public static bool operator <(BigInt left, BigInt right) => left.CompareTo(right) < 0;
    public static bool operator >(BigInt left, BigInt right) => left.CompareTo(right) > 0;
    public static bool operator <=(BigInt left, BigInt right) => left.CompareTo(right) <= 0;
    public static bool operator >=(BigInt left, BigInt right) => left.CompareTo(right) >= 0;

    /// <summary>Compares absolute values (mpz_cmpabs).</summary>
    public static int CompareAbs(BigInt a, BigInt b) => Nat.CompareNormalized(a.Limbs, a.LimbCount, b.Limbs, b.LimbCount);

    public static BigInt Min(BigInt a, BigInt b) => a.CompareTo(b) <= 0 ? a : b;
    public static BigInt Max(BigInt a, BigInt b) => a.CompareTo(b) >= 0 ? a : b;
}
