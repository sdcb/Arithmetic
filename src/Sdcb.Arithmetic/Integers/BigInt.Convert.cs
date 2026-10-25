using System.Diagnostics;
using System.Numerics;
using Sdcb.Arithmetic.Core;

namespace Sdcb.Arithmetic;

public readonly partial struct BigInt
{
    private const string Digits36 = "0123456789abcdefghijklmnopqrstuvwxyz";
    private const string Digits62 = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";

    /// <summary>Largest power of each base that fits in a ulong, and its exponent.</summary>
    private static readonly (ulong BigPow, int Digits)[] BaseInfo = BuildBaseInfo();

    private static (ulong, int)[] BuildBaseInfo()
    {
        var info = new (ulong, int)[63];
        for (int b = 2; b <= 62; b++)
        {
            ulong p = 1;
            int k = 0;
            while (p <= ulong.MaxValue / (ulong)b)
            {
                p *= (ulong)b;
                k++;
            }
            info[b] = (p, k);
        }
        return info;
    }

    private const int FormatBasecaseLimbs = 30;

    // powers[i] = bigPow^(2^i) per base, shared across calls; entries above PowCacheMaxLimbs are
    // recomputed per call so huge conversions do not pin their tables.
    private static readonly BigInt[]?[] s_powerCache = new BigInt[]?[63];
    private const int PowCacheMaxLimbs = 1 << 14;

    /// <summary>Returns at least <paramref name="count"/> entries of bigPow^(2^i) for base b.</summary>
    private static BigInt[] PowerTable(int b, int count)
    {
        BigInt[]? cached = Volatile.Read(ref s_powerCache[b]);
        if (cached is not null && cached.Length >= count)
        {
            return cached;
        }
        BigInt[] table = new BigInt[count];
        int have;
        if (cached is not null)
        {
            Array.Copy(cached, table, cached.Length);
            have = cached.Length;
        }
        else
        {
            table[0] = new BigInt(BaseInfo[b].BigPow);
            have = 1;
        }
        for (int i = have; i < count; i++)
        {
            table[i] = table[i - 1].Square();
        }
        int keep = 0;
        while (keep < count && table[keep].LimbCount <= PowCacheMaxLimbs)
        {
            keep++;
        }
        if (keep > (cached?.Length ?? 0))
        {
            Volatile.Write(ref s_powerCache[b], keep == count ? table : table[..keep]);
        }
        return table;
    }
    private const int ParseBasecaseDigits = 1000;

    // ---------------------------------------------------------------
    // Formatting
    // ---------------------------------------------------------------

    /// <summary>Formats in base 10.</summary>
    public override string ToString() => ToString(10);

    /// <summary>
    /// Formats in the given base (2..62). Bases up to 36 use lowercase digits;
    /// bases 37..62 use the 0-9, A-Z, a-z alphabet (like GMP).
    /// </summary>
    public string ToString(int toBase)
    {
        if (toBase < 2 || toBase > 62)
        {
            throw new ArgumentOutOfRangeException(nameof(toBase), "Base must be in 2..62.");
        }
        if (_size == 0)
        {
            return "0";
        }
        string alphabet = toBase <= 36 ? Digits36 : Digits62;
        string digits = BitOperations.IsPow2((uint)toBase)
            ? FormatPow2(toBase, alphabet)
            : FormatGeneral(toBase, alphabet);
        return _size < 0 ? "-" + digits : digits;
    }

    public string ToString(string? format, IFormatProvider? formatProvider)
    {
        if (string.IsNullOrEmpty(format) || format is "G" or "g" or "D" or "d" or "R" or "r")
        {
            return ToString(10);
        }
        if (format is "X")
        {
            return ToString(16).ToUpperInvariant();
        }
        if (format is "x")
        {
            return ToString(16);
        }
        throw new FormatException($"Unsupported format string '{format}' for BigInt.");
    }

    public bool TryFormat(Span<char> destination, out int charsWritten, ReadOnlySpan<char> format, IFormatProvider? provider)
    {
        string s = ToString(format.IsEmpty ? null : format.ToString(), provider);
        if (s.Length > destination.Length)
        {
            charsWritten = 0;
            return false;
        }
        s.CopyTo(destination);
        charsWritten = s.Length;
        return true;
    }

    private string FormatPow2(int b, string alphabet)
    {
        int bitsPerDigit = BitOperations.TrailingZeroCount((uint)b);
        long bits = GetBitLength();
        long nd = (bits + bitsPerDigit - 1) / bitsPerDigit;
        ReadOnlySpan<ulong> m = Limbs;
        return string.Create((int)nd, (Value: this, bitsPerDigit, alphabet), static (dest, state) =>
        {
            ReadOnlySpan<ulong> m = state.Value.Limbs;
            int bpd = state.bitsPerDigit;
            ulong mask = (1UL << bpd) - 1;
            long bitPos = 0;
            for (int i = dest.Length - 1; i >= 0; i--)
            {
                int limb = (int)(bitPos >> 6);
                int ofs = (int)(bitPos & 63);
                ulong v = m[limb] >> ofs;
                int got = 64 - ofs;
                if (got < bpd && limb + 1 < m.Length)
                {
                    v |= m[limb + 1] << got;
                }
                dest[i] = state.alphabet[(int)(v & mask)];
                bitPos += bpd;
            }
        });
    }

    private string FormatGeneral(int b, string alphabet)
    {
        (ulong bigPow, int chunk) = BaseInfo[b];
        BigInt abs = Abs();
        long bits = GetBitLength();
        // Upper bound on digit count.
        long ub = (long)(bits / Math.Log2(b)) + 2;

        if (LimbCount <= FormatBasecaseLimbs)
        {
            char[] buf = new char[ub];
            int written = FormatBasecase(buf, abs, b, alphabet, bigPow, chunk, pad: false);
            return new string(buf, buf.Length - written, written);
        }
        else
        {
            // Smallest level with bigPow^(2^level) > |value|.
            int level = 0;
            while ((double)(1L << level) * Math.Log2(bigPow) < bits - 1)
            {
                level++;
            }
            BigInt[] powers = PowerTable(b, level + 2);
            while (powers[level].CompareTo(abs) <= 0)
            {
                level++;
                if (level >= powers.Length)
                {
                    powers = PowerTable(b, level + 1);
                }
            }
            while (level > 0 && powers[level - 1].CompareTo(abs) > 0)
            {
                level--;
            }
            int width = chunk << level;
            char[] full = System.Buffers.ArrayPool<char>.Shared.Rent(width);
            FormatRec(full.AsSpan(0, width), abs, level, powers, b, alphabet, bigPow, chunk);
            // Trim leading zeros.
            int start = 0;
            while (start < width - 1 && full[start] == '0')
            {
                start++;
            }
            string s = new(full, start, width - start);
            System.Buffers.ArrayPool<char>.Shared.Return(full);
            return s;
        }
    }

    /// <summary>Writes v zero-padded to fill dest (dest.Length == chunk &lt;&lt; level).</summary>
    private static void FormatRec(Span<char> dest, BigInt v, int level, BigInt[] powers, int b, string alphabet, ulong bigPow, int chunk)
    {
        if (level == 0)
        {
            WriteChunk(dest, v._size == 0 ? 0 : v._d![0], b, alphabet);
            return;
        }
        if (v.LimbCount <= FormatBasecaseLimbs)
        {
            FormatBasecase(dest, v, b, alphabet, bigPow, chunk, pad: true);
            return;
        }
        int half = dest.Length >> 1;
        DivRem(v, powers[level - 1], out BigInt q, out BigInt r);
        FormatRec(dest.Slice(0, half), q, level - 1, powers, b, alphabet, bigPow, chunk);
        FormatRec(dest.Slice(half), r, level - 1, powers, b, alphabet, bigPow, chunk);
    }

    /// <summary>
    /// Repeated DivRemScalar conversion writing right-aligned into dest.
    /// When pad is true fills the entire span (leading zeros); returns chars written.
    /// </summary>
    private static int FormatBasecase(Span<char> dest, BigInt v, int b, string alphabet, ulong bigPow, int chunk, bool pad)
    {
        int n = v.LimbCount;
        if (n == 0)
        {
            if (pad)
            {
                dest.Fill('0');
                return dest.Length;
            }
            dest[^1] = '0';
            return 1;
        }
        Span<ulong> t = n <= 64 ? stackalloc ulong[n] : new ulong[n];
        v.Limbs.CopyTo(t);
        int pos = dest.Length;
        while (n > 1)
        {
            ulong rem = Nat.DivRemScalar(t, t, n, bigPow);
            n = Nat.TrimmedLength(t, n);
            WriteChunk(dest.Slice(pos - chunk, chunk), rem, b, alphabet);
            pos -= chunk;
        }
        // Final limb: write without forced padding.
        ulong last = t[0];
        if (b == 10)
        {
            int nd = 1;
            for (ulong p10 = 10; nd < 20 && last >= p10; p10 *= 10)
            {
                nd++;
            }
            pos -= nd;
            WriteDecimal(dest.Slice(pos, nd), last);
        }
        else
        {
            do
            {
                dest[--pos] = alphabet[(int)(last % (ulong)b)];
                last /= (ulong)b;
            }
            while (last != 0);
        }
        if (pad)
        {
            dest.Slice(0, pos).Fill('0');
            return dest.Length;
        }
        return dest.Length - pos;
    }

    private static void WriteChunk(Span<char> dest, ulong value, int b, string alphabet)
    {
        if (b == 10)
        {
            WriteDecimal(dest, value);
            return;
        }
        for (int i = dest.Length - 1; i >= 0; i--)
        {
            dest[i] = alphabet[(int)(value % (ulong)b)];
            value /= (ulong)b;
        }
        Debug.Assert(value == 0);
    }

    private const string DigitPairs =
        "00010203040506070809101112131415161718192021222324252627282930313233343536373839" +
        "40414243444546474849505152535455565758596061626364656667686970717273747576777879" +
        "8081828384858687888990919293949596979899";

    /// <summary>Writes value as exactly dest.Length decimal digits (zero-padded), two at a time.</summary>
    private static void WriteDecimal(Span<char> dest, ulong value)
    {
        int i = dest.Length;
        // Constant divisors compile to multiplications; split into 32-bit halves when possible.
        while (value >= 1UL << 32 && i >= 2)
        {
            ulong q = value / 100;
            int pair = (int)(value - q * 100) * 2;
            dest[--i] = DigitPairs[pair + 1];
            dest[--i] = DigitPairs[pair];
            value = q;
        }
        uint v = (uint)value;
        while (i >= 2)
        {
            uint q = v / 100;
            int pair = (int)(v - q * 100) * 2;
            dest[--i] = DigitPairs[pair + 1];
            dest[--i] = DigitPairs[pair];
            v = q;
        }
        if (i == 1)
        {
            dest[0] = (char)('0' + v % 10);
            v /= 10;
        }
        Debug.Assert(v == 0);
    }

    // ---------------------------------------------------------------
    // Parsing
    // ---------------------------------------------------------------

    public static BigInt Parse(string s) => Parse(s.AsSpan(), 10);
    public static BigInt Parse(string s, int fromBase) => Parse(s.AsSpan(), fromBase);
    public static BigInt Parse(ReadOnlySpan<char> s) => Parse(s, 10);
    public static BigInt Parse(string s, IFormatProvider? provider) => Parse(s.AsSpan(), 10);
    public static BigInt Parse(ReadOnlySpan<char> s, IFormatProvider? provider) => Parse(s, 10);

    public static bool TryParse(string? s, out BigInt result) => TryParse(s.AsSpan(), 10, out result);
    public static bool TryParse(string? s, IFormatProvider? provider, out BigInt result) => TryParse(s.AsSpan(), 10, out result);
    public static bool TryParse(ReadOnlySpan<char> s, IFormatProvider? provider, out BigInt result) => TryParse(s, 10, out result);
    public static bool TryParse(ReadOnlySpan<char> s, out BigInt result) => TryParse(s, 10, out result);

    /// <summary>
    /// Parses in the given base (2..62), or 0 to auto-detect the base from a
    /// 0x/0b/0o prefix or a leading 0 for octal (like mpz_set_str).
    /// Whitespace and '_' separators between digits are ignored.
    /// </summary>
    public static BigInt Parse(ReadOnlySpan<char> s, int fromBase)
    {
        if (!TryParse(s, fromBase, out BigInt result))
        {
            throw new FormatException("The value could not be parsed as a BigInt.");
        }
        return result;
    }

    public static bool TryParse(ReadOnlySpan<char> s, int fromBase, out BigInt result)
    {
        result = default;
        if (fromBase != 0 && (fromBase < 2 || fromBase > 62))
        {
            throw new ArgumentOutOfRangeException(nameof(fromBase), "Base must be 0 (auto) or in 2..62.");
        }
        s = s.Trim();
        if (s.IsEmpty)
        {
            return false;
        }
        int sign = 1;
        if (s[0] == '+' || s[0] == '-')
        {
            sign = s[0] == '-' ? -1 : 1;
            s = s.Slice(1);
        }
        if (fromBase == 0)
        {
            if (s.StartsWith("0x", StringComparison.OrdinalIgnoreCase))
            {
                fromBase = 16;
                s = s.Slice(2);
            }
            else if (s.StartsWith("0b", StringComparison.OrdinalIgnoreCase))
            {
                fromBase = 2;
                s = s.Slice(2);
            }
            else if (s.StartsWith("0o", StringComparison.OrdinalIgnoreCase))
            {
                fromBase = 8;
                s = s.Slice(2);
            }
            else if (s.Length > 1 && s[0] == '0')
            {
                fromBase = 8;
                s = s.Slice(1);
            }
            else
            {
                fromBase = 10;
            }
        }
        if (s.IsEmpty)
        {
            return false;
        }

        // Convert characters to digit values.
        byte[] digits = new byte[s.Length];
        int nd = 0;
        foreach (char c in s)
        {
            if (c == '_' || char.IsWhiteSpace(c))
            {
                continue;
            }
            int v = DigitValue(c, fromBase);
            if (v < 0)
            {
                return false;
            }
            digits[nd++] = (byte)v;
        }
        if (nd == 0)
        {
            return false;
        }

        ReadOnlySpan<byte> dspan = digits.AsSpan(0, nd);
        BigInt mag = BitOperations.IsPow2((uint)fromBase)
            ? ParsePow2(dspan, fromBase)
            : ParseGeneral(dspan, fromBase);
        result = sign < 0 ? -mag : mag;
        return true;
    }

    private static int DigitValue(char c, int b)
    {
        int v;
        if (b <= 36)
        {
            // Case-insensitive.
            v = c switch
            {
                >= '0' and <= '9' => c - '0',
                >= 'a' and <= 'z' => c - 'a' + 10,
                >= 'A' and <= 'Z' => c - 'A' + 10,
                _ => -1,
            };
        }
        else
        {
            v = c switch
            {
                >= '0' and <= '9' => c - '0',
                >= 'A' and <= 'Z' => c - 'A' + 10,
                >= 'a' and <= 'z' => c - 'a' + 36,
                _ => -1,
            };
        }
        return v < b ? v : -1;
    }

    private static BigInt ParsePow2(ReadOnlySpan<byte> digits, int b)
    {
        int bpd = BitOperations.TrailingZeroCount((uint)b);
        long totalBits = (long)digits.Length * bpd;
        int len = (int)((totalBits + 63) >> 6);
        ulong[] d = new ulong[len];
        long bitPos = 0;
        for (int i = digits.Length - 1; i >= 0; i--)
        {
            ulong v = digits[i];
            int limb = (int)(bitPos >> 6);
            int ofs = (int)(bitPos & 63);
            d[limb] |= v << ofs;
            if (ofs + bpd > 64)
            {
                d[limb + 1] |= v >> (64 - ofs);
            }
            bitPos += bpd;
        }
        return FromMagnitude(d, len, 1);
    }

    private static BigInt ParseGeneral(ReadOnlySpan<byte> digits, int b)
    {
        // Strip leading zeros.
        int start = 0;
        while (start < digits.Length && digits[start] == 0)
        {
            start++;
        }
        digits = digits.Slice(start);
        if (digits.IsEmpty)
        {
            return default;
        }
        if (digits.Length <= ParseBasecaseDigits)
        {
            return ParseBasecase(digits, b);
        }

        (ulong bigPow, int chunk) = BaseInfo[b];
        // powers[i] = bigPow^(2^i); find level so that chunk << level >= digits.Length.
        int levels = 0;
        long width = chunk;
        while (width < digits.Length)
        {
            levels++;
            width <<= 1;
        }
        return ParseRec(digits, PowerTable(b, levels + 1), levels, b, chunk);
    }

    private static BigInt ParseRec(ReadOnlySpan<byte> digits, BigInt[] powers, int level, int b, int chunk)
    {
        if (digits.Length <= ParseBasecaseDigits)
        {
            return ParseBasecase(digits, b);
        }
        while (level > 0 && ((long)chunk << (level - 1)) >= digits.Length)
        {
            level--;
        }
        Debug.Assert(level > 0);
        int lowLen = chunk << (level - 1);
        ReadOnlySpan<byte> high = digits.Slice(0, digits.Length - lowLen);
        ReadOnlySpan<byte> low = digits.Slice(digits.Length - lowLen);
        BigInt hi = ParseRec(high, powers, level - 1, b, chunk);
        BigInt lo = ParseRec(low, powers, level - 1, b, chunk);
        return hi * powers[level - 1] + lo;
    }

    private static BigInt ParseBasecase(ReadOnlySpan<byte> digits, int b)
    {
        (ulong bigPow, int chunk) = BaseInfo[b];
        int estLimbs = (int)((digits.Length * Math.Log2(b)) / 64) + 2;
        ulong[] buf = new ulong[estLimbs];
        int len = 0;

        int i = 0;
        int firstChunk = digits.Length % chunk;
        if (firstChunk == 0)
        {
            firstChunk = chunk;
        }
        while (i < digits.Length)
        {
            int take = i == 0 ? firstChunk : chunk;
            ulong chunkVal = 0;
            for (int j = 0; j < take; j++)
            {
                chunkVal = chunkVal * (ulong)b + digits[i + j];
            }
            ulong mulBy = i == 0 ? 1UL : bigPow;
            if (take != chunk && i != 0)
            {
                Debug.Assert(false, "only the first chunk may be short");
            }
            if (i != 0)
            {
                ulong carry = Nat.MulScalar(buf.AsSpan(0, len), buf.AsSpan(0, len), len, mulBy);
                if (carry != 0)
                {
                    buf[len++] = carry;
                }
            }
            ulong c2 = Nat.AddScalarInPlace(buf, len, chunkVal);
            if (c2 != 0)
            {
                buf[len++] = c2;
            }
            if (len == 0 && chunkVal != 0)
            {
                buf[0] = chunkVal;
                len = 1;
            }
            i += take;
        }
        return FromMagnitude(buf, len, 1);
    }

    // ---------------------------------------------------------------
    // Byte import/export (magnitude, like mpz_import/mpz_export)
    // ---------------------------------------------------------------

    /// <summary>Exports the absolute value as bytes (no sign); empty array for zero.</summary>
    public byte[] ToMagnitudeBytes(bool bigEndian = false)
    {
        if (_size == 0)
        {
            return [];
        }
        long bits = GetBitLength();
        int nb = (int)((bits + 7) >> 3);
        byte[] bytes = new byte[nb];
        ReadOnlySpan<ulong> m = Limbs;
        for (int i = 0; i < nb; i++)
        {
            byte v = (byte)(m[i >> 3] >> ((i & 7) * 8));
            bytes[bigEndian ? nb - 1 - i : i] = v;
        }
        return bytes;
    }

    /// <summary>Imports a magnitude from bytes with the given sign (mpz_import).</summary>
    public static BigInt FromMagnitudeBytes(ReadOnlySpan<byte> bytes, int sign = 1, bool bigEndian = false)
    {
        int nb = bytes.Length;
        if (nb == 0)
        {
            return default;
        }
        int len = (nb + 7) >> 3;
        ulong[] d = new ulong[len];
        for (int i = 0; i < nb; i++)
        {
            byte v = bytes[bigEndian ? nb - 1 - i : i];
            d[i >> 3] |= (ulong)v << ((i & 7) * 8);
        }
        return FromMagnitude(d, len, sign < 0 ? -1 : 1);
    }
}
