using System.Globalization;
using System.Text;
using Sdcb.Arithmetic.Core;

namespace Sdcb.Arithmetic;

// Formatting follows double's conventions:
//   ""/"R"/"G" (no precision): the shortest decimal that round-trips (Ryu's interval algorithm with
//   exact BigInt scaling), scientific notation when the decimal exponent is >= 33 or < -4 ("E+NN").
//   "G<n>": n significant digits, trailing zeros removed. "E<n>": n decimals, "E+ddd" exponent.
//   "F<n>" / "N<n>": n decimals (N with group separators). Precision rounding is exact, ties to even
//   (like double on .NET Core 3.0+). Parsing is correctly rounded from the exact decimal value.
public readonly partial struct Float128
{
    private const double Log10Of2 = 0.30102999566398119521;
    private const int GeneralDigits = 33; // FLT128_DIG: the threshold for scientific notation in "G"/"R"

    private static class Pow10Cache
    {
        public static readonly BigInt[] Values = Build(100);

        private static BigInt[] Build(int n)
        {
            BigInt[] t = new BigInt[n];
            t[0] = BigInt.One;
            for (int i = 1; i < n; i++)
            {
                t[i] = t[i - 1] * BigInt.Ten;
            }
            return t;
        }
    }

    private static BigInt Pow10(int k) =>
        k < Pow10Cache.Values.Length ? Pow10Cache.Values[k] : BigInt.Pow(BigInt.Ten, (ulong)k);

    private static readonly UInt128[] s_pow10u = BuildPow10U();

    private static UInt128[] BuildPow10U()
    {
        UInt128[] t = new UInt128[39];
        t[0] = UInt128.One;
        for (int i = 1; i < t.Length; i++)
        {
            t[i] = t[i - 1] * 10;
        }
        return t;
    }

    // ---------------------------------------------------------------
    // Shortest round-trip digits
    // ---------------------------------------------------------------

    /// <summary>floor(x * 2^e2 / 10^q) and whether the quotient is exact.</summary>
    private static UInt128 ScaleDown(UInt128 x, int e2, int q, out bool exact)
    {
        BigInt num = new(x);
        BigInt den = BigInt.One;
        if (q <= 0)
        {
            num *= Pow10(-q);
        }
        else
        {
            den = Pow10(q);
        }
        if (e2 >= 0)
        {
            num = num.ShiftLeft(e2);
        }
        else if (q <= 0)
        {
            exact = num.TrailingZeroCount() >= -e2;
            return (UInt128)num.ShiftRightTruncate(-e2);
        }
        else
        {
            den = den.ShiftLeft(-e2);
        }
        if (den.IsOne)
        {
            exact = true;
            return (UInt128)num;
        }
        BigInt.DivRem(num, den, out BigInt quo, out BigInt rem);
        exact = rem.IsZero;
        return (UInt128)quo;
    }

    /// <summary>The shortest decimal digits d with d * 10^exp10 rounding back to |value| (finite, nonzero).</summary>
    private static UInt128 ShortestDigits(UInt128 bits, out int exp10)
    {
        int ieeeExp = ExpOf(bits);
        UInt128 frac = FracOf(bits);
        UInt128 m2;
        int e2;
        if (ieeeExp == 0)
        {
            m2 = frac;
            e2 = 1 - ExponentBias - FractionBits - 2;
        }
        else
        {
            m2 = frac | Hidden;
            e2 = ieeeExp - ExponentBias - FractionBits - 2;
        }
        bool acceptBounds = (m2 & 1) == 0;
        UInt128 mv = m2 << 2;
        uint mmShift = frac != 0 || ieeeExp <= 1 ? 1u : 0u;
        UInt128 mp = mv + 2, mm = mv - 1 - mmShift;
        int bitLen = 128 - (int)UInt128.LeadingZeroCount(mp);
        int q = (int)Math.Floor((bitLen + e2) * Log10Of2) - 37;
        UInt128 vp = ScaleDown(mp, e2, q, out bool vpExact);
        UInt128 vr = ScaleDown(mv, e2, q, out bool vrTZ);
        UInt128 vm = ScaleDown(mm, e2, q, out bool vmTZ);
        if (vpExact && !acceptBounds)
        {
            vp--;
        }
        int removed = 0;
        int lastRemoved = 0;
        UInt128[] p10 = s_pow10u;
        for (int k = 16; k >= 1; k >>= 1)
        {
            UInt128 p = p10[k];
            while (true)
            {
                UInt128 vpk = vp / p, vmk = vm / p;
                if (vpk <= vmk)
                {
                    break;
                }
                UInt128 vrk = vr / p;
                UInt128 vrRem = vr - vrk * p;
                UInt128 low = p10[k - 1];
                vrTZ &= lastRemoved == 0 && vrRem % low == 0;
                lastRemoved = (int)(vrRem / low);
                vmTZ &= vm - vmk * p == 0;
                vp = vpk;
                vm = vmk;
                vr = vrk;
                removed += k;
            }
        }
        if (vmTZ && acceptBounds)
        {
            while (vm % 10 == 0)
            {
                vrTZ &= lastRemoved == 0;
                lastRemoved = (int)(vr % 10);
                vr /= 10;
                vp /= 10;
                vm /= 10;
                removed++;
            }
        }
        if (vrTZ && lastRemoved == 5 && (vr & 1) == 0)
        {
            lastRemoved = 4; // an exact tie rounds to even
        }
        bool up = (vr == vm && (!acceptBounds || !vmTZ)) || lastRemoved >= 5;
        exp10 = q + removed;
        return up ? vr + 1 : vr;
    }

    // ---------------------------------------------------------------
    // Exactly rounded digits
    // ---------------------------------------------------------------

    /// <summary>round(|value| * 10^k), ties to even (finite value).</summary>
    private static BigInt ScaledRound(UInt128 bits, long k)
    {
        int ieeeExp = ExpOf(bits);
        UInt128 frac = FracOf(bits);
        UInt128 m2 = ieeeExp == 0 ? frac : frac | Hidden;
        int e2 = (ieeeExp == 0 ? 1 : ieeeExp) - ExponentBias - FractionBits;
        BigInt num = new(m2);
        BigInt den = BigInt.One;
        if (k >= 0)
        {
            num *= Pow10((int)k);
        }
        else
        {
            den = Pow10((int)-k);
        }
        if (e2 >= 0)
        {
            num = num.ShiftLeft(e2);
        }
        else
        {
            den = den.ShiftLeft(-e2);
        }
        BigInt.DivRem(num, den, out BigInt q, out BigInt r);
        int c = r.IsZero ? -1 : (r << 1).CompareTo(den);
        if (c > 0 || (c == 0 && q.IsOdd))
        {
            q += BigInt.One;
        }
        return q;
    }

    /// <summary>|value| rounded to n significant digits (ties to even), with the decimal exponent of the first digit.</summary>
    private static string SignificantDigits(UInt128 bits, int n, out int exp10)
    {
        int e = (int)Math.Floor(ILogB(new Float128(bits & ~SignMask)) * Log10Of2);
        for (int attempt = 0; ; attempt++)
        {
            string s = ScaledRound(bits, (long)n - 1 - e).ToString();
            if (s.Length == n || attempt >= 3)
            {
                if (s.Length > n)
                {
                    e += s.Length - n; // 99..9 rounded up to 100..0
                    s = s[..n];
                }
                exp10 = e;
                return s;
            }
            e += s.Length - n;
        }
    }

    // ---------------------------------------------------------------
    // ToString / TryFormat
    // ---------------------------------------------------------------

    /// <summary>The shortest decimal string that parses back to exactly this value (like <see cref="double.ToString()"/>).</summary>
    public override string ToString() => Format(this, default, NumberFormatInfo.CurrentInfo);

    /// <summary>Formats with a standard format string: G[n], R, E[n], F[n] or N[n].</summary>
    public string ToString(string? format) => Format(this, format, NumberFormatInfo.CurrentInfo);

    /// <summary>Formats the shortest round-trip representation with culture-specific symbols.</summary>
    public string ToString(IFormatProvider? provider) => Format(this, default, NumberFormatInfo.GetInstance(provider));

    /// <summary>Formats with a standard format string (G[n], R, E[n], F[n], N[n]) and culture-specific symbols.</summary>
    public string ToString(string? format, IFormatProvider? provider) => Format(this, format, NumberFormatInfo.GetInstance(provider));

    /// <inheritdoc/>
    public bool TryFormat(Span<char> destination, out int charsWritten, ReadOnlySpan<char> format = default, IFormatProvider? provider = null)
    {
        string s = Format(this, format, NumberFormatInfo.GetInstance(provider));
        if (s.Length > destination.Length)
        {
            charsWritten = 0;
            return false;
        }
        s.CopyTo(destination);
        charsWritten = s.Length;
        return true;
    }

    /// <inheritdoc/>
    public bool TryFormat(Span<byte> utf8Destination, out int bytesWritten, ReadOnlySpan<char> format = default, IFormatProvider? provider = null) =>
        NumberText.TryWriteUtf8(Format(this, format, NumberFormatInfo.GetInstance(provider)), utf8Destination, out bytesWritten);

    private static string Format(Float128 value, ReadOnlySpan<char> format, NumberFormatInfo nfi)
    {
        char kind = 'G';
        int precision = -1;
        if (!format.IsEmpty)
        {
            kind = format[0];
            if (format.Length > 1)
            {
                if (!int.TryParse(format[1..], NumberStyles.None, CultureInfo.InvariantCulture, out precision) || precision > 999_999_999)
                {
                    throw new FormatException($"Format specifier '{format.ToString()}' is invalid.");
                }
            }
        }
        if (char.ToUpperInvariant(kind) is not ('G' or 'R' or 'E' or 'F' or 'N'))
        {
            throw new FormatException($"Format specifier '{format.ToString()}' is not supported by Float128.");
        }
        if (IsNaN(value))
        {
            return nfi.NaNSymbol;
        }
        bool neg = value.SignBit;
        if (IsInfinity(value))
        {
            return neg ? nfi.NegativeInfinitySymbol : nfi.PositiveInfinitySymbol;
        }
        StringBuilder sb = new(48);
        if (neg)
        {
            sb.Append(nfi.NegativeSign);
        }
        UInt128 bits = value._bits;
        bool zero = IsZero(value);
        switch (char.ToUpperInvariant(kind))
        {
            case 'R':
            case 'G':
            {
                string digits;
                int scale; // position of the decimal point relative to the first digit
                int maxDigits;
                if (zero)
                {
                    digits = "0";
                    scale = 1;
                    maxDigits = GeneralDigits;
                }
                else if (precision <= 0 || char.ToUpperInvariant(kind) == 'R')
                {
                    UInt128 d = ShortestDigits(bits, out int e10);
                    digits = d.ToString(CultureInfo.InvariantCulture);
                    scale = digits.Length + e10;
                    maxDigits = GeneralDigits;
                }
                else
                {
                    digits = SignificantDigits(bits, precision, out int e10);
                    scale = e10 + 1;
                    maxDigits = precision;
                }
                digits = digits.TrimEnd('0');
                if (digits.Length == 0)
                {
                    digits = "0";
                }
                if (scale > maxDigits || scale < -3)
                {
                    AppendScientific(sb, digits, scale - 1, kind == 'g' || kind == 'r' ? 'e' : 'E', 2, nfi);
                }
                else
                {
                    AppendFixed(sb, digits, scale, -1, nfi, group: false);
                }
                break;
            }
            case 'E':
            {
                int p = precision < 0 ? 6 : precision;
                string digits;
                int e10;
                if (zero)
                {
                    digits = new string('0', p + 1);
                    e10 = 0;
                }
                else
                {
                    digits = SignificantDigits(bits, p + 1, out e10);
                }
                sb.Append(digits[0]);
                if (p > 0)
                {
                    sb.Append(nfi.NumberDecimalSeparator).Append(digits, 1, p);
                }
                AppendExponent(sb, e10, kind == 'e' ? 'e' : 'E', 3, nfi);
                break;
            }
            default: // F, N
            {
                int p = precision < 0 ? nfi.NumberDecimalDigits : precision;
                string digits = zero ? "0" : ScaledRound(bits, p).ToString();
                // digits represent |value| * 10^p
                if (digits.Length <= p)
                {
                    digits = new string('0', p - digits.Length + 1) + digits;
                }
                AppendFixed(sb, digits, digits.Length - p, p, nfi, group: char.ToUpperInvariant(kind) == 'N');
                break;
            }
        }
        return sb.ToString();
    }

    /// <summary>Appends digits with the decimal point after <paramref name="scale"/> digits; decimals = -1 keeps all given digits.</summary>
    private static void AppendFixed(StringBuilder sb, string digits, int scale, int decimals, NumberFormatInfo nfi, bool group)
    {
        // integer part
        if (scale <= 0)
        {
            sb.Append('0');
        }
        else
        {
            int intLen = scale;
            int start = sb.Length;
            for (int i = 0; i < intLen; i++)
            {
                sb.Append(i < digits.Length ? digits[i] : '0');
            }
            if (group)
            {
                InsertGroups(sb, start, intLen, nfi);
            }
        }
        // fractional part
        int fracAvail = digits.Length - Math.Max(scale, 0);
        int leadingZeros = scale < 0 ? -scale : 0;
        int fracLen = decimals >= 0 ? decimals : leadingZeros + Math.Max(fracAvail, 0);
        if (fracLen > 0)
        {
            sb.Append(nfi.NumberDecimalSeparator);
            for (int i = 0; i < fracLen; i++)
            {
                int idx = i - leadingZeros + Math.Max(scale, 0);
                sb.Append(i < leadingZeros || idx >= digits.Length ? '0' : digits[idx]);
            }
        }
    }

    private static void InsertGroups(StringBuilder sb, int start, int length, NumberFormatInfo nfi)
    {
        int[] sizes = nfi.NumberGroupSizes;
        string sep = nfi.NumberGroupSeparator;
        if (sizes.Length == 0 || sep.Length == 0)
        {
            return;
        }
        // Walk from the right end of the integer part.
        int pos = start + length;
        int sizeIndex = 0;
        int size = sizes[0];
        while (size > 0 && pos - size > start)
        {
            pos -= size;
            sb.Insert(pos, sep);
            if (sizeIndex < sizes.Length - 1)
            {
                size = sizes[++sizeIndex];
            }
        }
    }

    private static void AppendScientific(StringBuilder sb, string digits, int exp, char e, int minExpDigits, NumberFormatInfo nfi)
    {
        sb.Append(digits[0]);
        if (digits.Length > 1)
        {
            sb.Append(nfi.NumberDecimalSeparator).Append(digits, 1, digits.Length - 1);
        }
        AppendExponent(sb, exp, e, minExpDigits, nfi);
    }

    private static void AppendExponent(StringBuilder sb, int exp, char e, int minDigits, NumberFormatInfo nfi)
    {
        sb.Append(e);
        sb.Append(exp < 0 ? nfi.NegativeSign : nfi.PositiveSign);
        string d = Math.Abs(exp).ToString(CultureInfo.InvariantCulture);
        if (d.Length < minDigits)
        {
            sb.Append('0', minDigits - d.Length);
        }
        sb.Append(d);
    }

    // ---------------------------------------------------------------
    // Parsing (correctly rounded from the exact decimal value)
    // ---------------------------------------------------------------

    private const NumberStyles DefaultStyle = NumberStyles.Float | NumberStyles.AllowThousands;

    public static Float128 Parse(string s) => Parse(s, DefaultStyle, null);

    public static Float128 Parse(string s, NumberStyles style) => Parse(s, style, null);

    public static Float128 Parse(string s, IFormatProvider? provider) => Parse(s, DefaultStyle, provider);

    public static Float128 Parse(string s, NumberStyles style, IFormatProvider? provider)
    {
        ArgumentNullException.ThrowIfNull(s);
        return Parse(s.AsSpan(), style, provider);
    }

    public static Float128 Parse(ReadOnlySpan<char> s, IFormatProvider? provider) => Parse(s, DefaultStyle, provider);

    public static Float128 Parse(ReadOnlySpan<char> s, NumberStyles style = DefaultStyle, IFormatProvider? provider = null) =>
        TryParse(s, style, provider, RoundingMode.NearestEven, out Float128 r)
            ? r
            : throw new FormatException("The input string was not in a correct format for Float128.");

    /// <summary>Parses (invariant culture, <see cref="NumberStyles.Float"/> | AllowThousands) with the given rounding.</summary>
    public static Float128 Parse(string s, RoundingMode rounding)
    {
        ArgumentNullException.ThrowIfNull(s);
        return TryParse(s, DefaultStyle, CultureInfo.InvariantCulture, rounding, out Float128 r)
            ? r
            : throw new FormatException("The input string was not in a correct format for Float128.");
    }

    public static bool TryParse(string? s, out Float128 result) => TryParse(s.AsSpan(), DefaultStyle, null, RoundingMode.NearestEven, out result);

    public static bool TryParse(ReadOnlySpan<char> s, out Float128 result) => TryParse(s, DefaultStyle, null, RoundingMode.NearestEven, out result);

    public static bool TryParse(string? s, IFormatProvider? provider, out Float128 result) =>
        TryParse(s.AsSpan(), DefaultStyle, provider, RoundingMode.NearestEven, out result);

    public static bool TryParse(ReadOnlySpan<char> s, IFormatProvider? provider, out Float128 result) =>
        TryParse(s, DefaultStyle, provider, RoundingMode.NearestEven, out result);

    public static bool TryParse(string? s, NumberStyles style, IFormatProvider? provider, out Float128 result) =>
        TryParse(s.AsSpan(), style, provider, RoundingMode.NearestEven, out result);

    public static bool TryParse(ReadOnlySpan<char> s, NumberStyles style, IFormatProvider? provider, out Float128 result) =>
        TryParse(s, style, provider, RoundingMode.NearestEven, out result);

    /// <summary>
    /// Parses a decimal number (culture symbols per <paramref name="provider"/>, NaN/Infinity symbols and the
    /// invariant nan/inf/infinity spellings), correctly rounded with <paramref name="rounding"/> from its exact value.
    /// Hex and binary specifiers are not supported.
    /// </summary>
    public static bool TryParse(ReadOnlySpan<char> s, NumberStyles style, IFormatProvider? provider, RoundingMode rounding, out Float128 result)
    {
        NumberText.ValidateStyle(style);
        if ((style & (NumberStyles.AllowHexSpecifier | NumberStyles.AllowBinarySpecifier)) != 0)
        {
            throw new ArgumentException("Hex and binary specifiers are not supported for Float128.", nameof(style));
        }
        if ((uint)rounding > (uint)RoundingMode.AwayFromZero)
        {
            throw new ArgumentOutOfRangeException(nameof(rounding));
        }
        result = default;
        NumberFormatInfo nfi = NumberFormatInfo.GetInstance(provider);
        ReadOnlySpan<char> t = NumberText.TrimForSymbols(s, style);
        if (TryParseSymbol(t, nfi, style, out result))
        {
            return true;
        }
        if (!NumberText.TryNormalize(s, style, provider, out string text, out _))
        {
            return false;
        }
        return TryParseNormalized(text, rounding, out result);
    }

    private static bool TryParseSymbol(ReadOnlySpan<char> t, NumberFormatInfo nfi, NumberStyles style, out Float128 result)
    {
        result = default;
        if (t.Equals(nfi.PositiveInfinitySymbol, StringComparison.OrdinalIgnoreCase))
        {
            result = PositiveInfinity;
            return true;
        }
        if (t.Equals(nfi.NegativeInfinitySymbol, StringComparison.OrdinalIgnoreCase))
        {
            result = NegativeInfinity;
            return true;
        }
        if (t.Equals(nfi.NaNSymbol, StringComparison.OrdinalIgnoreCase))
        {
            result = NaN;
            return true;
        }
        bool neg = false;
        ReadOnlySpan<char> w = t;
        if ((style & NumberStyles.AllowLeadingSign) != 0)
        {
            if (nfi.NegativeSign.Length > 0 && w.StartsWith(nfi.NegativeSign, StringComparison.Ordinal))
            {
                neg = true;
                w = w[nfi.NegativeSign.Length..];
            }
            else if (nfi.PositiveSign.Length > 0 && w.StartsWith(nfi.PositiveSign, StringComparison.Ordinal))
            {
                w = w[nfi.PositiveSign.Length..];
            }
            else if (w.Length > 0 && w[0] is '+' or '-')
            {
                neg = w[0] == '-';
                w = w[1..];
            }
        }
        if (w.Equals(nfi.NaNSymbol, StringComparison.OrdinalIgnoreCase) || w.Equals("nan", StringComparison.OrdinalIgnoreCase))
        {
            result = NaN;
            return true;
        }
        if (w.Equals(nfi.PositiveInfinitySymbol, StringComparison.OrdinalIgnoreCase) ||
            w.Equals("inf", StringComparison.OrdinalIgnoreCase) || w.Equals("infinity", StringComparison.OrdinalIgnoreCase) ||
            w.Equals("∞", StringComparison.Ordinal))
        {
            result = InfOf(neg);
            return true;
        }
        return false;
    }

    /// <summary>Parses the invariant form "[-]digits[.digits][e[-]digits]" produced by <see cref="NumberText.TryNormalize"/>.</summary>
    private static bool TryParseNormalized(ReadOnlySpan<char> s, RoundingMode rounding, out Float128 result)
    {
        result = default;
        bool neg = false;
        if (s.Length > 0 && s[0] == '-')
        {
            neg = true;
            s = s[1..];
        }
        int ePos = s.IndexOfAny('e', 'E');
        long exp10 = 0;
        if (ePos >= 0)
        {
            ReadOnlySpan<char> es = s[(ePos + 1)..];
            bool eneg = false;
            if (es.Length > 0 && es[0] is '+' or '-')
            {
                eneg = es[0] == '-';
                es = es[1..];
            }
            if (es.IsEmpty)
            {
                return false;
            }
            foreach (char c in es)
            {
                if (!char.IsAsciiDigit(c))
                {
                    return false;
                }
                if (exp10 < 1_000_000_000_000L)
                {
                    exp10 = exp10 * 10 + (c - '0');
                }
            }
            if (eneg)
            {
                exp10 = -exp10;
            }
            s = s[..ePos];
        }
        // Collect the significant digits (no leading zeros) and count fractional digits.
        char[]? rented = null;
        Span<char> buf = s.Length <= 256 ? stackalloc char[256] : (rented = System.Buffers.ArrayPool<char>.Shared.Rent(s.Length));
        try
        {
            int n = 0;
            bool seenDot = false, anyDigit = false;
            long fracDigits = 0;
            foreach (char c in s)
            {
                if (c == '.')
                {
                    if (seenDot)
                    {
                        return false;
                    }
                    seenDot = true;
                    continue;
                }
                if (!char.IsAsciiDigit(c))
                {
                    return false;
                }
                anyDigit = true;
                if (seenDot)
                {
                    fracDigits++;
                }
                if (n == 0 && c == '0')
                {
                    continue;
                }
                buf[n++] = c;
            }
            if (!anyDigit)
            {
                return false;
            }
            exp10 -= fracDigits;
            // Drop trailing zeros of the significand.
            while (n > 0 && buf[n - 1] == '0')
            {
                n--;
                exp10++;
            }
            if (n == 0)
            {
                result = ZeroOf(neg);
                return true;
            }
            long magnitude = n + exp10; // value in [10^(magnitude-1), 10^magnitude)
            if (magnitude > 4934)
            {
                result = Huge(neg, rounding);
                return true;
            }
            if (magnitude < -4967)
            {
                result = Tiny(neg, rounding);
                return true;
            }
            BigInt mant = ParseDigits(buf[..n]);
            if (exp10 >= 0)
            {
                result = RoundBig(neg, mant * BigInt.Pow(5, (ulong)exp10), exp10, sticky: false, rounding);
            }
            else
            {
                result = RoundRational(neg, mant, BigInt.Pow(5, (ulong)-exp10), exp10, rounding);
            }
            return true;
        }
        finally
        {
            if (rented is not null)
            {
                System.Buffers.ArrayPool<char>.Shared.Return(rented);
            }
        }
    }

    private static BigInt ParseDigits(ReadOnlySpan<char> digits)
    {
        if (digits.Length <= 38)
        {
            UInt128 v = UInt128.Zero;
            foreach (char c in digits)
            {
                v = v * 10 + (uint)(c - '0');
            }
            return new BigInt(v);
        }
        BigInt.TryParse(digits, 10, out BigInt m);
        return m;
    }

    // UTF-8

    public static Float128 Parse(ReadOnlySpan<byte> utf8Text, IFormatProvider? provider) => Parse(utf8Text, DefaultStyle, provider);

    public static bool TryParse(ReadOnlySpan<byte> utf8Text, IFormatProvider? provider, out Float128 result) =>
        TryParse(utf8Text, DefaultStyle, provider, out result);

    public static Float128 Parse(ReadOnlySpan<byte> utf8Text, NumberStyles style = DefaultStyle, IFormatProvider? provider = null) =>
        TryParse(utf8Text, style, provider, out Float128 r)
            ? r
            : throw new FormatException("The input string was not in a correct format for Float128.");

    public static bool TryParse(ReadOnlySpan<byte> utf8Text, NumberStyles style, IFormatProvider? provider, out Float128 result)
    {
        NumberText.ValidateStyle(style);
        return NumberText.TryParseUtf8(utf8Text, (ReadOnlySpan<char> s, out Float128 r) => TryParse(s, style, provider, RoundingMode.NearestEven, out r), out result);
    }
}
