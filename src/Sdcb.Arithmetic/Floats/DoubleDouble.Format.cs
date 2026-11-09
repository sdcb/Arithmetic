using System.Globalization;
using Sdcb.Arithmetic.Core;

namespace Sdcb.Arithmetic;

public readonly partial struct DoubleDouble
{
    /// <summary>Significant digits tried first by the round-trip formatter (enough for most values).</summary>
    private const int DefaultDigits = 32;

    // ---------------------------------------------------------------
    // Formatting (invariant culture, through the exact Real value)
    // ---------------------------------------------------------------

    /// <summary>
    /// The shortest decimal string, starting from 32 significant digits, that parses back to this
    /// exact (Hi, Lo) pair. Values whose words are far apart (e.g. 1 + 2^-200) need more digits.
    /// </summary>
    public override string ToString()
    {
        if (!double.IsFinite(_hi) || _hi == 0.0)
        {
            return ToString(DefaultDigits);
        }
        Real exact = ToRealExact();
        for (int nd = DefaultDigits; ; nd = nd < 40 ? nd + 1 : nd * 2)
        {
            string s = exact.ToString(nd);
            if (nd >= 1200 || (TryParseInvariant(s, out DoubleDouble back) && back.Equals(this)))
            {
                return s;
            }
        }
    }

    /// <summary>Formats with the given number of significant decimal digits (correctly rounded from the exact value).</summary>
    public string ToString(int significantDigits)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(significantDigits, 1);
        return ToRealExact().ToString(significantDigits);
    }

    /// <summary>
    /// Formats in the invariant culture. Supported formats: null, "", "G" or "R" (round-trip, see
    /// <see cref="ToString()"/>) and "G<i>n</i>" (<i>n</i> significant digits).
    /// </summary>
    public string ToString(string? format, IFormatProvider? formatProvider) => ToString(format.AsSpan());

    private string ToString(ReadOnlySpan<char> format)
    {
        if (format.IsEmpty || format.Equals("G", StringComparison.OrdinalIgnoreCase) || format.Equals("R", StringComparison.OrdinalIgnoreCase))
        {
            return ToString();
        }
        if ((format[0] == 'G' || format[0] == 'g') &&
            int.TryParse(format.Slice(1), NumberStyles.None, CultureInfo.InvariantCulture, out int digits) && digits > 0)
        {
            return ToString(digits);
        }
        throw new FormatException($"Unsupported format '{format.ToString()}' for DoubleDouble.");
    }

    public bool TryFormat(Span<char> destination, out int charsWritten, ReadOnlySpan<char> format, IFormatProvider? provider)
    {
        string s = ToString(format);
        if (s.Length > destination.Length)
        {
            charsWritten = 0;
            return false;
        }
        s.CopyTo(destination);
        charsWritten = s.Length;
        return true;
    }

    public bool TryFormat(Span<byte> utf8Destination, out int bytesWritten, ReadOnlySpan<char> format, IFormatProvider? provider) =>
        NumberText.TryWriteUtf8(ToString(format), utf8Destination, out bytesWritten);

    // ---------------------------------------------------------------
    // Parsing: the decimal is rounded to a Real with enough bits for every digit, then to the
    // nearest double-double (Hi = RN(v), Lo = RN(v - Hi)).
    // ---------------------------------------------------------------

    private static bool TryParseInvariant(ReadOnlySpan<char> s, out DoubleDouble result)
    {
        int prec = (int)Math.Min(4L * s.Length + 192, 1 << 20);
        if (Real.TryParse(s, prec, RealRounding.NearestEven, out Real r))
        {
            result = FromReal(r);
            return true;
        }
        result = default;
        return false;
    }

    /// <summary>Parses the invariant form "[±]digits[.digits][e±digits]", "NaN" or "[±]Infinity" (the provider is ignored, like <see cref="Real"/>).</summary>
    public static DoubleDouble Parse(string s, IFormatProvider? provider = null)
    {
        ArgumentNullException.ThrowIfNull(s);
        return Parse(s.AsSpan(), provider);
    }

    public static DoubleDouble Parse(ReadOnlySpan<char> s, IFormatProvider? provider = null) =>
        TryParseInvariant(s, out DoubleDouble r) ? r : throw new FormatException("The value could not be parsed as a DoubleDouble.");

    public static bool TryParse(string? s, IFormatProvider? provider, out DoubleDouble result) =>
        TryParseInvariant(s.AsSpan(), out result);

    public static bool TryParse(ReadOnlySpan<char> s, IFormatProvider? provider, out DoubleDouble result) =>
        TryParseInvariant(s, out result);

    /// <summary>Parses with .NET <see cref="NumberStyles"/> and the provider's symbols (hex/binary specifiers are not supported).</summary>
    public static DoubleDouble Parse(string s, NumberStyles style, IFormatProvider? provider)
    {
        ArgumentNullException.ThrowIfNull(s);
        return Parse(s.AsSpan(), style, provider);
    }

    public static DoubleDouble Parse(ReadOnlySpan<char> s, NumberStyles style, IFormatProvider? provider) =>
        TryParse(s, style, provider, out DoubleDouble r) ? r : throw new FormatException("The value could not be parsed as a DoubleDouble.");

    public static bool TryParse(string? s, NumberStyles style, IFormatProvider? provider, out DoubleDouble result) =>
        TryParse(s.AsSpan(), style, provider, out result);

    public static bool TryParse(ReadOnlySpan<char> s, NumberStyles style, IFormatProvider? provider, out DoubleDouble result)
    {
        NumberText.ValidateStyle(style);
        if ((style & (NumberStyles.AllowHexSpecifier | NumberStyles.AllowBinarySpecifier)) != 0)
        {
            throw new ArgumentException("Hex and binary specifiers are not supported for DoubleDouble.", nameof(style));
        }
        result = default;
        ReadOnlySpan<char> t = NumberText.TrimForSymbols(s, style);
        NumberFormatInfo nfi = NumberFormatInfo.GetInstance(provider);
        if (t.Equals(nfi.NaNSymbol, StringComparison.OrdinalIgnoreCase))
        {
            result = NaN;
            return true;
        }
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
        if (NumberText.TryNormalize(s, style, provider, out string text, out _))
        {
            return TryParseInvariant(text, out result);
        }
        ReadOnlySpan<char> word = t.Length > 0 && t[0] is '+' or '-' ? t.Slice(1) : t;
        if (word.Equals("nan", StringComparison.OrdinalIgnoreCase) || word.Equals("inf", StringComparison.OrdinalIgnoreCase) ||
            word.Equals("infinity", StringComparison.OrdinalIgnoreCase))
        {
            return TryParseInvariant(t, out result);
        }
        return false;
    }

    public static DoubleDouble Parse(ReadOnlySpan<byte> utf8Text, IFormatProvider? provider) =>
        TryParse(utf8Text, provider, out DoubleDouble r) ? r : throw new FormatException("The value could not be parsed as a DoubleDouble.");

    public static bool TryParse(ReadOnlySpan<byte> utf8Text, IFormatProvider? provider, out DoubleDouble result) =>
        NumberText.TryParseUtf8(utf8Text, (ReadOnlySpan<char> s, out DoubleDouble r) => TryParseInvariant(s, out r), out result);

    public static DoubleDouble Parse(ReadOnlySpan<byte> utf8Text, NumberStyles style, IFormatProvider? provider) =>
        TryParse(utf8Text, style, provider, out DoubleDouble r) ? r : throw new FormatException("The value could not be parsed as a DoubleDouble.");

    public static bool TryParse(ReadOnlySpan<byte> utf8Text, NumberStyles style, IFormatProvider? provider, out DoubleDouble result) =>
        NumberText.TryParseUtf8(utf8Text, (ReadOnlySpan<char> s, out DoubleDouble r) => TryParse(s, style, provider, out r), out result);
}
