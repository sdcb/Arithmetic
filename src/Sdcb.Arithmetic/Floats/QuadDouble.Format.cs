using System.Globalization;
using Sdcb.Arithmetic.Core;

namespace Sdcb.Arithmetic;

public readonly partial struct QuadDouble
{
    // ---------------------------------------------------------------
    // Formatting (exact value printed through BigFloat; invariant culture)
    // ---------------------------------------------------------------

    /// <summary>
    /// Formats with enough significant digits that <see cref="Parse(string, IFormatProvider?)"/> returns the
    /// same value: at least 67, and three more than the decimal width of the bits spanned by the components
    /// when they have gaps between them. Trailing zeros are trimmed, so exact short decimals stay short.
    /// </summary>
    private const int RoundTripDigits = 67;
    public override string ToString() => ToString(0);

    /// <summary>Formats with the given number of significant decimal digits (0 = enough to round-trip).</summary>
    public string ToString(int significantDigits)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(significantDigits);
        double x0 = _x0;
        if (double.IsNaN(x0))
        {
            return "NaN";
        }
        if (double.IsInfinity(x0))
        {
            return x0 > 0 ? "Infinity" : "-Infinity";
        }
        if (x0 == 0)
        {
            return double.IsNegative(x0) ? "-0" : "0";
        }
        int nd = significantDigits;
        if (nd == 0)
        {
            double last = _x3 != 0 ? _x3 : _x2 != 0 ? _x2 : _x1 != 0 ? _x1 : x0;
            int span = Math.ILogB(x0) - Math.Max(Math.ILogB(last) - 52, -1074) + 1;
            nd = Math.Max(RoundTripDigits, 3 + (int)Math.Ceiling(span * 0.3010299956639812));
        }
        string s = ToBigFloatExact(this).ToString(nd);
        return significantDigits == 0 && nd > RoundTripDigits ? PadDigits(s, nd) : s;
    }

    /// <summary>
    /// Restores trailing zeros trimmed by the formatter, so the digit count tells the parser how
    /// precise the text is (needed when the components span more than 212 bits).
    /// </summary>
    private static string PadDigits(string s, int nd)
    {
        int ePos = s.IndexOf('e');
        string mant = ePos < 0 ? s : s[..ePos];
        int count = 0;
        bool leading = true;
        foreach (char c in mant)
        {
            if (char.IsAsciiDigit(c))
            {
                leading &= c == '0';
                if (!leading)
                {
                    count++;
                }
            }
        }
        if (count >= nd)
        {
            return s;
        }
        string padded = (mant.Contains('.') ? mant : mant + ".") + new string('0', nd - count);
        return ePos < 0 ? padded : padded + s[ePos..];
    }

    /// <summary>
    /// Formats the value. <paramref name="format"/> may be null, empty, "G" or "R" (round-trip digits) or
    /// "G<i>n</i>" / "E<i>n</i>" (<i>n</i> significant digits); output is culture-invariant like <see cref="BigFloat"/>.
    /// </summary>
    public string ToString(string? format, IFormatProvider? formatProvider = null) => ToString(DigitsFromFormat(format));

    private static int DigitsFromFormat(ReadOnlySpan<char> format)
    {
        if (format.IsEmpty)
        {
            return 0;
        }
        char f = char.ToUpperInvariant(format[0]);
        if (f is not ('G' or 'R' or 'E'))
        {
            throw new FormatException($"Format '{format.ToString()}' is not supported by QuadDouble.");
        }
        if (format.Length == 1)
        {
            return 0;
        }
        if (f == 'R' || !int.TryParse(format.Slice(1), NumberStyles.None, CultureInfo.InvariantCulture, out int n) || n > 100000)
        {
            throw new FormatException($"Format '{format.ToString()}' is not supported by QuadDouble.");
        }
        return n;
    }

    public bool TryFormat(Span<char> destination, out int charsWritten, ReadOnlySpan<char> format = default, IFormatProvider? provider = null)
    {
        string s = ToString(DigitsFromFormat(format));
        if (s.Length > destination.Length)
        {
            charsWritten = 0;
            return false;
        }
        s.CopyTo(destination);
        charsWritten = s.Length;
        return true;
    }

    public bool TryFormat(Span<byte> utf8Destination, out int bytesWritten, ReadOnlySpan<char> format = default, IFormatProvider? provider = null) =>
        NumberText.TryWriteUtf8(ToString(DigitsFromFormat(format)), utf8Destination, out bytesWritten);

    // ---------------------------------------------------------------
    // Parsing (decimal text rounded through a BigFloat sized from the digit count)
    // ---------------------------------------------------------------

    /// <summary>
    /// Parses the invariant form <c>[+|-]digits[.digits][e[+|-]digits]</c>, <c>NaN</c>, <c>Infinity</c>/<c>inf</c>.
    /// The decimal is rounded once to max(214, ⌊(digits - 1)·log2 10⌋ - 2) bits, then split into QD components.
    /// </summary>
    private static bool TryParseInvariant(ReadOnlySpan<char> s, out QuadDouble result)
    {
        s = s.Trim();
        result = default;
        if (s.IsEmpty)
        {
            return false;
        }
        // Round once to the precision the digits carry (at least 214 bits, which every QD component
        // pattern of a 212-bit value fits): the decimal error of the text is then not mistaken for
        // extra low-order components, so ToString output round-trips exactly.
        int digits = 0;
        bool leading = true;
        foreach (char ch in s)
        {
            if (ch is 'e' or 'E')
            {
                break;
            }
            if (char.IsAsciiDigit(ch))
            {
                if (ch != '0')
                {
                    leading = false;
                }
                if (!leading)
                {
                    digits++;
                }
            }
        }
        long prec = Math.Max(214, (long)Math.Floor((digits - 1) * 3.321928094887362) - 2);
        if (!BigFloat.TryParse(s, (int)Math.Min(prec, 1 << 24), RoundingMode.NearestEven, out BigFloat r))
        {
            return false;
        }
        result = FromBigFloat(r);
        return true;
    }

    public static QuadDouble Parse(string s, IFormatProvider? provider = null)
    {
        ArgumentNullException.ThrowIfNull(s);
        return Parse(s.AsSpan(), provider);
    }

    public static QuadDouble Parse(ReadOnlySpan<char> s, IFormatProvider? provider = null) =>
        TryParseInvariant(s, out QuadDouble r) ? r : throw new FormatException("The value could not be parsed as a QuadDouble.");

    public static bool TryParse(string? s, IFormatProvider? provider, out QuadDouble result) =>
        TryParseInvariant(s.AsSpan(), out result);

    public static bool TryParse(ReadOnlySpan<char> s, IFormatProvider? provider, out QuadDouble result) =>
        TryParseInvariant(s, out result);

    public static bool TryParse(string? s, out QuadDouble result) => TryParseInvariant(s.AsSpan(), out result);

    public static bool TryParse(ReadOnlySpan<char> s, out QuadDouble result) => TryParseInvariant(s, out result);

    /// <summary>Parses with .NET <see cref="NumberStyles"/> and the provider's symbols (hex/binary specifiers are not supported).</summary>
    public static QuadDouble Parse(ReadOnlySpan<char> s, NumberStyles style, IFormatProvider? provider) =>
        TryParse(s, style, provider, out QuadDouble r) ? r : throw new FormatException("The value could not be parsed as a QuadDouble.");

    public static QuadDouble Parse(string s, NumberStyles style, IFormatProvider? provider)
    {
        ArgumentNullException.ThrowIfNull(s);
        return Parse(s.AsSpan(), style, provider);
    }

    public static bool TryParse(string? s, NumberStyles style, IFormatProvider? provider, out QuadDouble result) =>
        TryParse(s.AsSpan(), style, provider, out result);

    public static bool TryParse(ReadOnlySpan<char> s, NumberStyles style, IFormatProvider? provider, out QuadDouble result)
    {
        NumberText.ValidateStyle(style);
        if ((style & (NumberStyles.AllowHexSpecifier | NumberStyles.AllowBinarySpecifier)) != 0)
        {
            throw new ArgumentException("Hex and binary specifiers are not supported for QuadDouble.", nameof(style));
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

    public static QuadDouble Parse(ReadOnlySpan<byte> utf8Text, IFormatProvider? provider = null) =>
        TryParse(utf8Text, provider, out QuadDouble r) ? r : throw new FormatException("The value could not be parsed as a QuadDouble.");

    public static bool TryParse(ReadOnlySpan<byte> utf8Text, IFormatProvider? provider, out QuadDouble result) =>
        NumberText.TryParseUtf8(utf8Text, (ReadOnlySpan<char> s, out QuadDouble r) => TryParseInvariant(s, out r), out result);

    public static QuadDouble Parse(ReadOnlySpan<byte> utf8Text, NumberStyles style, IFormatProvider? provider) =>
        TryParse(utf8Text, style, provider, out QuadDouble r) ? r : throw new FormatException("The value could not be parsed as a QuadDouble.");

    public static bool TryParse(ReadOnlySpan<byte> utf8Text, NumberStyles style, IFormatProvider? provider, out QuadDouble result) =>
        NumberText.TryParseUtf8(utf8Text, (ReadOnlySpan<char> s, out QuadDouble r) => TryParse(s, style, provider, out r), out result);
}
