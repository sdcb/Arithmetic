using System.Buffers;
using System.Globalization;
using System.Text;

namespace Sdcb.Arithmetic.Core;

/// <summary>
/// Shared NumberStyles / UTF-8 plumbing for the generic-math parse and format entry points.
/// Culture-specific input is normalized to the invariant "[-]digits[.digits][e[-]digits]"
/// form that the type-specific parsers understand.
/// </summary>
internal static class NumberText
{
    internal enum Radix
    {
        Decimal,
        Hex,
        Binary,
    }

    private const NumberStyles WhiteStyles = NumberStyles.AllowLeadingWhite | NumberStyles.AllowTrailingWhite;

    public static void ValidateStyle(NumberStyles style)
    {
        const NumberStyles all = NumberStyles.Any | NumberStyles.AllowHexSpecifier | NumberStyles.AllowBinarySpecifier;
        if ((style & ~all) != 0)
        {
            throw new ArgumentException("Invalid NumberStyles value.", nameof(style));
        }
        bool hex = (style & NumberStyles.AllowHexSpecifier) != 0;
        bool bin = (style & NumberStyles.AllowBinarySpecifier) != 0;
        if ((hex && bin) || ((hex || bin) && (style & ~(WhiteStyles | NumberStyles.AllowHexSpecifier | NumberStyles.AllowBinarySpecifier)) != 0))
        {
            throw new ArgumentException("Hex and binary specifiers can only be combined with white-space styles.", nameof(style));
        }
    }

    private static bool IsWhite(char c) => c == ' ' || (c >= '\t' && c <= '\r');

    private static ReadOnlySpan<char> TrimWhite(ReadOnlySpan<char> s, NumberStyles style)
    {
        if ((style & NumberStyles.AllowLeadingWhite) != 0)
        {
            int i = 0;
            while (i < s.Length && IsWhite(s[i]))
            {
                i++;
            }
            s = s.Slice(i);
        }
        if ((style & NumberStyles.AllowTrailingWhite) != 0)
        {
            int j = s.Length;
            while (j > 0 && IsWhite(s[j - 1]))
            {
                j--;
            }
            s = s.Slice(0, j);
        }
        return s;
    }

    /// <summary>
    /// Validates <paramref name="s"/> against <paramref name="style"/> and produces either the
    /// invariant decimal form "[-]digits[.digits][e[-]digits]" or the bare hex/binary digit string.
    /// </summary>
    public static bool TryNormalize(ReadOnlySpan<char> s, NumberStyles style, IFormatProvider? provider, out string text, out Radix radix)
    {
        ValidateStyle(style);
        text = string.Empty;
        s = TrimWhite(s, style);
        if ((style & NumberStyles.AllowHexSpecifier) != 0)
        {
            radix = Radix.Hex;
            if (s.IsEmpty || s.ContainsAnyExcept("0123456789abcdefABCDEF"))
            {
                return false;
            }
            text = s.ToString();
            return true;
        }
        if ((style & NumberStyles.AllowBinarySpecifier) != 0)
        {
            radix = Radix.Binary;
            if (s.IsEmpty || s.ContainsAnyExcept("01"))
            {
                return false;
            }
            text = s.ToString();
            return true;
        }

        radix = Radix.Decimal;
        NumberFormatInfo nfi = NumberFormatInfo.GetInstance(provider);
        bool neg = false, signSeen = false;

        if ((style & NumberStyles.AllowParentheses) != 0 && s.Length >= 2 && s[0] == '(' && s[^1] == ')')
        {
            neg = true;
            signSeen = true;
            s = s.Slice(1, s.Length - 2);
        }
        if (!signSeen && (style & NumberStyles.AllowLeadingSign) != 0 && TryStripSign(ref s, nfi, leading: true, out bool n1))
        {
            neg = n1;
            signSeen = true;
        }
        if ((style & NumberStyles.AllowCurrencySymbol) != 0 && nfi.CurrencySymbol.Length > 0)
        {
            if (s.StartsWith(nfi.CurrencySymbol, StringComparison.Ordinal))
            {
                s = TrimWhite(s.Slice(nfi.CurrencySymbol.Length), style & NumberStyles.AllowLeadingWhite);
            }
            else if (s.EndsWith(nfi.CurrencySymbol, StringComparison.Ordinal))
            {
                s = TrimWhite(s.Slice(0, s.Length - nfi.CurrencySymbol.Length), style & NumberStyles.AllowTrailingWhite);
            }
        }
        if (!signSeen && (style & NumberStyles.AllowLeadingSign) != 0 && TryStripSign(ref s, nfi, leading: true, out bool n2))
        {
            neg = n2;
            signSeen = true;
        }
        if (!signSeen && (style & NumberStyles.AllowTrailingSign) != 0 && TryStripSign(ref s, nfi, leading: false, out bool n3))
        {
            neg = n3;
        }

        bool allowThousands = (style & NumberStyles.AllowThousands) != 0;
        bool allowDot = (style & NumberStyles.AllowDecimalPoint) != 0;
        bool allowExp = (style & NumberStyles.AllowExponent) != 0;
        string group = nfi.NumberGroupSeparator;
        string dec = nfi.NumberDecimalSeparator;

        var sb = new StringBuilder(s.Length + 1);
        if (neg)
        {
            sb.Append('-');
        }
        int digits = 0;
        bool seenDot = false;
        int i = 0;
        while (i < s.Length)
        {
            char c = s[i];
            if (char.IsAsciiDigit(c))
            {
                sb.Append(c);
                digits++;
                i++;
            }
            else if (allowDot && !seenDot && dec.Length > 0 && s.Slice(i).StartsWith(dec, StringComparison.Ordinal))
            {
                sb.Append('.');
                seenDot = true;
                i += dec.Length;
            }
            else if (allowThousands && !seenDot && digits > 0 && group.Length > 0 && s.Slice(i).StartsWith(group, StringComparison.Ordinal))
            {
                i += group.Length;
            }
            else if (allowExp && digits > 0 && (c == 'e' || c == 'E'))
            {
                i++;
                ReadOnlySpan<char> rest = s.Slice(i);
                bool expNeg = false;
                if (TryStripSign(ref rest, nfi, leading: true, out bool en))
                {
                    expNeg = en;
                }
                else if (!rest.IsEmpty && rest[0] is '+' or '-')
                {
                    expNeg = rest[0] == '-';
                    rest = rest.Slice(1);
                }
                if (rest.IsEmpty || rest.ContainsAnyExceptInRange('0', '9'))
                {
                    return false;
                }
                sb.Append('e');
                if (expNeg)
                {
                    sb.Append('-');
                }
                sb.Append(rest);
                i = s.Length;
            }
            else
            {
                return false;
            }
        }
        if (digits == 0)
        {
            return false;
        }
        text = sb.ToString();
        return true;
    }

    private static bool TryStripSign(ref ReadOnlySpan<char> s, NumberFormatInfo nfi, bool leading, out bool negative)
    {
        negative = false;
        string ns = nfi.NegativeSign, ps = nfi.PositiveSign;
        if (leading)
        {
            if (ns.Length > 0 && s.StartsWith(ns, StringComparison.Ordinal))
            {
                negative = true;
                s = s.Slice(ns.Length);
                return true;
            }
            if (ps.Length > 0 && s.StartsWith(ps, StringComparison.Ordinal))
            {
                s = s.Slice(ps.Length);
                return true;
            }
        }
        else
        {
            if (ns.Length > 0 && s.EndsWith(ns, StringComparison.Ordinal))
            {
                negative = true;
                s = s.Slice(0, s.Length - ns.Length);
                return true;
            }
            if (ps.Length > 0 && s.EndsWith(ps, StringComparison.Ordinal))
            {
                s = s.Slice(0, s.Length - ps.Length);
                return true;
            }
        }
        return false;
    }

    /// <summary>Strips leading/trailing white space allowed by the style (used for NaN/Infinity symbol checks).</summary>
    public static ReadOnlySpan<char> TrimForSymbols(ReadOnlySpan<char> s, NumberStyles style) => TrimWhite(s, style);

    // ---------------------------------------------------------------
    // UTF-8 helpers
    // ---------------------------------------------------------------

    public delegate bool CharParser<T>(ReadOnlySpan<char> s, out T result);

    public static bool TryParseUtf8<T>(ReadOnlySpan<byte> utf8, CharParser<T> parser, out T result)
    {
        int max = Encoding.UTF8.GetMaxCharCount(utf8.Length);
        char[]? rented = null;
        Span<char> buf = max <= 256 ? stackalloc char[256] : (rented = ArrayPool<char>.Shared.Rent(max));
        try
        {
            int n;
            try
            {
                n = Encoding.UTF8.GetChars(utf8, buf);
            }
            catch (DecoderFallbackException)
            {
                result = default!;
                return false;
            }
            return parser(buf.Slice(0, n), out result);
        }
        finally
        {
            if (rented is not null)
            {
                ArrayPool<char>.Shared.Return(rented);
            }
        }
    }

    public static bool TryWriteUtf8(string s, Span<byte> destination, out int bytesWritten)
    {
        if (Encoding.UTF8.GetByteCount(s) > destination.Length)
        {
            bytesWritten = 0;
            return false;
        }
        bytesWritten = Encoding.UTF8.GetBytes(s, destination);
        return true;
    }
}
