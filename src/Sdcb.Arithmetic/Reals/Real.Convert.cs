using System.Diagnostics;
using System.Globalization;

namespace Sdcb.Arithmetic;

public readonly partial struct Real
{
    // ---------------------------------------------------------------
    // Formatting (correctly rounded decimal, mpfr_get_str style)
    // ---------------------------------------------------------------

    /// <summary>
    /// Formats with just enough decimal digits to read the value back exactly
    /// (1 + ceil(p*log10(2)) significant digits).
    /// </summary>
    public override string ToString() => ToString(0, RealRounding.NearestEven);

    /// <summary>Formats with the given number of significant decimal digits (0 = enough for round-trip).</summary>
    public string ToString(int significantDigits, RealRounding rounding = RealRounding.NearestEven)
    {
        switch (_kind)
        {
            case Kind.NaN:
                return "NaN";
            case Kind.Inf:
                return IsNeg ? "-Infinity" : "Infinity";
            case Kind.Zero:
                return IsNeg ? "-0" : "0";
        }
        int nd = significantDigits > 0
            ? significantDigits
            : 1 + (int)Math.Ceiling(Precision * 0.3010299956639812);

        // Estimated decimal exponent E with v in [10^E, 10^(E+1)); corrected below.
        long msb = BinaryExponent - 1;
        long e10 = (long)Math.Floor(msb * 0.3010299956639812);

        string digits = DecimalDigits(nd, e10, rounding, out long e10Fixed);
        e10 = e10Fixed;

        string trimmed = digits.TrimEnd('0');
        if (trimmed.Length == 0)
        {
            trimmed = "0";
        }
        string sign = IsNeg ? "-" : "";
        if (e10 >= 0 && e10 <= 20)
        {
            if (e10 + 1 >= trimmed.Length)
            {
                return sign + trimmed + new string('0', (int)(e10 + 1 - trimmed.Length));
            }
            return sign + trimmed[..(int)(e10 + 1)] + "." + trimmed[(int)(e10 + 1)..];
        }
        if (e10 < 0 && e10 >= -5)
        {
            return sign + "0." + new string('0', (int)(-e10 - 1)) + trimmed;
        }
        string fracPart = trimmed.Length > 1 ? "." + trimmed[1..] : "";
        return $"{sign}{trimmed[0]}{fracPart}e{(e10 >= 0 ? "+" : "")}{e10}";
    }

    /// <summary>Correctly-rounded first nd decimal digits; fixes up the decimal exponent estimate.</summary>
    private string DecimalDigits(int nd, long e10, RealRounding rounding, out long e10Fixed)
    {
        for (int attempt = 0; ; attempt++)
        {
            BigInt i = ScaledDecimal(nd - 1 - e10, rounding);
            string s = i.ToString();
            if (s.Length == nd || attempt == 2)
            {
                e10Fixed = e10;
                if (s.Length > nd)
                {
                    // Rounding produced one extra digit (e.g. 999.. -> 1000..).
                    e10Fixed += s.Length - nd;
                    s = s[..nd];
                }
                return s;
            }
            e10 += s.Length - nd;
        }
    }

    /// <summary>Rounds |v| * 10^k to an integer under the value's sign and mode.</summary>
    private BigInt ScaledDecimal(long k, RealRounding rnd)
    {
        BigInt num = _man;
        BigInt den = BigInt.One;
        if (k >= 0)
        {
            num *= BigInt.Pow(BigInt.Ten, (ulong)k);
        }
        else
        {
            den = BigInt.Pow(BigInt.Ten, (ulong)(-k));
        }
        if (_exp >= 0)
        {
            num = num.ShiftLeft(_exp);
        }
        else
        {
            den = den.ShiftLeft(-_exp);
        }
        BigInt.DivRem(num, den, out BigInt q, out BigInt r);
        if (r.IsZero)
        {
            return q;
        }
        bool up = rnd switch
        {
            RealRounding.NearestEven => Compare2R(r, den, q),
            RealRounding.TowardZero => false,
            RealRounding.AwayFromZero => true,
            RealRounding.TowardPositive => !IsNeg,
            RealRounding.TowardNegative => IsNeg,
            _ => false,
        };
        return up ? q + BigInt.One : q;

        static bool Compare2R(BigInt r, BigInt den, BigInt q)
        {
            int c = (r << 1).CompareTo(den);
            return c > 0 || (c == 0 && q.IsOdd);
        }
    }

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

    // ---------------------------------------------------------------
    // Parsing (correctly rounded, mpfr_set_str for base 10)
    // ---------------------------------------------------------------

    public static Real Parse(string s, IFormatProvider? provider = null) => Parse(s.AsSpan(), provider);

    public static Real Parse(ReadOnlySpan<char> s, IFormatProvider? provider = null)
        => Parse(s, 0, RealRounding.NearestEven);

    public static Real Parse(ReadOnlySpan<char> s, int precision, RealRounding rounding = RealRounding.NearestEven)
    {
        if (!TryParse(s, precision, rounding, out Real result))
        {
            throw new FormatException("The value could not be parsed as a Real.");
        }
        return result;
    }

    public static bool TryParse(string? s, IFormatProvider? provider, out Real result)
        => TryParse(s.AsSpan(), 0, RealRounding.NearestEven, out result);

    public static bool TryParse(ReadOnlySpan<char> s, IFormatProvider? provider, out Real result)
        => TryParse(s, 0, RealRounding.NearestEven, out result);

    public static bool TryParse(ReadOnlySpan<char> s, int precision, RealRounding rounding, out Real result)
    {
        result = default;
        s = s.Trim();
        if (s.IsEmpty)
        {
            return false;
        }
        int prec = ResolvePrecision(precision, DefaultPrecision);

        bool neg = false;
        if (s[0] is '+' or '-')
        {
            neg = s[0] == '-';
            s = s.Slice(1);
        }
        if (s.Equals("nan", StringComparison.OrdinalIgnoreCase))
        {
            result = MakeNaN(prec);
            return true;
        }
        if (s.Equals("inf", StringComparison.OrdinalIgnoreCase) ||
            s.Equals("infinity", StringComparison.OrdinalIgnoreCase))
        {
            result = MakeInf(neg, prec);
            return true;
        }

        // [digits][.digits][e±exp]
        int ePos = s.IndexOfAny('e', 'E');
        long exp10 = 0;
        if (ePos >= 0)
        {
            if (!long.TryParse(s.Slice(ePos + 1), NumberStyles.AllowLeadingSign, CultureInfo.InvariantCulture, out exp10))
            {
                return false;
            }
            s = s.Slice(0, ePos);
        }
        int dot = s.IndexOf('.');
        BigInt mant;
        if (dot >= 0)
        {
            int fracDigits = s.Length - dot - 1;
            Span<char> merged = s.Length <= 512 ? stackalloc char[s.Length - 1] : new char[s.Length - 1];
            s.Slice(0, dot).CopyTo(merged);
            s.Slice(dot + 1).CopyTo(merged.Slice(dot));
            if (!BigInt.TryParse(merged, 10, out mant))
            {
                return false;
            }
            exp10 -= fracDigits;
        }
        else if (!BigInt.TryParse(s, 10, out mant))
        {
            return false;
        }

        if (mant.IsZero)
        {
            result = MakeZero(neg, prec);
            return true;
        }

        // value = mant * 10^exp10 = mant * 5^exp10 * 2^exp10
        if (exp10 >= 0)
        {
            BigInt m = mant * BigInt.Pow(5, (ulong)exp10);
            result = Round(m, exp10, sticky: false, neg, prec, rounding, out _);
        }
        else
        {
            BigInt den = BigInt.Pow(5, (ulong)(-exp10));
            long shift = Math.Max(0, prec + 3 + den.GetBitLength() - mant.GetBitLength());
            BigInt.DivRem(mant.ShiftLeft(shift), den, out BigInt q, out BigInt r);
            result = Round(q, exp10 - shift, sticky: !r.IsZero, neg, prec, rounding, out _);
        }
        return true;
    }

    // ---------------------------------------------------------------
    // Integer conversions
    // ---------------------------------------------------------------

    /// <summary>Converts to long, rounding per mode; throws on NaN/Inf/overflow (mpfr_get_si-ish).</summary>
    public long ToInt64(RealRounding rounding = RealRounding.TowardZero) => (long)ToBigInt(rounding);
}
