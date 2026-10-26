using System.Diagnostics;

namespace Sdcb.Arithmetic;

public readonly partial struct Real
{
    // ---------------------------------------------------------------
    // Mathematical constants with per-precision caching
    // (mpfr_const_pi / mpfr_const_log2 and EulerGamma's number).
    // Series are evaluated in fixed point with single-limb divisors,
    // giving O(n²/log q) with tiny constants — plenty fast to ~100k bits.
    // ---------------------------------------------------------------

    private static readonly object ConstLock = new();
    private static BigInt s_piCache;    // pi * 2^s_piBits
    private static int s_piBits;
    private static BigInt s_ln2Cache;
    private static int s_ln2Bits;
    private static BigInt s_eCache;
    private static int s_eBits;

    /// <summary>π correctly rounded to the given precision (mpfr_const_pi).</summary>
    public static Real Pi(int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = precision != 0 ? precision : DefaultPrecision;
        return Ziv(prec, rounding, wp =>
        {
            BigInt pi = PiFixed(wp);
            return (pi, (long)-wp, 64);
        });
    }

    /// <summary>ln 2 correctly rounded (mpfr_const_log2).</summary>
    public static Real Ln2(int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = precision != 0 ? precision : DefaultPrecision;
        return Ziv(prec, rounding, wp => (Ln2Fixed(wp), (long)-wp, 64));
    }

    /// <summary>EulerGamma's number e correctly rounded.</summary>
    public static Real E(int precision = 0, RealRounding rounding = RealRounding.NearestEven)
    {
        int prec = precision != 0 ? precision : DefaultPrecision;
        return Ziv(prec, rounding, wp => (EFixed(wp), (long)-wp, 64));
    }

    /// <summary>π at wp fractional bits (cached).</summary>
    internal static BigInt PiFixed(int wp)
    {
        lock (ConstLock)
        {
            if (s_piBits >= wp)
            {
                return s_piCache.ShiftRightTruncate(s_piBits - wp);
            }
        }
        int cp = wp + 64;
        BigInt pi;
        if (wp >= ConstAgmThreshold)
        {
            pi = PiFixedChudnovsky(cp);
        }
        else
        {
            // Machin: pi = 16*atan(1/5) - 4*atan(1/239)
            pi = (AtanRecipFixed(5, cp) << 4) - (AtanRecipFixed(239, cp) << 2);
        }
        lock (ConstLock)
        {
            if (cp > s_piBits)
            {
                s_piCache = pi;
                s_piBits = cp;
            }
        }
        return pi.ShiftRightTruncate(cp - wp);
    }

    /// <summary>ln2 at wp fractional bits (cached). ln2 = 2*atanh(1/3).</summary>
    internal static BigInt Ln2Fixed(int wp)
    {
        lock (ConstLock)
        {
            if (s_ln2Bits >= wp)
            {
                return s_ln2Cache.ShiftRightTruncate(s_ln2Bits - wp);
            }
        }
        int cp = wp + 64;
        BigInt v;
        if (wp >= ConstAgmThreshold)
        {
            // ln2 = (2/3) Σ 1/((2k+1) 9^k) by binary splitting (O(M(n) log n)).
            // Term ratio t_{k+1}/t_k = (2k+1) / (9 (2k+3)); t_0 = 1.
            long n = (cp + 80) * 10 / 32 + 2; // 9^k > 2^{cp+80}
            (BigInt _, BigInt q, BigInt t) = BinSplit(0, n, k => (2 * k + 1, 9 * (2 * k + 3)));
            v = (t << (cp + 1)) / (q * 3);
        }
        else
        {
            v = AtanhRecipFixed(3, cp) << 1;
        }
        lock (ConstLock)
        {
            if (cp > s_ln2Bits)
            {
                s_ln2Cache = v;
                s_ln2Bits = cp;
            }
        }
        return v.ShiftRightTruncate(cp - wp);
    }

    /// <summary>e at wp fractional bits (cached). e = Σ 1/k!.</summary>
    internal static BigInt EFixed(int wp)
    {
        lock (ConstLock)
        {
            if (s_eBits >= wp)
            {
                return s_eCache.ShiftRightTruncate(s_eBits - wp);
            }
        }
        int cp = wp + 64;
        BigInt sum;
        if (wp >= ConstAgmThreshold)
        {
            // e = Σ 1/k! by binary splitting: term ratio 1/(k+1), t_0 = 1.
            long n = 1;
            double lg = 0;
            while (-lg < cp + 80)
            {
                n++;
                lg += Math.Log2(n); // log2(N!) > cp + 80
            }
            (BigInt _, BigInt q, BigInt t) = BinSplit(0, n, static k => (1, k + 1));
            sum = (t << cp) / q;
        }
        else
        {
            BigInt term = BigInt.One << cp;
            sum = term;
            for (ulong k = 1; !term.IsZero; k++)
            {
                term = BigInt.DivRem(term, k, out _);
                sum += term;
            }
        }
        lock (ConstLock)
        {
            if (cp > s_eBits)
            {
                s_eCache = sum;
                s_eBits = cp;
            }
        }
        return sum.ShiftRightTruncate(cp - wp);
    }

    /// <summary>
    /// Binary splitting for Σ t_k with t_{k+1}/t_k = p_k/q_k over [a, b): returns products
    /// P = Π p, Q = Π q and T such that Σ (t_k / t_a) = T / Q.
    /// </summary>
    private static (BigInt P, BigInt Q, BigInt T) BinSplit(long a, long b, Func<long, (long p, long q)> ratio)
    {
        if (b - a == 1)
        {
            (long p, long q) = ratio(a);
            return (new BigInt(p), new BigInt(q), new BigInt(q));
        }
        long m = (a + b) / 2;
        (BigInt p1, BigInt q1, BigInt t1) = BinSplit(a, m, ratio);
        (BigInt p2, BigInt q2, BigInt t2) = BinSplit(m, b, ratio);
        return (p1 * p2, q1 * q2, t1 * q2 + p1 * t2);
    }

    /// <summary>
    /// Precision above which the constants use binary splitting (Chudnovsky for π)
    /// instead of the linear fixed-point series (configurable for tuning).
    /// </summary>
    internal static int ConstAgmThreshold = int.TryParse(Environment.GetEnvironmentVariable("SDCB_ARITHMETIC_CONSTAGM"), out int t) ? t : 2000;

    /// <summary>
    /// π at wp fractional bits by the Chudnovsky series with binary splitting:
    /// 1/π = 12 Σ (-1)^k (6k)! (13591409 + 545140134 k) / ((3k)! (k!)³ 640320^(3k+3/2)),
    /// i.e. π = 426880 √10005 · Q / T. Each term adds ~47.11 bits, so the cost is
    /// O(M(n) log² n) — several times faster than the Brent–Salamin AGM, whose ~log2(n)
    /// full-precision square roots dominate at high precision.
    /// </summary>
    private static BigInt PiFixedChudnovsky(int wp)
    {
        int p = wp + 32;
        long n = (long)(wp / 47.11) + 2;
        (_, BigInt q, BigInt t) = ChudnovskySplit(0, n, needP: false);
        Real r = Divide(FromBigInt(q, 0, p), FromBigInt(t, 0, p), p);
        Real s = Sqrt(FromBigInt(10005, 0, p), p);
        Real pi = Multiply(Multiply(r, s, p), FromBigInt(426880, 0, p), p);
        return pi.ToFixed(wp);
    }

    /// <summary>
    /// Chudnovsky binary splitting over [a, b): P = Π p_k, Q = Π q_k and T with
    /// Σ_{k∈[a,b)} term_k / term_a = T / Q (P is skipped for the rightmost branch).
    /// </summary>
    private static (BigInt P, BigInt Q, BigInt T) ChudnovskySplit(long a, long b, bool needP)
    {
        if (b - a == 1)
        {
            if (a == 0)
            {
                return (BigInt.One, BigInt.One, new BigInt(13591409));
            }
            BigInt p = (BigInt)(6 * a - 5) * (2 * a - 1) * (6 * a - 1);
            BigInt q = (BigInt)a * a * a * 10939058860032000L; // 640320³ / 24
            BigInt t = p * (13591409 + 545140134 * (BigInt)a);
            return (p, q, (a & 1) == 1 ? -t : t);
        }
        long m = (a + b) / 2;
        (BigInt p1, BigInt q1, BigInt t1) = ChudnovskySplit(a, m, needP: true);
        (BigInt p2, BigInt q2, BigInt t2) = ChudnovskySplit(m, b, needP);
        return (needP ? p1 * p2 : default, q1 * q2, t1 * q2 + p1 * t2);
    }
    /// <summary>atan(1/q) at wp fractional bits: Σ (-1)^i / ((2i+1) q^(2i+1)).</summary>
    internal static BigInt AtanRecipFixed(ulong q, int wp)
    {
        Debug.Assert(q >= 2);
        BigInt power = BigInt.DivRem(BigInt.One << wp, q, out _); // (1/q) * 2^wp
        ulong q2 = q * q;
        BigInt sum = power;
        bool subtract = true;
        for (ulong i = 1; ; i++, subtract = !subtract)
        {
            power = BigInt.DivRem(power, q2, out _);
            if (power.IsZero)
            {
                break;
            }
            BigInt term = BigInt.DivRem(power, 2 * i + 1, out _);
            sum = subtract ? sum - term : sum + term;
        }
        return sum;
    }

    /// <summary>atanh(1/q) at wp fractional bits: Σ 1 / ((2i+1) q^(2i+1)).</summary>
    internal static BigInt AtanhRecipFixed(ulong q, int wp)
    {
        BigInt power = BigInt.DivRem(BigInt.One << wp, q, out _);
        ulong q2 = q * q;
        BigInt sum = power;
        for (ulong i = 1; ; i++)
        {
            power = BigInt.DivRem(power, q2, out _);
            if (power.IsZero)
            {
                break;
            }
            sum += BigInt.DivRem(power, 2 * i + 1, out _);
        }
        return sum;
    }
}
