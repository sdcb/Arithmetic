using System.Buffers;
using System.Diagnostics;
using System.Runtime.CompilerServices;
using Sdcb.Arithmetic.Core;

namespace Sdcb.Arithmetic;

/// <summary>
/// Limb-level helpers that let Real work directly on the Nat limb kernels: intermediate
/// products, quotients and aligned sums live in pooled/stack scratch, and only the final
/// rounded significand is allocated.
/// </summary>
public readonly partial struct Real
{
    /// <summary>Scratch buffers up to this many limbs are stack allocated.</summary>
    private const int StackLimbs = 64;

    /// <summary>Precision (bits) from which Sqrt rounds a remainder-free approximate root.</summary>
    private static readonly int SqrtApprMinBits = int.TryParse(Environment.GetEnvironmentVariable("SDCB_ARITHMETIC_SQRT_APPR"), out int v) ? v : 300;

    /// <summary>A scratch buffer: the caller's stack span when large enough, else pooled. Dispose to return it.</summary>
    private readonly ref struct Scratch
    {
        private readonly ulong[]? _rented;
        public readonly Span<ulong> Span;

        public Scratch(Span<ulong> stack, int n)
        {
            if (n <= stack.Length)
            {
                _rented = null;
                Span = stack.Slice(0, n);
            }
            else
            {
                _rented = ArrayPool<ulong>.Shared.Rent(n);
                Span = _rented.AsSpan(0, n);
            }
        }

        public void Dispose()
        {
            if (_rented is not null)
            {
                ArrayPool<ulong>.Shared.Return(_rented);
            }
        }
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static bool TestBit(ReadOnlySpan<ulong> m, long bit)
        => ((m[(int)(bit >> 6)] >> (int)(bit & 63)) & 1) != 0;

    /// <summary>Number of limbs needed for a value of the given bit length.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static int LimbsFor(long bits) => checked((int)((bits + 63) >> 6));

    /// <summary>
    /// dst[0..dn) = src[0..sn) &lt;&lt; shift with zero fill below and above; dn must hold the
    /// shifted value.
    /// </summary>
    private static void ShiftLeftInto(Span<ulong> dst, int dn, ReadOnlySpan<ulong> src, int sn, long shift)
    {
        Debug.Assert(shift >= 0);
        int limbShift = checked((int)(shift >> 6)), bitShift = (int)(shift & 63);
        Debug.Assert(limbShift + sn <= dn);
        dst.Slice(0, limbShift).Clear();
        int top = limbShift + sn;
        ulong carry = 0;
        if (bitShift == 0)
        {
            src.Slice(0, sn).CopyTo(dst.Slice(limbShift));
        }
        else
        {
            carry = Nat.ShiftLeft(dst.Slice(limbShift), src, sn, bitShift);
        }
        if (top < dn)
        {
            dst[top] = carry;
            dst.Slice(top + 1, dn - top - 1).Clear();
        }
        else
        {
            Debug.Assert(carry == 0);
        }
    }

    /// <summary>Truncating right shift of a magnitude into a freshly allocated BigInt with the given sign.</summary>
    private static BigInt ShiftRightToBigInt(ReadOnlySpan<ulong> m, int mn, long shift, int sign)
    {
        long limbShiftL = shift >> 6;
        if (limbShiftL >= mn)
        {
            return default;
        }
        int limbShift = (int)limbShiftL, bitShift = (int)(shift & 63);
        int sn = mn - limbShift;
        ulong[] q = BigInt.NewLimbs(sn);
        if (bitShift == 0)
        {
            m.Slice(limbShift, sn).CopyTo(q);
        }
        else
        {
            Nat.ShiftRight(q, m.Slice(limbShift), sn, bitShift);
        }
        return BigInt.FromMagnitude(q, sn, sign);
    }

    /// <summary>r[0..an+bn) = |a| * |b| (dedicated squaring when a and b share storage).</summary>
    private static void MulInto(Span<ulong> r, BigInt a, BigInt b)
    {
        ReadOnlySpan<ulong> la = a.Limbs, lb = b.Limbs;
        if (la == lb)
        {
            Nat.Square(r, la, la.Length);
        }
        else if (la.Length >= lb.Length)
        {
            Nat.Mul(r, la, la.Length, lb, lb.Length);
        }
        else
        {
            Nat.Mul(r, lb, lb.Length, la, la.Length);
        }
    }

    /// <summary>
    /// The rounding core on a raw magnitude (high zero limbs allowed; value nonzero):
    /// sign * (m + θ) * 2^e, θ ∈ (0,1) iff sticky, rounded to prec bits with a single
    /// allocation for the result significand.
    /// </summary>
    internal static Real RoundLimbs(ReadOnlySpan<ulong> m, long e, bool sticky, bool negative, int prec, RealRounding rnd, out int ternary)
    {
        int mn = Nat.TrimmedLength(m, m.Length);
        Debug.Assert(mn > 0);
        if (prec == 0)
        {
            prec = DefaultPrecision;
        }
        long bits = Nat.BitLength(m, mn);
        if (bits <= 128 && prec <= 126)
        {
            UInt128 v = mn == 1 ? m[0] : new UInt128(m[1], m[0]);
            return RoundSmall(v, e, sticky, negative, prec, rnd, out ternary);
        }
        int qn = LimbsFor(prec);
        long drop = bits - prec;
        if (drop <= 0)
        {
            ulong[] w = BigInt.NewLimbs(qn);
            ShiftLeftInto(w, qn, m, mn, -drop);
            if (sticky)
            {
                return RoundSticky(w, qn, e + drop, negative, prec, rnd, out ternary);
            }
            ternary = 0;
            return new Real(new BigInt(w, qn), e + drop, prec, Kind.Normal, negative);
        }

        bool guard = TestBit(m, drop - 1);
        bool stickyLow = sticky || (drop > 1 && Nat.TrailingZeroCount(m, mn) < drop - 1);
        int limbShift = (int)(drop >> 6), bitShift = (int)(drop & 63);
        int sn = mn - limbShift;
        Debug.Assert(sn >= qn);
        ulong[] q = BigInt.NewLimbs(sn);
        if (bitShift == 0)
        {
            m.Slice(limbShift, sn).CopyTo(q);
        }
        else
        {
            Nat.ShiftRight(q, m.Slice(limbShift), sn, bitShift);
        }
        e += drop;
        ternary = 0;
        if (guard || stickyLow)
        {
            bool roundUp = rnd switch
            {
                RealRounding.NearestEven => guard && (stickyLow || (q[0] & 1) != 0),
                RealRounding.TowardZero => false,
                RealRounding.AwayFromZero => true,
                RealRounding.TowardPositive => !negative,
                RealRounding.TowardNegative => negative,
                _ => throw new ArgumentOutOfRangeException(nameof(rnd)),
            };
            if (roundUp)
            {
                IncrementSignificand(q, qn, prec, ref e);
                ternary = negative ? -1 : 1;
            }
            else
            {
                ternary = negative ? 1 : -1;
            }
        }
        return new Real(new BigInt(q, qn), e, prec, Kind.Normal, negative);
    }

    /// <summary>
    /// A prec-bit significand w followed by a sticky θ ∈ (0,1) with no guard bit: matches the
    /// BigInt rounding core, where nearest-even only rounds up when the guard bit is set.
    /// </summary>
    private static Real RoundSticky(ulong[] w, int qn, long e, bool negative, int prec, RealRounding rnd, out int ternary)
    {
        bool roundUp = rnd switch
        {
            RealRounding.NearestEven => false,
            RealRounding.TowardZero => false,
            RealRounding.AwayFromZero => true,
            RealRounding.TowardPositive => !negative,
            RealRounding.TowardNegative => negative,
            _ => throw new ArgumentOutOfRangeException(nameof(rnd)),
        };
        if (roundUp)
        {
            IncrementSignificand(w, qn, prec, ref e);
            ternary = negative ? -1 : 1;
        }
        else
        {
            ternary = negative ? 1 : -1;
        }
        return new Real(new BigInt(w, qn), e, prec, Kind.Normal, negative);
    }

    /// <summary>q[0..qn) += 1 for a prec-bit significand, renormalizing 2^prec to 2^(prec-1).</summary>
    private static void IncrementSignificand(Span<ulong> q, int qn, int prec, ref long e)
    {
        ulong carry = Nat.AddScalarInPlace(q, qn, 1);
        int topBits = prec & 63;
        if (carry != 0 || (topBits != 0 && (q[qn - 1] >> topBits) != 0))
        {
            q.Slice(0, qn).Clear();
            q[qn - 1] = 1UL << ((prec - 1) & 63);
            e++;
        }
    }

    /// <summary>
    /// CanRound on limbs without allocating: with F = m mod 2^drop and 1 &lt;= err &lt;= 2^62,
    /// rounding is ambiguous when F &lt; err, F &gt;= 2^drop - err, or (nearest) |F - 2^(drop-1)| &lt;= err.
    /// </summary>
    private static bool CanRoundLimbs(ReadOnlySpan<ulong> m, long drop, ulong err, RealRounding rnd)
    {
        Debug.Assert(drop >= 2 && err >= 1 && err <= 1UL << 62 && Nat.BitLength(m, m.Length) > drop);
        if (drop <= 65)
        {
            UInt128 low = m.Length > 1 ? new UInt128(m[1], m[0]) : m[0];
            UInt128 window = UInt128.One << (int)drop;
            UInt128 f = low & (window - 1);
            if (f < err || window - f <= err)
            {
                return false;
            }
            if (rnd == RealRounding.NearestEven)
            {
                UInt128 half = window >> 1;
                UInt128 diff = f >= half ? f - half : half - f;
                if (diff <= err)
                {
                    return false;
                }
            }
            return true;
        }
        // F = top * 2^(drop-1) + mid * 2^64 + l0, where mid spans bits [64, drop-1).
        ulong l0 = m[0];
        bool top = TestBit(m, drop - 1);
        bool midZero = BitRangeIs(m, 64, drop - 1, ones: false);
        bool midOnes = !midZero && BitRangeIs(m, 64, drop - 1, ones: true);
        ulong negErr = 0UL - err; // 2^64 - err
        if (!top && midZero && l0 < err)
        {
            return false;
        }
        if (top && midOnes && l0 >= negErr)
        {
            return false;
        }
        if (rnd == RealRounding.NearestEven)
        {
            if ((!top && midOnes && l0 >= negErr) || (top && midZero && l0 <= err))
            {
                return false;
            }
        }
        return true;
    }

    /// <summary>True when bits [lo, hi) of m are all zero (ones = false) or all one (ones = true).</summary>
    private static bool BitRangeIs(ReadOnlySpan<ulong> m, long lo, long hi, bool ones)
    {
        Debug.Assert(lo < hi && (hi - 1) >> 6 < m.Length);
        int i0 = (int)(lo >> 6), i1 = (int)((hi - 1) >> 6);
        ulong want = ones ? ulong.MaxValue : 0;
        for (int i = i0; i <= i1; i++)
        {
            ulong mask = ulong.MaxValue;
            if (i == i0)
            {
                mask &= ulong.MaxValue << (int)(lo & 63);
            }
            if (i == i1 && (hi & 63) != 0)
            {
                mask &= (1UL << (int)(hi & 63)) - 1;
            }
            if ((m[i] & mask) != (want & mask))
            {
                return false;
            }
        }
        return true;
    }

    // ---------------------------------------------------------------
    // Ring kernels on scratch buffers
    // ---------------------------------------------------------------

    /// <summary>Round(|a| * |b| * 2^e) with the product in scratch.</summary>
    [SkipLocalsInit]
    private static Real MulRound(BigInt a, BigInt b, long e, bool negative, int prec, RealRounding rnd, out int ternary)
    {
        int rn = a.LimbCount + b.LimbCount;
        using var s = new Scratch(stackalloc ulong[StackLimbs], rn);
        MulInto(s.Span, a, b);
        return RoundLimbs(s.Span, e, sticky: false, negative, prec, rnd, out ternary);
    }

    /// <summary>Round((|a| &lt;&lt; s) / |b| * 2^e) with numerator, quotient and remainder in scratch.</summary>
    [SkipLocalsInit]
    private static Real DivRound(BigInt a, long s, BigInt b, long e, bool negative, int prec, RealRounding rnd, out int ternary)
    {
        int an = a.LimbCount, dn = b.LimbCount;
        if (dn >= Nat.BZThreshold)
        {
            // One extra quotient limb lets the remainder-free quotient q̂ ∈ [q, q+2] decide the
            // rounding: the exact quotient lies in [q̂-2, q̂+1), so when q̂'s low limb L >= 3 it
            // never reaches a multiple of B and q̂ >> 64 with a sticky bit rounds identically.
            long s2 = s + 64;
            int nn2 = LimbsFor(a.GetBitLength() + s2);
            int qn2 = nn2 - dn + 1;
            using var sa = new Scratch(stackalloc ulong[StackLimbs], nn2 + qn2);
            Span<ulong> num2 = sa.Span.Slice(0, nn2), q2 = sa.Span.Slice(nn2, qn2);
            ShiftLeftInto(num2, nn2, a.Limbs, an, s2);
            Nat.DivApprox(q2, num2, nn2, b.Limbs, dn);
            if (q2[0] >= 3)
            {
                return RoundLimbs(q2.Slice(1), e, sticky: true, negative, prec, rnd, out ternary);
            }
        }
        int nn = LimbsFor(a.GetBitLength() + s);
        Debug.Assert(nn >= dn);
        int qn = nn - dn + 1;
        int total = nn + qn + dn;
        using var sc = new Scratch(stackalloc ulong[StackLimbs], total);
        Span<ulong> num = sc.Span.Slice(0, nn), q = sc.Span.Slice(nn, qn), r = sc.Span.Slice(nn + qn, dn);
        ShiftLeftInto(num, nn, a.Limbs, an, s);
        Nat.DivRem(q, r, num, nn, b.Limbs, dn);
        bool sticky = !Nat.IsZero(r, dn);
        return RoundLimbs(q, e, sticky, negative, prec, rnd, out ternary);
    }

    /// <summary>Round(sqrt(|a| &lt;&lt; s) * 2^e) with the shifted operand and root in scratch.</summary>
    [SkipLocalsInit]
    private static Real SqrtRound(BigInt a, long s, long e, int prec, RealRounding rnd, out int ternary)
    {
        if (prec >= SqrtApprMinBits)
        {
            // One extra root limb: S ∈ [r, r+1] for r = floor(root), so the exact root lies in
            // [S-1, S+1); with S's low limb L >= 2 it cannot reach a multiple of B and
            // S >> 64 with a sticky bit rounds identically, without the remainder.
            long s2 = s + 128;
            int nn2 = LimbsFor(a.GetBitLength() + s2);
            int sn2 = (nn2 + 1) / 2;
            using var sa = new Scratch(stackalloc ulong[StackLimbs], nn2 + sn2);
            Span<ulong> num2 = sa.Span.Slice(0, nn2), root2 = sa.Span.Slice(nn2, sn2);
            ShiftLeftInto(num2, nn2, a.Limbs, a.LimbCount, s2);
            Nat.SqrtApprox(root2, num2, nn2);
            if (root2[0] >= 2)
            {
                return RoundLimbs(root2.Slice(1), e, sticky: true, negative: false, prec, rnd, out ternary);
            }
        }
        int nn = LimbsFor(a.GetBitLength() + s);
        int sn = (nn + 1) / 2;
        using var sc = new Scratch(stackalloc ulong[StackLimbs], nn + sn + sn + 1);
        Span<ulong> num = sc.Span.Slice(0, nn), root = sc.Span.Slice(nn, sn), rem = sc.Span.Slice(nn + sn, sn + 1);
        ShiftLeftInto(num, nn, a.Limbs, a.LimbCount, s);
        int remSize = Nat.SqrtRem(root, rem, num, nn);
        return RoundLimbs(root, e, sticky: remSize != 0, negative: false, prec, rnd, out ternary);
    }

    /// <summary>
    /// Exact signed sum (-1)^aNeg |a| 2^ea + (-1)^bNeg |b| 2^eb, rounded; both operands are
    /// aligned into scratch (the caller bounds the exponent gap).
    /// </summary>
    [SkipLocalsInit]
    private static Real AddRound(BigInt a, long ea, bool aNeg, BigInt b, long eb, bool bNeg, int prec, RealRounding rnd, out int ternary)
    {
        long f = Math.Min(ea, eb);
        long sa = ea - f, sb = eb - f;
        int n = Math.Max(LimbsFor(a.GetBitLength() + sa), LimbsFor(b.GetBitLength() + sb)) + 1;
        using var sc = new Scratch(stackalloc ulong[StackLimbs], 2 * n);
        Span<ulong> xa = sc.Span.Slice(0, n), xb = sc.Span.Slice(n, n);
        ShiftLeftInto(xa, n, a.Limbs, a.LimbCount, sa);
        ShiftLeftInto(xb, n, b.Limbs, b.LimbCount, sb);
        if (aNeg == bNeg)
        {
            Nat.Add(xa, xa, xb, n);
            return RoundLimbs(xa, f, sticky: false, aNeg, prec, rnd, out ternary);
        }
        int c = Nat.Compare(xa, xb, n);
        if (c == 0)
        {
            ternary = 0;
            return MakeZero(rnd == RealRounding.TowardNegative, prec);
        }
        if (c > 0)
        {
            Nat.Sub(xa, xa, xb, n);
            return RoundLimbs(xa, f, sticky: false, aNeg, prec, rnd, out ternary);
        }
        Nat.Sub(xb, xb, xa, n);
        return RoundLimbs(xb, f, sticky: false, bNeg, prec, rnd, out ternary);
    }

    /// <summary>Round((|a| &lt;&lt; sh) ± θ) for the far-apart addition case: a sticky unit just below a's extended lsb.</summary>
    [SkipLocalsInit]
    private static Real AddStickyRound(BigInt a, long sh, long e, bool subtract, bool negative, int prec, RealRounding rnd, out int ternary)
    {
        int n = LimbsFor(a.GetBitLength() + sh);
        using var sc = new Scratch(stackalloc ulong[StackLimbs], n);
        ShiftLeftInto(sc.Span, n, a.Limbs, a.LimbCount, sh);
        if (subtract)
        {
            Nat.SubScalarInPlace(sc.Span, n, 1);
        }
        return RoundLimbs(sc.Span, e, sticky: true, negative, prec, rnd, out ternary);
    }

    // ---------------------------------------------------------------
    // Fixed-point kernels (value * 2^wp as a signed BigInt)
    // ---------------------------------------------------------------

    [SkipLocalsInit]
    private static BigInt FixMulCore(BigInt a, BigInt b, int wp)
    {
        if (a.IsZero || b.IsZero)
        {
            return default;
        }
        int sign = a.Sign * b.Sign;
        int rn = a.LimbCount + b.LimbCount;
        using var s = new Scratch(stackalloc ulong[StackLimbs], rn);
        MulInto(s.Span, a, b);
        return ShiftRightToBigInt(s.Span, rn, wp, sign);
    }

    [SkipLocalsInit]
    private static BigInt FixDivCore(BigInt a, BigInt b, int wp)
    {
        if (b.IsZero)
        {
            throw new DivideByZeroException();
        }
        if (a.IsZero)
        {
            return default;
        }
        int dn = b.LimbCount;
        int nn = LimbsFor(a.GetBitLength() + wp);
        if (nn < dn)
        {
            return a.ShiftLeft(wp) / b;
        }
        int qn = nn - dn + 1;
        ulong[] q = BigInt.NewLimbs(qn);
        using var sc = new Scratch(stackalloc ulong[StackLimbs], nn + dn);
        Span<ulong> num = sc.Span.Slice(0, nn), r = sc.Span.Slice(nn, dn);
        ShiftLeftInto(num, nn, a.Limbs, a.LimbCount, wp);
        Nat.DivRem(q, r, num, nn, b.Limbs, dn);
        return BigInt.FromMagnitude(q, qn, a.Sign * b.Sign);
    }

    [SkipLocalsInit]
    private static BigInt FixSqrtCore(BigInt a, int wp)
    {
        if (a.Sign <= 0)
        {
            return BigInt.Sqrt(a.ShiftLeft(wp));
        }
        int nn = LimbsFor(a.GetBitLength() + wp);
        int sn = (nn + 1) / 2;
        ulong[] root = BigInt.NewLimbs(sn);
        using var sc = new Scratch(stackalloc ulong[StackLimbs], nn);
        ShiftLeftInto(sc.Span, nn, a.Limbs, a.LimbCount, wp);
        Nat.SqrtRem(root, default, sc.Span, nn);
        return BigInt.FromMagnitude(root, sn, 1);
    }

    /// <summary>
    /// A mutable signed fixed-point accumulator for the Horner steps of SeriesRect:
    /// multiply by a word, truncating division by a word and signed addition all run in
    /// place on one pooled buffer, matching the BigInt truncation semantics exactly.
    /// </summary>
    private sealed class FixAccumulator : IDisposable
    {
        private ulong[] _buf;
        private int _n;
        private bool _neg;

        public FixAccumulator(int capacity)
        {
            _buf = ArrayPool<ulong>.Shared.Rent(Math.Max(4, capacity));
        }

        public bool IsZero => _n == 0;

        private void Ensure(int n)
        {
            if (_buf.Length < n)
            {
                ulong[] nb = ArrayPool<ulong>.Shared.Rent(Math.Max(n, 2 * _buf.Length));
                _buf.AsSpan(0, _n).CopyTo(nb);
                ArrayPool<ulong>.Shared.Return(_buf);
                _buf = nb;
            }
        }

        public void Set(BigInt v)
        {
            ReadOnlySpan<ulong> l = v.Limbs;
            Ensure(l.Length);
            l.CopyTo(_buf);
            _n = l.Length;
            _neg = v.Sign < 0;
        }

        public void MulSmall(ulong m)
        {
            if (_n == 0)
            {
                return;
            }
            Ensure(_n + 1);
            ulong hi = Nat.MulScalar(_buf, _buf, _n, m);
            if (hi != 0)
            {
                _buf[_n++] = hi;
            }
            else
            {
                _n = Nat.TrimmedLength(_buf, _n);
            }
        }

        public void DivSmall(ulong d)
        {
            if (_n == 0)
            {
                return;
            }
            Nat.DivRemScalar(_buf, _buf, _n, d);
            _n = Nat.TrimmedLength(_buf, _n);
        }

        public void Add(BigInt v)
        {
            ReadOnlySpan<ulong> l = v.Limbs;
            int vn = l.Length;
            if (vn == 0)
            {
                return;
            }
            bool vNeg = v.Sign < 0;
            if (_n == 0)
            {
                Set(v);
                return;
            }
            int n = Math.Max(_n, vn);
            Ensure(n + 1);
            Span<ulong> buf = _buf;
            if (_n < n)
            {
                buf.Slice(_n, n - _n).Clear();
            }
            if (vNeg == _neg)
            {
                ulong carry = Nat.AddInPlace(buf, n, l, vn);
                _n = n;
                if (carry != 0)
                {
                    buf[_n++] = carry;
                }
                return;
            }
            int c = Nat.CompareNormalized(buf, Nat.TrimmedLength(buf, n), l, vn);
            if (c == 0)
            {
                _n = 0;
                _neg = false;
                return;
            }
            if (c > 0)
            {
                Nat.SubInPlace(buf, n, l, vn);
            }
            else
            {
                // |v| > |acc|: acc = |v| - acc, taking v's sign.
                if (vn == n)
                {
                    Nat.Sub(buf, l, buf, n);
                }
                else
                {
                    using var sc = new Scratch(stackalloc ulong[StackLimbs], n);
                    l.CopyTo(sc.Span);
                    sc.Span.Slice(vn, n - vn).Clear();
                    Nat.Sub(buf, sc.Span, buf, n);
                }
                _neg = vNeg;
            }
            _n = Nat.TrimmedLength(buf, n);
        }

        public BigInt ToBigInt()
        {
            if (_n == 0)
            {
                return default;
            }
            ulong[] d = BigInt.NewLimbs(_n);
            _buf.AsSpan(0, _n).CopyTo(d);
            return BigInt.FromMagnitude(d, _n, _neg ? -1 : 1);
        }

        public void Dispose() => ArrayPool<ulong>.Shared.Return(_buf);
    }
}
