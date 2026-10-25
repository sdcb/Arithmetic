using System.Buffers;
using System.Diagnostics;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;

namespace Sdcb.Arithmetic.Core;

/// <summary>2x2 matrix with single-limb entries (each &lt; 2^63), determinant 1.</summary>
internal struct Matrix1
{
    public ulong U00, U01, U10, U11;
}

/// <summary>
/// A pair of multi-limb non-negative values (X, Y) updated as a row vector times
/// determinant-1 matrices. Used for the rows of <see cref="HalfGcdMatrix"/> and for the
/// gcdext cofactor. Buffers are zero beyond N.
/// </summary>
internal sealed class LimbRow
{
    public ulong[] X, Y;
    public int N;

    public LimbRow(int capacity, bool identityX)
    {
        X = new ulong[Math.Max(capacity, 2)];
        Y = new ulong[Math.Max(capacity, 2)];
        (identityX ? X : Y)[0] = 1;
        N = 1;
    }

    public void Ensure(int need)
    {
        if (X.Length < need)
        {
            int cap = Math.Max(need, X.Length * 3 / 2);
            Array.Resize(ref X, cap);
            Array.Resize(ref Y, cap);
        }
    }

    private void Renormalize(int len)
    {
        while (len > 1 && X[len - 1] == 0 && Y[len - 1] == 0)
        {
            len--;
        }
        N = len;
    }

    /// <summary>(X, Y) &lt;- (X, Y) * M.</summary>
    public void MulMatrix1(in Matrix1 m)
    {
        Ensure(N + 1);
        Nat.RowMul1(X, Y, N, m);
        Renormalize(N + 1);
    }

    /// <summary>Column col += q * column (1-col): col 0 means X += q*Y, col 1 means Y += q*X.</summary>
    public void UpdateQ(ReadOnlySpan<ulong> q, int col)
    {
        int qn = Nat.TrimmedLength(q, q.Length);
        if (qn == 0)
        {
            return;
        }
        ulong[] src = col == 0 ? Y : X;
        int sn = Nat.TrimmedLength(src, N);
        if (sn == 0)
        {
            return;
        }
        int pn = sn + qn;
        int len = Math.Max(N, pn);
        Ensure(len + 1);
        src = col == 0 ? Y : X;
        ulong[] dst = col == 0 ? X : Y;
        ulong carry;
        if (qn == 1)
        {
            carry = Nat.AddMulScalar(dst, src, sn, q[0]);
            carry = Nat.AddScalarInPlace(dst.AsSpan(sn), len + 1 - sn, carry);
        }
        else
        {
            ulong[] tmp = ArrayPool<ulong>.Shared.Rent(pn);
            Nat.MulAny(tmp, src, sn, q, qn);
            carry = Nat.AddInPlace(dst, len + 1, tmp, pn);
            ArrayPool<ulong>.Shared.Return(tmp);
        }
        Debug.Assert(carry == 0);
        Renormalize(len + 1);
    }

    /// <summary>(X, Y) &lt;- (X, Y) * M for a multi-limb matrix.</summary>
    public void MulMatrix(HalfGcdMatrix m)
    {
        m.Pad();
        int mn = m.N;
        int cap = N + mn + 2;
        ulong[] nx = new ulong[cap], ny = new ulong[cap];
        ulong[] tmp = ArrayPool<ulong>.Shared.Rent(cap);
        int xn = Nat.MulSum(nx, X, N, m.R0.X, mn, Y, N, m.R1.X, mn, tmp);
        int yn = Nat.MulSum(ny, X, N, m.R0.Y, mn, Y, N, m.R1.Y, mn, tmp);
        ArrayPool<ulong>.Shared.Return(tmp);
        X = nx;
        Y = ny;
        Renormalize(Math.Max(Math.Max(xn, yn), 1));
    }
}

/// <summary>
/// Reduction matrix M = (p00 p01; p10 p11), non-negative entries, determinant 1, such that
/// the original pair equals M times the reduced pair (Möller, 2008).
/// </summary>
internal sealed class HalfGcdMatrix : Nat.SubdivideHook
{
    public readonly LimbRow R0, R1;

    public HalfGcdMatrix(int n)
    {
        int cap = (n + 1) / 2 + 3;
        R0 = new LimbRow(cap, identityX: true);
        R1 = new LimbRow(cap, identityX: false);
    }

    public int N => Math.Max(R0.N, R1.N);

    public void MulMatrix1(in Matrix1 m)
    {
        R0.MulMatrix1(m);
        R1.MulMatrix1(m);
    }

    public void MulMatrix(HalfGcdMatrix m)
    {
        R0.MulMatrix(m);
        R1.MulMatrix(m);
    }

    /// <summary>Makes every entry buffer hold at least N valid (zero-padded) limbs.</summary>
    public void Pad()
    {
        int n = N;
        R0.Ensure(n);
        R1.Ensure(n);
    }

    public override void OnGcd(ReadOnlySpan<ulong> g, bool isA) => throw new InvalidOperationException();

    public override void OnQuotient(ReadOnlySpan<ulong> q, int col)
    {
        R0.UpdateQ(q, col);
        R1.UpdateQ(q, col);
    }
}

internal static partial class Nat
{
    /// <summary>Size (limbs) from which hgcd recurses instead of running Lehmer steps.</summary>
    internal static int HgcdThreshold = EnvInt("SDCB_ARITHMETIC_HGCD") ?? 200;

    /// <summary>Size (limbs) from which gcd/gcdext use the subquadratic hgcd reduction.</summary>
    internal static int GcdDcThreshold = EnvInt("SDCB_ARITHMETIC_GCD_DC") ?? 800;

    /// <summary>Callback for <see cref="GcdSubdivideStep"/>: receives the gcd or each quotient.</summary>
    internal abstract class SubdivideHook
    {
        /// <summary>The gcd was found; isA tells whether it is the current value of the caller's A (else B).</summary>
        public abstract void OnGcd(ReadOnlySpan<ulong> g, bool isA);

        /// <summary>Quotient q was applied: col 0 means B -= q*A, col 1 means A -= q*B.</summary>
        public abstract void OnQuotient(ReadOnlySpan<ulong> q, int col);
    }

    private sealed class GcdCapture : SubdivideHook
    {
        public ulong[]? G;
        public bool IsA;
        public LimbRow? Row;

        public override void OnGcd(ReadOnlySpan<ulong> g, bool isA)
        {
            G = g.ToArray();
            IsA = isA;
        }

        public override void OnQuotient(ReadOnlySpan<ulong> q, int col) => Row?.UpdateQ(q, col);
    }

    // ---------------------------------------------------------------
    // Small helpers
    // ---------------------------------------------------------------

    /// <summary>r = a*b for any operand order (either may be shorter); returns an+bn, or 0 if one is empty.</summary>
    internal static int MulAny(Span<ulong> r, ReadOnlySpan<ulong> a, int an, ReadOnlySpan<ulong> b, int bn)
    {
        if (an == 0 || bn == 0)
        {
            return 0;
        }
        if (an >= bn)
        {
            Mul(r, a, an, b, bn);
        }
        else
        {
            Mul(r, b, bn, a, an);
        }
        return an + bn;
    }

    /// <summary>dst = x*p + y*q (all non-negative); returns the normalized length. tmp needs the same capacity as dst.</summary>
    internal static int MulSum(Span<ulong> dst, ReadOnlySpan<ulong> x, int xn, ReadOnlySpan<ulong> p, int pn,
        ReadOnlySpan<ulong> y, int yn, ReadOnlySpan<ulong> q, int qn, Span<ulong> tmp)
    {
        xn = TrimmedLength(x, xn); pn = TrimmedLength(p, pn);
        yn = TrimmedLength(y, yn); qn = TrimmedLength(q, qn);
        int n1 = MulAny(dst, x, xn, p, pn);
        int n2 = MulAny(tmp, y, yn, q, qn);
        if (n2 == 0)
        {
            return TrimmedLength(dst, n1);
        }
        if (n1 == 0)
        {
            tmp.Slice(0, n2).CopyTo(dst);
            return TrimmedLength(dst, n2);
        }
        ulong c;
        int len;
        if (n1 >= n2)
        {
            c = Add(dst, dst, n1, tmp, n2);
            len = n1;
        }
        else
        {
            c = Add(dst, tmp, n2, dst, n1);
            len = n2;
        }
        dst[len] = c;
        return TrimmedLength(dst, len + 1);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void SubDD(ref ulong h, ref ulong l, ulong h2, ulong l2)
    {
        ulong nl = l - l2;
        h = h - h2 - (l < l2 ? 1UL : 0UL);
        l = nl;
    }

    /// <summary>(n1:n0) / (d1:d0) by shift-and-subtract (quotients are usually tiny); requires n1 &gt;= d1 &gt; 0.</summary>
    private static ulong Div2(out ulong r1, out ulong r0, ulong n1, ulong n0, ulong d1, ulong d0)
    {
        ulong q = 0;
        int dcnt = BitOperations.LeadingZeroCount(d1) - BitOperations.LeadingZeroCount(n1);
        Debug.Assert(dcnt >= 0);
        d1 = (d1 << dcnt) + (d0 >> 1 >> (63 - dcnt));
        d0 <<= dcnt;
        while (true)
        {
            ulong mask = n1 == d1 ? (n0 >= d0 ? ulong.MaxValue : 0) : (n1 > d1 ? ulong.MaxValue : 0);
            q = (q << 1) - mask;
            SubDD(ref n1, ref n0, mask & d1, mask & d0);
            if (dcnt-- == 0)
            {
                break;
            }
            d0 = (d1 << 63) | (d0 >> 1);
            d1 >>= 1;
        }
        r1 = n1;
        r0 = n0;
        return q;
    }

    // ---------------------------------------------------------------
    // hgcd2: double-limb Lehmer step (Möller's algorithm with Jebelean's condition)
    // ---------------------------------------------------------------

    /// <summary>
    /// Given the top two limbs of a and b, computes a reduction matrix M (entries &lt; 2^63)
    /// that is guaranteed to be a correct prefix of the Euclidean sequence of the full
    /// numbers. Returns false when no reduction is possible.
    /// </summary>
    internal static bool HalfGcd2Limb(ulong ah, ulong al, ulong bh, ulong bl, out Matrix1 m)
    {
        m = default;
        ulong u00, u01, u10, u11;
        const ulong Half = 1UL << 32;
        const ulong Limit1 = 1UL << 33;

        if (ah < 2 || bh < 2)
        {
            return false;
        }
        if (ah > bh || (ah == bh && al > bl))
        {
            SubDD(ref ah, ref al, bh, bl);
            if (ah < 2)
            {
                return false;
            }
            u00 = u01 = u11 = 1;
            u10 = 0;
        }
        else
        {
            SubDD(ref bh, ref bl, ah, al);
            if (bh < 2)
            {
                return false;
            }
            u00 = u10 = u11 = 1;
            u01 = 0;
        }

        bool doA = ah >= bh;
        bool singleStartsWithB;
        while (true)
        {
            if (doA)
            {
                Debug.Assert(ah >= bh);
                if (ah == bh)
                {
                    goto Done;
                }
                if (ah < Half)
                {
                    ah = (ah << 32) + (al >> 32);
                    bh = (bh << 32) + (bl >> 32);
                    singleStartsWithB = false;
                    break;
                }
                // a -= q b: multiply M from the right by (1 q; 0 1).
                SubDD(ref ah, ref al, bh, bl);
                if (ah < 2)
                {
                    goto Done;
                }
                if (ah <= bh)
                {
                    u01 += u00;
                    u11 += u10;
                }
                else
                {
                    ulong q = Div2(out ah, out al, ah, al, bh, bl);
                    if (ah < 2)
                    {
                        // a is too small, but q is correct.
                        u01 += q * u00;
                        u11 += q * u10;
                        goto Done;
                    }
                    q++;
                    u01 += q * u00;
                    u11 += q * u10;
                }
            }
            doA = true;

            Debug.Assert(bh >= ah);
            if (ah == bh)
            {
                goto Done;
            }
            if (bh < Half)
            {
                ah = (ah << 32) + (al >> 32);
                bh = (bh << 32) + (bl >> 32);
                singleStartsWithB = true;
                break;
            }
            // b -= q a: multiply M from the right by (1 0; q 1).
            SubDD(ref bh, ref bl, ah, al);
            if (bh < 2)
            {
                goto Done;
            }
            if (bh <= ah)
            {
                u00 += u01;
                u10 += u11;
            }
            else
            {
                ulong q = Div2(out bh, out bl, bh, bl, ah, al);
                if (bh < 2)
                {
                    u00 += q * u01;
                    u10 += q * u11;
                    goto Done;
                }
                q++;
                u00 += q * u01;
                u10 += q * u11;
            }
        }

        // Single-precision phase on the top 64 bits (the low half limb is discarded).
        bool doA1 = !singleStartsWithB;
        while (true)
        {
            if (doA1)
            {
                Debug.Assert(ah >= bh);
                ah -= bh;
                if (ah < Limit1)
                {
                    break;
                }
                if (ah <= bh)
                {
                    u01 += u00;
                    u11 += u10;
                }
                else
                {
                    ulong q = ah / bh;
                    ah -= q * bh;
                    if (ah < Limit1)
                    {
                        u01 += q * u00;
                        u11 += q * u10;
                        break;
                    }
                    q++;
                    u01 += q * u00;
                    u11 += q * u10;
                }
            }
            doA1 = true;

            Debug.Assert(bh >= ah);
            bh -= ah;
            if (bh < Limit1)
            {
                break;
            }
            if (bh <= ah)
            {
                u00 += u01;
                u10 += u11;
            }
            else
            {
                ulong q = bh / ah;
                bh -= q * ah;
                if (bh < Limit1)
                {
                    u00 += q * u01;
                    u10 += q * u11;
                    break;
                }
                q++;
                u00 += q * u01;
                u10 += q * u11;
            }
        }

    Done:
        m.U00 = u00; m.U01 = u01; m.U10 = u10; m.U11 = u11;
        return true;
    }

    // ---------------------------------------------------------------
    // Fused single-pass matrix kernels
    // ---------------------------------------------------------------

    /// <summary>
    /// In place (a; b) &lt;- M^-1 (a; b) = (u11 a - u01 b; u00 b - u10 a) over n limbs in
    /// one pass; both results are known to be non-negative. Returns the normalized size.
    /// </summary>
    internal static int Matrix22Mul1InverseVector(in Matrix1 m, Span<ulong> a, Span<ulong> b, int n)
    {
        if ((uint)n > (uint)a.Length || (uint)n > (uint)b.Length)
        {
            throw new ArgumentOutOfRangeException(nameof(n));
        }
        ref ulong ap = ref MemoryMarshal.GetReference(a);
        ref ulong bp = ref MemoryMarshal.GetReference(b);
        ulong u00 = m.U00, u01 = m.U01, u10 = m.U10, u11 = m.U11;
        long c0 = 0, c1 = 0;
        for (nint i = 0; i < n; i++)
        {
            ulong ai = Unsafe.Add(ref ap, i), bi = Unsafe.Add(ref bp, i);

            // 128-bit signed: u11*ai - u01*bi + c0 (|value| < 2^127).
            ulong h1 = BigMul(u11, ai, out ulong l1);
            ulong h2 = BigMul(u01, bi, out ulong l2);
            ulong lo = l1 - l2;
            ulong hi = h1 - h2 - (l1 < l2 ? 1UL : 0UL);
            ulong lo2 = lo + (ulong)c0;
            hi += (ulong)(c0 >> 63) + (lo2 < lo ? 1UL : 0UL);
            Unsafe.Add(ref ap, i) = lo2;
            c0 = (long)hi;

            h1 = BigMul(u00, bi, out l1);
            h2 = BigMul(u10, ai, out l2);
            lo = l1 - l2;
            hi = h1 - h2 - (l1 < l2 ? 1UL : 0UL);
            lo2 = lo + (ulong)c1;
            hi += (ulong)(c1 >> 63) + (lo2 < lo ? 1UL : 0UL);
            Unsafe.Add(ref bp, i) = lo2;
            c1 = (long)hi;
        }
        Debug.Assert(c0 == 0 && c1 == 0);
        while (n > 0 && Unsafe.Add(ref ap, n - 1) == 0 && Unsafe.Add(ref bp, n - 1) == 0)
        {
            n--;
        }
        return n;
    }

    /// <summary>In place (x, y) &lt;- (x u00 + y u10, x u01 + y u11) over n limbs; writes carries to x[n], y[n].</summary>
    internal static void RowMul1(Span<ulong> x, Span<ulong> y, int n, in Matrix1 m)
    {
        if ((uint)n >= (uint)x.Length || (uint)n >= (uint)y.Length)
        {
            throw new ArgumentOutOfRangeException(nameof(n));
        }
        ref ulong xp = ref MemoryMarshal.GetReference(x);
        ref ulong yp = ref MemoryMarshal.GetReference(y);
        ulong u00 = m.U00, u01 = m.U01, u10 = m.U10, u11 = m.U11;
        ulong cx = 0, cy = 0;
        for (nint i = 0; i < n; i++)
        {
            ulong xi = Unsafe.Add(ref xp, i), yi = Unsafe.Add(ref yp, i);
            // Entries < 2^63 keep each two-product sum plus carry below 2^128.
            ulong h1 = BigMul(xi, u00, out ulong l1);
            ulong h2 = BigMul(yi, u10, out ulong l2);
            ulong lo = l1 + l2;
            ulong hi = h1 + h2 + (lo < l1 ? 1UL : 0UL);
            ulong lo2 = lo + cx;
            Unsafe.Add(ref xp, i) = lo2;
            cx = hi + (lo2 < lo ? 1UL : 0UL);

            h1 = BigMul(xi, u01, out l1);
            h2 = BigMul(yi, u11, out l2);
            lo = l1 + l2;
            hi = h1 + h2 + (lo < l1 ? 1UL : 0UL);
            lo2 = lo + cy;
            Unsafe.Add(ref yp, i) = lo2;
            cy = hi + (lo2 < lo ? 1UL : 0UL);
        }
        Unsafe.Add(ref xp, n) = cx;
        Unsafe.Add(ref yp, n) = cy;
    }

    // ---------------------------------------------------------------
    // Subtraction + division step (used when hgcd2 cannot make progress)
    // ---------------------------------------------------------------

    /// <summary>
    /// One subtraction followed by one division on (a, b) of n limbs (the subdivision step of
    /// Möller's half-GCD). With s == 0, returns 0 once the gcd is found (reported via
    /// the hook). With s &gt; 0, never reduces either value to s limbs or fewer and returns 0
    /// when no such reduction is possible. a and b need n+1 limbs of storage.
    /// </summary>
    internal static int GcdSubdivideStep(Span<ulong> ap, Span<ulong> bp, int n, int s, SubdivideHook hook)
    {
        int an = TrimmedLength(ap, n), bn = TrimmedLength(bp, n);
        bool swapped = false;

        // Arrange a < b, subtract b -= a.
        if (an == bn)
        {
            int c = Compare(ap, bp, an);
            if (c == 0)
            {
                if (s == 0)
                {
                    hook.OnGcd(ap.Slice(0, an), isA: true);
                }
                return 0;
            }
            if (c > 0)
            {
                Swap(ref ap, ref bp, ref swapped);
            }
        }
        else if (an > bn)
        {
            Swap(ref ap, ref bp, ref swapped);
            (an, bn) = (bn, an);
        }
        if (an <= s)
        {
            if (s == 0)
            {
                hook.OnGcd(bp.Slice(0, bn), isA: swapped);
            }
            return 0;
        }

        ulong borrow = Sub(bp, bp, bn, ap, an);
        Debug.Assert(borrow == 0);
        bn = TrimmedLength(bp, bn);
        Debug.Assert(bn > 0);

        if (bn <= s)
        {
            // Undo the subtraction.
            ulong cy = Add(bp, ap, an, bp, bn);
            if (cy > 0)
            {
                bp[an] = cy;
            }
            return 0;
        }

        ReadOnlySpan<ulong> one = [1UL];
        if (an == bn)
        {
            int c = Compare(ap, bp, an);
            if (c == 0)
            {
                if (s > 0)
                {
                    hook.OnQuotient(one, swapped ? 1 : 0);
                }
                else
                {
                    hook.OnGcd(ap.Slice(0, an), isA: !swapped);
                }
                return 0;
            }
            hook.OnQuotient(one, swapped ? 1 : 0);
            if (c > 0)
            {
                Swap(ref ap, ref bp, ref swapped);
            }
        }
        else
        {
            hook.OnQuotient(one, swapped ? 1 : 0);
            if (an > bn)
            {
                Swap(ref ap, ref bp, ref swapped);
                (an, bn) = (bn, an);
            }
        }

        int qn = bn - an + 1;
        ulong[] pool = ArrayPool<ulong>.Shared.Rent(qn + an);
        Span<ulong> qp = pool.AsSpan(0, qn), rp = pool.AsSpan(qn, an);
        DivRem(qp, rp, bp, bn, ap, an);
        rp.CopyTo(bp);
        bn = TrimmedLength(bp, an);

        if (bn <= s)
        {
            if (s == 0)
            {
                hook.OnGcd(ap.Slice(0, an), isA: !swapped);
                ArrayPool<ulong>.Shared.Return(pool);
                return 0;
            }
            // The quotient is one too large: add back a and decrement it.
            if (bn > 0)
            {
                ulong cy = Add(bp, ap, an, bp, bn);
                if (cy != 0)
                {
                    bp[an++] = cy;
                }
            }
            else
            {
                ap.Slice(0, an).CopyTo(bp);
            }
            ulong b0 = SubScalarInPlace(qp, qn, 1);
            Debug.Assert(b0 == 0);
        }
        hook.OnQuotient(qp, swapped ? 1 : 0);
        ArrayPool<ulong>.Shared.Return(pool);
        return an;
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void Swap(ref Span<ulong> a, ref Span<ulong> b, ref bool swapped)
    {
        Span<ulong> t = a;
        a = b;
        b = t;
        swapped = !swapped;
    }

    // ---------------------------------------------------------------
    // hgcd: subquadratic half-gcd (Möller)
    // ---------------------------------------------------------------

    /// <summary>Extracts the normalized top two limbs of a and b at size n (n &gt;= 2).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void TopLimbs(ReadOnlySpan<ulong> a, ReadOnlySpan<ulong> b, int n, out ulong ah, out ulong al, out ulong bh, out ulong bl)
    {
        ulong mask = a[n - 1] | b[n - 1];
        Debug.Assert(mask != 0);
        ah = a[n - 1]; al = a[n - 2];
        bh = b[n - 1]; bl = b[n - 2];
        if ((long)mask >= 0)
        {
            int sh = BitOperations.LeadingZeroCount(mask);
            ulong a3 = n >= 3 ? a[n - 3] : 0, b3 = n >= 3 ? b[n - 3] : 0;
            ah = (ah << sh) | (al >> (64 - sh));
            al = (al << sh) | (a3 >> (64 - sh));
            bh = (bh << sh) | (bl >> (64 - sh));
            bl = (bl << sh) | (b3 >> (64 - sh));
        }
    }

    private static int HalfGcdStep(int n, Span<ulong> a, Span<ulong> b, int s, HalfGcdMatrix m)
    {
        ulong ah, al, bh, bl;
        if (n == s + 1)
        {
            if ((a[n - 1] | b[n - 1]) < 4)
            {
                return GcdSubdivideStep(a, b, n, s, m);
            }
            ah = a[n - 1]; al = a[n - 2];
            bh = b[n - 1]; bl = b[n - 2];
        }
        else
        {
            TopLimbs(a, b, n, out ah, out al, out bh, out bl);
        }
        if (HalfGcd2Limb(ah, al, bh, bl, out Matrix1 m1))
        {
            m.MulMatrix1(m1);
            return Matrix22Mul1InverseVector(m1, a, b, n);
        }
        return GcdSubdivideStep(a, b, n, s, m);
    }

    /// <summary>
    /// Reduces (a, b) of n limbs in place so that both stay above s = n/2 + 1 limbs, with
    /// the original pair equal to M times the reduced pair (M must start as identity).
    /// Returns the new size, or 0 when no reduction was possible. a, b need n+1 limbs.
    /// </summary>
    internal static int HalfGcd(Span<ulong> a, Span<ulong> b, int n, HalfGcdMatrix m)
    {
        int s = n / 2 + 1;
        if (n <= s)
        {
            return 0;
        }
        bool success = false;
        int nn;
        if (n >= HgcdThreshold)
        {
            int n2 = 3 * n / 4 + 1;
            int p = n / 2;
            nn = HalfGcd(a.Slice(p), b.Slice(p), n - p, m);
            if (nn > 0)
            {
                n = HalfGcdMatrixAdjust(m, p + nn, a, b, p);
                success = true;
            }
            while (n > n2)
            {
                nn = HalfGcdStep(n, a, b, s, m);
                if (nn == 0)
                {
                    return success ? n : 0;
                }
                n = nn;
                success = true;
            }
            if (n > s + 2)
            {
                p = 2 * s - n + 1;
                var m1 = new HalfGcdMatrix(n - p);
                nn = HalfGcd(a.Slice(p), b.Slice(p), n - p, m1);
                if (nn > 0)
                {
                    n = HalfGcdMatrixAdjust(m1, p + nn, a, b, p);
                    m.MulMatrix(m1);
                    success = true;
                }
            }
        }
        while (true)
        {
            nn = HalfGcdStep(n, a, b, s, m);
            if (nn == 0)
            {
                return success ? n : 0;
            }
            n = nn;
            success = true;
        }
    }

    /// <summary>
    /// Given a = α B^p + a_lo, b = β B^p + b_lo where (α, β) is the reduced top part (total
    /// size n), computes (a; b) &lt;- M^-1 (a; b). Returns the new size; storage n+1.
    /// </summary>
    internal static int HalfGcdMatrixAdjust(HalfGcdMatrix m, int n, Span<ulong> a, Span<ulong> b, int p)
    {
        m.Pad();
        int mn = m.N;
        Debug.Assert(p + mn < n);
        ulong[] pool = ArrayPool<ulong>.Shared.Rent(2 * (p + mn));
        Span<ulong> t0 = pool.AsSpan(0, p + mn), t1 = pool.AsSpan(p + mn, p + mn);
        ReadOnlySpan<ulong> p00 = m.R0.X.AsSpan(0, mn), p01 = m.R0.Y.AsSpan(0, mn);
        ReadOnlySpan<ulong> p10 = m.R1.X.AsSpan(0, mn), p11 = m.R1.Y.AsSpan(0, mn);

        // Values depending on a_lo first, before a is overwritten.
        MulAny(t0, p11, mn, a, p);
        MulAny(t1, p10, mn, a, p);

        // a = α B^p + p11 a_lo - p01 b_lo
        t0.Slice(0, p).CopyTo(a);
        ulong ah = Add(a.Slice(p), a.Slice(p), n - p, t0.Slice(p), mn);
        MulAny(t0, p01, mn, b, p);
        ulong cy = Sub(a, a, n, t0, p + mn);
        Debug.Assert(cy <= ah);
        ah -= cy;

        // b = β B^p + p00 b_lo - p10 a_lo
        MulAny(t0, p00, mn, b, p);
        t0.Slice(0, p).CopyTo(b);
        ulong bh = Add(b.Slice(p), b.Slice(p), n - p, t0.Slice(p), mn);
        cy = Sub(b, b, n, t1, p + mn);
        Debug.Assert(cy <= bh);
        bh -= cy;
        ArrayPool<ulong>.Shared.Return(pool);

        if (ah > 0 || bh > 0)
        {
            a[n] = ah;
            b[n] = bh;
            n++;
        }
        else if (a[n - 1] == 0 && b[n - 1] == 0)
        {
            n--;
        }
        Debug.Assert(a[n - 1] > 0 || b[n - 1] > 0);
        return n;
    }

    // ---------------------------------------------------------------
    // GCD drivers
    // ---------------------------------------------------------------

    /// <summary>
    /// Computes gcd(u, v) where both are stored in n limbs (at least one top limb nonzero,
    /// both values nonzero); u and v are destroyed and need n+1 limbs of storage. When row
    /// is given, it tracks the cofactor of the original u: returns the gcd and sets
    /// cofactor/negative so that g ≡ ±cofactor * u (mod v). With a row, may instead return
    /// null once both values fit one limb (left in u[0], v[0]); row then holds (p10, p11)
    /// with current u = p11 u0 - p01 v0 and current v = p00 v0 - p10 u0.
    /// </summary>
    internal static ulong[]? GcdCore(Span<ulong> u, Span<ulong> v, int n, LimbRow? row, out ulong[]? cofactor, out bool cofactorNegative)
    {
        var hook = new GcdCapture { Row = row };
        cofactor = null;
        cofactorNegative = false;
        while (n >= GcdDcThreshold)
        {
            int p = 2 * n / 3;
            var m = new HalfGcdMatrix(n - p);
            int nn = HalfGcd(u.Slice(p), v.Slice(p), n - p, m);
            if (nn > 0)
            {
                n = HalfGcdMatrixAdjust(m, p + nn, u, v, p);
                row?.MulMatrix(m);
            }
            else
            {
                n = GcdSubdivideStep(u, v, n, 0, hook);
                if (n == 0)
                {
                    goto Found;
                }
            }
        }
        while (n >= 2)
        {
            TopLimbs(u, v, n, out ulong uh, out ulong ul, out ulong vh, out ulong vl);
            if (HalfGcd2Limb(uh, ul, vh, vl, out Matrix1 m1))
            {
                n = Matrix22Mul1InverseVector(m1, u, v, n);
                row?.MulMatrix1(m1);
            }
            else
            {
                n = GcdSubdivideStep(u, v, n, 0, hook);
                if (n == 0)
                {
                    goto Found;
                }
            }
        }
        Debug.Assert(n == 1 && u[0] != 0 && v[0] != 0);
        if (row is null)
        {
            return [SmallPrimes.Gcd64(u[0], v[0])];
        }
        // Single-limb tail: the caller finishes with a word-level extended Euclid on u[0], v[0].
        return null;

    Found:
        if (row is not null)
        {
            // Original u0 = p00 A + p01 B, v0 = p10 A + p11 B (det 1), so
            // A = p11 u0 - p01 v0 and B = p00 v0 - p10 u0.
            cofactorNegative = !hook.IsA;
            ulong[] c = hook.IsA ? row.Y : row.X;
            cofactor = c.AsSpan(0, TrimmedLength(c, row.N)).ToArray();
        }
        return hook.G!;
    }
}
