using System.Buffers;
using System.Diagnostics;
using Sdcb.Arithmetic.Core;

namespace Sdcb.Arithmetic;

public readonly partial struct BigInt
{
    // ---------------------------------------------------------------
    // GCD / LCM (Lehmer's algorithm with 64-bit inner loop)
    // ---------------------------------------------------------------

    /// <summary>Greatest common divisor of |a| and |b| (mpz_gcd).</summary>
    public static BigInt Gcd(BigInt a, BigInt b)
    {
        a = a.Abs();
        b = b.Abs();
        if (a._size == 0)
        {
            return b;
        }
        if (b._size == 0)
        {
            return a;
        }
        if (CompareAbs(a, b) < 0)
        {
            (a, b) = (b, a);
        }
        if (b.LimbCount == 1)
        {
            ulong b0 = b._d![0];
            return new BigInt(SmallPrimes.Gcd64(b0, a.LimbCount == 1 ? a._d![0] : Nat.RemScalar(a.Limbs, a.LimbCount, b0)));
        }
        if (!PrepareGcdOperands(a, b, out ulong[] u, out ulong[] v, out int n, out BigInt early, out _))
        {
            return early;
        }
        ulong[] g = Nat.GcdCore(u, v, n, null, out _, out _)!;
        return FromMagnitude(g, g.Length, 1);
    }

    /// <summary>
    /// Brings |a| &gt; |b| &gt; 0 to equal-size working buffers (n limbs, n+1 storage) after one
    /// initial division when the sizes differ. Returns false with the result in early when
    /// that division already finishes (b divides a). The first buffer's cofactor equals the
    /// cofactor of a modulo b.
    /// </summary>
    private static bool PrepareGcdOperands(BigInt a, BigInt b, out ulong[] u, out ulong[] v, out int n, out BigInt early, out bool exact)
    {
        int an = a.LimbCount, bn = b.LimbCount;
        n = bn;
        u = new ulong[bn + 1];
        v = new ulong[bn + 1];
        b.Limbs.CopyTo(v);
        early = default;
        exact = false;
        if (an > bn)
        {
            ulong[] q = ArrayPool<ulong>.Shared.Rent(an - bn + 1);
            Nat.DivRem(q, u, a.Limbs, an, b.Limbs, bn);
            ArrayPool<ulong>.Shared.Return(q);
            if (Nat.IsZero(u, bn))
            {
                early = b;
                exact = true;
                return false;
            }
        }
        else
        {
            a.Limbs.CopyTo(u);
        }
        return true;
    }

    /// <summary>
    /// dst = p*A + q*B where A = a[0..an), B = b[0..bn), an &gt;= bn, |p|,|q| &lt; 2^62
    /// and the result is guaranteed non-negative. Returns the normalized length.
    /// dst must have an+1 limbs capacity.
    /// </summary>
    private static int LinearCombination(Span<ulong> dst, ReadOnlySpan<ulong> a, int an, ReadOnlySpan<ulong> b, int bn, long p, long q)
    {
        Debug.Assert(an >= bn);
        if (p == 0)
        {
            Debug.Assert(q >= 0);
            if (q == 0)
            {
                return 0;
            }
            dst[bn] = Nat.MulScalar(dst, b.Slice(0, bn), bn, (ulong)q);
            return Nat.TrimmedLength(dst, bn + 1);
        }
        if (q == 0)
        {
            Debug.Assert(p >= 0);
            dst[an] = Nat.MulScalar(dst, a.Slice(0, an), an, (ulong)p);
            return Nat.TrimmedLength(dst, an + 1);
        }
        if (p > 0 && q > 0)
        {
            dst[an] = Nat.MulScalar(dst, a.Slice(0, an), an, (ulong)p);
            ulong carry = Nat.AddMulScalar(dst.Slice(0, bn), b.Slice(0, bn), bn, (ulong)q);
            if (bn < an + 1)
            {
                Nat.AddScalarInPlace(dst.Slice(bn), an + 1 - bn, carry);
            }
            return Nat.TrimmedLength(dst, an + 1);
        }
        if (p > 0)
        {
            // p*A - |q|*B
            dst[an] = Nat.MulScalar(dst, a.Slice(0, an), an, (ulong)p);
            ulong borrow = Nat.SubMulScalar(dst.Slice(0, bn), b.Slice(0, bn), bn, (ulong)(-q));
            if (bn < an + 1)
            {
                borrow = Nat.SubScalarInPlace(dst.Slice(bn), an + 1 - bn, borrow);
            }
            Debug.Assert(borrow == 0);
            return Nat.TrimmedLength(dst, an + 1);
        }
        // q*B - |p|*A; non-negative implies an <= bn + 1 <= capacity.
        Debug.Assert(q > 0);
        dst[bn] = Nat.MulScalar(dst, b.Slice(0, bn), bn, (ulong)q);
        int len = bn + 1;
        Debug.Assert(an <= len);
        ulong borrow2 = Nat.SubMulScalar(dst.Slice(0, an), a.Slice(0, an), an, (ulong)(-p));
        if (an < len)
        {
            borrow2 = Nat.SubScalarInPlace(dst.Slice(an), len - an, borrow2);
        }
        Debug.Assert(borrow2 == 0);
        return Nat.TrimmedLength(dst, len);
    }

    /// <summary>
    /// Jacobi symbol bookkeeping along a Euclidean remainder sequence (u, v), u &gt;= v:
    /// the symbol equals Sign * (num/den) where den is odd and is u when DenFirst, else v.
    /// Lu/Lv hold u and v modulo 8, which is all reciprocity needs.
    /// </summary>
    private struct JacobiState
    {
        public int Sign;
        public bool DenFirst;
        public uint Lu, Lv;

        /// <summary>Applies the step u' = u - q*v and rotates the pair to (v, u').</summary>
        public void Step(ulong q)
        {
            uint lu2 = (uint)(Lu - q * Lv) & 7;
            if (DenFirst)
            {
                if ((Lv & 1) != 0)
                {
                    // (v/u) = ±(u/v) = ±(u'/v): denominator moves to v.
                    if ((Lu & 3) == 3 && (Lv & 3) == 3)
                    {
                        Sign = -Sign;
                    }
                }
                else
                {
                    // v even: (v/u') differs from (v/u) only through the (2/n) and
                    // reciprocity factors of v = 2w when v is exactly once even.
                    if ((Lv & 3) == 2)
                    {
                        uint w = (Lv >> 1) & 3;
                        if (TwoFactor(w, Lu) != TwoFactor(w, lu2))
                        {
                            Sign = -Sign;
                        }
                    }
                    DenFirst = false;
                    Lu = Lv;
                    Lv = lu2;
                    return;
                }
            }
            // Denominator is v, which becomes the first element of the rotated pair.
            DenFirst = true;
            Lu = Lv;
            Lv = lu2;
        }

        public void Swap()
        {
            DenFirst = !DenFirst;
            (Lu, Lv) = (Lv, Lu);
        }

        private static int TwoFactor(uint w, uint n)
        {
            n &= 7;
            int s = n is 3 or 5 ? -1 : 1;
            return (w & 3) == 3 && (n & 3) == 3 ? -s : s;
        }
    }

    /// <summary>
    /// Lehmer step on the top 62 bits of (u, v), u &gt;= v, na &gt;= 2: returns the cofactor
    /// matrix of the longest exactly-simulated run of Euclid quotients (identity when none).
    /// Signed 64-bit arithmetic stays overflow-free under the 2^61 cofactor cap.
    /// </summary>
    private static void LehmerMatrix(ReadOnlySpan<ulong> ua, int na, ReadOnlySpan<ulong> ub, int nb,
        out long ma, out long mb, out long mc, out long md, ref JacobiState js, bool trackJacobi)
    {
        long shift = Nat.BitLength(ua, na) - 62;
        Debug.Assert(shift > 0);
        long x = (long)ExtractBits64Span(ua, na, shift);
        long y = Nat.BitLength(ub, nb) <= shift ? 0 : (long)ExtractBits64Span(ub, nb, shift);
        Debug.Assert(x >= y && y >= 0);

        ma = 1; mb = 0; mc = 0; md = 1;
        const long Cap = 1L << 61;
        while (true)
        {
            long dy1 = y + mc, dy2 = y + md;
            if (dy1 <= 0 || dy2 <= 0)
            {
                break;
            }
            long q = (x + ma) / dy1;
            if (q != (x + mb) / dy2)
            {
                break;
            }
            Int128 nmc = (Int128)ma - (Int128)q * mc;
            Int128 nmd = (Int128)mb - (Int128)q * md;
            if (nmc > Cap || nmc < -Cap || nmd > Cap || nmd < -Cap)
            {
                break;
            }
            if (trackJacobi)
            {
                js.Step((ulong)q);
            }
            (ma, mc) = (mc, (long)nmc);
            (mb, md) = (md, (long)nmd);
            long t = x - q * y;
            x = y;
            y = t;
        }
    }

    /// <summary>Extracts 64 bits starting at bit index shift from a raw magnitude.</summary>
    private static ulong ExtractBits64Span(ReadOnlySpan<ulong> m, int n, long shift)
    {
        int limb = (int)(shift >> 6);
        int ofs = (int)(shift & 63);
        ulong lo = limb < n ? m[limb] : 0;
        ulong hi = limb + 1 < n ? m[limb + 1] : 0;
        return ofs == 0 ? lo : (lo >> ofs) | (hi << (64 - ofs));
    }

    /// <summary>GCD with a 64-bit value (mpz_gcd_ui).</summary>
    public static ulong Gcd(BigInt a, ulong b)
    {
        if (b == 0)
        {
            BigInt g = a.Abs();
            return g.LimbCount <= 1 ? (g._size == 0 ? 0 : g._d![0]) : throw new OverflowException("GCD does not fit in UInt64.");
        }
        if (a._size == 0)
        {
            return b;
        }
        return SmallPrimes.Gcd64(b, Nat.RemScalar(a.Limbs, a.LimbCount, b));
    }

    /// <summary>Least common multiple (mpz_lcm); result is non-negative.</summary>
    public static BigInt Lcm(BigInt a, BigInt b)
    {
        if (a._size == 0 || b._size == 0)
        {
            return default;
        }
        BigInt g = Gcd(a, b);
        return (a.Abs() / g) * b.Abs();
    }

    /// <summary>Multiplies by a small signed factor.</summary>
    private static BigInt MulSigned(BigInt v, long f)
    {
        if (f == 0 || v._size == 0)
        {
            return default;
        }
        ulong mag = f < 0 ? (ulong)(-f) : (ulong)f;
        int n = v.LimbCount;
        ulong[] d = new ulong[n + 1];
        d[n] = Nat.MulScalar(d, v.Limbs, n, mag);
        int sign = (f < 0) ^ (v._size < 0) ? -1 : 1;
        return FromMagnitude(d, n + 1, sign);
    }

    private static BigInt MulSigned128(BigInt v, Int128 f)
    {
        if (f == 0 || v._size == 0)
        {
            return default;
        }
        if (f == 1)
        {
            return v;
        }
        if (f >= long.MinValue + 1 && f <= long.MaxValue)
        {
            return MulSigned(v, (long)f);
        }
        return v * new BigInt(f);
    }

    // ---------------------------------------------------------------
    // Extended GCD / modular inverse
    // ---------------------------------------------------------------

    /// <summary>
    /// Extended GCD (mpz_gcdext): returns g = gcd(|a|, |b|) &gt;= 0 and cofactors
    /// with a*x + b*y == g.
    /// </summary>
    public static BigInt ExtendedGcd(BigInt a, BigInt b, out BigInt x, out BigInt y)
    {
        BigInt aa = a.Abs(), ba = b.Abs();
        if (ba._size == 0)
        {
            x = aa._size == 0 ? Zero : a._size < 0 ? MinusOne : One;
            y = Zero;
            return aa;
        }
        if (aa._size == 0)
        {
            x = Zero;
            y = b._size < 0 ? MinusOne : One;
            return ba;
        }

        int cmp = CompareAbs(aa, ba);
        if (cmp == 0)
        {
            x = Zero;
            y = b._size < 0 ? MinusOne : One;
            return aa;
        }

        // xa: cofactor of |a| modulo |b| (any representative).
        BigInt g, xa;
        if (aa.LimbCount == 1 && ba.LimbCount == 1)
        {
            ulong av = aa._d![0], bv = ba._d![0], p = av, r = bv;
            long m00 = 1, m10 = 0; // |cofactors| <= bv/g fit in 64 bits for these sizes
            while (r != 0)
            {
                ulong q = p / r;
                (p, r) = (r, p - q * r);
                (m00, m10) = (m10, m00 - (long)q * m10);
            }
            // Canonical range directly in 128-bit arithmetic.
            Int128 mod1 = bv / p, xs = m00 % mod1;
            if (xs < 0)
            {
                xs += mod1;
            }
            if (2 * xs > mod1)
            {
                xs -= mod1;
            }
            Int128 ys = ((Int128)p - (Int128)av * xs) / bv;
            x = new BigInt(a._size < 0 ? -xs : xs);
            y = new BigInt(b._size < 0 ? -ys : ys);
            return new BigInt(p);
        }
        else if (cmp > 0)
        {
            // gcd(|a| mod |b|, |b|): the cofactor of the remainder is that of |a|.
            if (!PrepareGcdOperands(aa, ba, out ulong[] u, out ulong[] v, out int n, out g, out _))
            {
                xa = Zero;
            }
            else
            {
                xa = GcdExtCore(u, v, n, out g);
            }
        }
        else
        {
            // Cofactor of |b| modulo |a| (from gcd(|b| mod |a|, |a|)), then that of |a|.
            if (!PrepareGcdOperands(ba, aa, out ulong[] u, out ulong[] v, out int n, out g, out _))
            {
                xa = One;
            }
            else
            {
                BigInt tb = GcdExtCore(u, v, n, out g);
                xa = DivExact(g - ba * tb, aa);
            }
        }

        // Canonical cofactors (mpz_gcdext): |x| < |b|/(2g), |y| < |a|/(2g), with x = sgn(a)
        // when |b| = 2g.
        BigInt mod = DivExact(ba, g);
        xa = Mod(xa, mod);
        if (xa.ShiftLeft(1) > mod)
        {
            xa -= mod;
        }
        BigInt yb = DivExact(g - aa * xa, ba);
        x = a._size < 0 ? -xa : xa;
        y = b._size < 0 ? -yb : yb;
        return g;
    }

    /// <summary>Runs the gcd core with cofactor tracking for the first buffer; returns that (signed) cofactor.</summary>
    private static BigInt GcdExtCore(ulong[] u, ulong[] v, int n, out BigInt g)
    {
        var row = new LimbRow(n + 2, identityX: false);
        ulong[]? gd = Nat.GcdCore(u, v, n, row, out ulong[]? c, out bool neg);
        if (gd is not null)
        {
            g = FromMagnitude(gd, gd.Length, 1);
            return c is null || c.Length == 0 ? Zero : FromMagnitude(c, c.Length, neg ? -1 : 1);
        }
        // Word-level extended Euclid on the single-limb pair: s*U + t*V = g, combined with
        // U = p11 u0 - p01 v0, V = p00 v0 - p10 u0, so the cofactor of u0 is s p11 - t p10.
        ulong p = u[0], r = v[0];
        Int128 s0 = 1, s1 = 0, t0 = 0, t1 = 1;
        while (r != 0)
        {
            ulong q = p / r;
            (p, r) = (r, p - q * r);
            (s0, s1) = (s1, s0 - q * s1);
            (t0, t1) = (t1, t0 - q * t1);
        }
        g = new BigInt(p);
        BigInt p10 = FromMagnitude(row.X, row.N, 1), p11 = FromMagnitude(row.Y, row.N, 1);
        return MulSigned128(p11, s0) - MulSigned128(p10, t0);
    }

    /// <summary>Modular inverse (mpz_invert): returns x with a*x ≡ 1 (mod m), 0 &lt;= x &lt; |m|. Throws when not invertible.</summary>
    public static BigInt ModInverse(BigInt a, BigInt m)
    {
        if (!TryModInverse(a, m, out BigInt inv))
        {
            throw new ArithmeticException("Value is not invertible modulo m.");
        }
        return inv;
    }

    /// <summary>Tries to compute the modular inverse; false when gcd(a, m) != 1.</summary>
    public static bool TryModInverse(BigInt a, BigInt m, out BigInt inverse)
    {
        inverse = default;
        BigInt mAbs = m.Abs();
        if (mAbs._size == 0)
        {
            return false;
        }
        if (mAbs.IsOne)
        {
            inverse = Zero;
            return true;
        }
        BigInt g = ExtendedGcd(a, mAbs, out BigInt x, out _);
        if (!g.IsOne)
        {
            return false;
        }
        inverse = Mod(x, mAbs);
        return true;
    }

    // ---------------------------------------------------------------
    // Jacobi / Legendre / Kronecker symbols
    // ---------------------------------------------------------------

    /// <summary>Jacobi symbol (a/n); n must be odd and positive (mpz_jacobi).</summary>
    public static int Jacobi(BigInt a, BigInt n)
    {
        if (n._size <= 0 || n.IsEven)
        {
            throw new ArgumentException("Jacobi symbol requires an odd positive denominator.", nameof(n));
        }
        return JacobiCore(Mod(a, n), n);
    }

    /// <summary>Legendre symbol (a/p) for odd prime p (mpz_legendre).</summary>
    public static int Legendre(BigInt a, BigInt p) => Jacobi(a, p);

    /// <summary>Kronecker symbol (a/n), the Jacobi symbol extended to all integers (mpz_kronecker).</summary>
    public static int Kronecker(BigInt a, BigInt n)
    {
        if (n._size == 0)
        {
            return a.IsOne || a.IsMinusOne ? 1 : 0;
        }
        int result = 1;
        if (n._size < 0)
        {
            n = -n;
            if (a._size < 0)
            {
                result = -result;
            }
        }
        if (n.IsEven)
        {
            if (a.IsEven)
            {
                return 0;
            }
            long twos = n.TrailingZeroCount();
            n = n.ShiftRightTruncate(twos);
            if ((twos & 1) != 0)
            {
                // (2/a) = 1 if a ≡ ±1 (mod 8), -1 if a ≡ ±3 (mod 8)
                ulong am8 = a._size < 0 ? (ulong)((8 - (long)RemainderAbs(a, 8)) & 7) : RemainderAbs(a, 8);
                if (am8 == 3 || am8 == 5)
                {
                    result = -result;
                }
            }
        }
        if (n.IsOne)
        {
            return result;
        }
        return result * JacobiCore(Mod(a, n), n);
    }

    /// <summary>Jacobi symbol for 0 &lt;= a &lt; n, n odd &gt;= 1, via Lehmer-accelerated Euclid with mod-8 reciprocity tracking.</summary>
    private static int JacobiCore(BigInt a, BigInt n)
    {
        if (a._size == 0)
        {
            return n.IsOne ? 1 : 0;
        }
        if (n.LimbCount == 1)
        {
            return Jacobi64(a._d![0], n._d![0]);
        }
        int cap = n.LimbCount + 1;
        ulong[] ua = new ulong[cap], ub = new ulong[cap], t1 = new ulong[cap], t2 = new ulong[cap];
        n.Limbs.CopyTo(ua);
        a.Limbs.CopyTo(ub);
        int na = n.LimbCount, nb = a.LimbCount;
        var js = new JacobiState { Sign = 1, DenFirst = true, Lu = (uint)ua[0] & 7, Lv = (uint)ub[0] & 7 };

        while (nb > 0)
        {
            if (na == 1)
            {
                ulong p = ua[0], r = ub[0];
                while (r != 0)
                {
                    ulong q = p / r;
                    js.Step(q);
                    (p, r) = (r, p - q * r);
                }
                return p == 1 ? js.Sign : 0;
            }
            long ma = 1, mb = 0, mc = 0, md = 1;
            JacobiState saved = js;
            if (nb > 1)
            {
                LehmerMatrix(ua, na, ub, nb, out ma, out mb, out mc, out md, ref js, trackJacobi: true);
            }
            if (mb == 0)
            {
                js = saved;
                Nat.DivRem(t2, t1, ua.AsSpan(0, na), na, ub.AsSpan(0, nb), nb);
                int nr = Nat.TrimmedLength(t1, nb);
                // Only q mod 8 matters for the low-bit bookkeeping.
                js.Step(t2[0]);
                (ua, ub, t1) = (ub, t1, ua);
                na = nb;
                nb = nr;
            }
            else
            {
                int n1 = LinearCombination(t1, ua, na, ub, nb, ma, mb);
                int n2 = LinearCombination(t2, ua, na, ub, nb, mc, md);
                (ua, t1) = (t1, ua);
                (ub, t2) = (t2, ub);
                na = n1;
                nb = n2;
                if (Nat.CompareNormalized(ua.AsSpan(0, na), na, ub.AsSpan(0, nb), nb) < 0)
                {
                    (ua, ub) = (ub, ua);
                    (na, nb) = (nb, na);
                    js.Swap();
                }
            }
        }
        return na == 1 && ua[0] == 1 ? js.Sign : 0;
    }

    /// <summary>Jacobi symbol on 64-bit values: 0 &lt;= a, n odd &gt;= 1.</summary>
    internal static int Jacobi64(ulong a, ulong n)
    {
        Debug.Assert((n & 1) != 0);
        int result = 1;
        a %= n;
        while (a != 0)
        {
            int twos = System.Numerics.BitOperations.TrailingZeroCount(a);
            a >>= twos;
            if ((twos & 1) != 0)
            {
                ulong nm8 = n & 7;
                if (nm8 == 3 || nm8 == 5)
                {
                    result = -result;
                }
            }
            if ((a & 3) == 3 && (n & 3) == 3)
            {
                result = -result;
            }
            (a, n) = (n % a, a);
        }
        return n == 1 ? result : 0;
    }
}
