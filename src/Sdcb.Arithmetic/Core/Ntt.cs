using System.Diagnostics;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.Arm;
using System.Runtime.Intrinsics.X86;

namespace Sdcb.Arithmetic.Core;

/// <summary>
/// Large-operand multiplication by number-theoretic transforms over three ~30-bit
/// primes, recombined with Garner's CRT.
/// Operands are cut into 32-bit coefficients; transforms use Harvey's lazy butterflies
/// with Shoup precomputed twiddles and run on Vector512 → Vector256 → Vector128 →
/// scalar, following <see cref="Isa"/>. Forward transforms are DIF (natural in,
/// bit-reversed out) and the inverse is DIT (bit-reversed in, natural out), so no
/// permutation pass is needed.
/// </summary>
internal static class Ntt
{
    /// <summary>Largest supported transform, in 32-bit coefficients (bounds twiddle-table memory).</summary>
    public const int MaxLog = 20;

    private static readonly Prime[] Primes = [new(998244353), new(985661441), new(754974721)];

    /// <summary>True when an (an x bn)-limb product fits a single transform.</summary>
    public static bool Fits(int an, int bn) => BitConverter.IsLittleEndian && 2L * (an + bn) <= 1L << MaxLog;

    /// <summary>r[0..an+bn) = a * b; b may alias a (squaring runs one forward transform).</summary>
    public static void Mul(Span<ulong> r, ReadOnlySpan<ulong> a, int an, ReadOnlySpan<ulong> b, int bn)
    {
        Debug.Assert(Fits(an, bn) && an >= 1 && bn >= 1);
        bool square = an == bn && Unsafe.AreSame(ref MemoryMarshal.GetReference(a), ref MemoryMarshal.GetReference(b));
        ReadOnlySpan<uint> a32 = MemoryMarshal.Cast<ulong, uint>(a.Slice(0, an));
        ReadOnlySpan<uint> b32 = MemoryMarshal.Cast<ulong, uint>(b.Slice(0, bn));
        int convLen = a32.Length + b32.Length - 1;
        // A cyclic transform of size n >= max(len) folds c[n+k] onto c[k]; up to WrapMax
        // wrapped coefficients are computed directly and subtracted, so a product just past
        // a power of two does not double the transform size.
        int log = 32 - BitOperations.LeadingZeroCount((uint)Math.Max(convLen - WrapMax, Math.Max(a32.Length, b32.Length)) - 1);
        if (log < 1)
        {
            log = 1;
        }
        int n = 1 << log;
        int wrap = Math.Max(0, convLen - n);
        UInt128[]? high = wrap > 0 ? HighCoefficients(a32, b32, n, convLen) : null;

        uint[] buf = System.Buffers.ArrayPool<uint>.Shared.Rent(4 * n);
        try
        {
            Span<uint> res = buf.AsSpan(0, 3 * n);
            Span<uint> fb = buf.AsSpan(3 * n, n);
            for (int k = 0; k < 3; k++)
            {
                Prime pr = Primes[k];
                Tables t = pr.GetTables(log);
                Span<uint> fa = res.Slice(k * n, n);
                Load(fa, a32, pr.P);
                Forward(fa, log, pr, t);
                if (square)
                {
                    Pointwise(fa, fa, pr);
                }
                else
                {
                    Load(fb, b32, pr.P);
                    Forward(fb, log, pr, t);
                    Pointwise(fa, fb, pr);
                }
                Inverse(fa, log, pr, t);
            }
            Recombine(MemoryMarshal.Cast<ulong, uint>(r.Slice(0, an + bn)), res, n, convLen, high);
        }
        finally
        {
            System.Buffers.ArrayPool<uint>.Shared.Return(buf);
        }
    }

    // ---------------------------------------------------------------
    // Prime field data
    // ---------------------------------------------------------------

    private sealed class Tables
    {
        // Layout: element [m + j] holds ω_{2m}^{±j} for m = 1..n/2, j < m, with Shoup companions.
        public readonly uint[] Fw, FwS, Iw, IwS;
        public readonly int Log;

        public Tables(Prime pr, int log)
        {
            Log = log;
            int n = 1 << log;
            Fw = new uint[n]; FwS = new uint[n]; Iw = new uint[n]; IwS = new uint[n];
            for (int m = 1; m < n; m <<= 1)
            {
                uint w = pr.RootOfUnity(2 * m);
                uint wi = Prime.PowMod(w, pr.P - 2, pr.P);
                ulong x = 1, y = 1;
                for (int j = 0; j < m; j++)
                {
                    Fw[m + j] = (uint)x;
                    FwS[m + j] = pr.Shoup((uint)x);
                    Iw[m + j] = (uint)y;
                    IwS[m + j] = pr.Shoup((uint)y);
                    x = x * w % pr.P;
                    y = y * wi % pr.P;
                }
            }
        }
    }

    private sealed class Prime
    {
        public readonly uint P;
        public readonly uint PInvNeg;   // -P^-1 mod 2^32 (Montgomery)
        public readonly uint Generator; // quadratic non-residue: generates the full 2-power subgroup
        private Tables? _tables;
        private readonly object _lock = new();

        public Prime(uint p)
        {
            P = p;
            uint inv = p;
            for (int i = 0; i < 5; i++)
            {
                inv *= 2 - p * inv;
            }
            PInvNeg = 0 - inv;
            uint g = 2;
            while (PowMod(g, (p - 1) / 2, p) != p - 1)
            {
                g++;
            }
            Generator = g;
        }

        public uint RootOfUnity(int order)
        {
            int v = BitOperations.TrailingZeroCount(P - 1);
            Debug.Assert(order <= 1 << v);
            uint w = PowMod(Generator, (P - 1) >> v, P); // order exactly 2^v
            for (int o = 1 << v; o > order; o >>= 1)
            {
                w = (uint)((ulong)w * w % P);
            }
            return w;
        }

        public uint Shoup(uint w) => (uint)(((ulong)w << 32) / P);

        public Tables GetTables(int log)
        {
            Tables? t = Volatile.Read(ref _tables);
            if (t is not null && t.Log >= log)
            {
                return t;
            }
            lock (_lock)
            {
                t = _tables;
                if (t is null || t.Log < log)
                {
                    t = new Tables(this, Math.Max(log, 10));
                    Volatile.Write(ref _tables, t);
                }
                return t;
            }
        }

        public static uint PowMod(uint b, uint e, uint p)
        {
            ulong r = 1, x = b % p;
            while (e != 0)
            {
                if ((e & 1) != 0)
                {
                    r = r * x % p;
                }
                x = x * x % p;
                e >>= 1;
            }
            return (uint)r;
        }
    }

    // ---------------------------------------------------------------
    // Scalar modular primitives (values kept lazily in [0, 2p) / [0, 4p))
    // ---------------------------------------------------------------

    /// <summary>x*w mod p in [0, 2p) for any x &lt; 2^32 (Shoup).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static uint ShoupMul(uint x, uint w, uint ws, uint p)
    {
        uint q = (uint)(((ulong)x * ws) >> 32);
        return x * w - q * p;
    }

    /// <summary>a*b*2^-32 mod p in [0, 2p) for a*b &lt; 2^32 * p (Montgomery).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static uint MontMul(uint a, uint b, uint p, uint pinvNeg)
    {
        ulong t = (ulong)a * b;
        uint m = (uint)t * pinvNeg;
        return (uint)((t + (ulong)m * p) >> 32);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static uint Reduce2p(uint x, uint p2) => Math.Min(x, x - p2);

    /// <summary>Copies coefficients reduced into [0, 2p), zero-padding the rest.</summary>
    private static void Load(Span<uint> dst, ReadOnlySpan<uint> src, uint p)
    {
        uint p2 = 2 * p;
        int i = 0;
        ref uint d = ref MemoryMarshal.GetReference(dst);
        ref uint s = ref MemoryMarshal.GetReference(src);
        if (Isa.Hw512)
        {
            Vector512<uint> vp2 = Vector512.Create(p2);
            for (; i + 16 <= src.Length; i += 16)
            {
                Vector512<uint> x = Vector512.LoadUnsafe(ref s, (nuint)i);
                x = Vector512.Min(x, x - vp2);
                Vector512.Min(x, x - vp2).StoreUnsafe(ref d, (nuint)i);
            }
        }
        if (Isa.Hw256)
        {
            Vector256<uint> vp2 = Vector256.Create(p2);
            for (; i + 8 <= src.Length; i += 8)
            {
                Vector256<uint> x = Vector256.LoadUnsafe(ref s, (nuint)i);
                x = Vector256.Min(x, x - vp2);
                Vector256.Min(x, x - vp2).StoreUnsafe(ref d, (nuint)i);
            }
        }
        if (Isa.Hw128)
        {
            Vector128<uint> vp2 = Vector128.Create(p2);
            for (; i + 4 <= src.Length; i += 4)
            {
                Vector128<uint> x = Vector128.LoadUnsafe(ref s, (nuint)i);
                x = Vector128.Min(x, x - vp2);
                Vector128.Min(x, x - vp2).StoreUnsafe(ref d, (nuint)i);
            }
        }
        for (; i < src.Length; i++)
        {
            Unsafe.Add(ref d, i) = Reduce2p(Reduce2p(Unsafe.Add(ref s, i), p2), p2);
        }
        dst.Slice(src.Length).Clear();
    }

    // ---------------------------------------------------------------
    // Transform drivers: pick the widest enabled vector tier.
    // ---------------------------------------------------------------

    private static void Forward(Span<uint> a, int log, Prime pr, Tables t)
    {
        if (Isa.Hw512 && log >= 5)
        {
            ForwardVec<Vector512<uint>, V512>(a, log, pr, t);
        }
        else if (Isa.Hw256 && log >= 4)
        {
            ForwardVec<Vector256<uint>, V256>(a, log, pr, t);
        }
        else if (Isa.Hw128Mul && log >= 3)
        {
            ForwardVec<Vector128<uint>, V128>(a, log, pr, t);
        }
        else
        {
            ForwardScalar(a, log, pr, t, 1 << log);
        }
    }

    private static void Inverse(Span<uint> a, int log, Prime pr, Tables t)
    {
        if (Isa.Hw512 && log >= 5)
        {
            InverseVec<Vector512<uint>, V512>(a, log, pr, t);
        }
        else if (Isa.Hw256 && log >= 4)
        {
            InverseVec<Vector256<uint>, V256>(a, log, pr, t);
        }
        else if (Isa.Hw128Mul && log >= 3)
        {
            InverseVec<Vector128<uint>, V128>(a, log, pr, t);
        }
        else
        {
            InverseScalar(a, log, pr, t, 1);
        }
    }

    private static void Pointwise(Span<uint> a, ReadOnlySpan<uint> b, Prime pr)
    {
        int i = 0, n = a.Length;
        ref uint ap = ref MemoryMarshal.GetReference(a);
        ref uint bp = ref MemoryMarshal.GetReference(b);
        if (Isa.Hw512)
        {
            i = PointwiseVec<Vector512<uint>, V512>(ref ap, ref bp, i, n, pr);
        }
        if (Isa.Hw256)
        {
            i = PointwiseVec<Vector256<uint>, V256>(ref ap, ref bp, i, n, pr);
        }
        if (Isa.Hw128Mul)
        {
            i = PointwiseVec<Vector128<uint>, V128>(ref ap, ref bp, i, n, pr);
        }
        for (; i < n; i++)
        {
            Unsafe.Add(ref ap, i) = MontMul(Unsafe.Add(ref ap, i), Unsafe.Add(ref bp, i), pr.P, pr.PInvNeg);
        }
    }

    // ---------------------------------------------------------------
    // Scalar transforms (also used for the smallest sizes).
    // ---------------------------------------------------------------

    /// <summary>DIF stages with half-size m from n/2 down to mStop (inclusive); inputs/outputs in [0, 2p).</summary>
    private static void ForwardScalar(Span<uint> a, int log, Prime pr, Tables t, int mTop)
    {
        uint p = pr.P, p2 = 2 * p;
        int n = 1 << log;
        ref uint d = ref MemoryMarshal.GetReference(a);
        for (int m = Math.Min(mTop, n) >> 1; m >= 1; m >>= 1)
        {
            for (int s = 0; s < n; s += 2 * m)
            {
                for (int j = 0; j < m; j++)
                {
                    ref uint xr = ref Unsafe.Add(ref d, s + j);
                    ref uint yr = ref Unsafe.Add(ref d, s + j + m);
                    uint x = xr, y = yr;
                    xr = Reduce2p(x + y, p2);
                    yr = ShoupMul(x - y + p2, t.Fw[m + j], t.FwS[m + j], p);
                }
            }
        }
    }

    /// <summary>DIT stages with half-size m from mStart up to n/2; inputs/outputs in [0, 4p).</summary>
    private static void InverseScalar(Span<uint> a, int log, Prime pr, Tables t, int mStart)
    {
        uint p = pr.P, p2 = 2 * p;
        int n = 1 << log;
        ref uint d = ref MemoryMarshal.GetReference(a);
        for (int m = mStart; m < n; m <<= 1)
        {
            for (int s = 0; s < n; s += 2 * m)
            {
                for (int j = 0; j < m; j++)
                {
                    ref uint xr = ref Unsafe.Add(ref d, s + j);
                    ref uint yr = ref Unsafe.Add(ref d, s + j + m);
                    uint x = Reduce2p(xr, p2);
                    uint tt = ShoupMul(yr, t.Iw[m + j], t.IwS[m + j], p);
                    xr = x + tt;
                    yr = x - tt + p2;
                }
            }
        }
    }

    // ---------------------------------------------------------------
    // Vector transforms, generic over the lane width.
    // ---------------------------------------------------------------

    private interface IVec<TV> where TV : struct
    {
        static abstract int W { get; }
        static abstract TV Load(ref uint p, nuint i);
        static abstract void Store(ref uint p, nuint i, TV v);
        static abstract TV Bc(uint x);
        static abstract TV Create(ReadOnlySpan<uint> s);
        static abstract TV Add(TV a, TV b);
        static abstract TV Sub(TV a, TV b);
        static abstract TV Min(TV a, TV b);
        static abstract TV MulLo(TV a, TV b);
        static abstract TV MulHi(TV a, TV b);
        static abstract TV MontMul(TV a, TV b, TV p, TV pinvNeg);
        /// <summary>Rearranges two consecutive vectors so that butterfly partners at distance m (&lt; W) share a lane; lane t of X has twiddle index t mod m.</summary>
        static abstract void Split(ref TV a, ref TV b, int m);
        static abstract void Merge(ref TV x, ref TV y, int m);
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static TV ShoupMulV<TV, TO>(TV x, TV w, TV ws, TV p) where TV : struct where TO : struct, IVec<TV>
        => TO.Sub(TO.MulLo(x, w), TO.MulLo(TO.MulHi(x, ws), p));

    private static void ForwardVec<TV, TO>(Span<uint> a, int log, Prime pr, Tables t) where TV : struct where TO : struct, IVec<TV>
    {
        int n = 1 << log, W = TO.W;
        TV p = TO.Bc(pr.P), p2 = TO.Bc(2 * pr.P);
        ref uint d = ref MemoryMarshal.GetReference(a);
        ref uint fw = ref MemoryMarshal.GetArrayDataReference(t.Fw);
        ref uint fws = ref MemoryMarshal.GetArrayDataReference(t.FwS);

        // Stages with m >= W: partners live in different vectors.
        for (int m = n >> 1; m >= W; m >>= 1)
        {
            for (int s = 0; s < n; s += 2 * m)
            {
                ref uint blk = ref Unsafe.Add(ref d, s);
                for (int j = 0; j < m; j += W)
                {
                    TV x = TO.Load(ref blk, (nuint)j), y = TO.Load(ref blk, (nuint)(j + m));
                    TV w = TO.Load(ref fw, (nuint)(m + j)), ws = TO.Load(ref fws, (nuint)(m + j));
                    TV sum = TO.Add(x, y);
                    TO.Store(ref blk, (nuint)j, TO.Min(sum, TO.Sub(sum, p2)));
                    TO.Store(ref blk, (nuint)(j + m), ShoupMulV<TV, TO>(TO.Add(TO.Sub(x, y), p2), w, ws, p));
                }
            }
        }

        // Stages with m < W: two vectors at a time, rearranged in registers.
        TV[] twv = new TV[2 * 8];
        int stages = 0;
        Span<uint> tmp = stackalloc uint[64];
        for (int m = W >> 1; m >= 1; m >>= 1, stages++)
        {
            for (int l = 0; l < W; l++)
            {
                tmp[l] = t.Fw[m + (l % m)];
                tmp[W + l] = t.FwS[m + (l % m)];
            }
            twv[2 * stages] = TO.Create(tmp.Slice(0, W));
            twv[2 * stages + 1] = TO.Create(tmp.Slice(W, W));
        }
        for (int s = 0; s < n; s += 2 * W)
        {
            TV va = TO.Load(ref d, (nuint)s), vb = TO.Load(ref d, (nuint)(s + W));
            int st = 0;
            for (int m = W >> 1; m >= 1; m >>= 1, st++)
            {
                TO.Split(ref va, ref vb, m);
                TV sum = TO.Add(va, vb);
                TV diff = TO.Add(TO.Sub(va, vb), p2);
                va = TO.Min(sum, TO.Sub(sum, p2));
                vb = m == 1 ? TO.Min(diff, TO.Sub(diff, p2)) : ShoupMulV<TV, TO>(diff, twv[2 * st], twv[2 * st + 1], p);
                TO.Merge(ref va, ref vb, m);
            }
            TO.Store(ref d, (nuint)s, va);
            TO.Store(ref d, (nuint)(s + W), vb);
        }
    }

    private static void InverseVec<TV, TO>(Span<uint> a, int log, Prime pr, Tables t) where TV : struct where TO : struct, IVec<TV>
    {
        int n = 1 << log, W = TO.W;
        TV p = TO.Bc(pr.P), p2 = TO.Bc(2 * pr.P);
        ref uint d = ref MemoryMarshal.GetReference(a);

        TV[] twv = new TV[2 * 8];
        Span<uint> tmp = stackalloc uint[64];
        int stages = 0;
        for (int m = 1; m < W; m <<= 1, stages++)
        {
            for (int l = 0; l < W; l++)
            {
                tmp[l] = t.Iw[m + (l % m)];
                tmp[W + l] = t.IwS[m + (l % m)];
            }
            twv[2 * stages] = TO.Create(tmp.Slice(0, W));
            twv[2 * stages + 1] = TO.Create(tmp.Slice(W, W));
        }
        for (int s = 0; s < n; s += 2 * W)
        {
            TV va = TO.Load(ref d, (nuint)s), vb = TO.Load(ref d, (nuint)(s + W));
            int st = 0;
            for (int m = 1; m < W; m <<= 1, st++)
            {
                TO.Split(ref va, ref vb, m);
                TV x = TO.Min(va, TO.Sub(va, p2));
                TV tt = m == 1 ? TO.Min(vb, TO.Sub(vb, p2)) : ShoupMulV<TV, TO>(vb, twv[2 * st], twv[2 * st + 1], p);
                va = TO.Add(x, tt);
                vb = TO.Add(TO.Sub(x, tt), p2);
                TO.Merge(ref va, ref vb, m);
            }
            TO.Store(ref d, (nuint)s, va);
            TO.Store(ref d, (nuint)(s + W), vb);
        }

        ref uint iw = ref MemoryMarshal.GetArrayDataReference(t.Iw);
        ref uint iws = ref MemoryMarshal.GetArrayDataReference(t.IwS);
        for (int m = W; m < n; m <<= 1)
        {
            for (int s = 0; s < n; s += 2 * m)
            {
                ref uint blk = ref Unsafe.Add(ref d, s);
                for (int j = 0; j < m; j += W)
                {
                    TV x = TO.Load(ref blk, (nuint)j), y = TO.Load(ref blk, (nuint)(j + m));
                    TV w = TO.Load(ref iw, (nuint)(m + j)), ws = TO.Load(ref iws, (nuint)(m + j));
                    x = TO.Min(x, TO.Sub(x, p2));
                    TV tt = ShoupMulV<TV, TO>(y, w, ws, p);
                    TO.Store(ref blk, (nuint)j, TO.Add(x, tt));
                    TO.Store(ref blk, (nuint)(j + m), TO.Add(TO.Sub(x, tt), p2));
                }
            }
        }
    }

    private static int PointwiseVec<TV, TO>(ref uint a, ref uint b, int i, int n, Prime pr) where TV : struct where TO : struct, IVec<TV>
    {
        TV p = TO.Bc(pr.P), pinv = TO.Bc(pr.PInvNeg);
        int W = TO.W;
        for (; i + W <= n; i += W)
        {
            TO.Store(ref a, (nuint)i, TO.MontMul(TO.Load(ref a, (nuint)i), TO.Load(ref b, (nuint)i), p, pinv));
        }
        return i;
    }

    // ---------------------------------------------------------------
    // CRT recombination (Garner) and carry propagation into 32-bit words.
    // ---------------------------------------------------------------

    private static void Recombine(Span<uint> r32, Span<uint> res, int n, int convLen, UInt128[]? high)
    {
        Prime p0 = Primes[0], p1 = Primes[1], p2 = Primes[2];
        uint P0 = p0.P, P1 = p1.P, P2 = p2.P;
        // Undo the Montgomery factor 2^-32 and the transform scale n: multiply by 2^32 / n.
        uint c0 = ScaleConst(p0, n), c1 = ScaleConst(p1, n), c2 = ScaleConst(p2, n);
        uint c0s = p0.Shoup(c0), c1s = p1.Shoup(c1), c2s = p2.Shoup(c2);
        uint inv01 = Prime.PowMod(P0 % P1, P1 - 2, P1), inv01s = p1.Shoup(inv01);
        uint p0m2 = P0 % P2, p0m2s = p2.Shoup(p0m2);
        uint inv012 = Prime.PowMod((uint)((ulong)P0 * P1 % P2), P2 - 2, P2), inv012s = p2.Shoup(inv012);
        ulong P01 = (ulong)P0 * P1;

        Span<uint> r0 = res.Slice(0, n), r1 = res.Slice(n, n), r2 = res.Slice(2 * n, n);
        var g = new GarnerConsts(c0, c0s, c1, c1s, c2, c2s, inv01, inv01s, p0m2, p0m2s, inv012, inv012s);
        int done = 0;
        ref uint g0 = ref MemoryMarshal.GetReference(r0);
        ref uint g1 = ref MemoryMarshal.GetReference(r1);
        ref uint g2 = ref MemoryMarshal.GetReference(r2);
        if (Isa.Hw512)
        {
            done = GarnerVec<Vector512<uint>, V512>(ref g0, ref g1, ref g2, done, Math.Min(convLen, n), g);
        }
        if (Isa.Hw256)
        {
            done = GarnerVec<Vector256<uint>, V256>(ref g0, ref g1, ref g2, done, Math.Min(convLen, n), g);
        }
        if (Isa.Hw128Mul)
        {
            done = GarnerVec<Vector128<uint>, V128>(ref g0, ref g1, ref g2, done, Math.Min(convLen, n), g);
        }
        for (int k = done; k < Math.Min(convLen, n); k++)
        {
            uint x0 = ShoupMul(r0[k], c0, c0s, P0);
            x0 = x0 >= P0 ? x0 - P0 : x0;
            uint x1 = ShoupMul(r1[k], c1, c1s, P1);
            x1 = x1 >= P1 ? x1 - P1 : x1;
            uint x2 = ShoupMul(r2[k], c2, c2s, P2);
            x2 = x2 >= P2 ? x2 - P2 : x2;
            // v1 = (x1 - x0) / p0 mod p1   (x0 < p0 < 2 p1)
            uint v1 = ShoupMul(x1 + 2 * P1 - x0, inv01, inv01s, P1);
            v1 = v1 >= P1 ? v1 - P1 : v1;
            // v2 = (x2 - x0 - v1 p0) / (p0 p1) mod p2
            uint a0 = x0 >= P2 ? x0 - P2 : x0;
            uint b0 = ShoupMul(v1, p0m2, p0m2s, P2);
            uint v2 = ShoupMul(x2 + 4 * P2 - a0 - b0, inv012, inv012s, P2);
            v2 = v2 >= P2 ? v2 - P2 : v2;
            r0[k] = x0;
            r1[k] = v1;
            r2[k] = v2;
        }

        // value_k = x0 + p0*v1 + p0*p1*v2, accumulated into 32-bit output words.
        ulong cl = 0, ch = 0;
        int outLen = r32.Length;
        for (int k = 0; k < outLen; k++)
        {
            if (k < convLen)
            {
                ulong lo, hi;
                if (k < n)
                {
                    lo = g0 + (ulong)g1 * P0;
                    hi = Nat.BigMul(P01, g2, out ulong mid);
                    ulong s = mid + lo;
                    hi += s < lo ? 1UL : 0UL;
                    lo = s;
                    if (high is not null && k < high.Length)
                    {
                        // Remove the wrapped-around coefficient c[k + n].
                        UInt128 v = new UInt128(hi, lo) - high[k];
                        hi = (ulong)(v >> 64);
                        lo = (ulong)v;
                    }
                    g0 = ref Unsafe.Add(ref g0, 1);
                    g1 = ref Unsafe.Add(ref g1, 1);
                    g2 = ref Unsafe.Add(ref g2, 1);
                }
                else
                {
                    UInt128 v = high![k - n];
                    hi = (ulong)(v >> 64);
                    lo = (ulong)v;
                }
                cl += lo;
                ch += hi + (cl < lo ? 1UL : 0UL);
            }
            r32[k] = (uint)cl;
            cl = (cl >> 32) | (ch << 32);
            ch >>= 32;
        }
        Debug.Assert(cl == 0 && ch == 0);
    }

    /// <summary>Largest number of wrapped coefficients corrected directly.</summary>
    private const int WrapMax = 96;

    /// <summary>Exact convolution coefficients c[n..convLen) (few terms each).</summary>
    private static UInt128[] HighCoefficients(ReadOnlySpan<uint> a, ReadOnlySpan<uint> b, int n, int convLen)
    {
        var high = new UInt128[convLen - n];
        for (int k = n; k < convLen; k++)
        {
            UInt128 s = 0;
            int iLo = Math.Max(0, k - b.Length + 1), iHi = Math.Min(k, a.Length - 1);
            for (int i = iLo; i <= iHi; i++)
            {
                s += (ulong)a[i] * b[k - i];
            }
            high[k - n] = s;
        }
        return high;
    }

    private readonly record struct GarnerConsts(uint C0, uint C0s, uint C1, uint C1s, uint C2, uint C2s,
        uint Inv01, uint Inv01s, uint P0m2, uint P0m2s, uint Inv012, uint Inv012s);

    private static int GarnerVec<TV, TO>(ref uint r0, ref uint r1, ref uint r2, int k, int count, GarnerConsts g)
        where TV : struct where TO : struct, IVec<TV>
    {
        uint P0 = Primes[0].P, P1 = Primes[1].P, P2 = Primes[2].P;
        TV p0 = TO.Bc(P0), p1 = TO.Bc(P1), p2 = TO.Bc(P2);
        TV c0 = TO.Bc(g.C0), c0s = TO.Bc(g.C0s), c1 = TO.Bc(g.C1), c1s = TO.Bc(g.C1s), c2 = TO.Bc(g.C2), c2s = TO.Bc(g.C2s);
        TV i01 = TO.Bc(g.Inv01), i01s = TO.Bc(g.Inv01s), m2 = TO.Bc(g.P0m2), m2s = TO.Bc(g.P0m2s), i012 = TO.Bc(g.Inv012), i012s = TO.Bc(g.Inv012s);
        TV twoP1 = TO.Bc(2 * P1), fourP2 = TO.Bc(4 * P2);
        int W = TO.W;
        for (; k + W <= count; k += W)
        {
            nuint i = (nuint)k;
            TV x0 = ShoupMulV<TV, TO>(TO.Load(ref r0, i), c0, c0s, p0);
            x0 = TO.Min(x0, TO.Sub(x0, p0));
            TV x1 = ShoupMulV<TV, TO>(TO.Load(ref r1, i), c1, c1s, p1);
            x1 = TO.Min(x1, TO.Sub(x1, p1));
            TV x2 = ShoupMulV<TV, TO>(TO.Load(ref r2, i), c2, c2s, p2);
            x2 = TO.Min(x2, TO.Sub(x2, p2));
            TV v1 = ShoupMulV<TV, TO>(TO.Sub(TO.Add(x1, twoP1), x0), i01, i01s, p1);
            v1 = TO.Min(v1, TO.Sub(v1, p1));
            TV a0 = TO.Min(x0, TO.Sub(x0, p2));
            TV b0 = ShoupMulV<TV, TO>(v1, m2, m2s, p2);
            TV v2 = ShoupMulV<TV, TO>(TO.Sub(TO.Sub(TO.Add(x2, fourP2), a0), b0), i012, i012s, p2);
            v2 = TO.Min(v2, TO.Sub(v2, p2));
            TO.Store(ref r0, i, x0);
            TO.Store(ref r1, i, v1);
            TO.Store(ref r2, i, v2);
        }
        return k;
    }

    private static uint ScaleConst(Prime pr, int n)
    {
        uint r = (uint)((1UL << 32) % pr.P);
        uint ninv = Prime.PowMod((uint)n % pr.P, pr.P - 2, pr.P);
        return (uint)((ulong)r * ninv % pr.P);
    }

    // ---------------------------------------------------------------
    // Lane-width implementations
    // ---------------------------------------------------------------

    private readonly struct V512 : IVec<Vector512<uint>>
    {
        private static readonly Vector512<uint> OddLanes = Vector512.Create(0xFFFFFFFF_00000000UL).AsUInt32();
        private static readonly Vector512<ulong> MergeLo4 = Vector512.Create(0UL, 1, 8, 9, 2, 3, 10, 11);
        private static readonly Vector512<ulong> MergeHi4 = Vector512.Create(4UL, 5, 12, 13, 6, 7, 14, 15);

        public static int W => 16;
        public static Vector512<uint> Load(ref uint p, nuint i) => Vector512.LoadUnsafe(ref p, i);
        public static void Store(ref uint p, nuint i, Vector512<uint> v) => v.StoreUnsafe(ref p, i);
        public static Vector512<uint> Bc(uint x) => Vector512.Create(x);
        public static Vector512<uint> Create(ReadOnlySpan<uint> s) => Vector512.Create(s);
        public static Vector512<uint> Add(Vector512<uint> a, Vector512<uint> b) => a + b;
        public static Vector512<uint> Sub(Vector512<uint> a, Vector512<uint> b) => a - b;
        public static Vector512<uint> Min(Vector512<uint> a, Vector512<uint> b) => Vector512.Min(a, b);
        public static Vector512<uint> MulLo(Vector512<uint> a, Vector512<uint> b) => a * b;

        public static Vector512<uint> MulHi(Vector512<uint> a, Vector512<uint> b)
        {
            Vector512<ulong> e = Avx512F.Multiply(a, b);
            Vector512<ulong> o = Avx512F.Multiply((a.AsUInt64() >>> 32).AsUInt32(), (b.AsUInt64() >>> 32).AsUInt32());
            return Vector512.ConditionalSelect(OddLanes, o.AsUInt32(), (e >>> 32).AsUInt32());
        }

        public static Vector512<uint> MontMul(Vector512<uint> a, Vector512<uint> b, Vector512<uint> p, Vector512<uint> pinv)
        {
            Vector512<ulong> te = Avx512F.Multiply(a, b);
            Vector512<ulong> to = Avx512F.Multiply((a.AsUInt64() >>> 32).AsUInt32(), (b.AsUInt64() >>> 32).AsUInt32());
            Vector512<ulong> qe = Avx512F.Multiply(Avx512F.Multiply(te.AsUInt32(), pinv).AsUInt32(), p);
            Vector512<ulong> qo = Avx512F.Multiply(Avx512F.Multiply(to.AsUInt32(), pinv).AsUInt32(), p);
            Vector512<ulong> re = (te + qe) >>> 32;
            Vector512<ulong> ro = to + qo;
            return Vector512.ConditionalSelect(OddLanes, ro.AsUInt32(), re.AsUInt32());
        }

        public static void Split(ref Vector512<uint> a, ref Vector512<uint> b, int m)
        {
            Vector512<uint> x, y;
            switch (m)
            {
                case 8:
                    x = Avx512F.Shuffle4x128(a, b, 0x44);
                    y = Avx512F.Shuffle4x128(a, b, 0xEE);
                    break;
                case 4:
                    x = Avx512F.Shuffle4x128(a, b, 0x88);
                    y = Avx512F.Shuffle4x128(a, b, 0xDD);
                    break;
                case 2:
                    x = Avx512F.UnpackLow(a.AsUInt64(), b.AsUInt64()).AsUInt32();
                    y = Avx512F.UnpackHigh(a.AsUInt64(), b.AsUInt64()).AsUInt32();
                    break;
                default:
                    x = Avx512F.Shuffle(a.AsSingle(), b.AsSingle(), 0x88).AsUInt32();
                    y = Avx512F.Shuffle(a.AsSingle(), b.AsSingle(), 0xDD).AsUInt32();
                    break;
            }
            a = x;
            b = y;
        }

        public static void Merge(ref Vector512<uint> x, ref Vector512<uint> y, int m)
        {
            Vector512<uint> a, b;
            switch (m)
            {
                case 8:
                    a = Avx512F.Shuffle4x128(x, y, 0x44);
                    b = Avx512F.Shuffle4x128(x, y, 0xEE);
                    break;
                case 4:
                    a = Avx512F.PermuteVar8x64x2(x.AsUInt64(), MergeLo4, y.AsUInt64()).AsUInt32();
                    b = Avx512F.PermuteVar8x64x2(x.AsUInt64(), MergeHi4, y.AsUInt64()).AsUInt32();
                    break;
                case 2:
                    a = Avx512F.UnpackLow(x.AsUInt64(), y.AsUInt64()).AsUInt32();
                    b = Avx512F.UnpackHigh(x.AsUInt64(), y.AsUInt64()).AsUInt32();
                    break;
                default:
                    a = Avx512F.UnpackLow(x, y);
                    b = Avx512F.UnpackHigh(x, y);
                    break;
            }
            x = a;
            y = b;
        }
    }

    private readonly struct V256 : IVec<Vector256<uint>>
    {
        private static readonly Vector256<uint> OddLanes = Vector256.Create(0xFFFFFFFF_00000000UL).AsUInt32();

        public static int W => 8;
        public static Vector256<uint> Load(ref uint p, nuint i) => Vector256.LoadUnsafe(ref p, i);
        public static void Store(ref uint p, nuint i, Vector256<uint> v) => v.StoreUnsafe(ref p, i);
        public static Vector256<uint> Bc(uint x) => Vector256.Create(x);
        public static Vector256<uint> Create(ReadOnlySpan<uint> s) => Vector256.Create(s);
        public static Vector256<uint> Add(Vector256<uint> a, Vector256<uint> b) => a + b;
        public static Vector256<uint> Sub(Vector256<uint> a, Vector256<uint> b) => a - b;
        public static Vector256<uint> Min(Vector256<uint> a, Vector256<uint> b) => Vector256.Min(a, b);
        public static Vector256<uint> MulLo(Vector256<uint> a, Vector256<uint> b) => a * b;

        public static Vector256<uint> MulHi(Vector256<uint> a, Vector256<uint> b)
        {
            Vector256<ulong> e = Avx2.Multiply(a, b);
            Vector256<ulong> o = Avx2.Multiply((a.AsUInt64() >>> 32).AsUInt32(), (b.AsUInt64() >>> 32).AsUInt32());
            return Vector256.ConditionalSelect(OddLanes, o.AsUInt32(), (e >>> 32).AsUInt32());
        }

        public static Vector256<uint> MontMul(Vector256<uint> a, Vector256<uint> b, Vector256<uint> p, Vector256<uint> pinv)
        {
            Vector256<ulong> te = Avx2.Multiply(a, b);
            Vector256<ulong> to = Avx2.Multiply((a.AsUInt64() >>> 32).AsUInt32(), (b.AsUInt64() >>> 32).AsUInt32());
            Vector256<ulong> qe = Avx2.Multiply(Avx2.Multiply(te.AsUInt32(), pinv).AsUInt32(), p);
            Vector256<ulong> qo = Avx2.Multiply(Avx2.Multiply(to.AsUInt32(), pinv).AsUInt32(), p);
            Vector256<ulong> re = (te + qe) >>> 32;
            Vector256<ulong> ro = to + qo;
            return Vector256.ConditionalSelect(OddLanes, ro.AsUInt32(), re.AsUInt32());
        }

        public static void Split(ref Vector256<uint> a, ref Vector256<uint> b, int m)
        {
            Vector256<uint> x, y;
            switch (m)
            {
                case 4:
                    x = Avx2.Permute2x128(a, b, 0x20);
                    y = Avx2.Permute2x128(a, b, 0x31);
                    break;
                case 2:
                    x = Avx2.UnpackLow(a.AsUInt64(), b.AsUInt64()).AsUInt32();
                    y = Avx2.UnpackHigh(a.AsUInt64(), b.AsUInt64()).AsUInt32();
                    break;
                default:
                    x = Avx.Shuffle(a.AsSingle(), b.AsSingle(), 0x88).AsUInt32();
                    y = Avx.Shuffle(a.AsSingle(), b.AsSingle(), 0xDD).AsUInt32();
                    break;
            }
            a = x;
            b = y;
        }

        public static void Merge(ref Vector256<uint> x, ref Vector256<uint> y, int m)
        {
            Vector256<uint> a, b;
            switch (m)
            {
                case 4:
                    a = Avx2.Permute2x128(x, y, 0x20);
                    b = Avx2.Permute2x128(x, y, 0x31);
                    break;
                case 2:
                    a = Avx2.UnpackLow(x.AsUInt64(), y.AsUInt64()).AsUInt32();
                    b = Avx2.UnpackHigh(x.AsUInt64(), y.AsUInt64()).AsUInt32();
                    break;
                default:
                    a = Avx2.UnpackLow(x, y);
                    b = Avx2.UnpackHigh(x, y);
                    break;
            }
            x = a;
            y = b;
        }
    }

    private readonly struct V128 : IVec<Vector128<uint>>
    {
        private static readonly Vector128<uint> OddLanes = Vector128.Create(0xFFFFFFFF_00000000UL).AsUInt32();

        public static int W => 4;
        public static Vector128<uint> Load(ref uint p, nuint i) => Vector128.LoadUnsafe(ref p, i);
        public static void Store(ref uint p, nuint i, Vector128<uint> v) => v.StoreUnsafe(ref p, i);
        public static Vector128<uint> Bc(uint x) => Vector128.Create(x);
        public static Vector128<uint> Create(ReadOnlySpan<uint> s) => Vector128.Create(s);
        public static Vector128<uint> Add(Vector128<uint> a, Vector128<uint> b) => a + b;
        public static Vector128<uint> Sub(Vector128<uint> a, Vector128<uint> b) => a - b;
        public static Vector128<uint> Min(Vector128<uint> a, Vector128<uint> b) => Vector128.Min(a, b);
        public static Vector128<uint> MulLo(Vector128<uint> a, Vector128<uint> b) => a * b;

        /// <summary>64-bit products of lanes 0 and 2 (x86 PMULUDQ semantics).</summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static Vector128<ulong> MulEven(Vector128<uint> a, Vector128<uint> b)
        {
            if (Sse2.IsSupported)
            {
                return Sse2.Multiply(a, b);
            }
            // Arm64: gather even lanes, then widening multiply.
            Vector64<uint> ae = AdvSimd.Arm64.UnzipEven(a, a).GetLower();
            Vector64<uint> be = AdvSimd.Arm64.UnzipEven(b, b).GetLower();
            return AdvSimd.MultiplyWideningLower(ae, be);
        }

        public static Vector128<uint> MulHi(Vector128<uint> a, Vector128<uint> b)
        {
            Vector128<ulong> e = MulEven(a, b);
            Vector128<ulong> o = MulEven((a.AsUInt64() >>> 32).AsUInt32(), (b.AsUInt64() >>> 32).AsUInt32());
            return Vector128.ConditionalSelect(OddLanes, o.AsUInt32(), (e >>> 32).AsUInt32());
        }

        public static Vector128<uint> MontMul(Vector128<uint> a, Vector128<uint> b, Vector128<uint> p, Vector128<uint> pinv)
        {
            Vector128<ulong> te = MulEven(a, b);
            Vector128<ulong> to = MulEven((a.AsUInt64() >>> 32).AsUInt32(), (b.AsUInt64() >>> 32).AsUInt32());
            Vector128<ulong> qe = MulEven((te.AsUInt32() * pinv), p);
            Vector128<ulong> qo = MulEven((to.AsUInt32() * pinv), p);
            Vector128<ulong> re = (te + qe) >>> 32;
            Vector128<ulong> ro = to + qo;
            return Vector128.ConditionalSelect(OddLanes, ro.AsUInt32(), re.AsUInt32());
        }

        public static void Split(ref Vector128<uint> a, ref Vector128<uint> b, int m)
        {
            Vector128<uint> x, y;
            if (m == 2)
            {
                // [a0 a1 b0 b1] / [a2 a3 b2 b3]
                if (Sse2.IsSupported)
                {
                    x = Sse2.UnpackLow(a.AsUInt64(), b.AsUInt64()).AsUInt32();
                    y = Sse2.UnpackHigh(a.AsUInt64(), b.AsUInt64()).AsUInt32();
                }
                else
                {
                    x = AdvSimd.Arm64.ZipLow(a.AsUInt64(), b.AsUInt64()).AsUInt32();
                    y = AdvSimd.Arm64.ZipHigh(a.AsUInt64(), b.AsUInt64()).AsUInt32();
                }
            }
            else if (Sse.IsSupported)
            {
                x = Sse.Shuffle(a.AsSingle(), b.AsSingle(), 0x88).AsUInt32();
                y = Sse.Shuffle(a.AsSingle(), b.AsSingle(), 0xDD).AsUInt32();
            }
            else
            {
                x = AdvSimd.Arm64.UnzipEven(a, b);
                y = AdvSimd.Arm64.UnzipOdd(a, b);
            }
            a = x;
            b = y;
        }

        public static void Merge(ref Vector128<uint> x, ref Vector128<uint> y, int m)
        {
            Vector128<uint> a, b;
            if (m == 2)
            {
                if (Sse2.IsSupported)
                {
                    a = Sse2.UnpackLow(x.AsUInt64(), y.AsUInt64()).AsUInt32();
                    b = Sse2.UnpackHigh(x.AsUInt64(), y.AsUInt64()).AsUInt32();
                }
                else
                {
                    a = AdvSimd.Arm64.ZipLow(x.AsUInt64(), y.AsUInt64()).AsUInt32();
                    b = AdvSimd.Arm64.ZipHigh(x.AsUInt64(), y.AsUInt64()).AsUInt32();
                }
            }
            else if (Sse2.IsSupported)
            {
                a = Sse2.UnpackLow(x, y);
                b = Sse2.UnpackHigh(x, y);
            }
            else
            {
                a = AdvSimd.Arm64.ZipLow(x, y);
                b = AdvSimd.Arm64.ZipHigh(x, y);
            }
            x = a;
            y = b;
        }
    }
}
