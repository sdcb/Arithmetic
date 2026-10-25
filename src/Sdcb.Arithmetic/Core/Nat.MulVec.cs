using System.Buffers;
using System.Diagnostics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.X86;

namespace Sdcb.Arithmetic.Core;

internal static partial class Nat
{
    /// <summary>Smaller-operand size (limbs) from which the AVX2 basecase beats the scalar one.</summary>
    internal static readonly int MulVecMin = EnvInt("SDCB_ARITHMETIC_MULVEC") ?? 8;

    /// <summary>Size (limbs) from which squaring uses the AVX2 product instead of the scalar squaring basecase.</summary>
    internal static readonly int SqrVecMin = EnvInt("SDCB_ARITHMETIC_SQRVEC") ?? 12;

    private const int VecPad = 16;   // zero digits on each side of the expanded operand
    private const int VecCols = 16;  // output digit columns per block (four Vector256<ulong>)

    /// <summary>
    /// Schoolbook product with AVX2 (vpmuludq). Operands are viewed as 32-bit digits;
    /// a is widened to one digit per 64-bit lane and every digit of b is split into 16-bit
    /// halves, so each 32x16-bit product is &lt; 2^48 and a column of up to 2^16 of them
    /// accumulates in a 64-bit lane without carries. Each block of 16 output columns is
    /// summed in registers and folded into r with one scalar carry pass.
    /// Requires an &gt;= bn &gt;= 1, bn &lt; 2^14 and r not overlapping a or b.
    /// </summary>
    internal static void MulSchoolbookVec256(Span<ulong> r, ReadOnlySpan<ulong> a, int an, ReadOnlySpan<ulong> b, int bn)
        => MulVec256(r, a, an, b, bn, an + bn);

    /// <summary>r[0..n) = (a * b) mod B^n for n-limb a and b (low half only, AVX2).</summary>
    internal static void MulLowVec256(Span<ulong> r, ReadOnlySpan<ulong> a, ReadOnlySpan<ulong> b, int n)
        => MulVec256(r, a, n, b, n, n);

    /// <summary>r[0..rn) = low rn limbs of a * b (rn &lt;= an + bn).</summary>
    [SkipLocalsInit]
    private static void MulVec256(Span<ulong> r, ReadOnlySpan<ulong> a, int an, ReadOnlySpan<ulong> b, int bn, int rn)
    {
        Debug.Assert(Avx2.IsSupported && an >= bn && bn >= 1 && bn < 1 << 14 && rn <= an + bn);
        int na = 2 * an, nb = 2 * bn;
        int need = na + 2 * VecPad + 2 * nb;
        ulong[]? rented = null;
        Span<ulong> buf = need <= 1024
            ? stackalloc ulong[1024]
            : (rented = ArrayPool<ulong>.Shared.Rent(need));

        ref ulong pad = ref MemoryMarshal.GetReference(buf);
        ref ulong ax = ref Unsafe.Add(ref pad, VecPad);           // ax[-VecPad .. na + VecPad)
        ref ulong bs = ref Unsafe.Add(ref ax, na + VecPad);       // (lo16, hi16) per digit of b
        buf.Slice(0, VecPad).Clear();
        buf.Slice(VecPad + na, VecPad).Clear();

        ReadOnlySpan<uint> a32 = MemoryMarshal.Cast<ulong, uint>(a.Slice(0, an));
        ref uint a32p = ref MemoryMarshal.GetReference(a32);
        nint i = 0;
        for (; i + 4 <= na; i += 4)
        {
            Avx2.ConvertToVector256Int64(Vector128.LoadUnsafe(ref a32p, (nuint)i)).AsUInt64().StoreUnsafe(ref ax, (nuint)i);
        }
        for (; i < na; i++)
        {
            Unsafe.Add(ref ax, i) = Unsafe.Add(ref a32p, i);
        }
        ReadOnlySpan<uint> b32 = MemoryMarshal.Cast<ulong, uint>(b.Slice(0, bn));
        for (int j = 0; j < nb; j++)
        {
            uint d = b32[j];
            Unsafe.Add(ref bs, 2 * j) = d & 0xFFFF;
            Unsafe.Add(ref bs, 2 * j + 1) = d >> 16;
        }

        Span<uint> r32 = MemoryMarshal.Cast<ulong, uint>(r.Slice(0, rn));
        int nr = 2 * rn;
        Span<ulong> fold = stackalloc ulong[2 * VecCols];
        ref ulong fd = ref MemoryMarshal.GetReference(fold);
        Vector256<ulong> m16 = Vector256.Create(0xFFFFUL), m32 = Vector256.Create(0xFFFFFFFFUL);
        ulong carry = 0, prev = 0;
        for (int k = 0; k < nr; k += VecCols)
        {
            // Columns k..k+15; lane t of column k+t takes a digit k+t-j.
            int jlo = Math.Max(0, k - (na - 1));
            int jhi = Math.Min(nb - 1, k + VecCols - 1);
            Vector256<ulong> l0 = default, l1 = default, l2 = default, l3 = default;
            Vector256<ulong> h0 = default, h1 = default, h2 = default, h3 = default;
            ref ulong ap = ref Unsafe.Add(ref ax, k - jlo);
            ref ulong bp = ref Unsafe.Add(ref bs, 2 * jlo);
            for (int j = jlo; j <= jhi; j++)
            {
                Vector256<uint> bl = Vector256.Create(bp).AsUInt32();
                Vector256<uint> bh = Vector256.Create(Unsafe.Add(ref bp, 1)).AsUInt32();
                Vector256<uint> x0 = Vector256.LoadUnsafe(ref ap).AsUInt32();
                Vector256<uint> x1 = Vector256.LoadUnsafe(ref ap, 4).AsUInt32();
                Vector256<uint> x2 = Vector256.LoadUnsafe(ref ap, 8).AsUInt32();
                Vector256<uint> x3 = Vector256.LoadUnsafe(ref ap, 12).AsUInt32();
                l0 += Avx2.Multiply(x0, bl);
                h0 += Avx2.Multiply(x0, bh);
                l1 += Avx2.Multiply(x1, bl);
                h1 += Avx2.Multiply(x1, bh);
                l2 += Avx2.Multiply(x2, bl);
                h2 += Avx2.Multiply(x2, bh);
                l3 += Avx2.Multiply(x3, bl);
                h3 += Avx2.Multiply(x3, bh);
                ap = ref Unsafe.Subtract(ref ap, 1);
                bp = ref Unsafe.Add(ref bp, 2);
            }
            // Column value l + h*2^16 = lo + hi*2^32 with lo < 2^33 at column t, hi at column t+1.
            FoldColumns(l0, h0, m16, m32, ref fd, 0);
            FoldColumns(l1, h1, m16, m32, ref fd, 4);
            FoldColumns(l2, h2, m16, m32, ref fd, 8);
            FoldColumns(l3, h3, m16, m32, ref fd, 12);
            int lim = Math.Min(VecCols, nr - k);
            for (int t = 0; t < lim; t++)
            {
                ulong c = Unsafe.Add(ref fd, t) + prev + carry;
                r32[k + t] = (uint)c;
                carry = c >> 32;
                prev = Unsafe.Add(ref fd, VecCols + t);
            }
        }
        Debug.Assert(rn < an + bn || (carry == 0 && prev == 0));
        if (rented is not null)
        {
            ArrayPool<ulong>.Shared.Return(rented);
        }
    }

    /// <summary>
    /// Lane masks for the squaring kernel: entry m (0..16) keeps lanes t &lt; m of a 16-column block.
    /// </summary>
    private static readonly ulong[] s_sqrMasks = BuildSqrMasks();

    private static ulong[] BuildSqrMasks()
    {
        ulong[] m = new ulong[(VecCols + 1) * VecCols];
        for (int k = 0; k <= VecCols; k++)
        {
            for (int t = 0; t < k; t++)
            {
                m[k * VecCols + t] = ulong.MaxValue;
            }
        }
        return m;
    }

    /// <summary>
    /// Schoolbook square with AVX2, same digit layout as <see cref="MulSchoolbookVec256"/>: only the
    /// products a_i a_j with i &lt; j are summed (doubled during the fold) and the diagonal squares
    /// are added in the scalar carry pass. Requires 1 &lt;= n &lt; 2^13 and r not overlapping a.
    /// </summary>
    [SkipLocalsInit]
    internal static void SquareSchoolbookVec256(Span<ulong> r, ReadOnlySpan<ulong> a, int n)
    {
        Debug.Assert(Avx2.IsSupported && n >= 1 && n < 1 << 13);
        int na = 2 * n;
        int need = na + 2 * VecPad + 2 * na;
        ulong[]? rented = null;
        Span<ulong> buf = need <= 1024
            ? stackalloc ulong[1024]
            : (rented = ArrayPool<ulong>.Shared.Rent(need));

        ref ulong pad = ref MemoryMarshal.GetReference(buf);
        ref ulong ax = ref Unsafe.Add(ref pad, VecPad);
        ref ulong bs = ref Unsafe.Add(ref ax, na + VecPad);
        buf.Slice(0, VecPad).Clear();
        buf.Slice(VecPad + na, VecPad).Clear();

        ReadOnlySpan<uint> a32 = MemoryMarshal.Cast<ulong, uint>(a.Slice(0, n));
        ref uint a32p = ref MemoryMarshal.GetReference(a32);
        nint i = 0;
        for (; i + 4 <= na; i += 4)
        {
            Avx2.ConvertToVector256Int64(Vector128.LoadUnsafe(ref a32p, (nuint)i)).AsUInt64().StoreUnsafe(ref ax, (nuint)i);
        }
        for (; i < na; i++)
        {
            Unsafe.Add(ref ax, i) = Unsafe.Add(ref a32p, i);
        }
        for (int j = 0; j < na; j++)
        {
            uint d = a32[j];
            Unsafe.Add(ref bs, 2 * j) = d & 0xFFFF;
            Unsafe.Add(ref bs, 2 * j + 1) = d >> 16;
        }

        Span<uint> r32 = MemoryMarshal.Cast<ulong, uint>(r.Slice(0, 2 * n));
        int nr = 2 * na;
        Span<ulong> fold = stackalloc ulong[2 * VecCols];
        ref ulong fd = ref MemoryMarshal.GetReference(fold);
        ref ulong masks = ref MemoryMarshal.GetArrayDataReference(s_sqrMasks);
        Vector256<ulong> m16 = Vector256.Create(0xFFFFUL), m32 = Vector256.Create(0xFFFFFFFFUL);
        ulong carry = 0, prev = 0;
        for (int k = 0; k < nr; k += VecCols)
        {
            // Column c = k+t takes a_{c-j} a_j for j > c-j, i.e. lanes t < 2j-k.
            int jlo = Math.Max(k / 2 + 1, k - (na - 1));
            int jhi = Math.Min(na - 1, k + VecCols - 1);
            // From jfull on every lane qualifies (2j - k >= 16).
            int jfull = Math.Max(jlo, (k + VecCols + 1) / 2);
            Vector256<ulong> l0 = default, l1 = default, l2 = default, l3 = default;
            Vector256<ulong> h0 = default, h1 = default, h2 = default, h3 = default;
            int j = jlo;
            for (; j < jfull && j <= jhi; j++)
            {
                ref ulong ap = ref Unsafe.Add(ref ax, k - j);
                ref ulong mk = ref Unsafe.Add(ref masks, (2 * j - k) * VecCols);
                Vector256<uint> bl = Vector256.Create(Unsafe.Add(ref bs, 2 * j)).AsUInt32();
                Vector256<uint> bh = Vector256.Create(Unsafe.Add(ref bs, 2 * j + 1)).AsUInt32();
                Vector256<uint> x0 = (Vector256.LoadUnsafe(ref ap) & Vector256.LoadUnsafe(ref mk)).AsUInt32();
                Vector256<uint> x1 = (Vector256.LoadUnsafe(ref ap, 4) & Vector256.LoadUnsafe(ref mk, 4)).AsUInt32();
                Vector256<uint> x2 = (Vector256.LoadUnsafe(ref ap, 8) & Vector256.LoadUnsafe(ref mk, 8)).AsUInt32();
                Vector256<uint> x3 = (Vector256.LoadUnsafe(ref ap, 12) & Vector256.LoadUnsafe(ref mk, 12)).AsUInt32();
                l0 += Avx2.Multiply(x0, bl);
                h0 += Avx2.Multiply(x0, bh);
                l1 += Avx2.Multiply(x1, bl);
                h1 += Avx2.Multiply(x1, bh);
                l2 += Avx2.Multiply(x2, bl);
                h2 += Avx2.Multiply(x2, bh);
                l3 += Avx2.Multiply(x3, bl);
                h3 += Avx2.Multiply(x3, bh);
            }
            if (j <= jhi)
            {
                ref ulong ap = ref Unsafe.Add(ref ax, k - j);
                ref ulong bp = ref Unsafe.Add(ref bs, 2 * j);
                for (; j <= jhi; j++)
                {
                    Vector256<uint> bl = Vector256.Create(bp).AsUInt32();
                    Vector256<uint> bh = Vector256.Create(Unsafe.Add(ref bp, 1)).AsUInt32();
                    Vector256<uint> x0 = Vector256.LoadUnsafe(ref ap).AsUInt32();
                    Vector256<uint> x1 = Vector256.LoadUnsafe(ref ap, 4).AsUInt32();
                    Vector256<uint> x2 = Vector256.LoadUnsafe(ref ap, 8).AsUInt32();
                    Vector256<uint> x3 = Vector256.LoadUnsafe(ref ap, 12).AsUInt32();
                    l0 += Avx2.Multiply(x0, bl);
                    h0 += Avx2.Multiply(x0, bh);
                    l1 += Avx2.Multiply(x1, bl);
                    h1 += Avx2.Multiply(x1, bh);
                    l2 += Avx2.Multiply(x2, bl);
                    h2 += Avx2.Multiply(x2, bh);
                    l3 += Avx2.Multiply(x3, bl);
                    h3 += Avx2.Multiply(x3, bh);
                    ap = ref Unsafe.Subtract(ref ap, 1);
                    bp = ref Unsafe.Add(ref bp, 2);
                }
            }
            // Double the off-diagonal sums (each lane stays below 2^63).
            FoldColumns(l0 << 1, h0 << 1, m16, m32, ref fd, 0);
            FoldColumns(l1 << 1, h1 << 1, m16, m32, ref fd, 4);
            FoldColumns(l2 << 1, h2 << 1, m16, m32, ref fd, 8);
            FoldColumns(l3 << 1, h3 << 1, m16, m32, ref fd, 12);
            int lim = Math.Min(VecCols, nr - k);
            for (int t = 0; t < lim; t++)
            {
                int c = k + t;
                ulong v = Unsafe.Add(ref fd, t) + prev + carry;
                ulong next = Unsafe.Add(ref fd, VecCols + t);
                if ((c & 1) == 0)
                {
                    ulong x = Unsafe.Add(ref ax, c >> 1);
                    ulong sq = x * x;
                    v += (uint)sq;
                    next += sq >> 32;
                }
                r32[c] = (uint)v;
                carry = v >> 32;
                prev = next;
            }
        }
        Debug.Assert(carry == 0 && prev == 0);
        if (rented is not null)
        {
            ArrayPool<ulong>.Shared.Return(rented);
        }
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static void FoldColumns(Vector256<ulong> l, Vector256<ulong> h, Vector256<ulong> m16, Vector256<ulong> m32, ref ulong fd, nuint t)
    {
        ((l & m32) + ((h & m16) << 16)).StoreUnsafe(ref fd, t);
        ((l >>> 32) + (h >>> 16)).StoreUnsafe(ref fd, VecCols + t);
    }
}
