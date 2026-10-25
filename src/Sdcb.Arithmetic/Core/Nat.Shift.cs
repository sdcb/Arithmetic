using System.Diagnostics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;

namespace Sdcb.Arithmetic.Core;

internal static partial class Nat
{
    /// <summary>
    /// r[0..n) = a[0..n) &lt;&lt; cnt (1 &lt;= cnt &lt;= 63). Returns the bits shifted out of the top.
    /// Operates high-to-low so r may alias a (or sit at a higher address).
    /// </summary>
    public static ulong ShiftLeft(Span<ulong> r, ReadOnlySpan<ulong> a, int n, int cnt)
    {
        Debug.Assert(cnt >= 1 && cnt <= 63 && n >= 1);
        if ((uint)n > (uint)r.Length || (uint)n > (uint)a.Length)
        {
            throw new ArgumentOutOfRangeException(nameof(n));
        }
        ref ulong rp = ref MemoryMarshal.GetReference(r);
        ref ulong ap = ref MemoryMarshal.GetReference(a);
        int tnc = 64 - cnt;
        ulong ret = Unsafe.Add(ref ap, n - 1) >> tnc;
        int i = n - 1; // next output index to produce (descending), i >= 1 uses a[i-1]
        if (Isa.Hw512)
        {
            for (; i >= 8; i -= 8)
            {
                nuint j = (nuint)(i - 7);
                Vector512<ulong> hi = Vector512.LoadUnsafe(ref ap, j), lo = Vector512.LoadUnsafe(ref ap, j - 1);
                (Vector512.ShiftLeft(hi, cnt) | Vector512.ShiftRightLogical(lo, tnc)).StoreUnsafe(ref rp, j);
            }
        }
        if (Isa.Hw256)
        {
            for (; i >= 4; i -= 4)
            {
                nuint j = (nuint)(i - 3);
                Vector256<ulong> hi = Vector256.LoadUnsafe(ref ap, j), lo = Vector256.LoadUnsafe(ref ap, j - 1);
                (Vector256.ShiftLeft(hi, cnt) | Vector256.ShiftRightLogical(lo, tnc)).StoreUnsafe(ref rp, j);
            }
        }
        if (Isa.Hw128)
        {
            for (; i >= 2; i -= 2)
            {
                nuint j = (nuint)(i - 1);
                Vector128<ulong> hi = Vector128.LoadUnsafe(ref ap, j), lo = Vector128.LoadUnsafe(ref ap, j - 1);
                (Vector128.ShiftLeft(hi, cnt) | Vector128.ShiftRightLogical(lo, tnc)).StoreUnsafe(ref rp, j);
            }
        }
        for (; i > 0; i--)
        {
            Unsafe.Add(ref rp, i) = (Unsafe.Add(ref ap, i) << cnt) | (Unsafe.Add(ref ap, i - 1) >> tnc);
        }
        rp = ap << cnt;
        return ret;
    }

    /// <summary>
    /// r[0..n) = a[0..n) &gt;&gt; cnt (1 &lt;= cnt &lt;= 63). Returns the bits shifted out of the
    /// bottom, placed in the most significant bits of the return value.
    /// Operates low-to-high so r may alias a (or sit at a lower address).
    /// </summary>
    public static ulong ShiftRight(Span<ulong> r, ReadOnlySpan<ulong> a, int n, int cnt)
    {
        Debug.Assert(cnt >= 1 && cnt <= 63 && n >= 1);
        if ((uint)n > (uint)r.Length || (uint)n > (uint)a.Length)
        {
            throw new ArgumentOutOfRangeException(nameof(n));
        }
        ref ulong rp = ref MemoryMarshal.GetReference(r);
        ref ulong ap = ref MemoryMarshal.GetReference(a);
        int tnc = 64 - cnt;
        ulong ret = ap << tnc;
        int last = n - 1; // outputs [0, last) use a[i+1]
        int i = 0;
        if (Isa.Hw512)
        {
            for (; i + 8 <= last; i += 8)
            {
                Vector512<ulong> lo = Vector512.LoadUnsafe(ref ap, (nuint)i), hi = Vector512.LoadUnsafe(ref ap, (nuint)i + 1);
                (Vector512.ShiftRightLogical(lo, cnt) | Vector512.ShiftLeft(hi, tnc)).StoreUnsafe(ref rp, (nuint)i);
            }
        }
        if (Isa.Hw256)
        {
            for (; i + 4 <= last; i += 4)
            {
                Vector256<ulong> lo = Vector256.LoadUnsafe(ref ap, (nuint)i), hi = Vector256.LoadUnsafe(ref ap, (nuint)i + 1);
                (Vector256.ShiftRightLogical(lo, cnt) | Vector256.ShiftLeft(hi, tnc)).StoreUnsafe(ref rp, (nuint)i);
            }
        }
        if (Isa.Hw128)
        {
            for (; i + 2 <= last; i += 2)
            {
                Vector128<ulong> lo = Vector128.LoadUnsafe(ref ap, (nuint)i), hi = Vector128.LoadUnsafe(ref ap, (nuint)i + 1);
                (Vector128.ShiftRightLogical(lo, cnt) | Vector128.ShiftLeft(hi, tnc)).StoreUnsafe(ref rp, (nuint)i);
            }
        }
        for (; i < last; i++)
        {
            Unsafe.Add(ref rp, i) = (Unsafe.Add(ref ap, i) >> cnt) | (Unsafe.Add(ref ap, i + 1) << tnc);
        }
        Unsafe.Add(ref rp, last) = Unsafe.Add(ref ap, last) >> cnt;
        return ret;
    }

    // ---------------------------------------------------------------
    // Comparison / normalization
    // ---------------------------------------------------------------

    /// <summary>Compares equal-length magnitudes. Returns -1/0/1.</summary>
    public static int Compare(ReadOnlySpan<ulong> a, ReadOnlySpan<ulong> b, int n)
    {
        for (int i = n - 1; i >= 0; i--)
        {
            ulong x = a[i], y = b[i];
            if (x != y)
            {
                return x < y ? -1 : 1;
            }
        }
        return 0;
    }

    /// <summary>Compares a[0..an) against b[0..bn) with an &gt;= bn. Returns -1/0/1.</summary>
    public static int CompareLonger(ReadOnlySpan<ulong> a, int an, ReadOnlySpan<ulong> b, int bn)
    {
        Debug.Assert(an >= bn);
        for (int i = an - 1; i >= bn; i--)
        {
            if (a[i] != 0)
            {
                return 1;
            }
        }
        return Compare(a, b, bn);
    }

    /// <summary>Compares fully-normalized magnitudes of possibly different lengths.</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int CompareNormalized(ReadOnlySpan<ulong> a, int an, ReadOnlySpan<ulong> b, int bn)
    {
        if (an != bn)
        {
            return an < bn ? -1 : 1;
        }
        return Compare(a, b, an);
    }

    /// <summary>Returns the number of significant limbs (strips high zero limbs).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static int TrimmedLength(ReadOnlySpan<ulong> a, int n)
    {
        while (n > 0 && a[n - 1] == 0)
        {
            n--;
        }
        return n;
    }

    /// <summary>Returns true if all n limbs are zero.</summary>
    public static bool IsZero(ReadOnlySpan<ulong> a, int n) => a.Slice(0, n).IndexOfAnyExcept(0UL) < 0;
}
