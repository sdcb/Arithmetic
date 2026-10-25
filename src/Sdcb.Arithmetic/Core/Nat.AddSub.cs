using System.Diagnostics;
using System.Runtime.CompilerServices;
using System.Runtime.InteropServices;
using System.Runtime.Intrinsics;

namespace Sdcb.Arithmetic.Core;

/// <summary>
/// Low-level multi-precision natural number kernels operating on little-endian
/// spans of 64-bit limbs (index 0 = least significant). The natural-number layer under BigInt, BigRational and Real.
/// All routines assume the caller has validated lengths; they are the hot inner
/// loops of the library.
/// </summary>
internal static partial class Nat
{
    // ---------------------------------------------------------------
    // Addition
    // ---------------------------------------------------------------

    /// <summary>r[0..n) = a[0..n) + b[0..n), returns carry (0 or 1). r may alias a or b.</summary>
    public static ulong Add(Span<ulong> r, ReadOnlySpan<ulong> a, ReadOnlySpan<ulong> b, int n)
    {
        Debug.Assert(n >= 0 && r.Length >= n && a.Length >= n && b.Length >= n);
        ref ulong rp = ref MemoryMarshal.GetReference(r);
        ref ulong ap = ref MemoryMarshal.GetReference(a);
        ref ulong bp = ref MemoryMarshal.GetReference(b);
        ulong carry = 0;
        int i = 0;
        if (n >= VecAddMin)
        {
            if (Isa.Hw512)
            {
                carry = AddVec512(ref rp, ref ap, ref bp, n, ref i, (uint)carry);
            }
            if (Isa.Hw256)
            {
                carry = AddVec256(ref rp, ref ap, ref bp, n, ref i, (uint)carry);
            }
        }
        // 4x unrolled carry chain.
        for (; i + 4 <= n; i += 4)
        {
            carry = AddWithCarry(Unsafe.Add(ref ap, i), Unsafe.Add(ref bp, i), carry, out Unsafe.Add(ref rp, i));
            carry = AddWithCarry(Unsafe.Add(ref ap, i + 1), Unsafe.Add(ref bp, i + 1), carry, out Unsafe.Add(ref rp, i + 1));
            carry = AddWithCarry(Unsafe.Add(ref ap, i + 2), Unsafe.Add(ref bp, i + 2), carry, out Unsafe.Add(ref rp, i + 2));
            carry = AddWithCarry(Unsafe.Add(ref ap, i + 3), Unsafe.Add(ref bp, i + 3), carry, out Unsafe.Add(ref rp, i + 3));
        }
        for (; i < n; i++)
        {
            carry = AddWithCarry(Unsafe.Add(ref ap, i), Unsafe.Add(ref bp, i), carry, out Unsafe.Add(ref rp, i));
        }
        return carry;
    }

    /// <summary>r[0..an) = a[0..an) + b[0..bn) with an &gt;= bn. Returns carry.</summary>
    public static ulong Add(Span<ulong> r, ReadOnlySpan<ulong> a, int an, ReadOnlySpan<ulong> b, int bn)
    {
        Debug.Assert(an >= bn);
        ulong carry = Add(r, a, b, bn);
        if (an > bn)
        {
            carry = AddScalar(r.Slice(bn), a.Slice(bn), an - bn, carry);
        }
        return carry;
    }

    /// <summary>r[0..n) = a[0..n) + b (single limb). Returns carry.</summary>
    public static ulong AddScalar(Span<ulong> r, ReadOnlySpan<ulong> a, int n, ulong b)
    {
        int i = 0;
        for (; i < n; i++)
        {
            ulong s = a[i] + b;
            r[i] = s;
            b = s < b ? 1UL : 0UL;
            if (b == 0) { i++; break; }
        }
        if (i < n)
        {
            a.Slice(i, n - i).CopyTo(r.Slice(i));
        }
        return b;
    }

    /// <summary>In-place r[0..n) += b (single limb). Returns carry.</summary>
    public static ulong AddScalarInPlace(Span<ulong> r, int n, ulong b)
    {
        for (int i = 0; i < n && b != 0; i++)
        {
            ulong s = r[i] + b;
            r[i] = s;
            b = s < b ? 1UL : 0UL;
        }
        return b;
    }

    /// <summary>In-place r[0..rn) += a[0..an), rn &gt;= an. Returns carry.</summary>
    public static ulong AddInPlace(Span<ulong> r, int rn, ReadOnlySpan<ulong> a, int an)
    {
        Debug.Assert(rn >= an);
        ulong carry = Add(r, r, a, an);
        if (carry != 0 && rn > an)
        {
            carry = AddScalarInPlace(r.Slice(an), rn - an, carry);
        }
        return carry;
    }

    // ---------------------------------------------------------------
    // Subtraction
    // ---------------------------------------------------------------

    /// <summary>r[0..n) = a[0..n) - b[0..n), returns borrow (0 or 1). r may alias.</summary>
    public static ulong Sub(Span<ulong> r, ReadOnlySpan<ulong> a, ReadOnlySpan<ulong> b, int n)
    {
        Debug.Assert(n >= 0 && r.Length >= n && a.Length >= n && b.Length >= n);
        ref ulong rp = ref MemoryMarshal.GetReference(r);
        ref ulong ap = ref MemoryMarshal.GetReference(a);
        ref ulong bp = ref MemoryMarshal.GetReference(b);
        ulong borrow = 0;
        int i = 0;
        if (n >= VecAddMin)
        {
            if (Isa.Hw512)
            {
                borrow = SubVec512(ref rp, ref ap, ref bp, n, ref i, (uint)borrow);
            }
            if (Isa.Hw256)
            {
                borrow = SubVec256(ref rp, ref ap, ref bp, n, ref i, (uint)borrow);
            }
        }
        for (; i + 4 <= n; i += 4)
        {
            borrow = SubWithBorrow(Unsafe.Add(ref ap, i), Unsafe.Add(ref bp, i), borrow, out Unsafe.Add(ref rp, i));
            borrow = SubWithBorrow(Unsafe.Add(ref ap, i + 1), Unsafe.Add(ref bp, i + 1), borrow, out Unsafe.Add(ref rp, i + 1));
            borrow = SubWithBorrow(Unsafe.Add(ref ap, i + 2), Unsafe.Add(ref bp, i + 2), borrow, out Unsafe.Add(ref rp, i + 2));
            borrow = SubWithBorrow(Unsafe.Add(ref ap, i + 3), Unsafe.Add(ref bp, i + 3), borrow, out Unsafe.Add(ref rp, i + 3));
        }
        for (; i < n; i++)
        {
            borrow = SubWithBorrow(Unsafe.Add(ref ap, i), Unsafe.Add(ref bp, i), borrow, out Unsafe.Add(ref rp, i));
        }
        return borrow;
    }

    /// <summary>r[0..an) = a[0..an) - b[0..bn), an &gt;= bn. Returns borrow.</summary>
    public static ulong Sub(Span<ulong> r, ReadOnlySpan<ulong> a, int an, ReadOnlySpan<ulong> b, int bn)
    {
        Debug.Assert(an >= bn);
        ulong borrow = Sub(r, a, b, bn);
        if (an > bn)
        {
            borrow = SubScalar(r.Slice(bn), a.Slice(bn), an - bn, borrow);
        }
        return borrow;
    }

    /// <summary>r[0..n) = a[0..n) - b (single limb). Returns borrow.</summary>
    public static ulong SubScalar(Span<ulong> r, ReadOnlySpan<ulong> a, int n, ulong b)
    {
        int i = 0;
        for (; i < n; i++)
        {
            ulong ai = a[i];
            r[i] = ai - b;
            b = ai < b ? 1UL : 0UL;
            if (b == 0) { i++; break; }
        }
        if (i < n)
        {
            a.Slice(i, n - i).CopyTo(r.Slice(i));
        }
        return b;
    }

    /// <summary>In-place r[0..n) -= b (single limb). Returns borrow.</summary>
    public static ulong SubScalarInPlace(Span<ulong> r, int n, ulong b)
    {
        for (int i = 0; i < n && b != 0; i++)
        {
            ulong ri = r[i];
            r[i] = ri - b;
            b = ri < b ? 1UL : 0UL;
        }
        return b;
    }

    /// <summary>In-place r[0..rn) -= a[0..an), rn &gt;= an. Returns borrow.</summary>
    public static ulong SubInPlace(Span<ulong> r, int rn, ReadOnlySpan<ulong> a, int an)
    {
        Debug.Assert(rn >= an);
        ulong borrow = Sub(r, r, a, an);
        if (borrow != 0 && rn > an)
        {
            borrow = SubScalarInPlace(r.Slice(an), rn - an, borrow);
        }
        return borrow;
    }

    // ---------------------------------------------------------------
    // Vector carry propagation. Per block: lane sums s, generate mask g
    // (lane wrapped) and propagate mask p (lane all-ones). The lanes that
    // receive an incoming carry are ((g << 1 | cin) + p) ^ p, and bit W of
    // that sum is the block carry-out: a W-lane carry-lookahead adder
    // computed with one scalar add.
    // ---------------------------------------------------------------

    /// <summary>Minimum limb count before the vector carry-propagation kernels pay off.</summary>
    private const int VecAddMin = 16;

    private static ulong AddVec512(ref ulong rp, ref ulong ap, ref ulong bp, int n, ref int index, uint carry)
    {
        Vector512<ulong> lanes = Vector512.Create(1UL, 2, 4, 8, 16, 32, 64, 128);
        nuint i = (nuint)index, end = (nuint)(n & ~7);
        for (; i < end; i += 8)
        {
            Vector512<ulong> va = Vector512.LoadUnsafe(ref ap, i);
            Vector512<ulong> s = va + Vector512.LoadUnsafe(ref bp, i);
            uint g = (uint)Vector512.LessThan(s, va).ExtractMostSignificantBits();
            uint p = (uint)Vector512.Equals(s, Vector512<ulong>.AllBitsSet).ExtractMostSignificantBits();
            uint x = ((g << 1) | carry) + p;
            carry = x >> 8;
            Vector512<ulong> inc = Vector512.Equals(Vector512.Create((ulong)(x ^ p)) & lanes, lanes);
            (s - inc).StoreUnsafe(ref rp, i);
        }
        index = (int)end;
        return carry;
    }

    private static ulong SubVec512(ref ulong rp, ref ulong ap, ref ulong bp, int n, ref int index, uint borrow)
    {
        Vector512<ulong> lanes = Vector512.Create(1UL, 2, 4, 8, 16, 32, 64, 128);
        nuint i = (nuint)index, end = (nuint)(n & ~7);
        for (; i < end; i += 8)
        {
            Vector512<ulong> va = Vector512.LoadUnsafe(ref ap, i);
            Vector512<ulong> vb = Vector512.LoadUnsafe(ref bp, i);
            Vector512<ulong> d = va - vb;
            uint g = (uint)Vector512.LessThan(va, vb).ExtractMostSignificantBits();
            uint p = (uint)Vector512.Equals(d, Vector512<ulong>.Zero).ExtractMostSignificantBits();
            uint x = ((g << 1) | borrow) + p;
            borrow = x >> 8;
            Vector512<ulong> dec = Vector512.Equals(Vector512.Create((ulong)(x ^ p)) & lanes, lanes);
            (d + dec).StoreUnsafe(ref rp, i);
        }
        index = (int)end;
        return borrow;
    }

    private static ulong AddVec256(ref ulong rp, ref ulong ap, ref ulong bp, int n, ref int index, uint carry)
    {
        Vector256<ulong> lanes = Vector256.Create(1UL, 2, 4, 8);
        nuint i = (nuint)index, end = (nuint)(n & ~3);
        for (; i < end; i += 4)
        {
            Vector256<ulong> va = Vector256.LoadUnsafe(ref ap, i);
            Vector256<ulong> s = va + Vector256.LoadUnsafe(ref bp, i);
            uint g = Vector256.LessThan(s, va).ExtractMostSignificantBits();
            uint p = Vector256.Equals(s, Vector256<ulong>.AllBitsSet).ExtractMostSignificantBits();
            uint x = ((g << 1) | carry) + p;
            carry = x >> 4;
            Vector256<ulong> inc = Vector256.Equals(Vector256.Create((ulong)(x ^ p)) & lanes, lanes);
            (s - inc).StoreUnsafe(ref rp, i);
        }
        index = (int)end;
        return carry;
    }

    private static ulong SubVec256(ref ulong rp, ref ulong ap, ref ulong bp, int n, ref int index, uint borrow)
    {
        Vector256<ulong> lanes = Vector256.Create(1UL, 2, 4, 8);
        nuint i = (nuint)index, end = (nuint)(n & ~3);
        for (; i < end; i += 4)
        {
            Vector256<ulong> va = Vector256.LoadUnsafe(ref ap, i);
            Vector256<ulong> vb = Vector256.LoadUnsafe(ref bp, i);
            Vector256<ulong> d = va - vb;
            uint g = Vector256.LessThan(va, vb).ExtractMostSignificantBits();
            uint p = Vector256.Equals(d, Vector256<ulong>.Zero).ExtractMostSignificantBits();
            uint x = ((g << 1) | borrow) + p;
            borrow = x >> 4;
            Vector256<ulong> dec = Vector256.Equals(Vector256.Create((ulong)(x ^ p)) & lanes, lanes);
            (d + dec).StoreUnsafe(ref rp, i);
        }
        index = (int)end;
        return borrow;
    }

    // ---------------------------------------------------------------
    // Primitive carry helpers
    // ---------------------------------------------------------------

    /// <summary>sum = a + b + carryIn (carryIn is 0 or 1); returns carry-out (0 or 1).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static ulong AddWithCarry(ulong a, ulong b, ulong carryIn, out ulong sum)
    {
        ulong s = a + b;
        ulong c = s < a ? 1UL : 0UL;
        s += carryIn;
        c += s < carryIn ? 1UL : 0UL;
        sum = s;
        return c;
    }

    /// <summary>diff = a - b - borrowIn (borrowIn is 0 or 1); returns borrow-out (0 or 1).</summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static ulong SubWithBorrow(ulong a, ulong b, ulong borrowIn, out ulong diff)
    {
        ulong d = a - b;
        ulong c = a < b ? 1UL : 0UL;
        ulong d2 = d - borrowIn;
        c += d < borrowIn ? 1UL : 0UL;
        diff = d2;
        return c;
    }
}
