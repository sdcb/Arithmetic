using System.Diagnostics;

namespace Sdcb.Arithmetic.Core;

/// <summary>
/// Montgomery multiplication context for an odd modulus, used by modular
/// exponentiation (Montgomery, 1985).
/// </summary>
internal sealed class MontgomeryContext
{
    private readonly ulong[] _m;      // modulus limbs
    private readonly int _n;          // limb count
    private readonly ulong _n0inv;    // -m^-1 mod 2^64
    private readonly ulong[] _rr;     // B^(2n) mod m (to enter Montgomery domain)
    private readonly ulong[] _rmod;   // B^n mod m (Montgomery form of 1)
    private readonly ulong[]? _minv;  // -m^-1 mod B^n, for REDC by multiplication
    private readonly ulong[]? _q;     // REDC scratch: n limbs
    private readonly ulong[]? _p;     // REDC scratch: 2n limbs

    /// <summary>Modulus size (limbs) from which REDC uses two multiplications instead of n AddMulScalar rows.</summary>
    private static readonly int s_redcMulMin = int.TryParse(Environment.GetEnvironmentVariable("SDCB_ARITHMETIC_REDC_MUL"), out int v) ? v : Isa.Hw256 ? 24 : 160;

    public int N => _n;
    public ReadOnlySpan<ulong> RMod => _rmod;

    public MontgomeryContext(ReadOnlySpan<ulong> modulus)
    {
        Debug.Assert(modulus.Length >= 1 && (modulus[0] & 1) != 0 && modulus[^1] != 0);
        _n = modulus.Length;
        _m = modulus.ToArray();
        _n0inv = unchecked((ulong)(-(long)Inv64(_m[0])));

        // rr = B^(2n) mod m ; rmod = B^n mod m
        _rr = PowerOfBMod(2 * _n);
        _rmod = PowerOfBMod(_n);
        if (_n >= s_redcMulMin)
        {
            _minv = NegInverseModB(_m, _n);
            _q = new ulong[_n];
            _p = new ulong[2 * _n];
        }
    }

    /// <summary>-m^-1 mod B^n by Newton-Hensel lifting (x' = x (2 - m x)).</summary>
    private static ulong[] NegInverseModB(ulong[] m, int n)
    {
        BigInt mm = BigInt.FromMagnitude((ulong[])m.Clone(), n, 1);
        BigInt x = new(Inv64(m[0]));
        for (long bits = 64; bits < 64L * n; bits *= 2)
        {
            long nb = Math.Min(2 * bits, 64L * n);
            BigInt mask = (BigInt.One << (int)nb) - BigInt.One;
            x = (x * (new BigInt(2) - ((mm & mask) * x & mask))) & mask;
        }
        BigInt neg = (BigInt.One << (64 * n)) - x;
        ulong[] r = new ulong[n];
        neg.Limbs.CopyTo(r);
        return r;
    }

    /// <summary>Computes B^k mod m via full division.</summary>
    private ulong[] PowerOfBMod(int k)
    {
        int an = k + 1;
        ulong[] a = new ulong[an];
        a[k] = 1;
        ulong[] q = new ulong[an - _n + 1];
        ulong[] r = new ulong[_n];
        Nat.DivRem(q, r, a, an, _m, _n);
        return r;
    }

    /// <summary>Returns x^-1 mod 2^64 for odd x (Newton iteration).</summary>
    internal static ulong Inv64(ulong x)
    {
        Debug.Assert((x & 1) != 0);
        ulong inv = x; // 3 correct bits
        inv *= 2 - x * inv; // 6
        inv *= 2 - x * inv; // 12
        inv *= 2 - x * inv; // 24
        inv *= 2 - x * inv; // 48
        inv *= 2 - x * inv; // 96 >= 64
        Debug.Assert(x * inv == 1);
        return inv;
    }

    /// <summary>
    /// REDC of the 2n-limb value t &lt; m B^n (t is destroyed); result (n limbs, &lt; m) into res.
    /// Each row's carry is parked in the limb it just cleared and all of them
    /// are added in one pass at the end, instead of propagating every row's carry.
    /// </summary>
    private void Reduce(Span<ulong> t, Span<ulong> res)
    {
        int n = _n;
        Debug.Assert(t.Length >= 2 * n);
        if (_minv is not null)
        {
            ReduceByMul(t, res);
            return;
        }
        for (int i = 0; i < n; i++)
        {
            ulong u = t[i] * _n0inv;
            t[i] = Nat.AddMulScalar(t.Slice(i, n), _m, n, u);
        }
        ulong cy = Nat.Add(res, t.Slice(n, n), t, n);
        if (cy != 0 || Nat.Compare(res, _m, n) >= 0)
        {
            Nat.Sub(res, res, _m, n);
        }
    }

    /// <summary>
    /// REDC with two products: q = t_lo * (-m^-1) mod B^n, then (t + q m) / B^n. The low halves
    /// of t and q m sum to exactly 0 or B^n, and only to 0 when t_lo is zero, so the carry into
    /// the high half is just (t_lo != 0).
    /// </summary>
    private void ReduceByMul(Span<ulong> t, Span<ulong> res)
    {
        int n = _n;
        Span<ulong> q = _q, p = _p;
        if (Isa.Hw256 && n <= Nat.MulKaratsubaThreshold)
        {
            Nat.MulLowVec256(q, t, _minv, n);
        }
        else
        {
            Nat.Mul(p, t, n, _minv, n);
            p.Slice(0, n).CopyTo(q);
        }
        Nat.Mul(p, q, n, _m, n);
        ulong cy = Nat.Add(res, t.Slice(n, n), p.Slice(n, n), n);
        if (!Nat.IsZero(t, n))
        {
            cy += Nat.AddScalarInPlace(res, n, 1);
        }
        if (cy != 0 || Nat.Compare(res, _m, n) >= 0)
        {
            Nat.Sub(res, res, _m, n);
        }
    }

    /// <summary>res = MontMul(x, y) = x*y*B^-n mod m. res may alias x or y.</summary>
    public void MontMul(Span<ulong> res, ReadOnlySpan<ulong> x, ReadOnlySpan<ulong> y, Span<ulong> t)
    {
        int n = _n;
        Nat.Mul(t, x, n, y, n);
        Reduce(t, res);
    }

    /// <summary>res = MontSqr(x). res may alias x.</summary>
    public void MontSqr(Span<ulong> res, ReadOnlySpan<ulong> x, Span<ulong> t)
    {
        int n = _n;
        Nat.Square(t, x, n);
        Reduce(t, res);
    }

    /// <summary>Enters the Montgomery domain: res = x * B^n mod m (x must be reduced mod m, n limbs zero-padded).</summary>
    public void ToMont(Span<ulong> res, ReadOnlySpan<ulong> x, Span<ulong> t) => MontMul(res, x, _rr, t);

    /// <summary>Leaves the Montgomery domain: res = x * B^-n mod m.</summary>
    public void FromMont(Span<ulong> res, ReadOnlySpan<ulong> x, Span<ulong> t)
    {
        int n = _n;
        t.Slice(0, 2 * n + 1).Clear();
        x.Slice(0, n).CopyTo(t);
        Reduce(t, res);
    }
}
