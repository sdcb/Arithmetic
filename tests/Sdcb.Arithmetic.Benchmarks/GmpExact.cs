using System.Runtime.InteropServices;

namespace Sdcb.Arithmetic.Benchmarks;

/// <summary>
/// Exact cross-check of BigInt against native GMP for the mpz API: random signs and sizes
/// (including subquadratic-algorithm sizes), results compared exactly. `-- gmp-exact [filter]`.
/// </summary>
internal static unsafe class GmpExact
{
    private const string L = "libgmp-10";
    private const CallingConvention C = CallingConvention.Cdecl;

    [DllImport(L, EntryPoint = "__gmpz_init", CallingConvention = C)] private static extern void Init(ref Gmp.MpzT x);
    [DllImport(L, EntryPoint = "__gmpz_tdiv_q", CallingConvention = C)] private static extern void TdivQ(ref Gmp.MpzT q, ref Gmp.MpzT n, ref Gmp.MpzT d);
    [DllImport(L, EntryPoint = "__gmpz_tdiv_r", CallingConvention = C)] private static extern void TdivR(ref Gmp.MpzT r, ref Gmp.MpzT n, ref Gmp.MpzT d);
    [DllImport(L, EntryPoint = "__gmpz_fdiv_q", CallingConvention = C)] private static extern void FdivQ(ref Gmp.MpzT q, ref Gmp.MpzT n, ref Gmp.MpzT d);
    [DllImport(L, EntryPoint = "__gmpz_fdiv_r", CallingConvention = C)] private static extern void FdivR(ref Gmp.MpzT r, ref Gmp.MpzT n, ref Gmp.MpzT d);
    [DllImport(L, EntryPoint = "__gmpz_cdiv_q", CallingConvention = C)] private static extern void CdivQ(ref Gmp.MpzT q, ref Gmp.MpzT n, ref Gmp.MpzT d);
    [DllImport(L, EntryPoint = "__gmpz_cdiv_r", CallingConvention = C)] private static extern void CdivR(ref Gmp.MpzT r, ref Gmp.MpzT n, ref Gmp.MpzT d);
    [DllImport(L, EntryPoint = "__gmpz_lcm", CallingConvention = C)] private static extern void Lcm(ref Gmp.MpzT r, ref Gmp.MpzT a, ref Gmp.MpzT b);
    [DllImport(L, EntryPoint = "__gmpz_sqrtrem", CallingConvention = C)] private static extern void SqrtRem(ref Gmp.MpzT r, ref Gmp.MpzT rem, ref Gmp.MpzT a);
    [DllImport(L, EntryPoint = "__gmpz_rootrem", CallingConvention = C)] private static extern void RootRem(ref Gmp.MpzT r, ref Gmp.MpzT rem, ref Gmp.MpzT a, uint k);
    [DllImport(L, EntryPoint = "__gmpz_perfect_square_p", CallingConvention = C)] private static extern int PerfectSquareP(ref Gmp.MpzT a);
    [DllImport(L, EntryPoint = "__gmpz_ior", CallingConvention = C)] private static extern void Ior(ref Gmp.MpzT r, ref Gmp.MpzT a, ref Gmp.MpzT b);
    [DllImport(L, EntryPoint = "__gmpz_com", CallingConvention = C)] private static extern void Com(ref Gmp.MpzT r, ref Gmp.MpzT a);
    [DllImport(L, EntryPoint = "__gmpz_fdiv_q_2exp", CallingConvention = C)] private static extern void FdivQ2Exp(ref Gmp.MpzT r, ref Gmp.MpzT a, uint b);
    [DllImport(L, EntryPoint = "__gmpz_cdiv_q_2exp", CallingConvention = C)] private static extern void CdivQ2Exp(ref Gmp.MpzT r, ref Gmp.MpzT a, uint b);
    [DllImport(L, EntryPoint = "__gmpz_tdiv_r_2exp", CallingConvention = C)] private static extern void TdivR2Exp(ref Gmp.MpzT r, ref Gmp.MpzT a, uint b);
    [DllImport(L, EntryPoint = "__gmpz_fdiv_r_2exp", CallingConvention = C)] private static extern void FdivR2Exp(ref Gmp.MpzT r, ref Gmp.MpzT a, uint b);
    [DllImport(L, EntryPoint = "__gmpz_cdiv_r_2exp", CallingConvention = C)] private static extern void CdivR2Exp(ref Gmp.MpzT r, ref Gmp.MpzT a, uint b);
    [DllImport(L, EntryPoint = "__gmpz_hamdist", CallingConvention = C)] private static extern uint Hamdist(ref Gmp.MpzT a, ref Gmp.MpzT b);
    [DllImport(L, EntryPoint = "__gmpz_scan0", CallingConvention = C)] private static extern uint Scan0(ref Gmp.MpzT a, uint start);
    [DllImport(L, EntryPoint = "__gmpz_scan1", CallingConvention = C)] private static extern uint Scan1(ref Gmp.MpzT a, uint start);
    [DllImport(L, EntryPoint = "__gmpz_jacobi", CallingConvention = C)] private static extern int Kronecker(ref Gmp.MpzT a, ref Gmp.MpzT b);
    [DllImport(L, EntryPoint = "__gmpz_prevprime", CallingConvention = C)] private static extern int PrevPrime(ref Gmp.MpzT r, ref Gmp.MpzT a);
    [DllImport(L, EntryPoint = "__gmpz_2fac_ui", CallingConvention = C)] private static extern void DoubleFac(ref Gmp.MpzT r, uint n);
    [DllImport(L, EntryPoint = "__gmpz_mfac_uiui", CallingConvention = C)] private static extern void MFac(ref Gmp.MpzT r, uint n, uint m);
    [DllImport(L, EntryPoint = "__gmpz_primorial_ui", CallingConvention = C)] private static extern void Primorial(ref Gmp.MpzT r, uint n);
    [DllImport(L, EntryPoint = "__gmpz_bin_ui", CallingConvention = C)] private static extern void BinUi(ref Gmp.MpzT r, ref Gmp.MpzT n, uint k);
    [DllImport(L, EntryPoint = "__gmpz_lucnum_ui", CallingConvention = C)] private static extern void LucNum(ref Gmp.MpzT r, uint n);
    [DllImport(L, EntryPoint = "__gmpz_remove", CallingConvention = C)] private static extern uint Remove(ref Gmp.MpzT r, ref Gmp.MpzT a, ref Gmp.MpzT f);
    [DllImport(L, EntryPoint = "__gmpz_congruent_p", CallingConvention = C)] private static extern int CongruentP(ref Gmp.MpzT a, ref Gmp.MpzT c, ref Gmp.MpzT d);
    [DllImport(L, EntryPoint = "__gmpz_divisible_p", CallingConvention = C)] private static extern int DivisibleP(ref Gmp.MpzT a, ref Gmp.MpzT d);
    [DllImport(L, EntryPoint = "__gmpz_get_d_2exp", CallingConvention = C)] private static extern double GetD2Exp(out int exp, ref Gmp.MpzT a);
    [DllImport(L, EntryPoint = "__gmpz_cmp_d", CallingConvention = C)] private static extern int CmpD(ref Gmp.MpzT a, double d);

    private static int s_checks, s_failures;
    private static string? s_filter;

    private static Gmp.MpzT Z(BigInt v)
    {
        Gmp.MpzT z = default;
        Init(ref z);
        Gmp.SetStr(ref z, v.ToString(16), 16);
        return z;
    }

    private static BigInt B(ref Gmp.MpzT z)
    {
        IntPtr s = Gmp.GetStr(IntPtr.Zero, 16, ref z);
        string hex = Marshal.PtrToStringAnsi(s)!;
        Gmp.FreeStr(s, (nuint)(hex.Length + 1));
        return BigInt.Parse(hex, 16);
    }

    private static void Eq(string op, BigInt ours, ref Gmp.MpzT z, string args)
    {
        s_checks++;
        BigInt g = B(ref z);
        if (g != ours)
        {
            Fail(op, args, $"gmp={g.ToString(16)} ours={ours.ToString(16)}");
        }
    }

    private static void EqL(string op, long ours, long gmp, string args)
    {
        s_checks++;
        if (ours != gmp)
        {
            Fail(op, args, $"gmp={gmp} ours={ours}");
        }
    }

    private static void Fail(string op, string args, string diff)
    {
        s_failures++;
        if (s_failures <= 60)
        {
            Console.WriteLine($"  MISMATCH {op} {args}: {diff}");
        }
    }

    private static bool Want(string name) => s_filter is null || name.Contains(s_filter, StringComparison.OrdinalIgnoreCase);

    private static BigInt Rand(Random r, int maxLimbs, bool allowNeg = true)
    {
        int limbs = r.Next(0, maxLimbs + 1);
        byte[] b = new byte[limbs * 8 + r.Next(0, 8)];
        r.NextBytes(b);
        if (b.Length > 16 && r.Next(4) == 0)
        {
            int s = r.Next(b.Length / 2), l = r.Next(b.Length - s);
            b.AsSpan(s, l).Fill(r.Next(2) == 0 ? (byte)0 : (byte)0xFF);
        }
        return BigInt.FromMagnitudeBytes(b, allowNeg && r.Next(2) == 0 ? -1 : 1);
    }

    private static int[] SizeMix => [1, 2, 3, 6, 12, 40, 120];

    public static void Run(string[] args)
    {
        s_filter = args.Length > 1 ? args[1] : null;
        var rng = new Random(2024);

        void Group(string name, int reps, Action<Random, int> body)
        {
            if (!Want(name))
            {
                return;
            }
            int before = s_failures;
            foreach (int size in SizeMix)
            {
                for (int i = 0; i < reps; i++)
                {
                    body(rng, size);
                }
            }
            Console.WriteLine($"{name,-14} {(s_failures == before ? "ok" : $"FAILED ({s_failures - before})")}");
        }

        Group("arith", 30, (r, n) =>
        {
            BigInt a = Rand(r, n), b = Rand(r, n);
            Gmp.MpzT za = Z(a), zb = Z(b), zr = default, zq = default;
            Init(ref zr); Init(ref zq);
            string s = $"a={a.ToString(16)} b={b.ToString(16)}";
            Gmp.Add(ref zr, ref za, ref zb); Eq("add", a + b, ref zr, s);
            Gmp.Sub(ref zr, ref za, ref zb); Eq("sub", a - b, ref zr, s);
            Gmp.Mul(ref zr, ref za, ref zb); Eq("mul", a * b, ref zr, s);
            Gmp.And(ref zr, ref za, ref zb); Eq("and", a & b, ref zr, s);
            Ior(ref zr, ref za, ref zb); Eq("ior", a | b, ref zr, s);
            Gmp.Xor(ref zr, ref za, ref zb); Eq("xor", a ^ b, ref zr, s);
            Com(ref zr, ref za); Eq("com", ~a, ref zr, s);
            if (!b.IsZero)
            {
                TdivQ(ref zq, ref za, ref zb); Eq("tdiv_q", a / b, ref zq, s);
                TdivR(ref zr, ref za, ref zb); Eq("tdiv_r", a % b, ref zr, s);
                BigInt.DivRem(a, b, out BigInt fq, out BigInt fr, IntegerRounding.Floor);
                FdivQ(ref zq, ref za, ref zb); Eq("fdiv_q", fq, ref zq, s);
                FdivR(ref zr, ref za, ref zb); Eq("fdiv_r", fr, ref zr, s);
                BigInt.DivRem(a, b, out BigInt cq, out BigInt cr, IntegerRounding.Ceiling);
                CdivQ(ref zq, ref za, ref zb); Eq("cdiv_q", cq, ref zq, s);
                CdivR(ref zr, ref za, ref zb); Eq("cdiv_r", cr, ref zr, s);
                BigInt prod = a * b;
                Gmp.MpzT zp = Z(prod);
                Gmp.DivExact(ref zq, ref zp, ref zb); Eq("divexact", BigInt.DivExact(prod, b), ref zq, s);
                EqL("divisible_p", BigInt.IsDivisible(a, b) ? 1 : 0, DivisibleP(ref za, ref zb) != 0 ? 1 : 0, s);
                BigInt c = Rand(r, n);
                Gmp.MpzT zc = Z(c);
                EqL("congruent_p", BigInt.IsCongruent(a, c, b) ? 1 : 0, CongruentP(ref za, ref zc, ref zb) != 0 ? 1 : 0, s);
                Gmp.Clear(ref zp); Gmp.Clear(ref zc);
            }
            Gmp.Clear(ref za); Gmp.Clear(ref zb); Gmp.Clear(ref zr); Gmp.Clear(ref zq);
        });

        Group("gcd", 20, (r, n) =>
        {
            BigInt f = Rand(r, Math.Max(1, n / 4), allowNeg: false) | BigInt.One;
            BigInt a = Rand(r, n) * (r.Next(2) == 0 ? f : BigInt.One), b = Rand(r, n) * f;
            Gmp.MpzT za = Z(a), zb = Z(b), zg = default, zs = default, zt = default;
            Init(ref zg); Init(ref zs); Init(ref zt);
            string s = $"a={a.ToString(16)} b={b.ToString(16)}";
            Gmp.Gcd(ref zg, ref za, ref zb); Eq("gcd", BigInt.Gcd(a, b), ref zg, s);
            Lcm(ref zg, ref za, ref zb); Eq("lcm", BigInt.Lcm(a, b), ref zg, s);
            BigInt g = BigInt.ExtendedGcd(a, b, out BigInt x, out BigInt y);
            Gmp.GcdExt(ref zg, ref zs, ref zt, ref za, ref zb);
            Eq("gcdext.g", g, ref zg, s);
            Eq("gcdext.s", x, ref zs, s);
            Eq("gcdext.t", y, ref zt, s);
            if (!b.IsZero && !b.IsOne && !b.IsMinusOne)
            {
                int gi = Gmp.Invert(ref zg, ref za, ref zb);
                bool ok = BigInt.TryModInverse(a, b, out BigInt inv);
                EqL("invert.ok", ok ? 1 : 0, gi != 0 ? 1 : 0, s);
                if (ok && gi != 0)
                {
                    Eq("invert", inv, ref zg, s);
                }
            }
            Gmp.Clear(ref za); Gmp.Clear(ref zb); Gmp.Clear(ref zg); Gmp.Clear(ref zs); Gmp.Clear(ref zt);
        });

        Group("powm", 6, (r, n) =>
        {
            int mn = Math.Min(n, 40);
            BigInt b = Rand(r, mn), e = Rand(r, Math.Min(n, 6), allowNeg: false), m = Rand(r, mn, allowNeg: false);
            if (m.IsZero)
            {
                return;
            }
            if (r.Next(2) == 0)
            {
                m |= BigInt.One;
            }
            Gmp.MpzT zb = Z(b), ze = Z(e), zm = Z(m), zr = default;
            Init(ref zr);
            Gmp.Powm(ref zr, ref zb, ref ze, ref zm);
            Eq("powm", BigInt.ModPow(b, e, m), ref zr, $"b={b.ToString(16)} e={e.ToString(16)} m={m.ToString(16)}");
            Gmp.Clear(ref zb); Gmp.Clear(ref ze); Gmp.Clear(ref zm); Gmp.Clear(ref zr);
        });

        Group("roots", 12, (r, n) =>
        {
            BigInt a = Rand(r, n, allowNeg: false);
            if (r.Next(3) == 0)
            {
                a = a * a;
            }
            Gmp.MpzT za = Z(a), zr = default, zm = default;
            Init(ref zr); Init(ref zm);
            string s = $"a={a.ToString(16)}";
            BigInt sq = BigInt.Sqrt(a, out BigInt rem);
            SqrtRem(ref zr, ref zm, ref za);
            Eq("sqrtrem.s", sq, ref zr, s);
            Eq("sqrtrem.r", rem, ref zm, s);
            EqL("perfect_square_p", BigInt.IsPerfectSquare(a) ? 1 : 0, PerfectSquareP(ref za) != 0 ? 1 : 0, s);
            EqL("perfect_power_p", BigInt.IsPerfectPower(a) ? 1 : 0, Gmp.PerfectPowerP(ref za) != 0 ? 1 : 0, s);
            uint k = (uint)r.Next(2, 9);
            BigInt rt = BigInt.Root(a, k, out BigInt rr);
            RootRem(ref zr, ref zm, ref za, k);
            Eq("rootrem.r", rt, ref zr, s + $" k={k}");
            Eq("rootrem.rem", rr, ref zm, s + $" k={k}");
            Gmp.Clear(ref za); Gmp.Clear(ref zr); Gmp.Clear(ref zm);
        });

        Group("bits", 20, (r, n) =>
        {
            BigInt a = Rand(r, n), b = Rand(r, n, allowNeg: false);
            Gmp.MpzT za = Z(a), zb = Z(b), zr = default;
            Init(ref zr);
            uint sh = (uint)r.Next(0, 64 * n + 70);
            string s = $"a={a.ToString(16)} sh={sh}";
            Gmp.Mul2Exp(ref zr, ref za, sh); Eq("mul_2exp", a << (int)sh, ref zr, s);
            Gmp.TdivQ2Exp(ref zr, ref za, sh); Eq("tdiv_q_2exp", a.ShiftRight(sh, IntegerRounding.Truncate), ref zr, s);
            FdivQ2Exp(ref zr, ref za, sh); Eq("fdiv_q_2exp", a.ShiftRight(sh, IntegerRounding.Floor), ref zr, s);
            CdivQ2Exp(ref zr, ref za, sh); Eq("cdiv_q_2exp", a.ShiftRight(sh, IntegerRounding.Ceiling), ref zr, s);
            TdivR2Exp(ref zr, ref za, sh); Eq("tdiv_r_2exp", a.RemainderPow2(sh, IntegerRounding.Truncate), ref zr, s);
            FdivR2Exp(ref zr, ref za, sh); Eq("fdiv_r_2exp", a.ModPow2(sh), ref zr, s);
            CdivR2Exp(ref zr, ref za, sh); Eq("cdiv_r_2exp", a.RemainderPow2(sh, IntegerRounding.Ceiling), ref zr, s);
            Gmp.MpzT zba = Z(a.Abs());
            EqL("popcount", a.Abs().PopCount(), Gmp.Popcount(ref zba) is var pc ? (long)(uint)pc : 0, s);
            EqL("hamdist", BigInt.HammingDistance(a.Abs(), b), Hamdist(ref zba, ref zb), s);
            uint st = (uint)r.Next(0, 64 * n + 8);
            long s0 = a.IndexOfClearBit(st), s1 = a.IndexOfSetBit(st);
            EqL("scan0", s0, Scan0(ref za, st) is var g0 && g0 == uint.MaxValue ? -1 : g0, s + $" st={st}");
            EqL("scan1", s1, Scan1(ref za, st) is var g1 && g1 == uint.MaxValue ? -1 : g1, s + $" st={st}");
            int e;
            double dm = GetD2Exp(out e, ref za);
            double om = a.ToDouble(out long oe);
            EqL("get_d_2exp", BitConverter.DoubleToInt64Bits(om), BitConverter.DoubleToInt64Bits(dm), s);
            EqL("get_d_2exp.e", oe, e, s);
            double dv = (r.NextDouble() - 0.5) * Math.Pow(2, r.Next(0, 64 * n + 2));
            EqL("cmp_d", Math.Sign(a.CompareTo(dv)), Math.Sign(CmpD(ref za, dv)), s + $" d={dv:R}");
            Gmp.Clear(ref za); Gmp.Clear(ref zb); Gmp.Clear(ref zr); Gmp.Clear(ref zba);
        });

        Group("symbols", 20, (r, n) =>
        {
            BigInt a = Rand(r, n), b = Rand(r, Math.Max(1, n), allowNeg: true);
            Gmp.MpzT za = Z(a), zb = Z(b);
            string s = $"a={a.ToString(16)} b={b.ToString(16)}";
            EqL("kronecker", BigInt.Kronecker(a, b), Kronecker(ref za, ref zb), s);
            BigInt odd = b.Abs() | BigInt.One;
            Gmp.MpzT zo = Z(odd);
            EqL("jacobi", BigInt.Jacobi(a, odd), Gmp.Jacobi(ref za, ref zo), s);
            Gmp.Clear(ref za); Gmp.Clear(ref zb); Gmp.Clear(ref zo);
        });

        Group("primes", 4, (r, n) =>
        {
            int pn = Math.Min(n, 12);
            BigInt a = Rand(r, pn, allowNeg: false);
            Gmp.MpzT za = Z(a), zr = default;
            Init(ref zr);
            string s = $"a={a.ToString(16)}";
            EqL("probab_prime_p", a.IsProbablePrime() ? 1 : 0, Gmp.ProbabPrimeP(ref za, 25) > 0 ? 1 : 0, s);
            Gmp.NextPrime(ref zr, ref za); Eq("nextprime", a.NextPrime(), ref zr, s);
            int gp = PrevPrime(ref zr, ref za);
            bool op = a.TryPrevPrime(out BigInt pp);
            EqL("prevprime.ok", op ? 1 : 0, gp != 0 ? 1 : 0, s);
            if (op && gp != 0)
            {
                Eq("prevprime", pp, ref zr, s);
            }
            Gmp.Clear(ref za); Gmp.Clear(ref zr);
        });

        Group("combinatorics", 3, (r, n) =>
        {
            Gmp.MpzT zr = default;
            Init(ref zr);
            uint k = (uint)r.Next(0, 40 * n + 10);
            Gmp.FacUi(ref zr, k); Eq("fac_ui", BigInt.Factorial(k), ref zr, $"n={k}");
            DoubleFac(ref zr, k); Eq("2fac_ui", BigInt.DoubleFactorial(k), ref zr, $"n={k}");
            uint m = (uint)r.Next(1, 7);
            MFac(ref zr, k, m); Eq("mfac_uiui", BigInt.MultiFactorial(k, m), ref zr, $"n={k} m={m}");
            Primorial(ref zr, k); Eq("primorial_ui", BigInt.Primorial(k), ref zr, $"n={k}");
            uint kk = (uint)r.Next(0, (int)k + 3);
            Gmp.BinUiUi(ref zr, k, kk); Eq("bin_uiui", BigInt.Binomial((ulong)k, kk), ref zr, $"n={k} k={kk}");
            BigInt nb = Rand(r, 2);
            Gmp.MpzT zn = Z(nb);
            uint kb = (uint)r.Next(0, 30);
            BinUi(ref zr, ref zn, kb); Eq("bin_ui", BigInt.Binomial(nb, kb), ref zr, $"n={nb} k={kb}");
            Gmp.FibUi(ref zr, k); Eq("fib_ui", BigInt.Fibonacci(k), ref zr, $"n={k}");
            LucNum(ref zr, k); Eq("lucnum_ui", BigInt.Lucas(k), ref zr, $"n={k}");
            Gmp.Clear(ref zr); Gmp.Clear(ref zn);
        });

        Group("convert", 10, (r, n) =>
        {
            BigInt a = Rand(r, n);
            Gmp.MpzT za = Z(a), zr = default;
            Init(ref zr);
            int bs = r.Next(2, 63);
            IntPtr p = Gmp.GetStr(IntPtr.Zero, bs, ref za);
            string gs = Marshal.PtrToStringAnsi(p)!;
            Gmp.FreeStr(p, (nuint)(gs.Length + 1));
            string os = a.ToString(bs);
            s_checks++;
            if (gs != os)
            {
                Fail("get_str", $"base={bs} a={a.ToString(16)}", $"gmp={gs} ours={os}");
            }
            Eq("set_str", BigInt.Parse(gs, bs), ref za, $"base={bs}");
            EqL("sizeinbase", a.CountDigits(bs), (long)Gmp.SizeInBase(ref za, bs), $"base={bs} a={a.ToString(16)}");
            BigInt f = new(r.Next(2, 12));
            Gmp.MpzT zf = Z(f);
            BigInt rm = BigInt.RemoveFactor(a * BigInt.Pow(f, (ulong)r.Next(0, 6)), f, out ulong cnt);
            Gmp.MpzT zx = Z(a * BigInt.Pow(f, (ulong)0));
            _ = zx;
            Gmp.MpzT zy = Z(rm * BigInt.Pow(f, cnt));
            uint gc = a.IsZero ? 0 : Remove(ref zr, ref zy, ref zf);
            if (!a.IsZero)
            {
                Eq("remove", rm, ref zr, $"a={a.ToString(16)} f={f}");
                EqL("remove.count", (long)cnt, gc, $"a={a.ToString(16)} f={f}");
            }
            Gmp.Clear(ref za); Gmp.Clear(ref zr); Gmp.Clear(ref zf); Gmp.Clear(ref zx); Gmp.Clear(ref zy);
        });

        // A few NTT/HGCD-sized operands.
        if (Want("large"))
        {
            int before = s_failures;
            foreach (int n in (int[])[900, 3000, 12000])
            {
                BigInt a = Rand(new Random(n), n), b = Rand(new Random(n + 1), n * 2 / 3 + 1);
                Gmp.MpzT za = Z(a), zb = Z(b), zr = default, zs = default, zt = default;
                Init(ref zr); Init(ref zs); Init(ref zt);
                string s = $"n={n}";
                Gmp.Mul(ref zr, ref za, ref zb); Eq("mul(large)", a * b, ref zr, s);
                if (!b.IsZero)
                {
                    TdivQ(ref zr, ref za, ref zb); Eq("tdiv_q(large)", a / b, ref zr, s);
                }
                Gmp.Gcd(ref zr, ref za, ref zb); Eq("gcd(large)", BigInt.Gcd(a, b), ref zr, s);
                BigInt g = BigInt.ExtendedGcd(a, b, out BigInt x, out BigInt y);
                Gmp.GcdExt(ref zr, ref zs, ref zt, ref za, ref zb);
                Eq("gcdext(large).s", x, ref zs, s);
                Eq("gcdext(large).t", y, ref zt, s);
                Gmp.MpzT zabs = Z(a.Abs());
                Gmp.Sqrt(ref zr, ref zabs); Eq("sqrt(large)", BigInt.Sqrt(a.Abs()), ref zr, s);
                Gmp.Clear(ref zabs);
                Gmp.Clear(ref za); Gmp.Clear(ref zb); Gmp.Clear(ref zr); Gmp.Clear(ref zs); Gmp.Clear(ref zt);
            }
            Console.WriteLine($"{"large",-14} {(s_failures == before ? "ok" : $"FAILED ({s_failures - before})")}");
        }

        Console.WriteLine($"gmp-exact: {s_checks} checks, {s_failures} failures");
        if (s_failures > 0)
        {
            Environment.ExitCode = 1;
        }
    }
}
