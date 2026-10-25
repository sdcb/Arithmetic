using System.Runtime.InteropServices;

namespace Sdcb.Arithmetic.Benchmarks;

/// <summary>
/// Bit-exact cross-check of Real against native MPFR: every function, all five rounding
/// modes, several precisions. Results are compared as (mantissa, exponent) pairs via
/// mpfr_get_z_2exp, plus NaN / signed infinity / signed zero. `-- exact [filter]`.
/// </summary>
internal static unsafe class MpfrExact
{
    private const string Lib = "libmpfr-6";

    [DllImport(Lib, EntryPoint = "mpfr_set_z_2exp", CallingConvention = CallingConvention.Cdecl)]
    private static extern int SetZ2Exp(ref Mpfr.MpfrT r, ref Gmp.MpzT z, int e, int rnd);

    [DllImport(Lib, EntryPoint = "mpfr_get_z_2exp", CallingConvention = CallingConvention.Cdecl)]
    private static extern int GetZ2Exp(ref Gmp.MpzT z, ref Mpfr.MpfrT x);

    [DllImport(Lib, EntryPoint = "mpfr_nan_p", CallingConvention = CallingConvention.Cdecl)]
    private static extern int NanP(ref Mpfr.MpfrT x);

    [DllImport(Lib, EntryPoint = "mpfr_inf_p", CallingConvention = CallingConvention.Cdecl)]
    private static extern int InfP(ref Mpfr.MpfrT x);

    [DllImport(Lib, EntryPoint = "mpfr_zero_p", CallingConvention = CallingConvention.Cdecl)]
    private static extern int ZeroP(ref Mpfr.MpfrT x);

    [DllImport(Lib, EntryPoint = "mpfr_signbit", CallingConvention = CallingConvention.Cdecl)]
    private static extern int SignBit(ref Mpfr.MpfrT x);

    [DllImport(Lib, EntryPoint = "mpfr_set_nan", CallingConvention = CallingConvention.Cdecl)]
    private static extern void SetNan(ref Mpfr.MpfrT x);

    [DllImport(Lib, EntryPoint = "mpfr_set_inf", CallingConvention = CallingConvention.Cdecl)]
    private static extern void SetInf(ref Mpfr.MpfrT x, int sign);

    [DllImport(Lib, EntryPoint = "mpfr_set_zero", CallingConvention = CallingConvention.Cdecl)]
    private static extern void SetZero(ref Mpfr.MpfrT x, int sign);

    [DllImport(Lib, EntryPoint = "mpfr_set_prec", CallingConvention = CallingConvention.Cdecl)]
    private static extern void SetPrec(ref Mpfr.MpfrT x, int prec);

    private static readonly IntPtr s_lib = NativeLibrary.Load(Lib, typeof(MpfrExact).Assembly, null);

    private delegate Real Unary(Real x, int prec, RealRounding rnd);
    private delegate Real Binary(Real x, Real y, int prec, RealRounding rnd);

    /// <summary>Input generator: returns a random Real of the given precision in the function's domain.</summary>
    private delegate Real Gen(Random rng, int prec);

    // string.GetHashCode() is randomized per process, which made the cross-check flaky.
    private static int StableSeed(string name)
    {
        int h = 17;
        foreach (char c in name)
        {
            h = unchecked(h * 31 + c);
        }
        return h & 0x7FFF;
    }

    private static int s_checks, s_failures;
    private static string? s_filter;

    public static void Run(string[] args)
    {
        // `exact [filter] [big]`: "big" adds multi-limb precisions that exercise the
        // pooled-scratch and long-field rounding paths (slow for the special functions).
        bool bigPrecs = args.Skip(1).Contains("big", StringComparer.OrdinalIgnoreCase);
        s_filter = args.Length > 1 && !args[1].Equals("big", StringComparison.OrdinalIgnoreCase) ? args[1] : null;
        int[] precs = bigPrecs ? [2, 7, 24, 53, 64, 113, 200, 640, 4097, 9000] : [2, 7, 24, 53, 64, 113, 200, 640];

        // Domains: generic (wide exponents), moderate (|x| up to ~2^8), unit interval, positive.
        Gen any = (r, p) => RandomReal(r, p, -40, 40, allowNeg: true);
        Gen mod = (r, p) => RandomReal(r, p, -20, 8, allowNeg: true);
        Gen pos = (r, p) => RandomReal(r, p, -40, 40, allowNeg: false);
        Gen unit = (r, p) => RandomReal(r, p, -30, 0, allowNeg: true);
        Gen big = (r, p) => RandomReal(r, p, -4, 12, allowNeg: true);
        Gen posMod = (r, p) => RandomReal(r, p, -20, 8, allowNeg: false);
        Gen geOne = (r, p) => Real.Add(Real.One, RandomReal(r, p, -30, 20, allowNeg: false), p);

        foreach (var (name, ours, native, gen) in UnaryCases(any, mod, pos, unit, big, posMod, geOne))
        {
            if (s_filter is not null && !name.Contains(s_filter, StringComparison.OrdinalIgnoreCase))
            {
                continue;
            }
            var f = (delegate* unmanaged[Cdecl]<Mpfr.MpfrT*, Mpfr.MpfrT*, int, int>)NativeLibrary.GetExport(s_lib, native);
            int before = s_failures;
            var rng = new Random(StableSeed(name));
            foreach (int prec in precs)
            {
                int count = prec > 300 ? 6 : 24;
                for (int i = 0; i < count; i++)
                {
                    Real x = gen(rng, rng.Next(2) == 0 ? prec : Math.Max(2, prec / 2 + rng.Next(prec)));
                    for (int rnd = 0; rnd < 5; rnd++)
                    {
                        Real o;
                        try
                        {
                            o = ours(x, prec, (RealRounding)rnd);
                        }
                        catch (Exception ex)
                        {
                            Fail(name, prec, rnd, $"x={Show(x)}", $"threw {ex.GetType().Name}: {ex.Message}");
                            continue;
                        }
                        Mpfr.MpfrT* mx = ToMpfr(x), mr = NewMpfr(prec);
                        f(mr, mx, rnd);
                        Check(name, prec, rnd, $"x={Show(x)}", o, mr);
                        Free(mx);
                        Free(mr);
                    }
                }
            }
            Console.WriteLine($"{name,-12} {(s_failures == before ? "ok" : $"FAILED ({s_failures - before})")}");
        }

        foreach (var (name, ours, native, gx, gy) in BinaryCases(any, mod, pos, posMod))
        {
            if (s_filter is not null && !name.Contains(s_filter, StringComparison.OrdinalIgnoreCase))
            {
                continue;
            }
            var f = (delegate* unmanaged[Cdecl]<Mpfr.MpfrT*, Mpfr.MpfrT*, Mpfr.MpfrT*, int, int>)NativeLibrary.GetExport(s_lib, native);
            int before = s_failures;
            var rng = new Random(StableSeed(name));
            foreach (int prec in precs)
            {
                int count = prec > 300 ? 6 : 24;
                for (int i = 0; i < count; i++)
                {
                    Real x = gx(rng, prec), y = gy(rng, rng.Next(2) == 0 ? prec : Math.Max(2, prec / 2 + rng.Next(prec)));
                    for (int rnd = 0; rnd < 5; rnd++)
                    {
                        Real o;
                        try
                        {
                            o = ours(x, y, prec, (RealRounding)rnd);
                        }
                        catch (Exception ex)
                        {
                            Fail(name, prec, rnd, $"x={Show(x)} y={Show(y)}", $"threw {ex.GetType().Name}: {ex.Message}");
                            continue;
                        }
                        Mpfr.MpfrT* mx = ToMpfr(x), my = ToMpfr(y), mr = NewMpfr(prec);
                        f(mr, mx, my, rnd);
                        Check(name, prec, rnd, $"x={Show(x)} y={Show(y)}", o, mr);
                        Free(mx);
                        Free(my);
                        Free(mr);
                    }
                }
            }
            Console.WriteLine($"{name,-12} {(s_failures == before ? "ok" : $"FAILED ({s_failures - before})")}");
        }
        Misc(precs);
        Console.WriteLine($"exact: {s_checks} checks, {s_failures} failures");
        if (s_failures > 0)
        {
            Environment.ExitCode = 1;
        }
    }

    private static bool Want(string name) => s_filter is null || name.Contains(s_filter, StringComparison.OrdinalIgnoreCase);

    /// <summary>Functions with integer arguments, extra outputs or vector inputs.</summary>
    private static void Misc(int[] precs)
    {
        var rng = new Random(77);
        var jn = (delegate* unmanaged[Cdecl]<Mpfr.MpfrT*, int, Mpfr.MpfrT*, int, int>)NativeLibrary.GetExport(s_lib, "mpfr_jn");
        var yn = (delegate* unmanaged[Cdecl]<Mpfr.MpfrT*, int, Mpfr.MpfrT*, int, int>)NativeLibrary.GetExport(s_lib, "mpfr_yn");
        var compound = (delegate* unmanaged[Cdecl]<Mpfr.MpfrT*, Mpfr.MpfrT*, int, int, int>)NativeLibrary.GetExport(s_lib, "mpfr_compound_si");
        var rootnsi = (delegate* unmanaged[Cdecl]<Mpfr.MpfrT*, Mpfr.MpfrT*, int, int, int>)NativeLibrary.GetExport(s_lib, "mpfr_rootn_si");
        var lgamma = (delegate* unmanaged[Cdecl]<Mpfr.MpfrT*, int*, Mpfr.MpfrT*, int, int>)NativeLibrary.GetExport(s_lib, "mpfr_lgamma");
        var remquo = (delegate* unmanaged[Cdecl]<Mpfr.MpfrT*, int*, Mpfr.MpfrT*, Mpfr.MpfrT*, int, int>)NativeLibrary.GetExport(s_lib, "mpfr_remquo");
        var fmodquo = (delegate* unmanaged[Cdecl]<Mpfr.MpfrT*, int*, Mpfr.MpfrT*, Mpfr.MpfrT*, int, int>)NativeLibrary.GetExport(s_lib, "mpfr_fmodquo");
        var modf = (delegate* unmanaged[Cdecl]<Mpfr.MpfrT*, Mpfr.MpfrT*, Mpfr.MpfrT*, int, int>)NativeLibrary.GetExport(s_lib, "mpfr_modf");
        var frexp = (delegate* unmanaged[Cdecl]<int*, Mpfr.MpfrT*, Mpfr.MpfrT*, int, int>)NativeLibrary.GetExport(s_lib, "mpfr_frexp");
        var euler = (delegate* unmanaged[Cdecl]<Mpfr.MpfrT*, int, int>)NativeLibrary.GetExport(s_lib, "mpfr_const_euler");
        var mpfrpi = (delegate* unmanaged[Cdecl]<Mpfr.MpfrT*, int, int>)NativeLibrary.GetExport(s_lib, "mpfr_const_pi");
        var mpfrlog2 = (delegate* unmanaged[Cdecl]<Mpfr.MpfrT*, int, int>)NativeLibrary.GetExport(s_lib, "mpfr_const_log2");
        var catalan = (delegate* unmanaged[Cdecl]<Mpfr.MpfrT*, int, int>)NativeLibrary.GetExport(s_lib, "mpfr_const_catalan");
        var facui = (delegate* unmanaged[Cdecl]<Mpfr.MpfrT*, uint, int, int>)NativeLibrary.GetExport(s_lib, "mpfr_fac_ui");
        var zetaui = (delegate* unmanaged[Cdecl]<Mpfr.MpfrT*, uint, int, int>)NativeLibrary.GetExport(s_lib, "mpfr_zeta_ui");
        var logui = (delegate* unmanaged[Cdecl]<Mpfr.MpfrT*, uint, int, int>)NativeLibrary.GetExport(s_lib, "mpfr_log_ui");
        var sum = (delegate* unmanaged[Cdecl]<Mpfr.MpfrT*, Mpfr.MpfrT**, uint, int, int>)NativeLibrary.GetExport(s_lib, "mpfr_sum");
        var dot = (delegate* unmanaged[Cdecl]<Mpfr.MpfrT*, Mpfr.MpfrT**, Mpfr.MpfrT**, uint, int, int>)NativeLibrary.GetExport(s_lib, "mpfr_dot");

        void Run(string name, Action<int, int> body)
        {
            if (!Want(name))
            {
                return;
            }
            int before = s_failures;
            foreach (int prec in precs)
            {
                for (int rnd = 0; rnd < 5; rnd++)
                {
                    body(prec, rnd);
                }
            }
            Console.WriteLine($"{name,-12} {(s_failures == before ? "ok" : $"FAILED ({s_failures - before})")}");
        }

        void Cmp(string name, int prec, int rnd, string args, Real ours, Mpfr.MpfrT* mr)
        {
            Check(name, prec, rnd, args, ours, mr);
        }

        Run("jn", (prec, rnd) =>
        {
            for (int i = 0; i < 3; i++)
            {
                int n = rng.Next(-6, 25);
                Real x = RandomReal(rng, prec, -10, 5, allowNeg: true);
                Mpfr.MpfrT* mx = ToMpfr(x), mr = NewMpfr(prec);
                jn(mr, n, mx, rnd);
                Cmp("jn", prec, rnd, $"n={n} x={Show(x)}", Real.BesselJ(n, x, prec, (RealRounding)rnd), mr);
                Free(mx); Free(mr);
            }
        });
        Run("yn", (prec, rnd) =>
        {
            for (int i = 0; i < 3; i++)
            {
                int n = rng.Next(-6, 25);
                Real x = RandomReal(rng, prec, -10, 5, allowNeg: false);
                Mpfr.MpfrT* mx = ToMpfr(x), mr = NewMpfr(prec);
                yn(mr, n, mx, rnd);
                Cmp("yn", prec, rnd, $"n={n} x={Show(x)}", Real.BesselY(n, x, prec, (RealRounding)rnd), mr);
                Free(mx); Free(mr);
            }
        });
        Run("compound", (prec, rnd) =>
        {
            for (int i = 0; i < 3; i++)
            {
                int n = rng.Next(-40, 40);
                Real x = Real.Add(Real.MinusOne, RandomReal(rng, prec, -20, 6, allowNeg: false), prec + 64);
                Mpfr.MpfrT* mx = ToMpfr(x), mr = NewMpfr(prec);
                compound(mr, mx, n, rnd);
                Cmp("compound", prec, rnd, $"n={n} x={Show(x)}", Real.Compound(x, n, prec, (RealRounding)rnd), mr);
                Free(mx); Free(mr);
            }
        });
        Run("rootn_si", (prec, rnd) =>
        {
            for (int i = 0; i < 3; i++)
            {
                int k = rng.Next(-9, 10);
                Real x = RandomReal(rng, prec, -30, 30, allowNeg: true);
                Mpfr.MpfrT* mx = ToMpfr(x), mr = NewMpfr(prec);
                rootnsi(mr, mx, k, rnd);
                Cmp("rootn_si", prec, rnd, $"k={k} x={Show(x)}", Real.RootN(x, (long)k, prec, (RealRounding)rnd), mr);
                Free(mx); Free(mr);
            }
        });
        Run("lgamma", (prec, rnd) =>
        {
            for (int i = 0; i < 3; i++)
            {
                Real x = RandomReal(rng, prec, -10, 4, allowNeg: true);
                Mpfr.MpfrT* mx = ToMpfr(x), mr = NewMpfr(prec);
                int sign;
                lgamma(mr, &sign, mx, rnd);
                Real o = Real.LogGamma(x, out int os, prec, (RealRounding)rnd);
                Cmp("lgamma", prec, rnd, $"x={Show(x)}", o, mr);
                if (os != sign && !o.IsNaN)
                {
                    Fail("lgamma", prec, rnd, $"x={Show(x)}", $"sign {os} != {sign}");
                }
                Free(mx); Free(mr);
            }
        });
        Run("remquo", (prec, rnd) =>
        {
            for (int i = 0; i < 3; i++)
            {
                Real x = RandomReal(rng, prec, -20, 40, allowNeg: true), y = RandomReal(rng, prec, -20, 20, allowNeg: true);
                Mpfr.MpfrT* mx = ToMpfr(x), my = ToMpfr(y), mr = NewMpfr(prec);
                int q;
                remquo(mr, &q, mx, my, rnd);
                Real o = Real.Ieee754Remainder(x, y, out long oq, prec, (RealRounding)rnd);
                Cmp("remquo", prec, rnd, $"x={Show(x)} y={Show(y)}", o, mr);
                // mpfr returns at least the low bits of a long (32 bits here): compare mod 2^30.
                if (((oq - q) & ((1L << 30) - 1)) != 0)
                {
                    Fail("remquo", prec, rnd, $"x={Show(x)} y={Show(y)}", $"quotient {oq} vs {q}");
                }
                fmodquo(mr, &q, mx, my, rnd);
                o = Real.TruncatedRemainder(x, y, out oq, prec, (RealRounding)rnd);
                Cmp("fmodquo", prec, rnd, $"x={Show(x)} y={Show(y)}", o, mr);
                if (((oq - q) & ((1L << 30) - 1)) != 0)
                {
                    Fail("fmodquo", prec, rnd, $"x={Show(x)} y={Show(y)}", $"quotient {oq} vs {q}");
                }
                Free(mx); Free(my); Free(mr);
            }
        });
        Run("modf/frexp", (prec, rnd) =>
        {
            for (int i = 0; i < 3; i++)
            {
                Real x = RandomReal(rng, prec + rng.Next(40), -20, 40, allowNeg: true);
                Mpfr.MpfrT* mx = ToMpfr(x), mi = NewMpfr(prec), mf = NewMpfr(prec);
                modf(mi, mf, mx, rnd);
                Real oi = Real.Truncate(x, out Real of, prec, (RealRounding)rnd);
                Cmp("modf.int", prec, rnd, $"x={Show(x)}", oi, mi);
                Cmp("modf.frac", prec, rnd, $"x={Show(x)}", of, mf);
                int e;
                frexp(&e, mi, mx, rnd);
                Real om = Real.Decompose(x, out long oe, prec, (RealRounding)rnd);
                Cmp("frexp", prec, rnd, $"x={Show(x)}", om, mi);
                if (oe != e)
                {
                    Fail("frexp", prec, rnd, $"x={Show(x)}", $"exp {oe} vs {e}");
                }
                Free(mx); Free(mi); Free(mf);
            }
        });
        Run("constants", (prec, rnd) =>
        {
            Mpfr.MpfrT* mr = NewMpfr(prec);
            euler(mr, rnd);
            Cmp("euler", prec, rnd, "", Real.EulerGamma(prec, (RealRounding)rnd), mr);
            catalan(mr, rnd);
            Cmp("catalan", prec, rnd, "", Real.CatalanConstant(prec, (RealRounding)rnd), mr);
            uint n = (uint)rng.Next(0, 300);
            facui(mr, n, rnd);
            Cmp("fac_ui", prec, rnd, $"n={n}", Real.Factorial(n, prec, (RealRounding)rnd), mr);
            zetaui(mr, n % 60, rnd);
            Cmp("zeta_ui", prec, rnd, $"n={n % 60}", Real.Zeta(n % 60, prec, (RealRounding)rnd), mr);
            uint m = (uint)rng.Next(0, int.MaxValue);
            logui(mr, m, rnd);
            Cmp("log_ui", prec, rnd, $"n={m}", Real.Log(m, prec, (RealRounding)rnd), mr);
            mpfrpi(mr, rnd);
            Cmp("const_pi", prec, rnd, "", Real.Pi(prec, (RealRounding)rnd), mr);
            mpfrlog2(mr, rnd);
            Cmp("const_log2", prec, rnd, "", Real.Ln2(prec, (RealRounding)rnd), mr);
            Free(mr);
        });
        Run("sum/dot", (prec, rnd) =>
        {
            int len = rng.Next(1, 12);
            var xs = new Real[len];
            var ys = new Real[len];
            var px = stackalloc Mpfr.MpfrT*[len];
            var py = stackalloc Mpfr.MpfrT*[len];
            for (int i = 0; i < len; i++)
            {
                xs[i] = RandomReal(rng, rng.Next(2, prec + 20), -60, 60, allowNeg: true);
                ys[i] = RandomReal(rng, rng.Next(2, prec + 20), -60, 60, allowNeg: true);
                if (i > 0 && rng.Next(4) == 0)
                {
                    xs[i] = -xs[i - 1]; // exact cancellation
                }
                px[i] = ToMpfr(xs[i]);
                py[i] = ToMpfr(ys[i]);
            }
            Mpfr.MpfrT* mr = NewMpfr(prec);
            sum(mr, px, (uint)len, rnd);
            Cmp("sum", prec, rnd, $"n={len}", Real.Sum(xs, prec, (RealRounding)rnd), mr);
            dot(mr, px, py, (uint)len, rnd);
            Cmp("dot", prec, rnd, $"n={len}", Real.Dot(xs, ys, prec, (RealRounding)rnd), mr);
            Free(mr);
            for (int i = 0; i < len; i++)
            {
                Free(px[i]);
                Free(py[i]);
            }
        });
    }

    private static IEnumerable<(string, Unary, string, Gen)> UnaryCases(Gen any, Gen mod, Gen pos, Gen unit, Gen big, Gen posMod, Gen geOne)
    {
        yield return ("sqrt", Real.Sqrt, "mpfr_sqrt", pos);
        yield return ("rec_sqrt", Real.ReciprocalSqrt, "mpfr_rec_sqrt", pos);
        yield return ("cbrt", Real.Cbrt, "mpfr_cbrt", any);
        yield return ("sqr", Real.Square, "mpfr_sqr", any);
        yield return ("exp", Real.Exp, "mpfr_exp", mod);
        yield return ("exp2", Real.Exp2, "mpfr_exp2", mod);
        yield return ("exp10", Real.Exp10, "mpfr_exp10", mod);
        yield return ("expm1", Real.ExpM1, "mpfr_expm1", mod);
        yield return ("log", Real.Log, "mpfr_log", pos);
        yield return ("log2", Real.Log2, "mpfr_log2", pos);
        yield return ("log10", Real.Log10, "mpfr_log10", pos);
        yield return ("log1p", Real.LogP1, "mpfr_log1p", pos);
        yield return ("sin", Real.Sin, "mpfr_sin", big);
        yield return ("cos", Real.Cos, "mpfr_cos", big);
        yield return ("tan", Real.Tan, "mpfr_tan", big);
        yield return ("asin", Real.Asin, "mpfr_asin", unit);
        yield return ("acos", Real.Acos, "mpfr_acos", unit);
        yield return ("atan", Real.Atan, "mpfr_atan", any);
        yield return ("sinh", Real.Sinh, "mpfr_sinh", mod);
        yield return ("cosh", Real.Cosh, "mpfr_cosh", mod);
        yield return ("tanh", Real.Tanh, "mpfr_tanh", mod);
        yield return ("asinh", Real.Asinh, "mpfr_asinh", any);
        yield return ("acosh", Real.Acosh, "mpfr_acosh", geOne);
        yield return ("atanh", Real.Atanh, "mpfr_atanh", unit);
        foreach (var c in ExtraUnaryCases(any, mod, pos, unit, big, posMod, geOne))
        {
            yield return c;
        }
    }

    private static IEnumerable<(string, Binary, string, Gen, Gen)> BinaryCases(Gen any, Gen mod, Gen pos, Gen posMod)
    {
        yield return ("add", Real.Add, "mpfr_add", any, any);
        yield return ("sub", Real.Subtract, "mpfr_sub", any, any);
        yield return ("mul", Real.Multiply, "mpfr_mul", any, any);
        yield return ("div", Real.Divide, "mpfr_div", any, any);
        yield return ("pow", Real.Pow, "mpfr_pow", pos, mod);
        yield return ("atan2", Real.Atan2, "mpfr_atan2", any, any);
        yield return ("hypot", Real.Hypot, "mpfr_hypot", any, any);
        yield return ("fmod", Real.TruncatedRemainder, "mpfr_fmod", any, any);
        yield return ("remainder", Real.Ieee754Remainder, "mpfr_remainder", any, any);
        yield return ("dim", Real.PositiveDifference, "mpfr_dim", any, any);
        foreach (var c in ExtraBinaryCases(any, mod, pos, posMod))
        {
            yield return c;
        }
    }

    // Functions added beyond the original surface register here.
    private static IEnumerable<(string, Unary, string, Gen)> ExtraUnaryCases(Gen any, Gen mod, Gen pos, Gen unit, Gen big, Gen posMod, Gen geOne)
    {
        Gen quarter = (r, p) => r.Next(3) == 0 ? Real.FromBigInt(new BigInt(r.Next(-64, 64)), 16).ShiftRight(2) : mod(r, p);
        Gen small = (r, p) => RandomReal(r, p, -20, 4, allowNeg: true);
        Gen gt = (r, p) => Real.Add(Real.MinusOne, RandomReal(r, p, -30, 20, allowNeg: false), p + 64);
        yield return ("sec", Real.Sec, "mpfr_sec", big);
        yield return ("csc", Real.Csc, "mpfr_csc", big);
        yield return ("cot", Real.Cot, "mpfr_cot", big);
        yield return ("sech", Real.Sech, "mpfr_sech", mod);
        yield return ("csch", Real.Csch, "mpfr_csch", mod);
        yield return ("coth", Real.Coth, "mpfr_coth", mod);
        yield return ("sinpi", Real.SinPi, "mpfr_sinpi", quarter);
        yield return ("cospi", Real.CosPi, "mpfr_cospi", quarter);
        yield return ("tanpi", Real.TanPi, "mpfr_tanpi", quarter);
        yield return ("exp2m1", Real.Exp2M1, "mpfr_exp2m1", quarter);
        yield return ("exp10m1", Real.Exp10M1, "mpfr_exp10m1", quarter);
        yield return ("log2p1", Real.Log2P1, "mpfr_log2p1", gt);
        yield return ("log10p1", Real.Log10P1, "mpfr_log10p1", gt);
        yield return ("gamma", Real.Gamma, "mpfr_gamma", small);
        yield return ("lngamma", Real.LogGamma, "mpfr_lngamma", small);
        yield return ("digamma", Real.Digamma, "mpfr_digamma", small);
        yield return ("erf", Real.Erf, "mpfr_erf", small);
        yield return ("erfc", Real.Erfc, "mpfr_erfc", small);
        yield return ("zeta", Real.Zeta, "mpfr_zeta", small);
        yield return ("eint", Real.ExponentialIntegral, "mpfr_eint", small);
        yield return ("li2", Real.Dilogarithm, "mpfr_li2", any);
        yield return ("j0", (x, p, r) => Real.BesselJ(0, x, p, r), "mpfr_j0", small);
        yield return ("j1", (x, p, r) => Real.BesselJ(1, x, p, r), "mpfr_j1", small);
        yield return ("y0", (x, p, r) => Real.BesselY(0, x, p, r), "mpfr_y0", (r, p) => RandomReal(r, p, -20, 4, allowNeg: false));
        yield return ("y1", (x, p, r) => Real.BesselY(1, x, p, r), "mpfr_y1", (r, p) => RandomReal(r, p, -20, 4, allowNeg: false));
        yield return ("ai", Real.AiryAi, "mpfr_ai", (r, p) => RandomReal(r, p, -20, 3, allowNeg: true));
    }

    private static IEnumerable<(string, Binary, string, Gen, Gen)> ExtraBinaryCases(Gen any, Gen mod, Gen pos, Gen posMod)
    {
        Gen small = (r, p) => RandomReal(r, p, -12, 4, allowNeg: true);
        yield return ("agm", Real.ArithmeticGeometricMean, "mpfr_agm", pos, pos);
        yield return ("beta", Real.Beta, "mpfr_beta", small, small);
        yield return ("powr", Real.Powr, "mpfr_powr", pos, mod);
    }

    private static void Check(string name, int prec, int rnd, string args, Real ours, Mpfr.MpfrT* mr)
    {
        string? diff = Diff(ours, ref *mr, prec);
        s_checks++;
        if (diff is not null)
        {
            Fail(name, prec, rnd, args, diff);
        }
    }

    private static string Show(Real x) => x.IsNumber
        ? $"{(x.SignBit ? "-" : "")}0x{x.Mantissa.ToString(16)}p{x.Exponent}[{x.Precision}]"
        : x.ToString();

    private static Mpfr.MpfrT* NewMpfr(int prec)
    {
        var p = (Mpfr.MpfrT*)NativeMemory.AllocZeroed((nuint)sizeof(Mpfr.MpfrT));
        Mpfr.Init2(ref *p, prec);
        return p;
    }

    private static void Free(Mpfr.MpfrT* p)
    {
        Mpfr.Clear(ref *p);
        NativeMemory.Free(p);
    }

    private static void Fail(string name, int prec, int rnd, string args, string diff)
    {
        s_failures++;
        if (s_failures <= 60)
        {
            Console.WriteLine($"  MISMATCH {name} prec={prec} rnd={rnd} {args}: {diff}");
        }
    }

    private static string? Diff(Real ours, ref Mpfr.MpfrT m, int prec)
    {
        if (NanP(ref m) != 0)
        {
            return ours.IsNaN ? null : $"mpfr=NaN ours={ours}";
        }
        if (ours.IsNaN)
        {
            return "ours=NaN";
        }
        bool mneg = SignBit(ref m) != 0;
        if (InfP(ref m) != 0)
        {
            return ours.IsInfinity && ours.SignBit == mneg ? null : $"mpfr={(mneg ? "-" : "+")}inf ours={ours}";
        }
        if (ZeroP(ref m) != 0)
        {
            return ours.IsZero && ours.SignBit == mneg ? null : $"mpfr={(mneg ? "-" : "+")}0 ours={ours}";
        }
        if (!ours.IsNumber)
        {
            return $"mpfr finite, ours={ours}";
        }
        Gmp.MpzT z = default;
        Gmp.Init2(ref z, (ulong)prec + 64);
        int e = GetZ2Exp(ref z, ref m);
        IntPtr s = Gmp.GetStr(IntPtr.Zero, 16, ref z);
        string hex = Marshal.PtrToStringAnsi(s)!;
        Gmp.FreeStr(s, (nuint)(hex.Length + 1));
        Gmp.Clear(ref z);
        BigInt mz = BigInt.Parse(hex, 16);
        BigInt om = ours.SignBit ? -ours.Mantissa : ours.Mantissa;
        if (ours.Precision != prec)
        {
            return $"precision {ours.Precision} != {prec}";
        }
        if (om == mz && ours.Exponent == e)
        {
            return null;
        }
        return $"mpfr={mz.ToString(16)}*2^{e} ours={om.ToString(16)}*2^{ours.Exponent}";
    }

    internal static Mpfr.MpfrT* ToMpfr(Real x)
    {
        Mpfr.MpfrT* p = NewMpfr(Math.Max(2, x.Precision));
        ref Mpfr.MpfrT m = ref *p;
        if (x.IsNaN)
        {
            SetNan(ref m);
        }
        else if (x.IsInfinity)
        {
            SetInf(ref m, x.SignBit ? -1 : 1);
        }
        else if (x.IsZero)
        {
            SetZero(ref m, x.SignBit ? -1 : 1);
        }
        else
        {
            BigInt man = x.SignBit ? -x.Mantissa : x.Mantissa;
            Gmp.MpzT z = default;
            Gmp.Init2(ref z, (ulong)x.Precision + 64);
            Gmp.SetStr(ref z, man.ToString(16), 16);
            int t = SetZ2Exp(ref m, ref z, (int)x.Exponent, 0);
            Gmp.Clear(ref z);
            if (t != 0)
            {
                throw new InvalidOperationException("inexact input transfer");
            }
        }
        return p;
    }

    internal static Real RandomReal(Random rng, int prec, int emin, int emax, bool allowNeg)
    {
        byte[] b = new byte[(prec + 7) / 8 + 1];
        rng.NextBytes(b);
        BigInt m = BigInt.FromMagnitudeBytes(b) | BigInt.One;
        long bits = m.GetBitLength();
        if (bits > prec)
        {
            m >>= (int)(bits - prec);
        }
        m = m.SetBit(prec - 1);
        // Occasionally exercise special mantissas (runs of ones / exact powers).
        int kind = rng.Next(12);
        if (kind == 0)
        {
            m = BigInt.One << (prec - 1);
        }
        else if (kind == 1)
        {
            m = (BigInt.One << prec) - 1;
        }
        int exp = rng.Next(emin, emax + 1);
        Real r = Real.FromBigInt(m, prec).ShiftLeft(exp - prec);
        return allowNeg && rng.Next(2) == 0 ? -r : r;
    }
}
