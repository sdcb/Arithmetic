using System.Runtime.InteropServices;

namespace Sdcb.Arithmetic.Benchmarks;

/// <summary>Minimal P/Invoke surface over the native GMP library (baseline for benchmarks).</summary>
public static unsafe class Gmp
{
    private const string Lib = "libgmp-10";

    [StructLayout(LayoutKind.Sequential)]
    public struct MpzT
    {
        public int Alloc;
        public int Size;
        public IntPtr D;
    }

    [DllImport(Lib, EntryPoint = "__gmpz_init2", CallingConvention = CallingConvention.Cdecl)]
    public static extern void Init2(ref MpzT x, ulong bits);

    [DllImport(Lib, EntryPoint = "__gmpz_clear", CallingConvention = CallingConvention.Cdecl)]
    public static extern void Clear(ref MpzT x);

    [DllImport(Lib, EntryPoint = "__gmpz_import", CallingConvention = CallingConvention.Cdecl)]
    public static extern void Import(ref MpzT x, nuint count, int order, nuint size, int endian, nuint nails, void* op);

    [DllImport(Lib, EntryPoint = "__gmpz_add", CallingConvention = CallingConvention.Cdecl)]
    public static extern void Add(ref MpzT r, ref MpzT a, ref MpzT b);

    [DllImport(Lib, EntryPoint = "__gmpz_sub", CallingConvention = CallingConvention.Cdecl)]
    public static extern void Sub(ref MpzT r, ref MpzT a, ref MpzT b);

    [DllImport(Lib, EntryPoint = "__gmpz_mul", CallingConvention = CallingConvention.Cdecl)]
    public static extern void Mul(ref MpzT r, ref MpzT a, ref MpzT b);

    [DllImport(Lib, EntryPoint = "__gmpz_tdiv_qr", CallingConvention = CallingConvention.Cdecl)]
    public static extern void TdivQr(ref MpzT q, ref MpzT r, ref MpzT n, ref MpzT d);

    [DllImport(Lib, EntryPoint = "__gmpz_gcd", CallingConvention = CallingConvention.Cdecl)]
    public static extern void Gcd(ref MpzT r, ref MpzT a, ref MpzT b);

    [DllImport(Lib, EntryPoint = "__gmpz_powm", CallingConvention = CallingConvention.Cdecl)]
    public static extern void Powm(ref MpzT r, ref MpzT b, ref MpzT e, ref MpzT m);

    [DllImport(Lib, EntryPoint = "__gmpz_get_str", CallingConvention = CallingConvention.Cdecl)]
    public static extern IntPtr GetStr(IntPtr str, int b, ref MpzT op);

    [DllImport(Lib, EntryPoint = "__gmpz_set_str", CallingConvention = CallingConvention.Cdecl)]
    public static extern int SetStr(ref MpzT r, [MarshalAs(UnmanagedType.LPStr)] string s, int b);

    [DllImport(Lib, EntryPoint = "__gmpz_probab_prime_p", CallingConvention = CallingConvention.Cdecl)]
    public static extern int ProbabPrimeP(ref MpzT n, int reps);

    [DllImport(Lib, EntryPoint = "__gmpz_sqrt", CallingConvention = CallingConvention.Cdecl)]
    public static extern void Sqrt(ref MpzT r, ref MpzT a);

    [DllImport(Lib, EntryPoint = "__gmpz_and", CallingConvention = CallingConvention.Cdecl)]
    public static extern void And(ref MpzT r, ref MpzT a, ref MpzT b);

    [DllImport(Lib, EntryPoint = "__gmpz_xor", CallingConvention = CallingConvention.Cdecl)]
    public static extern void Xor(ref MpzT r, ref MpzT a, ref MpzT b);

    [DllImport(Lib, EntryPoint = "__gmpz_mul_2exp", CallingConvention = CallingConvention.Cdecl)]
    public static extern void Mul2Exp(ref MpzT r, ref MpzT a, ulong bits);

    [DllImport(Lib, EntryPoint = "__gmpz_tdiv_q_2exp", CallingConvention = CallingConvention.Cdecl)]
    public static extern void TdivQ2Exp(ref MpzT r, ref MpzT a, ulong bits);

    [DllImport(Lib, EntryPoint = "__gmpz_popcount", CallingConvention = CallingConvention.Cdecl)]
    public static extern ulong Popcount(ref MpzT a);

    [DllImport(Lib, EntryPoint = "__gmpz_divexact", CallingConvention = CallingConvention.Cdecl)]
    public static extern void DivExact(ref MpzT q, ref MpzT n, ref MpzT d);

    [DllImport(Lib, EntryPoint = "__gmpz_gcdext", CallingConvention = CallingConvention.Cdecl)]
    public static extern void GcdExt(ref MpzT g, ref MpzT s, ref MpzT t, ref MpzT a, ref MpzT b);

    [DllImport(Lib, EntryPoint = "__gmpz_invert", CallingConvention = CallingConvention.Cdecl)]
    public static extern int Invert(ref MpzT r, ref MpzT a, ref MpzT m);

    [DllImport(Lib, EntryPoint = "__gmpz_fac_ui", CallingConvention = CallingConvention.Cdecl)]
    public static extern void FacUi(ref MpzT r, ulong n);

    [DllImport(Lib, EntryPoint = "__gmpz_nextprime", CallingConvention = CallingConvention.Cdecl)]
    public static extern void NextPrime(ref MpzT r, ref MpzT a);

    [DllImport(Lib, EntryPoint = "__gmpz_root", CallingConvention = CallingConvention.Cdecl)]
    public static extern int Root(ref MpzT r, ref MpzT a, ulong n);

    [DllImport(Lib, EntryPoint = "__gmpz_jacobi", CallingConvention = CallingConvention.Cdecl)]
    public static extern int Jacobi(ref MpzT a, ref MpzT b);

    [DllImport(Lib, EntryPoint = "__gmpz_fib_ui", CallingConvention = CallingConvention.Cdecl)]
    public static extern void FibUi(ref MpzT r, ulong n);

    [DllImport(Lib, EntryPoint = "__gmpz_bin_uiui", CallingConvention = CallingConvention.Cdecl)]
    public static extern void BinUiUi(ref MpzT r, ulong n, ulong k);

    [DllImport(Lib, EntryPoint = "__gmpz_perfect_power_p", CallingConvention = CallingConvention.Cdecl)]
    public static extern int PerfectPowerP(ref MpzT a);

    [DllImport(Lib, EntryPoint = "__gmpz_sizeinbase", CallingConvention = CallingConvention.Cdecl)]
    public static extern nuint SizeInBase(ref MpzT a, int b);

    [DllImport(Lib, EntryPoint = "__gmp_get_memory_functions", CallingConvention = CallingConvention.Cdecl)]
    private static extern void GetMemoryFunctions(out IntPtr alloc, out IntPtr realloc, out IntPtr free);

    /// <summary>Frees a string allocated by gmp using its registered free function.</summary>
    public static void FreeStr(IntPtr str, nuint len)
    {
        GetMemoryFunctions(out _, out _, out IntPtr free);
        var freeFn = (delegate* unmanaged[Cdecl]<IntPtr, nuint, void>)free;
        freeFn(str, len);
    }

    /// <summary>Creates an initialized mpz from limb data.</summary>
    public static MpzT FromLimbs(ReadOnlySpan<ulong> limbs)
    {
        MpzT z = default;
        Init2(ref z, (ulong)(limbs.Length + 1) * 64);
        fixed (ulong* p = limbs)
        {
            Import(ref z, (nuint)limbs.Length, -1, 8, 0, 0, p);
        }
        return z;
    }
}
