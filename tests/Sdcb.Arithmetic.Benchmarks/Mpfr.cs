using System.Runtime.InteropServices;

namespace Sdcb.Arithmetic.Benchmarks;

/// <summary>
/// Minimal P/Invoke surface over native MPFR (baseline for Real benchmarks).
/// mpfr_t is treated as an opaque 32-byte block (larger than the native struct)
/// and only manipulated through mpfr functions, avoiding layout assumptions.
/// </summary>
public static unsafe class Mpfr
{
    private const string Lib = "libmpfr-6";

    [StructLayout(LayoutKind.Sequential, Size = 32)]
    public struct MpfrT
    {
        private ulong _a, _b, _c, _d;
    }

    // Rounding: 0 = RNDN, 1 = RNDZ, 2 = RNDU, 3 = RNDD, 4 = RNDA
    [DllImport(Lib, EntryPoint = "mpfr_init2", CallingConvention = CallingConvention.Cdecl)]
    public static extern void Init2(ref MpfrT x, int prec);

    [DllImport(Lib, EntryPoint = "mpfr_clear", CallingConvention = CallingConvention.Cdecl)]
    public static extern void Clear(ref MpfrT x);

    [DllImport(Lib, EntryPoint = "mpfr_set_d", CallingConvention = CallingConvention.Cdecl)]
    public static extern int SetD(ref MpfrT x, double v, int rnd);

    [DllImport(Lib, EntryPoint = "mpfr_get_d", CallingConvention = CallingConvention.Cdecl)]
    public static extern double GetD(ref MpfrT x, int rnd);

    [DllImport(Lib, EntryPoint = "mpfr_set_str", CallingConvention = CallingConvention.Cdecl)]
    public static extern int SetStr(ref MpfrT x, [MarshalAs(UnmanagedType.LPStr)] string s, int b, int rnd);

    [DllImport(Lib, EntryPoint = "mpfr_get_str", CallingConvention = CallingConvention.Cdecl)]
    public static extern IntPtr GetStr(IntPtr str, out int exp, int b, nuint n, ref MpfrT op, int rnd);

    [DllImport(Lib, EntryPoint = "mpfr_free_str", CallingConvention = CallingConvention.Cdecl)]
    public static extern void FreeStr(IntPtr str);

    [DllImport(Lib, EntryPoint = "mpfr_add", CallingConvention = CallingConvention.Cdecl)]
    public static extern int Add(ref MpfrT r, ref MpfrT a, ref MpfrT b, int rnd);

    [DllImport(Lib, EntryPoint = "mpfr_sub", CallingConvention = CallingConvention.Cdecl)]
    public static extern int Sub(ref MpfrT r, ref MpfrT a, ref MpfrT b, int rnd);

    [DllImport(Lib, EntryPoint = "mpfr_mul", CallingConvention = CallingConvention.Cdecl)]
    public static extern int Mul(ref MpfrT r, ref MpfrT a, ref MpfrT b, int rnd);

    [DllImport(Lib, EntryPoint = "mpfr_div", CallingConvention = CallingConvention.Cdecl)]
    public static extern int Div(ref MpfrT r, ref MpfrT a, ref MpfrT b, int rnd);

    [DllImport(Lib, EntryPoint = "mpfr_sqrt", CallingConvention = CallingConvention.Cdecl)]
    public static extern int Sqrt(ref MpfrT r, ref MpfrT a, int rnd);

    [DllImport(Lib, EntryPoint = "mpfr_exp", CallingConvention = CallingConvention.Cdecl)]
    public static extern int Exp(ref MpfrT r, ref MpfrT a, int rnd);

    [DllImport(Lib, EntryPoint = "mpfr_log", CallingConvention = CallingConvention.Cdecl)]
    public static extern int Log(ref MpfrT r, ref MpfrT a, int rnd);

    [DllImport(Lib, EntryPoint = "mpfr_sin", CallingConvention = CallingConvention.Cdecl)]
    public static extern int Sin(ref MpfrT r, ref MpfrT a, int rnd);

    [DllImport(Lib, EntryPoint = "mpfr_cos", CallingConvention = CallingConvention.Cdecl)]
    public static extern int Cos(ref MpfrT r, ref MpfrT a, int rnd);

    [DllImport(Lib, EntryPoint = "mpfr_atan", CallingConvention = CallingConvention.Cdecl)]
    public static extern int Atan(ref MpfrT r, ref MpfrT a, int rnd);

    [DllImport(Lib, EntryPoint = "mpfr_pow", CallingConvention = CallingConvention.Cdecl)]
    public static extern int Pow(ref MpfrT r, ref MpfrT a, ref MpfrT b, int rnd);

    [DllImport(Lib, EntryPoint = "mpfr_const_pi", CallingConvention = CallingConvention.Cdecl)]
    public static extern int ConstPi(ref MpfrT r, int rnd);

    [DllImport(Lib, EntryPoint = "mpfr_free_cache", CallingConvention = CallingConvention.Cdecl)]
    public static extern void FreeCache();

    private static readonly IntPtr s_lib = NativeLibrary.Load(Lib, typeof(Mpfr).Assembly, null);

    /// <summary>Resolves an mpfr function with signature int f(mpfr_t r, const mpfr_t x, rnd).</summary>
    public static delegate* unmanaged[Cdecl]<ref MpfrT, ref MpfrT, int, int> Unary(string name)
        => (delegate* unmanaged[Cdecl]<ref MpfrT, ref MpfrT, int, int>)NativeLibrary.GetExport(s_lib, name);

    /// <summary>Resolves an mpfr function with signature int f(mpfr_t r, const mpfr_t x, const mpfr_t y, rnd).</summary>
    public static delegate* unmanaged[Cdecl]<ref MpfrT, ref MpfrT, ref MpfrT, int, int> Binary(string name)
        => (delegate* unmanaged[Cdecl]<ref MpfrT, ref MpfrT, ref MpfrT, int, int>)NativeLibrary.GetExport(s_lib, name);

    /// <summary>Resolves an mpfr function with signature int f(mpfr_t r, long n, const mpfr_t x, rnd) (jn/yn).</summary>
    public static delegate* unmanaged[Cdecl]<ref MpfrT, int, ref MpfrT, int, int> IntFirst(string name)
        => (delegate* unmanaged[Cdecl]<ref MpfrT, int, ref MpfrT, int, int>)NativeLibrary.GetExport(s_lib, name);

    /// <summary>Resolves an mpfr function with signature int f(mpfr_t r, const mpfr_t x, long n, rnd).</summary>
    public static delegate* unmanaged[Cdecl]<ref MpfrT, ref MpfrT, int, int, int> IntSecond(string name)
        => (delegate* unmanaged[Cdecl]<ref MpfrT, ref MpfrT, int, int, int>)NativeLibrary.GetExport(s_lib, name);

    /// <summary>Resolves an mpfr constant with signature int f(mpfr_t r, rnd).</summary>
    public static delegate* unmanaged[Cdecl]<ref MpfrT, int, int> Const(string name)
        => (delegate* unmanaged[Cdecl]<ref MpfrT, int, int>)NativeLibrary.GetExport(s_lib, name);

    /// <summary>Creates an initialized mpfr from a double at the given precision.</summary>
    public static MpfrT FromDouble(double v, int prec)
    {
        MpfrT t = default;
        Init2(ref t, prec);
        SetD(ref t, v, 0);
        return t;
    }

    /// <summary>Reads the decimal digits + exponent (0.digits × 10^exp convention).</summary>
    public static (string Digits, int Exp) GetDigits(ref MpfrT v, int digits)
    {
        IntPtr p = GetStr(IntPtr.Zero, out int e, 10, (nuint)digits, ref v, 0);
        string s = Marshal.PtrToStringAnsi(p)!;
        FreeStr(p);
        return (s, e);
    }
}
