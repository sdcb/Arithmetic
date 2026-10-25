using System.Runtime.Intrinsics;
using System.Runtime.Intrinsics.Arm;
using System.Runtime.Intrinsics.X86;

namespace Sdcb.Arithmetic.Core;

/// <summary>
/// Vector tiers used by the limb kernels (AVX-512 → AVX2 → 128-bit → scalar). Every member is
/// a JIT-time constant, so tiers that are not accelerated are removed from the generated code.
/// Tiers follow the runtime's own configuration: lower them with the standard .NET knobs, e.g.
/// <c>DOTNET_PreferredVectorBitWidth=128</c>, <c>DOTNET_EnableAVX2=0</c> or
/// <c>DOTNET_EnableHWIntrinsic=0</c>. All tiers produce bit-identical results.
/// </summary>
internal static class Isa
{
    // Vector512.IsHardwareAccelerated is false on parts where .NET disables 512-bit vectors by
    // default (frequency throttling), which we honor. Each tier also checks IsHardwareAccelerated
    // so DOTNET_PreferredVectorBitWidth lowers it even when the ISA itself stays supported.
    internal static bool Hw512 => Vector512.IsHardwareAccelerated && Avx512F.IsSupported && Avx512BW.IsSupported;
    internal static bool Hw256 => Vector256.IsHardwareAccelerated && Avx2.IsSupported;
    internal static bool Hw128 => Vector128.IsHardwareAccelerated;
    /// <summary>128-bit tier with the 32x32→64 lane multiplies the NTT needs (SSE4.1 or Arm64 AdvSimd).</summary>
    internal static bool Hw128Mul => Vector128.IsHardwareAccelerated && (Sse41.IsSupported || AdvSimd.Arm64.IsSupported);
}
