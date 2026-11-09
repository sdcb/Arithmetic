namespace Sdcb.Arithmetic;

public readonly partial struct DoubleDouble
{
    // ---------------------------------------------------------------
    // Hard-coded constants: each (C0, C1[, C2]) is the greedy split of the 500-bit BigFloat value
    // (C0 = RN(v), C1 = RN(v - C0), ...), generated with BigFloat and checked by the unit tests.
    // ---------------------------------------------------------------

    internal const double Pi0 = 3.141592653589793;
    internal const double Pi1 = 1.2246467991473532E-16;
    internal const double Tau0 = 6.283185307179586;
    internal const double Tau1 = 2.4492935982947064E-16;
    internal const double E0 = 2.718281828459045;
    internal const double E1 = 1.4456468917292502E-16;
    internal const double Ln2_0 = 0.6931471805599453;
    internal const double Ln2_1 = 2.3190468138462996E-17;
    internal const double Ln2_2 = 5.707708438416212E-34;
    internal const double Ln10_0 = 2.302585092994046;
    internal const double Ln10_1 = -2.1707562233822494E-16;
    internal const double Ln10_2 = -9.984262454465777E-33;
    internal const double PiOver2_0 = 1.5707963267948966;
    internal const double PiOver2_1 = 6.123233995736766E-17;
    internal const double PiOver2_2 = -1.4973849048591698E-33;
    internal const double PiOver4_0 = 0.7853981633974483;
    internal const double PiOver4_1 = 3.061616997868383E-17;
    internal const double ThreePiOver4_0 = 2.356194490192345;
    internal const double ThreePiOver4_1 = 9.184850993605148E-17;
    internal const double InvPi0 = 0.3183098861837907;
    internal const double InvPi1 = -1.9678676675182486E-17;
    internal const double TwoOverPi0 = 0.6366197723675814;
    internal const double Log2E0 = 1.4426950408889634;
    internal const double Log2E1 = 2.0355273740931033E-17;
    internal const double Log10E0 = 0.4342944819032518;
    internal const double Log10E1 = 1.098319650216765E-17;
    internal const double Log10_2_0 = 0.3010299956639812;
    internal const double Log10_2_1 = -2.8037281277851704E-18;
    internal const double Log2_10_0 = 3.321928094887362;
    internal const double Log2_10_1 = 1.661617516973592E-16;
    internal const double DegToRad0 = 0.017453292519943295;
    internal const double DegToRad1 = 2.9486522708701687E-19;
    internal const double RadToDeg0 = 57.29577951308232;
    internal const double RadToDeg1 = -1.9878495670576283E-15;

    // 1/k! as (hi, lo), k = 0..32.
    internal static ReadOnlySpan<double> InvFactHi =>
    [
        1, 1, 0.5, 0.16666666666666666, 0.041666666666666664, 0.008333333333333333, 0.001388888888888889,
        0.0001984126984126984, 2.48015873015873E-05, 2.7557319223985893E-06, 2.755731922398589E-07,
        2.505210838544172E-08, 2.08767569878681E-09, 1.6059043836821613E-10, 1.1470745597729725E-11,
        7.647163731819816E-13, 4.779477332387385E-14, 2.8114572543455206E-15, 1.5619206968586225E-16,
        8.22063524662433E-18, 4.110317623312165E-19, 1.9572941063391263E-20, 8.896791392450574E-22,
        3.868170170630684E-23, 1.6117375710961184E-24, 6.446950284384474E-26, 2.4795962632247976E-27,
        9.183689863795546E-29, 3.279889237069838E-30, 1.1309962886447716E-31, 3.7699876288159054E-33,
        1.216125041553518E-34, 3.8003907548547434E-36,
    ];

    internal static ReadOnlySpan<double> InvFactLo =>
    [
        0, 0, 0, 9.25185853854297E-18, 2.3129646346357427E-18, 1.1564823173178714E-19, -5.300543954373577E-20,
        1.7209558293420705E-22, 2.1511947866775882E-23, -1.858393274046472E-22, 2.3767714622250297E-23,
        -1.448814070935912E-24, -1.20734505911326E-25, 1.2585294588752098E-26, 2.0655512752830745E-28,
        7.03872877733453E-30, 4.399205485834081E-31, 1.6508842730861433E-31, 1.1910679660273754E-32,
        2.2141894119604265E-34, 1.4412973378659527E-36, -1.3643503830087908E-36, -7.911402614872376E-38,
        -8.843177655482344E-40, -3.6846573564509766E-41, -1.9330404233703465E-42, -1.2953730964765229E-43,
        1.4303150396787322E-45, 1.5117542744029879E-46, 1.0498015412959506E-47, 2.5870347832750324E-49,
        5.586290567888806E-51, 1.7457158024652518E-52,
    ];

    // ---------------------------------------------------------------
    // Public constants
    // ---------------------------------------------------------------

    /// <summary>0.</summary>
    public static DoubleDouble Zero => default;

    /// <summary>1.</summary>
    public static DoubleDouble One => new(1.0, 0.0, normalized: true);

    /// <summary>-1.</summary>
    public static DoubleDouble NegativeOne => new(-1.0, 0.0, normalized: true);

    /// <summary>π to 106 bits.</summary>
    public static DoubleDouble Pi => new(Pi0, Pi1, normalized: true);

    /// <summary>2π to 106 bits.</summary>
    public static DoubleDouble Tau => new(Tau0, Tau1, normalized: true);

    /// <summary>EulerGamma's number e to 106 bits.</summary>
    public static DoubleDouble E => new(E0, E1, normalized: true);

    /// <summary>ln 2 to 106 bits.</summary>
    public static DoubleDouble Ln2 => new(Ln2_0, Ln2_1, normalized: true);

    /// <summary>ln 10 to 106 bits.</summary>
    public static DoubleDouble Ln10 => new(Ln10_0, Ln10_1, normalized: true);

    /// <summary>
    /// The smallest positive value, 2^-1074 (same as <see cref="double.Epsilon"/>; .NET convention,
    /// shared by <see cref="QuadDouble.Epsilon"/> and <see cref="Float128.Epsilon"/>).
    /// For the relative precision use <see cref="MachineEpsilon"/>.
    /// </summary>
    public static DoubleDouble Epsilon => new(double.Epsilon, 0.0, normalized: true);

    /// <summary>
    /// The relative precision 2^-104, the same convention as the QD library's dd_real::_eps
    /// (twice the 2^-105 spacing of 106-bit values just above 1).
    /// </summary>
    public static DoubleDouble MachineEpsilon => new(4.930380657631324E-32, 0.0, normalized: true);

    /// <summary>
    /// The largest finite value: (<see cref="double.MaxValue"/>, 2^970·(1 - 2^-53)), the largest
    /// low word that still rounds back into the high word.
    /// </summary>
    public static DoubleDouble MaxValue => new(double.MaxValue, 9.979201547673598E+291, normalized: true);

    /// <summary>The most negative finite value, -<see cref="MaxValue"/>.</summary>
    public static DoubleDouble MinValue => new(double.MinValue, -9.979201547673598E+291, normalized: true);

    /// <summary>Not a number.</summary>
    public static DoubleDouble NaN => new(double.NaN, 0.0, normalized: true);

    /// <summary>+∞.</summary>
    public static DoubleDouble PositiveInfinity => new(double.PositiveInfinity, 0.0, normalized: true);

    /// <summary>-∞.</summary>
    public static DoubleDouble NegativeInfinity => new(double.NegativeInfinity, 0.0, normalized: true);
}
