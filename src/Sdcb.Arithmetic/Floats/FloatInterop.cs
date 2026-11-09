namespace Sdcb.Arithmetic;

// Conversions between the fixed-width float types.

public readonly partial struct QuadDouble
{
    /// <summary>Exact: a double-double is a quad-double with two zero tail words.</summary>
    public static implicit operator QuadDouble(DoubleDouble value) => new(value.Hi, value.Lo);

    /// <summary>Rounds to the leading ~106 bits (hi = x0, lo ≈ x1 + x2; not correctly rounded).</summary>
    public static explicit operator DoubleDouble(QuadDouble value) =>
        double.IsFinite(value.X0) && value.X1 != 0.0 ? new DoubleDouble(value.X0, value.X1 + value.X2) : (DoubleDouble)value.X0;

    /// <summary>
    /// Splits a binary128 value into up to three doubles (exact for every value whose exponent is
    /// in double's normal range). Values beyond double's range overflow to ±∞ or lose low bits to
    /// underflow, as with <c>(double)</c>.
    /// </summary>
    public static explicit operator QuadDouble(Float128 value)
    {
        double x0 = (double)value;
        if (!double.IsFinite(x0) || x0 == 0.0)
        {
            return new QuadDouble(x0);
        }
        Float128 r = value - x0;          // exact: the rounding residual of a 113-bit value
        double x1 = (double)r;
        double x2 = (double)(r - x1);     // exact as well
        return new QuadDouble(x0, x1, x2);
    }

    /// <summary>Correctly rounded (nearest-even) to binary128 via the exact BigFloat sum of the components.</summary>
    public static explicit operator Float128(QuadDouble value) => Float128.FromBigFloat(value);
}

public readonly partial struct Float128
{
    /// <summary>Correctly rounded (nearest-even): hi and lo are exact in binary128, so hi + lo rounds once.</summary>
    public static explicit operator Float128(DoubleDouble value) =>
        value.Lo == 0.0 ? (Float128)value.Hi : (Float128)value.Hi + (Float128)value.Lo;

    /// <summary>Rounds to the leading ~106 bits (hi = nearest double, lo = nearest double of the residual).</summary>
    public static explicit operator DoubleDouble(Float128 value)
    {
        double hi = (double)value;
        if (!double.IsFinite(hi) || hi == 0.0)
        {
            return hi;
        }
        return new DoubleDouble(hi, (double)(value - hi));
    }
}
