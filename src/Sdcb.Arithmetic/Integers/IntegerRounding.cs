namespace Sdcb.Arithmetic;

/// <summary>How an integer quotient (or right shift) is rounded when the division is inexact.</summary>
public enum IntegerRounding
{
    /// <summary>Round toward zero; the remainder has the sign of the dividend (C#'s <c>/</c> and <c>%</c>).</summary>
    Truncate = 0,
    /// <summary>Round toward negative infinity; the remainder has the sign of the divisor (Python's <c>//</c> and <c>%</c>).</summary>
    Floor = 1,
    /// <summary>Round toward positive infinity; the remainder has the opposite sign of the divisor.</summary>
    Ceiling = 2,
}

/// <summary>Result of a probabilistic primality test.</summary>
public enum Primality
{
    /// <summary>Definitely composite.</summary>
    Composite = 0,
    /// <summary>Probably prime (passed BPSW and the requested Miller-Rabin rounds; no counterexample is known).</summary>
    ProbablePrime = 1,
    /// <summary>Proven prime (small enough to be checked deterministically).</summary>
    Prime = 2,
}
