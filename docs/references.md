# References

`Sdcb.Arithmetic` is an independent C# implementation of published algorithms. It does not
contain GMP or MPFR source code. Native GMP and MPFR are used only as black-box references:
the `gmp-exact`/`exact` harnesses compare results, and `perf` compares timings.

Some comments name a GMP or MPFR function (for example "same role as `mpz_powm`") to say which
operation a method corresponds to. The algorithms themselves come from the literature below.

## General

- D. E. Knuth, *The Art of Computer Programming, Vol. 2: Seminumerical Algorithms*, 3rd ed., Addison-Wesley, 1997.
- R. P. Brent and P. Zimmermann, *Modern Computer Arithmetic*, Cambridge University Press, 2010.
- A. J. Menezes, P. C. van Oorschot and S. A. Vanstone, *Handbook of Applied Cryptography*, CRC Press, 1996.

## Integer kernels

| Area | Algorithm | Reference |
|---|---|---|
| Multiplication | Karatsuba | A. Karatsuba and Yu. Ofman, "Multiplication of many-digital numbers by automatic computers", *Doklady Akad. Nauk SSSR* 145, 1962. |
| Multiplication | Toom-3 and its interpolation sequence | A. L. Toom, 1963; S. A. Cook, PhD thesis, Harvard, 1966; M. Bodrato, "Towards optimal Toom-Cook multiplication for univariate and multivariate polynomials in characteristic 2 and 0", WAIFI 2007. |
| Multiplication | Number-theoretic transform over three primes + CRT | J. M. Pollard, "The fast Fourier transform in a finite field", *Math. Comp.* 25(114), 1971; H. L. Garner, "The residue number system", *IRE Trans. Electronic Computers* EC-8(2), 1959. |
| Division | Schoolbook division (Algorithm D) | Knuth, TAOCP Vol. 2, §4.3.1. |
| Division | Division by invariant integers (`invert_limb`-style reciprocal) | N. Möller and T. Granlund, "Improved division by invariant integers", *IEEE Trans. Computers* 60(2), 2011. |
| Division | Recursive divide-and-conquer division | C. Burnikel and J. Ziegler, "Fast recursive division", MPI-I-98-1-022, 1998. |
| Division | Newton reciprocal, approximate quotient | Brent and Zimmermann, *Modern Computer Arithmetic*, §1.4 and §3.4. |
| Exact division | Right-to-left exact division | T. Jebelean, "An algorithm for exact division", *J. Symbolic Computation* 15, 1993. |
| Exact division | Hensel (2-adic) division, divide-and-conquer | Brent and Zimmermann, *Modern Computer Arithmetic*, §1.4. |
| GCD | Lehmer's algorithm | D. H. Lehmer, "Euclid's algorithm for large numbers", *Amer. Math. Monthly* 45, 1938. |
| GCD | Double-digit Lehmer step | T. Jebelean, "A double-digit Lehmer-Euclid algorithm for finding the GCD of long integers", *J. Symbolic Computation* 19, 1995. |
| GCD | Subquadratic half-GCD | N. Möller, "On Schönhage's algorithm and subquadratic integer gcd computation", *Math. Comp.* 77, 2008. |
| Square root | Karatsuba square root | P. Zimmermann, "Karatsuba Square Root", INRIA Research Report RR-3805, 1999. |
| Modular arithmetic | Montgomery reduction | P. L. Montgomery, "Modular multiplication without trial division", *Math. Comp.* 44, 1985. |
| Primality | Miller-Rabin | G. L. Miller, *J. Comput. System Sci.* 13, 1976; M. O. Rabin, *J. Number Theory* 12, 1980. |
| Primality | BPSW (strong Lucas test) | R. Baillie and S. S. Wagstaff Jr., "Lucas pseudoprimes", *Math. Comp.* 35, 1980; C. Pomerance, J. L. Selfridge and S. S. Wagstaff Jr., "The pseudoprimes to 25·10⁹", *Math. Comp.* 35, 1980. |
| Combinatorics | Prime-swing factorial | P. Luschny, "Fast factorial functions" (online notes). |
| Combinatorics | Binomial via prime factorization | E. E. Kummer, *J. Reine Angew. Math.* 44, 1852 (Kummer's theorem). |
| Conversion | Divide-and-conquer radix conversion | Brent and Zimmermann, *Modern Computer Arithmetic*, §1.7; Knuth, TAOCP Vol. 2, §4.4. |

## Correctly rounded floating point (`Real`)

| Area | Reference |
|---|---|
| Correct rounding strategy (Ziv's loop) | A. Ziv, "Fast evaluation of elementary mathematical functions with correctly rounded last bit", *ACM TOMS* 17(3), 1991. |
| Semantics (rounding modes, special values, ternary) | L. Fousse, G. Hanrot, V. Lefèvre, P. Pélissier and P. Zimmermann, "MPFR: A multiple-precision binary floating-point library with correct rounding", *ACM TOMS* 33(2), 2007; IEEE Std 754-2019. |
| Series evaluation with rectangular splitting | M. S. Paterson and L. J. Stockmeyer, *SIAM J. Comput.* 2, 1973; D. M. Smith, "Efficient multiple-precision evaluation of elementary functions", *Math. Comp.* 52, 1989. |
| AGM logarithm and π | R. P. Brent, "Fast multiple-precision evaluation of elementary functions", *J. ACM* 23(2), 1976; E. Salamin, "Computation of π using arithmetic-geometric mean", *Math. Comp.* 30, 1976. |
| Binary splitting (e, ln 2) | B. Haible and T. Papanikolaou, "Fast multiprecision evaluation of series of rational numbers", ANTS-III, 1998. |
| Bernoulli numbers via tangent numbers | R. P. Brent and D. Harvey, "Fast computation of Bernoulli, tangent and secant numbers", *Computational and Analytical Mathematics*, Springer Proc. Math. Stat. 50, 2013. |
| Special functions (Stirling series, Euler-Maclaurin, Landen's identity, Hankel expansions, erf, Bessel, Airy, Ei, Li₂) | F. W. J. Olver et al. (eds.), *NIST Digital Library of Mathematical Functions*; M. Abramowitz and I. A. Stegun, *Handbook of Mathematical Functions*, 1964. |
| Normal random deviates | G. E. P. Box and M. E. Muller, "A note on the generation of random normal deviates", *Ann. Math. Statist.* 29, 1958. |
