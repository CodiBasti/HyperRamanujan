"""
Symbolic series expansions for the ellipse circumference in terms of
eccentricity e^2 and the symmetric parameter h = ((a-b)/(a+b))^2.
"""

from typing import List
import sympy as sp


def ivory_bessel_coefficient(n: int) -> sp.Rational:
    """
    Returns the n-th coefficient c_n in the Ivory-Bessel series:
    C / (pi * (a + b)) = sum_{n=0}^inf (binom(1/2, n))^2 * h^n
    """
    half = sp.Rational(1, 2)
    binom_val = sp.binomial(half, n)
    return sp.simplify(binom_val**2)


def ivory_bessel_series(order: int = 8) -> sp.Expr:
    """
    Returns the symbolic Taylor expansion of C / (pi * (a + b)) in terms of h
    up to O(h^order).
    """
    h = sp.Symbol("h", nonnegative=True)
    terms: List[sp.Expr] = [
        ivory_bessel_coefficient(n) * (h**n) for n in range(order)
    ]
    return sum(terms)


def get_ivory_bessel_coeffs(count: int = 10) -> List[sp.Rational]:
    """Returns a list of the first `count` exact rational coefficients of the h-series."""
    return [ivory_bessel_coefficient(n) for n in range(count)]
