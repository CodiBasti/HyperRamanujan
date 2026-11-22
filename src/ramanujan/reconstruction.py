"""
Reconstruction of Srinivasa Ramanujan's 1914 derivations for the circumference
of an ellipse ("Modular Equations and Approximations to pi", Section 16).

Explores three mathematical pathways capable of reproducing Ramanujan I and II
purely from analytical first principles:
1. Radical-Linear Ansatz (Undetermined Coefficients matching Ivory-Bessel series)
2. Quadratic Padé / Continued Fraction Reduction (Algebraic Approximation of 2F1)
3. Modular / Dedekind Eta / AGM Perturbation Approach
"""

from typing import Dict, Any
import sympy as sp
from src.core.series_expansions import ivory_bessel_series


def derive_ramanujan_I_symbolic() -> Dict[str, Any]:
    """
    Reconstructs Ramanujan's First Approximation:
    C / (pi * (a + b)) ≈ A - B * sqrt(1 - C * h)
    by matching the first 3 terms (h^0, h^1, h^2) of the Ivory-Bessel series:
    1 + h/4 + h^2/64 + ...
    """
    h = sp.Symbol("h", nonnegative=True)
    A, B, C_param = sp.symbols("A B C", real=True)

    ansatz = A - B * sp.sqrt(1 - C_param * h)
    target_series = ivory_bessel_series(order=4)

    # Expand ansatz up to h^2
    ansatz_series = sp.series(ansatz, h, 0, 3).removeO()
    diff_poly = sp.Poly(sp.expand(ansatz_series - ivory_bessel_series(order=3)), h)

    eqs = [sp.Eq(coeff, 0) for coeff in diff_poly.all_coeffs()]
    sol = sp.solve(eqs, (A, B, C_param), dict=True)[0]

    reconstructed = ansatz.subs(sol)
    error_series = sp.series(target_series - reconstructed, h, 0, 4)

    # Show geometric form in terms of a and b
    a, b = sp.symbols("a b", positive=True)
    h_ab = ((a - b) / (a + b)) ** 2
    geom_expr = sp.simplify(sp.pi * (a + b) * reconstructed.subs(h, h_ab))

    return {
        "ansatz": ansatz,
        "parameters": sol,
        "formula_h": reconstructed,
        "formula_ab": geom_expr,
        "leading_error_term": error_series,
    }


def derive_ramanujan_II_symbolic() -> Dict[str, Any]:
    """
    Reconstructs Ramanujan's Famous Second Approximation:
    C / (pi * (a + b)) ≈ 1 + (p * h) / (q + sqrt(4 - r * h))
    by matching the Ivory-Bessel expansion:
    1 + h/4 + h^2/64 + h^3/256 + 25*h^4/16384 + 49*h^5/65536 + ...

    Remarkable property ("Ramanujan's Miracle"):
    There are only 3 free parameters (p, q, r), which naturally match h^1, h^2, h^3.
    Yet for the unique rational solution (p=3, q=10, r=3), the h^4 coefficient ALSO
    vanishes identically (matching 25/16384 = 100/65536)!
    Why? Let's reveal the underlying quadratic equation for the function y(h) = (F(h) - 1)/h!
    """
    h = sp.Symbol("h", nonnegative=True)
    p, q, r = sp.symbols("p q r", real=True)

    ansatz = 1 + (p * h) / (q + sp.sqrt(4 - r * h))
    target_series = ivory_bessel_series(order=6)

    # Expand ansatz up to h^5
    ansatz_series = sp.series(ansatz, h, 0, 6).removeO()
    diff_expr = sp.expand(ansatz_series - target_series)

    # Extract coefficients for h^1, h^2, h^3
    eq1 = sp.Eq( diff_expr.coeff(h, 1), 0 )
    eq2 = sp.Eq( diff_expr.coeff(h, 2), 0 )
    eq3 = sp.Eq( diff_expr.coeff(h, 3), 0 )

    sol = sp.solve([eq1, eq2, eq3], (p, q, r), dict=True)[0]
    reconstructed = ansatz.subs(sol)

    # Check the h^4 and h^5 coefficients with this 3-parameter solution!
    reconstructed_series = sp.series(reconstructed, h, 0, 6).removeO()
    c4_exact = target_series.coeff(h, 4)
    c4_ramanujan = reconstructed_series.coeff(h, 4)
    c5_exact = target_series.coeff(h, 5)
    c5_ramanujan = reconstructed_series.coeff(h, 5)

    return {
        "ansatz": ansatz,
        "parameters": sol,
        "formula_h": reconstructed,
        "c4_exact": c4_exact,
        "c4_ramanujan": c4_ramanujan,
        "c4_match_exact": (c4_exact == c4_ramanujan),
        "c5_exact": c5_exact,
        "c5_ramanujan": c5_ramanujan,
        "c5_defect": c5_exact - c5_ramanujan,
    }


def explain_ramanujan_miracle_continued_fraction() -> Dict[str, Any]:
    """
    Explains HOW Ramanujan found 1 + 3h / (10 + sqrt(4 - 3h)) without brute-force!
    If we write the hypergeometric series as a Gauss Continued Fraction (or Stieltjes fraction):
    F(h) = 1 + (h/4) / (1 - (h/16) / (1 - (3h/16) / (1 - (3h/16) / (1 - (11h/48) / ...))))

    Look at the partial denominators:
    a_1 = 1/4
    a_2 = 1/16
    a_3 = 3/16
    a_4 = 3/16
    a_5 = 11/48
    When the tail of the continued fraction repeats periodically with period 1 at rate (3/16)*h,
    a periodic continued fraction T = 1 - (3h/16) / T satisfies T^2 - T + 3h/16 = 0,
    which solves to T = (1 + sqrt(1 - 3h/4)) / 2 = (2 + sqrt(4 - 3h)) / 4!
    Substituting this periodic tail T back into the 2nd level of the continued fraction yields
    Ramanujan's Second Approximation directly AND explains why 3 parameters match 4 orders!
    """
    h = sp.Symbol("h", nonnegative=True)
    target = ivory_bessel_series(order=7)

    # Compute J-fraction / S-fraction coefficients step by step:
    # F(h) - 1 = a1 * h / (1 - a2 * h / (1 - a3 * h / (1 - a4 * h / (1 - a5 * h / ...))))
    rem = target - 1
    cf_coeffs = []
    curr = rem
    for k in range(1, 6):
        # curr = a_k * h * (1 + O(h))
        s = sp.series(curr, h, 0, 7 - k).removeO()
        a_k = s.coeff(h, 1)
        cf_coeffs.append(a_k)
        # Next: curr = a_k * h / (1 - next_curr) => 1 - next_curr = a_k * h / curr
        next_curr = sp.series(1 - (a_k * h) / curr, h, 0, 7 - k).removeO()
        curr = next_curr

    # Notice cf_coeffs: [1/4, 1/16, 3/16, 3/16, 11/48]!
    # Because a_3 = 3/16 and a_4 = 3/16 ARE IDENTICAL, approximating the tail from level 3
    # onwards as a constant periodic continued fraction with parameter alpha = 3/16
    # automatically matches BOTH level 3 (h^3) AND level 4 (h^4)!
    T = (1 + sp.sqrt(1 - 4 * cf_coeffs[2] * h)) / 2
    ramanujan_from_cf = sp.simplify(1 + (cf_coeffs[0] * h) / (1 - (cf_coeffs[1] * h) / T))

    return {
        "continued_fraction_coeffs": cf_coeffs,
        "a3_equals_a4": cf_coeffs[2] == cf_coeffs[3],
        "periodic_tail_T": T,
        "reconstructed_from_periodic_cf": ramanujan_from_cf,
    }
