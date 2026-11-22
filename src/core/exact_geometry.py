"""
High-precision ground truth calculations for ellipse properties using mpmath.
Supports arbitrary precision for circumference, area, focal distance, and shape invariants.
"""

from typing import Dict
import mpmath as mp


def compute_h(a: mp.mpf, b: mp.mpf) -> mp.mpf:
    """Computes the Hudson/Ramanujan symmetry invariant h = ((a - b) / (a + b))^2."""
    a_mp, b_mp = mp.mpf(a), mp.mpf(b)
    if a_mp + b_mp == 0:
        return mp.mpf(0)
    return ((a_mp - b_mp) / (a_mp + b_mp)) ** 2


def compute_eccentricity(a: mp.mpf, b: mp.mpf) -> mp.mpf:
    """Computes the numerical eccentricity e = sqrt(1 - (b/a)^2) assuming a >= b."""
    a_mp, b_mp = mp.mpf(a), mp.mpf(b)
    if a_mp < b_mp:
        a_mp, b_mp = b_mp, a_mp
    if a_mp == 0:
        return mp.mpf(0)
    return mp.sqrt(1 - (b_mp / a_mp) ** 2)


def exact_circumference(a: mp.mpf, b: mp.mpf, dps: int = 50) -> mp.mpf:
    """
    Computes the exact circumference of an ellipse with semi-axes a and b
    to `dps` decimal places using the complete elliptic integral of the second kind:
    C = 4 * a * E(e^2) = pi * (a + b) * 2F1(-1/2, -1/2; 1; h)
    """
    with mp.workdps(dps):
        a_mp, b_mp = mp.mpf(a), mp.mpf(b)
        if a_mp < b_mp:
            a_mp, b_mp = b_mp, a_mp
        if a_mp == 0:
            return mp.mpf(0)
        if b_mp == 0:
            return 4 * a_mp
        m = 1 - (b_mp / a_mp) ** 2  # m = e^2
        return 4 * a_mp * mp.ellipe(m)


def exact_circumference_from_h(h: mp.mpf, dps: int = 50) -> mp.mpf:
    """
    Computes the normalized circumference C / (pi * (a + b)) directly from h
    using Gauss's hypergeometric function 2F1(-1/2, -1/2; 1; h).
    """
    with mp.workdps(dps):
        h_mp = mp.mpf(h)
        if h_mp == 1:
            return mp.mpf(4) / mp.pi
        return mp.hyp2f1(mp.mpf("-0.5"), mp.mpf("-0.5"), 1, h_mp)


def exact_ellipse_metrics(a: mp.mpf, b: mp.mpf, dps: int = 50) -> Dict[str, mp.mpf]:
    """
    Returns a dictionary of fundamental geometric properties of an ellipse
    (Useful for exploring relationships beyond circumference).
    """
    with mp.workdps(dps):
        a_mp, b_mp = mp.mpf(a), mp.mpf(b)
        if a_mp < b_mp:
            a_mp, b_mp = b_mp, a_mp
        h = compute_h(a_mp, b_mp)
        e = compute_eccentricity(a_mp, b_mp)
        c_linear = mp.sqrt(a_mp**2 - b_mp**2)
        area = mp.pi * a_mp * b_mp
        circ = exact_circumference(a_mp, b_mp, dps=dps)
        return {
            "a": a_mp,
            "b": b_mp,
            "h": h,
            "eccentricity": e,
            "linear_eccentricity_c": c_linear,
            "area": area,
            "circumference": circ,
            "semi_latus_rectum": (b_mp**2) / a_mp if a_mp != 0 else mp.mpf(0),
            "arithmetic_mean": (a_mp + b_mp) / 2,
            "geometric_mean": mp.sqrt(a_mp * b_mp),
            "quadratic_mean": mp.sqrt((a_mp**2 + b_mp**2) / 2),
        }
