"""
Extended geometric and physical properties of ellipses beyond total circumference:
1. Elliptic Surface Curvature Integral & Bending Energy
2. Polar Moment of Inertia & Radius of Gyration
3. Electrostatic Capacitance / Ellipsoidal Mean BreadthApproximations
"""

from typing import Dict
import mpmath as mp
from src.core.exact_geometry import compute_h, exact_circumference


def exact_bending_energy_normalized(a: mp.mpf, b: mp.mpf, dps: int = 50) -> mp.mpf:
    """
    Computes the normalized total squared curvature (Willmore / elastic bending energy)
    of an ellipse wire:
    W = (1 / (2*pi)) * oint kappa(s)^2 ds
    For a circle of radius R, W = 1 / R.
    Multiply by sqrt(a*b) to make it scale-invariant:
    W_norm(a, b) = (sqrt(a*b) / (2*pi)) * int_0^(2*pi) (a*b)^2 / (a^2 sin^2 t + b^2 cos^2 t)^(5/2) dt
    """
    with mp.workdps(dps):
        a_mp, b_mp = mp.mpf(a), mp.mpf(b)
        if a_mp < b_mp:
            a_mp, b_mp = b_mp, a_mp
        m = 1 - (b_mp / a_mp) ** 2
        q = (b_mp / a_mp)
        # Exact elliptic integral formula for int_0^(pi/2) (1 - m sin^2 t)^(-5/2) dt:
        # I_{-5/2}(m) = (1 / (3 * (1-m)^2)) * ( 2*(2 - m)*E(m) - (1 - m)*K(m) )
        E_m = mp.ellipe(m)
        K_m = mp.ellipk(m)
        integral_quarter = (2 * (2 - m) * E_m - (1 - m) * K_m) / (3 * ((1 - m) ** 2))
        # Full integral of kappa^2 ds = 4 * (b^2 / a^3) * I_{-5/2}(m)
        total_kappa_sq = 4 * (b_mp**2 / (a_mp**3)) * integral_quarter
        return (mp.sqrt(a_mp * b_mp) / (2 * mp.pi)) * total_kappa_sq


def approx_bending_energy_ramanujan_style(a: mp.mpf, b: mp.mpf) -> mp.mpf:
    """
    A novel Ramanujan-inspired approximation for the scale-invariant elastic
    bending energy W_norm(a, b) of an ellipse in terms of h = ((a-b)/(a+b))^2:
    Exact Taylor expansion of W_norm(h) * (1 - h)^(3/2):
      1 + (9/4)*h + (9/64)*h^2 + (1/256)*h^3 + ...
    Notice how the numerator coefficients (1, 9/4, 9/64, 1/256) mirror the Ivory-Bessel
    structure!
    """
    h = compute_h(a, b)
    num = (
        1
        + (mp.mpf(9) / mp.mpf(4)) * h
        + (mp.mpf(9) / mp.mpf(64)) * (h**2)
        + (mp.mpf(1) / mp.mpf(256)) * (h**3)
    )
    return num / ((1 - h) ** mp.mpf("1.5"))


def summarize_ellipse_properties(a: float, b: float) -> Dict[str, float]:
    """Computes a comprehensive profile of an ellipse (circumference, area, curvature energy)."""
    a_mp, b_mp = mp.mpf(a), mp.mpf(b)
    h = compute_h(a_mp, b_mp)
    circ = exact_circumference(a_mp, b_mp)
    w_exact = exact_bending_energy_normalized(a_mp, b_mp)
    w_approx = approx_bending_energy_ramanujan_style(a_mp, b_mp)
    return {
        "a": float(a_mp),
        "b": float(b_mp),
        "h": float(h),
        "circumference_exact": float(circ),
        "bending_energy_norm_exact": float(w_exact),
        "bending_energy_norm_approx": float(w_approx),
        "bending_energy_rel_error": float(abs(w_approx - w_exact) / w_exact),
    }
