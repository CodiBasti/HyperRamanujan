"""
Classical and historical approximations for the circumference of an ellipse.
All functions accept arbitrary precision numbers (mpmath.mpf or float)
and compute the normalized circumference C / (pi * (a + b)) as a function of h,
or full circumference C(a, b).
"""

from typing import Callable, Dict
import mpmath as mp


def approx_kepler(h: mp.mpf) -> mp.mpf:
    """
    Kepler's geometric-arithmetic approximation (1609) in normalized form:
    C / (pi * (a + b)) ≈ 1 - h/2  (since 2*sqrt(a*b)/(a+b) = sqrt(1 - h))
    """
    h_mp = mp.mpf(h)
    return (1 + mp.sqrt(1 - h_mp)) / 2


def approx_euler_quadratic(h: mp.mpf) -> mp.mpf:
    """
    Euler / Root-Mean-Square approximation (1773) in normalized form:
    C ≈ 2 * pi * sqrt((a^2 + b^2)/2) -> C / (pi * (a + b)) = sqrt(1 + h)
    """
    h_mp = mp.mpf(h)
    return mp.sqrt(1 + h_mp)


def approx_peano_sipos(h: mp.mpf) -> mp.mpf:
    """
    Sipos (1792) / Peano (1889) arithmetic-geometric mean combination:
    C ≈ pi * (3*(a+b)/2 - sqrt(a*b)) -> C / (pi*(a+b)) = 3/2 - sqrt(1-h)/2
    Matches h^0 and h^1.
    """
    h_mp = mp.mpf(h)
    return mp.mpf("1.5") - mp.sqrt(1 - h_mp) / 2


def approx_ramanujan_I(h: mp.mpf) -> mp.mpf:
    """
    Ramanujan's First Approximation (1914):
    C ≈ pi * (3*(a+b) - sqrt((3*a+b)*(a+3*b)))
    Since (3a+b)(a+3b) = 4*(a+b)^2 - (a-b)^2, in normalized form:
    C / (pi*(a+b)) = 3 - 2*sqrt(1 - h/4)
    Matches h^0, h^1, h^2!
    """
    h_mp = mp.mpf(h)
    return 3 - 2 * mp.sqrt(1 - mp.mpf("0.25") * h_mp)


def approx_hudson_pade(h: mp.mpf) -> mp.mpf:
    """
    Roger W. Hudson (1917) / Padé [1/1] approximation:
    C / (pi*(a+b)) ≈ (64 - 3*h^2) / (64 - 16*h) = 1 + (3*h)/(16 - 4*h)...
    Wait: Padé [1/1] of 1 + h/4 + h^2/64 is 1 + h / (4 - h/4) = (64 + 12*h) / (64 - 4*h).
    Hudson's famous formula is 0.25 * pi * (a+b) * (3*(1 + h/4) + 1/(1 - h/4)) = (64 - 3*h^2)/(64 - 16*h).
    Matches h^0, h^1, h^2, and gives c_3 = 1/256 (exact!).
    """
    h_mp = mp.mpf(h)
    return (64 - 3 * (h_mp**2)) / (64 - 16 * h_mp)


def approx_jacobsen_waadeland_pade22(h: mp.mpf) -> mp.mpf:
    """
    Jacobsen & Waadeland (1985) [2/2] Padé approximation:
    C / (pi*(a+b)) ≈ (256 - 48*h - 21*h^2) / (256 - 112*h + 3*h^2)
    Matches h^0, h^1, h^2, h^3, h^4 purely rationally (0 square roots!).
    """
    h_mp = mp.mpf(h)
    return (256 - 48 * h_mp - 21 * (h_mp**2)) / (256 - 112 * h_mp + 3 * (h_mp**2))


def approx_ramanujan_II(h: mp.mpf) -> mp.mpf:
    """
    Ramanujan's Famous Second Approximation (1914):
    C / (pi*(a+b)) ≈ 1 + 3*h / (10 + sqrt(4 - 3*h))
    Matches h^0, h^1, h^2, h^3, h^4! (5th order term is 190/262144 vs exact 196/262144, defect 3/131072).
    """
    h_mp = mp.mpf(h)
    return 1 + (3 * h_mp) / (10 + mp.sqrt(4 - 3 * h_mp))


def approx_cantrell_exact_end(h: mp.mpf) -> mp.mpf:
    """
    David W. Cantrell (2001) modification of Ramanujan II that is exact at h=1 (4/pi):
    Adds (4/pi - 14/11) * h^12.
    """
    h_mp = mp.mpf(h)
    ram2 = approx_ramanujan_II(h_mp)
    defect = (mp.mpf(4) / mp.pi) - (mp.mpf(14) / mp.mpf(11))
    return ram2 + defect * (h_mp**12)


def get_classical_approximations() -> Dict[str, Callable[[mp.mpf], mp.mpf]]:
    """Returns a registry of classical normalized circumference approximations f(h)."""
    return {
        "Kepler (1609)": approx_kepler,
        "Peano/Sipos (1792)": approx_peano_sipos,
        "Ramanujan I (1914)": approx_ramanujan_I,
        "Hudson [1/1] (1917)": approx_hudson_pade,
        "Jacobsen-Waadeland [2/2] (1985)": approx_jacobsen_waadeland_pade22,
        "Ramanujan II (1914)": approx_ramanujan_II,
        "Cantrell (2001)": approx_cantrell_exact_end,
    }
