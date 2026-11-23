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


def approx_matt_parker(h: mp.mpf) -> mp.mpf:
    """
    Matt Parker's 'Stand-up Maths' Approximation (2020):
    C ≈ pi * (53/3 * a + 717/35 * b - sqrt(269*a^2 + 667*a*b + 371*b^2))
    An asymmetric minimax modification of Ramanujan I optimized for moderate aspect
    ratios a/b in [1, 5] (h in [0, 0.44]), sacrificing circle exactness (37.01 ppm error at h=0).
    In terms of s = sqrt(h) = (a-b)/(a+b), with 2*a/(a+b) = 1+s and 2*b/(a+b) = 1-s:
    """
    h_mp = mp.mpf(h)
    s = mp.sqrt(h_mp)
    a_norm = 1 + s
    b_norm = 1 - s
    val = (
        (mp.mpf(53) / mp.mpf(3)) * a_norm
        + (mp.mpf(717) / mp.mpf(35)) * b_norm
        - mp.sqrt(269 * (a_norm**2) + 667 * a_norm * b_norm + 371 * (b_norm**2))
    )
    return val / 2


def approx_hudson_pade(h: mp.mpf) -> mp.mpf:
    """
    Ralph G. Hudson (1917) / Padé [2/1] approximation:
    C / (pi*(a+b)) ≈ (64 - 3*h^2) / (64 - 16*h).
    Hudson's famous formula is 0.25 * pi * (a+b) * (3*(1 + h/4) + 1/(1 - h/4)) = (64 - 3*h^2)/(64 - 16*h).
    Matches h^0, h^1, h^2, and h^3 (c_3 = 1/256 exact!).
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
    David W. Cantrell (2004) modification of Ramanujan II that is exact at h=1 (4/pi):
    Adds (4/pi - 14/11) * h^12.
    """
    h_mp = mp.mpf(h)
    ram2 = approx_ramanujan_II(h_mp)
    defect = (mp.mpf(4) / mp.pi) - (mp.mpf(14) / mp.mpf(11))
    return ram2 + defect * (h_mp**12)


def approx_ayala_raggi_r2_2exp(h: mp.mpf) -> mp.mpf:
    """
    Salvador E. Ayala-Raggi & Manuel Rendón-Marín (2025/2026, arXiv:2510.16191v2, AppliedMath):
    Transcendental exponential correction to Ramanujan II (`R2/2exp` with sigmoid damping):
    F(h) = F_RamII(h) / (1 - (A*exp(-B*(1-h)) + C*exp(-D*(1-h))) * sigma(h))
    where sigma(h) = 1 / (1 + exp(-60*(h - 0.33))).
    Uses 3 transcendental exp() evaluations + 1 sqrt() + 4 empirically fitted float64 constants.
    """
    h_mp = mp.mpf(h)
    ram2 = approx_ramanujan_II(h_mp)
    A = mp.mpf("3.37528e-4")
    B = mp.mpf("10.29662")
    C = mp.mpf("6.48093e-5")
    D = mp.mpf("40.89043")
    sigma = mp.mpf(1) / (1 + mp.exp(-60 * (h_mp - mp.mpf("0.33"))))
    corr = (A * mp.exp(-B * (1 - h_mp)) + C * mp.exp(-D * (1 - h_mp))) * sigma
    return ram2 / (1 - corr)


def approx_shankar_a1(h: mp.mpf) -> mp.mpf:
    """
    Uday Shankar (March 2026, arXiv:2603.03373, Attempt 1):
    Leading-order inner-radical perturbation of Ramanujan II matching through h^5:
    F(h) = 1 + 192*h / (640 + sqrt(16384 - 12288*h - 18*h^4))
         = 1 + 3*h / (10 + sqrt(4 - 3*h - (9/2048)*h^4)).
    """
    h_mp = mp.mpf(h)
    return 1 + (192 * h_mp) / (640 + mp.sqrt(16384 - 12288 * h_mp - 18 * (h_mp**4)))


def approx_shankar_a2(h: mp.mpf) -> mp.mpf:
    """
    Uday Shankar (March 2026, arXiv:2603.03373, Attempt 2):
    Depth-4 Stieltjes continued fraction closed with constant empirical tail rate 0.915/4:
    F(h) = ((128 - 48*h - 15*h^2) + (128 - 6*h^2)*S) / ((128 - 80*h + 3*h^2) + (128 - 32*h)*S),
    where S = sqrt(1 - 0.915*h). Matches through h^4 (near-h^5 via 0.915/4 ~= 11/48).
    """
    h_mp = mp.mpf(h)
    S = mp.sqrt(1 - mp.mpf("0.915") * h_mp)
    num = (128 - 48 * h_mp - 15 * (h_mp**2)) + (128 - 6 * (h_mp**2)) * S
    den = (128 - 80 * h_mp + 3 * (h_mp**2)) + (128 - 32 * h_mp) * S
    return num / den


def get_classical_approximations() -> Dict[str, Callable[[mp.mpf], mp.mpf]]:
    """Returns a registry of classical and modern literature approximations f(h)."""
    return {
        "Kepler (1609)": approx_kepler,
        "Peano/Sipos (1792)": approx_peano_sipos,
        "Ramanujan I (1914)": approx_ramanujan_I,
        "Matt Parker (2020)": approx_matt_parker,
        "Hudson [2/1] (1917)": approx_hudson_pade,
        "Jacobsen-Waadeland [2/2] (1985)": approx_jacobsen_waadeland_pade22,
        "Ramanujan II (1914)": approx_ramanujan_II,
        "Cantrell (2004)": approx_cantrell_exact_end,
        "Shankar A1 (2026)": approx_shankar_a1,
        "Shankar A2 (2026)": approx_shankar_a2,
        "Ayala-Raggi & Rendon-Marin (2026, R2/2exp)": approx_ayala_raggi_r2_2exp,
    }
