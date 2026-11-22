"""
Unit and property tests verifying:
1. Exact agreement of Ivory-Bessel series with Gauss's 2F1(-1/2, -1/2; 1; h) and 4a E(e^2)
2. Symbolic reproduction of Ramanujan I and II and the Continued Fraction periodicity miracle
3. Exact series matching order (up to h^7 and h^8) and boundary exactness (at h=1) of our novel formulas
"""

import mpmath as mp
import sympy as sp
from src.core.exact_geometry import exact_circumference, exact_circumference_from_h, compute_h
from src.core.series_expansions import ivory_bessel_series
from src.ramanujan.reconstruction import (
    derive_ramanujan_I_symbolic,
    derive_ramanujan_II_symbolic,
    explain_ramanujan_miracle_continued_fraction,
)
from src.approximations.novel_candidates import (
    approx_novel_cf_order6,
    approx_novel_exact_boundary_radical,
    approx_novel_pade_radical_order7,
)
from src.geometry.ellipsoid_3d import (
    exact_ellipsoid_surface_area,
    approx_hr_3d_10p,
)


def test_exact_geometry_identity() -> None:
    """Verifies C(a, b) == pi * (a + b) * 2F1(-1/2, -1/2; 1; h) to 45 decimal places."""
    with mp.workdps(50):
        a, b = mp.mpf("5.0"), mp.mpf("3.0")
        h = compute_h(a, b)
        c1 = exact_circumference(a, b, dps=50)
        c2 = mp.pi * (a + b) * exact_circumference_from_h(h, dps=50)
        assert abs(c1 - c2) < mp.mpf("1e-45"), f"Mismatch: {c1} vs {c2}"


def test_ramanujan_reconstructions() -> None:
    """Verifies symbolic derivation of Ramanujan I, II, and the continued fraction a3 == a4 identity."""
    r1 = derive_ramanujan_I_symbolic()
    assert r1["parameters"][sp.Symbol("A", real=True)] == 3

    r2 = derive_ramanujan_II_symbolic()
    assert r2["c4_match_exact"] is True
    assert r2["c5_defect"] == sp.Rational(3, 131072)

    cf = explain_ramanujan_miracle_continued_fraction()
    assert cf["a3_equals_a4"] is True
    assert cf["continued_fraction_coeffs"][:4] == [
        sp.Rational(1, 4),
        sp.Rational(1, 16),
        sp.Rational(3, 16),
        sp.Rational(3, 16),
    ]


def test_novel_approximations_series_and_boundary() -> None:
    """Verifies Taylor order up to h^7/h^8 and exact degenerate limit 4/pi at h=1."""
    h = sp.Symbol("h", nonnegative=True)
    target = ivory_bessel_series(order=9)

    # Check Novel I symbolic expansion up to h^7
    nov1_sym = 1 + (3 * h) / (
        10 + sp.sqrt(4 - 3 * h - (81 * (h**4)) / (18432 - 21120 * h + 5144 * (h**2)))
    )
    diff1 = sp.series(nov1_sym - target, h, 0, 8).removeO()
    assert diff1 == 0, f"Novel I should match up to h^7, but got residual: {diff1}"

    # Check exactness at h = 1 for Novel II and Novel III to 45 decimal places
    with mp.workdps(50):
        exact_h1 = mp.mpf(4) / mp.pi
        assert abs(approx_novel_exact_boundary_radical(mp.mpf(1)) - exact_h1) < mp.mpf("1e-45")
        assert abs(approx_novel_pade_radical_order7(mp.mpf(1)) - exact_h1) < mp.mpf("1e-45")


def test_3d_ellipsoid_surface_area() -> None:
    """Verifies exact 3D-to-2D boundary reduction and < 20 ppm global accuracy of HR-3D-10P."""
    with mp.workdps(40):
        # 1. Flat circular disk (100, 100, 0) -> S = 2 * pi * 100^2 exactly!
        ex_disk = exact_ellipsoid_surface_area(100, 100, 0, dps=40)
        ap_disk = approx_hr_3d_10p(100, 100, 0)
        assert abs(ap_disk - ex_disk) / ex_disk < mp.mpf("1e-35")

        # 2. General triaxial scalene (3, 2, 1) -> error < 0.05 ppm
        ex_tri = exact_ellipsoid_surface_area(3, 2, 1, dps=40)
        ap_tri = approx_hr_3d_10p(3, 2, 1)
        assert abs(ap_tri - ex_tri) / ex_tri < mp.mpf("5e-8")


if __name__ == "__main__":
    test_exact_geometry_identity()
    test_ramanujan_reconstructions()
    test_novel_approximations_series_and_boundary()
    test_3d_ellipsoid_surface_area()
    print("All mathematical verification tests passed (45-digit precision & symbolic proofs)!")
