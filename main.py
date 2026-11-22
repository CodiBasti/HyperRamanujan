"""
Main interactive CLI entry point for the Ellipse Circumference & Ramanujan Research Suite.
Runs:
1. Symbolic Reconstruction of Ramanujan I & II (including the Continued Fraction Periodicity Proof)
2. Discovery & Verification of 3 Novel High-Order Radical Approximations
3. High-Precision Numerical Benchmark across h in [0, 1]
4. Extended Ellipse Geometry Metrics (Elastic Bending Energy / Curvature Integral)
"""

from src.ramanujan.reconstruction import (
    derive_ramanujan_I_symbolic,
    derive_ramanujan_II_symbolic,
    explain_ramanujan_miracle_continued_fraction,
)
from src.benchmarks.error_analysis import print_benchmark_table
from src.benchmarks.pareto_efficiency import print_pareto_table
from src.geometry.extended_metrics import summarize_ellipse_properties
from src.geometry.ellipsoid_3d import print_3d_ellipsoid_benchmark


def main() -> None:
    print("=" * 88)
    print(" ELLIPSE RESEARCH SUITE: RAMANUJAN RECONSTRUCTION & NOVEL APPROXIMATIONS")
    print("=" * 88)

    print("\n[1] SYMBOLIC RECONSTRUCTION OF RAMANUJAN'S FORMULAS (1914)")
    print("-" * 88)
    r1 = derive_ramanujan_I_symbolic()
    print(f"Ramanujan I Ansatz      : C / (pi*(a+b)) ~= {r1['ansatz']}")
    print(f"Unique Series Solution  : {r1['parameters']} -> {r1['formula_h']}")
    print(f"Leading Error Term      : {r1['leading_error_term']}")

    r2 = derive_ramanujan_II_symbolic()
    print(f"\nRamanujan II Ansatz     : C / (pi*(a+b)) ~= {r2['ansatz']}")
    print(f"3-Parameter Solution    : {r2['parameters']} -> {r2['formula_h']}")
    print(f"The 'Ramanujan Miracle' : Matches h^4 coeff automatically! (Exact = {r2['c4_exact']}, RamII = {r2['c4_ramanujan']})")
    print(f"First Non-Zero Defect   : delta(c_5) = {r2['c5_exact']} - {r2['c5_ramanujan']} = {r2['c5_defect']} * h^5")

    cf = explain_ramanujan_miracle_continued_fraction()
    print("\nWHY DOES RAMANUJAN II MATCH 4 ORDERS WITH ONLY 3 PARAMETERS?")
    print(f"Stieltjes Continued Fraction coefficients of (2F1(-1/2, -1/2; 1; h) - 1):")
    print(f"  [a_1, a_2, a_3, a_4, a_5] = {cf['continued_fraction_coeffs']}")
    print(f"  Notice that a_3 = a_4 = 3/16 identically ({cf['a3_equals_a4']})!")
    print(f"  Replacing the tail from level 3 by a 1-periodic continued fraction T = 1 - (3h/16)/T")
    print(f"  yields the quadratic root T(h) = {cf['periodic_tail_T']} = (2 + sqrt(4 - 3h))/4,")
    print(f"  which collapses directly into 1 + 3h / (10 + sqrt(4 - 3h))!")

    print("\n[2] HIGH-PRECISION BENCHMARK: CLASSICAL VS. NOVEL APPROXIMATIONS (50 DPS)")
    print("-" * 88)
    print_benchmark_table()

    print("\n[3] PARETO-EFFICIENCY & COMPLEXITY RANKING ('WHAT IS THE BEST FORMULA?')")
    print("-" * 88)
    print_pareto_table()

    print("\n[4] EXTENDED ELLIPSE GEOMETRY: ELASTIC BENDING ENERGY (CURVATURE INTEGRAL)")
    print("-" * 88)
    for a, b in [(1.0, 1.0), (2.0, 1.0), (5.0, 1.0)]:
        props = summarize_ellipse_properties(a, b)
        print(
            f"a={a:.1f}, b={b:.1f} (h={props['h']:.4f}) -> "
            f"C_exact={props['circumference_exact']:.6f}, "
            f"W_norm_exact={props['bending_energy_norm_exact']:.6f}, "
            f"W_norm_approx={props['bending_energy_norm_approx']:.6f} "
            f"(Rel.Err={props['bending_energy_rel_error']:.2e})"
        )

    print("\n[5] 3D TRIAXIAL ELLIPSOID SURFACE AREA: THOMSEN (2004) VS. HYPER-RAMANUJAN 3D")
    print("-" * 88)
    print_3d_ellipsoid_benchmark()


if __name__ == "__main__":
    main()
