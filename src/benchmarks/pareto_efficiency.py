"""
Objective mathematical evaluation of 'The Best' approximation formula via
Pareto-Optimality and Computational Efficiency Scoring.

Balances 4 orthogonal metrics:
1. Computational Hardware Cost (CPU/FPU cycle weight: +, -, * = 1; / = 4; sqrt = 6; pow = 30)
2. Symbolic / Constant Description Complexity (number of constants & non-integer penalty)
3. Practical Asymptotic Precision (Guaranteed digits across moderate ellipses h in [0, 0.5])
4. Global Worst-Case Precision (Chebyshev L_infinity guaranteed digits across all h in [0, 1])
"""

from typing import Dict, Any, List
import mpmath as mp
from src.benchmarks.error_analysis import evaluate_approximations


# Hardware cycle weights & structural metadata for each formula
FORMULA_COMPLEXITY_SPEC: Dict[str, Dict[str, Any]] = {
    "Kepler (1609)": {
        "add_mul": 3, "div": 1, "sqrt": 1, "pow": 0, "transcendental_consts": 0, "taylor_order": 1
    },
    "Peano/Sipos (1792)": {
        "add_mul": 3, "div": 1, "sqrt": 1, "pow": 0, "transcendental_consts": 0, "taylor_order": 1
    },
    "Ramanujan I (1914)": {
        "add_mul": 4, "div": 0, "sqrt": 1, "pow": 0, "transcendental_consts": 0, "taylor_order": 2
    },
    "Hudson [1/1] (1917)": {
        "add_mul": 5, "div": 1, "sqrt": 0, "pow": 0, "transcendental_consts": 0, "taylor_order": 3
    },
    "Jacobsen-Waadeland [2/2] (1985)": {
        "add_mul": 9, "div": 1, "sqrt": 0, "pow": 0, "transcendental_consts": 0, "taylor_order": 4
    },
    "Ramanujan II (1914)": {
        "add_mul": 5, "div": 1, "sqrt": 1, "pow": 0, "transcendental_consts": 0, "taylor_order": 4
    },
    "Cantrell (2001)": {
        "add_mul": 11, "div": 1, "sqrt": 1, "pow": 0, "transcendental_consts": 1, "taylor_order": 4
    },
    "HR-7-Int (Order-7, pure integers)": {
        "add_mul": 12, "div": 2, "sqrt": 1, "pow": 0, "transcendental_consts": 0, "taylor_order": 7
    },
    "HR-6B (Order-6 + exact h=1)": {
        "add_mul": 12, "div": 2, "sqrt": 1, "pow": 0, "transcendental_consts": 1, "taylor_order": 6
    },
    "HR-8B (Order-8 + exact h=1)": {
        "add_mul": 16, "div": 2, "sqrt": 1, "pow": 0, "transcendental_consts": 1, "taylor_order": 8
    },
    "HR-8P (Order-8, -13/14 pole + exact h=1)": {
        "add_mul": 18, "div": 2, "sqrt": 1, "pow": 0, "transcendental_consts": 1, "taylor_order": 8
    },
    "HR-6M (Order-6, Minimax [3/1] + exact h=1)": {
        "add_mul": 16, "div": 2, "sqrt": 1, "pow": 0, "transcendental_consts": 3, "taylor_order": 6
    },
}


def compute_op_cost(spec: Dict[str, Any]) -> float:
    """
    Computes estimated FPU operational cost in normalized cycle units:
    +, -, * : 1.0 unit
    /       : 4.0 units
    sqrt    : 6.0 units
    pow/h^12: 12.0 units (or 4 squarings + multiplications, plus 3.0 per non-integer constant)
    """
    cycle_cost = (
        1.0 * spec["add_mul"]
        + 4.0 * spec["div"]
        + 6.0 * spec["sqrt"]
        + 12.0 * spec["pow"]
        + 2.5 * spec["transcendental_consts"]
    )
    return cycle_cost


def evaluate_pareto_efficiency() -> List[Dict[str, Any]]:
    """
    Evaluates all approximations on a dense grid of h in [0, 1] (including near-degenerate
    ellipses h = 0.98, 0.99, 0.999) and computes:
    - Digits at h = 0.1 (typical ellipse a/b ~ 1.92): D_typical = -log10(rel_err(0.1))
    - Digits worst-case in h in [0, 0.5] (a/b <= 5.83): D_practical
    - Digits worst-case globally in h in [0, 1]: D_global = -log10(max_rel_err)
    - Combined Efficiency Score = (0.4 * D_typical + 0.3 * D_practical + 0.3 * D_global) / sqrt(OpCost)
    """
    h_dense = [
        0.01, 0.05, 0.1, 0.2, 0.25, 0.35, 0.5,
        0.65, 0.75, 0.85, 0.9, 0.9425, 0.95, 0.969, 0.98, 0.99, 0.9945, 0.995, 0.998, 0.999, 0.9999, 0.99999, 1.0
    ]
    bench = evaluate_approximations(h_dense, dps=50)
    idx_0_1 = h_dense.index(0.1)
    idx_0_5 = h_dense.index(0.5)

    report = []
    for name, data in bench.items():
        spec = FORMULA_COMPLEXITY_SPEC[name]
        cost = compute_op_cost(spec)
        errs = data["relative_errors"]

        err_0_1 = max(errs[idx_0_1], 1e-16)
        err_pract = max(max(errs[: idx_0_5 + 1]), 1e-16)
        err_glob = max(max(errs), 1e-16)

        d_typ = float(-mp.log10(err_0_1))
        d_pract = float(-mp.log10(err_pract))
        d_glob = float(-mp.log10(err_glob))

        # Efficiency Index: Precision bits/digits gained per square-root of computational cost
        efficiency_score = (0.4 * d_typ + 0.3 * d_pract + 0.3 * d_glob) / (cost ** 0.5)

        report.append({
            "name": name,
            "taylor_order": spec["taylor_order"],
            "sqrt_count": spec["sqrt"],
            "op_cost": cost,
            "digits_at_h_0_1": d_typ,
            "digits_practical_max_h_0_5": d_pract,
            "digits_global_worst_case": d_glob,
            "max_global_ppm": err_glob * 1e6,
            "efficiency_score": efficiency_score,
        })

    return report


def print_pareto_table() -> None:
    """Prints the Pareto-Optimality & Efficiency Ranking table."""
    report = evaluate_pareto_efficiency()
    header = (
        f"| {'Formula / Method':<52} | {'Order':<5} | {'Cost':<5} | "
        f"{'Digits (h=0.1)':<14} | {'Digits (h<=0.5)':<15} | {'Min Digits [0,1]':<16} | "
        f"{'Max Err (ppm)':<13} | {'Eff. Score':<10} |"
    )
    sep = (
        "|" + "-" * 54 + "|" + "-" * 7 + "|" + "-" * 7 + "|"
        + "-" * 16 + "|" + "-" * 17 + "|" + "-" * 18 + "|"
        + "-" * 15 + "|" + "-" * 12 + "|"
    )
    print(header)
    print(sep)
    for row in report:
        print(
            f"| {row['name']:<52} | O(h^{row['taylor_order']}) | {row['op_cost']:<5.1f} | "
            f"{row['digits_at_h_0_1']:<14.2f} | {row['digits_practical_max_h_0_5']:<15.2f} | "
            f"{row['digits_global_worst_case']:<16.2f} | {row['max_global_ppm']:<13.2f} | "
            f"{row['efficiency_score']:<10.3f} |"
        )


if __name__ == "__main__":
    print_pareto_table()
