"""
Objective mathematical evaluation of 'The Best' approximation formula via:
1. Strict 2D Pareto-Frontier Dominance Analysis (OpCost vs. Guaranteed Precision — zero arbitrary weights!)
2. Information-Theoretic Efficiency (Bits of precision gained per FPU cycle unit)
   evaluated separately for:
   - Practical Ellipses (h <= 0.5, i.e. aspect ratio a/b <= 5.83): Eta_pract = D_typical / Cost
   - Global Worst-Case (h in [0, 1], all eccentricities up to flat line): Eta_global = D_global / Cost
   - Balanced Harmonic/Weighted Index: Eff_Score
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
    "Matt Parker (2020)": {
        "add_mul": 9, "div": 0, "sqrt": 1, "pow": 0, "transcendental_consts": 0, "taylor_order": 0
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
    pow     : 12.0 units
    non-int : 2.5 units (memory/LUT load penalty for irrational/fitted constants)
    """
    return (
        1.0 * spec["add_mul"]
        + 4.0 * spec["div"]
        + 6.0 * spec["sqrt"]
        + 12.0 * spec["pow"]
        + 2.5 * spec["transcendental_consts"]
    )


def evaluate_pareto_efficiency() -> List[Dict[str, Any]]:
    """
    Evaluates all approximations on a dense grid of h in [0, 1] (including h=0 and near-degenerate
    ellipses h = 0.9945, 0.998, 0.9999, 1.0) and determines:
    1. Strict 2D Pareto Frontier status (`is_pareto_optimal`):
       A formula is strictly Pareto-optimal if NO other formula has <= OpCost AND strictly higher
       precision in BOTH practical (h <= 0.5) and global (h in [0, 1]) regimes.
    2. Digits at h = 0.1 (D_typ), worst-case digits in h <= 0.5 (D_pract), and global worst-case digits (D_glob).
    """
    h_dense = [
        0.0, 0.01, 0.05, 0.1, 0.2, 0.25, 0.35, 0.5,
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

    # Determine strict multi-objective Pareto Frontier (Cost vs. D_pract vs. D_glob)
    for row in report:
        dominated = False
        for other in report:
            if other["name"] == row["name"]:
                continue
            if (
                other["op_cost"] <= row["op_cost"]
                and other["digits_practical_max_h_0_5"] >= row["digits_practical_max_h_0_5"]
                and other["digits_global_worst_case"] >= row["digits_global_worst_case"]
                and (
                    other["op_cost"] < row["op_cost"]
                    or other["digits_practical_max_h_0_5"] > row["digits_practical_max_h_0_5"]
                    or other["digits_global_worst_case"] > row["digits_global_worst_case"]
                )
            ):
                dominated = True
                break
        row["pareto_frontier"] = not dominated

    return report


def print_pareto_table() -> None:
    """Prints the Pareto-Optimality & Efficiency Ranking table."""
    report = evaluate_pareto_efficiency()
    header = (
        f"| {'Formula / Method':<44} | {'Order':<5} | {'Cost':<5} | "
        f"{'Digits (h=0.1)':<14} | {'Min (h<=0.5)':<12} | {'Min [0,1]':<9} | "
        f"{'Max Err (ppm)':<13} | {'Pareto?':<7} | {'Eff. Score':<10} |"
    )
    sep = (
        "|" + "-" * 46 + "|" + "-" * 7 + "|" + "-" * 7 + "|"
        + "-" * 16 + "|" + "-" * 14 + "|" + "-" * 11 + "|"
        + "-" * 15 + "|" + "-" * 9 + "|" + "-" * 12 + "|"
    )
    print(header)
    print(sep)
    for row in report:
        p_flag = "YES *" if row["pareto_frontier"] else "no"
        print(
            f"| {row['name']:<44} | O(h^{row['taylor_order']}) | {row['op_cost']:<5.1f} | "
            f"{row['digits_at_h_0_1']:<14.2f} | {row['digits_practical_max_h_0_5']:<12.2f} | "
            f"{row['digits_global_worst_case']:<9.2f} | {row['max_global_ppm']:<13.2f} | "
            f"{p_flag:<7} | {row['efficiency_score']:<10.3f} |"
        )


if __name__ == "__main__":
    print_pareto_table()
