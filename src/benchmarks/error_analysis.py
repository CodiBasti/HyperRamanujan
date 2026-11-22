"""
High-precision benchmarking suite comparing classical approximations
(Kepler, Euler, Ramanujan I, Hudson, Ramanujan II, Cantrell) against our newly
discovered continued-fraction and boundary-exact radical formulas across the entire
domain h in [0, 1].
"""

from typing import Dict, List, Any
import mpmath as mp
from src.core.exact_geometry import exact_circumference_from_h
from src.approximations.classical import get_classical_approximations
from src.approximations.novel_candidates import get_novel_approximations


def evaluate_approximations(
    h_values: List[float] = None, dps: int = 50
) -> Dict[str, Dict[str, Any]]:
    """
    Evaluates all registered approximations across a grid of h values
    and computes relative errors |C_approx - C_exact| / C_exact.
    """
    if h_values is None:
        h_values = [0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99, 1.0]

    with mp.workdps(dps):
        exact_vals = [exact_circumference_from_h(mp.mpf(str(h)), dps=dps) for h in h_values]
        all_models = {}
        all_models.update(get_classical_approximations())
        all_models.update(get_novel_approximations())

        results: Dict[str, Dict[str, Any]] = {}
        for name, func in all_models.items():
            rel_errors = []
            for h_val, ex_val in zip(h_values, exact_vals):
                h_mp = mp.mpf(str(h_val))
                ap_val = func(h_mp)
                rel_err = abs(ap_val - ex_val) / ex_val
                rel_errors.append(float(rel_err))

            results[name] = {
                "h_values": h_values,
                "relative_errors": rel_errors,
                "max_rel_error": max(rel_errors),
                "error_at_h_0_1": rel_errors[h_values.index(0.1)] if 0.1 in h_values else None,
                "error_at_h_0_5": rel_errors[h_values.index(0.5)] if 0.5 in h_values else None,
                "error_at_h_1_0": rel_errors[h_values.index(1.0)] if 1.0 in h_values else None,
            }
        return results


def print_benchmark_table() -> None:
    """Prints a formatted markdown table comparing all approximations."""
    h_grid = [0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0]
    results = evaluate_approximations(h_grid)

    header = (
        f"| {'Formula / Method':<48} | "
        + " | ".join([f"h={h:<4}" for h in h_grid])
        + " | Max Rel. Error |"
    )
    sep = "|" + "-" * 50 + "|" + "|".join(["-" * 10 for _ in h_grid]) + "|" + "-" * 16 + "|"
    print(header)
    print(sep)
    for name, data in results.items():
        errs = " | ".join([f"{e:.2e}" for e in data["relative_errors"]])
        print(f"| {name:<48} | {errs} | {data['max_rel_error']:.2e}     |")


if __name__ == "__main__":
    print_benchmark_table()
