"""
Generate publication-quality vector figures for the Hyper-Ramanujan (HR) paper.
Produces: paper/figures/error_comparison_2d.pdf (2-panel relative error & Pareto frontier).
"""

import os
import sys
from pathlib import Path
import mpmath as mp

# Ensure repository root is on sys.path when invoked as `python paper/generate_figures.py`
REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from src.core.exact_geometry import exact_circumference_from_h
from src.approximations.classical import (
    approx_ramanujan_II,
    approx_cantrell_exact_end,
    approx_shankar_a2,
    approx_ayala_raggi_r2_2exp,
)
from src.approximations.novel_candidates import (
    approx_novel_pade_pole_order8,
    approx_novel_minimax_radical,
)
from src.benchmarks.pareto_efficiency import evaluate_pareto_efficiency


def compute_curve_data():
    """Computes 50-digit relative error curves across h in [0.02, 0.9995]."""
    mp.mp.dps = 50
    # Dense grid with extra resolution near h -> 1
    hs = [mp.mpf(i) / 200 for i in range(4, 198)]  # 0.02 to 0.985
    hs += [mp.mpf("0.985") + mp.mpf(i) * mp.mpf("0.001") for i in range(1, 15)]  # up to 0.999
    hs += [mp.mpf("0.9995"), mp.mpf("0.9999")]

    models = [
        ("Ramanujan II (1914)", approx_ramanujan_II),
        ("Cantrell (2004)", approx_cantrell_exact_end),
        ("Shankar $A_2$ (2026)", approx_shankar_a2),
        ("Ayala-Raggi R2/2exp (2026)", approx_ayala_raggi_r2_2exp),
        ("Oertel HR-8P $\\star$ (2026)", approx_novel_pade_pole_order8),
        ("Oertel HR-6M (2026)", approx_novel_minimax_radical),
    ]

    exacts = [exact_circumference_from_h(h) for h in hs]
    h_floats = [float(h) for h in hs]
    curves = {}
    for label, fn in models:
        errs = []
        for h, ex in zip(hs, exacts):
            rel = abs(fn(h) - ex) / ex
            errs.append(max(float(rel), 1e-20))
        curves[label] = errs
    return h_floats, curves


def generate_error_figure(output_path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.family": "serif",
        "mathtext.fontset": "cm",
        "axes.labelsize": 10.5,
        "axes.titlesize": 11,
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "legend.fontsize": 8.5,
        "axes.linewidth": 0.8,
    })

    h_vals, curves = compute_curve_data()
    pareto_rows = evaluate_pareto_efficiency()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.6, 4.35), dpi=300)

    # --- Panel (a): Relative Error Curves ---
    styles = {
        "Ramanujan II (1914)":         dict(color="#6b7280", ls="--", lw=1.4, zorder=2),
        "Cantrell (2004)":             dict(color="#d97706", ls="-.", lw=1.5, zorder=3),
        "Shankar $A_2$ (2026)":        dict(color="#7c3aed", ls=":",  lw=1.8, zorder=4),
        "Ayala-Raggi R2/2exp (2026)":  dict(color="#0891b2", ls="--", lw=1.5, zorder=4),
        "Oertel HR-8P $\\star$ (2026)": dict(color="#1d4ed8", ls="-",  lw=2.2, zorder=6),
        "Oertel HR-6M (2026)":         dict(color="#059669", ls="-",  lw=1.9, zorder=5),
    }

    for label, errs in curves.items():
        ax1.semilogy(h_vals, errs, label=label, **styles[label])

    # Reference thresholds
    ax1.axhline(1e-6, color="#9ca3af", ls=":", lw=0.9, zorder=1)
    ax1.text(0.03, 1.4e-6, "1 ppm ($10^{-6}$)", fontsize=8, color="#4b5563")
    ax1.axhline(2.220446e-16, color="#dc2626", ls=":", lw=0.9, alpha=0.75, zorder=1)
    ax1.text(0.56, 3.6e-16, "IEEE-754 float64 $\\epsilon_{\\mathrm{mach}}$", fontsize=8, color="#b91c1c")

    ax1.set_xlim(0.0, 1.0)
    ax1.set_ylim(1e-18, 2e-3)
    ax1.set_xlabel(r"Shape parameter $h = ((a-b)/(a+b))^2$")
    ax1.set_ylabel(r"Relative error $\varepsilon_{\mathrm{rel}}(h)$")
    ax1.set_title("(a) 2D Relative Error Across Full Eccentricity Range")
    ax1.grid(True, which="major", ls=":", lw=0.5, alpha=0.55)
    ax1.legend(loc="upper left", frameon=True, framealpha=0.92, edgecolor="#d1d5db")

    # --- Panel (b): OpCost vs Global Max Error (ppm) ---
    short_names = {
        "Hudson [2/1] (1917)": ("Hudson [2/1]", (4, 5)),
        "Jacobsen-Waadeland [2/2] (1985)": ("JW [2/2]", (4, 5)),
        "Ramanujan II (1914)": ("Ram II", (4, 5)),
        "Shankar A1 (2026)": ("Shankar $A_1$", (4, 6)),
        "Cantrell (2004)": ("Cantrell", (-46, -12)),
        "Shankar A2 (2026)": ("Shankar $A_2$", (5, 4)),
        "Ayala-Raggi & Rendon-Marin (2026, R2/2exp)": ("R2/2exp\n(3 exp)", (-44, 6)),
        "HR-7-Int (Order-7, pure integers)": ("HR-7-Int", (5, 4)),
        "HR-6B (Order-6 + exact h=1)": ("HR-6B", (5, -10)),
        "HR-8B (Order-8 + exact h=1)": ("HR-8B", (5, 3)),
        "HR-8P (Order-8, -13/14 pole + exact h=1)": ("HR-8P $\\star$", (6, 4)),
        "HR-6M (Order-6, Minimax [3/1] + exact h=1)": ("HR-6M", (6, -6)),
    }

    # Plot Pareto envelope on (OpCost, Max Err ppm)
    pts_envelope = []
    for row in pareto_rows:
        name = row["name"]
        if name not in short_names:
            continue
        cost = row["op_cost"]
        err_ppm = row["max_global_ppm"]
        label, offset = short_names[name]
        is_hr = name.startswith("HR-")
        is_champ = "HR-8P" in name
        is_pareto = row["pareto_frontier"]

        if is_champ:
            ax2.scatter([cost], [err_ppm], s=140, marker="*", color="#1d4ed8", edgecolor="black", lw=0.6, zorder=7)
            ax2.annotate(label, (cost, err_ppm), xytext=offset, textcoords="offset points",
                         fontsize=8.5, fontweight="bold", color="#1d4ed8")
        elif is_hr:
            c_col = "#059669" if "HR-6M" in name else "#2563eb"
            ax2.scatter([cost], [err_ppm], s=48, marker="o", color=c_col, edgecolor="black", lw=0.5, zorder=6)
            ax2.annotate(label, (cost, err_ppm), xytext=offset, textcoords="offset points",
                         fontsize=8, fontweight="bold", color=c_col)
        elif is_pareto:
            ax2.scatter([cost], [err_ppm], s=42, marker="s", color="#6b7280", edgecolor="black", lw=0.5, zorder=5)
            ax2.annotate(label, (cost, err_ppm), xytext=offset, textcoords="offset points",
                         fontsize=8, color="#374151")
        else:
            ax2.scatter([cost], [err_ppm], s=46, marker="o", facecolors="none", edgecolors="#dc2626", lw=1.3, zorder=5)
            ax2.annotate(label, (cost, err_ppm), xytext=offset, textcoords="offset points",
                         fontsize=8, color="#b91c1c")

        if is_pareto:
            pts_envelope.append((cost, err_ppm))

    # Draw visual lower-left Pareto hull in (Cost, Max Err)
    pts_envelope.sort(key=lambda t: t[0])
    hull = []
    best_err = float("inf")
    for c, e in pts_envelope:
        if e < best_err:
            hull.append((c, e))
            best_err = e
    if len(hull) >= 2:
        ax2.plot([p[0] for p in hull], [p[1] for p in hull], color="#9ca3af", ls="--", lw=1.0, zorder=2, label="Global-error Pareto frontier")

    ax2.set_yscale("log")
    ax2.set_xlim(5, 96)
    ax2.set_ylim(0.15, 4500)
    ax2.set_xlabel(r"Hardware Operational Cost ($\mathrm{OpCost}$)")
    ax2.set_ylabel(r"Global Max Relative Error (ppm)")
    ax2.set_title("(b) Accuracy vs. Computational Cost (Lower-Left Is Better)")
    ax2.grid(True, which="major", ls=":", lw=0.5, alpha=0.55)

    # Compact legend for markers
    from matplotlib.lines import Line2D
    legend_handles = [
        Line2D([0], [0], marker="*", color="w", markerfacecolor="#1d4ed8", markeredgecolor="k", markersize=10, label="HR-8P (Pareto Champion)"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="#2563eb", markeredgecolor="k", markersize=6.5, label="Other HR Formulas (2026)"),
        Line2D([0], [0], marker="s", color="w", markerfacecolor="#6b7280", markeredgecolor="k", markersize=6, label="Classical Pareto-Optimal"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor="none", markeredgecolor="#dc2626", markeredgewidth=1.3, markersize=6.5, label="Strictly Dominated"),
    ]
    ax2.legend(handles=legend_handles, loc="upper right", frameon=True, framealpha=0.92, edgecolor="#d1d5db")

    fig.tight_layout(pad=1.3)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, format="pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] Saved vector figure to {output_path}")


if __name__ == "__main__":
    out = REPO_ROOT / "paper" / "figures" / "error_comparison_2d.pdf"
    if "--dry-run" in sys.argv:
        h_vals, curves = compute_curve_data()
        print(f"[DRY-RUN OK] Computed {len(h_vals)} grid points for {len(curves)} models.")
    else:
        generate_error_figure(out)
