"""Reproduce the analytic-circle and numerical-ellipse validation figure.

Run from any working directory:
    python /path/to/sinkhorn_curvature/code/make_figures.py

The ellipse inputs are produced by check_ellipse_bias.py. Circle eigenvalues
use the exact analytic Bessel formula, evaluated in double precision with
exponentially scaled modified Bessel functions. This is deterministic
mathematical validation, not an empirical dataset.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import platform
import tempfile

import numpy as np
import scipy
from scipy.special import ive
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "sinkhorn-matplotlib-cache"))
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


PLUM = "#76375F"
COPPER = "#BA6B3D"
TEAL = "#177D84"
CHARCOAL = "#3C3940"


def log_slope(epsilon: np.ndarray, error: np.ndarray) -> float:
    """Least-squares descriptive log-log slope; no sampling interpretation."""
    return float(np.polyfit(np.log(epsilon), np.log(error), 1)[0])


def circle_check() -> dict:
    radius = 1.0
    epsilon = np.geomspace(1e-4, 0.25, 100)
    kappa = 2 * radius**2 / epsilon
    denominator = ive(0, kappa)
    modes = []
    for mode in (1, 2, 3):
        p = ive(mode, kappa) / denominator
        # exp(-kappa) factors cancel in the ratio. The sinh expression
        # avoids a direct subtraction of 1-p**2.
        z = -np.log(p)
        metric = epsilon / (4 * np.sinh(z))
        limit = radius**2 / mode**2
        coefficient = (metric - limit) / epsilon
        rescaled_coefficient = mode**2 * coefficient
        linear_prediction = limit - epsilon / (4 * mode**2)
        remainder = metric - linear_prediction
        slope_mask = (epsilon >= 1e-3) & (epsilon <= 1e-2)
        slope = log_slope(epsilon[slope_mask], np.abs(remainder[slope_mask]))
        if not np.all((p > 0) & (p < 1) & (metric > 0)):
            raise ArithmeticError("Circle eigenvalues failed positivity bounds.")
        if abs(rescaled_coefficient[0] + 0.25) > 1e-3:
            raise ArithmeticError("Circle first-order coefficient failed its limit check.")
        if not (1.8 < slope < 2.2):
            raise ArithmeticError("Circle second-order remainder failed its order check.")
        modes.append({
            "mode": mode,
            "markov_eigenvalue": p.tolist(),
            "metric_eigenvalue": metric.tolist(),
            "limit_metric_eigenvalue": limit,
            "predicted_first_order_coefficient": -1 / (4 * mode**2),
            "measured_first_order_coefficient": coefficient.tolist(),
            "rescaled_first_order_coefficient": rescaled_coefficient.tolist(),
            "second_order_remainder": remainder.tolist(),
            "descriptive_second_order_log_slope": slope,
        })
    return {
        "radius": radius,
        "epsilon": epsilon.tolist(),
        "modes": modes,
        "formula": "p_k=I_k(2R^2/epsilon)/I_0(2R^2/epsilon); r_epsilon,k=epsilon*p_k/(2*(1-p_k^2))",
        "prediction": "r_epsilon,k=R^2/k^2-epsilon/(4k^2)+O(epsilon^2)",
        "evaluation": "Analytic Bessel formula evaluated in double precision with scipy.special.ive; not exact arithmetic.",
        "slope_epsilon_range": [1e-3, 1e-2],
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "matplotlib": matplotlib.__version__,
    }


def ellipse_check(results_dir: Path) -> dict:
    base = json.loads((results_dir / "ellipse_bias_results.json").read_text())
    refined = json.loads((results_dir / "ellipse_bias_refinement_1536.json").read_text())
    coarse = json.loads((results_dir / "ellipse_bias_refinement_768.json").read_text())
    rows = sorted(base["results"] + refined["results"], key=lambda row: row["epsilon"])
    epsilon = np.array([row["epsilon"] for row in rows])
    original_error = np.array([abs(row["metric"] - base["g0"]) for row in rows])
    corrected_error = np.array([abs(row["feasible_corrected_metric"] - base["g0"]) for row in rows])
    low = epsilon <= 0.02
    original_slope = log_slope(epsilon[low], original_error[low])
    corrected_slope = log_slope(epsilon[low], corrected_error[low])
    grid_error = max(abs(a["metric"]-b["metric"])
                     for a, b in zip(coarse["results"], refined["results"], strict=True))
    if not np.all((original_error > 0) & (corrected_error > 0) & (corrected_error < original_error)):
        raise ArithmeticError("Ellipse correction/error plotting check failed.")
    if grid_error > 1e-10:
        raise ArithmeticError("Ellipse refinement discrepancy is too large.")
    return {
        "epsilon": epsilon.tolist(),
        "original_absolute_error": original_error.tolist(),
        "corrected_absolute_error": corrected_error.tolist(),
        "limiting_metric": base["g0"],
        "descriptive_original_log_slope": original_slope,
        "descriptive_corrected_log_slope": corrected_slope,
        "slope_epsilon_range": [float(epsilon[low].min()), float(epsilon[low].max())],
        "maximum_768_vs_1536_grid_discrepancy": grid_error,
        "main_grid_size": base["n"],
        "small_epsilon_grid_size": refined["n"],
    }


def make_figure(circle: dict, ellipse: dict, figures_dir: Path) -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "legend.fontsize": 9,
        "axes.labelcolor": CHARCOAL,
        "text.color": CHARCOAL,
        "xtick.color": CHARCOAL,
        "ytick.color": CHARCOAL,
        "axes.edgecolor": CHARCOAL,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.facecolor": "white",
    })
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 4.65), layout="constrained")
    ax = axes[0]
    epsilon = np.array(circle["epsilon"])
    styles = [(PLUM, "-"), (COPPER, "--"), (TEAL, "-.")]
    for result, (colour, linestyle) in zip(circle["modes"], styles, strict=True):
        ax.semilogx(epsilon, result["rescaled_first_order_coefficient"],
                    color=colour, linestyle=linestyle, linewidth=2,
                    label=rf"Mode $k={result['mode']}$")
    ax.axhline(-0.25, color=CHARCOAL, linestyle=":", linewidth=1.4,
               label=r"Predicted limit $-1/4$")
    ax.set_title("(a) Circle: first-order coefficient", loc="left", pad=12)
    ax.set_xlabel(r"Regularisation parameter $\varepsilon$")
    ax.set_ylabel(r"$k^2(b_{\varepsilon,k}-R^2/k^2)/\varepsilon$")
    ax.set_xlim(epsilon.min(), epsilon.max())
    ax.legend(frameon=False, loc="lower left")
    ax.text(0.02, 0.88, r"Unit radius; analytic Bessel formula", transform=ax.transAxes,
            va="top", fontsize=9, color=CHARCOAL)

    ax = axes[1]
    epsilon = np.array(ellipse["epsilon"])
    ax.loglog(epsilon, ellipse["original_absolute_error"],
              color=PLUM, linewidth=2, marker="o", markersize=5,
              label=f"Original metric · slope {ellipse['descriptive_original_log_slope']:.2f}")
    ax.loglog(epsilon, ellipse["corrected_absolute_error"],
              color=COPPER, linewidth=2, marker="s", markersize=5,
              label=f"Positive correction · slope {ellipse['descriptive_corrected_log_slope']:.2f}")
    ax.set_title("(b) Ellipse: error after correction", loc="left", pad=12)
    ax.set_xlabel(r"Regularisation parameter $\varepsilon$")
    ax.set_ylabel("Absolute metric error")
    ax.set_xlim(epsilon.min()/1.15, epsilon.max()*1.15)
    ax.legend(frameon=False, loc="upper left")
    ax.text(0.98, 0.03, "Nonuniform density; varying curvature\nSlopes fitted for ε ≤ 0.02",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=9, color=CHARCOAL)
    for ax in axes:
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="both", which="major", color="#D9D3D7", linewidth=0.6, alpha=0.65)
        ax.set_axisbelow(True)
    figures_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(figures_dir / "sinkhorn_curvature_validation.pdf", bbox_inches="tight")
    fig.savefig(figures_dir / "sinkhorn_curvature_validation.png", dpi=240, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-dir", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    results_dir = args.project_dir / "results"
    figures_dir = args.project_dir / "figures"
    circle = circle_check()
    ellipse = ellipse_check(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    (results_dir / "circle_bessel_results.json").write_text(json.dumps(circle, indent=2)+"\n")
    (results_dir / "figure_validation_summary.json").write_text(json.dumps(ellipse, indent=2)+"\n")
    make_figure(circle, ellipse, figures_dir)
    print(json.dumps({
        "circle_smallest_epsilon": circle["epsilon"][0],
        "circle_scaled_coefficients_at_smallest_epsilon": [
            item["rescaled_first_order_coefficient"][0] for item in circle["modes"]],
        "circle_second_order_slopes": [item["descriptive_second_order_log_slope"] for item in circle["modes"]],
        "ellipse_original_slope": ellipse["descriptive_original_log_slope"],
        "ellipse_corrected_slope": ellipse["descriptive_corrected_log_slope"],
        "ellipse_grid_discrepancy": ellipse["maximum_768_vs_1536_grid_discrepancy"],
        "figures_directory": str(figures_dir),
    }, indent=2))


if __name__ == "__main__":
    main()
