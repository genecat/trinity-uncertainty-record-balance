"""Phase 4B robustness test for concentration-dispersion balance peaks.

This script reuses the Phase 4 record families and tests whether the observed
high-balance regimes persist under multiple alternative balance metrics rather
than only under the original

    S_balance = 1 - |C - D| / (C + D + eps)

definition. This is exploratory numerical research on a toy record-space model.
It does not establish a new physical law.
"""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt


X_MIN = -50.0
X_MAX = 50.0
N = 16384
HBAR = 1.0
TOLERANCE = 1e-8
EPS = 1e-12
SWEEP_POINTS = 41


def normalize_amplitude(amplitude: np.ndarray, dx: float) -> np.ndarray:
    """Normalize a complex amplitude so that integral |R(x)|^2 dx = 1."""
    norm = np.sqrt(np.sum(np.abs(amplitude) ** 2) * dx)
    if norm == 0.0:
        raise ValueError("Amplitude cannot be identically zero.")
    return amplitude / norm


def expectation_from_pdf(grid: np.ndarray, pdf: np.ndarray, spacing: float) -> tuple[float, float]:
    """Return the mean and variance for a normalized PDF on a uniform grid."""
    mean = np.sum(grid * pdf) * spacing
    variance = np.sum((grid - mean) ** 2 * pdf) * spacing
    variance = max(variance, 0.0)
    return mean, variance


def fourier_transform_record(record_x: np.ndarray, dx: float) -> tuple[np.ndarray, np.ndarray]:
    """Compute the continuous-convention Fourier transform and angular k grid."""
    k_grid = 2.0 * np.pi * np.fft.fftshift(np.fft.fftfreq(record_x.size, d=dx))
    record_k = (
        np.fft.fftshift(np.fft.fft(np.fft.ifftshift(record_x))) * dx / np.sqrt(2.0 * np.pi)
    )
    return k_grid, record_k


def gaussian_profile(x: np.ndarray, sigma: float, center: float = 0.0) -> np.ndarray:
    return np.exp(-((x - center) ** 2) / (4.0 * sigma**2))


def smooth_horizon_window(x: np.ndarray, boundary_scale: float, softness: float) -> np.ndarray:
    scaled = np.clip((np.abs(x) - boundary_scale) / softness, -700.0, 700.0)
    return 1.0 / (1.0 + np.exp(scaled))


def boundary_locked_gaussian_profile(
    x: np.ndarray,
    sigma: float,
    boundary_scale: float,
    softness: float,
) -> np.ndarray:
    return gaussian_profile(x, sigma=sigma) * smooth_horizon_window(x, boundary_scale, softness)


def overlap_record(
    x: np.ndarray,
    sigma: float,
    separation: float,
    left_weight: float,
    right_weight: float,
) -> np.ndarray:
    left = left_weight * gaussian_profile(x, sigma=sigma, center=-separation / 2.0)
    right = right_weight * gaussian_profile(x, sigma=sigma, center=separation / 2.0)
    return left + right


def oscillatory_selector_record(
    x: np.ndarray,
    sigma: float,
    beta: float,
    omega: float,
) -> np.ndarray:
    envelope = gaussian_profile(x, sigma=sigma)
    return envelope * np.exp(1j * beta * np.sin(omega * x))


def normalized_shannon_entropy(pdf: np.ndarray, spacing: float) -> float:
    """Compute grid-normalized Shannon entropy for a continuous PDF sample."""
    probabilities = pdf * spacing
    probabilities /= np.sum(probabilities)
    entropy = -np.sum(probabilities * np.log(probabilities + EPS))
    return float(entropy / np.log(pdf.size))


def evaluate_record(
    family: str,
    parameter_name: str,
    parameter_value: float,
    amplitude_x: np.ndarray,
    x_grid: np.ndarray,
    dx: float,
) -> dict[str, float | str]:
    record_x = normalize_amplitude(amplitude_x.astype(np.complex128), dx)
    pdf_x = np.abs(record_x) ** 2

    _, var_x = expectation_from_pdf(x_grid, pdf_x, dx)
    delta_x = np.sqrt(var_x)

    k_grid, record_k = fourier_transform_record(record_x, dx)
    dk = float(k_grid[1] - k_grid[0])
    pdf_k = np.abs(record_k) ** 2
    pdf_k /= np.sum(pdf_k) * dk

    _, var_k = expectation_from_pdf(k_grid, pdf_k, dk)
    delta_k = np.sqrt(var_k)
    uncertainty_product = delta_x * (HBAR * delta_k)
    ratio_to_hbar_over_2 = uncertainty_product / (HBAR / 2.0)

    concentration = 1.0 / (delta_x + EPS)
    dispersion = delta_k

    s_balance_original = 1.0 - abs(concentration - dispersion) / (concentration + dispersion + EPS)
    s_balance_log = float(np.exp(-abs(np.log((concentration + EPS) / (dispersion + EPS)))))
    s_balance_product_distance = 1.0 / ratio_to_hbar_over_2

    hx_norm = normalized_shannon_entropy(pdf_x, dx)
    hk_norm = normalized_shannon_entropy(pdf_k, dk)
    s_balance_entropy = 1.0 - abs(hx_norm - hk_norm)

    s_combined = float(np.mean([s_balance_original, s_balance_log, s_balance_entropy]))
    passed = uncertainty_product >= (HBAR / 2.0 - TOLERANCE)

    return {
        "family": family,
        "parameter_name": parameter_name,
        "parameter_value": float(parameter_value),
        "delta_x": delta_x,
        "delta_k": delta_k,
        "uncertainty_product": uncertainty_product,
        "ratio_to_hbar_over_2": ratio_to_hbar_over_2,
        "S_balance_original": s_balance_original,
        "S_balance_log": s_balance_log,
        "S_balance_product_distance": s_balance_product_distance,
        "Hx_norm": hx_norm,
        "Hk_norm": hk_norm,
        "S_balance_entropy": s_balance_entropy,
        "S_combined": s_combined,
        "pass_fail": "PASS" if passed else "FAIL",
    }


def run_sweeps(x_grid: np.ndarray, dx: float) -> list[dict[str, float | str]]:
    results: list[dict[str, float | str]] = []

    gaussian_sigmas = np.linspace(0.5, 12.0, SWEEP_POINTS)
    for sigma in gaussian_sigmas:
        results.append(
            evaluate_record(
                family="gaussian_width_sweep",
                parameter_name="sigma",
                parameter_value=sigma,
                amplitude_x=gaussian_profile(x_grid, sigma=sigma),
                x_grid=x_grid,
                dx=dx,
            )
        )

    softness_values = np.linspace(0.05, 5.0, SWEEP_POINTS)
    for softness in softness_values:
        results.append(
            evaluate_record(
                family="boundary_softness_sweep",
                parameter_name="softness",
                parameter_value=softness,
                amplitude_x=boundary_locked_gaussian_profile(
                    x_grid,
                    sigma=4.0,
                    boundary_scale=8.0,
                    softness=softness,
                ),
                x_grid=x_grid,
                dx=dx,
            )
        )

    separations = np.linspace(0.0, 20.0, SWEEP_POINTS)
    for separation in separations:
        results.append(
            evaluate_record(
                family="overlap_separation_sweep",
                parameter_name="separation_d",
                parameter_value=separation,
                amplitude_x=overlap_record(
                    x_grid,
                    sigma=2.0,
                    separation=separation,
                    left_weight=1.0,
                    right_weight=1.0,
                ),
                x_grid=x_grid,
                dx=dx,
            )
        )

    amplitude_ratios = np.linspace(0.1, 1.0, SWEEP_POINTS)
    for ratio in amplitude_ratios:
        results.append(
            evaluate_record(
                family="asymmetric_overlap_sweep",
                parameter_name="amplitude_ratio_r",
                parameter_value=ratio,
                amplitude_x=overlap_record(
                    x_grid,
                    sigma=2.0,
                    separation=8.0,
                    left_weight=1.0,
                    right_weight=ratio,
                ),
                x_grid=x_grid,
                dx=dx,
            )
        )

    k0_values = np.linspace(0.0, 5.0, SWEEP_POINTS)
    for k0 in k0_values:
        results.append(
            evaluate_record(
                family="selector_phase_sweep",
                parameter_name="k0",
                parameter_value=k0,
                amplitude_x=gaussian_profile(x_grid, sigma=4.0) * np.exp(1j * k0 * x_grid),
                x_grid=x_grid,
                dx=dx,
            )
        )

    beta_values = np.linspace(0.0, 5.0, SWEEP_POINTS)
    for beta in beta_values:
        results.append(
            evaluate_record(
                family="oscillatory_selector_sweep",
                parameter_name="beta",
                parameter_value=beta,
                amplitude_x=oscillatory_selector_record(
                    x_grid,
                    sigma=4.0,
                    beta=beta,
                    omega=1.5,
                ),
                x_grid=x_grid,
                dx=dx,
            )
        )

    return results


def save_results_csv(results: list[dict[str, float | str]], output_path: Path) -> None:
    fieldnames = [
        "family",
        "parameter_name",
        "parameter_value",
        "delta_x",
        "delta_k",
        "uncertainty_product",
        "ratio_to_hbar_over_2",
        "S_balance_original",
        "S_balance_log",
        "S_balance_product_distance",
        "Hx_norm",
        "Hk_norm",
        "S_balance_entropy",
        "S_combined",
        "pass_fail",
    ]
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            writer.writerow(result)


def group_by_family(results: list[dict[str, float | str]]) -> dict[str, list[dict[str, float | str]]]:
    grouped: dict[str, list[dict[str, float | str]]] = defaultdict(list)
    for result in results:
        grouped[str(result["family"])].append(result)
    return dict(grouped)


def max_row(rows: list[dict[str, float | str]], key: str) -> dict[str, float | str]:
    return max(rows, key=lambda row: float(row[key]))


def min_row(rows: list[dict[str, float | str]], key: str) -> dict[str, float | str]:
    return min(rows, key=lambda row: float(row[key]))


def compute_family_summary(results: list[dict[str, float | str]]) -> list[dict[str, float | str]]:
    summaries: list[dict[str, float | str]] = []
    grouped = group_by_family(results)

    for family in sorted(grouped):
        rows = grouped[family]
        best_original = max_row(rows, "S_balance_original")
        best_log = max_row(rows, "S_balance_log")
        best_entropy = max_row(rows, "S_balance_entropy")
        best_combined = max_row(rows, "S_combined")
        best_uncertainty = min_row(rows, "uncertainty_product")

        summaries.append(
            {
                "family": family,
                "max_S_balance_original": float(best_original["S_balance_original"]),
                "parameter_at_max_S_balance_original": float(best_original["parameter_value"]),
                "max_S_balance_log": float(best_log["S_balance_log"]),
                "parameter_at_max_S_balance_log": float(best_log["parameter_value"]),
                "max_S_balance_entropy": float(best_entropy["S_balance_entropy"]),
                "parameter_at_max_S_balance_entropy": float(best_entropy["parameter_value"]),
                "max_S_combined": float(best_combined["S_combined"]),
                "parameter_at_max_S_combined": float(best_combined["parameter_value"]),
                "min_uncertainty_product": float(best_uncertainty["uncertainty_product"]),
                "parameter_at_min_uncertainty_product": float(best_uncertainty["parameter_value"]),
            }
        )

    return summaries


def print_table(headers: list[str], rows: list[list[str]]) -> None:
    widths = [
        max(len(header), max(len(row[index]) for row in rows))
        for index, header in enumerate(headers)
    ]
    header_line = " | ".join(header.ljust(widths[index]) for index, header in enumerate(headers))
    separator = "-+-".join("-" * width for width in widths)
    print(header_line)
    print(separator)
    for row in rows:
        print(" | ".join(value.ljust(widths[index]) for index, value in enumerate(row)))


def print_family_summary(summaries: list[dict[str, float | str]]) -> None:
    headers = [
        "family",
        "max_S_balance_original",
        "parameter",
        "max_S_balance_log",
        "parameter",
        "max_S_balance_entropy",
        "parameter",
        "max_S_combined",
        "parameter",
        "min_uncertainty_product",
        "parameter",
    ]
    rows = []
    for summary in summaries:
        rows.append(
            [
                str(summary["family"]),
                f'{summary["max_S_balance_original"]:.6f}',
                f'{summary["parameter_at_max_S_balance_original"]:.6f}',
                f'{summary["max_S_balance_log"]:.6f}',
                f'{summary["parameter_at_max_S_balance_log"]:.6f}',
                f'{summary["max_S_balance_entropy"]:.6f}',
                f'{summary["parameter_at_max_S_balance_entropy"]:.6f}',
                f'{summary["max_S_combined"]:.6f}',
                f'{summary["parameter_at_max_S_combined"]:.6f}',
                f'{summary["min_uncertainty_product"]:.6f}',
                f'{summary["parameter_at_min_uncertainty_product"]:.6f}',
            ]
        )
    print_table(headers, rows)


def print_top_combined_records(results: list[dict[str, float | str]], top_n: int = 10) -> None:
    headers = [
        "family",
        "parameter_name",
        "parameter_value",
        "S_combined",
        "S_balance_original",
        "S_balance_log",
        "S_balance_entropy",
        "S_balance_product_distance",
        "uncertainty_product",
        "pass_fail",
    ]
    sorted_rows = sorted(results, key=lambda row: (-float(row["S_combined"]), float(row["uncertainty_product"])))
    rows = []
    for result in sorted_rows[:top_n]:
        rows.append(
            [
                str(result["family"]),
                str(result["parameter_name"]),
                f'{result["parameter_value"]:.6f}',
                f'{result["S_combined"]:.6f}',
                f'{result["S_balance_original"]:.6f}',
                f'{result["S_balance_log"]:.6f}',
                f'{result["S_balance_entropy"]:.6f}',
                f'{result["S_balance_product_distance"]:.6f}',
                f'{result["uncertainty_product"]:.6f}',
                str(result["pass_fail"]),
            ]
        )
    print_table(headers, rows)


def plot_best_combined_by_family(grouped: dict[str, list[dict[str, float | str]]], output_path: Path) -> None:
    family_names = sorted(grouped)
    best_rows = [max_row(grouped[family], "S_combined") for family in family_names]
    values = [float(row["S_combined"]) for row in best_rows]

    fig, ax = plt.subplots(figsize=(11.0, 5.8))
    bars = ax.bar(family_names, values, color="#2c7fb8")
    ax.set_ylabel("Best S_combined in family")
    ax.set_title("Robustness: best combined balance by family")
    ax.tick_params(axis="x", rotation=25)

    for bar, value, row in zip(bars, values, best_rows):
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            value,
            f"{value:.3f}\n@{float(row['parameter_value']):.2f}",
            ha="center",
            va="bottom",
            fontsize=8,
        )

    fig.tight_layout()
    fig.savefig(output_path, dpi=160)
    plt.close(fig)


def plot_scatter(
    grouped: dict[str, list[dict[str, float | str]]],
    x_key: str,
    y_key: str,
    xlabel: str,
    ylabel: str,
    title: str,
    output_path: Path,
) -> None:
    fig, ax = plt.subplots(figsize=(11.5, 6.8))
    for family in sorted(grouped):
        rows = grouped[family]
        x_values = [float(row[x_key]) for row in rows]
        y_values = [float(row[y_key]) for row in rows]
        ax.scatter(x_values, y_values, s=24, alpha=0.8, label=family)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_path, dpi=160)
    plt.close(fig)


def main() -> None:
    project_root = Path(__file__).resolve().parents[1]
    outputs_dir = project_root / "outputs"
    outputs_dir.mkdir(exist_ok=True)

    x_grid = np.linspace(X_MIN, X_MAX, N, endpoint=False)
    dx = float(x_grid[1] - x_grid[0])

    results = run_sweeps(x_grid, dx)
    summaries = compute_family_summary(results)
    grouped = group_by_family(results)

    save_results_csv(results, outputs_dir / "concentration_dispersion_robustness_results.csv")
    plot_best_combined_by_family(grouped, outputs_dir / "robustness_combined_balance_by_family.png")
    plot_scatter(
        grouped=grouped,
        x_key="S_balance_original",
        y_key="S_balance_log",
        xlabel="S_balance_original",
        ylabel="S_balance_log",
        title="Original vs log balance metrics",
        output_path=outputs_dir / "robustness_original_vs_log_balance.png",
    )
    plot_scatter(
        grouped=grouped,
        x_key="S_balance_entropy",
        y_key="S_combined",
        xlabel="S_balance_entropy",
        ylabel="S_combined",
        title="Entropy balance vs combined balance",
        output_path=outputs_dir / "robustness_entropy_vs_combined_balance.png",
    )
    plot_scatter(
        grouped=grouped,
        x_key="uncertainty_product",
        y_key="S_combined",
        xlabel="Uncertainty product",
        ylabel="S_combined",
        title="Uncertainty product vs combined balance",
        output_path=outputs_dir / "robustness_uncertainty_vs_combined_balance.png",
    )

    print("Family summary")
    print_family_summary(summaries)
    print()
    print("Top 10 records by S_combined")
    print_top_combined_records(results, top_n=10)


if __name__ == "__main__":
    main()
