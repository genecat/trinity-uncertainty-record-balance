"""Phase 2 Trinity/STO-inspired record-locking uncertainty test.

Interpretive mapping used in this toy model:
- R(x) is the record-localization face.
- R(k) is the action-gradient / momentum face.
- hbar remains the standard conversion factor p = hbar * k.

This script extends the Phase 1 uncertainty check with Trinity-inspired record
profiles that model boundary locking, scalar-clock expansion windows, selector
phase structure, and overlapping records. The goal is to test whether these
profiles preserve the standard Fourier uncertainty floor while generating
distinct concentration-dispersion behavior. It does not claim to derive hbar
from first principles.
"""

from __future__ import annotations

import csv
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


def normalize_amplitude(amplitude: np.ndarray, dx: float) -> np.ndarray:
    """Normalize a complex amplitude so integral |R(x)|^2 dx = 1."""
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
    """Smooth locking window used to model a softened record horizon."""
    scaled = np.clip((np.abs(x) - boundary_scale) / softness, -700.0, 700.0)
    return 1.0 / (1.0 + np.exp(scaled))


def boundary_locked_gaussian_profile(
    x: np.ndarray,
    sigma: float,
    boundary_scale: float,
    softness: float,
) -> np.ndarray:
    return gaussian_profile(x, sigma=sigma) * smooth_horizon_window(x, boundary_scale, softness)


def scalar_clock_expansion_record(
    x: np.ndarray,
    sigma_core: float,
    tail_strength: float,
    tail_scale: float,
) -> np.ndarray:
    core = gaussian_profile(x, sigma=sigma_core)
    # The weak tail is a toy expansion channel appended to a localized core.
    tail = tail_strength * np.exp(-np.abs(x) / tail_scale)
    return core + tail


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


def build_profiles(x: np.ndarray) -> list[tuple[str, np.ndarray]]:
    boundary_locked = boundary_locked_gaussian_profile(x, sigma=4.0, boundary_scale=10.0, softness=1.5)
    sharp_boundary_locked = boundary_locked_gaussian_profile(x, sigma=4.0, boundary_scale=10.0, softness=0.25)

    return [
        ("boundary_locked_gaussian", boundary_locked),
        ("sharp_boundary_locked_gaussian", sharp_boundary_locked),
        (
            "scalar_clock_expansion_record",
            scalar_clock_expansion_record(x, sigma_core=3.0, tail_strength=0.18, tail_scale=12.0),
        ),
        (
            "selector_phase_twisted_record",
            boundary_locked * np.exp(1j * 1.8 * x),
        ),
        ("two_time_overlap_record", overlap_record(x, sigma=1.5, separation=10.0, left_weight=1.0, right_weight=1.0)),
        ("asymmetric_overlap_record", overlap_record(x, sigma=1.5, separation=10.0, left_weight=1.0, right_weight=0.45)),
        ("oscillatory_selector_record", oscillatory_selector_record(x, sigma=4.0, beta=0.9, omega=0.7)),
    ]


def evaluate_profile(name: str, amplitude_x: np.ndarray, x_grid: np.ndarray, dx: float) -> dict[str, float | str]:
    record_x = normalize_amplitude(amplitude_x.astype(np.complex128), dx)
    pdf_x = np.abs(record_x) ** 2

    mean_x, var_x = expectation_from_pdf(x_grid, pdf_x, dx)
    delta_x = np.sqrt(var_x)

    k_grid, record_k = fourier_transform_record(record_x, dx)
    dk = float(k_grid[1] - k_grid[0])
    pdf_k = np.abs(record_k) ** 2

    # Renormalize after finite-grid FFT sampling so the k-space PDF integrates to 1.
    pdf_k /= np.sum(pdf_k) * dk

    mean_k, var_k = expectation_from_pdf(k_grid, pdf_k, dk)
    delta_k = np.sqrt(var_k)
    delta_p = HBAR * delta_k
    uncertainty_product = delta_x * delta_p
    ratio_to_hbar_over_2 = uncertainty_product / (HBAR / 2.0)

    # Toy balance diagnostic:
    # C measures localization concentration and D measures action-gradient dispersion.
    concentration = 1.0 / (delta_x + EPS)
    dispersion = delta_k
    s_balance = 1.0 - abs(concentration - dispersion) / (concentration + dispersion + EPS)
    passed = uncertainty_product >= (HBAR / 2.0 - TOLERANCE)

    return {
        "profile_name": name,
        "mean_x": mean_x,
        "mean_k": mean_k,
        "delta_x": delta_x,
        "delta_k": delta_k,
        "delta_p": delta_p,
        "uncertainty_product": uncertainty_product,
        "ratio_to_hbar_over_2": ratio_to_hbar_over_2,
        "C": concentration,
        "D": dispersion,
        "S_balance": s_balance,
        "pass_fail": "PASS" if passed else "FAIL",
    }


def print_results_table(results: list[dict[str, float | str]]) -> None:
    headers = [
        "profile_name",
        "mean_x",
        "mean_k",
        "delta_x",
        "delta_k",
        "uncertainty_product",
        "ratio_to_hbar_over_2",
        "C",
        "D",
        "S_balance",
        "pass_fail",
    ]
    rows = []
    for result in results:
        rows.append(
            [
                str(result["profile_name"]),
                f'{result["mean_x"]:.6f}',
                f'{result["mean_k"]:.6f}',
                f'{result["delta_x"]:.6f}',
                f'{result["delta_k"]:.6f}',
                f'{result["uncertainty_product"]:.6f}',
                f'{result["ratio_to_hbar_over_2"]:.6f}',
                f'{result["C"]:.6f}',
                f'{result["D"]:.6f}',
                f'{result["S_balance"]:.6f}',
                str(result["pass_fail"]),
            ]
        )

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


def save_results_csv(results: list[dict[str, float | str]], output_path: Path) -> None:
    fieldnames = [
        "profile_name",
        "mean_x",
        "mean_k",
        "delta_x",
        "delta_k",
        "delta_p",
        "uncertainty_product",
        "ratio_to_hbar_over_2",
        "C",
        "D",
        "S_balance",
        "pass_fail",
    ]
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            writer.writerow(result)


def save_bar_plot(
    results: list[dict[str, float | str]],
    key: str,
    ylabel: str,
    title: str,
    output_path: Path,
    reference_line: float | None = None,
    reference_label: str | None = None,
) -> None:
    profile_names = [str(result["profile_name"]) for result in results]
    values = [float(result[key]) for result in results]

    fig, ax = plt.subplots(figsize=(12, 6))
    bars = ax.bar(profile_names, values, color="#2c7fb8")
    if reference_line is not None:
        ax.axhline(reference_line, color="#d95f0e", linestyle="--", linewidth=2, label=reference_label)
        ax.legend()

    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.tick_params(axis="x", rotation=25)

    for bar, value in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            value,
            f"{value:.3f}",
            ha="center",
            va="bottom",
            fontsize=9,
        )

    fig.tight_layout()
    fig.savefig(output_path, dpi=160)
    plt.close(fig)


def main() -> None:
    project_root = Path(__file__).resolve().parents[1]
    outputs_dir = project_root / "outputs"
    outputs_dir.mkdir(exist_ok=True)

    x_grid = np.linspace(X_MIN, X_MAX, N, endpoint=False)
    dx = float(x_grid[1] - x_grid[0])

    results = []
    for name, amplitude_x in build_profiles(x_grid):
        results.append(evaluate_profile(name, amplitude_x, x_grid, dx))

    print_results_table(results)
    save_results_csv(results, outputs_dir / "trinity_record_locking_results.csv")
    save_bar_plot(
        results,
        key="uncertainty_product",
        ylabel="Uncertainty product (delta_x * delta_p)",
        title="Trinity-inspired uncertainty product by profile",
        output_path=outputs_dir / "trinity_uncertainty_products.png",
        reference_line=HBAR / 2.0,
        reference_label="hbar / 2",
    )
    save_bar_plot(
        results,
        key="S_balance",
        ylabel="Concentration-dispersion balance score",
        title="Trinity-inspired concentration-dispersion balance by profile",
        output_path=outputs_dir / "trinity_concentration_dispersion_balance.png",
    )


if __name__ == "__main__":
    main()
