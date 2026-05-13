"""Numerical uncertainty test for a 1D scalar-clock record-space toy model.

Trinity/STO interpretation for this script:
- R(x) represents the record-localization face.
- R(k) represents the action-gradient / momentum face.
- hbar is treated here as the standard conversion scale p = hbar * k.

This code verifies recovery of the standard Fourier uncertainty floor in a
controlled toy model. It does not claim a new derivation of hbar or new
quantum mechanics.
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
    """Compute a continuous-convention Fourier transform and angular k grid.

    numpy.fft.fftfreq returns frequencies in cycles per unit length. We convert
    them to angular wave numbers through k = 2*pi*f so that p = hbar*k is the
    standard momentum mapping used in the uncertainty relation.
    """
    k_grid = 2.0 * np.pi * np.fft.fftshift(np.fft.fftfreq(record_x.size, d=dx))
    record_k = (
        np.fft.fftshift(np.fft.fft(np.fft.ifftshift(record_x))) * dx / np.sqrt(2.0 * np.pi)
    )
    return k_grid, record_k


def gaussian_profile(x: np.ndarray, sigma: float) -> np.ndarray:
    return np.exp(-(x**2) / (4.0 * sigma**2))


def square_pulse_profile(x: np.ndarray, half_width: float) -> np.ndarray:
    return np.where(np.abs(x) <= half_width, 1.0, 0.0)


def exponential_decay_profile(x: np.ndarray, scale: float) -> np.ndarray:
    return np.exp(-np.abs(x) / scale)


def double_gaussian_profile(
    x: np.ndarray,
    sigma: float,
    separation: float,
    relative_phase: float = 0.0,
) -> np.ndarray:
    left = np.exp(-((x + separation / 2.0) ** 2) / (4.0 * sigma**2))
    right = np.exp(-((x - separation / 2.0) ** 2) / (4.0 * sigma**2)) * np.exp(1j * relative_phase)
    return left + right


def build_profiles(x: np.ndarray) -> list[tuple[str, np.ndarray]]:
    base_sigma = 4.0
    phase_k0 = 2.5
    return [
        ("Gaussian", gaussian_profile(x, sigma=base_sigma)),
        ("narrow Gaussian", gaussian_profile(x, sigma=1.0)),
        ("broad Gaussian", gaussian_profile(x, sigma=8.0)),
        ("square/windowed pulse", square_pulse_profile(x, half_width=3.0)),
        ("exponential decay", exponential_decay_profile(x, scale=3.0)),
        ("double Gaussian", double_gaussian_profile(x, sigma=1.2, separation=8.0)),
        ("phase-twisted Gaussian", gaussian_profile(x, sigma=base_sigma) * np.exp(1j * phase_k0 * x)),
    ]


def evaluate_profile(name: str, amplitude_x: np.ndarray, x_grid: np.ndarray, dx: float) -> dict[str, float | str]:
    record_x = normalize_amplitude(amplitude_x.astype(np.complex128), dx)
    pdf_x = np.abs(record_x) ** 2

    _, var_x = expectation_from_pdf(x_grid, pdf_x, dx)
    delta_x = np.sqrt(var_x)

    k_grid, record_k = fourier_transform_record(record_x, dx)
    dk = float(k_grid[1] - k_grid[0])
    pdf_k = np.abs(record_k) ** 2

    # Renormalize in k-space to keep the numerical integration exactly unitary
    # even after finite-grid truncation.
    pdf_k /= np.sum(pdf_k) * dk

    _, var_k = expectation_from_pdf(k_grid, pdf_k, dk)
    delta_k = np.sqrt(var_k)
    delta_p = HBAR * delta_k
    uncertainty_product = delta_x * delta_p
    ratio_to_hbar_over_2 = uncertainty_product / (HBAR / 2.0)
    passed = uncertainty_product >= (HBAR / 2.0 - TOLERANCE)

    return {
        "profile_name": name,
        "delta_x": delta_x,
        "delta_k": delta_k,
        "delta_p": delta_p,
        "uncertainty_product": uncertainty_product,
        "ratio_to_hbar_over_2": ratio_to_hbar_over_2,
        "pass_fail": "PASS" if passed else "FAIL",
    }


def print_results_table(results: list[dict[str, float | str]]) -> None:
    headers = [
        "profile_name",
        "delta_x",
        "delta_k",
        "delta_p",
        "uncertainty_product",
        "ratio_to_hbar_over_2",
        "pass_fail",
    ]
    rows = []
    for result in results:
        rows.append(
            [
                str(result["profile_name"]),
                f'{result["delta_x"]:.6f}',
                f'{result["delta_k"]:.6f}',
                f'{result["delta_p"]:.6f}',
                f'{result["uncertainty_product"]:.6f}',
                f'{result["ratio_to_hbar_over_2"]:.6f}',
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
        "delta_x",
        "delta_k",
        "delta_p",
        "uncertainty_product",
        "ratio_to_hbar_over_2",
        "pass_fail",
    ]
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            writer.writerow(result)


def save_uncertainty_plot(results: list[dict[str, float | str]], output_path: Path) -> None:
    profile_names = [str(result["profile_name"]) for result in results]
    uncertainty_products = [float(result["uncertainty_product"]) for result in results]

    fig, ax = plt.subplots(figsize=(11, 5.5))
    bars = ax.bar(profile_names, uncertainty_products, color="#2c7fb8")
    ax.axhline(HBAR / 2.0, color="#d95f0e", linestyle="--", linewidth=2, label="hbar / 2")
    ax.set_ylabel("Uncertainty product (delta_x * delta_p)")
    ax.set_title("Uncertainty product by record profile")
    ax.tick_params(axis="x", rotation=25)
    ax.legend()

    for bar, value in zip(bars, uncertainty_products):
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
    save_results_csv(results, outputs_dir / "uncertainty_results.csv")
    save_uncertainty_plot(results, outputs_dir / "uncertainty_products.png")


if __name__ == "__main__":
    main()
