"""Phase 3 canonical commutator test on the 1D record-space toy model.

Interpretive mapping used here:
- x_hat is the record-localization projection.
- p_hat = -i * hbar * d/dx is the temporal-action / action-gradient projection.

The commutator [x_hat, p_hat] tests whether projection order leaves the
standard phase-rotation residue i * hbar on the record amplitudes. This is a
standard numerical operator-algebra check on a toy record space, not a
first-principles derivation of Trinity/STO or of hbar.
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
PASS_THRESHOLD = 1e-2
WARN_THRESHOLD = 5e-2


def normalize_amplitude(amplitude: np.ndarray, dx: float) -> np.ndarray:
    """Normalize a complex amplitude so integral |psi(x)|^2 dx = 1."""
    norm = np.sqrt(np.sum(np.abs(amplitude) ** 2) * dx)
    if norm == 0.0:
        raise ValueError("Amplitude cannot be identically zero.")
    return amplitude / norm


def l2_norm(field: np.ndarray, dx: float) -> float:
    """Continuous L2 norm on the x grid."""
    return float(np.sqrt(np.sum(np.abs(field) ** 2) * dx))


def inner_product(left: np.ndarray, right: np.ndarray, dx: float) -> complex:
    """Continuous inner product <left|right> on the x grid."""
    return np.vdot(left, right) * dx


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


def overlap_record(x: np.ndarray, sigma: float, separation: float) -> np.ndarray:
    left = gaussian_profile(x, sigma=sigma, center=-separation / 2.0)
    right = gaussian_profile(x, sigma=sigma, center=separation / 2.0)
    return left + right


def oscillatory_selector_record(x: np.ndarray, sigma: float, beta: float, omega: float) -> np.ndarray:
    envelope = gaussian_profile(x, sigma=sigma)
    return envelope * np.exp(1j * beta * np.sin(omega * x))


def build_profiles(x: np.ndarray) -> list[tuple[str, np.ndarray]]:
    boundary_locked = boundary_locked_gaussian_profile(x, sigma=4.0, boundary_scale=10.0, softness=1.5)
    sharp_boundary_locked = boundary_locked_gaussian_profile(x, sigma=4.0, boundary_scale=10.0, softness=0.25)
    return [
        ("Gaussian", gaussian_profile(x, sigma=4.0)),
        ("boundary_locked_gaussian", boundary_locked),
        ("sharp_boundary_locked_gaussian", sharp_boundary_locked),
        ("selector_phase_twisted_record", boundary_locked * np.exp(1j * 1.8 * x)),
        ("two_time_overlap_record", overlap_record(x, sigma=1.5, separation=10.0)),
        ("oscillatory_selector_record", oscillatory_selector_record(x, sigma=4.0, beta=0.9, omega=0.7)),
    ]


def spectral_derivative(psi: np.ndarray, dx: float) -> np.ndarray:
    """Compute dpsi/dx using the FFT spectral derivative on the periodic grid."""
    k_grid = 2.0 * np.pi * np.fft.fftfreq(psi.size, d=dx)
    psi_hat = np.fft.fft(psi)
    return np.fft.ifft(1j * k_grid * psi_hat)


def x_hat(psi: np.ndarray, x_grid: np.ndarray) -> np.ndarray:
    return x_grid * psi


def p_hat(psi: np.ndarray, dx: float) -> np.ndarray:
    dpsi_dx = spectral_derivative(psi, dx)
    return -1j * HBAR * dpsi_dx


def classify_result(profile_name: str, relative_error: float) -> str:
    if relative_error < PASS_THRESHOLD:
        return "PASS"
    if profile_name == "sharp_boundary_locked_gaussian" and relative_error < WARN_THRESHOLD:
        return "WARN"
    return "FAIL"


def evaluate_profile(name: str, amplitude_x: np.ndarray, x_grid: np.ndarray, dx: float) -> dict[str, float | str]:
    psi = normalize_amplitude(amplitude_x.astype(np.complex128), dx)
    commutator_psi = x_hat(p_hat(psi, dx), x_grid) - p_hat(x_hat(psi, x_grid), dx)
    expected_psi = 1j * HBAR * psi
    residual = commutator_psi - expected_psi

    norm_expected = l2_norm(expected_psi, dx)
    norm_residual = l2_norm(residual, dx)
    relative_error = norm_residual / norm_expected

    normalized_overlap = inner_product(expected_psi, commutator_psi, dx) / (norm_expected**2)
    max_abs_residual = float(np.max(np.abs(residual)))

    return {
        "profile_name": name,
        "relative_error": relative_error,
        "overlap_real": float(np.real(normalized_overlap)),
        "overlap_imag": float(np.imag(normalized_overlap)),
        "max_abs_residual": max_abs_residual,
        "pass_fail": classify_result(name, relative_error),
    }


def print_results_table(results: list[dict[str, float | str]]) -> None:
    headers = [
        "profile_name",
        "relative_error",
        "overlap_real",
        "overlap_imag",
        "max_abs_residual",
        "pass_fail",
    ]
    rows = []
    for result in results:
        rows.append(
            [
                str(result["profile_name"]),
                f'{result["relative_error"]:.6e}',
                f'{result["overlap_real"]:.6f}',
                f'{result["overlap_imag"]:.6f}',
                f'{result["max_abs_residual"]:.6e}',
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
        "relative_error",
        "overlap_real",
        "overlap_imag",
        "max_abs_residual",
        "pass_fail",
    ]
    with output_path.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        for result in results:
            writer.writerow(result)


def save_relative_error_plot(results: list[dict[str, float | str]], output_path: Path) -> None:
    profile_names = [str(result["profile_name"]) for result in results]
    values = [float(result["relative_error"]) for result in results]

    fig, ax = plt.subplots(figsize=(11.5, 5.8))
    bars = ax.bar(profile_names, values, color="#2c7fb8")
    ax.axhline(PASS_THRESHOLD, color="#d95f0e", linestyle="--", linewidth=2, label="PASS threshold")
    ax.axhline(WARN_THRESHOLD, color="#756bb1", linestyle=":", linewidth=2, label="WARN threshold")
    ax.set_yscale("log")
    ax.set_ylabel("Relative commutator error")
    ax.set_title("Canonical commutator relative error by record profile")
    ax.tick_params(axis="x", rotation=25)
    ax.legend()

    for bar, value in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2.0,
            value,
            f"{value:.2e}",
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
    save_results_csv(results, outputs_dir / "record_operator_commutator_results.csv")
    save_relative_error_plot(results, outputs_dir / "commutator_relative_error.png")


if __name__ == "__main__":
    main()
