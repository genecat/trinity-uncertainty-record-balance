# Methods Summary

## Numerical Setup

All phases use the same core 1D toy-model discretization:

- spatial grid `x in [-50, 50]`
- `N = 16384` grid points
- `hbar = 1.0`
- normalized complex record amplitudes
- Fourier-space wave-number grid from `numpy.fft.fftfreq`

The shared interpretive mapping is:

- `R(x)` as a record-localization face
- `R(k)` as an action-gradient / momentum face
- `p = hbar * k` as the standard conversion rule

## Phase 1: Uncertainty Recovery

Phase 1 tests whether the toy record-space representation reproduces the standard Fourier uncertainty floor. The script evaluates several normalized profiles:

- Gaussian
- narrow Gaussian
- broad Gaussian
- square/windowed pulse
- exponential decay
- double Gaussian
- phase-twisted Gaussian

For each profile it computes:

- `delta_x`
- `delta_k`
- `delta_p = hbar * delta_k`
- `uncertainty_product = delta_x * delta_p`
- `ratio_to_hbar_over_2`

The numerical target is recovery of the standard lower bound `delta_x * delta_p >= hbar / 2`.

## Phase 2: Trinity-Inspired Record Profiles

Phase 2 adds structured profiles motivated by Trinity/STO language while keeping the numerical target modest: preserve the standard uncertainty floor while probing distinct spatial and Fourier-space organization.

Profiles tested:

- `boundary_locked_gaussian`
- `sharp_boundary_locked_gaussian`
- `scalar_clock_expansion_record`
- `selector_phase_twisted_record`
- `two_time_overlap_record`
- `asymmetric_overlap_record`
- `oscillatory_selector_record`

Phase 2 also introduces the toy concentration-dispersion diagnostic:

- `C = 1 / (delta_x + eps)`
- `D = delta_k`
- `S_balance = 1 - abs(C - D) / (C + D + eps)`

This score is not a standard quantum observable. It is only a toy structural diagnostic.

## Phase 3: Commutator Test

Phase 3 checks the canonical operator algebra on the same record-space grid.

Operators:

- `x_hat(psi) = x * psi`
- `p_hat(psi) = -i * hbar * d/dx`

The derivative is computed spectrally via FFT, and the test compares:

- `commutator_psi = x_hat(p_hat(psi)) - p_hat(x_hat(psi))`
- `expected_psi = i * hbar * psi`

Diagnostics include relative residual norm, normalized overlap with the expected result, and maximum pointwise residual. The goal is compatibility with the standard non-commuting operator structure.

## Phase 4: Concentration-Dispersion Sweep

Phase 4 scans six profile families using 41 samples each:

1. Gaussian width sweep
2. Boundary softness sweep
3. Overlap separation sweep
4. Asymmetric overlap sweep
5. Selector phase sweep
6. Oscillatory selector sweep

For each run it records uncertainty and toy structural diagnostics, then summarizes:

- family-level minimum and maximum uncertainty products
- maximum `S_balance`
- parameter values at the best-balance points
- the top-ranked records by `S_balance`

This phase isolates a reproducible difference between minimum uncertainty and maximum concentration-dispersion balance in the toy model.

## Phase 4B: Robustness Test

Phase 4B asks whether the high-balance regions found in Phase 4 are stable under alternative balance definitions. It reuses the same six families and 41-sample grids, then evaluates:

- `S_balance_original`
- `S_balance_log`
- `S_balance_product_distance`
- `S_balance_entropy`
- `S_combined`

`S_balance_product_distance` is treated as a control because it measures closeness to the Heisenberg floor rather than structural balance. The core robustness question is whether the same parameter windows stay prominent under:

- original concentration-dispersion matching
- logarithmic concentration-dispersion matching
- entropy-balance matching
- the combined structural score
