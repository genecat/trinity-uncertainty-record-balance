# A Record-Balance Wall in Scalar-Clock Toy Models:
# Concentration–Dispersion Structure Beyond Minimum Uncertainty

## Abstract

We study a one-dimensional scalar-clock record-space toy model designed to test a narrow question: whether a Trinity/STO-style interpretive mapping can remain compatible with standard quantum uncertainty and operator structure while also supporting an additional toy diagnostic of concentration-dispersion balance. The model treats `R(x)` as a record-localization face and `R(k)` as an action-gradient / momentum face, with `p = hbar * k` and `hbar = 1` in the numerical experiments. Phase 1 verifies recovery of the standard Fourier uncertainty floor for Gaussian and non-Gaussian records. Phase 2 introduces structured record-locking profiles, including boundary-locked, overlap, and oscillatory selector states, and shows that they remain above the standard uncertainty floor while generating distinct concentration-dispersion scores. Phase 3 verifies numerical recovery of the canonical commutator `[x_hat, p_hat] psi ≈ i hbar psi`. Phase 4 then sweeps parameterized record families and finds a reproducible separation between minimum-uncertainty states and maximum concentration-dispersion balance states. Gaussian states saturate the Heisenberg floor, but overlap-dominated records produce the highest toy structural balance scores. Phase 4B tests robustness under alternative balance metrics and shows that the overlap-separation family remains dominant under the combined structural score. We emphasize that this is exploratory numerical work on a toy model: it does not derive `hbar` from first principles, replace standard quantum mechanics, or prove Trinity/STO.

## Introduction

Many speculative interpretive programs begin by asking whether familiar quantum structure can be reproduced in a different conceptual language before stronger claims are attempted. The present work takes that limited approach. Rather than proposing a new mechanics, we test whether a simple scalar-clock record-space toy model can reproduce two standard benchmarks:

1. the Heisenberg uncertainty floor, and
2. the canonical position-momentum commutator.

Once those benchmarks are recovered numerically, we ask a second question that is explicitly exploratory rather than foundational: can a separate toy diagnostic reveal a structural regime that is not identical to minimum uncertainty?

This question is motivated by the idea that a record-space model might admit two distinguishable notions of "specialness":

- minimum uncertainty, associated with Gaussian saturation of the Fourier bound, and
- concentration-dispersion balance, associated with a more even structural matching between spatial concentration and momentum-space dispersion.

The present repository investigates that distinction numerically.

## Definitions

### Record-Space Variables

Let `R(x)` denote a normalized complex record amplitude on a finite 1D grid. The corresponding Fourier amplitude `R(k)` is computed numerically using FFT conventions consistent with `numpy.fft.fftfreq`, with angular wave number `k = 2 pi f`.

The numerical experiments fix:

- `x in [-50, 50]`
- `N = 16384`
- `hbar = 1`

### Standard Quantities

For each record we compute:

- `delta_x`
- `delta_k`
- `delta_p = hbar * delta_k`
- `uncertainty_product = delta_x * delta_p`
- `ratio_to_hbar_over_2 = uncertainty_product / (hbar / 2)`

### Toy Structural Diagnostic

The original concentration-dispersion balance score is:

- `C = 1 / (delta_x + eps)`
- `D = delta_k`
- `S_balance = 1 - abs(C - D) / (C + D + eps)`

This is not a standard quantum observable. It is a toy structural diagnostic only.

## Methods

### Phase 1: Uncertainty Recovery

We evaluate Gaussian, narrow Gaussian, broad Gaussian, square/windowed pulse, exponential decay, double Gaussian, and phase-twisted Gaussian profiles. The goal is to recover the standard uncertainty lower bound and verify that Gaussian-family states sit at or extremely near the minimum.

### Phase 2: Trinity-Inspired Record Profiles

We introduce structured profiles intended to mimic Trinity/STO-style record language without making any first-principles claim. These include:

- boundary-locked Gaussians,
- sharp boundary-locked Gaussians,
- scalar-clock expansion records,
- selector phase-twisted records,
- two-time overlap records,
- asymmetric overlap records,
- oscillatory selector records.

These profiles are used to test whether structure-rich records can remain compatible with the standard uncertainty floor while producing distinctive values of the toy balance score.

### Phase 3: Canonical Commutator Test

We define:

- `x_hat(psi) = x * psi`
- `p_hat(psi) = -i hbar d/dx`

with the derivative computed spectrally via FFT. For each test profile we compare:

- `[x_hat, p_hat] psi`
- `i hbar psi`

using residual norms, overlaps, and maximum pointwise residuals.

### Phase 4: Parameter Sweeps

We sweep six structured families using 41 samples each:

1. Gaussian width
2. Boundary softness
3. Overlap separation
4. Asymmetric overlap
5. Selector phase
6. Oscillatory selector strength

The purpose is to map regions of:

- minimum uncertainty,
- high concentration-dispersion balance,
- and possible toy-model balance walls or stable structural plateaus.

### Phase 4B: Robustness Metrics

To test whether Phase 4 is formula-specific, we compute:

- `S_balance_original`
- `S_balance_log`
- `S_balance_product_distance`
- `S_balance_entropy`
- `S_combined`

`S_balance_product_distance` is treated separately because it is constructed to reward minimum uncertainty rather than structural balance.

## Results

### Uncertainty Recovery

Phase 1 reproduces the standard uncertainty floor as expected. Gaussian, narrow Gaussian, broad Gaussian, and phase-twisted Gaussian states all land at or numerically indistinguishable from `delta_x delta_p = hbar / 2`. Non-Gaussian states such as square pulses, exponential decay records, and double Gaussians sit above the floor.

### Structured Record Profiles

Phase 2 shows that Trinity-inspired structured records preserve the uncertainty floor while displaying distinct toy balance values. In particular, overlap and asymmetry can raise the balance score even when the uncertainty product is no longer minimal.

### Commutator Compatibility

Phase 3 numerically recovers `[x_hat, p_hat] psi ≈ i hbar psi` with extremely small residual error for the smooth test profiles considered. This supports the limited claim that the toy record-space representation is compatible with the standard non-commuting operator structure.

### Balance-Wall Search

Phase 4 produces the main exploratory result of the repository. The highest `S_balance` values are not Gaussian minimum-uncertainty states. Instead, the strongest balance scores occur in overlap-dominated record families, especially the overlap-separation sweep near `d ≈ 9.5`. This is the key separation between:

- uncertainty minimization, which remains Gaussian-dominated, and
- concentration-dispersion balance, which appears structure-dominated in this toy model.

### Robustness

Phase 4B shows that the overlap-separation family remains the best performer under the combined structural score `S_combined`, which averages the original balance score, a logarithmic balance score, and an entropy-balance score. This suggests that the overlap-driven peak is not merely an artifact of the original `S_balance` definition.

## Discussion

The present results suggest a useful conceptual distinction. A record family can be optimal with respect to the Heisenberg floor without being optimal with respect to a toy structural balance criterion. Conversely, a highly structured record can score strongly on concentration-dispersion balance even while sitting well above the minimum uncertainty floor.

Within the limited language of this toy model, that separation resembles a "record-balance wall": a region where structured overlap records appear to achieve especially strong concentration-dispersion matching without collapsing back onto Gaussian minimum-uncertainty behavior.

This does not imply a new physical law. It does, however, suggest that if the Trinity/STO program is pursued further, balance diagnostics should be studied separately from uncertainty minimization rather than folded into it.

## Limitations

- The model is one-dimensional and highly simplified.
- The numerical domain is finite, with periodic FFT conventions.
- The balance metrics are toy diagnostics rather than standard observables.
- The broadest Gaussian states in the fixed `[-50, 50]` domain show small finite-box artifacts in the sweep phases.
- No first-principles derivation of quantum mechanics is attempted.

## Conclusion

This repository establishes three limited results. First, the scalar-clock record-space toy model reproduces the standard uncertainty floor for the tested profiles. Second, it reproduces the canonical position-momentum commutator numerically. Third, once those standard benchmarks are in place, it reveals a reproducible toy-model separation between minimum uncertainty and maximum concentration-dispersion balance. In the present experiments, Gaussian states remain the minimum-uncertainty states, while overlap-dominated structured records dominate the highest structural balance scores. That separation survives the Phase 4B robustness tests and therefore warrants further toy-model investigation, while remaining well inside a clearly bounded exploratory claim.

## References

- Placeholder for associated preprint.
- Placeholder for standard Fourier uncertainty references.
- Placeholder for numerical spectral-method references.
