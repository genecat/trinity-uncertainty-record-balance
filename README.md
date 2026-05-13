# trinity-uncertainty-record-balance

Numerical toy-model tests for the repository and draft paper:

**A Record-Balance Wall in Scalar-Clock Toy Models: Concentration–Dispersion Structure Beyond Minimum Uncertainty**

Paper/preprint DOI: `10.5281/zenodo.20172951`

Software archive DOI: `10.5281/zenodo.20172621`

## Short Description

This repository studies a 1D scalar-clock record-space toy model in which:

- `R(x)` is interpreted as a record-localization face.
- `R(k)` is interpreted as an action-gradient / momentum face.
- `p = hbar * k` uses the standard quantum conversion scale with `hbar = 1.0` in the numerical experiments.

The codebase is organized as a reproducible research package with five numerical phases:

1. Phase 1: recovery of the standard Heisenberg uncertainty floor.
2. Phase 2: Trinity/STO-inspired record-locking profiles.
3. Phase 3: canonical position-momentum commutator recovery.
4. Phase 4: concentration-dispersion balance sweeps.
5. Phase 4B: robustness tests across alternative balance metrics.

## Relation To The Paper

This repository is the computational companion to the draft manuscript:

**A Record-Balance Wall in Scalar-Clock Toy Models: Concentration–Dispersion Structure Beyond Minimum Uncertainty**

The current code is intended to support a public research release by:

- preserving the scripts used to generate the reported outputs,
- keeping CSV and PNG reproducibility artifacts under version control,
- documenting the claim boundaries explicitly,
- and providing a manuscript-ready summary of methods and results.

## Claim Boundary

This repository makes a narrow compatibility and exploratory-structure claim:

- It does **not** derive `hbar` from first principles.
- It does **not** replace standard quantum mechanics.
- It does **not** prove Trinity/STO.
- It **does** test whether a Trinity/STO-style interpretive mapping can reproduce standard uncertainty and commutator structure in a controlled toy model.
- It **does** identify a reproducible toy-model separation between minimum uncertainty and maximum concentration-dispersion balance in the present numerical setup.

## Installation

### Recommended Python Version

Python `3.11` or newer is recommended. The scripts use only:

- `numpy`
- `matplotlib`

### Setup

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## How To Run Each Phase

Run commands from the repository root.

### Phase 1: Uncertainty Recovery

```bash
python3 src/uncertainty_record_test.py
```

Outputs:

- `outputs/uncertainty_results.csv`
- `outputs/uncertainty_products.png`

### Phase 2: Trinity-Inspired Record Locking

```bash
python3 src/trinity_record_locking_test.py
```

Outputs:

- `outputs/trinity_record_locking_results.csv`
- `outputs/trinity_uncertainty_products.png`
- `outputs/trinity_concentration_dispersion_balance.png`

### Phase 3: Operator Commutator Test

```bash
python3 src/record_operator_commutator_test.py
```

Outputs:

- `outputs/record_operator_commutator_results.csv`
- `outputs/commutator_relative_error.png`

### Phase 4: Concentration-Dispersion Sweep

```bash
python3 src/concentration_dispersion_sweep.py
```

Outputs:

- `outputs/concentration_dispersion_sweep_results.csv`
- `outputs/sweep_uncertainty_vs_parameter.png`
- `outputs/sweep_s_balance_vs_parameter.png`
- `outputs/sweep_uncertainty_vs_s_balance.png`
- `outputs/sweep_ratio_vs_s_balance.png`

### Phase 4B: Balance Robustness Test

```bash
python3 src/concentration_dispersion_robustness_test.py
```

Outputs:

- `outputs/concentration_dispersion_robustness_results.csv`
- `outputs/robustness_combined_balance_by_family.png`
- `outputs/robustness_original_vs_log_balance.png`
- `outputs/robustness_entropy_vs_combined_balance.png`
- `outputs/robustness_uncertainty_vs_combined_balance.png`

## Expected Outputs

The high-level expected behavior is:

- Gaussian records saturate the standard uncertainty floor in Phase 1.
- Trinity-inspired structured profiles remain above the floor in Phase 2 while producing distinct toy balance scores.
- The canonical commutator is numerically recovered in Phase 3 with very small residual error for smooth profiles.
- Structured overlap profiles dominate the highest Phase 4 concentration-dispersion balance scores even though they are not minimum-uncertainty states.
- The same overlap-dominated high-balance region remains competitive in Phase 4B under multiple structural balance metrics.

## Reproducibility Notes

- All scripts use the same numerical box: `x in [-50, 50]`, `N = 16384`, `hbar = 1.0`.
- Fourier-space calculations use `numpy.fft.fftfreq`.
- CSV and PNG artifacts in `outputs/` are part of the published reproducibility package and are intentionally **not** ignored by `.gitignore`.
- A known finite-grid caveat appears for the very widest Gaussian states in the fixed `[-50, 50]` box. Those runs can drift slightly below the ideal uncertainty floor due to truncation and finite-domain effects rather than a new structural result.
- For a fuller reproducibility checklist, see [docs/reproducibility.md](/Users/eugenecatrambone/trinity-uncertainty-tests/docs/reproducibility.md).

## Repository Structure

```text
src/     numerical scripts for Phases 1-4B
outputs/ published reproducibility artifacts (CSVs and plots)
docs/    methods, claim boundaries, and reproducibility notes
paper/   manuscript draft, tables, and figure notes
```

## Citation

Citation metadata is provided in [CITATION.cff](/Users/eugenecatrambone/trinity-uncertainty-tests/CITATION.cff).

## Suggested Citation

Paper / preprint:

```text
Catrambone, Eugene. A Record-Balance Wall in Scalar-Clock Toy Models:
Concentration–Dispersion Structure Beyond Minimum Uncertainty. 2026.
Zenodo. https://doi.org/10.5281/zenodo.20172951
```

Software / repository:

```text
Catrambone, Eugene. A Record-Balance Wall in Scalar-Clock Toy Models:
Concentration–Dispersion Structure Beyond Minimum Uncertainty.
trinity-uncertainty-record-balance, 2026.
Zenodo. https://doi.org/10.5281/zenodo.20172621
```

Paper and software archive DOIs:

```text
Paper/preprint DOI: 10.5281/zenodo.20172951
Software DOI: 10.5281/zenodo.20172621
```

Please cite both the repository and the associated preprint when available.

## License

The code in this repository is released under the MIT License. See [LICENSE](/Users/eugenecatrambone/trinity-uncertainty-tests/LICENSE).

Research interpretation, manuscript language, and future publication terms may be updated separately as the preprint and repository mature.
