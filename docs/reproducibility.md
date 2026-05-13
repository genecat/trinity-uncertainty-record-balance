# Reproducibility Notes

Paper/preprint DOI: `10.5281/zenodo.20172951`

Archived software DOI: `10.5281/zenodo.20172621`

The Zenodo-tagged release `v0.1.1` is the archived software snapshot corresponding to DOI `10.5281/zenodo.20172621`. It includes archival metadata updates only and makes no methodological changes relative to `v0.1.0`.

## Recommended Python Version

Python `3.11` or newer is recommended.

## Dependencies

Install the project dependencies with:

```bash
pip install -r requirements.txt
```

The repository uses only:

- `numpy`
- `matplotlib`

## Commands To Reproduce All Outputs

Run these commands from the repository root:

```bash
python3 src/uncertainty_record_test.py
python3 src/trinity_record_locking_test.py
python3 src/record_operator_commutator_test.py
python3 src/concentration_dispersion_sweep.py
python3 src/concentration_dispersion_robustness_test.py
```

## Output Locations

All published reproducibility artifacts are written to `outputs/`.

CSV outputs:

- `outputs/uncertainty_results.csv`
- `outputs/trinity_record_locking_results.csv`
- `outputs/record_operator_commutator_results.csv`
- `outputs/concentration_dispersion_sweep_results.csv`
- `outputs/concentration_dispersion_robustness_results.csv`

Figure outputs:

- `outputs/uncertainty_products.png`
- `outputs/trinity_uncertainty_products.png`
- `outputs/trinity_concentration_dispersion_balance.png`
- `outputs/commutator_relative_error.png`
- `outputs/sweep_uncertainty_vs_parameter.png`
- `outputs/sweep_s_balance_vs_parameter.png`
- `outputs/sweep_uncertainty_vs_s_balance.png`
- `outputs/sweep_ratio_vs_s_balance.png`
- `outputs/robustness_combined_balance_by_family.png`
- `outputs/robustness_original_vs_log_balance.png`
- `outputs/robustness_entropy_vs_combined_balance.png`
- `outputs/robustness_uncertainty_vs_combined_balance.png`

## Finite-Grid Caveat

The widest Gaussian profiles in the fixed `x in [-50, 50]` box can drift slightly below the ideal uncertainty floor by a very small amount. In the current runs, that appears for the broadest entries in the Phase 4 and Phase 4B Gaussian-width sweeps.

This is best interpreted as a finite-domain truncation artifact:

- the Gaussian becomes too wide for the fixed numerical box,
- the tails are cut off,
- and the discrete approximation to the ideal continuum Gaussian is no longer exact.

This should not be interpreted as evidence against the standard uncertainty relation or as evidence for a new physical effect.

## Reproducibility Philosophy

The repository intentionally keeps CSV and PNG outputs under version control so that:

- readers can inspect the reported results without rerunning the code first,
- manuscript tables can be traced back to exact numerical artifacts,
- and later code revisions can be compared against the current baseline.
