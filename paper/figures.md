# Figure Inventory

## Phase 1

- `outputs/uncertainty_products.png`
  Shows the uncertainty product `delta_x * delta_p` for the baseline Phase 1 profiles. Gaussian-family states sit at the minimum, while non-Gaussian states rise above it.

## Phase 2

- `outputs/trinity_uncertainty_products.png`
  Shows the uncertainty product for the Trinity-inspired structured profiles.

- `outputs/trinity_concentration_dispersion_balance.png`
  Shows the toy concentration-dispersion balance score `S_balance` for the Phase 2 structured profiles.

## Phase 3

- `outputs/commutator_relative_error.png`
  Shows the relative residual error in the canonical commutator test `[x_hat, p_hat] psi ≈ i hbar psi` across the selected profile set.

## Phase 4

- `outputs/sweep_uncertainty_vs_parameter.png`
  Shows the uncertainty product as a function of the swept parameter for each Phase 4 family.

- `outputs/sweep_s_balance_vs_parameter.png`
  Shows the toy `S_balance` score as a function of the swept parameter for each Phase 4 family.

- `outputs/sweep_uncertainty_vs_s_balance.png`
  Scatter plot of uncertainty product versus `S_balance`, grouped by family. This is one of the clearest views of the separation between minimum uncertainty and maximum structural balance.

- `outputs/sweep_ratio_vs_s_balance.png`
  Scatter plot of `ratio_to_hbar_over_2` versus `S_balance`, grouped by family.

## Phase 4B

- `outputs/robustness_combined_balance_by_family.png`
  Shows the best `S_combined` value attained by each family in the robustness test.

- `outputs/robustness_original_vs_log_balance.png`
  Compares the original balance metric with the logarithmic balance metric across all sweep runs.

- `outputs/robustness_entropy_vs_combined_balance.png`
  Compares entropy-balance scores with the combined structural score.

- `outputs/robustness_uncertainty_vs_combined_balance.png`
  Shows how the combined structural score relates to the uncertainty product, highlighting that structural balance and minimum uncertainty are not the same criterion.
