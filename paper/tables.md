# Key Tables

## Phase 1: Uncertainty Recovery

| Profile | delta_x | delta_k | delta_p | uncertainty_product | ratio_to_hbar_over_2 | Status |
|---|---:|---:|---:|---:|---:|---|
| Gaussian | 4.000000 | 0.125000 | 0.125000 | 0.500000 | 1.000000 | PASS |
| narrow Gaussian | 1.000000 | 0.500000 | 0.500000 | 0.500000 | 1.000000 | PASS |
| broad Gaussian | 8.000000 | 0.062500 | 0.062500 | 0.500000 | 1.000000 | PASS |
| square/windowed pulse | 1.731979 | 8.701337 | 8.701337 | 15.070537 | 30.141073 | PASS |
| exponential decay | 2.121319 | 0.333533 | 0.333533 | 0.707530 | 1.415060 | PASS |
| double Gaussian | 4.168739 | 0.407655 | 0.407655 | 1.699406 | 3.398813 | PASS |
| phase-twisted Gaussian | 4.000000 | 0.125000 | 0.125000 | 0.500000 | 1.000000 | PASS |

## Phase 2: Trinity-Inspired Structured Profiles

| Profile | mean_x | mean_k | delta_x | delta_k | uncertainty_product | S_balance | Status |
|---|---:|---:|---:|---:|---:|---:|---|
| boundary_locked_gaussian | 0.000000 | 0.000000 | 3.562601 | 0.142851 | 0.508921 | 0.674549 | PASS |
| sharp_boundary_locked_gaussian | 0.000000 | 0.000000 | 3.789295 | 0.152446 | 0.577661 | 0.732301 | PASS |
| scalar_clock_expansion_record | -0.000000 | -0.000000 | 3.534695 | 0.152074 | 0.537536 | 0.699217 | PASS |
| selector_phase_twisted_record | -0.000000 | 1.800000 | 3.562601 | 0.142851 | 0.508921 | 0.674549 | PASS |
| two_time_overlap_record | -0.000000 | 0.000000 | 5.210924 | 0.326124 | 1.699406 | 0.740904 | PASS |
| asymmetric_overlap_record | -3.306441 | 0.000000 | 4.030548 | 0.327947 | 1.321807 | 0.861398 | PASS |
| oscillatory_selector_record | 0.000000 | 0.012500 | 4.000000 | 0.462514 | 1.850054 | 0.701741 | PASS |

## Phase 3: Canonical Commutator Test

| Profile | relative_error | overlap_real | overlap_imag | max_abs_residual | Status |
|---|---:|---:|---:|---:|---|
| Gaussian | 1.533237e-12 | 1.000000 | -0.000000 | 7.164712e-13 | PASS |
| boundary_locked_gaussian | 9.944468e-06 | 1.000000 | -0.000000 | 1.326703e-06 | PASS |
| sharp_boundary_locked_gaussian | 1.654165e-12 | 1.000000 | -0.000000 | 9.935992e-13 | PASS |
| selector_phase_twisted_record | 9.944779e-06 | 1.000000 | -0.000000 | 1.326744e-06 | PASS |
| two_time_overlap_record | 1.760499e-12 | 1.000000 | -0.000000 | 8.157717e-13 | PASS |
| oscillatory_selector_record | 1.884709e-12 | 1.000000 | -0.000000 | 8.345942e-13 | PASS |

## Phase 4: Family Summary

| Family | Runs | min_uncertainty_product | max_uncertainty_product | max_S_balance | parameter_at_max_S_balance | min_ratio_to_hbar_over_2 | Failures |
|---|---:|---:|---:|---:|---:|---:|---:|
| asymmetric_overlap_sweep | 41 | 0.529126 | 0.769168 | 0.869525 | 1.000000 | 1.058252 | 0 |
| boundary_softness_sweep | 41 | 0.500876 | 1.182727 | 0.916285 | 0.050000 | 1.001753 | 0 |
| gaussian_width_sweep | 41 | 0.499718 | 0.500000 | 0.666667 | 6.250000 | 0.999435 | 14 |
| oscillatory_selector_sweep | 41 | 0.500000 | 21.219095 | 0.920554 | 0.250000 | 1.000000 | 0 |
| overlap_separation_sweep | 41 | 0.500000 | 2.549386 | 0.980902 | 9.500000 | 1.000000 | 0 |
| selector_phase_sweep | 41 | 0.500000 | 0.500000 | 0.666667 | 3.125000 | 1.000000 | 0 |

## Phase 4: Top Balance Records

| Family | Parameter | Value | S_balance | uncertainty_product | ratio_to_hbar_over_2 |
|---|---|---:|---:|---:|---:|
| overlap_separation_sweep | separation_d | 9.500000 | 0.980902 | 1.038939 | 2.077878 |
| overlap_separation_sweep | separation_d | 9.000000 | 0.971349 | 0.944294 | 1.888588 |
| overlap_separation_sweep | separation_d | 10.000000 | 0.936954 | 1.134575 | 2.269151 |
| overlap_separation_sweep | separation_d | 8.500000 | 0.920870 | 0.853345 | 1.706689 |
| oscillatory_selector_sweep | beta | 0.250000 | 0.920554 | 1.172604 | 2.345208 |

## Phase 4B: Robustness Summary

| Family | max_S_original | param | max_S_log | param | max_S_entropy | param | max_S_combined | param | min_uncertainty_product | param |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| asymmetric_overlap_sweep | 0.869525 | 1.000000 | 0.769168 | 1.000000 | 0.530215 | 0.100000 | 0.689209 | 1.000000 | 0.529126 | 0.100000 |
| boundary_softness_sweep | 0.916285 | 0.050000 | 0.845504 | 0.050000 | 0.452790 | 0.050000 | 0.738193 | 0.050000 | 0.500876 | 2.896250 |
| gaussian_width_sweep | 0.666667 | 6.250000 | 0.500000 | 6.250000 | 0.831158 | 0.500000 | 0.665942 | 0.500000 | 0.499718 | 12.000000 |
| oscillatory_selector_sweep | 0.920554 | 0.250000 | 0.852803 | 0.250000 | 0.640800 | 4.625000 | 0.730811 | 0.250000 | 0.500000 | 0.000000 |
| overlap_separation_sweep | 0.980902 | 9.500000 | 0.962520 | 9.500000 | 0.545444 | 0.000000 | 0.790399 | 9.500000 | 0.500000 | 0.000000 |
| selector_phase_sweep | 0.666667 | 3.125000 | 0.500000 | 3.250000 | 0.402587 | 0.250000 | 0.523085 | 0.875000 | 0.500000 | 4.375000 |
