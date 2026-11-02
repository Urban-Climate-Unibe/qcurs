# archiv

Files that are no longer part of the package but are kept for reference.

- `R_0.2.0/`: the humidity levels under their 0.2.0 numbering (4 = saturation
  drift, 5 = stuck values, 6 = decoupling) and `load_packages.R`. Replaced in
  0.3.0 by `RH_QC_4_stuck_values`, `RH_QC_5_decoupling`,
  `RH_QC_6_saturation_drift`. Not loaded, not built (`.Rbuildignore`).
- `man_0.2.0/`: the Rd files that belonged to them.
- `vignettes/`: the former vignette, a copy of the R Markdown template in
  `inst/rmarkdown/templates/qc-workflow`; the template is the one to use.
- `R_0.4.2/`, `man_0.4.2/`, `tests_0.4.2/`: `T_QC_9_bias_shift` (0.4.0-0.4.2),
  its Rd and its tests. Removed in 0.5.0: on Thun-Steffisburg 2026 it objected
  only to site weather that lasted days (fog in the town with sun at Goldiwil,
  the cold lake in summer) and never to a sensor. The clock offset became
  level 9, the interpolation level 10. The archived tests use `shift_fixture()`
  from `tests/testthat/test-clock-offset.R`.
