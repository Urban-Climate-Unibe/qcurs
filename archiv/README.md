# archiv

Files that are no longer part of the package but are kept for reference.

- `R_0.2.0/`: the humidity levels under their 0.2.0 numbering (4 = saturation
  drift, 5 = stuck values, 6 = decoupling) and `load_packages.R`. Replaced in
  0.3.0 by `RH_QC_4_stuck_values`, `RH_QC_5_decoupling`,
  `RH_QC_6_saturation_drift`. Not loaded, not built (`.Rbuildignore`).
- `man_0.2.0/`: the Rd files that belonged to them.
- `vignettes/`: the former vignette, a copy of the R Markdown template in
  `inst/rmarkdown/templates/qc-workflow`; the template is the one to use.
