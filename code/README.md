# code/

This directory was reorganized (2026-07) out of a flat, 37-file layout into
four subdirectories, split by role:

- **`R/`** - Production function libraries, sourced by `_targets.R` via
  `tar_source("code/R")`. Everything the pipeline calls lives here:
  `depth_simulation.R`, `property_simulation.R` (both split out of the
  original `code/coweeta/sim-functions.R`), `soil_infill_functions.R`,
  `ssurgo` download logic (inside `soil_infill_functions.R`),
  `distribution_fitting.R`, and `bayesian_updating.R`.
- **`notebooks/`** - Exploratory / method-development R Markdown notebooks
  (distribution fitting, metalog methods, prototype driver notebooks). Not
  run by the pipeline; kept for their intellectual content and as potential
  paper material.
- **`coweeta/`** - The Coweeta Hydrologic Laboratory watershed case study,
  a self-contained sub-project.
- **`scratch/`** - Superseded, one-off, or broken scripts kept for historical
  reference only. Nothing here is sourced by the pipeline or any notebook.
- **`docs/`** - Non-analysis documentation (the original rrtools/workflowr
  project-scaffolding runbook).

## Deduplication done during the reorg

Several functions existed as 2-3 near-identical copies scattered across
notebooks. The following were extracted into a single canonical location in
`code/R/`, with the old inline definitions replaced by `source()` calls:

- `bayesian_update()` (grid-KDE Bayesian updating) -> `code/R/bayesian_updating.R`
- `generate_inverse_cdf_distribution()` and `calculate_summary_statistics()`
  (from `code/notebooks/soilgrids_simulation.Rmd`, the cleanest/most
  defensive versions found) -> `code/R/distribution_fitting.R`
- The full depth-simulation and correlated-property-simulation function
  families from `code/coweeta/sim-functions.R` -> split into
  `code/R/depth_simulation.R` and `code/R/property_simulation.R`
- `download_ssurgo_tabular()` and the soil-infill/imputation cascade were
  already consolidated in `soil_infill_functions.R`, which simply moved to
  `code/R/` and is now wired into `_targets.R` for the first time (it was
  previously written but never sourced by the pipeline)

Superseded copies (`code/scratch/local_functions.R`,
`code/scratch/local_functions_testing.R`, `code/scratch/bayesian_updating_normal.R`,
`code/scratch/spline_inv_cdf.R`) were moved to `scratch/` rather than deleted.

## Open TODOs (not addressed by this reorg - require scientific/methods judgment)

- **Metalog fitting internals** are reimplemented from scratch independently
  in `code/notebooks/metalog_functions.Rmd`, `code/notebooks/hwsd_wise_metalog_dist.Rmd`,
  and `code/notebooks/simulating_correlated_soil_properties_metalog.Rmd`. No
  canonical version has been picked - this requires judgment about which
  implementation is most correct, not a mechanical move.
- **`calculate_summary_statistics()`** has two more redefinitions, in
  `code/notebooks/hwsd_wise_metalog_dist.Rmd` and
  `code/notebooks/generating_empirical_parametric_distributions.Rmd`, that
  were left in place rather than merged into the canonical
  `code/R/distribution_fitting.R` version - they are not byte-identical
  (different signatures/behavior), so merging them is a methods decision.
- **`adjust_depthwise_property_GP_Quant()`** exists in two diverging versions:
  the one in `code/R/property_simulation.R` (used by the pipeline) and an
  earlier ancestor in `code/coweeta/depthwise-multivariate-profile-simulation.Rmd`
  with an added `sd_adjustment_factor` term. Not merged - see
  `code/coweeta/README.md`.
- **`simulate_cokey()`** in `code/R/property_simulation.R` appears to be an
  earlier, superseded implementation of the same logic as
  `simulate_cokey_generalized()` and may be dead code. Kept as-is pending a
  follow-up review.
- **Portability**: many notebooks (mostly in `code/notebooks/` and
  `code/coweeta/`) hardcode machine-specific absolute paths
  (`C:/R_Drive/...`, `C:/Python310/...`, `e:/gis_data/...`,
  `C:/LandPKS_API_SoilID-master/...`). The recommended pattern going forward
  is `here::here()` for in-repo paths and new entries in `config.yml` (which
  already holds `aoi_wkt`/`properties`) for paths to external data sources,
  rather than hardcoding them per-notebook.
- **No dependency lockfile.** `.binder/`'s Docker build was removed because it
  called `renv::restore()` against a lockfile that never existed. If
  reproducible-environment tooling becomes a priority, run `renv::init()` to
  generate `renv.lock` as a follow-up.
- **Known-broken code left flagged but unfixed** (see inline `# BROKEN:`
  comments): `code/scratch/compare-texture.R` sources a nonexistent
  `code/SOLUS/local-functions.R`; `code/coweeta/coweeta_soil_simulation.Rmd`
  calls an undefined `KSSL_VG_model()`; `code/coweeta/AWS-Map-Simulation.Rmd`'s
  GP-diagnostic section references undefined variables copy-pasted from
  `depthwise-multivariate-profile-simulation.Rmd`;
  `code/scratch/local_functions_testing.R` has a wrong-arity function call;
  `code/notebooks/simulating_correlated_soil_properties_metalog.Rmd` sources
  an external GitHub gist and calls `fit_dist()` on a non-dataframe.
- **`code/scratch/soil_infill_functions.py`** is the only non-R file in the
  repo and has no import/call relationship with any R file here. It likely
  belongs in a separate Python codebase (e.g. LandPKS SoilID) rather than
  this repo.
