# code/coweeta/

The Coweeta Hydrologic Laboratory watershed case study - a self-contained
application of the project's soil-property simulation methods to a specific
watershed.

- **`AWS-Map-Simulation.Rmd`** - main Coweeta available-water-storage
  simulation driver. Runs the depth/property simulation pipeline end-to-end
  and saves a key intermediate result (`data/coweeta-soil-sim-data.rda`). Its
  `source()` calls now point at `code/R/depth_simulation.R`,
  `code/R/property_simulation.R`, and `code/R/soil_infill_functions.R`
  (fixed during the code/ reorganization - it previously pointed at a
  `code/sim-functions.R` path that never resolved). The back half of this
  notebook (a GP-model diagnostic section) references variables that are
  never defined in this file - see the `# BROKEN:` comment at that chunk.
- **`SoilGrids_POLARIS_Download.Rmd`** - data-acquisition notebook (SoilGrids/
  POLARIS downloads via GDAL vsicurl) for Coweeta, which then pivots to an
  unrelated iSDAsoil STAC-catalog exploration for Rwanda/DRC.
- **`coweeta_soil_simulation.Rmd`** - the most complete, publication-oriented
  Coweeta pipeline (produces PDFs and rayshader 3D block diagrams). Duplicates
  significant simulation logic from `code/R/property_simulation.R` rather than
  reusing it, and calls an undefined `KSSL_VG_model()` function (flagged
  in-line, not fixed).
- **`depthwise-multivariate-profile-simulation.Rmd`** - pure method-development
  notebook (100% synthetic toy data) for the depth-wise Gaussian Process
  correction idea. This is the ancestor that the GP-diagnostic code in
  `AWS-Map-Simulation.Rmd` was copy-pasted from (without its setup code).

## Known, unresolved divergence

`adjust_depthwise_property_GP_Quant()` exists in two versions that have
diverged from a common ancestor:

- `code/R/property_simulation.R` (used by `_targets.R` and
  `AWS-Map-Simulation.Rmd`)
- `depthwise-multivariate-profile-simulation.Rmd` (uses
  `predict.GP(..., newdata = ...)` and an added `sd_adjustment_factor` term
  not present in the other version)

These were intentionally left unmerged during the code/ reorganization -
picking which behavior is correct is a scientific-methods decision, not a
file-move decision.
