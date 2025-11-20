# BCQRF Code Organization Guide

This document categorizes the codebase into **core implementation code** (reusable function libraries) versus **example/demonstration code** (workflows, notebooks, and case studies).

---

## Core Implementation Code (Function Libraries)

These are the foundational function libraries that implement the BCQRF methodology. They contain reusable functions that should be sourced/imported by other scripts.

### Primary Function Libraries

#### 1. `soil_infill_functions.R` (103 KB) ⭐ **MOST IMPORTANT**
**Purpose**: Primary data processing engine for soil property infilling

**Sections**:
- Section 1: Main workflow function (`infill_soil_property()`)
- Section 2: Horizon suitability functions (identify unsuitable horizons)
- Section 3: Data cleaning functions
- Section 4: Property configuration functions
- Section 5: Range value infilling functions
- Section 6: **6-tier hierarchical data recovery strategy**
- Section 7: Special functions (RFV handling, SSURGO download)

**Key Innovation**: Hierarchical infilling that respects pedological structure

**Use Case**: Import this when you need to fill missing soil property data from SSURGO or other sources

---

#### 2. `local_functions.R` (26 KB)
**Purpose**: Soil profile depth simulation and horizon boundary modeling

**Key Functions**:
- `get_aws_data_by_mukey()` - Query SSURGO by map unit key
- `tri_dist()` - Triangular distribution sampling
- `simulate_soil_profile_top_down()` - Top-down depth simulation
- `simulate_soil_profile_bottom_up()` - Bottom-up depth simulation
- `simulate_soil_profile_thickness()` - Monte Carlo thickness estimation (500 realizations)
- `simulate_and_perturb_soil_profiles()` - Profile perturbation using aqp::perturb()
- `query_osd_distinctness()` - OSD horizon boundary distinctness
- `infill_missing_distinctness()` - Fill missing boundary data

**Use Case**: Import when you need to simulate realistic soil horizon depths and boundaries

---

#### 3. `coweeta/sim-functions.R` (100 KB)
**Purpose**: Site-specific multivariate simulation functions for Coweeta LTER case study

**Key Functions**:
- `remove_organic_layer()` - Pre-processing to exclude O horizons
- `slice_and_aggregate_soil_data()` - Depth-weighted aggregation
- `sim_component_comp()` - Component composition simulation
- `simulate_correlated_triangular()` - Multivariate triangular sampling with correlations
- `adjust_depthwise_property_GP_Quant()` - Gaussian process adjustment
- `adjust_soil_data_parallel()` - Parallelized multi-property adjustment
- `simulate_soil_properties()` - Master simulation function
- `van_genuchten()` - Water retention model
- `simulate_vg_aws()` - Plant-available water storage simulation
- `calculate_aws_df()` - AWS aggregation to 100 cm

**Use Case**: Import for advanced multivariate simulation with spatial adjustment, especially for hydrological applications

---

#### 4. `soil_infill_functions.py` (34 KB)
**Purpose**: Python implementation of soil property infilling

**Key Functions**:
- `infill_soil_data(df)` - Main workflow for pandas DataFrames
- `infill_range_values(df)` - Range value infilling (_l and _h)
- `infill_parameter_range()` - Multi-strategy parameter infilling
- `infill_missing_texture_data()` - Texture-specific recovery

**Use Case**: Use when working in Python environments or integrating with Python-based ML pipelines

---

### Utility Function Libraries

#### 5. `spline_inv_cdf.R` (1.4 KB)
**Purpose**: Inverse CDF functions for spline-based quantile estimation

**Use Case**: Utility for converting between quantiles and values using spline interpolation

---

#### 6. `compare-texture.R` (6.1 KB)
**Purpose**: Soil texture classification utilities

**Use Case**: Classify soil texture classes from sand/silt/clay percentages

---

#### 7. `bayesian_updating_normal.R` (6.1 KB)
**Purpose**: Bayesian constraining framework implementations

**Methods**:
1. Normal distributions (analytical updates)
2. Biased likelihood with variance inflation
3. Non-normal distributions using kernel density estimation

**Use Case**: Constrain model predictions using empirical priors (Prior × Likelihood → Posterior)

---

## Example & Demonstration Code

These files demonstrate how to use the core functions, provide tutorials, and show complete workflows.

### Workflow Demonstration Scripts

#### 1. `assemble-pieces.R` (8.3 KB) - **Example Workflow**
**Purpose**: Demonstrates multi-source data integration workflow

**What it does**:
- Downloads soil data from SSURGO, SoilGrids, SOLUS
- Creates rasterized soil property maps
- Compares different data sources spatially
- Demonstrates `terra`, `soilDB` integration

**Type**: Executable example with specific geographic areas (CA, OK, IL, NE)

---

#### 2. `local_functions_testing.R` (30 KB) - **Testing/Development**
**Purpose**: Unit tests and validation for `local_functions.R`

**Type**: Testing and development script, not for production use

---

### R Markdown Notebooks (Tutorials & Demonstrations)

All `.Rmd` files are **demonstration notebooks** that illustrate methods and generate HTML reports.

#### Metalog Distribution Notebooks

##### 1. `Metalog_functions.Rmd` (27 KB)
**Demonstrates**:
- Metalog distribution fitting and sampling
- Basis function implementation (6 terms)
- `normalize()`, `denormalize()`, `fit_metalog()`, `recreate_metalog()`
- Using `rmetalog` package for bounded distributions

**Learning Goal**: Understand how to fit flexible bounded distributions to soil property data

---

##### 2. `Simulating correlated soil property values from empirical distributions modeled with metalog distributions using R.Rmd` (19 KB) - **KEY TUTORIAL**
**Demonstrates**:
- Complete multivariate simulation workflow
- Fitting metalog to each property independently
- Correlation matrix and Cholesky decomposition
- Simulating correlated soil properties
- Validation that correlations are preserved

**Learning Goal**: Master the complete workflow for multivariate correlated soil property simulation

---

##### 3. `KSSL_soil_property_correlation_matricies.Rmd` (16 KB)
**Demonstrates**:
- Computing empirical correlation matrices from KSSL database
- Compositional data analysis (ilr transformation for texture)
- Stratification by generalized horizon type (A, B, C, O, R, Cr)
- Database queries with SQLite

**Learning Goal**: Generate empirical correlation matrices for your own regions/soil types

---

##### 4. `Generating_Empirical_Parametric_Distributions.Rmd` (39 KB)
**Demonstrates**:
- Creating empirical distributions from summary statistics
- Iterative density adjustment to match statistical moments
- Piecewise-linear CDF construction

**Learning Goal**: Generate distributions when you only have summary stats (min, max, mean, mode, quantiles)

---

##### 5. `emprical_dist_testing.Rmd` (21 KB)
**Demonstrates**:
- Validation of empirical distribution generation
- Comparison with triangular and metalog distributions
- Statistical validation (mean, mode, quantiles)

**Learning Goal**: Understand trade-offs between distribution types

---

##### 6. `HWSD_WISE_Metalog_Dist.Rmd` (39 KB)
**Demonstrates**:
- Applying metalog methods to global soil databases (HWSD, WISE)
- Download → Extract → Fit metalog → Simulate workflow

**Learning Goal**: Scale methods to global datasets

---

#### Process & Workflow Notebooks

##### 7. `bcqrf_process.Rmd` (13 KB) - **END-TO-END WORKFLOW**
**Demonstrates**:
- Complete BCQRF implementation from start to finish
- Integration of data sources, infilling, simulation, validation
- Property lookup tables for SoilGrids, SOLUS, SSURGO

**Learning Goal**: See the complete BCQRF workflow in action

---

##### 8. `SoilGrids_simulation.Rmd` (56 KB)
**Demonstrates**:
- Processing SoilGrids soil profile object data
- Reshaping horizons to long format
- Extrapolating Q05/Q95 to min/max bounds
- Generating density-based distributions from quantile data

**Learning Goal**: Work with SoilGrids data and quantile-based predictions

---

##### 9. `constrained_QRF_ideas.Rmd` (18 KB)
**Demonstrates**:
- Theoretical framework for constrained quantile regression forest
- Implementation ideas and considerations

**Learning Goal**: Understand the QRF constraining methodology (conceptual, not fully implemented)

---

### Coweeta LTER Case Study Notebooks

These demonstrate the complete application of BCQRF methods to a real research site.

##### 1. `coweeta/coweeta_soil_simulation.Rmd` (30 KB) - **COMPLETE CASE STUDY**
**Demonstrates**:
- Multi-source data integration (SSURGO + RSS)
- SSURGO vs RSS comparison
- Compositional texture simulation with ilr transformation
- Multivariate simulation using `MASS::mvrnorm()`
- van Genuchten water retention parameters

**Learning Goal**: See how to apply all methods to a real research site

---

##### 2. `coweeta/depthwise-multivariate-profile-simulation.Rmd` (29 KB)
**Demonstrates**:
- Depth-wise (horizon-by-horizon) multivariate simulation
- Separate correlation structures for each depth interval

**Learning Goal**: Handle depth-dependent correlations

---

##### 3. `coweeta/SoilGrids_POLARIS_Download.Rmd` (21 KB)
**Demonstrates**:
- Data acquisition from SoilGrids and POLARIS systems
- Download, reproject, format workflow

**Learning Goal**: Acquire and pre-process external soil data sources

---

##### 4. `coweeta/AWS-Map-Simulation.Rmd` (26 KB)
**Demonstrates**:
- Plant-available water storage (AWS) simulation and mapping
- Hydrological modeling application
- Quantile-based mapping (5%, 50%, 95%)

**Learning Goal**: Apply simulated soil properties to hydrological applications

---

## Recommended Usage Patterns

### For Method Development/Research
Start with these **core libraries**:
1. `soil_infill_functions.R` - Data infilling
2. `local_functions.R` - Depth simulation
3. `bayesian_updating_normal.R` - Bayesian constraining

Source these in your own scripts as needed.

### For Learning the Methods
Follow this **tutorial sequence**:
1. `Metalog_functions.Rmd` - Learn metalog basics
2. `KSSL_soil_property_correlation_matricies.Rmd` - Generate correlation matrices
3. `Simulating correlated soil property values...Rmd` - Multivariate simulation
4. `bcqrf_process.Rmd` - See complete workflow
5. `coweeta/coweeta_soil_simulation.Rmd` - Complete case study

### For Adapting to Your Own Site
Use this **workflow**:
1. Start with `bcqrf_process.Rmd` as template
2. Import functions from `soil_infill_functions.R` and `local_functions.R`
3. Generate your own correlation matrices using `KSSL_soil_property_correlation_matricies.Rmd` approach
4. Follow `coweeta/coweeta_soil_simulation.Rmd` structure for your site
5. Adapt `coweeta/sim-functions.R` for site-specific needs

### For Python Users
1. Use `soil_infill_functions.py` for data infilling
2. Translate other functions to Python as needed
3. Can use R functions via `rpy2` if preferred

---

## Quick Reference: "Should I Modify This File?"

| File Type | Modify? | Purpose |
|-----------|---------|---------|
| `*_functions.R` | **YES** | Core libraries - add functions here |
| `*_functions.py` | **YES** | Core Python library - add functions here |
| `assemble-pieces.R` | **NO** | Example script - copy and adapt instead |
| `*_testing.R` | **NO** | Testing script - for validation only |
| `*.Rmd` | **NO** | Tutorials - copy chunks to your own notebooks |
| `coweeta/*.Rmd` | **NO** | Case study - use as template for your site |

---

## Function Import Examples

### R Example
```r
# Source core function libraries
source("code/soil_infill_functions.R")
source("code/local_functions.R")
source("code/bayesian_updating_normal.R")

# Now use the functions
my_data <- get_aws_data_by_mukey(c("123456", "789012"))
infilled_data <- infill_soil_property(my_data, "claytotal", get_default_property_config("clay"))
simulated_depths <- simulate_soil_profile_thickness(my_spc, n = 500)
```

### Python Example
```python
import pandas as pd
from code.soil_infill_functions import infill_soil_data, infill_range_values

# Load your data
df = pd.read_csv("my_soil_data.csv")

# Infill missing data
df_clean = infill_soil_data(df)
df_with_ranges = infill_range_values(df_clean)
```

---

## Summary Statistics

| Category | Count | Total Size |
|----------|-------|------------|
| **Core Function Libraries** | 7 files | 271 KB |
| **Example/Demo Scripts** | 2 files | 38 KB |
| **Tutorial Notebooks** | 9 files | 254 KB |
| **Case Study Notebooks** | 4 files | 106 KB |
| **Total** | 22 files | ~669 KB |

---

**Created**: 2025-11-20
**Purpose**: Guide users to understand code organization and usage patterns
