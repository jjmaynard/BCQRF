# Detailed Summary: BCQRF `/code/` Folder

## Overview

The **BCQRF** (Bayesian Constrained Quantile Regression Forest) `/code/` folder contains a comprehensive soil science research framework (~12,800 lines of code) that implements sophisticated methods for **soil property simulation, data infilling, and uncertainty quantification**. The codebase combines R and Python to bridge pedological knowledge, statistical methodology, and computational modeling.

---

## Primary Research Goals

1. **Download and process soil data** from multiple sources (SSURGO, SOLUS, SoilGrids, POLARIS, KSSL)
2. **Infill missing soil property data** using hierarchical, pedologically-aware strategies
3. **Simulate correlated soil properties** using metalog distributions with empirical correlations
4. **Apply Bayesian methods** to constrain model predictions with empirical data
5. **Generate realistic soil profile realizations** with depth variability and uncertainty quantification
6. **Produce spatially explicit predictions** for soil mapping and hydrological applications

---

## File Structure

### Main R Scripts (7 files)

- **`local_functions.R`** (26 KB) - Soil profile and depth simulation
- **`soil_infill_functions.R`** (103 KB) - **Primary data processing engine**
- **`soil_infill_functions.py`** (34 KB) - Python companion implementation
- **`bayesian_updating_normal.R`** (6 KB) - Bayesian constraining framework
- **`assemble-pieces.R`** (8 KB) - Multi-source workflow orchestration
- **`compare-texture.R`** (6 KB) - Texture classification utilities
- **`spline_inv_cdf.R`** (1 KB) - Inverse CDF functions

### R Markdown Notebooks (11 files)

**Metalog Distribution Notebooks:**
- `Metalog_functions.Rmd` - Basic metalog fitting and sampling
- `Simulating correlated soil property values from empirical distributions modeled with metalog distributions using R.Rmd` - Multivariate simulation workflow
- `KSSL_soil_property_correlation_matricies.Rmd` - Empirical correlation matrices
- `Generating_Empirical_Parametric_Distributions.Rmd` - Synthetic distributions from statistics

**Simulation Notebooks:**
- `SoilGrids_simulation.Rmd` - SoilGrids data processing
- `HWSD_WISE_Metalog_Dist.Rmd` - Global soil database applications
- `emprical_dist_testing.Rmd` - Distribution validation
- `constrained_QRF_ideas.Rmd` - Constrained QRF methodology

**Process Notebooks:**
- `bcqrf_process.Rmd` - End-to-end BCQRF workflow

### Coweeta Subdirectory

**Functions:**
- `coweeta/sim-functions.R` (100 KB) - Site-specific simulation functions

**Notebooks:**
- `coweeta/coweeta_soil_simulation.Rmd` - SSURGO vs RSS comparison
- `coweeta/depthwise-multivariate-profile-simulation.Rmd` - Depth-wise multivariate simulation
- `coweeta/SoilGrids_POLARIS_Download.Rmd` - Data acquisition
- `coweeta/AWS-Map-Simulation.Rmd` - Hydrological application

---

## Core Computational Methods

### 1. Hierarchical Soil Property Infilling

**Implementation:** `soil_infill_functions.R`

**6-tier priority system** for recovering missing data:

1. **Horizon Name Matching** - Find pedologically similar horizons (same designation)
2. **Depth-Weighted Infilling** - Use similar depth ranges with exponential distance weighting
3. **Within-Component Interpolation** - Fit vertical trends within profiles
4. **Cross-Component Interpolation** - Use other components in same map unit
5. **Related Property Estimation** - Apply soil science relationships (e.g., clay ↔ CEC)
6. **Statistical Fallback** - Group means as last resort

**Key Innovation**: Respects pedological structure while salvaging partial data, explicitly excluding unsuitable horizons (bedrock, cemented layers, organic).

#### Main Function Categories

**Section 1: Main Function**
- `infill_soil_property(df, property_name, property_config)` - Master workflow

**Section 2: Horizon Suitability**
- `identify_unsuitable_horizons()` - Detects R, Cr, O horizons
- `is_unsuitable_horizon_type()` - Pattern matching for cemented/restricted horizons
- `filter_unsuitable_horizons()` - Adds unsuitable_horizon flag

**Section 3: Data Cleaning**
- `clean_property_data()` - Validates and cleans soil property columns
- `apply_basic_range_limits()` - Enforces pedologically reasonable bounds

**Section 4: Property Configuration**
- `get_default_property_config()` - Configs for common properties
- `create_custom_property_config()` - Customizable specifications

**Section 5: Range Value Infilling**
- `infill_property_range_values()` - Primary range infilling
- `learn_property_ranges()` - Extracts variability patterns
- `get_property_contextual_ranges()` - Pedological knowledge

**Section 6: Property Data Recovery**
- `infill_missing_property_data()` - Orchestrates 6-tier strategy
- `horizon_name_property_infill()` - Tier 1
- `depth_weighted_property_infill()` - Tier 2
- `within_component_property_interpolation()` - Tier 3
- `cross_component_property_interpolation()` - Tier 4
- `related_property_estimation()` - Tier 5

**Section 7: Special Functions**
- `infill_rfv_property()` - Rock fragment volume handling
- `download_ssurgo_tabular()` - SSURGO download with restriction detection

---

### 2. Metalog Distribution Modeling

**Why Metalog?**
- Flexible 6-term basis functions accommodate diverse shapes
- Bounded distributions (0-100% for percentages)
- Efficient fitting via linear regression (not iterative optimization)
- Maintains marginal distributions under transformation

**Basis Functions:**
```r
1. log(y)
2. log(1-y)
3. y(1-y)log(y)
4. y(1-y)log(1-y)
5. y(1-y)[log(y)]²
6. y(1-y)[log(1-y)]²
```

**Key Functions:**
```r
normalize(x, min, max)           # Scale to (0,1)
denormalize(y, min, max)         # Scale back to original range
metalog_basis(y)                 # Calculate basis functions
fit_metalog(p, y)                # Fit via linear regression
recreate_metalog(coef, n)        # Generate samples
```

**Applications:**
- Fit metalog to each soil property independently
- Generate realistic bounded distributions from summary statistics
- Enable multivariate simulation with correlation preservation

**Implementation Files:**
- `Metalog_functions.Rmd` - Core functions
- `Simulating correlated soil property values...Rmd` - Multivariate use
- `Generating_Empirical_Parametric_Distributions.Rmd` - From summary stats

---

### 3. Multivariate Correlated Soil Property Simulation

**Workflow:**

1. **Fit metalog distribution** to each property independently
2. **Extract empirical correlation matrix** from KSSL database (stratified by horizon type)
3. **Apply Cholesky decomposition** to correlation matrix
4. **Generate uncorrelated samples** from metalog distributions
5. **Apply Cholesky transformation** to induce correlations
6. **Back-transform** to original scale

**Mathematical Framework:**
```r
# For each property i:
M_i = fit_metalog(quantiles_i)

# Empirical correlation matrix:
R = cor(properties_matrix)

# Cholesky decomposition:
L = chol(R)

# Uncorrelated samples:
U = matrix(recreate_metalog(M_i, n) for each i)

# Apply correlation:
X = U %*% L

# Result: correlated samples preserving marginal distributions
```

**Result**: Simulated soil properties with:
- Correct marginal distributions
- Correct empirical correlations
- Physically plausible combinations

**Implementation:** `Simulating correlated soil property values from empirical distributions modeled with metalog distributions using R.Rmd`

---

### 4. Compositional Data Analysis

**Challenge**: Sand/silt/clay percentages must sum to 100%

**Solution**: **Isometric log-ratio (ilr) transformation**

**Why ilr?**
- Transforms 3 compositional parts to 2 unconstrained dimensions
- Preserves distances and correlations
- Allows standard statistical operations
- Back-transformation preserves sum=100% constraint

**Workflow:**
```r
library(compositions)

# Transform to ilr space:
texture_ilr <- ilr(acomp(cbind(sand, silt, clay)))

# Calculate correlations in ilr space:
cor_matrix <- cor(texture_ilr)

# Simulate in ilr space:
sim_ilr <- mvrnorm(n, mu, Sigma)

# Back-transform to percentages:
sim_texture <- ilrInv(sim_ilr)
# sum(sim_texture) == 100 (guaranteed)
```

**Implementation:** `KSSL_soil_property_correlation_matricies.Rmd`

---

### 5. Depth Variability Simulation

**Implementation:** `local_functions.R`

**Method**: Monte Carlo simulation (500 realizations) using:
- **Triangular distributions** from SSURGO low/representative/high depth values
- **Top-down** simulations (start at surface, cumulate down)
- **Bottom-up** simulations (start at bedrock, cumulate up)
- **OSD distinctness codes** → boundary variability (Abrupt/Clear/Gradual/Diffuse)

**Key Functions:**

```r
tri_dist(n, a, b, c)  # Triangular distribution sampling
                      # a=min, b=max, c=mode

simulate_soil_profile_top_down(spc)
  # Start at surface, add horizon thicknesses sequentially

simulate_soil_profile_bottom_up(spc)
  # Start at bedrock, add horizons upward

simulate_soil_profile_thickness(spc, n=500)
  # Monte Carlo estimation of horizon thickness variability
  # Alternates top-down and bottom-up randomly

simulate_and_perturb_soil_profiles(spc, n=25)
  # Uses aqp::perturb() with thickness and boundary perturbation
  # Boundary perturbation based on OSD distinctness
```

**OSD Distinctness to Boundary Variability:**
| Distinctness | Boundary Thickness | Std Dev (cm) |
|--------------|-------------------|--------------|
| Abrupt       | < 2.5 cm          | 0.5          |
| Clear        | 2.5-6.5 cm        | 1.5          |
| Gradual      | 6.5-12.5 cm       | 3.0          |
| Diffuse      | > 12.5 cm         | 5.0          |

**Output**: ~25 simulated profiles per original profile with realistic depth uncertainty

---

### 6. Bayesian Constraining

**Framework**: Prior (empirical) × Likelihood (model) → Posterior (constrained)

**Mathematical Foundation:**
```
Posterior ∝ Prior(θ) × Likelihood(Data|θ)
P(θ|Data) ∝ P(Data|θ) × P(θ)
```

**Three Implementations** (bayesian_updating_normal.R):

#### 1. Normal Distributions (Analytical)
```r
# Prior: N(μ₀, σ₀²)
# Likelihood: N(μ₁, σ₁²)
# Posterior: N(μ_post, σ_post²)

σ_post² = 1 / (1/σ₀² + 1/σ₁²)  # Harmonic mean
μ_post = σ_post² × (μ₀/σ₀² + μ₁/σ₁²)  # Weighted average
```

#### 2. Biased Likelihood (Variance Inflation)
```r
# If model is biased, inflate variance:
σ₁² = σ₁² × inflation_factor  # e.g., 4x

# Reduces likelihood influence proportionally
# Prevents overconfident posterior
```

#### 3. Non-Normal Distributions (Kernel Density)
```r
# Estimate densities:
prior_density <- density(empirical_data)
likelihood_density <- density(model_predictions)

# Apply Bayes' rule pointwise:
posterior_density <- prior_density * likelihood_density

# Normalize:
posterior_density <- posterior_density / sum(posterior_density)

# Sample from posterior
```

**Application**: Constrain biased model predictions (SoilGrids, SOLUS QRF) using empirical SSURGO data as prior

---

### 7. Correlation Matrices by Horizon Type

**Implementation:** `KSSL_soil_property_correlation_matricies.Rmd`

**Data Processing:**

1. **Query SQLite KSSL database** for lab data
2. **Calculate bulk density** if missing: ρ_b = (1-porosity) × ρ_particle
3. **Compute rock fragment volume**: RFV = (ρ_particle - ρ_b) / ρ_particle × 100
4. **Extract properties**: sand, silt, clay, bulk density, water retention, CEC, pH, bases

**Compositional Analysis:**
- Sand/Silt/Clay → ilr transformation
- Calculate correlations in ilr space (2 dimensions)
- Store as correlation matrices

**Stratification:**
- Separate correlation matrices for each generalized horizon type:
  - **A horizon** - Surface, organic-rich
  - **B horizon** - Subsoil, clay accumulation
  - **C horizon** - Parent material
  - **O horizon** - Organic
  - **R horizon** - Bedrock
  - **Cr horizon** - Weathered bedrock

**Why Stratify?**
- Clay-CEC correlation stronger in B horizons
- Bulk density patterns differ by horizon
- Captures pedological differences

**Output**:
- One correlation matrix per generalized horizon
- Used in multivariate simulation via Cholesky decomposition

---

## Data Flow Architecture

```
SSURGO Database
    ↓
download_ssurgo_tabular()
    ↓
[Horizon restriction detection]
    ↓
infill_soil_property() [main workflow]
    │
    ├─→ [unsuitable horizon filtering] ← OSD distinctness data
    │
    ├─→ [range value infilling]
    │   └─→ [learns from data + contextual knowledge]
    │
    └─→ [missing data recovery] (6-tier priority)
        ├─→ [horizon matching]
        ├─→ [depth weighting]
        ├─→ [within-profile trends]
        ├─→ [cross-component analogs]
        ├─→ [property relationships]
        └─→ [group means fallback]
    ↓
[Clean soil property estimates with bounds]
    ↓
METALOG FITTING (per property, per horizon)
    ↓
CORRELATION MATRIX (empirical from KSSL)
    ↓
CHOLESKY DECOMPOSITION
    ↓
MULTIVARIATE SIMULATION
    ├─→ [uncorrelated metalog samples]
    └─→ [Cholesky-constrained correlation]
    ↓
[Simulated soil profiles with correlated properties]
    ↓
[Depth-wise simulations] ← local_functions.R simulations
    ↓
BAYESIAN CONSTRAINING
[Prior (empirical) × Likelihood (model) → Posterior]
    ↓
[Final constrained predictions]
```

---

## Coweeta LTER Case Study

**Purpose**: Complete applied workflow for Coweeta Hydrologic Laboratory site (North Carolina)

**Data Integration**:
- **SSURGO** - Ground-truthed survey data
- **RSS** - Rapid Soil Survey (areas lacking SSURGO)
- **SoilGrids** - Global gridded 3D soil properties
- **POLARIS** - Probabilistic Observatories for Land and Aquatic Research

**Workflow Architecture:**

1. **Download Phase** (`SoilGrids_POLARIS_Download.Rmd`)
   - Acquire raw data from multiple sources
   - Reproject to common CRS
   - Format for downstream analysis

2. **Processing Phase** (`coweeta_soil_simulation.Rmd`)
   - Split into SSURGO vs RSS components
   - Extract compositional texture data (sand/silt/clay)
   - Apply ilr transformation to preserve sum=100%
   - Calculate correlation matrices per generalized horizon

3. **Simulation Phase** (`depthwise-multivariate-profile-simulation.Rmd`)
   - Simulate (comppct_r × 10) realizations per component
   - Uses multivariate normal (`MASS::mvrnorm`) with Cholesky correlation
   - Back-transform ilr → texture percentages
   - Depth-specific correlation structures

4. **Application Phase** (`AWS-Map-Simulation.Rmd`)
   - Compute van Genuchten water retention parameters
   - Calculate plant-available water storage (AWS) to 100 cm
   - Summarize quantiles (5%, 50%, 95%) by map unit and depth
   - Generate maps for hydrological modeling

**Key Functions in `sim-functions.R`:**

| Function | Purpose |
|----------|---------|
| `remove_organic_layer()` | Pre-processing to exclude O horizons |
| `slice_and_aggregate_soil_data()` | Depth-weighted averaging to standard intervals |
| `sim_component_comp()` | Component composition simulation |
| `simulate_correlated_triangular()` | Multivariate triangular sampling |
| `adjust_depthwise_property_GP_Quant()` | Gaussian process adjustment |
| `adjust_soil_property_by_depth()` | Depth-aware property adjustment |
| `adjust_soil_data_parallel()` | Multi-property adjustment with parallelization |
| `simulate_soil_properties()` | Master simulation function |
| `van_genuchten()` | Water retention model |
| `simulate_vg_aws()` | Plant-available water simulation |
| `calculate_aws_df()` | AWS aggregation to 100 cm |

**Innovation**: SSURGO vs RSS comparison at same site reveals data source impacts on uncertainty quantification

**Output**: Quantile-based AWS maps (5%, 50%, 95%) for hydrological modeling and drought assessment

---

## Python Implementation

**File:** `soil_infill_functions.py` (34 KB)

**Purpose**: Data-driven robust infilling with pandas, numpy, scipy foundation

### Core Functions

#### `infill_soil_data(df)`
- Processes entire DataFrame with horizon-level filtering
- Identifies problematic horizons (missing texture data in top 50 cm)
- Attempts infilling or marks as incomplete
- Returns DataFrame with `texture_data_complete` flag

#### `infill_range_values(df)`
Robust `_l` and `_h` value infilling using:
- **Learned ranges**: Extract variability from complete data (by horizon, component, depth zone)
- **Contextual ranges**: Pedological knowledge lookup
- **Fallback ranges**: Conservative defaults

Parameter groups handled:
- **Texture** (sand, clay, silt) with sum-to-100% constraint
- **Bulk density** (0.8-2.2 g/cm³ range)
- **Water retention** (0-60% range)

#### `infill_parameter_range(df, param, config, group_name)`
Multi-strategy parameter infilling:
1. Learn from existing complete ranges
2. Apply contextual ranges
3. Infill missing `_l` values
4. Infill missing `_h` values
5. Enforce range constraints

#### Helper Functions
- `infill_missing_texture_data()` - Horizon-level texture recovery
- `learn_parameter_ranges()` - Extract variability patterns
- `get_contextual_ranges()` - Pedological knowledge
- `calculate_lower_bound()` - Context-aware lower estimation
- `calculate_upper_bound()` - Context-aware upper estimation
- `enforce_range_constraints()` - Final validation

**Design Philosophy**: Horizon-level processing rather than group exclusion allows salvage of partial data while flagging incomplete cases.

**Advantage over R**: Efficient pandas operations for large datasets, easier integration with Python-based ML workflows

---

## Major Dependencies

### R Packages

**Core Soil Science:**
- `soilDB` - SSURGO/NRCS database queries (SDA_query, mukey.wcs, fetchOSD)
- `aqp` - Soil Profile Collection (SPC) object model, perturb(), plotSPC()
- `terra` - Modern spatial raster/vector operations (replaces raster/sp)

**Statistical & Modeling:**
- `MASS` - `mvrnorm()` for multivariate simulation
- `compositions` - Compositional data analysis (ilr, acomp, ilrInv)
- `rmetalog` - Metalog distribution fitting and sampling
- `DBI`, `RSQLite` - Database connectivity
- `truncnorm` - Truncated normal distributions

**Data Manipulation:**
- `dplyr` - Grammar of data manipulation (filter, mutate, group_by, summarize)
- `tidyr` - Data reshaping (pivot_longer, pivot_wider)
- `stringr` - String manipulation
- `reshape2` - Legacy reshaping (melt/cast)

**Visualization & Reporting:**
- `ggplot2` - Publication-quality graphics
- `lattice`, `latticeExtra` - Trellis graphics for soil profiles
- `leaflet` - Interactive maps

**Utilities:**
- `here` - Portable file paths
- `roxygen2` - Documentation generation

### Python Packages

- `pandas` - DataFrame operations
- `numpy` - Numerical computing
- `scipy` - Scientific computing (optimization, statistics)

---

## File Dependency Graph

```
MAIN WORKFLOWS:
├── local_functions.R (26 KB)
│   └── Uses: aqp, soilDB
│   ├── → simulate_soil_profile_depth simulations (used in coweeta/)
│   └── → depth variability quantification
│
├── soil_infill_functions.R [103 KB - LARGEST]
│   ├── Requires: dplyr, terra, soilDB
│   ├── Section 1: Main entry point (infill_soil_property)
│   ├── Section 2-5: Data processing pipeline
│   ├── Section 6: Priority-based recovery (6-tier)
│   └── Section 7: Special case handling (RFV)
│
├── soil_infill_functions.py (34 KB)
│   └── Companion implementation in Python (pandas-based)
│
├── bayesian_updating_normal.R (6 KB)
│   └── Framework for prior×likelihood→posterior
│
└── METALOG SUBSYSTEM:
    ├── Metalog_functions.Rmd
    │   └── Basic metalog fitting (rmetalog package)
    ├── Simulating correlated soil property...Rmd
    │   ├── Fit metalog to each property
    │   ├── Cholesky decomposition
    │   └── Multivariate simulation
    ├── KSSL_soil_property_correlation_matricies.Rmd
    │   ├── Compute empirical correlations
    │   ├── By genhz (A, B, C, O, R, Cr)
    │   └── Using ilr for compositional texture
    └── Generating_Empirical_Parametric_Distributions.Rmd
        └── Create synthetic distributions from stats

COWEETA APPLICATION:
├── coweeta/sim-functions.R (100 KB)
│   ├── Multivariate simulation functions
│   ├── van Genuchten water retention
│   └── AWS calculation
├── coweeta/coweeta_soil_simulation.Rmd
│   ├── SSURGO vs RSS comparison
│   ├── Compositional texture via ilr
│   └── Multi-source integration
├── coweeta/depthwise-multivariate-profile-simulation.Rmd
│   └── Depth-specific correlation structure
├── coweeta/SoilGrids_POLARIS_Download.Rmd
│   └── Data acquisition
└── coweeta/AWS-Map-Simulation.Rmd
    └── Hydrological application

ORCHESTRATION:
├── assemble-pieces.R (8 KB)
│   ├── Downloads from multiple sources
│   ├── Creates comparison raster stacks
│   └── Uses terra, soilDB packages
└── bcqrf_process.Rmd (13 KB)
    └── End-to-end workflow demonstration
```

---

## Key Research Contributions

1. **Systematic Data Infilling** - 6-tier hierarchical approach respecting pedological structure and uncertainty

2. **Metalog-Based Multivariate Simulation** - Flexible bounded distributions with empirical correlations preserved via Cholesky decomposition

3. **Bayesian Constraining Framework** - Integration of empirical priors with model predictions to reduce bias and improve precision

4. **Compositional Data Methods** - Proper treatment of texture sum constraints using isometric log-ratio transformation

5. **Depth Variability Quantification** - Monte Carlo simulation of horizon boundaries using triangular distributions and OSD distinctness

6. **Applied Site Case Study** - Complete Coweeta workflow demonstrating method integration for hydrological applications

---

## Statistics Summary

| Metric | Count/Value |
|--------|-------------|
| **Total Lines of Code** | ~12,800 |
| **R Scripts** | 7 files |
| **R Markdown Notebooks** | 11 files |
| **Python Modules** | 1 file |
| **Functions** | 100+ named functions |
| **Largest File** | soil_infill_functions.R (103 KB) |
| **Data Sources** | SSURGO, KSSL, SoilGrids, POLARIS, RSS |
| **R Packages Used** | 20+ |
| **Python Packages Used** | 3 core packages |

---

## Scientific Applications

### Soil Mapping
- Spatially explicit predictions with uncertainty quantification
- Multi-source data fusion (SSURGO + SoilGrids + POLARIS)
- Gap-filling for incomplete survey areas

### Hydrological Modeling
- Plant-available water storage (AWS) estimation
- Water retention parameters (van Genuchten)
- Drought assessment and water balance modeling

### Uncertainty Quantification
- Probabilistic soil property distributions
- Quantile-based predictions (5%, 50%, 95%)
- Monte Carlo realizations for propagating uncertainty

### Data Fusion
- Integration of multiple soil databases
- Bayesian updating to constrain model predictions
- Comparison of data source impacts (SSURGO vs RSS)

### Environmental Modeling
- Realistic soil inputs for ecosystem models
- Correlated property generation for process-based models
- Depth-resolved soil profiles for vertical process modeling

---

## Computational Complexity

### Data Infilling
- **Time Complexity**: O(n × m²) where n = horizons, m = components per map unit
- **Space Complexity**: O(n × p) where p = properties
- **Parallelization**: Not implemented, but readily parallelizable by map unit

### Metalog Fitting
- **Time Complexity**: O(k³) per property where k = number of quantiles (typically 9)
- **Space Complexity**: O(k²)
- **Note**: Linear regression makes this very fast compared to iterative methods

### Multivariate Simulation
- **Time Complexity**: O(p³ + n×p²) where p = properties, n = realizations
  - O(p³) for Cholesky decomposition (one-time)
  - O(n×p²) for transformation of n samples
- **Space Complexity**: O(n×p)

### Monte Carlo Depth Simulation
- **Time Complexity**: O(r×h) where r = realizations (500), h = horizons
- **Space Complexity**: O(r×h)

---

## Future Development Directions

Based on the codebase structure and comments:

1. **Quantile Regression Forest Integration**
   - Full implementation of constrained QRF (currently conceptual in `constrained_QRF_ideas.Rmd`)
   - Spatial prediction incorporating covariates (elevation, aspect, geology, climate)

2. **Parallelization**
   - Implement parallel processing for map unit-level operations
   - GPU acceleration for large-scale simulations

3. **Web Interface**
   - Shiny application for interactive soil property simulation
   - User-friendly interface for non-programmers

4. **Machine Learning Integration**
   - Deep learning for soil property relationships
   - Automated feature selection for spatial prediction

5. **Expanded Data Sources**
   - Integration with additional global databases (WISE30sec, HWSD v2)
   - Real-time data updates from ongoing surveys

---

## Citation and Usage

This codebase supports reproducible soil science research. When using these methods, consider citing:

- **BCQRF Repository**: https://github.com/jjmaynard/BCQRF
- **Key Dependencies**: `soilDB`, `aqp`, `compositions`, `rmetalog`
- **Data Sources**: SSURGO (NRCS), KSSL, SoilGrids, POLARIS

For questions or contributions, see the repository README and CONTRIBUTING guidelines.

---

**Last Updated**: 2025-11-20
**Author**: Jonathan J Maynard, NMSU
**License**: MIT
