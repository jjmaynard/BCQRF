# Core Soil Property Infilling Functions - Part 1
# Modified to exclude R, Cr, and O horizons from infilling
# Required libraries
library(dplyr)
library(stringr)

# FUNCTION ORGANIZATION
# ---------------------
# 1. Main Functions
#    - downl(oad_ssurgo_tabular() - SSURGO data download, flag restrictions
#    - infill_soil_property() - Main entry point
#
# 2. Horizon Suitability Functions
#    - identify_unsuitable_horizons() - Detects R, Cr, O horizons
#    - is_unsuitable_horizon_type() - Detailed horizon type checking
#    - filter_unsuitable_horizons() - Adds flags and reports
#
# 3. Data Cleaning Functions
#    - clean_property_data() - Cleans and validates data
#    - apply_basic_range_limits() - Applies range validation
#
# 4. Property Configuration Functions
#    - get_default_property_config() - Default configs for common properties
#
# 5. Range Value Infilling Functions
#    - infill_property_range_values() - Main range infilling (clean and focused)
#    - learn_property_ranges() - Learns from existing data
#    - get_property_contextual_ranges() - Pedological knowledge
#    - calculate_property_lower_bound() - Lower bound calculations
#    - calculate_property_upper_bound() - Upper bound calculations
#    - get_contextual_spread() - Property-specific spreads
#
# 6. Property Data Recovery Functions
#    - infill_missing_property_data() - Main recovery strategy
#    - horizon_name_property_infill() - Horizon name matching
#    - standardize_horizon_name() - Name standardization
#    - calculate_horizon_similarity() - Similarity scoring
#    - depth_weighted_property_infill() - Depth-based infilling
#    - within_component_property_interpolation() - Within-component interpolation
#    - cross_component_property_interpolation() - Cross-component interpolation
#    - related_property_estimation() - Property relationship estimation
#    - calculate_depth_weighted_mean() - Statistical calculations
#
# 7. Special Functions
#    - infill_rfv_property() - Special RFV handling
#    - impute_rfv_values() - RFV-specific logic
#    - create_custom_property_config() - Custom configuration helper

# ============================================================================
# 1. MAIN FUNCTION
# ============================================================================

#' Download SSURGO Tabular Data with Horizon Restrictions
#'
#' This function downloads and processes SSURGO tabular data for a specified area of interest (AOI)
#' using its Well-Known Text (WKT) representation and a set of specified soil properties. It also
#' includes horizon restriction data to identify restrictive layers for soil property infilling.
#'
#' @param aoi_wkt Character. A Well-Known Text (WKT) representation of the area of interest.
#' @param properties Character vector. A vector of soil properties to download. Defaults to
#'   \code{c("db", "cec", "rfv", "clay", "ph", "sand", "silt", "soc")}.
#' @param include_restrictions Logical. Whether to include horizon restriction data. Defaults to TRUE.
#'
#' @return A data frame containing SSURGO tabular data with soil property measurements, restriction data,
#'   and associated metadata. When \code{include_restrictions = TRUE}, adds columns for identifying
#'   restrictive layers that should be excluded from soil property infilling.
#'
#' @details The AOI is first converted from WKT to a spatial object and reprojected to EPSG:5070.
#'   Then, the soil map unit keys (\code{mukey}) are retrieved using \code{soilDB::mukey.wcs}. An SQL query is
#'   constructed using the provided soil properties and submitted to the SSURGO database via \code{soilDB::SDA_query}.
#'   If rock fragment volume (\code{rfv}) is one of the requested properties, the function aggregates the low,
#'   representative, and high values per soil horizon.
#'
#'   Restriction data includes: horizon restrictions (reskind, resdept), cementation information,
#'   and other properties useful for identifying R, Cr, O, and other unsuitable horizons for infilling.
#'
#' @import dplyr
#' @import terra
#' @import soilDB
#' @export
download_ssurgo_tabular <- function(aoi_wkt, properties = c("db", "cec", "rfv", "clay", "ph", "sand", "silt", "soc"), include_restrictions = TRUE) {

  # Input validation
  if (missing(aoi_wkt) || is.null(aoi_wkt) || aoi_wkt == "") {
    stop("aoi_wkt cannot be empty or missing")
  }

  if (!is.character(properties) || length(properties) == 0) {
    stop("properties must be a non-empty character vector")
  }

  # Define lookup table for general soil properties and their corresponding SSURGO labels
  ssurgo_lookup <- data.frame(
    Property = c("db", "cec", "rfv", "clay", "ph", "sand", "silt", "soc"),
    SSURGO_Label_Low = c("dbovendry_l", "cec7_l", "chf.fragvol_l AS rfv_l", "claytotal_l", "ph1to1h2o_l", "sandtotal_l", "silttotal_l", "om_l"),
    SSURGO_Label_Rep = c("dbovendry_r", "cec7_r", "chf.fragvol_r AS rfv_r", "claytotal_r", "ph1to1h2o_r", "sandtotal_r", "silttotal_r", "om_r"),
    SSURGO_Label_High = c("dbovendry_h", "cec7_h", "chf.fragvol_h AS rfv_h", "claytotal_h", "ph1to1h2o_h", "sandtotal_h", "silttotal_h", "om_h"),
    stringsAsFactors = FALSE
  )

  # Validate requested properties
  invalid_props <- setdiff(properties, ssurgo_lookup$Property)
  if (length(invalid_props) > 0) {
    warning(paste("Invalid properties ignored:", paste(invalid_props, collapse = ", ")))
    properties <- intersect(properties, ssurgo_lookup$Property)
  }

  if (length(properties) == 0) {
    stop("No valid properties specified")
  }

  # Select the lookup table rows corresponding to the requested properties
  props <- ssurgo_lookup[ssurgo_lookup$Property %in% properties, ]

  # Convert AOI WKT to a spatial object and reproject to EPSG:5070 using terra functions
  tryCatch({
    aoi <- terra::vect(aoi_wkt, crs = "epsg:4326")
    aoi <- terra::project(aoi, "epsg:5070")
    aoi <- terra::as.polygons(terra::ext(aoi), crs = "epsg:5070")
  }, error = function(e) {
    stop(paste("Error processing AOI WKT:", e$message))
  })

  # Fetch SSURGO map unit keys (mukey) using soilDB
  tryCatch({
    mu <- soilDB::mukey.wcs(aoi = aoi, db = "gssurgo")
    mukey_list <- unique(terra::values(mu))
    mukey_list <- mukey_list[!is.na(mukey_list)]
  }, error = function(e) {
    stop(paste("Error fetching mukeys:", e$message))
  })

  if (length(mukey_list) == 0) {
    stop("No map unit keys found for the specified AOI")
  }

  # Format mukey(s) and property(s) for SQL query
  formatted_mukey <- paste0("(", paste0("'", as.integer(mukey_list), "'", collapse = ","), ")")
  ssurgo_variables <- c(props$SSURGO_Label_Low, props$SSURGO_Label_Rep, props$SSURGO_Label_High)
  formatted_properties <- paste(stats::na.omit(ssurgo_variables), collapse = ", ")

  # Add restriction-related fields if requested
  restriction_fields <- ""
  restriction_joins <- ""

  if (include_restrictions) {
    restriction_fields <- ",
      -- restriction data for identifying unsuitable horizons
      chf.fragkind, -- fragment kind (e.g., channers, flagstones)
      chf.fraghard, -- fragment hardness
      -- horizon designation and morphology for restriction detection
      cht.texcl, -- texture class from chtexture table
      cht.lieutex, -- in lieu texture for non-standard materials
      ch.desgnmaster, -- master horizon designation
      ch.desgnmasterprime, -- prime master horizon designation
      ch.desgnvert, -- vertical subdivision
      -- restriction layer information
      cr.reskind, -- restriction kind (e.g., 'Lithic bedrock', 'Paralithic bedrock')
      cr.resdept_l, cr.resdept_r, cr.resdept_h, -- restriction depth
      cr.resthk_l, cr.resthk_r, cr.resthk_h -- restriction thickness"

    restriction_joins <- "
      LEFT JOIN corestrictions AS cr ON co.cokey = cr.cokey
      LEFT JOIN chtexture AS cht ON ch.chkey = cht.chtkey"
  }

  # Construct the SQL query
  query <- sprintf(
    "
    SELECT
      -- contextual data
      co.mukey, co.cokey, compname, comppct_l, comppct_r, comppct_h,
      -- horizon morphology
      ch.chkey, hzname, hzdept_l, hzdept_r, hzdept_h, hzdepb_l, hzdepb_r, hzdepb_h, hzthk_l, hzthk_r, hzthk_h,
      -- soil property parameters
      %s%s
    FROM mapunit AS mu
      INNER JOIN component AS co ON mu.mukey = co.mukey
      INNER JOIN chorizon AS ch ON co.cokey = ch.cokey
      LEFT JOIN chfrags AS chf ON ch.chkey = chf.chkey%s
    WHERE mu.mukey IN %s
    ORDER BY co.cokey, ch.hzdept_r ASC;",
    formatted_properties, restriction_fields, restriction_joins, formatted_mukey
  )

  # Submit the query and fetch results using soilDB
  ssurgo_data <- NULL
  tryCatch({
    ssurgo_data <- soilDB::SDA_query(query)
  }, error = function(e) {
    stop(paste("Error executing SSURGO query:", e$message))
  })

  # Check if query returned valid data
  if (is.null(ssurgo_data)) {
    warning("SSURGO query returned NULL")
    return(data.frame())
  }

  if (!is.data.frame(ssurgo_data)) {
    warning("SSURGO query did not return a data frame")
    return(data.frame())
  }

  if (nrow(ssurgo_data) == 0) {
    warning("No data returned from SSURGO query")
    return(data.frame())
  }

  # If rock fragment volume (rfv) is requested, aggregate its values per soil horizon
  if ("rfv" %in% properties && all(c("rfv_l", "rfv_r", "rfv_h") %in% names(ssurgo_data))) {
    rfv_sum <- ssurgo_data %>%
      dplyr::group_by(chkey) %>%
      dplyr::summarise(
        rfv_l = sum(rfv_l, na.rm = TRUE),
        rfv_r = sum(rfv_r, na.rm = TRUE),
        rfv_h = sum(rfv_h, na.rm = TRUE),
        .groups = 'drop'
      )

    ssurgo_data <- ssurgo_data %>%
      dplyr::select(-dplyr::all_of(c("rfv_l", "rfv_r", "rfv_h"))) %>%
      dplyr::left_join(rfv_sum, by = "chkey") %>%
      dplyr::distinct()
  }

  # Add derived restriction indicators if restriction data was included
  if (include_restrictions && nrow(ssurgo_data) > 0) {
    ssurgo_data <- ssurgo_data %>%
      dplyr::mutate(
        # Flag for bedrock restrictions
        has_lithic_restriction = !is.na(reskind) & grepl("lithic|bedrock", reskind, ignore.case = TRUE),
        has_paralithic_restriction = !is.na(reskind) & grepl("paralithic", reskind, ignore.case = TRUE),

        # Flag for cemented horizons using texture class and in lieu texture
        is_cemented = (!is.na(texcl) & grepl("cemented|indurated", texcl, ignore.case = TRUE)) |
          (!is.na(lieutex) & grepl("cemented|indurated|bedrock|pavement", lieutex, ignore.case = TRUE)),

        # Combine horizon name patterns for restriction detection
        # This helps identify R, Cr, O horizons from hzname
        horizon_suggests_restriction = dplyr::case_when(
          is.na(hzname) ~ FALSE,
          grepl("^R$|^R[0-9]|^Cr|^O", hzname, ignore.case = FALSE) ~ TRUE,
          grepl("qm$|fragipan|duripan|petrocalcic|petrogypsic", hzname, ignore.case = TRUE) ~ TRUE,
          grepl("m$", hzname, ignore.case = FALSE) & nchar(hzname) > 1 ~ TRUE, # cemented horizons with 'm' suffix
          TRUE ~ FALSE
        ),

        # Master designation indicators
        master_designation_restriction = dplyr::case_when(
          is.na(desgnmaster) ~ FALSE,
          desgnmaster %in% c("R", "Cr", "O") ~ TRUE,
          TRUE ~ FALSE
        ),

        # Texture-based restriction indicators
        texture_suggests_restriction = dplyr::case_when(
          (!is.na(texcl) & grepl("cemented|indurated", texcl, ignore.case = TRUE)) ~ TRUE,
          (!is.na(lieutex) & grepl("cemented|indurated|bedrock|pavement", lieutex, ignore.case = TRUE)) ~ TRUE,
          TRUE ~ FALSE
        ),

        # Overall restriction flag for easy filtering
        # This can be used directly in the infilling functions
        is_potentially_restrictive = has_lithic_restriction |
          has_paralithic_restriction |
          is_cemented |
          horizon_suggests_restriction |
          master_designation_restriction |
          texture_suggests_restriction,

        # Add the unsuitable_horizon flag expected by infilling functions
        unsuitable_horizon = is_potentially_restrictive
      )
  } else if (!include_restrictions) {
    # Add basic unsuitable_horizon detection based on horizon names only
    ssurgo_data <- ssurgo_data %>%
      dplyr::mutate(
        unsuitable_horizon = dplyr::case_when(
          is.na(hzname) ~ FALSE,
          grepl("^R$|^R[0-9]|^Cr|^O", hzname, ignore.case = FALSE) ~ TRUE,
          grepl("qm$|fragipan|duripan|petrocalcic|petrogypsic", hzname, ignore.case = TRUE) ~ TRUE,
          grepl("m$", hzname, ignore.case = FALSE) & nchar(hzname) > 1 ~ TRUE,
          TRUE ~ FALSE
        )
      )
  }

  return(ssurgo_data)
}

infill_soil_property <- function(df, property_name, property_config = NULL) {
  #' Flexible function to infill missing data for any specified soil property
  #' Excludes R, Cr, and O horizons from infilling
  #'
  #' @param df Input soil data frame
  #' @param property_name Name of the property to infill (e.g., "sandtotal", "dbovendry")
  #' @param property_config Optional list with property-specific configuration
  #' @return Data frame with infilled property data

  # Validate inputs
  if (!is.character(property_name) || length(property_name) != 1) {
    stop("property_name must be a single character string")
  }

  r_col <- paste0(property_name, "_r")
  if (!r_col %in% names(df)) {
    stop(paste("Column", r_col, "not found in data frame"))
  }

  # Clean and validate the property column
  df <- clean_property_data(df, property_name)

  # Filter out problematic horizons before processing
  df <- filter_unsuitable_horizons(df)

  # Set default configuration if not provided
  if (is.null(property_config)) {
    property_config <- get_default_property_config(property_name)
  }

  # Step 1: Group by 'compname' and process each group (if column exists)
  if ("compname" %in% names(df)) {
    process_group <- function(group) {
      # Process each group individually with horizon-level filtering
      group <- group

      # Identify problematic horizons (missing data in top 50cm if depth info available)
      # BUT EXCLUDE R, Cr, and O horizons from being considered problematic
      if ("hzdepb_r" %in% names(group)) {
        top_50_mask <- group$hzdepb_r <= 50
        missing_data <- is.na(group[[r_col]])
        unsuitable_horizons <- identify_unsuitable_horizons(group)

        # Only consider horizons problematic if they're missing data, in top 50cm,
        # AND are suitable for infilling (not R, Cr, or O)
        problematic_horizons <- top_50_mask & missing_data & !unsuitable_horizons
      } else {
        # If no depth info, treat all missing as problematic except unsuitable horizons
        missing_data <- is.na(group[[r_col]])
        unsuitable_horizons <- identify_unsuitable_horizons(group)
        problematic_horizons <- missing_data & !unsuitable_horizons
      }

      if (any(problematic_horizons)) {
        # Try to infill missing data using various strategies
        group <- infill_missing_property_data(group, property_name, problematic_horizons, property_config)

        # Flag any remaining problematic horizons
        still_missing <- is.na(group[[r_col]])
        unsuitable_horizons <- identify_unsuitable_horizons(group)

        if ("hzdepb_r" %in% names(group)) {
          top_50_mask <- group$hzdepb_r <= 50
          # Data is complete if: not missing OR unsuitable for infilling OR not in top 50cm
          group$property_data_complete <- !(top_50_mask & still_missing & !unsuitable_horizons)
        } else {
          # Data is complete if: not missing OR unsuitable for infilling
          group$property_data_complete <- !(still_missing & !unsuitable_horizons)
        }
      } else {
        group$property_data_complete <- TRUE
      }

      return(group)
    }

    # Apply processing to all groups and combine results
    result_df <- df %>%
      group_by(compname) %>%
      group_modify(~ process_group(.x)) %>%
      ungroup()
  } else {
    # Process entire dataframe as one group
    missing_data <- is.na(df[[r_col]])
    unsuitable_horizons <- identify_unsuitable_horizons(df)
    problematic_horizons <- missing_data & !unsuitable_horizons

    if (any(problematic_horizons)) {
      result_df <- infill_missing_property_data(df, property_name, problematic_horizons, property_config)
      still_missing <- is.na(result_df[[r_col]])
      unsuitable_horizons <- identify_unsuitable_horizons(result_df)
      result_df$property_data_complete <- !(still_missing & !unsuitable_horizons)
    } else {
      result_df <- df
      result_df$property_data_complete <- TRUE
    }
  }

  # Infill the _l and _h values for all valid horizons (excluding unsuitable ones)
  result_df <- infill_property_range_values(result_df, property_name, property_config)

  return(result_df)
}

# ============================================================================
# 2. HORIZON SUITABILITY FUNCTIONS
# ============================================================================

identify_unsuitable_horizons <- function(df) {
  #' Identify horizons that should not be infilled (R, Cr, O horizons)
  #'
  #' @param df Input data frame
  #' @return Logical vector indicating unsuitable horizons

  if (!"hzname" %in% names(df)) {
    # If no horizon names, assume all are suitable
    return(rep(FALSE, nrow(df)))
  }

  unsuitable <- rep(FALSE, nrow(df))

  for (i in 1:nrow(df)) {
    hzname <- df$hzname[i]

    if (is.na(hzname)) {
      next  # Skip missing horizon names
    }

    # Clean and standardize horizon name
    hzname_clean <- toupper(trimws(as.character(hzname)))

    # Remove numbers and special characters for classification
    hzname_base <- str_replace_all(hzname_clean, "[0-9\\-\\.\\s]", "")

    # Check for unsuitable horizon types
    unsuitable[i] <- is_unsuitable_horizon_type(hzname_base)
  }

  return(unsuitable)
}

is_unsuitable_horizon_type <- function(hzname_base) {
  #' Check if a horizon type is unsuitable for standard soil property infilling
  #'
  #' @param hzname_base Cleaned horizon name base
  #' @return Logical indicating if horizon is unsuitable

  # O horizons (organic)
  if (substr(hzname_base, 1, 1) == "O") {
    return(TRUE)
  }

  # R horizons (bedrock)
  if (substr(hzname_base, 1, 1) == "R") {
    return(TRUE)
  }

  # Cr horizons (highly weathered/cemented parent material)
  if (substr(hzname_base, 1, 2) == "CR") {
    return(TRUE)
  }

  # Additional cemented horizons
  cemented_patterns <- c(
    "QM",    # Duripan
    "BQM",   # Cemented B horizon
    "CQM",   # Cemented C horizon
    "BTKQM", # Cemented argillic
    "BKQM",  # Cemented calcic
    "PETROCALCIC",
    "PETROGYPSIC",
    "DURIPAN",
    "FRAGIPAN"
  )

  for (pattern in cemented_patterns) {
    if (grepl(pattern, hzname_base)) {
      return(TRUE)
    }
  }

  # Check for 'm' suffix indicating cementation (e.g., Bkm, Btm)
  if (substr(hzname_base, nchar(hzname_base), nchar(hzname_base)) == "M" &&
      nchar(hzname_base) > 1) {
    return(TRUE)
  }

  return(FALSE)
}

filter_unsuitable_horizons <- function(df) {
  #' Add a flag for unsuitable horizons but don't remove them
  #' This preserves the data structure while marking horizons to skip
  #'
  #' @param df Input data frame
  #' @return Data frame with unsuitable_horizon flag

  df$unsuitable_horizon <- identify_unsuitable_horizons(df)

  # Optional: Print summary of unsuitable horizons found
  if (any(df$unsuitable_horizon) && "hzname" %in% names(df)) {
    unsuitable_types <- unique(df$hzname[df$unsuitable_horizon])
    unsuitable_types <- unsuitable_types[!is.na(unsuitable_types)]

    cat("Found", sum(df$unsuitable_horizon), "unsuitable horizons for infilling:\n")
    cat("Horizon types:", paste(unsuitable_types, collapse = ", "), "\n")
    cat("These horizons will be excluded from infilling but retained in the dataset.\n\n")
  }

  return(df)
}

# ============================================================================
# 3. DATA CLEANING FUNCTIONS
# ============================================================================

clean_property_data <- function(df, property_name) {
  #' Clean and validate property data, converting to numeric where needed
  #'
  #' @param df Input data frame
  #' @param property_name Name of the property to clean
  #' @return Data frame with cleaned property data

  property_cols <- paste0(property_name, c("_r", "_l", "_h"))

  for (col in property_cols) {
    if (col %in% names(df)) {
      # Convert to numeric, handling various input types
      original_vals <- df[[col]]

      # Handle factors by converting to character first
      if (is.factor(original_vals)) {
        df[[col]] <- as.numeric(as.character(original_vals))
      } else if (is.character(original_vals)) {
        # Clean character strings - remove common non-numeric characters
        cleaned_vals <- gsub("[^0-9.-]", "", original_vals)
        df[[col]] <- as.numeric(cleaned_vals)
      } else if (!is.numeric(original_vals)) {
        df[[col]] <- as.numeric(original_vals)
      }

      # Set infinite values to NA
      df[[col]][!is.finite(df[[col]])] <- NA

      # Apply basic range validation for known properties
      if (col == paste0(property_name, "_r")) {
        df[[col]] <- apply_basic_range_limits(df[[col]], property_name)
      }
    }
  }

  return(df)
}

apply_basic_range_limits <- function(values, property_name) {
  #' Apply basic range limits to catch obvious data errors
  #'
  #' @param values Numeric vector of values
  #' @param property_name Name of the property
  #' @return Vector with extreme outliers set to NA

  if (all(is.na(values))) return(values)

  # Define reasonable limits for common properties
  limits <- list(
    sandtotal = c(0, 100),
    claytotal = c(0, 100),
    silttotal = c(0, 100),
    dbovendry = c(0.3, 3.0),
    cec7 = c(0, 200),
    ph1to1h2o = c(2.5, 11.0),
    om = c(0, 100),
    rfv = c(0, 95)
  )

  if (property_name %in% names(limits)) {
    range_limits <- limits[[property_name]]
    values[values < range_limits[1] | values > range_limits[2]] <- NA
  }

  return(values)
}

# ============================================================================
# 4. PROPERTY CONFIGURATION FUNCTIONS
# ============================================================================

get_default_property_config <- function(property_name) {
  #' Get default configuration for common soil properties
  #'
  #' @param property_name Name of the soil property
  #' @return List with default configuration

  # Common texture properties
  if (property_name %in% c('sandtotal', 'claytotal', 'silttotal')) {
    return(list(
      type = 'texture',
      units = '%',
      typical_range = c(0, 100),
      fallback_range = 8,
      sum_constraint = 100,  # For texture properties
      related_properties = c('sandtotal', 'claytotal', 'silttotal')
    ))
  }

  # Bulk density
  else if (property_name == 'dbovendry') {
    return(list(
      type = 'bulk_density',
      units = 'g/cm³',
      typical_range = c(0.8, 2.2),
      fallback_range = 0.15,
      horizon_effects = TRUE
    ))
  }

  # Water retention
  else if (property_name %in% c('wthirdbar', 'wfifteenbar')) {
    return(list(
      type = 'water_retention',
      units = '%',
      typical_range = c(0, 60),
      fallback_range = 3,
      clay_dependent = TRUE,
      related_properties = c('wthirdbar', 'wfifteenbar')
    ))
  }

  # Rock fragment volume
  else if (property_name == 'rfv') {
    return(list(
      type = 'rock_fragments',
      units = '%',
      typical_range = c(0.01, 85),
      fallback_range = 5,
      zero_handling = 'special'  # Special case for zero/NA values
    ))
  }

  # Cation Exchange Capacity
  else if (property_name %in% c('cec7', 'cec82', 'ecec')) {
    return(list(
      type = 'cec',
      units = 'cmol(+)/kg',
      typical_range = c(0, 100),
      fallback_range = 5,
      clay_dependent = TRUE,
      om_dependent = TRUE,
      horizon_effects = TRUE,
      related_properties = c('claytotal', 'om')
    ))
  }

  # pH (various methods)
  else if (property_name %in% c('ph1to1h2o', 'ph01mcacl2', 'phh2o', 'phkcl', 'phnaf')) {
    return(list(
      type = 'ph',
      units = 'pH units',
      typical_range = c(3.0, 10.0),
      fallback_range = 0.5,
      horizon_effects = TRUE,
      parent_material_effects = TRUE,
      method_specific = TRUE
    ))
  }

  # Organic Matter/Carbon
  else if (property_name %in% c('om', 'oc', 'ompc', 'soc', 'ton', 'carbon_org')) {
    return(list(
      type = 'organic_matter',
      units = ifelse(property_name %in% c('om', 'ompc'), '%', 'g/kg'),
      typical_range = c(0, 50),
      fallback_range = 1.5,
      horizon_effects = TRUE,
      depth_dependent = TRUE,
      related_properties = c('cec7', 'ph1to1h2o')
    ))
  }

  # Generic property (fallback)
  else {
    return(list(
      type = 'generic',
      units = 'unknown',
      typical_range = NULL,
      fallback_range = 5,
      horizon_effects = FALSE
    ))
  }
}

# ============================================================================
# 5. RANGE VALUE INFILLING FUNCTIONS
# ============================================================================

infill_property_range_values <- function(df, property_name, property_config) {
  #' Infill _l and _h values for the specified property using data-driven approaches
  #' Excludes unsuitable horizons from range infilling
  #'
  #' @param df Input data frame
  #' @param property_name Name of the property
  #' @param property_config Property configuration list
  #' @return Data frame with infilled range values

  # Learn from existing complete ranges in the dataset (only from suitable horizons)
  learned_ranges <- learn_property_ranges(df, property_name, property_config)

  # Get contextual ranges based on property type and soil science knowledge
  context_ranges <- get_property_contextual_ranges(df, property_name, property_config)

  r_col <- paste0(property_name, "_r")
  l_col <- paste0(property_name, "_l")
  h_col <- paste0(property_name, "_h")

  # Infill missing _l values (only for suitable horizons)
  if (!l_col %in% names(df)) {
    df[[l_col]] <- NA
  }

  # Only infill suitable horizons
  suitable_mask <- if ("unsuitable_horizon" %in% names(df)) !df$unsuitable_horizon else rep(TRUE, nrow(df))
  missing_l_mask <- is.na(df[[l_col]]) & !is.na(df[[r_col]]) & suitable_mask

  if (any(missing_l_mask)) {
    df[missing_l_mask, l_col] <- apply(df[missing_l_mask, ], 1, function(row) {
      calculate_property_lower_bound(row, property_name, learned_ranges, context_ranges, property_config)
    })
  }

  # Infill missing _h values (only for suitable horizons)
  if (!h_col %in% names(df)) {
    df[[h_col]] <- NA
  }

  missing_h_mask <- is.na(df[[h_col]]) & !is.na(df[[r_col]]) & suitable_mask
  if (any(missing_h_mask)) {
    df[missing_h_mask, h_col] <- apply(df[missing_h_mask, ], 1, function(row) {
      calculate_property_upper_bound(row, property_name, learned_ranges, context_ranges, property_config)
    })
  }

  # Apply bounds checking (only to suitable horizons)
  if (!is.null(property_config$typical_range)) {
    df[[l_col]][suitable_mask] <- pmax(df[[l_col]][suitable_mask], property_config$typical_range[1], na.rm = TRUE)
    df[[h_col]][suitable_mask] <- pmin(df[[h_col]][suitable_mask], property_config$typical_range[2], na.rm = TRUE)
  } else {
    df[[l_col]][suitable_mask] <- pmax(df[[l_col]][suitable_mask], 0, na.rm = TRUE)  # Default: no negative values
  }

  # Ensure logical order: _l ≤ _r ≤ _h (with robust numeric conversion, only for suitable horizons)
  # Convert all columns to numeric first
  df[[l_col]] <- as.numeric(as.character(df[[l_col]]))
  df[[r_col]] <- as.numeric(as.character(df[[r_col]]))
  df[[h_col]] <- as.numeric(as.character(df[[h_col]]))

  df[[l_col]][suitable_mask] <- pmin(df[[l_col]][suitable_mask], df[[r_col]][suitable_mask], na.rm = TRUE)
  df[[h_col]][suitable_mask] <- pmax(df[[h_col]][suitable_mask], df[[r_col]][suitable_mask], na.rm = TRUE)

  return(df)
}

learn_property_ranges <- function(df, property_name, property_config) {
  #' Learn typical ranges from existing complete data for any property
  #' Only learns from suitable horizons (excludes R, Cr, O)
  #'
  #' @param df Input data frame
  #' @param property_name Name of the property
  #' @param property_config Property configuration
  #' @return List of learned ranges by context

  r_col <- paste0(property_name, "_r")
  l_col <- paste0(property_name, "_l")
  h_col <- paste0(property_name, "_h")

  # Ensure columns exist and convert to numeric
  for (col in c(r_col, l_col, h_col)) {
    if (col %in% names(df)) {
      df[[col]] <- as.numeric(as.character(df[[col]]))
    }
  }

  # Get complete records (have all three values and are numeric) AND are suitable horizons
  complete_mask <- !is.na(df[[r_col]]) & !is.na(df[[l_col]]) & !is.na(df[[h_col]]) &
    is.finite(df[[r_col]]) & is.finite(df[[l_col]]) & is.finite(df[[h_col]])

  # Exclude unsuitable horizons from learning
  if ("unsuitable_horizon" %in% names(df)) {
    suitable_mask <- !df$unsuitable_horizon
    complete_mask <- complete_mask & suitable_mask
  }

  complete_data <- df[complete_mask, ]

  if (nrow(complete_data) == 0) {
    fallback_range <- ifelse(is.null(property_config$fallback_range), 5, property_config$fallback_range)
    return(list(default_spread = as.numeric(fallback_range)))
  }

  # Calculate actual spreads
  lower_spreads <- complete_data[[r_col]] - complete_data[[l_col]]
  upper_spreads <- complete_data[[h_col]] - complete_data[[r_col]]

  ranges_by_context <- list()

  # Learn ranges by horizon type (if available) - only from suitable horizons
  if ('hzname' %in% names(complete_data)) {
    hznames <- unique(complete_data$hzname[!is.na(complete_data$hzname)])
    for (hzname in hznames) {
      hz_data <- complete_data[complete_data$hzname == hzname & !is.na(complete_data$hzname), ]
      if (nrow(hz_data) >= 3) {  # Need minimum sample size
        hz_lower <- hz_data[[r_col]] - hz_data[[l_col]]
        hz_upper <- hz_data[[h_col]] - hz_data[[r_col]]

        ranges_by_context[[paste0('hzname_', hzname)]] <- list(
          lower_spread_median = median(hz_lower, na.rm = TRUE),
          upper_spread_median = median(hz_upper, na.rm = TRUE),
          lower_spread_std = sd(hz_lower, na.rm = TRUE),
          upper_spread_std = sd(hz_upper, na.rm = TRUE),
          sample_size = nrow(hz_data)
        )
      }
    }
  }

  # Learn ranges by component group (if available)
  if ('compname' %in% names(complete_data)) {
    grps <- unique(complete_data$compname[!is.na(complete_data$compname)])
    for (grp in grps) {
      grp_data <- complete_data[complete_data$compname == grp & !is.na(complete_data$compname), ]
      if (nrow(grp_data) >= 5) {
        grp_lower <- grp_data[[r_col]] - grp_data[[l_col]]
        grp_upper <- grp_data[[h_col]] - grp_data[[r_col]]

        ranges_by_context[[paste0('compgrp_', grp)]] <- list(
          lower_spread_median = median(grp_lower, na.rm = TRUE),
          upper_spread_median = median(grp_upper, na.rm = TRUE),
          lower_spread_std = sd(grp_lower, na.rm = TRUE),
          upper_spread_std = sd(grp_upper, na.rm = TRUE),
          sample_size = nrow(grp_data)
        )
      }
    }
  }

  # Learn ranges by depth zone (if depth info available)
  if ('hzdepb_r' %in% names(complete_data)) {
    depth_zones <- list(
      surface = c(0, 30),
      subsurface = c(30, 100),
      deep = c(100, 200)
    )

    for (zone_name in names(depth_zones)) {
      zone_range <- depth_zones[[zone_name]]
      zone_data <- complete_data[complete_data$hzdepb_r > zone_range[1] &
                                   complete_data$hzdepb_r <= zone_range[2], ]

      if (nrow(zone_data) >= 3) {
        zone_lower <- zone_data[[r_col]] - zone_data[[l_col]]
        zone_upper <- zone_data[[h_col]] - zone_data[[r_col]]

        ranges_by_context[[paste0('depth_', zone_name)]] <- list(
          lower_spread_median = median(zone_lower, na.rm = TRUE),
          upper_spread_median = median(zone_upper, na.rm = TRUE),
          lower_spread_std = sd(zone_lower, na.rm = TRUE),
          upper_spread_std = sd(zone_upper, na.rm = TRUE),
          sample_size = nrow(zone_data)
        )
      }
    }
  }

  # Overall dataset statistics
  ranges_by_context[['overall']] <- list(
    lower_spread_median = median(lower_spreads, na.rm = TRUE),
    upper_spread_median = median(upper_spreads, na.rm = TRUE),
    lower_spread_std = sd(lower_spreads, na.rm = TRUE),
    upper_spread_std = sd(upper_spreads, na.rm = TRUE),
    sample_size = nrow(complete_data)
  )

  return(ranges_by_context)
}

get_property_contextual_ranges <- function(df, property_name, property_config) {
  #' Get property-specific contextual ranges based on soil science knowledge
  #'
  #' @param df Input data frame
  #' @param property_name Name of the property
  #' @param property_config Property configuration
  #' @return List of contextual ranges

  contextual_ranges <- list()

  # Property-specific knowledge
  if (property_config$type == 'texture') {
    contextual_ranges[[property_name]] <- switch(property_name,
                                                 sandtotal = list(typical_spread = 12, min_spread = 5, max_spread = 25),
                                                 claytotal = list(typical_spread = 8, min_spread = 3, max_spread = 20),
                                                 silttotal = list(typical_spread = 10, min_spread = 4, max_spread = 22),
                                                 list(typical_spread = 10, min_spread = 5, max_spread = 20)  # default for other texture
    )
  }

  else if (property_config$type == 'bulk_density') {
    contextual_ranges <- list(
      A_horizons = list(typical_spread = 0.15, range = c(0.8, 1.6)),
      B_horizons = list(typical_spread = 0.12, range = c(1.0, 1.8)),
      C_horizons = list(typical_spread = 0.10, range = c(1.2, 2.0))
    )
  }

  else if (property_config$type == 'water_retention') {
    contextual_ranges[[property_name]] <- switch(property_name,
                                                 wthirdbar = list(typical_spread = 3, clay_factor = 0.4),
                                                 wfifteenbar = list(typical_spread = 2, clay_factor = 0.3),
                                                 list(typical_spread = 2.5, clay_factor = 0.35)  # default
    )
  }

  else if (property_config$type == 'rock_fragments') {
    # RFV ranges depend on the actual value
    contextual_ranges[['rfv_ranges']] <- list(
      low = list(range = c(0.01, 5), typical_spread = 2),
      moderate = list(range = c(5, 15), typical_spread = 4),
      high = list(range = c(15, 35), typical_spread = 6),
      very_high = list(range = c(35, 60), typical_spread = 10),
      extreme = list(range = c(60, 85), typical_spread = 12)
    )
  }

  else {
    # Generic property - use simple percentage-based spread
    contextual_ranges[['generic']] <- list(
      typical_spread = property_config$fallback_range,
      percentage_based = TRUE
    )
  }

  return(contextual_ranges)
}

calculate_property_lower_bound <- function(row, property_name, learned_ranges, context_ranges, property_config) {
  #' Calculate appropriate lower bound for any property
  #'
  #' @param row Data frame row
  #' @param property_name Name of the property
  #' @param learned_ranges Learned ranges from data
  #' @param context_ranges Contextual ranges
  #' @param property_config Property configuration
  #' @return Calculated lower bound

  r_col <- paste0(property_name, "_r")
  r_value <- row[[r_col]]

  # Robust type checking and conversion
  if (is.null(r_value) || length(r_value) == 0) {
    return(NA)
  }

  # Convert to numeric if not already
  if (!is.numeric(r_value)) {
    r_value <- as.numeric(as.character(r_value))
  }

  # Check if conversion resulted in NA or if original was NA
  if (is.na(r_value) || !is.finite(r_value)) {
    return(NA)
  }

  # Priority order for range estimation
  spread_estimates <- list()

  # 1. Try horizon-specific learned range
  if ('hzname' %in% names(row) && !is.na(row[['hzname']])) {
    hz_key <- paste0("hzname_", row[['hzname']])
    if (hz_key %in% names(learned_ranges) && learned_ranges[[hz_key]]$sample_size >= 3) {
      spread <- learned_ranges[[hz_key]]$lower_spread_median
      if (!is.na(spread) && spread > 0) {
        spread_estimates <- append(spread_estimates, list(list(method = 'horizon_learned', spread = spread)))
      }
    }
  }

  # 2. Try component group learned range
  if ('compname' %in% names(row) && !is.na(row[['compname']])) {
    grp_key <- paste0("compgrp_", row[['compname']])
    if (grp_key %in% names(learned_ranges) && learned_ranges[[grp_key]]$sample_size >= 5) {
      spread <- learned_ranges[[grp_key]]$lower_spread_median
      if (!is.na(spread) && spread > 0) {
        spread_estimates <- append(spread_estimates, list(list(method = 'group_learned', spread = spread)))
      }
    }
  }

  # 3. Try depth zone learned range
  if ('hzdepb_r' %in% names(row) && !is.na(row[['hzdepb_r']])) {
    depth <- row[['hzdepb_r']]
    zone_key <- if (depth <= 30) 'depth_surface' else if (depth <= 100) 'depth_subsurface' else 'depth_deep'

    if (zone_key %in% names(learned_ranges) && learned_ranges[[zone_key]]$sample_size >= 3) {
      spread <- learned_ranges[[zone_key]]$lower_spread_median
      if (!is.na(spread) && spread > 0) {
        spread_estimates <- append(spread_estimates, list(list(method = 'depth_learned', spread = spread)))
      }
    }
  }

  # 4. Use overall dataset learned range
  if ('overall' %in% names(learned_ranges)) {
    spread <- learned_ranges[['overall']]$lower_spread_median
    if (!is.na(spread) && spread > 0) {
      spread_estimates <- append(spread_estimates, list(list(method = 'overall_learned', spread = spread)))
    }
  }

  # 5. Use contextual/pedological knowledge
  ctx_spread <- get_contextual_spread(row, property_name, context_ranges, property_config, 'lower')
  if (!is.na(ctx_spread)) {
    spread_estimates <- append(spread_estimates, list(list(method = 'contextual', spread = ctx_spread)))
  }

  # 6. Fallback to default
  if (length(spread_estimates) == 0) {
    fallback <- ifelse(is.null(property_config$fallback_range), 5, property_config$fallback_range)
    spread_estimates <- append(spread_estimates, list(list(method = 'fallback', spread = fallback)))
  }

  # Use the first (highest priority) estimate
  spread <- spread_estimates[[1]]$spread

  # Ensure spread is numeric and finite
  if (!is.numeric(spread) || !is.finite(spread) || spread < 0) {
    spread <- ifelse(is.null(property_config$fallback_range), 5, property_config$fallback_range)
    spread <- as.numeric(spread)
  }

  # Calculate lower bound
  lower_bound <- r_value - spread

  # Apply property-specific constraints
  if (!is.null(property_config$typical_range)) {
    lower_bound <- max(lower_bound, property_config$typical_range[1])
  } else {
    lower_bound <- max(lower_bound, 0)  # Default: never below 0
  }

  return(lower_bound)
}

calculate_property_upper_bound <- function(row, property_name, learned_ranges, context_ranges, property_config) {
  #' Calculate appropriate upper bound for any property
  #'
  #' @param row Data frame row
  #' @param property_name Name of the property
  #' @param learned_ranges Learned ranges from data
  #' @param context_ranges Contextual ranges
  #' @param property_config Property configuration
  #' @return Calculated upper bound

  r_col <- paste0(property_name, "_r")
  r_value <- row[[r_col]]

  # Robust type checking and conversion
  if (is.null(r_value) || length(r_value) == 0) {
    return(NA)
  }

  # Convert to numeric if not already
  if (!is.numeric(r_value)) {
    r_value <- as.numeric(as.character(r_value))
  }

  # Check if conversion resulted in NA or if original was NA
  if (is.na(r_value) || !is.finite(r_value)) {
    return(NA)
  }

  # Similar logic to lower bound but for upper spreads
  spread_estimates <- list()

  # 1. Horizon-specific learned range
  if ('hzname' %in% names(row) && !is.na(row[['hzname']])) {
    hz_key <- paste0("hzname_", row[['hzname']])
    if (hz_key %in% names(learned_ranges) && learned_ranges[[hz_key]]$sample_size >= 3) {
      spread <- learned_ranges[[hz_key]]$upper_spread_median
      if (!is.na(spread) && spread > 0) {
        spread_estimates <- append(spread_estimates, list(list(method = 'horizon_learned', spread = spread)))
      }
    }
  }

  # 2. Component group learned range
  if ('compname' %in% names(row) && !is.na(row[['compname']])) {
    grp_key <- paste0("compgrp_", row[['compname']])
    if (grp_key %in% names(learned_ranges) && learned_ranges[[grp_key]]$sample_size >= 5) {
      spread <- learned_ranges[[grp_key]]$upper_spread_median
      if (!is.na(spread) && spread > 0) {
        spread_estimates <- append(spread_estimates, list(list(method = 'group_learned', spread = spread)))
      }
    }
  }

  # 3. Depth zone learned range
  if ('hzdepb_r' %in% names(row) && !is.na(row[['hzdepb_r']])) {
    depth <- row[['hzdepb_r']]
    zone_key <- if (depth <= 30) 'depth_surface' else if (depth <= 100) 'depth_subsurface' else 'depth_deep'

    if (zone_key %in% names(learned_ranges) && learned_ranges[[zone_key]]$sample_size >= 3) {
      spread <- learned_ranges[[zone_key]]$upper_spread_median
      if (!is.na(spread) && spread > 0) {
        spread_estimates <- append(spread_estimates, list(list(method = 'depth_learned', spread = spread)))
      }
    }
  }

  # 4. Overall dataset learned range
  if ('overall' %in% names(learned_ranges)) {
    spread <- learned_ranges[['overall']]$upper_spread_median
    if (!is.na(spread) && spread > 0) {
      spread_estimates <- append(spread_estimates, list(list(method = 'overall_learned', spread = spread)))
    }
  }

  # 5. Contextual knowledge
  ctx_spread <- get_contextual_spread(row, property_name, context_ranges, property_config, 'upper')
  if (!is.na(ctx_spread)) {
    spread_estimates <- append(spread_estimates, list(list(method = 'contextual', spread = ctx_spread)))
  }

  # 6. Fallback
  if (length(spread_estimates) == 0) {
    fallback <- ifelse(is.null(property_config$fallback_range), 5, property_config$fallback_range)
    spread_estimates <- append(spread_estimates, list(list(method = 'fallback', spread = fallback)))
  }

  spread <- spread_estimates[[1]]$spread

  # Ensure spread is numeric and finite
  if (!is.numeric(spread) || !is.finite(spread) || spread < 0) {
    spread <- ifelse(is.null(property_config$fallback_range), 5, property_config$fallback_range)
    spread <- as.numeric(spread)
  }

  upper_bound <- r_value + spread

  # Apply property-specific upper limits
  if (!is.null(property_config$typical_range)) {
    upper_bound <- min(upper_bound, property_config$typical_range[2])
  }

  return(upper_bound)
}

get_contextual_spread <- function(row, property_name, context_ranges, property_config, bound_type) {
  #' Get contextual spread based on property type and soil characteristics
  #'
  #' @param row Data frame row
  #' @param property_name Name of the property
  #' @param context_ranges Contextual ranges
  #' @param property_config Property configuration
  #' @param bound_type 'lower' or 'upper'
  #' @return Contextual spread value

  r_col <- paste0(property_name, "_r")
  r_value <- row[[r_col]]

  # Convert to numeric if needed
  if (!is.numeric(r_value)) {
    r_value <- as.numeric(as.character(r_value))
  }

  # Return fallback if r_value is not valid
  if (is.na(r_value) || !is.finite(r_value)) {
    fallback <- ifelse(is.null(property_config$fallback_range), 5, property_config$fallback_range)
    return(as.numeric(fallback))
  }

  # Texture properties
  if (property_config$type == 'texture' && property_name %in% names(context_ranges)) {
    return(context_ranges[[property_name]]$typical_spread)
  }

  # Bulk density - varies by horizon type
  else if (property_config$type == 'bulk_density') {
    if ('hzname' %in% names(row) && !is.na(row[['hzname']])) {
      hzname <- toupper(as.character(row[['hzname']]))
      if (substr(hzname, 1, 1) == 'A' && 'A_horizons' %in% names(context_ranges)) {
        return(context_ranges[['A_horizons']]$typical_spread)
      } else if (substr(hzname, 1, 1) == 'B' && 'B_horizons' %in% names(context_ranges)) {
        return(context_ranges[['B_horizons']]$typical_spread)
      } else if (substr(hzname, 1, 1) == 'C' && 'C_horizons' %in% names(context_ranges)) {
        return(context_ranges[['C_horizons']]$typical_spread)
      }
    }
    return(property_config$fallback_range)
  }

  # Water retention - adjust for clay content
  else if (property_config$type == 'water_retention' && property_name %in% names(context_ranges)) {
    base_spread <- context_ranges[[property_name]]$typical_spread

    if ('claytotal_r' %in% names(row) && !is.na(row[['claytotal_r']])) {
      clay_content <- row[['claytotal_r']]
      clay_factor <- context_ranges[[property_name]]$clay_factor
      adjusted_spread <- base_spread + (clay_content * clay_factor / 100)
      return(adjusted_spread)
    }

    return(base_spread)
  }

  # CEC - varies by clay content and organic matter
  else if (property_config$type == 'cec') {
    base_spread <- property_config$fallback_range

    # Adjust based on clay content (higher clay = higher CEC variability)
    if ('claytotal_r' %in% names(row) && !is.na(row[['claytotal_r']])) {
      clay_content <- row[['claytotal_r']]
      clay_adjustment <- clay_content * 0.1  # 0.1 cmol/kg per % clay
      base_spread <- base_spread + clay_adjustment
    }

    # Adjust based on organic matter (higher OM = higher CEC variability)
    if ('om_r' %in% names(row) && !is.na(row[['om_r']])) {
      om_content <- row[['om_r']]
      om_adjustment <- om_content * 2  # 2 cmol/kg per % OM
      base_spread <- base_spread + om_adjustment
    }

    # Horizon-specific adjustments
    if ('hzname' %in% names(row) && !is.na(row[['hzname']])) {
      hzname <- toupper(as.character(row[['hzname']]))
      if (substr(hzname, 1, 1) == 'A') {
        base_spread <- base_spread * 1.3  # Surface horizons more variable
      } else if (substr(hzname, 1, 1) == 'C') {
        base_spread <- base_spread * 0.8  # Parent material more uniform
      }
    }

    return(base_spread)
  }

  # pH - varies by horizon and parent material
  else if (property_config$type == 'ph') {
    base_spread <- property_config$fallback_range

    # Horizon-specific adjustments
    if ('hzname' %in% names(row) && !is.na(row[['hzname']])) {
      hzname <- toupper(as.character(row[['hzname']]))
      if (substr(hzname, 1, 1) == 'A') {
        base_spread <- base_spread * 1.4  # Surface horizons more variable due to management
      } else if (grepl('BT', hzname)) {
        base_spread <- base_spread * 0.8  # Bt horizons often more uniform
      }
    }

    return(base_spread)
  }

  # Organic matter - strongly depth dependent
  else if (property_config$type == 'organic_matter') {
    base_spread <- property_config$fallback_range

    # Depth-dependent adjustments (surface horizons much more variable)
    if ('hzdept_r' %in% names(row) && !is.na(row[['hzdept_r']])) {
      depth <- row[['hzdept_r']]
      if (depth <= 15) {
        base_spread <- base_spread * 2.0  # Surface very variable
      } else if (depth <= 30) {
        base_spread <- base_spread * 1.5  # Subsurface moderately variable
      } else {
        base_spread <- base_spread * 0.6  # Deep horizons more uniform
      }
    }

    # Horizon-specific adjustments
    if ('hzname' %in% names(row) && !is.na(row[['hzname']])) {
      hzname <- toupper(as.character(row[['hzname']]))
      if (substr(hzname, 1, 1) == 'A') {
        base_spread <- base_spread * 1.5  # Surface mineral horizons highly variable
      }
      # Note: O horizons excluded from infilling, so no case needed
    }

    return(base_spread)
  }

  # Rock fragments - varies by value range
  else if (property_config$type == 'rock_fragments' && 'rfv_ranges' %in% names(context_ranges)) {
    rfv_ranges <- context_ranges[['rfv_ranges']]

    if (r_value <= 5) {
      return(rfv_ranges[['low']]$typical_spread)
    } else if (r_value <= 15) {
      return(rfv_ranges[['moderate']]$typical_spread)
    } else if (r_value <= 35) {
      return(rfv_ranges[['high']]$typical_spread)
    } else if (r_value <= 60) {
      return(rfv_ranges[['very_high']]$typical_spread)
    } else {
      return(rfv_ranges[['extreme']]$typical_spread)
    }
  }

  # Generic property - percentage-based or fixed spread
  else if (property_config$type == 'generic') {
    if ('generic' %in% names(context_ranges) && context_ranges[['generic']]$percentage_based) {
      return(r_value * 0.1)  # 10% of the value as default spread
    } else {
      return(property_config$fallback_range)
    }
  }

  # Fallback
  fallback <- ifelse(is.null(property_config$fallback_range), 5, property_config$fallback_range)
  return(as.numeric(fallback))
}

# ============================================================================
# 6. PROPERTY DATA RECOVERY FUNCTIONS
# ============================================================================

infill_missing_property_data <- function(group, property_name, problematic_mask, property_config) {
  #' Attempt to infill missing property data using various strategies
  #' Only infills suitable horizons (excludes R, Cr, O)
  #'
  #' @param group Data frame group to process
  #' @param property_name Name of the property to infill
  #' @param problematic_mask Logical vector for problematic horizons
  #' @param property_config Property configuration
  #' @return Data frame with infilled property data

  group <- group
  property_col <- paste0(property_name, "_r")

  # Ensure we have depth information for depth-based strategies
  has_depth <- all(c('hzdept_r', 'hzdepb_r') %in% names(group))

  # Strategy 1: PRIORITY - Match by horizon name (pedologically similar horizons)
  # Only among suitable horizons
  if ('hzname' %in% names(group)) {
    group <- horizon_name_property_infill(group, property_col, problematic_mask)
  }

  # Check if we still have missing data after horizon name matching
  still_missing <- is.na(group[[property_col]]) & problematic_mask
  if (!any(still_missing)) {
    return(group)  # All data recovered, no need for further strategies
  }

  # Strategy 2: Depth-weighted averaging from similar depth ranges (if depth available)
  # Only from suitable horizons
  if (has_depth) {
    group <- depth_weighted_property_infill(group, property_col, problematic_mask)
  }

  # Check progress after depth-weighted infill
  still_missing <- is.na(group[[property_col]]) & problematic_mask
  if (!any(still_missing)) {
    return(group)
  }

  # Strategy 3: Vertical interpolation within individual soil components (if depth available)
  # Only among suitable horizons
  if (has_depth) {
    group <- within_component_property_interpolation(group, property_col)
  }

  # Check progress after interpolation
  still_missing <- is.na(group[[property_col]]) & problematic_mask
  if (!any(still_missing)) {
    return(group)
  }

  # Strategy 4: Cross-component interpolation for similar depths (if depth available)
  # Only from suitable horizons
  if (has_depth) {
    group <- cross_component_property_interpolation(group, property_col)
  }

  # Strategy 5: Use related properties for estimation (property-specific)
  # Only for suitable horizons
  group <- related_property_estimation(group, property_name, property_config)

  # Strategy 6: Fallback to group statistics for any remaining missing values
  # Only using suitable horizons for statistics and only infilling suitable horizons
  remaining_missing <- is.na(group[[property_col]]) & problematic_mask
  if (any(remaining_missing)) {
    # Calculate statistics only from suitable horizons
    suitable_mask <- if ("unsuitable_horizon" %in% names(group)) !group$unsuitable_horizon else rep(TRUE, nrow(group))
    suitable_group <- group[suitable_mask, ]

    # Use depth-weighted group mean as last resort (if depth available)
    if (has_depth && nrow(suitable_group) > 0) {
      depth_weighted_mean <- calculate_depth_weighted_mean(suitable_group, property_col)
    } else if (nrow(suitable_group) > 0) {
      depth_weighted_mean <- mean(suitable_group[[property_col]], na.rm = TRUE)
    } else {
      depth_weighted_mean <- NA
    }

    if (!is.na(depth_weighted_mean)) {
      group[[property_col]][remaining_missing] <- depth_weighted_mean
    }
  }

  return(group)
}

horizon_name_property_infill <- function(group, property_col, problematic_mask) {
  #' Infill missing property data using horizons with matching or similar names
  #' Only uses suitable horizons as sources and only infills suitable horizons
  #'
  #' @param group Data frame group to process
  #' @param property_col Name of the property column to infill
  #' @param problematic_mask Logical vector for problematic horizons
  #' @return Data frame with infilled values

  # Input validation
  if (!is.data.frame(group) || nrow(group) == 0) {
    return(group)
  }

  if (!property_col %in% names(group)) {
    warning(paste("Property column", property_col, "not found in data"))
    return(group)
  }

  if (!"hzname" %in% names(group)) {
    warning("hzname column not found - cannot perform horizon name matching")
    return(group)
  }

  if (length(problematic_mask) != nrow(group)) {
    warning("problematic_mask length does not match number of rows")
    return(group)
  }

  # Initialize tracking column if not present
  if (!"infill_method" %in% names(group)) {
    group$infill_method <- ""
  }

  # Determine suitable horizons
  suitable_mask <- if ("unsuitable_horizon" %in% names(group)) {
    !group$unsuitable_horizon
  } else {
    rep(TRUE, nrow(group))
  }

  # Get indices that need infilling
  problematic_indices <- which(problematic_mask &
                                 is.na(group[[property_col]]) &
                                 !is.na(group$hzname))

  if (length(problematic_indices) == 0) {
    return(group)
  }

  # Pre-process horizon names for all rows to improve efficiency
  standardized_hznames <- sapply(group$hzname, function(x) {
    if (is.na(x)) return(NA_character_)
    standardize_horizon_name(x)
  })

  # Find potential source horizons (suitable, non-missing data, non-missing hzname)
  source_mask <- suitable_mask &
    !is.na(group[[property_col]]) &
    !is.na(standardized_hznames)

  source_indices <- which(source_mask)

  if (length(source_indices) == 0) {
    warning("No suitable source horizons found for horizon name matching")
    return(group)
  }

  # Process each problematic horizon
  for (idx in problematic_indices) {
    target_hzname_clean <- standardized_hznames[idx]

    if (is.na(target_hzname_clean) || target_hzname_clean == "") {
      next
    }

    # Find matching horizons from suitable sources
    matching_values <- numeric(0)
    similarity_scores <- numeric(0)

    for (source_idx in source_indices) {
      if (source_idx == idx) next  # Skip self

      source_hzname_clean <- standardized_hznames[source_idx]

      if (is.na(source_hzname_clean) || source_hzname_clean == "") {
        next
      }

      # Calculate similarity
      similarity <- calculate_horizon_similarity(target_hzname_clean, source_hzname_clean)

      if (similarity > 0.5) {  # Only use reasonably similar horizons
        matching_values <- c(matching_values, group[[property_col]][source_idx])
        similarity_scores <- c(similarity_scores, similarity)
      }
    }

    # Infill if matches found
    if (length(matching_values) > 0) {
      # Calculate weighted average
      if (length(matching_values) == 1) {
        infilled_value <- matching_values[1]
      } else {
        # Check for valid weights
        if (all(is.finite(similarity_scores)) && all(similarity_scores > 0)) {
          infilled_value <- weighted.mean(matching_values, similarity_scores)
        } else {
          infilled_value <- mean(matching_values)  # Fallback to simple mean
        }
      }

      # Assign infilled value
      group[[property_col]][idx] <- infilled_value

      # Update tracking information
      method_info <- paste0(property_col, ":hzname(", target_hzname_clean, ",n=",
                            length(matching_values), "); ")

      if (is.na(group$infill_method[idx]) || group$infill_method[idx] == "") {
        group$infill_method[idx] <- method_info
      } else {
        group$infill_method[idx] <- paste0(group$infill_method[idx], method_info)
      }
    }
  }

  return(group)
}

standardize_horizon_name <- function(hzname) {
  #' Standardize horizon names for better matching
  #'
  #' @param hzname Raw horizon name
  #' @return Cleaned and standardized horizon name

  if (is.na(hzname)) {
    return("")
  }

  # Convert to string and clean
  hzname <- toupper(trimws(as.character(hzname)))

  # Remove numbers at the end (layer designations)
  hzname <- str_replace(hzname, "\\d+$", "")

  # Remove common punctuation but keep important ones like '/'
  hzname <- str_replace_all(hzname, "[^\\w/]", "")

  return(hzname)
}

calculate_horizon_similarity <- function(hz1, hz2) {
  #' Calculate similarity between two horizon names
  #'
  #' @param hz1,hz2 Horizon names to compare
  #' @return Similarity score from 0 (no similarity) to 1 (identical)

  if (hz1 == hz2) {
    return(1.0)
  }

  if (hz1 == "" || hz2 == "") {
    return(0.0)
  }

  # Check for main horizon letter match
  main_hz1 <- ifelse(nchar(hz1) > 0, substr(hz1, 1, 1), "")
  main_hz2 <- ifelse(nchar(hz2) > 0, substr(hz2, 1, 1), "")

  if (main_hz1 == main_hz2) {
    base_score <- 0.8  # Same main horizon type

    # Bonus for additional character matches
    chars1 <- strsplit(hz1, "")[[1]]
    chars2 <- strsplit(hz2, "")[[1]]
    common_chars <- intersect(chars1, chars2)
    bonus <- length(common_chars) / max(length(chars1), length(chars2)) * 0.2

    return(min(1.0, base_score + bonus))
  }

  # Check for related horizons (pedologically similar)
  related_groups <- list(
    c('A', 'AP', 'AE'),  # Surface horizons
    c('E', 'EB', 'BE'),  # Eluvial horizons
    c('B', 'BT', 'BW', 'BC', 'BS'),  # Illuvial/subsurface horizons
    c('C', 'CB', 'CR'),  # Parent material - NOTE: CR included for similarity but will be excluded from infilling
    c('O', 'OA', 'OE')   # Organic horizons - NOTE: O included for similarity but will be excluded from infilling
  )

  for (group in related_groups) {
    if (main_hz1 %in% group && main_hz2 %in% group) {
      return(0.6)  # Related horizon types
    }
  }

  # Handle transitional horizons like 'B/C', 'A/E'
  if (grepl("/", hz1) || grepl("/", hz2)) {
    parts1 <- ifelse(grepl("/", hz1), strsplit(hz1, "/")[[1]], hz1)
    parts2 <- ifelse(grepl("/", hz2), strsplit(hz2, "/")[[1]], hz2)

    max_sim <- 0
    for (p1 in parts1) {
      for (p2 in parts2) {
        sim <- calculate_horizon_similarity(trimws(p1), trimws(p2))
        max_sim <- max(max_sim, sim)
      }
    }

    return(max_sim * 0.8)  # Slight penalty for composite matching
  }

  return(0.0)  # No similarity
}

depth_weighted_property_infill <- function(group, property_col, problematic_mask) {
  #' Fill missing values using depth-weighted averages from similar depth ranges
  #' Only uses suitable horizons as sources
  #'
  #' @param group Data frame group to process
  #' @param property_col Name of the property column to infill
  #' @param problematic_mask Logical vector for problematic horizons
  #' @return Data frame with infilled values

  # Input validation
  if (!is.data.frame(group) || nrow(group) == 0) {
    return(group)
  }

  if (!property_col %in% names(group)) {
    warning(paste("Property column", property_col, "not found in data"))
    return(group)
  }

  if (length(problematic_mask) != nrow(group)) {
    warning("problematic_mask length does not match number of rows")
    return(group)
  }

  # Check for required depth columns
  required_depth_cols <- c("hzdept_r", "hzdepb_r")
  missing_depth_cols <- setdiff(required_depth_cols, names(group))

  if (length(missing_depth_cols) > 0) {
    warning(paste("Missing required depth columns:", paste(missing_depth_cols, collapse = ", ")))
    return(group)
  }

  # Initialize tracking column if not present
  if (!"infill_method" %in% names(group)) {
    group$infill_method <- ""
  }

  # Determine suitable horizons
  suitable_mask <- if ("unsuitable_horizon" %in% names(group)) {
    !group$unsuitable_horizon
  } else {
    rep(TRUE, nrow(group))
  }

  # Get indices that need infilling
  problematic_indices <- which(problematic_mask &
                                 is.na(group[[property_col]]) &
                                 !is.na(group$hzdept_r) &
                                 !is.na(group$hzdepb_r))

  if (length(problematic_indices) == 0) {
    return(group)
  }

  # Pre-identify potential source horizons
  source_mask <- suitable_mask &
    !is.na(group[[property_col]]) &
    !is.na(group$hzdept_r) &
    !is.na(group$hzdepb_r) &
    is.finite(group[[property_col]]) &
    is.finite(group$hzdept_r) &
    is.finite(group$hzdepb_r)

  source_indices <- which(source_mask)

  if (length(source_indices) == 0) {
    warning("No suitable source horizons found for depth-weighted infilling")
    return(group)
  }

  # Pre-calculate source horizon properties for efficiency
  source_tops <- group$hzdept_r[source_indices]
  source_bottoms <- group$hzdepb_r[source_indices]
  source_mids <- (source_tops + source_bottoms) / 2
  source_values <- group[[property_col]][source_indices]

  # Process each problematic horizon
  for (idx in problematic_indices) {
    target_top <- group$hzdept_r[idx]
    target_bottom <- group$hzdepb_r[idx]

    # Validate target depth values
    if (!is.finite(target_top) || !is.finite(target_bottom) || target_top >= target_bottom) {
      next
    }

    target_mid <- (target_top + target_bottom) / 2

    # Calculate overlaps and similarities for all source horizons
    overlaps <- pmax(0, pmin(target_bottom, source_bottoms) - pmax(target_top, source_tops))
    depth_diffs <- abs(target_mid - source_mids)
    depth_similarities <- 1 / (1 + depth_diffs)

    # Apply tolerance filter (20cm tolerance for depth difference)
    tolerance_mask <- (overlaps > 0) | (depth_diffs < 20)

    if (!any(tolerance_mask)) {
      next  # No horizons within tolerance
    }

    # Filter to horizons within tolerance, excluding self
    valid_source_mask <- tolerance_mask & (source_indices != idx)

    if (!any(valid_source_mask)) {
      next  # No valid sources after filtering
    }

    # Get filtered values and calculate weights
    filtered_overlaps <- overlaps[valid_source_mask]
    filtered_similarities <- depth_similarities[valid_source_mask]
    filtered_values <- source_values[valid_source_mask]

    # Calculate combined weights (overlap + depth similarity)
    weights <- filtered_overlaps + filtered_similarities

    # Ensure weights are positive and finite
    valid_weight_mask <- is.finite(weights) & weights > 0

    if (!any(valid_weight_mask)) {
      next  # No valid weights
    }

    # Final filtering
    final_weights <- weights[valid_weight_mask]
    final_values <- filtered_values[valid_weight_mask]

    # Calculate weighted average
    if (length(final_values) > 0) {
      if (length(final_values) == 1) {
        weighted_avg <- final_values[1]
      } else {
        weighted_avg <- weighted.mean(final_values, final_weights)
      }

      # Validate result
      if (is.finite(weighted_avg)) {
        group[[property_col]][idx] <- weighted_avg

        # Update tracking information
        depth_range <- paste0(round(target_top, 1), "-", round(target_bottom, 1), "cm")
        method_info <- paste0(property_col, ":depth_wt(", depth_range, ",n=",
                              length(final_values), "); ")

        if (is.na(group$infill_method[idx]) || group$infill_method[idx] == "") {
          group$infill_method[idx] <- method_info
        } else {
          group$infill_method[idx] <- paste0(group$infill_method[idx], method_info)
        }
      }
    }
  }

  return(group)
}

within_component_property_interpolation <- function(group, property_col) {
  #' Interpolate missing values within individual soil components based on depth
  #' Only uses suitable horizons for interpolation
  #'
  #' @param group Data frame group to process
  #' @param property_col Name of the property column to infill
  #' @return Data frame with interpolated values

  group <- group
  suitable_mask <- if ("unsuitable_horizon" %in% names(group)) !group$unsuitable_horizon else rep(TRUE, nrow(group))

  # Group by individual soil component if component ID exists
  if ('cokey' %in% names(group)) {
    component_groups <- split(group, group$cokey)
  } else {
    # Fallback: treat entire group as one component
    component_groups <- list(group)
  }

  processed_components <- list()

  for (comp_data in component_groups) {
    comp_data <- comp_data[order(comp_data$hzdepb_r), ]
    comp_suitable_mask <- if ("unsuitable_horizon" %in% names(comp_data)) !comp_data$unsuitable_horizon else rep(TRUE, nrow(comp_data))

    # Only interpolate among suitable horizons
    suitable_comp_data <- comp_data[comp_suitable_mask, ]

    if (any(is.na(suitable_comp_data[[property_col]])) && any(!is.na(suitable_comp_data[[property_col]]))) {
      # Interpolate missing values based on depth using approx
      valid_indices <- which(!is.na(suitable_comp_data[[property_col]]))
      if (length(valid_indices) >= 2) {
        # Create interpolation for suitable horizons only
        missing_suitable_indices <- which(is.na(suitable_comp_data[[property_col]]))

        interpolated_values <- approx(x = suitable_comp_data$hzdepb_r[valid_indices],
                                      y = suitable_comp_data[[property_col]][valid_indices],
                                      xout = suitable_comp_data$hzdepb_r[missing_suitable_indices],
                                      rule = 2)$y

        # Map back to original component data
        suitable_comp_data[[property_col]][missing_suitable_indices] <- interpolated_values
        comp_data[comp_suitable_mask, ] <- suitable_comp_data
      }
    }

    processed_components <- append(processed_components, list(comp_data))
  }

  if (length(processed_components) > 1) {
    return(do.call(rbind, processed_components))
  } else {
    return(processed_components[[1]])
  }
}

cross_component_property_interpolation <- function(group, property_col) {
  #' Use data from other components at similar depths to fill missing values
  #' Only uses suitable horizons as sources
  #'
  #' @param group Data frame group to process
  #' @param property_col Name of the property column to infill
  #' @return Data frame with infilled values

  # Input validation
  if (!is.data.frame(group) || nrow(group) == 0) {
    return(group)
  }

  if (!property_col %in% names(group)) {
    warning(paste("Property column", property_col, "not found in data"))
    return(group)
  }

  # Check for required depth columns
  required_depth_cols <- c("hzdept_r", "hzdepb_r")
  missing_depth_cols <- setdiff(required_depth_cols, names(group))

  if (length(missing_depth_cols) > 0) {
    warning(paste("Missing required depth columns:", paste(missing_depth_cols, collapse = ", ")))
    return(group)
  }

  # Initialize tracking column if not present
  if (!"infill_method" %in% names(group)) {
    group$infill_method <- ""
  }

  # Determine suitable horizons
  suitable_mask <- if ("unsuitable_horizon" %in% names(group)) {
    !group$unsuitable_horizon
  } else {
    rep(TRUE, nrow(group))
  }

  # Only try to infill suitable horizons with missing data and valid depths
  missing_mask <- is.na(group[[property_col]]) &
    suitable_mask &
    !is.na(group$hzdept_r) &
    !is.na(group$hzdepb_r) &
    is.finite(group$hzdept_r) &
    is.finite(group$hzdepb_r)

  if (!any(missing_mask)) {
    return(group)
  }

  # Pre-identify potential source horizons
  source_mask <- suitable_mask &
    !is.na(group[[property_col]]) &
    !is.na(group$hzdept_r) &
    !is.na(group$hzdepb_r) &
    is.finite(group[[property_col]]) &
    is.finite(group$hzdept_r) &
    is.finite(group$hzdepb_r)

  source_indices <- which(source_mask)

  if (length(source_indices) == 0) {
    warning("No suitable source horizons found for cross-component interpolation")
    return(group)
  }

  # Pre-calculate source horizon properties for efficiency
  source_tops <- group$hzdept_r[source_indices]
  source_bottoms <- group$hzdepb_r[source_indices]
  source_mids <- (source_tops + source_bottoms) / 2
  source_values <- group[[property_col]][source_indices]

  # Validate source midpoints
  valid_source_mask <- is.finite(source_mids)
  if (!any(valid_source_mask)) {
    warning("No valid source horizon midpoints found")
    return(group)
  }

  # Filter sources to only valid ones
  source_indices <- source_indices[valid_source_mask]
  source_mids <- source_mids[valid_source_mask]
  source_values <- source_values[valid_source_mask]

  missing_indices <- which(missing_mask)
  depth_tolerance <- 15  # 15cm tolerance for similar depths

  # Process each missing value
  for (idx in missing_indices) {
    target_top <- group$hzdept_r[idx]
    target_bottom <- group$hzdepb_r[idx]

    # Validate target depths
    if (!is.finite(target_top) || !is.finite(target_bottom) || target_top >= target_bottom) {
      next
    }

    target_depth_mid <- (target_top + target_bottom) / 2

    if (!is.finite(target_depth_mid)) {
      next
    }

    # Calculate depth differences for all source horizons
    depth_diffs <- abs(target_depth_mid - source_mids)

    # Find horizons within tolerance, excluding self
    within_tolerance <- (depth_diffs <= depth_tolerance) & (source_indices != idx)

    if (!any(within_tolerance)) {
      next  # No horizons within tolerance
    }

    # Get matching values and calculate weights
    matching_values <- source_values[within_tolerance]
    matching_diffs <- depth_diffs[within_tolerance]

    # Calculate weights (closer depths get higher weight)
    weights <- 1 / (1 + matching_diffs)

    # Validate weights and values
    valid_match_mask <- is.finite(weights) &
      is.finite(matching_values) &
      weights > 0

    if (!any(valid_match_mask)) {
      next  # No valid matches
    }

    # Filter to valid matches
    final_values <- matching_values[valid_match_mask]
    final_weights <- weights[valid_match_mask]

    # Calculate weighted average
    if (length(final_values) > 0) {
      if (length(final_values) == 1) {
        weighted_value <- final_values[1]
      } else {
        weighted_value <- weighted.mean(final_values, final_weights)
      }

      # Validate result
      if (is.finite(weighted_value)) {
        group[[property_col]][idx] <- weighted_value

        # Update tracking information
        depth_info <- paste0(round(target_depth_mid, 1), "cm")
        method_info <- paste0(property_col, ":cross_comp(", depth_info, ",n=",
                              length(final_values), ",tol=", depth_tolerance, "cm); ")

        if (is.na(group$infill_method[idx]) || group$infill_method[idx] == "") {
          group$infill_method[idx] <- method_info
        } else {
          group$infill_method[idx] <- paste0(group$infill_method[idx], method_info)
        }
      }
    }
  }

  return(group)
}

related_property_estimation <- function(group, property_name, property_config) {
  #' Use related properties to estimate missing values (property-specific logic)
  #' Only estimates for suitable horizons
  #'
  #' @param group Data frame group to process
  #' @param property_name Name of the property to estimate
  #' @param property_config Property configuration
  #' @return Data frame with estimated values

  group <- group
  property_col <- paste0(property_name, "_r")
  suitable_mask <- if ("unsuitable_horizon" %in% names(group)) !group$unsuitable_horizon else rep(TRUE, nrow(group))

  # Only proceed if we have related properties defined
  if (is.null(property_config$related_properties)) {
    return(group)
  }

  missing_mask <- is.na(group[[property_col]]) & suitable_mask  # Only estimate for suitable horizons
  if (!any(missing_mask)) return(group)

  # Texture properties - use sum constraint (sand + clay + silt ≈ 100%)
  if (property_config$type == 'texture') {
    texture_cols <- paste0(property_config$related_properties, "_r")
    available_cols <- intersect(texture_cols, names(group))

    if (length(available_cols) >= 2) {  # Need at least 2 other texture components
      for (idx in which(missing_mask)) {
        other_values <- group[idx, available_cols]
        if (sum(!is.na(other_values)) >= 2) {
          sum_others <- sum(other_values, na.rm = TRUE)
          estimated_value <- max(0, min(100, 100 - sum_others))  # Constrain to 0-100%
          group[[property_col]][idx] <- estimated_value
        }
      }
    }
  }

  # Water retention - use clay relationship if available
  else if (property_config$type == 'water_retention' && 'claytotal_r' %in% names(group)) {
    # Simple empirical relationships (these could be made more sophisticated)
    clay_col <- 'claytotal_r'

    for (idx in which(missing_mask)) {
      clay_content <- group[[clay_col]][idx]
      if (!is.na(clay_content)) {
        if (property_name == 'wthirdbar') {
          # Rough approximation: field capacity ≈ 0.3 * clay + 10
          estimated_value <- max(0, min(60, 0.3 * clay_content + 10))
        } else if (property_name == 'wfifteenbar') {
          # Rough approximation: wilting point ≈ 0.4 * clay + 2
          estimated_value <- max(0, min(40, 0.4 * clay_content + 2))
        } else {
          next  # Unknown water retention property
        }
        group[[property_col]][idx] <- estimated_value
      }
    }
  }

  # CEC - use clay and organic matter relationships
  else if (property_config$type == 'cec') {
    for (idx in which(missing_mask)) {
      clay <- ifelse('claytotal_r' %in% names(group), group$claytotal_r[idx], NA)
      om <- ifelse('om_r' %in% names(group), group$om_r[idx], NA)

      if (!is.na(clay) || !is.na(om)) {
        # Base CEC estimation using typical relationships
        estimated_cec <- 0

        if (!is.na(clay)) {
          # Clay contribution: varies by clay type (using moderate activity assumption)
          estimated_cec <- estimated_cec + (clay * 0.5)  # 0.5 cmol/kg per % clay
        }

        if (!is.na(om)) {
          # Organic matter contribution: ~200 cmol/kg per % OM
          estimated_cec <- estimated_cec + (om * 20)
        }

        # Minimum CEC even for sandy soils
        estimated_cec <- max(2, estimated_cec)

        # Constrain to reasonable range
        estimated_value <- max(0, min(100, estimated_cec))
        group[[property_col]][idx] <- estimated_value
      }
    }
  }

  # pH - use landscape position and parent material if available
  else if (property_config$type == 'ph') {
    for (idx in which(missing_mask)) {
      # Default to slightly acidic (common in many agricultural soils)
      estimated_ph <- 6.2

      # Adjust based on horizon if available
      if ('hzname' %in% names(group) && !is.na(group$hzname[idx])) {
        hzname <- toupper(as.character(group$hzname[idx]))
        if (substr(hzname, 1, 1) == 'A') {
          estimated_ph <- estimated_ph - 0.3  # Surface horizons often more acidic
        } else if (substr(hzname, 1, 1) == 'C') {
          estimated_ph <- estimated_ph + 0.2  # Parent material often less acidic
        }
      }

      # Adjust based on organic matter if available (high OM often more acidic)
      if ('om_r' %in% names(group) && !is.na(group$om_r[idx])) {
        om_content <- group$om_r[idx]
        if (om_content > 5) {
          estimated_ph <- estimated_ph - 0.4  # High OM soils often more acidic
        }
      }

      # Constrain to reasonable range
      estimated_value <- max(3.0, min(10.0, estimated_ph))
      group[[property_col]][idx] <- estimated_value
    }
  }

  # Organic matter - use depth and horizon relationships
  else if (property_config$type == 'organic_matter') {
    for (idx in which(missing_mask)) {
      # Base estimation on depth if available
      estimated_om <- 2.0  # Default moderate OM

      if ('hzdept_r' %in% names(group) && !is.na(group$hzdept_r[idx])) {
        depth <- group$hzdept_r[idx]
        if (depth <= 15) {
          estimated_om <- 3.5  # Surface horizon
        } else if (depth <= 30) {
          estimated_om <- 1.8  # Shallow subsurface
        } else if (depth <= 50) {
          estimated_om <- 0.8  # Deep subsurface
        } else {
          estimated_om <- 0.3  # Very deep
        }
      }

      # Adjust based on horizon type if available
      if ('hzname' %in% names(group) && !is.na(group$hzname[idx])) {
        hzname <- toupper(as.character(group$hzname[idx]))
        if (substr(hzname, 1, 1) == 'A') {
          estimated_om <- estimated_om * 1.5  # Surface mineral horizon
        } else if (substr(hzname, 1, 1) == 'C') {
          estimated_om <- 0.2  # Parent material
        }
        # Note: O horizons are excluded from infilling, so no case for them here
      }

      # Adjust based on texture if available (clay holds more OM)
      if ('claytotal_r' %in% names(group) && !is.na(group$claytotal_r[idx])) {
        clay_content <- group$claytotal_r[idx]
        if (clay_content > 35) {
          estimated_om <- estimated_om * 1.3  # High clay soils retain more OM
        } else if (clay_content < 15) {
          estimated_om <- estimated_om * 0.7  # Sandy soils retain less OM
        }
      }

      # Constrain to reasonable range
      estimated_value <- max(0, min(50, estimated_om))
      group[[property_col]][idx] <- estimated_value
    }
  }

  # Bulk density - use texture relationship if available
  else if (property_config$type == 'bulk_density' && any(c('sandtotal_r', 'claytotal_r') %in% names(group))) {
    for (idx in which(missing_mask)) {
      # Simple relationship: higher clay = lower bulk density, higher sand = higher bulk density
      clay <- ifelse('claytotal_r' %in% names(group), group$claytotal_r[idx], NA)
      sand <- ifelse('sandtotal_r' %in% names(group), group$sandtotal_r[idx], NA)

      if (!is.na(clay) || !is.na(sand)) {
        # Base bulk density around 1.4 g/cm³, adjust for texture
        base_bd <- 1.4

        if (!is.na(clay)) {
          base_bd <- base_bd - (clay - 20) * 0.01  # Decrease with higher clay
        }
        if (!is.na(sand)) {
          base_bd <- base_bd + (sand - 50) * 0.005  # Increase with higher sand
        }

        # Constrain to reasonable range
        estimated_value <- max(0.8, min(2.2, base_bd))
        group[[property_col]][idx] <- estimated_value
      }
    }
  }

  return(group)
}

calculate_depth_weighted_mean <- function(group, property_col) {
  #' Calculate depth-weighted mean for a property, giving more weight to thicker horizons
  #' Only uses suitable horizons for calculation
  #'
  #' @param group Data frame group
  #' @param property_col Name of the property column
  #' @return Depth-weighted mean value

  suitable_mask <- if ("unsuitable_horizon" %in% names(group)) !group$unsuitable_horizon else rep(TRUE, nrow(group))
  valid_data <- group[!is.na(group[[property_col]]) & suitable_mask, ]

  if (nrow(valid_data) == 0) {
    return(NA)
  }

  # Calculate horizon thickness as weight (if depth info available)
  if (all(c('hzdept_r', 'hzdepb_r') %in% names(valid_data))) {
    thicknesses <- valid_data$hzdepb_r - valid_data$hzdept_r
    values <- valid_data[[property_col]]

    if (sum(thicknesses) == 0) {
      return(mean(values))
    }

    return(weighted.mean(values, thicknesses))
  } else {
    # Fallback to simple mean if no depth info
    return(mean(valid_data[[property_col]]))
  }
}

# ============================================================================
# 7. SPECIAL FUNCTIONS
# ============================================================================

infill_rfv_property <- function(df) {
  #' Special function for RFV infilling that handles zero/NA cases
  #' Only processes suitable horizons
  #'
  #' @param df Input data frame
  #' @return Data frame with infilled RFV values

  # Apply row-wise RFV imputation only to suitable horizons
  suitable_mask <- if ("unsuitable_horizon" %in% names(df)) !df$unsuitable_horizon else rep(TRUE, nrow(df))

  result_df <- df %>%
    rowwise() %>%
    do({
      current_row_num <- as.numeric(rownames(.))
      if (length(current_row_num) == 0 || current_row_num > length(suitable_mask)) {
        # Handle edge case where row number can't be determined
        if ("unsuitable_horizon" %in% names(.) && .$unsuitable_horizon) {
          .  # Return unchanged for unsuitable horizons
        } else {
          impute_rfv_values(.)
        }
      } else if (suitable_mask[current_row_num]) {
        impute_rfv_values(.)
      } else {
        .  # Return unchanged for unsuitable horizons
      }
    }) %>%
    ungroup()

  return(result_df)
}

impute_rfv_values <- function(row) {
  #' Improved RFV (Rock Fragment Volume) imputation with better logic and safety checks
  #' Only for suitable horizons
  #'
  #' @param row Data frame row (as named vector or list)
  #' @return Modified row with imputed RFV values

  # Convert row to list for easier manipulation
  if (is.data.frame(row)) {
    row <- as.list(row[1, ])
  }

  # Handle missing RFV or zero RFV - both mean "no significant rock fragments"
  if (is.na(row[["rfv_r"]]) || row[["rfv_r"]] == 0) {
    row[["rfv_r"]] <- 0.02
    row[["rfv_l"]] <- 0.01
    row[["rfv_h"]] <- 0.03
    return(as.data.frame(row))
  }

  rfv_r <- row[["rfv_r"]]

  # Ensure rfv_r is within reasonable bounds (0.01-85%)
  rfv_r <- max(0.01, min(rfv_r, 85))
  row[["rfv_r"]] <- rfv_r

  # Determine appropriate range based on rfv_r value and context
  if (rfv_r >= 0.01 && rfv_r <= 5) {
    lower_spread <- min(rfv_r - 0.01, 2)
    upper_spread <- 4
  } else if (rfv_r > 5 && rfv_r <= 15) {
    lower_spread <- 3
    upper_spread <- 5
  } else if (rfv_r > 15 && rfv_r <= 35) {
    lower_spread <- 5
    upper_spread <- 8
  } else if (rfv_r > 35 && rfv_r <= 60) {
    lower_spread <- 8
    upper_spread <- 12
  } else {  # rfv_r > 60
    lower_spread <- 10
    upper_spread <- min(15, 85 - rfv_r)
  }

  # Adjust spreads based on horizon context if available
  if ('hzname' %in% names(row) && !is.na(row[['hzname']])) {
    hzname <- toupper(as.character(row[['hzname']]))

    if (substr(hzname, 1, 1) == 'A' && rfv_r > 0.01) {
      upper_spread <- upper_spread * 1.3  # Surface horizons more variable
    } else if (grepl('BT', hzname)) {
      lower_spread <- lower_spread * 0.8  # Bt horizons often more uniform
      upper_spread <- upper_spread * 0.8
    } else if (substr(hzname, 1, 1) == 'C') {
      lower_spread <- lower_spread * 1.2  # Parent material more variable
      upper_spread <- upper_spread * 1.2
    }
  }

  # Calculate bounds with safety checks
  rfv_l <- max(0.01, rfv_r - lower_spread)
  rfv_h <- min(85, rfv_r + upper_spread)

  # Only update _l and _h if they're missing
  if (is.na(row[["rfv_l"]])) {
    row[["rfv_l"]] <- rfv_l
  }

  if (is.na(row[["rfv_h"]])) {
    row[["rfv_h"]] <- rfv_h
  }

  # Final consistency check
  if (!is.na(row[["rfv_l"]]) && !is.na(row[["rfv_h"]])) {
    row[["rfv_l"]] <- min(row[["rfv_l"]], rfv_r)
    row[["rfv_h"]] <- max(row[["rfv_h"]], rfv_r)

    if (row[["rfv_l"]] <= 0) {
      row[["rfv_l"]] <- 0.01
    }

    if (row[["rfv_h"]] - row[["rfv_l"]] < 0.02) {
      spread <- max(0.02, rfv_r * 0.2)
      row[["rfv_l"]] <- max(0.01, rfv_r - spread/2)
      row[["rfv_h"]] <- min(85, rfv_r + spread/2)
    }
  }

  return(as.data.frame(row))
}

create_custom_property_config <- function(property_name,
                                          property_type = "generic",
                                          units = "unknown",
                                          typical_range = NULL,
                                          fallback_range = 5,
                                          related_properties = NULL,
                                          special_options = NULL) {
  #' Create a custom property configuration
  #'
  #' @param property_name Name of the property
  #' @param property_type Type: "texture", "bulk_density", "water_retention", "rock_fragments", "cec", "ph", "organic_matter", "generic"
  #' @param units Units of measurement
  #' @param typical_range Numeric vector c(min, max) of typical values
  #' @param fallback_range Default spread for range calculation
  #' @param related_properties Character vector of related property names
  #' @param special_options List of additional options
  #' @return Property configuration list

  # Validate inputs
  if (!is.character(property_name) || length(property_name) != 1) {
    stop("property_name must be a single character string")
  }

  valid_types <- c("texture", "bulk_density", "water_retention", "rock_fragments",
                   "cec", "ph", "organic_matter", "generic")
  if (!property_type %in% valid_types) {
    warning(paste("property_type should be one of:", paste(valid_types, collapse = ", ")))
  }

  # Start with base configuration
  config <- list(
    type = property_type,
    units = units,
    fallback_range = fallback_range
  )

  # Add typical range if provided
  if (!is.null(typical_range)) {
    if (!is.numeric(typical_range) || length(typical_range) != 2) {
      stop("typical_range must be a numeric vector of length 2: c(min, max)")
    }
    if (typical_range[1] >= typical_range[2]) {
      stop("typical_range[1] must be less than typical_range[2]")
    }
    config$typical_range <- typical_range
  }

  # Add related properties if provided
  if (!is.null(related_properties)) {
    if (!is.character(related_properties)) {
      stop("related_properties must be a character vector")
    }
    config$related_properties <- related_properties
  }

  # Add type-specific defaults based on property_type
  if (property_type == "texture") {
    config$sum_constraint <- 100
    if (is.null(config$typical_range)) {
      config$typical_range <- c(0, 100)
    }
    if (is.null(config$related_properties)) {
      config$related_properties <- c('sandtotal', 'claytotal', 'silttotal')
    }
  }

  else if (property_type == "bulk_density") {
    config$horizon_effects <- TRUE
    if (is.null(config$typical_range)) {
      config$typical_range <- c(0.8, 2.2)
    }
  }

  else if (property_type == "water_retention") {
    config$clay_dependent <- TRUE
    if (is.null(config$typical_range)) {
      config$typical_range <- c(0, 60)
    }
    if (is.null(config$related_properties)) {
      config$related_properties <- c('wthirdbar', 'wfifteenbar')
    }
  }

  else if (property_type == "rock_fragments") {
    config$zero_handling <- 'special'
    if (is.null(config$typical_range)) {
      config$typical_range <- c(0.01, 85)
    }
  }

  else if (property_type == "cec") {
    config$clay_dependent <- TRUE
    config$om_dependent <- TRUE
    config$horizon_effects <- TRUE
    if (is.null(config$typical_range)) {
      config$typical_range <- c(0, 100)
    }
    if (is.null(config$related_properties)) {
      config$related_properties <- c('claytotal', 'om')
    }
  }

  else if (property_type == "ph") {
    config$horizon_effects <- TRUE
    config$parent_material_effects <- TRUE
    config$method_specific <- TRUE
    if (is.null(config$typical_range)) {
      config$typical_range <- c(3.0, 10.0)
    }
  }

  else if (property_type == "organic_matter") {
    config$horizon_effects <- TRUE
    config$depth_dependent <- TRUE
    if (is.null(config$typical_range)) {
      config$typical_range <- c(0, 50)
    }
    if (is.null(config$related_properties)) {
      config$related_properties <- c('cec7', 'ph1to1h2o')
    }
  }

  # Add any special options provided
  if (!is.null(special_options)) {
    if (!is.list(special_options)) {
      stop("special_options must be a list")
    }
    config <- c(config, special_options)
  }

  return(config)
}

# ============================================================================
# HELPER FUNCTIONS FOR SPECIAL OPERATIONS
# ============================================================================

validate_property_config <- function(config, property_name) {
  #' Validate a property configuration object
  #'
  #' @param config Property configuration list
  #' @param property_name Name of the property (for error messages)
  #' @return TRUE if valid, stops with error if invalid

  if (!is.list(config)) {
    stop(paste("Configuration for", property_name, "must be a list"))
  }

  required_fields <- c("type", "units", "fallback_range")
  missing_fields <- setdiff(required_fields, names(config))
  if (length(missing_fields) > 0) {
    stop(paste("Configuration for", property_name, "missing required fields:",
               paste(missing_fields, collapse = ", ")))
  }

  if (!is.numeric(config$fallback_range) || config$fallback_range <= 0) {
    stop(paste("fallback_range for", property_name, "must be a positive number"))
  }

  if (!is.null(config$typical_range)) {
    if (!is.numeric(config$typical_range) || length(config$typical_range) != 2) {
      stop(paste("typical_range for", property_name, "must be a numeric vector of length 2"))
    }
    if (config$typical_range[1] >= config$typical_range[2]) {
      stop(paste("typical_range for", property_name, "must have min < max"))
    }
  }

  return(TRUE)
}

get_rfv_range_category <- function(rfv_value) {
  #' Categorize RFV values into standard ranges
  #'
  #' @param rfv_value Rock fragment volume percentage
  #' @return Character string indicating category

  if (is.na(rfv_value) || rfv_value <= 0) {
    return("none")
  } else if (rfv_value <= 5) {
    return("low")
  } else if (rfv_value <= 15) {
    return("moderate")
  } else if (rfv_value <= 35) {
    return("high")
  } else if (rfv_value <= 60) {
    return("very_high")
  } else {
    return("extreme")
  }
}

apply_property_constraints <- function(values, property_config) {
  #' Apply property-specific constraints to values
  #'
  #' @param values Numeric vector of values to constrain
  #' @param property_config Property configuration
  #' @return Constrained values

  if (is.null(values) || all(is.na(values))) {
    return(values)
  }

  # Apply typical range constraints if available
  if (!is.null(property_config$typical_range)) {
    values <- pmax(values, property_config$typical_range[1], na.rm = TRUE)
    values <- pmin(values, property_config$typical_range[2], na.rm = TRUE)
  }

  # Apply specific constraints based on property type
  if (property_config$type == "texture") {
    values <- pmax(values, 0, na.rm = TRUE)
    values <- pmin(values, 100, na.rm = TRUE)
  } else if (property_config$type == "ph") {
    values <- pmax(values, 0, na.rm = TRUE)
    values <- pmin(values, 14, na.rm = TRUE)
  } else if (property_config$type == "rock_fragments") {
    values <- pmax(values, 0.01, na.rm = TRUE)
    values <- pmin(values, 95, na.rm = TRUE)
  } else {
    # Generic constraint: no negative values unless explicitly allowed
    values <- pmax(values, 0, na.rm = TRUE)
  }

  return(values)
}


# ============================================================================
# SYSTEM DOCUMENTATION AND TECHNICAL NOTES
# ============================================================================

# FUNCTION ORGANIZATION
# ---------------------
# 1. Main Function
#    - infill_soil_property() - Main entry point
#
# 2. Horizon Suitability Functions
#    - identify_unsuitable_horizons() - Detects R, Cr, O horizons
#    - is_unsuitable_horizon_type() - Detailed horizon type checking
#    - filter_unsuitable_horizons() - Adds flags and reports
#
# 3. Data Cleaning Functions
#    - clean_property_data() - Cleans and validates data
#    - apply_basic_range_limits() - Applies range validation
#
# 4. Property Configuration Functions
#    - get_default_property_config() - Default configs for common properties
#
# 5. Range Value Infilling Functions
#    - infill_property_range_values() - Main range infilling (clean and focused)
#    - learn_property_ranges() - Learns from existing data
#    - get_property_contextual_ranges() - Pedological knowledge
#    - calculate_property_lower_bound() - Lower bound calculations
#    - calculate_property_upper_bound() - Upper bound calculations
#    - get_contextual_spread() - Property-specific spreads
#
# 6. Property Data Recovery Functions
#    - infill_missing_property_data() - Main recovery strategy
#    - horizon_name_property_infill() - Horizon name matching
#    - standardize_horizon_name() - Name standardization
#    - calculate_horizon_similarity() - Similarity scoring
#    - depth_weighted_property_infill() - Depth-based infilling
#    - within_component_property_interpolation() - Within-component interpolation
#    - cross_component_property_interpolation() - Cross-component interpolation
#    - related_property_estimation() - Property relationship estimation
#    - calculate_depth_weighted_mean() - Statistical calculations
#
# 7. Special Functions
#    - infill_rfv_property() - Special RFV handling
#    - impute_rfv_values() - RFV-specific logic
#    - create_custom_property_config() - Custom configuration helper

# HORIZON EXCLUSION SYSTEM
# -------------------------
# Automatically excludes these horizon types from infilling:
# - R horizons: Bedrock
# - Cr horizons: Highly weathered/cemented parent material
# - O horizons: Organic horizons
# - Cemented horizons: QM, BQM, CQM, fragipan, duripan, etc.
# - 'm' suffix: Bkm, Btm (cemented horizons)

# RECOVERY STRATEGY HIERARCHY
# ---------------------------
# 1. Horizon name matching (Priority 1): Pedologically similar horizons
# 2. Depth-weighted averaging (Priority 2): Similar depth ranges with overlap
# 3. Within-component interpolation (Priority 3): Vertical trends in profiles
# 4. Cross-component interpolation (Priority 4): Similar depths, different components
# 5. Related property estimation (Priority 5): Soil science relationships
# 6. Statistical fallback (Priority 6): Group means as last resort

# RANGE CALCULATION APPROACH
# ---------------------------
# Hierarchical approach for _l and _h values:
# 1. Learned ranges from existing complete data (by horizon, component group, depth zone)
# 2. Contextual/pedological knowledge (property-specific)
# 3. Fallback values (conservative defaults)
# Always maintains logical order: _l ≤ _r ≤ _h

# PROPERTY-SPECIFIC FEATURES
# ---------------------------
# Texture: Sum constraint (sand + clay + silt ≈ 100%), cross-property estimation
# Bulk Density: Horizon-dependent (A < B < C), texture relationships
# CEC: Clay + organic matter dependencies, horizon effects
# pH: Horizon and parent material effects, depth relationships
# Water Retention: Clay-dependent empirical relationships
# RFV: Special zero/NA handling, value-dependent range calculations
# Organic Matter: Depth-dependent, horizon-specific adjustments

# DATA QUALITY FEATURES
# ----------------------
# - Robust type conversion (factor → character → numeric)
# - Range validation and outlier detection
# - Consistency checking (_l ≤ _r ≤ _h)
# - Method tracking for transparency
# - Comprehensive error handling
# - Validation functions for configurations

# EXTENSIBILITY
# -------------
# - Custom property configurations via create_custom_property_config()
# - Support for new property types
# - Flexible special options system
# - Modular function organization
# - Easy integration with existing workflows

# PERFORMANCE CONSIDERATIONS
# ---------------------------
# - Progressive strategy application with early termination
# - Efficient grouping and vectorization where possible
# - Memory-conscious processing for large datasets
# - Optional progress reporting for long-running operations




# ============================================================================
# COMPREHENSIVE USAGE EXAMPLES - SOIL PROPERTY INFILLING SYSTEM
# ============================================================================

# Complete examples covering all 4 parts of the soil infilling functions
# All functions automatically exclude R, Cr, and O horizons

# Basic Usage Examples
# --------------------

# # 1. BASIC SINGLE PROPERTY INFILLING
# # Basic usage - automatically excludes R, Cr, O horizons
# infilled_data <- infill_soil_property(soil_data, "sandtotal")
#
# # Check which horizons were excluded from infilling
# excluded_horizons <- infilled_data[infilled_data$unsuitable_horizon == TRUE, c("hzname", "cokey")]
# print(excluded_horizons)
#
# # 2. MULTIPLE PROPERTY PROCESSING
# # Process multiple properties - all will skip unsuitable horizons
# properties <- c("sandtotal", "claytotal", "dbovendry", "cec7", "ph1to1h2o")
# for (prop in properties) {
#   infilled_data <- infill_soil_property(infilled_data, prop)
# }
#
# # Check unsuitable horizon summary
# excluded_summary <- table(infilled_data$hzname[infilled_data$unsuitable_horizon == TRUE])
# print(excluded_summary)
#
# # Advanced Range Infilling Examples
# # ---------------------------------
#
# # 3. RANGE VALUE ANALYSIS
# # Basic range infilling for sand content
# sand_config <- get_default_property_config("sandtotal")
# soil_data_with_ranges <- infill_property_range_values(soil_data, "sandtotal", sand_config)
#
# # Check learned ranges from dataset
# learned_sand_ranges <- learn_property_ranges(soil_data, "sandtotal", sand_config)
# print(learned_sand_ranges$overall)  # Overall dataset statistics
#
# # Get contextual ranges for bulk density
# bd_config <- get_default_property_config("dbovendry")
# bd_context <- get_property_contextual_ranges(soil_data, "dbovendry", bd_config)
# print(bd_context)
#
# # Property Recovery Strategy Examples
# # -----------------------------------
#
# # 4. RECOVERY STRATEGY TESTING
# # Main recovery strategy for sand content
# problematic_mask <- is.na(soil_data$sandtotal_r) & !soil_data$unsuitable_horizon
# recovered_data <- infill_missing_property_data(soil_data, "sandtotal",
#                                                problematic_mask, sand_config)
#
# # Test horizon name matching
# horizon_matched <- horizon_name_property_infill(soil_data, "sandtotal_r", problematic_mask)
#
# # Check horizon similarity scores
# hz1 <- "Bt1"
# hz2 <- "Bt2"
# similarity <- calculate_horizon_similarity(standardize_horizon_name(hz1),
#                                            standardize_horizon_name(hz2))
# print(paste("Similarity between", hz1, "and", hz2, ":", similarity))
#
# # Special Functions Examples
# # --------------------------
#
# # 5. RFV SPECIAL HANDLING
# # Process RFV with special zero/NA handling
# soil_data_rfv <- infill_rfv_property(soil_data)
#
# # Check RFV categories
# soil_data_rfv$rfv_category <- sapply(soil_data_rfv$rfv_r, get_rfv_range_category)
# table(soil_data_rfv$rfv_category)
#
# # 6. CUSTOM PROPERTY CONFIGURATIONS
# # Create config for custom nitrogen property
# nitrogen_config <- create_custom_property_config(
#   property_name = "nitrogen_total",
#   property_type = "generic",
#   units = "mg/kg",
#   typical_range = c(0, 50),
#   fallback_range = 2,
#   special_options = list(
#     measurement_method = "kjeldahl",
#     detection_limit = 0.1
#   )
# )
#
# # Validate the configuration
# validate_property_config(nitrogen_config, "nitrogen_total")
#
# # Custom chemical property
# phosphorus_config <- create_custom_property_config(
#   property_name = "phosphorus_available",
#   property_type = "generic",
#   units = "ppm",
#   typical_range = c(0, 200),
#   fallback_range = 5,
#   related_properties = c("ph1to1h2o", "claytotal"),
#   special_options = list(
#     extraction_method = "bray1",
#     horizon_dependent = TRUE
#   )
# )
#
# # Texture property with sum constraint
# custom_sand_config <- create_custom_property_config(
#   property_name = "sand_fine",
#   property_type = "texture",
#   units = "%",
#   typical_range = c(0, 80),
#   fallback_range = 6,
#   related_properties = c("sand_coarse", "sand_medium", "sand_fine"),
#   special_options = list(
#     particle_size_range = "0.1-0.25mm",
#     sum_constraint = 100
#   )
# )
#
# # Complete Workflow Examples
# # --------------------------
#
# # 7. COMPREHENSIVE SOIL DATASET PROCESSING
# # Complete workflow for processing a soil dataset
# process_complete_soil_dataset <- function(soil_data) {
#
#   # Standard physical properties
#   physical_properties <- c("sandtotal", "claytotal", "silttotal", "dbovendry", "rfv")
#
#   # Chemical properties
#   chemical_properties <- c("cec7", "ph1to1h2o", "om")
#
#   # Water retention properties
#   water_properties <- c("wthirdbar", "wfifteenbar")
#
#   # Process all properties
#   all_properties <- c(physical_properties, chemical_properties, water_properties)
#
#   result_data <- soil_data
#   processing_log <- list()
#
#   for (prop in all_properties) {
#     cat("Processing", prop, "...\n")
#
#     # Special handling for RFV
#     if (prop == "rfv") {
#       result_data <- infill_rfv_property(result_data)
#     } else {
#       # Standard processing
#       if (paste0(prop, "_r") %in% names(result_data)) {
#         result_data <- infill_soil_property(result_data, prop)
#
#         # Track processing results
#         missing_before <- sum(is.na(soil_data[[paste0(prop, "_r")]]))
#         missing_after <- sum(is.na(result_data[[paste0(prop, "_r")]]))
#         processing_log[[prop]] <- list(
#           missing_before = missing_before,
#           missing_after = missing_after,
#           filled = missing_before - missing_after
#         )
#       }
#     }
#   }
#
#   # Summary report
#   cat("\n=== PROCESSING SUMMARY ===\n")
#   total_unsuitable <- sum(result_data$unsuitable_horizon)
#   cat("Unsuitable horizons excluded:", total_unsuitable, "\n")
#
#   for (prop in names(processing_log)) {
#     log_entry <- processing_log[[prop]]
#     cat(sprintf("%s: %d missing → %d missing (%d filled)\n",
#                 prop, log_entry$missing_before, log_entry$missing_after, log_entry$filled))
#   }
#
#   return(result_data)
# }
#
# # 8. QUALITY CONTROL AND VALIDATION
# # Quality control functions
# check_data_quality <- function(processed_data) {
#
#   # Check for logical consistency in texture
#   if (all(c("sandtotal_r", "claytotal_r", "silttotal_r") %in% names(processed_data))) {
#     texture_sums <- processed_data$sandtotal_r + processed_data$claytotal_r + processed_data$silttotal_r
#     texture_issues <- abs(texture_sums - 100) > 5  # Allow 5% tolerance
#
#     if (any(texture_issues, na.rm = TRUE)) {
#       cat("Warning:", sum(texture_issues, na.rm = TRUE), "records with texture sum issues\n")
#     }
#   }
#
#   # Check range consistency (_l ≤ _r ≤ _h)
#   properties_to_check <- c("sandtotal", "claytotal", "silttotal", "dbovendry", "cec7", "ph1to1h2o")
#
#   for (prop in properties_to_check) {
#     r_col <- paste0(prop, "_r")
#     l_col <- paste0(prop, "_l")
#     h_col <- paste0(prop, "_h")
#
#     if (all(c(r_col, l_col, h_col) %in% names(processed_data))) {
#       # Check _l ≤ _r ≤ _h
#       l_issues <- processed_data[[l_col]] > processed_data[[r_col]]
#       h_issues <- processed_data[[h_col]] < processed_data[[r_col]]
#
#       if (any(l_issues, na.rm = TRUE) || any(h_issues, na.rm = TRUE)) {
#         cat("Warning:", prop, "has", sum(l_issues | h_issues, na.rm = TRUE), "range consistency issues\n")
#       }
#     }
#   }
# }

