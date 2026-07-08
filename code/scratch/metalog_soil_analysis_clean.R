# =============================================================================
# SIMPLIFIED METALOG SOIL ANALYSIS - WORKING VERSION
# This version works with base R and handles the key field-level analysis
# =============================================================================

# Try to load packages safely
safe_require <- function(package) {
  suppressWarnings(require(package, character.only = TRUE, quietly = TRUE))
}

# Package availability flags
HAS_DPLYR <- safe_require("dplyr")
HAS_GGPLOT2 <- safe_require("ggplot2")
HAS_RMETALOG <- safe_require("rmetalog")

# =============================================================================
# IMPROVED METALOG FUNCTIONS - MANUAL IMPLEMENTATION
# These functions provide more robust metalog fitting and evaluation
# =============================================================================

# Helper function: Quantile function for metalog
quantile_metalog <- function(a, y, term, bounds, boundedness) {
  f <- (y - 0.5)
  l <- log(y / (1 - y))
  x <- a[1] + a[2] * l + a[3] * f * l
  if (term > 3) {
    x <- x + a[4] * f
  }
  o <- 2
  e <- 2
  if (term > 4) {
    for (i in 5:term) {
      if (i %% 2 == 0) {
        x <- x + a[i] * f^e * l
        e <- e + 1
      }
      if (i %% 2 != 0) {
        x <- x + a[i] * f^o
        o <- o + 1
      }
    }
  }
  if (boundedness == "sl") {
    x <- bounds[1] + exp(x)
  }
  if (boundedness == "su") {
    x <- bounds[2] - exp(-x)
  }
  if (boundedness == "b") {
    x <- (bounds[1] + bounds[2] * exp(x)) / (1 + exp(x))
  }
  return(x)
}

# Helper function: PDF (density function) for metalog
pdf_metalog <- function(a, y, term, bounds, boundedness) {
  d <- y * (1 - y)
  f <- (y - 0.5)
  l <- log(y / (1 - y))
  x <- a[2] / d
  if (a[3] != 0) {
    x <- x + a[3] * ((f / d) + l)
  }
  if (term > 3) {
    x <- x + a[4]
  }
  e <- 1
  o <- 1
  if (term > 4) {
    for (i in 5:term) {
      if (i %% 2 != 0) {
        x <- x + ((o + 1) * a[i] * f^o)
        o <- o + 1
      }
      if (i %% 2 == 0) {
        x <- x + a[i] * (((f^(e + 1)) / d) + (e + 1) * (f^e) * l)
        e <- e + 1
      }
    }
  }
  x <- (x^(-1))
  if (boundedness != "u") {
    M <- quantile_metalog(a, y, term, bounds = bounds, boundedness = "u")
  }
  if (boundedness == "sl") {
    x <- x * exp(-M)
  }
  if (boundedness == "su") {
    x <- x * exp(M)
  }
  if (boundedness == "b") {
    x <- (x * (1 + exp(M))^2) / ((bounds[2] - bounds[1]) * exp(M))
  }
  return(x)
}

# Newton's method to find probability for given quantile
newtons_method <- function(coefficients, q, term, bounds, boundedness) {
  alpha_step <- 0.1  # Increased from 0.01 for faster convergence
  err <- 1e-06  # Slightly relaxed from 1e-07
  temp_err <- 0.1
  y_now <- 0.5
  i <- 1
  max_iter <- 1000  # Reduced from 10000 - if not converging by then, it won't
  
  while (temp_err > err && i <= max_iter) {
    # Quantile function
    quantile_val <- quantile_metalog(coefficients, y_now, term, bounds, boundedness)
    
    # PDF (density function)
    density_val <- pdf_metalog(coefficients, y_now, term, bounds, boundedness)
    
    # Newton's method iteration
    frist_function <- quantile_val - q
    derv_function <- density_val
    
    # Check for zero or invalid derivative
    if(is.na(derv_function) || derv_function == 0 || !is.finite(derv_function)) {
      return(NA)  # Silent failure, don't spam warnings
    }
    
    # Adaptive step size - reduce if overshooting
    if(i > 1 && temp_err > 0.01) {
      alpha_step <- 0.05
    }
    
    y_next <- y_now - alpha_step * (frist_function / derv_function)  # FIXED: Division not multiplication
    temp_err <- abs(y_next - y_now)
    
    # Ensure y stays within valid bounds
    y_next <- max(min(y_next, 0.99999), 1e-06)
    y_now <- y_next
    
    i <- i + 1
  }
  
  # Return NA if didn't converge, but silently
  if(i > max_iter) {
    return(NA)
  }
  
  return(y_now)
}

# Manual metalog quantile function
manual_qmetalog <- function(coefficients, y, term, boundedness = "u", bounds = c(NA, NA)) {
  # Input validation
  if (max(y) >= 1 | min(y) <= 0) {
    stop("Error: y must be a positive numeric vector between 0 and 1")
  }
  if (term < 2 | term%%1 != 0 | length(term) > 1) {
    stop("Error: term must be a single positive integer >= 2")
  }
  
  # Initialize the basis functions
  Y <- data.frame(y1 = rep(1, length(y)))         # y1 = constant 1
  Y$y2 <- log(y / (1 - y))                       # y2 = log(y / (1 - y))
  
  if (term > 2) {
    Y$y3 <- (y - 0.5) * Y$y2                     # y3 = (y - 0.5) * y2
  }
  if (term > 3) {
    Y$y4 <- (y - 0.5)                            # y4 = (y - 0.5)
  }
  if (term > 4) {
    for (i in 5:term) {
      y_col <- paste0("y", i)
      if (i %% 2 != 0) {  # Odd terms
        Y[[y_col]] <- Y$y4^(i %/% 2)
      } else {            # Even terms
        z_col <- paste0("y", i - 1)
        Y[[y_col]] <- Y$y2 * Y[[z_col]]
      }
    }
  }
  
  # Convert to matrix
  Y <- as.matrix(Y)
  
  # Apply the coefficients
  s <- Y %*% coefficients[1:term]
  
  # Handle boundedness
  if (boundedness == "sl") {
    s <- bounds[1] + exp(s)
  } else if (boundedness == "su") {
    s <- bounds[2] - exp(-s)
  } else if (boundedness == "b") {
    s <- (bounds[1] + (bounds[2] * exp(s))) / (1 + exp(s))
  }
  
  return(as.numeric(s))
}

# Manual metalog density function
dmetalog_manual <- function(coefficients, q, term = 3, bounds = c(NA, NA), boundedness = "u") {
  tryCatch({
    # Compute probabilities (y) for each quantile using Newton's method
    probabilities <- sapply(q, function(x) {
      result <- newtons_method(coefficients = coefficients, q = x, term = term, bounds = bounds, boundedness = boundedness)
      if(is.na(result)) return(NA)
      return(result)
    })
    
    # Check how many failed
    n_failed <- sum(is.na(probabilities))
    if(n_failed > 0) {
      pct_failed <- (n_failed / length(probabilities)) * 100
      # Only warn if more than 10% failed
      if(pct_failed > 10) {
        warning(paste0(round(pct_failed, 1), "% of probabilities could not be computed via Newton's method"))
      }
      # If too many failed, return all NA
      if(pct_failed > 50) {
        return(rep(NA, length(q)))
      }
    }
    
    # Compute the density (PDF) for each probability
    densities <- sapply(probabilities, function(y) {
      if(is.na(y)) return(NA)
      pdf_metalog(coefficients = coefficients, y = y, term = term, bounds = bounds, boundedness = boundedness)
    })
    
    return(densities)
  }, error = function(e) {
    warning(paste("dmetalog_manual failed:", e$message))
    return(rep(NA, length(q)))
  })
}

#' Load C_datasets Study SOC Data (Working version)
#' @return data.frame with SOC stock measurements
load_c_datasets_soc_data <- function() {
  
  cat("Loading C_datasets study SOC data...\n")
  
  # Check for data files
  locations_file <- "data/C_datasets_Study_Data/locations.csv"
  measurements_file <- "data/C_datasets_Study_Data/measurements.csv"
  
  if(!file.exists(locations_file) || !file.exists(measurements_file)) {
    stop("Data files not found in data/C_datasets_Study_Data/")
  }
  
  # Load data
  locations <- read.csv(locations_file, stringsAsFactors = FALSE)
  measurements <- read.csv(measurements_file, stringsAsFactors = FALSE)
  
  cat("Loaded", nrow(locations), "sampling locations\n")
  cat("Loaded", nrow(measurements), "SOC measurements\n")
  
  # Merge data
  combined_data <- merge(measurements, locations, by = "location_id", all.x = TRUE)
  
  # Extract field IDs
  combined_data$field_id <- gsub("([0-9]+[A-Z]+)[0-9]{4}$", "\\1", combined_data$location_id)
  no_match <- combined_data$field_id == combined_data$location_id
  if(any(no_match)) {
    combined_data$field_id[no_match] <- gsub("([0-9]+[A-Z])[0-9]+$", "\\1", combined_data$location_id[no_match])
  }
  
  combined_data$field_code <- paste(combined_data$site, combined_data$field_id, sep = "_")
  
  cat("Identified", length(unique(combined_data$field_id)), "unique fields\n")
  cat("Field codes:", paste(head(sort(unique(combined_data$field_id)), 10), collapse = ", "), "...\n")
  
  # Standardize column names
  names(combined_data) <- tolower(names(combined_data))
  
  # Calculate SOC stocks for 0-30 cm
  # SOC stock (Mg C/ha) = SOCc (%) × BD (g/cm³) × depth (cm) × 100
  if("socc" %in% names(combined_data) && "bd" %in% names(combined_data)) {
    combined_data$depth_thickness <- combined_data$sample_depth_max - combined_data$sample_depth_min
    combined_data$soc_stock <- combined_data$socc * combined_data$bd * combined_data$depth_thickness * 100
    
    # Get 0-30 cm data
    soc_0_30 <- combined_data[combined_data$sample_depth_max <= 30, ]
    
    if(nrow(soc_0_30) > 0) {
      # Aggregate by location
      soc_summary <- aggregate(soc_stock ~ location_id, data = soc_0_30, FUN = sum, na.rm = TRUE)
      names(soc_summary)[2] <- "soc_stock_0_30"
      
      # Get unique location info - check what columns are available
      cat("Available columns in combined_data:", paste(names(combined_data), collapse = ", "), "\n")
      
      # Get columns that actually exist
      available_cols <- intersect(c("location_id", "site", "field_id", "field_code"), names(combined_data))
      locations_unique <- combined_data[!duplicated(combined_data$location_id), available_cols]
      
      # Merge results
      final_data <- merge(soc_summary, locations_unique, by = "location_id")
      
      cat("Calculated SOC stocks (0-30 cm) for", nrow(final_data), "locations\n")
      cat("Available columns:", paste(names(final_data), collapse = ", "), "\n")
      
      return(final_data)
    }
  }
  
  cat("Warning: Could not calculate SOC stocks\n")
  return(NULL)
}

# Placeholder functions for compatibility
download_scan_data <- function(...) { return(NULL) }
download_kssl_data <- function(...) { return(NULL) }
process_scan_moisture_data <- function(...) { return(NULL) }
process_kssl_data <- function(...) { return(NULL) }
load_combined_soil_data <- function(...) { return(load_c_datasets_soc_data()) }
#' Load comprehensive soil data from multiple sources
#' @param include_c_datasets Logical, include C_datasets study data
#' @param include_scan Logical, include SCAN network data
#' @param include_kssl Logical, include KSSL laboratory data
#' @param scan_sites Vector of SCAN site codes
#' @param kssl_mlras Vector of MLRA codes for KSSL
#' @param kssl_series Vector of soil series for KSSL
#' @return Combined soil property data frame
load_comprehensive_soil_data <- function(include_c_datasets = TRUE,
                                        include_scan = FALSE,
                                        include_kssl = FALSE,
                                        scan_sites = NULL,
                                        kssl_mlras = NULL,
                                        kssl_series = NULL) {
  
  all_data <- data.frame()
  
  # Load C_datasets if requested
  if(include_c_datasets) {
    cat("Loading C_datasets SOC data...\n")
    c_data <- load_c_datasets_soc_data()
    if(!is.null(c_data) && nrow(c_data) > 0) {
      c_data$data_source <- "c_datasets"
      all_data <- rbind(all_data, c_data)
    }
  }
  
  # Load SCAN data if requested
  if(include_scan) {
    cat("Loading SCAN network data...\n")
    scan_data <- load_scan_data(scan_sites)
    if(!is.null(scan_data) && nrow(scan_data) > 0) {
      scan_data$data_source <- "scan"
      all_data <- rbind(all_data, scan_data)
    }
  }
  
  # Load KSSL data if requested
  if(include_kssl) {
    cat("Loading KSSL laboratory data...\n")
    kssl_data <- load_kssl_data(kssl_mlras, kssl_series)
    if(!is.null(kssl_data) && nrow(kssl_data) > 0) {
      kssl_data$data_source <- "kssl"
      all_data <- rbind(all_data, kssl_data)
    }
  }
  
  if(nrow(all_data) == 0) {
    stop("No soil data could be loaded from any source")
  }
  
  return(all_data)
}

#' Load combined soil data (simplified version for fallback)
#' @param include_scan Logical, include SCAN data
#' @param scan_sites Vector of SCAN site codes
#' @return Combined soil data
load_combined_soil_data <- function(include_scan = TRUE, scan_sites = NULL) {
  
  # Try C_datasets first
  c_data <- load_c_datasets_soc_data()
  
  # If SCAN requested, try to add it
  if(include_scan) {
    scan_data <- load_scan_data(scan_sites)
    if(!is.null(scan_data) && nrow(scan_data) > 0) {
      scan_data$data_source <- "scan"
      c_data$data_source <- "c_datasets"
      
      # Try to combine if possible
      tryCatch({
        return(rbind(c_data, scan_data))
      }, error = function(e) {
        cat("Could not combine SCAN and C_datasets data, using C_datasets only\n")
        return(c_data)
      })
    }
  }
  
  return(c_data)
}

#' Load SCAN network soil moisture and temperature data using soilDB
#' @param site_codes Vector of SCAN site codes 
#' @return Data frame with soil measurements
load_scan_data <- function(site_codes = c(2002, 2113, 2129)) {
  
  if(!safe_require("soilDB")) {
    cat("soilDB package not available for SCAN data\n")
    return(NULL)
  }
  
  tryCatch({
    cat("Downloading SCAN data for sites:", paste(site_codes, collapse = ", "), "\n")
    
    # Use soilDB to fetch SCAN data - returns a LIST not a data frame
    scan_data <- soilDB::fetchSCAN(site.code = site_codes, 
                                   year = 2020:2023,
                                   report = 'SOIL')
    
    if(is.null(scan_data)) {
      cat("No SCAN data retrieved\n")
      return(NULL)
    }
    
    # List of reports to process and their corresponding property names
    reports_to_process <- list(
      SMS = "moisture",  # Soil Moisture
      STO = "temperature" # Soil Temperature
    )
    
    all_processed_data <- data.frame()
    
    for (report_code in names(reports_to_process)) {
      property_name <- reports_to_process[[report_code]]
      
      cat("Processing report:", report_code, "for property:", property_name, "\n")
      
      if (!is.null(scan_data[[report_code]]) && nrow(scan_data[[report_code]]) > 0) {
        
        report_df <- scan_data[[report_code]]
        cat("  Retrieved", nrow(report_df), "observations for", report_code, "\n")
        
        # Clean and filter the data
        valid_data <- report_df[!is.na(report_df$value) & is.finite(report_df$value), ]
        
        # Additional filtering for moisture (0-100%)
        if (report_code == "SMS") {
          valid_data <- valid_data[valid_data$value >= 0 & valid_data$value <= 100, ]
        }
        
        if (nrow(valid_data) > 0) {
          # Convert to our standard format
          processed_data <- data.frame(
            site_id = paste0("SCAN_", valid_data$Site),
            location_id = paste0("SCAN_", valid_data$Site, "_", valid_data$depth, "cm"),
            property_type = property_name,
            property_value = valid_data$value,
            depth_cm = valid_data$depth,
            measurement_date = valid_data$Date,
            stringsAsFactors = FALSE
          )
          
          all_processed_data <- rbind(all_processed_data, processed_data)
          cat("  Processed", nrow(processed_data), "valid observations for", property_name, "\n")
        } else {
          cat("  No usable data found for", report_code, "after filtering\n")
        }
      } else {
        cat("  No data found in SCAN response for report:", report_code, "\n")
      }
    }
    
    if (nrow(all_processed_data) > 0) {
      return(all_processed_data)
    }
    
    cat("No usable SCAN data found after processing all reports\n")
    return(NULL)
    
  }, error = function(e) {
    cat("Error loading SCAN data:", e$message, "\n")
    return(NULL)
  })
}#' Load KSSL laboratory data using soilDB
#' @param mlra_codes Vector of MLRA codes
#' @param series_names Vector of soil series names
#' @return Data frame with lab measurements
load_kssl_data <- function(mlra_codes = c("103", "104", "107"), 
                          series_names = c("clarion", "nicollet", "webster")) {
  
  if(!safe_require("soilDB")) {
    cat("soilDB package not available for KSSL data\n")
    return(NULL)
  }
  
  if(!safe_require("aqp")) {
    cat("aqp package not available for SoilProfileCollection handling\n")
    return(NULL)
  }
  
  tryCatch({
    cat("Downloading KSSL data for MLRAs:", paste(mlra_codes, collapse = ", "), "\n")
    cat("Soil series:", paste(series_names, collapse = ", "), "\n")
    
    # Use soilDB to fetch KSSL data - returns SoilProfileCollection
    kssl_result <- soilDB::fetchKSSL(mlra = mlra_codes, 
                                     series = series_names)
    
    if(is.null(kssl_result)) {
      cat("KSSL result is NULL\n")
      return(NULL)
    }
    
    cat("KSSL returned SoilProfileCollection with", length(kssl_result), "profiles\n")
    
    # Extract horizon and site data from SoilProfileCollection
    horizons_df <- as.data.frame(aqp::horizons(kssl_result))
    sites_df <- as.data.frame(aqp::site(kssl_result))
    
    cat("Extracted", nrow(horizons_df), "horizons and", nrow(sites_df), "sites\n")
    
    # Merge site and horizon data
    combined <- merge(horizons_df, sites_df, by = "pedon_key", all.x = TRUE)
    cat("Successfully merged data:", nrow(combined), "rows\n")
    
    if(nrow(combined) == 0) {
      cat("No combined data available\n")
      return(NULL)
    }
    
    # Process key soil properties based on the actual KSSL structure
    processed_data <- data.frame()
    
    cat("Processing soil properties...\n")
    cat("Available horizon columns:", paste(head(names(horizons_df), 20), collapse = ", "), "...\n")
    
    # Updated property mappings based on actual KSSL column names
    properties <- list(
      clay = "clay",
      silt = "silt", 
      sand = "sand",
      organic_carbon = "oc",
      ph_h2o = "ph_h2o",
      ph_cacl2 = "ph_cacl2",
      bulk_density = "db_13b",  # 1/3 bar bulk density
      water_retention_15 = "w15l2",  # 15 bar water content
      caco3 = "caco3"  # Calcium carbonate
    )
    
    cat("Searching for properties...\n")
    
    for(prop_name in names(properties)) {
      prop_col <- properties[[prop_name]]
      
      if(prop_col %in% names(combined)) {
        cat("Processing", prop_name, "from column:", prop_col, "\n")
        
        # Get valid data indices
        valid_indices <- !is.na(combined[[prop_col]]) & is.finite(combined[[prop_col]])
        n_valid <- sum(valid_indices)
        
        cat("  Valid values:", n_valid, "out of", nrow(combined), "\n")
        
        if(n_valid > 10) {
          # Extract data for valid indices only
          values <- combined[[prop_col]][valid_indices]
          
          # Create location identifiers
          site_ids <- combined$pedon_key[valid_indices]
          location_ids <- combined$labsampnum[valid_indices]
          depths_top <- combined$hzn_top[valid_indices]
          depths_bot <- combined$hzn_bot[valid_indices]
          mlras <- combined$mlra[valid_indices]
          series <- combined$taxonname[valid_indices]
          
          # Create processed data frame
          prop_data <- data.frame(
            site_id = paste0("KSSL_", site_ids),
            location_id = paste0("KSSL_", location_ids),
            property_type = prop_name,
            property_value = values,
            depth_top = depths_top,
            depth_bot = depths_bot,
            depth_cm = (depths_top + depths_bot) / 2,  # midpoint depth
            mlra = mlras,
            series = series,
            stringsAsFactors = FALSE
          )
          
          processed_data <- rbind(processed_data, prop_data)
          cat("  Added", nrow(prop_data), "records for", prop_name, "\n")
        } else {
          cat("  Insufficient valid data for", prop_name, "\n")
        }
      } else {
        cat("  Column", prop_col, "not found for", prop_name, "\n")
      }
    }
    
    if(nrow(processed_data) > 0) {
      cat("Processed", nrow(processed_data), "KSSL measurements\n")
      cat("Properties found:", paste(unique(processed_data$property_type), collapse = ", "), "\n")
      return(processed_data)
    } else {
      cat("No usable KSSL data found\n")
      return(NULL)
    }
    
  }, error = function(e) {
    cat("Error loading KSSL data:", e$message, "\n")
    return(NULL)
  })
}

fit_metalog_distribution <- function(values, property_name, bounds = NULL, max_terms = 15) {
  if(!HAS_RMETALOG) {
    cat("rmetalog package not available\n")
    return(NULL)
  }
  
  # Remove missing values
  values <- values[!is.na(values) & is.finite(values)]
  
  if(length(values) < 10) {
    cat("Insufficient data for metalog fitting (n =", length(values), ")\n")
    return(NULL)
  }
  
  cat("Fitting metalog distribution for", property_name, "(n =", length(values), ")\n")
  
  # PRE-FLIGHT CHECK: Skip if data characteristics predict timeout
  data_range <- max(values) - min(values)
  data_mean <- mean(values)
  range_mean_ratio <- data_range / data_mean
  
  # Empirical thresholds from IL-BR and IL-RT testing:
  # - High variability (range/mean > 1.0) AND small sample (< 200) = likely timeout
  # - OR very high skew (|skew| > 1.0)
  skip_metalog <- FALSE
  skip_reason <- ""
  
  if(range_mean_ratio > 1.0 && length(values) < 200) {
    skip_metalog <- TRUE
    skip_reason <- paste0("High variability (range/mean = ", round(range_mean_ratio, 2), 
                         ") with small sample (n = ", length(values), ")")
  }
  
  if(!skip_metalog) {
    skew_val <- tryCatch(moments::skewness(values), error = function(e) 0)
    if(abs(skew_val) > 1.2) {
      skip_metalog <- TRUE
      skip_reason <- paste0("High skewness (", round(skew_val, 2), ")")
    }
  }
  
  if(skip_metalog) {
    cat("  SKIPPING metalog fit:", skip_reason, "\n")
    cat("  Empirical testing shows this combination likely causes timeouts\n")
    return(NULL)
  }
  
  tryCatch({
    # Determine bounds based on data
    # UPDATED: Use fully bounded metalog as primary approach based on empirical testing
    # Bounded metalogs show 0.07-0.66% AIC improvement and 11-33% KS improvement in 4/5 sites
    if(is.null(bounds)) {
      data_min <- min(values)
      data_max <- max(values)
      data_range <- data_max - data_min
      
      # Use fully bounded with data-driven bounds (10% buffer)
      lower_bound <- max(0, data_min - 0.1 * data_range)
      upper_bound <- data_max + 0.1 * data_range
      bounds <- c(lower_bound, upper_bound)
      boundedness <- "b"
    } else {
      if(is.na(bounds[2])) {
        boundedness <- "sl"
      } else if(is.na(bounds[1])) {
        boundedness <- "su"
      } else {
        boundedness <- "b"
      }
    }
    
    # Try to fit metalog using rmetalog
    cat("  Attempting rmetalog fit with bounds:", bounds, "boundedness:", boundedness, "\n")
    metalog_fit <- tryCatch({
      # Primary attempt: bounded with term_limit = 4 (optimal based on testing)
      # TIMEOUT PROTECTION: Prevent infinite loops in rmetalog
      cat("  Calling rmetalog::metalog with term_limit = 4, step_len = 0.01\n")
      
      setTimeLimit(cpu = 10, elapsed = 10, transient = TRUE)
      result <- rmetalog::metalog(values, term_limit = 4, bounds = bounds, 
                                   boundedness = boundedness, step_len = 0.01)
      setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
      
      cat("  Primary bounded fit succeeded\n")
      result
    }, error = function(e1) {
      setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
      
      if(grepl("time limit|timeout", e1$message, ignore.case = TRUE)) {
        cat("  Primary fit TIMEOUT (>10 sec) - trying with fewer terms\n")
      } else {
        cat("  Primary fit failed:", e1$message, "\n")
      }
      
      # Fallback 1: Try with term_limit = 3
      tryCatch({
        cat("  Trying bounded fit with term_limit = 3\n")
        setTimeLimit(cpu = 10, elapsed = 10, transient = TRUE)
        result <- rmetalog::metalog(values, term_limit = 3, bounds = bounds, 
                                     boundedness = boundedness, step_len = 0.01)
        setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
        cat("  Bounded fit with 3 terms succeeded\n")
        result
      }, error = function(e2) {
        setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
        cat("  Bounded fit with 3 terms failed\n")
        
        # Fallback 2: Try with term_limit = 2
        tryCatch({
          cat("  Trying bounded fit with term_limit = 2\n")
          setTimeLimit(cpu = 10, elapsed = 10, transient = TRUE)
          result <- rmetalog::metalog(values, term_limit = 2, bounds = bounds, 
                                       boundedness = boundedness, step_len = 0.01)
          setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
          cat("  Bounded fit with 2 terms succeeded\n")
          result
        }, error = function(e3) {
          setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
          cat("  All bounded fitting attempts failed/timed out\n")
          return(NULL)
        })
      })
    })
    
    if(is.null(metalog_fit)) {
      cat("  WARNING: Metalog fitting failed completely for this site. Skipping metalog analysis.\n")
      cat("  Continuing with traditional distributions only...\n\n")
      next  # Skip to next site
    }
    
    if(!is.null(metalog_fit)) {
      # Find best valid term BY MINIMUM AIC (not maximum term count)
      if("Validation" %in% names(metalog_fit) && "term" %in% names(metalog_fit$Validation)) {
        valid_terms <- metalog_fit$Validation$term
        
        # DATA CHARACTERISTIC CHECK: Predict problematic cases
      data_range <- max(values) - min(values)
      data_mean <- mean(values)
      range_mean_ratio <- data_range / data_mean
      high_variability <- range_mean_ratio > 1.2
      large_sample <- length(values) > 200
      
      # Skip custom AIC calculation if data characteristics predict timeout issues
      skip_custom_aic <- high_variability && large_sample
      
      if(skip_custom_aic) {
        cat("  High variability detected (range/mean =", round(range_mean_ratio, 2), ")\n")
        cat("  Skipping custom AIC calculation, using minimum validated term count\n")
        best_term <- min(valid_terms)
        best_aic <- NA
        cat("  Selected", best_term, "terms (conservative choice)\n")
      } else {
        # Calculate AIC for each valid term count with timeout protection
        cat("  Calculating AIC for each term count to select optimal...\n")
        term_aics <- numeric(length(valid_terms))
        names(term_aics) <- valid_terms
        aic_timeouts <- 0
        
        for(i in seq_along(valid_terms)) {
          k <- valid_terms[i]
          cat("    Testing", k, "terms...")
          col_name <- paste0("M", k)
          
          # TIMEOUT PROTECTION: Use setTimeLimit for each AIC calculation
          aic_result <- tryCatch({
            setTimeLimit(cpu = 5, elapsed = 5, transient = TRUE)
            
            result <- if(col_name %in% colnames(metalog_fit$M)) {
              # Calculate AIC using M matrix approach with error handling
              p_seq <- metalog_fit$M[, "y"]
              Q_seq <- metalog_fit$M[, col_name]
              
              # Sort for proper differentiation
              sort_idx <- order(Q_seq)
              q_sorted <- Q_seq[sort_idx]
              p_sorted <- p_seq[sort_idx]
              
              # Calculate derivative dQ/dp numerically at midpoints
              dq <- diff(q_sorted)
              dp <- diff(p_sorted)
              
              # Avoid division by zero
              valid_dp <- abs(dp) > 1e-10
              if(sum(valid_dp) < length(dp) * 0.5) {
                # Too many invalid derivatives
                Inf
              } else {
                dq_dp <- dq[valid_dp] / dp[valid_dp]
                q_mid <- (q_sorted[-1][valid_dp] + q_sorted[-length(q_sorted)][valid_dp]) / 2
                pdf_mid <- 1 / dq_dp
                
                # Interpolate to get PDF at observed values
                densities <- approx(x = q_mid, y = pdf_mid, xout = values, rule = 2)$y
                
                # Check for valid densities
                valid_densities <- !is.na(densities) & densities > 0 & is.finite(densities)
                pct_valid <- sum(valid_densities) / length(densities) * 100
                
                if(pct_valid < 80) {
                  # Insufficient valid densities
                  Inf
                } else {
                  log_lik <- sum(log(densities[valid_densities]))
                  if(!is.finite(log_lik)) {
                    Inf
                  } else {
                    aic_val <- 2 * k - 2 * log_lik
                    cat("    ", k, "terms: AIC =", round(aic_val, 2), 
                        " (", round(pct_valid, 1), "% valid densities)\n")
                    aic_val
                  }
                }
              }
            } else {
              Inf
            }
            
            setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
            result
          }, error = function(e) {
            setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
            if(grepl("time limit|timeout", e$message, ignore.case = TRUE)) {
              cat("    ", k, "terms: TIMEOUT (>5 sec)\n")
              aic_timeouts <<- aic_timeouts + 1
            } else {
              cat("    ", k, "terms: AIC calculation failed (", e$message, ")\n")
            }
            Inf
          })
          
          term_aics[i] <- aic_result
          
          # FALLBACK: If 2+ timeouts, stop trying and use minimum terms
          if(aic_timeouts >= 2) {
            cat("  Multiple timeouts detected - falling back to minimum term count\n")
            break
          }
        }
        
        # Select term count with MINIMUM AIC
        # Filter out Inf values (failed calculations)
        valid_aic_idx <- which(is.finite(term_aics))
        
        # FALLBACK: If no valid AICs or multiple timeouts, use minimum validated terms
        if(length(valid_aic_idx) == 0 || aic_timeouts >= 2) {
          cat("  AIC calculation unreliable - using minimum validated term count\n")
          best_term <- min(valid_terms)
          best_aic <- NA
          cat("  Selected", best_term, "terms (conservative fallback)\n")
        } else {
          best_idx <- valid_aic_idx[which.min(term_aics[valid_aic_idx])]
          best_term <- valid_terms[best_idx]
          best_aic <- term_aics[best_idx]
          cat("  Selected", best_term, "terms (minimum AIC =", round(best_aic, 2), ")\n")
        }
      }
        
        # Extract coefficients for the best term
        coeff_col <- paste0("a", best_term)
        if(coeff_col %in% colnames(metalog_fit$A)) {
          coefficients <- metalog_fit$A[, coeff_col]
          
          return(list(
            fit = metalog_fit, 
            best_terms = best_term, 
            feasible = TRUE, 
            property = property_name,
            coefficients = coefficients,
            bounds = bounds,
            boundedness = boundedness,
            aic = best_aic  # Store the calculated AIC
          ))
        }
      }
      
      # Fallback for different package structures
      return(list(
        fit = metalog_fit, 
        best_terms = 3, 
        feasible = TRUE, 
        property = property_name,
        coefficients = NULL,
        bounds = bounds,
        boundedness = boundedness,
        aic = NA
      ))
    }
    return(NULL)
  }, error = function(e) {
    cat("Error fitting metalog:", e$message, "\n")
    return(NULL)
  })
}

fit_traditional_distributions <- function(values, property_name) {
  
  # Remove missing values
  values <- values[!is.na(values) & is.finite(values)]
  
  if(length(values) < 10) {
    return(NULL)
  }
  
  cat("Fitting traditional distributions for", property_name, "(n =", length(values), ")\n")
  
  results <- list()
  
  # Load fitdistrplus if available
  if(!safe_require("fitdistrplus")) {
    cat("fitdistrplus package not available, fitting simple normal only\n")
    normal_fit <- list(
      distribution = "normal",
      parameters = list(mean = mean(values), sd = sd(values)),
      aic = 2*2 + 2*length(values)*log(sd(values)) + length(values)*log(2*pi) + length(values)
    )
    return(list(normal = normal_fit))
  }
  
  # Normal distribution
  tryCatch({
    fit_norm <- fitdistrplus::fitdist(values, "norm")
    results$normal <- fit_norm
    cat("  Successfully fitted normal distribution\n")
  }, error = function(e) {
    cat("  Failed to fit normal distribution:", e$message, "\n")
  })
  
  # Log-normal distribution (for positive values)
  if(all(values > 0)) {
    tryCatch({
      fit_lnorm <- fitdistrplus::fitdist(values, "lnorm")
      results$lognormal <- fit_lnorm
      cat("  Successfully fitted log-normal distribution\n")
    }, error = function(e) {
      cat("  Failed to fit log-normal distribution:", e$message, "\n")
    })
  }
  
  # Gamma distribution (for positive values)
  if(all(values > 0)) {
    tryCatch({
      fit_gamma <- tryCatch({
        # Try MLE with starting values first
        mean_val <- mean(values)
        var_val <- var(values)
        fitdistrplus::fitdist(values, "gamma", 
                              start = list(shape = mean_val^2/var_val, 
                                           rate = mean_val/var_val),
                              optim.method = "BFGS")
      }, error = function(e) {
        # Fall back to MME if MLE fails
        fitdistrplus::fitdist(values, "gamma", method = "mme")
      })
      results$gamma <- fit_gamma
      cat("  Successfully fitted gamma distribution\n")
    }, error = function(e) {
      cat("  Failed to fit gamma distribution:", e$message, "\n")
    })
  }
  
  # Beta distribution (for values between 0 and 1 or percentages)
  if(property_name %in% c("clay", "silt", "sand", "moisture") || 
     (all(values >= 0) && all(values <= 1))) {
    tryCatch({
      if(max(values) > 1) {
        # Convert percentages to proportions
        scaled_values <- values / 100
      } else {
        scaled_values <- values
      }
      
      # Adjust for boundary values
      scaled_values[scaled_values <= 0] <- 0.001
      scaled_values[scaled_values >= 1] <- 0.999
      
      # Calculate method of moments starting values for beta
      m <- mean(scaled_values)
      v <- var(scaled_values)
      
      # Check if variance is valid for beta (must be < m*(1-m))
      if(v < m * (1 - m)) {
        shape1_start <- m * ((m * (1 - m) / v) - 1)
        shape2_start <- (1 - m) * ((m * (1 - m) / v) - 1)
        
        # Try MLE with L-BFGS-B (better for bounded parameters)
        fit_beta <- fitdistrplus::fitdist(scaled_values, "beta",
                                          start = list(shape1 = shape1_start, 
                                                       shape2 = shape2_start),
                                          optim.method = "L-BFGS-B",
                                          lower = c(shape1 = 0.001, shape2 = 0.001))
      } else {
        # Fall back to method of moments if variance too large
        fit_beta <- fitdistrplus::fitdist(scaled_values, "beta", 
                                          method = "mme")
      }
      
      results$beta <- fit_beta
      cat("  Successfully fitted beta distribution\n")
    }, error = function(e) {
      cat("  Failed to fit beta distribution:", e$message, "\n")
    })
  }
  
  cat("  Successfully fitted", length(results), "traditional distributions\n")
  
  return(results)
}

compare_distribution_performance <- function(values, metalog_fit, traditional_fits) {
  
  if(is.null(metalog_fit) || is.null(traditional_fits)) {
    return(NULL)
  }
  
  results <- data.frame(
    distribution = character(),
    aic = numeric(),
    bic = numeric(),
    ks_statistic = numeric(),
    ks_p_value = numeric(),
    terms = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Clean values
  values <- values[!is.na(values) & is.finite(values)]
  n <- length(values)
  
  # Traditional distributions
  for(dist_name in names(traditional_fits)) {
    fit <- traditional_fits[[dist_name]]
    
    tryCatch({
      # Get AIC and BIC from fitdist object
      if("aic" %in% names(fit)) {
        aic_val <- fit$aic
        bic_val <- fit$aic + (log(n) - 2) * length(fit$estimate)
      } else {
        # Fallback for simple fits
        aic_val <- fit$aic
        bic_val <- aic_val + log(n) * 2
      }
      
      # Kolmogorov-Smirnov test
      ks_test <- tryCatch({
        if(dist_name == "beta" && max(values) > 1) {
          # For beta distribution with percentage data
          scaled_values <- pmax(0.001, pmin(0.999, values / 100))
          ks.test(scaled_values, function(x) do.call(paste0("p", fit$distname), 
                                                     c(list(q = x), as.list(fit$estimate))))
        } else if("distname" %in% names(fit) && "estimate" %in% names(fit)) {
          ks.test(values, function(x) do.call(paste0("p", fit$distname), 
                                              c(list(q = x), as.list(fit$estimate))))
        } else {
          # Simple normal case
          ks.test(values, "pnorm", mean(values), sd(values))
        }
      }, error = function(e) {
        list(statistic = NA, p.value = NA)
      })
      
      # Determine number of parameters
      n_params <- if("estimate" %in% names(fit)) length(fit$estimate) else 2
      
      results <- rbind(results, data.frame(
        distribution = dist_name,
        aic = aic_val,
        bic = bic_val,
        ks_statistic = ks_test$statistic,
        ks_p_value = ks_test$p.value,
        terms = n_params,
        stringsAsFactors = FALSE
      ))
      
    }, error = function(e) {
      cat("Error processing", dist_name, "distribution:", e$message, "\n")
    })
  }
  
  # Metalog distribution
  if(!is.null(metalog_fit) && metalog_fit$feasible) {
    tryCatch({
      cat("  Starting metalog performance calculation...\n")
      n <- length(values)
      metalog_obj <- metalog_fit$fit
      k <- metalog_fit$best_terms
      cat("  Using", k, "terms for", n, "observations\n")
      
      # Use the stored M matrix (quantiles) and Y vector (probabilities) from rmetalog
      # This avoids Newton's method convergence issues entirely
      
      # Calculate AIC using numerical differentiation
      cat("  Calculating AIC...\n")
      metalog_aic <- tryCatch({
        # Check if AIC was already calculated during term selection
        if(!is.null(metalog_fit$aic) && !is.na(metalog_fit$aic) && is.finite(metalog_fit$aic)) {
          cat("  Using pre-calculated AIC:", round(metalog_fit$aic, 2), "\n")
          metalog_fit$aic
        } else {
          # Calculate AIC using M matrix approach
          # Get the quantile column for best term (columns are named M2, M3, M4, etc.)
          col_name <- paste0("M", k)
        
        if(col_name %in% colnames(metalog_obj$M)) {
          # Extract stored quantiles and probabilities
          # The M matrix has columns m2, M2, m3, M3, ..., and "y" (probability)
          stored_quantiles <- metalog_obj$M[, col_name]
          stored_probs <- metalog_obj$M[, "y"]  # Probability is stored in M matrix, not Y
          
          # For each data value, approximate its PDF using numerical differentiation
          # PDF(x) ≈ 1 / (dQ/dp) where Q is quantile function
          
          # Sort stored quantiles for interpolation
          sort_idx <- order(stored_quantiles)
          q_sorted <- stored_quantiles[sort_idx]
          p_sorted <- stored_probs[sort_idx]
          
          # Calculate derivative dQ/dp numerically
          dq <- diff(q_sorted)
          dp <- diff(p_sorted)
          dq_dp <- dq / dp  # Derivative of quantile function
          
          # PDF is reciprocal: f(x) = 1 / (dQ/dp)
          # Associate each derivative with midpoint
          q_mid <- (q_sorted[-1] + q_sorted[-length(q_sorted)]) / 2
          pdf_mid <- 1 / dq_dp
          
          # Interpolate to get PDF at observed values
          densities <- approx(x = q_mid, y = pdf_mid, xout = values, rule = 2)$y
          
          # Check for valid densities
          valid_densities <- !is.na(densities) & densities > 0 & is.finite(densities)
          pct_valid <- sum(valid_densities) / length(densities) * 100
          
          if(pct_valid < 50) {
            cat("  Metalog AIC calculation: only", round(pct_valid, 1), "% valid densities - skipping\n")
            NA_real_
          } else {
            # Calculate log-likelihood using valid densities
            log_lik <- sum(log(densities[valid_densities]))
            if(is.finite(log_lik)) {
              aic_metalog <- 2 * k - 2 * log_lik
              cat("  Metalog AIC calculated successfully:", round(aic_metalog, 2), 
                  "(", round(pct_valid, 1), "% valid densities)\n")
              aic_metalog
            } else {
              cat("  Metalog AIC: infinite log-likelihood\n")
              NA_real_
            }
          }
          } else {
            cat("  Column", col_name, "not found in M matrix\n")
            cat("  Available columns:", paste(colnames(metalog_obj$M), collapse = ", "), "\n")
            NA_real_
          }
        }
      }, error = function(e) {
        cat("  Metalog AIC calculation failed:", e$message, "\n")
        NA_real_
      })
      
      metalog_bic <- if(!is.na(metalog_aic)) metalog_aic + log(n) * k else NA_real_
      
      # Calculate KS statistic using stored quantiles for CDF
      cat("  Calculating KS statistic...\n")
      ks_test_metalog <- tryCatch({
        col_name <- paste0("M", k)
        
        if(col_name %in% colnames(metalog_obj$M)) {
          stored_quantiles <- metalog_obj$M[, col_name]
          stored_probs <- metalog_obj$M[, "y"]  # Probability is in M matrix
          
          cat("    Creating CDF function with", length(stored_quantiles), "points\n")
          
          # Create CDF function using interpolation
          metalog_cdf <- function(x) {
            result <- approx(x = stored_quantiles, y = stored_probs, xout = x, rule = 2)$y
            return(result)
          }
          
          # Test CDF function before running KS test
          cat("    Testing CDF function...\n")
          test_cdf <- tryCatch({
            test_vals <- quantile(values, probs = c(0.25, 0.5, 0.75))
            test_result <- metalog_cdf(test_vals)
            cat("    CDF test successful\n")
            TRUE
          }, error = function(e) {
            cat("    CDF test failed:", e$message, "\n")
            FALSE
          })
          
          if(test_cdf) {
            cat("    Running KS test...\n")
            ks_result <- ks.test(values, metalog_cdf)
            cat("    KS test complete\n")
            ks_result
          } else {
            list(statistic = NA, p.value = NA)
          }
        } else {
          cat("    Column", col_name, "not found\n")
          list(statistic = NA, p.value = NA)
        }
      }, error = function(e) {
        cat("  Metalog KS test failed:", e$message, "\n")
        list(statistic = NA, p.value = NA)
      })
      
      cat("  Building results data frame...\n")
      
      results <- rbind(results, data.frame(
        distribution = "metalog",
        aic = metalog_aic,
        bic = metalog_bic,
        ks_statistic = ks_test_metalog$statistic,
        ks_p_value = ks_test_metalog$p.value,
        terms = k,
        stringsAsFactors = FALSE
      ))
      
      cat("  Metalog results added successfully\n")
      
    }, error = function(e) {
      cat("Error calculating metalog performance metrics:", e$message, "\n")
      cat("Traceback:", paste(sys.calls(), collapse = "\n"), "\n")
    })
  }
  
  return(results)
}

cat("Simplified metalog soil analysis functions loaded successfully\n")

#' Load C_datasets with specific property focus
#' @param property_filter Character indicating which property to focus on
#' @return Data frame in long format with property_type and property_value
load_c_datasets_property_data <- function(property_filter = "bulk_density") {
  
  cat("Loading C_datasets with property filter:", property_filter, "\n")
  
  # Load the measurements data
  measurements_file <- "data/C_datasets_Study_Data/measurements.csv"
  if(!file.exists(measurements_file)) {
    stop("C_datasets measurements.csv not found")
  }
  
  measurements <- read.csv(measurements_file)
  cat("Loaded", nrow(measurements), "measurement records\n")
  
  # Create property-specific datasets
  if(property_filter == "bulk_density") {
    # Focus on bulk density
    result_data <- measurements %>%
      filter(!is.na(BD)) %>%
      select(site, location_id, sample_depth_min, sample_depth_max, 
             property_type = "BD", property_value = BD) %>%
      mutate(property_type = "bulk_density_0_30")
      
  } else if(property_filter == "soc_concentration") {
    # Focus on SOC concentration
    result_data <- measurements %>%
      filter(!is.na(SOCc)) %>%
      select(site, location_id, sample_depth_min, sample_depth_max, 
             property_type = "SOCc", property_value = SOCc) %>%
      mutate(property_type = "soc_concentration_0_30")
      
  } else if(property_filter == "depth_gradient") {
    # Focus on depth patterns - use SOC but analyze by depth
    result_data <- measurements %>%
      filter(!is.na(SOCc)) %>%
      mutate(depth_category = ifelse(sample_depth_max <= 15, "shallow", 
                                   ifelse(sample_depth_max <= 30, "medium", "deep"))) %>%
      select(site, location_id, sample_depth_min, sample_depth_max, 
             property_type = "depth_category", property_value = SOCc) %>%
      mutate(property_type = paste0("soc_depth_", property_type))
      
  } else if(property_filter == "texture_simulation") {
    # REMOVED: Simulated texture data not allowed for empirical analysis
    stop("Texture simulation removed - only real empirical data permitted")
      
  } else if(property_filter == "integrated") {
    # Multiple properties combined
    bd_data <- measurements %>%
      filter(!is.na(BD)) %>%
      select(site, location_id, sample_depth_min, sample_depth_max, 
             property_value = BD) %>%
      mutate(property_type = "bulk_density")
    
    soc_data <- measurements %>%
      filter(!is.na(SOCc)) %>%
      select(site, location_id, sample_depth_min, sample_depth_max, 
             property_value = SOCc) %>%
      mutate(property_type = "soc_concentration")
    
    result_data <- rbind(bd_data, soc_data)
    
  } else {
    # Default to SOC stock
    result_data <- measurements %>%
      mutate(thickness = sample_depth_max - sample_depth_min,
             soc_stock_kg_m2 = SOCc * BD * thickness / 100) %>%
      filter(!is.na(soc_stock_kg_m2)) %>%
      select(site, location_id, sample_depth_min, sample_depth_max, 
             property_type = "soc_stock", property_value = soc_stock_kg_m2) %>%
      mutate(property_type = "soc_stock_0_30")
  }
  
  cat("Created property dataset with", nrow(result_data), "records for", property_filter, "\n")
  cat("Property types:", paste(unique(result_data$property_type), collapse = ", "), "\n")
  cat("Value range:", round(min(result_data$property_value, na.rm=T), 3), 
      "to", round(max(result_data$property_value, na.rm=T), 3), "\n")
  
  return(result_data)
}
