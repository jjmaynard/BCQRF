#' @title List of Simulation Functions
#'
#' @description This document provides an overview of the key functions used in soil property and profile depth simulations.
#'
#' @section Function List:
#'
#' X -  Moved to soil_infill_functions `impute_rfv_values` - Imputes missing or zero RFV (rock fragment volume) values using default or derived values.
#' 2. `remove_organic_layer` - Removes organic horizons (O horizons) and adjusts mineral soil depths accordingly.
#' X -  Moved to soil_infill_functions `infill_soil_data` - Fills missing values in soil texture, bulk density, and water retention data.
#' 4. `slice_and_aggregate_soil_data` - Slices soil data into 1 cm increments and aggregates by standard depth intervals.
#' 5. `tri_dist` - Samples random values from a triangular distribution given minimum, mode, and maximum.
#' 6. `sim_component_comp` - Simulates soil component composition using a triangular distribution.
#' 7. `simulate_correlated_triangular` - Generates correlated random samples from multiple triangular distributions.
#' 8. `adjust_depthwise_property_GP_Quant` - Adjusts simulated soil property values using Gaussian Process (GP) predictions.
#' 9. `adjust_soil_property_by_depth` - Adjusts soil property values for a given soil component across depth layers.
#' 10. `adjust_soil_data` - Adjusts multiple soil properties using GP modeling and reshapes data.
#' 11. `adjust_soil_data_parallel` - Parallelized version of `adjust_soil_data` for faster processing.
#' 12. `simulate_soil_properties` - Simulates soil properties using correlation matrices and Monte Carlo methods.
#' 13. `calculate_mode` - Computes the mode (most frequent value) from a numeric vector.
#' 14. `simulate_cokey` - Simulates soil properties for a given soil component (`cokey`).
#' 15. `van_genuchten` - Computes soil water retention based on the van Genuchten model.
#' 16. `simulate_vg_aws` - Simulates plant-available water storage (AWS) using van Genuchten parameters.
#' 17. `calculate_aws_df` - Calculates AWS in the top 100 cm of soil based on simulated soil properties.
#' 18. `get_aws_data_by_mukey` - Queries SSURGO soil data for specified map unit keys (mukeys).
#' 19. `query_osd_distinctness` - Retrieves and converts OSD horizon distinctness codes to boundary offset values.
#' 20. `infill_missing_distinctness` - Fills missing distinctness values based on generalized soil horizons.
#' 21. `infill_missing_depth_variability` - Imputes missing soil depth values with representative ±2 cm adjustments.
#' 22. `simulate_soil_profile_top_down` - Simulates soil profile horizon depths using a top-down approach.
#' 23. `simulate_soil_profile_bottom_up` - Simulates soil profile horizon depths using a bottom-up approach.
#' 24. `simulate_soil_profile_thickness` - Estimates horizon thickness variability via multiple simulations.
#' 25. `simulate_and_perturb_soil_profiles` - Generates perturbed soil profiles by adjusting thickness and boundary distinctness.
#' 26. `simulate_profile_depths_by_collection` - Runs depth simulations for all profiles in a soil profile collection.
#' 27. `simulate_profile_depths_by_collection_parallel` - Parallelized version of `simulate_profile_depths_by_collection`.
#' 28. `simulate_profile_depths_by_mukey` - Retrieves SSURGO data and simulates soil profile depths for a given map unit key (`mukey`).
#' 29. `evaluate_simulated_depths` - Evaluates if simulated soil profile depths fall outside specified depth ranges.
#'
#' @note This list includes functions for soil property simulations, depth variability modeling, and Monte Carlo-based uncertainty analysis.




#' Remove Organic Layers and Adjust Depths (Including Depth Estimates)
#'
#' This function removes rows with organic horizons (where 'hzname' contains a capital 'O')
#' and adjusts the depths ('hzdept_r', 'hzdepb_r', 'hzdept_l', 'hzdepb_l', 'hzdept_h', 'hzdepb_h')
#' within each group (grouped by 'compname'). The depth adjustments ensure that the remaining horizons'
#' depths are recalculated based on the removal of the organic layers, preserving the original horizon thicknesses.
#' If the '_l' or '_h' values are missing or the columns do not exist, only the '_r' values will be adjusted.
#'
#' @param df A data frame containing soil horizon data, with columns 'hzname', 'hzdept_r', 'hzdepb_r',
#'           and optionally 'hzdept_l', 'hzdepb_l', 'hzdept_h', 'hzdepb_h', and 'compname'.
#' @return A data frame with organic layers removed and depth adjustments made within each group.
#' @export
remove_organic_layer <- function(df) {

  # Define a helper function to process each group.
  process_group <- function(group) {
    # Filter out rows where 'hzname' contains 'O' (Organic layers)
    filtered_group <- group[!grepl("O", group$hzname), ]

    # If the group is empty after filtering, return it directly.
    if (nrow(filtered_group) == 0) {
      return(filtered_group)
    }

    # Extract the first value of hzdept_r to use as a baseline.
    initial_hzdept_r <- filtered_group$hzdept_r[1]

    # Subtract the baseline from all depth columns using dplyr::mutate() and dplyr::across()
    filtered_group <- filtered_group |>
      dplyr::mutate(dplyr::across(dplyr::starts_with("hzdep"), ~ . - initial_hzdept_r))

    # Adjust 'hzdept_l' and 'hzdept_h' to 0 if 'hzdept_r' is 0.
    for (i in seq_len(nrow(filtered_group))) {
      if (filtered_group$hzdept_r[i] == 0) {
        if ("hzdept_l" %in% names(filtered_group)) {
          if (!is.na(filtered_group$hzdept_l[i])) {
            filtered_group$hzdept_l[i] <- 0
          }
        }
        if ("hzdept_h" %in% names(filtered_group)) {
          if (!is.na(filtered_group$hzdept_h[i])) {
            filtered_group$hzdept_h[i] <- 0
          }
        }
      }
    }

    # Calculate thickness for the 'hzdept_r' and 'hzdepb_r' columns.
    thickness_r <- filtered_group$hzdepb_r - filtered_group$hzdept_r

    # Adjust the 'hzdept_r' and 'hzdepb_r' depths.
    filtered_group$hzdept_r <- cumsum(c(0, thickness_r[-length(thickness_r)]))
    filtered_group$hzdepb_r <- cumsum(thickness_r)

    # If 'hzdept_l' and 'hzdepb_l' exist and are not all missing, adjust their depths.
    if (all(c("hzdept_l", "hzdepb_l") %in% names(filtered_group))) {
      if (!all(is.na(filtered_group$hzdept_l))) {
        thickness_l <- filtered_group$hzdepb_l - filtered_group$hzdept_l
        filtered_group$hzdept_l <- cumsum(c(0, thickness_l[-length(thickness_l)]))
        filtered_group$hzdepb_l <- cumsum(thickness_l)
      }
    }

    # If 'hzdept_h' and 'hzdepb_h' exist and are not all missing, adjust their depths.
    if (all(c("hzdept_h", "hzdepb_h") %in% names(filtered_group))) {
      if (!all(is.na(filtered_group$hzdept_h))) {
        thickness_h <- filtered_group$hzdepb_h - filtered_group$hzdept_h
        filtered_group$hzdept_h <- cumsum(c(0, thickness_h[-length(thickness_h)]))
        filtered_group$hzdepb_h <- cumsum(thickness_h)
      }
    }

    return(filtered_group)
  }

  # Group the data by 'cokey', apply process_group to each group, and then ungroup the result.
  result <- df |>
    dplyr::group_by(cokey) |>
    dplyr::group_modify(~ process_group(.x)) |>
    dplyr::ungroup()

  return(result)
}

#' Slice and Aggregate Soil Data at Specified Depth Intervals
#'
#' This function slices a data frame with soil data into 1 cm increments based on
#' depth ranges provided in the 'hzdept_r' and 'hzdepb_r' columns, and calculates
#' mean values for each depth increment across all other numeric data columns.
#'
#' @param df A data frame where each row represents a soil sample with 'hzdept_r'
#'        and 'hzdepb_r' columns representing depth ranges.
#' @param depth_ranges A list of vectors, where each vector contains two elements
#'        representing the top and bottom of a depth range (e.g., list(c(0, 30), c(30, 100))).
#'        Default is list(c(0, 30), c(30, 100)).
#' @return A data frame with depth ranges and mean values of soil properties for each range.
#' @export
slice_and_aggregate_soil_data <- function(df, depth_ranges = list(c(0, 30), c(30, 100))) {

  # Validate input
  if (!all(c("hzdept_r", "hzdepb_r") %in% colnames(df))) {
    stop("Input data frame must contain 'hzdept_r' and 'hzdepb_r' columns")
  }

  # Select numeric columns for aggregation, excluding the depth range columns
  data_columns <- df |>
    dplyr::select(dplyr::where(is.numeric)) |>
    dplyr::select(-hzdept_r, -hzdepb_r) |>
    colnames()

  # Generate a data frame for each 1 cm increment within each row's depth range
  rows_list <- list()
  for (i in 1:nrow(df)) {
    row <- df[i, ]
    # Convert the selected columns to a named list
    row_data <- as.list(row[data_columns])

    # Create increments for the depth range
    for (depth in seq(row$hzdept_r, row$hzdepb_r - 1)) {
      # Add the 'Depth' value to the row data
      row_data$Depth <- depth
      # Append the modified row data to the list
      rows_list <- append(rows_list, list(as.data.frame(row_data, stringsAsFactors = FALSE)))
    }
  }

  # Combine the list of data frames into a single data frame
  aggregated_data <- do.call(rbind, rows_list)

  # Convert the Depth column to numeric
  aggregated_data$Depth <- as.numeric(aggregated_data$Depth)

  # Process each depth range
  results <- list()
  for (i in seq_along(depth_ranges)) {
    range <- depth_ranges[[i]]
    top <- range[1]
    bottom <- range[2]

    # Subset data for the current depth range
    # Use >= for the first depth (to include 0) and < for the bottom (to avoid overlap)
    if (top == 0) {
      subset <- aggregated_data |>
        dplyr::select(dplyr::where(is.numeric)) |>
        dplyr::filter(Depth >= top & Depth < bottom)
    } else {
      subset <- aggregated_data |>
        dplyr::select(dplyr::where(is.numeric)) |>
        dplyr::filter(Depth >= top & Depth < bottom)
    }

    # Only process if there's data in this range
    if (nrow(subset) > 0) {
      # Calculate the mean for each column in the subset
      actual_max_depth <- max(subset$Depth)
      mean_values <- colMeans(subset, na.rm = TRUE)
      mean_values["hzdept_r"] <- top

      # Set bottom depth to actual max depth if it's less than the target bottom
      if (actual_max_depth + 1 < bottom) {  # +1 because we use < instead of <=
        mean_values["hzdepb_r"] <- actual_max_depth + 1
      } else {
        mean_values["hzdepb_r"] <- bottom
      }

      # Include the compname (from the first row of aggregated_data) if it exists
      if ("compname" %in% colnames(aggregated_data)) {
        result_row <- c(compname = aggregated_data$compname[1], mean_values)
      } else {
        result_row <- mean_values
      }

      results <- append(results, list(result_row))
    } else {
      # Create a row with NA values for ranges with no data
      if ("compname" %in% colnames(aggregated_data)) {
        na_row <- c(compname = aggregated_data$compname[1],
                    rep(NA, length(data_columns) + 1),  # +1 for Depth column
                    hzdept_r = top,
                    hzdepb_r = bottom)
        names(na_row) <- c("compname", data_columns, "Depth", "hzdept_r", "hzdepb_r")
      } else {
        na_row <- c(rep(NA, length(data_columns) + 1),  # +1 for Depth column
                    hzdept_r = top,
                    hzdepb_r = bottom)
        names(na_row) <- c(data_columns, "Depth", "hzdept_r", "hzdepb_r")
      }

      results <- append(results, list(na_row))
    }
  }

  # Convert results to a data frame
  if (length(results) > 0) {
    result_df <- do.call(rbind, results) |> as.data.frame()

    # Convert numeric columns back to numeric (they may have become character during rbind)
    numeric_cols <- c(data_columns, "Depth", "hzdept_r", "hzdepb_r")
    for (col in numeric_cols) {
      if (col %in% colnames(result_df)) {
        result_df[[col]] <- as.numeric(result_df[[col]])
      }
    }
  } else {
    result_df <- data.frame()
  }

  return(result_df)
}


#' Triangular Distribution Sampling Function, code from the 'triangle' R package
#'
#' This function generates random samples from a triangular distribution
#' defined by the minimum value (`a`), maximum value (`b`), and the mode (`c`).
#'
#' @param n Integer, number of random samples to generate (default is 1).
#'          If `n` is a vector, the length of the vector is used.
#' @param a Numeric, the minimum value of the distribution (default is 0).
#' @param b Numeric, the maximum value of the distribution (default is 1).
#' @param c Numeric, the mode or most likely value of the distribution.
#'          It defaults to the midpoint `(a + b) / 2` if not provided.
#' @return A vector of random samples from the triangular distribution.
#'         If invalid inputs are provided, a vector of `NaN` is returned.
#' @examples
#' # Generate 10 random samples from a triangular distribution
#' samples <- tri_dist(10, a = 0, b = 1, c = 0.5)
#' print(samples)
#' @export
tri_dist <- function(n = 1, a = 0, b = 1, c = (a + b)/2) {

  # If n is a vector, take its length as the number of samples
  if (length(n) > 1)
    n <- length(n)

  # Ensure n is a valid number and positive
  if (n < 1 | is.na(n))
    stop(paste("invalid argument: n =", n))

  # Ensure n is an integer (in case it isn't)
  n <- floor(n)

  # Check for any missing (NA) or invalid parameters in a, b, or c
  if (any(is.na(c(a, b, c))))
    return(rep(NaN, times = n))

  # Check if the mode (c) is within the range [a, b]
  if (a > c | b < c)
    return(rep(NaN, times = n))

  # Check for infinite values in a, b, or c
  if (any(is.infinite(c(a, b, c))))
    return(rep(NaN, times = n))

  # Generate n uniform random numbers between 0 and 1
  p <- runif(n)

  # Determine which values of p should be mapped to the left side of the triangle (a to c)
  # and which values should be mapped to the right side (c to b)
  if (a != c) {
    # For general case where mode (c) is not equal to the minimum (a)
    i <- which((a + sqrt(p * (b - a) * (c - a))) <= c)
    j <- which((b - sqrt((1 - p) * (b - a) * (b - c))) > c)
  } else {
    # Special case when a == c (the distribution starts from the mode)
    i <- which((a + sqrt(p * (b - a) * (c - a))) < c)
    j <- which((b - sqrt((1 - p) * (b - a) * (b - c))) >= c)
  }

  # Calculate the left side of the triangular distribution for indices in i
  if (length(i) != 0)
    p[i] <- a + sqrt(p[i] * (b - a) * (c - a))

  # Calculate the right side of the triangular distribution for indices in j
  if (length(j) != 0)
    p[j] <- b - sqrt((1 - p[j]) * (b - a) * (b - c))

  # Return the vector of random samples
  return(p)
}


#' Simulate Soil Component Composition
#'
#' This function simulates soil component compositions from a given dataframe
#' that includes the low, mode (representative), and high values for each component.
#' It uses a triangular distribution to generate random values for each component
#' and adjusts the results to ensure the total sum is 100%.
#'
#' @param data Dataframe containing the soil component information. It must contain
#'        the columns: 'compname' (component names), 'comppct_r' (mode/representative value),
#'        'comppct_l' (low value), and 'comppct_h' (high value).
#' @param n_simulations Integer, number of simulations to run (default is 1000).
#' @return A dataframe where each row corresponds to a simulated composition
#'         and each column represents one component. The percentages in each row sum to 100.
#'
#' @examples
#' # Create an example dataframe
#' data_table <- data.frame(
#'   compname = c("Component1", "Component2", "Component3", "Component4"),
#'   comppct_r = c(30, 40, 20, 10),  # Mode (representative values)
#'   comppct_l = c(25, 35, 15, 5),   # Low values
#'   comppct_h = c(35, 45, 25, 15)   # High values
#' )
#' # Simulate 1000 compositions
#' simulated_df <- sim_component_comp(data_table, n_simulations = 1000)
#' head(simulated_df)
#' @export
sim_component_comp <- function(data, n_simulations = 1000) {

  data <- data |>
    dplyr::select(mukey, cokey, compname, comppct_l, comppct_r, comppct_h) |>
    dplyr::distinct()

  # Fill missing comppct_l and comppct_h based on comppct_r
  data$comppct_l[is.na(data$comppct_l)] <- data$comppct_r[is.na(data$comppct_l)] - 2
  data$comppct_h[is.na(data$comppct_h)] <- data$comppct_r[is.na(data$comppct_h)] + 2

  n_components <- nrow(data)  # Number of components
  compositions <- matrix(NA, nrow = n_simulations, ncol = n_components)

  # Generate n_simulations random samples from triangular distributions for each component
  for (i in 1:n_components) {
    low  <- data$comppct_l[i]
    mode <- data$comppct_r[i]
    high <- data$comppct_h[i]

    compositions[, i] <- tri_dist(n_simulations, a = low, b = high, c = mode)
  }

  # Normalize across each column (note: the comment indicates rows, but the code aggregates columns)
  sim_compositions <- apply(compositions, 2, function(col) {
    round(sum(col) / 100, 0)
  })

  data$sim_comppct <- sim_compositions

  return(data)
}


#' Simulate Correlated Samples from Triangular Distributions
#'
#' This function generates correlated random samples from multiple triangular distributions
#' with specified lower limits, modes, and upper limits. The correlations between the variables
#' are controlled by the input correlation matrix. The function uses a Cholesky decomposition
#' to induce the specified correlations between the generated samples.
#'
#' @param n Integer, the number of samples to generate.
#' @param params List of parameter sets for the triangular distributions, where each element
#'        is a vector of three values:
#'        - Lower limit (`a`), Mode (`b`), and Upper limit (`c`) for the triangular distribution.
#'        For example: `params = list(c(a1, b1, c1), c(a2, b2, c2), ...)`.
#' @param correlation_matrix A square matrix specifying the desired correlations between
#'        the variables. It must be positive semi-definite and of size equal to the number of triangular distributions.
#' @param random_seed Optional integer, to set a specific random seed for reproducibility.
#' @return A matrix of correlated samples where each column corresponds to samples from one triangular distribution.
#' @details
#' The function first generates uncorrelated standard normal variables and then uses a Cholesky
#' decomposition of the correlation matrix to introduce correlations. The correlated standard
#' normal variables are then transformed into uniform variables using the cumulative distribution
#' function (CDF) of the normal distribution. Finally, these uniform variables are transformed
#' into samples from triangular distributions using the inverse CDF of the triangular distribution.
#'
#' @examples
#' # Define parameters for three triangular distributions
#' params <- list(c(0, 5, 10),  # Triangle distribution with limits 0, 10 and mode 5
#'                c(1, 4, 6),   # Triangle distribution with limits 1, 6 and mode 4
#'                c(2, 3, 5))   # Triangle distribution with limits 2, 5 and mode 3
#'
#' # Define a correlation matrix
#' correlation_matrix <- matrix(c(1, 0.5, 0.3,
#'                                0.5, 1, 0.2,
#'                                0.3, 0.2, 1), nrow = 3)
#'
#' # Simulate 1000 correlated samples
#' samples <- simulate_correlated_triangular(1000, params, correlation_matrix, random_seed = 42)
#' head(samples)
#' @export
simulate_correlated_triangular <- function(n, params, correlation_matrix, random_seed = NULL) {

  # Optionally set the seed for reproducibility (base R function)
  if (!is.null(random_seed)) {
    set.seed(random_seed)
  }

  # Generate uncorrelated standard normal variables using MASS::mvrnorm (non‐base function)
  uncorrelated_normal <- MASS::mvrnorm(n, mu = rep(0, length(params)), Sigma = diag(length(params)))

  # Cholesky decomposition of the correlation matrix (base R function: chol)
  L <- chol(correlation_matrix)

  # Compute correlated normal variables by multiplying the uncorrelated normals with L
  correlated_normal <- uncorrelated_normal %*% L

  # Initialize a matrix to store the samples for each triangular distribution (base R function: matrix)
  samples <- matrix(NA, nrow = n, ncol = length(params))

  # Loop through each triangular distribution
  for (i in seq_along(params)) {
    a <- params[[i]][1]  # Lower limit of the triangular distribution
    b <- params[[i]][2]  # Mode (peak) of the triangular distribution
    c <- params[[i]][3]  # Upper limit of the triangular distribution

    # Get the correlated normal variable for this triangular distribution
    normal_var <- correlated_normal[, i]

    # Transform the normal variable to a uniform distribution using the CDF of the normal distribution (base R: pnorm)
    u <- pnorm(normal_var)

    # Handle the case of a degenerate triangular distribution (division by zero)
    if (c == a) {
      samples[, i] <- rep(a, n)
    } else {
      # Determine which values fall on the left side of the triangle
      condition <- u <= (b - a) / (c - a)

      # For values satisfying the condition (left side of the triangle)
      samples[condition, i] <- a + sqrt(u[condition] * (c - a) * (b - a))

      # For values not satisfying the condition (right side of the triangle)
      samples[!condition, i] <- c - sqrt((1 - u[!condition]) * (c - a) * (c - b))
    }
  }

  # Return the matrix of correlated samples
  return(samples)
}


#' Adjust Depthwise Property Values Using GP Predictions
#'
#' This function adjusts simulated property values across different depths to align with the relative changes predicted by a Gaussian Process (GP) model. The adjustment ensures that the simulated values reflect the trends indicated by the GP model while maintaining their statistical properties.
#'
#' @param simulated_values A numeric matrix of simulated property values, where each row corresponds to a depth and each column corresponds to a simulation.
#' @param gp_model A fitted Gaussian Process model used to predict property means at different depths.
#' @param depths A numeric vector of depth values corresponding to the rows of `simulated_values`.
#'
#' @return A numeric matrix of adjusted simulated property values with the same dimensions as `simulated_values`.
#'
#' @examples
#' # Assuming 'simulated_values' is your matrix of simulated data,
#' # 'gp_model' is your fitted GP model, and 'depths' is your depth vector:
#' adjusted_values <- adjust_depthwise_property_GP_Quant(simulated_values, gp_model, depths)
#'
#' @export
adjust_depthwise_property_GP_Quant <- function(simulated_values, gp_model, depths) {

  # Get GP-predicted means for each depth using the GPfit package's predict.GP function
  gp_pred_means <- GPfit::predict.GP(gp_model, xnew = as.matrix(depths))$Y_hat

  # Initialize adjusted values matrix with the original simulated values
  adjusted_values <- simulated_values
  surface_values <- simulated_values[1, ]

  # Compute the empirical cumulative distribution function (ECDF) of the surface depth from simulated values
  surface_ecdf <- ecdf(surface_values)
  surface_quantiles <- surface_ecdf(surface_values)

  # Loop through each depth starting from the second one
  for (i in 2:length(depths)) {

    # Extract simulated values at the previous and current depths
    simulated_prev <- adjusted_values[i - 1, ]
    simulated_curr <- adjusted_values[i, ]

    # Calculate the ratio of GP-predicted means between current and previous depths
    gp_mean_ratio <- gp_pred_means[i] / gp_pred_means[i - 1]

    # Initialize a vector to store adjusted simulated values for the current depth
    adjusted_curr <- numeric(length(simulated_curr))

    # Adjust each simulated value at the current depth
    for (j in 1:length(simulated_curr)) {
      # Get the quantile for the simulated value at the previous depth
      q <- surface_quantiles[j]

      # Find the corresponding value at the same quantile in the current depth's distribution
      quantile_value <- quantile(simulated_curr[j], probs = q)

      # Adjust the current value based on the GP mean ratio and SD adjustment factor
      adjusted_curr[j] <- quantile_value + (simulated_prev[j] * gp_mean_ratio - quantile_value)
    }

    # Calculate ECDF for the adjusted current depth values
    ecdf_corrected <- ecdf(adjusted_curr)

    # Adjust corrected values to match original ECDF at the current depth
    corrected_adjusted <- quantile(simulated_values[i, ], probs = ecdf_corrected(adjusted_curr))

    # Store the adjusted simulated values back into the matrix
    adjusted_values[i, ] <- corrected_adjusted
  }

  # Return the matrix of adjusted simulated values
  return(adjusted_values)
}


#' Process Soil Properties for a 'cokey' Group
#'
#' This function processes each specified soil property within a 'cokey' group by adjusting the property values across depths using GP modeling.
#'
#' @param data_group A data frame representing a group of data for a specific 'cokey'.
#' @param properties A character vector of soil property names to be processed.
#'
#' @return A data frame with adjusted properties for the 'cokey' group.
#' @export
adjust_soil_property_by_depth <- function(data_group, properties) {

  # Initialize an empty list to store results for each property
  adjusted_data_list <- list()

  # Loop through each property and apply the correction
  for (property in properties) {
    # Select the relevant columns and reshape for the current property
    property_matrix <- data_group |>
      dplyr::select(hzdept_r, dplyr::all_of(property), simulation_number) |>
      tidyr::pivot_wider(
        names_from = simulation_number,
        values_from = dplyr::all_of(property),
        values_fill = NA
      ) |>
      dplyr::arrange(hzdept_r)

    # Extract depth levels and property matrix
    prop_depths <- property_matrix$hzdept_r
    property_matrix <- property_matrix |>
      dplyr::select(-hzdept_r) |>
      as.matrix()

    # Compute mean values at each depth for the current property
    mean_values <- rowMeans(property_matrix, na.rm = TRUE)

    # Scale the depth values to (0, 1)
    scaled_prop_depths <- scale(prop_depths, center = min(prop_depths), scale = max(prop_depths) - min(prop_depths))

    # Fit Gaussian Process model (from GPfit package)
    gp_model <- GPfit::GP_fit(X = as.matrix(scaled_prop_depths), Y = mean_values)

    # Adjust the original matrix based on GP predictions.
    # Note: The function 'adjust_depthwise_property_GP_Quant' must be defined in your environment.
    adjusted_matrix <- adjust_depthwise_property_GP_Quant(property_matrix, gp_model, scaled_prop_depths)

    # Convert adjusted matrix back to a data frame and rename columns
    adjusted_df <- as.data.frame(adjusted_matrix)
    colnames(adjusted_df) <- paste0(property, "_sim", seq_len(ncol(adjusted_df)))

    # Add depth column back using dplyr::bind_cols
    adjusted_df <- dplyr::bind_cols(hzdept_r = prop_depths, adjusted_df)

    # Store the adjusted data for this property in the list
    adjusted_data_list[[property]] <- adjusted_df
  }

  # Combine all adjusted properties into a single data frame by joining on 'hzdept_r'
  adjusted_data_combined <- adjusted_data_list |>
    purrr::reduce(dplyr::full_join, by = "hzdept_r")

  # Return the combined data with all adjusted properties for the current 'cokey' group
  return(adjusted_data_combined)
}


#' Process Soil Data with Gaussian Process Modeling
#'
#' This function processes soil data by adjusting specified soil properties for each unique 'cokey' group.
#' It applies Gaussian Process (GP) modeling to adjust property values across depths, reshapes the data,
#' and adds any missing columns from the original dataset.
#'
#' @param data A data frame containing soil simulation data with columns including 'cokey', 'hzdept_r', 'simulation_number', and the soil properties.
#' @param properties A character vector of soil property names to be processed.
#'
#' @return A data frame with adjusted soil properties, reshaped to long format and including all necessary columns.
#'
#' @details
#' The function performs the following steps:
#' 1. Groups the data by 'cokey' and processes each group separately using 'adjust_soil_property_by_depth'.
#' 2. Reshapes the data back to long format.
#' 3. Adds any missing columns from the original data.
#'
#' @import dplyr
#' @import tidyr
#' @import GPfit
#' @import purrr
#'
#' @examples
#' properties <- c("bulk_density_third_bar", "water_retention_third_bar", "water_retention_15_bar",
#'                 "rfv", "sand_total", "silt_total", "clay_total")
#' final_processed_data <- adjust_soil_data(ssurgo_soil_sim_df, properties)
#' @export
adjust_soil_data <- function(data, properties) {

  # Apply the function to each 'cokey' group in the data
  processed_data <- data |>
    dplyr::group_by(cokey) |>
    dplyr::group_modify(~ adjust_soil_property_by_depth(.x, properties)) |>
    dplyr::ungroup()

  # Reshape the data using pivot_longer to move from wide to long format
  processed_data_long <- processed_data |>
    tidyr::pivot_longer(
      cols = -c(cokey, hzdept_r),
      names_to = c(".value", "simulation_number"),
      names_pattern = "(.*)_sim(\\d+)",
      values_drop_na = TRUE
    ) |>
    dplyr::mutate(simulation_number = as.integer(simulation_number)) |>
    dplyr::arrange(cokey, hzdept_r, simulation_number)

  # Select the columns missing from 'processed_data_long' in the original data
  missing_columns <- data |>
    dplyr::select(cokey, hzdept_r, simulation_number, compname, mukey, hzdepb_r, unique_id) |>
    dplyr::distinct()

  # Join the missing columns with 'processed_data_long'
  final_data <- processed_data_long |>
    dplyr::left_join(missing_columns, by = c("cokey", "hzdept_r", "simulation_number"))

  # Return the final processed data
  return(final_data)
}


#' Process Soil Data with Gaussian Process Modeling in Parallel
#'
#' This function processes soil data by adjusting specified soil properties for each unique 'cokey' group.
#' It applies Gaussian Process (GP) modeling to adjust property values across depths, reshapes the data,
#' and adds any missing columns from the original dataset.
#'
#' @param data A data frame containing soil simulation data with columns including 'cokey', 'hzdept_r', 'simulation_number', and the soil properties.
#' @param properties A character vector of soil property names to be processed.
#'
#' @return A data frame with adjusted soil properties, reshaped to long format and including all necessary columns.
#'
#' @details
#' The function performs the following steps:
#' 1. Groups the data by 'cokey' and processes each group separately using 'adjust_soil_property_by_depth' in parallel.
#' 2. Reshapes the data back to long format.
#' 3. Adds any missing columns from the original data.
#'
#' @import dplyr
#' @import tidyr
#' @import GPfit
#' @import purrr
#' @import furrr
#' @import future
#'
#' @examples
#' properties <- c("bulk_density_third_bar", "water_retention_third_bar", "water_retention_15_bar",
#'                 "rfv", "sand_total", "silt_total", "clay_total")
#' final_processed_data <- adjust_soil_data_parallel(ssurgo_soil_sim_df, properties)
#' @export
adjust_soil_data_parallel <- function(data, properties) {

  # Set up the parallel backend using future package's plan
  future::plan(future::multisession)  # Use multisession plan; adjust as needed

  # Apply the function to each 'cokey' group in the data in parallel
  processed_data <- data |>
    dplyr::group_by(cokey) |>
    tidyr::nest() |>
    dplyr::mutate(processed = furrr::future_map_dfr(data, ~ adjust_soil_property_by_depth(.x, properties))) |>
    dplyr::select(-data) |>
    tidyr::unnest(cols = processed) |>
    dplyr::ungroup()

  # Reshape the data using pivot_longer to move from wide to long format
  processed_data_long <- processed_data |>
    tidyr::pivot_longer(
      cols = -c(cokey, hzdept_r),
      names_to = c(".value", "simulation_number"),
      names_pattern = "(.*)_sim(\\d+)",
      values_drop_na = TRUE
    ) |>
    dplyr::mutate(simulation_number = as.integer(simulation_number)) |>
    dplyr::arrange(cokey, hzdept_r, simulation_number)

  # Select the columns missing from 'processed_data_long' in the original data
  missing_columns <- data |>
    dplyr::select(cokey, hzdept_r, simulation_number, compname, mukey, hzdepb_r, unique_id) |>
    dplyr::distinct()

  # Join the missing columns with 'processed_data_long'
  final_data <- processed_data_long |>
    dplyr::left_join(missing_columns, by = c("cokey", "hzdept_r", "simulation_number"))

  # Return the final processed data
  return(final_data)
}


#' Simulate Soil Properties
#'
#' This function simulates soil properties based on input data from a soil profile dataset.
#' It calculates a local soil property correlation matrix and simulates soil properties.
#'
#' @param data A data frame containing input soil horizon data.
#' @param global_cor_matrices A named list of global correlation matrices.
#' @param global_cor_txt_matrices A named list of global text correlation matrices.
#'
#' @return A data frame with simulated soil properties.
#'
#' @import dplyr
#' @import tidyr
#' @import GPfit
#' @import purrr
#' @import Hmisc
#' @import Matrix
#' @export
simulate_soil_properties <- function(data, global_cor_matrices, global_cor_txt_matrices) {

  # --- Step 1: Subset and Infill Data ---

  # Define the columns of interest.
  all_sim_columns <- c(
    "compname", "mukey", "cokey", "hzname", "sim_comppct", "hzdept_r", "hzdepb_r",
    "dbovendry_l", "dbovendry_r", "dbovendry_h",
    "wthirdbar_l", "wthirdbar_r", "wthirdbar_h",
    "wfifteenbar_l", "wfifteenbar_r", "wfifteenbar_h",
    "sandtotal_l", "sandtotal_r", "sandtotal_h",
    "silttotal_l", "silttotal_r", "silttotal_h",
    "claytotal_l", "claytotal_r", "claytotal_h",
    "rfv_l", "rfv_r", "rfv_h",
    "ph1to1h2o_l", "ph1to1h2o_r", "ph1to1h2o_h",
    "cec7_l", "cec7_r", "cec7_h",
    "om_l", "om_r", "om_h"
  )

  # These are the simulation data columns that also map to the lookup table later.
  sim_data_columns <- c(
    "hzdept_r", "hzdepb_r",
    "dbovendry_l", "dbovendry_r", "dbovendry_h",
    "wthirdbar_l", "wthirdbar_r", "wthirdbar_h",
    "wfifteenbar_l", "wfifteenbar_r", "wfifteenbar_h",
    "sandtotal_l", "sandtotal_r", "sandtotal_h",
    "silttotal_l", "silttotal_r", "silttotal_h",
    "claytotal_l", "claytotal_r", "claytotal_h",
    "rfv_l", "rfv_r", "rfv_h",
    "ph1to1h2o_l", "ph1to1h2o_r", "ph1to1h2o_h",
    "cec7_l", "cec7_r", "cec7_h",
    "om_l", "om_r", "om_h"
  )

  # Keep only the columns available in 'data'
  sim_columns <- intersect(all_sim_columns, colnames(data))
  sim_data_columns <- intersect(sim_data_columns, colnames(data))

  # Subset and infill the data (assuming infill_soil_data() is defined elsewhere)
  sim_data <- data %>%
    dplyr::select(dplyr::all_of(sim_columns)) %>%
    infill_soil_data()

  # Create or recode the genetic horizon identifier ('genhz')
  if (!"genhz" %in% colnames(sim_data)) {
    sim_data <- sim_data %>%
      dplyr::mutate(genhz = sub("^[0-9]+", "", hzname))  # remove any leading digits

    sim_data$genhz <- generalizeHz(
      sim_data$genhz,
      new = c("O", "A", "B", "C", "Cr", "R"),
      pattern = c("O", "^A", "^B", "^C", "^Cr", "^R"),
      ordered = TRUE
    )
  }

  # --- Step 2: Prepare Global Correlation Matrices ---

  # Lookup table to map simulation columns to correlation matrix variable names.
  lookup_table <- data.frame(
    sim_columns = c(
      "dbovendry_r", "wthirdbar_r", "wfifteenbar_r",
      "claytotal_r", "sandtotal_r", "rfv_r",
      "ph1to1h2o_r", "cec7_r", "om_r"
    ),
    cor_names = c(
      "db", "wr_3b", "wr_15b",
      "ilr1", "ilr2", "rfv",
      "ph", "cec", "soc"
    ),
    stringsAsFactors = FALSE
  )

  # Determine the correlation variables to include.
  cor_vars <- lookup_table$cor_names[lookup_table$sim_columns %in% sim_data_columns]

  # Subset each global correlation matrix to only the required variables.
  global_cor_matrices <- lapply(global_cor_matrices, function(mat) {
    mat[cor_vars, cor_vars, drop = FALSE]
  })

  # --- Step 3: Calculate Local Correlation Matrices ---

  # Select simulation properties (columns ending with "_r") along with genhz.
  rep_columns <- sim_data %>%
    dplyr::select(dplyr::ends_with("_r"), genhz)

  # Stop if there is not more than one soil property to simulate.
  if(ncol(rep_columns %>% dplyr::select(dplyr::ends_with("_r"))) < 2) {
    stop("More than one property is required for simulation")
  }

  # If all three soil texture columns exist, compute the ILR transformation.
  texture_cols <- c("sandtotal_r", "silttotal_r", "claytotal_r")
  if (all(texture_cols %in% colnames(sim_data))) {
    ilr_site <- sim_data %>%
      dplyr::select(dplyr::all_of(texture_cols)) %>%
      ilr() %>%    # assumes ilr() is defined elsewhere
      as.matrix()

    # Add ILR components to the simulation data.
    rep_columns$ilr1 <- ilr_site[, 1]
    rep_columns$ilr2 <- ilr_site[, 2]

    # Remove the original texture columns from rep_columns.
    rep_columns <- rep_columns %>%
      dplyr::select(-dplyr::all_of(texture_cols))
  }

  # Replace any zero values in the simulation columns (except excluded ones) with 0.01.
  exclude_cols <- c("hzdept_r", "hzdepb_r", "genhz")
  include_cols <- setdiff(names(rep_columns), exclude_cols)
  rep_columns <- rep_columns %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(include_cols), ~ ifelse(. == 0, 0.01, .)))

  # Ensure 'genhz' is a factor.
  rep_columns$genhz <- factor(rep_columns$genhz)

  # Helper function to check positive definiteness.
  is_positive_definite <- function(mat) {
    tryCatch({
      chol(mat)
      TRUE
    }, error = function(e) FALSE)
  }

  # Initialize lists to hold local correlation matrices.
  correlation_matrices <- list()
  txt_correlation_matrices_local <- list()

  # Loop through each unique genetic horizon.
  for (gh in unique(rep_columns$genhz)) {
    rep_columns_filtered <- rep_columns %>%
      dplyr::filter(genhz == gh, wthirdbar_r != 0.01)

    correlation_data <- rep_columns_filtered %>%
      dplyr::select(dplyr::all_of(include_cols))

    # Calculate local correlation if enough observations exist;
    # otherwise, fall back to the corresponding global matrix.
    if (nrow(correlation_data) > 4) {
      local_cor <- Hmisc::rcorr(as.matrix(correlation_data), type = "pearson")$r
    } else {
      local_cor <- global_cor_matrices[[as.character(gh)]]
    }

    # Force symmetry.
    local_cor <- (local_cor + t(local_cor)) / 2

    # Ensure the matrix is positive definite.
    if (!is_positive_definite(local_cor)) {
      if (any(is.na(local_cor)) || any(is.infinite(local_cor))) {
        local_cor <- global_cor_matrices[[as.character(gh)]]
      } else {
        local_cor <- as.matrix(Matrix::nearPD(local_cor)$mat)
      }
    }

    # Store the matrices using the genetic horizon as key.
    correlation_matrices[[as.character(gh)]] <- local_cor
    txt_correlation_matrices_local[[as.character(gh)]] <- global_cor_txt_matrices[[as.character(gh)]]
  }

  # --- Step 4: Simulate by cokey ---

  result_list <- list()
  unique_cokeys <- unique(sim_data$cokey)

  for (ck in unique_cokeys) {
    sim_subset <- sim_data %>% dplyr::filter(cokey == ck)

    # The function simulate_cokey_generalized is assumed to be defined elsewhere.
    sim_result <- simulate_cokey_generalized(sim_subset, correlation_matrices, txt_correlation_matrices_local)

    result_list[[as.character(ck)]] <- sim_result
  }

  # Combine and return the simulation results.
  sim_data_df <- do.call(rbind, result_list)
  return(sim_data_df)
}


#' Calculate Mode
#'
#' This function calculates the mode (i.e., the most frequently occurring value) of a given vector.
#'
#' @param x A vector of values.
#' @return The mode of the vector.
#' @examples
#' calculate_mode(c(1, 2, 2, 3, 3, 3, 4))
#' @export
calculate_mode <- function(x) {
  freq_table <- table(x)                # Create a frequency table
  mode_value <- as.numeric(names(freq_table)[which.max(freq_table)])  # Find the most frequent value
  return(mode_value)
}


#' Simulate Soil Properties for a Specific Cokey (Generalized)
#'
#' This function simulates a flexible set of soil properties for each row in a given cokey group.
#' If the necessary columns for texture (sand, silt, clay) are present, it also simulates texture
#' data using compositional transformations. The joint probability distributions are built from
#' correlation matrices corresponding to the genetic horizon (genhz).
#'
#' @param sim_cokey A data frame containing soil data for a specific cokey group.
#' @param correlation_matrices A list of correlation matrices keyed by genetic horizon (genhz) used for soil property simulation.
#' @param txt_correlation_matrices A list of correlation matrices keyed by genetic horizon (genhz) for soil texture (sand, silt, clay).
#'
#' @return A data frame containing simulated soil properties for the cokey group, including metadata
#' (compname, mukey, cokey, depth, simulation number, and a unique identifier).
#'
#' @import dplyr
#' @import compositions
#'
#' @examples
#' \dontrun{
#'   # Suppose 'sim_cokey_data' has columns for bulk density, pH, and texture.
#'   # correlation_matrices has correlation info for these properties by genhz,
#'   # and txt_correlation_matrices has correlation info for texture.
#'   sim_result <- simulate_cokey_generalized(sim_cokey_data, correlation_matrices, txt_correlation_matrices)
#' }
#' @export
simulate_cokey_generalized <- function(sim_cokey, correlation_matrices, txt_correlation_matrices = NULL) {

  # The fixed order in which properties must appear (rows/cols in correlation matrix)
  param_order <- c("db", "wr_3b", "wr_15b", "ilr1", "ilr2", "rfv", "ph", "cec", "soc")

  # Columns required for texture (sand, silt, clay)
  texture_cols_required <- c("sandtotal_l", "sandtotal_r", "sandtotal_h",
                             "silttotal_l", "silttotal_r", "silttotal_h",
                             "claytotal_l", "claytotal_r", "claytotal_h")

  # Helper to extract c(l, r, h) if all present
  get_param_set <- function(row, prefix) {
    lcol <- paste0(prefix, "_l")
    rcol <- paste0(prefix, "_r")
    hcol <- paste0(prefix, "_h")
    if (all(c(lcol, rcol, hcol) %in% names(row))) {
      return(c(row[[lcol]], row[[rcol]], row[[hcol]]))
    }
    return(NULL)
  }

  # Initialize
  sim_data_out <- list()

  for (i in seq_len(nrow(sim_cokey))) {
    row <- sim_cokey[i, ]
    genhz_val <- as.character(row$genhz)

    # Identify the local correlation matrix
    local_corr <- correlation_matrices[[genhz_val]]

    # Check if texture is fully present
    has_texture <- all(texture_cols_required %in% names(row))

    # 1) If texture is present, do the triangular simulation for sand/silt/clay
    #    get ilr1/ilr2 c(l,r,h).
    ilr1_lrh <- NULL
    ilr2_lrh <- NULL

    if (has_texture && !is.null(txt_correlation_matrices)) {
      txt_corr <- txt_correlation_matrices[[genhz_val]]
      # Triangular params for sand/silt/clay
      params_txt <- list(
        c(row$sandtotal_l, row$sandtotal_r, row$sandtotal_h),
        c(row$silttotal_l, row$silttotal_r, row$silttotal_h),
        c(row$claytotal_l, row$claytotal_r, row$claytotal_h)
      )
      # Simulate texture
      sim_txt <- compositions::acomp(
        simulate_correlated_triangular(as.integer(row$sim_comppct), params_txt, txt_corr)
      )
      # ILR transform
      sim_txt_ilr <- compositions::ilr(sim_txt)
      ilr1_vals <- sim_txt_ilr[, 1]
      ilr2_vals <- sim_txt_ilr[, 2]
      ilr1_lrh <- c(min(ilr1_vals), calculate_mode(ilr1_vals), max(ilr1_vals))
      ilr2_lrh <- c(min(ilr2_vals), calculate_mode(ilr2_vals), max(ilr2_vals))
    }

    # 2) Build a named list of property parameter sets in the EXACT param_order
    param_list <- vector("list", length(param_order))
    names(param_list) <- param_order

    # Fill in each property if present:

    # db
    db_set <- get_param_set(row, "dbovendry")
    if (!is.null(db_set)) param_list[["db"]] <- db_set

    # wr_3b
    wr3_set <- get_param_set(row, "wthirdbar")
    if (!is.null(wr3_set)) param_list[["wr_3b"]] <- wr3_set

    # wr_15b
    wr15_set <- get_param_set(row, "wfifteenbar")
    if (!is.null(wr15_set)) param_list[["wr_15b"]] <- wr15_set

    # ilr1 / ilr2 from texture
    if (has_texture && !is.null(ilr1_lrh) && !is.null(ilr2_lrh)) {
      param_list[["ilr1"]] <- ilr1_lrh
      param_list[["ilr2"]] <- ilr2_lrh
    }

    # rfv
    rfv_set <- get_param_set(row, "rfv")
    if (!is.null(rfv_set)) param_list[["rfv"]] <- rfv_set

    # ph
    ph_set <- get_param_set(row, "ph1to1h2o")
    if (!is.null(ph_set)) param_list[["ph"]] <- ph_set

    # cec
    cec_set <- get_param_set(row, "cec7")
    if (!is.null(cec_set)) param_list[["cec"]] <- cec_set

    # soc / om
    om_set <- get_param_set(row, "om")
    if (!is.null(om_set)) param_list[["soc"]] <- om_set

    # Remove NULL items
    param_list <- param_list[!vapply(param_list, is.null, logical(1))]

    # If no valid params, skip
    if (!length(param_list)) {
      message("No recognized properties for row ", i, ". Skipping.")
      next
    }

    # Convert param_list to the final structure for simulation
    params_for_sim <- unname(param_list)  # remove names so simulate_correlated_triangular sees a list

    # Subset correlation matrix to the columns we actually have
    keep_cols <- names(param_list)  # e.g. c("db","wr_3b","wr_15b","ilr1","ilr2","rfv","ph","cec","soc") if present
    local_corr_sub <- local_corr[keep_cols, keep_cols, drop = FALSE]

    # 3) Now run the simulation in the correct order
    tryCatch({
      n_sim <- as.integer(row$sim_comppct)
      sim_data <- simulate_correlated_triangular(
        n = n_sim,
        params = params_for_sim,
        correlation_matrix = local_corr_sub
      )
      sim_data <- as.data.frame(sim_data)

      # Assign property names based on param_list
      colnames(sim_data) <- names(param_list)

      # If we had water retention columns, do scaling
      if ("wr_3b" %in% names(param_list)) {
        sim_data[["wr_3b"]] <- sim_data[["wr_3b"]] / 100
      }
      if ("wr_15b" %in% names(param_list)) {
        sim_data[["wr_15b"]] <- sim_data[["wr_15b"]] / 100
      }

      # If texture was present, do inverse ILR
      if (has_texture && all(c("ilr1", "ilr2") %in% names(param_list))) {
        sim_txt <- compositions::ilrInv(sim_data[, c("ilr1", "ilr2")])
        sim_txt <- as.data.frame(sim_txt)
        colnames(sim_txt) <- c("sand_total", "silt_total", "clay_total")
        sim_txt <- sim_txt * 100

        # Remove ilr1 / ilr2
        sim_data <- sim_data[, setdiff(names(sim_data), c("ilr1", "ilr2"))]

        # Combine with final texture columns
        sim_data <- cbind(sim_data, sim_txt)
      }

      # Attach metadata
      sim_data$compname <- row$compname
      sim_data$mukey    <- row$mukey
      sim_data$cokey    <- row$cokey
      sim_data$hzdept_r <- row$hzdept_r
      sim_data$hzdepb_r <- row$hzdepb_r

      # Add simulation number & unique id
      sim_data$simulation_number <- seq_len(nrow(sim_data))
      sim_data$unique_id <- paste0(row$cokey, "-", sprintf("%02d", sim_data$simulation_number))

      sim_data_out <- append(sim_data_out, list(sim_data))
    }, error = function(e) {
      message("Error for row ", i, ": ", e$message)
    })
  }

  # Combine results
  final_df <- do.call(rbind, sim_data_out)
  return(final_df)
}


#' Simulate Soil Properties for a Specific Cokey
#'
#' This function simulates soil property values for each row within a given cokey group.
#' It performs the following steps depending on input data:
#'   1. If `texture` data is present:
#'      a.  Simulates sand, silt, and clay percentages using triangular distributions and converts the results
#'          using compositional transformations (via acomp and ilr).
#'      b.  Extracts lower, modal, and upper values from the ilr-transformed texture data.
#'   2. Constructs parameters for additional properties (bulk density, water retention at two pressure levels, and rfv)
#'      and simulates these properties using correlated triangular distributions.
#'   3. If `texture` data is present:
#'          a. Performs an inverse ilr transformation on the simulated texture data, scales the results,
#'          and assembles the complete simulated dataset for the cokey group.
#'
#' @param sim_cokey A data frame containing soil data for a specific cokey group.
#' @param correlation_matrices A list of correlation matrices keyed by genetic horizon (genhz) used for soil property simulation.
#' @param texture_correlation_matrices A list of correlation matrices keyed by genetic horizon (genhz) for soil texture (sand, silt, clay).
#'
#' @return A data frame containing simulated soil properties for the cokey group, including metadata (compname, mukey, cokey,
#'         depth, simulation number, and a unique identifier).
#'
#' @import dplyr
#' @import compositions
#'
#' @examples
#' \dontrun{
#'   sim_result <- simulate_cokey(sim_cokey_data, correlation_matrices, texture_correlation_matrix)
#' }
#' @export
simulate_cokey <- function(sim_cokey, correlation_matrices, txt_correlation_matrices) {
  # Initialize an empty list to store the simulated data for each row
  sim_data_out <- list()

  for (index in seq_len(nrow(sim_cokey))) {
    row <- sim_cokey[index, ]
    local_correlation_matrix <- correlation_matrices[[as.character(row$genhz)]]

    # If texture is present
    if (any(c('sandtotal_r', 'silttotal_r', 'claytotal_r') %in% names(row))) {
      texture_correlation_matrix <- txt_correlation_matrices[[as.character(row$genhz)]]
      # Step 2a: Simulate sand, silt, clay percentages
      params_txt <- list(
        c(row$sandtotal_l, row$sandtotal_r, row$sandtotal_h),
        c(row$silttotal_l, row$silttotal_r, row$silttotal_h),
        c(row$claytotal_l, row$claytotal_r, row$claytotal_h)
      )

      # Step 2b: Convert simulated texture data using acomp and ilr transformation
      simulated_txt <- compositions::acomp(simulate_correlated_triangular(as.integer(row$sim_comppct), params_txt, texture_correlation_matrix))
      simulated_txt_ilr <- compositions::ilr(simulated_txt)

      # Step 2c: Extract lower, modal, and upper values for ilr1 and ilr2
      ilr1_values <- simulated_txt_ilr[, 1]
      ilr2_values <- simulated_txt_ilr[, 2]

      ilr1_l <- min(ilr1_values)
      ilr1_r <- calculate_mode(ilr1_values)
      ilr1_h <- max(ilr1_values)
      ilr2_l <- min(ilr2_values)
      ilr2_r <- calculate_mode(ilr2_values)
      ilr2_h <- max(ilr2_values)
    }

    params <- list(
      c(row$dbovendry_l, row$dbovendry_r, row$dbovendry_h),
      c(row$wthirdbar_l, row$wthirdbar_r, row$wthirdbar_h),
      c(row$wfifteenbar_l, row$wfifteenbar_r, row$wfifteenbar_h),
      c(ilr1_l, ilr1_r, ilr1_h),
      c(ilr2_l, ilr2_r, ilr2_h),
      c(row$rfv_l, row$rfv_r, row$rfv_h)
    )

    # Step 2d: Simulate all properties and perform inverse ilr transformation
    tryCatch({
      n_sim <- as.integer(row$sim_comppct)
      sim_data <- simulate_correlated_triangular(
        n = n_sim,
        params = params,
        correlation_matrix = local_correlation_matrix
      )

      sim_data <- data.frame(sim_data)
      colnames(sim_data) <- c("bulk_density_third_bar", "water_retention_third_bar", "water_retention_15_bar", "ilr1", "ilr2", "rfv")

      sim_data$water_retention_third_bar <- sim_data$water_retention_third_bar / 100
      sim_data$water_retention_15_bar <- sim_data$water_retention_15_bar / 100

      sim_txt <- compositions::ilrInv(sim_data[, c("ilr1", "ilr2")])
      sim_txt <- data.frame(sim_txt)
      colnames(sim_txt) <- c("sand_total", "silt_total", "clay_total")
      sim_txt <- sim_txt * 100

      multi_sim <- cbind(sim_data[, !base::colnames(sim_data) %in% c("ilr1", "ilr2")], sim_txt)
      multi_sim$compname <- row$compname
      multi_sim$mukey <- row$mukey
      multi_sim$cokey <- row$cokey
      multi_sim$hzdept_r <- row$hzdept_r
      multi_sim$hzdepb_r <- row$hzdepb_r

      # Add simulation number
      multi_sim$simulation_number <- seq_len(nrow(multi_sim))

      # Create unique identifier: 'cokey' plus "_" and the simulation number
      multi_sim$unique_id <- paste0(row$cokey, "-", sprintf("%02d", multi_sim$simulation_number))

      sim_data_out <- base::append(sim_data_out, list(multi_sim))
    }, error = function(e) {
      base::cat("Error: ", e$message, "\n")
    })
  }

  # Concatenate simulated values into a single data frame
  sim_data_df <- dplyr::bind_rows(sim_data_out) |>
    base::replace(is.na(.), NA)

  return(sim_data_df)
}


#' van Genuchten Soil Water Retention Model
#'
#' This function computes the volumetric water content (`theta`) for a given matric potential (`h`)
#' using the van Genuchten equation. The equation models the soil water retention curve, describing
#' the relationship between soil water content and soil matric potential.
#'
#' @param h Numeric, the matric potential (pressure head) in cmH2O. Often negative for unsaturated soils
#'          (e.g., -33 cmH2O for field capacity, -1500 cmH2O for wilting point).
#' @param alpha Numeric, the van Genuchten parameter related to the inverse of the air-entry value.
#'              Units are in 1/cm.
#' @param n Numeric, the van Genuchten shape parameter related to pore size distribution. Must be greater than 1.
#' @param theta_r Numeric, the residual water content (minimum water content).
#' @param theta_s Numeric, the saturated water content (maximum water content).
#'
#' @return Numeric, the predicted volumetric water content (`theta`) at the given matric potential (`h`).
#'
#' @details The van Genuchten model equation is given by:
#'   \deqn{\theta(h) = \theta_r + \frac{(\theta_s - \theta_r)}{(1 + (|\alpha h|)^n)^{m}}}
#' where \eqn{m = 1 - \frac{1}{n}}. This equation models the water retention curve in soils,
#' describing how water content decreases with decreasing matric potential.
#'
#' @examples
#' # Example: Calculate water content at field capacity (-33 kPa) using van Genuchten parameters
#' h_fc <- -33 * 10.19716  # Field capacity matric potential in cmH2O (conversion from kPa)
#' alpha <- 0.08           # Example alpha (1/cm)
#' n <- 1.2                # Example n parameter (dimensionless)
#' theta_r <- 0.05         # Residual water content
#' theta_s <- 0.4          # Saturated water content
#' theta_fc <- van_genuchten(h_fc, alpha, n, theta_r, theta_s)
#' print(theta_fc)
#' @export
van_genuchten <- function(h, alpha, n, theta_r, theta_s) {
  m <- 1 - (1 / n)
  as.numeric(theta_r + (theta_s - theta_r) / ((1 + (abs(alpha * h))^n)^m))
}


#' Simulate Available Water Holding Capacity (AWHC) using van Genuchten Parameters from the ROSETTA Model
#'
#' This function simulates the Available Water Holding Capacity (AWHC) for soil components based on the variability
#' of van Genuchten parameters estimated from the ROSETTA model. The parameters for the van Genuchten model
#' (alpha, n, theta_r, theta_s) are sampled using Monte Carlo simulations based on the provided means and standard
#' deviations for each parameter. The function calculates water retention at specific matric potentials (Field Capacity
#' and Permanent Wilting Point) and computes AWHC as the difference between these two points.
#'
#' @param data Dataframe containing soil parameters for each component. It must include the following columns:
#'        - `alpha`: The van Genuchten alpha parameter (mean).
#'        - `sd_alpha`: The standard deviation for alpha.
#'        - `npar`: The van Genuchten n parameter (mean).
#'        - `sd_npar`: The standard deviation for n.
#'        - `theta_r`: The residual water content (mean).
#'        - `sd_theta_r`: The standard deviation for theta_r.
#'        - `theta_s`: The saturated water content (mean).
#'        - `sd_theta_s`: The standard deviation for theta_s.
#'        - `compname`: The component name.
#'        - `hzname`: The horizon name.
#'        - `layerID`: A unique identifier for the soil layer (used as part of the key for simulation results).
#' @param n_simulations Integer, number of Monte Carlo simulations to perform (default is 100).
#'
#' @return A list where each element contains simulated water retention parameters for a given soil component.
#'         Each entry is a dataframe with the van Genuchten parameters and the corresponding water content at
#'         Field Capacity (FC), Permanent Wilting Point (PWP), and Available Water Holding Capacity (AWHC).
#'
#' @import dplyr
#'
#' @examples
#' \dontrun{
#'   data <- data.frame(
#'     compname = c("Component1", "Component2"),
#'     hzname = c("A", "Bt"),
#'     layerID = c("Component1_A", "Component2_Bt"),
#'     alpha = c(-1.84, -1.73),
#'     sd_alpha = c(0.06, 0.06),
#'     npar = c(0.13, 0.11),
#'     sd_npar = c(0.01, 0.01),
#'     theta_r = c(0.05, 0.07),
#'     sd_theta_r = c(0.01, 0.015),
#'     theta_s = c(0.4, 0.45),
#'     sd_theta_s = c(0.02, 0.025)
#'   )
#'   results <- simulate_vg_aws(data, n_simulations = 1000)
#'   head(results[["Component1_A"]])
#' }
#' @export
simulate_vg_aws <- function(data, n_simulations = 100) {

  # Initialize result list to store the simulations for each component
  simulation_results <- list()

  # Iterate over each row in the data
  for (i in seq_len(nrow(data))) {

    # Extract van Genuchten parameters and their uncertainties for the current component
    alpha_mean    <- data$alpha[i]
    alpha_sd      <- data$sd_alpha[i]
    n_mean        <- data$npar[i]
    n_sd          <- data$sd_npar[i]
    theta_r_mean  <- data$theta_r[i]
    theta_r_sd    <- data$sd_theta_r[i]
    theta_s_mean  <- data$theta_s[i]
    theta_s_sd    <- data$sd_theta_s[i]

    # If any of the parameters are NA, skip this row
    if (is.na(alpha_mean) | is.na(n_mean) | is.na(theta_r_mean) | is.na(theta_s_mean)) {
      next
    }

    # Generate random samples for the van Genuchten parameters based on the provided means and standard deviations
    set.seed(123)  # Set seed for reproducibility
    alpha_samples   <- stats::rnorm(n_simulations, mean = alpha_mean, sd = alpha_sd)
    n_samples       <- stats::rnorm(n_simulations, mean = n_mean, sd = n_sd)
    theta_r_samples <- stats::rnorm(n_simulations, mean = theta_r_mean, sd = theta_r_sd)
    theta_s_samples <- stats::rnorm(n_simulations, mean = theta_s_mean, sd = theta_s_sd)

    # Convert matric potentials from kPa to cmH2O (1 kPa = 10.19716 cmH2O)
    h_fc  <- -33 * 10.19716    # Field capacity (FC) matric potential in cmH2O
    h_pwp <- -1500 * 10.19716  # Permanent wilting point (PWP) matric potential in cmH2O

    # Create a dataframe to hold the sampled parameters
    results <- data.frame(
      alpha   = 10^(alpha_samples),   # Convert alpha back from log scale
      n       = 10^(n_samples),       # Convert npar back from log scale
      theta_r = theta_r_samples,
      theta_s = theta_s_samples,
      sim_num = seq_len(n_simulations)
    )

    # Calculate water retention at FC and PWP for each set of sampled parameters using the van Genuchten function
    results <- results |>
      dplyr::rowwise() |>
      dplyr::mutate(
        theta_fc  = van_genuchten(h_fc, alpha, n, theta_r, theta_s),   # Water content at Field Capacity
        theta_pwp = van_genuchten(h_pwp, alpha, n, theta_r, theta_s)    # Water content at Permanent Wilting Point
      ) |>
      dplyr::ungroup()

    # Calculate Available Water Holding Capacity (AWHC) as the difference between FC and PWP
    results$AWHC <- results$theta_fc - results$theta_pwp

    # Store the results in the list with a key combining component name and layerID
    simulation_results[[paste(data$layerID[i], i, sep = "_")]] <- results
  }

  # Return the list of simulation results for each component
  return(simulation_results)
}


#' Calculate Available Water Storage Data Frame
#'
#' This function calculates available water storage (AWHC) in the top 100 cm of the soil profile
#' based on simulated soil properties. It first runs the ROSETTA model using soil properties,
#' creates a unique layer ID, then simulates water storage using van Genuchten parameters.
#' The function groups the data by depth intervals and computes average AWHC for standard soil depth intervals.
#'
#' @param sim_data_df A data frame containing simulated soil properties.
#'
#' @return A data frame of AWHC for standard soil depth intervals.
#'
#' @import dplyr
#' @import tidyr
#' @import aqp
#' @import soilDB
#'
#' @examples
#' \dontrun{
#'   sim_aws <- calculate_aws_df(simulated_soil_data)
#'   head(sim_aws)
#' }
#' @export
calculate_aws_df <- function(sim_data_df) {
  # Run Rosetta using soil properties
  variables <- c(
    "sand_total", "silt_total", "clay_total",
    "bulk_density_third_bar", "water_retention_third_bar", "water_retention_15_bar"
  )
  rosetta_data <- soilDB::ROSETTA(sim_data_df, vars = variables, v = "3", include.sd = TRUE)

  # Create layerID by combining compname and hzdept_r
  sim_data_df$layerID <- paste(sim_data_df$compname, sim_data_df$hzdept_r, sep = "_")
  rosetta_data$layerID <- sim_data_df$layerID

  # Simulate available water storage using van Genuchten parameters
  aws <- simulate_vg_aws(rosetta_data)

  # Use Map to add new columns to each data frame in the list
  aws <- Map(function(df, top_value, bottom_value, compname_value, cokey_value) {
    df$top <- top_value
    df$bottom <- bottom_value
    df$compname <- compname_value
    df$cokey <- cokey_value
    return(df)
  }, aws, sim_data_df$hzdept_r, sim_data_df$hzdepb_r, sim_data_df$compname, sim_data_df$cokey)

  # Combine the list of data frames into one
  sim_aws_df <- dplyr::bind_rows(aws)

  # Assign depth structure to the data frame using aqp
  aqp::depths(sim_aws_df) <- cokey ~ hzdept_r + hzdepb_r

  # Compute slab summaries over specified depth intervals (e.g., 0, 5, 15, 30, 60, 100 cm)
  # Using sim_aws_df as the input to the slab function
  prof.slab <- aqp::slab(sim_aws_df,
                         fm = cokey ~ AWHC,
                         slab.structure = c(0, 5, 15, 30, 60, 100),
                         slab.fun = aqp::mean_na)

  # Reshape the slab output to wide format
  prof.slab.1 <- prof.slab |>
    dplyr::select(-contributing_fraction) |>
    tidyr::pivot_wider(names_from = variable, values_from = value)

  # Clean up missing values if necessary (this line may be adjusted based on desired behavior)
  prof.slab.1[is.na(prof.slab.1)] <- NA

  return(prof.slab.1)
}


#' Query SSURGO Database for Soil Data by Mukey
#'
#' This function queries the NRCS SSURGO database for specific soil attributes based on the provided map unit key (mukey).
#'
#' @param mukeys A vector of string or numeric map unit keys (mukeys) representing the map units to retrieve data for.
#'
#' @return A data frame with soil horizon data for the given mukeys, including aggregated rock fragment data.
#'
#' @import dplyr
#' @import soilDB
#' @export
get_aws_data_by_mukey <- function(mukeys) {
  # Step 1: Format mukey(s) for SQL
  formatted_mukey <- paste0("(", paste0("'", mukeys, "'", collapse = ","), ")")

  # Step 2: Construct the SQL query
  query <- sprintf("
    SELECT
      -- contextual data
      co.mukey, co.cokey, compname, comppct_l, comppct_r, comppct_h,
      -- horizon morphology
      ch.chkey, hzname, hzdept_l, hzdept_r, hzdept_h, hzdepb_l, hzdepb_r, hzdepb_h, hzthk_l, hzthk_r, hzthk_h,
      -- soil texture parameters
      sandtotal_l, sandtotal_r, sandtotal_h,
      silttotal_l, silttotal_r, silttotal_h,
      claytotal_l, claytotal_r, claytotal_h,
      -- bulk density and water retention parameters
      dbovendry_l, dbovendry_r, dbovendry_h,
      wthirdbar_l, wthirdbar_r, wthirdbar_h,
      wfifteenbar_l, wfifteenbar_r, wfifteenbar_h,
      -- rock fragment volume from chfrags table
      chf.fragvol_l AS rfv_l, chf.fragvol_r AS rfv_r, chf.fragvol_h AS rfv_h
    FROM mapunit AS mu
      INNER JOIN component AS co ON mu.mukey = co.mukey
      INNER JOIN chorizon AS ch ON co.cokey = ch.cokey
      LEFT JOIN chfrags AS chf ON ch.chkey = chf.chkey
    WHERE mu.mukey IN %s
    ORDER BY co.cokey, ch.hzdept_r ASC;", formatted_mukey)

  # Step 3: Submit the query and fetch results
  result <- soilDB::SDA_query(query)

  # Step 4: Group by `chkey` and sum rock fragment values (rfv_l, rfv_r, rfv_h)
  rfv_sum <- result |>
    dplyr::group_by(chkey) |>
    dplyr::summarise(
      rfv_l = sum(rfv_l, na.rm = TRUE),
      rfv_r = sum(rfv_r, na.rm = TRUE),
      rfv_h = sum(rfv_h, na.rm = TRUE)
    ) |>
    dplyr::ungroup()

  # Step 5: Merge the summed values back to the main dataset and remove duplicates
  result <- result |>
    dplyr::select(-dplyr::all_of(c("rfv_l", "rfv_r", "rfv_h"))) |>
    dplyr::left_join(rfv_sum, by = "chkey") |>
    dplyr::distinct()

  # Step 6: Return the final result
  return(result)
}


#--------------------------------------------------------------------------------------------------------
# Functions to run soil profile depth simulations

#' Query OSD Data and Convert Horizon Distinctness to Offset
#'
#' This function retrieves Official Series Description (OSD) data for a given list of soil series (component names),
#' extracts horizon-related fields such as id, top, bottom, hzname, and distinctness, and converts the horizon
#' distinctness codes into offset values using the `hzDistinctnessCodeToOffset` function.
#'
#' @param horizon_data Data frame containing SSURGO horizon data.
#'
#' @return A data frame containing the following columns:
#'         - `id`: Unique identifier for the soil profile.
#'         - `hzname`: The name of the soil horizon.
#'         - `distinctness`: The distinctness code for the horizon (e.g., A, C, G, D).
#'         - `genhz`: The generalized horizon designation.
#'         - `bound_sd`: The calculated offset value from the distinctness code.
#'
#' @import dplyr
#' @import soilDB
#'
#' @examples
#' \dontrun{
#'   soils <- c('amador', 'pentz', 'pardee', 'auburn', 'loafercreek', 'millvilla')
#'   result <- query_osd_distinctness(soils)
#'   head(result)
#' }
#' @export
query_osd_distinctness <- function(horizon_data) {
  # Subset SSURGO horizon data to keep only compname and hzname
  ssurgo_horizon_data <- horizon_data |>
    dplyr::select(compname, hzname)

  # Ensure 'genhz' column exists in ssurgo_horizon_data
  if (!"genhz" %in% colnames(ssurgo_horizon_data)) {
    ssurgo_horizon_data$genhz <- generalizeHz(
      ssurgo_horizon_data$hzname,
      new = c('O', 'A', 'B', 'C', 'Cr', 'R'),
      pattern = c('O', '^A', '^B', '^C', '^Cr', '^R'),
      ordered = TRUE
    )
  }

  # Rename 'compname' to 'id'
  colnames(ssurgo_horizon_data)[colnames(ssurgo_horizon_data) == "compname"] <- "id"

  # Fetch OSD data for the list of soil component names using soilDB
  s <- soilDB::fetchOSD(unique(horizon_data$compname))

  # Check if the necessary fields are available in the OSD horizon data
  if ("distinctness" %in% names(s@horizons) &&
      "top" %in% names(s@horizons) &&
      "bottom" %in% names(s@horizons) &&
      "hzname" %in% names(s@horizons) &&
      "id" %in% names(s@horizons)) {

    # Extract relevant fields: id, hzname, distinctness
    osd_horizon_data <- s@horizons[, c("id", "hzname", "distinctness")]

    # Generate generalized horizon codes for the OSD data
    osd_horizon_data$genhz <- generalizeHz(
      osd_horizon_data$hzname,
      new = c('O', 'A', 'B', 'C', 'Cr', 'R'),
      pattern = c('O', '^\\d*A', '^\\d*B', '^\\d*C', '^\\d*Cr', '^\\d*R'),
      ordered = TRUE
    )

    # Identify missing genhz values in ssurgo_horizon_data that are not in osd_horizon_data
    missing_genhz <- setdiff(ssurgo_horizon_data$genhz, osd_horizon_data$genhz)

    # If there are missing genhz values, append them with NA for distinctness
    if (length(missing_genhz) > 0) {
      missing_rows <- ssurgo_horizon_data[ssurgo_horizon_data$genhz %in% missing_genhz, c("id", "hzname", "genhz")]
      missing_rows$distinctness <- NA
      osd_horizon_data <- base::rbind(osd_horizon_data, missing_rows)
      message("Added missing genhz values: ", paste(missing_genhz, collapse = ", "))
    } else {
      message("No missing genhz values found.")
    }

    # Infill missing distinctness values using a custom function
    osd_horizon_data <- infill_missing_distinctness(osd_horizon_data)

    # Convert the distinctness codes to offset values
    osd_horizon_data$bound_sd <- hzDistinctnessCodeToOffset(osd_horizon_data$distinctness)

    # Rename the columns for clarity
    names(osd_horizon_data) <- c("id", "hzname", "distinctness", "genhz", "bound_sd")

    return(osd_horizon_data)
  } else {
    stop("Necessary fields (distinctness, id, top, bottom, hzname) are not available in the OSD data.")
  }
}


#' Infill Missing Distinctness Values for Horizons
#'
#' This function fills in missing `distinctness` values for common horizon names (O, A, B, C, R, Cr)
#' based on typical boundary characteristics observed in the field. It assigns reasonable default
#' values for common horizon types.
#'
#' ## Default Boundary Distinctness Assignments:
#' - **O Horizons (Organic Layers):** diffuse
#' - **A Horizons (Topsoil):** clear
#' - **B Horizons (Subsoil):** gradual
#' - **C Horizons (Parent Material):** gradual
#' - **Cr Horizons (Weathered Bedrock):** gradual
#' - **R Horizons (Bedrock):** abrupt
#'
#' @param horizon_data A data frame containing horizon data. It must include `hzname` and `distinctness` columns.
#'
#' @return A data frame with missing `distinctness` values infilled based on the horizon name or generalized horizon group.
#'
#' @examples
#' \dontrun{
#'   # Example horizon data
#'   df <- data.frame(
#'     hzname = c("A", "B", "C", "R", "O", "Cr"),
#'     distinctness = c(NA, "gradual", NA, NA, "diffuse", NA),
#'     stringsAsFactors = FALSE
#'   )
#'   df <- infill_missing_distinctness(df)
#'   print(df)
#' }
#' @export
infill_missing_distinctness <- function(horizon_data) {
  # List of default boundary distinctness values for common horizons
  default_distinctness <- list(
    "O"  = "diffuse",  # Organic horizon
    "A"  = "clear",    # Topsoil
    "E"  = "clear",    # Subsoil
    "B"  = "gradual",  # Subsoil
    "C"  = "gradual",  # Parent material
    "Cr" = "gradual",  # Weathered bedrock
    "R"  = "abrupt"    # Bedrock
  )

  # Assign a generalized horizon group (genhz) if it's not already present
  if (!"genhz" %in% colnames(horizon_data)) {
    horizon_data$genhz <- aqp::generalizeHz(
      horizon_data$hzname,
      new = c('O', 'A', 'E', 'B', 'C', 'Cr', 'R'),
      pattern = c('O', '^\\d*A', '^\\d*E', '^\\d*B', '^\\d*C', '^\\d*Cr', '^\\d*R'),
      ordered = TRUE
    )
  }

  # Function to assign default distinctness based on horizon name or generalized horizon group
  fill_default_distinctness <- function(hzname, genhz, distinctness) {
    if (!is.na(distinctness)) {
      return(distinctness)  # If distinctness exists, return it
    }
    # Use default distinctness based on genhz first, then hzname if available
    if (!is.na(genhz) && genhz %in% names(default_distinctness)) {
      return(default_distinctness[[genhz]])
    } else if (!is.na(hzname) && hzname %in% names(default_distinctness)) {
      return(default_distinctness[[hzname]])
    } else {
      return(NA)  # Return NA if no match is found
    }
  }

  # Apply the default distinctness infilling using mapply
  horizon_data$distinctness <- mapply(
    fill_default_distinctness,
    horizon_data$hzname,
    horizon_data$genhz,
    horizon_data$distinctness
  )

  return(horizon_data)
}


#' Infill Missing Depth Variability
#'
#' This function infills missing depth variability values in a horizon data frame by replacing missing
#' top and bottom depth bounds with the representative depth value plus or minus 2, ensuring that the
#' resulting values are not negative.
#'
#' @param horizon_data A data frame containing horizon depth data. It must include the following columns:
#'   - `hzdept_r`: Representative top depth.
#'   - `hzdept_l`: Top low value (may be missing).
#'   - `hzdept_h`: Top high value (may be missing).
#'   - `hzdepb_r`: Representative bottom depth.
#'   - `hzdepb_l`: Bottom low value (may be missing).
#'   - `hzdepb_h`: Bottom high value (may be missing).
#'
#' @return A data frame with missing depth variability values filled in. Missing top low values are replaced with
#'   `hzdept_r - 2` and missing top high values with `hzdept_r + 2`. Similarly, missing bottom low values are replaced
#'   with `hzdepb_r - 2` and missing bottom high values with `hzdepb_r + 2`. All values are ensured to be non-negative.
#'
#' @examples
#' \dontrun{
#'   # Example usage:
#'   data <- data.frame(
#'     hzdept_r = c(10, 20),
#'     hzdept_l = c(NA, 18),
#'     hzdept_h = c(NA, NA),
#'     hzdepb_r = c(30, 40),
#'     hzdepb_l = c(NA, NA),
#'     hzdepb_h = c(NA, 42)
#'   )
#'   infill_missing_depth_variability(data)
#' }
#' @export
infill_missing_depth_variability <- function(horizon_data) {
  # Infill missing top low values (hzdept_l) with hzdept_r - 2, ensuring non-negative results
  horizon_data$hzdept_l <- ifelse(
    is.na(horizon_data$hzdept_l),
    pmax(horizon_data$hzdept_r - 2, 0),
    horizon_data$hzdept_l
  )

  # Infill missing top high values (hzdept_h) with hzdept_r + 2
  horizon_data$hzdept_h <- ifelse(
    is.na(horizon_data$hzdept_h),
    pmax(horizon_data$hzdept_r + 2, 0),
    horizon_data$hzdept_h
  )

  # Infill missing bottom low values (hzdepb_l) with hzdepb_r - 2
  horizon_data$hzdepb_l <- ifelse(
    is.na(horizon_data$hzdepb_l),
    pmax(horizon_data$hzdepb_r - 2, 0),
    horizon_data$hzdepb_l
  )

  # Infill missing bottom high values (hzdepb_h) with hzdepb_r + 2
  horizon_data$hzdepb_h <- ifelse(
    is.na(horizon_data$hzdepb_h),
    pmax(horizon_data$hzdepb_r + 2, 0),
    horizon_data$hzdepb_h
  )

  return(horizon_data)
}


#' Top-down Simulation of Soil Profile
#'
#' This function simulates a soil profile by calculating the top and bottom depths for each horizon
#' in a top-down manner. The first horizon always starts at 0 cm, and the top of each subsequent
#' horizon is set equal to the bottom depth of the previous horizon. The bottom depth of each horizon
#' is simulated using a triangular distribution via the `tri_dist` function. The triangular distribution
#' is defined by a lower bound (`hzdepb_l`), a representative value (`hzdepb_r`), and an upper bound (`hzdepb_h`).
#'
#' @param horizon_data A data frame containing horizon data. It should include the following columns:
#'   - `hzname`: Name of the horizon.
#'   - `hzdepb_l`: Lower bound for the bottom depth of the horizon.
#'   - `hzdepb_r`: Representative bottom depth of the horizon.
#'   - `hzdepb_h`: Upper bound for the bottom depth of the horizon.
#'
#' @return A data frame with the following columns:
#'   - `hzname`: The name of the horizon.
#'   - `top`: The simulated top depth of the horizon.
#'   - `bottom`: The simulated bottom depth of the horizon.
#'
#' @examples
#' \dontrun{
#'   # Example horizon data
#'   horizon_data <- data.frame(
#'     hzname = c("A", "B", "C"),
#'     hzdepb_l = c(20, 40, 60),
#'     hzdepb_r = c(25, 45, 65),
#'     hzdepb_h = c(30, 50, 70)
#'   )
#'   profile <- simulate_soil_profile_top_down(horizon_data)
#'   print(profile)
#' }
#' @export
simulate_soil_profile_top_down <- function(horizon_data) {
  n <- nrow(horizon_data)

  # Initialize a result data frame with hzname, top, and bottom columns
  result <- data.frame(
    hzname = horizon_data$hzname,
    top = rep(NA, n),
    bottom = rep(NA, n)
  )

  # Set the top of the first horizon to 0
  result$top[1] <- 0

  for (i in 1:n) {
    if (i > 1) {
      # Set the top of the current horizon to be the bottom of the previous horizon
      result$top[i] <- result$bottom[i - 1]
    }

    # Determine bounds for bottom depth simulation
    bottom_low  <- max(horizon_data$hzdepb_l[i], result$top[i])
    bottom_high <- max(bottom_low, horizon_data$hzdepb_h[i])
    bottom_mode <- min(max(horizon_data$hzdepb_r[i], bottom_low), bottom_high)

    # Generate a single random sample using the tri_dist function
    bottom_depth <- tri_dist(n = 1, a = bottom_low, b = bottom_high, c = bottom_mode)
    bottom_depth <- round(bottom_depth)

    # Ensure the bottom depth is not less than the top depth
    if (bottom_depth < result$top[i]) {
      bottom_depth <- result$top[i]
    }
    result$bottom[i] <- bottom_depth
  }

  return(result)
}


#' Bottom-up Simulation of Soil Profile
#'
#' This function simulates a soil profile from the bottom upward by calculating the top and bottom
#' depths for each horizon. The simulation starts from the bottom horizon and proceeds upward,
#' setting the bottom depth of each horizon to the top depth of the horizon immediately below.
#' For each horizon, the bottom and top depths are simulated using a triangular distribution via
#' the `tri_dist` function, with parameters defined by the lower, representative, and upper depth
#' values. The function ensures that there are no gaps between horizons and that the top of the
#' uppermost horizon is set to 0.
#'
#' @param horizon_data A data frame containing horizon data. The data frame must include the following columns:
#'   - `hzname`: The name of the horizon.
#'   - `hzdept_l`: The lower bound for the top depth of the horizon.
#'   - `hzdept_r`: The representative top depth of the horizon.
#'   - `hzdept_h`: The upper bound for the top depth of the horizon.
#'   - `hzdepb_l`: The lower bound for the bottom depth of the horizon.
#'   - `hzdepb_r`: The representative bottom depth of the horizon.
#'   - `hzdepb_h`: The upper bound for the bottom depth of the horizon.
#'
#' @return A data frame with the following columns:
#'   - `hzname`: The name of the horizon.
#'   - `top`: The simulated top depth of the horizon.
#'   - `bottom`: The simulated bottom depth of the horizon.
#'
#' @examples
#' \dontrun{
#'   # Example horizon data
#'   horizon_data <- data.frame(
#'     hzname = c("A", "B", "C"),
#'     hzdept_l = c(0, 20, 40),
#'     hzdept_r = c(0, 25, 45),
#'     hzdept_h = c(0, 30, 50),
#'     hzdepb_l = c(20, 40, 60),
#'     hzdepb_r = c(25, 45, 65),
#'     hzdepb_h = c(30, 50, 70)
#'   )
#'   profile <- simulate_soil_profile_bottom_up(horizon_data)
#'   print(profile)
#' }
#' @export
simulate_soil_profile_bottom_up <- function(horizon_data) {
  n <- nrow(horizon_data)

  # Initialize a result data frame with hzname, top, and bottom columns
  result <- data.frame(
    hzname = horizon_data$hzname,
    top = rep(NA, n),
    bottom = rep(NA, n)
  )

  # Process horizons from bottom to top
  for (i in n:1) {
    # Simulate bottom depth using tri_dist
    bottom_low  <- horizon_data$hzdepb_l[i]
    bottom_high <- horizon_data$hzdepb_h[i]
    bottom_mode <- horizon_data$hzdepb_r[i]

    # Generate a single random sample using tri_dist
    bottom_depth <- tri_dist(n = 1, a = bottom_low, b = bottom_high, c = bottom_mode)
    bottom_depth <- round(bottom_depth)

    if (i == n) {
      # For the bottom horizon, simulate top depth
      top_low  <- horizon_data$hzdept_l[i]
      top_high <- min(horizon_data$hzdept_h[i], bottom_depth)
      top_mode <- min(max(horizon_data$hzdept_r[i], top_low), top_high)

      # Generate a single random sample using tri_dist
      top_depth <- tri_dist(n = 1, a = top_low, b = top_high, c = top_mode)
      top_depth <- round(top_depth)

      # Ensure top depth is not greater than bottom depth
      if (top_depth > bottom_depth) {
        top_depth <- bottom_depth
      }
    } else {
      # For other horizons, set bottom depth to the top depth of the horizon below
      bottom_depth <- result$top[i + 1]

      # Simulate top depth using tri_dist
      top_low  <- horizon_data$hzdept_l[i]
      top_high <- min(horizon_data$hzdept_h[i], bottom_depth)
      top_mode <- min(max(horizon_data$hzdept_r[i], top_low), top_high)

      # Generate a single random sample using tri_dist
      top_depth <- tri_dist(n = 1, a = top_low, b = top_high, c = top_mode)
      top_depth <- round(top_depth)

      # Ensure top depth is not greater than bottom depth
      if (top_depth > bottom_depth) {
        top_depth <- bottom_depth
      }
    }

    # Assign the simulated depths to the result
    result$top[i]    <- top_depth
    result$bottom[i] <- bottom_depth
  }

  # Ensure the top depth of the uppermost horizon is 0
  result$top[1] <- 0

  # Adjust the bottom depth of the first horizon if necessary
  if (result$bottom[1] < result$top[1]) {
    result$bottom[1] <- result$top[1]
  }

  # Ensure no gaps between horizons
  for (i in 1:(n - 1)) {
    if (result$bottom[i] != result$top[i + 1]) {
      result$bottom[i] <- result$top[i + 1]
    }
  }

  return(result)
}


#' Simulate Soil Profile Thickness
#'
#' This function simulates the thickness of soil horizons for a given soil profile using both top-down and bottom-up
#' simulation methods. It performs multiple simulations and calculates the standard deviation of horizon thickness
#' (i.e., the difference between representative bottom and top depths) across the simulations. The summarized results
#' include the horizon name, representative top and bottom depths, and the standard deviation of thickness for each horizon.
#'
#' @param horizon_data A dataframe containing horizon data. It must include the following columns:
#'        - `hzname`: The name of each horizon (e.g., A, B, C).
#'        - `hzdept_r`: The representative top depth of each horizon.
#'        - `hzdepb_r`: The representative bottom depth of each horizon.
#'        - `hzdept_l`: The low estimate for the top depth of each horizon.
#'        - `hzdept_h`: The high estimate for the top depth of each horizon.
#'        - `hzdepb_l`: The low estimate for the bottom depth of each horizon.
#'        - `hzdepb_h`: The high estimate for the bottom depth of each horizon.
#' @param n_simulations Integer, number of simulations to perform (default is 500).
#'
#' @return A dataframe containing the following columns:
#'         - `hzname`: The horizon name.
#'         - `top`: The representative top depth of the horizon (from `hzdept_r`).
#'         - `bottom`: The representative bottom depth of the horizon (from `hzdepb_r`).
#'         - `thickness_sd`: The standard deviation of the horizon thickness (bottom - top) across the simulations.
#'
#' @examples
#' \dontrun{
#'   # Example horizon data with representative and bound estimates for depths
#'   horizon_data <- data.frame(
#'     hzname = c("A", "B", "C"),
#'     hzdept_r = c(0, 20, 35),
#'     hzdepb_r = c(20, 35, 50),
#'     hzdept_l = c(0, 15, 30),
#'     hzdept_h = c(0, 25, 40),
#'     hzdepb_l = c(15, 30, 45),
#'     hzdepb_h = c(25, 40, 60)
#'   )
#'   set.seed(123)
#'   summarized_results <- simulate_soil_profile_thickness(horizon_data, n_simulations = 500)
#'   print(summarized_results)
#' }
#'
#' @import dplyr
#' @export
simulate_soil_profile_thickness <- function(horizon_data, n_simulations = 500) {
  # Step 1: Infill missing depth variability values
  horizon_data <- infill_missing_depth_variability(horizon_data)

  results_list <- list()

  # Loop over the number of simulations
  for (sim in 1:n_simulations) {
    # Randomly select the simulation direction: top-down or bottom-up
    if (stats::runif(1) < 0.5) {
      # Simulate profile using the top-down method
      result <- simulate_soil_profile_top_down(horizon_data)
      method <- "Top-Down"
    } else {
      # Simulate profile using the bottom-up method
      result <- simulate_soil_profile_bottom_up(horizon_data)
      method <- "Bottom-Up"
    }

    # Calculate thickness for each horizon as the difference between bottom and top depths
    result$Thickness <- result$bottom - result$top

    # Record the simulation method and simulation number
    result$Method <- method
    result$Simulation <- sim

    # Store the current simulation result
    results_list[[sim]] <- result
  }

  # Combine all simulation results into a single dataframe
  all_results <- do.call(rbind, results_list)

  # Reorder columns for clarity: Simulation, Method, hzname, top, bottom, Thickness
  all_results <- all_results[, c("Simulation", "Method", "hzname", "top", "bottom", "Thickness")]

  # Group results by horizon name and calculate the standard deviation of thickness for each horizon
  horizon_sd <- all_results %>%
    dplyr::group_by(hzname) %>%
    dplyr::summarise(thickness_sd = round(sd(Thickness), 2)) %>%
    dplyr::ungroup()

  # Merge the representative top and bottom depths from horizon_data
  summarized_results <- merge(
    horizon_data[, c("hzname", "hzdept_r", "hzdepb_r")],
    horizon_sd,
    by = "hzname"
  )

  # Rename columns for consistency: hzdept_r -> top, hzdepb_r -> bottom
  summarized_results <- summarized_results %>%
    dplyr::rename(
      top = hzdept_r,
      bottom = hzdepb_r
    )

  # Return the summarized results
  return(summarized_results)
}


#' Simulate and Perturb Soil Profiles
#'
#' This function takes a SoilProfileCollection object and perturbs its horizons to simulate variability
#' in soil profile thickness and boundary depths. The simulation is performed in two steps:
#' first, horizon thickness is perturbed based on simulated thickness variability (using a top-down or bottom-up
#' simulation method), and then the boundary depths are perturbed using distinctness-derived offsets.
#' If the profile contains only one horizon, the function simply replicates the profile.
#'
#' @param soil_profile A SoilProfileCollection object containing soil profile data.
#'
#' @return A SoilProfileCollection object with perturbed horizon depths.
#'
#' @import dplyr
#' @import aqp
#' @export
simulate_and_perturb_soil_profiles <- function(soil_profile) {
  # Step 1: Extract horizons from the soil profile and select relevant columns
  soil_profile@horizons <- aqp::horizons(soil_profile) %>%
    dplyr::select(mukey, id, cokey, compname, hzname, sim_comppct,
                  hzdept_l, hzdept_r, hzdept_h, hzdepb_l, hzdepb_r, hzdepb_h, hzthk_l, sim_comppct)

  horizon_data <- aqp::horizons(soil_profile)
  horizon_data$genhz <- aqp::generalizeHz(
    horizon_data$hzname,
    new = c('O','A','B','C','Cr','R'),
    pattern = c('O', '^\\d*A','^\\d*B','^\\d*C','^\\d*Cr','^\\d*R'),
    ordered = TRUE
  )
  n_simulations <- unique(horizon_data$sim_comppct)  # Set number of simulations based on sim_comppct

  # Check if there's only one horizon (e.g., R horizon only)
  if (nrow(horizon_data) == 1) {
    message("Only one horizon in the profile. Skipping perturbation for this profile.")

    # Return the original soil profile replicated n_simulations times (no perturbation needed)
    replicated_profiles <- list()
    for (i in 1:n_simulations) {
      # Create a copy of the soil profile and assign a unique profile ID
      new_profile <- soil_profile
      aqp::profile_id(new_profile) <- paste0(aqp::profile_id(new_profile), "_sim_", i)
      replicated_profiles[[i]] <- new_profile
    }
    simulated_profiles <- aqp::combine(replicated_profiles)
    return(simulated_profiles)
  }

  # Step 1 (continued): Simulate soil profile thickness
  simulated_thickness <- simulate_soil_profile_thickness(horizon_data, n_simulations)
  simulated_thickness$genhz <- aqp::generalizeHz(
    simulated_thickness$hzname,
    new = c('O','A','B','C','Cr','R'),
    pattern = c('O', '^\\d*A','^\\d*B','^\\d*C','^\\d*Cr','^\\d*R'),
    ordered = TRUE
  )

  # Step 2: Query OSD distinctness and get bound_sd values
  distinctness_data <- query_osd_distinctness(horizon_data)
  bound_lut <- distinctness_data %>%
    dplyr::group_by(id, genhz) %>%
    dplyr::summarise(bound_sd = mean(bound_sd, na.rm = TRUE)) %>%
    dplyr::ungroup()

  # Step 3: Merge the simulated thickness data and distinctness data
  combined_data <- merge(simulated_thickness, bound_lut, by = "genhz") %>%
    dplyr::arrange(top)

  # Step 4: Add thickness standard deviation and bound_sd to the soil profile horizons
  horizon_data <- horizon_data %>%
    dplyr::left_join(combined_data %>% dplyr::select(hzname, thickness_sd, bound_sd), by = "hzname")
  soil_profile@horizons <- horizon_data %>%
    dplyr::select(mukey, id, cokey, compname, hzname, sim_comppct, hzdept_r, hzdepb_r, thickness_sd, bound_sd)

  # Step 3: First perturbation (horizon thickness)
  perturbed_profiles_thickness <- aqp::perturb(
    soil_profile,
    n = n_simulations,
    thickness.attr = "thickness_sd"
  )

  # Step 4: Second perturbation (boundary distinctness)
  list_of_profiles <- list()
  for (i in seq_along(perturbed_profiles_thickness)) {
    single_profile <- perturbed_profiles_thickness[i]
    perturbed_profile <- aqp::perturb(
      single_profile,
      n = 1,
      boundary.attr = "bound_sd",
      id = single_profile@site$id
    )
    list_of_profiles[[i]] <- perturbed_profile
  }

  simulated_profiles <- aqp::combine(list_of_profiles)

  # Step 5: Identify out-of-range horizons and adjust depths
  adjust_out_of_range_profiles <- function(simulated_profiles, horizon_data) {
    sim_horizons <- aqp::horizons(simulated_profiles) %>%
      dplyr::select(id, mukey, cokey, compname, hzname, top = hzdept_r, bottom = hzdepb_r)

    merged_data <- merge(
      sim_horizons,
      horizon_data[, c("hzname", "hzdept_l", "hzdepb_h")],
      by = "hzname",
      all.x = TRUE
    )

    merged_data <- merged_data %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        top = dplyr::if_else(is.na(hzdept_l), top, pmax(top, hzdept_l)),
        bottom = dplyr::if_else(is.na(hzdepb_h), bottom, pmin(bottom, hzdepb_h))
      ) %>%
      dplyr::ungroup()

    simulated_profiles_hzdata <- aqp::horizons(simulated_profiles) %>%
      dplyr::left_join(merged_data %>% dplyr::select(hzname, id, top, bottom), by = c("id", "hzname"))
    simulated_profiles@horizons <- simulated_profiles_hzdata

    return(simulated_profiles)
  }

  simulated_profiles <- adjust_out_of_range_profiles(simulated_profiles, horizon_data)

  # Step 6: Return the adjusted simulated profiles
  return(simulated_profiles)
}


#' Run Soil Profile Depth Simulations for a SoilProfileCollection
#'
#' This function simulates soil profile depths for each profile in a SoilProfileCollection by applying
#' the `simulate_and_perturb_soil_profiles` function to each individual profile. The simulated profiles
#' from all individual profiles are then combined into a single SoilProfileCollection. A seed can be set for reproducibility.
#'
#' @param soil_collection A SoilProfileCollection object containing soil profile data.
#' @param seed An integer value for setting the random seed (default is 123) for reproducible simulations.
#'
#' @return A SoilProfileCollection object containing the simulated soil profiles.
#'
#' @examples
#' \dontrun{
#'   # Assuming 'soil_collection' is a valid SoilProfileCollection object
#'   simulated_collection <- simulate_profile_depths_by_collection(soil_collection, seed = 123)
#' }
#'
#' @import aqp
#' @export
simulate_profile_depths_by_collection <- function(soil_collection, seed = 123) {

  # Check if input is a valid SoilProfileCollection
  if (!inherits(soil_collection, "SoilProfileCollection")) {
    stop("Input must be a SoilProfileCollection.")
  }

  # Set seed for reproducibility
  set.seed(seed)

  # Initialize an empty list to store all simulated profiles
  all_simulated_profiles <- list()

  # Loop over each profile in the SoilProfileCollection
  for (i in seq_along(soil_collection)) {
    # Extract the current soil profile
    soil_profile <- soil_collection[i, ]

    # Simulate and perturb the current soil profile
    simulated_profiles <- simulate_and_perturb_soil_profiles(soil_profile)

    # Append the simulated profiles to the list
    all_simulated_profiles[[i]] <- simulated_profiles
  }

  # Combine all simulated profiles into a single SoilProfileCollection
  combined_simulated_profiles <- aqp::combine(all_simulated_profiles)

  # Return the combined simulated profiles
  return(combined_simulated_profiles)
}


#' Run Soil Profile Depth Simulations in Parallel for a SoilProfileCollection
#'
#' This function simulates soil profile depths for each profile in a SoilProfileCollection in parallel
#' using the 'future.apply' package. It sets up a parallel backend with the specified number of cores,
#' processes each soil profile individually using the `simulate_and_perturb_soil_profiles` function,
#' and combines the results into a single SoilProfileCollection. A seed is used for reproducibility.
#'
#' @param soil_collection A SoilProfileCollection object containing soil profile data.
#' @param seed An integer value to set the random seed for reproducibility (default is 123).
#' @param n_cores Integer, number of cores to use for parallel processing (default is 6).
#'
#' @return A SoilProfileCollection object containing the simulated and perturbed soil profiles.
#'
#' @examples
#' \dontrun{
#'   simulated_profiles <- simulate_profile_depths_by_collection_parallel(soil_collection, seed = 123, n_cores = 6)
#'   print(simulated_profiles)
#' }
#'
#' @import future
#' @import future.apply
#' @import aqp
#' @export
simulate_profile_depths_by_collection_parallel <- function(soil_collection, seed = 123, n_cores = 6) {

  # Set up the parallel backend using the future package
  future::plan(future::multisession, workers = n_cores)

  # Check if input is a valid SoilProfileCollection
  if (!inherits(soil_collection, "SoilProfileCollection")) {
    stop("Input must be a SoilProfileCollection.")
  }

  # Set seed for reproducibility
  set.seed(seed)

  # Define a function to simulate a single soil profile
  simulate_single_profile <- function(i) {
    soil_profile <- soil_collection[i, ]  # Select the ith profile
    cat("Processing profile ID:", soil_profile$id, "\n")
    simulate_and_perturb_soil_profiles(soil_profile)
  }

  # Run the simulation in parallel using future_lapply
  result <- tryCatch({
    all_simulated_profiles <- future.apply::future_lapply(seq_along(soil_collection), simulate_single_profile)

    # Revert to sequential processing
    future::plan(future::sequential)

    # Combine all simulated profiles into a single SoilProfileCollection
    combined_simulated_profiles <- aqp::combine(all_simulated_profiles)

    return(combined_simulated_profiles)
  }, error = function(e) {
    cat("Error in processing profiles:\n")
    cat("Error message:", e$message, "\n")
    # Optionally return NULL to indicate failure
    return(NULL)
  })

  return(result)
}


#' Simulate Profile Depths by Mukey
#'
#' This function queries soil data for a given map unit key (mukey) using the SSURGO database,
#' sets up the corresponding soil profile data, and then simulates and perturbs the soil profiles
#' for that mukey. It performs a specified number of simulations per profile (using the
#' `simulate_and_perturb_soil_profiles` function) and combines the results into a single
#' SoilProfileCollection. A seed can be set for reproducibility.
#'
#' @param mukey A string or numeric value representing the map unit key to query.
#' @param n_simulations Integer, the number of simulations to perform for each profile (default is 100).
#' @param seed An integer to set the random seed for reproducibility (default is 123).
#'
#' @return A SoilProfileCollection object containing the simulated and perturbed soil profiles.
#'
#' @examples
#' \dontrun{
#'   # Query and simulate profiles for a specific mukey
#'   simulated_profiles <- simulate_profile_depths_by_mukey("123456", n_simulations = 100, seed = 123)
#'   print(simulated_profiles)
#' }
#'
#' @import aqp
#' @export
simulate_profile_depths_by_mukey <- function(mukey, n_simulations = 100, seed = 123) {
  # Step 1: Query the data based on the mukey
  mu_data <- get_aws_data_by_mukey(mukeys = mukey)

  if (is.null(mu_data)) {
    stop("No data found for the provided mukey.")
  }

  # Step 2: Set up the soil profile data
  set.seed(seed)  # For reproducibility
  aqp::depths(mu_data) <- compname ~ hzdept_r + hzdepb_r
  aqp::hzdesgnname(mu_data) <- "hzname"

  # Initialize an empty list to store all simulated profiles
  all_simulated_profiles <- list()

  # Step 3: Loop over each profile in mu_data
  for (i in seq_along(mu_data)) {
    # Select the current profile
    soil_profile <- mu_data[i, ]

    # Step 4: Simulate n_simulations for the current soil profile
    simulated_profiles <- simulate_and_perturb_soil_profiles(soil_profile, n_simulations = n_simulations)

    # Append the simulated profiles to the list
    all_simulated_profiles[[i]] <- simulated_profiles
  }

  # Step 5: Combine all simulated profiles into a single SoilProfileCollection
  combined_simulated_profiles <- aqp::combine(all_simulated_profiles)

  # Step 6: Return the combined simulated profiles
  return(combined_simulated_profiles)
}


#' Evaluate Simulated Profile Depths
#'
#' This function evaluates whether the simulated soil profile depths fall outside the specified range.
#' It extracts the simulated top and bottom depths from a SoilProfileCollection, merges these with the
#' original horizon data, and flags horizons where the simulated top is less than the lower bound or the
#' simulated bottom exceeds the upper bound.
#'
#' @param simulated_profiles A SoilProfileCollection object containing simulated soil profiles.
#' @param horizon_data A data frame with original horizon data, including columns:
#'        - `hzname`: Horizon name.
#'        - `hzdept_l`: Lower bound for the top depth.
#'        - `hzdepb_h`: Upper bound for the bottom depth.
#'
#' @return A data frame containing rows (horizons) where the simulated depths are out of range.
#'
#' @examples
#' \dontrun{
#'   out_of_range <- evaluate_simulated_depths(simulated_profiles, horizon_data)
#'   head(out_of_range)
#' }
#'
#' @import dplyr
#' @import aqp
#' @export
evaluate_simulated_depths <- function(simulated_profiles, horizon_data) {
  # Extract horizons from the simulated profiles
  sim_horizons <- aqp::horizons(simulated_profiles) %>%
    dplyr::select(mukey, cokey, compname, hzname, top = hzdept_r, bottom = hzdepb_r)

  # Merge simulated horizons with original horizon data for comparison
  merged_data <- merge(
    sim_horizons,
    horizon_data[, c("hzname", "hzdept_l", "hzdepb_h")],
    by = "hzname",
    all.x = TRUE
  )

  # Check if simulated top and bottom depths fall outside the specified range
  merged_data$out_of_range <- with(merged_data, (top < hzdept_l) | (bottom > hzdepb_h))

  # Return rows where simulated depths are out of range
  out_of_range_data <- merged_data[merged_data$out_of_range == TRUE, ]

  return(out_of_range_data)
}


#' aws_mukey_quant <- function(sim_aws_df) {
#'   #' Need to wrap in function in form that can be applied to each mukey.
#'   #' @return A list with available water storage prediction interval (5,50,5,PIW90) for mukey.
#'   aws_grouped <- sim_aws_df %>% group_by(top)
#'   data_len_depth <- summarise(aws_grouped, depth_len = n())
#'   sim_aws_df <- left_join(sim_aws_df, data_len_depth, by = "top")
#'
#'   # Step 1b: Reshape aws data by mukey
#'   aws_grouped_bottom <- sim_aws_df %>% group_by(bottom)
#'   aws_quant_list <- summarise(aws_grouped_bottom, aws_quant = quantile(sim_aws_df, probs = c(0.05, 0.50, 0.95)))
#'
#'   # Step 1c: Rename and group data
#'   aws05 <- calculate_aws_quant(aws_quant_list, "0.05")
#'   aws95 <- calculate_aws_quant(aws_quant_list, "0.95")
#'   aws50 <- calculate_aws_quant(aws_quant_list, "0.50")
#'
#'   # Width of 90th prediction interval for available water storage (aws_PIW90)
#'   aws_PIW90 <- round(aws95$aws0.95_100 - aws05$aws0.05_100, 2)
#'
#'   return(c(aws_PIW90, aws05, aws95, aws50))
#'
#'
#' # Calculate Available Water Storage (AWS) for Region of Interest (ROI)
#' calculate_aws_quant <- function(df, quantile) {
#'   total <- sum(df[[quantile]] * df$depth * df$n)
#'   return(data.frame(!!paste0("aws", quantile, "_100") := total))
#' }
