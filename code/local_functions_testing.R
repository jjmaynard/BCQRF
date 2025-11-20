






# Define the function
get_aws_data_by_mukey <- function(mukeys) {
  #' Query SSURGO database for soil data by mukey
  #'
  #' This function queries the NRCS SSURGO database for specific soil attributes based on the provided map unit key (mukey).
  #'
  #' @param mukeys A vector of string or numeric map unit keys (mukeys) representing the map units to retrieve data for.
  #'
  #' @return A data frame with soil horizon data for the given mukeys, including aggregated rock fragment data.

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
    -- tables of interest
    FROM mapunit AS mu
    -- implicit filtering
    INNER JOIN component AS co ON mu.mukey = co.mukey
    INNER JOIN chorizon AS ch ON co.cokey = ch.cokey
    -- join with chfrags table using the correct chkey
    LEFT JOIN chfrags AS chf ON ch.chkey = chf.chkey
    -- filter by provided mukey
    WHERE mu.mukey IN %s
    -- order for visual inspection
    ORDER BY co.cokey, ch.hzdept_r ASC;", formatted_mukey)

  # Step 3: Submit the query and fetch results
  result <- SDA_query(query)

  # Step 4: Group by `chkey` and sum rock fragment values (rfv_l, rfv_r, rfv_h)
  rfv_sum <- result %>%
    dplyr::group_by(chkey) %>%
    summarise(
      rfv_l = sum(rfv_l, na.rm = TRUE),
      rfv_r = sum(rfv_r, na.rm = TRUE),
      rfv_h = sum(rfv_h, na.rm = TRUE)
    ) %>%
    ungroup()

  # Step 5: Merge the summed values back to the main dataset and remove duplicates
  result <- result %>%
    dplyr::select(-c(rfv_l, rfv_r, rfv_h)) %>%  # Remove original rfv columns
    left_join(rfv_sum, by = "chkey") %>%  # Merge summed values
    distinct()  # Remove any duplicates

  # Step 6: Return the final result
  return(result)
}

# Query OSD Data and Convert Horizon Distinctness to Offset
query_osd_distinctness <- function(soils) {
  #'
  #' This function retrieves Official Series Description (OSD) data for a given list of soil series (component names),
  #' extracts horizon-related fields such as id, top, bottom, hzname, and distinctness, and converts the horizon
  #' distinctness codes into offset values using the `hzDistinctnessCodeToOffset` function.
  #'
  #' The function returns a data frame containing the extracted horizon data along with the distinctness code and the
  #' corresponding offset value.
  #'
  #' @param soils Character vector of soil series (component names) for which to fetch OSD data. These names
  #'              correspond to the Official Series Descriptions in the OSD database.
  #'
  #' @return A data frame containing the following columns:
  #'         - `id`: Unique identifier for the soil profile.
  #'         - `top`: The top depth of the horizon (in cm).
  #'         - `bottom`: The bottom depth of the horizon (in cm).
  #'         - `hzname`: The name of the soil horizon.
  #'         - `distinctness`: The distinctness code for the horizon (e.g., A, C, G, D).
  #'         - `bound_sd`: The calculated SD value from the distinctness code.
  #'
  #' @examples
  #' # Example soil series to query
  #' soils <- c('amador', 'pentz', 'pardee', 'auburn', 'loafercreek', 'millvilla')
  #'
  #' # Query OSD data and convert distinctness codes to offsets
  #' result <- query_osd_distinctness(soils)
  #'
  #' # View the result
  #' head(result)
  #'
  #' @export

  # Fetch OSD data for the list of soil component names
  s <- fetchOSD(soils)

  # Check if the necessary fields are available in the horizon data
  if ("distinctness" %in% names(s@horizons) &&
      "top" %in% names(s@horizons) &&
      "bottom" %in% names(s@horizons) &&
      "hzname" %in% names(s@horizons) &&
      "id" %in% names(s@horizons)) {

    # Extract relevant fields: id, top, bottom, hzname, distinctness
    horizon_data <- s@horizons[, c("id", "top", "bottom", "hzname", "distinctness")]
    horizon_data$genhz <- generalizeHz(
      horizon_data$hzname,
      new = c('O','A','B','C', 'Cr','R'),
      pattern = c('O', '^\\d*A','^\\d*B','^\\d*C', '^\\d*Cr', '^\\d*R'),
      ordered = TRUE
    )
    horizon_data <- infill_missing_distinctness(horizon_data)

    # Convert the distinctness codes to offset values using hzDistinctnessCodeToOffset
    horizon_data$bound_sd <- hzDistinctnessCodeToOffset(horizon_data$distinctness)

    # Rename the columns for clarity
    names(horizon_data) <- c("id", "top", "bottom", "hzname", "distinctness", "genhz", "bound_sd")

    return(horizon_data)
  } else {
    stop("Necessary fields (distinctness, id, top, bottom, hzname) are not available in the OSD data.")
  }
}


infill_missing_distinctness <- function(horizon_data) {
  #' Infill Missing Distinctness Values for Horizons
  #'
  #' This function fills in missing `distinctness` values for common horizon names (O, A, B, C, R, Cr)
  #' based on typical boundary characteristics observed in the field. It assigns reasonable default
  #' values for common horizon types.
  #'
  #' ## Default Boundary Distinctness Assignments:
  #' - **O Horizons (Organic Layers):** Diffuse
  #' - **A Horizons (Topsoil):** Clear
  #' - **B Horizons (Subsoil):** Gradual
  #' - **C Horizons (Parent Material):** Gradual
  #' - **Cr Horizons (Weathered Bedrock):** Gradual
  #' - **R Horizons (Bedrock):** Abrupt
  #'
  #' @param horizon_data A data frame containing horizon data. It must include `hzname` and `distinctness` columns.
  #'
  #' @return A data frame with missing `distinctness` values infilled based on the horizon name or generalized horizon group.

  # List of default boundary distinctness values for common horizons
  default_distinctness <- list(
    "O" = "Diffuse",      # Organic horizon
    "A" = "Clear",        # Topsoil
    "B" = "Gradual",      # Subsoil
    "C" = "Gradual",      # Parent material
    "Cr" = "Gradual",     # Weathered bedrock
    "R" = "Abrupt"        # Bedrock
  )

  # Assign a generalized horizon group (genhz) if it's not already present
  if (!"genhz" %in% colnames(horizon_data)) {
    horizon_data$genhz <- generalizeHz(
      horizon_data$hzname,
      new = c('O','A','B','C','Cr','R'),
      pattern = c('O', '^A','^B','^C','^Cr', '^R'),
      ordered = TRUE
    )
  }

  # Function to assign default distinctness based on horizon name or generalized horizon group
  fill_default_distinctness <- function(hzname, genhz, distinctness) {
    if (!is.na(distinctness)) {
      return(distinctness)  # If distinctness exists, return it
    }

    # Use default distinctness based on hzname or genhz
    if (!is.na(genhz) && genhz %in% names(default_distinctness)) {
      return(default_distinctness[[genhz]])
    } else if (!is.na(hzname) && hzname %in% names(default_distinctness)) {
      return(default_distinctness[[hzname]])
    } else {
      return(NA)  # Return NA if no match is found
    }
  }

  # Apply the default distinctness infilling
  horizon_data$distinctness <- mapply(
    fill_default_distinctness,
    horizon_data$hzname,
    horizon_data$genhz,
    horizon_data$distinctness
  )

  return(horizon_data)
}


# Triangular Distribution Sampling Function, code from the 'triangle' R package
tri_dist <- function(n = 1, a = 0, b = 1, c = (a + b)/2) {
  #' Triangular Distribution Sampling
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


# Function to infill missing values with representative values ± 2, ensuring no negative values
infill_missing_depth_variability <- function(horizon_data) {
  # Infill missing top low values (hzdept_l) with hzdept_r - 2
  horizon_data$hzdept_l <- ifelse(is.na(horizon_data$hzdept_l),
                                  pmax(horizon_data$hzdept_r - 2, 0),  # Ensure no negative values
                                  horizon_data$hzdept_l)

  # Infill missing top high values (hzdept_h) with hzdept_r + 2
  horizon_data$hzdept_h <- ifelse(is.na(horizon_data$hzdept_h),
                                  pmax(horizon_data$hzdept_r + 2, 0),
                                  horizon_data$hzdept_h)

  # Infill missing bottom low values (hzdepb_l) with hzdepb_r - 2
  horizon_data$hzdepb_l <- ifelse(is.na(horizon_data$hzdepb_l),
                                  pmax(horizon_data$hzdepb_r - 2, 0),
                                  horizon_data$hzdepb_l)

  # Infill missing bottom high values (hzdepb_h) with hzdepb_r + 2
  horizon_data$hzdepb_h <- ifelse(is.na(horizon_data$hzdepb_h),
                                  pmax(horizon_data$hzdepb_r + 2, 0),
                                  horizon_data$hzdepb_h)

  return(horizon_data)
}


# Top-down simulation function
simulate_soil_profile_top_down <- function(horizon_data) {
  n <- nrow(horizon_data)

  # Initialize a result data frame with hzname, top, and bottom columns
  result <- data.frame(
    hzname = horizon_data$hzname,
    top = rep(NA, n),
    bottom = rep(NA, n)
  )

  # Initialize the top of the first horizon
  result$top[1] <- 0

  for (i in 1:n) {
    if (i > 1) {
      # Set the top of the current horizon equal to the bottom of the previous horizon
      result$top[i] <- result$bottom[i - 1]
    }

    # Simulate bottom depth using tri_dist instead of simulate_triangular
    bottom_low <- max(horizon_data$hzdepb_l[i], result$top[i])
    bottom_high <- max(bottom_low, horizon_data$hzdepb_h[i])
    bottom_mode <- min(max(horizon_data$hzdepb_r[i], bottom_low), bottom_high)

    # Generate a single random sample using tri_dist
    bottom_depth <- tri_dist(n = 1, a = bottom_low, b = bottom_high, c = bottom_mode)
    bottom_depth <- round(bottom_depth)

    # Ensure bottom depth is at least the top depth
    if (bottom_depth < result$top[i]) {
      bottom_depth <- result$top[i]
    }
    result$bottom[i] <- bottom_depth
  }

  return(result)
}



# Bottom-up simulation function
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
    bottom_low <- horizon_data$hzdepb_l[i]
    bottom_high <- horizon_data$hzdepb_h[i]
    bottom_mode <- horizon_data$hzdepb_r[i]

    # Generate a single random sample using tri_dist
    bottom_depth <- tri_dist(n = 1, a = bottom_low, b = bottom_high, c = bottom_mode)
    bottom_depth <- round(bottom_depth)

    # For the bottom horizon, simulate top depth
    if (i == n) {
      top_low <- horizon_data$hzdept_l[i]
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
      # For other horizons, set bottom depth to top depth of the horizon below
      bottom_depth <- result$top[i + 1]

      # Simulate top depth using tri_dist
      top_low <- horizon_data$hzdept_l[i]
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
    result$top[i] <- top_depth
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
#' This function simulates the thickness of soil horizons for a given soil profile, using top-down and bottom-up
#' methods, and calculates the standard deviation of horizon thickness across multiple simulations. The function
#' returns summarized results that include the horizon name, representative top and bottom depths, and the standard
#' deviation of thickness for each horizon.
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
#'         - `sd_thickness`: The standard deviation of the horizon thickness across the simulations.
#'
#' @examples
#' # Example horizon data (with representative top and bottom depths)
#' horizon_data <- data.frame(
#'   hzname = c("A", "B", "C"),
#'   hzdept_r = c(0, 20, 35),  # Representative top depths
#'   hzdepb_r = c(20, 35, 50),  # Representative bottom depths
#'   hzdept_l = c(0, 15, 30),   # Low top depths
#'   hzdept_h = c(0, 25, 40),   # High top depths
#'   hzdepb_l = c(15, 30, 45),  # Low bottom depths
#'   hzdepb_h = c(25, 40, 60)   # High bottom depths
#' )
#'
#' # Simulate 500 soil profiles with thickness calculations
#' set.seed(123)
#' n_simulations <- 500
#' summarized_results <- simulate_soil_profile_thickness(horizon_data, n_simulations)
#'
#' # View the summarized results
#' print(summarized_results)
#'
#' @export
simulate_soil_profile_thickness <- function(horizon_data, n_simulations = 500) {
  # Step 1: Infill missing values
  horizon_data <- infill_missing_depth_variability(horizon_data)

  results_list <- list()

  # Loop over the number of simulations
  for (sim in 1:n_simulations) {
    # Randomly select the simulation direction
    if (runif(1) < 0.5) {
      # Simulate top-down direction
      result <- simulate_soil_profile_top_down(horizon_data)
      method <- "Top-Down"
    } else {
      # Simulate bottom-up direction
      result <- simulate_soil_profile_bottom_up(horizon_data)
      method <- "Bottom-Up"
    }

    # Calculate the thickness for each horizon as the difference between bottom and top depths
    result$Thickness <- result$bottom - result$top

    # Add the method used and the simulation number to the result
    result$Method <- method
    result$Simulation <- sim

    # Store the result of the current simulation in the results list
    results_list[[sim]] <- result
  }

  # Combine all simulation results into a single data frame
  all_results <- do.call(rbind, results_list)

  # Reorder columns for clarity: Simulation, Method, hzname, top, bottom, Thickness
  all_results <- all_results[, c("Simulation", "Method", "hzname", "top", "bottom", "Thickness")]

  # Group the results by hzname and calculate the standard deviation of thickness for each horizon
  horizon_sd <- all_results %>%
    group_by(hzname) %>%
    summarise(
      thickness_sd = sd(Thickness) %>% round(2)         # Calculate the standard deviation of thickness
    )

  # Merge the representative top and bottom depths from horizon_data
  summarized_results <- merge(horizon_data[, c("hzname", "hzdept_r", "hzdepb_r")], horizon_sd,
                              by.x = "hzname", by.y = "hzname")

  # Rename hzdept_r and hzdepb_r to top and bottom for consistency
  summarized_results <- summarized_results %>%
    rename(
      top = hzdept_r,
      bottom = hzdepb_r
    )

  # Return the summarized results
  return(summarized_results)
}


simulate_and_perturb_soil_profiles <- function(soil_profile, n_simulations = 25) {

  # Step 0: format data
  # Soil profile data
  horizon_data <- soil_profile@horizons
  # OSD soil component names
  osd_soils <- unique(horizons(soil_profile)$compname)

  # Step 1: Simulate soil profile thickness
  simulated_thickness <- simulate_soil_profile_thickness(horizon_data, n_simulations)
  simulated_thickness$genhz <- generalizeHz(
    horizon_data$hzname,
    new = c('O','A','B','C', 'Cr','R'),
    pattern = c('O', '^\\d*A','^\\d*B','^\\d*C', '^\\d*Cr', '^\\d*R'),
    ordered = TRUE
  )

  # Step 2: Query OSD distinctness and get bound_sd values
  distinctness_data <- query_osd_distinctness(osd_soils)
  bound_lut <- distinctness_data %>%
    group_by(id, genhz) %>%
    summarise(bound_sd = mean(bound_sd, na.rm = TRUE)) %>%  # Calculate the mean of 'bound_sd'
    ungroup()  # Remove grouping

  # Step 3: Merge the simulated thickness data and distinctness data
  combined_data <- merge(simulated_thickness, bound_lut, by = "genhz")

  # Step 4: Add sd_thickness (thickness.attr) and bound_sd (boundary.attr) to the soil profile collection
  horizons(soil_profile)$thickness_sd <- combined_data$thickness_sd
  horizons(soil_profile)$bound_sd <- combined_data$bound_sd

  # Step 3: First perturbation (horizon thickness)
  perturbed_profiles_thickness <- aqp::perturb(
    soil_profile,
    n = n_simulations,
    thickness.attr = "thickness_sd", max.depth=250
  )
  horizons(perturbed_profiles_thickness)$bound_sd <- horizons(soil_profile)$bound_sd

  # Step 4: Second perturbation (boundary distinctness)
  # Initialize an empty list to store profiles
  list_of_profiles <- list()

  for (i in 1:length(perturbed_profiles_thickness)) {
    # Extract individual profile
    single_profile <- perturbed_profiles_thickness[i]

    # Perturb boundary depths
    perturbed_profile <- aqp::perturb(
      single_profile,
      n = 1,
      boundary.attr = "bound_sd", max.depth=250
    )

    # Store the perturbed profile
    list_of_profiles[[i]] <- perturbed_profile
  }

  # Combine all perturbed profiles into a SoilProfileCollection
  simulated_profiles <- aqp::combine(list_of_profiles)

  # Step 5: Identify out-of-range horizons and adjust depths
  adjust_out_of_range_profiles <- function(simulated_profiles, horizon_data) {
    sim_horizons <- horizons(simulated_profiles) %>% dplyr::select(mukey,cokey,compname, hzname, top = hzdept_r, bottom = hzdepb_r)

    # Merge simulated horizons with the original horizon data for comparison
    merged_data <- merge(
      sim_horizons,
      horizon_data[, c("hzname", "hzdept_l", "hzdepb_h")],
      by = "hzname",
      all.x = TRUE
    )

    # Adjust depths if they are out of range
    merged_data$top <- pmax(merged_data$top, merged_data$hzdept_l)  # Adjust if top is too shallow
    merged_data$bottom <- pmin(merged_data$bottom, merged_data$hzdepb_h)  # Adjust if bottom is too deep

    # Update the adjusted data back to the profile
    horizons(simulated_profiles)$top <- merged_data$top
    horizons(simulated_profiles)$bottom <- merged_data$bottom
    return(simulated_profiles)
  }

  # Apply the adjustment for out-of-range horizons
  simulated_profiles_adjusted <- adjust_out_of_range_profiles(simulated_profiles, horizon_data)

  # Step 6: Return the adjusted simulated profiles
  return(simulated_profiles_adjusted)

}

simulate_profile_depths_by_mukey <- function(mukey, n_simulations = 100, seed = 123) {

  # Step 1: Query the data based on the mukey
  mu_data <- get_aws_data_by_mukey(mukeys = mukey)

  if (is.null(mu_data)) {
    stop("No data found for the provided mukey.")
  }

  # Step 2: Set up the soil profile data
  set.seed(seed)  # For reproducibility
  depths(mu_data) <- compname ~ hzdept_r + hzdepb_r
  hzdesgnname(mu_data) <- "hzname"

  # Initialize an empty list to store all simulated profiles
  all_simulated_profiles <- list()

  # Step 3: Loop over each profile in mu_data
  for (i in seq_along(mu_data)) {
    # Select each profile
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

simulate_profile_depths_by_collection <- function(soil_collection, n_simulations = 100, seed = 123) {
  # Load necessary libraries
  library(aqp)

  # Check if input is a valid SoilProfileCollection
  if (!inherits(soil_collection, "SoilProfileCollection")) {
    stop("Input must be a SoilProfileCollection.")
  }

  # Step 1: Set up soil profile data for reproducibility
  set.seed(seed)

  # Initialize an empty list to store all simulated profiles
  all_simulated_profiles <- list()

  # Step 2: Loop over each profile in the SoilProfileCollection
  for (i in seq_along(soil_collection)) {
    # Select each profile
    soil_profile <- soil_collection[i, ]

    # Step 3: Simulate n_simulations for the current soil profile
    simulated_profiles <- simulate_and_perturb_soil_profiles(soil_profile, n_simulations = n_simulations)

    # Append the simulated profiles to the list
    all_simulated_profiles[[i]] <- simulated_profiles
  }

  # Step 4: Combine all simulated profiles into a single SoilProfileCollection
  combined_simulated_profiles <- aqp::combine(all_simulated_profiles)

  # Step 5: Return the combined simulated profiles
  return(combined_simulated_profiles)
}
#----------------------------------------------------------------------------------------------------------------------
# Example Dataset

mukey <- 545835
combined_profiles <- simulate_profile_depths_by_mukey(mukey, n_simulations = 20)

# tighter margins
par(mar = c(0, 0, 0, 1))
# profile sketches with customization of style
plotSPC(combined_profiles, print.id = FALSE, max.depth = 250, cex.names = 0.66, width = 0.35, name.style = 'center-center', depth.axis = list(style = 'compact', line = -2), color='claytotal_r')

#---------------------------------------------------------------------------------------------------------------------------



#Additional code for testing
# mukeys for Coweeta
mukeys <- c(
  545800, 545801, 545803, 545805, 545806, 545807, 545811, 545813, 545814, 545815,
  545818, 545819, 545820, 545829, 545830, 545831, 545835, 545836, 545837, 545838,
  545842, 545843, 545853, 545854, 545855, 545857, 545859, 545860, 545861, 545863,
  545874, 545875, 545876, 545878, 545882, 545885, 545886, 545887
)


ssurgo_data <- get_aws_data_by_mukey(mukeys)
# Simulate 500 soil profiles with thickness calculations
set.seed(123)
n_simulations <- 20

# Load sample data and convert it into a SoilProfileCollection
ssurgo_data_sub <- ssurgo_data %>% dplyr::filter(mukey == 545838)
depths(ssurgo_data_sub) <- compname ~ hzdept_r + hzdepb_r

# Select a profile to use as the basis for simulation
s <- ssurgo_data_sub[1,]
hzdesgnname(s) <- "hzname"

# Example horizon data and OSD soils
horizon_data <- s@horizons

# Example OSD soil component names
osd_soils <- unique(horizons(s)$compname)

# Step 6: Simulate 25 new profiles
simulated_profiles <- simulate_and_perturb_soil_profiles(s, horizon_data, osd_soils, n_simulations = 100)

# View the simulated profiles
plot(simulated_profiles)


perturbed_profiles_thickness <- aqp::perturb(
  soil_profile,
  n = n_simulations,
  thickness.attr = "thickness_sd"
)

perturbed_profiles_boundary <- aqp::perturb(
  soil_profile,
  n = n_simulations,
  boundary.attr = "bound_sd"
)

# tighter margins
par(mar = c(0, 0, 0, 1))
# profile sketches with customization of style
plotSPC(simulated_profiles, print.id = FALSE, max.depth = 250, cex.names = 0.66, width = 0.35, name.style = 'center-center', depth.axis = list(style = 'compact', line = -2), color='claytotal_r')

# profile sketches with customization of style
plotSPC(perturbed_profiles_thickness, print.id = FALSE, max.depth = 250, cex.names = 0.66, width = 0.35, name.style = 'center-center', depth.axis = list(style = 'compact', line = -2), color='claytotal_r')

# profile sketches with customization of style
plotSPC(perturbed_profiles_boundary, print.id = FALSE, max.depth = 250, cex.names = 0.66, width = 0.35, name.style = 'center-center', depth.axis = list(style = 'compact', line = -2), color='claytotal_r')


range(simulated_profiles@horizons$hzdepb_r)
range(perturbed_profiles_boundary@horizons$hzdepb_r)
range(perturbed_profiles_thickness@horizons$hzdepb_r)




# Function to evaluate if the simulated profile depths fall outside the specified range
evaluate_simulated_depths <- function(simulated_profiles, horizon_data) {

  # Extract horizons from simulated profiles
  sim_horizons <- horizons(simulated_profiles) %>% dplyr::select(mukey,cokey,compname, hzname, top = hzdept_r, bottom = hzdepb_r)

  # Merge simulated horizons with the original horizon data for comparison
  merged_data <- merge(
    sim_horizons,
    horizon_data[, c("hzname", "hzdept_l", "hzdepb_h")],
    by = "hzname",
    all.x = TRUE
  )

  # Check if simulated top and bottom depths fall outside the specified range
  merged_data$out_of_range <- with(merged_data, {
    (top < hzdept_l) | (bottom > hzdepb_h)
  })
  simulated_profiles$out_of_range == merged_data$out_of_range
# Return rows where simulated depths are out of range
out_of_range_data <- merged_data[merged_data$out_of_range == TRUE, ]

return(out_of_range_data)
}


# Assuming you have simulated profiles and horizon data
# Apply the evaluation function
out_of_range_horizons <- evaluate_simulated_depths(simulated_profiles, horizon_data)

# View the results
print(out_of_range_horizons)


sim1 <- aqp::perturb(soil_profile, n = 20, boundary.attr = "bound_sd")
sim2 <- aqp::perturb(soil_profile, n = 20, thickness.attr = "thickness_sd")
sim3 <- aqp::perturb(soil_profile, n = 20, boundary.attr = "thickness_sd")
sim4 <- aqp::perturb(soil_profile, n = 20, thickness.attr = "bound_sd")



