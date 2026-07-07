#' @title Percentile-Based Distribution Simulation Helpers
#'
#' @description Canonical `generate_inverse_cdf_distribution()` and
#' `calculate_summary_statistics()` functions, extracted during the code/
#' reorganization from code/notebooks/soilgrids_simulation.Rmd (the cleanest of
#' several parallel "simulate from percentiles" implementations found across the
#' project, and the most defensively-coded of three `calculate_summary_statistics()`
#' redefinitions). Other, non-canonical "simulate from percentiles" approaches
#' (a monotonic-spline inverse-CDF in code/scratch/spline_inv_cdf.R, and
#' density/KDE- and metalog-based approaches in code/notebooks/) are left in place
#' as parallel exploratory methods rather than merged in, since choosing among them
#' is a scientific-methods decision, not a reorg decision.

#' Build a piecewise-linear inverse-CDF sampler from percentile columns
#'
#' Given a one-row data frame of percentile values (columns named like "P0", "P5",
#' "P50", "P95", "P100"), builds a piecewise-linear inverse-CDF via `approxfun()`
#' and draws `n` samples from it via inverse-transform sampling.
#'
#' @param quantile_df A data frame with at least one row and percentile-named columns.
#' @param percentile_cols Character vector of candidate percentile column names to use.
#' @param n Number of samples to draw.
#' @return A numeric vector of `n` simulated values.
#' @export
generate_inverse_cdf_distribution <- function(quantile_df,
                                              percentile_cols = c("P0","P5","P50","P95","P100"),
                                              n = 1000) {
  # 1) Verify percentile columns are present and valid
  present_cols <- percentile_cols[percentile_cols %in% names(quantile_df)]
  if (length(present_cols) == 0) {
    stop("No valid percentile columns found in the data frame.")
  }

  # Remove any columns with NA or -1 in the first row (adapt to your own rules)
  valid_cols <- present_cols[
    !is.na(quantile_df[1, present_cols]) &
      quantile_df[1, present_cols] != -1
  ]
  if (length(valid_cols) < 2) {
    stop("Need at least 2 valid quantile columns to construct inverse CDF.")
  }

  # 2) Extract the numeric values (the actual data) for each valid percentile col
  x_vals <- unlist(quantile_df[1, valid_cols], use.names = FALSE)

  # 3) Convert column names like "P0","P5","P50","P95","P100" into probabilities 0, 0.05, 0.5, 0.95, 1.0
  #    This step assumes that your column names always start with 'P' and the rest is the number.
  #    If your naming is different, adapt accordingly.
  get_prob <- function(colname) {
    # remove 'P' and convert to numeric, e.g. "5" => 5. Then divide by 100 => 0.05
    prob <- as.numeric(sub("P", "", colname)) / 100
    return(prob)
  }
  cdf_vals <- sapply(valid_cols, get_prob)

  # 4) Sort them in ascending order, so we have a proper CDF
  #    (In case the user gave them out of order)
  sort_idx <- order(cdf_vals)
  cdf_vals <- cdf_vals[sort_idx]
  x_vals   <- x_vals[sort_idx]

  # 5) Build a piecewise-linear inverse CDF function using approxfun
  #    We'll invert cdf_vals -> x_vals
  inv_cdf <- approxfun(
    x = cdf_vals,
    y = x_vals,
    method = "linear",
    rule   = 2   # rule=2: out-of-range values are clamped to the boundary
  )

  # 6) Sample from uniform(0,1) and map through the inverse CDF
  u <- runif(n, min = 0, max = 1)
  samples <- inv_cdf(u)

  return(samples)
}

#' Compute summary statistics (mean, SD, CV, percentiles, quartiles) for a numeric vector
#'
#' @param data A numeric vector.
#' @param percentile_probs Numeric vector of percentile probabilities (0-1) to compute.
#' @return A one-row data frame of summary statistics, including dynamically-named
#'   percentile columns (e.g. P10, P20, ...).
#' @export
calculate_summary_statistics <- function(data, percentile_probs = seq(0.1, 0.9, by = 0.1)) {
    # Ensure the input is a numeric vector
    if (!is.numeric(data)) {
        stop("Data must be a numeric vector.")
    }
    if (!is.numeric(percentile_probs) || any(percentile_probs < 0 | percentile_probs > 1)) {
        stop("Percentile probabilities must be numeric values between 0 and 1.")
    }

    # Calculating basic statistics
    num <- length(data)
    mean_value <- mean(data)
    std_dev <- sd(data)
    cv <- (std_dev / mean_value) * 100
    median_value <- median(data)
    mad_value <- mad(data)
    min_value <- min(data)
    max_value <- max(data)
    variance_value <- var(data)
    se <- std_dev / sqrt(num)

    # Calculating percentiles
    percentiles <- quantile(data, probs = percentile_probs, na.rm = TRUE)

    # Calculating quartiles
    quartiles <- quantile(data, probs = c(0.25, 0.75), na.rm = TRUE)
    quart1 <- quartiles[1]
    quart3 <- quartiles[2]

    # Creating a named vector of percentiles (e.g., P05, P50, P95)
    percentile_names <- paste0("P", sprintf("%02d", as.integer(percentile_probs * 100)))
    percentile_values <- as.list(percentiles)
    names(percentile_values) <- percentile_names

    # Combining results into a single data frame
    summary_df <- data.frame(
        Num = num,
        Mean = mean_value,
        STD = std_dev,
        CV = cv,
        Median = median_value,
        MAD = mad_value,
        Min = min_value,
        Max = max_value,
        Var = variance_value,
        Quart1 = quart1,
        Quart3 = quart3,
        SE = se
    )

    # Append percentile columns dynamically
    summary_df <- cbind(summary_df, as.data.frame(t(percentile_values)))

    return(summary_df)
}
