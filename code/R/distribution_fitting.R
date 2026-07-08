#' @title Percentile-Based Distribution Simulation
#'
#' @description Consolidates several "simulate a distribution from summary
#' percentiles" approaches that previously existed as separate, one-off
#' functions scattered across code/notebooks/ and code/scratch/ (a
#' piecewise-linear inverse-CDF, a monotonic-spline inverse-CDF, a KDE/
#' truncated-normal approach, a metalog fit, and beta/normal parametric
#' fits). These are genuinely different statistical approaches (different
#' smoothness, tail-extrapolation, and distributional assumptions), so
#' rather than picking one as "correct", `simulate_from_percentiles()`
#' exposes all of them behind a single `method` argument with a shared
#' input-validation front end. `compare_percentile_methods()` and
#' `validate_percentile_methods_synthetic()` (below) support comparing
#' methods against each other and, where a known ground-truth distribution
#' is available, against reconstruction accuracy.
#'
#' `generate_inverse_cdf_distribution()` is kept as a thin backwards-compatible
#' wrapper around `simulate_from_percentiles(..., method = "linear_cdf")`
#' since it is called directly from code/notebooks/soilgrids_simulation.Rmd.

#' Parse percentile columns (e.g. "P0","P5","P50") from a one-row data frame
#'
#' Shared validation/extraction front end used by every method in
#' `simulate_from_percentiles()`. Column names must be of the form "P<num>"
#' where <num> is the percentile (0-100). Columns that are NA or -1 (the
#' sentinel used elsewhere in this project for "not available") are dropped.
#'
#' @param quantile_df A data frame with at least one row and percentile-named columns.
#' @param percentile_cols Character vector of candidate percentile column names to use.
#' @return A list with `probs` (sorted probabilities, 0-1) and `values` (matching values).
#' @keywords internal
extract_percentile_pairs <- function(quantile_df, percentile_cols) {
  present_cols <- percentile_cols[percentile_cols %in% names(quantile_df)]
  if (length(present_cols) == 0) {
    stop("No valid percentile columns found in the data frame.")
  }

  valid_cols <- present_cols[
    !is.na(quantile_df[1, present_cols]) &
      quantile_df[1, present_cols] != -1
  ]
  if (length(valid_cols) < 2) {
    stop("Need at least 2 valid quantile columns.")
  }

  values <- unlist(quantile_df[1, valid_cols], use.names = FALSE)
  probs <- as.numeric(sub("P", "", valid_cols)) / 100

  sort_idx <- order(probs)
  list(probs = probs[sort_idx], values = values[sort_idx])
}

#' Piecewise-linear inverse-CDF sampler (nonparametric, exact at knots, clamped tails)
#' @keywords internal
sim_linear_cdf <- function(probs, values, n) {
  inv_cdf <- approxfun(x = probs, y = values, method = "linear", rule = 2)
  inv_cdf(runif(n, min = 0, max = 1))
}

#' Monotonic-spline inverse-CDF sampler (smoother than linear, still exact at knots)
#'
#' @param bounds Optional length-2 vector giving the value at prob 0 and prob 1.
#'   If NULL, extends the observed value range by 5% on each side.
#' @keywords internal
sim_spline <- function(probs, values, n, bounds = NULL) {
  if (is.null(bounds)) {
    pad <- 0.05 * diff(range(values))
    if (pad == 0) pad <- 1
    bounds <- c(min(values) - pad, max(values) + pad)
  }
  # Only pad in the 0/1 endpoints if they aren't already supplied (e.g. via
  # "P0"/"P100" columns) - otherwise splinefun() sees duplicate x values.
  has_p0 <- probs[1] == 0
  has_p1 <- probs[length(probs)] == 1
  ext_probs <- c(if (!has_p0) 0, probs, if (!has_p1) 1)
  ext_values <- c(if (!has_p0) bounds[1], values, if (!has_p1) bounds[2])
  inv_cdf <- splinefun(ext_probs, ext_values, method = "monoH.FC")
  inv_cdf(runif(n, min = 0, max = 1))
}

#' KDE sampler: seeds truncated-normal draws around each percentile, then
#' resamples from their kernel density estimate
#'
#' Generalizes the position-based (not column-name-based) logic previously
#' duplicated across code/notebooks/generating_empirical_parametric_distributions.Rmd
#' and code/notebooks/hwsd_wise_metalog_dist.Rmd.
#'
#' @param sample_size Number of truncated-normal draws used to build the KDE
#'   (higher = smoother density estimate, at more compute cost).
#' @keywords internal
sim_kde <- function(probs, values, n, sample_size = 10000) {
  if (!requireNamespace("truncnorm", quietly = TRUE)) {
    stop("Package 'truncnorm' is required for method = 'kde'.")
  }
  k <- length(values)
  dense_points <- numeric(0)
  for (i in seq_len(k)) {
    sd_i <- if (i == 1 || i == k) {
      abs(values[i] - values[min(k, i + 1)]) / 3
    } else {
      min(abs(values[i] - values[i - 1]), abs(values[i] - values[i + 1])) / 2
    }
    sd_i <- max(sd_i, 0.1)
    lower <- if (i == 1) values[1] else -Inf
    upper <- if (i == k) values[k] else Inf
    dense_points <- c(
      dense_points,
      truncnorm::rtruncnorm(ceiling(sample_size / k), a = lower, b = upper, mean = values[i], sd = sd_i)
    )
  }

  d <- density(dense_points, adjust = 1, kernel = "gaussian", n = sample_size * 10)
  in_range <- d$x >= min(values) & d$x <= max(values)
  d$x <- d$x[in_range]
  d$y <- d$y[in_range]
  d$y <- d$y / sum(d$y)

  sample(d$x, n, replace = TRUE, prob = d$y)
}

#' Metalog distribution sampler (fits a flexible, monotonic quantile function)
#'
#' Uses the public rmetalog package API (`metalog()` / `rmetalog()`) rather
#' than a hand-rolled reimplementation of its internals.
#'
#' @section KNOWN RELIABILITY ISSUES: `rmetalog::metalog()`'s internal
#'   fitting can be numerically unstable. Separately, this project's own
#'   prior experience (see code/scratch/metalog_soil_analysis_clean.R)
#'   found it can *hang* on high-variability-with-small-sample or
#'   high-skew data, which is why this function wraps each fit attempt in
#'   `setTimeLimit()` plus a fallback ladder that retries with
#'   progressively fewer terms.
#'
#'   Separately, while developing this function (R 4.5.2, rmetalog 1.0.3),
#'   `rmetalog::metalog()` calls made via `Rscript -e '<one-liner>'`
#'   reproducibly segfaulted (9/9 attempts, on the package's own textbook
#'   example and other inputs) - but the identical calls made from a
#'   sourced `.R` file, both called directly and through this function,
#'   succeeded 18/18 times. So the segfault looks like an artifact of that
#'   specific invocation style rather than a defect triggered by normal
#'   usage (sourced scripts, knitting an Rmd, running the `targets`
#'   pipeline) - but it was only characterized in this one environment. A
#'   segfault kills the R process outright and cannot be caught by
#'   `tryCatch()`, so if you see one, no amount of defensive code here can
#'   fully protect against it. Because the full risk profile isn't
#'   characterized across every calling context, "metalog" is left out of
#'   the default `methods` list in `compare_percentile_methods()` and
#'   `validate_percentile_methods_synthetic()` below - opt in explicitly,
#'   and if you ever see a crash, test in an isolated subprocess (e.g.
#'   `callr::r(function() rmetalog::metalog(...))`) before trusting it in
#'   a long-running pipeline.
#'
#' @param boundedness One of "u" (unbounded), "sl" (semi-bounded below),
#'   "su" (semi-bounded above), "b" (bounded). See `?rmetalog::metalog`.
#' @param bounds Length-2 vector of (lower, upper) bounds. Required unless
#'   boundedness = "u".
#' @param term Number of metalog terms to try first (higher = more flexible
#'   fit, but more prone to overfitting sparse percentiles and more prone to
#'   the instability described above). Capped at the number of valid
#'   percentiles and at 9. If fitting fails or times out, retries with
#'   progressively fewer terms down to 2 before giving up.
#' @param time_limit_sec Per-attempt wall-clock timeout in seconds, passed to
#'   `setTimeLimit()`. Does not protect against a segfault (see above).
#' @keywords internal
sim_metalog <- function(probs, values, n, bounds = NULL, boundedness = "u", term = 9, time_limit_sec = 10) {
  if (!requireNamespace("rmetalog", quietly = TRUE)) {
    stop("Package 'rmetalog' is required for method = 'metalog'.")
  }
  # p = 0 or 1 is mathematically invalid for metalog's logit-based fit
  # (log(p / (1 - p)) is undefined at the boundary) - drop them rather than
  # let rmetalog fail on them unpredictably.
  interior <- probs > 0 & probs < 1
  probs <- probs[interior]
  values <- values[interior]
  term_limit <- min(length(values), 9, term)
  if (term_limit < 2) {
    stop("Need at least 2 valid percentiles strictly between 0 and 1 to fit a metalog distribution.")
  }
  if (boundedness != "u" && is.null(bounds)) {
    stop("bounds = c(lower, upper) is required when boundedness != 'u'.")
  }
  bounds <- if (is.null(bounds)) c(0, 1) else bounds

  # Fallback ladder: try the requested term count, then progressively fewer
  # terms, since a lower-order fit is less prone to the numerical
  # instability described above. Each attempt is wrapped in setTimeLimit()
  # so a hang is aborted rather than left to run indefinitely.
  fit <- NULL
  for (k in seq(term_limit, 2, by = -1)) {
    fit <- tryCatch({
      setTimeLimit(cpu = time_limit_sec, elapsed = time_limit_sec, transient = TRUE)
      on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE), add = TRUE)
      rmetalog::metalog(
        values,
        bounds = bounds,
        boundedness = boundedness,
        term_limit = k,
        probs = probs
      )
    }, error = function(e) {
      warning(sprintf("metalog fit with %d terms failed/timed out: %s", k, conditionMessage(e)))
      NULL
    })
    if (!is.null(fit)) {
      term_limit <- k
      break
    }
  }
  if (is.null(fit)) {
    stop("Metalog fitting failed for all term counts from ", term, " down to 2.")
  }

  rmetalog::rmetalog(fit, n = n, term = term_limit)
}

#' Beta distribution sampler (fits a Beta distribution to percentiles rescaled to [0,1])
#' @param bounds Length-2 vector of (lower, upper) support bounds. Required.
#' @keywords internal
sim_beta <- function(probs, values, n, bounds) {
  if (!requireNamespace("fitdistrplus", quietly = TRUE)) {
    stop("Package 'fitdistrplus' is required for method = 'beta'.")
  }
  if (missing(bounds) || is.null(bounds)) {
    stop("bounds = c(lower, upper) is required for method = 'beta'.")
  }
  if (length(values) < 3) {
    stop("At least 3 percentiles are required to robustly fit the Beta distribution.")
  }

  values_scaled <- (values - bounds[1]) / (bounds[2] - bounds[1])
  fit <- tryCatch(
    fitdistrplus::fitdist(values_scaled, "beta"),
    error = function(e) stop("Error fitting Beta distribution: ", conditionMessage(e))
  )
  samples <- rbeta(n, shape1 = fit$estimate[1], shape2 = fit$estimate[2])
  samples * (bounds[2] - bounds[1]) + bounds[1]
}

#' Normal distribution sampler (P50 as mean; SD from the widest available percentile pair)
#' @keywords internal
sim_normal <- function(probs, values, n) {
  if (length(values) < 3) {
    stop("At least 3 percentiles are required to robustly fit the normal distribution.")
  }
  mean_idx <- which(probs == 0.5)
  if (length(mean_idx) == 0) {
    stop("Missing median (P50) value for normal distribution fitting.")
  }
  mean_val <- values[mean_idx]

  # Use the widest available (non-0/1) probability pair to estimate SD, generalizing
  # the original hardcoded P5/P95 (whose z-score span is exactly qnorm(0.95)-qnorm(0.05)).
  tail_idx <- probs > 0 & probs < 1 & probs != 0.5
  if (sum(tail_idx) < 2) {
    stop("Need at least two non-median percentiles to estimate the standard deviation.")
  }
  p_lo <- min(probs[tail_idx]); p_hi <- max(probs[tail_idx])
  x_lo <- values[probs == p_lo][1]; x_hi <- values[probs == p_hi][1]
  std_val <- (x_hi - x_lo) / (qnorm(p_hi) - qnorm(p_lo))

  rnorm(n, mean = mean_val, sd = std_val)
}

#' Simulate values from a distribution reconstructed from summary percentiles
#'
#' Dispatches to one of several methods for turning a handful of known
#' percentiles into simulated draws. See the file-level `@description` above
#' for how the methods differ and when to prefer one over another.
#'
#' @param quantile_df A data frame with at least one row and percentile-named
#'   columns (e.g. "P0", "P5", "P50", "P95", "P100").
#' @param method One of "linear_cdf" (default), "spline", "kde", "metalog",
#'   "beta", "normal".
#' @param percentile_cols Character vector of candidate percentile column names.
#' @param n Number of samples to draw.
#' @param bounds Optional length-2 vector of (lower, upper) bounds. Required
#'   for method = "beta"; used as padding bounds for "spline" if given;
#'   required for "metalog" unless `boundedness = "u"`.
#' @param ... Additional method-specific arguments (`sample_size` for "kde";
#'   `boundedness`/`term` for "metalog").
#' @return A numeric vector of `n` simulated values.
#' @export
simulate_from_percentiles <- function(quantile_df,
                                       method = c("linear_cdf", "spline", "kde", "metalog", "beta", "normal"),
                                       percentile_cols = c("P0", "P5", "P50", "P95", "P100"),
                                       n = 1000,
                                       bounds = NULL,
                                       ...) {
  method <- match.arg(method)
  pairs <- extract_percentile_pairs(quantile_df, percentile_cols)

  # Methods take different extra arguments (e.g. metalog's boundedness/term
  # vs kde's sample_size). Passing a shared `...` straight through would
  # error on whichever method doesn't declare a given name (this bit us in
  # testing: passing boundedness/term for metalog broke kde in the same
  # compare_percentile_methods() call). Instead, filter `...` down to just
  # the names each target function actually declares.
  extra_args <- list(...)
  filter_args <- function(fn) extra_args[names(extra_args) %in% names(formals(fn))]

  switch(method,
    linear_cdf = sim_linear_cdf(pairs$probs, pairs$values, n),
    spline     = do.call(sim_spline, c(list(pairs$probs, pairs$values, n, bounds = bounds), filter_args(sim_spline))),
    kde        = do.call(sim_kde, c(list(pairs$probs, pairs$values, n), filter_args(sim_kde))),
    metalog    = do.call(sim_metalog, c(list(pairs$probs, pairs$values, n, bounds = bounds), filter_args(sim_metalog))),
    beta       = sim_beta(pairs$probs, pairs$values, n, bounds = bounds),
    normal     = sim_normal(pairs$probs, pairs$values, n)
  )
}

#' Build a piecewise-linear inverse-CDF sampler from percentile columns
#'
#' Backwards-compatible wrapper around
#' `simulate_from_percentiles(..., method = "linear_cdf")`. Kept because
#' code/notebooks/soilgrids_simulation.Rmd calls this function by name.
#'
#' @inheritParams simulate_from_percentiles
#' @export
generate_inverse_cdf_distribution <- function(quantile_df,
                                              percentile_cols = c("P0","P5","P50","P95","P100"),
                                              n = 1000) {
  simulate_from_percentiles(quantile_df, method = "linear_cdf", percentile_cols = percentile_cols, n = n)
}

#' Run several percentile-reconstruction methods on the same data and compare them
#'
#' Draws `n` samples from each requested method and returns both the raw
#' samples and a summary-statistics table (via `calculate_summary_statistics()`)
#' side by side, so the shape/spread of each reconstruction can be compared
#' by eye or by downstream tests (e.g. `ks.test()` between pairs of methods).
#'
#' @inheritParams simulate_from_percentiles
#' @param methods Character vector of methods to run (see
#'   `simulate_from_percentiles`). "metalog" is deliberately excluded from
#'   the default - see the KNOWN CRASH RISK note on `sim_metalog()` above
#'   before opting into it.
#' @return A list with `samples` (named list of numeric vectors, one per
#'   method) and `summary` (data frame of summary statistics, one row per method).
#' @export
compare_percentile_methods <- function(quantile_df,
                                        methods = c("linear_cdf", "spline", "kde"),
                                        percentile_cols = c("P0", "P5", "P50", "P95", "P100"),
                                        n = 1000,
                                        bounds = NULL,
                                        ...) {
  samples <- setNames(
    lapply(methods, function(m) {
      tryCatch(
        simulate_from_percentiles(quantile_df, method = m, percentile_cols = percentile_cols, n = n, bounds = bounds, ...),
        error = function(e) {
          warning(sprintf("method '%s' failed: %s", m, conditionMessage(e)))
          NA_real_
        }
      )
    }),
    methods
  )

  summary_df <- do.call(rbind, lapply(methods, function(m) {
    s <- samples[[m]]
    if (length(s) == 1 && is.na(s)) return(NULL)
    cbind(method = m, calculate_summary_statistics(s))
  }))
  rownames(summary_df) <- NULL

  list(samples = samples, summary = summary_df)
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

#' Evaluate percentile-reconstruction methods against a known ground-truth distribution
#'
#' For each replicate: draws a large sample from `true_rng`, computes the
#' requested percentiles from it (simulating "what we'd observe" as summary
#' data), reconstructs a distribution from those percentiles using each
#' method, then scores the reconstruction against the true draws via the
#' Kolmogorov-Smirnov statistic and relative mean/SD error. Averages scores
#' across replicates.
#'
#' @param true_rng A function(n) that draws n samples from the ground-truth distribution.
#' @param percentile_probs Numeric vector of probabilities (0-1) to reveal as "known" percentiles.
#' @param methods Character vector of methods to evaluate. "metalog" is
#'   deliberately excluded from the default - see the KNOWN CRASH RISK note
#'   on `sim_metalog()` above before opting into it.
#' @param n_true Number of ground-truth draws per replicate (used both to
#'   derive the revealed percentiles and as the KS-test reference sample).
#' @param n_sim Number of samples to reconstruct per method per replicate.
#' @param n_reps Number of independent replicates to average over.
#' @param bounds Optional bounds passed through to methods that use them.
#' @return A data frame with one row per method: mean KS statistic, mean
#'   absolute relative mean error, and mean absolute relative SD error across replicates.
#' @export
validate_percentile_methods_synthetic <- function(true_rng,
                                                    percentile_probs = c(0, 0.05, 0.5, 0.95, 1),
                                                    methods = c("linear_cdf", "spline", "kde"),
                                                    n_true = 5000,
                                                    n_sim = 1000,
                                                    n_reps = 20,
                                                    bounds = NULL) {
  percentile_cols <- paste0("P", round(percentile_probs * 100))

  results <- lapply(seq_len(n_reps), function(rep_i) {
    truth <- true_rng(n_true)
    percentile_values <- quantile(truth, probs = percentile_probs, names = FALSE)
    quantile_df <- as.data.frame(setNames(as.list(percentile_values), percentile_cols))

    rows <- lapply(methods, function(m) {
      recon <- tryCatch(
        simulate_from_percentiles(quantile_df, method = m, percentile_cols = percentile_cols, n = n_sim, bounds = bounds),
        error = function(e) NULL
      )
      if (is.null(recon) || length(recon) == 0 || all(is.na(recon))) {
        return(data.frame(method = m, ks_stat = NA_real_, mean_rel_err = NA_real_, sd_rel_err = NA_real_))
      }
      ks <- suppressWarnings(ks.test(recon, truth)$statistic)
      mean_rel_err <- abs(mean(recon) - mean(truth)) / abs(mean(truth))
      sd_rel_err <- abs(sd(recon) - sd(truth)) / sd(truth)
      data.frame(method = m, ks_stat = unname(ks), mean_rel_err = mean_rel_err, sd_rel_err = sd_rel_err)
    })
    do.call(rbind, rows)
  })

  combined <- do.call(rbind, results)
  agg <- aggregate(cbind(ks_stat, mean_rel_err, sd_rel_err) ~ method, data = combined, FUN = mean, na.rm = TRUE)
  agg[order(agg$ks_stat), ]
}
