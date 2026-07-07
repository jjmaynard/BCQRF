#' @title Grid-Based Bayesian Updating
#'
#' @description Canonical `bayesian_update()` function, extracted during the code/
#' reorganization. This exact function previously existed as three verbatim (or
#' near-verbatim) copies in code/scratch/bayesian_updating_normal.R (as inline,
#' non-function demo code), code/notebooks/constrained_qrf_ideas.Rmd, and
#' code/notebooks/soilgrids_simulation.Rmd. Those notebooks now source() this file
#' instead of redefining the function.

#' Combine a prior and likelihood distribution into a posterior via grid-based KDE
#'
#' Estimates prior and likelihood densities over a shared value grid via kernel
#' density estimation, multiplies them per Bayes' rule, and samples from the
#' resulting (normalized) posterior.
#'
#' @param prior_distribution Numeric vector of prior samples.
#' @param likelihood_distribution Numeric vector of likelihood samples.
#' @param grid_range Optional length-2 vector giving the grid's min/max. If NULL,
#'   derived from the combined range of both inputs, padded by 1.
#' @param grid_resolution Step size of the evaluation grid.
#' @return A numeric vector of 1000 samples drawn from the posterior distribution.
#' @export
bayesian_update <- function(prior_distribution, likelihood_distribution, grid_range = NULL, grid_resolution = 0.01) {
  # Determine the range of the grid dynamically if not provided
  if (is.null(grid_range)) {
    grid_min <- min(c(prior_distribution, likelihood_distribution)) - 1
    grid_max <- max(c(prior_distribution, likelihood_distribution)) + 1
  } else {
    grid_min <- grid_range[1]
    grid_max <- grid_range[2]
  }

  # Define the value grid
  value_grid <- seq(grid_min, grid_max, by = grid_resolution)

  # Estimate densities for the prior and likelihood
  prior_density <- density(prior_distribution, bw = "nrd0", from = grid_min, to = grid_max, n = length(value_grid))
  likelihood_density <- density(likelihood_distribution, bw = "nrd0", from = grid_min, to = grid_max, n = length(value_grid))

  # Interpolate densities to match the grid
  prior_prob <- approxfun(prior_density$x, prior_density$y)(value_grid)
  likelihood_prob <- approxfun(likelihood_density$x, likelihood_density$y)(value_grid)

  # Normalize densities to ensure they sum to 1 (if not already normalized)
  prior_prob <- prior_prob / sum(prior_prob)
  likelihood_prob <- likelihood_prob / sum(likelihood_prob)

  # Compute the posterior using Bayes' rule: posterior proportional to prior * likelihood
  posterior_prob <- prior_prob * likelihood_prob
  posterior_prob <- posterior_prob / sum(posterior_prob)  # Normalize posterior

  # Sample from the posterior
  posterior_samples <- sample(value_grid, size = 1000, prob = posterior_prob, replace = TRUE)

  # Return the posterior samples as the constrained distribution
  return(posterior_samples)
}
