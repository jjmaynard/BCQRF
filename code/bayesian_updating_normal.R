# Load necessary library
library(ggplot2)

# Define the parameters for the prior distribution (based on empirical data)
prior_mean <- 3.5
prior_sd <- 1.0

# Define the parameters for the likelihood distribution (based on model predictions)
likelihood_mean <- 4.0
likelihood_sd <- 0.3

# Bayesian updating for the posterior mean and standard deviation
# Formula for Bayesian updating when both prior and likelihood are normal distributions:
# Posterior mean = (prior_mean / prior_variance + likelihood_mean / likelihood_variance) /
#                  (1 / prior_variance + 1 / likelihood_variance)
# Posterior variance = 1 / (1 / prior_variance + 1 / likelihood_variance)

# Calculate variances
prior_var <- prior_sd^2
likelihood_var <- likelihood_sd^2

# Calculate posterior mean and variance
posterior_mean <- (prior_mean / prior_var + likelihood_mean / likelihood_var) /
  (1 / prior_var + 1 / likelihood_var)
posterior_var <- 1 / (1 / prior_var + 1 / likelihood_var)
posterior_sd <- sqrt(posterior_var)

# Print posterior mean and standard deviation
cat("Posterior Mean:", posterior_mean, "\n")
cat("Posterior SD:", posterior_sd, "\n")

# Generate samples for visualization
set.seed(42)  # for reproducibility
n_samples <- 1000
prior_samples <- rnorm(n_samples, mean = prior_mean, sd = prior_sd)
likelihood_samples <- rnorm(n_samples, mean = likelihood_mean, sd = likelihood_sd)
posterior_samples <- rnorm(n_samples, mean = posterior_mean, sd = posterior_sd)

# Combine into a data frame for plotting
data <- data.frame(
  Value = c(prior_samples, likelihood_samples, posterior_samples),
  Distribution = factor(rep(c("Prior (Empirical)", "Likelihood (Model)", "Posterior (Constrained)"),
                            each = n_samples))
)

# Plot the distributions
ggplot(data, aes(x = Value, fill = Distribution)) +
  geom_density(alpha = 0.5) +
  labs(title = "Bayesian Updating to Constrain Model Predictions",
       x = "Soil Organic Carbon (%)",
       y = "Density") +
  theme_minimal() +
  scale_fill_manual(values = c("skyblue", "salmon", "lightgreen"))







# Load necessary library
library(ggplot2)

# Define the empirical prior distribution
prior_mean <- 3.5
prior_sd <- 1.0

# Define the model predictions (biased likelihood)
biased_likelihood_mean <- 6.0  # Significantly biased mean
biased_likelihood_sd <- 0.3    # High precision (low variance)

# Inflate the likelihood variance to reduce its influence
inflation_factor <- 4
adjusted_likelihood_sd <- sqrt(biased_likelihood_sd^2 * inflation_factor)

# Bayesian updating
prior_var <- prior_sd^2
adjusted_likelihood_var <- adjusted_likelihood_sd^2

posterior_mean <- (prior_mean / prior_var + biased_likelihood_mean / adjusted_likelihood_var) /
  (1 / prior_var + 1 / adjusted_likelihood_var)
posterior_var <- 1 / (1 / prior_var + 1 / adjusted_likelihood_var)
posterior_sd <- sqrt(posterior_var)

# Truncate the posterior within the empirical range
posterior_samples <- rnorm(1000, mean = posterior_mean, sd = posterior_sd)
empirical_min <- prior_mean - 3 * prior_sd
empirical_max <- prior_mean + 3 * prior_sd
posterior_samples_truncated <- posterior_samples[posterior_samples >= empirical_min &
                                                   posterior_samples <= empirical_max]

# Combine into a data frame for visualization
data <- data.frame(
  Value = c(rnorm(984, mean = prior_mean, sd = prior_sd),
            rnorm(984, mean = biased_likelihood_mean, sd = biased_likelihood_sd),
            posterior_samples_truncated),
  Distribution = factor(rep(c("Prior (Empirical)",
                              "Likelihood (Biased Model)",
                              "Posterior (Constrained)"), each = length(posterior_samples_truncated)))
)

# Plot the distributions
ggplot(data, aes(x = Value, fill = Distribution)) +
  geom_density(alpha = 0.5) +
  labs(title = "Bayesian Updating with Biased Predictions",
       x = "Soil Organic Carbon (%)",
       y = "Density") +
  theme_minimal() +
  scale_fill_manual(values = c("skyblue", "salmon", "lightgreen"))




# Load necessary libraries
library(ggplot2)
library(dplyr)

# Define the empirical (prior) and model prediction (likelihood) distributions
# Simulate non-normal distributions for demonstration purposes
set.seed(42)
empirical_prior <- c(rnorm(500, mean = 3.5, sd = 0.8), runif(200, min = 1.5, max = 2.5))
model_predictions <- c(rnorm(700, mean = 6.0, sd = 0.4), rnorm(300, mean = 5.5, sd = 0.3))

# Define a grid of values for evaluation
value_grid <- seq(0, 8, by = 0.01)

# Estimate densities for the prior and likelihood
prior_density <- density(empirical_prior, bw = "nrd0", from = 0, to = 8, n = length(value_grid))
likelihood_density <- density(model_predictions, bw = "nrd0", from = 0, to = 8, n = length(value_grid))

# Interpolate densities to match the grid
prior_prob <- approxfun(prior_density$x, prior_density$y)(value_grid)
likelihood_prob <- approxfun(likelihood_density$x, likelihood_density$y)(value_grid)

# Normalize densities to ensure they sum to 1 (if not already normalized)
prior_prob <- prior_prob / sum(prior_prob)
likelihood_prob <- likelihood_prob / sum(likelihood_prob)

# Compute the posterior using Bayes' rule: posterior ∝ prior * likelihood
posterior_prob <- prior_prob * likelihood_prob
posterior_prob <- posterior_prob / sum(posterior_prob)  # Normalize posterior

# Create a data frame for visualization
posterior_samples <- sample(value_grid, size = 1000, prob = posterior_prob, replace = TRUE)

data <- data.frame(
  Value = c(empirical_prior, model_predictions, posterior_samples),
  Distribution = factor(c(
    rep("Prior (SSURGO-Mapunit)", length(empirical_prior)),
    rep("Likelihood (SOLUS-QRF predictions at selected pixel)", length(model_predictions)),
    rep("Posterior (SSURGO Constrained SOLUS)", length(posterior_samples))
  ))
)

# Plot the distributions
ggplot(data, aes(x = Value, fill = Distribution)) +
  geom_density(alpha = 0.5) +
  labs(title = "Bayesian Updating",
       x = "Soil Organic Carbon (%)",
       y = "Density") +
  theme_minimal() +
  scale_fill_manual(values = c("skyblue", "salmon", "lightgreen")) +
  scale_x_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 1))
