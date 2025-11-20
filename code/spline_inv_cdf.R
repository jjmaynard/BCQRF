# Load required package
library(splines)  # For spline interpolation

# Step 1: Define given percentiles and corresponding quantile values
percentiles <- c(0.05,  0.50,  0.95)  # Given percentiles
values <- c(10,  35,  70)  # Corresponding values from empirical data

# Step 2: Create an interpolated CDF function
# Extend the range slightly to cover full [0,1] probability space
extended_percentiles <- c(0, percentiles, 1)
extended_values <- c(min(values) - 5, values, max(values) + 5)  # Extend min/max slightly

# Use spline interpolation to create an approximate inverse CDF
inv_cdf <- splinefun(extended_percentiles, extended_values, method = "monoH.FC")

# Step 3: Perform inverse transform sampling
set.seed(123)  # For reproducibility
uniform_samples <- runif(10000)  # Generate 1000 random numbers from U(0,1)
empirical_samples <- inv_cdf(uniform_samples)  # Apply inverse CDF

# Step 4: Plot the results
hist(empirical_samples, breaks = 30, probability = TRUE, col = "lightblue",
     main = "Histogram of 1000 Samples from Empirical Distribution",
     xlab = "Sampled Values")
library(ggplot2)
ggplot(data.frame(x = empirical_samples), aes(x = x)) +
  geom_density(fill = "blue", alpha = 0.5) +
  labs(title = "Density Plot of Empirical Distribution",
       x = "Sampled Values",
       y = "Density") +
  theme_minimal()
