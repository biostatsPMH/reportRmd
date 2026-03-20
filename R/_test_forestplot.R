# Test Script for Improved Forest Plot Functions ----

# This script demonstrates the improved forestplot functions with:
# 1. Native pipe usage
# 2. Proper variable ordering (by max OR, then by factor levels within variables)
# 3. Integration of m_summary for robust data handling
# 4. Combined UV/MV plotting in single function

# Load required libraries ----
library(ggplot2)
library(dplyr)
library(scales)

# Source the improved functions ----
# source("forestplot_improved.R")
# source("autosum.R")  # For m_summary and helper functions

# Example with simulated data ----
set.seed(123)
n <- 200

# Create example dataset ----
example_data <- data.frame(
  # Response variable ----
  outcome = rbinom(n, 1, 0.4),

  # Categorical predictor with factor levels ----
  treatment = factor(
    sample(c("Control", "Treatment A", "Treatment B"), n, replace = TRUE),
    levels = c("Control", "Treatment A", "Treatment B")
  ),

  # Another categorical predictor ----
  stage = factor(
    sample(c("Early", "Advanced"), n, replace = TRUE),
    levels = c("Early", "Advanced")
  ),

  # Continuous predictor ----
  age = rnorm(n, 60, 10),

  # Binary predictor ----
  sex = factor(sample(c("Female", "Male"), n, replace = TRUE))
)

# Add some missing data to demonstrate N differences ----
example_data$age[sample(1:n, 10)] <- NA

# Fit multivariable model ----
mv_model <- glm(
  outcome ~ treatment + stage + age + sex,
  data = example_data,
  family = binomial()
)

# Display model summary ----
cat("\n=== Model Summary ===\n")
print(summary(mv_model))

# Example 1: Basic MV forest plot (adjusted estimates only) ----
cat("\n=== Example 1: Adjusted Estimates Only ===\n")
load_all()
p1 <- forestplotMV(
  model = mv_model,
  data = example_data,
  include_unadjusted = FALSE,
  showRef = TRUE,
  logScale = TRUE
)
print(p1)

# Example 2: MV forest plot with unadjusted estimates ----
cat("\n=== Example 2: Both Adjusted and Unadjusted Estimates ===\n")
p2 <- forestplotMV(
  model = mv_model,
  data = example_data,
  include_unadjusted = TRUE,
  showRef = TRUE,
  logScale = TRUE
)
print(p2)

# Example 3: UV forest plot (backwards compatibility) ----
cat("\n=== Example 3: Univariable Forest Plot (Backwards Compatible) ===\n")
p3 <- forestplotUV(
  response = "outcome",
  covs = c("treatment", "stage", "age", "sex"),
  data = example_data,
  family = binomial(),
  showRef = TRUE,
  logScale = TRUE
)
print(p3)

# Example 4: Hide reference levels ----
cat("\n=== Example 4: Hide Reference Levels ===\n")
p4 <- forestplotMV(
  model = mv_model,
  data = example_data,
  include_unadjusted = FALSE,
  showRef = FALSE,
  logScale = TRUE
)
print(p4)

# Example 5: Linear scale instead of log scale ----
cat("\n=== Example 5: Linear Scale ===\n")
p5 <- forestplotMV(
  model = mv_model,
  data = example_data,
  include_unadjusted = TRUE,
  showRef = TRUE,
  logScale = FALSE
)
print(p5)

# Example 6: Custom colors ----
cat("\n=== Example 6: Custom Colors ===\n")
p6 <- forestplotMV(
  model = mv_model,
  data = example_data,
  include_unadjusted = TRUE,
  colours = c("#4169E1", "#000000", "#DC143C"),  # Blue, Black, Crimson
  showRef = TRUE,
  logScale = TRUE
)
print(p6)

# Demonstration of ordering ----
cat("\n=== Ordering Demonstration ===\n")
cat("Variables are ordered by their maximum OR (highest to lowest)\n")
cat("Within each variable:\n")
cat("  - Reference level appears first (top)\n")
cat("  - Other levels follow factor order (for factors)\n")
cat("  - Or ordered by OR (for continuous/non-factors)\n\n")

# Show the deprecated function warning ----
cat("\n=== Deprecated Function Warning ===\n")
suppressWarnings({
  try(forestplotUVMV(p3, p1), silent = FALSE)
})

cat("\n=== All tests complete ===\n")
