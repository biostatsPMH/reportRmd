# ggkmcif3 Testing Scripts
# Comprehensive test suite for all parameters in ggkmcif3Parameters

# Load required libraries
library(survival)
library(ggplot2)
library(dplyr)
library(cowplot)

# Create sample datasets for testing ----

# Sample KM data (time-to-event with censoring)
set.seed(123)
n <- 200
km_data <- data.frame(
  time = rexp(n, rate = 0.1) |> round(2),
  status = rbinom(n, 1, 0.7),  # 70% event rate
  treatment = factor(sample(c("Control", "Treatment"), n, replace = TRUE)),
  age_group = factor(sample(c("Young", "Old"), n, replace = TRUE)),
  age_cont = rnorm(n, 60, 15) |> round(1),
  sex = factor(sample(c("Male", "Female"), n, replace = TRUE))
)

# Sample CIF data (competing risks)
cif_data <- data.frame(
  time = rexp(n, rate = 0.1) |> round(2),
  status = sample(0:2, n, replace = TRUE, prob = c(0.3, 0.4, 0.3)),  # 0=censor, 1=event1, 2=event2
  treatment = factor(sample(c("Control", "Treatment"), n, replace = TRUE)),
  age_group = factor(sample(c("Young", "Old"), n, replace = TRUE)),
  risk_score = rnorm(n, 0, 1) |> round(2)
)

# Quick Test
data("pembrolizumab")
data=pembrolizumab;cov="cohort";response=c("os_time","os_status")
ggkmcif3(response,cov,data,censor.size=2.5,censor.stroke=7)

ggkmcif3(response,cov,data,times=seq(0,42,6))
ggkmcif3(response,cov,data,times=seq(0,42,6),xlim=c(0,30))

# Test Helper Functions ----

#' Run a single test case
#' @param test_name Character string describing the test
#' @param data Data frame to use
#' @param response Character vector with time and status columns
#' @param test_args List of arguments to pass to ggkmcif3
run_test <- function(test_name, data, response, test_args = list()) {
  cat("\n", "=", "\n")
  cat("Testing:", test_name, "\n")
  cat("=", "\n")

  tryCatch({
    args <- c(list(response = response, data = data), test_args)
    return(args)
    result <- do.call(ggkmcif3, args)
    plot(result)
    cat("✓ Test passed successfully\n")
    return(result)
  }, error = function(e) {
    cat("✗ Test failed with error:", e$message, "\n")
    return(NULL)
  })
}

# Basic Functionality Tests ----

test_basic_functionality <- function() {
  cat("\n", "#", "\n")
  cat("# BASIC FUNCTIONALITY TESTS\n")
  cat("#", "\n")

  # Test 1: Basic KM plot
args <-  run_test("Basic KM plot - no covariates",
           km_data,
           c("time", "status"))
do.call(ggkmcif3,args)
ggkmcif3(data=km_data,response = c("time", "status"),cov="age_group")

  # Test 2: Basic KM plot with covariate
  run_test("Basic KM plot - with treatment covariate",
           km_data,
           c("time", "status"),
           list(cov = "treatment"))

  # Test 3: Basic CIF plot
  run_test("Basic CIF plot - no covariates",
           cif_data,
           c("time", "status"),
           list(type = "CIF"))

  # Test 4: Basic CIF plot with covariate
  run_test("Basic CIF plot - with treatment covariate",
           cif_data,
           c("time", "status"),
           list(type = "CIF", cov = "treatment"))

  # Test 5: Auto-detection of plot type
  run_test("Auto-detect KM type",
           km_data,
           c("time", "status"),
           list(cov = "treatment"))

  run_test("Auto-detect CIF type",
           cif_data,
           c("time", "status"),
           list(cov = "treatment"))
}

# Plot Customization Tests ----

test_plot_customization <- function() {
  cat("\n", "#", "\n")
  cat("# PLOT CUSTOMIZATION TESTS\n")
  cat("#", "\n")

  # Test colors
  run_test("Custom colors",
           km_data,
           c("time", "status"),
           list(cov = "treatment", col = c("red", "blue")))

  # Test line types
  run_test("Custom line types",
           km_data,
           c("time", "status"),
           list(cov = "treatment", linetype = c(1, 2)))

  # Test labels
  run_test("Custom labels",
           km_data,
           c("time", "status"),
           list(cov = "treatment",
                xlab = "Follow-up Time (years)",
                ylab = "Overall Survival",
                main = "Survival Analysis by Treatment"))

  # Test strata labels
  run_test("Custom strata labels",
           km_data,
           c("time", "status"),
           list(cov = "treatment",
                stratalabs = c("Standard Care", "New Treatment"),
                strataname = "Treatment Group"))

  # Test sizing parameters
  run_test("Custom sizing",
           km_data,
           c("time", "status"),
           list(cov = "treatment",
                fsize = 14,
                lsize = 1.2,
                censor.size = 1.0))
}

# Statistical Features Tests ----

test_statistical_features <- function() {
  cat("\n", "#", "\n")
  cat("# STATISTICAL FEATURES TESTS\n")
  cat("#", "\n")

  # Test p-values
  run_test("P-values enabled",
           km_data,
           c("time", "status"),
           list(cov = "treatment", pval = TRUE))

  run_test("P-values disabled",
           km_data,
           c("time", "status"),
           list(cov = "treatment", pval = FALSE))

  # Test confidence curves
  run_test("Confidence curves",
           km_data,
           c("time", "status"),
           list(cov = "treatment",
                conf.curves = TRUE,
                conf.type = "log"))

  # Test different confidence types
  conf_types <- c("plain", "log", "log-log", "logit")
  for(conf_type in conf_types) {
    run_test(paste("Confidence type:", conf_type),
             km_data,
             c("time", "status"),
             list(cov = "treatment",
                  conf.curves = TRUE,
                  conf.type = conf_type))
  }

  # Test hazard ratios
  run_test("Hazard ratios",
           km_data,
           c("time", "status"),
           list(cov = "treatment",
                HR = TRUE,
                HR_pval = TRUE))
}

# Median and Set Time Tests ----

test_median_and_set_time <- function() {
  cat("\n", "#", "\n")
  cat("# MEDIAN AND SET TIME TESTS\n")
  cat("#", "\n")

  # Test median features
  run_test("Median text",
           km_data,
           c("time", "status"),
           list(cov = "treatment",
                median.text = TRUE,
                median.CI = TRUE))

  run_test("Median lines",
           km_data,
           c("time", "status"),
           list(cov = "treatment",
                median.lines = TRUE,
                median.lsize = 1.5))

  # Test set time features
  run_test("Set time analysis - single time",
           km_data,
           c("time", "status"),
           list(cov = "treatment",
                set.time = 5,
                set.time.text = "5-year survival",
                set.time.line = TRUE,
                set.time.CI = TRUE))

  run_test("Set time analysis - multiple times",
           km_data,
           c("time", "status"),
           list(cov = "treatment",
                set.time = c(2, 5, 10),
                set.time.text = "Survival at key timepoints",
                set.time.line = TRUE))
}

# Risk Table Tests ----

test_risk_table <- function() {
  cat("\n", "#", "\n")
  cat("# RISK TABLE TESTS\n")
  cat("#", "\n")

  # Test risk table enabled/disabled
  run_test("Risk table enabled",
           km_data,
           c("time", "status"),
           list(cov = "treatment", table = TRUE))

  run_test("Risk table disabled",
           km_data,
           c("time", "status"),
           list(cov = "treatment", table = FALSE))

  # Test risk table customization
  run_test("Custom risk table labels",
           km_data,
           c("time", "status"),
           list(cov = "treatment",
                stratalabs.table = c("Ctrl", "Trt"),
                strataname.table = "Group",
                Numbers_at_risk_text = "Number at risk"))

  # Test risk table height
  run_test("Custom risk table height",
           km_data,
           c("time", "status"),
           list(cov = "treatment",
                tbl.height = 0.3))

  # Test with custom times
  run_test("Custom time points",
           km_data,
           c("time", "status"),
           list(cov = "treatment",
                times = c(0, 2, 4, 6, 8, 10)))
}

# CIF-Specific Tests ----

test_cif_specific <- function() {
  cat("\n", "#", "\n")
  cat("# CIF-SPECIFIC TESTS\n")
  cat("#", "\n")

  # Test different events
  run_test("CIF - Event 1",
           cif_data,
           c("time", "status"),
           list(type = "CIF",
                cov = "treatment",
                plot.event = 1))

  run_test("CIF - Event 2",
           cif_data,
           c("time", "status"),
           list(type = "CIF",
                cov = "treatment",
                plot.event = 2))

  # Test multiple events
  run_test("CIF - Multiple events",
           cif_data,
           c("time", "status"),
           list(type = "CIF",
                cov = "treatment",
                plot.event = c(1, 2),
                event = "col"))

  # Test flipped CIF
  run_test("CIF - Flipped curves",
           cif_data,
           c("time", "status"),
           list(type = "CIF",
                cov = "treatment",
                flip.CIF = TRUE))

  # Test event labels
  run_test("CIF - Custom event labels",
           cif_data,
           c("time", "status"),
           list(type = "CIF",
                cov = "treatment",
                plot.event = c(1, 2),
                eventlabs = c("Death", "Relapse"),
                event.name = "Event Type"))
}

# Continuous Covariate Tests ----

test_continuous_covariates <- function() {
  cat("\n", "#", "\n")
  cat("# CONTINUOUS COVARIATE TESTS\n")
  cat("#", "\n")

  # Test continuous covariate with median cut
  run_test("Continuous covariate - median cut",
           km_data,
           c("time", "status"),
           list(cov = "age_cont"))

  # Test continuous covariate with custom cut
  run_test("Continuous covariate - custom cut",
           km_data,
           c("time", "status"),
           list(cov = "age_cont",
                cut = 65))

  # Test continuous covariate with custom labels
  run_test("Continuous covariate - custom labels",
           km_data,
           c("time", "status"),
           list(cov = "age_cont",
                cut = 65,
                stratalabs = c("Younger", "Older"),
                strataname = "Age Group"))
}

# Advanced Formatting Tests ----

test_advanced_formatting <- function() {
  cat("\n", "#", "\n")
  cat("# ADVANCED FORMATTING TESTS\n")
  cat("#", "\n")

  # Test legend positioning
  legend_positions <- c("left", "top", "right", "bottom", "none")
  for(pos in legend_positions) {
    run_test(paste("Legend position:", pos),
             km_data,
             c("time", "status"),
             list(cov = "treatment",
                  legend.pos = pos))
  }

  # Test custom legend position with coordinates
  run_test("Legend position - coordinates",
           km_data,
           c("time", "status"),
           list(cov = "treatment",
                legend.pos = c(0.8, 0.8)))

  # Test p-value positioning
  run_test("Custom p-value position",
           km_data,
           c("time", "status"),
           list(cov = "treatment",
                pval = TRUE,
                pval.pos = c(5, 0.2)))

  # Test axis limits
  run_test("Custom axis limits",
           km_data,
           c("time", "status"),
           list(cov = "treatment",
                xlim = c(0, 15),
                ylim = c(0.3, 1)))

  # Test digit precision
  run_test("Custom digit precision",
           km_data,
           c("time", "status"),
           list(cov = "treatment",
                HR = TRUE,
                HR_pval = TRUE,
                median.text = TRUE,
                HR.digits = 3,
                HR.pval.digits = 4,
                pval.digits = 4,
                median.digits = 2))
}

# Censor Mark Tests ----

test_censor_marks <- function() {
  cat("\n", "#", "\n")
  cat("# CENSOR MARK TESTS\n")
  cat("#", "\n")

  # Test censor marks enabled
  run_test("Censor marks - enabled",
           km_data,
           c("time", "status"),
           list(cov = "treatment",
                censor.marks = TRUE))

  # Test censor marks disabled
  run_test("Censor marks - disabled",
           km_data,
           c("time", "status"),
           list(cov = "treatment",
                censor.marks = FALSE))

  # Test custom censor mark styling
  run_test("Censor marks - custom styling",
           km_data,
           c("time", "status"),
           list(cov = "treatment",
                censor.marks = TRUE,
                censor.size = 2.0,
                censor.stroke = 2.5))
}

# Complex Scenario Tests ----

test_complex_scenarios <- function() {
  cat("\n", "#", "\n")
  cat("# COMPLEX SCENARIO TESTS\n")
  cat("#", "\n")

  # Test 1: Publication-ready KM plot
  run_test("Publication-ready KM plot",
           km_data,
           c("time", "status"),
           list(cov = "treatment",
                conf.curves = TRUE,
                conf.type = "log",
                pval = TRUE,
                median.text = TRUE,
                median.lines = TRUE,
                median.CI = TRUE,
                set.time = c(2, 5),
                set.time.text = "Key timepoints",
                set.time.line = TRUE,
                set.time.CI = TRUE,
                stratalabs = c("Standard of Care", "Experimental Treatment"),
                strataname = "Treatment Arm",
                xlab = "Time (years)",
                ylab = "Overall Survival Probability",
                main = "Kaplan-Meier Survival Analysis",
                fsize = 12,
                col = c("#E31A1C", "#1F78B4")))

  # Test 2: Complex CIF plot
  run_test("Complex CIF plot",
           cif_data,
           c("time", "status"),
           list(type = "CIF",
                cov = "treatment",
                plot.event = c(1, 2),
                event = "col",
                conf.curves = TRUE,
                pval = TRUE,
                eventlabs = c("Disease Progression", "Death"),
                event.name = "Competing Events",
                stratalabs = c("Control Group", "Treatment Group"),
                col = c("#FF7F00", "#33A02C", "#6A3D9A", "#B15928"),
                xlab = "Follow-up Time (months)",
                main = "Cumulative Incidence Analysis"))

  # Test 3: Multi-group comparison
  run_test("Multi-group comparison",
           km_data,
           c("time", "status"),
           list(cov = "age_group",
                pval = TRUE,
                conf.curves = TRUE,
                median.text = TRUE,
                table = TRUE,
                col = c("#1B9E77", "#D95F02", "#7570B3")))

  # Test 4: No risk table, custom positioning
  run_test("No risk table with custom text positioning",
           km_data,
           c("time", "status"),
           list(cov = "treatment",
                table = FALSE,
                pval = TRUE,
                median.text = TRUE,
                pval.pos = c(8, 0.9),
                median.pos = c(2, 0.3),
                legend.pos = "bottom"))
}

# Edge Case Tests ----

test_edge_cases <- function() {
  cat("\n", "#", "\n")
  cat("# EDGE CASE TESTS\n")
  cat("#", "\n")

  # Test with small dataset
  small_data <- km_data |> slice_head(n = 20)
  run_test("Small dataset",
           small_data,
           c("time", "status"),
           list(cov = "treatment"))

  # Test with no events
  no_events_data <- km_data |> mutate(status = 0)
  run_test("No events (all censored)",
           no_events_data,
           c("time", "status"),
           list(cov = "treatment"))

  # Test with single group
  single_group <- km_data |> filter(treatment == "Control")
  run_test("Single group only",
           single_group,
           c("time", "status"))

  # Test returns parameter
  result <- run_test("Returns parameter",
                     km_data,
                     c("time", "status"),
                     list(cov = "treatment",
                          returns = TRUE))

  if(!is.null(result)) {
    cat("✓ Returns object structure validated\n")
    print(str(result))
  }
}

# Performance and Stress Tests ----

test_performance <- function() {
  cat("\n", "#", "\n")
  cat("# PERFORMANCE TESTS\n")
  cat("#", "\n")

  # Test with larger dataset
  set.seed(456)
  large_data <- data.frame(
    time = rexp(1000, rate = 0.05),
    status = rbinom(1000, 1, 0.6),
    group = factor(sample(letters[1:5], 1000, replace = TRUE))
  )

  start_time <- Sys.time()
  run_test("Large dataset (n=1000)",
           large_data,
           c("time", "status"),
           list(cov = "group"))
  end_time <- Sys.time()
  cat("Execution time:", difftime(end_time, start_time, units = "secs"), "seconds\n")

  # Test with many time points
  run_test("Many time points",
           km_data,
           c("time", "status"),
           list(cov = "treatment",
                times = seq(0, 20, by = 0.5)))
}

# Main Test Runner ----

run_all_tests <- function() {
  cat("Starting comprehensive ggkmcif3 testing suite...\n")
  cat("Date:", Sys.time(), "\n")

  # Run all test categories
  test_basic_functionality()
  test_plot_customization()
  test_statistical_features()
  test_median_and_set_time()
  test_risk_table()
  test_cif_specific()
  test_continuous_covariates()
  test_advanced_formatting()
  test_censor_marks()
  test_complex_scenarios()
  test_edge_cases()
  test_performance()

  cat("\n", "=", "\n")
  cat("Testing suite completed!\n")
  cat("=", "\n")
}

# Parameter Validation Tests ----

test_parameter_validation <- function() {
  cat("\n", "#", "\n")
  cat("# PARAMETER VALIDATION TESTS\n")
  cat("#", "\n")

  # Test all ggkmcif3Parameters defaults
  default_params <- ggkmcif3Parameters(strataname = "test")

  cat("Testing all default parameters from ggkmcif3Parameters:\n")
  param_names <- names(default_params)

  for(param in param_names) {
    if(param %in% c("returns", "print.n.missing")) next  # Skip logical params that might cause issues

    test_args <- list()
    test_args[[param]] <- default_params[[param]]
    test_args$cov <- "treatment"

    run_test(paste("Parameter:", param),
             km_data,
             c("time", "status"),
             test_args)
  }
}

# Interactive Test Runner ----

run_specific_tests <- function(test_categories = NULL) {
  available_tests <- list(
    "basic" = test_basic_functionality,
    "customization" = test_plot_customization,
    "statistical" = test_statistical_features,
    "median" = test_median_and_set_time,
    "risktable" = test_risk_table,
    "cif" = test_cif_specific,
    "continuous" = test_continuous_covariates,
    "formatting" = test_advanced_formatting,
    "censor" = test_censor_marks,
    "complex" = test_complex_scenarios,
    "edge" = test_edge_cases,
    "performance" = test_performance,
    "validation" = test_parameter_validation
  )

  if(is.null(test_categories)) {
    cat("Available test categories:\n")
    for(name in names(available_tests)) {
      cat("-", name, "\n")
    }
    cat("\nUse run_specific_tests(c('basic', 'customization')) to run specific tests\n")
    cat("Use run_all_tests() to run everything\n")
    return(invisible())
  }

  for(category in test_categories) {
    if(category %in% names(available_tests)) {
      available_tests[[category]]()
    } else {
      cat("Unknown test category:", category, "\n")
    }
  }
}

run_all_tests()
