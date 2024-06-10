library(reportRmd)
data("pembrolizumab")

test_that("rm_compact_summary works without a grouping variable", {
  output <- rm_compactsummary(pembrolizumab, xvars = c("age", "change_ctdna_group"), tableOnly = TRUE)
  expect_equal(output[["Covariate"]], c("age Median (Min, Max)", "change_ctdna_group n(%)"))
  expect_equal(output[["Full Sample (n=94)"]], c("59.1 (21.1, 81.8)", "40 (55)"))
  expect_equal(output[["Missing"]], c(0, 21))
})

test_that("rm_compact_summary works with a grouping variable", {
  output <- rm_compactsummary(pembrolizumab, xvars = c("age", "change_ctdna_group"), grp = "sex", tableOnly = TRUE)
  expect_equal(output[["Covariate"]], c("age Median (Min, Max)", "change_ctdna_group n(%)"))
  expect_equal(output[["Full Sample (n=94)"]], c("59.1 (21.1, 81.8)", "40 (55)"))
  expect_equal(output[["Female (n=58)"]], c("56.6 (34.1, 78.2)", "21 (52)"))
  expect_equal(output[["Male (n=36)"]], c("61.2 (21.1, 81.8)", "19 (58)"))
  expect_equal(output[["p-value"]], c("0.30", "0.84"))
  expect_equal(output[["Missing"]], c(0, 21))
})

test_that("rm_compact_summary works with non-logical use_mean and vector digits", {
  output <- rm_compactsummary(data = pembrolizumab, xvars = c("age", "change_ctdna_group", "l_size"), grp = "sex", use_mean = "age", tableOnly = TRUE, digits = c("age" = 2, "l_size" = 3), digits.cat = 1, iqr = TRUE, show.tests = TRUE)
  expect_equal(output[["Covariate"]], c("age Mean (sd)", "change_ctdna_group n(%)", "l_size Median (Q1, Q3)"))
  expect_equal(output[["Full Sample (n=94)"]], c("57.86 (12.75)", "40 (54.8)", "73.500 (49.250, 108.750)"))
  expect_equal(output[["Female (n=58)"]], c("56.95 (12.59)", "21 (52.5)", "68.000 (44.250, 97.750)"))
  expect_equal(output[["Male (n=36)"]], c("59.32 (13.05)", "19 (57.6)", "93.000 (65.500, 121.000)"))
  expect_equal(output[["p-value"]], c("0.39", "0.84", "0.066"))
  expect_equal(output[["Missing"]], c(0, 21, 0))
  expect_equal(output[["pTest"]], c("t-test", "ChiSq", "Wilcoxon Rank Sum"))
})

test_that("rm_compactsummary works with mixed use_mean and digits", {
  output <- rm_compactsummary(data = pembrolizumab, xvars = c("tmb", "l_size", "baseline_ctdna", "orr"), grp = "sex", use_mean = c("tmb", "l_size"), tableOnly = TRUE, digits = c("tmb" = 3, "baseline_ctdna" = 2), pvalue = TRUE, show.tests = FALSE)
  expect_equal(ncol(output), 5)
  expect_equal(output[["Covariate"]], c("tmb Mean (sd)", "l_size Mean (sd)", "baseline_ctdna Median (Min, Max)", "orr n(%)"))
  expect_equal(output[["Full Sample (n=94)"]], c("0.911 (0.969)", "87.9 (59.6)", "86.03 (0.00, 4475.01)", "78 (83)"))
  expect_equal(output[["p-value"]], c("0.76", "0.23", "0.97", "0.18"))
})

test_that("rm_compactsummary works with effSize = TRUE, show.tests = TRUE", {
  output <- rm_compactsummary(data = pembrolizumab, xvars = c("age", "change_ctdna_group"), grp = "cohort", use_mean = "age", tableOnly = TRUE, digits = 2, effSize = TRUE, show.tests = TRUE)
  expect_equal(output[["Covariate"]], c("age Mean (sd)", "change_ctdna_group n(%)"))
  expect_equal(output[["Full Sample (n=94)"]], c("57.86 (12.75)", "40 (55)"))
  expect_equal(output[["A (n=16)"]], c("62.88 (6.08)", "8 (57)"))
  expect_equal(output[["B (n=18)"]], c("56.29 (13.99)", "7 (64)"))
  expect_equal(output[["C (n=18)"]], c("57.89 (10.83)", "5 (50)"))
  expect_equal(output[["D (n=12)"]], c("63.84 (9.98)", "2 (20)"))
  expect_equal(output[["E (n=30)"]], c("53.70 (15.26)", "18 (64)"))
  expect_equal(output[["p-value"]], c("0.069", "0.18"))
  expect_equal(output[["Effect Size (95% CI)"]], c("0.05 (0.00, 0.18)", "0.30 (0.23, 0.52)"))
  expect_equal(output[["Missing"]], c(0, 21))
  expect_equal(output[["pTest"]], c("ANOVA", "Fisher Exact"))
  expect_equal(output[["effStat"]], c("Omega Sq", "Cramer's V"))
})

test_that("rm_compactsummary works with effSize = TRUE, show.tests = TRUE, but pvalue = FALSE", {
  output <- rm_compactsummary(data = pembrolizumab, xvars = c("age", "change_ctdna_group"), grp = "cohort", use_mean = "age", tableOnly = TRUE, digits = 2, pvalue = FALSE, effSize = TRUE, show.tests = TRUE)
  expect_equal(ncol(output), 8)
  expect_equal(output[["Effect Size (95% CI)"]], NULL)
})

test_that("rm_compactsummary works with effSize = TRUE and pvalue = TRUE, but show.tests = FALSE", {
  output <- rm_compactsummary(data = pembrolizumab, xvars = c("tmb", "l_size"), grp = "cohort", tableOnly = TRUE, digits = 2, pvalue = TRUE, effSize = TRUE, show.tests = FALSE)
  expect_equal(ncol(output), 9)
  expect_equal(output[["Effect Size (95% CI)"]], c("0.19 (0.07, 0.41)", "0.31 (0.14, 0.58)"))
  expect_equal(output[["pTest"]], NULL)
  expect_equal(output[["effStat"]], NULL)
})

## Extra Tests
# check sd == 0
temp <- data.frame()
temp <- pembrolizumab[, "age"]
temp[["age"]] <- rep(19, nrow(pembrolizumab))
test_that("rm_compactsummary works with sd == 0", {
  output <- rm_compactsummary(data = temp, xvars = "age", use_mean = "age", tableOnly = TRUE)
  expect_equal(output[["Full Sample (n=94)"]], "19.0 (0.0)")
})

