# Shared simulated data ------------------------------------------------------

sim_uv_data <- function(n = 300, seed = 11) {
  set.seed(seed)
  d <- data.frame(
    age = rnorm(n, 60, 10),
    tmb = rexp(n, 1),
    sex = factor(sample(c("M", "F"), n, TRUE)),
    grp = factor(sample(c("A", "B", "C"), n, TRUE)),
    chr = sample(c("x", "y"), n, TRUE),
    dt = as.Date("2020-01-01") + sample(1:1000, n, TRUE),
    stringsAsFactors = FALSE
  )
  d$y <- 2 * d$age + 5 * (d$sex == "M") + rnorm(n, 0, 20)
  d$bin <- rbinom(n, 1, stats::plogis(scale(d$y)[, 1]))
  d$cnt <- rpois(n, exp(1 + 0.01 * d$age))
  d$ordy <- factor(sample(c("lo", "mid", "hi"), n, TRUE),
                   levels = c("lo", "mid", "hi"), ordered = TRUE)
  d$time <- rexp(n, 0.05)
  d$status <- rbinom(n, 1, 0.7)
  d$one <- factor("only")
  d
}

est_of <- function(tab, variable, col = 2) {
  x <- tab[[col]][tab[[1]] == variable]
  as.numeric(sub("^\\s*([-0-9.e]+).*", "\\1", x[1]))
}


# Estimates -------------------------------------------------------------------

test_that("linear estimates match lm", {
  d <- sim_uv_data()
  tab <- rm_uvsum(response = "y", covs = c("age", "sex"), data = d, tableOnly = TRUE)
  expect_equal(est_of(tab, "age"), unname(round(coef(lm(y ~ age, data = d))[2], 2)))
  expect_equal(est_of(tab, "M"), unname(round(coef(lm(y ~ sex, data = d))[2], 2)))
})

test_that("logistic and poisson estimates are exponentiated", {
  d <- sim_uv_data()
  lg <- rm_uvsum(response = "bin", covs = "age", data = d, tableOnly = TRUE)
  expect_match(names(lg)[2], "OR")
  expect_equal(est_of(lg, "age"),
               unname(round(exp(coef(glm(bin ~ age, data = d, family = binomial))[2]), 2)))

  ps <- rm_uvsum(response = "cnt", covs = "age", data = d, tableOnly = TRUE,
                 type = "poisson")
  expect_match(names(ps)[2], "RR")
  expect_equal(est_of(ps, "age"),
               unname(round(exp(coef(glm(cnt ~ age, data = d, family = poisson))[2]), 2)))
})

test_that("two response variables give a survival model", {
  skip_if_not_installed("survival")
  d <- sim_uv_data()
  tab <- suppressWarnings(
    rm_uvsum(response = c("time", "status"), covs = "age", data = d, tableOnly = TRUE))
  expect_match(names(tab)[2], "HR")
  expect_true("Event" %in% names(tab))
})

test_that("CIwidth is honoured", {
  d <- sim_uv_data()
  tab <- rm_uvsum(response = "y", covs = "age", data = d, tableOnly = TRUE, CIwidth = 0.90)
  expect_match(names(tab)[2], "90%CI", fixed = TRUE)
  ci <- round(confint(lm(y ~ age, data = d), level = 0.90)[2, ], 2)
  cell <- tab[[2]][tab[[1]] == "age"]
  expect_equal(as.numeric(regmatches(cell, gregexpr("[-0-9.]+", cell))[[1]][2:3]),
               unname(ci))
})


# Covariates that have to be dropped ------------------------------------------
# All three of these used to reach the model-fitting code and error, because
# the reduced covariate list was not passed on.

test_that("date covariates are dropped, not fit", {
  d <- sim_uv_data()
  expect_message(
    rm_uvsum(response = "y", covs = c("age", "dt"), data = d, tableOnly = TRUE),
    "Dates can not be used")
  tab <- suppressMessages(
    rm_uvsum(response = "y", covs = c("age", "dt"), data = d, tableOnly = TRUE))
  expect_false("dt" %in% tab[[1]])
  expect_true("age" %in% tab[[1]])
})

test_that("the response is dropped if it also appears in covs", {
  d <- sim_uv_data()
  expect_warning(
    tab <- rm_uvsum(response = "y", covs = c("y", "age"), data = d, tableOnly = TRUE),
    "is the response")
  expect_false("y" %in% tab[[1]])
  expect_true("age" %in% tab[[1]])
})

test_that("covariates with a single observed value are dropped", {
  d <- sim_uv_data()
  expect_warning(
    tab <- rm_uvsum(response = "y", covs = c("age", "one"), data = d, tableOnly = TRUE),
    "only one unique value")
  expect_false("one" %in% tab[[1]])
  expect_equal(nrow(tab), 1L)
})

test_that("dropping every covariate gives a clear error", {
  d <- sim_uv_data()
  expect_error(
    suppressWarnings(rm_uvsum(response = "y", covs = "one", data = d, tableOnly = TRUE)),
    "No covariates remain")
})


# reflevel --------------------------------------------------------------------

test_that("reflevel reorders an ordered response", {
  # relevel() refuses ordered factors, so this used to be a hard error
  skip_if_not_installed("MASS")
  d <- sim_uv_data()
  expect_no_error(suppressMessages(
    rm_uvsum(response = "ordy", covs = "age", data = d,
             tableOnly = TRUE, reflevel = "hi")))
  expect_equal(levels(reportRmd:::set_ref_level(d$ordy, "hi")),
               c("hi", "lo", "mid"))
  expect_true(is.ordered(reportRmd:::set_ref_level(d$ordy, "hi")))
  # an unordered factor keeps its class too
  expect_false(is.ordered(reportRmd:::set_ref_level(d$sex, "M")))
  expect_equal(levels(reportRmd:::set_ref_level(d$sex, "M")), c("M", "F"))
})

test_that("an invalid reflevel is reported clearly", {
  d <- sim_uv_data()
  expect_error(
    rm_uvsum(response = "ordy", covs = "age", data = d, tableOnly = TRUE, reflevel = "zzz"),
    "not a level of the response")
})


# removeInf -------------------------------------------------------------------

test_that("removeInf blanks unstable estimates", {
  d <- sim_uv_data()
  d$sep <- d$bin   # perfect separation

  keep <- suppressWarnings(
    rm_uvsum(response = "bin", covs = c("sep", "age"), data = d, tableOnly = TRUE))
  drop <- suppressWarnings(suppressMessages(
    rm_uvsum(response = "bin", covs = c("sep", "age"), data = d,
             tableOnly = TRUE, removeInf = TRUE)))

  expect_match(keep[[2]][keep[[1]] == "sep"], "Inf")
  expect_true(is.na(drop[[2]][drop[[1]] == "sep"]))
  expect_true(is.na(drop[["p-value"]][drop[[1]] == "sep"]))
  # the covariate is still listed, and other covariates are untouched
  expect_true("sep" %in% drop[[1]])
  expect_equal(keep[[2]][keep[[1]] == "age"], drop[[2]][drop[[1]] == "age"])
})


# Other options ---------------------------------------------------------------

test_that("character and ordered covariates are handled", {
  d <- sim_uv_data()
  expect_no_error(rm_uvsum(response = "y", covs = "chr", data = d, tableOnly = TRUE))
  expect_no_error(rm_uvsum(response = "y", covs = "ordy", data = d, tableOnly = TRUE))
})

test_that("a tibble gives the same table as a data frame", {
  skip_if_not_installed("dplyr")
  d <- sim_uv_data()
  a <- rm_uvsum(response = "y", covs = c("age", "sex"), data = d, tableOnly = TRUE)
  b <- rm_uvsum(response = "y", covs = c("age", "sex"),
                data = dplyr::as_tibble(d), tableOnly = TRUE)
  expect_equal(a, b, ignore_attr = TRUE)
})

test_that("missing response values are dropped with a warning", {
  d <- sim_uv_data()
  d$y[1:20] <- NA
  expect_warning(
    tab <- suppressMessages(
      rm_uvsum(response = "y", covs = "age", data = d, tableOnly = TRUE)),
    "missing response")
  expect_equal(tab$N[1], nrow(d) - 20)
})

test_that("p.adjust adjusts and unformattedp stays numeric", {
  d <- sim_uv_data()
  a <- rm_uvsum(response = "y", covs = c("age", "sex", "grp"), data = d,
                tableOnly = TRUE, unformattedp = TRUE)
  b <- rm_uvsum(response = "y", covs = c("age", "sex", "grp"), data = d,
                tableOnly = TRUE, unformattedp = TRUE, p.adjust = "bonferroni")
  expect_type(a[["p-value"]], "double")
  expect_true(all(b[["p-value"]] >= a[["p-value"]], na.rm = TRUE))
})

test_that("showN and showEvent drop their columns", {
  d <- sim_uv_data()
  tab <- rm_uvsum(response = "bin", covs = "age", data = d, tableOnly = TRUE,
                  showN = FALSE, showEvent = FALSE)
  expect_false("N" %in% names(tab))
  expect_false("Event" %in% names(tab))
})

test_that("returnModels returns one model per covariate", {
  d <- sim_uv_data()
  m <- rm_uvsum(response = "y", covs = c("age", "sex"), data = d, returnModels = TRUE)
  expect_named(m, c("age", "sex"))
})

test_that("dimchk attribute matches the returned table", {
  d <- sim_uv_data()
  tab <- rm_uvsum(response = "y", covs = c("age", "grp"), data = d, tableOnly = TRUE)
  expect_equal(attr(tab, "dimchk"), dim(tab))
  expect_equal(attr(tab, "to_indent"), which(!tab[[1]] %in% c("age", "grp")))
})
