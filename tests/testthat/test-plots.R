skip_if_not_installed("survival")
skip_if_not_installed("ggplot2")

sim_plot_data <- function(n = 250, seed = 41) {
  set.seed(seed)
  d <- data.frame(
    age = rnorm(n, 60, 10),
    sex = factor(sample(c("M", "F"), n, TRUE)),
    grp = factor(sample(c("A", "B", "C"), n, TRUE))
  )
  d$y <- 2 * d$age + 5 * (d$sex == "M") + rnorm(n, 0, 20)
  d$bin <- rbinom(n, 1, stats::plogis(scale(d$y)[, 1]))
  d$time <- rexp(n, 0.05)
  d$status <- rbinom(n, 1, 0.7)
  d$cmprsk <- sample(0:2, n, TRUE)
  d
}

quietly <- function(expr) suppressWarnings(suppressMessages(expr))


# ggkmcif2 --------------------------------------------------------------------
# These all aborted on ggplot2 < 3.5.0, because theme elements added in 3.5.0
# were set unconditionally. DESCRIPTION now declares a floor of 3.4.0 and the
# 3.5.0-only elements are applied only when they exist.

test_that("ggkmcif2 builds a plot", {
  d <- sim_plot_data()
  expect_s3_class(
    quietly(ggkmcif2(response = c("time", "status"), cov = "sex", data = d)),
    "ggplot")
})

test_that("ggkmcif2 works without a covariate", {
  d <- sim_plot_data()
  expect_s3_class(
    quietly(ggkmcif2(response = c("time", "status"), data = d)),
    "ggplot")
})

test_that("a numeric legend position works on any supported ggplot2", {
  d <- sim_plot_data()
  expect_s3_class(
    quietly(ggkmcif2(response = c("time", "status"), cov = "sex", data = d,
                     legend.pos = c(0.8, 0.8))),
    "ggplot")
})

test_that("the plot renders to a device", {
  d <- sim_plot_data()
  p <- quietly(ggkmcif2(response = c("time", "status"), cov = "sex", data = d))
  f <- tempfile(fileext = ".png")
  on.exit(unlink(f), add = TRUE)
  expect_no_error(quietly(ggplot2::ggsave(f, p, width = 6, height = 5)))
  expect_true(file.exists(f))
})

test_that("the type argument is case insensitive", {
  d <- sim_plot_data()
  for (ty in c("KM", "km", "CIF", "cif")) {
    expect_s3_class(
      quietly(ggkmcif2(response = c("time", "status"), cov = "sex",
                       data = d, type = ty)),
      "ggplot")
  }
})

test_that("an unrecognised type is still rejected", {
  d <- sim_plot_data()
  expect_error(
    quietly(ggkmcif2(response = c("time", "status"), cov = "sex",
                     data = d, type = "nonsense")),
    "KM or CIF")
})

test_that("competing risks data defaults to CIF", {
  d <- sim_plot_data()
  expect_s3_class(
    quietly(ggkmcif2(response = c("time", "cmprsk"), cov = "sex", data = d)),
    "ggplot")
})

test_that("optional plot elements do not break the build", {
  d <- sim_plot_data()
  expect_s3_class(
    quietly(ggkmcif2(response = c("time", "status"), cov = "sex", data = d,
                     table = FALSE)),
    "ggplot")
  expect_s3_class(
    quietly(ggkmcif2(response = c("time", "status"), cov = "sex", data = d,
                     conf.curves = TRUE)),
    "ggplot")
})


# Forest plots ----------------------------------------------------------------

test_that("forestplotUV builds a plot", {
  d <- sim_plot_data()
  expect_s3_class(
    quietly(forestplotUV(response = "bin", covs = c("age", "sex"),
                         data = d, family = "binomial")),
    "ggplot")
})

test_that("forestplotMV builds a plot, with and without unadjusted estimates", {
  d <- sim_plot_data()
  m <- glm(bin ~ age + sex, data = d, family = binomial)
  expect_s3_class(quietly(forestplotMV(m)), "ggplot")
  expect_s3_class(quietly(forestplotMV(m, include_unadjusted = TRUE)), "ggplot")
})

test_that("forestplotMV handles a formula held in a variable", {
  d <- sim_plot_data()
  f <- bin ~ age + sex
  expect_s3_class(
    quietly(forestplotMV(glm(f, data = d, family = binomial),
                         include_unadjusted = TRUE)),
    "ggplot")
})

test_that("forestplotMV honours conf.level with unadjusted estimates", {
  # conf.level was passed to rm_uvsum as an unevaluated symbol and failed to
  # resolve in that function's frame
  d <- sim_plot_data()
  m <- glm(bin ~ age + sex, data = d, family = binomial)
  expect_s3_class(
    quietly(forestplotMV(m, include_unadjusted = TRUE, conf.level = 0.90)),
    "ggplot")
})

test_that("forestplotUVMV combines the two", {
  d <- sim_plot_data()
  uv <- quietly(forestplotUV(response = "bin", covs = c("age", "sex"),
                             data = d, family = "binomial"))
  mv <- quietly(forestplotMV(glm(bin ~ age + sex, data = d, family = binomial)))
  expect_s3_class(quietly(forestplotUVMV(uv, mv)), "ggplot")
})

test_that("plotuv builds a panel", {
  d <- sim_plot_data()
  expect_no_error(
    quietly(plotuv(response = "bin", covs = c("age", "sex"), data = d)))
})
