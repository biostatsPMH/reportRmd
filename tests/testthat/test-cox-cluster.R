test_that("clustered univariate coxph uses robust variances and reports clusters", {
  data("pembrolizumab")
  tab <- rm_uvsum(
    response = c("os_time", "os_status"),
    covs = c("age", "sex"),
    data = pembrolizumab,
    id = "cohort",
    tableOnly = TRUE
  )
  expect_true(any(grepl("clusters", names(tab))))
  m <- rm_uvsum(
    response = c("os_time", "os_status"),
    covs = "age",
    data = pembrolizumab,
    id = "cohort",
    returnModels = TRUE
  )
  expect_false(is.null(m$age$naive.var))
})

test_that("clustered coxph global p-values use robust Wald", {
  data("pembrolizumab")
  fit <- survival::coxph(
    survival::Surv(os_time, os_status) ~ age + cohort,
    data = pembrolizumab,
    cluster = id
  )
  g <- gp(fit)
  expect_equal(attr(g, "global_p"), "Robust Wald test")
  expect_false(is.na(g$global_p[g$var == "cohort"]))
})

test_that("include_unadjusted carries clustering and N columns", {
  data("pembrolizumab")
  fit <- survival::coxph(
    survival::Surv(os_time, os_status) ~ age + sex,
    data = pembrolizumab,
    cluster = cohort
  )
  tab <- rm_mvsum(fit, include_unadjusted = TRUE, tableOnly = TRUE)
  expect_true(any(grepl("^Unadjusted", names(tab))))
  expect_true(any(grepl("^N", names(tab))))
})
