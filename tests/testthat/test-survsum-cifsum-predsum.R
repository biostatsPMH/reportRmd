skip_if_not_installed("survival")

sim_tte_data <- function(n = 250, seed = 51) {
  set.seed(seed)
  d <- data.frame(
    age = rnorm(n, 60, 10),
    sex = factor(sample(c("F", "M"), n, TRUE)),
    grp = factor(sample(c("A", "B", "C"), n, TRUE))
  )
  d$time <- round(rexp(n, 0.05), 1)
  d$status <- rbinom(n, 1, 0.65)
  d$cmp <- sample(0:2, n, TRUE, c(.35, .40, .25))
  d$grpNA <- d$grp
  d$grpNA[1:20] <- NA
  d$truth <- rbinom(n, 1, stats::plogis(scale(d$age)[, 1]))
  d$pred <- stats::plogis(scale(d$age)[, 1] + rnorm(n, 0, 0.5))
  d
}

quietly <- function(expr) suppressWarnings(suppressMessages(expr))


# rm_survsum -----------------------------------------------------------------

test_that("rm_survsum returns one row per group plus the log-rank block", {
  d <- sim_tte_data()
  tab <- quietly(rm_survsum(data = d, time = "time", status = "status",
                            group = "grp", survtimes = c(12, 24),
                            survtimeunit = "month", tableOnly = TRUE))
  expect_true(all(levels(d$grp) %in% tab[[1]]))
  expect_true(any(vapply(tab, function(col) any(grepl("Log Rank", col)), logical(1))))
  expect_true(any(grepl("12month", names(tab))))
  expect_true(any(grepl("24month", names(tab))))

  ev <- tab[["Events/Total"]][tab[[1]] == "A"]
  expect_equal(ev, paste0(sum(d$status[d$grp == "A"]), "/", sum(d$grp == "A")))
})

test_that("rm_survsum works without a grouping variable", {
  d <- sim_tte_data()
  tab <- quietly(rm_survsum(data = d, time = "time", status = "status",
                            survtimes = 12, survtimeunit = "month",
                            tableOnly = TRUE))
  expect_equal(nrow(tab), 1L)
  expect_false("Group" %in% names(tab))
})

test_that("rm_survsum accepts a group name held in a variable", {
  d <- sim_tte_data()
  g <- "grp"
  a <- quietly(rm_survsum(data = d, time = "time", status = "status", group = g,
                          survtimes = 12, survtimeunit = "month", tableOnly = TRUE))
  b <- quietly(rm_survsum(data = d, time = "time", status = "status", group = "grp",
                          survtimes = 12, survtimeunit = "month", tableOnly = TRUE))
  expect_equal(a, b, ignore_attr = TRUE)
})

test_that("rm_survsum drops missing group values with a message", {
  d <- sim_tte_data()
  expect_message(
    rm_survsum(data = d, time = "time", status = "status", group = "grpNA",
               survtimes = 12, survtimeunit = "month", tableOnly = TRUE),
    "missing data were removed")
})

test_that("rm_survsum honours CIwidth", {
  d <- sim_tte_data()
  tab <- quietly(rm_survsum(data = d, time = "time", status = "status",
                            group = "grp", survtimes = 12, survtimeunit = "month",
                            CIwidth = 0.90, tableOnly = TRUE))
  expect_true(any(grepl("90% CI", names(tab))))
  expect_false(any(grepl("95% CI", names(tab))))
})

test_that("rm_survsum requires the time unit when survtimes are given", {
  d <- sim_tte_data()
  expect_error(
    rm_survsum(data = d, time = "time", status = "status", survtimes = 12,
               tableOnly = TRUE),
    "survtimeunit")
})


# rm_cifsum ------------------------------------------------------------------

test_that("rm_cifsum returns one row per stratum plus Gray's test", {
  skip_if_not_installed("cmprsk")
  d <- sim_tte_data()
  tab <- quietly(rm_cifsum(data = d, time = "time", status = "cmp",
                           group = "grp", eventtimes = c(12, 24),
                           eventtimeunit = "month", tableOnly = TRUE))
  expect_true(all(levels(d$grp) %in% tab[[1]]))
  expect_true(any(vapply(tab, function(col) any(grepl("Gray", col)), logical(1))))
  expect_true(any(grepl("12month", names(tab))))
})

test_that("rm_cifsum works without a grouping variable", {
  skip_if_not_installed("cmprsk")
  d <- sim_tte_data()
  tab <- quietly(rm_cifsum(data = d, time = "time", status = "cmp",
                           eventtimes = 12, eventtimeunit = "month",
                           tableOnly = TRUE))
  expect_equal(nrow(tab), 1L)
})

test_that("rm_cifsum accepts a group name held in a variable", {
  skip_if_not_installed("cmprsk")
  d <- sim_tte_data()
  g <- "grp"
  a <- quietly(rm_cifsum(data = d, time = "time", status = "cmp", group = g,
                         eventtimes = 12, eventtimeunit = "month", tableOnly = TRUE))
  b <- quietly(rm_cifsum(data = d, time = "time", status = "cmp", group = "grp",
                         eventtimes = 12, eventtimeunit = "month", tableOnly = TRUE))
  expect_equal(a, b, ignore_attr = TRUE)
})


# predsum --------------------------------------------------------------------

test_that("predsum reports the standard diagnostic statistics", {
  d <- sim_tte_data()
  tab <- quietly(predsum(true = d$truth, predicted = d$pred, cutoff = 0.5,
                         tableOnly = TRUE))
  for (s in c("Prevalence", "Sensitivity", "Specificity", "PPV", "NPV",
              "Accuracy", "AUC")) {
    expect_true(s %in% tab[[1]], info = s)
  }
  expect_equal(as.numeric(tab[[2]][tab[[1]] == "Sample Size"]), length(d$truth))
  expect_equal(as.numeric(tab[[2]][tab[[1]] == "Events"]), sum(d$truth))
})

test_that("predsum sensitivity and specificity match a hand-built table", {
  d <- sim_tte_data()
  tab <- quietly(predsum(true = d$truth, predicted = d$pred, cutoff = 0.5,
                         tableOnly = TRUE))
  cls <- as.integer(d$pred >= 0.5)
  sens <- sum(cls == 1 & d$truth == 1) / sum(d$truth == 1)
  spec <- sum(cls == 0 & d$truth == 0) / sum(d$truth == 0)
  num <- function(x) as.numeric(sub("^\\s*([0-9.]+).*", "\\1", x))
  expect_equal(num(tab[[2]][tab[[1]] == "Sensitivity"]), round(sens, 2))
  expect_equal(num(tab[[2]][tab[[1]] == "Specificity"]), round(spec, 2))
})

test_that("predsum requires a cutoff for continuous predictions", {
  d <- sim_tte_data()
  expect_error(predsum(true = d$truth, predicted = d$pred, tableOnly = TRUE),
               "cutoff")
})

test_that("predsum accepts a factor outcome with an explicit positive level", {
  d <- sim_tte_data()
  f <- factor(ifelse(d$truth == 1, "Yes", "No"))
  tab <- quietly(predsum(true = f, predicted = d$pred, cutoff = 0.5,
                         positive = "Yes", tableOnly = TRUE))
  expect_equal(as.numeric(tab[[2]][tab[[1]] == "Events"]), sum(f == "Yes"))
})


# predsum: reproducible bootstrap CIs ----------------------------------------

.pd <- sim_tte_data()   # built once: regenerating it would re-seed the RNG
boot_vals <- function(...) {
  quietly(predsum(true = .pd$truth, predicted = .pd$pred, cutoff = 0.5,
                  tableOnly = TRUE, nboot = 200,
                  output = c("brier", "f1", "aucpr"), ...))[[2]]
}

test_that("bootstrap intervals are identical across calls", {
  out <- replicate(3, paste(boot_vals(), collapse = "|"))
  expect_equal(length(unique(out)), 1L)
})

test_that("seed and seed = NULL behave as documented", {
  expect_false(identical(paste(boot_vals(seed = 1), collapse = "|"),
                         paste(boot_vals(seed = 99), collapse = "|")))
  out <- replicate(4, paste(boot_vals(seed = NULL), collapse = "|"))
  expect_gt(length(unique(out)), 1L)
})

test_that("predsum leaves the calling session's random state alone", {
  set.seed(42); before <- runif(3)
  set.seed(42); invisible(boot_vals())
  expect_equal(runif(3), before)
})
