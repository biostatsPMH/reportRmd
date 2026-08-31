# Shared simulated data ------------------------------------------------------

sim_sum_data <- function(n = 300, seed = 21) {
  set.seed(seed)
  d <- data.frame(
    age = rnorm(n, 60, 10),
    bmi = rnorm(n, 27, 4),
    sex = factor(sample(c("M", "F"), n, TRUE)),
    grp = factor(sample(c("A", "B", "C"), n, TRUE)),
    chr = sample(c("x", "y"), n, TRUE),
    lgl = sample(c(TRUE, FALSE), n, TRUE),
    cnt = rpois(n, 3),
    my_var.x = rnorm(n),
    stringsAsFactors = FALSE
  )
  d$miss <- d$age
  d$miss[1:40] <- NA
  d$catNA <- factor(sample(c("p", "q", NA), n, TRUE, c(.45, .45, .10)))
  d$one <- factor("only")
  d$allna <- NA
  d$unusedlv <- factor(sample(c("a", "b"), n, TRUE), levels = c("a", "b", "never"))
  d$dt <- as.Date("2020-01-01") + sample(1:900, n, TRUE)
  d
}


# ---------------------------------------------------------------------------
# rm_covsum
# ---------------------------------------------------------------------------

test_that("summary statistics match the data", {
  d <- sim_sum_data()
  tab <- rm_covsum(data = d, covs = c("age", "sex"), tableOnly = TRUE)
  expect_match(tab[[2]][tab[[1]] == "Mean (sd)"][1],
               sprintf("^%s ", format(round(mean(d$age), 1), nsmall = 1)))
  n_m <- sum(d$sex == "M")
  expect_match(tab[[2]][tab[[1]] == "M"],
               paste0("^", n_m, " \\(", round(100 * n_m / nrow(d)), "\\)"))
})

test_that("group columns carry the right denominators", {
  d <- sim_sum_data()
  tab <- rm_covsum(data = d, maincov = "sex", covs = "grp", tableOnly = TRUE)
  expect_true(any(grepl(paste0("F \\(n=", sum(d$sex == "F"), "\\)"), names(tab))))
  expect_true(any(grepl(paste0("M \\(n=", sum(d$sex == "M"), "\\)"), names(tab))))
})

test_that("the grouping variable is dropped from covs", {
  d <- sim_sum_data()
  expect_warning(
    tab <- rm_covsum(data = d, maincov = "sex", covs = c("sex", "age"), tableOnly = TRUE),
    "grouping variable")
  expect_false("sex" %in% tab[[1]])
  expect_true("age" %in% tab[[1]])
})

test_that("dropping every covariate gives a clear error", {
  d <- sim_sum_data()
  expect_error(
    suppressWarnings(rm_covsum(data = d, maincov = "sex", covs = "sex", tableOnly = TRUE)),
    "No covariates remain")
})

test_that("excludeLevels leaves a non-numeric p-value entry alone", {
  d <- sim_sum_data()
  # "excl" is not a p-value; coercing it used to raise a coercion warning
  expect_no_warning(
    tab <- rm_covsum(data = d, maincov = "sex", covs = "unusedlv",
                     excludeLevels = list(unusedlv = "b"), tableOnly = TRUE))
  expect_true("excl" %in% tab[["p-value"]])
})

test_that("excludeLevels survives a p-value adjustment", {
  d <- sim_sum_data()
  expect_no_warning(
    tab <- rm_covsum(data = d, maincov = "sex", covs = c("unusedlv", "age", "bmi"),
                     excludeLevels = list(unusedlv = "b"),
                     p.adjust = "bonferroni", tableOnly = TRUE))
  expect_true("excl" %in% tab[["p-value"]])
  # real p-values are still adjusted
  raw <- rm_covsum(data = d, maincov = "sex", covs = c("age", "bmi"),
                   tableOnly = TRUE, unformattedp = TRUE)
  adj <- rm_covsum(data = d, maincov = "sex", covs = c("age", "bmi"),
                   p.adjust = "bonferroni", tableOnly = TRUE, unformattedp = TRUE)
  expect_true(all(as.numeric(adj[["p-value"]]) >= as.numeric(raw[["p-value"]]),
                  na.rm = TRUE))
})

test_that("include_missing controls NA levels of the grouping variable only", {
  d <- sim_sum_data()
  d$sex[1:30] <- NA
  drop <- rm_covsum(data = d, maincov = "sex", covs = "age", tableOnly = TRUE)
  keep <- rm_covsum(data = d, maincov = "sex", covs = "age",
                    include_missing = TRUE, tableOnly = TRUE)
  expect_true(any(grepl("n=270", names(drop))))
  expect_gt(ncol(keep), ncol(drop))
})

test_that("covariate missingness is always reported", {
  d <- sim_sum_data()
  tab <- rm_covsum(data = d, maincov = "sex", covs = "miss", tableOnly = TRUE)
  expect_true("Missing" %in% tab[[1]])
  expect_equal(as.numeric(tab[[2]][tab[[1]] == "Missing"]), sum(is.na(d$miss)))
})

test_that("character, logical and single-level covariates work", {
  d <- sim_sum_data()
  expect_no_error(rm_covsum(data = d, maincov = "sex", covs = "chr", tableOnly = TRUE))
  expect_no_error(rm_covsum(data = d, maincov = "sex", covs = "lgl", tableOnly = TRUE))
  expect_no_error(rm_covsum(data = d, maincov = "sex", covs = "one", tableOnly = TRUE))
})

test_that("dropLevels controls empty factor levels", {
  d <- sim_sum_data()
  keep <- rm_covsum(data = d, maincov = "sex", covs = "unusedlv",
                    dropLevels = FALSE, tableOnly = TRUE)
  drop <- rm_covsum(data = d, maincov = "sex", covs = "unusedlv",
                    dropLevels = TRUE, tableOnly = TRUE)
  expect_true("never" %in% keep[[1]])
  expect_false("never" %in% drop[[1]])
})

test_that("display options behave", {
  d <- sim_sum_data()
  expect_false("p-value" %in%
    names(rm_covsum(data = d, maincov = "sex", covs = "age",
                    pvalue = FALSE, tableOnly = TRUE)))
  expect_false(any(grepl("Full Sample",
    names(rm_covsum(data = d, maincov = "sex", covs = "age",
                    full = FALSE, tableOnly = TRUE)))))
  expect_true("StatTest" %in%
    names(rm_covsum(data = d, maincov = "sex", covs = "age",
                    show.tests = TRUE, tableOnly = TRUE)))
  expect_true("Effect Size" %in%
    names(rm_covsum(data = d, maincov = "sex", covs = "age",
                    effSize = TRUE, tableOnly = TRUE)))
})

test_that("nicenames replaces dots and underscores", {
  d <- sim_sum_data()
  tab <- rm_covsum(data = d, maincov = "sex", covs = "my_var.x", tableOnly = TRUE)
  expect_true("my var x" %in% tab[[1]])
})

test_that("a tibble gives the same table as a data frame", {
  skip_if_not_installed("dplyr")
  d <- sim_sum_data()
  a <- rm_covsum(data = d, maincov = "sex", covs = c("age", "grp"), tableOnly = TRUE)
  b <- rm_covsum(data = dplyr::as_tibble(d), maincov = "sex",
                 covs = c("age", "grp"), tableOnly = TRUE)
  expect_equal(a, b, ignore_attr = TRUE)
})


# ---------------------------------------------------------------------------
# rm_compactsum
# ---------------------------------------------------------------------------

test_that("the grouping variable is dropped from xvars", {
  d <- sim_sum_data()
  # the reduced list used not to reach the summarising loop, so the variable
  # was still rendered, with empty group cells
  expect_warning(
    tab <- rm_compactsum(data = d, xvars = c("sex", "age"), grp = "sex", tableOnly = TRUE),
    "grouping variable")
  expect_equal(nrow(tab), 1L)
  expect_false(any(grepl("^sex", tab[[1]])))
  expect_false(anyNA(tab[[3]]))
})

test_that("all-missing variables are dropped from xvars", {
  d <- sim_sum_data()
  expect_warning(
    tab <- rm_compactsum(data = d, xvars = c("age", "allna"), grp = "sex", tableOnly = TRUE),
    "only missing data")
  expect_false("allna" %in% tab[[1]])
  expect_equal(nrow(tab), 1L)
})

test_that("dropping every variable gives a clear error", {
  d <- sim_sum_data()
  expect_error(
    suppressWarnings(rm_compactsum(data = d, xvars = "sex", grp = "sex", tableOnly = TRUE)),
    "No variables remain")
})

test_that("a single-level factor does not raise an internal warning", {
  d <- sim_sum_data()
  # a 1 x k contingency table drops to a vector; min(dim(NULL)) used to warn
  expect_no_warning(
    rm_compactsum(data = d, xvars = c("age", "one"), grp = "sex", tableOnly = TRUE))
})

test_that("binary variables are collapsed to a single row", {
  d <- sim_sum_data()
  tab <- rm_compactsum(data = d, xvars = "sex", tableOnly = TRUE)
  expect_equal(nrow(tab), 1L)
  expect_equal(tab[[1]], "sex - M")
  expect_match(tab[[2]], paste0("^", sum(d$sex == "M"), " \\("))
})

test_that("use_mean switches between mean and median summaries", {
  d <- sim_sum_data()
  med <- rm_compactsum(data = d, xvars = "age", grp = "sex", tableOnly = TRUE)
  mn <- rm_compactsum(data = d, xvars = "age", grp = "sex",
                      use_mean = TRUE, tableOnly = TRUE)
  expect_match(med[[2]], "-")                                   # median (Q1-Q3)
  expect_match(mn[[2]], sprintf("^%s ", format(round(mean(d$age), 1), nsmall = 1)))

  # a character vector selects which variables use the mean
  mix <- rm_compactsum(data = d, xvars = c("age", "bmi"), grp = "sex",
                       use_mean = "age", tableOnly = TRUE)
  bmi_med <- rm_compactsum(data = d, xvars = "bmi", grp = "sex", tableOnly = TRUE)
  expect_equal(mix[[2]][1], mn[[2]][1])
  expect_equal(mix[[2]][2], bmi_med[[2]][1])
})

test_that("use_mean rejects variables that are not in xvars", {
  d <- sim_sum_data()
  expect_error(
    rm_compactsum(data = d, xvars = "age", grp = "sex", use_mean = "bmi",
                  tableOnly = TRUE),
    "use_mean")
})

test_that("a missing grouping value is dropped with a message", {
  d <- sim_sum_data()
  d$sex[1:20] <- NA
  expect_message(
    rm_compactsum(data = d, xvars = "age", grp = "sex", tableOnly = TRUE),
    "observations missing sex")
})

test_that("missing counts appear in their own column", {
  d <- sim_sum_data()
  tab <- rm_compactsum(data = d, xvars = "miss", grp = "sex", tableOnly = TRUE)
  expect_true("Missing" %in% names(tab))
  expect_equal(as.numeric(tab[["Missing"]][1]), sum(is.na(d$miss)))
})

test_that("display options behave", {
  d <- sim_sum_data()
  expect_true("pTest" %in%
    names(rm_compactsum(data = d, xvars = "age", grp = "sex",
                        show.tests = TRUE, tableOnly = TRUE)))
  expect_true(any(grepl("Effect Size",
    names(rm_compactsum(data = d, xvars = "age", grp = "sex",
                        effSize = TRUE, tableOnly = TRUE)))))
  expect_equal(nrow(rm_compactsum(data = d, xvars = "age", grp = "sex",
                                  all.stats = TRUE, tableOnly = TRUE)), 4L)
  expect_type(rm_compactsum(data = d, xvars = "age", grp = "sex",
                            unformattedp = TRUE, tableOnly = TRUE)[["p-value"]],
              "double")
})

test_that("dates are summarised but not tested", {
  d <- sim_sum_data()
  expect_message(
    tab <- rm_compactsum(data = d, xvars = c("age", "dt"), grp = "sex", tableOnly = TRUE),
    "date variables")
  expect_true(is.na(tab[["p-value"]][tab[[1]] == "dt"]))
})

test_that("p.adjust adjusts upward", {
  d <- sim_sum_data()
  raw <- rm_compactsum(data = d, xvars = c("age", "bmi", "cnt"), grp = "sex",
                       tableOnly = TRUE, unformattedp = TRUE)
  adj <- rm_compactsum(data = d, xvars = c("age", "bmi", "cnt"), grp = "sex",
                       p.adjust = "bonferroni", tableOnly = TRUE, unformattedp = TRUE)
  expect_true(all(adj[["p-value"]] >= raw[["p-value"]], na.rm = TRUE))
})

test_that("a tibble gives the same table as a data frame", {
  skip_if_not_installed("dplyr")
  d <- sim_sum_data()
  a <- rm_compactsum(data = d, xvars = c("age", "grp"), grp = "sex", tableOnly = TRUE)
  b <- rm_compactsum(data = dplyr::as_tibble(d), xvars = c("age", "grp"),
                     grp = "sex", tableOnly = TRUE)
  expect_equal(a, b, ignore_attr = TRUE)
})


# Regression: variable-selection styles ---------------------------------------

test_that("column names held in a variable do not raise a deprecation warning", {
  d <- sim_sum_data()
  v <- c("age", "grp")
  g <- "sex"
  expect_no_warning(rm_covsum(data = d, maincov = g, covs = v, tableOnly = TRUE))
  expect_no_warning(rm_compactsum(data = d, xvars = v, grp = g, tableOnly = TRUE))
  # literal strings, bare names and tidyselect helpers all still work
  expect_no_error(rm_compactsum(data = d, xvars = c(age, grp), grp = sex, tableOnly = TRUE))
  expect_no_error(rm_covsum(data = d, maincov = "sex",
                            covs = dplyr::starts_with("b"), tableOnly = TRUE))
})

test_that("an unknown column name is reported by name", {
  d <- sim_sum_data()
  expect_error(rm_covsum(data = d, maincov = "sex", covs = c("age", "nope"),
                         tableOnly = TRUE),
               "nope")
})

test_that("two-level ordered factors are treated as binary", {
  d <- sim_sum_data()
  d$ord2 <- factor(sample(c("Low", "High"), nrow(d), TRUE),
                   levels = c("Low", "High"), ordered = TRUE)
  # class() is c("ordered", "factor"); the old test warned and missed this
  expect_no_warning(
    tab <- rm_compactsum(data = d, xvars = "ord2", grp = "sex", tableOnly = TRUE))
  expect_equal(nrow(tab), 1L)
  expect_equal(tab[[1]], "ord2 - High")
})


# Reproducible categorical tests ---------------------------------------------

sim_mc_data <- function(n = 400, seed = 3) {
  set.seed(seed)
  data.frame(
    g = factor(sample(LETTERS[1:6], n, TRUE)),
    v = factor(sample(letters[1:12], n, TRUE))
  )
}

pv <- function(tab) as.numeric(stats::na.omit(tab[["p-value"]]))

test_that("a Monte Carlo p-value is identical across calls", {
  # a 12 x 6 table defeats the exact test, so both functions fall back to
  # simulation; without a fixed seed the p-value moved on every call
  d <- sim_mc_data()
  cs <- replicate(3, pv(suppressMessages(
    rm_covsum(data = d, maincov = "g", covs = "v",
              tableOnly = TRUE, unformattedp = TRUE))))
  ks <- replicate(3, pv(suppressMessages(
    rm_compactsum(data = d, xvars = "v", grp = "g",
                  tableOnly = TRUE, unformattedp = TRUE))))
  expect_equal(length(unique(cs)), 1L)
  expect_equal(length(unique(ks)), 1L)
})

test_that("rm_covsum and rm_compactsum give the same Monte Carlo p-value", {
  d <- sim_mc_data()
  cs <- pv(suppressMessages(rm_covsum(data = d, maincov = "g", covs = "v",
                                      tableOnly = TRUE, unformattedp = TRUE)))
  ks <- pv(suppressMessages(rm_compactsum(data = d, xvars = "v", grp = "g",
                                          tableOnly = TRUE, unformattedp = TRUE)))
  expect_equal(cs, ks)
})

test_that("seed and B are honoured", {
  d <- sim_mc_data()
  a <- pv(suppressMessages(rm_covsum(data = d, maincov = "g", covs = "v",
                                     tableOnly = TRUE, unformattedp = TRUE, seed = 1)))
  b <- pv(suppressMessages(rm_covsum(data = d, maincov = "g", covs = "v",
                                     tableOnly = TRUE, unformattedp = TRUE, seed = 99)))
  expect_false(isTRUE(all.equal(a, b)))

  # B is reported and reaches fisher.test: 1/(B+1) is the smallest attainable
  # p-value, so the value must be a multiple of 1/(B+1)
  small <- suppressMessages(rm_covsum(data = d, maincov = "g", covs = "v",
                                      tableOnly = TRUE, unformattedp = TRUE, B = 200))
  expect_equal(round(pv(small) * 201) / 201, pv(small), tolerance = 1e-8)
})

test_that("seed = NULL restores the varying behaviour", {
  d <- sim_mc_data()
  out <- replicate(4, pv(suppressMessages(
    rm_covsum(data = d, maincov = "g", covs = "v",
              tableOnly = TRUE, unformattedp = TRUE, seed = NULL))))
  expect_gt(length(unique(out)), 1L)
})

test_that("the calling session's random state is left untouched", {
  d <- sim_mc_data()
  set.seed(42); before <- runif(3)
  set.seed(42); invisible(suppressMessages(
    rm_covsum(data = d, maincov = "g", covs = "v", tableOnly = TRUE)))
  expect_equal(runif(3), before)
})

test_that("the reported test names the number of replicates", {
  d <- sim_mc_data()
  cs <- suppressMessages(rm_covsum(data = d, maincov = "g", covs = "v",
                                   show.tests = TRUE, tableOnly = TRUE))
  ks <- suppressMessages(rm_compactsum(data = d, xvars = "v", grp = "g",
                                       show.tests = TRUE, tableOnly = TRUE))
  expect_match(cs[["StatTest"]][1], "B=1000", fixed = TRUE)
  expect_match(ks[["pTest"]][1], "B=1000", fixed = TRUE)
})


# Expected rather than observed counts ----------------------------------------

test_that("the exact test is triggered by expected, not observed, counts", {
  # 3 observed in a cell but all expected counts comfortably above 5
  set.seed(11)
  d <- data.frame(
    g = factor(rep(c("A", "B"), each = 60)),
    v = factor(c(rep(c("x", "y", "z"), times = c(3, 27, 30)),
                 rep(c("x", "y", "z"), times = c(20, 20, 20))))
  )
  expect_true(all(suppressWarnings(chisq.test(table(d$v, d$g))$expected) >= 5))
  cs <- suppressMessages(rm_covsum(data = d, maincov = "g", covs = "v",
                                   show.tests = TRUE, tableOnly = TRUE))
  expect_equal(cs[["StatTest"]][1], "Chi Sq")
})

test_that("rm_covsum and rm_compactsum choose the same test", {
  set.seed(101)
  n <- 220
  d <- data.frame(
    arm = factor(sample(c("Control", "Experimental"), n, TRUE)),
    ecog = factor(sample(0:3, n, TRUE, c(.35, .40, .20, .05))),
    site = factor(sample(paste0("Site ", LETTERS[1:8]), n, TRUE)),
    sex = factor(sample(c("F", "M"), n, TRUE))
  )
  cs <- suppressMessages(rm_covsum(data = d, maincov = "arm",
                                   covs = c("ecog", "site"),
                                   show.tests = TRUE, tableOnly = TRUE))
  ks <- suppressMessages(rm_compactsum(data = d, xvars = c("ecog", "site"),
                                       grp = "arm", show.tests = TRUE,
                                       tableOnly = TRUE))
  for (v in c("ecog", "site")) {
    expect_equal(
      as.numeric(cs[["p-value"]][cs[[1]] == v]),
      as.numeric(ks[["p-value"]][ks[[1]] == v]),
      tolerance = 1e-6
    )
  }
})
