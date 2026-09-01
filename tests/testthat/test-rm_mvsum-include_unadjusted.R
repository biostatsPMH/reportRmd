# Shared simulated data ------------------------------------------------------

sim_mv_data <- function(n = 400, seed = 1) {
  set.seed(seed)
  d <- data.frame(
    age = rnorm(n, 60, 10),
    tmb = rexp(n, 1),
    bmi = rnorm(n, 27, 4),
    sex = factor(sample(c("M", "F"), n, TRUE)),
    grp = factor(sample(c("A", "B", "C"), n, TRUE))
  )
  d$y <- 2 * d$age + 5 * (d$sex == "M") + rnorm(n, 0, 20)
  d$bin <- rbinom(n, 1, stats::plogis(scale(d$y)[, 1]))
  d$cnt <- rpois(n, exp(1 + 0.01 * d$age))
  d$time <- rexp(n, 0.05)
  d$status <- rbinom(n, 1, 0.7)
  d
}

# Pull the point estimate out of a "1.23 (0.99, 1.50)" cell
est_of <- function(tab, variable, col) {
  x <- tab[[col]][tab[[1]] == variable]
  as.numeric(sub("^\\s*([-0-9.e]+).*", "\\1", x[1]))
}

unadj_col <- function(tab) grep("^Unadjusted .*CI", names(tab), value = TRUE)[1]
adj_col   <- function(tab) grep("^Adjusted .*CI",   names(tab), value = TRUE)[1]


# 1. The unadjusted column really is univariate -------------------------------

test_that("unadjusted estimates match separately fitted univariate models", {
  d <- sim_mv_data()
  tab <- rm_mvsum(lm(y ~ age + sex + grp, data = d),
                  tableOnly = TRUE, include_unadjusted = TRUE)

  expect_equal(est_of(tab, "age", unadj_col(tab)),
               unname(round(coef(lm(y ~ age, data = d))[2], 2)))
  expect_equal(est_of(tab, "M", unadj_col(tab)),
               unname(round(coef(lm(y ~ sex, data = d))[2], 2)))
  expect_equal(est_of(tab, "age", adj_col(tab)),
               unname(round(coef(lm(y ~ age + sex + grp, data = d))[2], 2)))
  # adjusted and unadjusted must not be identical here
  expect_false(isTRUE(all.equal(est_of(tab, "M", unadj_col(tab)),
                                est_of(tab, "M", adj_col(tab)))))
})

test_that("include_unadjusted adds columns but not rows", {
  d <- sim_mv_data()
  for (fm in list(y ~ age + sex + grp, y ~ age * sex, y ~ age + tmb + bmi)) {
    m <- lm(fm, data = d)
    a <- rm_mvsum(m, tableOnly = TRUE)
    b <- rm_mvsum(m, tableOnly = TRUE, include_unadjusted = TRUE)
    expect_equal(nrow(a), nrow(b), info = deparse(fm))
    expect_equal(a[[1]], b[[1]], info = deparse(fm))
  }
})

test_that("adjusted column is unchanged by include_unadjusted", {
  d <- sim_mv_data()
  m <- lm(y ~ age + sex + grp, data = d)
  a <- rm_mvsum(m, tableOnly = TRUE)
  b <- rm_mvsum(m, tableOnly = TRUE, include_unadjusted = TRUE)
  expect_equal(a[["Estimate(95%CI)"]], b[[adj_col(b)]])
  expect_equal(a[["p-value"]], b[["Adjusted p-value"]])
  expect_equal(a$N, b$N)
})


# 2. coxph: response is a Surv() object ---------------------------------------

test_that("coxph works without the deprecated data argument", {
  skip_if_not_installed("survival")
  d <- sim_mv_data()
  m <- survival::coxph(survival::Surv(time, status) ~ age + sex + grp, data = d)

  tab <- expect_no_error(
    rm_mvsum(m, tableOnly = TRUE, include_unadjusted = TRUE))
  expect_true(unadj_col(tab) %in% names(tab))
  expect_false(all(is.na(tab[[unadj_col(tab)]])))

  uv <- survival::coxph(survival::Surv(time, status) ~ age, data = d)
  expect_equal(est_of(tab, "age", unadj_col(tab)),
               unname(round(exp(coef(uv))[1], 2)))
})

test_that("coxph gives the same answer with and without the data argument", {
  skip_if_not_installed("survival")
  d <- sim_mv_data()
  m <- survival::coxph(survival::Surv(time, status) ~ age + sex, data = d)
  a <- rm_mvsum(m, tableOnly = TRUE, include_unadjusted = TRUE)
  b <- suppressWarnings(
    rm_mvsum(m, data = d, tableOnly = TRUE, include_unadjusted = TRUE))
  expect_equal(a, b, ignore_attr = TRUE)
})


# 3. CIwidth other than 0.95 --------------------------------------------------

test_that("CIwidth is honoured on both sides of the table", {
  d <- sim_mv_data()
  tab <- expect_no_error(
    rm_mvsum(lm(y ~ age + sex, data = d), tableOnly = TRUE,
             include_unadjusted = TRUE, CIwidth = 0.90))

  expect_match(unadj_col(tab), "90%CI", fixed = TRUE)
  expect_match(adj_col(tab), "90%CI", fixed = TRUE)

  ci90 <- round(confint(lm(y ~ age, data = d), level = 0.90)[2, ], 2)
  cell <- tab[[unadj_col(tab)]][tab[[1]] == "age"]
  expect_equal(as.numeric(regmatches(cell,
    gregexpr("[-0-9.]+", cell))[[1]][2:3]), unname(ci90))
})

test_that("unadjusted CI width differs from the 95% default", {
  d <- sim_mv_data()
  m <- lm(y ~ age + sex, data = d)
  t90 <- rm_mvsum(m, tableOnly = TRUE, include_unadjusted = TRUE, CIwidth = 0.90)
  t95 <- rm_mvsum(m, tableOnly = TRUE, include_unadjusted = TRUE, CIwidth = 0.95)
  expect_false(identical(t90[[unadj_col(t90)]], t95[[unadj_col(t95)]]))
})


# 4. How the formula was written should not matter ----------------------------

test_that("a formula held in a variable works", {
  d <- sim_mv_data()
  f <- y ~ age + sex
  a <- expect_no_error(
    rm_mvsum(lm(f, data = d), tableOnly = TRUE, include_unadjusted = TRUE))
  b <- rm_mvsum(lm(y ~ age + sex, data = d),
                tableOnly = TRUE, include_unadjusted = TRUE)
  expect_equal(a, b, ignore_attr = TRUE)
})

test_that("a dot formula works", {
  d <- sim_mv_data()[, c("y", "age", "sex")]
  a <- expect_no_error(
    rm_mvsum(lm(y ~ ., data = d), tableOnly = TRUE, include_unadjusted = TRUE))
  b <- rm_mvsum(lm(y ~ age + sex, data = d),
                tableOnly = TRUE, include_unadjusted = TRUE)
  expect_equal(a, b, ignore_attr = TRUE)
})

test_that("update()d models work", {
  d <- sim_mv_data()
  m <- update(lm(y ~ age, data = d), . ~ . + sex)
  expect_no_error(rm_mvsum(m, tableOnly = TRUE, include_unadjusted = TRUE))
})


# 5. Transformed predictors ---------------------------------------------------

test_that("transformed predictors get a univariate estimate on the same scale", {
  d <- sim_mv_data()
  tab <- expect_no_error(
    rm_mvsum(lm(y ~ log(tmb) + sex, data = d),
             tableOnly = TRUE, include_unadjusted = TRUE))

  expect_true("log(tmb)" %in% tab[[1]])
  expect_equal(est_of(tab, "log(tmb)", unadj_col(tab)),
               unname(round(coef(lm(y ~ log(tmb), data = d))[2], 2)))
})

test_that("I() terms get a univariate estimate", {
  d <- sim_mv_data()
  tab <- expect_no_error(rm_mvsum(lm(y ~ I(age^2) + sex, data = d),
                                  tableOnly = TRUE, include_unadjusted = TRUE))
  expect_equal(est_of(tab, "I(age^2)", unadj_col(tab)),
               unname(round(coef(lm(y ~ I(age^2), data = d))[2], 2)))
})

test_that("include_unadjusted does not add failures for matrix-valued terms", {
  # poly() breaks rm_mvsum in getVarLevels() regardless of include_unadjusted.
  # This is a separate, pre-existing issue; the test records that
  # include_unadjusted neither causes nor masks it.
  d <- sim_mv_data()
  m <- lm(y ~ poly(age, 2) + sex, data = d)
  base_fails <- inherits(try(rm_mvsum(m, tableOnly = TRUE), silent = TRUE),
                         "try-error")
  uv_fails <- inherits(
    try(rm_mvsum(m, tableOnly = TRUE, include_unadjusted = TRUE), silent = TRUE),
    "try-error")
  expect_equal(uv_fails, base_fails)
})


# 6. Covariates appearing only inside an interaction --------------------------

test_that("a continuous term used only in an interaction is still reported", {
  d <- sim_mv_data()
  tab <- rm_mvsum(lm(y ~ grp + age:sex, data = d),
                  tableOnly = TRUE, include_unadjusted = TRUE)

  expect_true("age" %in% tab[[1]])   # was silently dropped
  expect_true("sex" %in% tab[[1]])
  expect_equal(est_of(tab, "age", unadj_col(tab)),
               unname(round(coef(lm(y ~ age, data = d))[2], 2)))
  # appended rows carry no adjusted estimate
  expect_true(is.na(tab[[adj_col(tab)]][tab[[1]] == "age"][1]))
})

test_that("continuous x continuous interactions are recognised", {
  d <- sim_mv_data()
  tab <- rm_mvsum(lm(y ~ age + tmb + age:bmi, data = d),
                  tableOnly = TRUE, include_unadjusted = TRUE)

  expect_true("bmi" %in% tab[[1]])   # was never appended
  expect_equal(est_of(tab, "bmi", unadj_col(tab)),
               unname(round(coef(lm(y ~ bmi, data = d))[2], 2)))
})

test_that("a fully specified interaction adds no extra rows", {
  d <- sim_mv_data()
  m <- lm(y ~ age * sex, data = d)
  expect_equal(nrow(rm_mvsum(m, tableOnly = TRUE)),
               nrow(rm_mvsum(m, tableOnly = TRUE, include_unadjusted = TRUE)))
})


# 7. Arguments passed through to rm_uvsum -------------------------------------

test_that("digits applies to the unadjusted column", {
  d <- sim_mv_data()
  tab <- rm_mvsum(lm(y ~ age + sex, data = d), tableOnly = TRUE,
                  include_unadjusted = TRUE, digits = 4)
  n_dp <- function(x) nchar(sub("^[-0-9]+\\.", "", sub(" .*$", "", x)))
  expect_equal(n_dp(tab[[unadj_col(tab)]][tab[[1]] == "age"]),
               n_dp(tab[[adj_col(tab)]][tab[[1]] == "age"]))
})

test_that("whichp applies to both columns", {
  d <- sim_mv_data()
  tab <- rm_mvsum(lm(y ~ age + sex + grp, data = d), tableOnly = TRUE,
                  include_unadjusted = TRUE, whichp = "global")
  # global p-values sit on the header row, level rows are blank
  i <- which(tab[[1]] == "grp")
  expect_false(is.na(tab[["Unadjusted p-value"]][i]))
  expect_true(is.na(tab[["Unadjusted p-value"]][which(tab[[1]] == "B")]))
})

test_that("unformattedp returns numeric p-values in both columns", {
  d <- sim_mv_data()
  tab <- rm_mvsum(lm(y ~ age + sex, data = d), tableOnly = TRUE,
                  include_unadjusted = TRUE, unformattedp = TRUE)
  expect_type(tab[["Unadjusted p-value"]], "double")
  expect_type(tab[["Adjusted p-value"]], "double")
})

test_that("p.adjust affects only the adjusted column", {
  d <- sim_mv_data()
  m <- lm(y ~ age + sex + grp, data = d)
  a <- rm_mvsum(m, tableOnly = TRUE, include_unadjusted = TRUE,
                unformattedp = TRUE)
  b <- rm_mvsum(m, tableOnly = TRUE, include_unadjusted = TRUE,
                unformattedp = TRUE, p.adjust = "bonferroni")
  expect_equal(a[["Unadjusted p-value"]], b[["Unadjusted p-value"]])
  expect_true(all(b[["Adjusted p-value"]] >= a[["Adjusted p-value"]], na.rm = TRUE))
})


# 8. Joining the two tables ---------------------------------------------------

test_that("factors sharing level names are not confused", {
  d <- sim_mv_data()
  d$f1 <- factor(rep(c("Low", "High"), length.out = nrow(d)))
  d$f2 <- factor(rep(c("Low", "High"), each = 2, length.out = nrow(d)))
  d$y <- d$y + 10 * (d$f2 == "Low")

  tab <- rm_mvsum(lm(y ~ f1 + f2 + age, data = d),
                  tableOnly = TRUE, include_unadjusted = TRUE)
  low <- which(tab[[1]] == "Low")
  expect_length(low, 2)
  expect_equal(as.numeric(sub("^\\s*([-0-9.]+).*", "\\1",
                              tab[[unadj_col(tab)]][low[1]])),
               unname(round(coef(lm(y ~ f1, data = d))[2], 2)))
  expect_equal(as.numeric(sub("^\\s*([-0-9.]+).*", "\\1",
                              tab[[unadj_col(tab)]][low[2]])),
               unname(round(coef(lm(y ~ f2, data = d))[2], 2)))
})

test_that("factor levels containing a colon are handled", {
  d <- sim_mv_data()
  levels(d$grp) <- c("A:1", "B:2", "C:3")
  tab <- rm_mvsum(lm(y ~ age + grp, data = d),
                  tableOnly = TRUE, include_unadjusted = TRUE)
  expect_false(is.na(tab[[unadj_col(tab)]][tab[[1]] == "B:2"]))
})

test_that("a continuous covariate placed after a factor keeps its own estimate", {
  d <- sim_mv_data()
  tab <- rm_mvsum(lm(y ~ grp + sex + age, data = d),
                  tableOnly = TRUE, include_unadjusted = TRUE)
  expect_equal(est_of(tab, "age", unadj_col(tab)),
               unname(round(coef(lm(y ~ age, data = d))[2], 2)))
})


# 9. Missing data -------------------------------------------------------------

test_that("univariate models use the multivariable complete-case sample", {
  d <- sim_mv_data()
  d$tmb[1:100] <- NA
  m <- lm(y ~ age + tmb + sex, data = d)
  tab <- rm_mvsum(m, tableOnly = TRUE, include_unadjusted = TRUE)

  cc <- model.frame(m)
  expect_equal(est_of(tab, "age", unadj_col(tab)),
               unname(round(coef(lm(y ~ age, data = cc))[2], 2)))
  expect_equal(unique(tab$N[tab[[1]] == "age"]), nrow(cc))
})


# 10. Model families ----------------------------------------------------------

test_that("glm binomial reports odds ratios in both columns", {
  d <- sim_mv_data()
  tab <- rm_mvsum(glm(bin ~ age + sex, data = d, family = binomial),
                  tableOnly = TRUE, include_unadjusted = TRUE)
  expect_match(unadj_col(tab), "OR")
  expect_equal(est_of(tab, "age", unadj_col(tab)),
               unname(round(exp(coef(glm(bin ~ age, data = d,
                                         family = binomial))[2]), 2)))
})

test_that("glm poisson reports rate ratios in both columns", {
  d <- sim_mv_data()
  tab <- rm_mvsum(glm(cnt ~ age + sex, data = d, family = poisson),
                  tableOnly = TRUE, include_unadjusted = TRUE)
  expect_match(unadj_col(tab), "RR")
  expect_equal(est_of(tab, "age", unadj_col(tab)),
               unname(round(exp(coef(glm(cnt ~ age, data = d,
                                         family = poisson))[2]), 2)))
})


# 11. Display options ---------------------------------------------------------

test_that("showN and showEvent still drop their columns", {
  d <- sim_mv_data()
  tab <- rm_mvsum(glm(bin ~ age + sex, data = d, family = binomial),
                  tableOnly = TRUE, include_unadjusted = TRUE,
                  showN = FALSE, showEvent = FALSE)
  expect_false("N" %in% names(tab))
  expect_false("Event" %in% names(tab))
})

test_that("column order is covariate, unadjusted pair, adjusted pair", {
  d <- sim_mv_data()
  tab <- rm_mvsum(lm(y ~ age + sex, data = d),
                  tableOnly = TRUE, include_unadjusted = TRUE)
  expect_equal(names(tab)[1:5],
               c("Covariate", unadj_col(tab), "Unadjusted p-value",
                 adj_col(tab), "Adjusted p-value"))
})

test_that("indent and bold attributes cover appended rows", {
  d <- sim_mv_data()
  tab <- rm_mvsum(lm(y ~ grp + age:sex, data = d),
                  tableOnly = TRUE, include_unadjusted = TRUE)
  hdr <- which(tab[[1]] %in% c("grp", "age:sex", "age", "sex"))
  expect_equal(sort(setdiff(seq_len(nrow(tab)), attr(tab, "to_indent"))),
               sort(hdr))
  expect_true(all(hdr %in% attr(tab, "bold_cells")[, 1]))
  expect_equal(attr(tab, "dimchk"), dim(tab))
})

test_that("the formatted table renders", {
  d <- sim_mv_data()
  expect_no_error(rm_mvsum(lm(y ~ age + sex, data = d),
                           include_unadjusted = TRUE))
})
