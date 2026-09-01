# Simulated analysis files with the kinds of quirks real study data carries
suppressMessages(library(reportRmd))

# 1. A tidy-ish oncology trial file -----------------------------------------
sim_trial <- function(n = 220, seed = 101) {
  set.seed(seed)
  d <- data.frame(
    patient_id = sprintf("PT-%03d", seq_len(n)),
    age_at_dx = round(rnorm(n, 63, 11), 1),
    sex = factor(sample(c("Female", "Male"), n, TRUE, c(.46, .54))),
    ecog = factor(sample(0:3, n, TRUE, c(.35, .40, .20, .05))),
    stage = factor(
      sample(c("I", "II", "III", "IV"), n, TRUE, c(.15, .25, .35, .25)),
      levels = c("I", "II", "III", "IV"),
      ordered = TRUE
    ),
    arm = factor(sample(c("Control", "Experimental"), n, TRUE)),
    site = factor(sample(paste0("Site ", LETTERS[1:8]), n, TRUE)),
    bmi = round(rnorm(n, 27, 5), 1),
    ldh = round(rlnorm(n, 5.4, .5)),
    smoker = sample(
      c("Never", "Former", "Current", NA),
      n,
      TRUE,
      c(.4, .35, .2, .05)
    ),
    stringsAsFactors = FALSE
  )
  # realistic missingness: biomarker measured in a subset only
  d$biomarker <- round(rlnorm(n, 1.2, .8), 2)
  d$biomarker[sample(n, round(.35 * n))] <- NA
  lp <- -3 +
    .03 * d$age_at_dx +
    .6 * (d$arm == "Experimental") +
    .5 * as.numeric(d$stage)
  d$response <- rbinom(n, 1, stats::plogis(lp - mean(lp)))
  d$os_time <- round(rexp(n, .02 + .004 * as.numeric(d$stage)), 1)
  d$os_status <- rbinom(n, 1, .65)
  d$pfs_months <- round(d$os_time * runif(n, .3, 1), 1)
  d
}

# 2. The same study as it usually arrives: types not yet cleaned -------------
sim_messy <- function(n = 220, seed = 102) {
  d <- sim_trial(n, seed)
  d$ecog <- as.character(d$ecog) # numeric-looking factor
  d$sex <- ifelse(d$sex == "Female", "F", "M") # character, not factor
  d$response <- ifelse(d$response == 1, "Yes", "No") # character outcome
  d$stage <- as.character(d$stage) # loses the ordering
  d$treatment_start_date <- as.Date("2019-01-01") + sample(0:1200, n, TRUE)
  d$`Prior lines of therapy` <- sample(0:4, n, TRUE)
  d$almost_constant <- factor(c(rep("No", n - 3), rep("Yes", 3)))
  d
}

# 3. A rare-disease file: small n, sparse cells ------------------------------
sim_small <- function(n = 34, seed = 103) {
  set.seed(seed)
  data.frame(
    age = round(rnorm(n, 55, 14), 1),
    sex = factor(sample(c("F", "M"), n, TRUE)),
    subtype = factor(sample(c("A", "B", "C"), n, TRUE, c(.6, .3, .1))),
    mutation = factor(sample(c("Wildtype", "Mutant"), n, TRUE, c(.85, .15))),
    resp = rbinom(n, 1, .4),
    time = round(rexp(n, .05), 1),
    status = rbinom(n, 1, .6)
  )
}

# 4. Repeated measures ------------------------------------------------------
sim_long <- function(n_pt = 60, k = 4, seed = 104) {
  set.seed(seed)
  d <- data.frame(
    id = rep(seq_len(n_pt), each = k),
    visit = factor(rep(paste0("V", seq_len(k)), n_pt)),
    arm = factor(rep(sample(c("A", "B"), n_pt, TRUE), each = k))
  )
  d$age <- rep(round(rnorm(n_pt, 60, 9), 1), each = k)
  d$score <- 50 +
    rep(rnorm(n_pt, 0, 8), each = k) +
    3 * as.numeric(d$visit) +
    4 * (d$arm == "B") +
    rnorm(nrow(d), 0, 5)
  d
}
