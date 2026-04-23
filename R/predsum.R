
#' Summarise Performance of a Prediction Model
#'
#' Standardised reporting of binary prediction model performance. Returns a
#' data frame of metrics (point estimate + confidence interval) and optionally
#' renders it via \code{reportRmd::outTable()}.
#'
#' @param true Vector of true class labels. Must be 0/1 or a two-level factor.
#' @param predicted Vector of predicted values, either dichotomous (0/1) or
#'   numeric scores. If scores, \code{cutoff} must be supplied.
#' @param cutoff Numeric cut-off applied as \code{predicted >= cutoff} to
#'   define the positive class. Required when \code{predicted} is non-binary
#'   and ignored (with a message) otherwise.
#' @param positive The value of \code{true} treated as the positive class.
#'   Defaults to \code{1} when \code{true} is 0/1, or to the second level
#'   when \code{true} is a factor.
#' @param ci Confidence interval method: \code{"clopper-pearson"} (default;
#'   Delong is used for AUC) or \code{"bootstrap"} (percentile).
#' @param ci_width Nominal coverage, default \code{0.95}.
#' @param ci_separate logical. If \code{FALSE} (default) estimate and CI are
#'   combined in a single \code{Value (95\%CI)} column. If \code{TRUE}, the
#'   output has four columns: \code{Statistic}, \code{Estimate}, \code{lower},
#'   \code{upper}.
#' @param nboot Number of bootstrap replicates, default \code{2000}.
#' @param digits Number of digits for rounding, default \code{2}.
#' @param output One of \code{"epi"} (default), \code{"data-science"},
#'   \code{"all"}, or a character vector of statistic names. Recognised names
#'   (case-insensitive): \code{prevalence, sensitivity, specificity, ppv, npv,
#'   accuracy, auc, aucpr, f1, brier, tp, fp, tn, fn, n, events, positive,
#'   cutoff}. The aliases \code{recall} and \code{precision} are accepted and,
#'   when supplied, are reflected in the row labels. Preset names cannot be
#'   combined with individual statistics.
#' @param tableOnly logical. If \code{TRUE} the data frame is returned
#'   silently; if \code{FALSE} (default) the data frame is rendered via
#'   \code{reportRmd::outTable()} and returned invisibly.
#'
#' @details
#' Clopper-Pearson intervals are used for proportion-based metrics
#' (prevalence, sensitivity, specificity, PPV, NPV, accuracy). AUC uses the
#' Delong variance estimator on the logit scale. Brier score, F1 and AUC-PR
#' have no closed-form interval and are obtained by percentile bootstrap;
#' when \code{ci = "clopper-pearson"} this is done automatically and a
#' message is emitted.
#'
#' When \code{predicted} is dichotomous, AUC, AUC-PR, Brier and the cut-off row are
#' omitted.
#'
#' @references
#' Delong ER, Delong DM, Clarke-Pearson DL (1988). Comparing the areas under
#' two or more correlated receiver operating characteristic curves: a
#' nonparametric approach. \emph{Biometrics} 44(3):837-845.
#'
#' Clopper CJ, Pearson ES (1934). The use of confidence or fiducial limits
#' illustrated in the case of the binomial. \emph{Biometrika} 26(4):404-413.
#'
#' @examples
#' # Simulate true labels and model scores
#' set.seed(1)
#' n <- 300
#' true  <- rbinom(n, 1, 0.3)
#' score <- plogis(rnorm(n, mean = ifelse(true == 1, 0.9, -0.9)))
#'
#' # Epidemiology-style summary from probability scores
#' predsum(true, score, cutoff = 0.5, tableOnly = TRUE)
#'
#' # Data-science style (Recall, Precision, F1, Brier, AUC-PR)
#' predsum(true, score, cutoff = 0.5,
#'         output = "data-science", tableOnly = TRUE)
#'
#' # Custom selection with separate CI columns
#' predsum(true, score, cutoff = 0.5,
#'         output = c("sensitivity", "specificity", "auc"),
#'         ci_separate = TRUE, tableOnly = TRUE)
#'
#' # Dichotomous predictions: AUC, AUC-PR, Brier and cut-off are omitted
#' pred_class <- as.integer(score >= 0.5)
#' predsum(true, pred_class, tableOnly = TRUE)
#'
#' @return A data frame of performance statistics.
#' @export
predsum <- function(true,
                    predicted,
                    cutoff = NULL,
                    positive = NULL,
                    ci = c("clopper-pearson", "bootstrap"),
                    ci_width = 0.95,
                    ci_separate = FALSE,
                    nboot = 2000,
                    digits = 2,
                    output = "epi",
                    tableOnly = FALSE) {

  ci <- match.arg(ci)
  alpha <- 1 - ci_width


  # Validate and coerce true
  if (is.factor(true)) {
    if (nlevels(true) != 2) stop("`true` must have exactly two levels.")
    if (is.null(positive)) positive <- levels(true)[2]
    if (!positive %in% levels(true)) stop("`positive` not a level of `true`.")
    true_bin <- as.integer(true == positive)
  } else {
    u <- unique(stats::na.omit(true))
    if (!all(u %in% c(0, 1))) stop("`true` must be 0/1 or a two-level factor.")
    if (is.null(positive)) positive <- 1
    true_bin <- as.integer(true == positive)
  }

  # Determine predicted type and validate cutoff
  pred_num <- as.numeric(predicted)
  is_dichot <- all(stats::na.omit(pred_num) %in% c(0, 1))
  if (!is_dichot) {
    if (is.null(cutoff))
      stop("`cutoff` must be supplied when `predicted` is not dichotomous.")
    if (!is.numeric(cutoff) || length(cutoff) != 1L || !is.finite(cutoff))
      stop("`cutoff` must be a single finite numeric value.")
  }
  if (is_dichot && !is.null(cutoff)) {
    message("`predicted` is already dichotomous; `cutoff` ignored.")
    cutoff <- NULL
  }

  # Drop NAs
  keep <- !is.na(true_bin) & !is.na(pred_num)
  if (any(!keep)) message(sum(!keep), " row(s) with NA dropped.")
  true_bin <- true_bin[keep]
  pred_num <- pred_num[keep]
  if (length(true_bin) == 0L)
    stop("No non-missing observations remaining after NA removal.")

  # Derive class predictions
  if (is_dichot) {
    pred_class <- as.integer(pred_num)
    scores <- NULL
  } else {
    pred_class <- as.integer(pred_num >= cutoff)
    scores <- pred_num
  }
  n <- length(true_bin)

  # Point estimates via shared helper
  pe <- .compute_metrics(true_bin, pred_class, scores)
  TP <- pe$tp; FP <- pe$fp; TN <- pe$tn; FN <- pe$fn
  prev <- pe$prevalence; sens <- pe$sensitivity; spec <- pe$specificity
  ppv  <- pe$ppv; npv  <- pe$npv; acc <- pe$accuracy
  f1_val <- pe$f1; brier <- pe$brier
  auc_val <- pe$auc; aucpr_val <- pe$aucpr

  # Clopper-Pearson helper
  cp <- function(x, m) {
    if (is.na(x) || m == 0) return(c(NA_real_, NA_real_))
    lb <- if (x == 0) 0 else stats::qbeta(alpha / 2, x, m - x + 1)
    ub <- if (x == m) 1 else stats::qbeta(1 - alpha / 2, x + 1, m - x)
    c(lb, ub)
  }

  # Quantile helper for bootstrap percentile CIs
  qf <- function(v) stats::quantile(v, c(alpha / 2, 1 - alpha / 2),
                                    na.rm = TRUE, names = FALSE)

  # Compute CIs
  if (ci == "clopper-pearson") {
    prev_ci <- cp(sum(true_bin), n)
    sens_ci <- cp(TP, TP + FN)
    spec_ci <- cp(TN, TN + FP)
    ppv_ci  <- cp(TP, TP + FP)
    npv_ci  <- cp(TN, TN + FN)
    acc_ci  <- cp(TP + TN, n)
    auc_ci  <- if (!is.null(scores)) .delong_ci(scores, true_bin, alpha)
    else c(NA_real_, NA_real_)

    # Shared bootstrap resamples for F1, Brier, AUC-PR
    message("Bootstrapping CIs for Brier, F1 and AUC-PR (", nboot, " replicates).")
    boot_mat <- matrix(NA_real_, nrow = nboot, ncol = 3,
                       dimnames = list(NULL, c("f1", "brier", "aucpr")))
    for (b in seq_len(nboot)) {
      idx <- sample.int(n, n, replace = TRUE)
      sc  <- if (!is.null(scores)) scores[idx] else NULL
      m   <- .compute_metrics(true_bin[idx], pred_class[idx], sc)
      boot_mat[b, "f1"]    <- m$f1
      boot_mat[b, "brier"] <- m$brier
      boot_mat[b, "aucpr"] <- m$aucpr
    }
    f1_ci    <- qf(boot_mat[, "f1"])
    brier_ci <- if (!is.null(scores)) qf(boot_mat[, "brier"]) else c(NA_real_, NA_real_)
    aucpr_ci <- if (!is.null(scores)) qf(boot_mat[, "aucpr"]) else c(NA_real_, NA_real_)

  } else {
    # Full percentile bootstrap via shared helper
    cols <- c("prev","sens","spec","ppv","npv","acc","f1","brier","auc","aucpr")
    boot_mat <- matrix(NA_real_, nrow = nboot, ncol = length(cols),
                       dimnames = list(NULL, cols))
    for (b in seq_len(nboot)) {
      idx <- sample.int(n, n, replace = TRUE)
      sc  <- if (!is.null(scores)) scores[idx] else NULL
      m   <- .compute_metrics(true_bin[idx], pred_class[idx], sc)
      boot_mat[b, ] <- c(m$prevalence, m$sensitivity, m$specificity,
                         m$ppv, m$npv, m$accuracy, m$f1, m$brier,
                         m$auc, m$aucpr)
    }
    prev_ci  <- qf(boot_mat[, "prev"])
    sens_ci  <- qf(boot_mat[, "sens"])
    spec_ci  <- qf(boot_mat[, "spec"])
    ppv_ci   <- qf(boot_mat[, "ppv"])
    npv_ci   <- qf(boot_mat[, "npv"])
    acc_ci   <- qf(boot_mat[, "acc"])
    f1_ci    <- qf(boot_mat[, "f1"])
    brier_ci <- if (!is.null(scores)) qf(boot_mat[, "brier"]) else c(NA_real_, NA_real_)
    auc_ci   <- if (!is.null(scores)) qf(boot_mat[, "auc"])   else c(NA_real_, NA_real_)
    aucpr_ci <- if (!is.null(scores)) qf(boot_mat[, "aucpr"]) else c(NA_real_, NA_real_)
  }

  # Assemble metric registry
  metrics <- list(
    prevalence  = list(label = "Prevalence",  est = prev,      ci = prev_ci,  count = FALSE),
    sensitivity = list(label = "Sensitivity", est = sens,      ci = sens_ci,  count = FALSE),
    specificity = list(label = "Specificity", est = spec,      ci = spec_ci,  count = FALSE),
    ppv         = list(label = "PPV",         est = ppv,       ci = ppv_ci,   count = FALSE),
    npv         = list(label = "NPV",         est = npv,       ci = npv_ci,   count = FALSE),
    accuracy    = list(label = "Accuracy",    est = acc,       ci = acc_ci,   count = FALSE),
    auc         = list(label = "AUC",         est = auc_val,   ci = auc_ci,   count = FALSE),
    aucpr       = list(label = "AUC-PR",      est = aucpr_val, ci = aucpr_ci, count = FALSE),
    f1          = list(label = "F1",          est = f1_val,    ci = f1_ci,    count = FALSE),
    brier       = list(label = "Brier",       est = brier,     ci = brier_ci, count = FALSE),
    tp          = list(label = "TP",          est = TP,        ci = NULL,     count = TRUE),
    fp          = list(label = "FP",          est = FP,        ci = NULL,     count = TRUE),
    tn          = list(label = "TN",          est = TN,        ci = NULL,     count = TRUE),
    fn          = list(label = "FN",          est = FN,        ci = NULL,     count = TRUE),
    n           = list(label = "Sample Size", est = n,         ci = NULL,     count = TRUE),
    events      = list(label = "Events",      est = sum(true_bin), ci = NULL, count = TRUE),
    positive    = list(label = "Positive Category", est = positive, ci = NULL, count = TRUE),
    cutoff      = list(label = "Cut-off",     est = if (is.null(cutoff)) NA else cutoff,
                       ci = NULL, count = TRUE)
  )

  # Preset key sets
  epi_keys <- c("prevalence","sensitivity","specificity","ppv","npv","accuracy",
                "auc","n","events","positive","cutoff")
  ds_keys  <- c("prevalence","sensitivity","specificity","ppv","npv","accuracy",
                "aucpr","f1","brier","n","events","positive","cutoff")
  all_keys <- c("prevalence","sensitivity","specificity","ppv","npv","accuracy",
                "auc","aucpr","f1","brier","tp","fp","tn","fn",
                "n","events","positive","cutoff")
  preset_names <- c("epi", "data-science", "all")

  # Resolve output selection
  if (length(output) == 1L && output %in% preset_names) {
    keys <- switch(output,
                   "epi"          = epi_keys,
                   "data-science" = ds_keys,
                   "all"          = all_keys)
    if (output == "data-science") {
      metrics$sensitivity$label <- "Recall"
      metrics$ppv$label         <- "Precision"
    }
  } else {
    if (any(output %in% preset_names)) {
      stop("Preset names ('epi', 'data-science', 'all') cannot be combined ",
           "with individual statistics.")
    }
    req <- tolower(output)
    alias_map <- c(recall = "sensitivity", precision = "ppv")
    canon <- ifelse(req %in% names(alias_map), alias_map[req], req)
    bad <- setdiff(canon, names(metrics))
    if (length(bad)) stop("Unknown statistic(s): ", paste(bad, collapse = ", "))
    if ("recall" %in% req)    metrics$sensitivity$label <- "Recall"
    if ("precision" %in% req) metrics$ppv$label         <- "Precision"
    keys <- canon
  }

  # Drop AUC / AUC-PR / Brier / cut-off when predictions dichotomous
  if (is.null(scores)) keys <- setdiff(keys, c("auc", "aucpr", "brier", "cutoff"))

  # Format rows
  fmt_num <- function(x) formatC(round(x, digits), format = "f", digits = digits)

  rows <- lapply(keys, function(k) {
    m <- metrics[[k]]
    if (isTRUE(m$count)) {
      est_chr <- as.character(m$est)
      if (ci_separate) {
        data.frame(Statistic = m$label, Estimate = est_chr,
                   lower = NA_character_, upper = NA_character_,
                   stringsAsFactors = FALSE)
      } else {
        data.frame(Statistic = m$label, Value = est_chr,
                   stringsAsFactors = FALSE)
      }
    } else {
      est_chr <- if (is.na(m$est)) NA_character_ else fmt_num(m$est)
      lb <- m$ci[1]; ub <- m$ci[2]
      if (ci_separate) {
        data.frame(Statistic = m$label,
                   Estimate  = est_chr,
                   lower     = if (is.na(lb)) NA_character_ else fmt_num(lb),
                   upper     = if (is.na(ub)) NA_character_ else fmt_num(ub),
                   stringsAsFactors = FALSE)
      } else {
        val <- if (is.na(m$est)) NA_character_
        else if (is.na(lb) || is.na(ub)) est_chr
        else paste0(est_chr, " (", fmt_num(lb), ", ", fmt_num(ub), ")")
        data.frame(Statistic = m$label, Value = val,
                   stringsAsFactors = FALSE)
      }
    }
  })

  df <- do.call(rbind, rows)
  if (!ci_separate) {
    names(df)[2] <- paste0("Value (", round(ci_width * 100), "%CI)")
  }

  # Render or return
  if (!tableOnly) {
    if (!requireNamespace("reportRmd", quietly = TRUE)) {
      stop("Package 'reportRmd' required when tableOnly = FALSE.")
    }
    print(reportRmd::outTable(df))
    return(invisible(df))
  }
  df
}


# Core metric computation shared by point estimates and bootstrap
.compute_metrics <- function(true_bin, pred_class, scores) {
  tp <- sum(pred_class == 1 & true_bin == 1)
  fp <- sum(pred_class == 1 & true_bin == 0)
  tn <- sum(pred_class == 0 & true_bin == 0)
  fn <- sum(pred_class == 0 & true_bin == 1)
  nn <- length(true_bin)
  prev  <- if (nn > 0) mean(true_bin) else NA_real_
  sens  <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  spec  <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
  ppv   <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
  npv   <- if ((tn + fn) > 0) tn / (tn + fn) else NA_real_
  acc   <- if (nn > 0) (tp + tn) / nn else NA_real_
  f1    <- if (!is.na(sens) && !is.na(ppv) && (sens + ppv) > 0)
    2 * sens * ppv / (sens + ppv) else NA_real_
  brier <- if (!is.null(scores)) mean((scores - true_bin)^2) else NA_real_
  auc   <- if (!is.null(scores)) .fast_auc(scores, true_bin) else NA_real_
  aucpr <- if (!is.null(scores)) .auc_pr(scores, true_bin)   else NA_real_
  list(tp = tp, fp = fp, tn = tn, fn = fn, n = nn,
       prevalence = prev, sensitivity = sens, specificity = spec,
       ppv = ppv, npv = npv, accuracy = acc, f1 = f1,
       brier = brier, auc = auc, aucpr = aucpr)
}


# Fast AUC via Wilcoxon rank-sum
.fast_auc <- function(scores, y) {
  x1 <- scores[y == 1]; x2 <- scores[y == 0]
  n1 <- length(x1); n2 <- length(x2)
  if (n1 == 0 || n2 == 0) return(NA_real_)
  r <- rank(c(x1, x2))
  (sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2) / (n1 * n2)
}


# Area under the precision-recall curve (trapezoidal)
.auc_pr <- function(scores, y) {
  P <- sum(y == 1)
  if (P == 0 || length(scores) == 0) return(NA_real_)
  ord <- order(scores, decreasing = TRUE)
  s <- scores[ord]; yy <- y[ord]
  tp <- cumsum(yy == 1)
  fp <- cumsum(yy == 0)
  is_last_of_tie <- c(s[-1] != s[-length(s)], TRUE)
  tp <- tp[is_last_of_tie]; fp <- fp[is_last_of_tie]
  recall    <- tp / P
  precision <- tp / (tp + fp)
  recall    <- c(0, recall)
  precision <- c(precision[1], precision)
  sum((recall[-1] - recall[-length(recall)]) *
        (precision[-1] + precision[-length(precision)]) / 2)
}


# Delong (1988) variance for a single AUC; logit-transformed CI
.delong_ci <- function(scores, y, alpha) {
  x1 <- scores[y == 1]; x0 <- scores[y == 0]
  n1 <- length(x1); n0 <- length(x0)
  if (n1 < 2 || n0 < 2) return(c(NA_real_, NA_real_))
  r_all <- rank(c(x1, x0))
  r_pos <- rank(x1); r_neg <- rank(x0)
  V10 <- (r_all[seq_len(n1)] - r_pos) / n0
  V01 <- 1 - (r_all[(n1 + 1):(n1 + n0)] - r_neg) / n1
  auc <- mean(V10)
  v <- stats::var(V10) / n1 + stats::var(V01) / n0
  if (!is.finite(v) || v <= 0) return(c(NA_real_, NA_real_))
  auc_c    <- min(max(auc, 1e-8), 1 - 1e-8)
  logit    <- log(auc_c / (1 - auc_c))
  se_logit <- sqrt(v) / (auc_c * (1 - auc_c))
  z  <- stats::qnorm(1 - alpha / 2)
  lb <- stats::plogis(logit - z * se_logit)
  ub <- stats::plogis(logit + z * se_logit)
  c(lb, ub)
}


