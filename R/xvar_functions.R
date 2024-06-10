#' Create a summary table for an individual covariate
#'
#' @param xvar character with the name of covariate to include in table
#' @param data dataframe containing data
#' @param grp character with the name of the grouping variable
#' @param covTitle character with the name of the covariate (predictor) column.
#'   The default is to leave this empty for output or, for table only output to
#'   use the column name 'Covariate'
#' @param digits numeric specifying the number of digits for summarizing mean
#'   data. Otherwise, can specify for individual covariates using a vector of
#'   digits where each element is named using the covariate name. If a covariate
#'   is not in the vector the default will be used for it (default is 1).
#'   See examples
#' @param digits.cat numeric specifying the numer of digits for the proportions
#'   when summarizing categorical data (default is 0)
#' @param iqr logical indicating if you want to display the interquartile range
#'   (Q1, Q3) as opposed to (min, max) in the summary for continuous variables
#' @param all.stats logical indicating if all summary statistics (Q1, Q3 + min,
#'   max on a separate line) should be displayed. Overrides iqr
#' @param pvalue logical indicating if you want p-values included in the table
#' @param effSize logical indicating if you want effect sizes and their 95%
#'   confidence intervals included in the table. Can only be obtained if pvalue
#'   is also requested. Effect sizes calculated include Cramer's V for
#'   categorical variables, and Cohen's d, Wilcoxon r, Epsilon-squared, or
#'   Omega-squared for numeric/continuous variables
#' @param show.tests logical indicating if the type of statistical test and
#'   effect size (if effSize = TRUE) used should be shown in a column beside the
#'   p-values. Ignored if pvalue = FALSE
#' @param percentage choice of how percentages are presented, either column
#'   (default) or row
#' @return A data frame is returned
#'
#' @keywords internal

xvar_function <- function(xvar, data, grp, covTitle = "", digits = 1, digits.cat = 0, iqr = FALSE, all.stats = FALSE, pvalue = TRUE, effSize = FALSE, show.tests = FALSE, percentage = "col") {
  UseMethod("xvar_function", xvar)
}

#' Helper function: xvar_function.default
#'
#' @param ...
#'
#' @keywords internal
xvar_function.default <- function(xvar, ...) {
  stop("No default method for xvar_function. The xvar class must be known")
}

xvar_function.rm_binary <- function(xvar, data, grp, covTitle = "", digits = 1, digits.cat = 0, iqr = FALSE, all.stats = FALSE, pvalue = TRUE, effSize = FALSE, show.tests = FALSE, percentage = "col") {
  if (!pvalue) {
    show.tests = FALSE
  }
  ## *** why change the class?
  class(xvar) <- "character"
  df <- data.frame(Covariate = paste0(xvar, " n(%)"))
  if (covTitle == "") {
    colnames(df) <- " "
  }
  else {
    colnames(df) <- covTitle
  }
  x_var <- data[[xvar]]
  if (percentage == "row") {
    df[, paste0("Full Sample (n=", nrow(data), ")")] <- as.character(sum(x_var, na.rm = TRUE))
  }
  else {
    df[, paste0("Full Sample (n=", nrow(data), ")")] <- paste0(sum(x_var, na.rm = TRUE), " (", format(round((100*sum(x_var, na.rm = TRUE) / (nrow(data) - sum(is.na(x_var)))), digits.cat), nsmall = digits.cat), ")")
  }

  if (!missing(grp)) {
    group_var <- data[[grp]]
    grp <- as.factor(grp)
    n_levels <- nlevels(group_var)
    grp_levels <- levels(group_var)
    columns <- lapply(as.list(grp_levels), binary_xvar_helper, data = data, xvar = xvar, grp = grp, digits.cat = digits.cat, percentage = percentage)
    i = 1
    for (grp_level in grp_levels) {
      sub <- subset(data, group_var == grp_level)
      title <- paste0(grp_level, " (n=", nrow(sub), ")")
      df[1, title] <- columns[i]
      i = i + 1
    }
    cont_table <- table(x_var, group_var)
    cont_table <- cont_table[rowSums(cont_table) > 0, colSums(cont_table) > 0]
    if (is.null(dim(cont_table)) || min(dim(cont_table)) < 2) {
      df[1, "Missing"] <- sum(is.na(x_var), na.rm = TRUE)
      return(df)
    }
    if (!any(chi.test.rm(cont_table)$expected < 5)) {
      chisq_test <- chi.test.rm(cont_table)
      df[1, "p-value"] <- chisq_test$p.value
      if (effSize) {
        output <- calc_CramerV(chisq_test)
        df[1, "Effect Size (95% CI)"] <- psthr(output)
      }
      df[1, "Missing"] <- sum(is.na(x_var))
      if (show.tests) {
        df[1, "pTest"] <- "ChiSq"
        if (effSize) {
          df[1, "effStat"] <- "Cramer's V"
        }
      }
    }
    else if (any(chi.test.rm(cont_table)$expected < 5)) {
      fisher_test <- fisher.test.rm(cont_table)
      df[1, "p-value"] <-fisher_test$p.value
      if (effSize) {
        output <- calc_CramerV(fisher_test)
        df[1, "Effect Size (95% CI)"] <- psthr(output)
      }
      if (show.tests) {
        df[1, "pTest"] <- "Fisher Exact"
        if (effSize) {
          df[1, "effStat"] <- "Cramer's V"
        }
      }
    }
  }
  df[1, "Missing"] <- sum(is.na(x_var))
  return(df)
}

xvar_function.rm_mean <- function(xvar, data, grp, covTitle = "", digits = 1, digits.cat = 0, iqr = FALSE, all.stats = FALSE, pvalue = TRUE, effSize = FALSE, show.tests = FALSE, percentage = "col") {
  if (!pvalue) {
    show.tests = FALSE
  }
  class(xvar) <- "character"
  df <- data.frame(Covariate = paste0(xvar, " Mean (sd)"))
  if (covTitle == "") {
    colnames(df) <- " "
  }
  else {
    colnames(df) <- covTitle
  }
  x_var <- data[[xvar]]
  df[, paste0("Full Sample (n=", nrow(data), ")")] <- paste0(format(round(mean(x_var, na.rm = TRUE), digits), nsmall = digits), " (", format(round(sd(x_var, na.rm = TRUE), digits), nsmall = digits), ")")

  if (!missing(grp)) {
    group_var <- data[[grp]]
    grp_columns <- lapply(levels(group_var), mean_by_grp, data = data, xvar = xvar, grp = grp, digits = digits)
    i = 1
    for (level in levels(group_var)) {
      sub <- subset(data, group_var == level)
      title <- paste0(level, " (n=", nrow(sub), ")")
      df[, title] <- grp_columns[i]
      i = i+1
    }
    no_na_data <- na.omit(data.frame(x_var = data[[xvar]], group_var = data[[grp]]))
    if (length(unique(na.omit(group_var))) < 2 | length(unique(no_na_data$group_var)) < 2) {
      df[, "Missing"] <- sum(is.na(x_var))
      return(df)
    }

    if (length(levels(group_var)) == 2) {
      t_test <- t.test.rm(x_var, group_var)
      df[, "p-value"] <- t_test$p.value
      N <- nrow(data)
      if (effSize) {
        output <- calc_cohenD(t_test)
        df[, "Effect Size (95% CI)"] <- psthr(output)
      }
      df[, "Missing"] <- sum(is.na(x_var))
      if (show.tests) {
        df[, "pTest"] <- "t-test"
        if (effSize) {
          df[, "effStat"] <- "Cohen's d"
        }
      }
    }
    else if (length(levels(group_var)) > 2) {
      anova_test <- stats::aov(x_var ~ group_var)
      df[, "p-value"] <- summary(anova_test)[[1]][["Pr(>F)"]][1]
      if (effSize) {
        output <- calc_omegaSq(anova_test)
        df[, "Effect Size (95% CI)"] <- psthr(output)
      }
      df[, "Missing"] <- sum(is.na(x_var))
      if (show.tests) {
        df[, "pTest"] <- "ANOVA"
        if (effSize) {
          df[, "effStat"] <- "Omega Sq"
        }
      }
    }
  }
  df[, "Missing"] <- sum(is.na(x_var))
  return(df)
}

xvar_function.rm_median <- function(xvar, data, grp, covTitle = "", digits = 1, digits.cat = 0, iqr = FALSE, all.stats = FALSE, pvalue = TRUE, effSize = FALSE, show.tests = FALSE, percentage = "col") {
  if (!pvalue) {
    show.tests = FALSE
  }
  class(xvar) <- "character"
  x_var <- data[[xvar]]
  if (all.stats) {
    df <- data.frame(Covariate = c(xvar, "  Median (Q1, Q3)", "  Range (min, max)"))
    if (covTitle == "") {
      colnames(df) <- " "
    }
    else {
      colnames(df) <- covTitle
    }
    bracket_iqr <- paste0("(", format(round(stats::quantile(x_var, na.rm = TRUE, prob = 0.25), digits), nsmall = digits), ", ", format(round(stats::quantile(x_var, na.rm = TRUE, prob = 0.75), digits), nsmall = digits), ")")
    bracket_range <- paste0("(", format(round(min(x_var, na.rm = TRUE), digits), nsmall = digits), ", ", format(round(max(x_var, na.rm = TRUE), digits), nsmall = digits), ")")
    df[2, paste0("Full Sample (n=", nrow(data), ")")] <- paste0(format(round(median(x_var, na.rm = TRUE), digits), nsmall = digits), " ", bracket_iqr)
    df[3, paste0("Full Sample (n=", nrow(data), ")")] <- bracket_range
  }
  else {
    df <- data.frame(Covariate = paste0(xvar, ifelse(!iqr, " Median (Min, Max)", " Median (Q1, Q3)")))
    if (covTitle == "") {
      colnames(df) <- " "
    }
    else {
      colnames(df) <- covTitle
    }
    bracket <- ifelse(!iqr, paste0("(", format(round(min(x_var, na.rm = TRUE), digits), nsmall = digits), ", ", format(round(max(x_var, na.rm = TRUE), digits), nsmall = digits), ")"), paste0("(", format(round(stats::quantile(x_var, na.rm = TRUE, prob = 0.25), digits), nsmall = digits), ", ", format(round(stats::quantile(x_var, na.rm = TRUE, prob = 0.75), digits), nsmall = digits), ")"))
    df[, paste0("Full Sample (n=", nrow(data), ")")] <- paste0(format(round(median(x_var, na.rm = TRUE), digits), nsmall = digits), " ", bracket)
  }

  if (!missing(grp)) {
    group_var <- data[[grp]]
    if (all.stats) {
      grp_med <- lapply(as.list(levels(group_var)), median_by_grp, data = data, xvar = xvar, grp = grp, iqr = T, digits = digits)
      grp_range <- lapply(as.list(levels(group_var)), median_by_grp, data = data, xvar = xvar, grp = grp, iqr = F, digits = digits, range_only = T)
      i = 1
      for (level in levels(group_var)) {
        sub <- subset(data, group_var == level)
        title <- paste0(level, " (n=", nrow(sub), ")")
        df[2, title] <- grp_med[i]
        df[3, title] <- grp_range[i]
        i = i+1
      }
    }
    else {
      grp_columns <- lapply(as.list(levels(group_var)), median_by_grp, data = data, xvar = xvar, grp = grp, iqr = iqr, digits = digits)
      i = 1
      for (level in levels(group_var)) {
        sub <- subset(data, group_var == level)
        title <- paste0(level, " (n=", nrow(sub), ")")
        df[, title] <- grp_columns[i]
        i = i+1
      }
    }
    no_na_data <- na.omit(data.frame(x_var = data[[xvar]], group_var = data[[grp]]))
    if (length(unique(na.omit(group_var))) < 2 | length(unique(no_na_data$group_var)) < 2) {
      df[1, "Missing"] <- sum(is.na(x_var))
      return(df)
    }
    if (length(unique(group_var)) == 2) {
      wilcox_test <- wilcox.test.rm(x_var, group_var)
      df[1, "p-value"] <- wilcox_test$p.value
      if (effSize) {
        output <- calc_WilcoxonR(wilcox_test)
        df[1, "Effect Size (95% CI)"] <- psthr(output)
      }
      df[1, "Missing"] <- sum(is.na(x_var))
      if (show.tests) {
        df[1, "pTest"] <- "Wilcoxon Rank Sum"
        if (effSize) {
          df[1, "effStat"] <- "Wilcoxon r"
        }
      }
    }
    else if (length(levels(group_var)) > 2) {
      kruskal_test <- kruskal.test.rm(x_var, group_var)
      df[1, "p-value"] <- kruskal_test$p.value
      if (effSize) {
        output <- calc_epsilonSq(kruskal_test)
        df[1, "Effect Size (95% CI)"] <- psthr(output)
      }
      df[1, "Missing"] <- sum(is.na(x_var))
      if (show.tests) {
        df[1, "pTest"] <- "Kruskal Wallis"
        if (effSize) {
          df[1, "effStat"] <- "Epsilon sq"
        }
      }
    }
  }
  df[1, "Missing"] <- sum(is.na(x_var))
  return(df)
}

xvar_function.rm_categorical <- function(xvar, data, grp, covTitle = "", digits = 1, digits.cat = 0, iqr = FALSE, all.stats = FALSE, pvalue = TRUE, effSize = FALSE, show.tests = FALSE, percentage = "col") {
  if (!pvalue) {
    show.tests = FALSE
  }
  class(xvar) <- "character"
  x_var <- data[[xvar]]
  rows <- c(paste0(xvar, " n(%)"))
  for (xvar_level in levels(x_var)) {
    rows <- append(rows, xvar_level)
  }
  df <- data.frame(Covariate = rows)
  if (covTitle == "") {
    colnames(df) <- " "
  }
  else {
    colnames(df) <- covTitle
  }
  i = 2
  for (xvar_level in levels(x_var)) {
    xvar_subset <- subset(data, x_var == xvar_level)
    if (percentage == "row") {
      df[i, paste0("Full Sample (n=", nrow(data), ")")] <- as.character(nrow(xvar_subset))
    }
    else {
      df[i, paste0("Full Sample (n=", nrow(data), ")")] <- paste0(nrow(xvar_subset), " (", format(round((100*nrow(xvar_subset) / (nrow(data) - sum(is.na(x_var)))), digits.cat), nsmall = digits.cat), ")")
    }
    i = i + 1
  }
  if (!missing(grp)) {
    group_var <- data[[grp]]
    i = 2
    for (xvar_level in levels(x_var)) {
      columns <- categ_xvar_helper(xvar_level, data, xvar, grp, digits.cat, percentage)
      df[i, paste0("Full Sample (n=", nrow(data), ")")] <- columns[2, "Full Sample"]
      df[i, columns[1, levels(group_var)]] <- columns[2, levels(group_var)]
      i = i + 1
    }
    cont_table <- table(x_var, group_var)
    cont_table <- cont_table[rowSums(cont_table) > 0, colSums(cont_table) > 0]
    if (is.null(dim(cont_table)) | min(dim(cont_table)) < 2) {
      df[1, "Missing"] <- sum(is.na(x_var))
      return(df)
    }

    if (!any(chi.test.rm(cont_table)$expected < 5)) {
      chisq_test <- chi.test.rm(cont_table)
      df[1, "p-value"] <- chisq_test$p.value
      if (effSize) {
        output <- calc_CramerV(chisq_test)
        df[1, "Effect Size (95% CI)"] <- psthr(output)
      }
      df[1, "Missing"] <- sum(is.na(x_var))
      if (show.tests) {
        df[1, "pTest"] <- "ChiSq"
        if (effSize) {
          df[1, "effStat"] <- "Cramer's V"
          if (effSize) {
            df[1, "effStat"] <- "Cramer's V"
          }
        }
      }
    }

    else if (any(chi.test.rm(cont_table)$expected < 5)) {
      fisher_test <- fisher.test.rm(cont_table)
      df[1, "p-value"] <- fisher_test$p.value
      if (effSize) {
        output <- calc_CramerV(fisher_test)
        df[1, "Effect Size (95% CI)"] <- psthr(output)
      }
      if (show.tests) {
        df[1, "pTest"] <- "Fisher Exact"
        if (effSize) {
          df[1, "effStat"] <- "Cramer's V"
          if (effSize) {
            df[1, "effStat"] <- "Cramer's V"
          }
        }
      }
    }
  }
  df[1, "Missing"] <- sum(is.na(x_var))
  return(df)
}

xvar_function.rm_two_level <- function(xvar, data, grp, covTitle = "", digits = 1, digits.cat = 0, iqr = FALSE, all.stats = FALSE, pvalue = TRUE, effSize = FALSE, show.tests = FALSE, percentage = "col") {
  if (!pvalue) {
    show.tests = FALSE
  }
  class(xvar) <- "character"
  temp <- data.frame()
  x_var <- data[[xvar]]
  unique_levels <- unique(x_var)
  unique_levels <- sort(unique_levels)
  binary_column <- ifelse(x_var == unique_levels[1], 0, 1)

  if (!missing(grp)) {
    temp <- data[, grp]
    temp[[xvar]] <- binary_column
  }
  else {
    temp <- data[, xvar]
    temp[[xvar]] <- binary_column
  }
  df <- data.frame(Covariate = paste0(xvar, " n(%)"))
  if (covTitle == "") {
    colnames(df) <- " "
  }
  else {
    colnames(df) <- covTitle
  }
  x_var <- temp[[xvar]]
  if (percentage == "row") {
    df[, paste0("Full Sample (n=", nrow(temp), ")")] <- as.character(sum(x_var, na.rm = TRUE))
  }
  else {
    df[, paste0("Full Sample (n=", nrow(temp), ")")] <- paste0(sum(x_var, na.rm = TRUE), " (", format(round((100*sum(x_var, na.rm = TRUE) / (nrow(data) - sum(is.na(x_var)))), digits.cat), nsmall = digits.cat), ")")
  }

  if (!missing(grp)) {
    group_var <- temp[[grp]]
    grp <- as.factor(grp)
    n_levels <- nlevels(group_var)
    grp_levels <- levels(group_var)
    columns <- lapply(as.list(grp_levels), binary_xvar_helper, temp, xvar = xvar, grp = grp, digits.cat = digits.cat, percentage = percentage)
    i = 1
    for (grp_level in grp_levels) {
      sub <- subset(temp, group_var == grp_level)
      title <- paste0(grp_level, " (n=", nrow(sub), ")")
      df[1, title] <- columns[i]
      i = i + 1
    }
    cont_table <- table(x_var, group_var)
    cont_table <- cont_table[rowSums(cont_table) > 0, colSums(cont_table) > 0]
    if (is.null(dim(cont_table)) | min(dim(cont_table)) < 2) {
      df[1, "Missing"] <- sum(is.na(x_var))
      return(df)
    }
    if (!any(chi.test.rm(cont_table)$expected < 5)) {
      chisq_test <- chi.test.rm(cont_table)
      df[1, "p-value"] <- chisq_test$p.value
      if (effSize) {
        output <- calc_CramerV(chisq_test)
        df[1, "Effect Size (95% CI)"] <- psthr(output)
      }
      df[1, "Missing"] <- sum(is.na(x_var))
      if (show.tests) {
        df[1, "pTest"] <- "ChiSq"
        if (effSize) {
          df[1, "effStat"] <- "Cramer's V"
        }
      }
    }
    else if (any(chi.test.rm(cont_table)$expected < 5)) {
      fisher_test <- fisher.test.rm(cont_table)
      df[1, "p-value"] <-fisher_test$p.value
      if (effSize) {
        output <- calc_CramerV(fisher_test)
        df[1, "Effect Size (95% CI)"] <- psthr(output)
      }
      if (show.tests) {
        df[1, "pTest"] <- "Fisher Exact"
        if (effSize) {
          df[1, "effStat"] <- "Cramer's V"
        }
      }
    }
  }
  df[1, "Missing"] <- sum(is.na(x_var))
  return(df)
}
