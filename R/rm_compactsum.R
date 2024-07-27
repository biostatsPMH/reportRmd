#' Output a compact summary table
#'
#' Returns a data frame corresponding to a descriptive table.
#'
#' Comparisons for categorical variables default to chi-square tests, but if
#' there are counts of <5 then the Fisher Exact test will be used. For grouping
#' variables with two levels, either t-tests (mean) or wilcoxon tests (median)
#' will be used for numerical variables. Otherwise, ANOVA (mean) or kruskal-
#' wallis tests will be used. The statistical test used can be displayed by
#' specifying show.tests = TRUE. Statistical tests and effect sizes for grp and/
#' or xvars with less than 2 counts in any level will not be shown.
#'
#'
#'
#' @param data dataframe containing data
#' @param xvars character vector with the names of covariates to include in
#'   table
#' @param grp character with the name of the grouping variable
#' @param use_mean logical indicating whether mean and standard deviation will
#'   be returned for continuous variables instead of median. Otherwise, can
#'   specify for individual variables using a character vector containing the
#'   names of covariates to return mean and sd for (if use_mean is not supplied,
#'   all covariates will have median summaries). See examples
#' @param caption character containing table caption (default is no caption)
#' @param tableOnly logical, if TRUE then a dataframe is returned, otherwise a
#'   formatted printed object is returned (default is FALSE)
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
#' @param nicenames logical indicating if you want to replace . and _ in strings
#'.  with a space
#' @param iqr logical indicating if you want to display the interquartile range
#'   (Q1, Q3) as opposed to (min, max) in the summary for continuous variables
#' @param all.stats logical indicating if all summary statistics (Q1, Q3 + min,
#'   max on a separate line) should be displayed. Overrides iqr
#' @param pvalue logical indicating if you want p-values included in the table
#' @param effSize logical indicating if you want effect sizes and their 95%
#'   confidence intervals included in the table. Effect sizes calculated include
#'   Cramer's V for categorical variables, and Cohen's d, Wilcoxon r,
#'   Epsilon-squared, or Omega-squared for numeric/continuous variables
#' @param p.adjust p-adjustments to be performed
#' @param unformattedp logical indicating if you would like the p-value to be
#'   returned unformatted (ie. not rounded or prefixed with '<'). Best used with
#'   tableOnly = T and outTable function. See examples
#' @param show.tests logical indicating if the type of statistical test and
#'   effect size (if effSize = TRUE) used should be shown in a column beside the
#'   p-values.
#' @param full logical indicating if you want the full sample included in the
#'   table, ignored if grp is not specified
#' @param percentage choice of how percentages are presented, either column
#'   (default) or row
#' @returns A character vector of the table source code, unless tableOnly = TRUE
#'   in which case a data frame is returned. The output has the following
#'   attribute:
#'
#'   * "description", which describes what is included in the
#'   output table and the type of statistical summary for each covariate. When
#'   applicable, the types of statistical tests used will be included. If
#'   effSize = TRUE, the effect sizes for each covariate will also be mentioned.
#'
#' @returns The "Missing" column of the output table refers to the total missing
#'   values for each xvar across all groups, rather than missingness within the
#'   grp variable. If the number of missing values for all xvars specified is 0,
#'   the "Missing" column will be removed from the table.
#' @references Smithson, M. (2002). Noncentral Confidence Intervals for
#'   Standardized Effect Sizes. (07/140 ed., Vol.
#'   140). SAGE Publications. \url{https://doi.org/10.4135/9781412983761.n4}
#' @references Steiger, J. H. (2004). Beyond the F Test: Effect Size Confidence
#'   Intervals and Tests of Close Fit in the Analysis of Variance and Contrast
#'   Analysis. Psychological Methods, 9(2), 164–182.
#'   \url{https://doi.org/10.1037/1082-989X.9.2.164}
#' @references Kelley, T. L. (1935). An Unbiased Correlation Ratio Measure.
#'   Proceedings of the National Academy of Sciences - PNAS, 21(9), 554–559.
#'   \url{https://doi.org/10.1073/pnas.21.9.554}
#' @references Okada, K. (2013). Is Omega Squared Less Biased? A Comparison of
#'   Three Major Effect Size Indices in One-Way ANOVA. Behavior Research
#'   Methods, 40(2), 129-147.
#' @references Breslow, N. (1970). A generalized Kruskal-Wallis test for
#' comparing K samples subject to unequal patterns of censorship. Biometrika,
#' 57(3), 579-594.
#' @references FRITZ, C. O., MORRIS, P. E., & RICHLER, J. J. (2012). Effect Size
#' Estimates: Current Use, Calculations, and Interpretation. Journal of
#' Experimental Psychology. General, 141(1), 2–18.
#' \url{https://doi.org/10.1037/a0024338}
#'
#' @examples
#' data("pembrolizumab")
#' rm_compactsum(data = pembrolizumab, xvars = c("age",
#' "change_ctdna_group", "l_size", "pdl1"), grp = "sex", use_mean = "age",
#' digits = c("age" = 2, "l_size" = 3), digits.cat = 1, iqr = TRUE,
#' show.tests = TRUE)
#'
#' # To show effect sizes
#' rm_compactsum(data = pembrolizumab, xvars = c("age",
#' "change_ctdna_group"), grp = "sex", use_mean = "age", digits = 2,
#' effSize = TRUE, show.tests = TRUE)
#'
#' # To return unformatted p-values
#' rm_compactsum(data = pembrolizumab, xvars = c("l_size",
#' "change_ctdna_group"), grp = "cohort", effSize = TRUE, unformattedp = TRUE)
#'
#' # To view self-generated description
#' summary_tab <- rm_compactsum(data=pembrolizumab, xvars =
#' c("change_ctdna_group", "orr", "age"), grp = "cohort", effSize = TRUE,
#' show.tests = TRUE)
#' attr(summary_tab, "description")
#'
#' @export
rm_compactsum <- function(data, xvars, grp, use_mean, caption = NULL, tableOnly = FALSE, covTitle = "", digits = 1, digits.cat = 0,  nicenames = TRUE, iqr = FALSE, all.stats = FALSE, pvalue = TRUE, effSize = FALSE, p.adjust = "none", unformattedp = FALSE, show.tests = FALSE, full = TRUE, percentage = "col") {
  if (missing(data))
    stop("data is a required argument")
  if (missing(xvars))
    stop("xvars is a required argument")
  if (!inherits(data, "data.frame"))
    stop("data must be supplied as a data frame.")
  if (!inherits(xvars, "character"))
    stop("xvars must be supplied as a character vector or string indicating variables in data")
  missing_vars = setdiff(xvars, names(data))
  if (length(missing_vars) > 0) {
    stop(paste0("These xvars are not in the data: '", paste0(missing_vars, collapse = "', '"), "'"))
  }
  if (missing(use_mean)) {
    use_mean <- FALSE
  }
  if (!missing(grp)) {
    if (!inherits(grp, "character") | length(grp) > 1)
      stop("grp must be supplied as a string indicating a variable in data")
  }
  if (!((is.logical(use_mean) & length(use_mean) == 1) | (is.character(use_mean) & length(use_mean) > 0))) {
    stop("use_mean must be a character vector or a logical")
  }
  if (!is.numeric(digits)) {
    stop("digits must be a single numeric or a vector of digits")
  }
  if (!(is.numeric(digits.cat) & length(digits.cat) == 1)) {
    stop("digits.cat must be a single numeric")
  }
  if (!(percentage %in% c("col", "row"))) {
    stop("percentage argument must be either 'row' or 'col'")
  }
  argList <- as.list(match.call(expand.dots = TRUE)[-1])
  argsToPass <- intersect(names(formals(xvar_function)), names(argList))
  argsToPass <- setdiff(argsToPass,"xvars")
  args <- argList[argsToPass]

  dt <- as.name(args$data)

  xvars <- unique(xvars)

  for (xvar in xvars) {
    if (is.logical(data[[xvar]]) || is.character(data[[xvar]])) {
      data[[xvar]] <- as.factor(data[[xvar]])
      args$data <- data
    }
  }
  if (!missing(grp)) {
    if (!(grp %in% names(data))) {
      stop("grp is not in the data")
    }
    if (grp %in% xvars){
      warning(paste(grp,'is the grouping variable and cannot appear as a covariate. \n',
                    'It is omitted from the output.'))
      xvars <- setdiff(xvars, grp)
    }
    if (is.logical(data[[grp]]) | is.character(data[[grp]]) | (is.numeric(data[[grp]]) & length(unique(data[[grp]])) <= 5)) {
      data[[grp]] <- as.factor(data[[grp]])
      args$data <- data
    }
    else if (is.numeric(data[[grp]]) & length(unique(data[[grp]])) > 5) {
      stop("Convert grp to a factor")
    }
  }
  if (!missing(grp)) {
    grp_missing <- length(which(is.na(data[[grp]])))
    if (grp_missing > 0) {
      message(paste0("There are ", grp_missing, " missing cases for grouping variable '", grp, "'."))
    }
  }
  if (!missing(grp) & (effSize | show.tests | pvalue)) {
    grp_tab <- table(data[[grp]])
    if (any(grp_tab < 2)) {
      warning("Small counts in '", grp, "'. No statistical tests or effect size will be reported.")
      args$effSize = FALSE
      args$pvalue = FALSE
      args$show.tests = FALSE
    }
    else if (any(grp_tab < 5)) {
      warning("Small sample size in '", grp, "' group may lead to unstable effect sizes.")
    }
  }
  for (xvar in xvars) {
    if (inherits(data[[xvar]], "Date") || inherits(data[[xvar]], "POSIXt")) {
      xvars <- setdiff(xvars, xvar)
      warning(paste("date variable", xvar, "will be ignored"))
    }
  }
  for (xvar in xvars) {
    if (is.character(data[[xvar]]) | is.logical(data[[xvar]])) {
      data[[xvar]] <- as.factor(data[[xvar]])
      args$data <- data
    }
  }
  ignored_xvars <- c()
  if (!(missing(use_mean))) {
    if (!is.logical(use_mean)) {
      for (xvar in use_mean) {
        if (!(xvar %in% names(data))) {
          stop(paste0("variable '", xvar, "' in use_mean is not in data '", dt, "'"))
        }
        if (!(xvar %in% xvars)) {
          stop(paste0("variable '", xvar, "' in use_mean is not in xvars"))
        }
        if (is.factor(data[[xvar]]) | is_binary(data[[xvar]]) | grepl(class(data[[xvar]]),"factor")) {
          ignored_xvars <- c(ignored_xvars, xvar)
        }
      }
      if (length(ignored_xvars) > 0) {
        warning(paste("use_mean will be ignored for non-numerical xvars:", paste0(ignored_xvars, collapse = ", ")))
      }
    }
  }
  ignored_xvars <- c()
  if (!(missing(digits))) {
    if (length(digits) > 1) {
      for (xvar in names(digits)) {
        if (!(xvar %in% names(data))) {
          stop(paste0("variable '", xvar, "' in digits is not in data '", dt, "'"))
        }
        if (!(xvar %in% xvars)) {
          stop(paste0("variable '", xvar, "' in digits is not in xvars"))
        }
        if (is.factor(data[[xvar]]) | is_binary(data[[xvar]]) | grepl(class(data[[xvar]]),"factor")) {
          ignored_xvars <- c(ignored_xvars, xvar)
        }
      }
      if (length(ignored_xvars) > 0) {
        warning(paste("digits will be ignored for non-numerical xvars:", paste0(ignored_xvars, collapse = ", ")))
      }
    }
  }
  if (tableOnly) {
    args$covTitle <- "Covariate"
  }

  output_list <- NULL
  for (xvar in xvars) {
    if (grepl(class(data[[xvar]]),"factor") & length(unique(na.omit(data[[xvar]]))) == 2) {
      class(xvar) <- c(class(xvar),"rm_two_level")
    }
    else if (inherits(data[[xvar]],"factor")) {
      class(xvar) <- c(class(xvar),"rm_categorical")
    }
    else if (is.numeric(data[[xvar]]) && is_binary(data[[xvar]])) {
      class(xvar) <- c(class(xvar), "rm_binary")
    }
    else if (is.numeric(data[[xvar]])) {
      if (identical(use_mean, FALSE)) {
        class(xvar) <- c(class(xvar),"rm_median")
      }
      else {
        if (grepl(class(use_mean), "character")) {
          if (xvar %in% use_mean) {
            class(xvar) <- c(class(xvar),"rm_mean")
          }
          else {
            class(xvar) <- c(class(xvar),"rm_median")
          }
        }
        else if (grepl(class(use_mean), "logical")) {
          if (identical(use_mean, TRUE)) {
            class(xvar) <- c(class(xvar),"rm_mean")
          }
          else {
            class(xvar) <- c(class(xvar),"rm_median")
          }
        }
        else {
          stop("use_mean must be a logical or character vector")
        }
      }
    }
    else {
      stop("xvar must be binary, numeric, or categorical")
    }
    if (length(names(digits)) >= 1) {
      if (xvar %in% names(digits)) {
        args$digits = digits[[xvar]]
      }
      else {
        args$digits = 1
      }
    }
    else if (length(digits) == 1 && is.numeric(digits)) {
      args$digits = digits
    }
    args$xvar = xvar
    output_list[[xvar]] <- do.call(xvar_function, args)
  }
  result <-dplyr::bind_rows(output_list)
  if (all(na.omit(result[["Missing"]]) == 0)) {
    result <- result[, -which(names(result) == "Missing")]
  }
  if (!full) {
    result <- result[, -grep("^Full Sample", names(result))]
  }
  if ("p-value" %in% colnames(result)) {
    if (!pvalue) {
      result <- result[, -which(names(result) == "p-value")]
    }
    else {
      method <- p.adjust
      result[["p-value"]] <- p.adjust(result[["p-value"]], method = method)
      if (!unformattedp) {
        result[["p-value"]] <- formatp(result[["p-value"]])
      }
    }
  }
  lbl <- c()
  for (xvar in xvars) {
    if (inherits(data[[xvar]],"factor") & length(unique(na.omit(data[[xvar]]))) > 2) {
      lbl <- c(lbl, xvar, levels(data[[xvar]]))
    }
    else {
      lbl <- c(lbl, xvar)
    }
  }
  if (nicenames) {
    if (typeof(args$data) == "symbol") {
    result[, 1] <- replaceLbl(args$data, lbl)
    }
    else {
      result[, 1] <- replaceLbl(dt, lbl)
    }
  }
  result[, 1] <- paste0(result[, 1],result$`disp`)
  result$`disp` <- NULL
  attr(result, "description") <- generate_description(xvars, output_list)
  if (tableOnly) {
    return(result)
  }
  to_indent <- c()
  n <- 0
  for (xvar in xvars) {
    n = n+1
    if (inherits(data[[xvar]],"factor") & length(unique(na.omit(data[[xvar]]))) > 2) {
      to_indent <- c(to_indent, n + (1:length(unique(na.omit(data[[xvar]])))))
      n = n + length(unique(na.omit(data[[xvar]])))
    }
  }
  nicetable <- outTable(result, caption = caption, nicenames = nicenames, to_indent = to_indent)
  attr(nicetable, "description") <- generate_description(xvars, output_list)
  return(nicetable)
}

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
#'   confidence intervals included in the table. Effect sizes calculated include
#'   Cramer's V for categorical variables, and Cohen's d, Wilcoxon r,
#'   Epsilon-squared, or Omega-squared for numeric/continuous variables
#' @param show.tests logical indicating if the type of statistical test and
#'   effect size (if effSize = TRUE) used should be shown in a column beside the
#'   p-values.
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
  if (!(pvalue | effSize)) {
    show.tests = FALSE
  }
  ## *** why change the class?
  class(xvar) <- "character"
  df <- data.frame(Covariate = xvar)
  df[["disp"]] <- " n (%)"
  if (covTitle == "") {
    names(df$`Covariate`) <- " "
  }
  else {
    names(df) <- covTitle
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
    no_na <- subset(data, !is.na(data[[xvar]]))
    no_na_tab <- table(no_na[[grp]])
    if (any(no_na_tab < 2)) {
      effSize <- FALSE
      pvalue <- FALSE
      show.tests <- FALSE
    }
    if (pvalue | effSize | show.tests) {
      if (!any(chi.test.rm(cont_table)$expected < 5)) {
        chisq_test <- chi.test.rm(cont_table)
        df[1, "p-value"] <- chisq_test$p.value
        if (effSize) {
          output <- calc_CramerV(chisq_test)
          df[1, "Effect Size (95% CI)"] <- format_delta(output)
          attr(df, "eff_size") <- "Cramer's V"
        }
        df[1, "Missing"] <- sum(is.na(x_var))
        if (show.tests & pvalue) {
          df[1, "pTest"] <- "ChiSq"
        }
        if (pvalue) {
          attr(df, "stat_test") <- "chi-square test"
        }
        if (show.tests & effSize) {
          df[1, "effStat"] <- "Cramer's V"
        }
      }
      else if (any(chi.test.rm(cont_table)$expected < 5)) {
        fisher_test <- fisher.test.rm(cont_table)
        df[1, "p-value"] <-fisher_test$p.value
        if (effSize) {
          output <- calc_CramerV(fisher_test)
          df[1, "Effect Size (95% CI)"] <- format_delta(output)
          attr(df, "eff_size") <- "Cramer's V"
        }
        if (show.tests & pvalue) {
          df[1, "pTest"] <- "Fisher Exact"
        }
        if (pvalue) {
          attr(df, "stat_test") <- "Fisher's Exact Test"
        }
        if (show.tests & effSize) {
          df[1, "effStat"] <- "Cramer's V"
        }
      }
    }
  }
  df[1, "Missing"] <- sum(is.na(x_var))
  attr(df, "stat_sum") <- "counts (%)"
  return(df)
}

xvar_function.rm_mean <- function(xvar, data, grp, covTitle = "", digits = 1, digits.cat = 0, iqr = FALSE, all.stats = FALSE, pvalue = TRUE, effSize = FALSE, show.tests = FALSE, percentage = "col") {
  if (!(pvalue | effSize)) {
    show.tests = FALSE
  }
  class(xvar) <- "character"
  df <- data.frame(Covariate = xvar)
  df[["disp"]] <- " Mean (sd)"
  if (covTitle == "") {
    names(df$`Covariate`) <- " "
  }
  else {
    names(df$`Covariate`) <- covTitle
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
    no_na <- subset(data, !is.na(data[[xvar]]))
    no_na_tab <- table(no_na[[grp]])
    if (any(no_na_tab < 2)) {
      effSize <- FALSE
      pvalue <- FALSE
      show.tests <- FALSE
    }
    if (pvalue | effSize | show.tests) {
      if (length(levels(group_var)) == 2) {
        t_test <- t.test.rm(x_var, group_var)
        df[, "p-value"] <- t_test$p.value
        N <- nrow(data)
        if (effSize) {
          output <- calc_cohenD(t_test)
          df[, "Effect Size (95% CI)"] <- format_delta(output)
          attr(df, "eff_size") <- "Cohen's d"
        }
        df[, "Missing"] <- sum(is.na(x_var))
        if (show.tests & pvalue) {
          df[, "pTest"] <- "t-test"
        }
        if (pvalue) {
          attr(df, "stat_test") <- "independent t-test"
        }
        if (show.tests & effSize) {
          df[, "effStat"] <- "Cohen's d"
        }
      }
      else if (length(levels(group_var)) > 2) {
        anova_test <- stats::aov(x_var ~ group_var)
        df[, "p-value"] <- summary(anova_test)[[1]][["Pr(>F)"]][1]
        if (effSize) {
          output <- calc_omegaSq(anova_test)
          df[, "Effect Size (95% CI)"] <- format_delta(output)
          attr(df, "eff_size") <- "Omega Squared"
        }
        df[, "Missing"] <- sum(is.na(x_var))
        if (show.tests & pvalue) {
          df[, "pTest"] <- "ANOVA"
        }
        if (pvalue) {
          attr(df, "stat_test") <- "ANOVA"
        }
        if (show.tests & effSize) {
          df[, "effStat"] <- "Omega Sq"
        }
      }
    }
  }
  df[, "Missing"] <- sum(is.na(x_var))
  attr(df, "stat_sum") <- "mean (sd)"
  return(df)
}

xvar_function.rm_median <- function(xvar, data, grp, covTitle = "", digits = 1, digits.cat = 0, iqr = FALSE, all.stats = FALSE, pvalue = TRUE, effSize = FALSE, show.tests = FALSE, percentage = "col") {
  if (!(pvalue | effSize)) {
    show.tests = FALSE
  }
  class(xvar) <- "character"
  x_var <- data[[xvar]]
  if (all.stats) {
    df <- data.frame(Covariate = c(xvar, "  Median (Q1, Q3)", "  Range (min, max)"))
    if (covTitle == "") {
      names(df$`Covariate`) <- " "
    }
    else {
      names(df$`Covariate`) <- covTitle
    }
    bracket_iqr <- paste0("(", format(round(stats::quantile(x_var, na.rm = TRUE, prob = 0.25), digits), nsmall = digits), ", ", format(round(stats::quantile(x_var, na.rm = TRUE, prob = 0.75), digits), nsmall = digits), ")")
    bracket_range <- paste0("(", format(round(min(x_var, na.rm = TRUE), digits), nsmall = digits), ", ", format(round(max(x_var, na.rm = TRUE), digits), nsmall = digits), ")")
    df[2, paste0("Full Sample (n=", nrow(data), ")")] <- paste0(format(round(median(x_var, na.rm = TRUE), digits), nsmall = digits), " ", bracket_iqr)
    df[3, paste0("Full Sample (n=", nrow(data), ")")] <- bracket_range
  }
  else {
    df <- data.frame(Covariate = xvar)
    df[["disp"]] <- ifelse(!iqr, " Median (Min, Max)", " Median (Q1, Q3)")
    if (covTitle == "") {
      names(df$`Covariate`) <- " "
    }
    else {
      names(df$`Covariate`) <- covTitle
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
    no_na <- subset(data, !is.na(data[[xvar]]))
    no_na_tab <- table(no_na[[grp]])
    if (any(no_na_tab < 2)) {
      effSize <- FALSE
      pvalue <- FALSE
      show.tests <- FALSE
    }
    if (pvalue | effSize | show.tests) {
      if (length(unique(group_var)) == 2) {
        wilcox_test <- wilcox.test.rm(x_var, group_var)
        df[1, "p-value"] <- wilcox_test$p.value
        if (effSize) {
          output <- calc_WilcoxonR(wilcox_test)
          df[1, "Effect Size (95% CI)"] <- format_delta(output)
          attr(df, "eff_size") <- "Wilcoxon R"
        }
        df[1, "Missing"] <- sum(is.na(x_var))
        if (show.tests & pvalue) {
          df[1, "pTest"] <- "Wilcoxon Rank Sum"
        }
        if (pvalue) {
          attr(df, "stat_test") <- "Wilcoxon rank sum test"
        }
        if (show.tests & effSize) {
          df[1, "effStat"] <- "Wilcoxon r"
        }
      }
      else if (length(levels(group_var)) > 2) {
        kruskal_test <- kruskal.test.rm(x_var, group_var)
        df[1, "p-value"] <- kruskal_test$p.value
        if (effSize) {
          output <- calc_epsilonSq(kruskal_test)
          df[1, "Effect Size (95% CI)"] <- format_delta(output)
          attr(df, "eff_size") <- "Epsilon-squared"
        }
        df[1, "Missing"] <- sum(is.na(x_var))
        if (show.tests & pvalue) {
          df[1, "pTest"] <- "Kruskal Wallis"
        }
        if (pvalue) {
          attr(df, "stat_test") <- "Kruskal Wallis Test"
        }
        if (show.tests & effSize) {
          df[1, "effStat"] <- "Epsilon sq"
        }
      }
    }
  }
  df[1, "Missing"] <- sum(is.na(x_var))
  if (iqr) {
    attr(df, "stat_sum") <- "median (IQR)"
  }
  else {
    attr(df, "stat_sum") <- "median (min/max)"
  }
  return(df)
}

xvar_function.rm_categorical <- function(xvar, data, grp, covTitle = "", digits = 1, digits.cat = 0, iqr = FALSE, all.stats = FALSE, pvalue = TRUE, effSize = FALSE, show.tests = FALSE, percentage = "col") {
  if (!(pvalue | effSize)) {
    show.tests = FALSE
  }
  class(xvar) <- "character"
  x_var <- data[[xvar]]
  rows <- c(paste0(xvar))
  for (xvar_level in levels(x_var)) {
    rows <- append(rows, xvar_level)
  }
  df <- data.frame(Covariate = rows, "disp" = rep("", length(rows)))
  df[1, "disp"] <-  " n (%)"
  if (covTitle == "") {
    names(df$`Covariate`) <- " "
  }
  else {
    names(df$`Covariate`) <- covTitle
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
    no_na <- subset(data, !is.na(data[[xvar]]))
    no_na_tab <- table(no_na[[grp]])
    if (any(no_na_tab < 2)) {
      effSize <- FALSE
      pvalue <- FALSE
      show.tests <- FALSE
    }
    if (pvalue | effSize | show.tests) {
      if (!any(chi.test.rm(cont_table)$expected < 5)) {
        chisq_test <- chi.test.rm(cont_table)
        df[1, "p-value"] <- chisq_test$p.value
        if (effSize) {
          output <- calc_CramerV(chisq_test)
          df[1, "Effect Size (95% CI)"] <- format_delta(output)
          attr(df, "eff_size") <- "Cramer's V"
        }
        df[1, "Missing"] <- sum(is.na(x_var))
        if (show.tests & pvalue) {
          df[1, "pTest"] <- "ChiSq"
        }
        if (pvalue) {
          attr(df, "stat_test") <- "chi-square test"
        }
        if (show.tests & effSize) {
          df[1, "effStat"] <- "Cramer's V"
        }
      }

      else if (any(chi.test.rm(cont_table)$expected < 5)) {
        fisher_test <- fisher.test.rm(cont_table)
        df[1, "p-value"] <- fisher_test$p.value
        if (effSize) {
          output <- calc_CramerV(fisher_test)
          df[1, "Effect Size (95% CI)"] <- format_delta(output)
          attr(df, "eff_size") <- "Cramer's V"
        }
        if (show.tests & pvalue) {
          df[1, "pTest"] <- "Fisher Exact"
        }
        if (pvalue) {
          attr(df, "stat_test") <- "Fisher's Exact Test"
        }
        if (show.tests & effSize) {
          df[1, "effStat"] <- "Cramer's V"
        }
      }
    }
  }
  attr(df, "stat_sum") <- "counts (%)"
  df[1, "Missing"] <- sum(is.na(x_var))

  return(df)
}

xvar_function.rm_two_level <- function(xvar, data, grp, covTitle = "", digits = 1, digits.cat = 0, iqr = FALSE, all.stats = FALSE, pvalue = TRUE, effSize = FALSE, show.tests = FALSE, percentage = "col") {
  if (!(pvalue | effSize)) {
    show.tests = FALSE
  }
  class(xvar) <- "character"
  temp <- data.frame()
  x_var <- data[[xvar]]
  unique_levels <- unique(x_var)
  unique_levels <- sort(unique_levels)
  show_level <- as.character(unique_levels[2])
  binary_column <- ifelse(x_var == unique_levels[1], 0, 1)
  if (!missing(grp)) {
    temp <- subset(data, select = grp)
    temp[[xvar]] <- binary_column
  }
  else {
    temp <- subset(data, select = xvar)
    temp[[xvar]] <- binary_column
  }
  df <- data.frame(Covariate = xvar)
  df[["disp"]] <-  paste0(" - ", show_level, " (n (%))")
  if (covTitle == "") {
    names(df$`Covariate`) <- " "
  }
  else {
    names(df$`Covariate`) <- covTitle
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
    no_na <- subset(data, !is.na(data[[xvar]]))
    no_na_tab <- table(no_na[[grp]])
    if (any(no_na_tab < 2)) {
      effSize <- FALSE
      pvalue <- FALSE
      show.tests <- FALSE
    }
    if (pvalue | effSize | show.tests) {
      if (!any(chi.test.rm(cont_table)$expected < 5)) {
        chisq_test <- chi.test.rm(cont_table)
        df[1, "p-value"] <- chisq_test$p.value
        if (effSize) {
          output <- calc_CramerV(chisq_test)
          df[1, "Effect Size (95% CI)"] <- format_delta(output)
          attr(df, "eff_size") <- "Cramer's V"
        }
        df[1, "Missing"] <- sum(is.na(x_var))
        if (show.tests & pvalue) {
          df[1, "pTest"] <- "ChiSq"
        }
        if (pvalue) {
        attr(df, "stat_test") <- "chi-square test"
        }
        if (show.tests & effSize) {
          df[1, "effStat"] <- "Cramer's V"
        }
      }
      else if (any(chi.test.rm(cont_table)$expected < 5)) {
        fisher_test <- fisher.test.rm(cont_table)
        df[1, "p-value"] <-fisher_test$p.value
        if (effSize) {
          output <- calc_CramerV(fisher_test)
          df[1, "Effect Size (95% CI)"] <- format_delta(output)
          attr(df, "eff_size") <- "Cramer's V"
        }
        if (show.tests & pvalue) {
          df[1, "pTest"] <- "Fisher Exact"
        }
        if (pvalue) {
          attr(df, "stat_test") <- "Fisher's Exact Test"
        }
        if (show.tests & effSize) {
          df[1, "effStat"] <- "Cramer's V"
        }
      }
    }
  }
  df[1, "Missing"] <- sum(is.na(x_var))
  attr(df, "stat_sum") <- "counts (%)"
  return(df)
}

delta_CI <- function(htest,CIwidth=0.95){
  # determine what test result is being passed in
  if (inherits(htest,"aov")){
    class(htest) <- c(class(htest),"rm_F")
  } else if (inherits(htest,"htest")){
    if (grepl("Chi",htest$method)){
      htest <- uncorrectedChi(htest)
      class(htest) <- c(class(htest),"rm_chi")
    } else if (grepl("Kruskal",htest$method)) {
      class(htest) <- c(class(htest),"rm_chi")
    } else if (grepl("Wilcoxon rank sum",htest$method)) {
      class(htest) <- c(class(htest),"rm_chi")
    } else if (grepl("t-test",htest$method)) {
      class(htest) <- c(class(htest),"rm_t")
    } else if (grepl("Fisher",htest$method)) {
      class(htest) <- c(class(htest),"rm_chi")
    }}  else stop("htest must be the output of a call to t.test, chisq.test, kruskal.test, aov or fisher.test.rm")
  dL <- delta_l(htest,CIwidth)
  dU <- delta_u(htest,CIwidth)
  ci <- c(dL,dU)
  attr(ci,"confidence") <- CIwidth
  return(ci)
}

delta_l <- function(htest,CIwidth){
  UseMethod("delta_l")
}

delta_l.default <- function(...){
  stop("No default method for delta_l. The test type must be known")
}

delta_l.rm_F <- function(htest,CIwidth) {
  hsum <- summary(htest)[[1]]
  Fstat <- hsum[["F value"]][1]
  nu1 <- hsum[["Df"]][1]
  nu2 <- hsum[["Df"]][2]
  ulim <- 1 - (1-CIwidth)/2
  if (stats::pf(Fstat,nu1,nu2) < ulim) return(0)

  lc <- c(.001,Fstat/2,Fstat)
  while(stats::pf(Fstat,nu1,nu2,lc[1])<ulim) {
    lc <- c(lc[1]/4,lc[1],lc[3])
  }
  while(stats::pf(Fstat,nu1,nu2,lc[3])>ulim) {
    lc <- c(lc[1],lc[3],lc[3]+Fstat)
  }
  diff <- 1
  while(diff > .00001) {
    if(stats::pf(Fstat,nu1,nu2,lc[2])<ulim){
      lc <- c(lc[1],(lc[1]+lc[2])/2,lc[2])
    }  else lc <- c(lc[2],(lc[2]+lc[3])/2,lc[3])
    diff <- abs(stats::pf(Fstat,nu1,nu2,lc[2]) - ulim)
  }
  names(lc) <- NULL
  return(lc[2])
}

delta_l.rm_chi <- function(htest,CIwidth) {
  chival <- htest$statistic; df <- htest$parameter
  ulim <- 1 - (1-CIwidth)/2
  if (stats::pchisq(chival,df)< ulim) return(0)

  lc <- c(.001,chival/2,chival)
  while(stats::pchisq(chival,df,lc[1])<ulim) {
    if(stats::pchisq(chival,df)<ulim)
      return(0)
    lc <- c(lc[1]/4,lc[1],lc[3])
  }
  diff <- 1
  while(diff > .00001) {
    if(stats::pchisq(chival,df,lc[2])<ulim){
      lc <- c(lc[1],(lc[1]+lc[2])/2,lc[2])
    }  else lc <- c(lc[2],(lc[2]+lc[3])/2,lc[3])
    diff <- abs(stats::pchisq(chival,df,lc[2]) - ulim)
  }
  names(lc) <- NULL
  return(lc[2])
}

delta_l.rm_t <- function(htest,CIwidth) {
  tval <- abs(htest$statistic); df <- htest$parameter
  ulim <- 1 - (1-CIwidth)/2
  if (stats::pt(tval,df) < ulim) return(0)

  lc <- c(-tval,tval/2,tval)
  while(stats::pt(tval,df,lc[1])<ulim) {
    lc <- c(lc[1]-tval,lc[1],lc[3])
  }
  diff <- 1
  while(diff > .00001) {
    if(stats::pt(tval,df,lc[2])<ulim)
      lc <- c(lc[1],(lc[1]+lc[2])/2,lc[2])
    else lc <- c(lc[2],(lc[2]+lc[3])/2,lc[3])
    diff <- abs(stats::pt(tval,df,lc[2]) - ulim)
  }
  names(lc) <- NULL
  return(lc[2])
}

delta_u <- function(htest,CIwidth){
  UseMethod("delta_u")
}

delta_u.default <- function(...){
  stop("No default method for delta_u. The test type must be known")
}

delta_u.rm_F <- function(htest,CIwidth) {
  hsum <- summary(htest)[[1]]
  Fstat <- hsum[["F value"]][1]
  nu1 <- hsum[["Df"]][1]
  nu2 <- hsum[["Df"]][2]
  alpha <- 1-CIwidth
  if (stats::pf(Fstat,nu1,nu2)<alpha/2) return(0)
  uc <- c(Fstat,2*Fstat,3*Fstat)
  llim <- (1-CIwidth)/2
  while(stats::pf(Fstat,nu1,nu2,uc[1])<llim) {
    uc <- c(uc[1]/4,uc[1],uc[3])
  }
  while(stats::pf(Fstat,nu1,nu2,uc[3])>llim) {
    uc <- c(uc[1],uc[3],uc[3]+Fstat)
  }
  diff <- 1
  while(diff > .00001) {
    if(stats::pf(Fstat,nu1,nu2,uc[2])<llim)
      uc <- c(uc[1],(uc[1]+uc[2])/2,uc[2])
    else uc <- c(uc[2],(uc[2]+uc[3])/2,uc[3])
    diff <- abs(stats::pf(Fstat,nu1,nu2,uc[2]) - llim)
    lcdf <- stats::pf(Fstat,nu1,nu2,uc[2])
  }
  names(uc) <- NULL
  return(uc[2])
}

delta_u.rm_chi <- function(htest,CIwidth) {
  alpha <- 1-CIwidth
  chival <- htest$statistic; df <- htest$parameter
  if (stats::pchisq(chival,df)<alpha/2) return(0)
  uc <- c(chival,2*chival,3*chival)
  llim <- (1-CIwidth)/2
  while(stats::pchisq(chival,df,uc[1])<llim) {
    uc <- c(uc[1]/4,uc[1],uc[3])
  }
  while(stats::pchisq(chival,df,uc[3])>llim) {
    uc <- c(uc[1],uc[3],uc[3]+chival)
  }
  diff <- 1
  while(diff > .00001) {
    if(stats::pchisq(chival,df,uc[2])<llim)
      uc <- c(uc[1],(uc[1]+uc[2])/2,uc[2])
    else uc <- c(uc[2],(uc[2]+uc[3])/2,uc[3])
    diff <- abs(pchisq(chival,df,uc[2]) - llim)
    lcdf <- stats::pchisq(chival,df,uc[2])
  }
  names(uc) <- NULL
  return(uc[2])
}

delta_u.rm_t <- function(htest,CIwidth) {
  alpha <- 1-CIwidth
  tval <- abs(htest$statistic); df <- htest$parameter
  if (stats::pt(tval,df)<alpha/2) return(0)
  llim <- (1-CIwidth)/2
  uc <- c(tval,1.5*tval,2*tval)
  while(stats::pt(tval,df,uc[3])>llim) {
    uc <- c(uc[1],uc[3],uc[3]+tval)
  }
  diff <- 1
  while(diff > .00001) {
    if(stats::pt(tval,df,uc[2])<llim)
      uc <- c(uc[1],(uc[1]+uc[2])/2,uc[2])
    else uc <- c(uc[2],(uc[2]+uc[3])/2,uc[3])
    diff <- abs(stats::pt(tval,df,uc[2]) - llim)
    lcdf <- stats::pt(tval,df,uc[2])
  }
  names(uc) <- NULL
  return(uc[2])
}

uncorrectedChi <- function(x) {
  stopifnot(inherits(x, "htest") )
  stopifnot("X-squared" %in% names(x[["statistic"]]))
  return(stats::chisq.test(x$observed,correct = FALSE))
}

fisher.test.rm <- function(x,...){
  chi.out <- suppressWarnings(stats::chisq.test(x,...))
  rtn <- stats::fisher.test(x,...)
  rtn$observed <- chi.out$observed
  rtn$statistic <- chi.out$statistic
  rtn$parameter <- chi.out$parameter
  rtn$k <- min(dim(chi.out$observed), na.rm = TRUE)
  rtn$N <- sum(chi.out$observed)
  rtn$method <- paste(rtn$method,"with additional chisq.test arguments")
  class(rtn) <- c(class(rtn),"fisher.test.rm")
  return(rtn)
}

chi.test.rm <- function(x,...){
  chi_test <- suppressWarnings(stats::chisq.test(x,...))
  chi_test$k <- min(dim(chi_test$observed), na.rm = TRUE)
  chi_test$N <- sum(chi_test$observed)
  class(chi_test) <- c(class(chi_test),"chi.test.rm")
  return(chi_test)
}

t.test.rm <- function(xvar,grp){
  t_test <- stats::t.test(xvar~grp)
  n <- unlist(lapply(levels(grp),function(g){
    length(na.omit(xvar[grp==g]))
  }))
  names(n) <- levels(grp)
  t_test$n <- n
  t_test$xvar <- xvar
  t_test$grp <- grp
  class(t_test) <- c(class(t_test),"t.test.rm")
  return(t_test)
}

wilcox.test.rm <- function(xvar,grp){
  xg <- na.omit(cbind(xvar,grp))
  x <- xg[xg[,2]==1,1]
  y <- xg[xg[,2]==2,1]
  wilcox_test <- suppressWarnings(stats::wilcox.test(x,y))
  Ua <- wilcox_test$statistic
  Ub <- suppressWarnings(stats::wilcox.test(y,x)$statistic)
  U1 = min(Ua,Ub)
  U2 = max(Ua,Ub)
  n1 <- length(na.omit(x))
  n2 <- length(na.omit(y))
  wilcox_test$U=min(U1,U2)
  Z = (U1+.5-(U1+U2)/2)/sqrt((n1*n2*(n1+n2+1)/12))
  n <- unlist(lapply(levels(grp),function(g){
    length(na.omit(xvar[grp==g]))
  }))
  names(n) <- levels(grp)
  wilcox_test$n <- n
  wilcox_test$Z <- Z
  wilcox_test$statistic=Z^2
  wilcox_test$parameter=1
  wilcox_test$xvar <- xvar
  wilcox_test$grp <- grp
  class(wilcox_test) <- c(class(wilcox_test),"wilcox.test.rm")
  return(wilcox_test)
}

kruskal.test.rm <- function(xvar,grp){
  kruskal_test <- stats::kruskal.test(xvar~grp)
  n <- length(na.omit(xvar))
  kruskal_test$n <- n
  kruskal_test$xvar <- xvar
  kruskal_test$grp <- grp
  class(kruskal_test) <- c(class(kruskal_test),"kruskal.test.rm")
  return(kruskal_test)
}

anova_toOmegaSq <- function(anova_summary){
  tbl <- anova_summary[[1]]
  num <- tbl[1,2]-tbl[1,1]*tbl[2,3]
  den <- tbl[1,2]+tbl[2,2]+tbl[2,3]
  return(num/den)
}

chi_toCramer <- function(chisq_test){
  obs <- chisq_test$observed
  V = sqrt(chisq_test$statistic/(sum(obs)*(min(dim(obs))-1)))
  return(V)
}

t_toCohen <-function(t_test){
  if (!inherits(t_test,"t.test.rm")) stop("t_test must be a function returned from t.test.rm")
  n1 <- t_test$n[1]
  n2 <- t_test$n[2]
  cohen=abs(t_test$statistic*sqrt((n1+n2)/(n1*n2)))
  return(cohen)
}

chi_toEpsilonSq <- function(kruskal_test){
  n <- kruskal_test$n
  epsilonSq <- kruskal_test$statistic/(n-1)
}

wilcox_effSize <- function(wilcox_test){
  abs(wilcox_test$Z/sqrt(sum(wilcox_test$n)))
}

lambda_toOmegaSq <- function(lambda,N){
  lambda/(lambda+N)
}

lambda_toCramer <- function(lambda,N,k,df){
  sqrt((lambda+df)/(N*(k-1)))
}

lambda_toEpsilon <- function(lamba,N,df){
  (lamba+df)/((N^2-1)/(N+1))
}

lambda_toCohen <- function(lambda,N1,N2){
  lambda/sqrt((N1*N2)/(N1+N2))
}

lambda_toWilcoxR <- function(lambda,N){
  sqrt(lambda/N)
}

calc_CramerV <- function(chisq_test, CIwidth = 0.95) {
  cramer <- chi_toCramer(chisq_test)
  df_cont <- data.frame(chisq_test$observed)
  df <- dplyr::bind_rows(lapply(1:nrow(df_cont),function(x){
    data.frame(df_cont[rep(x,df_cont$Freq[x]),-ncol(df_cont)])
  }))

  bs_cramer <- function(data,indices){
    dt <- data[indices,]
    new_chi <- chi.test.rm(dt[[1]],dt[[2]])

    chi_toCramer(new_chi)
  }
  b_cramer <- boot::boot(df,bs_cramer,R=1000)
  b_ci <- boot::boot.ci(b_cramer,conf = CIwidth,type="perc")
  output = c("cramer v"=cramer,lower=b_ci$perc[4],upper=b_ci$perc[5])
  return(output)
}

calc_cohenD <- function(t_test, CIwidth = 0.95) {
  cohen <- t_toCohen(t_test)
  df <- data.frame(xvar = t_test$xvar, grp = t_test$grp)

  bs_cohen <- function(data,indices){
    dt <- data[indices,]
    new_t <- tryCatch(t.test.rm(dt$xvar,dt$grp))
    t_toCohen(new_t)
  }

  b_cohen <- boot::boot(df,bs_cohen,R=1000)
  b_ci <- boot::boot.ci(b_cohen,conf = CIwidth,type="perc")
  output = c("cohen d"=cohen,lower=b_ci$perc[4],upper=b_ci$perc[5])
  return(output)
}

calc_WilcoxonR <- function(wilcox_test, CIwidth = 0.95) {
  r <- wilcox_effSize(wilcox_test)
  df <- data.frame(xvar = wilcox_test$xvar, grp = wilcox_test$grp)

  bs_r <- function(data,indices){
    dt <- data[indices,]
    new_wilcox <- wilcox.test.rm(dt$xvar, dt$grp)
    wilcox_effSize(new_wilcox)
  }
  b_r <- boot::boot(df,bs_r,R=1000)
  b_ci <- boot::boot.ci(b_r,conf = CIwidth,type="perc")
  output = c("wilcoxon r"=r,lower=b_ci$perc[4],upper=b_ci$perc[5])
  return(output)
}

calc_epsilonSq <- function(kruskal_test, CIwidth = 0.95) {
  eps <- chi_toEpsilonSq(kruskal_test)
  df <- data.frame(xvar = kruskal_test$xvar, grp = kruskal_test$grp)

  bs_epsilon <- function(data,indices){
    dt <- data[indices,]
    new_kruskal <- kruskal.test.rm(dt$xvar,dt$grp)

    chi_toEpsilonSq(new_kruskal)
  }
  b_epsilon <- boot::boot(df,bs_epsilon,R=1000)
  b_ci <- boot::boot.ci(b_epsilon,conf = CIwidth,type="perc")
  output = c("epsilon sq"=eps,lower=b_ci$perc[4],upper=b_ci$perc[5])
  return(output)
}

calc_omegaSq <- function(anova_test, CIwidth = 0.95){
  summary_anova <- summary(anova_test)
  omega <- anova_toOmegaSq(summary_anova)

  bs_omega <- function(data,indices){
    dt <- data[indices,]
    new_anova <- stats::update(anova_test,data=dt)
    new_summary <- summary(new_anova)
    anova_toOmegaSq(new_summary)
  }
  b_omega <- boot::boot(anova_test$model,bs_omega,R=1000)
  b_ci <- boot::boot.ci(b_omega,conf = CIwidth,type="perc")
  output = c("omega squared"=omega,lower=b_ci$perc[4],upper=b_ci$perc[5])
  return(output)
}

pstprn <- reportRmd:::pstprn

format_delta <- function(x,digits=2){
  out <- sapply(x, function(x){

    format(round(x,digits),nsmall=digits)

  })

  out <- unlist(out,use.names = FALSE)
  ci <- paste0(out[1], " (", out[2], ", ", out[3], ")")
  return(ci)
}

