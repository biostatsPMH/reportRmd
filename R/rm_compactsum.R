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
#'   p-values. Ignored if pvalue = FALSE
#' @param full logical indicating if you want the full sample included in the
#'   table, ignored if grp is not specified
#' @param percentage choice of how percentages are presented, either column
#'   (default) or row
#' @returns A character vector of the table source code, unless tableOnly = TRUE
#'   in which case a data frame is returned. The output has the following
#'   attribute:
#'
#'   * "paragraph", which describes what is included in the
#'   output table and the type of statistical summary for each covariate. When
#'   applicable, the types of statistical tests used will be included. If
#'   effSize = TRUE, the effect sizes for each covariate will also be mentioned.
#'
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
#' @export
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
  # args$covTitle = covTitle
  # args$digits = digits
  # args$digits.cat = digits.cat
  # args$iqr = iqr
  # args$all.stats = all.stats
  # args$pvalue = pvalue
  # args$effSize = effSize
  # args$show.tests = show.tests
  # args$percentage = percentage

  dt <- as.name(args$data)
  if (!missing(grp)) {
    if (!(grp %in% names(data))) {
      stop("grp is not in the data")
    }
    if (grp %in% xvars){
      warning(paste(grp,'is the grouping variable and can not appear as a covariate. \n',
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
    #### !!!! in xvar function, if effSIze then fdo table on grp and get num of ppl in each group. if any are fewer than 5 print warning of unstable effsize
  }
  ### NE Instead NA
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

  if (!pvalue) {
    show.tests <- FALSE
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
  # print(output_list)
  # descr_sum <- list()
  # for (xvar in xvars) {
  #   x_sum <- attr(output_list[[xvar]], "stat_sum")
  #   if (!is.null(x_sum)) {
  #     if (!(x_sum %in% names(descr_sum))) {
  #       descr_sum[[x_sum]] <- xvar
  #     }
  #     else {
  #       descr_sum[[x_sum]] <- c(descr_sum[[x_sum]], xvar)
  #     }
  #   }
  # }
  # descr_eff <- list()
  # for (xvar in xvars) {
  #   x_eff <- attr(output_list[[xvar]], "eff_size")
  #   if (!is.null(x_eff)) {
  #     if (!(x_eff %in% names(descr_eff))) {
  #       descr_eff[[x_eff]] <- xvar
  #     }
  #     else {
  #       descr_eff[[x_eff]] <- c(descr_eff[[x_eff]], xvar)
  #     }
  #   }
  # }
  # for (x_eff in names(descr_eff)) {
  #   print(paste0("The '", paste(descr_eff[[x_eff]], collapse = "', '"), "' variables use ", x_eff, " for effect size."))
  # }
  # for (x_sum in names(descr_sum)) {
  #   print(paste(x_sum, "statistical summaries are displayed for '", paste(descr_sum[[x_sum]], collapse = "', '"), "' variables."))
  # }
  # # descr <-


  result <-dplyr::bind_rows(output_list)
  if (all(result[["Missing"]] == 0))
    result <- result[, -which(names(result) == "Missing")]
  if (!full) {
    result <- result[, -2]
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
  attr(result, "paragraph") <- generate_paragraph(xvars, output_list)
  if (tableOnly) {
    return(result)
  }
  nicetable <- outTable(result, caption = caption, nicenames = nicenames)
  attr(nicetable, "paragraph") <- generate_paragraph(xvars, output_list)
  return(nicetable)
}
