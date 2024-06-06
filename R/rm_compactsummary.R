#' Output a compact summary table
#'
#' Returns a data frame corresponding to a descriptive table.
#'
#' Comparisons for categorical variables default to chi-square tests, but if
#' there are counts of <5 then the Fisher Exact test will be used. For grouping
#' variables with two levels, either t-tests (mean) or wilcoxon tests (median)
#' will be used for numerical variables. Otherwise, ANOVA (mean) or kruskal-
#' wallis tests will be used. The statistical test used can be displayed by
#' specifying show.tests = TRUE.
#'
#'
#' @param data dataframe containing data
#' @param xvars character vector with the names of covariates to include in
#'   table
#' @param grp character with the name of the grouping variable
#' @param use_mean logical indicating whether mean and standard deviation will
#'   be returned for continuous variables instead of median. Otherwise, can
#'   specify for individual variabels using a character vector containing the
#'   names of covariates to return mean and sd for (default is FALSE).
#'   See examples
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
#'   confidence intervals included in the table. Can only be obtained if pvalue
#'   is also requested. Effect sizes calculated include Cramer's V for
#'   categorical variables, and Cohen's d, Wilcoxon r, Epsilon-squared, or
#'   Omega-squared for numeric/continuous variables
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
#' @return A character vector of the table source code, unless tableOnly = TRUE
#' in which case a data frame is returned
#'
#' @references Smithson, M. (2002). Noncentral CIwidthidence Intervals for
#'   Standardized Effect Sizes. In CIwidthidence Intervals (07/140 ed., Vol.
#'   140). SAGE Publications. https://doi.org/10.4135/9781412983761.n4
#' @references Steiger, J. H. (2004). Beyond the F Test: Effect Size Confidence
#'   Intervals and Tests of Close Fit in the Analysis of Variance and Contrast
#'   Analysis. Psychological Methods, 9(2), 164–182.
#'   https://doi.org/10.1037/1082-989X.9.2.164
#' @references Kelley, T. L. (1935). An Unbiased Correlation Ratio Measure.
#'   Proceedings of the National Academy of Sciences - PNAS, 21(9), 554–559.
#'   https://doi.org/10.1073/pnas.21.9.554
#' @references Okada, K. (2013). Is Omega Squared Less Biased? A Comparison of
#'   Three Major Effect Size Indices in One-Way ANOVA. Behavior Research
#'   Methods, 40(2), 129-147.
#' @references Breslow, N. (1970). A generalized Kruskal-Wallis test for
#' comparing K samples subject to unequal patterns of censorship. Biometrika,
#' 57(3), 579-594.
#' @references FRITZ, C. O., MORRIS, P. E., & RICHLER, J. J. (2012). Effect Size
#' Estimates: Current Use, Calculations, and Interpretation. Journal of
#' Experimental Psychology. General, 141(1), 2–18.
#' https://doi.org/10.1037/a0024338
#'
#' @export
#' @examples
#' data("pembrolizumab")
#' rm_compactsummary(data = pembrolizumab, xvars = c("age",
#' "change_ctdna_group", "l_size", "pdl1"), grp = "sex", use_mean = "age",
#' digits = c("age" = 2, "l_size" = 3), digits.cat = 1, iqr = TRUE,
#' show.tests = TRUE)
#'
#' # To Show Effect Sizes
#' rm_compactsummary(data = pembrolizumab, xvars = c("age",
#' "change_ctdna_group"), grp = "cohort", use_mean = "age", digits = 2,
#' effSize = TRUE, show.tests = TRUE)
#'
#' # To return unformatted p-values
#' rm_compactsummary(data = pembrolizumab, xvars = c("l_size",
#' "change_ctdna_group"), grp = "cohort", effSize = TRUE, unformattedp = TRUE)
#'
#' @export
rm_compactsummary <- function(data, xvars, grp, use_mean, caption = NULL, tableOnly = FALSE, covTitle = "", digits = 1, digits.cat = 0,  nicenames = TRUE, iqr = FALSE, all.stats = FALSE, pvalue = TRUE, effSize = FALSE, p.adjust = "none", unformattedp = FALSE, show.tests = FALSE, full = TRUE, percentage = "col") {
  argList <- as.list(match.call(expand.dots = TRUE)[-1])
  argsToPass <- intersect(names(formals(xvar_function)), names(argList))
  argsToPass <- setdiff(argsToPass,"xvars")
  args <- argList[argsToPass]

  if (!pvalue) {
    args$effSize <- FALSE
    args$show.tests <- FALSE
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
    else if (is.numeric(data[[xvar]]) && is_binary(na.omit(data[[xvar]]))) {
      class(xvar) <- c(class(xvar), "rm_binary")
    }
    else if (is.numeric(data[[xvar]])) {
      if (missing(use_mean)) {
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
  if (tableOnly) {
    return(result)
  }
  return(outTable(result, caption = caption, nicenames = nicenames))
}
