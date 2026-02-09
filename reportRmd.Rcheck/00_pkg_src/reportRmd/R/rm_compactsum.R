#'Output a compact summary table
#'
#'Outputs a table formatted for pdf, word or html output with summary statistics
#'
#'Comparisons for categorical variables default to chi-square tests, but if
#'there are counts of <5 then the Fisher Exact test will be used. For grouping
#'variables with two levels, either t-tests (mean) or wilcoxon tests (median)
#'will be used for numerical variables. Otherwise, ANOVA (mean) or Kruskal-
#'Wallis tests will be used. The statistical test used can be displayed by
#'specifying show.tests = TRUE. Statistical tests and effect sizes for grp and/
#'or xvars with less than 2 counts in any level will not be shown.
#'
#'Effect sizes are calculated as Cohen d for between group differences if the
#'variable is summarised with the mean, otherwise Wilcoxon R if summarised with
#'a median. Cramer's V is used for categorical variables, omega is used for
#'differences in means among more than two groups and epsilon for differences in
#'medians among more than two groups. Confidence intervals are calculated using
#'bootstrapping.
#'
#'tidyselect can only be used for xvars and grp arguments. Additional arguments
#'(digits, use_mean) must be passed in using characters if variable names are
#'used.
#'
#'@param data dataframe containing data
#'@param xvars character vector with the names of covariates to include in table
#'@param grp character with the name of the grouping variable
#'@param use_mean logical indicating whether mean and standard deviation will be
#'  returned for continuous variables instead of median. Otherwise, can specify
#'  for individual variables using a character vector containing the names of
#'  covariates to return mean and sd for (if use_mean is not supplied, all
#'  covariates will have median summaries). See examples.
#'@param caption character containing table caption (default is no caption)
#'@param tableOnly logical, if TRUE then a dataframe is returned, otherwise a
#'  formatted printed object is returned (default is FALSE)
#'@param covTitle character with the name of the covariate (predictor) column.
#'  The default is to leave this empty for output or, for table only output to
#'  use the column name 'Covariate'
#'@param digits numeric specifying the number of digits for summarizing mean
#'  data. Digits can be specified for individual variables using a named vector
#'  in the format digits=c("var1"=2,"var2"=3). If a variable is not in the
#'  vector the default will be used for it (default is 1). See examples
#'@param digits.cat numeric specifying the number of digits for the proportions
#'  when summarizing categorical data (default is 0)
#'@param nicenames logical indicating if you want to replace . and _ in strings
#'  .  with a space
#'@param iqr logical indicating if you want to display the interquartile range
#'  (Q1-Q3) as opposed to (min-max) in the summary for continuous variables
#'@param all.stats logical indicating if all summary statistics (Q1, Q3 + min,
#'  max on a separate line) should be displayed. Overrides iqr
#'@param pvalue logical indicating if you want p-values included in the table
#'@param effSize logical indicating if you want effect sizes and their 95%
#'  confidence intervals included in the table. Effect sizes calculated include
#'  Cramer's V for categorical variables, and Cohen's d, Wilcoxon r,
#'  Epsilon-squared, or Omega-squared for numeric/continuous variables
#'@param p.adjust p-adjustments to be performed
#'@param unformattedp logical indicating if you would like the p-value to be
#'  returned unformatted (ie. not rounded or prefixed with '<'). Best used with
#'  tableOnly = T and outTable function. See examples
#'@param show.sumstats logical indicating if the type of statistical summary
#'  (mean, median, etc) used should be shown.
#'@param show.tests logical indicating if the type of statistical test and
#'  effect size (if effSize = TRUE) used should be shown in a column beside the
#'  p-values.
#'@param full logical indicating if you want the full sample included in the
#'  table, ignored if grp is not specified
#'@param percentage choice of how percentages are presented, either column
#'  (default) or row
#'@returns A character vector of the table source code, unless tableOnly = TRUE
#'  in which case a data frame is returned. The output has the following
#'  attribute:
#'
#'   * "description", which describes what is included in the
#'  output table and the type of statistical summary for each covariate. When
#'  applicable, the types of statistical tests used will be included. If effSize
#'  = TRUE, the effect sizes for each covariate will also be mentioned.
#'
#'@references Smithson, M. (2002). Noncentral Confidence Intervals for
#'  Standardized Effect Sizes. (07/140 ed., Vol. 140). SAGE Publications.
#'  \doi{10.4135/9781412983761.n4}
#'@references Steiger, J. H. (2004). Beyond the F Test: Effect Size Confidence
#'  Intervals and Tests of Close Fit in the Analysis of Variance and Contrast
#'  Analysis. Psychological Methods, 9(2), 164–182.
#'  \doi{10.1037/1082-989X.9.2.164}
#'@references Kelley, T. L. (1935). An Unbiased Correlation Ratio Measure.
#'  Proceedings of the National Academy of Sciences - PNAS, 21(9), 554–559.
#'  \doi{10.1073/pnas.21.9.554}
#'@references Okada, K. (2013). Is Omega Squared Less Biased? A Comparison of
#'  Three Major Effect Size Indices in One-Way ANOVA. Behavior Research Methods,
#'  40(2), 129-147.
#'@references Breslow, N. (1970). A generalized Kruskal-Wallis test for
#'  comparing K samples subject to unequal patterns of censorship. Biometrika,
#'  57(3), 579-594.
#'@references FRITZ, C. O., MORRIS, P. E., & RICHLER, J. J. (2012). Effect Size
#'  Estimates: Current Use, Calculations, and Interpretation. Journal of
#'  Experimental Psychology. General, 141(1), 2–18. \doi{10.1037/a0024338}
#'
#' @examples
#' data("pembrolizumab")
#' rm_compactsum(data = pembrolizumab, xvars = c("age",
#' "change_ctdna_group", "l_size", "pdl1"), grp = "sex", use_mean = "age",
#' digits = c("age" = 2, "l_size" = 3), digits.cat = 1, iqr = TRUE,
#' show.tests = TRUE)
#'
#' # Other Examples (not run)
#' ## Include the summary statistic in the variable column
#' #rm_compactsum(data = pembrolizumab, xvars = c("age",
#' #"change_ctdna_group"), grp = "sex", use_mean = "age", show.sumstats=TRUE)
#'
#' ## To show effect sizes
#' #rm_compactsum(data = pembrolizumab, xvars = c("age",
#' #"change_ctdna_group"), grp = "sex", use_mean = "age", digits = 2,
#' #effSize = TRUE, show.tests = TRUE)
#'
#' ## To return unformatted p-values
#' #rm_compactsum(data = pembrolizumab, xvars = c("l_size",
#' #"change_ctdna_group"), grp = "cohort", effSize = TRUE, unformattedp = TRUE)
#'
#' ## Using tidyselect
#' #pembrolizumab |> rm_compactsum(xvars = c(age, sex, pdl1), grp = cohort,
#' #effSize = TRUE)
#'
#'@export
rm_compactsum <- function(data, xvars, grp, use_mean, caption = NULL, tableOnly = FALSE, covTitle = "", digits = 1, digits.cat = 0,  nicenames = TRUE, iqr = TRUE, all.stats = FALSE, pvalue = TRUE, effSize = FALSE, p.adjust = "none", unformattedp = FALSE, show.sumstats =FALSE,show.tests = FALSE, full = TRUE, percentage = "col") {
  if (missing(data))
    stop("data is a required argument")
  if (!inherits(data, "data.frame"))
    stop("data must be supplied as a data frame.")
  if (missing(xvars))
    stop("xvars is a required argument")

  x_vars <- tidyselect::eval_select(expr = tidyselect::enquo(xvars), data = data[unique(names(data))],
                                    allow_rename = FALSE)
  x_vars <- names(x_vars)
  if (!missing(grp)) {
    grouping_var <- tidyselect::eval_select(expr = tidyselect::enquo(grp), data = data[unique(names(data))],
                                            allow_rename = FALSE)
    grouping_var <- names(grouping_var)
  }
  else {
    grouping_var <- NULL
  }
  argList <- as.list(match.call(expand.dots = TRUE)[-1])
  argsToPass <- intersect(names(formals(xvar_function)), names(argList))
  argsToPass <- setdiff(argsToPass,"xvars")
  args <- argList[argsToPass]
  xvars <- x_vars
  dt_msg <- FALSE

  if (!missing(grp) && length(grouping_var)>1) stop("Only one grouping variable is allowed")

  if (all.stats) {
    use_mean <- FALSE
  }

  missing_vars = setdiff(c(grouping_var,x_vars), names(data))
  if (length(missing_vars) > 0) {
    stop(paste0("These xvars are not in the data: '", paste0(missing_vars, collapse = "', '"), "'"))
  }
  all_miss <- sapply(x_vars, function(x) all(is.na(data[[x]])))
  all_miss <- names(all_miss)[all_miss]
  if (length(all_miss)>0){
    warning(paste("The following variables have only missing data and will be omitted from the output:",paste0(all_miss,collapse=", ")))
    x_vars <- setdiff(x_vars,all_miss)
  }
  if (!missing(grp)) {
    if (!(grouping_var %in% names(data))) {
      stop(paste("grp variable",grouping_var,"is not in the data"))
    }
    if (grouping_var %in% x_vars){
      warning(paste(grouping_var,'is the grouping variable and can not appear as a covariate. \n',
                    'It is omitted from the output.'))
      x_vars <- setdiff(x_vars, grouping_var)
    }
    data[[grouping_var]] <- factor(data[[grouping_var]])
  }
  if (!missing(grp) & (effSize | show.tests | pvalue)) {
    grp_tab <- table(data[[grouping_var]])
    if (any(grp_tab < 2)) {
      warning("Small counts in '", grouping_var, "'. No statistical tests or effect size will be reported.")
      args$effSize = FALSE
      args$pvalue = FALSE
      args$show.tests = FALSE
    }
    else if (any(grp_tab < 5)) {
      warning("Small sample size in '", grouping_var, "' group may lead to unstable effect sizes.")
    }
  }
  if (missing(use_mean)) {
    use_mean <- FALSE
  } else {
    if (!inherits(use_mean,"logical")){
      if (!all(use_mean %in% names(data))){
        stop("use_mean must be logical or a character vector of variables in data")
      }
      if (!all(use_mean %in% x_vars)){
        stop("use_mean must be logical or a character vector of variables specified in x_vars")
      }
      for (xvar in use_mean){
        if (!is.numeric(data[[xvar]]))
          stop(paste0("use_mean can only be specified for numeric columns. Variable ",xvar," is not numeric."))
      }
    }
  }
  if (!is.numeric(digits)) {
    stop("digits must be a single numeric or a vector of digits")
  }
  if (length(digits)>1){
    if (!all(names(digits) %in% names(data)))
      stop(paste("The following variables were specified in digits but are not found in the data:", setdiff(names(digits),names(data))))
  }
  if (!(is.numeric(digits.cat) & length(digits.cat) == 1)) {
    stop("digits.cat must be a single numeric value")
  }
  if (!(percentage %in% c("col", "row"))) {
    stop("percentage argument must be either 'row' or 'col'")
  }

  data <- dplyr::select(data, dplyr::all_of(c(grouping_var,x_vars)))
  data_lbls <- extract_labels(data)

  if (!missing(grp) ) {
    n_na <- sum(is.na(data[[grouping_var]]))
    if (n_na>0) message(paste("There are",n_na,"observations missing",grouping_var,"\nThese will be removed from the data."))
    miss_data <- data
    data <- data[!is.na(data[[grouping_var]]),]
  }

  for (xvar in xvars) {
    if ( inherits(data[[xvar]], "POSIXt")) {
      data[[xvar]] <- as.Date(data[[xvar]])
    }
    if (inherits(data[[xvar]], "Date"))  {
      dt_msg <- TRUE
    }
    if (inherits(data[[xvar]],"difftime")){
      data[[xvar]] <- as.numeric(data[[xvar]])
    }
    if (is.character(data[[xvar]]) | is.logical(data[[xvar]])) {
      data[[xvar]] <- as.factor(data[[xvar]])
    }
    if (inherits(data[[xvar]],"haven_labelled")){
      lbl <- attr(data[[xvar]],"labels")
      attributes(data[[xvar]]) <- NULL
      newx <- factor(data[[xvar]],levels=lbl,labels=names(lbl))
    }
  }
  if (dt_msg) message("no statistical tests will be applied to date variables, date variables will be summarised with median")

  args$grp <- grouping_var
  args$data <- data
  if (tableOnly) {
    if (covTitle=="") args$covTitle <- "Covariate"  }

  output_list <- NULL

  for (xvar in xvars) {
    class(xvar) <- c(class(xvar),assign_method(data,xvar,use_mean))
    args$digits <- assign_digits(xvar,digits)
    args$xvar = xvar
    output_list[[xvar]] <- do.call(xvar_function, args)
  }
  result <- dplyr::bind_rows(output_list,.id="var")
  if (nrow(result)==0) return()

  if (all(na.omit(result[["Missing"]]) == 0))
    result <- result[, -which(names(result) == "Missing")]

  if ("p-value" %in% colnames(result)) {
    if (!pvalue) {
      result <- result[, -which(names(result) == "p-value")]
    }
    else {
      method <- p.adjust
      result[["p-value"]] <- p.adjust(result[["p-value"]], method = method)
      small_p <- which(result$`p-value`<.05)
      if (!unformattedp) {
        result[["p-value"]] <- formatp(result[["p-value"]])
      }
    }
  }
  to_indent <- which(result$Covariate != result$var & result$disp=="")
  if (all(result$disp=="")) {
     bold_cells <- cbind(which(result$Covariate == result$var ), rep(1, length(which(result$Covariate == result$var))))
  } else   bold_cells <- cbind(which(result$Covariate == result$var & result$disp!="" ), rep(1, length(which(result$Covariate == result$var))))
  result$var <- NULL
  lbl <- result[, 1]
  data <- set_labels(data,data_lbls)
  if (nicenames) {
    result[, 1] <- replaceLbl(data, lbl)
  }
  if ("ref" %in% names(result)) {
    result[, 1] <- paste0(result[, 1],ifelse(is.na(result$`ref`),"",result$`ref`))
    result$`ref` <- NULL
  }
  if (show.sumstats)  result[, 1] <- paste0(result[, 1],result$`disp`)
  result$`disp` <- NULL
  if (covTitle != "") {
    names(result)[1] <- covTitle
  }
  else {
    names(result)[1] <- ""
  }
  if (!full) {
    result <- result[, -2]
  }
  attr(result, "description") <- generate_description(xvars, output_list)
  if (!missing(grp)){
    if (n_na>0) {
      attr(result, "description") <- paste(n_na,"entries were missing data on",grouping_var,"and were removed prior to analysis.",attr(result, "description"))
    }
  }
  attr(result,'to_indent') <- to_indent
  attr(result,'bold_cells') <- bold_cells
  attr(result, "dimchk") <- dim(result)

  if (tableOnly) {
    if (names(result)[1]=="") names(result)[1] <- "Covariate"
    return(result)
  }
  if ("p-value" %in% colnames(result)) {
    p_col <- which(names(result) == "p-value")
    bold_cells <- rbind(bold_cells, cbind(small_p, rep(p_col, length(small_p))))
  }
  if (interactive()) cat(attr(result, "description"))
  outTable(result, caption = caption, nicenames = nicenames, to_indent = to_indent, bold_cells = bold_cells)
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
#' @param digits.cat numeric specifying the number of digits for the proportions
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
#'   p-values. Ignored if pvalue = FALSE
#' @param percentage choice of how percentages are presented, either column
#'   (default) or row
#' @return A data frame is returned
#'
#' @keywords internal

xvar_function <- function(xvar, data, grp, covTitle = "", digits = 1, digits.cat = 0, iqr = TRUE, all.stats = FALSE, pvalue = TRUE, effSize = FALSE, show.tests = FALSE, percentage = "col") {
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

xvar_function.rm_date <- function(xvar, data, grp, covTitle = "", digits = 1, digits.cat = 0, iqr = TRUE, all.stats = FALSE, pvalue = TRUE, effSize = FALSE, show.tests = FALSE, percentage = "col") {
  # no testing implemented for dates, all stats not run, just summarised with median and either iqr or range
  pvalue=FALSE; effSize=FALSE; show.tests = FALSE

  class(xvar) <- "character"
  x_var <- data[[xvar]]
  df <- data.frame(Covariate = xvar)
  df[["disp"]] <- ifelse(!iqr, " Median (Min-Max)", " Median (Q1-Q3)")

  date_sum <- function(x_var){
    var_sum <- as.character(summary(x_var))
    if (iqr){
      x_iqr <- var_sum[c(2,5)]
      bracket <- paste0("(", x_iqr[1], " to ",x_iqr[2], ")")
    } else {
      x_rng <- var_sum[c(1,6)]
      bracket <- paste0("(", x_rng[1], " to ",x_rng[2], ")")
    }
    return(paste0(var_sum[3], " ", bracket))
  }
  df[, paste0("Full Sample (n=", nrow(data), ")")] <- date_sum(x_var)

  if (!missing(grp)) {
    group_var <- factor(data[[grp]])

    grp_columns <- lapply(levels(group_var), function(g) date_sum(x_var[group_var==g]))
    i = 1
    for (level in levels(group_var)) {
      sub <- subset(data, group_var == level)
      title <- paste0(level, " (n=", nrow(sub), ")")
      df[, title] <- grp_columns[i]
      i = i+1
    }

    no_na_data <- na.omit(data.frame(x_var = x_var, group_var = group_var))
    if (length(unique(na.omit(group_var))) < 2 | length(unique(no_na_data$group_var)) < 2) {
      df[1, "Missing"] <- sum(is.na(x_var))
      return(df)
    }

  }
  df[1, "Missing"] <- sum(is.na(x_var))
  if (all.stats){
    df$Covariate <- df$disp
    df <- dplyr::bind_rows(data.frame(Covariate = xvar),df)
  }
  if (iqr) {
    attr(df, "stat_sum") <- "median (IQR)"
  }  else {
    attr(df, "stat_sum") <- "median (min/max)"
  }
  return(df)
}


xvar_function.rm_binary <- function(xvar, data, grp, covTitle = "", digits = 1, digits.cat = 0, iqr = TRUE, all.stats = FALSE, pvalue = TRUE, effSize = FALSE, show.tests = FALSE, percentage = "col") {
  if (!(pvalue | effSize)) {
    show.tests = FALSE
  }
  class(xvar) <- "character"
  df <- data.frame(Covariate = xvar)
  df[["disp"]] <- " n (%)"
  x_var <- data[[xvar]]
  if (percentage == "row") {
    df[, paste0("Full Sample (n=", nrow(data), ")")] <- as.character(sum(x_var, na.rm = TRUE))
  }
  else {
    df[, paste0("Full Sample (n=", nrow(data), ")")] <- paste0(sum(x_var, na.rm = TRUE), " (", format(round((100*sum(x_var, na.rm = TRUE) / (nrow(data) - sum(is.na(x_var)))), digits.cat), nsmall = digits.cat), "%)")
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
    if (pvalue | effSize | show.tests) {
      if (!any(chi_test_rm(cont_table)$expected < 5)) {
        chisq_test <- chi_test_rm(cont_table)
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
      else if (any(chi_test_rm(cont_table)$expected < 5)) {
        fisher_test <- fisher_test_rm(cont_table)
        df[1, "p-value"] <-fisher_test$p.value
        if (effSize) {
          output <- calc_CramerV(fisher_test)
          df[1, "Effect Size (95% CI)"] <- format_delta(output)
          attr(df, "eff_size") <- "Cramer's V"
        }
        if (show.tests & pvalue) {
          df[1, "pTest"] <- fisher_test$p_type
        }
        if (pvalue) {
          attr(df, "stat_test") <- fisher_test$stat_test
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

xvar_function.rm_mean <- function(xvar, data, grp, covTitle = "", digits = 1, digits.cat = 0, iqr = TRUE, all.stats = FALSE, pvalue = TRUE, effSize = FALSE, show.tests = FALSE, percentage = "col") {
  if (!(pvalue | effSize)) {
    show.tests = FALSE
  }
  class(xvar) <- "character"
  df <- data.frame(Covariate = xvar)
  df[["disp"]] <- " Mean (sd)"
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
    if (pvalue | effSize | show.tests) {
      if (length(unique(no_na_data$group_var)) == 2) {
        t_test <- t_test_rm(x_var, group_var)
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
      else if (length(unique(no_na_data$group_var)) > 2) {
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
      } else { # fewer than two groups, no comparisons
        if (interactive()) message(paste("Only one group with non-missing data for",xvar,"no statistical tests performed."))
      }
    }
  }
  df[, "Missing"] <- sum(is.na(x_var))
  attr(df, "stat_sum") <- "mean (sd)"
  return(df)
}

xvar_function.rm_median <- function(xvar, data, grp, covTitle = "", digits = 1, digits.cat = 0, iqr = TRUE, all.stats = FALSE, pvalue = TRUE, effSize = FALSE, show.tests = FALSE, percentage = "col") {
  if (!(pvalue | effSize)) {
    show.tests = FALSE
  }
  class(xvar) <- "character"
  x_var <- data[[xvar]]
  if (all.stats) {
    df <- data.frame(Covariate = c(xvar, "  Mean (sd)", "  Median (Q1-Q3)", "  Range (min-max)"))
    df[["disp"]] <- ""
    display_mean <- paste0(format(round(mean(x_var, na.rm = TRUE), digits), nsmall = digits), " (", format(round(sd(x_var, na.rm = TRUE), digits), nsmall = digits), ")")
    x_iqr <- stats::quantile(x_var, na.rm = TRUE, prob = c(0.25,0.75))
    bracket_iqr <- paste0("(", format(round(x_iqr[1], digits), nsmall = digits), "-", ifelse(x_iqr[2]<0," (",""),format(round(x_iqr[2], digits), nsmall = digits), ifelse(x_iqr[2]<0,")",""), ")")
    x_rng <- range(x_var, na.rm = TRUE)
    bracket_range <- paste0("(", format(round(x_rng[1], digits), nsmall = digits), "-", ifelse(x_rng[2]<0," (",""), format(round(x_rng[2], digits), nsmall = digits), ifelse(x_rng[2]<0,")",""), ")")
    df[2, paste0("Full Sample (n=", nrow(data), ")")] <- display_mean
    df[3, paste0("Full Sample (n=", nrow(data), ")")] <- paste0(format(round(median(x_var, na.rm = TRUE), digits), nsmall = digits), " ", bracket_iqr)
    df[4, paste0("Full Sample (n=", nrow(data), ")")] <- bracket_range
  }
  else {
    df <- data.frame(Covariate = xvar)
    df[["disp"]] <- ifelse(!iqr, " Median (Min-Max)", " Median (Q1-Q3)")
    if (iqr){
      x_iqr <- stats::quantile(x_var, na.rm = TRUE, prob = c(0.25,0.75))
      bracket <- paste0("(", format(round(x_iqr[1], digits), nsmall = digits), "-", ifelse(x_iqr[2]<0,"(",""), format(round(x_iqr[2], digits), nsmall = digits), ifelse(x_iqr[2]<0,")",""), ")")
    } else {
      x_rng <- range(x_var, na.rm = TRUE)
      bracket <- paste0("(", format(round(x_rng[1], digits), nsmall = digits), "-", ifelse(x_rng[2]<0,"(",""), format(round(x_rng[2], digits), nsmall = digits), ifelse(x_rng[2]<0,")",""), ")")
    }
    df[, paste0("Full Sample (n=", nrow(data), ")")] <- paste0(format(round(median(x_var, na.rm = TRUE), digits), nsmall = digits), " ", bracket)
  }

  if (!missing(grp)) {
    group_var <- data[[grp]]
    if (all.stats) {
      grp_mean <- lapply(as.list(levels(group_var)), mean_by_grp, data = data, xvar = xvar, grp = grp, digits = digits)
      grp_med <- lapply(as.list(levels(group_var)), median_by_grp, data = data, xvar = xvar, grp = grp, iqr = T, digits = digits)
      grp_range <- lapply(as.list(levels(group_var)), median_by_grp, data = data, xvar = xvar, grp = grp, iqr = F, digits = digits, range_only = T)
      i = 1
      for (level in levels(group_var)) {
        sub <- subset(data, group_var == level)
        title <- paste0(level, " (n=", nrow(sub), ")")
        df[2, title] <- grp_mean[i]
        df[3, title] <- grp_med[i]
        df[4, title] <- grp_range[i]
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
    no_na_data <- na.omit(data.frame(x_var = x_var, group_var = group_var))
    if (length(unique(na.omit(group_var))) < 2 | length(unique(no_na_data$group_var)) < 2) {
      df[1, "Missing"] <- sum(is.na(x_var))
      return(df)
    }
    no_na <- subset(data, !is.na(data[[xvar]]))
    no_na_tab <- table(no_na[[grp]])
    if (pvalue | effSize | show.tests) {
      if (length(unique(no_na_data$group_var)) == 2) {
        wilcox_test <- wilcox_test_rm(x_var, group_var)
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
      else if (length(unique(no_na_data$group_var)) > 2) {
        kruskal_test <- kruskal_test_rm(x_var, group_var)
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
      } else { # fewer than two groups, no comparisons
        if (interactive()) message(paste("Only one group with non-missing data for",xvar,"no statistical tests performed."))
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

xvar_function.rm_categorical <- function(xvar, data, grp, covTitle = "", digits = 1, digits.cat = 0, iqr = TRUE, all.stats = FALSE, pvalue = TRUE, effSize = FALSE, show.tests = FALSE, percentage = "col") {
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
  attr(df, "stat_sum") <- "counts (%)"
  i = 2
  for (xvar_level in levels(x_var)) {
    xvar_subset <- subset(data, x_var == xvar_level)
    if (percentage == "row") {
      df[i, paste0("Full Sample (n=", nrow(data), ")")] <- as.character(nrow(xvar_subset))
    }
    else {
      df[i, paste0("Full Sample (n=", nrow(data), ")")] <- paste0(nrow(xvar_subset), " (", format(round((100*nrow(xvar_subset) / (nrow(data) - sum(is.na(x_var)))), digits.cat), nsmall = digits.cat), "%)")
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
      if (interactive()) message(paste("Only non-missing data for one group, no statistical test for",xvar,"performed."))
      return(df)
    }
    no_na <- subset(data, !is.na(data[[xvar]]))
    no_na_tab <- table(no_na[[grp]])
    # if (any(no_na_tab < 2)) {
    #   effSize <- FALSE
    #   pvalue <- FALSE
    #   show.tests <- FALSE
    # }
    if (pvalue | effSize | show.tests) {
      if (!any(chi_test_rm(cont_table)$expected < 5)) {
        chisq_test <- chi_test_rm(cont_table)
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
      else if (any(chi_test_rm(cont_table)$expected < 5)) {
        fisher_test <- fisher_test_rm(cont_table)
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
          attr(df, "stat_test") <- "Fisher's exact test"
        }
        if (show.tests & effSize) {
          df[1, "effStat"] <- "Cramer's V"
        }
      }
    }
  }
  df[1, "Missing"] <- sum(is.na(x_var))
  return(df)
}

xvar_function.rm_two_level <- function(xvar, data, grp, covTitle = "", digits = 1, digits.cat = 0, iqr = TRUE, all.stats = FALSE, pvalue = TRUE, effSize = FALSE, show.tests = FALSE, percentage = "col") {
  if (!(pvalue | effSize)) {
    show.tests = FALSE
  }
  class(xvar) <- "character"
  temp <- data.frame()
  x_var <- data[[xvar]]

  if (length(unique(na.omit(data[[xvar]]))) == 2) {
    unique_levels <- unique(x_var)
    unique_levels <- sort(unique_levels)
    binary_column <- ifelse(x_var == unique_levels[1], 0, 1)
    level_shown <- unique_levels[2]
  }
  else { # only one unique value in x_var
    level_shown <- unique(x_var)[1]
    binary_column <- ifelse(x_var == level_shown, 1, 0)
  }

  if (!missing(grp)) {
    temp <- data[, grp,drop=FALSE]
    temp[[xvar]] <- binary_column
  }
  else {
    temp <- data[, xvar,drop=FALSE]
    temp[[xvar]] <- binary_column
  }
  df <- data.frame(Covariate = xvar)
  df[["ref"]] <-  paste(" -", level_shown)
  df[["disp"]] <-  "n (%)"
  x_var <- temp[[xvar]]
  if (percentage == "row") {
    df[, paste0("Full Sample (n=", nrow(temp), ")")] <- as.character(sum(x_var, na.rm = TRUE))
  }
  else {
    df[, paste0("Full Sample (n=", nrow(temp), ")")] <- paste0(sum(x_var, na.rm = TRUE), " (", format(round((100*sum(x_var, na.rm = TRUE) / (nrow(data) - sum(is.na(x_var)))), digits.cat), nsmall = digits.cat), "%)")
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
      if (!any(chi_test_rm(cont_table)$expected < 5)) {
        chisq_test <- chi_test_rm(cont_table)
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
      else if (any(chi_test_rm(cont_table)$expected < 5)) {
        fisher_test <- fisher_test_rm(cont_table)
        df[1, "p-value"] <-fisher_test$p.value
        if (effSize) {
          output <- calc_CramerV(fisher_test)
          df[1, "Effect Size (95% CI)"] <- format_delta(output)
          attr(df, "eff_size") <- "Cramer's V"
        }
        if (show.tests & pvalue) {
          df[1, "pTest"] <- fisher_test$p_type
        }
        if (pvalue) {
          attr(df, "stat_test") <- fisher_test$stat_test
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

assign_method <- function(data,xvar,use_mean){
  if (grepl(class(data[[xvar]]),"factor") & length(unique(na.omit(data[[xvar]]))) <= 2) {
    return("rm_two_level")
  }
  else if (inherits(data[[xvar]],"factor")) {
    return("rm_categorical")
  }
  else if (inherits(data[[xvar]],"Date") ) {
    return( "rm_date")
  }
  else if (is.numeric(data[[xvar]]) && is_binary(data[[xvar]])) {
    return( "rm_binary")
  }
  else if (is.numeric(data[[xvar]])) {
    if (identical(use_mean, FALSE)) {
      return("rm_median")
    }
    else {
      if (grepl(class(use_mean), "character")) {
        if (xvar %in% use_mean) {
          return("rm_mean")
        }
        else {
          return("rm_median")
        }
      }
      else if (grepl(class(use_mean), "logical")) {
        if (identical(use_mean, TRUE)) {
          return("rm_mean")
        }
        else {
          return("rm_median")
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
}

assign_digits <- function(xvar,digits){
  if (length(digits)==1) return(digits)
  if (xvar %in% names(digits)) return(digits[[xvar]])
  # otherwise default to 1
  return(1)
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
    }}  else stop("htest must be the output of a call to t.test, chisq.test, kruskal.test, aov or fisher_test_rm")
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

fisher_test_rm <- function(x,...){
  chi.out <- suppressWarnings(stats::chisq.test(x,...))
  rtn <- try(stats::fisher.test(x,...),silent=T)
  if (inherits(rtn,"try-error")){
    message("Using MC sim. Use set.seed() prior to function for reproducible results.")
    rtn <- stats::fisher.test(x, simulate.p.value = T)
    rtn$p_type <- "MC sim"
    rtn$stat_test <- "Fisher's exact test with Monte Carlo simulation"
  } else {
    rtn$p_type <- "Fisher Exact"
    rtn$stat_test <- "Fisher's exact test"
  }
  rtn$observed <- chi.out$observed
  rtn$statistic <- chi.out$statistic
  rtn$parameter <- chi.out$parameter
  rtn$k <- min(dim(chi.out$observed), na.rm = TRUE)
  rtn$N <- sum(chi.out$observed)
  rtn$method <- paste(rtn$method,"with additional chisq.test arguments")
  class(rtn) <- c(class(rtn),"fisher_test_rm")
  return(rtn)
}

chi_test_rm <- function(x,...){
  chi_test <- suppressWarnings(stats::chisq.test(x,...))
  chi_test$k <- min(dim(chi_test$observed), na.rm = TRUE)
  chi_test$N <- sum(chi_test$observed)
  class(chi_test) <- c(class(chi_test),"chi_test_rm")
  return(chi_test)
}

t_test_rm <- function(xvar,grp){
  t_test <- stats::t.test(xvar~grp)
  n <- unlist(lapply(levels(grp),function(g){
    length(na.omit(xvar[grp==g]))
  }))
  names(n) <- levels(grp)
  t_test$n <- n
  t_test$xvar <- xvar
  t_test$grp <- grp
  class(t_test) <- c(class(t_test),"t_test_rm")
  return(t_test)
}

wilcox_test_rm <- function(xvar,grp){
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
  class(wilcox_test) <- c(class(wilcox_test),"wilcox_test_rm")
  return(wilcox_test)
}

kruskal_test_rm <- function(xvar,grp){
  kruskal_test <- stats::kruskal.test(xvar~grp)
  n <- length(na.omit(xvar))
  kruskal_test$n <- n
  kruskal_test$xvar <- xvar
  kruskal_test$grp <- grp
  class(kruskal_test) <- c(class(kruskal_test),"kruskal_test_rm")
  return(kruskal_test)
}

anova_toOmegaSq <- function(anova_summary){
  tbl <- anova_summary[[1]]
  num <- tbl[1,2]-tbl[1,1]*tbl[2,3]
  den <- tbl[1,2]+tbl[2,2]+tbl[2,3]
  return(num/den)
}

anova_toEpsilonSq <- function(anova_summary){
  tbl <- anova_summary[[1]]
  num <- tbl[1,2]-tbl[1,1]*tbl[2,3]
  den <- tbl[1,2]+tbl[2,2]
  return(num/den)
}

chi_toCramer <- function(chisq_test){
  obs <- chisq_test$observed
  V = sqrt(chisq_test$statistic/(sum(obs)*(min(dim(obs))-1)))
  return(V)
}

t_toCohen <-function(t_test){
  if (!inherits(t_test,"t_test_rm")) stop("t_test must be a function returned from t_test_rm")
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
    new_chi <- chi_test_rm(dt[[1]],df[[2]])

    chi_toCramer(new_chi)
  }
  b_cramer <- boot::boot(df,bs_cramer,R=1000)
  b_ci <- boot::boot.ci(b_cramer,conf = CIwidth,type="basic")
  output = c("cramer v"=cramer,lower=b_ci$basic[4],upper=b_ci$basic[5])
  return(output)
}

calc_cohenD <- function(t_test, CIwidth = 0.95) {
  cohen <- t_toCohen(t_test)
  df <- data.frame(xvar = t_test$xvar, grp = t_test$grp)

  bs_cohen <- function(data,indices){
    dt <- data[indices,]
    new_t <- tryCatch(t_test_rm(dt$xvar,dt$grp))
    t_toCohen(new_t)
  }

  b_cohen <- boot::boot(df,bs_cohen,R=1000)
  b_ci <- boot::boot.ci(b_cohen,conf = CIwidth,type="basic")
  output = c("cohen d"=cohen,lower=b_ci$basic[4],upper=b_ci$basic[5])
  return(output)
}

calc_WilcoxonR <- function(wilcox_test, CIwidth = 0.95) {
  r <- wilcox_effSize(wilcox_test)
  df <- data.frame(xvar = wilcox_test$xvar, grp = wilcox_test$grp)

  bs_r <- function(data,indices){
    dt <- data[indices,]
    new_wilcox <- wilcox_test_rm(dt$xvar, dt$grp)
    wilcox_effSize(new_wilcox)
  }
  b_r <- boot::boot(df,bs_r,R=1000)
  b_ci <- boot::boot.ci(b_r,conf = CIwidth,type="basic")
  output = c("wilcoxon r"=r,lower=b_ci$basic[4],upper=b_ci$basic[5])
  return(output)
}

calc_epsilonSq <- function(kruskal_test, CIwidth = 0.95) {
  eps <- chi_toEpsilonSq(kruskal_test)
  df <- data.frame(xvar = kruskal_test$xvar, grp = kruskal_test$grp)

  bs_epsilon <- function(data,indices){
    dt <- data[indices,]
    new_kruskal <- kruskal_test_rm(dt$xvar,dt$grp)

    chi_toEpsilonSq(new_kruskal)
  }
  b_epsilon <- boot::boot(df,bs_epsilon,R=1000)
  b_ci <- boot::boot.ci(b_epsilon,conf = CIwidth,type="basic")
  output = c("epsilon sq"=eps,lower=b_ci$basic[4],upper=b_ci$basic[5])
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
  b_ci <- boot::boot.ci(b_omega,conf = CIwidth,type="basic")
  output = c("omega squared"=omega,lower=b_ci$basic[4],upper=b_ci$basic[5])
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



mean_by_grp <- function(grp_level, data, xvar, grp, digits = 1) {
  x_var <- data[[xvar]]
  group_var <- data[[grp]]
  subset_grp <- subset(data, group_var == grp_level)
  new_xvar <- subset_grp[[xvar]]

  if (all(is.na((subset_grp[[xvar]])))) {
    return("NA")
  }
  return(paste0(format(round(mean(new_xvar, na.rm = TRUE), digits), nsmall = digits), " (", format(round(sd(new_xvar, na.rm = TRUE), digits), nsmall = digits), ")"))
}

median_by_grp <- function(grp_level, data, xvar, grp, iqr = FALSE, digits = 1, range_only = F) {
  x_var <- data[[xvar]]
  group_var <- data[[grp]]
  subset_grp <- subset(data, group_var == grp_level)
  new_xvar <- subset_grp[[xvar]]
  if (all(is.na((subset_grp[[xvar]])))) {
    return("NA")
  }
  if (iqr) {
    x_iqr <- stats::quantile(new_xvar, na.rm = TRUE, prob = c(0.25,0.75))
    bracket <- paste0("(", format(round(x_iqr[1], digits), nsmall = digits), "-", ifelse(x_iqr[2]<0,"(",""),format(round(x_iqr[2], digits), nsmall = digits), ifelse(x_iqr[2]<0,")",""),")")
  }
  else {
    x_rng <- range(new_xvar,na.rm=T)
    bracket <- paste0("(", format(round(x_rng[1], digits), nsmall = digits), "-", ifelse(x_rng[2]<0,"(",""),format(round(x_rng[2], digits), nsmall = digits), ifelse(x_rng[2]<0,")",""),")")
  }
  if (range_only && !iqr) {
    return(bracket)
  }
  return(paste0(format(round(median(new_xvar, na.rm = TRUE), digits), nsmall = digits), " ", bracket))
}

categ_xvar_helper <- function(xvar_level, data, xvar, grp, digits.cat = 0, percentage = "col",pct_symbol=TRUE) {
  if (pct_symbol) r_brckt <- "%)" else r_brckt <- ")"
  x_var <- data[[xvar]]
  group_var <- data[[grp]]
  nlevels_grp <- length(levels(group_var))
  to_return <- matrix(nrow = 2, ncol = nlevels_grp + 1)
  colnames(to_return) <- c("Full Sample", levels(group_var))

  subset_xvar <- subset(data, x_var == xvar_level)
  missing_xvar <- sum(is.na(x_var))
  for (grp_level in levels(group_var)) {
    subset_grp <- subset(data, group_var == grp_level)
    missing_grp <- sum(is.na(subset_grp[[xvar]]))
    subset_xvar_grp <- subset(data, group_var == grp_level & x_var == xvar_level)
    to_return[1, grp_level] <- paste0(grp_level, " (n=", nrow(subset_grp), ")")
    if (percentage == "row") {
      to_return[2, grp_level] <- paste0(nrow(subset_xvar_grp), " (", ifelse((nrow(subset_grp) - missing_grp) == 0, 0, format(round((100*nrow(subset_xvar_grp) / (nrow(subset_xvar) - missing_xvar)), digits.cat), nsmall = digits.cat)), r_brckt)
    }
    else {
      to_return[2, grp_level] <- paste0(nrow(subset_xvar_grp), " (", ifelse((nrow(subset_grp) - missing_grp) == 0, 0, format(round((100*nrow(subset_xvar_grp) / (nrow(subset_grp) - missing_grp)), digits.cat), nsmall = digits.cat)), r_brckt)
    }
    if (all(is.na((subset_grp[[xvar]])))) {
      to_return[2, grp_level] <- "NA"
    }
  }
  to_return[1, "Full Sample"] <- paste0("Full Sample (n=", nrow(data), ")" )
  if (percentage == "row") {
    to_return[2, "Full Sample"] <- as.character(nrow(subset_xvar))
  }
  else {
    to_return[2, "Full Sample"] <- paste0(nrow(subset_xvar), " (", format(round((100*nrow(subset_xvar) / (nrow(data) - missing_xvar)), digits.cat), nsmall = digits.cat), r_brckt)
  }
  return(to_return)
}

binary_xvar_helper <- function(grp_level, data, xvar, grp, digits.cat = 0, percentage = "col",pct_symbol=TRUE) {
  if (pct_symbol) r_brckt <- "%)" else r_brckt <- ")"
  x_var <- data[[xvar]]
  group_var <- data[[grp]]
  subset_grp <- subset(data, group_var == grp_level)
  if (all(is.na((subset_grp[[xvar]])))) {
    return("NA")
  }
  num_missing <- sum(is.na(subset_grp[[xvar]]))
  if (percentage == "row") {
    bracket <- paste0(" (", format(round(100*sum(subset_grp[[xvar]], na.rm = TRUE) / sum(x_var, na.rm = TRUE), digits.cat), nsmall = digits.cat), r_brckt)
  }
  else {
    bracket <- paste0(" (", format(round(100*sum(subset_grp[[xvar]], na.rm = TRUE) / (nrow(subset_grp) - num_missing), digits.cat), nsmall = digits.cat), r_brckt)
  }
  return(paste0(sum(subset_grp[[xvar]], na.rm = TRUE), bracket))
}

is_binary <- function(x) {
  attributes(x) <- NULL
  all(unique(na.omit(x)) %in% c(0, 1))
}

format_strings <- function(variables) {
  if (length(variables) == 1) {
    return(variables)
  } else if (length(variables) == 2) {
    return(paste(variables[1], "and", variables[2]))
  } else {
    last_var <- variables[length(variables)]
    first_vars <- paste(variables[-length(variables)], collapse = ", ")
    return(paste0(first_vars, ", and ", last_var))
  }
}

generate_description <- function(xvars, tables) {
  descr_sum <- list()
  for (xvar in xvars) {
    x_sum <- attr(tables[[xvar]], "stat_sum")
    if (!is.null(x_sum)) {
      if (!(x_sum %in% names(descr_sum))) {
        descr_sum[[x_sum]] <- xvar
      }
      else {
        descr_sum[[x_sum]] <- c(descr_sum[[x_sum]], xvar)
      }
    }
  }
  descr_eff <- list()
  for (xvar in xvars) {
    x_eff <- attr(tables[[xvar]], "eff_size")
    if (!is.null(x_eff)) {
      if (!(x_eff %in% names(descr_eff))) {
        descr_eff[[x_eff]] <- xvar
      }
      else {
        descr_eff[[x_eff]] <- c(descr_eff[[x_eff]], xvar)
      }
    }
  }
  descr_stat <- list()
  for (xvar in xvars) {
    x_stat <- attr(tables[[xvar]], "stat_test")
    if (!is.null(x_stat)) {
      if (!(x_stat %in% names(descr_stat))) {
        descr_stat[[x_stat]] <- xvar
      }
      else {
        descr_stat[[x_stat]] <- c(descr_stat[[x_stat]], xvar)
      }
    }
  }
  ret <- ""
  if (length(names(descr_sum)) > 0) {
    sum_sent <- c()
    for (x_sum in names(descr_sum)) {
      summ<- paste0(x_sum, " for ", format_strings(descr_sum[[x_sum]]))
      sum_sent <- c(sum_sent, summ)
    }
    sum_sent <- paste0("Descriptive statistics were calculated as ", format_strings(sum_sent), ".")
    ret <- paste(ret, sum_sent)
  }

  if (length(names(descr_stat)) > 0) {
    stat_sent <- c()
    for (x_stat in names(descr_stat)) {
      stat <- paste0(x_stat, " for ", format_strings(descr_stat[[x_stat]]))
      stat_sent <- c(stat_sent, stat)
    }
    stat_sent <- paste0("Between group comparisons were made using ", format_strings(stat_sent), ".")
    ret <- paste(ret, stat_sent)
  }

  if (length(names(descr_eff)) > 0) {
    eff_sent <- c()
    for (x_eff in names(descr_eff)) {
      eff <- paste0(x_eff, " for ", format_strings(descr_eff[[x_eff]]))
      eff_sent <- c(eff_sent, eff)
    }
    eff_sent <- paste0("Reported effect sizes are ", format_strings(eff_sent), ". 1000 bootstrap samples were used to calculate effect size confidence intervals.")
    ret <- paste(ret, eff_sent)
  }
  ret <- trimws(ret, "both")
  return(ret)
}
