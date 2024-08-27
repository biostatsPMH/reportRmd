# Based on code in pivot_longer
selectMyVars <- function(data, grouping_var,x_vars,other_args=FALSE,...) {
  grouping_var <- tidyselect::eval_select(expr = enquo(grouping_var), data = data[unique(names(data))],
                                          allow_rename = FALSE)
  grp <- names(grouping_var)
  if (length(grp)>1) stop("Only one grouping variable is allowed")

  x_vars <- tidyselect::eval_select(expr = enquo(x_vars), data = data[unique(names(data))],
                                    allow_rename = FALSE)
  x_vars <- names(x_vars)
  cat("Grouping var: ",grp,"\n")
  cat("x variables: ",paste(x_vars,collapse = ", "),"\n")

  argList <- as.list(match.call()[-1])
  print(argList)
  print(all(sapply(argList$x_vars[-1], is.character)))
}

# pembrolizumab |> selectMyVars(grouping_var = cohort,x_vars = c(l_size,pdl1,age,change_ctdna_group))
# pembrolizumab |> selectMyVars(grouping_var = "cohort",x_vars = c("l_size",'pdl1','age','change_ctdna_group'))

rm_tidycompact <- function(data, xvars, grp, use_mean, caption = NULL, tableOnly = FALSE, covTitle = "", digits = 1, digits.cat = 0,  nicenames = TRUE, iqr = FALSE, all.stats = FALSE, pvalue = TRUE, effSize = FALSE, p.adjust = "none", unformattedp = FALSE, show.tests = FALSE, full = TRUE, percentage = "col") {
  if (missing(data))
    stop("data is a required argument")
  if (missing(xvars))
    stop("xvars is a required argument")
  if (!inherits(data, "data.frame"))
    stop("data must be supplied as a data frame.")
  # if (!inherits(xvars, "character"))
  #   stop("xvars must be supplied as a character vector or string indicating variables in data")

  grouping_var <- tidyselect::eval_select(expr = enquo(grp), data = data[unique(names(data))],
                                          allow_rename = FALSE)
  grouping_var <- names(grouping_var)
  if (length(grouping_var)>1) stop("Only one grouping variable is allowed")

  x_vars <- tidyselect::eval_select(expr = enquo(xvars), data = data[unique(names(data))],
                                    allow_rename = FALSE)
  x_vars <- names(x_vars)

  missing_vars = setdiff(x_vars, names(data))
  if (length(missing_vars) > 0) {
    stop(paste0("These xvars are not in the data: '", paste0(missing_vars, collapse = "', '"), "'"))
  }
  if (missing(use_mean)) {
    use_mean <- FALSE
  }
  # if (!missing(grp)) {
  #   if (!inherits(grp, "character") | length(grp) > 1)
  #     stop("grp must be supplied as a string indicating a variable in data")
  # }
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
  if (!all(sapply(argList$xvars[-1], is.character))) {
    argList$xvars <- x_vars
  }
  if (!is.character(argList$grp)) {
    argList$grp <- grouping_var
  }
  argsToPass <- intersect(names(formals(xvar_function)), names(argList))
  argsToPass <- setdiff(argsToPass,"xvars")
  args <- argList[argsToPass]
  grp <- grouping_var
  xvars <- x_vars
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
  if (all(na.omit(result[["Missing"]]) == 0))
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
  lbl <- result[, 1]
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
  if (covTitle != "") {
    names(result)[1] <- covTitle
  }
  else {
    names(result)[1] <- ""
  }
  attr(result, "description") <- generate_description(xvars, output_list)
  if (tableOnly) {
    return(result)
  }
  nicetable <- outTable(result, caption = caption, nicenames = nicenames)
  attr(nicetable, "description") <- generate_description(xvars, output_list)
  return(nicetable)
}
# pembrolizumab |> rm_tidycompact(xvars = c(age, change_ctdna_group, l_size, pdl1),
#                                 grp = sex, use_mean = "age", digits = c("age" = 2, "l_size" = 3),
#                                 digits.cat = 1, iqr = TRUE, show.tests = TRUE)
#
# pembrolizumab |> rm_tidycompact(xvars = c("age", "change_ctdna_group", "l_size", "pdl1"),
#                                 grp = "sex", use_mean = "age", digits = c("age" = 2, "l_size" = 3),
#                                 digits.cat = 1, iqr = TRUE, show.tests = TRUE)
