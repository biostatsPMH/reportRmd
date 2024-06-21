
mean_by_grp <- function(grp_level, data, xvar, grp, digits = 1) {
  x_var <- data[[xvar]]
  group_var <- data[[grp]]
  subset_grp <- subset(data, group_var == grp_level)
  new_xvar <- subset_grp[[xvar]]

  if (all(is.na((subset_grp[[xvar]])))) {
    return("NE")
  }
  return(paste0(format(round(mean(new_xvar, na.rm = TRUE), digits), nsmall = digits), " (", format(round(sd(new_xvar, na.rm = TRUE), digits), nsmall = digits), ")"))
}

median_by_grp <- function(grp_level, data, xvar, grp, iqr = FALSE, digits = 1, range_only = F) {
  x_var <- data[[xvar]]
  group_var <- data[[grp]]
  subset_grp <- subset(data, group_var == grp_level)
  new_xvar <- subset_grp[[xvar]]
  if (all(is.na((subset_grp[[xvar]])))) {
    return("NE")
  }
  if (iqr) {
    bracket <- paste0("(", format(round(stats::quantile(new_xvar, na.rm = TRUE, prob = 0.25), digits), nsmall = digits), ", ", format(round(stats::quantile(new_xvar, na.rm = TRUE, prob = 0.75), digits), nsmall = digits), ")")
  }
  else {
    bracket <- paste0("(", format(round(min(new_xvar, na.rm = TRUE), digits), nsmall = digits), ", ", format(round(max(new_xvar, na.rm = TRUE), digits), nsmall = digits), ")")
  }
  if (range_only && !iqr) {
    return(bracket)
  }
  return(paste0(format(round(median(new_xvar, na.rm = TRUE), digits), nsmall = digits), " ", bracket))
}

categ_xvar_helper <- function(xvar_level, data, xvar, grp, digits.cat = 0, percentage = "col") {
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
      to_return[2, grp_level] <- paste0(nrow(subset_xvar_grp), " (", ifelse((nrow(subset_grp) - missing_grp) == 0, 0, format(round((100*nrow(subset_xvar_grp) / (nrow(subset_xvar) - missing_xvar)), digits.cat), nsmall = digits.cat)), ")")
    }
    else {
      to_return[2, grp_level] <- paste0(nrow(subset_xvar_grp), " (", ifelse((nrow(subset_grp) - missing_grp) == 0, 0, format(round((100*nrow(subset_xvar_grp) / (nrow(subset_grp) - missing_grp)), digits.cat), nsmall = digits.cat)), ")")
    }
    if (all(is.na((subset_grp[[xvar]])))) {
      to_return[2, grp_level] <- "NE"
    }
  }
  to_return[1, "Full Sample"] <- paste0("Full Sample (n=", nrow(data), ")" )
  if (percentage == "row") {
    to_return[2, "Full Sample"] <- as.character(nrow(subset_xvar))
  }
  else {
    to_return[2, "Full Sample"] <- paste0(nrow(subset_xvar), " (", format(round((100*nrow(subset_xvar) / (nrow(data) - missing_xvar)), digits.cat), nsmall = digits.cat), ")")
  }
  return(to_return)
}

binary_xvar_helper <- function(grp_level, data, xvar, grp, digits.cat = 0, percentage = "col") {
  x_var <- data[[xvar]]
  group_var <- data[[grp]]
  subset_grp <- subset(data, group_var == grp_level)
  if (all(is.na((subset_grp[[xvar]])))) {
    return("NE")
  }
  num_missing <- sum(is.na(subset_grp[[xvar]]))
  if (percentage == "row") {
    bracket <- paste0(" (", format(round(100*sum(subset_grp[[xvar]], na.rm = TRUE) / sum(x_var, na.rm = TRUE), digits.cat), nsmall = digits.cat), ")")
  }
  else {
    bracket <- paste0(" (", format(round(100*sum(subset_grp[[xvar]], na.rm = TRUE) / (nrow(subset_grp) - num_missing), digits.cat), nsmall = digits.cat), ")")
  }
  return(paste0(sum(subset_grp[[xvar]], na.rm = TRUE), bracket))
}

is_binary <- function(x) all(unique(na.omit(x)) %in% c(0, 1))

generate_paragraph <- function(xvars, tables) {
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
  intro <- paste0("This table provides statistical summaries for covariates '", paste0(xvars, collapse = "', '"), "'.")
  ret <- intro
  if (length(names(descr_sum)) > 0) {
    for (x_sum in names(descr_sum)) {
      sum_stat <- paste0(x_sum, " statistical summaries are displayed for '", paste0(descr_sum[[x_sum]], collapse = "', '"), "' variable(s).")
      ret <- paste(ret, sum_stat)
    }
  }
  for (x_stat in names(descr_stat)) {
    test <- paste0("The '", x_stat, "' was used for covariate(s) '", paste(descr_stat[[x_stat]], collapse = "', '"), "'.")
    ret <- paste(ret, test)
  }
  if (length(names(descr_eff)) > 0) {
    for (x_eff in names(descr_eff)) {
      eff <- paste0("The '", paste(descr_eff[[x_eff]], collapse = "', '"), "' variable(s) use ", x_eff, " for effect size.")
      ret <- paste(ret, eff)
    }
  }
  return(ret)
}

