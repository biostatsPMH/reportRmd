
mean_by_grp <- function(grp_level, data, xvar, grp, digits = 1) {
  x_var <- data[[xvar]]
  group_var <- data[[grp]]
  subset_grp <- subset(data, group_var == grp_level)
  new_xvar <- subset_grp[[xvar]]

  num_na <- sum(is.na(new_xvar))
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
  num_na <- sum(is.na(new_xvar))
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
      to_return[2, grp_level] <- "NA"
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
    return("NA")
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

is_binary <- function(x) all(unique(x) %in% c(0, 1))


