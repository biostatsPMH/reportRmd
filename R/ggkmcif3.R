# Utility function for default values
`%||%` <- function(x, y) if (is.null(x)) y else x

# Internal helper: format "value(lower-upper)" CI string
# @keywords internal
# @noRd
format_ci <- function(value, lower, upper, digits) {
  paste0(round_sprintf(value, digits = digits),
         "(", round_sprintf(lower, digits = digits),
         "-", round_sprintf(upper, digits = digits), ")")
}

# Data Preparation Functions ----

#' Validate and prepare input data
#' @param data Input dataframe
#' @param response Character vector with time and status column names
#' @param cov Covariate column name (optional)
#' @param print.n.missing Whether to print missing data message
validate_and_prepare_data <- function(data, response, cov = NULL, print.n.missing = TRUE) {
  data <- data.frame(data)

  # Identify missing data
  missing_cols <- if (!is.null(cov)) c(response, cov) else response
  remove <- rowSums(is.na(data[, missing_cols])) > 0

  # Print missing data message
  if (print.n.missing && sum(remove) > 0) {
    msg <- if (sum(remove) == 1) {
      "1 observation with missing data was removed."
    } else {
      paste(sum(remove), "observations with missing data were removed.")
    }
    message(msg)
  }

  data[!remove, ]
}

#' Process covariate variable (factor conversion, numeric cutoffs)
#' @param data Input data
#' @param cov Covariate column name
#' @param cut Numeric cutoff for continuous variables
#' @param stratalabs Custom strata labels
process_covariate <- function(data, cov, cut = NULL, stratalabs = NULL) {
  if (is.null(cov)) return(list(data = data, stratalabs = NULL, numeric = FALSE))

  # Coerce to factor if needed
  if (!is.factor(data[[cov]]) && !is.numeric(data[[cov]])) {
    message("Coercing the cov variable to factor")
    data[[cov]] <- factor(data[[cov]])
  }

  numeric <- is.numeric(data[[cov]])

  # Handle numeric covariates
  if (numeric) {
    if (is.null(cut)) cut <- median(data[[cov]], na.rm = TRUE)

    le_code <- "\u2264"
    levels_labels <- c(paste0(le_code, round(cut, 2)), paste0(">", round(cut, 2)))

    data[[cov]] <- factor(
      ifelse(data[[cov]] <= cut, levels_labels[1], levels_labels[2]),
      levels = levels_labels
    )

    if (is.null(stratalabs)) stratalabs <- levels_labels
  }

  # Drop unused levels and set labels
  data[[cov]] <- droplevels(data[[cov]])
  if (is.null(stratalabs)) stratalabs <- levels(data[[cov]])

  list(data = data, stratalabs = stratalabs, numeric = numeric)
}

# Survival Analysis Functions ----

#' Fit Kaplan-Meier survival curves
#' @param data Input data
#' @param response Time and status variables
#' @param cov Covariate (optional)
#' @param conf.type Confidence interval type
fit_km_model <- function(data, response, cov = NULL, conf.type = "log") {
  formula_str <- paste("survival::Surv(", response[1], ",", response[2], ")")
  rhs <- if (is.null(cov)) "1" else cov
  formula <- as.formula(paste(formula_str, "~", rhs))

  sfit <- survival::survfit(formula, data = data, conf.type = conf.type)

  # Calculate p-value if multiple groups
  pval <- if (!is.null(cov)) {
    sdiff <- survival::survdiff(formula, data = data)
    pchisq(sdiff$chisq, length(sdiff$n) - 1, lower.tail = FALSE)
  } else NA

  list(sfit = sfit, pval = pval)
}

#' Fit Competing Risks Model
#'
#' @param data Input dataframe
#' @param response Character vector with time and status column names
#' @param cov Covariate column name (optional)
#' @return List containing fit object and group separator
fit_cif_model <- function(data, response, cov = NULL) {
  if (is.null(cov)) {
    # Single group analysis
    invisible(utils::capture.output(
      fit <- cmprsk::cuminc(data[[response[1]]], data[[response[2]]])
    ))
    gsep <- " "
  } else {
    # Multiple group analysis
    newgpvar <- paste0(data[[cov]], ":")
    newgpvar <- factor(newgpvar, levels = paste0(levels(data[[cov]]), ":"))
    invisible(utils::capture.output(
      fit <- cmprsk::cuminc(data[[response[1]]], data[[response[2]]], newgpvar)
    ))
    gsep <- ": "
  }

  list(fit = fit, gsep = gsep)
}


# CIF-specific Functions ----

#' Calculate Median Time to Event for CIF
#'
#' @param fit CIF fit object
#' @param event_name Name of the event in fit object
#' @return Median time to event
calculate_cif_median <- function(fit, event_name) {
  estall <- fit[[event_name]]$est
  timall <- fit[[event_name]]$time
  timall[estall >= 0.5][1]
}

#' Process CIF Median Values
#'
#' @param fit CIF fit object
#' @param plot.event Events to plot
#' @param stratalabs Strata labels
#' @param median.lines Whether to calculate for median lines
#' @param median.text Whether to add median text
#' @param median.digits Number of digits for median
#' @param multiple_lines Whether there are multiple strata
#' @return List with updated stratalabs and median values
process_cif_medians <- function(fit, plot.event, stratalabs,
                                median.lines = FALSE, median.text = FALSE,
                                median.digits = 3, multiple_lines = FALSE) {

  median_vals <- NULL
  median_txt <- NULL

  if (!median.lines && !median.text) {
    return(list(stratalabs = stratalabs, median_vals = NULL, median_txt = NULL))
  }

  # Get event names that match plot.event
  last_character <- substr(names(fit), nchar(names(fit)), nchar(names(fit)))
  get_values <- names(fit)[last_character %in% plot.event]

  if (length(get_values) > 0) {
    median_vals <- sapply(get_values, function(x) calculate_cif_median(fit, x))

    if (!multiple_lines) {
      median_txt <- round_sprintf(median_vals, digits = median.digits)
    } else if (median.text && length(plot.event) == 1) {
      stratalabs <- paste(stratalabs, ", Median=",
                          round_sprintf(median_vals, digits = median.digits))
    }
  }

  list(stratalabs = stratalabs, median_vals = median_vals, median_txt = median_txt)
}

#' Calculate CIF Estimates at Specific Time Points
#'
#' @param fit CIF fit object
#' @param plot.event Events to plot
#' @param set.time Time points to evaluate
#' @param set.time.CI Whether to include confidence intervals
#' @param set.time.digits Number of digits
#' @param multiple_lines Whether there are multiple strata
#' @return Data frame with time-specific estimates
calculate_cif_timepoints <- function(fit, plot.event, set.time,
                                     set.time.CI = FALSE, set.time.digits = 3,
                                     multiple_lines = FALSE) {

  last_character <- substr(names(fit), nchar(names(fit)), nchar(names(fit)))
  get_values <- names(fit)[last_character %in% plot.event]

  set.surv.text <- NULL
  set.surv <- NULL

  for (time_i in set.time) {
    z <- qnorm(1 - (1 - 0.95) / 2)

    val <- cmprsk::timepoints(fit, times = time_i)$est[
      rownames(cmprsk::timepoints(fit, times = time_i)$est) %in% get_values, ]
    var <- cmprsk::timepoints(fit, times = time_i)$var[
      rownames(cmprsk::timepoints(fit, times = time_i)$est) %in% get_values, ]

    lower <- val^exp(-z * sqrt(var) / (val * log(val)))
    upper <- val^exp(z * sqrt(var) / (val * log(val)))

    keep_sum <- data.frame(value = val, time = time_i)
    keep_sum$set.CI <- if (set.time.CI) {
      format_ci(val, lower, upper, set.time.digits)
    } else {
      round_sprintf(val, digits = set.time.digits)
    }

    if (is.null(set.surv.text)) {
      set.surv.text <- paste0(time_i, " ", keep_sum$set.CI)
    } else {
      set.surv.text <- paste0(set.surv.text, ",", time_i, " ", keep_sum$set.CI)
    }

    set.surv <- rbind(keep_sum, set.surv)
  }

  list(set.surv = set.surv, set.surv.text = set.surv.text)
}

#' Process CIF Time-Specific Estimates
#'
#' @param fit CIF fit object
#' @param plot.event Events to plot
#' @param stratalabs Strata labels
#' @param set.time.text Text label for time points
#' @param set.time Time points to evaluate
#' @param set.time.line Whether to add lines
#' @param set.time.CI Whether to include confidence intervals
#' @param set.time.digits Number of digits
#' @param multiple_lines Whether there are multiple strata
#' @return List with updated stratalabs and time-specific estimates
process_cif_timepoints <- function(fit, plot.event, stratalabs,
                                   set.time.text = NULL, set.time = NULL,
                                   set.time.line = FALSE, set.time.CI = FALSE,
                                   set.time.digits = 3, multiple_lines = FALSE) {

  if (is.null(set.time.text) && !set.time.line) {
    return(list(stratalabs = stratalabs, set.surv = NULL, set.surv.text = NULL))
  }

  if (is.null(set.time)) set.time <- 5

  result <- calculate_cif_timepoints(fit, plot.event, set.time,
                                     set.time.CI, set.time.digits,
                                     multiple_lines)

  # Add to stratalabs if appropriate
  if (!is.null(set.time.text)) {
    if (multiple_lines && length(plot.event) == 1) {
      stratalabs <- paste0(stratalabs, ", ", result$set.surv.text)
    } else if (!multiple_lines) {
      # For single group, format the text appropriately
      formatted_text <- sapply(set.time, function(t) {
        idx <- which(result$set.surv$time == t)
        if (length(idx) > 0) {
          paste0(t, " ", set.time.text, "=", result$set.surv$set.CI[idx])
        }
      })
      result$set.surv.text <- paste(formatted_text, collapse = ",")
    }
  }

  list(stratalabs = stratalabs,
       set.surv = result$set.surv,
       set.surv.text = result$set.surv.text)
}

#' Create Survival Fit for Risk Table in CIF
#'
#' @param data Input data
#' @param response Time and status variables
#' @param cov Covariate (optional)
#' @return Survival fit object for risk table
create_cif_risk_table_sfit <- function(data, response, cov = NULL) {
  # Create temporary data where all events are treated the same
  temp <- data
  temp[[response[2]]][temp[[response[2]]] > 0] <- 1

  formula_str <- paste("survival::Surv(", response[1], ",", response[2], ")")
  rhs <- if (is.null(cov)) "1" else cov
  formula <- as.formula(paste(formula_str, "~", rhs))

  survival::survfit(formula, data = temp)
}




# HR and Statistical Functions ----

# Internal helper: append HR values and p-values to strata labels given a fitted model.
# @param compact If TRUE, use compact formatting (no spaces) for plot annotations.
# @keywords internal
# @noRd
append_hr_labels <- function(stratalabs, fit, HR, HR_pval, HR.digits,
                             HR.pval.digits, compact = FALSE) {
  hr_prefix <- if (compact) "HR=" else "HR = "
  HR_vals <- paste0(hr_prefix, sapply(seq(length(stratalabs) - 1), function(i) {
    psthr(summary(fit)$conf.int[i, c(1, 3, 4)], y = HR.digits, compact = compact)
  }))
  if (HR) stratalabs[-1] <- paste(stratalabs[-1], HR_vals)
  if (HR_pval) {
    stratalabs[-1] <- paste(stratalabs[-1],
                            sapply(summary(fit)$coef[, 5], lpvalue2, digits = HR.pval.digits))
  }
  stratalabs[1] <- paste(stratalabs[1], "REF")
  stratalabs
}

#' Add hazard ratios to strata labels
#' @param stratalabs Original strata labels
#' @param data Input data
#' @param response Time and status variables
#' @param cov Covariate
#' @param type Model type ("KM" or "CIF")
#' @param plot.event Event for CIF
#' @param HR Whether to include HR
#' @param HR_pval Whether to include HR p-value
#' @param HR.digits Number of digits for HR
#' @param HR.pval.digits Number of digits for HR p-value
add_km_hazard_ratios <- function(stratalabs, data, response, cov, type, plot.event = 1,
                              HR = FALSE, HR_pval = FALSE, HR.digits = 2, HR.pval.digits = 3) {

  if (!HR && !HR_pval) return(stratalabs)
  if (is.null(cov)) return(stratalabs)

  multiple_lines <- length(stratalabs) > 1
  if (!multiple_lines) return(stratalabs)

  if (type == "KM") {
    coxfit <- survival::coxph(
      as.formula(paste(paste("survival::Surv(", response[1], ",", response[2], ")", sep = ""), "~", cov, sep = "")),
      data = data
    )
    stratalabs <- append_hr_labels(stratalabs, coxfit, HR, HR_pval, HR.digits, HR.pval.digits)

  } else if (type == "CIF" && length(plot.event) == 1 && plot.event[1] == 1) {
    crrfit <- crrRx(
      as.formula(paste(paste(response, collapse = "+"), "~", cov, sep = "")),
      data = data
    )
    stratalabs <- append_hr_labels(stratalabs, crrfit, HR, HR_pval, HR.digits, HR.pval.digits)
  }

  stratalabs
}

#' Add CIF Hazard Ratios to Strata Labels
#'
#' @param stratalabs Original strata labels
#' @param data Input data
#' @param response Time and status variables
#' @param cov Covariate
#' @param plot.event Event for CIF (must be 1 for HR calculation)
#' @param HR Whether to include HR
#' @param HR_pval Whether to include HR p-value
#' @param HR.digits Number of digits for HR
#' @param HR.pval.digits Number of digits for HR p-value
#' @return Updated strata labels
add_cif_hazard_ratios <- function(stratalabs, data, response, cov,
                                  plot.event = 1, HR = FALSE, HR_pval = FALSE,
                                  HR.digits = 2, HR.pval.digits = 3) {

  if (!HR && !HR_pval) return(stratalabs)
  if (is.null(cov)) return(stratalabs)
  if (length(stratalabs) <= 1) return(stratalabs)
  if (length(plot.event) != 1 || plot.event[1] != 1) return(stratalabs)

  # Fit competing risks regression
  crrfit <- crrRx(
    as.formula(paste(paste(response, collapse = "+"), "~", cov)),
    data = data
  )

  append_hr_labels(stratalabs, crrfit, HR, HR_pval, HR.digits, HR.pval.digits,
                   compact = TRUE)
}
# Summary Statistics Functions ----

#' Calculate median survival times and add to labels
#' @param sfit Survival fit object (for KM)
#' @param fit CIF fit object
#' @param stratalabs Strata labels
#' @param type Model type
#' @param plot.event Events to plot (for CIF)
#' @param median.text Whether to add median text
#' @param median.CI Whether to include CI
#' @param median.digits Number of digits
calculate_and_add_median_times <- function(sfit = NULL, fit = NULL, stratalabs, type = "KM",
                                           plot.event = 1, median.text = FALSE,
                                           median.CI = FALSE, median.digits = 3) {

  if (!median.text || type == "CIF") return(list(stratalabs = stratalabs, median_vals = NULL))

  multiple_lines <- length(stratalabs) > 1

  if (type == "KM") {
    if (multiple_lines) {
      median_vals <- summary(sfit)$table[, "median"]
      median_lower <- summary(sfit)$table[, "0.95LCL"]
      median_upper <- summary(sfit)$table[, "0.95UCL"]

      med_txt <- if (median.CI) {
        format_ci(median_vals, median_lower, median_upper, median.digits)
      } else round_sprintf(median_vals, digits = median.digits)

      stratalabs <- paste(stratalabs, ", Median=", med_txt)

    } else {
      median_vals <- summary(sfit)$table["median"]
      median_lower <- summary(sfit)$table["0.95LCL"]
      median_upper <- summary(sfit)$table["0.95UCL"]

      median_txt <- if (median.CI) {
        format_ci(median_vals, median_lower, median_upper, median.digits)
      } else round_sprintf(median_vals, digits = median.digits)
    }
  }

  list(stratalabs = stratalabs, median_vals = median_vals %||% NULL)
}

#' Calculate survival/CIF at specific time points and add to labels
#' @param sfit Survival fit object (for KM)
#' @param fit CIF fit object
#' @param stratalabs Strata labels
#' @param type Model type
#' @param plot.event Events to plot (for CIF)
#' @param set.time.text Text label for time points
#' @param set.time Time points to evaluate
#' @param set.time.line boolean to specify if you want the survival added as
#'   lines to the plot at a specified point
#' @param set.time.CI Whether to include CI
#' @param set.time.digits Number of digits
calculate_and_add_time_specific_estimates <- function(sfit = NULL, fit = NULL, stratalabs,
                                                      type = "KM", plot.event = 1,
                                                      set.time.text = NULL, set.time = NULL,
                                                      set.time.line = FALSE,
                                                      set.time.CI = FALSE, set.time.digits = 3) {

  if (is.null(set.time.text) && !set.time.line) {
    return(list(stratalabs = stratalabs, set.surv = NULL))
  }

  if (is.null(set.time)) set.time <- 5

  multiple_lines <- length(stratalabs) > 1

  if (type == "KM") {
    if (multiple_lines) {
      final.set.text <- NULL
      set.surv <- NULL

      for (time_i in set.time) {
        sum <- summary(sfit, time = c(0, time_i))
        df_text <- data.frame(strata = sum$strata, time = sum$time)
        df_text$set.CI <- if (set.time.CI) {
          format_ci(sum$surv, sum$lower, sum$upper, set.time.digits)
        } else round_sprintf(sum$surv, digits = set.time.digits)

        dup_strata <- sum$strata[duplicated(sum$strata)]
        keep_sum <- df_text[!(sum$strata %in% dup_strata) | sum$time != 0, ]
        keep_sum$set.CI[keep_sum$time == 0] <- NA
        set.surv <- rbind(set.surv, keep_sum)

        if (is.null(final.set.text)) {
          final.set.text <- paste0(time_i, " ", set.time.text, "=", keep_sum$set.CI)
        } else {
          final.set.text <- paste0(final.set.text, ",", time_i, " ", set.time.text, "=", keep_sum$set.CI)
        }
      }

      if (!is.null(set.time.text)) {
        stratalabs <- paste0(stratalabs, ", ", final.set.text)
      }

    } else {
      set.surv.text <- NULL
      set.surv <- NULL

      for (time_i in set.time) {
        sum <- summary(sfit, time = c(0, time_i))
        df_text <- data.frame(time = sum$time)
        df_text$set.CI <- if (set.time.CI) {
          format_ci(sum$surv, sum$lower, sum$upper, set.time.digits)
        } else round_sprintf(sum$surv, digits = set.time.digits)

        keep_sum <- df_text
        if (length(sum$surv) > 1) keep_sum <- df_text[sum$time != 0, ]
        keep_sum$set.CI[keep_sum$time == 0] <- NA
        set.surv <- rbind(set.surv, keep_sum)

        if (is.null(set.surv.text)) {
          set.surv.text <- paste0(time_i, " ", set.time.text, "=", keep_sum$set.CI)
        } else {
          set.surv.text <- paste0(set.surv.text, ",", time_i, " ", set.time.text, "=", keep_sum$set.CI)
        }
      }
    }

  } else if (type == "CIF") {
    # Handle CIF time-specific estimates
    if (multiple_lines) {
      set.surv.text <- NULL
      set.surv <- NULL

      for (time_i in set.time) {
        z <- qnorm(1 - (1 - 0.95)/2)
        last_character <- substr(names(fit), nchar(names(fit)), nchar(names(fit)))
        get_values <- names(fit)[last_character %in% plot.event]

        val <- cmprsk::timepoints(fit, times = time_i)$est[rownames(cmprsk::timepoints(fit, times = time_i)$est) %in% get_values, ]
        var <- cmprsk::timepoints(fit, times = time_i)$var[rownames(cmprsk::timepoints(fit, times = time_i)$est) %in% get_values, ]
        lower <- val^exp(-z * sqrt(var)/(val * log(val)))
        upper <- val^exp(z * sqrt(var)/(val * log(val)))

        keep_sum <- data.frame(val, time_i)
        names(keep_sum)[1] <- "value"
        keep_sum$set.CI <- if (set.time.CI) {
          format_ci(val, lower, upper, set.time.digits)
        } else round_sprintf(val, digits = set.time.digits)

        if (is.null(set.surv.text)) {
          set.surv.text <- paste0(time_i, " ", set.time.text, "=", keep_sum$set.CI)
        } else {
          set.surv.text <- paste0(set.surv.text, ",", time_i, " ", set.time.text, "=", keep_sum$set.CI)
        }
        set.surv <- rbind(keep_sum, set.surv)
      }

      if (!is.null(set.time.text) && length(plot.event) == 1) {
        stratalabs <- paste0(stratalabs, ", ", set.surv.text)
      }
    }
  }

  list(stratalabs = stratalabs, set.surv = set.surv)
}

# Data Frame Creation Functions ----

#' Create plotting dataframe for KM curves
#' @param sfit Survival fit object
#' @param stratalabs Strata labels
#' @param conf.curves Whether to include confidence intervals
#' @param conf.type Confidence interval type
create_km_dataframe <- function(sfit, stratalabs, conf.curves = FALSE, conf.type = "log") {
  multiple_lines <- length(stratalabs) > 1

  df <- data.frame(
    time = sfit$time,
    n.risk = sfit$n.risk,
    n.censor = sfit$n.censor,
    n.event = sfit$n.event,
    surv = sfit$surv,
    strata = if (multiple_lines) summary(sfit, censored = TRUE)$strata else factor("All")
  )

  # Add confidence intervals if requested
  if (conf.type != "none" && conf.curves) {
    df$upper <- sfit$upper
    df$lower <- sfit$lower
  } else {
    df$upper <- NA
    df$lower <- NA
  }

  # Set strata labels
  levels(df$strata) <- stratalabs

  # Add time zero point
  zeros <- data.frame(
    time = 0, surv = 1, upper = 1, lower = 1,
    strata = if (multiple_lines) levels(df$strata) else factor("All")
  )

  zeros <- match_column_types(zeros,df)
  dplyr::bind_rows(zeros, df) |>
    transform(strata = factor(strata, levels = stratalabs))
}

#' Create CIF Plotting Data Frame
#'
#' @param fit CIF fit object
#' @param gsep Group separator from fit_cif_model
#' @param plot.event Events to plot
#' @param stratalabs Strata labels
#' @param conf.type Confidence interval type
#' @param flip.CIF Whether to flip the CIF curve
#' @param eventlabs Event labels
#' @param cov Covariate for proper factor levels
#' @param data Original data (for factor levels)
#' @return Data frame for plotting
create_cif_dataframe <- function(fit, gsep, plot.event, stratalabs,
                                 conf.type = "log", flip.CIF = FALSE,
                                 eventlabs = NULL, cov = NULL, data = NULL) {

  # Extract test if present
  if (!is.null(fit$Tests)) {
    test <- fit$Tests
    fit <- fit[names(fit) != "Tests"]
    attr(fit, "test") <- test
  }

  # Process fit into dataframe
  fit2 <- lapply(fit, `[`, 1:3)
  gnames <- names(fit2)

  fit2_list <- lapply(seq_along(gnames), function(ind) {
    df <- as.data.frame(fit2[[ind]])
    df$name <- gnames[ind]
    df
  })

  df <- do.call(rbind, fit2_list)

  # Parse event and strata information
  df$event <- sapply(strsplit(df$name, split = gsep), `[`, 2)
  df$strata <- sapply(strsplit(df$name, split = gsep), `[`, 1)

  # Set strata factor levels
  if (!is.null(cov) && !is.null(data)) {
    df$strata <- factor(df$strata, levels = levels(data[[cov]]))
    levels(df$strata) <- stratalabs
    df$strata <- factor(df$strata, levels = stratalabs)
  } else {
    df$strata <- "ALL"
  }

  # Rename and calculate additional columns
  names(df)[names(df) == "est"] <- "surv"
  df$std <- sqrt(df$var)

  # Filter for requested events
  df <- df[df$event %in% plot.event, ]

  # Initialize confidence intervals
  df$upper <- NA
  df$lower <- NA

  # Calculate confidence intervals if requested
  if (conf.type == "log") {
    z <- qnorm(1 - (1 - 0.95) / 2)
    df$lower <- df$surv^exp(-z * sqrt(df$var) / (df$surv * log(df$surv)))
    df$upper <- df$surv^exp(z * sqrt(df$var) / (df$surv * log(df$surv)))
  } else if (conf.type != "none") {
    message("Only log confidence intervals available for CIF")
  }

  # Flip CIF if requested
  if (flip.CIF) {
    df$surv <- 1 - df$surv
    df$upper <- 1 - df$upper
    df$lower <- 1 - df$lower
  }

  # Apply event labels if provided
  if (!is.null(eventlabs)) {
    df$event <- factor(df$event)
    levels(df$event) <- eventlabs
  }

  # Preserve any test information
  if (exists("test")) {
    attr(df, "test") <- test
  }

  df
}
# Statistical Test Functions ----

#' Extract Gray's test results from CIF fit
#' @param fit CIF fit object or dataframe with test attribute
#' @param plot.event Events to test
extract_grays_test <- function(fit, plot.event = 1) {
  test <- NULL

  if (is.data.frame(fit)) {
    test <- attr(fit, "test")
  } else if ("Tests" %in% names(fit)) {
    test <- fit$Tests
  }

  if (is.null(test)) return(NA)

  if (length(plot.event) == 1) {
    test[rownames(test) == plot.event, 2]
  } else {
    test[, 2]
  }
}

# Plotting Helper Functions ----

#' Get current ggplot2 theme base size
#' @noRd
get_current_base_size <- function(default = 11) {
  current_theme <- ggplot2::theme_get()
  base_size <- current_theme$text$size
  if (is.null(base_size)) base_size <- default
  return(base_size)
}

#' Create base ggplot for survival curves
#' @param df Plotting dataframe
#' @param type Plot type ("KM" or "CIF")
#' @param xlab x axis label
#' @param ylab y axis label
#' @param multiple_lines Whether multiple strata
#' @param plot.event Events to plot
#' @param event How to distinguish events ("col" or "linetype")
#' @param lsize Line size
#' @param fsize Font size
#' @param col colours vector
#' @param linetype Line types vector
#' @param legend.pos Legend position
#' @param legend.title Legend title
#' @param times Time breaks
#' @param ylim Y-axis limits
#' @param xlim X-axis limits
#' @param main Plot title
create_base_plot <- function(df, type, xlab = "Time", ylab = "Survival Probability",
                             multiple_lines, plot.event = 1, event = "col", lsize = 0.5,
                             fsize, col, linetype = NULL, legend.pos = "bottom",
                             legend.title = NULL,
                             times, ylim = c(0, 1), xlim = NULL, main = NULL) {

  # Determine maximum x value
  maxxval <- max(times)
  maxxlim <- if (is.null(xlim)) maxxval else xlim[2]

  # Create base plot
  if (type == "CIF" && length(plot.event) > 1) {
    if (event == "linetype") {
      p <- ggplot2::ggplot(df) +
        ggplot2::geom_step(ggplot2::aes(time, surv, colour = strata, linetype = event), linewidth = lsize)
      if (!is.null(legend.title)) p <- p + ggplot2::guides(linetype=guide_legend(legend.title))
    } else if (event == "col") {
      p <- ggplot2::ggplot(df) +
        ggplot2::geom_step(ggplot2::aes(time, surv, colour = event, linetype = strata), linewidth = lsize)
      if (!is.null(legend.title)) p <- p + ggplot2::guides(colour=guide_legend(legend.title))
    }
  } else {
    p <- ggplot2::ggplot(df) +
      ggplot2::geom_step(ggplot2::aes(time, surv, group = strata, linetype = linetype, col = strata), linewidth = lsize)
    if (!is.null(legend.title)) p <- p + ggplot2::guides(colour=guide_legend(legend.title))
  }

  # Apply theme and scales
  p <- p +
    ggplot2::theme_classic(base_size = fsize) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(margin = ggplot2::margin(t = 0), vjust = 1),
      axis.text.x.top = ggplot2::element_text(margin = ggplot2::margin(b = 0), vjust = 0),
      axis.text.y = ggplot2::element_text(margin = ggplot2::margin(r = 0), hjust = 1),
      axis.text.y.right = ggplot2::element_text(margin = ggplot2::margin(l = 0), hjust = 0),
      axis.title.x.bottom = ggplot2::element_text(vjust = 4),
      panel.grid.minor = ggplot2::element_blank(),
      legend.key = ggplot2::element_rect(colour = "transparent", fill = "transparent"),
      legend.key.spacing.y = unit(-.25, "lines"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      legend.background = ggplot2::element_blank()
    ) +
    ggplot2::scale_x_continuous(paste0("\n", xlab), breaks = times, limits = c(0, maxxval)) +
    ggplot2::coord_cartesian(xlim = c(0, maxxlim)) +
    ggplot2::scale_y_continuous(paste0(ylab, "\n"), limits = ylim)
  if (is.null(legend.title)){
    p <- p+ ggplot2::theme(legend.title = element_blank())
  } else p <- p+ ggplot2::theme(legend.title = element_text(margin = margin(b = 0,l=5),hjust = 0.5))
  # Add title if specified
  if (!is.null(main)) {
    p <- p + ggplot2::ggtitle(main) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = fsize))
  }
  # position the legend
  if (is.numeric(legend.pos)) {
    p <- p +
      ggplot2::theme(
        legend.position = "inside",
        legend.position.inside = legend.pos,
        legend.justification.inside = legend.pos)
  } else {
    p <- p + ggplot2::theme(legend.position = legend.pos)
  }
  p
}

#' Add confidence bands to plot
#' @param p ggplot object
#' @param df Plotting dataframe
#' @param type Plot type
#' @param event Event distinction method
#' @param plot.event Events to plot
add_confidence_bands <- function(p, df, type, event = "col", plot.event = 1) {

  conf_data <- df[!is.na(df$upper) & !is.na(df$lower), ]
  if (nrow(conf_data) == 0) return(p)

  if (type == "KM") {
    p + ggplot2::geom_ribbon(
      data = conf_data,
      ggplot2::aes(x = time, fill = strata, ymin = lower, ymax = upper),
      inherit.aes = FALSE, alpha = 0.2, show.legend = FALSE
    )
  } else if (type == "CIF") {
    if (event == "linetype" || length(plot.event) == 1) {
      for (evnt in unique(df$event)) {
        event_data <- conf_data[conf_data$event == evnt, ]
        if (nrow(event_data) > 0) {
          p <- p + ggplot2::geom_ribbon(
            data = event_data,
            ggplot2::aes(x = time, fill = strata, ymin = lower, ymax = upper),
            inherit.aes = FALSE, alpha = 0.2, show.legend = FALSE
          )
        }
      }
    } else if (event == "col" && length(plot.event) > 1) {
      for (stra in unique(df$strata)) {
        strata_data <- conf_data[conf_data$strata == stra, ]
        if (nrow(strata_data) > 0) {
          p <- p + ggplot2::geom_ribbon(
            data = strata_data,
            ggplot2::aes(x = time, fill = event, ymin = lower, ymax = upper),
            inherit.aes = FALSE, alpha = 0.2, show.legend = FALSE
          )
        }
      }
    }
    p
  }
}

#' Add censor marks to KM plot
#' @param p ggplot object
#' @param df Plotting dataframe
#' @param censor.size Size of censor marks
#' @param censor.stroke Stroke of censor marks
#' @param shape Shape for censor marks, defaults to "|", but can use an
#'   character or standard geom_point shapes (0-24)
add_censor_marks <- function(p, df,censor.size = 0.5, censor.stroke = 1.5, shape="|") {

  if (is.null(shape)) return(p)

  censor_data <- subset(df, n.censor > 0)
  if (nrow(censor_data) == 0) return(p)

  if (is.numeric(shape)){
    if (!shape %in% 0:24){
      message("censor.shape must be a number 0-24 or a character")
      shape = 3
    }
    if (missing(censor.size)) censor.size <-2.5
  }
  p + ggplot2::geom_point(
    data = censor_data,
    ggplot2::aes(x = time, y = surv, group = strata, col = strata),
    shape = shape, size = censor.size, stroke = censor.stroke, show.legend = FALSE)
}

#' Add statistical test results to plot
#' @param p ggplot object
#' @param type Plot type
#' @param multiple_lines Whether multiple strata
#' @param pval_result P-value
#' @param pval.pos Position for p-value text
#' @param times Time breaks
#' @param xlim X-axis limits
#' @param ylim Y-axis limits
#' @param psize Text size for p-value
#' @param pval.digits Number of digits for p-value
#' @param plot.event Events being plotted
#' @param eventlabs Event labels
add_statistical_tests <- function(p, type, multiple_lines, pval_result, pval.pos = NULL,
                                  times, xlim, ylim, psize = 3.5, pval.digits = 3,
                                  plot.event = 1, eventlabs = NULL) {

  if (!multiple_lines || is.na(pval_result)) return(p)

  if (type == "KM") {
    pvaltxt <- paste(lpvalue2(pval_result, pval.digits), "(Log Rank)")

    if (is.null(pval.pos)) {
      p + ggplot2::annotate("text", x = Inf, y = Inf,
                            hjust = 1, vjust = 1,
                            label = pvaltxt, size = psize)
      # p + ggplot2::annotate("text", x = 0.85 * xlim[2], y = ylim[1],
      #                       label = pvaltxt, size = psize)
    } else {
      p + ggplot2::annotate("text", x = pval.pos[1], y = pval.pos[2],
                            label = pvaltxt, size = psize)
    }
  } else if (type == "CIF") {
    if (length(plot.event) == 1) {
      pvaltxt <- paste(lpvalue2(pval_result, pval.digits), "(Gray's test)")

      if (is.null(pval.pos)) {
        p + ggplot2::annotate("text", x = Inf, y = -Inf,
                              hjust = 1.1, vjust = -0.5,
                              label = pvaltxt, size = psize)
        # p + ggplot2::annotate("text", x = 0.85 * xlim[2], y = ylim[1],
        #                       label = pvaltxt, size = psize)
      } else {
        p + ggplot2::annotate("text", x = pval.pos[1], y = pval.pos[2],
                              label = pvaltxt, size = psize)
      }
    } else {
      pvaltxt <- sapply(pval_result, lpvalue2, pval.digits)
      pvaltxt <- c("Gray's test", paste(eventlabs %||% paste("Event", plot.event), pvaltxt))

      if (is.null(pval.pos)) {
        p + ggplot2::annotate("text", x = 0.85 * xlim[2],
                              y = c(0.12, 0.08, 0.04), label = pvaltxt, size = psize)
      } else {
        p + ggplot2::annotate("text", x = pval.pos[1],
                              y = c(pval.pos[2], pval.pos[2] - 0.04, pval.pos[2] - 0.08),
                              label = pvaltxt, size = psize)
      }
    }
  }
}

#' Add median survival text to plot
#' @param p ggplot object
#' @param type Plot type
#' @param multiple_lines Whether multiple strata
#' @param median_txt Median text
#' @param median.pos Position for median text
#' @param times Time breaks
#' @param ylim Y-axis limits
#' @param median.size Text size
#' @param plot.event Events being plotted
#' @param eventlabs Event labels
add_median_text <- function(p, type, multiple_lines, median_txt, median.pos = NULL,
                            times, ylim, median.size = 3, plot.event = 1, eventlabs = NULL) {

  if (multiple_lines || is.null(median_txt)) return(p)

  if (length(plot.event) == 2 && !is.null(eventlabs)) {
    median_txt <- paste(paste0(eventlabs, ":", median_txt), collapse = "\n")
  }

  if (is.null(median.pos)) {
    if (type == "KM") {
      median.pos <- c(0.1 * max(times), ylim[1])
    } else {
      median.pos <- c(0.1 * max(times), ylim[2] * 0.95)
    }
  }

  p + ggplot2::annotate("text", x = median.pos[1], y = median.pos[2],
                        label = median_txt, size = median.size)
}

#' Add survival at set time text to plot
#' @param p ggplot object
#' @param type Plot type
#' @param multiple_lines Whether multiple strata
#' @param set.surv.text Survival text
#' @param set.pos Position for survival text
#' @param times Time breaks
#' @param ylim Y-axis limits
#' @param set.size Text size
#' @param plot.event Events being plotted
#' @param eventlabs Event labels
add_set_time_text <- function(p, type, multiple_lines, set.surv.text, set.pos = NULL,
                              times, ylim, set.size = 3, plot.event = 1, eventlabs = NULL) {

  if (multiple_lines || is.null(set.surv.text)) return(p)

  if (length(plot.event) == 2 && !is.null(eventlabs)) {
    set.surv.text <- paste(paste0(eventlabs, ":", set.surv.text), collapse = "\n")
  }

  if (is.null(set.pos)) {
    if (type == "KM") {
      set.pos <- c(0.1 * max(times), ylim[1] + 0.1)
    } else {
      set.pos <- c(0.1 * max(times), ylim[2] * 0.85)
    }
  }

  p + ggplot2::annotate("text", x = set.pos[1], y = set.pos[2],
                        label = set.surv.text, size = set.size)
}

#' Add median survival lines to plot
#' @param p ggplot object
#' @param median_vals Median values
#' @param median.lsize Line size for median lines
add_median_lines <- function(p, median_vals, median.lsize = 1) {

  if (is.null(median_vals)) return(p)

  temp <- data.frame(x = median_vals, y = 0.5)
  p <- p + ggplot2::geom_segment(data = temp, ggplot2::aes(x = x, xend = x, y = 0, yend = 0.5),
                                 lty = 2, lwd = median.lsize)

  temp2 <- data.frame(x = max(median_vals, na.rm = TRUE), y = 0.5)
  p + ggplot2::geom_segment(data = temp2, ggplot2::aes(x = 0, xend = x, y = 0.5, yend = 0.5),
                            lty = 2, lwd = median.lsize)
}

#' Add survival at set time lines to plot
#' @param p ggplot object
#' @param set.surv Data frame with survival at set times
#' @param set.lsize Line size
add_set_time_lines <- function(p, set.surv, set.lsize = 1) {

  if (is.null(set.surv)) return(p)

  set.surv$y <- as.numeric(sub(" *\\(.*", "", set.surv$set.CI))
  set.surv$time[is.na(set.surv$y)] <- NA

  p <- p + ggplot2::geom_segment(data = set.surv, ggplot2::aes(x = 0, xend = time, y = y, yend = y),
                                 lty = 2, lwd = set.lsize)
  p + ggplot2::geom_segment(data = set.surv, ggplot2::aes(x = time, xend = time, y = 0, yend = y),
                            lty = 2, lwd = set.lsize)
}

#' Apply colour and linetype scales
#' @param p ggplot object
#' @param col colours vector
#' @param linetype Line types vector
#' @param stratalabs Strata labels
#' @param eventlabs Event labels
#' @param multiple_lines Whether multiple strata
#' @param plot.event Events being plotted
#' @param event How events are distinguished
apply_scales_and_guides <- function(p, col, linetype = NULL, stratalabs, eventlabs = NULL,
                                    multiple_lines, plot.event = 1, event = "col") {

  # Determine labels for colour and linetype
  if (multiple_lines && (length(plot.event) == 1 || event == "linetype")) {
    col_labs <- stratalabs
  } else {
    col_labs <- eventlabs %||% paste("Event", plot.event)
  }

  if (multiple_lines && (length(plot.event) == 1 || event == "col")) {
    linetype_labs <- stratalabs
  } else {
    linetype_labs <- eventlabs %||% paste("Event", plot.event)
  }

  # Apply scales
  p <- p +
    ggplot2::scale_colour_manual(values = col, labels = col_labs) +
    ggplot2::scale_fill_manual(values = col, labels = col_labs)

  if (!is.null(linetype)) {
    p <- p + ggplot2::scale_linetype_manual(labels = linetype_labs, values = linetype)
  } else {
    p <- p + ggplot2::scale_linetype_discrete(labels = linetype_labs)
  }

  # Handle guide display for multi-event plots
  if (event == "linetype" && length(plot.event) == 2 && !multiple_lines) {
    p <- p + ggplot2::guides(col = "none")
  }

  if (event == "col" && length(plot.event) == 2 && !multiple_lines) {
    p <- p + ggplot2::guides(linetype = "none")
  }

  p
}

#' Create risk table for survival plot
#' @param sfit Survival fit object
#' @param times Time points for risk table
#' @param xlim X-axis limits
#' @param stratalabs Strata labels
#' @param stratalabs.table Table-specific strata labels
#' @param strataname.table Table strata name
#' @param Numbers_at_risk_text Text for numbers at risk
#' @param multiple_lines Whether multiple strata
#' @param col colours for strata
#' @param fsize Font size
#' @param nsize Number size in table
create_risk_table <- function(sfit, times, xlim,stratalabs, stratalabs.table = NULL,
                              strataname.table = "", Numbers_at_risk_text = "At Risk",
                              multiple_lines = TRUE, col = NULL, fsize = 12, nsize = 3) {

  # Determine maximum x value
  maxxval <- max(times)
  maxxlim <- if (is.null(xlim)) maxxval else xlim[2]

    if (is.null(stratalabs.table)) stratalabs.table <- stratalabs

  sfit.summary <- summary(sfit, times = times, extend = TRUE)

  risk.data <- data.frame(
    strata = if (multiple_lines) sfit.summary$strata else factor("All"),
    time = sfit.summary$time,
    n.risk = sfit.summary$n.risk
  ) |>
    transform(strata = factor(strata, levels = rev(levels(strata))))

  yticklabs <- unname(rev(stratalabs.table))

  if (is.null(col)) {
    cols1 <- scales::hue_pal()(length(stratalabs))
  } else {
    cols1 <- col
  }

  p <- ggplot2::ggplot(risk.data, ggplot2::aes(x = time, y = strata, label = format(n.risk, nsmall = 0))) +
    ggplot2::labs(x="",y="")+
    ggplot2::geom_text(hjust = "middle", vjust = "center", size = nsize) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(breaks = times) +
    ggplot2::coord_cartesian(xlim = c(0, maxxlim)) +
    ggplot2::scale_y_discrete(breaks = as.character(levels(risk.data$strata)), labels = yticklabs) +
    ggplot2::theme(
      legend.position = "none",
      text = ggplot2::element_text(size = fsize),
      plot.margin = ggplot2::unit(c(-0.5, 0.4, 0.1, 0.2), "lines"),
      axis.title.x = ggplot2::element_text(size = fsize, vjust = 1),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(face = "bold", hjust = 1, colour = rev(cols1))
    )
  if (!is.null(Numbers_at_risk_text)) {
    p <- p +
      ggtitle(Numbers_at_risk_text) +
      theme(plot.title = element_text(size = fsize,
                                      face = "bold",
                                      margin = margin(b = .5)))
  }
  # Handle large dash replacement if needed
  if (all(as.character(yticklabs) == "-")) {
    p <- .set_large_dash_as_ytext(p)
  }

  p
}

# Utility Functions ----
#' @noRd
process_legend_pos <- function(legend.pos) {
  if (is.character(legend.pos) && length(legend.pos) == 1) {
    # Use switch for string positions
    switch(legend.pos,
           "left" = "left",
           "top" = "top",
           "right" = "right",
           "bottom" = "bottom",
           "none" = "none",
           stop("Invalid legend position string. Must be one of: 'left', 'top', 'right', 'bottom', 'none'"))
  } else if (is.numeric(legend.pos) && length(legend.pos) == 2) {
    # Return numeric vector as-is for custom positioning
    legend.pos
  } else {
    stop("legend.pos must be either a character string or a numeric vector of length 2")
  }
}

# NOTE: break_function() is now defined in helper.R
# This eliminates code duplication across files
# Main Plotting Function ----

#' Main plotting function (refactored with separate CIF functions)
#' @param response Character vector with time and status column names
#' @param cov Covariate column name (optional)
#' @param data Input dataframe
#' @param type Plot type ("KM" or "CIF", auto-detected if NULL)
#' @param pval Whether to show p-values
#' @param conf.curves Whether to show confidence bands
#' @param table Whether to include risk table
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param col colours vector
#' @param times Numeric vector of times for the x-axis
#' @param plot.event Events to plot
#' @param returns Whether to return list with plot and at risk table
#' @param ... Additional arguments
#' @export
ggkmcif2 <- function(response, cov = NULL, data, pval = TRUE,
                     conf.curves = FALSE, table = TRUE, xlab = "Time",
                     ylab = NULL, col = NULL, times = NULL, type = NULL,
                     plot.event = 1, returns = FALSE, ...) {

  # Process arguments and set defaults
  mainArgs <- as.list(match.call(expand.dots = TRUE)[-1])
  toAdd <- mainArgs[setdiff(names(mainArgs), ls())]
  # Capture the parent environment before entering lapply & eval all
  parentEnv <- parent.frame()
  # toAdd <- lapply(toAdd, function(arg) {
  #   if (inherits(arg, "call")) arg <- eval(arg)
  #   return(arg)
  # })
  toAdd <- toAdd |>
    lapply(function(arg) eval(arg, envir = parentEnv))
  list2env(toAdd, environment())





  # Add default parameters
  if (!("strataname" %in% names(mainArgs))) {
    strataname <- ifelse(is.null(cov), "", nicename(cov))
  }
  if (!("legend.title" %in% names(mainArgs))) {
    legend.title <- ifelse(is.null(strataname), "", strataname)
  }

  defaultExtraArgs <- ggkmcif2Parameters(strataname = strataname)
  argsToAdd <- defaultExtraArgs[setdiff(names(defaultExtraArgs), names(mainArgs))]
  argsToAdd <- lapply(argsToAdd, function(arg) {
    if (inherits(arg, "call")) arg <- eval(arg)
    return(arg)
  })
  list2env(argsToAdd, environment())
  event <- event[1]

  # Auto-detect plot type
  if (is.null(type)) {
    if (length(unique(na.omit(data[[response[2]]]))) > 3) stop("More than 3 unique values detected in the event type. For KM curves there should be two unique values and for CIF either two, or three if there is a competing risk")
    type <- if (length(unique(na.omit(data[[response[2]]]))) < 3) "KM" else "CIF"
  }
  # Check & set legend position
  if (!("legend.pos" %in% names(mainArgs))){
    if (type == "KM") legend.pos <- c(0,0) else legend.pos <- c(0,1)
  } else legend.pos <- process_legend_pos(legend.pos)


  # Set default y-axis label
  if (type == "KM") {
    if (is.null(ylab)) ylab <- "Survival Probability"
  } else if (type == "CIF") {
    if (is.null(ylab)) ylab <- "Probability of an Event"
    # Disable median for CIF
    median.lines <- FALSE
    median.text <- FALSE
  } else {
    stop("Type must be either KM or CIF")
  }

  # Handle flip.CIF settings
  if (flip.CIF && type == "CIF") {
    median.text <- FALSE
    median.lines <- FALSE
    set.time.text <- NULL
    set.time.line <- FALSE
    set.time <- NULL
  }

  # Validate and prepare data
  data <- validate_and_prepare_data(data, response, cov, print.n.missing)

  # Process covariate
  cov_result <- process_covariate(data, cov, cut, stratalabs)
  data <- cov_result$data
  stratalabs <- cov_result$stratalabs %||% "All"
  multiple_lines <- !is.null(cov)

  # Set font size
  if (missing(fsize)) fsize <- get_current_base_size()
  if (missing(psize)) psize <- (fsize * 0.8) / .pt

  # Process colours and line types
  if (is.null(stratalabs.table)) stratalabs.table <- stratalabs

  if (!is.null(col) && !is.null(cov) && (event != "col" || length(plot.event) == 1)) {
    col <- rep(col, length.out = length(levels(data[[cov]])))
  }
  if (!is.null(linetype) && !is.null(cov) && (event != "linetype" || length(plot.event) == 1)) {
    linetype <- rep(linetype, length.out = length(levels(data[[cov]])))
  }

  if (is.null(col)) {
    col_length <- if (length(plot.event) == 1 || event != "col") {
      ifelse(is.null(cov), 1, length(levels(data[[cov]])))
    } else 2
    col <- color_palette_surv_ggplot(col_length)
  }

  # Initialize common variables
  median_vals <- NULL
  median_txt <- NULL
  set.surv <- NULL
  set.surv.text <- NULL
  sfit <- NULL

  # Fit models and create data frames based on type
  if (type == "KM") {
    # KM processing
    km_result <- fit_km_model(data, response, cov, conf.type)
    sfit <- km_result$sfit

    # Add hazard ratios to labels
    stratalabs <- add_km_hazard_ratios(stratalabs, data, response, cov, type, plot.event,
                                    HR, HR_pval, HR.digits, HR.pval.digits)

    # Add median times to labels
    median_result <- calculate_and_add_median_times(sfit, NULL, stratalabs, type, plot.event,
                                                    median.text, median.CI, median.digits)
    stratalabs <- median_result$stratalabs
    median_vals <- median_result$median_vals

    # Add time-specific estimates to labels
    time_result <- calculate_and_add_time_specific_estimates(sfit, NULL, stratalabs, type, plot.event,
                                                             set.time.text, set.time, set.time.line,
                                                             set.time.CI, set.time.digits)
    stratalabs <- time_result$stratalabs
    set.surv <- time_result$set.surv

    # Create plotting dataframe
    df <- create_km_dataframe(sfit, stratalabs, conf.curves, conf.type)

    # Store single-group median text for later use
    if (!multiple_lines && median.text) {
      median_txt <- paste0("Median=", round_sprintf(median_vals, median.digits))
    }

    # Store single-group set time text
    if (!multiple_lines && !is.null(set.time.text)) {
      set.surv.text <- paste(sapply(set.time, function(t) {
        idx <- which(set.surv$time == t)
        if (length(idx) > 0) {
          paste0(t, " ", set.time.text, "=", set.surv$set.CI[idx])
        }
      }), collapse = ",")
    }

  } else { # CIF

    # Fit CIF model
    cif_result <- fit_cif_model(data, response, cov)
    fit <- cif_result$fit
    gsep <- cif_result$gsep

    # Add hazard ratios for CIF
    stratalabs <- add_cif_hazard_ratios(stratalabs, data, response, cov,
                                        plot.event, HR, HR_pval, HR.digits, HR.pval.digits)

    # Process median values
    median_result <- process_cif_medians(fit, plot.event, stratalabs,
                                         median.lines, median.text, median.digits,
                                         multiple_lines)
    stratalabs <- median_result$stratalabs
    median_vals <- median_result$median_vals
    median_txt <- median_result$median_txt

    # Process time-specific estimates
    timepoint_result <- process_cif_timepoints(fit, plot.event, stratalabs,
                                               set.time.text, set.time, set.time.line,
                                               set.time.CI, set.time.digits, multiple_lines)
    stratalabs <- timepoint_result$stratalabs
    set.surv <- timepoint_result$set.surv
    set.surv.text <- timepoint_result$set.surv.text

    # Create plotting dataframe
    df <- create_cif_dataframe(fit, gsep, plot.event, stratalabs,
                               conf.type, flip.CIF, eventlabs, cov, data)

    # Create survival fit for risk table if needed
    if (table) {
      sfit <- create_cif_risk_table_sfit(data, response, cov)
    }
  }


  # Set time breaks
  m <- max(nchar(stratalabs))
  maxxval <- max(df$time, if (!is.null(times)) times[length(times)] else 0)
  maxxlim <- if (is.null(xlim)) maxxval else xlim[2]
  if (is.null(times)) times <- break_function(maxxlim)
  if (is.null(xlim)) xlim <- c(min(times), max(times))


  # Create base plot
  p <- create_base_plot(df, type, xlab, ylab, multiple_lines, plot.event, event,
                        lsize, fsize, col, linetype, legend.pos, legend.title,
                        times, ylim, xlim, main)

  # Add confidence bands
  if (conf.curves) {
    p <- add_confidence_bands(p, df, type, event, plot.event)
  }

  # Add censor marks for KM
  if (censor.marks && type == "KM") {
    p <- add_censor_marks(p, df, censor.size, censor.stroke, censor.symbol)
  }

  # Add statistical tests
  if (pval && multiple_lines) {
    if (type == "KM") {
      p <- add_statistical_tests(p, type, multiple_lines, km_result$pval, pval.pos,
                                 times, xlim, ylim, psize, pval.digits, plot.event, eventlabs)
    } else if (type == "CIF") {
      gray_pval <- extract_grays_test(df, plot.event)
      p <- add_statistical_tests(p, type, multiple_lines, gray_pval, pval.pos,
                                 times, xlim, ylim, psize, pval.digits, plot.event, eventlabs)
    }
  }

  # Add median text for single group
  if (median.text && !multiple_lines && !is.null(median_txt)) {
    p <- add_median_text(p, type, multiple_lines, median_txt, median.pos,
                         times, ylim, median.size, plot.event, eventlabs)
  }

  # Add set time text for single group
  if (!is.null(set.time.text) && !multiple_lines && !is.null(set.surv.text)) {
    p <- add_set_time_text(p, type, multiple_lines, set.surv.text, set.pos,
                           times, ylim, set.size, plot.event, eventlabs)
  }

  # Add median lines
  if (median.lines && !is.null(median_vals)) {
    p <- add_median_lines(p, median_vals, median.lsize)
  }

  # Add set time lines
  if (set.time.line && !is.null(set.surv)) {
    p <- add_set_time_lines(p, set.surv, set.lsize)
  }

  # Apply colour and linetype scales
  p <- apply_scales_and_guides(p, col, linetype, stratalabs, eventlabs,
                               multiple_lines, plot.event, event)

  # Create combined plot with risk table
  if (table && !is.null(sfit)) {
    data.table <- create_risk_table(sfit, times, xlim,
                                    stratalabs, stratalabs.table,
                                    strataname.table, Numbers_at_risk_text,
                                    multiple_lines, col, fsize, nsize)

    # Set table height
    line_size_in_inches <- fsize/72
    dev_height_inches <- dev.size("in")[2]
    p_risk_lines <- length(unique(stratalabs))+2
    p2_risk_height <- line_size_in_inches*p_risk_lines
    p1_height <-  dev_height_inches - p2_risk_height
    rel_height = c(1,(dev_height_inches-p1_height)/p1_height)


    p_combined <- cowplot::plot_grid(p, data.table, nrow = 2, ncol = 1,
                                     rel_heights = rel_height,
                                     align = 'v', axis = 'lr')

    if (returns) {
      return(list(plot = p, table = data.table, combined = p_combined))
    } else {
      return(p_combined)
    }
  } else {
    if (returns) {
      return(list(plot = p))
    } else {
      return(p)
    }
  }
}

#' Additional parameters passed to ggkmcif2
#'
#' This section documents the additional parameters for \link{ggkmcif2}.
#'
#' @param HR boolean to specify if you want hazard ratios included in the plot
#' @param HR_pval boolean to specify if you want HR p-values in the plot
#' @param conf.type One of "none"(the default), "plain", "log" , "log-log" or
#'   "logit". Only enough of the string to uniquely identify it is necessary.
#'   The first option causes confidence intervals not to be generated. The
#'   second causes the standard intervals curve +- k *se(curve), where k is
#'   determined from conf.int. The log option calculates intervals based on the
#'   cumulative hazard or log(survival). The log-log option bases the intervals
#'   on the log hazard or log(-log(survival)), and the logit option on
#'   log(survival/(1-survival)).
#' @param main String corresponding to main title. When NULL uses Kaplan-Meier
#'   Plot s, and "Cumulative Incidence Plot for CIF"
#'
#' @param stratalabs string corresponding to the labels of the covariate, when
#'   NULL will use the levels of the covariate
#' @param strataname String of the covariate name default is  nicename(cov)
#' @param stratalabs.table String corresponding to the levels of the covariate
#'   for the number at risk table, when NULL will use the levels of the
#'   covariate. Can use a string of "-" when the labels are long
#' @param strataname.table String of the covariate name for the number at risk
#'   table default is  nicename(cov
#'
#' @param median.text boolean to specify if you want the median values added to
#'   the legend (or as added text if there are no covariates), for KM only
#' @param median.lines boolean to specify if you want the median values added as
#'   lines to the plot, for KM only
#' @param median.CI boolean to specify if you want the 95\% confidence interval
#'   with the median text (Only for KM)
#' @param set.time.text string for the text to add survival at a specified time
#'   (eg. year OS)
#' @param set.time.line boolean to specify if you want the survival added as
#'   lines to the plot at a specified point
#' @param set.time Numeric values of the specific time of interest, default is 5
#'   (Multiple values can be entered)
#' @param set.time.CI boolean to specify if you want the 95\% confidence
#'   interval with the set time text
#'
#' @param censor.marks logical value. If TRUE, includes censor marks (only for
#'   KM curves)
#' @param censor.size size of censor marks, default is 3
#' @param censor.stroke stroke of censor marks, default is 1.5
#' @param censor.symbol either a character or a number 0-24 specifying the
#'   ggplot shape to be used as the censor symbol
#' @param fsize font size
#' @param nsize font size for numbers in the numbers at risk table
#' @param lsize line size
#' @param psize size of the pvalue
#' @param median.size size of the median text (Only when there are no
#'   covariates)
#' @param median.pos vector of length 2 corresponding to the median position
#'   (Only when there are no covariates)
#' @param median.lsize line size of the median lines
#' @param set.size size of the survival at a set time text (Only when there are
#'   no covariates)
#' @param set.pos  vector of length 2 corresponding to the survival at a set
#'   point position (Only when there are no covariates)
#' @param set.lsize line size of the survival at set points
#' @param ylim vector of length 2 corresponding to limits of y-axis. Default to
#'   NULL
#' @param linetype vector of line types; default is solid for all lines
#' @param xlim  vector of length 2 corresponding to limits of x-axis. Default to
#'   NULL
#' @param legend.pos A string corresponding to the legend position
#'   ("left","top", "right", "bottom", "none") or a numeric vector specifying
#'   the internal coordinates of the plot ie c(0.5,.0.5) for the centre of the
#'   plot.
#' @param legend.title a string for the title of the legend, defaults to
#'   strataname
#' @param pval.pos  vector of length 2 corresponding to the p-value position
#' @param event String specifying if the event should be mapped to the colour,
#'   or linetype when plotting both events to colour = "col", line type
#' @param flip.CIF boolean to flip the CIF curve to start at 1
#' @param cut numeric value indicating where to divide a continuous covariate
#'   (default is the median)
#' @param eventlabs String corresponding to the event type names
#' @param event.name String corresponding to the label of the event types
#' @param Numbers_at_risk_text String for the label of the number at risk, set
#'   Numbers_at_risk_text=NULL to remove
#' @param tbl.height Height of the at risk table, relative to plot. To set the
#'   table to half the height of the plot use tbl.height = 0.5
#' @param HR.digits Number of digits printed of the  hazard ratio
#' @param HR.pval.digits Number of digits printed of the hazard ratio pvalue
#' @param pval.digits Number of digits printed of the Gray's/log rank pvalue
#'
#' @param median.digits Number of digits printed of the median pvalue
#' @param set.time.digits Number of digits printed of the probability at a
#'   specified time
#' @param print.n.missing Logical, should the number of missing be shown !Needs
#'   to be checked
#' @param returns Logical, if TRUE a list contain the plot and at risk table is
#'   returned
#'
#' @name ggkmcif2Parameters
#' @export
ggkmcif2Parameters <- function(HR = FALSE, HR_pval = FALSE, conf.type = "log",
                               main = NULL, stratalabs = NULL, strataname,
                               stratalabs.table = NULL, strataname.table = strataname,
                               median.text = FALSE, median.lines = FALSE, median.CI = FALSE,
                               set.time.text = NULL, set.time.line = FALSE, set.time = 5,
                               set.time.CI = FALSE,
                               censor.marks = TRUE, censor.size = 2,
                               censor.stroke = 1.5, censor.symbol ="|",
                               fsize, nsize = 3, lsize = .7, psize = 3.5,
                               median.size = 3, median.pos = NULL, median.lsize = 1, set.size = 3,
                               set.pos = NULL, set.lsize = 1, ylim = c(0, 1),
                               linetype = NULL, xlim = NULL, legend.pos,
                               legend.title = strataname,
                               pval.pos = NULL,
                               event = c("col", "linetype"), flip.CIF = FALSE,
                               cut = NULL, eventlabs = NULL, event.name = NULL,
                               Numbers_at_risk_text = "At risk",
                               tbl.height = NULL,
                               HR.digits = 2, HR.pval.digits = 3, pval.digits = 3, median.digits = 3,
                               set.time.digits = 3, print.n.missing = TRUE, returns = FALSE) {
  return(as.list(environment(), all = TRUE))
}

# Additional Required Functions ----

#' Function to extract ggplot colours
#' @noRd
.extract_ggplot_colours <- function(p, grp.levels) {
  # Extract colours from ggplot object
  g <- ggplot_build(p)
  colours <- unique(g$data[[1]]$colour)
  if (length(colours) < length(grp.levels)) {
    colours <- rep(colours, length.out = length(grp.levels))
  }
  colours[1:length(grp.levels)]
}

#' Function to set large dash as y-axis text
#' @noRd
.set_large_dash_as_ytext <- function(plot) {
  # Handle large dash formatting for y-axis text
  plot + ggplot2::theme(axis.text.y = ggplot2::element_text(family = "mono"))
}

calculate_text_height_base <- function(text = "Ag", font_size = 11, units = "inches") {
  height_inches <- strheight(text, cex = font_size / 12, units = "inches")

  switch(units,
         "cm" = height_inches * 2.54,
         "mm" = height_inches * 25.4,
         "inches" = height_inches,
         stop("Supported units: 'cm', 'mm', 'inches'"))
}
