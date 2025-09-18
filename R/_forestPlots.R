# Helper function to order data by risk while keeping variable levels together
.order_by_risk <- function(tab) {
  # For each variable, find the maximum absolute log odds ratio
  # This handles both OR > 1 and OR < 1 scenarios
  var_max_risk <- aggregate(estimate, by = list(var.name = tab$var.name),
                            FUN = function(x) max(abs(log(x)), na.rm = TRUE))
  names(var_max_risk)[2] <- "max_log_or"

  # Merge back with original data
  tab <- merge(tab, var_max_risk, by = "var.name", all.x = TRUE)

  # Order by maximum risk within variable, then by variable name, then by level order
  tab <- tab[order(-tab$max_log_or, tab$var.name, -tab$level.order), ]

  # Remove the temporary column
  tab$max_log_or <- NULL

  return(tab)
}
.safe_numeric <- function(x) {
  # Handle empty strings and non-numeric values
  x[x == "" | is.na(x)] <- NA
  suppressWarnings(as.numeric(x))
}

#' Shared forest plot creation function
#'
#' Creates a forest plot from standardized data. This is an internal function used by
#' both forestplotUV and forestplotMV to ensure consistent plotting behavior.
#'
#' @param tab A data.frame containing the standardized forest plot data. Must contain
#'   the following columns:
#'   \describe{
#'     \item{variable}{Character. The variable name (e.g., "age", "sex")}
#'     \item{var.name}{Character. The forward-filled variable name for grouping}
#'     \item{level.name}{Character. The level name for categorical variables (NA for continuous)}
#'     \item{estimate}{Numeric. The point estimate (OR, RR, HR, etc.)}
#'     \item{estimate.label}{Character. Formatted estimate with CI (e.g., "1.23 (0.98-1.45)" or "Reference")}
#'     \item{conf.low}{Numeric. Lower confidence interval bound}
#'     \item{conf.high}{Numeric. Upper confidence interval bound}
#'     \item{N}{Numeric. Sample size for each group (optional, used if showN=TRUE)}
#'     \item{Event}{Numeric. Number of events for each group (optional, used if showEvent=TRUE)}
#'   }
#'   Additional columns may be present but are not used by this function.
#' @param x_lab Character. Label for the x-axis (e.g., "Odds Ratio", "Hazard Ratio")
#' @param conf.level Numeric. Confidence level for intervals (default 0.95)
#' @param orderByRisk Logical. Whether to order variables by risk (currently not implemented)
#' @param colours Character vector or "default". Colors for estimates <1, =1, >1
#' @param showEst Logical. Whether to show estimates in y-axis labels
#' @param rmRef Logical. Whether to remove reference category rows
#' @param logScale Logical. Whether to use log scale for x-axis
#' @param nxTicks Numeric. Number of tick marks for log scale
#' @param showN Logical. Whether to show sample sizes on secondary y-axis
#' @param showEvent Logical. Whether to show event counts on secondary y-axis
#'
#' @return A ggplot2 object representing the forest plot
#'
#' @details
#' The tab data.frame should be prepared using either .prepare_uv_data() or
#' .prepare_mv_data() to ensure proper formatting. The estimate.label column
#' should contain "1.0 (Reference)" for reference categories and formatted
#' estimates like "1.23 (0.98-1.45)" for other categories.
#'
#' @keywords internal
.create_forest_plot <- function(tab, x_lab, conf.level = 0.95, orderByRisk = TRUE,
                                colours = "default", showEst = TRUE, rmRef = FALSE,
                                logScale = getOption("reportRmd.logScale", TRUE),
                                nxTicks = 5, showN = TRUE, showEvent = TRUE) {

  # Remove reference rows if requested
  if (rmRef) {
    tab <- tab[setdiff(1:nrow(tab), which(tab$estimate.label == "1.0 (Reference)")), ]
  }

  yvals <- 1:nrow(tab)
  tab$estimate.label <- ifelse(is.na(tab$estimate.label), "", tab$estimate.label)
  tab$estimate.label <- ifelse(tab$estimate.label == "1.0 (Reference)", "(Reference)", tab$estimate.label)

  # Create y-axis labels
  if (showEst) {
    yLabels <- data.frame(
      y.pos = yvals,
      labels = ifelse(is.na(tab$level.name),
                      paste(tab$variable, ": ", tab$estimate.label, sep = ""),
                      paste(tab$level.name, ": ", tab$estimate.label, sep = ""))
    )
  } else {
    yLabels <- data.frame(
      y.pos = yvals,
      labels = ifelse(is.na(tab$level.name),
                      tab$variable,
                      ifelse(tab$estimate.label == "(Reference)",
                             paste(tab$level.name, ": ", tab$estimate.label, sep = ""),
                             tab$level.name))
    )
  }

  yLabels$labels <- gsub("_", " ", yLabels$labels)
  yLabels <- yLabels[!is.na(yLabels$y.pos), ]

  # Prepare plot data
  tab$x.val <- ifelse(tab$estimate.label == "(Reference)", 1, tab$estimate)
  tab$y.val <- yLabels$y.pos
  tab$colour <- ifelse(tab$x.val < 1, "a", ifelse(tab$x.val == 1, "b", "c"))

  # Set colors
  if (length(colours) == 1) {
    colours <- c(a = "#006B3C", b = "black", c = "#FF0800")
  } else {
    names(colours) <- c("a", "b", "c")
  }
  colours <- colours[sort(unique(tab$colour))]

  # Create secondary axis
  Axis <- scale_y_continuous(breaks = yLabels$y.pos, labels = yLabels$labels)
  themeSecAxis <- NULL

  if (showN & !showEvent) {
    Axis <- scale_y_continuous(
      breaks = yLabels$y.pos,
      labels = yLabels$labels,
      sec.axis = dup_axis(breaks = yLabels$y.pos, labels = tab$N, name = "N")
    )
    themeSecAxis <- theme(axis.title.y.right = element_text(angle = 0, hjust = 0.5, vjust = 0.5))
  }

  if (showEvent & "Event" %in% names(tab)) {
    Axis <- scale_y_continuous(
      breaks = yLabels$y.pos,
      labels = yLabels$labels,
      sec.axis = dup_axis(breaks = yLabels$y.pos,
                          labels = paste(tab$N, tab$Event, sep = " : "),
                          name = "N : Event")
    )
    themeSecAxis <- theme(axis.title.y.right = element_text(angle = 270, hjust = 0.5, vjust = 0.5))
  }

  # Create the plot
  suppressWarnings({
    tryCatch({
      p <- ggplot(tab, aes_(x = ~x.val, y = ~y.val, colour = ~colour)) +
        geom_point(na.rm = TRUE, size = 2) +
        geom_errorbarh(aes_(xmin = ~conf.low, xmax = ~conf.high),
                       height = 0, size = 0.9, na.rm = TRUE) +
        geom_vline(xintercept = 1) +
        labs(y = "", x = x_lab) +
        guides(colour = "none") +
        Axis +
        scale_colour_manual(values = colours) +
        theme_bw() +
        theme(
          axis.text.y = element_text(
            face = ifelse(tab$variable == tab$var.name | is.na(tab$var.name), "bold", "plain"),
            hjust = 0
          ),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.ticks = element_blank()
        ) +
        themeSecAxis

      if (logScale) {
        p + scale_x_log10(breaks = scales::log_breaks(n = nxTicks))
      } else {
        p
      }
    }, error = function(e) { NULL })
  })
}

# Helper function to prepare UV data
.prepare_uv_data <- function(response, covs, data, digits, id, corstr, family,
                             conf.level, showN, showEvent, orderByRisk = TRUE) {
  tab <- uvsum2(response, covs, data, digits = digits, id = id, corstr = corstr,
                family = family, type = NULL, gee = FALSE, strata = 1,
                nicenames = FALSE, showN = showN, showEvent = showEvent,
                CIwidth = conf.level, reflevel = NULL, returnModels = FALSE)

  # Standardize column names and format
  tab$estimate.label <- tab[, 2]
  tab$estimate.label[which(tab$estimate.label == "Reference")] <- "1.0 (Reference)"

  # Extract numeric values more safely
  estimate_text <- gsub(" .*", "", tab[, 2])
  estimate_text[estimate_text == "Reference"] <- "1"
  tab$estimate <- .safe_numeric(estimate_text)

  # Extract confidence intervals more safely
  ci_text <- gsub("\\(([^()]*)\\)|.", "\\1", tab[, 2])
  conf_low_text <- gsub(",.*", "", ci_text)
  conf_high_text <- gsub("^\\S*\\s+", "", ci_text)

  tab$conf.low <- .safe_numeric(conf_low_text)
  tab$conf.high <- .safe_numeric(conf_high_text)

  tab$level.name <- tab[, 1]
  tab$var.name <- NA
  tab$var.name[which(tab$Covariate %in% covs)] <- tab$level.name[which(tab$Covariate %in% covs)]

  # Forward fill variable names
  y <- tab$var.name
  y_forward_fill <- fillNAs(y)
  tab <- cbind(tab, y, y_forward_fill)
  tab$var.name <- tab$y_forward_fill
  tab$level.order <- sequence(rle(tab$var.name)$lengths)

  # Order by risk if requested (same as MV function)
  if (orderByRisk) {
    tab <- tab[order(rank(tab$estimate), tab$var.name), ]
  }

  # Add variable order
  dt <- as.data.frame(unique(tab$var.name))
  colnames(dt) <- "var.name"
  dt$var.order <- 1:nrow(dt) + 1
  tab <- merge(tab, dt, by = "var.name", all = TRUE)
  tab <- tab[order(tab$var.order, -tab$level.order), ]

  # Standardize remaining columns
  tab$p.value <- tab$"p-value"
  tab$p.label <- paste(format(round(as.numeric(tab$p.value), 3), nsmall = 3), sep = "")
  tab$variable <- tab$Variable

  # Select final columns
  vars <- c("variable", "var.name", "level.name", "level.order", "estimate",
            "p.label", "p.value", "conf.low", "conf.high", "var.order",
            "estimate.label", "N", "Event")
  if (!("Event" %in% names(tab))) {
    vars <- setdiff(vars, "Event")
  }
  tab <- tab[, vars]
  as.data.frame(tab)
}

# Helper function to prepare MV data
.prepare_mv_data <- function(model, data, digits, conf.level, showN, showEvent, orderByRisk = TRUE) {
  tab <- mvsum(model, data, digits = digits, markup = FALSE, sanitize = FALSE,
               nicenames = FALSE, showN = showN, showEvent = showEvent,
               CIwidth = conf.level)

  # Standardize column names and format
  tab$estimate.label <- tab[, 2]
  tab$estimate.label[which(tab$estimate.label == "Reference")] <- "1.0 (Reference)"

  # Extract numeric values more safely
  estimate_text <- gsub(" .*", "", tab[, 2])
  estimate_text[estimate_text == "Reference"] <- "1"
  tab$estimate <- .safe_numeric(estimate_text)

  # Extract confidence intervals more safely
  ci_text <- gsub("\\(([^()]*)\\)|.", "\\1", tab[, 2])
  conf_low_text <- gsub(",.*", "", ci_text)
  conf_high_text <- gsub("^\\S*\\s+", "", ci_text)

  tab$conf.low <- .safe_numeric(conf_low_text)
  tab$conf.high <- .safe_numeric(conf_high_text)

  tab$level.name <- tab[, 1]
  tab$var.name <- NA

  # Get covariates from model
  covs <- colnames(model$model)
  tab$var.name[which(tab$Covariate %in% covs)] <- tab$level.name[which(tab$Covariate %in% covs)]

  # Forward fill variable names
  y <- tab$var.name
  y_forward_fill <- fillNAs(y)
  tab <- cbind(tab, y, y_forward_fill)
  tab$var.name <- tab$y_forward_fill
  tab$level.order <- sequence(rle(tab$var.name)$lengths)

  # Order by risk if requested (keeping variable levels together)
  if (orderByRisk) {
    tab <- .order_by_risk(tab)
  }

  # Add variable order AFTER risk ordering
  dt <- as.data.frame(unique(tab$var.name))
  colnames(dt) <- "var.name"
  dt$var.order <- 1:nrow(dt)
  tab <- merge(tab, dt, by = "var.name", all = TRUE)

  # Final sort to ensure proper display order
  if (orderByRisk) {
    # If ordered by risk, maintain that order
    tab <- tab[order(tab$var.order, tab$level.order), ]
  } else {
    # Default alphabetical ordering
    tab <- tab[order(tab$var.order, -tab$level.order), ]
  }

  # Standardize remaining columns
  tab$p.value <- tab$"p-value"
  tab$p.label <- paste(format(round(as.numeric(tab$p.value), 3), nsmall = 3), sep = "")
  tab$variable <- tab$Covariate

  # Select final columns
  vars <- c("variable", "var.name", "level.name", "level.order", "estimate",
            "p.label", "p.value", "conf.low", "conf.high", "var.order",
            "estimate.label", "N", "Event")
  if (!("Event" %in% names(tab))) {
    vars <- setdiff(vars, "Event")
  }
  tab <- tab[, vars]
  as.data.frame(tab)
}

#' Create an univariable forest plot using ggplot2
#'
#' This function will send and take log or logistic regression fit from glm or geeglm
#' from uvsum function, and display the OR or RR for each variable on the appropriate log scale.
#'
#' @param response character vector with names of columns to use for response
#' @param covs character vector with names of columns to use for covariates
#' @param data dataframe containing your data
#' @param model fitted model object
#' @param id character vector which identifies clusters. Only used for geeglm
#' @param corstr character string specifying the correlation structure. Only
#'   used for geeglm. The following are permitted: '"independence"',
#'   '"exchangeable"', '"ar1"', '"unstructured"' and '"userdefined"'
#' @param family description of the error distribution and link function to be
#'   used in the model. Only used for geeglm
#' @param digits number of digits to round to
#' @param conf.level controls the width of the confidence interval
#' @param orderByRisk logical, should the plot be ordered by risk
#' @param colours can specify colours for risks less than, 1 and greater than
#'   1.0. Default is red, black, green
#' @param showEst logical, should the risks be displayed on the plot in text
#' @param rmRef logical, should the reference levels be removed for the plot?
#' @param logScale logical, should OR/RR be shown on log scale, defaults to
#'   TRUE, or reportRmd.logScale if set. See https://doi.org/10.1093/aje/kwr156 for why you may prefer a
#'   linear scale.
#' @param nxTicks Number of tick marks supplied to the log_breaks function to
#'   produce
#' @param showN Show number of observations per variable and category
#' @param showEvent Show number of events per variable and category
#' @import ggplot2
#' @importFrom scales log_breaks
#' @keywords plot
#' @return a plot object
#' @export
#' @examples
#' data("pembrolizumab")
#' forestplotUV(response="orr", covs=c("change_ctdna_group", "sex", "age", "l_size"),
#' data=pembrolizumab, family='binomial')
forestplotUV2 <- function(response, covs, data, id = NULL, corstr = NULL,
                         model = "glm", family = NULL, digits = getOption("reportRmd.digits", 2),
                         conf.level = 0.95, orderByRisk = TRUE, colours = "default",
                         showEst = TRUE, rmRef = FALSE,
                         logScale = getOption("reportRmd.logScale", TRUE),
                         nxTicks = 5, showN = TRUE, showEvent = TRUE) {

  # Prepare data
  tab <- .prepare_uv_data(response, covs, data, digits, id, corstr, family,
                          conf.level, showN, showEvent, orderByRisk)

  # Determine x-axis label
  x_lab <- "Unadjusted Odds Ratio"  # Default
  if (inherits(model, "glm")) {
    if (model$family$link == "log") {
      x_lab <- "Unadjusted Relative Risk"
    } else if (model$family$link == "logit") {
      x_lab <- "Unadjusted Odds Ratio"
    }
  }

  # Create and return the plot
  .create_forest_plot(tab, x_lab, conf.level, orderByRisk, colours, showEst,
                      rmRef, logScale, nxTicks, showN, showEvent)
}

#' Create a multivariable forest plot using ggplot2
#'
#' This function will send and take log or logistic regression fit from glm or geeglm
#' from mvsum function, and display the OR or RR for each variable on the appropriate log scale.
#'
#' @param model an object output from the glm or geeglm function, must be from a logistic
#'   regression
#' @param data dataframe containing your data
#' @param conf.level controls the width of the confidence interval
#' @param orderByRisk logical, should the plot be ordered by risk
#' @param colours can specify colours for risks less than, 1 and greater than
#'   1.0. Default is red, black, green
#' @param showEst logical, should the risks be displayed on the plot in text
#' @param rmRef logical, should the reference levels be removed for the plot?
#' @param digits number of digits to use displaying estimates
#' @param logScale logical, should OR/RR be shown on log scale, defaults to
#'   TRUE, or reportRmd.logScale if set. See https://doi.org/10.1093/aje/kwr156 for why you may prefer a
#'   linear scale.
#' @param nxTicks Number of tick marks supplied to the log_breaks function to
#'   produce
#' @param showN Show number of observations per variable and category
#' @param showEvent Show number of events per variable and category
#' @import ggplot2
#' @importFrom scales log_breaks
#' @keywords plot
#' @return a plot object
#' @export
#' @examples
#' data("pembrolizumab")
#' glm_fit = glm(orr~change_ctdna_group+sex+age+l_size,
#' data=pembrolizumab,family = 'binomial')
#' forestplotMV(glm_fit)
forestplotMV2 <- function(model, data, conf.level = 0.95, orderByRisk = TRUE,
                         colours = "default", showEst = TRUE, rmRef = FALSE,
                         digits = getOption("reportRmd.digits", 2),
                         logScale = getOption("reportRmd.logScale", TRUE),
                         nxTicks = 5, showN = TRUE, showEvent = TRUE) {

  # Prepare data
  tab <- .prepare_mv_data(model, data, digits, conf.level, showN, showEvent, orderByRisk)

  # Determine x-axis label
  x_lab <- "Adjusted Odds Ratio"  # Default
  if (inherits(model, "glm")) {
    if (model$family$link == "log") {
      x_lab <- "Adjusted Relative Risk"
    } else if (model$family$link == "logit") {
      x_lab <- "Adjusted Odds Ratio"
    }
  }

  # Create and return the plot
  .create_forest_plot(tab, x_lab, conf.level, orderByRisk, colours, showEst,
                      rmRef, logScale, nxTicks, showN, showEvent)
}
