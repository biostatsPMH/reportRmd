forestplotMV_updated <- function(model, data, conf.level = 0.95, orderByRisk = TRUE,
                         colours = "default", showEst = TRUE, rmRef = FALSE,
                         digits = getOption("reportRmd.digits", 2),
                         logScale = getOption("reportRmd.logScale", TRUE),
                         nxTicks = 5, showN = TRUE, showEvent = TRUE) {

  if (inherits(model, "glm")) {
    if (model$family$link == "log") {
      x_lab <- "Adjusted Relative Risk"
    }
    else if (model$family$link == "logit") {
      x_lab <- "Adjusted Odds Ratio"
    }
    else stop("model must be a logit or log link fit")
  }
  else {
    x_lab <- "Adjusted Odds Ratio"
  }

  ###################################
  # Get data from m_summary
  tab <- m_summary(model, CIwidth = conf.level, digits = digits, for_plot = TRUE)

  # Handle ordering (replicate original logic exactly)
  if (orderByRisk) {
    # Original: tab <- tab[order(rank(tab$estimate), tab$var.name), ]
    tab$estimate_for_order <- ifelse(is.na(tab$est), 1, tab$est)
    tab <- tab[order(rank(tab$estimate_for_order), tab$var), ]

    # Original variable ordering logic
    dt <- as.data.frame(unique(tab$var))
    colnames(dt) <- "var"
    dt$var.order <- 1:nrow(dt)
    dt$var.order <- dt$var.order + 1
    tab <- merge(tab, dt, by = "var", all = TRUE)

    # Create level.order (original: sequence(rle(tab$var.name)$lengths))
    tab$level.order <- sequence(rle(tab$var)$lengths)

    # Final ordering: original was tab$var.order, -tab$level.order
    tab <- tab[order(tab$var.order, -tab$level.order), ]
  } else {
    # Use original ordering from m_summary
    tab <- tab[order(tab$ord), ]
  }

  # Remove reference categories if requested
  if (rmRef) {
    tab <- tab[!(tab$Est_CI == "Reference"), ]
  }

  # Check if Events column exists
  if (!("Events" %in% names(tab))) {
    showEvent <- FALSE
  }

  ###################################
  # Set up plotting data with correct y-axis ordering
  yvals <- nrow(tab):1  # Reverse order so first row appears at top

  # Create y-axis labels
  if (showEst) {
    yLabels <- data.frame(
      y.pos = yvals,
      labels = ifelse(is.na(tab$header),
                      ifelse(is.na(tab$lvl),
                             paste(tab$var, ": ", tab$Est_CI, sep = ""),
                             paste(tab$lvl, ": ", tab$Est_CI, sep = "")),
                      tab$var)
    )
  } else {
    yLabels <- data.frame(
      y.pos = yvals,
      labels = ifelse(is.na(tab$header),
                      ifelse(is.na(tab$lvl), tab$var, tab$lvl),
                      tab$var)
    )
  }

  yLabels$labels <- gsub("_", " ", yLabels$labels)
  yLabels <- yLabels[!is.na(yLabels$y.pos), ]

  # Set up plot coordinates and colors
  tab$x.val <- ifelse(tab$Est_CI == "Reference", 1, tab$est)
  tab$y.val <- yvals
  tab$colour <- ifelse(tab$x.val < 1, "a", ifelse(tab$x.val == 1, "b", "c"))

  # Set up colors
  if (length(colours) == 1) {
    colours <- c(a = "#006B3C", b = "black", c = "#FF0800")
  } else {
    names(colours) <- c("a", "b", "c")
  }

  # Set up axes based on showN and showEvent parameters
  if (showN & !showEvent) {
    Axis <- scale_y_continuous(
      breaks = yLabels$y.pos,
      labels = yLabels$labels,
      sec.axis = dup_axis(breaks = yLabels$y.pos,
                          labels = tab$n, name = "N")
    )
    themeSecAxis <- theme(axis.title.y.right = element_text(angle = 0, hjust = 0.5, vjust = 0.5))
  }

  if (showEvent && "Events" %in% names(tab)) {
    Axis <- scale_y_continuous(
      breaks = yLabels$y.pos,
      labels = yLabels$labels,
      sec.axis = dup_axis(breaks = yLabels$y.pos,
                          labels = paste(tab$n, tab$Events, sep = " : "),
                          name = "N : Event")
    )
    themeSecAxis <- theme(axis.title.y.right = element_text(angle = 270, hjust = 0.5, vjust = 0.5))
  }

  if (!showN & !showEvent) {
    Axis <- scale_y_continuous(breaks = yLabels$y.pos, labels = yLabels$labels)
    themeSecAxis <- NULL
  }

  colours <- colours[sort(unique(tab$colour))]

  # Create the plot
  suppressWarnings({
    tryCatch({
      p <- ggplot(tab, aes_(x = ~x.val, y = ~y.val, colour = ~colour)) +
        geom_point(na.rm = TRUE, size = 2) +
        geom_errorbarh(aes_(xmin = ~lwr, xmax = ~upr),
                       height = 0, size = 0.9, na.rm = TRUE) +
        geom_vline(xintercept = 1) +
        labs(y = "", x = x_lab) +
        guides(colour = "none") +
        Axis +
        scale_colour_manual(values = colours) +
        theme_bw() +
        theme(axis.text.y = element_text(
          face = ifelse(!is.na(tab$header), "bold", "plain"),
          hjust = 0),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.ticks = element_blank()) +
        themeSecAxis

      if (logScale) {
        p + scale_x_log10(breaks = scales::log_breaks(n = nxTicks))
      } else {
        p
      }
    }, error = function(e) {})
  })
}
