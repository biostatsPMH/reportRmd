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
#' @noRd
#' @examples
#' data("pembrolizumab")
#' forestplotUV(response="orr", covs=c("change_ctdna_group", "sex", "age", "l_size"),
#' data=pembrolizumab, family='binomial')
forestplotUV = function (response, covs, data, id = NULL, corstr = NULL,
                         model = "glm", family = NULL, digits = getOption("reportRmd.digits",2), conf.level = 0.95,
                         orderByRisk = TRUE, colours = "default", showEst = TRUE, rmRef = FALSE,
                         logScale=getOption("reportRmd.logScale",TRUE), nxTicks = 5, showN = TRUE, showEvent = TRUE)
{
  if (inherits(model, "glm")) {
    if (model$family$link == "log") {
      x_lab = "Unadjusted Relative Risk"
    }
    else if (model$family$link == "logit") {
      x_lab = "Unadjusted Odds Ratio"
    }
    else stop("model must be a logit or log link fit")
  }
  else {
    x_lab = "Unadjusted Odds Ratio"
  }
  tab = uvsum2(response, covs, data, digits = digits, id = NULL, corstr = NULL,
               family = NULL, type = NULL, gee = FALSE, strata = 1,
               nicenames = F, showN = TRUE, showEvent = TRUE,
               CIwidth = conf.level, reflevel = NULL, returnModels = FALSE)
  tab$estimate.label <- tab[,2];
  tab$estimate.label[which(tab$estimate.label == "Reference")] <- "1.0 (Reference)";
  tab$estimate <- as.numeric(gsub(" .*", "", tab[,2]));
  tab$conf.low <- as.numeric(gsub(",.*", "", gsub("\\(([^()]*)\\)|.", "\\1", tab[,2])));
  tab$conf.high <- as.numeric(gsub("^\\S*\\s+", "", gsub("\\(([^()]*)\\)|.", "\\1", tab[,2])));
  tab$level.name <- tab[,1];
  tab$var.name <- NA;
  tab$var.name[which(tab$Covariate %in% covs)] <- tab$level.name[which(tab$Covariate %in% covs)];
  y <- tab$var.name;
  y_forward_fill <- fillNAs(y);
  tab <- cbind(tab, y, y_forward_fill);
  tab$var.name <- tab$y_forward_fill;
  tab$level.order <- sequence(rle(tab$var.name)$lengths)
  tab <- tab[order(rank(tab$estimate), tab$var.name), ];
  dt <- as.data.frame(unique(tab$var.name));
  colnames(dt) <- "var.name";
  dt$var.order <- 1:nrow(dt);
  dt$var.order <- dt$var.order + 1;
  tab <- merge(tab, dt, by = "var.name", all = T);
  tab <- tab[order(tab$var.order, -tab$level.order), ];
  tab$p.value <- tab$"p-value";
  tab$p.label <- paste(format(round(as.numeric(tab$p.value), 3), nsmall=3), sep="");
  tab$variable <- tab$Variable;
  tab <- tab[, c("variable", "var.name", "level.name", "level.order", "estimate", "p.label", "p.value", "conf.low", "conf.high", "var.order", "estimate.label", "N", "Event")];
  tab <- as.data.frame(tab);
  if (rmRef)
    tab = tab[setdiff(1:nrow(tab), which(tab$estimate.label ==
                                           "1.0 (Reference)")), ]
  yvals = 1:nrow(tab)
  tab$estimate.label = ifelse(is.na(tab$estimate.label), "",
                              tab$estimate.label)
  tab$estimate.label = ifelse(tab$estimate.label == "1.0 (Reference)",
                              "(Reference)", tab$estimate.label)
  if (showEst) {
    yLabels = data.frame(y.pos = yvals, labels = ifelse(is.na(tab$level.name),
                                                        paste(tab$variable, ": ", tab$estimate.label, sep=""),
                                                        paste(tab$level.name, ": ", tab$estimate.label, sep="")))
  }
  else {
    yLabels = data.frame(y.pos = yvals, labels = ifelse(is.na(tab$level.name),
                                                        tab$variable, ifelse(tab$estimate.label == "(Reference)",
                                                                             paste(tab$level.name, ": ", tab$estimate.label, sep=""), tab$level.name)))
  }
  yLabels$labels <- gsub("_", " ", yLabels$labels)
  yLabels <- yLabels[!is.na(yLabels$y.pos), ]
  tab$x.val = ifelse(tab$estimate.label == "(Reference)", 1,
                     tab$estimate)
  tab$y.val = yLabels$y.pos
  tab$colour <- ifelse(tab$x.val < 1, "a", ifelse(tab$x.val ==
                                                    1, "b", "c"))
  if (length(colours) == 1) {
    colours = c(a = "#006B3C", b = "black", c = "#FF0800")
  }
  else {
    names(colours) = c("a", "b", "c")
  }
  if (showN & !showEvent) {
    Axis = scale_y_continuous(breaks = yLabels$y.pos,
                              labels = yLabels$labels, sec.axis = dup_axis(breaks = yLabels$y.pos,
                                                                           labels = tab$N, name = "N"))
    themeSecAxis = theme(axis.title.y.right = element_text(angle = 0, hjust = 0.5, vjust = 0.5))
  }
  if (showEvent) {
    Axis = scale_y_continuous(breaks = yLabels$y.pos,
                              labels = yLabels$labels, sec.axis = dup_axis(breaks = yLabels$y.pos,
                                                                           labels = paste(tab$N, tab$Event, sep=" : "), name = "N : Event"))
    themeSecAxis = theme(axis.title.y.right = element_text(angle = 270, hjust = 0.5, vjust = 0.5))
  }
  if (!showN & !showEvent) {
    Axis = scale_y_continuous(breaks = yLabels$y.pos,
                              labels = yLabels$labels)
    themeSecAxis = NULL;
  }
  colours <- colours[sort(unique(tab$colour))]
  suppressWarnings({tryCatch({
    p = ggplot(tab, aes_(x = ~x.val, y = ~y.val, colour = ~colour)) +
      geom_point(na.rm = TRUE, size = 2) +
      geom_errorbarh(aes_(xmin = ~conf.low,
                          xmax = ~conf.high), height = 0, size = 0.9, na.rm = TRUE) +
      geom_vline(xintercept = 1) + labs(y = "", x = x_lab) +
      guides(colour = "none") + Axis + scale_colour_manual(values = colours) +
      theme_bw() +
      theme(axis.text.y = element_text(face = ifelse(tab$variable ==
                                                       tab$var.name | is.na(tab$var.name), "bold", "plain"),
                                       hjust = 0),
            panel.grid.major.y = element_blank(),
            panel.grid.minor.y = element_blank(),
            axis.ticks = element_blank()) +
      themeSecAxis
    if (logScale)
      p + scale_x_log10(breaks = scales::log_breaks(n = nxTicks))
    else p
  }, error=function(e){})})
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
#' @noRd
#' @examples
#' data("pembrolizumab")
#' glm_fit = glm(orr~change_ctdna_group+sex+age+l_size,
#' data=pembrolizumab,family = 'binomial')
#' forestplotMV(glm_fit)
forestplotMV = function (model, data,conf.level = 0.95, orderByRisk = TRUE,
                         colours = "default", showEst = TRUE, rmRef = FALSE,
                         digits=getOption("reportRmd.digits",2),
                         logScale=getOption("reportRmd.logScale",TRUE), nxTicks = 5, showN = TRUE, showEvent = TRUE)
{
  if (inherits(model, "glm")) {
    if (model$family$link == "log") {
      x_lab = "Adjusted Relative Risk"
    }
    else if (model$family$link == "logit") {
      x_lab = "Adjusted Odds Ratio"
    }
    else stop("model must be a logit or log link fit")
  }
  else {
    x_lab = "Adjusted Odds Ratio"
  }
  tab = mvsum(model, data, digits = digits, markup = F, sanitize = F,
              nicenames = F, showN = TRUE, showEvent = showEvent, CIwidth = conf.level)
  tab$estimate.label <- tab[,2];
  tab$estimate.label[which(tab$estimate.label == "Reference")] <- "1.0 (Reference)";
  tab$estimate <- as.numeric(gsub(" .*", "", tab[,2]));
  tab$conf.low <- as.numeric(gsub(",.*", "", gsub("\\(([^()]*)\\)|.", "\\1", tab[,2])));
  tab$conf.high <- as.numeric(gsub("^\\S*\\s+", "", gsub("\\(([^()]*)\\)|.", "\\1", tab[,2])));
  tab$level.name <- tab[,1];
  tab$var.name <- NA;
  covs <- colnames(model$model);
  tab$var.name[which(tab$Covariate %in% covs)] <- tab$level.name[which(tab$Covariate %in% covs)];
  y <- tab$var.name;
  y_forward_fill <- fillNAs(y);
  tab <- cbind(tab, y, y_forward_fill);
  tab$var.name <- tab$y_forward_fill;
  tab$level.order <- sequence(rle(tab$var.name)$lengths)
  tab <- tab[order(rank(tab$estimate), tab$var.name), ];
  dt <- as.data.frame(unique(tab$var.name));
  colnames(dt) <- "var.name";
  dt$var.order <- 1:nrow(dt);
  dt$var.order <- dt$var.order + 1;
  tab <- merge(tab, dt, by = "var.name", all = T);
  tab <- tab[order(tab$var.order, -tab$level.order), ];
  tab$p.value <- tab$"p-value";
  tab$p.label <- paste(format(round(as.numeric(tab$p.value), 3), nsmall=3), sep="");
  tab$variable <- tab$Covariate;
  vars <- c("variable", "var.name", "level.name", "level.order", "estimate", "p.label", "p.value", "conf.low", "conf.high", "var.order", "estimate.label", "N", "Event")
  if (!('Event' %in% names(tab))) {
    showEvent <- FALSE
    vars <- setdiff(vars,"Event")
  }
  tab <- tab[, vars];
  tab <- as.data.frame(tab);

  if (rmRef)
    tab = tab[setdiff(1:nrow(tab), which(tab$estimate.label ==
                                           "1.0 (Reference)")), ]
  yvals = 1:nrow(tab)
  tab$estimate.label = ifelse(is.na(tab$estimate.label), "",
                              tab$estimate.label)
  tab$estimate.label = ifelse(tab$estimate.label == "1.0 (Reference)",
                              "(Reference)", tab$estimate.label)
  if (showEst) {
    yLabels = data.frame(y.pos = yvals, labels = ifelse(is.na(tab$level.name),
                                                        paste(tab$variable, ": ", tab$estimate.label, sep=""),
                                                        paste(tab$level.name, ": ", tab$estimate.label, sep="")))
  }
  else {
    yLabels = data.frame(y.pos = yvals, labels = ifelse(is.na(tab$level.name),
                                                        tab$variable, ifelse(tab$estimate.label == "(Reference)",
                                                                             paste(tab$level.name, ": ", tab$estimate.label, sep=""), tab$level.name)))
  }
  yLabels$labels <- gsub("_", " ", yLabels$labels)
  yLabels <- yLabels[!is.na(yLabels$y.pos), ]
  tab$x.val = ifelse(tab$estimate.label == "(Reference)", 1,
                     tab$estimate)
  tab$y.val = yLabels$y.pos
  tab$colour <- ifelse(tab$x.val < 1, "a", ifelse(tab$x.val ==
                                                    1, "b", "c"))
  if (length(colours) == 1) {
    colours = c(a = "#006B3C", b = "black", c = "#FF0800")
  }
  else {
    names(colours) = c("a", "b", "c")
  }
  if (showN & !showEvent) {
    Axis = scale_y_continuous(breaks = yLabels$y.pos,
                              labels = yLabels$labels, sec.axis = dup_axis(breaks = yLabels$y.pos,
                                                                           labels = tab$N, name = "N"))
    themeSecAxis = theme(axis.title.y.right = element_text(angle = 0, hjust = 0.5, vjust = 0.5))
  }
  if (showEvent) {
    Axis = scale_y_continuous(breaks = yLabels$y.pos,
                              labels = yLabels$labels, sec.axis = dup_axis(breaks = yLabels$y.pos,
                                                                           labels = paste(tab$N, tab$Event, sep=" : "), name = "N : Event"))
    themeSecAxis = theme(axis.title.y.right = element_text(angle = 270, hjust = 0.5, vjust = 0.5))
  }
  if (!showN & !showEvent) {
    Axis = scale_y_continuous(breaks = yLabels$y.pos,
                              labels = yLabels$labels)
    themeSecAxis = NULL;
  }
  colours <- colours[sort(unique(tab$colour))]
  suppressWarnings({tryCatch({
    p = ggplot(tab, aes_(x = ~x.val, y = ~y.val, colour = ~colour)) +
      geom_point(na.rm = TRUE, size = 2) + geom_errorbarh(aes_(xmin = ~conf.low,
                                                               xmax = ~conf.high), height = 0, size = 0.9, na.rm = TRUE) +
      geom_vline(xintercept = 1) + labs(y = "", x = x_lab) +
      guides(colour = "none") + Axis + scale_colour_manual(values = colours) +
      theme_bw() + theme(axis.text.y = element_text(face = ifelse(tab$variable ==
                                                                    tab$var.name | is.na(tab$var.name), "bold", "plain"),
                                                    hjust = 0), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
                         axis.ticks = element_blank()) + themeSecAxis
    if (logScale)
      p + scale_x_log10(breaks = scales::log_breaks(n = nxTicks))
    else p
  }, error=function(e){})})
}


#' Combine an univariable and multivariable forest plot using ggplot2
#'
#' This function will take log or logistic regression fit forest plot output
#' from forestplotUV and forestplotMV functions and display the combined
#' adjusted and unadjusted OR or RR for each variable on the appropriate
#' log scale. Please note that total N and reference-level N is taken from
#' unadjusted model.
#'
#' @param UVmodel an UV model object output from the forestplotUV function
#' @param MVmodel a MV model object output from the forestplotMV function
#' @param model fitted model object
#' @param family description of the error distribution and link function to be
#'   used in the model. Only used for geeglm
#' @param digits number of digits to round to
#' @param orderByRisk logical, should the plot be ordered by risk
#' @param colours can specify colours for risks less than, 1 and greater than
#'   1.0. Default is red, black, green
#' @param showEst logical, should the risks be displayed on the plot in text
#' @param rmRef logical, should the reference levels be removed for the plot?
#' @param logScale logical, should OR/RR be shown on log scale, defaults to
#'   TRUE. See https://doi.org/10.1093/aje/kwr156 for why you may prefer a
#'   linear scale.
#' @param nxTicks Number of tick marks supplied to the log_breaks function to
#'   produce
#' @param showN Show number of observations per variable and category
#' @param showEvent Show number of events per variable and category
#' @import ggplot2
#' @importFrom scales log_breaks
#' @keywords plot
#' @return a plot object
#' @noRd
#' @examples
#' data("pembrolizumab")
#' UVp = forestplotUV(response="orr", covs=c("change_ctdna_group", "sex", "age",
#' "l_size"), data=pembrolizumab, family='binomial')
#' MVp = forestplotMV(glm(orr~change_ctdna_group+sex+age+l_size,
#' data=pembrolizumab,family = 'binomial'))
#' forestplotUVMV(UVp, MVp)
forestplotUVMV = function (UVmodel, MVmodel, model = "glm",
                           family = NULL, digits=getOption("reportRmd.digits",2), orderByRisk = TRUE, colours = "default",
                           showEst = TRUE, rmRef = FALSE, logScale = FALSE, nxTicks = 5,
                           showN = TRUE, showEvent = TRUE)
{
  if (inherits(model, "glm")) {
    if (model$family$link == "log") {
      x_lab = "Relative Risk"
    }
    else if (model$family$link == "logit") {
      x_lab = "Odds Ratio"
    }
    else stop("model must be a logit or log link fit")
  }
  else {
    x_lab = "Odds Ratio"
  }
  UVmodel$data$type <- "Unadjusted";
  MVmodel$data$type <- "Adjusted";
  tab <- rbind(UVmodel$data, MVmodel$data);
  tab$IsVarName <- paste(tab$estimate.label, tab$type, sep="");
  tab[which(tab$estimate.label != "(Reference)" & is.na(tab$estimate)), ]$IsVarName <- "Y";
  tab <- tab[!duplicated(tab[ , c("var.name", "variable", "IsVarName")]),];
  tab$Reference <- paste(tab$estimate.label, tab$type, sep="");
  tab[which(tab$estimate.label == "(Reference)"), ]$Reference <- "Y";
  tab <- tab[!duplicated(tab[, c("var.name", "estimate.label", "Reference")]),];
  y <- tab$var.name;
  y_forward_fill <- fillNAs(y);
  tab <- cbind(tab, y, y_forward_fill);
  tab$var.name <- tab$y_forward_fill;
  tab$level.order <- sequence(rle(tab$var.name)$lengths)
  tab <- tab[order(rank(tab$estimate), tab$var.name), ];
  dt <- as.data.frame(unique(tab$var.name));
  colnames(dt) <- "var.name";
  dt$var.order <- 1:nrow(dt);
  dt$var.order <- dt$var.order + 1;
  tab <- tab[, -which(colnames(tab) %in% c("var.order"))];
  tab <- merge(tab, dt, by = "var.name", all = T);
  tab <- tab[order(tab$var.order, tab$level.order, tab$type, decreasing = c(F, F, T), method="radix"), ];
  tab <- tab[, c("variable", "var.name", "level.name", "level.order", "estimate", "p.label", "p.value", "conf.low", "conf.high", "var.order", "estimate.label", "N", "Event", "type")];
  tab <- as.data.frame(tab);
  if (rmRef)
    tab = tab[setdiff(1:nrow(tab), which(tab$estimate.label ==
                                           "1.0 (Reference)")), ]
  yvals = 1:nrow(tab)
  tab$estimate.label = ifelse(is.na(tab$estimate.label), "",
                              tab$estimate.label)
  tab$estimate.label = ifelse(tab$estimate.label == "1.0 (Reference)",
                              "(Reference)", tab$estimate.label)
  if (showEst) {
    yLabels = data.frame(y.pos = yvals, labels = ifelse(is.na(tab$level.name),
                                                        paste(tab$variable, ": ", tab$estimate.label, sep=""),
                                                        paste(tab$level.name, ": ", tab$estimate.label, sep="")))
  }
  else {
    yLabels = data.frame(y.pos = yvals, labels = ifelse(is.na(tab$level.name),
                                                        tab$variable, ifelse(tab$estimate.label == "(Reference)",
                                                                             paste(tab$level.name, ": ", tab$estimate.label, sep=""), tab$level.name)))
  }
  yLabels$labels <- gsub("_", " ", yLabels$labels)
  yLabels <- yLabels[!is.na(yLabels$y.pos), ]
  tab$x.val = ifelse(tab$estimate.label == "(Reference)", 1,
                     tab$estimate)
  tab$y.val = yLabels$y.pos
  tab$colour <- ifelse(tab$x.val < 1, "a", ifelse(tab$x.val ==
                                                    1, "b", "c"))
  if (length(colours) == 1) {
    colours = c(a = "#006B3C", b = "black", c = "#FF0800")
  }
  else {
    names(colours) = c("a", "b", "c")
  }
  if (showN & !showEvent) {
    Axis = scale_y_continuous(breaks = yLabels$y.pos,
                              labels = yLabels$labels, sec.axis = dup_axis(breaks = yLabels$y.pos,
                                                                           labels = tab$N, name = "N"))
    themeSecAxis = theme(axis.title.y.right = element_text(angle = 0, hjust = 0.5, vjust = 0.5))
    warning(paste("Total N and reference-level N is taken from unadjusted model."))
  }
  if (showEvent) {
    Axis = scale_y_continuous(breaks = yLabels$y.pos,
                              labels = yLabels$labels, sec.axis = dup_axis(breaks = yLabels$y.pos,
                                                                           labels = paste(tab$N, tab$Event, sep=" : "), name = "N : Event"))
    themeSecAxis = theme(axis.title.y.right = element_text(angle = 270, hjust = 0.5, vjust = 0.5))
    warning(paste("Total N and reference-level N is taken from unadjusted model."))
  }
  if (!showN & !showEvent) {
    Axis = scale_y_continuous(breaks = yLabels$y.pos,
                              labels = yLabels$labels)
    themeSecAxis = NULL;
  }
  colours <- colours[sort(unique(tab$colour))]
  suppressWarnings({tryCatch({
    p = ggplot(tab, aes_(x = ~x.val, y = ~y.val, colour = ~colour, linetype = ~type, shape = ~type)) +
      geom_point(na.rm = TRUE, size = 2) + geom_errorbarh(aes_(xmin = ~conf.low,
                                                               xmax = ~conf.high), height = 0, size = 0.9, na.rm = TRUE) +
      geom_vline(xintercept = 1) + labs(y = "", x = x_lab) +
      guides(colour = "none") + Axis + scale_colour_manual(values = colours) +
      scale_shape_manual(values = c(15, 16)) +
      theme_bw() + theme(axis.text.y = element_text(face = ifelse(tab$variable ==
                                                                    tab$var.name | is.na(tab$var.name), "bold", "plain"),
                                                    hjust = 0), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
      theme(legend.title=element_blank(), legend.margin = margin(-8, 0, 0, 0),
            legend.spacing.x = unit(2, "mm"), legend.spacing.y = unit(0, "mm"),
            legend.position = "bottom", axis.ticks = element_blank()) +
      themeSecAxis
    if (logScale)
      p + scale_x_log10(breaks = scales::log_breaks(n = nxTicks))
    else p
  }, error=function(e){})})
}
