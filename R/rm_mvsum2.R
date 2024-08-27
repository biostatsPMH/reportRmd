#'Format a regression model nicely for 'Rmarkdown'
#'
#'Multivariable (or univariate) regression models are re-formatted for reporting
#'and a global p-value is added for the evaluation of factor variables.
#'
#'Global p-values are likelihood ratio tests for lm, glm and polr models. For
#'lme models an attempt is made to re-fit the model using ML and if,successful
#'LRT is used to obtain a global p-value. For coxph models the model is re-run
#'without robust variances with and without each variable and a LRT is
#'presented. If unsuccessful a Wald p-value is returned. For GEE and CRR models
#'Wald global p-values are returned. For negative binomial models a deviance
#'test is used.
#'
#'If the variance inflation factor is requested (VIF=T) then a generalised VIF
#'will be calculated in the same manner as the car package.
#'
#'If the MASS package is loaded, profile likelihood confidence intervals will be
#'calculated; otherwise Wald confidence intervals will be calculated.
#'Users should look for the message "Waiting for profiling to be done...", which
#'indicates that profile likelihoods are calculated.
#'
#'The number of decimals places to display the statistics can be changed with
#'digits, but this will not change the display of p-values. If more significant
#'digits are required for p-values then use tableOnly=TRUE and format as
#'desired.
#'@param model model fit
#'@param data data that model was fit on (an attempt will be made to extract
#'  this from the model)
#'@param digits number of digits to round estimates to, does not affect p-values
#'@param covTitle character with the names of the covariate (predictor) column.
#'  The default is to leave this empty for output or, for table only output to
#'  use the column name 'Covariate'.
#'@param showN boolean indicating sample sizes should be shown for each
#'  comparison, can be useful for interactions
#'@param showEvent boolean indicating if number of events should be shown. Only
#'  available for logistic.
#'@param CIwidth width for confidence intervals, defaults to 0.95
#'@param vif boolean indicating if the variance inflation factor should be
#'  included. See details
#'@param whichp string indicating whether you want to display p-values for
#'  levels within categorical data ("levels"), global p values ("global"), or
#'  both ("both"). Irrelevant for continuous predictors.
#'@param caption table caption
#'@param tableOnly boolean indicating if unformatted table should be returned
#'@param p.adjust p-adjustments to be performed (Global p-values only)
#'@param unformattedp boolean indicating if you would like the p-value to be
#'  returned unformatted (ie not rounded or prefixed with '<'). Should be used
#'  in conjuction with the digits argument.
#'@param nicenames boolean indicating if you want to replace . and _ in strings
#'  with a space
#'@param chunk_label Deprecated, previously used in Word to allow
#'  cross-referencing, this should now be done at the chunk level.
#'@param fontsize PDF/HTML output only, manually set the table fontsize
#'@return A character vector of the table source code, unless tableOnly=TRUE in
#'  which case a data frame is returned
#'@export
#'@references John Fox & Georges Monette (1992) Generalized Collinearity
#'  Diagnostics, Journal of the American Statistical Association, 87:417,
#'  178-183, \doi{10.1080/01621459.1992.10475190}
#'@references  John Fox and Sanford Weisberg (2019). An {R} Companion to Applied
#'  Regression, Third Edition. Thousand Oaks CA: Sage.
#' @examples
#' data("pembrolizumab")
#' glm_fit = glm(change_ctdna_group~sex:age+baseline_ctdna+l_size,
#' data=pembrolizumab,family = 'binomial')
#' rm_mvsum(glm_fit)
#'
#' #linear model with p-value adjustment
#' lm_fit=lm(baseline_ctdna~age+sex+l_size+tmb,data=pembrolizumab)
#' rm_mvsum(lm_fit,p.adjust = "bonferroni")
#' #Coxph
#' require(survival)
#' res.cox <- coxph(Surv(os_time, os_status) ~ sex+age+l_size+tmb, data = pembrolizumab)
#' rm_mvsum(res.cox, vif=TRUE)
rm_mvsum2 <- function(model, data, digits=getOption("reportRmd.digits",2),covTitle='',showN=TRUE,showEvent=TRUE,CIwidth=0.95, vif=FALSE,
                     whichp=c("levels","global","both"),
                     caption=NULL,tableOnly=FALSE,p.adjust='none',unformattedp=FALSE,nicenames = TRUE,chunk_label, fontsize){
  if (unformattedp) formatp <- function(x) {as.numeric(x)}
  whichp <- match.arg(whichp)

  lifecycle::deprecate_soft("0.2.0","rm_mvsum2(data)")
  lifecycle::deprecate_soft("0.2.0","rm_mvsum2(chunk_label)")
  if (any(is.na(model$coefficients))) stop(paste0('rm_mvsum2 cannot run when any model coeffcients are NA.\nThe following model coefficients could not be estimated:\n',
                                                  paste(names(model$coefficients)[is.na(model$coefficients)],collapse = ", "),
                                                  "\nPlease re-fit a valid model prior to reporting. Do you need to run droplevels?"))
  # get the table
  tab <- m_summary(model, CIwidth = CIwidth, digits = digits, vif = vif, whichp = whichp, for_plot = FALSE)
  att_tab <- attributes(tab)
  to_indent <- which(!(tab[["Variable"]] %in% attr(model$terms, "term.labels")))
  bold_cells <- cbind(which(tab[["Variable"]] %in% attr(model$terms, "term.labels")), rep(1, length(which(tab[["Variable"]] %in% attr(model$terms, "term.labels")))))
  attr(tab, "to_indent") <- to_indent
  attr(tab, "bold_cells") <- bold_cells


  # perform p-value adjustment across variable-level p-values remove factor p-values
  # if ("Global p-value" %in% names(tab)){
  #   if(p.adjust!='none') {
  #     raw_p <- ifelse(is.na(tab[["Global p-value"]]),tab[["p-value"]],tab[["Global p-value"]])
  #     raw_p[!att_tab$varID] <- NA
  #     p_sig <- suppressWarnings(stats::p.adjust(raw_p,method=p.adjust))
  #     message('Global p-values were adjusted according to the ',p.adjust,' method. Factor level p-values have been removed.')
  #     tab[["raw p-value"]]<-formatp(raw_p)
  #   } else{
  #     p_sig <- ifelse(is.na(tab[["Global p-value"]]),tab[["p-value"]],tab[["Global p-value"]])
  #   }
  #   tab[["p-value"]]  <- sapply(p_sig,formatp)
  #   tab <- tab[,grep("Global p-value",names(tab),invert = T)]
  # } else {
  #   raw_p <- tab[["p-value"]]
  #   p_sig <- suppressWarnings(stats::p.adjust(raw_p,method=p.adjust))
  #   tab[["p-value"]] <- sapply(p_sig,formatp)
  # }
  # to_bold_p <- which(as.numeric(tab[["p-value"]])<.05)


  # if (length(to_bold_p)>0)  bold_cells <- rbind(bold_cells,
  #                                               matrix(cbind(to_bold_p, which(names(tab)=='p-value')),ncol=2))


  # formatting pval column
  method <- p.adjust
  tab[["p-value"]] <- p.adjust(tab[["p-value"]], method = method)
  if (!unformattedp) {
    tab[["p-value"]] <- formatp(tab[["p-value"]])
  }
  p_col <- (which(names(tab) == "p-value"))
  bold_cells <- rbind(bold_cells, cbind(which(as.numeric(gsub("[^0-9\\.]", "", tab[["p-value"]])) < 0.05), rep(p_col, length(which(as.numeric(gsub("[^0-9\\.]", "", tab[["p-value"]])) < 0.05)))))
  if (nrow(bold_cells) < 1) {
    bold_cells <- NULL
  }

  if (nicenames){
    attr(tab,"termnames") <- tab$Variable
    tab$Variable <- replaceLbl(model$call[["data"]], tab$Variable)
  }

  names(tab)[1] <- covTitle
  for (a in setdiff(names(att_tab),names(attributes(tab)))) attr(tab,a) <- att_tab[[a]]
  if (tableOnly){
    if (names(tab)[1]=='') names(tab)[1] <- 'Covariate'
    attr(tab, 'to_indent') <- to_indent
    attr(tab,'bold_cells') <- bold_cells
    attr(tab,'dimchk') <- dim(tab)
    return(tab)
  }
  argL <- list(tab=tab,to_indent=to_indent,bold_cells = bold_cells,
               caption=caption, digits = digits,
               chunk_label=ifelse(missing(chunk_label),'NOLABELTOADD',chunk_label))
  if (!missing(fontsize)) argL[['fontsize']] <- fontsize
  do.call(outTable, argL)
}
