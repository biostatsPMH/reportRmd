# OLD rm_uvsum, rm_mvsum for reference
#'Output several univariate models nicely in a single table
#'
#'A table with the model parameters from running separate univariate models on
#'each covariate. For factors with more than two levels a Global p-value is
#'returned.
#'
#'Global p-values are likelihood ratio tests for lm, glm and polr models. For
#'lme models an attempt is made to re-fit the model using ML and if,successful
#'LRT is used to obtain a global p-value. For coxph models the model is re-run
#'without robust variances with and without each variable and a LRT is
#'presented. If unsuccessful a Wald p-value is returned. For GEE and CRR models
#'Wald global p-values are returned.
#'
#'The number of decimals places to display the statistics can be changed with
#'digits, but this will not change the display of p-values. If more significant
#'digits are required for p-values then use tableOnly=TRUE and format as
#'desired.
#'@param response string vector with name of response
#'@param covs character vector with the names of columns to fit univariate
#'  models to
#'@param data dataframe containing data
#'@param digits number of digits to round estimates and CI to. Does not affect
#'  p-values.
#'@param covTitle character with the names of the covariate (predictor) column.
#'  The default is to leave this empty for output or, for table only output to
#'  use the column name 'Covariate'.
#'@param caption character containing table caption (default is no caption)
#'@param tableOnly boolean indicating if unformatted table should be returned
#'@param removeInf boolean indicating if infinite estimates should be removed
#'  from the table
#' @param p.adjust p-adjustments to be performed. Uses the
#'  [p.adjust] function from base R
#'@param unformattedp boolean indicating if you would like the p-value to be
#'  returned unformatted (ie not rounded or prefixed with '<'). Should be used
#'  in conjunction with the digits argument.
#'@param whichp string indicating whether you want to display p-values for
#'  levels within categorical data ("levels"), global p values ("global"), or
#'  both ("both"). Irrelevant for continuous predictors.
#'@param chunk_label only used if output is to Word to allow cross-referencing
#'@param  gee boolean indicating if gee models should be fit to account for
#'  correlated observations. If TRUE then the id argument must specify the
#'  column in the data which indicates the correlated clusters.
#'@param id character vector which identifies clusters. Only used for geeglm
#'@param corstr character string specifying the correlation structure. Only used
#'  for geeglm. The following are permitted: '"independence"', '"exchangeable"',
#'  '"ar1"', '"unstructured"' and '"userdefined"'
#'@param family description of the error distribution and link function to be
#'  used in the model. Only used for geeglm
#'@param type string indicating the type of univariate model to fit. The
#'  function will try and guess what type you want based on your response. If
#'  you want to override this you can manually specify the type. Options include
#'  "linear", "logistic", "poisson",coxph", "crr", "boxcox", "ordinal", "geeglm"
#'@param offset string specifying the offset term to be used for Poisson or
#'  negative binomial regression. Example: offset="log(follow_up)"
#'@param strata character vector of covariates to stratify by. Only used for
#'  coxph and crr
#'@param nicenames boolean indicating if you want to replace . and _ in strings
#'  with a space
#'@param showN boolean indicating if you want to show sample sizes
#'@param showEvent boolean indicating if you want to show number of events. Only
#'  available for logistic.
#'@param CIwidth width of confidence interval, default is 0.95
#'@param reflevel manual specification of the reference level. Only used for
#'  ordinal regression This will allow you to see which model is not fitting if
#'  the function throws an error
#'@param returnModels boolean indicating if a list of fitted models should be
#'  returned. If this is TRUE then the models will be returned, but the output
#'  will be suppressed. In addition to the model elements a data element will be
#'  appended to each model so that the fitted data can be examined, if
#'  necessary. See Details
#'@param fontsize PDF/HTML output only, manually set the table fontsize
#'@param forceWald `r lifecycle::badge("deprecated")` `forceWald = TRUE` is no longer supported; this function will always use profile likelihoods as per the inclusion of the MASS confidence intervals into base from from R 4.4.0
#' @seealso \code{\link{covsum}},\code{\link{fisher.test}},
#'   \code{\link{chisq.test}}, \code{\link{wilcox.test}},
#'   \code{\link{kruskal.test}}, \code{\link{anova}}, \code{\link[rstatix:cramer_v]{rstatix::cramer_v}},
#'   \code{\link[rstatix:eta_squared]{rstatix:eta_squared}}, and \code{\link{outTable}}
#'@return A character vector of the table source code, unless tableOnly=TRUE in
#'  which case a data frame is returned
#'@export
#' @examples
#' # Examples are for demonstration and are not meaningful
#' # Coxph model with 90% CI
#' data("pembrolizumab")
#' rm_uvsum(response = c('os_time','os_status'),
#' covs=c('age','sex','baseline_ctdna','l_size','change_ctdna_group'),
#' data=pembrolizumab,CIwidth=.9)
#'
#' # Linear model with default 95% CI
#' rm_uvsum(response = 'baseline_ctdna',
#' covs=c('age','sex','l_size','pdl1','tmb'),
#' data=pembrolizumab)
#'
#' # Logistic model with default 95% CI
#' rm_uvsum(response = 'os_status',
#' covs=c('age','sex','l_size','pdl1','tmb'),
#' data=pembrolizumab,family = binomial)

#' # Poisson models returned as model list
#' mList <- rm_uvsum(response = 'baseline_ctdna',
#' covs=c('age','sex','l_size','pdl1','tmb'),
#' data=pembrolizumab, returnModels=TRUE)
#' #'
#' # GEE on correlated outcomes
#' data("ctDNA")
#' rm_uvsum(response = 'size_change',
#' covs=c('time','ctdna_status'),
#' gee=TRUE,
#' id='id', corstr="exchangeable",
#' family=gaussian("identity"),
#' data=ctDNA,showN=TRUE)
rm_uvsum <- function(response, covs , data , digits=getOption("reportRmd.digits",2), covTitle='',caption=NULL,
                     tableOnly=FALSE,removeInf=FALSE,p.adjust='none',unformattedp=FALSE,
                     whichp=c("levels","global","both"),
                     chunk_label,
                     gee=FALSE,id = NULL,corstr = NULL,family = NULL,type = NULL,
                     offset=NULL,
                     strata = 1,
                     nicenames = TRUE,showN=TRUE,showEvent=TRUE,CIwidth = 0.95,
                     reflevel=NULL,returnModels=FALSE,fontsize,forceWald=FALSE){

  if (missing(data)) stop('data is a required argument')
  if (missing(covs)) stop('covs is a required argument') else covs <- unique(covs)
  if (missing(response)) stop('response is a required argument')
  if (length(response)>2) stop('The response must be a single outcome for linear, logistic and ordinal models or must specify the time and event status variables for survival models.')
  if (!inherits(data,'data.frame')) stop('data must be supplied as a data frame.')
  if (!inherits(covs,'character')) stop('covs must be supplied as a character vector or string indicating variables in data')
  if (!missing(id)) if (!is.null(id)) if(!inherits(id,"character") | length(id)>1) stop("id must be specified as a character referring to a variable name, id='varname'")
  if (is.null(strata)) strata <- 1
  if ( is.na(strata) | strata=="") strata <- 1
  missing_vars = na.omit(setdiff(c(response, covs,id,ifelse(strata==1,NA,strata)), names(data)))
  if (length(missing_vars) > 0) stop(paste("These variables are not in the data:\n",
                                           paste0(missing_vars,collapse=csep())))
  if (strata==1) nm <- c(response,covs) else nm <- c(strata,response,covs)
  if (!all(names(data[,nm])==names(data.frame(data[,nm])))) stop('Non-standard variable names detected.\n Try converting data with new_data <- data.frame(data) \n then use new variable names in rm_uvsum.' )
  if (missing(forceWald)) forceWald = getOption("reportRmd.forceWald",FALSE)
  argList <- as.list(match.call()[-1])
  df_nm <- matchdata(argList$data)
  whichp <- match.arg(whichp)
  for (v in covs) {
    if (inherits(data[[v]], c("character", "ordered"))) data[[v]] <- factor(data[[v]], ordered = F)
    if (inherits(data[[v]],c('Date','POSIXt'))) {
      covs <- setdiff(covs,v)
      message(paste('Dates can not be used as predictors, try creating a time variable.\n The variable',v,'does not appear in the table.'))
    }

    df <- na.omit(data[,c(response,v)])
    if (v %in% response){
      warning(paste(v,'is the response and can not appear in the covariate.\n',
                    'It is omitted from the output.'))
      covs <- setdiff(covs,v)
    }
    if (length(unique(df[[v]]))==1) {
      warning(paste(v,'has only one unique value for non-missing response.\n',
                    'It is omitted from the output.'))
      covs <- setdiff(covs,v)
    }
  }

  if (unformattedp) formatp <- function (x,...){x}
  # get the table
  rtn <- uvsum(response,covs,data,digits=digits,markup = FALSE,sanitize=FALSE,
               gee=gee,id = id, offset=offset,
               corstr = corstr,family = family,type = type,strata = strata,
               nicenames = FALSE,showN = showN,showEvent = showEvent,
               CIwidth = CIwidth,reflevel=reflevel,returnModels=returnModels,forceWald = forceWald)
  if (returnModels) tab <- rtn[[1]] else tab <- rtn
  att_tab <- attributes(tab)
  cap_warn <- character(0)
  if (removeInf){
    # Do not display unstable estimates
    inf_values =  grep('Inf',tab[,2])
    if (length(inf_values)>0){
      if ('Global p-values' %in% names(tab)) to_hide <-2:4 else to_hide <-2:3
      tab[inf_values,to_hide] <-NA
      cap_warn <- paste0(cap_warn,ifelse(identical(cap_warn,character(0)),'',', '),
                         'Covariates with unstable estimates:',
                         paste(tab$Covariate[inf_values],collapse=','),'.')
    }
  }
  # if an adjustment was made, add this to the cap_warn text
  if (p.adjust!='none') cap_warn <- paste0(cap_warn,'. Global p-values were adjusted according to the ',p.adjust,' method. Factor level p-values have been removed.')

  to_indent <- which(!attr(tab,"varID"))
  to_bold_name <- which(attr(tab,"varID"))
  bold_cells <- arrayInd(to_bold_name, dim(tab))

  if (nicenames){
    attr(tab,"termnames") <- tab$Covariate
    tab$Covariate <- replaceLbl(df_nm, tab$Covariate)
  }
  # decide which p-values to keep
  if (whichp=="levels"){
    if ("Global p-value" %in% names(tab)) tab[["Global p-value"]] <- NULL
  } else if (whichp=="global"){
    if ("Global p-value" %in% names(tab)) {
      tab$`p-value` <- sapply(1:nrow(tab), function(x) {
        ifelse(att_tab$varID[x],
               ifelse(tab[["Global p-value"]][x]!="",tab[["Global p-value"]][x],tab$`p-value`[x]),"")
      })
      tab[["Global p-value"]] <- NULL
    }
  } # if both then leave as is

  if ("Global p-value" %in% names(tab)){
    tab[["Global p-value"]][which(tab[["Global p-value"]]==''|tab[["Global p-value"]]=='NA')] <-NA
    if(p.adjust!='none') {
      raw_p <- ifelse(is.na(tab[["Global p-value"]]),tab[["p-value"]],tab[["Global p-value"]])
      raw_p[!att_tab$varID] <- NA
      p_sig <- suppressWarnings(stats::p.adjust(raw_p,method=p.adjust))
      message('Global p-values were adjusted according to the ',p.adjust,' method. Factor level p-values have been removed.')
      tab[["raw p-value"]]<-formatp(raw_p)
    } else{
      p_sig <- ifelse(is.na(tab[["Global p-value"]]),tab[["p-value"]],tab[["Global p-value"]])
    }
    tab[["p-value"]]  <- sapply(p_sig,formatp)
    tab <- tab[,grep("Global p-value",names(tab),invert = T)]
  } else {
    raw_p <- tab[["p-value"]]
    p_sig <- suppressWarnings(stats::p.adjust(raw_p,method=p.adjust))
    tab[["p-value"]] <- sapply(p_sig,formatp)
  }
  to_bold_p <- which(as.numeric(p_sig)<.05)

  if (length(to_bold_p)>0) bold_cells <- rbind(bold_cells,
                                               matrix(cbind(to_bold_p, which(names(tab)=='p-value')),ncol=2))

  names(tab)[1] <-covTitle
  for (a in setdiff(names(att_tab),names(attributes(tab)))) attr(tab,a) <- att_tab[[a]]
  if (tableOnly){
    if (names(tab)[1]=='') names(tab)[1]<- 'Covariate'
    if (length(cap_warn)>0) message(cap_warn)
    attr(tab,"data") <- df_nm
    attr(tab,"data call") <- deparse1(argList$data)
    attr(tab, 'to_indent') <- to_indent
    attr(tab,'bold_cells') <- bold_cells
    attr(tab,'dimchk') <- dim(tab)
    return(tab)
  }
  if (returnModels) return (rtn$models)
  argL <- list(tab=tab, digits = digits,
               to_indent=to_indent,bold_cells=bold_cells,
               caption=caption,
               chunk_label=ifelse(missing(chunk_label),'NOLABELTOADD',chunk_label))
  if (!missing(fontsize)) argL[['fontsize']] <- fontsize
  do.call(outTable, argL)

}




#' Format a regression model nicely for 'Rmarkdown'
#'
#' Multivariable (or univariate) regression models are re-formatted for
#' reporting and a global p-value is added for the evaluation of factor
#' variables.
#'
#' Global p-values are likelihood ratio tests for lm, glm and polr models. For
#' lme models an attempt is made to re-fit the model using ML and if,successful
#' LRT is used to obtain a global p-value. For coxph models the model is re-run
#' without robust variances with and without each variable and a LRT is
#' presented. If unsuccessful a Wald p-value is returned. For GEE and CRR models
#' Wald global p-values are returned. For negative binomial models a deviance
#' test is used.
#'
#' If the variance inflation factor is requested (VIF=T) then a generalised VIF
#' will be calculated in the same manner as the car package.
#'
#' The number of decimals places to display the statistics can be changed with
#' digits, but this will not change the display of p-values. If more significant
#' digits are required for p-values then use tableOnly=TRUE and format as
#' desired.
#' @param model model fit
#' @param data data that model was fit on (an attempt will be made to extract
#'   this from the model)
#' @param digits number of digits to round estimates to, does not affect
#'   p-values
#' @param covTitle character with the names of the covariate (predictor) column.
#'   The default is to leave this empty for output or, for table only output to
#'   use the column name 'Covariate'.
#' @param showN boolean indicating sample sizes should be shown for each
#'   comparison, can be useful for interactions
#' @param showEvent boolean indicating if number of events should be shown. Only
#'   available for logistic.
#' @param CIwidth width for confidence intervals, defaults to 0.95
#' @param vif boolean indicating if the variance inflation factor should be
#'   included. See details
#'@param whichp string indicating whether you want to display p-values for
#'  levels within categorical data ("levels"), global p values ("global"), or
#'  both ("both"). Irrelevant for continuous predictors.
#' @param caption table caption
#' @param tableOnly boolean indicating if unformatted table should be returned
#' @param p.adjust p-adjustments to be performed. Uses the
#'  [p.adjust] function from base R
#' @param unformattedp boolean indicating if you would like the p-value to be
#'   returned unformatted (ie not rounded or prefixed with '<'). Should be used
#'   in conjuction with the digits argument.
#' @param nicenames boolean indicating if you want to replace . and _ in strings
#'   with a space
#' @param chunk_label only used if output is to Word to allow cross-referencing
#' @param fontsize PDF/HTML output only, manually set the table fontsize
#' @return A character vector of the table source code, unless tableOnly=TRUE in
#'   which case a data frame is returned
#' @export
#' @references John Fox & Georges Monette (1992) Generalized Collinearity
#'   Diagnostics, Journal of the American Statistical Association, 87:417,
#'   178-183, \doi{10.1080/01621459.1992.10475190}
#' @references  John Fox and Sanford Weisberg (2019). An {R} Companion to
#'   Applied Regression, Third Edition. Thousand Oaks CA: Sage.
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
rm_mvsum <- function(model, data, digits=getOption("reportRmd.digits",2),covTitle='',showN=TRUE,showEvent=TRUE,CIwidth=0.95, vif=FALSE,
                     whichp=c("levels","global","both"),
                     caption=NULL,tableOnly=FALSE,p.adjust='none',unformattedp=FALSE,nicenames = TRUE,chunk_label, fontsize){
  if (unformattedp) formatp <- function(x) {as.numeric(x)}
  whichp <- match.arg(whichp)

  # get the table
  tab <- mvsum(model=model,data=data,digits=digits,markup = FALSE,
               sanitize = FALSE, nicenames = FALSE,showN=showN,showEvent=showEvent,CIwidth = CIwidth,vif=vif)
  att_tab <- attributes(tab)
  if ("Global p-value" %in% names(tab)){
    tab[["Global p-value"]][which(tab[["Global p-value"]]==''|tab[["Global p-value"]]=='NA')] <-NA
  }
  to_indent <- which(!attr(tab,"varID"))
  to_bold_name <- which(attr(tab,"varID"))
  bold_cells <- arrayInd(to_bold_name, dim(tab))

  # decide which p-values to keep
  if (whichp=="levels"){
    if ("Global p-value" %in% names(tab)) tab[["Global p-value"]] <- NULL
  } else if (whichp=="global"){
    if ("Global p-value" %in% names(tab)) {
      tab$`p-value` <- sapply(1:nrow(tab), function(x) {
        ifelse(att_tab$varID[x],
               ifelse(tab[["Global p-value"]][x]!="",tab[["Global p-value"]][x],tab$`p-value`[x]),"")
      })
      tab[["Global p-value"]] <- NULL
    }
  } # if both then leave as is

  # perform p-value adjustment across variable-level p-values remove factor p-values
  if ("Global p-value" %in% names(tab)){
    if(p.adjust!='none') {
      raw_p <- ifelse(is.na(tab[["Global p-value"]]),tab[["p-value"]],tab[["Global p-value"]])
      raw_p[!att_tab$varID] <- NA
      p_sig <- suppressWarnings(stats::p.adjust(raw_p,method=p.adjust))
      message('Global p-values were adjusted according to the ',p.adjust,' method. Factor level p-values have been removed.')
      tab[["raw p-value"]]<-formatp(raw_p)
    } else{
      p_sig <- ifelse(is.na(tab[["Global p-value"]]),tab[["p-value"]],tab[["Global p-value"]])
    }
    tab[["p-value"]]  <- sapply(p_sig,formatp)
    tab <- tab[,grep("Global p-value",names(tab),invert = T)]
  } else {
    raw_p <- tab[["p-value"]]
    p_sig <- suppressWarnings(stats::p.adjust(raw_p,method=p.adjust))
    tab[["p-value"]] <- sapply(p_sig,formatp)
  }
  to_bold_p <- which(as.numeric(p_sig)<.05)

  if (length(to_bold_p)>0)  bold_cells <- rbind(bold_cells,
                                                matrix(cbind(to_bold_p, which(names(tab)=='p-value')),ncol=2))


  if (nicenames){
    attr(tab,"termnames") <- tab$Covariate
    tab$Covariate <- replaceLbl(att_tab$data, tab$Covariate)
  }
  names(tab)[1] <-covTitle
  for (a in setdiff(names(att_tab),names(attributes(tab)))) attr(tab,a) <- att_tab[[a]]
  if (tableOnly){
    if (names(tab)[1]=='') names(tab)[1]<- 'Covariate'
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
