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
#'
#'tidyselect can only be used for response and covs variables. Additional
#'arguments must be passed in using characters
#'
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
#'@param p.adjust p-adjustments to be performed (Global p-values only)
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
#'@param forceWald boolean indicating if Wald confidence intervals should be
#'  used instead of profile likelihood. This is not recommended, but can speed
#'  up computations. To use throughout a document use
#'  options(reportRmd.forceWald=TRUE)
#'@seealso
#'\code{\link{uvsum}},\code{\link{lm}},\code{\link{glm}},\code{\link{crr}},
#'\code{\link{coxph}}, \code{\link{lme}},\code{\link{geeglm}},\code{\link{polr}}
#'@return A character vector of the table source code, unless tableOnly=TRUE in
#'  which case a data frame is returned
#'@export
#' @examples
#' # Examples are for demonstration and are not meaningful
#' # Coxph model with 90% CI
#' data("pembrolizumab")
#' rm_uvsum2(response = c('os_time','os_status'),
#' covs=c('age','sex','baseline_ctdna','l_size','change_ctdna_group'),
#' data=pembrolizumab,CIwidth=.9)
#'
#' # Linear model with default 95% CI
#' rm_uvsum2(response = 'baseline_ctdna',
#' covs=c('age','sex','l_size','pdl1','tmb'),
#' data=pembrolizumab)
#'
#' # Logistic model with default 95% CI
#' rm_uvsum2(response = 'os_status',
#' covs=c('age','sex','l_size','pdl1','tmb'),
#' data=pembrolizumab,family = binomial)

#' # Poisson models returned as model list
#' mList <- rm_uvsum2(response = 'baseline_ctdna',
#' covs=c('age','sex','l_size','pdl1','tmb'),
#' data=pembrolizumab, returnModels=TRUE)
#' #'
#' # GEE on correlated outcomes
#' data("ctDNA")
#' rm_uvsum2(response = 'size_change',
#' covs=c('time','ctdna_status'),
#' gee=TRUE,
#' id='id', corstr="exchangeable",
#' family=gaussian("identity"),
#' data=ctDNA,showN=TRUE)
#'
#' # Using tidyselect
#' pembrolizumab |> rm_uvsum2(response = sex,
#' covs = c(age, cohort))
rm_uvsum2 <- function(response, covs , data , digits=getOption("reportRmd.digits",2), covTitle='',caption=NULL,
                      tableOnly=FALSE,removeInf=FALSE,p.adjust='none',unformattedp=FALSE,
                      whichp=c("levels","global","both"),
                      chunk_label,
                      gee=FALSE,id = NULL,corstr = NULL,family = NULL,type = NULL,
                      offset=NULL,
                      strata = 1,
                      nicenames = TRUE,showN=TRUE,showEvent=TRUE,CIwidth = 0.95,
                      reflevel=NULL,returnModels=FALSE,fontsize,forceWald){
  response_var <- tidyselect::eval_select(expr = enquo(response), data = data[unique(names(data))],
                                          allow_rename = FALSE)
  response_var <- names(response_var)

  x_vars <- tidyselect::eval_select(expr = enquo(covs), data = data[unique(names(data))],
                                    allow_rename = FALSE)
  x_vars <- names(x_vars)
  if (missing(data)) stop('data is a required argument')
  if (missing(covs)) stop('covs is a required argument') else covs <- unique(x_vars)
  if (missing(response)) stop('response is a required argument')
  empty <- NULL
  if ("" %in% c(strata, type, offset, id)) {
    args <- list(strata = strata, type = type, offset = offset, id = id)
    empty <- names(args)[which(args == "")]
    for (var in empty) {
      fun <- get(match.call()[[1]])
      assign(var, formals(fun)[[var]])
    }
    warning(paste0("empty string arguments "), paste(empty, collapse = ", "), " will be ignored")
  }

  response <- response_var
  covs <- x_vars
  if (length(response)>2) stop('The response must be a single outcome for linear, logistic and ordinal models or must specify the time and event status variables for survival models.')
  if (!inherits(data,'data.frame')) stop('data must be supplied as a data frame.')
  # if (!inherits(covs,'character')) stop('covs must be supplied as a character vector or string indicating variables in data')
  missing_vars = na.omit(setdiff(c(response, covs,id,ifelse(strata==1,NA,strata)), names(data)))
  if (length(missing_vars) > 0) stop(paste("These variables are not in the data:\n",
                                           paste0(missing_vars,collapse=csep())))
  if (strata==1) nm <- c(response,covs) else nm <- c(strata,response,covs)
  if (!all(names(data[,nm])==names(data.frame(data[,nm])))) stop('Non-standard variable names detected.\n Try converting data with new_data <- data.frame(data) \n then use new variable names in rm_uvsum.' )
  if (missing(forceWald)) forceWald = getOption("reportRmd.forceWald",FALSE)

  argList <- as.list(match.call()[-1])
  if (!all(sapply(argList$covs[-1], is.character))) {
    argList$covs <- x_vars
  }
  if (!is.character(argList$response)) {
    argList$response <- response_var
  }
  if ("tableOnly" %in% names(argList)) {
    argList[["tableOnly"]] <- NULL
  }
  if (!is.null(empty)) {
    for (var in empty) {
      if (is.null(eval(as.name(var)))) {
        argList[var] <- list(NULL)
      }
      else {
        argList[[var]] <- eval(as.name(var))
      }
    }
  }

  nomiss <- nrow(na.omit(data[,response,drop=FALSE]))
  if (nrow(data)!=nomiss) warning(paste("Cases with missing response data have been removed.\n",
                                          nrow(data)-nomiss,"case(s) removed."))
  ## TO DO - split out testing for this function,
  ## and data changes for the uvsum2 function
  ## need to do this, because when we run the models
  ## we and summarise the data from the same function
  ## that contains the changed data
  ## -------------------------------------------------
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

  # get the table
  tab <- do.call(uvsum2,argList)
  # rtn <- uvsum2(response,covs,data,digits=digits,
  #               gee=gee,id = id, offset=offset,
  #               corstr = corstr,family = family,type = type,strata = strata,
  #               nicenames = FALSE,showN = showN,showEvent = showEvent,
  #               CIwidth = CIwidth,reflevel=reflevel,returnModels=returnModels,forceWald = forceWald)

  # If user specifies return models, don't format a table, just return a list of models
  to_indent <- attr(tab, "to_indent")
  bold_cells <- attr(tab, "bold_cells")
  if (returnModels) return (tab$models)


  att_tab <- attributes(tab)


  # Adjust / format pvalues based on unformattedp, p.adjust

  # Add variable labels
  # may be easier to use extract_labels on data frame

  # need to get bold/indent attributes

  # formatting pval column
  method <- p.adjust
  tab[["p-value"]] <- p.adjust(tab[["p-value"]], method = method)
  if (!unformattedp) {
    tab[["p-value"]] <- formatp(tab[["p-value"]])
  }

  # changing UB to Inf, LB to 0
  tab[, 2] <- sapply(tab[, 2], process_ci)

  p_col <- (which(names(tab) == "p-value"))
  bold_cells <- rbind(bold_cells, cbind(which(as.numeric(gsub("[^0-9\\.]", "", tab[["p-value"]])) < 0.05), rep(p_col, length(which(as.numeric(gsub("[^0-9\\.]", "", tab[["p-value"]])) < 0.05)))))
  if (nrow(bold_cells) < 1) {
    bold_cells <- NULL
  }

  names(tab)[1] <- covTitle
  lbl <- tab[, 1]
  if (nicenames) {
    tab[, 1] <- replaceLbl(argList$data, lbl)
  }
  argL <- list(tab=tab, digits = digits,
               to_indent=to_indent,bold_cells=bold_cells
  )
  # Add attributes if returning a table
  for (a in setdiff(names(att_tab),names(attributes(tab)))) attr(tab,a) <- att_tab[[a]]
  if (tableOnly){
    if (names(tab)[1]=='') names(tab)[1]<- 'Covariate'
    attr(tab,"data call") <- deparse1(argList$data)
    attr(tab, 'to_indent') <- to_indent
    attr(tab,'bold_cells') <- bold_cells
    attr(tab,'dimchk') <- dim(tab)
    return(tab)
  }

  do.call(outTable, argL)

}

uvsum2 <- function (response, covs, data, digits=getOption("reportRmd.digits",2),id = NULL, corstr = NULL, family = NULL,
                    type = NULL, offset=NULL, gee=FALSE,strata = 1, nicenames = TRUE,
                    showN = TRUE, showEvent = TRUE, CIwidth = 0.95, reflevel=NULL,returnModels=FALSE,forceWald, whichp = "level")
{
  argList <- as.list(match.call()[-1])
  if (missing(forceWald)) forceWald = getOption("reportRmd.forceWald",FALSE)
  if (inherits(data[[response[1]]],"character")) data[[response[1]]] <- factor(data[[response[1]]])
  if (!inherits(strata,"numeric")) {
    strataVar = strata
    strata <- sapply(strata, function(stra) {
      paste("strata(", stra, ")", sep = "")
    })
  }  else {
    strataVar <- ""
    strata <- ""
  }
  if (length(response)==1) {
    if (sum(is.na(data[[response]]))>0) message(paste(sum(is.na(data[[response]])),"observations with missing outcome removed."))
    data <- subset(data,!is.na(data[[response]]))
  } else {
    if (sum(is.na(data[[response[1]]])|is.na(data[[response[2]]]))>0) message(paste(sum(is.na(data[[response[1]]])|is.na(data[[response[2]]])),"observations with missing outcome removed."))
    data <- subset(data,!(is.na(data[[response[1]]])|is.na(data[[response[2]]])))
  }
  # Set family if user specifies type
  if (!is.null(type)) {
    if (length(response)==1 & (type %in% c('coxph','crr')))
      stop('Please specify two variables in the response for survival models. \nExample: response=c("time","status")')
    if (length(response)==2 & !(type %in% c('coxph','crr')))
      stop('Response can only be of length one for non-survival models.')
    if (type == "logistic") {
      if (is.null(family)) family='binomial'
    }
    else if (type == "poisson") {
      if (all(data[[response]]==as.integer(data[[response]]))){
        data[[response]]=as.integer(data[[response]])
      }
      else {
        stop('Poisson regression requires an integer response.')
      }
      if (is.null(family)) family='poisson'
    }
    else if (type == "negbin") {
      if (all(data[[response]]==as.integer(data[[response]]))){
        data[[response]]=as.integer(data[[response]])
      }
      else {
        stop('Negative binomial regression requires an integer response.')
      }
      if (!is.null(family)) message('For negative binomial regression currently only the log link is implemented.')
    }
    else if (type == "linear" | type == "boxcox") {
      if (is.null(family)) family='gaussian'
    }
    else if (type == "ordinal") {
      if (!inherits(data[[response[1]]],c("factor","ordered"))) {
        warning("Response variable is not a factor, will be converted to an ordered factor")
        data[[response]] <- factor(data[[response]],
                                   ordered = T)
      }
      if (!is.null(reflevel)) {
        data[[response]] <- stats::relevel(data[[response]],
                                           ref = reflevel)
      }
    }
    else {
      stop("type must be either coxph, logistic, linear, poisson, negbin, boxcox, crr, ordinal (or NULL)")
    }
  }    else {
    if (length(response) == 2) {
      # Check that responses are numeric
      for (i in 1:2) if (!is.numeric(data[[response[i]]])) stop('Both response variables must be numeric')
      if (length(unique(na.omit(data[[response[2]]]))) < 3) {
        type <- "coxph"
      }
      else {
        type <- "crr"
      }
    } else if (length(unique(na.omit(data[[response]]))) == 2) {
      type <- "logistic"
      family="binomial"
    } else if (inherits(data[[response[1]]],"ordered")) {
      type <- "ordinal"
      if (!is.null(reflevel)) {
        data[[response]] <- stats::relevel(data[[response]],
                                           ref = reflevel)
      }
    } else if (inherits(data[[response[1]]],"integer")) {
      type <- "poisson"
      family="poisson"
    } else {
      if (!inherits(data[[response[1]]],"numeric")) stop('Response variable must be numeric')
      type <- "linear"
      beta <- "Estimate"
      family='gaussian'
    }
  }
  # Do some more model checking --------
  if (forceWald) confint <- confint.default
  if (strata != "" & type != "coxph") {
    stop("strata can only be used with coxph")
  }
  if (!is.null(id)){
    if (!(gee | type =='coxph')) {
      warning('id argument will be ignored. This is used only for survival strata or clustering in GEE. To run a GEE model set gee=TRUE.')
    }
  }
  if (!is.null(offset) & !(type %in% c('poisson','negbin'))) {
    warning('Offset terms only used for Poisson and negative binomial regression.\nOffset term will be ignored.')
  }
  if (!is.null(corstr)){
    if (! (gee | type =='coxph')) {
      warning('id argument will be ignored. This is used only for survival strata or clustering in GEE. To run a GEE model set gee=TRUE.')
    }
  }
  if (!is.null(offset)){
    ovars <- unlist(strsplit(offset,"[^a-zA-Z_]"))
    if(length(intersect(names(data),ovars))==0){
      stop(paste('Variable names in the offset term contains special characters. \nPlease remove special characters, except "_" from the variable name and re-fit.\n',
                 'offset =',offset))
    } else ovars <- intersect(names(data),ovars)
  } else ovars <- NULL
  if (gee){
    if (!type %in% c('linear','logistic','poisson')) stop('GEE models currently only implemented for Poisson, logistic or linear regression.')
    if (is.null(id)) stop('The id argument must be set for gee models to indicate clusters.')
    if (is.null(corstr)) stop ('You must provide correlation structure (i.e. corstr="independence") for GEE models.')
  }
  # Assign the class to response --------
  data("uvmodels", envir=environment())

  tp_merge <- merge(data.frame(type=type,family=ifelse(is.null(family),NA,family),gee=gee),uvmodels,all.x = T)
  if (nrow(tp_merge)>1) stop("Can not detect regression type, try specifying family")
  if (is.na(tp_merge$autoreg_class) | tp_merge$autoreg_class=="raise_error"){
    stop(paste("Can not detect valid regression type for the combination of type =",type,
               "family =",family, "and gee=",gee))
  }
  class(response) <- tp_merge$autoreg_class
  if (returnModels) modelList <- NULL
  modelList <- NULL
  for (cov in covs) {
    modelList[[cov]] <- autoreg(response, data, cov, id, strata, family, offset, corstr)
  }

  summaryList <- NULL
  summaryList <- lapply(modelList,m_summary,digits= digits, CIwidth=CIwidth, vif = FALSE,whichp="level", for_plot = FALSE)
  summaryList <- dplyr::bind_rows(summaryList)

  # Other stuff here!

  if (!showN) {
    summaryList[["N"]] <- NULL
  }
  if (!showEvent) {
    summaryList[["Event"]] <- NULL
  }
  to_indent <- which(!(summaryList[["Variable"]] %in% covs))
  bold_cells <- cbind(which(summaryList[["Variable"]] %in% covs), rep(1, length(which(summaryList[["Variable"]] %in% covs))))
  attr(summaryList, "to_indent") <- to_indent
  attr(summaryList, "bold_cells") <- bold_cells
  if (returnModels) return(list(summaryList,models=modelList)) else return(summaryList)
}
