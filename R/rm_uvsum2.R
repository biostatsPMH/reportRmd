rm_uvsum2 <- function(response, covs , data , digits=getOption("reportRmd.digits",2), covTitle='',caption=NULL,
                     tableOnly=FALSE,removeInf=FALSE,p.adjust='none',unformattedp=FALSE,
                     whichp=c("levels","global","both"),
                     chunk_label,
                     gee=FALSE,id = NULL,corstr = NULL,family = NULL,type = NULL,
                     offset=NULL,
                     strata = 1,
                     nicenames = TRUE,showN=TRUE,showEvent=TRUE,CIwidth = 0.95,
                     reflevel=NULL,returnModels=FALSE,fontsize,forceWald){

  if (missing(data)) stop('data is a required argument')
  if (missing(covs)) stop('covs is a required argument') else covs <- unique(covs)
  if (missing(response)) stop('response is a required argument')
  if (length(response)>2) stop('The response must be a single outcome for linear, logistic and ordinal models or must specify the time and event status variables for survival models.')
  if (!inherits(data,'data.frame')) stop('data must be supplied as a data frame.')
  if (!inherits(covs,'character')) stop('covs must be supplied as a character vector or string indicating variables in data')
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
  rtn <- uvsum2(response,covs,data,digits=digits,markup = FALSE,sanitize=FALSE,
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

uvsum2 <- function (response, covs, data, digits=getOption("reportRmd.digits",2),id = NULL, corstr = NULL, family = NULL,
                   type = NULL, offset=NULL, gee=FALSE,strata = 1, nicenames = TRUE,
                   showN = TRUE, showEvent = TRUE, CIwidth = 0.95, reflevel=NULL,returnModels=FALSE,forceWald)
{
  if (missing(forceWald)) forceWald = getOption("reportRmd.forceWald",FALSE)
  if (!sanitize)  sanitizestr <- identity
  if (!nicenames) nicename <- identity
  if (inherits(data[[response[1]]],"character")) data[[response[1]]] <- factor(data[[response[1]]])
  if (!inherits(strata,"numeric")) {
    strataVar = strata
    strata <- sapply(strata, function(stra) {
      paste("strata(", stra, ")", sep = "")
    })
  }
  else {
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
  }
  else {
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

  if (forceWald) confint <- confint.default
  if (strata != "" & type != "coxph") {
    stop("strata can only be used with coxph")
  }
  if (!is.null(id)){
    if (! (gee | type =='coxph')) {
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

  if (returnModels) modelList <- NULL
  modelList <- c()
  for (cov in covs) {
    cov_model <- autoreg(response, data, cov, id, strata, family, offset, corstr)
    modelList <- c(modelList, cov_model)
  }
  summaryList <- NULL
  for (model in modelList) {
    cov_summary <- modelsum(model, digits, CIwidth, vif, whichp, ...)
    summaryList <- c(summaryList, cov_summary)
  }
  summaryList <- bindrows(summaryList)




  table <- lapply(out, function(x) {
    return(x[[1]])
  })
  varID <- do.call("c",lapply(table,function(x){
    return(stats::setNames(c(TRUE,rep(FALSE,nrow(x)-1)),x[,1]))
  }))
  table <- do.call("rbind", lapply(table, data.frame,
                                   stringsAsFactors = FALSE))
  colName <- c("Covariate", sanitizestr(beta),
               "p-value", "Global p-value")
  if (showN) colName <- c(colName,"N")
  if (showEvent & type == "logistic") colName <- c(colName,"Event")
  colnames(table) <- colName
  table[,"Global p-value"] <- ifelse(table[,'p-value']=='',table[,"Global p-value"],'')
  if (all(table[,"Global p-value"]=='')) table <- table[, -which(colnames(table)=="Global p-value")]
  colnames(table) <- sapply(colnames(table), lbld)
  attr(table,"varID") <- varID
  if (returnModels) return(list(table,models=modelList)) else return(table)
}
