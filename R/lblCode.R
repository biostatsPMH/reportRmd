rRmd.env <- new.env(parent = emptyenv())
rRmd.env$varInfo <- NULL

#' Clear all variable labels associated with data frames
#'
#' @export
#' @examples
#' # Clear all variable data frames and check
#' clearVariableLabels()
#' getVariableLabels()
clearVariableLabels <- function(){
 rRmd.env$varInfo <- NULL
}

#' Set variable labels
#'
#' Assigns default or data frame specific tables containing variable labels
#'
#' The variable labels need to be stored in two-column data frames. The first
#' column must contain the variable name and the second column the variable
#' label. The column names are not used.
#'
#' If a data frame specific table of labels is provided it will be used for that
#' data frame only and labels contained in the default table will not be used.
#' If no labels are supplied for a variable then the variable name will be used
#' with '_' replaced by a space.
#' @param default a data frame containing two columns, variable names and
#'   variable labels to be used as default variable labels.
#' @param ... pairs of data.frame, variable label pairs in the format
#'   data=labels to provide variable labels specific to a data frame
#' @export
#' @examples
#' # create a data frame contain variable labels and set as default
#' varLabels <- data.frame(vars=c('Var1','Var2','Var3'),lbls=c('Age','Sex','Score'))
#' setVariableLabels(varLabels)
#'
#' # add a separate table of data labels for the ctDNA data
#' data(ctDNA)
#' ctDNA_names <- data.frame(var=names(ctDNA),
#' label=c('Patient ID',
#'         'Study Cohort',
#'         'Change in ctDNA since baseline',
#'         'Number of weeks on treatment',
#'         'Percentage change in tumour measurement'))
#' setVariableLabels(ctDNA=ctDNA_names)
#'
#' # Check that both the default and ctDNA label tables are set
#' getVariableLabels()
setVariableLabels <- function(default,...){
  pars <- as.list(match.call()[-1])
  if (length(pars)==0) stop('No variable labels specified.')
  old_pars <- rRmd.env$varInfo
  for (x in names(pars)){
    if (x %in% names(old_pars))   message(paste(x,'variable names will be updated.' ))
  }
  #check that all labels exist as data frames and that if more than one table is supplied it is named
  if (any(names(pars)=="")) stop("Only the default table may be unnamed. \nTo supply data frame specific labels please provide in the format data=labels.\nSee examples.")
  for (df in pars){
    if (!exists(as.character(df))) stop(paste(df,'does not exist. Please specify a valid data frame for labels. \nSee examples.'))
    if (!inherits(eval(df),'data.frame')) stop(paste(df,'is not a data frame. Variable labels must be specified in a data frame with two columns: variable names and variable labels.'))
    if (ncol(eval(df))!=2) stop(paste(df,'must be a data frame with two columns: variable names and variable labels.'))
  }
  rm(df)
  for (nm in setdiff(names(pars),'default')){
    if (!exists(as.character(nm))) stop(paste(nm,' does not exist. Please specify a valid data frame. \nSee examples.'))
    if (!inherits(get0(nm),'data.frame')) stop(paste(nm,'is not a data frame. Please specify the variable labels in the format data=labels.\nSee examples.'))
  }
  if (names(pars)[1]!='default' & !('default' %in% names(old_pars))) message('There is no default variable label table specified. \nVariable labels will only be matched to the specified table(s).')

  new_pars <- c(old_pars[setdiff(names(old_pars),names(pars))], pars)
  rRmd.env$varInfo <- new_pars
}


#' Show all variable label data frames
#'
#' @export
#' @examples
#' # add a separate table of data labels for the ctDNA data
#' data(ctDNA)
#' ctDNA_names <- data.frame(var=names(ctDNA),
#' label=c('Patient ID',
#'         'Study Cohort',
#'         'Change in ctDNA since baseline',
#'         'Number of weeks on treatment',
#'         'Percentage change in tumour measurement'))
#' setVariableLabels(ctDNA=ctDNA_names)
#'
#' # Check that both the default and ctDNA label tables are set
#' getVariableLabels()
getVariableLabels <-function(){
  vl <- rRmd.env$varInfo
  vL <- data.frame(data.frame = names(vl),
                   variable.names = as.character(vl))
  return(vL)
}


#' Extract variable labels from labelled data and set
#'
#' Retrieve labels set using the haven, expss or sjlabelled and set the variable
#' label table.
#'
#' By default the variable names will be associated with the data frame. To make
#' the variable labels the default table specify default=TRUE.
#'
#' If no labels are supplied for a variable then the variable name will be used
#' with '_' replaced by a space.
#' @param data data frame with labelled variables to be extracted.
#' @param default logical, should the variable labels be used as the default
#'   labels for all data frames, default is FALSE.
#' @export
#' @examples
#' # example code (not run)
#' # library(sjlabelled)
#' # data(efc)
#' # extractLabels(efc)
extractLabels <- function(data,default=FALSE){
  argList <- match.call()[-1]
  if (!inherits(data,'data.frame')) stop('data must be a data frame with labelled variables to extract.')
  lblv <- lapply(data, function(x){
    typ <- grep('label',class(x),value=TRUE)
    if (length(typ)==0) {
      lbl <- attr(x,"label")
    } else if (typ %in% c("labelled","haven_labelled")) lbl <- attr(x,"label")
    return(ifelse(is.null(lbl),NA,lbl))
  })
  if (all(is.na(unlist(lblv)))) stop('No labelled variables detected in data.\nCurrently this function can extract labels created by the haven, expss and sjlabelled')
  lbldf <- data.frame(var=names(unlist(lblv)),lbl=unlist(lblv))
  rownames(lbldf) <- NULL
  lbldf$lbl <- ifelse(is.na(lbldf$lbl),lbldf$var,lbldf$lbl)
  dn <- paste0(as.character(argList$data),  "_names")
  cl <- paste(dn,"<<-lbldf")
  eval(parse(text = cl))
  cl <- paste("setVariableLabels(",ifelse(default,"default",as.character(match.call()[-1])),"=",dn,")")
  eval(parse(text = cl))
  message(paste(dn,'created and added to list of variable tables'))
  return(get0(dn))
}

getVL <- function(data_name){
  if (!exists(data_name)) stop(paste(data_name),'does not exist.')
  var_info <- rRmd.env$varInfo
  vl <- NULL
  if (!is.null(var_info)){
    if ('default' %in% names(var_info)) vl <- var_info[['default']]
    if (data_name %in% names(var_info)) vl <- var_info[[data_name]]
  }
  if (!is.null(vl)){
    if (!vldTbl(vl,data_name)) return(NULL)
    vL <- eval(vl)
    if (ncol(vL)!=2)    message(paste(as.character(vl),'has been altered. The first two columns are assumed to be variable names and labels.'))
    names(vL)[1:2] <- c('var','lbl')
    return(vL)
  } else return(NULL)
}


vldTbl <- function(vL, dn){
  if (!exists(dn)) {
    message(paste(dn," does not exist"))
    return(FALSE)
  }
  if (!inherits(dn,'character')) {
    message('data table must be specified as a character string.')
    return(FALSE)
  }
  if (!inherits(vL,'name')) {
    message('variable names must be specified as a data.frame.')
    return(FALSE)
  }
  if (!exists(vL)){
    message(paste(as.character(vL)," does not exist, assigning default variable labels"))
    return(FALSE)
  }
  df <- get0(dn)
  if (!inherits(df,'data.frame')) {
    message(paste(dn,'is not a data frame'))
    return(FALSE)
  }
  if (!inherits(eval(vL),'data.frame')) {
    message(paste(as.character(vL),'is not a data frame, assigning default variable labels'))
    return(FALSE)
  }
  return(TRUE)
}

# return variable labels associated with variables
replaceLbl <- function(dn,cv){
  if(length(dn)>1) dn <- dn[2]

  if (!inherits(dn,'character')) stop('data table must be specified as a character string.')
  if (!inherits(cv,'character')) stop('variable name must be specified as a character string.')
  lbl <- getVL(dn)
  vl <- data.frame(var=cv,ord=1:length(cv))
  if (!is.null(lbl)){
    cvnew <- merge(vl,lbl,all.x=T)
    cvnew <- cvnew[order(cvnew$ord),]
  } else {
    cvnew <- vl
    cvnew$lbl <- NA
  }
  cvnew$lbl <- ifelse(is.na(cvnew$lbl),nicename(cvnew$var),cvnew$lbl)
  return(cvnew$lbl)
}
