
clearVariableLabels <- function(){
  options(reportRmd.v_info = NULL)
}

setVariableLabels <- function(default,...){
  pars <- as.list(match.call()[-1])
  if (length(pars)==0) stop('No variable labels specified.')
  old_pars <- getOption("reportRmd.v_info")
  for (x in names(pars)){
    if (x %in% names(old_pars))   message(paste(x,'variable names will be updated.' ))
  }
  #check that all labels exist at data frames and that if more than one table is supplied it is named
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
  if (names(pars)[1]!='default' & !('default' %in% names(old_pars))) message('There is no default variable label table specified. Variable labels will only be matched to the specified table(s).')
  new_pars <- c(old_pars[setdiff(names(old_pars),names(pars))], pars)
  options(reportRmd.v_info = new_pars)

}


getVariableLabels <-function(){
  vl <- getOption("reportRmd.v_info")
  vL <- data.frame(data.frame = names(vl),
                   variable.names = as.character(vl))
  return(vL)
}

test_fun <- function(data){
#  return(match.call())
  pars <- as.list(match.call()[-1])
  dn <- as.character(pars$data)
  getVL(dn)
}

# retrieve labels for a data frame
getVL <- function(data_name){
  dn <- data_name
   var_info <- getOption("reportRmd.v_info",NULL)
   vl <- NULL
   if (!is.null(var_info)){
     if ('default' %in% names(var_info)) vl <- var_info[['default']]
     if (dn %in% names(var_info)) vl <- var_info[[dn]]
   }
   if (!is.null(vl)){
     vL <- eval(vl)
     names(vL) <- c('var','lbl')
      return(vL)
    } else return(NULL)
}


extractLabels <- function(data){
  if (!inherits(data,'data.frame')) stop('data must be a data frame with labelled variables to extract.')
  lblv <- lapply(data, function(x){
    typ <- grep('label',class(x),value=TRUE)
    if (length(typ)==0) {
      lbl <- attr(x,"label")
    } else if (typ %in% c("labelled","haven_labelled")) lbl <- attr(x,"label")
    return(ifelse(is.null(lbl),NA,lbl))
  })
  if (all(is.na(unlist(lblv)))) stop('No labelled variables detected in data.\nCurrently this function can extract labels created by the haven, expss')
  lbldf <- data.frame(var=names(unlist(lblv)),lbl=unlist(lblv))
  rownames(lbldf) <- NULL
  lbldf$lbl <- ifelse(is.na(lbldf$lbl),lbldf$var,lbldf$lbl)
  dn <- paste0(as.character(match.call()[-1]),  "_names")
  eval(parse(text = paste(dn,"<<-lbldf")))
  message(paste(dn,'created and added to list of variable tables'))
  cl <- paste("setVariableLabels(",as.character(match.call()[-1]),"=",dn,")")
  eval(parse(text = cl))
  return(get0(dn))
}

replaceLbl <- function(dn,cv){
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
