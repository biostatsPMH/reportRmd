devtools::load_all()


varNames <- data.frame(var=names(pembrolizumab),
                       label=c('Patient ID','Age at study entry','Patient Sex','Study Cohort','Target lesion size at baseline','PD L1 percent','log of tumour size','Baseline ctDNA','Did ctDNA increase or decrease from baseline to cycle 3','Objective Response','Clinical Beneficial Response','Overall survival status', 'Overall survival time in months','Progression free survival status','Progression free survival time in months'))


ctDNA_names <- data.frame(var=names(ctDNA),
                          label=c('Patient ID',
                                  'Study Cohort',
                                  'Change in ctDNA since baseline',
                                  'Number of weeks on treatment',
                                  'Percentage change in tumour measurement'))


getOption("reportRmd.variable_info")

clearVariableLabels <- function(){
  options(reportRmd.v_info = NULL)
}

setVariableLabels <- function(default,...){
  pars <- as.list(match.call()[-1])
#  return(pars)
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
  for (nm in setdiff(names(pars),'default')){
    if (!exists(as.character(nm))) stop(paste(nm,' does not exist. Please specify a valid data frame. \nSee examples.'))
    if (!inherits(get0(nm),'data.frame')) stop(paste(nm,'is not a data frame. Please specify the variable labels in the format data=labels.\nSee examples.'))
  }
  if (names(pars)[1]!='default' & !('default' %in% names(old_pars))) message('There is no default variable label table specified. Variable labels will only be matched to the specified table(s).')
  new_pars <- c(old_pars[setdiff(names(old_pars),names(pars))], pars)
  options(reportRmd.v_info = new_pars)

}

pars <- setVariableLabels(varNames,ctDNA=ctDNA_names)

getVariableLabels <-function(){
  vl <- getOption("reportRmd.v_info")
  vL <- data.frame(data.frame = names(vl),
                   variable.names = as.character(vl))

}
test_fun <- function(data){
#  return(match.call())
  getVL(match.call())
}

getVL <- function(call){
   pars <- as.list(call[-1])
   dn <- as.character(pars$data)
   var_info <- getOption("reportRmd.v_info",NULL)
   vl <- NULL
   if (!is.null(var_info)){
     if ('default' %in% names(var_info)) vl <- var_info[['default']]
     if (dn %in% names(var_info)) vl <- var_info[[dn]]
   }
    vL <- data.frame(var=names(get0(dn)))
   if (!is.null(vl)){
     vl1 <- eval(vl)
     names(vl1) <- c('var','lbl')
    vL <- merge(vL,vl1,all.x = T)
    } else vL$lbl <- NA
    vL$lbl <- ifelse(is.na(vL$lbl),vL$var,vL$lbl)
  return(vL)
}

extractLabels <- function(data){
  if (!inherits(data,'data.frame')) stop('data must be a data frame with labelled variables to extract.')
  lblv <- lapply(data, function(x){
    typ <- grep('label',class(x),value=TRUE)
    if (length(typ)==0) return(NA)
    if (typ=="haven_labelled") lbl <- attr(x,"label")
    if (typ=="labelled") lbl <- attr(x,"label")
    return(lbl)
  })
  if (all(is.na(unlist(lblv)))) stop('No labelled variables detected in data.\nCurrently this function can extract labels from the haven')
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
extractLabels(mtcars)

dn<-'mtcars_names'
get0(dn)
mtcars_names
setVariableLabels(varNames)
test_fun(ctDNA)

setVariableLabels(varNames,ctDNA=ctDNA_names)
 test_fun(ctDNA)

setVariableLabels(ctDNA=ctDNA_names)
test_fun(pembrolizumab)
getOption("reportRmd.v_info")

call <- test_fun(pembrolizumab)


clearVariableLabels()
setVariableLabels(varNames)
getOption("reportRmd.v_info")
setVariableLabels(ctDNA=ctDNA_names)
getOption("reportRmd.v_info")


setVariableLabels(ctDNA=ctDNA_names)
test_fun(pembrolizumab)
getOption("reportRmd.v_info")

library(haven)
df <- read_sav("C:/Users/lavery/OneDrive - UHN/Jones/CancerRelatedFatigue/ARCC breast_FINAL.sav")

class(df)
class(df[[1]])
df2 <- df[1:5,]
is_labelled(df2)
class(df2$MaritalStatus)

data=df2
getLabels <- function(data){
  if (!inherits(data,'data.frame')) stop('data must be a data frame with labelled variables to extract.')
  lblv <- lapply(data, function(x){
    typ <- grep('label',class(x),value=TRUE)
    if (length(typ)==0) return(NA)
    if (typ=="haven_labelled") lbl <- attr(x,"label")

  })
  lbldf <- data.frame(var=names(unlist(lblv)),lbl=unlist(lblv))
  rownames(lbldf) <- NULL
  lbldf$lbl <- ifelse(is.na(lbldf$lbl),lbldf$var,lbldf$lbl)

}

library(expss)
data(mtcars)
mtcars = apply_labels(mtcars,
                      mpg = "Miles/(US) gallon",
                      cyl = "Number of cylinders",
                      disp = "Displacement (cu.in.)",
                      hp = "Gross horsepower",
                      drat = "Rear axle ratio",
                      wt = "Weight (1000 lbs)",
                      qsec = "1/4 mile time",
                      vs = "Engine",
                      vs = c("V-engine" = 0,
                             "Straight engine" = 1),
                      am = "Transmission",
                      am = c("Automatic" = 0,
                             "Manual"=1),
                      gear = "Number of forward gears",
                      carb = "Number of carburetors"
)

