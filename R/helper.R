xcn <- function(v){
  sapply(v, function(x){
    colHead <- toupper(x)
    if (nchar(colHead)>1){
      l1 = substr(colHead,1,1)
      l2 = substr(colHead,2,2)
      rtn <- 26*which(LETTERS==l1)+which(LETTERS==l2)
    } else {
      rtn <- which(LETTERS==colHead)
    }
  })}

is.error <- function(x) inherits(x, "try-error")
csep<-function(){return(", ")}

.negloglik.boxcox <- function (lambda.val, data, xmat, lik.method = "ML")
{
  if (length(lambda.val) == 2) {
    data <- data + lambda.val[2]
    lambda <- lambda.val[1]
  }
  else lambda <- lambda.val
  lambda <- unname(lambda)
  n <- length(data)
  beta.size <- ncol(xmat)
  if (isTRUE(all.equal(unname(lambda), 0)))
    yt <- log(data)
  else yt <- ((data^lambda) - 1)/lambda
  beta <- solve(crossprod(xmat), crossprod(xmat, yt))
  ss <- sum((drop(yt) - drop(xmat %*% beta))^2)
  if (lik.method == "ML")
    neglik <- (n/2) * log(ss) - ((lambda - 1) * sum(log(data)))
  if (lik.method == "RML") {
    xx <- crossprod(xmat)
    if (length(as.vector(xx)) == 1)
      choldet <- 0.5 * log(xx)
    else choldet <- sum(log(diag(chol(xx))))
    neglik <- ((n - beta.size)/2) * log(ss) + choldet - ((lambda -
                                                            1) * sum(log(data)))
  }
  if (mode(neglik) != "numeric")
    neglik <- Inf
  return(drop(neglik))
}



#' Round retaining digits
#'
#' Round retaining digits
#'@param x a vector
#'@param digits numeric
#'@keywords helper
niceNum <- function(x,digits=2){
  rndx = sapply(x, function(x) {
    if(is.na(x)) return(x)
    if(is.null(x)) return(x)
    format(round(as.numeric(x),digits),nsmall=digits)})
  return(gsub(" ","",rndx))
}



#' Paste with parentheses
#'
#' Paste with parentheses
#'
#'@param x a vector
#'@keywords helper
pstprn<-function(x){paste(x[1]," (",paste(x[-1],collapse=csep()),")",sep="")}

#' Round and paste with parentheses
#'
#' Round and paste with parentheses
#'
#' @param x a numeric vector
#' @param y integer corresponding to the number of digits to round by
#'@keywords helper
psthr<- function (x, y = 2)
{
  x <- sapply(x, function(x) {
    ifelse(abs(x) < 0.01 | abs(x) > 1000, format(x, scientific = TRUE,digits = y), format(round(x, y),nsmall = y))
  })
  pstprn(x)
}

covnm<-function(betanames,call){
  sapply(betanames,function(betaname){

    # indx<-which(sapply(call,function(cov){charmatch(cov,betaname)})==1)
    indx=which(sapply(call,function(cov)grepl(cov,betaname,fixed=TRUE))) ## changed on Feb 21, 2019 for checkings
    if(length(indx)==1) return(call[indx])
    #If one  facorname is a subset of another
    indx2<-which.max(sapply(call[indx],nchar))
    if(length(indx2)==1) return(call[indx[indx2]])
    indx3<-which(sapply(call[indx2],function(c){substr(betaname,1,nchar(c))==c}))
    if(length(indx3)==1)  return(call[indx[indx2[indx3]]])
  })
}

alleql<-function(x,y){
  !any((x==y)==F)
}


betaindx<-function(x){
  i=1
  out<-1
  result<-NULL
  while(TRUE){
    if(i+1>length(x)){
      result<-c(result,list(out))
      return(result)
    }
    else if(alleql(x[[i+1]],x[[i]])){
      out<-c(out,i+1)
    }
    else{
      result<-c(result,list(out))
      out<-i+1
    }
    i=i+1
  }
}


#' Capitalize a string
#'
#' Capitalize a string
#'
#' @param x string
#' @keywords helper
cap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}

#' Lean strings for printing
#'
#' Returns strings with . and _ replaced by a space. This is nice when printing
#' column names of your dataframe in a report
#' @param strings vector of strings to give a nice name
#' @param check_numbers boolean indicating if numbers with decimals should be
#'   checked for and retained.
#' @keywords helper
nicename <-function (strings,check_numbers=TRUE)
{
  out <- sapply(strings, function(x) {
    original_x <- x
    x <- chartr(".", " ", x)
    x <- chartr("_", " ", x)
    if(check_numbers){
      p.positions <- gregexpr(pattern ='\\d\\.[0-9]+',original_x)[[1]]+1
      for(pos in p.positions){
        substr(x,pos,pos) <- '.'
      }

    }
    x <- gsub(" +", " ", x)
    return(x)
  })
  return(out)
}


#' Formats p-values
#'
#' Returns <0.001 if pvalue is <0.001. Else rounds the pvalue to specified
#' significant digits
#'
#' @param x an integer
#' @param digits the number of significant digits to return
#' @keywords helper
pvalue<-function(x,digits){
  if(is.na(x)|inherits(x,"character")) return(x)
  else if (x<=0.001) return("<0.001")
  else return(signif(x,digits))
}

#' Specific p-value formatting
#'
#' If p < 0.001 returns "<0.001", if p < 0.01 returns p to 3 decimal places
#' otherwise returns p to 2 decimal places
#' @param pvalues a vector of p values
#' @keywords helper
formatp<- function(pvalues){
  p_out <- sapply(pvalues, function(x){
    xsig <-suppressWarnings(as.numeric(x))
    fmtX <- ifelse(xsig<0.001,"<0.001",
                   ifelse(xsig<0.1,format(round(xsig,3),nsmall=3),
                          format(round(xsig,2),nsmall=2)))
    x <- ifelse(x=='excl','excl',fmtX)
  })
  p_out = unname(p_out)
  return(p_out)
}


sanitize <- function(str) {
  result <- str
  result <- gsub("\\\\", "SANITIZE.BACKSLASH", result)
  result <- gsub("$", "\\$", result, fixed = TRUE)
  result <- gsub(">", "$>$", result, fixed = TRUE)
  result <- gsub("<", "$<$", result, fixed = TRUE)
  result <- gsub("|", "$|$", result, fixed = TRUE)
  result <- gsub("{", "\\{", result, fixed = TRUE)
  result <- gsub("}", "\\}", result, fixed = TRUE)
  result <- gsub("%", "\\%", result, fixed = TRUE)
  result <- gsub("&", "\\&", result, fixed = TRUE)
  result <- gsub("_", "\\_", result, fixed = TRUE)
  result <- gsub("#", "\\#", result, fixed = TRUE)
  result <- gsub("^", "\\verb|^|", result, fixed = TRUE)
  result <- gsub("~", "\\~{}", result, fixed = TRUE)
  result <- gsub("SANITIZE.BACKSLASH", "$\\backslash$",
                 result, fixed = TRUE)
  return(result)
}

#'Sanitizes strings to not break LaTeX
#'
#'Strings with special characters will break LaTeX if returned 'asis' by knitr.
#'This happens every time we use one of the main reportRmd functions. We first
#'sanitize our strings with this function to stop LaTeX from breaking.
#'
#'@param str a vector of strings to sanitize
#' @keywords helper
sanitizestr<-function(str){
  as.vector(sapply(str,function(char){sanitize(char)}))
}

#' Bold strings in LaTeX
#'
#'@param strings A vector of strings to bold.
#' @keywords helper
lbld<-function(strings){sapply(strings,function(x){
  if(is.null(x)) return(x)
  if(is.na(x)) return(x)
  return(paste("\\textbf{",x,"}",sep=""))})}

#' Bold strings in HTML
#'
#'@param strings A vector of strings to bold.
#' @keywords helper
hbld<-function(strings){sapply(strings,function(x){
  if(is.null(x)) return(x)
  if(is.na(x)) return(x)
  return(paste('<span style="font-weight: bold;">',x,"</span>",sep=""))})}


ltgt <- function(x){
  # replace < and > with their html entities
  sapply(x,function(x){
    z <- try(gsub("<",'&lt;',x),silent=T)
    if (inherits(z,"try-error")){
      # Warn user and try to convert
      warning(paste("The following string contains non-ASCII charcters and may not display properly:\n",x))
      z <- try(gsub("<",'&lt;',iconv(x, to = "ASCII", sub = "")),silent=T)
      if (inherits(z,"try-error")) return(NA)
    }
    z <- gsub(">",'&gt;',z)
    return(z)
  })
}

#' Replace dollar signs with html for proper HTML output
#'
#'@param s a character vector
#'@keywords helper
rmds <- function(s){
  sapply(s,function(x){
#    x <- gsub("<0.001",'&lt;0.001',x)
    x <- ltgt(x)
    gsub("[$]",'<span style="display: inline">&#36</span>',x)
  })
}



#'Add spaces to strings in LaTeX
#'
#'Add spaces to strings in LaTeX. Returns appends ~~~ before the string
#'
#'@param x string
#'@keywords helper
addspace<-function(x){
  paste("~~~",x,sep="")
}
#' Formats p-values for LaTeX
#'
#' Returns <0.001 if pvalue is <0.001. Else rounds the pvalue to specified significant digits. Will bold the p-value if it is <= 0.05
#' @param x an integer
#' @param sigdigits number of significant digit to report
#'@keywords helper
lpvalue <- function (x, sigdigits = 2)
{
  if (is.na(x) | inherits(x,"character") )
    return(x)
  else if (x <= 0.001)
    return("\\textbf{$<$0.001}")
  else if (x <= 0.1)
    x = format(round(x, 3), nsmall = 3)
  else x = format(round(x, sigdigits), nsmall = sigdigits)
  if (x <= 0.05)
    return(paste("\\textbf{", x, "}", sep = ""))
  else return(x)
}


removedollar<-function(x){
  colnms<-strsplit(x,":")
  indx<-unlist(lapply(colnms,function(colnm) sapply(colnm, function(coln) regexpr("$",coln,fixed=T)[1]+1)))
  if(length(unique(indx))==1){
    if(unique(indx)!=0) x<-unlist(lapply(colnms,function(colnm) paste(substring(colnm,indx[1]),collapse=":")))
  }
  return(x)
}

modelmatrix<-function(f,data=NULL){
  k<-as.character(f)
  y<-NULL
  if(!length(k)%in%c(2,3)) stop("formula not properly formed")
  if(length(k)==3) {
    f<-stats::as.formula(paste("~",k[2],"+",k[3],sep=""))
    y<-stats::model.matrix(as.formula(paste("~",k[2],sep="")),data)[,-1,drop=F]}
  x<-model.matrix(f,data)[,-1,drop=F]
  colnames(x)<-removedollar(colnames(x))
  if(!is.null(y)){
    return(list(x[,1:ncol(y),drop=F],x[,(ncol(y)+1):ncol(x),drop=F]))
  }else{
    return(x)
  }}

nicecall <- function(model_call) {
  call_str <- deparse1(model_call)
  call_str <- gsub("[\"]","'",call_str)
  return(call_str)
}
matchdata <- function(dataArg){
  df_str <- as.character(dataArg)
  if (length(df_str)>1) df_str = df_str[2]
  no_fnc <- gsub("[A-Za-z]+[(]","",df_str)
  txt_bts <- unlist(strsplit(no_fnc,split = "[^A-Za-z0-9_.]"))
  txt_bts <- txt_bts[txt_bts!=""]
  obj <- intersect(txt_bts,ls(name=".GlobalEnv"))
  if (length(obj)>0){
    dfInd <-sapply(obj,function(x)inherits(get0(x),'data.frame'))
    df <- obj[dfInd]
    if (length(df)>1){
      message("Multiple data objects found in function call")
      return(NULL)
    } else return(df)
  } else {
    message("Data object could not be extracted from function call")
    return(NULL)
  }
}

matchcovariate=function(betanames,ucall){
  out=as.vector(sapply(betanames,function(betaname){
    splitbetaname=unlist(strsplit(betaname,":",fixed=T))
    out=sapply(splitbetaname,function(bname){
      bname=gsub(" ","",bname) # added 14 Dec 2020 to allow matching with centred variables
      #indx=which(sapply(ucall,function(cov)charmatch(cov,bname))==1)
      indx=which(sapply(ucall,function(cov)grepl(cov,bname,fixed=TRUE))) ## changed on Feb 21, 2019 for checkings
      if(length(indx)==1)return(indx)
      #If one  facorname is a subset of another
      indx2<-which.max(sapply(ucall[indx],nchar))
      if(length(indx2)==1) return(indx[indx2])
      indx3<-which(sapply(ucall[indx2],function(c){substr(betaname,1,nchar(c))==c}))
      if(length(indx3)==1)  return(ucall[indx[indx2[indx3]]])
      return(-1)
    })
    if(-1 %in% out) return(-1)
    result=0
    n=length(out)
    for(i in 1:length(out)){
      result=result+out[i]*100^(n-1)
      n=n-1
    }
    return(result)}))
  if(-1 %in% out) return(-1)
  return (out)
}

# Adapted From the CAR package to compute VIF
GVIF <- function(model){
  v <- vcov(model)
  ind <- attr(model.matrix(model), "assign")
  if (0 %in% ind) {
    v <- v[-1, -1]
    ind <- ind[-1]
  }
  xvar <- labels(terms(model))
  if (length(xvar)<2) {
    return(data.frame(Covariate=xvar,VIF=NA))
  }

  R <- stats::cov2cor(v)
  detR <- det(R)
  result <- matrix(0, length(xvar), 2)
  for (var in 1:length(xvar)) {
    terms <- which(ind == var)
    result[var, 1] <- det(as.matrix(R[terms, terms])) * det(as.matrix(R[-terms,
                                                                        -terms]))/detR
    result[var, 2] <- length(terms)
  }
  if (all(result[, 2] == 1)){
    rtn <- result[, 1]
  } else rtn <- result[, 1]^(1/(2 * result[, 2]))
  data.frame(Covariate=xvar,VIF=rtn)
}

# (ggsurv) ---------------------------------------------------------

round_sprintf <- function(value,digits){
  sprintf( paste0("%.",digits,"f"), round(value,digits))
}

pstprn0 <- function (x)
{
  paste0(x[1], "(", paste0(x[-1], collapse = ","),
         ")", sep = "")
}

psthr0 <- function (x, digits = 2)
{
  x <- sapply(x, function(x) {
    ifelse(abs(x) < 0.01 | abs(x) > 1000, format(x, scientific = TRUE,
                                                 digits = digits), round_sprintf(x, digits))
  })
  pstprn0(x)
}

break_function <- function(xmax){

  xmax_length <- ifelse(xmax>1,nchar(round(xmax)),round(abs(log10(xmax))))

  byx <- if(xmax>1) {round(xmax/10,digits = 2-xmax_length)
  }else round(xmax/10,digits = xmax_length+1)

  breaks <- seq(0,xmax,by=byx)
  if(max(breaks)<byx) breaks <- c(breaks,max(breaks)+byx)
  return(breaks)
}

lpvalue2 <- function (x,digits)
{
  if (is.na(x) |    inherits(x,"character") )
    return(x)
  else if (x < 10^-(digits))
    return(paste0("p < ",10^-(digits)))
  else return(paste0("p = ",round_sprintf(x, digits)))

}

.extract_ggplot_colors <- function(p, grp.levels){
  g <- ggplot_build(p)
  .cols <- unlist(unique(g$data[[1]]["colour"]))
  if(!is.null(grp.levels)){
    if(length(.cols)==1) .cols <- rep(.cols, length(grp.levels))
    names(.cols) <- grp.levels
  }
  .cols
}

.set_large_dash_as_ytext <- function(ggp){
  ggp + theme(axis.text.y = element_text(size = 50, vjust = 0.35),
              axis.ticks.y = element_blank())
}

##This function is used by the survfit package
survfit_confint <- function(p, se, logse=TRUE, conf.type, conf.int=0.95,
                            selow, ulimit=TRUE) {
  zval <- qnorm(1- (1-conf.int)/2, 0,1)
  if (missing(selow)) scale <- 1.0
  else scale <- ifelse(selow==0, 1.0, selow/se)  # avoid 0/0 at the origin
  if (!logse) se <- ifelse(se==0, 0, se/p)   # se of log(survival) = log(p)

  if (conf.type=='plain') {
    se2 <- se* p * zval  # matches equation 4.3.1 in Klein & Moeschberger
    if (ulimit) list(lower= pmax(p -se2*scale, 0), upper = pmin(p + se2, 1))
    else  list(lower= pmax(p -se2*scale, 0), upper = p + se2)
  }
  else if (conf.type=='log') {
    #avoid some "log(0)" messages
    xx <- ifelse(p==0, NA, p)
    se2 <- zval* se
    temp1 <- exp(log(xx) - se2*scale)
    temp2 <- exp(log(xx) + se2)
    if (ulimit) list(lower= temp1, upper= pmin(temp2, 1))
    else  list(lower= temp1, upper= temp2)
  }
  else if (conf.type=='log-log') {
    xx <- ifelse(p==0 | p==1, NA, p)
    se2 <- zval * se/log(xx)
    temp1 <- exp(-exp(log(-log(xx)) - se2*scale))
    temp2 <- exp(-exp(log(-log(xx)) + se2))
    list(lower = temp1 , upper = temp2)
  }
  else if (conf.type=='logit') {
    xx <- ifelse(p==0, NA, p)  # avoid log(0) messages
    se2 <- zval * se *(1 + xx/(1-xx))

    temp1 <- 1- 1/(1+exp(log(p/(1-p)) - se2*scale))
    temp2 <- 1- 1/(1+exp(log(p/(1-p)) + se2))
    list(lower = temp1, upper=temp2)
  }
  else if (conf.type=="arcsin") {
    xx <- ifelse(p==0, NA, p)
    se2 <- .5 *zval*se * sqrt(xx/(1-xx))
    list(lower= (sin(pmax(0, asin(sqrt(xx)) - se2*scale)))^2,
         upper= (sin(pmin(pi/2, asin(sqrt(xx)) + se2)))^2)
  }
  else stop("invalid conf.int type")
}


color_palette_surv_ggplot <- function(length){
  if(length==1) return("black")
  if(length==2) return(c("#D53E4F","#3288BD"))
  if(length==3) return(c("#D53E4F","#ABDDA4","#3288BD"))
  if(length==4) return(c("#D53E4F","#FDAE61","#ABDDA4","#3288BD"))
  if(length==5) return(c("#D53E4F","#FDAE61","#FEE08B","#ABDDA4","#3288BD"))
  if(length==6) return(c("#D53E4F","#FDAE61","#FEE08B","#ABDDA4","#66C2A5","#3288BD"))
  if(length==7) return(c("#D53E4F","#F46D43","#FDAE61","#FEE08B","#ABDDA4","#66C2A5","#3288BD"))
  if(length==8) return(c("#D53E4F","#F46D43","#FDAE61","#FEE08B","#ABDDA4","#66C2A5","#3288BD","#5E4FA2"))
  if(length==9) return(c("#9E0142","#D53E4F","#F46D43","#FDAE61","#FEE08B","#ABDDA4","#66C2A5","#3288BD","#5E4FA2"))
  if(length==10) return(c("black","#9E0142","#D53E4F","#F46D43","#FDAE61","#FEE08B","#ABDDA4","#66C2A5","#3288BD","#5E4FA2"))
  if(length>10) {message("10 colours maximum in default")}
  return(rep(c("black","#9E0142","#D53E4F","#F46D43","#FDAE61","#FEE08B","#ABDDA4","#66C2A5","#3288BD","#5E4FA2"),length.out=length))
}


# (forestplot2) ---------------------------------------------------------
format_glm = function(glm_fit,conf.level = 0.95,digits=c(2,3),orderByRisk=TRUE){
  if (! class(glm_fit)[1] %in% c('glm','geeglm','polr')) stop('Only objects of class glm, geeglm and polr are accepted.')

  #extracting ORs and p values
  Z = stats::qnorm(1-(1-conf.level)/2)
  tab <- as.data.frame(summary(glm_fit)$coefficients)
  tab <- cbind(variable= rownames(tab),tab)
  rownames(tab) <- NULL

  if (class(glm_fit)[1] %in% c("glm", "geeglm")){
    names(tab) =  c("variable","estimate",  "std.error" ,"statistic", "p.value")
    tab = tab[-which(tab$variable=='(Intercept)'),]
  }  else {
    names(tab) =  c("variable","estimate",  "std.error" ,"statistic")
    tab$coef.type = ifelse(grepl("[|]",tab$variable),"scale","coefficient")
    tab <- tab[tab$coef.type=='coefficient',]
    tab$p.value = stats::pnorm(abs(tab$statistic),lower.tail = FALSE) * 2
  }

  tab$conf.low=exp(tab$estimate-Z*tab$std.error)
  tab$conf.high=exp(tab$estimate+Z*tab$std.error)
  tab$estimate = exp(tab$estimate)
  tab$estimate.label = paste0(niceNum(tab$estimate), ' (',niceNum(tab$conf.low),', ',niceNum(tab$conf.high),')')


  tab$p.label = ifelse(tab$p.value<0.001, '<0.001', niceNum(tab$p.value,digits[2]))
  names(tab)[1] = 'variable'

  tab = tab[,c('variable', 'estimate', 'p.label', 'p.value', 'conf.low', 'conf.high')]


  if (orderByRisk){
    tab$var.order = rank(tab$estimate)
  } else{
    tab$var.order = 1:nrow(tab)
  }

  # Extract the reference levels if needed
  if (length(glm_fit$xlevels)!=0){
    ref_levels <- NULL
    for (i in seq_along(glm_fit$xlevels)){
      ref_levels <- rbind(ref_levels,
                          data.frame(var.name=rep(names(glm_fit$xlevels)[i],length(glm_fit$xlevels[[i]])+1),
                                     level.name = c(names(glm_fit$xlevels)[i],glm_fit$xlevels[[i]]),
                                     level.order=1:(length(glm_fit$xlevels[[i]])+1),
                                     variable=paste0(names(glm_fit$xlevels)[i],c('',glm_fit$xlevels[[i]]))))
    }


    tab = merge(ref_levels, tab, by='variable',all = T)

    tab$estimate.label = ifelse(is.na(tab$estimate), '1.0 (Reference)',
                                paste0(niceNum(tab$estimate), ' (',niceNum(tab$conf.low),', ',niceNum(tab$conf.high),')'))

    varOrders <- tapply(X = tab$var.order,
                        INDEX=tab$var.name,
                        FUN = function(x) min(x,na.rm=T))
    varOrderLookup <- data.frame(var.name=names(varOrders),var.order=varOrders)


    varOrderLookup <- stats::na.omit(tab[,c("var.name","var.order")])

    for (i in 1:nrow(varOrderLookup)){
      tab$var.order[tab$var.name==varOrderLookup$var.name[i]] <- varOrderLookup$var.order[i]
    }

    tab$estimate.label = ifelse(tab$level.name %in% names(glm_fit$xlevels),NA_character_,tab$estimate.label)
    tab[order(tab$var.order,tab$level.order,decreasing=c(F,T)),]
  } else {
    tab$estimate.label = paste0(niceNum(tab$estimate), ' (',niceNum(tab$conf.low),', ',niceNum(tab$conf.high),')')
    tab$level.order=1
    tab$var.name=tab$variable
    tab$level.name=tab$variable
    tab[order(tab$var.order),]
  }

}

# New function to strip centering from a covariate
getvarname = function(betaname){
  sapply(betaname,function(x){
    x = gsub('I[(]','',x)
    x = gsub('[-+].*','',x)
    x = trimws(x)
    return(x)
  })
}

lbl_count <- function(y){
  q75 <- summary(y)[5]
  return(data.frame(y=max(y),  label=paste('n =',length(y))))
}

betaWithCI <-function(betaname,CIwidth=0.95){
  paste0(betaname,"(",100*CIwidth,"%CI)")
}

niceStr <- function (strings)
{
  out <- sapply(strings, function(x) {
    x <- chartr('/',' ',x)
    x <- chartr(".", " ", x)
    x <- chartr("_", " ", x)
    return(x)
  })
  return(out)
}

wrp_lbl <- function(x,width = 10){
  x <- niceStr(x)
  #  strwrap(x,width = width) # doesn't work nicely with spaces
  lst <- strwrap(x,width = width,simplify = F)
  for (i in seq_along(lst)) lst[[i]] <- paste(lst[[i]],collapse='\n')
  unlist(lst)
}


label_wrap_reportRx <- function (width = 25, multi_line = TRUE) {
  fun <- function(labels) {
    labels <- label_value(labels, multi_line = multi_line)
    lapply(labels, function(x) {
      x <- niceStr(x)
      x <- strwrap(x, width = width, simplify = FALSE)
      vapply(x, paste, character(1), collapse = "\n")
    })
  }
  structure(fun, class = "labeller")
}






reportRx_pal <- function(
    direction = 1
) {

  function(n) {
    if (n>10) warning('Ten colour maximum, colours will be recycled.')

    colour_list <- color_palette_surv_ggplot(n)

    colour_list <- unname(unlist(colour_list))
    if (direction >= 0) colour_list else rev(colour_list)
  }
}

scale_colour_reportRx <- function(
    direction = 1,
    ...
) {
  ggplot2::discrete_scale(
    aesthetics = c("colour","fill"),
    scale_name = "reportRx",
    reportRx_pal( direction),
    ...
  )
}

fillNAs <- function(x) {
  ind = which(!is.na(x))
  if(is.na(x[1]))
    ind = c(1,ind)
  rep(x[ind], times = diff(c(ind, length(x) + 1) ))
}


#' Extract Function and Package Information from Current Document
#'
#' This function analyses the current file (an R script, Rmd or qmd file) to
#' extract information about all functions called within the code, identifies
#' their associated packages, and returns a summary of packages used with
#' version and citation information.
#'
#' @description The function automatically detects the current R script file
#' (works best in RStudio), parses the code to identify function calls,
#' determines which packages they belong to, and creates a summary of all
#' non-base R packages used in the script. It handles both namespace-qualified
#' function calls (e.g., `dplyr::filter`) and regular function calls, while
#' filtering out base R functions and control structures.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{package_name}{Character. Name of the package}
#'   \item{functions_called}{Character. Comma-separated list of functions called from this package}
#'   \item{package_version}{Character. Version number of the installed package}
#'   \item{package_citation}{Character. Formatted citation for the package}
#' }
#'
#' @note
#' \itemize{
#'   \item Works best when run from RStudio with an active source file
#'   \item Requires that referenced packages are already loaded/installed
#'   \item Will not detect functions called through indirect methods (e.g., `do.call()`)
#' }
#'
#' @examples
#' \dontrun{
#' # Run this function from within an R script to analyze its dependencies
#' package_info <- extract_package_details()
#' print(package_info)
#' }
#'
#' @seealso \code{\link[utils]{getAnywhere}},
#' \code{\link[utils]{packageVersion}}, \code{\link[utils]{citation}}
#'
#' @importFrom utils getAnywhere packageVersion citation
#' @importFrom  rstudioapi isAvailable getSourceEditorContext
#' @importFrom dplyr distinct group_by summarise mutate
#' @export
extract_package_details <- function() {

  # Get the current file path
  # Try different methods to detect the current file
  get_current_file <- function() {
    # Try rstudioapi first (works in RStudio)
    if (requireNamespace("rstudioapi", quietly = TRUE)) {
      if (rstudioapi::isAvailable()) {
        file_path <- rstudioapi::getSourceEditorContext()$path
        if (!is.null(file_path) && file_path != "") {
          return(file_path)
        }
      }
    }
    # Try to get the currently executing script
    # This works when source() is used
    if (exists("ofile", envir = parent.frame())) {
      return(get("ofile", envir = parent.frame()))
    }

    # Try sys.frame approach
    for (i in sys.nframe():1) {
      if (exists("ofile", envir = sys.frame(i))) {
        return(get("ofile", envir = sys.frame(i)))
      }
    }

    # If all else fails, prompt user
    stop("Could not detect current file. Please ensure you're running this from RStudio with an active file, or use rstudioapi package.")
  }

  try(file_path <- get_current_file(),silent = T)
  if (inherits("try-error",file_path)) stop("Current file can not be identified.")


  # Read the file content
  text <- readLines(file_path, warn = FALSE) |>
    paste(collapse = "\n")

  # Remove HTML comments (<!-- -->)
  text <- gsub("<!--.*?-->", "", text, perl = TRUE)

  # Split into lines for processing
  lines <- strsplit(text, "\n")[[1]]

  # Remove comments (everything after # on each line)
  # But preserve # inside strings
  remove_comments <- function(line) {
    # Simple approach: find # not inside quotes
    in_single_quote <- FALSE
    in_double_quote <- FALSE
    chars <- strsplit(line, "")[[1]]
    result <- character()

    for (i in seq_along(chars)) {
      char <- chars[i]

      # Check for escape sequences
      if (i > 1 && chars[i-1] == "\\") {
        result <- c(result, char)
        next
      }

      # Toggle quote states
      if (char == "'" && !in_double_quote) {
        in_single_quote <- !in_single_quote
      } else if (char == '"' && !in_single_quote) {
        in_double_quote <- !in_double_quote
      } else if (char == "#" && !in_single_quote && !in_double_quote) {
        # Found a comment outside of quotes
        break
      }

      result <- c(result, char)
    }

    paste(result, collapse = "")
  }

  # Apply comment removal to each line
  lines <- sapply(lines, remove_comments, USE.NAMES = FALSE)

  # Combine back into single text
  text <- paste(lines, collapse = "\n")

  # First, remove function definitions to avoid false positives
  # Match "function(" only when preceded by space, =, <-, or start of line
  text_no_func_def <- gsub("(^|\\s|=|<-)function\\s*\\(", "\\1FUNCTION_DEFINITION(", text, perl = TRUE)

  # Initialize results list
  all_functions <- list()

  # Pattern 1: Namespace-qualified function calls (package::function or package:::function)
  ns_pattern <- "\\b([a-zA-Z][a-zA-Z0-9.]*)(:::|::)([a-zA-Z_][a-zA-Z0-9._]*)\\s*\\("
  ns_matches <- gregexpr(ns_pattern, text_no_func_def, perl = TRUE)
  ns_calls <- regmatches(text_no_func_def, ns_matches)[[1]]

  if (length(ns_calls) > 0) {
    # Extract package and function names
    for (call in ns_calls) {
      # Remove the opening parenthesis and whitespace
      call_clean <- gsub("\\s*\\($", "", call)

      # Split by :: or :::
      parts <- strsplit(call_clean, ":::|::")[[1]]
      if (length(parts) == 2) {
        all_functions[[length(all_functions) + 1]] <- list(
          function_name = parts[2],
          package_name = parts[1],
          is_namespaced = TRUE
        )
      }
    }
  }

  # Pattern 2: Regular function calls (not namespace-qualified)
  # Exclude calls that are preceded by :: or :::
  regular_pattern <- "(?<!:)(?<!::|:::)\\b([a-zA-Z_][a-zA-Z0-9._]*)\\s*\\("
  regular_matches <- gregexpr(regular_pattern, text_no_func_def, perl = TRUE)
  regular_calls <- regmatches(text_no_func_def, regular_matches)[[1]]

  if (length(regular_calls) > 0) {
    # Extract just the function names (remove the trailing "(")
    function_names <- gsub("\\s*\\($", "", regular_calls)

    # Remove some common non-function patterns
    control_structures <- c("if", "for", "while", "repeat", "function",
                            "case_when","c","list","data.frame","tibble",
                            "FUNCTION_DEFINITION", "switch")
    function_names <- function_names[!function_names %in% control_structures]

    # Add to list
    for (func_name in function_names) {
      all_functions[[length(all_functions) + 1]] <- list(
        function_name = func_name,
        package_name = NA,
        is_namespaced = FALSE
      )
    }
  }

  # Convert to data frame for easier manipulation
  if (length(all_functions) == 0) {
    return(data.frame(
      function_name = character(),
      package_name = character(),
      package_version = character(),
      stringsAsFactors = FALSE
    ))
  }

  functions_df <- do.call(rbind, lapply(all_functions, as.data.frame))

  # Get unique function-package combinations
  functions_df <- functions_df |> dplyr::distinct()

  # Check if not namespaced functions are base R function, and if so, remove
  is_base_R <- function(func_name){
    sapply(func_name, function(f){
      exists(f, mode = "function", envir = baseenv())
    })
  }

  functions_df <- functions_df |>
    dplyr::mutate(
      package_name = ifelse(is.na(package_name),ifelse(is_base_R(function_name),"base",NA),package_name)) |>
    dplyr::filter(!grepl("base",package_name)) |>
    dplyr::filter(function_name != "extract_package_details")

  # Get package information for non-namespaced functions
  get_function_package <- function(func_names) {
    sapply(func_names,function(func_name) {
      # find the function
      where_found <- getAnywhere(func_name)
      if (length(where_found$where)>0){
        if ( ".GlobalEnv" %in% where_found$where) return("GlobalEnv")
        pkg_list <- unique(gsub("package[:]|namespace[:]","",
                                grep("package|namespace",where_found$where,value = T)))
        return(pkg_list[1])
      } else return("Unknown")
    },USE.NAMES=F,simplify=T)
  }

  functions_df <- functions_df |>
    dplyr::mutate(
      package_name = ifelse(is.na(package_name),
                            get_function_package(function_name),package_name)) |>
    dplyr::filter(!grepl("Unknown",package_name))

  # Function to get package version
  get_package_version <- function(pkg_name) {
    sapply(pkg_name, function(pn){
      tryCatch(
        as.character(packageVersion(pn)),
        error = function(e) NA_character_
      )
    })
  }
  # Function to get package citation
  get_package_citation <- function(pkg_name) {
    sapply(pkg_name, function(pn){
      tryCatch(
        format(citation(pn),style = "text"),
        error = function(e) NA_character_
      )
    })
  }
  packages_df <- functions_df |>
    dplyr::group_by(package_name) |>
    dplyr::summarise(functions_called = paste(unique(function_name),collapse=", ")) |>
    dplyr::ungroup() |>
    dplyr::filter(package_name !="GlobalEnv") |>
    dplyr::add_row(package_name ="utils") |>
    dplyr::mutate(package_version = get_package_version(package_name),
                  package_citation = get_package_citation(package_name)) |>
    dplyr::group_by(package_citation) |>
    dplyr::slice_tail(n=1) |>
    dplyr::ungroup() |>
    dplyr::mutate(package_name = gsub("utils","R",package_name))
  ord <- c((1:nrow(packages_df))[-which(packages_df$package_name=="R")],which(packages_df$package_name=="R"))
  packages_df <-packages_df[ord,]
  return(packages_df)
}

