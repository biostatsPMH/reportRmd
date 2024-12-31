mterms <- function(model) {
  UseMethod("mterms", model)
}
mterms.default <- function(model){
  names(model$coefficients)[!grepl("intercept",
                                   names(model$coefficients),ignore.case = T)]
}

mterms.lmerModLmerTest <- function(model){
  if (requireNamespace("nlme", quietly = TRUE)) {
  names(nlme::fixef(model))[!grepl("intercept",
                              names(nlme::fixef(model)),ignore.case = T)]
  } else stop("Summarising mixed effects models requires the nlme package be installed")
}

getVarLevels <- function(model){
  ord <- NULL
  nt <- function(str) length(strsplit(str,":")[[1]])-1
  vrs<-try(attr(model$terms,"term.labels"),silent = T)
  if (inherits(vrs,"try-error"))  vrs<-try(attr(terms(model),"term.labels"))
  if (inherits(vrs,"try-error")) stop("Model terms could not be found.")
  if (any(sapply(vrs,nt)>1)) stop("Summary functions will not work with three-way interactions.")
  terms <- mterms(model)
  lvls<-setdiff(terms,vrs)
  df <- data.frame(terms=terms)
  df$var <- ifelse(df$term %in% vrs,df$term,NA)
  df$lvl <- ifelse(df$term %in% lvls,df$term,NA)
  int_terms <- grep("[:]",vrs,value=T)
  for (v in int_terms){
    vr <- unlist(strsplit(v, ":"))
    p <- paste0("^", vr[1], ".*:", vr[2], ".*$")
    vind <- grep(p, df$terms)
    if (all(is.na(df$var[vind]))) df$var[vind] <- v
  }
  if (any(is.na(df$var))){
    vind <- lapply(vrs,function(v) which(grepl(paste0("^",v),df$terms)|grepl(paste0("[:]",v),df$terms)))
    names(vind) <- vrs
    vind <- vind[unlist(lapply(vind,length))>0]
    vr <- character(max(unlist(vind)))
    for (name in names(vind)) {
      indices <- vind[[name]]
      for (index in indices) {
        if (vr[index] == "") {
          vr[index] <- name
        } else {
          vr[index] <- paste0(vr[index],":", name)
        }
      }
    }
    vr2 <- ifelse(grepl("[:]",df$terms),vr,gsub("[:].*","",vr))
    df$var <- ifelse(is.na(df$var),vr2,df$var)
    lvl2 <- mapply(function(v,l){
      v=strsplit(v,"[:]")[[1]]
      l=strsplit(l,"[:]")[[1]]
      paste0(mapply(function(v,l) sub(v,"",l),v,l,USE.NAMES = F),collapse = ":")
    }, df$var,df$lvl,USE.NAMES=F)
    df$lvl=sub("^:","",lvl2)
  }
  df$lvl <- ifelse(df$lvl=="NA:NA",NA,df$lvl)
  df$lvl <- ifelse(df$lvl=="NA",NA,df$lvl)
  df$ord <- 1:nrow(df)

  # add sample size and events
  md <- get_model_data(model)
  if (is.null(md)){
    warning("Model data could not be extracted, simplified summary provided")
  } else {
    events = get_event_counts(model)
    if (!is.null(events)) ed <- md[events==1,] else ed <- NULL
  }
  df$v1 <- sub(":.*","",df$var)
  df$v2 <- sub(".*:","",df$var)
  df$v2 <- ifelse(df$v1==df$v2,NA,df$v2)

  vcls <- sapply(na.omit(unique(c(df$v1,df$v2))), function(x) ifelse(is.numeric(md[[x]]),"numeric","factor"))
  int_terms <- unique(grep("[:]",df$var,value=T))
  freq1 <- NULL
  for (i in int_terms){
    vcl <- vcls[strsplit(i,":")[[1]]]
    if (all(vcl=="factor")){
      if (!is.null(md)){
        ntbl <- eval(parse(text=paste("data.frame(with(md,table(",sub(":",',',i),")))")))
        ntbl$var <- i
        ntbl$lvl <- paste0(ntbl[,1],":",ntbl[,2])
        if (!is.null(ed)){
          etbl <- eval(parse(text=paste("data.frame(with(ed,table(",sub(":",',',i),")))")))
          etbl$var <- i
          etbl$lvl <- paste0(ntbl[,1],":",ntbl[,2])
          names(etbl) <- sub("Freq","Events",names(etbl))
          ntbl <- merge(ntbl,etbl)
        }
      } else {ntbl <- data.frame(var=i)}

    } else if (any(vcl=="factor")){
      if (!is.null(md)){
        ntbl <- eval(parse(text=paste("data.frame(with(md,table(",names(vcl)[vcl=="factor"],")))")))
        ntbl$var <- i
        ntbl$lvl <- ntbl[,1]
        if (!is.null(ed)){
          etbl <- eval(parse(text=paste("data.frame(with(ed,table(",names(vcl)[vcl=="factor"],")))")))
          etbl$var <- i
          etbl$lvl <- etbl[,1]
          names(etbl) <- sub("Freq","Events",names(etbl))
          ntbl <- merge(ntbl,etbl)
        }
      } else {ntbl <- data.frame(var=i)}
    } else{
      if (!is.null(md)){
        ntbl <- data.frame(var=i,lvl=NA,Freq=nrow(md))
        if (!is.null(ed)){
          etbl <- data.frame(var=i,lvl=NA,Events=nrow(ed))
          ntbl <- merge(ntbl,etbl)
        }
      } else {ntbl <- data.frame(var=i)}
    }
    freq1 <- dplyr::bind_rows(freq1,ntbl)
  }
  freq2 <- NULL
  for (i in setdiff(na.omit(unique(df$var)),int_terms)){
    vcl <- vcls[i]
    if (vcl=="factor"){
      if (!is.null(md)){
        ntbl <- eval(parse(text=paste("data.frame(with(md,table(",i,")))")))
        ntbl$var <- i
        ntbl$lvl <- ntbl[,1]
        if (!is.null(ed)){
          etbl <- eval(parse(text=paste("data.frame(with(ed,table(",i,")))")))
          etbl$var <- i
          etbl$lvl <- etbl[,1]
          names(etbl) <- sub("Freq","Events",names(etbl))
          ntbl <- merge(ntbl,etbl)
        }
      } else {ntbl <- data.frame(var=i)}
    } else {
      if (!is.null(md)){
        ntbl <- data.frame(var=i,lvl=NA,Freq=nrow(md))
        if (!is.null(ed)){
          etbl <- data.frame(var=i,lvl=NA,Events=nrow(ed))
          ntbl <- merge(ntbl,etbl)
        }
      } else {ntbl <- data.frame(var=i)}
    }
    freq2 <- dplyr::bind_rows(freq2,ntbl)
  }
  freq <- dplyr::bind_rows(freq1,freq2)
  freq$n <- freq$Freq
  freq <- freq[,intersect(names(freq),c("var","lvl","n","Events")),drop=FALSE]
  out <- suppressMessages(dplyr::full_join(df,freq))
  out <- dplyr::arrange(out,ord)
  vord <- data.frame(var=na.omit(unique(out$var)),
                     vord=1:length(na.omit(unique(out$var))))
  out <- suppressMessages(dplyr::full_join(out,vord))
  out$ord <- ifelse(is.na(out$ord),0,out$ord)
  out <- dplyr::arrange(out,vord,ord)
  out$ref <- out$ord==0
  extra_lvl <- xtr_lvls(out$ref )
  rtn <- out[setdiff(1:nrow(out),extra_lvl),
             intersect(c("var","lvl","n","ref","Events","terms"),names(out))]
  return(rtn)
}

xtr_lvls <- function(vec){
  nx_vec <- c(vec[-1],FALSE)
  pr_vec <- c(FALSE,vec[-length(vec)])
  cnst_T <- which(vec & nx_vec | vec & pr_vec)
  return(cnst_T[c(0,diff(cnst_T))==1])
}


