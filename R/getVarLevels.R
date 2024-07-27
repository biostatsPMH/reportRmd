getVarLevels <- function(model){
  nt <- function(str) length(strsplit(str,":")[[1]])-1
  terms<-names(model$coefficients)[!grepl("intercept",names(model$coefficients),ignore.case = T)]
  vrs<-attr(model$terms,"term.labels")
  if (any(sapply(vrs,nt)>1)) stop("Summary functions will not work with three-way interactions.")
  lvls<-setdiff(names(model$coefficients),vrs)
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
  df$v1 <- sub(":.*","",df$var)
  df$v2 <- sub(".*:","",df$var)
  df$v2 <- ifelse(df$v1==df$v2,NA,df$v2)
  data= get_model_data(model)
  vcls <- sapply(na.omit(unique(c(df$v1,df$v2))), function(x) ifelse(inherits(data[[x]],"numeric"),"numeric","factor"))
  int_terms <- unique(grep("[:]",df$var,value=T))
  freq1 <- NULL
  for (i in int_terms){
    vcl <- vcls[strsplit(i,":")[[1]]]
    if (all(vcl=="factor")){
      ntbl <- eval(parse(text=paste("data.frame(with(data,table(",sub(":",',',i),")))")))
      ntbl$var <- i
      ntbl$lvl <- paste0(ntbl[,1],":",ntbl[,2])
    } else if (any(vcl=="factor")){
      ntbl <- eval(parse(text=paste("data.frame(with(data,table(",names(vcl)[vcl=="factor"],")))")))
      ntbl$var <- i
      ntbl$lvl <- ntbl[,1]
    } else{
      ntbl <- data.frame(var=i,lvl=NA,Freq=nrow(data))
    }
    freq1 <- dplyr::bind_rows(freq1,ntbl)
  }
  freq2 <- NULL
  for (i in setdiff(na.omit(unique(df$var)),int_terms)){
    vcl <- vcls[i]
    if (vcl=="factor"){
      ntbl <- eval(parse(text=paste("data.frame(with(data,table(",i,")))")))
      ntbl$var <- i
      ntbl$lvl <- ntbl[,1]
    } else ntbl <- data.frame(var=i,lvl=NA,Freq=nrow(data))
    freq2 <- dplyr::bind_rows(freq2,ntbl)
  }
  freq <- dplyr::bind_rows(freq1,freq2)
  freq$n <- freq$Freq
  freq <- freq[,c("var","lvl","n")]
  out <- dplyr::full_join(df,freq)
  out <- dplyr::arrange(out,ord)
  vord <- data.frame(var=na.omit(unique(out$var)),
                                 vord=1:length(na.omit(unique(out$var))))
  out <- dplyr::full_join(out,vord)
  out$ord <- ifelse(is.na(out$ord),0,out$ord)
  out <- dplyr::arrange(out,vord,ord)
  out$ref <- out$ord==0
  extra_lvl <- xtr_zrs(!out$ref )
  out <- out[!extra_lvl,c("terms","var","lvl","n","ref")]
}

xtr_zrs <- function(vec) {
  # Create a logical vector where consecutive zeroes are marked
  is_zero <- vec == 0
  is_xtr_zero <- is_zero & c(FALSE, is_zero[-length(is_zero)])

  return(is_xtr_zero & !c(FALSE, !is_xtr_zero[-length(is_xtr_zero)]))
}




# getVarLevels <- function(model){
#   terms<-names(model$coefficients)[!grepl("intercept",names(model$coefficients),ignore.case = T)]
#   vrs<-attr(model$terms,"term.labels")
#   lvls<-setdiff(names(model$coefficients),vrs)
#   df <- data.frame(terms=terms)
#   df$var <- ifelse(df$term %in% vrs,df$term,NA)
#   df$lvl <- ifelse(df$term %in% lvls,df$term,NA)
#   if (any(is.na(df$var))){
#     nt <- function(str) length(strsplit(str,":")[[1]])-1
#     vind <- lapply(vrs,function(v) which(grepl(paste0("^",v),df$terms)|grepl(paste0("[:]",v),df$terms)))
#     names(vind) <- vrs
#     vind <- vind[unlist(lapply(vind,length))>0]
#     vr <- character(max(unlist(vind)))
#     for (name in names(vind)) {
#       indices <- vind[[name]]
#       for (index in indices) {
#         if (vr[index] == "") {
#           vr[index] <- name
#         } else {
#           vr[index] <- paste0(vr[index],":", name)
#         }
#       }
#     }
#     vr2 <- ifelse(grepl("[:]",df$terms),vr,gsub("[:].*","",vr))
#     df$var <-vr2
#     lvl2 <- mapply(function(v,l){
#       v=strsplit(v,"[:]")[[1]]
#       l=strsplit(l,"[:]")[[1]]
#       paste0(mapply(function(v,l) sub(v,"",l),v,l,USE.NAMES = F),collapse = ":")
#     }, df$var,df$lvl,USE.NAMES=F)
#     df$lvl=lvl2
#   }
#   df
# }
