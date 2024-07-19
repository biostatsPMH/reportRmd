getVarLevels <- function(model){
  terms<-names(model$coefficients)[!grepl("intercept",names(model$coefficients),ignore.case = T)]
  vrs<-attr(model$terms,"term.labels")
  lvls<-setdiff(names(model$coefficients),vrs)
  df <- data.frame(terms=terms)
  df$var <- ifelse(df$term %in% vrs,df$term,NA)
  df$lvl <- ifelse(df$term %in% lvls,df$term,NA)
  if (any(is.na(df$var))){
    nt <- function(str) length(strsplit(str,":")[[1]])-1
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
    df$var <-vr2
    lvl2 <- mapply(function(v,l){
      v=strsplit(v,"[:]")[[1]]
      l=strsplit(l,"[:]")[[1]]
      paste0(mapply(function(v,l) sub(v,"",l),v,l,USE.NAMES = F),collapse = ":")
    }, df$var,df$lvl,USE.NAMES=F)
    df$lvl=lvl2
  }
  df
}
