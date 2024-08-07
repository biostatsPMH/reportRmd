# Combind model components and variable levels and sample sizes

m_summary <- function(model,CIwidth=.95,digits=2,vif = FALSE,whichp="level", for_plot = FALSE){

  m_coeff <- coeffSum(model,CIwidth,digits)
  if (any(!is.na(m_coeff$lwr) & !is.na(m_coeff$upr) & (m_coeff$lwr == m_coeff$upr))) message("Zero-width confidence interval detected. Check predictor units.")
  lvls <- getVarLevels(model)
  lvls$ord  <- 1:nrow(lvls)

  cs <- merge(lvls,m_coeff,all.x = T)

  cs <- cs[order(cs$ord),]
  rownames(cs) <- NULL
  # add variable header rows
  cs$header <- NA
  hdr_vars <- names(table(cs$var))[table(cs$var)>1]
  for (v in hdr_vars){
    hdr_rw <- sort(which(cs$var %in% v))[1]
    cs <- dplyr::add_row(cs,.before = hdr_rw,var=v,header=TRUE)
  }
  for (var in unique(na.omit(cs[["var"]]))) {
    if (var %in% cs[["var"]] & length(which(cs[["var"]] == var)) == 3) {
      p <- cs[max(which(cs[["var"]] == var)), "p_value"]
      cs[min(which(cs[["var"]] == var)), "p_value"] <- p
      cs[max(which(cs[["var"]] == var)), "p_value"] <- NA
    }
  }
  for (i in 1:nrow(cs)) {
    if (is.na(cs[i, "terms"]) & !is.na(cs[i, "header"])) {
      cs[i, "terms"] <- cs[i, "var"]
    }
    else if (is.na(cs[i, "terms"])) {
      cs[i, "terms"] <- paste0(cs[i, "var"], cs[i, "lvl"])
    }
  }

  if (vif){
    VIF <- try(GVIF(model),silent = TRUE)
    names(VIF)[1] <- "terms"
    cs <- dplyr::full_join(cs,VIF, by = join_by(terms))
  }

  if (whichp!="level"){
    global_p <- gp(model)
    colnames(global_p) <- c("terms", "global_p")
    cs <- dplyr::full_join(cs,global_p, by = dplyr::join_by(terms))
  }
  for (i in 1:nrow(cs)) {
    if (!is.na(cs[i, "header"])) {
      v <- cs[i, "var"]
      n <- sum(cs[which(cs[, "var"] == v), "n"], na.rm = T)
      if ("Events" %in% colnames(cs)) {
        e <- sum(cs[which(cs[, "var"] == v), "Events"], na.rm = T)
        cs[i, c("n", "Events")] <- c(n, e)
      }
      else {
        cs[i, "n"] <- n
      }
    }
  }
  cs$ord <- 1:nrow(cs)
  attr(cs,'estLabel') <- attr(m_coeff,'estLabel')

  if (for_plot) {
    return(cs)
  }

  estLbl <- attr(cs, "estLabel")
  var_col <- c()
  for (i in 1:nrow(cs)) {
    if (!is.na(cs[i, "ref"]) & cs[i, "ref"] == TRUE) {
      cs[i, "Est_CI"] <- "Reference"
    }
    if (!is.na(cs[i, "header"])) {
      var_col <- c(var_col, cs[i, "var"])
    }
    else if (length(which(cs[, "var"] == cs[i, "var"])) == 1) {
      var_col <- c(var_col, cs[i, "var"])
    }
    else {
      var_col <- c(var_col, cs[i, "lvl"])
    }
  }
  cs <- cbind(data.frame(Variable = var_col), cs)

  # if (for_plot) {
  #   return(cs)
  # }


  if (whichp == "both") {
    for (i in 1:nrow(cs)) {
      if (!is.na(cs[i, "header"]) & (length(which(cs$var == cs[i, "var"])) > 3)) {
        cs[i, "p_value"] <- cs[i, "global_p"]
      }
    }
  }
  else if (whichp == "global") {
    cs[["p_value"]] <- cs[["global_p"]]
  }

  #
  #   lbl <- c()
  #   for (i in 1:nrow(mcoeff)) {
  #     if (is.na(mcoeff[i, "lvl"])) {
  #       lbl <- c()
  #     }
  #   }


  cols_to_keep <- c("Variable", "Est_CI", "p_value", "n")
  new_colnames <- c("Variable", estLbl, "p-value", "N")
  if ("Events" %in% colnames(cs)) {
    cols_to_keep <- c(cols_to_keep, "Events")
    new_colnames <- c(new_colnames, "Event")
  }
  if (vif) {
    cols_to_keep <- c(cols_to_keep, "VIF")
    new_colnames <- c(new_colnames, "VIF")
  }
  # print(cols_to_keep)
  # print(cs)
  # print(new_colnames)
  cs <- cs[, cols_to_keep]
  colnames(cs) <- new_colnames

  # #remove any columns that only have NA values
  # for (col in colnames(mcoeff)) {
  #   if (all(is.na(mcoeff[[col]]))) {
  #     mcoeff[[col]] <- NULL
  #   }
  # }


  # if (nicenames) {
  #   mcoeff[["Variable"]] <- replaceLbl(mcoeff, "Variable")
  # }




  return(cs)
}

# Extract model components ------------
coeffSum <- function(model,CIwidth=.95,digits=2) {
  CIwidth=CIwidth;digits=digits
  UseMethod("coeffSum",model)
}

coeffSum.lme <- function(model,CIwidth=.95,digits=2) {
  ms <- data.frame(summary(model)$tTable)
  pt <- 1-(1-CIwidth)/2
  ci <- data.frame(terms=rownames(ms),
                   lwr=ms$Value-qt(pt,df=ms$DF)*ms$Std.Error,
                   upr=ms$Value+qt(pt,df=ms$DF)*ms$Std.Error)
  cs <- data.frame(
    terms=rownames(ms),
    est=ms$Value,
    p_value = ms$p.value
  )
  cs <- merge(cs,ci,all.x = T)
  cs$Est_CI <- apply(cs[,c('est','lwr','upr')],MARGIN = 1,function(x) psthr(x,digits))
  attr(cs,'estLabel') <- betaWithCI("Estimate",CIwidth)
  return(cs)
}

coeffSum.lmerMod <- function(model,CIwidth=.95,digits=2) {
  stop("Method not implemented for lmer fit from lme4,\nre-fit model using lmeTest package.")
}

coeffSum.lmerModLmerTest <- function(model,CIwidth=.95,digits=2) {
  ms <- data.frame(summary(model)$coefficients)
  pt <- 1-(1-CIwidth)/2
  df <- sw_df(model)
  ms$t_value <- ms$Estimate/ms$Std..Error
  ms$p_value <- 2*pt(abs(ms$t_value),df,lower.tail=FALSE)
  ci <- data.frame(terms=rownames(ms),
                   lwr=ms$Estimate-qt(pt,df)*ms$Std..Error,
                   upr=ms$Estimate+qt(pt,df)*ms$Std..Error)
  cs <- data.frame(
    terms=rownames(ms),
    est=ms$Estimate,
    p_value= ms$p_value
  )
  cs <- merge(cs,ci,all.x = T)
  cs$Est_CI <- apply(cs[,c('est','lwr','upr')],MARGIN = 1,function(x) psthr(x,digits))
  attr(cs,'estLabel') <- betaWithCI("Estimate",CIwidth)
  return(cs)
}

# this code taken from the lmerTest package to compute the Satterthwaite df
sw_df <- function(model){
  qform <- function(x,A){
    sum(x * (A %*% x))
  }
  L <- diag(length(fixef(model)))
  var_con <- qform(L, vcov(model))
  grad_var_con <- vapply(model@Jac_list, function(x) qform(L,x), numeric(1L))
  satt_denom <- qform(grad_var_con, model@vcov_varpar)
  df <- drop(2 * var_con^2/satt_denom)
  return(df)
}

coeffSum.default <- function(model,CIwidth=.95,digits=2) {
  ms <- summary(model)$coefficients
  ci <- as.data.frame(confint(model,level = CIwidth))
  names(ci) <- c("lwr","upr")
  ci$terms <- rownames(ci)
  rownames(ci) <- NULL
  cs <- data.frame(
    terms=rownames(ms),
    est=ms[,1],
    p_value = ms[,4]
  )
  cs <- merge(cs,ci,all.x = T)
  cs$Est_CI <- apply(cs[,c('est','lwr','upr')],MARGIN = 1,function(x) psthr(x,digits))
  attr(cs,'estLabel') <- betaWithCI("Estimate",CIwidth)
  return(cs)
}

coeffSum.crrRx <- function(model,CIwidth=.95,digits=2) {
  ms <- data.frame(model$coeffTbl)
  ci <- try(exp(confint(model,level = CIwidth)),silent = T)
  if (!inherits(ci,"try-error")){
    if (!inherits(ci,"matrix")) {
      ci <- matrix(ci,ncol=2)
      rownames(ci) <- rownames(ms)[1]
    }
    ci <- data.frame(ci)
    names(ci) <- c("lwr","upr")
    ci$terms <- rownames(ci)
    rownames(ci) <- NULL
  }   else {
    Z_mult = qnorm(1 - (1 - CIwidth)/2)
    ci <- data.frame(lwr=exp(ms[, 1] - Z_mult * ms[, 3]),
                     upr=exp(ms[, 1] + Z_mult * ms[, 3]))
    ci$terms <- rownames(ms)
  }
  cs <- data.frame(
    terms=rownames(ms),
    est=ms[,2],
    p_value = ms[,5]
  )
  cs <- merge(cs,ci,all.x = T)
  cs$Est_CI <- apply(cs[,c('est','lwr','upr')],MARGIN = 1,function(x) psthr(x,digits))
  attr(cs,'estLabel') <- betaWithCI("HR",CIwidth)
  return(cs)
}

coeffSum.coxph <- function(model,CIwidth=.95,digits=2) {
  ms <- summary(model)$coefficients
  ci <- exp(as.data.frame(confint(model,level = CIwidth)))
  names(ci) <- c("lwr","upr")
  ci$terms <- rownames(ci)
  rownames(ci) <- NULL
  cs <- data.frame(
    terms=rownames(ms),
    est=ms[,2],
    p_value = ms[,5]
  )
  cs <- merge(cs,ci,all.x = T)
  cs$Est_CI <- apply(cs[,c('est','lwr','upr')],MARGIN = 1,function(x) psthr(x,digits))
  attr(cs,'estLabel') <- betaWithCI("HR",CIwidth)
  return(cs)
}

coeffSum.glm <- function(model,CIwidth=.95,digits=2) {
  ms <- data.frame(summary(model)$coefficients)
  if (model$family$link %in% c("logit","log")){
    ci <- as.data.frame(exp(confint(model,level = CIwidth)))
    ms$Estimate <- exp(ms$Estimate)
  } else {
    ci <- as.data.frame(confint(model,level = CIwidth))
  }
  names(ci) <- c("lwr","upr")
  ci$terms <- rownames(ci)
  rownames(ci) <- NULL

  cs <- data.frame(
    terms=rownames(ms),
    est=ms[,1],
    p_value = ms[,4]
  )

  cs <- merge(cs,ci,all.x = T)
  cs$Est_CI <- apply(cs[,c('est','lwr','upr')],MARGIN = 1,function(x) psthr(x,digits))
  if (model$family$link == "logit"){
    attr(cs,'estLabel') <- betaWithCI("OR",CIwidth)
  } else if (model$family$link == "log"){
    attr(cs,'estLabel') <- betaWithCI("RR",CIwidth)
  } else {
    attr(cs,'estLabel') <- betaWithCI("Estimate",CIwidth)
  }
  return(cs)
}

coeffSum.polr <- function(model,CIwidth=.95,digits=2) {
  ms <- summary(model)$coefficients
  ci <- try(exp(confint(model,level = CIwidth)),silent = T)
  if (!inherits(ci,"try-error")){
    if (!inherits(ci,"matrix")) {
      ci <- matrix(ci,ncol=2)
      rownames(ci) <- rownames(ms)[1]
    }
    ci <- data.frame(ci)
    names(ci) <- c("lwr","upr")
    ci$terms <- rownames(ci)
    rownames(ci) <- NULL
  }   else {
    m <- summary(model)$coefficients
    Z_mult = qnorm(1 - (1 - CIwidth)/2)
    ci <- data.frame(lwr=exp(m[, 1] - Z_mult * m[, 2]),
                     upr=exp(m[, 1] + Z_mult * m[, 2]))
    ci$terms <- rownames(ci)
  }
  pvalues = stats::pnorm(abs(ms[, "Value"]/ms[,"Std. Error"]), lower.tail = FALSE) * 2

  cs <- data.frame(
    terms=rownames(ms),
    est=ms[,2],
    p_value = pvalues
  )

  cs <- merge(cs,ci,all.x = T)
  cs <- cs[!grepl("[|]",cs$terms),]
  cs$Est_CI <- apply(cs[,c('est','lwr','upr')],MARGIN = 1,function(x) psthr(x,digits))
  attr(cs,'estLabel') <- betaWithCI("OR",CIwidth)
  return(cs)
}


# Extract data from a fitted model ---------------
get_model_data <- function(model){
  UseMethod("get_model_data", model)
}

get_model_data.default <- function(model){
  return(NULL)
}
get_model_data.lm <- function(model){
  return(model$model)
}
get_model_data.lme <- function(model){
  return(model$data)
}
get_model_data.lmerMod <- function(model){
  return(model.frame(model))
}
get_model_data.lmerModLmerTest <- function(model){
  return(model.frame(model))
}
get_model_data.crrRx <- function(model){
  return(model$model)
}
get_model_data.polr <- function(model){
  return(model$model)
}
get_model_data.coxph <- function(model){
  if (is.null(model$data)){
  df <- try(stats::model.frame(model$call$formula,
                               eval(parse(text = paste("data=",
                                                       deparse(model$call$data))))), silent = TRUE)
  } else (df <- model$data)
  if (inherits(df,'try-error')) {
    warning ("Model data could not be extracted")
    return(NULL)
  }
  return(df)
}

# may need to add other methods

# Extract event counts from a fitted model ---------------
get_event_counts <- function(model){
  UseMethod("get_event_counts",model)
}

get_event_counts.default <- function(model){
  return(NULL)
}
get_event_counts.coxph <- function(model){
  md <- get_model_data(model)
  y <- md[[1]]
  if (ncol(y)==2) return(y[,2])
  if (any(grepl("[+]",y))){
    st <- ifelse(grepl("[+]",y),0,1)
    return(st)
  }
  return(NULL)
}
get_event_counts.crrRx <- function(model){
  md <- get_model_data(model)
  if (is.null(md)) return(NULL)
  return(md[[2]])
}
get_event_counts.glm <- function(model){
  if (model$family$family=="binomial"|model$family$family=="quasibinomial"){
    return(model$y)
  }
}
# Calculate a global p-value for categorical variables --------
gp <- function(model) {
  UseMethod("gp", model)
}
gp.default <- function(model,CIwidth=.95,digits=2) { # lm, negbin
  terms <- attr(model$terms, "term.labels")
  globalpvalue <- drop1(model,scope=terms,test = "Chisq")
  gp <- data.frame(var=rownames(globalpvalue)[-1],
                   global_p = globalpvalue[-1,5])
  attr(gp,"global_p") <-"LRT"
  return(gp)
}


gp.coxph <- function(model,CIwidth=.95,digits=2) {
  terms <- attr(model$terms, "term.labels")
  globalpvalue <- drop1(model,scope=terms,test="Chisq")
  gp <- data.frame(var=terms,
                   global_p = globalpvalue[["Pr(>Chi)"]][-1])
  attr(gp,"global_p") <-"LRT"
  return(gp)
}

gp.crrRx <- function(model,CIwidth=.95,digits=2) {
  terms <- strsplit(trimws(gsub(".*~","", deparse(model$call[[1]]))),"[+]")[[1]]
  terms <- sapply(terms,trimws)
  gp_vals <- data.frame(var=terms,
                        global_p = NA)
  rownames(gp_vals) <- NULL
  if (length(terms)>1){
  for (t in terms){
    x <- ifelse(length(setdiff(terms,t))>0,setdiff(terms,t),1)
    eval(parse(text = paste('m2 <-try(crrRx(',paste(paste(setdiff(names(model$model),terms),collapse = "+"),
                                                    "~", x, sep = ""),
                            ',data = model$model))')))

    if (!inherits(m2,"try-error")) {
      degf <- length(grep(t,names(model$coef)))
      gp <- pchisq(2*(model$loglik-m2$loglik),degf)
    } else gp <- NA
    gp_vals$global_p[which(gp_vals$var==t)] <- gp
    attr(gp_vals,"global_p") <-"LRT"
  }}
  return(gp_vals)
}

gp.glm <- function(model,CIwidth=.95,digits=2) {
  terms <- attr(model$terms, "term.labels")
  globalpvalue <- drop1(model,scope=terms,test="LRT")
  gp <- data.frame(var=rownames(globalpvalue)[-1],
                   global_p = globalpvalue[-1,5])
  attr(gp,"global_p") <-"LRT"
  return(gp)
}
gp.lme <- function(model,CIwidth=.95,digits=2) {
  terms <- attr(model$terms, "term.labels")
  globalpvalue <- drop1(update(model,method="ML"),scope=terms,test = "Chisq")
  gp <- data.frame(var=rownames(globalpvalue)[-1],
                   global_p = globalpvalue[-1,5])
  attr(gp,"global_p") <-"LRT"
  return(gp)
}
gp.lmerMod <- function(model,CIwidth=.95,digits=2) {
  terms <- attr(terms(model), "term.labels")
  globalpvalue <- drop1(update(model),scope=terms,test = "Chisq")
  gp <- data.frame(var=rownames(globalpvalue)[-1],
                   global_p = globalpvalue$`Pr(Chi)`[-1])
  attr(gp,"global_p") <-"LRT"
  return(gp)
}
gp.lmerModLmerTest <- function(model,CIwidth=.95,digits=2) {
  terms <- attr(terms(model), "term.labels")
  globalpvalue <- drop1(update(model),scope=terms,test = "Chisq")
  gp <- data.frame(var=rownames(globalpvalue),
                   global_p = globalpvalue$`Pr(>F)`)
  attr(gp,"global_p") <-"LRT"
  return(gp)
}

gp.polr <- function(model,CIwidth=.95,digits=2) {
  terms <- attr(model$terms, "term.labels")
  globalpvalue <- drop1(model,scope=terms,test="Chisq")
  gp <- data.frame(var=rownames(globalpvalue)[-1],
                   global_p = globalpvalue[["Pr(>Chi)"]][-1])
  attr(gp,"global_p") <-"LRT"
  return(gp)
}

gp.gee <- function(model,CIwidth=.95,digits=2) {
  terms <- attr(model$terms, "term.labels")
  terms <- sapply(terms,trimws)
  gp_vals <- data.frame(var=terms,
                        global_p = NA)
  rownames(gp_vals) <- NULL
  for (t in terms){
    covariateindex <- grep(paste0("^",t),names(model$coefficients))
    gp <- try(aod::wald.test(b = model$coefficients[covariateindex],
                             Sigma = (model$geese$vbeta)[covariateindex, covariateindex],
                             Terms = seq_len(length(model$coefficients[covariateindex])))$result$chi2[3],
              silent = T)
    if (inherits(gp,"try-error")) gp <- NA
    gp_vals$global_p[which(gp_vals$var==t)] <- gp
  }
  attr(gp_vals,"global_p") <-"Wald test"
  return(gp_vals)
}
