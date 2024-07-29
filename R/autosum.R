

model.summary <- function(model,digits=2,CIwidth = 0.95, whichp = FALSE, ...){
  mcoeff <- coeffSum(model,CIwidth,digits)
  # check units - if any lwr==upr issue warning
  if (any(mcoeff$lwr==mcoeff$upr)) message("Zero-width confidence interval detected. Check predictor units.")

  # Basically, we want to take mcoeff and
  terms <- attr(model$terms, "term.labels")
  if (all(mcoeff$Term %in% terms)){
    # No categorical variables, no need to add any reference data
  } else {
    # Categorical variables
    mcoeff$order <- 1:nrow(mcoeff) # keep track of the order of the coefficients
    cat_vars <- NULL
    for (catVar in setdiff(terms,mcoeff$Term)){
      vpos <- sapply(paste0(catVar,model$xlevels[[catVar]]),function(x) {
        p = which(mcoeff$Term==x)
        if (length(p)==0) p = NA
        return(p)})
      vpos <- vpos[!is.na(vpos)]
      reg_lvls <- gsub(catVar,"",mcoeff$Term[vpos])
      ref_lvl <- setdiff(model$xlevels[[catVar]],reg_lvls)

      for (i in 2:ncol(df)) data.frame(table(df[[1]],df[,i]))

      catInfo <- data.frame(Term = mcoeff$Term[vpos],
                            Variable= catVar,
                            var_level=reg_lvls,
                            ref_level=ref_lvl)
      cat_vars <- dplyr::bind_rows(cat_vars,catInfo)
    }
    # add to the data frame
    mcoeff <- merge(mcoeff,cat_vars,all.x=T)
  }
  mcoeff$Variable[is.na(mcoeff$Variable)] <- mcoeff$Term[is.na(mcoeff$Variable)]


  tpos <- as.numeric(factor(mcoeff$Variable))
  vars <- sapply(tpos,function(x) terms[x])
  # This needs to calculate "global" p-values for categorical variables
  # works for linear models - need to test all the others!
  # it adds the global-p-value to the dataframe
  if (!all(mcoeff$Term==vars)){
    mcoeff$variable <- vars
    gobal_p <- gp(model)
    mcoeff <- merge(mcoeff,gobal_p,all.x = TRUE,sort=FALSE)
  }
  mcoeff <- mcoeff[order(mcoeff$order),]
  # return(mcoeff)
  # From here - need to add the reference levels as rows to the table
  for (var in unique(mcoeff$variable)) {
    ref_row <- data.frame(var_level = mcoeff$ref_level[1], Est_CI = "Reference", Term = paste0(mcoeff$Variable[1], mcoeff$ref_level[1]), variable = mcoeff$variable[1], Variable = mcoeff$Variable[1], ref_level = mcoeff$ref_level[1])
    mcoeff <- dplyr::bind_rows(ref_row, mcoeff)
    subset_var <- subset(mcoeff, variable == var)
    if (length(unique(subset_var$Term)) > 2) {
      if (whichp == "global" | whichp == "both") {
        var_row <- data.frame(var_level = mcoeff$Variable[1], p_value = mcoeff$global_p[1])
      }
      else { # whichp == "level"
        var_row <- data.frame(var_level = mcoeff$Variable[1])
      }
    }
    else { # two or less levels

      non_ref <- setdiff(subset_var$var_level, mcoeff$ref_level)
      # print(non_ref)
      # print(which(mcoeff$var_level == non_ref))

      mcoeff[which(mcoeff$var_level == non_ref), "p_value"] <- NA
      var_row <- data.frame(var_level = mcoeff$Variable[1], p_value = as.numeric(mcoeff$global_p))
    }
    mcoeff <- dplyr::bind_rows(var_row, mcoeff)
  }


  # reference_row <- c(mcoeff$ref_level[1], "Reference", rep(NA, ncol(mcoeff) - 2))
  # mcoeff <- dplyr::bind_rows(var_row, ref_row, mcoeff)
  # mcoeff <- select(mcoeff, var_level, Est_CI, p_value, est, lwr, upr, global_p)
  View(mcoeff)
}

# Combind model components and variable levels and sample sizes

m_summary <- function(model,CIwidth=.95,digits=2,vif = FALSE,whichp="level"){
  m_coeff <- coeffSum(model,CIwidth,digits)
  lvls <- getVarLevels(model)
  lvls$ord  <- 1:nrow(lvls)

  cs <- merge(lvls,m_coeff,all.x = T)
  if (vif){
    VIF <- try(GVIF(model),silent = TRUE)
    names(VIF)[1] <- "var"
    cs <- full_join(cs,VIF)
  }

  if (whichp!="level"){
    global_p <- gp(model)
    #print(global_p)
    cs <- merge(cs,global_p,all = T)
  }
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
  return(cs)
}

# Extract model components ------------
coeffSum <- function(model,CIwidth=.95,digits=2) {
  CIwidth=CIwidth;digits=digits
  UseMethod("coeffSum",model)
}

coeffSum.default <- function(model,CIwidth=.95,digits=2) {
  ms <- summary(model)$coefficients
  ci <- as.data.frame(confint(model,level = CIwidth))
  names(ci) <- c("lwr","upr")
  ci$terms <- rownames(ci)
  rownames(ci) <- NULL
  cs <- data.frame(
    terms=rownames(ms),
    est=ms[,2],
    p_value = ms[,4]
  )
  cs <- merge(cs,ci,all.x = T)
  cs$Est_CI <- apply(cs[,c('est','lwr','upr')],MARGIN = 1,function(x) psthr(x,digits))
  if (inherits(model,"coxph")){
    attr(cs,'estLabel') <- betaWithCI("HR",CIwidth)
  } else {
    attr(cs,'estLabel') <- betaWithCI("Estimate",CIwidth)
  }

  return(cs)
}

coeffSum.crrRx <- function(model,CIwidth=.95,digits=2) {
  ms <- model$coeffTbl
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
    ci$terms <- rownames(ci)
  }
  cs <- data.frame(
    terms=rownames(ms),
    est=ms[,2],
    p_value = ms[,4]
  )
  cs <- merge(cs,ci,all.x = T)
  cs$Est_CI <- apply(cs[,c('est','lwr','upr')],MARGIN = 1,function(x) psthr(x,digits))
  attr(cs,'estLabel') <- betaWithCI("HR",CIwidth)
  return(cs)
}

coeffSum.glm <- function(model,CIwidth=.95,digits=2) {
  ms <- summary(model)$coefficients
  if (model$family$link %in% c("logit","log")){
    ci <- as.data.frame(exp(confint(model,level = CIwidth)))
  } else {
    ci <- as.data.frame(confint(model,level = CIwidth))
  }
  names(ci) <- c("lwr","upr")
  ci$terms <- rownames(ci)
  rownames(ci) <- NULL

  cs <- data.frame(
    terms=rownames(ms),
    est=ms[,2],
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
get_model_data.crrRx <- function(model){
  return(model$model)
}
get_model_data.polr <- function(model){
  return(model$model)
}
get_model_data.coxph <- function(model){
  df <- try(stats::model.frame(model$call$formula,
                               eval(parse(text = paste("data=",
                                                       deparse(model$call$data))))), silent = TRUE)
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

gp.polr <- function(model,CIwidth=.95,digits=2) {
  terms <- attr(model$terms, "term.labels")
  globalpvalue <- drop1(model,scope=terms,test="ChiSq")
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
