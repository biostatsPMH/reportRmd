

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
    mcoeff <- bind_rows(ref_row, mcoeff)
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
    mcoeff <- bind_rows(var_row, mcoeff)
  }


  # reference_row <- c(mcoeff$ref_level[1], "Reference", rep(NA, ncol(mcoeff) - 2))
  # mcoeff <- bind_rows(var_row, ref_row, mcoeff)
  # mcoeff <- select(mcoeff, var_level, Est_CI, p_value, est, lwr, upr, global_p)
  View(mcoeff)
}

# Extract model components ------------
coeffSum <- function(model,CIwidth=.95,digits=2,...) {
  CIwidth=CIwidth;digits=digits
  UseMethod("coeffSum",model)
}

coeffSum.default <- function(model,CIwidth=.95,digits=2,...) {
  ms <- summary(model)$coefficients
  ci <- confint(model,level=CIwidth)
  cs <- data.frame(
    Term=rownames(ms),
    est=ms[,1],
    p_value = ms[,4],
    lwr=ci[,1],
    upr=ci[,2]
  )
  rownames(cs) <- NULL

  cs$Est_CI <- apply(cs[,c('est','lwr','upr')],MARGIN = 1,function(x) psthr(x,digits))
  attr(cs,'estLabel') <- betaWithCI("Estimate",CIwidth)
  return(cs[-1,])
}
coeffSum.geeglm <- function(model,CIwidth=.95,digits=2,...) {
  ms <- summary(model)$coefficients
  if (grepl("log",model$family$link)){
    estFun <- exp
    ci_mult <- stats::qnorm(1 - (1 - CIwidth)/2)
  } else {
    estFun <- identity
    ci_mult <- stats::qt(1 - (1 - CIwidth)/2,model$df.residual)
  }
  ci <- cbind(estFun(ms[,1]-ci_mult*ms[,2]),estFun(ms[,1]+ci_mult*ms[,2]))
  cs <- data.frame(
    Term=rownames(ms),
    est=estFun(ms[,1]),
    p_value = ms[,4],
    lwr=ci[,1],
    upr=ci[,2]
  )
  rownames(cs) <- NULL

  cs$Est_CI <- apply(cs[,c('est','lwr','upr')],MARGIN = 1,function(x) psthr(x,digits))
  beta <- ifelse(model$family$family=="gaussian","Estimate",
                 ifelse(model$family$family=="binomial","OR",
                        ifelse(model$family$family=="poisson","RR","GEE Estimate")))
  attr(cs,'estLabel') <- betaWithCI(beta,CIwidth)
  return(cs[-1,])
}

coeffSum.glm <- function(model,CIwidth=.95,digits=2,...) {
  ms <- summary(model)$coefficients
  ci <- confint(model,level=CIwidth)
  if (grepl("log",model$family$link)) estFun <- exp else estFun <- identity
  ci <- estFun(ci)
  cs <- data.frame(
    Term=rownames(ms),
    est=estFun(ms[,1]),
    p_value = ms[,4],
    lwr=ci[,1],
    upr=ci[,2]
  )
  rownames(cs) <- NULL
  cs$Est_CI <- apply(cs[,c('est','lwr','upr')],MARGIN = 1,function(x) psthr(x,digits))
  beta <- ifelse(model$family$family=="gaussian","Estimate",
                 ifelse(model$family$family=="binomial","OR",
                        ifelse(model$family$family=="poisson","RR","GLM Estimate")))
  attr(cs,'estLabel') <- betaWithCI(beta,CIwidth)
  return(cs[-1,])
}

coeffSum.negbin <- function(model,CIwidth=.95,digits=2,...) {
  ms <- summary(model)$coefficients
  ci <- exp(confint(model,level=CIwidth))
  cs <- data.frame(
    Term=rownames(ms),
    est=exp(ms[,1]),
    p_value = ms[,4],
    lwr=ci[,1],
    upr=ci[,2]
  )
  rownames(cs) <- NULL
  cs$Est_CI <- apply(cs[,c('est','lwr','upr')],MARGIN = 1,function(x) psthr(x,digits))
  attr(cs,'estLabel') <- betaWithCI("RR",CIwidth)
  return(cs[-1,])
}

coeffSum.coxph <- function(model,CIwidth=.95,digits=2,...) {
  ms <- summary(model)$coefficients
  ci <- exp(confint(model,level = CIwidth))
  cs <- data.frame(
    Term=rownames(ms),
    est=exp(ms[,1]),
    p_value = ms[,5],
    lwr = ci[,1],
    upr=ci[,2]
  )
  rownames(cs) <- NULL
  cs$Est_CI <- apply(cs[,c('est','lwr','upr')],MARGIN = 1,function(x) psthr(x,digits))
  attr(cs,'estLabel') <- betaWithCI("HR",CIwidth)
  return(cs)
}

coeffSum.crr <- function(model,CIwidth=.95,digits=2,...) {
  out <- summary(model, conf.int = CIwidth)
  ms <- out$coef
  ci <- out$conf.int
  cs <- data.frame(
    Term=rownames(ms),
    est=ms[,1],
    p_value = ms[,5],
    lwr = ci[,3],
    upr = ci[,4]
  )
  rownames(cs) <- NULL
  cs$Est_CI <- apply(cs[,c('est','lwr','upr')],MARGIN = 1,function(x) psthr(x,digits))
  attr(cs,'estLabel') <- betaWithCI("HR",CIwidth)
  return(cs)
}

coeffSum.lme <- function(model,CIwidth=.95,digits=2,...) {
  ms <- summary(model)$tTable
  t_mult <- qt(1 - (1 - CIwidth)/2,ms[,3])
  cs <- data.frame(
    Term=rownames(ms),
    est=ms[,1],
    p_value = ms[,5],
    lwr = ms[,1]-t_mult*ms[,2],
    upr = ms[,1]+t_mult*ms[,2]
  )
  rownames(cs) <- NULL
  cs$Est_CI <- apply(cs[,c('est','lwr','upr')],MARGIN = 1,function(x) psthr(x,digits))
  attr(cs,'estLabel') <- betaWithCI("Estimate",CIwidth)
  return(cs)
}

coeffSum.polr <- function(model,CIwidth=.95,digits=2,...) {
  ms <- summary(model)$coefficients
  ci <- matrix(exp(confint(model,level=CIwidth)),ncol = 2)
  cs <- data.frame(
    Term=rownames(ms),
    est=exp(ms[,1]),
    p_value = stats::pt(abs(ms[,3]),model$df.residual, lower.tail = FALSE)*2,
    lwr=ci[,1],
    upr=ci[,2]
  )
  rownames(cs) <- NULL
  cs$Est_CI <- apply(cs[,c('est','lwr','upr')],MARGIN = 1,function(x) psthr(x,digits))
  cs <- cs[cs$Term %in% names(model$coefficients),]
  attr(cs,'estLabel') <- betaWithCI("OR",CIwidth)
  return(cs)
}



# Calculate a global p-value for categorical variables --------
gp <- function(model) {
  UseMethod("gp", model)
}
gp.default <- function(model,CIwidth=.95,digits=2) { # lm, negbin
  terms <- attr(model$terms, "term.labels")
  globalpvalue <- drop1(model,scope=terms,test = "Chisq")
  gp <- data.frame(Variable=rownames(globalpvalue)[-1],
                   global_p = globalpvalue[-1,5])
  attr(gp,"global_p") <-"LRT"
  return(gp)
}


gp.coxph <- function(model,CIwidth=.95,digits=2) {
  terms <- attr(model$terms, "term.labels")
  globalpvalue <- drop1(model,scope=terms,test="Chisq")
  gp <- data.frame(Variable=terms,
                   global_p = globalpvalue[["Pr(>Chi)"]][-1])
  attr(gp,"global_p") <-"LRT"
  return(gp)
}

gp.crr <- function(model,CIwidth=.95,digits=2) {
  terms <- strsplit(trimws(gsub(".*~","", deparse(model$call[[1]]))),"[+]")[[1]]
  terms <- sapply(terms,trimws)
  gp_vals <- data.frame(Variable=terms,
                        global_p = NA)
  rownames(gp_vals) <- NULL
  for (t in terms){
    eval(parse(text = paste('m2 <-try(crrRx(',paste(paste(setdiff(names(model$data),terms),collapse = "+"),
                                                    "~", setdiff(terms,t), sep = ""),
                            ',data = model$data))')))

    if (!inherits(m2,"try-error")) {
      degf <- length(grep(t,names(model$coef)))
      gp <- pchisq(2*(model$loglik-m2$loglik),degf)
    } else gp <- NA
    gp_vals$global_p[which(gp_vals$Variable==t)] <- gp
    attr(gp_vals,"global_p") <-"LRT"
  }
  return(gp_vals)
}

gp.glm <- function(model,CIwidth=.95,digits=2) {
  terms <- attr(model$terms, "term.labels")
  globalpvalue <- drop1(model,scope=terms,test="LRT")
  gp <- data.frame(Variable=rownames(globalpvalue)[-1],
                   global_p = globalpvalue[-1,5])
  attr(gp,"global_p") <-"LRT"
  return(gp)
}
gp.lme <- function(model,CIwidth=.95,digits=2) {
  terms <- attr(model$terms, "term.labels")
  globalpvalue <- drop1(update(model,method="ML"),scope=terms,test = "Chisq")
  gp <- data.frame(Variable=rownames(globalpvalue)[-1],
                   global_p = globalpvalue[-1,5])
  attr(gp,"global_p") <-"LRT"
  return(gp)
}

gp.polr <- function(model,CIwidth=.95,digits=2) {
  terms <- attr(model$terms, "term.labels")
  globalpvalue <- drop1(model,scope=terms,test="ChiSq")
  gp <- data.frame(Variable=rownames(globalpvalue)[-1],
                   global_p = globalpvalue[["Pr(>Chi)"]][-1])
  attr(gp,"global_p") <-"LRT"
  return(gp)
}

gp.gee <- function(model,CIwidth=.95,digits=2) {
  terms <- attr(model$terms, "term.labels")
  terms <- sapply(terms,trimws)
  gp_vals <- data.frame(Variable=terms,
                        global_p = NA)
  rownames(gp_vals) <- NULL
  for (t in terms){
    covariateindex <- grep(paste0("^",t),names(model$coefficients))
    gp <- try(aod::wald.test(b = model$coefficients[covariateindex],
                             Sigma = (model$geese$vbeta)[covariateindex, covariateindex],
                             Terms = seq_len(length(model$coefficients[covariateindex])))$result$chi2[3],
              silent = T)
    if (inherits(gp,"try-error")) gp <- NA
    gp_vals$global_p[which(gp_vals$Variable==t)] <- gp
  }
  attr(gp_vals,"global_p") <-"Wald test"
  return(gp_vals)
}


