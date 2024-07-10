#
#
# model.summary <- function(model,digits=2,CIwidth = 0.95, ...){
#   mcoeff <- coeffSum(model,CIwidth,digits)
#   # check units - if any lwr==upr issue warning
#   if (any(mcoeff$lwr==mcoeff$upr)) message("Zero-width confidence interval detected. Check predictor units.")
#
#   terms <- attr(model$terms, "term.labels")
#   if (all(mcoeff$Term %in% terms)){
#     # No categorical variables
#   } else {
#     # Categorical variables
#     for (catVar in setdiff(terms,mcoeff$Term)){
#       vpos <- sapply(paste0(catVar,model$xlevels[[catVar]]),function(x) {
#         p = which(mcoeff$Term==x)
#         if (length(p)==0) p = NA
#         return(p)})
#       vpos <- vpos[!is.na(vpos)]
#       reg_lvls <- gsub(catVar,"",mcoeff$Term[vpos])
#       ref_lvl <- setdiff(model$xlevels[[catVar]],reg_lvls)
#
#       catInfo <- data.frame(Term = mcoeff$Term[vpos],
#                             Variable= catVar,
#                             var_level=reg_lvls,
#                             ref_level=ref_lvl)
#     }
#
#   }
#   tpos <- model$assign[-1]
#   vars <- sapply(tpos,function(x) terms[x])
#   if (!all(mcoeff$Term==vars)){
#     mcoeff$variable <- vars
#     drop_p <- drop1(model,scope=terms,test = "Chisq")
#     gp <- data.frame(variable=rownames(drop_p)[-1],
#                      global_p = drop_p[-1,5])
#     mcoeff$order <- 1:nrow(mcoeff)
#     return <- merge(mcoeff,gp,all.x = TRUE,sort=FALSE)
#     return <- return[order(return$order),]
#     return$level <- mapply(function(term,var){
#       sub(var,'',term)
#     },return$Term,return$variable)
#     for (term in terms) {
#       reg_lvls <- gsub(term,"",names(model$coefficients))[varInd]
#       ref_lvl <- setdiff(unique(model$model[[term]]),reg_lvls)
#       print(ref_lvl)
#     }
#
#   }
# }

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
gp <- function(model, reduced_model, ...) {
  UseMethod("gp", model,reduced_model,...)
}
gp.default <- function(model,CIwidth=.95,digits=2,...) { # lm, negbin
  terms <- attr(model$terms, "term.labels")
  globalpvalue <- drop1(model,scope=terms,test = "Chisq")
}
gp.glm <- function(model,CIwidth=.95,digits=2,...) {
  terms <- attr(model$terms, "term.labels")
  globalpvalue <- drop1(model,scope=terms,test="LRT")
}
gp.lme <- function(model,CIwidth=.95,digits=2,...) {
  terms <- attr(model$terms, "term.labels")
  globalpvalue <- drop1(update(model,method="ML"),scope=terms,test = "Chisq")
}
gp.gee <- function(model,CIwidth=.95,digits=2,...) {
  globalpvalue <- try(aod::wald.test(b = model$coefficients[covariateindex],
                                     Sigma = (model$geese$vbeta)[covariateindex, covariateindex],
                                     Terms = seq_len(length(model$coefficients[covariateindex])))$result$chi2[3],
                      silent = T)
  if (inherits(globalpvalue,"try-error")) return(NULL)

}

wald_gp <- function(){
  globalpvalue <- try(aod::wald.test(b = model$coef$fixed[covariateindex],
                                     Sigma = vcov(model)[covariateindex, covariateindex],
                                     Terms = seq_along(covariateindex))$result$chi2[3],silent = T)

}
anova_gp <- function(model,terms,...){
  gp_vals <- data.frame(Term=terms,global_p=NA)
  for (t in terms){
    reduced_model <-try( update(model,as.formula(paste(".~",paste(setdiff(terms,t),collapse = "+")))),silent = T)
    if (!inherits(reduced_model,"try-error")) {
      print(attr(reduced_model$terms, "term.labels"))
      gp <-try(stats::na.omit(anova(model,reduced_model)[,"Pr(Chi)"]))
    } else gp <- NA
    gp_vals$global_p[which(gp_vals$Term==t)] <- gp
  }
  return(gp_vals)
}
