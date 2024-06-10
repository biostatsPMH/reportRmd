# Cramer's V CI
# References (need to place in wrapper function):
# Smithson, M. (2002). Noncentral CIwidthidence Intervals for Standardized Effect Sizes. In CIwidthidence Intervals (07/140 ed., Vol. 140). SAGE Publications. https://doi.org/10.4135/9781412983761.n4
# Steiger, J. H. (2004). Beyond the F Test: Effect Size Confidence Intervals and Tests of Close Fit in the Analysis of Variance and Contrast Analysis. Psychological Methods, 9(2), 164–182. https://doi.org/10.1037/1082-989X.9.2.164
# ref for epsilonSq: Kelley, T. L. (1935). An Unbiased Correlation Ratio Measure. Proceedings of the National Academy of Sciences - PNAS, 21(9), 554–559. https://doi.org/10.1073/pnas.21.9.554
# Ref to support omegaSq as best ANOVA effect size:
# Okada, K. (2013). Is Omega Squared Less Biased? A Comparison of Three Major Effect Size Indices in One-Way ANOVA. Behavior Research Methods, 40(2), 129-147.
# Ref for the KW
# Breslow, N. (1970). A generalized Kruskal-Wallis test for comparing K samples subject to unequal patterns of censorship. Biometrika, 57(3), 579-594.
# Ref for wilxon test
# FRITZ, C. O., MORRIS, P. E., & RICHLER, J. J. (2012). Effect Size Estimates: Current Use, Calculations, and Interpretation. Journal of Experimental Psychology. General, 141(1), 2–18. https://doi.org/10.1037/a0024338

# Website for ncp calculation and example:
# https://methods-sagepub-com.myaccess.library.utoronto.ca/book/confidence-intervals/n4.xml

# General Process:
# Identity the appropriate effect size transformation from the test statistic
# Cramer's V is from ChiSq,
#CIs noncentral chi-sq/t noncentrality parameter (ncp):

# DELTA (NCP) FUNCTIONS ----------------------------------------------------------------

# This function will return the NCP CI limits for the F, chi-sq and t distributions,
# for computing CIs of Cramer's V, Cohens d,  omega squared and epsilon squared for the Kruskal Wallis
# see beyond the F statistic for the rationale & algorithm

delta_CI <- function(htest,CIwidth=0.95){
  # determine what test result is being passed in
  if (inherits(htest,"aov")){
    class(htest) <- c(class(htest),"rm_F")
  } else if (inherits(htest,"htest")){
    if (grepl("Chi",htest$method)){
      htest <- uncorrectedChi(htest)
      class(htest) <- c(class(htest),"rm_chi")
    } else if (grepl("Kruskal",htest$method)) {
      class(htest) <- c(class(htest),"rm_chi")
    } else if (grepl("Wilcoxon rank sum",htest$method)) {
      class(htest) <- c(class(htest),"rm_chi")
    } else if (grepl("t-test",htest$method)) {
      class(htest) <- c(class(htest),"rm_t")
    } else if (grepl("Fisher",htest$method)) {
      class(htest) <- c(class(htest),"rm_chi")
    }}  else stop("htest must be the output of a call to t.test, chisq.test, kruskal.test, aov or fisher.test.rm")
  dL <- delta_l(htest,CIwidth)
  dU <- delta_u(htest,CIwidth)
  ci <- c(dL,dU)
  attr(ci,"confidence") <- CIwidth
  return(ci)
}

delta_l <- function(htest,CIwidth){
  UseMethod("delta_l")
}

delta_l.default <- function(...){
  stop("No default method for delta_l. The test type must be known")
}

delta_l.rm_F <- function(htest,CIwidth) {
  hsum <- summary(htest)[[1]]
  Fstat <- hsum[["F value"]][1]
  nu1 <- hsum[["Df"]][1]
  nu2 <- hsum[["Df"]][2]
  ulim <- 1 - (1-CIwidth)/2
  if (stats::pf(Fstat,nu1,nu2) < ulim) return(0)

  #	This first part finds a lower and upper value from which to start.
  lc <- c(.001,Fstat/2,Fstat)
  #	This first part finds a lower value from which to start.
  while(stats::pf(Fstat,nu1,nu2,lc[1])<ulim) {
    lc <- c(lc[1]/4,lc[1],lc[3])
  }
  while(stats::pf(Fstat,nu1,nu2,lc[3])>ulim) {
    lc <- c(lc[1],lc[3],lc[3]+Fstat)
  }
  #	This next part finds the lower limit for the ncp.
  diff <- 1
  while(diff > .00001) {
    if(stats::pf(Fstat,nu1,nu2,lc[2])<ulim){
      lc <- c(lc[1],(lc[1]+lc[2])/2,lc[2])
    }  else lc <- c(lc[2],(lc[2]+lc[3])/2,lc[3])
    diff <- abs(stats::pf(Fstat,nu1,nu2,lc[2]) - ulim)
  }
  names(lc) <- NULL
  return(lc[2])
}

delta_l.rm_chi <- function(htest,CIwidth) {
  chival <- htest$statistic; df <- htest$parameter
  ulim <- 1 - (1-CIwidth)/2
  if (stats::pchisq(chival,df)< ulim) return(0)

  #	This first part finds a lower value from which to start.
  lc <- c(.001,chival/2,chival)
  while(stats::pchisq(chival,df,lc[1])<ulim) {
    if(stats::pchisq(chival,df)<ulim)
      return(0)
    lc <- c(lc[1]/4,lc[1],lc[3])
  }
  #	This next part finds the lower limit for the ncp.
  diff <- 1
  while(diff > .00001) {
    if(stats::pchisq(chival,df,lc[2])<ulim){
      lc <- c(lc[1],(lc[1]+lc[2])/2,lc[2])
    }  else lc <- c(lc[2],(lc[2]+lc[3])/2,lc[3])
    diff <- abs(stats::pchisq(chival,df,lc[2]) - ulim)
  }
  names(lc) <- NULL
  return(lc[2])
}

delta_l.rm_t <- function(htest,CIwidth) {
  tval <- abs(htest$statistic); df <- htest$parameter
  ulim <- 1 - (1-CIwidth)/2
  if (stats::pt(tval,df) < ulim) return(0)

  #	This first part finds a lower and upper value from which to start.
  lc <- c(-tval,tval/2,tval)
  while(stats::pt(tval,df,lc[1])<ulim) {
    lc <- c(lc[1]-tval,lc[1],lc[3])
  }
  #	This next part finds the lower limit for the ncp.
  diff <- 1
  while(diff > .00001) {
    if(stats::pt(tval,df,lc[2])<ulim)
      lc <- c(lc[1],(lc[1]+lc[2])/2,lc[2])
    else lc <- c(lc[2],(lc[2]+lc[3])/2,lc[3])
    diff <- abs(stats::pt(tval,df,lc[2]) - ulim)
  }
  names(lc) <- NULL
  return(lc[2])
}

delta_u <- function(htest,CIwidth){
  UseMethod("delta_u")
}

delta_u.default <- function(...){
  stop("No default method for delta_u. The test type must be known")
}

delta_u.rm_F <- function(htest,CIwidth) {
  hsum <- summary(htest)[[1]]
  Fstat <- hsum[["F value"]][1]
  nu1 <- hsum[["Df"]][1]
  nu2 <- hsum[["Df"]][2]
  alpha <- 1-CIwidth
  if (stats::pf(Fstat,nu1,nu2)<alpha/2) return(0)
  #	This first part finds upper and lower starting values.
  uc <- c(Fstat,2*Fstat,3*Fstat)
  llim <- (1-CIwidth)/2
  while(stats::pf(Fstat,nu1,nu2,uc[1])<llim) {
    uc <- c(uc[1]/4,uc[1],uc[3])
  }
  while(stats::pf(Fstat,nu1,nu2,uc[3])>llim) {
    uc <- c(uc[1],uc[3],uc[3]+Fstat)
  }
  #	This next part finds the upper limit for the ncp.
  diff <- 1
  while(diff > .00001) {
    if(stats::pf(Fstat,nu1,nu2,uc[2])<llim)
      uc <- c(uc[1],(uc[1]+uc[2])/2,uc[2])
    else uc <- c(uc[2],(uc[2]+uc[3])/2,uc[3])
    diff <- abs(stats::pf(Fstat,nu1,nu2,uc[2]) - llim)
    lcdf <- stats::pf(Fstat,nu1,nu2,uc[2])
  }
  names(uc) <- NULL
  return(uc[2])
}

delta_u.rm_chi <- function(htest,CIwidth) {
  alpha <- 1-CIwidth
  chival <- htest$statistic; df <- htest$parameter
  if (stats::pchisq(chival,df)<alpha/2) return(0)
  #	This first part finds upper and lower starting values.
  uc <- c(chival,2*chival,3*chival)
  llim <- (1-CIwidth)/2
  while(stats::pchisq(chival,df,uc[1])<llim) {
    uc <- c(uc[1]/4,uc[1],uc[3])
  }
  while(stats::pchisq(chival,df,uc[3])>llim) {
    uc <- c(uc[1],uc[3],uc[3]+chival)
  }
  #	This next part finds the upper limit for the ncp.
  diff <- 1
  while(diff > .00001) {
    if(stats::pchisq(chival,df,uc[2])<llim)
      uc <- c(uc[1],(uc[1]+uc[2])/2,uc[2])
    else uc <- c(uc[2],(uc[2]+uc[3])/2,uc[3])
    diff <- abs(pchisq(chival,df,uc[2]) - llim)
    lcdf <- stats::pchisq(chival,df,uc[2])
  }
  names(uc) <- NULL
  return(uc[2])
}

delta_u.rm_t <- function(htest,CIwidth) {
  alpha <- 1-CIwidth
  tval <- abs(htest$statistic); df <- htest$parameter
  if (stats::pt(tval,df)<alpha/2) return(0)
  llim <- (1-CIwidth)/2
  #	This first part finds an upper value from which to start.
  uc <- c(tval,1.5*tval,2*tval)
  while(stats::pt(tval,df,uc[3])>llim) {
    uc <- c(uc[1],uc[3],uc[3]+tval)
  }
  #	This next part finds the upper limit for the ncp.
  diff <- 1
  while(diff > .00001) {
    if(stats::pt(tval,df,uc[2])<llim)
      uc <- c(uc[1],(uc[1]+uc[2])/2,uc[2])
    else uc <- c(uc[2],(uc[2]+uc[3])/2,uc[3])
    diff <- abs(stats::pt(tval,df,uc[2]) - llim)
    lcdf <- stats::pt(tval,df,uc[2])
  }
  names(uc) <- NULL
  return(uc[2])
}

# HELPER FUNCTIONS ----------------------------------------------------------------
# Function to make sure that the chi-square was not run with continuity correction
uncorrectedChi <- function(x) {
  stopifnot(inherits(x, "htest") )
  stopifnot("X-squared" %in% names(x[["statistic"]]))
  return(stats::chisq.test(x$observed,correct = FALSE))
}

# Function to preserve components of the data for effect size estimates
# Note that this can't be called as a formula - it needs xvar and group to be
# passed as separate variables

# Function to add the observed counts to the output of the fisher.test

fisher.test.rm <- function(x,...){
  chi.out <- suppressWarnings(stats::chisq.test(x,...))
  rtn <- stats::fisher.test(x,...)
  rtn$observed <- chi.out$observed
  rtn$statistic <- chi.out$statistic
  rtn$parameter <- chi.out$parameter
  rtn$k <- min(dim(chi.out$observed), na.rm = TRUE)
  rtn$N <- sum(chi.out$observed)
  rtn$method <- paste(rtn$method,"with additional chisq.test arguments")
  class(rtn) <- c(class(rtn),"fisher.test.rm")
  return(rtn)
}

chi.test.rm <- function(x,...){
  chi_test <- suppressWarnings(stats::chisq.test(x,...))
  chi_test$k <- min(dim(chi_test$observed), na.rm = TRUE)
  chi_test$N <- sum(chi_test$observed)
  class(chi_test) <- c(class(chi_test),"chi.test.rm")
  return(chi_test)
}

t.test.rm <- function(xvar,grp){
  t_test <- stats::t.test(xvar~grp)
  n <- unlist(lapply(levels(grp),function(g){
    length(na.omit(xvar[grp==g]))
  }))
  names(n) <- levels(grp)
  t_test$n <- n
  class(t_test) <- c(class(t_test),"t.test.rm")
  return(t_test)
}

# But we don't know the distribution of the resulting effect size under the
# alternative hypothesis

wilcox.test.rm <- function(xvar,grp){
  xg <- na.omit(cbind(xvar,grp))
  x <- xg[xg[,2]==1,1]
  y <- xg[xg[,2]==2,1]
  wilcox_test <- suppressWarnings(stats::wilcox.test(x,y))
  Ua <- wilcox_test$statistic
  Ub <- suppressWarnings(stats::wilcox.test(y,x)$statistic)
  U1 = min(Ua,Ub)
  U2 = max(Ua,Ub)
  n1 <- length(na.omit(x))
  n2 <- length(na.omit(y))
  # na <- length(na.omit(x))
  # nb <- length(na.omit(y))
  # mu <- na*nb/2
  # su = sqrt(na*nb*(na+nb+1)/12)
  # Z = abs(min(U1,U2)-mu)/su
  wilcox_test$U=min(U1,U2)
  Z = (U1+.5-(U1+U2)/2)/sqrt((n1*n2*(n1+n2+1)/12))
  n <- unlist(lapply(levels(grp),function(g){
    length(na.omit(xvar[grp==g]))
  }))
  names(n) <- levels(grp)
  wilcox_test$n <- n
  wilcox_test$Z <- Z
  # for use with the delta_CI
  wilcox_test$statistic=Z^2
  wilcox_test$parameter=1
  wilcox_test$xvar <- xvar
  wilcox_test$grp <- grp
  class(wilcox_test) <- c(class(wilcox_test),"wilcox.test.rm")
  return(wilcox_test)
}

kruskal.test.rm <- function(xvar,grp){
  kruskal_test <- stats::kruskal.test(xvar~grp)
  n <- length(na.omit(xvar))
  kruskal_test$n <- n
  class(kruskal_test) <- c(class(kruskal_test),"kruskal.test.rm")
  return(kruskal_test)
}

# Conversion Functions  ----------------------------------------------------------------

# First you need to write the conversion functions - these should just convert a statistic to an effect size
# # I don't think this is correct actually
# FtoOmegaSq <- function(Fstat,nu1,nu2){
#   n=nu1+nu2+1
#   k=nu1+1
#   omegaSq=(Fstat-1)/(Fstat+((n-k+1)/(k-1)))
#   # omega squared must be in the range 0,1
#   omegaSq[omegaSq<0] <- 0
#   omegaSq[omegaSq>1] <- 1
#   return(omegaSq)
# }

# From Okada2013
anova_toOmegaSq <- function(anova_summary){
  tbl <- anova_summary[[1]]
  num <- tbl[1,2]-tbl[1,1]*tbl[2,3]
  den <- tbl[1,2]+tbl[2,2]+tbl[2,3]
  return(num/den)
}

# From Smithson
chi_toCramer <- function(chisq_test){
  obs <- chisq_test$observed
  V = sqrt(chisq_test$statistic/(sum(obs)*(min(dim(obs))-1)))
  return(V)
}

# From Smithson
t_toCohen <-function(t_test){
  if (!inherits(t_test,"t.test.rm")) stop("t_test must be a function returned from t.test.rm")
  n1 <- t_test$n[1]
  n2 <- t_test$n[2]
  cohen=abs(t_test$statistic*sqrt((n1+n2)/(n1*n2)))
  return(cohen)
  # This is approximate and actually assumes equal sample sizes in each group
#  abs(2*tStat$statistic/sqrt(N))
}

# REF Tomczak, M., & Tomczak, E. (2014). The need to report effect size estimates revisited. An overview of some recommended measures of effect size. Trends in Sport Sciences, 21(1), 19-25.
# This is for the Kruskal Wallis chi squared
chi_toEpsilonSq <- function(kruskal_test){
  n <- kruskal_test$n
  epsilonSq <- kruskal_test$statistic/(n-1)
}

# This is for eta sq from Wilcoxon R
wilcox_effSize <- function(wilcox_test){
  abs(wilcox_test$Z/sqrt(sum(wilcox_test$n)))
}

# Example usage
# Assume t_test_result is the result of a call to t.test
# t_test_result <- t.test(group1, group2)
# cohen_d_value <- cohen_d(t_test_result)
# print(cohen_d_value)


# And a corresponding function to convert the noon-centrality parameter (from delta_ci) to the effect size
# here, lambda is the non-centrality parameter and N is the total sample size
# This functions help us convert the CI on the non-central distribution to a
# CI on the effect size
# From Steiger
lambda_toOmegaSq <- function(lambda,N){
  lambda/(lambda+N)
}

# From Smithson
lambda_toCramer <- function(lambda,N,k,df){
  sqrt((lambda+df)/(N*(k-1)))
}

lambda_toEpsilon <- function(lamba,N,df){
  (lamba+df)/((N^2-1)/(N+1))
}

lambda_toCohen <- function(lambda,N1,N2){
  lambda/sqrt((N1*N2)/(N1+N2))
}

# Based on the fact that eta squared is Z2/N and Z2 has a chi square distribution
lambda_toWilcoxR <- function(lambda,N){
  sqrt(lambda/N)
}

# Effect Size (with CI) functions ----------------------------------------------------
# Then, we can write a function to compute the effect size and CI

calc_omegaSq <- function(anova_test){
  summary_anova <- summary(anova_test)
  # First calculate the effect size
  omega <- anova_toOmegaSq(summary_anova)

  # use boot for bootstrap re-sampling
  output = c("omega squared"=omega,lower=eff_ci[1],upper=eff_ci[2])
  return(output)
}

# calc_omegaSq <- function(anova_test){
#   summary_anova <- summary(anova_test)
#   # First calculate the effect size
#   omega <- anova_toOmegaSq(summary_anova)
#   # then the delta CI  - this is the CI on the non-central F, not the CI on omega
#   delta_ci <- delta_CI(anova_test)
#   # then, we need to also transform the delta_ci to the effect size scale
#   eff_ci <- lambda_toOmegaSq(delta_ci,sum(summary_anova[[1]]$Df)+1)
#   eff_ci[eff_ci<0] <-0
#   # return the effect size, lower bound, upper bound
#   output = c("omega squared"=omega,lower=eff_ci[1],upper=eff_ci[2])
#   return(output)
# }

calc_CramerV <- function(chisq_test){
  cramer <- chi_toCramer(chisq_test)
  delta_ci <- delta_CI(chisq_test)
  eff_ci <- lambda_toCramer(delta_ci,sum(chisq_test$observed),min(dim(chisq_test$observed), na.rm = TRUE),chisq_test$parameter)
  eff_ci[eff_ci<0] <-0
  # return the effect size, lower bound, upper bound
  output = c("Cramer V"=cramer,lower=eff_ci[1],upper=eff_ci[2])
  return(output)
}

calc_cohenD <- function(t_test){
  cohen <-  t_toCohen(t_test)
  delta_ci <- delta_CI(t_test)
  eff_ci <- lambda_toCohen(delta_ci,t_test$n[1],t_test$n[2])
  # return the effect size, lower bound, upper bound
  eff_ci[eff_ci<0] <-0
  output = c("Cohen d"=cohen,lower=eff_ci[1],upper=eff_ci[2])
  return(output)
}

# This is for the Wilcoxon test
calc_WilcoxonR <- function(wilcox_test){
  r <- wilcox_effSize(wilcox_test)
  delta_ci <- delta_CI(wilcox_test)
  eff_ci <- lambda_toWilcoxR(delta_ci,sum(wilcox_test$n))
  eff_ci[eff_ci<0] <-0
  # return the effect size, lower bound, upper bound
  output = c("Wilcoxn R"=r,lower=eff_ci[1],upper=eff_ci[2])
  return(output)
}

# # WRONG This is for the Wilcoxon test
# calc_etaSq <- function(wilcox_test){
#   eta2 <- chi_toetaSq(wilcox_test)
#   delta_ci <- delta_CI(wilcox_test)
#   eff_ci <- lambda_toEtaSq(delta_ci,sum(wilcox_test$n))
#   eff_ci[eff_ci<0] <-0
#   # return the effect size, lower bound, upper bound
#   output = c("eta squared"=eta2,lower=eff_ci[1],upper=eff_ci[2])
#   return(output)
# }

# This is for the Kruskal Wallis and, for now, the Wilcoxon test
calc_epsilonSq <- function(kruskal_test){
  # The wilcoxon is a special case of the krukal wallis for two groups
  if (inherits(kruskal_test,"wilcox.test.rm")){
    kruskal_test <- kruskal.test.rm(kruskal_test$xvar,kruskal_test$grp)
  }
  epsilonSq <- chi_toEpsilonSq(kruskal_test)
  delta_ci <- delta_CI(kruskal_test)
  eff_ci <- lambda_toEpsilon(delta_ci,kruskal_test$n,kruskal_test$parameter)
  # return the effect size, lower bound, upper bound
  eff_ci[eff_ci<0] <-0
  output = c("Epsilon Squared"=epsilonSq,lower=eff_ci[1],upper=eff_ci[2])
  return(output)
}




