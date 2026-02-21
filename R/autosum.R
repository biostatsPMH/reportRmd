#'Output a table for multivariate or univariate regression models
#'
#'A dataframe corresponding to a univariate or multivariate regression
#'table. If for_plot = TRUE, estimates and confidence interval bounds will
#'also be displayed separately for easy plotting.
#'
#'Global p-values are likelihood ratio tests for lm, glm and polr models. For
#'lme models an attempt is made to re-fit the model using ML and if,successful
#'LRT is used to obtain a global p-value. For coxph models the model is re-run
#'without robust variances with and without each variable and a LRT is
#'presented. If unsuccessful a Wald p-value is returned. For GEE and CRR models
#'Wald global p-values are returned. For negative binomial models a deviance
#'test is used.
#'
#'If the variance inflation factor is requested (VIF=TRUE) then a generalised VIF
#'will be calculated in the same manner as the car package.
#'
#' As of R 4.4.0 the likelihood profiles are included in base R.
#'
#'The number of decimals places to display the statistics can be changed with
#'digits, but this will not change the display of p-values. If more significant
#'digits are required for p-values then use tableOnly=TRUE and format as
#'desired.
#'@param model model fit
#'@param CIwidth width for confidence intervals, defaults to 0.95
#'@param digits number of digits to round estimates to, does not affect p-values
#'@param vif boolean indicating if the variance inflation factor should be
#'  included. See details
#'@param whichp string indicating whether you want to display p-values for
#'  levels within categorical data ("levels"), global p values ("global"), or
#'  both ("both"). Irrelevant for continuous predictors. When for_plot = TRUE,
#'  global p values will be displayed in a separate column from p values.
#'  If whichp = "levels", global p values will not be included in the outputted
#'  table.
#'@param for_plot boolean indicating whether or not the function will be used
#'  for plotting. Default is FALSE
#'@keywords internal
#'@examples
#' \dontrun{data("pembrolizumab")
#' uv_lm <- lm(age~sex,data=pembrolizumab)
#' m_summary(uv_lm, digits = 3, for_plot = FALSE)
#'
#' mv_binom <- glm(orr~age+sex+cohort,family = 'binomial',data = pembrolizumab)
#' m_summary(mv_binom, whichp = "both", for_plot = TRUE)}
m_summary <- function(model,CIwidth=.95,digits=2,vif = FALSE,whichp="levels",
                      for_plot = FALSE){

  m_coeff <- coeffSum(model,CIwidth,digits)
  m_coeff$Est_CI <- apply(m_coeff[,c('est','lwr','upr')],MARGIN = 1,function(x) psthr(x,digits))

  if (any(!is.na(m_coeff$lwr) & !is.na(m_coeff$upr) & (m_coeff$lwr == m_coeff$upr))) message("Zero-width confidence interval detected. Check predictor units.")

  # For multiply imputed models, swap to first fit for downstream functions
  if (inherits(model, "mira")) {
    model <- model$analyses[[1]]
    vif <- FALSE
    message("VIF is not available for multiply imputed models.")
    if (whichp != "levels") {
      whichp <- "levels"
      message("Global p-values are not yet supported for multiply imputed models.")
    }
  }

  lvls <- getVarLevels(model)
  lvls$ord  <- 1:nrow(lvls)

  cs <- merge(lvls,m_coeff,all.x = TRUE)

  cs <- cs[order(cs$ord),]
  rownames(cs) <- NULL
  # add variable header rows
  cs$header <- NA
  hdr_vars <- names(table(cs$var))[table(cs$var)>1]
  for (v in hdr_vars){
    hdr_rw <- sort(which(cs$var %in% v))[1]
    cs <- dplyr::add_row(cs,.before = hdr_rw,var=v,header=TRUE)
  }
  for (i in 1:nrow(cs)) {
    if (is.na(cs[i, "terms"]) & !is.na(cs[i, "header"])) {
      cs[i, "terms"] <- cs[i, "var"]
    }
    else if (is.na(cs[i, "terms"])) {
      cs[i, "terms"] <- paste0(cs[i, "var"], cs[i, "lvl"])
    }
  }
  terms <- NULL
  if (vif){
    VIF <- try(GVIF(model),silent = TRUE)
    if (inherits(VIF,"try-error") || all(is.infinite(VIF$VIF)) || all(is.na(VIF$VIF))){
      message("Variance inflation factor can not be calculated")
      vif <- FALSE
    } else {
      names(VIF)[1] <- "terms"
      VIF$VIF <- round(VIF$VIF,digits)
      cs <- suppressMessages(dplyr::full_join(cs,VIF, by = dplyr::join_by(terms)))
    }
  }

  if (whichp!="levels"){
    global_p <- gp(model)
    if (!all(is.na(global_p))){
      colnames(global_p) <- c("terms", "global_p")
      cs <- suppressMessages(dplyr::left_join(cs,global_p, by = dplyr::join_by(terms)))
    }
  }

  model_n <- nrow(get_model_data(model))
  model_events <- {
    ev <- get_event_counts(model)
    if (!is.null(ev)) sum(ev == 1) else NULL
  }
  for (i in 1:nrow(cs)) {
    if (!is.na(cs[i, "header"])) {
      v <- cs[i, "var"]
      n <- sum(cs[which(cs[, "var"] == v), "n"], na.rm = TRUE)
      # For interaction terms without main effects, child rows have no per-level
      # N (all NA), so the sum is 0. Fall back to the model's total N.
      if (n == 0 && !is.null(model_n)) n <- model_n
      if ("Events" %in% colnames(cs)) {
        e <- sum(cs[which(cs[, "var"] == v), "Events"], na.rm = TRUE)
        if (e == 0 && !is.null(model_events)) e <- model_events
        cs[i, c("n", "Events")] <- c(n, e)
      }
      else {
        cs[i, "n"] <- n
      }
    }
  }
  cs$ord <- 1:nrow(cs)
  attr(cs,'estLabel') <- attr(m_coeff,'estLabel')
  estLbl <- attr(cs, "estLabel")
  varIDs <- unique(cs[,"var"])

  if (for_plot) {
    return(cs)
  }
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
  if (whichp == "both" & ("global_p" %in% names(cs))) {
    for (i in 1:nrow(cs)) {
      if (!is.na(cs[i, "header"])) {
        cs[i, "p_value"] <- cs[i, "global_p"]
      }
    }
    cs[["global_p"]] <- NULL
  }
  else if (whichp == "global" ) {
    if ("global_p" %in% names(cs)){
      cs[["p_value"]] <- cs[["global_p"]]
      cs[["global_p"]] <- NULL
    } else {
      cs[["p_value"]] <- NA
    }
  }
  cols_to_keep <- na.omit(c("Variable", "Est_CI", "p_value",
                    ifelse("global_p" %in% names(cs),"global_p",NA), "n"))
  new_colnames <- na.omit(c("Variable", estLbl, "p-value",
                    ifelse("global_p" %in% names(cs),"Global p-value",NA),"N"))
  if ("Events" %in% colnames(cs)) {
    cols_to_keep <- c(cols_to_keep, "Events")
    new_colnames <- c(new_colnames, "Event")
  }
  if (vif) {
    cols_to_keep <- c(cols_to_keep, "VIF")
    new_colnames <- c(new_colnames, "VIF")
  }
  cs <- cs[, cols_to_keep]
  colnames(cs) <- new_colnames
  attr(cs, "estLabel") <- estLbl
  attr(cs, "covs") <- varIDs
  return(cs)
}

process_ci <- function(ci_string, digits = 2) {
  # Check if the string contains a valid confidence interval in the format (LB, UB)
  if (grepl("\\([-0-9.eE.NA, ]+\\)", ci_string)) {
    # Extract the lower bound (LB) and upper bound (UB) within the parentheses
    ci_values <- sub(".*\\(([-0-9.eE]+|NA), ([-0-9.eE]+|NA|Inf)\\).*", "\\1,\\2", ci_string)
    ci_values <- unlist(strsplit(ci_values, ","))
    # Convert to numeric or handle "NA" and "Inf" as needed
    lower_bound <- ifelse(ci_values[1] == "NA", NA, suppressWarnings(format(as.numeric(ci_values[1]), nsmall = digits)))
    upper_bound <- ifelse(ci_values[2] %in% c("NA", "Inf"), ci_values[2], suppressWarnings(format(as.numeric(ci_values[2]), nsmall = digits)))
    # Apply rounding and conditions
    if (!is.na(lower_bound)) {
      if (round(as.numeric(lower_bound), digits) == 0) lower_bound <- 0
      else if (abs(as.numeric(lower_bound)) < 0.000001) lower_bound <- 0
      # else if (lower_bound < -10000000) lower_bound <- -Inf
    }
    if (upper_bound == "NA") {
      upper_bound <- Inf
    }
    else if (as.numeric(upper_bound) > 10000000) {
      upper_bound <- Inf
    }
    else if (abs(as.numeric(upper_bound)) < 0.000001) {
      upper_bound <- 0
    }

    # Reconstruct the CI string with the updated bounds
    new_ci <- paste0("(", lower_bound, ", ", upper_bound, ")")
    updated_string <- sub("\\([-0-9.eE, NAInf ]+\\)", new_ci, ci_string)

    return(updated_string)
  } else {
    return(ci_string)  # Return the original string if it doesn't contain a valid CI
  }
}

# Extract model components ------------
#' Extract model coefficient summary
#'
#' S3 generic to extract formatted coefficient summaries from fitted models.
#' @param model A fitted model object.
#' @param CIwidth Confidence interval width (default 0.95).
#' @param digits Number of digits for rounding (default 2).
#' @return A data frame of formatted model coefficients.
#' @keywords internal
#' @export
coeffSum <- function(model,CIwidth=.95,digits=2) {
  CIwidth=CIwidth;digits=digits
  UseMethod("coeffSum",model)
}

#' @export
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
  cs <- merge(cs,ci,all.x = TRUE)

  attr(cs,'estLabel') <- betaWithCI("Estimate",CIwidth)
  return(cs)
}

#' @export
coeffSum.lmerMod <- function(model,CIwidth=.95,digits=2) {
  stop("Method not implemented for lmer fit from lme4,\nre-fit model using lmeTest package.")
}

#' @export
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
  cs <- merge(cs,ci,all.x = TRUE)

  attr(cs,'estLabel') <- betaWithCI("Estimate",CIwidth)
  return(cs)
}

# this code taken from the lmerTest package to compute the Satterthwaite df
sw_df <- function(model){
  qform <- function(x,A){
    sum(x * (A %*% x))
  }
  if (requireNamespace("nlme", quietly = TRUE)) {
    L <- diag(length(nlme::fixef(model)))
  } else stop("Summarising mixed effects models requires the nlme package be installed")


  var_con <- qform(L, vcov(model))
  grad_var_con <- vapply(model@Jac_list, function(x) qform(L,x), numeric(1L))
  satt_denom <- qform(grad_var_con, model@vcov_varpar)
  df <- drop(2 * var_con^2/satt_denom)
  return(df)
}

#' @export
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
  cs <- merge(cs,ci,all.x = TRUE)

  attr(cs,'estLabel') <- betaWithCI("Estimate",CIwidth)
  return(cs)
}

#' @export
coeffSum.crrRx <- function(model,CIwidth=.95,digits=2) {
  ms <- data.frame(model$coeffTbl)
  ci <- try(exp(confint(model,level = CIwidth)),silent = TRUE)
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
  cs <- merge(cs,ci,all.x = TRUE)

  attr(cs,'estLabel') <- betaWithCI("HR",CIwidth)
  return(cs)
}

#' @export
coeffSum.tidycrr <- function(model,CIwidth=.95,digits=2) {
  ms <- model$tidy
  Z_mult = qnorm(1 - (1 - CIwidth)/2)
  cs <- data.frame(terms=ms$term,
                   est=exp(ms$estimate),
                   lwr=exp(ms$estimate-Z_mult*ms$std.error),
                   upr=exp(ms$estimate+Z_mult*ms$std.error),
                   p_value = ms$p.value)

  attr(cs,'estLabel') <- betaWithCI("HR",CIwidth)
  return(cs)
}

#' @export
coeffSum.coxph <- function(model,CIwidth=.95,digits=2) {
  ms <- summary(model)$coefficients
  ci <- exp(as.data.frame(confint(model,level = CIwidth)))
  names(ci) <- c("lwr","upr")
  ci$terms <- rownames(ci)
  rownames(ci) <- NULL
  cs <- data.frame(
    terms=rownames(ms),
    est=ms[,"exp(coef)"],
    p_value = ms[,"Pr(>|z|)"]
  )
  cs <- merge(cs,ci,all.x = TRUE)
  # test for PH assumption
  zph <- try(survival::cox.zph(model),silent = TRUE)
  if (!inherits(zph,"try-error")){
    if (any(zph$table[,"p"]<.05)){
      zph_var <- setdiff(rownames(zph$table)[which(zph$table[,"p"]<.05)],"GLOBAL")
      global_zph <- which(rownames(zph$table)=="GLOBAL")
      if (length(global_zph)>0 & nrow(zph$table)>2){
        if (zph$table[global_zph,"p"]<.05) warning("Global Cox PH assumption may be violated.")
      }
      if (length(zph_var)>0) warning(paste("Warning Cox PH assumption may be violated for these variables:",
                                           paste(setdiff(zph_var,"GLOBAL"),collapse=",")))
    }
  }
  attr(cs,'estLabel') <- betaWithCI("HR",CIwidth)
  return(cs)
}

#' @export
coeffSum.glm <- function(model,CIwidth=.95,digits=2) {
  ms <- data.frame(summary(model)$coefficients)
  if (model$family$link %in% c("logit","log")){
    ci <- suppressMessages(try(as.data.frame(exp(confint(model,level = CIwidth))),silent=TRUE))
    if (inherits(ci,"try-error")){
      ci <- as.data.frame(exp(confint.default(model,level = CIwidth)))
    }
    ms$Estimate <- exp(ms$Estimate)
  } else {
    ci <- suppressMessages(try(as.data.frame(confint(model,level = CIwidth)),silent = TRUE))
    if (inherits(ci,"try-error")){
      ci <- as.data.frame(confint.default(model,level = CIwidth))
    }}
  names(ci) <- c("lwr","upr")
  ci$terms <- rownames(ci)
  rownames(ci) <- NULL

  cs <- data.frame(
    terms=rownames(ms),
    est=ms[,1],
    p_value = ms[,4]
  )

  cs <- merge(cs,ci,all.x = TRUE)

  if (model$family$link == "logit"){
    attr(cs,'estLabel') <- betaWithCI("OR",CIwidth)
  } else if (model$family$link == "log"){
    attr(cs,'estLabel') <- betaWithCI("RR",CIwidth)
  } else {
    attr(cs,'estLabel') <- betaWithCI("Estimate",CIwidth)
  }
  return(cs)
}

#' @export
coeffSum.polr <- function(model,CIwidth=.95,digits=2) {
  ms <- summary(model)$coefficients
  ci <- try(exp(confint(model,level = CIwidth)),silent = TRUE)
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
    est=exp(ms[,1]),
    p_value = pvalues
  )

  cs <- merge(cs,ci,all.x = TRUE)
  cs <- cs[!grepl("[|]",cs$terms),]

  attr(cs,'estLabel') <- betaWithCI("OR",CIwidth)
  return(cs)
}

#' @export
coeffSum.mira <- function(model, CIwidth = 0.95, digits = 2) {
  if (!requireNamespace("mice", quietly = TRUE))
    stop("The mice package is required for multiply imputed model summaries.")

  fit1 <- model$analyses[[1]]

  # --- Validate event count consistency across imputations ---
  events <- lapply(model$analyses, function(f) get_event_counts(f))
  if (!is.null(events[[1]])) {
    event_sums <- sapply(events, sum, na.rm = TRUE)
    if (length(unique(event_sums)) > 1) {
      stop(
        "Event counts are not consistent across imputations.\n",
        "  Event totals per imputation: ",
        paste(event_sums, collapse = ", "), "\n",
        "This may indicate that the outcome variable was imputed, ",
        "which is not supported.\n",
        "Please ensure the outcome is complete or excluded from the imputation model."
      )
    }
  }

  # --- Pool estimates via Rubin's rules ---
  pooled <- mice::pool(model)

  # Detect if exponentiation is needed
  exponentiate <- FALSE
  if (inherits(fit1, "glm") && fit1$family$link %in% c("logit", "log")) {
    exponentiate <- TRUE
  } else if (inherits(fit1, "coxph")) {
    exponentiate <- TRUE
  }

  ps <- summary(pooled, conf.int = TRUE, conf.level = CIwidth,
                exponentiate = exponentiate)

  # Build output in standard coeffSum format
  ci_cols <- grep("%", names(ps), value = TRUE)
  cs <- data.frame(
    terms   = ps$term,
    est     = ps$estimate,
    lwr     = ps[[ci_cols[1]]],
    upr     = ps[[ci_cols[2]]],
    p_value = ps$p.value,
    stringsAsFactors = FALSE
  )

  # Remove intercept row (consistent with other coeffSum methods
  # where intercept is dropped during merge with getVarLevels)
  cs <- cs[!grepl("intercept", cs$terms, ignore.case = TRUE), ]

  # Set estimate label based on model type
  if (inherits(fit1, "glm")) {
    if (fit1$family$link == "logit") {
      attr(cs, "estLabel") <- betaWithCI("OR", CIwidth)
    } else if (fit1$family$link == "log") {
      attr(cs, "estLabel") <- betaWithCI("RR", CIwidth)
    } else {
      attr(cs, "estLabel") <- betaWithCI("Estimate", CIwidth)
    }
  } else if (inherits(fit1, "coxph")) {
    attr(cs, "estLabel") <- betaWithCI("HR", CIwidth)
  } else {
    attr(cs, "estLabel") <- betaWithCI("Estimate", CIwidth)
  }

  return(cs)
}

# Extract coefficients from a fitted model ---------------
#' Extract model coefficients
#'
#' S3 generic to extract raw coefficients from a fitted model.
#' @param model A fitted model object.
#' @return A named numeric vector of model coefficients.
#' @keywords internal
#' @export
get_model_coef <- function(model){
  UseMethod("get_model_coef", model)
}

#' @export
get_model_coef.default <- function(model){
  return(model$coefficients)
}

#' @export
get_model_coef.tidycrr <- function(model){
  return(model$coefs)
}

# Extract data from a fitted model ---------------
#' Extract data from a fitted model
#'
#' S3 generic to extract the data frame used to fit a model.
#' @param model A fitted model object.
#' @return A data frame, or NULL if data cannot be extracted.
#' @keywords internal
#' @export
get_model_data <- function(model){
  UseMethod("get_model_data", model)
}

#' @export
get_model_data.default <- function(model){
  return(NULL)
}
#' @export
get_model_data.lm <- function(model){
  return(model$model)
}
#' @export
get_model_data.lme <- function(model){
  return(model$data)
}
#' @export
get_model_data.lmerMod <- function(model){
  return(model.frame(model))
}
#' @export
get_model_data.lmerModLmerTest <- get_model_data.lmerMod
#' @export
get_model_data.crrRx <- function(model){
  return(model$model)
}
#' @export
get_model_data.tidycrr <- function(model){
  return(model$data)
}
#' @export
get_model_data.polr <- function(model){
  return(model$model)
}
#' @export
get_model_data.coxph <- function(model){
  if (is.null(model$data)){
    df <- try(stats::model.frame(model$call$formula,
                                 eval(parse(text = paste("data=",
                                                         deparse(model$call$data))))), silent = TRUE)
  } else { df <- model$data }
  if (inherits(df,'try-error')) {
    warning ("Model data could not be extracted")
    return(NULL)
  }
  return(df)
}

# may need to add other methods

# Extract event counts from a fitted model ---------------
#' Extract event counts from a fitted model
#'
#' S3 generic to extract event and sample size counts from a fitted model.
#' @param model A fitted model object.
#' @return A named list with event count information, or NULL.
#' @keywords internal
#' @export
get_event_counts <- function(model){
  UseMethod("get_event_counts",model)
}

#' @export
get_event_counts.default <- function(model){
  return(NULL)
}
#' @export
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
#' @export
get_event_counts.crrRx <- function(model){
  md <- get_model_data(model)
  if (is.null(md)) return(NULL)
  return(md[[2]])
}
#' @export
get_event_counts.tidycrr <- function(model){
  md <- get_model_data(model)
  if (is.null(md)) return(NULL)
  stat_var <- names(model$blueprint$ptypes$outcomes)[2]
  return(as.numeric(md[[stat_var]])-1)
}

#' @export
get_event_counts.glm <- function(model){
  if (model$family$family=="binomial"|model$family$family=="quasibinomial"){
    return(model$y)
  }
}

# S3 methods for multiply imputed (mira) objects ---------------
#' @export
get_model_coef.mira <- function(model) {
  return(model$analyses[[1]]$coefficients)
}

#' @export
get_model_data.mira <- function(model) {
  return(get_model_data(model$analyses[[1]]))
}

#' @export
get_event_counts.mira <- function(model) {
  return(get_event_counts(model$analyses[[1]]))
}

# Calculate a global p-value for categorical variables --------
#' Calculate global p-values for categorical variables
#'
#' S3 generic to compute global (Type II/III) p-values for categorical
#' predictors in a fitted model.
#' @param model A fitted model object.
#' @param ... Additional arguments passed to methods.
#' @return A data frame with columns \code{var} and \code{global_p}.
#' @keywords internal
#' @export
gp <- function(model, ...) {
  UseMethod("gp", model)
}

# Shared helper for drop1-based global p-value extraction.
# extract_pvals: function(globalpvalue, terms) -> data.frame(var, global_p)
gp_from_drop1 <- function(model, terms, test = "Chisq",
                           extract_pvals = NULL) {
  globalpvalue <- try(stats::drop1(model, scope = terms, test = test),
                      silent = TRUE)
  if (inherits(globalpvalue, "try-error")) {
    message("Global p values could not be evaluated.")
    return(NA)
  }
  if (!is.null(extract_pvals)) {
    gp <- extract_pvals(globalpvalue, terms)
  } else {
    gp <- data.frame(var = rownames(globalpvalue)[-1],
                     global_p = globalpvalue[-1, 5])
  }
  attr(gp, "global_p") <- "LRT"
  return(gp)
}

#' @export
gp.default <- function(model, ...) { # lm
  data <- model$model
  model <- stats::update(model, data = data)
  terms <- attr(model$terms, "term.labels")
  return(gp_from_drop1(model, terms, test = "Chisq"))
}

#' @export
gp.negbin <- function(model, ...) {
  data <- model$model
  if (any(grepl("offset\\(",names(data)))){
    call_str_vc <- as.character(model$call)
    offset_term <- grep("offset\\(",call_str_vc,value=TRUE)
    offset_str <- gsub("[)]","",gsub(".*offset\\(","",offset_term))
    names(data)[grep("offset\\(",names(data))] <- offset_str
  }
  model <- stats::update(model, data = data)
  terms <- attr(model$terms, "term.labels")
  return(gp_from_drop1(model, terms, test = "Chisq"))
}

#' @export
gp.coxph <- function(model, ...) {
  terms <- attr(model$terms, "term.labels")
  return(gp_from_drop1(model, terms, test = "Chisq",
    extract_pvals = function(gp_result, terms) {
      data.frame(var = terms,
                 global_p = gp_result[["Pr(>Chi)"]][-1])
    }))
}

# Shared helper for competing-risks LRT global p-values (crrRx / tidycrr).
# fit_reduced_fn(model, rhs, data) should return a fitted model or try-error.
# get_loglik_fn(model) should return the log-likelihood scalar.
gp_crr_lrt <- function(model, terms, data, lhs, fit_reduced_fn, get_loglik_fn) {
  gp_vals <- data.frame(var = terms, global_p = NA)
  rownames(gp_vals) <- NULL
  if (length(terms) > 1) {
    for (t in terms) {
      x <- if (length(setdiff(terms, t)) > 0) setdiff(terms, t) else "1"
      m2 <- fit_reduced_fn(model, paste(x, collapse = "+"), data)
      if (!inherits(m2, "try-error")) {
        degf <- length(grep(t, names(model$coef)))
        gp <- pchisq(2 * (get_loglik_fn(model) - get_loglik_fn(m2)), degf)
      } else gp <- NA
      gp_vals$global_p[which(gp_vals$var == t)] <- gp
    }
  }
  attr(gp_vals, "global_p") <- "LRT"
  return(gp_vals)
}

#' @export
gp.crrRx <- function(model, ...) {
  terms <- attr(model$terms, "term.labels")
  lhs <- paste(names(model$model)[1:2], collapse = "+")
  gp_crr_lrt(model, terms, model$model, lhs,
    fit_reduced_fn = function(mod, rhs, md) {
      try(eval(parse(text = paste0("crrRx(", lhs, "~", rhs, ", data = md)"))))
    },
    get_loglik_fn = function(m) m$loglik)
}

#' @export
gp.tidycrr <- function(model, ...) {
  terms <- names(model$blueprint$ptypes$predictors)
  if (!requireNamespace("tidycmprsk", quietly = TRUE)) {
    message("tidycmprsk package not installed. Global p values could not be evaluated.")
    gp_vals <- data.frame(var = terms, global_p = NA)
    attr(gp_vals, "global_p") <- "LRT"
    return(gp_vals)
  }
  md <- get_model_data(model)
  lhs <- gsub("~.*", "", deparse(model$formula))
  gp_crr_lrt(model, terms, md, lhs,
    fit_reduced_fn = function(mod, rhs, md) {
      try(eval(parse(text = paste0("tidycmprsk::crr(", lhs, "~", rhs, ", data = md)"))))
    },
    get_loglik_fn = function(m) m$cmprsk$loglik)
}

#' @export
gp.glm <- function(model, ...) {
  data <- model$model
  if ("(offset)" %in% names(data)){
    call_str <- paste(deparse(model$call),collapse="")
    call_str_vc <- as.character(model$call)
    num_commas <- length(unlist(gregexpr(",", call_str)))
    offset_str <- call_str_vc[num_commas+2]
    names(data)[which(names(data)=="(offset)")] <- offset_str
  }
  model <- stats::update(model, data = data)
  terms <- attr(model$terms, "term.labels")
  return(gp_from_drop1(model, terms, test = "LRT"))
}

#' @export
gp.lme <- function(model, ...) {
  terms <- attr(model$terms, "term.labels")
  model <- stats::update(model, method = "ML")
  return(gp_from_drop1(model, terms, test = "Chisq"))
}

#' @export
gp.lmerMod <- function(model, ...) {
  terms <- attr(terms(model), "term.labels")
  model <- stats::update(model)
  return(gp_from_drop1(model, terms, test = "Chisq",
    extract_pvals = function(gp_result, terms) {
      data.frame(var = rownames(gp_result)[-1],
                 global_p = gp_result$`Pr(Chi)`[-1])
    }))
}

#' @export
gp.lmerModLmerTest <- function(model, ...) {
  terms <- attr(terms(model), "term.labels")
  model <- stats::update(model)
  return(gp_from_drop1(model, terms, test = "Chisq",
    extract_pvals = function(gp_result, terms) {
      data.frame(var = rownames(gp_result),
                 global_p = gp_result$`Pr(>F)`)
    }))
}

#' @export
gp.polr <- function(model, ...) {
  data <- model$model
  model <- stats::update(model, data = data)
  terms <- attr(model$terms, "term.labels")
  return(gp_from_drop1(model, terms, test = "Chisq",
    extract_pvals = function(gp_result, terms) {
      data.frame(var = rownames(gp_result)[-1],
                 global_p = gp_result[["Pr(>Chi)"]][-1])
    }))
}


#' @export
gp.geeglm <- function(model, ...) {
  terms <- attr(model$terms, "term.labels")
  terms <- sapply(terms,trimws)
  gp_vals <- data.frame(var=terms,
                        global_p = NA)
  rownames(gp_vals) <- NULL
  msg <- FALSE
  for (t in terms){
    covariateindex <- grep(paste0("^",t),names(model$coefficients))
    gp <- try(aod::wald.test(b = model$coefficients[covariateindex],
                             Sigma = (model$geese$vbeta)[covariateindex, covariateindex],
                             Terms = seq_len(length(model$coefficients[covariateindex])))$result$chi2[3],
              silent = TRUE)
    if (inherits(gp,"try-error")) {
      msg <- TRUE
      gp <- NA
    }
    gp_vals$global_p[which(gp_vals$var==t)] <- gp
  }
  if (msg) message("Global p values could not be evaluated.")
  attr(gp_vals,"global_p") <-"Wald test"
  return(gp_vals)
}
