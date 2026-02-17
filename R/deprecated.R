#'Get univariate summary dataframe
#'
#'Returns a dataframe corresponding to a univariate regression table
#'
#'Univariate summaries for a number of covariates, the type of model can be
#'specified. If unspecified the function will guess the appropriate model based
#'on the response variable.
#'
#'Confidence intervals are extracted using confint where possible. Otherwise
#'Student t distribution is used for linear models and the Normal distribution
#'is used for proportions.
#'
#'returnModels can be used to return a list of the univariate models, which will
#'be the same length as covs. The data used to run each model will include all
#'cases with observations on the response and covariate. For gee models the data
#'are re-ordered so that the ids appear sequentially and proper estimates are
#'given.
#' @seealso [uvsum2()] for the current implementation
#' @return Same as `uvsum2()`
#' @noRd
uvsum <- function(..., markup, sanitize, forceWald) {
  lifecycle::deprecate_warn(
    when = "0.1.2",
    what = "uvsum()",
    with = "uvsum2()",
    details = c(
      "Please update your code to use uvsum2():",
      "uvsum2(response, covs, data, ...)",
      "uvsum will be removed in future versions of reportRmd"
    )
  )

  # Handle deprecated parameters with specific warnings
  if (!missing(markup)) {
    lifecycle::deprecate_warn("0.2.0", "uvsum(markup)")
  }
  if (!missing(sanitize)) {
    lifecycle::deprecate_warn("0.2.0", "uvsum(sanitize)")
  }
  if (!missing(forceWald)) {
    lifecycle::deprecate_warn("0.2.0", "uvsum(forceWald)")
  }

  # Remove deprecated parameters and call uvsum2
  dots <- list(...)
  valid_args <- dots[!names(dots) %in% c("markup", "sanitize", "forceWald")]
  do.call(uvsum2, valid_args)
}

#' Get multivariate summary dataframe
#'
#' Returns a dataframe with the model summary and global p-value for multi-level
#' variables.
#'
#' Global p-values are likelihood ratio tests for lm, glm and polr models. For
#' lme models an attempt is made to re-fit the model using ML and if,successful
#' LRT is used to obtain a global p-value. For coxph models the model is re-run
#' without robust variances with and without each variable and a LRT is
#' presented. If unsuccessful a Wald p-value is returned. For GEE and CRR models
#' Wald global p-values are returned.
#'
#' If the variance inflation factor is requested (VIF=TRUE) then a generalised VIF
#' will be calculated in the same manner as the car package.
#'
#' VIF for competing risk models is computed by fitting a linear model with a
#' dependent variable comprised of the sum of the model independent variables
#' and then calculating VIF from this linear model.
#'
#' @param model fitted model object
#' @param data dataframe containing data
#' @param digits number of digits to round to
#' @param showN boolean indicating sample sizes should be shown for each
#'   comparison, can be useful for interactions
#' @param showEvent boolean indicating if number of events should be shown. Only
#'   available for logistic.
#' @param markup boolean indicating if you want latex markup
#' @param sanitize boolean indicating if you want to sanitize all strings to not
#'   break LaTeX
#' @param nicenames boolean indicating if you want to replace . and _ in strings
#'   with a space.
#' @param CIwidth width for confidence intervals, defaults to 0.95
#' @param vif boolean indicating if the variance inflation factor should be
#'   included. See details
#' @keywords dataframe
#' @importFrom stats na.omit formula model.frame anova qnorm vcov setNames getCall
#' @importFrom utils capture.output
#' @references John Fox & Georges Monette (1992) Generalized Collinearity
#'   Diagnostics, Journal of the American Statistical Association, 87:417,
#'   178-183, DOI: 10.1080/01621459.1992.10475190
#' @references  John Fox and Sanford Weisberg (2019). An {R} Companion to
#'   Applied Regression, Third Edition. Thousand Oaks CA: Sage.
mvsum <- function (model, data, digits=getOption("reportRmd.digits",2), showN = TRUE, showEvent = TRUE, markup = TRUE, sanitize = TRUE, nicenames = TRUE,
                   CIwidth = 0.95, vif=TRUE){
  lifecycle::deprecate_soft("0.2.0","covsum(markup)")
  lifecycle::deprecate_soft("0.2.0","covsum(sanitize)")
  lifecycle::deprecate_warn(
    when = "0.1.2",
    what = "mvsum()",
    with = "modelsummary()",
    details = c(
      "Please update your code to use modelsummary():",
      "mvsum will be removed in future versions of reportRmd"
    )
  )

  if (any(is.na(model$coefficients))) stop(paste0('rm_mvsum can not run when any model coeffcients are NA.\nThe following model coefficients could not be estimated:\n',
                                                  paste(names(model$coefficients)[is.na(model$coefficients)],collapse = ", "),
                                                  "\nPlease re-fit a valid model prior to reporting. Do you need to run droplevels?"))
  if (!markup) {
    lbld <- identity
    addspace <- identity
    lpvalue <- identity
  }
  if (!sanitize)
    sanitizestr <- identity
  if (!nicenames)
    nicename <- identity
  if (inherits(model,c("lm", "lme", "multinom",
                       "survreg", "polr"))) {
    call <- Reduce(paste,
                   deparse(stats::formula(model$terms),
                           width.cutoff = 500))
  }  else if (inherits(model,c("crr"))) {
    call <- paste(deparse(model$formula), collapse = "")
  }  else call <- paste(deparse(model$formula), collapse = "")
  call <- unlist(strsplit(call, "~", fixed = TRUE))[2]
  call <- unlist(strsplit(call, ",", fixed = TRUE))[1]
  if (substr(call, nchar(call), nchar(call)) == "\"")
    call <- substr(call, 1, nchar(call) - 1)
  call <- unlist(strsplit(call, "\"", fixed = TRUE))[1]
  call <- unlist(strsplit(call, "+", fixed = TRUE))
  call <- unlist(strsplit(call, "*", fixed = TRUE))
  call <- unlist(strsplit(call, ":", fixed = TRUE))
  call <- unique(call)
  call <- call[which(is.na(sapply(call, function(cov) {
    charmatch("strata(", cov)
  })) == TRUE)]
  call <- gsub("\\s", "", call)
  type <- class(model)[1]
  if (!isTRUE(model$family$link) && !isTRUE(model$family$link %in% c("log", "logit"))){
    showEvent = FALSE
  }
  if (type == "lm") {
    betanames <- attributes(summary(model)$coef)$dimnames[[1]][-1]
    beta <- "Estimate"
    expnt = FALSE
    ss_data <- model$model
  }
  else if (type == "polr") {
    expnt = TRUE
    betanames <- names(model$coefficients)
    beta <- "OR"
    ss_data <- model$model
  }
  else if (type == "lme") {
    expnt = FALSE
    betanames <- names(model$coef$fixed)[-1]
    beta <- "Estimate"
    ss_data <- model$data
  }
  else if (type == "glm") {
    if (model$family$link == "logit") {
      beta <- "OR"
      expnt = TRUE
    } else if (model$family$link == "log") {
      beta <- "RR"
      expnt = TRUE
    } else {
      beta <- "Estimate"
      expnt = FALSE
    }
    if ( model$family$family=="poisson") showEvent <- FALSE
    betanames <- names(model$coef)[-1]
    ss_data <- model$model
  }
  else if (type == "negbin") {
    betanames <- attributes(summary(model)$coef)$dimnames[[1]][-1]
    beta <- "RR"
    expnt = TRUE
    ss_data <- model$model
    showEvent <- FALSE
  }
  else if (type == "geeglm") {
    if (model$family$link == "logit") {
      beta <- "OR"
      expnt = TRUE
    } else if (model$family$link == "log") {
      beta <- "RR"
      expnt = TRUE
    } else {
      beta <- "Estimate"
      expnt = FALSE
    }
    betanames <- attributes(summary(model)$coef)$row.names[-1]
    if ( model$family$family=="poisson") showEvent <- FALSE
    ss_data <- model$model
  }
  else if (type == "coxph" | type == "crr") {
    beta <- "HR"
    expnt = TRUE
    betanames <- attributes(summary(model)$coef)$dimnames[[1]]
    ss_data <- try(stats::model.frame(model$call$formula, eval(parse(text = paste("data=",
                                                                                  deparse(model$call$data))))), silent = TRUE)
    if (inherits(ss_data,'try-error') & type == "crr") ss_data <- try(model$model)
  }
  else {
    stop("type must be either polr, coxph, glm, lm, geeglm, crr, lme, negbin (or NULL)")
  }
  if (inherits(ss_data,"data.frame")) {
    if ('(weights)' %in% names(ss_data))
      names(ss_data)<- gsub('[(]weights[)]',as.character(model$call[['weights']]),names(ss_data))
    if (any(grepl('offset[(]',names(ss_data)))){
      ot <- which(grepl('offset[(]',names(ss_data)))
      vn <- gsub('[)]','',gsub('offset[(]',"",names(ss_data)[ot]))
      ss_data[[vn]] <- ss_data[,ot]
    }
    data <- ss_data
  } else if (type=='crr'){
    if (missing(data)){
      stop("Data can not be derived from model, data argument must be supplied.")
    } else if (model$n!=nrow(data)) {
      if (showN) stop('For crr models, the supplied data frame can contain only non-missing data.\n Either set showN = FALSE or run na.omit() on a data frame containing only model variables.')
    }
  } else if (type=='coxph'){
    if (missing(data)) stop("Data can not be derived from model, data argument must be supplied.")
    data <- na.omit(data[,c(dimnames(model$y)[[2]],betanames)])
  } else {
    stop("Data can not be derived from model, check model object.")
  }
  beta = betaWithCI(beta, CIwidth)
  ucall = unique(call)
  if (length(setdiff(ucall,names(data)))>0) stop('Currently this function is only implemented to work with standard variable names.\n Try converting the data to a standard data.frame with data.frame(data) and re-running the model to use rm_mvsum.')
  indx = try(matchcovariate(betanames, ucall),silent = TRUE)
  if (is.error(indx)) stop('This function not yet implemented for complex function calls. Try re-specifying the model.')
  for (v in ucall) {
    if (inherits(data[[v]], "character"))
      data[[v]] <- factor(data[[v]])
  }
  if (min(indx) == -1)
    stop("Factor name + level name is the same as another factor name. Please change. Will fix this issue in future.")
  y <- betaindx(indx)
  if (type %in% c("lm", "glm", "negbin","geeglm", "lme")) {
    y <- lapply(y, function(x) {
      x + 1
    })
    betanames <- c("intercept", betanames)
  }
  out <- lapply(y, function(covariateindex) {
    betaname <- betanames[covariateindex]
    betaname <- strsplit(betaname, ":", fixed = TRUE)
    oldcovname <- covnm(betaname[[1]], call)
    oldcovname <- getvarname(oldcovname)
    oldcovname <- paste(oldcovname,collapse = ":")
    levelnameslist <- lapply(betaname, function(level) {
      mapply(function(lvl, cn) {
        result <- ifelse(length(grep(paste0(cn, cn),
                                     lvl)) > 0, unlist(sub(paste0(cn, cn), cn, lvl)),
                         unlist(sub(cn, "", lvl)))
        out <- ifelse(result == "", cn, result)
      }, level, oldcovname)
    })
    levelnames <- unlist(lapply(levelnameslist, function(x) paste(x,
                                                                  collapse = ":")))
    covariatename <- oldcovname
    reference = NULL
    title = NULL
    body = NULL
    if (type == "lme") {
      globalpvalue <- NA
      f <- paste0('. ~ . -',oldcovname)
      if ( length(f)==1){
        m_small <- try(stats::update(model,paste0('. ~ . -',oldcovname),data=data,method='ML'),silent=TRUE)
        if (!is.error(m_small)){
          m_new <- stats::update(model,method='ML')
          globalpvalue <- try(as.vector(stats::na.omit(anova(m_small,m_new)[,"p-value"])),silent=TRUE) # LRT
        }
      }
      if (is.na(globalpvalue)| is.error(globalpvalue)) {
        globalpvalue <- try(aod::wald.test(b = model$coef$fixed[covariateindex],
                                           Sigma = vcov(model)[covariateindex, covariateindex],
                                           Terms = seq_along(covariateindex))$result$chi2[3],silent = TRUE)
      }
    } else if (type  =='negbin'){
      m_small <- try(stats::update(model,paste0('. ~ . -',oldcovname),data=data),silent = TRUE)
      globalpvalue <- try(as.vector(stats::na.omit(anova(m_small,model)[,"Pr(Chi)"])),silent = TRUE)
    } else if (type  =='glm'){
      m_small <- try(stats::update(model,paste0('. ~ . -',oldcovname),data=data),silent = TRUE)
      globalpvalue <- try(as.vector(stats::na.omit(anova(m_small,model,test='LRT')[,"Pr(>Chi)"])),silent = TRUE)
    } else if (type == "polr") {
      m_small <- try(stats::update(model,paste0('. ~ . -',oldcovname),data=data),silent=TRUE)
      globalpvalue <- try(as.vector(stats::na.omit(anova(m_small,model)[,"Pr(Chi)"])),silent=TRUE)
    } else if (type == "crr" ) { # Leave as Wald Test
      globalpvalue <- try(aod::wald.test(b = model$coef[covariateindex],
                                         Sigma = model$var[covariateindex, covariateindex],
                                         Terms = seq_along(covariateindex))$result$chi2[3],
                          silent = TRUE)
    } else if (type=='geeglm'){ # Leave as Wald Test
      globalpvalue <- try(aod::wald.test(b = model$coefficients[covariateindex],
                                         Sigma = (model$geese$vbeta)[covariateindex, covariateindex],
                                         Terms = seq_len(length(model$coefficients[covariateindex])))$result$chi2[3],
                          silent = TRUE)

    } else if (type=='coxph') {
      m_data <- data
      names(m_data)[1] <- 'y'
      m_full <- try(stats::update(model,as.formula('y ~ . '),data=m_data),silent=TRUE)
      m_small <- try(stats::update(model,paste0('y ~ . -',oldcovname),data=m_data),silent=TRUE)
      gp_aov <- try(anova(m_small,m_full),silent = TRUE)

      if (inherits(gp_aov,'try-error')) globalpvalue <- gp_aov else globalpvalue <- as.vector(stats::na.omit(gp_aov[,4]))

    } else {
      m_small <- try(stats::update(model,paste0('. ~ . -',oldcovname),data=data),silent=TRUE)
      globalpvalue <- try(as.vector(stats::na.omit(anova(m_small,model)[,"Pr(>F)"])),silent = TRUE)
    }
    if (is.error(globalpvalue)) globalpvalue <- "NA"
    if (length(globalpvalue)==0) globalpvalue <- "NA"
    if (!identical(lpvalue,identity)) globalpvalue <- lpvalue(globalpvalue,digits)
    if (type == "coxph" | type == "crr") {
      hazardratio <- c(apply(matrix(summary(model, conf.int = CIwidth)$conf.int[covariateindex,
                                                                                c(1, 3, 4)], ncol = 3), 1, psthr,digits))
      pvalues <- c(sapply(summary(model)$coef[covariateindex,
                                              5], lpvalue))
    }
    else if (type %in% c('glm','negbin') & expnt) {
      m <- summary(model, conf.int = CIwidth)$coefficients
      Z_mult = qnorm(1 - (1 - CIwidth)/2)
      hazardratio <- apply(cbind(exp(m[covariateindex, 1]),
                                 exp(m[covariateindex, 1] - Z_mult * m[covariateindex, 2]),
                                 exp(m[covariateindex, 1] + Z_mult * m[covariateindex, 2])), 1, psthr,digits)
      pvalues <- c(sapply(m[covariateindex, 4], lpvalue))
    }
    else if (type == "geeglm" & expnt) {
      m <- summary(model, conf.int = CIwidth)$coefficients
      Z_mult = qnorm(1 - (1 - CIwidth)/2)
      hazardratio <- apply(cbind(exp(m[covariateindex, 1]),
                                 exp(m[covariateindex, 1] - Z_mult * m[covariateindex,2]),
                                 exp(m[covariateindex, 1] + Z_mult * m[covariateindex, 2])), 1, psthr,digits)
      pvalues <- c(sapply(m[covariateindex, 4], lpvalue))
    }
    else if (type == "polr") {
      m <- summary(model)$coefficients
      Z_mult = qnorm(1 - (1 - CIwidth)/2)
      hazardratio <- apply(cbind(exp(m[covariateindex,1]),
                                 exp(m[covariateindex, 1] - Z_mult * m[covariateindex, 2]),
                                 exp(m[covariateindex, 1] + Z_mult * m[covariateindex, 2])), 1, psthr,digits)
      pvalues = stats::pnorm(abs(m[covariateindex, "Value"]/m[covariateindex,
                                                              "Std. Error"]), lower.tail = FALSE) * 2
      pvalues <- c(sapply(pvalues, lpvalue))
    }
    else if (type == "lm" | type == "glm" & !expnt) {
      T_mult = abs(stats::qt((1 - CIwidth)/2, model$df.residual))
      m <- summary(model, conf.int = CIwidth)$coefficients
      hazardratio <- apply(cbind(m[covariateindex, "Estimate"],
                                 m[covariateindex, "Estimate"] - T_mult * m[covariateindex, "Std. Error"],
                                 m[covariateindex, "Estimate"] + T_mult * m[covariateindex, "Std. Error"]), 1, psthr,digits)
      pvalues <- sapply(m[covariateindex, 4], lpvalue)
    }
    else if (type == "geeglm" & !expnt) {
      T_mult = abs(stats::qt((1 - CIwidth)/2, model$df.residual))
      m <- summary(model, conf.int = CIwidth)$coefficients
      hazardratio <- apply(cbind(m[covariateindex, "Estimate"],
                                 m[covariateindex, "Estimate"] - T_mult * m[covariateindex, "Std.err"],
                                 m[covariateindex, "Estimate"] + T_mult * m[covariateindex, "Std.err"]), 1, psthr,digits)
      pvalues <- sapply(m[covariateindex, 4], lpvalue)
    }
    else if (type == "lme") {
      T_mult = abs(stats::qt((1 - CIwidth)/2, summary(model)$fixDF$X))[covariateindex]
      m <- summary(model, conf.int = CIwidth)$tTable
      hazardratio <- apply(cbind(m[covariateindex, 1],
                                 m[covariateindex, 1] - T_mult * m[covariateindex, 2],
                                 m[covariateindex, 1] + T_mult * m[covariateindex,2]), 1, psthr,digits)
      pvalues <- c(sapply(m[covariateindex, 5], lpvalue))
    }
    if (length(betaname[[1]]) == 1) {
      if (!inherits(data[[oldcovname]],"factor")) {
        title <- c(covariatename, hazardratio,pvalues, globalpvalue)
      }     else if (length(levelnames) == 1) {
        title <- c(covariatename, "", pvalues,globalpvalue)
        if (!is.null(data))
          reference <- c(addspace(sanitizestr(names(table(data[,
                                                               which(names(data) == oldcovname)]))[1])),
                         "Reference", "", "")
        body <- c(levelnames, hazardratio, "",
                  "")
      }      else {
        if (!is.null(data)) {
          reference <- c(addspace(sanitizestr(names(table(data[,
                                                               which(names(data) == oldcovname)]))[1])),
                         "Reference", "", "")
        }
        title <- c(covariatename, "", "",
                   globalpvalue)
        body <- cbind(levelnames, hazardratio, pvalues,
                      rep("", length(levelnames)))
      }
    }    else {
      if (length(levelnames) != 1) {
        title <- c(covariatename, "", "",
                   globalpvalue)
        body <- cbind(levelnames, hazardratio, pvalues,
                      rep("", length(levelnames)))
      }      else {
        title <- c(covariatename, hazardratio, pvalues,
                   globalpvalue)

      }
    }
    out <- rbind(title, reference, body)
    if (out[1, 2] == "") {
      if (length(grep(":", title[1])) > 0) {
        ss_N = unlist(lapply(levelnameslist,
                             function(level) {
                               N <- mapply(function(cn, lvl) {
                                 if (cn == lvl) {
                                   nrow(data)
                                 } else {
                                   sum(data[[cn]] == sub(cn,"",lvl))
                                 }
                               }, unlist(strsplit(oldcovname,":")), level)
                               return(min(N))
                             }))
      }
      else {
        ss_N = as.vector(table(data[[oldcovname]]))
      }
      ss_N <- c(nrow(data),ss_N) # Add in the total for the variable
    }
    else {
      ss_N = nrow(data)
    }
    out <- cbind(out, ss_N)
    if (showEvent){
      if (out[1, 2] == "") {
        if (length(grep(":", title[1])) > 0) {
          ss_Event = unlist(lapply(levelnameslist,
                                   function(level) {
                                     Event <- mapply(function(cn, lvl) {
                                       if (cn == lvl) {
                                         nrow(ss_data[which(ss_data[,1] %in% c(1, levels(ss_data[,1])[2])),])
                                       } else {
                                         sum(ss_data[which(ss_data[,1] %in% c(1, levels(ss_data[,1])[2])),][[cn]] == sub(cn,"",lvl))
                                       }
                                     }, unlist(strsplit(oldcovname,":")), level)
                                     return(min(Event))
                                   }))
        }
        else {
          ss_Event = as.vector(table(ss_data[which(ss_data[,1] %in% c(1, levels(ss_data[,1])[2])),][[oldcovname]]))
        }
        ss_Event <- c(nrow(ss_data[which(ss_data[,1] %in% c(1, levels(ss_data[,1])[2])),]),ss_Event) # Add in the total for the variable
      }
      else {
        ss_Event = nrow(ss_data[which(ss_data[,1] %in% c(1, levels(ss_data[,1])[2])),])
      }
      out <- cbind(out, ss_Event)
    }
    rownames(out) <- NULL
    colnames(out) <- NULL
    return(list(out, nrow(out)))
  })

  table <- lapply(out, function(x) {
    return(x[[1]])
  })
  varID <- do.call("c",lapply(table,function(x){
    return(stats::setNames(c(TRUE,rep(FALSE,nrow(x)-1)),x[,1]))
  }))
  index <- unlist(lapply(out, function(x) {
    return(x[[2]])
  }))
  table <- do.call("rbind", lapply(table, data.frame,
                                   stringsAsFactors = FALSE))
  if(length(names(table))==5){
    colnames(table) <- c("Covariate", sanitizestr(beta), "p-value",
                         "Global p-value","N")
  } else  colnames(table) <- c("Covariate", sanitizestr(beta), "p-value",
                               "Global p-value","N","Event")
  table[,"Global p-value"] <- ifelse(table[,'p-value']=='',table[,"Global p-value"],'')
  if (all(table[,"Global p-value"]=='')) table <- table[, -which(colnames(table)=="Global p-value")]
  if (!showN) table <- table[, setdiff(colnames(table),"N")]
  if (!showEvent) table <- table[, setdiff(colnames(table),"Event")]
  if (vif) {
    if (type %in% c('geeglm','lme','negbin')){
      message('VIF not yet implemented for negative binomial, mixed effects or GEE models.')
    } else {
      if (type=='crr'){
        xnm <- intersect(names(data),names(model$coef))
        data$y <- rowSums(data[,xnm],na.rm = TRUE)+stats::rnorm(nrow(data),0,2)
        mvif <- lm(formula = paste('y~',paste(xnm,collapse = '+')),data=data)
        VIF <- try(GVIF(mvif),silent = TRUE)
      } else VIF <- try(GVIF(model),silent = TRUE)
      if (!inherits(VIF,'try-error')) {
        if (nrow(VIF)>1){
          vifcol <- character(nrow(table))
          ind <- match(VIF$Covariate,table$Covariate)
          for (x in 1:length(ind)) vifcol[ind[x]] <- niceNum(VIF$VIF[x],digits = digits)
          table <- cbind(table,VIF=vifcol)
        }
      } else warning('VIF could not be computed for the model.')
    }}
  if (nicenames) table[,1] <- nicename(table[,1])
  colnames(table) <- sapply(colnames(table), lbld)
  attr(table,'covs') <- ucall
  attr(table,"varID") <- varID
  #mc <- paste(utils::capture.output(model$call),collapse="")
  dataArg <- stats::getCall(model)$data
  #dn <- sub(pattern=".*data = (\\w+).*",replacement = "\\1",x=mc)
  dn <- matchdata(dataArg)
  if (is.null(dn)){
    warning('Model data not found. No variable labels will be assigned to variables.')
  } else  {
    attr(table,"data") <- dn
    attr(table,"data call") <- dataArg
    attr(table,"model call") <- nicecall(model$call)
  }
  return(table)
}

#' Create a forest plot using ggplot2 (DEPRECATED)
#'
#'#' @description
#' **Deprecated**: Please use [forestplotMV()] instead.
#'
#' `r lifecycle::badge("deprecated")`
#'
#' This function will be removed in a future version.
#'
#' This function will accept a log or logistic regression fit from glm or
#' geeglm, and display the OR or RR for each variable on the appropriate log
#' scale.
#'
#' @param model an object output from the glm or geeglm function, must be from a
#'   logistic regression
#' @param conf.level controls the width of the confidence interval
#' @param orderByRisk logical, should the plot be ordered by risk
#' @param colours can specify colours for risks less than, 1 and greater than
#'   1.0. Default is red, black, green
#' @param showEst logical, should the risks be displayed on the plot in text
#' @param rmRef logical, should the reference levels be removed for the plot?
#' @param logScale logical, should OR/RR be shown on log scale, defaults to
#'   TRUE, or reportRmd.logScale if set. See https://doi.org/10.1093/aje/kwr156
#'   for why you may prefer a linear scale.
#' @param nxTicks Number of tick marks supplied to the log_breaks function to
#'   produce
#' @importFrom lifecycle deprecate_warn
#' @import ggplot2
#' @importFrom scales log_breaks
#' @keywords plot
#' @return a plot object
#' @export
forestplot2 = function(model,conf.level=0.95,orderByRisk=TRUE,colours='default',showEst=TRUE,rmRef=FALSE,logScale=getOption("reportRmd.logScale",TRUE),nxTicks=5){

  # Warn
  lifecycle::deprecate_warn(
    when = "0.1.2",
    what = "forestplot2()",
    with = "forestplotMV()"
  )


  # Call the new function
  forestplotMV(
    model = model,
    conf.level = conf.level,
    colours = colours,
    showEst = showEst,
    showRef = !rmRef,  # Invert: rmRef=TRUE means showRef=FALSE
    logScale = logScale,
    nxTicks = nxTicks
  )
}

# Summary functions for plots --------------------------

# Survival Curves --------------------------------------------------------------


#' Create Kaplan-Meier or cumulative incidence plots
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `ggkmcif()` was deprecated in version 0.1.2 and will be removed in a future version.
#' Please use [ggkmcif2()] instead.
#'
#' @inheritParams ggkmcif2
#'
#' @return See [ggkmcif2()] for return value details
#'
#' @keywords internal
#' @export
ggkmcif <- function(response,cov=NULL,data,type=NULL,
                    pval = TRUE,HR=FALSE,HR_pval=FALSE,conf.curves=FALSE,conf.type = "log",table = TRUE,
                    times = NULL,xlab = "Time",ylab=NULL ,
                    main = NULL,stratalabs = NULL,strataname = nicename(cov),
                    stratalabs.table=NULL,strataname.table=strataname,
                    median.text=FALSE,median.lines=FALSE,median.CI=FALSE,
                    set.time.text=NULL,set.time.line=FALSE,set.time=5,set.time.CI=FALSE,
                    censor.marks = TRUE,censor.size = 3,censor.stroke = 1.5,
                    fsize = 10, nsize = 3, lsize = 1, psize = 3.5,
                    median.size=3,median.pos=NULL,median.lsize=1,
                    set.size=3,set.pos=NULL,set.lsize=1,
                    ylim=c(0,1), col=NULL,linetype=NULL, xlim=NULL,
                    legend.pos = NULL,  pval.pos=NULL,plot.event=1,event=c("col","linetype"),flip.CIF =FALSE,
                    cut=NULL,eventlabs=NULL,event.name=NULL,Numbers_at_risk_text="Numbers at risk",
                    HR.digits = 2,HR.pval.digits=3, pval.digits=3,
                    median.digits=3,set.time.digits=3,returns = FALSE,print.n.missing=TRUE){

  lifecycle::deprecate_warn(
    when = "0.1.2",
    what = "ggkmcif()",
    with = "ggkmcif2()",
    details = c(
      "Please update your code to use ggkmcif2()",
      "ggkmcif() has been deprecated."
    ))

  # Get all supplied arguments
  args_list <- as.list(match.call())[-1]

  # Call ggkmcif2 with supplied arguments
  do.call(ggkmcif2, args_list)
}

modify_ggkmcif <- function(list_gg){
  lifecycle::deprecate_stop("0.1.0","modify_ggkmcif()","ggkmcif2()",
                            details = "ggkmcif has been replaced by ggkmcif2 which uses cowplot to export plotting elements, modify_ggkmcif has been deprecated.")
}

#' combine components of a call to ggkmci
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `ggkmcif()` was deprecated in version 0.1.2 and will be removed in a future version.
#' @param list_gg A list of ggplot objects from `ggkmcif()`. (Deprecated)
#' Please use [ggkmcif2()] instead.
ggkmcif_paste <- function(list_gg){
  lifecycle::deprecate_stop("0.1.0","ggkmcif_paste()","ggkmcif2()",
                            details = "ggkmcif has been replaced by ggkmcif2 which uses cowplot to export plotting elements, ggkmcif_paste has been deprecated.")

}
