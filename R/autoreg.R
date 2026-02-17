
# Fitting Functions ---------------------------
#' Fit a regression model based on response type
#'
#' S3 generic that dispatches to the appropriate model fitting function
#' based on the class of the response variable.
#' @param response The response variable (used for S3 dispatch).
#' @param data A data frame containing the model data.
#' @param x_var Character string of the predictor variable name.
#' @param id Optional subject ID for GEE models.
#' @param strata Optional stratification variable.
#' @param family Optional family for GLM models.
#' @param offset Optional offset term.
#' @param corstr Correlation structure for GEE models (default "independence").
#' @return A fitted model object.
#' @keywords internal
#' @export
autoreg <- function(response,data,x_var,id=NULL,strata="",family=NULL,offset=NULL, corstr = "independence"){
  UseMethod("autoreg",response)
}

# Internal helper: fit a model using eval(parse(...)) with a simple formula
# @keywords internal
# @noRd
autoreg_simple_fit <- function(fit_call, response, x_var, data, extra_args = "") {
  m2 <- NULL
  f <- paste(response, "~", x_var, sep = "")
  eval(parse(text = paste0("m2 <- ", fit_call, "(", f,
                            ifelse(nzchar(extra_args), paste0(",", extra_args), ""),
                            ", data = data)")))
  return(m2)
}

#' @export
autoreg.rm_default <- function(response,data,x_var,id=NULL,strata="",family=NULL,offset=NULL, corstr = "independence"){
  autoreg_simple_fit("lm", response, x_var, data)
}

#' @export
autoreg.rm_lm <- autoreg.rm_default

#' @export
autoreg.rm_coxph <- function(response,data,x_var,id=NULL,strata="",family=NULL,offset=NULL, corstr = "independence"){
  m2 <- NULL
  if (all(data[[response[2]]]==0)) stop('No events observed, can\'t fit a Cox model.')
  if (all(data[[response[2]]]==1)) message(paste('All participants with non-missing',x_var,'experienced the event.'))
  f <- paste(paste("survival::Surv(",
                   response[1], ",", response[2], ")",
                   sep = ""), "~", x_var, ifelse(strata ==
                                                   "", "", "+"), paste(strata,
                                                                       collapse = "+"), sep = "")
  if (is.null(id)) {
    eval(parse(text = paste('m2 <- survival::coxph(formula=as.formula(',f,'), data = data)')))
  } else{
    eval(parse(text = paste('m2 <- survival::coxph(formula=as.formula(',f,'),id =',id,', data = data)')))
  }
  m2$data <- try(stats::model.frame(m2$call$formula,data=data), silent = TRUE)
  if (inherits(m2$data,"try-error")) m2$data <- NULL
  return(m2)
}

#' @export
autoreg.rm_crr <- function(response,data,x_var,id=NULL,strata="",family=NULL,offset=NULL, corstr = "independence"){
  autoreg_simple_fit("crrRx", paste(response, collapse = "+"), x_var, data)
}

#' @export
autoreg.rm_glm <- function(response,data,x_var,id=NULL,strata="",family=NULL,offset=NULL, corstr = "independence"){
  m2 <- NULL
  eval(parse(text = paste("m2 <- glm(",paste(response, "~",x_var, sep = ""),
                          ",family = ",family,",",
                          ifelse(is.null(offset),"",paste("offset=",offset,",")),
                          "data = data)")))
  return(m2)
}

#' @export
autoreg.rm_gee <-function(response,data,x_var,id,strata="",family=NULL,offset=NULL, corstr = "independence"){
  # Check if geepack is available ----
  if (!requireNamespace("geepack", quietly = TRUE)) {
    stop("Package 'geepack' is required for GEE models. ",
         "Please install it with: install.packages('geepack')",
         call. = FALSE)
  }
  m2 <- NULL
  idf <- as.numeric(as.factor(data[[id]]))
  class(response) <- "character"
  if (inherits(data[[response]],"factor")) data[[response]] <- as.numeric(data[[response]])-1
  eval(parse(text = paste0("m2 <- geepack::geeglm(",paste(response, "~",x_var, sep = ""),
                           ",family = ",family,",",
                           ifelse(is.null(offset),"",paste("offset=",offset,",")),
                           "data = data, id = idf, corstr = '",corstr,"')")))
  return(m2)
}


#' @export
autoreg.rm_negbin <-function(response,data,x_var,id=NULL,strata="",family=NULL,offset=NULL, corstr = "independence"){
  # Check if MASS is available ----
  if (!requireNamespace("MASS", quietly = TRUE)) {
    stop("Package 'MASS' is required for negative binomial models. ",
         "Please install it with: install.packages('MASS')",
         call. = FALSE)
  }
  m2 <- NULL
  f <- paste(response, "~",x_var,
             ifelse(is.null(offset),"",paste0("+offset(",offset,")")),
             sep = "")
  eval(parse(text = paste("m2 <- MASS::glm.nb(",f,
                          ",link = log",",",
                          "data = data)")))
  return(m2)
}

#' @export
autoreg.rm_boxcox <-function(response,data,x_var,id=NULL,strata="",family=NULL,offset=NULL, corstr = "independence"){
  autoreg_simple_fit("boxcoxfitRx", response, x_var, data)
}

#' @export
autoreg.rm_ordinal <-function(response,data,x_var,id=NULL,strata="",family=NULL,offset=NULL, corstr = "independence"){
  m2 <- NULL
  eval(parse(text = paste('m2 = MASS::polr(data = data,',
                          paste(response,"~", x_var, sep = ""),
                          ',method = "logistic",Hess = TRUE)')))
  return(m2)
}
