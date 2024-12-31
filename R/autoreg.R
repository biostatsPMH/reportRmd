
# Fitting Functions ---------------------------
# default is linear
autoreg <- function(response,data,x_var,id=NULL,strata="",family=NULL,offset=NULL, corstr = "independence"){
  UseMethod("autoreg",response)
}

autoreg.rm_default <- function(response,data,x_var,id=NULL,strata="",family=NULL,offset=NULL, corstr = "independence"){
  m2 <- NULL
  eval(parse(text = paste('m2 <- lm(',
                          paste(response, "~",x_var, sep = ""),
                          ',data = data)')))
  return(m2)
}

autoreg.rm_lm <- function(response,data,x_var,id=NULL,strata="",family=NULL,offset=NULL, corstr = "independence"){
  m2 <- NULL
  eval(parse(text = paste('m2 <- lm(',
                          paste(response, "~",x_var, sep = ""),
                          ',data = data)')))
  return(m2)
}

autoreg.rm_coxph <- function(response,data,x_var,id=NULL,strata="",family=NULL,offset=NULL, corstr = "independence"){
  m2 <- NULL
  if (all(data[[response[2]]]==0)) stop('No events observed, can\'t fit a Cox model.')
  if (all(data[[response[2]]]==1)) stop(paste('All participants with non-missing',x_var,'experienced the event. \nConsider continuous regression.'))
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

autoreg.rm_crr <- function(response,data,x_var,id=NULL,strata="",family=NULL,offset=NULL, corstr = "independence"){
  m2 <- NULL
  eval(parse(text = paste('m2 <- crrRx(',paste(paste(response,collapse = "+"),
                                               "~", x_var, sep = ""),
                          ',data = data)')))
  return(m2)
}

autoreg.rm_glm <- function(response,data,x_var,id=NULL,strata="",family=NULL,offset=NULL, corstr = "independence"){
  m2 <- NULL
  eval(parse(text = paste("m2 <- glm(",paste(response, "~",x_var, sep = ""),
                          ",family = ",family,",",
                          ifelse(is.null(offset),"",paste("offset=",offset,",")),
                          "data = data)")))
  return(m2)
}

autoreg.rm_gee <-function(response,data,x_var,id,strata="",family=NULL,offset=NULL, corstr = "independence"){
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


autoreg.rm_negbin <-function(response,data,x_var,id=NULL,strata="",family=NULL,offset=NULL, corstr = "independence"){
  m2 <- NULL
  f <- paste(response, "~",x_var,
             ifelse(is.null(offset),"",paste0("+offset(",offset,")")),
             sep = "")
  eval(parse(text = paste("m2 <- MASS::glm.nb(",f,
                          ",link = log",",",
                          "data = data)")))
  return(m2)
}
autoreg.rm_boxcox <-function(response,data,x_var,id=NULL,strata="",family=NULL,offset=NULL, corstr = "independence"){
  m2 <- NULL
  eval(parse(text = paste('m2 <- boxcoxfitRx(',
                          paste(response,"~", x_var, sep = ""),
                          ',data = data)')))

  return(m2)
}

autoreg.rm_ordinal <-function(response,data,x_var,id=NULL,strata="",family=NULL,offset=NULL, corstr = "independence"){
  m2 <- NULL
  eval(parse(text = paste('m2 = MASS::polr(data = data,',
                          paste(response,"~", x_var, sep = ""),
                          ',method = "logistic",Hess = TRUE)')))
  return(m2)
}
