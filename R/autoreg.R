# Function to determine the appropriate regression type based on the response variable
derive_type <- function(data,response,gee) {
  # list2env(attributes(data),envir = environment())
  if (length(response)==2){
    status_var <- data[[response[2]]]
    if (length(unique(na.omit(data[[response[2]]])))==3){
      return("rm_crr")
    }
    if (length(unique(na.omit(data[[response[2]]])))==2){
      return("rm_cox")
    }
  }
  var <- data[[response]]
  if (is.numeric(var)) {
    if (all(var %in% c(0,1))) {
      if (gee) return(c("rm_gee","rm_binomial")) else return(c("rm_glm","rm_binomial"))
    }
    if (length(unique(na.omit(var)))==2){
      message("Binary response detected.\nTo use binomial regression recode data in 0/1 format or change to factor.\nLinear regression will be performed.")
      if (gee) return(c("rm_gee","rm_gaussian")) else return("rm_continuous")
    }
    if (is.integer(var)){
      message("Integer response detected. Poisson regression will be performed, set the type argument to change this.")
      if (gee) return(c("rm_gee","rm_poisson")) else return(c("rm_glm","rm_poisson"))
    }
  }
  if (is.factor(var)) {
    n_levels <- length(levels(var))
    if (n_levels == 2) {
      if (gee) return(c("rm_gee","rm_binomial")) else return(c("rm_glm","rm_binomial"))
    }
    if (inherits(var,"ordered")){
      return("rm_ordinal")
    }
    stop("Unordered factor response detected. Multinomial regression not yet implemented.\nTo use ordinal regression set the response class to an ordered factor.")
  }
  return("unknown")
}



# Fitting Functions ---------------------------

# default is linear
autoreg <- function(response,data,x_var,id=NULL,strata="",family=NULL,offset=NULL, corstr = "independence"){
  UseMethod("autoreg",response)
}

autoreg.rm_default <- function(response,data,x_var,id=NULL,strata="",family=NULL,offset=NULL, corstr = "independence"){
  eval(parse(text = paste('m2 <- lm(',
                          paste(response, "~",x_var, sep = ""),
                          ',data = data)')))
  return(m2)
}

autoreg.rm_lm <- function(response,data,x_var,id=NULL,strata="",family=NULL,offset=NULL, corstr = "independence"){
  eval(parse(text = paste('m2 <- lm(',
                          paste(response, "~",x_var, sep = ""),
                          ',data = data)')))
  return(m2)
}

autoreg.rm_coxph <- function(response,data,x_var,id=NULL,strata="",family=NULL,offset=NULL, corstr = "independence"){
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
  return(m2)
}

autoreg.rm_crr <- function(response,data,x_var,id=NULL,strata="",family=NULL,offset=NULL, corstr = "independence"){
  eval(parse(text = paste('m2 <- crrRx(',paste(paste(response,collapse = "+"),
                                               "~", x_var, sep = ""),
                          ',data = data)')))
  m2$data <- data
  return(m2)
}

autoreg.rm_glm <- function(response,data,x_var,id=NULL,strata="",family=NULL,offset=NULL, corstr = "independence"){
  eval(parse(text = paste("m2 <- glm(",paste(response, "~",x_var, sep = ""),
                          ",family = ",family,",",
                          ifelse(is.null(offset),"",paste("offset=",offset,",")),
                          "data = data)")))
  return(m2)
}

autoreg.rm_gee <-function(response,data,x_var,id=NULL,strata="",family=NULL,offset=NULL, corstr = "independence"){
  eval(parse(text = paste0("m2 <- geepack::geeglm(",paste(response, "~",x_var, sep = ""),
                           ",family = ",family,",",
                           ifelse(is.null(offset),"",paste("offset=",offset,",")),
                           "data = data, id = id, corstr = '",corstr,"')")))
  return(m2)
}

autoreg.rm_negbin <-function(response,data,x_var,id=NULL,strata="",family=NULL,offset=NULL, corstr = "independence"){
  f <- paste(response, "~",x_var,
             ifelse(is.null(offset),"",paste0("+offset(",offset,")")),
             sep = "")
  eval(parse(text = paste("m2 <- MASS::glm.nb(",f,
                          ",link = log",",",
                          "data = data)")))
  return(m2)
}
autoreg.rm_boxcox <-function(response,data,x_var,id=NULL,strata="",family=NULL,offset=NULL, corstr = "independence"){
  eval(parse(text = paste('m2 <- boxcoxfitRx(',
                          paste(response,"~", x_var, sep = ""),
                          ',data = data)')))

  return(m2)
}
autoreg.rm_ordinal <-function(response,data,x_var,id=NULL,strata="",family=NULL,offset=NULL, corstr = "independence"){
  eval(parse(text = paste('m2 = MASS::polr(data = data,',
                          paste(response,"~", x_var, sep = ""),
                          ',method = "logistic",Hess = TRUE)')))
  return(m2)
}


