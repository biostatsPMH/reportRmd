# # Function to determine the appropriate regression type based on the response variable
derive_type <- function(data,response,gee,type,family) {

  if (!is.null(type)){
    if (!type %in% uvmodels$type ) stop(paste("Type must be one of ",
                                              paste0(sort(unique(uvmodels$type)),collapse = ", ")))
    tp_merge <- merge(data.frame(type=type,family=family,gee=gee),uvmodels,all.x = T)
    if (nrow(tp_merge)>1) stop("Can not detect regression type, try specifying family")
    if (is.na(tp_merge$autoreg_class)){
      stop("Can not detect valid regression type")
    } else (return(tp_merge$autoreg_class))
  }
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
      if (gee) {
        rtn <- "rm_gee"
        attr(rtn,"family") = "binomial"
        return(rtn)
      } else {
        rtn <- "rm_glm"
        attr(rtn,"family") = "binomial"
        return(rtn)
      }
      if (length(unique(na.omit(var)))==2){
        message("Binary response detected.\nTo use binomial regression recode data in 0/1 format or change to factor.\nLinear regression will be performed.")
        if (gee) {
          rtn <- "rm_gee"
          attr(rtn,"family") = "gaussian"
          return(rtn)
        } else return("rm_lm")
      }
      if (is.integer(var)){
        message("Integer response detected. Poisson regression will be performed, set the type argument to change this.")
        if (gee) {
          rtn <- "rm_gee"
          attr(rtn,"family") = "poisson"
          return(rtn)
        } else {
          rtn <- "rm_glm"
          attr(rtn,"family") = "poisson"
          return(rtn)
        }
      }
    }
  }
  if (is.factor(var)) {
    n_levels <- length(levels(var))
    if (n_levels == 2) {
      if (gee) {
        rtn <- "rm_gee"
        attr(rtn,"family") = "binomial"
        return(rtn)
      } else {
        rtn <- "rm_glm"
        attr(rtn,"family") = "binomial"
        return(rtn)
      }
    }
    if (inherits(var,"ordered")){
      return("rm_ordinal")
    }
    stop("Unordered factor response detected. Multinomial regression not yet implemented.\nTo use ordinal regression set the response class to an ordered factor.")
  }
  return("unknown")
}
