# need to add whichp
uv_test <- function(response, covs, data, digits=getOption("reportRmd.digits",2),id = NULL, corstr = NULL, family = NULL,
                    type = NULL, offset=NULL, gee=FALSE,strata = 1, markup = TRUE, sanitize = TRUE, nicenames = TRUE,
                    showN = TRUE, showEvent = TRUE, CIwidth = 0.95, reflevel=NULL,returnModels=FALSE,forceWald, vif = T, whichp="both") {


}

uv_model <- function(modelList, vif = T, CIwidth = 0.95, digits = 2, whichp = "both", returnModels = F, tableOnly = F) {
  summaryList <- list()
  #models must be stored in a list, not vector
  for (i in 1:length(modelList)) {
    cov_summary <- modelsum(modelList[[i]], digits, CIwidth, vif, whichp)
    summaryList[[i]] <- cov_summary
    i = i+1
  }
  output <- bind_rows(summaryList)
  cols_to_keep <- c()
  if ("VIF" %in% colnames(output)) {
    cols_to_keep <- c(cols_to_keep, "VIF")
  }
  if ("events" %in% colnames(output)) {
    cols_to_keep <- c(cols_to_keep, "events")
  }
  output <- output[, c("Variable", "Est_CI", "p_value", "n", cols_to_keep)]

  if (returnModels) return(list(output,models=modelList))

  if (tableOnly) {
    return(output)
  }
  #return nicely formatted table
  return(output)
}
