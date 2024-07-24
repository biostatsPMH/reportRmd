# need to add whichp
uv_test <- function(response, covs, data, digits=getOption("reportRmd.digits",2),id = NULL, corstr = NULL, family = NULL,
                    type = NULL, offset=NULL, gee=FALSE,strata = 1, markup = TRUE, sanitize = TRUE, nicenames = TRUE,
                    showN = TRUE, showEvent = TRUE, CIwidth = 0.95, reflevel=NULL,returnModels=FALSE,forceWald, vif = T,whichp="both") {


  # assigned_class <- derive_type(data, response, gee)
  # class(response) <- c(class(response), assigned_class)
  modelList <- c()
  for (i in 1:length(covs)) {
    cov_model <- autoreg(response, data, x_var = covs[i], id, strata, family, offset, corstr)
    modelList[[i]] <- cov_model
  }
  summaryList <- c()
  for (model in modelList) {
    cov_summary <- modelsum(model, digits, CIwidth, vif, whichp)
    summaryList <- c(summaryList, cov_summary)
  }
  output <- bind_rows(summaryList)
  return(output)
}

uv_model <- function(modelList, vif = T, CIwidth = 0.95, digits = 2, whichp = "both") {
  summaryList <- list()
  for (i in 1:length(modelList)) {
    cov_summary <- modelsum(modelList[[i]], digits, CIwidth, vif, whichp)
    summaryList[[i]] <- cov_summary
    i = i+1
  }
  # for (model in modelList) {
  #   cov_summary <- modelsum(model, digits, CIwidth, vif, whichp)
  #   summaryList <- c(summaryList, cov_summary)
  # }
  output <- bind_rows(summaryList)
  return(output)
}
