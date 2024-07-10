
model_sum <- function(xvar, response, model, data, digits=2,CIwidth = 0.95) {
  if (missing(model)) {
    print("model missing")
  }
  mcoeff <- coeffSum(model, CIwidth, digits)
  tab <- data.frame(Covariate = xvar)
  i = 2
  for (level in levels(data[[xvar]])) { #this needs to be changed
    tab[i, ] <- level
    i = i + 1
  }
  if ("rm_lm" %in% class(response)) {
    tab[[paste0("Estimate (", CIwidth * 100, "% CI)")]] <- c("", "Reference", mcoeff$Est_CI)
  }
  tab[1, "p-value"] <- mcoeff$p_value
  print(which(is.na(tab[["p-value"]])))
  tab[which(is.na(tab[["p-value"]])), "p-value"] <- ""
  # tab[["N"]] <- NA
  tab[1, "N"] <- length(data[[xvar]]) #do I include NAs???
  print(tab[["N"]])
  i = 2
  for (level in levels(data[[xvar]])) { #this needs to be changed
    print(level)
    tab[i, "N"] <- nrow(subset(data, data[[xvar]] == level))
    i = i + 1
  }
  print(tab)
}
