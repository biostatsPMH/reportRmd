source("autosum.R")

modelsum <- function(model, digits = 2, CIwidth = 0.95, whichp = FALSE, ...) {
  if (!(whichp %in% c(FALSE, "level", "global", "both"))) {
    stop("argument whichp must be FALSE, or one of 'level', 'global', or 'both'")
  }
  mcoeff <- coeffSum(model,CIwidth,digits)
  print(mcoeff)

  # check units - if any lwr==upr issue warning
  if (any(mcoeff$lwr==mcoeff$upr)) message("Zero-width confidence interval detected. Check predictor units.")

  terms <- attr(model$terms, "term.labels") # this is the list of which predictor (?) variables were used for the model
  if (all(mcoeff$Term %in% terms)) {
    # no categorical variables
  }
  else {
    mcoeff$order <- 1:nrow(mcoeff) # keep track of the order of the coefficients
    cat_vars <- NULL
    for (catVar in setdiff(terms,mcoeff$Term)){
      vpos <- sapply(paste0(catVar,model$xlevels[[catVar]]),function(x) {
        p = which(mcoeff$Term==x)
        if (length(p)==0) p = NA
        return(p)})
      vpos <- vpos[!is.na(vpos)]
      reg_lvls <- gsub(catVar,"",mcoeff$Term[vpos])
      ref_lvl <- setdiff(model$xlevels[[catVar]],reg_lvls)


      print(length(mcoeff$Term[vpos]))
      print(length(catVar))
      print(length(reg_lvls))
      print(length(ref_lvl))
      catInfo <- data.frame(Term = mcoeff$Term[vpos],
                            Variable= catVar,
                            var_level=reg_lvls,
                            ref_level=ref_lvl)
      cat_vars <- dplyr::bind_rows(cat_vars,catInfo)
    }
    mcoeff <- merge(mcoeff,cat_vars,all.x=T)
  }
  mcoeff$Variable[is.na(mcoeff$Variable)] <- mcoeff$Term[is.na(mcoeff$Variable)]

  tpos <- as.numeric(factor(mcoeff$Variable))
  vars <- sapply(tpos,function(x) terms[x])

  if (!all(mcoeff$Term==vars)){
    mcoeff$variable <- vars #why do we need this column?
    gobal_p <- gp(model)
    mcoeff <- merge(mcoeff,gobal_p,all.x = TRUE,sort=FALSE)
    print(mcoeff)
  }
  mcoeff <- mcoeff[order(mcoeff$order),]
  # From here - need to add the reference levels as rows to the table

  for (var in mcoeff$Variable) {
    if (var %in% setdiff(terms,mcoeff$Term)) {
      #find where its first term is
      first_t <- min(which(mcoeff$Variable == var))
      print(first_t)
      ref <- mcoeff[first_t, "ref_level"]
      print("ref below")
      print(ref)
      var_row <- data.frame(Term = var, Variable = var)
      ref_row <- data.frame(Term = paste0(var, ref), Variable = var, Est_CI = "Reference")
      print(nrow(ref_row))
      mcoeff <- bind_rows(mcoeff[1:(first_t - 1), ],
                        var_row, ref_row,
                        mcoeff[- (1:(first_t - 1)), ])
      rownames(mcoeff) <- 1:nrow(mcoeff)
      View(mcoeff)
    }
  }

  # View(mcoeff)
}
