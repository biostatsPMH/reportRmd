

model.summary <- function(model,digits=2,CIwidth = 0.95, whichp = FALSE, ...){
  mcoeff <- coeffSum(model,CIwidth,digits)
  # check units - if any lwr==upr issue warning
  if (any(mcoeff$lwr==mcoeff$upr)) message("Zero-width confidence interval detected. Check predictor units.")

  # Basically, we want to take mcoeff and
  terms <- attr(model$terms, "term.labels")
  if (all(mcoeff$Term %in% terms)){
    # No categorical variables, no need to add any reference data
  } else {
    # Categorical variables
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

      for (i in 2:ncol(df)) data.frame(table(df[[1]],df[,i]))

      catInfo <- data.frame(Term = mcoeff$Term[vpos],
                            Variable= catVar,
                            var_level=reg_lvls,
                            ref_level=ref_lvl)
      cat_vars <- dplyr::bind_rows(cat_vars,catInfo)
    }
    # add to the data frame
    mcoeff <- merge(mcoeff,cat_vars,all.x=T)
  }
  mcoeff$Variable[is.na(mcoeff$Variable)] <- mcoeff$Term[is.na(mcoeff$Variable)]


  tpos <- as.numeric(factor(mcoeff$Variable))
  vars <- sapply(tpos,function(x) terms[x])
  # This needs to calculate "global" p-values for categorical variables
  # works for linear models - need to test all the others!
  # it adds the global-p-value to the dataframe
  if (!all(mcoeff$Term==vars)){
    mcoeff$variable <- vars
    gobal_p <- gp(model)
    mcoeff <- merge(mcoeff,gobal_p,all.x = TRUE,sort=FALSE)
  }
  mcoeff <- mcoeff[order(mcoeff$order),]
  # return(mcoeff)
  # From here - need to add the reference levels as rows to the table
  for (var in unique(mcoeff$variable)) {
    ref_row <- data.frame(var_level = mcoeff$ref_level[1], Est_CI = "Reference", Term = paste0(mcoeff$Variable[1], mcoeff$ref_level[1]), variable = mcoeff$variable[1], Variable = mcoeff$Variable[1], ref_level = mcoeff$ref_level[1])
    mcoeff <- dplyr::bind_rows(ref_row, mcoeff)
    subset_var <- subset(mcoeff, variable == var)
    if (length(unique(subset_var$Term)) > 2) {
      if (whichp == "global" | whichp == "both") {
        var_row <- data.frame(var_level = mcoeff$Variable[1], p_value = mcoeff$global_p[1])
      }
      else { # whichp == "level"
        var_row <- data.frame(var_level = mcoeff$Variable[1])
      }
    }
    else { # two or less levels

      non_ref <- setdiff(subset_var$var_level, mcoeff$ref_level)
      # print(non_ref)
      # print(which(mcoeff$var_level == non_ref))

      mcoeff[which(mcoeff$var_level == non_ref), "p_value"] <- NA
      var_row <- data.frame(var_level = mcoeff$Variable[1], p_value = as.numeric(mcoeff$global_p))
    }
    mcoeff <- dplyr::bind_rows(var_row, mcoeff)
  }


  # reference_row <- c(mcoeff$ref_level[1], "Reference", rep(NA, ncol(mcoeff) - 2))
  # mcoeff <- dplyr::bind_rows(var_row, ref_row, mcoeff)
  # mcoeff <- select(mcoeff, var_level, Est_CI, p_value, est, lwr, upr, global_p)
  View(mcoeff)
}

# Extract model components ------------
coeffSum <- function(model,CIwidth=.95,digits=2,vif=FALSE,whichp="level") {
  CIwidth=CIwidth;digits=digits
  UseMethod("coeffSum",model)
}

coeffSum.default <- function(model,CIwidth=.95,digits=2,vif = FALSE,whichp="level") {
  ms <- summary(model)$coefficients
  ci <- as.data.frame(exp(confint(model,level = CIwidth)))
  names(ci) <- c("lb","ub")
  ci$terms <- rownames(ci)
  rownames(ci) <- NULL
  lvls <- getVarLevels(model)

  cs <- data.frame(
    terms=rownames(ms),
    est=ms[,2],
    p_value = ms[,4]
  ) |> dplyr::left_join(lvls)
  rownames(cs) <- NULL
  cs <- dplyr::left_join(cs,ci)

  xvars <- model$model[,-1,drop=FALSE]
  var_types <- attr(model$terms, "dataClasses")

  events_ss <- lapply(names(xvars),function(v){
    if (var_types[v] == "numeric") return(data.frame(Variable=v,n=nrow(xvars), var = v))
    if (var_types[[v]] == "factor") {
      d2 <-data.frame(table(xvars[[v]]))
      names(d2) <- c("lvl","n")
      d2$var=v
      d2$Variable=paste0(v,d2$lvl)
      print("d2 below")
      print(d2)
      d3 <- data.frame(Variable=v,n=nrow(xvars),var=v)
      print("d3 below")
      print(d3)
      bind_rows(d3,d2)
    }})
  events_ss <- bind_rows(events_ss)
  print(events_ss)
  cs$Est_CI <- apply(cs[,c('est','lwr','upr')],MARGIN = 1,function(x) psthr(x,digits))
  attr(cs,'estLabel') <- betaWithCI("HR",CIwidth)

  events_ss$order <- 1:nrow(events_ss)
  cs <- cs[-1, ]
  cs <- merge(cs,events_ss, all= T, by = "Variable")
  print("mereged cs below")
  print(cs)

  if (vif){
    VIF <- try(GVIF(model),silent = TRUE)
    names(VIF)[1] <- "var"
    cs <- full_join(cs,VIF, by = join_by("Variable" == "var"))
  }

  if (whichp!="level"){
    global_p <- gp(model)
    print(global_p)
    cs <- merge(cs,global_p,all = T)
  }
  for (var in unique(na.omit(cs[["var"]]))) {
    if (var %in% cs[["var"]] & length(which(cs[["var"]] == var)) == 3) {
      p <- cs[max(which(cs[["var"]] == var)), "p_value"]
      cs[min(which(cs[["var"]] == var)), "p_value"] <- p
      cs[max(which(cs[["var"]] == var)), "p_value"] <- NA
    }
  }
  cs <- cs[order(cs$order), ]
  return(cs)
}


#   ms <- summary(model)$coefficients
#   ci <- confint(model,level=CIwidth)
#
#   var_types <- attr(model$terms,"dataClasses")
#   m_df <- model$model
#   refs_df <- data.frame(terms = c(), var = c(), lvl = c())
#   ss <-lapply(names(m_df)[-1],function(v){
#     if (var_types[[v]]=="numeric") return(data.frame(Variable=v,n=nrow(m_df)))
#     if (var_types[[v]]=="factor") {
#       tab <- table(m_df[[v]])
#       for (lev in names(tab)) {
#         if (!(paste0(v, lev) %in% attr(model$coefficients, "names"))) {
#           df <- data.frame(terms = paste0(v, lev), var = v, lvl = lev)
#           refs_df <<- dplyr::bind_rows(refs_df, df)
#         }
#       }
#       names(tab) <- paste0(v, names(tab))
#       d <- data.frame(tab)
#       names(d) <- c("Variable","n")
#       return(d)
#     }
#   })
#   ss <- dplyr::bind_rows(ss)
#
#   cs <- data.frame(
#     Term=rownames(ms),
#     est=ms[,1],
#     p_value = ms[,4],
#     lwr=ci[,1],
#     upr=ci[,2]
#   )
#   rownames(cs) <- NULL
#
#   cs$Est_CI <- apply(cs[,c('est','lwr','upr')],MARGIN = 1,function(x) psthr(x,digits))
#   attr(cs,'estLabel') <- betaWithCI("Estimate",CIwidth)
#   cs <- cs[-1,]
#   term_col <- cs[["Term"]]
#   ex <- getVarLevels(model)
#   ex <- dplyr::bind_rows(ex, refs_df)
#   cs <- full_join(cs, ex, by = c("Term" = "terms"))
#
#   cs <- full_join(ss, cs, by = c("Variable" = "Term"))
#   i = 1
#   for (var in cs[["Variable"]]) {
#     if (!(var %in% term_col)) {
#       cs[i, "Est_CI"] <- "Reference"
#     }
#     i = i+1
#   }
#   old_cs <- cs
#   for (v in setdiff(attr(model$terms, "term.labels"),cs[["Variable"]])) {
#     first_t <- min(which(cs[["var"]] == v))
#     var_row <- data.frame(Variable = v, n = sum(subset(cs, var == v)[["n"]]),
#                           events = sum(subset(cs, var == v)[["events"]]),
#                           var = v, lvl = "NA")
#     if (first_t == 1) {
#       cs <- dplyr::bind_rows(var_row, cs)
#     }
#     else {
#       cs <- dplyr::bind_rows(cs[1:(first_t - 1), ], var_row, cs[- (1:(first_t - 1)), ])
#     }
#   }
#   cs <- cs[, names(old_cs)]
#
#   if (vif) {
#     VIF <- try(GVIF(model),silent = TRUE)
#     if (!inherits(VIF,'try-error')) {
#       if (nrow(VIF)>1){
#         cs <- full_join(cs, VIF, by = c("Variable" = "Covariate"))
#         # vifcol <- character(nrow(cs))
#         # ind <- match(VIF$Covariate,cs$`Variable`)
#
#         # for (x in 1:length(ind)) vifcol[ind[x]] <- niceNum(VIF$VIF[x],digits = digits)
#         # table <- cbind(table,VIF=vifcol)
#       }
#     } else warning('VIF could not be computed for the model.')
#   }
#
#   return(cs)
# }
# coeffSum.geeglm <- function(model,CIwidth=.95,digits=2,vif = FALSE,whichp="level") {
#   ms <- summary(model)$coefficients
#   if (grepl("log",model$family$link)){
#     estFun <- exp
#     ci_mult <- stats::qnorm(1 - (1 - CIwidth)/2)
#   } else {
#     estFun <- identity
#     ci_mult <- stats::qt(1 - (1 - CIwidth)/2,model$df.residual)
#   }
#   ci <- cbind(estFun(ms[,1]-ci_mult*ms[,2]),estFun(ms[,1]+ci_mult*ms[,2]))
#
#   refs_df <- data.frame(terms = c(), var = c(), lvl = c())
#   var_types <- attr(model$terms,"dataClasses")
#   m_df <- model$model
#   ss <-lapply(names(m_df)[-1],function(v){
#     if (var_types[[v]]=="numeric") return(data.frame(Variable=v,n=nrow(m_df)))
#     if (var_types[[v]]=="factor") {
#       tab <- table(m_df[[v]])
#       for (lev in names(tab)) {
#         if (!(paste0(v, lev) %in% attr(model$coefficients, "names"))) {
#           df <- data.frame(terms = paste0(v, lev), var = v, lvl = lev)
#           refs_df <<- dplyr::bind_rows(refs_df, df)
#         }
#       }
#       names(tab) <- paste0(v, names(tab))
#       d <- data.frame(tab)
#       names(d) <- c("Variable","n")
#       return(d)
#     }
#   })
#   ss <- dplyr::bind_rows(ss)
#
#   if (model$family[1]=="binomial"){
#     events <- lapply(names(m_df)[-1],function(v){
#       if (var_types[[v]]=="numeric") return(data.frame(Variable=v,events=sum(model$y)))
#       if (var_types[[v]]=="factor") {
#         tab <- table(model$y,m_df[[v]])
#         colnames(tab) <- paste0(v, names(tab[1,]))
#         d <- data.frame(tab) |>
#           dplyr::filter(Var1==1) |>
#           dplyr::select(-Var1)
#         names(d) <- c("Variable","events")
#         return(d)
#       }
#     })
#     events <- dplyr::bind_rows(events)
#     counts <- merge(ss, events, sort = FALSE)
#   }
#   else {
#     counts <- ss
#   }
#   cs <- data.frame(
#     Term=rownames(ms),
#     est=estFun(ms[,1]),
#     p_value = ms[,4],
#     lwr=ci[,1],
#     upr=ci[,2]
#   )
#   rownames(cs) <- NULL
#
#   cs$Est_CI <- apply(cs[,c('est','lwr','upr')],MARGIN = 1,function(x) psthr(x,digits))
#   beta <- ifelse(model$family$family=="gaussian","Estimate",
#                  ifelse(model$family$family=="binomial","OR",
#                         ifelse(model$family$family=="poisson","RR","GEE Estimate")))
#   attr(cs,'estLabel') <- betaWithCI(beta,CIwidth)
#
#   cs <- cs[-1,]
#   term_col <- cs[["Term"]]
#
#   ex <- getVarLevels(model)
#   ex <- dplyr::bind_rows(ex, refs_df)
#   cs <- full_join(cs, ex, by = c("Term" = "terms"))
#
#
#   cs <- full_join(counts, cs, by = c("Variable" = "Term"))
#   #adding reference levels in:
#   i = 1
#   for (var in cs[["Variable"]]) {
#     if (!(var %in% term_col)) {
#       cs[i, "Est_CI"] <- "Reference"
#     }
#     i = i+1
#   }
#   return(cs)
#}

coeffSum.glm <- function(model,CIwidth=.95,digits=2,vif = FALSE,whichp="level") {
  ms <- summary(model)$coefficients
  ci <- exp(confint(model,level = CIwidth))

  status <- model$model[[1]]
  xvars <- model$model[,-1,drop=FALSE]
  var_types <- attr(model$terms, "dataClasses")

  events_ss <- lapply(names(xvars),function(v){
    if (var_types[v] == "numeric") return(data.frame(Variable=v,events=sum(status),n=nrow(xvars), var = v))
    if (var_types[[v]] == "factor") {
      d1 <- data.frame(table(status,xvars[[v]]))  |>
        dplyr::filter(status==1)   |>
        dplyr::select(-status)
      names(d1) <- c("lvl","events")
      d1$var=v
      d2 <-data.frame(table(xvars[[v]]))
      names(d2) <- c("lvl","n")
      d2$var=v
      d <- merge(d1,d2)
      d$Variable=paste0(v,d$lvl)
      d3 <- data.frame(Variable=v,events=sum(status),n=nrow(xvars),var=v)
      bind_rows(d3,d)
    }})
  events_ss <- bind_rows(events_ss)

  cs <- data.frame(
    Variable=rownames(ms),
    est=ms[,2],
    p_value = ms[,4],
    lwr = ci[,1],
    upr = ci[,2]
  )
  rownames(cs) <- NULL
  cs$Est_CI <- apply(cs[,c('est','lwr','upr')],MARGIN = 1,function(x) psthr(x,digits))
  attr(cs,'estLabel') <- betaWithCI("HR",CIwidth)
  cs <- cs[-1,]
  events_ss$order <- 1:nrow(events_ss)

  cs <- full_join(cs,events_ss,by = "Variable")

  if (vif){
    VIF <- try(GVIF(model),silent = TRUE)
    names(VIF)[1] <- "var"
    cs <- full_join(cs,VIF, by = join_by("Variable" == "var"))
  }

  if (whichp!="level"){
    global_p <- gp(model)
    cs <- merge(cs,global_p,all = T)
  }
  for (var in unique(na.omit(cs[["var"]]))) {
    if (var %in% cs[["var"]] & length(which(cs[["var"]] == var)) == 3) {
      p <- cs[max(which(cs[["var"]] == var)), "p_value"]
      cs[min(which(cs[["var"]] == var)), "p_value"] <- p
      cs[max(which(cs[["var"]] == var)), "p_value"] <- NA
    }
  }
  cs <- cs[order(cs$order), ]
  return(cs)
  # ms <- summary(model)$coefficients
  # ci <- confint(model,level=CIwidth)
  # if (grepl("log",model$family$link)) estFun <- exp else estFun <- identity
  # ci <- estFun(ci)
  #
  # # For Clarina:
  # # Calculate the sample size
  # refs_df <- data.frame(terms = c(), var = c(), lvl = c())
  # var_types <- attr(model$terms,"dataClasses")
  # m_df <- model$model
  # ss <-lapply(names(m_df)[-1],function(v){
  #   if (var_types[[v]]=="numeric") return(data.frame(Variable=v,n=nrow(m_df)))
  #   if (var_types[[v]]=="factor") {
  #
  #     tab <- table(m_df[[v]])
  #
  #     for (lev in names(tab)) {
  #       if (!(paste0(v, lev) %in% attr(model$coefficients, "names"))) {
  #         df <- data.frame(terms = paste0(v, lev), var = v, lvl = lev)
  #         refs_df <<- dplyr::bind_rows(refs_df, df)
  #       }
  #     }
  #
  #     names(tab) <- paste0(v, names(tab))
  #     # print(names(table(m_df[[v]])))
  #     d <- data.frame(tab)
  #     names(d) <- c("Variable","n")
  #     return(d)
  #   }
  # })
  # ss <- dplyr::bind_rows(ss)
  #
  # # ss needs to be added to cs below
  #
  # # Calculate the number of events for binomial models
  # if (model$family[1]=="binomial"){
  #   events <- lapply(names(m_df)[-1],function(v){
  #     if (var_types[[v]]=="numeric") return(data.frame(Variable=v,events=sum(model$y)))
  #     if (var_types[[v]]=="factor") {
  #       tab <- table(model$y,m_df[[v]])
  #       colnames(tab) <- paste0(v, names(tab[1,]))
  #       d <- data.frame(tab) |>
  #         dplyr::filter(Var1==1) |>
  #         dplyr::select(-Var1)
  #       names(d) <- c("Variable","events")
  #       return(d)
  #     }
  #   })
  #   events <- dplyr::bind_rows(events)
  #   counts <- merge(ss, events, sort = FALSE)
  # }
  # else {
  #   counts <- ss
  # }
  #
  # # event counts need to be added as well - we also need to do this for cox_ph models, and gee models with family=binomial
  #
  # cs <- data.frame(
  #   Term=rownames(ms),
  #   est=estFun(ms[,1]),
  #   p_value = ms[,4],
  #   lwr=ci[,1],
  #   upr=ci[,2]
  # )
  # rownames(cs) <- NULL
  # cs$Est_CI <- apply(cs[,c('est','lwr','upr')],MARGIN = 1,function(x) psthr(x,digits))
  # beta <- ifelse(model$family$family=="gaussian","Estimate",
  #                ifelse(model$family$family=="binomial","OR",
  #                       ifelse(model$family$family=="poisson","RR","GLM Estimate")))
  # attr(cs,'estLabel') <- betaWithCI(beta,CIwidth)
  # cs <- cs[-1,]
  #
  # term_col <- cs[["Term"]]
  #
  # ex <- getVarLevels(model)
  # ex <- dplyr::bind_rows(ex, refs_df)
  # cs <- full_join(cs, ex, by = c("Term" = "terms"))
  #
  # print(cs)
  # cs <- full_join(counts, cs, by = c("Variable" = "Term"))
  # #adding reference levels in:
  # i = 1
  # for (var in cs[["Variable"]]) {
  #   if (!(var %in% term_col)) {
  #     print(var)
  #     cs[i, "Est_CI"] <- "Reference"
  #   }
  #   i = i+1
  # }
  # print(cs)
  # old_cs <- cs
  # for (v in setdiff(attr(model$terms, "term.labels"),cs[["Variable"]])) {
  #   first_t <- min(which(cs[["var"]] == v))
  #   var_row <- data.frame(Variable = v, n = sum(subset(cs, var == v)[["n"]]),
  #                         events = sum(subset(cs, var == v)[["events"]]),
  #                         var = v, lvl = "NA")
  #   if (first_t == 1) {
  #     cs <- dplyr::bind_rows(var_row, cs)
  #   }
  #   else {
  #     print(cs[1:(8 - 1), ])
  #     cs <- dplyr::bind_rows(cs[1:(first_t - 1), ], var_row, cs[- (1:(first_t - 1)), ])
  #   }
  # }
  # cs <- cs[, names(old_cs)]
  #
  # if (vif) {
  #   VIF <- try(GVIF(model),silent = TRUE)
  #   if (!inherits(VIF,'try-error')) {
  #     if (nrow(VIF)>1){
  #       cs <- full_join(cs, VIF, by = c("Variable" = "Covariate"))
  #       # vifcol <- character(nrow(cs))
  #       # ind <- match(VIF$Covariate,cs$`Variable`)
  #
  #       # for (x in 1:length(ind)) vifcol[ind[x]] <- niceNum(VIF$VIF[x],digits = digits)
  #       # table <- cbind(table,VIF=vifcol)
  #     }
  #   } else warning('VIF could not be computed for the model.')
  # }
  #
  # return(cs)
}

coeffSum.negbin <- function(model,CIwidth=.95,digits=2,vif = FALSE,whichp="level") {
  ms <- summary(model)$coefficients
  ci <- exp(confint(model,level = CIwidth))

  status <- model$model[[1]]
  xvars <- model$model[,-1,drop=FALSE]
  var_types <- attr(model$terms, "dataClasses")

  events_ss <- lapply(names(xvars),function(v){
    if (var_types[v] == "numeric") return(data.frame(Variable=v,n=nrow(xvars), var = v))
    if (var_types[[v]] == "factor") {
      d2 <-data.frame(table(xvars[[v]]))
      names(d2) <- c("lvl","n")
      d2$var=v
      d2$Variable=paste0(v,d2$lvl)
      d3 <- data.frame(Variable=v,n=nrow(xvars),var=v)
      bind_rows(d3,d2)
    }})
  events_ss <- bind_rows(events_ss)
  cs <- data.frame(
    Variable=rownames(ms),
    est=ms[,2],
    p_value = ms[,4],
    lwr = ci[,1],
    upr = ci[,2]
  )
  rownames(cs) <- NULL
  cs$Est_CI <- apply(cs[,c('est','lwr','upr')],MARGIN = 1,function(x) psthr(x,digits))
  attr(cs,'estLabel') <- betaWithCI("HR",CIwidth)

  events_ss$order <- 1:nrow(events_ss)
  cs <- cs[-1, ]
  cs <- merge(cs,events_ss, all= T, by = "Variable")

  if (vif){
    VIF <- try(GVIF(model),silent = TRUE)
    names(VIF)[1] <- "var"
    cs <- full_join(cs,VIF, by = join_by("Variable" == "var"))
  }

  if (whichp!="level"){
    global_p <- gp(model)
    cs <- merge(cs,global_p,all = T)
  }
  for (var in unique(na.omit(cs[["var"]]))) {
    if (var %in% cs[["var"]] & length(which(cs[["var"]] == var)) == 3) {
      p <- cs[max(which(cs[["var"]] == var)), "p_value"]
      cs[min(which(cs[["var"]] == var)), "p_value"] <- p
      cs[max(which(cs[["var"]] == var)), "p_value"] <- NA
    }
  }
  cs <- cs[order(cs$order), ]
  return(cs)
}

coeffSum.coxph <- function(model,CIwidth=.95,digits=2,vif = FALSE,whichp="level") {
  ms <- summary(model)$coefficients
  ci <- exp(confint(model,level = CIwidth))

  status <- model$model[[1]]
  xvars <- model$model[,-1,drop=FALSE]
  var_types <- attr(model$terms, "dataClasses")

  if (!is.null(model$model)) {
    events_ss <- lapply(names(xvars),function(v){
    if (var_types[v] == "numeric") return(data.frame(Variable=v,events=sum(status),n=nrow(xvars), var = v))
    if (var_types[[v]] == "factor") {
      d1 <- data.frame(table(status,xvars[[v]]))  |>
        dplyr::filter(status==1)   |>
        dplyr::select(-status)
      names(d1) <- c("lvl","events")
      d1$var=v
      d2 <-data.frame(table(xvars[[v]]))
      names(d2) <- c("lvl","n")
      d2$var=v
      d <- merge(d1,d2)
      d$Variable=paste0(v,d$lvl)
      d3 <- data.frame(Variable=v,events=sum(status),n=nrow(xvars),var=v)
      merge(d3,d,all = T)
    }})
  }

  cs <- data.frame(
    Variable=rownames(ms),
    est=ms[,2],
    p_value = ms[,5],
    lwr = ci[,1],
    upr = ci[,2]
  )
  rownames(cs) <- NULL
  cs$Est_CI <- apply(cs[,c('est','lwr','upr')],MARGIN = 1,function(x) psthr(x,digits))
  attr(cs,'estLabel') <- betaWithCI("HR",CIwidth)

  cs <- merge(cs,events_ss,all = T)

  if (vif){
    VIF <- try(GVIF(model),silent = TRUE)
    names(VIF)[1] <- "var"
    cs <- merge(cs,VIF)
  }

  if (whichp!="level"){
    global_p <- gp(model)
    cs <- merge(cs,global_p,all = T)
  }

  return(cs)
}

coeffSum.crr <- function(model,CIwidth=.95,digits=2,vif = FALSE,whichp="level") {
  out <- summary(model, conf.int = CIwidth)
  ms <- out$coef
  ci <- out$conf.int

  status <- model$model[[2]]
  xvars <- model$model[,-(1:2),drop=FALSE]
  var_types <- attr(model$terms, "dataClasses")

  events_ss <- lapply(names(xvars),function(v){
    if (var_types[v] == "numeric") return(data.frame(Variable=v,events=sum(status),n=nrow(xvars)))
    if (var_types[[v]] == "factor") {
      d1 <- data.frame(table(status,xvars[[v]]))  |>
        dplyr::filter(status==1)   |>
        dplyr::select(-status)
      names(d1) <- c("lvl","events")
      d1$var=v
      d2 <-data.frame(table(xvars[[v]]))
      names(d2) <- c("lvl","n")
      d2$var=v
      d <- merge(d1,d2)
      d$Variable=paste0(v,d$lvl)
      d3 <- data.frame(Variable=v,events=sum(status),n=nrow(xvars),var=v)
      merge(d3,d,all = T)
    }})

  cs <- data.frame(
    Variable=rownames(ms),
    est=ms[,1],
    p_value = ms[,5],
    lwr = ci[,3],
    upr = ci[,4]
  )
  rownames(cs) <- NULL
  cs$Est_CI <- apply(cs[,c('est','lwr','upr')],MARGIN = 1,function(x) psthr(x,digits))
  attr(cs,'estLabel') <- betaWithCI("HR",CIwidth)

  cs <- merge(cs,events_ss,all = T)

  if (vif){
    df <- data.frame(y=stats::rnorm(nrow(model$model),0,2))
    df <- cbind(df,xvars)
    mvif <- lm(y~.,data=df)
    VIF <- try(GVIF(mvif),silent = TRUE)
    names(VIF)[1] <- "var"
    cs <- merge(cs,VIF)
  }

  if (whichp!="level"){
    global_p <- gp(model)
    cs <- merge(cs,global_p,all = T)
  }
  return(cs)
}

coeffSum.lme <- function(model,CIwidth=.95,digits=2,vif = FALSE,whichp="level") {
  ms <- summary(model)$tTable
  t_mult <- qt(1 - (1 - CIwidth)/2,ms[,3])

  var_types <- attr(model$terms,"dataClasses")
  m_df <- model$model
  ss <-lapply(names(m_df)[-1],function(v){
    if (var_types[[v]]=="numeric") return(data.frame(Variable=v,n=nrow(m_df)))
    if (var_types[[v]]=="factor") {
      tab <- table(m_df[[v]])
      names(tab) <- paste0(v, names(tab))
      d <- data.frame(tab)
      names(d) <- c("Variable","n")
      return(d)
    }
  })
  ss <- dplyr::bind_rows(ss)

  cs <- data.frame(
    Term=rownames(ms),
    est=ms[,1],
    p_value = ms[,5],
    lwr = ms[,1]-t_mult*ms[,2],
    upr = ms[,1]+t_mult*ms[,2]
  )
  rownames(cs) <- NULL
  cs$Est_CI <- apply(cs[,c('est','lwr','upr')],MARGIN = 1,function(x) psthr(x,digits))
  attr(cs,'estLabel') <- betaWithCI("Estimate",CIwidth)

  cs <- cs[-1,]
  term_col <- cs[["Term"]]
  cs <- full_join(ss, cs, by = c("Variable" = "Term"))
  #adding reference levels in:
  i = 1
  for (var in cs[["Variable"]]) {
    if (!(var %in% term_col)) {
      cs[i, "Est_CI"] <- "Reference"
    }
    i = i+1
  }
  return(cs)
}

coeffSum.polr <- function(model,CIwidth=.95,digits=2,vif = FALSE,whichp="level") {
  ms <- summary(model)$coefficients
  ci <- exp(confint(model,level = CIwidth))

  status <- model$model[[1]]
  xvars <- model$model[,-1,drop=FALSE]
  var_types <- attr(model$terms, "dataClasses")

  events_ss <- lapply(names(xvars),function(v){
    if (var_types[v] == "numeric") return(data.frame(Variable=v,n=nrow(xvars)))
    if (var_types[[v]] == "factor") {
      d2 <-data.frame(table(xvars[[v]]))
      names(d2) <- c("lvl","n")
      d2$var=v
      d2$Variable=paste0(v,d2$lvl)
      d3 <- data.frame(Variable=v,n=nrow(xvars),var=v)
      bind_rows(d3,d2)
    }})
  events_ss <- bind_rows(events_ss)
  cs <- data.frame(
    Variable=rownames(ms),
    est=ms[,2],
    p_value = ms[,4],
    lwr = ci[,1],
    upr = ci[,2]
  )
  rownames(cs) <- NULL
  cs$Est_CI <- apply(cs[,c('est','lwr','upr')],MARGIN = 1,function(x) psthr(x,digits))
  attr(cs,'estLabel') <- betaWithCI("HR",CIwidth)

  events_ss$order <- 1:nrow(events_ss)
  cs <- cs[-1, ]
  cs <- merge(cs,events_ss, all= T, by = "Variable")

  if (vif){
    VIF <- try(GVIF(model),silent = TRUE)
    names(VIF)[1] <- "var"
    cs <- full_join(cs,VIF, by = join_by("Variable" == "var"))
  }

  if (whichp!="level"){
    global_p <- gp(model)
    print(global_p)
    cs <- merge(cs,global_p,all = T)
  }
  cs <- cs[order(cs$order), ]
  return(cs)
}



# Extract data from a fitted model ---------------
get_model_data <- function(model){
  UseMethod("get_model_data", model)
}

get_model_data.default <- function(model){
  return(model$model)
}
# may need to add other methods


# Calculate a global p-value for categorical variables --------
gp <- function(model) {
  UseMethod("gp", model)
}
gp.default <- function(model,CIwidth=.95,digits=2) { # lm, negbin
  terms <- attr(model$terms, "term.labels")
  globalpvalue <- drop1(model,scope=terms,test = "Chisq")
  gp <- data.frame(Variable=rownames(globalpvalue)[-1],
                   global_p = globalpvalue[-1,5])
  attr(gp,"global_p") <-"LRT"
  return(gp)
}


gp.coxph <- function(model,CIwidth=.95,digits=2) {
  terms <- attr(model$terms, "term.labels")
  globalpvalue <- drop1(model,scope=terms,test="Chisq")
  gp <- data.frame(Variable=terms,
                   global_p = globalpvalue[["Pr(>Chi)"]][-1])
  attr(gp,"global_p") <-"LRT"
  return(gp)
}

gp.crr <- function(model,CIwidth=.95,digits=2) {
  terms <- strsplit(trimws(gsub(".*~","", deparse(model$call[[1]]))),"[+]")[[1]]
  terms <- sapply(terms,trimws)
  gp_vals <- data.frame(Variable=terms,
                        global_p = NA)
  rownames(gp_vals) <- NULL
  if (length(terms)>1){
  for (t in terms){
    x <- ifelse(length(setdiff(terms,t))>0,setdiff(terms,t),1)
    eval(parse(text = paste('m2 <-try(crrRx(',paste(paste(setdiff(names(model$model),terms),collapse = "+"),
                                                    "~", x, sep = ""),
                            ',data = model$model))')))

    if (!inherits(m2,"try-error")) {
      degf <- length(grep(t,names(model$coef)))
      gp <- pchisq(2*(model$loglik-m2$loglik),degf)
    } else gp <- NA
    gp_vals$global_p[which(gp_vals$Variable==t)] <- gp
    attr(gp_vals,"global_p") <-"LRT"
  }}
  return(gp_vals)
}

gp.glm <- function(model,CIwidth=.95,digits=2) {
  terms <- attr(model$terms, "term.labels")
  globalpvalue <- drop1(model,scope=terms,test="LRT")
  gp <- data.frame(Variable=rownames(globalpvalue)[-1],
                   global_p = globalpvalue[-1,5])
  attr(gp,"global_p") <-"LRT"
  return(gp)
}
gp.lme <- function(model,CIwidth=.95,digits=2) {
  terms <- attr(model$terms, "term.labels")
  globalpvalue <- drop1(update(model,method="ML"),scope=terms,test = "Chisq")
  gp <- data.frame(Variable=rownames(globalpvalue)[-1],
                   global_p = globalpvalue[-1,5])
  attr(gp,"global_p") <-"LRT"
  return(gp)
}

gp.polr <- function(model,CIwidth=.95,digits=2) {
  terms <- attr(model$terms, "term.labels")
  globalpvalue <- drop1(model,scope=terms,test="ChiSq")
  gp <- data.frame(Variable=rownames(globalpvalue)[-1],
                   global_p = globalpvalue[["Pr(>Chi)"]][-1])
  attr(gp,"global_p") <-"LRT"
  return(gp)
}

gp.gee <- function(model,CIwidth=.95,digits=2) {
  terms <- attr(model$terms, "term.labels")
  terms <- sapply(terms,trimws)
  gp_vals <- data.frame(Variable=terms,
                        global_p = NA)
  rownames(gp_vals) <- NULL
  for (t in terms){
    covariateindex <- grep(paste0("^",t),names(model$coefficients))
    gp <- try(aod::wald.test(b = model$coefficients[covariateindex],
                             Sigma = (model$geese$vbeta)[covariateindex, covariateindex],
                             Terms = seq_len(length(model$coefficients[covariateindex])))$result$chi2[3],
              silent = T)
    if (inherits(gp,"try-error")) gp <- NA
    gp_vals$global_p[which(gp_vals$Variable==t)] <- gp
  }
  attr(gp_vals,"global_p") <-"Wald test"
  return(gp_vals)
}
