#' Extract model term names
#'
#' S3 generic to extract predictor term names from a fitted model,
#' excluding the intercept.
#' @param model A fitted model object.
#' @return A character vector of term names.
#' @keywords internal
#' @export
mterms <- function(model) {
  UseMethod("mterms", model)
}
#' @export
mterms.default <- function(model){
  nms <- names(model$coefficients)
  # Exclude intercept and inestimable (NA) coefficients
  nms[!grepl("intercept", nms, ignore.case = TRUE) & !is.na(model$coefficients)]
}

#' @export
mterms.tidycrr <- function(model){
  model$tidy$term
}

#' @export
mterms.lme <- function(model){
  if (requireNamespace("nlme", quietly = TRUE)) {
  names(nlme::fixef(model))[!grepl("intercept",
                              names(nlme::fixef(model)),ignore.case = TRUE)]
  } else stop("Summarising mixed effects models requires the nlme package be installed")
}

#' @export
mterms.lmerModLmerTest <- function(model){
  if (requireNamespace("nlme", quietly = TRUE)) {
  names(nlme::fixef(model))[!grepl("intercept",
                              names(nlme::fixef(model)),ignore.case = TRUE)]
  } else stop("Summarising mixed effects models requires the nlme package be installed")
}

getVarLevels <- function(model){
  ord <- NULL
  nt <- function(str) length(strsplit(str,":")[[1]])-1
  vrs<-try(attr(model$terms,"term.labels"),silent = TRUE)
  if (inherits(vrs,"try-error"))  vrs<-try(attr(terms(model),"term.labels"))
  if (is.null(vrs) &inherits(model,"tidycrr")) vrs <- try(names(model$blueprint$ptypes$predictors))
  if (inherits(vrs,"try-error")) stop("Model terms could not be found.")
  if (any(sapply(vrs,nt)>1)) stop("Summary functions will not work with three-way interactions.")
  terms <- mterms(model)
  lvls<-setdiff(terms,vrs)
  df <- data.frame(terms=terms)
  df$var <- ifelse(df$term %in% vrs,df$term,NA)
  df$lvl <- ifelse(df$term %in% lvls,df$term,NA)
  int_terms <- grep("[:]",vrs,value=TRUE)
  for (v in int_terms){
    vr <- unlist(strsplit(v, ":"))
    p <- paste0("^", vr[1], ".*:", vr[2], ".*$")
    vind <- grep(p, df$terms)
    if (all(is.na(df$var[vind]))) df$var[vind] <- v
  }
  if (any(is.na(df$var))){
    vind <- lapply(vrs,function(v) which(grepl(paste0("^",v),df$terms)|grepl(paste0("[:]",v),df$terms)))
    names(vind) <- vrs
    vind <- vind[unlist(lapply(vind,length))>0]
    if (length(unlist(vind)) == 0) return(df)
    vr <- character(max(unlist(vind)))
    for (name in names(vind)) {
      indices <- vind[[name]]
      for (index in indices) {
        if (vr[index] == "") {
          vr[index] <- name
        } else {
          vr[index] <- paste0(vr[index],":", name)
        }
      }
    }
    vr2 <- ifelse(grepl("[:]",df$terms),vr,gsub("[:].*","",vr))
    df$var <- ifelse(is.na(df$var),vr2,df$var)
    lvl2 <- mapply(function(v,l){
      v=strsplit(v,"[:]")[[1]]
      l=strsplit(l,"[:]")[[1]]
      paste0(mapply(function(v,l) sub(v,"",l),v,l,USE.NAMES = FALSE),collapse = ":")
    }, df$var,df$lvl,USE.NAMES=FALSE)
    df$lvl=sub(":$","",sub("^:","",lvl2))
  }
  df$lvl <- ifelse(df$lvl=="NA:NA",NA,df$lvl)
  df$lvl <- ifelse(df$lvl=="NA",NA,df$lvl)
  df$ord <- 1:nrow(df)

  # add sample size and events
  md <- get_model_data(model)
  if (is.null(md)){
    warning("Model data could not be extracted, simplified summary provided")
  } else {
    events = get_event_counts(model)
    if (!is.null(events)) ed <- md[events==1,] else ed <- NULL
  }
  df$v1 <- sub(":.*","",df$var)
  df$v2 <- sub(".*:","",df$var)
  df$v2 <- ifelse(df$v1==df$v2,NA,df$v2)

  if (is.null(md)){
    vcls <- setNames(rep("factor", length(na.omit(unique(c(df$v1,df$v2))))),
                     na.omit(unique(c(df$v1,df$v2))))
  } else {
    vcls <- sapply(na.omit(unique(c(df$v1,df$v2))), function(x) ifelse(is.numeric(md[[x]]),"numeric","factor"))
  }
  int_terms <- unique(grep("[:]",df$var,value=TRUE))

  # Helper: build a frequency/events table for a single variable
  build_freq_row <- function(varname, vcl, md, ed) {
    if (is.null(md)) return(data.frame(var = varname))
    # Determine which factor variables to tabulate
    fct_vars <- names(vcl)[vcl == "factor"]
    if (length(fct_vars) == 0) {
      # All numeric: just report total N
      ntbl <- data.frame(var = varname, lvl = NA, Freq = nrow(md))
      if (!is.null(ed)) ntbl <- merge(ntbl, data.frame(var = varname, lvl = NA, Events = nrow(ed)))
    } else {
      tbl_vars <- paste(fct_vars, collapse = ",")
      ntbl <- eval(parse(text = paste0("data.frame(with(md, table(", tbl_vars, ")))")))
      ntbl$var <- varname
      if (length(fct_vars) == 1) {
        ntbl$lvl <- ntbl[, 1]
      } else {
        ntbl$lvl <- paste0(ntbl[, 1], ":", ntbl[, 2])
      }
      if (!is.null(ed)) {
        etbl <- eval(parse(text = paste0("data.frame(with(ed, table(", tbl_vars, ")))")))
        etbl$var <- varname
        if (length(fct_vars) == 1) {
          etbl$lvl <- etbl[, 1]
        } else {
          etbl$lvl <- paste0(etbl[, 1], ":", etbl[, 2])
        }
        names(etbl) <- sub("Freq", "Events", names(etbl))
        ntbl <- merge(ntbl, etbl)
      }
    }
    return(ntbl)
  }

  freq_list <- list()
  for (i in int_terms) {
    vcl <- vcls[strsplit(i, ":")[[1]]]
    freq_list[[length(freq_list) + 1]] <- build_freq_row(i, vcl, md, ed)
  }
  for (i in setdiff(na.omit(unique(df$var)), int_terms)) {
    vcl <- vcls[i]
    freq_list[[length(freq_list) + 1]] <- build_freq_row(i, vcl, md, ed)
  }
  freq <- dplyr::bind_rows(freq_list)
  freq$n <- freq$Freq
  freq <- freq[,intersect(names(freq),c("var","lvl","n","Events")),drop=FALSE]
  out <- suppressMessages(dplyr::full_join(df,freq))
  out <- dplyr::arrange(out,ord)
  vord <- data.frame(var=na.omit(unique(out$var)),
                     vord=1:length(na.omit(unique(out$var))))
  out <- suppressMessages(dplyr::full_join(out,vord))
  out$ord <- ifelse(is.na(out$ord),0,out$ord)
  out <- dplyr::arrange(out,vord,ord)
  out$ref <- out$ord==0

  # Remove spurious reference rows for interaction terms without main effects.
  # When A:B is in the model but A (a factor) is not a main effect, R estimates
  # all levels of A in the interaction — there is no dropped reference level.
  # The frequency table join creates a ref row for the first factor level, but
  # it is not meaningful and should be removed.
  # Before removing, transfer N/Events to the coefficient rows by matching
  # factor levels within the coefficient lvl strings.
  main_effects <- vrs[!grepl(":", vrs)]
  for (int_var in int_terms) {
    int_components <- strsplit(int_var, ":")[[1]]
    factor_components <- int_components[int_components %in% names(vcls)[vcls == "factor"]]
    # Check if any factor component lacks a main effect
    has_no_main <- any(!factor_components %in% main_effects)
    if (has_no_main) {
      # Frequency-only rows (from full_join, ord==0) have N/Events per level
      freq_rows <- which(out$var == int_var & out$ref == TRUE)
      coeff_rows <- which(out$var == int_var & out$ref == FALSE)
      if (length(freq_rows) > 0 && length(coeff_rows) > 0) {
        # Match frequency rows to coefficient rows by factor level
        for (fi in freq_rows) {
          freq_lvl <- out$lvl[fi]  # e.g., "Female" or "Female:A"
          # Find coefficient row whose lvl contains this factor level
          for (ci in coeff_rows) {
            coeff_lvl <- out$lvl[ci]  # e.g., "sexFemale:age"
            if (!is.na(coeff_lvl) && !is.na(freq_lvl) && grepl(freq_lvl, coeff_lvl, fixed = TRUE)) {
              out$n[ci] <- out$n[fi]
              if ("Events" %in% names(out)) out$Events[ci] <- out$Events[fi]
              break
            }
          }
        }
        out <- out[-freq_rows, , drop = FALSE]
      }
    }
  }

  extra_lvl <- xtr_lvls(out$ref )
  rtn <- out[setdiff(1:nrow(out),extra_lvl),
             intersect(c("var","lvl","n","ref","Events","terms"),names(out))]
  return(rtn)
}

xtr_lvls <- function(vec){
  nx_vec <- c(vec[-1],FALSE)
  pr_vec <- c(FALSE,vec[-length(vec)])
  cnst_T <- which(vec & nx_vec | vec & pr_vec)
  return(cnst_T[c(0,diff(cnst_T))==1])
}


