#'Format a regression model nicely for 'Rmarkdown'
#'
#'Multivariable (or univariate) regression models are re-formatted for reporting
#'and a global p-value is added for the evaluation of factor variables.
#'
#'Global p-values are likelihood ratio tests for lm, glm and polr models. For
#'lme models an attempt is made to re-fit the model using ML and if successful
#'LRT is used to obtain a global p-value. For lmer models (lme4), if the
#'lmerTest package is installed, Satterthwaite-based p-values and F-test global
#'p-values are used; otherwise Wald z-based p-values and chi-squared LRT global
#'p-values are returned. For glmer models (lme4), Wald z-based p-values are
#'used with chi-squared LRT global p-values. Estimates are exponentiated for
#'binomial (OR) and poisson/negative binomial (RR) families. For coxph models
#'the model is re-run without robust variances with and without each variable
#'and a LRT is presented. If unsuccessful a Wald p-value is returned. For GEE
#'and CRR models Wald global p-values are returned. For negative binomial
#'models a deviance test is used.
#'
#'If the variance inflation factor is requested (VIF=TRUE, default) then a
#'generalised VIF will be calculated in the same manner as the car package.
#'
#'As of version 0.1.1 if global p-values are requested they will be included in
#'the p-value column.
#'
#' As of R 4.4.0 profile likelihood confidence intervals will
#'be calculated automatically and there is no longer an option to force Wald
#'tests.
#'
#'The number of decimals places to display the statistics can be changed with
#'digits, but this will not change the display of p-values. If more significant
#'digits are required for p-values then use tableOnly=TRUE and format as
#'desired.
#'@param model model fit
#'@param data `r lifecycle::badge("deprecated")` data is no longer required as
#'  it is extracted from the model. This argument will be removed in the future
#'@param digits number of digits to round estimates to, does not affect p-values
#'@param covTitle character with the names of the covariate (predictor) column.
#'  The default is to leave this empty for output or, for table only output to
#'  use the column name 'Covariate'.
#'@param showN boolean indicating sample sizes should be shown for each
#'  comparison, can be useful for interactions
#'@param showEvent boolean indicating if number of events should be shown. Only
#'  available for logistic.
#'@param CIwidth width for confidence intervals, defaults to 0.95
#'@param vif boolean indicating if the variance inflation factor should be
#'  included. See details
#'@param whichp string indicating whether you want to display p-values for
#'  levels within categorical data ("levels"), global p values ("global"), or
#'  both ("both"). Irrelevant for continuous predictors.
#'@param caption table caption
#'@param tableOnly boolean indicating if unformatted table should be returned
#'@param p.adjust p-adjustments to be performed. Uses the [p.adjust] function
#'  from base R
#'@param unformattedp boolean indicating if you would like the p-value to be
#'  returned unformatted (ie not rounded or prefixed with '<'). Should be used
#'  in conjunction with the digits argument.
#'@param nicenames boolean indicating if you want to replace . and _ in strings
#'  with a space
#' @param include_unadjusted Logical. If TRUE, includes univariate estimates
#'   alongside multivariable estimates. Default is FALSE. The univariate models
#'   are fit on the same complete-case sample as the multivariable model, so
#'   that the two sets of estimates are directly comparable and the N column
#'   applies to both. Any p-value adjustment requested via `p.adjust` is applied
#'   to the multivariable p-values only. Covariates that appear only inside an
#'   interaction term are appended to the end of the table with an unadjusted
#'   estimate and no adjusted estimate.
#'@param chunk_label Deprecated, previously used in Word to allow
#'  cross-referencing, this should now be done at the chunk level.
#'@param fontsize PDF/HTML output only, manually set the table fontsize
#'@return A character vector of the table source code, unless tableOnly=TRUE in
#'  which case a data frame is returned
#'@export
#'@references John Fox & Georges Monette (1992) Generalized Collinearity
#'  Diagnostics, Journal of the American Statistical Association, 87:417,
#'  178-183, \doi{10.1080/01621459.1992.10475190}
#'@references  John Fox and Sanford Weisberg (2019). An {R} Companion to Applied
#'  Regression, Third Edition. Thousand Oaks CA: Sage.
#' @examples
#' data("pembrolizumab")
#' glm_fit = glm(change_ctdna_group~sex:age+baseline_ctdna+l_size,
#' data=pembrolizumab,family = 'binomial')
#' rm_mvsum(glm_fit)
#'
#' #linear model with p-value adjustment
#' lm_fit=lm(baseline_ctdna~age+sex+l_size+tmb,data=pembrolizumab)
#' rm_mvsum(lm_fit,p.adjust = "bonferroni")
#' #Coxph
#' require(survival)
#' res.cox <- coxph(Surv(os_time, os_status) ~ sex+age+l_size+tmb, data = pembrolizumab)
#' rm_mvsum(res.cox, vif=TRUE)
#'
#' # lmer (lme4 mixed effects model) - single random intercept
#' \donttest{
#' if (require(lme4)){
#' lmer_fit <- lme4::lmer(age ~ sex + pdl1 + (1|cohort), data = pembrolizumab)
#' rm_mvsum(lmer_fit)
#' }
#' }
#'
#' # lmer with multiple random effects and global p-values
#' \donttest{
#' if (require(lme4) && require(geepack)){
#' data(dietox, package = "geepack")
#' dietox$Cu <- as.factor(dietox$Cu)
#' lmer_fit2 <- lme4::lmer(Weight ~ Cu + Time + (1|Pig) + (1|Litter), data = dietox)
#' rm_mvsum(lmer_fit2, whichp = "both")
#' }
#' }
#'
#' # glmer (binomial mixed effects model) - odds ratios
#' \donttest{
#' if (require(lme4)){
#' data(cbpp, package = "lme4")
#' glmer_fit <- lme4::glmer(cbind(incidence, size - incidence) ~ period + (1|herd),
#'   data = cbpp, family = binomial)
#' rm_mvsum(glmer_fit)
#' }
#' }
#'
#' # glmer.nb (negative binomial mixed effects model) - rate ratios
#' \donttest{
#' if (require(lme4) && require(geepack)){
#' data(dietox, package = "geepack")
#' dietox$Cu <- as.factor(dietox$Cu)
#' nb_fit <- lme4::glmer.nb(Weight ~ Cu + Time + (1|Pig), data = dietox)
#' rm_mvsum(nb_fit, whichp = "both")
#' }
#' }
rm_mvsum <- function(model, data, digits=getOption("reportRmd.digits",2),covTitle='',showN=TRUE,showEvent=TRUE,CIwidth=0.95, vif=TRUE,
                     whichp=c("levels","global","both"),
                     caption=NULL,tableOnly=FALSE,p.adjust='none',unformattedp=FALSE,nicenames = TRUE,include_unadjusted=FALSE,
                     chunk_label, fontsize){
  if (unformattedp) formatp <- function(x) {as.numeric(x)}
  whichp <- match.arg(whichp)

  # Handle multiply imputed (mira) models
  is_mira <- inherits(model, "mira")
  if (is_mira) {
    if (!requireNamespace("mice", quietly = TRUE))
      stop("The mice package is required for multiply imputed model summaries.")
    fit1 <- model$analyses[[1]]
  }

  if (!missing(data)) lifecycle::deprecate_soft("0.1.1","rm_mvsum(data)")
  if (!missing(chunk_label)) lifecycle::deprecate_soft("0.1.1","rm_mvsum(chunk_label)")
  model_coef <- get_model_coef(model)
  if (any(is.na(model_coef))) warning(paste0('The following model coefficients could not be estimated and are excluded from the table:\n',
                                             paste(names(model_coef)[is.na(model_coef)],collapse = ", "),
                                             "\nConsider re-fitting the model or running droplevels."))
  # get the table
  tab <- m_summary(model, CIwidth = CIwidth, digits = digits, vif = vif, whichp = whichp, for_plot = FALSE)
  if (include_unadjusted && is_mira) {
    message("Unadjusted estimates are not supported for multiply imputed models.")
    include_unadjusted <- FALSE
  }
  if (include_unadjusted && inherits(model, "logistf")) {
    message("Unadjusted estimates are not supported for logistf models. ",
            "Univariate models would be fit with ordinary (unpenalized) ",
            "logistic regression, which is not comparable to the penalized ",
            "multivariable estimates and may fail under separation.")
    include_unadjusted <- FALSE
  }
  extra_terms <- character(0)
  if (include_unadjusted) {
    m_sum <- m_summary(model, CIwidth = CIwidth, digits = digits, vif = vif, whichp = whichp, for_plot = TRUE)

    # For the univariate models we need the original data with individual
    # columns, not the model frame which may hold Surv() objects or
    # transformed terms
    uv_data <- if (!missing(data)) data else get_uv_data(model)
    ma <- get_model_args(model, data = uv_data)
    uv <- prepare_uv_terms(model, ma, uv_data)

    tabUV <- rm_uvsum(response = uv$response, covs = uv$covs, data = uv$data,
                      digits = digits, CIwidth = CIwidth, whichp = whichp,
                      showEvent = showEvent,
                      tableOnly = TRUE, nicenames = FALSE, unformattedp = unformattedp)
    tab <- combine_uv_mv(tabUV, m_sum, tab,
                         term_labels = ma$term_labels,
                         uv_terms = uv$covs, uv_labels = uv$labels)
    extra_terms <- attr(tab, "extra_terms")
  }
  if (!showN) {
    rmc <- grep("^N(\\s|$)", names(tab))
    if (length(rmc)>0) tab <- tab[,-rmc ]
  }
  if (!showEvent) {
    rmc <- grep("^Event(\\s|$)", names(tab))
    if (length(rmc)>0) tab <- tab[,-rmc ]
  }
  att_tab <- attributes(tab)
  model_terms <- if (is_mira) fit1$terms else tryCatch(model$terms, error = function(e) terms(model))
  header_terms <- c(attr(model_terms, "term.labels"), extra_terms)
  is_header_row <- tab[["Variable"]] %in% header_terms
  to_indent <- which(!is_header_row)
  bold_cells <- cbind(which(is_header_row), rep(1, sum(is_header_row)))
  attr(tab, "to_indent") <- to_indent
  attr(tab, "bold_cells") <- bold_cells

  # Find all estimate columns (contain "95%CI")
  est_cols <- grep("\\([0-9.]+%CI\\)", names(tab), value = TRUE)

  # Temporarily remove already-formatted unadjusted p-values before
  # format_bold_pvalues, which would corrupt string values like "<0.001"
  unadj_p_col <- NULL
  if (include_unadjusted && "Unadjusted p-value" %in% names(tab)) {
    unadj_p_col <- tab[["Unadjusted p-value"]]
    tab[["Unadjusted p-value"]] <- NULL
  }

  # Format and bold p-values
  pv <- format_bold_pvalues(tab, bold_cells,
                            unformattedp = unformattedp, p.adjust = p.adjust)
  tab <- pv$tab; bold_cells <- pv$bold_cells

  # Restore unadjusted p-values in the correct position (after unadjusted estimate)
  if (!is.null(unadj_p_col)) {
    unadj_est_col <- grep("^Unadjusted.*\\([0-9.]+%CI\\)", names(tab), value = TRUE)[1]
    insert_pos <- if (!is.na(unadj_est_col)) which(names(tab) == unadj_est_col) + 1 else 2
    tab <- data.frame(
      tab[, seq_len(insert_pos - 1), drop = FALSE],
      `Unadjusted p-value` = unadj_p_col,
      tab[, insert_pos:ncol(tab), drop = FALSE],
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  # changing UB to Inf, LB to 0 for all estimate columns
  for (ecol in est_cols) {
    tab[[ecol]] <- sapply(tab[[ecol]], process_ci)
  }

  if (nicenames){
    attr(tab,"termnames") <- tab$Variable
    md <- try(get_model_data(if (is_mira) fit1 else model))
    if (inherits(md, "try-error")) {
      warning("Unable to extract data from model, using variable names")
    } else {
      tab$Variable <- replaceLbl(md, tab$Variable)
    }
  }

  names(tab)[1] <- covTitle
  for (a in setdiff(names(att_tab),names(attributes(tab)))) attr(tab,a) <- att_tab[[a]]
  if (tableOnly){
    if (names(tab)[1]=='') names(tab)[1] <- 'Covariate'
    attr(tab, 'to_indent') <- to_indent
    attr(tab,'bold_cells') <- bold_cells
    attr(tab,'dimchk') <- dim(tab)
    return(tab)
  }
  argL <- list(tab=tab,to_indent=to_indent,bold_cells = bold_cells,
               caption=caption, digits = digits,
               chunk_label=ifelse(missing(chunk_label),'NOLABELTOADD',chunk_label))
  if (!missing(fontsize)) argL[['fontsize']] <- fontsize
  do.call(outTable, argL)
}

#' Extract response and predictor names from a fitted model
#'
#' Uses the model's terms object rather than parsing the call, so that
#' formulas held in variables, `.` formulas and transformed terms are all
#' handled.
#'
#' @param model a fitted model
#' @param data optional data frame used to expand a `.` formula
#' @return a list with elements `response`, `predictors` and `term_labels`
#' @keywords internal
get_model_args <- function(model, data = NULL) {

  f <- tryCatch(stats::formula(model), error = function(e) NULL)
  tt <- tryCatch(stats::terms(model), error = function(e) NULL)
  if (is.null(tt) && !is.null(f)) {
    tt <- tryCatch(
      if (is.null(data)) stats::terms(f) else stats::terms(f, data = data),
      error = function(e) NULL)
  }

  # Fall back to parsing the call for model classes with no terms/formula method
  if (is.null(tt)) {
    av <- as.character(model$call)
    av_f <- av[which(grepl("~", av))]
    if (length(av_f) == 0)
      stop("Unable to extract the model formula; supply the data argument.")
    f <- stats::as.formula(av_f[1])
    tt <- stats::terms(f)
  }
  if (is.null(f)) f <- stats::formula(tt)

  term_labels <- attr(tt, "term.labels")
  # drop random effect terms, e.g. (1 | id)
  term_labels <- trimws(term_labels[!grepl("\\|", term_labels)])

  resp_idx <- attr(tt, "response")
  response_vars <- if (!is.null(resp_idx) && resp_idx > 0) all.vars(f[[2]]) else character(0)

  predictor_vars <- unique(trimws(unlist(strsplit(term_labels, "[:*]"))))

  list(
    response = response_vars,
    predictors = predictor_vars,
    term_labels = term_labels
  )
}


#' Retrieve the data used to fit a model
#'
#' Returns the original data frame supplied to the model call where possible,
#' restricted to the rows actually used in the fit. The model frame is only
#' used as a fallback because it can hold composite columns (e.g. `Surv()`)
#' and transformed terms rather than the original variables.
#'
#' @param model a fitted model
#' @return a data frame, or NULL if the data can not be recovered
#' @keywords internal
get_uv_data <- function(model) {

  mf <- tryCatch(get_model_data(model), error = function(e) NULL)

  cl <- tryCatch(stats::getCall(model), error = function(e) NULL)
  orig <- NULL
  if (!is.null(cl) && !is.null(cl[["data"]])) {
    env <- tryCatch(environment(stats::formula(model)), error = function(e) NULL)
    if (is.null(env) || !is.environment(env)) env <- parent.frame()
    orig <- tryCatch(eval(cl[["data"]], envir = env), error = function(e) NULL)
  }

  if (!is.data.frame(orig)) return(mf)

  # restrict the original data to the rows used in the fit
  if (is.data.frame(mf) && nrow(mf) < nrow(orig)) {
    keep <- rownames(orig) %in% rownames(mf)
    if (sum(keep) == nrow(mf)) orig <- orig[keep, , drop = FALSE]
  }
  orig
}


#' Assemble the data and covariate names needed for the univariate models
#'
#' Terms that are not plain columns of the data (e.g. `log(x)`, `poly(x, 2)`)
#' are copied across from the model frame under a syntactically valid name so
#' that [rm_uvsum] can fit them, and are mapped back to their original label
#' afterwards.
#'
#' @param model a fitted model
#' @param ma the output of [get_model_args]
#' @param uv_data data frame returned by [get_uv_data]
#' @return a list with `data`, `covs`, `labels` and `response`
#' @keywords internal
prepare_uv_terms <- function(model, ma, uv_data) {

  if (!is.data.frame(uv_data))
    stop("Unable to extract the data used to fit the model; supply the data argument.")

  mf <- tryCatch(stats::model.frame(model), error = function(e) NULL)
  if (!is.data.frame(mf)) mf <- tryCatch(get_model_data(model), error = function(e) NULL)

  # Expand a composite response (e.g. Surv(time, status)) if the original
  # columns are not available
  response <- ma$response
  missing_resp <- setdiff(response, names(uv_data))
  if (length(missing_resp) > 0 && is.data.frame(mf) && ncol(mf) > 0) {
    y <- mf[[1]]
    if (inherits(y, "Surv") && nrow(uv_data) == nrow(mf)) {
      ym <- as.matrix(y)
      if (ncol(ym) == length(response)) {
        for (i in seq_along(response)) uv_data[[response[i]]] <- ym[, i]
      }
    }
  }
  missing_resp <- setdiff(response, names(uv_data))
  if (length(missing_resp) > 0)
    stop("Response variable(s) not found in the model data: ",
         paste(missing_resp, collapse = ", "))

  covs <- character(0)
  labels <- character(0)
  dropped <- character(0)

  for (cp in ma$predictors) {
    if (cp %in% names(uv_data)) {
      covs <- c(covs, cp); labels <- c(labels, cp); next
    }
    if (is.data.frame(mf) && cp %in% names(mf) && nrow(mf) == nrow(uv_data) &&
        is.null(dim(mf[[cp]]))) {
      nm <- make.names(cp)
      while (nm %in% c(names(uv_data), covs)) nm <- paste0(nm, ".")
      uv_data[[nm]] <- mf[[cp]]
      covs <- c(covs, nm); labels <- c(labels, cp); next
    }
    dropped <- c(dropped, cp)
  }

  if (length(dropped) > 0)
    warning("Unadjusted estimates could not be computed for: ",
            paste(dropped, collapse = ", "))
  if (length(covs) == 0)
    stop("None of the model predictors could be matched to columns in the data.")

  list(data = uv_data, covs = covs, labels = labels, response = response)
}


#' Attach variable/level keys to a summary table
#'
#' Rows whose label matches a model term are the term itself (a continuous
#' covariate, or the header of a categorical or interaction block); all other
#' rows are levels belonging to the most recent term. Deriving the key this way
#' avoids mis-assigning continuous covariates, which have no header row.
#'
#' @param x a summary table
#' @param varcol name of the column holding the labels
#' @param terms character vector of model term labels
#' @return `x` with `var` and `lvl` columns added
#' @keywords internal
add_var_lvl_keys <- function(x, varcol, terms) {
  lab <- as.character(x[[varcol]])
  # compare ignoring whitespace, e.g. "poly(age, 2)" vs "poly(age,2)", but
  # keep the term label itself so that it can be matched against the data
  m <- match(gsub("\\s", "", lab), gsub("\\s", "", terms))
  is_term <- !is.na(m)
  # the first row is always a term, even if the label can not be matched
  if (length(is_term) > 0) is_term[1] <- TRUE
  x$var <- ifelse(is_term, ifelse(is.na(m), lab, terms[m]), NA_character_)
  x <- tidyr::fill(x, "var", .direction = "down")
  x$lvl <- ifelse(is_term, NA_character_, lab)
  x
}


combine_uv_mv <- function(tabUV, m_sum, tabMV, term_labels = NULL,
                          uv_terms = NULL, uv_labels = NULL) {

  ci_pat <- "\\([0-9.]+%CI\\)"

  # Detect estimate column name dynamically
  est_col_uv <- grep(ci_pat, names(tabUV), value = TRUE)[1]
  est_col_mv <- grep(ci_pat, names(tabMV), value = TRUE)[1]

  if (is.na(est_col_uv)) stop("Cannot find estimate column in tabUV")
  if (is.na(est_col_mv)) stop("Cannot find estimate column in tabMV")

  # Standardize column names for internal processing
  tabUV_work <- tabUV
  names(tabUV_work)[names(tabUV_work) == est_col_uv] <- "Est_CI"

  tabMV_work <- tabMV
  names(tabMV_work)[names(tabMV_work) == est_col_mv] <- "Est_CI"

  if (is.null(term_labels)) term_labels <- unique(stats::na.omit(tabMV_work$Variable))

  # Map the sanitized names used for the univariate fits back to the model's
  # term labels so that the two tables can be joined
  if (!is.null(uv_terms) && !is.null(uv_labels)) {
    m <- match(tabUV_work$Covariate, uv_terms)
    tabUV_work$Covariate[!is.na(m)] <- uv_labels[stats::na.omit(m)]
  }
  uv_term_labels <- if (!is.null(uv_labels)) uv_labels else term_labels

  # 1. Reconstruct (var, lvl) keys for both tables
  tabUV2 <- add_var_lvl_keys(tabUV_work, "Covariate", uv_term_labels)

  out <- tabMV_work |>
    dplyr::mutate(mv_order = dplyr::row_number()) |>
    add_var_lvl_keys("Variable", term_labels) |>
    dplyr::left_join(
      tabUV2 |>
        dplyr::select(dplyr::any_of(c("var", "lvl", "Est_CI", "p-value"))) |>
        dplyr::rename(
          `Unadjusted Est_CI` = "Est_CI",
          `Unadjusted p-value` = "p-value"
        ),
      by = c("var", "lvl"),
      relationship = "many-to-one"
    ) |>
    dplyr::rename(
      `Adjusted Est_CI` = "Est_CI",
      `Adjusted p-value` = "p-value"
    ) |>
    dplyr::arrange(.data$mv_order) |>
    dplyr::select(dplyr::any_of(c(
      "Variable", "var",
      "Unadjusted Est_CI", "Unadjusted p-value",
      "Adjusted Est_CI", "Adjusted p-value",
      "N", "Event", "VIF"
    )))

  # 2. Identify main effects that appear only inside interaction terms
  main_terms <- term_labels[!grepl(":", term_labels)]
  interaction_terms <- term_labels[grepl(":", term_labels)]
  uv_only_vars <- character(0)
  if (length(interaction_terms) > 0) {
    uv_only_vars <- setdiff(
      unique(unlist(strsplit(interaction_terms, ":"))),
      main_terms)
    uv_only_vars <- uv_only_vars[uv_only_vars %in% tabUV2$var]
  }

  # Remove 'var' column before appending UV-only rows
  out <- out |> dplyr::select(-"var")
  out_cols <- names(out)
  extra_terms <- character(0)

  # 3. Append UV-only main effects at the end
  if (length(uv_only_vars) > 0) {
    uv_rows_to_add <- list()

    blank_row <- function(v) {
      r <- stats::setNames(as.list(rep(NA, length(out_cols))), out_cols)
      r$Variable <- v
      r
    }

    for (v in uv_only_vars) {
      uv_rows <- tabUV2[which(tabUV2$var == v), , drop = FALSE]
      if (nrow(uv_rows) == 0) next
      extra_terms <- c(extra_terms, v)

      level_rows <- uv_rows[!is.na(uv_rows$lvl), , drop = FALSE]
      term_row <- uv_rows[is.na(uv_rows$lvl), , drop = FALSE]

      if (nrow(level_rows) > 0) {
        # Categorical: header row then one row per level
        new_header <- blank_row(v)
        if (nrow(term_row) > 0 && "N" %in% names(term_row)) new_header$N <- term_row$N[1]
        uv_rows_to_add[[length(uv_rows_to_add) + 1]] <-
          as.data.frame(new_header, check.names = FALSE, stringsAsFactors = FALSE)

        for (k in seq_len(nrow(level_rows))) {
          lv_row <- level_rows[k, ]
          new_row <- blank_row(lv_row$lvl)
          new_row$`Unadjusted Est_CI` <- lv_row$Est_CI
          new_row$`Unadjusted p-value` <- lv_row$`p-value`
          if ("N" %in% names(lv_row)) new_row$N <- lv_row$N
          if ("Event" %in% out_cols && "Event" %in% names(lv_row)) new_row$Event <- lv_row$Event
          uv_rows_to_add[[length(uv_rows_to_add) + 1]] <-
            as.data.frame(new_row, check.names = FALSE, stringsAsFactors = FALSE)
        }
      } else if (nrow(term_row) > 0) {
        # Continuous: single row
        uv_row <- term_row[1, ]
        new_row <- blank_row(v)
        new_row$`Unadjusted Est_CI` <- uv_row$Est_CI
        new_row$`Unadjusted p-value` <- uv_row$`p-value`
        if ("N" %in% names(uv_row)) new_row$N <- uv_row$N
        if ("Event" %in% out_cols && "Event" %in% names(uv_row)) new_row$Event <- uv_row$Event
        uv_rows_to_add[[length(uv_rows_to_add) + 1]] <-
          as.data.frame(new_row, check.names = FALSE, stringsAsFactors = FALSE)
      }
    }

    if (length(uv_rows_to_add) > 0) {
      out <- dplyr::bind_rows(out, uv_rows_to_add)
    }
  }

  # Rename columns back to original format
  names(out)[names(out) == "Unadjusted Est_CI"] <-
    paste0("Unadjusted ", est_col_mv)
  names(out)[names(out) == "Adjusted Est_CI"] <-
    paste0("Adjusted ", est_col_mv)

  attr(out, "extra_terms") <- extra_terms
  out
}
