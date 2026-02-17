#'Format a regression model nicely for 'Rmarkdown'
#'
#'Multivariable (or univariate) regression models are re-formatted for reporting
#'and a global p-value is added for the evaluation of factor variables.
#'
#'Global p-values are likelihood ratio tests for lm, glm and polr models. For
#'lme models an attempt is made to re-fit the model using ML and if,successful
#'LRT is used to obtain a global p-value. For coxph models the model is re-run
#'without robust variances with and without each variable and a LRT is
#'presented. If unsuccessful a Wald p-value is returned. For GEE and CRR models
#'Wald global p-values are returned. For negative binomial models a deviance
#'test is used.
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
#'   alongside multivariable estimates. Default is FALSE.
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
  if (any(is.na(model_coef))) stop(paste0('rm_mvsum cannot run when any model coeffcients are NA.\nThe following model coefficients could not be estimated:\n',
                                          paste(names(model_coef)[is.na(model_coef)],collapse = ", "),
                                          "\nPlease re-fit a valid model prior to reporting. Do you need to run droplevels?"))
  # get the table
  tab <- m_summary(model, CIwidth = CIwidth, digits = digits, vif = vif, whichp = whichp, for_plot = FALSE)
  if (include_unadjusted && is_mira) {
    message("Unadjusted estimates are not supported for multiply imputed models.")
    include_unadjusted <- FALSE
  }
  if (include_unadjusted) {
    m_sum <- m_summary(model, CIwidth = CIwidth, digits = digits, vif = vif, whichp = whichp, for_plot = TRUE)
    ma <- get_model_args(model)

    # For univariate models, we need the original data with individual columns
    # not the model frame which may have Surv() objects for survival models
    uv_data <- if (!missing(data)) data else get_model_data(model)

    tabUV <- rm_uvsum(response = ma$response, covs = ma$predictors, data = uv_data,
                      tableOnly = TRUE, nicenames = FALSE, unformattedp = unformattedp)
    tab <- combine_uv_mv(tabUV, m_sum, tab)
  }
  if (!showN) {
    rmc <- which(names(tab)=="N")
    if (length(rmc)>0) tab <- tab[,-rmc ]
  }
  if (!showEvent) {
    rmc <- which(names(tab)=="Event")
    if (length(rmc)>0) tab <- tab[,-rmc ]
  }
  att_tab <- attributes(tab)
  model_terms <- if (is_mira) fit1$terms else model$terms
  to_indent <- which(!(tab[["Variable"]] %in% attr(model_terms, "term.labels")))
  bold_cells <- cbind(which(tab[["Variable"]] %in% attr(model_terms, "term.labels")), rep(1, length(which(tab[["Variable"]] %in% attr(model_terms, "term.labels")))))
  attr(tab, "to_indent") <- to_indent
  attr(tab, "bold_cells") <- bold_cells

  # Find all estimate columns (contain "95%CI")
  est_cols <- grep("\\(95%CI\\)", names(tab), value = TRUE)

  # Format and bold p-values
  pv <- format_bold_pvalues(tab, bold_cells,
                            unformattedp = unformattedp, p.adjust = p.adjust)
  tab <- pv$tab; bold_cells <- pv$bold_cells

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

get_model_args <- function(model) {
  av <- as.character(model$call)
  av_f <- av[which(grepl("~",av))]
  f <- as.formula(av_f)

  #  response variables (strip functions like Surv)
  response_vars <- all.vars(f[[2]])

  #  predictor variables (ignore interactions)
  term_labels <- attr(terms(f), "term.labels")

  predictor_vars <- unique(unlist(
    strsplit(term_labels, "[:*]")
  ))

  list(
    response = response_vars,
    predictors = predictor_vars
  )
}


combine_uv_mv <- function(tabUV, m_sum, tabMV) {

  # Detect estimate column name dynamically
  est_col_uv <- grep("\\(95%CI\\)", names(tabUV), value = TRUE)[1]
  est_col_mv <- grep("\\(95%CI\\)", names(tabMV), value = TRUE)[1]

  if (is.na(est_col_uv)) stop("Cannot find estimate column in tabUV")
  if (is.na(est_col_mv)) stop("Cannot find estimate column in tabMV")

  # Standardize column names for internal processing
  tabUV_work <- tabUV
  names(tabUV_work)[names(tabUV_work) == est_col_uv] <- "Est_CI"

  tabMV_work <- tabMV
  names(tabMV_work)[names(tabMV_work) == est_col_mv] <- "Est_CI"

  # 1. Reconstruct (var, lvl) for tabUV using its structure
  tabUV2 <- tabUV_work |>
    dplyr::mutate(
      is_header = (is.na(Est_CI) | Est_CI == "<NA>") &
                  (is.na(`p-value`) | `p-value` == "<NA>"),
      var = dplyr::if_else(is_header, Covariate, NA_character_),
      lvl = dplyr::if_else(!is_header, Covariate, NA_character_)
    ) |>
    tidyr::fill(var, .direction = "down")

  # 2. Attach (var, lvl) to MV table and merge with UV
  # Preserve original tabMV row order
  out <- tabMV_work |>
    dplyr::mutate(
      mv_order = dplyr::row_number(),
      is_header = is.na(Est_CI) | Est_CI == "<NA>",
      var = ifelse(is_header, Variable, NA_character_)
    ) |>
    tidyr::fill(var, .direction = "down") |>
    dplyr::mutate(
      lvl = ifelse(is_header, NA_character_, Variable)
    ) |>
    dplyr::left_join(
      tabUV2 |>
        dplyr::select(dplyr::any_of(c("var", "lvl", "Est_CI", "p-value"))) |>
        dplyr::rename(
          `Unadjusted Est_CI` = Est_CI,
          `Unadjusted p-value` = `p-value`
        ),
      by = c("var", "lvl")
    ) |>
    dplyr::rename(
      `Adjusted Est_CI` = Est_CI,
      `Adjusted p-value` = `p-value`
    ) |>
    dplyr::arrange(mv_order) |>
    dplyr::select(dplyr::any_of(c(
      "Variable", "var",
      "Unadjusted Est_CI", "Unadjusted p-value",
      "Adjusted Est_CI", "Adjusted p-value",
      "N", "Event", "VIF"
    )))

  # 3. Identify UV-only main effects from interactions
  mv_main_vars <- unique(out$var[!grepl(":", out$var)])
  out_cols <- names(out)
  uv_only_vars <- c()

  # Find all interaction terms in MV model and extract their components
  interaction_terms <- unique(out$var[grepl(":", out$var)])
  for (int_term in interaction_terms) {
    components <- unlist(strsplit(int_term, ":"))
    for (comp in components) {
      if (!(comp %in% mv_main_vars) && !(comp %in% uv_only_vars)) {
        uv_only_vars <- c(uv_only_vars, comp)
      }
    }
  }

  # Remove 'var' column before appending UV-only rows
  out <- out |> dplyr::select(-var)

  # 4. Append UV-only main effects at the end
  if (length(uv_only_vars) > 0) {
    uv_rows_to_add <- list()

    for (v in uv_only_vars) {
      uv_rows <- tabUV2 |> dplyr::filter(var == v)
      if (nrow(uv_rows) == 0) next

      has_header <- any(uv_rows$is_header)

      if (has_header) {
        # Categorical: add header then levels
        header_row <- uv_rows |> dplyr::filter(is_header)
        level_rows <- uv_rows |> dplyr::filter(!is_header)

        new_header <- as.list(rep(NA, length(out_cols)))
        names(new_header) <- out_cols
        new_header$Variable <- v
        new_header$N <- header_row$N[1]
        uv_rows_to_add[[length(uv_rows_to_add) + 1]] <- as.data.frame(new_header)

        for (k in seq_len(nrow(level_rows))) {
          lv_row <- level_rows[k, ]
          new_row <- as.list(rep(NA, length(out_cols)))
          names(new_row) <- out_cols
          new_row$Variable <- lv_row$lvl
          new_row$`Unadjusted Est_CI` <- lv_row$Est_CI
          new_row$`Unadjusted p-value` <- lv_row$`p-value`
          new_row$N <- lv_row$N
          if ("Event" %in% names(lv_row)) new_row$Event <- lv_row$Event
          uv_rows_to_add[[length(uv_rows_to_add) + 1]] <- as.data.frame(new_row)
        }
      } else {
        # Continuous: single row
        uv_row <- uv_rows[1, ]
        new_row <- as.list(rep(NA, length(out_cols)))
        names(new_row) <- out_cols
        new_row$Variable <- v
        new_row$`Unadjusted Est_CI` <- uv_row$Est_CI
        new_row$`Unadjusted p-value` <- uv_row$`p-value`
        new_row$N <- uv_row$N
        if ("Event" %in% names(uv_row)) new_row$Event <- uv_row$Event
        uv_rows_to_add[[length(uv_rows_to_add) + 1]] <- as.data.frame(new_row)
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

  out
}
