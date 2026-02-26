# Model Registry for reportRmd Package ----
# Internal constants defining model types and their S3 dispatch classes

#' Get model class from model specifications
#'
#' Internal function that maps model type, family, and GEE usage to the
#' appropriate S3 class for autoreg() dispatch.
#'
#' @param type Model type (linear, logistic, poisson, negbin, ordinal, boxcox, coxph, crr)
#' @param family Model family (gaussian, binomial, poisson, or NULL)
#' @param gee Logical indicating if GEE model
#' @return Character string of S3 class name
#' @keywords internal
get_model_class <- function(type, family = NULL, gee = FALSE) {

  # Model registry: maps (type, family, gee) -> (class, beta_label)
  # This replaces the uvmodels data object with a more maintainable structure

  registry <- list(
    # Standard GLM/LM models ----
    list(type = "linear",   family = "gaussian", gee = FALSE,
         class = "rm_lm",     beta = "Estimate"),
    list(type = "logistic", family = "binomial", gee = FALSE,
         class = "rm_glm",    beta = "OR"),
    list(type = "poisson",  family = "poisson",  gee = FALSE,
         class = "rm_glm",    beta = "RR"),
    list(type = "negbin",   family = NA,         gee = FALSE,
         class = "rm_negbin", beta = "RR"),
    list(type = "ordinal",  family = NA,         gee = FALSE,
         class = "rm_ordinal", beta = "OR"),
    list(type = "boxcox",   family = "gaussian", gee = FALSE,
         class = "rm_boxcox", beta = "Estimate"),

    # Survival models ----
    list(type = "coxph",    family = NA,         gee = FALSE,
         class = "rm_coxph",  beta = "HR"),
    list(type = "crr",      family = NA,         gee = FALSE,
         class = "rm_crr",    beta = "HR"),

    # GEE models (limited support) ----
    list(type = "linear",   family = "gaussian", gee = TRUE,
         class = "rm_gee",    beta = "Estimate"),
    list(type = "logistic", family = "binomial", gee = TRUE,
         class = "rm_gee",    beta = "OR"),
    list(type = "poisson",  family = "poisson",  gee = TRUE,
         class = "rm_gee",    beta = "RR")
  )

  # Normalize NA handling
  family_normalized <- if (is.null(family)) NA_character_ else as.character(family)

  # Find matching entry
  matches <- vapply(registry, function(entry) {
    type_match <- entry$type == type
    family_match <- (is.na(entry$family) && is.na(family_normalized)) ||
                    (!is.na(entry$family) && !is.na(family_normalized) &&
                     entry$family == family_normalized)
    gee_match <- entry$gee == gee

    type_match && family_match && gee_match
  }, logical(1))

  if (!any(matches)) {
    # Check if it's an unsupported GEE combination
    if (gee && type %in% c("negbin", "ordinal", "boxcox", "coxph", "crr")) {
      stop(paste0(
        "GEE models not supported for type '", type, "'. ",
        "GEE is only available for linear, logistic, and poisson regression."
      ))
    }

    stop(paste0(
      "Could not find matching model class for: ",
      "type = '", type, "', ",
      "family = '", ifelse(is.null(family), "NULL", family), "', ",
      "gee = ", gee, ". ",
      "Please check your model specification."
    ))
  }

  matched_entry <- registry[[which(matches)[1]]]

  return(list(
    class = matched_entry$class,
    beta_label = matched_entry$beta
  ))
}

#' Get beta label for model type
#'
#' Returns the appropriate coefficient label (OR, HR, RR, Estimate)
#' for a given model class.
#'
#' @param model_class S3 class name (e.g., "rm_glm", "rm_coxph")
#' @return Character string for beta label
#' @keywords internal
get_beta_label <- function(model_class) {
  labels <- c(
    rm_lm = "Estimate",
    rm_glm = "OR",      # Default for rm_glm, context-dependent
    rm_negbin = "RR",
    rm_ordinal = "OR",
    rm_boxcox = "Estimate",
    rm_coxph = "HR",
    rm_crr = "HR",
    rm_gee = "Estimate" # Default for rm_gee, context-dependent
  )

  labels[model_class]
}
