# Helper Functions for reportRmd Package ----

# Column/String Utilities ----

#' Convert Excel column letters to numbers
#'
#' Converts Excel-style column letters (A, B, ..., Z, AA, AB, ...) to
#' numeric column indices.
#'
#' @param v Vector of column letter codes
#' @return Numeric vector of column indices
#' @keywords internal
#' @examples
#' \dontrun{
#' xcn("A")   # Returns 1
#' xcn("Z")   # Returns 26
#' xcn("AA")  # Returns 27
#' }
xcn <- function(v) {
  sapply(v, function(x) {
    col_head <- toupper(x)
    if (nchar(col_head) > 1) {
      l1 <- substr(col_head, 1, 1)
      l2 <- substr(col_head, 2, 2)
      rtn <- 26 * which(LETTERS == l1) + which(LETTERS == l2)
    } else {
      rtn <- which(LETTERS == col_head)
    }
    return(rtn)
  })
}

#' Check if object is an error
#'
#' @param x Object to check
#' @return Logical indicating if x inherits from "try-error"
#' @keywords internal
is.error <- function(x) {
  inherits(x, "try-error")
}

#' Column separator for paste operations
#'
#' @return Character string ", "
#' @keywords internal
csep <- function() {
  return(", ")
}

.negloglik.boxcox <- function (lambda.val, data, xmat, lik.method = "ML")
{
  if (length(lambda.val) == 2) {
    data <- data + lambda.val[2]
    lambda <- lambda.val[1]
  }
  else lambda <- lambda.val
  lambda <- unname(lambda)
  n <- length(data)
  beta.size <- ncol(xmat)
  if (isTRUE(all.equal(unname(lambda), 0)))
    yt <- log(data)
  else yt <- ((data^lambda) - 1)/lambda
  beta <- solve(crossprod(xmat), crossprod(xmat, yt))
  ss <- sum((drop(yt) - drop(xmat %*% beta))^2)
  if (lik.method == "ML")
    neglik <- (n/2) * log(ss) - ((lambda - 1) * sum(log(data)))
  if (lik.method == "RML") {
    xx <- crossprod(xmat)
    if (length(as.vector(xx)) == 1)
      choldet <- 0.5 * log(xx)
    else choldet <- sum(log(diag(chol(xx))))
    neglik <- ((n - beta.size)/2) * log(ss) + choldet - ((lambda -
                                                            1) * sum(log(data)))
  }
  if (mode(neglik) != "numeric")
    neglik <- Inf
  return(drop(neglik))
}

# REMOVED: Duplicate fillNAs() definition
# The vectorized version at line 710 is kept for better performance



# Number and String Formatting Functions ----

#' Round numbers with trailing zeros
#'
#' Rounds numeric values to specified decimal places while preserving
#' trailing zeros (e.g., 1.5 -> "1.50" with digits=2).
#'
#' @param x Numeric vector to round
#' @param digits Number of decimal places (default 2)
#' @return Character vector of rounded numbers with trailing zeros
#' @keywords helper
#' @examples
#' \dontrun{
#' niceNum(c(1.5, 2.345, NA), digits = 2)
#' # Returns: c("1.50", "2.35", NA)
#' }
niceNum <- function(x, digits = 2) {
  rndx <- sapply(x, function(x) {
    if (is.na(x)) return(x)
    if (is.null(x)) return(x)
    format(round(as.numeric(x), digits), nsmall = digits)
  })
  return(gsub(" ", "", rndx))
}



#' Paste vector elements with parentheses
#'
#' Formats a vector as "first (second, third, ...)" where remaining
#' elements are comma-separated inside parentheses.
#'
#' @param x Vector of values (first element shown separately)
#' @return Character string with first element followed by remaining elements in parentheses
#' @keywords helper
#' @examples
#' \dontrun{
#' pstprn(c(1.5, 1.2, 1.8))
#' # Returns: "1.5 (1.2, 1.8)"
#' }
pstprn <- function(x) {
  paste0(x[1], " (", paste(x[-1], collapse = csep()), ")")
}

#' Round and paste with parentheses (smart formatting)
#'
#' Rounds numeric values and formats as "value (lower, upper)" with
#' intelligent formatting:
#' \itemize{
#'   \item Values with |x| < 0.01 or |x| > 1000: scientific notation
#'   \item Other values: standard rounding with trailing zeros
#' }
#'
#' @param x Numeric vector to round and format
#' @param y Number of digits/significant figures (default 2)
#' @return Character string with first element followed by remaining elements in parentheses
#' @keywords helper
#' @examples
#' \dontrun{
#' psthr(c(1.234, 1.123, 1.345), 2)
#' # Returns: "1.23 (1.12, 1.35)"
#' psthr(c(0.001, 0.0005, 0.0015), 2)
#' # Returns: "1.0e-03 (5.0e-04, 1.5e-03)"
#' }
psthr <- function(x, y = 2) {
  x <- sapply(x, function(x) {
    ifelse(abs(x) < 0.01 | abs(x) > 1000,
           format(x, scientific = TRUE, digits = y),
           format(round(x, y), nsmall = y))
  })
  pstprn(x)
}

# Model Coefficient and Term Matching Functions ----

#' Match coefficient names to covariate names
#'
#' Matches model coefficient names (with factor levels) back to original
#' covariate names from the model call. Handles cases where one factor
#' name is a subset of another.
#'
#' @param betanames Vector of coefficient names from model
#' @param call Vector of covariate names from model formula
#' @return Vector of matched covariate names
#' @keywords internal
covnm <- function(betanames, call) {
  sapply(betanames, function(betaname) {
    # Find indices where call elements are found in betaname
    # Changed from charmatch to grepl on Feb 21, 2019
    indx <- which(sapply(call, function(cov) grepl(cov, betaname, fixed = TRUE)))

    if (length(indx) == 1) return(call[indx])

    # If one factor name is a subset of another, choose longest match
    indx2 <- which.max(sapply(call[indx], nchar))
    if (length(indx2) == 1) return(call[indx[indx2]])

    # Check if betaname starts with the matched covariate
    indx3 <- which(sapply(call[indx2], function(c) {
      substr(betaname, 1, nchar(c)) == c
    }))
    if (length(indx3) == 1) return(call[indx[indx2[indx3]]])
  })
}

#' Check if all elements are equal
#'
#' @param x First vector/list
#' @param y Second vector/list
#' @return Logical indicating if all elements are equal
#' @keywords internal
alleql <- function(x, y) {
  !any((x == y) == FALSE)
}

#' Group sequential elements with same values
#'
#' Creates index groups for sequential elements that have the same values.
#' Used for grouping factor levels in model coefficients.
#'
#' @param x List of vectors to group
#' @return List of index vectors indicating groups
#' @keywords internal
betaindx <- function(x) {
  i <- 1
  out <- 1
  result <- NULL

  while (TRUE) {
    if (i + 1 > length(x)) {
      result <- c(result, list(out))
      return(result)
    } else if (alleql(x[[i + 1]], x[[i]])) {
      out <- c(out, i + 1)
    } else {
      result <- c(result, list(out))
      out <- i + 1
    }
    i <- i + 1
  }
}


#' Capitalize a string
#'
#' Capitalize a string
#'
#' @param x string
#' @keywords helper
cap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}

#' Lean strings for printing
#'
#' Returns strings with . and _ replaced by a space. This is nice when printing
#' column names of your dataframe in a report
#' @param strings vector of strings to give a nice name
#' @param check_numbers boolean indicating if numbers with decimals should be
#'   checked for and retained.
#' @keywords helper
nicename <-function (strings,check_numbers=TRUE)
{
  out <- sapply(strings, function(x) {
    original_x <- x
    x <- chartr(".", " ", x)
    x <- chartr("_", " ", x)
    if(check_numbers){
      p.positions <- gregexpr(pattern ='\\d\\.[0-9]+',original_x)[[1]]+1
      for(pos in p.positions){
        substr(x,pos,pos) <- '.'
      }

    }
    x <- gsub(" +", " ", x)
    return(x)
  })
  return(out)
}


# NOTE: pvalue() function removed - it was not used anywhere in the codebase
# Use formatp() for general p-value formatting instead

# P-value Formatting Functions ----
# These functions format p-values for different output contexts

#' Format p-values for tables (HTML/LaTeX)
#'
#' Standard p-value formatting for table output:
#' - p < 0.001: returns "<0.001"
#' - p < 0.1: returns p to 3 decimal places
#' - p >= 0.1: returns p to 2 decimal places
#'
#' This is the main p-value formatting function used throughout the package.
#' Used by: rm_compactsum(), rm_mvsum(), rm_uvsum()
#'
#' @param pvalues Vector of p-values (numeric or character)
#' @return Character vector of formatted p-values
#' @keywords helper
#' @examples
#' \dontrun{
#' formatp(c(0.0001, 0.045, 0.123))
#' # Returns: c("<0.001", "0.045", "0.12")
#' }
formatp <- function(pvalues) {
  p_out <- sapply(pvalues, function(x) {
    xsig <- suppressWarnings(as.numeric(x))
    fmt_x <- ifelse(xsig < 0.001, "<0.001",
                    ifelse(xsig < 0.1,
                           format(round(xsig, 3), nsmall = 3),
                           format(round(xsig, 2), nsmall = 2)))
    x <- ifelse(x == "excl", "excl", fmt_x)
    return(x)
  })
  p_out <- unname(p_out)
  return(p_out)
}


# LaTeX and HTML Formatting Functions ----

#' Sanitize strings for LaTeX output
#'
#' Escapes special LaTeX characters to prevent compilation errors.
#' This is an internal function called by sanitizestr().
#'
#' Special characters handled: dollar sign, ampersand, percent, hash,
#' underscore, braces, tilde, caret, angle brackets, pipe, backslash
#'
#' @param str Character string to sanitize
#' @return Sanitized string safe for LaTeX
#' @keywords internal
sanitize <- function(str) {
  result <- str
  # Temporarily replace backslash to avoid double-escaping
  result <- gsub("\\\\", "SANITIZE.BACKSLASH", result)
  result <- gsub("$", "\\$", result, fixed = TRUE)
  result <- gsub(">", "$>$", result, fixed = TRUE)
  result <- gsub("<", "$<$", result, fixed = TRUE)
  result <- gsub("|", "$|$", result, fixed = TRUE)
  result <- gsub("{", "\\{", result, fixed = TRUE)
  result <- gsub("}", "\\}", result, fixed = TRUE)
  result <- gsub("%", "\\%", result, fixed = TRUE)
  result <- gsub("&", "\\&", result, fixed = TRUE)
  result <- gsub("_", "\\_", result, fixed = TRUE)
  result <- gsub("#", "\\#", result, fixed = TRUE)
  result <- gsub("^", "\\verb|^|", result, fixed = TRUE)
  result <- gsub("~", "\\~{}", result, fixed = TRUE)
  # Restore backslash
  result <- gsub("SANITIZE.BACKSLASH", "$\\backslash$", result, fixed = TRUE)
  return(result)
}

#' Sanitize string vector for LaTeX
#'
#' Vectorized wrapper around sanitize() that processes multiple strings.
#' Strings with special characters will break LaTeX if returned 'asis' by
#' knitr. This function is called by main reportRmd functions before output
#' to ensure LaTeX-safe strings.
#'
#' @param str Vector of strings to sanitize
#' @return Vector of LaTeX-safe strings
#' @keywords helper
#' @seealso \code{\link{sanitize}}
sanitizestr <- function(str) {
  as.vector(sapply(str, function(char) {
    sanitize(char)
  }, USE.NAMES = FALSE))
}

#' Bold strings for LaTeX output
#'
#' Wraps strings in LaTeX bold formatting (\\textbf{}).
#'
#' @param strings Vector of strings to bold
#' @return Vector of strings wrapped in \\textbf{}
#' @keywords helper
lbld <- function(strings) {
  sapply(strings, function(x) {
    if (is.null(x)) return(x)
    if (is.na(x)) return(x)
    return(paste0("\\textbf{", x, "}"))
  }, USE.NAMES = FALSE)
}

#' Bold strings for HTML output
#'
#' Wraps strings in HTML bold formatting using inline CSS.
#'
#' @param strings Vector of strings to bold
#' @return Vector of strings wrapped in HTML bold span
#' @keywords helper
hbld <- function(strings) {
  sapply(strings, function(x) {
    if (is.null(x)) return(x)
    if (is.na(x)) return(x)
    return(paste0('<span style="font-weight: bold;">', x, "</span>"))
  }, USE.NAMES = FALSE)
}

#' Escape less-than and greater-than signs for HTML
#'
#' Replaces < and > with their HTML entities (&lt; and &gt;).
#' Handles non-ASCII characters by attempting conversion.
#'
#' @param x Vector of strings to process
#' @return Vector with < and > replaced by HTML entities
#' @keywords internal
ltgt <- function(x) {
  sapply(x, function(x) {
    z <- try(gsub("<", "&lt;", x), silent = TRUE)
    if (inherits(z, "try-error")) {
      # Warn user and try to convert non-ASCII characters
      warning(paste0(
        "The following string contains non-ASCII characters ",
        "and may not display properly:\n", x
      ))
      z <- try(gsub("<", "&lt;", iconv(x, to = "ASCII", sub = "")),
               silent = TRUE)
      if (inherits(z, "try-error")) return(NA)
    }
    z <- gsub(">", "&gt;", z)
    return(z)
  }, USE.NAMES = FALSE)
}

#' Replace dollar signs for HTML output
#'
#' Escapes dollar signs and < > symbols for proper HTML rendering.
#' Converts $ to HTML entity and calls ltgt() for angle brackets.
#'
#' @param s Character vector to process
#' @return Character vector with escaped special characters
#' @keywords helper
rmds <- function(s) {
  sapply(s, function(x) {
    x <- ltgt(x)  # Handle < and >
    gsub("[$]", '<span style="display: inline">&#36</span>', x)
  }, USE.NAMES = FALSE)
}

#' Add LaTeX spacing before string
#'
#' Prepends "~~~" (LaTeX non-breaking spaces) before string.
#' Used for indentation in LaTeX tables.
#'
#' @param x Character string
#' @return String with "~~~" prepended
#' @keywords helper
addspace <- function(x) {
  paste0("~~~", x)
}
#' Format p-values for LaTeX output with bolding
#'
#' Formats p-values for LaTeX documents with automatic bolding of
#' significant values (p <= 0.05).
#'
#' Formatting rules:
#' \itemize{
#'   \item p <= 0.001: returns bold "< 0.001" in LaTeX format
#'   \item p <= 0.05: returns p rounded to specified digits, bolded
#'   \item p > 0.05: returns p rounded to specified digits, not bolded
#'   \item p <= 0.1: uses 3 decimal places
#'   \item p > 0.1: uses sigdigits decimal places
#' }
#'
#' Used by: rm_covsum() and related LaTeX output functions in main.R
#'
#' @param x Numeric p-value
#' @param sigdigits Number of significant digits to report (default 2)
#' @return Character string with LaTeX formatting
#' @keywords helper
#' @examples
#' \dontrun{
#' lpvalue(0.0001)  # Returns bold < 0.001 in LaTeX
#' lpvalue(0.023)   # Returns: "\\textbf{0.023}"
#' lpvalue(0.156)   # Returns: "0.16"
#' }
lpvalue <- function(x, sigdigits = 2) {
  if (is.na(x) || inherits(x, "character")) {
    return(x)
  } else if (x <= 0.001) {
    return("\\textbf{$<$0.001}")
  } else if (x <= 0.1) {
    x <- format(round(x, 3), nsmall = 3)
  } else {
    x <- format(round(x, sigdigits), nsmall = sigdigits)
  }

  if (as.numeric(x) <= 0.05) {
    return(paste0("\\textbf{", x, "}"))
  } else {
    return(x)
  }
}


#' Remove dollar sign prefix from column names
#'
#' Removes common prefix (before $) from interaction term column names.
#' Used to clean up model matrix column names.
#'
#' @param x Character vector of column names
#' @return Character vector with dollar sign prefixes removed
#' @keywords internal
removedollar <- function(x) {
  colnms <- strsplit(x, ":")
  indx <- unlist(lapply(colnms, function(colnm) {
    sapply(colnm, function(coln) {
      regexpr("$", coln, fixed = TRUE)[1] + 1
    })
  }))

  if (length(unique(indx)) == 1) {
    if (unique(indx) != 0) {
      x <- unlist(lapply(colnms, function(colnm) {
        paste(substring(colnm, indx[1]), collapse = ":")
      }))
    }
  }
  return(x)
}

#' Create model matrix from formula
#'
#' Extracts model matrix from formula, optionally separating response
#' variables from predictors. Removes intercept column and cleans
#' column names.
#'
#' @param f Formula object
#' @param data Data frame (optional)
#' @return Model matrix or list of matrices (y and x) if response present
#' @keywords internal
modelmatrix <- function(f, data = NULL) {
  k <- as.character(f)
  y <- NULL

  if (!length(k) %in% c(2, 3)) {
    stop("Formula not properly formed")
  }

  if (length(k) == 3) {
    # Formula has response variable
    f <- stats::as.formula(paste0("~", k[2], "+", k[3]))
    y <- stats::model.matrix(
      stats::as.formula(paste0("~", k[2])),
      data
    )[, -1, drop = FALSE]
  }

  x <- stats::model.matrix(f, data)[, -1, drop = FALSE]
  colnames(x) <- removedollar(colnames(x))

  if (!is.null(y)) {
    return(list(
      x[, 1:ncol(y), drop = FALSE],
      x[, (ncol(y) + 1):ncol(x), drop = FALSE]
    ))
  } else {
    return(x)
  }
}

#' Format model call as clean string
#'
#' Converts model call to string with single quotes instead of double quotes.
#'
#' @param model_call Call object from model
#' @return Character string of formatted call
#' @keywords internal
nicecall <- function(model_call) {
  call_str <- base::deparse1(model_call)
  call_str <- gsub('[\"]', "'", call_str)
  return(call_str)
}

#' Extract data frame name from model call argument
#'
#' Attempts to identify the data frame object used in a model call
#' by searching the global environment.
#'
#' @param dataArg Data argument from model call
#' @return Character string of data frame name, or NULL if not found
#' @keywords internal
matchdata <- function(dataArg) {
  df_str <- as.character(dataArg)
  if (length(df_str) > 1) df_str <- df_str[2]

  # Remove function calls
  no_fnc <- gsub("[A-Za-z]+[(]", "", df_str)

  # Extract valid R object names
  txt_bts <- unlist(strsplit(no_fnc, split = "[^A-Za-z0-9_.]"))
  txt_bts <- txt_bts[txt_bts != ""]

  # Find matching objects in global environment
  obj <- intersect(txt_bts, ls(name = ".GlobalEnv"))

  if (length(obj) > 0) {
    # Check which are data frames
    df_ind <- sapply(obj, function(x) {
      inherits(get0(x), "data.frame")
    })
    df <- obj[df_ind]

    if (length(df) > 1) {
      message("Multiple data objects found in function call")
      return(NULL)
    } else {
      return(df)
    }
  } else {
    message("Data object could not be extracted from function call")
    return(NULL)
  }
}

#' Match coefficient names to covariate indices
#'
#' Matches model coefficient names (including interactions) to original
#' covariate indices from the model call. Handles centered variables and
#' interaction terms. Returns a numeric encoding for sorting.
#'
#' Changes:
#' - Feb 21, 2019: Changed from charmatch to grepl for more reliable matching
#' - Dec 14, 2020: Added space removal to support centered variables
#'
#' @param betanames Vector of coefficient names from model
#' @param ucall Vector of unique covariate names from model call
#' @return Numeric vector of encoded covariate indices, or -1 if matching fails
#' @keywords internal
matchcovariate <- function(betanames, ucall) {
  out <- as.vector(sapply(betanames, function(betaname) {
    # Split interaction terms
    splitbetaname <- unlist(strsplit(betaname, ":", fixed = TRUE))

    out <- sapply(splitbetaname, function(bname) {
      # Remove spaces to handle centered variables (e.g., "I(x - mean(x))")
      bname <- gsub(" ", "", bname)

      # Find matching covariate indices
      indx <- which(sapply(ucall, function(cov) {
        grepl(cov, bname, fixed = TRUE)
      }))

      if (length(indx) == 1) return(indx)

      # If one factor name is a subset of another, choose longest match
      indx2 <- which.max(sapply(ucall[indx], nchar))
      if (length(indx2) == 1) return(indx[indx2])

      # Check if betaname starts with matched covariate
      indx3 <- which(sapply(ucall[indx2], function(c) {
        substr(betaname, 1, nchar(c)) == c
      }))
      if (length(indx3) == 1) return(ucall[indx[indx2[indx3]]])

      return(-1)
    })

    if (-1 %in% out) return(-1)

    # Encode indices as base-100 number for sorting
    result <- 0
    n <- length(out)
    for (i in seq_along(out)) {
      result <- result + out[i] * 100^(n - 1)
      n <- n - 1
    }
    return(result)
  }))

  if (-1 %in% out) return(-1)
  return(out)
}

# Statistical Utility Functions ----

#' Compute Generalized Variance Inflation Factor (GVIF)
#'
#' Calculates GVIF for model terms to detect multicollinearity.
#' Adapted from the car package.
#'
#' For terms with degrees of freedom > 1 (e.g., factor variables),
#' returns GVIF^(1/(2*df)) for comparability with standard VIF.
#'
#' @param model Fitted model object (lm, glm, coxph, etc.)
#' @return Data frame with columns: Covariate (term name) and VIF (value)
#' @keywords internal
#' @references
#' Fox, J. and Weisberg, S. (2019). An R Companion to Applied Regression,
#' Third Edition. Thousand Oaks CA: Sage.
GVIF <- function(model) {
  v <- stats::vcov(model)
  ind <- attr(stats::model.matrix(model), "assign")

  # Remove intercept if present
  if (0 %in% ind) {
    v <- v[-1, -1]
    ind <- ind[-1]
  }

  xvar <- labels(stats::terms(model))

  # Return NA for models with < 2 predictors
  if (length(xvar) < 2) {
    return(data.frame(Covariate = xvar, VIF = NA))
  }

  R <- stats::cov2cor(v)
  detR <- det(R)
  result <- matrix(0, length(xvar), 2)

  for (var in seq_along(xvar)) {
    terms <- which(ind == var)
    result[var, 1] <- det(as.matrix(R[terms, terms])) *
      det(as.matrix(R[-terms, -terms])) / detR
    result[var, 2] <- length(terms)
  }

  # Adjust GVIF for terms with df > 1
  if (all(result[, 2] == 1)) {
    rtn <- result[, 1]
  } else {
    rtn <- result[, 1]^(1 / (2 * result[, 2]))
  }

  data.frame(Covariate = xvar, VIF = rtn)
}

# (ggsurv) ---------------------------------------------------------

# Survival Analysis Formatting Functions ----
# These variants use no spaces and are specific to survival/ggkmcif plots

#' Round with sprintf formatting
#'
#' Rounds values using sprintf for precise decimal formatting.
#' Used internally by psthr0().
#'
#' @param value Numeric value to round
#' @param digits Number of decimal places
#' @return Character string with exactly 'digits' decimal places
#' @keywords internal
round_sprintf <- function(value, digits) {
  sprintf(paste0("%.", digits, "f"), round(value, digits))
}

#' Paste with parentheses (no spaces)
#'
#' Similar to pstprn() but without spaces. Used for compact
#' formatting in survival plots.
#'
#' @param x Vector of values
#' @return Character string in format "x[1](x[2],x[3],...)"
#' @keywords internal
pstprn0 <- function(x) {
  paste0(x[1], "(", paste0(x[-1], collapse = ","), ")")
}

#' Round and paste with parentheses (no spaces, for plots)
#'
#' Compact version of psthr() without spaces. Used specifically
#' for survival plot annotations where space is limited.
#'
#' @param x Numeric vector to round and format
#' @param digits Number of decimal places (default 2)
#' @return Character string in format "x[1](x[2],x[3])"
#' @keywords internal
psthr0 <- function(x, digits = 2) {
  x <- sapply(x, function(x) {
    ifelse(abs(x) < 0.01 | abs(x) > 1000,
           format(x, scientific = TRUE, digits = digits),
           round_sprintf(x, digits))
  })
  pstprn0(x)
}

#' Custom break function for axis scaling
#'
#' Creates axis breaks based on maximum value with custom logic for
#' determining break intervals based on magnitude.
#'
#' @param xmax Maximum value for axis
#' @return Numeric vector of break positions
#' @keywords internal
break_function_custom <- function(xmax) {
  xmax_length <- ifelse(xmax > 1,
                        nchar(round(xmax)),
                        round(abs(log10(xmax))))

  byx <- if (xmax > 1) {
    round(xmax / 10, digits = 2 - xmax_length)
  } else {
    round(xmax / 10, digits = xmax_length + 1)
  }

  breaks <- seq(0, xmax, by = byx)
  if (max(breaks) < byx) breaks <- c(breaks, max(breaks) + byx)
  return(breaks)
}

#' Standard break function using pretty()
#'
#' Creates nice axis breaks using R's built-in pretty() function.
#' This is the standard approach used by ggkmcif functions.
#'
#' @param x Maximum value for axis
#' @param n Desired number of intervals (default 5)
#' @return Numeric vector of break positions
#' @keywords internal
break_function <- function(x, n = 5) {
  pretty(c(0, x), n = n)
}

#' Format p-values for plot annotations
#'
#' Formats p-values specifically for display in plots (e.g., survival curves).
#' Returns formatted string with "p = " or "p < " prefix.
#'
#' Formatting rules:
#' \itemize{
#'   \item p < 10^-digits: returns "p < threshold" (e.g., "p < 0.001")
#'   \item p >= threshold: returns "p = value" rounded to specified digits
#' }
#'
#' Used by: ggkmcif2() for survival curve annotations in main.R and ggkmcif3.R
#'
#' @param x Numeric p-value
#' @param digits Number of decimal places to display (default from context)
#' @return Character string with "p = " or "p < " prefix
#' @keywords helper
#' @examples
#' \dontrun{
#' lpvalue2(0.0001, 3)  # Returns: "p < 0.001"
#' lpvalue2(0.0456, 3)  # Returns: "p = 0.046"
#' }
lpvalue2 <- function(x, digits) {
  if (is.na(x) || inherits(x, "character")) {
    return(x)
  } else if (x < 10^-(digits)) {
    return(paste0("p < ", 10^-(digits)))
  } else {
    return(paste0("p = ", round_sprintf(x, digits)))
  }
}

.extract_ggplot_colors <- function(p, grp.levels){
  g <- ggplot2::ggplot_build(p)
  .cols <- unlist(unique(g$data[[1]]["colour"]))
  if(!is.null(grp.levels)){
    if(length(.cols)==1) .cols <- rep(.cols, length(grp.levels))
    names(.cols) <- grp.levels
  }
  .cols
}

.set_large_dash_as_ytext <- function(ggp){
  ggp + ggplot2::theme(
    axis.text.y = ggplot2::element_text(size = 50, vjust = 0.35),
    axis.ticks.y = ggplot2::element_blank()
  )
}

##This function is used by the survfit package
survfit_confint <- function(p, se, logse=TRUE, conf.type, conf.int=0.95,
                            selow, ulimit=TRUE) {
  zval <- stats::qnorm(1- (1-conf.int)/2, 0, 1)
  if (missing(selow)) scale <- 1.0
  else scale <- ifelse(selow==0, 1.0, selow/se)  # avoid 0/0 at the origin
  if (!logse) se <- ifelse(se==0, 0, se/p)   # se of log(survival) = log(p)

  if (conf.type=='plain') {
    se2 <- se* p * zval  # matches equation 4.3.1 in Klein & Moeschberger
    if (ulimit) list(lower= pmax(p -se2*scale, 0), upper = pmin(p + se2, 1))
    else  list(lower= pmax(p -se2*scale, 0), upper = p + se2)
  }
  else if (conf.type=='log') {
    #avoid some "log(0)" messages
    xx <- ifelse(p==0, NA, p)
    se2 <- zval* se
    temp1 <- exp(log(xx) - se2*scale)
    temp2 <- exp(log(xx) + se2)
    if (ulimit) list(lower= temp1, upper= pmin(temp2, 1))
    else  list(lower= temp1, upper= temp2)
  }
  else if (conf.type=='log-log') {
    xx <- ifelse(p==0 | p==1, NA, p)
    se2 <- zval * se/log(xx)
    temp1 <- exp(-exp(log(-log(xx)) - se2*scale))
    temp2 <- exp(-exp(log(-log(xx)) + se2))
    list(lower = temp1 , upper = temp2)
  }
  else if (conf.type=='logit') {
    xx <- ifelse(p==0, NA, p)  # avoid log(0) messages
    se2 <- zval * se *(1 + xx/(1-xx))

    temp1 <- 1- 1/(1+exp(log(p/(1-p)) - se2*scale))
    temp2 <- 1- 1/(1+exp(log(p/(1-p)) + se2))
    list(lower = temp1, upper=temp2)
  }
  else if (conf.type=="arcsin") {
    xx <- ifelse(p==0, NA, p)
    se2 <- .5 *zval*se * sqrt(xx/(1-xx))
    list(lower= (sin(pmax(0, asin(sqrt(xx)) - se2*scale)))^2,
         upper= (sin(pmin(pi/2, asin(sqrt(xx)) + se2)))^2)
  }
  else stop("invalid conf.int type")
}


color_palette_surv_ggplot <- function(length){
  if(length==1) return("black")
  if(length==2) return(c("#D53E4F","#3288BD"))
  if(length==3) return(c("#D53E4F","#ABDDA4","#3288BD"))
  if(length==4) return(c("#D53E4F","#FDAE61","#ABDDA4","#3288BD"))
  if(length==5) return(c("#D53E4F","#FDAE61","#FEE08B","#ABDDA4","#3288BD"))
  if(length==6) return(c("#D53E4F","#FDAE61","#FEE08B","#ABDDA4","#66C2A5","#3288BD"))
  if(length==7) return(c("#D53E4F","#F46D43","#FDAE61","#FEE08B","#ABDDA4","#66C2A5","#3288BD"))
  if(length==8) return(c("#D53E4F","#F46D43","#FDAE61","#FEE08B","#ABDDA4","#66C2A5","#3288BD","#5E4FA2"))
  if(length==9) return(c("#9E0142","#D53E4F","#F46D43","#FDAE61","#FEE08B","#ABDDA4","#66C2A5","#3288BD","#5E4FA2"))
  if(length==10) return(c("black","#9E0142","#D53E4F","#F46D43","#FDAE61","#FEE08B","#ABDDA4","#66C2A5","#3288BD","#5E4FA2"))
  if(length>10) {message("10 colours maximum in default")}
  return(rep(c("black","#9E0142","#D53E4F","#F46D43","#FDAE61","#FEE08B","#ABDDA4","#66C2A5","#3288BD","#5E4FA2"),length.out=length))
}


# (forestplot2) ---------------------------------------------------------
format_glm = function(glm_fit,conf.level = 0.95,digits=c(2,3),orderByRisk=TRUE){
  if (! class(glm_fit)[1] %in% c('glm','geeglm','polr')) stop('Only objects of class glm, geeglm and polr are accepted.')

  #extracting ORs and p values
  Z = stats::qnorm(1-(1-conf.level)/2)
  tab <- as.data.frame(summary(glm_fit)$coefficients)
  tab <- cbind(variable= rownames(tab),tab)
  rownames(tab) <- NULL

  if (class(glm_fit)[1] %in% c("glm", "geeglm")){
    names(tab) =  c("variable","estimate",  "std.error" ,"statistic", "p.value")
    tab = tab[-which(tab$variable=='(Intercept)'),]
  }  else {
    names(tab) =  c("variable","estimate",  "std.error" ,"statistic")
    tab$coef.type = ifelse(grepl("[|]",tab$variable),"scale","coefficient")
    tab <- tab[tab$coef.type=='coefficient',]
    tab$p.value = stats::pnorm(abs(tab$statistic),lower.tail = FALSE) * 2
  }

  tab$conf.low=exp(tab$estimate-Z*tab$std.error)
  tab$conf.high=exp(tab$estimate+Z*tab$std.error)
  tab$estimate = exp(tab$estimate)
  tab$estimate.label = paste0(niceNum(tab$estimate), ' (',niceNum(tab$conf.low),', ',niceNum(tab$conf.high),')')


  tab$p.label = ifelse(tab$p.value<0.001, '<0.001', niceNum(tab$p.value,digits[2]))
  names(tab)[1] = 'variable'

  tab = tab[,c('variable', 'estimate', 'p.label', 'p.value', 'conf.low', 'conf.high')]


  if (orderByRisk){
    tab$var.order = rank(tab$estimate)
  } else{
    tab$var.order = 1:nrow(tab)
  }

  # Extract the reference levels if needed
  if (length(glm_fit$xlevels)!=0){
    ref_levels <- NULL
    for (i in seq_along(glm_fit$xlevels)){
      ref_levels <- rbind(ref_levels,
                          data.frame(var.name=rep(names(glm_fit$xlevels)[i],length(glm_fit$xlevels[[i]])+1),
                                     level.name = c(names(glm_fit$xlevels)[i],glm_fit$xlevels[[i]]),
                                     level.order=1:(length(glm_fit$xlevels[[i]])+1),
                                     variable=paste0(names(glm_fit$xlevels)[i],c('',glm_fit$xlevels[[i]]))))
    }


    tab = merge(ref_levels, tab, by='variable',all = T)

    tab$estimate.label = ifelse(is.na(tab$estimate), '1.0 (Reference)',
                                paste0(niceNum(tab$estimate), ' (',niceNum(tab$conf.low),', ',niceNum(tab$conf.high),')'))

    varOrders <- tapply(X = tab$var.order,
                        INDEX=tab$var.name,
                        FUN = function(x) min(x,na.rm=T))
    varOrderLookup <- data.frame(var.name=names(varOrders),var.order=varOrders)


    varOrderLookup <- stats::na.omit(tab[,c("var.name","var.order")])

    for (i in 1:nrow(varOrderLookup)){
      tab$var.order[tab$var.name==varOrderLookup$var.name[i]] <- varOrderLookup$var.order[i]
    }

    tab$estimate.label = ifelse(tab$level.name %in% names(glm_fit$xlevels),NA_character_,tab$estimate.label)
    tab[order(tab$var.order,tab$level.order,decreasing=c(F,T)),]
  } else {
    tab$estimate.label = paste0(niceNum(tab$estimate), ' (',niceNum(tab$conf.low),', ',niceNum(tab$conf.high),')')
    tab$level.order=1
    tab$var.name=tab$variable
    tab$level.name=tab$variable
    tab[order(tab$var.order),]
  }

}

# New function to strip centering from a covariate
getvarname = function(betaname){
  sapply(betaname,function(x){
    x = gsub('I[(]','',x)
    x = gsub('[-+].*','',x)
    x = trimws(x)
    return(x)
  })
}

lbl_count <- function(y){
  q75 <- summary(y)[5]
  return(data.frame(y=max(y),  label=paste('n =',length(y))))
}

betaWithCI <-function(betaname,CIwidth=0.95){
  paste0(betaname,"(",100*CIwidth,"%CI)")
}

niceStr <- function (strings)
{
  out <- sapply(strings, function(x) {
    x <- chartr('/',' ',x)
    x <- chartr(".", " ", x)
    x <- chartr("_", " ", x)
    return(x)
  })
  return(out)
}

wrp_lbl <- function(x,width = 10){
  x <- niceStr(x)
  #  strwrap(x,width = width) # doesn't work nicely with spaces
  lst <- strwrap(x,width = width,simplify = F)
  for (i in seq_along(lst)) lst[[i]] <- paste(lst[[i]],collapse='\n')
  unlist(lst)
}


label_wrap_reportRx <- function (width = 25, multi_line = TRUE) {
  fun <- function(labels) {
    labels <- ggplot2::label_value(labels, multi_line = multi_line)
    lapply(labels, function(x) {
      x <- niceStr(x)
      x <- strwrap(x, width = width, simplify = FALSE)
      vapply(x, paste, character(1), collapse = "\n")
    })
  }
  structure(fun, class = "labeller")
}






reportRx_pal <- function(
    direction = 1
) {

  function(n) {
    if (n>10) warning('Ten colour maximum, colours will be recycled.')

    colour_list <- color_palette_surv_ggplot(n)

    colour_list <- unname(unlist(colour_list))
    if (direction >= 0) colour_list else rev(colour_list)
  }
}

scale_colour_reportRx <- function(
    direction = 1,
    ...
) {
  ggplot2::discrete_scale(
    aesthetics = c("colour","fill"),
    scale_name = "reportRx",
    reportRx_pal( direction),
    ...
  )
}

#' Forward fill NA values (vectorized implementation)
#'
#' Efficiently fills NA values by carrying forward the last non-NA value.
#' Uses vectorized operations for better performance than loop-based approaches.
#'
#' @param x Vector with NAs to fill
#' @return Vector with NAs filled
#' @keywords internal
#' @examples
#' \dontrun{
#' fillNAs(c(1, NA, NA, 2, NA, 3))  # Returns: c(1, 1, 1, 2, 2, 3)
#' }
fillNAs <- function(x) {
  ind <- which(!is.na(x))
  if (is.na(x[1])) {
    ind <- c(1, ind)
  }
  rep(x[ind], times = diff(c(ind, length(x) + 1)))
}


#' Extract Function and Package Information from Current Document
#'
#' This function analyses the current file (an R script, Rmd or qmd file) to
#' extract information about all functions called within the code, identifies
#' their associated packages, and returns a summary of packages used with
#' version and citation information.
#'
#' @description The function automatically detects the current R script file
#' (works best in RStudio), parses the code to identify function calls,
#' determines which packages they belong to, and creates a summary of all
#' non-base R packages used in the script. It handles both namespace-qualified
#' function calls (e.g., `dplyr::filter`) and regular function calls, while
#' filtering out base R functions and control structures.
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{package_name}{Character. Name of the package}
#'   \item{functions_called}{Character. Comma-separated list of functions called from this package}
#'   \item{package_version}{Character. Version number of the installed package}
#'   \item{package_citation}{Character. Formatted citation for the package}
#' }
#'
#' @note
#' \itemize{
#'   \item Works best when run from RStudio with an active source file
#'   \item Requires that referenced packages are already loaded/installed
#'   \item Will not detect functions called through indirect methods (e.g., `do.call()`)
#' }
#'
#' @examples
#' \dontrun{
#' # Run this function from within an R script to analyze its dependencies
#' package_info <- extract_package_details()
#' print(package_info)
#' }
#'
#' @seealso \code{\link[utils]{getAnywhere}},
#' \code{\link[utils]{packageVersion}}, \code{\link[utils]{citation}}
#'
#' @importFrom utils getAnywhere packageVersion citation
#' @importFrom  rstudioapi isAvailable getSourceEditorContext
#' @importFrom dplyr distinct group_by summarise mutate
#' @export
extract_package_details <- function() {

  # Get the current file path
  # Try different methods to detect the current file
  get_current_file <- function() {
    # Try rstudioapi first (works in RStudio)
    if (requireNamespace("rstudioapi", quietly = TRUE)) {
      if (rstudioapi::isAvailable()) {
        file_path <- rstudioapi::getSourceEditorContext()$path
        if (!is.null(file_path) && file_path != "") {
          return(file_path)
        }
      }
    }
    # Try to get the currently executing script
    # This works when source() is used
    if (exists("ofile", envir = parent.frame())) {
      return(get("ofile", envir = parent.frame()))
    }

    # Try sys.frame approach
    for (i in sys.nframe():1) {
      if (exists("ofile", envir = sys.frame(i))) {
        return(get("ofile", envir = sys.frame(i)))
      }
    }

    # If all else fails, prompt user
    stop("Could not detect current file. Please ensure you're running this from RStudio with an active file, or use rstudioapi package.")
  }

  try(file_path <- get_current_file(),silent = T)
  if (inherits("try-error",file_path)) stop("Current file can not be identified.")


  # Read the file content
  text <- readLines(file_path, warn = FALSE) |>
    paste(collapse = "\n")

  # Remove HTML comments (<!-- -->)
  text <- gsub("<!--.*?-->", "", text, perl = TRUE)

  # Split into lines for processing
  lines <- strsplit(text, "\n")[[1]]

  # Remove comments (everything after # on each line)
  # But preserve # inside strings
  remove_comments <- function(line) {
    # Simple approach: find # not inside quotes
    in_single_quote <- FALSE
    in_double_quote <- FALSE
    chars <- strsplit(line, "")[[1]]
    result <- character()

    for (i in seq_along(chars)) {
      char <- chars[i]

      # Check for escape sequences
      if (i > 1 && chars[i-1] == "\\") {
        result <- c(result, char)
        next
      }

      # Toggle quote states
      if (char == "'" && !in_double_quote) {
        in_single_quote <- !in_single_quote
      } else if (char == '"' && !in_single_quote) {
        in_double_quote <- !in_double_quote
      } else if (char == "#" && !in_single_quote && !in_double_quote) {
        # Found a comment outside of quotes
        break
      }

      result <- c(result, char)
    }

    paste(result, collapse = "")
  }

  # Apply comment removal to each line
  lines <- sapply(lines, remove_comments, USE.NAMES = FALSE)

  # Combine back into single text
  text <- paste(lines, collapse = "\n")

  # First, remove function definitions to avoid false positives
  # Match "function(" only when preceded by space, =, <-, or start of line
  text_no_func_def <- gsub("(^|\\s|=|<-)function\\s*\\(", "\\1FUNCTION_DEFINITION(", text, perl = TRUE)

  # Initialize results list
  all_functions <- list()

  # Pattern 1: Namespace-qualified function calls (package::function or package:::function)
  ns_pattern <- "\\b([a-zA-Z][a-zA-Z0-9.]*)(:::|::)([a-zA-Z_][a-zA-Z0-9._]*)\\s*\\("
  ns_matches <- gregexpr(ns_pattern, text_no_func_def, perl = TRUE)
  ns_calls <- regmatches(text_no_func_def, ns_matches)[[1]]

  if (length(ns_calls) > 0) {
    # Extract package and function names
    for (call in ns_calls) {
      # Remove the opening parenthesis and whitespace
      call_clean <- gsub("\\s*\\($", "", call)

      # Split by :: or :::
      parts <- strsplit(call_clean, ":::|::")[[1]]
      if (length(parts) == 2) {
        all_functions[[length(all_functions) + 1]] <- list(
          function_name = parts[2],
          package_name = parts[1],
          is_namespaced = TRUE
        )
      }
    }
  }

  # Pattern 2: Regular function calls (not namespace-qualified)
  # Exclude calls that are preceded by :: or :::
  regular_pattern <- "(?<!:)(?<!::|:::)\\b([a-zA-Z_][a-zA-Z0-9._]*)\\s*\\("
  regular_matches <- gregexpr(regular_pattern, text_no_func_def, perl = TRUE)
  regular_calls <- regmatches(text_no_func_def, regular_matches)[[1]]

  if (length(regular_calls) > 0) {
    # Extract just the function names (remove the trailing "(")
    function_names <- gsub("\\s*\\($", "", regular_calls)

    # Remove some common non-function patterns
    control_structures <- c("if", "for", "while", "repeat", "function",
                            "case_when","c","list","data.frame","tibble",
                            "FUNCTION_DEFINITION", "switch")
    function_names <- function_names[!function_names %in% control_structures]

    # Add to list
    for (func_name in function_names) {
      all_functions[[length(all_functions) + 1]] <- list(
        function_name = func_name,
        package_name = NA,
        is_namespaced = FALSE
      )
    }
  }

  # Convert to data frame for easier manipulation
  if (length(all_functions) == 0) {
    return(data.frame(
      function_name = character(),
      package_name = character(),
      package_version = character(),
      stringsAsFactors = FALSE
    ))
  }

  functions_df <- do.call(rbind, lapply(all_functions, as.data.frame))

  # Get unique function-package combinations
  functions_df <- functions_df |> dplyr::distinct()

  # Check if not namespaced functions are base R function, and if so, remove
  is_base_R <- function(func_name){
    sapply(func_name, function(f){
      exists(f, mode = "function", envir = baseenv())
    })
  }

  functions_df <- functions_df |>
    dplyr::mutate(
      package_name = ifelse(is.na(package_name),ifelse(is_base_R(function_name),"base",NA),package_name)) |>
    dplyr::filter(!grepl("base",package_name)) |>
    dplyr::filter(function_name != "extract_package_details")

  # Get package information for non-namespaced functions
  get_function_package <- function(func_names) {
    sapply(func_names,function(func_name) {
      # find the function
      where_found <- getAnywhere(func_name)
      if (length(where_found$where)>0){
        if ( ".GlobalEnv" %in% where_found$where) return("GlobalEnv")
        pkg_list <- unique(gsub("package[:]|namespace[:]","",
                                grep("package|namespace",where_found$where,value = T)))
        return(pkg_list[1])
      } else return("Unknown")
    },USE.NAMES=F,simplify=T)
  }

  functions_df <- functions_df |>
    dplyr::mutate(
      package_name = ifelse(is.na(package_name),
                            get_function_package(function_name),package_name)) |>
    dplyr::filter(!grepl("Unknown",package_name))

  # Function to get package version
  get_package_version <- function(pkg_name) {
    sapply(pkg_name, function(pn){
      tryCatch(
        as.character(packageVersion(pn)),
        error = function(e) NA_character_
      )
    })
  }
  # Function to get package citation
  get_package_citation <- function(pkg_name) {
    sapply(pkg_name, function(pn){
      tryCatch(
        format(citation(pn),style = "text"),
        error = function(e) NA_character_
      )
    })
  }
  packages_df <- functions_df |>
    dplyr::group_by(package_name) |>
    dplyr::summarise(functions_called = paste(unique(function_name),collapse=", ")) |>
    dplyr::ungroup() |>
    dplyr::filter(package_name !="GlobalEnv") |>
    dplyr::add_row(package_name ="utils") |>
    dplyr::mutate(package_version = get_package_version(package_name),
                  package_citation = get_package_citation(package_name)) |>
    dplyr::group_by(package_citation) |>
    dplyr::slice_tail(n=1) |>
    dplyr::ungroup() |>
    dplyr::mutate(package_name = gsub("utils","R",package_name))
  ord <- c((1:nrow(packages_df))[-which(packages_df$package_name=="R")],which(packages_df$package_name=="R"))
  packages_df <-packages_df[ord,]
  return(packages_df)
}

