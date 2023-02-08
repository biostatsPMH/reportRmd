#' Retrieve columns number from spreadsheet columns specified as unquoted letters
#' @param ... unquoted excel column headers (i.e. excelCol(A,CG,AA)) separated by commas
#' @importFrom rlang as_string
#' @returns a numeric vector corresponding to columns in a spreadsheet
#' @export
#' @examples
#' ## Find the column numbers for excel columns AB, CE and BB
#' excelCol(AB,CE,bb)
excelCol<- function(...){
  args <- as.list(match.call())[-1]
  args <-unname(unlist(lapply(args,function(x) {rlang::as_string(x)})))
  if (sum(unlist(lapply(args, function(x) grepl("[^A-Za-z]",x))))>0) {
    stop('Only valid Excel column names can be supplied, separated by commas.')
  }
  rtn<-sapply(args, function(x){
    colHead <- toupper(x)
    if (nchar(colHead)>1){
      l1 = substr(colHead,1,1)
      l2 = substr(colHead,2,2)
      rtn <- 26*which(LETTERS==l1)+which(LETTERS==l2)
    } else {
      rtn <- which(LETTERS==colHead)
    }
  })
  names(rtn) <- toupper(names(rtn))
  return(rtn)
}
#' fit box cox transformed linear model
#'
#' Wrapper function to fit fine and gray competing risk model using function crr
#' from package cmprsk
#'
#' @param f formula for the model. Currently the formula only works by using the
#'   name of the column in a dataframe. It does not work by using $ or []
#'   notation.
#' @param data dataframe containing data
#' @param lambda boolean indicating if you want to output the lamda used in the
#'   boxcox transformation. If so the function will return a list of length 2
#'   with the model as the first element and a vector of length 2 as the second.
#' @returns a list containing the linear model (lm) object and, if requested,
#'   lambda
#' @keywords model
boxcoxfitRx<-function(f,data,lambda=FALSE){
  x<-as.character(f)[3]
  y<-as.character(f)[2]
  time<- gsub("\\s","",unlist(strsplit(y,"+",fixed=TRUE))[1])
  covs<-removedollar(x)
  tempindexboxcoxfitRx<-seq_len(nrow(data))
  df1<-data.frame(tempindexboxcoxfitRx,data)
  f2<-as.formula(paste("tempindexboxcoxfitRx+",y,"~",x))
  temp<-modelmatrix(f2,df1)
  ff<-list(temp[[1]][,-1,drop=F],temp[[2,drop=F]])
  temp<-temp[[1]][,1,drop=F]
  lambda1<-unlist(unlist(geoR_boxcoxfit(ff[[1]],ff[[2]],lambda2=TRUE))[1:2])
  ff[[1]]<-((ff[[1]]+lambda1[2])^lambda1[1]-1)/lambda1[1]
  df<-merge(df1,temp,by="tempindexboxcoxfitRx")[,-1,drop=F]
  df[,time]<-ff[[1]]
  out<-lm(f,data=df)
  out$call<-paste("~",covs)
  if(lambda)  return(list(out,lambda1))
  return(out)
}

#' Parameter Estimation for the Box-Cox Transformation
#'
#' This function is copied from the geoR package which has been removed from the
#' CRAN repository.
#'
#' For more information see:
#' https://cran.r-project.org/web/packages/geoR/index.html
#' @param object a vector with the data
#' @param xmat a matrix with covariates values. Defaults to rep(1, length(y)).
#' @param lambda numerical value(s) for the transformation parameter lambda.
#'   Used as the initial value in the function for parameter estimation. If not
#'   provided default values are assumed. If multiple values are passed the one
#'   with highest likelihood is used as initial value.
#' @param lambda2 ogical or numerical value(s) of the additional transformation
#'   (see DETAILS below). Defaults to NULL. If TRUE this parameter is also
#'   estimated and the initial value is set to the absolute value of the minimum
#'   data. A numerical value is provided it is used as the initial value.
#'   Multiple values are allowed as for lambda.
#' @param add.to.data a constant value to be added to the data.
#' @importFrom stats optim optimize
geoR_boxcoxfit <- function (object, xmat, lambda, lambda2 = NULL, add.to.data = 0)
{
  call.fc <- match.call()
  data <- object + add.to.data
  if (is.null(lambda2) && any(data <= 0))
    stop("Transformation requires positive data")
  data <- as.vector(data)
  n <- length(data)
  if (missing(xmat))
    xmat <- rep(1, n)
  xmat <- as.matrix(xmat)
  if (any(xmat[, 1] != 1))
    xmat <- cbind(1, xmat)
  xmat <- xmat[!is.na(data), , drop = FALSE]
  data <- data[!is.na(data)]
  n <- length(data)
  beta.size <- ncol(xmat)
  if (nrow(xmat) != length(data))
    stop("xmat and data have incompatible lengths")
  lik.method <- "ML"
  if (all(data > 0))
    absmin <- 0
  else absmin <- abs(min(data)) + 1e-05 * diff(range(data))
  if (!is.null(lambda2)) {
    if (missing(lambda))
      lambda.ini <- seq(-2, 2, by = 0.2)
    else lambda.ini <- lambda
    lambda2.ini <- 0
    if (isTRUE(lambda2))
      lambda2.ini <- absmin
    else if (mode(lambda2) == "numeric")
      lambda2.ini <- lambda2
    lambdas.ini <- as.matrix(expand.grid(lambda.ini, lambda2.ini))
    if (length(as.matrix(lambdas.ini)) > 2) {
      lamlik <- apply(lambdas.ini, 1, .negloglik.boxcox,
                      data = data + absmin, xmat = xmat, lik.method = lik.method)
      lambdas.ini <- lambdas.ini[which(lamlik == min(lamlik)),
      ]
    }
    lambdas.ini <- unname(drop(lambdas.ini))
    lik.lambda <- stats::optim(par = lambdas.ini, fn = .negloglik.boxcox,
                               method = "L-BFGS-B", lower = c(-Inf, absmin), data = data,
                               xmat = xmat, lik.method = lik.method)
  }
  else {
    lik.lambda <- stats::optimize(.negloglik.boxcox, interval = c(-5,
                                                                  5), data = data, xmat = xmat, lik.method = lik.method)
    lik.lambda <- list(par = lik.lambda$minimum, value = lik.lambda$objective,
                       convergence = 0, message = "function optimize used")
  }
  lambda.fit <- lik.lambda$par
  if (length(lambda.fit) == 1)
    lambda.fit <- c(lambda.fit, 0)
  data <- data + lambda.fit[2]
  if (isTRUE(all.equal(unname(lambda.fit[1]), 0)))
    yt <- log(data)
  else yt <- ((data^lambda.fit[1]) - 1)/lambda.fit[1]
  beta <- solve(crossprod(xmat), crossprod(xmat, yt))
  mu <- drop(xmat %*% beta)
  sigmasq <- sum((yt - mu)^2)/n
  if (lik.method == "ML")
    loglik <- drop((-(n/2) * (log(2 * pi) + log(sigmasq) +
                                1)) + (lambda.fit[1] - 1) * sum(log(data)))
  temp <- 1 + lambda.fit[1] * mu
  fitted.y <- ((temp^((1/lambda.fit[1]) - 2)) * (temp^2 + ((1 -
                                                              lambda.fit[1])/2) * sigmasq))
  variance.y <- (temp^((2/lambda.fit[1]) - 2)) * sigmasq
  if (beta.size == 1) {
    fitted.y <- unique(fitted.y)
    variance.y <- unique(fitted.y)
  }
  beta <- drop(beta)
  if (length(beta) > 1)
    names(beta) <- paste("beta", 0:(beta.size - 1), sep = "")
  if (length(lik.lambda$par) == 1)
    lambda.fit <- lambda.fit[1]
  if (length(lik.lambda$par) == 2)
    names(lambda.fit) <- c("lambda", "lambda2")
  res <- list(lambda = lambda.fit, beta.normal = drop(beta),
              sigmasq.normal = sigmasq, loglik = loglik, optim.results = lik.lambda)
  res$call <- call.fc
  oldClass(res) <- "boxcoxfit"
  return(res)
}


#' fit crr model
#'
#' Wrapper function to fit fine and gray competing risk model using function crr
#' from package cmprsk
#'
#' @param f formula for the model. Currently the formula only works by using the
#'   name of the column in a dataframe. It does not work by using $ or []
#'   notation.
#' @param data dataframe containing data
#' @keywords model
#' @returns a competing risk model with the call appended to the list
#' @importFrom cmprsk crr
#' @seealso \code{\link{crr}}
#' @examples
#' # From the crr help file:
#' set.seed(10)
#' ftime <- rexp(200)
#' fstatus <- sample(0:2,200,replace=TRUE)
#' cov <- matrix(runif(600),nrow=200)
#' dimnames(cov)[[2]] <- c('x1','x2','x3')
#' df <- data.frame(ftime,fstatus,cov)
#' m1 <- crrRx(as.formula('ftime+fstatus~x1+x2+x3'),df)
#' # Nicely output to report:
#' rm_mvsum(m1,data=df,showN = TRUE,vif=TRUE)
#' @export
crrRx<-function(f,data){
  k<-as.character(f)[3]
  covs<-removedollar(k)
  ff<-modelmatrix(f,data)
  m1<-cmprsk::crr(ff[[1]][,1],ff[[1]][,2],ff[[2]])
  m1$call<-paste("~",covs)
  return(m1)
}




#' Get covariate summary dataframe
#'
#' Returns a dataframe corresponding to a descriptive table.
#'
#' Comparisons for categorical variables default to chi-square tests, but if
#' there are counts of <5 then the Fisher Exact test will be used and if this is
#' unsuccessful then a second attempt will be made computing p-values using MC
#' simulation. If testcont='ANOVA' then the t-test with unequal variance will be
#' used for two groups and an ANOVA will be used for three or more. The
#' statistical test used can be displayed by specifying show.tests=TRUE.
#'
#' The number of decimals places to display the statistics can be changed with
#' digits, but this will not change the display of p-values. If more significant
#' digits are required for p-values then use tableOnly=TRUE and format as
#' desired.
#'
#' @param data dataframe containing data
#' @param covs character vector with the names of columns to include in table
#' @param maincov covariate to stratify table by
#' @param digits number of digits for summarizing mean data, does not affect
#'   p-values
#' @param numobs named list overriding the number of people you expect to have
#'   the covariate
#' @param markup boolean indicating if you want latex markup
#' @param sanitize boolean indicating if you want to sanitize all strings to not
#'   break LaTeX
#' @param nicenames boolean indicating if you want to replace . and _ in strings
#'   with a space
#' @param IQR boolean indicating if you want to display the inter quantile range
#'   (Q1,Q3) as opposed to (min,max) in the summary for continuous variables
#' @param all.stats boolean indicating if all summary statistics (Q1,Q3 +
#'   min,max on a separate line) should be displayed. Overrides IQR.
#' @param pvalue boolean indicating if you want p-values included in the table
#' @param effSize boolean indicating if you want effect sizes included in the
#'   table. Can only be obtained if pvalue is also requested. Effect sizes are
#'   calculated with the rstatix package using Cramer V for categorical and Eta
#'   Squared for continuous covariates.
#' @param show.tests boolean indicating if the type of statistical used should
#'   be shown in a column beside the pvalues. Ignored if pvalue=FALSE.
#' @param dropLevels logical, indicating if empty factor levels be dropped from
#'   the output, default is TRUE.
#' @param excludeLevels a named list of covariate levels to exclude from
#'   statistical tests in the form list(varname =c('level1','level2')). These
#'   levels will be excluded from association tests, but not the table. This can
#'   be useful for levels where there is a logical skip (ie not missing, but not
#'   presented). Ignored if pvalue=FALSE.
#' @param full boolean indicating if you want the full sample included in the
#'   table, ignored if maincov is NULL
#' @param digits.cat number of digits for the proportions when summarizing
#'   categorical data (default: 0)
#' @param testcont test of choice for continuous variables,one of
#'   \emph{rank-sum} (default) or \emph{ANOVA}
#' @param testcat test of choice for categorical variables,one of
#'   \emph{Chi-squared} (default) or \emph{Fisher}
#' @param include_missing Option to include NA values of maincov. NAs will not
#'   be included in statistical tests
#' @param percentage choice of how percentages are presented ,one of
#'   \emph{column} (default) or \emph{row}
#' @keywords dataframe
#' @importFrom stats lm sd
#' @importFrom rstatix cramer_v eta_squared
#' @seealso \code{\link{fisher.test}},\code{\link{chisq.test}},
#'   \code{\link{wilcox.test}},\code{\link{kruskal.test}},and
#'   \code{\link{anova}}
covsum <- function (data, covs, maincov = NULL, digits = 1, numobs = NULL,
                    markup = TRUE, sanitize = TRUE, nicenames = TRUE, IQR = FALSE,
                    all.stats = FALSE, pvalue = TRUE, effSize = FALSE, show.tests = FALSE, dropLevels = TRUE,
                    excludeLevels = NULL, full = TRUE, digits.cat = 0, testcont = c("rank-sum test",
                                                                                    "ANOVA"), testcat = c("Chi-squared", "Fisher"), include_missing = FALSE,
                    percentage = c("column", "row"))
{
  if (missing(data))
    stop("data is a required argument")
  if (missing(covs))
    stop("covs is a required argument")
  if (!pvalue & effSize)
    stop("effSize can only be specified when pvalue = TRUE")
  if (!inherits(data, "data.frame"))
    stop("data must be supplied as a data frame.")
  if (!inherits(covs, "character"))
    stop("covs must be supplied as a character vector or string indicating variables in data")
  if (!is.null(maincov)) {
    if (!inherits(maincov, "character") | length(maincov) >
        1)
      stop("maincov must be supplied as a string indicating a variable in data")
  }
  missing_vars = setdiff(c(maincov, covs), names(data))
  if (length(missing_vars) > 0) {
    stop(paste("These covariates are not in the data:", paste0(missing_vars,
                                                               collapse = csep())))
  }
  for (v in c(maincov, covs)) {
    if (inherits(data[[v]], "logical"))
      data[[v]] <- factor(data[[v]])
    if (inherits(data[[v]], "character"))
      data[[v]] <- factor(data[[v]])
    if (inherits(data[[v]], c("Date", "POSIXt"))) {
      covs <- setdiff(covs, v)
      message(paste("Dates can not be summarised in this version of reportRmd.\n The variable",
                    v, "does not appear in the table."))
    }
  }
  if (!all(names(data[, c(maincov, covs)]) == names(data.frame(data[,
                                                                    c(maincov, covs)]))))
    warning("Non-standard variable names may cause problems. Check output.")
  missing_testcat <- missing(testcat)
  testcont <- match.arg(testcont)
  testcat <- match.arg(testcat)
  percentage <- match.arg(percentage)
  if (!pvalue) {
    show.tests <- FALSE
    excludeLevels <- NULL
  }
  if (!markup) {
    lbld <- identity
    addspace <- identity
    lpvalue <- identity
  }
  digits <- as.integer(digits)
  digits.cat <- as.integer(digits.cat)
  if (digits < 0)
    stop("parameter 'digits' cannot be negative!")
  if (digits.cat < 0)
    stop("parameter 'digits.cat' cannot be negative!")
  if (!sanitize)
    sanitizestr <- identity
  if (!nicenames)
    nicename <- identity
  if (dropLevels)
    data <- droplevels(data)
  if (!is.null(maincov)) {
    if (include_missing == FALSE)
      data <- data[!is.na(data[[maincov]]), ]
    levels <- names(table(data[[maincov]], useNA = "ifany"))
    levels <- c(list(levels), as.list(levels))
    if (length(na.omit(unique(data[[maincov]]))) == 1) {
      warning(paste("Only one value of the main covariate exists, show data for maincov =",
                    na.omit(unique(data[[maincov]]))))
      maincov <- NULL
      full = TRUE
      levels <- "NOMAINCOVNULLNA"
    }
  }
  else {
    full = TRUE
    levels <- "NOMAINCOVNULLNA"
  }
  N = nrow(data)
  if (!is.null(maincov)) {
    nmaincov <- c(sum(table(data[[maincov]], useNA = "ifany")),
                  table(data[[maincov]], useNA = "ifany"))
  }
  else {
    nmaincov <- N
    p <- NULL
  }
  out <- lapply(covs, function(cov) {
    ismiss = F
    n <- sum(table(data[[cov]]))
    if (!is.null(excludeLevels[[cov]])) {
      excludeLevel = excludeLevels[[cov]]
    }
    else excludeLevel = ""
    factornames <- NULL
    if (is.null(numobs[[cov]]))
      numobs[[cov]] <- nmaincov
    if (numobs[[cov]][1] - n > 0) {
      ismiss = T
      factornames <- c(factornames, "Missing")
    }
    if (is.factor(data[[cov]])) {
      factornames <- c(levels(data[[cov]]), factornames)
      if (!is.null(maincov)) {
        if (pvalue) {
          pdata = data[!(data[[cov]] %in% excludeLevel),
          ]
          lowcounts <- sum(table(pdata[[maincov]], pdata[[cov]],
                                 exclude = excludeLevel) < 5) > 0
          if (!missing_testcat & testcat == "Chi-squared" &
              lowcounts)
            warning(paste("Low counts are present in",
                          cov, "variable consider Fisher's test."),
                    call. = F)
          if ((missing_testcat & lowcounts) | testcat ==
              "Fisher") {
            p_type <- "Fisher Exact"
            p <- try(stats::fisher.test(pdata[[maincov]],
                                        pdata[[cov]])$p.value, silent = TRUE)
            if (effSize) {
              e_type <- "Cramer"
              e <- try(rstatix::cramer_v(table(pdata[[cov]], pdata[[maincov]])), silent = TRUE)
            }
            if (is.error(p)) {
              message("Using MC sim. Use set.seed() prior to function for reproducible results.")
              p <- try(stats::fisher.test(pdata[[maincov]],
                                          pdata[[cov]], simulate.p.value = T)$p.value,
                       silent = TRUE)
              p_type <- "MC sim"
              if (effSize) {
                e_type <- "Cramer"
                e <- try(rstatix::cramer_v(table(pdata[[cov]], pdata[[maincov]])), silent = TRUE)
              }
            }
          }
          else {
            p_type = "Chi Sq"
            p = try(stats::chisq.test(pdata[[maincov]],
                                      pdata[[cov]])$p.value, silent = TRUE)
            if (effSize) {
              e_type <- "Cramer"
              e <- try(rstatix::cramer_v(table(pdata[[cov]], pdata[[maincov]])), silent = TRUE)
            }
          }
          if (is.error(p))
            p <- NA
          p <- lpvalue(p)
          if (effSize) {
            if (is.error(e))
              e <- NA
            e <- lpvalue(e)
          }
        }
      }
      if (percentage == "column") {
        onetbl <- mapply(function(sublevel, N) {
          missing <- NULL
          if (is.na(sublevel[1]) | sublevel[1] != "NOMAINCOVNULLNA") {
            subdata <- subset(data, subset = data[[maincov]] %in%
                                sublevel)
          }
          else {
            subdata <- data
          }
          table <- table(subdata[[cov]])
          tbl <- table(subdata[[cov]])
          n <- sum(tbl)
          prop <- round(tbl/n, 2 + digits.cat) * 100
          prop <- sapply(prop, function(x) {
            if (!is.nan(x)) {
              x
            }
            else {
              0
            }
          })
          prop.fmt <- sprintf(paste0("%.", digits.cat,
                                     "f"), prop)
          tbl <- mapply(function(num, prop) {
            paste(num, " (", prop, ")", sep = "")
          }, tbl, prop.fmt)
          if (ismiss)
            missing <- N - n
          tbl <- c(tbl, lbld(missing))
          return(tbl)
        }, levels, numobs[[cov]])
      }
      if (percentage == "row") {
        onetbl <- mapply(function(sublevel, N) {
          missing <- NULL
          if (is.na(sublevel[1]) | sublevel[1] != "NOMAINCOVNULLNA") {
            subdata <- subset(data, subset = data[[maincov]] %in%
                                sublevel)
          }
          else {
            subdata <- data
          }
          table <- table(subdata[[cov]])
          tbl <- table(subdata[[cov]])
          n <- sum(tbl)
          if (ismiss)
            missing <- N - n
          tbl <- c(tbl, lbld(missing))
          return(tbl)
        }, levels, numobs[[cov]])
        if (ismiss) {
          if (dim(onetbl)[1] > 2) {
            onetbl[-nrow(onetbl), -1] <- t(apply(onetbl[-nrow(onetbl),
                                                        -1], 1, function(x) {
                                                          x <- as.numeric(x)
                                                          prop <- round(x/sum(x), 2 + digits.cat) *
                                                            100
                                                          prop.fmt <- sprintf(paste0("%.", digits.cat,
                                                                                     "f"), prop)
                                                          return(paste(x, " (", prop.fmt, ")", sep = ""))
                                                        }))
          }
          else {
            onetbl[-nrow(onetbl), -1] <- (function(x) {
              x <- as.numeric(x)
              prop <- round(x/sum(x), 2 + digits.cat) *
                100
              prop.fmt <- sprintf(paste0("%.", digits.cat,
                                         "f"), prop)
              return(paste(x, " (", prop.fmt, ")", sep = ""))
            })(onetbl[-nrow(onetbl), -1])
          }
        }
        else {
          if (!is.null(dim(onetbl))) {
            onetbl[, -1] <- t(apply(onetbl[, -1], 1,
                                    function(x) {
                                      x <- as.numeric(x)
                                      prop <- round(x/sum(x), 2 + digits.cat) *
                                        100
                                      prop.fmt <- sprintf(paste0("%.", digits.cat,
                                                                 "f"), prop)
                                      return(paste(x, " (", prop.fmt, ")",
                                                   sep = ""))
                                    }))
          }
          else {
            onetbl[-1] <- (function(x) {
              x <- as.numeric(x)
              prop <- round(x/sum(x), 2 + digits.cat) *
                100
              prop.fmt <- sprintf(paste0("%.", digits.cat,
                                         "f"), prop)
              return(paste(x, " (", prop.fmt, ")", sep = ""))
            })(onetbl[-1])
          }
        }
      }
    }
    else {
      if (all.stats) {
        factornames <- c("Mean (sd)", "Median (Q1,Q3)",
                         "Range (min, max)", factornames)
      }
      else factornames <- c("Mean (sd)", ifelse(IQR, "Median (Q1,Q3)",
                                                "Median (Min,Max)"), factornames)
      if (!is.null(maincov)) {
        if (pvalue) {
          if (testcont[1] == "rank-sum test") {
            if (length(unique(data[[maincov]])) == 2) {
              p_type = "Wilcoxon Rank Sum"
              p <- try(stats::wilcox.test(data[[cov]] ~
                                            data[[maincov]])$p.value, silent = T)
              if (effSize) {
                e_type <- "Eta sq"
                e <- try(rstatix::eta_squared(stats::aov(data[[cov]]~data[[maincov]], data=data)), silent = TRUE)
              }
            }
            else {
              p_type = "Kruskal Wallis"
              p <- try(stats::kruskal.test(data[[cov]] ~
                                             data[[maincov]])$p.value, silent = T)
              if (effSize) {
                e_type <- "Eta sq"
                e <- try(rstatix::eta_squared(stats::aov(data[[cov]]~data[[maincov]], data=data)), silent = TRUE)
              }
            }
          }
          else {
            if (length(unique(data[[maincov]])) == 2) {
              p_type = "t-test"
              p <- try(stats::t.test(data[[cov]] ~ data[[maincov]])$p.value,
                       silent = TRUE)
              if (effSize) {
                e_type <- "Eta sq"
                e <- try(rstatix::eta_squared(stats::aov(data[[cov]]~data[[maincov]], data=data)), silent = TRUE)
              }
            }
            else {
              p_type = "ANOVA"
              p <- try(stats::anova(stats::lm(data[[cov]] ~
                                                data[[maincov]]))[5][[1]][1], silent = TRUE)
              if (effSize) {
                e_type <- "Eta sq"
                e <- try(rstatix::eta_squared(stats::aov(data[[cov]]~data[[maincov]], data=data)), silent = TRUE)
              }
            }
          }
          if (is.error(p))
            p <- NA
          p <- lpvalue(p)
          if (effSize) {
            if (is.error(e))
              e <- NA
            e <- lpvalue(e)
          }
        }
      }
      onetbl <- mapply(function(sublevel, N) {
        missing <- NULL
        if (is.na(sublevel[1]) | sublevel[1] != "NOMAINCOVNULLNA") {
          subdata <- subset(data, subset = data[[maincov]] %in%
                              sublevel)
        }
        else {
          subdata <- data
        }
        if (ismiss) {
          n <- sum(table(subdata[[cov]]))
          missing <- N - n
        }
        sumCov <- round(summary(subdata[[cov]]), digits)
        if (sumCov[4] == "NaN") {
          meansd <- ""
          mmm <- ""
          if (all.stats)
            mmm = c("", "")
        }
        else {
          meansd <- paste(niceNum(sumCov["Mean"], digits),
                          " (", niceNum(sd(subdata[[cov]], na.rm = T),
                                        digits), ")", sep = "")
          mmm <- if (IQR | all.stats) {
            if (all(c(sumCov["Median"], sumCov["1st Qu."],
                      sumCov["3rd Qu."]) == floor(c(sumCov["Median"],
                                                    sumCov["1st Qu."], sumCov["3rd Qu."])))) {
              paste(sumCov["Median"], " (", sumCov["1st Qu."],
                    csep(), sumCov["3rd Qu."], ")", sep = "")
            }
            else {
              paste(niceNum(sumCov["Median"], digits),
                    " (", niceNum(sumCov["1st Qu."], digits),
                    csep(), niceNum(sumCov["3rd Qu."], digits),
                    ")", sep = "")
            }
          }
          else {
            if (all(c(sumCov["Median"], sumCov["Min."],
                      sumCov["Max."]) == floor(c(sumCov["Median"],
                                                 sumCov["Min."], sumCov["Max."])))) {
              paste(sumCov["Median"], " (", sumCov["Min."],
                    csep(), sumCov["Max."], ")", sep = "")
            }
            else {
              paste(niceNum(sumCov["Median"], digits),
                    " (", niceNum(sumCov["Min."], digits),
                    csep(), niceNum(sumCov["Max."], digits),
                    ")", sep = "")
            }
          }
          if (all.stats) {
            mmm <- c(mmm, if (all(c(sumCov["Min."], sumCov["Max."]) ==
                                  floor(c(sumCov["Min."], sumCov["Max."])))) {
              paste("(", sumCov["Min."], csep(), sumCov["Max."],
                    ")", sep = "")
            } else {
              paste("(", niceNum(sumCov["Min."], digits),
                    csep(), niceNum(sumCov["Max."], digits),
                    ")", sep = "")
            })
          }
        }
        tbl <- c(meansd, mmm, lbld(missing))
        return(tbl)
      }, levels, numobs[[cov]])
    }
    factornames <- addspace(sanitizestr(nicename(factornames)))
    if (is.null(nrow(onetbl))) {
      onetbl <- matrix(data = onetbl, ncol = length(onetbl),
                       nrow = 1)
    }
    onetbl <- cbind(factornames, onetbl)
    if (!is.null(maincov)) {
      onetbl <- rbind(c(lbld(sanitizestr(nicename(cov))),
                        rep("", length(levels[[1]]) + 1)), onetbl)
      if (pvalue) {
        p_NA = rep("", nrow(onetbl) - 1)
        p_NA[levels(data[[cov]]) %in% excludeLevel] <- "excl"
        onetbl <- cbind(onetbl, c(p, p_NA))
      }
      if (effSize) {
        e_NA = rep("", nrow(onetbl) - 1)
        e_NA[levels(data[[cov]]) %in% excludeLevel] <- "excl"
        onetbl <- cbind(onetbl, c(e, e_NA))
      }
      if (show.tests & effSize) {
        onetbl <- cbind(onetbl, c(paste(p_type, ", ", e_type, sep=""), rep("", nrow(onetbl) -
                                                                             1)))
      }
      if (show.tests & !effSize) {
        onetbl <- cbind(onetbl, c(paste(p_type), rep("", nrow(onetbl) -
                                                       1)))
      }
    }
    else {
      onetbl <- rbind(c(lbld(sanitizestr(nicename(cov))),
                        ""), onetbl)
    }
    rownames(onetbl) <- NULL
    colnames(onetbl) <- NULL
    return(onetbl)
  })
  table <- do.call("rbind", lapply(out, data.frame, stringsAsFactors = FALSE))
  table <- data.frame(apply(table, 2, unlist), stringsAsFactors = FALSE)
  rownames(table) <- NULL
  if (!is.null(maincov)) {
    colnm_table <- c("Covariate", paste("Full Sample (n=",
                                        N, ")", sep = ""), mapply(function(x, y) {
                                          paste(x, " (n=", y, ")", sep = "")
                                        }, names(table(data[[maincov]], useNA = "ifany")), table(data[[maincov]],
                                                                                                 useNA = "ifany")))
    if (pvalue)
      colnm_table <- c(colnm_table, "p-value")
    if (effSize)
      colnm_table <- c(colnm_table, "Effect Size")
    if (show.tests)
      colnm_table <- c(colnm_table, "StatTest")
    colnames(table) <- colnm_table
  }
  else {
    colnames(table) <- c("Covariate", paste("n=", N, sep = ""))
  }
  colnames(table) <- sanitizestr(colnames(table))
  if (!full)
    table <- table[, -2]
  return(table)
}

#' Get univariate summary dataframe
#'
#' Returns a dataframe corresponding to a univariate regression table
#'
#' Univariate summaries for a number of covariates, the type of model can be
#' specified. If unspecified the function will guess the appropriate model based
#' on the response variable.
#'
#' Confidence intervals are extracted using confint where possible. Otherwise
#' Student t distribution is used for linear models and the Normal distribution
#' is used for proportions.
#'
#' returnModels can be used to return a list of the univariate models, which
#' will be the same length as covs. The data used to run each model will include
#' all cases with observations on the response and covariate. For gee models the
#' data are re-ordered so that the ids appear sequentially and proper estimates
#' are given.
#' @param response string vector with name of response
#' @param covs character vector with the names of columns to fit univariate
#'   models to
#' @param data dataframe containing data
#' @param digits number of digits to round to
#' @param id character vector which identifies clusters. Used for GEE and coxph
#'   models.
#' @param corstr character string specifying the correlation structure. Only
#'   used for geeglm. The following are permitted: '"independence"',
#'   '"exchangeable"', '"ar1"', '"unstructured"' and '"userdefined"'
#' @param family specify details of the model used. This argument does not need
#'   to be specified and should be used with caution. By default, gaussian
#'   errors are used for linear models, the binomial family with logit link is
#'   used for logistic regression and poisson with log link is used for poisson
#'   regression. This can be specified with the type argument, or will be
#'   inferred from the data type. See \code{\link{family}}. Ignored for ordinal
#'   and survival regression and if the type argument is not explicitly
#'   specified.
#' @param type string indicating he type of univariate model to fit. The
#'   function will try and guess what type you want based on your response. If
#'   you want to override this you can manually specify the type. Options
#'   include "linear", "logistic", "poisson", coxph", "crr", "boxcox" and
#'   "ordinal"
#' @param  gee boolean indicating if gee models should be fit to account for
#'   correlated observations. If TRUE then the id argument must specify the
#'   column in the data which indicates the correlated clusters.
#' @param strata character vector of covariates to stratify by. Only used for
#'   coxph and crr
#' @param markup boolean indicating if you want latex markup
#' @param sanitize boolean indicating if you want to sanitize all strings to not
#'   break LaTeX
#' @param nicenames boolean indicating if you want to replace . and _ in strings
#'   with a space
#' @param showN boolean indicating if you want to show sample sizes
#' @param CIwidth width of confidence interval, default is 0.95
#' @param reflevel manual specification of the reference level. Only used for
#'   ordinal. This may allow you to debug if the function throws an error.
#' @param returnModels boolean indicating if a list of fitted models should be
#'   returned.
#' @seealso
#' \code{\link{lm}},\code{\link{glm}},\code{\link{crr}},\code{\link{coxph}},
#' \code{\link{lme}},\code{\link{geeglm}},\code{\link{polr}}
#' @keywords dataframe
#' @importFrom MASS polr
#' @importFrom survival coxph Surv
#' @importFrom aod wald.test
#' @importFrom geepack geeglm
#' @importFrom stats na.omit as.formula anova glm lm qnorm qt confint
uvsum <- function (response, covs, data, digits=2,id = NULL, corstr = NULL, family = NULL,
                   type = NULL, gee=FALSE,strata = 1, markup = TRUE, sanitize = TRUE, nicenames = TRUE,
                   showN = TRUE, CIwidth = 0.95, reflevel=NULL,returnModels=FALSE)
{

  if (!markup) {
    lbld <- identity
    addspace <- identity
    lpvalue <- identity
  }
  if (!sanitize)  sanitizestr <- identity
  if (!nicenames) nicename <- identity
  if (inherits(data[[response[1]]],"character")) data[[response[1]]] <- factor(data[[response[1]]])
  if (!inherits(strata,"numeric")) {
    strataVar = strata
    strata <- sapply(strata, function(stra) {
      paste("strata(", stra, ")", sep = "")
    })
  }
  else {
    strataVar <- ""
    strata <- ""
  }
  if (!is.null(type)) {
    if (length(response)==1 & (type %in% c('coxph','crr')))
      stop('Please specify two variables in the response for survival models. \nExample: response=c("time","status")')
    if (length(response)==2 & !(type %in% c('coxph','crr')))
      stop('Response can only be of length one for non-survival models.')
    if (type == "logistic") {
      beta <- "OR"
      if (is.null(family)) family='binomial'
    }
    else if (type == "poisson") {
      if (all(data[[response]]==as.integer(data[[response]]))){
        data[[response]]=as.integer(data[[response]])
      }
      else {
        stop('Poisson regression requires an integer response.')
      }
      beta <- "RR"
      if (is.null(family)) family='poisson'
    }
    else if (type == "linear" | type == "boxcox") {
      beta <- "Estimate"
      if (is.null(family)) family='gaussian'
    }
    else if (type == "coxph" | type == "crr") {
      beta <- "HR"
    }
    else if (type == "ordinal") {
      if (!inherits(data[[response[1]]],c("factor","ordered"))) {
        warning("Response variable is not a factor, will be converted to an ordered factor")
        data[[response]] <- factor(data[[response]],
                                   ordered = T)
      }
      if (!is.null(reflevel)) {
        data[[response]] <- stats::relevel(data[[response]],
                                           ref = reflevel)
      }
      beta <- "OR"
    }
    else {
      stop("type must be either coxph, logistic, linear, poisson, boxcox, crr, ordinal (or NULL)")
    }
  }
  else {
    if (length(response) == 2) {
      # Check that responses are numeric
      for (i in 1:2) if (!is.numeric(data[[response[i]]])) stop('Both response variables must be numeric')
      if (length(unique(data[[response[2]]])) < 3) {
        type <- "coxph"
      }
      else {
        type <- "crr"
      }
      beta <- "HR"
    }
    else if (length(na.omit(unique(data[[response]]))) == 2) {
      type <- "logistic"
      beta <- "OR"
      family="binomial"
    }
    else if (inherits(data[[response[1]]],"ordered")) {
      type <- "ordinal"
      beta <- "OR"
      if (!is.null(reflevel)) {
        data[[response]] <- stats::relevel(data[[response]],
                                           ref = reflevel)
      }
    }
    else if (inherits(data[[response[1]]],"integer")) {
      type <- "poisson"
      beta <- "RR"
      family="poisson"
    }
    else {
      if (!inherits(data[[response[1]]],"numeric")) stop('Response variable must be numeric')
      type <- "linear"
      beta <- "Estimate"
      family='gaussian'
    }
  }
  beta = betaWithCI(beta, CIwidth)
  if (strata != "" & type != "coxph") {
    stop("strata can only be used with coxph")
  }
  if (!is.null(id)){
    if (! (gee | type =='coxph')) {
      warning('id argument will be ignored. This is used only for survival strata or clustering in GEE. To run a GEE model set gee=TRUE.')
    }
  }
  if (!is.null(corstr)){
    if (! (gee | type =='coxph')) {
      warning('id argument will be ignored. This is used only for survival strata or clustering in GEE. To run a GEE model set gee=TRUE.')
    }
  }
  if (gee){
    if (!type %in% c('linear','logistic')) stop('GEE models currently only implemented for logistic or linear regression.')
    if (is.null(id)) stop('The id argument must be set for gee models to indicate clusters.')
    if (is.null(corstr)) stop ('You must provide correlation structure (i.e. corstr="independence") for GEE models.')
  }
  if (returnModels) modelList <- NULL
  out <- lapply(covs, function(x_var) {
    data <- data[,intersect(c(response, x_var, strataVar,id),names(data))]
    data <- stats::na.omit(data)
    m2 <- NULL
    if (gee){
      data <- data[order(data[[id]]),]
      idf <- as.numeric(as.factor(data[[id]]))
      data$idf <- idf
    }
    if (inherits(data[[x_var]],c("ordered", "factor"))) {
      data[[x_var]] = droplevels(data[[x_var]])
    }
    if (is.factor(data[[x_var]])) {
      x_var_str <- x_var
      levelnames = sapply(sapply(sapply(levels(data[[x_var]]),
                                        nicename), sanitizestr), addspace)
      x_var_str <- lbld(sanitizestr(nicename(x_var)))
      title <- NULL
      body <- NULL
    } else x_var_str <- lbld(sanitizestr(nicename(x_var)))

    if (type == "coxph") {
      f <- paste(paste("survival::Surv(",
                       response[1], ",", response[2], ")",
                       sep = ""), "~", x_var, ifelse(strata ==
                                                       "", "", "+"), paste(strata,
                                                                           collapse = "+"), sep = "")
      if (is.null(id)) {
        eval(parse(text = paste('m2 <- survival::coxph(formula=as.formula(',f,'), data = data)')))
      } else{
        eval(parse(text = paste('m2 <- survival::coxph(formula=as.formula(',f,'),id =',id,', data = data)')))
      }
      m <- summary(m2,conf.int = CIwidth)
      hr <- m$conf.int[, c(1, 3, 4)]
      pvals <- m$coefficients[,"Pr(>|z|)"]
      globalpvalue <- m$logtest['pvalue']

    }
    else if (type == "crr") {
      eval(parse(text = paste('m2 <- crrRx(',paste(paste(response,collapse = "+"),
                                                   "~", x_var, sep = ""),
                              ',data = data)')))
      m <- summary(m2,conf.int = CIwidth)
      hr <- m$conf.int[, c(1, 3, 4)]
      pvals <- m$coef[,5]
      globalpvalue <- try(aod::wald.test(b = m2$coef,
                                         Sigma = m2$var, Terms = seq_len(length(m2$coef)))$result$chi2[3])

    }
    else if (type %in% c("logistic","poisson")) {
      if (gee){
        eval(parse(text = paste0("m2 <- geepack::geeglm(",paste(response, "~",x_var, sep = ""),
                                 ",family = ",family,",",
                                 "data = data, id = idf, corstr = '",corstr,"')")))
        globalpvalue <- try(aod::wald.test(b = m2$coefficients[-1],
                                           Sigma = (m2$geese$vbeta)[-1, -1], Terms = seq_len(length(m2$coefficients[-1])))$result$chi2[3])
        m <- summary(m2)$coefficients
        Zmult = stats::qnorm(1 - (1 - CIwidth)/2)
        hr <- cbind(exp(m[,1]),exp(m[, 1] - Zmult * m[, 2]),
                    exp(m[,1] + Zmult * m[, 2]))
        pvals <- m[-1,"Pr(>|W|)"]
      }
      else{
        eval(parse(text = paste("m2 <- glm(",paste(response, "~",x_var, sep = ""),
                                ",family = ",family,",",
                                "data = data)")))
        m2_null <- stats::update(m2,formula=as.formula(paste0(response,'~1')),data=m2$model)
        globalpvalue <- try(as.vector(stats::na.omit(anova(m2_null,m2,test="LRT")[,"Pr(>Chi)"]))) # LRT
        m <- summary(m2)$coefficients
        hr <- cbind(exp(m[,1]),exp(confint(m2,level=CIwidth)[,]))
        pvals <- m[-1,"Pr(>|z|)"]

      }
      hr <- hr[-1,]
    }
    else if (type %in% c("linear", "boxcox")) {
      if (gee){
        eval(parse(text = paste0("m2 <- geepack::geeglm(",
                                 paste(response, "~",x_var, sep = ""),
                                 ",data = data, id = idf, corstr = '",corstr,
                                 "', family = ",family,")")))
        m <- summary(m2)$coefficients
        globalpvalue <- try(aod::wald.test(b = m2$coefficients[-1],
                                           Sigma = vcov(m2)[-1, -1], Terms = seq_len(length(m2$coefficients[-1])))$result$chi2[3],silent = T)
        Tmult = stats::qt(1 - (1 - CIwidth)/2, m2$df.residual)
        hr <- cbind(m[,1], m[, 1] - Tmult * m[, 2],
                    m[, 1] +Tmult * m[, 2])
        pvals <- m[-1,"Pr(>|W|)"]
      } else {
        if (type =='linear') {
          eval(parse(text = paste('m2 <- lm(',
                                  paste(response, "~",x_var, sep = ""),
                                  ',data = data)')))
        }
        else {
          eval(parse(text = paste('m2 <- boxcoxfitRx(',
                                  paste(response,"~", x_var, sep = ""),
                                  ',data = data)')))

        }
        m2_null <- lm(formula=as.formula(paste0(response,'~1')),data=m2$model)
        globalpvalue <- try(as.vector(stats::na.omit(anova(m2_null,m2,test="LRT")[,"Pr(>Chi)"])),silent = T) # LRT
        m <- summary(m2)$coefficients
        hr <- cbind(m[,1], confint(m2,level=CIwidth))
        pvals <- m[-1,4]
      }
      hr <- hr[-1,]
    }
    else if (type == "ordinal") {
      eval(parse(text = paste('m2 = MASS::polr(data = data,',
                              paste(response,"~", x_var, sep = ""),
                              ',method = "logistic",Hess = TRUE)')))
      m <- data.frame(summary(m2)$coef)
      m <- m[grep(x_var, rownames(summary(m2)$coef)),]

      m2_null <- stats::update(m2,data=m2$model,formula=as.formula(paste0(response,'~1' )))
      globalpvalue <- try(as.vector(stats::na.omit(anova(m2_null,m2)[,"Pr(Chi)"])))

      pvals <- stats::pt(abs(m[,3]),m2$df.residual, lower.tail = FALSE)*2
      if (length(pvals)>1){
        hr <- cbind(exp(m[,1]),exp(confint(m2,level=CIwidth)))
      } else {
        hr <- c(exp(m[,1]),exp(confint(m2,level=CIwidth)))
      }
    }
    hrmat <- matrix(hr,ncol = 3)
    if (is.error(globalpvalue))  globalpvalue <- "NA"
    if (is.factor(data[[x_var]])){
      hazardratio <- c("Reference", apply(hrmat, 1, psthr,digits))

      if (length(pvals)>1){
        pvalue <- c("", sapply(pvals,lpvalue))
        title <- c(x_var_str, "", "", lpvalue(globalpvalue))

      } else {
        pvalue <- sapply(pvals,lpvalue)
        title <- c(x_var_str, "", pvalue, lpvalue(globalpvalue))
      }
      if (length(levelnames) == 2) {
        body <- cbind(levelnames, hazardratio, c("",
                                                 ""), c("", ""))
      }
      else {
        body <- cbind(levelnames, hazardratio, pvalue,
                      rep("", length(levelnames)))
      }
      out <- rbind(title, body)
    } else {
      out <- matrix(c(x_var_str,
                      psthr(hr,digits),
                      lpvalue(pvals),
                      lpvalue(globalpvalue)),
                    ncol = 4)

    }
    if (showN) {
      n_by_level = nrow(data)
      if (is.factor(data[[x_var]])) {
        n_by_level = c(n_by_level, as.vector(table(data[[x_var]])))
      }
      out <- cbind(out, n_by_level)
    }
    if (returnModels) {
      m2$data <- data
      modelList[[x_var]] <<- m2
    }
    rownames(out) <- NULL
    colnames(out) <- NULL
    return(list(out, nrow(out)))
  })

  table <- lapply(out, function(x) {
    return(x[[1]])
  })
  table <- do.call("rbind", lapply(table, data.frame,
                                   stringsAsFactors = FALSE))
  colName <- c("Covariate", sanitizestr(beta),
               "p-value", "Global p-value")
  if (showN) colName <- c(colName,"N")
  colnames(table) <- colName
  table[,"Global p-value"] <- ifelse(table[,'p-value']=='',table[,"Global p-value"],'')
  if (all(table[,"Global p-value"]=='')) table <- table[, -which(colnames(table)=="Global p-value")]
  colnames(table) <- sapply(colnames(table), lbld)
  if (returnModels) return(list(table,models=modelList)) else return(table)
}

#' Get multivariate summary dataframe
#'
#' Returns a dataframe with the model summary and global p-value for multi-level
#' variables.
#'
#' Global p-values are likelihood ratio tests for lm, glm and polr models. For
#' lme models an attempt is made to re-fit the model using ML and if,successful
#' LRT is used to obtain a global p-value. For coxph models the model is re-run
#' without robust variances with and without each variable and a LRT is
#' presented. If unsuccessful a Wald p-value is returned. For GEE and CRR models
#' Wald global p-values are returned.
#'
#' If the variance inflation factor is requested (VIF=T) then a generalised VIF
#' will be calculated in the same manner as the car package.
#'
#' VIF for competing risk models is computed by fitting a linear model with a
#' dependent variable comprised of the sum of the model independent variables
#' and then calculating VIF from this linear model.
#'
#' @param model fitted model object
#' @param data dataframe containing data
#' @param digits number of digits to round to
#' @param showN boolean indicating sample sizes should be shown for each
#'   comparison, can be useful for interactions
#' @param markup boolean indicating if you want latex markup
#' @param sanitize boolean indicating if you want to sanitize all strings to not
#'   break LaTeX
#' @param nicenames boolean indicating if you want to replace . and _ in strings
#'   with a space.
#' @param CIwidth width for confidence intervals, defaults to 0.95
#' @param vif boolean indicating if the variance inflation factor should be
#'   included. See details
#' @keywords dataframe
#' @importFrom stats na.omit formula model.frame anova qnorm vcov
#' @references John Fox & Georges Monette (1992) Generalized Collinearity
#'   Diagnostics, Journal of the American Statistical Association, 87:417,
#'   178-183, DOI: 10.1080/01621459.1992.10475190
#' @references  John Fox and Sanford Weisberg (2019). An {R} Companion to
#'   Applied Regression, Third Edition. Thousand Oaks CA: Sage. URL:
#'   https://socialsciences.mcmaster.ca/jfox/Books/Companion
mvsum <- function (model, data, digits=2, showN = FALSE, markup = TRUE, sanitize = TRUE, nicenames = TRUE,
                   CIwidth = 0.95, vif=TRUE)
{
  if (!markup) {
    lbld <- identity
    addspace <- identity
    lpvalue <- identity
  }
  if (!sanitize)
    sanitizestr <- identity
  if (!nicenames)
    nicename <- identity
  if (inherits(model,c("lm", "lme", "multinom",
                       "survreg", "polr"))) {
    call <- Reduce(paste,
                   deparse(stats::formula(model$terms),
                           width.cutoff = 500))
  }  else if (inherits(model,c("crr"))) {
    call <- paste(deparse(model$call), collapse = "")
  }  else call <- paste(deparse(model$formula), collapse = "")
  call <- unlist(strsplit(call, "~", fixed = T))[2]
  call <- unlist(strsplit(call, ",", fixed = T))[1]
  if (substr(call, nchar(call), nchar(call)) == "\"")
    call <- substr(call, 1, nchar(call) - 1)
  call <- unlist(strsplit(call, "\"", fixed = T))[1]
  call <- unlist(strsplit(call, "+", fixed = T))
  call <- unlist(strsplit(call, "*", fixed = T))
  call <- unlist(strsplit(call, ":", fixed = T))
  call <- unique(call)
  call <- call[which(is.na(sapply(call, function(cov) {
    charmatch("strata(", cov)
  })) == T)]
  call <- gsub("\\s", "", call)
  type <- class(model)[1]
  if (type == "lm") {
    betanames <- attributes(summary(model)$coef)$dimnames[[1]][-1]
    beta <- "Estimate"
    expnt = FALSE
    ss_data <- model$model
  }
  else if (type == "polr") {
    expnt = TRUE
    betanames <- names(model$coefficients)
    beta <- "OR"
    ss_data <- model$model
  }
  else if (type == "lme") {
    expnt = FALSE
    betanames <- names(model$coef$fixed)[-1]
    beta <- "Estimate"
    ss_data <- model$data
  }
  else if (type == "glm") {
    if (model$family$link == "logit") {
      beta <- "OR"
      expnt = TRUE
    } else if (model$family$link == "log") {
      beta <- "RR"
      expnt = TRUE
    } else {
      beta <- "Estimate"
      expnt = FALSE
    }
    betanames <- names(model$coef)[-1]
    ss_data <- model$model
  }
  else if (type == "geeglm") {
    if (model$family$link == "logit") {
      beta <- "OR"
      expnt = TRUE
    } else if (model$family$link == "log") {
      beta <- "RR"
      expnt = TRUE
    } else {
      beta <- "Estimate"
      expnt = FALSE
    }
    betanames <- attributes(summary(model)$coef)$row.names[-1]
    ss_data <- model$model
  }
  else if (type == "coxph" | type == "crr") {
    beta <- "HR"
    expnt = TRUE
    # betanames <- attributes(model$terms)$term.labels
    # if (is.null(betanames))    betanames <- attributes(summary(model)$coef)$dimnames[[1]]
    betanames <- attributes(summary(model)$coef)$dimnames[[1]]
    ss_data <- try(stats::model.frame(model$call$formula, eval(parse(text = paste("data=",
                                                                                  deparse(model$call$data))))), silent = TRUE)
  }
  else {
    stop("type must be either polr, coxph, glm, lm, geeglm, crr, lme (or NULL)")
  }
  if (inherits(ss_data,"data.frame")) {
    if ('(weights)' %in% names(ss_data))
      names(ss_data)<- gsub('[(]weights[)]',as.character(model$call[['weights']]),names(ss_data))
    data <- ss_data
  } else if (type=='crr'){
    if (missing(data)){
      stop("Data can not be derived from model, data argument must be supplied.")
    } else if (model$n!=nrow(data)) {
      if (showN) stop('For crr models, the supplied data frame can contain only non-missing data.\n Either set showN = FALSE or run na.omit() on a data frame containing only model variables.')
    }
  } else if (type=='coxph'){
    if (missing(data)) stop("Data can not be derived from model, data argument must be supplied.")
    data <- na.omit(data[,c(dimnames(model$y)[[2]],betanames)])
  } else {
    stop("Data can not be derived from model, check model object.")
  }
  beta = betaWithCI(beta, CIwidth)
  ucall = unique(call)
  if (length(setdiff(ucall,names(data)))>0) stop('Currently this function is only implemented to work with standard variable names.\n Try converting the data to a standard data.frame with data.frame(data) and re-running the model to use rm_mvsum.')
  indx = try(matchcovariate(betanames, ucall),silent = T)
  if (is.error(indx)) stop('This function not yet implemented for complex function calls. Try re-specifying the model.')
  for (v in ucall) {
    if (inherits(data[[v]], "character"))
      data[[v]] <- factor(data[[v]])
  }
  if (min(indx) == -1)
    stop("Factor name + level name is the same as another factor name. Please change. Will fix this issue later")
  y <- betaindx(indx)
  if (type %in% c("lm", "glm", "geeglm", "lme")) {
    y <- lapply(y, function(x) {
      x + 1
    })
    betanames <- c("intercept", betanames)
  }
  out <- lapply(y, function(covariateindex) {
    betaname <- betanames[covariateindex]
    betaname <- strsplit(betaname, ":", fixed = T)
    oldcovname <- covnm(betaname[[1]], call)
    oldcovname <- getvarname(oldcovname)
    oldcovname <- paste(oldcovname,collapse = ":")
    levelnameslist <- lapply(betaname, function(level) {
      mapply(function(lvl, cn) {
        result <- ifelse(length(grep(paste0(cn, cn),
                                     lvl)) > 0, unlist(sub(paste0(cn, cn), cn, lvl)),
                         unlist(sub(cn, "", lvl)))
        out <- ifelse(result == "", cn, result)
      }, level, oldcovname)
    })
    levelnames <- unlist(lapply(levelnameslist, function(x) paste(x,
                                                                  collapse = ":")))
    covariatename <- oldcovname
    reference = NULL
    title = NULL
    body = NULL
    if (type == "lme") {
      globalpvalue <- NA
      f <- paste0('. ~ . -',oldcovname)
      if ( length(f)==1){
        m_small <- try(stats::update(model,paste0('. ~ . -',oldcovname),data=data,method='ML'),silent=TRUE)
        if (!is.error(m_small)){
          m_new <- stats::update(model,method='ML')
          globalpvalue <- try(as.vector(stats::na.omit(anova(m_small,m_new)[,"p-value"])),silent=T) # LRT
        }
      }
      if (is.na(globalpvalue)| is.error(globalpvalue)) {
        globalpvalue <- try(aod::wald.test(b = model$coef$fixed[covariateindex],
                                           Sigma = vcov(model)[covariateindex, covariateindex],
                                           Terms = seq_along(covariateindex))$result$chi2[3],silent = T)
      }
    } else if (type=='glm'){
      m_small <- try(stats::update(model,paste0('. ~ . -',oldcovname),data=data),silent = T)
      globalpvalue <- try(as.vector(stats::na.omit(anova(m_small,model,test='LRT')[,"Pr(>Chi)"])),silent = T)
    } else if (type == "polr") {
      m_small <- try(stats::update(model,paste0('. ~ . -',oldcovname),data=data),silent=TRUE)
      globalpvalue <- try(as.vector(stats::na.omit(anova(m_small,model)[,"Pr(Chi)"])),silent=TRUE)
    } else if (type == "crr" ) { # Leave as Wald Test
      globalpvalue <- try(aod::wald.test(b = model$coef[covariateindex],
                                         Sigma = model$var[covariateindex, covariateindex],
                                         Terms = seq_along(covariateindex))$result$chi2[3],
                          silent = T)
    } else if (type=='geeglm'){ # Leave as Wald Test
      globalpvalue <- try(aod::wald.test(b = model$coefficients[covariateindex],
                                         Sigma = (model$geese$vbeta)[covariateindex, covariateindex],
                                         Terms = seq_len(length(model$coefficients[covariateindex])))$result$chi2[3],
                          silent = T)

    } else if (type=='coxph') {
      m_data <- data
      names(m_data)[1] <- 'y'
      # R 4.2.2 breaks this
      # if (ncol(m_data)>2){
      #   # run models without robust variances
      #   m_full <- try(survival::coxph(y~.,data = m_data,robust=FALSE),silent=TRUE)
      #   m_small <- try(survival::coxph(y~.,data = m_data[,-which(names(m_data)==oldcovname)],robust=FALSE),silent=TRUE)
      #   globalpvalue <- try(as.vector(stats::na.omit(anova(m_small,m_full)[,"Pr(>|Chi|)"])),silent = T)
      # } else {
      #   globalpvalue <- try(as.vector(stats::na.omit(anova(model)[,"Pr(>|Chi|)"])),silent = T)
      # }
      # New code 19 Dec 2022
      # run models without robust variances
      m_full <- try(survival::coxph(y~.,data = m_data,robust=FALSE),silent=TRUE)
      m_small <- try(survival::coxph(y~.,data = m_data[,-which(names(m_data)==oldcovname)],robust=FALSE),silent=TRUE)
      gp_aov <- try(anova(m_small,m_full),silent = T)
      if (inherits(gp_aov,'try-error')) globalpvalue <- gp_aov else globalpvalue <- as.vector(stats::na.omit(gp_aov[,4]))

    } else {
      m_small <- try(stats::update(model,paste0('. ~ . -',oldcovname),data=data),silent=TRUE)
      globalpvalue <- try(as.vector(stats::na.omit(anova(m_small,model)[,"Pr(>F)"])),silent = T)
    }
    if (is.error(globalpvalue)) globalpvalue <- "NA"
    if (!identical(lpvalue,identity)) globalpvalue <- lpvalue(globalpvalue,digits)
    if (type == "coxph" | type == "crr") {
      hazardratio <- c(apply(matrix(summary(model, conf.int = CIwidth)$conf.int[covariateindex,
                                                                                c(1, 3, 4)], ncol = 3), 1, psthr,digits))
      pvalues <- c(sapply(summary(model)$coef[covariateindex,
                                              5], lpvalue))
    }
    else if (type == "glm" & expnt) {
      m <- summary(model, conf.int = CIwidth)$coefficients
      Z_mult = qnorm(1 - (1 - CIwidth)/2)
      hazardratio <- apply(cbind(exp(m[covariateindex, 1]),
                                 exp(m[covariateindex, 1] - Z_mult * m[covariateindex, 2]),
                                 exp(m[covariateindex, 1] + Z_mult * m[covariateindex, 2])), 1, psthr,digits)
      pvalues <- c(sapply(m[covariateindex, 4], lpvalue))
    }
    else if (type == "geeglm" & expnt) {
      m <- summary(model, conf.int = CIwidth)$coefficients
      Z_mult = qnorm(1 - (1 - CIwidth)/2)
      hazardratio <- apply(cbind(exp(m[covariateindex, 1]),
                                 exp(m[covariateindex, 1] - Z_mult * m[covariateindex,2]),
                                 exp(m[covariateindex, 1] + Z_mult * m[covariateindex, 2])), 1, psthr,digits)
      pvalues <- c(sapply(m[covariateindex, 4], lpvalue))
    }
    else if (type == "polr") {
      m <- summary(model)$coefficients
      Z_mult = qnorm(1 - (1 - CIwidth)/2)
      hazardratio <- apply(cbind(exp(m[covariateindex,1]),
                                 exp(m[covariateindex, 1] - Z_mult * m[covariateindex, 2]),
                                 exp(m[covariateindex, 1] + Z_mult * m[covariateindex, 2])), 1, psthr,digits)
      pvalues = stats::pnorm(abs(m[covariateindex, "Value"]/m[covariateindex,
                                                              "Std. Error"]), lower.tail = FALSE) * 2
      pvalues <- c(sapply(pvalues, lpvalue))
    }
    else if (type == "lm" | type == "glm" & !expnt) {
      T_mult = abs(stats::qt((1 - CIwidth)/2, model$df.residual))
      m <- summary(model, conf.int = CIwidth)$coefficients
      hazardratio <- apply(cbind(m[covariateindex, "Estimate"],
                                 m[covariateindex, "Estimate"] - T_mult * m[covariateindex, "Std. Error"],
                                 m[covariateindex, "Estimate"] + T_mult * m[covariateindex, "Std. Error"]), 1, psthr,digits)
      pvalues <- sapply(m[covariateindex, 4], lpvalue)
    }
    else if (type == "geeglm" & !expnt) {
      T_mult = abs(stats::qt((1 - CIwidth)/2, model$df.residual))
      m <- summary(model, conf.int = CIwidth)$coefficients
      hazardratio <- apply(cbind(m[covariateindex, "Estimate"],
                                 m[covariateindex, "Estimate"] - T_mult * m[covariateindex, "Std.err"],
                                 m[covariateindex, "Estimate"] + T_mult * m[covariateindex, "Std.err"]), 1, psthr,digits)
      pvalues <- sapply(m[covariateindex, 4], lpvalue)
    }
    else if (type == "lme") {
      T_mult = abs(stats::qt((1 - CIwidth)/2, summary(model)$fixDF$X))[covariateindex]
      m <- summary(model, conf.int = CIwidth)$tTable
      hazardratio <- apply(cbind(m[covariateindex, 1],
                                 m[covariateindex, 1] - T_mult * m[covariateindex, 2],
                                 m[covariateindex, 1] + T_mult * m[covariateindex,2]), 1, psthr,digits)
      pvalues <- c(sapply(m[covariateindex, 5], lpvalue))
    }
    if (length(betaname[[1]]) == 1) {
      if (!is.factor(data[, oldcovname])) {
        # title <- c(nicename(covariatename), hazardratio,pvalues, globalpvalue)
        title <- c(covariatename, hazardratio,pvalues, globalpvalue)
      }     else if (length(levelnames) == 1) {
        title <- c(covariatename, "", pvalues,globalpvalue)
        if (!is.null(data))
          reference <- c(addspace(sanitizestr(names(table(data[,
                                                               which(names(data) == oldcovname)]))[1])),
                         "Reference", "", "")
        body <- c(levelnames, hazardratio, "",
                  "")
      }      else {
        if (!is.null(data)) {
          reference <- c(addspace(sanitizestr(names(table(data[,
                                                               which(names(data) == oldcovname)]))[1])),
                         "Reference", "", "")
        }
        title <- c(covariatename, "", "",
                   globalpvalue)
        body <- cbind(levelnames, hazardratio, pvalues,
                      rep("", length(levelnames)))
      }
    }    else {
      if (length(levelnames) != 1) {
        title <- c(covariatename, "", "",
                   globalpvalue)
        body <- cbind(levelnames, hazardratio, pvalues,
                      rep("", length(levelnames)))
      }      else {
        title <- c(covariatename, hazardratio, pvalues,
                   globalpvalue)

      }
    }
    out <- rbind(title, reference, body)
    if (out[1, 2] == "") {
      if (length(grep(":", title[1])) > 0) {
        ss_N = unlist(lapply(levelnameslist,
                                   function(level) {
                                     N <- mapply(function(cn, lvl) {
                                       if (cn == lvl) {
                                         nrow(data)
                                       } else {
                                         sum(data[[cn]] == lvl)
                                       }
                                     }, oldcovname, level)
                                     return(min(N))
                                   }))
      }
      else {
        ss_N = as.vector(table(data[[oldcovname]]))
      }
      ss_N <- c(nrow(data),ss_N) # Add in the total for the variable
    }
    else {
      ss_N = nrow(data)
    }
    out <- cbind(out, ss_N)
    rownames(out) <- NULL
    colnames(out) <- NULL
    return(list(out, nrow(out)))
  })
  table <- lapply(out, function(x) {
    return(x[[1]])
  })
  index <- unlist(lapply(out, function(x) {
    return(x[[2]])
  }))
  table <- lapply(out, function(x) {
    return(x[[1]])
  })
  index <- unlist(lapply(out, function(x) {
    return(x[[2]])
  }))
  table <- do.call("rbind", lapply(table, data.frame,
                                   stringsAsFactors = FALSE))
  colName = c("Covariate", sanitizestr(beta), "p-value",
              "Global p-value","N")
  colnames(table) <- colName
  table[,"Global p-value"] <- ifelse(table[,'p-value']=='',table[,"Global p-value"],'')
  if (all(table[,"Global p-value"]=='')) table <- table[, -which(colnames(table)=="Global p-value")]
  if (!showN) table <- table[, -which(colnames(table)=="N")]
  if (vif) {
    if (type=='geeglm'|type=='lme'){
        message('VIF not yet implemented for mixed effects/GEE models.')
    } else {
      if (type=='crr'){
        xnm <- intersect(names(data),names(model$coef))
        data$y <- rowSums(data[,xnm],na.rm = TRUE)+stats::rnorm(nrow(data),0,2)
        mvif <- lm(formula = paste('y~',paste(xnm,collapse = '+')),data=data)
        VIF <- try(GVIF(mvif),silent = TRUE)
      } else VIF <- try(GVIF(model),silent = TRUE)
      if (!inherits(VIF,'try-error')) {
        if (nrow(VIF)>1){
        vifcol <- character(nrow(table))
        ind <- match(VIF$Covariate,table$Covariate)
        for (x in 1:length(ind)) vifcol[ind[x]] <- niceNum(VIF$VIF[x],digits = digits)
        table <- cbind(table,VIF=vifcol)
      }
        } else warning('VIF could not be computed for the model.')
    }}
  if (nicenames) table[,1] <- nicename(table[,1])
  colnames(table) <- sapply(colnames(table), lbld)
  attr(table,'covs') <- ucall
  return(table)
}




#' Create a forest plot using ggplot2
#'
#' This function will accept a log or logistic regression fit from glm or geeglm, and
#' display the OR or RR for each variable on the appropriate log scale.
#'
#' @param model an object output from the glm or geeglm function, must be from a logistic
#'   regression
#' @param conf.level controls the width of the confidence interval
#' @param orderByRisk logical, should the plot be ordered by risk
#' @param colours can specify colours for risks less than, 1 and greater than
#'   1.0. Default is red, black, green
#' @param showEst logical, should the risks be displayed on the plot in text
#' @param rmRef logical, should the reference levels be removed for the plot?
#' @param logScale logical, should OR/RR be shown on log scale, defaults to
#'   TRUE. See https://doi.org/10.1093/aje/kwr156 for why you may prefer a
#'   linear scale.
#' @param nxTicks Number of tick marks supplied to the log_breaks function to
#'   produce
#' @import ggplot2
#' @importFrom scales log_breaks
#' @keywords plot
#' @return a plot object
#' @export
#' @examples
#' glm_fit = glm(orr~change_ctdna_group+sex+age+l_size,
#' data=pembrolizumab,family = 'binomial')
#' forestplot2(glm_fit)
forestplot2 = function(model,conf.level=0.95,orderByRisk=TRUE,colours='default',showEst=TRUE,rmRef=FALSE,logScale=TRUE,nxTicks=5){

  if (inherits(model,'glm')){
    if(model$family$link=='log'){
      x_lab = 'Relative Risk'
    } else if (model$family$link=='logit'){
      x_lab='Odds Ratio'
    } else stop('model must be a logit or log link fit')
  } else {
    x_lab='Odds Ratio'
  }

  tab = format_glm(model,conf.level = conf.level,orderByRisk=orderByRisk)
  if (rmRef) tab = tab[setdiff(1:nrow(tab),which(tab$estimate.label=='1.0 (Reference)')),]

  yvals=1:nrow(tab)
  tab$estimate.label = ifelse(is.na(tab$estimate.label),'',tab$estimate.label)
  tab$estimate.label = ifelse(tab$estimate.label == '1.0 (Reference)','(Reference)',tab$estimate.label)

  if (showEst){
    yLabels = data.frame(y.pos=yvals,
                         labels=ifelse(is.na(tab$level.name),
                                       paste(tab$variable,tab$estimate.label),
                                       paste(tab$level.name,tab$estimate.label)))
  } else {
    yLabels = data.frame(y.pos=yvals,
                         labels=ifelse(is.na(tab$level.name),
                                       tab$variable,
                                       ifelse(tab$estimate.label == '(Reference)',
                                              paste(tab$level.name,tab$estimate.label),
                                              tab$level.name)
                         ))
  }
  yLabels$labels <- gsub('_',' ',yLabels$labels)
  yLabels <- yLabels[!is.na(yLabels$y.pos),]
  # TODO: add indents using '  ' to category labels, omit RR estimates?, make hjust=0

  tab$x.val = ifelse(tab$estimate.label == '(Reference)',1,tab$estimate)
  tab$y.val = yLabels$y.pos

  # set colours
  tab$colour <- ifelse(tab$x.val<1,'a',ifelse(tab$x.val==1,'b','c'))

  if (colours=='default'){
    colours = c(a='red',b='black',c='darkgreen')
  }  else {
    names(colours) = c('a','b','c')
  }

  # ensure that colours are always red, black, green
  colours <- colours[sort(unique(tab$colour))]

  p = ggplot(tab, aes_(x=~x.val,y=~y.val,colour=~colour))+
    geom_point(na.rm=TRUE,size=2) +
    geom_errorbarh(aes_(xmin = ~conf.low, xmax = ~conf.high),
                   height  = 0,
                   size   = 0.9,
                   na.rm=TRUE) +
    geom_vline(xintercept = 1.0) +
    labs(y='',x=x_lab) +
    guides(colour='none')+
    scale_y_continuous(breaks = yLabels$y.pos,labels=yLabels$labels) +
    scale_colour_manual(values=colours)+
    theme_bw() +
    theme(axis.text.y = element_text(face=ifelse(tab$variable ==tab$var.name | is.na(tab$var.name),"bold","plain"),
                                     hjust = 0),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank())

  if (logScale) p + scale_x_log10(breaks=scales::log_breaks(n=nxTicks)) else p
}


#' Create an univariable forest plot using ggplot2
#'
#' This function will send and take log or logistic regression fit from glm or geeglm
#' from uvsum function, and display the OR or RR for each variable on the appropriate log scale.
#'
#' @param response character vector with names of columns to use for response
#' @param covs character vector with names of columns to use for covariates
#' @param data dataframe containing your data
#' @param conf.level controls the width of the confidence interval
#' @param orderByRisk logical, should the plot be ordered by risk
#' @param colours can specify colours for risks less than, 1 and greater than
#'   1.0. Default is red, black, green
#' @param showEst logical, should the risks be displayed on the plot in text
#' @param rmRef logical, should the reference levels be removed for the plot?
#' @param logScale logical, should OR/RR be shown on log scale, defaults to
#'   TRUE. See https://doi.org/10.1093/aje/kwr156 for why you may prefer a
#'   linear scale.
#' @param nxTicks Number of tick marks supplied to the log_breaks function to
#'   produce
#' @param showN Show number of observations per variable and category
#' @import ggplot2
#' @importFrom scales log_breaks
#' @importFrom data.table ":="
#' @keywords plot
#' @return a plot object
#' @export
#' @examples
#' forestplotUV(response="orr", covs=c("change_ctdna_group", "sex", "age", "l_size"), data=pembrolizumab, family='binomial')
forestplotUV = function (response, covs, data, id = NULL, corstr = NULL, model = "glm",
                          family = NULL, digits = 2, conf.level = 0.95, orderByRisk = TRUE, colours = "default",
                          showEst = TRUE, rmRef = FALSE, logScale = FALSE, nxTicks = 5, showN = TRUE)
{
  if (inherits(model, "glm")) {
    if (model$family$link == "log") {
      x_lab = "Unadjusted Relative Risk"
    }
    else if (model$family$link == "logit") {
      x_lab = "Unadjusted Odds Ratio"
    }
    else stop("model must be a logit or log link fit")
  }
  else {
    x_lab = "Unadjusted Odds Ratio"
  }
  #tab = format_glm(model, conf.level = conf.level, orderByRisk = orderByRisk)
  ###################################
  tab = uvsum(response, covs, data, digits = 2, id = NULL, corstr = NULL,
                          family = NULL, type = NULL, gee = FALSE, strata = 1, markup = F,
                          sanitize = F, nicenames = F, showN = TRUE, CIwidth = conf.level,
                          reflevel = NULL, returnModels = FALSE)
  tab$estimate.label <- tab[,2];
  tab$estimate.label[which(tab$estimate.label == "Reference")] <- "1.0 (Reference)";
  tab$estimate <- as.numeric(gsub(" .*", "", tab[,2]));
  tab$conf.low <- as.numeric(gsub(",.*", "", gsub("\\(([^()]*)\\)|.", "\\1", tab[,2])));
  tab$conf.high <- as.numeric(gsub("^\\S*\\s+", "", gsub("\\(([^()]*)\\)|.", "\\1", tab[,2])));
  tab$level.name <- tab[,1];
  tab$var.name <- NA;
  tab$var.name[which(tab$Covariate %in% covs)] <- tab$level.name[which(tab$Covariate %in% covs)];
  dt <- data.table(y = tab$var.name )
  dt[, y_forward_fill := y[1], .(cumsum(!is.na(y)))]
  tab <- cbind(tab, dt)
  tab$var.name <- tab$y_forward_fill;
  tab$level.order <- sequence(rle(tab$var.name)$lengths)
  tab <- tab[order(rank(tab$estimate), tab$var.name), ];
  dt <- data.table(var.name = unique(tab$var.name));
  dt$var.order <- 1:nrow(dt);
  dt$var.order <- dt$var.order + 1;
  tab <- merge(tab, dt, by = "var.name", all = T);
  tab <- tab[order(tab$var.order, -tab$level.order), ];
  tab$p.value <- tab$"p-value";
  tab$p.label <- paste(format(round(as.numeric(tab$p.value), 3), nsmall=3), sep="");
  tab$variable <- tab$Covariate;
  tab <- tab[, c("variable", "var.name", "level.name", "level.order", "estimate", "p.label", "p.value", "conf.low", "conf.high", "var.order", "estimate.label", "N")];
  tab <- as.data.frame(tab);
  ###################################
  if (rmRef)
    tab = tab[setdiff(1:nrow(tab), which(tab$estimate.label ==
                                           "1.0 (Reference)")), ]
  yvals = 1:nrow(tab)
  tab$estimate.label = ifelse(is.na(tab$estimate.label), "",
                              tab$estimate.label)
  tab$estimate.label = ifelse(tab$estimate.label == "1.0 (Reference)",
                              "(Reference)", tab$estimate.label)
  if (showEst) {
    yLabels = data.frame(y.pos = yvals, labels = ifelse(is.na(tab$level.name),
                                                        paste(tab$variable, tab$estimate.label), paste(tab$level.name,
                                                                                                       tab$estimate.label)))
  }
  else {
    yLabels = data.frame(y.pos = yvals, labels = ifelse(is.na(tab$level.name),
                                                        tab$variable, ifelse(tab$estimate.label == "(Reference)",
                                                                             paste(tab$level.name, tab$estimate.label), tab$level.name)))
  }
  yLabels$labels <- gsub("_", " ", yLabels$labels)
  yLabels <- yLabels[!is.na(yLabels$y.pos), ]
  tab$x.val = ifelse(tab$estimate.label == "(Reference)", 1,
                     tab$estimate)
  tab$y.val = yLabels$y.pos
  tab$colour <- ifelse(tab$x.val < 1, "a", ifelse(tab$x.val ==
                                                    1, "b", "c"))
  if (colours == "default") {
    colours = c(a = "red", b = "black", c = "darkgreen")
  }
  else {
    names(colours) = c("a", "b", "c")
  }
  if (showN) {
    Axis = scale_y_continuous(breaks = yLabels$y.pos,
                              labels = yLabels$labels, sec.axis = dup_axis(breaks = yLabels$y.pos,
                                                                           labels = tab$N, name = "N"))
  }
  else {
    Axis = scale_y_continuous(breaks = yLabels$y.pos,
                              labels = yLabels$labels)
  }
  colours <- colours[sort(unique(tab$colour))]
  suppressWarnings({tryCatch({
    p = ggplot(tab, aes_(x = ~x.val, y = ~y.val, colour = ~colour)) +
      geom_point(na.rm = TRUE, size = 2) + geom_errorbarh(aes_(xmin = ~conf.low,
                                                               xmax = ~conf.high), height = 0, size = 0.9, na.rm = TRUE) +
      geom_vline(xintercept = 1) + labs(y = "", x = x_lab) +
      guides(colour = "none") + Axis + scale_colour_manual(values = colours) +
      theme_bw() + theme(axis.text.y = element_text(face = ifelse(tab$variable ==
                                                                    tab$var.name | is.na(tab$var.name), "bold", "plain"),
                                                    hjust = 0), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
                         axis.ticks = element_blank(),
                         axis.title.y.right = element_text(angle = 0, hjust = 0.5, vjust = 0.5) )
    if (logScale)
      p + scale_x_log10(breaks = scales::log_breaks(n = nxTicks))
    else p
  }, error=function(e){})})
}


#' Create a multivariable forest plot using ggplot2
#'
#' This function will send and take log or logistic regression fit from glm or geeglm
#' from mvsum function, and display the OR or RR for each variable on the appropriate log scale.
#'
#' @param model an object output from the glm or geeglm function, must be from a logistic
#'   regression
#' @param conf.level controls the width of the confidence interval
#' @param orderByRisk logical, should the plot be ordered by risk
#' @param colours can specify colours for risks less than, 1 and greater than
#'   1.0. Default is red, black, green
#' @param showEst logical, should the risks be displayed on the plot in text
#' @param rmRef logical, should the reference levels be removed for the plot?
#' @param logScale logical, should OR/RR be shown on log scale, defaults to
#'   TRUE. See https://doi.org/10.1093/aje/kwr156 for why you may prefer a
#'   linear scale.
#' @param nxTicks Number of tick marks supplied to the log_breaks function to
#'   produce
#' @param showN Show number of observations per variable and category
#' @import ggplot2
#' @importFrom scales log_breaks
#' @importFrom data.table ":="
#' @keywords plot
#' @return a plot object
#' @export
#' @examples
#' glm_fit = glm(orr~change_ctdna_group+sex+age+l_size,
#' data=pembrolizumab,family = 'binomial')
#' forestplotMV(glm_fit)
forestplotMV = function (model, conf.level = 0.95, orderByRisk = TRUE, colours = "default",
                          showEst = TRUE, rmRef = FALSE, logScale = FALSE, nxTicks = 5, showN = TRUE)
{
  if (inherits(model, "glm")) {
    if (model$family$link == "log") {
      x_lab = "Adjusted Relative Risk"
    }
    else if (model$family$link == "logit") {
      x_lab = "Adjusted Odds Ratio"
    }
    else stop("model must be a logit or log link fit")
  }
  else {
    x_lab = "Adjusted Odds Ratio"
  }
  #tab = format_glm(model, conf.level = conf.level, orderByRisk = orderByRisk)
  ###################################
  tab = mvsum(model, data, digits = 2,
                          markup = F, sanitize = F, nicenames = F, showN = TRUE, CIwidth = conf.level)
  tab$estimate.label <- tab[,2];
  tab$estimate.label[which(tab$estimate.label == "Reference")] <- "1.0 (Reference)";
  tab$estimate <- as.numeric(gsub(" .*", "", tab[,2]));
  tab$conf.low <- as.numeric(gsub(",.*", "", gsub("\\(([^()]*)\\)|.", "\\1", tab[,2])));
  tab$conf.high <- as.numeric(gsub("^\\S*\\s+", "", gsub("\\(([^()]*)\\)|.", "\\1", tab[,2])));
  tab$level.name <- tab[,1];
  tab$var.name <- NA;
  covs <- colnames(model$model);
  tab$var.name[which(tab$Covariate %in% covs)] <- tab$level.name[which(tab$Covariate %in% covs)];
  dt <- data.table(y = tab$var.name )
  dt[, y_forward_fill := y[1], .(cumsum(!is.na(y)))]
  tab <- cbind(tab, dt)
  tab$var.name <- tab$y_forward_fill;
  tab$level.order <- sequence(rle(tab$var.name)$lengths)
  tab <- tab[order(rank(tab$estimate), tab$var.name), ];
  dt <- data.table(var.name = unique(tab$var.name));
  dt$var.order <- 1:nrow(dt);
  dt$var.order <- dt$var.order + 1;
  tab <- merge(tab, dt, by = "var.name", all = T);
  tab <- tab[order(tab$var.order, -tab$level.order), ];
  tab$p.value <- tab$"p-value";
  tab$p.label <- paste(format(round(as.numeric(tab$p.value), 3), nsmall=3), sep="");
  tab$variable <- tab$Covariate;
  tab <- tab[, c("variable", "var.name", "level.name", "level.order", "estimate", "p.label", "p.value", "conf.low", "conf.high", "var.order", "estimate.label", "N")];
  tab <- as.data.frame(tab);
  ###################################
  if (rmRef)
    tab = tab[setdiff(1:nrow(tab), which(tab$estimate.label ==
                                           "1.0 (Reference)")), ]
  yvals = 1:nrow(tab)
  tab$estimate.label = ifelse(is.na(tab$estimate.label), "",
                              tab$estimate.label)
  tab$estimate.label = ifelse(tab$estimate.label == "1.0 (Reference)",
                              "(Reference)", tab$estimate.label)
  if (showEst) {
    yLabels = data.frame(y.pos = yvals, labels = ifelse(is.na(tab$level.name),
                                                        paste(tab$variable, tab$estimate.label), paste(tab$level.name,
                                                                                                       tab$estimate.label)))
  }
  else {
    yLabels = data.frame(y.pos = yvals, labels = ifelse(is.na(tab$level.name),
                                                        tab$variable, ifelse(tab$estimate.label == "(Reference)",
                                                                             paste(tab$level.name, tab$estimate.label), tab$level.name)))
  }
  yLabels$labels <- gsub("_", " ", yLabels$labels)
  yLabels <- yLabels[!is.na(yLabels$y.pos), ]
  tab$x.val = ifelse(tab$estimate.label == "(Reference)", 1,
                     tab$estimate)
  tab$y.val = yLabels$y.pos
  tab$colour <- ifelse(tab$x.val < 1, "a", ifelse(tab$x.val ==
                                                    1, "b", "c"))
  if (colours == "default") {
    colours = c(a = "red", b = "black", c = "darkgreen")
  }
  else {
    names(colours) = c("a", "b", "c")
  }
  if (showN) {
    Axis = scale_y_continuous(breaks = yLabels$y.pos,
                              labels = yLabels$labels, sec.axis = dup_axis(breaks = yLabels$y.pos,
                                                                           labels = tab$N, name = "N"))
  }
  else {
    Axis = scale_y_continuous(breaks = yLabels$y.pos,
                              labels = yLabels$labels)
  }
  colours <- colours[sort(unique(tab$colour))]
  suppressWarnings({tryCatch({
    p = ggplot(tab, aes_(x = ~x.val, y = ~y.val, colour = ~colour)) +
      geom_point(na.rm = TRUE, size = 2) + geom_errorbarh(aes_(xmin = ~conf.low,
                                                               xmax = ~conf.high), height = 0, size = 0.9, na.rm = TRUE) +
      geom_vline(xintercept = 1) + labs(y = "", x = x_lab) +
      guides(colour = "none") + Axis + scale_colour_manual(values = colours) +
      theme_bw() + theme(axis.text.y = element_text(face = ifelse(tab$variable ==
                                                                    tab$var.name | is.na(tab$var.name), "bold", "plain"),
                                                    hjust = 0), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
                         axis.ticks = element_blank(),
                         axis.title.y.right = element_text(angle = 0, hjust = 0.5, vjust = 0.5) )
    if (logScale)
      p + scale_x_log10(breaks = scales::log_breaks(n = nxTicks))
    else p
  }, error=function(e){})})
}


#' Combine an univariable and multivariable forest plot using ggplot2
#'
#' This function will take log or logistic regression fit forest plot output
#' from forestplotUV and forestplotMV functions and display the combined
#' adjusted and unadjusted OR or RR for each variable on the appropriate
#' log scale. Please note that total N and reference-level N is taken from
#' unadjusted model.
#'
#' @param UVmodel an UV model object output from the forestplotUV function
#' @param MVmodel a MV model object output from the forestplotMV function
#' @param orderByRisk logical, should the plot be ordered by risk
#' @param colours can specify colours for risks less than, 1 and greater than
#'   1.0. Default is red, black, green
#' @param showEst logical, should the risks be displayed on the plot in text
#' @param rmRef logical, should the reference levels be removed for the plot?
#' @param logScale logical, should OR/RR be shown on log scale, defaults to
#'   TRUE. See https://doi.org/10.1093/aje/kwr156 for why you may prefer a
#'   linear scale.
#' @param nxTicks Number of tick marks supplied to the log_breaks function to
#'   produce
#' @param showN Show number of observations per variable and category
#' @import ggplot2
#' @importFrom scales log_breaks
#' @importFrom data.table ":="
#' @keywords plot
#' @return a plot object
#' @export
#' @examples
#' UVp = forestplotUV(response="orr", covs=c("change_ctdna_group", "sex", "age",
#' "l_size"), data=pembrolizumab, family='binomial')
#' MVp = forestplotMV(glm(orr~change_ctdna_group+sex+age+l_size,
#' data=pembrolizumab,family = 'binomial'))
#' forestplotUVMV(UVp, MVp)
forestplotUVMV = function (UVmodel, MVmodel, model = "glm",
                            family = NULL, digits = 2, orderByRisk = TRUE, colours = "default",
                            showEst = TRUE, rmRef = FALSE, logScale = FALSE, nxTicks = 5, showN = TRUE)
{
  if (inherits(model, "glm")) {
    if (model$family$link == "log") {
      x_lab = "Relative Risk"
    }
    else if (model$family$link == "logit") {
      x_lab = "Odds Ratio"
    }
    else stop("model must be a logit or log link fit")
  }
  else {
    x_lab = "Odds Ratio"
  }
  #tab = format_glm(model, conf.level = conf.level, orderByRisk = orderByRisk)
  ###################################
  UVmodel$data$type <- "Unadjusted";
  MVmodel$data$type <- "Adjusted";
  tab <- rbind(UVmodel$data, MVmodel$data);
  tab <- tab[!duplicated(tab[ , c("var.name", "variable")]),];
  tab$Reference <- paste(tab$estimate.label, tab$type, sep="");
  tab[which(tab$estimate.label == "(Reference)"), ]$Reference <- "Y";
  tab <- tab[!duplicated(tab[, c("var.name", "estimate.label", "Reference")]),];
  dt <- data.table(y = tab$var.name )
  dt[, y_forward_fill := y[1], .(cumsum(!is.na(y)))]
  tab <- cbind(tab, dt)
  tab$var.name <- tab$y_forward_fill;
  tab$level.order <- sequence(rle(tab$var.name)$lengths)
  tab <- tab[order(rank(tab$estimate), tab$var.name), ];
  dt <- data.table(var.name = unique(tab$var.name));
  dt$var.order <- 1:nrow(dt);
  dt$var.order <- dt$var.order + 1;
  tab <- tab[, -which(colnames(tab) %in% c("var.order"))];
  tab <- merge(tab, dt, by = "var.name", all = T);
  tab <- tab[order(tab$var.order, tab$level.order, tab$type, decreasing = c(F, F, T), method="radix"), ];
  tab <- tab[, c("variable", "var.name", "level.name", "level.order", "estimate", "p.label", "p.value", "conf.low", "conf.high", "var.order", "estimate.label", "N", "type")];
  tab <- as.data.frame(tab);
  ###################################
  if (rmRef)
    tab = tab[setdiff(1:nrow(tab), which(tab$estimate.label ==
                                           "1.0 (Reference)")), ]
  yvals = 1:nrow(tab)
  tab$estimate.label = ifelse(is.na(tab$estimate.label), "",
                              tab$estimate.label)
  tab$estimate.label = ifelse(tab$estimate.label == "1.0 (Reference)",
                              "(Reference)", tab$estimate.label)
  if (showEst) {
    yLabels = data.frame(y.pos = yvals, labels = ifelse(is.na(tab$level.name),
                                                        paste(tab$variable, tab$estimate.label), paste(tab$level.name,
                                                                                                       tab$estimate.label)))
  }
  else {
    yLabels = data.frame(y.pos = yvals, labels = ifelse(is.na(tab$level.name),
                                                        tab$variable, ifelse(tab$estimate.label == "(Reference)",
                                                                             paste(tab$level.name, tab$estimate.label), tab$level.name)))
  }
  yLabels$labels <- gsub("_", " ", yLabels$labels)
  yLabels <- yLabels[!is.na(yLabels$y.pos), ]
  tab$x.val = ifelse(tab$estimate.label == "(Reference)", 1,
                     tab$estimate)
  tab$y.val = yLabels$y.pos
  tab$colour <- ifelse(tab$x.val < 1, "a", ifelse(tab$x.val ==
                                                    1, "b", "c"))
  if (colours == "default") {
    colours = c(a = "red", b = "black", c = "darkgreen")
  }
  else {
    names(colours) = c("a", "b", "c")
  }
  if (showN) {
    Axis = scale_y_continuous(breaks = yLabels$y.pos,
                              labels = yLabels$labels, sec.axis = dup_axis(breaks = yLabels$y.pos,
                                                                           labels = tab$N, name = "N"))
    warning(paste("Total N and reference-level N is taken from unadjusted model."))
  }
  else {
    Axis = scale_y_continuous(breaks = yLabels$y.pos,
                              labels = yLabels$labels)
  }
  colours <- colours[sort(unique(tab$colour))]
  suppressWarnings({tryCatch({
    p = ggplot(tab, aes_(x = ~x.val, y = ~y.val, colour = ~colour, linetype = ~type)) +
      geom_point(na.rm = TRUE, size = 2) + geom_errorbarh(aes_(xmin = ~conf.low,
                                                               xmax = ~conf.high), height = 0, size = 0.9, na.rm = TRUE) +
      geom_vline(xintercept = 1) + labs(y = "", x = x_lab) +
      guides(colour = "none") + Axis + scale_colour_manual(values = colours) +
      theme_bw() + theme(axis.text.y = element_text(face = ifelse(tab$variable ==
                                                                    tab$var.name | is.na(tab$var.name), "bold", "plain"),
                                                    hjust = 0), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
      theme(legend.title=element_blank(), legend.margin = margin(-8, 0, 0, 0),
            legend.spacing.x = unit(2, "mm"), legend.spacing.y = unit(0, "mm"),
            legend.position = "bottom", axis.ticks = element_blank(),
            axis.title.y.right = element_text(angle = 0, hjust = 0.5, vjust = 0.5) )
    if (logScale)
      p + scale_x_log10(breaks = scales::log_breaks(n = nxTicks))
    else p
  }, error=function(e){})})
}


#' Plot multiple bivariate relationships in a single plot
#'
#' This function is designed to accompany \code{\link{uvsum}} as a means of
#' visualising the results, and uses similar syntax.
#'
#' Plots are displayed as follows: If response is continuous For a numeric
#' predictor scatterplot For a categorical predictor: If 20+ observations
#' available boxplot, otherwise dotplot with median line If response is a factor
#' For a numeric predictor: If 20+ observations available boxplot, otherwise
#' dotplot with median line For a categorical predictor barplot Response
#' variables are shown on the ordinate (y-axis) and covariates on the abscissa
#' (x-axis)
#'
#' @param response character vector with names of columns to use for response
#' @param covs character vector with names of columns to use for covariates
#' @param data dataframe containing your data
#' @param showN boolean indicating whether sample sizes should be shown on the
#'   plots
#' @param showPoints boolean indicating whether individual data points should be
#'   shown when n>20 in a category
#' @param na.rm boolean indicating whether na values should be shown or removed
#' @param response_title character value with title of the plot
#' @param return_plotlist boolean indicating that the list of plots should be
#'   returned instead of a plot, useful for applying changes to the plot, see
#'   details
#' @param ncol the number of columns of plots to be display in the ggarrange
#'   call, defaults to 2
#' @param p_margins sets the TRBL margins of the individual plots, defaults to
#'   c(0,0.2,1,.2)
#' @param bpThreshold Default is 20, if there are fewer than 20 observations in
#'   a category then dotplots, as opposed to boxplots are shown.
#' @param mixed should a mix of dotplots and boxplots be shown based on sample
#'   size? If false then all categories will be shown as either dotplots, or
#'   boxplots according the bpThreshold and the smallest category size
#' @keywords plot
#' @returns a list containing plots for each variable in covs
#' @importFrom ggplot2 ggplot aes_string geom_boxplot geom_point geom_text
#'   stat_summary scale_x_discrete stat theme labs .data
#' @importFrom ggpubr ggarrange
#' @importFrom stats median
#' @return a plot object
#' @export
#' @examples
#' ## Run multiple univariate analyses on the pembrolizumab dataset to predict cbr and
#' ## then visualise the relationships.
#' rm_uvsum(data=pembrolizumab,
#' response='cbr',covs=c('age','sex','l_size','baseline_ctdna'))
#' plotuv(data=pembrolizumab,  response='cbr',
#' covs=c('age','sex','l_size','baseline_ctdna'),showN=TRUE)
#' @seealso \code{\link{ggplot}} and \code{\link{ggarrange}}
plotuv <- function(response,covs,data,showN=FALSE,showPoints=TRUE,na.rm=TRUE,
                   response_title=NULL,return_plotlist=FALSE,ncol=2,p_margins=c(0,0.2,1,.2),
                   bpThreshold=20,mixed=TRUE){
  for (v in c(response,covs)){
    if (!v %in% names(data)) stop(paste(v,'is not a variable in data.'))
    if (inherits(data[[v]],'character')) data[[v]] <- factor(data[[v]])
  }

  if (is.null(response_title)) response_title = response
  response_title = niceStr(response_title)
  plist <- NULL
  if (inherits(data[[response]],c('factor','ordered'))){
    use_common_legend = TRUE
    # ensure that all levels have the same colours for all plots
    lvls <- NULL
    for (x_var in covs){
      t <-table(data[[response]][!is.na(data[[x_var]])])
      lvls <- unique(c(lvls,names(t)[which(t>0)]))
    }
    levels(data[[response]]) <- c(levels(data[[response]])[which(levels(data[[response]])%in% lvls)],rep(NA,length(levels(data[[response]]))-length(lvls)))
    niceStr(levels(data[[response]]))
    lvlCol <- reportRx_pal()(length(levels(data[[response]])))
    names(lvlCol) = levels(data[[response]])
    for (x_var in covs){
      flip=FALSE
      # remove missing data, if requested
      if (na.rm) pdata = stats::na.omit(data[,c(response,x_var)]) else pdata = data[,c(response,x_var)]

      if (inherits(pdata[[x_var]],'numeric')){
        if (all(table(pdata[[response]])<bpThreshold) | (any(table(pdata[[response]])<bpThreshold) & !mixed)){
          p<-ggplot(data=pdata, aes(x=.data[[response]],y=.data[[x_var]],fill=.data[[response]]),colour=.data[[response]]) +
            geom_dotplot(binaxis = "y",stackdir = "center",dotsize = .8  ) +
            stat_summary(fun = median, fun.min = median, fun.max = median,
                         geom = "crossbar", width = 0.5) +
            coord_flip()
          flip = TRUE
        } else{
          if (any(table(pdata[[response]])<bpThreshold)){
            message(paste('Boxplots not shown for categories with fewer than', bpThreshold ,'observations.'))
          }
          pdata$alpha <- factor(ifelse(pdata[[response]] %in% names(table(pdata[[response]]))[table(pdata[[response]])<bpThreshold],'light','regular'),
                                levels=c('light','regular'))
          pdata$lty <- factor(ifelse(pdata[[response]] %in% names(table(pdata[[response]]))[table(pdata[[response]])<bpThreshold],'0','1'),
                              levels = c('0','1'))
          black_points <- pdata[!pdata[[response]] %in% names(table(pdata[[response]]))[table(pdata[[response]])<bpThreshold],]
          coloured_points <- pdata[pdata[[response]] %in% names(table(pdata[[response]]))[table(pdata[[response]])<bpThreshold],]
          p <- ggplot(data=pdata, aes(y=.data[[response]],x=.data[[x_var]],fill=.data[[response]])) +
            geom_boxplot(aes(alpha=.data[['alpha']],linetype=.data[['lty']]),outlier.shape = NA)  +
            scale_alpha_manual(breaks=c('light','regular'),values=c(0,1)) +
            scale_linetype_manual(breaks=c('0','1'),values = c(0,1))
          if (showPoints) {
            p <- p +geom_jitter(data=coloured_points,aes(colour=.data[[response]]), alpha=0.9)
            p <- p +geom_jitter(data= black_points,
                                color="black", size=0.4, alpha=0.9)
          }
          if (showN){
            p <-  p+
              stat_summary(aes(x=min(.data[[x_var]])),geom='label',vjust=-0.5,hjust=0,fun.data = lbl_count,label.size=0,fill='white',label.padding = unit(0.15, "lines"),alpha=.8)
          }

        }
        p<- p+
          theme(axis.text.y=element_blank(),
                axis.ticks.y = element_blank())
      } else {  # x_var is categorical
        p <- ggplot(data=pdata, aes(x=.data[[x_var]],fill=.data[[response]])) +
          geom_bar(position=position_dodge()) +
          scale_x_discrete(labels= function(x) wrp_lbl(x))
        if (showN){
          p <- p +
            geom_text(aes(label=stat(count)),position = position_dodge(width = 1),stat='count',vjust=1)
        }
        if (length(unique(pdata[[x_var]]))>8){
          p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
        }
      }
      p <- p  +
        theme_bw() +
        theme(
          plot.title = element_text(size=10),
          plot.margin = unit(p_margins, "lines")) +
        guides(alpha='none',linetype='none',colour='none')+
        scale_colour_manual(values=lvlCol)+
        scale_fill_manual(values=lvlCol)

      if (flip){
        plist[[x_var]] <- p +
          labs(y=niceStr(x_var),x='',fill=response_title)
      } else {
        plist[[x_var]] <- p +
          labs(x=niceStr(x_var),y='',fill=response_title)

      }
    }

  } else{ # Response is numeric
    use_common_legend = FALSE # colours have different meanings, indicated on x axis
    for (x_var in covs){
      # remove missing data, if requested
      if (na.rm) pdata = stats::na.omit(data[,c(response,x_var)]) else pdata = data[,c(response,x_var)]

      if (inherits(pdata[[x_var]],'numeric')){
        p <- ggplot(data=pdata, aes(y=.data[[response]],x=.data[[x_var]])) +
          geom_point()
      } else{
        if (all(table(pdata[[x_var]])<bpThreshold)){
          p <- ggplot(data=pdata, aes(x=.data[[x_var]],y=.data[[response]],fill=.data[[x_var]]),colour=.data[[x_var]]) +
            geom_dotplot(binaxis = "y",stackdir = "center" ,dotsize = .8) +
            stat_summary(fun = median, fun.min = median, fun.max = median,
                         geom = "crossbar", width = 0.5)+
            scale_x_discrete(labels= function(x) wrp_lbl(x))
        } else{
          if (any(table(pdata[[x_var]])<bpThreshold)){
            message(paste('Boxplots not shown for categories with fewer than', bpThreshold ,'observations.'))
          }
          pdata$alpha <- factor(ifelse(pdata[[x_var]] %in% names(table(pdata[[x_var]]))[table(pdata[[x_var]])<bpThreshold],'light','regular'),
                                c('light','regular'))
          pdata$lty <- factor(ifelse(pdata[[x_var]] %in% names(table(pdata[[x_var]]))[table(pdata[[x_var]])<bpThreshold],'0','1'),
                              levels=c('0','1'))
          black_points <- pdata[!pdata[[x_var]] %in% names(table(pdata[[x_var]]))[table(pdata[[x_var]])<bpThreshold],]
          coloured_points <- pdata[pdata[[x_var]] %in% names(table(pdata[[x_var]]))[table(pdata[[x_var]])<bpThreshold],]
          p <- ggplot(data=pdata, aes(x=.data[[x_var]],y=.data[[response]],fill=.data[[x_var]])) +
            geom_boxplot(aes(alpha=.data[['alpha']],linetype=.data[['lty']]),outlier.shape = NA)  +
            scale_alpha_manual(breaks=c('light','regular'),values=c(0,1)) +
            scale_linetype_manual(breaks=c('0','1'),values = c(0,1))+
            scale_x_discrete(labels= function(x) wrp_lbl(x))
          if (showPoints) {
            p <- p +geom_jitter(data=coloured_points,aes(colour=.data[[x_var]]), alpha=0.9)
            p <- p +geom_jitter(data= black_points,
                                color="black", size=0.4, alpha=0.9)
          }
          if (showN){
            p <-  p+
              stat_summary(aes(y=max(.data[[response]])),geom='label',fun.data = lbl_count,label.size=0,fill='white',label.padding = unit(0.15, "lines"),alpha=.8)
          }

        }
      }
      plist[[x_var]] <- p  +
        theme_bw() +
        theme(
          legend.position = 'none',
          plot.title = element_text(size=9),
          plot.margin = unit(p_margins, "lines")) +
        labs(x=niceStr(x_var),y=niceStr(response_title)) +
        guides(colour='none',linetype='none',alpha='none')+
        scale_colour_reportRx()

    }
  }
  # if the first plot doesn't have all the levels, take the legend from a plot that does
  if (inherits(data[[response]],c('factor','ordered'))){
    lvls_miss<-sapply(covs,function(x) length(setdiff(names(lvlCol),unique(data[[response]][!is.na(data[[x]])]))))
    if (lvls_miss[1]>0) legend.grob <- ggpubr::get_legend(plist[[which(lvls_miss==0)[1]]]) else legend.grob <- NULL
  } else legend.grob <- NULL
  if (return_plotlist){
    return(plist)
  } else{   suppressMessages(ggpubr::ggarrange(plotlist=plist,
                                               common.legend = use_common_legend,
                                               ncol=ncol,
                                               nrow=ceiling(length(plist)/ncol),
                                               legend.grob = legend.grob))
  }}

# Rmarkdown Reporting --------------------------------------------------------------

#' Print tables to PDF/Latex HTML or Word
#'
#' Output the table nicely to whatever format is appropriate. This is the output
#' function used by the rm_* printing functions.
#'
#' Entire rows can be bolded, or specific cells. Currently indentation refers to
#' the first column only. By default, underscores in column names are converted
#' to spaces. To disable this set rm_ to FALSE
#'
#' @param tab a table to format
#' @param row.names a string specifying the column name to assign to the
#'   rownames. If NULL (the default) then rownames are removed.
#' @param to_indent numeric vector indicating which rows to indent in the first
#'   column.
#' @param  bold_headers boolean indicating if the column headers should be
#'   bolded
#' @param rows_bold numeric vector indicating which rows to bold
#' @param bold_cells array indices indicating which cells to bold. These will be
#'   in addition to rows bolded by rows_bold.
#' @param caption table caption
#' @param digits number of digits to round numeric columns to, wither a single
#'   number or a vector corresponding to the number of numeric columns in tab
#' @param align string specifying column alignment, defaults to left alignment
#'   of the first column and right alignment of all other columns. The align
#'   argument accepts a single string with 'l' for left, 'c' for centre and 'r'
#'   for right, with no separations. For example, to set the left column to be
#'   centred, the middle column right-aligned and the right column left aligned
#'   use: align='crl'
#' @param applyAttributes boolean indicating if the function should use
#'   to_indent and bold_cells formatting attributes. This will only work
#'   properly if the dimensions of the table output from rm_covsum, rm_uvsum etc
#'   haven't changed.
#' @param keep.rownames should the row names be included in the output
#' @param fontsize PDF/HTML output only, manually set the table fontsize
#' @param chunk_label only used knitting to Word docs to allow cross-referencing
#' @return A character vector of the table source code, unless tableOnly=TRUE in
#'   which case a data frame is returned
#' @export
#' @examples
#' # To make custom changes or change the fontsize in PDF/HTML
#' tab <- rm_covsum(data=pembrolizumab,maincov = 'change_ctdna_group',
#' covs=c('age','sex','pdl1','tmb','l_size'),show.tests=TRUE,tableOnly = TRUE)
#' outTable(tab, fontsize=7)
#'
#' # To bold columns with the variable names
#'  rows_bold <- c(1,4,7,10,13)
#'  outTable(tab,rows_bold = rows_bold)
#'
#'  # To bold the estimates for male/female
#'  bold_cells <- as.matrix(expand.grid(5:6,1:ncol(tab)))
#'  outTable(tab,bold_cells= bold_cells)
outTable <- function(tab,row.names=NULL,to_indent=numeric(0),bold_headers=TRUE,
                     rows_bold=numeric(0),bold_cells=NULL,caption=NULL,digits,align,
                     applyAttributes=TRUE,keep.rownames=FALSE,fontsize,chunk_label){

  # strip tibble aspects
  tab=as.data.frame(tab)
  if (!is.null(row.names)) {
    tab <- cbind(rownames(tab),tab)
    names(tab)[1] <- row.names
  }
  rownames(tab) <- NULL

  # define column alignment
  if (missing(align)){
    alignSpec = paste(c('l',rep('r',ncol(tab)-1)),collapse = '',sep='')
  } else{
    alignSpec = gsub('[^lrc]+','',paste(align,collapse=''))
    alignSpec = substr(alignSpec,1,ncol(tab))
    if (nchar(alignSpec)<ncol(tab)) {
      lastchar = substr(alignSpec,nchar(alignSpec),nchar(alignSpec))
      alignSpec <- paste0(alignSpec,paste(rep(lastchar,ncol(tab)-nchar(alignSpec)),collapse=''))
    }
    if (!identical(alignSpec, align)){ warning(paste0('Argument align did not conform to expectations, align="',alignSpec,'" used instead'))}
  }
  # round and format numeric columns if digits is specified
  if (!missing(digits)){
    coltypes <- unlist(lapply(tab, class))
    if (any(coltypes=='numeric')){
      numCols <- names(coltypes)[coltypes=='numeric']
      colRound <- cbind(numCols,digits)
      colDigits <- as.numeric(colRound[,2])
      names(colDigits) <- colRound[,1]
      for (v in numCols) tab[[v]] <- sapply(tab[[v]],function(x) niceNum(x,digits=colDigits[v]))
    }
  }
  out_fmt = ifelse(is.null(knitr::pandoc_to()),'html',
                   ifelse(knitr::pandoc_to(c('doc','docx')),'doc',
                          ifelse(knitr::is_latex_output(),'latex','html')))

  chunk_label = ifelse(missing(chunk_label),'NOLABELTOADD',chunk_label)

  if (applyAttributes){
    if (!is.null(attr(tab,'dimchk'))){
      if (all(attr(tab,'dimchk')==dim(tab))){
        if (!is.null(attr(tab,'to_indent'))) to_indent <- attr(tab,'to_indent')
        if (!is.null(attr(tab,'bold_cells'))) bold_cells <- attr(tab,'bold_cells')
      }
    }
  }
  if (is.null(to_indent)) to_indent = numeric(0)
  to_indent = as.vector(to_indent)

  if (length(rows_bold)>0){
    arrInd <- as.matrix(expand.grid(rows_bold,1:ncol(tab)))
    bold_cells <- rbind(bold_cells,arrInd)
    dimnames(bold_cells) <- NULL
    bold_cells <- bold_cells[!duplicated(bold_cells),]
  }
  if (!is.null(bold_cells)){
    bold_cells <- bold_cells[!duplicated(bold_cells),,drop=FALSE]
    bold_cells <- bold_cells[!is.na(tab[bold_cells]),,drop=FALSE]
  }
  if (out_fmt=='doc'){
    caption = if (!is.null(caption)) {ifelse(chunk_label=='NOLABELTOADD',caption,paste0('(\\#tab:',chunk_label,')',caption))}
    tab[is.na(tab)] <-'&nbsp;' # This is necessary to assign the 'Compact' style to empty cells
    tab[tab==''] <-'&nbsp;'

    tab[[1]][to_indent] <- sapply(tab[[1]][to_indent],function(x) paste('&nbsp;&nbsp;',x))
    pander::pander(tab,
                   caption=caption,
                   emphasize.strong.cells=bold_cells,
                   split.table=Inf, split.cells=15,
                   justify = alignSpec)


  } else {  # For PDF, HTML
    # set NA to empty in kable
    oldop <- options()
    on.exit(options(oldop))
    options(knitr.kable.NA = '')
    if (out_fmt=='latex') {
      names(tab) <- sanitize(names(tab))
      if(!is.null(caption)) caption <- sanitize(caption)
      for (v in 1:ncol(tab)) tab[[v]] <- sanitize(tab[[v]])
      if (!is.null(bold_cells)) tab[bold_cells] <- sapply(tab[bold_cells],function(x) lbld(x))
    }
    if (out_fmt=='html') {
      if (!is.null(bold_cells)) tab[bold_cells] <- sapply(tab[bold_cells],function(x) hbld(x))
      for (v in 1:ncol(tab)) tab[[v]] <- rmds(tab[[v]])
    }
    #    names(tab) <- sanitize(names(tab))
    # This may not work as expected if a small table is split over pages
    # better to also repeat table headers
    if (nrow(tab)>30){
      kout <- knitr::kable(tab, format = out_fmt,
                           escape = FALSE,
                           booktabs=TRUE,
                           longtable=TRUE,
                           linesep='',
                           caption=caption,
                           align =alignSpec)
      kout <- kableExtra::kable_styling(kout,latex_options = c('repeat_header'))
    } else {
      kout <- knitr::kable(tab, format = out_fmt,
                           escape = FALSE,
                           booktabs=TRUE,
                           longtable=FALSE,
                           linesep='',
                           caption=caption,
                           align = alignSpec)
    }
    kout <- kableExtra::add_indent(kout,positions = to_indent)
    if (!missing(fontsize)){
      kout <- kableExtra::kable_styling(kout,font_size = fontsize)
    }
    if (out_fmt=='html'){
      kout <- kableExtra::kable_styling(kout,full_width = T)

    }
    kout
  }
}

#' Combine two table columns into a single column with levels of one nested
#' within levels of the other.
#'
#' This function accepts a data frame (via the data argument) and combines two
#' columns into a single column with values from the head_col serving as headers
#' and values of the to_col displayed underneath each header. The resulting
#' table is then passed to outTable for printing and output, to use the grouped
#' table as a data frame specify tableOnly=TRUE. By default the headers will be
#' bolded and the remaining values indented.
#'
#' Note that it is possible to combine multiple tables (more than two) with this
#' function.
#'
#' @param data dataframe
#' @param head_col character value specifying the column name with the headers
#' @param to_col character value specifying the column name to add the headers
#'   into
#' @param colHeader character with the desired name of the first column.
#'   The default is to leave this empty for output or, for table only output to
#'   use the column name 'col1'.
#' @param caption table caption
#' @param indent Boolean should the original values in the to_col be indented
#' @param boldheaders Boolean should the header column values be bolded
#' @param hdr_prefix character value that will prefix headers
#' @param hdr_suffix character value that will suffix headers
#' @param digits number of digits to round numeric columns to, wither a single
#'   number or a vector corresponding to the number of numeric columns
#' @param tableOnly boolean indicating if the table should be formatted for
#'   printing or returned as a data frame
#' @param fontsize PDF/HTML output only, manually set the table fontsize
#' @return A character vector of the table source code, unless tableOnly=TRUE in
#'   which case a data frame is returned
#' @export
#' @examples
#' ## Investigate models to predict baseline ctDNA and tumour size and display together
#' ## (not clinically useful!)
#' fit1 <- lm(baseline_ctdna~age+l_size+pdl1,data=pembrolizumab)
#' m1 <- rm_mvsum(fit1,tableOnly=TRUE)
#' m1$Response = 'ctDNA'
#' fit2 <- lm(l_size~age+baseline_ctdna+pdl1,data=pembrolizumab)
#' m2 <- rm_mvsum(fit2,tableOnly=TRUE)
#' m2$Response = 'Tumour Size'
#' rbind(m1,m2)
#' nestTable(rbind(m1,m2),head_col='Response',to_col='Covariate')
nestTable <- function(data,head_col,to_col,colHeader ='',caption=NULL,indent=TRUE,boldheaders=TRUE,hdr_prefix='',hdr_suffix='',digits=2,tableOnly=FALSE,fontsize){

  # strip any grouped data or tibble properties
  if (inherits(data,'data.frame')){
    colNames <- names(data)
    data <- data.frame(data)
  } else stop('data must be a data.frame')
  if (length(which(names(data)==head_col))==0) stop ('head_col must be a string specifying a variable in data')
  if (length(which(names(data)==to_col))==0) stop ('to_col must be a string specifying a variable in data')
  # re-order columns so that the head_col and to_col appear to the left
  colOrd <- c(which(names(data)==head_col),which(names(data)==to_col),
              setdiff(1:ncol(data),c(which(names(data)==head_col),which(names(data)==to_col))))
  data <- data[,colOrd]
  # ensure that the data are sorted by the header column and to column in the order they first appear
  # necessary if there is a misplaced row
  data[[head_col]] <- factor(data[[head_col]],levels=unique(data[[head_col]]),ordered = T)
  data[[to_col]] <- factor(data[[to_col]],levels=unique(data[[to_col]]),ordered = T)
  data <- data[order(data[[head_col]],data[[to_col]]),]
  data[[head_col]] <- as.character(data[[head_col]])
  data[[to_col]] <- as.character(data[[to_col]])
  new_row = data[1,]

  # round and format numeric columns if digits is specified
  if (!missing(digits)){
    coltypes <- unlist(lapply(data, class))
    numCols <- names(coltypes)[coltypes=='numeric']
    if (length(numCols)>0){
      colRound <- cbind(numCols,digits)
      colDigits <- as.numeric(colRound[,2])
      names(colDigits) <- colRound[,1]
      for (v in numCols) data[[v]] <- sapply(data[[v]],function(x) niceNum(x,digits=colDigits[v]))
    }
  }

  for (i in 1:ncol(new_row)) new_row[1,i] <- NA
  new_headers = unique(data[[head_col]])
  repeat{
    header_index = which(!duplicated(data[[head_col]]) & !is.na(data[[head_col]]))[1]
    new_row[[to_col]] <- data[[head_col]][header_index]

    if (header_index>1){
      data = rbind(data[1:(header_index-1),],new_row,data[(header_index):nrow(data),])
    } else {
      data = rbind(new_row,data)
    }

    data[[head_col]][data[[head_col]]==new_row[[to_col]]] <- NA
    if (sum(is.na(data[[head_col]]))==nrow(data)) break
  }
  header_rows <- which(data[[to_col]] %in% new_headers)
  to_indent <- which(!(data[[to_col]] %in% new_headers))
  if (!indent) to_indent <- numeric(0)

  data[[to_col]][header_rows] <- paste0(hdr_prefix,data[[to_col]][header_rows],hdr_suffix)

  data <- data[,setdiff(names(data),head_col)]
  names(data) <- c(colHeader,setdiff(colNames,c(to_col,head_col)))
  if (tableOnly){
    if (names(data)[1]=='')  names(data)[1] <- 'Col1'
    return(data)
  }
  if (boldheaders) rows_bold = header_rows else rows_bold=numeric(0)
  argL <- list(tab=data,to_indent=to_indent,rows_bold=rows_bold,caption=caption)
  if (!missing(fontsize)) argL[['fontsize']] <- fontsize
  do.call(outTable, argL)
}


#' Outputs a descriptive covariate table
#'
#' Returns a data frame corresponding to a descriptive table.
#'
#' Comparisons for categorical variables default to chi-square tests, but if
#' there are counts of <5 then the Fisher Exact test will be used and if this is
#' unsuccessful then a second attempt will be made computing p-values using MC
#' simulation. If testcont='ANOVA' then the t-test with unequal variance will be
#' used for two groups and an ANOVA will be used for three or more. The
#' statistical test used can be displayed by specifying show.tests=TRUE.
#'
#' Effect size can be obtained when p-value is requested.
#'
#' Further formatting options are available using tableOnly=TRUE and outputting
#' the table with a call to outTable.
#'
#' @param data dataframe containing data
#' @param covs character vector with the names of columns to include in table
#' @param maincov covariate to stratify table by
#' @param caption character containing table caption (default is no caption)
#' @param tableOnly Logical, if TRUE then a dataframe is returned, otherwise a
#'   formatted printed object is returned (default).
#' @param covTitle character with the names of the covariate (predictor) column.
#'   The default is to leave this empty for output or, for table only output to
#'   use the column name 'Covariate'.
#' @param digits number of digits for summarizing mean data
#' @param digits.cat number of digits for the proportions when summarizing
#'   categorical data (default: 0)
#' @param nicenames boolean indicating if you want to replace . and _ in strings
#'   with a space
#' @param IQR boolean indicating if you want to display the inter quantile range
#'   (Q1,Q3) as opposed to (min,max) in the summary for continuous variables
#' @param all.stats boolean indicating if all summary statistics (Q1,Q3 +
#'   min,max on a separate line) should be displayed. Overrides IQR.
#' @param pvalue boolean indicating if you want p-values included in the table
#' @param effSize boolean indicating if you want effect sizes included in the
#'   table. Can only be obtained if pvalue is also requested. Effect sizes are
#'   calculated with the rstatix package using Cramer V for categorical and Eta
#'   Squared for continuous covariates.
#' @param unformattedp boolean indicating if you would like the p-value to be
#'   returned unformatted (ie not rounded or prefixed with '<'). Best used with
#'   tableOnly = T and outTable function. See examples.
#' @param show.tests boolean indicating if the type of statistical used should
#'   be shown in a column beside the pvalues. Ignored if pvalue=FALSE.
#' @param testcont test of choice for continuous variables,one of
#'   \emph{rank-sum} (default) or \emph{ANOVA}
#' @param testcat test of choice for categorical variables,one of
#'   \emph{Chi-squared} (default) or \emph{Fisher}
#' @param full boolean indicating if you want the full sample included in the
#'   table, ignored if maincov is NULL
#' @param include_missing Option to include NA values of maincov. NAs will not
#'   be included in statistical tests
#' @param percentage choice of how percentages are presented, one of
#'   \emph{column} (default) or \emph{row}
#' @param dropLevels logical, indicating if empty factor levels be dropped from
#'   the output, default is TRUE.
#' @param excludeLevels a named list of covariate levels to exclude from
#'   statistical tests in the form list(varname =c('level1','level2')). These
#'   levels will be excluded from association tests, but not the table. This can
#'   be useful for levels where there is a logical skip (ie not missing, but not
#'   presented). Ignored if pvalue=FALSE.
#' @param numobs named list overriding the number of people you expect to have
#'   the covariate
#' @param fontsize PDF/HTML output only, manually set the table fontsize
#' @param chunk_label only used if output is to Word to allow cross-referencing
#' @keywords dataframe
#' @return A character vector of the table source code, unless tableOnly=TRUE in
#'   which case a data frame is returned
#' @export
#' @seealso \code{\link{covsum}},\code{\link{fisher.test}},
#'   \code{\link{chisq.test}}, \code{\link{wilcox.test}},
#'   \code{\link{kruskal.test}}, \code{\link{anova}}, \code{\link{cramer_v}},
#'   \code{\link{eta_squared}}, and \code{\link{outTable}}
#' @examples
#' rm_covsum(data=pembrolizumab, maincov = 'orr',
#' covs=c('age','sex','pdl1','tmb','l_size','change_ctdna_group'),
#' show.tests=TRUE)
#'
#' # To Show Effect Sizes
#' rm_covsum(data=pembrolizumab, maincov = 'orr',
#' covs=c('age','sex'),
#' effSize=TRUE)
#'
#' # To make custom changes or change the fontsize in PDF/HTML
#' tab <- rm_covsum(data=pembrolizumab,maincov = 'change_ctdna_group',
#' covs=c('age','sex','pdl1','tmb','l_size'),show.tests=TRUE,tableOnly = TRUE)
#' outTable(tab, fontsize=7)
#'
#' # To return unformatted p-values
#' tab <- rm_covsum(data=pembrolizumab, maincov = 'orr',
#' covs=c('age','sex','pdl1','tmb','l_size','change_ctdna_group'),
#' show.tests=TRUE,unformattedp=TRUE,tableOnly=TRUE)
#' outTable(tab,digits=5)
#' outTable(tab,digits=5, applyAttributes=FALSE) # remove bold/indent
rm_covsum <- function (data, covs, maincov = NULL, caption = NULL, tableOnly = FALSE,
                       covTitle = "", digits = 1, digits.cat = 0, nicenames = TRUE,
                       IQR = FALSE, all.stats = FALSE, pvalue = TRUE, effSize = FALSE, unformattedp = FALSE,
                       show.tests = FALSE, testcont = c("rank-sum test", "ANOVA"),
                       testcat = c("Chi-squared", "Fisher"), full = TRUE, include_missing = FALSE,
                       percentage = c("column", "row"), dropLevels = TRUE, excludeLevels = NULL,
                       numobs = NULL, fontsize,chunk_label)
{
  if (unformattedp)
    formatp <- function(x) {
      as.numeric(x)
    }
  argList <- as.list(match.call(expand.dots = TRUE)[-1])
  argsToPass <- intersect(names(formals(covsum)), names(argList))
  covsumArgs <- argList[names(argList) %in% argsToPass]
  covsumArgs[["markup"]] <- FALSE
  covsumArgs[["sanitize"]] <- FALSE
  tab <- do.call(covsum, covsumArgs)
  if (nicenames)
    output_var_names <- gsub("[_.]", " ", covs)
  else output_var_names <- covs
  vI <- unlist(sapply(output_var_names, function(x) which(x ==
                                                            tab[[1]])[1]))
  to_indent <- setdiff(1:nrow(tab), vI)
  to_bold_name <- vI
  if (nicenames)
    tab$Covariate <- gsub("[_.]", " ", tab$Covariate)
  names(tab)[1] <- covTitle
  bold_cells <- arrayInd(to_bold_name, dim(tab))
  if ("p-value" %in% names(tab)) {
    to_bold_p <- which(as.numeric(tab[["p-value"]]) < 0.05)
    p_vals <- tab[["p-value"]]
    new_p <- sapply(p_vals, formatp)
    tab[["p-value"]] <- new_p
    if (length(to_bold_p) > 0)
      bold_cells <- rbind(bold_cells, matrix(cbind(to_bold_p,
                                                   which(names(tab) == "p-value")), ncol = 2))
  }
  if ("Effect Size" %in% names(tab)) {
    e_vals <- tab[["Effect Size"]]
    new_e <- sapply(e_vals, formatp)
    tab[["Effect Size"]] <- new_e
  }
  if (tableOnly) {
    if (names(tab)[1] == "")
      names(tab)[1] <- "Covariate"
    attr(tab, "to_indent") <- to_indent
    attr(tab, "bold_cells") <- bold_cells
    attr(tab, "dimchk") <- dim(tab)
    return(tab)
  }
  argL <- list(tab = tab, to_indent = to_indent, bold_cells = bold_cells,
           caption = caption, chunk_label = ifelse(missing(chunk_label),
                                                   "NOLABELTOADD", chunk_label))
  if (!missing(fontsize)) argL[['fontsize']] <- fontsize
  do.call(outTable, argL)
}

#' Output several univariate models nicely in a single table
#'
#' A table with the model parameters from running separate univariate models on
#' each covariate. For factors with more than two levels a Global p-value is
#' returned.
#'
#' Global p-values are likelihood ratio tests for lm, glm and polr models. For
#' lme models an attempt is made to re-fit the model using ML and if,successful
#' LRT is used to obtain a global p-value. For coxph models the model is re-run
#' without robust variances with and without each variable and a LRT is
#' presented. If unsuccessful a Wald p-value is returned. For GEE and CRR models
#' Wald global p-values are returned.
#'
#' The number of decimals places to display the statistics can be changed with
#' digits, but this will not change the display of p-values. If more significant
#' digits are required for p-values then use tableOnly=TRUE and format as
#' desired.
#' @param response string vector with name of response
#' @param covs character vector with the names of columns to fit univariate
#'   models to
#' @param data dataframe containing data
#' @param digits number of digits to round estimates and CI to. Does not affect
#'   p-values.
#' @param covTitle character with the names of the covariate (predictor) column.
#'   The default is to leave this empty for output or, for table only output to
#'   use the column name 'Covariate'.
#' @param caption character containing table caption (default is no caption)
#' @param tableOnly boolean indicating if unformatted table should be returned
#' @param removeInf boolean indicating if infinite estimates should be removed
#'   from the table
#' @param p.adjust p-adjustments to be performed (Global p-values only)
#' @param unformattedp boolean indicating if you would like the p-value to be
#'   returned unformatted (ie not rounded or prefixed with '<'). Should be used
#'   in conjunction with the digits argument.
#' @param chunk_label only used if output is to Word to allow cross-referencing
#' @param  gee boolean indicating if gee models should be fit to account for
#'   correlated observations. If TRUE then the id argument must specify the
#'   column in the data which indicates the correlated clusters.
#' @param id character vector which identifies clusters. Only used for geeglm
#' @param corstr character string specifying the correlation structure. Only
#'   used for geeglm. The following are permitted: '"independence"',
#'   '"exchangeable"', '"ar1"', '"unstructured"' and '"userdefined"'
#' @param family description of the error distribution and link function to be
#'   used in the model. Only used for geeglm
#' @param type string indicating the type of univariate model to fit. The
#'   function will try and guess what type you want based on your response. If
#'   you want to override this you can manually specify the type. Options
#'   include "linear", "logistic", "poisson",coxph", "crr", "boxcox", "ordinal",
#'   "geeglm"
#' @param strata character vector of covariates to stratify by. Only used for
#'   coxph and crr
#' @param nicenames boolean indicating if you want to replace . and _ in strings
#'   with a space
#' @param showN boolean indicating if you want to show sample sizes
#' @param CIwidth width of confidence interval, default is 0.95
#' @param reflevel manual specification of the reference level. Only used for
#'   ordinal regression This will allow you to see which model is not fitting if
#'   the function throws an error
#' @param returnModels boolean indicating if a list of fitted models should be
#'   returned. If this is TRUE then the models will be returned, but the output
#'   will be suppressed. In addition to the model elements a data element will
#'   be appended to each model so that the fitted data can be examined, if
#'   necessary. See Details
#' @param fontsize PDF/HTML output only, manually set the table fontsize
#' @seealso
#' \code{\link{uvsum}},\code{\link{lm}},\code{\link{glm}},\code{\link{crr}},
#' \code{\link{coxph}},
#' \code{\link{lme}},\code{\link{geeglm}},\code{\link{polr}}
#' @return A character vector of the table source code, unless tableOnly=TRUE in
#'   which case a data frame is returned
#' @export
#' @examples
#' # Examples are for demonstration and are not meaningful
#' # Coxph model with 90% CI
#' rm_uvsum(response = c('os_time','os_status'),
#' covs=c('age','sex','baseline_ctdna','l_size','change_ctdna_group'),
#' data=pembrolizumab,CIwidth=.9)
#'
#' # Linear model with default 95% CI
#' rm_uvsum(response = 'baseline_ctdna',
#' covs=c('age','sex','l_size','pdl1','tmb'),
#' data=pembrolizumab)
#'
#' # Logistic model with default 95% CI
#' rm_uvsum(response = 'os_status',
#' covs=c('age','sex','l_size','pdl1','tmb'),
#' data=pembrolizumab,family = binomial)

#' # Poisson models returned as model list
#' mList <- rm_uvsum(response = 'baseline_ctdna',
#' covs=c('age','sex','l_size','pdl1','tmb'),
#' data=pembrolizumab, returnModels=TRUE)
#' # mList$sex$data # will expose the modelled data
#'
#' # GEE on correlated outcomes
#' rm_uvsum(response = 'size_change',
#' covs=c('time','ctdna_status'),
#' gee=TRUE,
#' id='id', corstr="exchangeable",
#' family=gaussian("identity"),
#' data=ctDNA,showN=TRUE)
rm_uvsum <- function(response, covs , data , digits=2, covTitle='',caption=NULL,
                     tableOnly=FALSE,removeInf=FALSE,p.adjust='none',unformattedp=FALSE,
                     chunk_label,
                     gee=FALSE,id = NULL,corstr = NULL,family = NULL,type = NULL,
                     strata = 1,
                     nicenames = TRUE,showN = TRUE,CIwidth = 0.95,
                     reflevel=NULL,returnModels=FALSE,fontsize){

  if (missing(data)) stop('data is a required argument')
  if (missing(covs)) stop('covs is a required argument')
  if (missing(response)) stop('response is a required argument')
  if (length(response)>2) stop('The response must be a single outcome for linear, logistic and ordinal models or must specify the time and event status variables for survival models.')
  if (!inherits(data,'data.frame')) stop('data must be supplied as a data frame.')
  if (!inherits(covs,'character')) stop('covs must be supplied as a character vector or string indicating variables in data')
  missing_vars = na.omit(setdiff(c(response, covs,id,ifelse(strata==1,NA,strata)), names(data)))
  if (length(missing_vars) > 0) stop(paste("These variables are not in the data:\n",
                                           paste0(missing_vars,collapse=csep())))
  if (strata==1) nm <- c(response,covs) else nm <- c(strata,response,covs)
  if (!all(names(data[,nm])==names(data.frame(data[,nm])))) stop('Non-standard variable names detected.\n Try converting data with new_data <- data.frame(data) \n then use new variable names in rm_uvsum.' )

  for (v in covs) {
    if (inherits(data[[v]], c("character", "ordered"))) data[[v]] <- factor(data[[v]], ordered = F)
    if (inherits(data[[v]],c('Date','POSIXt'))) {
      covs <- setdiff(covs,v)
      message(paste('Dates can not be used as predictors, try creating a time variable.\n The variable',v,'does not appear in the table.'))
    }

    df <- na.omit(data[,c(response,v)])
    if (v %in% response){
      warning(paste(v,'is the response and can not appear in the covariate.\n',
                    'It is omitted from the output.'))
      covs <- setdiff(covs,v)
    }
    if (length(unique(df[[v]]))==1) {
      warning(paste(v,'has only one unique value for non-missing response.\n',
                    'It is omitted from the output.'))
      covs <- setdiff(covs,v)
    }
  }

  if (unformattedp) formatp <- function (x,...){x}
  # get the table
  rtn <- uvsum(response,covs,data,digits=digits,markup = FALSE,sanitize=FALSE,
               gee=gee,id = id,
               corstr = corstr,family = family,type = type,strata = strata,
               nicenames = nicenames,showN = showN,
               CIwidth = CIwidth,reflevel=reflevel,returnModels=returnModels)
  if (returnModels) tab <- rtn[[1]] else tab <- rtn
  cap_warn <- character(0)
  if (removeInf){
    # Do not display unstable estimates
    inf_values =  grep('Inf',tab[,2])
    if (length(inf_values)>0){
      if ('Global p-values' %in% names(tab)) to_hide <-2:4 else to_hide <-2:3
      tab[inf_values,to_hide] <-NA
      cap_warn <- paste0(cap_warn,ifelse(identical(cap_warn,character(0)),'',', '),
                         'Covariates with unstable estimates:',
                         paste(tab$Covariate[inf_values],collapse=','),'.')
    }
  }
  # if an adjustment was made, add this to the cap_warn text
  if (p.adjust!='none') cap_warn <- paste0(cap_warn,'. Global p-values were adjusted according to the ',p.adjust,' method.')

  if (nicenames) output_var_names <- gsub('[_.]',' ',covs) else output_var_names <- covs
#  vI <- unlist(sapply(output_var_names, function (x) grep(x,tab[[1]])[1]))
  vI <- unlist(sapply(output_var_names, function (x) which(x==tab[[1]])[1]))
  to_indent <- setdiff(1:nrow(tab),vI)
  to_bold_name <- vI
  if (nicenames) tab$Covariate <- gsub('[_.]',' ',tab$Covariate)
  bold_cells <- arrayInd(to_bold_name, dim(tab))

  if ("Global p-value" %in% names(tab)){
    tab[["Global p-value"]][which(tab[["Global p-value"]]==''|tab[["Global p-value"]]=='NA')] <-NA
  }

  # perform p-value adjustment across all p-values
  if ("Global p-value" %in% names(tab)){
    raw_p <- ifelse(is.na(tab[["Global p-value"]]),tab[["p-value"]],tab[["Global p-value"]])
    p_sig <- suppressWarnings(stats::p.adjust(raw_p,method=p.adjust))
    to_bold_p <- which(p_sig<0.05)
    p_sig <- sapply(p_sig,formatp)
    gp_vals <- which(!is.na(tab[["Global p-value"]]))
    tab[["Global p-value"]][gp_vals]  <- p_sig[gp_vals]
    p_vals <- which(tab[["p-value"]]!='' & !is.na(tab[["p-value"]]))
    tab[["p-value"]][p_vals]  <- p_sig[p_vals]

    if (length(to_bold_p)>0) bold_cells <- rbind(bold_cells,
                                                 matrix(cbind(to_bold_p, which(names(tab)=='Global p-value')),ncol=2))

  } else {
    raw_p <- tab[["p-value"]]
    p_sig <- suppressWarnings(stats::p.adjust(raw_p,method=p.adjust))
    tab[["p-value"]] <- sapply(p_sig,formatp)
  }
  if(p.adjust!='none') tab[["raw p-value"]]<-formatp(raw_p)

  to_bold_p <- which(p_sig<.05)

  if (length(to_bold_p)>0) bold_cells <- rbind(bold_cells,
                                               matrix(cbind(to_bold_p, which(names(tab)=='p-value')),ncol=2))

  names(tab)[1] <-covTitle

  if (tableOnly){
    if (names(tab)[1]=='') names(tab)[1]<- 'Covariate'
    if (length(cap_warn)>0) message(cap_warn)
    attr(tab, 'to_indent') <- to_indent
    attr(tab,'bold_cells') <- bold_cells
    attr(tab,'dimchk') <- dim(tab)
    return(tab)
  }
  if (returnModels) return (rtn$models)
  argL <- list(tab=tab, digits = digits,
           to_indent=to_indent,bold_cells=bold_cells,
           caption=caption,
           chunk_label=ifelse(missing(chunk_label),'NOLABELTOADD',chunk_label))
  if (!missing(fontsize)) argL[['fontsize']] <- fontsize
  do.call(outTable, argL)

}


#' Format a regression model nicely for 'Rmarkdown'
#'
#' Multivariable (or univariate) regression models are re-formatted for
#' reporting and a global p-value is added for the evaluation of factor
#' variables.
#'
#' Global p-values are likelihood ratio tests for lm, glm and polr models. For
#' lme models an attempt is made to re-fit the model using ML and if,successful
#' LRT is used to obtain a global p-value. For coxph models the model is re-run
#' without robust variances with and without each variable and a LRT is
#' presented. If unsuccessful a Wald p-value is returned. For GEE and CRR models
#' Wald global p-values are returned.
#'
#' If the variance inflation factor is requested (VIF=T) then a generalised VIF
#' will be calculated in the same manner as the car package.
#'
#' The number of decimals places to display the statistics can be changed with
#' digits, but this will not change the display of p-values. If more significant
#' digits are required for p-values then use tableOnly=TRUE and format as
#' desired.
#' @param model model fit
#' @param data data that model was fit on (an attempt will be made to extract
#'   this from the model)
#' @param digits number of digits to round estimates to, does not affect
#'   p-values
#' @param covTitle character with the names of the covariate (predictor) column.
#'   The default is to leave this empty for output or, for table only output to
#'   use the column name 'Covariate'.
#' @param showN boolean indicating sample sizes should be shown for each
#'   comparison, can be useful for interactions
#' @param CIwidth width for confidence intervals, defaults to 0.95
#' @param vif boolean indicating if the variance inflation factor should be
#'   included. See details
#' @param caption table caption
#' @param tableOnly boolean indicating if unformatted table should be returned
#' @param p.adjust p-adjustments to be performed (Global p-values only)
#' @param unformattedp boolean indicating if you would like the p-value to be
#'   returned unformatted (ie not rounded or prefixed with '<'). Should be used
#'   in conjuction with the digits argument.
#' @param nicenames boolean indicating if you want to replace . and _ in strings
#'   with a space
#' @param chunk_label only used if output is to Word to allow cross-referencing
#' @param fontsize PDF/HTML output only, manually set the table fontsize
#' @return A character vector of the table source code, unless tableOnly=TRUE in
#'   which case a data frame is returned
#' @export
#' @references John Fox & Georges Monette (1992) Generalized Collinearity
#'   Diagnostics, Journal of the American Statistical Association, 87:417,
#'   178-183, DOI: 10.1080/01621459.1992.10475190
#' @references  John Fox and Sanford Weisberg (2019). An {R} Companion to
#'   Applied Regression, Third Edition. Thousand Oaks CA: Sage. URL:
#'   https://socialsciences.mcmaster.ca/jfox/Books/Companion
#' @examples
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
rm_mvsum <- function(model, data, digits=2,covTitle='',showN=FALSE,CIwidth=0.95, vif=TRUE,
                     caption=NULL,tableOnly=FALSE,p.adjust='none',unformattedp=FALSE,nicenames = TRUE,chunk_label, fontsize){
  if (unformattedp) formatp <- function(x) {as.numeric(x)}
  # get the table
  tab <- mvsum(model=model,data=data,digits=digits,markup = FALSE,
               sanitize = FALSE, nicenames = FALSE,showN=showN,CIwidth = CIwidth,vif=vif)

  to_indent <- setdiff(1:nrow(tab),
                       sapply(attr(tab,'covs'),function(x) grep(x,tab$Covariate)[1],
                              USE.NAMES = FALSE))

  if ("Global p-value" %in% names(tab)){
    tab[["Global p-value"]][which(tab[["Global p-value"]]==''|tab[["Global p-value"]]=='NA')] <-NA
    to_indent <- setdiff(to_indent,which(!is.na(tab[["Global p-value"]])))
  }
  to_bold_name <- setdiff(1:nrow(tab),to_indent)
  bold_cells <- arrayInd(to_bold_name, dim(tab))

  # perform p-value adjustment across all p-values
  if ("Global p-value" %in% names(tab)){
    raw_p <- ifelse(is.na(tab[["Global p-value"]]),tab[["p-value"]],tab[["Global p-value"]])
    p_sig <- suppressWarnings(stats::p.adjust(raw_p,method=p.adjust))
    p_sig <- sapply(p_sig,formatp)
    tab[["Global p-value"]][!is.na(tab[["Global p-value"]])]  <- p_sig[!is.na(tab[["Global p-value"]])]
    tab[["p-value"]][is.na(tab[["Global p-value"]])]  <- p_sig[is.na(tab[["Global p-value"]])]
    to_bold_p <- which(p_sig<0.05)
    if (length(to_bold_p)>0) bold_cells <- rbind(bold_cells,
                                                 matrix(cbind(to_bold_p, which(names(tab)=='Global p-value')),ncol=2))

  } else {
    raw_p <- tab[["p-value"]]
    p_sig <- suppressWarnings(stats::p.adjust(raw_p,method=p.adjust))
    tab[["p-value"]] <- sapply(p_sig,formatp)
  }
  if(p.adjust!='none') tab[["raw p-value"]]<-formatp(raw_p)

  to_bold_p <- which(p_sig<.05)

  if (length(to_bold_p)>0)  bold_cells <- rbind(bold_cells,
                                                matrix(cbind(to_bold_p, which(names(tab)=='p-value')),ncol=2))


  if (nicenames) tab$Covariate <- gsub('[_.]',' ',tab$Covariate)
  names(tab)[1] <-covTitle
  if (tableOnly){
    if (names(tab)[1]=='') names(tab)[1]<- 'Covariate'
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

#' Combine univariate and multivariable regression tables
#'
#' This function will combine rm_uvsum and rm_mvsum outputs into a single table.
#' The tableOnly argument must be set to TRUE when tables to be combined are
#' created. The resulting table will be in the same order as the uvsum table
#' and will contain the same columns as the uvsum and mvsum tables, but the
#' p-values will be combined into a single column. There must be a variable overlapping
#' between the uvsum and mvsum tables and all variables in the mvsum table
#' must also appear in the uvsum table.
#'
#'
#' @param uvsumTable Output from rm_uvsum, with tableOnly=TRUE
#' @param mvsumTable  Output from rm_mvsum, with tableOnly=TRUE
#' @param covTitle character with the names of the covariate (predictor) column.
#'   The default is to leave this empty for output or, for table only output to
#'   use the column name 'Covariate'.
#' @param vif boolean indicating if the variance inflation factor should be
#'   shown if present in the mvsumTable. Default is FALSE.
#' @param caption table caption
#' @param tableOnly boolean indicating if unformatted table should be returned
#' @param chunk_label only used if output is to Word to allow cross-referencing
#' @param fontsize PDF/HTML output only, manually set the table fontsize
#' @seealso
#'   \code{\link{rm_uvsum}},\code{\link{rm_mvsum}}
#' @return A character vector of the table source code, unless tableOnly=TRUE in
#'   which case a data frame is returned
#' @export
#' @examples
#' require(survival)
#'
#' uvTab <- rm_uvsum(response = c('os_time','os_status'),
#' covs=c('age','sex','baseline_ctdna','l_size','change_ctdna_group'),
#' data=pembrolizumab,tableOnly=TRUE)
#' mv_surv_fit <- coxph(Surv(os_time,os_status)~age+sex+
#' baseline_ctdna+l_size+change_ctdna_group, data=pembrolizumab)
#' uvTab <- rm_mvsum(mv_surv_fit)
#'
#' #linear model
#' uvtab<-rm_uvsum(response = 'baseline_ctdna',
#' covs=c('age','sex','l_size','pdl1','tmb'),
#' data=pembrolizumab,tableOnly=TRUE)
#' lm_fit=lm(baseline_ctdna~age+sex+l_size+tmb,data=pembrolizumab)
#' mvtab<-rm_mvsum(lm_fit,tableOnly = TRUE)
#' rm_uv_mv(uvtab,mvtab,tableOnly=TRUE)
#'
#' #logistic model
#' uvtab<-rm_uvsum(response = 'os_status',
#' covs=c('age','sex','l_size','pdl1','tmb'),
#' data=pembrolizumab,family = binomial,tableOnly=TRUE)
#' logis_fit<-glm(os_status~age+sex+l_size+pdl1+tmb,data = pembrolizumab,family = 'binomial')
#' mvtab<-rm_mvsum(logis_fit,tableOnly = TRUE)
#' rm_uv_mv(uvtab,mvtab,tableOnly=TRUE)
rm_uv_mv <- function(uvsumTable,mvsumTable,covTitle='',vif=FALSE,caption=NULL,tableOnly=FALSE,chunk_label,fontsize){
  # Check that tables are data frames and not kable objects
  if (!inherits(uvsumTable,'data.frame')) stop('uvsumTable must be a data.frame. Did you forget to specify tableOnly=TRUE?')
  if (!inherits(mvsumTable,'data.frame')) stop('mvsumTable must be a data.frame. Did you forget to specify tableOnly=TRUE?')
  # Check that the first columns have the same name
  if (names(uvsumTable)[1] != names(mvsumTable)[1]) stop('The covariate columns must have the same name in both tables')
  # Check that there is overlap between the variables
  if (length(intersect(uvsumTable[,1],mvsumTable[,1]))==0) stop('There are no overlaping variables between the models, tables couldn\'t be combined.')
  # Check that all the variables in the multivariate model are in the univariate model
  if (length(setdiff(mvsumTable[,1],uvsumTable[,1]))>0) {
    stop(paste('The following variables were not in the univariate model:',paste0(setdiff(mvsumTable[,1],uvsumTable[,1]),collapse=", "),
               '\nRun uvsum with all the variables in the multivariable model.'))
  }
  if (is.null(attr(uvsumTable,'to_indent'))) {
    warning('Please re-generate the uvsumTable to correct variable indenting')
    to_indent <- numeric(0)
  } else{
    to_indent <- attr(uvsumTable,'to_indent')
  }
  if (is.null(attr(mvsumTable,'covs'))){
    stop('Please re-generate the mvsumTable to identify the covariates')
  } else {
    attr(mvsumTable,'indented') <- which(!nicename(mvsumTable[,1]) %in% nicename(attr(mvsumTable,'covs')))
  }

  x <- lapply(list(uvsumTable,mvsumTable), function(t) {
    p_cols <- grep('p-value',names(t))
    # add a column for the variable name
    vname <- character(nrow(t))
    vname[setdiff(1:nrow(t),attr(t,'to_indent'))] <- t[,1][setdiff(1:nrow(t),attr(t,'to_indent'))]
    for (i in 1:nrow(t)) vname[i] <- ifelse(vname[i]=='',vname[i-1],vname[i])
    if ('Global p-value' %in% names(t)){
      t[['Global p-value']][t[['Global p-value']]==''] <- NA
      p_var <- ifelse(is.na(t[['Global p-value']]),t[['p-value']] ,t[['Global p-value']])
    } else p_var <- t[['p-value']]
    p_var <- ifelse(p_var=='',NA,p_var)
    t$p <- p_var
    t$var_level <- paste(vname,t[,1],sep='_')
    return(t[,setdiff(1:ncol(t),p_cols)])
  })
  x[[1]]$varOrder = 1:nrow(x[[1]])
  names(x[[1]])[2] <- paste('Unadjusted',names(x[[1]])[2])
  names(x[[2]])[2] <- paste('Adjusted',names(x[[2]])[2])
  vifind <- which(names(x[[2]])=='VIF')
  if (length(vifind)>0){
    if (!vif) x[[2]] <- x[[2]][,-vifind] else x[[2]] <- x[[2]][,c(setdiff(1:ncol(x[[2]]),vifind),vifind)]
  }
  for (vn in setdiff(names(x[[2]])[3:ncol(x[[2]])],c('VIF','var_level'))) names(x[[2]]) <- gsub(vn, paste(vn,'(adj)'),names(x[[2]]))
  out <- merge(x[[1]],x[[2]],by='var_level',all=TRUE)
  out <- out[,-which(names(out)=='var_level')]
  out <- out[,-grep('[.]y',names(out))]
  names(out) <- gsub('[.]x','',names(out))
  out <- out[order(out$varOrder),-which(names(out)=='varOrder')]

  names(out)[1] <-covTitle
  if (tableOnly){
    if (names(out)[1]=='') names(out)[1]<- 'Covariate'
    return(out)
  }

  to_bold_name <- setdiff(1:nrow(out),to_indent)
  bold_cells <- arrayInd(to_bold_name, dim(out))

  to_bold_p <- which(out[["p"]]<.05 & !is.na(out[["p"]]))
  if (length(to_bold_p)>0) bold_cells <- rbind(bold_cells,
                                               matrix(cbind(to_bold_p, which(names(out)=='p')),ncol=2))
  to_bold_p <- which(out[["p (adj)"]]<.05 & !is.na(out[["p (adj)"]]))
  if (length(to_bold_p)>0) bold_cells <- rbind(bold_cells,
                                               matrix(cbind(to_bold_p, which(names(out)=='p (adj)')),ncol=2))

  argL <- list(tab=out,to_indent=to_indent,bold_cells = bold_cells,
           caption=caption)
  if (!missing(fontsize)) argL[['fontsize']] <- fontsize
  do.call(outTable, argL)
}



# Survival Curves --------------------------------------------------------------


#' Plot KM and CIF curves with ggplot
#'
#' This function will plot a KM or CIF curve with option to add the number at
#' risk. You can specify if you want confidence bands, the hazard ratio, and
#' pvalues, as well as the units of time used.
#'
#' Note that for proper pdf output of special characters the following code
#' needs to be included in the first chunk of the rmd
#' knitr::opts_chunk$set(dev="cairo_pdf")
#'
#' @param response character vector with names of columns to use for response
#' @param cov String specifying the column name of stratification variable
#' @param data dataframe containing your data
#' @param type string indicating he type of univariate model to fit. The
#'   function will try and guess what type you want based on your response. If
#'   you want to override this you can manually specify the type. Options
#'   include "KM", and ,"CIF"
#'
#' @param pval boolean to specify if you want p-values in the plot (Log Rank
#'   test for KM and Gray's test for CIF)
#' @param HR boolean to specify if you want hazard ratios included in the plot
#' @param HR_pval boolean to specify if you want HR p-values in the plot
#' @param conf.curves boolean to specify if you want confidence interval bands
#' @param conf.type One of "none"(the default), "plain", "log" , "log-log" or
#'   "logit". Only enough of the string to uniquely identify it is necessary.
#'   The first option causes confidence intervals not to be generated. The
#'   second causes the standard intervals curve +- k *se(curve), where k is
#'   determined from conf.int. The log option calculates intervals based on the
#'   cumulative hazard or log(survival). The log-log option bases the intervals
#'   on the log hazard or log(-log(survival)), and the logit option on
#'   log(survival/(1-survival)).
#' @param table Logical value. If TRUE, includes the number at risk table
#' @param times Numeric vector of times for the x-axis
#' @param xlab String corresponding to xlabel. By default is "Time"
#' @param ylab String corresponding to ylabel. When NULL uses "Survival
#'   probability" for KM cuves, and "Probability of an event" for CIF
#' @param main String corresponding to main title. When NULL uses Kaplan-Meier
#'   Plot s, and "Cumulative Incidence Plot for CIF"
#'
#' @param stratalabs string corresponding to the labels of the covariate, when
#'   NULL will use the levels of the covariate
#' @param strataname String of the covariate name default is  nicename(cov)
#' @param stratalabs.table String corresponding to the levels of the covariate
#'   for the number at risk table, when NULL will use the levels of the
#'   covariate. Can use a string of "-" when the labels are long
#' @param strataname.table String of the covariate name for the number at risk
#'   table default is  nicename(cov
#'
#' @param median.text boolean to specify if you want the median values added to
#'   the legend (or as added text if there are no covariates), for KM only
#' @param median.lines boolean to specify if you want the median values added as
#'   lines to the plot, for KM only
#' @param median.CI boolean to specify if you want the 95\% confidence interval
#'   with the median text (Only for KM)
#' @param set.time.text string for the text to add survival at a specified time
#'   (eg. year OS)
#' @param set.time.line boolean to specify if you want the survival added as
#'   lines to the plot at a specified point
#' @param set.time Numeric values of the specific time of interest, default is 5
#'   (Multiple values can be entered)
#' @param set.time.CI boolean to specify if you want the 95\% confidence
#'   interval with the set time text
#'
#' @param censor.marks logical value. If TRUE, includes censor marks (only for
#'   KM curves)
#' @param censor.size size of censor marks, default is 3
#' @param censor.stroke stroke of censor marks, default is 1.5
#' @param fsize font size
#' @param nsize font size for numbers in the numbers at risk table
#' @param lsize line size
#' @param psize size of the pvalue
#' @param median.size size of the median text (Only when there are no
#'   covariates)
#' @param median.pos vector of length 2 corresponding to the median position
#'   (Only when there are no covariates)
#' @param median.lsize line size of the median lines
#' @param set.size size of the survival at a set time text (Only when there are
#'   no covariates)
#' @param set.pos  vector of length 2 corresponding to the survival at a set
#'   point position (Only when there are no covariates)
#' @param set.lsize line size of the survival at set points
#' @param ylim vector of length 2 corresponding to limits of y-axis. Default to
#'   NULL
#' @param col vector of colours
#' @param linetype vector of line types
#' @param xlim  vector of length 2 corresponding to limits of x-axis. Default to
#'   NULL
#' @param legend.pos Can be either a string corresponding to the legend position
#'   ("left","top", "right", "bottom", "none") or a vector of length 2
#'   corresponding to the legend position (uses normalized units (ie the
#'   c(0.5,0.5) is the middle of the plot))
#' @param pval.pos  vector of length 2 corresponding to the p-value position
#' @param plot.event  Which event(s) to plot (1,2, or c(1,2))
#' @param event String specifying if the event should be mapped to the colour,
#'   or linetype when plotting both events to colour = "col", line type
#' @param flip.CIF boolean to flip the CIF curve to start at 1
#' @param cut numeric value indicating where to divide a continuous covariate
#'   (default is the median)
#' @param eventlabs String corresponding to the event type names
#' @param event.name String corresponding to the label of the event types
#' @param Numbers_at_risk_text String for the label of the number at risk
#' @param HR.digits Number of digits printed of the  hazard ratio
#' @param HR.pval.digits Number of digits printed of the hazard ratio pvalue
#' @param pval.digits Number of digits printed of the Gray's/log rank pvalue
#'
#' @param median.digits Number of digits printed of the median pvalue
#' @param set.time.digits Number of digits printed of the probability at a
#'   specified time
#' @param print.n.missing Logical, should the number of missing be shown !Needs
#'   to be checked
#' @param returns Logical value returns a list with all ggplot objects in a list
#' @importFrom stats median qnorm as.formula pchisq model.matrix time
#' @importFrom cmprsk cuminc
#' @importFrom survival survfit survdiff
#' @importFrom ggplot2 ggplot
#' @importFrom gridExtra grid.arrange
#' @examples
#' # Simple plot without confidence intervals
#' ggkmcif(response = c('os_time','os_status'),
#' cov='cohort',
#' data=pembrolizumab)
#'
#' # Plot with median survival time
#' ggkmcif(response = c('os_time','os_status'),
#' cov='sex',
#' data=pembrolizumab,
#' median.text = TRUE,median.lines=TRUE,conf.curves=TRUE)
#'
#' # Plot with specified survival times and log-log CI
#' ggkmcif(response = c('os_time','os_status'),
#' cov='sex',
#' data=pembrolizumab,
#' median.text = FALSE,set.time.text = 'mo OS',
#' set.time = c(12,24),conf.type = 'log-log',conf.curves=TRUE)
#'
#' # KM plot with 95% CI and censor marks
#' ggkmcif(c('os_time','os_status'),'sex',data = pembrolizumab, type = 'KM',
#' HR=TRUE, HR_pval = TRUE, conf.curves = TRUE,conf.type='log-log',
#' set.time.CI = TRUE, censor.marks=TRUE)
#' @return Nothing is returned unless returns = TRUE is used. With returns =
#'   TRUE,  if table=TRUE (the default) a table style graphic with survival plot
#'   and number at risk table is returned. Otherwise a plot with the survival
#'   curves is returned.
#' @export
ggkmcif <- function(response,cov=NULL,data,type=NULL,
                    pval = TRUE,HR=FALSE,HR_pval=FALSE,conf.curves=FALSE,conf.type = "log",table = TRUE,
                    times = NULL,xlab = "Time",ylab=NULL ,
                    main = NULL,stratalabs = NULL,strataname = nicename(cov),
                    stratalabs.table=NULL,strataname.table=strataname,
                    median.text=FALSE,median.lines=FALSE,median.CI=FALSE,
                    set.time.text=NULL,set.time.line=FALSE,set.time=5,set.time.CI=FALSE,
                    censor.marks = TRUE,censor.size = 3,censor.stroke = 1.5,
                    fsize = 10, nsize = 3, lsize = 1, psize = 3.5,
                    median.size=3,median.pos=NULL,median.lsize=1,
                    set.size=3,set.pos=NULL,set.lsize=1,
                    ylim=c(0,1), col=NULL,linetype=NULL, xlim=NULL,
                    legend.pos = NULL,  pval.pos=NULL,plot.event=1,event=c("col","linetype"),flip.CIF =FALSE,
                    cut=NULL,eventlabs=NULL,event.name=NULL,Numbers_at_risk_text="Numbers at risk",
                    HR.digits = 2,HR.pval.digits=3, pval.digits=3,
                    median.digits=3,set.time.digits=3,returns = FALSE,print.n.missing=TRUE){
  event <- match.arg(event)

  # to enable the function to work with tibbles
  data <- data.frame(data)

  ##Removing all missing data
  #rowSums(is.na(final[ , 5:6])) == 0,
  if(!is.null(cov)) remove <- rowSums(is.na(data[ ,c(response,cov)])) > 0
  if(is.null(cov)) remove <- rowSums(is.na(data[ ,c(response)])) > 0

  data <- data[!remove,]

  if (print.n.missing==T & sum(remove)==1 ) {
    message('1 observation with missing data was removed.')
  } else if(print.n.missing==T & sum(remove) > 0) message(paste(sum(remove),'observations with missing data were removed.'))


  if (!is.factor(data[,cov])  & !is.numeric(data[,cov]) & !is.null(cov)){ message("Coercing the cov variable to factor"); data[,cov] <- factor(data[,cov])}

  if (is.numeric(data[,cov])){
    numeric = T
    if(is.null(cut)) cut <- median(data[,cov],na.rm = T)
    data[,cov] <- factor(ifelse(data[,cov]<=cut,paste0("<=",round(cut,2)),paste0(">",round(cut,2))),levels = c(paste0("<=",round(cut,2)),paste0(">",round(cut,2))))
    out_fmt = ifelse(is.null(knitr::pandoc_to()),'html',
                     ifelse(knitr::pandoc_to(c('doc','docx')),'doc',
                            ifelse(knitr::is_latex_output(),'latex','html')))
    #  le_code <- ifelse(out_fmt=='latex',"\\le","\u2264")
    le_code <- "\u2264"
    if(is.null(stratalabs)) stratalabs = c(paste0(le_code,round(cut,2)),paste0(">",round(cut,2)))
  }

  data[,cov] <- droplevels(data[,cov])
  if(is.null(stratalabs)) stratalabs <- levels(data[,cov])
  if(is.null(stratalabs.table)) stratalabs.table <- stratalabs

  if(!is.null(col) & !is.null(cov) & (event != "col" | length(plot.event)==1)) col <- rep(col,length.out=length(levels(data[,cov])))
  if(!is.null(linetype) & !is.null(cov) & (event != "linetype" | length(plot.event)==1)) linetype <- rep(linetype,length.out=length(levels(data[,cov])))


  if(is.null(col)){
    if(length(plot.event)==1|event != "col"){
      col_length <- ifelse(is.null(cov),1,length(levels(data[,cov])))
    }else col_length = 2

    col <- color_palette_surv_ggplot(col_length)
  }

  # Specifying the type of plot ----------------------------------------------


  #Specifying KM or CIF & is.null(type)
  if (length(unique(data[, response[2]]))< 3 & is.null(type)) type = "KM"
  if (length(unique(data[, response[2]]))>= 3 & is.null(type)) type = "CIF"

  ##REMOVING MEDIAN FOR CIF
  if(type=="CIF") {
    median.lines = FALSE
    median.text=FALSE
  }


  if(type=="KM") {
    if(is.null(main)) main <- "Kaplan-Meier Plot"
    if(is.null(ylab)) ylab = "Survival Probability"
  }else if(type=="CIF"){
    #For the number at risk calcuation
    #data_sfit[,response[2]][data_sfit[,response[2]]!=0] <- 1
    if(is.null(main)) main <- "Cumulative Incidence Plot"
    if(is.null(ylab)) ylab = "Probability of an Event"
  }else(stop("Type must be either KM or CIF"))


  if(flip.CIF==T & type=="CIF") {
    median.text=F
    median.lines=F
    set.time.text=NULL
    set.time.line=F
    set.time=NULL
  }

  # Labels ------------------------------------------------------------------
  multiple_lines <- !is.null(cov)


  # HR and p-val cox----------------------------------------------------------------------
  if(type=="KM" & multiple_lines & (HR|HR_pval)){
    coxfit <- survival::coxph(as.formula(paste(paste("survival::Surv(", response[1],
                                                     ",", response[2], ")", sep = ""), "~", cov,
                                               sep = "")), data = data)


    HR_vals <- paste0("HR=",sapply(seq(length(stratalabs)-1),function(i){
      return(psthr0(summary(coxfit)$conf.int[i,c(1, 3, 4)],digits=HR.digits))
    }))

    if(HR)stratalabs[-1] <- paste(stratalabs[-1],HR_vals)
    if(HR_pval) stratalabs[-1] <- paste(stratalabs[-1],sapply(summary(coxfit)$coef[,5],lpvalue2,digits=HR.pval.digits))
    stratalabs[1] <- paste(stratalabs[1],"REF")
  }


  # HR and p-val crr --------------------------------------------------------
  if(type=="CIF" & multiple_lines & (HR|HR_pval)&length(plot.event==1) & plot.event[1]==1){
    crrfit <- crrRx(as.formula(paste(paste(response,
                                           collapse = "+"), "~", cov, sep = "")),
                    data = data)

    HR_vals <- paste0("HR=",sapply(seq(length(stratalabs)-1),function(i){
      return(psthr0(summary(crrfit)$conf.int[i,c(1, 3, 4)],digits=HR.digits))
    }))

    if(HR)stratalabs[-1] <- paste(stratalabs[-1],HR_vals)
    if(HR_pval) stratalabs[-1] <- paste(stratalabs[-1],sapply(summary(crrfit)$coef[,5],lpvalue2,digits=HR.pval.digits))
    stratalabs[1] <- paste(stratalabs[1],"REF")

  }
  # Model fitting KM and creating a dataframe of times--------------------------------------------------------
  if(type=="KM"){
    if(!multiple_lines){

      sfit <- survival::survfit(as.formula(paste(paste("survival::Surv(", response[1],
                                                       ",", response[2], ")", sep = ""), "~", 1,
                                                 sep = "")), data = data,conf.type=conf.type)

      if(median.lines==T|median.text==TRUE){
        median_vals <- summary(sfit)$table['median']
        median_lower <- summary(sfit)$table['0.95LCL']
        median_upper <- summary(sfit)$table['0.95UCL']

        median_txt <- if(median.CI==TRUE) {paste0(round_sprintf(median_vals,digits=median.digits),"(",round_sprintf(median_lower,digits=median.digits),"-",
                                                  round_sprintf(median_upper,digits=median.digits),")")
        }else round_sprintf(median_vals,digits=median.digits)
      }
      if(!is.null(set.time.text)|set.time.line==TRUE){

        set.surv.text = NULL
        set.surv = NULL
        for(time_i in set.time)
        {
          sum <- summary(sfit,time=c(0,time_i)) ##The zero is to make sure all levels are included
          df_text <- data.frame(time=sum$time)
          df_text$set.CI<- if(set.time.CI==TRUE) {paste0(round_sprintf(sum$surv,digits=set.time.digits),"(",round_sprintf(sum$lower,digits=set.time.digits),"-",
                                                         round_sprintf(sum$upper,digits=set.time.digits),")")
          }else round_sprintf(sum$surv,digits=set.time.digits)
          keep_sum <- df_text
          if(length(sum$surv)>1) keep_sum <- df_text[sum$time!=0,]
          keep_sum$set.CI[keep_sum$time==0] <- NA
          set.surv <- rbind(set.surv,keep_sum)

          if(is.null(set.surv.text)) {set.surv.text <- paste0(time_i," ",set.time.text,"=",keep_sum$set.CI)
          }else  set.surv.text <- paste0(set.surv.text,",",time_i," ",set.time.text,"=",keep_sum$set.CI)

        }

      }
      stratalabs <- "All"

    }else{
      sfit <- survival::survfit(as.formula(paste(paste("survival::Surv(", response[1],
                                                       ",", response[2], ")", sep = ""), "~", cov,
                                                 sep = "")), data = data,conf.type=conf.type)

      if(median.lines==T|median.text==TRUE) {
        median_vals <- summary(sfit)$table[,'median']
        median_lower <- summary(sfit)$table[,'0.95LCL']
        median_upper <- summary(sfit)$table[,'0.95UCL']

        med_txt <- if(median.CI==TRUE) {paste0(round_sprintf(median_vals,digits=median.digits),"(",round_sprintf(median_lower,digits=median.digits),"-",
                                               round_sprintf(median_upper,digits=median.digits),")")
        }else round_sprintf(median_vals,digits=median.digits)

        if(median.text==TRUE) stratalabs <- paste(stratalabs,", Median=",med_txt)
      }

      if(!is.null(set.time.text)|set.time.line==TRUE) {

        final.set.text = NULL
        set.surv = NULL
        for(time_i in set.time)
        {
          sum <- summary(sfit,time=c(0,time_i)) ##The zero is to make sure all levels are included
          df_text <- data.frame(strat=sum$strata,time=sum$time)
          df_text$set.CI<- if(set.time.CI==TRUE) {paste0(round_sprintf(sum$surv,digits=set.time.digits),"(",round_sprintf(sum$lower,digits=set.time.digits),"-",
                                                         round_sprintf(sum$upper,digits=set.time.digits),")")
          }else round_sprintf(sum$surv,digits=set.time.digits)
          dup_strata = sum$strata[duplicated(sum$strata)]
          keep_sum <- df_text[!(sum$strata %in% dup_strata) | sum$time!=0,]
          keep_sum$set.CI[keep_sum$time==0] <- NA

          set.surv <- rbind(set.surv,keep_sum)


          if(is.null(final.set.text)) {final.set.text <- paste0(time_i," ",set.time.text,"=",keep_sum$set.CI)
          }else  final.set.text <- paste0(final.set.text,",",time_i," ",set.time.text,"=",keep_sum$set.CI)

        }


        if(!is.null(set.time.text)) stratalabs <- paste0(stratalabs,", ",final.set.text)
      }
    }
    df <- NULL
    df <- data.frame(time = sfit$time,
                     n.risk = sfit$n.risk,  n.censor = sfit$n.censor,
                     n.event = sfit$n.event, surv = sfit$surv,
                     strata = if(multiple_lines){
                       summary(sfit, censored = T)$strata
                     }else factor("All"),
                     upper = if(conf.type != "none" & conf.curves==TRUE){
                       sfit$upper
                     }else factor(NA),
                     lower = if(conf.type != "none"& conf.curves==TRUE){
                       sfit$lower
                     }else factor(NA))
    levels(df$strata) <- stratalabs
    zeros <- data.frame(time = 0, surv = 1,
                        strata = if(multiple_lines){
                          levels(df$strata)
                        }else factor("All"),
                        upper = 1, lower = 1)
    df <- plyr::rbind.fill(zeros, df) # Forcing the curves to start at 1

    df$strata <- factor(df$strata,levels=stratalabs)
  }


  # Model fitting CIF -------------------------------------------------------
  if(type=="CIF"){

    ### Functions for getting the median values
    median_time_to_event <- function(name){
      estall=get(name,fit)$est
      timall=get(name,fit)$time
      return(timall[estall>=0.5][1])

    }

    if(!multiple_lines){

      invisible(utils::capture.output(fit <-  cuminc( data[,response[1]],data[,response[2]] )))
      stratalabs <- " "
      gsep = " "

      last_character <- substr(names(fit), nchar(names(fit)), nchar(names(fit)))
      get_values <- names(fit)[last_character %in% plot.event]
      if(median.lines==T|median.text==TRUE) {
        median_vals = sapply(get_values, median_time_to_event)
        median_txt <- round_sprintf(median_vals,digits=median.digits)

      }
      if(!is.null(set.time.text)|set.time.line==TRUE) {

        set.surv.text = NULL
        set.surv = NULL
        for(time_i in set.time)
        {
          z <- qnorm(1-(1-0.95)/2)
          val <- cmprsk::timepoints(fit,times=time_i)$est[rownames(cmprsk::timepoints(fit,times=time_i)$est) %in% get_values,]
          var <- cmprsk::timepoints(fit,times=time_i)$var[rownames(cmprsk::timepoints(fit,times=time_i)$est) %in% get_values,]

          lower <- val ^ exp(-z*sqrt(var)/(val*log(val)))
          upper <- val ^ exp(z*sqrt(var)/(val*log(val)))

          keep_sum <- data.frame(val,time_i)
          names(keep_sum)[1] <- 'value'
          keep_sum$set.CI<- if(set.time.CI==TRUE) {paste0(round_sprintf(val,digits=set.time.digits),"(",round_sprintf(lower,digits=set.time.digits),"-",
                                                          round_sprintf(upper,digits=set.time.digits),")")
          }else round_sprintf(val,digits=set.time.digits)

          if(is.null(set.surv.text)) {set.surv.text <- paste0(time_i," ",set.time.text,"=",keep_sum$set.CI)
          }else  set.surv.text <- paste0(set.surv.text,",",time_i," ",set.time.text,"=",keep_sum$set.CI)


          set.surv <- rbind(keep_sum,set.surv)

        }
      }

      if(table){ #Sfit is for the numbers at risk so both events are counted the same way

        temp <- data
        temp[,response[2]][temp[,response[2]] > 0] <- 1
        sfit <- survival::survfit(as.formula(paste(paste("Surv(", response[1],
                                                         ",", response[2], ")", sep = ""), "~", 1,
                                                   sep = "")), data = temp)



      }

    }else{
      newgpvar <- paste0(data[,cov],":")
      newgpvar <- factor(newgpvar, levels = paste0(levels(data[,cov]),":") )
      invisible(utils::capture.output(fit <- cuminc(data[,response[1]],data[,response[2]], newgpvar)))
      gsep = ": "


      last_character <- substr(names(fit), nchar(names(fit)), nchar(names(fit)))
      get_values <- names(fit)[last_character %in% plot.event]

      if(median.lines==T|median.text==TRUE) median_vals <- sapply(get_values, median_time_to_event)
      if(median.text==T & length(plot.event)==1) stratalabs <- paste(stratalabs,", Median=",round_sprintf(median_vals,digits=median.digits))

      if(!is.null(set.time.text)|set.time.line==TRUE) {
        set.surv.text = NULL
        set.surv = NULL
        for(time_i in set.time)
        {
          z <- qnorm(1-(1-0.95)/2)
          val <- cmprsk::timepoints(fit,times=time_i)$est[rownames(cmprsk::timepoints(fit,times=time_i)$est) %in% get_values,]
          var <- cmprsk::timepoints(fit,times=time_i)$var[rownames(cmprsk::timepoints(fit,times=time_i)$est) %in% get_values,]

          lower <- val ^ exp(-z*sqrt(var)/(val*log(val)))
          upper <- val ^ exp(z*sqrt(var)/(val*log(val)))

          keep_sum <- data.frame(val,time_i)
          names(keep_sum)[1] <- 'value'
          keep_sum$set.CI<- if(set.time.CI==TRUE) {paste0(round_sprintf(val,digits=set.time.digits),"(",round_sprintf(lower,digits=set.time.digits),"-",
                                                          round_sprintf(upper,digits=set.time.digits),")")
          }else round_sprintf(val,digits=set.time.digits)

          if(is.null(set.surv.text)) {set.surv.text <- paste0(time_i," ",set.time.text,"=",keep_sum$set.CI)
          }else  set.surv.text <- paste0(set.surv.text,",",time_i," ",set.time.text,"=",keep_sum$set.CI)


          set.surv <- rbind(keep_sum,set.surv)

        }


      }
      if(!is.null(set.time.text)& length(plot.event)==1) stratalabs <- paste0(stratalabs,", ",set.surv.text)

      if(table){ #Sfit is for the numbers at risk so both events are counted the same way

        temp <- data
        temp[,response[2]][temp[,response[2]] > 0] <- 1
        sfit <- survival::survfit(as.formula(paste(paste("Surv(", response[1],
                                                         ",", response[2], ")", sep = ""), "~", cov,
                                                   sep = "")), data = temp)

      }

    }

    ##Setting up the data frame
    if (!is.null(fit$Tests)){
      test <- fit$Test
      fit <- fit[names(fit) != "Tests"]
    }
    fit2 <- lapply(fit, `[`, 1:3)
    gnames <- names(fit2)

    fit2_list <- lapply(seq_along(gnames), function(ind) {
      df <- as.data.frame(fit2[[ind]])
      df$name <- gnames[ind]
      df
    })

    df <- do.call(rbind, fit2_list)
    df$event <- sapply(strsplit(df$name, split = gsep), `[`, 2)
    df$strata <- sapply(strsplit(df$name, split = gsep), `[`, 1)
    df$strata <- factor(df$strata, levels = levels(data[,cov]) )
    levels(df$strata) <- stratalabs

    if(multiple_lines){df$strata <- factor(df$strata,levels=stratalabs)
    }else df$strata <- "ALL"

    df$std <- std <- sqrt(df$var)

    names(df)[names(df)=="est"] <- 'surv'

    df <- df[df$event %in% plot.event,]
    df$upper <- NA
    df$lower <- NA


    if(conf.type!="log" & conf.type!="none") message("Only log confidence intervals avaliable for CIF")


    if(conf.type=="log") {
      z <- qnorm(1-(1-0.95)/2)
      df$lower <- df$surv ^ exp(-z*sqrt(df$var)/(df$surv*log(df$surv)))
      df$upper <- df$surv ^ exp(z*sqrt(df$var)/(df$surv*log(df$surv)))
    }

    if(flip.CIF) {
      df$surv <- 1-df$surv
      df$upper <- 1-df$upper
      df$lower <- 1-df$lower
    }

    if(!is.null(eventlabs)) {
      df$event <- factor(df$event)
      levels(df$event) <- eventlabs
    }else eventlabs <- c(1,2)



  }


  # axis and legend positions --------------------------------------------------------------------
  m <- max(nchar(stratalabs))
  maxxval = max(df$time,times[length(times)])
  if( is.null(xlim) ){
    maxxlim =  maxxval
  }else {
    if( length(xlim)!=2 | !is.numeric(xlim) ) stop("xlim must be a numeric vector of length 2.")
    maxxlim = xlim[2]
  }


  linetype_name <- col_name <- strataname

  leg.pos <- legend.pos
  d <- length(levels(df$strata))

  if(!multiple_lines & (type=="KM" | (type=="CIF" & length(plot.event)==1))){
    leg.pos <- 'none'
  }else if(is.null(legend.pos)) {
    if(type=="CIF" & flip.CIF==F){ leg.pos <- c(min(0.1+m/200,0.5), 0.9-d*0.05)
    }else leg.pos <- c(min(0.1+m/200,0.5), 0.05+d*0.05)
  }


  if(is.null(times)) times <- break_function(maxxlim)

  # plotting ----------------------------------------------------------------

  if(type=="CIF" & length(plot.event)>1){
    if(is.null(event.name)) event.name <- 'event'
    if(event=='linetype') {p <- ggplot(df)+ geom_step( aes(time, surv, color = strata,linetype=event),size = lsize); linetype_name = event.name}
    if(event=='col') {p <- ggplot(df)+ geom_step( aes(time, surv, color = event,linetype=strata),size = lsize); col_name=event.name}
  }else{
    p <- ggplot(df) + geom_step(aes(time, surv, group = strata,linetype = strata, col=strata), size = lsize)
  }


  # Confidence intervals to the plot ----------------------------------------


  if(conf.type!="none" & conf.curves==TRUE){
    if(type=="KM")p <- p +  geom_ribbon(data=df[!is.na(df$upper) & !is.na(df$lower),], aes(x=time, fill=strata, ymin = lower, ymax = upper), inherit.aes = FALSE, alpha = 0.2, show.legend = F)

    if(type=="CIF" & (event == "linetype"|length(plot.event)==1)){
      for( evnt in unique(df$event) ){
        p <- p + geom_ribbon(data=df[!is.na(df$upper) & !is.na(df$lower) & df$event == evnt,], aes(x=time, fill=strata, ymin = lower, ymax = upper), inherit.aes = FALSE, alpha = 0.2, show.legend = F)
      }
    }

    if(type=="CIF" & event == "col" & length(plot.event)>1){
      for( stra in unique(df$strata) ){
        p <- p + geom_ribbon(data=df[!is.na(df$upper) & !is.na(df$lower) & df$strata == stra,], aes(x=time, fill=event, ymin = lower, ymax = upper), inherit.aes = FALSE, alpha = 0.2, show.legend = F)
      }
    }

  }

  # Modifying axis, titles, fonts, and themeing -------------------------------


  p <- p+
    theme_classic(base_size=fsize) +
    theme(
      axis.text.x = element_text(margin = margin(t = 0), vjust = 1),
      axis.text.x.top = element_text(margin = margin(b = 0), vjust = 0),
      axis.text.y = element_text(margin = margin(r = 0), hjust = 1),
      axis.text.y.right = element_text(margin = margin(l = 0), hjust = 0),
      axis.title.x.bottom = element_text(vjust = 4))+
    scale_x_continuous(paste0("\n",xlab), breaks = times,
                       limits = c(0, maxxval)) +
    coord_cartesian(xlim=c(0,maxxlim)) + ### changes the actual plotted limits if needed
    scale_y_continuous(paste0(ylab,"\n"), limits = ylim) +
    theme(panel.grid.minor = element_blank()) +
    theme(legend.key = element_rect(colour = "transparent", fill = "transparent"),
          legend.background=element_blank(),
          legend.position = leg.pos) +
    labs(linetype = linetype_name, col=col_name) +
    ggtitle(main) +
    theme(plot.title = element_text(face="bold",hjust = 0.5,size = fsize))


  # censoring ---------------------------------------------------------------


  if( censor.marks & type=="KM") p <- p + geom_point(data = subset(df, n.censor>0),
                                                     aes(x=time, y=surv, group = strata, col=strata),
                                                     shape=3, size=censor.size, stroke = censor.stroke,show.legend =F)


  # Log rank p-val ----------------------------------------------------------
  if(pval & type=="KM" & multiple_lines) {
    sdiff <- survival::survdiff(eval(sfit$call$formula), data = eval(sfit$call$data))
    pval <- pchisq(sdiff$chisq, length(sdiff$n)-1, lower.tail = FALSE)
    pvaltxt <- lpvalue2(pval,pval.digits)
    pvaltxt <- paste(pvaltxt,"(Log Rank)")
    #pvaltxt <- paste("p =", signif(pval, 3), "(Log Rank)")
    if(is.null(pval.pos)){
      p <- p + annotate("text", x = 0.85 * max(times), y = ylim[1], label = pvaltxt, size = psize)
    }else p <- p + annotate("text", x = pval.pos[1], y = pval.pos[2], label = pvaltxt, size = psize)
  }


  # Gray's pvalue -----------------------------------------------------------

  if( !is.null(cov) & pval & type=="CIF") {
    if(length(plot.event)==1 ){
      test <- test[rownames(test)==plot.event,]
      pval <- test[2]
      pvaltxt <- lpvalue2(pval,pval.digits)
      pvaltxt <- paste(pvaltxt,"(Gray's test)")
      if(is.null(pval.pos)){
        p <- p + annotate("text", x = 0.85 * max(times), y = ylim[1], label = pvaltxt, size = psize)
      }else p <- p + annotate("text", x = pval.pos[1], y = c(pval.pos[2]), label = pvaltxt, size = psize)
    }else{



      pval <- test[,2]
      pvaltxt <- sapply(pval,lpvalue2,pval.digits)
      pvaltxt <- c("Gray's test",paste(eventlabs, pvaltxt))
      if(is.null(pval.pos)){

        p <- p + annotate("text", x = 0.85 * max(df$time), y = c(0.12,0.08,0.04), label = pvaltxt, size = psize)
      }else p <- p + annotate("text", x = pval.pos[1], y = c(pval.pos[2],pval.pos[2]-0.04, pval.pos[2]-0.08), label = pvaltxt, size = psize)
    }
  }


  # median and set time with one level text ---------------------------------------------------

  if(median.text &!multiple_lines) {
    median_txt <- paste0("Median=",median_txt)
    if(length(plot.event)==2) median_txt <- paste(paste0(eventlabs,":",median_txt),collapse="\n")


    if(is.null(median.pos) & type=="KM")  median.pos <- c( 0.1 * max(times),ylim[1])
    if(is.null(median.pos) & type=="CIF")  median.pos <- c( 0.1 * max(times),ylim[2]*.95)

    p <- p + annotate("text", x =median.pos[1], y = median.pos[2], label = median_txt, size = median.size)

  }

  if(!is.null(set.time.text) &!multiple_lines) {
    #set.surv.text <- paste0(set.time,"",set.time.text,"=",round_sprintf(set.surv,digits=set.time.digits))
    if(length(plot.event)==2) set.surv.text <- paste(paste0(eventlabs,":",set.surv.text),collapse="\n")

    if(is.null(set.pos) & type=="KM")  set.pos <- c( 0.1 * max(times),ylim[1]+0.1)
    if(is.null(set.pos) & type=="CIF")  set.pos <- c( 0.1 * max(times),ylim[2]*.85)


    p <- p + annotate("text", x =set.pos[1], y = set.pos[2], label = set.surv.text, size = set.size)

  }


  # median and set time lines ------------------------------------------------------------

  if(median.lines){
    temp <- data.frame(x=median_vals,y=0.5)
    p <- p + geom_segment(data=temp,aes(x=x,xend=x,y=0,yend=0.5),lty=2,lwd=median.lsize)
    temp2 <- data.frame(x=max(median_vals,na.rm=TRUE),y=0.5)
    p <- p + geom_segment(data=temp2,aes(x=0,xend=x,y=0.5,yend=0.5),lty=2,lwd=median.lsize)
  }

  if(set.time.line){

    set.surv$y <- as.numeric(sub(" *\\(.*", "", set.surv$set.CI))
    set.surv$time[is.na(set.surv$y)] <- NA

    p <- p + geom_segment(data=set.surv,aes(x=0,xend=time,y=y,yend=y),lty=2,lwd=set.lsize)
    p <- p + geom_segment(data=set.surv,aes(x=time,xend=time,y=0,yend=y),lty=2,lwd=set.lsize)
  }


  # Colour, linetype and fill -----------------------------------------------


  ##Getting the colour labels
  if(multiple_lines==T & (length(plot.event)==1|event=="linetype")){
    col_labs <- stratalabs
  }else col_labs <- eventlabs

  p <- p + scale_colour_manual(values = col,labels=col_labs)
  p <- p + scale_fill_manual(values = col,labels=col_labs)

  ##Getting the linetype labels
  if(multiple_lines==T & (length(plot.event)==1|event=="col")){
    linetype_labs <- stratalabs
  }else linetype_labs <- eventlabs


  if(!is.null(linetype)) p <- p + scale_linetype_manual(labels=linetype_labs,values = linetype)
  if(is.null(linetype)) p <- p + scale_linetype_discrete(labels=linetype_labs)


  if(event == "linetype" & length(plot.event)==2 & multiple_lines==F){

    p <- p + guides(col = F)
  }

  if(event == "col" & length(plot.event)==2 & multiple_lines==F){

    p <- p + guides(linetype = F)
  }


  # Numbers at risk ---------------------------------------------------------

  if(table) {
    ## Create table graphic to include at-risk numbers

    blank.pic <- ggplot(df, aes(time, surv)) +
      geom_blank() +
      theme_bw() +
      scale_x_continuous(breaks = times, limits = c(0, maxxval)) +
      theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
            axis.title.x = element_blank(), axis.title.y = element_blank(),
            axis.ticks = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(), panel.background = element_rect(fill = "transparent"))

    sfit.summary <- summary(sfit, times = times, extend = TRUE)
    risk.data <- data.frame(strata = if(multiple_lines){
      sfit.summary$strata
    }else factor("All"),
    time = sfit.summary$time,
    n.risk = sfit.summary$n.risk)
    # if risk and event do paste0(sfit.summary$n.risk, "(",sfit.summary$n.events,")")
    risk.data$strata <- factor(risk.data$strata, levels=rev(levels(risk.data$strata)))

    if(!is.null(col)){cols1 <- col
    }else{cols1 <- .extract_ggplot_colors(p, grp.levels = stratalabs)}



    if(multiple_lines==F|(event == "col" & length(plot.event)>1)) cols1 = "black"


    n_strata <- length(stratalabs)
    #yticklabs <- rep("-", n_strata)


    yticklabs <- unname(rev(stratalabs.table))
    #strataname.table = ifelse(n_strata==1, "\n", strataname)

    data.table <- ggplot(risk.data, aes(x = time, y = strata, label = format(n.risk, nsmall = 0))) +
      geom_text(hjust="middle", vjust="center", size = nsize) +
      theme_bw() +
      scale_x_continuous(Numbers_at_risk_text, breaks = times, limits = c(0,maxxval)) +
      coord_cartesian(xlim=c(0,maxxlim)) + ### changes the actual plotted limits if needed
      theme(legend.position = "none") +
      theme(text = element_text(size = fsize)) + ##increase font size
      theme(plot.margin = unit(c(-1.5, 1, 0.1, 0.2), "lines"))+
      scale_y_discrete(strataname.table, breaks = as.character(levels(risk.data$strata)), labels = yticklabs)

    data.table <- data.table +suppressWarnings(theme(axis.title.x = element_text(size = fsize, vjust = 1), panel.grid.major = element_blank(),
                                                     panel.grid.minor = element_blank(), panel.border = element_blank(),
                                                     axis.text.x = element_blank(), axis.ticks = element_blank(),
                                                     axis.text.y = element_text(face = "bold", hjust = 1,colour = rev(cols1))))




    if(all(as.character(yticklabs)=="-")) data.table <- .set_large_dash_as_ytext(data.table)



    gA <- ggplotGrob(p)
    gB <- ggplotGrob(blank.pic)
    gC <- ggplotGrob(data.table)
    maxWidth = grid::unit.pmax(gA$widths[2:5], gC$widths[2:5])
    gA$widths[2:5] <- as.list(maxWidth)
    gB$widths[2:5] <- as.list(maxWidth)
    gC$widths[2:5] <- as.list(maxWidth)

    gridExtra::grid.arrange(gA, gB, gC,
                            clip = FALSE, nrow = 3, ncol = 1,
                            heights = unit(c(2, .1, .25), c("null", "null", "null")))

    if(returns) {
      if(table==TRUE){
        a <- list(p,blank.pic,data.table)
      } else a <- p

      return(a)
    }
  } else(return(p))


}

modify_ggkmcif <- function(list_gg){
  .Deprecated("ggkmcif_paste")
  ggkmcif_paste(list_gg)
}

#' Plot KM and CIF curves with ggplot
#'
#' This function puts together a survival curve, and a number at risk table
#'
#' @param list_gg list containing the results of ggkmcif
#' @return a gtable with three elements, the survival curve, a spacer and the
#'   number at risk table
#' @export
#' @examples
#' plot <- ggkmcif(response=c('pfs_time','pfs_status'),
#' data=pembrolizumab,returns = TRUE)
#'
#' # Highlighting a section of the curve
#' plot[[1]] <- plot[[1]] +
#' ggplot2::geom_rect(xmin=4,xmax=8,ymin=0.15,ymax=0.4,alpha=0.01,fill='yellow')
#'
#' # Putting the curve back together
#' ggkmcif_paste(plot)
ggkmcif_paste <- function(list_gg){
  gA <- ggplotGrob(list_gg[[1]])
  gB <- ggplotGrob(list_gg[[2]])
  gC <- ggplotGrob(list_gg[[3]])
  maxWidth = grid::unit.pmax(gA$widths[2:5], gC$widths[2:5])
  gA$widths[2:5] <- as.list(maxWidth)
  gB$widths[2:5] <- as.list(maxWidth)
  gC$widths[2:5] <- as.list(maxWidth)

  gridExtra::grid.arrange(gA, gB, gC,
                          clip = FALSE, nrow = 3, ncol = 1,
                          heights = unit(c(2, .1, .25), c("null", "null", "null")))
}

# Survival Summaries --------------------------------------------------------------

#' Display event counts, expected event counts and logrank test of differences
#'
#' This is a wrapper function around the survdiff function to display overall
#' event rates and group-specific rates along with the log-rank test of a
#' difference in survival between groups in a single table suitable for markdown
#' output. Median survival times are included by default but can be removed
#' setting median=FALSE
#' @param data data frame containing survival data
#' @param time string indicating survival time variable
#' @param status string indicating event status variable
#' @param covs character vector indicating variables to group observations by
#' @param strata string indicating the variable to stratify observations by
#' @param includeVarNames boolean indicating if the variable names should be
#'   included in the output table, default is FALSE
#' @param digits the number of digits in the survival rate
#' @param showCols character vector indicating which of the optional columns to
#'   display, defaults to c('N','Observed','Expected')
#' @param CIwidth width of the median survival estimates, default is 95%
#' @param conf.type type of confidence interval see \code{\link{survfit}} for
#'   details. Default is 'log'.
#' @param caption table caption
#' @param tableOnly should a dataframe or a formatted object be returned
#' @param fontsize PDF/HTML output only, manually set the table fontsize
#' @importFrom  survival survdiff Surv strata
#' @seealso \code{\link{survdiff}}
#' @examples
#' #' # Differences between sex
#' rm_survdiff(data=pembrolizumab,time='os_time',status='os_status',
#' covs='sex',digits=1)
#'
#' # Differences between sex, stratified by cohort
#' rm_survdiff(data=pembrolizumab,time='os_time',status='os_status',
#' covs='sex',strata='cohort',digits=1)

#' # Differences between sex/cohort groups
#' rm_survdiff(data=pembrolizumab,time='os_time',status='os_status',
#' covs=c('sex','cohort'),digits=1)
#' @return A character vector of the survival table source code, unless tableOnly=TRUE in
#'   which case a data frame is returned
#' @export
rm_survdiff <- function(data,time,status,covs,strata,includeVarNames=FALSE,
                        digits=1,showCols=c('N','Observed','Expected'),CIwidth=0.95,
                        conf.type='log',caption=NULL,tableOnly=FALSE,fontsize){
  if (missing(data)) stop('data is a required argument')
  if (missing(time)) stop('time is a required argument')
  if (missing(status)) stop('status is a required argument')
  if (missing(covs)) stop('covs is a required argument')
  if  (!inherits(data,'data.frame')) stop('data must be supplied as a data frame')
  if( !inherits(time,'character' )| length(time)>1)
    stop('time must be supplied as a string indicating a variable in data')
  if (!inherits(status,'character' )| length(status)>1)
    stop('status must be supplied as a string indicating a variable in data')
  if (!inherits(covs,'character' ))
    stop('covs must be supplied as a character vector or string indicating variables in data')
  missing_vars = na.omit(setdiff(c(time, status,covs), names(data)))
  if (length(missing_vars) > 0) stop(paste("These variables are not in the data:\n",
                                           paste0(missing_vars,collapse=csep())))
  if (missing(strata)){
    s_f <- ''
    lr_txt <- ''
  } else {
    if (!inherits(strata,'character' ) | length(strata)>1 | !strata %in% names(data)) stop('strata must be supplied as a string indicating a variable in data')
    if (length(intersect(covs,strata))>0) stop('A variable can appear in covs or strata but not both.')
    s_f <- paste0('+strata(',strata,')')
    lr_txt <-paste('stratified by',strata)
  }
  lr_test = survival::survdiff(as.formula(paste0("survival::Surv(",time,',',status,') ~',
                                                 paste0(covs,collapse = '+'),s_f)),
                               data = data)

  gfit <- survival::survfit(as.formula(paste0("survival::Surv(",time,',',status,') ~',
                                              paste0(covs,collapse = '+'))),
                            data = data)
  gtbl <- t(summary(gfit)$table)
  ofit <- survival::survfit(as.formula(paste0("survival::Surv(",time,',',status,') ~1')),
                            data = data,conf.type=conf.type,conf.int=CIwidth)
  otbl <- summary(ofit)$table
  otbl <- data.frame(Overall=otbl)
  mtbl <- data.frame(t(cbind(otbl,gtbl)))
  m_CI <- apply(mtbl[,grep('median|LCL|UCL',names(mtbl))],1,function(x) psthr(x,y=digits))
  m_CI_nm <-paste0('Median (',CIwidth*100,'%CI)')
  if (inherits(lr_test$obs,'numeric')) {
    ob <- lr_test$obs
    ex <- lr_test$exp
    df <- length(lr_test$obs)-1
  } else {
    ob <- rowSums(lr_test$obs)
    ex <- rowSums(lr_test$exp)
    df <- nrow(lr_test$obs)-1
  }
  lr_data <-data.frame(group=names(lr_test$n),
                       N=as.numeric(lr_test$n),
                       Observed=ob,
                       Expected=niceNum(ex,1))
  lr_data <- rbind(c(group='Overall',
                     N=otbl['records','Overall'],
                     Observed = otbl['events','Overall'],
                     Expected=NA),
                   lr_data)
  lr_data <- cbind(lr_data,m_CI)
  names(lr_data) <- gsub('m_CI',m_CI_nm,names(lr_data))
  rownames(lr_data) <- NULL
  lr_data <- rbind(lr_data,c('Log Rank Test',NA,NA,NA,
                             paste('ChiSq =',niceNum(lr_test$chisq,1),'on',df,'df')))
  lr_data <- rbind(lr_data,c(lr_txt,NA,NA,NA,paste('p-value =',formatp(pchisq(lr_test$chisq,df,lower.tail = FALSE)))))
  if (!includeVarNames){
    for (v in covs) lr_data$group <- gsub(paste0(v,'='),'',lr_data$group)
  }
  lr_data <- lr_data[,setdiff(names(lr_data),setdiff(c('N','Observed','Expected'),showCols))]
  if (tableOnly) return(lr_data)
  to_indent <- setdiff(1:nrow(lr_data),c(1,nrow(lr_data),nrow(lr_data)-1))
  argL <- list(lr_data,to_indent=to_indent,caption=caption,align=c('lrrrr'))
  if (!missing(fontsize)) argL[['fontsize']] <- fontsize
  do.call(outTable, argL)


}

#' Summarise survival data by group
#'
#' Displays event counts, median survival time and survival rates at specified
#' times points for the entire cohort and by group. The logrank test of
#' differences in survival curves is displayed.
#'
#' This summary table is supplied for simple group comparisons only. To examine
#' differences in groups with stratification see \code{\link{rm_survdiff}}. To
#' summarise differences in survival rates controlling for covariates see
#' \code{\link{rm_survtime}}.
#'
#' @param data data frame containing survival data
#' @param time string indicating survival time variable
#' @param status string indicating event status variable
#' @param group string or character vector indicating the variable(s) to group
#'   observations by. If this is left as NULL (the default) then summaries are
#'   provided for the entire cohort.
#' @param survtimes numeric vector specifying when survival probabilities should
#'   be calculated.
#' @param survtimeunit unit of time to suffix to the time column label if
#'   survival probabilities are requested, should be plural
#' @param survtimesLbls if supplied, a vector the same length as survtimes with
#'   descriptions (useful for displaying years with data provided in months)
#' @param CIwidth width of the survival probabilities, default is 95%
#' @param unformattedp boolean indicating if you would like the p-value to be
#'   returned unformatted (ie not rounded or prefixed with '<'). Should be used
#'   in conjunction with the digits argument.
#' @param conf.type type of confidence interval see \code{\link{survfit}} for
#'   details. Default is 'log'.
#' @param na.action default is to omit missing values, but can be set to throw
#'   and error using na.action='na.fail'
#' @param showCounts boolean indicating if the at risk, events and censored
#'   columns should be output, default is TRUE
#' @param digits the number of digits in the survival rate, default is 2.
#' @param caption table caption for markdown output
#' @param tableOnly should a dataframe or a formatted object be returned
#' @param fontsize PDF/HTML output only, manually set the table fontsize
#' @importFrom  survival survfit Surv
#' @seealso \code{\link{survfit}}
#' @return A character vector of the survival table source code, unless tableOnly=TRUE in
#'   which case a data frame is returned
#' @export
#' @examples
#' # Simple median survival table
#' rm_survsum(data=pembrolizumab,time='os_time',status='os_status')
#'
#' # Survival table with yearly survival rates
#' rm_survsum(data=pembrolizumab,time='os_time',status='os_status',
#' survtimes=c(12,24),survtimesLbls=1:2, survtimeunit='yr')
#'
#' #Median survival by group
#' rm_survsum(data=pembrolizumab,time='os_time',status='os_status',group='sex')
#'
#' # Survival Summary by cohort, displayed in years
#' rm_survsum(data=pembrolizumab,time='os_time',status='os_status',
#' group="cohort",survtimes=seq(12,72,12),
#' survtimesLbls=seq(1,6,1),
#' survtimeunit='years')
#'
#' # Survival Summary by Sex and ctDNA group
#' rm_survsum(data=pembrolizumab,time='os_time',status='os_status',
#' group=c('sex','change_ctdna_group'),survtimes=c(12,24),survtimeunit='mo')
#'
rm_survsum <- function(data,time,status,group=NULL,survtimes=NULL,
                       survtimeunit,survtimesLbls=NULL,CIwidth=0.95,unformattedp=FALSE,
                       conf.type='log',
                       na.action='na.omit',showCounts=TRUE,digits=2,caption=NULL,tableOnly=FALSE,fontsize){
  if (missing(data)) stop('data is a required argument')
  if (missing(time)) stop('time is a required argument')
  if (missing(survtimeunit)) if (!is.null(survtimes)) stop('survtimeunit must be specified if survtimes are set. Example survtimeunit="year"')
  if (!is.null(survtimes)) timelbl <- paste0('Time (',survtimeunit,')')
  if  (!inherits(data,'data.frame')) stop('data must be supplied as a data frame')
  if( !inherits(time,'character')| length(time)>1)
    stop('time must be supplied as a string indicating a variable in data')
  if (!inherits(status,'character')| length(status)>1)
    stop('status must be supplied as a string indicating a variable in data')
  if (!is.null(survtimes)) if (is.null(survtimesLbls))  survtimesLbls <- survtimes
  if (length(survtimesLbls)!=length(survtimes)) stop('If supplied, the survtimesLbls vector must be the same length as survtime')
  if (unformattedp) formatp <- function (x,...){x}

  data <- data.frame(data)
  n_cols <- c('n.risk','n.event','n.censor')

  missing_vars = na.omit(setdiff(c(time, status,group), names(data)))
  if (length(missing_vars) > 0) stop(paste("These variables are not in the data:\n",
                                           paste0(missing_vars,collapse=csep())))

  sfit <- survival::survfit(as.formula(paste0("survival::Surv(",time,',',status,') ~',ifelse(is.null(group),"1",paste(group,collapse='+')))),
                            data = data,conf.type=conf.type,conf.int=CIwidth)
  nFit <- sum(sfit$n,na.rm=TRUE)
  if (nrow(data)-nFit==1) {
    message('1 observation with missing data was removed.')
  } else if (nrow(data)>nFit) message(paste(nrow(data)-nFit ,'observations with missing data were removed.'))

  if (!is.null(survtimes)) {
    ofit <- summary(sfit,times=survtimes,extend=!is.null(survtimes))
    colsToExtract <- which(names(sfit) %in% c("strata","time","sr","surv","lower","upper"))
    tb <- data.frame(do.call(cbind,ofit[colsToExtract]))
    tb$sr <- apply(tb[,c("surv","lower","upper")], 1, psthr,digits)
    if (!is.null(group)){
      tb$strata <- factor(tb$strata,levels=unique(as.numeric(ofit$strata)),labels = levels(ofit$strata))
      w <- matrix(nrow=length(unique(tb$strata)),ncol=length(unique(tb$time)),
                  dimnames=list(unique(tb$strata),unique(tb$time)))
      for (i in 1:nrow(tb)) w[which(rownames(w)==tb$strata[i]),which(colnames(w)==tb$time[i])] <- tb$sr[i]
    } else{
      w <- matrix(nrow=1,ncol=length(unique(tb$time)))
      colnames(w) <- unique(tb$time)
      for (i in 1:nrow(tb)) w[1,which(colnames(w)==tb$time[i])] <- tb$sr[i]
    }
  } else w <- NULL

  mtbl <- summary(sfit)$table
  if (inherits(mtbl,'matrix')){
    m_CI <- apply(mtbl[,grep('median|LCL|UCL',colnames(mtbl))],1,function(x) psthr(x,y=digits))
    lr <- survival::survdiff(as.formula(paste0("survival::Surv(",time,',',status,') ~',
                                               paste0(group,collapse = '+'))),data = data)
    nt <- paste0(lr$obs,'/',lr$n)
    w <- cbind(nt,m_CI,w)
    df <- length(lr$obs)-1
    gl <- rownames(w)
    for (v in group) gl <- gsub(paste0(v,'='),"",gl)
    tab <- cbind(gl,data.frame(w))
    rownames(tab) <- NULL
    nm <- c(paste(group,collapse = ','),'Events/Total',paste0('Median (',CIwidth*100,'%CI)'))
    if (!is.null(survtimes)) nm <- c(nm,paste0(survtimesLbls,survtimeunit," (", CIwidth*100,"% CI)"))
    names(tab) <- nm
    tab <- rbind(tab,c(rep("",length(survtimes)),'Log Rank Test','ChiSq',paste(niceNum(lr$chisq,1),'on',df,'df')))
    tab <- rbind(tab,c(rep("",length(survtimes)+1),'p-value',formatp(pchisq(lr$chisq,df,lower.tail = FALSE))))
  } else{
    m_CI <- psthr(mtbl[grep('median|LCL|UCL',names(mtbl))],y=digits)
    nt <-paste0(mtbl["events"],"/",mtbl["n.start"])
    w <- cbind(nt,m_CI,w)
    tab <- data.frame(w)
    rownames(tab) <- NULL
    nm <- c('Events/Total',paste0('Median (',CIwidth*100,'%CI)'))
    if (!is.null(survtimes)) nm <- c(nm,paste0(survtimesLbls,survtimeunit," (", CIwidth*100,"% CI)"))
    names(tab) <- nm

  }

  if (tableOnly){
    return(tab)
  }
  argL <- list(tab,caption=caption,
           align = paste0('l',paste0(rep('r',ncol(tab)-1),collapse = '')))
  if (!missing(fontsize)) argL[['fontsize']] <- fontsize
  do.call(outTable, argL)

}

#' Display survival rates and events for specified times
#'
#' This is a wrapper for the survfit function to output a tidy display for
#' reporting. Either Kaplan Meier or Cox Proportional Hazards models may be used
#' to estimate the survival probabilities.
#'
#' If covariates are supplied then a Cox proportional hazards model is fit for
#' the entire cohort and each strata. Otherwise the default is for Kaplan-Meier
#' estimates. Setting type = 'PH' will force a proportional hazards model.
#' @param data data frame containing survival data
#' @param time string indicating survival time variable
#' @param status string indicating event status variable
#' @param covs character vector with the names of variables to adjust for in
#'   coxph fit
#' @param strata string indicating the variable to group observations by. If
#'   this is left as NULL (the default) then event counts and survival rates are
#'   provided for the entire cohort.
#' @param type survival function, if no covs are specified defaults to
#'   Kaplan-Meier, otherwise the Cox PH model is fit. Use type='PH' to fit a Cox
#'   PH model with no covariates.
#' @param survtimes numeric vector specifying when survival probabilities should
#'   be calculated.
#' @param survtimeunit unit of time to suffix to the time column label if
#'   survival probabilities are requested, should be plural
#' @param strata.prefix character value describing the grouping variable
#' @param survtimesLbls if supplied, a vector the same length as survtimes with
#'   descriptions (useful for displaying years with data provided in months)
#' @param showCols character vector specifying which of the optional columns to
#'   display, defaults to c('At Risk','Events','Censored')
#' @param CIwidth width of the survival probabilities, default is 95%
#' @param conf.type type of confidence interval see \code{\link{survfit}} for
#'   details. Default is 'log'.
#' @param na.action default is to omit missing values, but can be set to throw
#'   and error using na.action='na.fail'
#' @param showCounts boolean indicating if the at risk, events and censored
#'   columns should be output, default is TRUE
#' @param digits the number of digits in the survival rate, default is 2.
#' @param caption table caption for markdown output
#' @param tableOnly should a dataframe or a formatted object be returned
#' @param fontsize PDF/HTML output only, manually set the table fontsize
#' @importFrom  survival survfit Surv
#' @seealso \code{\link{survfit}}
#' @return A character vector of the survival table source code, unless tableOnly=TRUE in
#'   which case a data frame is returned
#' @export
#' @examples
#' # Kaplan-Mieir survival probabilities with time displayed in years
#' rm_survtime(data=pembrolizumab,time='os_time',status='os_status',
#' strata="cohort",type='KM',survtimes=seq(12,72,12),
#' survtimesLbls=seq(1,6,1),
#' survtimeunit='years')
#'
#' # Cox Proportional Hazards survivial probabilities
#' rm_survtime(data=pembrolizumab,time='os_time',status='os_status',
#' strata="cohort",type='PH',survtimes=seq(12,72,12),survtimeunit='months')
#'
#' # Cox Proportional Hazards survivial probabilities controlling for age
#' rm_survtime(data=pembrolizumab,time='os_time',status='os_status',
#' covs='age',strata="cohort",survtimes=seq(12,72,12),survtimeunit='months')
#'
rm_survtime <- function(data,time,status,covs=NULL,strata=NULL,type='KM',survtimes,
                        survtimeunit,strata.prefix=NULL,survtimesLbls=NULL,
                        showCols=c('At Risk','Events','Censored'),CIwidth=0.95,conf.type='log',
                        na.action='na.omit',showCounts=TRUE,digits=2,caption=NULL,tableOnly=FALSE,fontsize){
  if (missing(data)) stop('data is a required argument')
  if (missing(time)) stop('time is a required argument')
  if (missing(status)) stop('status is a required argument')
  if (missing(survtimes)) stop('survtimes must be specified as a numeric vector')
  if (missing(survtimeunit)) timelbl <- 'Time' else timelbl <- paste0('Time (',survtimeunit,')')
  if  (!inherits(data,'data.frame')) stop('data must be supplied as a data frame')
  if( !inherits(time,'character')| length(time)>1)
    stop('time must be supplied as a string indicating a variable in data')
  if (!inherits(status,'character')| length(status)>1)
    stop('status must be supplied as a string indicating a variable in data')
  if (!missing(covs)){
    if (!inherits(covs,'character'))
      stop('covs must be supplied as a character vector or string indicating variables in data')
  }
  if (is.null(survtimesLbls)) survtimesLbls = survtimes
  if (length(survtimesLbls)!=length(survtimes)) stop('If supplied, the survtimesLbls vector must be the same length as survtime')
  if (!missing(strata) & !missing(covs)){
    if (length(intersect(covs,strata))>0) stop('A variable can appear in covs or strata but not both.')

  }
  timelbl <- paste0('Time (',survtimeunit,')')
  data <- data.frame(data)
  n_cols <- c('n.risk','n.event','n.censor')

  missing_vars = na.omit(setdiff(c(time, status,covs,strata), names(data)))
  if (length(missing_vars) > 0) stop(paste("These variables are not in the data:\n",
                                           paste0(missing_vars,collapse=csep())))

  if (!is.null(covs)) type <- 'PH'
  if (is.null(covs) & type =='PH') covs <- '1'

  survtime_est <- NULL
  if (type=='KM'){
    sfit <- survival::survfit(as.formula(paste0("survival::Surv(",time,',',status,') ~1')),
                              data = data,conf.type=conf.type,conf.int=CIwidth)
  }  else{ # CoxPH
    sfit<- survival::survfit(survival::coxph(as.formula(paste0("Surv(",time,',',status,') ~',
                                                               paste(covs,collapse='+'))),
                                             data = data),conf.type=conf.type,conf.int=CIwidth)
  }
  colsToExtract <- which(names(sfit) %in% c("time","n.risk","n.event","n.censor","sr","surv","lower","upper"))
  if (nrow(data)-sfit$n==1) {
    message('1 observation with missing data was removed.')
  } else if (nrow(data)>sfit$n) message(paste(nrow(data)-sfit$n ,'observations with missing data were removed.'))


  ofit <- summary(sfit,times=survtimes,extend=!is.null(survtimes))
  ofit <- data.frame(do.call(cbind,ofit[colsToExtract]))
  ofit$N <- sfit$n
  survtime_est[['Overall']] <- ofit
  header <- data.frame(matrix(nrow = 1,ncol=ncol(ofit)))
  names(header) <- names(ofit)

  if (!is.null(strata)){
    if (inherits(data[[strata]],c('ordered','factor'))){
      levelnames <- levels(data[[strata]])
      levelnames <- levelnames[levelnames %in% unique(data[[strata]])]
    } else  levelnames <- unique(data[[strata]])

    for (g in levelnames){
      if (type=='KM'){
        sfit<- survival::survfit(as.formula(paste0("survival::Surv(",time,',',status,') ~1')),
                                 data = data[data[[strata]]==g,],conf.type=conf.type,conf.int=CIwidth)
      }
      else{ # CoxPH
        sfit<- survival::survfit(survival::coxph(as.formula(paste0("Surv(",time,',',status,') ~',
                                                                   paste(covs,sep='+'))),
                                                 data = data[data[[strata]]==g,]),conf.type=conf.type,conf.int=CIwidth)
      }
      gfit <- summary(sfit,times=survtimes,extend=!is.null(survtimes))
      gfit <- data.frame(do.call(cbind,gfit[colsToExtract]))
      gfit$N <- sfit$n
      survtime_est[[as.character(g)]] <- gfit
    }}

  out <- lapply(names(survtime_est),function(x){
    xtbl <- survtime_est[[x]]
    xtbl$sr <- apply(xtbl[,c("surv","lower","upper")], 1, psthr,digits)
    hdtxt <- ifelse(x=='Overall',x,
                    ifelse(is.null(strata.prefix),x,
                           paste(strata.prefix,x,sep='=')))

    header[1,] <- c(hdtxt,rep('',ncol(header)-1))
    header[1,n_cols[1]] = xtbl$N[1]
    header$sr <- ''
    xnew <- rbind(header,xtbl)

    return(xnew[,c(names(xtbl)[1],n_cols,'sr')])
  })
  tab <- do.call(rbind, out)
  rownames(tab) <- NULL
  tab_times <- tab
  kp <- apply(tab[,n_cols],1,function(x){
    sum(as.numeric(x),na.rm=TRUE) > 0
  })

  tab <- tab[kp,]
  names(tab) <- gsub('time',timelbl,names(tab))
  names(tab) <- gsub('sr',paste0("Survival Rate (", CIwidth*100,"\\% CI)"),names(tab))
  names(tab) <- gsub('n.event','Events',names(tab))
  names(tab) <- gsub('n.censor','Censored',names(tab))
  names(tab) <- gsub('n.risk','At Risk',names(tab))
  tab <- tab[,setdiff(names(tab),setdiff(c('At Risk','Events','Censored'),showCols))]

  if (tableOnly){
    return(tab)
  }
  to_indent  <- which(tab[[timelbl]] %in% survtimes)
  argL <- list(tab,to_indent=to_indent,caption=caption,
           align = ifelse(showCounts,'lrrrr','lr'))
  if (!missing(fontsize)) argL[['fontsize']] <- fontsize
  do.call(outTable, argL)

}

