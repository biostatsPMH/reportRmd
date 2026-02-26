#' Retrieve columns number from spreadsheet columns specified as unquoted letters
#' @param ... unquoted excel column headers (i.e. excelCol(A,CG,AA)) separated by commas
#' @importFrom rlang as_string
#' @returns a numeric vector corresponding to columns in a spreadsheet
#' @export
#' @examples
#' ## Find the column numbers for excel columns AB, CE and BB
#' excelCol(AB,CE,bb)
#' ## Get the columns between A and K and Z
#' excelCol(A-K,Z)
excelCol<- function(...){
  args <- as.list(match.call())[-1]
  len <- lapply(args,nchar)
  if (any(unlist(len)>2)) stop('Columns names must be between A and ZZ. \nVariable names can not be used in this function.')
  if (any(grepl('-',args))) {
    rtn<-unname(unlist(lapply(args,function(x) {
      x1 <- gsub("-",",",x)
      x1 <- x1[which(x1!=",")]
      if (sum(grepl("[^A-Za-z]",x1))>0)
        stop('Only valid Excel column names can be supplied, separated by commas or hyphens.')
      x2 <- sapply(x1, xcn)
      if (length(x2)>1) x2 <- seq(x2[1],x2[2],1)
      return(x2)
    })))
  } else{
    args <-unname(unlist(lapply(args,function(x) {rlang::as_string(x)})))
    if (sum(unlist(lapply(args, function(x) grepl("[^A-Za-z]",x))))>0) {
      stop('Only valid Excel column names can be supplied, separated by commas or hyphens.')
    }
    rtn <- xcn(args)
    names(rtn) <- toupper(names(rtn))
  }
  return(rtn)
}



#' Retrieve spreadsheet column letter-names from columns indices
#'
#' Creates a vector of spreadsheet-style letter-names corresponding to column numbers
#'
#' This is the inverse function of excelCol
#'
#' @param columnIndices vector of integer column indices
#' @returns a character vector corresponding to the spreadsheet column headings
#' @export
#' @examples
#' ## Find the column numbers for excel columns AB, CE and BB
#' colIndices <- excelCol(AB,CE,bb)
#' ## Go back to the column names
#' excelColLetters(colIndices)
excelColLetters<- function(columnIndices){
  if (!is.integer(as.integer(columnIndices))) stop('columnIndices must be a vector of numeric column indices')
  out <- sapply(columnIndices,function(x){
    rmd <- x%%26
    quo <- x%/% 26
    if (rmd==0){
      rmd <-26
      quo <- quo-1
    }
    fl <- ifelse(quo>0,LETTERS[quo],"")
    sl <- ifelse(rmd>0,LETTERS[rmd],"")
    return(paste0(fl,sl))
  })
  names(out) <- columnIndices
  return(out)
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
  }      else {
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
#' @seealso \code{\link[cmprsk:crr]{cmprsk::crr}}
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
  argList <- as.list(match.call()[-1])
  k<-as.character(f)[3]
  covs<-removedollar(k)
  ff<-modelmatrix(f,data)
  m1<-cmprsk::crr(ff[[1]][,1],ff[[1]][,2],ff[[2]])
  m2 <- summary(m1)
  m1$formula <- paste("~",covs)
  m1$terms <- covs
  covstr <- extract_terms(covs)
  attr(m1$terms,"term.labels") <- covstr
  names(covstr) <- covstr
  attr(m1$terms,"dataClasses") <- sapply(covstr,function(x)class(data[[x]]))

  vrs <- sapply(covstr, function(x){
    strsplit(x,":")[[1]]
  })
  xvr <- unique(unlist(vrs))
  m1$model <- na.omit(data[,c(colnames(ff[[1]]),intersect(names(data),xvr))])
  attr(m1$model,"terms") <- paste(paste(colnames(ff[[1]]),collapse = "+"),
                                  "~", covstr, sep = "")

  m1$coeffTbl <- m2$coef
  m1$coefficients <- m1$coef
  m1$call<-as.call(list(f,data=argList$data))
  # m1$data <- data
  attr(m1,"class") <- c(attr(m1,"class"),"crrRx")
  return(m1)
}

extract_terms <- function(terms) {
  # Initialize an empty vector to store the terms
  ts <- c()
  pp <- unlist(strsplit(terms, "\\+"))
  for (part in pp) {
    part <- trimws(part)
    if (grepl("\\*", part)) {
      vars <- unlist(strsplit(part, "\\*"))
      var1 <- trimws(vars[1])
      var2 <- trimws(vars[2])
      ts <- c(ts, var1, var2, paste(var1, var2, sep = ":"))
    } else {
      ts <- c(ts, part)
    }
  }
  return(ts)
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
#'   table. Can only be obtained if pvalue is also requested. Effect sizes
#'   calculated include Cramer's V for categorical variables, Cohen's d,
#'   Wilcoxon r, or Eta-squared for numeric/continuous variables.
#' @param show.tests boolean indicating if the type of statistical test and
#'   effect size used should be shown in a column beside the pvalues. Ignored if
#'   pvalue=FALSE.
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
#' @importFrom stats lm sd setNames aov
#' @importFrom rstatix cramer_v eta_squared
#' @seealso \code{\link{fisher.test}},\code{\link{chisq.test}},
#'   \code{\link{wilcox.test}},\code{\link{kruskal.test}},and
#'   \code{\link{anova}}
#' @references Ellis, P.D. (2010) The essential guide to effect sizes:
#'   statistical power, meta-analysis, and the interpretation of research
#'   results. Cambridge: Cambridge University
#'   Press.\doi{10.1017/CBO9780511761676}
#' @references Lakens, D. (2013)  Calculating and reporting effect sizes to
#'   facilitate cumulative science: a practical primer for t-tests and ANOVAs.
#'   Frontiers in Psychology, 4; 863:1-12. \doi{10.3389/fpsyg.2013.00863}
covsum <- function (data, covs, maincov = NULL, digits = 1, numobs = NULL,
                    markup = FALSE, sanitize = FALSE, nicenames = TRUE, IQR = FALSE,
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
    if (inherits(data[[v]], "difftime"))
      data[[v]] <- as.numeric(data[[v]])
    if (inherits(data[[v]], "logical"))
      data[[v]] <- factor(data[[v]])
    if (inherits(data[[v]], "character"))
      data[[v]] <- factor(data[[v]])
    if (inherits(data[[v]], c("POSIXt"))) {
      covs <- setdiff(covs, v)
      message(paste("POSIXt can not be summarised in this version of reportRmd.\n The variable",
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
    else excludeLevel = NA
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
              e_type <- "Cramer's V"
              e <- try(sqrt((stats::chisq.test(pdata[[maincov]],
                                               pdata[[cov]])$statistic/N)/
                              (min(length(unique(pdata[[maincov]]))-1,length(unique(pdata[[cov]]))-1))), silent = TRUE)
            }
            if (is.error(p)) {
              message("Using MC sim. Use set.seed() prior to function for reproducible results.")
              p <- try(stats::fisher.test(pdata[[maincov]],
                                          pdata[[cov]], simulate.p.value = T)$p.value,
                       silent = TRUE)
              p_type <- "MC sim"
              if (effSize) {
                e_type <- "Cramer's V"
                e <- try(sqrt((stats::chisq.test(pdata[[maincov]],
                                                 pdata[[cov]])$statistic/N)/
                                (min(length(unique(pdata[[maincov]]))-1,length(unique(pdata[[cov]]))-1))), silent = TRUE)
              }
            }
          }
          else {
            p_type = "Chi Sq"
            p = try(stats::chisq.test(pdata[[maincov]],
                                      pdata[[cov]])$p.value, silent = TRUE)
            if (effSize) {
              e_type <- "Cramer's V"
              e <- try(sqrt((stats::chisq.test(pdata[[maincov]],
                                               pdata[[cov]])$statistic/N)/
                              (min(length(unique(pdata[[maincov]]))-1,length(unique(pdata[[cov]]))-1))), silent = TRUE)
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
                e_type <- "Wilcoxon r"
                e <- try(ifelse(is.finite(qnorm(stats::wilcox.test(data[[cov]]~data[[maincov]], data=data)$p.value/2)),
                                abs(qnorm(stats::wilcox.test(data[[cov]]~data[[maincov]], data=data)$p.value/2))/sqrt(N),
                                abs(qnorm(0.0001/2))/sqrt(N)),
                         silent = TRUE)
              }
            }
            else {
              p_type = "Kruskal Wallis"
              p <- try(stats::kruskal.test(data[[cov]] ~
                                             data[[maincov]])$p.value, silent = T)
              if (effSize) {
                e_type <- "Eta sq"
                e <- try(summary(stats::aov(data[[cov]]~data[[maincov]]))[[1]]$`Sum Sq`[1]/
                           (summary(stats::aov(data[[cov]]~data[[maincov]]))[[1]]$`Sum Sq`[1]+summary(stats::aov(data[[cov]]~data[[maincov]]))[[1]]$`Sum Sq`[2]))
              }
            }
          }
          else {
            if (length(unique(data[[maincov]])) == 2) {
              p_type = "t-test"
              p <- try(stats::t.test(data[[cov]] ~ data[[maincov]])$p.value,
                       silent = TRUE)
              if (effSize) {
                e_type <- "Cohen's d"
                e <- try(abs(2*stats::t.test(data[[cov]] ~ data[[maincov]])$statistic/sqrt(N)), silent = TRUE)
              }
            }
            else {
              p_type = "ANOVA"
              p <- try(stats::anova(stats::lm(data[[cov]] ~
                                                data[[maincov]]))[5][[1]][1], silent = TRUE)
              if (effSize) {
                e_type <- "Eta sq"
                e <- try(summary(stats::aov(data[[cov]]~data[[maincov]]))[[1]]$`Sum Sq`[1]/
                           (summary(stats::aov(data[[cov]]~data[[maincov]]))[[1]]$`Sum Sq`[1]+summary(stats::aov(data[[cov]]~data[[maincov]]))[[1]]$`Sum Sq`[2]))
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
        sumCov <- try(round(summary(subdata[[cov]]), digits))

        if (sumCov[4] == "NaN") {
          meansd <- ""
          mmm <- ""
          if (all.stats)
            mmm = c("", "")
        }
        else if (inherits(subdata[[cov]],"Date")) {
          meansd <- paste(as.Date(floor(as.numeric(niceNum(sumCov["Mean"], digits))), origin="1970-01-01"),
                          " (", niceNum(sd(subdata[[cov]], na.rm = T),
                                        digits), " days)", sep = "")
          mmm <- if (IQR | all.stats) {
            if (all(as.Date(c(sumCov["Median"], sumCov["1st Qu."],
                              sumCov["3rd Qu."]), origin="1970-01-01") == as.Date(floor(c(sumCov["Median"],
                                                                                          sumCov["1st Qu."], sumCov["3rd Qu."])), origin="1970-01-01"))) {
              paste(as.Date(as.numeric(sumCov["Median"]), origin="1970-01-01"), " (", as.Date(as.numeric(sumCov["1st Qu."]), origin="1970-01-01"),
                    csep(), as.Date(as.numeric(sumCov["3rd Qu."]), origin="1970-01-01"), ")", sep = "")
            }
            else {
              paste(as.Date(floor(as.numeric(niceNum(sumCov["Median"], digits))), origin="1970-01-01"),
                    " (", as.Date(floor(as.numeric(niceNum(sumCov["1st Qu."], digits))), origin="1970-01-01"),
                    csep(), as.Date(floor(as.numeric(niceNum(sumCov["3rd Qu."], digits))), origin="1970-01-01"),
                    ")", sep = "")
            }
          }
          else {
            if (all(as.Date(c(sumCov["Median"], sumCov["Min."],
                              sumCov["Max."]), origin="1970-01-01") == as.Date(floor(c(sumCov["Median"],
                                                                                       sumCov["Min."], sumCov["Max."])), origin="1970-01-01"))) {
              paste(as.Date(as.numeric(sumCov["Median"]), origin="1970-01-01"), " (", as.Date(as.numeric(sumCov["Min."]), origin="1970-01-01"),
                    csep(), as.Date(as.numeric(sumCov["Max."]), origin="1970-01-01"), ")", sep = "")
            }
            else {
              paste(as.Date(as.numeric(niceNum(sumCov["Median"], digits)), origin="1970-01-01"),
                    " (", as.Date(as.numeric(niceNum(sumCov["Min."], digits)), origin="1970-01-01"),
                    csep(), as.Date(as.numeric(niceNum(sumCov["Max."], digits)), origin="1970-01-01"),
                    ")", sep = "")
            }
          }
          if (all.stats) {
            mmm <- c(mmm, if (all(as.Date(c(sumCov["Min."], sumCov["Max."]), origin="1970-01-01") ==
                                  as.Date(as.numeric(floor(c(sumCov["Min."], sumCov["Max."]))), origin="1970-01-01"))) {
              paste("(", as.Date(as.numeric(sumCov["Min."]), origin="1970-01-01"), csep(), as.Date(as.numeric(sumCov["Max."]), origin="1970-01-01"),
                    ")", sep = "")
            } else {
              paste("(", as.Date(as.numeric(niceNum(sumCov["Min."], digits)), origin="1970-01-01"),
                    csep(), as.Date(as.numeric(niceNum(sumCov["Max."], digits)), origin="1970-01-01"),
                    ")", sep = "")
            })
          }
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
  varID <- do.call("c",lapply(out,function(x){
    return(stats::setNames(c(TRUE,rep(FALSE,nrow(x)-1)),x[,1]))
  }))
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
  attr(table,"varID") <- varID
  return(table)
}


#' Extract and prepare forest plot data from m_summary output
#'
#' @param summary_output output from m_summary with for_plot = TRUE
#' @param model_type character, "adjusted" or "unadjusted"
#' @param digits number of digits for rounding
#' @keywords internal
prepare_forest_data <- function(summary_output, model_type = "adjusted", digits = 2) {

  # Extract key columns from m_summary output
  tab <- summary_output

  # Remove header rows - we'll handle variable names in labels instead
  tab <- tab[is.na(tab$header) | tab$header == FALSE, ]

  # Create estimate label
  tab$estimate.label <- ifelse(
    !is.na(tab$ref) & tab$ref == TRUE,
    "Reference",
    tab$Est_CI
  )

  # Extract numeric estimate and confidence bounds
  tab$estimate <- tab$est
  tab$conf.low <- tab$lwr
  tab$conf.high <- tab$upr

  # Set level name
  tab$level.name <- tab$lvl

  # Set variable name (for grouping)
  tab$var.name <- tab$var

  # Add model type
  tab$type <- model_type

  # Keep relevant columns
  cols <- c("var", "var.name", "level.name", "estimate", "estimate.label",
            "conf.low", "conf.high", "p_value", "n", "type", "ref")

  # Add Events if present
  if ("Events" %in% names(tab)) {
    cols <- c(cols, "Events")
  }

  tab <- tab[, cols]

  # Rename for consistency
  names(tab)[names(tab) == "var"] <- "variable"
  names(tab)[names(tab) == "p_value"] <- "p.value"
  names(tab)[names(tab) == "n"] <- "N"
  if ("Events" %in% names(tab)) {
    names(tab)[names(tab) == "Events"] <- "Event"
  }

  return(tab)
}

#' Order forest plot data by risk and factor levels
#'
#' @param tab data frame prepared by prepare_forest_data
#' @param model original model object (for factor level info)
#' @param data original data
#' @keywords internal
order_forest_data <- function(tab, model, data) {

  # Initialize level_order column for all rows
  tab$level_order <- NA

  # Get model variables and their classes
  model_frame <- model.frame(model)
  model_vars <- names(model_frame)[-1]  # Exclude response

  # For each variable, determine ordering
  for (var in unique(tab$variable)) {
    var_rows_idx <- which(tab$variable == var)

    # Check if variable is a factor in original data
    if (var %in% names(data) && is.factor(data[[var]])) {
      # Order by factor levels
      factor_levels <- levels(data[[var]])

      # Reference level gets order 0
      ref_idx <- var_rows_idx[which(!is.na(tab$ref[var_rows_idx]) &
                                      tab$ref[var_rows_idx] == TRUE)]
      if (length(ref_idx) > 0) {
        tab$level_order[ref_idx] <- 0
      }

      # Other levels get order based on factor
      for (i in seq_along(factor_levels)) {
        level_idx <- var_rows_idx[which(tab$level.name[var_rows_idx] == factor_levels[i])]
        if (length(level_idx) > 0) {
          # Assign order to all matching indices that don't already have one
          needs_order <- is.na(tab$level_order[level_idx])
          if (any(needs_order)) {
            tab$level_order[level_idx[needs_order]] <- i
          }
        }
      }

    } else {
      # For continuous or non-factor variables, order by estimate
      var_estimates <- tab$estimate[var_rows_idx]
      tab$level_order[var_rows_idx] <- rank(var_estimates, na.last = TRUE)

      # But put reference first if it exists
      ref_idx <- var_rows_idx[which(!is.na(tab$ref[var_rows_idx]) &
                                      tab$ref[var_rows_idx] == TRUE)]
      if (length(ref_idx) > 0) {
        tab$level_order[ref_idx] <- 0
      }
    }
  }

  # Calculate max estimate per variable for between-variable ordering
  var_max_est <- tab |>
    dplyr::group_by(variable) |>
    dplyr::summarise(max_est = max(estimate, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(max_est)) |>
    dplyr::mutate(var_order = dplyr::row_number())

  # Merge var_order back
  tab <- tab |>
    dplyr::left_join(var_max_est[, c("variable", "var_order")], by = "variable")

  # Sort by var_order (descending max OR), then level_order (ascending)
  tab <- tab |>
    dplyr::arrange(var_order, level_order)

  # Mark first level for each variable (for label formatting)
  tab$is_first_level <- FALSE
  for (var in unique(tab$variable)) {
    var_rows <- which(tab$variable == var)
    if (length(var_rows) > 0) {
      tab$is_first_level[var_rows[1]] <- TRUE
    }
  }

  return(tab)
}


#' Create a multivariable forest plot using ggplot2
#'
#' This function creates forest plots from fitted regression models, with optional
#' inclusion of unadjusted estimates. It uses m_summary for robust data extraction
#' and properly handles factor level ordering and reference levels.
#'
#' @param model an object output from the glm or geeglm function, must be from a
#'   logistic or log-link regression
#' @param data dataframe containing your data (required if include_unadjusted = TRUE)
#' @param include_unadjusted logical, should unadjusted estimates be included?
#'   Default is FALSE
#' @param conf.level controls the width of the confidence interval (default 0.95)
#' @param colours can specify colours for risks less than, equal to, and greater
#'   than 1.0. Default is green, black, red
#' @param showEst logical, should the risks be displayed on the plot in text?
#'   Default is TRUE
#' @param showRef logical, should reference levels be shown? Default is TRUE
#' @param digits number of digits to use displaying estimates (default 2)
#' @param logScale logical, should OR/RR be shown on log scale? Defaults to TRUE.
#'   See https://doi.org/10.1093/aje/kwr156 for why you may prefer a linear scale
#' @param nxTicks Number of tick marks for x-axis (default 5)
#' @param showN Show number of observations per variable and category (default TRUE)
#' @param showEvent Show number of events per variable and category (default TRUE)
#' @param xlim Numeric vector of length 2 specifying x-axis limits (default c(0.2, 5))
#' @import ggplot2
#' @importFrom scales log_breaks
#' @keywords plot
#' @return a ggplot object
#' @export
#' @examples
#' data("pembrolizumab")
#' glm_fit <- glm(orr ~ change_ctdna_group + sex + age + l_size,
#'                data = pembrolizumab, family = 'binomial')
#'
#' # Adjusted only
#' forestplotMV(glm_fit, data = pembrolizumab)
#'
#' # Both adjusted and unadjusted
#' forestplotMV(glm_fit, data = pembrolizumab, include_unadjusted = TRUE)
forestplotMV <- function(model, data = NULL, include_unadjusted = FALSE,
                         conf.level = 0.95, colours = "default",
                         showEst = TRUE, showRef = TRUE,
                         digits = getOption("reportRmd.digits", 2),
                         logScale = getOption("reportRmd.logScale", TRUE),
                         nxTicks = 5, showN = TRUE, showEvent = TRUE,
                         xlim = c(0.2, 5)) {

  # Determine axis label based on model family
  if (inherits(model, "glm")) {
    if (model$family$link == "log") {
      x_lab_adj <- "Adjusted Relative Risk"
      x_lab_unadj <- "Unadjusted Relative Risk"
    } else if (model$family$link == "logit") {
      x_lab_adj <- "Adjusted Odds Ratio"
      x_lab_unadj <- "Unadjusted Odds Ratio"
    } else {
      stop("model must be a logit or log link fit")
    }
  } else {
    x_lab_adj <- "Adjusted Odds Ratio"
    x_lab_unadj <- "Unadjusted Odds Ratio"
  }

  # Set x-axis label
  x_lab <- if (include_unadjusted) {
    gsub("Adjusted |Unadjusted ", "", x_lab_adj)
  } else {
    x_lab_adj
  }

  # Extract data from model if not provided (needed for order_forest_data)
  if (is.null(data)) {
    data <- get_model_data(model)
  }

  # Get adjusted estimates using m_summary
  tab_mv <- m_summary(model, CIwidth = conf.level, digits = digits,
                      whichp = "levels", for_plot = TRUE)
  tab_mv <- prepare_forest_data(tab_mv, model_type = "Adjusted", digits = digits)

  # Get model N for comparison
  mv_n <- nrow(model.frame(model))

  # If unadjusted estimates requested, fit UV models
  if (include_unadjusted) {
    message("Fitting univariate models for each predictor")

    # Extract predictor names from model
    model_formula <- formula(model)
    response_var <- all.vars(model_formula)[1]
    predictor_vars <- attr(terms(model), "term.labels")

    # Fit UV model for each predictor
    uv_tabs <- list()
    for (pred in predictor_vars) {
      # Create UV formula
      uv_formula <- as.formula(paste(response_var, "~", pred))

      # Fit UV model with same family
      if (inherits(model, "glm")) {
        uv_model <- glm(uv_formula, data = data, family = model$family)
      } else {
        stop("Only glm models currently supported for unadjusted estimates")
      }

      # Get UV summary
      uv_tab <- m_summary(uv_model, CIwidth = conf.level, digits = digits,
                          whichp = "levels", for_plot = TRUE)
      uv_tabs[[pred]] <- prepare_forest_data(uv_tab, model_type = "Unadjusted",
                                             digits = digits)
    }

    # Combine UV tabs
    tab_uv <- do.call(rbind, uv_tabs)

    # Check for sample size differences and message user
    uv_n <- max(tab_uv$N, na.rm = TRUE)
    if (uv_n != mv_n) {
      message(sprintf(
        "Note: Adjusted model N=%d may differ from unadjusted model N=%d due to missing data in covariates",
        mv_n, uv_n
      ))
    }

    # Combine MV and UV data
    tab <- rbind(tab_mv, tab_uv)

    # Remove duplicate reference levels (keep only one)
    ref_rows <- which(!is.na(tab$ref) & tab$ref == TRUE)
    if (length(ref_rows) > 0) {
      # Group by variable and level.name to find duplicates
      ref_data <- tab[ref_rows, ]
      ref_data$row_idx <- ref_rows

      # Keep only first reference per variable+level combination
      duplicated_refs <- duplicated(paste(ref_data$variable, ref_data$level.name))
      rows_to_remove <- ref_data$row_idx[duplicated_refs]

      if (length(rows_to_remove) > 0) {
        tab <- tab[-rows_to_remove, ]
      }
    }

  } else {
    tab <- tab_mv
  }

  # Order data by risk
  tab <- order_forest_data(tab, model, data)

  # Remove reference levels if requested
  if (!showRef) {
    tab <- tab[is.na(tab$ref) | tab$ref == FALSE, ]
  }

  # Create y position values with tighter within-variable spacing
  # This makes levels within a variable "hang together" visually
  tab$y.val <- NA
  current_y <- nrow(tab)

  for (i in 1:nrow(tab)) {
    tab$y.val[i] <- current_y

    # Check if next row is a different variable (create larger gap)
    if (i < nrow(tab)) {
      if (tab$variable[i] != tab$variable[i + 1]) {
        # Larger gap between variables
        current_y <- current_y - 1.5
      } else {
        # Smaller gap within variable
        current_y <- current_y - 0.6
      }
    }
  }

  yvals <- tab$y.val

  # Format estimate labels
  tab$estimate.label[is.na(tab$estimate.label)] <- ""
  tab$estimate.label[tab$estimate.label == "Reference"] <- "(Reference)"

  # Create y-axis labels with variable name on first level
  if (showEst) {
    yLabels <- data.frame(
      y.pos = yvals,
      labels = ifelse(
        tab$is_first_level,
        # First level: show variable and level name
        paste(tab$variable, ": ", tab$level.name, ": ", tab$estimate.label, sep = ""),
        ifelse(
          !is.na(tab$type) & tab$type == "Unadjusted",
          # Unadjusted: show (Unadjusted) instead of level name
          paste("      (Unadjusted): ", tab$estimate.label, sep = ""),
          # Other levels: show level name
          paste("   ", tab$level.name, ": ", tab$estimate.label, sep = "")
        )
      )
    )
  } else {
    yLabels <- data.frame(
      y.pos = yvals,
      labels = ifelse(
        tab$is_first_level,
        paste(tab$variable, ": ", tab$level.name, sep = ""),
        ifelse(
          !is.na(tab$type) & tab$type == "Unadjusted",
          "      (Unadjusted)",
          paste("   ", tab$level.name, sep = "")
        )
      )
    )
  }

  # Clean up labels
  yLabels$labels <- gsub("_", " ", yLabels$labels)
  yLabels$labels <- gsub(": NA", "", yLabels$labels)
  yLabels <- yLabels[!is.na(yLabels$y.pos), ]

  # Set x values for plotting
  tab$x.val <- ifelse(tab$estimate.label == "(Reference)", 1, tab$estimate)

  # Handle extreme confidence intervals with xlim
  tab$ci_low_arrow <- FALSE
  tab$ci_high_arrow <- FALSE
  tab$show_point <- TRUE  # Flag to control point visibility

  if (!is.null(xlim) && length(xlim) == 2) {
    # Identify CIs outside limits OR NA for arrow indicators ----
    # NA bounds indicate unstable CIs (e.g., perfect separation) - show arrow in that direction
    tab$ci_low_arrow <- is.na(tab$conf.low) | (!is.na(tab$conf.low) & tab$conf.low < xlim[1])
    tab$ci_high_arrow <- is.na(tab$conf.high) | (!is.na(tab$conf.high) & tab$conf.high > xlim[2])

    # Clip CI bounds to limits ----
    tab$conf.low_plot <- pmax(tab$conf.low, xlim[1], na.rm = TRUE)
    tab$conf.high_plot <- pmin(tab$conf.high, xlim[2], na.rm = TRUE)

    # Don't show point if estimate is outside limits ----
    tab$show_point <- (tab$x.val >= xlim[1] & tab$x.val <= xlim[2]) |
      (tab$estimate.label == "(Reference)")

    # Set x position for points (only used when show_point = TRUE) ----
    tab$x.val_plot <- tab$x.val
  } else {
    # No clipping ----
    tab$conf.low_plot <- tab$conf.low
    tab$conf.high_plot <- tab$conf.high
    tab$x.val_plot <- tab$x.val
  }

  # Remove CIs from reference levels ----
  is_reference <- !is.na(tab$ref) & tab$ref == TRUE
  tab$conf.low_plot[is_reference] <- NA
  tab$conf.high_plot[is_reference] <- NA
  tab$ci_low_arrow[is_reference] <- FALSE
  tab$ci_high_arrow[is_reference] <- FALSE

  # Don't show arrows if estimate itself is NA ----
  tab$ci_low_arrow[is.na(tab$estimate)] <- FALSE
  tab$ci_high_arrow[is.na(tab$estimate)] <- FALSE

  # Assign colours based on risk
  tab$colour <- ifelse(
    tab$x.val < 1, "a",
    ifelse(tab$x.val == 1, "b", "c")
  )

  # Set colour palette
  if (length(colours) == 1) {
    colours <- c(a = "#006B3C", b = "black", c = "#FF0800")
  } else {
    names(colours) <- c("a", "b", "c")
  }
  colours <- colours[sort(unique(tab$colour))]

  # Create secondary axis for N and Event
  if (showN & !showEvent) {
    Axis <- scale_y_continuous(
      breaks = yLabels$y.pos,
      labels = yLabels$labels,
      sec.axis = dup_axis(
        breaks = yLabels$y.pos,
        labels = tab$N,
        name = "N"
      )
    )
    themeSecAxis <- theme(
      axis.title.y.right = element_text(angle = 0, hjust = 0.5, vjust = 0.5)
    )
  } else if (showEvent) {
    # Check if Event column exists
    if ("Event" %in% names(tab)) {
      Axis <- scale_y_continuous(
        breaks = yLabels$y.pos,
        labels = yLabels$labels,
        sec.axis = dup_axis(
          breaks = yLabels$y.pos,
          labels = paste(tab$N, tab$Event, sep = " : "),
          name = "N : Event"
        )
      )
      themeSecAxis <- theme(
        axis.title.y.right = element_text(angle = 270, hjust = 0.5, vjust = 0.5)
      )
    } else {
      # Fall back to just N if Event not available
      Axis <- scale_y_continuous(
        breaks = yLabels$y.pos,
        labels = yLabels$labels,
        sec.axis = dup_axis(
          breaks = yLabels$y.pos,
          labels = tab$N,
          name = "N"
        )
      )
      themeSecAxis <- theme(
        axis.title.y.right = element_text(angle = 0, hjust = 0.5, vjust = 0.5)
      )
    }
  } else {
    Axis <- scale_y_continuous(
      breaks = yLabels$y.pos,
      labels = yLabels$labels
    )
    themeSecAxis <- NULL
  }

  # Create base plot
  if (include_unadjusted) {
    # Plot with shape and linetype for UV vs MV
    p <- ggplot(tab, aes(x = x.val_plot, y = y.val, colour = colour,
                         linetype = type, shape = type)) +
      geom_point(data = subset(tab, show_point), na.rm = TRUE, size = 2) +
      geom_errorbarh(
        aes(xmin = conf.low_plot, xmax = conf.high_plot),
        height = 0, size = 0.9, na.rm = TRUE
      ) +
      # Add arrows for CIs extending beyond limits
      geom_segment(
        data = subset(tab, ci_low_arrow),
        aes(x = xlim[1], xend = xlim[1] * 0.95, y = y.val, yend = y.val,
            colour = colour, linetype = type),
        arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
        size = 0.9, na.rm = TRUE, show.legend = FALSE
      ) +
      geom_segment(
        data = subset(tab, ci_high_arrow),
        aes(x = xlim[2], xend = xlim[2] * 1.05, y = y.val, yend = y.val,
            colour = colour, linetype = type),
        arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
        size = 0.9, na.rm = TRUE, show.legend = FALSE
      ) +
      scale_shape_manual(values = c(15, 16)) +
      scale_linetype_manual(values = c("solid", "dashed")) +
      guides(
        shape = guide_legend(title = NULL, nrow = 1),
        linetype = guide_legend(title = NULL, nrow = 1)
      )
  } else {
    # Plot without shape/linetype
    p <- ggplot(tab, aes(x = x.val_plot, y = y.val, colour = colour)) +
      geom_point(data = subset(tab, show_point), na.rm = TRUE, size = 2) +
      geom_errorbarh(
        aes(xmin = conf.low_plot, xmax = conf.high_plot),
        height = 0, size = 0.9, na.rm = TRUE
      ) +
      # Add arrows for CIs extending beyond limits
      geom_segment(
        data = subset(tab, ci_low_arrow),
        aes(x = xlim[1], xend = xlim[1] * 0.95, y = y.val, yend = y.val, colour = colour),
        arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
        size = 0.9, na.rm = TRUE, show.legend = FALSE
      ) +
      geom_segment(
        data = subset(tab, ci_high_arrow),
        aes(x = xlim[2], xend = xlim[2] * 1.05, y = y.val, yend = y.val, colour = colour),
        arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
        size = 0.9, na.rm = TRUE, show.legend = FALSE
      )
  }

  # Add common plot elements
  p <- p +
    geom_vline(xintercept = 1) +
    labs(y = "", x = x_lab) +
    guides(colour = "none") +
    Axis +
    scale_colour_manual(values = colours) +
    theme_bw() +
    theme(
      axis.text.y = element_text(
        face = ifelse(tab$is_first_level, "bold", "plain"),
        hjust = 0
      ),
      legend.title = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.ticks = element_blank()
    )

  # Add secondary axis theme if needed
  if (!is.null(themeSecAxis)) {
    p <- p + themeSecAxis
  }

  # Add log scale if requested
  if (logScale) {
    # Use xlim if provided, otherwise calculate from data
    if (!is.null(xlim) && length(xlim) == 2) {
      x_min <- xlim[1]
      x_max <- xlim[2]
    } else {
      # Calculate appropriate limits to show all data
      # Handle infinite or extreme CI bounds
      finite_low <- tab$conf.low[is.finite(tab$conf.low)]
      finite_high <- tab$conf.high[is.finite(tab$conf.high)]

      # If we have finite values, use them
      if (length(finite_low) > 0 && length(finite_high) > 0) {
        x_min <- min(finite_low, na.rm = TRUE)
        x_max <- max(finite_high, na.rm = TRUE)

        # Cap extreme values to keep plot readable
        # If max is more than 100x the min, cap it
        if (x_max / x_min > 100) {
          x_max <- x_min * 100
          message("Note: Very wide confidence intervals detected. X-axis capped for readability.")
        }
      } else {
        # Fallback to estimate range if CIs are all infinite
        finite_est <- tab$estimate[is.finite(tab$estimate)]
        if (length(finite_est) > 0) {
          x_min <- min(finite_est, na.rm = TRUE) * 0.5
          x_max <- max(finite_est, na.rm = TRUE) * 2
        } else {
          # Last resort: default range
          x_min <- 0.1
          x_max <- 10
        }
      }
    }

    # Expand limits slightly to ensure all labels show
    x_range <- log10(x_max) - log10(x_min)
    x_expand <- x_range * 0.1

    p <- p +
      scale_x_log10(
        breaks = scales::log_breaks(n = nxTicks),
        limits = c(10^(log10(x_min) - x_expand), 10^(log10(x_max) + x_expand)),
        expand = expansion(mult = 0.05),
        oob = scales::squish  # Squish out-of-bounds values instead of removing
      )
  }

  return(p)
}


#' Create an univariable forest plot using ggplot2
#'
#' This function creates forest plots from univariable regression models.
#' It is maintained for backwards compatibility. For new code, consider using
#' forestplotMV() which can now handle both adjusted and unadjusted estimates.
#'
#' @param response character vector with names of columns to use for response
#' @param covs character vector with names of columns to use for covariates
#' @param data dataframe containing your data
#' @param model fitted model object (default "glm")
#' @param id character vector which identifies clusters. Only used for geeglm
#' @param corstr character string specifying the correlation structure. Only
#'   used for geeglm
#' @param family description of the error distribution and link function to be
#'   used in the model
#' @param digits number of digits to round to (default 2)
#' @param conf.level controls the width of the confidence interval (default 0.95)
#' @param colours can specify colours for risks less than, equal to, and greater
#'   than 1.0. Default is green, black, red
#' @param showEst logical, should the risks be displayed on the plot in text?
#'   Default is TRUE
#' @param showRef logical, should reference levels be shown? Default is TRUE
#' @param logScale logical, should OR/RR be shown on log scale? Defaults to TRUE
#' @param nxTicks Number of tick marks for x-axis (default 5)
#' @param showN Show number of observations per variable and category (default TRUE)
#' @param showEvent Show number of events per variable and category (default TRUE)
#' @param xlim numeric vector of length 2 specifying x-axis limits. Default is
#'   c(0.2, 5). Confidence intervals extending beyond these limits will be
#'   shown with arrows. Set to NULL to disable limit constraining.
#' @import ggplot2
#' @importFrom scales log_breaks
#' @keywords plot
#' @return a ggplot object
#' @export
#' @examples
#' data("pembrolizumab")
#' forestplotUV(response = "orr",
#'              covs = c("change_ctdna_group", "sex", "age", "l_size"),
#'              data = pembrolizumab, family = 'binomial')
forestplotUV <- function(response, covs, data, model = "glm", id = NULL,
                         corstr = NULL, family = NULL,
                         digits = getOption("reportRmd.digits", 2),
                         conf.level = 0.95, colours = "default",
                         showEst = TRUE, showRef = TRUE,
                         logScale = getOption("reportRmd.logScale", TRUE),
                         nxTicks = 5, showN = TRUE, showEvent = TRUE,
                         xlim = c(0.2, 5)) {

  # Determine axis label
  if (inherits(model, "glm") || model == "glm") {
    if (!is.null(family)) {
      # Handle both string and family object inputs
      if (is.character(family)) {
        # Family passed as string
        if (family == "binomial") {
          x_lab <- "Unadjusted Odds Ratio"
        } else {
          x_lab <- "Unadjusted Odds Ratio"
        }
      } else if (is.function(family) || is.list(family)) {
        # Family passed as function or family object
        fam_obj <- if (is.function(family)) family() else family
        if (fam_obj$link == "logit") {
          x_lab <- "Unadjusted Odds Ratio"
        } else if (fam_obj$link == "log") {
          x_lab <- "Unadjusted Relative Risk"
        } else {
          x_lab <- "Unadjusted Odds Ratio"
        }
      } else {
        x_lab <- "Unadjusted Odds Ratio"
      }
    } else {
      x_lab <- "Unadjusted Odds Ratio"
    }
  } else {
    x_lab <- "Unadjusted Odds Ratio"
  }

  # Fit UV models for each covariate
  uv_tabs <- list()

  for (cov in covs) {
    # Create formula
    uv_formula <- as.formula(paste(response, "~", cov))

    # Fit model
    if (model == "glm" || inherits(model, "glm")) {
      if (!is.null(family)) {
        uv_model <- glm(uv_formula, data = data, family = family)
      } else {
        uv_model <- glm(uv_formula, data = data)
      }
    } else {
      stop("Currently only glm models are supported")
    }

    # Get summary
    uv_tab <- m_summary(uv_model, CIwidth = conf.level, digits = digits,
                        whichp = "levels", for_plot = TRUE)
    uv_tabs[[cov]] <- prepare_forest_data(uv_tab, model_type = "Unadjusted",
                                          digits = digits)
  }

  # Combine all UV models
  tab <- do.call(rbind, uv_tabs)

  # Create a pseudo-model for ordering purposes
  # We'll use the first UV model as template
  pseudo_formula <- as.formula(paste(response, "~", paste(covs, collapse = " + ")))
  if (!is.null(family)) {
    pseudo_model <- glm(pseudo_formula, data = data, family = family)
  } else {
    pseudo_model <- glm(pseudo_formula, data = data)
  }

  # Order data by risk
  tab <- order_forest_data(tab, pseudo_model, data)

  # Remove reference levels if requested
  if (!showRef) {
    tab <- tab[is.na(tab$ref) | tab$ref == FALSE, ]
  }

  # Create y position values with tighter within-variable spacing
  # This makes levels within a variable "hang together" visually
  tab$y.val <- NA
  current_y <- nrow(tab)

  for (i in 1:nrow(tab)) {
    tab$y.val[i] <- current_y

    # Check if next row is a different variable (create larger gap)
    if (i < nrow(tab)) {
      if (tab$variable[i] != tab$variable[i + 1]) {
        # Larger gap between variables
        current_y <- current_y - 1.5
      } else {
        # Smaller gap within variable
        current_y <- current_y - 0.6
      }
    }
  }

  yvals <- tab$y.val

  # Format estimate labels
  tab$estimate.label[is.na(tab$estimate.label)] <- ""
  tab$estimate.label[tab$estimate.label == "Reference"] <- "(Reference)"

  # Create y-axis labels with variable name on first level
  if (showEst) {
    yLabels <- data.frame(
      y.pos = yvals,
      labels = ifelse(
        tab$is_first_level,
        paste(tab$variable, ": ", tab$level.name, ": ", tab$estimate.label, sep = ""),
        paste("   ", tab$level.name, ": ", tab$estimate.label, sep = "")
      )
    )
  } else {
    yLabels <- data.frame(
      y.pos = yvals,
      labels = ifelse(
        tab$is_first_level,
        paste(tab$variable, ": ", tab$level.name, sep = ""),
        paste("   ", tab$level.name, sep = "")
      )
    )
  }

  # Clean up labels
  yLabels$labels <- gsub("_", " ", yLabels$labels)
  yLabels$labels <- gsub(": NA", "", yLabels$labels)
  yLabels <- yLabels[!is.na(yLabels$y.pos), ]

  # Set x values for plotting
  tab$x.val <- ifelse(tab$estimate.label == "(Reference)", 1, tab$estimate)

  # Handle extreme confidence intervals with xlim
  tab$ci_low_arrow <- FALSE
  tab$ci_high_arrow <- FALSE
  tab$show_point <- TRUE  # Flag to control point visibility

  if (!is.null(xlim) && length(xlim) == 2) {
    # Identify CIs outside limits OR NA for arrow indicators ----
    # NA bounds indicate unstable CIs (e.g., perfect separation) - show arrow in that direction
    tab$ci_low_arrow <- is.na(tab$conf.low) | (!is.na(tab$conf.low) & tab$conf.low < xlim[1])
    tab$ci_high_arrow <- is.na(tab$conf.high) | (!is.na(tab$conf.high) & tab$conf.high > xlim[2])

    # Clip CI bounds to limits ----
    tab$conf.low_plot <- pmax(tab$conf.low, xlim[1], na.rm = TRUE)
    tab$conf.high_plot <- pmin(tab$conf.high, xlim[2], na.rm = TRUE)

    # Don't show point if estimate is outside limits ----
    tab$show_point <- (tab$x.val >= xlim[1] & tab$x.val <= xlim[2]) |
      (tab$estimate.label == "(Reference)")

    # Set x position for points (only used when show_point = TRUE) ----
    tab$x.val_plot <- tab$x.val
  } else {
    # No clipping ----
    tab$conf.low_plot <- tab$conf.low
    tab$conf.high_plot <- tab$conf.high
    tab$x.val_plot <- tab$x.val
  }

  # Remove CIs from reference levels ----
  is_reference <- !is.na(tab$ref) & tab$ref == TRUE
  tab$conf.low_plot[is_reference] <- NA
  tab$conf.high_plot[is_reference] <- NA
  tab$ci_low_arrow[is_reference] <- FALSE
  tab$ci_high_arrow[is_reference] <- FALSE

  # Don't show arrows if estimate itself is NA ----
  tab$ci_low_arrow[is.na(tab$estimate)] <- FALSE
  tab$ci_high_arrow[is.na(tab$estimate)] <- FALSE

  # Assign colours based on risk
  tab$colour <- ifelse(
    tab$x.val < 1, "a",
    ifelse(tab$x.val == 1, "b", "c")
  )

  # Set colour palette
  if (length(colours) == 1) {
    colours <- c(a = "#006B3C", b = "black", c = "#FF0800")
  } else {
    names(colours) <- c("a", "b", "c")
  }
  colours <- colours[sort(unique(tab$colour))]

  # Create secondary axis for N and Event
  if (showN & !showEvent) {
    Axis <- scale_y_continuous(
      breaks = yLabels$y.pos,
      labels = yLabels$labels,
      sec.axis = dup_axis(
        breaks = yLabels$y.pos,
        labels = tab$N,
        name = "N"
      )
    )
    themeSecAxis <- theme(
      axis.title.y.right = element_text(angle = 0, hjust = 0.5, vjust = 0.5)
    )
  } else if (showEvent) {
    # Check if Event column exists
    if ("Event" %in% names(tab)) {
      Axis <- scale_y_continuous(
        breaks = yLabels$y.pos,
        labels = yLabels$labels,
        sec.axis = dup_axis(
          breaks = yLabels$y.pos,
          labels = paste(tab$N, tab$Event, sep = " : "),
          name = "N : Event"
        )
      )
      themeSecAxis <- theme(
        axis.title.y.right = element_text(angle = 270, hjust = 0.5, vjust = 0.5)
      )
    } else {
      # Fall back to just N if Event not available
      Axis <- scale_y_continuous(
        breaks = yLabels$y.pos,
        labels = yLabels$labels,
        sec.axis = dup_axis(
          breaks = yLabels$y.pos,
          labels = tab$N,
          name = "N"
        )
      )
      themeSecAxis <- theme(
        axis.title.y.right = element_text(angle = 0, hjust = 0.5, vjust = 0.5)
      )
    }
  } else {
    Axis <- scale_y_continuous(
      breaks = yLabels$y.pos,
      labels = yLabels$labels
    )
    themeSecAxis <- NULL
  }

  # Create plot
  p <- ggplot(tab, aes(x = x.val_plot, y = y.val, colour = colour)) +
    geom_point(data = subset(tab, show_point), na.rm = TRUE, size = 2) +
    geom_errorbarh(
      aes(xmin = conf.low_plot, xmax = conf.high_plot),
      height = 0, size = 0.9, na.rm = TRUE
    ) +
    # Add arrows for CIs extending beyond limits
    geom_segment(
      data = subset(tab, ci_low_arrow),
      aes(x = xlim[1], xend = xlim[1] * 0.95, y = y.val, yend = y.val, colour = colour),
      arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
      size = 0.9, na.rm = TRUE, show.legend = FALSE
    ) +
    geom_segment(
      data = subset(tab, ci_high_arrow),
      aes(x = xlim[2], xend = xlim[2] * 1.05, y = y.val, yend = y.val, colour = colour),
      arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
      size = 0.9, na.rm = TRUE, show.legend = FALSE
    ) +
    geom_vline(xintercept = 1) +
    labs(y = "", x = x_lab) +
    guides(colour = "none") +
    Axis +
    scale_colour_manual(values = colours) +
    theme_bw() +
    theme(
      axis.text.y = element_text(
        face = ifelse(tab$is_first_level, "bold", "plain"),
        hjust = 0
      ),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      axis.ticks = element_blank()
    )

  # Add secondary axis theme if needed
  if (!is.null(themeSecAxis)) {
    p <- p + themeSecAxis
  }

  # Add log scale if requested
  if (logScale) {
    # Use xlim if provided, otherwise calculate from data
    if (!is.null(xlim) && length(xlim) == 2) {
      x_min <- xlim[1]
      x_max <- xlim[2]
    } else {
      # Calculate appropriate limits to show all data
      # Handle infinite or extreme CI bounds
      finite_low <- tab$conf.low[is.finite(tab$conf.low)]
      finite_high <- tab$conf.high[is.finite(tab$conf.high)]

      # If we have finite values, use them
      if (length(finite_low) > 0 && length(finite_high) > 0) {
        x_min <- min(finite_low, na.rm = TRUE)
        x_max <- max(finite_high, na.rm = TRUE)

        # Cap extreme values to keep plot readable
        # If max is more than 100x the min, cap it
        if (x_max / x_min > 100) {
          x_max <- x_min * 100
          message("Note: Very wide confidence intervals detected. X-axis capped for readability.")
        }
      } else {
        # Fallback to estimate range if CIs are all infinite
        finite_est <- tab$estimate[is.finite(tab$estimate)]
        if (length(finite_est) > 0) {
          x_min <- min(finite_est, na.rm = TRUE) * 0.5
          x_max <- max(finite_est, na.rm = TRUE) * 2
        } else {
          # Last resort: default range
          x_min <- 0.1
          x_max <- 10
        }
      }
    }

    # Expand limits slightly to ensure all labels show
    x_range <- log10(x_max) - log10(x_min)
    x_expand <- x_range * 0.1

    p <- p +
      scale_x_log10(
        breaks = scales::log_breaks(n = nxTicks),
        limits = c(10^(log10(x_min) - x_expand), 10^(log10(x_max) + x_expand)),
        expand = expansion(mult = 0.05),
        oob = scales::squish  # Squish out-of-bounds values instead of removing
      )
  }

  return(p)
}


#' Combine univariable and multivariable forest plot (DEPRECATED)
#'
#' This function is deprecated. Please use forestplotMV() with
#' include_unadjusted = TRUE instead.
#'
#' @param UVmodel an UV model object output from the forestplotUV function
#' @param MVmodel a MV model object output from the forestplotMV function
#' @param ... additional arguments (ignored)
#' @export
forestplotUVMV <- function(UVmodel, MVmodel, ...) {
  .Deprecated(
    "forestplotMV",
    msg = paste(
      "forestplotUVMV is deprecated.",
      "Please use forestplotMV(model, data, include_unadjusted = TRUE) instead."
    )
  )

  # Return the MV plot for backwards compatibility
  return(MVmodel)
}



#'Plot multiple bivariate relationships in a single plot
#'
#'This function is designed to accompany \code{\link{rm_uvsum}} as a means of
#'visualising the results, and uses similar syntax.
#'
#'Plots are displayed as follows: If response is continuous For a numeric
#'predictor scatterplot For a categorical predictor: If 20+ observations
#'available boxplot, otherwise dotplot with median line If response is a factor
#'For a numeric predictor: If 20+ observations available boxplot, otherwise
#'dotplot with median line For a categorical predictor barplot Response
#'variables are shown on the ordinate (y-axis) and covariates on the abscissa
#'(x-axis)
#'
#'Variable names are replaced by their labels if available, or by tidy versions
#'if not. Set use_labels=FALSE to use the variable names.
#'
#'@param response character vector with names of columns to use for response
#'@param covs character vector with names of columns to use for covariates
#'@param data dataframe containing your data
#'@param showN boolean indicating whether sample sizes should be shown on the
#'  plots
#'@param showPoints boolean indicating whether individual data points should be
#'  shown when n>20 in a category
#'@param na.rm boolean indicating whether na values should be shown or removed
#'@param response_title character value with title of the plot
#'@param return_plotlist boolean indicating that the list of plots should be
#'  returned instead of a plot, useful for applying changes to the plot, see
#'  details
#'@param ncol the number of columns of plots to be display in the ggarrange
#'  call, defaults to 2
#'@param p_margins sets the TRBL margins of the individual plots, defaults to
#'  c(0,0.2,1,.2)
#'@param bpThreshold Default is 20, if there are fewer than 20 observations in a
#'  category then dotplots, as opposed to boxplots are shown.
#'@param mixed should a mix of dotplots and boxplots be shown based on sample
#'  size? If false then all categories will be shown as either dotplots, or
#'  boxplots according the bpThreshold and the smallest category size
#'@param violin Show violin plots instead of boxplots. This will override
#'  bpThreshold and mixed.
#'@param position for categorical variables how should barplots be presented.
#'  Default is "dodge" IF stack is TRUE then n will not be shown.
#'@param use_labels boolean, default is true if the variables have label
#'  attributes this will be shown in the plot instead of the variable names, or
#'  if there are no labels then tidy versions of the variable names will be
#'  used. If use_labels=FALSE the variable names will be used.
#'@keywords plot
#'@returns a list containing plots for each variable in covs
#'@importFrom ggplot2 ggplot aes_string geom_boxplot geom_point geom_text
#'  stat_summary scale_x_discrete stat theme labs .data
#'@importFrom ggpubr ggarrange
#'@importFrom stats median
#'@export
#' @examples
#' ## Run multiple univariate analyses on the pembrolizumab dataset to predict cbr and
#' ## then visualise the relationships.
#' data("pembrolizumab")
#' rm_uvsum(data=pembrolizumab,
#' response='cbr',covs=c('age','sex','l_size','baseline_ctdna'))
#' plotuv(data=pembrolizumab,  response='cbr',
#' covs=c('age','sex','l_size','baseline_ctdna'),showN=TRUE)
#'@seealso \code{\link[ggplot2:ggplot]{ggplot2::ggplot}} and
#'  \code{\link[ggpubr:ggarrange]{ggpubr::ggarrange}}
#'  \code{\link{replace_plot_labels}}
plotuv <- function(response,covs,data,showN=FALSE,showPoints=TRUE,na.rm=TRUE,
                   response_title=NULL,return_plotlist=FALSE,ncol=2,p_margins=c(0,0.2,1,.2),
                   bpThreshold=20,mixed=TRUE,violin=FALSE,position=c("dodge","stack","fill"),use_labels=TRUE) {

  if (missing(response) & !is.null(response_title)) {
    warning("response_title will be ignored because no response variable was provided")
  }
  if (violin){
    showPoints=FALSE
  }
  if (missing(response)) {
    for (v in covs){
      if (!v %in% names(data)) stop(paste(v,'is not a variable in data.'))
      if (inherits(data[[v]],'character')) data[[v]] <- factor(data[[v]])
    }
    plist <- NULL

    for (x_var in covs) {
      if (na.rm) {
        pdata <- stats::na.omit(data[, x_var])
      }
      else {
        pdata <- data[, x_var]
      }
      if (inherits(pdata[[x_var]], "numeric")) { # if x_var is numeric
        p <- ggplot(data = pdata, aes(x = .data[[x_var]]))
        if (violin) {
          p <- p + geom_violin()
        }
        else {
          p <- p + geom_boxplot(outlier.shape = 16, outlier.colour = "red", na.rm = TRUE)
        }
      }
      else {  # x_var is categorical
        p <- ggplot(data = pdata, aes(x = .data[[x_var]])) + geom_bar(fill = "white", colour = "black") + scale_x_discrete(labels = function(x) wrp_lbl(x))
        if (showN){
          p <- p + geom_text(aes(label = stat(count)), position = position_dodge(width = 1), stat = 'count', vjust = 2, size = 2.5)
        }
        if (length(unique(pdata[[x_var]]))>8){
          p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
        }
      }
      plist[[x_var]] <- p + labs(x=niceStr(x_var),y='',fill=response_title)
    }
    # replace variable names with variable labels
    if (use_labels){
      plist <- lapply(plist,function(x) replace_plot_labels(x))
    }

    if (return_plotlist){
      return(plist)
    }
    else {
      suppressMessages(ggpubr::ggarrange(plotlist=plist,
                                         ncol=ncol,
                                         nrow=ceiling(length(plist)/ncol)))
    }
  }
  else {
    for (v in c(response,covs)){
      if (!v %in% names(data)) stop(paste(v,'is not a variable in data.'))
      if (inherits(data[[v]],'character')) data[[v]] <- factor(data[[v]])
    }
    position = match.arg(position)
    if (position=="stack"){
      bar_position=ggplot2::position_stack
      showN=FALSE
    } else if (position=="dodge"){
      bar_position=ggplot2::position_dodge
    } else if (position=="fill"){
      bar_position=ggplot2::position_fill
      showN=FALSE
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
            p <- ggplot(data=pdata, aes(y=.data[[response]],x=.data[[x_var]],fill=.data[[response]]))
            if (violin) {
              p <- p + geom_violin(aes(alpha=.data[['alpha']],linetype=.data[['lty']]))
            } else{
              p <- p + geom_boxplot(aes(alpha=.data[['alpha']],linetype=.data[['lty']]),outlier.shape = NA)
            }
            p <- p +
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
            geom_bar(position=bar_position()) +
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
            p <- ggplot(data=pdata, aes(x=.data[[x_var]],y=.data[[response]],fill=.data[[x_var]]))
            if (violin) {
              p <- p + geom_violin(aes(alpha=.data[['alpha']],linetype=.data[['lty']]))
            } else{
              p <- p + geom_boxplot(aes(alpha=.data[['alpha']],linetype=.data[['lty']]),outlier.shape = NA)
            }
            p <- p + scale_alpha_manual(breaks=c('light','regular'),values=c(0,1)) +
              scale_linetype_manual(breaks=c('0','1'),values = c(0,1))+
              scale_x_discrete(labels= function(x) wrp_lbl(x))
            p <- p +geom_point(data=coloured_points,aes(colour=.data[[x_var]]),
                               position=position_jitterdodge(),alpha=0.9)

            if (showPoints) {
              p <- p +geom_jitter(data= black_points,
                                  color="black", size=0.4, alpha=0.9)
            }
            if (showN){
              p <-  p+
                stat_summary(aes(y=max(.data[[response]])*1.05,vjust=0),geom='label',fun.data = lbl_count,label.size=0,fill='white',label.padding = unit(0.15, "lines"),alpha=.8)+
                geom_blank(aes( y=max(.data[[response]])*1.1))

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
    # replace variable names with variable labels
    if (use_labels){
      plist <- lapply(plist,function(x) replace_plot_labels(x))
    }
    if (return_plotlist){
      return(plist)
    } else{   suppressMessages(ggpubr::ggarrange(plotlist=plist,
                                                 common.legend = use_common_legend,
                                                 ncol=ncol,
                                                 nrow=ceiling(length(plist)/ncol),
                                                 legend.grob = legend.grob))
    }
  }
}


# Rmarkdown Reporting ----
#' Print tables to PDF/Latex HTML or Word
#'
#' Output the table nicely to whatever format is appropriate. This is the output
#' function used by the rm_* printing functions.
#'
#' Entire rows can be bolded, or specific cells. Currently indentation refers to
#' the first column only. By default, underscores in column names are converted
#' to spaces. To disable this set nicenames to FALSE
#'
#' @param tab a table to format
#' @param row.names a string specifying the column name to assign to the
#'   rownames. If NULL (the default) then rownames are removed.
#' @param to_indent numeric vector indicating which rows to indent in the first
#'   column.
#' @param bold_headers boolean indicating if the column headers should be
#'   bolded
#' @param rows_bold numeric vector indicating which rows to bold
#' @param bold_cells array indices indicating which cells to bold. These will be
#'   in addition to rows bolded by rows_bold.
#' @param caption table caption
#' @param digits number of digits to round numeric columns to, either a single
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
#' @param nicenames boolean indicating if you want to replace . and _ in strings
#'   with a space
#' @param fontsize PDF/HTML output only, manually set the table fontsize
#' @param chunk_label only used knitting to Word docs to allow cross-referencing
#' @param format if specified ('html','latex') will override the
#'   global pandoc setting
#' @return A character vector of the table source code, unless tableOnly=TRUE in
#'   which case a data frame is returned
#' @export
outTable <- function(tab, row.names = NULL, to_indent = numeric(0), bold_headers = TRUE,
                     rows_bold = numeric(0), bold_cells = NULL, caption = NULL,
                     digits = getOption("reportRmd.digits", 2), align,
                     applyAttributes = TRUE, keep.rownames = FALSE, nicenames = TRUE,
                     fontsize, chunk_label, format = NULL) {

  # Input validation
  if (!inherits(tab, "data.frame")) stop("tab must be a data frame")
  if (nrow(tab) == 0) return(NULL)

  # Strip tibble aspects and convert to data.frame
  tab <- as.data.frame(tab)

  # Handle row names
  if (!is.null(row.names)) {
    tab <- cbind(rownames(tab), tab)
    names(tab)[1] <- row.names
  }
  rownames(tab) <- NULL

  # Store original column names before any processing
  original_names <- names(tab)

  # Validate bold_cells specification
  if (!is.null(bold_cells)) {
    if (nrow(bold_cells) == 0) {
      bold_cells <- NULL
    } else {
      if (max(bold_cells[, 1]) > nrow(tab) | max(bold_cells[, 2]) > ncol(tab)) {
        message("Cell bolding incorrectly specified, no cells will be bolded")
        bold_cells <- NULL
      }
    }
  }

  # Define column alignment
  if (missing(align)) {
    alignSpec <- paste(c('l', rep('r', ncol(tab) - 1)), collapse = '')
  } else {
    alignSpec <- gsub('[^lrc]+', '', paste(align, collapse = ''))
    alignSpec <- substr(alignSpec, 1, ncol(tab))
    if (nchar(alignSpec) < ncol(tab)) {
      lastchar <- substr(alignSpec, nchar(alignSpec), nchar(alignSpec))
      alignSpec <- paste0(alignSpec, paste(rep(lastchar, ncol(tab) - nchar(alignSpec)), collapse = ''))
    }
    if (!identical(alignSpec, align)) {
      warning(paste0('Argument align did not conform to expectations, align="', alignSpec, '" used instead'))
    }
  }

  # Round and format numeric columns
  if (!missing(digits)) {
    numCols <- names(tab)[sapply(tab, inherits, 'numeric')]
    if (length(numCols) > 0) {
      colRound <- cbind(numCols, digits)
      colDigits <- as.numeric(colRound[, 2])
      names(colDigits) <- colRound[, 1]
      for (v in numCols) {
        tab[[v]] <- sapply(tab[[v]], function(x) niceNum(x, digits = colDigits[v]))
      }
    }
  }

  # Determine output format
  out_fmt <- if (!is.null(format) && format %in% c('html', 'latex')) {
    format
  } else if (is.null(knitr::pandoc_to())) {
    'html'
  } else if (knitr::pandoc_to(c('doc', 'docx'))) {
    'doc'
  } else if (knitr::is_latex_output()) {
    'latex'
  } else {
    'html'
  }

  # Handle chunk label
  chunk_label <- if (missing(chunk_label)) 'NOLABELTOADD' else chunk_label

  # Apply stored attributes if requested
  if (applyAttributes && !is.null(attr(tab, 'dimchk'))) {
    if (all(attr(tab, 'dimchk') == dim(tab))) {
      if (!is.null(attr(tab, 'to_indent'))) to_indent <- attr(tab, 'to_indent')
      if (!is.null(attr(tab, 'bold_cells'))) bold_cells <- attr(tab, 'bold_cells')
    }
  }

  # Ensure to_indent is a vector
  to_indent <- as.vector(to_indent)
  if (is.null(to_indent)) to_indent <- numeric(0)

  # Apply nice names if requested
  if (nicenames) {
    names(tab) <- nicename(names(tab))
  }

  # Combine rows_bold with bold_cells
  if (length(rows_bold) > 0) {
    arrInd <- as.matrix(expand.grid(rows_bold, 1:ncol(tab)))
    bold_cells <- rbind(bold_cells, arrInd)
    dimnames(bold_cells) <- NULL
    bold_cells <- bold_cells[!duplicated(bold_cells), , drop = FALSE]
  }

  # Clean up bold_cells
  if (!is.null(bold_cells)) {
    bold_cells <- bold_cells[!duplicated(bold_cells), , drop = FALSE]
    bold_cells <- bold_cells[!is.na(tab[bold_cells]), , drop = FALSE]
  }

  # Word/Doc output
  if (out_fmt == 'doc') {
    caption <- if (!is.null(caption)) {
      if (chunk_label == 'NOLABELTOADD') {
        caption
      } else {
        paste0('(\\#tab:', chunk_label, ')', caption)
      }
    }

    # Handle Date columns (convert to character to avoid errors)
    date_cols <- names(tab)[sapply(tab, inherits, "Date")]
    if (length(date_cols) > 0) {
      for (v in date_cols) {
        tab[[v]] <- as.character(tab[[v]])
      }
    }

    # Preserve column names before conversion
    col_names <- names(tab)

    # Convert all columns to character
    tab <- tab |>
      lapply(as.character) |>
      as.data.frame(stringsAsFactors = FALSE, check.names = FALSE)

    # Restore column names explicitly
    names(tab) <- col_names

    # Replace NA and empty strings with non-breaking space
    tab[is.na(tab)] <- '&nbsp;'
    tab[tab == ''] <- '&nbsp;'

    # Apply indentation
    if (length(to_indent) > 0) {
      tab[[1]][to_indent] <- paste('&nbsp;&nbsp;', tab[[1]][to_indent])
    }

    # Output with pander
    pander::pander(
      tab,
      caption = caption,
      emphasize.strong.cells = bold_cells,
      split.table = Inf,
      split.cells = 15,
      justify = alignSpec
    )

  } else {
    # PDF/HTML output

    # Set NA to empty in kable
    oldop <- options()
    on.exit(options(oldop))
    options(knitr.kable.NA = '')

    # Format-specific text processing
    if (out_fmt == 'latex') {
      names(tab) <- sanitize(names(tab))
      if (!is.null(caption)) caption <- sanitize(caption)
      for (v in 1:ncol(tab)) tab[[v]] <- sanitize(tab[[v]])
      if (!is.null(bold_cells)) tab[bold_cells] <- sapply(tab[bold_cells], function(x) lbld(x))
    }

    if (out_fmt == 'html') {
      names(tab) <- ltgt(names(tab))
      for (v in 1:ncol(tab)) tab[[v]] <- rmds(tab[[v]])
      if (!is.null(bold_cells)) tab[bold_cells] <- sapply(tab[bold_cells], function(x) hbld(x))
    }

    # Determine if long table needed
    long_table <- nrow(tab) > 30

    # Create kable output
    kout <- knitr::kable(
      tab,
      format = out_fmt,
      escape = FALSE,
      booktabs = TRUE,
      longtable = long_table,
      linesep = '',
      caption = caption,
      align = alignSpec
    )

    # Apply styling based on format
    if (out_fmt == "html") {
      kout <- kout |>
        kableExtra::kable_styling(fixed_thead = TRUE) |>
        kableExtra::kable_styling(full_width = TRUE)
    } else {
      kout <- kout |>
        kableExtra::kable_styling(latex_options = c('repeat_header'))
    }

    # Apply indentation
    kout <- kableExtra::add_indent(kout, positions = to_indent)

    # Apply font size if specified
    if (!missing(fontsize)) {
      kout <- kableExtra::kable_styling(kout, font_size = fontsize)
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
#' data(pembrolizumab)
#' fit1 <- lm(baseline_ctdna~age+l_size+pdl1,data=pembrolizumab)
#' m1 <- rm_mvsum(fit1,tableOnly=TRUE)
#' m1$Response = 'ctDNA'
#' fit2 <- lm(l_size~age+baseline_ctdna+pdl1,data=pembrolizumab)
#' m2 <- rm_mvsum(fit2,tableOnly=TRUE)
#' m2$Response = 'Tumour Size'
#' nestTable(rbind(m1,m2),head_col='Response',to_col='Covariate')
nestTable <- function(data,head_col,to_col,colHeader ='',caption=NULL,indent=TRUE,boldheaders=TRUE,hdr_prefix='',hdr_suffix='',digits=getOption("reportRmd.digits",2),tableOnly=FALSE,fontsize){

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

  data <- data[,setdiff(names(data),head_col),drop=FALSE]
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

#' Output a scrollable table
#'
#' This function accepts the output of a aa call to knitr::kable or
#' reportRmd::outTable and, if the output format is html, will produce a
#' scrollable table. Otherwise a regular table will be output for pandoc/latex
#'
#' @param knitrTable output from a call to knitr::kable or outTable
#' @param pixelHeight the height of the scroll box in pixels, default is 500
#' @importFrom kableExtra scroll_box
#' @examples
#' data("pembrolizumab")
#' tab <- rm_covsum(data=pembrolizumab,maincov = 'change_ctdna_group',
#' covs=c('age','cohort','sex','pdl1','tmb','l_size'),full=F)
#' scrolling_table(tab,pixelHeight=300)
#' @export
scrolling_table <- function(knitrTable,pixelHeight=500){

  out_fmt = ifelse(is.null(knitr::pandoc_to()),'html',
                   ifelse(knitr::pandoc_to(c('doc','docx')),'doc',
                          ifelse(knitr::is_latex_output(),'latex','html')))
  if (out_fmt %in% c('doc','latex')) return(knitrTable)

  if (!inherits(knitrTable,"knitr_kable")) stop("This function requires a knitr_kable object.\nTry running reportRmd::outTable prior to use.")

  kableExtra::scroll_box(
    knitrTable,
    height = paste0(pixelHeight,"px;"),
    box_css = "border: 1px solid #ddd; padding: 5px; ")
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
#' A newer version of this function is \link{rm_compactsum} which is more
#' flexible and displays fewer rows of output.
#'
#' @param data dataframe containing data
#' @param covs Character vector of covariate names to summarize.
#'   Can also be specified using the \code{xvars} alias.
#' @param maincov Character string specifying the grouping variable.
#'   Can also be specified using the \code{grp} alias.
#' @param xvars Alias for \code{covs}. Accepts either unquoted variable names
#'   (e.g., \code{c(age, sex)}) or character strings (e.g., \code{c("age", "sex")}).
#' @param grp Alias for \code{maincov}. Accepts either an unquoted variable name
#'   (e.g., \code{sex}) or a character string (e.g., \code{"sex"}).
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
#'   table. Can only be obtained if pvalue is also requested. Effect sizes
#'   calculated include Cramer's V for categorical variables, Cohen's d,
#'   Wilcoxon r, or Eta-squared for numeric/continuous variables.
#' @param p.adjust p-adjustments to be performed. Uses the
#'  [p.adjust] function from base R
#' @param unformattedp boolean indicating if you would like the p-value to be
#'   returned unformatted (ie not rounded or prefixed with '<'). Best used with
#'   tableOnly = T and outTable function. See examples.
#' @param show.tests boolean indicating if the type of statistical test and
#'   effect size used should be shown in a column beside the pvalues. Ignored if
#'   pvalue=FALSE.
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
#'   \code{\link{kruskal.test}}, \code{\link{anova}}, \code{\link[rstatix:cramer_v]{rstatix::cramer_v}},
#'   \code{\link[rstatix:eta_squared]{rstatix:eta_squared}}, and \code{\link{outTable}}
#' @references Ellis, P.D. (2010) The essential guide to effect sizes:
#' statistical power, meta-analysis, and the interpretation of research
#' results. Cambridge: Cambridge University Press.\doi{10.1017/CBO9780511761676}
#' @references Lakens, D. (2013)  Calculating and reporting effect sizes to
#' facilitate cumulative science: a practical primer for t-tests and ANOVAs.
#' Frontiers in Psychology, 4; 863:1-12. \doi{10.3389/fpsyg.2013.00863}
#' @examples
#' data("pembrolizumab")
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
rm_covsum <- function (data, covs=NULL, maincov = NULL, caption = NULL, tableOnly = FALSE,
                       covTitle = "", digits = 1, digits.cat = 0, nicenames = TRUE,
                       IQR = FALSE, all.stats = FALSE, pvalue = TRUE, effSize = FALSE,
                       p.adjust='none', unformattedp = FALSE,
                       show.tests = FALSE, testcont = c("rank-sum test", "ANOVA"),
                       testcat = c("Chi-squared", "Fisher"), full = TRUE, include_missing = FALSE,
                       percentage = c("column", "row"), dropLevels = TRUE, excludeLevels = NULL,
                       numobs = NULL, fontsize,chunk_label, xvars=NULL,grp=NULL)
{
  if (missing(data))
    stop("data is a required argument")
  if (!inherits(data, "data.frame"))
    stop("data must be supplied as a data frame.")

  # Handle grp alias for maincov ----
  if (!is.null(substitute(grp)) && is.null(maincov)) {
    grp_expr <- substitute(grp)
    if (is.character(grp)) {
      # Already a string ----
      maincov <- grp
    } else {
      # Convert symbol to string ----
      maincov <- deparse(grp_expr)
    }
  }

  # Handle xvars alias for covs ----
  if (!is.null(substitute(xvars)) && is.null(covs)) {
    xvars_expr <- substitute(xvars)

    if (is.character(xvars)) {
      # Already strings ----
      covs <- xvars
    } else if (is.call(xvars_expr) && identical(xvars_expr[[1]], as.name("c"))) {
      # Extract from c(...) - could be symbols or strings ----
      covs <- sapply(as.list(xvars_expr)[-1], function(x) {
        if (is.character(eval(x, parent.frame()))) {
          eval(x, parent.frame())
        } else {
          deparse(x)
        }
      })
    } else {
      # Single variable ----
      covs <- deparse(xvars_expr)
    }
  }

  # Check that we have required arguments ----
  if (is.null(covs)) {
    stop("Either 'covs' or 'xvars' must be provided")
  }

  argList <- as.list(match.call(expand.dots = TRUE)[-1])
  df_nm <- deparse(argList$data)
  argsToPass <- intersect(names(formals(covsum)), names(argList))
  covsumArgs <- argList[names(argList) %in% argsToPass]
  covsumArgs[["data"]] <- data  # required to avoid error if data object shares the name of an R function
  covsumArgs[["markup"]] <- FALSE
  covsumArgs[["sanitize"]] <- FALSE
  covsumArgs[["nicenames"]] <- FALSE
  covsumArgs[["covs"]] <- covs
  covsumArgs[["maincov"]] <- maincov

  tab <- do.call(covsum, covsumArgs)
  Sys.sleep(1)
  to_indent <- which(!attr(tab,"varID"))
  to_bold_name <- which(attr(tab,"varID"))
  bold_cells <- arrayInd(to_bold_name, dim(tab))

  if (nicenames) tab$Covariate <- replaceLbl(data, tab$Covariate)
  names(tab)[1] <- covTitle
  if ("p-value" %in% names(tab)) {
    if (p.adjust!='none'){
      #tab[["p (unadjusted)"]] <- tab[["p-value"]]
      unadjusted_p <- tab[["p-value"]]
      tab[["p-value"]] <- p.adjust(unadjusted_p,method=p.adjust)
    }
    to_bold_p <- which(as.numeric(tab[["p-value"]]) < 0.05)
    if (!unformattedp){
      #      if (!is.null(tab[["p (unadjusted)"]])) tab[["p (unadjusted)"]] <- sapply(tab[["p (unadjusted)"]],formatp)
      tab[["p-value"]] <- sapply(tab[["p-value"]],formatp)
    }
    if (length(to_bold_p) > 0)
      bold_cells <- rbind(bold_cells, matrix(cbind(to_bold_p,
                                                   which(names(tab) == "p-value")), ncol = 2))
    tab[["p-value"]] <- sapply(tab[["p-value"]],function(x) ifelse(nchar(x)==0,NA,x))
  }
  if ("Effect Size" %in% names(tab)) {
    e_vals <- tab[["Effect Size"]]
    new_e <- sapply(e_vals, formatp)
    tab[["Effect Size"]] <- new_e
  }
  if (tableOnly) {
    if (names(tab)[1] == "")
      names(tab)[1] <- "Covariate"
    attr(tab,"data") <- df_nm
    attr(tab,"data call") <- deparse1(argList$data)
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


#' Combine univariate and multivariable regression tables
#'
#' This function will combine rm_uvsum and rm_mvsum outputs into a single table.
#' The tableOnly argument must be set to TRUE when tables to be combined are
#' created. The resulting table will be in the same order as the uvsum table and
#' will contain the same columns as the uvsum and mvsum tables, but the p-values
#' will be combined into a single column. There must be a variable overlapping
#' between the uvsum and mvsum tables and all variables in the mvsum table must
#' also appear in the uvsum table.
#'
#'
#' @param uvsumTable Output from rm_uvsum, with tableOnly=TRUE
#' @param mvsumTable  Output from rm_mvsum, with tableOnly=TRUE
#' @param covTitle character with the names of the covariate (predictor) column.
#'   The default is to leave this empty for output or, for table only output to
#'   use the column name 'Covariate'.
#' @param vif boolean indicating if the variance inflation factor should be
#'   shown if present in the mvsumTable. Default is FALSE.
#' @param showN boolean indicating if sample sizes should be displayed.
#' @param showEvent boolean indicating if number of events (dichotomous
#'   outcomes) should be displayed.
#' @param caption table caption
#' @param tableOnly boolean indicating if unformatted table should be returned
#' @param chunk_label only used if output is to Word to allow cross-referencing
#' @param fontsize PDF/HTML output only, manually set the table fontsize
#' @seealso \code{\link{rm_uvsum}},\code{\link{rm_mvsum}}
#' @return A character vector of the table source code, unless tableOnly=TRUE in
#'   which case a data frame is returned
#' @export
#' @examples
#' require(survival)
#' data("pembrolizumab")
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
rm_uv_mv <- function(uvsumTable,mvsumTable,covTitle='',vif=FALSE,showN=FALSE, showEvent=FALSE,caption=NULL,tableOnly=FALSE,chunk_label,fontsize){
  # Check that tables are data frames and not kable objects
  if (!inherits(uvsumTable,'data.frame')) stop('uvsumTable must be a data.frame. Did you forget to specify tableOnly=TRUE?')
  if (!inherits(mvsumTable,'data.frame')) stop('mvsumTable must be a data.frame. Did you forget to specify tableOnly=TRUE?')
  # Check that the first columns have the same name
  if (names(uvsumTable)[1] != names(mvsumTable)[1]) stop('The covariate columns must have the same name in both tables')
  # check that variable label use is consistent
  if (is.null(attr(uvsumTable,"termlabels")) != is.null(attr(mvsumTable,"termlabels"))) stop("Both tables must either use variable labels or variable names\nRe-run summaries with the same data set.")
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
  if (showN){
    nc <- unlist(lapply(x, function(z){ length(which(names(z)=="N"))}))
    if (any(nc==0))    warning('To show sample size run rm_uvsum, rm_mvsum with showN=T')
  }
  x[[1]]$varOrder = 1:nrow(x[[1]])
  if (!showN){
    x <-lapply(x, function(z){
      nc <- which(names(z)=="N")
      if (length(nc>0)) z <- z[,-nc]
      return(z)
    })
  }
  if (!showEvent){
    x <-lapply(x, function(z){
      nc <- which(names(z)=="Event")
      if (length(nc>0)) z <- z[,-nc]
      return(z)
    })
  }
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



# Summary functions for plots --------------------------

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
#'   probability" for KM curves, and "Probability of an event" for CIF
#' @param main String corresponding to main title. When NULL uses Kaplan-Meier
#'   Plot s, and "Cumulative Incidence Plot for CIF"
#'
#' @param stratalabs string corresponding to the labels of the covariate, when
#'   NULL will use the levels of the covariate
#' @param strataname String of the covariate name default is nicename(cov)
#' @param stratalabs.table String corresponding to the levels of the covariate
#'   for the number at risk table, when NULL will use the levels of the
#'   covariate. Can use a string of "-" when the labels are long
#' @param strataname.table String of the covariate name for the number at risk
#'   table default is nicename(cov)
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
#' data("pembrolizumab")
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

  lifecycle::deprecate_warn(
    when = "0.1.2",
    what = "ggkmcif()",
    with = "ggkmcif2()",
    details = c(
      "Please update your code to use ggkmcif2()",
      "ggkmcif will be removed in future versions of reportRmd"
    )
  )

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
    # try dplyr::bind_rows

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


    if(returns) {
      if(table==TRUE){
        a <- list(p,blank.pic,data.table)
      } else a <- p

      return(a)
    } else {
      gridExtra::grid.arrange(gA, gB, gC,
                              clip = FALSE, nrow = 3, ncol = 1,
                              heights = unit(c(2, .1, .25), c("null", "null", "null")))

    }
  } else(return(p))


}

modify_ggkmcif <- function(list_gg){
  lifecycle::deprecate_warn("0.1.0","modify_ggkmcif","ggkmcif2()")
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
#' data("pembrolizumab")
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
  lifecycle::deprecate_warn("0.1.0","ggkmcif_paste()","ggkmcif2()")

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
#' difference in survival between groups in a single table suitable for
#' markdown output. Median survival times are included by default but can be
#' removed setting median=FALSE
#' @param data data frame containing survival data
#' @param time string indicating survival time variable
#' @param status string indicating event status variable
#' @param covs character vector indicating variables to group observations by
#' @param strata string indicating the variable to stratify observations by
#' @param includeVarNames boolean indicating if the variable names should be
#'   included in the output table, default is FALSE
#' @param digits the number of digits in the survival rate
#' @param showCols character vector indicating which of the optional columns
#'   to display, defaults to c('N','Observed','Expected')
#' @param CIwidth width of the median survival estimates, default is 95%
#' @param conf.type type of confidence interval see
#'   \code{\link[survival:survfit]{survival::survfit}} for details. Default is 'log'.
#' @param caption table caption
#' @param tableOnly should a dataframe or a formatted object be returned
#' @param fontsize PDF/HTML output only, manually set the table fontsize
#' @param unformattedp boolean indicating if you would like the p-value to be
#'   returned unformatted (ie not rounded or prefixed with '<'). Should be used
#'   in conjunction with the digits argument.
#' @importFrom  survival survdiff Surv strata
#' @seealso \code{\link[survival:survdiff]{survival::survdiff}}
#' @examples
#' #' # Differences between sex
#' data("pembrolizumab")
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
                        conf.type='log',caption=NULL,tableOnly=FALSE,fontsize,unformattedp=FALSE){
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
  if (unformattedp)
    formatp <- function(x, ...) {x}
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
#' @param conf.type type of confidence interval see \code{\link[survival:survfit]{survival::survfit}} for
#'   details. Default is 'log'.
#' @param na.action default is to omit missing values, but can be set to throw
#'   and error using na.action='na.fail'
#' @param showCounts boolean indicating if the at risk, events and censored
#'   columns should be output; default is TRUE
#' @param showLogrank boolean indicating if the log-rank test statistic and
#'   p-value should be output; default is TRUE
#' @param eventProb boolean indicating if event probabilities, rather than
#'   survival probabilities, should be displayed; default is FALSE
#' @param digits the number of digits in the survival rate, default is 2, unless
#'   the reportRmd.digits option is set
#' @param caption table caption for markdown output
#' @param tableOnly should a dataframe or a formatted object be returned
#' @param fontsize PDF/HTML output only, manually set the table fontsize
#' @importFrom  survival survfit Surv
#' @seealso \code{\link[survival:survfit]{survival::survfit}}
#' @return A character vector of the survival table source code, unless
#'   tableOnly=TRUE in which case a data frame is returned
#' @export
#' @examples
#' # Simple median survival table
#' data("pembrolizumab")
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
rm_survsum <- function (data, time, status, group = NULL, survtimes = NULL,
                        survtimeunit, survtimesLbls = NULL, CIwidth = 0.95, unformattedp = FALSE,
                        conf.type = "log", na.action = "na.omit", showCounts = TRUE, showLogrank = TRUE, eventProb = FALSE,
                        digits = getOption("reportRmd.digits",2), caption = NULL, tableOnly = FALSE,fontsize)
{
  if (missing(data))
    stop("data is a required argument")
  if (missing(time))
    stop("time is a required argument")
  if (missing(survtimeunit))
    if (!is.null(survtimes))
      stop("survtimeunit must be specified if survtimes are set. Example survtimeunit=\"year\"")
  if (!is.null(survtimes))
    timelbl <- paste0("Time (", survtimeunit, ")")
  if (!inherits(data, "data.frame"))
    stop("data must be supplied as a data frame")
  if (!inherits(time, "character") | length(time) > 1)
    stop("time must be supplied as a string indicating a variable in data")
  if (!inherits(status, "character") | length(status) > 1)
    stop("status must be supplied as a string indicating a variable in data")
  if (!is.null(survtimes))
    if (is.null(survtimesLbls))
      survtimesLbls <- survtimes
  if (length(survtimesLbls) != length(survtimes))
    stop("If supplied, the survtimesLbls vector must be the same length as survtime")
  if (unformattedp)
    formatp <- function(x, ...) {x}
  data <- data.frame(data)
  n_cols <- c("n.risk", "n.event", "n.censor")
  missing_vars = na.omit(setdiff(c(time, status, group), names(data)))
  if (length(missing_vars) > 0)
    stop(paste("These variables are not in the data:\n",
               paste0(missing_vars, collapse = csep())))
  sfit <- survival::survfit(as.formula(paste0("survival::Surv(",
                                              time, ",", status, ") ~", ifelse(is.null(group), "1",
                                                                               paste(group, collapse = "+")))), data = data, conf.type = conf.type,
                            conf.int = CIwidth)
  nFit <- sum(sfit$n, na.rm = TRUE)
  if (nrow(data) - nFit == 1) {
    message("1 observation with missing data was removed.")
  }
  else if (nrow(data) > nFit)
    message(paste(nrow(data) - nFit, "observations with missing data were removed."))
  if (!is.null(survtimes)) {
    ofit <- summary(sfit, times = survtimes, extend = !is.null(survtimes))
    colsToExtract <- which(names(sfit) %in% c("strata",
                                              "time", "sr", "surv", "lower", "upper"))
    tb <- data.frame(do.call(cbind, ofit[colsToExtract]))
    if (eventProb == F){
      tb$sr <- apply(tb[, c("surv", "lower", "upper")], 1,
                     psthr, digits)
    }
    if (eventProb == T){
      tb$surv <- 1-tb$surv
      tb$lower_temp <- tb$lower; tb$upper_temp <- tb$upper
      tb$upper <- 1-tb$lower_temp
      tb$lower <- 1-tb$upper_temp
      tb$sr <- apply(tb[, c("surv", "lower", "upper")], 1,
                     psthr, digits)
    }
    if (!is.null(group)) {
      tb$strata <- factor(tb$strata, levels = unique(as.numeric(ofit$strata)),
                          labels = levels(ofit$strata))
      w <- matrix(nrow = length(unique(tb$strata)), ncol = length(unique(tb$time)),
                  dimnames = list(unique(tb$strata), unique(tb$time)))
      for (i in 1:nrow(tb)) w[which(rownames(w) == tb$strata[i]),
                              which(colnames(w) == tb$time[i])] <- tb$sr[i]
    }
    else {
      w <- matrix(nrow = 1, ncol = length(unique(tb$time)))
      colnames(w) <- unique(tb$time)
      for (i in 1:nrow(tb)) w[1, which(colnames(w) ==
                                         tb$time[i])] <- tb$sr[i]
    }
  }
  else w <- NULL
  mtbl <- summary(sfit)$table
  if (inherits(mtbl, "matrix")) {
    m_CI <- apply(mtbl[, grep("median|LCL|UCL", colnames(mtbl))],
                  1, function(x) psthr(x, y = digits))
    m_CI <- gsub("NA","Not Estimable",m_CI)
    lr <- survival::survdiff(as.formula(paste0("survival::Surv(",
                                               time, ",", status, ") ~", paste0(group, collapse = "+"))),
                             data = data)
    nt <- paste0(lr$obs, "/", lr$n)
    w <- cbind(nt, m_CI, w)
    df <- length(lr$obs) - 1
    gl <- rownames(w)
    for (v in group) gl <- gsub(paste0(v, "="), "", gl)
    tab <- cbind(gl, data.frame(w))
    rownames(tab) <- NULL
    nm <- c("Group", "Events/Total",
            paste0("Median (", CIwidth * 100, "%CI)"))
    if (!is.null(survtimes))
      nm <- c(nm, paste0(survtimesLbls, survtimeunit,
                         " (", CIwidth * 100, "% CI)"))
    names(tab) <- nm
    if (showLogrank == TRUE){
      tab <- rbind(tab, c(rep("", length(survtimes)), "Log Rank Test",
                          "ChiSq", paste(niceNum(lr$chisq, 1), "on", df, "df")))
      tab <- rbind(tab, c(rep("", length(survtimes) + 1),
                          "p-value", formatp(pchisq(lr$chisq, df, lower.tail = FALSE))))
    }
  }
  else {
    m_CI <- psthr(mtbl[grep("median|LCL|UCL", names(mtbl))],
                  y = digits)
    m_CI <- gsub("NA","Not Estimable",m_CI)
    nt <- paste0(mtbl["events"], "/", mtbl["n.start"])
    w <- cbind(nt, m_CI, w)
    tab <- data.frame(w)
    rownames(tab) <- NULL
    nm <- c("Events/Total", paste0("Median (", CIwidth *
                                     100, "%CI)"))
    if (!is.null(survtimes))
      nm <- c(nm, paste0(survtimesLbls, survtimeunit,
                         " (", CIwidth * 100, "% CI)"))
    names(tab) <- nm
  }
  if (tableOnly) {
    return(tab)
  }
  outTable(tab, caption = caption, align = paste0("l", paste0(rep("r",
                                                                  ncol(tab) - 1), collapse = "")))
}






#' Summarize cumulative incidence by group
#'
#' Displays event counts and event rates at specified time points for the entire
#' cohort and by group. Gray's test of differences in cumulative incidence is
#' displayed.
#'
#' @param data data frame containing survival data
#' @param time string indicating survival time variable
#' @param status string indicating event status variable; must have at least 3
#'   levels, e.g. 0 = censor, 1 = event, 2 = competing risk
#' @param group string or character vector indicating the variable to group
#'   observations by
#' @param eventcode numerical variable indicating event, default is 1
#' @param cencode numerical variable indicating censored observation, default is
#'   0
#' @param eventtimes numeric vector specifying when event probabilities should
#'   be calculated
#' @param eventtimeunit unit of time to suffix to the time column label if event
#'   probabilities are requested, should be plural
#' @param eventtimeLbls if supplied, a vector the same length as eventtimes with
#'   descriptions (useful for displaying years with data provided in months)
#' @param CIwidth width of the event probabilities, default is 95%
#' @param unformattedp boolean indicating if you would like the p-value to be
#'   returned unformatted (ie not rounded or prefixed with '<'). Should be used
#'   in conjunction with the digits argument.
#' @param na.action default is to omit missing values, but can be set to throw
#'   and error using na.action='na.fail'
#' @param showCounts boolean indicating if the at risk, events and censored
#'   columns should be output, default is TRUE
#' @param showGraystest boolean indicating Gray's test should be included in the
#'   final table, default is TRUE
#' @param digits the number of digits to report in the event probabilities,
#'   default is 2.
#' @param caption table caption for markdown output
#' @param tableOnly should a dataframe or a formatted object be returned
#'
#' @return A character vector of the event table source code, unless
#'   tableOnly=TRUE in which case a data frame is returned
#' @export
#'
#' @examples
#' library(survival)
#' data(pbc)
#'
#' # Event probabilities at various time points with replacement time labels
#' rm_cifsum(data=pbc,time='time',status='status',
#' eventtimes=c(1825,3650),eventtimeLbls=c(5,10),eventtimeunit='yr')
#'
#' # Event probabilities by one group
#' rm_cifsum(data=pbc,time='time',status='status',group='trt',
#' eventtimes=c(1825,3650),eventtimeunit='day')
#'
#'
#' # Event probabilities by multiple groups
#' rm_cifsum(data=pbc,time='time',status='status',group=c('trt','sex'),
#' eventtimes=c(1825,3650),eventtimeunit='day')
#'
rm_cifsum <- function (data, time, status, group = NULL, eventcode = 1, cencode = 0, eventtimes,
                       eventtimeunit, eventtimeLbls = NULL, CIwidth = 0.95, unformattedp = FALSE,
                       na.action = "na.omit", showCounts = TRUE, showGraystest = TRUE,
                       digits = 2, caption = NULL, tableOnly = FALSE
) {

  if (missing(data))
    stop("data is a required argument")
  if (missing(time))
    stop("time is a required argument")
  if (missing(eventtimes))
    stop("eventtimes is a required argument")
  if (missing(eventtimeunit))
    stop("eventtimeunit is a required argument. Example eventtimeunit=\"year\"")
  if (!is.null(eventtimes))
    timelbl <- paste0("Time (", eventtimeunit, ")")
  if (!inherits(data, "data.frame"))
    stop("data must be supplied as a data frame")
  if (!inherits(time, "character") | length(time) > 1)
    stop("time must be supplied as a string indicating a variable in data")
  if (!inherits(status, "character") | length(status) > 1)
    stop("status must be supplied as a string indicating a variable in data")
  if (!is.null(eventtimes))
    if (is.null(eventtimeLbls))
      eventtimeLbls <- eventtimes
  if (length(eventtimeLbls) != length(eventtimes))
    stop("If supplied, the eventtimeLbls vector must be the same length as eventtimes")
  if (length(unique(data[[status]]))==2)
    stop("Only two unique statuses exist in the data. Consider using rm_survsum when competing risks are absent.")
  if (unformattedp)
    formatp <- function(x, ...) {
      x
    }
  missing_vars = na.omit(setdiff(c(time, status, group), names(data)))
  if (length(missing_vars) > 0)
    stop(paste("These variables are not in the data:\n",
               paste0(missing_vars, collapse = csep())))

  data <- data.frame(data)
  n.full <- nrow(data)
  data <- na.omit(data[,c(time, status, group)])
  n.non.missing <- nrow(data)
  if (n.full - n.non.missing == 1) {
    message("1 observation with missing data was removed.")
  }
  else if (n.full > n.non.missing)
    message(paste(n.full - n.non.missing, "observations with missing data were removed."))


  if (!is.null(group)){
    if (length(group)==1){
      data[[group]] <- as.factor(data[[group]])
      fit <- cmprsk::cuminc(ftime = data[[time]], fstatus = data[[status]], group = data[[group]], cencode = 0)
    }
    if (length(group)>1){
      data[["group.combined"]] <- as.factor(apply(data[,group], 1, paste0, collapse=", "))
      fit <- cmprsk::cuminc(ftime = data[[time]], fstatus = data[[status]], group = data[["group.combined"]], cencode = 0)
      group <- "group.combined"
    }
  }
  else{
    fit <- cmprsk::cuminc(ftime = data[[time]], fstatus = data[[status]], cencode = 0)
  }


  event.comb <- data.frame()
  for (i in 1:length(eventtimes)){
    dat <- data.frame(cmprsk::timepoints(fit, eventtimes[i]))
    dat2 <- dat[which(as.numeric(substr(rownames(dat), nchar(rownames(dat))-1, nchar(rownames(dat)))) == eventcode),]
    names(dat2) <- c("cif", "var")
    dat2$lower <- dat2$cif ** exp((-1*abs(qnorm((1-CIwidth)/2))*sqrt(dat2$var))/(dat2$cif*log(dat2$cif)))
    dat2$upper <- dat2$cif ** exp((abs(qnorm((1-CIwidth)/2))*sqrt(dat2$var))/(dat2$cif*log(dat2$cif)))
    if (!is.null(group)){
      dat2$strata <- levels(data[[group]])
    }
    dat2$time <- rep(eventtimes[i], nrow(dat2))

    event.comb <- rbind(event.comb, dat2)
  }

  event.comb$sr <- apply(event.comb[, c("cif", "lower", "upper")], 1, psthr, digits)

  if (!is.null(group)) {
    w <- matrix(nrow = length(unique(event.comb$strata)), ncol = length(unique(event.comb$time)),
                dimnames = list(unique(event.comb$strata), unique(event.comb$time)))
    for (i in 1:nrow(event.comb)) w[which(rownames(w) == event.comb$strata[i]),
                                    which(colnames(w) == event.comb$time[i])] <- event.comb$sr[i]

    if (showCounts == T){
      num.event <- vector()
      num.total <- vector()
      for(i in 1:length(levels(data[[group]]))){
        datai<- data[data[[group]]==levels(data[[group]])[i] & !is.na(data[[group]]),]

        num.event.i <- nrow(datai[!is.na(datai[[time]]) & datai[[status]]==eventcode,])
        num.total.i <- nrow(datai[!is.na(datai[[time]] & !is.na(datai[[status]])),])

        num.event <- c(num.event, num.event.i)
        num.total <- c(num.total, num.total.i)
      }
      event.total <- paste0(num.event, "/", num.total)
      tab <- data.frame(cbind(levels(data[[group]]), event.total, w))
      names(tab) <- c("Strata", "Event/Total", paste0(eventtimeLbls, eventtimeunit, " (", CIwidth * 100, "% CI)"))
    }
    else{
      tab <- data.frame(cbind(levels(data[[group]]), w))
      names(tab) <- c("Strata", paste0(eventtimeLbls, eventtimeunit, " (", CIwidth * 100, "% CI)"))
    }
  }
  else {
    if (showCounts == T){
      w <- matrix(nrow = 1, ncol = length(unique(event.comb$time)))
      colnames(w) <- unique(event.comb$time)
      for (i in 1:nrow(event.comb)) w[1, which(colnames(w) == event.comb$time[i])] <- event.comb$sr[i]

      datai <- data[which(!is.na(data[[time]]) & !is.na(data[[status]])),]
      num.event <- nrow(datai[!is.na(datai[[time]]) & datai[[status]]==eventcode,])
      num.total <- nrow(datai)
      event.total <- paste0(num.event, "/", num.total)
      tab <- data.frame(cbind(event.total, w))
      names(tab) <- c("Event/Total", paste0(eventtimeLbls, eventtimeunit, " (", CIwidth * 100, "% CI)"))
    }
    else {
      tab <- data.frame(w)
      names(tab) <- c(paste0(eventtimeLbls, eventtimeunit, " (", CIwidth * 100, "% CI)"))
    }
  }

  if (!is.null(group) & showGraystest == TRUE){
    gray <- fit$Tests[which(as.numeric(rownames(fit$Tests))==eventcode),]
    tab <- rbind(tab, c(rep("", ncol(tab)-3), "Gray's Test",
                        "ChiSq", paste(niceNum(gray[1], 1), "on", round(gray[3]), "df")))
    tab <- rbind(tab, c(rep("", ncol(tab)-2),
                        "p-value", formatp(gray[2])))
  }
  if (tableOnly) {
    return(tab)
  }
  outTable(tab, caption = caption, align = paste0("l", paste0(rep("r",
                                                                  ncol(tab) - 1), collapse = "")))


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
#' @param conf.type type of confidence interval see \code{\link[survival:survfit]{survival::survfit}} for
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
#' @seealso \code{\link[survival:survfit]{survival::survfit}}
#' @return A character vector of the survival table source code, unless tableOnly=TRUE in
#'   which case a data frame is returned
#' @export
#' @examples
#' # Kaplan-Mieir survival probabilities with time displayed in years
#' data("pembrolizumab")
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
                        na.action='na.omit',showCounts=TRUE,digits=getOption("reportRmd.digits",2),caption=NULL,tableOnly=FALSE,fontsize){
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

