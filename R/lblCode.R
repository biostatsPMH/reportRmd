
#' Clear variable labels from a data frame
#'
#' This function will remove all label attributes from variables in the data.
#'
#' To change or remove individual labels use set_labels or set_var_labels
#' @param data the data frame to remove labels from
#' @export
#' @examples
#' # Set a few variable labels for ctDNA
#' ctDNA <- ctDNA |> set_var_labels(
#'    ctdna_status="detectable ctDNA",
#'   cohort="A cohort label")
#' # Clear all variable data frames and check
#' clear_labels(ctDNA)
clear_labels <- function(data){
  for (v in names(data)) {
    v_att <- attributes(data[[v]])
    v_lbl_att <- setdiff(unique(c("label",grep('label',names(v_att),value=TRUE))),"labels")
    for (a in v_lbl_att) attr(data[[v]],a) <- NULL
  }
  return(data)
}


#' Extract variable labels from labelled data frame
#'
#' Extract variable labels from data and return a data frame with labels
#'
#' All variable names will be returned, even those with no labels.
#' If the label attribute has length greater than one the values will be
#' concatenated and returned as a single string separated by sep
#' @param data the data frame to extract labels from
#' @param sep character used to separate multiple labels, defaults to "_"
#' @export
#' @examples
#' # Set a few variable labels for ctDNA
#' ctDNA <- ctDNA |> set_var_labels(
#'    ctdna_status="detectable ctDNA",
#'   cohort="A cohort label")
#' # Extract labels
#' extract_labels(ctDNA)
extract_labels <- function(data,sep="_"){
  lbls <-sapply(names(data), function(v){
    lbl = paste(attr(data[[v]],"label"),sep=sep)
    ifelse(length(lbl)==0,NA,lbl)
  })
  rtn <- data.frame(variable=names(data),label=lbls)
  rownames(rtn) <- NULL
  return(rtn)
}

#' Set variable labels
#'
#' Set variable labels for a data frame using name-label pairs.
#'
#' If no label is provided for a variable then the existing label will not be
#' changed. To remove a label set the label to NA.
#' @param data data frame containing variables to be labelled
#' @param ... Name-label pairs the name gives the name of the column in the
#'   output and the label is a character vector of length one.
#' @export
#' @examples
#' # set labels using name-label pairs
#' # and return labelled data frame
#' ctDNA |> set_var_labels(
#'    ctdna_status="detectable ctDNA",
#'   cohort="A cohort label")
set_var_labels = function (data, ...) {
  args <- as.list(match.call())[-1]
  dnm <- as.character(args[1])
  varLbls <- args[-1]
  for (i in seq_along(varLbls)){
    v <- names(varLbls)[i]
    l <- varLbls[i]
    if (!v %in% names(data)) {
      message(paste(v,'not a variable in',dnm,"\nLabel not added."))
    } else {
      if (length(l)>1) message(paste("Label for",v,'has more than one element.\n Only the first will be used.'))
      attr(data[[v]], "label") <- as.character(l[1])
    }
  }
  return(data)
}


#' Set variable labels
#'
#' Assign variable labels to a data.frame from a lookup table.
#'
#' Useful if variable labels have been imported from a data dictionary. The
#' first column in names_labels must contain the variable name and the second
#' column the variable label. The column names are not used.
#'
#' If no label is provided then the existing label will not be changed. To
#' remove a label set the label to NA.
#'
#' @param data data frame to be labelled
#' @param names_labels data frame with column 1 containing variable names from
#'   data and column 2 containing variable labels. Other columns will be
#'   ignored.
#' @export
#' @examples
#' # create data frame with labels
#' lbls <- data.frame(c1=c('cohort','size_change'),
#' c2=c('Cancer cohort','Change in tumour size'))
#' # set labels and return labelled data frame
#' set_labels(ctDNA,lbls)
set_labels <- function(data,names_labels){
  if (!inherits(data,"data.frame")) {
    stop("data must be a data frame")
  }
  if (!inherits(names_labels,"data.frame")) {
    stop("names_labels must be a data frame")
  }
  if (ncol(names_labels)<2) stop("names_labels must be a data frame with at least two columns")

  if (ncol(names_labels)>2) message("The names_labels data frame contains more than two columns.\nVariable names will be taken from the first column and variable labels from the second column.")
  for (v in 1:ncol(names_labels)) names_labels[[v]] <- as.character(names_labels[[v]] )

  varIndx <- which (names_labels[[1]] %in% names(data))
  v_lbls <- names_labels[varIndx,]

  for (i in 1:nrow(v_lbls)) attr(data[[v_lbls[[1]][i]]], "label") <- v_lbls[[2]][i]
  return(data)
}

# return variable labels associated with variables
replaceLbl <- function(data_arg,cv){
  if (!inherits(data_arg,"character")) dn <- paste(deparse1(data_arg),collapse="") else dn <- data_arg
  if (!inherits(dn,'character')) stop('data table must be specified as a character string.')
  if (!inherits(cv,'character')) stop('variable name must be specified as a character string.')
  lbl <- extract_labels(get0(dn))
  vl <- data.frame(variable=cv,ord=1:length(cv))
  if (!is.null(lbl)){
    cvnew <- merge(vl,lbl,all.x=T)
    cvnew <- cvnew[order(cvnew$ord),]
    cvnew <- cvnew[!duplicated(cvnew),]
  } else {
    cvnew <- vl
    cvnew$label <- NA
  }
  cvnew$label <- ifelse(is.na(cvnew$label),nicename(cvnew$variable),cvnew$label)
  return(cvnew$label)
}


