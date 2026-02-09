#'  Survival data 
#'  
#'  Survival status and ctDNA levels for patients receiving pembrolizumab 
#'  
#' @usage data('pembrolizumab')
#' @format A data frame with 94 rows and 15 variables:
#' \describe{ 
#'   \item{id}{Patient ID}
#'   \item{age}{Age at study entry}
#'   \item{sex}{Patient Sex}
#'   \item{cohort}{Study Cohort}
#'   \item{l_size}{Target lesion size at baseline}
#'   \item{pdl1}{PD L1 percent}
#'   \item{tmb}{log of TMB}
#'   \item{baseline_ctdna}{Baseline ctDNA}
#'   \item{change_ctdna_group}{Did ctDNA increase or decrease from baseline to cycle 3}
#'   \item{orr}{Objective Response}
#'   \item{cbr}{Clinical Beneficial Response}
#'   \item{os_status}{Overall survival status}
#'   \item{os_time}{Overall survival time in months}
#'   \item{pfs_status}{Progression free survival status}
#'   \item{pfs_time}{Progression free survival time in months}
#' } 
#' @source \url{https://www.nature.com/articles/s43018-020-0096-5} 
"pembrolizumab" 


#'  Tumour size change over time 
#'  
#'  Longitudinal changes in tumour size since baseline for patients by changes in ctDNA status (clearance, decrease or increase) since baseline. 
#'  
#' @usage data('ctDNA')
#' @format A data frame with 270 rows and 5 variables:
#' \describe{ 
#'   \item{id}{Patient ID}
#'   \item{cohort}{Study Cohort}
#'   \item{ctdna_status}{Change in ctDNA since baseline}
#'   \item{time}{Number of weeks on treatment}
#'   \item{size_change}{Percentage change in tumour measurement}
#' } 
#' @source \url{https://www.nature.com/articles/s43018-020-0096-5} 
"ctDNA" 


# NOTE: uvmodels data object is deprecated and no longer used
# Model type mapping is now handled by get_model_class() in model_registry.R
# This provides better maintainability and eliminates the need for a separate data file
# The uvmodels.rda file can be removed from data/ directory in a future release

