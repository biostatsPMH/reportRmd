#'  ctDNA 
#'  
#'  Circulating tumour DNA levels 
#'  
#'  ctDNA levels and RECIST status for patients with solid tumor patients treated with pembrolizumab 
#'  
#' @format A data frame with 186 rows and 9 variables:
#' \describe{ 
#'   \item{id}{Patient ID}
#'   \item{cohort}{Study Cohort: A = Squamous cell carcinoma of soft pallate, B = Triple negative breast cancer, C = Ovarian, high grade serous, D = Melanoma, E = Other Solid Tumor}
#'   \item{log10_change_ctdna}{Transformed value of the change in cDNA over the interval}
#'   \item{interval_start_week}{Start of interval measured since initiation of therapy}
#'   \item{interval_end_week}{End of interval measured since initiation of therapy}
#'   \item{interval_recist}{Recist 1.1 status}
#'   \item{reason_off_trial}{Reason for removal from trial}
#'   \item{os_status}{Overall survival status}
#'   \item{on_treatment}{Was patient on treatment at end of trial}
#' } 
#' @source \url{https://www.nature.com/articles/s43018-020-0096-5} 
"ctDNA" 


#'  pembrolizumab 
#'  
#'  Survival data 
#'  
#'  Survival status and ctDNA levels for patients receiving pembrolizumab 
#'  
#' @format A data frame with 94 rows and 15 variables:
#' \describe{ 
#'   \item{id}{Patient ID}
#'   \item{age}{Age at study entry}
#'   \item{sex}{Patient Sex}
#'   \item{cohort}{Study Cohort: A = Squamous cell carcinoma of soft pallate, B = Triple negative breast cancer, C = Ovarian, high grade serous, D = Melanoma, E = Other Solid Tumor}
#'   \item{l_size}{Target lesion size at baseline}
#'   \item{pdl1}{PD L1 percent}
#'   \item{tmb}{log of TMB}
#'   \item{baseline_ctdna}{Baseline ctDNA}
#'   \item{change_ctdna_group}{Did ctDNA increase or decrease from baseline to cycle 3}
#'   \item{orr}{Objective Response}
#'   \item{cbr}{Clinical Beneficial Response}
#'   \item{os_status}{Overall survival status, 0 = alive, 1 = deceased}
#'   \item{os_time}{Overall survival time in months}
#'   \item{pfs_status}{Progression free survival status, 0 = progression free, 1 = progressed}
#'   \item{pfs_time}{Progression free survival time in months}
#' } 
#' @source \url{https://www.nature.com/articles/s43018-020-0096-5} 
"pembrolizumab" 


