# The data for the package comes from the Bratman paper in Nature Cancer
# Bratman, S.V., Yang, S.Y.C., Iafolla, M.A.J. et al. 
# Personalized circulating tumor DNA analysis as a predictive biomarker in 
# solid tumor patients treated with pembrolizumab. 
# Nat Cancer 1, 873–881 (2020). https://doi.org/10.1038/s43018-020-0096-5
# Manuscript & Data links: https://www.nature.com/articles/s43018-020-0096-5
library(tidyverse)
data_path <- 'data-raw/PMH_data.xlsx'

# Figure 2 data becomes the ctDNA dataset (some variables)

# # This was run initally and the pembrolizumab sheet was placed into the Excel
# # Figure 3 data + Supplemental Table 3 data 
# # is merged to become the pembrolizumab data
# fig3 <- readxl::read_excel('data-raw/PMH_data.xlsx',sheet = 'Fig3')
# fig3 <- fig3 %>%
#   select(!c(baseline_ctDNA,C3_ctDNA,change_ctDNA))
# st3 <- readxl::read_excel('data-raw/PMH_data.xlsx',sheet = 'ST3')
# pembrolizumab <- full_join(fig3,st3)
# pembrolizumab %>% write_csv('pembrolizumab.csv')

sheets <- readxl::excel_sheets(data_path)
url <- 'https://www.nature.com/articles/s43018-020-0096-5'

# Save the data information to R/data.R
# Import File Descriptions and corresponding data file

for (f in grep('_var',sheets,value=T)){
  var_info <- readxl::read_excel('data-raw/PMH_data.xlsx',sheet = f)
  datafile <- readxl::read_excel('data-raw/PMH_data.xlsx',sheet = gsub('_var.*','',f))
  
  # Extract the dataset name, description and details
  name <- var_info$Description[which(var_info$Variable=='NAME')]
  if (length(name)==0) {
    sink()
    stop('Each variable information sheet needs to specify the Name, Description and Details of the Dataset.')
  }

  details <- var_info$Description[var_info$Variable=='DETAILS']
  description <- var_info$Description[var_info$Variable=='DESCRIPTION']
  var_info <- var_info %>%
    filter(!Variable %in% c('NAME','DESCRIPTION','DETAILS'))
  
    datafile <- datafile %>%
    select(all_of(var_info$Variable))
  
  # Ensure all variable names are in lowercase 
  # & that spaces are converted to '_'
  names(datafile) <- tolower(names(datafile))
  names(datafile) <- gsub(' ','_',names(datafile))
  var_info$Variable <- tolower(var_info$Variable)
  var_info$Variable <- gsub(' ','_',var_info$Variable)
  
  
  # Write the Documentation File
  sink('R/data.R',append = !(f==grep('_var',sheets,value=T)[1]))
  cat("#' ",name,'\n')
  cat("#' ",'\n')
  cat("#' ",description,'\n')
  cat("#' ",'\n')
  cat("#' ",details,'\n')
  cat("#' ",'\n')
  
  cat("#' @format A data frame with ",nrow(datafile)," rows and ",ncol(datafile)," variables:\n",sep="")
  cat("#' \\describe{ \n",sep="")
  for (i in 1:nrow(var_info)){
    cat("#'   \\item{",var_info$Variable[i],"}{",var_info$Description[i],"}\n",sep="")
  }
  cat("#' } \n")
  cat("#' @source \\url{",url,"} \n",sep="")
  cat('"',name,'"',' \n',sep="")
  cat("\n\n")
  # Close the data documentation file
  sink()
  
  # Change all character variables to factors
  for (v in setdiff(names(datafile),'id')){
      if ('character' %in% class(datafile[[v]]))  datafile[[v]] <- factor(datafile[[v]])
    }
  # Output the Data
  eval(parse(text = paste(name,' <- datafile')))
  eval(parse(text = paste0("usethis::use_data(",name,", overwrite = TRUE)")))
}


