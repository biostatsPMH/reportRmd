pkgname <- "reportRmd"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('reportRmd')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("clear_labels")
### * clear_labels

flush(stderr()); flush(stdout())

### Name: clear_labels
### Title: Clear variable labels from a data frame
### Aliases: clear_labels

### ** Examples

# Set a few variable labels for ctDNA
data("ctDNA")
ctDNA <- ctDNA |> set_var_labels(
   ctdna_status="detectable ctDNA",
  cohort="A cohort label")
# Clear all variable data frames and check
clear_labels(ctDNA)



cleanEx()
nameEx("crrRx")
### * crrRx

flush(stderr()); flush(stdout())

### Name: crrRx
### Title: fit crr model
### Aliases: crrRx
### Keywords: model

### ** Examples

# From the crr help file:
set.seed(10)
ftime <- rexp(200)
fstatus <- sample(0:2,200,replace=TRUE)
cov <- matrix(runif(600),nrow=200)
dimnames(cov)[[2]] <- c('x1','x2','x3')
df <- data.frame(ftime,fstatus,cov)
m1 <- crrRx(as.formula('ftime+fstatus~x1+x2+x3'),df)
# Nicely output to report:
rm_mvsum(m1,data=df,showN = TRUE,vif=TRUE)



cleanEx()
nameEx("excelCol")
### * excelCol

flush(stderr()); flush(stdout())

### Name: excelCol
### Title: Retrieve columns number from spreadsheet columns specified as
###   unquoted letters
### Aliases: excelCol

### ** Examples

## Find the column numbers for excel columns AB, CE and BB
excelCol(AB,CE,bb)
## Get the columns between A and K and Z
excelCol(A-K,Z)



cleanEx()
nameEx("excelColLetters")
### * excelColLetters

flush(stderr()); flush(stdout())

### Name: excelColLetters
### Title: Retrieve spreadsheet column letter-names from columns indices
### Aliases: excelColLetters

### ** Examples

## Find the column numbers for excel columns AB, CE and BB
colIndices <- excelCol(AB,CE,bb)
## Go back to the column names
excelColLetters(colIndices)



cleanEx()
nameEx("extract_labels")
### * extract_labels

flush(stderr()); flush(stdout())

### Name: extract_labels
### Title: Extract variable labels from labelled data frame
### Aliases: extract_labels

### ** Examples

# Set a few variable labels for ctDNA
data("ctDNA")
ctDNA <- ctDNA |> set_var_labels(
   ctdna_status="detectable ctDNA",
  cohort="A cohort label")
# Extract labels
extract_labels(ctDNA)



cleanEx()
nameEx("extract_package_details")
### * extract_package_details

flush(stderr()); flush(stdout())

### Name: extract_package_details
### Title: Extract Function and Package Information from Current Document
### Aliases: extract_package_details

### ** Examples

## Not run: 
##D # Run this function from within an R script to analyze its dependencies
##D package_info <- extract_package_details()
##D print(package_info)
## End(Not run)




cleanEx()
nameEx("fillNAs")
### * fillNAs

flush(stderr()); flush(stdout())

### Name: fillNAs
### Title: Forward fill NA values (vectorized implementation)
### Aliases: fillNAs
### Keywords: internal

### ** Examples

## Not run: 
##D fillNAs(c(1, NA, NA, 2, NA, 3))  # Returns: c(1, 1, 1, 2, 2, 3)
## End(Not run)



cleanEx()
nameEx("forestplotMV")
### * forestplotMV

flush(stderr()); flush(stdout())

### Name: forestplotMV
### Title: Create a multivariable forest plot using ggplot2
### Aliases: forestplotMV
### Keywords: plot

### ** Examples

data("pembrolizumab")
glm_fit <- glm(orr ~ change_ctdna_group + sex + age + l_size,
               data = pembrolizumab, family = 'binomial')

# Adjusted only
forestplotMV(glm_fit, data = pembrolizumab)

# Both adjusted and unadjusted
forestplotMV(glm_fit, data = pembrolizumab, include_unadjusted = TRUE)



cleanEx()
nameEx("forestplotUV")
### * forestplotUV

flush(stderr()); flush(stdout())

### Name: forestplotUV
### Title: Create an univariable forest plot using ggplot2
### Aliases: forestplotUV
### Keywords: plot

### ** Examples

data("pembrolizumab")
forestplotUV(response = "orr",
             covs = c("change_ctdna_group", "sex", "age", "l_size"),
             data = pembrolizumab, family = 'binomial')



cleanEx()
nameEx("formatp")
### * formatp

flush(stderr()); flush(stdout())

### Name: formatp
### Title: Format p-values for tables (HTML/LaTeX)
### Aliases: formatp
### Keywords: helper

### ** Examples

## Not run: 
##D formatp(c(0.0001, 0.045, 0.123))
##D # Returns: c("<0.001", "0.045", "0.12")
## End(Not run)



cleanEx()
nameEx("ggkmcif")
### * ggkmcif

flush(stderr()); flush(stdout())

### Name: ggkmcif
### Title: Plot KM and CIF curves with ggplot
### Aliases: ggkmcif

### ** Examples

data("pembrolizumab")
# Simple plot without confidence intervals
ggkmcif(response = c('os_time','os_status'),
cov='cohort',
data=pembrolizumab)

# Plot with median survival time
ggkmcif(response = c('os_time','os_status'),
cov='sex',
data=pembrolizumab,
median.text = TRUE,median.lines=TRUE,conf.curves=TRUE)

# Plot with specified survival times and log-log CI
ggkmcif(response = c('os_time','os_status'),
cov='sex',
data=pembrolizumab,
median.text = FALSE,set.time.text = 'mo OS',
set.time = c(12,24),conf.type = 'log-log',conf.curves=TRUE)

# KM plot with 95% CI and censor marks
ggkmcif(c('os_time','os_status'),'sex',data = pembrolizumab, type = 'KM',
HR=TRUE, HR_pval = TRUE, conf.curves = TRUE,conf.type='log-log',
set.time.CI = TRUE, censor.marks=TRUE)



cleanEx()
nameEx("ggkmcif2_2025")
### * ggkmcif2_2025

flush(stderr()); flush(stdout())

### Name: ggkmcif2_2025
### Title: Plot KM and CIF curves with ggplot
### Aliases: ggkmcif2_2025

### ** Examples

# Simple plot without confidence intervals
data("pembrolizumab")
ggkmcif2(response = c('os_time','os_status'),
cov='cohort',
data=pembrolizumab)

# Plot with median survival time
ggkmcif2(response = c('os_time','os_status'),
cov='sex',
data=pembrolizumab,
median.text = TRUE,median.lines=TRUE,conf.curves=TRUE)

# Plot with specified survival times and log-log CI
ggkmcif2(response = c('os_time','os_status'),
cov='sex',
data=pembrolizumab,
median.text = FALSE,set.time.text = 'mo OS',
set.time = c(12,24),conf.type = 'log-log',conf.curves=TRUE)

# KM plot with 95% CI and censor marks
ggkmcif2(c('os_time','os_status'),'sex',data = pembrolizumab, type = 'KM',
HR=TRUE, HR_pval = TRUE, conf.curves = TRUE,conf.type='log-log',
set.time.CI = TRUE, censor.marks=TRUE)



cleanEx()
nameEx("ggkmcif_paste")
### * ggkmcif_paste

flush(stderr()); flush(stdout())

### Name: ggkmcif_paste
### Title: Plot KM and CIF curves with ggplot
### Aliases: ggkmcif_paste

### ** Examples

data("pembrolizumab")
plot <- ggkmcif(response=c('pfs_time','pfs_status'),
data=pembrolizumab,returns = TRUE)

# Highlighting a section of the curve
plot[[1]] <- plot[[1]] +
ggplot2::geom_rect(xmin=4,xmax=8,ymin=0.15,ymax=0.4,alpha=0.01,fill='yellow')

# Putting the curve back together
ggkmcif_paste(plot)



cleanEx()
nameEx("lpvalue2")
### * lpvalue2

flush(stderr()); flush(stdout())

### Name: lpvalue2
### Title: Format p-values for plot annotations
### Aliases: lpvalue2
### Keywords: helper

### ** Examples

## Not run: 
##D lpvalue2(0.0001, 3)  # Returns: "p < 0.001"
##D lpvalue2(0.0456, 3)  # Returns: "p = 0.046"
## End(Not run)



cleanEx()
nameEx("m_summary")
### * m_summary

flush(stderr()); flush(stdout())

### Name: m_summary
### Title: Output a table for multivariate or univariate regression models
### Aliases: m_summary
### Keywords: internal

### ** Examples

## Not run: 
##D data("pembrolizumab")
##D uv_lm <- lm(age~sex,data=pembrolizumab)
##D m_summary(uv_lm, digits = 3, for_plot = FALSE)
##D 
##D mv_binom <- glm(orr~age+sex+cohort,family = 'binomial',data = pembrolizumab)
##D m_summary(mv_binom, whichp = "both", for_plot = TRUE)
## End(Not run)



cleanEx()
nameEx("nestTable")
### * nestTable

flush(stderr()); flush(stdout())

### Name: nestTable
### Title: Combine two table columns into a single column with levels of
###   one nested within levels of the other.
### Aliases: nestTable

### ** Examples

## Investigate models to predict baseline ctDNA and tumour size and display together
## (not clinically useful!)
data(pembrolizumab)
fit1 <- lm(baseline_ctdna~age+l_size+pdl1,data=pembrolizumab)
m1 <- rm_mvsum(fit1,tableOnly=TRUE)
m1$Response = 'ctDNA'
fit2 <- lm(l_size~age+baseline_ctdna+pdl1,data=pembrolizumab)
m2 <- rm_mvsum(fit2,tableOnly=TRUE)
m2$Response = 'Tumour Size'
nestTable(rbind(m1,m2),head_col='Response',to_col='Covariate')



cleanEx()
nameEx("niceNum")
### * niceNum

flush(stderr()); flush(stdout())

### Name: niceNum
### Title: Round numbers with trailing zeros
### Aliases: niceNum
### Keywords: helper

### ** Examples

## Not run: 
##D niceNum(c(1.5, 2.345, NA), digits = 2)
##D # Returns: c("1.50", "2.35", NA)
## End(Not run)



cleanEx()
nameEx("plotuv")
### * plotuv

flush(stderr()); flush(stdout())

### Name: plotuv
### Title: Plot multiple bivariate relationships in a single plot
### Aliases: plotuv
### Keywords: plot

### ** Examples

## Run multiple univariate analyses on the pembrolizumab dataset to predict cbr and
## then visualise the relationships.
data("pembrolizumab")
rm_uvsum(data=pembrolizumab,
response='cbr',covs=c('age','sex','l_size','baseline_ctdna'))
plotuv(data=pembrolizumab,  response='cbr',
covs=c('age','sex','l_size','baseline_ctdna'),showN=TRUE)



cleanEx()
nameEx("psthr")
### * psthr

flush(stderr()); flush(stdout())

### Name: psthr
### Title: Round and paste with parentheses (smart formatting)
### Aliases: psthr
### Keywords: helper

### ** Examples

## Not run: 
##D psthr(c(1.234, 1.123, 1.345), 2)
##D # Returns: "1.23 (1.12, 1.35)"
##D psthr(c(0.001, 0.0005, 0.0015), 2)
##D # Returns: "1.0e-03 (5.0e-04, 1.5e-03)"
## End(Not run)



cleanEx()
nameEx("pstprn")
### * pstprn

flush(stderr()); flush(stdout())

### Name: pstprn
### Title: Paste vector elements with parentheses
### Aliases: pstprn
### Keywords: helper

### ** Examples

## Not run: 
##D pstprn(c(1.5, 1.2, 1.8))
##D # Returns: "1.5 (1.2, 1.8)"
## End(Not run)



cleanEx()
nameEx("replace_plot_labels")
### * replace_plot_labels

flush(stderr()); flush(stdout())

### Name: replace_plot_labels
### Title: Replace variable names with labels in ggplot
### Aliases: replace_plot_labels

### ** Examples

## Not run: 
##D data("pembrolizumab")
##D p <- ggplot(pembrolizumab,aes(x=change_ctdna_group,y=baseline_ctdna)) +
##D geom_boxplot()
##D replace_plot_labels(p)
##D pembrolizumab <- set_var_labels(pembrolizumab,
##D change_ctdna_group="Change in ctDNA group")
##D p <- ggplot(pembrolizumab,aes(x=change_ctdna_group,y=baseline_ctdna)) +
##D geom_boxplot()
##D replace_plot_labels(p)
##D # Can also be used with a pipe, but expression needs to be wrapped in a brace
##D (ggplot(pembrolizumab,aes(x=change_ctdna_group,y=baseline_ctdna)) +
##D geom_boxplot()) |> replace_plot_labels()
## End(Not run)



cleanEx()
nameEx("rm_cifsum")
### * rm_cifsum

flush(stderr()); flush(stdout())

### Name: rm_cifsum
### Title: Summarize cumulative incidence by group
### Aliases: rm_cifsum

### ** Examples

library(survival)
data(pbc)

# Event probabilities at various time points with replacement time labels
rm_cifsum(data=pbc,time='time',status='status',
eventtimes=c(1825,3650),eventtimeLbls=c(5,10),eventtimeunit='yr')

# Event probabilities by one group
rm_cifsum(data=pbc,time='time',status='status',group='trt',
eventtimes=c(1825,3650),eventtimeunit='day')


# Event probabilities by multiple groups
rm_cifsum(data=pbc,time='time',status='status',group=c('trt','sex'),
eventtimes=c(1825,3650),eventtimeunit='day')




cleanEx()
nameEx("rm_compactsum")
### * rm_compactsum

flush(stderr()); flush(stdout())

### Name: rm_compactsum
### Title: Output a compact summary table
### Aliases: rm_compactsum

### ** Examples

data("pembrolizumab")
rm_compactsum(data = pembrolizumab, xvars = c("age",
"change_ctdna_group", "l_size", "pdl1"), grp = "sex", use_mean = "age",
digits = c("age" = 2, "l_size" = 3), digits.cat = 1, iqr = TRUE,
show.tests = TRUE)

# Other Examples (not run)
## Include the summary statistic in the variable column
#rm_compactsum(data = pembrolizumab, xvars = c("age",
#"change_ctdna_group"), grp = "sex", use_mean = "age", show.sumstats=TRUE)

## To show effect sizes
#rm_compactsum(data = pembrolizumab, xvars = c("age",
#"change_ctdna_group"), grp = "sex", use_mean = "age", digits = 2,
#effSize = TRUE, show.tests = TRUE)

## To return unformatted p-values
#rm_compactsum(data = pembrolizumab, xvars = c("l_size",
#"change_ctdna_group"), grp = "cohort", effSize = TRUE, unformattedp = TRUE)

## Using tidyselect
#pembrolizumab |> rm_compactsum(xvars = c(age, sex, pdl1), grp = cohort,
#effSize = TRUE)




cleanEx()
nameEx("rm_covsum")
### * rm_covsum

flush(stderr()); flush(stdout())

### Name: rm_covsum
### Title: Outputs a descriptive covariate table
### Aliases: rm_covsum
### Keywords: dataframe

### ** Examples

data("pembrolizumab")
rm_covsum(data=pembrolizumab, maincov = 'orr',
covs=c('age','sex','pdl1','tmb','l_size','change_ctdna_group'),
show.tests=TRUE)

# To Show Effect Sizes
rm_covsum(data=pembrolizumab, maincov = 'orr',
covs=c('age','sex'),
effSize=TRUE)

# To make custom changes or change the fontsize in PDF/HTML
tab <- rm_covsum(data=pembrolizumab,maincov = 'change_ctdna_group',
covs=c('age','sex','pdl1','tmb','l_size'),show.tests=TRUE,tableOnly = TRUE)
outTable(tab, fontsize=7)

# To return unformatted p-values
tab <- rm_covsum(data=pembrolizumab, maincov = 'orr',
covs=c('age','sex','pdl1','tmb','l_size','change_ctdna_group'),
show.tests=TRUE,unformattedp=TRUE,tableOnly=TRUE)
outTable(tab,digits=5)
outTable(tab,digits=5, applyAttributes=FALSE) # remove bold/indent



cleanEx()
nameEx("rm_mvsum")
### * rm_mvsum

flush(stderr()); flush(stdout())

### Name: rm_mvsum
### Title: Format a regression model nicely for 'Rmarkdown'
### Aliases: rm_mvsum

### ** Examples

data("pembrolizumab")
glm_fit = glm(change_ctdna_group~sex:age+baseline_ctdna+l_size,
data=pembrolizumab,family = 'binomial')
rm_mvsum(glm_fit)

#linear model with p-value adjustment
lm_fit=lm(baseline_ctdna~age+sex+l_size+tmb,data=pembrolizumab)
rm_mvsum(lm_fit,p.adjust = "bonferroni")
#Coxph
require(survival)
res.cox <- coxph(Surv(os_time, os_status) ~ sex+age+l_size+tmb, data = pembrolizumab)
rm_mvsum(res.cox, vif=TRUE)



cleanEx()
nameEx("rm_survdiff")
### * rm_survdiff

flush(stderr()); flush(stdout())

### Name: rm_survdiff
### Title: Display event counts, expected event counts and logrank test of
###   differences
### Aliases: rm_survdiff

### ** Examples

#' # Differences between sex
data("pembrolizumab")
rm_survdiff(data=pembrolizumab,time='os_time',status='os_status',
covs='sex',digits=1)

# Differences between sex, stratified by cohort
rm_survdiff(data=pembrolizumab,time='os_time',status='os_status',
covs='sex',strata='cohort',digits=1)
# Differences between sex/cohort groups
rm_survdiff(data=pembrolizumab,time='os_time',status='os_status',
covs=c('sex','cohort'),digits=1)



cleanEx()
nameEx("rm_survsum")
### * rm_survsum

flush(stderr()); flush(stdout())

### Name: rm_survsum
### Title: Summarise survival data by group
### Aliases: rm_survsum

### ** Examples

# Simple median survival table
data("pembrolizumab")
rm_survsum(data=pembrolizumab,time='os_time',status='os_status')

# Survival table with yearly survival rates
rm_survsum(data=pembrolizumab,time='os_time',status='os_status',
survtimes=c(12,24),survtimesLbls=1:2, survtimeunit='yr')

#Median survival by group
rm_survsum(data=pembrolizumab,time='os_time',status='os_status',group='sex')

# Survival Summary by cohort, displayed in years
rm_survsum(data=pembrolizumab,time='os_time',status='os_status',
group="cohort",survtimes=seq(12,72,12),
survtimesLbls=seq(1,6,1),
survtimeunit='years')

# Survival Summary by Sex and ctDNA group
rm_survsum(data=pembrolizumab,time='os_time',status='os_status',
group=c('sex','change_ctdna_group'),survtimes=c(12,24),survtimeunit='mo')




cleanEx()
nameEx("rm_survtime")
### * rm_survtime

flush(stderr()); flush(stdout())

### Name: rm_survtime
### Title: Display survival rates and events for specified times
### Aliases: rm_survtime

### ** Examples

# Kaplan-Mieir survival probabilities with time displayed in years
data("pembrolizumab")
rm_survtime(data=pembrolizumab,time='os_time',status='os_status',
strata="cohort",type='KM',survtimes=seq(12,72,12),
survtimesLbls=seq(1,6,1),
survtimeunit='years')

# Cox Proportional Hazards survivial probabilities
rm_survtime(data=pembrolizumab,time='os_time',status='os_status',
strata="cohort",type='PH',survtimes=seq(12,72,12),survtimeunit='months')

# Cox Proportional Hazards survivial probabilities controlling for age
rm_survtime(data=pembrolizumab,time='os_time',status='os_status',
covs='age',strata="cohort",survtimes=seq(12,72,12),survtimeunit='months')




cleanEx()
nameEx("rm_uv_mv")
### * rm_uv_mv

flush(stderr()); flush(stdout())

### Name: rm_uv_mv
### Title: Combine univariate and multivariable regression tables
### Aliases: rm_uv_mv

### ** Examples

require(survival)
data("pembrolizumab")
uvTab <- rm_uvsum(response = c('os_time','os_status'),
covs=c('age','sex','baseline_ctdna','l_size','change_ctdna_group'),
data=pembrolizumab,tableOnly=TRUE)
mv_surv_fit <- coxph(Surv(os_time,os_status)~age+sex+
baseline_ctdna+l_size+change_ctdna_group, data=pembrolizumab)
uvTab <- rm_mvsum(mv_surv_fit)

#linear model
uvtab<-rm_uvsum(response = 'baseline_ctdna',
covs=c('age','sex','l_size','pdl1','tmb'),
data=pembrolizumab,tableOnly=TRUE)
lm_fit=lm(baseline_ctdna~age+sex+l_size+tmb,data=pembrolizumab)
mvtab<-rm_mvsum(lm_fit,tableOnly = TRUE)
rm_uv_mv(uvtab,mvtab,tableOnly=TRUE)

#logistic model
uvtab<-rm_uvsum(response = 'os_status',
covs=c('age','sex','l_size','pdl1','tmb'),
data=pembrolizumab,family = binomial,tableOnly=TRUE)
logis_fit<-glm(os_status~age+sex+l_size+pdl1+tmb,data = pembrolizumab,family = 'binomial')
mvtab<-rm_mvsum(logis_fit,tableOnly = TRUE)
rm_uv_mv(uvtab,mvtab,tableOnly=TRUE)



cleanEx()
nameEx("rm_uvsum")
### * rm_uvsum

flush(stderr()); flush(stdout())

### Name: rm_uvsum
### Title: Output several univariate models nicely in a single table
### Aliases: rm_uvsum

### ** Examples

# Examples are for demonstration and are not meaningful
# Coxph model with 90% CI
data("pembrolizumab")
rm_uvsum(response = c('os_time','os_status'),
covs=c('age','sex','baseline_ctdna','l_size','change_ctdna_group'),
data=pembrolizumab,CIwidth=.9)

# Linear model with default 95% CI
rm_uvsum(response = 'baseline_ctdna',
covs=c('age','sex','l_size','pdl1','tmb'),
data=pembrolizumab)

# Logistic model with default 95% CI
rm_uvsum(response = 'os_status',
covs=c('age','sex','l_size','pdl1','tmb'),
data=pembrolizumab,family = binomial)
# Poisson models returned as model list
mList <- rm_uvsum(response = 'baseline_ctdna',
covs=c('age','sex','l_size','pdl1','tmb'),
data=pembrolizumab, returnModels=TRUE)
#'
# GEE on correlated outcomes
data("ctDNA")
rm_uvsum(response = 'size_change',
covs=c('time','ctdna_status'),
gee=TRUE,
id='id', corstr="exchangeable",
family=gaussian("identity"),
data=ctDNA,showN=TRUE)

# Using tidyselect
pembrolizumab |> rm_uvsum(response = sex,
covs = c(age, cohort))



cleanEx()
nameEx("scrolling_table")
### * scrolling_table

flush(stderr()); flush(stdout())

### Name: scrolling_table
### Title: Output a scrollable table
### Aliases: scrolling_table

### ** Examples

data("pembrolizumab")
tab <- rm_covsum(data=pembrolizumab,maincov = 'change_ctdna_group',
covs=c('age','cohort','sex','pdl1','tmb','l_size'),full=F)
scrolling_table(tab,pixelHeight=300)



cleanEx()
nameEx("set_labels")
### * set_labels

flush(stderr()); flush(stdout())

### Name: set_labels
### Title: Set variable labels
### Aliases: set_labels

### ** Examples

data("ctDNA")
# create data frame with labels
lbls <- data.frame(c1=c('cohort','size_change'),
c2=c('Cancer cohort','Change in tumour size'))
# set labels and return labelled data frame
set_labels(ctDNA,lbls)



cleanEx()
nameEx("set_var_labels")
### * set_var_labels

flush(stderr()); flush(stdout())

### Name: set_var_labels
### Title: Set variable labels
### Aliases: set_var_labels

### ** Examples

# set labels using name-label pairs
# and return labelled data frame
data("ctDNA")
ctDNA |> set_var_labels(
   ctdna_status="detectable ctDNA",
  cohort="A cohort label")



cleanEx()
nameEx("xcn")
### * xcn

flush(stderr()); flush(stdout())

### Name: xcn
### Title: Convert Excel column letters to numbers
### Aliases: xcn
### Keywords: internal

### ** Examples

## Not run: 
##D xcn("A")   # Returns 1
##D xcn("Z")   # Returns 26
##D xcn("AA")  # Returns 27
## End(Not run)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
