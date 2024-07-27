# Test Code
library(devtools)
# load_all()
data("pembrolizumab")
head(pembrolizumab)
library(dplyr)

CIwidth=.95;digits=2

# Linear
# Fitting the Models -----------------
response = "age"
class(response) <-c(class(response),"rm_lm")
x_var <- "sex"
id=NULL;strata="";family=NULL;offset=NULL
lm_fit <- autoreg(response,data=pembrolizumab,x_var)
coeffSum(lm_fit)
mv_lm <- lm(pdl1 ~ age+sex+cohort,data = pembrolizumab)
mv_lm_int <- lm(pdl1 ~ age*sex+sex*cohort+age*l_size,data = pembrolizumab)
model <- mv_lm_int
coeffSum(mv_lm_int) # Extra print statements

# Binomial_fit -------------------
response="orr"
class(response) <-c(class(response),"rm_glm")
binom_fit <- autoreg(response,data=pembrolizumab,x_var,family="binomial")
model <- binom_fit
coeffSum(binom_fit)
mv_binom <- glm(orr~age+sex+cohort,family = 'binomial',data = pembrolizumab)
rm_mvsum(mv_binom)
coeffSum(mv_binom) # Not working - but Poisson is

# Poisson -------------------
pembrolizumab$int_var <- rpois(n=nrow(pembrolizumab),lambda = 1)
response="int_var"
class(response) <-c(class(response),"rm_glm")
pois_fit <- autoreg(response,data=pembrolizumab,x_var,family="poisson",offset=NULL)
coeffSum(pois_fit)
mv_pois_fit <-  glm(formula = int_var ~ age+sex+cohort, family = poisson, data = pembrolizumab)
coeffSum(mv_pois_fit)

rm_mvsum(mv_pois_fit)

# Negative Binomial -------------------
response="int_var"
class(response) <-c(class(response),"rm_negbin")
neg_fit <- autoreg(response,data=pembrolizumab,x_var,family="poisson",offset=NULL)
coeffSum(neg_fit)
class(neg_fit)
mv_nb <- MASS::glm.nb(int_var~age+sex+cohort,data=pembrolizumab,link=log)
rm_mvsum(mv_nb)
coeffSum(mv_nb)

# Ordinal -------------------
pembrolizumab$ord_var <- factor(ifelse(pembrolizumab$int_var>2,2,pembrolizumab$int_var),ordered = T)
response="ord_var"
class(response) <-c(class(response),"rm_ordinal")
ord_fit <- autoreg(response,data=pembrolizumab,x_var,family=NULL,offset=NULL)
coeffSum(ord_fit)
mv_ord <- MASS::polr(ord_var ~ age+ sex+cohort, data = pembrolizumab, Hess = TRUE,
                     method = "logistic")
coeffSum(mv_ord) # Error
rm_mvsum(mv_ord)

# Cox PH -------------------
response = c("os_time","os_status")
class(response) <-c(class(response),"rm_coxph")
cox_fit <- autoreg(response,data=pembrolizumab,x_var="sex",family=NULL,offset=NULL,id=NULL,strata = "")
model <- cox_fit
coeffSum(cox_fit,whichp="both") # global-p doesn't do anything yet
library(survival)
mv_cox <- survival::coxph(Surv(os_time,os_status) ~ age+sex+cohort, data = pembrolizumab)
rm_mvsum(mv_cox)
coeffSum(mv_cox) # Doesn't work

# CRR -------------------
pembrolizumab$os_status2 <- pembrolizumab$os_status
pembrolizumab$os_status2[sample(1:nrow(pembrolizumab),10,replace = F)] <-2
response = c("os_time","os_status2")
class(response) <-c(class(response),"rm_crr")
crr_fit <- autoreg(response,data=pembrolizumab,x_var="sex",family=NULL,offset=NULL,id=NULL,strata = "")
coeffSum(crr_fit,whichp = "both") # will return NA for only 1 variable, that's fine

mv_crr <- crrRx(as.formula('os_time+os_status2~age+sex+cohort'),data=pembrolizumab)
coeffSum(mv_crr) # funny ordering
rm_mvsum(mv_crr,data = pembrolizumab) # This errors, but example below does not

# From the crr help file - this works
set.seed(10)
ftime <- rexp(200)
fstatus <- sample(0:2,200,replace=TRUE)
cov <- matrix(runif(600),nrow=200)
dimnames(cov)[[2]] <- c('x1','x2','x3')
cov2 <- data.frame(cov) |>
  dplyr::mutate(x4=cut(x3,seq(0,1,.25)))
df <- data.frame(ftime,fstatus,cov)
m1 <- crrRx(as.formula('ftime+fstatus~x1+x2+x3'),df)
rm_mvsum(m1,data=df,showN = T)
coeffSum(m1) # this works

df2 <- data.frame(ftime,fstatus,cov2)
m2 <- crrRx(as.formula('ftime+fstatus~x1+x2+x4'),df2)
rm_mvsum(m2,data=df2,showN = T)
coeffSum(m2)


# GEE Model - binomial  - both should work now
response="orr"; x_var <- "sex"
class(response) <-c(class(response),"rm_gee")

gee_out <- autoreg(response,data=pembrolizumab,x_var,family="binomial",id="id")

pembrolizumab$orr2 <- ifelse(pembrolizumab$orr=="CR/PR",1,0)
response="orr2"; x_var <- "sex"
class(response) <-c(class(response),"rm_gee")

gee_out <- autoreg(response,data=pembrolizumab,x_var,family="binomial",id="id")



# Summarising the model output -----------------
lm_sum <- coeffSum(lm_fit,CIwidth=0.95)
class(lm_sum)
typeof(lm_sum)
coeffSum(binom_fit, CIwidth=0.95)
coeffSum(pois_fit, CIwidth=0.95)
coeffSum(neg_fit, CIwidth=0.95)
coeffSum(ord_fit, CIwidth=0.95)
coeffSum(cox_fit, CIwidth=0.95)
coeffSum(gee_out,CIwidth=0.95)
coeffSum(crr_fit,CIwidth=0.95)

# for testing global p
model <- survival::coxph(formula = as.formula(survival::Surv(os_time,
                                                             os_status) ~ age+sex+cohort), data = pembrolizumab)

model <- MASS::polr(formula = ord_var ~ age+sex+cohort, data = pembrolizumab, Hess = TRUE,
           method = "logistic")

idf <- as.numeric(factor(pembrolizumab$id))
model <- geepack::geeglm(orr2 ~ age+sex+cohort, family = binomial,
                         data = pembrolizumab,
                         id = idf, corstr = "independence")


crr_mod <- crrRx(f = (os_time + os_status ~ sex + age + cohort), data = pembrolizumab)

pembrolizumab$os_status2 <- pembrolizumab$os_status
pembrolizumab$os_status2[sample(1:nrow(pembrolizumab),10,replace = F)] <-2
response = c("os_time","os_status2")
class(response) <-c(class(response),"rm_crr")
crr_uv <- autoreg(response,data=pembrolizumab,x_var="sex",family=NULL,offset=NULL,id=NULL,strata = "")
coeffSum.crr(crr_uv)
coeffSum(crr_uv)

coeffSum.crr(crr_mod)
