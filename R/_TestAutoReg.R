# Test Code
library(devtools)
# load_all()
data("pembrolizumab")
head(pembrolizumab)

# Linear

# Fitting the Models -----------------
response = "age"
class(response) <-c(class(response),"rm_lm")
x_var <- "sex"
id=NULL;strata="";family=NULL;offset=NULL
lm_fit <- autoreg(response,data=pembrolizumab,x_var)
coeffSum(lm_fit)

# Binomialm_fit# Binomial
response="orr"
class(response) <-c(class(response),"rm_glm")
binom_fit <- autoreg(response,data=pembrolizumab,x_var,family="binomial")
model <- binom_fit
coeffSum(binom_fit)

# Poisson
pembrolizumab$int_var <- rpois(n=nrow(pembrolizumab),lambda = 1)
response="int_var"
class(response) <-c(class(response),"rm_glm")
pois_fit <- autoreg(response,data=pembrolizumab,x_var,family="poisson",offset=NULL)
coeffSum(pois_fit)

mv_pois_fit <-  glm(formula = int_var ~ age+cohort, family = poisson, data = pembrolizumab)
coeffSum(mv_pois_fit)

# Negative Binomial
response="int_var"
class(response) <-c(class(response),"rm_negbin")
neg_fit <- autoreg(response,data=pembrolizumab,x_var,family="poisson",offset=NULL)
coeffSum(neg_fit)
class(neg_fit)

# Ordinal
pembrolizumab$ord_var <- factor(ifelse(pembrolizumab$int_var>2,2,pembrolizumab$int_var),ordered = T)
response="ord_var"
class(response) <-c(class(response),"rm_ordinal")
ord_fit <- autoreg(response,data=pembrolizumab,x_var,family=NULL,offset=NULL)
coeffSum(ord_fit)

# Cox PH -
response = c("os_time","os_status")
class(response) <-c(class(response),"rm_coxph")
cox_fit <- autoreg(response,data=pembrolizumab,x_var="sex",family=NULL,offset=NULL,id=NULL,strata = "")
model <- cox_fit
coeffSum(cox_fit,whichp="both") # global-p doesn't do anything yet

# CRR
pembrolizumab$os_status2 <- pembrolizumab$os_status
pembrolizumab$os_status2[sample(1:nrow(pembrolizumab),10,replace = F)] <-2
response = c("os_time","os_status2")
class(response) <-c(class(response),"rm_crr")
crr_fit <- autoreg(response,data=pembrolizumab,x_var="sex",family=NULL,offset=NULL,id=NULL,strata = "")
coeffSum.crr(crr_fit,whichp = "both") # will return NA for only 1 variable, that's fine


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
