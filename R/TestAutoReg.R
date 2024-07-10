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

# Binomial
response="orr"
class(response) <-c(class(response),"rm_glm")
binom_fit <- autoreg(response,data=pembrolizumab,x_var,family="binomial")


# Poisson
pembrolizumab$int_var <- rpois(n=nrow(pembrolizumab),lambda = 1)
response="int_var"
class(response) <-c(class(response),"rm_glm")
pois_fit <- autoreg(response,data=pembrolizumab,x_var,family="poisson",offset=NULL)


# Negative Binomial
response="int_var"
class(response) <-c(class(response),"rm_negbin")
neg_fit <- autoreg(response,data=pembrolizumab,x_var,family="poisson",offset=NULL)
neg_fit
class(neg_fit)

# Ordinal
pembrolizumab$ord_var <- factor(ifelse(pembrolizumab$int_var>2,2,pembrolizumab$int_var),ordered = T)
response="ord_var"
class(response) <-c(class(response),"rm_ordinal")
ord_fit <- autoreg(response,data=pembrolizumab,x_var,family=NULL,offset=NULL)

# Cox PH -
response = c("os_time","os_status")
class(response) <-c(class(response),"rm_coxph")
cox_fit <- autoreg(response,data=pembrolizumab,x_var="sex",family=NULL,offset=NULL,id=NULL,strata = "")

# GEE Model - binomial

response="orr"
class(response) <-c(class(response),"rm_gee")

gee_out <- autoreg(response,data=pembrolizumab,x_var,family="binomial",id=id)
gee_out

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



# TO DO:
# - test all the different autoreg functions
# - test all the coeffSum functions
# write a function to accept data, response, covs and use the derive_type function to assign a class to the response
