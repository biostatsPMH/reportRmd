
# Testing for getVarLevels and coeffSum ------------

data("pembrolizumab")
devtools::load_all()
## Linear -----
uv_lm <- lm(age~sex,data=pembrolizumab)
coeffSum(uv_lm)
getVarLevels(uv_lm)
m_summary(uv_lm)

auto_lm <- autoreg.rm_lm("age",data=pembrolizumab,"sex")
coeffSum(auto_lm)
getVarLevels(auto_lm)
m_summary(auto_lm)

mv_lm <- lm(pdl1 ~ age+sex+cohort,data = pembrolizumab)
gp(mv_lm)
coeffSum(mv_lm)
getVarLevels(mv_lm)
m_summary(mv_lm)
mv_lm2 <- lm(pdl1 ~ age*sex+sex*cohort+age*l_size,data = pembrolizumab)
gp(mv_lm2)
coeffSum(mv_lm2)
getVarLevels(mv_lm2)
m_summary(mv_lm2)

uvsum_lm <- rm_uvsum("age", covs = c("sex", "l_size"), data = pembrolizumab)
uvsum2_lm <- rm_uvsum2("age", covs = c("sex", "l_size"), data = pembrolizumab)
## Binomial -----
uv_binom <- glm(formula=as.formula(orr~sex),data=pembrolizumab,family="binomial")
coeffSum(uv_binom)
getVarLevels(uv_binom)
m_summary(uv_binom)

auto_binom <- autoreg.rm_glm("orr",data=pembrolizumab,"sex",family="binomial")
coeffSum(auto_binom)
getVarLevels(auto_binom)
m_summary(auto_binom)

mv_binom <- glm(orr~age+sex+cohort,family = 'binomial',data = pembrolizumab)
gp(mv_binom)
coeffSum(mv_binom)
getVarLevels(mv_binom)
m_summary(mv_binom)

mv_binom2 <- glm(orr~age:sex+cohort,family = 'binomial',data = pembrolizumab)
gp(mv_binom2)
coeffSum(mv_binom2)
getVarLevels(mv_binom2)
m_summary(mv_binom2)

uvsum_binom <- rm_uvsum("orr", covs = c("age", "sex", "cohort"), data = pembrolizumab, family = "binomial")
uvsum2_binom <- rm_uvsum2("orr", covs = c("age", "sex", "cohort"), data = pembrolizumab, family = "binomial")

## Poisson -----
pembrolizumab$int_var <- rpois(n=nrow(pembrolizumab),lambda = 1)
uv_pois <- glm(formula=as.formula(int_var~age),data=pembrolizumab, family="poisson",offset=NULL)
coeffSum(uv_pois)
getVarLevels(uv_pois)
m_summary(uv_pois)

auto_pois <- autoreg.rm_glm("int_var",data=pembrolizumab,"age",family="poisson",offset=NULL)
coeffSum(auto_pois)
getVarLevels(auto_pois)
m_summary(auto_pois)

mv_pois <- glm(formula = int_var ~ age+sex+cohort, family = poisson, data = pembrolizumab)
gp(mv_pois)
coeffSum(mv_pois)
getVarLevels(mv_pois)
m_summary(mv_pois)

mv_pois2 <-glm(formula = int_var ~ age+sex:cohort, family = poisson, data = pembrolizumab)
gp(mv_pois2)
coeffSum(mv_pois2)

getVarLevels(mv_pois2)
m_summary(mv_pois2)

uvsum_pois <- rm_uvsum("int_var", covs = c("age", "sex", "cohort"), data = pembrolizumab, family = "poisson",offset=NULL)
uvsum2_pois <- rm_uvsum2("int_var", covs = c("age", "sex", "cohort"), data = pembrolizumab, family = "poisson",offset=NULL)
## Negative Binomial -----
uv_negbin <- MASS::glm.nb(formula=as.formula(int_var~sex),data=pembrolizumab)
coeffSum(uv_negbin)
getVarLevels(uv_negbin)
m_summary(uv_negbin)

auto_negbin <- autoreg.rm_negbin("int_var",data=pembrolizumab,"sex",offset=NULL)
coeffSum(auto_negbin)
getVarLevels(auto_negbin)
m_summary(auto_negbin)

mv_negbin <- MASS::glm.nb(int_var~age+sex+cohort,data=pembrolizumab,link=log)
gp(mv_negbin)
coeffSum(mv_negbin)
getVarLevels(mv_negbin)
m_summary(mv_negbin)

mv_negbin2 <- MASS::glm.nb(int_var~age*sex+cohort,data=pembrolizumab,link=log)
gp(mv_negbin2)
coeffSum(mv_negbin2)
getVarLevels(mv_negbin2)
m_summary(mv_negbin2)

uvsum_negbin <- rm_uvsum("int_var", covs = c("age", "sex", "cohort"), data = pembrolizumab, offset=NULL, type = "negbin")
uvsum2_negbin <- rm_uvsum2("int_var", covs = c("age", "sex", "cohort"), data = pembrolizumab, offset=NULL, type = "negbin")

## Ordinal -----
pembrolizumab$ord_var <- factor(ifelse(pembrolizumab$int_var>2,2,pembrolizumab$int_var),ordered = TRUE)

uv_ord <- MASS::polr(formula=as.formula(ord_var~sex),data=pembrolizumab)
coeffSum(uv_ord)
getVarLevels(uv_ord)
m_summary(uv_ord)

auto_ord <- autoreg.rm_ordinal("ord_var",data=pembrolizumab, "sex")
coeffSum(auto_ord)
getVarLevels(auto_ord)
m_summary(auto_ord)

mv_ord <- MASS::polr(ord_var ~ age+ sex+cohort, data = pembrolizumab, Hess = TRUE,
                     method = "logistic")
gp(mv_ord)
coeffSum(mv_ord)
getVarLevels(mv_ord)
m_summary(mv_ord)


mv_ord2 <- MASS::polr(ord_var ~ age:sex+cohort, data = pembrolizumab, Hess = TRUE,
                      method = "logistic")
gp(mv_ord2)
coeffSum(mv_ord2)
getVarLevels(mv_ord2)
m_summary(mv_ord2)

uvsum_ord <- rm_uvsum("ord_var", covs = c("age", "sex", "cohort"), data = pembrolizumab)
uvsum2_ord <- rm_uvsum2("ord_var", covs = c("age", "sex", "cohort"), data = pembrolizumab)


## Cox PH -----
uv_cox <- coxph(Surv(os_time, os_status) ~ sex, data = pembrolizumab)
coeffSum(uv_cox)
getVarLevels(uv_cox)
m_summary(uv_cox)

auto_cox <- autoreg.rm_coxph(response=c("os_time","os_status"),data=pembrolizumab,"sex",family=NULL,offset=NULL,id=NULL,strata = "")
coeffSum(auto_cox)
getVarLevels(auto_cox)
m_summary(auto_cox)

mv_cox <- survival::coxph(Surv(os_time,os_status) ~ age+sex+cohort, data = pembrolizumab)
gp(mv_cox)
coeffSum(mv_cox)
getVarLevels(mv_cox)
m_summary(mv_cox)

mv_cox2 <- survival::coxph(Surv(os_time,os_status) ~ age:sex+cohort, data = pembrolizumab)
model <- survival::coxph(Surv(os_time,os_status) ~ age:sex+cohort, data = pembrolizumab)
gp(mv_cox2)
coeffSum(mv_cox2)
getVarLevels(mv_cox2)
m_summary(mv_cox2)


#
uvsum_cox <- rm_uvsum(response=c("os_time","os_status"), covs = c("age", "sex", "cohort"), data = pembrolizumab,family=NULL,offset=NULL,id=NULL,strata = "")
uvsum2_cox <- rm_uvsum2(response=c("os_time","os_status"), covs = c("age", "sex", "cohort"), data = pembrolizumab,family=NULL,offset=NULL,id=NULL,strata = "")
uvsum_cox <- rm_uvsum(response=c("os_time","os_status"), covs = c("age", "sex", "cohort"), data = pembrolizumab,family=NULL,offset=NULL,id=NULL)
uvsum2_cox <- rm_uvsum2(response=c("os_time","os_status"), covs = c("age", "sex", "cohort"), data = pembrolizumab,family=NULL,offset=NULL,id=NULL)

## CRR -----
pembrolizumab$os_status2 <- pembrolizumab$os_status
pembrolizumab$os_status2[sample(1:nrow(pembrolizumab),10,replace = FALSE)] <-2

uv_crr <- crrRx(as.formula('os_time+os_status2~ cohort'),data=pembrolizumab)
coeffSum(uv_crr)
getVarLevels(uv_crr)
m_summary(uv_crr)

auto_crr <- autoreg.rm_crr(response=c("os_time","os_status2"), data=pembrolizumab, "sex")
coeffSum(auto_crr)
getVarLevels(auto_crr)
m_summary(auto_crr)

mv_crr <-crrRx(as.formula('os_time+os_status2~age+sex+cohort'),data=pembrolizumab)
gp(mv_crr)
coeffSum(mv_crr)
getVarLevels(mv_crr)
m_summary(mv_crr)

mv_crr2 <-crrRx(as.formula('os_time+os_status2~age*sex+age:cohort'),data=pembrolizumab)
gp(mv_crr2)
coeffSum(mv_crr2)
getVarLevels(mv_crr2)
m_summary(mv_crr2)

uvsum_crr <- rm_uvsum(response=c("os_time","os_status2"), covs = c("age", "sex", "cohort"), data = pembrolizumab)
uvsum2_crr <- rm_uvsum2(response=c("os_time","os_status2"), covs = c("age", "sex", "cohort"), data = pembrolizumab)

## GEE  (NEED MORE MODELS HERE) -----
pembrolizumab$orr2 <- ifelse(pembrolizumab$orr=="CR/PR",1,0)
uv_gee <- geepack::geeglm(formula=as.formula(orr2~sex),data=pembrolizumab,id=pembrolizumab$id)
coeffSum(uv_gee)
getVarLevels(uv_gee)
m_summary(uv_gee)

auto_gee <- autoreg.rm_gee("orr",data=pembrolizumab, "sex",id="id")
coeffSum(auto_gee)
getVarLevels(auto_gee)


mv_gee <- geepack::geeglm(orr2 ~ age+sex+cohort, family = binomial,
                          data = pembrolizumab,
                          id = pembrolizumab$id, corstr = "independence")
mv_gee_glm <- glm(orr2 ~ age+sex+cohort, family = binomial,
                          data = pembrolizumab)
rm_mvsum2(mv_gee,vif=T)
rm_mvsum2(mv_gee_glm,vif=T)

mv_gee <- geepack::geeglm(orr2 ~ age+sex+cohort, family = binomial,
                          data = pembrolizumab,
                          id = 1:nrow(pembrolizumab), corstr = "independence")
rm_mvsum2(mv_gee,vif=T)

gp(mv_gee)
coeffSum(mv_gee)
getVarLevels(mv_gee)
mv_gee2 <- geepack::geeglm(orr2 ~ age:sex+cohort, family = binomial,
                          data = pembrolizumab,
                          id = pembrolizumab$id, corstr = "independence")
gp(mv_gee2)
coeffSum(mv_gee2)
getVarLevels(mv_gee2)

uvsum_gee <- rm_uvsum(response="orr", covs = c("age", "sex", "cohort"), data = pembrolizumab)
uvsum2_gee <- rm_uvsum2(response="orr", covs = c("age", "sex", "cohort"), data = pembrolizumab)


uvsum2(response = c('os_time','os_status'),
                  covs=c('age','sex','baseline_ctdna','l_size','change_ctdna_group'),
                  data=pembrolizumab,CIwidth=.9)


uvsum2(response = c('os_time','os_status2'),
       covs=c('age','sex','baseline_ctdna','l_size','change_ctdna_group'),
       data=pembrolizumab,CIwidth=.9)

rm_uvsum2(response=c("os_time","os_status"), covs = c("age", "sex", "cohort"), data = pembrolizumab,family=NULL,offset=NULL,id="",strata = "")


# works fine
data("pembrolizumab")
rm_compactsum(data = pembrolizumab, xvars = c("age",
                                              "change_ctdna_group", "l_size", "pdl1","cohort"), grp = "sex", use_mean = "age",
              digits = c("age" = 2, "l_size" = 3), digits.cat = 1)

# but not if everyone is in the same groups
pembrolizumab$change_ctdna_group2 = "changed"
rm_compactsum(data = pembrolizumab, xvars = c("age",
                                              "change_ctdna_group2", "l_size", "pdl1"), grp = "sex", use_mean = "age",
              digits = c("age" = 2, "l_size" = 3), digits.cat = 1)

# This looks good
rm_compactsum(data = pembrolizumab, xvars = c("age",
                                              "change_ctdna_group", "l_size", "pdl1","cohort"), grp = "sex", use_mean = "age",
              digits = c("age" = 2, "l_size" = 1), digits.cat = 1,all.stats = F)

#  So do all these
rm_compactsum(data = pembrolizumab, xvars = c("age",
                                              "change_ctdna_group", "l_size", "pdl1","cohort"), grp = "sex", use_mean = "age",
              digits = c("age" = 2, "l_size" = 1), digits.cat = 1,all.stats = T)

rm_compactsum(data = pembrolizumab, xvars = c("age"), grp = "sex", use_mean = "age",all.stats = T)


rm_compactsum(data = pembrolizumab, xvars = c("age"), grp = "sex", all.stats = T)

rm_compactsum(data = pembrolizumab, xvars = c("age"), grp = "sex", all.stats = T)

rm_compactsum(data = pembrolizumab, xvars = c("age","l_size"), grp = "sex", all.stats = T)

rm_compactsum(data = pembrolizumab, xvars = c("cohort"), grp = "sex", all.stats = T)

rm_compactsum(data = pembrolizumab, xvars = c("age","cohort"), grp = "sex", all.stats = T)

# This gives the same results for me for the p-values ?
uvsum_ord <- rm_uvsum("ord_var", covs = c("age", "sex", "cohort"), data = pembrolizumab)
uvsum2_ord <- rm_uvsum2("ord_var", covs = c("age", "sex", "cohort"), data = pembrolizumab)

mv_binom <- glm(orr~age+sex+cohort,family = 'binomial',data = pembrolizumab)
# This is performing regular CI (shouldn't!!)
rm_mvsum(mv_binom)
exp(confint.default(mv_binom))
# This is performing Profile Likelihood CI - which is what we want
rm_mvsum2(mv_binom)
exp(confint(mv_binom))

mv_binom2 <- glm(orr~age:sex+cohort,family = 'binomial',data = pembrolizumab)
# This is performing regular CI (shouldn't!!)
rm_mvsum(mv_binom2)
exp(confint.default(mv_binom2))
# This is performing Profile Likelihood CI - which is what we want
rm_mvsum2(mv_binom2)
exp(confint(mv_binom2))

mv_lm2 <- lm(pdl1 ~ age*sex+sex*cohort+age*l_size,data = pembrolizumab)
rm_mvsum(mv_lm2)


rm_mvsum(model)
rm_mvsum2(model)


# GEE Models
# Note - these don't work with factors for responses like regular glm, response variable has to be numeric
mv_gee <- geepack::geeglm(os_status ~ age+sex+cohort,
                          family = binomial,
                          data = pembrolizumab,
                          id = pembrolizumab$id, corstr = "independence")

rm_mvsum(mv_gee)
rm_mvsum2(mv_gee)
model <- mv_gee


## LM GEE fitting
## To do: add models for binomial & poisson and test against rm_mvsum
mv_lm_gee <- geepack::geeglm(pdl1 ~ age+sex+cohort,data = pembrolizumab,
                             id=pembrolizumab$id, # always need to specify for gee models
                             family = gaussian)
rm_mvsum2(mv_lm_gee)

# geeglm will fail if the interactions are too complex
mv_lm2_gee <- geepack::geeglm(pdl1 ~ age*sex+cohort+l_size,data = pembrolizumab,
                              id=pembrolizumab$id,
                              family = gaussian)
rm_mvsum2(mv_lm2_gee)

data("pembrolizumab")

tab <- pembrolizumab |>
  rm_compactsum(xvars = c(pdl1,age,change_ctdna_group),grp="sex")

attributes(tab)
rm_mvsum(mv_lm2_gee, showN = T, showEvent = T, vif = T, whichp = "both")
rm_mvsum2(mv_lm2_gee, showN = T, showEvent = T, vif = T, whichp = "both")


mv_binom2_gee <- geepack::geeglm(orr2~age:sex+cohort,family = 'binomial',data = pembrolizumab,id=1:nrow(pembrolizumab))

rm_mvsum(mv_binom2_gee, showN = T, showEvent = T, vif = T, whichp = "both")
rm_mvsum2(mv_binom2_gee, showN = T, showEvent = T, vif = T, whichp = "both")


lung[,"sex",drop=FALSE]
pembrolizumab[,"sex"]


# THIS GIVES AN ERROR
pembrolizumab$Counts <- rpois(nrow(pembrolizumab),lambda = 3)
pembrolizumab$length_followup <- rnorm(nrow(pembrolizumab),mean = 72,sd=3)

rm_uvsum2(data=pembrolizumab, response='Counts', type='negbin',
         covs=c('age','cohort'),offset = "log(length_followup)")



model <- MASS::glm.nb(Counts ~ age + offset(log(length_followup)), link = log,
             data = pembrolizumab)
rm_mvsum(model)
rm_mvsum2(model)
names(model$model)


# Sep 2024 Errors from Katherine
data("pembrolizumab")
rm_compactsum(
  data = pembrolizumab,
  xvars = "age")
rm_compactsum(
  data = pembrolizumab,
  xvars = "age",
  all.stats = T
)
