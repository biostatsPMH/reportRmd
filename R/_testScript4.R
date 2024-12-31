
# Testing for getVarLevels and coeffSum ------------
data("pembrolizumab")

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
m_summary(mv_binom,forceWald = F)
m_summary(mv_binom,forceWald = T)
tab <- m_summary(mv_binom)
rm_mvsum(mv_binom,showN=F,showEvent=F)
mv_binom2 <- glm(orr~age:sex+cohort,family = 'binomial',data = pembrolizumab)
gp(mv_binom2)
coeffSum(mv_binom2)
getVarLevels(mv_binom2)
m_summary(mv_binom2)

uvsum_binom <- rm_uvsum("orr", covs = c("age", "sex", "cohort"), data = pembrolizumab, family = "binomial")
rm_uvsum("orr", covs = c("age", "sex", "cohort"), data = pembrolizumab, family = "binomial",forceWald = T)

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
rm_mvsum2(mv_crr)

mv_crr2 <-crrRx(as.formula('os_time+os_status2~age*sex+age:cohort'),data=pembrolizumab)
gp(mv_crr2)
coeffSum(mv_crr2)
getVarLevels(mv_crr2)
m_summary(mv_crr2)
rm_mvsum2(mv_crr2,whichp = "both")

uvsum_crr <- rm_uvsum(response=c("os_time","os_status2"), covs = c("age", "sex", "cohort"), data = pembrolizumab)
uvsum2_crr <- rm_uvsum2(response=c("os_time","os_status2"), covs = c("age", "sex", "cohort"), data = pembrolizumab)

# GEE ---------------
pembrolizumab$orr2 <- ifelse(pembrolizumab$orr=="CR/PR",1,0)

mv_gee2 <- geepack::geeglm(orr2 ~ age:sex+cohort, family = binomial,
                          data = pembrolizumab,
                          id = pembrolizumab$id, corstr = "independence")
gp(mv_gee2)
coeffSum(mv_gee2)
getVarLevels(mv_gee2)

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
                                              "change_ctdna_group", "l_size", "pdl1","cohort"),
              grp = "sex", use_mean = "age",
              digits = c("age" = 2, "l_size" = 3), digits.cat = 1)

rm_compactsum(data = pembrolizumab, xvars = c("age",
                                              "change_ctdna_group", "l_size", "pdl1","cohort"),
              grp = "sex", use_mean = T, percentage = "col",
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
rm_mvsum(mv_binom2,whichp = "both")
exp(confint.default(mv_binom2))
# This is performing Profile Likelihood CI - which is what we want
rm_mvsum2(mv_binom2,whichp = "both")
exp(confint(mv_binom2))

rm_mvsum2(mv_pois,whichp = "both")
exp(confint(mv_pois))

mv_lm2 <- lm(pdl1 ~ age*sex+sex*cohort+age*l_size,data = pembrolizumab)
rm_mvsum(mv_lm,whichp = "both")


rm_mvsum(model)
rm_mvsum2(model,whichp="both")


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
rm_mvsum(mv_lm_gee)
rm_mvsum2(mv_lm_gee)

uv_lm_gee <- rm_uvsum2("pdl1", covs = c("age", "sex", "cohort"),data = pembrolizumab,
                             id=pembrolizumab$id, # always need to specify for gee models
                             family = gaussian, type = "gee")


# geeglm will fail if the interactions are too complex
mv_lm2_gee <- geepack::geeglm(pdl1 ~ age*sex+cohort+l_size,data = pembrolizumab,
                              id=pembrolizumab$id,
                              family = gaussian)
rm_mvsum(mv_lm2_gee)
rm_mvsum2(mv_lm2_gee)


## Poisson
mv_pois_gee <- geepack::geeglm(pdl1 ~ age+sex+cohort,data = pembrolizumab,
                             id=pembrolizumab$id, # always need to specify for gee models
                             family = "poisson")
rm_mvsum(mv_pois_gee)
rm_mvsum2(mv_pois_gee)

uv_pois_gee <- rm_uvsum2("pdl1", covs = c("age", "sex", "cohort"),data = pembrolizumab,
                      id=pembrolizumab$id, # always need to specify for gee models
                      family = poisson, type = "gee")

mv_pois2_gee <- geepack::geeglm(pdl1 ~ age*sex+cohort+l_size,data = pembrolizumab,
                              id=pembrolizumab$id,
                              family = "poisson")
rm_mvsum(mv_pois2_gee)
rm_mvsum2(mv_pois2_gee)

## Binomial
pembrolizumab$binary <- rbinom(nrow(pembrolizumab), 1, 0.5)
mv_binom_gee <- geepack::geeglm(binary ~ age+sex+cohort,data = pembrolizumab,
                               id=pembrolizumab$id, # always need to specify for gee models
                               family = "binomial")
rm_mvsum(mv_binom_gee)
rm_mvsum2(mv_binom_gee)

uv_binom_gee <- rm_uvsum2("pdl1", covs = c("age", "sex", "cohort"),data = pembrolizumab,
                        id=pembrolizumab$id, # always need to specify for gee models
                        family = binomial, type = "gee")

mv_binom2_gee <- geepack::geeglm(binary ~ age*sex+cohort+l_size,data = pembrolizumab,
                                id=pembrolizumab$id,
                                family = "binomial")
# mv_binom2_gee <- geepack::geeglm(orr~age:sex+cohort,data = pembrolizumab,
#                                  id=pembrolizumab$id,
#                                  family = "binomial")
rm_mvsum(mv_binom2_gee)
rm_mvsum2(mv_binom2_gee)


# errors from spreadsheet
# all ok
rm_uvsum(response=c("os_time","os_status"), covs = c("age", "sex", "cohort"), data = pembrolizumab,family=NULL,offset=NULL,id=NULL,strata = "")
rm_uvsum(response=c("os_time","os_status"), covs = c("age", "sex", "cohort"), data = pembrolizumab,family=NULL,offset=NULL,id=NULL,strata = NULL)
rm_uvsum(response=c("os_time","os_status"), covs = c("age", "sex", "cohort"), data = pembrolizumab,family=NULL,offset=NULL,id=NULL,strata = NA)
rm_uvsum2(response=c("os_time","os_status"), covs = c("age", "sex", "cohort"), data = pembrolizumab,family=NULL,offset=NULL,id=NULL,strata = NA)
rm_uvsum2(response=c("os_time","os_status"), covs = c("age", "sex", "cohort"), data = pembrolizumab,family=NULL,offset=NULL,id=NULL,strata = NULL)
rm_uvsum(response=c("os_time","os_status"), covs = c("age", "sex", "cohort"), data = pembrolizumab,family=NULL,offset=NULL,id=NULL)
rm_uvsum2(response=c("os_time","os_status"), covs = c("age", "sex", "cohort"), data = pembrolizumab,family=NULL,offset=NULL,id=NULL)
rm_uvsum2(response=c("os_time","os_status"), covs = c("age", "sex", "cohort"), data = pembrolizumab,family=NULL,offset=NULL,id=NULL,strata = "")
rm_uvsum(response=c("os_time","os_status"), covs = c("age", "sex", "cohort"), data = pembrolizumab,family=NULL,offset=NULL,id=NULL,strata = "")
mv_crr <-crrRx(as.formula('os_time+os_status2~age+sex+cohort'),data=pembrolizumab)

rm_mvsum(mv_crr) # working now

mv_crr2 <-crrRx(as.formula('os_time+os_status2~age*sex+age:cohort'),data=pembrolizumab)
rm_mvsum2(mv_crr2)

rm_mvsum2(mv_lm2) # this is good

rm_mvsum(mv_binom2, showN = T, showEvent = T, vif = T, whichp = "both",tableOnly = T)

rm_mvsum2(mv_binom2, showN = T, showEvent = T, vif = T, whichp = "both")
summary(mv_binom2) # This is always correct
m_summary(mv_binom2) # This has wrong p-values for age:sex
coeffSum(mv_binom2) # This is right
gp(mv_binom2) # This is right

mv_pois2
rm_mvsum2(mv_pois2)
rm_mvsum(mv_negbin, showN = T, showEvent = T, vif = T, whichp = "both")
rm_mvsum2(mv_negbin, showN = T, showEvent = T, vif = T, whichp = "both")
rm_mvsum(mv_negbin2, showN = T, showEvent = T, vif = T, whichp = "both")
rm_mvsum2(mv_negbin2, showN = T, showEvent = T, vif = T, whichp = "both")

rm_mvsum(mv_ord, showN = T, showEvent = T, vif = T, whichp = "both",tableOnly = T)
rm_mvsum2(mv_ord, showN = T, showEvent = T, vif = T, whichp = "both")
coeffSum(mv_ord)
gp(mv_ord)
rm_mvsum(mv_ord2, showN = T, showEvent = T, vif = T, whichp = "both")
rm_mvsum2(mv_ord2, showN = T, showEvent = T, vif = T, whichp = "both")


rm_mvsum(mv_cox2, showN = T, showEvent = T, vif = T, whichp = "both",tableOnly = T)
rm_mvsum2(mv_cox2, showN = T, showEvent = T, vif = T, whichp = "both")
coeffSum(mv_cox2)
gp(mv_cox2)

rm_mvsum(mv_crr, showN = T, showEvent = T, vif = T, whichp = "both")

rm_mvsum2(mv_crr, showN = T, showEvent = T, vif = T, whichp = "both")
rm_mvsum2(mv_crr,whichp = "both")


rm_mvsum(mv_lm_gee, showN = T, showEvent = T, vif = T, whichp = "both")
rm_mvsum2(mv_lm_gee, showN = T, showEvent = T, vif = T, whichp = "both")

uv_lm_gee <- rm_uvsum("pdl1", covs = c("age", "sex", "cohort"),data = pembrolizumab,
                             id=pembrolizumab$id,
                             family =binomial, type="gee",corstr = "exchangeable")

uv_binom_gee <- rm_uvsum("pdl1", covs = c("age", "sex", "cohort"),data = pembrolizumab,
                        id=pembrolizumab$id, # always need to specify for gee models
                        family = binomial, type = "gee")
gp.crrRx(mv_crr2)

rm_uvsum2(response = c('os_time','os_status'),
           covs=c('age','sex','baseline_ctdna','l_size','change_ctdna_group'),
           data=pembrolizumab,CIwidth=.9,returnModels = F)

model <- survival::coxph(formula = as.formula(survival::Surv(os_time,
                                                             os_status) ~ change_ctdna_group), data = pembrolizumab)
gp(model)

pembrolizumab$change_ctdna_group
rm_uvsum2(response = c('os_time','os_status'),
          covs=c('change_ctdna_group'),
          data=pembrolizumab,CIwidth=.9)


#gee package
data("warpbreaks")
mv_gee <- gee::gee(breaks ~ tension, id=wool, data=warpbreaks, corstr="exchangeable")
class(mv_gee)
gp(mv_gee)
rm_mvsum2(mv_gee)
response = "breaks"
data=warpbreaks[sample(1:nrow(warpbreaks)),]
x_var="tension"
id="wool"
strata="";family="gaussian";offset=NULL; corstr = "independence"

test <- rm_compactsum(data = pembrolizumab, xvars = c("age",
                                              "change_ctdna_group2", "l_size", "pdl1"), grp = "sex", use_mean = "age",
              digits = c("age" = 2, "l_size" = 3), digits.cat = 1,tableOnly = T)
test[1,1] <-"**age** Mean (sd)"
test
outTable(test)


rm_mvsum(mv_binom2, showN = T, showEvent = T, vif = T, whichp = "both",table=T)
rm_mvsum2(mv_binom2, showN = T, showEvent = T, vif = T, whichp = "both")
coeffSum(mv_binom2)
gp(mv_binom2)


rm_uvsum2(response = c('os_time','os_status'),covs=c('age','sex','baseline_ctdna','l_size','change_ctdna_group'),
          data=pembrolizumab,CIwidth=.9)

row.names=NULL;to_indent=numeric(0);bold_headers=TRUE;
rows_bold=numeric(0);bold_cells=NULL;caption=NULL;digits=getOption("reportRmd.digits",2);
applyAttributes=TRUE;keep.rownames=FALSE; nicenames=TRUE;format=NULL

tab <- rm_mvsum2(mv_binom2, showN = T, showEvent = T, vif = T, whichp = "global",tableOnly = T)


out <- rm_mvsum2(mv_binom2, showN = T, showEvent = T, vif = T, whichp = "global")
rm_mvsum2(mv_binom2, showN = T, showEvent = T, vif = T, whichp = "levels")
m_summary(mv_binom2, whichp = "global")
gp(mv_ord)


drop1.default <- function (object, scope, scale = 0, test = c("none", "Chisq"),
          k = 2, trace = FALSE, ...)
{

  tl <- attr(terms(object), "term.labels")
  if (missing(scope))
    scope <- drop.scope(object)
  else {
    if (!is.character(scope))
      scope <- attr(terms(update.formula(object, scope)),
                    "term.labels")
    if (!all(match(scope, tl, 0L) > 0L))
      stop("scope is not a subset of term labels")
  }
  ns <- length(scope)
  ans <- matrix(nrow = ns + 1L, ncol = 2L, dimnames = list(c("<none>",
                                                             scope), c("df", "AIC")))
  ans[1, ] <- extractAIC(object, scale, k = k, ...)
  n0 <- nobs(object, use.fallback = TRUE)
  env <- environment(formula(object))
  for (i in seq_len(ns)) {
    tt <- scope[i]
    if (trace > 1) {
      cat("trying -", tt, "\n", sep = "")
      flush.console()
    }
    nfit <- update(object, as.formula(paste("~ . -", tt)),
                   evaluate = FALSE)
#    nfit <- eval(nfit, envir = env)
    nfit <- eval(nfit)
    ans[i + 1L, ] <- extractAIC(nfit, scale, k = k, ...)
    nnew <- nobs(nfit, use.fallback = TRUE)
    if (all(is.finite(c(n0, nnew))) && nnew != n0)
      stop("number of rows in use has changed: remove missing values?")
  }
  dfs <- ans[1L, 1L] - ans[, 1L]
  dfs[1L] <- NA
  aod <- data.frame(Df = dfs, AIC = ans[, 2])
  test <- match.arg(test)
  if (test == "Chisq") {
    dev <- ans[, 2L] - k * ans[, 1L]
    dev <- dev - dev[1L]
    dev[1L] <- NA
    nas <- !is.na(dev)
    P <- dev
    P[nas] <- safe_pchisq(dev[nas], dfs[nas], lower.tail = FALSE)
    aod[, c("LRT", "Pr(>Chi)")] <- list(dev, P)
  }
  head <- c("Single term deletions", "\nModel:", deparse(formula(object)),
            if (scale > 0) paste("\nscale: ", format(scale), "\n"))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}

# HIV data error
load("/Users/lisaavery/Library/CloudStorage/OneDrive-Personal/UofT HIV in Motion_NIH/combined_dataset.rda")

rm_compactsum(data=combined_df,
              xvars =c('Age', 'Sex', 'Gender', 'YearHIVDx', 'ARVS', 'CD4C', 'Aerobic_Exercise', 'Strength_Exercise', 'ExStatus', 'ExDays', 'MOS_EIS', 'MOS_TANGS', 'MOS_AFFS', 'MOS_PSI', 'MOS_Total', 'PHQ8_TOTAL','Pain', 'Mental_Health', 'Admin_Mode'),
              grp = "Study",pvalue=F)


xvars =c('Age', 'Sex', 'Gender', 'YearHIVDx', 'ARVS', 'CD4C', 'Aerobic_Exercise', 'Strength_Exercise', 'ExStatus', 'ExDays', 'MOS_EIS', 'MOS_TANGS', 'MOS_AFFS', 'MOS_PSI', 'MOS_Total', 'PHQ8_TOTAL','Pain', 'Mental_Health', 'Admin_Mode')

for (v in xvars){
  print(v)
  rm_compactsum(data=combined_df,
                xvars = v,
                grp = "Study",pvalue=F)
}

data <- combined_df
xvar <- "ARVS"
x = combined_df$ARVS
attributes(x) <- NULL
as.character(x)
is_binary <- function(x) {
  attributes(x) <- NULL
  all(unique(na.omit(x)) %in% c(0, 1))
}


lbl <- attr(data[[xvar]],"labels")
attributes(data[[xvar]]) <- NULL
newx <- factor(data[[xvar]],levels=lbl,labels=names(lbl))

table(newx,data[[xvar]])


all_pred <- c("Age","BMI",
              "Pregestational_diabetes","Thyroid_disorder","Chronic_hypertension",
              "Gravida","Para","GA_at_delivery_index_pregnancy",
              "More_than_1_previous_PPROM","Prior_Term_Pregnancy",
              "PPIP_Normal","PPIP_Chorio","PPIP_MVM_FVM",
              "ASA","Progesterone","ATB_suppression","Antepartum_bleed",
              "Ureaplasma_Mycoplasma_at_12_14_weeks","Ureaplasma_Mycoplasma_at_16_18_weeks",
              "Negative_UU_after_treatment","persistance_UU_after_treatment",
              "Positive_BV_swab","Need_for_cerclage",
              "cx_short_16","cx_short_20","cx_short_24",
              "CL_at_16_weeks","CL_at_20_weeks","CL_at_24_weeks")

pred_vars <- c("Age","BMI",
               "Pregestational_diabetes","Thyroid_disorder","Chronic_hypertension",
               "Gravida","Para","GA_at_delivery_index_pregnancy",
               "More_than_1_previous_PPROM","Prior_Term_Pregnancy",
               "PPIP_Normal","PPIP_Chorio","PPIP_MVM_FVM",
               "ASA","Progesterone","ATB_suppression","Antepartum_bleed",
               "Ureaplasma_Mycoplasma_at_12_14_weeks","Ureaplasma_Mycoplasma_at_16_18_weeks",
               "Negative_UU_after_treatment","persistance_UU_after_treatment",
               "Positive_BV_swab","Need_for_cerclage")

load(file="test_ws.rda")

rm_uvsum2(data=df_g1,response="Second_PPROM",
                                      covs=pred_vars,
                                      showEvent = F,tableOnly = T)

# Test plots ---------
 data("pembrolizumab")
 p <- ggplot(pembrolizumab,aes(x=change_ctdna_group,y=baseline_ctdna)) +
 geom_boxplot()
 replace_plot_labels(p)
 pembrolizumab <- set_var_labels(pembrolizumab,
 change_ctdna_group="Change in ctDNA group")
 p <- ggplot(pembrolizumab,aes(x=change_ctdna_group,y=baseline_ctdna)) +
 geom_boxplot()
 replace_plot_labels(p)
 # Can also be used with a pipe, but expression needs to be wrapped in a brace
 (ggplot(pembrolizumab,aes(x=change_ctdna_group,y=baseline_ctdna)) +
 geom_boxplot()) |> replace_plot_labels()

 p <- ggplot(pembrolizumab,aes(x=change_ctdna_group,y=baseline_ctdna)) +
   geom_boxplot() +
   facet_wrap(~sex)
 replace_plot_labels(p)


 # Test uv mv ----------
  require(survival)
  data("pembrolizumab")
  uvTab <- rm_uvsum2(response = c('os_time','os_status'),
  covs=c('age','sex','baseline_ctdna','l_size','change_ctdna_group'),
  data=pembrolizumab,tableOnly=TRUE)
  mv_surv_fit <- coxph(Surv(os_time,os_status)~age+sex+
  baseline_ctdna+l_size+change_ctdna_group, data=pembrolizumab)
  mvTab <- rm_mvsum2(mv_surv_fit,tableOnly = TRUE)
  rm_uv_mv(uvTab,mvTab)

  # linear model
  uvtab<-rm_uvsum(response = 'baseline_ctdna',
  covs=c('age','sex','l_size','pdl1','tmb'),
  data=pembrolizumab,tableOnly=TRUE)
  # pdl1 not included in model
  lm_fit=lm(baseline_ctdna~age+sex+l_size+tmb,data=pembrolizumab)
  mvtab<-rm_mvsum(lm_fit,tableOnly = TRUE)
  rm_uv_mv(uvtab,mvtab)

  #logistic model
  uvtab<-rm_uvsum(response = 'os_status',
  covs=c('age','sex','l_size','pdl1','tmb'),
  data=pembrolizumab,family = binomial,tableOnly=TRUE)
  logis_fit<-glm(os_status~age+sex+l_size+pdl1+tmb,data = pembrolizumab,family = 'binomial')
  mvtab<-rm_mvsum(logis_fit,tableOnly = TRUE)
  rm_uv_mv(uvtab,mvtab,tableOnly=TRUE)


call_str <- deparse(model$call)
call_str_vc <- as.character(model$call)

# Count the number of commas
num_commas <- length(unlist(gregexpr(",", call_str)))
offset_str <- call_str_vc[num_commas+2]
