#-------------------------------------------------------------------------------------
# Data Definition
data(cancer, package='survival')
data(dietox,package='geepack')

lung <- cancer
lung$Status=factor(lung$status-1)
lung$Sex = factor(lung$sex,labels = c('Male','Female'))
lung$AgeGroup = factor(cut(lung$age, breaks=seq(0,100,10)), labels = c('30s','40s','50s','60s','70s','80s'))
lung$OneLevelFactor = factor(x='one level')
lung = lung[order(lung$Status),]

lung$x_null = rnorm(nrow(lung))
lung$x_pred = c(rnorm(sum(lung$Status==0),0,1),
                rnorm(sum(lung$Status==1),1,1))


#data(dietox)
dietox$Cu     <- as.factor(dietox$Cu)

test_row=data.frame(type=c(1,NA,1,1,2,2,2,2,3,3),age=c(5,NA,10,24,35,12,45,34,NA,12),group=c('A','A',NA,'A','A','B','I','C','A','B'),
                group2=c('A','A','A','A','A','A','A','A','A','A'),group3=c('A','A','A','A','A','A','A','A','A',NA),group4=c('A','B','A','A','B','B','B','B',NA,NA))

lung_missing <- lung
set.seed(123)
lung_missing$ph.ecog[sample(c(1:228),50)] <- NA
lung_missing$AgeGroup[sample(c(1:228),50)] <- NA

#-------------------------------------------------------------------------------------
# NOTES
# As of R 4.4.0 profile likelihoods are automatically calculated for glm models
# This is true even for linear models fit with the glm function
# These tests have been updated 30 Dec 2024 to reflect these changes
#-------------------------------------------------------------------------------------

test_that("covsum calculates correctly with no maincov", {
  output = covsum(data=lung,
                  covs=c('Status','Sex','wt.loss','OneLevelFactor'))
  expect_equal(names(output) , c("Covariate",'n=228'))
  expect_equal(output$Covariate , c("Status","0","1","Sex","Male","Female","wt loss","Mean (sd)","Median (Min,Max)","Missing","OneLevelFactor","one level"))
  expect_equal(output$`n=228`,c("","63 (28)","165 (72)","","138 (61)","90 (39)","","9.8 (13.1)","7 (-24, 68)","14","","228 (100)"))

})

test_that("covsum calculates correctly with maincov", {
  output = covsum(data=lung,
                  maincov='Status',
                  covs=c('Sex','wt.loss','OneLevelFactor'))
  expect_equal(names(output) ,c("Covariate","Full Sample (n=228)","0 (n=63)","1 (n=165)","p-value") )
  expect_equal(output$Covariate , c("Sex","Male","Female","wt loss","Mean (sd)","Median (Min,Max)","Missing","OneLevelFactor","one level"))
  expect_equal(output[,3],c("","26 (41)","37 (59)","","9.1 (12.9)","4 (-10, 49)","1","","63 (100)"))

})

test_that("covsum rounds variables correctly", {
  output = covsum(data=lung,
                  covs=c('Status','Sex','wt.loss','OneLevelFactor','x_pred'),
                  digits=3,
                  digits.cat = 1)
  mean_sd = paste0(format(round(mean(lung$x_pred,na.rm = T),3),nsmall=3), " (",format(round(sd(lung$x_pred,na.rm = T),3),nsmall=3),")")
  med_min_max = paste0(format(round(median(lung$x_pred,na.rm = T),3),nsmall=3), " (",format(round(min(lung$x_pred,na.rm = T),3),nsmall=3),", ",format(round(max(lung$x_pred,na.rm = T),3),nsmall=3),")")
  expect_equal(names(output) , c("Covariate",'n=228'))
  expect_equal(output$Covariate , c("Status","0","1","Sex","Male","Female","wt loss","Mean (sd)","Median (Min,Max)","Missing","OneLevelFactor","one level","x pred","Mean (sd)",'Median (Min,Max)'))
  expect_equal(output$`n=228`,c("","63 (27.6)","165 (72.4)","","138 (60.5)","90 (39.5)","","9.832 (13.140)","7 (-24, 68)","14","","228 (100.0)","",mean_sd,med_min_max))
  output = covsum(data=lung,
                  covs=c('Status','Sex','wt.loss','OneLevelFactor','x_pred'),
                  digits=3,
                  digits.cat = 1,
                  IQR = T)
  med_Q1_Q3 = paste0(format(round(median(lung$x_pred,na.rm = T),3),nsmall=3), " (",format(round(quantile(lung$x_pred,.25,na.rm = T),3),nsmall=3),", ",format(round(quantile(lung$x_pred,.75,na.rm = T),3),nsmall=3),")")
  expect_equal(names(output) , c("Covariate",'n=228'))
  expect_equal(output$Covariate , c("Status","0","1","Sex","Male","Female","wt loss","Mean (sd)","Median (Q1,Q3)","Missing","OneLevelFactor","one level","x pred","Mean (sd)",'Median (Q1,Q3)'))
  expect_equal(output$`n=228`,c("","63 (27.6)","165 (72.4)","","138 (60.5)","90 (39.5)","","9.832 (13.140)","7.000 (0.000, 15.750)","14","","228 (100.0)","",mean_sd,med_Q1_Q3))

})

test_that("covsum calculates rows correctly", {
  output = covsum(data=test_row,
                  maincov='type',
                  covs=c('group','group2','group3','group4'),
                  percentage = 'row',testcat = 'Fisher')
  expect_equal(names(output) ,c("Covariate","Full Sample (n=9)","1 (n=3)","2 (n=4)",'3 (n=2)',"p-value") )
  expect_equal(output$Covariate , c("group", "A", "B", "C", "I", "Missing", "group2", "A", "group3",
                                    "A", "Missing", "group4", "A", "B", "Missing"))
  expect_equal(output[,3],c("", "2 (50)", "0 (0)", "0 (0)", "0 (0)", "1", "", "3 (33)",
                            "", "3 (38)", "0", "", "3 (100)", "0 (0)", "0"))

  expect_equal(output[,4],c("", "1 (25)", "1 (50)", "1 (100)", "1 (100)", "0", "", "4 (44)",
                            "", "4 (50)", "0", "", "0 (0)", "4 (100)", "0"))
  expect_equal(output[,5],c("", "1 (25)", "1 (50)", "0 (0)", "0 (0)", "0", "", "2 (22)",
                            "", "1 (12)", "1", "", "0 (0)", "0 (0)", "2"))

})

test_that("covsum includes missing correctly", {
  output = covsum(data=lung_missing,maincov = 'ph.ecog',
                  covs=c('Sex','AgeGroup','age','meal.cal'),include_missing = TRUE,pvalue = FALSE)

  expect_equal(names(output) ,c("Covariate", "Full Sample (n=228)", "0 (n=49)", "1 (n=86)",
                                "2 (n=41)", "3 (n=1)", "NA (n=51)") )

  expect_equal(output[,2],c("", "138 (61)", "90 (39)", "", "2 (1)", "17 (10)", "58 (33)",
                            "65 (37)", "34 (19)", "2 (1)", "50", "", "62.4 (9.1)",
                            "63 (39, 82)", "", "928.8 (402.2)", "975 (96, 2600)", "47"))

  expect_equal(output[,4],c("", "53 (62)", "33 (38)", "", "1 (1)", "6 (9)", "21 (31)",
                            "28 (42)", "11 (16)", "0 (0)", "19", "", "61.7 (9.2)",
                            "64 (40, 80)", "", "940.9 (343.6)", "975 (280, 2450)", "21"))

})

test_that("covsum includes missing correctly when presenting row percentages", {
  output = covsum(data=lung_missing,maincov = 'ph.ecog',
                  covs=c('Sex','AgeGroup','age','meal.cal'),include_missing = TRUE,pvalue = FALSE,percentage = 'row',digits.cat = 3)

  expect_equal(names(output) ,c("Covariate", "Full Sample (n=228)", "0 (n=49)", "1 (n=86)",
                                "2 (n=41)", "3 (n=1)", "NA (n=51)") )

  expect_equal(output[,2],c("", "138", "90", "", "2", "17", "58", "65", "34", "2", "50",
                            "", "62.4 (9.1)", "63 (39, 82)", "", "928.8 (402.2)", "975 (96, 2600)",
                            "47"))

  expect_equal(output[,4],c("", "53 (38.406)", "33 (36.667)", "", "1 (50.000)", "6 (35.294)",
                            "21 (36.207)", "28 (43.077)", "11 (32.353)", "0 (0.000)", "19",
                            "", "61.7 (9.2)", "64 (40, 80)", "", "940.9 (343.6)", "975 (280, 2450)",
                            "21"))

})

# # This runs ok, but not in the rmd check for some reason
# test_that("uvsum2 logistic regression CIS are correct",{
#   digits = 2 # TODO: add function flexibility
#   expected <- list(c("1.20 (0.94, 1.55)", "2.65 (1.98, 3.65)"),
#                    c("1.20 (0.89, 1.63)", "2.65 (1.87, 3.89)"),
#                    c("1.20 (0.81, 1.80)", "2.65 (1.69, 4.42)"),
#                    c("1.20 (0.79, 1.86)", "2.65 (1.63, 4.65)"))
#   ci_widths <- c(0.9,.95,.99,.995)
#   for (i in seq_along(ci_widths)){
#     ci_width <- ci_widths[i]
# #    print(ci_width)
#     # set.seed(1)
#     m1 = glm(Status~x_null,data=lung,family='binomial')
#     x1=summary(m1)$coefficients
#     m2 = glm(Status~x_pred,data=lung,family='binomial')
#     x2=summary(m2)$coefficients
#     digits <- 2
# # #    DONE: updated the wald CIs to likelihood profile CIs?
# #     expected = c(  paste(format(round(exp(x1[2,1]),digits),nsmall=digits),
# #                          paste0("(",paste0(format(round(exp(confint(m1,level=ci_width)[2,]),digits),nsmall=digits),collapse=", "),")")),
# #                    paste(format(round(exp(x2[2,1]),digits),nsmall=digits),
# #                           paste0("(",paste0(format(round(exp(confint(m2,level=ci_width)[2,]),digits),nsmall=digits),collapse=", "),")")))
# #     print(expected)
#     # OLD - WALD CIs
#     # expected = c(  paste(format(round(exp(x1[2,1]),digits),nsmall=digits),
#     #                      paste0("(",paste0(format(round(c(exp(x1[2,1]-qnorm(1-(1-ci_width)/2)*x1[2,2]),exp(x1[2,1]+qnorm(1-(1-ci_width)/2)*x1[2,2])),digits),nsmall=digits),collapse=", "),")")),
#     #                paste(format(round(exp(x2[2,1]),digits),nsmall=digits),
#     #                      paste0("(",paste0(format(round(c(exp(x2[2,1]-qnorm(1-(1-ci_width)/2)*x2[2,2]),exp(x2[2,1]+qnorm(1-(1-ci_width)/2)*x2[2,2])),digits),nsmall=digits),collapse=", "),")")))
#     set.seed(1)
#     output = uvsum2(response = 'Status',
#                    covs=c('x_null','x_pred'),
#                    data=lung,
#                    type='logistic',
#                    CIwidth = ci_width)
#
#     expect_equal(output[,2],expected[[i]])
#     expect_equal(names(output)[2],paste0("OR(",ci_width*100,"%CI)"))
#
#   }
#
# })

test_that("uvsum2 linear regression CIS are correct",{
  digits = 2 # TODO: add function flexibility
  for (ci_width in c(0.9,.95,.99,.995)){
    m1 = lm(age~wt.loss,data=lung)
    x1=summary(m1)$coefficients
    m2 = lm(age~Sex,data=lung)
    x2=summary(m2)$coefficients
    expected = c(  paste(format(round(x1[2,1],digits),nsmall=digits),
                         paste0("(",paste0(format(round(c(x1[2,1]-qt(1-(1-ci_width)/2,m1$df.residual)*x1[2,2],x1[2,1]+qt(1-(1-ci_width)/2,m1$df.residual)*x1[2,2]),digits),nsmall=digits),collapse=", "),")")),
                   paste(format(round(x2[2,1],digits),nsmall=digits),
                         paste0("(",paste0(format(round(c(x2[2,1]-qt(1-(1-ci_width)/2,m2$df.residual)*x2[2,2],x2[2,1]+qt(1-(1-ci_width)/2,m2$df.residual)*x2[2,2]),digits),nsmall=digits),collapse=", "),")")))
    expected = gsub(',  ',', ',expected)
    output = uvsum2(response = 'age',
                   covs=c('wt.loss','Sex'),
                   data=lung,
                   type='linear',
                   CIwidth = ci_width)

    expect_equal(output[c(1,4),2],expected)
    expect_equal(names(output)[2],paste0("Estimate(",ci_width*100,"%CI)"))
  }

})

test_that("rm_mvsum outputs lm models correctly",{
  fit_lm = lm(wt.loss~age+Sex,data=lung)
    output = rm_mvsum(fit_lm,showN = FALSE,tableOnly = TRUE)
    covs = c('age','Sex','Male','Female')
    est = c('0.03 (-0.16, 0.23)',NA,'Reference','-3.38 (-7.00, 0.25)')
    names = c('Covariate','Estimate(95%CI)','p-value','VIF')
    expect_equal(output[,1],covs)
    expect_equal(output[,2],est)
    expect_equal(names(output),names)
})

test_that("rm_mvsum outputs glm linear models correctly",{
  fit_glm = glm(wt.loss~age+Sex,data=lung)
  output = rm_mvsum(fit_glm,CIwidth = .99,showN = F,tableOnly = T,vif=T)
  names = c('Covariate','Estimate(99%CI)','p-value','VIF')
  covs = c('age','Sex','Male','Female')
  est = c('0.03 (-0.22, 0.29)',NA,'Reference','-3.38 (-8.11, 1.36)')
  expect_equal(output[,1],covs)
  expect_equal(output[,2],est)
  expect_equal(names(output),names)
})

test_that("rm_mvsum outputs glm binomial models correctly",{
  fit_glm = glm(Status~age+Sex,data=lung,family='binomial')
  output = rm_mvsum(fit_glm,vif=FALSE,showN = F,showEvent = F,tableOnly = T)
  names = c('Covariate','OR(95%CI)','p-value')
  covs = c('age','Sex','Male','Female')
  est = c('1.03 (1.00, 1.07)',NA,'Reference','0.35 (0.19, 0.64)')
  expect_equal(output[,1],covs)
  expect_equal(output[,2],est)
  expect_equal(names(output),names)
})

test_that("rm_mvsum outputs glm poisson models correctly",{
  fit_glm = glm(status~age+Sex,data=lung,family='poisson')
  output = rm_mvsum(fit_glm,vif=FALSE,showN = F,showEvent = F,tableOnly = T)
  names = c('Covariate','RR(95%CI)','p-value')
  covs = c('age','Sex','Male','Female')
  est = c('1.00 (0.99, 1.01)',NA,'Reference','0.88 (0.72, 1.09)')
  expect_equal(output[,1],covs)
  expect_equal(output[,2],est)
  expect_equal(names(output),names)
})

test_that("uvsum2 outputs geeglm models correctly",{
  # mf1 <- formula(Weight ~ Cu)
  # gee1 <- geeglm(mf1, data=dietox, id=Pig, family=gaussian("identity"), corstr="ar1")
  # mf2 <- formula(Weight ~ Time)
  # gee2 <- geeglm(mf2, data=dietox, id=Pig, family=gaussian("identity"), corstr="ar1")
  output = uvsum2(response = 'Weight',covs=c('Cu','Time'), gee=T,
                 data=dietox, id='Pig', family=gaussian("identity"), corstr="ar1",whichp='both')
  names = c('Variable','Estimate(95%CI)','p-value','N')
  covs = c('Cu','Cu000','Cu035','Cu175','Time')
  est = c(NA,"Reference","-0.49 (-3.50, 2.52)","1.78 (-1.89, 5.45)", "6.73 (6.58, 6.88)")
  expect_equal(output[,1],covs)
  expect_equal(output[,2],est)
  expect_equal(names(output),names)
})

test_that("rm_mvsum outputs geeglm models correctly",{
  mf1 <- formula(Weight ~ Cu+Time)
  gee1 <- geepack::geeglm(mf1, data=dietox, id=Pig, family=gaussian("identity"), corstr="ar1")
  output = rm_mvsum(gee1,data=dietox,showN = T,tableOnly = T,whichp='both',vif=F)
  names = c('Covariate','Estimate(95%CI)','p-value','N')
  covs = c('Cu','Cu000','Cu035','Cu175','Time')
  est = c(NA,"Reference","-0.47 (-3.29, 2.35)","1.21 (-2.30, 4.71)", "6.73 (6.58, 6.88)")
  expect_equal(output[,1],covs)
  expect_equal(output[,2],est)
  expect_equal(names(output),names)
})

# # Uncomment this if you need to ensure the tests are being run
# test_that("this script is being executed",{
#   expect_equal("This script was run","YES!")
# })
# TODO: Finishing writing CI test scripts for uvsum2 (boxcos, coxph, crr)
# Write script to check sample size calculations

# TODO: test_that('CI works correctly ')

#lung$sex <- ifelse(lung$sex==1,"F","M")
# lung$crr_status <- lung$status*1
# lung$crr_status[seq(1,200,by=2)] <- 2

#mvsum2(crrRx(time+crr_status~age+sex,data=lung),data=lung)
#mvsum2(coxph(Surv(time,status)~age+sex,data=lung))

#uvsum2(response = c('time','crr_status'),covs=c('age','sex'),data=lung)
#uvsum2(response = c('time','status'),covs=c('age','sex'),data=lung)

