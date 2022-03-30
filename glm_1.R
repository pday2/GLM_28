# GLM - Homicide Victims
# Group 28
# Due: 30 March 2022
# Peter Day r0866276
# Augustin Kuntz r0879499
# Qusai A M Iwidat r0726650
# Liam Stein r0876581
# Kylan Young r0711789
# Martha Efeti Tondo r0886202

# Change your setwd
# if (Sys.info()[1] == "Windows") {
#   setwd("C:/Users/peter/My Tresors/Documentsacer/KULeuven/GLM/Project")
# } else {
#   setwd("/home/muddy/Tresors/Documentsacer/KULeuven/GLM/Project")
# }
setwd("C:/Workdir/GLM/GLM_28")

library("tidyverse")
library("dplyr")
library("ggplot2")
library("car") # For Anova
library("aod") # for wald.test
library("MASS") # for glm.nb - negative binomial GLM
library("DHARMa") #simulateResiduals
library("pscl") # for zeroinfl model
# glm.RR from https://rpubs.com/kaz_yos/poisson
glm.RR <- function(GLM.RESULT, digits = 3) {
  if (GLM.RESULT$family$family == "binomial") {
    LABEL <- "OR"
  } else if (GLM.RESULT$family$family == "poisson") {
    LABEL <- "RR"
  } else {
    stop("Not logistic or Poisson model")
  }
  COEF      <- stats::coef(GLM.RESULT)
  CONFINT   <- stats::confint(GLM.RESULT)
  TABLE     <- cbind(coef=COEF, CONFINT)
  TABLE.EXP <- round(exp(TABLE), digits)
  colnames(TABLE.EXP)[1] <- LABEL
  TABLE.EXP
}

hosmerlem = function(y, yhat, g=10) {
  fcutyhat = cut(yhat,breaks = quantile(yhat, probs=seq(0,1, 1/g)), include.lowest=TRUE)
  obs = xtabs(cbind(1 - y, y) ~ cutyhat)
  expect = xtabs(cbind(1 - yhat, yhat) ~ cutyhat)
  chisq = sum((obs - expect)^2/expect)
  P = 1 - pchisq(chisq, g - 2)
  return(list(chisq=chisq,p.value=P))
}


# 1. resp: The number of victims the respondent knows
# 2. race: The race of the respondent (black or white)
data <- read.csv("homicide.csv", header = TRUE, stringsAsFactors = TRUE)
data <- data %>% dplyr::select(c(resp, race))
# Exploratory
summary(data)
mean(data$resp) # 0.144
var(data$resp) # 0.295 - Overdispersion may be an issue as Var(y) > E(y)
ggplot(data=data, aes(x=resp)) + geom_histogram(binwidth=1) 
ggplot(data, aes(x=resp, fill=race)) + geom_bar(position='dodge', stat='count')
# counts for each race by number of known homicides
(tab <- table(data$resp, data$race))
table(data)
round(prop.table(tab), 3)
#data[which(data$race=="black"),]
#data[which(data$race=="white"),]
#data$resp <- as.factor(data$resp)

# 1. Fit a Poisson model
poisFit <- glm(resp ~ race, data = data, family = poisson(link = "log"))
summary(poisFit)
poisFit$coefficients

# Relative Risk see slide 450
cbind(Est=poisFit$coefficients, confint(poisFit))
# Risk Ratio slide 462
# 2. Risk ratio and the corresponding confidence interval.
round(exp(cbind(RR=poisFit$coefficients, confint(poisFit))), 3)

# The average number of known homicide victims for Whites is
#     0.177 times that of Blacks


# or use glm.RR - same numbers as above
glm.RR(poisFit)
# 3. Calculate the ratio of the means of the response for each race
#   (mean response for black/mean response for white)
mrwhite <- mean(predict.glm(poisFit, data[which(data$race=="white"),]))
mrblack <- mean(predict.glm(poisFit, data[which(data$race=="black"),]))
ratioBW <- mrblack/mrwhite
ratioBW #ratio = 0.27

# 4. Calculate the predictions of the models for each rate (or race???)
# Predict case per person (n = 1) 
exp(predict(poisFit, newdata = data.frame(race = "white", n = 1)))
exp(predict(poisFit, newdata = data.frame(race = "black", n = 1)))

# 5. Analyze the GOF of the model (deviance test, over-dispersion...)

## GOF tests
# Pearson test
X2=sum(residuals(poisFit, type = "pearson")^2)
n=dim(data)[1]
p=length(coef(poisFit))
data.frame(X2s=X2,pvalue=(1-pchisq(X2,n-p))) #p-value not good
# X2 = 2279, pval = 0 - large lack of fit, model not good

# Deviance tests
Dev=summary(poisFit)$deviance
df=summary(poisFit)$df.residual
data.frame(Dev=Dev, df=df, pvalue=(1-pchisq(Dev,df))) #P-value good, deviation is relatively close to df
# dev=844 < df=1306, p-value=1 => FtR H0 model is good fit???? hmmmmm


# Likelihood Ratio Test
Anova(poisFit, test="LR", type=3) #Looks good, significant difference between races

# Wald Test
Anova(poisFit, test="Wald", type=3) #Looks good
# or
wald.test(Sigma=vcov(poisFit), b=coef(poisFit), Terms = 1:2)
# p-value < 0.001 => evidence for association between race and resp, good fit



#Residuals diagnostics
sim.pois <- simulateResiduals(poisFit,plot=T)
hist(sim.pois)
testUniformity(sim.pois)
testDispersion(sim.pois) #Overdispersion issue

# Hosmer-Lemeshow test
# HLPois <- hosmerlem(y=data$resp, yhat=fitted(poisFit), g=10)
# HL test returns error for all g except g=1 which results in other errors
# Error in cut.default(yhat, breaks = quantile(yhat, probs = seq(0, 1, 1/g)),
#   : 'breaks' are not unique 

# 6. Fit a negative binomial model and get estimated model based variances
#    (per race) for the counts. Compare them with the observed 

# Fit NB model
nbFit <- glm.nb(resp ~ race, data )
summary(nbFit)
nbFit$coefficients #less effect for racewhite
nbFit$theta

# Model based variances
var_black_nb <- exp(nbFit$coefficients[1]) + (1/nbFit$theta)*exp(nbFit$coefficients[1])^2
var_white_nb <- exp(nbFit$coefficients[1] + nbFit$coefficients[2]) +(1/nbFit$theta)*exp(nbFit$coefficients[1] + nbFit$coefficients[2])^2

# Observed variances
var_black_obs <- var(data[which(data$race == "black"),1])
var_white_obs <- var(data[which(data$race == "white"),1])

# comparing the variances
results = data.frame(Race = c("black","white"), NB_Variance = c(var_black_nb,var_white_nb), Obs_Variance = c(var_black_obs,var_white_obs) )
print(results)


# 7. Fit a Quasi-likelihood model

QLM <-glm(resp ~ race, data = data, family =quasipoisson)
summary(QLM)


# 8. Conclusion
aicp <- AIC(poisFit) # 1121
aicn <- AIC(nbFit) # 1001 - nb better than pois
aicq <- AIC(QLM) # fails because q-l does not have a likelihood
llp <- logLik(poisFit) # -558
lln <- logLik(nbFit) # -497 again nb better than pois
llq <- logLik(QLM) # fails because q-l does not have a likelihood
X2p <- sum(residuals(poisFit, type = "pearson")^2)
p=length(coef(poisFit))
data.frame(X2s=X2p,pvalue=(1-pchisq(X2p,n-p)))
X2n <- sum(residuals(nbFit, type = "pearson")^2)
p=length(coef(poisFit))
data.frame(X2s=X2n,pvalue=(1-pchisq(X2n,n-p)))
X2q <- sum(residuals(QLM, type = "pearson")^2)
p=length(coef(poisFit))
data.frame(X2s=X2q,pvalue=(1-pchisq(X2q,n-p)))
summPois <- summary(poisFit)
summNB <- summary(nbFit)
summQLM <- summary(QLM)
# phi - scale or estimate of phi
summPois$dispersion
summNB$dispersion
summQLM$dispersion
# deviance(fit): Returns the residual deviance D(y, μ̂) for the fitted glm
# deviance closer to 0 is better
devp <- deviance(poisFit)
devn <- deviance(nbFit)
devq <- deviance(QLM)

tab <- matrix(round(c(aicp, aicn, aicq,
                llp, lln, llq, 
                devp, devn, devq,
                X2p, X2n, X2q), 3), ncol = 4, byrow = FALSE)
colnames(tab) <- c('AIC','logLike','dev', 'Chi2')
rownames(tab) <- c('Poisson','NegBin','QLM')
tab <- as.table(tab)
tab


# Just want to check Zero Inflated Poisson since there are many 0 responses
zip <- zeroinfl(resp ~ race | race, data = data)
summary(zip)
zinb <- zeroinfl(resp ~ race | race, data = data, dist = "negbin")
summary(zinb)
aiczip <- AIC(zip) # 998.7 only just slightly better than NegBin at 1001.8
aiczin <- AIC(zinb) # 999.0


