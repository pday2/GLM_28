# GLM - Homicide Victims
# Group 28
# Due: 30 March 2022
# Peter Day r0866276
# Augustin Kuntz r0879499
# Qusai A M Iwidat r0726650
# Liam Stein r0876581
# Kylan Young r0711789
# Martha Efeti Tondo r0886202


library("tidyverse")
library("dplyr")
library("ggplot2")
library("car") # For Anova
library("aod") # for wald.test
library("MASS") # for glm.nb - negative binomial GLM

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

# Change your setwd
if (Sys.info()[1] == "Windows") {
  setwd("C:/Users/peter/My Tresors/Documentsacer/KULeuven/GLM/Project")
} else {
  setwd("/home/muddy/Tresors/Documentsacer/KULeuven/GLM/Project")   
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

(lambda_B <- exp(poisFit$coefficients[[1]]))
# hmmmmmmmmmmm not quite???????????????????
(lambda_W <- exp(poisFit$coefficients[[1]] + poisFit$coefficients[[2]]))
# I think this needs fixing!!!!!!!!!!!!!!!!
(ratioBW <- lambda_B/lambda_W)


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
data.frame(X2s=X2,pvalue=(1-pchisq(X2,n-p)))

# Deviance tests
Dev=summary(poisFit)$deviance
df=summary(poisFit)$df.residual
data.frame(Dev=Dev, df=df, pvalue=(1-pchisq(Dev,df)))

# Likelihood Ratio Test
Anova(poisFit, test="LR", type=3)

# Wald Test
Anova(poisFit, test="Wald", type=3)
# or
wald.test(Sigma=vcov(poisFit), b=coef(poisFit), Terms = 1:2)


# From https://rpubs.com/kaz_yos/poisson
#      Goodness of fit test If the residual deviance is close enough to the 
#      residual degrees of freedom, it is a good fit. It can be tested by 
#      Chi-squared test.
cbind(residual.deviance           = deviance(poisFit),
     residual.degrees.of.freedom = df.residual(poisFit),
     chisq.p.value               = pchisq(deviance(poisFit), df.residual(poisFit), lower = F)
)


# 6. Fit a negative binomial model and get estimated model based variances
#    (per race) for the counts. Compare them with the observed 

nbFit <- glm.nb(resp ~ race, data = data)
summary(nbFit)

# 7. Fit a Quasi-likelihood model

QLM <-glm(resp ~ race, data = data, family =quasipoisson)
summary(QLM)









library(DHARMa)
# citation("DHARMa")
# https://cran.microsoft.com/web/packages/DHARMa/vignettes/DHARMa.html


# glm.nb in library("MASS")


library("pscl") # for zeroinfl model
