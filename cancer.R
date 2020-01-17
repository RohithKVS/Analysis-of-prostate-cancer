#The following is the R code for building a linear model

#read the data
prostate = read.csv("prostate_cancer(1).csv",sep=",",header=T)

#As vesinv is a qualitative variable, we must convert it into a factor so that R treats it as a qualitative variable
prostate$vesinv=as.factor(prostate$vesinv)

#Get boxplot of the response variable
boxplot(prostate$psa,main="Boxplot of the original PSA")

#As we have outliers in the response, we use natural log transformation of response
prostate$psa=log(prostate$psa)
#Get boxplot of the transformed response
boxplot(prostate$psa,main="Boxplot of the transformed PSA")

#Look at the relationship between response and each predictor
plot(prostate$cancervol,prostate$psa, main="Plot of PSA against cancervol")
fit1 <- lm(prostate$psa ~ prostate$cancervol, data = prostate)
abline(fit1)
plot(prostate$weight,prostate$psa, main="Plot of PSA against weight")
fit2 <- lm(prostate$psa ~ prostate$weight, data = prostate)
abline(fit2)
plot(prostate$age,prostate$psa, main="Plot of PSA against age")
fit3 <- lm(prostate$psa ~ prostate$age, data = prostate)
abline(fit3)
plot(prostate$benpros,prostate$psa, main="Plot of PSA against benpros")
fit4 <- lm(prostate$psa ~ prostate$benpros, data = prostate)
abline(fit4)
plot(prostate$vesinv,prostate$psa, main="Plot of PSA against vesinv")
fit5 <- lm(prostate$psa ~ prostate$vesinv, data = prostate)
abline(fit5)
plot(prostate$capspen,prostate$psa, main="Plot of PSA against capspen")
fit6 <- lm(prostate$psa ~ prostate$capspen, data = prostate)
abline(fit6)
plot(prostate$gleason,prostate$psa, main="Plot of PSA against gleason")
fit7 <- lm(prostate$psa ~ prostate$gleason, data = prostate)
abline(fit7)

#Do a forward selection AIC to get a model and look at the model.
fit8.forward <- step(lm(prostate$psa ~ 1, data = prostate),scope = list(upper = ~prostate$cancervol + prostate$weight + prostate$age + prostate$benpros + prostate$vesinv + prostate$capspen + prostate$gleason),direction = "forward")
summary(fit8.forward)

#Do a backward selection AIC to get a model and look at the model.
fit9.backward <- step(lm(prostate$psa ~ prostate$cancervol + prostate$weight + prostate$age + prostate$benpros + prostate$vesinv + prostate$capspen + prostate$gleason, data = prostate),scope = list(lower = ~1), direction = "backward")
#Look at the model
summary(fit9.backward)

#Do both the selections to get a model and look at the model.
fit10.both <- step(lm(prostate$psa ~ 1, data = prostate),scope = list(lower = ~1, upper = ~prostate$cancervol + prostate$weight + prostate$age + prostate$benpros + prostate$vesinv + prostate$capspen + prostate$gleason),direction = "both")
summary(fit10.both)

#Build a new model based on the predictors provided by AIC and look at the model.
fit11 = lm(prostate$psa ~ prostate$cancervol + prostate$gleason + prostate$benpros + prostate$vesinv, data = prostate)
summary(fit11)

#Build a model based on all the predictors and look at the model.
fit12 = lm(prostate$psa ~ + prostate$cancervol + prostate$gleason + prostate$benpros + prostate$capspen + prostate$vesinv + prostate$age + prostate$weight, data = prostate)
summary(fit12)

#This function checks which model is better
anova(fit11,fit12)

#Check the most frequent category of the qualitative variable(vesinv)
if(length(subset(prostate$vesinv, prostate$vesinv == 0)) > length(subset(prostate$vesinv, prostate$vesinv == 1)))
{
  z1=0
}else{
  z1=1
}

#The below steps perform model diagnostics of the model

#residual plot
plot(fitted(fit11), resid(fit11),main="Residual plot of the model")
abline(h = 0)

# normal QQ plot
qqnorm(resid(fit11),main="QQ plot of the model")
qqline(resid(fit11))

#Get the estimates from the final model
estimate=fit11$coefficients

#Estimate the psa level whose quantitative predictors are at the sample means of the variables and qualitative predictors are at the most frequent category.
#As we have transformed the response variable, we must convert it back to the original response. This is done by taking exponential
psa_estimate=exp(estimate["(Intercept)"]+(estimate["prostate$cancervol"]*mean(prostate$cancervol))+(estimate["prostate$gleason"]*mean(prostate$gleason))+(estimate["prostate$benpros"]*mean(prostate$benpros))+(estimate["prostate$vesinv1"]*z1))

#Display the estimate
psa_estimate
