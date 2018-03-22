library(mice)
library(lme4)
library(mvtnorm)
library(lmerTest)

Daily <- read.csv("C:/Users/Josie/Dropbox/DesignedMissingness/Daily with Tox.csv")

names(Daily)[1] <- c('PID')
Daily <- Daily[,c('PID', 'Day', 'Q1', 'Q2', 'Q3', 'Q4', 'Q5', 'Q6', 'Q7',
                  'MissedDose', 'ZEduc', 'ZIncom45', 'Gender','Age', 
                  'ZAlcTox','ZCESDFU', 'ZAUDIT', 'DrinkYN')]

comp <- complete(mice(Daily, m=1, seed=817236),1)
comp$Age <- as.numeric(comp$Age)
comp$Day <- comp$Day-1
rm(Daily)

comp$Age <- as.numeric(scale(comp$Age))
comp$Q1 <- as.factor(comp$Q1)
comp$Q2 <- as.factor(comp$Q2)
comp$Q3 <- as.factor(comp$Q3)
comp$Q4 <- as.factor(comp$Q4)
comp$Q5 <- as.factor(comp$Q5)
comp$Q6 <- as.factor(comp$Q6)
comp$Q7 <- as.integer(comp$Q7)
comp$DrinkYN <- as.factor(comp$DrinkYN)
comp$MissedDose <- as.factor(comp$MissedDose)
comp$Gender <- as.factor(comp$Gender)

fitMD <- glm(data=comp, MissedDose ~ DrinkYN + ZAlcTox + Day, family=binomial(link=logit))
fitranMD <- glmer(data=comp, MissedDose ~ DrinkYN + ZAlcTox + Day + (1|PID), family=binomial(link=logit))

#Find fit to simulate
fit1 <- glm(data=comp, DrinkYN~ZIncom45+Age+ZAlcTox+ZAUDIT, family=binomial(link=logit))
fit2<- glm(data=comp, Q1~ZEduc+ZIncom45+ZAlcTox+ZCESDFU, family=binomial(link=logit))
fit3<- glm(data=comp, Q2~Q1+ZEduc+ZIncom45+Age+ZAlcTox
           +ZCESDFU+ZAUDIT+Gender, family=binomial(link=logit))
fit4<- glm(data=comp, Q3~Q1+Q2+ZEduc+ZIncom45+ZAlcTox
           +ZCESDFU+Gender, family=binomial(link=logit))
fit5<- glm(data=comp, Q4~Q1+Q2+Q3+ZEduc+ZAlcTox
           +ZCESDFU+ZAUDIT+Gender, family=binomial(link=logit))
fit6<- glm(data=comp, Q5~Q1+Q3+Q4+ZEduc+ZIncom45+ZAlcTox
           +ZAUDIT, family=binomial(link=logit))
fit7<- glm(data=comp, Q6~Q3+Q5+DrinkYN+ZEduc+ZIncom45+Age+ZAlcTox
           +ZCESDFU+ZAUDIT, family=binomial(link=logit))

