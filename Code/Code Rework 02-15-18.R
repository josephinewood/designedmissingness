
# Install needed packages, if necessary
# install.packages(c('lme4','mvtnorm','lmerTest'))

# Set working directory
setwd("~/Dropbox/Research/DesignedMissingness/")

# Load necessary libraries
library(mice)
library(lme4)
library(mvtnorm)
library(lmerTest)

# Read in Daily With Tox CSV data file
Daily <- read.csv('~/Dropbox/Research/DesignedMissingness/Daily with Tox.csv')

# Retain desired variables
names(Daily)[1] <- c('PID')
Daily <- Daily[,c('PID', 'Day', 'Q1', 'Q2', 'Q3', 'Q4', 'Q5', 'Q6', 'Q7'
                  ,'MissedDose'
                  , 'ZEduc', 'ZIncom45', 'Gender', 'Age', 'ZAlcTox'
                  ,'ZCESDFU', 'ZAUDIT', 'DrinkYN')]

# Impute missing values in new dataset, 'comp'...Daily has 4.74% missing (sum(is.na(Daily))/prod(dim(Daily)))
comp <- complete(mice(Daily, m=1, seed=817236),1)
# Convert Age to data type 'double'
comp$Age <- as.numeric(comp$Age)
# Shift days back 1 day
comp$Day <- comp$Day-1
# Remove the Daily dataset from the environment
rm(Daily)

#comp$Age <- as.numeric(scale(comp$Age)) # why???
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

# Correlation Matrix
# Use model.matrix to convert factor variables to dummy encoding
# Matches OG matrix!!!
corrmat <- cor(model.matrix(~.-1,data=comp[,c('ZEduc', 'ZIncom45', 'Age', 'ZAlcTox', 'ZCESDFU', 'ZAUDIT')]))
round(corrmat, 2) # display the matrix, rounding values to 2 decimal places

# Fit logistic models to each variable to use as foundation for simulation
# OG author built models that were seemingly random? New models are the target~all others
# Drink_fit <- glm(data=comp, DrinkYN~., family=binomial(link=logit)) # warning: fitted probs numerically 0 or 1 occurred
# 
# Q1_fit <- glm(data=comp, Q1~., family=binomial(link=logit))
# 
# Q2_fit <- glm(data=comp, Q2~., family=binomial(link=logit))
# 
# Q3_fit <- glm(data=comp, Q3~., family=binomial(link=logit))
# 
# Q4_fit <- glm(data=comp, Q4~., family=binomial(link=logit))
# 
# Q5_fit <- glm(data=comp, Q5~., family=binomial(link=logit))
# 
# Q6_fit <- glm(data=comp, Q6~., family=binomial(link=logit))

# OG stopped with Q6...but don't we want to simulate all of the variables? Nvm done with the corrmat
# Q7_fit <- lm(data=comp, Q7~.)
# 
# MD_fit <- lm(data=comp, MissedDose~.)
# 
# Edu_fit <- lm(data=comp, ZEduc~.)
# 
# Inc_fit <- lm(data=comp, ZIncom45~.)
# 
# Gender_fit <- lm(data=comp, Gender~.)
# 
# Age_fit <- lm(data=comp, Q6~.)
# 
# Alc_fit <- lm(data=comp, ZAlcTox~.)
# 
# Cesdfu_fit <- lm(data=comp, ZCESDFU~.)
# 
# Audit_fit <- lm(data=comp, ZAUDIT~.)

# for (n in names(comp[,-1])){
#   n <- as.name(n)
#   nam <- paste(n, "_fit", sep = "")
#   mod <- glm(data=comp,n~.-1,family=binomial(link=logit))
#   assign(nam, mod)
# }

# OG fit models for the questions using variables simulated via corrmat
fit1 <- glm(data=comp, DrinkYN~ZIncom45+Age+ZAlcTox+ZAUDIT, family=binomial(link=logit))
fit2<- glm(data=comp, Q1~ZEduc+ZIncom45+ZAlcTox+ZCESDFU, family=binomial(link=logit))
fit3<- glm(data=comp, Q2~Q1+ZEduc+ZIncom45+Age+ZAlcTox+ZCESDFU+ZAUDIT+Gender, family=binomial(link=logit))
fit4<- glm(data=comp, Q3~Q1+Q2+ZEduc+ZIncom45+ZAlcTox+ZCESDFU+Gender, family=binomial(link=logit))
fit5<- glm(data=comp, Q4~Q1+Q2+Q3+ZEduc+ZAlcTox+ZCESDFU+ZAUDIT+Gender, family=binomial(link=logit))
fit6<- glm(data=comp, Q5~Q1+Q3+Q4+ZEduc+ZIncom45+ZAlcTox+ZAUDIT, family=binomial(link=logit))
fit7<- glm(data=comp, Q6~Q3+Q5+DrinkYN+ZEduc+ZIncom45+Age+ZAlcTox+ZCESDFU+ZAUDIT, family=binomial(link=logit))

# Simulate
set.seed(1234)
sim.data <- list()
for (i in 1:200) {
  X <- as.data.frame(rmvnorm(60, mean=rep(0,6), sigma=corrmat))
  names(X) <- c('ZEduc', 'ZIncom45', 'Age', 'ZAlcTox', 'ZCESDFU', 'ZAUDIT')
  X$PID <- 1:60
  X$Gender <- as.factor(rbinom(60, 1, prob=0.5))
  X <- X[rep(row.names(X), 45),] #Repeat each row 45 times, all current variables are time invariant (ZAlcTox, too??)
  X <- X[order(X$PID),] #Order by PID
  X$Day <- rep(0:44, 60) #Create a variable containing day variable 0-44

  # Create Variables
  X$DrinkYN <- 1/(1+exp(-1*predict(fit1, newdata = X)))
  X$DrinkYN <- as.factor(rbinom(2700, 1, prob = X$DrinkYN))
  X$Q1 <- 1/(1+exp(-1*predict(fit2, newdata=X)))
  X$Q1 <- as.factor(rbinom(2700, 1, prob = X$Q1))
  X$Q2 <- 1/(1+exp(-1*predict(fit3, newdata = X)))
  X$Q2 <- as.factor(rbinom(2700,1,prob=X$Q2))
  X$Q3 <- 1/(1+exp(-1*predict(fit4, newdata = X)))
  X$Q3 <- as.factor(rbinom(2700,1,prob=X$Q3))
  X$Q4 <- 1/(1+exp(-1*predict(fit5, newdata = X)))
  X$Q4 <- as.factor(rbinom(2700,1,prob=X$Q4))
  X$Q5 <- 1/(1+exp(-1*predict(fit6, newdata = X)))
  X$Q5 <- as.factor(rbinom(2700,1,prob=X$Q5))
  X$Q6 <- 1/(1+exp(-1*predict(fit7, newdata = X)))
  X$Q6 <- as.factor(rbinom(2700,1,prob=X$Q6))
  X$Q7 <- rpois(2700, 2) # lambda = 2...mean of 2. why???
  X[X$DrinkYN==0,c('Q7')] <- 0
  
  # Make Missed Dose
  # This is where issues come in!!!
  # Closest I could get to these numbers was:
  # MD <- glm(data=comp,MissedDose~ZEduc+ZIncom45+Age+ZAlcTox+ZCESDFU+ZAUDIT,family=binomial(link=logit))
  # MD <- glm(data=comp,MissedDose~DrinkYN+ZAlcTox+Day,family=binomial(link=logit))
  MD <- glmer(MissedDose ~ DrinkYN*Day + ZAlcTox + (1|PID), data = comp, family=binomial(link=logit))
  # Use a predict function instead of coding in the values
  int <- -2.3
  b1 <- 0.56
  b2 <- 1.28
  b3 <- 0.2
  b4 <- 0
  sigmaPID <- 1
  
  X$groupErr <- rnorm(60, mean=0, sd=sigmaPID)[X$PID]
  X$groupErr2 <- rnorm(60, mean=0, sd=sigmaPID)[X$PID]
  X$MissedDose <- int +X$groupErr+b1*as.numeric(as.character(X$DrinkYN))+b2*X$ZAlcTox+b3*X$Day+b4*as.numeric(as.character(X$DrinkYN))*X$Day
  X$MissedDose <- 1/(1+exp(-X$MissedDose))
  X$MissedDose <- rbinom(2700, 1, prob=X$MissedDose)
  
  sim.data[[i]] <- X
}
