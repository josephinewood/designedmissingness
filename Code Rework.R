
# Install needed packages, if necessary
# install.packages(c('lme4','mvtnorm','lmerTest'))

# Set working directory
setwd("~/GitHub/designedmissingness/")

# Load necessary libraries
library(mice)
library(lme4)
library(mvtnorm)
library(lmerTest)

# Read in Daily With Tox CSV data file
Daily <- read.csv('~/GitHub/designedmissingness/Daily with Tox.csv')

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
for (i in 1:5) {
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
  int <- MD@beta[1]
  b1 <- MD@beta[2] # DrinkYN1
  b2 <- MD@beta[3] # Day
  b3 <- MD@beta[4] # ZAlcTox
  b4 <- MD@beta[5] # DrinkYN1:Day
  sigmaPID <- MD@theta
  
  X$groupErr <- rnorm(60, mean=0, sd=sigmaPID)[X$PID]
  X$groupErr2 <- rnorm(60, mean=0, sd=sigmaPID)[X$PID]
  X$MissedDose <- predict(MD,X,allow.new.levels=TRUE)
  X$MissedDose <- 1/(1+exp(-X$MissedDose))
  X$MissedDose <- rbinom(2700, 1, prob=X$MissedDose)
  
  sim.data[[i]] <- X
}

#Remove unnecessary before analysis
rm(comp, sigma, X, b1, b2,b3,b4, fit1, fit2, fit3, fit4, fit5, fit6, fit7,
   i, int, sigmaPID)

l1 <- sim.data[[1]]

## Create MISSINGNESS
rremove <- function(nrem, x) { # this literally just randomly removes nrem columns. Like, the entire column. Why???
  id <- sample(length(x), nrem)
  x[id] <- NA
  x
}

grps <- list(X='Q1',A=c('Q2','Q3'),B=c('Q4','Q5'),C=c('Q6','Q7'))

split.form <- function(set,grps){
  set$block[sample(1:nrow(set),nrow(set),FALSE)] <- c('A','B','C')
  set[set$block == "A",grps$A] <- NA
  set[set$block == "B",grps$B] <- NA
  set[set$block == "C",grps$C] <- NA
  set
}


wave.des <- function(set, nmiss, ...) { # nmiss is the number of days each person would miss? 11?? Like, miss a fourth of the days?
  id <- replicate(60, sample(2:44, size=nmiss))
  for (i in 1:60) {
    set[set$PID==i & set$Day %in% id[,i], c('Q1','Q2','Q3','Q4','Q5','Q6','Q7')] <- NA
  }
  set
}

# alt.split <- function(set,n,bound=0){
#   set[sample(set$Day <= bound,n),c('Q1','Q2','Q3','Q4','Q5','Q6','Q7')] <- t(apply(set[sample(set$Day <= bound,n),],split.form))
#   # above won't work, is applying split.form() to different sampled rows??
#   # maybe assign og sampled rows to variable and then apply split.form() to those and then reassign to set?
#                                                                         
# }

# Step 1 - Apply PM methods to simulated datasets
set.seed(917236)
split.pm <- lapply(sim.data,function(x){split.form(x,grps)})
wave.pm <- lapply(sim.data,function(x){wave.des(x,11)})

# Step 2 - Impute simulated datasets, for both methods
split.imp <- lapply(split.pm, function(x) {mice(x, m=5)})
split.data <- lapply(split.imp,function(x) {complete(x,1)})

wave.imp <- lapply(wave.pm, function(x) {mice(x, m=5)})
wave.data <- lapply(wave.imp,function(x) {complete(x,1)})

# for(j in c(2,4,6)) { # 2, 4, 6? Why?
#   new.data <-lapply(sim.data, function(x) {split.form(x, grps)})
#   split.imp[[j/2]] <- lapply(new.data, function(x) {mice(x, m=5)})
# }

# l1.split <- split.form(l1, grps)
# l1.split.imp <- mice(l1, m=5)
# l1.split.imp2 <- mice(l1, m=50)

# imp2 <- list()
# set.seed(917236)
# for(j in c(2,4,6)) {
#   new.data <-lapply(sim.data, function(x) {split.form(x, n=j, bound=2)})
#   imp2[[j/2]] <- lapply(new.data, function(x) {mice(x, m=5)})
# } 

# wave.imp <- list()
# set.seed(918273)
# for(i in c(11, 22, 33)) {
#   new.data <- lapply(sim.data, function(x) {wave.des(set=x, nmiss=i)})
#   wave.imp[[i/11]] <- lapply(new.data, function(x) {mice(x, m=5)})
# }

# Step 3 - Build models from imputed simulation data for both methods
set.seed(12345)
splitmods <- list()
for (i in 1:length(split.data)){
  dat <- data.frame(split.data[i])
  model <- glm(MissedDose ~ DrinkYN*Day + ZAlcTox, data = dat, family=binomial(link=logit))
  # No "+ (1|PID)" because the data doesn't have PID, it indexes the rows as 1,1.1,1.2,...,2,2.1,2.2,...Change this?
  splitmods[i] <- data.frame(model$coefficients)
}

set.seed(12345)
wavemods <- list()
for (i in 1:length(wave.data)){
  dat <- data.frame(wave.data[i])
  model <- glm(MissedDose ~ DrinkYN*Day + ZAlcTox, data = dat, family=binomial(link=logit))
  # No "+ (1|PID)" because the data doesn't have PID, it indexes the rows as 1,1.1,1.2,...,2,2.1,2.2,...Change this?
  wavemods[i] <- data.frame(model$coefficients)
}

# Step 4 - Calculate mean parameter estimates for both methods
split.means
for (i in 1:5){ # there are 5 parameter estimates
  
}

analyze.comp <- function(df, ...) {
  fit <- glmer(data=df, MissedDose~DrinkYN+ZAlcTox+Day+(1|PID), family=binomial(link=logit),
               control=glmerControl(optimizer='bobyqa', optCtrl = list(maxfun=2000)))
  g0.est <- summary(fit)$coefficients[1,1]
  g0.se <- summary(fit)$coefficients[1,2]
  g0.L <- g0.est - 1.96*g0.se
  g0.U <- g0.est + 1.96*g0.se
  b1.est <- summary(fit)$coefficients[2,1]
  b1.se <- summary(fit)$coefficients[2,2]
  b1.L <- b1.est - 1.96*b1.se
  b1.U <- b1.est + 1.96*b1.se
  b2.est <- summary(fit)$coefficients[3,1]
  b2.se <- summary(fit)$coefficients[3,2]
  b2.L <- b2.est - 1.96*b2.se
  b2.U <- b2.est + 1.96*b2.se
  b3.est <- summary(fit)$coefficients[4,1]
  b3.se <- summary(fit)$coefficients[4,2]
  b3.L <- b3.est - 1.96*b3.se
  b3.U <- b3.est + 1.96*b3.se
  c(g0.est, g0.L, g0.U, b1.est, b1.L, b1.U, b2.est, b2.L, b2.U, b3.est, b3.L, b3.U)
}

library(beepr)
beep(2, expr= full.res <- lapply(sim.data, function(x) {analyze.comp(x)}))


mean(unlist(lapply(full.res, function(x) {abs(x[10]-(0.2))})))
mean(unlist(lapply(full.res, function(x) {abs(x[10]-(0.2))/0.2})))          # What is this???
mean(unlist(lapply(full.res, function(x) {(x[10]-(0.2))**2})))
mean(unlist(lapply(full.res, function(x) {(x[11]<0.2)&(x[12]>0.2)})))
mean(unlist(lapply(full.res, function(x) {x[12]-x[11]})))


mi.ml <- function(set,...) {
  fit <- with(set, glmer(MissedDose~DrinkYN+ZAlcTox+Day+(1 |PID), family=binomial(link=logit),
                         control=glmerControl(optimizer='bobyqa', optCtrl = list(maxfun=2000))))
  p <- pool(fit)
  res <- summary(p)
  c(res[1,1], res[1,6], res[1,7], res[1,9], res[2,1], res[2,6], res[2,7], res[2,9],
    res[3,1], res[3,6], res[3,7], res[3,9], res[4,1], res[4,6], res[4,7], res[4,9])
}

fit <- with(split.imp, glmer(MissedDose~DrinkYN+ZAlcTox+Day+(1 |PID), family=binomial(link=logit),
                       control=glmerControl(optimizer='bobyqa', optCtrl = list(maxfun=2000))))
fit2 <- with(imp2, glmer(MissedDose~DrinkYN+ZAlcTox+Day+(1 |PID), family=binomial(link=logit),
                         control=glmerControl(optimizer='bobyqa', optCtrl = list(maxfun=2000))))
p <- pool(fit)
p2 <- pool(fit2)


sf2.low.res <- lapply(imp2[[1]], function(x) {mi.ml(x)})
sf2.med.res <- lapply(imp2[[2]], function(x) {mi.ml(x)})
sf2.high.res <-lapply(imp2[[3]], function(x) {mi.ml(x)})
rm(imp2)

wd.low.res <- lapply(imp3[[1]], function(x) {mi.ml(x)})
wd.med.res <- lapply(imp3[[2]], function(x) {mi.ml(x)})
wd.high.res <-lapply(imp3[[3]], function(x) {mi.ml(x)})
rm(imp3)



mean(unlist(lapply(sf2.low.res, function(x) {abs(x[9]-(0.64))})))
mean(unlist(lapply(sf2.med.res, function(x) {abs(x[9]-(0.64))})))
mean(unlist(lapply(sf2.high.res, function(x) {abs(x[9]-(0.64))})))

mean(unlist(lapply(sf2.low.res, function(x) {abs(x[9]-(0.64))/0.64})))
mean(unlist(lapply(sf2.med.res, function(x) {abs(x[9]-(0.64))/0.64})))
mean(unlist(lapply(sf2.high.res, function(x) {abs(x[9]-(0.64))/0.64})))

mean(unlist(lapply(sf2.low.res, function(x) {(x[9]-(0.64))**2})))
mean(unlist(lapply(sf2.med.res, function(x) {(x[9]-(0.64))**2})))
mean(unlist(lapply(sf2.high.res, function(x) {(x[9]-(0.64))**2})))

mean(unlist(lapply(sf2.low.res, function(x) {(x[10]<0.64)&(x[11]>0.64)})))
mean(unlist(lapply(sf2.med.res, function(x) {(x[10]<0.64)&(x[11]>0.64)})))
mean(unlist(lapply(sf2.high.res, function(x) {(x[10]<0.64)&(x[11]>0.64)})))

mean(unlist(lapply(sf2.low.res, function(x) {x[11]-x[10]})))
mean(unlist(lapply(sf2.med.res, function(x) {x[11]-x[10]})))
mean(unlist(lapply(sf2.high.res, function(x) {x[11]-x[10]})))

mean(unlist(lapply(sf2.low.res, function(x) {x[12]})))
mean(unlist(lapply(sf2.med.res, function(x) {x[12]})))
mean(unlist(lapply(sf2.high.res, function(x) {x[12]})))



##SIMULATE WITH LESS INTER-SURVEY RELATIONSHIP

#Find fit to simulate
fit1 <- glm(data=comp, DrinkYN~ZIncom45+Age+ZAlcTox+ZAUDIT, family=binomial(link=logit))
fit2<- glm(data=comp, Q1~ZEduc+ZIncom45+ZAlcTox+ZCESDFU, family=binomial(link=logit))
fit3<- glm(data=comp, Q2~ZEduc+ZIncom45+Age+ZAlcTox
           +ZCESDFU+ZAUDIT+Gender, family=binomial(link=logit))
fit4<- glm(data=comp, Q3~ZEduc+ZIncom45+ZAlcTox
           +ZCESDFU+Gender, family=binomial(link=logit))
fit5<- glm(data=comp, Q4~ZEduc+ZAlcTox
           +ZCESDFU+ZAUDIT+Gender, family=binomial(link=logit))
fit6<- glm(data=comp, Q5~ZEduc+ZIncom45+ZAlcTox
           +ZAUDIT, family=binomial(link=logit))
fit7<- glm(data=comp, Q6~ZEduc+ZIncom45+Age+ZAlcTox
           +ZCESDFU+ZAUDIT, family=binomial(link=logit))

#Prep for simulate
sigma <- matrix(c(1.0, 0.3, 0.0,-0.1, 0.1,-0.2,
                  0.3, 1.0,-0.1,-0.1, 0.0, 0.1,
                  0.0,-0.1, 1.0, 0.0,-0.1, 0.0,
                  -0.1,-0.1, 0.0, 1.0, 0.1, 0.0,
                  0.1, 0.0,-0.1, 0.1, 1.0, 0.4,
                  -0.2, 0.1, 0.0, 0.0, 0.4, 1.0), nrow=6)


#Simulate
set.seed(2982364)
sim.data <- list()
for (i in 1:200) {
  X <- as.data.frame(rmvnorm(60, mean=rep(0,6), sigma=sigma))
  names(X) <- c('ZEduc', 'ZIncom45', 'Age', 'ZAlcTox', 'ZCESDFU', 'ZAUDIT')
  X$PID <- 1:60
  X$Gender <- as.factor(rbinom(60, 1, prob=0.5))
  X <- X[rep(row.names(X), 45),] #Repeat each row 45 times, all curent variables are time invarient
  X <- X[order(X$PID),] #Order by PID
  X$Day <- rep(0:44, 60) #Create a variable containing day variable 0-44
  
  
  #Create Variables
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
  X$Q7 <- rpois(2700, 2)
  X[X$DrinkYN==0,c('Q7')] <- 0
  
  #Make Missed Dose
  int <- -2.3
  b1 <- 0.56
  b2 <- 1.28
  b3 <- 0.2
  sigmaPID <- 1
  
  X$groupErr <- rnorm(60, mean=0, sd=sigmaPID)[X$PID]
  X$groupErr2 <- rnorm(60, mean=0, sd=sigmaPID)[X$PID]
  X$MissedDose <- int +b1*as.numeric(as.character(X$DrinkYN))+b2*X$ZAlcTox+b3*X$Day+X$groupErr
  X$MissedDose <- 1/(1+exp(-X$MissedDose))
  X$MissedDose <- rbinom(2700, 1, prob=X$MissedDose)
  
  sim.data[[i]] <- X
}


#Remove unnecessary before analysis
rm(comp, sigma, X, b1, b2,b3, fit1, fit2, fit3, fit4, fit5, fit6, fit7,
   i, int, sigmaPID)



full.res <- lapply(sim.data, function(x) {analyze.comp(x)})

mean(unlist(lapply(full.res, function(x) {abs(x[10]-(0.2))})))
mean(unlist(lapply(full.res, function(x) {abs(x[10]-(0.2))/0.2})))
mean(unlist(lapply(full.res, function(x) {(x[10]-(0.2))**2})))
mean(unlist(lapply(full.res, function(x) {(x[11]<0.2)&(x[12]>0.2)})))
mean(unlist(lapply(full.res, function(x) {x[12]-x[10]})))





##SIMULATE WITH NO INTER-SURVEY RELATIONSHIP


#Prep for simulate
sigma <- matrix(c(1.0, 0.3, 0.0,-0.1, 0.1,-0.2,
                  0.3, 1.0,-0.1,-0.1, 0.0, 0.1,
                  0.0,-0.1, 1.0, 0.0,-0.1, 0.0,
                  -0.1,-0.1, 0.0, 1.0, 0.1, 0.0,
                  0.1, 0.0,-0.1, 0.1, 1.0, 0.4,
                  -0.2, 0.1, 0.0, 0.0, 0.4, 1.0), nrow=6)


#Simulate
set.seed(2982364)
sim.data <- list()
for (i in 1:200) {
  X <- as.data.frame(rmvnorm(60, mean=rep(0,6), sigma=sigma))
  names(X) <- c('ZEduc', 'ZIncom45', 'Age', 'ZAlcTox', 'ZCESDFU', 'ZAUDIT')
  X$PID <- 1:60
  X$Gender <- as.factor(rbinom(60, 1, prob=0.5))
  X <- X[rep(row.names(X), 45),] #Repeat each row 45 times, all curent variables are time invarient
  X <- X[order(X$PID),] #Order by PID
  X$Day <- rep(0:44, 60) #Create a variable containing day variable 0-44
  
  
  #Create Variables
  X$DrinkYN <- as.factor(rbinom(2700, 1, prob = 0.22))
  X$Q1 <- as.factor(rbinom(2700, 1, prob = 0.17))
  X$Q2 <- as.factor(rbinom(2700,1,prob=0.15))
  X$Q3 <- as.factor(rbinom(2700,1,prob=0.24))
  X$Q4 <- as.factor(rbinom(2700,1,prob=0.22))
  X$Q5 <- as.factor(rbinom(2700,1,prob=0.12))
  X$Q6 <- as.factor(rbinom(2700,1,prob=0.18))
  X$Q7 <- rpois(2700, 2)
  X[X$DrinkYN==0,c('Q7')] <- 0
  
  #Make Missed Dose
  int <- -2.3
  b1 <- 0.56
  b2 <- 1.28
  b3 <- 0.2
  sigmaPID <- 1
  
  X$groupErr <- rnorm(60, mean=0, sd=sigmaPID)[X$PID]
  X$groupErr2 <- rnorm(60, mean=0, sd=sigmaPID)[X$PID]
  X$MissedDose <- int +b1*as.numeric(as.character(X$DrinkYN))+b2*X$ZAlcTox+b3*X$Day+X$groupErr
  X$MissedDose <- 1/(1+exp(-X$MissedDose))
  X$MissedDose <- rbinom(2700, 1, prob=X$MissedDose)
  
  sim.data[[i]] <- X
}


#Remove unnecessary before analysis
rm(comp, sigma, X, b1, b2,b3
   i, int, sigmaPID)



full.res <- lapply(sim.data, function(x) {analyze.comp(x)})

mean(unlist(lapply(full.res, function(x) {abs(x[10]-(0.2))})))
mean(unlist(lapply(full.res, function(x) {abs(x[10]-(0.2))/0.2})))
mean(unlist(lapply(full.res, function(x) {(x[10]-(0.2))**2})))
mean(unlist(lapply(full.res, function(x) {(x[11]<0.2)&(x[12]>0.2)})))
mean(unlist(lapply(full.res, function(x) {x[12]-x[10]})))