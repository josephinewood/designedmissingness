
# Install needed packages, if necessary
# install.packages(c('lme4','mvtnorm','lmerTest'))

# Set working directory
setwd("~/GitHub/designedmissingness/")
setwd("~/Dropbox/designedmissingnessGit/")

# Load necessary libraries
library(mice)
library(lme4)
library(mvtnorm)
library(lmerTest)

# Read in Daily With Tox CSV data file
Daily <- read.csv('~/GitHub/designedmissingness/Daily with Tox.csv')
#Daily <- read.csv('Daily with Tox.csv')

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
comp$Q7 <- as.factor(comp$DrinkYN)
comp$Q8 <- as.integer(comp$Q7)
comp$MissedDose <- as.factor(comp$MissedDose)
comp$Gender <- as.factor(comp$Gender)

cols.keep <- c('PID', 'Day', 'Q1', 'Q2', 'Q3', 'Q4', 'Q5', 'Q6', 'Q7', 'Q8', 'MissedDose', 'ZEduc', 'ZIncom45', 'Gender', 'Age', 'ZAlcTox', 'ZCESDFU', 'ZAUDIT')

comp <- comp[,cols.keep]

no.pid <- 60

pop.mod <- glmer(MissedDose ~ Q7 + ZAlcTox + Day + Q7*Day + (1|PID), data = comp, family=binomial(link=logit))

pop.mod.coeffs <- summary(pop.mod)$coef[,1:2]

# Correlation Matrix
corrmat <- cor(model.matrix(~.-1,data=comp[,c('ZEduc', 'ZIncom45', 'Age', 'ZAlcTox', 'ZCESDFU', 'ZAUDIT')]))
round(corrmat, 2) # display the matrix, rounding values to 2 decimal places

# OG fit models for the questions using variables simulated via corrmat
fit1 <- glm(data=comp, Q7~ZIncom45+Age+ZAlcTox+ZAUDIT, family=binomial(link=logit))
fit2<- glm(data=comp, Q1~ZEduc+ZIncom45+ZAlcTox+ZCESDFU, family=binomial(link=logit))
fit3<- glm(data=comp, Q2~Q1+ZEduc+ZIncom45+Age+ZAlcTox+ZCESDFU+ZAUDIT+Gender, family=binomial(link=logit))
fit4<- glm(data=comp, Q3~Q1+Q2+ZEduc+ZIncom45+ZAlcTox+ZCESDFU+Gender, family=binomial(link=logit))
fit5<- glm(data=comp, Q4~Q1+Q2+Q3+ZEduc+ZAlcTox+ZCESDFU+ZAUDIT+Gender, family=binomial(link=logit))
fit6<- glm(data=comp, Q5~Q1+Q3+Q4+ZEduc+ZIncom45+ZAlcTox+ZAUDIT, family=binomial(link=logit))
fit7<- glm(data=comp, Q6~Q3+Q5+Q7+ZEduc+ZIncom45+Age+ZAlcTox+ZCESDFU+ZAUDIT, family=binomial(link=logit))

# Simulate
set.seed(1234)
no.sim = 5 # number of simulated datasets
sim.data <- list()
for (i in 1:no.sim) {
  X <- as.data.frame(rmvnorm(60, mean=rep(0,6), sigma=corrmat))
  names(X) <- c('ZEduc', 'ZIncom45', 'Age', 'ZAlcTox', 'ZCESDFU', 'ZAUDIT')
  X$PID <- 1:60
  X$Gender <- as.factor(rbinom(60, 1, prob=0.5))
  X <- X[rep(row.names(X), 45),] #Repeat each row 45 times, all current variables are time invariant (ZAlcTox, too??)
  X <- X[order(X$PID),] #Order by PID
  X$Day <- rep(0:44, 60) #Create a variable containing day variable 0-44

  # Create Variables
  X$Q7 <- 1/(1+exp(-1*predict(fit1, newdata = X)))
  X$Q7 <- as.factor(rbinom(2700, 1, prob = X$Q7))
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
  X$Q8 <- rpois(2700, 2) # lambda = 2...mean of 2. why???
  X[X$Q7==0,c('Q8')] <- 0
  
  # Make Missed Dose
  MD <- glmer(MissedDose ~ Q7 + ZAlcTox + Day + Q7*Day + (1|PID), data = comp, family=binomial(link=logit))
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

## Create MISSINGNESS
# rremove <- function(nrem, x) { # this literally just randomly removes nrem columns. Like, the entire column. Why???
#   id <- sample(length(x), nrem)
#   x[id] <- NA
#   x
# }

grps <- list(X=c('Q1','Q8'),A=c('Q2','Q3'),B=c('Q4','Q5'),C=c('Q6','Q7'))

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
    set[set$PID==i & set$Day %in% id[,i], c('Q1','Q2','Q3','Q4','Q5','Q6','Q7','Q8')] <- NA
  }
  set
}

# Step 1 - Apply PM methods to simulated datasets
set.seed(917236)
split.pm <- lapply(sim.data,function(x){split.form(x,grps)})
wave.pm <- lapply(sim.data,function(x){wave.des(x,11)})

no.imp = 5 # number of imputed datasets

# Step 2 - Impute simulated datasets, for both methods

impute <- function(pm){ # pass split.data
  imp <- lapply(pm, function(x) {mice(x, m=no.imp)}) # loop through the split.imps in the complete function
  res.data <- list()
  for (i in 1:no.sim){
    res.data[[i]] <- list()
    for (j in 1:no.imp){
      res.data[[i]][[j]] <- complete(imp[[i]],j)
    }
  }
  res.data
}

split.data <- impute(split.pm)
wave.data <- impute(wave.pm)

# Step 3 - Build models from imputed simulation data for both methods
# Result: a list of lists of model coefficients for each method
# Coefficient order: Intercept, DrinkYN2, Day, ZAlcTox, and DrinkYN2:Day

mod.maker <- function(dat){ # passing split.data[[i]]
  modbetas <- list()
  modse <- list()
  for (j in 1:length(dat)){
    df <- data.frame(dat[[j]])
    model <- glmer(MissedDose ~ Q7 + ZAlcTox + Day + Q7*Day + (1|PID), data = df, family=binomial(link=logit))
    modbetas[[j]] <- summary(model)$coef[,1] # only take the coefficients of the variables, not the intercept
    modse[[j]] <- summary(model)$coef[,2]
  }
  return(list(betas=modbetas,se=modse))
}

pm.mods <- function(dat){
  set.seed(12345)
  resmods <- list()
  for (i in 1:no.sim){
    resmods[i] <- list(mod.maker(dat[[i]]))
  }
  resmods
}

splitmods <- pm.mods(split.data)
wavemods <- pm.mods(wave.data)



TD <- function(mods){
  intervals <- list()
  for (i in 1:no.sim){
    Q <- apply(do.call(rbind,mods[[i]]$betas),2,mean)
    B <- apply(do.call(rbind,mods[[i]]$betas),2,var)
    W <- apply(do.call(rbind,mods[[i]]$se)^2,2,mean)
    
    #Combioned standard error
    Tm <- sqrt((1+1/no.imp)*B + W)
    df <- (no.imp - 1) * (1 + 1/no.imp*W/B)^2
    
    intervals[[i]] <- cbind(Q ,lower = Q - qt(0.975,df)*Tm, upper = Q + qt(0.975,df)*Tm)
  }
  intervals
}

split.intervals <- TD(splitmods)
wave.intervals <- TD(wavemods)

coverage <- function(intervals){
  fits <- list() # lists of lists of whether or not a parameter estimate fell in the sim's CI
  for (i in 1:no.sim){
    incl <- list()
    for (j in 1:5){
      incl[[j]] <- (pop.mod.coeffs[j] >= intervals[[i]][j,2]) & (pop.mod.coeffs[j] <= intervals[[i]][j,3])
    }
    fits[[i]] <- incl
  }
  
  fit.mat <- matrix(unlist(allsplit.fits),5)
  cover <- apply(splitmat,2,mean)
  return(list(fit.mat,cover))
}

split.cover <- coverage(split.intervals)
wave.cover <- coverage(wave.intervals)


# BIAS #

bias <- function(intervals){
  results <- list()
  numerator <- list()
  Qdiffs <- list()
  for (j in 1:no.imp){
    Qdiffs[[j]] <- list(NA,5)
    for (i in 1:no.sim){
      Qdiffs[[j]][i] <- unlist(intervals[[i]][j,1]) - unlist(pop.mod.coeffs[j,1])
    }
    numerator[j] <- Reduce("+",Qdiffs[[j]])
    results[j] <- unlist(numerator[j])/no.sim
  }
  results
}

split.bias <- bias(split.intervals)
wave.bias <- bias(wave.intervals)

# MSE #

mse <- function(intervals){
  results <- list()
  numerator <- list()
  Qdiffs2 <- list()
  for (j in 1:no.imp){
    Qdiffs2[[j]] <- list(NA,5)
    for (i in 1:no.sim){
      Qdiffs2[[j]][i] <- (unlist(intervals[[i]][j,1]) - unlist(pop.mod.coeffs[j,1]))^2
    }
    numerator[j] <- Reduce("+",Qdiffs2[[j]])
    results[j] <- unlist(numerator[j])/no.sim
  }
  results
}

split.mse <- mse(split.intervals)
wave.mse <- mse(wave.intervals)

# (1+(1/D))B/TD

# FMI #

fmi <- function(mods){
  for (i in 1:no.sim){
    B <- apply(do.call(rbind,mods[[i]]$betas),2,var)
    W <- apply(do.call(rbind,mods[[i]]$se)^2,2,mean)
    
    # Combined variance
    TD <- (1 + 1/no.imp)*B + W
    FMI <- ((1 + (1/no.imp))*B)/TD
  }
  FMI
}

split.fmi <- fmi(splitmods)
wave.fmi <- fmi(wavemods)

# COMPARISON TABLE #
# coverage, bias, mse, fmi for each variable, each method

split.table <- data.frame(split.cover[2],cbind(split.bias,split.mse),split.fmi)
colnames(split.table) <- c('Coverage','Bias','MSE','FMI')
split.table

wave.table <- data.frame(wave.cover[2],cbind(wave.bias,wave.mse),wave.fmi)
colnames(wave.table) <- c('Coverage','Bias','MSE','FMI')
wave.table





























