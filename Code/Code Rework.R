# #nohup R --vanilla CMD BATCH /home/gmatthews1/designedMissingness/simulationPar.R /home/gmatthews1/designedMissingness/simulationPar.Rout &
# 
# 
# nohup R --vanilla CMD BATCH '--args no.pid=60 corr.scale=1 howmuch="low" nmiss=11 nrem=2' /home/gmatthews1/designedMissingness/simulationPar.R /home/gmatthews1/designedMissingness/simulationPar60_1_low.Rout &
# nohup R --vanilla CMD BATCH '--args no.pid=120 corr.scale=1 howmuch="low" nmiss=11 nrem=2' /home/gmatthews1/designedMissingness/simulationPar.R /home/gmatthews1/designedMissingness/simulationPar120_1_low.Rout &
#   
#   nohup R --vanilla CMD BATCH '--args no.pid=60 corr.scale=1 howmuch="med" nmiss=22 nrem=4' /home/gmatthews1/designedMissingness/simulationPar.R /home/gmatthews1/designedMissingness/simulationPar60_1_med.Rout &
#   nohup R --vanilla CMD BATCH '--args no.pid=120 corr.scale=1 howmuch="med" nmiss=22 nrem=4' /home/gmatthews1/designedMissingness/simulationPar.R /home/gmatthews1/designedMissingness/simulationPar60_1_med.Rout &
#   
#   nohup R --vanilla CMD BATCH '--args no.pid=60 corr.scale=1 howmuch="high" nmiss=33 nrem=6' /home/gmatthews1/designedMissingness/simulationPar.R /home/gmatthews1/designedMissingness/simulationPar60_1_high.Rout &
#   nohup R --vanilla CMD BATCH '--args no.pid=120 corr.scale=1 howmuch="high" nmiss=33 nrem=6' /home/gmatthews1/designedMissingness/simulationPar.R /home/gmatthews1/designedMissingness/simulationPar60_1_high.Rout &

#hightime == 1
nohup R --vanilla CMD BATCH '--args no.pid=60 corr.scale=1 howmuch="low" nmiss=11 nrem=2 hightime=1' /home/gmatthews1/designedMissingness/simulationPar.R /home/gmatthews1/designedMissingness/simulationPar60_1_low_highttime.Rout &
nohup R --vanilla CMD BATCH '--args no.pid=120 corr.scale=1 howmuch="low" nmiss=11 nrem=2 hightime=1' /home/gmatthews1/designedMissingness/simulationPar.R /home/gmatthews1/designedMissingness/simulationPar120_1_low_highttime.Rout &

  nohup R --vanilla CMD BATCH '--args no.pid=60 corr.scale=1 howmuch="med" nmiss=22 nrem=4 hightime=1' /home/gmatthews1/designedMissingness/simulationPar.R /home/gmatthews1/designedMissingness/simulationPar60_1_med_highttime.Rout &
  nohup R --vanilla CMD BATCH '--args no.pid=120 corr.scale=1 howmuch="med" nmiss=22 nrem=4 hightime=1' /home/gmatthews1/designedMissingness/simulationPar.R /home/gmatthews1/designedMissingness/simulationPar120_1_med_highttime.Rout &

  nohup R --vanilla CMD BATCH '--args no.pid=60 corr.scale=1 howmuch="high" nmiss=33 nrem=6 hightime=1' /home/gmatthews1/designedMissingness/simulationPar.R /home/gmatthews1/designedMissingness/simulationPar60_1_high_highttime.Rout &
  nohup R --vanilla CMD BATCH '--args no.pid=120 corr.scale=1 howmuch="high" nmiss=33 nrem=6 hightime=1' /home/gmatthews1/designedMissingness/simulationPar.R /home/gmatthews1/designedMissingness/simulationPar120_1_high_highttime.Rout &

  
  
nohup R --vanilla CMD BATCH '--args no.pid=60 corr.scale=0.5 howmuch="low" nmiss=11 nrem=2' /home/gmatthews1/designedMissingness/simulationPar.R /home/gmatthews1/designedMissingness/simulationPar60_1_low.Rout &
nohup R --vanilla CMD BATCH '--args no.pid=120 corr.scale=0.5 howmuch="low" nmiss=11 nrem=2' /home/gmatthews1/designedMissingness/simulationPar.R /home/gmatthews1/designedMissingness/simulationPar120_1_low.Rout &

  nohup R --vanilla CMD BATCH '--args no.pid=60 corr.scale=0.5 howmuch="med" nmiss=22 nrem=4' /home/gmatthews1/designedMissingness/simulationPar.R /home/gmatthews1/designedMissingness/simulationPar60_1_med.Rout &
  nohup R --vanilla CMD BATCH '--args no.pid=120 corr.scale=0.5 howmuch="med" nmiss=22 nrem=4' /home/gmatthews1/designedMissingness/simulationPar.R /home/gmatthews1/designedMissingness/simulationPar60_1_med.Rout &

  nohup R --vanilla CMD BATCH '--args no.pid=60 corr.scale=0.5 howmuch="high" nmiss=33 nrem=6' /home/gmatthews1/designedMissingness/simulationPar.R /home/gmatthews1/designedMissingness/simulationPar60_1_high.Rout &
  nohup R --vanilla CMD BATCH '--args no.pid=120 corr.scale=0.5 howmuch="high" nmiss=33 nrem=6' /home/gmatthews1/designedMissingness/simulationPar.R /home/gmatthews1/designedMissingness/simulationPar60_1_high.Rout &

# Load necessary libraries
library(mice)
library(lme4)
library(mvtnorm)
library(lmerTest)
library(parallel)

#Passing arguments from the command line. 
#arguments are no.pid
args=(commandArgs(TRUE))

eval(parse(text=args[[1]]))
print(eval(parse(text=args[[1]])))

for(i in 1:length(args)){
  eval(parse(text=args[[i]]))
}

#Settings: 
#sampling size: 60 and 120 subjects
#Correlation settings: regular and low (1 and 0.5)
#Missingness amount: low, medium, high
#Wave missingness levels are nmiss = 11, 22, and 33.  
#Split form missingness levels are nrem = 2, 4, 6

start <- Sys.time()

no.sim = 1000 # number of simulated datasets
no.imp = 5 # number of imputed datasets

# no.pid <- 60
# corr.scale <- 1



# Install needed packages, if necessary
# install.packages(c('lme4','mvtnorm','lmerTest'))

# Set working directory
#setwd("~/GitHub/designedmissingness/")
#setwd("~/Dropbox/designedmissingnessGit/")



# Read in Daily With Tox CSV data file
#Daily <- read.csv('~/GitHub/designedmissingness/Daily with Tox.csv')
#Daily <- read.csv('Daily with Tox.csv')
Daily <- read.csv('/home/gmatthews1/designedMissingness/Daily with Tox.csv')

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
#We re-assign drinkYN to be Q7 to match the survey description in the manuscript?  
comp$Q1 <- as.factor(comp$Q1)
comp$Q2 <- as.factor(comp$Q2)
comp$Q3 <- as.factor(comp$Q3)
comp$Q4 <- as.factor(comp$Q4)
comp$Q5 <- as.factor(comp$Q5)
comp$Q6 <- as.factor(comp$Q6)
comp$Q8 <- as.integer(comp$Q7)
comp$Q7 <- as.factor(comp$DrinkYN)
comp$MissedDose <- as.factor(comp$MissedDose)
comp$Gender <- as.factor(comp$Gender)

cols.keep <- c('PID', 'Day', 'Q1', 'Q2', 'Q3', 'Q4', 'Q5', 'Q6', 'Q7', 'Q8', 'MissedDose', 'ZEduc', 'ZIncom45', 'Gender', 'Age', 'ZAlcTox', 'ZCESDFU', 'ZAUDIT')

comp <- comp[,cols.keep]



pop.mod <- glmer(MissedDose ~ Q7 + ZAlcTox + Day + (1|PID), data = comp, family=binomial(link=logit))

pop.mod.coeffs <- summary(pop.mod)$coef[,1:2]

no.var <- nrow(pop.mod.coeffs)

if (hightime == 1){
  pop.mod.coeffs[,1] <- pop.mod.coeffs[,1] * c(2,0.5,0.5,1)
}

# Correlation Matrix
corrmat <- cor(model.matrix(~.-1,data=comp[,c('ZEduc', 'ZIncom45', 'Age', 'ZAlcTox', 'ZCESDFU', 'ZAUDIT')]))
round(corrmat, 2) # display the matrix, rounding values to 2 decimal places

#Create low, medium, and high correlation settings.  
corrmat[lower.tri(corrmat)] <- corr.scale*corrmat[lower.tri(corrmat)] 
corrmat[upper.tri(corrmat)] <- corr.scale*corrmat[upper.tri(corrmat)] 

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

sim.data <- list()
for (i in 1:no.sim) {print(i)
  X <- as.data.frame(rmvnorm(no.pid, mean=rep(0,6), sigma=corrmat))
  names(X) <- c('ZEduc', 'ZIncom45', 'Age', 'ZAlcTox', 'ZCESDFU', 'ZAUDIT')
  X$PID <- 1:no.pid
  X$Gender <- as.factor(rbinom(no.pid, 1, prob=0.5))
  X <- X[rep(row.names(X), 45),] #Repeat each row 45 times, all current variables are time invariant (ZAlcTox, too??)
  X <- X[order(X$PID),] #Order by PID
  X$Day <- rep(0:44, no.pid) #Create a variable containing day variable 0-44

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
  MD <- glmer(MissedDose ~ Q7 + ZAlcTox + Day + (1|PID), data = comp, family=binomial(link=logit))
  sigmaPID <- MD@theta
  
  if (hightime == 1){
    MD@beta <- MD@beta * c(2,0.5,0.5,1)
  }
  
  X$groupErr <- rnorm(no.pid, mean=0, sd=sigmaPID)[X$PID]
  X$groupErr2 <- rnorm(no.pid, mean=0, sd=sigmaPID)[X$PID]
  X$MissedDose <- predict(MD,X,allow.new.levels=TRUE)
  X$MissedDose <- 1/(1+exp(-X$MissedDose))
  X$MissedDose <- rbinom(2700, 1, prob=X$MissedDose)
  
  sim.data[[i]] <- X
}

mod.maker.nomiss <- function(dat){ #
  modbetas <- list()
  modse <- list()
  for (i in 1:length(dat)){
    try(model <- glmer(MissedDose ~ Q7 + ZAlcTox + Day + (1|PID), data = dat[[i]], family=binomial(link=logit)))
    try(modbetas[[i]] <- summary(model)$coef[,1]) # only take the coefficients of the variables, not the intercept
    try(modse[[i]] <- summary(model)$coef[,2])
  }
  return(list(betas=modbetas,se=modse))
}

#Get results for the data with no missingess.
fullmods <- mod.maker.nomiss(sim.data)

#Remove unnecessary before analysis
rm(comp, sigma, X, b1, b2,b3,b4, fit1, fit2, fit3, fit4, fit5, fit6, fit7,
   i, int, sigmaPID)

save.image(paste0("/home/gmatthews1/designedMissingness/simResults20180509_",no.pid,"_",corr.scale,"_",howmuch,"_hightime",hightime,".RData"))

## Create MISSINGNESS
# rremove <- function(nrem, x) { # this literally just randomly removes nrem columns. Like, the entire column. Why???
#   id <- sample(length(x), nrem)
#   x[id] <- NA
#   x
# }

#grps <- list(X=c('Q1','Q8'),A=c('Q2','Q3'),B=c('Q4','Q5'),C=c('Q6','Q7'))

# split.form <- function(set,grps){
#   
#    for (i in 1:no.pid){
#      set[set$PID == i,grps$A]
#    }
#   
#   set$block[sample(1:nrow(set),nrow(set),FALSE)] <- c('A','B','C')
#   set[set$block == "A",grps$A] <- NA
#   set[set$block == "B",grps$B] <- NA
#   set[set$block == "C",grps$C] <- NA
#   set
# }

split.form <- function(set, nrem = 2, bound = 0){
  
set[set$Day >= bound, paste0("Q",1:8)] <- t(apply(set[set$Day >= bound,paste0("Q",1:8)], 1, function(x){
    ind <- sample(1:8, nrem)
    x[ind] <- NA
    as.numeric(x)
  }))
    
  return(set)

}


wave.des <- function(set, nmiss, ...) { # nmiss is the number of days each person would miss? 11?? Like, miss a fourth of the days?
  id <- replicate(no.pid, sample(2:44, size=nmiss))
  for (i in 1:no.pid) {
    set[set$PID==i & set$Day %in% id[,i], c('Q1','Q2','Q3','Q4','Q5','Q6','Q7','Q8')] <- NA
  }
  set
}

# Step 1 - Apply PM methods to simulated datasets
set.seed(917236)
split.pm <- mclapply(sim.data,function(x){split.form(x,nrem = nrem, bound = 0)},mc.cores = 12)
altered.split.pm <- mclapply(sim.data,function(x){split.form(x,nrem = nrem, bound = 2)},mc.cores = 12)
wave.pm <- mclapply(sim.data,function(x){wave.des(x,nmiss = nmiss)},mc.cores = 12)

save.image(paste0("/home/gmatthews1/designedMissingness/simResults20180509_",no.pid,"_",corr.scale,"_",howmuch,".RData"))

# Step 2 - Impute simulated datasets, for both methods

impute <- function(pm){ # pass split.data
  imp <- mclapply(pm, function(x) {mice(x, m=no.imp)},mc.cores = 18) # loop through the split.imps in the complete function
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
altered.split.data <- impute(altered.split.pm)
wave.data <- impute(wave.pm)

save.image(paste0("/home/gmatthews1/designedMissingness/simResults20180509_",no.pid,"_",corr.scale,"_",howmuch,"_hightime",hightime,".RData"))

# Step 3 - Build models from imputed simulation data for both methods
# Result: a list of lists of model coefficients for each method
# Coefficient order: Intercept, DrinkYN2, Day, ZAlcTox, and DrinkYN2:Day


#This should really be changed to run in parallel.  
#Loops are inefficient.  
mod.maker <- function(dat){ # passing split.data[[i]]
  modbetas <- list()
  modse <- list()
  for (j in 1:length(dat)){
    df <- data.frame(dat[[j]])
    try(model <- glmer(MissedDose ~ Q7 + ZAlcTox + Day + (1|PID), data = df, family=binomial(link=logit)))
    try(modbetas[[j]] <- summary(model)$coef[,1]) # only take the coefficients of the variables, not the intercept
    try(modse[[j]] <- summary(model)$coef[,2])
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
alteredsplitmods <- pm.mods(altered.split.data)
wavemods <- pm.mods(wave.data)

save.image(paste0("/home/gmatthews1/designedMissingness/simResults20180509_",no.pid,"_",corr.scale,"_",howmuch,"_hightime",hightime,".RData"))

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
altered.split.intervals <- TD(alteredsplitmods)
wave.intervals <- TD(wavemods)

coverage <- function(intervals){
  #How I would have written it.  Results are exactly the same.  
  # out <- c()
  # cov <- lapply(intervals, function(x){
  #   for (j in 1:no.var){
  #     out[j] <- (pop.mod.coeffs[j,1] >= x[j,2]) & (pop.mod.coeffs[j] <= x[j,3])
  #   }
  #   out
  # })
  # 
  # apply(do.call(rbind,cov),2,mean)
  
  fits <- list() # lists of lists of whether or not a parameter estimate fell in the sim's CI
  for (i in 1:no.sim){
    incl <- list()
    for (j in 1:no.var){
      incl[[j]] <- (pop.mod.coeffs[j] >= intervals[[i]][j,2]) & (pop.mod.coeffs[j] <= intervals[[i]][j,3])
    }
    fits[[i]] <- incl
  }
  
  fit.mat <- matrix(unlist(fits),no.var)
  cover <- apply(fit.mat,1,mean)
  return(list(fit.mat,cover))
}

split.cover <- coverage(split.intervals)
altered.split.cover <- coverage(altered.split.intervals)
wave.cover <- coverage(wave.intervals)

# CI Length #

cilength <- function(intervals){
  length <- list() # lists of lists of whether or not a parameter estimate fell in the sim's CI
  for (i in 1:no.sim){
    len <- list()
    for (j in 1:no.var){
      len[[j]] <-  intervals[[i]][j,3] - intervals[[i]][j,2]
    }
    length[[i]] <- len
  }
  
  fit.mat <- matrix(unlist(length),no.var)
  out <- apply(fit.mat,1,mean)
  return(list(fit.mat,out))
}

split.length <- cilength(split.intervals)
altered.split.length <- cilength(altered.split.intervals)
wave.length <- cilength(wave.intervals)


# BIAS #

bias <- function(intervals){
  results <- list()
  numerator <- list()
  Qdiffs <- list()
  #j-th variable
  for (j in 1:nrow(intervals[[1]])){
    Qdiffs[[j]] <- list(NA,nrow(intervals[[1]]))
    for (i in 1:no.sim){
      Qdiffs[[j]][i] <- unlist(intervals[[i]][j,1]) - unlist(pop.mod.coeffs[j,1])
    }
    
  results[j] <- mean(unlist(Qdiffs[[j]]))
    
  }
  results
}

split.bias <- bias(split.intervals)
altered.split.bias <- bias(altered.split.intervals)
wave.bias <- bias(wave.intervals)

# BIAS #

pctbias <- function(intervals){
  results <- list()
  numerator <- list()
  Qdiffs <- list()
  #j-th variable
  for (j in 1:nrow(intervals[[1]])){
    Qdiffs[[j]] <- list(NA,nrow(intervals[[1]]))
    for (i in 1:no.sim){
      Qdiffs[[j]][i] <- (unlist(intervals[[i]][j,1]) - unlist(pop.mod.coeffs[j,1])) / (unlist(pop.mod.coeffs[j,1]))
    }
    
    results[j] <- mean(unlist(Qdiffs[[j]]))
    
  }
  results
}

split.pctbias <- pctbias(split.intervals)
altered.split.pctbias <- pctbias(altered.split.intervals)
wave.pctbias <- pctbias(wave.intervals)

# MSE #

mse <- function(intervals){
  results <- c()
  numerator <- list()
  Qdiffs2 <- list()
  for (j in 1:no.var){
    Qdiffs2[[j]] <- rep(NA,no.sim)
    for (i in 1:no.sim){
      Qdiffs2[[j]][i] <- (unlist(intervals[[i]][j,1]) - unlist(pop.mod.coeffs[j,1]))^2
    }
    
  results[j] <- mean(Qdiffs2[[j]])
    
  }
  return(results)
}

split.mse <- mse(split.intervals)
altered.split.mse <- mse(altered.split.intervals)
wave.mse <- mse(wave.intervals)

# (1+(1/D))B/TD

# FMI #

#Note: Remove the rows where FMI is 0 likely due to non-vergence. 
fmi <- function(mods){
  FMI <- matrix(NA, nrow = no.sim, ncol = length(mods[[1]]$betas[[1]]))
  for (i in 1:no.sim){
    B <- apply(do.call(rbind,mods[[i]]$betas),2,var)
    W <- apply(do.call(rbind,mods[[i]]$se)^2,2,mean)
    
    # Combined variance
    TD <- (1 + 1/no.imp)*B + W
    FMI[i,] <- ((1 + (1/no.imp))*B)/TD
  }
  return(apply(FMI, 2, mean))
}

split.fmi <- fmi(splitmods)
altered.split.fmi <- fmi(alteredsplitmods)
wave.fmi <- fmi(wavemods)

# COMPARISON TABLE #
# coverage, bias, mse, fmi for each variable, each method

split.table <- data.frame(split.cover[2],cbind(split.bias,split.mse),split.fmi,split.length,split.pctbias)
colnames(split.table) <- c('Coverage','Bias','MSE','FMI','CI Length','Pct Bias')
split.table

altered.split.table <- data.frame(altered.split.cover[2],cbind(altered.split.bias,altered.split.mse),altered.split.fmi,altered.split.length,altered.split.pctbias)
colnames(altered.split.table) <- c('Coverage','Bias','MSE','FMI','CI Length','Pct Bias')
altered.split.table

wave.table <- data.frame(wave.cover[2],cbind(wave.bias,wave.mse),wave.fmi,wave.length,wave.pctbias)
colnames(wave.table) <- c('Coverage','Bias','MSE','FMI','CI Length','Pct Bias')
wave.table



save.image(paste0("/home/gmatthews1/designedMissingness/simResults20180509_",no.pid,"_",corr.scale,"_",howmuch,"_hightime",hightime,".RData"))
d <- list(split.table,altered.split.table,wave.table)
save(d, file = paste0("/home/gmatthews1/designedMissingness/tables20180509_",no.pid,"_",corr.scale,"_",howmuch,"_hightime",hightime,".RData"))


end <- Sys.time()
print(end-start)



























