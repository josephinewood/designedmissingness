
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
comp$Q7 <- as.factor(comp$DrinkYN)
comp$Q8 <- as.integer(comp$Q7)
comp$MissedDose <- as.factor(comp$MissedDose)
comp$Gender <- as.factor(comp$Gender)

cols.keep <- c('PID', 'Day', 'Q1', 'Q2', 'Q3', 'Q4', 'Q5', 'Q6', 'Q7', 'Q8', 'MissedDose', 'ZEduc', 'ZIncom45', 'Gender', 'Age', 'ZAlcTox', 'ZCESDFU', 'ZAUDIT')

comp <- comp[,cols.keep]

no.pid <- 60

pop.mod <- glm(data=comp,MissedDose~Q7+ZAlcTox+Day+Q7*Day,family=binomial(link=logit))

# pop.pids <- list()
pop.mod.coeffs <- coef(pop.mod)[2:5] # if glmer object: take the first row of coefficients, excluding the intercept, the betas are the same for each PID but the intercept varies.


# sim.pids <- list()
# sim.pid.means <- list()
# 
# for (k in 1:no.pid){  # only 59 PIDs instead of 60...is it using PID = 1 as a baseline and that intercept is 0? the other PID list have 60...oh is that because we've simulated another person? I think this is the cause. Is that an issue? Maybe just don't compare results to the original population, just to their simulated parent data set?
#     pop.pids[[k]]<- coef(pop.mod)[[1]][[1]][k]
# }

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
  # This is where issues come in!!!
  # Closest I could get to these numbers was:
  # MD <- glm(data=comp,MissedDose~ZEduc+ZIncom45+Age+ZAlcTox+ZCESDFU+ZAUDIT,family=binomial(link=logit))
  # MD <- glm(data=comp,MissedDose~Q7+ZAlcTox+Day+Q7*Day,family=binomial(link=logit))
  MD <- glmer(MissedDose ~ Q7 + ZAlcTox + Day + Q7*Day + (1|PID), data = comp, family=binomial(link=logit))
  # Use a predict function instead of coding in the values
  # int <- MD@beta[1]
  # b1 <- MD@beta[2] # DrinkYN1
  # b2 <- MD@beta[3] # Day
  # b3 <- MD@beta[4] # ZAlcTox
  # b4 <- MD@beta[5] # DrinkYN1:Day
  sigmaPID <- MD@theta
  
  X$groupErr <- rnorm(60, mean=0, sd=sigmaPID)[X$PID]
  X$groupErr2 <- rnorm(60, mean=0, sd=sigmaPID)[X$PID]
  X$MissedDose <- predict(MD,X,allow.new.levels=TRUE)
  X$MissedDose <- 1/(1+exp(-X$MissedDose))
  X$MissedDose <- rbinom(2700, 1, prob=X$MissedDose)
  
  sim.data[[i]] <- X
}

sim.mods <- list()
sim.mod.intercepts <- list()
sim.mod.coeffs <- list()

for (i in 1:no.sim){
  sim.mods[[i]] <- glmer(MissedDose ~ Q7 + ZAlcTox + Day + Q7*Day + (1|PID), data = as.data.frame(sim.data[i]), family=binomial(link=logit))
  #sim.mod.intercepts[[i]] <- coef(sim.mods[[i]])[[1]][1]
  sim.mod.coeffs[[i]] <- coef(sim.mods[[i]])[[1]][1,2:5] # take the first row and only the last 4 columns to get coefficients for everything but the intercepts
}

# sim.pids <- list()
# sim.pid.means <- list()
# 
# for (k in 1:no.pid){
#   for (i in 1:no.sim){
#       sim.pids[[k]][[i]] <- sim.mod.intercepts[[i]][[1]][[1]][k] # subscripting thing and no. of items not multiple of replacement length warning
#   }
# }
# 
# for (k in 1:no.pid){
#   sim.pid.means[k] <- mean(sim.pids[[k]])   # means for wave and split are very similar but not the same!
# }

#Remove unnecessary before analysis
rm(comp, sigma, X, b1, b2,b3,b4, fit1, fit2, fit3, fit4, fit5, fit6, fit7,
   i, int, sigmaPID)

l1 <- sim.data[[1]]

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

no.imp = 5 # number of imputed datasets

# Step 2 - Impute simulated datasets, for both methods
split.imp <- lapply(split.pm, function(x) {mice(x, m=no.imp)}) # loop through the split.imps in the complete function
split.data <- list()
for (i in 1:no.sim){
  split.data[[i]] <- list()
  for (j in 1:no.imp){
  split.data[[i]][[j]] <- complete(split.imp[[i]],j)
  }
}

# s1 <- data.frame(split.data[1][1])
# s3 <- data.frame(split.data[3][1]) # very similar, but I guess that is obvious as the missingness isn't a lot
# w1 <- data.frame(wave.data[1][1])


wave.imp <- lapply(wave.pm, function(x) {mice(x, m=no.imp)})
wave.data <- list()
for (i in 1:no.sim){
  wave.data[[i]] <- list()
  for (j in 1:no.imp){
    wave.data[[i]][[j]] <- complete(wave.imp[[i]],j)
  }
}

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
# Result: a list of lists of model coefficients for each method
# Coefficient order: Intercept, DrinkYN2, Day, ZAlcTox, and DrinkYN2:Day

mod.maker <- function(dat){ # passing split.data[[i]]
  mods <- list()
  for (j in 1:length(dat)){
    df <- data.frame(dat[[j]])
    model <- glmer(MissedDose ~ Q7 + ZAlcTox + Day + Q7*Day + (1|PID), data = df, family=binomial(link=logit))
    mods[[j]] <- coef(model)[[1]][1,2:5] # only take the coefficients of the variables, not the intercept
  }
  mods
}
# df <- data.frame(split.data[[1]][[1]])
# model <- glmer(MissedDose ~ DrinkYN*Day + Q1 + Q2 + Q3 + Q4 + Q5 + Q6 + Q7 + (1|PID), data = df, family=binomial(link="logit"))
# mod <- mod.maker(split.data[[1]]) # <-- every person has their own intercept??? the coefficients are the same except ints

set.seed(12345)
splitmods <- list()
for (i in 1:no.sim){
  splitmods[i] <- list(mod.maker(split.data[[i]])) # fails to converge
}

set.seed(12345)
wavemods <- list()
for (i in 1:no.sim){
  wavemods[i] <- list(mod.maker(wave.data[[i]])) # also fails to converge
}

# Step 4 - Calculate mean parameter estimates for both methods
# compare means of coefficents overall and then find means of each person's intercept and compare those to original???
split.coeff.means <- list()
wave.coeff.means <- list()

for (i in 1:no.sim){
  for (j in 1:no.imp){
    split.coeff.means[[i]] <- list(colMeans(as.data.frame(do.call(rbind, splitmods[[i]])))) # using split.coeff.means[[i]] throws an error about the subscript being out of bounds, but if run with split.coeff.means[i] and THEN double brackets, it works???
  }
}

for (i in 1:no.sim){
  for (j in 1:no.imp){
    wave.coeff.means[[i]] <- list(colMeans(as.data.frame(do.call(rbind, wavemods[[i]])))) # same here
  }
}

# split.pids <- list()
# wave.pids <- list()
# 
# split.pid.means <- list()
# wave.pid.means <- list()
# 
# 
# for (k in 1:no.pid){
#   for (i in 1:no.sim){
#     for (j in 1:no.imp){
#       split.pids[[k]][i][j] <- splitmods[[i]][[j]][[1]][[1]][k] # same subscripting thing and no. of items not multiple of replacement length warning...results in 5 intercepts per PID, not 25 as expected, I think it takes the first intercept for each set of 5 imputation sets
#     }
#   }
# }
# 
# for (k in 1:no.pid){
#   for (i in 1:no.sim){
#     for (j in 1:no.imp){
#       wave.pids[[k]][i][j] <- wavemods[[i]][[j]][[1]][[1]][k] # same subscripting thing and no. of items not multiple of replacement length warning...results in 5 intercepts per PID, not 25 as expected
#     }
#   }
# }
# 
# for (k in 1:no.pid){
#   split.pid.means[k] <- mean(split.pids[[k]])   # means for wave and split are very similar but not the same!
# }
# 
# for (k in 1:no.pid){
#   wave.pid.means[k] <- mean(wave.pids[[k]])
# }
# 
# 
# # Compare means to sim datasets and original true model??? 
# 
# 
# # PIDs: Is difference significant??? Do a permutation test?!
# # Coefficients: is the imputed missing estimate within the CI of the simulated/true estimate? That's what OG did
# 
# pid.means <- cbind(sim.pid.means,split.pid.means,wave.pid.means)
# # coeff.means <- cbind(sim.mod.coeffs,split.coeff.means,wave.coeff.means)
# 
# 
# sim.split.pid.means <- list()
# 
# for (k in 1:no.pid){
#   sim.split.pid.means[k] <- sim.pid.means[[k]]-split.pid.means[[k]]
# }

# Std. Error = summary(sim.mods[[1]])$coefficients[i,2] where i is the row of the variable you care about
sim.allstd.errors <- list() # USE IF CARE ABOUT ALL VARIABELS, NOT JUST Q7
for (i in 1:no.sim){
   for (j in 2:5){ # number of variables is 5
     sim.allstd.errors[[i]][j-1] <- summary(sim.mods[[i]])$coefficients[j,2] # subscripting error again
   }
 }

sim.allcis <- list()    # FIX IF NEED CIs FOR ALL VARIABLES, NOT JUST Q7
for (i in 1:no.sim){
   allci <- list()
   for (j in 1:4){
     Est <- sim.mod.coeffs[[i]][j]
     L <- c(sim.mod.coeffs[[i]][j] - sim.allstd.errors[[i]][j])
     U <- c(sim.mod.coeffs[[i]][j] + sim.allstd.errors[[i]][j])
     allci[[j]] <- cbind(Est,L,U)
   }
   sim.allcis[[i]] <- allci
}

# USE IF CARE ABOUT ALL OF THE VARIABLES
allsplit.fits <- list() # lists of lists of whether or not a parameter estimate fell in the sim's CI
for (i in 1:no.sim){
  incl <- list()
  for (j in 1:4){
    incl[[j]] <- (split.coeff.means[[i]][[1]][j] >= sim.allcis[[i]][[j]][[2]]) & (split.coeff.means[[i]][[1]][j] <= sim.allcis[[i]][[j]][[3]])
  }
  allsplit.fits[[i]] <- incl
}

# USE IF ONLY CARE ABOUT Q7 

sim.std.errors <- list() # USE IF ONLY CARE ABOUT Q7
for (i in 1:no.sim){
  sim.std.errors[i] <- summary(sim.mods[[i]])$coefficients[2,2] # subscripting error again
}

sim.cis <- list()
for (i in 1:no.sim){
  Est <- sim.mod.coeffs[[i]][[1]]
  L <- c(sim.mod.coeffs[[i]][[1]] - 1.96*sim.std.errors[[i]])
  U <- c(sim.mod.coeffs[[i]][[1]] + 1.96*sim.std.errors[[i]])
  sim.cis[[i]] <- cbind(Est,L,U)
}

split.fits <- list() # lists of lists of whether or not a parameter estimate fell in the sim's CI
for (i in 1:no.sim){
  split.fits[[i]] <- (split.coeff.means[[i]][[1]][1] >= sim.cis[[i]][2]) & (split.coeff.means[[i]][[1]][1] <= sim.cis[[i]][3])
}

wave.fits <- list() # lists of lists of whether or not a parameter estimate fell in the sim's CI
for (i in 1:no.sim){
  wave.fits[[i]] <- (wave.coeff.means[[i]][[1]][1] >= sim.cis[[i]][2]) & (wave.coeff.means[[i]][[1]][1] <= sim.cis[[i]][3])
}


analyze.comp <- function(df, ...) {
  fit <- glm(data=df, MissedDose~Q7+ZAlcTox+Day+Q7*Day, family=binomial(link=logit))
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

# library(beepr)
# beep(2, expr= 

full.res <- lapply(sim.data, function(x) {analyze.comp(x)})


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