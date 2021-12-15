# clear workspace
rm(list=ls())


## load library ##

#library(jagsUI)  #once loaded, XQuartz will be used to view data sets, so wait until ready to run model

library(lme4)
library(car)
library(plyr)
library(dplyr)
library(data.table) ## for fread
library(ggplot2)
library(gridExtra)
library(ggpubr)

############# set working directory ###########

#setwd("~/Desktop/Penn State/Classwork/CEC_Project/Chapter 2/R Scripts/Concentration Analyses/Concentration Analyses - Part 1")

####################################################################


#########################################################################  
####################### OVARY TISSUE ####################################
#########################################################################


###### ###### ###### ######
## load contaminant data ###
###### ###### ###### ###### 

dat1 <- read.csv('All_Compounds_Cleaned.csv')

## rename Site ID to Site and  Tissue Type to Tissue
names(dat1)[4]<-"Site"
names(dat1)[5]<-"Tissue"

## subset out compounds of interest and their corresponding detection limits (DLs), and Site and Tissue

dat <- dat1 %>%
  select('Site',
         'Tissue',
         'Total_PCB',
         'Total_PCB_DL',
         'Hexachlorobenzene','Hexachlorobenzene_DL',
         'p.p..DDE','p.p..DDE_DL',
         'Galaxolide','Galaxolide_DL',
         'Naphthalene','Naphthalene_DL',
         'Phenanthrene','Phenanthrene_DL',
         'Methyl.Triclosan','Methyl_Triclosan_DL',
         'Mirex','Mirex_DL',
         'Pyrene','Pyrene_DL')


###### ###### ###### ######  
## load land use data #### 
###### ###### ###### ######  

land_use <- fread('Upstream_Immediate_totals.csv')
unique(land_use$Site_ID)

## order by site
land_use <- land_use[order(Site_ID)]

## subset orginal data for columns of interest
land_use_2 <- land_use[,c('Site_ID','Scale','Developed_total', 'BMP_Intensity', 'Cultivated_total')]

## create columns in contaminant data set for each landuse type
dat$PCT.AG <- 1:186
dat$PCT.DEV <- 1:186
dat$BMP <- 1:186

###### ###### 
#### AG ###### 
###### ###### 

## fill in columns (in contaminant data set) with appropriate agricultural landuse value for each site

dat$PCT.AG <- ifelse(dat$Site=='ANT', as.numeric((land_use_2[2,5])),dat$PCT.AG) #ANT
dat$PCT.AG <- ifelse(dat$Site=='BE', as.numeric((land_use_2[4,5])),dat$PCT.AG) #BE
dat$PCT.AG <- ifelse(dat$Site=='CHIL', as.numeric((land_use_2[6,5])),dat$PCT.AG) #CHIL
dat$PCT.AG <- ifelse(dat$Site=='LOY', as.numeric((land_use_2[8,5])),dat$PCT.AG) #LOY
dat$PCT.AG <- ifelse(dat$Site=='MON', as.numeric((land_use_2[10,5])),dat$PCT.AG) #MON
dat$PCT.AG <- ifelse(dat$Site=='PC', as.numeric((land_use_2[12,5])),dat$PCT.AG) #PC
dat$PCT.AG <- ifelse(dat$Site=='POT', as.numeric((land_use_2[14,5])),dat$PCT.AG) #POT
dat$PCT.AG <- ifelse(dat$Site=='SUSHA', as.numeric((land_use_2[16,5])),dat$PCT.AG) #SUSHA
dat$PCT.AG <- ifelse(dat$Site=='WBM', as.numeric((land_use_2[18,5])),dat$PCT.AG) #WMB
dat$PCT.AG <- ifelse(dat$Site=='WBMC', as.numeric((land_use_2[20,5])),dat$PCT.AG) #WBMC
dat$PCT.AG <- ifelse(dat$Site=='WYA', as.numeric((land_use_2[22,5])),dat$PCT.AG) #WYA

## change values to numeric
dat$PCT.AG <- as.numeric(dat$PCT.AG)

## standardize and logit-transform 
dat$ag.logit.z <- as.numeric(scale(logit(dat$PCT.AG)))
#hist(dat$ag.logit.z)

###### ###### 
#### DEV ######  
###### ###### 

## fill in columns (in contaminant data set) with appropriate developed landuse value for each site

dat$PCT.DEV <- ifelse(dat$Site=='ANT', as.numeric((land_use_2[2,3])),dat$PCT.DEV) #ANT
dat$PCT.DEV <- ifelse(dat$Site=='BE', as.numeric((land_use_2[4,3])),dat$PCT.DEV) #BE
dat$PCT.DEV <- ifelse(dat$Site=='CHIL', as.numeric((land_use_2[6,3])),dat$PCT.DEV) #CHIL
dat$PCT.DEV <- ifelse(dat$Site=='LOY', as.numeric((land_use_2[8,3])),dat$PCT.DEV) #LOY
dat$PCT.DEV <- ifelse(dat$Site=='MON', as.numeric((land_use_2[10,3])),dat$PCT.DEV) #MON
dat$PCT.DEV <- ifelse(dat$Site=='PC', as.numeric((land_use_2[12,3])),dat$PCT.DEV) #PC
dat$PCT.DEV <- ifelse(dat$Site=='POT', as.numeric((land_use_2[14,3])),dat$PCT.DEV) #POT
dat$PCT.DEV <- ifelse(dat$Site=='SUSHA', as.numeric((land_use_2[16,3])),dat$PCT.DEV) #SUSHA
dat$PCT.DEV <- ifelse(dat$Site=='WBM', as.numeric((land_use_2[18,3])),dat$PCT.DEV) #WMB
dat$PCT.DEV <- ifelse(dat$Site=='WBMC', as.numeric((land_use_2[20,3])),dat$PCT.DEV) #WBMC
dat$PCT.DEV <- ifelse(dat$Site=='WYA', as.numeric((land_use_2[22,3])),dat$PCT.DEV) #WYA

## change values to numeric
dat$PCT.DEV <- as.numeric(dat$PCT.DEV)

## standardize and logit-transform 
dat$dev.logit.z <- as.numeric(scale(logit(dat$PCT.DEV)))
#hist(dat$dev.logit.z)

###### ######
#### BMPs ### 
###### ######

## fill in columns (in contaminant data set) with appropriate BMP intensity value for each site

dat$BMP <- ifelse(dat$Site=='ANT', as.numeric((land_use_2[1,4])),dat$BMP) #ANT
dat$BMP <- ifelse(dat$Site=='BE', as.numeric((land_use_2[3,4])),dat$BMP) #BE
dat$BMP <- ifelse(dat$Site=='CHIL', as.numeric((land_use_2[5,4])),dat$BMP) #CHIL
dat$BMP <- ifelse(dat$Site=='LOY', as.numeric((land_use_2[7,4])),dat$BMP) #LOY
dat$BMP <- ifelse(dat$Site=='MON', as.numeric((land_use_2[9,4])),dat$BMP) #MON
dat$BMP <- ifelse(dat$Site=='PC', as.numeric((land_use_2[11,4])),dat$BMP) #PC
dat$BMP <- ifelse(dat$Site=='POT', as.numeric((land_use_2[13,4])),dat$BMP) #POT
dat$BMP <- ifelse(dat$Site=='SUSHA', as.numeric((land_use_2[15,4])),dat$BMP) #SUSHA
dat$BMP <- ifelse(dat$Site=='WBM', as.numeric((land_use_2[17,4])),dat$BMP) #WMB
dat$BMP <- ifelse(dat$Site=='WBMC', as.numeric((land_use_2[19,4])),dat$BMP) #WBMC
dat$BMP <- ifelse(dat$Site=='WYA', as.numeric((land_use_2[21,4])),dat$BMP) #WYA

## change values to numeric
dat$BMP <- as.numeric(dat$BMP)

## standardize and logit-transform

dat$BMP.z <- as.numeric(scale(dat$BMP))
#hist(dat$BMP.logit.z)

#### subset for tissue of interest ###### 

dat <- subset(dat, Tissue =="Ovary")

###### ###### ###### ###### ###### ###### ###### ###### ### 
###### ###### FOR-LOOP FOR OVARY TISSUE MODEL ###### ###### 
###### ###### ###### ###### ###### ###### ###### ###### ### 

## save number of columns for future use in model
dat_length <- ncol(dat)

## save column names for future use in model
dat_col_names <- colnames(dat)

## create empty matrix to fill in with probility of a positive effect and estimated effects
Immed_Ovary_Prob_Pos <- matrix(nrow =(dat_length-8), ncol = 3)
Immed_Ovary_effect <- matrix(nrow =(dat_length-8), ncol = 3)

for (colIdx in seq(from=3, to=dat_length-7, by=2)) { ### BEGIN FOR-LOOP (looping through contaminants)
  
  # print to keep track of model progress
  print("NEW LOOP ITERATION.  WORKING ON:")
  print(dat_col_names[colIdx])
  
  # Set the appropriate Censoring Limit for contaminant
  dat$censorLimitVec <- dat[,colIdx+1]
  
  # Must tell JAGS which observations are ABOVE Censoring limit
  dat$isNotCensored <- (dat[,colIdx] > dat$censorLimitVec) # Returns TRUE if not below DL and NA if below DL
  
  # Convert NA in dat$isNotCensored to FALSE
  dat$isNotCensored <- ifelse(is.na(dat$isNotCensored), FALSE, dat$isNotCensored)
  
  #################################################################
  ########## BUGS CODE ############################################
  #################################################################
  
  ## copy of data set
  jags.dat.2 <- copy(dat)
  
  ## load jagsUI package ##
  library(jagsUI)
  
  #################################################
  ###################  AG LANDUSE  ################
  #################################################
  
  ## print info to know which loop model is on
  print('CULTIVATED IMMEDIATE OVARY')
  print(dat_col_names[colIdx])
  
  # Define the model in the BUGS language and write a text file
  sink("vary.int.fixed.slope.txt")
  cat("
    model {
    
    
    # Likelihood: 
    # Level-1 of the model
    for (i in 1:n){ 
    isNotCensored[i] ~ dinterval(y[i] , censorLimitVec[i] )
    y[i] ~ dnorm(mu[i], tau)               
    mu[i] <- alpha[group[i]]   
    } 
    
    
    # Level-2 of the model
    for(j in 1:J){
    alpha[j] ~ dnorm(mu.alpha[j],tau.alpha)
    mu.alpha[j] <- alpha.gamma + beta.gamma *ag[j]
    }
    
    
    # Priors and derived quantities
    sigma ~ dunif(0, 100)
    alpha.gamma ~ dnorm(0, 0.0001)
    beta.gamma ~ dnorm(0,0.0001)
    sigma.alpha ~ dunif(0,100)
   
   # Derived quantities
    tau <- pow(sigma,-2) # precision

    tau.alpha <- pow(sigma.alpha,-2) # precision
    
    } # end model
    ",fill = TRUE)
  sink()
  
  # Number of groups (i.e., different sites)
  J <- length(unique(jags.dat.2$Site))
  
  # Site indicator
  jags.dat.2$G <- as.numeric(as.factor(as.numeric(as.factor(jags.dat.2$Site))))
  
  #predictor scaled and log transformed for model fit
  ag <- unique(jags.dat.2$ag.logit.z)
  
  #log transform concentrations
  Y <- log(dat[,colIdx])
  n <- length(Y)
  
  # Load data
  data <- list(y = Y, n = n,
               ag=ag,
               group = jags.dat.2$G,
               J = J,
               censorLimitVec = log(jags.dat.2$censorLimitVec),
               isNotCensored = as.numeric(jags.dat.2$isNotCensored) )
  
  ### Initial values ###
  yInit <- rep( NA , nrow(jags.dat.2) )
  ## set to half the censor limit
  yInit[is.na(Y)] <- jags.dat.2$censorLimitVec[is.na(Y)]*0.5
  
  # Log-transform initial values for censored data
  yInit <- log(yInit)
  
  # Check initial values
  head(yInit)
  
  # Initial values
  inits <- function (){
    list (alpha.gamma = rnorm(1),
          sigma=runif(1),
          beta.gamma=rnorm(1),
          sigma.alpha=runif(1),
          y=yInit)
    
  }
  
  # Parameters to be monitored
  parameters <- c("alpha.gamma",
                  "sigma",
                  "beta.gamma",
                  "sigma.alpha",
                  "mu.alpha",
                  "alpha")
  
  # MCMC settings
  ni <- 30000
  nt <- 1
  nb <- 10000
  nc <- 3
  
  # Call JAGS from R 
  out <- jags(data, inits, parameters, "vary.int.fixed.slope.txt", n.chains = nc, 
              n.thin = nt, n.iter = ni, n.burnin = nb)
  
  # Summarize posteriors
  print(out, dig = 3)
  
  ### check for convergence ###
  
  #traceplot(out)

  #convergence and rhat <1.1 
  
  # See what max Rhat value is
  max(out$summary[, c("Rhat")])
  
  ### check for significance ###
  print('Quantiles:')
  print('AG IMMEDIATE OVARY')
  print(dat_col_names[colIdx])
  print('Mu.Alpha:')
  print(quantile(out$sims.list$mu.alpha, c(0.025, 0.975)))
  print('Alpha:')
  print(quantile(out$sims.list$alpha, c(0.025, 0.975)))
  print('Beta.gamma:')
  print(quantile(out$sims.list$beta.gamma, c(0.025, 0.975)))
  
  #overlaps 0 = relationship is NOT significant 
  
  
  ########## Save output for graphing later ############
  
  ## extract mu.alpha matrix
  mu.alpha.ag <- out$sims.list$mu.alpha
  ## convert to dataframe (df)
  mu.alpha.ag <- as.data.frame(mu.alpha.ag)
  ## grab mean of each column and credible intervals
  meanPopAve.ag <- apply(mu.alpha.ag, 2, mean)
  upperCI.PopAve.ag <- apply(mu.alpha.ag, 2, quantile, probs=c(0.975) )
  lowerCI.PopAve.ag <- apply(mu.alpha.ag, 2, quantile, probs=c(0.025) )
  ## combine into df for plotting
  df.ag <- rbind(meanPopAve.ag, upperCI.PopAve.ag,lowerCI.PopAve.ag,ag)
  df.ag <- as.data.frame(t(df.ag))
  
  ## extract Alpha matrix
  alpha.ag <- out$sims.list$alpha
  ## convert to df
  alpha.ag <- as.data.frame(alpha.ag)
  ## grab mean of each column (mean site concentration) and credible intervals
  meanSiteAve.ag <- apply(alpha.ag, 2, mean)
  upperCI.SiteAve.ag <- apply(alpha.ag, 2, quantile, probs=c(0.975) )
  lowerCI.SiteAve.ag <- apply(alpha.ag, 2, quantile, probs=c(0.025) )
  ## combine into df for plotting
  df.ag.2 <- rbind(meanSiteAve.ag, upperCI.SiteAve.ag,lowerCI.SiteAve.ag,ag)
  df.ag.2 <- as.data.frame(t(df.ag.2))
  
  ##### Probability of a positive effect #####
  
  ## pull out beta.gammas
  beta.gamma.ag <- out$sims.list$beta.gamma
  beta.gamma.ag <- as.data.frame(beta.gamma.ag)
  
  #create function to look at the mean of only the positive relationships of the posterior iterations
  probs.fun <- function(x){mean(x>0)}
  
  #apply function and make into dataframe
  
  probs.ag <- apply(beta.gamma.ag,2,probs.fun)
  print('Prob Pos Effect AG')
  print(dat_col_names[colIdx])
  print(probs.ag)

  ## fill in probability of a positive matrix
  Immed_Ovary_Prob_Pos[colIdx-2,1] <- round(probs.ag,2)
  
  ## fill in values with mean post dist (95% CI) to empty matrix 
  avg.beta.gamma.ag <- mean(out$sims.list$beta.gamma)
  beta.lower.ag <- quantile(out$sims.list$beta.gamma, 0.025)
  beta.upper.ag <- quantile(out$sims.list$beta.gamma, 0.975)
  Immed_Ovary_effect[colIdx-2,1] <- paste(round(avg.beta.gamma.ag,2),"[",round(beta.lower.ag,2),",", round(beta.upper.ag,2),"]")
  
  
  #################################################
  ################### DEV LANDUSE #################  
  #################################################
  
  ## print to keep track of model progress
  print('DEVELOPED IMMEDIATE OVARY')
  print(dat_col_names[colIdx])
  
  # Define the model in the BUGS language and write a text file
  sink("vary.int.fixed.slope.txt")
  cat("
    model {
    
    
    # Likelihood: 
    # Level-1 of the model
    for (i in 1:n){ 
    isNotCensored[i] ~ dinterval(y[i] , censorLimitVec[i] )
    y[i] ~ dnorm(mu[i], tau)               
    mu[i] <- alpha[group[i]]   
    } 
    
    
    # Level-2 of the model
    for(j in 1:J){
    alpha[j] ~ dnorm(mu.alpha[j],tau.alpha)
    mu.alpha[j] <- alpha.gamma + beta.gamma *dev[j]
    }
    
    
    # Priors and derived quantities
    sigma ~ dunif(0, 100)
    alpha.gamma ~ dnorm(0, 0.0001)
    beta.gamma ~ dnorm(0,0.0001)
    sigma.alpha ~ dunif(0,100)
   
   # Derived quantities
    tau <- pow(sigma,-2) # precision

    tau.alpha <- pow(sigma.alpha,-2) # precision
    
    } # end model
    ",fill = TRUE)
  sink()
  
  # Number of groups (different sites)
  J <- length(unique(jags.dat.2$Site))
  
  # Site indicator
  jags.dat.2$G <- as.numeric(as.factor(as.numeric(as.factor(jags.dat.2$Site))))
  
  # Predictor scaled and log transformed for model fit
  dev <- unique(jags.dat.2$dev.logit.z)
  
  # log transform concentrations
  Y <- log(dat[,colIdx])
  # number of observations
  n <- length(Y)
  
  # Load data
  data <- list(y = Y, n = n,
               dev=dev,
               group = jags.dat.2$G,
               J = J,
               censorLimitVec = log(jags.dat.2$censorLimitVec),
               isNotCensored = as.numeric(jags.dat.2$isNotCensored) )
  
  ### Initial values ###
  yInit <- rep( NA , nrow(jags.dat.2) )
  ## set to half the censor limit
  yInit[is.na(Y)] <- jags.dat.2$censorLimitVec[is.na(Y)]*0.5
  
  # Log-transform initial values for censored data
  yInit <- log(yInit)
  
  # Check initial values
  head(yInit)
  
  # Initial values
  inits <- function (){
    list (alpha.gamma = rnorm(1),
          sigma=runif(1),
          beta.gamma=rnorm(1),
          sigma.alpha=runif(1),
          y=yInit)
    
  }
  
  # Parameters to be monitored
  parameters <- c("alpha.gamma",
                  "sigma",
                  "beta.gamma",
                  "sigma.alpha",
                  "mu.alpha",
                  "alpha")
  
  # MCMC settings
  ni <- 30000
  nt <- 1
  nb <- 10000
  nc <- 3
  
  # Call JAGS from R 
  out <- jags(data, inits, parameters, "vary.int.fixed.slope.txt", n.chains = nc, 
              n.thin = nt, n.iter = ni, n.burnin = nb)
  
  # Summarize posteriors
  print(out, dig = 3)
  
  ### check for convergence ###
  
  #traceplot(out)
  
  #convergence and rhat <1.1 
  # See what max Rhat value is
  max(out$summary[, c("Rhat")])
  
  ### check for significance ###
  print('Quantiles:')
  print('DEVELOPED IMMEDIATE OVARY')
  print(dat_col_names[colIdx])
  print('Mu.Alpha:')
  print(quantile(out$sims.list$mu.alpha, c(0.025, 0.975)))
  print('Beta.gamma:')
  print(quantile(out$sims.list$beta.gamma, c(0.025, 0.975)))
  
  #overlaps 0 = relationship is NOT significant 
  
  ########## Save output for graphing later ############
  
  ## extract mu.alpha matrix
  mu.alpha.dev <- out$sims.list$mu.alpha
  ## convert to df
  mu.alpha.dev <- as.data.frame(mu.alpha.dev)
  ## grab mean of each column and credible intervals
  meanPopAve.dev <- apply(mu.alpha.dev, 2, mean)
  upperCI.PopAve.dev <- apply(mu.alpha.dev, 2, quantile, probs=c(0.975) )
  lowerCI.PopAve.dev <- apply(mu.alpha.dev, 2, quantile, probs=c(0.025) )
  ## combine into df for plotting
  df.dev <- rbind(meanPopAve.dev, upperCI.PopAve.dev,lowerCI.PopAve.dev,dev)
  df.dev <- as.data.frame(t(df.dev))
  
  ## extract Alpha matrix
  alpha.dev <- out$sims.list$alpha
  ## convert to df
  alpha.dev <- as.data.frame(alpha.dev)
  ## grab mean of each column (mean site concentration) and credible intervals
  meanSiteAve.dev <- apply(alpha.dev, 2, mean)
  upperCI.SiteAve.dev <- apply(alpha.dev, 2, quantile, probs=c(0.975) )
  lowerCI.SiteAve.dev <- apply(alpha.dev, 2, quantile, probs=c(0.025) )
  ## combine into df for plotting
  df.dev.2 <- rbind(meanSiteAve.dev, upperCI.SiteAve.dev,lowerCI.SiteAve.dev,dev)
  df.dev.2 <- as.data.frame(t(df.dev.2))
  
  ##### Probability of a positive effect #####
  ## pull out beta.gammas
  beta.gamma.dev <- out$sims.list$beta.gamma
  beta.gamma.dev <- as.data.frame(beta.gamma.dev)
  
  #create function to look at the mean of only the positive relationships of the posterior iterations
  probs.fun <- function(x){mean(x>0)}
  #apply function and make into dataframe
  
  probs.dev <- apply(beta.gamma.dev,2,probs.fun)
  print('Prob Pos Effect DEV')
  print(dat_col_names[colIdx])
  print(probs.dev)
  
  ## fill in  prob pos matrix
  Immed_Ovary_Prob_Pos[colIdx-2,2] <- round(probs.dev,2)
  
  ## fill in values with mean post dist (95% CI) to empty matrix 
  avg.beta.gamma.dev <- mean(out$sims.list$beta.gamma)
  beta.lower.dev <- quantile(out$sims.list$beta.gamma, 0.025)
  beta.upper.dev <- quantile(out$sims.list$beta.gamma, 0.975)
  Immed_Ovary_effect[colIdx-2,2] <- paste(round(avg.beta.gamma.dev,2),"[",round(beta.lower.dev,2),",", round(beta.upper.dev,2),"]")
  
  #####################################
  ########## BMP INTENSITY ###########
  #####################################
  
  ## print to track of model progress
  print('BMPs IMMEDIATE OVARY')
  print(dat_col_names[colIdx])
  
  # Define the model in the BUGS language and write a text file
  sink("vary.int.fixed.slope.txt")
  cat("
    model {
    
    
    # Likelihood: 
    # Level-1 of the model
    for (i in 1:n){ 
    isNotCensored[i] ~ dinterval(y[i] , censorLimitVec[i] )
    y[i] ~ dnorm(mu[i], tau)               
    mu[i] <- alpha[group[i]]   
    } 
    
    
    # Level-2 of the model
    for(j in 1:J){
    alpha[j] ~ dnorm(mu.alpha[j],tau.alpha)
    mu.alpha[j] <- alpha.gamma + beta.gamma *BMP[j]
    }
    
    
    # Priors and derived quantities
    sigma ~ dunif(0, 100)
    alpha.gamma ~ dnorm(0, 0.0001)
    beta.gamma ~ dnorm(0,0.0001)
    sigma.alpha ~ dunif(0,100)
   
   # Derived quantities
    tau <- pow(sigma,-2) # precision

    tau.alpha <- pow(sigma.alpha,-2) # precision
    
    } # end model
    ",fill = TRUE)
  sink()
  
  # Number of groups (different sites)
  J <- length(unique(jags.dat.2$Site))
  
  # Site indicator
  jags.dat.2$G <- as.numeric(as.factor(as.numeric(as.factor(jags.dat.2$Site))))
  
   
  # Predictor scaled and log transformed for model fit
  BMP <- unique(jags.dat.2$BMP.z)
  
  # log transform concentrations
  Y <- log(dat[,colIdx])
  
  # number of observations
  n <- length(Y)
  
  # Load data
  data <- list(y = Y, n = n,
               BMP=BMP,
               group = jags.dat.2$G,
               J = J,
               censorLimitVec = log(jags.dat.2$censorLimitVec),
               isNotCensored = as.numeric(jags.dat.2$isNotCensored) )
  
  ### Initial values ###
  yInit <- rep( NA , nrow(jags.dat.2) )
  ## set to half the censor limit
  yInit[is.na(Y)] <- jags.dat.2$censorLimitVec[is.na(Y)]*0.5
  
  # Log-transform initial values for censored data
  yInit <- log(yInit)
  
  # check initial values
  head(yInit)
  
  # Initial values
  inits <- function (){
    list (alpha.gamma = rnorm(1),
          sigma=runif(1),
          beta.gamma=rnorm(1),
          sigma.alpha=runif(1),
          y=yInit)
    
  }
  
  # Parameters monitored
  parameters <- c("alpha.gamma",
                  "sigma",
                  "beta.gamma",
                  "sigma.alpha",
                  "mu.alpha",
                  "alpha")
  
  # MCMC settings
  ni <- 30000
  nt <- 1
  nb <- 10000
  nc <- 3
  
  # Call JAGS from R 
  out <- jags(data, inits, parameters, "vary.int.fixed.slope.txt", n.chains = nc, 
              n.thin = nt, n.iter = ni, n.burnin = nb)
  
  # Summarize posteriors
  print(out, dig = 3)
  
  ### check for convergence ###
  
  #traceplot(out)
  
  #convergence and rhat <1.1 
  # See what max Rhat value is
  max(out$summary[, c("Rhat")])
  
  ### check for significance ###
  print('Quantiles:')
  print('BMPs IMMEDIATE OVARY')
  print(dat_col_names[colIdx])
  print('Mu.alpha:')
  print(quantile(out$sims.list$mu.alpha, c(0.025, 0.975)))
  print('Beta.gamma:')
  print(quantile(out$sims.list$beta.gamma, c(0.025, 0.975)))
  
  #overlaps 0 = relatitonship is NOT significant 
  
  ########## Save output for graphing later ############
  
  ## extract mu.alpha matrix
  mu.alpha.BMP <- out$sims.list$mu.alpha
  ## convert to df
  mu.alpha.BMP <- as.data.frame(mu.alpha.BMP)
  ## grab mean of each column and credible intervals
  meanPopAve.BMP <- apply(mu.alpha.BMP, 2, mean)
  upperCI.PopAve.BMP <- apply(mu.alpha.BMP, 2, quantile, probs=c(0.975) )
  lowerCI.PopAve.BMP <- apply(mu.alpha.BMP, 2, quantile, probs=c(0.025) )
  ## combine into df for plotting
  df.BMP <- rbind(meanPopAve.BMP, upperCI.PopAve.BMP,lowerCI.PopAve.BMP,BMP)
  df.BMP <- as.data.frame(t(df.BMP))
  
  ## extract Alpha matrix
  alpha.BMP <- out$sims.list$alpha
  ## convert to df
  alpha.BMP <- as.data.frame(alpha.BMP)
  ## grab mean of each column (mean site concentration) and credible intervals
  meanSiteAve.BMP <- apply(alpha.BMP, 2, mean)
  upperCI.SiteAve.BMP <- apply(alpha.BMP, 2, quantile, probs=c(0.975) )
  lowerCI.SiteAve.BMP <- apply(alpha.BMP, 2, quantile, probs=c(0.025) )
  ## combine into df for plotting
  df.BMP.2 <- rbind(meanSiteAve.BMP, upperCI.SiteAve.BMP,lowerCI.SiteAve.BMP,BMP)
  df.BMP.2 <- as.data.frame(t(df.BMP.2))
  
  ##### Probability of a positive effect #####
  ## pull out beta.gammas
  beta.gamma.BMP <- out$sims.list$beta.gamma
  beta.gamma.BMP <- as.data.frame(beta.gamma.BMP)
  
  #create function to look at the mean of only the positive relationships of the posterior iterations
  probs.fun <- function(x){mean(x>0)}
  #apply function and make into dataframe
  
  probs.BMP <- apply(beta.gamma.BMP,2,probs.fun)
  print('Prob Pos Effect BMPs')
  print(dat_col_names[colIdx])
  print(probs.BMP)
  
  ## fill in prob pos matrix
  Immed_Ovary_Prob_Pos[colIdx-2,3] <- round(probs.BMP,2)
  
  ## fill in values with mean post dist (95% CI) to empty matrix
  avg.beta.gamma.BMP <- mean(out$sims.list$beta.gamma)
  beta.lower.BMP <- quantile(out$sims.list$beta.gamma, 0.025)
  beta.upper.BMP <- quantile(out$sims.list$beta.gamma, 0.975)
  Immed_Ovary_effect[colIdx-2,3] <- paste(round(avg.beta.gamma.BMP,2),"[",round(beta.lower.BMP,2),",", round(beta.upper.BMP,2),"]")
  
  
  ############################################
  ###### Graph all 3 model outputs and #######
  ###### put all 3 figs into one panel ######
  ############################################
  
  
  ### ag ###
  p.ag<- ggplot() +
    geom_ribbon(df.ag,mapping = aes(x=ag, ymin=lowerCI.PopAve.ag, ymax=upperCI.PopAve.ag),fill = "grey90")+
    geom_line(df.ag, mapping =  aes(x=ag, y=meanPopAve.ag)) +
    geom_point(df.ag.2, mapping = aes(x=ag, y=meanSiteAve.ag), size = 2)+
    geom_errorbar(df.ag.2,mapping = aes(x=ag, ymin=lowerCI.SiteAve.ag, ymax=upperCI.SiteAve.ag), width=0.05) +
    theme_classic() +
    xlab("Standardized, Logit[Agricultural Land]") +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size=10, vjust = 0.1))
  
  
  ### dev ###
  p.dev<- ggplot() + 
    geom_ribbon(df.dev,mapping = aes(x=dev, ymin=lowerCI.PopAve.dev, ymax=upperCI.PopAve.dev),fill = "grey90")+
    geom_line(df.dev, mapping =  aes(x=dev, y=meanPopAve.dev)) +
    geom_point(df.dev.2, mapping = aes(x=dev, y=meanSiteAve.dev), size =2)+
    geom_errorbar(df.dev.2,mapping = aes(x=dev, ymin=lowerCI.SiteAve.dev, ymax=upperCI.SiteAve.dev), width=0.05) +
    theme_classic() +
    xlab("Standardized, Logit[Developed Land]") +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size=10, vjust = 0.1))
  
  
  ### BMPs ###
  p.BMP<- ggplot() + 
    geom_ribbon(df.BMP,mapping = aes(x=BMP, ymin=lowerCI.PopAve.BMP, ymax=upperCI.PopAve.BMP), fill = "grey90") +
    geom_line(df.BMP, mapping =  aes(x=BMP, y=meanPopAve.BMP)) +
    geom_point(df.BMP.2, mapping = aes(x=BMP, y=meanSiteAve.BMP), size = 2)+
    geom_errorbar(df.BMP.2,mapping = aes(x=BMP, ymin=lowerCI.SiteAve.BMP, ymax=upperCI.SiteAve.BMP), width=0.05) +
    theme_classic() +
    xlab("Standardized [BMP Intensity]") +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size=10, vjust = 0.1))
  
  
  ### arrange all figures onto one plot ###
  name_figure <- paste(dat_col_names[colIdx],"immed_ovary.pdf",sep="_")
  chemical <- paste(dat_col_names[colIdx], "Immediate Catchment")
  
  p.all <- grid.arrange(
    p.ag,
    p.dev,
    p.BMP,
    ncol = 1,
    left = text_grob("Mean Site Concentration log(ng/g)",
                     size =12, hjust=0.5, vjust=0.3, rot=90))
  
 ggsave(filename = name_figure,
        plot = p.all,
        height = 8,
        width = 6.5,
        units = "in")
  
} ### END FOR-LOOP (of contaminants)






## re-name columns and rows and remove NAs
colnames(Immed_Ovary_Prob_Pos) <- c("Agricultural", "Developed","BMPs")
rownames(Immed_Ovary_Prob_Pos) <- dat_col_names[3:20]
Immed_Ovary_Prob_Pos <- na.omit(Immed_Ovary_Prob_Pos)


## re-name columns and rows and remove NAs
colnames(Immed_Ovary_effect) <- c("Agricultural", "Developed","BMPs")
rownames(Immed_Ovary_effect) <- dat_col_names[3:20]
Immed_Ovary_effect <- na.omit(Immed_Ovary_effect)

## export df 

write.csv(Immed_Ovary_Prob_Pos, "Immed_Ovary_Prob_Pos_summary.csv")
write.csv(Immed_Ovary_effect, "Immed_Ovary_effect_summary.csv")














  

##################################################################################  
####################### JUVENILE (YOY) TISSUE ####################################
##################################################################################

# clear workspace
 rm(list=ls())

##### NOTE: work flow below is the same as above but for juvenile data set #####

###### ###### ###### ######
## load contaminant data ###
###### ###### ###### ###### 

dat1 <- read.csv('All_Compounds_Cleaned.csv')

#rename Site ID to Site and  Tissue Type to Tissue
names(dat1)[4]<-"Site"
names(dat1)[5]<-"Tissue"

## subset out compounds of interest and their corresponding detection limits (DLs), and Site and Tissue
dat <- dat1 %>%
  select('Site',
         'Tissue',
         'Total_PCB',
         'Total_PCB_DL',
         'Hexachlorobenzene','Hexachlorobenzene_DL',
         'p.p..DDE','p.p..DDE_DL',
         'Galaxolide','Galaxolide_DL',
         'Naphthalene','Naphthalene_DL',
         'Phenanthrene','Phenanthrene_DL',
         'Methyl.Triclosan','Methyl_Triclosan_DL',
         'Mirex','Mirex_DL',
         'Pyrene','Pyrene_DL')


###### ###### ###### ######  
## load land use data ###### 
###### ###### ###### ###### 

land_use <- fread('Upstream_Immediate_totals.csv')
unique(land_use$Site_ID)

## order by Site_ID
land_use <- land_use[order(Site_ID)]

## pull out columns of interest
land_use_2 <- land_use[,c('Site_ID','Scale','Developed_total', 'BMP_Intensity', 'Cultivated_total')]

## create columns for landuse data in contaminant data set
dat$PCT.AG <- 1:186
dat$PCT.DEV <- 1:186
dat$BMP <- 1:186

###### ###### 
#### AG ###### 
###### ###### 

## fill in columns (in contaminant data set) with appropriate agricultural landuse value for each site

dat$PCT.AG <- ifelse(dat$Site=='ANT', as.numeric((land_use_2[2,5])),dat$PCT.AG) #ANT
dat$PCT.AG <- ifelse(dat$Site=='BE', as.numeric((land_use_2[4,5])),dat$PCT.AG) #BE
dat$PCT.AG <- ifelse(dat$Site=='CHIL', as.numeric((land_use_2[6,5])),dat$PCT.AG) #CHIL
dat$PCT.AG <- ifelse(dat$Site=='LOY', as.numeric((land_use_2[8,5])),dat$PCT.AG) #LOY
dat$PCT.AG <- ifelse(dat$Site=='MON', as.numeric((land_use_2[10,5])),dat$PCT.AG) #MON
dat$PCT.AG <- ifelse(dat$Site=='PC', as.numeric((land_use_2[12,5])),dat$PCT.AG) #PC
dat$PCT.AG <- ifelse(dat$Site=='POT', as.numeric((land_use_2[14,5])),dat$PCT.AG) #POT
dat$PCT.AG <- ifelse(dat$Site=='SUSHA', as.numeric((land_use_2[16,5])),dat$PCT.AG) #SUSHA
dat$PCT.AG <- ifelse(dat$Site=='WBM', as.numeric((land_use_2[18,5])),dat$PCT.AG) #WMB
dat$PCT.AG <- ifelse(dat$Site=='WBMC', as.numeric((land_use_2[20,5])),dat$PCT.AG) #WBMC
dat$PCT.AG <- ifelse(dat$Site=='WYA', as.numeric((land_use_2[22,5])),dat$PCT.AG) #WYA

dat$PCT.AG <- as.numeric(dat$PCT.AG)

#### standardize and logit-transfrom ###### 

dat$ag.logit.z <- as.numeric(scale(logit(dat$PCT.AG)))
#hist(dat$ag.logit.z)

###### ###### 
#### DEV ######  
###### ###### 

## fill in columns (in contaminant data set) with appropriate developed landuse value for each site

dat$PCT.DEV <- ifelse(dat$Site=='ANT', as.numeric((land_use_2[2,3])),dat$PCT.DEV) #ANT
dat$PCT.DEV <- ifelse(dat$Site=='BE', as.numeric((land_use_2[4,3])),dat$PCT.DEV) #BE
dat$PCT.DEV <- ifelse(dat$Site=='CHIL', as.numeric((land_use_2[6,3])),dat$PCT.DEV) #CHIL
dat$PCT.DEV <- ifelse(dat$Site=='LOY', as.numeric((land_use_2[8,3])),dat$PCT.DEV) #LOY
dat$PCT.DEV <- ifelse(dat$Site=='MON', as.numeric((land_use_2[10,3])),dat$PCT.DEV) #MON
dat$PCT.DEV <- ifelse(dat$Site=='PC', as.numeric((land_use_2[12,3])),dat$PCT.DEV) #PC
dat$PCT.DEV <- ifelse(dat$Site=='POT', as.numeric((land_use_2[14,3])),dat$PCT.DEV) #POT
dat$PCT.DEV <- ifelse(dat$Site=='SUSHA', as.numeric((land_use_2[16,3])),dat$PCT.DEV) #SUSHA
dat$PCT.DEV <- ifelse(dat$Site=='WBM', as.numeric((land_use_2[18,3])),dat$PCT.DEV) #WMB
dat$PCT.DEV <- ifelse(dat$Site=='WBMC', as.numeric((land_use_2[20,3])),dat$PCT.DEV) #WBMC
dat$PCT.DEV <- ifelse(dat$Site=='WYA', as.numeric((land_use_2[22,3])),dat$PCT.DEV) #WYA

dat$PCT.DEV <- as.numeric(dat$PCT.DEV)

#### standardize and logit-transfrom ###### 

dat$dev.logit.z <- as.numeric(scale(logit(dat$PCT.DEV)))
#hist(dat$dev.logit.z)

###### ###### ###### 
#### BMPS ###### 
###### ###### ###### 

## fill in columns (in contaminant data set) with appropriate BMP intensity value for each site

dat$BMP <- ifelse(dat$Site=='ANT', as.numeric((land_use_2[1,4])),dat$BMP) #ANT
dat$BMP <- ifelse(dat$Site=='BE', as.numeric((land_use_2[3,4])),dat$BMP) #BE
dat$BMP <- ifelse(dat$Site=='CHIL', as.numeric((land_use_2[5,4])),dat$BMP) #CHIL
dat$BMP <- ifelse(dat$Site=='LOY', as.numeric((land_use_2[7,4])),dat$BMP) #LOY
dat$BMP <- ifelse(dat$Site=='MON', as.numeric((land_use_2[9,4])),dat$BMP) #MON
dat$BMP <- ifelse(dat$Site=='PC', as.numeric((land_use_2[11,4])),dat$BMP) #PC
dat$BMP <- ifelse(dat$Site=='POT', as.numeric((land_use_2[13,4])),dat$BMP) #POT
dat$BMP <- ifelse(dat$Site=='SUSHA', as.numeric((land_use_2[15,4])),dat$BMP) #SUSHA
dat$BMP <- ifelse(dat$Site=='WBM', as.numeric((land_use_2[17,4])),dat$BMP) #WMB
dat$BMP <- ifelse(dat$Site=='WBMC', as.numeric((land_use_2[19,4])),dat$BMP) #WBMC
dat$BMP <- ifelse(dat$Site=='WYA', as.numeric((land_use_2[21,4])),dat$BMP) #WYA

dat$BMP <- as.numeric(dat$BMP)

#### standardize and logit-transfrom ###### 

dat$BMP.z <- as.numeric(scale(dat$BMP))
#hist(dat$BMP.logit.z)

#### subset for tissue of interest ###### 

dat <- subset(dat, Tissue =="YOY")

###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ######
###### ###### ###### ###### FOR-LOOP FOR JUVENILE (YOY) TISSUE MODEL ###### ###### 
###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### ###### 

## pull out columns for future use
dat_length <- ncol(dat)

## pull out columne names for future use
dat_col_names <- colnames(dat)

## create matrix to fill in with the probability of a positive effect and estimated effects
Immed_YOY_Prob_Pos <- matrix(nrow = (dat_length-8), ncol = 3)
Immed_YOY_effect <- matrix(nrow=(dat_length-8), ncol = 3)


for (colIdx in seq(from=3, to=dat_length-7, by=2)) { ### BEGIN FOR-LOOP (looping through contaminants)
  
  ## print to keep track of model progress
  print("NEW LOOP ITERATION.  WORKING ON:")
  print(dat_col_names[colIdx])
  
  # Set the appropriate Censoring Limit for contaminant 
  dat$censorLimitVec <- dat[,colIdx+1]
  
  # Must tell JAGS which observations are ABOVE Censoring limit
  dat$isNotCensored <- (dat[,colIdx] > dat$censorLimitVec) # Returns TRUE if not below DL and NA if below DL
  
  # Convert NA in dat$isNotCensored to FALSE
  dat$isNotCensored <- ifelse(is.na(dat$isNotCensored), FALSE, dat$isNotCensored)
  
  #################################################################
  ########## BUGS CODE ############################################
  #################################################################
  
  ## copy of data set
  jags.dat.2 <- copy(dat)
  
  ## load jagsUI package ##
  library(jagsUI)
  
  ####################################################
  ################### AG LANDUSE  ####################
  ####################################################
  
  print('AG JUVENILE')
  print(dat_col_names[colIdx])
  
  # Define the model in the BUGS language and write a text file
  sink("vary.int.fixed.slope.txt")
  cat("
    model {
    
    
    # Likelihood: 
    # Level-1 of the model
    for (i in 1:n){ 
    isNotCensored[i] ~ dinterval(y[i] , censorLimitVec[i] )
    y[i] ~ dnorm(mu[i], tau)               
    mu[i] <- alpha[group[i]]   
    } 
    
    
    # Level-2 of the model
    for(j in 1:J){
    alpha[j] ~ dnorm(mu.alpha[j],tau.alpha)
    mu.alpha[j] <- alpha.gamma + beta.gamma *ag[j]
    }
    
    
    # Priors and derived quantities
    sigma ~ dunif(0, 100)
    alpha.gamma ~ dnorm(0, 0.0001)
    beta.gamma ~ dnorm(0,0.0001)
    sigma.alpha ~ dunif(0,100)
   
   # Derived quantities
    tau <- pow(sigma,-2) # precision

    tau.alpha <- pow(sigma.alpha,-2) # precision
    
    } # end model
    ",fill = TRUE)
  sink()
  
  # Number of groups (different sites)
  J <- length(unique(jags.dat.2$Site))
  
  # Site indicator
  jags.dat.2$G <- as.numeric(as.factor(as.numeric(as.factor(jags.dat.2$Site))))
  
  #predictor scaled and log transformed for model fit
  ag <- unique(jags.dat.2$ag.logit.z)
  
  #log transform concentrations
  Y <- log(dat[,colIdx])
  n <- length(Y)
  
  # Load data
  data <- list(y = Y, n = n,
               ag=ag,
               group = jags.dat.2$G,
               J = J,
               censorLimitVec = log(jags.dat.2$censorLimitVec),
               isNotCensored = as.numeric(jags.dat.2$isNotCensored) )
  
  ### Initial values ###
  yInit <- rep( NA , nrow(jags.dat.2) )
  ## set to half the censor limit
  yInit[is.na(Y)] <- jags.dat.2$censorLimitVec[is.na(Y)]*0.5
  
  # Log-transform initial values for censored data
  yInit <- log(yInit)
  
  head(yInit)
  # Initial values
  inits <- function (){
    list (alpha.gamma = rnorm(1),
          sigma=runif(1),
          beta.gamma=rnorm(1),
          sigma.alpha=runif(1),
          y=yInit)
    
  }
  
  # Parameters to be monitored
  parameters <- c("alpha.gamma",
                  "sigma",
                  "beta.gamma",
                  "sigma.alpha",
                  "mu.alpha",
                  "alpha")
  
  # MCMC settings
  ni <- 30000
  nt <- 1
  nb <- 10000
  nc <- 3
  
  # Call JAGS from R 
  out <- jags(data, inits, parameters, "vary.int.fixed.slope.txt", n.chains = nc, 
              n.thin = nt, n.iter = ni, n.burnin = nb)
  
  # Summarize posteriors
  print(out, dig = 3)
  
  ### check for convergence ###
  #traceplot(out)
  
  #convergence and rhat <1.1 
  
  # See what max Rhat value is
  max(out$summary[, c("Rhat")])
  
  ### check for significance ###
  print('Quantiles:')
  print('CULTIVATED IMMEDIATE YOY')
  print(dat_col_names[colIdx])
  print('Mu.Alpha:')
  print(quantile(out$sims.list$mu.alpha, c(0.025, 0.975)))
  print('Beta.gamma:')
  print(quantile(out$sims.list$beta.gamma, c(0.025, 0.975)))
  
  #overlaps 0 = relatitonship is NOT significant 
  
  
  ########## Save output for graphing later ############
  
  
  ## extract mu.alpha matrix
  mu.alpha.ag <- out$sims.list$mu.alpha
  ## convert to df
  mu.alpha.ag <- as.data.frame(mu.alpha.ag)
  ## grab mean of each column and credible intervals
  meanPopAve.ag <- apply(mu.alpha.ag, 2, mean)
  upperCI.PopAve.ag <- apply(mu.alpha.ag, 2, quantile, probs=c(0.975) )
  lowerCI.PopAve.ag <- apply(mu.alpha.ag, 2, quantile, probs=c(0.025) )
  ## combine into df for plotting
  df.ag <- rbind(meanPopAve.ag, upperCI.PopAve.ag,lowerCI.PopAve.ag,ag)
  df.ag <- as.data.frame(t(df.ag))
  
  ## extract Alpha matrix
  alpha.ag <- out$sims.list$alpha
  ## convert to df
  alpha.ag <- as.data.frame(alpha.ag)
  ## grab mean of each column (mean site concentration) and credible intervals
  meanSiteAve.ag <- apply(alpha.ag, 2, mean)
  upperCI.SiteAve.ag <- apply(alpha.ag, 2, quantile, probs=c(0.975) )
  lowerCI.SiteAve.ag <- apply(alpha.ag, 2, quantile, probs=c(0.025) )
  ## combine into df for plotting
  df.ag.2 <- rbind(meanSiteAve.ag, upperCI.SiteAve.ag,lowerCI.SiteAve.ag,ag)
  df.ag.2 <- as.data.frame(t(df.ag.2))
  
  ##### Probability of a positive effect #####
  
  ## pull out beta.gammas
  beta.gamma.ag <- out$sims.list$beta.gamma
  beta.gamma.ag <- as.data.frame(beta.gamma.ag)
  
  #create function to look at the mean of only the positive relationships of the posterior iterations
  probs.fun <- function(x){mean(x>0)}
  #apply function and make into dataframe
  
  probs.ag <- apply(beta.gamma.ag,2,probs.fun)
  print('Prob Pos Effect AG')
  print(dat_col_names[colIdx])
  print(probs.ag)
  
  ## fill in matrix with prob pos
  Immed_YOY_Prob_Pos[colIdx-2,1] <- round(probs.ag,2)
  
  ## fill in values with mean post dist (95% CI) to empty matrix 
  avg.beta.gamma.ag <- mean(out$sims.list$beta.gamma)
  beta.lower.ag <- quantile(out$sims.list$beta.gamma, 0.025)
  beta.upper.ag <- quantile(out$sims.list$beta.gamma, 0.975)
  Immed_YOY_effect[colIdx-2,1] <- paste(round(avg.beta.gamma.ag,2),"[",round(beta.lower.ag,2),",", round(beta.upper.ag,2),"]")
  
  ##################################################
  ################### DEV LANDUSE ##################
  ##################################################
  
  # print to keep track of model progress
  print('DEVELOPED JUVENILE')
  print(dat_col_names[colIdx])
  
  # Define the model in the BUGS language and write a text file
  sink("vary.int.fixed.slope.txt")
  cat("
    model {
    
    
    # Likelihood: 
    # Level-1 of the model
    for (i in 1:n){ 
    isNotCensored[i] ~ dinterval(y[i] , censorLimitVec[i] )
    y[i] ~ dnorm(mu[i], tau)               
    mu[i] <- alpha[group[i]]   
    } 
    
    
    # Level-2 of the model
    for(j in 1:J){
    alpha[j] ~ dnorm(mu.alpha[j],tau.alpha)
    mu.alpha[j] <- alpha.gamma + beta.gamma *dev[j]
    }
    
    
    # Priors and derived quantities
    sigma ~ dunif(0, 100)
    alpha.gamma ~ dnorm(0, 0.0001)
    beta.gamma ~ dnorm(0,0.0001)
    sigma.alpha ~ dunif(0,100)
   
   # Derived quantities
    tau <- pow(sigma,-2) # precision

    tau.alpha <- pow(sigma.alpha,-2) # precision
    
    } # end model
    ",fill = TRUE)
  sink()
  
  # Number of groups (different sites)
  J <- length(unique(jags.dat.2$Site))
  
  # Site indicator
  jags.dat.2$G <- as.numeric(as.factor(as.numeric(as.factor(jags.dat.2$Site))))
  
  #predictor scaled and log transformed for model fit
  dev <- unique(jags.dat.2$dev.logit.z)
  
  #log transform concentrations
  Y <- log(dat[,colIdx])
  n <- length(Y)
  
  # Load data
  data <- list(y = Y, n = n,
               dev=dev,
               group = jags.dat.2$G,
               J = J,
               censorLimitVec = log(jags.dat.2$censorLimitVec),
               isNotCensored = as.numeric(jags.dat.2$isNotCensored) )
  
  ### Initial values ###
  yInit <- rep( NA , nrow(jags.dat.2) )
  ## set to half the censor limit
  yInit[is.na(Y)] <- jags.dat.2$censorLimitVec[is.na(Y)]*0.5
  
  # Log-transform initial values for censored data
  yInit <- log(yInit)
  
  head(yInit)
  
  # Initial values
  inits <- function (){
    list (alpha.gamma = rnorm(1),
          sigma=runif(1),
          beta.gamma=rnorm(1),
          sigma.alpha=runif(1),
          y=yInit)
    
  }
  
  # Parameters to be monitored
  parameters <- c("alpha.gamma",
                  "sigma",
                  "beta.gamma",
                  "sigma.alpha",
                  "mu.alpha",
                  "alpha")
  
  # MCMC settings
  ni <- 30000
  nt <- 1
  nb <- 10000
  nc <- 3
  
  # Call JAGS from R 
  out <- jags(data, inits, parameters, "vary.int.fixed.slope.txt", n.chains = nc, 
              n.thin = nt, n.iter = ni, n.burnin = nb)
  
  # Summarize posteriors
  print(out, dig = 3)
  
  ### check for convergence ###
  
  #traceplot(out)
  
  #convergence and rhat <1.1 
  # See what max Rhat value is
  max(out$summary[, c("Rhat")])
  
  ### check for significance ###
  print('Quantiles:')
  print('DEVELOPED IMMEDIATE YOY')
  print(dat_col_names[colIdx])
  print('Mu.Alpha:')
  print(quantile(out$sims.list$mu.alpha, c(0.025, 0.975)))
  print('Beta.gamma:')
  print(quantile(out$sims.list$beta.gamma, c(0.025, 0.975)))
  
  #overlaps 0 = relatitonship is NOT significant 
  
  ########## Save output for graphing later ############
  
  ## extract mu.alpha matrix
  mu.alpha.dev <- out$sims.list$mu.alpha
  ## convert to dataframe
  mu.alpha.dev <- as.data.frame(mu.alpha.dev)
  ## grab mean of each column and credible intervals
  meanPopAve.dev <- apply(mu.alpha.dev, 2, mean)
  upperCI.PopAve.dev <- apply(mu.alpha.dev, 2, quantile, probs=c(0.975) )
  lowerCI.PopAve.dev <- apply(mu.alpha.dev, 2, quantile, probs=c(0.025) )
  ## combine into df for plotting
  df.dev <- rbind(meanPopAve.dev, upperCI.PopAve.dev,lowerCI.PopAve.dev,dev)
  df.dev <- as.data.frame(t(df.dev))
  
  ## extract Alpha matrix
  alpha.dev <- out$sims.list$alpha
  ## convert to df
  alpha.dev <- as.data.frame(alpha.dev)
  ## grab mean of each column (mean site concentration) and credible intervals
  meanSiteAve.dev <- apply(alpha.dev, 2, mean)
  upperCI.SiteAve.dev <- apply(alpha.dev, 2, quantile, probs=c(0.975) )
  lowerCI.SiteAve.dev <- apply(alpha.dev, 2, quantile, probs=c(0.025) )
  ## combine into df for plotting
  df.dev.2 <- rbind(meanSiteAve.dev, upperCI.SiteAve.dev,lowerCI.SiteAve.dev,dev)
  df.dev.2 <- as.data.frame(t(df.dev.2))
  
  ##### Probability of a positive effect #####
  
  ## pull out beta.gammas
  beta.gamma.dev <- out$sims.list$beta.gamma
  beta.gamma.dev <- as.data.frame(beta.gamma.dev)
  
  #create function to look at the mean of only the positive relationships of the posterior iterations
  probs.fun <- function(x){mean(x>0)}
  #apply function and make into dataframe
  
  probs.dev <- apply(beta.gamma.dev,2,probs.fun)
  print('Prob Pos Effect DEV')
  print(dat_col_names[colIdx])
  print(probs.dev)
  
  ## fill in matrix with prob pos
  Immed_YOY_Prob_Pos[colIdx-2,2] <- round(probs.dev,2)
  
  ## fill in values with mean post dist (95% CI) to empty matrix 
  avg.beta.gamma.dev <- mean(out$sims.list$beta.gamma)
  beta.lower.dev <- quantile(out$sims.list$beta.gamma, 0.025)
  beta.upper.dev <- quantile(out$sims.list$beta.gamma, 0.975)
  Immed_YOY_effect[colIdx-2,2] <- paste(round(avg.beta.gamma.dev,2),"[",round(beta.lower.dev,2),",", round(beta.upper.dev,2),"]")
  
  #####################################################
  ################### BMP INTENSITY ###################
  #####################################################
  
  # print to keep track of model progress
  print('BMPs JUVENILE')
  print(dat_col_names[colIdx])
  
  # Define the model in the BUGS language and write a text file
  sink("vary.int.fixed.slope.txt")
  cat("
    model {
    
    
    # Likelihood: 
    # Level-1 of the model
    for (i in 1:n){ 
    isNotCensored[i] ~ dinterval(y[i] , censorLimitVec[i] )
    y[i] ~ dnorm(mu[i], tau)               
    mu[i] <- alpha[group[i]]   
    } 
    
    
    # Level-2 of the model
    for(j in 1:J){
    alpha[j] ~ dnorm(mu.alpha[j],tau.alpha)
    mu.alpha[j] <- alpha.gamma + beta.gamma *BMP[j]
    }
    
    
    # Priors and derived quantities
    sigma ~ dunif(0, 100)
    alpha.gamma ~ dnorm(0, 0.0001)
    beta.gamma ~ dnorm(0,0.0001)
    sigma.alpha ~ dunif(0,100)
   
   # Derived quantities
    tau <- pow(sigma,-2) # precision

    tau.alpha <- pow(sigma.alpha,-2) # precision
    
    } # end model
    ",fill = TRUE)
  sink()
  
  # Number of groups (different sites)
  J <- length(unique(jags.dat.2$Site))
  
  # Site indicator
  jags.dat.2$G <- as.numeric(as.factor(as.numeric(as.factor(jags.dat.2$Site))))
  
   
  # predictor scaled and log transformed for model fit
  BMP <- unique(jags.dat.2$BMP.z)
   
  
  #log transform concentrations
  Y <- log(dat[,colIdx])
  n <- length(Y)
  
  # Load data
  data <- list(y = Y, n = n,
               BMP=BMP,
               group = jags.dat.2$G,
               J = J,
               censorLimitVec = log(jags.dat.2$censorLimitVec),
               isNotCensored = as.numeric(jags.dat.2$isNotCensored) )
  
  ### Initial values ###
  yInit <- rep( NA , nrow(jags.dat.2) )
  ## set to half the censor limit
  yInit[is.na(Y)] <- jags.dat.2$censorLimitVec[is.na(Y)]*0.5
  
  # Log-transform initial values for censored data
  yInit <- log(yInit)
  
  head(yInit)
  
  # Initial values
  inits <- function (){
    list (alpha.gamma = rnorm(1),
          sigma=runif(1),
          beta.gamma=rnorm(1),
          sigma.alpha=runif(1),
          y=yInit)
    
  }
  
  # Parameters to be monitored
  parameters <- c("alpha.gamma",
                  "sigma",
                  "beta.gamma",
                  "sigma.alpha",
                  "mu.alpha",
                  "alpha")
  
  # MCMC settings
  ni <- 30000
  nt <- 1
  nb <- 10000
  nc <- 3
  
  # Call JAGS from R 
  out <- jags(data, inits, parameters, "vary.int.fixed.slope.txt", n.chains = nc, 
              n.thin = nt, n.iter = ni, n.burnin = nb)
  
  # Summarize posteriors
  print(out, dig = 3)
  
  ### check for convergence ###
  
  #traceplot(out)
  
  #convergence and rhat <1.1 
  # See what max Rhat value is
  max(out$summary[, c("Rhat")])
  
  ### check for significance ###
  print('Quantiles:')
  print('BMPs IMMEDIATE YOY')
  print(dat_col_names[colIdx])
  print('Mu.alpha:')
  print(quantile(out$sims.list$mu.alpha, c(0.025, 0.975)))
  print('Beta.gamma:')
  print(quantile(out$sims.list$beta.gamma, c(0.025, 0.975)))
  
  #overlaps 0 = relatitonship is NOT significant 
  
  ########## Save output for graphing later ############
  
  ## extract mu.alpha matrix
  mu.alpha.BMP <- out$sims.list$mu.alpha
  ## convert to df
  mu.alpha.BMP <- as.data.frame(mu.alpha.BMP)
  ## grab mean of each column and credible intervals
  meanPopAve.BMP <- apply(mu.alpha.BMP, 2, mean)
  upperCI.PopAve.BMP <- apply(mu.alpha.BMP, 2, quantile, probs=c(0.975) )
  lowerCI.PopAve.BMP <- apply(mu.alpha.BMP, 2, quantile, probs=c(0.025) )
  ## combine into df for plotting
  df.BMP <- rbind(meanPopAve.BMP, upperCI.PopAve.BMP,lowerCI.PopAve.BMP,BMP)
  df.BMP <- as.data.frame(t(df.BMP))
  
  ## extract Alpha matrix
  alpha.BMP <- out$sims.list$alpha
  ## convert to df
  alpha.BMP <- as.data.frame(alpha.BMP)
  ## grab mean of each column (mean site concentration) and credible intervals
  meanSiteAve.BMP <- apply(alpha.BMP, 2, mean)
  upperCI.SiteAve.BMP <- apply(alpha.BMP, 2, quantile, probs=c(0.975) )
  lowerCI.SiteAve.BMP <- apply(alpha.BMP, 2, quantile, probs=c(0.025) )
  ## combine into df for plotting
  df.BMP.2 <- rbind(meanSiteAve.BMP, upperCI.SiteAve.BMP,lowerCI.SiteAve.BMP,BMP)
  df.BMP.2 <- as.data.frame(t(df.BMP.2))
  
  ##### Probability of a positive effect #####
  
  ## pull out beta.gammas
  beta.gamma.BMP <- out$sims.list$beta.gamma
  beta.gamma.BMP <- as.data.frame(beta.gamma.BMP)
  
  #create function to look at the mean of only the positive relationships of the posterior iterations
  probs.fun <- function(x){mean(x>0)}
  #apply function and make into dataframe
  
  probs.BMP <- apply(beta.gamma.BMP,2,probs.fun)
  print('Prob Pos Effect BMPs')
  print(dat_col_names[colIdx])
  print(probs.BMP)
  
  ## fill in matrix with prob pos
  Immed_YOY_Prob_Pos[colIdx-2,3] <- round(probs.BMP,2)
  
  ## fill in values with mean post dist (95% CI) to empty matrix
  avg.beta.gamma.BMP <- mean(out$sims.list$beta.gamma)
  beta.lower.BMP <- quantile(out$sims.list$beta.gamma, 0.025)
  beta.upper.BMP <- quantile(out$sims.list$beta.gamma, 0.975)
  Immed_YOY_effect[colIdx-2,3] <- paste(round(avg.beta.gamma.BMP,2),"[",round(beta.lower.BMP,2),",", round(beta.upper.BMP,2),"]")
  
  
  ############################################
  ###### Graph all 3 model outputs and #######
  ###### put all 3 figs into one panel ######
  ############################################
  
  ### ag ###
  p.ag<- ggplot() +
    geom_ribbon(df.ag,mapping = aes(x=ag, ymin=lowerCI.PopAve.ag, ymax=upperCI.PopAve.ag), fill = "grey90")+
    geom_line(df.ag, mapping =  aes(x=ag, y=meanPopAve.ag)) +
    geom_point(df.ag.2, mapping = aes(x=ag, y=meanSiteAve.ag), size =2)+
    geom_errorbar(df.ag.2,mapping = aes(x=ag, ymin=lowerCI.SiteAve.ag, ymax=upperCI.SiteAve.ag), width=0.05) +
    theme_classic() +
    xlab("Standardized, Logit[Agricultural Land]") +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size=10, vjust = 0.1))
  
  
  ### dev ###
  p.dev<- ggplot() + 
    geom_ribbon(df.dev,mapping = aes(x=dev, ymin=lowerCI.PopAve.dev, ymax=upperCI.PopAve.dev), fill = "grey90")+
    geom_line(df.dev, mapping =  aes(x=dev, y=meanPopAve.dev)) +
    geom_point(df.dev.2, mapping = aes(x=dev, y=meanSiteAve.dev), size =2)+
    geom_errorbar(df.dev.2,mapping = aes(x=dev, ymin=lowerCI.SiteAve.dev, ymax=upperCI.SiteAve.dev), width=0.05) +
    theme_classic() +
    xlab("Standardized, Logit[Developed Land]") +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size=10, vjust = 0.1))
  
  
  ### BMPs ###
  p.BMP <- ggplot() + 
    geom_ribbon(df.BMP,mapping = aes(x=BMP, ymin=lowerCI.PopAve.BMP, ymax=upperCI.PopAve.BMP), fill = "grey90") +
    geom_line(df.BMP, mapping =  aes(x=BMP, y=meanPopAve.BMP)) +
    geom_point(df.BMP.2, mapping = aes(x=BMP, y=meanSiteAve.BMP), size =2)+
    geom_errorbar(df.BMP.2,mapping = aes(x=BMP, ymin=lowerCI.SiteAve.BMP, ymax=upperCI.SiteAve.BMP), width=0.05) +
    theme_classic() +
    xlab("Standardized [BMP Intensity]") +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size=10, vjust = 0.1))
  
  
  ### arrange all figures onto one plot ###
  name_figure <- paste(dat_col_names[colIdx],"immed_YOY.pdf",sep="_")
  chemical <- paste(dat_col_names[colIdx], "Immediate Catchment")
  
  p.all <- grid.arrange(
    p.ag,
    p.dev,
    p.BMP,
    ncol = 1,
    left = text_grob("Mean Site Concentration log(ng/g)",
                     size =12, hjust=0.5, vjust=0.3, rot=90))
  
  ggsave(filename = name_figure,
         plot = p.all,
         height = 8,
         width = 6.5,
         units = "in")
  
} ### END FOR-LOOP (of contaminants)






## re-name columns and rows and remove NAs
colnames(Immed_YOY_Prob_Pos) <- c("Agricultural", "Developed","BMPs")
rownames(Immed_YOY_Prob_Pos) <- dat_col_names[3:20]
Immed_YOY_Prob_Pos <- na.omit(Immed_YOY_Prob_Pos)

## re-name columns and rows and remove NAs
colnames(Immed_YOY_effect) <- c("Agricultural", "Developed","BMPs")
rownames(Immed_YOY_effect) <- dat_col_names[3:20]
Immed_YOY_effect <- na.omit(Immed_YOY_effect)

## export df
write.csv(Immed_YOY_Prob_Pos, "Immed_YOY_Prob_Pos_summary.csv")
write.csv(Immed_YOY_effect, "Immed_YOY_effect_summary.csv")



