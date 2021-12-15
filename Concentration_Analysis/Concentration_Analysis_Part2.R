# clear workspace
rm(list=ls())


## load library ##

#library(jagsUI) ## once loaded, XQuartz will be used to view data sets, so wait until ready to run model

library(MCMCpack)
library(lme4)
library(car)
library(plyr)
library(dplyr)
library(data.table) ## for fread
library(ggplot2)
library(gridExtra)
library(ggpubr)
#library(jagsUI)

############# set working directory ###########

#setwd("~/Desktop/Penn State/Classwork/CEC_Project/Chapter 2/R Scripts/Concentration Analyses/Concentration Analysis - Part 2")

####################################################################



###### ###### ###### ######
## load contaminant data ###
###### ###### ###### ###### 

dat1 <- read.csv('All_Compounds_Cleaned.csv')

#rename Site ID to Site and Tissue Type to Tissue
names(dat1)[4]<-"Site"
names(dat1)[5]<-"Tissue"

## subset out compounds of interest and their corresponding detection limits (DLs), and Site and Tissue
dat <- dat1 %>% select('Site',
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

dat$Site <- as.factor(dat$Site)
dat$Tissue <- as.factor(dat$Tissue)

## order to see which site have both YOY and  Ovary data
dat <- dat[order(dat$Site),]

## filter for sites that have both ovary and YOY tissue data ##
## i.e., WBM, POT, ANT, CHIl, BE, PC
dat.2 <- dat[dat$Site ==  "ANT"| dat$Site ==  "BE"| dat$Site == "CHIL" | dat$Site =="POT" | dat$Site == "PC"| dat$Site == "WBM",] 



###### ###### ###### ###### ###### ###### ###### ###### ###
###### ###### FOR-LOOP FOR TISSUE COMPARISON MODEL ###### # 
###### ###### ###### ###### ###### ###### ###### ###### ### 

## save number of columns for future use in model
dat_length <- ncol(dat.2)

## save column names for future use in model
dat_col_names <- colnames(dat.2)

## create empty matrices to fill in with model results for plotting and tables ##

name <- seq(from=3, to= 20, by=2)
Immed_effect_beta <- matrix(nrow =54, ncol = 3)
colnames(Immed_effect_beta) <- c("Mean.Tissue.Diff", "Lower.CI", "Upper.CI")
rownames(Immed_effect_beta) <- c(rep("Total_PCB", 6),
                                  rep("Hexachlorobenzene", 6),
                                  rep("p.p..DDE", 6),
                                  rep("Galaxolide", 6),
                                  rep("Naphthalene", 6),
                                  rep("Phenanthrene", 6),
                                  rep("Methyl.Triclosan", 6),
                                  rep("Mirex", 6),
                                  rep("Pyrene", 6))

Beta.Table <- matrix(nrow =9, ncol = 6)
colnames(Beta.Table) <- c("ANT","BE","CHIL","PC","POT", "WBM")
rownames(Beta.Table) <- dat_col_names[name]

Immed_effect_beta.2 <- matrix(nrow =9, ncol =3)
rownames(Immed_effect_beta.2) <- dat_col_names[name]
colnames(Immed_effect_beta.2) <- c("Average.Diff","Lower.Diff","Upper.Diff")
  
Beta.Table.2 <- matrix(nrow =9, ncol = 1)
rownames(Beta.Table.2) <- dat_col_names[name]
colnames(Beta.Table.2) <- "Beta.Average"

Immed_effect_alpha <- matrix(nrow =54, ncol = 6)
colnames(Immed_effect_alpha) <- c("Ovary.Mean", "Ovary.Lower.CI", "Ovary.Upper.CI","Juvenile.Mean", "Juvenile.Lower.CI", "Juvenile.Upper.CI")
rownames(Immed_effect_alpha) <- c(rep("Total_PCB", 6),
                             rep("Hexachlorobenzene", 6),
                             rep("p.p..DDE", 6),
                             rep("Galaxolide", 6),
                             rep("Naphthalene", 6),
                             rep("Phenanthrene", 6),
                             rep("Methyl.Triclosan", 6),
                             rep("Mirex", 6),
                             rep("Pyrene", 6))

Alpha.Table <- matrix(nrow =9, ncol = 12)
colnames(Alpha.Table) <- c("ANT.Ovary",
                           "BE.Ovary",
                           "CHIL.Ovary",
                           "PC.Ovary",
                           "POT.Ovary",
                           "WBM.Ovary",
                           "ANT.Juvenile",
                           "BE.Juvenile",
                           "CHIL.Juvenile",
                           "PC.Juvenile",
                           "POT.Juvenile",
                           "WBM.Juvenile")
rownames(Alpha.Table) <- dat_col_names[name]

Alpha.Table.2 <- matrix(nrow =9, ncol = 2)
rownames(Alpha.Table.2) <- dat_col_names[name]
colnames(Alpha.Table.2) <- c("Ovary", "Juvenile")

Immed_effect_alpha.2 <- matrix(nrow =9, ncol = 6)
rownames(Immed_effect_alpha.2) <- dat_col_names[name]
colnames(Immed_effect_alpha.2) <- c("Ovary.Mean",
                                    "Ovary.Lower",
                                    "Ovary.Upper",
                                    "Juvenile.Mean",
                                    "Juvenile.Lower",
                                    "Juvenile.Upper")

## matrix to print r-hats for checking convergence of each model

r_hats <- matrix(nrow = 9, ncol = 1)
rownames(r_hats) <- c("Total_PCB",
                      "Hexachlorobenzene",
                      "p.p..DDE",
                      "Galaxolide",
                      "Naphthalene",
                      "Phenanthrene",
                      "Methyl.Triclosan",
                      "Mirex",
                      "Pyrene")

colnames(r_hats) <- "R Hats"

##for tables
rowIdx <-1 
##for plotting
rowIdx2 <- 1

for (colIdx in seq(from=3, to=dat_length, by=2)) {  ### BEGIN FOR-LOOP
  
  ## print to track model progression
  print("NEW LOOP ITERATION.  WORKING ON:")
  print(dat_col_names[colIdx])
  
  # Set the appropriate Censoring Limit for contaminant
  dat.2$censorLimitVec <- dat.2[,colIdx+1]
  
  # Must tell JAGS which observations are ABOVE Censoring limit
  dat.2$isNotCensored <- (dat.2[,colIdx] > dat.2$censorLimitVec) # Returns TRUE if not below DL and NA if below DL
  
  # Convert NA in dat$isNotCensored to FALSE
  dat.2$isNotCensored <- ifelse(is.na(dat.2$isNotCensored), FALSE, dat.2$isNotCensored)
  
  #################################################################
  ########## BUGS CODE ############################################
  #################################################################
  
  ## copy of data set
  jags.dat.2 <- copy(dat.2)
  
  ## load jagsUI package ##
  library(jagsUI)
  
  ####################################################################################
  ################### MODEL (WITHOUT LANDUSE) #########################################
  ####################################################################################
  
  print('CULTIVATED IMMEDIATE BOTH')
  print(dat_col_names[colIdx])
  
  # Define the model in the BUGS language and write a text file
  sink("vary.int.vary.slope.txt")
  cat("
    model {
    
    
    # Likelihood: 
    # Level-1 of the model
    for (i in 1:n){ 
    isNotCensored[i] ~ dinterval(y[i] , censorLimitVec[i] )
    y[i] ~ dnorm(mu[i], tau)               
    mu[i] <- alpha[group[i]] + beta[group[i]] * tissue[i]         
    } 
    
    
    # Level-2 of the model
    for(j in 1:J){
    alpha[j] <- BB[j,1]
    beta[j] <- BB[j,2]
    
    BB[j,1:K] ~ dmnorm(BB.hat[j,], Tau.B[,]) 
    
    BB.hat[j,1] <- mu.alpha 
    BB.hat[j,2] <- mu.beta
    }
    
    
    # Priors and derived quantities
    sigma ~ dunif(0, 100)
    tau <- pow(sigma,-2) # precision
    sigma2 <- pow(sigma,2)
    
    mu.alpha ~ dnorm(0, 0.0001)
    mu.beta ~ dnorm(0, 0.0001)
    
    for(j in 1:J){
    yoy.site.mean[j] <- exp(alpha[j])
    ovary.site.mean[j] <- exp(alpha[j] + beta[j])
    Ovary.effect[j] <- beta[j]
    site.mean.diff[j] <- ovary.site.mean[j] - yoy.site.mean[j]
    }
    
    overall.ovary <- mu.alpha + mu.beta
    overall.yoy <- mu.alpha
    overall.ovary.effect <- mu.beta
    
    ### Model variance-covariance
    Tau.B[1:K,1:K] ~ dwish(W[,], df)
    df <- K+1
    Sigma.B[1:K,1:K] <- inverse(Tau.B[,])
    for (k in 1:K){
    for (k.prime in 1:K){
    rho.B[k,k.prime] <- Sigma.B[k,k.prime]/
    sqrt(Sigma.B[k,k]*Sigma.B[k.prime,k.prime])
    }
    sigma.B[k] <- sqrt(Sigma.B[k,k])
    }

    
    } # end model
    ",fill = TRUE)
  sink()
  
  # Number of groups (different sites)
  J <- length(unique(jags.dat.2$Site))
  
  # Site indicator
  jags.dat.2$G <- as.numeric(as.factor(as.numeric(as.factor(jags.dat.2$Site))))
  
  ## change tissue types to 1's (Ovary) and  0's (YOY)
  tissue <- as.numeric(dat.2$Tissue)
  tissue <- ifelse(tissue==1, 1,0)
  
  ## number of varying parameters
  
  K <- 2
  
  W <- diag(K)
  
  #log transform concentrations
  Y <- log(dat.2[,colIdx])
  n <- length(Y)
  
  # Load data
  data <- list(y = Y, n = n,
               group = jags.dat.2$G,
               J = J,
               tissue = tissue,
               censorLimitVec = log(jags.dat.2$censorLimitVec),
               isNotCensored = as.numeric(jags.dat.2$isNotCensored),
               K=K,
               W=W)
  
  ### Initial values ###
  yInit <- rep( NA , nrow(jags.dat.2) )
  
  ## set to half the censor limit
  yInit[is.na(Y)] <- jags.dat.2$censorLimitVec[is.na(Y)]*0.5
  
  # Log-transform initial values for censored data
  yInit <- log(yInit)
  
  head(yInit)
  
  # Initial values
  inits <- function (){
    list (mu.alpha = rnorm(1),
          mu.beta=rnorm(1),
          sigma=runif(1),
          BB=matrix(rnorm(J*K),nrow=J,ncol=K),
          Tau.B=rwish(K+1,diag(K)),y=yInit)
    
  }
  
  # Parameters monitored
  parameters <- c("alpha",
                  "beta",
                  "mu.alpha",
                  "mu.beta",
                  "BB",
                  "sigma",
                  "Sigma.B",
                  "yoy.site.mean",
                  "ovary.site.mean",
                  "Ovary.effect",
                  "BB.hat",
                  "overall.ovary",
                  "overall.yoy",
                  "overall.ovary.effect",
                  "site.mean.diff")
  
  # MCMC settings
  ni <- 800000
  nt <- 50
  nb <- 50000
  nc <- 3
  
  # Call JAGS from R 
  out <- jags(data, inits, parameters, "vary.int.vary.slope.txt", n.chains = nc, 
              n.thin = nt, n.iter = ni, n.burnin = nb)
  
  # Summarize posteriors
  print(out, dig = 3)
  
  ### check for convergence ###
  
  #traceplot(out)
  
  ### check for significance ###
  print('Quantiles:')
  print(dat_col_names[colIdx])
  
  #convergence and rhat <1.1 
  
  # See what max Rhat value is
  r_hats[rowIdx] <- max(out$summary[, c("Rhat")])
  
  ####### pulling out information for plotting and tables ########
  
  #### Save output for graphing later #######
  
  ################################
  ## extract Juvenile site mean ##
  ################################ 
  
  ## exponentiated in model --  back transform from log transformation 
  
  alpha.yoy <- as.data.frame(out$sims.list$yoy.site.mean)
  ## grab mean of each column (mean site concentration) and credible intervals
  meanSiteAve.yoy <- round(apply(alpha.yoy, 2, mean),2)
  upperCI.SiteAve.yoy <- round(apply(alpha.yoy, 2, quantile, probs=c(0.975)),2)
  lowerCI.SiteAve.yoy <- round(apply(alpha.yoy, 2, quantile, probs=c(0.025)),2)
  ## combine into df for plotting
  df.alpha.yoy <- cbind(meanSiteAve.yoy, lowerCI.SiteAve.yoy, upperCI.SiteAve.yoy)
  
 
  #############################
  ## extract Ovary site mean ##
  #############################
  
  ## exponentiated in model --  back transform from log transformation 
  
  alpha.ovary <- as.data.frame(out$sims.list$ovary.site.mean)
  ## grab mean of each column (mean site concentration) and credible intervals
  meanSiteAve.ovary <- round(apply(alpha.ovary, 2, mean),2)
  upperCI.SiteAve.ovary <- round(apply(alpha.ovary, 2, quantile, probs=c(0.975)),2)
  lowerCI.SiteAve.ovary <- round(apply(alpha.ovary, 2, quantile, probs=c(0.025)),2)
  ## combine into df for plotting
  df.alpha.ovary <- cbind(meanSiteAve.ovary,lowerCI.SiteAve.ovary, upperCI.SiteAve.ovary)
  
  Alpha.Table[rowIdx, 1:6] <- paste(round(meanSiteAve.ovary,2),"[",round(lowerCI.SiteAve.ovary,2),",", round(upperCI.SiteAve.ovary,2),"]")
  Alpha.Table[rowIdx, 7:12] <- paste(round(meanSiteAve.yoy,2),"[",round(lowerCI.SiteAve.yoy,2),",", round(upperCI.SiteAve.yoy,2),"]")
  
  Immed_effect_alpha[(rowIdx2):(rowIdx2+5), 1:3] <- df.alpha.ovary[,]
  Immed_effect_alpha[(rowIdx2):(rowIdx2+5), 4:6] <- df.alpha.yoy[,]
  
  
  #######################################################################################
  ## extract grand mean (mu.alpha + mu.beta) - ovary (mean concentration across all sites)
  #######################################################################################
  
  alpha.ovary.2 <- as.data.frame(exp(out$sims.list$overall.ovary))
  ## grab mean of each column and credible intervals
  meanSiteAve.ovary.2 <- round(apply(alpha.ovary.2, 2, mean),2) #change to median for naphthalene
  upperCI.SiteAve.ovary.2 <- round(apply(alpha.ovary.2, 2, quantile, probs=c(0.975)),2)
  lowerCI.SiteAve.ovary.2 <- round(apply(alpha.ovary.2, 2, quantile, probs=c(0.025)),2)
  ## combine into df for plotting
  df.alpha.ovary.2 <- cbind(meanSiteAve.ovary.2,lowerCI.SiteAve.ovary.2, upperCI.SiteAve.ovary.2)
  
  
  
  #######################################################################################
  ## extract grand mean (mu.alpha's) - juvenile (mean concentration across all sites)
  #######################################################################################
  
  alpha.yoy.2 <- as.data.frame(exp(out$sims.list$overall.yoy))
  ## grab mean of each column and credible intervals
  meanSiteAve.yoy.2 <- round(apply(alpha.yoy.2, 2, mean),2) #change to median for naphthalene
  upperCI.SiteAve.yoy.2 <- round(apply(alpha.yoy.2, 2, quantile, probs=c(0.975)),2)
  lowerCI.SiteAve.yoy.2 <- round(apply(alpha.yoy.2, 2, quantile, probs=c(0.025)),2)
  ## combine into df for plotting
  df.alpha.yoy.2 <- cbind(meanSiteAve.yoy.2, lowerCI.SiteAve.yoy.2, upperCI.SiteAve.yoy.2)
  
  Alpha.Table.2[rowIdx,1] <- paste(round(meanSiteAve.ovary.2,2),"[",round(lowerCI.SiteAve.ovary.2,2),",", round(upperCI.SiteAve.ovary.2,2),"]")
  Alpha.Table.2[rowIdx,2] <- paste(round(meanSiteAve.yoy.2,2),"[",round(lowerCI.SiteAve.yoy.2,2),",", round(upperCI.SiteAve.yoy.2,2),"]")
  
  Immed_effect_alpha.2[rowIdx,1:3] <- df.alpha.ovary.2[,]
  Immed_effect_alpha.2[rowIdx,4:6] <- df.alpha.yoy.2[,]
  
  
  
  #######################################################################################
  ## extract mean site difference in ovary and juvenile tissue concentrations
  #######################################################################################
  
  beta <- as.data.frame(out$sims.list$site.mean.diff)
  ## grab mean of each column and credible intervals
  meanTissue.diff <- round(apply(beta, 2, mean),0)
  upperTissue.diff <- round(apply(beta, 2, quantile, probs=c(0.975)),1)
  lowerTissue.diff <- round(apply(beta, 2, quantile, probs=c(0.025)),1)
  ## combine into df for plotting
  df.beta <- cbind(meanTissue.diff, lowerTissue.diff, upperTissue.diff)

  Immed_effect_beta[(rowIdx2):(rowIdx2+5),] <- df.beta[,]
  Beta.Table[rowIdx,] <- paste(round(meanTissue.diff,1),"[",round(lowerTissue.diff,1),",", round(upperTissue.diff,1),"]")

  #######################################################################################
  ## Extract grand mean difference (overall difference between ovary and tissue means across all sites)
  #######################################################################################
  
  beta.2 <- as.data.frame(exp((out$sims.list$overall.ovary)-(out$sims.list$overall.yoy)))
  ## grab mean of each column and credible intervals
  meanTissue.diff.2 <- round(apply(beta.2, 2, mean),2) #change to median for naphthalene
  upperTissue.diff.2 <- round(apply(beta.2, 2, quantile, probs=c(0.975)),2)
  lowerTissue.diff.2 <- round(apply(beta.2, 2, quantile, probs=c(0.025)),2)
  ## combine into df for plotting
  df.beta.2 <- cbind(meanTissue.diff.2, lowerTissue.diff.2, upperTissue.diff.2)
  
  Beta.Table.2[rowIdx,] <- paste(round(meanTissue.diff.2,1),"[",round(lowerTissue.diff.2,1),",", round(upperTissue.diff.2,1),"]")
  Immed_effect_beta.2[rowIdx,] <- df.beta.2[,]
  
  rowIdx <- (rowIdx+1)
  rowIdx2 <- (rowIdx2+6)
  

} ### END FOR-LOOP


### finalize and save table info ### 

write.csv(Alpha.Table, "Alpha.Table.csv")
write.csv(Alpha.Table.2, "Alpha.Table.2.csv")
write.csv(Beta.Table, "Beta.Table.csv")
write.csv(Beta.Table.2, "Beta.Table.2.csv")

write.csv(Immed_effect_alpha, "Immed_effect_alpha.csv")
write.csv(Immed_effect_alpha.2, "Immed_effect_alpha.2.csv")
write.csv(Immed_effect_beta, "Immed_effect_beta.csv")
write.csv(Immed_effect_beta.2, "Immed_effect_beta.2.csv")


Nap.alpha.2 <- Alpha.Table.2[5,]
Nap.alpha <- Alpha.Table[5,]
Nap.beta <-Beta.Table[5,]
Nap.beta.2 <-Beta.Table.2[5,]

write.csv(Nap.alpha.2, "Nap.alpha.2.csv")
write.csv(Nap.alpha, "Nap.alpha.csv")
write.csv(Nap.beta, "Nap.beta.csv")
write.csv(Nap.beta.2, "Nap.beta.2.csv")

Nap.Immed_effect_alpha <- Immed_effect_alpha[25:30,]
Nap.Immed_effect_alpha.2 <- Immed_effect_alpha.2[5,]
Nap.Immed_effect_beta <-Immed_effect_beta[25:30,]
Nap.Immed_effect_beta.2 <-Immed_effect_beta.2[5,]

write.csv(Nap.Immed_effect_alpha.2, "Nap.Immed_effect_alpha.2.csv")
write.csv(Nap.Immed_effect_alpha, "Nap.Immed_effect_alpha.csv")
write.csv(Nap.Immed_effect_beta, "Nap.Immed_effect_beta.csv")
write.csv(Nap.Immed_effect_beta.2, "Nap.Immed_effect_beta.2.csv")
