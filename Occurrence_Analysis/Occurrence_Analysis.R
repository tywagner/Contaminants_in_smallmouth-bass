# clear workspace
rm(list=ls())

##load libraries
library(data.table) #for fread
library(jagsUI)
library(lme4)
library(MCMCpack)
library(ggplot2)
library(ggpubr) #for ggarrange
#library(ggeffects) #for predicted probabilities


## set working directory ##
#setwd("~/Desktop/Penn State/Classwork/CEC_Project/Chapter 2/R Scripts/Occurrence Analysis/Logistic Regression - jags")

################################
###### load juvenile data ######
################################

dat.juv <- fread("dat.final.p1.juv.csv")

#str(dat.juv)

#change names of columns
colnames(dat.juv)[5] <- "Ag" ## agricultural landuse
colnames(dat.juv)[7] <- "Dev" ## developed landuse
colnames(dat.juv)[8] <- "PCT.BMP"
colnames(dat.juv)[9] <- "BMP" ## BMP Intensity
colnames(dat.juv)[14] <- "ppdde"
colnames(dat.juv)[49] <- "methylsalicylate"
colnames(dat.juv)[26] <- "TBP"

dat.juv$Site <- as.factor(dat.juv$Site)

#remove weird first column
dat.juv <- dat.juv[,c(2:52)]

#remove compounds that don't meet the 20% frequency of detection threshold
#keep:hepta, ppdde(?),galax, chlorophyros,TBP,nap,phen,pyrene,caffeine,idole
#isophorone,methylsal
dat.juv <- dat.juv[,c(1:9,12:13,15,18,25,35:37,41,45:46,48)]


##############################
###### load ovary data ###### 
##############################

dat.ovary <- fread("dat.final.p1.ovary.csv")

#change names of columns
colnames(dat.ovary)[5] <- "Ag" ## agricultural landuse
colnames(dat.ovary)[7] <- "Dev" ## developed landuse
colnames(dat.ovary)[8] <- "PCT.BMP" 
colnames(dat.ovary)[9] <- "BMP" ## BMP Intensity
colnames(dat.ovary)[17] <- "ppdde"
colnames(dat.ovary)[17] <- "opDDT"
colnames(dat.ovary)[27] <- "diethylphthalate"
colnames(dat.ovary)[28] <- "DEHP"
colnames(dat.ovary)[51] <- "methyltriclosan"
colnames(dat.ovary)[34] <- "methylnaphthalene"

# change Site to factor
dat.ovary$Site <- as.factor(dat.ovary$Site)

#remove first column
dat.ovary <- dat.ovary[,c(2:56)]

#remove compounds that don't meet the 20% frequency of detection threshold
#keep: hexa, hepta, octochlor, ppdde,opDDT,mirex, galax, diethylphal,
#DEHP,methylnap,anthracene,nap,phen, benzos, methyltric
dat.ovary <- dat.ovary[,c(1:9,12:14,16:19,26:27,33:34,36:37,41,50)]



############################
###### load "Both" Data ######
############################
#Data for analysis with all sites that had both juvenile and ovary tissue data#

dat.both <- fread("dat.final.p2.csv")

#change names of columns
colnames(dat.both)[5] <- "Ag" ## agricultural landuse
colnames(dat.both)[7] <- "Dev" ## developed landuse
colnames(dat.both)[8] <- "PCT.BMP" 
colnames(dat.both)[9] <- "BMP" ## BMP Intensity
colnames(dat.both)[14] <- "ppdde"
colnames(dat.both)[21] <- "diethylphthalate"
colnames(dat.both)[42] <- "methyltriclosan"

# change Site and Tissue to factors
dat.both$Site <- as.factor(dat.both$Site)
dat.both$Tissue <- as.factor(dat.both$Tissue)

# indicate that ovary tissue is the reference cell (i.e., ovary = 1; juvenile = 0)
dat.both$Tissue <- ifelse(dat.both$Tissue == 'Ovary', 1, 0)

# remove first column
dat.both <- dat.both[,c(2:44)]

#remove compounds that don't meet the 20% frequency of detection threshold
#keep: hepta, ppdde, galax, diethylphal,nap,phen,pyrene,isophorone,methyltric
dat.both <- dat.both[,c(1:9,12:14,20,27:29,38,41)]


####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### 

####### ####### ####### ####### 
####### BOTH Tissues Model ######
####### ####### ####### ####### 

## 9 compounds
## analysis with all sites that had both juvenile and ovary tissue data

####### ####### ####### ####### ####### ####### ####### ####### 
## matrices to fill in with estimates and credible intervals ##
####### ####### ####### ####### ####### ####### ####### ####### 

## Overall Ag ##
overall.both.ag <- matrix(ncol= 6,nrow = 9)
overall.both.ag <- as.data.frame(overall.both.ag)
rownames(overall.both.ag) <- c("Heptachlor",
                              "ppdde",
                              "Galaxolide",
                              "Diethyl phthalate",
                              "Naphthalene",
                              "Phenanthrene",
                              "Pyrene",
                              "Isophorone",
                              "Methyl triclosan")
overall.both.ag[,7] <- c("Heptachlor",
                          "ppdde",
                          "Galaxolide",
                          "Diethyl phthalate",
                          "Naphthalene",
                          "Phenanthrene",
                          "Pyrene",
                          "Isophorone",
                          "Methyl triclosan")
colnames(overall.both.ag) <- c("Ovary",
                              "Ovary.Upper",
                              "Ovary.Lower",
                              "Juvenile",
                              "Juvenile.Upper",
                              "Juvenile.Lower",
                              "Chem.Names")
## for manuscript

overall.both.ag.2 <- matrix(ncol= 2,nrow = 9)
overall.both.ag.2 <- as.data.frame(overall.both.ag.2)
rownames(overall.both.ag.2) <- c("Heptachlor",
                               "ppdde",
                               "Galaxolide",
                               "Diethyl phthalate",
                               "Naphthalene",
                               "Phenanthrene",
                               "Pyrene",
                               "Isophorone",
                               "Methyl triclosan")
colnames(overall.both.ag.2) <- c("Ovary",
                               "Juvenile")



## Overall DEV ##
overall.both.dev <- matrix(ncol= 7,nrow = 9)
overall.both.dev <- as.data.frame(overall.both.dev)
rownames(overall.both.dev) <- c("Heptachlor",
                               "ppdde",
                               "Galaxolide",
                               "Diethyl phthalate",
                               "Naphthalene",
                               "Phenanthrene",
                               "Pyrene",
                               "Isophorone",
                               "Methyl triclosan")
overall.both.dev[,7] <- c("Heptachlor",
                          "ppdde",
                          "Galaxolide",
                          "Diethyl phthalate",
                          "Naphthalene",
                          "Phenanthrene",
                          "Pyrene",
                          "Isophorone",
                          "Methyl triclosan")
colnames(overall.both.dev) <- c("Ovary",
                               "Ovary.Upper",
                               "Ovary.Lower",
                               "Juvenile",
                               "Juvenile.Upper",
                               "Juvenile.Lower",
                               "Chem.Names")

## for manuscript

overall.both.dev.2 <- matrix(ncol= 2,nrow = 9)
overall.both.dev.2 <- as.data.frame(overall.both.dev.2)
rownames(overall.both.dev.2) <- c("Heptachlor",
                                 "ppdde",
                                 "Galaxolide",
                                 "Diethyl phthalate",
                                 "Naphthalene",
                                 "Phenanthrene",
                                 "Pyrene",
                                 "Isophorone",
                                 "Methyl triclosan")
colnames(overall.both.dev.2) <- c("Ovary",
                                 "Juvenile")

## Overall BMP ##
overall.both.BMP <- matrix(ncol= 7,nrow = 9)
overall.both.BMP <- as.data.frame(overall.both.BMP)
rownames(overall.both.BMP) <- c("Heptachlor",
                                "ppdde",
                                "Galaxolide",
                                "Diethyl phthalate",
                                "Naphthalene",
                                "Phenanthrene",
                                "Pyrene",
                                "Isophorone",
                                "Methyl triclosan")
overall.both.BMP[,7] <- c("Heptachlor",
                       "ppdde",
                       "Galaxolide",
                       "Diethyl phthalate",
                       "Naphthalene",
                       "Phenanthrene",
                       "Pyrene",
                       "Isophorone",
                       "Methyl triclosan")
colnames(overall.both.BMP) <- c("Ovary",
                                "Ovary.Upper",
                                "Ovary.Lower",
                                "Juvenile",
                                "Juvenile.Upper",
                                "Juvenile.Lower",
                                "Chem.Names")

## for manuscript

overall.both.BMP.2 <- matrix(ncol= 2,nrow = 9)
overall.both.BMP.2 <- as.data.frame(overall.both.BMP.2)
rownames(overall.both.BMP.2) <- c("Heptachlor",
                                  "ppdde",
                                  "Galaxolide",
                                  "Diethyl phthalate",
                                  "Naphthalene",
                                  "Phenanthrene",
                                  "Pyrene",
                                  "Isophorone",
                                  "Methyl triclosan")
colnames(overall.both.BMP.2) <- c("Ovary",
                                  "Juvenile")

## Overall Difference in Occurrence ##
overall.difference <- matrix(ncol= 10,nrow = 9)
overall.difference <- as.data.frame(overall.difference)
rownames(overall.difference) <- c("Heptachlor",
                                "ppdde",
                                "Galaxolide",
                                "Diethyl phthalate",
                                "Naphthalene",
                                "Phenanthrene",
                                "Pyrene",
                                "Isophorone",
                                "Methyl triclosan")
overall.difference[,10] <- c("Heptachlor",
                          "ppdde",
                          "Galaxolide",
                          "Diethyl phthalate",
                          "Naphthalene",
                          "Phenanthrene",
                          "Pyrene",
                          "Isophorone",
                          "Methyl triclosan")
colnames(overall.difference) <- c("Ag.Diff",
                                "Ag.Diff.Upper",
                                "Ag.Diff.Lower",
                                "Dev.Diff",
                                "Dev.Diff.Upper",
                                "Dev.Diff.Lower",
                                "BMP.Diff",
                                "BMP.Diff.Upper",
                                "BMP.Diff.Lower",
                                "Chem.Names")

## for manuscript

overall.difference.2 <- matrix(ncol= 3,nrow = 9)
overall.difference.2 <- as.data.frame(overall.difference.2)
rownames(overall.difference.2) <- c("Heptachlor",
                                  "ppdde",
                                  "Galaxolide",
                                  "Diethyl phthalate",
                                  "Naphthalene",
                                  "Phenanthrene",
                                  "Pyrene",
                                  "Isophorone",
                                  "Methyl triclosan")
colnames(overall.difference.2) <- c("Ag.Diff","Dev.Diff","BMP.Diff")

####### ####### 
## BY SITES ####### 
####### ####### 

## Ag.Sites.YOY ##
sites.YOY.ag <- matrix(ncol= 19,nrow = 9)
sites.YOY.ag <- as.data.frame(sites.YOY.ag)
rownames(sites.YOY.ag) <- c("Heptachlor",
                               "ppdde",
                               "Galaxolide",
                               "Diethyl phthalate",
                               "Naphthalene",
                               "Phenanthrene",
                               "Pyrene",
                               "Isophorone",
                               "Methyl triclosan")
sites.YOY.ag[,19] <- c("Heptachlor",
                         "ppdde",
                         "Galaxolide",
                         "Diethyl phthalate",
                         "Naphthalene",
                         "Phenanthrene",
                         "Pyrene",
                         "Isophorone",
                         "Methyl triclosan")
colnames(sites.YOY.ag) <- c("ANT",
                            "ANT.Upper",
                            "ANT.Lower",
                            "BE",
                            "BE.Upper",
                            "BE.Lower",
                            "CHIL",
                            "CHIL.Upper",
                            "CHIL.Lower",
                            "PC",
                            "PC.Upper",
                            "PC.Lower",
                            "POT",
                            "POT.Upper",
                            "POT.Lower",
                            "WBM",
                            "WBM.Upper",
                            "WBM.Lower",
                            "Chem.Names")

##for manuscript
sites.YOY.ag.2 <- matrix(ncol= 6,nrow = 9)
sites.YOY.ag.2 <- as.data.frame(sites.YOY.ag.2)
rownames(sites.YOY.ag.2) <- c("Heptachlor",
                            "ppdde",
                            "Galaxolide",
                            "Diethyl phthalate",
                            "Naphthalene",
                            "Phenanthrene",
                            "Pyrene",
                            "Isophorone",
                            "Methyl triclosan")
colnames(sites.YOY.ag.2) <- c("ANT",
                            "BE",
                            "CHIL",
                            "PC",
                            "POT",
                            "WBM")

## Ag.Sites.OVARY ##
sites.OVARY.ag <- matrix(ncol= 19,nrow = 9)
sites.OVARY.ag <- as.data.frame(sites.OVARY.ag)
rownames(sites.OVARY.ag) <- c("Heptachlor",
                            "ppdde",
                            "Galaxolide",
                            "Diethyl phthalate",
                            "Naphthalene",
                            "Phenanthrene",
                            "Pyrene",
                            "Isophorone",
                            "Methyl triclosan")
sites.OVARY.ag[,19] <- c("Heptachlor",
                        "ppdde",
                        "Galaxolide",
                        "Diethyl phthalate",
                        "Naphthalene",
                        "Phenanthrene",
                        "Pyrene",
                        "Isophorone",
                        "Methyl triclosan")
colnames(sites.OVARY.ag) <- c("ANT",
                            "ANT.Upper",
                            "ANT.Lower",
                            "BE",
                            "BE.Upper",
                            "BE.Lower",
                            "CHIL",
                            "CHIL.Upper",
                            "CHIL.Lower",
                            "PC",
                            "PC.Upper",
                            "PC.Lower",
                            "POT",
                            "POT.Upper",
                            "POT.Lower",
                            "WBM",
                            "WBM.Upper",
                            "WBM.Lower",
                            "Chem.Names")

##for manuscript
sites.OVARY.ag.2 <- matrix(ncol= 6,nrow = 9)
sites.OVARY.ag.2 <- as.data.frame(sites.OVARY.ag.2)
rownames(sites.OVARY.ag.2) <- c("Heptachlor",
                              "ppdde",
                              "Galaxolide",
                              "Diethyl phthalate",
                              "Naphthalene",
                              "Phenanthrene",
                              "Pyrene",
                              "Isophorone",
                              "Methyl triclosan")
colnames(sites.OVARY.ag.2) <- c("ANT",
                              "BE",
                              "CHIL",
                              "PC",
                              "POT",
                              "WBM")

## Dev.Sites.YOY ##
sites.YOY.dev <- matrix(ncol= 19,nrow = 9)
sites.YOY.dev <- as.data.frame(sites.YOY.dev)
rownames(sites.YOY.dev) <- c("Heptachlor",
                            "ppdde",
                            "Galaxolide",
                            "Diethyl phthalate",
                            "Naphthalene",
                            "Phenanthrene",
                            "Pyrene",
                            "Isophorone",
                            "Methyl triclosan")
sites.YOY.dev[,19] <- c("Heptachlor",
                          "ppdde",
                          "Galaxolide",
                          "Diethyl phthalate",
                          "Naphthalene",
                          "Phenanthrene",
                          "Pyrene",
                          "Isophorone",
                          "Methyl triclosan")
colnames(sites.YOY.dev) <- c("ANT",
                            "ANT.Upper",
                            "ANT.Lower",
                            "BE",
                            "BE.Upper",
                            "BE.Lower",
                            "CHIL",
                            "CHIL.Upper",
                            "CHIL.Lower",
                            "PC",
                            "PC.Upper",
                            "PC.Lower",
                            "POT",
                            "POT.Upper",
                            "POT.Lower",
                            "WBM",
                            "WBM.Upper",
                            "WBM.Lower",
                            "Chem.Names")

##for manuscript
sites.YOY.dev.2 <- matrix(ncol= 6,nrow = 9)
sites.YOY.dev.2 <- as.data.frame(sites.YOY.dev.2)
rownames(sites.YOY.dev.2) <- c("Heptachlor",
                                "ppdde",
                                "Galaxolide",
                                "Diethyl phthalate",
                                "Naphthalene",
                                "Phenanthrene",
                                "Pyrene",
                                "Isophorone",
                                "Methyl triclosan")
colnames(sites.YOY.dev.2) <- c("ANT",
                                "BE",
                                "CHIL",
                                "PC",
                                "POT",
                                "WBM")

##  Dev.Sites.OVARY ##
sites.OVARY.dev <- matrix(ncol= 19,nrow = 9)
sites.OVARY.dev <- as.data.frame(sites.OVARY.dev)
rownames(sites.OVARY.dev) <- c("Heptachlor",
                              "ppdde",
                              "Galaxolide",
                              "Diethyl phthalate",
                              "Naphthalene",
                              "Phenanthrene",
                              "Pyrene",
                              "Isophorone",
                              "Methyl triclosan")
sites.OVARY.dev[,19] <- c("Heptachlor",
                        "ppdde",
                        "Galaxolide",
                        "Diethyl phthalate",
                        "Naphthalene",
                        "Phenanthrene",
                        "Pyrene",
                        "Isophorone",
                        "Methyl triclosan")
colnames(sites.OVARY.dev) <- c("ANT",
                              "ANT.Upper",
                              "ANT.Lower",
                              "BE",
                              "BE.Upper",
                              "BE.Lower",
                              "CHIL",
                              "CHIL.Upper",
                              "CHIL.Lower",
                              "PC",
                              "PC.Upper",
                              "PC.Lower",
                              "POT",
                              "POT.Upper",
                              "POT.Lower",
                              "WBM",
                              "WBM.Upper",
                              "WBM.Lower",
                              "Chem.Names")

##for manuscript
sites.OVARY.dev.2 <- matrix(ncol= 6,nrow = 9)
sites.OVARY.dev.2 <- as.data.frame(sites.OVARY.dev.2)
rownames(sites.OVARY.dev.2) <- c("Heptachlor",
                               "ppdde",
                               "Galaxolide",
                               "Diethyl phthalate",
                               "Naphthalene",
                               "Phenanthrene",
                               "Pyrene",
                               "Isophorone",
                               "Methyl triclosan")
colnames(sites.OVARY.dev.2) <- c("ANT",
                               "BE",
                               "CHIL",
                               "PC",
                               "POT",
                               "WBM")


##  BMP.Sites.YOY ##
sites.YOY.BMP <- matrix(ncol= 19,nrow = 9)
sites.YOY.BMP <- as.data.frame(sites.YOY.BMP)
rownames(sites.YOY.BMP) <- c("Heptachlor",
                             "ppdde",
                             "Galaxolide",
                             "Diethyl phthalate",
                             "Naphthalene",
                             "Phenanthrene",
                             "Pyrene",
                             "Isophorone",
                             "Methyl triclosan")
sites.YOY.BMP[,19] <- c("Heptachlor",
                          "ppdde",
                          "Galaxolide",
                          "Diethyl phthalate",
                          "Naphthalene",
                          "Phenanthrene",
                          "Pyrene",
                          "Isophorone",
                          "Methyl triclosan")
colnames(sites.YOY.BMP) <- c("ANT",
                             "ANT.Upper",
                             "ANT.Lower",
                             "BE",
                             "BE.Upper",
                             "BE.Lower",
                             "CHIL",
                             "CHIL.Upper",
                             "CHIL.Lower",
                             "PC",
                             "PC.Upper",
                             "PC.Lower",
                             "POT",
                             "POT.Upper",
                             "POT.Lower",
                             "WBM",
                             "WBM.Upper",
                             "WBM.Lower",
                             "Chem.Names")

##for manuscript
sites.YOY.BMP.2 <- matrix(ncol= 6,nrow = 9)
sites.YOY.BMP.2 <- as.data.frame(sites.YOY.BMP.2)
rownames(sites.YOY.BMP.2) <- c("Heptachlor",
                                 "ppdde",
                                 "Galaxolide",
                                 "Diethyl phthalate",
                                 "Naphthalene",
                                 "Phenanthrene",
                                 "Pyrene",
                                 "Isophorone",
                                 "Methyl triclosan")
colnames(sites.OVARY.dev.2) <- c("ANT",
                                 "BE",
                                 "CHIL",
                                 "PC",
                                 "POT",
                                 "WBM")


## BMP.Sites.OVARY ##
sites.OVARY.BMP <- matrix(ncol= 19,nrow = 9)
sites.OVARY.BMP <- as.data.frame(sites.OVARY.BMP)
rownames(sites.OVARY.BMP) <- c("Heptachlor",
                               "ppdde",
                               "Galaxolide",
                               "Diethyl phthalate",
                               "Naphthalene",
                               "Phenanthrene",
                               "Pyrene",
                               "Isophorone",
                               "Methyl triclosan")
sites.OVARY.BMP[,19] <- c("Heptachlor",
                           "ppdde",
                           "Galaxolide",
                           "Diethyl phthalate",
                           "Naphthalene",
                           "Phenanthrene",
                           "Pyrene",
                           "Isophorone",
                           "Methyl triclosan")
colnames(sites.OVARY.BMP) <- c("ANT",
                               "ANT.Upper",
                               "ANT.Lower",
                               "BE",
                               "BE.Upper",
                               "BE.Lower",
                               "CHIL",
                               "CHIL.Upper",
                               "CHIL.Lower",
                               "PC",
                               "PC.Upper",
                               "PC.Lower",
                               "POT",
                               "POT.Upper",
                               "POT.Lower",
                               "WBM",
                               "WBM.Upper",
                               "WBM.Lower",
                               "Chem.Names")

##for manuscript
sites.OVARY.BMP.2 <- matrix(ncol= 6,nrow = 9)
sites.OVARY.BMP.2 <- as.data.frame(sites.OVARY.BMP.2)
rownames(sites.OVARY.BMP.2) <- c("Heptachlor",
                               "ppdde",
                               "Galaxolide",
                               "Diethyl phthalate",
                               "Naphthalene",
                               "Phenanthrene",
                               "Pyrene",
                               "Isophorone",
                               "Methyl triclosan")
colnames(sites.OVARY.BMP.2) <- c("ANT",
                                 "BE",
                                 "CHIL",
                                 "PC",
                                 "POT",
                                 "WBM")


## Differences in Occurrance Sites - Ag ##
sites.diff.ag <- matrix(ncol= 19,nrow = 9)
sites.diff.ag <- as.data.frame(sites.diff.ag)
rownames(sites.diff.ag) <- c("Heptachlor",
                               "ppdde",
                               "Galaxolide",
                               "Diethyl phthalate",
                               "Naphthalene",
                               "Phenanthrene",
                               "Pyrene",
                               "Isophorone",
                               "Methyl triclosan")
sites.diff.ag[,19] <- c("Heptachlor",
                          "ppdde",
                          "Galaxolide",
                          "Diethyl phthalate",
                          "Naphthalene",
                          "Phenanthrene",
                          "Pyrene",
                          "Isophorone",
                          "Methyl triclosan")
colnames(sites.diff.ag) <- c("ANT",
                               "ANT.Upper",
                               "ANT.Lower",
                               "BE",
                               "BE.Upper",
                               "BE.Lower",
                               "CHIL",
                               "CHIL.Upper",
                               "CHIL.Lower",
                               "PC",
                               "PC.Upper",
                               "PC.Lower",
                               "POT",
                               "POT.Upper",
                               "POT.Lower",
                               "WBM",
                               "WBM.Upper",
                               "WBM.Lower",
                               "Chem.Names")

##for manuscript
sites.diff.ag.2 <- matrix(ncol= 6,nrow = 9)
sites.diff.ag.2 <- as.data.frame(sites.diff.ag.2)
rownames(sites.diff.ag.2) <- c("Heptachlor",
                               "ppdde",
                               "Galaxolide",
                               "Diethyl phthalate",
                               "Naphthalene",
                               "Phenanthrene",
                               "Pyrene",
                               "Isophorone",
                               "Methyl triclosan")
colnames(sites.diff.ag.2) <- c("ANT",
                                 "BE",
                                 "CHIL",
                                 "PC",
                                 "POT",
                                 "WBM")

## Differences in Occurrance Sites - Dev ##
sites.diff.dev <- matrix(ncol= 19,nrow = 9)
sites.diff.dev <- as.data.frame(sites.diff.dev)
rownames(sites.diff.dev) <- c("Heptachlor",
                             "ppdde",
                             "Galaxolide",
                             "Diethyl phthalate",
                             "Naphthalene",
                             "Phenanthrene",
                             "Pyrene",
                             "Isophorone",
                             "Methyl triclosan")
sites.diff.dev[,19] <- c("Heptachlor",
                        "ppdde",
                        "Galaxolide",
                        "Diethyl phthalate",
                        "Naphthalene",
                        "Phenanthrene",
                        "Pyrene",
                        "Isophorone",
                        "Methyl triclosan")
colnames(sites.diff.dev) <- c("ANT",
                             "ANT.Upper",
                             "ANT.Lower",
                             "BE",
                             "BE.Upper",
                             "BE.Lower",
                             "CHIL",
                             "CHIL.Upper",
                             "CHIL.Lower",
                             "PC",
                             "PC.Upper",
                             "PC.Lower",
                             "POT",
                             "POT.Upper",
                             "POT.Lower",
                             "WBM",
                             "WBM.Upper",
                             "WBM.Lower",
                             "Chem.Names")

##for manuscript
sites.diff.dev.2 <- matrix(ncol= 6,nrow = 9)
sites.diff.dev.2 <- as.data.frame(sites.diff.dev.2)
rownames(sites.diff.dev.2) <- c("Heptachlor",
                               "ppdde",
                               "Galaxolide",
                               "Diethyl phthalate",
                               "Naphthalene",
                               "Phenanthrene",
                               "Pyrene",
                               "Isophorone",
                               "Methyl triclosan")
colnames(sites.diff.dev.2) <- c("ANT",
                               "BE",
                               "CHIL",
                               "PC",
                               "POT",
                               "WBM")

## Differences in Occurrance Sites - BMP ##
sites.diff.BMP <- matrix(ncol= 19,nrow = 9)
sites.diff.BMP <- as.data.frame(sites.diff.BMP)
rownames(sites.diff.BMP) <- c("Heptachlor",
                              "ppdde",
                              "Galaxolide",
                              "Diethyl phthalate",
                              "Naphthalene",
                              "Phenanthrene",
                              "Pyrene",
                              "Isophorone",
                              "Methyl triclosan")
sites.diff.BMP[,19] <- c("Heptachlor",
                         "ppdde",
                         "Galaxolide",
                         "Diethyl phthalate",
                         "Naphthalene",
                         "Phenanthrene",
                         "Pyrene",
                         "Isophorone",
                         "Methyl triclosan")
colnames(sites.diff.BMP) <- c("ANT",
                              "ANT.Upper",
                              "ANT.Lower",
                              "BE",
                              "BE.Upper",
                              "BE.Lower",
                              "CHIL",
                              "CHIL.Upper",
                              "CHIL.Lower",
                              "PC",
                              "PC.Upper",
                              "PC.Lower",
                              "POT",
                              "POT.Upper",
                              "POT.Lower",
                              "WBM",
                              "WBM.Upper",
                              "WBM.Lower",
                              "Chem.Names")

##for manuscript
sites.diff.BMP.2 <- matrix(ncol= 6,nrow = 9)
sites.diff.BMP.2 <- as.data.frame(sites.diff.BMP.2)
rownames(sites.diff.BMP.2) <- c("Heptachlor",
                                "ppdde",
                                "Galaxolide",
                                "Diethyl phthalate",
                                "Naphthalene",
                                "Phenanthrene",
                                "Pyrene",
                                "Isophorone",
                                "Methyl triclosan")
colnames(sites.diff.BMP.2) <- c("ANT",
                                "BE",
                                "CHIL",
                                "PC",
                                "POT",
                                "WBM")


dat.col.names <- colnames(dat.both)
dat.both <- data.frame(dat.both)
chem.names <- c("Heptachlor",
                "ppdde",
                "Galaxolide",
                "Diethyl phthalate",
                "Naphthalene",
                "Phenanthrene",
                "Pyrene",
                "Isophorone",
                "Methyl triclosan")

#################################################################
########## BUGS CODE ############################################
#################################################################

rowIdx <- 1
for (colIdx in 10:18){ ## FOR-LOOP START: Loop through contaminants
  
  for(landcol in seq(from =5 , to =9, by =2)){ ## FOR-LOOP START: Loop through different landuse/BMP intensity
  
    ## keep track of model loops  
  print("NEW LOOP ITERATION.  WORKING ON:")
  print(dat.col.names[colIdx])
  print(dat.col.names[landcol])
  
  
# Define the model in the BUGS language and write a text file
sink("tissue.model.txt")
cat("
model {
 
 
# Likelihood:
# Level-1 of the model
for (i in 1:n){
   y[i] ~ dbin(p[i],1)  # distributional assumption
   p[i] <- exp(lp[i])/(1+exp(lp[i])) # logit link function
   lp[i] <- alpha[group[i]] + beta[group[i]] * tissue[i] # linear predictor   
  }
 
 
# Level-2 of the model
for(j in 1:J){
  alpha[j] <- BB[j,1]
  beta[j] <- BB[j,2]
  
  BB[j,1:K] ~ dmnorm(BB.hat[j,], Tau.B[,])
  
  BB.hat[j,1] <- mu.alpha + alpha.gamma*landuse[j]
  BB.hat[j,2] <- mu.beta  + beta.gamma*landuse[j]
} 

   
    
    # Priors and derived quantities
    sigma ~ dunif(0, 100)
    tau <- pow(sigma,-2) # precision
    sigma2 <- pow(sigma,2)
    
    mu.alpha ~ dnorm(0, 0.0001)
    mu.beta ~ dnorm(0, 0.0001)
    alpha.gamma ~ dnorm(0, 0.0001)
    beta.gamma~ dnorm(0, 0.0001)
    
    for(j in 1:J){
    yoy.site[j] <- ilogit(alpha[j]) # TW added
    ovary.site[j] <- ilogit(alpha[j] + beta[j]) # TW added
    Ovary.effect[j] <- ilogit(beta[j])
    site.tissue.diff.prob[j] <- (ovary.site[j] - yoy.site[j]) 
    }
    
    overall.ovary <- ilogit(mu.alpha + mu.beta)
    overall.yoy <- ilogit(mu.alpha)
    overall.ovary.effect <- ilogit(mu.beta)
    overall.tissue.diff.prob <- (overall.ovary - overall.yoy) # TW modified
    
    
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

# Number of sites
J <- length(unique(dat.both$Site))

# Site indicator
G <- as.numeric(as.factor(as.numeric(as.factor(dat.both$Site))))

## tissue
tissue <- dat.both$Tissue

## landuse
landuse <- unique(dat.both[,landcol])

# Number of parameters
K <- 2

# Create identity matrix for Wishart dist'n
#!!!!!!!Number of parameters to estimate (K)

W <- diag(K)

Y <- dat.both[,colIdx]
n <- length(Y)

# Load data
data <- list(y = Y,
             group = G,
             n = n,
             J = J,
             tissue = tissue,
             landuse=landuse,
             K=K,
             W=W )


# Initial values
inits <- function (){
  list (mu.alpha = rnorm(1),
        mu.beta=rnorm(1),
        alpha.gamma=rnorm(1),
        beta.gamma = rnorm(1),
        sigma=runif(1), #may not need this? in other code
        BB=matrix(rnorm(J*K),nrow=J,ncol=K),
        Tau.B=rwish(K+1,diag(K)) )
}


# Parameters monitored
parameters <- c("alpha",
                "beta",
                "mu.alpha",
                "mu.beta",
                "alpha.gamma",
                "BB",
                "sigma",
                "Sigma.B",
                "yoy.site",
                "ovary.site",
                "Ovary.effect",
                "BB.hat",
                "overall.ovary",
                "overall.yoy",
                "overall.ovary.effect",
                "site.tissue.diff.prob",
                "overall.tissue.diff.prob")


# MCMC settings
ni <- 100000
nt <- 1
nb <- 30000
nc <- 3



start.time = Sys.time()         # Set timer 
# Call JAGS from R 

out1 <- jags(data, inits, parameters, "tissue.model.txt", n.chains = nc, 
             n.thin = nt, n.iter = ni, n.burnin = nb)

# 
end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 
# Calculate computation time


# Summarize posteriors
print(out1, dig = 3)

#determine if rhat values are <1.1 and save output in csv's to review after
rhat.1 <- data.frame(out1$Rhat)

rhat.final <- matrix(NA, nrow = 6, ncol = 21)

for(col in 1:ncol(rhat.1)){ ## loop over columns to fill in matrices
  for(row in 1:nrow(rhat.1)){ ## loop over rows to fill in matrices
    
rhat.final[,] <- if(rhat.1[row, col] >1.1){
  paste("True")
  } else {
    paste("False")
  }
    
  }
  
}

path.rhats <- "./rhats"

if(landcol == 5){ ## ag loop
file.name = paste("rhat_ag", chem.names[rowIdx], ".csv", sep="")
write.csv(rhat.final, file.path(path.rhats,file = file.name))

} else if (landcol ==7){ ## dev loop
  file.name = paste("rhat_dev", chem.names[rowIdx], ".csv", sep="")
  write.csv(rhat.final, file.path(path.rhats,file = file.name))
  
} else{ ## BMP loop
  file.name = paste("rhat_BMP", chem.names[rowIdx], ".csv", sep="")
  write.csv(rhat.final, file.path(path.rhats,file = file.name))
}

# sum1 <- out$BUGSoutput$summary
# write.csv(sum1,'summary.csv',row.names=T)
# str(out)

###### ###### ###### ###### 
###### fill in tables ###### 
###### ###### ###### ###### 

###### OVERALL (mu alphas and mu alphas + mu betas) occurrence estimates ######

if(landcol== 5){  ## ag loop
       #ovary
       overall.both.ag[rowIdx,1] <- round(data.frame(mean(out1$sims.list$overall.ovary)),2)
       overall.both.ag[rowIdx,2] <- round(data.frame(quantile(out1$sims.list$overall.ovary, probs=c(0.975))),2)
       overall.both.ag[rowIdx,3] <- round(data.frame(quantile(out1$sims.list$overall.ovary, probs=c(0.025))),2)
       #juvenile
       overall.both.ag[rowIdx,4] <- round(data.frame(mean(out1$sims.list$overall.yoy)),2)
       overall.both.ag[rowIdx,5] <- round(data.frame(quantile(out1$sims.list$overall.yoy, probs=c(0.975))),2)
       overall.both.ag[rowIdx,6] <- round(data.frame(quantile(out1$sims.list$overall.yoy, probs=c(0.025))),2)
       }else if(landcol == 7){ ## dev loop
              #ovary
              overall.both.dev[rowIdx,1] <- round(data.frame(mean(out1$sims.list$overall.ovary)),2)
              overall.both.dev[rowIdx,2] <- round(data.frame(quantile(out1$sims.list$overall.ovary, probs=c(0.975))),2)
              overall.both.dev[rowIdx,3] <- round(data.frame(quantile(out1$sims.list$overall.ovary, probs=c(0.025))),2)
              #juvenile
              overall.both.dev[rowIdx,4] <- round(data.frame(mean(out1$sims.list$overall.yoy)),2)
              overall.both.dev[rowIdx,5] <- round(data.frame(quantile(out1$sims.list$overall.yoy, probs=c(0.975))),2)
              overall.both.dev[rowIdx,6] <- round(data.frame(quantile(out1$sims.list$overall.yoy, probs=c(0.025))),2)
              }else{ ## BMP loop
                #ovary
              overall.both.BMP[rowIdx,1] <- round(data.frame(mean(out1$sims.list$overall.ovary)),2)
              overall.both.BMP[rowIdx,2] <- round(data.frame(quantile(out1$sims.list$overall.ovary, probs=c(0.975))),2)
              overall.both.BMP[rowIdx,3] <- round(data.frame(quantile(out1$sims.list$overall.ovary, probs=c(0.025))),2)
              #juvenile
              overall.both.BMP[rowIdx,4] <- round(data.frame(mean(out1$sims.list$overall.yoy)),2)
              overall.both.BMP[rowIdx,5] <- round(data.frame(quantile(out1$sims.list$overall.yoy, probs=c(0.975))),2)
              overall.both.BMP[rowIdx,6] <- round(data.frame(quantile(out1$sims.list$overall.yoy, probs=c(0.025))),2)
              }

###### OVERALL (mu betas) difference in occurrence ######

if(landcol== 5){  ## ag loop
  #ovary
  overall.difference[rowIdx,1] <- round(data.frame(mean(out1$sims.list$overall.tissue.diff.prob)),2)
  overall.difference[rowIdx,2] <- round(data.frame(quantile(out1$sims.list$overall.tissue.diff.prob, probs=c(0.975))),2)
  overall.difference[rowIdx,3] <- round(data.frame(quantile(out1$sims.list$overall.tissue.diff.prob, probs=c(0.025))),2)
  
}else if(landcol == 7){ ## dev loop
  #ovary
  overall.difference[rowIdx,4] <- round(data.frame(mean(out1$sims.list$overall.tissue.diff.prob)),2)
  overall.difference[rowIdx,5] <- round(data.frame(quantile(out1$sims.list$overall.tissue.diff.prob, probs=c(0.975))),2)
  overall.difference[rowIdx,6] <- round(data.frame(quantile(out1$sims.list$overall.tissue.diff.prob, probs=c(0.025))),2)
  
}else{ ## BMP loop
  #ovary
  overall.difference[rowIdx,7] <- round(data.frame(mean(out1$sims.list$overall.tissue.diff.prob)),2)
  overall.difference[rowIdx,8] <- round(data.frame(quantile(out1$sims.list$overall.tissue.diff.prob, probs=c(0.975))),2)
  overall.difference[rowIdx,9] <- round(data.frame(quantile(out1$sims.list$overall.tissue.diff.prob, probs=c(0.025))),2)
  
}
###### By Site (alphas and alphas + betas) occurrence estimates ######      

## Ovary

if(landcol== 5){ ## ag loop
  #ag
  sites.OVARY.ag[rowIdx,c(1,4,7,10,13,16)] <- t(round(data.frame(apply(out1$sims.list$ovary.site, 2, mean)),2))
  sites.OVARY.ag[rowIdx,c(2,5,8,11,14,17)] <- t(round(data.frame(apply(out1$sims.list$ovary.site, 2, quantile, probs=c(0.975))),2))
  sites.OVARY.ag[rowIdx,c(3,6,9,12,15,18)] <- t(round(data.frame(apply(out1$sims.list$ovary.site, 2, quantile, probs=c(0.025))),2))
}else if(landcol== 7){ ## dev loop
#dev
sites.OVARY.dev[rowIdx,c(1,4,7,10,13,16)] <- t(round(data.frame(apply(out1$sims.list$ovary.site, 2, mean)),2))
sites.OVARY.dev[rowIdx,c(2,5,8,11,14,17)] <- t(round(data.frame(apply(out1$sims.list$ovary.site, 2, quantile, probs=c(0.975))),2))
sites.OVARY.dev[rowIdx,c(3,6,9,12,15,18)] <- t(round(data.frame(apply(out1$sims.list$ovary.site, 2, quantile, probs=c(0.025))),2))
}else{ ## BMP loop
#BMP
sites.OVARY.BMP[rowIdx,c(1,4,7,10,13,16)] <- t(round(data.frame(apply(out1$sims.list$ovary.site, 2, mean)),2))
sites.OVARY.BMP[rowIdx,c(2,5,8,11,14,17)] <- t(round(data.frame(apply(out1$sims.list$ovary.site, 2, quantile, probs=c(0.975))),2))
sites.OVARY.BMP[rowIdx,c(3,6,9,12,15,18)] <- t(round(data.frame(apply(out1$sims.list$ovary.site, 2, quantile, probs=c(0.025))),2))
}

## Juvenile

if(landcol== 5){ ## ag loop
  #ag
  sites.YOY.ag[rowIdx,c(1,4,7,10,13,16)] <- t(round(data.frame(apply(out1$sims.list$yoy.site, 2, mean)),2))
  sites.YOY.ag[rowIdx,c(2,5,8,11,14,17)] <- t(round(data.frame(apply(out1$sims.list$yoy.site, 2, quantile, probs=c(0.975))),2))
  sites.YOY.ag[rowIdx,c(3,6,9,12,15,18)] <- t(round(data.frame(apply(out1$sims.list$yoy.site, 2, quantile, probs=c(0.025))),2))
}else if(landcol== 7){ ## dev loop
  #dev
  sites.YOY.dev[rowIdx,c(1,4,7,10,13,16)] <- t(round(data.frame(apply(out1$sims.list$yoy.site, 2, mean)),2))
  sites.YOY.dev[rowIdx,c(2,5,8,11,14,17)] <- t(round(data.frame(apply(out1$sims.list$yoy.site, 2, quantile, probs=c(0.975))),2))
  sites.YOY.dev[rowIdx,c(3,6,9,12,15,18)] <- t(round(data.frame(apply(out1$sims.list$yoy.site, 2, quantile, probs=c(0.025))),2))
}else{ ## BMP loop
  #BMP
  sites.YOY.BMP[rowIdx,c(1,4,7,10,13,16)] <- t(round(data.frame(apply(out1$sims.list$yoy.site, 2, mean)),2))
  sites.YOY.BMP[rowIdx,c(2,5,8,11,14,17)] <- t(round(data.frame(apply(out1$sims.list$yoy.site, 2, quantile, probs=c(0.975))),2))
  sites.YOY.BMP[rowIdx,c(3,6,9,12,15,18)] <- t(round(data.frame(apply(out1$sims.list$yoy.site, 2, quantile, probs=c(0.025))),2))
}  


###### By Site (betas) differences in occurrence ######      

## Differences

if(landcol== 5){
  #ag
  sites.diff.ag[rowIdx,c(1,4,7,10,13,16)] <- t(round(data.frame(apply(out1$sims.list$site.tissue.diff.prob, 2, mean)),2))
  sites.diff.ag[rowIdx,c(2,5,8,11,14,17)] <- t(round(data.frame(apply(out1$sims.list$site.tissue.diff.prob, 2, quantile, probs=c(0.975))),2))
  sites.diff.ag[rowIdx,c(3,6,9,12,15,18)] <- t(round(data.frame(apply(out1$sims.list$site.tissue.diff.prob, 2, quantile, probs=c(0.025))),2))
}else if(landcol== 7){
  #dev
  sites.diff.dev[rowIdx,c(1,4,7,10,13,16)] <- t(round(data.frame(apply(out1$sims.list$site.tissue.diff.prob, 2, mean)),2))
  sites.diff.dev[rowIdx,c(2,5,8,11,14,17)] <- t(round(data.frame(apply(out1$sims.list$site.tissue.diff.prob, 2, quantile, probs=c(0.975))),2))
  sites.diff.dev[rowIdx,c(3,6,9,12,15,18)] <- t(round(data.frame(apply(out1$sims.list$site.tissue.diff.prob, 2, quantile, probs=c(0.025))),2))
}else{
  #BMP
  sites.diff.BMP[rowIdx,c(1,4,7,10,13,16)] <- t(round(data.frame(apply(out1$sims.list$site.tissue.diff.prob, 2, mean)),2))
  sites.diff.BMP[rowIdx,c(2,5,8,11,14,17)] <- t(round(data.frame(apply(out1$sims.list$site.tissue.diff.prob, 2, quantile, probs=c(0.975))),2))
  sites.diff.BMP[rowIdx,c(3,6,9,12,15,18)] <- t(round(data.frame(apply(out1$sims.list$site.tissue.diff.prob, 2, quantile, probs=c(0.025))),2))
}


} ## END OF FOR-LOOP: Landuse
  rowIdx <- rowIdx+1
} ## END OF FOR-LOOP: Contaminant

##### ##### ##### ##### ##### ##### ##### 
##### END OF MODEL AND TABLE FOR LOOPS ##### 
##### ##### ##### ##### ##### ##### ##### 

## combine columns for manuscript 

#### overall occurrence estimates ####

seq.mean <- c(1,4)
seq.upper <- c(2,5)
seq.lower <- c(3,6)

for(rowIdx in 1:9){ ## loop over rows to fill in table
  for(column in 1:2){ ## loop over columns to fill in table
    #ag
    overall.both.ag.2[rowIdx,column] <- paste(overall.both.ag[rowIdx,seq.mean[column]],"[",overall.both.ag[rowIdx,seq.lower[column]],",", overall.both.ag[rowIdx,seq.upper[column]],"]")
    #dev
    overall.both.dev.2[rowIdx,column] <- paste(overall.both.dev[rowIdx,seq.mean[column]],"[",overall.both.dev[rowIdx,seq.lower[column]],",", overall.both.dev[rowIdx,seq.upper[column]],"]")
    #BMP
    overall.both.BMP.2[rowIdx,column] <- paste(overall.both.BMP[rowIdx,seq.mean[column]],"[",overall.both.BMP[rowIdx,seq.lower[column]],",", overall.both.BMP[rowIdx,seq.upper[column]],"]")
    
  }
}

#### overall differences in occurrence ####

seq.mean <- c(1,4,7)
seq.upper <- c(2,5,8)
seq.lower <- c(3,6,9)

for(rowIdx in 1:9){
  for(column in 1:3){
    overall.difference.2[rowIdx,column] <- paste(overall.difference[rowIdx,seq.mean[column]],"[",overall.difference[rowIdx,seq.lower[column]],",", overall.difference[rowIdx,seq.upper[column]],"]")
  }
}

#### occurrence estimates by site ####
seq.mean <- c(1,4,7,10,13,16)
seq.upper <- c(2,5,8,11,14,17)
seq.lower <- c(3,6,9,12,15,18)

## ovary
for(rowIdx in 1:9){ ## loop over rows to fill in table
  for(column in 1:6){ ## loop over columns to fill in table
    #ag
    sites.OVARY.ag.2[rowIdx,column] <- paste(sites.OVARY.ag[rowIdx,seq.mean[column]],"[",sites.OVARY.ag[rowIdx,seq.lower[column]],",", sites.OVARY.ag[rowIdx,seq.upper[column]],"]")
    #dev
    sites.OVARY.dev.2[rowIdx,column] <- paste(sites.OVARY.dev[rowIdx,seq.mean[column]],"[",sites.OVARY.dev[rowIdx,seq.lower[column]],",", sites.OVARY.dev[rowIdx,seq.upper[column]],"]")
    #BMP
    sites.OVARY.BMP.2[rowIdx,column] <- paste(sites.OVARY.BMP[rowIdx,seq.mean[column]],"[",sites.OVARY.BMP[rowIdx,seq.lower[column]],",", sites.OVARY.BMP[rowIdx,seq.upper[column]],"]")
    
  }
}

## juveniles
for(rowIdx in 1:9){ ## loop over rows to fill in table
  for(column in 1:6){ ## loop over columns to fill in table
    #ag
    sites.YOY.ag.2[rowIdx,column] <- paste(sites.YOY.ag[rowIdx,seq.mean[column]],"[",sites.YOY.ag[rowIdx,seq.lower[column]],",", sites.YOY.ag[rowIdx,seq.upper[column]],"]")
    #dev
    sites.YOY.dev.2[rowIdx,column] <- paste(sites.YOY.dev[rowIdx,seq.mean[column]],"[",sites.YOY.dev[rowIdx,seq.lower[column]],",", sites.YOY.dev[rowIdx,seq.upper[column]],"]")
    #BMP
    sites.YOY.BMP.2[rowIdx,column] <- paste(sites.YOY.BMP[rowIdx,seq.mean[column]],"[",sites.YOY.BMP[rowIdx,seq.lower[column]],",", sites.YOY.BMP[rowIdx,seq.upper[column]],"]")
 
    }
}

#### by site differences in occurrence ####
seq.mean <- c(1,4,7,10,13,16)
seq.upper <- c(2,5,8,11,14,17)
seq.lower <- c(3,6,9,12,15,18)

for(rowIdx in 1:9){
  for(column in 1:6){
    #ag
    sites.diff.ag.2[rowIdx,column] <- paste(sites.diff.ag[rowIdx,seq.mean[column]],"[",sites.diff.ag[rowIdx,seq.lower[column]],",", sites.diff.ag[rowIdx,seq.upper[column]],"]")
    #dev
    sites.diff.dev.2[rowIdx,column] <- paste(sites.diff.dev[rowIdx,seq.mean[column]],"[",sites.diff.dev[rowIdx,seq.lower[column]],",", sites.diff.dev[rowIdx,seq.upper[column]],"]")
    #BMP
    sites.diff.BMP.2[rowIdx,column] <- paste(sites.diff.BMP[rowIdx,seq.mean[column]],"[",sites.diff.BMP[rowIdx,seq.lower[column]],",", sites.diff.BMP[rowIdx,seq.upper[column]],"]")
    
  }
}



## Save all tables ##

path.overall.tables <- "./both.overall.tables"

## overall 
write.csv(overall.both.ag,file.path(path.overall.tables,"overall.both.ag.csv"))
write.csv(overall.both.dev,file.path(path.overall.tables,"overall.both.dev.csv"))
write.csv(overall.both.BMP,file.path(path.overall.tables,"overall.both.BMP.csv"))

write.csv(overall.both.ag.2,file.path(path.overall.tables,"overall.both.ag.2.csv"))
write.csv(overall.both.dev.2,file.path(path.overall.tables,"overall.both.dev.2.csv"))
write.csv(overall.both.BMP.2,file.path(path.overall.tables,"overall.both.BMP.2.csv"))

write.csv(overall.difference,file.path(path.overall.tables,"overall.difference.csv"))
write.csv(overall.difference.2,file.path(path.overall.tables,"overall.difference.2.csv"))

## by site

path.sites.tables <- "./both.sites.tables"

write.csv(sites.OVARY.ag,file.path(path.sites.tables, "sites.OVARY.ag.csv"))
write.csv(sites.OVARY.dev,file.path(path.sites.tables, "sites.OVARY.dev.csv"))
write.csv(sites.OVARY.BMP,file.path(path.sites.tables, "sites.OVARY.BMP.csv"))

write.csv(sites.OVARY.ag.2,file.path(path.sites.tables, "sites.OVARY.ag.2.csv"))
write.csv(sites.OVARY.dev.2,file.path(path.sites.tables, "sites.OVARY.dev.2.csv"))
write.csv(sites.OVARY.BMP.2,file.path(path.sites.tables, "sites.OVARY.BMP.2.csv"))

write.csv(sites.YOY.ag,file.path(path.sites.tables, "sites.YOY.ag.csv"))
write.csv(sites.YOY.dev,file.path(path.sites.tables, "sites.YOY.dev.csv"))
write.csv(sites.YOY.BMP,file.path(path.sites.tables, "sites.YOY.BMP.csv"))

write.csv(sites.YOY.ag.2,file.path(path.sites.tables, "sites.YOY.ag.2.csv"))
write.csv(sites.YOY.dev.2,file.path(path.sites.tables, "sites.YOY.dev.2.csv"))
write.csv(sites.YOY.BMP.2,file.path(path.sites.tables, "sites.YOY.BMP.2.csv"))

write.csv(sites.diff.ag,file.path(path.sites.tables, "sites.diff.ag.csv"))
write.csv(sites.diff.dev,file.path(path.sites.tables, "sites.diff.dev.csv"))
write.csv(sites.diff.BMP,file.path(path.sites.tables, "sites.diff.BMP.csv"))

write.csv(sites.diff.ag.2,file.path(path.sites.tables, "sites.diff.ag.2.csv"))
write.csv(sites.diff.dev.2,file.path(path.sites.tables, "sites.diff.dev.2.csv"))
write.csv(sites.diff.BMP.2,file.path(path.sites.tables, "sites.diff.BMP.2.csv"))


##### ##### ##### ###
##### Plotting ##### 
##### ##### ##### ##


##### Ag #####

## create long data set to use for plotting
sites.ag.combo <- as.data.frame(matrix(NA,nrow = 108, ncol=6))
colnames(sites.ag.combo) <- c("Chem","Tissue","Site","Mean","Up", "Low")
sites.ag.combo[,1] <- c(rep("Heptachlor",12),
                        rep("ppdde",12),
                        rep("Galaxolide",12),
                        rep("Diethyl phthalate",12),
                        rep("Naphthalene",12),
                        rep("Phenanthrene",12),
                        rep("Pyrene",12),
                        rep("Isophorone",12),
                        rep("Methyl triclosan",12))
sites.ag.combo[,2] <- rep(c("Ovary","Juvenile"),54)
sites.ag.combo[,3] <- rep(c(rep("Antietam",2),
                           rep("Bald Eagle",2),
                           rep("Chillisquaque",2),
                           rep("Pine Creek",2),
                           rep("Potomac",2),
                           rep("West Branch Mahantango",2)),9)

row.num.ovary <- 1
row.num.juv <- 2

for(j in 1:9){
for(i in seq(from=1, to=18, by=3)){
  sites.ag.combo[row.num.ovary,4:6] <- sites.OVARY.ag[j,c(i, i+1, i+2)]
  sites.ag.combo[row.num.juv,4:6] <- sites.YOY.ag[j,c(i, i+1, i+2)]
  row.num.ovary <- row.num.ovary + 2
  row.num.juv <- row.num.juv + 2
}
}


seq.1 <- matrix(NA, nrow =9, ncol = 12)
seq.1[1,] <- c(1:12)
seq.1[2,] <- c(13:24)
seq.1[3,] <- c(25:36)
seq.1[4,] <- c(37:48)
seq.1[5,] <- c(49:60)
seq.1[6,] <- c(61:72)
seq.1[7,] <- c(73:84)
seq.1[8,] <- c(85:96)
seq.1[9,] <- c(97:108)

chem.names <- c("Heptachlor",
                "ppdde",
                "Galaxolide",
                "Diethyl phthalate",
                "Naphthalene",
                "Phenanthrene",
                "Pyrene",
                "Isophorone",
                "Methyl triclosan")

for(rows in 1:9){
  chem <- chem.names[rows]
  seq.2 <- seq.1[rows,]
  
p <- ggplot(sites.ag.combo[seq.2,], aes(x= Tissue, y= Mean))+
  facet_wrap(~Site,nrow = 2, scales = 'free',strip.position = "top")+
  geom_point(size=1.8)+
  geom_errorbar(aes(x = Tissue, ymin = Low, ymax= Up), width = 0.1, size = 0.8)+
  ylim(0,1)+
  theme_classic() +
  ylab("Predicted Probability of Occurence")+
  xlab("Tissue Type")+
  theme(axis.text.y = element_text(size = 12,
                                   hjust = 0.5,
                                   vjust = 1),
        plot.margin = unit(c(2,1,2,2), "lines"),
        axis.text.x = element_text(angle = 45,
                                   vjust = 1,
                                   hjust = 1,
                                   size = 12),
        strip.background = element_blank(),
        strip.text = element_text(size = 12),
        axis.title.x = element_text(size = 16,
                                    hjust = 0.5,
                                    vjust = -3),
        axis.title.y = element_text(size = 16,
                                    hjust = 0.5,
                                    vjust = 5,
                                    angle = 90),
        strip.placement = "outside",
        panel.spacing = unit(2, "lines"))

file.name = paste("p.ag.occurrence_", chem, ".jpeg", sep="")
ggsave(filename = file.name, p, width = 22, height = 18, units = "cm", dpi = 1000)

}


##### Dev #####

## create long data set to use for plotting
sites.dev.combo <- as.data.frame(matrix(NA,nrow = 108, ncol=6))
colnames(sites.dev.combo) <- c("Chem","Tissue","Site","Mean","Up", "Low")
sites.dev.combo[,1] <- c(rep("Heptachlor",12),
                        rep("ppdde",12),
                        rep("Galaxolide",12),
                        rep("Diethyl phthalate",12),
                        rep("Naphthalene",12),
                        rep("Phenanthrene",12),
                        rep("Pyrene",12),
                        rep("Isophorone",12),
                        rep("Methyl triclosan",12))
sites.dev.combo[,2] <- rep(c("Ovary","Juvenile"),54)
sites.dev.combo[,3] <- rep(c(rep("Antietam",2),
                            rep("Bald Eagle",2),
                            rep("Chillisquaque",2),
                            rep("Pine Creek",2),
                            rep("Potomac",2),
                            rep("West Branch Mahantango",2)),9)

row.num.ovary <- 1
row.num.juv <- 2

for(j in 1:9){
  for(i in seq(from=1, to=18, by=3)){
    sites.dev.combo[row.num.ovary,4:6] <- sites.OVARY.dev[j,c(i, i+1, i+2)]
    sites.dev.combo[row.num.juv,4:6] <- sites.YOY.dev[j,c(i, i+1, i+2)]
    row.num.ovary <- row.num.ovary + 2
    row.num.juv <- row.num.juv + 2
  }
}


seq.1 <- matrix(NA, nrow =9, ncol = 12)
seq.1[1,] <- c(1:12)
seq.1[2,] <- c(13:24)
seq.1[3,] <- c(25:36)
seq.1[4,] <- c(37:48)
seq.1[5,] <- c(49:60)
seq.1[6,] <- c(61:72)
seq.1[7,] <- c(73:84)
seq.1[8,] <- c(85:96)
seq.1[9,] <- c(97:108)

chem.names <- c("Heptachlor",
                "ppdde",
                "Galaxolide",
                "Diethyl phthalate",
                "Naphthalene",
                "Phenanthrene",
                "Pyrene",
                "Isophorone",
                "Methyl triclosan")

for(rows in 1:9){
  chem <- chem.names[rows]
  seq.2 <- seq.1[rows,]
  
  p <- ggplot(sites.dev.combo[seq.2,], aes(x= Tissue, y= Mean))+
    facet_wrap(~Site,nrow = 2, scales = 'free',strip.position = "top")+
    geom_point(size=1.8)+
    geom_errorbar(aes(x = Tissue, ymin = Low, ymax= Up), width = 0.1, size = 0.8)+
    ylim(0,1)+
    theme_classic() +
    ylab("Predicted Probability of Occurence")+
    xlab("Tissue Type")+
    theme(axis.text.y = element_text(size = 12,
                                     hjust = 0.5,
                                     vjust = 1),
          plot.margin = unit(c(2,1,2,2), "lines"),
          axis.text.x = element_text(angle = 45,
                                     vjust = 1,
                                     hjust = 1,
                                     size = 12),
          strip.background = element_blank(),
          strip.text = element_text(size = 12),
          axis.title.x = element_text(size = 16,
                                      hjust = 0.5,
                                      vjust = -3),
          axis.title.y = element_text(size = 16,
                                      hjust = 0.5,
                                      vjust = 5,
                                      angle = 90),
          strip.placement = "outside",
          panel.spacing = unit(2, "lines"))
  
  file.name = paste("p.dev.occurrence_", chem, ".jpeg", sep="")
  ggsave(filename = file.name, p, width = 22, height = 18, units = "cm", dpi = 1000)
  
}


##### BMP #####

## create long data set to use for plotting
sites.BMP.combo <- as.data.frame(matrix(NA,nrow = 108, ncol=6))
colnames(sites.BMP.combo) <- c("Chem","Tissue","Site","Mean","Up", "Low")
sites.BMP.combo[,1] <- c(rep("Heptachlor",12),
                        rep("ppdde",12),
                        rep("Galaxolide",12),
                        rep("Diethyl phthalate",12),
                        rep("Naphthalene",12),
                        rep("Phenanthrene",12),
                        rep("Pyrene",12),
                        rep("Isophorone",12),
                        rep("Methyl triclosan",12))
sites.BMP.combo[,2] <- rep(c("Ovary","Juvenile"),54)
sites.BMP.combo[,3] <- rep(c(rep("Antietam",2),
                            rep("Bald Eagle",2),
                            rep("Chillisquaque",2),
                            rep("Pine Creek",2),
                            rep("Potomac",2),
                            rep("West Branch Mahantango",2)),9)

row.num.ovary <- 1
row.num.juv <- 2

for(j in 1:9){
  for(i in seq(from=1, to=18, by=3)){
    sites.BMP.combo[row.num.ovary,4:6] <- sites.OVARY.BMP[j,c(i, i+1, i+2)]
    sites.BMP.combo[row.num.juv,4:6] <- sites.YOY.BMP[j,c(i, i+1, i+2)]
    row.num.ovary <- row.num.ovary + 2
    row.num.juv <- row.num.juv + 2
  }
}


seq.1 <- matrix(NA, nrow =9, ncol = 12)
seq.1[1,] <- c(1:12)
seq.1[2,] <- c(13:24)
seq.1[3,] <- c(25:36)
seq.1[4,] <- c(37:48)
seq.1[5,] <- c(49:60)
seq.1[6,] <- c(61:72)
seq.1[7,] <- c(73:84)
seq.1[8,] <- c(85:96)
seq.1[9,] <- c(97:108)

chem.names <- c("Heptachlor",
                "ppdde",
                "Galaxolide",
                "Diethyl phthalate",
                "Naphthalene",
                "Phenanthrene",
                "Pyrene",
                "Isophorone",
                "Methyl triclosan")

for(rows in 1:9){
  chem <- chem.names[rows]
  seq.2 <- seq.1[rows,]
  
  p <- ggplot(sites.BMP.combo[seq.2,], aes(x= Tissue, y= Mean))+
    facet_wrap(~Site,nrow = 2, scales = 'free',strip.position = "top")+
    geom_point(size=1.8)+
    geom_errorbar(aes(x = Tissue, ymin = Low, ymax= Up), width = 0.1, size = 0.8)+
    ylim(0,1)+
    theme_classic() +
    ylab("Predicted Probability of Occurence")+
    xlab("Tissue Type")+
    theme(axis.text.y = element_text(size = 12,
                                     hjust = 0.5,
                                     vjust = 1),
          plot.margin = unit(c(2,1,2,2), "lines"),
          axis.text.x = element_text(angle = 45,
                                     vjust = 1,
                                     hjust = 1,
                                     size = 12),
          strip.background = element_blank(),
          strip.text = element_text(size = 12),
          axis.title.x = element_text(size = 16,
                                      hjust = 0.5,
                                      vjust = -3),
          axis.title.y = element_text(size = 16,
                                      hjust = 0.5,
                                      vjust = 5,
                                      angle = 90),
          strip.placement = "outside",
          panel.spacing = unit(2, "lines"))
  
  file.name = paste("p.BMP.occurrence_", chem, ".jpeg", sep="")
  ggsave(filename = file.name, p, width = 22, height = 18, units = "cm", dpi = 1000)
  
}






####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### 

####### ####### ####### #######  ###
####### Juvenile Tissue Model ######
####### ####### ####### ####### ####

####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### 

## 12 compounds 

##################################################
### Matrices to fill with estimates from model ###
##################################################

## Overall ##
overall.occurrence.juv <- matrix(ncol=9, nrow = 12)
overall.occurrence.juv <- as.data.frame(overall.occurrence.juv)
rownames(overall.occurrence.juv) <- c("Heptachlor",
                                      "ppdde",
                                      "Galaxolide",
                                      "Chlorpyrifos",
                                      "TBP",
                                      "Naphthalene",
                                      "Phenanthrene",
                                      "Pyrene",
                                      "Caffeine",
                                      "Indole",
                                      "Isophorone",
                                      "Methyl salicylate")
colnames(overall.occurrence.juv) <- c("Ag",
                                        "Ag.UP",
                                        "Ag.LOW",
                                        "Dev",
                                        "Dev.UP",
                                        "Dev.LOW",
                                        "BMP",
                                        "BMP.UP",
                                        "BMP.LOW")
## for table in manuscript ##
overall.occurrence.juv.2 <- matrix(ncol=3, nrow = 12)
overall.occurrence.juv.2 <- as.data.frame(overall.occurrence.juv.2)
rownames(overall.occurrence.juv.2) <- c("Heptachlor",
                                        "ppdde",
                                        "Galaxolide",
                                        "Chlorpyrifos",
                                        "TBP",
                                        "Naphthalene",
                                        "Phenanthrene",
                                        "Pyrene",
                                        "Caffeine",
                                        "Indole",
                                        "Isophorone",
                                        "Methyl salicylate")
colnames(overall.occurrence.juv.2) <- c("Ag",
                                          "Dev",
                                          "BMP")
################
## By sites ##
################

## overall ##
site.occurrence.juv <- array(NA, c(12,9,7)) #contaminants, land-use (means, CIs), sites
colnames(site.occurrence.juv) <- c("Ag",
                                     "Ag.UP",
                                     "Ag.LOW",
                                     "Dev",
                                     "Dev.UP",
                                     "Dev.LOW",
                                     "BMP",
                                     "BMP.UP",
                                     "BMP.LOW")
rownames(site.occurrence.juv) <- c("Heptachlor",
                                   "ppdde",
                                   "Galaxolide",
                                   "Chlorpyrifos",
                                   "TBP",
                                   "Naphthalene",
                                   "Phenanthrene",
                                   "Pyrene",
                                   "Caffeine",
                                   "Indole",
                                   "Isophorone",
                                   "Methyl salicylate")

## for table in manuscript ##
site.occurrence.juv.2 <- array(NA, c(12,3,7)) #contaminants, land-use (means, CIs), sites
colnames(site.occurrence.juv.2) <- c("Ag",
                                       "Dev",
                                       "BMP")
rownames(site.occurrence.juv.2) <- c("Heptachlor",
                                     "ppdde",
                                     "Galaxolide",
                                     "Chlorpyrifos",
                                     "TBP",
                                     "Naphthalene",
                                     "Phenanthrene",
                                     "Pyrene",
                                     "Caffeine",
                                     "Indole",
                                     "Isophorone",
                                     "Methyl salicylate")
####################
## land-use effects
####################

## overall ##

landuse.effect.occurrence.juv <- matrix(ncol=9, nrow = 12)
landuse.effect.occurrence.juv <- as.data.frame(landuse.effect.occurrence.juv)
rownames(landuse.effect.occurrence.juv) <- c("Heptachlor",
                                             "ppdde",
                                             "Galaxolide",
                                             "Chlorpyrifos",
                                             "TBP",
                                             "Naphthalene",
                                             "Phenanthrene",
                                             "Pyrene",
                                             "Caffeine",
                                             "Indole",
                                             "Isophorone",
                                             "Methyl salicylate")
colnames(landuse.effect.occurrence.juv) <- c("Ag",
                                               "Ag.UP",
                                               "Ag.LOW",
                                               "Dev",
                                               "Dev.UP",
                                               "Dev.LOW",
                                               "BMP",
                                               "BMP.UP",
                                               "BMP.LOW")

## for table in manuscript ##

landuse.effect.occurrence.juv.2 <- matrix(ncol=3, nrow = 12)
landuse.effect.occurrence.juv.2 <- as.data.frame(landuse.effect.occurrence.juv.2)
rownames(landuse.effect.occurrence.juv.2) <- c("Heptachlor",
                                               "ppdde",
                                               "Galaxolide",
                                               "Chlorpyrifos",
                                               "TBP",
                                               "Naphthalene",
                                               "Phenanthrene",
                                               "Pyrene",
                                               "Caffeine",
                                               "Indole",
                                               "Isophorone",
                                               "Methyl salicylate")
colnames(landuse.effect.occurrence.juv.2) <- c("Ag",
                                                 "Dev",
                                                 "BMP")




## Matrices to fill in for plotting panel plots ##

#lines
probPredPopAve.Ag <- array(NA, c(210000,100,12) )
probPredPopAve.Dev <- array(NA, c(210000,100,12) )
probPredPopAve.BMP <- array(NA, c(210000,100,12) )

#means
meanProbPopAve.Ag <- matrix(NA, nrow = 12, ncol = 100)
meanProbPopAve.Dev <- matrix(NA, nrow = 12, ncol = 100)
meanProbPopAve.BMP <- matrix(NA, nrow = 12, ncol = 100)

#upper
upperCI.PopAve.Ag <- matrix(NA, nrow = 12, ncol = 100)
upperCI.PopAve.Dev <- matrix(NA, nrow = 12, ncol = 100)
upperCI.PopAve.BMP <- matrix(NA, nrow = 12, ncol = 100)

#lower
lowerCIA.PopAve.Ag <- matrix(NA, nrow = 12, ncol = 100)
lowerCIA.PopAve.Dev <- matrix(NA, nrow = 12, ncol = 100)
lowerCIA.PopAve.BMP <- matrix(NA, nrow = 12, ncol = 100)


dat.col.names <- colnames(dat.juv)
dat.juv <- data.frame(dat.juv)
chem.names <- dat.col.names[c(10:21)]
compIdx <- 1
rowIdx <- 1

for(colIdx in 10:21){  ## FOR-LOOP START: Loop through contaminants
  for(landcol in seq(from = 4, to=8, by=2)){ ## FOR-LOOP START: Loop through landuse types
    
    # print information to keep track of model progress
    print("NEW LOOP ITERATION.  WORKING ON:")
    print(dat.col.names[colIdx])
    print(dat.col.names[landcol])
    
# Define the model in the BUGS language and write a text file
sink("model.juv.txt")
cat("
model {
 
 
# Likelihood:
# Level-1 of the model
for (i in 1:n){
   y[i] ~ dbin(p[i],1)  # distributional assumption
   p[i] <- exp(lp[i])/(1+exp(lp[i])) # logit link function
   lp[i] <- alpha[group[i]]  # linear predictor   
  }

# Level-2 of the model
    for(j in 1:J){
    alpha[j] ~ dnorm(mu.alpha[j],tau.alpha)
    mu.alpha[j] <- alpha.gamma + beta.gamma *landuse[j]
    }
    
    
    # Priors and derived quantities
    sigma ~ dunif(0, 100)
    alpha.gamma ~ dnorm(0, 0.0001)
    beta.gamma ~ dnorm(0,0.0001)
    sigma.alpha ~ dunif(0,100)
   
   overall.mean <- alpha.gamma
   site.mean <- alpha
   landuse.effect <- beta.gamma
   
   
   # Derived quantities
    tau <- pow(sigma,-2) # precision

    tau.alpha <- pow(sigma.alpha,-2) # precision
    
    } # end model
    ",fill = TRUE)
sink()

# Number of sites
J <- length(unique(dat.juv$Site))

# Site indicator
G <- as.numeric(as.factor(as.numeric(as.factor(dat.juv$Site))))

#predictor scaled and log transformed for model fit
landuse <- unique(dat.juv[,landcol])

#log transform concentrations
Y <- dat.juv[,colIdx]
n <- length(Y)

# Load data
data <- list(y = Y,
             n = n,
             landuse=landuse,
             group = G,
             J = J )

# Initial values
inits <- function (){
  list (alpha.gamma = rnorm(1),
        sigma=runif(1),
        beta.gamma=rnorm(1),
        sigma.alpha=runif(1))
}

# Parameters monitored
parameters <- c("alpha.gamma",
                "sigma",
                "beta.gamma",
                "sigma.alpha",
                "mu.alpha",
                "alpha",
                "overall.mean",
                "site.mean",
                "landuse.effect")

# MCMC settings
ni <- 100000
nt <- 1
nb <- 30000
nc <- 3

start.time = Sys.time()         # Set timer 
# Call JAGS from R 

out1 <- jags(data, inits, parameters, "model.juv.txt", n.chains = nc, 
             n.thin = nt, n.iter = ni, n.burnin = nb)

# 
end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 
# Calculate computation time


# Summarize posteriors
print(out1, dig = 3)

#determine if rhat values are <1.1 and save output in csv's to review after
rhat.1 <- data.frame(out1$Rhat)

rhat.final <- matrix(NA, nrow = 7, ncol = 10)

for(col in ncol(rhat.1)){ ## loop through columns to fill in matrix
  for(row in 1:nrow(rhat.1)){ ## loop through rows to fill in matrix
    
    rhat.final[,] <- if(rhat.1[row, col] >1.1){
      paste("True")
    } else {
      paste("False")
    }
    
  }
  
}

path.rhats <- "./rhats"

if(landcol == 4){
  file.name = paste("rhat_juv_ag", chem.names[rowIdx], ".csv", sep="")
  write.csv(rhat.final, file.path(path.rhats,file = file.name))
  
} else if (landcol ==6){
  file.name = paste("rhat_juv_dev", chem.names[rowIdx], ".csv", sep="")
  write.csv(rhat.final, file.path(path.rhats,file = file.name))
} else{
  file.name = paste("rhat_juv_BMP", chem.names[rowIdx], ".csv", sep="")
  write.csv(rhat.final, file.path(path.rhats,file = file.name))
}


# sum1 <- out$BUGSoutput$summary
# write.csv(sum1,'summary.csv',row.names=T)
# str(out)

############################
### Tables of Estimates ###
############################

#####################################################################
## Overall Estimated Occurrence in Juvenile Tissue (alphas.gammas) ##
#####################################################################

#everything separate
if(landcol == 4){ ## ag loop
  overall.occurrence.juv[compIdx,1] <- round(data.frame(mean(plogis(out1$sims.list$overall.mean))),2)
  overall.occurrence.juv[compIdx,2] <- round(data.frame(quantile(plogis(out1$sims.list$overall.mean),probs=c(0.975))),2)
  overall.occurrence.juv[compIdx,3] <- round(data.frame(quantile(plogis(out1$sims.list$overall.mean), probs=c(0.025))),2)
} else if (landcol == 6){ ## dev loop
  
  overall.occurrence.juv[compIdx,4] <- round(data.frame(mean(plogis(out1$sims.list$overall.mean))),2)
  overall.occurrence.juv[compIdx,5] <- round(data.frame(quantile(plogis(out1$sims.list$overall.mean),probs=c(0.975))),2)
  overall.occurrence.juv[compIdx,6] <- round(data.frame(quantile(plogis(out1$sims.list$overall.mean),probs=c(0.025))),2)
  
} else { ## BMP loop
  overall.occurrence.juv[compIdx,7] <- round(data.frame(mean(plogis(out1$sims.list$overall.mean))),2)
  overall.occurrence.juv[compIdx,8] <- round(data.frame(quantile(plogis(out1$sims.list$overall.mean),probs=c(0.975))),2)
  overall.occurrence.juv[compIdx,9] <- round(data.frame(quantile(plogis(out1$sims.list$overall.mean),probs=c(0.025))),2)
}
#everything pasted into one column (for each land-use type)
#ag
overall.occurrence.juv.2[compIdx,1] <- paste(round(overall.occurrence.juv[compIdx,1],2),"[",round(overall.occurrence.juv[compIdx,3],2),",", round(overall.occurrence.juv[compIdx,2],2),"]")
#dev
overall.occurrence.juv.2[compIdx,2] <- paste(round(overall.occurrence.juv[compIdx,4],2),"[",round(overall.occurrence.juv[compIdx,6],2),",", round(overall.occurrence.juv[compIdx,5],2),"]")
#BMP
overall.occurrence.juv.2[compIdx,3] <- paste(round(overall.occurrence.juv[compIdx,7],2),"[",round(overall.occurrence.juv[compIdx,9],2),",", round(overall.occurrence.juv[compIdx,8],2),"]")

############################################################
## Estimated Occurrence in Juvenile Tissue By Site (alphas)
############################################################

for(site in 1:7){
    #ag
  if(landcol == 4){
    site.occurrence.juv[compIdx,1,site] <- round(data.frame(apply(plogis(out1$sims.list$site.mean), 2, mean)),2)[site,]
    site.occurrence.juv[compIdx,2,site] <- round(data.frame(apply(plogis(out1$sims.list$site.mean), 2, quantile, probs=c(0.975))),2)[site,]
    site.occurrence.juv[compIdx,3,site] <- round(data.frame(apply(plogis(out1$sims.list$site.mean), 2, quantile, probs=c(0.025))),2)[site,]
  } else if(landcol == 6){
    #dev
    site.occurrence.juv[compIdx,4,site] <- round(data.frame(apply(plogis(out1$sims.list$site.mean), 2, mean)),2)[site,]
    site.occurrence.juv[compIdx,5,site] <- round(data.frame(apply(plogis(out1$sims.list$site.mean), 2, quantile, probs=c(0.975))),2)[site,]
    site.occurrence.juv[compIdx,6,site] <- round(data.frame(apply(plogis(out1$sims.list$site.mean), 2, quantile, probs=c(0.025))),2)[site,]
  }else{
    #BMP
    site.occurrence.juv[compIdx,7,site] <- round(data.frame(apply(plogis(out1$sims.list$site.mean), 2, mean)),2)[site,]
    site.occurrence.juv[compIdx,8,site] <- round(data.frame(apply(plogis(out1$sims.list$site.mean), 2, quantile, probs=c(0.975))),2)[site,]
    site.occurrence.juv[compIdx,9,site] <- round(data.frame(apply(plogis(out1$sims.list$site.mean), 2, quantile, probs=c(0.025))),2)[site,]
  }
}

## everything summarized into one column each (for each land-use type)

for(site in 1:7){
  #ag
  site.occurrence.juv.2[compIdx,1,site] <- paste(round(site.occurrence.juv[compIdx,1,site],2),"[",round(site.occurrence.juv[compIdx,3,site],2),",", round(site.occurrence.juv[compIdx,2,site],2),"]")
  #dev
  site.occurrence.juv.2[compIdx,2,site] <- paste(round(site.occurrence.juv[compIdx,4,site],2),"[",round(site.occurrence.juv[compIdx,6,site],2),",", round(site.occurrence.juv[compIdx,5,site],2),"]")
  #BMP
  site.occurrence.juv.2[compIdx,3,site] <- paste(round(site.occurrence.juv[compIdx,7,site],2),"[",round(site.occurrence.juv[compIdx,9,site],2),",", round(site.occurrence.juv[compIdx,8,site],2),"]")
}

################################################
## Land-use effect on occurrence (beta.gammas)
################################################

#everything separate
if(landcol== 4){
  #ag
  landuse.effect.occurrence.juv[compIdx,1] <- round(data.frame(mean(out1$sims.list$beta.gamma)),2)
  landuse.effect.occurrence.juv[compIdx,2] <- round(data.frame(quantile(out1$sims.list$beta.gamma, probs=c(0.975))),2)
  landuse.effect.occurrence.juv[compIdx,3] <- round(data.frame(quantile(out1$sims.list$beta.gamma,probs=c(0.025))),2)
} else if (landcol == 6){
  #dev
  landuse.effect.occurrence.juv[compIdx,4] <- round(data.frame(mean(out1$sims.list$beta.gamma)),2)
  landuse.effect.occurrence.juv[compIdx,5] <- round(data.frame(quantile(out1$sims.list$beta.gamma, probs=c(0.975))),2)
  landuse.effect.occurrence.juv[compIdx,6] <- round(data.frame(quantile(out1$sims.list$beta.gamma,probs=c(0.025))),2)
  
} else {
  #BMP
  landuse.effect.occurrence.juv[compIdx,7] <- round(data.frame(mean(out1$sims.list$beta.gamma)),2)
  landuse.effect.occurrence.juv[compIdx,8] <- round(data.frame(quantile(out1$sims.list$beta.gamma, probs=c(0.975))),2)
  landuse.effect.occurrence.juv[compIdx,9] <- round(data.frame(quantile(out1$sims.list$beta.gamma,probs=c(0.025))),2)
}
#everything pasted into one column (for each land-use type)
#ag
landuse.effect.occurrence.juv.2[compIdx,1] <- paste(round(landuse.effect.occurrence.juv[compIdx,1],2),"[",round(landuse.effect.occurrence.juv[compIdx,3],2),",", round(landuse.effect.occurrence.juv[compIdx,2],2),"]")
#dev
landuse.effect.occurrence.juv.2[compIdx,2] <- paste(round(landuse.effect.occurrence.juv[compIdx,4],2),"[",round(landuse.effect.occurrence.juv[compIdx,6],2),",", round(landuse.effect.occurrence.juv[compIdx,5],2),"]")
#BMP
landuse.effect.occurrence.juv.2[compIdx,3] <- paste(round(landuse.effect.occurrence.juv[compIdx,7],2),"[",round(landuse.effect.occurrence.juv[compIdx,9],2),",", round(landuse.effect.occurrence.juv[compIdx,8],2),"]")

################################################
#### pull out estimates for plotting later #### 
################################################

if(landcol== 4){ ## ag loop
  # Population-average parameters
  ests <- out1$summary[c("alpha.gamma","beta.gamma"),1]
  
  # Fake predictor
  X3 <- seq(min(unique(dat.ovary$Ag)), max(unique(dat.ovary$Ag)), length=100)
  
  
  # Linear predictor for population average effect
  lpredPopAve <- ests[1] + ests[2]*X3 
  plot(plogis(lpredPopAve)~X3, type ="l")
  
  # Number of simulations
  N <- out1$mcmc.info$n.samples
  
  linPredPopAve <- matrix(NA, ncol=length(X3), nrow=N ) #container for predicted values
  
  mu.A <- out1$sims.list$alpha.gamma
  mu.B <- out1$sims.list$beta.gamma
  
  for(i in 1:N){
    for(t in 1:length(X3)){
      linPredPopAve[i,t] <- mu.A[i] + mu.B[i]*X3[t] 
    }
  }
  
  ## store output for panel plot later ##
  probPredPopAve.Ag[,,compIdx] <- plogis(linPredPopAve)
  meanProbPopAve.Ag[compIdx,] <- apply(probPredPopAve.Ag[,,compIdx], 2, mean)
  upperCI.PopAve.Ag[compIdx,] <- apply(probPredPopAve.Ag[,,compIdx], 2, quantile, probs=c(0.975) )
  lowerCIA.PopAve.Ag[compIdx,] <- apply(probPredPopAve.Ag[,,compIdx], 2, quantile, probs=c(0.025) )
  
}else if(landcol == 6){ ## dev loop
  
  # Population-average parameters
  ests <- out1$summary[c("alpha.gamma","beta.gamma"),1]
  
  # Fake predictor
  X3 <- seq(min(unique(dat.ovary$Dev)), max(unique(dat.ovary$Dev)), length=100)
  
  # Linear predictor for population average effect
  lpredPopAve <- ests[1] + ests[2]*X3 
  plot(plogis(lpredPopAve)~X3, type ="l")
  
  # Number of simulations
  N <- out1$mcmc.info$n.samples
  
  linPredPopAve <- matrix(NA, ncol=length(X3), nrow=N ) #container for predicted values
  
  mu.A <- out1$sims.list$alpha.gamma
  mu.B <- out1$sims.list$beta.gamma
  
  for(i in 1:N){
    for(t in 1:length(X3)){
      linPredPopAve[i,t] <- mu.A[i] + mu.B[i]*X3[t] 
    }
  }
  
  probPredPopAve.Dev[,,compIdx] <- plogis(linPredPopAve)
  meanProbPopAve.Dev[compIdx,] <- apply(probPredPopAve.Dev[,,compIdx], 2, mean)
  upperCI.PopAve.Dev[compIdx,] <- apply(probPredPopAve.Dev[,,compIdx], 2, quantile, probs=c(0.975) )
  lowerCIA.PopAve.Dev[compIdx,] <- apply(probPredPopAve.Dev[,,compIdx], 2, quantile, probs=c(0.025) )
  
} else{ ## BMP loop
  # Population-average parameters
  ests <- out1$summary[c("alpha.gamma","beta.gamma"),1]
  
  # Fake predictor
  X3 <- seq(min(unique(dat.ovary$BMP)), max(unique(dat.ovary$BMP)), length=100)
  
  # Linear predictor for population average effect
  lpredPopAve <- ests[1] + ests[2]*X3 
  plot(plogis(lpredPopAve)~X3, type ="l")
  
  # Number of simulations
  N <- out1$mcmc.info$n.samples
  
  linPredPopAve <- matrix(NA, ncol=length(X3), nrow=N ) #container for predicted values
  
  mu.A <- out1$sims.list$alpha.gamma
  mu.B <- out1$sims.list$beta.gamma
  
  for(i in 1:N){
    for(t in 1:length(X3)){
      linPredPopAve[i,t] <- mu.A[i] + mu.B[i]*X3[t] 
    }
  }
  
  probPredPopAve.BMP[,,compIdx] <- plogis(linPredPopAve)
  meanProbPopAve.BMP[compIdx,] <- apply(probPredPopAve.BMP[,,compIdx], 2, mean)
  upperCI.PopAve.BMP[compIdx,] <- apply(probPredPopAve.BMP[,,compIdx], 2, quantile, probs=c(0.975) )
  lowerCIA.PopAve.BMP[compIdx,] <- apply(probPredPopAve.BMP[,,compIdx], 2, quantile, probs=c(0.025) )
}

}  ## FOR-LOOP END: Landuse
  
  compIdx <- compIdx+1
  rowIdx <- rowIdx+1
  
} ## FOR-LOOP END: Contaminants



##### ##### ##### ##### ##### ##### ##### 
##### END OF MODEL AND TABLE FOR LOOPS ##### 
##### ##### ##### ##### ##### ##### ##### 

## Save all tables ##

write.csv(overall.occurrence.juv,"overall.occurrence.juv.csv")
write.csv(overall.occurrence.juv.2, "overall.occurrence.juv.2.csv")
write.csv(site.occurrence.juv, "site.occurrence.juv.csv")
write.csv(site.occurrence.juv.2, "site.occurrence.juv.2.csv")
write.csv(landuse.effect.occurrence.juv, "landuse.effect.occurrence.juv.csv")
write.csv(landuse.effect.occurrence.juv.2, "landuse.effect.occurrence.juv.2.csv")

## Panel Plots ##
################# PLOT ############################

############
#### Ag ####
############

tiff("Occurence_Ag_Juv.tiff", 
     height=6, width=6, units='in', res=450, compression='lzw')
def.par <- par(no.readonly = TRUE) 		

size.labels = 1
size.text = 1
axissize <- 0.8
x.label = 'Standardized[logit(Agricultural Land)]'
y.label = 'Predicted Probability of Occurrence'


nf <- layout(matrix(c(1:12),nrow=4,ncol=3,byrow=TRUE),  TRUE) 
layout.show(nf)
par(mar=c(0.0,0.1,0.1,0.1),oma=c(3,3,0,1),mai=c(0.0,0.05,0.05,0) )

# Compound-specific plots for Each Land Use

Ymin <- 0
Ymax <- 1

x <- seq(min(unique(dat.ovary$Ag)), max(unique(dat.ovary$Ag)), length=100)

y <- runif(100, min = 0, max = 1)

ylabnums <- c(0.0,0.2,0.4,0.8,1.0)

for(i in 1:12){
  
  plot(x,y, ylim=c(Ymin, Ymax), xlim=c(min(x), max(x)), axes=F, ylab='', xlab='', type='n')
  
  if( i <= 9){
    axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, labels=F ) 
  } else {
    axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01)
  }	
  
  if( i ==1 | i ==4 | i ==7 | i ==10){
    axis(side=2,cex.axis=axissize , mgp=c(0,0.3,0),tck= -0.01, las=1,at=format(pretty(ylabnums), digits=2), 
         labels=format(pretty(ylabnums), digits=2))
  } else {
    axis(side=2,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, labels=F)
  }	
  
  #axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01)
  axis(side=2,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, las=1, labels=F)
  # 
  # 
  i.for <- order(x)
  i.back <- order(x, decreasing = TRUE )
  x.polygon <- c(x[i.for] , x[i.back] )
  Lower <- lowerCIA.PopAve.Ag[i,]
  Upper <- upperCI.PopAve.Ag[i,]
  y.polygon <- c( Lower[i.for] , Upper[i.back] )
  polygon( x.polygon , y.polygon , col = "lightgray" , border = NA )
  # 
  lines(x,meanProbPopAve.Ag[i,], lwd = 1, col="black", lty = 1)
  # 
  text(-2, 0.9, i, cex=1)
  #text(0.6,0.9,'(7021)',cex=0.8)
  # 
  mtext(y.label, line = 1.5, side = 2, cex = size.text,outer=T)
  mtext(x.label, line = 1.5, side = 1, cex = size.text, outer=T)
  # 
  box()
}
par(def.par)
dev.off()

##############
#### Dev #####
##############

tiff("Occurence_Dev_Juv.tiff", 
     height=6, width=6, units='in', res=450, compression='lzw')
def.par <- par(no.readonly = TRUE) 		

size.labels = 1
size.text = 1
axissize <- 0.8
x.label = 'Standardized[logit(Developed Land)]'
y.label = 'Predicted Probability of Occurrence'


nf <- layout(matrix(c(1:12),nrow=4,ncol=3,byrow=TRUE),  TRUE) 
layout.show(nf)
par(mar=c(0.0,0.1,0.1,0.1),oma=c(3,3,0,1),mai=c(0.0,0.05,0.05,0) )

# Compound-specific plots for Each Land Use

Ymin <- 0
Ymax <- 1

x <- seq(min(unique(dat.ovary$Dev)), max(unique(dat.ovary$Dev)), length=100)

y <- runif(100, min = 0, max = 1)

ylabnums <- c(0.0,0.2,0.4,0.8,1.0)

for(i in 1:12){
  
  plot(x,y, ylim=c(Ymin, Ymax), xlim=c(min(x), max(x)), axes=F, ylab='', xlab='', type='n')
  
  if( i <= 9){
    axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, labels=F ) 
  } else {
    axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01)
  }	
  
  if( i ==1 | i ==4 | i ==7 | i ==10){
    axis(side=2,cex.axis=axissize , mgp=c(0,0.3,0),tck= -0.01, las=1,at=format(pretty(ylabnums), digits=2), 
         labels=format(pretty(ylabnums), digits=2))
  } else {
    axis(side=2,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, labels=F)
  }	
  
  #axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01)
  axis(side=2,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, las=1, labels=F)
  # 
  # 
  i.for <- order(x)
  i.back <- order(x, decreasing = TRUE )
  x.polygon <- c(x[i.for] , x[i.back] )
  Lower <- lowerCIA.PopAve.Dev[i,]
  Upper <- upperCI.PopAve.Dev[i,]
  y.polygon <- c( Lower[i.for] , Upper[i.back] )
  polygon( x.polygon , y.polygon , col = "lightgray" , border = NA )
  # 
  lines(x,meanProbPopAve.Dev[i,], lwd = 1, col="black", lty = 1)
  # 
  text(-1, 0.9, i, cex=1)
  #text(0.6,0.9,'(7021)',cex=0.8)
  # 
  mtext(y.label, line = 1.5, side = 2, cex = size.text,outer=T)
  mtext(x.label, line = 1.5, side = 1, cex = size.text, outer=T)
  # 
  box()
}
par(def.par)
dev.off()

##############
#### BMP #####
##############

tiff("Occurence_BMP_Juv.tiff", 
     height=6, width=6, units='in', res=450, compression='lzw')
def.par <- par(no.readonly = TRUE) 		

size.labels = 1
size.text = 1
axissize <- 0.8
x.label = 'Standardized[BMP]'
y.label = 'Predicted Probability of Occurrence'


nf <- layout(matrix(c(1:12),nrow=4,ncol=3,byrow=TRUE),  TRUE) 
layout.show(nf)
par(mar=c(0.0,0.1,0.1,0.1),oma=c(3,3,0,1),mai=c(0.0,0.05,0.05,0) )

# Compound-specific plots for Each Land Use

Ymin <- 0
Ymax <- 1

x <- seq(min(unique(dat.ovary$BMP)), max(unique(dat.ovary$BMP)), length=100)
y <- runif(100, min = 0, max = 1)

ylabnums <- c(0.0,0.2,0.4,0.8,1.0)

for(i in 1:12){
  
  plot(x,y, ylim=c(Ymin, Ymax), xlim=c(min(x), max(x)), axes=F, ylab='', xlab='', type='n')
  
  if( i <= 9){
    axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, labels=F ) 
  } else {
    axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01)
  }	
  
  if( i ==1 | i ==4 | i ==7 | i ==10){
    axis(side=2,cex.axis=axissize , mgp=c(0,0.3,0),tck= -0.01, las=1,at=format(pretty(ylabnums), digits=2), 
         labels=format(pretty(ylabnums), digits=2))
  } else {
    axis(side=2,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, labels=F)
  }	
  
  #axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01)
  axis(side=2,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, las=1, labels=F)
  # 
  # 
  i.for <- order(x)
  i.back <- order(x, decreasing = TRUE )
  x.polygon <- c(x[i.for] , x[i.back] )
  Lower <- lowerCIA.PopAve.BMP[i,]
  Upper <- upperCI.PopAve.BMP[i,]
  y.polygon <- c( Lower[i.for] , Upper[i.back] )
  polygon( x.polygon , y.polygon , col = "lightgray" , border = NA )
  # 
  lines(x,meanProbPopAve.BMP[i,], lwd = 1, col="black", lty = 1)
  # 
  text(4, 0.9, i, cex=1)
  #text(0.6,0.9,'(7021)',cex=0.8)
  # 
  mtext(y.label, line = 1.5, side = 2, cex = size.text,outer=T)
  mtext(x.label, line = 1.5, side = 1, cex = size.text, outer=T)
  # 
  box()
}
par(def.par)
dev.off()

####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### 

####### ####### ####### ####### 
####### Ovary Tissue Model ####
####### ####### ####### #######

####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### ####### 

# 15 compounds

#############################################
### Matrices to fill with model estimates ###
#############################################

## Overall ##
overall.occurrence.ovary <- matrix(ncol=9, nrow = 15)
overall.occurrence.ovary <- as.data.frame(overall.occurrence.ovary)
rownames(overall.occurrence.ovary) <- c("Hexachlorobenzene",
                                        "Heptachlor",
                               "Octachlorostyrene",
                               "ppdde",
                               "o,p'-DDT",
                               "Mirex",
                               "Galaxolide",
                               "Diethyl phthalate",
                               "DEHP",
                               "2-Methyl naphthalene",
                               "Anthracene",
                               "Naphthalene",
                               "Phenanthrene",
                               "Benzophenone",
                               "Methyl Triclosan")
colnames(overall.occurrence.ovary) <- c("Ag",
                               "Ag.UP",
                               "Ag.LOW",
                               "Dev",
                               "Dev.UP",
                               "Dev.LOW",
                               "BMP",
                               "BMP.UP",
                               "BMP.LOW")
## for table in manuscript ##
overall.occurrence.ovary.2 <- matrix(ncol=3, nrow = 15)
overall.occurrence.ovary.2 <- as.data.frame(overall.occurrence.ovary.2)
rownames(overall.occurrence.ovary.2) <- c("Hexachlorobenzene",
                                          "Heptachlor",
                                          "Octachlorostyrene",
                                          "ppdde",
                                          "o,p'-DDT",
                                          "Mirex",
                                          "Galaxolide",
                                          "Diethyl phthalate",
                                          "DEHP",
                                          "2-Methyl naphthalene",
                                          "Anthracene",
                                          "Naphthalene",
                                          "Phenanthrene",
                                          "Benzophenone",
                                          "Methyl Triclosan")
colnames(overall.occurrence.ovary.2) <- c("Ag",
                                          "Dev",
                                          "BMP")
################
## By sites ##
################

## overall ##
site.occurrence.ovary <- array(NA, c(15,9,10)) #contaminants, land-use (means, CIs), sites
colnames(site.occurrence.ovary) <- c("Ag",
                                     "Ag.UP",
                                     "Ag.LOW",
                                     "Dev",
                                     "Dev.UP",
                                     "Dev.LOW",
                                     "BMP",
                                     "BMP.UP",
                                     "BMP.LOW")
rownames(site.occurrence.ovary) <- c("Hexachlorobenzene",
                                     "Heptachlor",
                                     "Octachlorostyrene",
                                     "ppdde",
                                     "o,p'-DDT",
                                     "Mirex",
                                     "Galaxolide",
                                     "Diethyl phthalate",
                                     "DEHP",
                                     "2-Methyl naphthalene",
                                     "Anthracene",
                                     "Naphthalene",
                                     "Phenanthrene",
                                     "Benzophenone",
                                     "Methyl Triclosan")

## for table in manuscript ##
site.occurrence.ovary.2 <- array(NA, c(15,3,10)) #contaminants, land-use (means, CIs), sites
colnames(site.occurrence.ovary.2) <- c("Ag",
                                     "Dev",
                                     "BMP")
rownames(site.occurrence.ovary.2) <- c("Hexachlorobenzene",
                                       "Heptachlor",
                                       "Octachlorostyrene",
                                       "ppdde",
                                       "o,p'-DDT",
                                       "Mirex",
                                       "Galaxolide",
                                       "Diethyl phthalate",
                                       "DEHP",
                                       "2-Methyl naphthalene",
                                       "Anthracene",
                                       "Naphthalene",
                                       "Phenanthrene",
                                       "Benzophenone",
                                       "Methyl Triclosan")
####################
## land-use effects
####################

## overall ##

landuse.effect.occurrence.ovary <- matrix(ncol=9, nrow = 15)
landuse.effect.occurrence.ovary <- as.data.frame(landuse.effect.occurrence.ovary)
rownames(landuse.effect.occurrence.ovary) <- c("Hexachlorobenzene",
                                               "Heptachlor",
                                               "Octachlorostyrene",
                                               "ppdde",
                                               "o,p'-DDT",
                                               "Mirex",
                                               "Galaxolide",
                                               "Diethyl phthalate",
                                               "DEHP",
                                               "2-Methyl naphthalene",
                                               "Anthracene",
                                               "Naphthalene",
                                               "Phenanthrene",
                                               "Benzophenone",
                                               "Methyl Triclosan")
colnames(landuse.effect.occurrence.ovary) <- c("Ag",
                                        "Ag.UP",
                                        "Ag.LOW",
                                        "Dev",
                                        "Dev.UP",
                                        "Dev.LOW",
                                        "BMP",
                                        "BMP.UP",
                                        "BMP.LOW")

## for table in manuscript ##

landuse.effect.occurrence.ovary.2 <- matrix(ncol=3, nrow = 15)
landuse.effect.occurrence.ovary.2 <- as.data.frame(landuse.effect.occurrence.ovary.2)
rownames(landuse.effect.occurrence.ovary.2) <- c("Hexachlorobenzene",
                                                 "Heptachlor",
                                                 "Octachlorostyrene",
                                                 "ppdde",
                                                 "o,p'-DDT",
                                                 "Mirex",
                                                 "Galaxolide",
                                                 "Diethyl phthalate",
                                                 "DEHP",
                                                 "2-Methyl naphthalene",
                                                 "Anthracene",
                                                 "Naphthalene",
                                                 "Phenanthrene",
                                                 "Benzophenone",
                                                 "Methyl Triclosan")
colnames(landuse.effect.occurrence.ovary.2) <- c("Ag",
                                          "Dev",
                                          "BMP")

####################################################
## Matrices to fill in for plotting panel plots ##
####################################################

#lines
probPredPopAve.Ag <- array(NA, c(210000,100,15) )
probPredPopAve.Dev <- array(NA, c(210000,100,15) )
probPredPopAve.BMP <- array(NA, c(210000,100,15) )

#means
meanProbPopAve.Ag <- matrix(NA, nrow = 15, ncol = 100)
meanProbPopAve.Dev <- matrix(NA, nrow = 15, ncol = 100)
meanProbPopAve.BMP <- matrix(NA, nrow = 15, ncol = 100)

#upper
upperCI.PopAve.Ag <- matrix(NA, nrow = 15, ncol = 100)
upperCI.PopAve.Dev <- matrix(NA, nrow = 15, ncol = 100)
upperCI.PopAve.BMP <- matrix(NA, nrow = 15, ncol = 100)

#lower
lowerCIA.PopAve.Ag <- matrix(NA, nrow = 15, ncol = 100)
lowerCIA.PopAve.Dev <- matrix(NA, nrow = 15, ncol = 100)
lowerCIA.PopAve.BMP <- matrix(NA, nrow = 15, ncol = 100)

dat.col.names <- colnames(dat.ovary)
chem.names <- dat.col.names[10:24]
dat.ovary <- data.frame(dat.ovary)

compIdx <- 1
rowIdx <- 1

for(colIdx in 10:24){  ## FOR-LOOP START: Loop through contaminants
  
  for(landcol in seq(from = 4, to=8, by=2)){ ## FOR-LOOP START: Loop through landuse types
    
    print("NEW LOOP ITERATION.  WORKING ON:")
    print(dat.col.names[colIdx])
    print(dat.col.names[landcol])
    
    
  # Define the model in the BUGS language and write a text file
  sink("model.ovary.txt")
  cat("
model {
 
 
# Likelihood:
# Level-1 of the model
for (i in 1:n){
   y[i] ~ dbin(p[i],1)  # distributional assumption
   p[i] <- exp(lp[i])/(1+exp(lp[i])) # logit link function
   lp[i] <- alpha[group[i]]  # linear predictor   
  }

# Level-2 of the model
    for(j in 1:J){
    alpha[j] ~ dnorm(mu.alpha[j],tau.alpha)
    mu.alpha[j] <- alpha.gamma + beta.gamma *landuse[j]
    }
    
    
    # Priors and derived quantities
    sigma ~ dunif(0, 100)
    alpha.gamma ~ dnorm(0, 0.0001)
    beta.gamma ~ dnorm(0,0.0001)
    sigma.alpha ~ dunif(0,100)
    
   overall.mean <- alpha.gamma
   site.mean <- alpha
   landuse.effect <- beta.gamma
   
   # Derived quantities
    tau <- pow(sigma,-2) # precision

    tau.alpha <- pow(sigma.alpha,-2) # precision
    
    } # end model
    ",fill = TRUE)
  sink()
  
  # Number of sites
  J <- length(unique(dat.ovary$Site))
  
  # Site indicator
  G <- as.numeric(as.factor(as.numeric(as.factor(dat.ovary$Site))))
  
  #predictor scaled and log transformed for model fit
  landuse <- unique(dat.ovary[,landcol])
  
  #log transform concentrations
  Y <- dat.ovary[,colIdx]
  n <- length(Y)
  
  # Load data
  data <- list(y = Y,
               n = n,
               landuse=landuse,
               group = G,
               J = J )
  
  # Initial values
  inits <- function (){
    list (alpha.gamma = rnorm(1),
          sigma=runif(1),
          beta.gamma=rnorm(1),
          sigma.alpha=runif(1))
  }
  
  # Parameters monitored
  parameters <- c("alpha.gamma",
                  "sigma",
                  "beta.gamma",
                  "sigma.alpha",
                  "mu.alpha",
                  "alpha",
                  "overall.mean",
                  "site.mean",
                  "landuse.effect")
  
  # MCMC settings
  ni <- 100000
  nt <- 1
  nb <- 30000
  nc <- 3
  
  start.time = Sys.time()         # Set timer 
  # Call JAGS from R 
  
  out1 <- jags(data, inits, parameters, "model.ovary.txt", n.chains = nc, 
               n.thin = nt, n.iter = ni, n.burnin = nb)
  
  # 
  end.time = Sys.time()
  elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
  cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 
  # Calculate computation time
  
  
  # Summarize posteriors
  print(out1, dig = 3)
  
  #determine if rhat values are <1.1 and save output in csv's to review after
  rhat.1 <- data.frame(out1$Rhat)

  rhat.final <- matrix(NA, nrow = 10, ncol = 10)
  
  for(col in ncol(rhat.1)){ ## loop through columns to fill in matrix
    for(row in 1:nrow(rhat.1)){ ## loop through rows to fill in matrix
      
      rhat.final[,] <- if(rhat.1[row, col] >1.1){
        paste("True")
      } else {
        paste("False")
      }
      
    }
    
  }
  
  
  path.rhats <- "./rhats"
  
  if(landcol == 4){ ## ag loop
    file.name = paste("rhat_ovary_ag", chem.names[rowIdx], ".csv", sep="")
    write.csv(rhat.final, file.path(path.rhats,file = file.name))
    
  } else if (landcol ==6){ ## dev loop
    file.name = paste("rhat_ovary_dev", chem.names[rowIdx], ".csv", sep="")
    write.csv(rhat.final, file.path(path.rhats,file = file.name))
  } else{ ## BMP loops
    file.name = paste("rhat_ovary_BMP", chem.names[rowIdx], ".csv", sep="")
    write.csv(rhat.final, file.path(path.rhats,file = file.name))
  }
  
  
  ## check convergence via Rhat values ##
  # See what max Rhat value is (<1.1)
  
  max(out1$summary[, c("Rhat")])
  
  # sum1 <- out$BUGSoutput$summary
  # write.csv(sum1,'summary.csv',row.names=T)
  # str(out)
  
  
  ############################
  ### Tables of Estimates ###
  ############################
  
  ########################################################
  ## Overall Estimated Occurrence in Ovary Tissue (mu.alphas)
  ####################################################
  dim(plogis(out1$sims.list$overall.mean))
  #everything separate
  if(landcol == 4){ ## ag loop
  overall.occurrence.ovary[compIdx,1] <- round(data.frame(mean(plogis(out1$sims.list$overall.mean))),2)
  overall.occurrence.ovary[compIdx,2] <- round(data.frame(quantile(plogis(out1$sims.list$overall.mean),probs=c(0.975))),2)
  overall.occurrence.ovary[compIdx,3] <- round(data.frame(quantile(plogis(out1$sims.list$overall.mean),probs=c(0.025))),2)
  } else if (landcol == 6){ ## dev loop
    
  overall.occurrence.ovary[compIdx,4] <- round(data.frame(mean(plogis(out1$sims.list$overall.mean))),2)
  overall.occurrence.ovary[compIdx,5] <- round(data.frame(quantile(plogis(out1$sims.list$overall.mean),probs=c(0.975))),2)
  overall.occurrence.ovary[compIdx,6] <- round(data.frame(quantile(plogis(out1$sims.list$overall.mean),probs=c(0.025))),2)
  
  } else { ## BMP loop
    overall.occurrence.ovary[compIdx,7] <- round(data.frame(mean(plogis(out1$sims.list$overall.mean))),2)
    overall.occurrence.ovary[compIdx,8] <- round(data.frame(quantile(plogis(out1$sims.list$overall.mean),probs=c(0.975))),2)
    overall.occurrence.ovary[compIdx,9] <- round(data.frame(quantile(plogis(out1$sims.list$overall.mean),probs=c(0.025))),2)
  }
  #everything pasted into one column (for each land-use type)
  #ag
  overall.occurrence.ovary.2[compIdx,1] <- paste(round(overall.occurrence.ovary[compIdx,1],2),"[",round(overall.occurrence.ovary[compIdx,3],2),",", round(overall.occurrence.ovary[compIdx,2],2),"]")
  #dev
  overall.occurrence.ovary.2[compIdx,2] <- paste(round(overall.occurrence.ovary[compIdx,4],2),"[",round(overall.occurrence.ovary[compIdx,6],2),",", round(overall.occurrence.ovary[compIdx,5],2),"]")
  #BMP
  overall.occurrence.ovary.2[compIdx,3] <- paste(round(overall.occurrence.ovary[compIdx,7],2),"[",round(overall.occurrence.ovary[compIdx,9],2),",", round(overall.occurrence.ovary[compIdx,8],2),"]")
  
  ############################################################
  ## Estimated Occurrence in Ovary Tissue By Site (alphas)
  ############################################################
  
  for(site in 1:10){
  #ag
  if(landcol == 4){
  site.occurrence.ovary[compIdx,1,site] <- round(data.frame(apply(plogis(out1$sims.list$site.mean), 2, mean)),2)[site,]
  site.occurrence.ovary[compIdx,2,site] <- round(data.frame(apply(plogis(out1$sims.list$site.mean), 2, quantile, probs=c(0.975))),2)[site,]
  site.occurrence.ovary[compIdx,3,site] <- round(data.frame(apply(plogis(out1$sims.list$site.mean), 2, quantile, probs=c(0.025))),2)[site,]
  } else if(landcol == 6){
  #dev
  site.occurrence.ovary[compIdx,4,site] <- round(data.frame(apply(plogis(out1$sims.list$site.mean), 2, mean)),2)[site,]
  site.occurrence.ovary[compIdx,5,site] <- round(data.frame(apply(plogis(out1$sims.list$site.mean), 2, quantile, probs=c(0.975))),2)[site,]
  site.occurrence.ovary[compIdx,6,site] <- round(data.frame(apply(plogis(out1$sims.list$site.mean), 2, quantile, probs=c(0.025))),2)[site,]
  }else{
  #BMP
  site.occurrence.ovary[compIdx,7,site] <- round(data.frame(apply(plogis(out1$sims.list$site.mean), 2, mean)),2)[site,]
  site.occurrence.ovary[compIdx,8,site] <- round(data.frame(apply(plogis(out1$sims.list$site.mean), 2, quantile, probs=c(0.975))),2)[site,]
  site.occurrence.ovary[compIdx,9,site] <- round(data.frame(apply(plogis(out1$sims.list$site.mean), 2, quantile, probs=c(0.025))),2)[site,]
  }
  }

  ## everything summarized into one column each (for each land-use type)
  
  for(site in 1:10){
  #ag
  site.occurrence.ovary.2[compIdx,1,site] <- paste(round(site.occurrence.ovary[compIdx,1,site],2),"[",round(site.occurrence.ovary[compIdx,3,site],2),",", round(site.occurrence.ovary[compIdx,2,site],2),"]")
  #dev
  site.occurrence.ovary.2[compIdx,2,site] <- paste(round(site.occurrence.ovary[compIdx,4,site],2),"[",round(site.occurrence.ovary[compIdx,6,site],2),",", round(site.occurrence.ovary[compIdx,5,site],2),"]")
  #BMP
  site.occurrence.ovary.2[compIdx,3,site] <- paste(round(site.occurrence.ovary[compIdx,7,site],2),"[",round(site.occurrence.ovary[compIdx,9,site],2),",", round(site.occurrence.ovary[compIdx,8,site],2),"]")
  }
  ################################################
  ## Land-use effect on occurrence (beta.gammas)
  ################################################
  
  #everything separate
  if(landcol== 4){
    #ag
    landuse.effect.occurrence.ovary[compIdx,1] <- round(data.frame(mean(out1$sims.list$beta.gamma)),2)
    landuse.effect.occurrence.ovary[compIdx,2] <- round(data.frame(quantile(out1$sims.list$beta.gamma, probs=c(0.975))),2)
    landuse.effect.occurrence.ovary[compIdx,3] <- round(data.frame(quantile(out1$sims.list$beta.gamma,probs=c(0.025))),2)
  } else if (landcol == 6){
    #dev
    landuse.effect.occurrence.ovary[compIdx,4] <- round(data.frame(mean(out1$sims.list$beta.gamma)),2)
    landuse.effect.occurrence.ovary[compIdx,5] <- round(data.frame(quantile(out1$sims.list$beta.gamma, probs=c(0.975))),2)
    landuse.effect.occurrence.ovary[compIdx,6] <- round(data.frame(quantile(out1$sims.list$beta.gamma,probs=c(0.025))),2)
    
  } else {
    #BMP
    landuse.effect.occurrence.ovary[compIdx,7] <- round(data.frame(mean(plogis(out1$sims.list$beta.gamma))),2)
    landuse.effect.occurrence.ovary[compIdx,8] <- round(data.frame(quantile(out1$sims.list$beta.gamma, probs=c(0.975))),2)
    landuse.effect.occurrence.ovary[compIdx,9] <- round(data.frame(quantile(out1$sims.list$beta.gamma,probs=c(0.025))),2)
  }
  #everything pasted into one column (for each land-use type)
  #ag
  landuse.effect.occurrence.ovary.2[compIdx,1] <- paste(round(landuse.effect.occurrence.ovary[compIdx,1],2),"[",round(landuse.effect.occurrence.ovary[compIdx,3],2),",", round(landuse.effect.occurrence.ovary[compIdx,2],2),"]")
  #dev
  landuse.effect.occurrence.ovary.2[compIdx,2] <- paste(round(landuse.effect.occurrence.ovary[compIdx,4],2),"[",round(landuse.effect.occurrence.ovary[compIdx,6],2),",", round(landuse.effect.occurrence.ovary[compIdx,5],2),"]")
  #BMP
  landuse.effect.occurrence.ovary.2[compIdx,3] <- paste(round(landuse.effect.occurrence.ovary[compIdx,7],2),"[",round(landuse.effect.occurrence.ovary[compIdx,9],2),",", round(landuse.effect.occurrence.ovary[compIdx,8],2),"]")
  
  
  ################################################
  #### Pull out estimates for plotting later ####  
  ################################################
  
  if(landcol== 4){ ##ag loop
  # Population-average parameters
  ests <- out1$summary[c("alpha.gamma","beta.gamma"),1]
  
  # Fake predictor
  X3 <- seq(min(unique(dat.ovary$Ag)), max(unique(dat.ovary$Ag)), length=100)
  
  
  # Linear predictor for population average effect
  lpredPopAve <- ests[1] + ests[2]*X3 
  plot(plogis(lpredPopAve)~X3, type ="l")
  
  # Number of simulations
  N <- out1$mcmc.info$n.samples
  
  linPredPopAve <- matrix(NA, ncol=length(X3), nrow=N ) #container for predicted values
  
  mu.A <- out1$sims.list$alpha.gamma
  mu.B <- out1$sims.list$beta.gamma
  
  for(i in 1:N){
    for(t in 1:length(X3)){
      linPredPopAve[i,t] <- mu.A[i] + mu.B[i]*X3[t] 
    }
  }
  
  ## store output for panel plot later ##
  probPredPopAve.Ag[,,compIdx] <- plogis(linPredPopAve)
  meanProbPopAve.Ag[compIdx,] <- apply(probPredPopAve.Ag[,,compIdx], 2, mean)
  upperCI.PopAve.Ag[compIdx,] <- apply(probPredPopAve.Ag[,,compIdx], 2, quantile, probs=c(0.975) )
  lowerCIA.PopAve.Ag[compIdx,] <- apply(probPredPopAve.Ag[,,compIdx], 2, quantile, probs=c(0.025) )
 
  } else if(landcol == 6){ ## dev loop
    
    # Population-average parameters
    ests <- out1$summary[c("alpha.gamma","beta.gamma"),1]
    
    # Fake predictor
    X3 <- seq(min(unique(dat.ovary$Dev)), max(unique(dat.ovary$Dev)), length=100)
    
    # Linear predictor for population average effect
    lpredPopAve <- ests[1] + ests[2]*X3 
    plot(plogis(lpredPopAve)~X3, type ="l")
    
    # Number of simulations
    N <- out1$mcmc.info$n.samples
    
    linPredPopAve <- matrix(NA, ncol=length(X3), nrow=N ) #container for predicted values
    
    mu.A <- out1$sims.list$alpha.gamma
    mu.B <- out1$sims.list$beta.gamma
    
    for(i in 1:N){
      for(t in 1:length(X3)){
        linPredPopAve[i,t] <- mu.A[i] + mu.B[i]*X3[t] 
      }
    }
    
    probPredPopAve.Dev[,,compIdx] <- plogis(linPredPopAve)
    meanProbPopAve.Dev[compIdx,] <- apply(probPredPopAve.Dev[,,compIdx], 2, mean)
    upperCI.PopAve.Dev[compIdx,] <- apply(probPredPopAve.Dev[,,compIdx], 2, quantile, probs=c(0.975) )
    lowerCIA.PopAve.Dev[compIdx,] <- apply(probPredPopAve.Dev[,,compIdx], 2, quantile, probs=c(0.025) )
    
  } else{ ## BMP loop
    # Population-average parameters
    ests <- out1$summary[c("alpha.gamma","beta.gamma"),1]
    
    # Fake predictor
    X3 <- seq(min(unique(dat.ovary$BMP)), max(unique(dat.ovary$BMP)), length=100)
    
    # Linear predictor for population average effect
    lpredPopAve <- ests[1] + ests[2]*X3 
    plot(plogis(lpredPopAve)~X3, type ="l")
    
    # Number of simulations
    N <- out1$mcmc.info$n.samples
    
    linPredPopAve <- matrix(NA, ncol=length(X3), nrow=N ) #container for predicted values
    
    mu.A <- out1$sims.list$alpha.gamma
    mu.B <- out1$sims.list$beta.gamma
    
    for(i in 1:N){
      for(t in 1:length(X3)){
        linPredPopAve[i,t] <- mu.A[i] + mu.B[i]*X3[t] 
      }
    }
    
    probPredPopAve.BMP[,,compIdx] <- plogis(linPredPopAve)
    meanProbPopAve.BMP[compIdx,] <- apply(probPredPopAve.BMP[,,compIdx], 2, mean)
    upperCI.PopAve.BMP[compIdx,] <- apply(probPredPopAve.BMP[,,compIdx], 2, quantile, probs=c(0.975) )
    lowerCIA.PopAve.BMP[compIdx,] <- apply(probPredPopAve.BMP[,,compIdx], 2, quantile, probs=c(0.025) )
  }
  
  }  ## END OF FOR-LOOP: Landuse
  
  compIdx <- compIdx+1 
  rowIdx <- rowIdx+1
  
}   ## END OF FOR-LOOP: Contaminants


##### ##### ##### ##### ##### ##### ##### 
##### END OF MODEL AND TABLE FOR LOOPS ##### 
##### ##### ##### ##### ##### ##### ##### 

## Save all tables ##

write.csv(overall.occurrence.ovary,"overall.occurrence.ovary.csv")
write.csv(overall.occurrence.ovary.2, "overall.occurrence.ovary.2.csv")
write.csv(site.occurrence.ovary, "site.occurrence.ovary.csv")
write.csv(site.occurrence.ovary.2, "site.occurrence.ovary.2.csv")
write.csv(landuse.effect.occurrence.ovary, "landuse.effect.occurrence.ovary.csv")
write.csv(landuse.effect.occurrence.ovary.2, "landuse.effect.occurrence.ovary.2.csv")


## Panel Plots ##
################# PLOT ############################

############
#### Ag ####
############

  tiff("Occurence_Ag_Ovary.tiff", 
       height=6, width=6, units='in', res=450, compression='lzw')
  def.par <- par(no.readonly = TRUE) 		
  
  size.labels = 1
  size.text = 1
  axissize <- 0.8
  x.label = 'Standardized[logit(Agricultural Land)]'
  y.label = 'Predicted Probability of Occurrence'
  
  
  nf <- layout(matrix(c(1:15),nrow=5,ncol=3,byrow=TRUE),  TRUE) 
  layout.show(nf)
  par(mar=c(0.0,0.1,0.1,0.1),oma=c(3,3,0,1),mai=c(0.0,0.05,0.05,0) )

  # Compound-specific plots for Each Land Use
  
  Ymin <- 0
  Ymax <- 1
  
  x <- seq(min(unique(dat.ovary$Ag)), max(unique(dat.ovary$Ag)), length=100)
  
  y <- runif(100, min = 0, max = 1)
  
  ylabnums <- c(0.0,0.2,0.4,0.8,1.0)
  
  for(i in 1:15){
    
    plot(x,y, ylim=c(Ymin, Ymax), xlim=c(min(x), max(x)), axes=F, ylab='', xlab='', type='n')
    
    if( i <= 12){
      axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, labels=F ) 
    } else {
      axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01)
    }	
    
    if( i ==1 | i ==4 | i ==7 | i ==10 | i ==13){
      axis(side=2,cex.axis=axissize , mgp=c(0,0.3,0),tck= -0.01, las=1,at=format(pretty(ylabnums), digits=2), 
           labels=format(pretty(ylabnums), digits=2))
    } else {
      axis(side=2,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, labels=F)
    }	
  
  #axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01)
  axis(side=2,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, las=1, labels=F)
  # 
  # 
  i.for <- order(x)
  i.back <- order(x, decreasing = TRUE )
  x.polygon <- c(x[i.for] , x[i.back] )
  Lower <- lowerCIA.PopAve.Ag[i,]
  Upper <- upperCI.PopAve.Ag[i,]
  y.polygon <- c( Lower[i.for] , Upper[i.back] )
  polygon( x.polygon , y.polygon , col = "lightgray" , border = NA )
  # 
  lines(x,meanProbPopAve.Ag[i,], lwd = 1, col="black", lty = 1)
  # 
  text(-2, 0.9, i, cex=1)
  #text(0.6,0.9,'(7021)',cex=0.8)
  # 
  mtext(y.label, line = 1.5, side = 2, cex = size.text,outer=T)
  mtext(x.label, line = 1.5, side = 1, cex = size.text, outer=T)
  # 
  box()
  }
  par(def.par)
  dev.off()
  
  ##############
  #### Dev #####
  ##############
  
  tiff("Occurence_Dev_Ovary.tiff", 
       height=6, width=6, units='in', res=450, compression='lzw')
  def.par <- par(no.readonly = TRUE) 		
  
  size.labels = 1
  size.text = 1
  axissize <- 0.8
  x.label = 'Standardized[logit(Developed Land)]'
  y.label = 'Predicted Probability of Occurrence'
  
  
  nf <- layout(matrix(c(1:15),nrow=5,ncol=3,byrow=TRUE),  TRUE) 
  layout.show(nf)
  par(mar=c(0.0,0.1,0.1,0.1),oma=c(3,3,0,1),mai=c(0.0,0.05,0.05,0) )
  
  # Compound-specific plots for Each Land Use
  
  Ymin <- 0
  Ymax <- 1
  
  x <- seq(min(unique(dat.ovary$Dev)), max(unique(dat.ovary$Dev)), length=100)
  
  y <- runif(100, min = 0, max = 1)
  
  ylabnums <- c(0.0,0.2,0.4,0.8,1.0)
  
  for(i in 1:15){
    
    plot(x,y, ylim=c(Ymin, Ymax), xlim=c(min(x), max(x)), axes=F, ylab='', xlab='', type='n')
    
    if( i <= 12){
      axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, labels=F ) 
    } else {
      axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01)
    }	
    
    if( i ==1 | i ==4 | i ==7 | i ==10 | i ==13){
      axis(side=2,cex.axis=axissize , mgp=c(0,0.3,0),tck= -0.01, las=1,at=format(pretty(ylabnums), digits=2), 
           labels=format(pretty(ylabnums), digits=2))
    } else {
      axis(side=2,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, labels=F)
    }	
    
    #axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01)
    axis(side=2,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, las=1, labels=F)
    # 
    # 
    i.for <- order(x)
    i.back <- order(x, decreasing = TRUE )
    x.polygon <- c(x[i.for] , x[i.back] )
    Lower <- lowerCIA.PopAve.Dev[i,]
    Upper <- upperCI.PopAve.Dev[i,]
    y.polygon <- c( Lower[i.for] , Upper[i.back] )
    polygon( x.polygon , y.polygon , col = "lightgray" , border = NA )
    # 
    lines(x,meanProbPopAve.Dev[i,], lwd = 1, col="black", lty = 1)
    # 
    text(-1, 0.9, i, cex=1)
    #text(0.6,0.9,'(7021)',cex=0.8)
    # 
    mtext(y.label, line = 1.5, side = 2, cex = size.text,outer=T)
    mtext(x.label, line = 1.5, side = 1, cex = size.text, outer=T)
    # 
    box()
  }
  par(def.par)
  dev.off()
  
  ##############
  #### BMP #####
  ##############
  
  tiff("Occurence_BMP_Ovary.tiff", 
       height=6, width=6, units='in', res=450, compression='lzw')
  def.par <- par(no.readonly = TRUE) 		
  
  size.labels = 1
  size.text = 1
  axissize <- 0.8
  x.label = 'Standardized[BMP]'
  y.label = 'Predicted Probability of Occurrence'
  
  
  nf <- layout(matrix(c(1:15),nrow=5,ncol=3,byrow=TRUE),  TRUE) 
  layout.show(nf)
  par(mar=c(0.0,0.1,0.1,0.1),oma=c(3,3,0,1),mai=c(0.0,0.05,0.05,0) )
  
  # Compound-specific plots for Each Land Use
  
  Ymin <- 0
  Ymax <- 1
  
  x <- seq(min(unique(dat.ovary$BMP)), max(unique(dat.ovary$BMP)), length=100)
  
  y <- runif(100, min = 0, max = 1)
  
  ylabnums <- c(0.0,0.2,0.4,0.8,1.0)
  
  for(i in 1:15){
    
    plot(x,y, ylim=c(Ymin, Ymax), xlim=c(min(x), max(x)), axes=F, ylab='', xlab='', type='n')
    
    if( i <= 12){
      axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, labels=F ) 
    } else {
      axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01)
    }	
    
    if( i ==1 | i ==4 | i ==7 | i ==10 | i ==13){
      axis(side=2,cex.axis=axissize , mgp=c(0,0.3,0),tck= -0.01, las=1,at=format(pretty(ylabnums), digits=2), 
           labels=format(pretty(ylabnums), digits=2))
    } else {
      axis(side=2,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, labels=F)
    }	
    
    #axis(side=1,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01)
    axis(side=2,cex.axis=axissize , mgp=c(1,0,0),tck= -0.01, las=1, labels=F)
    # 
    # 
    i.for <- order(x)
    i.back <- order(x, decreasing = TRUE )
    x.polygon <- c(x[i.for] , x[i.back] )
    Lower <- lowerCIA.PopAve.BMP[i,]
    Upper <- upperCI.PopAve.BMP[i,]
    y.polygon <- c( Lower[i.for] , Upper[i.back] )
    polygon( x.polygon , y.polygon , col = "lightgray" , border = NA )
    # 
    lines(x,meanProbPopAve.BMP[i,], lwd = 1, col="black", lty = 1)
    # 
    text(4, 0.9, i, cex=1)
    #text(0.6,0.9,'(7021)',cex=0.8)
    # 
    mtext(y.label, line = 1.5, side = 2, cex = size.text,outer=T)
    mtext(x.label, line = 1.5, side = 1, cex = size.text, outer=T)
    # 
    box()
  }
  par(def.par)
  dev.off()
  
  