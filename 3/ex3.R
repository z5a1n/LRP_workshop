rm(list=ls())
library(ggplot2)
library(reshape2)
source("functions.R") # contains functions "survivorship_F" and "MSYcalc"

# Read in data

Data <- readRDS("3/Ex3_data.rda")

WAA <- Data[[1]] # weight-at-age 1968-2020
MAT <- Data[[2]] # maturity-at-age 1968-2020
VUL <- Data[[3]] # vulnerability-at-age 1968-2020

D <- Data[[4]] # Model estimated parameters by year 

########################################################################################################################
# Move h/R0/phi0 parameterization of B-H SRR to a/b parameterization and calculate "unfished equilibrium SSB"
########################################################################################################################

# Model assumes a Beverton-Holt SRR with steepness h = 0.75
h <- 0.75 # model assumed steepness
R0 <- 3.505041 # model estimated equilibrium unfished recruitment [provided]
M <- 0.35 # model assumed natural mortality rate

# mean unfished spawning biomass per recruit (phi0) over fist 5 years (5 years = mean generation time)
phi0_5<-mean(D$phi0[1:5])

# B-H a and b. Here we assume that the stock recruitment relationship is estimate using phi0 from the first 5 years
BHa <- 4*h/(phi0_5*(1-h)) # estimated Beverton-Holt a
BHb <- 1/R0*(BHa-1/phi0_5) # estimated Beverton-Holt b

#"virgin" unfished SSB0
SSB0 <- R0/(BHa-BHb*R0)

########################################################################################################################
# Calculate SSBmsy using WAA, MAT, VUL from first 5 years and last 10 years [Note: 53 years of data]
# Note: MSYcalc function can be used to estimate MSY reference points for various time periods. 
########################################################################################################################

# MSYcalc: function from source.R returns a list with Fmsy, msy, SSBmsy from M=M, waa=weight-at-age, mat=maturity-at-age, sel=vulnerability-at-age, Beverton-Holt a and b
# MSYcalc <- function(M,waa,mat,sel,a,b)

calc_years1_5 <- MSYcalc(M=0.35,waa=apply(WAA[1:5,],2,mean), mat=apply(MAT[1:5,],2,mean), sel=apply(VUL[1:5,],2,mean),a=BHa, b=BHb)
calc_last10years <- MSYcalc(M=0.35,waa=apply(WAA[44:53,],2,mean), mat=apply(MAT[44:53,],2,mean), sel=apply(VUL[44:53,],2,mean),a=BHa, b=BHb)

calc_years1_5$Fmsy
calc_years1_5$SSBmsy
calc_years1_5$msy

########################################################################################################################
# Plots 
########################################################################################################################

# Plots for Weight-at-age, Maturity-at-age, and vulnerability-at-age
ggplot(melt(cbind(WAA,Year=1968:2020),id.vars = "Year"),aes(x=Year,y=value,group=variable,color=variable)) + geom_path() + theme_classic() + labs(x="Year", y="Weight (kg)") + labs(color = "Age")
ggplot(melt(cbind(MAT,Year=1968:2020),id.vars = "Year"),aes(x=Year,y=value,group=variable,color=variable)) + geom_path() + theme_classic() + labs(x="Year", y="Maturity") + labs(color = "Age")
ggplot(melt(cbind(VUL,Year=1968:2020),id.vars = "Year"),aes(x=Year,y=value,group=variable,color=variable)) + geom_path() + theme_classic() + labs(x="Year", y="Vulnerability") + labs(color = "Age")

########################################################################################################################
# Plots from data frame D
########################################################################################################################

ggplot(D[!is.na(D$Rec),],aes(y=Rec,x=SSB,label=Year)) + geom_point() +
  theme_classic() + labs(x="SSB (kt)", y="Recruitment (10^9)") + expand_limits(y=0) + expand_limits(x=0) +
  geom_text(mapping=aes(y=Rec,x=SSB,label=Year),nudge_y = 0.5,size=2) 

#plot SRR and observations
ggplot() + geom_point(D[!is.na(D$Rec),],mapping=aes(y=Rec,x=SSB)) +
  theme_classic() + labs(x="SSB (kt)", y="Recruitment (10^9)") + expand_limits(y=0) + expand_limits(x=0) +
  geom_function(fun=function(x) BHa*x/(1+BHb*x)) 

#Plot Historical Recruitment
ggplot(D[!is.na(D$Rec),],aes(y=Rec,x=Year)) + geom_path() + theme_classic() + labs(x="Year", y="Recruitment (10^9)") + expand_limits(y=0) 

#Plot Historical Recruitment Deviations
ggplot(D[!is.na(D$rec_dev),],aes(y=log(rec_dev),x=Year)) + geom_path() + theme_classic() + labs(x="Year", y="log(Recruitment Deviations)") + expand_limits(y=0) 

#Plot Historical SSB 
ggplot(D,aes(y=SSB,x=Year)) + geom_path() + theme_classic() + labs(x="Year", y="SSB (kt)") + expand_limits(y=0)

#Plot Historical Catch
ggplot(D,aes(y=Catch,x=Year)) + geom_path() + theme_classic() + labs(x="Year", y="Catch (kt)") + expand_limits(y=0)

#Plot F
ggplot(D,aes(y=f,x=Year)) + geom_path() + theme_classic() + labs(x="Year", y="F") + expand_limits(y=0) 

#Plot dynamic SSB0
ggplot(D) + 
  geom_path(mapping=aes(y=SSB,x=Year)) +
  geom_path(mapping=aes(y=dSSB0b,x=Year),color="blue") +
  geom_path(mapping=aes(y=dSSB0a,x=Year),color="red") +
  geom_path(mapping=aes(y=dSSB0g,x=Year),color="purple") +
  geom_path(mapping=aes(y=dSSB0m,x=Year),color="grey") +
  theme_classic() + labs(x="Year", y="SSB (kt)") + expand_limits(y=0) +
  geom_hline(yintercept=SSB0, linetype="dashed", color = "green") 

#Plot SSBmsy with WAA,MAT,VUL from different time periods
ggplot(D) + 
  geom_path(mapping=aes(y=SSB,x=Year)) +
  theme_classic() + labs(x="Year", y="SSB (kt)") + expand_limits(y=0) +
  geom_hline(yintercept=calc_years1_5$SSBmsy, linetype="dashed", color = "red") +
  geom_hline(yintercept=calc_last10years$SSBmsy, linetype="dashed", color = "blue") 
  
# Plot acoustic index 1999-2020
  #Add x yr moving average (example 3 years)
  D$MA_Index <- NA
  D$MA_Index[D$Year%in%2001:2020] <- apply(cbind(D$Acoustic_Index[D$Year%in%2001:2020],D$Acoustic_Index[D$Year%in%2000:2019],D$Acoustic_Index[D$Year%in%1999:2018]),1,mean)

  #Add loess smoother (example span = 0.5)
  lsmooth <- loess(Acoustic_Index ~ Year,data=D,span=0.5)  
  D$lowess_Index <- NA
  D$lowess_Index[D$Year%in%1999:2020] <- predict(lsmooth)
  
  ggplot(D[!is.na(D$Acoustic_Index),]) + geom_path(mapping=aes(y=Acoustic_Index,x=Year)) + theme_classic() + labs(x="Year", y="Acoustic SSB (kt)") + expand_limits(y=0,x=1968) +
    geom_path(data=D[!is.na(D$MA_Index),],mapping=aes(y=MA_Index,x=Year),color="red") +
    geom_path(data=D[!is.na(D$lowess_Index),],mapping=aes(y=lowess_Index,x=Year),color="blue") 
  
  