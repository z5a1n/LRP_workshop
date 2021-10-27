rm(list=ls())
library(ggplot2)
library(reshape2)
library(cowplot)
source("functions.R") # contains functions "survivorship_F" and "MSYcalc"

# Read in data

Data <- readRDS("4/Ex4_data.rda")

WAA <- Data$WAA # weight-at-age 1968-2020
MAT <- Data$MAT # maturity-at-age 1968-2020
VUL65 <- as.data.frame(Data$VUL65) # vulnerability-at-age 1968-2020 for model with h = 0.65
VUL95 <- as.data.frame(Data$VUL95) # vulnerability-at-age 1968-2020 for model with h = 0.95
D65 <- Data$D65 # Model estimated parameters by year for model with h = 0.65
D95 <- Data$D95 # Model estimated parameters by year for model with h = 0.95

VUL75 <- as.data.frame(Data$VUL75) # vulnerability-at-age 1968-2020 for model with h = 0.75
D75 <- Data$D75 # Model estimated parameters by year for model with h = 0.75

########################################################################################################################
# Move h/R0/phi0 parameterization of B-H SRR to a/b parameterization and calculate "unfished equilibrium SSB"
########################################################################################################################

#############
# h = 0.65
#############

# Model assumes a Beverton-Holt SRR with steepness h = 0.65
h_65 <- 0.65 # model assumed steepness
R0_65 <- 3.567227 # model estimated equilibrium unfished recruitment [provided]
M <- 0.35 # model assumed natural mortality rate

# mean unfished spawning biomass per recruit (phi0) over fist 5 years (5 years = mean generation time)
phi0_5_65<-mean(D65$phi0[1:5])

# B-H a and b. Here we assume that the stock recruitment relationship is estimate using phi0 from the first 5 years
BHa_65 <- 4*h_65/(phi0_5_65*(1-h_65)) # estimated Beverton-Holt a
BHb_65 <- 1/R0_65*(BHa_65-1/phi0_5_65) # estimated Beverton-Holt b

#"virgin" unfished SSB0
SSB0_65<- R0_65/(BHa_65-BHb_65*R0_65)

#############
# h = 0.95
#############

# Model assumes a Beverton-Holt SRR with steepness h = 0.95
h_95 <- 0.95 # model assumed steepness
R0_95 <- 3.377007 # model estimated equilibrium unfished recruitment [provided]
M <- 0.35 # model assumed natural mortality rate

# mean unfished spawning biomass per recruit (phi0) over fist 5 years (5 years = mean generation time)
phi0_5_95<-mean(D95$phi0[1:5])

# B-H a and b. Here we assume that the stock recruitment relationship is estimate using phi0 from the first 5 years
BHa_95 <- 4*h_95/(phi0_5_95*(1-h_95)) # estimated Beverton-Holt a
BHb_95 <- 1/R0_95*(BHa_95-1/phi0_5_95) # estimated Beverton-Holt b

#"virgin" unfished SSB0
SSB0_95<- R0_95/(BHa_95-BHb_95*R0_95)

#############
# h = 0.75
#############

# Model assumes a Beverton-Holt SRR with steepness h = 0.95
h_75 <- 0.75 # model assumed steepness
R0_75 <- 3.505041 # model estimated equilibrium unfished recruitment [provided]
M <- 0.35 # model assumed natural mortality rate

# mean unfished spawning biomass per recruit (phi0) over fist 5 years (5 years = mean generation time)
phi0_5_75<-mean(D75$phi0[1:5])

# B-H a and b. Here we assume that the stock recruitment relationship is estimate using phi0 from the first 5 years
BHa_75 <- 4*h_75/(phi0_5_75*(1-h_75)) # estimated Beverton-Holt a
BHb_75 <- 1/R0_75*(BHa_75-1/phi0_5_75) # estimated Beverton-Holt b

#"virgin" unfished SSB0
SSB0_75<- R0_75/(BHa_75-BHb_75*R0_75)

########################################################################################################################
# Calculate SSBmsy using WAA, MAT, VUL from first 5 years and last 10 years [Note: 53 years of data]
# Note: MSYcalc function can be used to estimate MSY reference points for various time periods. 
########################################################################################################################

# MSYcalc: function from source.R returns a list with Fmsy, msy, SSBmsy from M=M, waa=weight-at-age, mat=maturity-at-age, sel=vulnerability-at-age, Beverton-Holt a and b
# MSYcalc <- function(M,waa,mat,sel,a,b)

calc65_years1_5 <- MSYcalc(M=0.35,waa=apply(WAA[1:5,],2,mean), mat=apply(MAT[1:5,],2,mean), sel=apply(VUL65[1:5,],2,mean),a=BHa_65, b=BHb_65)
calc95_years1_5 <- MSYcalc(M=0.35,waa=apply(WAA[1:5,],2,mean), mat=apply(MAT[1:5,],2,mean), sel=apply(VUL95[1:5,],2,mean),a=BHa_95, b=BHb_95)
calc75_years1_5 <- MSYcalc(M=0.35,waa=apply(WAA[1:5,],2,mean), mat=apply(MAT[1:5,],2,mean), sel=apply(VUL75[1:5,],2,mean),a=BHa_75, b=BHb_75)

calc65_last10years <- MSYcalc(M=0.35,waa=apply(WAA[44:53,],2,mean), mat=apply(MAT[44:53,],2,mean), sel=apply(VUL65[44:53,],2,mean),a=BHa_65, b=BHb_65)
calc95_last10years <- MSYcalc(M=0.35,waa=apply(WAA[44:53,],2,mean), mat=apply(MAT[44:53,],2,mean), sel=apply(VUL95[44:53,],2,mean),a=BHa_95, b=BHb_95)
calc75_last10years <- MSYcalc(M=0.35,waa=apply(WAA[44:53,],2,mean), mat=apply(MAT[44:53,],2,mean), sel=apply(VUL75[44:53,],2,mean),a=BHa_75, b=BHb_75)

########################################################################################################################
# Plots 
########################################################################################################################

# Plots for Weight-at-age, Maturity-at-age, and vulnerability-at-age
ggplot(melt(cbind(WAA,Year=1968:2020),id.vars = "Year"),aes(x=Year,y=value,group=variable,color=variable)) + geom_path() + theme_classic() + labs(x="Year", y="Weight (kg)") + labs(color = "Age")
ggplot(melt(cbind(MAT,Year=1968:2020),id.vars = "Year"),aes(x=Year,y=value,group=variable,color=variable)) + geom_path() + theme_classic() + labs(x="Year", y="Maturity") + labs(color = "Age")
ggplot(melt(cbind(VUL65,Year=1968:2020),id.vars = "Year"),aes(x=Year,y=value,group=variable,color=variable)) + geom_path() + theme_classic() + labs(x="Year", y="Vulnerability (Model h = 0.65)") + labs(color = "Age")
ggplot(melt(cbind(VUL95,Year=1968:2020),id.vars = "Year"),aes(x=Year,y=value,group=variable,color=variable)) + geom_path() + theme_classic() + labs(x="Year", y="Vulnerability (Model h = 0.95)") + labs(color = "Age")

#plot SRRs and observations

p1 <- ggplot(D65[!is.na(D65$Rec),],aes(y=Rec,x=SSB,label=Year)) + geom_point() + theme_classic() + labs(x="SSB (kt)", y="Recruitment (10^9)") + expand_limits(y=13) + expand_limits(x=0) + ggtitle("h=0.65") + geom_function(fun=function(x) BHa_65*x/(1+BHb_65*x),color="red") 
p2 <- ggplot(D95[!is.na(D95$Rec),],aes(y=Rec,x=SSB,label=Year)) + geom_point() + theme_classic() + labs(x="SSB (kt)", y="Recruitment (10^9)") + expand_limits(y=13) + expand_limits(x=0) + ggtitle("h=0.95") + geom_function(fun=function(x) BHa_95*x/(1+BHb_95*x),color="blue") 
p3 <- ggplot(D75[!is.na(D95$Rec),],aes(y=Rec,x=SSB,label=Year)) + geom_point() + theme_classic() + labs(x="SSB (kt)", y="Recruitment (10^9)") + expand_limits(y=13) + expand_limits(x=0) + ggtitle("h=0.75") + geom_function(fun=function(x) BHa_75*x/(1+BHb_75*x)) 
p4 <- ggplot(D75[!is.na(D95$Rec),],aes(y=Rec,x=SSB,label=Year)) + theme_classic() + labs(x="SSB (kt)", y="Recruitment (10^9)") + expand_limits(y=13) + expand_limits(x=0) + ggtitle("h=0.65 red, h=0.95 blue, h=0.75 black") + 
  geom_function(fun=function(x) BHa_75*x/(1+BHb_75*x)) +
  geom_function(fun=function(x) BHa_65*x/(1+BHb_65*x),color="red") +
  geom_function(fun=function(x) BHa_95*x/(1+BHb_95*x),color="blue") 

plot_grid(p1,p2,p3,p4)

#Plot Historical SSB 
ggplot(D75,aes(y=SSB,x=Year)) + geom_path() + theme_classic() + labs(x="Year", y="SSB (kt)") + expand_limits(y=0) + ggtitle("h=0.65 red, h=0.95 blue, h=0.75 black") +
  geom_path(data=D65,mapping=aes(y=SSB,x=Year),color="red") +
  geom_path(data=D95,mapping=aes(y=SSB,x=Year),color="blue") 

#Plot Historical Catch
ggplot(D75,aes(y=Catch,x=Year)) + geom_path() + theme_classic() + labs(x="Year", y="Catch (kt)") + expand_limits(y=0)

#Plot F
ggplot(D75,aes(y=f,x=Year)) + geom_path() + theme_classic() + labs(x="Year", y="F") + expand_limits(y=0) + ggtitle("h=0.65 red, h=0.95 blue, h=0.75 black") +
  geom_path(data=D65,mapping=aes(y=f,x=Year),color="red") +
  geom_path(data=D95,mapping=aes(y=f,x=Year),color="blue") 

#Plot SSB0 (one dynamic option and virgin SSB0 plotted)
ggplot(D75,aes(y=SSB,x=Year)) + geom_path() + theme_classic() + labs(x="Year", y="SSB (kt)") + expand_limits(y=0) + ggtitle("h=0.65 red, h=0.95 blue, h=0.75 black") +
  geom_path(data=D65,mapping=aes(y=SSB,x=Year),color="red") +
  geom_path(data=D95,mapping=aes(y=SSB,x=Year),color="blue") +

  geom_path(data=D75, mapping=aes(y=dSSB0a,x=Year),linetype=3) +
  geom_path(data=D65, mapping=aes(y=dSSB0a,x=Year),color="red",linetype=3) +
  geom_path(data=D95, mapping=aes(y=dSSB0a,x=Year),color="blue",linetype=3) +
  
  geom_hline(yintercept=SSB0_75,linetype="dashed") +
  geom_hline(yintercept=SSB0_65,linetype="dashed",color="red") +
  geom_hline(yintercept=SSB0_95,linetype="dashed",color="blue") 

#Plot SSBmsy (assuming WAA, MAT, VUL from first 5 years)
ggplot(D75,aes(y=SSB,x=Year)) + geom_path() + theme_classic() + labs(x="Year", y="SSB (kt)") + expand_limits(y=0) + ggtitle("h=0.65 red, h=0.95 blue, h=0.75 black") +
  geom_path(data=D65,mapping=aes(y=SSB,x=Year),color="red") +
  geom_path(data=D95,mapping=aes(y=SSB,x=Year),color="blue") +
  geom_hline(yintercept=calc75_years1_5$SSBmsy,linetype="dashed") +
  geom_hline(yintercept=calc65_years1_5$SSBmsy,linetype="dashed",color="red") +
  geom_hline(yintercept=calc95_years1_5$SSBmsy,linetype="dashed",color="blue") 
  
#Plot SSBmsy (assuming WAA, MAT, VUL from last 10 years)
ggplot(D75,aes(y=SSB,x=Year)) + geom_path() + theme_classic() + labs(x="Year", y="SSB (kt)") + expand_limits(y=0) + ggtitle("h=0.65 red, h=0.95 blue, h=0.75 black") +
  geom_path(data=D65,mapping=aes(y=SSB,x=Year),color="red") +
  geom_path(data=D95,mapping=aes(y=SSB,x=Year),color="blue") +
  geom_hline(yintercept=calc75_last10years$SSBmsy,linetype="dashed") +
  geom_hline(yintercept=calc65_last10years$SSBmsy,linetype="dashed",color="red") +
  geom_hline(yintercept=calc95_last10years$SSBmsy,linetype="dashed",color="blue") 

# Plot acoustic index 1999-2020
  #Add x yr moving average (example 3 years)
  D75$MA_Index <- NA
  D75$MA_Index[D75$Year%in%2001:2020] <- apply(cbind(D75$Acoustic_Index[D75$Year%in%2001:2020],D75$Acoustic_Index[D75$Year%in%2000:2019],D75$Acoustic_Index[D75$Year%in%1999:2018]),1,mean)
  
  #Add loess smoother (example span = 0.5)
  lsmooth <- loess(Acoustic_Index ~ Year,data=D75,span=0.5)  
  D75$lowess_Index <- NA
  D75$lowess_Index[D75$Year%in%1999:2020] <- predict(lsmooth)
    
ggplot(D75[!is.na(D75$Acoustic_Index),]) + geom_path(mapping=aes(y=Acoustic_Index,x=Year)) + theme_classic() + labs(x="Year", y="Acoustic SSB (kt)") + expand_limits(y=0,x=1968) +
  geom_path(data=D75[!is.na(D75$MA_Index),],mapping=aes(y=MA_Index,x=Year),color="red") +
  geom_path(data=D75[!is.na(D75$lowess_Index),],mapping=aes(y=lowess_Index,x=Year),color="blue") 
  




