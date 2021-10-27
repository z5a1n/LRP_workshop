rm(list=ls())
library(ggplot2)
source("functions.R") # contains functions "survivorship_F" and "MSYcalc"

# Read in the two data sets

AA <- read.csv("2/Ex2_At_Age_Data.csv")
D <- read.csv("2/Ex2_Data.csv")

########################################################################################################################
# Move h/R0/phi0 parameterization of B-H SRR to a/b parameterization and calculate "unfished equilibrium SSB"
########################################################################################################################

# Model assumes a Beverton-Holt SRR with steepness h = 0.75
h <- 0.75 # model assumed steepness
R0 <- 3.505041 # model estimated equilibrium unfished recruitment [provided]
M <- 0.35 # model assumed natural mortality rate

# survivorship_F: function from functions.R (calculates survivorship-at-age from f=F (default is F=0), M=M, n_ages = number of age classes (including 0 and plus group), sel=vulnerability-at-age when F !=0)
# survivorship_F <- function(f=0,M,n_ages,sel)
l_age <- survivorship_F(M=M,n_ages=length(AA$w_age)) # survivorship-at-age (unfished, F=0)
phi0 <- sum(l_age*AA$w_age*AA$m_age) # unfished spawning biomass  per recruit

# B-H a and b
BHa <- 4*h/(phi0*(1-h)) # estimated Beverton-Holt a
BHb <- 1/R0*(BHa-1/phi0) # estimated Beverton-Holt b

# Unfished equilibrum SSB0
SSB0 <- R0/(BHa-BHb*R0)

########################################################################################################################
# Calculate equilibrium SSBmsy
########################################################################################################################

# MSYcalc: function from functions.R returns a list with Fmsy, msy, SSBmsy from M=M, waa=weight-at-age, mat=maturity-at-age, sel=vulnerability-at-age, Beverton-Holt a and b
# MSYcalc <- function(M,waa,mat,sel,a,b)

calc <- MSYcalc(M=0.35,waa=AA$w_age,mat=AA$m_age,sel=AA$v_age,a=BHa,b=BHb)

Fmsy <- calc$Fmsy
msy <- calc$msy
SSBmsy <- calc$SSBmsy


########################################################################################################################
# Plots from data frame AA
########################################################################################################################

# Plots for Weight-at-age, Maturity-at-age, and vulnerability-at-age
ggplot(AA,aes(y=w_age/1000,x=age)) + geom_path() + theme_classic() + labs(x="Age", y="Weight (kg)") + expand_limits(y=0) + expand_limits(x=0) 
ggplot(AA,aes(y=m_age,x=age)) + geom_path() + theme_classic() + labs(x="Age", y="Maturity") + expand_limits(y=0) + expand_limits(x=0) 
ggplot(AA,aes(y=v_age,x=age)) + geom_path() + theme_classic() + labs(x="Age", y="Vulnerability") + expand_limits(y=0) + expand_limits(x=0) 

########################################################################################################################
# Plots from data frame D
########################################################################################################################
#plot SR observations
ggplot(D[!is.na(D$Rec),],aes(y=Rec,x=SSB,label=Year)) + geom_point() + theme_classic() + labs(x="SSB (kt)", y="Recruitment (10^9)") + expand_limits(y=0) + expand_limits(x=0) +
  geom_text(mapping=aes(y=Rec,x=SSB,label=Year),nudge_y = 0.5,size=2) 
#plot SRR and observations
ggplot() + geom_point(D[!is.na(D$Rec),],mapping=aes(y=Rec,x=SSB)) + theme_classic() + labs(x="SSB (kt)", y="Recruitment (10^9)") + expand_limits(y=0) + expand_limits(x=0) +
  geom_function(fun=function(x) BHa*x/(1+BHb*x)) 
#Plot Historical Recruitment
ggplot(D[!is.na(D$Rec),],aes(y=Rec,x=Year)) + geom_path() + theme_classic() + labs(x="Year", y="Recruitment (10^9)") + expand_limits(y=0) 
#Plot Historical SSB 
ggplot(D,aes(y=SSB,x=Year)) + geom_path() +theme_classic() + labs(x="Year", y="SSB (kt)") + expand_limits(y=0)
#Plot Historical Catch
ggplot(D,aes(y=Catch,x=Year)) + geom_path() + theme_classic() + labs(x="Year", y="Catch (kt)") + expand_limits(y=0)
#Plot F
ggplot(D,aes(y=f,x=Year)) + geom_path() + theme_classic() + labs(x="Year", y="F") + expand_limits(y=0) 

#Plot Equilibrium SSB0, Equilibrium SSBmsy
ggplot(D) + 
  geom_path(mapping=aes(y=SSB,x=Year)) +
  theme_classic() + labs(x="Year", y="SSB (kt)") + expand_limits(y=0) +
  geom_hline(yintercept=SSB0, linetype="dashed", color = "blue") + 
  geom_hline(yintercept=SSBmsy, color = "green") 


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

