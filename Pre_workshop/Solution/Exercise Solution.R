# The calculation assume that the system is at equilibrium. Mean at-age parameters are provided.

rm(list=ls())
library(ggplot2)
source("functions.R") # contains functions "survivorship_F" and "MSYcalc"

# Read in the data sets

D <- read.csv("Pre_workshop/inputs.csv")

# at-age variables with 11+ group

w_age <- D$Weight*1000 # weight-at-age in g; ages 1 to 11+
m_age <- D$Maturity # maturity; ages 1 to 11+
v_age <- D$Vulnerability # vulnerability or selectivity; ages 1 to 11+

# other Parameters

M <- D$M[1] # assumed natural mortality rate
h <- D$h[1] # assumed steepness of the Beverton-Holt stock recruitment relationship
R0 <- D$R0[1] # model estimated unfished equilibrium recruitment (units = billions of recruits)

########################################################################################################################
# Move h/R0/phi0 parameterization of B-H SRR to a/b parameterization 
########################################################################################################################

# survivorship_F: 
# function from functions.R (calculates survivorship-at-age from f=F, M=M, waa=weight-at-age, mat=maturity-at-age, sel=vulnerability-at-age
# for unfished survivorship provide just M, waam mat 
# survivorship_F <- function(f=0,M,waa,mat,sel)
l_age <- survivorship_F(M=M,waa=w_age,mat=m_age) # unfished survivorship-at-age (F=0)
phi0 <- sum(l_age*w_age*m_age) # unfished spawning biomass per recruit

# B-H a and b
BHa <- 4*h/(phi0*(1-h)) # estimated Beverton-Holt a
BHb <- 1/R0*(BHa-1/phi0) # estimated Beverton-Holt b


########################################################################################################################
# Calculate equilibrium SSBmsy
########################################################################################################################

# MSYcalc: 
# function from functions.R returns a list with Fmsy, msy, SSBmsy from:
# M=M, waa=weight-at-age, mat=maturity-at-age, sel=vulnerability-at-age, Beverton-Holt a and b
# MSYcalc <- function(M,waa,mat,sel,a,b)

calc <- MSYcalc(M=M,waa=w_age,mat=m_age,sel=v_age,a=BHa,b=BHb) #biomass units (g X billion) = kt

########################################################################################################################
# Solutions
########################################################################################################################

#1. Equilibrium SSBMSY

  SSBMSY = calc$SSBmsy

#2. Equilibrium SSB0

  SSB0 = R0/(BHa-BHb*R0)

#3. Equilibrium SSB at F40% (assuming no stock recruitment relationship)

  # set up a data frame with various F, calculate biomass per recruit (phiF) and SPR for various F
  DF <- data.frame(f=seq(0.01,3,0.01),phiF=NA,SPR=NA)

  for(i in 1:nrow(DF)){
    l_f_age <- survivorship_F(f=DF$f[i],M=M,waa=w_age,mat=m_age,sel=v_age,message=F)
    DF$phiF[i] <- sum(l_f_age*w_age*m_age)
  }

  # Define SPR
  DF$SPR <- DF$phiF/phi0

  # Find F with SPR closest to 40%
  DF$abs_SPR_minus40 <- abs(DF$SPR - 0.4)
  f_at_SPR_40 <- DF$f[which(DF$abs_SPR_minus40 == min(DF$abs_SPR_minus40))]

  # Calculate Equilibrium SSB at F = "f_at_SPR_40" by multiplying biomass per recruit (phiF) by number of recruits
  Eq_SSB_at_SPR_40 = DF$phiF[DF$f==f_at_SPR_40] * R0

#4. SSB at 50% R0
  
  SSB_R50 = 0.5*R0/(BHa-BHb*0.5*R0)
