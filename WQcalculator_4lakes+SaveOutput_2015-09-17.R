# Calculate total loads and water quality for four lakes given direct drainage loads
# SRC 2015-09-10

rm(list = ls())
graphics.off()

# Functions +++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Estimate Chl from TP using the Filstrup model
# See Filstrup_Demo+MendotaData_2015-02-08.R for validation
# Paper is:
# Filstrup, C.T., T. Wagner, P.A. Soranno, E.H. Stanley,
# C.A. Stow, K.E. Webster and J.A. Downing. 2014.
# Regional variability among nonlinear chlorophyll-phosphorus
# relationships in lakes. Limnol. Oceanogr.59: 1691-1703.
filstrup = function(TP) {
  # Parameters from Table 3 of Filstrup et al.
  A = 2.14 # upper asymptote
  B = 0.49 # inflection point
  C = 3.21 # rate of increase
  D = 0.08 # lower asymptote
  TP.GM = 1.1287 # grand mean log10(TP) from Filstrup 2015-02-08
  lTP = log10(1000*TP) # log10 of TP in ug/L
  lChl = D + ( (A-D)/(1 + exp(-C*(lTP-TP.GM-B))) )  # Equation 1
  Chl = 10^lChl  # Chl in ug/L
  return(Chl)
}

# Forward transform
# Parameter must lie between 0 and Pmax; transform to unbounded parameter for fitting
FT = function(P,Pmax) {
  Ps = tan( (pi*P/Pmax) - (pi/2))
  return(Ps)
}

# Back transform
# Unbounded parameter to a value between 0 and Pmax
BT = function(Ps,Pmax) {
  P = (Pmax/pi)*( (pi/2) + atan(Ps))
}

# Function to Simulate P and export through one year
# This version does book-keeping on P0 each year, rather
# than starting over at a new value each year.
Ppred.v2 = function(par,P0,L,Pmax,dt,Nstep,NY) {
  S = BT(par[1],Pmax)
  H = BT(par[2],Pmax)
  W = BT(par[3],Pmax)
  P1hat = rep(0,NY)
  Ehat = rep(0,NY)
  # Compute estimated P1 and Export for each year
  for(iy in 1:NY) {
    Pstart = ifelse(iy==1,P0,P1hat[iy-1])
    Pt = rep(Pstart,Nstep)
    Et = rep(0,Nstep)
    # Integrate by Euler
    for(i in 2:Nstep){
      #Xt = H*(Pt[i-1] + W*L[iy])
      Xt = H*Pt[i-1] + W*L[iy]
      Et[i] = Et[i-1] + Xt*dt
      Pt[i] = Pt[i-1] + (L[iy]-S*Pt[i-1]- Xt)*dt
    } # end integration loop within one year
    P1hat[iy] = Pt[Nstep]
    Ehat[iy] = Et[Nstep]
  }  # end loop over years
  outlist = list(P1hat,Ehat)
  return(outlist)
}

# END OF FUNCTIONS ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Constants for all scenarios 

# Set Daphnia pulicaria to absent (0) or present (1) for Mendota and Monona
# Notes: (1) Presence/absence of Daphnia pulicaria will depend on time and WSC scenario
#        (2) Effect of Daphnia pulicaria has not been measured in Waubesa & Kegonsa
Dpul.Me = 0 # Daphnia pulicaria in Mendota
Dpul.Mo = 0 # Daphnia pulicaria in Monona

# Toggle stochastic within-year effects off (0) or on (1)
stoch = 0

# Constants for within-year simulation of P dynamics
Pmax = 1 # maximum value of balance model parameters, for the transform function
Nstep = 30 # number of time steps per year
dt = 1/Nstep  # time step for within-year dynamics


# <<<<<<<<<<<<<<<<
# Read the direct drainage loads from THMB
# Variable names are scenario code dot lake

cname=c('year','AR.Me','AR.Mo','AR.Wa','AR.Ke','AI.Me','AI.Mo','AI.Wa','AI.Ke',
        'NW.Me','NW.Mo','NW.Wa','NW.Ke','CC.Me','CC.Mo','CC.Wa','CC.Ke')
Raw = read.csv('LakeP.csv',header=F,col.names=cname)

Year = Raw$year
NY = length(Year) # number of years to simulate

DDloads = as.matrix(Raw[2:17])  # matrix of loads
Ploads = DDloads  # matrix to hold the total P loads (for Mendota only, same as DD loads)

# >>>>>>>>>>>>>>>>

# Load the data files needed to complete the calculations for the four lakes ****************

# Load Water Quality Regressions
# Predictors for Log summer TP:
#   Mendota input and Dpul
#   Monona output and Dpul
#   Waubesa and Kegonsa: outputs
# Predictors for Log summer Secchi transparency:
#   Mendota input and Dpul
#   Monona output and Dpul
#   Waubesa and Kegonsa P1 (NOvember 1 P after the growing season)
#save(LTPreg.Me,LTPreg.Mo,LTPreg.Wa,LTPreg.Ke,
#     LSecreg.Me,LSecreg.Mo,LSecreg.Wa,LSecreg.Ke,
#     file = 'WaterQualityRegressions.Rdata')
load(file = 'WaterQualityRegressions.Rdata')

# Parameters from 'BalanceModelFit_Mendota+2sigmas.R'
# par.t.Me is the transformed parameters needed to simulate the model
# P0bar.Me is the mean P0
# Shat.Me is the sedimentation coefficient
# Hhat.Me is the hydrologic outflow coefficient for P0
# What.Me is the hydrologic outflow coefficient for L
# sigmas are model errors for P and export respectively
# Save line:
#save(par.t.Me,P0bar.Me,Shat.Me,Hhat.Me,What.Me,sigmaP.Me,sigmaX.Me, file='Pbalance.Me.Rdata')
load(file='Pbalance.Me.Rdata')

# Parameters from 'BalanceModelFit_Monona+2sigmas.R'
# par.t.Mo is the transformed parameters needed to simulate the model
# P0bar.Mo is the mean P0
# Shat.Mo is the sedimentation coefficient
# Hhat.Mo is the hydrologic outflow coefficient for P0
# What.Mo is the hydrologic outflow coefficient for L
# sigmas are model errors for P and export respectively
# Save line:
#save(par.t.Mo,P0bar.Mo,Shat.Mo,Hhat.Mo,What.Mo,sigmaP.Mo,sigmaX.Mo, file='Pbalance.Mo.Rdata')
load(file='Pbalance.Mo.Rdata')

# Parameters from 'BalanceModelFit_Waubesa+1s+2sigmas.R'
# par.t.Wa is the transformed parameters needed to simulate the model
# P0bar.Wa is the mean P0
# Shat.Wa is the sedimentation coefficient
# Hhat.Wa is the hydrologic outflow coefficient for P0
# What.Wa is the hydrologic outflow coefficient for L
# sigmas are model errors for P and export respectively
# Save line:
#save(par.t.Wa,P0bar.Wa,Shat.Wa,Hhat.Wa,What.Wa,sigmaP.Wa,sigmaX.Wa, file='Pbalance.Wa.Rdata')
load(file='Pbalance.Wa.Rdata')

# Parameters from 'BalanceModelFit_Kegonsa+1s+2sigmas.R'
# par.t.Ke is the transformed parameters needed to simulate the model
# P0bar.Ke is the mean P0
# Shat.Ke is the sedimentation coefficient
# Hhat.Ke is the hydrologic outflow coefficient for P0
# What.Ke is the hydrologic outflow coefficient for L
# sigmas are model errors for P and export respectively
# Save line:
#save(par.t.Ke,P0bar.Ke,Shat.Ke,Hhat.Ke,What.Ke,sigmaP.Ke,sigmaX.Ke, file='Pbalance.Ke.Rdata')
load(file='Pbalance.Ke.Rdata')

# End loading of the additional data files ********************************************************

# Start loop over all 16 simulations  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# DDloads is NY*16 matrix of direct drainage loads
# Scenario sequence AR, AI, NW, CC; lake sequence Mo, Me, Wa, Ke
# Names in simulation code:
# AnnLoad.Me is load to Mendota (both total and direct drainage)
# DDload.Mo is direct drainage load to Monona
# DDload.Wa is direct drainage load to Waubesa
# DDload.Ke is direct drainage load to Kegonsa

# Matrices for output
TPmat = matrix(0,nr=NY,nc=16)
CHLmat = matrix(0,nr=NY,nc=16)
SECmat = matrix(0,nr=NY,nc=16)

for(iscen in 1:4)  { # start loop over all 4 scenarios
  # lake indices for output matrix columns
  iMe = 4*(iscen-1)+1
  iMo = 4*(iscen-1)+2
  iWa = 4*(iscen-1)+3
  iKe = 4*(iscen-1)+4
  # Direct drainage loads for scenario
  AnnLoad.Me = DDloads[,iMe]
  DDload.Mo = DDloads[,iMo]
  DDload.Wa = DDloads[,iWa]
  DDload.Ke = DDloads[,iKe]

# Compute Mendota Budget and WQ *******************************************************************

# Simulate within-year dynamics to get Nov P and Export
Sim.Me = Ppred.v2(par.t.Me,P0bar.Me,AnnLoad.Me,Pmax,dt,Nstep,NY)
PNovsim.Me = Sim.Me[[1]]
XPsim.Me = Sim.Me[[2]]

# Compute summer TP
LTPb0 = LTPreg.Me$coefficients[1] # intercept
LTPb1 = LTPreg.Me$coefficients[2] # effect of input
LTPb2 = LTPreg.Me$coefficients[3] # effect of D. pul.
LTPdf = LTPreg.Me$df.residual
LTPse = sd(LTPreg.Me$residuals,na.rm=T)/sqrt(LTPdf)
noise = stoch*LTPse*rt(NY,df=LTPdf)
TPsummer.Me = exp( LTPb0 + LTPb1*AnnLoad.Me + LTPb2*Dpul.Me + noise )
TPmat[,iMe] = TPsummer.Me

# Compute summer Chl
Chlsummer.Me = filstrup(TPsummer.Me)
CHLmat[,iMe] = Chlsummer.Me

# Compute summer Secchi transparency
LSecb0 = LSecreg.Me$coefficients[1] # intercept
LSecb1 = LSecreg.Me$coefficients[2] # effect of input
LSecb2 = LSecreg.Me$coefficients[3] # effect of D. pul.
LSecdf = LSecreg.Me$df.residual
LSecse = sd(LSecreg.Me$residuals,na.rm=T)/sqrt(LSecdf)
noise = stoch*LSecse*rt(NY,df=LSecdf)
Secsummer.Me = exp( LSecb0 + LSecb1*AnnLoad.Me + LSecb2*Dpul.Me + noise)
SECmat[,iMe] = Secsummer.Me

# End Mendota, Compute Monona Budget and WQ ***********************************

# Annual load is Direct Drainage plus river transfer
AnnLoad.Mo = DDload.Mo + XPsim.Me

# Run within-year simulations to get P in Nov and Export
# Simulate within-year dynamics to get Nov P and Export
Sim.Mo = Ppred.v2(par.t.Mo,P0bar.Mo,AnnLoad.Mo,Pmax,dt,Nstep,NY)
PNovsim.Mo = Sim.Mo[[1]]
XPsim.Mo = Sim.Mo[[2]]

# Calculate summer TP using regression
LTPb0 = LTPreg.Mo$coefficients[1] # intercept
LTPb1 = LTPreg.Mo$coefficients[2] # effect of output
LTPb2 = LTPreg.Mo$coefficients[3] # effect of D. pul.
LTPdf = LTPreg.Mo$df.residual
LTPse = sd(LTPreg.Mo$residuals,na.rm=T)/sqrt(LTPdf)
noise = stoch*LTPse*rt(NY,df=LTPdf)
TPsummer.Mo = exp( LTPb0 + LTPb1*XPsim.Mo + LTPb2*Dpul.Mo + noise )
TPmat[,iMo] =TPsummer.Mo

# Compute summer Chl
Chlsummer.Mo = filstrup(TPsummer.Mo)
CHLmat[,iMo] = Chlsummer.Mo

# Compute summer Secchi transparency
LSecb0 = LSecreg.Mo$coefficients[1] # intercept
LSecb1 = LSecreg.Mo$coefficients[2] # effect of input
LSecb2 = LSecreg.Mo$coefficients[3] # effect of D. pul.
LSecdf = LSecreg.Mo$df.residual
LSecse = sd(LSecreg.Mo$residuals,na.rm=T)/sqrt(LSecdf)
noise = stoch*LSecse*rt(NY,df=LSecdf)
Secsummer.Mo = exp( LSecb0 + LSecb1*XPsim.Mo + LSecb2*Dpul.Mo + noise)
SECmat[,iMo] = Secsummer.Mo

# End Monona, Compute Waubesa Budget and WQ ***********************************

# Annual load is Direct Drainage plus river transfer
AnnLoad.Wa = DDload.Wa + XPsim.Mo

# Run within-year simulations to get P in Nov and Export
# Simulate within-year dynamics to get Nov P and Export
Sim.Wa = Ppred.v2(par.t.Wa,P0bar.Wa,AnnLoad.Wa,Pmax,dt,Nstep,NY)
PNovsim.Wa = Sim.Wa[[1]]
XPsim.Wa = Sim.Wa[[2]]

# Calculate summer TP using regression
LTPb0 = LTPreg.Wa$coefficients[1] # intercept
LTPb1 = LTPreg.Wa$coefficients[2] # effect of output
LTPdf = LTPreg.Wa$df.residual
LTPse = sd(LTPreg.Wa$residuals,na.rm=T)/sqrt(LTPdf)
noise = stoch*LTPse*rt(NY,df=LTPdf)
TPsummer.Wa = exp( LTPb0 + LTPb1*XPsim.Wa + noise)
TPmat[,iWa] = TPsummer.Wa

# Compute summer Chl
Chlsummer.Wa = filstrup(TPsummer.Wa)
CHLmat[,iWa] = Chlsummer.Wa

# Calculate summer Secchi transparency
LSecb0 = LSecreg.Wa$coefficients[1] # intercept
LSecb1 = LSecreg.Wa$coefficients[2] # effect of P1 (Nov. P mass)
LSecdf = LSecreg.Wa$df.residual
LSecse = sd(LSecreg.Wa$residuals,na.rm=T)/sqrt(LSecdf)
noise = stoch*LSecse*rt(NY,df=LSecdf)
Secsummer.Wa = exp( LSecb0 + LSecb1*PNovsim.Wa + noise)
SECmat[,iWa] = Secsummer.Wa

# End Waubesa, Compute Kegonsa Budget and TP ***********************************

# Annual load is Direct Drainage plus river transfer
AnnLoad.Ke = DDload.Ke + XPsim.Wa

# Run within-year simulations to get P in Nov and Export
# Simulate within-year dynamics to get Nov P and Export
Sim.Ke = Ppred.v2(par.t.Ke,P0bar.Ke,AnnLoad.Ke,Pmax,dt,Nstep,NY)
PNovsim.Ke = Sim.Ke[[1]]
XPsim.Ke = Sim.Ke[[2]]

# Calculate summer TP using regression
LTPb0 = LTPreg.Ke$coefficients[1]
LTPb1 = LTPreg.Ke$coefficients[2]
LTPdf = LTPreg.Ke$df.residual
LTPse = sd(LTPreg.Ke$residuals,na.rm=T)/sqrt(LTPdf)
noise = stoch*LTPse*rt(NY,df=LTPdf)
TPsummer.Ke = exp( LTPb0 + LTPb1*XPsim.Ke + noise)
TPmat[,iKe] = TPsummer.Ke

# Compute summer Chl
Chlsummer.Ke = filstrup(TPsummer.Ke)
CHLmat[,iKe] = Chlsummer.Ke

# Calculate summer Secchi transparency
LSecb0 = LSecreg.Ke$coefficients[1] # intercept
LSecb1 = LSecreg.Ke$coefficients[2] # effect of P1 (Nov. P mass)
LSecdf = LSecreg.Ke$df.residual
LSecse = sd(LSecreg.Ke$residuals,na.rm=T)/sqrt(LSecdf)
noise = stoch*LSecse*rt(NY,df=LSecdf)
Secsummer.Ke = exp( LSecb0 + LSecb1*PNovsim.Ke + noise)
SECmat[,iKe] = Secsummer.Ke

# END COMPUTATION OF BUDGETS AND WATER QUALITY FOR THE FOUR LAKES ************************
} # End loop over scenarios >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


# <<<<<<<<<<<<<<<<<<<<<<<<<

# Generate plots
# Color scheme:
# Nested Watershed = red
# Accelerated Innovation = gold or tan
# Connected communities = blue
# Abandonment & renewal = green

# Smoothed plots
library('zoo')

# Plots for TP

lTPmat = log(TPmat)
lTProll = rollapply(lTPmat,width=9,FUN=mean,na.rm=T,na.pad=T,by.column=T)
Yroll = exp(lTProll)

windows()
par(mfrow=c(2,2),cex.axis=1.4,cex.lab=1.4,cex.main=1.4,font.main=1,mar=c(4, 4.2, 2.5, 2) + 0.1)
yrange = range(Yroll,na.rm=T)
plot(Year,Yroll[,1],type='l',lwd=2,col='forestgreen',ylim=yrange,log='y',
     xlab='Year',ylab='Summer [TP], mg/L',main='Mendota')
points(Year,Yroll[,5],type='l',lwd=2,col='darkgoldenrod1')
points(Year,Yroll[,9],type='l',lwd=2,col='red')
points(Year,Yroll[,13],type='l',lwd=2,col='blue')

#ymin = yrange[1]
#yrange = c(ymin,0.1)
plot(Year,Yroll[,2],type='l',lwd=2,col='forestgreen',ylim=yrange,log='y',
     xlab='Year',ylab='Summer [TP], mg/L',main='Monona')
points(Year,Yroll[,6],type='l',lwd=2,col='darkgoldenrod1')
points(Year,Yroll[,10],type='l',lwd=2,col='red')
points(Year,Yroll[,14],type='l',lwd=2,col='blue')

legend('topleft',
       legend=c('AI','AR','CC','NW'),bty='n',text.font=2,
       lty=c(1,1,1,1),lwd=c(2,2,2,2),seg.len=1,
       col=c('darkgoldenrod1','forestgreen','blue','red') )

#yrange = c(ymin,0.2)
plot(Year,Yroll[,3],type='l',lwd=2,col='forestgreen',ylim=yrange,log='y',
     xlab='Year',ylab='Summer [TP], mg/L',main='Waubesa')
points(Year,Yroll[,7],type='l',lwd=2,col='darkgoldenrod1')
points(Year,Yroll[,11],type='l',lwd=2,col='red')
points(Year,Yroll[,15],type='l',lwd=2,col='blue')

#yrange = c(ymin,0.25)
plot(Year,Yroll[,4],type='l',lwd=2,col='forestgreen',ylim=yrange,log='y',
     xlab='Year',ylab='Summer [TP], mg/L',main='Kegonsa')
points(Year,Yroll[,8],type='l',lwd=2,col='darkgoldenrod1')
points(Year,Yroll[,12],type='l',lwd=2,col='red')
points(Year,Yroll[,16],type='l',lwd=2,col='blue')

# Plots for Chlorophyll

# smooth in log units
#lCHLmat = log(CHLmat)
#lCHLroll = rollapply(lCHLmat,width=9,FUN=mean,na.rm=T,na.pad=T,by.column=T)
#Yroll = exp(lCHLroll)

# smooth in natural units
Yroll = rollapply(CHLmat,width=9,FUN=mean,na.rm=T,na.pad=T,by.column=T) # DUMMY for testing

windows()
par(mfrow=c(2,2),cex.axis=1.4,cex.lab=1.4,cex.main=1.4,font.main=1,mar=c(4, 4.2, 2.5, 2) + 0.1)
yrange = range(Yroll,na.rm=T)
plot(Year,Yroll[,1],type='l',lwd=2,col='forestgreen',ylim=yrange,log='y',
     xlab='Year',ylab='Summer [Chl], ug/L',main='Mendota')
points(Year,Yroll[,5],type='l',lwd=2,col='darkgoldenrod1')
points(Year,Yroll[,9],type='l',lwd=2,col='red')
points(Year,Yroll[,13],type='l',lwd=2,col='blue')

#ymin = yrange[1]
#yrange = c(ymin,0.1)
plot(Year,Yroll[,2],type='l',lwd=2,col='forestgreen',ylim=yrange,log='y',
     xlab='Year',ylab='Summer [Chl], ug/L',main='Monona')
points(Year,Yroll[,6],type='l',lwd=2,col='darkgoldenrod1')
points(Year,Yroll[,10],type='l',lwd=2,col='red')
points(Year,Yroll[,14],type='l',lwd=2,col='blue')

legend('topleft',
       legend=c('AI','AR','CC','NW'),bty='n',text.font=2,
       lty=c(1,1,1,1),lwd=c(2,2,2,2),seg.len=1,
       col=c('darkgoldenrod1','forestgreen','blue','red') )

#yrange = c(ymin,0.2)
plot(Year,Yroll[,3],type='l',lwd=2,col='forestgreen',ylim=yrange,log='y',
     xlab='Year',ylab='Summer [Chl], ug/L',main='Waubesa')
points(Year,Yroll[,7],type='l',lwd=2,col='darkgoldenrod1')
points(Year,Yroll[,11],type='l',lwd=2,col='red')
points(Year,Yroll[,15],type='l',lwd=2,col='blue')

#yrange = c(ymin,0.25)
plot(Year,Yroll[,4],type='l',lwd=2,col='forestgreen',ylim=yrange,log='y',
     xlab='Year',ylab='Summer [Chl], ug/L',main='Kegonsa')
points(Year,Yroll[,8],type='l',lwd=2,col='darkgoldenrod1')
points(Year,Yroll[,12],type='l',lwd=2,col='red')
points(Year,Yroll[,16],type='l',lwd=2,col='blue')

# Plots for Secchi

# smooth in natural units
Yroll = rollapply(SECmat,width=9,FUN=mean,na.rm=T,na.pad=T,by.column=T) # DUMMY for testing

windows()
par(mfrow=c(2,2),cex.axis=1.4,cex.lab=1.4,cex.main=1.4,font.main=1,mar=c(4, 4.2, 2.5, 2) + 0.1)
yrange = range(Yroll[,1],Yroll[,5],Yroll[,9],Yroll[,13],na.rm=T)
plot(Year,Yroll[,1],type='l',lwd=2,col='forestgreen',ylim=yrange,log='y',
     xlab='Year',ylab='Summer Secchi, m',main='Mendota')
points(Year,Yroll[,5],type='l',lwd=2,col='darkgoldenrod1')
points(Year,Yroll[,9],type='l',lwd=2,col='red')
points(Year,Yroll[,13],type='l',lwd=2,col='blue')

legend('bottomleft',
       legend=c('AI','AR','CC','NW'),bty='n',text.font=2,
       lty=c(1,1,1,1),lwd=c(2,2,2,2),seg.len=1,
       col=c('darkgoldenrod1','forestgreen','blue','red') )

yrange = range(Yroll[,2],Yroll[,6],Yroll[,10],Yroll[,14],na.rm=T)
plot(Year,Yroll[,2],type='l',lwd=2,col='forestgreen',ylim=yrange,log='y',
     xlab='Year',ylab='Summer Secchi, m',main='Monona')
points(Year,Yroll[,6],type='l',lwd=2,col='darkgoldenrod1')
points(Year,Yroll[,10],type='l',lwd=2,col='red')
points(Year,Yroll[,14],type='l',lwd=2,col='blue')

yrange = range(Yroll[,3],Yroll[,7],Yroll[,11],Yroll[,15],na.rm=T)
plot(Year,Yroll[,3],type='l',lwd=2,col='forestgreen',ylim=yrange,log='y',
     xlab='Year',ylab='Summer Secchi, m',main='Waubesa')
points(Year,Yroll[,7],type='l',lwd=2,col='darkgoldenrod1')
points(Year,Yroll[,11],type='l',lwd=2,col='red')
points(Year,Yroll[,15],type='l',lwd=2,col='blue')

yrange = range(Yroll[,4],Yroll[,8],Yroll[,12],Yroll[,16],na.rm=T)
plot(Year,Yroll[,4],type='l',lwd=2,col='forestgreen',ylim=yrange,log='y',
     xlab='Year',ylab='Summer Secchi, m',main='Kegonsa')
points(Year,Yroll[,8],type='l',lwd=2,col='darkgoldenrod1')
points(Year,Yroll[,12],type='l',lwd=2,col='red')
points(Year,Yroll[,16],type='l',lwd=2,col='blue')

# Save output
# NY is number of years, Year is 4-digit year
# DDloads is the annual direct drainage loads of P to the lakes
# CHLmat is summer Chl
# SECmat is summer Secchi
# TPmat is summer TP
save(NY,Year,DDloads,CHLmat,SECmat,TPmat,file='WQoutput.Rdata')
