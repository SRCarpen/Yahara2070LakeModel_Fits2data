# Demonstration of Filstrup et al. model
# to predict Chl in ug/L from TP in ug/L
# SRC 2015-02-07

# Model is equation 1
# Our region is the Rock River region (RR)
# Parameter values were taken from Table 3 for RR
#
# The paper citation is:
# Filstrup, C.T., T. Wagner, P.A. Soranno, E.H. Stanley,
# C.A. Stow, K.E. Webster and J.A. Downing. 2014.
# Regional variability among nonlinear chlorophyll-phosphorus
# relationships in lakes. Limnol. Oceanogr.59: 1691-1703.

rm(list = ls())  # clear memory
graphics.off()  # clear graphics windows

# Look at Mendota Chl (despite warts)
# This dataset is from \Yahara\YLP_model5 in July 2010
# Read in the data set
# Save line from WQ-daily_JulAugSep_MakeData_V1.r
# Save data for model fitting
#save(P1Nov,P0Nov,LNov,XPNov,Secchi,Chl,TP,DRP,Dpul,
# file='Nov+JulyAugSepDailyWQ.Rdata')

load('Nov+JulyAugSepDailyWQ.Rdata')
Chl = ifelse(Chl==0,NA,Chl)

# Convert to ug/L, Rename Mendota Chl and TP, discard the rest
Chl.Me = Chl
TP.Me = TP*1000
rm(P1Nov,P0Nov,LNov,XPNov,Secchi,Chl,TP,DRP,Dpul)

# TP gradient (probably much wider than we will see in data or simulations)
NTP = 100
lTP = seq(-1,3,length.out=NTP) # gradient of log TP
TP = 10^lTP # gradient of TP in ug/L
lnTP = log(TP)  # natural log TP

# Parameters from Table 3 of Filstrup et al.
A = 2.14 # upper asymptote
B = 0.49 # inflection point
C = 3.21 # rate of increase
D = 0.08 # lower asymptote

# Note there is a typo 2 lines below eq. 2 in the Filstrup et al. paper
# The symbol 'TP' in the equations has to be log TP
# Otherwise the equation does not make any sense
# Correction confirmed by email with Christopher Filstrup on 8 Feb 2015
# In addition, Filstrup points out that TP is centered around
# the grand mean of log10(TP) which is 1.1287
TP.GM = 1.1287 # grand mean log10(TP) from Filstrup 2015-02-08
lChl = D + ( (A-D)/(1 + exp(-C*(lTP-TP.GM-B))) )  # Equation 1
Chl = 10^lChl  # Chl in ug/L

windows()
par(mfrow=c(1,1),mar=c(5, 4.3, 4, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
plot(TP,Chl,type='l',lwd=2,col='blue',log='xy',
     xlab = 'TP', ylab = 'Chl')
points(TP,TP,type='l',lty=2)

# Mendota data
windows()
par(mfrow=c(1,1),mar=c(5, 4.3, 4, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
plot(TP.Me,Chl.Me,type='p',pch=19,col='forestgreen',log='xy',
     xlab = 'TP', ylab = 'Chl')

# Overlay
windows()
par(mfrow=c(1,1),mar=c(5, 4.3, 4, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
plot(TP,Chl,type='l',lwd=2,col='blue',log='xy',xlim=c(5,500),
     xlab = 'TP', ylab = 'Chl',
     main='Filstrup et al. model blue line, Mendota data green points')
points(TP.Me,Chl.Me,type='p',pch=19,col='forestgreen')

