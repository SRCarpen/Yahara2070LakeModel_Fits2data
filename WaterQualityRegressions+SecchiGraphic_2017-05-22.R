# Compute regressions for water quality versus P budget components
# SRC 29 August 2014

rm(list = ls())
graphics.off()

# Load the complete database from PBudgets_LogMeanWQ_4Lakes_2013-01-26.r
# Variables saved are Year (1976-2008) and
# for each of the four lakes, Mendota-Monona-Waubesa-Kegonsa:
# Total Load, Load from Land, P0, P1,flushing coef., sedimentation coef., Export
# AND
# Mean Log TP, Mean Log Secchi, Mean Log DRP, Proportion Detectable DRP,
# AND presence/absence of D. pulicaria for Mendota and Monona only
# Save line was
# save(YR7608,
#      Me_in,M.1,P0.Me,P1.Me,f.1,s.1,Me.out,
#      LTP.Me,LSecchi.Me,LDRP.Me,DetectDRP.Me,Dpul.Me,
#      Mo_in,M.2,P0.Mo,P1.Mo,f.2,s.2,Mo.out,
#      LTP.Mo,LSecchi.Mo,LDRP.Mo,DetectDRP.Mo,Dpul.Mo,
#      Wa_in,M.3,P0.Wa,P1.Wa,f.3,s.3,Wa.out,
#      LTP.Wa,LSecchi.Wa,LDRP.Wa,DetectDRP.Wa,
#      Ke_in,M.4,P0.Ke,P1.Ke,f.4,s.4,Ke.out,
#      LTP.Ke,LSecchi.Ke,LDRP.Ke,DetectDRP.Ke,
#      file='PBudgets+MeanLogWQ-4Lakes.Rdata')
load(file='PBudgets+MeanLogWQ-4Lakes.Rdata')

LTPreg.Me = lm(LTP.Me ~ Me_in + Dpul.Me)
print(summary(LTPreg.Me))

LTPreg.Mo = lm(LTP.Mo ~  Mo.out + Dpul.Mo)
print(summary(LTPreg.Mo))

LTPreg.Wa = lm(LTP.Wa ~ Wa.out)
print(summary(LTPreg.Wa))

LTPreg.Ke = lm(LTP.Ke ~ Ke.out)
print(summary(LTPreg.Ke))

LSecreg.Me = lm(LSecchi.Me ~ Me_in + Dpul.Me)
print(summary(LSecreg.Me))

LSecreg.Mo = lm(LSecchi.Mo ~ Mo.out + Dpul.Mo)
print(summary(LSecreg.Mo))

LSecreg.Wa = lm(LSecchi.Wa ~ P1.Wa)
print(summary(LSecreg.Wa))

LSecreg.Ke = lm(LSecchi.Ke ~ P1.Ke)
print(summary(LSecreg.Ke))

# Save Water Quality Regressions
# Predictors for Log summer TP:
#   Mendota input and Dpul
#   Monona output and Dpul
#   Waubesa and Kegonsa: outputs
# Predictors for Log summer Secchi transparency:
#   Mendota input and Dpul
#   Monona output and Dpul
#   Waubesa and Kegonsa P1 (NOvember 1 P after the growing season)
save(LTPreg.Me,LTPreg.Mo,LTPreg.Wa,LTPreg.Ke,
     LSecreg.Me,LSecreg.Mo,LSecreg.Wa,LSecreg.Ke,
     file = 'WaterQualityRegressions.Rdata')

# Plot Secchi results
xrange = range(LSecreg.Me$fitted.values,LSecreg.Mo$fitted.values,
               LSecreg.Wa$fitted.values,LSecreg.Ke$fitted.values)
yrange = range(LSecchi.Me,LSecchi.Mo,LSecchi.Wa,LSecchi.Ke,na.rm=T)
windows()
par(cex.axis=1.5,cex.lab=1.5)
plot(LSecreg.Me$fitted.values,LSecchi.Me,type='p',pch=19,col='purple',
     xlim=xrange,ylim=yrange,xlab='Predicted log Secchi',ylab='Observed log Secchi',
     main='Summer means')
points(LSecreg.Mo$fitted.values,LSecchi.Mo,type='p',pch=19,col='blue')
points(LSecreg.Wa$fitted.values[5:33],LSecchi.Wa[5:33],type='p',pch=19,col='forestgreen')
points(LSecreg.Ke$fitted.values[5:33],LSecchi.Ke[5:33],type='p',pch=19,col='red')
abline(a=0,b=1,lty=2,lwd=2)
legend('topleft',legend = c('Mendota','Monona','Waubesa','Kegonsa'),
       bty='n',pch = c(19,19,19,19),
       col=c('purple','blue','forestgreen','red'))

