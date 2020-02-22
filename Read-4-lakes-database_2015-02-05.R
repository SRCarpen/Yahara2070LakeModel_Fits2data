# Script to read Yahara 4 lakes annualized data
# S.R. Carpenter, 2015-02-05

rm(list = ls())  # clear memory

# Load the complete database for P Mass Budgets and Lake TP
# Original source is PBudgets_LogMeanWQ_4Lakes_2013-01-26.r
# P mass budgets and summer TP extracted by Extract_Pdata_Only_2013-08-05.R

# Variables saved are Year (1976-2008) and
# for each of the four lakes, Mendota-Monona-Waubesa-Kegonsa:
# Total Load, Load from Land, P0, P1, Export, Mean summer TP, Daphnia for Me and Mo only

# Save line was
# save(YR7608,
#      L.Me,LandLoad.Me,P0.Me,P1.Me,XP.Me,LTP.Me,Dpul.Me,
#      L.Mo,LandLoad.Mo,P0.Mo,P1.Mo,XP.Mo,LTP.Mo,Dpul.Mo,
#      L.Wa,LandLoad.Wa,P0.Wa,P1.Wa,XP.Wa,LTP.Wa,
#      L.Ke,LandLoad.Ke,P0.Ke,P1.Ke,XP.Ke,LTP.Ke,
#      file='PMassBudgets+LogTP-4Lakes.Rdata')
load(file='PMassBudgets+LogTP-4Lakes.Rdata')
