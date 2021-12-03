#------------------------------------------------------------------------------
# The code is to reproduce the figure 5-6
# Note: Prior to running this code, the code for figure 3 should be run
#------------------------------------------------------------------------------

## Loading required R packages
library(EnvStats)
library(readxl)     
library(FME)        
library(mrgsolve)
library(dplyr)
library(minpack.lm)  
library(ggplot2)
library(tidyr)
library(ggprism)    
library(patchwork)  

# Input the PBPK model
source (file = "GenPBPK.R") # Loading the generic PBPK model code
source (file = "GetPC.R")   # Loading the tissue composition based tissue:plasma model
source (file = 'Pars.R')    # Loading the parameters

## Load Model
mod <- mcode_cache("pbpk", GenricPBPK) 

## Set up the seed and Monte Carlo iterations numbers
set.seed (123)
N = 1000

# ------------------------------------------------------------------------------
# Penicillin G (PG) simulation scenario
# Label dose (IU): 3000 units/lb bwt/day (https://www.etoolsage.com/converter/IU_Converter.asp)
# Label dose (mg/kg): 5 repeat dose 6.5 mg/kg; (lb = 0.454 kg) via IM injection
# Extra-label dose: 5 repeat dose 32.5 mg/kg (5*label dose) via IM injection
# Tolerance of PG: 0.05 ug/g for tissues in cattle (USFDA); 
# Action limit of PG: 0.025 ug/g for swine (FSIS, 2013);
#-------------------------------------------------------------------------------

## Calculation of idata
## PG
idata_C_PG <- MCsim (Pars_C_PG)
idata_S_PG <- MCsim (Pars_S_PG)


## Scenario A (label use): 6.5 mg/kg x 5 via IM injection
## Scenario B (extra lab use): 32.5 (5x label dose) mg/kg x 5 via IM injection

## Cattle
Sim_C_A_PG <-pred (pars = Pars_C_PG, drug = "PG", idata = idata_C_PG, 
                   tinterval = 24, Dose = 6.5, 
                   Dtimes = 5, route = 'im')

Sim_C_B_PG <-pred (pars = Pars_C_PG,drug = "PG", idata = idata_C_PG, 
                   tinterval = 24, Dose = 32.5, 
                   Dtimes = 5, route = 'im')

## Swine
Sim_S_A_PG <-pred (pars = Pars_S_PG, drug = "PG", idata = idata_S_PG, 
                   tinterval = 24, Dose = 6.5, 
                   Dtimes = 5, route = 'im')

Sim_S_B_PG <-pred (pars = Pars_S_PG,drug = "PG", idata = idata_S_PG, 
                   tinterval = 24, Dose = 32.5, 
                   Dtimes = 5, route = 'im')



#-------------------------------------------------------------------------------
## Ggplot
## Cattle

p_PG_C_A1_L<- MCplot(Sim_C_A_PG, target = "Liver", TOL = 0.05, tdoses = 5, Y.limit = 1e-5)  
p_PG_C_A1_K<- MCplot(Sim_C_A_PG, target = "Kidney",TOL = 0.05, tdoses = 5, Y.limit = 1e-5, color = 'bisque3')  
p_PG_C_A1_M<- MCplot(Sim_C_A_PG, target = "Muscle",TOL = 0.05, tdoses = 5, Y.limit = 1e-5, color = 'yellow')  

p1_PG <-p_PG_C_A1_L+p_PG_C_A1_K+p_PG_C_A1_M


## Swine
p_PG_S_A1_L<- MCplot(Sim_S_A_PG, target = "Liver" , TOL = 0.025, tdoses = 5, Y.limit = 1e-4)  
p_PG_S_A1_K<- MCplot(Sim_S_A_PG, target = "Kidney", TOL = 0.025, tdoses = 5, Y.limit = 1e-4, color = 'bisque3')  
p_PG_S_A1_M<- MCplot(Sim_S_A_PG, target = "Muscle", TOL = 0.025, tdoses = 5, Y.limit = 1e-4, color = 'yellow')  

p2_PG <-p_PG_S_A1_L + p_PG_S_A1_K + p_PG_S_A1_M


## Estimated the Withdraw intervals
PG_WDI_C_A<-Sim_C_A_PG %>% group_by (Time) %>% summarise (
    M99 = quantile(CM, probs  = 0.99),
    L99 = quantile(CL, probs  = 0.99),
    K99 = quantile(CK, probs  = 0.99)
)

PG_WDI_C_B<-Sim_C_B_PG %>% group_by (Time) %>% summarise (
    M99 = quantile(CM, probs  = 0.99),
    L99 = quantile(CL, probs  = 0.99),
    K99 = quantile(CK, probs  = 0.99)
)


PG_WDI_S_A<-Sim_S_A_PG %>% group_by (Time) %>% summarise (
    M99 = quantile(CM, probs  = 0.99),
    L99 = quantile(CL, probs  = 0.99),
    K99 = quantile(CK, probs  = 0.99)
)

PG_WDI_S_B<-Sim_S_B_PG %>% group_by (Time) %>% summarise (
    M99 = quantile(CM, probs  = 0.99),
    L99 = quantile(CL, probs  = 0.99),
    K99 = quantile(CK, probs  = 0.99)
)



## Tolerance of PG 50 ppb (0.05 ug/g) in Edible Tissue in cattle
WDIs_PG_C_L_A <- PG_WDI_C_A %>% filter(round(L99,2) <= 0.05 & Time >= 24*5)%>%mutate(WDIs = (Time/24)-5)%>%select(WDIs)%>%min()
WDIs_PG_C_K_A <- PG_WDI_C_A %>% filter(round(K99,2) <= 0.05 & Time >= 24*5)%>%mutate(WDIs = (Time/24)-5)%>%select(WDIs)%>%min()
WDIs_PG_C_M_A <- PG_WDI_C_A %>% filter(round(M99,2) <= 0.05 & Time >= 24*5)%>%mutate(WDIs = (Time/24)-5)%>%select(WDIs)%>%min()

WDIs_PG_C_L_B <- PG_WDI_C_B %>% filter(round(L99,2) <= 0.05 & Time >= 24*5)%>%mutate(WDIs = (Time/24)-5)%>%select(WDIs)%>%min()
WDIs_PG_C_K_B <- PG_WDI_C_B %>% filter(round(K99,2) <= 0.05 & Time >= 24*5)%>%mutate(WDIs = (Time/24)-5)%>%select(WDIs)%>%min()
WDIs_PG_C_M_B <- PG_WDI_C_B %>% filter(round(M99,2) <= 0.05 & Time >= 24*5)%>%mutate(WDIs = (Time/24)-5)%>%select(WDIs)%>%min()

## Tolerance of PG 25 ppb (0.025 ug/g) in Edible Tissue in swine
WDIs_PG_S_L_A <- PG_WDI_S_A %>% filter(round(L99,2) <= 0.025 & Time >= 24*5)%>%mutate(WDIs = (Time/24)-5)%>%select(WDIs)%>%min()
WDIs_PG_S_K_A <- PG_WDI_S_A %>% filter(round(K99,2) <= 0.025 & Time >= 24*5)%>%mutate(WDIs = (Time/24)-5)%>%select(WDIs)%>%min()
WDIs_PG_S_M_A <- PG_WDI_S_A %>% filter(round(M99,2) <= 0.025 & Time >= 24*5)%>%mutate(WDIs = (Time/24)-5)%>%select(WDIs)%>%min()

WDIs_PG_S_L_B <- PG_WDI_S_B %>% filter(round(L99,2) <= 0.025 & Time >= 24*5)%>%mutate(WDIs = (Time/24)-5)%>%select(WDIs)%>%min()
WDIs_PG_S_K_B <- PG_WDI_S_B %>% filter(round(K99,2) <= 0.025 & Time >= 24*5)%>%mutate(WDIs = (Time/24)-5)%>%select(WDIs)%>%min()
WDIs_PG_S_M_B <- PG_WDI_S_B %>% filter(round(M99,2) <= 0.025 & Time >= 24*5)%>%mutate(WDIs = (Time/24)-5)%>%select(WDIs)%>%min()

# ------------------------------------------------------------------------------
# Flunixin (FLU) simulation scenario 
# Label dose for cattle: 2.2 mg/kg of 3 repeated IV injections
# Label dose for swine: 2.2 mg/kg of single IM injections
# Extra-label use for Cattle: 2.2 mg/kg of 3 repeated IM
# Extra-label use for Swine: 2.2  mg/kg of 3  repeated  IM  injections   
# Tolerance of FLU for Cattle: 0.125 (ug/g) for liver; 0.025 (ug/g) for muscle
# Tolerance of FLU for Swine: 0.03 (ug/g) for liver; 0.025 (ug/g) for muscle
#-------------------------------------------------------------------------------
## Calculation of idata
## FLU

idata_C_FLU <- MCsim (Pars_C_FLU)
idata_S_FLU <- MCsim (Pars_S_FLU)


## Cattle
## Scenario A (Label dose) 2.2 mg/kg x 3 via IV injection
## Scenario B (Extra label) 2.2 mg/kg x 3 via IM injection

## Label use
Sim_C_A_FLU <-pred (pars      = Pars_C_FLU, 
                    drug      = "FLU", 
                    idata     = idata_C_FLU, 
                    tinterval = 24, 
                    Dose      = 2.2, 
                    Dtimes    = 3, 
                    route     = 'iv')

## Extra label use
Sim_C_B_FLU <-pred (pars      = Pars_C_FLU, 
                    drug      = "FLU", 
                    idata     = idata_C_FLU, 
                    tinterval = 24, 
                    Dose      = 2.2, 
                    Dtimes    = 3, 
                    route     = 'im')


## Swine
## Scenario A (Label dose): 2.2 mg/kg x 1 via IM injection
## Scenario B: 2.2 mg/kg x 3 via IM injection

## Label use
Sim_S_A_FLU <-pred (pars      = Pars_S_FLU, 
                    drug      = "FLU", 
                    idata     = idata_S_FLU, 
                    tinterval = 24, 
                    Dose      = 2.2, 
                    Dtimes    = 1, 
                    route     = 'im')

## Extra label use
Sim_S_B_FLU <-pred (pars      = Pars_S_FLU, 
                    drug      = "FLU", 
                    idata     = idata_S_FLU, 
                    tinterval = 24, 
                    Dose      = 2.2, 
                    Dtimes    = 3, 
                    route     = 'im')


#--------------------------------_------------------
## Ggplot
## Cattle
p_FLU_C_A1_L<- MCplot(Sim_C_A_FLU, target = "Liver", tdoses = 3, TOL  = 0.125,Y.limit=1e-6) 
p_FLU_C_A1_K<- MCplot(Sim_C_A_FLU, target = "Kidney", tdoses = 3, TOL = 0.125, color = 'bisque3',Y.limit=1e-6) 
p_FLU_C_A1_M<- MCplot(Sim_C_A_FLU, target = "Muscle", tdoses = 3, TOL = 0.025, color = 'yellow',Y.limit=1e-6) 

p1_FLU<-(p_FLU_C_A1_L + p_FLU_C_A1_K + p_FLU_C_A1_M)

## Swine
p_FLU_S_A1_L<- MCplot(Sim_S_A_FLU, target = "Liver", tdoses = 1, TOL =0.03, Y.limit=1e-4) 
p_FLU_S_A1_K<- MCplot(Sim_S_A_FLU, target = "Kidney", tdoses = 1, TOL=0.03,color = 'bisque3',Y.limit=1e-4) 
p_FLU_S_A1_M<- MCplot(Sim_S_A_FLU, target = "Muscle", tdoses = 1, TOL=0.025,color = 'yellow',Y.limit=1e-4) 

p2_FLU<-(p_FLU_S_A1_L + p_FLU_S_A1_K + p_FLU_S_A1_M)


## Estimated the Withdraw intervals
FLU_WDI_C_A<-Sim_C_A_FLU %>% group_by (Time) %>% summarise (
    M99 = quantile(CM, probs  = 0.99),
    L99 = quantile(CL, probs  = 0.99),
    K99 = quantile(CK, probs  = 0.99)
)

FLU_WDI_C_B<-Sim_C_B_FLU %>% group_by (Time) %>% summarise (
    M99 = quantile(CM, probs  = 0.99),
    L99 = quantile(CL, probs  = 0.99),
    K99 = quantile(CK, probs  = 0.99)
)

## Swine
FLU_WDI_S_A<-Sim_S_A_FLU %>% group_by (Time) %>% summarise (
    M99 = quantile(CM, probs  = 0.99),
    L99 = quantile(CL, probs  = 0.99),
    K99 = quantile(CK, probs  = 0.99)
)

FLU_WDI_S_B<-Sim_S_B_FLU %>% group_by (Time) %>% summarise (
    M99 = quantile(CM, probs  = 0.99),
    L99 = quantile(CL, probs  = 0.99),
    K99 = quantile(CK, probs  = 0.99)
)

## Tolerance (TOL) of flunixin for liver is 0.125 µg/g, and for muscle is 0.025 µg/g in cattle.
## TOL for liver (0.125 µg/g) was used for kidney and plasma
WDIs_FLU_C_L_A <- FLU_WDI_C_A %>% filter(round(L99,3) <= 0.125 & Time >=24*3)%>%mutate(WDIs = (Time/24)-3)%>%select(WDIs)%>%min()
WDIs_FLU_C_K_A <- FLU_WDI_C_A %>% filter(round(K99,3) <= 0.125 & Time >=24*3)%>%mutate(WDIs = (Time/24)-3)%>%select(WDIs)%>%min()
WDIs_FLU_C_M_A <- FLU_WDI_C_A %>% filter(round(M99,3) <= 0.025 & Time >=24*3)%>%mutate(WDIs = (Time/24)-3)%>%select(WDIs)%>%min()

WDIs_FLU_C_L_B <- FLU_WDI_C_B %>% filter(round(L99,3) <= 0.125 & Time >=24*3)%>%mutate(WDIs = (Time/24)-3)%>%select(WDIs)%>%min()
WDIs_FLU_C_K_B <- FLU_WDI_C_B %>% filter(round(K99,3) <= 0.125 & Time >=24*3)%>%mutate(WDIs = (Time/24)-3)%>%select(WDIs)%>%min()
WDIs_FLU_C_M_B <- FLU_WDI_C_B %>% filter(round(M99,3) <= 0.025 & Time >=24*3)%>%mutate(WDIs = (Time/24)-3)%>%select(WDIs)%>%min()


## The TOL of flunixin in swine is 0.030 µg/g for liver, and 0.025 µg/g for muscle
## TOL for liver (0.030 µg/g in swine) was used for kidney and plasma

WDIs_FLU_S_L_A <- FLU_WDI_S_A %>% filter(round(L99,3) <= 0.030& Time >=24*1)%>%mutate(WDIs = (Time/24)-1)%>%select(WDIs)%>%min()
WDIs_FLU_S_K_A <- FLU_WDI_S_A %>% filter(round(K99,3) <= 0.030& Time >=24*1)%>%mutate(WDIs = (Time/24)-1)%>%select(WDIs)%>%min()
WDIs_FLU_S_M_A <- FLU_WDI_S_A %>% filter(round(M99,3) <= 0.025& Time >=24*1)%>%mutate(WDIs = (Time/24)-1)%>%select(WDIs)%>%min()

WDIs_FLU_S_L_B <- FLU_WDI_S_B %>% filter(round(L99,3) <= 0.030& Time >=24*3)%>%mutate(WDIs = (Time/24)-3)%>%select(WDIs)%>%min()
WDIs_FLU_S_K_B <- FLU_WDI_S_B %>% filter(round(K99,3) <= 0.030& Time >=24*3)%>%mutate(WDIs = (Time/24)-3)%>%select(WDIs)%>%min()
WDIs_FLU_S_M_B <- FLU_WDI_S_B %>% filter(round(M99,3) <= 0.025& Time >=24*3)%>%mutate(WDIs = (Time/24)-3)%>%select(WDIs)%>%min()



# ------------------------------------------------------------------------------
# Florfenicol (FLO) and simulation scenario 
# Label dose for cattle: 20 mg/kg of 2 repeated IM injections at 48 hours interval
# Alternative label dose for cattle: 40 mg/kg of single SC injection at 24 hours interval
# Label dose for swine: 100 ppm in water for 5 days via oral administration
# Extra-label use for Cattle: 20 mg/kg of 3 repeated IM injections at 48-hours interval
# Alternative extra-label dose for cattle: 40 mg/kg of 3 repeated SC injections at 96-hours interval
# Extra-label use for Swine: 40  mg/kg of 2  repeated  SC  injections  
# Tolerance in US (FLOA) in cattle: 3.7 (ug/g) in liver; 0.3 (ug/g) in muscle (FDA, 2017)
# Tolerance in China and EU (FLO+FLOA): 3 (ug/g) in liver; 0.2 (ug/g) in muscle; 0.3 (ug/g) in kidney;
# MRLs for FLO+FLOA in swine in US: 2.5 (ug/g) in liver; 0.2 (ug/g) in muscle;
# MRLs for FLO+FLOA in swine in China and EU: 2 (ug/g) in liver; 0.3 (ug/g) in muscle; 0.5 (ug/g) in Kidney
#-------------------------------------------------------------------------------
## Calculation of idata
## FLO
idata_C_FLO <- MCsim (Pars_C_FLO)
idata_S_FLO <- MCsim (Pars_S_FLO)



## Cattle
## Scenario A (label use): 20 mg/kg x 2 via IM injection (48 hrs interval)
## Scenario B (label use): 40 mg/kg x 1 via SC injection (24 hrs interval)
## Scenario C (Extra-label use): 20 mg/kg x 3 via IM injection (48 hrs interval)
## Scenario D (Extra-label use): 40 mg/kg x 3 via SC injection (96 hrs interval)

Sim_C_A_FLO <-pred (pars = Pars_C_FLO, 
                    drug = "FLO", idata = idata_C_FLO, 
                    tinterval = 48, Dose = 20, Dtimes = 2, route = 'im')

Sim_C_B_FLO <-pred (pars = Pars_C_FLO, 
                    drug = "FLO", idata = idata_C_FLO, 
                    tinterval = 24, Dose = 40, Dtimes = 1, route = 'sc')

Sim_C_C_FLO <-pred (pars = Pars_C_FLO, 
                    drug = "FLO", idata = idata_C_FLO, 
                    tinterval = 48, Dose = 20, Dtimes = 3, route = 'im')

Sim_C_D_FLO <-pred (pars = Pars_C_FLO, 
                    drug = "FLO", idata = idata_C_FLO, 
                    tinterval = 96, Dose = 40, Dtimes = 3, route = 'sc')
## Swine
## Scenario A (label use): 30 mg/kg x 5 via PO (24 hrs interval)
## 100 mg/l (in water) * 10 l/day (water consumption from Thacker (2001)) = 1000 mg/day
## 1000 mg/day / 70 kg (swine body weight) = 14 mg/kg/day

## Scenario B (extra-label use): 40 mg/kg x 2 via SC  (96 hr interval)

Sim_S_A_FLO <-pred (pars = Pars_S_FLO, 
                    drug = "FLO", idata = idata_S_FLO, 
                    tinterval = 24, Dose = 14, Dtimes = 5, route = 'po')

Sim_S_B_FLO <-pred (pars = Pars_S_FLO, 
                    drug = "FLO", idata = idata_S_FLO, 
                    tinterval = 96, Dose = 40, Dtimes = 2, route = 'sc')



#-------------------------------------------------------------------------------
## Ggplot
## Cattle
## Cattle
p_FLO_C_A1_L<- MCplot(Sim_C_A_FLO, target = "Liver1", tinterval = 48, tdoses = 2,  TOL=3.7,Y.limit = 1e-6) 
p_FLO_C_A1_K<- MCplot(Sim_C_A_FLO, target = "Kidney1", tinterval = 48, tdoses = 2, TOL=3.7, color = 'bisque3',Y.limit=1e-6) 
p_FLO_C_A1_M<- MCplot(Sim_C_A_FLO, target = "Muscle1", tinterval = 48, tdoses = 2, TOL=0.3, color = 'yellow',Y.limit=1e-6) 

p1_FLO<-(p_FLO_C_A1_L + p_FLO_C_A1_K + p_FLO_C_A1_M)

## Swine
p_FLO_S_A1_L<- MCplot(Sim_S_A_FLO, target = "Liver1",  tdoses = 5, TOL=2.5,Y.limit =1e-6) 
p_FLO_S_A1_K<- MCplot(Sim_S_A_FLO, target = "Kidney1", tdoses = 5, TOL=2.5,color = 'bisque3',Y.limit=1e-6) 
p_FLO_S_A1_M<- MCplot(Sim_S_A_FLO, target = "Muscle1", tdoses = 5, TOL=0.2,color = 'yellow',Y.limit=1e-6) 

p2_FLO<-(p_FLO_S_A1_L + p_FLO_S_A1_K + p_FLO_S_A1_M)

#p1_FLO /p2_FLO 


## Estimated the Withdraw intervals
FLO_WDI_C_A<-Sim_C_A_FLO %>% group_by (Time) %>% summarise (
    M99 = quantile(CM_Met, probs  = 0.99),
    L99 = quantile(CL_Met, probs  = 0.99),
    K99 = quantile(CK_Met, probs  = 0.99)
)

FLO_WDI_C_B<-Sim_C_B_FLO %>% group_by (Time) %>% summarise (
    M99 = quantile(CM_Met, probs  = 0.99),
    L99 = quantile(CL_Met, probs  = 0.99),
    K99 = quantile(CK_Met, probs  = 0.99)
)

FLO_WDI_C_C<-Sim_C_C_FLO %>% group_by (Time) %>% summarise (
    M99 = quantile(CM_Met, probs  = 0.99),
    L99 = quantile(CL_Met, probs  = 0.99),
    K99 = quantile(CK_Met, probs  = 0.99)
)

FLO_WDI_C_D<-Sim_C_D_FLO %>% group_by (Time) %>% summarise (
    M99 = quantile(CM_Met, probs  = 0.99),
    L99 = quantile(CL_Met, probs  = 0.99),
    K99 = quantile(CK_Met, probs  = 0.99)
)

## Swine
FLO_WDI_S_A<-Sim_S_A_FLO %>% group_by (Time) %>% summarise (
    M99 = quantile(CM_Met, probs  = 0.99),
    L99 = quantile(CL_Met, probs  = 0.99),
    K99 = quantile(CK_Met, probs  = 0.99)
)

FLO_WDI_S_B<-Sim_S_B_FLO %>% group_by (Time) %>% summarise (
    M99 = quantile(CM_Met, probs  = 0.99),
    L99 = quantile(CL_Met, probs  = 0.99),
    K99 = quantile(CK_Met, probs  = 0.99)
)

## Tolerance (TOL) of Florfenicol for liver is 3.7 µg/g  and for muscle is 0.3 µg/g  in cattle.
## TOL for liver was used for kidney 
WDIs_FLO_C_L_A <- FLO_WDI_C_A %>% filter(round(L99,3) <= 3.7 & Time >=24*4)%>%mutate(WDIs = (Time/24)-4)%>%select(WDIs)%>%min()
WDIs_FLO_C_K_A <- FLO_WDI_C_A %>% filter(round(K99,3) <= 3.7 & Time >=24*4)%>%mutate(WDIs = (Time/24)-4)%>%select(WDIs)%>%min()
WDIs_FLO_C_M_A <- FLO_WDI_C_A %>% filter(round(M99,3) <= 0.3 & Time >=24*4)%>%mutate(WDIs = (Time/24)-4)%>%select(WDIs)%>%min()

WDIs_FLO_C_L_B <- FLO_WDI_C_B %>% filter(round(L99,3) <= 3.7 & Time >=24*1)%>%mutate(WDIs = (Time/24)-1)%>%select(WDIs)%>%min()
WDIs_FLO_C_K_B <- FLO_WDI_C_B %>% filter(round(K99,3) <= 3.7 & Time >=24*1)%>%mutate(WDIs = (Time/24)-1)%>%select(WDIs)%>%min()
WDIs_FLO_C_M_B <- FLO_WDI_C_B %>% filter(round(M99,3) <= 0.3 & Time >=24*1)%>%mutate(WDIs = (Time/24)-1)%>%select(WDIs)%>%min()

WDIs_FLO_C_L_C <- FLO_WDI_C_C %>% filter(round(L99,3) <= 3.7 & Time >=24*6)%>%mutate(WDIs = (Time/24)-6)%>%select(WDIs)%>%min()
WDIs_FLO_C_K_C <- FLO_WDI_C_C %>% filter(round(K99,3) <= 3.7 & Time >=24*6)%>%mutate(WDIs = (Time/24)-6)%>%select(WDIs)%>%min()
WDIs_FLO_C_M_C <- FLO_WDI_C_C %>% filter(round(M99,3) <= 0.3 & Time >=24*6)%>%mutate(WDIs = (Time/24)-6)%>%select(WDIs)%>%min()

WDIs_FLO_C_L_D <- FLO_WDI_C_D %>% filter(round(L99,3) <= 3.7 & Time >=24*12)%>%mutate(WDIs = (Time/24)-12)%>%select(WDIs)%>%min()
WDIs_FLO_C_K_D <- FLO_WDI_C_D %>% filter(round(K99,3) <= 3.7 & Time >=24*12)%>%mutate(WDIs = (Time/24)-12)%>%select(WDIs)%>%min()
WDIs_FLO_C_M_D <- FLO_WDI_C_D %>% filter(round(M99,3) <= 0.3 & Time >=24*12)%>%mutate(WDIs = (Time/24)-12)%>%select(WDIs)%>%min()

## The TOL of Florfenicol in swine is 2.5 µg/g for liver, and 0.2 µg/g for muscle
## TOL for liver was used for kidney and plasma
WDIs_FLO_S_L_A <- FLO_WDI_S_A %>% filter(round(L99,3) <= 2.5 & Time >=24*5)%>%mutate(WDIs = (Time/24)-5)%>%select(WDIs)%>%min()
WDIs_FLO_S_K_A <- FLO_WDI_S_A %>% filter(round(K99,3) <= 2.5 & Time >=24*5)%>%mutate(WDIs = (Time/24)-5)%>%select(WDIs)%>%min()
WDIs_FLO_S_M_A <- FLO_WDI_S_A %>% filter(round(M99,3) <= 0.2 & Time >=24*5)%>%mutate(WDIs = (Time/24)-5)%>%select(WDIs)%>%min()

WDIs_FLO_S_L_B <- FLO_WDI_S_B %>% filter(round(L99,3) <= 2.5 & Time >=24*8)%>%mutate(WDIs = (Time/24)-8)%>%select(WDIs)%>%min()
WDIs_FLO_S_K_B <- FLO_WDI_S_B %>% filter(round(K99,3) <= 2.5 & Time >=24*8)%>%mutate(WDIs = (Time/24)-8)%>%select(WDIs)%>%min()
WDIs_FLO_S_M_B <- FLO_WDI_S_B %>% filter(round(M99,3) <= 0.2 & Time >=24*8)%>%mutate(WDIs = (Time/24)-8)%>%select(WDIs)%>%min()


##///////////////////////////////////////////////////////////////////////////////
## Plot


ggsave("Fig 5.tiff",scale = 1.5,
       plot = p1_FLU/p1_FLO/p1_PG,
       width = 25, height = 15, units = "cm", dpi=320)


ggsave("Fig 6.tiff",scale = 1.5,
       plot = p2_FLU/p2_FLO/p2_PG,
       width = 25, height = 15, units = "cm", dpi=320)







