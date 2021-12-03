#------------------------------------------------------------------------------
# The code is to reproduce the Table S4-S7
#------------------------------------------------------------------------------

## Loading required R package
library(FME)
library(mrgsolve)
library(dplyr)
library(minpack.lm)  
library(ggplot2)
library(tidyr)
library(patchwork)

## Input the PBPK model code and requring code
source (file = "GenPBPK.R")  # Loading the generic PBPK model code
source (file = "CirPlot.R")  # Loading the circleplot code
source (file = 'Pars.R')     # Loading the parameters

## Load Model
mod <- mcode_cache("pbpk", GenricPBPK) #refer to mcode function in mrgsolve user guide 3.1.2 Inline

## MC simulation

### Define the prediction function
pred <- function (pars, drug, tinterval, Dose, tdose, route) {
  
  ## Get out of log domain
  parsinput <- exp(pars)
  BW <- parsinput["BW"]
  
  Fracsc = ifelse(is.na(parsinput["Fracsc"]), 1, parsinput["Fracsc"]) 
  Fracim = ifelse(is.na(parsinput["Fracim"]), 1, parsinput["Fracim"]) 
  
  ## Exposure scenarios
  BW          = BW
  tinterval   = tinterval
  TDOSE       = tdose
  Frac        = switch(route, "im"= Fracim, "sc" = Fracsc)
  MW          = switch(drug, "PG" = 334.4, "FLO" =358.21, "FLU" = 296.24)
  MW1         = switch(drug, "PG" = 334.4, "FLO" = 247.28, "FLU" = 312.24)
  Dose        = Dose
  DOSEfast    = Dose*BW*Frac/MW              
  DOSEslow    = Dose*BW*(1-Frac)/MW  
  DOSE        = Dose*BW/MW  
  
  ## 
  if (route == "im") {
    ev_1 <- ev (ID   = 1, amt  = DOSEfast, ii = tinterval, 
                addl = TDOSE - 1, cmt  = "Amtsiteim", replicate = FALSE)
    ev_2 <- ev (ID   = 1, amt  = DOSEslow, ii = tinterval, 
                addl = TDOSE - 1, cmt  = "ADOSEim", replicate = FALSE)
    ev_3 <- ev (ID   = 1, amt  = DOSEfast+DOSEslow, ii = tinterval, 
                addl = TDOSE - 1, cmt  = "ADOSE", replicate = FALSE)
    ex <- ev_1 + ev_2 + ev_3
    
  }
  
  if (route == "sc") {
    ev_1 <- ev (ID   = 1, amt  = DOSEfast, ii = tinterval, 
                addl = TDOSE - 1, cmt  = "Amtsitesc", replicate = FALSE)
    ev_2 <- ev (ID   = 1, amt  = DOSEslow, ii = tinterval, 
                addl = TDOSE - 1, cmt  = "ADOSEsc", replicate = FALSE)
    
    ev_3 <- ev (ID   = 1, amt  = DOSEfast+DOSEslow, ii = tinterval, 
                addl = TDOSE - 1, cmt  = "ADOSE", replicate = FALSE)
    
    ex <- ev_1 + ev_2 + ev_3
  }
  
  if (route == "iv") {
    ev_1 <- ev (ID   = 1, amt  = DOSE, ii = tinterval, tinf = 0.01,
                addl = TDOSE - 1, cmt  = "APlas_free", replicate = FALSE)
    ev_2 <- ev (ID   = 1, amt  = DOSE, ii = tinterval, tinf = 0.01,
                addl = TDOSE - 1, cmt  = "ADOSE",  replicate = FALSE)
    ex <- ev_1+ev_2
    
  }
  
  if (route == "po") {
    ev_1 <- ev (ID   = 1, amt  = DOSE, ii = tinterval, 
                addl = TDOSE - 1, cmt  = "AST", replicate = FALSE)
    ev_2 <- ev (ID   = 1, amt  = DOSE, ii = tinterval, 
                addl = TDOSE - 1, cmt  = "ADOSE", replicate = FALSE)
    ex <- ev_1+ev_2
  }
  
  
  tsamp  = tgrid(0, tinterval*(TDOSE - 1) + 24*30, 0.1) 
  
  ## Simulation
  out <- mod %>% param (parsinput) %>% 
    update(atol = 1E-10, rtol = 1E-6, maxsteps = 50000) %>%
    mrgsim_d (data = ex, tgrid = tsamp)
  
  
  outdf = cbind.data.frame( Time  = out$time,
                            CP    = out$Plasma*MW,
                            CL    = out$Liver*MW,
                            CK    = out$Kidney*MW,
                            CM    = out$Muscle*MW,
                            CF    = out$Fat*MW,
                            AUCCP = out$AUC_CP*MW,
                            AUCCM = out$AUC_CM*MW,
                            AUCCL = out$AUC_CL*MW,
                            AUCCK = out$AUC_CK*MW,
                            AUCCP1 = out$AUC_CP1*MW1,
                            AUCCM1 = out$AUC_CM1*MW1,
                            AUCCL1 = out$AUC_CL1*MW1,
                            AUCCK1 = out$AUC_CK1*MW1)
  return (outdf)
  
}



## Define the sensitivity function
NSC_func <- function (pars, drug, tinterval, Dose, tdose, route, AUC="AUCCM") {
  n <- length(pars)
  NSC_AUC     <- matrix(NA, nrow = length(pars) , ncol = 1)
  
  for (i in 1:n) {
    if(drug == "FLU" & i <= 2) {
      pars.new     <- pars %>% replace(i, log(exp((pars[i]))*0.99))
      delta.pars1  <- (1-exp(pars.new))/exp(pars[i])
      new          <- pred(pars.new, drug, tinterval, Dose, tdose, route)
      ori          <- pred(pars, drug, tinterval, Dose, tdose, route)
      delta.pars   <- exp(pars[i])/(exp(pars[i])*-0.99)
      
    } else { 
    pars.new     <- pars %>% replace(i, log(exp((pars[i]))*1.01))
    new          <- pred(pars.new, drug, tinterval, Dose, tdose, route)
    ori          <- pred(pars, drug, tinterval, Dose, tdose, route)
    delta.pars   <- exp(pars[i])/(exp(pars[i])*0.01)
    }
    
    ## Estimated the AUC
    AUC.new   =  new[new$Time==max(new$Time),]-new[new$Time==tdose*tinterval,]
    AUC.ori   =  ori[ori$Time==max(ori$Time),]-ori[ori$Time==tdose*tinterval,]
    
    delta.AUC   =  (AUC.new - AUC.ori)%>%select(AUC)
    AUC.ori_sel =  AUC.ori%>%select(AUC)
    
    NSC_AUC [i, 1]   <- as.numeric((delta.AUC/AUC.ori_sel) * delta.pars)
    
  }
  
  return (list(NSC_AUC = NSC_AUC))
}

## Define the NSC 
# ------------------------------------------------------------------------------
# Penicillin G (PG) simulation scenario
# Label dose (IU): 3000 units/lb bwt/day (https://www.etoolsage.com/converter/IU_Converter.asp)
# Label dose (mg/kg): 5 repeat dose 6.5 mg/kg; (lb = 0.454 kg) via IM injection for swine and cattle
#-------------------------------------------------------------------------------

## Sensitivity function for Cattle
## Scenario: 6.5 mg/kg x 5 via IM injection
NSC_C_PG_L <- data.frame(NSC_func (log(Pars_C_PG), tinterval = 24, 
                                   tdose = 5, Dose = 6.5, AUC="AUCCL",
                                   route = "im", drug = 'PG'))


NSC_C_PG_K <- data.frame(NSC_func (log(Pars_C_PG), tinterval = 24, 
                                   tdose = 5, Dose = 6.5, AUC="AUCCK",
                                   route = "im", drug = 'PG'))


NSC_C_PG_M <- data.frame(NSC_func (log(Pars_C_PG), tinterval = 24, 
                                 tdose = 5, Dose = 6.5, AUC="AUCCM",
                                 route = "im", drug = 'PG'))



## Sensitivity function for Swine
## Scenario: 6.5 mg/kg x 5 via IM injection
NSC_S_PG_L <- data.frame(NSC_func (log(Pars_S_PG), tinterval = 24, 
                                 tdose = 5, Dose = 6.5, AUC="AUCCL",
                                 route = "im", drug = 'PG'))

NSC_S_PG_K <- data.frame(NSC_func (log(Pars_S_PG), tinterval = 24, 
                                 tdose = 5, Dose = 6.5, AUC="AUCCK",
                                 route = "im", drug = 'PG'))

NSC_S_PG_M <- data.frame(NSC_func (log(Pars_S_PG), tinterval = 24, 
                                 tdose = 5, Dose = 6.5, AUC="AUCCM",
                                 route = "im", drug = 'PG'))


## Name the row 
rownames(NSC_C_PG_L) = rownames(NSC_C_PG_K) = rownames(NSC_C_PG_M) = names(Pars_C_PG)
rownames(NSC_S_PG_L) = rownames(NSC_S_PG_K) = rownames(NSC_S_PG_M) = names(Pars_S_PG)

NSC_C_PG <- cbind.data.frame(NSC_C_PG_L, NSC_C_PG_K, NSC_C_PG_M)

NSC_S_PG <- cbind.data.frame(NSC_S_PG_L, NSC_S_PG_K, NSC_S_PG_M)


colnames(NSC_C_PG) = colnames(NSC_S_PG)=c("NSC_AUCL","NSC_AUCK","NSC_AUCM")

write.csv(NSC_C_PG, file = 'NSC_C_PG.csv')
write.csv(NSC_S_PG, file = 'NSC_S_PG.csv')


# ------------------------------------------------------------------------------
# Flunixin (FLU) simulation scenario 
# Label dose for cattle: 2.2 mg/kg of 3 repeated IV injections
# Label dose for swine: 2.2 mg/kg of single IM injections
#-------------------------------------------------------------------------------

## Sensitivity function for Cattle
## Scenario: 2.2 mg/kg x 3 via IV injection
NSC_Pars_C_FLU <- Pars_C_FLU [names(Pars_C_FLU) %in% c("Fracsc","Fracim")==FALSE]

NSC_C_FLU_L <- data.frame(NSC_func (log(NSC_Pars_C_FLU), tinterval = 24, 
                                  tdose = 3, Dose = 2.2, AUC="AUCCL",
                                  route = "iv", drug = 'FLU'))

NSC_C_FLU_K <- data.frame(NSC_func (log(NSC_Pars_C_FLU), tinterval = 24, 
                                  tdose = 3, Dose = 2.2, AUC="AUCCK",
                                  route = "iv", drug = 'FLU'))

NSC_C_FLU_M <- data.frame(NSC_func (log(NSC_Pars_C_FLU), tinterval = 24, 
                                  tdose = 3, Dose = 2.2, AUC="AUCCM",
                                  route = "iv", drug = 'FLU'))

NSC_C_FLU_L1 <- data.frame(NSC_func (log(NSC_Pars_C_FLU), tinterval = 24, 
                                    tdose = 3, Dose = 2.2, AUC="AUCCL1",
                                    route = "iv", drug = 'FLU'))

NSC_C_FLU_K1 <- data.frame(NSC_func (log(NSC_Pars_C_FLU), tinterval = 24, 
                                    tdose = 3, Dose = 2.2, AUC="AUCCK1",
                                    route = "iv", drug = 'FLU'))

NSC_C_FLU_M1 <- data.frame(NSC_func (log(NSC_Pars_C_FLU), tinterval = 24, 
                                    tdose = 3, Dose = 2.2, AUC="AUCCM1",
                                    route = "iv", drug = 'FLU'))
## Sensitivity function for Swine
## Scenario: 2.2 mg/kg x 1 via IM injection
NSC_Pars_S_FLU <- Pars_S_FLU [names(Pars_S_FLU) %in% c("Fracsc","Fracim")==FALSE]

NSC_S_FLU_L <- data.frame(NSC_func (log(NSC_Pars_S_FLU), tinterval = 24, 
                                  tdose = 1, Dose = 2.2, AUC="AUCCL",
                                  route = "im", drug = 'FLU'))

NSC_S_FLU_K <- data.frame(NSC_func (log(NSC_Pars_S_FLU), tinterval = 24, 
                                  tdose = 1, Dose = 2.2, AUC="AUCCK",
                                  route = "im", drug = 'FLU'))

NSC_S_FLU_M <- data.frame(NSC_func (log(NSC_Pars_S_FLU), tinterval = 24, 
                                  tdose = 1, Dose = 2.2, AUC="AUCCM",
                                  route = "im", drug = 'FLU'))

NSC_S_FLU_L1 <- data.frame(NSC_func (log(NSC_Pars_S_FLU), tinterval = 24, 
                                    tdose = 1, Dose = 2.2, AUC="AUCCL1",
                                    route = "im", drug = 'FLU'))

NSC_S_FLU_K1 <- data.frame(NSC_func (log(NSC_Pars_S_FLU), tinterval = 24, 
                                    tdose = 1, Dose = 2.2, AUC="AUCCK1",
                                    route = "im", drug = 'FLU'))

NSC_S_FLU_M1 <- data.frame(NSC_func (log(NSC_Pars_S_FLU), tinterval = 24, 
                                    tdose = 1, Dose = 2.2, AUC="AUCCM1",
                                    route = "im", drug = 'FLU'))


## Name the row 
rownames(NSC_C_FLU_L)  = rownames(NSC_C_FLU_K)  = rownames(NSC_C_FLU_M)  = names(NSC_Pars_C_FLU)
rownames(NSC_S_FLU_L)  = rownames(NSC_S_FLU_K)  = rownames(NSC_S_FLU_M)  = names(NSC_Pars_S_FLU)
rownames(NSC_C_FLU_L1) = rownames(NSC_C_FLU_K1) = rownames(NSC_C_FLU_M1) = names(NSC_Pars_C_FLU)
rownames(NSC_S_FLU_L1) = rownames(NSC_S_FLU_K1) = rownames(NSC_S_FLU_M1) = names(NSC_Pars_S_FLU)

NSC_C_FLU <- cbind.data.frame(NSC_C_FLU_L, NSC_C_FLU_K, NSC_C_FLU_M,
                             NSC_C_FLU_L1, NSC_C_FLU_K1, NSC_C_FLU_M1)

NSC_S_FLU <- cbind.data.frame(NSC_S_FLU_L, NSC_S_FLU_K, NSC_S_FLU_M,
                             NSC_S_FLU_L1,NSC_S_FLU_K1, NSC_S_FLU_M1)


colnames(NSC_C_FLU) = colnames(NSC_S_FLU)=c("NSC_AUCL","NSC_AUCK","NSC_AUCM",
                                          "NSC_AUCL1","NSC_AUCK1","NSC_AUCM1")


write.csv(NSC_C_FLU, file = 'NSC_C_FLU.csv')
write.csv(NSC_S_FLU, file = 'NSC_S_FLU.csv')


# ------------------------------------------------------------------------------
# Florfenicol (FLO) and simulation scenario 
# Label dose for cattle: 20 mg/kg of 2 repeated IM injections at 48 hours interval
# Label dose for swine: 100 ppm in water for 5 days via oral administration
#-------------------------------------------------------------------------------

## Sensitivity function for Cattle
## Scenario: 20 mg/kg x 2 via IM injection (48 hrs interval)
NSC_C_FLO_L <- data.frame(NSC_func (log(Pars_C_FLO), tinterval = 48, 
                                  tdose = 2, Dose = 20, AUC="AUCCL",
                                  route = "im", drug = 'FLO'))

NSC_C_FLO_K <- data.frame(NSC_func (log(Pars_C_FLO), tinterval = 48, 
                                  tdose = 2, Dose = 20, AUC="AUCCK",
                                  route = "im", drug = 'FLO'))

NSC_C_FLO_M <- data.frame(NSC_func (log(Pars_C_FLO), tinterval = 48, 
                                  tdose = 2, Dose = 20, AUC="AUCCM",
                                  route = "im", drug = 'FLO'))

NSC_C_FLO_L1 <- data.frame(NSC_func (log(Pars_C_FLO), tinterval = 48, 
                                    tdose = 2, Dose = 20, AUC="AUCCL1",
                                    route = "im", drug = 'FLO'))

NSC_C_FLO_K1 <- data.frame(NSC_func (log(Pars_C_FLO), tinterval = 48, 
                                    tdose = 2, Dose = 20, AUC="AUCCK1",
                                    route = "im", drug = 'FLO'))

NSC_C_FLO_M1 <- data.frame(NSC_func (log(Pars_C_FLO), tinterval = 48, 
                                    tdose = 2, Dose = 20, AUC="AUCCM1",
                                    route = "im", drug = 'FLO'))


## Sensitivity function for Swine
## Swine
## Scenario A (label use): 30 mg/kg x 5 via PO (24 hrs interval)
## 100 mg/l (in water) * 10 l/day (water consumption from Thacker (2001)) = 1000 mg/day
## 1000 mg/day / 70 kg (swine body weight) = 14 mg/kg/day
NSC_S_FLO_L <- data.frame(NSC_func (log(Pars_S_FLO), tinterval = 24, 
                                  tdose = 5, Dose = 14, AUC="AUCCL",
                                  route = "po", drug = 'FLO'))

NSC_S_FLO_K <- data.frame(NSC_func (log(Pars_S_FLO), tinterval = 24, 
                                  tdose = 5, Dose = 14, AUC="AUCCK",
                                  route = "po", drug = 'FLO'))

NSC_S_FLO_M <- data.frame(NSC_func (log(Pars_S_FLO), tinterval = 24, 
                                  tdose = 5, Dose = 14, AUC="AUCCM",
                                  route = "po", drug = 'FLO'))

NSC_S_FLO_L1 <- data.frame(NSC_func (log(Pars_S_FLO), tinterval = 24, 
                                    tdose = 5, Dose = 14, AUC="AUCCL1",
                                    route = "po", drug = 'FLO'))

NSC_S_FLO_K1 <- data.frame(NSC_func (log(Pars_S_FLO), tinterval = 24, 
                                    tdose = 5, Dose = 14, AUC="AUCCK1",
                                    route = "po", drug = 'FLO'))

NSC_S_FLO_M1 <- data.frame(NSC_func (log(Pars_S_FLO), tinterval = 24, 
                                    tdose = 5, Dose = 14, AUC="AUCCM1",
                                    route = "po", drug = 'FLO'))

## Name the row 
rownames(NSC_C_FLO_L) = rownames(NSC_C_FLO_K) = rownames(NSC_C_FLO_M) = names(Pars_C_FLO)
rownames(NSC_S_FLO_L) = rownames(NSC_S_FLO_K) = rownames(NSC_S_FLO_M) = names(Pars_S_FLO)
rownames(NSC_C_FLO_L1) = rownames(NSC_C_FLO_K1) = rownames(NSC_C_FLO_M1) = names(Pars_C_FLO)
rownames(NSC_S_FLO_L1) = rownames(NSC_S_FLO_K1) = rownames(NSC_S_FLO_M1) = names(Pars_S_FLO)

NSC_C_FLO <- cbind.data.frame(NSC_C_FLO_L, NSC_C_FLO_K, NSC_C_FLO_M,
                              NSC_C_FLO_L1, NSC_C_FLO_K1, NSC_C_FLO_M1)

NSC_S_FLO <- cbind.data.frame(NSC_S_FLO_L, NSC_S_FLO_K, NSC_S_FLO_M,
                              NSC_S_FLO_L1,NSC_S_FLO_K1, NSC_S_FLO_M1)


colnames(NSC_C_FLO) = colnames(NSC_S_FLO)=c("NSC_AUCL","NSC_AUCK","NSC_AUCM",
                                            "NSC_AUCL1","NSC_AUCK1","NSC_AUCM1")


write.csv(NSC_C_FLO, file = 'NSC_C_FLO.csv')
write.csv(NSC_S_FLO, file = 'NSC_S_FLO.csv')


