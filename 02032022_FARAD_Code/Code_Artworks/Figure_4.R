#------------------------------------------------------------------------------
# The code is to reproduce the figure 4
#------------------------------------------------------------------------------

## Loading required R package
library(FME)
library(mrgsolve)
library(dplyr)
library(minpack.lm)  ## R-package for model fitting
library(ggplot2)
library(tidyr)
library(patchwork)

## Input the PBPK model code and cirplot code
source (file = "GenPBPK.R") # Loading the generic PBPK model code
source (file = "CirPlot.R")  # Loading the circleplot code
source (file = 'Pars.R')  # Loading the parameters

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
# Label dose (mg/kg): 5 repeat dose 6.5 mg/kg; (lb = 0.454 kg) via IM injection for cattle and swine
# Extra-label dose: 5 repeat dose 32.5 mg/kg (5*label dose) via IM injection
# Tolerance of PG: 0.05 ug/g for tissues in cattle (USFDA); 
# Action limit of PG: 0.025 ug/g for swine (FSIS, 2013);
#-------------------------------------------------------------------------------


## Sensitivity function for Cattle
## Scenario: 6.5 mg/kg x 5 via IM injection

NSC_C_PG <- data.frame(NSC_func (log(Pars_C_PG), tinterval = 24, 
                                    tdose = 5, Dose = 6.5, AUC="AUCCM",
                                    route = "im", drug = 'PG'))

## Sensitivity function for Swine
## Scenario: 6.5 mg/kg x 5 via IM injection
NSC_S_PG <- data.frame(NSC_func (log(Pars_S_PG), tinterval = 24, 
                                    tdose = 5, Dose = 6.5, AUC="AUCCM",
                                    route = "im", drug = 'PG'))



rownames(NSC_C_PG) = names(Pars_C_PG)
rownames(NSC_S_PG) = names(Pars_S_PG)
colnames(NSC_C_PG) = colnames(NSC_S_PG) = "NSC"

NSC_C_PG<-NSC_C_PG %>% mutate(Pars = rownames(NSC_C_PG),
                              Drug = "PG",
                              Species = "C")

NSC_S_PG<-NSC_S_PG %>% mutate(Pars = rownames(NSC_S_PG),
                              Drug = "PG",
                              Species = "S")

unique (NSC_C_PG %>% filter (abs(NSC) >= 0.1))
unique (NSC_S_PG %>% filter (abs(NSC) >= 0.1))


# ------------------------------------------------------------------------------
# Flunixin (FLU) simulation scenario 
# Label dose for cattle: 2.2 mg/kg of 3 repeated IV injections
# Label dose for swine: 2.2 mg/kg of single IM injections
# Extra-label use for Cattle: 2.2 mg/kg of 3 repeated IM
# Extra-label use for Swine: 2.2  mg/kg of 3  repeated  IM  injections   
# Tolerance of FLU for Cattle: 0.125 (ug/g) for liver; 0.025 (ug/g) for muscle
# Tolerance of FLU for Swine: 0.03 (ug/g) for liver; 0.025 (ug/g) for muscle
# Dosing regimen: 2.2 mg/kg for 3 repeated dose via IM injection
#-------------------------------------------------------------------------------

## Sensitivity function for Cattle
## Scenario: 2.2 mg/kg x 3 via IV injection
NSC_Pars_C_FLU <- Pars_C_FLU [names(Pars_C_FLU) %in% c("Fracsc","Fracim")==FALSE]


NSC_C_FLU <- data.frame(NSC_func (log(NSC_Pars_C_FLU), tinterval = 24, 
                                    tdose = 3, Dose = 2.2, AUC="AUCCM",
                                    route = "iv", drug = 'FLU'))


## Sensitivity function for Swine
## Scenario: 2.2 mg/kg x 1 via IM injection
NSC_Pars_S_FLU <- Pars_S_FLU [names(Pars_S_FLU) %in% c("Fracsc","Fracim")==FALSE]

NSC_S_FLU <- data.frame(NSC_func (log(NSC_Pars_S_FLU), tinterval = 24, 
                                    tdose = 1, Dose = 2.2, AUC="AUCCM",
                                    route = "im", drug = 'FLU'))


rownames(NSC_C_FLU) = names(NSC_Pars_C_FLU)
rownames(NSC_S_FLU) = names(NSC_Pars_S_FLU)
colnames(NSC_S_FLU) = colnames(NSC_C_FLU) = "NSC"

NSC_C_FLU<-NSC_C_FLU %>% mutate(Pars = rownames(NSC_C_FLU),
                                Drug = "FLU",
                                Species = "C")

NSC_S_FLU<-NSC_S_FLU %>% mutate(Pars = rownames(NSC_S_FLU),
                                Drug = "FLU",
                                Species = "S")

unique (NSC_C_FLU %>% filter (abs(NSC) >= 0.1))
unique (NSC_S_FLU %>% filter (abs(NSC) >= 0.1))


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
# Dosing regimen: 20 mg/kg with 2 repeated IM injection and 48-hrs interval
#-------------------------------------------------------------------------------

## Sensitivity function for Cattle
## Scenario: 20 mg/kg x 2 via IM injection (48 hrs interval)
NSC_C_FLO <- data.frame(NSC_func (log(Pars_C_FLO), tinterval = 48, 
                                  tdose = 2, Dose = 20, AUC="AUCCM",
                                  route = "im", drug = 'FLO'))

## Sensitivity function for Swine
## Scenario (label use): 30 mg/kg x 5 via PO (24 hrs interval)
## 100 mg/l (in water) * 10 l/day (water consumption from Thacker (2001)) = 1000 mg/day
## 1000 mg/day / 70 kg (swine body weight) = 14 mg/kg/day
NSC_S_FLO <- data.frame(NSC_func (log(Pars_S_FLO), tinterval = 24, 
                                  tdose = 5, Dose = 14, AUC="AUCCM",
                                  route = "po", drug = 'FLO'))


rownames(NSC_C_FLO) = names(Pars_C_FLO)
rownames(NSC_S_FLO) = names(Pars_S_FLO)

colnames(NSC_C_FLO) = colnames(NSC_S_FLO) = "NSC"

NSC_C_FLO<-NSC_C_FLO %>% mutate(Pars = rownames(NSC_C_FLO),
                                Drug = "FLO",
                                Species = "C")

NSC_S_FLO<-NSC_S_FLO %>% mutate(Pars = rownames(NSC_S_FLO),
                                Drug = "FLO",
                                Species = "S")

unique (NSC_C_FLO %>% filter (abs(NSC) >= 0.1))
unique (NSC_S_FLO %>% filter (abs(NSC) >= 0.1))



#-------------------------------------------------------------------------------
## Make the circle plot
## Combine all NSC data sets
NSC <- rbind.data.frame(NSC_C_PG, NSC_C_FLU, NSC_C_FLO,
                        NSC_S_PG, NSC_S_FLU, NSC_S_FLO)

## Reorder the column
NSC <- NSC %>% relocate(c("Pars", "Drug", "Species"), .before = "NSC")

## 
NSC_H <- unique (NSC %>% group_by (Drug) %>%
                   filter (abs(NSC) >= 0.1))

Pdata_C <- NSC_H %>% filter(Species == "C") %>% mutate(value = abs(NSC)*100) %>% 
           select(Pars, Drug, value)

Pdata_C$Pars[c(4:7,21:26,36:40,41)]<-c("QKC","QMC","VLC","VRestC",
                                    "QLC","QKC","QMC","QRestC","VMC","VRestC",
                                    "QLC","QKC","QMC","VLC","VFC","VRestC")

Pdata_S <- NSC_H %>% filter(Species == "S") %>% mutate(value = abs(NSC)*100) %>% 
           select(Pars, Drug, value)

Pdata_S$Pars[c(3:5,14:16,25:32)]<-c("QLC","QKC","QMC",
                                    "QKC", "QMC","QRestC",
                                    "QLC","QKC","QMC", "QFC","QRestC",
                                    "VFC","VMC","VRest")



## Make plot
#Pdata_C<-Pdata_C%>%filter(!grepl("Frac", Pars))
#Pdata_S<-Pdata_S%>%filter(!grepl("Frac", Pars))

p1 <- Cirplot (Pdata_C) + scale_fill_brewer(palette = "Dark2") #scale_fill_manual(values=c("red", "blue", "green"))
p2 <- Cirplot (Pdata_S) + scale_fill_brewer(palette = "Dark2") #scale_fill_manual(values=c("red", "blue", "green"))

p1+p2

ggsave("Fig 4.tiff",scale = 1.1,
       plot = p1+p2,
       width = 30, height = 25, units = "cm", dpi=320)




