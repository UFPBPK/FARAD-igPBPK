#------------------------------------------------------------------------------
# The code is to reproduce the figure S3
# Note: Prior to running this code, the code for figure 3 should be run
#------------------------------------------------------------------------------

## Source code for Monte Carlo simulation
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

# ------------------------------------------------------------------------------
# Florfenicol (FLO) and simulation scenario 
# Label dose: 
# Extra-label dose: 
# Tolerance in US (FLOA) in cattle: 3.7 (ug/g) in liver; 0.3 (ug/g) in muscle (FDA, 2017)
# Tolerance in China and EU (FLO+FLOA): 3 (ug/g) in liver; 0.2 (ug/g) in muscle; 0.3 (ug/g) in kidney;
# MRLs for FLO+FLOA in swine in US: 2.5 (ug/g) in liver; 0.2 (ug/g) in muscle;
# MRLs for FLO+FLOA in swine in China and EU: 2 (ug/g) in liver; 0.3 (ug/g) in muscle; 0.5 (ug/g) in Kidney
# Dosing regimen: 20 mg/kg with 2 repeated IM injection and 48-hrs interval
#-------------------------------------------------------------------------------
## Calculation of idata
## FLO
Data_FLO <- read.csv(file = "Data_FLO.csv")
head(Data_FLO)


## Read the dataset and later used in model calibration
Obs_A1  <- Data_FLO %>% filter(Study ==1)
Obs_A2  <- Data_FLO %>% filter(Study ==2)
Obs_A3  <- Data_FLO %>% filter(Study ==3)
Obs_A4  <- Data_FLO %>% filter(Study ==4)
Obs_A5  <- Data_FLO %>% filter(Study ==5)
Obs_A6  <- Data_FLO %>% filter(Study ==6)
Obs_A7  <- Data_FLO %>% filter(Study ==7)
Obs_A8  <- Data_FLO %>% filter(Study ==8)
Obs_A9  <- Data_FLO %>% filter(Study ==9)
Obs_A10 <- Data_FLO %>% filter(Study ==10)
Obs_A11 <- Data_FLO %>% filter(Study ==11)
Obs_A12 <- Data_FLO %>% filter(Study ==12)
Obs_A13 <- Data_FLO %>% filter(Study ==13)
Obs_A14 <- Data_FLO %>% filter(Study ==14)
Obs_A15 <- Data_FLO %>% filter(Study ==15)
Obs_A16 <- Data_FLO %>% filter(Study ==16)
Obs_A17 <- Data_FLO %>% filter(Study ==17)
Obs_A18 <- Data_FLO %>% filter(Study ==18)
Obs_A19 <- Data_FLO %>% filter(Study ==19)
Obs_A20 <- Data_FLO %>% filter(Study ==20)
Obs_A21 <- Data_FLO %>% filter(Study ==21)
Obs_A22 <- Data_FLO %>% filter(Study ==22)


Data_Cal <-rbind.data.frame(
    Obs_A1_1 <- Obs_A1 %>% filter(Dose == "22.5")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A1_1'),
    Obs_A1_2 <- Obs_A1 %>% filter(Dose == "30")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A1_2'),
    Obs_A1_3 <- Obs_A1 %>% filter(Dose == "15")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A1_3'),
    Obs_A2_1 <- Obs_A2 %>% filter(Route == 'IV')%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A2_1'),
    Obs_A2_2 <- Obs_A2 %>% filter(Route == 'IM')%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A2_2'),
    Obs_A2_3 <- Obs_A2 %>% filter(Route == 'PO')%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A2_3'),
    Obs_A3_1 <- Obs_A3 %>% filter(Dose == 5)%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A3_1'),
    Obs_A3_2 <- Obs_A3 %>% filter(Dose == 20)%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A3_2'),
    Obs_A4_1 <- Obs_A4 %>% filter(Matrix == 'Plasma')%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A4_1'),
    Obs_A4_2 <- Obs_A4 %>% filter(Matrix == 'Liver')%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A4_2'),
    Obs_A4_3 <- Obs_A4 %>% filter(Matrix == 'Kidney')%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A4_3'),
    Obs_A4_4 <- Obs_A4 %>% filter(Matrix == 'Muscle')%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A4_4'),
    Obs_A5_1 <- Obs_A5 %>% filter(Route == 'IV')%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A5_1'),
    Obs_A5_2 <- Obs_A5 %>% filter(Route == 'IM')%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A5_2'),
    Obs_A5_3 <- Obs_A5 %>% filter(Route == 'PO')%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A5_3'),
    Obs_A6_1 <- Obs_A6 %>% filter(Matrix == 'Liver')%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A6_1'),
    Obs_A6_2 <- Obs_A6 %>% filter(Matrix == 'Kidney')%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A6_2'),
    Obs_A6_3 <- Obs_A6 %>% filter(Matrix == 'Muscle')%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A6_3'),
    #Obs_A6_4 <- Obs_A6 %>% filter(Matrix == 'Fat')%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A6_4'),
    Obs_A7_1 <- Obs_A7 %>% filter(Matrix == 'Liver')%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A7_1'),
    Obs_A7_2 <- Obs_A7 %>% filter(Matrix == 'Kidney')%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A7_2'),
    Obs_A7_3 <- Obs_A7 %>% filter(Matrix == 'Muscle')%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A7_3'),
    #Obs_A7_4 <- Obs_A7 %>% filter(Matrix == 'Fat')%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A7_4'),
    Obs_A8_1 <- Obs_A8 %>% filter(Route == 'IM')%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A8_1'),
    Obs_A8_2 <- Obs_A8 %>% filter(Route == 'PO' & Dose == 250)%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A8_2'),
    Obs_A9_1 <- Obs_A9 %>% filter(Route == 'IM')%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A9_1'),
    Obs_A10_1 <- Obs_A10 %>% filter(Route == 'IV')%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A10_1'),
    Obs_A10_2 <- Obs_A10 %>% filter(Route == 'IM')%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A10_2'),
    Obs_A11_1 <- Obs_A11 %>% filter(Route == "IV")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A11_1'),
    Obs_A11_2 <- Obs_A11 %>% filter(Route == "PO")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A11_2'),
    Obs_A12_1 <- Obs_A12 %>% select(Time = Time, Conc = Conc)%>%mutate(Study = 'A12_1'),
    Obs_A13_1 <- Obs_A13 %>% filter(Route == "IV")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A13_1'),
    Obs_A13_2 <- Obs_A13 %>% filter(Route == "IM")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A13_2'),
    Obs_A14_1 <- Obs_A14 %>% filter(Route == "IM")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A14_1'),
    Obs_A14_2 <- Obs_A14 %>% filter(Route == "SC")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A14_2'),
    #Obs_A15_1 <- Obs_A15 %>% filter(Route == "SC")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A15_1'),
    Obs_A16_1 <- Obs_A16 %>% filter(Route == "IV")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A16_1'),
    Obs_A17_1 <- Obs_A17 %>% filter(Route == "IV")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A17_1'),
    Obs_A18_1 <- Obs_A18 %>% filter(Route == "PO")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A18_1'),
    Obs_A19_1 <- Obs_A19 %>% filter(Route == "IV")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A19_1'),
    Obs_A20_1 <- Obs_A20 %>% filter(Matrix == "Liver") %>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A20_1'),
    Obs_A20_2 <- Obs_A20 %>% filter(Matrix == "Kidney")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A20_2'),
    Obs_A20_3 <- Obs_A20 %>% filter(Matrix == "Muscle")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A20_3'),
    Obs_A21_1 <- Obs_A21 %>% filter(Matrix == "Liver") %>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A21_1'),
    Obs_A21_2 <- Obs_A21 %>% filter(Matrix == "Muscle")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A21_2'),
    Obs_A22_1 <- Obs_A22 %>% filter(Matrix == "Liver")%>%select(Time = Time, Conc = Conc)%>%mutate(Study = 'A22_1')
    
)
    
##
pred <- function (Vpars, BW, tinterval = 24, Dose, Dtimes, route) {
    
    ## Get out of log domain
    parsinput <- Vpars
    parsinput["BW"] <- BW
    Fracsc = ifelse(is.na(parsinput["Fracsc"]), 1, parsinput["Fracsc"]) 
    Fracim = ifelse(is.na(parsinput["Fracim"]), 1, parsinput["Fracim"]) 
    
    ## Exposure scenarios
    BW          = BW
    tinterval   = tinterval
    TDOSE       = Dtimes
    Frac        = switch(route, "im"= Fracim, "sc" = Fracsc)
    MW          = 296.24
    MW1         = 312.24
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
    
    if (route == "oral") {
        ev_1 <- ev (ID   = 1, amt  = DOSE, ii = tinterval, 
                    addl = TDOSE - 1, cmt  = "AST", replicate = FALSE)
        ev_2 <- ev (ID   = 1, amt  = DOSE, ii = tinterval, 
                    addl = TDOSE - 1, cmt  = "ADOSE", replicate = FALSE)
        ex <- ev_1+ev_2
    }
    
    
    tsamp  = tgrid(0, tinterval*(TDOSE - 1) + 24*60, 0.1) 
    
    ## Simulation
    out <- mod %>% param (parsinput) %>% 
        update(atol = 1E-10, rtol = 1E-5, maxsteps = 50000) %>%
        mrgsim_d (data = ex, tgrid = tsamp)
    
    
    outdf = cbind.data.frame( Time  = out$time,
                              CP    = out$Plasma*MW,
                              CP_Met= out$P_Met*MW1,
                              CL    = out$Liver*MW,
                              CL_Met= out$L_Met*MW1,
                              CK    = out$Kidney*MW,
                              CK_Met= out$K_Met*MW1,
                              CM    = out$Muscle*MW,
                              CM_Met= out$M_Met*MW1,
                              CF    = out$Fat*MW)
    return (outdf)
    
}

## Cattle
# Prediction
Out_Cal     <- rbind.data.frame(
    
    # Prediction
    out_A1_1    <- pred (Vpars = Pars_S_FLO, BW = 36.5, Dose = 22.5, Dtimes = 1, route = "im")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A1_1'),
    out_A1_2    <- pred (Vpars = Pars_S_FLO, BW = 36.5, Dose = 30, Dtimes = 1, route = "im")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A1_2'),
    out_A1_3    <- pred (Vpars = Pars_S_FLO, BW = 36.5, Dose = 15, Dtimes = 2, tinterval = 48,route = "im")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A1_3'),
    out_A2_1    <- pred (Vpars = Pars_S_FLO, BW = 27, Dose = 20, Dtimes = 1, route = "iv")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A2_1'),
    out_A2_2    <- pred (Vpars = Pars_S_FLO, BW = 27, Dose = 20, Dtimes = 1, route = "im")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A2_2'),
    out_A2_3    <- pred (Vpars = Pars_S_FLO, BW = 27, Dose = 20, Dtimes = 1, route = "oral")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A2_3'),
    out_A3_1    <- pred (Vpars = Pars_S_FLO, BW = 36.5, Dose = 5, Dtimes = 1, route = "im")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A3_1'),
    out_A3_2    <- pred (Vpars = Pars_S_FLO, BW = 36.5, Dose = 20, Dtimes = 1, route = "im")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A3_2'),
    out_A4_1    <- pred (Vpars = Pars_S_FLO, BW = 36.5, Dose = 20, Dtimes = 1, route = "im")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A4_1'),
    out_A4_2    <- pred (Vpars = Pars_S_FLO, BW = 36.5, Dose = 20, Dtimes = 1, route = "im")%>%
                   select(Time = Time, Conc = CL)%>%mutate(Study = 'A4_2'),
    out_A4_3    <- pred (Vpars = Pars_S_FLO, BW = 36.5, Dose = 20, Dtimes = 1, route = "im")%>%
                   select(Time = Time, Conc = CK)%>%mutate(Study = 'A4_3'),
    out_A4_4    <- pred (Vpars = Pars_S_FLO, BW = 36.5, Dose = 20, Dtimes = 1, route = "im")%>%
                   select(Time = Time, Conc = CM)%>%mutate(Study = 'A4_4'),
    out_A5_1    <- pred (Vpars = Pars_S_FLO, BW = 36.5, Dose = 20, Dtimes = 1, route = "iv")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A5_1'),
    out_A5_2    <- pred (Vpars = Pars_S_FLO, BW = 36.5, Dose = 20, Dtimes = 1, route = "im")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A5_2'),
    out_A5_3    <- pred (Vpars = Pars_S_FLO, BW = 36.5, Dose = 20, Dtimes = 1, route = "oral")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A5_3'),
    out_A6_1    <- pred (Vpars = Pars_S_FLO, BW = 59, Dose = 20, Dtimes = 5, route = "oral")%>%
                   select(Time = Time, Conc = CL_Met)%>%mutate(Study = 'A6_1'),
    out_A6_2    <- pred (Vpars = Pars_S_FLO, BW = 59, Dose = 20, Dtimes = 5, route = "oral")%>%
                   select(Time = Time, Conc = CK_Met)%>%mutate(Study = 'A6_2'),
    out_A6_3    <- pred (Vpars = Pars_S_FLO, BW = 59, Dose = 20, Dtimes = 5, route = "oral")%>%
                   select(Time = Time, Conc = CM_Met)%>%mutate(Study = 'A6_3'),
    # out_A6_4    <- pred (Vpars = Pars_S_FLO, BW = 59, Dose = 20, Dtimes = 5, route = "oral")%>%
    #                select(Time = Time, Conc = CF_Met)%>%mutate(Study = 'A6_4'),
    out_A7_1    <- pred (Vpars = Pars_S_FLO, BW = 59, Dose = 10, Dtimes = 5, route = "oral")%>%
                   select(Time = Time, Conc = CL_Met)%>%mutate(Study = 'A7_1'),
    out_A7_2    <- pred (Vpars = Pars_S_FLO, BW = 59, Dose = 10, Dtimes = 5, route = "oral")%>%
                   select(Time = Time, Conc = CK_Met)%>%mutate(Study = 'A7_2'),
    out_A7_3    <- pred (Vpars = Pars_S_FLO, BW = 59, Dose = 10, Dtimes = 5, route = "oral")%>%
                   select(Time = Time, Conc = CM_Met)%>%mutate(Study = 'A7_3'),
    # out_A7_4    <- pred (Vpars = Pars_S_FLO, BW = 59, Dose = 10, Dtimes = 5, route = "oral")%>%
    #                select(Time = Time, Conc = CF_Met)%>%mutate(Study = 'A7_4'),
    out_A8_1    <- pred (Vpars = Pars_S_FLO, BW = 40, Dose = 15, Dtimes = 1, route = "im")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A8_1'),
    out_A8_2    <- pred (Vpars = Pars_S_FLO, BW = 40, Dose = 15, Dtimes = 3, route = "oral")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A8_2'),
    out_A9_1    <- pred (Vpars = Pars_S_FLO, BW = 33.6, Dose = 20, Dtimes = 1, route = "im")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A9_1'),
    out_A10_1   <- pred (Vpars = Pars_S_FLO, BW = 20, Dose = 30, Dtimes = 1, route = "iv")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A10_1'),
    out_A10_2   <- pred (Vpars = Pars_S_FLO, BW = 20, Dose = 30, Dtimes = 1, route = "im")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A10_2'),
    out_A11_1   <- pred (Vpars = Pars_C_FLO,  BW = 64, Dose = 22, Dtimes = 1, route = "iv")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A11_1'),
    out_A11_2   <- pred (Vpars = Pars_C_FLO, BW = 64, Dose = 22, Dtimes = 1, route = "oral")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A11_2'),
    out_A12_1   <- pred (Vpars = Pars_C_FLO, BW = 204, Dose = 40, Dtimes = 1, route = "sc")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A12_1'),
    out_A13_1   <- pred (Vpars = Pars_C_FLO, BW = 113, Dose = 20, Dtimes = 1, route = "iv")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A13_1'),
    out_A13_2   <- pred (Vpars = Pars_C_FLO, BW = 113, Dose = 20, Dtimes = 1, route = "im")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A13_2'),
    out_A14_1   <- pred (Vpars = Pars_C_FLO, BW = 563, Dose = 40, Dtimes = 1, route = "im")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A14_1'),
    out_A14_2   <- pred (Vpars = Pars_C_FLO, BW = 563, Dose = 40, Dtimes = 1, route = "sc")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A14_2'),
    # out_A15_1   <- pred (Vpars = Pars_C_FLO, BW = 700, Dose = 40, Dtimes = 1, route = "sc")%>%
    #                select(Time = Time, Conc = CP)%>%mutate(Study = 'A15_1'),
    out_A16_1   <- pred (Vpars = Pars_C_FLO, BW = 688, Dose = 2.2, Dtimes = 1, route = "iv")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A16_1'),
    out_A17_1   <- pred (Vpars = Pars_C_FLO, BW = 290, Dose = 20, Dtimes = 1, route = "iv")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A17_1'),
    out_A18_1   <- pred (Vpars = Pars_C_FLO, BW = 79.5, Dose = 20, Dtimes = 1, route = "oral")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A18_1'),
    out_A19_1   <- pred (Vpars = Pars_C_FLO, BW = 550, Dose = 50, Dtimes = 1, route = "iv")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A19_1'),
    out_A20_1   <- pred (Vpars = Pars_C_FLO, BW = 228, Dose = 40, Dtimes = 1, route = "sc")%>%
                   select(Time = Time, Conc = CL_Met)%>%mutate(Study = 'A20_1'),
    out_A20_2   <- pred (Vpars = Pars_C_FLO, BW = 228, Dose = 40, Dtimes = 1, route = "sc")%>%
                   select(Time = Time, Conc = CK_Met)%>%mutate(Study = 'A20_2'),
    out_A20_3   <- pred (Vpars = Pars_C_FLO, BW = 228, Dose = 40, Dtimes = 1, route = "sc")%>%
                   select(Time = Time, Conc = CM_Met)%>%mutate(Study = 'A20_3'),
    out_A21_1   <- pred (Vpars = Pars_C_FLO, BW = 302, Dose = 40, Dtimes = 1, route = "sc")%>%
                   select(Time = Time, Conc = CL_Met)%>%mutate(Study = 'A21_1'),
    out_A21_2   <- pred (Vpars = Pars_C_FLO, BW = 302, Dose = 40, Dtimes = 1, route = "sc")%>%
                   select(Time = Time, Conc = CM_Met)%>%mutate(Study = 'A21_2'),
    out_A22_1   <- pred (Vpars = Pars_C_FLO, BW = 292, Dose = 40, Dtimes = 1, route = "sc")%>%
                   select(Time = Time, Conc = CL_Met)%>%mutate(Study = 'A22_1')
    
)
    
    
###/////////////////////////////////////////////////////////////////////////////
## Plot
Levels <- c('A1_1', 'A1_2', 'A1_3', 'A2_1','A2_2', 
            'A2_3','A3_1', 'A3_2','A4_1', 'A4_2',
            'A4_3', 'A4_4', 'A5_1', 'A5_2','A5_3', 
            'A6_1','A6_2', 'A6_3', 'A7_1','A7_2', 
            'A7_3', 'A8_1', 'A8_2','A9_1', 'A10_1',
            'A10_2', 'A11_1','A11_2', 'A12_1','A13_1', 
            'A13_2', 'A14_1', 'A14_2', 'A16_1',
            'A17_1','A18_1','A19_1','A20_1','A20_2',
            'A20_3','A21_1','A21_2','A22_1')

Data_Cal$Study <- factor(Data_Cal$Study, levels=Levels)
Out_Cal$Study <- factor(Out_Cal$Study, levels=Levels)

p<-ggplot(data = Data_Cal, aes(x=Time, y=Conc)) + 
    geom_point(shape = 21, colour = "black", fill = "white", size = 1.5, stroke = 1.2) + 
    geom_line(data = data.frame(Out_Cal), aes(x = Time, y = Conc ), 
              linetype = 2, size = 0.6, colour = "black") +
    
    scale_y_log10 (labels = function(x) format(x, scientific = TRUE),limits = c(0.0001, NA))  + 
    facet_wrap(~Study, scales = "free",ncol=5) +
    theme_bw() +
    theme (
        plot.background         = element_rect (fill="White"),
        text                    = element_text (family = "Times"),   # text front (Time new roman)
        strip.background        = element_blank(),
        strip.text              = element_blank(),
        panel.grid.major        = element_blank(),
        panel.grid.minor        = element_blank(),
        axis.text               = element_text (size   = 10, colour = "black", face = "bold"),    # tick labels along axes
        axis.title              = element_text (size   = 10, colour = "black", face = "bold"),   # label of axes
        legend.position         ='none') +
    labs (x = "",  y = "")

    
    
## Export the figure
ggsave("Fig S2.tiff",scale = 1.5,
       plot = p,
       width = 25, height = 15, units = "cm", dpi=320)





