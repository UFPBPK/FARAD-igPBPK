#------------------------------------------------------------------------------
# The code is to reproduce the Figure 3
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
source (file = 'Pars.R')    # Loading the other functions used in this R script

## Load Model
mod <- mcode_cache("pbpk", GenricPBPK) #refer to mcode function in mrgsolve user guide 3.1.2 Inline


# ------------------------------------------------------------------------------
# Penicillin G (PG) simulation scenario
# Label dose (IU): 3000 units/lb bwt/day (https://www.etoolsage.com/converter/IU_Converter.asp)
# Label dose (mg/kg): 6.5 mg/kg; (lb = 0.454 kg)x 5 via IM injection
# Extra-label dose: 32.5 mg/kg (5*label dose) and 65 mg/kg
# Tolerance of PG: 0.05 ug/g for tissues in cattle (USFDA);
# Action limit of PG: 25 ng/g for swine (FSIS, 2013);
# LOD of PG: 1.8 ng/g in kidney; 0.7 in muscle for swine
# Dosing regimen: 3.96 for 5 repeated dose via IM injection
#-------------------------------------------------------------------------------

# Input dataset
Data_PG <- read.csv(file = "Data_PG.csv")
head(Data_PG)

## Study 7; Papich et al., 1993; IM, SC; P; Dose: 20, 65, 65 mg/kg/day
Obs_A1    <- Data_PG %>% filter(Study == 1)
Obs_A2    <- Data_PG %>% filter(Study == 2)
Obs_A3    <- Data_PG %>% filter(Study == 3)
Obs_A4    <- Data_PG %>% filter(Study == 4)
Obs_A5    <- Data_PG %>% filter(Study == 5)
Obs_A6    <- Data_PG %>% filter(Study == 6)
Obs_A7    <- Data_PG %>% filter(Study == 7)
Obs_A8    <- Data_PG %>% filter(Study == 8)
Obs_A9    <- Data_PG %>% filter(Study == 9)
Obs_A10   <- Data_PG %>% filter(Study == 10)
Obs_A11   <- Data_PG %>% filter(Study == 11)

Data_Cal <-rbind.data.frame(
    Obs_A1_1  <- Obs_A1 %>% filter(Route == "IM")%>%
                 select(Time = Time, Conc = Conc)%>%mutate(Study = 'A1_1', Species = 'Swine'),
    
    Obs_A1_2  <- Obs_A1 %>% filter(Route == "SC")%>%
                 select(Time = Time, Conc = Conc)%>%mutate(Study = 'A1_2', Species = 'Swine'),
    
    Obs_A2_1  <- Obs_A2 %>% filter(Dose == "14.9" & Matrix == 'P')%>%
                 select(Time = Time, Conc = Conc)%>%mutate(Study = 'A2_1', Species = 'Swine'),
    
    Obs_A2_2  <- Obs_A2 %>% filter(Dose == "14.9" & Matrix == 'K')%>%
                 select(Time = Time, Conc = Conc)%>%mutate(Study = 'A2_2', Species = 'Swine'),
    
    Obs_A2_3  <- Obs_A2 %>% filter(Dose == "14.9" & Matrix == 'M')%>%
                 select(Time = Time, Conc = Conc)%>%mutate(Study = 'A2_3', Species = 'Swine'),
    
    Obs_A2_4  <- Obs_A2 %>% filter(Dose == "65.4" & Matrix == 'P')%>%
                 select(Time = Time, Conc = Conc)%>%mutate(Study = 'A2_4', Species = 'Swine'),
    
    Obs_A2_5  <- Obs_A2 %>% filter(Dose == "65.4" & Matrix == 'L')%>%
                 select(Time = Time, Conc = Conc)%>%mutate(Study = 'A2_5', Species = 'Swine'),
    
    Obs_A2_6  <- Obs_A2 %>% filter(Dose == "65.4" & Matrix == 'K')%>%
                 select(Time = Time, Conc = Conc)%>%mutate(Study = 'A2_6', Species = 'Swine'),
    
    Obs_A2_7  <- Obs_A2 %>% filter(Dose == "65.4" & Matrix == 'M')%>%
                 select(Time = Time, Conc = Conc)%>%mutate(Study = 'A2_7', Species = 'Swine'),
    
    Obs_A3_1  <- Obs_A3 %>% filter(Matrix == "P")%>%
                 select(Time = Time, Conc = Conc)%>%mutate(Study = 'A3_1', Species = 'Swine'),
    
    Obs_A3_2  <- Obs_A3 %>% filter(Matrix == "M")%>%
                 select(Time = Time, Conc = Conc)%>%mutate(Study = 'A3_2', Species = 'Swine'),
    
    Obs_A3_3  <- Obs_A3 %>% filter(Matrix == "K")%>%
                 select(Time = Time, Conc = Conc)%>%mutate(Study = 'A3_3', Species = 'Swine'),
    
    # Obs_A4_1  <- Obs_A4 %>% filter(Matrix == 'P')%>%
    #              select(Time = Time, Conc = Conc)%>%mutate(Study = 'A4_1'),
    # 
    # Obs_A4_2  <- Obs_A4 %>% filter(Matrix == 'M')%>%
    #              select(Time = Time, Conc = Conc)%>%mutate(Study = 'A4_2'),
    # 
    # Obs_A4_3  <- Obs_A4 %>% filter(Matrix == 'K')%>%
    #              select(Time = Time, Conc = Conc)%>%mutate(Study = 'A4_3'),
    # 
    # Obs_A4_4  <- Obs_A4 %>% filter(Matrix == 'L')%>%
    #              select(Time = Time, Conc = Conc)%>%mutate(Study = 'A4_4'),
    
    Obs_A5_1  <- Obs_A5 %>% filter(Matrix == 'P')%>%
                 select(Time = Time, Conc = Conc)%>%mutate(Study = 'A5_1', Species = 'Swine'),
    
    Obs_A5_2  <- Obs_A5 %>% filter(Matrix == 'M')%>%
                 select(Time = Time, Conc = Conc)%>%mutate(Study = 'A5_2', Species = 'Swine'),
    
    Obs_A5_3  <- Obs_A5 %>% filter(Matrix == 'K')%>%
                 select(Time = Time, Conc = Conc)%>%mutate(Study = 'A5_3', Species = 'Swine'),
    
    Obs_A6_1  <- Obs_A6 %>% filter(Matrix == 'P' & Dose == 6.5)%>%
                 select(Time = Time, Conc = Conc)%>%mutate(Study = 'A6_1', Species = 'Swine'),
    
    Obs_A6_2  <- Obs_A6 %>% filter(Matrix == 'P' & Dose == 32.5)%>%
                 select(Time = Time, Conc = Conc)%>%mutate(Study = 'A6_2', Species = 'Swine'),
    
    Obs_A6_3  <- Obs_A6 %>% filter(Matrix == 'M' & Dose == 6.5)%>%
                 select(Time = Time, Conc = Conc)%>%mutate(Study = 'A6_3', Species = 'Swine'),
    
    Obs_A6_4  <- Obs_A6 %>% filter(Matrix == 'M' & Dose == 32.5)%>%
                 select(Time = Time, Conc = Conc)%>%mutate(Study = 'A6_4', Species = 'Swine'),
    
    Obs_A6_5  <- Obs_A6 %>% filter(Matrix == 'L' & Dose == 32.5)%>%
                 select(Time = Time, Conc = Conc)%>%mutate(Study = 'A6_5', Species = 'Swine'),
    
    Obs_A6_6  <- Obs_A6 %>% filter(Matrix == 'K' & Dose == 6.5)%>%
                 select(Time = Time, Conc = Conc)%>%mutate(Study = 'A6_6', Species = 'Swine'),
    
    Obs_A6_7  <- Obs_A6 %>% filter(Matrix == 'K' & Dose == 32.5)%>%
                 select(Time = Time, Conc = Conc)%>%mutate(Study = 'A6_7', Species = 'Swine'),
    
    Obs_A7_1  <- Obs_A7 %>% filter(Route == 'IM' & Dose == "23.7" & Repeat == "5")%>%
                 select(Time = Time, Conc = Conc)%>%mutate(Study = 'A7_1', Species = 'Cattle'),
   
    Obs_A7_2  <- Obs_A7 %>% filter(Route == 'IM' & Dose == "65.4" & Repeat == "5")%>%
                 select(Time = Time, Conc = Conc)%>%mutate(Study = 'A7_2', Species = 'Cattle'),
   
    Obs_A7_3  <- Obs_A7 %>% filter(Route == 'IM' & Dose == "65.4" & Repeat == "1")%>%
                 select(Time = Time, Conc = Conc)%>%mutate(Study = 'A7_3', Species = 'Cattle'),
   
    Obs_A7_4  <- Obs_A7 %>% filter(Route == 'SC'& Dose == "65.4")%>%
                 select(Time = Time, Conc = Conc)%>%mutate(Study = 'A7_4', Species = 'Cattle'),
   
    Obs_A8_1  <- Obs_A8 %>% filter(Route == 'IM' & Dose == "23.7" & Matrix == "P")%>%
                 select(Time = Time, Conc = Conc) %>% mutate(Study = 'A8_1', Species = 'Cattle'),
   
    Obs_A8_2  <- Obs_A8 %>% filter(Route == 'IM' & Dose == "23.7" & Matrix == "L")%>%
                 select(Time = Time, Conc = Conc)%>% mutate(Study = 'A8_2', Species = 'Cattle'),
   
    Obs_A8_3  <- Obs_A8 %>% filter(Route == 'IM' & Dose == "23.7" & Matrix == "K")%>%
                 select(Time = Time, Conc = Conc)%>% mutate(Study = 'A8_3', Species = 'Cattle'),
   
    Obs_A8_4  <- Obs_A8 %>% filter(Route == 'IM' & Dose == "23.7" & Matrix == "M")%>%
                select(Time = Time, Conc = Conc)%>% mutate(Study = 'A8_4', Species = 'Cattle'),
   
    Obs_A8_5  <- Obs_A8 %>% filter(Route == 'IM' & Dose == "65.4" & Matrix == "P")%>%
                select(Time = Time, Conc = Conc)%>% mutate(Study = 'A8_5', Species = 'Cattle'),
   
    Obs_A8_6  <- Obs_A8 %>% filter(Route == 'IM' & Dose == "65.4" & Matrix == "L")%>%
                select(Time = Time, Conc = Conc)%>% mutate(Study = 'A8_6', Species = 'Cattle'),
   
    Obs_A8_7  <- Obs_A8 %>% filter(Route == 'IM' & Dose == "65.4" & Matrix == "K")%>%
                select(Time = Time, Conc = Conc)%>% mutate(Study = 'A8_7', Species = 'Cattle'),
   
    Obs_A8_8  <- Obs_A8 %>% filter(Route == 'IM' & Dose == "65.4" & Matrix == "M")%>%
                select(Time = Time, Conc = Conc)%>% mutate(Study = 'A8_8', Species = 'Cattle'),
    
    Obs_A9_1  <- Obs_A9 %>% select(Time = Time, Conc = Conc)%>% mutate(Study = 'A9', Species = 'Cattle'),
   
    Obs_A10_1 <- Obs_A10 %>% filter(Matrix == 'P')%>%
                select(Time = Time, Conc = Conc)%>% mutate(Study = 'A10_1', Species = 'Cattle'),
   
    Obs_A10_2 <- Obs_A10 %>% filter(Matrix == 'K')%>%
                select(Time = Time, Conc = Conc)%>% mutate(Study = 'A10_2', Species = 'Cattle'),
    
    Obs_A11_1   <- Obs_A11 %>% select(Time = Time, Conc =  Conc)%>% mutate(Study = 'A11', Species = 'Cattle')
    
  )



## Calculation of idata
## PG
pred <- function (Vpars, BW, tinterval = 24, Dose, Dtimes, route) {
    
    ## Get out of log domain
    parsinput <- Vpars
    parsinput["BW"] <- BW
    Fracsc = ifelse(is.na(parsinput["Fracsc"]), 0.5, parsinput["Fracsc"]) 
    Fracim = ifelse(is.na(parsinput["Fracim"]), 0.5, parsinput["Fracim"]) 
    
    ## Exposure scenarios
    BW          = BW
    tinterval   = tinterval
    TDOSE       = Dtimes
    Frac        = switch(route, "im"= Fracim, "sc" = Fracsc)
    MW          = 334.4
    MW1         = 334.4
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
    
    
    tsamp  = tgrid(0, tinterval*(TDOSE - 1) + 24*30, 0.1) 
    
    ## Simulation
    out <- mod %>% param (parsinput) %>% 
        update(atol = 1E-6, rtol = 1E-3, maxsteps = 50000) %>%
        mrgsim_d (data = ex, tgrid = tsamp)
    
    
    outdf = cbind.data.frame( Time  = out$time,
                              CP    = out$Plasma*MW,
                              CL    = out$Liver*MW,
                              CK    = out$Kidney*MW,
                              CM    = out$Muscle*MW,
                              CF    = out$Fat*MW)
    return (outdf)
    
}

## Cattle
# Prediction

Out_Cal     <- rbind.data.frame(
out_A1_1    <- pred (Vpars = Pars_S_PG, BW = 3.3, Dose = 99, Dtimes = 1, route = "im")%>%
               select(Time = Time, Conc = CP)%>%mutate(Study = 'A1_1', Species = 'Swine'),

out_A1_2    <- pred (Vpars = Pars_S_PG, BW = 3.6, Dose = 99, Dtimes = 1, route = "sc")%>%
               select(Time = Time, Conc = CP)%>%mutate(Study = 'A1_2', Species = 'Swine'),

out_A2_1    <- pred (Vpars = Pars_S_PG, BW = 100, Dose = 14.9, Dtimes = 3, route = "im")%>%
               select(Time = Time, Conc = CP)%>%mutate(Study = 'A2_1', Species = 'Swine'),

out_A2_2    <- pred (Vpars = Pars_S_PG, BW = 100, Dose = 14.9, Dtimes = 3, route = "im")%>%
               select(Time = Time, Conc = CK)%>%mutate(Study = 'A2_2', Species = 'Swine'),

out_A2_3    <- pred (Vpars = Pars_S_PG, BW = 100, Dose = 14.9, Dtimes = 3, route = "im")%>%
               select(Time = Time, Conc = CM)%>%mutate(Study = 'A2_3', Species = 'Swine'),

out_A2_4    <- pred (Vpars = Pars_S_PG, BW = 100, Dose = 65.4, Dtimes = 5, route = "im")%>%
               select(Time = Time, Conc = CP)%>%mutate(Study = 'A2_4', Species = 'Swine'),

out_A2_5    <- pred (Vpars = Pars_S_PG, BW = 100, Dose = 65.4, Dtimes = 5, route = "im")%>%
               select(Time = Time, Conc = CL)%>%mutate(Study = 'A2_5', Species = 'Swine'),

out_A2_6    <- pred (Vpars = Pars_S_PG, BW = 100, Dose = 65.4, Dtimes = 5, route = "im")%>%
               select(Time = Time, Conc = CK)%>%mutate(Study = 'A2_6', Species = 'Swine'),

out_A2_7    <- pred (Vpars = Pars_S_PG, BW = 100, Dose = 65.4, Dtimes = 5, route = "im")%>%
               select(Time = Time, Conc = CM)%>%mutate(Study = 'A2_7', Species = 'Swine'),

out_A3_1    <- pred (Vpars = Pars_S_PG, BW = 100, Dose = 14.9, Dtimes = 3, route = "im")%>%
               select(Time = Time, Conc = CP)%>%mutate(Study = 'A3_1', Species = 'Swine'),

out_A3_2    <- pred (Vpars = Pars_S_PG, BW = 100, Dose = 14.9, Dtimes = 3, route = "im")%>%
               select(Time = Time, Conc = CM)%>%mutate(Study = 'A3_2', Species = 'Swine'),

out_A3_3    <- pred (Vpars = Pars_S_PG, BW = 100, Dose = 14.9, Dtimes = 3, route = "im")%>%
               select(Time = Time, Conc = CK)%>%mutate(Study = 'A3_3', Species = 'Swine'),

# out_A4_1    <- pred (Vpars = Pars_S_PG, BW = 209, Dose = 32.5, Dtimes = 1, route = "im")%>%
#                select(Time = Time, Conc = CM)%>%mutate(Study = 'A4_1'),
# 
# out_A4_2    <- pred (Vpars = Pars_S_PG, BW = 209, Dose = 32.5, Dtimes = 1, route = "im")%>%
#                select(Time = Time, Conc = CK)%>%mutate(Study = 'A4_2'),
# 
# out_A4_3    <- pred (Vpars = Pars_S_PG, BW = 209, Dose = 32.5, Dtimes = 1, route = "im")%>%
#                select(Time = Time, Conc = CL)%>%mutate(Study = 'A4_3'),
# 
# out_A4_4    <- pred (Vpars = Pars_S_PG, BW = 209, Dose = 32.5, Dtimes = 1, route = "im")%>%
#                select(Time = Time, Conc = CP)%>%mutate(Study = 'A4_4'),

out_A5_1    <- pred (Vpars = Pars_S_PG, BW = 228, Dose = 32.5, Dtimes = 3, route = "im")%>%
               select(Time = Time, Conc = CP)%>%mutate(Study = 'A5_1', Species = 'Swine'),

out_A5_2    <- pred (Vpars = Pars_S_PG, BW = 228, Dose = 32.5, Dtimes = 3, route = "im")%>%
               select(Time = Time, Conc = CM)%>%mutate(Study = 'A5_2', Species = 'Swine'),

out_A5_3    <- pred (Vpars = Pars_S_PG, BW = 228, Dose = 32.5, Dtimes = 3, route = "im")%>%
               select(Time = Time, Conc = CK)%>%mutate(Study = 'A5_3', Species = 'Swine'),

out_A6_1    <- pred (Vpars = Pars_S_PG, BW = 238, Dose = 6.5, Dtimes = 3, route = "im") %>%
               select(Time = Time, Conc = CP)%>%mutate(Study = 'A6_1', Species = 'Swine'),

out_A6_2    <- pred (Vpars = Pars_S_PG, BW = 238, Dose = 32.5, Dtimes = 3, route = "im")%>%
               select(Time = Time, Conc = CP)%>%mutate(Study = 'A6_2', Species = 'Swine'),

out_A6_3    <- pred (Vpars = Pars_S_PG, BW = 238, Dose = 6.5, Dtimes = 3, route = "im") %>%
               select(Time = Time, Conc = CM)%>%mutate(Study = 'A6_3', Species = 'Swine'),

out_A6_4    <- pred (Vpars = Pars_S_PG, BW = 238, Dose = 32.5, Dtimes = 3, route = "im")%>%
               select(Time = Time, Conc = CM)%>%mutate(Study = 'A6_4', Species = 'Swine'),

out_A6_5    <- pred (Vpars = Pars_S_PG, BW = 238, Dose = 32.5, Dtimes = 3, route = "im")%>%
               select(Time = Time, Conc = CL)%>%mutate(Study = 'A6_5', Species = 'Swine'),

out_A6_6    <- pred (Vpars = Pars_S_PG, BW = 238, Dose = 6.5, Dtimes = 3, route = "im") %>%
               select(Time = Time, Conc = CK)%>%mutate(Study = 'A6_6', Species = 'Swine'),

out_A6_7    <- pred (Vpars = Pars_S_PG, BW = 238, Dose = 32.5, Dtimes = 3, route = "im")%>%
               select(Time = Time, Conc = CK)%>%mutate(Study = 'A6_7', Species = 'Swine'),

out_A7_1    <- pred (Vpars = Pars_C_PG, BW = 480, Dose = 23.7, Dtimes = 5, route = "im")%>%
               select(Time = Time, Conc = CP)%>%mutate(Study = 'A7_1' , Species = 'Cattle'),

out_A7_2    <- pred (Vpars = Pars_C_PG, BW = 480, Dose = 65.4, Dtimes = 5, route = "im")%>%
               select(Time = Time, Conc = CP)%>%mutate(Study = 'A7_2' , Species = 'Cattle'),

out_A7_3    <- pred (Vpars = Pars_C_PG, BW = 480, Dose = 65.4, Dtimes = 1, route = "im")%>%
               select(Time = Time, Conc = CP)%>%mutate(Study = 'A7_3' , Species = 'Cattle'),

out_A7_4    <- pred (Vpars = Pars_C_PG, BW = 480, Dose = 65.4, Dtimes = 1, route = "sc")%>%
               select(Time = Time, Conc = CP)%>%mutate(Study = 'A7_4' , Species = 'Cattle'),

out_A8_1    <- pred (Vpars = Pars_C_PG, BW = 485, Dose = 23.7, Dtimes = 5, route = "im")%>%
               select(Time = Time, Conc = CP)%>%mutate(Study = 'A8_1' , Species = 'Cattle'),

out_A8_2    <- pred (Vpars = Pars_C_PG, BW = 485, Dose = 23.7, Dtimes = 5, route = "im")%>%
               select(Time = Time, Conc = CL)%>%mutate(Study = 'A8_2' , Species = 'Cattle'),

out_A8_3    <- pred (Vpars = Pars_C_PG, BW = 485, Dose = 23.7, Dtimes = 5, route = "im")%>%
               select(Time = Time, Conc = CK)%>%mutate(Study = 'A8_3' , Species = 'Cattle'),

out_A8_4    <- pred (Vpars = Pars_C_PG, BW = 485, Dose = 23.7, Dtimes = 5, route = "im")%>%
               select(Time = Time, Conc = CM)%>%mutate(Study = 'A8_4' , Species = 'Cattle'),

out_A8_5    <- pred (Vpars = Pars_C_PG, BW = 485, Dose = 65.4, Dtimes = 5, route = "im")%>%
               select(Time = Time, Conc = CP)%>%mutate(Study = 'A8_5' , Species = 'Cattle'),

out_A8_6    <- pred (Vpars = Pars_C_PG, BW = 485, Dose = 65.4, Dtimes = 5, route = "im")%>%
               select(Time = Time, Conc = CL)%>%mutate(Study = 'A8_6' , Species = 'Cattle'),

out_A8_7    <- pred (Vpars = Pars_C_PG, BW = 485, Dose = 65.4, Dtimes = 5, route = "im")%>%
               select(Time = Time, Conc = CK)%>%mutate(Study = 'A8_7' , Species = 'Cattle'),

out_A8_8    <- pred (Vpars = Pars_C_PG, BW = 485, Dose = 65.4, Dtimes = 5, route = "im")%>%
               select(Time = Time, Conc = CM)%>%mutate(Study = 'A8_8' , Species = 'Cattle'),

out_A9      <- pred (Vpars = Pars_C_PG, BW = 93, Dose = 8.9, Dtimes = 1, route = "sc")%>%
               select(Time = Time, Conc = CP)%>%mutate(Study = 'A9' , Species = 'Cattle'),

out_A10_1   <- pred (Vpars = Pars_C_PG, BW = 262, Dose = 6.9, Dtimes = 3, tinterval = 12, route = "im")%>%
               select(Time = Time, Conc = CP)%>%mutate(Study = 'A10_1' , Species = 'Cattle'),

out_A10_2   <- pred (Vpars = Pars_C_PG, BW = 262, Dose = 6.9, Dtimes = 3, tinterval = 12, route = "im")%>%
               select(Time = Time, Conc = CK)%>%mutate(Study = 'A10_2' , Species = 'Cattle'),

out_A11    <- pred (Vpars = Pars_C_PG, BW = 624, Dose = 20.8, Dtimes = 1, route = "im")%>%
              select(Time = Time, Conc = CP)%>%mutate(Study = 'A11' , Species = 'Cattle')

)



###/////////////////////////////////////////////////////////////////////////////
## Plot
Levels <- c('A1_1', 'A1_2','A2_1','A2_2','A2_3', 'A2_4', 
            'A2_5', 'A2_6','A2_7','A3_1', 'A3_2', 'A3_3', 
            'A5_1', 'A5_2','A5_3', 'A6_1','A6_2','A6_3', 
            'A6_4', 'A6_5','A6_6', 'A6_7','A7_1','A7_2', 
            'A7_3', 'A7_4','A8_1', 'A8_2','A8_3','A8_4', 
            'A8_5', 'A8_6','A8_7', 'A8_8','A9_1','A10_1',
            'A10_2','A11')

Data_Cal$Study <- factor(Data_Cal$Study, levels=Levels)
Out_Cal$Study <- factor(Out_Cal$Study, levels=Levels)

p<-ggplot(data = Data_Cal, aes(x=Time, y=Conc)) + 
    geom_point(shape = 21, colour = "red", fill = "white", size = 1.5, stroke = 1.2) + 
    geom_line(data = data.frame(Out_Cal), aes(x = Time, y = Conc,linetype = Species), 
              size = 0.6, colour = "black") +
   
    scale_y_log10 (labels = function(x) format(x, scientific = TRUE))  + 
    facet_wrap(~Study, scales = "free",ncol=6) +
    theme_bw() +
    theme (
        panel.background = element_rect(fill = "white"),
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
ggsave("Fig 3.tiff",scale = 1.5,
       plot = p,
       width = 25, height = 15, units = "cm", dpi=320)


















