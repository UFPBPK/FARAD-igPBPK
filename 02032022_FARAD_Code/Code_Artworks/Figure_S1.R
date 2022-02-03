#------------------------------------------------------------------------------
# The code is to reproduce the figure S1
# Note: Prior to running this code, the code for figure 3 should be run
#------------------------------------------------------------------------------

## Source code for Monte Carlo simulation
library(EnvStats)
library(readxl)     ## R-package for importing excel file to R
library(FME)        ## R-package for the function of 'modFit' and 'modCost'
library(mrgsolve)
library(dplyr)
library(minpack.lm)  ## R-package for model fitting
library(ggplot2)
library(tidyr)
library(ggprism)    ## R-package for theme of 'GraphPad Prism'
library(patchwork)  ## R-package for combination of seperate ggplot

# Input the PBPK model
source (file = "GenPBPK.R") # Loading the generic PBPK model code
source (file = "GetPC.R")   # Loading the tissue composition based tissue:plasma model
source (file = 'Pars.R')  # Loading the parameters

## Load Model
mod <- mcode_cache("pbpk", GenricPBPK) #refer to mcode function in mrgsolve user guide 3.1.2 Inline



# ------------------------------------------------------------------------------
# Flunixin (FLU) simulation scenario 
# Label dose: 2.2 mg/kg of 3 repeated IV injections
# Extra-label use for Cattle: 2.2 mg/kg of 3 repeated IM
# Extra-label use for Swine: 2.2  mg/kg (single) or  3  repeated  IM  injections   
# Tolerance of FLU for Cattle: 0.125 (ug/g) for liver; 0.025 (ug/g) for muscle
# Tolerance of FLU for Swine: 0.03 (ug/g) for liver; 0.025 (ug/g) for muscle
# Dosing regimen: 2.2 mg/kg for 3 repeated dose via IM injection
#-------------------------------------------------------------------------------
## Calculation of idata
## FLU
# Input dataset
Data_FLU<- read.csv(file = "Data_FLU.csv")
head(Data_FLU)

## Read the dataset and later used in model calibration
Obs_A1   <- Data_FLU %>% filter(Study == 1) 
Obs_A2   <- Data_FLU %>% filter(Study == 2)
Obs_A3   <- Data_FLU %>% filter(Study == 3)
Obs_A4   <- Data_FLU %>% filter(Study == 4)
Obs_A5   <- Data_FLU %>% filter(Study == 5)
Obs_A6   <- Data_FLU %>% filter(Study == 6)
Obs_A7   <- Data_FLU %>% filter(Study == 7) 
Obs_A8   <- Data_FLU %>% filter(Study == 8)
Obs_A9   <- Data_FLU %>% filter(Study == 9)
Obs_A10  <- Data_FLU %>% filter(Study == 10)
Obs_A11  <- Data_FLU %>% filter(Study == 11)
Obs_A12  <- Data_FLU %>% filter(Study == 12)
Obs_A13  <- Data_FLU %>% filter(Study == 13)
Obs_A14  <- Data_FLU %>% filter(Study == 14)

Data_Cal <-rbind.data.frame(
    
    Obs_A1_1 <- Obs_A1%>%filter(Matrix == "P")%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A1_1'),
    Obs_A1_2 <- Obs_A1%>%filter(Matrix == "P-Met")%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A1_2'),
    Obs_A2_1 <- Obs_A2%>%filter(Route == "IV")%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A2_1'),
    Obs_A2_2 <- Obs_A2%>%filter(Route == "IM")%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A2_2'),
    Obs_A2_3 <- Obs_A2%>%filter(Route == "PO")%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A2_3'),
    Obs_A3_1 <- Obs_A3%>%filter(Matrix == "L")%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A3_1'),
    Obs_A3_2 <- Obs_A3%>%filter(Matrix == "K")%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A3_2'),
    Obs_A3_3 <- Obs_A3%>%filter(Matrix == "M")%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A3_3'),
    Obs_A3_4 <- Obs_A3%>%filter(Matrix == "F")%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A3_4'),
    Obs_A4_1 <- Obs_A4%>%filter(Matrix == "P")%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A4_1'),
    Obs_A5_1 <- Obs_A5%>%filter(Matrix == "L")%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A5_1'),
    Obs_A5_2 <- Obs_A5%>%filter(Matrix == "M")%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A5_2'),
    Obs_A5_3 <- Obs_A5%>%filter(Matrix == "K")%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A5_3'),
    Obs_A6_1 <- Obs_A6%>%filter(Matrix == "P"&Comp. == 'FLU')%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A6_1'),
    Obs_A6_2 <- Obs_A6%>%filter(Matrix == "P"&Comp. == '5OH')%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A6_2'),
    Obs_A6_3 <- Obs_A6%>%filter(Matrix == "L"&Comp. == 'FLU')%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A6_3'),
    Obs_A6_4 <- Obs_A6%>%filter(Matrix == "L"&Comp. == '5OH')%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A6_4'),
    Obs_A6_5 <- Obs_A6%>%filter(Matrix == "K"&Comp. == 'FLU')%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A6_5'),
    Obs_A6_6 <- Obs_A6%>%filter(Matrix == "K"&Comp. == '5OH')%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A6_6'),
    Obs_A7_1 <- Obs_A7%>% filter(Route == "IM")%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A7_1'),
    Obs_A7_2 <- Obs_A7%>%filter(Route == "PO")%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A7_2'),
    Obs_A8_1 <- Obs_A8%>%filter(Matrix == "L")%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A8_1'),
    Obs_A8_2 <- Obs_A8%>%filter(Matrix == "M")%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A8_2'),
    Obs_A8_3 <- Obs_A8%>%filter(Matrix == "K")%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A8_3'),
    Obs_A8_4 <- Obs_A8%>%filter(Matrix == "F")%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A8_4'),
    Obs_A9_1 <- Obs_A9%>%filter(Route == "IV")%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A9_1'),
    Obs_A9_2 <- Obs_A9%>%filter(Route == "SC")%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A9_2'),
    Obs_A10_1 <- Obs_A10%>% filter(Route == "IV")%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A10_1'),
    Obs_A10_2 <- Obs_A10%>% filter(Route == "IM")%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A10_2'),
    Obs_A11_1 <- Obs_A11%>% filter(Matrix == "P")%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A11_1'),
    Obs_A11_2 <- Obs_A11%>% filter(Matrix == "L")%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A11_2'),
    Obs_A11_3 <- Obs_A11%>% filter(Matrix == "M")%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A11_3'),
    Obs_A11_4 <- Obs_A11%>% filter(Matrix == "K")%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A11_4'),
    Obs_A11_5 <- Obs_A11%>% filter(Matrix == "L-Met")%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A11_5'),
    Obs_A12_1 <- Obs_A12%>% filter(Study ==12 & Matrix == "P-Met")%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A12_1'),
    Obs_A13_1 <- Obs_A13 %>% filter(Study ==13)%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A13_1'),
    Obs_A14_1 <- Obs_A14 %>% filter(Study ==14)%>%select(Time = Time, Conc = Conc.)%>%mutate(Study = 'A14_1')
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
    
    
    tsamp  = tgrid(0, tinterval*(TDOSE - 1) + 24*30, 0.1) 
    
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
    out_A1_1    <- pred (Vpars = Pars_S_FLU, BW = 40.15, Dose = 3,   Dtimes = 1, route = "iv")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A1_1'),
    out_A1_2    <- pred (Vpars = Pars_S_FLU, BW = 40.15, Dose = 3,   Dtimes = 1, route = "iv")%>%
                   select(Time = Time, Conc = CP_Met)%>%mutate(Study = 'A1_2'),
    out_A2_1    <- pred (Vpars = Pars_S_FLU, BW = 168,   Dose = 2.2, Dtimes = 1, route = "iv")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A2_1'),
    out_A2_2    <- pred (Vpars = Pars_S_FLU, BW = 168,   Dose = 2.2, Dtimes = 1, route = "im")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A2_2'),
    out_A2_3    <- pred (Vpars = Pars_S_FLU, BW = 168,   Dose = 2.2, Dtimes = 1, route = "oral")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A2_3'),
    out_A3_1    <- pred (Vpars = Pars_S_FLU, BW = 40,    Dose = 2.2, Dtimes = 3, route = "im")%>%
                   select(Time = Time, Conc = CL)%>%mutate(Study = 'A3_1'),
    out_A3_2    <- pred (Vpars = Pars_S_FLU, BW = 40,    Dose = 2.2, Dtimes = 3, route = "im")%>%
                   select(Time = Time, Conc = CK)%>%mutate(Study = 'A3_2'),
    out_A3_3    <- pred (Vpars = Pars_S_FLU, BW = 40,    Dose = 2.2, Dtimes = 3, route = "im")%>%
                   select(Time = Time, Conc = CM)%>%mutate(Study = 'A3_3'),
    out_A3_4    <- pred (Vpars = Pars_S_FLU, BW = 40,    Dose = 2.2, Dtimes = 3, route = "im")%>%
                   select(Time = Time, Conc = CF)%>%mutate(Study = 'A3_4'),
    out_A4_1    <- pred (Vpars = Pars_S_FLU, BW = 26.5,  Dose = 2,   Dtimes = 1, route = "iv")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A4_1'),
    out_A5_1    <- pred (Vpars = Pars_S_FLU, BW = 26.5,  Dose = 2.4, Dtimes = 3, route = "im")%>%
                   select(Time = Time, Conc = CL)%>%mutate(Study = 'A5_1'),
    out_A5_2    <- pred (Vpars = Pars_S_FLU, BW = 26.5,  Dose = 2.4, Dtimes = 3, route = "im")%>%
                   select(Time = Time, Conc = CM)%>%mutate(Study = 'A5_2'),
    out_A5_3    <- pred (Vpars = Pars_S_FLU, BW = 26.5,  Dose = 2.4, Dtimes = 3, route = "im")%>%
                   select(Time = Time, Conc = CK)%>%mutate(Study = 'A5_3'),
    out_A6_1    <- pred (Vpars = Pars_S_FLU, BW = 128,   Dose = 2.2, Dtimes = 1, route = "im")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A6_1'),
    out_A6_2    <- pred (Vpars = Pars_S_FLU, BW = 128,   Dose = 2.2, Dtimes = 1, route = "im")%>%
                   select(Time = Time, Conc = CP_Met)%>%mutate(Study = 'A6_2'),
    out_A6_3    <- pred (Vpars = Pars_S_FLU, BW = 128,   Dose = 2.2, Dtimes = 1, route = "im")%>%
                   select(Time = Time, Conc = CL)%>%mutate(Study = 'A6_3'),
    out_A6_4    <- pred (Vpars = Pars_S_FLU, BW = 128,   Dose = 2.2, Dtimes = 1, route = "im")%>%
                   select(Time = Time, Conc = CL_Met)%>%mutate(Study = 'A6_4'),
    out_A6_5    <- pred (Vpars = Pars_S_FLU, BW = 128,   Dose = 2.2, Dtimes = 1, route = "im")%>%
                   select(Time = Time, Conc = CK)%>%mutate(Study = 'A6_5'),
    out_A6_6    <- pred (Vpars = Pars_S_FLU, BW = 128,   Dose = 2.2, Dtimes = 1, route = "im")%>%
                   select(Time = Time, Conc = CK_Met)%>%mutate(Study = 'A6_6'),
    out_A7_1    <- pred (Vpars = Pars_S_FLU, BW = 3.4, Dose = 2.2,   Dtimes = 1, route = "im")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A7_1'),
    out_A7_2    <- pred (Vpars = Pars_S_FLU, BW = 3.4, Dose = 2.2, Dtimes = 1, route = "iv")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A7_2'),
    out_A8_1    <- pred (Vpars = Pars_C_FLU, BW = 127, Dose = 2.2, Dtimes = 3, route = "iv")%>%
                   select(Time = Time, Conc = CL)%>%mutate(Study = 'A8_1'),
    out_A8_2    <- pred (Vpars = Pars_C_FLU, BW = 127, Dose = 2.2, Dtimes = 3, route = "iv")%>%
                   select(Time = Time, Conc = CM)%>%mutate(Study = 'A8_2'),
    out_A8_3    <- pred (Vpars = Pars_C_FLU, BW = 127, Dose = 2.2, Dtimes = 3, route = "iv")%>%
                   select(Time = Time, Conc = CK)%>%mutate(Study = 'A8_3'),
    out_A8_4    <- pred (Vpars = Pars_C_FLU, BW = 127, Dose = 2.2, Dtimes = 3, route = "iv")%>%
                   select(Time = Time, Conc = CF)%>%mutate(Study = 'A8_4'),
    out_A9_1    <- pred (Vpars = Pars_C_FLU, BW = 288, Dose = 2.2, Dtimes = 1, route = "iv")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A9_1'),
    out_A9_2    <- pred (Vpars = Pars_C_FLU, BW = 288, Dose = 2.2, Dtimes = 1, route = "sc")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A9_2'),
    out_A10_1   <- pred (Vpars = Pars_C_FLU, BW = 462, Dose = 2.2, Dtimes = 1, route = "iv")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A10_1'),
    out_A10_2   <- pred (Vpars = Pars_C_FLU, BW = 462, Dose = 2.2, Dtimes = 1, route = "im")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A10_2'),
    out_A11_1   <- pred (Vpars = Pars_C_FLU, BW = 53.3, Dose = 2.2, Dtimes = 3, route = "iv")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A11_1'),
    out_A11_2   <- pred (Vpars = Pars_C_FLU, BW = 53.3, Dose = 2.2, Dtimes = 3, route = "iv")%>%
                   select(Time = Time, Conc = CL)%>%mutate(Study = 'A11_2'),
    out_A11_3   <- pred (Vpars = Pars_C_FLU, BW = 53.3, Dose = 2.2, Dtimes = 3, route = "iv")%>%
                   select(Time = Time, Conc = CM)%>%mutate(Study = 'A11_3'),
    out_A11_4   <- pred (Vpars = Pars_C_FLU, BW = 53.3, Dose = 2.2, Dtimes = 3, route = "iv")%>%
                   select(Time = Time, Conc = CK)%>%mutate(Study = 'A11_4'),
    out_A11_5   <- pred (Vpars = Pars_C_FLU, BW = 53.3, Dose = 2.2, Dtimes = 3, route = "iv")%>%
                   select(Time = Time, Conc = CL_Met)%>%mutate(Study = 'A11_5'),
    out_A12_1   <- pred (Vpars = Pars_C_FLU, BW = 570, Dose = 2.2, Dtimes = 4, route = "iv")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A12_1'),
    out_A13_1   <- pred (Vpars = Pars_C_FLU, BW = 60.2, Dose = 2.2, Dtimes = 1, route = "iv")%>%
                    select(Time = Time, Conc = CP)%>%mutate(Study = 'A13_1'),
    out_A14_1   <- pred (Vpars = Pars_C_FLU, BW = 562, Dose = 2.2, Dtimes = 1, route = "iv")%>%
                   select(Time = Time, Conc = CP)%>%mutate(Study = 'A14_1')
)
    
    
###/////////////////////////////////////////////////////////////////////////////
## Plot
Levels <- c('A1_1', 'A1_2', 'A2_1', 'A2_2','A2_3', 
            'A3_1','A3_2', 'A3_3','A3_4', 'A4_1',
            'A5_1', 'A5_2', 'A5_3', 'A6_1','A6_2', 
            'A6_3','A6_4', 'A6_5','A6_6', 'A7_1',
            'A7_2', 'A8_1', 'A8_2', 'A8_3','A8_4', 
            'A9_1','A9_2', 'A10_1','A10_2', 'A11_1',
            'A11_2', 'A11_3', 'A11_4', 'A11_5','A12_1', 
            'A13_1','A14_1')

Data_Cal$Study <- factor(Data_Cal$Study, levels=Levels)
Out_Cal$Study <- factor(Out_Cal$Study, levels=Levels)


p<-ggplot(data = Data_Cal, aes(x=Time, y=Conc)) + 
    geom_point(shape = 21, colour = "black", fill = "white", size = 1.5, stroke = 1.2) + 
    geom_line(data = data.frame(Out_Cal), aes(x = Time, y = Conc ), 
              linetype = 2, size = 0.6, colour = "black") +
    
    scale_y_log10 (labels = function(x) format(x, scientific = TRUE))  + 
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




ggsave("Fig S1.tiff",scale = 1.5,
       plot = p,
       width = 25, height = 15, units = "cm", dpi=320)

    
    











