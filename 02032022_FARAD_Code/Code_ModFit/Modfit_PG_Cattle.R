## Loading required R package
library(FME)
library(mrgsolve)
library(dplyr)
library(minpack.lm)  
library(ggplot2)
library(tidyr)

# Input the PBPK model
source (file = "GenPBPK.R") # Loading the generic PBPK model code
source (file = "GetPC.R")   # Loading the tissue composition based tissue:plasma model

## Load Model
mod <- mcode_cache("pbpk", GenricPBPK) 

## Input the Cattle parameters
## Physiological parameters (Please refer to Table 2 for the full name of parameters)

Phy_pars_Cattle <- c(
    BW          = 300,      
    Htc         = 0.378,     
    QCC         = 5.45,    
    QLCa        = 0.44,	                      
    QKCa        = 0.11,	                   
    QMCa        = 0.28,	                     
    QFCa        = 0.08,
    QRestCa     = (1-0.44-0.11-0.28-0.08),
    VLCa        = 0.0122,                                  
    VKCa        = 0.0021,                             
    VFCa        = 0.1218,                                   
    VMCa        = 0.361,                                 
    VbloodCa    = 0.0399,
    VRestCa     = (1-0.0122-0.0021-0.1218-0.361-0.0399)
)

## Chemical-specific parameters for florfenicol (Please refer to Table 3 for the full name of parameters)
Chme_pars_PG <- c(
    fR          = 0.5,              
    fR1         = 0.5,
    Fracsc       = 0.5,
    Fracim       = 0.65,
    Kunabs      = 0.01,                     
    Kim         = 0.07,                             
    Ksc         = 0.25,                         
    Kabs        = 1.91,
    Kdissim     = 0.007,     
    Kdisssc     = 0.005,      
    KehcC        = 0.05,
    KmetC       = 0.05,                            
    KurineC     = 1.4,           
    KbileC      = 0.05           
)   



## Estiamte the partition coefficient
dat_tiss_Cattle <- read.csv(file = "Tiss_Cattle.csv")

## Chemical characterization of FLO
## Basic information; logP, pka collected from: https://pubchem.ncbi.nlm.nih.gov/compound/Penicillin-g#section=Solubility
## Abbreviations: 
## ogP: the log-scale n-octanol:buffer partition coefficient;--- 
## pka: acid disscoiation constant; 
## BP is the in vitro blood:plasma ratio;--
## fup: fraction unbound in plasma;---
## dat: input the tissue composition data;---
## type: chemical type (#1-neutral; #2-monoprotic acid; #3-monoprotic base; #4-diprotic acid; #5-diprotic base--
## #6-monoprotic acid monoprotic base (acid comes first); #7-triprotic acid; #8-triprotic base; ---
## #9-diprotic acid monoprotic base (first two are acid);  #10-diprotic base monoprotic acid (first one is acid))

PC_all_Cattle<-PCcoef_all(logP = 1.83, 
                          pKa  = 2.74, 
                          fup  = (1-0.366), 
                          BP   = 1, 
                          type = 2, #3-monoprotic base
                          dat  = dat_tiss_Cattle)


PC_all_Cattle<- do.call(rbind,PC_all_Cattle)

PC_all_Cattle <- PC_all_Cattle%>%as.data.frame()%>%mutate(Method = rownames(PC_all_Cattle)) %>% 
    filter(Method != "pcoeff_PKsim") %>%
    summarise(PF = mean(as.numeric(Kp_ad)), 
              PL = mean(as.numeric(Kp_li)), 
              PK = mean(as.numeric(Kp_ki)), 
              PM = mean(as.numeric(Kp_mu)),
              PRest = mean(as.numeric(Kp_sp, Kp_sk, Kp_lu, Kp_he, Kp_br, Kp_bo)))



PC_Cattle_PG <- c(PL = PC_all_Cattle$PL,
                  PL1= PC_all_Cattle$PL,
                  PK = PC_all_Cattle$PK, 
                  PK1= PC_all_Cattle$PK, 
                  PM = PC_all_Cattle$PM,
                  PM1= PC_all_Cattle$PM, 
                  PF = 0.04,#PC_all_Cattle$PF,
                  PRest1 = PC_all_Cattle$PRest,
                  PRest = PC_all_Cattle$PRest)


## Defined the parameters
VarPars_Cattle <- c(Chme_pars_PG, PC_Cattle_PG)   ## Variable parameters 
FixPars_Cattle <- c(Phy_pars_Cattle, VarPars_Cattle) ## Fixed parameters

##
##////////////////////////////////////////////////////////////////////////////////
# Input dataset
## Model fitting
Data_PG <- read.csv(file = "Data_PG.csv")
head(Data_PG)

## Study 7; Papich et al., 1993; IM, SC; P; Dose: 20, 65, 65 mg/kg/day
Obs_A7      <- Data_PG %>% filter(Study ==7)
Obs_A7_IM_1 <- Obs_A7 %>% filter(Route == 'IM' & Dose == "23.7" & Repeat == "5")%>%select(Time = Time, CP = Conc)
Obs_A7_IM_2 <- Obs_A7 %>% filter(Route == 'IM' & Dose == "65.4" & Repeat == "5")%>%select(Time = Time, CP = Conc)
Obs_A7_IM_3 <- Obs_A7 %>% filter(Route == 'IM' & Dose == "65.4" & Repeat == "1")%>%select(Time = Time, CP = Conc)
Obs_A7_SC   <- Obs_A7 %>% filter(Route == 'SC'& Dose == "65.4")%>%select(Time = Time, CP = Conc)

## Study 8;Korsrud et al., 1993;  IM; P, L, K, M; Dose: 24, 65 mg/kg/day
Obs_A8 <- Data_PG %>% filter(Study ==8)
Obs_A8_IM_1_P <- Obs_A8 %>% filter(Route == 'IM' & Dose == "23.7" & Matrix == "P")%>%select(Time = Time, CP = Conc)
Obs_A8_IM_1_L <- Obs_A8 %>% filter(Route == 'IM' & Dose == "23.7" & Matrix == "L")%>%select(Time = Time, CL = Conc)
Obs_A8_IM_1_K <- Obs_A8 %>% filter(Route == 'IM' & Dose == "23.7" & Matrix == "K")%>%select(Time = Time, CK = Conc)
Obs_A8_IM_1_M <- Obs_A8 %>% filter(Route == 'IM' & Dose == "23.7" & Matrix == "M")%>%select(Time = Time, CM = Conc)
Obs_A8_IM_2_P <- Obs_A8 %>% filter(Route == 'IM' & Dose == "65.4" & Matrix == "P")%>%select(Time = Time, CP = Conc)
Obs_A8_IM_2_L <- Obs_A8 %>% filter(Route == 'IM' & Dose == "65.4" & Matrix == "L")%>%select(Time = Time, CL = Conc)
Obs_A8_IM_2_K <- Obs_A8 %>% filter(Route == 'IM' & Dose == "65.4" & Matrix == "K")%>%select(Time = Time, CK = Conc)
Obs_A8_IM_2_M <- Obs_A8 %>% filter(Route == 'IM' & Dose == "65.4" & Matrix == "M")%>%select(Time = Time, CM = Conc)

## Study 9; Trollder et al. 1986;  SC; Dose: 7.4 mg/kg/day; P
Obs_A9   <- Data_PG %>% filter(Study ==9)
Obs_A9_SC <- Obs_A9 %>% select(Time = Time, CP = Conc)

## Study 10; Chiesa et al., 2006; IM; Dose: 7 mg/kg/day; K,P; Repeat: 3;
Obs_A10 <- Data_PG %>% filter(Study ==10)
Obs_A10_P <- Obs_A10 %>% filter(Matrix == 'P')%>%select(Time = Time, CP = Conc)
Obs_A10_K <- Obs_A10 %>% filter(Matrix == 'K')%>%select(Time = Time, CK = Conc)

## Study 11; Chiesa et al., 2006; IM; Dose: 7 mg/kg/day; K,P; Repeat: 3;
Obs_A11 <- Data_PG %>% filter(Study ==11)
Obs_A11_P <- Obs_A11 %>% filter(Matrix == 'P')%>%select(Time = Time, CP = Conc)



### Define the prediction function
pred <- function (fixpars, Vpars, BW, tinterval = 24, Dose, Dtimes, route) {
    
    ## Get out of log domain
    parsinput <- exp(Vpars) 
    fixpars["BW"] <- BW
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
    out <- mod %>% param (fixpars) %>% param (parsinput) %>% 
        update(atol = 1E-10, rtol = 1E-6, maxsteps = 50000) %>%
        mrgsim_d (data = ex, tgrid = tsamp)
    
    
    outdf = cbind.data.frame( Time  = out$time,
                              CP    = out$Plasma*MW,
                              CL    = out$Liver*MW,
                              CK    = out$Kidney*MW,
                              CM    = out$Muscle*MW,
                              CF    = out$Fat*MW)
    return (outdf)
    
}

## Define cost function

Cost_Cattle <- function (pars, w) {
    
    # Prediction
    out_A7_IM_1  <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 480, Dose = 23.7, Dtimes = 5, route = "im")
    out_A7_IM_2  <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 480, Dose = 65.4, Dtimes = 5, route = "im")
    out_A7_IM_3  <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 480, Dose = 65.4, Dtimes = 1, route = "im")
    out_A7_SC    <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 480, Dose = 65.4, Dtimes = 1, route = "sc")
    out_A8_IM_1  <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 485, Dose = 23.7, Dtimes = 5, route = "im")
    out_A8_IM_2  <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 485, Dose = 65.4, Dtimes = 5, route = "im")
    out_A9_SC    <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 93, Dose = 8.9, Dtimes = 1, route = "sc")
    #out_A10_IM    <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 262, Dose = 6.9, Dtimes = 3, tinterval = 12, route = "im")
    
    
    cost<- modCost  (model = out_A7_IM_1, obs = Obs_A7_IM_1, x ="Time", weight = w)
    cost<- modCost  (model = out_A7_IM_2, obs = Obs_A7_IM_2, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A7_IM_3, obs = Obs_A7_IM_3, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A7_SC, obs = Obs_A7_SC, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A8_IM_1, obs = Obs_A8_IM_1_P, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A8_IM_1, obs = Obs_A8_IM_1_L, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A8_IM_1, obs = Obs_A8_IM_1_K, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A8_IM_1, obs = Obs_A8_IM_1_M, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A8_IM_2, obs = Obs_A8_IM_2_P, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A8_IM_2, obs = Obs_A8_IM_2_L, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A8_IM_2, obs = Obs_A8_IM_2_K, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A8_IM_2, obs = Obs_A8_IM_2_M, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A9_SC, obs = Obs_A9_SC, x ="Time", cost = cost, weight = w)
    # cost<- modCost  (model = out_A10_IM, obs = Obs_A10_P, x ="Time", cost = cost, weight = w)
    # cost<- modCost  (model = out_A10_IM, obs = Obs_A10_K, x ="Time", cost = cost, weight = w)
    # 
    return(cost)
}

Cost_Cattle(log(VarPars_Cattle), w="mean")

## Set up the initial parameter list

theta_init <- c(
    fR          = 0.62,                  
    Fracsc      = 0.70,
    Fracim      = 0.60,
    Kunabs      = 0.5,                    
    Kim         = 0.1,                             
    Ksc         = 0.25,                             
    Kabs        = 1.91,
    KehcC       = 0.001,
    KmetC       = 0.2,                            
    Kdissim     = 0.007,     
    Kdisssc     = 0.005,    
    KurineC     = 0.8,            
    KbileC      = 0.6,          
    PL          = 1,
    PM          = 0.3,
    PF          = 0.04,
    PK          = 2,
    PRest        = 7.49
)

## Sensitivity function (FME) 
## Check the sensitive parameters in the model
Sns_Cattle <- sensFun(func = Cost_Cattle, w = "mean",
                      parms = log(theta_init), varscale = 1)

Sen_1 <- summary(Sns_Cattle)
plot(summary(Sns_Cattle))

## Selected senstive parameters; 
theta_Cattle <- theta_init[abs(Sen_1$Min) >1*mean(abs(Sen_1$Min))]
theta_Cattle

##
theta_Cattle <- c(
    fR          = 0.62,                  
    Fracsc      = 0.70,
    Fracim      = 0.60,
    Kunabs      = 0.5,                    
    Kim         = 0.1,                             
    Ksc         = 0.25,                             
    #Kabs        = 1.91,
    KehcC       = 0.001,
    KmetC       = 0.2,                            
    Kdissim     = 0.007,     
    #Kdisssc     = 0.005,    
    KurineC     = 0.8,            
    KbileC      = 0.6,          
    PL          = 1,
    PM          = 0.3,
    #PF          = 0.04,
    PK          = 2,
    PRest        = 7.49
)


# ## Selected sensitive parameters; 

Fit_Cattle_PG <- modFit(f = Cost_Cattle, p = log(theta_Cattle), w = "mean",method ="Marq",
                        control = nls.lm.control(nprint = 1))
summary(Fit_Cattle_PG) 


exp(Fit_Cattle_PG$par)

###
Ccattle_PG<-Cost_Cattle(Fit_Cattle_PG$par, w="mean")

## Make plot
PDat_Cattle <- cbind.data.frame (OBS = Ccattle_PG$residuals$obs,
                                 PRE = Ccattle_PG$residuals$mod,
                                 RES = Ccattle_PG$residuals$res)

PDat_Cattle <- PDat_Cattle %>% mutate (Log.OBS = log(OBS,10), Log.PRE = log(PRE,10), Species = "Cattle")

fit_Cattle <- lm(Log.PRE ~ Log.OBS, data = PDat_Cattle)
summary(fit_Cattle)

PlotDat_Cattle <- PDat_Cattle %>% mutate(prediction = predict(fit_Cattle), OPR = PRE/OBS)


p1 <-
    ggplot(PlotDat_Cattle, aes(Log.PRE, Log.OBS)) +
    geom_point  (aes(shape   = as.factor(Species)), colour = "steelblue4", size = 4)  +
    geom_abline (intercept = 0,
                 slope     = 1,
                 color     ="steelblue4", size = 1, alpha = 0.8) +
    annotation_logticks() +
    scale_y_continuous(limits = c(-3,1), labels = scales::math_format(10^.x))+
    scale_x_continuous(limits = c(-3,1),labels = scales::math_format(10^.x))
p1

## Save the fitting results
saveRDS(Fit_Cattle_PG, file = 'Fit_Cattle_PG.rds')
saveRDS(Ccattle_PG, file = 'Ccattle_PG.rds')



#///////////////////////////////////////////////////////////////////////////////
## Model evaluation
### Make cost function
Cost_Cattle_eva <- function (pars, w) {
    
    out_A11   <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 624, Dose = 20.8, Dtimes = 1, route = "im")
    out_A10_IM    <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 262, Dose = 6.9, Dtimes = 3, tinterval = 12, route = "im")
    
    # Prediction
    cost<- modCost  (model = out_A11, obs = Obs_A11_P, x ="Time", weight = w)
    cost<- modCost  (model = out_A10_IM, obs = Obs_A10_P, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A10_IM, obs = Obs_A10_K, x ="Time", cost = cost, weight = w)
    
    return(cost)
}


Cost_Cattle_eva (Fit_Cattle_PG$par, w="mean")


## Make the data for the plot
Cattle_PG_eva  <-Cost_Cattle_eva(Fit_Cattle_PG$par, w="mean")

## Make a plot
PDat_eva   <- cbind.data.frame (OBS = Cattle_PG_eva$residuals$obs,
                                PRE = Cattle_PG_eva$residuals$mod,
                                RES = Cattle_PG_eva$residuals$res)

PDat_eva <- PDat_eva %>% mutate (Log.OBS = log(OBS, 10), Log.PRE = log(PRE, 10), Species = "Swine")


## Estimate the r-squared values
fit_eva <- lm(Log.OBS ~ Log.PRE, data = PDat_eva)
summary(fit_eva)

##
PlotDat_eva <- PDat_eva %>% mutate(prediction = predict(fit_eva), OPR = PRE/OBS)


p2 <-
    ggplot(PlotDat_eva, aes(Log.PRE, Log.OBS)) +
    geom_point  (aes(shape   = as.factor(Species)), colour = "steelblue4", size = 4)  +
    geom_abline (intercept = 0,
                 slope     = 1,
                 color     ="steelblue4", size = 1, alpha = 0.8) +
    annotation_logticks() +
    scale_y_continuous(limits = c(-4,2), labels = scales::math_format(10^.x))+
    scale_x_continuous(limits = c(-4,2),labels = scales::math_format(10^.x))

p2

# Save the fitting results
saveRDS(Cattle_PG_eva, file = 'Cattle_PG_eva.rds')







