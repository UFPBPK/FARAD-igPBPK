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


## Chemical-specific parameters for flunixin (Please refer to Table 3 for the full name of parameters)
Chme_pars_FLU <- c(
    fR          = 0.05,              
    fR1         = 0.01,
    Kunabs      = 0.5,                      
    Kabs        = 0.4,
    Kim         = 1.00,                            
    Ksc         = 0.40,                           
    KmetC       = 0.020,                            
    KehcC       = 0.05,     
    KurineC     = 0.10,              
    KurineC1    = 0.10,
    KbileC      = 0.50,             
    KbileC1     = 0.10, 
    Kdissim     = 1e-5,     
    Kdisssc     = 1e-5      
)   


## Loading tissue-composition data
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

PC_all_Cattle<-PCcoef_all(logP = 4.17, 
                          pKa  = 1.88, 
                          fup  = (1-0.95), 
                          BP   = 1, 
                          type = 2, 
                          dat  = dat_tiss_Cattle)

PC_all_Cattle<- do.call(rbind,PC_all_Cattle)

PC_all_Cattle <- PC_all_Cattle%>%as.data.frame()%>%mutate(Method = rownames(PC_all_Cattle)) %>% 
    filter(Method != "pcoeff_PKsim") %>%
    summarise(PF = mean(as.numeric(Kp_ad)), 
              PL = mean(as.numeric(Kp_li)), 
              PK = mean(as.numeric(Kp_ki)), 
              PM = mean(as.numeric(Kp_mu)),
              PRest = mean(as.numeric(Kp_sp, Kp_sk, Kp_lu, Kp_he, Kp_br, Kp_bo)))




PC_Cattle_FLU <- c(PL = PC_all_Cattle$PL,
                   PL1= PC_all_Cattle$PL,
                   PK = 4, #PC_all_Cattle$PK, 
                   PK1= 4, #PC_all_Cattle$PK, 
                   PM = 0.5,#PC_all_Cattle$PM,
                   PM1= PC_all_Cattle$PM, 
                   PF = PC_all_Cattle$PF,
                   PRest1 = PC_all_Cattle$PRest,
                   PRest = PC_all_Cattle$PRest)



## Defined the parameters
VarPars_Cattle <- c(Chme_pars_FLU, PC_Cattle_FLU)   ## Variable parameters 
FixPars_Cattle <- c(Phy_pars_Cattle, VarPars_Cattle) ## Fixed parameters




##---------------------------------------------------------------------
## Model fitting
data_flu <- read.csv(file = "Data_FLU.csv")
head(data_flu)


## Read the data set and later used in model calibration
## Study 8: FDA (1998) 
Obs_A8     <- data_flu %>% filter(Study ==8)
Obs_A8_L   <- Obs_A8 %>% filter(Matrix == "L")%>%select(Time = Time, CL = Conc.)
Obs_A8_M   <- Obs_A8 %>% filter(Matrix == "M")%>%select(Time = Time, CM = Conc.)
Obs_A8_K   <- Obs_A8 %>% filter(Matrix == "K")%>%select(Time = Time, CK = Conc.)
Obs_A8_F   <- Obs_A8 %>% filter(Matrix == "F")%>%select(Time = Time, CF = Conc.)

## Study 9: Shelver et al. (2013) 
Obs_A9     <- data_flu %>% filter(Study ==9)
Obs_A9_IV  <- Obs_A9 %>% filter(Route == "IV")%>%select(Time = Time, CP = Conc.)
Obs_A9_SC  <- Obs_A9 %>% filter(Route == "SC")%>%select(Time = Time, CP = Conc.)

## Study 10: Odensvik and Johansson (1995)
Obs_A10    <- data_flu %>% filter(Study ==10)
Obs_A10_IV <- Obs_A10 %>% filter(Route == "IV")%>%select(Time = Time, CP = Conc.)
Obs_A10_IM <- Obs_A10 %>% filter(Route == "IM")%>%select(Time = Time, CP = Conc.)

## Study 11: Kissell et al. (2016) for evaluation
Obs_A11 <- data_flu %>% filter(Study ==11)
Obs_A11_P <- Obs_A11 %>% filter(Matrix == "P")%>%select(Time = Time, CP = Conc.)
Obs_A11_L <- Obs_A11 %>% filter(Matrix == "L")%>%select(Time = Time, CL = Conc.)
Obs_A11_M <- Obs_A11 %>% filter(Matrix == "M")%>%select(Time = Time, CM = Conc.)
Obs_A11_K <- Obs_A11 %>% filter(Matrix == "K")%>%select(Time = Time, CK = Conc.)
Obs_A11_LMet <- Obs_A11 %>% filter(Matrix == "L-Met")%>%select(Time = Time, CL_Met = Conc.)

## Study 12: Jaroszewski et al. (2008)
Obs_A12 <- data_flu %>% filter(Study ==12 & Matrix == "P")%>%select(Time = Time, CP = Conc.)
Obs_A12_PMet <- data_flu %>% filter(Study ==12 & Matrix == "P-Met")%>%select(Time = Time, CP_Met = Conc.)

## Study 13: Kleinhenz et al. (2016)
Obs_A13 <- data_flu %>% filter(Study ==13)%>%select(Time = Time, CP = Conc.)

## Study 14: Odensvik (1995)
Obs_A14 <- data_flu %>% filter(Study ==14)%>%select(Time = Time, CP = Conc.)



## Prediction function
pred <- function (fixpars, Vpars, BW, tinterval = 24, Dose, Dtimes, route) {
    
    ## Get out of log domain
    parsinput <- exp(Vpars) 
    fixpars["BW"] <- BW
    Fracsc = ifelse(is.na(parsinput["Fracsc"]), 1, parsinput["Fracsc"]) 
    Fracim = ifelse(is.na(parsinput["Fracim"]), 1, parsinput["Fracim"]) 
    
    ## Make the dosing regiment
    BW          = BW
    tinterval   = tinterval
    TDOSE       = Dtimes
    Frac        = switch(route, "im"= Fracim, "sc" = Fracsc)
    MW          = 296.24
    MW1         = 312.24
    DOSEfast    = Dose*BW*Frac/MW              
    DOSEslow    = Dose*BW*(1-Frac)/MW  
    DOSE        = Dose*BW/MW  
    
    ## Create the exposure scenario
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

## Create the Cost function
Cost_Cattle <- function (pars, w) {
    
    # Prediction
    out_A8    <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 127, Dose = 2.2, Dtimes = 3, route = "iv")
    out_A9_IV <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 288, Dose = 2.2, Dtimes = 1, route = "iv")
    out_A9_SC <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 288, Dose = 2.2, Dtimes = 1, route = "sc")
    out_A10_IV<- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 462, Dose = 2.2, Dtimes = 1, route = "iv")
    out_A10_IM<- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 462, Dose = 2.2, Dtimes = 1, route = "im")
    out_A13   <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 60.2, Dose = 2.2, Dtimes = 1, route = "iv")
    out_A14   <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 562, Dose = 2.2, Dtimes = 1, route = "iv")
    out_A12   <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 570, Dose = 2.2, Dtimes = 4, route = "iv")
    
    
    ## Cost function
    cost<- modCost  (model = out_A8, obs = Obs_A8_L, x ="Time", weight = w)
    cost<- modCost  (model = out_A8, obs = Obs_A8_M, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A8, obs = Obs_A8_K, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A8, obs = Obs_A8_F, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A9_IV, obs = Obs_A9_IV, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A9_SC, obs = Obs_A9_SC, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A10_IV, obs = Obs_A10_IV, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A10_IM, obs = Obs_A10_IM, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A13, obs = Obs_A13, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A14, obs = Obs_A14, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A12, obs = Obs_A12, x ="Time",  cost = cost, weight = w)
    cost<- modCost  (model = out_A12, obs = Obs_A12_PMet, x ="Time", cost = cost, weight = w)
    
    return(cost)
}


##
Cost_Cattle(log(VarPars_Cattle), w="mean")


## Initial parameters
theta_init <- c(
    fR          = 0.1,             
    fR1         = 0.01,
    Kunabs      = 0.8,                      
    Kabs        = 0.4,
    Kim         = 0.5,
    Ksc         = 0.40,
    KmetC       = 0.02,
    KehcC       = 0.01,  
    KurineC     = 0.5,        
    KurineC1    = 0.05,    
    KbileC      = 0.5,        
    KbileC1     = 0.5, 
    PK          = 4,
    PK1         = 4,
    PF          = 0.6,
    PL          = 2,
    PL1         = 2,
    PRest1      = 8,
    PRest       = 8,
    PM          = 0.5,
    PM1         = 0.5)


## Sensitivity function (FROM FME package) 
## Check the sensitive parameters in the model
Sns_Cattle <- sensFun(func = Cost_Cattle, w = "mean",
                      parms = log(theta_init), varscale = 1)

Sen_1      <- summary(Sns_Cattle)
plot(summary(Sns_Cattle))


## Selected sensitive parameters; 
theta_Cattle <- theta_init[abs(Sen_1$Mean) >1.1*mean(abs(Sen_1$Mean))]
theta_Cattle

## Customized parameters
theta_Cattle <- c(
    fR          = 0.1,             
    fR1         = 0.01,
    Kunabs      = 0.8,                      
    #Kabs        = 0.4,
    #Kim         = 0.5,                          
    #Ksc         = 0.40,                         
    KmetC       = 0.02,                        
    KehcC       = 0.01,  
    KurineC     = 0.5,        
    KurineC1    = 0.05,    
    KbileC      = 0.5,        
    KbileC1     = 0.5, 
    PK          = 4,
    PK1         = 4,
    PF          = 0.6,
    PL          = 2,
    #PL1         = 2,
    PRest1      = 8,
    PRest       = 8,
    PM          = 0.5)
    #PM1         = 0.5)

# ## Selected sensitive parameters; 
Fit_Cattle_FLU <- modFit(f       = Cost_Cattle, 
                         p       = log(theta_Cattle), 
                         w       = "mean",
                         method  = "Marq",
                         control = nls.lm.control(nprint = 1))

summary(Fit_Cattle_FLU) 


exp(Fit_Cattle_FLU$par)


## Create the plot
Ccattle_FLU<-Cost_Cattle(Fit_Cattle_FLU$par, w="mean")

a<-data.frame(Ccattle_FLU$residuals)

# Make a plot
PDat_Cattle <- cbind.data.frame (OBS = Ccattle_FLU$residuals$obs,
                          PRE = Ccattle_FLU$residuals$mod,
                          RES = Ccattle_FLU$residuals$res)

PDat_Cattle <- PDat_Cattle %>% mutate (Log.OBS = log(OBS,10), Log.PRE = log(PRE,10), Species = "Cattle")

fit_Cattle <- lm(Log.OBS ~ Log.PRE, data = PDat_Cattle)
summary(fit_Cattle)

PlotDat_Cattle <- PDat_Cattle %>% mutate(prediction = predict(fit_Cattle), OPR = PRE/OBS)


p1 <-
    ggplot(PlotDat_Cattle, aes(Log.PRE, Log.OBS)) +
    geom_point  (aes(shape   = as.factor(Species)), colour = "steelblue4", size = 4)  +
    geom_abline (intercept = 0,
                 slope     = 1,
                 color     ="steelblue4", size = 1, alpha = 0.8) +
    annotation_logticks() +
    scale_y_continuous(limits = c(-4,3),labels = scales::math_format(10^.x))+
    scale_x_continuous(limits = c(-4,3),labels = scales::math_format(10^.x))
p1

## Save the fitting results
saveRDS(Fit_Cattle_FLU, file = 'Fit_Cattle_FLU.rds')
saveRDS(Ccattle_FLU, file = 'Ccattle_FLU.rds')


##//////////////////////////////////////////////////////////////////////////////
## Model evaluation
## Study 11: Kissell et al. (2016)

### Make cost function
Cost_Cattle_eva <- function (pars, w) {
    
    # Prediction
    out_A11_1 <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 53.3, Dose = 2.2, Dtimes = 3, route = "iv")
    out_A11_2 <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 53.3, Dose = 2.2, Dtimes = 3, route = "iv")
    
    cost<- modCost  (model = out_A11_1, obs = Obs_A11_P, x ="Time", weight = w)
    cost<- modCost  (model = out_A11_2, obs = Obs_A11_L, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A11_2, obs = Obs_A11_M, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A11_2, obs = Obs_A11_K, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A11_2, obs = Obs_A11_LMet, x ="Time", cost = cost, weight = w)

    
    return(cost)
}


Cost_Cattle_eva (Fit_Cattle_FLU$par, w="mean")

## Make the data for the plot
Cattle_FLU_eva <-Cost_Cattle_eva (Fit_Cattle_FLU$par, w="mean")


## Make a plot
PDat_eva   <- cbind.data.frame (OBS = Cattle_FLU_eva$residuals$obs,
                            PRE = Cattle_FLU_eva$residuals$mod,
                            RES = Cattle_FLU_eva$residuals$res)

PDat_eva <- PDat_eva %>% mutate (Log.OBS = log(OBS, 10), Log.PRE = log(PRE, 10), Species = "Cattle")


## Estimate the r-squared values
fit_eva <- lm(Log.OBS ~ Log.PRE, data = PDat_eva)
summary(fit_eva)

## Plot the results of model evaluation
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

## Save the fitting results
saveRDS(Cattle_FLU_eva, file = 'Cattle_FLU_eva.rds')






