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
Chme_pars_FLO <- c(
    fR          = 0.80,              
    fR1         = 0.71,
    Fracsc      = 0.60,
    Fracim      = 0.45,
    Kabs        = 0.4,
    KehcC       = 0.05,
    Kunabs      = 5e-5,                       
    Kim         = 0.1618,                             
    Ksc         = 0.128,                            
    KmetC       = 0.025,                             
    KurineC     = 0.05102,       
    KurineC1    = 4.5675E-2	,  
    KbileC      = 1E-2,            
    KbileC1     = 1E-2,     
    Kdissim     = 0.0104,     
    Kdisssc     = 0.0054      
)   

## Predict the partition coefficient (PC)
dat_tiss_Cattle <- read.csv(file = "Tiss_Cattle.csv")

## Chemical characterization of FLO
## Basic information; logP, pka collected from: https://www.agilent.com/cs/library/applications/5990-3615EN.pdf 
## Abbreviations: 
## ogP: the log-scale n-octanol:buffer partition coefficient;--- 
## pka: acid disscoiation constant; 
## BP is the in vitro blood:plasma ratio;--
## fup: fraction unbound in plasma;---
## dat: input the tissue composition data;---
## type: chemical type (#1-neutral; #2-monoprotic acid; #3-monoprotic base; #4-diprotic acid; #5-diprotic base--
## #6-monoprotic acid monoprotic base (acid comes first); #7-triprotic acid; #8-triprotic base; ---
## #9-diprotic acid monoprotic base (first two are acid);  #10-diprotic base monoprotic acid (first one is acid))

PC_all_Cattle_1<-PCcoef_all(logP = 0.39, 
                          pKa  = 9.87, 
                          fup  = (1-0.19), 
                          BP   = 1, 
                          type = 3, #3-monoprotic base
                          dat  = dat_tiss_Cattle)


PC_all_Cattle_2<-PCcoef_all(logP = 0.73, 
                          pKa  = 10.6, 
                          fup  = (1-0.29), 
                          BP   = 1, 
                          type = 3, #3-monoprotic base
                          dat  = dat_tiss_Cattle)


PC_all_Cattle_1<- do.call(rbind,PC_all_Cattle_1)

PC_all_Cattle_1 <- PC_all_Cattle_1%>%as.data.frame()%>%mutate(Method = rownames(PC_all_Cattle_1)) %>% 
    filter(Method != "pcoeff_PKsim") %>%
    summarise(PF = mean(as.numeric(Kp_ad)), 
              PL = mean(as.numeric(Kp_li)), 
              PK = mean(as.numeric(Kp_ki)), 
              PM = mean(as.numeric(Kp_mu)),
              PRest = mean(as.numeric(Kp_sp, Kp_sk, Kp_lu, Kp_he, Kp_br, Kp_bo)))



PC_Cattle_FLO <- c(PL = PC_all_Cattle_1$PL,
                   PL1= 8.74,
                   PK = PC_all_Cattle_1$PK, 
                   PK1= 1.3, 
                   PM = PC_all_Cattle_1$PM,
                   PM1= 0.27, 
                   PF = PC_all_Cattle_1$PF,
                   PRest1 = 0.1332, # yang et al. (2019)
                   PRest  = 0.1332)

## Defined the parameters
VarPars_Cattle <- c(Chme_pars_FLO, PC_Cattle_FLO)   ## Variable parameters 
FixPars_Cattle <- c(Phy_pars_Cattle, VarPars_Cattle) ## Fixed parameters


##---------------------------------------------------------------------
## Model fitting
Data_FLO <- read.csv(file = "Data_FLO.csv")
head(Data_FLO)

## Read the dataset and later used in model calibration
## Study 11; Varma et al., 1986; IV, PO, Matrix = serum/ plasma; Dose = 22 mg/kg/day
Obs_A11 <- Data_FLO %>% filter(Study ==11)
Obs_A11_IV <- Obs_A11 %>% filter(Route == "IV")%>%select(Time = Time, CP = Conc)
Obs_A11_PO <- Obs_A11 %>% filter(Route == "PO")%>%select(Time = Time, CP = Conc)

## Study 12: Sidhu et al., 2014
Obs_A12 <- Data_FLO %>% filter(Study ==12)
Obs_A12_SC <- Obs_A12 %>%select(Time = Time, CP = Conc)

## Study 13: Lobell et al., 1994
Obs_A13    <- Data_FLO %>% filter(Study ==13)
Obs_A13_IV <- Obs_A13 %>% filter(Route == "IV")%>%select(Time = Time, CP = Conc)
Obs_A13_IM <- Obs_A13 %>% filter(Route == "IM")%>%select(Time = Time, CP = Conc)

## Study 14: Lacroix et al., 2011
Obs_A14    <- Data_FLO %>% filter(Study ==14)
Obs_A14_IM <- Obs_A14 %>% filter(Route == "IM")%>%select(Time = Time, CP = Conc)
Obs_A14_SC <- Obs_A14 %>% filter(Route == "SC")%>%select(Time = Time, CP = Conc)

## Study 15: Kawalek et al., 2016
Obs_A15    <- Data_FLO %>% filter(Study ==15)
Obs_A15_SC <- Obs_A15 %>% filter(Route == "SC")%>%select(Time = Time, CP = Conc)

## Study 16: Gilliam et al., 2008
Obs_A16    <- Data_FLO %>% filter(Study ==16)
Obs_A16_IV <- Obs_A16 %>% filter(Route == "IV")%>%select(Time = Time, CP = Conc)

## Study 17: de Craene et al., 1997
Obs_A17    <- Data_FLO %>% filter(Study ==17)
Obs_A17_IV <- Obs_A17 %>% filter(Route == "IV")%>%select(Time = Time, CP = Conc)

## Study 18: Croubels et al., 2006
Obs_A18    <- Data_FLO %>% filter(Study ==18)
Obs_A18_PO <- Obs_A18 %>% filter(Route == "PO")%>%select(Time = Time, CP = Conc)

## Study 19: Bretzlaff et al., 1987
Obs_A19    <- Data_FLO %>% filter(Study ==19)
Obs_A19_IV <- Obs_A19 %>% filter(Route == "IV")%>%select(Time = Time, CP = Conc)

## Study 20; Schering-Plough Animal Health Corporation, 2008
Obs_A20   <- Data_FLO %>% filter(Study == 20)
Obs_A20_L <- Obs_A20 %>% filter(Matrix == "Liver") %>%select(Time = Time, CL_Met = Conc)
Obs_A20_K <- Obs_A20 %>% filter(Matrix == "Kidney")%>%select(Time = Time, CK_Met = Conc)
Obs_A20_M <- Obs_A20 %>% filter(Matrix == "Muscle")%>%select(Time = Time, CM_Met = Conc)

## Study 21; Laboratories, Ltd., 2015
Obs_A21   <- Data_FLO %>% filter(Study == 21)
Obs_A21_L <- Obs_A21 %>% filter(Matrix == "Liver") %>%select(Time = Time, CL_Met = Conc)
Obs_A21_M <- Obs_A21 %>% filter(Matrix == "Muscle")%>%select(Time = Time, CM_Met = Conc)

## Study 22; Intervet, Inc., 2009
Obs_A22   <- Data_FLO %>% filter(Study == 22)
Obs_A22_L <- Obs_A22 %>% filter(Matrix == "Liver")%>%select(Time = Time, CL_Met = Conc)


##
###
pred <- function (fixpars, Vpars, BW, tinterval = 24, Dose, Dtimes, route) {
    
    ## Get out of log domain
    parsinput <- exp(Vpars) 
    fixpars["BW"] <- BW
    Fracsc = ifelse(is.na(parsinput["Fracsc"]), 0.60, parsinput["Fracsc"]) 
    Fracim = ifelse(is.na(parsinput["Fracim"]), 0.45, parsinput["Fracim"]) 
    
    ## Exposure scenarios
    BW          = BW
    tinterval   = tinterval
    TDOSE       = Dtimes
    Frac        = switch(route, "im"= Fracim, "sc" = Fracsc)
    MW          = 358.21 
    MW1         = 247.28
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
## If you run the cost function, you will find the Warning message listed below:
## In regularize.values(x, y, ties, missing(ties)) :collapsing to unique 'x' values
## It casued by the function approx in later R 3.6.0 version; But it doesn't influence the results

Cost_Cattle <- function (pars, w) {
    
    # Prediction
    out_A11_IV  <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 64, Dose = 22, Dtimes = 1, route = "iv")
    out_A11_PO  <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 64, Dose = 22, Dtimes = 1, route = "oral")
    out_A12_SC  <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 204, Dose = 40, Dtimes = 1, route = "sc")
    out_A13_IV  <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 113, Dose = 20, Dtimes = 1, route = "iv")
    out_A13_IM  <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 113, Dose = 20, Dtimes = 1, route = "im")
    out_A14_IM  <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 563, Dose = 40, Dtimes = 1, route = "im")
    out_A14_SC  <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 563, Dose = 40, Dtimes = 1, route = "sc")
    #out_A15_SC  <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 700, Dose = 40, Dtimes = 1, route = "sc")
    out_A16_IV  <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 688, Dose = 2.2, Dtimes = 1, route = "iv")
    out_A17_IV  <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 290, Dose = 20, Dtimes = 1, route = "iv")
    out_A18_PO  <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 79.5, Dose = 20, Dtimes = 1, route = "oral")
    out_A19_IV  <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 550, Dose = 50, Dtimes = 1, route = "iv")
    out_A20_SC  <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 228, Dose = 40, Dtimes = 1, route = "sc")
    out_A21_SC  <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 302, Dose = 40, Dtimes = 1, route = "sc")
    out_A22_SC  <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 292, Dose = 40, Dtimes = 1, route = "sc")
    
    
    cost<- modCost  (model = out_A11_IV, obs = Obs_A11_IV, x ="Time", weight = w)
    cost<- modCost  (model = out_A11_PO, obs = Obs_A11_PO, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A12_SC, obs = Obs_A12_SC, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A13_IV, obs = Obs_A13_IV, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A13_IM, obs = Obs_A13_IM, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A14_IM, obs = Obs_A14_IM, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A14_SC, obs = Obs_A14_SC, x ="Time", cost = cost, weight = w)
    #cost<- modCost  (model = out_A15_SC, obs = Obs_A15_SC, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A16_IV, obs = Obs_A16_IV, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A17_IV, obs = Obs_A17_IV, x ="Time", cost = cost, weight = w)
    #cost<- modCost  (model = out_A18_PO, obs = Obs_A18_PO, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A19_IV, obs = Obs_A19_IV, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A20_SC, obs = Obs_A20_L, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A20_SC, obs = Obs_A20_K, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A20_SC, obs = Obs_A20_M, x ="Time", cost = cost, weight = w)
    #cost<- modCost  (model = out_A21_SC, obs = Obs_A21_L, x ="Time", cost = cost, weight = w)
    #cost<- modCost  (model = out_A21_SC, obs = Obs_A21_M, x ="Time", cost = cost, weight = w)
    #cost<- modCost  (model = out_A22_SC, obs = Obs_A22_L, x ="Time", cost = cost, weight = w)
    
    return(cost)
}

Cost_Cattle(log(VarPars_Cattle), w="mean")

## Defined the initial parameter list

theta_init <- c(
    Fracsc    = 0.26,
    Fracim    = 0.45,
    fR        = 0.80,              
    fR1       = 0.71,
    Kunabs    = 0.01,
    Kabs      = 1.4,
    KehcC     = 0.15,
    KmetC     = 0.25,
    Kmet1C    = 0.01,
    KurineC   = 0.3,
    KurineC1  = 1E-3,
    KbileC    = 1E-1,
    KbileC1   = 1E-1,
    Kim       = 0.1618,                             
    Ksc       = 0.128,
    Kdissim   = 0.0104,
    Kdisssc   = 0.0054,
    PL        = 2.2,
    PL1       = 8.74,
    PM        = 1.33,
    PM1       = 1.33,
    PF        = 0.615,
    PK        = 0.9,
    PK1       = 0.9,
    PRest     = 0.133,
    PRest1    = 0.121
)

## Sensitivity function (FME) 
## Check the senstive parameters in the model
Sns_Cattle <- sensFun(func = Cost_Cattle, w = "mean",
                      parms = log(theta_init), varscale = 1)

Sen_1 <- summary(Sns_Cattle)
plot(summary(Sns_Cattle))



## Selected senstive parameters; 
theta <- theta_init[abs(Sen_1$Min) >1.2*mean(abs(Sen_1$Min))]
theta

## customize parameter
theta <- c(
           Fracsc     = 0.66,
           Fracim     = 0.45,
           fR         = 0.80,              
           #fR1       = 0.71,
           Kunabs     = 0.01,
           #Kabs      = 1.4,
           KehcC     = 0.15,
           KmetC     = 0.25,
           #Kmet1C    = 0.01,
           KurineC   = 0.3,
           KurineC1  = 1E-3,
           KbileC    = 1E-1,
           KbileC1   = 1E-1,
           #Kim       = 0.1618,                             
           #Ksc       = 0.128,
           #Kdissim   = 0.0104,
           Kdisssc   = 0.0054,
           PL        = 2.2,
           PL1       = 8.74,
           #PM        = 1.33,
           #PM1       = 1.33,
           PF        = 0.615,
           PK        = 0.9,
           #PK1       = 0.9,
           PRest     = 0.133,
           PRest1    = 0.121)

# ## Selected sensitive parameters; 

Fit_Cattle_FLO <- modFit(f = Cost_Cattle, p = log(theta), w = "mean",method ="Marq",
              control = nls.lm.control(nprint = 1))

summary(Fit_Cattle_FLO) 

exp(Fit_Cattle_FLO$par)

###
Ccattle_FLO<-Cost_Cattle(Fit_Cattle_FLO$par, w="mean")

PDat <- cbind.data.frame (OBS = Ccattle_FLO$residuals$obs,
                          PRE = Ccattle_FLO$residuals$mod,
                          RES = Ccattle_FLO$residuals$res)

PDat <- PDat %>% mutate (Log.OBS = log(OBS, 10), Log.PRE = log(PRE, 10), Species = "Cattle")

fit <- lm(Log.OBS ~ Log.PRE, data = PDat)
summary(fit)

PlotDat <- PDat %>% mutate(prediction = predict(fit), OPR = PRE/OBS)


p1 <-
    ggplot(PlotDat, aes(Log.PRE, Log.OBS)) +
    geom_point  (aes(shape   = as.factor(Species)), colour = "steelblue4", size = 4)  +
    geom_abline (intercept = 0,
                 slope     = 1,
                 color     ="steelblue4", size = 1, alpha = 0.8) +
    annotation_logticks() +
    scale_y_continuous(limits = c(-2,2), labels = scales::math_format(10^.x))+
    scale_x_continuous(limits = c(-2,2),labels = scales::math_format(10^.x))
p1

# Save the fitting results
saveRDS(Fit_Cattle_FLO, file = 'Fit_Cattle_FLO.rds')
saveRDS(Ccattle_FLO, file = 'Ccattle_FLO.rds')

##//////////////////////////////////////////////////////////////////////////////
## Model evaluation
## Study 11: Jaroszewski et al. (2008)

### Make cost function
Cost_Cattle_eva <- function (pars, w) {
    
    # Prediction
    out_A18_PO  <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 79.5, Dose = 20, Dtimes = 1, route = "oral")
    out_A22_SC  <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 292, Dose = 40, Dtimes = 1, route = "sc")
    out_A21_SC  <- pred (fixpars = FixPars_Cattle, Vpars = pars, BW = 302, Dose = 40, Dtimes = 1, route = "sc")
    
    cost<- modCost  (model = out_A18_PO, obs = Obs_A18_PO, x ="Time", weight = w)
    cost<- modCost  (model = out_A22_SC, obs = Obs_A22_L, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A21_SC, obs = Obs_A21_L, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A21_SC, obs = Obs_A21_M, x ="Time", cost = cost, weight = w)
    
    return(cost)
}


Cost_Cattle_eva (Fit_Cattle_FLO$par, w="mean")

## Make the data for the plot
Cattle_FLO_eva <-Cost_Cattle_eva (Fit_Cattle_FLO$par, w="mean")


## Make a plot
PDat_eva   <- cbind.data.frame (OBS = Cattle_FLO_eva$residuals$obs,
                                PRE = Cattle_FLO_eva$residuals$mod,
                                RES = Cattle_FLO_eva$residuals$res)

PDat_eva <- PDat_eva %>% mutate (Log.OBS = log(OBS, 10), Log.PRE = log(PRE, 10), Species = "Cattle")


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
saveRDS(Cattle_FLO_eva, file = 'Cattle_FLO_eva.rds')








