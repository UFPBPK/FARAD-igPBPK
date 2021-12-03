## Loading required R package
library(FME)
library(mrgsolve)
library(dplyr)
library(minpack.lm)  ## R-package for model fitting
library(ggplot2)
library(tidyr)

# Input the PBPK model
source (file = "GenPBPK.R") # Loading the generic PBPK model code
source (file = "GetPC.R")   # Loading the tissue composition based tissue:plasma model

## Load Model
mod <- mcode_cache("pbpk", GenricPBPK) #refer to mcode function in mrgsolve user guide 3.1.2 Inline


## Input the Cattle parameters
## Physiological parameters (Please refer to Table 2 for the full name of parameters)
Phy_pars_Swine <- c(
  BW          = 33.182,   
  Htc         = 0.412,     
  QCC         = 8.7,    
  QLCa        = 0.273,	                        
  QKCa        = 0.114,	                    
  QMCa        = 0.342,	                       
  QFCa        = 0.128,
  QRestCa     = (1-0.273-0.114-0.342-0.128),
  VLCa        = 0.0204,                                   
  VKCa        = 0.0037,                                  
  VFCa        = 0.1544,                                     
  VMCa        = 0.3632,                                    
  VbloodCa    = 0.0412,
  VRestCa     = (1-0.0204-0.0037-0.1544-0.3632-0.0412)
)

## Chemical-specific parameters (Please refer to Table 3 for the full name of parameters)
Chme_pars_FLO <- c(
    fR         = 0.81,               
    fR1        = 0.79,
    Fracsc      = 0.60,
    Fracim      = 0.45,
    Kfeces      = 0.01,                 
    Kim         = 0.1618,                              
    Ksc         = 0.128,                            
    Kabs        = 0.4,
    Kunabs      = 0.01,
    KmetC       = 0.25,                            
    KurineC     = 0.05102,            
    KurineC1    = 4.5675E-3	, 
    KbileC      = 1E-3,         
    KbileC1     = 1E-3,     
    Kdissim     = 0.0104,     
    Kdisssc     = 0.0054      
)   

##
## Estiamte the partition coefficient
dat_tiss_Swine <- read.csv(file = "Tiss_Pig.csv")

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

PC_all_Swine<-PCcoef_all(logP = -0.12, 
                          pKa  = 9.61, 
                          fup  = (1-0.19),  
                          BP   = 1, 
                          type = 3, #3-monoprotic base
                          dat  = dat_tiss_Swine)


PC_all_Swine<- do.call(rbind,PC_all_Swine)

PC_all_Swine <- PC_all_Swine%>%as.data.frame()%>%mutate(Method = rownames(PC_all_Swine)) %>% 
    filter(Method != "pcoeff_PKsim") %>%
    summarise(PF = mean(as.numeric(Kp_ad)), 
              PL = median(as.numeric(Kp_li)), 
              PK = mean(as.numeric(Kp_ki)), 
              PM = mean(as.numeric(Kp_mu)),
              PRest = mean(as.numeric(Kp_sp, Kp_sk, Kp_lu, Kp_he, Kp_br, Kp_bo)))



PC_Swine_FLO <- c(PL = 0.2,#PC_all_Swine$PL,#1.3,
                  PL1= 8.74,
                  PK = PC_all_Swine$PK,#0.9, 
                  PK1= 1.3, 
                  PM = PC_all_Swine$PM,#0.9,
                  PM1= 0.27, 
                  PF = 0.2,#PC_all_Swine$PF,#2.7,
                  PRest1 = 0.01,
                  PRest = 0.01)#PC_all_Swine$PRest)#1.1332)


## Defined the parameters
VarPars_Swine <- c(Chme_pars_FLO, PC_Swine_FLO)   ## Variable parameters 
FixPars_Swine <- c(Phy_pars_Swine, VarPars_Swine) ## Fixed parameters

##---------------------------------------------------------------------
## Model fitting
Data_FLO <- read.csv(file = "Data_FLO.csv")
head(Data_FLO)

## Read the dataset and later used in model calibration
## B1: Study 1; Emberchts et al., 2013; IM, Matrix = serum/ plasma; Dose: 22 mg/kg/day; Rpeat: 1,1,2
Obs_B1 <- Data_FLO %>% filter(Study == 1)
Obs_B1_IM_1 <- Obs_B1 %>% filter(Dose == "22.5" & Time >=96)%>%select(Time = Time, CP = Conc)
Obs_B1_IM_2 <- Obs_B1 %>% filter(Dose == "30"& Time >=96)%>%select(Time = Time, CP = Conc)
Obs_B1_IM_3 <- Obs_B1 %>% filter(Dose == "15"& Time >=96)%>%select(Time = Time, CP = Conc)

## B2: Study 2; Jiang et al., 2006; IV, IM, PO; Matrix: plasma; Dose: 20 mg/kg/day
Obs_B2 <- Data_FLO %>% filter(Study ==2)
Obs_B2_IV <- Obs_B2 %>% filter(Route == 'IV')%>%select(Time = Time, CP = Conc)
Obs_B2_IM <- Obs_B2 %>% filter(Route == 'IM')%>%select(Time = Time, CP = Conc)
Obs_B2_PO <- Obs_B2 %>% filter(Route == 'PO')%>%select(Time = Time, CP = Conc)

## B3: Study 3; Kim et al., 2008; IM; Plasma; Dose: 5, 20 mg/kg/day
Obs_B3 <- Data_FLO %>% filter(Study ==3)
Obs_B3_IM_1 <- Obs_B3 %>% filter(Dose == 5)%>%select(Time = Time, CP = Conc)
Obs_B3_IM_2 <- Obs_B3 %>% filter(Dose == 20)%>%select(Time = Time, CP = Conc)

## B4: Study 4; Lin et al., 2002; IM; K, M, L, P; Dose: 20 mg/kg/day
Obs_B4 <- Data_FLO %>% filter(Study ==4)
Obs_B4_P <- Obs_B4 %>% filter(Matrix == 'Plasma')%>%select(Time = Time, CP = Conc)
Obs_B4_L <- Obs_B4 %>% filter(Matrix == 'Liver')%>%select(Time = Time, CL = Conc)
Obs_B4_K <- Obs_B4 %>% filter(Matrix == 'Kidney')%>%select(Time = Time, CK = Conc)
Obs_B4_M <- Obs_B4 %>% filter(Matrix == 'Muscle')%>%select(Time = Time, CM = Conc)

## B5: Study 5; Liu et al., 2003; IV, IM, PO; P; Dose: 20 mg/kg/day
Obs_B5 <- Data_FLO %>% filter(Study ==5)
Obs_B5_IV <- Obs_B5 %>% filter(Route == 'IV')%>%select(Time = Time, CP = Conc)
Obs_B5_IM <- Obs_B5 %>% filter(Route == 'IM')%>%select(Time = Time, CP = Conc)
Obs_B5_PO <- Obs_B5 %>% filter(Route == 'PO')%>%select(Time = Time, CP = Conc)

## B6: Study 6; Report, 2002; PO; Dose: 20 mg/kg/day; L, K, M, F; FLOA
Obs_B6 <- Data_FLO %>% filter(Study ==6)
Obs_B6_L <- Obs_B6 %>% filter(Matrix == 'Liver')%>%select(Time = Time, CL_Met = Conc)
Obs_B6_K <- Obs_B6 %>% filter(Matrix == 'Kidney')%>%select(Time = Time, CK_Met = Conc)
Obs_B6_M <- Obs_B6 %>% filter(Matrix == 'Muscle')%>%select(Time = Time, CM_Met = Conc)
Obs_B6_F <- Obs_B6 %>% filter(Matrix == 'Fat')%>%select(Time = Time, CF_Met = Conc)

## B7: Study 7; Report, 2006; PO; Dose: 10 mg/kg/day; L, K, M, F; Repeat; FLOA
Obs_B7 <- Data_FLO %>% filter(Study ==7)
Obs_B7_L <- Obs_B7 %>% filter(Matrix == 'Liver')%>%select(Time = Time, CL_Met = Conc)
Obs_B7_K <- Obs_B7 %>% filter(Matrix == 'Kidney')%>%select(Time = Time, CK_Met = Conc)
Obs_B7_M <- Obs_B7 %>% filter(Matrix == 'Muscle')%>%select(Time = Time, CM_Met = Conc)
Obs_B7_F <- Obs_B7 %>% filter(Matrix == 'Fat')%>%select(Time = Time, CF_Met = Conc)

## B8: Study 8; Voorspoels et al., 1999; IM, PO; Dose: 15 mg/kg/day; P; Repeat: 1, 1
Obs_B8 <- Data_FLO %>% filter(Study ==8)
Obs_B8_IM <- Obs_B8 %>% filter(Route == 'IM')%>%select(Time = Time, CP = Conc)
Obs_B8_PO <- Obs_B8 %>% filter(Route == 'PO' & Dose == 250)%>%select(Time = Time, CP = Conc)

## B9: Study 9; Zhang et al., 2016; IM; Dose: 20 mg/kg/day; P;
Obs_B9 <- Data_FLO %>% filter(Study ==9)
Obs_B9_IM <- Obs_B9 %>% filter(Route == 'IM')%>%select(Time = Time, CP = Conc)

## B10: Study 10; Lei et al., 2018; IV, IM; Dose: 30 mg/kg/day; P;
Obs_B10 <- Data_FLO %>% filter(Study ==10)
Obs_B10_IV <- Obs_B10 %>% filter(Route == 'IV')%>%select(Time = Time, CP = Conc)
Obs_B10_IM <- Obs_B10 %>% filter(Route == 'IM')%>%select(Time = Time, CP = Conc)


###
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
  
  
  tsamp  = tgrid(0, tinterval*(TDOSE - 1) + 24*20, 0.1) 
  
  ## Simulation
  out <- mod %>% param (fixpars) %>% param (parsinput) %>% 
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


## Create the Cost function
## If you run the cost function, you will find the Warning message listed below:
## In regularize.values(x, y, ties, missing(ties)) :collapsing to unique 'x' values
## It casued by the function approx in later R 3.6.0 version; But it doesn't influence the results


Cost_Swine <- function (pars, w) {
    
    # Prediction
    out_B1_IM_1  <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 36.5, Dose = 22.5, Dtimes = 1, route = "im")
    out_B1_IM_2  <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 36.5, Dose = 30, Dtimes = 1, route = "im")
    out_B1_IM_3  <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 36.5, Dose = 15, Dtimes = 2, tinterval = 48,route = "im")
    out_B2_IV    <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 27, Dose = 20, Dtimes = 1, route = "iv")
    out_B2_IM    <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 27, Dose = 20, Dtimes = 1, route = "im")
    out_B2_PO    <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 27, Dose = 20, Dtimes = 1, route = "oral")
    out_B3_IM_1  <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 36.5, Dose = 5, Dtimes = 1, route = "im")
    out_B3_IM_2  <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 36.5, Dose = 20, Dtimes = 1, route = "im")
    out_B4_IM    <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 36.5, Dose = 20, Dtimes = 1, route = "im")
    out_B5_IV    <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 36.5, Dose = 20, Dtimes = 1, route = "iv")
    out_B5_IM    <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 36.5, Dose = 20, Dtimes = 1, route = "im")
    out_B5_PO    <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 36.5, Dose = 20, Dtimes = 1, route = "oral")
    out_B6_PO    <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 59, Dose = 20, Dtimes = 5, route = "oral")
    out_B8_IM    <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 40, Dose = 15, Dtimes = 1, route = "im")
    out_B8_PO    <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 40, Dose = 15, Dtimes = 3, route = "oral")
    out_B9_IM    <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 33.6, Dose = 20, Dtimes = 1, route = "im")

    cost<- modCost  (model = out_B1_IM_1, obs = Obs_B1_IM_1, x ="Time", weight = w)
    cost<- modCost  (model = out_B1_IM_2, obs = Obs_B1_IM_2, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_B1_IM_3, obs = Obs_B1_IM_3, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_B2_IV, obs = Obs_B2_IV, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_B2_IM, obs = Obs_B2_IM, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_B2_PO, obs = Obs_B2_PO, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_B3_IM_1, obs = Obs_B3_IM_1, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_B3_IM_2, obs = Obs_B3_IM_2, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_B4_IM, obs = Obs_B4_P, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_B4_IM, obs = Obs_B4_L, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_B4_IM, obs = Obs_B4_K, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_B4_IM, obs = Obs_B4_M, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_B5_IV, obs = Obs_B5_IV, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_B5_IM, obs = Obs_B5_IM, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_B5_PO, obs = Obs_B5_PO, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_B8_IM, obs = Obs_B8_IM, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_B8_PO, obs = Obs_B8_PO, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_B9_IM, obs = Obs_B9_IM, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_B6_PO, obs = Obs_B6_L, x ="Time", cost = cost,weight = w)
    cost<- modCost  (model = out_B6_PO, obs = Obs_B6_K, x ="Time", cost = cost,weight = w)
    cost<- modCost  (model = out_B6_PO, obs = Obs_B6_M, x ="Time", cost = cost,weight = w)

    
    return(cost)
}

#Cost_Swine(log(VarPars_Swine), w="mean")
## Define the initial paraemter list
theta_init <- c(
  Fracsc      = 0.60,
  Fracim      = 0.65,
  fR          = 0.81,               
  fR1         = 0.79,
  Kim         = 0.16,                          
  Ksc         = 0.128,                           
  Kabs        = 1.9,
  KehcC        = 0.15,
  Kunabs      = 0.1,
  KmetC       = 0.025,                          
  Kmet1C      = 0.001,
  KurineC     = 0.51,          
  KurineC1    = 4.5675E-3, 
  KbileC      = 1E-2,            
  KbileC1     = 1E-3,      
  Kdissim     = 0.0104,    
  Kdisssc     = 0.0054,   
  PL          = 0.2,
  PL1         = 8.79,
  PM          = 0.7,
  PM1         = 2.27, 
  PF          = 0.29,
  PK          = 0.68,
  PK1         = 8.30,
  PRest       = 0.01
)

## Senstivity function (FME) 
## Check the senstive parameters in the model
Sns_Swine <- sensFun(func = Cost_Swine, w = "mean",
                      parms = log(theta_init), varscale = 1)

Sen_1 <- summary(Sns_Swine)
plot(summary(Sns_Swine))


## Selected senstive parameters;
theta <- theta_init[abs(Sen_1$Min) >1.1*mean(abs(Sen_1$Min))]
theta

theta <- c(
           #Fracsc      = 0.60,
           Fracim      = 0.65,
           fR          = 0.80,               
           fR1         = 0.70,
           Kim         = 0.16,                          
           #Ksc         = 0.128,                           
           #Kabs        = 1.9,
           KehcC        = 0.15,
           Kunabs      = 0.1,
           KmetC       = 0.025,                          
           #Kmet1C      = 0.001,
           KurineC     = 0.1,          
           KurineC1    = 4.5675E-3, 
           KbileC      = 1E-2,            
           KbileC1     = 1E-3,      
           Kdissim     = 0.0104,    
           #Kdisssc     = 0.0054,   
           PL          = 0.2,
           PL1         = 8.79,
           #PM          = 0.7,
           PM1         = 2.27, 
           PF          = 0.29,
           #PK          = 0.68,
           PK1         = 8.30,
           PRest       = 1.13,
           PRest1      = 0.01)

## Cost_Swine(log(theta), w="mean")

## Selected sensitive parameters; 
Fit_Swine_FLO <- modFit(f = Cost_Swine, p = log(theta), w = "mean", 
                        method ="Marq", # Nelder-Mead 
                        control = nls.lm.control(nprint = 1))

summary(Fit_Swine_FLO) 


exp(Fit_Swine_FLO$par)


###
CSwine_FLO<-Cost_Swine(Fit_Swine_FLO$par, w="mean")
a<-as.data.frame(CSwine_FLO$residuals)


PDat <- cbind.data.frame (OBS = CSwine_FLO$residuals$obs,
                          PRE = CSwine_FLO$residuals$mod,
                          RES = CSwine_FLO$residuals$res)

PDat <- PDat %>% mutate (Log.OBS = log(OBS,10), Log.PRE = log(PRE,10), Species = "Swine")

fit <- lm(Log.PRE ~ Log.OBS, data = PDat)
summary(fit)

PlotDat <- PDat %>% mutate(prediction = predict(fit), OPR = PRE/OBS)


p1 <-
    ggplot(PlotDat, aes(Log.PRE, Log.OBS)) +
    geom_point  (aes(shape   = as.factor(Species)), colour = "steelblue4", size = 4)  +
    geom_abline (intercept = 0,
                 slope     = 1,
                 color     ="steelblue4", size = 1, alpha = 0.8) +
    annotation_logticks() +
    scale_y_continuous(limits = c(-3,3), labels = scales::math_format(10^.x))+
    scale_x_continuous(limits = c(-3,3),labels = scales::math_format(10^.x))


p1

# Save the fitting results
saveRDS(Fit_Swine_FLO, file = 'Fit_Swine_FLO.rds')
saveRDS(CSwine_FLO, file = 'CSwine_FLO.rds')


##//////////////////////////////////////////////////////////////////////////////
## Model evaluation
## Study 11: Jaroszewski et al. (2008)

### Make cost function
Cost_Swine_eva <- function (pars, w) {
  
  # Prediction
  out_B7_PO    <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 59, Dose = 10, Dtimes = 5, route = "oral")
  out_B10_IV    <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 20, Dose = 30, Dtimes = 1, route = "iv")
  out_B10_IM    <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 20, Dose = 30, Dtimes = 1, route = "im")
  
  cost<- modCost  (model = out_B7_PO, obs = Obs_B7_L, x ="Time", weight = w)
  cost<- modCost  (model = out_B7_PO, obs = Obs_B7_K, x ="Time", cost = cost, weight = w)
  cost<- modCost  (model = out_B7_PO, obs = Obs_B7_M, x ="Time", cost = cost, weight = w)
  cost<- modCost  (model = out_B10_IV, obs = Obs_B10_IV, x ="Time", cost = cost, weight = w)
  cost<- modCost  (model = out_B10_IM, obs = Obs_B10_IM, x ="Time", cost = cost, weight = w)
  
  return(cost)
}


Cost_Swine_eva (Fit_Swine_FLO$par, w="mean")

## Make the data for the plot
Swine_FLO_eva <-Cost_Swine_eva (Fit_Swine_FLO$par, w="mean")


## Make a plot
PDat_eva   <- cbind.data.frame (OBS = Swine_FLO_eva$residuals$obs,
                                PRE = Swine_FLO_eva$residuals$mod,
                                RES = Swine_FLO_eva$residuals$res)

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
saveRDS(Swine_FLO_eva, file = 'Swine_FLO_eva.rds')










