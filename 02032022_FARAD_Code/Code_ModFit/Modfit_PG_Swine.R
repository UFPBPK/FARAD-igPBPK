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


## Input the swine parameters
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
Chme_pars_PG <- c(
    fR          = 0.84,                   
    fR1         = 0.84,
    Fracsc      = 0.5,
    Fracim      = 0.65,
    Kunabs      = 0.01,                         
    Kim         = 0.07,                              
    Ksc         = 0.25,                             
    Kabs        = 1.91,
    Kdissim     = 0.007,     
    Kdisssc     = 0.005,      
    KehcC       = 0.05,
    KmetC       = 0.05,                             
    KurineC     = 1.4,           
    KbileC      = 0.5
)   

##
## Estiamte the partition coefficient
dat_tiss_Swine <- read.csv(file = "Tiss_Pig.csv")

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

PC_all_Swine<-PCcoef_all(logP = 1.83, 
                         pKa  = 2.74, 
                         fup  = (1-0.366), 
                         BP   = 1, 
                         type = 2, #3-monoprotic base
                         dat  = dat_tiss_Swine)


PC_all_Swine<- do.call(rbind,PC_all_Swine)

PC_all_Swine <- PC_all_Swine%>%as.data.frame()%>%mutate(Method = rownames(PC_all_Swine)) %>% 
    filter(Method != "pcoeff_PKsim") %>%
    summarise(PF = mean(as.numeric(Kp_ad)), 
              PL = mean(as.numeric(Kp_li)), 
              PK = mean(as.numeric(Kp_ki)), 
              PM = mean(as.numeric(Kp_mu)),
              PRest = mean(as.numeric(Kp_sp, Kp_sk, Kp_lu, Kp_he, Kp_br, Kp_bo)))




PC_Swine_PG <- c(PL = PC_all_Swine$PL,
                 PL1= PC_all_Swine$PL,
                 PK = PC_all_Swine$PK, 
                 PK1= PC_all_Swine$PK, 
                 PM = PC_all_Swine$PM,
                 PM1= PC_all_Swine$PM, 
                 PF = PC_all_Swine$PF,
                 PRest1 = PC_all_Swine$PRest,
                 PRest = PC_all_Swine$PRest)

## Defined the parameters
VarPars_Swine <- c(Chme_pars_PG, PC_Swine_PG)   ## Variable parameters 
FixPars_Swine <- c(Phy_pars_Swine, VarPars_Swine) ## Fixed parameters

##---------------------------------------------------------------------
## Model fitting
Data_PG <- read.csv(file = "Data_PG.csv")
head(Data_PG)

## Read the dataset and later used in model calibration
## Study 1: Ranheim et al., 2002; IM,SC; Matrix = serum/ plasma; Dose: 15.5, 11.7 mg/kg/day; Repeat: 1,1
Obs_A1    <- Data_PG %>% filter(Study == 1)
Obs_A1_IM <- Obs_A1 %>% filter(Route == "IM")%>%select(Time = Time, CP = Conc)
Obs_A1_SC <- Obs_A1 %>% filter(Route == "SC")%>%select(Time = Time, CP = Conc)

## Study 2: Korsrud et al., 1998; IM; Matrix: P, L, K, M, F; Dose: 15, 65 mg/kg/day; Rpeat: 3
Obs_A2_1    <- Data_PG %>% filter(Study ==2 & Dose == "14.9")
Obs_A2_2    <- Data_PG %>% filter(Study ==2 & Dose == "65.4")
Obs_A2_1_P  <- Obs_A2_1 %>% filter(Matrix == 'P')%>%select(Time = Time, CP = Conc)
Obs_A2_1_K  <- Obs_A2_1 %>% filter(Matrix == 'K')%>%select(Time = Time, CK = Conc)
Obs_A2_1_M  <- Obs_A2_1 %>% filter(Matrix == 'M')%>%select(Time = Time, CM = Conc)
Obs_A2_2_P  <- Obs_A2_2 %>% filter(Matrix == 'P')%>%select(Time = Time, CP = Conc)
Obs_A2_2_L  <- Obs_A2_2 %>% filter(Matrix == 'L')%>%select(Time = Time, CL = Conc)
Obs_A2_2_K  <- Obs_A2_2 %>% filter(Matrix == 'K')%>%select(Time = Time, CK = Conc)
Obs_A2_2_M  <- Obs_A2_2 %>% filter(Matrix == 'M')%>%select(Time = Time, CM = Conc)

## Study 3: Korsrud et al., 1996; IM; Plasma; Dose: 15 mg/kg/day; Rpeat: 3
Obs_A3    <- Data_PG %>% filter(Study ==3)
Obs_A3_P  <- Obs_A3 %>% filter(Matrix == "P")%>%select(Time = Time, CP = Conc)
Obs_A3_M  <- Obs_A3 %>% filter(Matrix == "M")%>%select(Time = Time, CM = Conc)
Obs_A3_K  <- Obs_A3 %>% filter(Matrix == "K")%>%select(Time = Time, CK = Conc)

## Study 4: Apley et al., 2009: Study 4; IM; P,M,L,K,F; Dose: 32.5 mg/kg/day; Repeat: 1
# Obs_A4    <- Data_PG %>% filter(Study == 4)
# Obs_A4_P  <- Obs_A4 %>% filter(Matrix == 'P')%>%select(Time = Time, CP = Conc)
# Obs_A4_M  <- Obs_A4 %>% filter(Matrix == 'M')%>%select(Time = Time, CM = Conc)
# Obs_A4_K  <- Obs_A4 %>% filter(Matrix == 'K')%>%select(Time = Time, CK = Conc)
# Obs_A4_L  <- Obs_A4 %>% filter(Matrix == 'L')%>%select(Time = Time, CL = Conc)

## Study 5:Lupton et al., 2014; IM; P,M,K; Dose: 32.5 mg/kg/day; Repeat: 3
Obs_A5    <- Data_PG %>% filter(Study == 5)
Obs_A5_P  <- Obs_A5 %>% filter(Matrix == 'P')%>%select(Time = Time, CP = Conc)
Obs_A5_M  <- Obs_A5 %>% filter(Matrix == 'M')%>%select(Time = Time, CM = Conc)
Obs_A5_K  <- Obs_A5 %>% filter(Matrix == 'K')%>%select(Time = Time, CK = Conc)


## Study 6:Li et al., 2019; 
Obs_A6    <- Data_PG %>% filter(Study == 6)
Obs_A6_P_1  <- Obs_A6 %>% filter(Matrix == 'P' & Dose == 6.5)%>%select(Time = Time, CP = Conc)
Obs_A6_P_2  <- Obs_A6 %>% filter(Matrix == 'P'& Dose == 32.5)%>%select(Time = Time, CP = Conc)
Obs_A6_M_1  <- Obs_A6 %>% filter(Matrix == 'M'& Dose == 6.5)%>%select(Time = Time, CM = Conc)
Obs_A6_M_2  <- Obs_A6 %>% filter(Matrix == 'M'& Dose == 32.5)%>%select(Time = Time, CM = Conc)
Obs_A6_L_2  <- Obs_A6 %>% filter(Matrix == 'L'& Dose == 32.5)%>%select(Time = Time, CL = Conc)
Obs_A6_K_1  <- Obs_A6 %>% filter(Matrix == 'K'& Dose == 6.5)%>%select(Time = Time, CK = Conc)
Obs_A6_K_2  <- Obs_A6 %>% filter(Matrix == 'K'& Dose == 32.5)%>%select(Time = Time, CK = Conc)


###
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
        update(atol = 1E-12, rtol = 1E-6, maxsteps = 50000) %>%
        mrgsim_d (data = ex, tgrid = tsamp)
    
    
    outdf = cbind.data.frame( Time  = out$time,
                              CP    = out$Plasma*MW,
                              CL    = out$Liver*MW,
                              CK    = out$Kidney*MW,
                              CM    = out$Muscle*MW,
                              CF    = out$Fat*MW)
    return (outdf)
    
}

## Create the Cost function
## If you run the cost function, you will find the Warning message listed below:
## In regularize.values(x, y, ties, missing(ties)) :collapsing to unique 'x' values
## It casued by the function approx in later R 3.6.0 version; But it doesn't influence the results


Cost_Swine <- function (pars, w) {
    
    # Prediction
    out_A1_IM   <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 3.3, Dose = 99, Dtimes = 1, route = "im")
    out_A1_SC   <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 3.6, Dose = 99, Dtimes = 1, route = "sc")
    out_A2_IM_1 <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 100, Dose = 14.9, Dtimes = 3, route = "im")
    out_A2_IM_2 <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 100, Dose = 65.4, Dtimes = 5, route = "im")
    out_A3_IM   <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 100, Dose = 14.9, Dtimes = 3, route = "im")
    out_A5_IM   <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 228, Dose = 32.5, Dtimes = 3, route = "im")
   
    
    cost<- modCost  (model = out_A1_IM, obs = Obs_A1_IM, x ="Time", weight = w)
    cost<- modCost  (model = out_A1_SC, obs = Obs_A1_SC, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A2_IM_1, obs = Obs_A2_1_P, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A2_IM_1, obs = Obs_A2_1_K, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A2_IM_1, obs = Obs_A2_1_M, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A2_IM_2, obs = Obs_A2_2_P, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A2_IM_2, obs = Obs_A2_2_L, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A2_IM_2, obs = Obs_A2_2_K, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A2_IM_2, obs = Obs_A2_2_M, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A3_IM, obs = Obs_A3_P, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A3_IM, obs = Obs_A3_M, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A3_IM, obs = Obs_A3_K, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A5_IM, obs = Obs_A5_P, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A5_IM, obs = Obs_A5_M, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A5_IM, obs = Obs_A5_K, x ="Time", cost = cost, weight = w)

    return(cost)
}

Cost_Swine(log(VarPars_Swine), w="mean")

## Sensitivity function (FME) 
## Check the sensitive parameters in the model
# Sns_Swine <- sensFun(func = Cost_Swine, w = "mean",
#                       parms = log(VarPars_Swine), varscale = 1)
# 
# Sen_1 <- summary(Sns_Swine)
# plot(summary(Sns_Swine))

# ## Selected senstive parameters;
# theta_Swine <- VarPars_Swine[abs(Sen_1$Mean) >1.2*mean(abs(Sen_1$Mean))]
# theta_Swine

###
theta_Swine <- c(
    fR          = 0.64,                   
    #fR1        = 0.834,
    #Fracsc      = 0.50,
    Fracim      = 0.6,
    Kunabs      = 0.8,                       
    Kim         = 0.07,                              
    #Ksc        = 0.02,                               
    #Kabs       = 0.40,
    KehcC       = 0.01,
    KmetC       = 0.5,                            
    #Kdissim     = 0.007,     
    #Kdisssc    = 0.005,    
    KurineC     = 1.4,              
    KbileC      = 0.5,           
    PL          = 0.08,
    PM          = 0.15,
    #PF          = 0.35,
    PK          = 1.2,
    PRest       = 0.479)


# ## Selected sensitive parameters; 


Fit_Swine_PG <- modFit(f = Cost_Swine, 
                       p = log(theta_Swine), 
                       w = "mean",
                       method ="Marq",
                       lower = -Inf, upper = log(theta_Swine*10),
                       control = nls.lm.control(nprint = 1))

summary(Fit_Swine_PG) 


exp(Fit_Swine_PG$par)

###
CSwine_PG<-Cost_Swine(Fit_Swine_PG$par, w="mean")

## Make a plot
PDat_Swine <- cbind.data.frame (
    OBS = CSwine_PG$residuals$obs,
    PRE = CSwine_PG$residuals$mod,
    RES = CSwine_PG$residuals$res)

PDat_Swine <- PDat_Swine %>% filter(OBS>0)%>%mutate (Log.OBS = log(OBS,10), Log.PRE = log(PRE,10), Species = "Swine")

fit_Swine <- lm(Log.PRE ~ Log.OBS, data = PDat_Swine)
summary(fit_Swine)

PlotDat_Swine <- PDat_Swine %>% mutate(prediction = predict(fit_Swine), OPR = PRE/OBS)
    

p1 <-
    ggplot(PlotDat_Swine, aes(Log.PRE, Log.OBS)) +
    geom_point  (aes(shape   = as.factor(Species)), colour = "steelblue4", size = 4)  +
    geom_abline (intercept = 0,
                 slope     = 1,
                 color     ="steelblue4", size = 1, alpha = 0.8) +
    annotation_logticks() +
    scale_y_continuous(limits = c(-3,1), labels = scales::math_format(10^.x))+
    scale_x_continuous(limits = c(-3,1),labels = scales::math_format(10^.x))
p1

## Save the fitting results
saveRDS(Fit_Swine_PG, file = 'Fit_Swine_PG.rds')
saveRDS(CSwine_PG, file = 'CSwine_PG.rds')

#///////////////////////////////////////////////////////////////////////////////
## Model evaluation
### Make cost function
Cost_Swine_eva <- function (pars, w) {
    
    out_A6_IM_1   <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 238, Dose = 6.5, Dtimes = 3, route = "im")
    out_A6_IM_2   <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 238, Dose = 32.5, Dtimes = 3, route = "im")

    cost<- modCost  (model = out_A6_IM_1, obs = Obs_A6_P_1, x ="Time", cost = cost,weight = w)
    cost<- modCost  (model = out_A6_IM_1, obs = Obs_A6_M_1, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A6_IM_1, obs = Obs_A6_K_1, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A6_IM_2, obs = Obs_A6_P_2, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A6_IM_2, obs = Obs_A6_M_2, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A6_IM_2, obs = Obs_A6_K_2, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A6_IM_2, obs = Obs_A6_L_2, x ="Time", cost = cost, weight = w)
    
    
    return(cost)
}


Cost_Swine_eva (Fit_Swine_PG$par, w="mean")


## Make the data for the plot
CSwine_PG_eva  <-Cost_Swine_eva(Fit_Swine_PG$par, w="mean")

## Make a plot
PDat_eva   <- cbind.data.frame (OBS = CSwine_PG_eva$residuals$obs,
                                PRE = CSwine_PG_eva$residuals$mod,
                                RES = CSwine_PG_eva$residuals$res)

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
saveRDS(CSwine_PG_eva, file = 'CSwine_PG_eva.rds')












