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
Chme_pars_FLU <- c(
    fR          = 0.05,                  
    fR1         = 0.01,
    Kunabs      = 0.5,                         
    Kabs        = 0.4,
    Kim         = 1.00,                            
    Ksc         = 0.40,                           
    KmetC       = 0.25,                             
    KehcC       = 0.15,     
    KurineC     = 0.10,             
    KurineC1    = 0.10,
    KbileC      = 0.10,           
    KbileC1     = 0.10, 
    Kdissim     = 1e-5,     
    Kdisssc     = 1e-5      
)   

##
## Estiamte the partition coefficient
## Abbreviations: 
## ogP: the log-scale n-octanol:buffer partition coefficient;--- 
## pka: acid disscoiation constant; 
## BP is the in vitro blood:plasma ratio;--
## fup: fraction unbound in plasma;---
## dat: input the tissue composition data;---
## type: chemical type (#1-neutral; #2-monoprotic acid; #3-monoprotic base; #4-diprotic acid; #5-diprotic base--
## #6-monoprotic acid monoprotic base (acid comes first); #7-triprotic acid; #8-triprotic base; ---
## #9-diprotic acid monoprotic base (first two are acid);  #10-diprotic base monoprotic acid (first one is acid))

dat_tiss_Swine <- read.csv(file = "Tiss_Pig.csv")

PC_all_Swine<-PCcoef_all(logP = 4.17, 
                         pKa  = 1.88, 
                         fup  = (1-0.95), 
                         BP   = 1, 
                         type = 2, 
                         dat  = dat_tiss_Swine)

PC_all_Swine<- do.call(rbind,PC_all_Swine)

PC_all_Swine <- PC_all_Swine%>%as.data.frame()%>%mutate(Method = rownames(PC_all_Swine)) %>% 
    filter(Method != "pcoeff_PKsim") %>%
    summarise(PF = mean(as.numeric(Kp_ad)), 
              PL = mean(as.numeric(Kp_li)), 
              PK = mean(as.numeric(Kp_ki)), 
              PM = mean(as.numeric(Kp_mu)),
              PRest = mean(as.numeric(Kp_sp, Kp_sk, Kp_lu, Kp_he, Kp_br, Kp_bo)))



PC_Swine_FLU <- c(PL = 10.52,#PC_all_Swine$PL,
                  PL1= 9.26,#PC_all_Swine$PL,
                  PM = PC_all_Swine$PM, 
                  PM1= PC_all_Swine$PM, 
                  PK = 4,#PC_all_Swine$PM,
                  PK1= 4,#PC_all_Swine$PM, 
                  PF = 0.6,#PC_all_Swine$PF,
                  PRest1 = 5,##PC_all_Swine$PRest,
                  PRest = 8)#PC_all_Swine$PRest)


##
VarPars_Swine <- c(Chme_pars_FLU, PC_Swine_FLU)   ## Variable parameters 
FixPars_Swine <- c(PC_Swine_FLU, VarPars_Swine) ## Fixed parameters



# Loading Swine data
## Swine data
data_flu <- read.csv(file = "Data_FLU.csv")
head(data_flu)

## Study 1: Howard et al. (2014)
Obs_A1      <- data_flu %>% filter(Study ==1) 
Obs_A1_P    <- Obs_A1%>% filter(Matrix == "P") %>%select(Time = Time, CP = Conc.)
Obs_A1_PMet <- Obs_A1 %>% filter(Matrix == "P-Met") %>%select(Time = Time, CP_Met = Conc.)

## Study 2: Pairis-Garcial et al. (2010)
Obs_A2    <- data_flu %>% filter(Study ==2)
Obs_A2_IV <- Obs_A2%>% filter(Route == "IV")%>%select(Time = Time, CP = Conc.)
Obs_A2_IM <- Obs_A2%>% filter(Route == "IM")%>%select(Time = Time, CP = Conc.)
Obs_A2_PO <- Obs_A2%>% filter(Route == "PO")%>%select(Time = Time, CP = Conc.)

## Study 3: FDA (2005)
Obs_A3   <- data_flu %>% filter(Study == 3)
Obs_A3_L <- Obs_A3 %>% filter(Matrix == "L")%>%select(Time = Time, CL = Conc.)
Obs_A3_K <- Obs_A3 %>% filter(Matrix == "K")%>%select(Time = Time, CK = Conc.)
Obs_A3_M <- Obs_A3 %>% filter(Matrix == "M")%>%select(Time = Time, CM = Conc.)
Obs_A3_F <- Obs_A3 %>% filter(Matrix == "F")%>%select(Time = Time, CF = Conc.)

## Study 4: FDA (2005)
Obs_A4   <- data_flu %>% filter(Study == 4)
Obs_A4_P <- Obs_A4%>% filter(Matrix == "P") %>%select(Time = Time, CP = Conc.)

## Study 5: EMEA (1999)
Obs_A5   <- data_flu %>% filter(Study == 5)
Obs_A5_L <- Obs_A5 %>% filter(Matrix == "L")%>%select(Time = Time, CL = Conc.)
Obs_A5_M <- Obs_A5 %>% filter(Matrix == "M")%>%select(Time = Time, CM = Conc.)
Obs_A5_K <- Obs_A5 %>% filter(Matrix == "K")%>%select(Time = Time, CK = Conc.)

## Study 6: Bates et al. (2020)
Obs_A6      <- data_flu %>% filter(Study == 6)
Obs_A6_P    <- Obs_A6 %>% filter(Matrix == "P" & Comp. == 'FLU')%>%select(Time = Time, CP = Conc.)
Obs_A6_PMet <- Obs_A6 %>% filter(Matrix == "P" & Comp. == '5OH')%>%select(Time = Time, CP_Met = Conc.)
Obs_A6_L    <- Obs_A6 %>% filter(Matrix == "L"& Comp. == 'FLU')%>%select(Time = Time, CL = Conc.)
Obs_A6_LMet <- Obs_A6 %>% filter(Matrix == "L"& Comp. == '5OH')%>%select(Time = Time, CL_Met = Conc.)
Obs_A6_K    <- Obs_A6 %>% filter(Matrix == "K"& Comp. == 'FLU')%>%select(Time = Time, CK = Conc.)
Obs_A6_KMet <- Obs_A6 %>% filter(Matrix == "K"& Comp. == '5OH')%>%select(Time = Time, CK_Met = Conc.)


###
pred <- function (fixpars, Vpars, BW, tinterval = 24, Dose, Dtimes, route) {
    
    ## Get out of log domain
    parsinput <- exp(Vpars) 
    fixpars["BW"] <- BW
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

### Make cost function
Cost_Swine <- function (pars, w) {
    
    # Prediction
    out_A1    <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 40.15, Dose = 3,   Dtimes = 1, route = "iv")
    out_A2_IV <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 168,   Dose = 2.2, Dtimes = 1, route = "iv")
    out_A2_IM <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 168,   Dose = 2.2, Dtimes = 1, route = "im")
    out_A2_PO <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 168,   Dose = 2.2, Dtimes = 1, route = "oral")
    out_A3    <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 40,    Dose = 2.2, Dtimes = 3, route = "im")
    out_A4    <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 26.5,  Dose = 2,   Dtimes = 1, route = "iv")
    out_A5    <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 26.5,  Dose = 2.4, Dtimes = 3, route = "im")
    out_A6    <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 128,   Dose = 2.2, Dtimes = 1, route = "im")
    
    
    cost<- modCost  (model = out_A1,    obs = Obs_A1_P,    x ="Time", weight = w)
    cost<- modCost  (model = out_A1,    obs = Obs_A1_PMet, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A2_IV, obs = Obs_A2_IV,   x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A2_IM, obs = Obs_A2_IM,   x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A2_PO, obs = Obs_A2_PO,   x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A3,    obs = Obs_A3_L,    x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A3,    obs = Obs_A3_M,    x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A3,    obs = Obs_A3_K,    x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A3,    obs = Obs_A3_F,    x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A4,    obs = Obs_A4_P,    x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A5,    obs = Obs_A5_L,    x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A5,    obs = Obs_A5_M,    x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A5,    obs = Obs_A5_K,    x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A6,    obs = Obs_A6_P,    x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A6,    obs = Obs_A6_PMet, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A6,    obs = Obs_A6_L,    x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A6,    obs = Obs_A6_LMet, x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A6,    obs = Obs_A6_K,    x ="Time", cost = cost, weight = w)
    cost<- modCost  (model = out_A6,    obs = Obs_A6_KMet, x ="Time", cost = cost, weight = w)

    
    return(cost)
}

Cost_Swine(log(VarPars_Swine), w="mean")

## Defined the inital parameter list
theta_init <- c(
  fR          = 0.05,                  
  fR1         = 0.01,
  Kunabs      = 0.5,                     
  Kabs        = 0.4,
  Kim         = 1,                              
  Ksc         = 0.40,                         
  KmetC       = 0.02,                          
  KehcC       = 0.015, 
  KurineC     = 0.1,          
  KurineC1    = 0.01,   
  KbileC      = 0.1,          
  KbileC1     = 0.04,  
  PK          = 4,
  PK1         = 4,
  PF          = 0.6,
  PL          = 10.52,
  PL1         = 9.26,
  PRest1      = 5,
  PRest       = 8,
  PM1         = 0.4,
  PM          = 0.4
)


## Sensitivity function (FME)
## Check the Sensitivity parameters in the model
Sns_Swine <- sensFun(func = Cost_Swine ,
                     parms = log(theta_init), varscale = 1, w = "mean")

Sen_2 <- summary(Sns_Swine)
plot(summary(Sns_Swine))


# ## Selected senstive parameters;
theta <- theta_init[abs(Sen_2$Mean) > 1.2*mean(abs(Sen_2$Mean))]
theta


## Customized parameter value
theta_swine <- c(
  fR          = 0.05,                  
  #fR1         = 0.01,
  #Kunabs      = 0.5,                     
  #Kabs        = 0.4,
  #Kim         = 1,                              
  #Ksc         = 0.40,                         
  KmetC       = 0.02,                          
  KehcC       = 0.15, 
  KurineC     = 0.1,          
  KurineC1    = 0.01,   
  KbileC      = 0.1,          
  KbileC1     = 0.04,  
  PK          = 8,
  #PK1         = 4,
  #PF          = 0.6,
  PL          = 10.52,
  #PL1         = 9.26,
  #PRest1      = 5,
  PRest       = 8,
  #PM1         = 0.4,
  PM          = 0.4
  )
  

##
## Selected sensitive parameters; 

Fit_Swine_FLU <- modFit(f = Cost_Swine, 
                        p = log(theta_swine),
                        lower = -Inf, upper = log(theta_swine*10),
                        method ="Marq",w = "mean", control = nls.lm.control(nprint = 1))

summary(Fit_Swine_FLU) 


exp(Fit_Swine_FLU$par)

## Make the data for the plot
CSwine_FLU <-Cost_Swine(Fit_Swine_FLU$par, w="mean")

## Make a plot
PDat   <- cbind.data.frame (OBS = CSwine_FLU$residuals$obs,
                          PRE = CSwine_FLU$residuals$mod,
                          RES = CSwine_FLU$residuals$res)

PDat <- PDat %>% mutate (Log.OBS = log(OBS, 10), Log.PRE = log(PRE, 10), Species = "Swine")


## Estimate the r-squared values
fit <- lm(Log.OBS ~ Log.PRE, data = PDat)
summary(fit)

##
PlotDat <- PDat %>% mutate(prediction = predict(fit), OPR = PRE/OBS)


p1 <-
    ggplot(PlotDat, aes(Log.PRE, Log.OBS)) +
    geom_point  (aes(shape   = as.factor(Species)), colour = "steelblue4", size = 4)  +
    geom_abline (intercept = 0,
                 slope     = 1,
                 color     ="steelblue4", size = 1, alpha = 0.8) +
    annotation_logticks() +
    scale_y_continuous(limits = c(-4,2), labels = scales::math_format(10^.x))+
    scale_x_continuous(limits = c(-4,2),labels = scales::math_format(10^.x))

p1


## Save the fitting results
saveRDS(Fit_Swine_FLU, file = 'Fit_Swine_FLU.rds')
saveRDS(CSwine_FLU, file = 'CSwine_FLU.rds')

##//////////////////////////////////////////////////////////////////////////////
## Model evaluation
## Study 7: Kittrell et al. (2020)
Obs_A7       <- data_flu %>% filter(Study ==7) 
Obs_A7_IM    <- Obs_A7%>% filter(Route == "IM") %>%select(Time = Time, CP = Conc.)
Obs_A7_PO    <- Obs_A7 %>% filter(Route == "PO") %>%select(Time = Time, CP = Conc.)

### Make cost function
Cost_Swine_eva <- function (pars, w) {
  
  # Prediction
  out_A7_IM    <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 3.4, Dose = 2.2,   Dtimes = 1, route = "im")
  out_A7_PO    <- pred (fixpars = FixPars_Swine, Vpars = pars, BW = 3.4, Dose = 2.2, Dtimes = 1, route = "iv")
  
  
  cost<- modCost  (model = out_A7_IM,    obs = Obs_A7_IM, x ="Time", weight = w)
  
  cost<- modCost  (model = out_A7_PO,    obs = Obs_A7_PO, x ="Time", cost = cost, weight = w)
  
  
  return(cost)
}


Cost_Swine_eva (Fit_Swine_FLU$par, w="mean")


## Make the data for the plot
CSwine_FLU_eva  <-Cost_Swine_eva(Fit_Swine_FLU$par, w="mean")

## Make a plot
PDat_eva   <- cbind.data.frame (OBS = CSwine_FLU_eva$residuals$obs,
                            PRE = CSwine_FLU_eva$residuals$mod,
                            RES = CSwine_FLU_eva$residuals$res)

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

## Save the fitting results
saveRDS(CSwine_FLU_eva, file = 'CSwine_FLU_eva.rds')











