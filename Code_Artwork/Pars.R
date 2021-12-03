## Loading R-packages
library(ggplot2)
library(extrafont)
library(ggh4x)

## Input the rds file
Fit_Swine_PG = readRDS(file = 'Fit_Swine_PG.rds')
Fit_Cattle_PG = readRDS(file = 'Fit_Cattle_PG.rds')
Fit_Swine_FLU = readRDS(file = 'Fit_Swine_FLU.rds')
Fit_Cattle_FLU = readRDS(file = 'Fit_Cattle_FLU.rds')
Fit_Swine_FLO = readRDS(file = 'Fit_Swine_FLO.rds')
Fit_Cattle_FLO = readRDS(file = 'Fit_Cattle_FLO.rds')

## Select the sensitive parameters by using sensitivity analysis
## Define the parameters; These parameters were initial values---
## --will be instead with fitted parameters
## Input the Cattle parameters
##-----------------------------------------------------------
# PG
Pars_C_PG <- c(
    fR          = 0.5,              
    fR1         = 0.5,
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
    VRestCa     = (1-0.0122-0.0021-0.1218-0.361-0.0399),
    Fracsc      = 0.50,
    Fracim      = 0.65,
    Kabs        = 1.91,
    Kunabs      = 0.01,                       
    Kim         = 0.07,                               
    Ksc         = 0.25,                              
    KehcC       = 0.05,
    KmetC       = 0.05,                          
    Kdissim     = 0.007,     
    Kdisssc     = 0.005,    
    KurineC     = 1.4,              
    KbileC      = 0.05,
    PL          = 1.21,
    PM          = 1.35,
    PF          = 0.04,
    PK          = 1.25,
    PRest       = 0.8       
)

## Input the Swine parameters
## Physiological parameters
Pars_S_PG <- c(
    fR          = 0.634,               
    fR1         = 0.634,
    BW          = 92,   
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
    VRestCa     = (1-0.0204-0.0037-0.1544-0.3632-0.0412),
    Fracsc      = 0.50,
    Fracim      = 0.65,
    Kunabs      = 0.01,                         
    Kim         = 0.07,                              
    Ksc         = 0.25,                               
    Kabs        = 1.91,
    KehcC       = 0.05,
    KmetC       = 0.05,                             
    Kdissim     = 0.007,     
    Kdisssc     = 0.005,   
    KurineC     = 1.4,             
    KbileC      = 0.05,
    PL          = 1.46,
    PM          = 0.88,
    PF          = 0.249,
    PK          = 1.068,
    PRest       = 0.73
)

#-------------------------------------------------------------------------------
## FLU
Pars_C_FLU <- c(
    ## Fixed or physiological parameters
    fR          = 0.05,              
    fR1         = 0.01,
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
    VRestCa     = (1-0.0122-0.0021-0.1218-0.361-0.0399),
    Fracsc       = 1,
    Fracim       = 1,
    ## Chemical-specific parameters
    Kunabs      = 0.5,                       
    Kabs        = 0.4,
    Kim         = 1,                               
    Ksc         = 0.40,                            
    KmetC       = 0.020,                             
    KehcC       = 0.05,      
    KurineC     = 0.1,              
    KurineC1    = 0.1,
    KbileC      = 0.5,           
    KbileC1     = 0.1, 
    PK          = 2.87,
    PK1         = 2.87,
    PF          = 0.16,
    PL          = 3.11,
    PL1         = 3.11,
    PRest1      = 1.65,
    PRest       = 1.65,
    PM          = 2.96,
    PM1         = 2.96
)

## Swine
Pars_S_FLU <- c(
    ## Fixed or physiological parameters
    fR          = 0.05,                  
    fR1         = 0.01,
    BW          = 92,   
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
    VRestCa     = (1-0.0204-0.0037-0.1544-0.3632-0.0412),
    Fracsc      = 1,
    Fracim      = 1,
    ## Chemical-specific parameters
    Kunabs      = 0.5,                        
    Kabs        = 1.9,
    Kim         = 0.5,                             
    Ksc         = 0.40,                             
    KmetC       = 0.020,                           
    KehcC       = 0.15,    
    KurineC     = 0.050,           
    KurineC1    = 0.010,
    KbileC      = 0.010,             
    KbileC1     = 0.010, 
    Kdissim     = 1e-5,     
    Kdisssc     = 1e-5,      
    PK          = 4,
    PK1         = 4,
    PF          = 0.6,
    PL          = 10.52,
    PL1         = 9.26,
    PRest1      = 5,
    PRest       = 8,
    PM          = 1.57,
    PM1         = 1.57
)

#-------------------------------------------------------------------------------
## FLO
## Cattle
Pars_C_FLO <- c(
    ## Fixed or physiological parameters
    fR          = 0.80,              
    fR1         = 0.71,
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
    VRestCa     = (1-0.0122-0.0021-0.1218-0.361-0.0399),
    Fracsc      = 0.60,
    Fracim      = 0.45,
    ## Chemical-specific parameters
    Kunabs      = 0.5,                       
    Kabs        = 1.91,
    KehcC       = 0.05,
    KmetC       = 0.0025,
    KurineC     = 0.05102,
    KurineC1    = 4.5675E-2,
    KbileC      = 1E-2,
    KbileC1     = 1E-2,
    Kim         = 0.1618,                             
    Ksc         = 0.128,
    Kdissim     = 0.0104,
    Kdisssc     = 0.0054,
    PL          = 2.2,
    PL1         = 8.74,
    PM          = 1.33,
    PM1         = 0.9,
    PF          = 0.615,
    PK          = 0.9,
    PK1         = 1.30,
    PRest       = 1.31,
    PRest1      = 0.01
)

## Swine
Pars_S_FLO <- c(
    ## Fixed or physiological parameters
    fR          = 0.81,               
    fR1         = 0.79,
    BW          = 92,   
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
    VRestCa     = (1-0.0204-0.0037-0.1544-0.3632-0.0412),
    Fracsc      = 0.50,
    Fracim      = 0.65,
    ## Chemical-specific parameters
    Kunabs      = 0.01,                        
    Kabs        = 0.4,
    Kim         = 0.1618,                             
    Ksc         = 0.128,                             
    KehcC       = 0.05,
    KmetC       = 0.25,                            
    KurineC     = 0.05102,           
    KurineC1    = 4.5675E-3, 
    KbileC      = 1E-3,      
    KbileC1     = 1E-3,     
    Kdissim     = 0.0104,   
    Kdisssc     = 0.0054,  
    PL          = 0.2,
    PL1         = 8.795916,
    PM          = 1.3,
    PM1         = 0.27, 
    PF          = 0.29,
    PK          = 3.47,
    PK1         = 1.3,
    PRest       = 0.01,
    PRest1      = 0.01
)

## Define SD list for physiological parameters
Pars_SD_Cattle <- c(
   QCC         = 1.47,    
   QLCa        = 0.25,	              
   QKCa        = 0.08,	              
   QMCa        = 0.09,	                      
   QFCa        = 0.024,
   QRestCa     = 0.027,
   VLCa        = 0.0018,                                     
   VKCa        = 0.0005,                                  
   VFCa        = 0.0506,                                    
   VMCa        = 0.117,                                  
   VbloodCa    = 0.0068,
   VRestCa     = 0.1389
)

Pars_SD_Swine <- c(
   QCC         = 1.62,    
   QLCa        = 0.082,	              
   QKCa        = 0.032,	              
   QMCa        = 0.306,	                      
   QFCa        = 0.038,
   QRestCa     = 0.042,
   VLCa        = 0.0033,                                     
   VKCa        = 0.0011,                                  
   VFCa        = 0.0265,                                    
   VMCa        = 0.0266,                                  
   VbloodCa    = 0.0046,
   VRestCa     = 0.125
)


## Extract the best optimized parameters and instead the initial parameters with them
## PG
C.par.PG <- exp(Fit_Cattle_PG$par)
S.par.PG <- exp(Fit_Swine_PG$par)

Pars_C_PG [names(C.par.PG)]  <- as.numeric(C.par.PG) 
Pars_S_PG [names(S.par.PG)]  <- as.numeric(S.par.PG) 

## FLU
C.par.FLU <- exp(Fit_Cattle_FLU$par)
S.par.FLU <- exp(Fit_Swine_FLU$par)

Pars_C_FLU [names(C.par.FLU)]  <- as.numeric(C.par.FLU) 
Pars_S_FLU [names(S.par.FLU)]  <- as.numeric(S.par.FLU) 

## FLO
C.par.FLO <- exp(Fit_Cattle_FLO$par)
S.par.FLO <- exp(Fit_Swine_FLO$par)

Pars_C_FLO [names(C.par.FLO)]  <- as.numeric(C.par.FLO) 
Pars_S_FLO [names(S.par.FLO)]  <- as.numeric(S.par.FLO) 

## Define the Monte Carlo simulation function to output "idata"
MCsim <- function(pars) {
    
  if (pars ["BW"] >100) {sd_pars <- Pars_SD_Cattle
  }else{sd_pars <- Pars_SD_Swine}
  
    ## Define the vector
    CV   <-vector()
    M    <- vector()
    Dist <- vector()
    
    ## Make a loop
    for (i in 1:length(pars)){
        
        CV[i]  = switch(substr(names(pars)[i],1,1),
                        "f"   = 0.1,   ## CV of percentage of plasma protein: 0.1 
                        "B"   = 0.2,   ## CV of BW 
                        "H"   = 0.1,   ## CV of htc: 0.1
                        "Q"   = as.numeric(sd_pars[i-4])/as.numeric(pars[i]),   ## CV of blood flow: 0.3
                        "V"   = as.numeric(sd_pars[i-4])/as.numeric(pars[i]),   ## CV of tissue volume: 0.3
                        "K"   = 0.3,   ## CV of kinetic constant: 0.3
                        "P"   = 0.2,   ## CV of partition coefficient: 0.2 
                        0.1)
        
        Dist[i] = switch(substr(names(pars)[i],1,1),
                         "f"   = 2, ## 1 indicate the normal distribution
                         "B"   = 1, ## 2 indicate the log-normal distribution
                         "H"   = 1,
                         "Q"   = 2,
                         "V"   = 1,
                         "K"   = 2,
                         "P"   = 2,
                         1)
        
        M[i] = as.numeric(pars[i])
        
    }
    
    Dat <- cbind.data.frame(Mean = M, CV = CV, Dist = Dist) %>% 
        mutate(SD    = M*CV, 
               logM  = log(M/sqrt(1+(SD/M)^2)),
               logSD = sqrt(log(1+(SD/M)^2)))
    
    ## Assign the distribution type to parameters: 1: normal; 2: lognormal
    rownames(Dat)<-names(pars)
    Data1 <- Dat%>% filter(Dist == 1)
    Data2 <- Dat%>% filter(Dist == 2)
    
    mc1<-list()
    mc2<-list()
    
    for(i in 1:dim(Data1)[1]) {
        mc1[[i]]<-rnormTrunc (N, 
                              min  = qnorm(0.025, mean = Data1$Mean[i], sd = Data1$SD[i]), 
                              max  = qnorm(0.975, mean = Data1$Mean[i], sd = Data1$SD[i]),  
                              mean = Data1$Mean[i], 
                              sd   = Data1$SD[i])
    }
    
    for(i in 1:dim(Data2)[1]) {
        mc2[[i]]<-rlnormTrunc (N, 
                               min = qlnorm(0.025, mean = Data2$logM[i], sd = Data2$logSD[i]), 
                               max = qlnorm(0.975, mean = Data2$logM[i], sd = Data2$logSD[i]), 
                               meanlog = Data2$logM[i], 
                               sdlog   = Data2$logSD[i])
    }
    
    mc1 <-do.call(cbind,mc1)
    mc2 <-do.call(cbind,mc2)
    
    colnames(mc1) <- rownames(Data1)
    colnames(mc2) <- rownames(Data2)
    idata <- cbind(ID = 1:N, mc1, mc2) 
    
    return(as.data.frame(idata))
}


## Define the prediction function
pred <- function (pars, drug, idata, tinterval = 24, Dose, Dtimes, route) {
    
    ## Exposure scenarios
    if (is.null(idata$BW) == TRUE) {BW = rep(pars["BW"], N)} else {BW =idata$BW}
    
    if (route == "im") {
        if (is.null(idata$Fracim) == TRUE) {Frac = rep(as.numeric(pars["Fracim"]), N)} 
        else {Frac =idata$Fracim}
    } else if (route == "sc") {
        if (is.null(idata$Fracsc) == TRUE) {Frac = rep(as.numeric(pars["Fracsc"]), N)} 
        else {Frac = idata$Fracsc}
    } else {Frac = 1}
    
    tinterval   = tinterval
    TDOSE       = Dtimes
    MW          = switch(drug, "PG" = 334.4, "FLO" = 358.21, "FLU" = 296.24)
    MW1         = switch(drug, "PG" = 334.4, "FLO" = 247.28, "FLU" = 312.24)
    DOSEfast    = Dose*BW*Frac/MW              
    DOSEslow    = Dose*BW*(1-Frac)/MW  
    DOSE        = Dose*BW/MW  
    
    ## 
    if (route == "im") {
        ev_1 <- ev (ID   = 1:N, amt  = DOSEfast, ii = tinterval, 
                    addl = TDOSE - 1, cmt  = "Amtsiteim", replicate = FALSE)
        ev_2 <- ev (ID   = 1:N, amt  = DOSEslow, ii = tinterval, 
                    addl = TDOSE - 1, cmt  = "ADOSEim", replicate = FALSE)
        
        ex <- ev_1 + ev_2 
        
    }
    
    if (route == "sc") {
        ev_1 <- ev (ID   = 1:N, amt  = DOSEfast, ii = tinterval, 
                    addl = TDOSE - 1, cmt  = "Amtsitesc", replicate = FALSE)
        ev_2 <- ev (ID   = 1:N, amt  = DOSEslow, ii = tinterval, 
                    addl = TDOSE - 1, cmt  = "ADOSEsc", replicate = FALSE)
        
        ex <- ev_1 + ev_2 
    }
    
    if (route == "iv") {
        ev_1 <- ev (ID   = 1:N, amt  = DOSE, ii = tinterval, tinf = 0.01,
                    addl = TDOSE - 1, cmt  = "APlas_free", replicate = FALSE)
        
        ex <- ev_1 
        
    }
    
    if (route == "po") {
        ev_1 <- ev (ID   = 1:N, amt  = DOSE, ii = tinterval, tinf = 0.01,
                    addl = TDOSE - 1, cmt  = "AST", replicate = FALSE)
        
        ex <- ev_1
    }
    
    
    tsamp  = tgrid(0, tinterval*(TDOSE - 1) + 24*60, 0.5) 
    
    ## Simulation
    out <- mod %>% param(pars)%>%
        update(atol = 1E-10, rtol = 1E-5, maxsteps = 50000) %>%
        mrgsim (idata = idata, ev = ex, tgrid = tsamp)
    
    
    outdf = cbind.data.frame( ID    = out$ID,
                              Time  = out$time,
                              CP    = out$Plasma*MW,
                              CP_Met= out$P_Met*MW1,
                              CL    = out$Liver*MW,
                              CL_Met= out$L_Met*MW1,
                              CL_Plus=out$Liver*MW+out$L_Met*MW1,
                              CK    = out$Kidney*MW,
                              CK_Met= out$K_Met*MW1,
                              CK_Plus=out$Kidney*MW+out$K_Met*MW1,
                              CM    = out$Muscle*MW,
                              CM_Met= out$M_Met*MW1,
                              CM_Plus=out$Muscle*MW+out$M_Met*MW1,
                              CF    = out$Fat*MW)
    
    
    return (outdf)
    
}
windowsFonts(Times=windowsFont("Times New Roman"))

## Plot theme
ptheme<-theme (
    plot.background         = element_rect (fill="White"),
    text                    = element_text (family = "Times"),   # text front (Time new roman)
    panel.border = element_rect(colour = "black", fill=NA, size=2.5),
    #axis.line.x = element_line(size = 2, linetype = "solid", colour = "black"),
    #axis.line.y = element_line(size = 2, linetype = "solid", colour = "black"),
    panel.background        = element_rect (fill="White"),
    panel.grid.major        = element_blank(),
    panel.grid.minor        = element_blank(),
    ggh4x.axis.ticks.length.minor = rel(0.5),
    axis.text               = element_text (size   = 15, colour = "black", face = "bold"),    # tick labels along axes
    axis.title              = element_text (size   = 18, colour = "black", face = "bold"),
    axis.ticks.length.x = unit(.25, "cm"),
    axis.ticks = element_line(size = 1),# label of axes
    legend.position='none') 

## Define the plot function
###---------------------------------------------------
MCplot <- function (data, target, TOL, tdoses, tinterval=24, color="gray", Y.limit=1e-3) {
    
    PDat <- data %>% group_by(Time)%>%
        summarise( Time = Time/24,
                   Median = switch(target,
                                   "Plasma"  = median(CP),
                                   "Plasma1" = median(CP_Met),
                                   "Liver"   = median(CL),
                                   "Liver1"  = median(CL_Met),
                                   "Liver+"  = median(CL_Plus),
                                   "Kidney"  = median(CK),
                                   "Kidney1" = median(CK_Met),
                                   "Kidney+" = median(CK_Plus),
                                   "Muscle"  = median(CM),
                                   "Muscle1" = median(CM_Met),
                                   "Muscle+" = median(CM_Plus)),
                   
                   Ub     = switch(target,
                                   "Plasma"  = quantile(CP, probs  = 0.99),
                                   "Plasma1" = quantile(CP_Met), probs  = 0.99,
                                   "Liver"   = quantile(CL, probs  = 0.99),
                                   "Liver1"  = quantile(CL_Met, probs  = 0.99),
                                   "Liver+"  = quantile(CL_Plus, probs  = 0.99),
                                   "Kidney"  = quantile(CK, probs  = 0.99),
                                   "Kidney1" = quantile(CK_Met, probs  = 0.99),
                                   "Kidney+" = quantile(CK_Plus, probs  = 0.99),
                                   "Muscle"  = quantile(CM, probs  = 0.99),
                                   "Muscle1" = quantile(CM_Met, probs  = 0.99),
                                   "Muscle+" = quantile(CM_Plus, probs  = 0.99)),
                   
                   Lb     = switch(target,
                                   "Plasma"  = quantile(CP, probs  = 0.01),
                                   "Plasma1" = quantile(CP_Met), probs  = 0.01,
                                   "Liver"   = quantile(CL, probs  = 0.01),
                                   "Liver1"  = quantile(CL_Met, probs  = 0.01),
                                   "Liver+"  = quantile(CL_Plus, probs  = 0.01),
                                   "Kidney"  = quantile(CK, probs  = 0.01),
                                   "Kidney1" = quantile(CK_Met, probs  = 0.01),
                                   "Kidney+" = quantile(CK_Plus, probs  = 0.01),
                                   "Muscle"  = quantile(CM, probs  = 0.01),
                                   "Muscle1" = quantile(CM_Met, probs  = 0.01),
                                   "Muscle+" = quantile(CM_Plus, probs  = 0.01)),
                   .groups = 'drop')
    

    x.intercept <- PDat %>%filter(Time >= ((tdoses*tinterval/24)))%>% filter (Ub <= TOL) %>%select(Time)%>%min()
    endtime     <- x.intercept + (tdoses*tinterval/24)+5
    #Y.intercept <- (PDat %>% filter (Time >=endtime)%>%select(Lb))%>% max()

    
    p <-ggplot (PDat, aes(Time, Median))+ 
        geom_ribbon(aes(ymin = Lb, ymax = Ub), fill =color)+
        geom_line  (aes(x = Time, y = Ub), colour = "black",  lwd = 1,linetype=2) +
        geom_line  (aes(x = Time, y = Lb), colour = "black",  lwd = 1,linetype=2) +
        geom_line  ( lwd=1, colour = "black") +
        scale_x_continuous (    breaks = pretty(0:endtime, n = 5),
                                minor_breaks = pretty(0:endtime, n=endtime),
                                label  = pretty((0:endtime), n = 5)-((tdoses*tinterval/24)-1),
                                limits = c(0, endtime)) +
        scale_y_log10( 
            labels = function(x) format(x, scientific = TRUE),
            limits = c(Y.limit, NA))+
        annotation_logticks(short = unit(2,"mm"),
                            mid = unit(3,"mm"),
                            long = unit(3,"mm"), size = 0.6,
                            sides = "l") +
        
        labs(x="",y="",face="bold") 
    
  
        p1<-p+
            # geom_vline(aes(xintercept = x.intercept),
            #              size = 1, color = "red", linetype = 2,
            #              show.legend = F)+
            geom_line(aes(y = TOL),color = 'black',size = 1, linetype = 'twodash', show.legend = F)+
            ptheme#+ # add tolerance line
    
    
    return(p1)
}
    
    
    







