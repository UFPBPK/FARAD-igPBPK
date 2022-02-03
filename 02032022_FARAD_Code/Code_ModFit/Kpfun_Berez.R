##---------------------------------------------------------------------------------------------------
# This function generates a list of the partition coefficients based on the BEREZHKOVSKIY, LEONID M. (2004) equations
# The code was modfied from the Metrum research group: DOI: 10.1002/jps.20073
# Author: Wei-Chun Chou
# Adivisor: Zhoumeng Lin
# Date: 2020/10/23
##----------------------------------------------------------------------------------------------------

## Loading required packages
.libPaths("C:/Users/weichunc/OneDrive - Kansas State University/Documents/R/win-library/3.6") 
library(dplyr)

## The function calcuate the tissue:plasma partition coefficients based on the method---
## of the BEREZHKOVSKIY, LEONID M. (2004) equations: DOI: 10.1002/jps.20073

## Abbreviations: logP: the log-scale n-octanol:buffer partition coefficient;--- 
## pka: acid disscoiation constant; 
## BP is the Blood to plasma concentration ratio;--
## fup: fraction unbound in plasma;---
## dat: input the tissue composition data;---
## type: chemical type (#1-neutral; #2-monoprotic acid; #3-monoprotic base; #4-diprotic acid; #5-diprotic base--
## #6-monoprotic acid monoprotic base (acid comes first); #7-triprotic acid; #8-triprotic base; ---
## #9-diprotic acid monoprotic base (first two are acid);  #10-diprotic base monoprotic acid (first one is acid))

Kp_Berez <- function (logP, pKa, fup, BP = 1, type = 1, dat) {
    
    ## Exctract the tissue composition ---
    ## by filtering the data to exclude the "Plasma", "Adipose", "Red blood cells (RBCs)"
    dat_all <- dat %>% filter (!tissue %in% c("Plasma", "Adipose", "RBCs"))
    
    ## Create a empty vector to save the partition coefficient for each organs
    n <- length (dat$tissue)
    Kp_all <- vector(mode = "numeric", length = n) 
    
    Vwp  <- dat$f_water[dat$tissue == "Plasma"] # Extract the fraction volume of water in plasma
    Vnlp <- dat$f_n_l[dat$tissue == "Plasma"]   # Extract the fraction volume of fractional volume of neutral lipids in plasma
    Vphp <- dat$f_pl[dat$tissue == "Plasma"]    # Extract the fraction volume of fractional volume of phospholipids in plasma
    
    
    ## Extract the tissue composition ---
    ## by filtering the data to exclude the "Plasma", and "RBCs" 
    dat_tis <- dat %>% filter (!tissue %in% c("Plasma", "RBCs"))
    
    ## Extract the f_water, f_n_l and f_pl of tissues, except for the adipose
    Vwt   <- dat_tis$f_water[dat_tis$tissue != 'Adipose'] # Extract the fraction volume of water in tissues
    Vnlt  <- dat_tis$f_n_l[dat_tis$tissue != 'Adipose'] # Extract the fraction volume of fractional volume of neutral lipids in tissues
    Vpht  <- dat_tis$f_pl[dat_tis$tissue != 'Adipose'] # Extract the fraction volume of fractional volume of phospholipids in tissues
    
    ## Extract the f_water, f_n_l and f_pl of adipose
    Vwad   <- dat_tis$f_water[dat_tis$tissue == 'Adipose']
    Vnlad  <- dat_tis$f_n_l[dat_tis$tissue == 'Adipose']
    Vphad  <- dat_tis$f_pl[dat_tis$tissue == 'Adipose']    
    
    ## Estiamte the the D_star: "olive oil:buffer partition coefficient" of both the--
    ## nonionized and ionized species at pH 7.4
    
    pH <- dat$pH[dat$tissue == "Adipose"]
    logD <- 1.115*logP-1.35 # log10-scale olive oil:buffer (water) partition coefficient of nonionized species
    
    ## Estiamte the logD_star based on different type of chemicals
    ## Equations are based on the Pulin and Theil, 2001: DOI: 10.1002/1520-6017(200104)90:4<436::aid-jps1002>3.0.co;2-p 
    
    logD_star <- switch(type,
                        #1-neutral
                        logD,
                        #2-monoprotic acid
                        logD-log10(1+10^(pH-pKa)),
                        #3-monoprotic base
                        logD-log10(1+10^(pKa-pH)), 
                        #4-diprotic acid
                        logD-log10(1+10^(2*pH-pKa[1]-pKa[2])),
                        #5-diprotic base
                        logD-log10(1+10^(pKa[1]+pKa[2]-2*pH)), 
                        #6-monoprotic acid monoprotic base (acid comes first)
                        logD-log10(1+10^(pKa[2]-pKa[1])),  
                        #7-triprotic acid
                        logD-log10(1+10^(3*pH-pKa[1]-pKa[2]-pKa[3])),  
                        #8-triprotic base
                        logD-log10(1+10^(pKa[1]+pKa[2]+pKa[3]-3*pH)),  
                        #9-diprotic acid monoprotic base (first two are acid)
                        logD-log10(1+10^(pH-pKa[1]-pKa[2]+pKa[3])), 
                        #10-diprotic base monoprotic acid (first one is acid)
                        logD-log10(1+10^(pKa[2]+pKa[3]-pKa[1]-pH)))
    
    
    ## Estimate the partition coefficient (Kp) of adipose (Kp_ad) and tissues (Kp_tis) 
    fut <- 1/(1+((1-fup)/fup)*0.5) # the unbound fraction of drug/chemical in tissues (fut) or plasma (fup)
    D_star <- 10^logD_star
    Kp_ad <- ((D_star*(Vnlad+0.3*Vphad)+((Vwad/fut)+0.7*Vphad))/(D_star*(Vnlp+0.3*Vphp)+((Vwp/fup)+0.7*Vphp)))
    
    P <- 10^logP  # P: n-octanol:buffer partition coefficient of the nonionized species at pH 7.4,
    Kp_tis <- ((P*(Vnlt+0.3*Vpht)+((Vwt/fut)+0.7*Vpht))/(P*(Vnlp+0.3*Vphp)+((Vwp/fup)+0.7*Vphp))) 
    
    ## Name the Kp
    ## tolower function is used to translate characters vector to small case
    nms_all <- dat_all$tissue %>% substr(1,2) %>% tolower()
    nms_all <- paste("Kp", nms_all, sep ="_")
    nms <- c("Kp_ad",nms_all)
    Kp <- as.list(c(Kp_ad, Kp_tis))
    names(Kp) <- nms
    
    return (Kp)
    
}









