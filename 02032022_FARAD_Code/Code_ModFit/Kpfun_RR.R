##---------------------------------------------------------------------------------------------------
# This function generates a list of the partition coefficients based on the Rodgers and Rowland method
# The code was modfied from the Metrum research group: DOI: https://doi.org/10.1124/dmd.120.090498
# Author: Wei-Chun Chou
# Adivisor: Zhoumeng Lin
# Date: 2020/10/26
##----------------------------------------------------------------------------------------------------


## Loading required packages
library(dplyr)

## The function calcuate the tissue:plasma partition coefficients based on the method---
## of the Rodgers and Rowland 2005: https://doi.org/10.1002/jps.20322 ---
## and Rodgers and Rowland 2006: https://doi.org/10.1002/jps.20502

## Abbreviations: logP: the log-scale n-octanol:buffer partition coefficient;--- 
## pka: acid disscoiation constant; 
## BP is the Blood to plasma concentration ratio;--
## fup: fraction unbound in plasma;---
## dat: input the tissue composition data;---
## type: chemical type (#1-neutral; #2-monoprotic acid; #3-monoprotic base; #4-diprotic acid; #5-diprotic base--
## #6-monoprotic acid monoprotic base (acid comes first); #7-triprotic acid; #8-triprotic base; ---
## #9-diprotic acid monoprotic base (first two are acid);  #10-diprotic base monoprotic acid (first one is acid))

Kp_RR <- function (logP, pKa = 0, fup, BP = 1, type = 1, dat) {
    
    ## Exctract the tissue composition ---
    ## by filtering the data to exclude the "Plasma", "Adipose", "Red blood cells (RBCs)"
    dat_all <- dat %>% filter (!tissue %in% c("Plasma", "Adipose", "RBCs"))
    
    ## by filtering the data to get the tissue composition for "Adipose"
    dat_ad <- dat %>% filter (tissue == "Adipose")
    ## by filtering the data to get the tissue composition for "RBCs"
    dat_rbc <- dat %>% filter(tissue == "RBCs") 
    ## by filtering the data to get the tissue composition for "Plasma"
    dat_plas <- dat %>% filter(tissue == "Plasma") 
    
    pH_IW   <- 7           ## pH of intracellular tissue water
    pH_P    <- 7.4             ## pH of plasma
    pH_RBC  <- 7.22            ## pH of blood cells
    P       <- 10^(logP)       ## octonal:water partition coeff
    logP_OW <- 1.115*logP - 1.35 ## logP_OW is the log10-scale olive oil:buffer (water) partition coefficient of nonionized species
    P_OW    <- 10^(logP_OW) ## same as "D" in PT and Berez method; olive oil:buffer (water) partition coefficient of nonionized species
    Ka      <- 10^(-pKa)      ## Acid dissociation constant
    HCT     <- 0.45          ## hematocrit
    
    #Calculate Kp values
    Kpu_bc <- (HCT - 1 + BP)/(HCT*fup) ## Blood cell to plasma water concentration ratio
    
    X <- switch(type,
                #1-neutral
                0,   
                #2-monoprotic acid
                10^(pH_IW-pKa),
                #3-monoprotic base
                10^(pKa-pH_IW),
                #4-diprotic acid
                10^(pH_IW-pKa[1])+10^(2*pH_IW-pKa[1]-pKa[2]),
                #5-diprotic base
                10^(pKa[2]-pH_IW)+10^(pKa[1]+pKa[2]-2*pH_IW), 
                #6-monoprotic acid monoprotic base (acid comes first)
                10^(pKa[2]-pH_IW)+10^(pH_IW-pKa[1]),  
                #7-triprotic acid
                10^(pH_IW-pKa[1])+10^(2*pH_IW-pKa[1]-pKa[2])+10^(3*pH_IW-pKa[1]-pKa[2]-pKa[3]),  
                #8-triprotic base
                10^(pKa[3]-pH_IW)+10^(pKa[3]+pKa[2]-2*pH_IW)+10^(pKa[1]+pKa[2]+pKa[3]-3*pH_IW),  
                #9-diprotic acid monoprotic base (first two are acid)
                10^(pKa[3]-pH_IW)+10^(pH_IW-pKa[1])+10^(2*pH_IW-pKa[1]-pKa[2]), 
                #10-diprotic base monoprotic acid (first one is acid)
                10^(pH_IW-pKa[1])+10^(pKa[3]-pH_IW)+10^(pKa[2]+pKa[3]-2*pH_IW))       
    
    Y <- switch(type,
                #1-neutral
                0,   
                #2-monoprotic acid
                10^(pH_P-pKa),
                #3-monoprotic base
                10^(pKa-pH_P), 
                #4-diprotic acid
                10^(pH_P-pKa[1])+10^(2*pH_P-pKa[1]-pKa[2]),
                #5-diprotic base
                10^(pKa[2]-pH_P)+10^(pKa[1]+pKa[2]-2*pH_P), 
                #6-monoprotic acid monoprotic base (acid comes first)
                10^(pKa[2]-pH_P)+10^(pH_P-pKa[1]),  
                #7-triprotic acid
                10^(pH_P-pKa[1])+10^(2*pH_P-pKa[1]-pKa[2])+10^(3*pH_P-pKa[1]-pKa[2]-pKa[3]),  
                #8-triprotic base
                10^(pKa[3]-pH_P)+10^(pKa[3]+pka[2]-2*pH_P)+10^(pKa[1]+pKa[2]+pKa[3]-3*pH_P),  
                #9-diprotic acid monoprotic base (first two are acid)
                10^(pKa[3]-pH_P)+10^(pH_P-pKa[1])+10^(2*pH_P-pKa[1]-pKa[2]), 
                #10-diprotic base monoprotic acid (first one is acid)
                10^(pH_P-pKa[1])+10^(pKa[3]-pH_P)+10^(pKa[2]+pKa[3]-2*pH_P))       
    
    Z <- switch(type,
                #1-neutral
                1,   
                #2-monoprotic acid
                1,
                #3-monoprotic base
                10^(pKa-pH_RBC), 
                #4-diprotic acid
                1,
                #5-diprotic base
                10^(pKa[2]-pH_RBC)+10^(pKa[1]+pKa[2]-2*pH_RBC), 
                #6-monoprotic acid monoprotic base (acid comes first)
                10^(pKa[2]-pH_RBC)+10^(pH_RBC-pKa[1]),  
                #7-triprotic acid
                1,  
                #8-triprotic base
                10^(pKa[3]-pH_RBC)+10^(pKa[3]+pka[2]-2*pH_RBC)+10^(pKa[1]+pKa[2]+pKa[3]-3*pH_RBC),  
                #9-diprotic acid monoprotic base (first two are acid)
                10^(pKa[3]-pH_RBC)+10^(pH_RBC-pKa[1])+10^(2*pH_RBC-pKa[1]-pKa[2]), 
                #10-diprotic base monoprotic acid (first one is acid)
                10^(pH_RBC-pKa[1])+10^(pKa[3]-pH_RBC)+10^(pKa[2]+pKa[3]-2*pH_RBC)) 
    
    
    Ka_PR <- (1/fup - 1 - (P*dat_plas$f_n_l + (0.3*P + 0.7)*dat_plas$f_n_pl)/(1+Y))
    Ka_AP <- (Kpu_bc - (1 + Z)/(1 + Y)*dat_rbc$f_iw - (P*dat_rbc$f_n_l + (0.3*P + 0.7)*dat_rbc$f_n_pl)/(1 + Y)) * (1 + Y)/dat_rbc$f_a_pl/Z
    
    
    # Assign the moderate to strong bases type_calc=1 and everything else type_calc=2 
    type_calc <- ifelse((type==3 & pKa[1]>7) | (type==5 & pKa[1] >7) | (type==6 & pKa[2] > 7) | (type==8 & pKa[1] > 7) | (type==9 & pKa[3]>7) | (type==10 & pKa[2]>7), 1,2)
    
    # Re-assign the neutrals type_calc=3
    if(type==1){type_calc=3}  #neutrals
    
    
    # Multiply by fup to get Kp rather than Kpu
    if(type_calc==1){  #moderate to strong bases
        Kp_tis <- (dat_all$f_ew + ((1 + X)/(1 + Y))*dat_all$f_iw + ((P*dat_all$f_n_l + (0.3*P + 0.7)*dat_all$f_n_pl))/(1 + Y) + (Ka_AP*dat_all$f_a_pl*X)/(1 + Y))*fup  #non lipid
        Kp_ad <- (dat_ad$f_ew + ((1 + X)/(1 + Y))*dat_ad$f_iw + ((P_OW*dat_ad$f_n_l + (0.3*P_OW + 0.7)*dat_ad$f_n_pl))/(1 + Y) + (Ka_AP*dat_ad$f_a_pl*X)/(1 + Y))*fup  #lipid
    }else if(type_calc==2){   #acidic and zwitterions
        Kp_tis <- (dat_all$f_ew + ((1 + X)/(1 + Y))*dat_all$f_iw + ((P*dat_all$f_n_l + (0.3*P + 0.7)*dat_all$f_n_pl))/(1 + Y) + (Ka_PR*dat_all$AR*X)/(1 + Y))*fup  #non lipid
        Kp_ad <- (dat_ad$f_ew + ((1 + X)/(1 + Y))*dat_ad$f_iw + ((P_OW*dat_ad$f_n_l + (0.3*P_OW + 0.7)*dat_ad$f_n_pl))/(1 + Y) + (Ka_PR*dat_ad$AR*X)/(1 + Y))*fup #lipid
    }else{  #neutrals
        Kp_tis <- (dat_all$f_ew + ((1 + X)/(1 + Y))*dat_all$f_iw + ((P*dat_all$f_n_l + (0.3*P + 0.7)*dat_all$f_n_pl))/(1 + Y) + (Ka_PR*dat_all$LR*X)/(1 + Y))*fup  #non lipid
        Kp_ad <- (dat_ad$f_ew + ((1 + X)/(1 + Y))*dat_ad$f_iw + ((P_OW*dat_ad$f_n_l + (0.3*P_OW + 0.7)*dat_ad$f_n_pl))/(1 + Y) + (Ka_PR*dat_ad$LR*X)/(1 + Y))*fup  #lipid
    }
    
    
    ## Name the Kp
    ## tolower function is used to translate characters vector to small case
    nms_all <- dat_all$tissue %>% substr(1,2) %>% tolower()
    nms_all <- paste("Kp", nms_all, sep ="_")
    nms <- c("Kp_ad",nms_all)
    Kp <- as.list(c(Kp_ad, Kp_tis))
    names(Kp) <- nms
    
    return (Kp)
    
}