##---------------------------------------------------------------------------------------------------
# This function generates a list of the partition coefficients based on the Schmitt method
# The code was modfied from the Metrum research group: DOI: https://doi.org/10.1124/dmd.120.090498
# Author: Wei-Chun Chou
# Adivisor: Zhoumeng Lin
# Date: 2020/10/26
##----------------------------------------------------------------------------------------------------


## Loading required packages
library(dplyr)

## The function calcuate the tissue:plasma partition coefficients based on the method---
## of the Schmitt 2008: DOI: 10.1016/j.tiv.2007.09.010
## some modified method is based on the method of Pearce et al., 2017; DOI: 10.1007/s10928-017-9548-7; used in httk package

## Abbreviations: logP: the log-scale n-octanol:buffer partition coefficient;--- 
## pka: acid disscoiation constant; BP is the in vitro blood:plasma ratio;--
## fup: fraction unbound in plasma;---
## dat: input the tissue composition data;---
## type: chemical type (#1-neutral; #2-monoprotic acid; #3-monoprotic base; #4-diprotic acid; #5-diprotic base--
## #6-monoprotic acid monoprotic base (acid comes first); #7-triprotic acid; #8-triprotic base; ---
## #9-diprotic acid monoprotic base (first two are acid);  #10-diprotic base monoprotic acid (first one is acid))

Kp_Schmit <- function(logP, pKa, fup, type = 1, dat){
    
    ## Exctract the tissue composition ---
    ## by filtering the data to exclude the "Plasma", "Red blood cells (RBCs)"
    dat_all <- dat %>% filter (!tissue %in% c("Plasma", "RBCs"))
    
    # logMA is the log of membrane affinity = phosphatidylcholine:water (neutral phospholipid:water) partition coefficient;
    # The estimation of logMA was derved from the Yun & Edginton (2013); difference with Schmitt, 2009
    
    logMA <- 1.24 + 0.304*logP  #in case we don't have a direct logMA; equation from Yun & Edginton (2013)
    K_n_pl <- 10^logMA    #neutral phospholipids:water partition coefficient
    
    # K_protein: protein:water partition; Schmitt, 2008
    # K_protein <- 0.163+0.0221*K_n_pl from schmitt, 2008; 
    K_protein <- ((0.81 + 0.11 * K_n_pl)/24.92)*5 # Here we used the equaltion lised in Utsey et al., 2020

    pH <- dat_all$pH
    alpha <- 1e-3  #ratio between ditribution coefficient at given pH (D) and that in neutral form (D0)
    logP_OW <- 1.115*logP - 1.35 ## logP_OW is the log10-scale olive oil:buffer (water) partition coefficient of nonionized species
    P_OW    <- 10^(logP_OW) ## same as "D" in PT and Berez method; olive oil:buffer (water) partition coefficient of nonionized species
    
    # Estimated the olive oil:buffer (water) partition coefficient at different pH based on Henderson-Hasselbalch Equations
    W <- switch(type,
                #1-neutral
                0,
                #2-monoprotic acid
                10^(pH-pKa),
                #3-monoprotic base
                10^(pKa-pH),
                #4-diprotic acid
                10^(pH-pKa[1])+10^(2*pH-pKa[1]-pKa[2]), 
                #5-diprotic base
                10^(pKa[2]-pH)+10^(pKa[1]+pKa[2]-2*pH), 
                #6-monoprotic acid monoprotic base (acid comes first)
                10^(pKa[2]-pH)+10^(pH-pKa[1]),
                #7-triprotic acid
                10^(pH-pKa[1])+10^(2*pH-pKa[1]-pKa[2])+10^(3*pH-pKa[1]-pKa[2]-pKa[3]),
                #8-triprotic base
                10^(pKa[1]+pKa[2]+pKa[3]-3*pH))
    
    
    ## Calculation of the neutral lipid:water partition coefficient (K_n_l) and ---
    ## acid phospholipid:water partition ceofficient (K_a_pl)
    ## neutral phospolipid:water partition coefficient (K_n_pl)
    ## In the assumption of schmitt (2008), the Do:w is equal to K_n_l
    if(type==1 | type==2 | type==4 | type==7){ # neutral, monoprotic acid, diprotic acid, triprotic acid
        K_n_l  <- P_OW*(((1 - alpha)/(1 + W)) + alpha) # eq. 13 from schmitt, 2008
        K_a_pl <- K_n_pl*((1/(1 + W)) + 0.05*(1-(1/(1 + W)))) # different with the eq. 17 from schmitt, 2008; error in the schmitt paper
    }
    else if(type==3){ # monoprotic base
        K_n_l  <- P_OW*(((1-alpha)/(1+W))+alpha) 
        K_a_pl <- K_n_pl*((1/(1 + W)) + 20*(1-(1/(1 + W)))) 
    }
    else if(type==5){ # diprotic base 
        F1 <- (1/(1+10^(pKa[1]-pH)))
        F2 <- (1/(1+10^(pKa[2]-pH)))
        K_n_l <- P_OW*(F1*F2 + alpha*((1-F1)*F2 + F1*(1-F2)) + (1-F1)*(1-F2))
        K_a_pl <- K_n_pl*(F1*F2 + 20*((1-F1)*F2 + F1*(1-F2)) + (1-F1)*(1-F2))
    }
    else if(type==6){ # monoprotic acid monoprotic base (acid comes first)
        F1 <- (1/(1+10^(pH-pKa[1])))
        F2 <- (1/(1+10^(pKa[2]-pH)))
        K_n_l <- P_OW*(F1*F2 + alpha*((1-F1)*F2 + F1*(1-F2)) + (1-F1)*(1-F2))
        K_a_pl <- K_n_pl*(F1*F2 + 0.05*(1-F1)*F2 + 20*(1-(F2))*F1 + (1-F1)*(1-F2))
    }
    else if(type==8){ # triprotic base 
        F1 <- (1/(1+10^(pKa[1]-pH)))
        F2 <- (1/(1+10^(pKa[2]-pH)))
        F3 <- (1/(1+10^(pKa[3]-pH)))
        K_n_l <- P_OW*(F1*F2*F3 + alpha*((1-F1)*F2*F3 + F1*(1-F2)*F3 + F1*F2*(1-F3) + (1-F1)*(1-F2)*F3 + (1-F1)*F2*(1-F3) + F1*(1-F2)*(1-F3) + (1-F1)*(1-F2)*(1-F3)))
        K_a_pl <- K_n_pl*(F1*F2*F3 + 20*((1-F1)*F2*F3 + F1*(1-F2)*F3 + F1*F2*(1-F3)) + (1-F1)*(1-F2)*F3 + (1-F1)*F2*(1-F3) + F1*(1-F2)*(1-F3) + (1-F1)*(1-F2)*(1-F3))
    }
    else if(type==9){ # diprotic acid monoprotic base (first two are acid)
        F1 <- (1/(1+10^(pH-pKa[1])))
        F2 <- (1/(1+10^(pH-pKa[2])))
        F3 <- (1/(1+10^(pKa[3]-pH)))
        K_n_l <- P_OW*(F1*F2*F3 + alpha*((1-F1)*F2*F3 + F1*(1-F2)*F3 + F1*F2*(1-F3) + (1-F1)*(1-F2)*F3 + (1-F1)*(F2)*(1-F3) + F1*(1-F2)*(1-F3) + (1-F1)*(1-F2)*(1-F3)))
        K_a_pl <- K_n_pl*(F1*F2*F3 + 0.05*((1-F1)*F2*F3 + F1*(1-F2)*F3 + (1-F1)*(1-F2)*F3) + 20*F1*F2*(1-F3) + (1-F1)*F2*(1-F3) + F1*(1-F2)*(1-F3) + (1-F1)*(1-F2)*(1-F3))
    }
    else if(type==10){ # diprotic base monoprotic acid (first one is acid)
        F1 <- (1/(1+10^(pH-pKa[1])))
        F2 <- (1/(1+10^(pKa[2]-pH)))
        F3 <- (1/(1+10^(pKa[3]-pH)))
        K_n_l <- P_OW*(F1*F2*F3 + alpha*((1-F1)*F2*F3 + F1*(1-F2)*F3 + F1*F2*(1-F3) + (1-F1)*(1-F2)*F3 + (1-F1)*(F2)*(1-F3) + F1*(1-F2)*(1-F3) + (1-F1)*(1-F2)*(1-F3)))
        K_a_pl <- K_n_pl*(F1*F2*F3 + 0.05*(1-F1)*F2*F3 + 20*(F1*(1-F2)*F3 + F1*F2*(1-F3) + F1*(1-F2)*(1-F3)) + (1-F1)*F2*(1-F3) + (1-F1)*(1-F2)*F3 + (1-F1)*(1-F2)*(1-F3))
    }
    
    kp <- (dat_all$f_water+(K_n_l*dat_all$f_n_l)+(K_n_pl*dat_all$f_n_pl)+(K_a_pl*dat_all$f_a_pl)+(K_protein*dat_all$f_proteins))*fup
    
    dat2 <- data.frame(tissue=dat_all$tissue, Kp=kp)
    name <- dat2$tissue %>% substr(1,2) %>% tolower()
    name <- paste("Kp", name, sep ="_")
    Kp <- as.list(dat2$Kp)
    names(Kp) <- name
    
    return(Kp)
}
    