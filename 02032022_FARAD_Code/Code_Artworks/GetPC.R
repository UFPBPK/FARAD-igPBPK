##---------------------------------------------------------------------------------------------------
# This function acquire the partition coefficients through several function
# The code was modfied from the Metrum research group: DOI: https://doi.org/10.1124/dmd.120.090498
# Author: Wei-Chun Chou
# Adivisor: Zhoumeng Lin
# Date: 2020/10/26
##----------------------------------------------------------------------------------------------------

## Loading required packages
library(dplyr)

## Abbreviations: logP: the log-scale n-octanol:buffer partition coefficient;--- 
## pka: acid disscoiation constant; 
## BP is the in vitro blood:plasma ratio;--
## fup: fraction unbound in plasma;---
## dat: input the tissue composition data;---
## type: chemical type (#1-neutral; #2-monoprotic acid; #3-monoprotic base; #4-diprotic acid; #5-diprotic base--
## #6-monoprotic acid monoprotic base (acid comes first); #7-triprotic acid; #8-triprotic base; ---
## #9-diprotic acid monoprotic base (first two are acid);  #10-diprotic base monoprotic acid (first one is acid))

PCcoef <- function(logP, pKa, fup, BP, type = 2, pred="P&T", dat){ ## deafult type  = 2 (acid); Poulin and Theil method 
    if(pred=="P&T"){
        source("Kpfun_PT.R")
        pcoeff <- Kp_PT(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, dat=dat)
    }else if(pred=="Berez"){  #Berezhkovskiy 
        source("Kpfun_Berez.R")
        pcoeff <- Kp_Berez(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, dat=dat)
    }else if(pred == "pksim"){  #standard PK-Sim, Willmann et al. 2008
        source("Kpfun_Pksim.R")  
        pcoeff <- Kp_Pksim(logP=logP, fup=fup, dat=dat)
    }else if(pred == "Schmitt"){  #Schmitt, Walter 2008
        source("Kpfun_Schmitt.R")
        pcoeff <- Kp_Schmit(logP=logP, pKa=pKa, fup=fup, type=type, dat=dat)
    }else if(pred == "RR"){#Rodgers and Rowland 2006
        source("Kpfun_RR.R")  
        pcoeff <- Kp_RR(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, dat=dat)
    }else {
        source("Kpfun_Qppr.R")  
        pcoeff <- Kp_Qppr(logP=logP, fup=fup, dat=dat)
    }
    return(pcoeff)
}

PCcoef_all <- function(logP, pKa, fup, BP, type = 2, dat){ ## deafult type  = 2 (acid); Poulin and Theil method 
        source("Kpfun_PT.R")        
        source("Kpfun_Berez.R")        
        source("Kpfun_Pksim.R")         
        source("Kpfun_Schmitt.R")        
        source("Kpfun_RR.R")          
        source("Kpfun_Qppr.R")  

        pcoeff_PT <- Kp_PT(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, dat=dat)
        pcoeff_Br <- Kp_Berez(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, dat=dat)
        pcoeff_PKsim <- Kp_Pksim(logP=logP, fup=fup, dat=dat)
        pcoeff_Schmit <- Kp_Schmit(logP=logP, pKa=pKa, fup=fup, type=type, dat=dat)
        pcoeff_RR <- Kp_RR(logP=logP, pKa=pKa, fup=fup, BP=BP, type=type, dat=dat)
        pcoeff_Qppr <- Kp_Qppr(logP=logP, fup=fup, dat=dat)
      
    return(list(pcoeff_PT= pcoeff_PT, 
                pcoeff_Br = pcoeff_Br,
                pcoeff_RR = pcoeff_RR,
                pcoeff_PKsim=pcoeff_PKsim, 
                pcoeff_Schmit=pcoeff_Schmit))
}
