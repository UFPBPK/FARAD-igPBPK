##---------------------------------------------------------------------------------------------------
# This function generates a list of the partition coefficients based on the PKSim method
# The code was modfied from the Metrum research group: DOI: https://doi.org/10.1124/dmd.120.090498
# Author: Wei-Chun Chou
# Adivisor: Zhoumeng Lin
# Date: 2020/10/26
##----------------------------------------------------------------------------------------------------


## Loading required packages
library(dplyr)

## The function calcuate the tissue:plasma partition coefficients based on the method---
## of the PKSim: DOI: https://doi.org/10.1517/17425255.1.1.159

## Abbreviations: 
## logP: the log-scale n-octanol:buffer partition coefficient;--- 
## pka: acid disscoiation constant; BP is the in vitro blood:plasma ratio;--
## dat: input the tissue composition data;---

Kp_Pksim <- function(logP, fup, dat){
    
    ## Exctract the tissue composition ---
    ## by filtering the data to exclude the "Plasma", "Red blood cells (RBCs)"
    dat_all <- dat %>% filter (!tissue %in% c("Plasma", "RBCs"))
    
    # logMA is the log of membrane affinity = phosphatidylcholine:water (neutral phospholipid:water) partition coefficient;
    # The estimation of logMA was derved from the Yun & Edginton (2013); difference with Schmitt, 2009
    
    logMA <- 1.24 + 0.304*logP  #in case we don't have a direct logMA; equation from Yun & Edginton (2013)
    K_n_pl <- 10^logMA    #neutral phospholipids:water partition coefficient
    
    # K_protein: protein:water partition; Schmitt, 2008
    # K_protein <- 0.163+0.0221*K_n_pl from schmitt, 2008; 
    K_protein <- ((0.81 + 0.11 * K_n_pl)/24.92)*5 # Here we used the equaltion lised in Utsey et al., 2020; From PK-Sim (very similar value to the other method)
    
    kp <- (dat_all$f_water + (K_n_pl*dat_all$f_lipids) + (K_protein*dat_all$f_proteins))*fup
    #denom <- 0.945 + (10^logMA*0.00575) + (0.93*fup)  #plasma fractions
    #kp <- kp/denom  #according to Willmann et al. (2005)
    dat2 <- data.frame(tissue=dat_all$tissue, Kp=kp)
    name <- dat2$tissue %>% substr(1,2) %>% tolower()
    name <- paste("Kp", name, sep ="_")
    Kp <- as.list(dat2$Kp)
    names(Kp) <- name
    
    return(Kp)
}
    
    
    
    
    