##---------------------------------------------------------------------------------------------------
# This function generates a list of the partition coefficients based on the DeJongh method
# The code based on DeJongh et al., Arch Toxicol (1997): DOI: 10.1007/s002040050463
# Author: Wei-Chun Chou
# Adivisor: Zhoumeng Lin
# Date: 2020/10/26
##----------------------------------------------------------------------------------------------------

## Loading required packages
library(dplyr)

## Abbreviations: 
## logP: the log-scale n-octanol:buffer partition coefficient;--- 
## pka: acid disscoiation constant; BP is the in vitro blood:plasma ratio;--
## dat: input the tissue composition data;---

Kp_Qppr <- function(logP, fup, dat){
    
    ## Exctract the tissue composition ---
    ## by filtering the data to exclude the "Plasma", "Red blood cells (RBCs)"
    dat_all <- dat %>% filter (!tissue %in% c("Plasma", "RBCs"))
    dat_b   <- dat %>% filter (tissue ==c("Plasma"))
        
    ## Water (Vw) and lipid (Vl) fractions of tissues in rats. Values are normalized to Vw + Vl=1
    Vw_all = (dat_all$f_water)/(dat_all$f_water + dat_all$f_lipids)
    Vl_all = 1 - Vw_all
    Vw_b   = (dat_b$f_water)/(dat_b$f_water + dat_b$f_lipids)
    Vl_b   = 1 - Vw_b
    
    ## All pararmeter values obtained from the Table 2 in DeJongh et al. (1997)
    ## A and B is collander-type coefficients
    
    A <- c(0.7, NA, 0.48, NA, 0.57, NA, 0.44, NA, 0.29, NA, NA)
    B <- c(-0.02, NA, -0.21, NA, -0.19, NA, -0.19, NA, -0.55, NA, NA)
    
    ## Klw_t is the lipid-water partition coefficient of a compound within a specific tissue
    ## Klw_b is the lipid-water partition coefficient of a compound within blood
    Klw_b = Klw_t = (10^logP)^A
    
    ## Calculated the partition coefficient (PC) of tissues by QPPR model
    Kp <- (Vl_all[which(!is.na(A))]*Klw_t[which(!is.na(A))] + Vw_all[which(!is.na(A))])/(Vl_b*Klw_b[which(!is.na(A))] + Vw_b) + B[which(!is.na(B))]
    name <- dat_all$tissue[which(!is.na(A))] %>% substr(1,2) %>% tolower()
    name <- paste("Kp", name, sep ="_")
    names(Kp) <- name
    
    return(Kp)
    
}
    
    
    