#------------------------------------------------------------------------------
# The code is to reproduce the figure 2
#------------------------------------------------------------------------------

## Loading R packages
library(ggplot2)
library(dplyr)
library(ggprism) ## This package provided you the ?GraphPad Prism? like plot
library(patchwork) ## This package is used to add plots together


## Loading RDS files; These RDS files produced from model calibration
Ccattle_FLO  <- readRDS(file = 'Ccattle_FLO.rds')
CSwine_FLO   <- readRDS(file = 'CSwine_FLO.rds')
Ccattle_FLU  <- readRDS(file = 'Ccattle_FLU.rds')
CSwine_FLU   <- readRDS(file = 'CSwine_FLU.rds')
Ccattle_PG   <- readRDS(file = 'Ccattle_PG.rds')
CSwine_PG    <- readRDS(file = 'CSwine_PG.rds')

## Loading RDS files; These RDS files produced from model evaluation
Ccattle_FLO_eva  <- readRDS(file = 'Cattle_FLO_eva.rds')
CSwine_FLO_eva   <- readRDS(file = 'CSwine_FLO_eva.rds')
Ccattle_FLU_eva  <- readRDS(file = 'Cattle_FLU_eva.rds')
CSwine_FLU_eva   <- readRDS(file = 'CSwine_FLU_eva.rds')
Ccattle_PG_eva   <- readRDS(file = 'Cattle_PG_eva.rds')
CSwine_PG_eva    <- readRDS(file = 'CSwine_PG_eva.rds')


## Make the table for plot the figure for 'Florfenicol'
## plot data of cattle for Florfenicol (FLO)
PDat_C_FLO <- cbind.data.frame (OBS = Ccattle_FLO$residuals$obs, 
                                PRE = Ccattle_FLO$residuals$mod,
                                RES = Ccattle_FLO$residuals$res,
                                Species = "Cattle")

## Model evaluation 
PDat_C_FLO_eva <- cbind.data.frame (OBS = Ccattle_FLO_eva$residuals$obs, 
                                    PRE = Ccattle_FLO_eva$residuals$mod,
                                    RES = Ccattle_FLO_eva$residuals$res,
                                    Species = "Cattle")

## plot data of swine for Florfenicol (FLO)
PDat_S_FLO <- cbind.data.frame (OBS = CSwine_FLO$residuals$obs, 
                                PRE = CSwine_FLO$residuals$mod,
                                RES = CSwine_FLO$residuals$res,
                                Species = "Swine")
## Model evaluation 
PDat_S_FLO_eva <- cbind.data.frame (OBS = CSwine_FLO_eva$residuals$obs, 
                                PRE = CSwine_FLO_eva$residuals$mod,
                                RES = CSwine_FLO_eva$residuals$res,
                                Species = "Swine")

## plot data of cattle for Penicillin G (PG)
PDat_C_PG <- cbind.data.frame (OBS  = Ccattle_PG$residuals$obs, 
                               PRE = Ccattle_PG$residuals$mod,
                               RES = Ccattle_PG$residuals$res,
                               Species = "Cattle")

## Model evaluation 
PDat_C_PG_eva <- cbind.data.frame (OBS  = Ccattle_PG_eva$residuals$obs, 
                               PRE = Ccattle_PG_eva$residuals$mod,
                               RES = Ccattle_PG_eva$residuals$res,
                               Species = "Cattle")

## plot data of swine for Penicillin G  (PG)
PDat_S_PG <- cbind.data.frame (OBS = CSwine_PG$residuals$obs, 
                               PRE = CSwine_PG$residuals$mod,
                               RES = CSwine_PG$residuals$res,
                               Species = "Swine")

## Model evaluation 
PDat_S_PG_eva <- cbind.data.frame (OBS = CSwine_PG_eva$residuals$obs, 
                               PRE = CSwine_PG_eva$residuals$mod,
                               RES = CSwine_PG_eva$residuals$res,
                               Species = "Swine")

## plot data of cattle for Flunxin (FLU)
PDat_C_FLU <- cbind.data.frame (OBS  = Ccattle_FLU$residuals$obs, 
                                PRE = Ccattle_FLU$residuals$mod,
                                RES = Ccattle_FLU$residuals$res,
                                Species = "Cattle")

## Model evaluation 
PDat_C_FLU_eva <- cbind.data.frame (OBS  = Ccattle_FLU_eva$residuals$obs, 
                                PRE = Ccattle_FLU_eva$residuals$mod,
                                RES = Ccattle_FLU_eva$residuals$res,
                                Species = "Cattle")

## plot data of swine for Flunxin (FLU)
PDat_S_FLU <- cbind.data.frame (OBS = CSwine_FLU$residuals$obs, 
                                PRE = CSwine_FLU$residuals$mod,
                                RES = CSwine_FLU$residuals$res,
                                Species = "Swine")

## Model evaluation 
PDat_S_FLU_eva <- cbind.data.frame (OBS = CSwine_FLU_eva$residuals$obs, 
                                PRE = CSwine_FLU_eva$residuals$mod,
                                RES = CSwine_FLU_eva$residuals$res,
                                Species = "Swine")

## combine the data of cattle and swine
## FLO
PDat_FLO <- rbind(PDat_C_FLO, PDat_S_FLO)
PDat_FLO <- PDat_FLO %>% mutate (Log.OBS = log(OBS,10), Log.PRE = log(PRE,10))

## Model evaluation 
PDat_FLO_eva <- rbind(PDat_C_FLO_eva, PDat_S_FLO_eva)
PDat_FLO_eva <- PDat_FLO_eva %>% mutate (Log.OBS = log(OBS,10), Log.PRE = log(PRE,10))

## PG
PDat_PG <- rbind(PDat_C_PG, PDat_S_PG)
PDat_PG <- PDat_PG %>% mutate (Log.OBS = log(OBS,10), Log.PRE = log(PRE,10))


## Model evaluation 
PDat_PG_eva  <- rbind(PDat_C_PG_eva , PDat_S_PG_eva)
PDat_PG_eva  <- PDat_PG_eva  %>% mutate (Log.OBS = log(OBS,10), Log.PRE = log(PRE,10))

## FLU
PDat_FLU<- rbind(PDat_C_FLU, PDat_S_FLU)
PDat_FLU <- PDat_FLU %>% mutate (Log.OBS = log(OBS,10), Log.PRE = log(PRE,10))

## Model evaluation 
PDat_FLU_eva  <- rbind(PDat_C_FLU_eva , PDat_S_FLU_eva )
PDat_FLU_eva  <- PDat_FLU_eva  %>% mutate (Log.OBS = log(OBS,10), Log.PRE = log(PRE,10))


## Estimate the R-squared using linear regression
fit_FLO <- lm(Log.PRE ~ Log.OBS, data = PDat_FLO)
summary(fit_FLO)

fit_PG <- lm(Log.PRE ~ Log.OBS, data = PDat_PG)
summary(fit_PG)

fit_FLU <- lm(Log.PRE ~ Log.OBS , data = PDat_FLU)
summary(fit_FLU)


fit_FLO_eva <- lm(Log.PRE ~ Log.OBS, data = PDat_FLO_eva)
summary(fit_FLO_eva)

fit_PG_eva <- lm(Log.PRE ~ Log.OBS, data = PDat_PG_eva)
summary(fit_PG_eva)

fit_FLU_eva <- lm(Log.PRE ~ Log.OBS , data = PDat_FLU_eva)
summary(fit_FLU_eva)

## Estimate the MAPE
# z <- (sum(PDat_FLO$OBS/PDat_FLO$PRE>2))/(dim(PDat_FLO)[1])
# z_eva <- (sum(PDat_FLO_eva$OBS/PDat_FLO_eva$PRE>2))/(dim(PDat_FLO_eva)[1])
# 
# MAPE_FLO <- mean(PDat_FLU$OBS/PDat_FLU$PRE)
# MAPE_FLO <- mean(PDat_PG$OBS/PDat_PG$PRE)
# 


## Add the observed-to-prediction ratios column
PlotDat_FLO <- PDat_FLO %>% mutate(prediction = predict(fit_FLO), OPR = PRE/OBS)
PlotDat_PG  <- PDat_PG %>% mutate(prediction = predict(fit_PG), OPR = PRE/OBS)
PlotDat_FLU <- PDat_FLU %>% mutate(prediction = predict(fit_FLU), OPR = PRE/OBS)

PlotDat_FLO_eva  <- PDat_FLO_eva  %>% mutate(prediction = predict(fit_FLO_eva ), OPR = PRE/OBS)
PlotDat_PG_eva   <- PDat_PG_eva  %>% mutate(prediction = predict(fit_PG_eva ), OPR = PRE/OBS)
PlotDat_FLU_eva  <- PDat_FLU_eva  %>% mutate(prediction = predict(fit_FLU_eva), OPR = PRE/OBS)

## Plot theme
ptheme<-theme (
        plot.background         = element_rect (fill="White"),
        text                    = element_text (family = "Times"),   # text front (Time new roman)
        panel.border            = element_rect (colour = "black", fill=NA, size=2),
        panel.background        = element_rect (fill="White"),
        panel.grid.major        = element_blank(),
        panel.grid.minor        = element_blank(),
        axis.text               = element_text (size   = 15, colour = "black", face = "bold"),    # tick labels along axes
        axis.title              = element_text (size   = 18, colour = "black", face = "bold"),   # label of axes
        legend.position='none') 


## GGplot
## FLO

## Plot for model calibration of FLO
p_FLO <- 
    ggplot(PlotDat_FLO, aes(Log.PRE, Log.OBS)) + 
    geom_point  (aes(colour = factor(Species),  shape = factor(Species)), size = 4)+
    scale_color_manual(values = c('#999999', '#999999')) +
    #scale_alpha_manual(values = c(.7, 0.5))+
    geom_abline (intercept = 0, 
                 slope     = 1,
                 color     ="black", size = 1, alpha = 0.8, linetype = 2) +
    annotation_logticks() +
    scale_y_continuous(limits = c(-3,2), labels = scales::math_format(10^.x))+
    scale_x_continuous(limits = c(-3,2),labels = scales::math_format(10^.x)) +
    ptheme + labs (x = "", y = "")


## Plot for model evaluation of FLO
p_FLO_eva <- 
    ggplot(PlotDat_FLO_eva, aes(Log.PRE, Log.OBS)) + 
    geom_point  (aes(colour = factor(Species),  shape = factor(Species)), size = 4)+
    scale_color_manual(values = c('#999999', '#999999')) +
    #scale_alpha_manual(values = c(.7, 0.5))+
    geom_abline (intercept = 0, 
                 slope     = 1,
                 color     ="black", size = 1, alpha = 0.8, linetype = 2) +
    annotation_logticks() +
    scale_y_continuous(limits = c(-3,2), labels = scales::math_format(10^.x))+
    scale_x_continuous(limits = c(-3,2),labels = scales::math_format(10^.x)) +
    ptheme + labs (x = "", y = "")



## PG
## Plot for model calibration of PG
p_PG <- 
    ggplot(PlotDat_PG, aes(Log.PRE, Log.OBS)) + 
    geom_point  (aes(colour = factor(Species),  shape = factor(Species)), size = 4)+
    scale_color_manual(values = c('blue', 'blue')) +
    #scale_alpha_manual(values = c(.8, 0.5))+
    geom_abline (intercept = 0, 
                 slope     = 1,
                 color     ="black", size = 1, alpha = 0.8, linetype = 2) +
    annotation_logticks() +
    scale_y_continuous(limits = c(-3,2), labels = scales::math_format(10^.x))+
    scale_x_continuous(limits = c(-3,2),labels = scales::math_format(10^.x)) +
    ptheme + labs (x = "", y = "")



## Plot for model evaluation of PG
p_PG_eva  <- 
    ggplot(PlotDat_PG_eva , aes(Log.PRE, Log.OBS)) + 
    geom_point  (aes(colour = factor(Species),  shape = factor(Species)), size = 4)+
    scale_color_manual(values = c('blue', 'blue')) +
    #scale_alpha_manual(values = c(.8, 0.5))+
    geom_abline (intercept = 0, 
                 slope     = 1,
                 color     ="black", size = 1, alpha = 0.8, linetype = 2) +
    annotation_logticks() +
    scale_y_continuous(limits = c(-3,2), labels = scales::math_format(10^.x))+
    scale_x_continuous(limits = c(-3,2),labels = scales::math_format(10^.x)) +
    ptheme + labs (x = "", y = "")



## FLU
## Plot for model calibration of FLU
p_FLU <- 
    ggplot(PlotDat_FLU, aes(Log.PRE, Log.OBS)) + 
    geom_point  (aes(colour = factor(Species),  shape = factor(Species)), size = 4)+
    scale_color_manual(values = c('red', 'red')) +
    #scale_alpha_manual(values = c(.8, 0.5))+
    geom_abline (intercept = 0, 
                 slope     = 1,
                 color     ="black", size = 1, alpha = 0.8, linetype = 2) +
    annotation_logticks() +
    scale_y_continuous(limits = c(-3,2), labels = scales::math_format(10^.x))+
    scale_x_continuous(limits = c(-3,2),labels = scales::math_format(10^.x)) +
    ptheme + labs (x = "", y = "")

## Plot for model evaluation of FLU
p_FLU_eva <- 
    ggplot(PlotDat_FLU_eva, aes(Log.PRE, Log.OBS)) + 
    geom_point  (aes(colour = factor(Species),  shape = factor(Species)), size = 4)+
    scale_color_manual(values = c('red', 'red')) +
    #scale_alpha_manual(values = c(.8, 0.5))+
    geom_abline (intercept = 0, 
                 slope     = 1,
                 color     ="black", size = 1, alpha = 0.8, linetype = 2) +
    annotation_logticks() +
    scale_y_continuous(limits = c(-3,2), labels = scales::math_format(10^.x))+
    scale_x_continuous(limits = c(-3,2),labels = scales::math_format(10^.x)) +
    ptheme + labs (x = "", y = "")



## Export figures

ggsave("Fig 2a.tiff",scale = 1.5,
       plot = p_FLU + p_FLO + p_PG ,
       width = 30, height = 10, units = "cm", dpi=320)


ggsave("Fig 2b.tiff",scale = 1.5,
       plot = p_FLU_eva + p_FLO_eva + p_PG_eva,
       width = 30, height = 10, units = "cm", dpi=320)
















