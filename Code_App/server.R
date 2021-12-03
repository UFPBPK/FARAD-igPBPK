
## Loading required r package
library(shiny)          # R package for shiny App
library(shinydashboard) # R package for shiny App UI structure
library(shinyjs)        # R package for hide and show a shiny element
library(mrgsolve)       # R package for PBPK model
library(ggplot2)        # R package for plot output
library(truncnorm)      # R package for Truncated normal distribution
library(EnvStats)       # R Package for Environmental Statistics, Including US EPA Guidance
library(rmarkdown)      # R package for output report
library(dplyr)          # R package for dataframe manupulation
library(knitr)          # R package for output reportsource("function_2.R")
library(ggplot2)
library(extrafont)
library(ggh4x)



#//////////////////////////////////////////////////////////////////////////
server <- function(input, output) {
  animal    = isolate(input$animal)
  drug      = isolate(input$drug)
  dose      = isolate(input$doselevel)
  tdoses    = isolate(input$numdose)  
  tinterval = isolate(input$tinterval)
  nsim      = isolate(input$N)
  route     = isolate(input$route)
  target    = isolate(input$target)
  TOL <- switch(target,
                'Liver' = isolate(input$tolerance_Liver),
                'Muscle'= isolate(input$tolerance_Muscle),
                'Kidney' = isolate(input$tolerance_Muscle))
  
  
  observeEvent(input$reset,{reset("Therapeutic")})
  observeEvent(input$reset1,{reset("All_parameters")})
  
  
  ## Load Model
  mod <- mcode_cache("pbpk", GenricPBPK) #refer to mcode function in mrgsolve user guide 3.1.2 Inline
  
  ##
  
  
  
  ## Set Up Simulation Subjects
  r1 <- reactive({
    # if (input$MRL == "US") 
    # {TOL = isolate(input$tolerance)
    # end_time = isolate(input$simu_time)}
    # if (input$MRL == "EU") 
    # {TOL = isolate(input$tolerance)
    #end_time = isolate(input$simu_time_EU)}
    
    
    # Dosing, repeated doses
    input$action | input$action1
    
    animal    = isolate(input$animal)
    drug      = isolate(input$drug)
    target    = isolate(input$target)
    N         = isolate(input$N) 

    if (animal == 'Cattle') {
      if (drug == 'Flunixin' | drug =='5-OH flunixin') {
        tdoses    = isolate(input$numdose)  
        tinterval = isolate(input$tinterval)
        pars <- Pars_C_FLU
        
        TOL <- switch(target,
                      'Liver' = isolate(input$tolerance_Liver),
                      'Muscle'= isolate(input$tolerance_Muscle),
                      'Kidney' = isolate(input$tolerance_Muscle))
        
      } else if (drug == 'Florfenicol' | drug =='Florfenicol amine') {
        tdoses    = isolate(input$numdose_3)  
        tinterval = isolate(input$tinterval_3)
        route     = isolate(input$route_3)
        
        pars <- Pars_C_FLO
        TOL <- switch(target,
                      'Liver' = isolate(input$tolerance_Liver_3),
                      'Muscle'= isolate(input$tolerance_Muscle_3),
                      'Kidney' = isolate(input$tolerance_Muscle_3))
        
      } else {
        tdoses    = isolate(input$numdose_5)  
        tinterval = isolate(input$tinterval_5)
        route     = isolate(input$route_5)
        
        pars <- Pars_C_PG
        TOL <- switch(target,
                      'Liver' = isolate(input$tolerance_Liver_5),
                      'Muscle'= isolate(input$tolerance_Muscle_5),
                      'Kidney' = isolate(input$tolerance_Muscle_5))
      }
    }
    
    if (animal == 'Swine') {
      if (drug == 'Flunixin' | drug =='5-OH flunixin') {
        tdoses    = isolate(input$numdose_2)  
        tinterval = isolate(input$tinterval_2)
        route     = isolate(input$route_2)
        
        pars <- Pars_S_FLU
        TOL <- switch(target,
                      'Liver' = isolate(input$tolerance_Liver_2),
                      'Muscle'= isolate(input$tolerance_Muscle_2),
                      'Kidney' = isolate(input$tolerance_Muscle_2))
        
      } else if (drug == 'Florfenicol' | drug =='Florfenicol amine') {
        tdoses    = isolate(input$numdose_4)  
        tinterval = isolate(input$tinterval_4)
        route     = isolate(input$route_4)
        
        pars <- Pars_S_FLO
        TOL <- switch(target,
                      'Liver' = isolate(input$tolerance_Liver_4),
                      'Muscle'= isolate(input$tolerance_Muscle_4),
                      'Kidney' = isolate(input$tolerance_Muscle_4))
        
      } else {
        tdoses    = isolate(input$numdose_6)  
        tinterval = isolate(input$tinterval_6)
        route     = isolate(input$route_6)
        
        pars <- Pars_S_PG
        TOL <- switch(target,
                      'Liver' = isolate(input$tolerance_Liver_6),
                      'Muscle'= isolate(input$tolerance_Muscle_6),
                      'Kidney' = isolate(input$tolerance_Muscle_6))
      }
    }
    
    #{Parameters for Various Exposure Scenarios}
    
    
        
    pars2 <- c(
      BW          = input$BW,       
      QCC         = input$QCC,    
      QLCa        = input$QLC,	                        
      QKCa        = input$QKC,	                    
      QMCa        = input$QMC,	                     
      QFCa        = input$QFC,                 
      QRestCa     = input$QRestC,                
      VLCa        = input$VLC,                                  
      VKCa        = input$VKC,                                   
      VFCa        = input$VFC,                                     
      VMCa        = input$VMC,                                    
      VbloodCa    = input$VbC,    
      VRestCa     = input$VRestC,    
      Fracsc      = input$Fracsc,
      Fracim      = input$Fracim,
      fR          = input$fR,
      Kabs        = input$Kabs,
      Kunabs      = input$Kunabs,                       
      Kim         = input$Kim,                               
      Ksc         = input$Ksc,                              
      Kdissim     = input$Kdissim,     
      Kdisssc     = input$Kdisssc,    
      KehcC       = input$KehcC,
      KmetC       = input$KmetC,                          
      KurineC     = input$KurineC,              
      KbileC      = input$KbileC,
      PL          = input$PL,
      PM          = input$PM,
      PF          = input$PF,
      PK          = input$PK,
      PRest       = input$PRest           
    )
    
    pars [names(pars2)] <- pars2
    
    #N = nsim # 1000, Number of iterations in the Monte Carlo simulation
    idata<-MCsim(N=input$N, pars)
    
    set.seed(830) 
    PDat<-MC.pred (pars = pars, drug = drug, idata = idata,
                   tinterval = tinterval, dose = dose, dtimes = tdoses, route = route)
    
    
    
    return(PDat)
  })
  
  
  
  
  
  ## Parent compound plot
  r2 <-  reactive({
    
    input$action | input$action1
    animal    = isolate(input$animal)
    drug      = isolate(input$drug)
    target    = isolate(input$target)
    data      = r1()

    if (animal == 'Cattle') {
      if (drug == 'Flunixin' | drug =='5-OH flunixin') {
        tdoses    = isolate(input$numdose)  
        tinterval = isolate(input$tinterval)
        
        TOL <- switch(target,
                    'Liver' = isolate(input$tolerance_Liver),
                    'Muscle'= isolate(input$tolerance_Muscle),
                    'Kidney' = isolate(input$tolerance_Muscle),
                    'Plasma' = isolate(input$tolerance_Muscle))
      
      } else if (drug == 'Florfenicol' | drug =='Florfenicol amine') {
        tdoses    = isolate(input$numdose_3)  
        tinterval = isolate(input$tinterval_3)
        
        
        TOL <- switch(target,
                      'Liver' = isolate(input$tolerance_Liver_3),
                      'Muscle'= isolate(input$tolerance_Muscle_3),
                      'Kidney' = isolate(input$tolerance_Muscle_3),
                      'Plasma' = isolate(input$tolerance_Muscle_3))
        
      } else {
        tdoses    = isolate(input$numdose_5)  
        tinterval = isolate(input$tinterval_5)
        
        TOL <- switch(target,
                      'Liver' = isolate(input$tolerance_Liver_5),
                      'Muscle'= isolate(input$tolerance_Muscle_5),
                      'Kidney' = isolate(input$tolerance_Muscle_5),
                      'Plasma' = isolate(input$tolerance_Muscle_5))
      }
    }
    
    if (animal == 'Swine') {
      if (drug == 'Flunixin' | drug =='5-OH flunixin') {
        tdoses    = isolate(input$numdose_2)  
        tinterval = isolate(input$tinterval_2)
        

        TOL <- switch(target,
                      'Liver' = isolate(input$tolerance_Liver_2),
                      'Muscle'= isolate(input$tolerance_Muscle_2),
                      'Kidney' = isolate(input$tolerance_Muscle_2),
                      'Plasma' = isolate(input$tolerance_Muscle_2))
        
      } else if (drug == 'Florfenicol' | drug =='Florfenicol amine') {
        tdoses    = isolate(input$numdose_4)  
        tinterval = isolate(input$tinterval_4)
        
        
        TOL <- switch(target,
                      'Liver' = isolate(input$tolerance_Liver_4),
                      'Muscle'= isolate(input$tolerance_Muscle_4),
                      'Kidney' = isolate(input$tolerance_Muscle_4),
                      'Plasma' = isolate(input$tolerance_Muscle_4))
        
      } else {
        tdoses    = isolate(input$numdose_6)  
        tinterval = isolate(input$tinterval_6)
        
        
        TOL <- switch(target,
                      'Liver' = isolate(input$tolerance_Liver_6),
                      'Muscle'= isolate(input$tolerance_Muscle_6),
                      'Kidney' = isolate(input$tolerance_Muscle_6),
                      'Plasma' = isolate(input$tolerance_Muscle_6))
      }
    }
    
    
    
    target1 <- switch(drug,
                     "Flunixin"  = target,
                     "Florfenicol" = target,
                     "Penicillin G" = target,
                     "5-OH flunixin" = paste0(target,"1"),
                     "Florfenicol amine"= paste0(target,"1"))
  
    p <- switch(target,
                Liver  = wdplot(data,target=target1),
                Plasma = wdplot(data,target=target1),
                Kidney = wdplot(data,target=target1),
                Muscle = wdplot(data,target=target1),
                Fat    = wdplot(data,target=target1))
   
    
    end_time    <- data %>% filter_(paste0(target1,'.ub','<=',TOL, ' & Time >(tdoses*tinterval/24-1)'))%>%select(Time)%>%min()
    endtime     <-end_time + tdoses*tinterval/24
    x.intercept = (data%>%filter_(paste0(target1,'.ub','<=',TOL, ' & Time >(tdoses*tinterval/24-1)'))%>%select(Time)%>%min())
    
    # end_time <- data %>% filter_(paste0(target1,'.lb','>=',0.001, ' & Time >0'))%>%select(Time)%>%max()
    # endtime <-end_time + tdoses*tinterval/24
    # 
    
    p<-p+
       scale_x_continuous (breaks = pretty(0:endtime),
      #                     minor_breaks = pretty(0:endtime, n=endtime),
                           label  = pretty((0:endtime), n = 5)-(tdoses*tinterval/24-1), 
                           limits = c(NA, endtime), expand = c(0,0), guide = "axis_minor") +       
      #  
      geom_line(aes(y = TOL),color = 'black',size = 0.5, linetype = 'twodash', show.legend = F)+ #+ # add tolerance line
      #scale_y_log10(limits = c(1e-15, NA))+
      scale_y_log10(#breaks = (10^(-3):NA),
      labels = scales::trans_format("log10", scales::math_format(10^.x)),
      limits = c(10^-6, NA))+
      annotation_logticks(short = unit(2,"mm"),
                          mid = unit(3,"mm"),
                          long = unit(4,"mm"),size = 1,
                          sides = "l") +
      geom_vline(aes(xintercept = x.intercept), 
                 size = 1, color = "red", linetype = 2,
                 show.legend = F) +
      labs(x="Time (Day)",y=paste('Concentration in', target, "(mg/kg or ppm)"),face="bold") +
     
      ptheme 
    
    return(p)
  })
  
  
  r3 <-  reactive({
    input$action | input$action1
    animal    = isolate(input$animal)
    drug      = isolate(input$drug)
    target    = isolate(input$target)
    data      = r1()
    
    
    if (animal == 'Cattle') {
      if (drug == 'Flunixin' | drug =='5-OH flunixin') {
        tdoses    = isolate(input$numdose)  
        tinterval = isolate(input$tinterval)
        
        TOL_liver  = isolate(input$tolerance_Liver)
        TOL_Muscle = isolate(input$tolerance_Muscle)
        TOL_Kidney = isolate(input$tolerance_Muscle)
        TOL_Plasma = isolate(input$tolerance_Muscle)
  
      } else if (drug == 'Florfenicol' | drug =='Florfenicol amine') {
        tdoses    = isolate(input$numdose_3)  
        tinterval = isolate(input$tinterval_3)
        
        TOL_liver  = isolate(input$tolerance_Liver_3)
        TOL_Muscle = isolate(input$tolerance_Muscle_3)
        TOL_Kidney = isolate(input$tolerance_Muscle_3)
        TOL_Plasma = isolate(input$tolerance_Muscle_3)
        
       
      } else {
        tdoses    = isolate(input$numdose_5)  
        tinterval = isolate(input$tinterval_5)
        
        TOL_liver  = isolate(input$tolerance_Liver_5)
        TOL_Muscle = isolate(input$tolerance_Muscle_5)
        TOL_Kidney = isolate(input$tolerance_Muscle_5)
        TOL_Plasma = isolate(input$tolerance_Muscle_5)
        
      }
    }
    
    if (animal == 'Swine') {
      if (drug == 'Flunixin' | drug =='5-OH flunixin') {
        tdoses    = isolate(input$numdose_2)  
        tinterval = isolate(input$tinterval_2)
        
        
        TOL_liver  = isolate(input$tolerance_Liver_2)
        TOL_Muscle = isolate(input$tolerance_Muscle_2)
        TOL_Kidney = isolate(input$tolerance_Muscle_2)
        TOL_Plasma = isolate(input$tolerance_Muscle_2)
        
        
      } else if (drug == 'Florfenicol' | drug =='Florfenicol amine') {
        tdoses    = isolate(input$numdose_4)  
        tinterval = isolate(input$tinterval_4)
        
        
        TOL_liver  = isolate(input$tolerance_Liver_4)
        TOL_Muscle = isolate(input$tolerance_Muscle_4)
        TOL_Kidney = isolate(input$tolerance_Muscle_4)
        TOL_Plasma = isolate(input$tolerance_Muscle_4)
        
      } else {
        tdoses    = isolate(input$numdose_6)  
        tinterval = isolate(input$tinterval_6)
        
        
        TOL_liver  = isolate(input$tolerance_Liver_6)
        TOL_Muscle = isolate(input$tolerance_Muscle_6)
        TOL_Kidney = isolate(input$tolerance_Muscle_6)
        TOL_Plasma = isolate(input$tolerance_Muscle_6)
      }
    }
    
  
    
    plot.fun <-function (p, target, TOL) {
      
      end_time    <- data %>% filter_(paste0(target,'.ub','<=',TOL, ' & Time >(tdoses*tinterval/24-1)'))%>%select(Time)%>%min()
      endtime     <- end_time + tdoses*tinterval/24
      x.intercept = (data%>%filter_(paste0(target,'.ub','<=',TOL, ' & Time >(tdoses*tinterval/24-1)'))%>%select(Time)%>%min())
      x.intercept1 = (data%>%select_(paste0(target,'.ub')))%>%max()
      
      p1<-p + 
        scale_x_continuous (    breaks = pretty(0:endtime, n = 5),
        #                         minor_breaks = pretty(0:endtime, n=endtime),
                                 label  = pretty((0:endtime), n = 5)-((tdoses*tinterval/24)-1),
                                limits = c(0, endtime)) +
        geom_line(aes(y = TOL),color = 'black',size = 0.5, linetype = 'twodash', show.legend = F)+ #+ # add tolerance line
        scale_y_log10(#breaks = (10^(-3):NA),
          labels = scales::trans_format("log10", scales::math_format(10^.x)),
          limits = c(10^-6, NA))+
        annotation_logticks(short = unit(2,"mm"),
                            mid = unit(3,"mm"),
                            long = unit(4,"mm"),size = 1,
                            sides = "l") +
        labs(x="Time (Day)",y=paste('Concentration in', target, "(mg/kg or ppm)"),face="bold") +
        
        ptheme 

      
      if (x.intercept1 < TOL) {
      p1<-p1} else { 
        p1<-p1+  
        geom_vline(aes(xintercept = x.intercept), 
                   size = 1, color = "red", linetype = 2,
                   show.legend = F)}
        
      return(p1)
    }
    
    if (drug == 'Florfenicol' | drug =='Flunixin' | drug == 'Penicillin G') {
    
    Plot.liver   <- plot.fun(wdplot(data,target="Liver"), target="Liver",  TOL_liver)
    Plot.kidney  <- plot.fun(wdplot(data,target="Kidney"),target="Kidney", TOL_Kidney)
    Plot.Muscle  <- plot.fun(wdplot(data,target="Muscle"),target="Muscle", TOL_Muscle)
    Plot.Plasma  <- plot.fun(wdplot(data,target="Plasma"),target="Plasma", TOL_Plasma)
    
    } else {
    Plot.liver  <- plot.fun(wdplot(data,target="Liver1"),target="Liver1", TOL_liver)
    Plot.kidney <- plot.fun(wdplot(data,target="Kidney1"),target="Kidney1", TOL_Kidney)
    Plot.Muscle <- plot.fun(wdplot(data,target="Muscle1"),target="Muscle1", TOL_Muscle)
    Plot.Plasma <- plot.fun(wdplot(data,target="Plasma1"),target="Plasma1",TOL_Plasma)
    }
      
    
    return(list(Plot.liver = Plot.liver , 
                Plot.kidney = Plot.kidney,
                Plot.Muscle = Plot.Muscle,
                Plot.Plasma = Plot.Plasma))
  })
  
  ##////////////////////////////////////////////////////////////////////////////
  ## Output 
    
  output$title1 <- renderText({
    input$action | input$action1
    
    target    = isolate(input$target)
    drug      = isolate(input$drug)
    
    return(paste(input$drug, 'concentraitons in liver'))
    
    
  })
  
  output$title2 <- renderText({
    input$action | input$action1
    
    target    = isolate(input$target)
    drug      = isolate(input$drug)
    
    target    = isolate(input$target)
    drug      = isolate(input$drug)
    
    return(paste(drug, 'concentraitons in kidney'))
    
    
  })
  
  output$title3 <- renderText({
    input$action | input$action1
    
    target    = isolate(input$target)
    drug      = isolate(input$drug)
    
    return(paste(input$drug, 'concentraitons in muscle'))
    
    
  }) 
  
  output$title4 <- renderText({
    input$action | input$action1
    
    target    = isolate(input$target)
    drug      = isolate(input$drug)
    
    return(paste(drug, 'concentraitons in plasma'))
    
    
  })
  
  
  
  
  
  output$wdtplot <- renderPlot({
    withProgress(message = 'Creating plot', detail='Please Wait...', value = 10, 
                 expr={input$action | input$action1
                   a=r2()
                   return(a)}) 
    
  }, width = 800)
  
  output$wdplot_liver <- renderPlot({
    withProgress(message = 'Creating plot', detail='Please Wait...', value = 10, 
                 expr={input$action | input$action1
                   a=(r3()[[1]])
                   a })
    
    
    
  })
  
  
  
  output$wdplot_kidney <- renderPlot({
    withProgress(message = 'Creating plot', detail='Please Wait...', value = 10, 
                 expr={input$action | input$action1
                   a=(r3()[[2]])
                   a })
    
    
    
  })
  
  
  output$wdplot_muscle <- renderPlot({
    withProgress(message = 'Creating plot', detail='Please Wait...', value = 10, 
                 expr={input$action | input$action1
                   a=(r3()[[3]])
                   a })
    
    
    
  })
  
  

  output$wdplot_plasma <- renderPlot({
    withProgress(message = 'Creating plot', detail='Please Wait...', value = 10, 
                 expr={input$action | input$action1
                   a=(r3()[[4]])
                   a })
    
    
    
  })
  
  
  
  
  
  
  output$info <- renderText({
    input$action | input$action1

    animal    = isolate(input$animal)
    drug      = isolate(input$drug)
    #dose      = isolate(input$doselevel)
    #tdoses    = isolate(input$numdose)  
    #tinterval = isolate(input$tinterval)
    target    = isolate(input$target)
    #route     = isolate(input$route)
    data      = r1()

    if (animal == 'Cattle') {
      if (drug == 'Flunixin' | drug =='5-OH flunixin') {
        tdoses    = isolate(input$numdose)  
        tinterval = isolate(input$tinterval)
        
        TOL <- switch(target,
                      'Liver' = isolate(input$tolerance_Liver),
                      'Muscle'= isolate(input$tolerance_Muscle),
                      'Kidney' = isolate(input$tolerance_Muscle),
                      'Plasma' = isolate(input$tolerance_Muscle))
        
      } else if (drug == 'Florfenicol' | drug =='Florfenicol amine') {
        tdoses    = isolate(input$numdose_3)  
        tinterval = isolate(input$tinterval_3)
        
        
        TOL <- switch(target,
                      'Liver' = isolate(input$tolerance_Liver_3),
                      'Muscle'= isolate(input$tolerance_Muscle_3),
                      'Kidney' = isolate(input$tolerance_Muscle_3),
                      'Plasma' = isolate(input$tolerance_Muscle_3))
        
      } else {
        tdoses    = isolate(input$numdose_5)  
        tinterval = isolate(input$tinterval_5)
        
        TOL <- switch(target,
                      'Liver' = isolate(input$tolerance_Liver_5),
                      'Muscle'= isolate(input$tolerance_Muscle_5),
                      'Kidney' = isolate(input$tolerance_Muscle_5),
                      'Plasma' = isolate(input$tolerance_Muscle_5))
      }
    }
    
    if (animal == 'Swine') {
      if (drug == 'Flunixin' | drug =='5-OH flunixin') {
        tdoses    = isolate(input$numdose_2)  
        tinterval = isolate(input$tinterval_2)
        
        
        TOL <- switch(target,
                      'Liver' = isolate(input$tolerance_Liver_2),
                      'Muscle'= isolate(input$tolerance_Muscle_2),
                      'Kidney' = isolate(input$tolerance_Muscle_2),
                      'Plasma' = isolate(input$tolerance_Muscle_2))
        
      } else if (drug == 'Florfenicol' | drug =='Florfenicol amine') {
        tdoses    = isolate(input$numdose_4)  
        tinterval = isolate(input$tinterval_4)
        
        
        TOL <- switch(target,
                      'Liver' = isolate(input$tolerance_Liver_4),
                      'Muscle'= isolate(input$tolerance_Muscle_4),
                      'Kidney' = isolate(input$tolerance_Muscle_4),
                      'Plasma' = isolate(input$tolerance_Muscle_4))
        
      } else {
        tdoses    = isolate(input$numdose_6)  
        tinterval = isolate(input$tinterval_6)
        
        
        TOL <- switch(target,
                      'Liver' = isolate(input$tolerance_Liver_6),
                      'Muscle'= isolate(input$tolerance_Muscle_6),
                      'Kidney' = isolate(input$tolerance_Muscle_6),
                      'Plasma' = isolate(input$tolerance_Muscle_6))
      }
    }
    
    
    target1 <- switch(drug,
                     "Flunixin"  = target,
                     "Florfenicol" = target,
                     "Penicillin G" = target,
                     "5-OH flunixin" = paste0(target,"1"),
                     "Florfenicol amine"= paste0(target,"1"))

    a<-round(data %>% filter_(paste0(target1,'.ub','<=',TOL, ' & Time >(tdoses*tinterval/24-1)'))%>%select(Time)%>%min(),2)
    a1=(a-tdoses*tinterval/24)+1
    
    
    b1 <- paste("Figure: Results of PBPK model with Monte Carlo simulation for", drug, "concentrations in", target1, "for", animal)
    b2 <- paste("following dosing regiment:", route, "administration at", dose,"mg/kg for", tdoses, "doses with", tinterval, "intervals.")
    b3 <- paste("Each of the simulations was run for 1000 iterations. The median (black line), 1th and 99th (black dashed lines) percentiles of simulated results were plotted.")
    b4 <- paste("The tolerance is shown on each of panels using the horizontal dashed line", "(Tolearance for liver is", TOL,"µg/g)")
    b5 <- paste("The intersection between the x-axis and the black vertical line indicates the withdraw intervals (WDIs).")

    
    c <- paste("The recommended withdrawal interval (WDI) is",round(a1, 2))
    
    info <-paste(b1, b2,"", "Note:", b3, b4, b5, sep = "\n")   
    
    return(info)
  })
  
  output$info2 <- renderText({
    input$action | input$action1
    
    animal    = isolate(input$animal)
    drug      = isolate(input$drug)
    N         = isolate(input$N)
    dose      = isolate(input$doselevel)
    tdoses    = isolate(input$numdose)  
    tinterval = isolate(input$tinterval)
    target    = isolate(input$target)
    route     = isolate(input$route)
    data      = r1()
    
    
    if (animal == 'Cattle') {
      if (drug == 'Flunixin' | drug =='5-OH flunixin') {
        tdoses    = isolate(input$numdose)  
        tinterval = isolate(input$tinterval)
        
        TOL <- switch(target,
                      'Liver' = isolate(input$tolerance_Liver),
                      'Muscle'= isolate(input$tolerance_Muscle),
                      'Kidney' = isolate(input$tolerance_Muscle),
                      'Plasma' = isolate(input$tolerance_Muscle))
        
      } else if (drug == 'Florfenicol' | drug =='Florfenicol amine') {
        tdoses    = isolate(input$numdose_3)  
        tinterval = isolate(input$tinterval_3)
        
        
        TOL <- switch(target,
                      'Liver' = isolate(input$tolerance_Liver_3),
                      'Muscle'= isolate(input$tolerance_Muscle_3),
                      'Kidney' = isolate(input$tolerance_Muscle_3),
                      'Plasma' = isolate(input$tolerance_Muscle_3))
        
      } else {
        tdoses    = isolate(input$numdose_5)  
        tinterval = isolate(input$tinterval_5)
        
        TOL <- switch(target,
                      'Liver' = isolate(input$tolerance_Liver_5),
                      'Muscle'= isolate(input$tolerance_Muscle_5),
                      'Kidney' = isolate(input$tolerance_Muscle_5),
                      'Plasma' = isolate(input$tolerance_Muscle_5))
      }
    }
    
    if (animal == 'Swine') {
      if (drug == 'Flunixin' | drug =='5-OH flunixin') {
        tdoses    = isolate(input$numdose_2)  
        tinterval = isolate(input$tinterval_2)
        
        
        TOL <- switch(target,
                      'Liver' = isolate(input$tolerance_Liver_2),
                      'Muscle'= isolate(input$tolerance_Muscle_2),
                      'Kidney' = isolate(input$tolerance_Muscle_2),
                      'Plasma' = isolate(input$tolerance_Muscle_2))
        
      } else if (drug == 'Florfenicol' | drug =='Florfenicol amine') {
        tdoses    = isolate(input$numdose_4)  
        tinterval = isolate(input$tinterval_4)
        
        
        TOL <- switch(target,
                      'Liver' = isolate(input$tolerance_Liver_4),
                      'Muscle'= isolate(input$tolerance_Muscle_4),
                      'Kidney' = isolate(input$tolerance_Muscle_4),
                      'Plasma' = isolate(input$tolerance_Muscle_4))
        
      } else {
        tdoses    = isolate(input$numdose_6)  
        tinterval = isolate(input$tinterval_6)
        
        
        TOL <- switch(target,
                      'Liver' = isolate(input$tolerance_Liver_6),
                      'Muscle'= isolate(input$tolerance_Muscle_6),
                      'Kidney' = isolate(input$tolerance_Muscle_6),
                      'Plasma' = isolate(input$tolerance_Muscle_6))
      }
    }
    
    
    target <- switch(drug,
                     "Flunixin"  = target,
                     "Florfenicol" = target,
                     "Penicillin G" = target,
                     "5-OH flunixin" = paste0(target,"1"),
                     "Florfenicol amine"= paste0(target,"1"))
    
    a<-round(data %>% filter_(paste0(target,'.ub','<=',TOL, ' & Time >(tdoses*tinterval/24-1)'))%>%select(Time)%>%min(),2)
    a1=(a-tdoses*tinterval/24)+1
    info <- paste("The recommended withdrawal interval (WDI) is",round(a1, 2))
    
    
    return(info)
    
  }) 
  
  ###/////////////////////////////////////// Dowload figures
  ## Download the plot on the main panel (wdplot)
  
 
  
  output$downloadwdt <- downloadHandler(
    filename = 'concentration.jpg',
    content = function(file) {
      ggsave(file, r2())
    },
    contentType = 'image/jpg'
  )
  
  output$downl <- downloadHandler(
    filename = paste(input$drug, 'concentraitons in liver.jpg'),
    content = function(file) {
      ggsave(file, plot=r3()[[1]])
    },
    contentType = 'image/jpg'
  )
  
  output$downl1 <- downloadHandler(
    filename = "Metabolites concentrations in liver.jpg",
    content = function(file) {
      ggsave(file, r3()[[2]])
    },
    contentType = 'image/jpg'
  )
  
  output$downk <- downloadHandler(
    filename =  paste(input$drug, 'concentraitons in kidney.jpg'),
    content = function(file) {
      ggsave(file, r3()[[3]])
    },
    contentType = 'image/jpg'
  )
  
  output$downk1 <- downloadHandler(
    filename = "Metabolites Concentrations in kidney.jpg",
    content = function(file) {
      ggsave(file, r3()[[4]])
    },
    contentType = 'image/jpg'
  )
  

  output$downm <- downloadHandler(
    filename =  "Parent compound concentrations in muscle.jpg",
    content = function(file) {
      ggsave(file, r3()[[5]])
    },
    contentType = 'image/jpg'
  )
  
  output$downm1 <- downloadHandler(
    filename = "Metabolites concentrations in muscle.jpg",
    content = function(file) {
      ggsave(file, r3()[[6]])
    },
    contentType = 'image/jpg'
  )

  ## Download table------------------------------------------------------------
  output$table1 <- renderTable ({
    datatable <- r1()
  })
  
  output$downloadtable <- downloadHandler(
    filename = 'table.csv',
    content = function(file) {
      write.csv(format(r1(), digits = 4),file)
    }
  )
  
## End-------------------------------------------------------------------------  
  
  
    
## Download Report-------------------------------------------------------------
  output$downloadreport <- downloadHandler(
    filename = function() {
      paste('my-report', sep = '.', switch(
        input$format, PDF = 'pdf', HTML = 'html', Word = 'docx'
      ))
    },

    content = function(file) {

        src <- normalizePath('report3.Rmd')

      # temporarily switch to the temp dir, in case you do not have write
      # permission to the current working directory
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      file.copy(src, 'report3.Rmd', overwrite = TRUE)

      library(rmarkdown)
      out <- render('report3.Rmd', switch(
        input$format,
        PDF = pdf_document(), HTML = html_document(), Word = word_document()
      ))
      file.rename(out, file)
    }
  )
  
 ##---------------------------------------------------------------------------- 
  
  
 ## Report output //////////////////////////////////////////////////////// 
  get_report <- reactive({
    input$action | input$action1
    
    animal    = isolate(input$animal)
    drug      = isolate(input$drug)
    target    = isolate(input$target)
    N         = isolate(input$N)
    data_tissue <- r1()
    
    if (animal == 'Cattle') {
      if (drug == 'Flunixin' | drug =='5-OH flunixin') {
        tdoses    = isolate(input$numdose)  
        tinterval = isolate(input$tinterval)
        dose      = isolate(input$doselevel)
        route     = isolate(input$route)
        
        TOL_liver  = isolate(input$tolerance_Liver)
        TOL_Muscle = isolate(input$tolerance_Muscle)
        TOL_Kidney = isolate(input$tolerance_Muscle)
        TOL_Plasma = isolate(input$tolerance_Muscle)
        
        
      } else if (drug == 'Florfenicol' | drug =='Florfenicol amine') {
        tdoses    = isolate(input$numdose_3)  
        tinterval = isolate(input$tinterval_3)
        dose      = isolate(input$doselevel_3)
        route     = isolate(input$route_3)
        
        
        TOL_liver  = isolate(input$tolerance_Liver_3)
        TOL_Muscle = isolate(input$tolerance_Muscle_3)
        TOL_Kidney = isolate(input$tolerance_Muscle_3)
        TOL_Plasma = isolate(input$tolerance_Muscle_3)
        
      } else {
        tdoses    = isolate(input$numdose_5)  
        tinterval = isolate(input$tinterval_5)
        dose      = isolate(input$doselevel_5)
        route     = isolate(input$route_5)
        
        TOL_liver  = isolate(input$tolerance_Liver_5)
        TOL_Muscle = isolate(input$tolerance_Muscle_5)
        TOL_Kidney = isolate(input$tolerance_Muscle_5)
        TOL_Plasma = isolate(input$tolerance_Muscle_5)
        
      }
    }
    
    if (animal == 'Swine') {
      if (drug == 'Flunixin' | drug =='5-OH flunixin') {
        tdoses    = isolate(input$numdose_2)  
        tinterval = isolate(input$tinterval_2)
        dose      = isolate(input$doselevel_2)
        route     = isolate(input$route_2)
        
        TOL_liver  = isolate(input$tolerance_Liver_2)
        TOL_Muscle = isolate(input$tolerance_Muscle_2)
        TOL_Kidney = isolate(input$tolerance_Muscle_2)
        TOL_Plasma = isolate(input$tolerance_Muscle_2)
        
      } else if (drug == 'Florfenicol' | drug =='Florfenicol amine') {
        tdoses    = isolate(input$numdose_4)  
        tinterval = isolate(input$tinterval_4)
        dose      = isolate(input$doselevel_4)
        route     = isolate(input$route_4)
        
        
        TOL_liver  = isolate(input$tolerance_Liver_4)
        TOL_Muscle = isolate(input$tolerance_Muscle_4)
        TOL_Kidney = isolate(input$tolerance_Muscle_4)
        TOL_Plasma = isolate(input$tolerance_Muscle_4)
        
      } else {
        tdoses    = isolate(input$numdose_6)  
        tinterval = isolate(input$tinterval_6)
        dose      = isolate(input$doselevel_6)
        route     = isolate(input$route_6)
        
        
        TOL_liver  = isolate(input$tolerance_Liver_6)
        TOL_Muscle = isolate(input$tolerance_Muscle_6)
        TOL_Kidney = isolate(input$tolerance_Muscle_6)
        TOL_Plasma = isolate(input$tolerance_Muscle_6)
      }
    }
    
    a = round(data_tissue%>%filter_(paste0('Liver','.ub','<=',TOL_liver, ' & Time >(tdoses*tinterval/24-1)'))%>%select(Time)%>%min())
    b = round(data_tissue%>%filter_(paste0('Muscle','.ub','<=',TOL_Muscle, ' & Time >(tdoses*tinterval/24-1)'))%>%select(Time)%>%min())
    c = round(data_tissue%>%filter_(paste0('Kidney','.ub','<=',TOL_Kidney, ' & Time >(tdoses*tinterval/24-1)'))%>%select(Time)%>%min())
    d = round(data_tissue%>%filter_(paste0('Plasma','.ub','<=',TOL_Plasma, ' & Time >(tdoses*tinterval/24-1)'))%>%select(Time)%>%min())
    
    WDI_data_liver   = (a-tdoses*tinterval/24)+1
    WDI_data_muscle  = (b-tdoses*tinterval/24)+1 
    WDI_data_kidney  = (c-tdoses*tinterval/24)+1 
    WDI_data_plasma  = (c-tdoses*tinterval/24)+1 
    
    

    
    return(list(animal, drug, target, N,dose,route,tdoses,tinterval,
                TOL_liver, TOL_Muscle, TOL_Kidney, TOL_Plasma, 
                WDI_data_liver, WDI_data_muscle, WDI_data_kidney, WDI_data_plasma))
  })
  
  
## End of report output //////////////////////////////////////////////////////// 
  
  
  
  
  
  
  
  
}