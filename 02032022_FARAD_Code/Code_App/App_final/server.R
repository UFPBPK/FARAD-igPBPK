server <- function(input, output) {
  
  
  # DEFINE SETS-----------------------------------------------------------------
  
  observeEvent(input$reset,{reset("Therapeutic")})
  observeEvent(input$reset1,{reset("All_parameters")})
  
  
  ## Load mrgsolve pbpk Model
  mod <- mcode_cache("pbpk", GenricPBPK) 
  
  
  ## Set Up Simulation Subjects
  set_select <- reactive({
    
    input$action | input$action1
    
    isolate({
      animal = input$animal
      drug   = input$drug
      N      = input$N
      
      if (animal == 'Cattle') {
        
        if (drug == 'Flunixin') {
          tdoses    = input$ndose_BC_FLU  
          tinterval = input$int_BC_FLU
          target    = input$tar_BC_FLU
          route     = input$ro_BC_FLU
          dose      = input$dose_BC_FLU
          pars      = Pars_C_FLU
          
          
        } else if (drug == 'Florfenicol') {
          tdoses    = input$ndose_BC_FLO  
          tinterval = input$int_BC_FLO
          target    = input$tar_BC_FLO
          route     = input$ro_BC_FLO
          dose      = input$dose_BC_FLO
          pars      = Pars_C_FLO
          
          
        } else {
          tdoses    = input$ndose_BC_PG  
          tinterval = input$int_BC_PG 
          target    = input$tar_BC_PG 
          route     = input$ro_BC_PG 
          dose      = input$dose_BC_PG
          pars      = Pars_C_PG
          
          
          
        }
      }
      
      
      if (animal == 'Swine') {
        if (drug == 'Flunixin') {
          tdoses    = input$ndose_SW_FLU  
          tinterval = input$int_SW_FLU
          target    = input$tar_SW_FLU
          route     = input$ro_SW_FLU
          dose      = input$dose_SW_FLU
          pars      = Pars_S_FLU
          
          
          
        } else if (drug == 'Florfenicol') {
          tdoses    = input$ndose_SW_FLO  
          tinterval = input$int_SW_FLO
          target    = input$tar_SW_FLO
          route     = input$ro_SW_FLO
          dose      = input$dose_SW_FLO
          pars      = Pars_S_FLO
          
        } else {
          tdoses    = input$ndose_SW_PG  
          tinterval = input$int_SW_PG 
          target    = input$tar_SW_PG 
          route     = input$ro_SW_PG 
          dose      = input$dose_SW_PG
          pars      = Pars_S_PG
          
        }
      }
      
      if (animal == 'Cattle' & drug == 'Flunixin') {
        TOL = c (0.125, 0.125, 0.025)
      } else if (animal == 'Cattle' & drug == 'Florfenicol') {
        TOL = c (3.7, 3.7, 0.3)
      } else if (animal == 'Cattle' & drug == 'Penicillin G') {
        TOL = c (0.05, 0.05, 0.05)
      } else if (animal == 'Swine' & drug == 'Flunixin') {
        TOL = c (0.03, 0.03, 0.025)
      } else if (animal == 'Swine' & drug == 'Florfenicol') {
        TOL = c (2.5, 2.5, 0.2)
      } else if (animal == 'Swine' & drug == 'Penicillin G') {
        TOL = c (0.025, 0.025, 0.025)
      }
      
      TOL <- switch(target,
                    'Liver'  = TOL[1],
                    'Muscle' = TOL[3],
                    'Kidney' = TOL[2],
                    'Plasma' = TOL[1])
      
    })   
    return(list(animal = animal, drug = drug, tdoses = tdoses, tinterval =tinterval, 
                target = target, route = route, N = N, 
                TOL = TOL, dose = dose, pars=pars))
    
    
  })
  
  
  data_1 <- reactive({
    
    # input$action | input$action1
    # isolate({drug = input$drug})
    
    setSelect <- set_select()
    
    pars      <- setSelect$pars
    drug      <- setSelect$drug
    tdoses    <- setSelect$tdoses
    tinterval <- setSelect$tinterval
    target    <- setSelect$target
    route     <- setSelect$route
    dose      <- setSelect$dose
    N         <- setSelect$N
    
    #{Parameters for Various Exposure Scenarios}
    
    
    #N = nsim # 1000, Number of iterations in the Monte Carlo simulation
    set.seed(123) 
    idata<-MCsim(N=N, pars)
    
    
    PDat<-MC.pred (N=N, pars = pars, drug = drug, idata = idata,
                   tinterval = tinterval, dose = dose, dtimes = tdoses, route = route)
    
    
    return(PDat)
  })
  
  
  data_2 <- reactive({
    
    input$action | input$action1
    isolate({drug = input$drug})
    
    setSelect <- set_select()
    
    pars      <- setSelect$pars
    
    tdoses    <- setSelect$tdoses
    tinterval <- setSelect$tinterval
    target    <- setSelect$target
    route     <- setSelect$route
    dose      <- setSelect$dose
    N         <- setSelect$N
    
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
    
    pars [names(pars)] <- pars
    
    #N = nsim # 1000, Number of iterations in the Monte Carlo simulation
    set.seed(123) 
    idata<-MCsim(N=N, pars)
    
    
    PDat<-MC.pred (N=N, pars = pars, drug = drug, idata = idata,
                   tinterval = tinterval, dose = dose, dtimes = tdoses, route = route)
    
    
    return(PDat)
  })
  
  
  ## Parent compound plot-------------------------------------------------------
  plot_output <-  reactive({
    
    # input$action | input$action1
    # isolate({drug = input$drug
    # })
    
    setSelect <- set_select()
    data      <- data_1()
    drug      <- setSelect$drug
    animal    <- setSelect$animal
    tinterval <- setSelect$tinterval
    tdoses    <- setSelect$tdoses
    target    <- setSelect$target
    TOL       <- setSelect$TOL
    
    target1 <- switch(drug,
                      "Flunixin"         = target,
                      "Florfenicol"      = paste0(target,"1"), # Market residual is Florfenicol amine
                      "Penicillin G"     = target)
    
    p <- switch(target,
                Liver  = wdplot(data,target=target1),
                Plasma = wdplot(data,target=target1),
                Kidney = wdplot(data,target=target1),
                Muscle = wdplot(data,target=target1),
                Fat    = wdplot(data,target=target1))
    
    
    
    
    
    end_time    <- data %>% filter_(paste0('round(',target1,'.99',',3)','<=',TOL, 
                                           ' & Time > ((tdoses-1)*tinterval)'))%>%mutate(WDI = Time/24)%>%select(WDI)%>%min()
    
    endtime     <- end_time + tdoses*tinterval/24 
    
    x.intercept <- end_time 
    
    p<-p+
      scale_x_continuous (breaks       = pretty(0:endtime, n = 5),
                          minor_breaks = pretty(0:endtime, n=endtime),
                          label        = pretty((0:endtime), n = 5)-(((tdoses-1)*tinterval/24)),
                          limits       = c(NA, endtime), guide = "axis_minor") +
      
      geom_line(aes(y = TOL),color = 'black', size = 0.5, linetype = 'twodash', show.legend = F)+ #+ # add tolerance line
      scale_y_log10(
        labels = scales::trans_format("log10", scales::math_format(10^.x)),
        limits = c(10^-6, NA)) +
      annotation_logticks(short  = unit(2,"mm"),
                          mid    = unit(3,"mm"),
                          long   = unit(4,"mm"),size = 1,
                          sides  = "l") +
      geom_vline(aes(xintercept  = x.intercept), 
                 size = 1, color = "red", linetype = 2,
                 show.legend = F) +
      labs(x="Time (Day)",y=paste('Concentration in', target, "(mg/kg or ppm)"),face="bold") +
      
      ptheme 
    
    return(p)
    
  }) # END OF PLOT_OUTPUT FUNCTION
  
  # PLOT FUNCTION FOR PLOT_OUTPUT()
  output$plot_1 <- renderPlot({
    withProgress(message = 'Creating plot', 
                 detail='Please Wait...', 
                 value = 10, 
                 expr={input$action | input$action1
                   a=plot_output()
                   
                   return(a)}) 
  })
  
  
  ## PLOT FOR THE TAB OF MODEL PRAMETERS-----------------------------------------  
  plot_output_2 <- reactive({
    
    # input$action | input$action1
    # isolate({
    #   drug = input$drug
    #   animal = input$animal
    # })
    
    setSelect <- set_select()
    drug      <- setSelect$drug
    animal    <- setSelect$animal
    data      <- data_1()
    tinterval <- setSelect$tinterval
    tdoses    <- setSelect$tdoses
    target    <- setSelect$target
    
    #})
    ## PLOT FUNCTION ----------------------------------------------------------
    
    plot.fun <-function (p, target, TOL) {
      
      end_time    <- data %>% filter_(paste0('round(',target,'.99',',2)','<=',TOL, 
                              ' & Time > ((tdoses-1)*tinterval)'))%>%select(Time)%>%min()
      
      endtime     <- end_time/24 + tdoses*tinterval/24 
      
      x.intercept <- end_time/24 

      p1<-p+
        scale_x_continuous (breaks       = pretty(0:endtime, n = 5),
                            minor_breaks = pretty(0:endtime, n=endtime),
                            label        = pretty((0:endtime), n = 5)-((tdoses*tinterval/24)-tinterval/24),
                            limits       = c(NA, endtime), guide = "axis_minor") +
        
        geom_line(aes(y = TOL),color = 'black', size = 0.5, linetype = 'twodash', show.legend = F)+ #+ # add tolerance line
        scale_y_log10(
          labels = scales::trans_format("log10", scales::math_format(10^.x)),
          limits = c(10^-6, NA)) +
        annotation_logticks(short  = unit(2,"mm"),
                            mid    = unit(3,"mm"),
                            long   = unit(4,"mm"),size = 1,
                            sides  = "l") +
        geom_vline(aes(xintercept  = x.intercept), 
                   size = 1, color = "red", linetype = 2,
                   show.legend = F) +
        labs(x="Time (Day)",y=paste('Concentration in', target, "(mg/kg or ppm)"),face="bold") +
        
        ptheme 
      
      

      return(p1)
    }
    
    # END OF PLOT FUNCTION-----------------------------------------------------
    
    
    if (animal == 'Cattle' & drug == 'Flunixin') {
      TOL = c (0.125, 0.125, 0.025)
    } else if (animal == 'Cattle' & drug == 'Florfenicol') {
      TOL = c (3.7, 3.7, 0.3)
    } else if (animal == 'Cattle' & drug == 'Penicillin G') {
      TOL = c (0.05, 0.05, 0.05)
    } else if (animal == 'Swine' & drug == 'Flunixin') {
      TOL = c (0.03, 0.03, 0.025)
    } else if (animal == 'Swine' & drug == 'Florfenicol') {
      TOL = c (2.5, 2.5, 0.2)
    } else if (animal == 'Swine' & drug == 'Penicillin G') {
      TOL = c (0.025, 0.025, 0.025)
    }
    
    
    
    if (drug =='Flunixin' | drug == 'Penicillin G') {
      
      Plot.liver   <- plot.fun(p=wdplot(data, target="Liver"), target="Liver",  TOL=TOL[1])
      Plot.kidney  <- plot.fun(p=wdplot(data, target="Kidney"),target="Kidney", TOL=TOL[2])
      Plot.Muscle  <- plot.fun(p=wdplot(data, target="Muscle"),target="Muscle", TOL=TOL[3])
      Plot.Plasma  <- plot.fun(p=wdplot(data, target="Plasma"),target="Plasma", TOL=TOL[1])
      
    } else if (drug == "Florfenicol"){
      Plot.liver  <- plot.fun(wdplot(data,target="Liver1"), target="Liver1",  TOL=TOL[1])
      Plot.kidney <- plot.fun(wdplot(data,target="Kidney1"),target="Kidney1", TOL=TOL[2])
      Plot.Muscle <- plot.fun(wdplot(data,target="Muscle1"),target="Muscle1", TOL=TOL[3])
      Plot.Plasma <- plot.fun(wdplot(data,target="Plasma1"),target="Plasma1", TOL=TOL[1])
    }
    
    
    return(list(Plot.liver = Plot.liver , 
                Plot.kidney = Plot.kidney,
                Plot.Muscle = Plot.Muscle,
                Plot.Plasma = Plot.Plasma))
  })
  
  
  # PLOT FUNCTION FOR PLOT_OUTPUT()---------------------------------------------
  

  
  output$plot_liver <- renderPlot({
    withProgress(message = 'Creating plot', detail='Please Wait...', 
                 value = 10, 
                 expr={input$action | input$action1
                   a=(plot_output_2()$Plot.liver)
                   return(a) })
    
    
    
  })
  
  
  
  output$plot_kidney <- renderPlot({
    withProgress(message = 'Creating plot', detail='Please Wait...', value = 10, 
                 expr={input$action | input$action1
                   a=(plot_output_2()$Plot.kidney)
                   return(a) })
    
    
    
  })
  
  
  output$plot_muscle <- renderPlot({
    withProgress(message = 'Creating plot', detail='Please Wait...', value = 10, 
                 expr={input$action | input$action1
                   a=(plot_output_2()$Plot.Muscle)
                   return(a) })
    
    
    
  })
  
  
  
  output$plot_plasma <- renderPlot({
    withProgress(message = 'Creating plot', detail='Please Wait...', value = 10, 
                 expr={input$action | input$action1
                   a=(plot_output_2()$Plot.Plasma)
                   return(a) })
    
    
    
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
  
  
  
  
  output$info <- renderText({
    
    # input$action | input$action1
    # isolate({
    #   drug   = input$drug 
    #   animal = input$animal
    # })
    
    setSelect <- set_select()
    drug      <- setSelect$drug
    animal    <- setSelect$animal
    pars      <- setSelect$pars
    tdoses    <- setSelect$tdoses
    tinterval <- setSelect$tinterval
    target    <- setSelect$target
    route     <- setSelect$route
    dose      <- setSelect$dose
    N         <- setSelect$N
    TOL       <- setSelect$TOL
    
    
    drug1   <- switch(drug,
                      "Flunixin"         = "Flunixin",
                      "Florfenicol"      = "Florfenicol amine",
                      "Penicillin G"     = "Penicillin G")
    
    
    b1 <- paste("Figure: Results of PBPK model with Monte Carlo simulation for", drug, "concentrations in", target, "for", animal)
    b2 <- paste("following dosing regiment:", route, "administration at", dose,"mg/kg for", tdoses, "doses with", tinterval, "intervals.")
    b3 <- paste("The median, 1th and 99th percentiles of simulated results were plotted.")
    b4 <- paste("The tolerance is shown on each of panels using the horizontal dashed line", "(Tolearance for liver is", TOL,"µg/g)")
    b5 <- paste("The intersection between the x-axis and the black vertical line indicates the withdraw intervals (WDIs).")
    
    
    
    info <-paste(b1, b2,"", "Note:", b3, b4, b5, sep = "\n")   
    
    return(info)
  })
  
  output$info2 <- renderText({
    
    # input$action | input$action1
    # isolate({
    #   drug = input$drug 
    #   animal = input$animal
    # })
    # 
    setSelect <- set_select()
    
    data      <- data_1()
    animal    <- setSelect$animal
    drug      <- setSelect$drug
    pars      <- setSelect$pars
    tdoses    <- setSelect$tdoses
    tinterval <- setSelect$tinterval
    target    <- setSelect$target
    route     <- setSelect$route
    dose      <- setSelect$dose
    N         <- setSelect$N
    TOL       <- setSelect$TOL
    
    target1 <- switch(drug,
                      "Flunixin"         = target,
                      "Florfenicol"      = paste0(target,"1"),
                      "Penicillin G"     = target)
    
    
    a<-data %>% filter_(paste0('round(',target1,'.99',',3)','<=',TOL, ' & Time > ((tdoses-1)*tinterval)'))%>%
               mutate(WDI = Time/24-((tdoses-1)*tinterval/24)) %>%select(WDI)%>%min()
    
    a1 = round(a, 2)
    
    info <- paste("The recommended withdrawal interval (WDI) is",a1)
    
    
    return(info)
    
  }) 
  
  ###/////////////////////////////////////// Dowload figures
  ## Download the plot on the main panel (wdplot)
  
  
  
  output$downloadwdt <- downloadHandler(
    filename = 'concentration.jpg',
    content = function(file) {
      ggsave(file, plot_output(), width = 40, height = 20, units = "cm")
    },
    contentType = 'image/jpg'
  )
  
  output$downl <- downloadHandler(
    filename = paste(input$drug, 'concentraitons in liver.jpg'),
    content = function(file) {
      ggsave(file, plot_output_2()[[1]])
    },
    contentType = 'image/jpg'
  )
  
  output$downl1 <- downloadHandler(
    filename = "Metabolites concentrations in liver.jpg",
    content = function(file) {
      ggsave(file, plot_output_2()[[2]])
    },
    contentType = 'image/jpg'
  )
  
  output$downk <- downloadHandler(
    filename =  paste(input$drug, 'concentraitons in kidney.jpg'),
    content = function(file) {
      ggsave(file, plot_output_2()[[3]])
    },
    contentType = 'image/jpg'
  )
  
  output$downk1 <- downloadHandler(
    filename = "Metabolites Concentrations in kidney.jpg",
    content = function(file) {
      ggsave(file, plot_output_2()[[4]])
    },
    contentType = 'image/jpg'
  )
  
  
  output$downm <- downloadHandler(
    filename =  "Parent compound concentrations in muscle.jpg",
    content = function(file) {
      ggsave(file, plot_output_2()[[5]])
    },
    contentType = 'image/jpg'
  )
  
  output$downm1 <- downloadHandler(
    filename = "Metabolites concentrations in muscle.jpg",
    content = function(file) {
      ggsave(file, plot_output_2()[[6]])
    },
    contentType = 'image/jpg'
  )
  
  ## Download table------------------------------------------------------------
  
  table <- reactive({
    setSelect <- set_select()

    data      <- data.frame(data_1())
    animal    <- setSelect$animal
    drug      <- setSelect$drug

    if (drug =='Flunixin' | drug == 'Penicillin G') {

      a<-data %>% select(Liver.01 = Liver.01, Liver.50 = Liver.50, Liver.99 = Liver.99,
                         Kidney.01 = Kidney.01, Kidney.50 = Kidney.50, Kidney.99 = Kidney.99,
                         Muscle.01 = Muscle.01,  Muscle.50 = Muscle.50, Muscle.99 = Muscle.99,
                         Plasma.01 = Plasma.01,  Plasma.50 = Plasma.50, Plasma.99 = Plasma.99)
     
    } else {
      a<-data %>% select(Liver1.01 = Liver1.01, Liver1.50 = Liver1.50, Liver1.99 = Liver1.99,
                         Kidney1.01 = Kidney1.01, Kidney1.50 = Kidney1.50, Kidney1.99 = Kidney1.99,
                         Muscle1.01 = Muscle1.01,  Muscle1.50 = Muscle1.50, Muscle1.99 = Muscle1.99,
                         Plasma1.01 = Plasma1.01,  Plasma1.50 = Plasma1.50, Plasma1.99 = Plasma1.99)
      
    }

      return(a)
  })
    
  
  output$table1 <- renderTable ({
    datatable <- table()
  })
  
  
  output$downloadtable <- downloadHandler(
    filename = 'table.csv',
    content = function(file) {
      write.csv(format(data.frame(table()), digits = 4),file)
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
    
    isolate({
      drug   = input$drug
      animal = input$animal
    })
    
    setSelect <- set_select()
    
    pars      <- setSelect$pars
    tdoses    <- setSelect$tdoses
    tinterval <- setSelect$tinterval
    target    <- setSelect$target
    route     <- setSelect$route
    dose      <- setSelect$dose
    N         <- setSelect$N
    data      <- data_1()
    
    
    if (animal == 'Cattle' & drug == 'Flunixin') {
      TOL = c (0.125, 0.125, 0.025)
    } else if (animal == 'Cattle' & drug == 'Florfenicol') {
      TOL = c (3.7, 3.7, 0.3)
    } else if (animal == 'Cattle' & drug == 'Penicillin G') {
      TOL = c (0.05, 0.05, 0.05)
    } else if (animal == 'Swine' & drug == 'Flunixin') {
      TOL = c (0.03, 0.03, 0.025)
    } else if (animal == 'Swine' & drug == 'Florfenicol') {
      TOL = c (2.5, 2.5, 0.2)
    } else if (animal == 'Swine' & drug == 'Penicillin G') {
      TOL = c (0.025, 0.025, 0.025)
    }
    
    
    if (drug == "Florfenicol") {
      
      a_target  = 'Liver1'
      b_target  = 'Muscle1'
      c_target  = 'Kidney1'
      d_target  = 'Plasma1'
      
    } else {
      a_target  = 'Liver'
      b_target  = 'Muscle'
      c_target  = 'Kidney'
      d_target  = 'Plasma'
    }
    
    
    a = data%>%filter_(paste0('round(',a_target,'.99',',3)','<=',TOL[1], ' & Time >(tdoses-1)*tinterval'))%>%select(Time)%>%min()
    b = data%>%filter_(paste0('round(',b_target,'.99',',3)','<=',TOL[3], ' & Time >(tdoses-1)*tinterval'))%>%select(Time)%>%min()
    c = data%>%filter_(paste0('round(',c_target,'.99',',3)','<=',TOL[2], ' & Time >(tdoses-1)*tinterval'))%>%select(Time)%>%min()
    d = data%>%filter_(paste0('round(',d_target,'.99',',3)','<=',TOL[1], ' & Time >(tdoses-1)*tinterval'))%>%select(Time)%>%min()
    
    WDI_data_liver   = round((a/24-((tdoses-1)*tinterval/24)), 2)
    WDI_data_muscle  = round((b/24-((tdoses-1)*tinterval/24)), 2) 
    WDI_data_kidney  = round((c/24-((tdoses-1)*tinterval/24)), 2) 
    WDI_data_plasma  = round((d/24-((tdoses-1)*tinterval/24)), 2) 
    
    
    
    
    return(list(animal, drug, target, N,dose,route,tdoses,tinterval,
                TOL[1], TOL[3], TOL[2], TOL[1], 
                WDI_data_liver, WDI_data_muscle, WDI_data_kidney, WDI_data_plasma))
  })
  
  
  ## End of report output //////////////////////////////////////////////////////// 
  
  
}
