
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

## Loading mrgsolve code
source("Function.R")


## Sidebar: make your main item
sidebar <- dashboardSidebar(
  sidebarMenu(menuItem("Home", tabName = 'plot', icon = icon("home"), selected = TRUE),
              menuItem('Model Parameters', tabName = 'modelparameter', icon = icon('line-chart')),
              menuItem("Output Table", tabName = 'table', icon = icon("table")),
              menuItem("Model Structure", tabName = 'modelstructure', icon = icon("area-chart")),
              menuItem("Code", tabName = 'code', icon = icon("code")),
              # menuItem("Tutorial", tabName = "tutorial", icon = icon("mortar-board")),
              menuItem("Output Report", tabName = "report", icon = icon("file-pdf-o")),
              menuItem("About", tabName = 'about', icon = icon("question-circle"))
  ))


## Body: the content of menuItem
body <- dashboardBody(
  tabItems( # to match upt the menuItem with tabItems 
    
    ## First item content
    tabItem( 
      tabName = 'plot', # the tabName have to match with the tabName of "menuItem"
      sidebarLayout( # sidebarLayout: layout a sidebar (sidebarPanel) and main area (mainPanel)
        sidebarPanel(width = 2,
                     fluidRow(
                       useShinyjs(),
                       tags$h4(strong("Parameters for Therapeutic Scenario"),style = "font-size:20px;"),
                       div(
                         id = "Therapeutic",
                         selectInput(inputId = 'animal', label = 'Species', 
                                     choices = c('Cattle','Swine'), width = 200, selected = "Cattle"
                         ),
                         selectInput(inputId = 'drug', label = 'Drugs', 
                                     choices = c('Flunixin','Florfenicol','Penicillin G',
                                                 '5-OH flunixin', 'Florfenicol amine'), width = 200, selected = "Flunixin"  
                         ),
                         
                         selectInput(inputId = 'target', label = 'Target tissue', 
                                     choices = c('Plasma','Liver','Kidney','Muscle'), width = 200, selected = 'Liver'
                         ),
                         
                         numericInput(inputId = 'simu_time', label = 'Simulation time after last administration', value = 150, 
                                      min = 1, max = 200, step = 1, width = 200
                         ),
                         numericInput(inputId = 'N', label = 'Number of animials', value = 500, min = 1, step = 1, width = 200
                         ),
                         
                         conditionalPanel(condition = "input.drug == 'Flunixin' && input.animal == 'Cattle'",
                                          splitLayout(
                                            selectInput(inputId = 'route', label = 'Administration routes', 
                                                        choices = c('IM', 'IV', 'SC','PO'), width = 200, selected = 'im'
                                            )),
                                          numericInput(inputId = 'doselevel', label = 'Dose level (mg/kg)', 
                                                       value = 2.2, min = 0,  step = 0.01, width = 200
                                          ),
                                          numericInput(inputId = 'tinterval', label = 'Dose interval (h)', value = 24, 
                                                       min = 1, max = 24*10, step = 1, width = 200
                                          ),
                                          numericInput(inputId = 'numdose', label = 'Number of administrations', value = 3, 
                                                       min = 1, max = 5*10, step = 1, width = 200
                                          ),
                                          numericInput(inputId = 'tolerance_Liver', label = 'TL in Liver (mg/kg or ppm)', 
                                                       value = 0.125, min = 0.01, step = 0.01, width = 200
                                          ),
                                          numericInput(inputId = 'tolerance_Muscle', label = 'TL in Muscle (mg/kg or ppm)', 
                                                       value = 0.025, min = 0.01, step = 0.01, width = 200
                                          ), 
                                          numericInput(inputId = 'tolerance_Kidney', label = 'TL in Kidney (mg/kg or ppm)', 
                                                       value = 0.125,  min = 0.01, step = 0.01, width = 200
                                          )
                         ),
                         
                         conditionalPanel(condition = "input.drug == 'Flunixin' && input.animal == 'Swine'",
                                          splitLayout(
                                            selectInput(inputId = 'route_2', label = 'Administration routes', 
                                                        choices = c('IM', 'IV', 'SC','PO'), width = 200, selected = 'im'
                                            )),
                                          numericInput(inputId = 'doselevel_2', label = 'Dose level (mg/kg)', 
                                                       value = 2.2, min = 0,  step = 0.01, width = 200
                                          ),
                                          numericInput(inputId = 'tinterval_2', label = 'Dose interval (h)', 
                                                       value = 24, min = 1, max = 24*10, step = 1, width = 200
                                          ),
                                          numericInput(inputId = 'numdose_2', label = 'Number of administrations', 
                                                       value = 3, min = 1, max = 5*10, step = 1, width = 200
                                          ),
                                          numericInput(inputId = 'tolerance_Liver_2', label = 'TL in Liver (mg/kg or ppm)', 
                                                       value = 0.03, min = 0.01, step = 0.01, width = 200
                                          ),
                                          numericInput(inputId = 'tolerance_Muscle_2', label = 'TL in Muscle (mg/kg or ppm)', 
                                                       value = 0.025, min = 0.01, step = 0.01, width = 200
                                          ), 
                                          numericInput(inputId = 'tolerance_Kidney_2', label = 'TL in Kidney (mg/kg or ppm)', 
                                                       value = 0.03, min = 0.01, step = 0.01, width = 200
                                          )
                         ),
                         
                         conditionalPanel(condition = "input.drug == '5-OH flunixin' && input.animal == 'Cattle'",
                                          splitLayout(
                                            selectInput(inputId = 'route', label = 'Administration routes', 
                                                        choices = c('IM', 'IV', 'SC','PO'), width = 200, selected = 'im'
                                            )),
                                          numericInput(inputId = 'doselevel', label = 'Dose level (mg/kg)', 
                                                       value = 2.2, min = 0,  step = 0.01, width = 200
                                          ),
                                          numericInput(inputId = 'tinterval', label = 'Dose interval (h)', value = 24, 
                                                       min = 1, max = 24*10, step = 1, width = 200
                                          ),
                                          numericInput(inputId = 'numdose', label = 'Number of administrations', value = 3, 
                                                       min = 1, max = 5*10, step = 1, width = 200
                                          ),
                                          numericInput(inputId = 'tolerance_Liver', label = 'TL in Liver (mg/kg or ppm)', value = 0.125, 
                                                       min = 0.01, step = 0.01, width = 200
                                          ),
                                          numericInput(inputId = 'tolerance_Muscle', label = 'TL in Muscle (mg/kg or ppm)', value = 0.025, 
                                                       min = 0.01, step = 0.01, width = 200
                                          ), 
                                          numericInput(inputId = 'tolerance_Kidney', label = 'TL in Kidney (mg/kg or ppm)', value = 0.125, 
                                                       min = 0.01, step = 0.01, width = 200
                                          )
                         ),
                         
                         conditionalPanel(condition = "input.drug == '5-OH flunixin' && input.animal == 'Swine'",
                                          splitLayout(
                                            selectInput(inputId = 'route_2', label = 'Administration routes', 
                                                        choices = c('IM', 'IV', 'SC','PO'), width = 200, selected = 'im'
                                            )),
                                          numericInput(inputId = 'doselevel_2', label = 'Dose level (mg/kg)', 
                                                       value = 2.2, min = 0,  step = 0.01, width = 200
                                          ),
                                          numericInput(inputId = 'tinterval_2', label = 'Dose interval (h)', value = 24, 
                                                       min = 1, max = 24*10, step = 1, width = 200
                                          ),
                                          numericInput(inputId = 'numdose_2', label = 'Number of administrations', value = 3, 
                                                       min = 1, max = 5*10, step = 1, width = 200
                                          ),
                                          numericInput(inputId = 'tolerance_Liver_2', label = 'TL in Liver (mg/kg or ppm)', value = 0.03, 
                                                       min = 0.01, step = 0.01, width = 200
                                          ),
                                          numericInput(inputId = 'tolerance_Muscle_2', label = 'TL in Muscle (mg/kg or ppm)', value = 0.025, 
                                                       min = 0.01, step = 0.01, width = 200
                                          ), 
                                          numericInput(inputId = 'tolerance_Kidney_2', label = 'TL in Kidney (mg/kg or ppm)', value = 0.03, 
                                                       min = 0.01, step = 0.01, width = 200
                                          )
                         ),
                         
                         conditionalPanel(condition = "input.drug == 'Florfenicol' && input.animal == 'Cattle'",
                                          splitLayout(
                                            selectInput(inputId = 'route_3', label = 'Administration routes', 
                                                        choices = c('IM', 'IV', 'SC','PO'), width = 200, selected = 'im'
                                            )),
                                          numericInput(inputId = 'doselevel_3', label = 'Dose level (mg/kg)', 
                                                       value = 20, min = 0,  step = 0.01, width = 200
                                          ),
                                          numericInput(inputId = 'tinterval_3', label = 'Dose interval (h)', 
                                                       value = 48, min = 1, max = 24*10, step = 1, width = 200
                                          ),
                                          numericInput(inputId = 'numdose_3', label = 'Number of administrations', 
                                                       value = 2, min = 1, max = 5*10, step = 1, width = 200
                                          ),
                                          numericInput(inputId = 'tolerance_Liver_3', label = 'TL in Liver (mg/kg or ppm)', 
                                                       value = 3.7, min = 0.01, step = 0.01, width = 200
                                          ),
                                          numericInput(inputId = 'tolerance_Muscle_3', label = 'TL in Muscle (mg/kg or ppm)', 
                                                       value = 0.3, min = 0.01, step = 0.01, width = 200
                                          ), 
                                          numericInput(inputId = 'tolerance_Kidney_3', label = 'TL in Kidney (mg/kg or ppm)', 
                                                       value = 0.3, min = 0.01, step = 0.01, width = 200
                                          )
                         ),
                         
                         conditionalPanel(condition = "input.drug == 'Florfenicol' && input.animal == 'Swine'",
                                          splitLayout(
                                            selectInput(inputId = 'route_4', label = 'Administration routes', 
                                                        choices = c('IM', 'IV', 'SC','PO'), width = 200, selected = 'PO'
                                            )),
                                          numericInput(inputId = 'doselevel_4', label = 'Dose level (mg/kg)', 
                                                       value = 10, min = 0,  step = 0.01, width = 200
                                          ),
                                          numericInput(inputId = 'tinterval_4', label = 'Dose interval (h)', 
                                                       value = 24, min = 1, max = 24*10, step = 1, width = 200
                                          ),
                                          numericInput(inputId = 'numdose_4', label = 'Number of administrations', 
                                                       value = 5, min = 1, max = 5*10, step = 1, width = 200
                                          ),
                                          numericInput(inputId = 'tolerance_Liver_4', label = 'TL in Liver (mg/kg or ppm)', 
                                                       value = 2.5, min = 0.01, step = 0.01, width = 200
                                          ),
                                          numericInput(inputId = 'tolerance_Muscle_4', label = 'TL in Muscle (mg/kg or ppm)', 
                                                       value = 0.2, min = 0.01, step = 0.01, width = 200
                                          ), 
                                          numericInput(inputId = 'tolerance_Kidney_4', label = 'TL in Kidney (mg/kg or ppm)', 
                                                       value = 0.2, min = 0.01, step = 0.01, width = 200
                                          )
                         ),
                         
                         conditionalPanel(condition = "input.drug == 'Florfenicol amine' && input.animal == 'Cattle'",
                                          splitLayout(
                                            selectInput(inputId = 'route_3', label = 'Administration routes', 
                                                        choices = c('IM', 'IV', 'SC','PO'), width = 200, selected = 'im'
                                            )),
                                          numericInput(inputId = 'doselevel_3', label = 'Dose level (mg/kg)', 
                                                       value = 20, min = 0,  step = 0.01, width = 200
                                          ),
                                          numericInput(inputId = 'tinterval_3', label = 'Dose interval (h)', 
                                                       value = 48, min = 1, max = 24*10, step = 1, width = 200
                                          ),
                                          numericInput(inputId = 'numdose_3', label = 'Number of administrations', 
                                                       value = 2, min = 1, max = 5*10, step = 1, width = 200
                                          ),
                                          numericInput(inputId = 'tolerance_Liver_3', label = 'TL in Liver (mg/kg or ppm)', 
                                                       value = 3.7, min = 0.01, step = 0.01, width = 200
                                          ),
                                          numericInput(inputId = 'tolerance_Muscle_3', label = 'TL in Muscle (mg/kg or ppm)', 
                                                       value = 0.3, min = 0.01, step = 0.01, width = 200
                                          ), 
                                          numericInput(inputId = 'tolerance_Kidney_3', label = 'TL in Kidney (mg/kg or ppm)', 
                                                       value = 0.3, min = 0.01, step = 0.01, width = 200
                                          )
                         ),
                         
                         conditionalPanel(condition = "input.drug == 'Florfenicol amine' && input.animal == 'Swine'",
                                          splitLayout(
                                            selectInput(inputId = 'route_4', label = 'Administration routes', 
                                                        choices = c('IM', 'IV', 'SC','PO'), width = 200, selected = 'PO'
                                            )),
                                          numericInput(inputId = 'doselevel_4', label = 'Dose level (mg/kg)', 
                                                       value = 10, min = 0,  step = 0.01, width = 200
                                          ),
                                          numericInput(inputId = 'tinterval_4', label = 'Dose interval (h)', 
                                                       value = 24, min = 1, max = 24*10, step = 1, width = 200
                                          ),
                                          numericInput(inputId = 'numdose_4', label = 'Number of administrations', 
                                                       value = 5, min = 1, max = 5*10, step = 1, width = 200
                                          ),
                                          numericInput(inputId = 'tolerance_Liver_4', label = 'TL in Liver (mg/kg or ppm)', 
                                                       value = 3.7, min = 0.01, step = 0.01, width = 200
                                          ),
                                          numericInput(inputId = 'tolerance_Muscle_4', label = 'TL in Muscle (mg/kg or ppm)', 
                                                       value = 0.3, min = 0.01, step = 0.01, width = 200
                                          ), 
                                          numericInput(inputId = 'tolerance_Kidney_4', label = 'TL in Kidney (mg/kg or ppm)', 
                                                       value = 0.3, min = 0.01, step = 0.01, width = 200
                                          )
                         ),
                         conditionalPanel(condition = "input.drug == 'Penicillin G' && input.animal == 'Cattle'",
                                          splitLayout(
                                            selectInput(inputId = 'route_5', label = 'Administration routes', 
                                                        choices = c('IM', 'IV', 'SC','PO'), width = 200, selected = 'im'
                                            )),
                                          numericInput(inputId = 'doselevel_5', label = 'Dose level (mg/kg)', 
                                                       value = 6.5, min = 0,  step = 0.01, width = 200
                                          ),
                                          numericInput(inputId = 'tinterval_5', label = 'Dose interval (h)', 
                                                       value = 24, min = 1, max = 24*10, step = 1, width = 200
                                          ),
                                          numericInput(inputId = 'numdose_5', label = 'Number of administrations', 
                                                       value = 5, min = 1, max = 5*10, step = 1, width = 200
                                          ),
                                          numericInput(inputId = 'tolerance_Liver_5', label = 'TL in Liver (mg/kg or ppm)', 
                                                       value = 0.05, min = 0.01, step = 0.01, width = 200
                                          ),
                                          numericInput(inputId = 'tolerance_Muscle_5', label = 'TL in Muscle (mg/kg or ppm)', 
                                                       value = 0.05, min = 0.01, step = 0.01, width = 200
                                          ), 
                                          numericInput(inputId = 'tolerance_Kidney_5', label = 'TL in Kidney (mg/kg or ppm)', 
                                                       value = 0.05, min = 0.01, step = 0.01, width = 200
                                          )
                         ),
                         
                         conditionalPanel(condition = "input.drug == 'Penicillin G' && input.animal == 'Swine'",
                                          splitLayout(
                                            selectInput(inputId = 'route_6', label = 'Administration routes', 
                                                        choices = c('IM', 'IV', 'SC','PO'), width = 200, selected = 'im'
                                            )),
                                          numericInput(inputId = 'doselevel_6', label = 'Dose level (mg/kg)', 
                                                       value = 6.5, min = 0,  step = 0.01, width = 200
                                          ),
                                          numericInput(inputId = 'tinterval_6', label = 'Dose interval (h)', 
                                                       value = 24, min = 1, max = 24*10, step = 1, width = 200
                                          ),
                                          numericInput(inputId = 'numdose_6', label = 'Number of administrations', 
                                                       value = 5, min = 1, max = 5*10, step = 1, width = 200
                                          ),
                                          numericInput(inputId = 'tolerance_Liver_6', label = 'TL in Liver (mg/kg or ppm)', 
                                                       value = 0.05, min = 0.01, step = 0.01, width = 200
                                          ),
                                          numericInput(inputId = 'tolerance_Muscle_6', label = 'TL in Muscle (mg/kg or ppm)', 
                                                       value = 0.05, min = 0.01, step = 0.01, width = 200
                                          ), 
                                          numericInput(inputId = 'tolerance_Kidney_6', label = 'TL in Kidney (mg/kg or ppm)', 
                                                       value = 0.05, min = 0.01, step = 0.01, width = 200
                                          )
                         ),
                         
                         conditionalPanel(condition = "input.drug == 'Penicillin G'",
                                          tags$strong(tags$a(href="https://www.accessdata.fda.gov/scripts/cdrh/cfdocs/cfcfr/CFRSearch.cfm?fr=522.1696a", 'Tolerance by FDA', target="_blank"
                                          ))),
                         conditionalPanel(condition = "input.drug == 'Flunixin'",
                                          tags$strong(tags$a(href="https://www.accessdata.fda.gov/scripts/cdrh/cfdocs/cfcfr/CFRSearch.cfm?fr=556.286", 'Tolerance by FDA', target="_blank"
                                          ))),
                         
                         conditionalPanel(condition = "input.drug == 'Florfenicol'",
                                          tags$strong(tags$a(href="https://www.accessdata.fda.gov/scripts/cdrh/cfdocs/cfcfr/CFRSearch.cfm?CFRPart=556&showFR=1", 'Tolerance by FDA', target="_blank"
                                          )),
                                          tags$strong(tags$p("Use of florfenicol in cattle and swine is extralabel in the USA"))
                                          
                         ),
                         
                         tags$p(" "),
                         
                         actionButton(inputId = 'action',label = 'Apply Changes',icon("paper-plane"), style="width:200px"),
                         
                         actionButton(inputId = 'reset', label = 'Default Values',icon("refresh"), style="width:200px") 
                       ))
        ), # end of sidebarPanel
        
        mainPanel(width = 8,
                  fixedRow(
                    box(title = tags$h1(strong('Extralabel Withdrawal Interval Plot'), style = "font-size:25px;"),
                        plotOutput(outputId = 'wdtplot', click = "plot_click"#, width="80%",
                                   #dblclick = "plot_dblclick",
                                   #hover = "plot_hover",
                                   #brush = "plot_brush"
                        ), height = 900, width = 3000,
                        br(), tags$h1(strong(verbatimTextOutput("info")),style = "font-size:20px;"), 
                        #conditionalPanel(condition = "input.animal == 'Cattle'",
                                         tags$strong(verbatimTextOutput("info2")),
                                         #tags$strong("Labelled dose: 22 mg/kg Oral dose via feed continuously for 7 - 14 days"), br(),br(),
                                         #tags$strong("Labelled WDT: 5 day"), br(), br(),
                                         #tags$strong("Tolerance (mg/kg or ppm): Muscle = 2,  Liver = 6, Kidney = 12"),br(), br(),
                                         #tags$strong(tags$a(href="https://www.accessdata.fda.gov/scripts/cdrh/cfdocs/cfcfr/cfrsearch.cfm?fr=556.500", 'Click here to see florfenicol tolerances in edible tissues by FDA', target="_blank"
                                         #))
                                         
                        #),
                        #conditionalPanel(condition = "input.animal == 'Swine'",
                                         #tags$strong(textOutput("info3")),
                                         #tags$p("Labelled WDT (day): 5"),
                                         #tags$p("Labelled dose (mg/kg): 6.5"),
                                         #tags$p("Tolerance (ppm or ug/g): 2 for Muscle"),
                                         #tags$p("Tolerance (ppm or ug/g): 6 for Liver"),
                                         #tags$p("Tolerance (ppm or ug/g): 12 for Kidney")
                        #),
                        
                        tags$strong("Disclaimer: The generic PBPK model is calibrated for the FLU, FLU an PG in cattle or swine"),br(),
                        tags$strong("Abbreviations: FLU- Flunxin, FLO-Florfenicol, FLOA- Florefenicol Amine", "TOL-Tolerance"),br(), br(), 
                        downloadButton('downloadwdt','Download')
                    )
                  ) # end of fixedRow
        ) # end of mainPanel
      )), # end of the tabItem of "plot"
    
    
    ## Second item content
    tabItem(
      tabName = 'modelparameter',
      sidebarLayout(
        sidebarPanel(width = 2,
                     fluidRow(
                       useShinyjs(),
                       tabsetPanel(id = 'tabs', type = 'tabs',
                                   tabPanel(title = 'Physiological parameters'),
                                   tabPanel(title = 'Chemical-specific parameters')),br(),
                       div(
                         id = "all_parameters",
                         actionButton(inputId = 'action1',label = 'Apply Changes',icon("paper-plane"), style="width:200px"),
                         
                         actionButton(inputId = 'reset1', label = 'Default Values',icon("refresh"), style="width:200px") ,
                         br(),
                         br(),
                         
                         conditionalPanel(condition = "input.tabs == 'Physiological parameters'&& input.animal == 'Cattle'",
                                          splitLayout(
                                            numericInput(inputId = 'BW', label = 'BW', value = 250, 
                                                         step = 0.01, width = 100),
                                            
                                            numericInput(inputId = 'QCC', label = 'QCC', value = 5.97,
                                                         step = 0.01, width = 100)
                                          ),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'QLC', label = 'QLC', value = 0.405,
                                                         step = 0.01, width = 100),
                                            
                                            numericInput(inputId = 'QKC', label = 'QKC', value = 0.09,
                                                         step = 0.01, width = 100)
                                          ),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'QMC', label = 'QMC', value = 0.18,
                                                         step = 0.01, width = 100),
                                            
                                            numericInput(inputId = 'QFC', label = 'QFC', value = 0.08,
                                                         step = 0.01, width = 100)
                                          ),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'QRestC', label = 'QRestC', value = 0.245,
                                                         step = 0.01, width = 100),
                                            
                                            numericInput(inputId = 'VLC', label = 'VLC', value = 0.014,
                                                         step = 0.001, width = 100)
                                            
                                          ),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'VKC', label = 'VKC', value = 0.0025,
                                                         step = 0.01, width = 100),
                                            
                                            numericInput(inputId = 'VFC', label = 'VFC', value = 0.15,
                                                         step = 0.01, width = 100)
                                          ),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'VMC', label = 'VMC', value = 0.27,
                                                         step = 0.01, width = 100),
                                            
                                            numericInput(inputId = 'VbloodC', label = 'VbC', value = 0.040,
                                                         step = 0.001, width = 100)
                                            
                                          ),
                                          
                                          
                                          splitLayout(
                                            numericInput(inputId = 'VRestC', label = 'VRestC', value = 0.5235,
                                                         step = 0.01, width = 100)
                                            
                                          )
                         ), 
                         
                         conditionalPanel(condition = "input.tabs == 'Physiological parameters'&& input.animal == 'Swine'",
                                          splitLayout(
                                            numericInput(inputId = 'BW', label = 'BW', value = 33, 
                                                         step = 0.01, width = 100
                                            )),
                                          splitLayout(
                                            numericInput(inputId = 'QCC', label = 'QCC', value = 8.543,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'QLC', label = 'QLC', value = 0.273,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'QKC', label = 'QKC', value = 0.116,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'QMC', label = 'QMC', value = 0.293,
                                                         step = 0.01, width = 100
                                                         
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'QFC', label = 'QFC', value = 0.128,
                                                         step = 0.01, width = 100
                                                         
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'QRestC', label = 'QRestC', value = 0.190,
                                                         step = 0.01, width = 100
                                                         
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'VLC', label = 'VLC', value = 0.023,
                                                         step = 0.001, width = 100
                                                         
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'VKC', label = 'VKC', value = 0.0045,
                                                         step = 0.01, width = 100
                                                         
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'VFC', label = 'VFC', value = 0.235,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'VMC', label = 'VMC', value = 0.355,
                                                         step = 0.01, width = 100
                                                         
                                            )),
                                          
                                          
                                          splitLayout(
                                            numericInput(inputId = 'VbloodC', label = 'VbC', value = 0.050,
                                                         step = 0.001, width = 100
                                                         
                                            )),
                                          
                                          
                                          splitLayout(
                                            numericInput(inputId = 'VRestC', label = 'VRestC', value = 0.333,
                                                         step = 0.01, width = 100
                                                         
                                            ))), 
                         
                         conditionalPanel(condition = "input.tabs == 'Chemical-specific parameters'&& input.animal == 'Cattle'&& input.drug == 'Flunixin'",
                                          splitLayout(
                                            numericInput(inputId = 'Fracim', label = 'Fracim', value = 0.50,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'Fracsc', label = 'Fracsc', value = 0.65,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'fR', label = 'fR', value = 0.05,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'PL', label = 'PL', value = 3.06,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'PK', label = 'PK', value = 6.5,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'PM', label = 'PM', value = 0.287,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'PF', label = 'PF', value = 0.16,
                                                         step = 0.0001, width = 100
                                            )),
                                          
                                          
                                          splitLayout(
                                            numericInput(inputId = 'PRest', label = 'PRest', value = 6.7,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          
                                          splitLayout(
                                            numericInput(inputId = 'Kim', label = 'Kim', value = 0.26 ,
                                                         step = 0.001, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'Ksc', label = 'Ksc', value = 0.40,
                                                         
                                                         step = 0.001, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'Kdissim', label = 'Kdissim', value = 1e-5,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'Kdisssc', label = 'Kdisssc', value = 1e-5,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'KehcC', label = 'KehcC', value = 8.87e-5,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'KmetC', label = 'KmetC', value = 0.20,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'KbileC', label = 'KbileC', value = 0.010,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'KurineC', label = 'KurineC', value = 0.5,
                                                         step = 0.001, width = 100
                                                         
                                            ))),
                         
                         conditionalPanel(condition = "input.tabs == 'Chemical-specific parameters'&& input.animal == 'Cattle'&& input.drug == 'Florfenicol'",
                                          splitLayout(
                                            numericInput(inputId = 'Fracim', label = 'Fracim', value = 0.57,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'Fracsc', label = 'Fracsc', value = 0.33,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'fR', label = 'fR', value = 0.80,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'PL', label = 'PL', value = 2.2,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'PK', label = 'PK', value = 2.52,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'PM', label = 'PM', value = 1.30,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'PF', label = 'PF', value = 0.60,
                                                         step = 0.0001, width = 100
                                            )),
                                          
                                          
                                          splitLayout(
                                            numericInput(inputId = 'PRest', label = 'PRest', value = 0.11,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          
                                          splitLayout(
                                            numericInput(inputId = 'Kim', label = 'Kim', value = 0.15 ,
                                                         step = 0.001, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'Ksc', label = 'Ksc', value = 0.09,
                                                         
                                                         step = 0.001, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'Kdissim', label = 'Kdissim', value = 0.01,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'Kdisssc', label = 'Kdisssc', value = 0.0054,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'KehcC', label = 'KehcC', value = 0.05,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'KmetC', label = 'KmetC', value = 0.0025,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'KbileC', label = 'KbileC', value = 1e-2,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'KurineC', label = 'KurineC', value = 0.05,
                                                         step = 0.001, width = 100
                                                         
                                            ))),
                         
                         conditionalPanel(condition = "input.tabs == 'Chemical-specific parameters'&& input.animal == 'Cattle'&& input.drug == 'Penicillin G'",
                                          splitLayout(
                                            numericInput(inputId = 'Fracim', label = 'Fracim', value = 0.775,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'Fracsc', label = 'Fracsc', value = 0.209,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'fR', label = 'fR', value = 0.50,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'PL', label = 'PL', value = 3.04,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'PK', label = 'PK', value = 4.96,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'PM', label = 'PM', value = 0.09,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'PF', label = 'PF', value = 0.152,
                                                         step = 0.0001, width = 100
                                            )),
                                          
                                          
                                          splitLayout(
                                            numericInput(inputId = 'PRest', label = 'PRest', value = 0.81,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          
                                          splitLayout(
                                            numericInput(inputId = 'Kim', label = 'Kim', value = 0.05 ,
                                                         step = 0.001, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'Ksc', label = 'Ksc', value = 0.25,
                                                         
                                                         step = 0.001, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'Kdissim', label = 'Kdissim', value = 0.007,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'Kdisssc', label = 'Kdisssc', value = 0.005,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'KehcC', label = 'KehcC', value = 0.09,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'KmetC', label = 'KmetC', value = 0.73,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'KbileC', label = 'KbileC', value = 0.129,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'KurineC', label = 'KurineC', value = 0.83,
                                                         step = 0.001, width = 100
                                                         
                                            ))),
                         
                         
                         conditionalPanel(condition = "input.tabs == 'Chemical-specific parameters'&& input.animal == 'Swine'&& input.drug == 'Flunixin'",
                                          splitLayout(
                                            numericInput(inputId = 'Fracim', label = 'Fracim', value = 0.65,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'Fracsc', label = 'Fracsc', value = 0.50,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'fR', label = 'fR', value = 0.05,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'PL', label = 'PL', value = 1.08,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'PK', label = 'PK', value = 3.03,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'PM', label = 'PM', value = 0.09,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'PF', label = 'PF', value = 0.11,
                                                         step = 0.0001, width = 100
                                            )),
                                          
                                          
                                          splitLayout(
                                            numericInput(inputId = 'PRest', label = 'PRest', value = 12.6,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          
                                          splitLayout(
                                            numericInput(inputId = 'Kim', label = 'Kim', value = 1 ,
                                                         step = 0.001, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'Ksc', label = 'Ksc', value = 0.4,
                                                         
                                                         step = 0.001, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'Kdissim', label = 'Kdissim', value = 1e-5,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'Kdisssc', label = 'Kdisssc', value = 1e-5,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'KehcC', label = 'KehcC', value = 0.15,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'KmetC', label = 'KmetC', value = 0.02,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'KbileC', label = 'KbileC', value = 0.01,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'KurineC', label = 'KurineC', value = 0.05,
                                                         step = 0.001, width = 100
                                                         
                                            ))),
                         
                         conditionalPanel(condition = "input.tabs == 'Chemical-specific parameters'&& input.animal == 'Swine'&& input.drug == 'Florfenicol'",
                                          splitLayout(
                                            numericInput(inputId = 'Fracim', label = 'Fracim', value = 0.61,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'Fracsc', label = 'Fracsc', value = 0.50,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'fR', label = 'fR', value = 0.81,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'PL', label = 'PL', value = 0.12,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'PK', label = 'PK', value = 0.68,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'PM', label = 'PM', value = 0.7,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'PF', label = 'PF', value = 0.29,
                                                         step = 0.0001, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'PRest', label = 'PRest', value = 1.31,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          
                                          splitLayout(
                                            numericInput(inputId = 'Kim', label = 'Kim', value = 0.11 ,
                                                         step = 0.001, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'Ksc', label = 'Ksc', value = 0.128,
                                                         
                                                         step = 0.001, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'Kdissim', label = 'Kdissim', value = 0.02,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'Kdisssc', label = 'Kdisssc', value = 0.65,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'KehcC', label = 'KehcC', value = 0.23,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'KmetC', label = 'KmetC', value = 0.02,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'KbileC', label = 'KbileC', value = 0.005,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'KurineC', label = 'KurineC', value = 0.496,
                                                         step = 0.001, width = 100
                                                         
                                            ))),
                         
                         conditionalPanel(condition = "input.tabs == 'Chemical-specific parameters'&& input.animal == 'Swine'&& input.drug == 'Penicillin G'",
                                          splitLayout(
                                            numericInput(inputId = 'Fracim', label = 'Fracim', value = 0.67,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'Fracsc', label = 'Fracsc', value = 0.50,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'fR', label = 'fR', value = 0.634,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'PL', label = 'PL', value = 0.03,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'PK', label = 'PK', value = 0.84,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'PM', label = 'PM', value = 0.07,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'PF', label = 'PF', value = 0.06,
                                                         step = 0.0001, width = 100
                                            )),
                                          
                                          
                                          splitLayout(
                                            numericInput(inputId = 'PRest', label = 'PRest', value = 0.73,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          
                                          splitLayout(
                                            numericInput(inputId = 'Kim', label = 'Kim', value = 0.07 ,
                                                         step = 0.001, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'Ksc', label = 'Ksc', value = 0.25,
                                                         
                                                         step = 0.001, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'Kdissim', label = 'Kdissim', value = 0.007,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'Kdisssc', label = 'Kdisssc', value = 0.005,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'KehcC', label = 'KehcC', value = 0.463,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'KmetC', label = 'KmetC', value = 0.612,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'KbileC', label = 'KbileC', value = 0.89,
                                                         step = 0.01, width = 100
                                            )),
                                          
                                          splitLayout(
                                            numericInput(inputId = 'KurineC', label = 'KurineC', value = 1.19,
                                                         step = 0.001, width = 100
                                            )))
                         
                       )
                     )),
        
        mainPanel(width = 10,
                  tabItem(tabName = 'modelparameter',
                          fluidRow(
                            box(downloadButton('downl','Download'),
                                title = textOutput("title1"), if(FALSE) {'FLU Concentrations in liver'}, status = "primary", solidHeader = TRUE,
                                plotOutput(outputId = 'wdplot_liver', height = 350)
                            ),
                            
                            
                            box(downloadButton('downk','Download'),
                                title = textOutput("title2"), if(FALSE) {'FLU Concentrations in kidney'},status = "primary", solidHeader = TRUE,
                                plotOutput(outputId = 'wdplot_kidney', height = 350)
                            ),
                            
                            
                            box(downloadButton('downm','Download'),
                                title = textOutput("title3"), if(FALSE){'FLU Concentrations in muscle'},status = "primary", solidHeader = TRUE,
                                plotOutput(outputId = 'wdplot_muscle', height = 350)
                            ),
                            
                            box(downloadButton('downm1','Download'),
                                title = textOutput("title4"), if(FALSE){'FLU Concentrations in plasma'},status = "primary", solidHeader = TRUE,
                                plotOutput(outputId = 'wdplot_plasma', height = 350)
                            )
                          )
                  )
        )
        
      )
      
    ),
    
    ##
    tabItem(
      tabName = 'table',
      box(width = "100%", height = "100%", status = "primary", solidHeader = TRUE, title = "Concentration in tissues",
          downloadButton('downloadtable','Download'),
          tableOutput('table1')
      )
    ), # end ot the tabItem of "table"
    
    
    
    ## Third item content
    tabItem(
      tabName = 'modelstructure',
      fluidRow(
        box(width = "100%", status = "primary",title= "Model Structure", align="center", 
            solidHeader = TRUE, tags$img(src = 'ModStructure1.png', width = '70%',height = '120%'),
            tags$h2(strong("
                   Fig. 1. A schematic diagram of the physiologically based pharmacokinetic (PBPK) model for Flunixin, 
                   Florfenicol and Penicillin G in cattle and swine. Four different administration routes
                   including intraPlasma (IV), subcuteneous (SC), intramuscular (IM) injections and oral (PO) administration are presented in the
                   model. 
                   "), align = 'left', style = "font-family: 'Times', serif; font-weight: 500; 
                   font-size: 20px; line-height: 1; color: #404040;")
        )
      )
    ),
    
    tabItem(
      tabName = 'about',
      fluidRow(
        box(width = 10,status = "primary",title= "About", 
            solidHeader = TRUE, tags$img(src = 'farad.png', height = "50px"),
            tags$br(),
            tags$strong("
                   The Extralabel Withdrawal Interval Simulator is based on the generic PBPK model for Flunixin, Florfenicol and Penicillin G
                   in cattle and swine.
                   "), 
            
            tags$br(),
            tags$strong("
                   The project is supported by the US Department of Agriculture National Institute of Food and Agriculture for the Food Animal Residue
                   Avoidance Databank (FARAD) Program (http://www.farad.org/).
                   "),tags$br(), tags$strong("Contact Info:",tags$br(),"Zhoumeng Lin, BMed, PhD, DABT",tags$br(),"linzhoumeng@ufl.edu",tags$br(),
                                             "(352) 273-9188")
        )
      )
    ),
    
    
    tabItem(
      tabName = 'report',
      box(width = 8, status = "primary", solidHeader = TRUE, title="Output Report", 
          radioButtons('format', 'Document format', c('HTML', 'Word', 'PDF'),
                       inline = TRUE),
          downloadButton('downloadreport', "Download Report")
      )
    ),
    
    tabItem(
      tabName = 'code',
      box( width = NULL, status = "primary", solidHeader = TRUE, title="PBPK model",
           downloadButton('downloadcode', 'Download'),
           br(),br(),
           pre(includeText("GenPBPK.R"))
      ))
    
    
    
  ))



# The dashboardPage() function expects three components: 
# header, sidebar, and body            
ui<-dashboardPage( 
  skin = 'yellow',
  dashboardHeader(title = "Withdrawal Interval Simulator"),
  sidebar,
  body
)






