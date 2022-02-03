## Loading global class, which containing required function
source("Global.R")
source("ui_design.R")

# The dashboardPage() function expects three components: 
# header, sidebar, and body            
ui<-dashboardPage( 
  skin = 'yellow',
  dashboardHeader(title = "Withdrawal Interval Simulator"),
  
  ## Sidebar: make your main item
  sidebar <- dashboardSidebar(
    width = 200,
    
    tags$head( 
      tags$style(HTML(".main-sidebar { font-size: 18px;}")) #change the font size to 15
    ),
    
    sidebarMenu(menuItem("Plot", tabName = 'plot', icon = icon("line-chart"), selected = TRUE),
                menuItem('Model Parameters', tabName = 'modelparameter', icon = icon('adjust')),
                menuItem("Output Table", tabName = 'table', icon = icon("table")),
                menuItem("Model Structure", tabName = 'modelstructure', icon = icon("area-chart")),
                menuItem("Code", tabName = 'code', icon = icon("code")),
                menuItem("Output Report", tabName = "report", icon = icon("file-pdf-o")),
                menuItem("About", tabName = 'about', icon = icon("question-circle"))
    )),
  
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
                         tags$h4(strong("Parameters for Therapeutic Scenario"),style = "font-size:15px;"),
                         div(
                           id = "Therapeutic",
                           selectInput(inputId = 'animal', label = 'Species', 
                                       choices = c('Cattle','Swine'), width = 200, selected = "Cattle"),
                           selectInput(inputId = 'drug', label = 'Drugs', 
                                       choices = c('Flunixin','Florfenicol','Penicillin G'), width = 200, selected = "Flunixin"),
                           
                           numericInput(inputId = 'simu_time', label = 'Simulation time after last administration', value = 150, 
                                        min = 1, max = 200, step = 1, width = 200),
                           
                           numericInput(inputId = 'N', label = 'Number of animials', value = 1000, 
                                        min = 1, step = 1, width = 200),
                           
                           conditionalPanel(
                             condition = "input.drug == 'Flunixin' && input.animal == 'Cattle'",
                             SidebarDesign (a = 'BC', b = 'FLU')
                           ),
                           
                           conditionalPanel(
                             condition = "input.drug == 'Flunixin' && input.animal == 'Swine'",
                             SidebarDesign (a = 'SW', b = 'FLU')
                           ),
                           
                           conditionalPanel(
                             condition = "input.drug == 'Florfenicol' && input.animal == 'Cattle'",
                             SidebarDesign (a = 'BC', b = 'FLO')
                           ),
                           
                           conditionalPanel(
                             condition = "input.drug == 'Florfenicol' && input.animal == 'Swine'",
                             SidebarDesign (a = 'SW', b = 'FLO')
                           ),
                           
                           conditionalPanel(
                             condition = "input.drug == 'Penicillin G' && input.animal == 'Cattle'",
                             SidebarDesign (a = 'BC', b = 'PG')
                           ),
                           
                           conditionalPanel(
                             condition = "input.drug == 'Penicillin G' && input.animal == 'Swine'",
                             SidebarDesign (a = 'SW', b = 'PG')
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
          
          mainPanel(width = 9,
                    fixedRow(
                      box(title = tags$h1(strong('Extralabel Withdrawal Interval Plot'), style = "font-size:25px;"),
                          plotOutput("plot_1", width = "90%"), height = 900, width = 80, 
                          br(), tags$h1(strong(verbatimTextOutput("info")),style = "font-size:20px;"), 
                          tags$strong(verbatimTextOutput("info2")),
                          
                          downloadButton('downloadwdt','Download')
                      )
                    ) 
          ) 
        )
      ), # end of the tabItem of "plot"
      
      
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
                           
                           ## SidebarMP_1 is a function from the file 'ui_design.R'
                           conditionalPanel(condition = "input.tabs == 'Physiological parameters'&& input.animal == 'Cattle'",
                                            SidebarMP_1 (a = 'BC')
                           ), 
                           
                           conditionalPanel(condition = "input.tabs == 'Physiological parameters'&& input.animal == 'Swine'",
                                            SidebarMP_1 (a = 'SW')
                           ), 
                           
                           ## SidebarMP_2 is a function from the file 'ui_design.R'
                           conditionalPanel(condition = "input.tabs == 'Chemical-specific parameters'&& input.animal == 'Cattle'&& input.drug == 'Flunixin'",
                                            SidebarMP_2 (a = 'BC', b = 'FLU')
                           ),
                           
                           conditionalPanel(condition = "input.tabs == 'Chemical-specific parameters'&& input.animal == 'Cattle'&& input.drug == 'Florfenicol'",
                                            SidebarMP_2 (a = 'BC', b = 'FLO')
                           ),
                           
                           conditionalPanel(condition = "input.tabs == 'Chemical-specific parameters'&& input.animal == 'Cattle'&& input.drug == 'Penicillin G'",
                                            SidebarMP_2 (a = 'BC', b = 'PG')
                           ),
                           
                           conditionalPanel(condition = "input.tabs == 'Chemical-specific parameters'&& input.animal == 'Swine'&& input.drug == 'Flunixin'",
                                            SidebarMP_2 (a = 'SW', b = 'FLU')
                           ),
                           
                           conditionalPanel(condition = "input.tabs == 'Chemical-specific parameters'&& input.animal == 'Swine'&& input.drug == 'Florfenicol'",
                                            SidebarMP_2 (a = 'SW', b = 'FLO')
                           ),
                           
                           conditionalPanel(condition = "input.tabs == 'Chemical-specific parameters'&& input.animal == 'Swine'&& input.drug == 'Penicillin G'",
                                            SidebarMP_2 (a = 'SW', b = 'PG')
                           )
                           
                         )
                       )),
          
          mainPanel(width = 10,
                    tabItem(tabName = 'modelparameter',
                            fluidRow(
                              box(downloadButton('downl','Download'),
                                  title = textOutput("title1"), if(FALSE) {'FLU Concentrations in liver'}, status = "primary", solidHeader = TRUE,
                                  plotOutput(outputId = 'plot_liver', height = 350)
                              ),
                              
                              
                              box(downloadButton('downk','Download'),
                                  title = textOutput("title2"), if(FALSE) {'FLU Concentrations in kidney'},status = "primary", solidHeader = TRUE,
                                  plotOutput(outputId = 'plot_kidney', height = 350)
                              ),
                              
                              
                              box(downloadButton('downm','Download'),
                                  title = textOutput("title3"), if(FALSE){'FLU Concentrations in muscle'},status = "primary", solidHeader = TRUE,
                                  plotOutput(outputId = 'plot_muscle', height = 350)
                              ),
                              
                              box(downloadButton('downm1','Download'),
                                  title = textOutput("title4"), if(FALSE){'FLU Concentrations in plasma'},status = "primary", solidHeader = TRUE,
                                  plotOutput(outputId = 'plot_plasma', height = 350)
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
              solidHeader = TRUE, tags$img(src = 'ModStructure.png', width = '70%',height = '120%'),
              tags$h2(strong("
                   Fig. 1. A schematic diagram of the physiologically based pharmacokinetic (PBPK) model for Flunixin, 
                   Florfenicol and Penicillin G in cattle and swine. Four different administration routes
                   including intraPlasma (IV), subcuteneous (SC), intramuscular (IM) injections and oral administration are presented in the
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
)