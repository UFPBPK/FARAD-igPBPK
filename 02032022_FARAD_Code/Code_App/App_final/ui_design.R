# SIDEBAR_DESIGN---------------------------------------------------------------
SidebarDesign <- function(a = 'BC', b = 'FLU') {
   
  if (a == 'BC' & b == 'FLU') {
    value = c ('IV','Liver',2.2, 24, 3, 0.125, 0.125, 0.025)
  } else if (a == 'SW' & b == 'FLU') {
    value = c ('IM','Liver',2.2, 24, 1, 0.03, 0.03, 0.025)
  } else if (a == 'BC' & b == 'FLO') {
    value = c ('IM','Liver',20, 48, 2, 3.7, 3.7, 0.3)
  } else if (a == 'SW' & b == 'FLO') {
    value = c ('PO','Liver',14, 24, 5, 2.5, 2.5, 0.2)
  } else if (a == 'BC' & b == 'PG') {
    value = c ('IM','Liver',6.5, 24, 5, 0.05, 0.05, 0.05)
  } else if (a == 'SW' & b == 'PG') {
    value = c ('IM','Liver',6.5, 24, 5, 0.025, 0.025, 0.025)
  }

    
  
  div(
      selectInput(inputId = paste0("ro_", a, "_", b), label = 'Administration routes', 
                choices = c('IM', 'IV', 'SC','PO'), width = 200, selected = value[1]),
      
      selectInput(inputId = paste0("tar_", a, "_", b), label = 'Target tissue', 
                choices = c('Plasma','Liver','Kidney','Muscle'), width = 200, selected = value[2]),
      
      numericInput(inputId = paste0("dose_", a, "_", b), label = 'Dose level (mg/kg)', 
                 value = value[3], min = 0,  step = 0.01, width = 200),
      
      numericInput(inputId = paste0("int_", a, "_", b), label = 'Dose interval (h)', 
                   value = value[4], min = 1, max = 24*10, step = 1, width = 200),
      
      numericInput(inputId = paste0("ndose_", a, "_", b), label = 'Number of administrations', 
                   value = value[5], min = 1, max = 5*10, step = 1, width = 200),
      
      conditionalPanel (
        condition = paste0("input.","tar_", a, "_", b, "==", "'Liver'"),
        numericInput(inputId = 'tol_Liver', label =  'Tolerance in Liver (mg/kg or ppm)',
                 value   = value[6], min = 0.01, step = 0.01, width = 200)
        ),
      
      conditionalPanel (
        condition = paste0("input.","tar_", a, "_", b, "==", "'Kidney'"),
        numericInput(inputId = 'tol_Kidney', label =  'Tolerance in Kidney (mg/kg or ppm)',
                     value   = value[7], min = 0.01, step = 0.01, width = 200)
        ),
      
      conditionalPanel (
        condition = paste0("input.","tar_", a, "_", b, "==", "'Muscle'"),
        numericInput(inputId = 'tol_Muscle', label =  'Tolerance in Muscle (mg/kg or ppm)',
                     value   = value[8], min = 0.01, step = 0.01, width = 200)
        )
  )

}

## SIDEBAR FOR THE CONTENT OF MODEL PARAMETERS---------------------------------
SidebarMP_1 <- function(a = 'BC') {

if (a == 'BC') {
  Pars <- Pars_C_PG

} else if (a == 'SW') {
  Pars <- Pars_S_PG
  
}

div ( 
  splitLayout(
    numericInput(inputId = 'BW', label = 'BW', 
                 value = Pars[3], step = 0.01, width = 100),
    
    numericInput(inputId = 'QCC', label = 'QCC', 
                 value = Pars[5], step = 0.01, width = 100)
    ),
  
  splitLayout(
    numericInput(inputId = 'QLC', label = 'QLC', 
                 value = Pars[6], step = 0.01, width = 100),
    
    numericInput(inputId = 'QKC', label = 'QKC', 
                 value = Pars[7], step = 0.01, width = 100)
    ),
  
  splitLayout(
    numericInput(inputId = 'QMC', label = 'QMC', 
                 value = Pars[8], step = 0.01, width = 100),
    
    numericInput(inputId = 'QFC', label = 'QFC', 
                 value = Pars[9], step = 0.01, width = 100)
    ),
  
  splitLayout(
    numericInput(inputId = 'QRestC', label = 'QRestC', 
                 value = Pars[10], step = 0.01, width = 100),
  
    numericInput(inputId = 'VLC', label = 'VLC', 
                 value = Pars[11], step = 0.001, width = 100)
    ),
  
  splitLayout(
    numericInput(inputId = 'VKC', label = 'VKC', 
                 value = Pars[12], step = 0.01, width = 100),
    
    numericInput(inputId = 'VFC', label = 'VFC', 
                 value = Pars[13], step = 0.01, width = 100)
    ),
  
  splitLayout(
    numericInput(inputId = 'VMC', label = 'VMC', 
                 value = Pars[14], step = 0.01, width = 100),
    
    numericInput(inputId = 'VbloodC', label = 'VbC', 
                 value = Pars[15], step = 0.001, width = 100)
    ),
  
  splitLayout(
    numericInput(inputId = 'VRestC', label = 'VRestC', 
                 value = Pars[16], step = 0.01, width = 100)
    )
  )
 
}
  
  
SidebarMP_2 <- function(a = 'BC', b = 'PG') {
  
  if (a == 'BC') {
    if (b == 'PG') {
    Pars <- Pars_C_PG [c(1,2,17:33)]
    } else if (b == 'FLU' | b == '5-OH') {
      Pars <- Pars_C_FLU [c(1,2,17:33)]
    } else if (b == 'FLO' | b == 'FLOA') {
      Pars <- Pars_C_FLO [c(1,2,17:33)]
    }
  } else if (a == 'SW') {
    if (b == 'PG') {
      Pars <- Pars_S_PG [c(1,2,17:33)]
    } else if (b == 'FLU' | b == '5-OH') {
      Pars <- Pars_S_FLU [c(1,2,17:33)]
    } else if (b == 'FLO' | b == 'FLOA') {
      Pars <- Pars_S_FLO [c(1,2,17:33)]
    }
  }
  
  
  div(
    splitLayout(
      numericInput(inputId = 'Fracim', label = 'Fracim', 
                   value = Pars[4], step = 0.01, width = 100),
    
      numericInput(inputId = 'Fracsc', label = 'Fracsc', 
                   value = Pars[3], step = 0.01, width = 100)
      ),
    
    splitLayout(
      numericInput(inputId = 'fR', label = 'fR', 
                   value = Pars[1], step = 0.01, width = 100),
    
      numericInput(inputId = 'PL', label = 'PL', 
                   value = Pars[15], step = 0.01, width = 100)
      ),
    
    splitLayout(
      numericInput(inputId = 'PK', label = 'PK', 
                   value = Pars[18], step = 0.01, width = 100),
    
      numericInput(inputId = 'PM', label = 'PM', 
                   value = Pars[16], step = 0.01, width = 100)
      ),
    
    splitLayout(
      numericInput(inputId = 'PF', label = 'PF', 
                   value = Pars[17], step = 0.0001, width = 100),

      numericInput(inputId = 'PRest', label = 'PRest', 
                   value = Pars[19], step = 0.01, width = 100)
      ),
    
    
    splitLayout(
      numericInput(inputId = 'Kim', label = 'Kim', 
                   value = Pars[7], step = 0.001, width = 100),
    
      numericInput(inputId = 'Ksc', label = 'Ksc', 
                   value = Pars[8], step = 0.001, width = 100)
      ),
    
    splitLayout(
      numericInput(inputId = 'Kdissim', label = 'Kdissim', 
                   value = Pars[11], step = 0.01, width = 100),
    
      numericInput(inputId = 'Kdisssc', label = 'Kdisssc', 
                   value = Pars[12], step = 0.01, width = 100)
      ),
    
    splitLayout(
      numericInput(inputId = 'KehcC', label = 'KehcC', 
                   value = Pars[9], step = 0.01, width = 100),
    
      numericInput(inputId = 'KmetC', label = 'KmetC', 
                   value = Pars[10], step = 0.01, width = 100)
      ),
    
    splitLayout(
      numericInput(inputId = 'KbileC', label = 'KbileC', 
                   value = Pars[14], step = 0.01, width = 100),
    
      numericInput(inputId = 'KurineC', label = 'KurineC', 
                   value = Pars[13], step = 0.001, width = 100)
      )
  )
}
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  




