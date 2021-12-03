##---------------------------------------------------------------------------------------------------
# This mrgsolve  code is developed for the generic PBPK model 
# to simulate the concentrations of Flunixin (FLU), Penicillin G (PG), Florfenicol (FLO) 
# for swine and beef cattle
#----------------------------------------------------------------------------------------------------

GenricPBPK <- '
$PROB @annotated
# Generic PBPK model for Flunixin (FLU), Penicillin G (PG), Florfenicol (FLO)  
- Author    : Wei-Chun Chou
- Adivisor  : Zhoumeng Lin
- Date      : Oct, 2021
- Strucutre : GI tract, Muscle, rest of body, fat, kidneys, liver, venous and arterial blood
- Default physiological parameters were used from swine, chemical-specific parsmeters was used for FLU
- Default physiological parameters value obtained from Lin et al. (2020)
- Default chemical-specific parameters value obtained from Li et al. (2019)

//-----------------------------------------------------------------------------------------------------
$PARAM @annotated
BW          : 33.182   : Body Weight (kg)                                        , (Lin et al. 2020)
Htc         : 0.412    : Hematocrit                                              , (Lin et al. 2020)
QCC         : 8.700    : L/h/kg, Cardiac Output                                  , (Lin et al. 2020)
QLCa        : 0.273	   : unitless, Fraction of blood flow to the liver           , (Lin et al. 2020)
QKCa        : 0.114	   : Fraction of blood flow to the kidneys                   , (Lin et al. 2020)
QMCa        : 0.342	   : Fraction of blood flow to the muscle                    , (Lin et al. 2020)
QFCa        : 0.128	   : Fraction of blood flow to the fat                       , (Lin et al. 2020)
QRestCa     : 0.143    : Fraction of blood flow to the rest of body              , Value from (total sum equals to 1) (1-QLC-QKC-QFC-QMC)
VLCa        : 0.020    : Fractional liver tissue                                 , (Lin et al. 2020)
VKCa        : 0.0037   : Fractional kidney tissue                                , (Lin et al. 2020)
VFCa        : 0.154    : Fractional fat tissue                                   , (Lin et al. 2020)
VMCa        : 0.363    : Fractional muscle tissue                                , (Lin et al. 2020)
VbloodCa    : 0.041    : Total blood volume                                      , Total sum of venous and arterial blood)
VRestCa     : 0.417    : Fractional rest of body                                 , Value calculated from 1-VLC-VKC-VFC-VMC-VbloodC
PL          : 10.52    : Liver/plasma partition coefficient                      , (Li et al. 2019)   
PK          : 4.00     : Kidney/plasma partition coefficient                     , (Li et al. 2019) 
PM          : 0.50     : Muscle/plasma partition coefficient                     , (Li et al. 2019)  
PF          : 0.60     : Fat/plasma partition coefficient                        , (Li et al. 2019)  
PRest       : 8.00     : Rest of body/plasma partition coefficient               , (Li et al. 2019)  
PL1         : 10.52    : unitless,PL for metabolites                             , (Li et al. 2019) 
PK1         : 4        : unitless,PK for metabolites                             , (Li et al. 2019) 
PM1         : 0.5      : unitless,PM for metabolites                             , (Li et al. 2019)
PRest1      : 8        : unitless, PRest for metabolites, default value is assumed to be the same as the parent compound
fR          : 0.05     : unitless, Percentage of drug NOT bound to plasma proteins         
fR1         : 0.01     : unitless, Percentage of drug mtabolites NOT bound to plasma proteins
Kfeces      : 0.01     : 1/h, fecal elimination rate constant                     
Kim         : 1.00     : 1/h, Absorption rate constant of IM administration                       
Ksc         : 0.40     : 1/h, Absorption rate constant of SC administration                    
KmetC       : 0.20     : 1/h/kg, metabolic rate constant                         , 
Kmet1C      : 0        : 1/h/kg, metabolic rate constant for other metabolic pathway
KurineC     : 0.10     : L/h/kg, urine elimination rate constant for drug        ,  
KurineC1    : 0.10     : L/h/kg, urine elimination rate constant for metabolites , 
KbileC      : 0.10     : L/h/kg, biliary secretion rate constant for drug        , 
KbileC1     : 0.10     : L/h/kg, biliary secretion rate constant for metabolites , 
Kdissim     : 1e-5     : 1/h, Absorption rate constant of IM administration
Kdisssc     : 1e-5     : 1/h, Absorption rate constant of SC administration
GE          : 0.182    : 1/h, Gastric emptying time  (Yang et al., 2013)
Kabs        : 0.40     : 1/h, Rate of absorption of durg from small intestine to liver    
Kunabs      : 5e-1     : 1/h, Rate of unabsorbed dose to appear in feces 
KehcC       : 0.05     : 1/h/kg, rate constant for the enterohepatic circulation

//--------------------------------------------------------------------------------------------------------------
$MAIN
double sumQ       = QLCa + QKCa + QMCa + QFCa + QRestCa;             // Sum up cardiac output fraction
double sumV       = VLCa + VKCa + VMCa + VFCa + VbloodCa + VRestCa;  // Sum up the tissue volumes
double QLC        = QLCa/sumQ;                                       // Adjusted blood flow rate fraction to liver
double QKC        = QKCa/sumQ;                                       // Adjusted blood flow rate fraction to kidney
double QMC        = QMCa/sumQ;                                       // Adjusted blood flow rate fraction to muscle
double QFC        = QFCa/sumQ;                                       // Adjusted blood flow rate fraction to fat
double QRestC     = 1-QLC-QKC-QMC-QFC;                               // Adjusted blood flow rate fraction to rest of body
double QRest1C    = 1-QLC-QKC-QMC;                                   // Adjusted blood flow rate fraction to rest of body for metabolites
double VLC        = VLCa/sumV;                                       // Adjusted fraction of tissue volume of liver
double VKC        = VKCa/sumV;                                       // Adjusted fraction of tissue volume of kidney
double VMC        = VMCa/sumV;                                       // Adjusted fraction of tissue volume of muscle
double VFC        = VFCa/sumV;                                       // Adjusted fraction of tissue volume of fat
double VbloodC    = VbloodCa/sumV;                                   // Adjusted fraction of tissue volume of blood
double VRestC     = 1-VLC-VbloodC-VKC-VMC-VFC;                       // Adjusted fraction of tissue volume of rest of body
double VRest1C    = 1-VLC-VbloodC-VKC-VMC;                           // Adjusted fraction of tissue volume of rest of body for metabolites
double Free       = fR;                                              // Percentage of drug not bound to plasma protein for parent compound
double Free1      = fR1;                                             // Percentage of drug not bound to plasma protein for metabolites
double QC         = QCC*BW;                                          // Cardiac output
double QL         = QLC*QC;                                          // Blood flows in Liver
double QK         = QKC*QC;                                          // Blood flows in Kidney
double QF         = QFC*QC;                                          // Blood flows in Fat
double QM         = QMC*QC;                                          // Blood flows in Muscle
double QRest      = QRestC*QC;                                       // Blood flows in Rest of body (parent compound submodel)
double QRest1     = QRest1C*QC;                                      // Blood flows in Rest of body (metabolites submodel)
double VL         = VLC*BW;                                          // Tissues vloume of Liver
double VK         = VKC*BW;                                          // Tissues vloume of Kidney
double VF         = VFC*BW;                                          // Tissues vloume of Fat
double VM         = VMC*BW;                                          // Tissues vloume of Muscle
double VRest      = VRestC*BW;                                       // Tissues vloume of Rest of body (parent compound submodel)
double VRest1     = VRest1C*BW;                                      // BTissues vloume of Rest of body (metabolites submodel)
double Vblood     = VbloodC*BW;                                      // BTissues vloume of lood
double VPlas      = Vblood*(1-Htc);                                  // Tissues vloume of Plasma
double Kmet       = KmetC*BW;                                        // Metabolic rate constant for major metabolites 
double Kmet1      = Kmet1C*BW;                                       // Metabolic rate constant for other metabolites 
double Kurine     = KurineC*BW;                                      // Urinary elimination rate constant for parent compound
double Kurine1    = KurineC1*BW;                                     // Urinary elimination rate constant for metabolites
double Kbile      = KbileC*BW;                                       // Biliary elimination rate constant for parent compound
double Kbile1     = KbileC1*BW;                                      // Biliary elimination rate constant for metabolites
double Kehc       = KehcC*BW;                                        // Rate constant for enterohepatic circulation

//---------------------------------------------------------------------------------------
$INIT @annotated

ADOSE           : 0   : mg, Amount of input dose; virtual compartment
ADOSEim         : 0   : mg, Amount of input dose of intramuscular injection
ADOSEsc         : 0   : mg, Amount of input dose of subcutaneous injection
Absorbim        : 0   : mg, Amount of drug absorbed by intramuscular injection     
Amtsiteim       : 0   : mg, amount of the drug at the intramuscular injection site     
Absorbsc        : 0   : mg, Amount of drug absorbed by subcutaneous injection    
Amtsitesc       : 0   : mg, Amount of drug at the subcutaneous injection site
APlas_free      : 0   : mg, Amount of unbound drug in the plasma compartment  
APlas_free1     : 0   : mg, Amount of unbound metabolites in the plasma compartment  
AL              : 0   : mg, Amount of drug in the liver compartment    
AL1             : 0   : mg, Amount of drug metabolites in the liver compartment 
AK              : 0   : mg, Amount of drug in the kidney compartment   
AK1             : 0   : mg, Amount of drug metabolites in the kidney compartment 
AM              : 0   : mg, Amount of drug in the muscle compartment   
AM1             : 0   : mg, Amount of drug metabolites in the muscle compartment   
AF              : 0   : mg, Amount of drug in the fat compartment
ARest           : 0   : mg, Amount of drug in the rest of body compartment
ARest1          : 0   : mg, Amount of drug metabolites in the rest of body compartment
Amet            : 0   : mg, Amount of drug metabolites produced in liver
Amet1           : 0   : mg, Amount of drug metabolites produced in liver   
Aurine          : 0   : mg, Amount of drug excretion by the urine   
Aurine1         : 0   : mg, Amount of drug metabolites excretion by the urine   
Abile           : 0   : mg, Amount of drug through biliary secretion
Abile1          : 0   : mg, Amount of drug metabolites through biliary secretion
Afeces          : 0   : mg, Amount of drug excretion by feces   
AST             : 0   : mg, Amount of durg in the stomach
ASI             : 0   : mg, Amount of durg in the small intestine
AabsSI          : 0   : mg, Amount of drug absorbed by small intestine 
AUCCP           : 0   : mg/L*hr, Area under curve of drug in the plasma
AUCCP1          : 0   : mg/L*hr, Area under curve of drug metabolites in the plasma
AUCCL           : 0   : mg/L*hr, Area under curve of drug in the liver
AUCCL1          : 0   : mg/L*hr, Area under curve of drug metabolites in the liver
AUCCK           : 0   : mg/L*hr, Area under curve of drug in the kidney
AUCCK1          : 0   : mg/L*hr, Area under curve of drug metabolites in the kidney
AUCCM           : 0   : mg/L*hr, Area under curve of drug in the mammary
AUCCM1          : 0   : mg/L*hr, Area under curve of drug metabolites in the mammary

//----------------------------------------------------------------------------------------
$ODE //to define model differential equations 

// Dosing for IM and SC intramuscular 
//-- The chaning rate of the amount of dose via IM  
double Rim          = Kim*Amtsiteim;          // Rim, drug absorption rate of intramuscular route (mmol/h)
double Rsiteim      = -Rim + Kdissim*ADOSEim; // Rsiteim, changing rate of drug at the intramuscular injection site (mmol/h) 

//-- The chaning rate of the amount of dose via SC
double Rsc          = Ksc*Amtsitesc;          // Rsc, drug absorption rate of subcutaneous route (mmol /h)
double Rsitesc      = -Rsc + Kdisssc*ADOSEsc; // Rsitesc, changing rate of drug at the subcutaneous injection site (mmol/h)

//-- The amount of paraent compound in dosing compartment (virtual) via IM
dxdt_Absorbim       = Rim;                    // Absorbim, amount of drug absorbed by intramuscular injection (mmol)
dxdt_Amtsiteim      = Rsiteim;                // Amtsiteim, amount of the drug at the intramuscular injection site (mmol)
dxdt_ADOSEim        = -Kdissim*ADOSEim;

//-- The amount of paraent compound in dosing compartment (virtual) via SC
dxdt_Absorbsc       = Rsc;                    // Absorbsc, amount of drug absorbed by subcutaneous injection (mmol)
dxdt_Amtsitesc      = Rsitesc;                // Amtsitesc, amount of the drug at the subcutaneous injection site (mmol)
dxdt_ADOSEsc        = -Kdisssc*ADOSEsc;

// Concentrations of paraent/metabolized compounds in the tissues/capillary blood (mmol/L)
double CPlas_free   = APlas_free/VPlas;       // Concentration of unbound parent drug in the blood 
double CPlas_free1  = APlas_free1/VPlas;      // Concentration of unbound metabolited drug in the blood 
double CPlas        = CPlas_free/Free;        // Concentration of total drug in the blood 
double CPlas1       = CPlas_free1/Free1;      // Concentration of total drug metabolites in the plasma 
double CL           = AL/VL;                  // Concentration of parent drug in the tissue of liver 
double CL1          = AL1/VL;                 // Concentration of parent drug in the tissue of liver 
double CVL          = CL/PL;                  // Concentration of parent drug in the capillary blood of liver (mmol/L)
double CVL1         = CL1/PL1;                // Concentration of parent drug in the capillary blood of liver (mmol/L)
double CK           = AK/VK;                  // Concentration of parent drug in the tissue of kidney 
double CK1          = AK1/VK;                 // Concentration of parent drug in the tissue of kidney 
double CVK          = CK/PK;                  // Concentration of parent drug in the capillary blood of kidney (mmol/L)
double CVK1         = CK1/PK1;                // Concentration of parent drug in the capillary blood of kidney (mmol/L)
double CM           = AM/VM;                  // Concentration of parent drug in the tissue of muscle (mmol/L)
double CM1          = AM1/VM;                 // Concentration of parent drug in the tissue of muscle (mmol/L)
double CVM          = CM/PM;                  // Concentration of parent drug in the capillary blood of muscle (mmol/L)
double CVM1         = CM1/PM1;                // Concentration of parent drug in the capillary blood of muscle (mmol/L)
double CF           = AF/VF;                  // Concentration of parent drug in the tissue of fat (mmol/L)
double CVF          = CF/PF;                  // Concentration of parent drug in the capillary blood of  fat (mmol/L)
double CRest        = ARest/VRest;            // Crest drug concentration in the rest of the body (mmol/L)
double CVRest       = CRest/PRest;            // Concentration of parent drug in the capillary blood of the rest of body
double CRest1       = ARest1/VRest1;          // Crest drug concentration in the rest of the body (mmol/L)
double CVRest1      = CRest1/PRest1;          // Concentration of parent drug in the capillary blood of the rest of body
double CV           = (QL*CVL + QK*CVK + QF*CVF + QM*CVM + QRest*CVRest)/QC; // CV, the drug concentration in the blood (mmol/L)
double CV1          = (QL*CVL1 + QK*CVK1 + QM*CVM1 + QRest1*CVRest1)/QC; // CV, the drug concentration in the blood (mmol/L)


// Paraent compounds in blood compartment
//-- The chaning rate of the amount of parent/metabolized compounds in blood   
double RPlas_free   = QC*(CV-CPlas)*Free + Rsc + Rim;
double RPlas_free1  = QC*(CV1-CPlas1)*Free1;

//-- The amount and AUC of parent/metabolized compounds in blood 
dxdt_APlas_free     = RPlas_free;
dxdt_APlas_free1    = RPlas_free1;

// Paraent compounds in liver compartment
//-- The chaning rate of the amount of dose  
double Rmet         = Kmet*AL;               // Rmet, the metabolic rate from Flunixin to 5OH Flunixin in liver (mmol/h)
double Rmet1        = Kmet1*AL;         
double Rbile        = Kbile*CVL;             //Rbile, the biliary secretion rate of FLU (mmol/h)
double Rbile1       = Kbile1*CVL1;
double RL           = QL*(CPlas - CVL)*Free - Rmet - Rmet1-Rbile + Kabs*ASI + Kehc*ASI; // RL, the changing rate of the amount of FLU in liver (mmol/h)
double RL1          = QL*(CPlas1 - CVL1)*Free1 + Rmet - Rbile1;

//-- The amount and AUC of parent/metabolized compounds in liver compartment 
dxdt_Amet           = Rmet; 
dxdt_Amet1          = Rmet1;
dxdt_Abile          = Rbile;  
dxdt_Abile1         = Rbile1;
dxdt_AL             = RL; 
dxdt_AL1            = RL1;


// GI compartments and Enterohepatic Circulation of FLU
//-- The chaning rate of the amount of dose  
double RST          = -GE*AST;                                                    // mg/h,   Rate of chagne in stomach caomprtment
double RSI          = Rbile + Rbile1 + GE*AST - Kabs*ASI - Kehc*ASI - Kunabs*ASI; // mg/h,   Rate of chagne in small intestine caomprtment
double RabsSI       = Kabs*ASI;                                                   // mg/h,   Rate of absorption in the small intestine

//-- The amount and AUC of parent/metabolized compounds in liver compartment 
dxdt_AST            = RST;                                                      // mg, Amount in Stomach
dxdt_ASI            = RSI;                                                      // mg, Amount in small intestine
dxdt_AabsSI         = RabsSI;                                                   // mg, Amount absorbed in the small intestine


// Feces compartment
double Rfeces       = Kunabs*ASI;    // mg/h,   Rate of change in feces compartment
dxdt_Afeces         = Rfeces; 


// kidney compartment, flow-limited model
//-- The chaning rate of the amount of dose  
double Rurine       = Kurine*CVK;    //Rurine, the urinary secretion rate of FLU (mmol/h)
double Rurine1      = Kurine1*CVK1;    //Rurine, the urinary secretion rate of FLU (mmol/h)
double RK           = QK*(CPlas - CVK)*Free - Rurine; // RK, the changing rate of the amount of drug in kidney (mmol/h)
double RK1          = QK*(CPlas1 - CVK1)*Free1 - Rurine1;

dxdt_Aurine         = Rurine;
dxdt_Aurine1        = Rurine1;
dxdt_AK             = RK;            // AK, amount of drug in kidney (mmol)
dxdt_AK1            = RK1;

// FLU in muscle compartment, flow-limited model
double RM           = QM*(CPlas - CVM)*Free; //RM, the changing rate of the amount of drug in muscle (mmol/h)  
double RM1          = QM*(CPlas1 - CVM1)*Free1;
dxdt_AM             = RM; //AM, amount of the drug in muscle (mmol)
dxdt_AM1            = RM1;


// FLU in fat compartment, flow-limited model
double RF           = QF*(CPlas - CVF)*Free; // RF, the changing rate of the amount of drug in fat (mmol/h)
dxdt_AF             = RF;     // AF, amount of the drug in fat (mmol)

// FLU in the compartment of rest of body, flow-limited model
double Rrest        = QRest*(CPlas - CVRest)*Free; // Rrest, the changing rate of the amount of drug in the rest of the body (mmol/h)
double Rrest1       = QRest1*(CPlas1 - CVRest1)*Free1;

dxdt_ARest          = Rrest;  // Arest, amount of the drug in the rest of the body (mmol)
dxdt_ARest1         = Rrest1; // Arest, amount of the drug in the rest of the body (mmol)
dxdt_ADOSE          = 0;

// {AUC equation}
dxdt_AUCCP          = CPlas;
dxdt_AUCCP1         = CPlas1;
dxdt_AUCCL          = CL;
dxdt_AUCCL1         = CL1;
dxdt_AUCCK          = CK;
dxdt_AUCCK1         = CK1;
dxdt_AUCCM          = CM;
dxdt_AUCCM1         = CM1;


// {Mass balance equations}
double Qbal         = QC-QM-QRest-QF-QK-QL;
double Input        = Amtsiteim + Amtsitesc + ADOSEim + ADOSEsc;
double Atissue      = APlas_free + AM + ARest + AF + AK + AL + AST + ASI;
double Aloss        = Aurine + Amet + Afeces ;
double Atotal       = Atissue + Aloss + Input;
double Bal          = ADOSE  - Atotal; 
double Qbal1        = QC - QL - QK - QM - QRest1;
double Atissue1     = APlas_free1 + AM1 + AK1 + AL1 + ARest1 + Abile1 + Aurine1;
double Input1       = Amet;
double Bal1         = Input1 - Atissue1;


//------------------------------------------------------------------------------
$TABLE 
capture Plasma      = CPlas;
capture P_Met       = CPlas1;
capture Liver       = CL;
capture L_Met       = CL1;
capture Kidney      = CK;
capture K_Met       = CK1;
capture Muscle      = CM;
capture M_Met       = CM1;
capture Fat         = CF;
capture AUC_CP      = AUCCP; 
capture AUC_CP1     = AUCCP1;
capture AUC_CL      = AUCCL;
capture AUC_CL1     = AUCCL1;
capture AUC_CK      = AUCCK;
capture AUC_CK1     = AUCCK1;
capture AUC_CM      = AUCCM;
capture AUC_CM1     = AUCCM1;
'














