##---------------------------------------------------------------------------------------------------
# This mrgsolve code is developed for the generic model for PG, FLU, FLO
# Author   : Wei-Chun Chou
# Advisor  : Zhoumeng Lin
# Date     : 2020/10/26
##----------------------------------------------------------------------------------------------------

GenricPBPK <-'
$PROB @annotated
# Generic PBPK model for Flunixin (Flu), Penicillin G (PG), Florfenicol (Flo)  
- Author    : Wei-Chun Chou
- Adivisor  : Zhoumeng Lin
- Date      : Oct, 2020
- Strucutre : GI tract, Muscle, rest of body, fat, kidneys, liver, venous and arterial blood
- Default physiological parameters was used from swine, chemical-specific parsmeters was used for Flu 

//-----------------------------------------------------------------------------------------------------
$PARAM @annotated
BW          : 33.182   : Body Weight (kg)                                        , value from Li et al. (2017)
Htc         : 0.35     : Hematocrit                                              , value from Estienne et al, 2019
QCC         : 8.543    : L/h/kg, Cardiac Output                                  , value from Li et al. (2017)
QLCa        : 0.273	   : Fraction of blood flow to the liver                     , value from Li et al. (2017)
QKCa        : 0.116	   : Fraction of blood flow to the kidneys                   , value from Li et al. (2017)
QMCa        : 0.293	   : Fraction of blood flow to the muscle                    , value from Li et al. (2017)
QFCa        : 0.128	   : Fraction of blood flow to the fat                       , value from Li et al. (2017)
QRestCa     : 0.190    : Fraction of blood flow to the rest of body              , value from (total sum equals to 1) (1-QLC-QKC-QFC-QMC)
VLCa        : 0.023    : Fractional liver tissue                                 , (Li et al. 2017)
VKCa        : 0.0045   : Fractional kidney tissue                                , (Li et al. 2017)
VFCa        : 0.235    : Fractional fat tissue                                   , (Li et al. 2017)
VMCa        : 0.355    : Fractional muscle tissue                                , (Li et al. 2017)
VbloodCa    : 0.050    : Total blood volume                                      , total sum of venous and arterial blood)
VRestCa     : 0.322    : Fractional rest of body                                 , value calculated from 1-VLC-VKC-VFC-VMC-VbloodC
PL          : 10.52    : Liver/plasma partition coefficient                      , model fitting   
PK          : 4.00     : Kidney/plasma partition coefficient                     , model fitting 
PM          : 0.50     : Muscle/plasma partition coefficient                     , model fitting 
PF          : 0.60     : Fat/plasma partition coefficient                        , model fitting 
PRest       : 8.00     : Rest of body/plasma partition coefficient               , model fitting 
PL1         : 10.52    : PL for metabolites                                      , default value is assumed to be the same as the parent compound
PK1         : 4        : PK for metabolites                                      , default value is assumed to be the same as the parent compound
PM1         : 0.5      : PM for metabolites                                      , default value is assumed to be the same as the parent compound
PRest1      : 8        : PRest for metabolites                                   , default value is assumed to be the same as the parent compound
Fub         : 0.95     : Percentage of drug bound to plasma proteins             , 
Fub1        : 0.99     : Percentage of drug mtabolites bound to plasma proteins  ,
Kfeces      : 0.01     : 1/h, fecal elimination rate constant                    , model fitting 
Kim         : 1.00     : 1/h, IM absorption rate constant                        , model fitting 
Ksc         : 0.40     : 1/h, SC absorption rate constant                        , model fitting 
KmetC       : 0.20     : 1/h/kg, metabolic rate constant                         , model fitting
Kmet1C      : 0        :
KurineC     : 0.10     : L/h/kg, urine elimination rate constant for drug        , model fitting 
KurineC1    : 0.10     : L/h/kg, urine elimination rate constant for metabolites , model fitting
KbileC      : 0.10     : L/h/kg, biliary secretion rate constant for drug        , model fitting
KbileC1     : 0.10     : L/h/kg, biliary secretion rate constant for metabolites , model fitting
Kdissim     : 1e-5     : 1/h, 
Kdisssc     : 1e-5     : 1/h, 
GE          : 0.182    : 1/h,  Gastric emptying time  (Yang et al., 2013)
Kabs        : 1.91     : 1/(h, Rate of absorption of durg from small intestine to liver    
Kunabs      : 5e-5     : 1/(h*BW^-0.25), Rate of unabsorbed dose to appear in feces (initial value assumed the same as PFOA (7.05e-5) from Worley and Fisher, 2015 and then re-fitting)
KehcC       : 0.05     : rate constant for the enterohepatic circulation


//-----------------------------------------------------------------------------------------------
$MAIN

double sumQ       = QLCa + QKCa + QMCa + QFCa + QRestCa;             // sum up cardiac output fraction
double sumV       = VLCa + VKCa + VMCa + VFCa + VbloodCa + VRestCa; // sum up the tissue volumes
double QLC        = QLCa/sumQ;        // adjusted blood flow rate fraction to liver
double QKC        = QKCa/sumQ;        // adjusted blood flow rate fraction to kidney
double QMC        = QMCa/sumQ;        // adjusted blood flow rate fraction to muscle
double QFC        = QFCa/sumQ;        // adjusted blood flow rate fraction to fat
double QRestC     = QRestCa/sumQ;     // adjusted blood flow rate fraction to rest of body
double QRest1C    = 1-QLC-QKC-QMC;     // adjusted blood flow rate fraction to rest of body

double VLC        = VLCa/sumV;       // adjusted fraction of tissue volume of liver
double VKC        = VKCa/sumV;       // adjusted fraction of tissue volume of kidney
double VMC        = VMCa/sumV;       // adjusted fraction of tissue volume of muscle
double VFC        = VFCa/sumV;       // adjusted fraction of tissue volume of fat
double VbloodC    = VbloodCa/sumV;
double VRestC     = VRestCa/sumV;    // adjusted fraction of tissue volume of rest of body
double VRest1C    = 1-VLC-VbloodC-VKC-VMC;
double Free       = 1-Fub;            // Percentage of drug not bound to plasma protein for paranet compound
double Free1      = 1-Fub1;            // Percentage of drug not bound to plasma protein for paranet compound

double QC         = QCC*BW;       // Cardiac output
double QL         = QLC*QC;       // Liver
double QK         = QKC*QC;       // Kidney
double QF         = QFC*QC;       // Fat
double QM         = QMC*QC;       // Muscle
double QRest      = QRestC*QC;    // Rest of body
double QRest1     = QRest1C*QC;    // Rest of body

double VL         = VLC*BW;     // Liver
double VK         = VKC*BW;     // Kidney
double VF         = VFC*BW;     // Fat
double VM         = VMC*BW;     // Muscle
double VRest      = VRestC*BW;          // Rest of body
double VRest1     = VRest1C*BW;          // Rest of body

double Vblood     = VbloodC*BW;         // Blood
double VPlas      = Vblood*(1-Htc);    // Plasma

double Kmet       = KmetC*BW;    // Metabolic rate constant from FLU to 5OH FLU
double Kmet1      = Kmet1C*BW;
double Kurine     = KurineC*BW;  // Urinary elimination rate constant for FLU
double Kurine1    = KurineC1*BW;
double Kbile      = KbileC*BW;   // Biliary elimination rate constant for FLU
double Kbile1     = KbileC1*BW;
double Kehc       = KehcC*BW;


//---------------------------------------------------------------------------------------
$INIT @annotated

ADOSE           : 0   : mg, Amount of input dose; virtual compartment
ADOSEim         : 0   : mg, Amount of input dose of intramuscular injection
ADOSEsc         : 0   : mg, Amount of input dose of subcutaneous injection
Absorbim        : 0   : mg, Amount of drug absorbed by intramuscular injection     
Amtsiteim       : 0   : mg, amount of the drug at the intramuscular injection site     
Absorbsc        : 0   : mg, Amount of drug absorbed by subcutaneous injection    
Amtsitesc       : 0   : mg, Amount of drug at the subcutaneous injection site
APlas           : 0   : mg, Amount of unbound drug in the plasma compartment  
APlas1          : 0   : mg, Amount of unbound metabolites in the plasma compartment  
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
$ODE // {to define Ordinary Differential Equations} 

// Dosing for IM and SC route
//-- The chaning rate of the amount of dose via IM  
double Rim          = Kim*Amtsiteim;           // Rim, drug absorption rate of intramuscular route (mmol/h)
double Rsiteim      = -Rim + Kdissim*ADOSEim;  // Rsiteim, changing rate of drug at the intramuscular injection site (mmol/h) 

//-- The chaning rate of the amount of dose via SC
double Rsc          = Ksc*Amtsitesc;           // Rsc, drug absorption rate of subcutaneous route (mmol /h)
double Rsitesc      = -Rsc + Kdisssc*ADOSEsc;  // Rsitesc, changing rate of drug at the subcutaneous injection site (mmol/h)

//-- The amount of paraent compound in dosing compartment (virtual) via IM
dxdt_Absorbim       = Rim;                     // Absorbim, amount of drug absorbed by intramuscular injection (mmol)
dxdt_Amtsiteim      = Rsiteim;                 // Amtsiteim, amount of the drug at the intramuscular injection site (mmol)
dxdt_ADOSEim        = -Kdissim*ADOSEim;

//-- The amount of paraent compound in dosing compartment (virtual) via SC
dxdt_Absorbsc       = Rsc;                     // Absorbsc, amount of drug absorbed by subcutaneous injection (mmol)
dxdt_Amtsitesc      = Rsitesc;                 // Amtsitesc, amount of the drug at the subcutaneous injection site (mmol)
dxdt_ADOSEsc        = -Kdisssc*ADOSEsc;

// Concentrations of paraent/metabolized compounds in the tissues/capillary blood (mmol/L)
double CPlas        = APlas/VPlas;           // Concentration of total drug in the blood 
double CPlas1       = APlas1/VPlas;           // Concentration of total drug metabolites in the plasma 
double CPlas_free   = CPlas*Free;             // Concentration of unbound drug in the blood 
double CPlas_free1  = CPlas1*Free1;
double CL           = AL/VL;               // Concentration of parent drug in the tissue of liver 
double CL1          = AL1/VL;               // Concentration of parent drug in the tissue of liver 
double CVL          = CL/PL;               // Concentration of parent drug in the capillary blood of liver (mmol/L)
double CVL1         = CL1/PL1;               // Concentration of parent drug in the capillary blood of liver (mmol/L)
double CK           = AK/VK;               // Concentration of parent drug in the tissue of kidney 
double CK1          = AK1/VK;               // Concentration of parent drug in the tissue of kidney 
double CVK          = CK/PK;               // Concentration of parent drug in the capillary blood of kidney (mmol/L)
double CVK1         = CK1/PK1;               // Concentration of parent drug in the capillary blood of kidney (mmol/L)
double CM           = AM/VM;               // Concentration of parent drug in the tissue of muscle (mmol/L)
double CM1          = AM1/VM;               // Concentration of parent drug in the tissue of muscle (mmol/L)
double CVM          = CM/PM;               // Concentration of parent drug in the capillary blood of muscle (mmol/L)
double CVM1         = CM1/PM1;               // Concentration of parent drug in the capillary blood of muscle (mmol/L)
double CF           = AF/VF;               // Concentration of parent drug in the tissue of fat (mmol/L)
double CVF          = CF/PF;               // Concentration of parent drug in the capillary blood of  fat (mmol/L)
double CRest        = ARest/VRest;         // Crest drug concentration in the rest of the body (mmol/L)
double CVRest       = CRest/PRest;         // Concentration of parent drug in the capillary blood of the rest of body
double CRest1       = ARest1/VRest1;         // Crest drug concentration in the rest of the body (mmol/L)
double CVRest1      = CRest1/PRest1;         // Concentration of parent drug in the capillary blood of the rest of body
double CV           = (QL*CVL + QK*CVK + QF*CVF + QM*CVM + QRest*CVRest + Rsc + Rim)/QC; // CV, the drug concentration in the blood (mmol/L)
double CV1          = (QL*CVL1 + QK*CVK1 + QM*CVM1 + QRest1*CVRest1)/QC; // CV, the drug concentration in the blood (mmol/L)


// Paraent compounds in blood compartment
//-- The chaning rate of the amount of parent/metabolized compounds in blood   
double RPlas        = QC*(CV-CPlas_free) ;
double RPlas1       = QC*(CV1-CPlas_free1);

//-- The amount and AUC of parent/metabolized compounds in blood 
dxdt_APlas          = RPlas;
dxdt_APlas1         = RPlas1;

// Paraent compounds in liver compartment
//-- The chaning rate of the amount of dose  
double Rmet         = Kmet*AL;     // Rmet, the metabolic rate from Flunixin to 5OH Flunixin in liver (mmol/h)
double Rmet1        = Kmet1*AL;
double Rbile        = Kbile*CVL;   //Rbile, the biliary secretion rate of FLU (mmol/h)
double Rbile1       = Kbile1*CVL1;
double RL           = QL*(CPlas_free - CVL) - Rmet - Rmet1- Rbile + Kabs*ASI + Kehc*ASI; // RL, the changing rate of the amount of FLU in liver (mmol/h)
double RL1          = QL*(CPlas_free1 - CVL1) + Rmet - Rbile1;

//-- The amount and AUC of parent/metabolized compounds in liver compartment 
dxdt_Amet           = Rmet;
dxdt_Amet1          = Rmet1;
dxdt_Abile          = Rbile;  
dxdt_Abile1         = Rbile1;
dxdt_AL             = RL; 
dxdt_AL1            = RL1;


// GI compartments and Enterohepatic Circulation of FLU
//-- The chaning rate of the amount of dose  
double RST          = -GE*AST;                                                             // mg/h,   Rate of chagne in stomach caomprtment
double RSI          = Rbile + Rbile1 + GE*AST - Kabs*ASI - Kehc*ASI - Kunabs*ASI;                                                // mg/h,   Rate of chagne in small intestine caomprtment
double RabsSI       = Kabs*ASI;                                                                   // mg/h,   Rate of absorption in the small intestine

//-- The amount and AUC of parent/metabolized compounds in liver compartment 
dxdt_AST            = RST; // mg, Amount in Stomach
dxdt_ASI            = RSI;                                                                             // mg,     Amount in small intestine
dxdt_AabsSI         = RabsSI;                                                                       // mg,     Amount absorbed in the small intestine


// Feces compartment
double Rfeces       = Kunabs*ASI;    // mg/h,   Rate of change in feces compartment
dxdt_Afeces         = Rfeces; 


// kidney compartment, flow-limited model
//-- The chaning rate of the amount of dose  
double Rurine       = Kurine*CVK;    //Rurine, the urinary secretion rate of FLU (mmol/h)
double Rurine1      = Kurine1*CVK1;    //Rurine, the urinary secretion rate of FLU (mmol/h)
double RK           = QK*(CPlas_free - CVK) - Rurine; // RK, the changing rate of the amount of drug in kidney (mmol/h)
double RK1          = QK*(CPlas_free1 - CVK1) - Rurine1;

dxdt_Aurine         = Rurine;
dxdt_Aurine1        = Rurine1;
dxdt_AK             = RK;            // AK, amount of drug in kidney (mmol)
dxdt_AK1            = RK1;

// FLU in muscle compartment, flow-limited model
double RM           = QM*(CPlas_free - CVM); //RM, the changing rate of the amount of drug in muscle (mmol/h)  
double RM1          = QM*(CPlas_free1 - CVM1);
dxdt_AM             = RM; //AM, amount of the drug in muscle (mmol)
dxdt_AM1            = RM1;


// FLU in fat compartment, flow-limited model
double RF           = QF*(CPlas_free - CVF); // RF, the changing rate of the amount of drug in fat (mmol/h)
dxdt_AF             = RF; // AF, amount of the drug in fat (mmol)

// FLU in the compartment of rest of body, flow-limited model
double Rrest        = QRest*(CPlas_free - CVRest); // Rrest, the changing rate of the amount of drug in the rest of the body (mmol/h)
double Rrest1       = QRest1*(CPlas_free1 - CVRest1);

dxdt_ARest          = Rrest; // Arest, amount of the drug in the rest of the body (mmol)
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
double Atissue      = APlas + AM + ARest + AF + AK + AL + AST + ASI;
double Aloss        = Aurine + Amet + Afeces ;
double Atotal       = Atissue + Aloss + Input;
double Bal          = ADOSE  - Atotal; 

double Qbal1        = QC - QL - QK - QM - QRest1;
double Atissue1     = APlas1 + AM1 + AK1 + AL1 + ARest1 + Abile1 + Aurine1;
double Input1       = Amet;
double Bal1         = Input1 - Atissue1;


// {Output table}
$TABLE 
capture Plasma      = CPlas;
capture P_Met       = CPlas_free1;
capture Liver       = CL;
capture L_Met       = CL1;
capture Kidney      = CK;
capture K_Met       = CK1;
capture Muscle      = CM;
capture M_Met       = CM1;
capture Fat         = CF;
capture Mbal        = Bal;
capture Mbal1       = Bal1;
capture AUC_CP      = AUCCP; 
capture AUC_CP1     = AUCCP1;
capture AUC_CL      = AUCCL;
capture AUC_CL1     = AUCCL1;
capture AUC_CK      = AUCCK;
capture AUC_CK1     = AUCCK1;
capture AUC_CM      = AUCCM;
capture AUC_CM1     = AUCCM1;
'













