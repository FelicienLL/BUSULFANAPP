$PARAM TVCL = 5.85, TVV = 18.8
$PARAM BW_CL = 0.834, BW_V = 0.927, DAY234_CLV = 0.942
$PARAM SEX_V = 1, MALIGN_V = 1, MALIGN_CL = 1, GSTA1_CL = 1, 
$PARAM ETA1 = 0, ETA2 = 0 // IIV 
$PARAM ETA3 = 0, ETA4 = 0, ETA5 = 0, ETA6 = 0, ETA7 = 0, ETA8 = 0, ETA9 = 0, ETA10 = 0, ETA11 = 0, ETA12 = 0 // IOV 

$PARAM @annotated @covariates
BW    : 25.0 : Body weight (kg)
DAY   :  1.0 : Day (1 or more, or 0)
DIAG  : 0.0  : Diagnosis (0 malign, 1 non-malign)
SEX   : 0.0  : Sex (0 male, 1 female)
GSTA1 : 0.0  : GSTA1 variant carrier (0 no, 1 yes)
  
$OMEGA @block @name IIV
0.0678 0.0358 0.0266
$OMEGA @block @name IOV_DAY0
0.0198 0.0132 0.0192
$OMEGA @block @name IOV_DAY1
0.0198 0.0132 0.0192
$OMEGA @block @name IOV_DAY2
0.0198 0.0132 0.0192
$OMEGA @block @name IOV_DAY3
0.0198 0.0132 0.0192
$OMEGA @block @name IOV_DAY4
0.0198 0.0132 0.0192

$SIGMA 0.00371 0

$PKMODEL cmt = "CENT"
$MAIN

double MALIGN = 0 ;
if(DIAG == 0) MALIGN = 1 ; 
if(DIAG == 1) MALIGN = 0 ; 

double DAY234 = 0.0 ;
if(DAY > 1.0) DAY234 = 1.0 ; 

double DAY1 = 0 ; 
double DAY2 = 0 ; 
double DAY3 = 0 ; 
double DAY4 = 0 ; 
double DAY5 = 0 ; 
if(DAY == 1.0) DAY1 = 1 ;
if(DAY == 2.0) DAY2 = 1 ;
if(DAY == 3.0) DAY3 = 1 ;
if(DAY == 4.0) DAY4 = 1 ;
if(DAY == 5.0) DAY5 = 1 ;

double IIV_CL = ETA1 + ETA(1) ; 
double IIV_V  = ETA2 + ETA(2) ; 
double IOV_CL = DAY1 * (ETA3 + ETA(3)) + DAY2 * (ETA5 + ETA(5)) + DAY3 * (ETA7 + ETA(7)) + DAY4 * (ETA9 + ETA(9))   + DAY5 * (ETA11 + ETA(11)) ; 
double IOV_V  = DAY1 * (ETA4 + ETA(4)) + DAY2 * (ETA6 + ETA(6)) + DAY3 * (ETA8 + ETA(8)) + DAY4 * (ETA10 + ETA(10)) + DAY5 * (ETA12 + ETA(12)) ; 
double COV_CL = pow(BW/25.0, BW_CL) * pow(DAY234_CLV, DAY234) * pow(MALIGN_CL, MALIGN) * pow(GSTA1_CL, GSTA1) ; 
double COV_V  = pow(BW/25.0, BW_V)  * pow(DAY234_CLV, DAY234) * pow(MALIGN_V,  MALIGN) * pow(SEX_V,    SEX  ) ; 

double CL = TVCL * COV_CL * exp(IIV_CL + IOV_CL) ; 
double V  = TVV  * COV_V  * exp(IIV_V  + IOV_V ) ; 

$TABLE
double DV = 1000 * CENT/V * (1 + EPS(1)) + EPS(2) ; 
if(DV < 0.0) DV = 0; 
$CAPTURE DV CL COV_CL V
