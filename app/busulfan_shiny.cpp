$PARAM TVCL = 5.85, TVV = 18.8
$PARAM BW_CL = 0.834, BW_V = 0.927, DAY234_CLV = 0.942
//$PARAM SEX_V = 1, MALIGN_V = 1, MALIGN_CL = 1, GSTA1_CL = 1, 
$PARAM ETA1 = 0, ETA2 = 0 // IIV 
$PARAM ETA3 = 0, ETA4 = 0, ETA5 = 0, ETA6 = 0, ETA7 = 0, ETA8 = 0, ETA9 = 0, ETA10 = 0, ETA11 = 0, ETA12 = 0 // IOV 

$PARAM @annotated @covariates
BW    : 25.0 : Body weight (kg)
DAY   :  1.0 : Day (0 to 4)
// SEX   :  0.0 : Sex (0 male, 1 female)
// GSTA1 :  0.0 : GSTA1 genotype (0 WT, 1 HT or HM mutation)
// MALIGN:  0.0 : Malignancy (0 no-cancer, 1 cancer)

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

double COV_CL = pow(BW/25.0, BW_CL) ; 
double COV_V  = pow(BW/25.0, BW_CL) ; 

double IIV_CL = ETA1 + ETA(1) ; 
double IIV_V  = ETA2 + ETA(2) ; 

double DAY0_CL =              exp(ETA3  + ETA(3))  ;
double DAY0_V  =              exp(ETA4  + ETA(4))  ;
double DAY1_CL =              exp(ETA5  + ETA(5))  ;
double DAY1_V  =              exp(ETA6  + ETA(6))  ;
double DAY2_CL = DAY234_CLV * exp(ETA7  + ETA(7))  ;
double DAY2_V  = DAY234_CLV * exp(ETA8  + ETA(8))  ;
double DAY3_CL = DAY234_CLV * exp(ETA9  + ETA(9))  ;
double DAY3_V  = DAY234_CLV * exp(ETA10 + ETA(10)) ;
double DAY4_CL = DAY234_CLV * exp(ETA11 + ETA(11)) ;
double DAY4_V  = DAY234_CLV * exp(ETA12 + ETA(12)) ;

double DAY0 = 0 ; 
double DAY1 = 0 ; 
double DAY2 = 0 ; 
double DAY3 = 0 ; 
double DAY4 = 0 ; 
if(DAY == 0.0) DAY0 = 1 ;
if(DAY == 1.0) DAY1 = 1 ;
if(DAY == 2.0) DAY2 = 1 ;
if(DAY == 3.0) DAY3 = 1 ;
if(DAY == 4.0) DAY4 = 1 ;

double IOV_CL = DAY0 * DAY0_CL + DAY1 * DAY1_CL + DAY2 * DAY2_CL + DAY3 * DAY3_CL + DAY4 * DAY4_CL ; 
double IOV_V  = DAY0 * DAY0_V  + DAY1 * DAY1_V  + DAY2 * DAY2_V  + DAY3 * DAY3_V  + DAY4 * DAY4_V  ; 

double CL  = TVCL * COV_CL * exp(IIV_CL) * IOV_CL  ; 
double V   = TVV  * COV_V  * exp(IIV_V ) * IOV_V   ; 

$TABLE
double DV = 1000 * CENT/V * (1 + EPS(1)) + EPS(2) ; 
if(DV < 0.0) DV = 0; 
$CAPTURE DV CL V 
