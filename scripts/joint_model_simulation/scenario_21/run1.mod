;; SCENARIO 21
$PROBLEM    JOINT MODEL SCENARIO 21
$INPUT      ID TIME AMT ADDL II RATE DV MDV EVID CMT FLAG AGE WT EGFR
ALB LDH ULDH SEX RACE PS TT TUS PDL GRP TRT
$DATA      small_dataset1.csv IGNORE=@
$ABBREVIATED COMRES=17
;exclusion criteria of patients
$SUBROUTINE ADVAN6 TOL=5
$MODEL      COMP=(CENTRAL,DEFDOSE,DEFOBS) COMP=(PERIPHERAL)
COMP=(EFFECT) COMP=(CHSURV) COMP=(CHTRTDUR)
$PK
;PK MODEL 
;V1SEX-DEFINITION START
IF(SEX .EQ. 0) V1SEX = 1
IF(SEX .EQ. 1) V1SEX = EXP(THETA(7))
;V1SEX-DEFINITION END
;V1-RELATION START
V1COV = V1SEX
;V1-REALATION END
;CLEGFR-DEFINITION START
CLEGFR = (EGFR/90)**THETA(8)
;CLEGFR-DEFINITION END
;CLWT-DEFINITION START
CLWT = ((WT/80)**THETA(5)) 
;CLWT-DEFINITION ENDS
;CLSEX-DEFINITION START
IF(SEX .EQ. 0) CLSEX = 1
IF(SEX .EQ. 1) CLSEX = EXP(THETA(9))
;CLSEX-DEFINITION END
;CLPS-DEFINITION START
IF(PS .EQ. 0) CLPS = 1
IF(PS .EQ. 1) CLPS = EXP(THETA(10))
;CLPS-DEFINITION END
;CLRAAS-DEFINITION START (race)
IF(RACE .EQ. 1) THEN
CLRAAS = EXP(THETA(11))
ELSE
CLRAAS = 1
ENDIF
;CLRASS-DEFINITION END
;CLALB-DEFINITION START
CLALB = (ALB/4)**THETA(12)
;CLALB-DEFINITION END
;CLBLDH-DEFINITION START (baseline LDH)
CLBLDH = (LDH/200)**THETA(13)
;CLBLDH-DEFINITION END
;CLPOPOTH-DEFINITION START
IF(TT .EQ. 1) CLPOPOTH = THETA(14)
;CLPOPOTH-DEFINITION END
;CL-RELATION START
CLCOV = CLEGFR*CLSEX*CLPS*CLRAAS*CLALB*CLBLDH*CLPOPOTH*CLWT
;CL-RELATION END
;EMAX FUNCTION FOR DEFINING TIME DEPENDANT CLEARANCE 
;EMAXPS-DEFINITION START (EFFECT OF PERFORMANCE STATUS ON EMAX)
IF(PS .EQ. 0) EMAXPS = 0
IF(PS .EQ. 1) EMAXPS = THETA(15)
;EMAXPS-DEFINITION END
;EMAX-RELATION START
EMAXCOV = EMAXPS
;EMAX-RELATION END
;EMAX PARAMETERS
TVEMAX = THETA(16)
T50 = THETA(17)
HILL = THETA(18)
TVCOVEMAX = TVEMAX + EMAXCOV
EMAX = TVCOVEMAX + ETA(4)
;EMAX FUNTION
CLTD = EXP((EMAX*(TIME**HILL))/(T50**HILL+TIME**HILL))
;TIME DEPENDENCY END
;***PK PARAMETERS****
;TTIME DEPENDANT CLEARANCE
TVCL1 = THETA(1)
TVCL = CLTD*CLCOV*TVCL1
CL = TVCL*EXP(ETA(1))
;STABLE CLEARANCE (INITIAL CLEARANCE) IN mL/h USED FOR THE HAZARD FUNCTION
CLSTAB = TVCL1*CLCOV*EXP(ETA(1))*1000
TVEXV1V2 = THETA(6)
TVV1 = THETA(2)*((WT/80)**TVEXV1V2)
TVV1 = V1COV*TVV1
V1 = TVV1*EXP(ETA(2))
TVQ = THETA(3)*((WT/80)**THETA(5))
Q = TVQ
TVV2 = THETA(4)*((WT/80)**TVEXV1V2)
V2 = TVV2*EXP(ETA(3))
S1 = V1
REP = IREP
;****TURN OF THE PK MODEL****
;IF (NEWIND .NE. 2) COUNT = 0
;F1 = 1
;IF (COUNT .GE. 1) F1 = 0 
;***** CONCENTRATION IN EFFECT COMPARTMENT****
KE0 = THETA(31)
;*****SURVIVAL MODEL******
;*****COVARIATE VALUES****
BASHAZ = THETA(19) * EXP(ETA(5))
COVSEX = THETA(20) * SEX
COVPS  = THETA(21) * PS
COVPDL = THETA(22) * PDL
COVTUS = THETA(23) * TUS
COVAGE = THETA(24) * AGE
COVALB = THETA(25) * ALB
COVLDH = THETA(26) * ULDH
COVWT  = THETA(27) * WT
COVTRT = THETA(28) * TRT
COVCONC = THETA(29) 
COVEFF = COVSEX + COVPS + COVPDL + COVTUS + COVAGE + COVALB + COVWT + COVTRT 
;*****TREATMENT DURATION MODEL*****
LAM = THETA(30)
;; THE SIMULATION PART FOR TTE SIMULATIONS ;;
IF (ICALL.EQ.4) THEN ; The event time sim $problem
IF (NEWIND.EQ.0) THEN ; Only for the first record
COM(6) = 1 ; Reset simulation ID counter
COM(4) = 26280 ; Set max time/censoring time
ENDIF
IF (NEWIND.EQ.1) THEN ; For every new ind except first in dataset
ICOUNT = COM(6) + 1 ; Update individual counter over simulations
COM(6) = ICOUNT
ENDIF
;; THE SURVIVAL SIMULATIONS
IF (NEWIND.NE.2) THEN ; For every new individual
CALL RANDOM(2,R)
COM(3) = -1 ; Variable for survival at event time
COM(2) = R ; Store the random number
COM(1) = -1 ; Variable for the event time
COM(7) = 0 ; Individual event counter
CALL RANDOM(3,R)
COM(10) = -1 ; Variable for treatment duration at event time
COM(9) = R ; Store the random number
COM(8) = -1 ; Variable for the event time
COM(11) = 0 ; Individual event counte
COM(12) = 1 ; bioavailability
COM(13) = 0 ; concentration in V1 at a time of treatment discontinuation
COM(14) = 0 ; concentartion in V1 at the time of death
COM(15) = 0 ; concentration in effect compartment  at the time of treatment discontinuation
COM(16) = 0 ; concentration in effect compartment  at the time of death
COM(17) = 0 ; concentration in effect compartment at current time
F1 = 1
ENDIF
ENDIF
;---------MTIME for increasing precision in $DES --------
IF (ICALL.EQ.4) THEN
IF (TIME.EQ.0) TEMP=0
TEMP=TEMP+24
MTIME(1)=TEMP
MTDIFF=1
ENDIF
F1=COM(12) ; the bioavailability is set to zero after the treatment is discontinued
$DES
;******PK MODEL*******
K10 = CL/V1
K12 = Q/V1
K21 = Q/V2
DADT(1) = -K10*A(1) - K12*A(1) + K21*A(2)
DADT(2) =  K12*A(1) - K21*A(2)
C1 = A(1)/V1
;**** EFFECT COMPARTMENT*****
DADT(3) = KE0*A(1) - KE0*A(3)
COM(17) = A(3)/V1 ; concentration in effect compartment
IF(GRP .EQ. 1) COM(17) = 0 ; CONCENTRATION IS ZERO IF A PATIENT IS IN CONTROL GROUP
;*****HAZARD FUNCTION******
DADT(4) = (BASHAZ*EXP(COVEFF + COVCONC*COM(17)))/24
SUR = EXP(-A(4))
IF(COM(2).GT.SUR.AND.COM(1).EQ.-1) THEN ; If event save event time in COM(1)
COM(1) = T
COM(3) = SUR
COM(14) = A(1)/V1 ; CONCENTRATION IN V1 AT THE TIME OF DEATH
COM(16) = A(3)/V1 ; CONCENTRATION IN EFFECT COMPARTMENT AT THE TIME OF DEATH
ENDIF
;*****TREATMENT DURATION FUNCTION*****
DADT(5) = LAM/24
DUR = EXP(-A(5))
IF(COM(9).GT.DUR.AND.COM(8).EQ.-1) THEN ; If event save event time in COM(8)
COM(8) = T
COM(10) = DUR 
COM(12) = 0 ; variable determening the bioavailability
COM(13) = A(1)/V1 ; CONCENTRATION IN V1 AT THE TIME OF TREATMENT DISCONTINUATION
COM(15) = A(3)/V1 ; CONCENTRATION IN EFFECT COMPARTMENT AT THE TIME OF TREATMENT DISCONTINUATION
ENDIF
$ERROR
"FIRST
"@CHARACTER(LEN=100)::FMT ! Define FORMAT string for writing dataset
IF (FLAG.EQ.2) THEN
EP1   = EPS(1)
IPRED = F
IRES  = DV - IPRED
W     = F
IWRES = IRES/W
Y = IPRED * EXP(EPS(1))
ENDIF
;******SURVIVAL MODEL******
;concentration in effect compartment
CE = A(3)/V1
IF(GRP .EQ. 1) CE = 0
;; SURVIVAL TTE MODEL
CHZ = A(4)
SURX = EXP(-CHZ)
IF (COM(1).GT.COM(4)) THEN ;IF T > ENDTIME, T=ENDTIME
; Check survival again at endtime
IF (COM(2).GT.SURX) THEN
COM(1) = COM(4)
ELSE
COM(1) = -1 ;Integrated too far, reset event
ENDIF
ENDIF
EVT = COM(1) ; Save Event time
RNM = COM(2) ; Save random number, just for debugging
ENDTIME = COM(4) ; Endtime of study
;; TREATMENT DURATION TTE MODEL
CUHAZ = A(5)
DURX = EXP(-CUHAZ)
IF (COM(8).GT.COM(4)) THEN ;IF T > ENDTIME, T=ENDTIME
; Check survival again at endtime
IF (COM(9).GT.DURX) THEN
COM(8) = COM(4)
ELSE
COM(8) = -1 ;Integrated too far, reset event
ENDIF
ENDIF
EVT1 = COM(8) ; Save Event time
RND = COM(9) ; Save random number, just for debugging
ENDTIME1 = COM(4) ; Endtime of study
; ADD RTTE, DV TO OUTPUT, SET DV=0 IF NO EVENT OR CENSORED, DV=1 IF EVENT, RTTE = 1 IF EVENT OR CENSORED
ICOUNT = COM(6)+(IREP-1)*NINDR ;NINDR The  number  of  individual records in the data set containing an observation record
ITER = IREP
; Define the format of the output file
"LAST
"FMT='(E13.7,10(1XE13.7))' ! The output FORMAT
" !Write all events
"  IF (NEWIND.EQ.0) THEN !Open file at first record
"   OPEN (99, FILE = 'simtab1.dat', POSITION='APPEND')
"       IF (IREP.EQ.1) THEN !Write header for 1st subproblem
"           WRITE (99,'(A,9(1XA))') 'ID','DV','TIME','RTTE','SURX/DURX','ICOUNT','TRT', 'CONCEFF','CONCV1', 'FLAG'
"       ENDIF
"  ENDIF
" IF (EVT1.NE.-1) THEN !If an EVENT
"   MYDV=1
"   RTTE=1
"   TMDV=0
"   COM(11) = COM(11) + 1 !Update Event counter
"   WRITE (99,FMT) ID,MYDV,EVT1,RTTE,DURX,ICOUNT,TRT,COM(15),COM(13),1.0
"   COM(8) = -1 !Reset Event time variable
"   COM(9) = 0 !Reset Random variable
"   COM(10) = -1 !Reset treatemend duration variable
"   COM(12) = 0 !bioavailability after the treatment is stoped
"   ENDIF
" IF (EVT.NE.-1) THEN !If an EVENT
"   MYDV=1
"   RTTE=1
"   TMDV=0
"   COM(7) = COM(7) + 1 !Update Event counter
"   WRITE (99,FMT) ID,MYDV,EVT,RTTE,SURX,ICOUNT,TRT,COM(16),COM(14),2.0
"   COM(1) = -1 !Reset Event time variable
"   COM(2) = 0 !Reset Random variable
"   COM(3) = -1 !Reset survival variable
"ENDIF   
" IF (LIREC.EQ.NDREC.AND.COM(7).EQ.0) THEN !Right Censoring 
"   MYDV = 0
"   TMDV = 0
"   RTTE = 1
"   TMP=COM(4)
"   WRITE (99,FMT) ID,MYDV,TMP,RTTE,SURX,ICOUNT,TRT,COM(17),C1,3.0
" ENDIF
" IF (NDREC.EQ.LIREC.AND.NIREC.EQ.NINDR) THEN !Last record for last individual
"   CLOSE(99) ! Close File pointer
" ENDIF
$THETA  
;*****PK PARAMETERS*****
0.0121 ; CL l/h
4.39 ; V1 L
0.0396 ; Q l/h
2.59 ; V2 L
0.489 ; EXCLQ BW effect on CL and Q
0.621 ; EXV1V2 BW effect on V
-0.187 ; V1SEX
0.153 ; CLEGFR
-0.181 ; CLSEX
0.126 ; CLPS performance status ecog 1 effect on CL
-0.116 ; CLRAASE race asian
-0.861 ; CLALB
0.287 ; CLBLDH BASELINE LDH
1 ; CLPOPOTH esophageal squamous cell cancer patient population
-0.109 ; EMAXPS performance status ecog 1 effect on EMAX
-0.387 ; TVEMAX maximum change in CL
1400 ; T50 time to 50pct change in CL
2.12 ; HILL coef.
;******HAZARD FUNCTIOM PARAMETERS******
0.021 ; baseline hazard?
0.0296 ; sex female vs male
-0.0976 ; performance status 1:0
-0.0932 ; pdl1 status positive vs negative
-0.06452 ; tumor size cm
0.0072 ; baseline age yr
0.1776 ; baseline albumin g/dl
0.5362 ; baseline LDH xULN
-0.0275 ; baseline weight
-1.6222 ; GROUP 1=1,GROUP 2 =0
-0.428 ; concentration in effect compartment
;*****TREATMENT DURATION*****
0.008887 ; LAMBDA value log(2)/(median treatment duration)
;*****EFFECT COMPARTMENT*****
0.0000536 ; eliminatin rate constant of the effect compartment
$OMEGA  BLOCK(2) FIX
0.0728  ;  IIV on CL
0.0352 0.0896  ; COV IIV CL:VC, IIV on VC
$OMEGA  0.261  FIX  ;  IIV on VP
0.042  FIX  ; IIV on EMAX
0  FIX  ; place holder from surv model
$SIGMA  0.219  FIX  ;         RV
$SIMULATION (154078) (739076 UNIFORM) (520709 UNIFORM) ONLYSIM NSUB=1 ;
;$TABLE      ID TIME CL Q V1 V2 ETA1 ETA2 ETA3 ETA4 EP1 IPRED PRED FLAG
;GRP FILE=pkdata.tbl NOPRINT ONEHEADER
;$TABLE      ID TIME SUR DUR CE FLAG GRP F1 COM(12)
;NOPRINT ONEHEADER FILE=survdata.tbl
