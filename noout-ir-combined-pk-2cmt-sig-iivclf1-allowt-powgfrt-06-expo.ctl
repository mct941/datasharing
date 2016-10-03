
;; -----------------------------------------------------------------------------
;; Name: /TMX-67/TMX-67_301/PK-IR/sims/noout-ir-combined-pk-2cmt-sig-iivclf1-allowt-powgfrt-06-expo.ctl
;; Wed, 15 Jun 2016 14:45:58 -0400 (sbihorel)
;; Written for use with NONMEM Version 7
;; -----------------------------------------------------------------------------
;; Last Modification (by user):
;;
;; -----------------------------------------------------------------------------
;; PURPOSE: Predict steady-state exposures based upon final parameter estimates
;;          from model ../refine/noout-ir-combined-pk-2cmt-sig-iivclf1-allowt-powgfrt-lin-06.ctl
;; -----------------------------------------------------------------------------

$PROBLEM noout-ir-combined-pk-2cmt-sig-iivclf1-allowt-powgfrt-06-expo

$INPUT ID NUM ONUM TIME TSPD CMT AMT RATE SS II 
 UDOSE EVID MDV DV TYPE DOSE STUDY AGE SEX RACE 
 ETHN SMK ALC BMI WT GFRT BUA SGPT PPI RFCAT
 BETA1 BETA2 BETA3
         
$DATA ../Data/sshell-3.csv
  IGNORE=@
  
$SUBROUTINES ADVAN6 TRANS1 TOL9

$MODEL
  COMP=(DEPOT DEFDOSE)
  COMP=(CENTRAL DEFOBS)
  COMP=(PERIPH)
  COMP=(AUC)
  
$THETA
  9.26448     ;--th1- CL: Apparent Clearance (L/h)
  42.3031     ;--th2- Vc: Apparent Central Volume (L)
  2.3417      ;--th3- Q: Apparent Distribution Clearance (L/h)
  59.3907     ;--th4- Vp: Apparent Peripheral Volume (L)
  3.59597     ;--th5- Ka: Absorption Rate Constant (1/h)
  0.594134    ;--th6- D1: Absorption Zero-order Duration (h)  
  1 FIX       ;--th7- F1: Relative Bioavailability (-)
  0.75 FIX    ;--th8- POW1: Allometric Power on CL and Q (-)
  1 FIX       ;--th9- POW2: Allometric Power on Vc and Vp (-)
  0.37142     ;--th10- CL: Power of Glomerular Filtration Effect (min/mL)
  
$OMEGA 
  0 FIX       ;--eta1- IIV in CL [exp]
  0 FIX       ;--eta2- IIV in F1 [exp]
  0 FIX       ;--eta3- IIV in Ka [exp]
   
$SIGMA 
  0.392885    ;--eps1- RV [ccv]

$PK 

  CALLFL=-2
  ; Median covariate value
  MWT = 101.3
  MGFRT = 66.4
  
  ; Covariate block
  FW1 = (WT/MWT)**THETA(8)
  FW2 = (WT/MWT)**THETA(9)
  COV1 = (GFRT/MGFRT)**THETA(10)
  
  TVCL = THETA(1)*FW1*COV1
  TVV2 = THETA(2)*FW2
  TVQ = THETA(3)*FW1
  TVV3 = THETA(4)*FW2
  TVKA = THETA(5)
  TVD1 = THETA(6)
  TVF1 = THETA(7)
  
  CL = TVCL*EXP(0*ETA(1)+BETA1)
  V2 = TVV2
  Q = TVQ
  V3 = TVV3
  KA = TVKA*EXP(0*ETA(3)+BETA3)
  D1 = TVD1
  F1 = TVF1*EXP(0*ETA(2)+BETA2)*DOSE
  
  ; AMT [mg]
  ; V2  [L]
  ; DV  [ng/mL or ug/L]
  
  S2=V2/1000
  
  K20 = CL/V2
  K23 = Q/V2
  K32 = Q/V3
  
$DES
  AON=0
  IF (TIME>=23.9) AON=1

  DADT(1) = -KA*A(1)
  DADT(2) =  KA*A(1) - (K20+K23)*A(2) + K32*A(3)
  DADT(3) = K23*A(2) - K32*A(3)

  ; type=-1 and 0 are dose records
  ; type=2 evid=2 records to obtain  Cmax,Cmin, and AUC(0-24)

  ; Definitions
  ; cmax:  max concentration between doses
  ; tmax:  time of cmax relative to time of previous dose
  ; cmin:  min concentration between doses
  ; tcmin: time of cmin relative to time of previous dose
  ; pauc:  cumulative auc for the previous dose
  ; prvd:  T (time) of previous dose
  ; cauc:  cumulative auc
  ; auct:  auctau defined as cauc-pauc
  ; cave:  auctau/tau

  ; Capture current concentration (ng/mL) and cauc (ng*hr/mL) (after last dose)
  CP = A(2)/S2
  CAUC = A(4)

  ; Initialize all variables at first dose
  ; CAUC is initialized by AUC compartment
  IF (T.LE.0) THEN
    CMAX=0
    TMAX=0
    CMIN=999999      ;number higher than any possible concentration
    CAVE=0
    PAUC=0
    PRVD=0
    AUCT=0
    CTAU=0
  ENDIF

  ; Re-initialize at other dose records (type=0)
  IF(TYPE.EQ.0)THEN
    CMAX=0
    TMAX=0
    CMIN=999999
    CAVE=0
    PAUC=CAUC
    PRVD=T
  ENDIF

  ; Compute cmax by comparing current predicted concentration to 
  ; previous value of Cmax. If higher, reset Cmax. If = or lower, retain Cmax
  ; this will find the latest highest concentration in the dose interval
  ; if there are 2 equal peaks it will assign the later one
  IF (CP.GE.CMAX) THEN
    CMAX=CP
    TMAX=T-PRVD
  ELSE
    CMAX=CMAX
    TMAX=TMAX
  ENDIF

  ; Compute cmin by comparing current predicted concentration to 
  ; previous value of Cmin. If lower, reset Cmin. If = or higher, retain Cmin
  ; this will find the earliest minimum concentration in the dose interval
  ; if there are 2 equal mins it will assign the earlier one
  IF(CP.LT.CMIN)THEN
     CMIN=CP
     TCMIN=T-PRVD
  ELSE
     CMIN=CMIN
     TCMIN=TCMIN
  ENDIF
  
  ; Calculate AUC(0-tau)
  ; current cumulative AUC - cumulative AUC of previous dose
  ; at the end of the dose interval(type=2 records) this will be AUC(0-24)
  ; AUCT and Cave are not the correct values until the end of dose interval
  ; (type=2 records)
  AUCT=CAUC-PAUC

  IF (T.NE.PRVD) THEN
    CAVE=AUCT/(T-PRVD) 
  ELSE
    CAVE=0             ;this will happen for dose records
  ENDIF;

  ; CP unit is ng/mL and time unit is h so AUC unit is ng*h/mL
  DADT(4) = AON*CP
      
$ERROR
  ; Flag dose records (they do not contribute to MVOF)
  ; the flag value will be added to the predicted to prevent
  ; the log(0) for dose records. 
  
  DFLAG=0
  IF(AMT.NE.0) DFLAG=1

  ; If a predicted concentration on sample record is 0
  ; flag so it can be set to as small a value as possible to prevent log(0)

  CFLAG=0
  CONC=A(2)*1000/V2
  IF (DFLAG.EQ.0 .AND. CONC.EQ.0) CFLAG=1
  IPRED=CONC+DFLAG+CFLAG*10**(-1*16)
  IRES=DV-IPRED
  W=IPRED
  IWRES=IRES/W
  
  Y=IPRED+W*EPS(1)
  

$SIMULATION (1234) NSUB=1 ONLYSIM

$TABLE NUM ONUM ID TIME TSPD EVID MDV IPRED CL KA F1 BETA1 BETA2 BETA3 CMAX TMAX AUCT
  CMIN CAVE STUDY DOSE AGE SEX RACE ETHN SMK ALC BMI WT GFRT BUA SGPT PPI RFCAT 
  TYPE NOPRINT ONEHEADER
  FILE=noout-ir-combined-pk-2cmt-sig-iivclf1-allowt-powgfrt-06-expo.tbl
  

