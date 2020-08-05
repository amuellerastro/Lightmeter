;added double keyword to avg and constants

function robust_mean_mod,Y,CUT,Sigma,Num_Rej,GoodInd=GoodInd
;+
; NAME:
;    Robust_Mean 
;
; PURPOSE:
;    Outlier-resistant determination of the mean and standard deviation.
;
; EXPLANATION:
;    Robust_Mean trims away outliers using the median and the median
;    absolute deviation.    An approximation formula is used to correct for
;    the trunction caused by trimming away outliers
;
; CALLING SEQUENCE:
;    mean = Robust_Mean( VECTOR, Sigma_CUT, Sigma_Mean, Num_RejECTED)
;
; INPUT ARGUMENT:
;       VECTOR    = Vector to average
;       Sigma_CUT = Data more than this number of standard deviations from the
;               median is ignored. Suggested values: 2.0 and up.
;
; OUTPUT ARGUMENT:
;       Mean  = the mean of the input vector, numeric scalar
;
; KEYWORDS:
;
;       GoodInd = The indices of the values not rejected
;
; OPTIONAL OUTPUTS:
;Sigma_Mean = the approximate standard deviation of the mean, numeric
;            scalar.  This is the Sigma of the distribution divided by sqrt(N-1)
;            where N is the number of unrejected points. The larger
;            SIGMA_CUT, the more accurate. It will tend to underestimate the
;            true uncertainty of the mean, and this may become significant for
;            cuts of 2.0 or less.
;       Num_RejECTED = the number of points trimmed, integer scalar
;
; EXAMPLE:
;       IDL> a = randomn(seed, 10000)    ;Normal distribution with 10000 pts
;       IDL> Robust_Mean,a, 3, mean, meansig, num    ;3 Sigma clipping   
;       IDL> print, mean, meansig,num
;
;       The mean should be near 0, and meansig should be near 0.01 ( =
;        1/sqrt(10000) ).    
; PROCEDURES USED:
;       AVG() - compute simple mean
; REVISION HISTORY:
;       Written, H. Freudenreich, STX, 1989; Second iteration added 5/91.
;       Use MEDIAN(/EVEN)    W. Landsman   April 2002
;       Correct conditional test, higher order truncation correction formula
;                R. Arendt/W. Landsman   June 2002
;       New truncation formula for sigma H. Freudenriech  July 2002
;-  
On_Error,2
 if N_params() LT 2 then begin
     print,'Syntax - Robust_Mean(Vector, Sigma_cut, [ Sigma_mean, '
     print,'                                  Num_Rejected ])'
     return,1
 endif  
 Npts    = N_Elements(Y)
 YMed    = MEDIAN(Y,/EVEN)
 AbsDev  = ABS(Y-YMed)
 MedAbsDev = MEDIAN(AbsDev,/EVEN)/0.6745d0
 IF MedAbsDev LT 1.0D-24 THEN MedAbsDev = AVG(AbsDev, /double)/0.8d0  
 Cutoff    = Cut*MedAbsDev  
 GoodInd = WHERE( AbsDev LE Cutoff, Num_Good )
 GoodPts = Y[ GoodInd ]
 Mean    = AVG( GoodPts, /double )
 Sigma   = SQRT( TOTAL((GoodPts-Mean)^2.d0)/Num_Good )
 Num_Rej = Npts - Num_Good ; Compenate Sigma for truncation (formula by HF):
 SC = Cut > 1.0d0
 IF SC LE 4.50d0 THEN $
    SIGMA=SIGMA/(-0.15405d0+0.90723d0*SC-0.23584d0*SC^2.d0+0.020142d0*SC^3.d0)  
 Cutoff = Cut*Sigma  
 GoodInd = WHERE( AbsDev LE Cutoff, Num_Good )
 GoodPts = Y[ GoodInd ]
 mean    = AVG( GoodPts, /double )
 Sigma   = SQRT( TOTAL((GoodPts-mean)^2.d0)/Num_Good )
 Num_Rej = Npts - Num_Good ; Fixed bug (should check for SC not Sigma) & add higher order correction
 SC = Cut > 1.0
 IF SC LE 4.50 THEN $
    SIGMA=SIGMA/(-0.15405+0.90723*SC-0.23584*SC^2.d0+0.020142*SC^3.d0) 
; Now the standard deviation of the mean:
 Sigma = Sigma/SQRT(Npts-1.0d0)  
return,mean
 END





