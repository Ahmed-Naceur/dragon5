*DECK DETPLAT
      SUBROUTINE DETPLAT(DT,RESPON,REF,KC,PDD,REGUL,PDF,PDT,NDDG,
     +                   AP,BP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Computes platinum detector reading with delay time consideration
* and normalized to a reference value.
*
*Copyright:
* Copyright (C) 2010 Ecole Polytechnique de Montreal.
*
*Author(s): 
* E. Varin
*
*Parameters: 
* DT       
* RESPON   
* REF      
* KC       
* PDD      
* REGUL    
* PDF      
* PDT      
* NDDG     
* AP       
* BP       
*
*-----------------------------------------------------------------------
*
      INTEGER  NDDG,I
      REAL RESPON,REF,DT
      REAL AEP,AP(NDDG),BP(NDDG),PDT(NDDG),PDF(NDDG+1),PDD(NDDG),PDO
      LOGICAL  REGUL
*
*     FLUX AT DETECTOR SITE
*     ---------------------
      IF ( KC.EQ.0) THEN
*     REFERENCE CALCULATION
*
          PDO = RESPON
          IF (REGUL) THEN
             CALL DETFIL(PDO,RESPON,RESPON,25.,DT,1)
          ENDIF
          IF (PDO.EQ.0.0) THEN
            PDO = 1.
            WRITE(6,*)'===> DETECTOR HAS AN INITIAL ZERO VALUE <==='
          ENDIF
          REF = PDO
      ENDIF
*
*     NORMALIZATION TO REFERENCE VALUE
*
      RESPON = RESPON/REF
*
      DO 10 I = 1 , NDDG
*
        AEP   = DT*PDT(I)
        AP(I) = 1.0
*
        IF (AEP.GT.20.0) THEN
           AEP = 20.0
        ENDIF
        IF (AEP.GT.1.0E-6) THEN
           AP(I) = EXP ( -AEP )
        ENDIF
        BP(I) = 1.0 - AP(I)
*
10    CONTINUE
*
      IF (KC.EQ.0) THEN
*
*     INITIALISATION AT REFERENCE CALCULATION
*     ---------------------------------------
*
         DO 20 I = 1 , NDDG
           PDD(I) = PDF(I+1) * RESPON
20       CONTINUE
      ENDIF
*
*     CALCULATION OF DETECTOR RESPONSE
*     --------------------------------
*
      PDO = 0.0
      PDO = PDF(1)*RESPON
      DO 30 I = 1 , NDDG
           PDD(I) = AP(I)*PDD(I) + BP(I)*PDF(I+1)*RESPON
           PDO    = PDO + PDD(I)
 30   CONTINUE
*
      RESPON = PDO
*
      RETURN
      END
