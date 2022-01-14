*DECK QIJCMP
      SUBROUTINE QIJCMP(NREG,NSOUT,NPIJ,NGRP,NCOR,VOLSUR,SIGTAL,DPR,
     > NPSYS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compression of PIJ matrices in symmetric format.
*
*Copyright:
* Copyright (C) 1994 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* NREG    number of zones for geometry.
* NSOUT   number of surfaces for geometry. 
* NPIJ    number of probabilities in one group.
* NGRP    number of groups.
* NCOR    maximum number of corners.
* VOLSUR  volumes.
* SIGTAL  materials and albedos.
* NPSYS   non-converged energy group indices.
*
*Parameters: input/output
* DPR     collision probabilities on input and
*         compress probability matrix at output.
*
*Comments:
*  Format of compress probability matrix
*    NPLEN=(NREG+NSOUT+2)*(NREG+NSOUT+1)/2
*    IND(I,J)=MAX(I+NSOUT+1,J+NSOUT+1)
*            *(MAX(I+NSOUT+1,J+NSOUT+1)-1)/2
*            +MIN(I+NSOUT+1,J+NSOUT+1)
*    IS=-NSOUT,-1; JS=-NSOUT,IS; I=IND(IS,JS)
*      PROB(I)=VOLSUR(IS)*PSS(IS,JS)
*    IV=1,NREG; JS=-NSOUT,-1;    I=IND(IV,JS)
*      SIGT(IV).GT.0.0
*        PROB(I)=SIGT(IV)*VOLSUR(IV)*PVS(IV,JS)
*      SIGT(IV).EQ.0.0
*        PROB(I)=VOLSUR(IV)*PVS(IV,JS)
*     IV=1,NREG; JV=1,IV;       I=IND(IV,JV)
*       SIGT(IV).GT.0.0 AND SIGT(JV).GT.0.0
*         PROB(I)=SIGT(IV)*SIGT(JV)*VOLSUR(IV)*PVV(IV,JV)
*       SIGT(IV).GT.0.0 AND SIGT(JV).EQ.0.0
*         PROB(I)=SIGT(IV)*VOLSUR(IV)*PVV(IV,JV)
*       SIGT(IV).EQ.0.0 AND SIGT(JV).GT.0.0
*         PROB(I)=SIGT(JV)*VOLSUR(IV)*PVV(IV,JV)
*       SIGT(IV).EQ.0.0 AND SIGT(JV).EQ.0.0
*         PROB(I)=VOLSUR(IV)*PVV(IV,JV)
*
*-----------------------------------------------------------------------
*
      IMPLICIT         NONE
      INTEGER          NREG,NSOUT,NPIJ,NGRP,NCOR,NPSYS(NGRP)
      INTEGER          IPR,IL,JL,IG,INDIJ
      REAL             VOLSUR(NSOUT:NREG),SIGTAL(NSOUT:NREG,NGRP),ZERO
      DOUBLE PRECISION DPR(NPIJ,NGRP),ZCOR,ZCOR1,DZERO
      PARAMETER      ( ZERO=0.0, DZERO=0.0D0 )
C----
C  SYMMETRIZE AND STORE IN PROB
C----
      INDIJ= 0
      DO 5 IL = 1, NREG-NSOUT+1
         INDIJ= INDIJ + IL
         DO 1 IG= 1, NGRP
           IF(NPSYS(IG).NE.0)
     >     DPR(INDIJ,IG)= DPR(INDIJ,IG) + DPR(INDIJ,IG)
    1    CONTINUE
    5 CONTINUE
      IF( NCOR.NE.1 )THEN
         ZCOR1= 1.0D0/DBLE(NCOR)
         ZCOR=  1.0D0/DBLE(NCOR*NCOR)
         INDIJ= 0
         DO 35 IL    = NSOUT, NREG
            IF( IL.GT.0 ) ZCOR= ZCOR1
            DO 25 JL = NSOUT, IL
               INDIJ= INDIJ + 1
               IF( JL.GT.0 ) ZCOR= 1.0D0
               DO 15 IG= 1, NGRP
                  IF(NPSYS(IG).NE.0)
     >            DPR(INDIJ,IG)= ZCOR * DPR(INDIJ,IG)
   15          CONTINUE
   25       CONTINUE
   35    CONTINUE
      ENDIF
      IPR=-((1-NSOUT)*NSOUT)/2
      DO 80 IL= NSOUT,-1
         IPR= IPR+1
         DO 70 IG= 1, NGRP
            IF(NPSYS(IG).NE.0) DPR(IPR,IG)= DBLE(VOLSUR(IL))
   70    CONTINUE
   80 CONTINUE
      IPR= IPR+1
      DO 90 IG= 1, NGRP
         DPR(IPR,IG)= DZERO
   90 CONTINUE
      DO 110 IL= 1,NREG
         IPR= IPR-NSOUT+IL
         DO 100 IG= 1, NGRP
            IF(NPSYS(IG).EQ.0) GO TO 100
            IF( SIGTAL(IL,IG).EQ.ZERO )THEN
               DPR(IPR,IG)= DBLE(VOLSUR(IL))
            ELSE
               DPR(IPR,IG)= DBLE(VOLSUR(IL)*SIGTAL(IL,IG))
            ENDIF
  100    CONTINUE
  110 CONTINUE
C
      RETURN
      END
