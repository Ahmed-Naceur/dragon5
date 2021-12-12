*DECK PIJCMP
      SUBROUTINE PIJCMP(NREG,NSOUT,NCOR,DPR,VOLSUR,LPIJK,PROB)
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
* NCOR    maximum number of corners.
* DPR     collision probabilities.
* VOLSUR  volumes.
* LPIJK   pijk flag.
*
*Parameters: output
* PROB    compress probability matrix.
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
      IMPLICIT NONE
*----
* VARIABLES
*----
      INTEGER          NREG,NSOUT,NCOR,IPR,IL,JL,IVOL,IUN
      DOUBLE PRECISION DPR(-NSOUT:NREG,-NSOUT:NREG),PROB(*),ZCOR,ZCOR1
      REAL             VOLSUR(-NSOUT:NREG),COEF
      LOGICAL          LPIJK
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: WORK
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(WORK((NREG+NSOUT+2)*(NREG+NSOUT+1)/2))
*----
*  SYMMETRIZE AND STORE IN PROB
*----
      IPR= 0
      DO  150 IL    = -NSOUT, NREG
         DO  160 JL = -NSOUT, IL
            IPR= IPR+1
            WORK(IPR)= DPR(IL,JL) + DPR(JL,IL)
  160    CONTINUE
  150 CONTINUE
      IF( NCOR.EQ.1 )THEN
         IPR= 0
         DO  250 IL    = -NSOUT, NREG
            DO  260 JL = -NSOUT, IL
               IPR= IPR+1
               PROB(IPR)= WORK(IPR)
  260       CONTINUE
  250    CONTINUE
      ELSE
         IPR= 0
         ZCOR1= 1.0D0/DBLE(NCOR)
         ZCOR=  1.0D0/DBLE(NCOR*NCOR)
         DO  251 IL    = -NSOUT, NREG
            IF( IL.GT.0 ) ZCOR= ZCOR1
            DO  261 JL = -NSOUT, IL
               IPR= IPR+1
               IF( JL.GT.0 ) ZCOR= 1.0D0
               PROB(IPR)= ZCOR * WORK(IPR)
  261       CONTINUE
  251    CONTINUE
      ENDIF
*----
*  CHARGE VOLUMES IN THE PROB MATRIX
*----
      COEF=1.0
      IVOL= NSOUT*(NSOUT+1)/2
      DO 300 IUN= -NSOUT, NREG
         IF( IUN.LE.0 )THEN
            IVOL= IVOL+1
            IF(LPIJK) COEF= 3./4.
         ELSE
            IVOL= IVOL+NSOUT+IUN
            IF(LPIJK) COEF= 2./3.
         ENDIF
         IF( PROB(IVOL).NE.0.0 )THEN
            CALL XABORT( 'PIJCMP: UNEXPECTED VALUE IN PROB MATRIX' )
         ENDIF
         PROB(IVOL) = VOLSUR(IUN)*COEF
  300 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(WORK)
      RETURN
      END
