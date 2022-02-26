*DECK PIJRDG
      SUBROUTINE PIJRDG(IPRT,NREG,NSOUT,SIGTAL,PROB)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Diagonal normalization of collision probs (CP).
*
*Copyright:
* Copyright (C) 1991 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* IPRT    print level.
* NREG    number of zones for geometry.
* NSOUT   number of surfaces for geometry.
* SIGTAL  albedo-sigt vector.
*
*Parameters: input/output
* PROB    CP matrix for all types.
*
*References:
* R. Roy and G. Marleau,
* Normalization Techniques for CP Matrices,
* CONF/PHYSOR-90, Marseille/France, V 2, P IX-40 (1990).
*-----------------------------------------------------------------------
*
      IMPLICIT   NONE
      INTEGER    IPRT,NREG,NSOUT,IR,JR,IPRB,IPRF,IUNK,JUNK,IVOL
      REAL       SIGTAL(-NSOUT:NREG)
      INTEGER    IUNOUT,IPRINT
      PARAMETER (IUNOUT=6, IPRINT=4)
      DOUBLE PRECISION PROB(*),BILAN
      IPRB= 0
      IUNK= 0
      IVOL= NSOUT*(NSOUT+1)/2
*
*     RENORMALIZE ALL DIAGONAL ELEMENTS OF MATRIX *PROB*
      IF( IPRT.GE.IPRINT )THEN
         WRITE(IUNOUT,9000) 'PIJRDG'
      ENDIF
      DO 100 IR = -NSOUT, NREG
         IUNK= IUNK+1
         BILAN=0.0
         DO 10 JR= -NSOUT, IR-1
            IPRB= IPRB+1
            IF( JR.LT.0.OR.SIGTAL(JR).GT.0.0 )THEN
               BILAN=BILAN + PROB(IPRB)
            ENDIF
   10    CONTINUE
         IPRB= IPRB+1
         IPRF= IPRB
         JUNK= IUNK
         DO 20 JR=  IR+1 , NREG
            IPRF= IPRF+JUNK
            JUNK= JUNK+1
            IF( JR.LT.0.OR.SIGTAL(JR).GT.0.0 )THEN
               BILAN=BILAN + PROB(IPRF)
            ENDIF
   20    CONTINUE
         IF( IR.LT.0 )THEN
            IVOL= IVOL+1
            PROB(IPRB)= PROB(IVOL)-BILAN
            IF( IPRT.GE.IPRINT )THEN
              WRITE(IUNOUT,9001) -IR,BILAN
            ENDIF
         ELSEIF( IR.GT.0 )THEN
            IVOL= IVOL+IUNK-1
            IF( SIGTAL(IR).GT.0.0 )THEN
*
*              VOIDS ARE NOT BE RENORMALIZED
               PROB(IPRB)= PROB(IVOL)-BILAN
              IF( IPRT.GE.IPRINT )THEN
                WRITE(IUNOUT,9002) IR,BILAN
              ENDIF
            ENDIF
         ELSE
            IVOL= IVOL+1
         ENDIF
 100  CONTINUE
*
9000  FORMAT('Diagonal correction factors for CP in ',A6)
9001  FORMAT('Surface ',I10,5X,E15.6)
9002  FORMAT('Region  ',I10,5X,E15.6)
      RETURN
      END
