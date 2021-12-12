*DECK EDIBHX
      SUBROUTINE EDIBHX (MAXPTS,IPTRK,NREG,IMERGE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Reset merging indices for the double heterogeneity option (Bihet).
*
*Copyright:
* Copyright (C) 2014 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* MAXPTS  allocated storage for arrays of dimension NREG.
* IPTRK   pointer to the tracking LCM object (L_TRACK signature).
*
*Parameters: input/output
* NREG    number of volumes in the macro geometry on input and
*         number of volumes in the composite geometry at output.
* IMERGE  merging indices in the macro geometry on input and
*         merging indices in the composite geometry at output.
*
*-----------------------------------------------------------------------
*
      USE         GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK
      INTEGER MAXPTS,NREG,IMERGE(MAXPTS)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40)
      INTEGER ISTATE(NSTATE)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IBI,NS
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: MIXGR
*
      CALL LCMSIX(IPTRK,'BIHET',1)
      CALL LCMGET(IPTRK,'PARAM',ISTATE)
      IR1=ISTATE(1)
      IR2=ISTATE(2)
      NREG2=ISTATE(3)
      NG=ISTATE(4)
      NSMAX=ISTATE(5)
      ALLOCATE(IBI(NREG2),NS(NG),MIXGR(NSMAX,NG,IR2-IR1))
      CALL LCMGET(IPTRK,'IBI',IBI)
      CALL LCMGET(IPTRK,'NS',NS)
      CALL LCMGET(IPTRK,'MIXGR',MIXGR)
      CALL LCMSIX(IPTRK,' ',2)
      NREG=NREG2
      DO 20 IKK=1,NREG2
      IF(IBI(IKK).GT.IR1) THEN
         I=IBI(IKK)-IR1
         DO 15 J=1,NG
         DO 10 K=1,NS(J)
         IF(MIXGR(K,J,I).NE.0) THEN
            NREG=NREG+1
            IMERGE(NREG)=IMERGE(IKK)
         ENDIF
   10    CONTINUE
   15    CONTINUE
      ENDIF
   20 CONTINUE
      DEALLOCATE(MIXGR,NS,IBI)
      IF(NREG.GT.MAXPTS) CALL XABORT('EDIBHX: MAXPTS IS TOO SMALL.')
      RETURN
      END
