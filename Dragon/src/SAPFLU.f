*DECK SAPFLU
      SUBROUTINE SAPFLU(IMPX,NCALS,IPSAP,IPFLUX,IPDEPL,NGA,NRT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To recover the flux of the reference calculation.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IMPX    print parameter.
* NCALS   index of the elementary calculation.
* IPSAP   pointer to the Saphyb.
* IPFLUX  pointer to the reference flux (L_FLUX signature).
* IPDEPL  pointer to the burnup object (L_BURNUP signature).
* NGA     number of groups in the reference calculation.
* NRT     number of unknowns per group in the reference calculation.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPSAP,IPFLUX,IPDEPL
      INTEGER IMPX,NCALS,NGA,NRT
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPFLUX
      REAL BIRRAD(2)
      CHARACTER TEXT12*12,HSMG*131
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FLXREF
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(FLXREF(NRT,NGA))
*----
*  RECOVER THE FLUX NORMATIZATION FACTOR.
*----
      IF(C_ASSOCIATED(IPDEPL)) THEN
         CALL LCMGET(IPDEPL,'BURNUP-IRRAD',BIRRAD)
         BURN=BIRRAD(1)
         CALL LCMLEN(IPDEPL,'FLUX-NORM',ILONG,ITYLCM)
         IF(ILONG.EQ.0) THEN
            WRITE(HSMG,'(40HSAPFLU: THE ''FLUX-NORM'' RECORD IS NOT SE,
     1      20HT FOR BURNUP STEP AT,E12.5,14H MW-DAY/TONNE.)') BURN
            CALL XABORT(HSMG)
         ENDIF
         CALL LCMGET(IPDEPL,'FLUX-NORM',FNORM)
         IF(IMPX.GT.0) WRITE(6,100) FNORM,BURN
      ELSE
         FNORM=1.0
         IF(IMPX.GT.0) WRITE(6,110)
      ENDIF
*
      JPFLUX=LCMGID(IPFLUX,'FLUX')
      DO 20 IGR=1,NGA
      CALL LCMGDL(JPFLUX,IGR,FLXREF(1,IGR))
      DO 10 IRT=1,NRT
      FLXREF(IRT,IGR)=FLXREF(IRT,IGR)*FNORM*1.0E13
   10 CONTINUE
   20 CONTINUE
*
      WRITE(TEXT12,'(''calc'',I8)') NCALS
      CALL LCMSIX(IPSAP,TEXT12,1)
      CALL LCMSIX(IPSAP,'divers',1)
      CALL LCMPUT(IPSAP,'FLXREF',NRT*NGA,2,FLXREF)
      CALL LCMSIX(IPSAP,' ',2)
      CALL LCMSIX(IPSAP,' ',2)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(FLXREF)
      RETURN
*
  100 FORMAT(45H SAPFLU: NORMALIZE THE FLUX WITH THE FACTOR =,1P,E12.5,
     1 26H TAKEN FROM BURNUP STEP AT,E12.5,14H MW-DAY/TONNE.)
  110 FORMAT(36H SAPFLU: THE FLUX IS NOT NORMALIZED.)
      END
