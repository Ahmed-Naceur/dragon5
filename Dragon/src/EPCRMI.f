*DECK EPCRMI
      SUBROUTINE EPCRMI(IPMIC,IPRINT,NIS,NBISO,NMIXT,NIFISS,
     >                  NAMISO,NISOU,IDVF,IDMF)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Cross reference variance isotopes and MICROLIB isotopes.
*
*Copyright:
* Copyright (C) 2009 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):
* G. Marleau
*
*Parameters: input
* IPMIC   pointer to microlib.
* IPRINT  print level.
* NIS     number of isotopes on EPC.
* NBISO   number of isotopes on MICROLIB.
* NMIXT   number of mixtures on MICROLIB.
* NIFISS  number of fissiles isotopes on MICROLIB.
*
*Parameters: output
* NAMISO  array containing the isotope names.
* NISOU   MICROLIB isotopes used.
* IDVF    variance isotopes to analyze and fission id.
* IDMF    MICROLIB isotopes to analyze and fission id.
*
*-----------------------------------------------------------------------
*
      USE              GANLIB
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      TYPE(C_PTR)      IPMIC
      INTEGER          IPRINT,NIS,NBISO,NMIXT,NIFISS
      INTEGER          NAMISO(3,NIS),NISOU(3,NBISO),
     >                 IDVF(2,NIS),IDMF(2,NBISO)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='EPCRMI')
      INTEGER          ILCMUP,ILCMDN
      PARAMETER       (ILCMUP=1,ILCMDN=2)
*----
*  Local variables
*----
      INTEGER          IPRTL,NBIU,ISO,JSO,IFI
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NISON,FID,FNM
*----
*  Scratch storage allocation
*   NISON   MICROLIB isotopes reference names
*   FID     MICROLIB fissile id
*   FNM     MICROLIB fissile name
*----
      ALLOCATE(NISON(3,NBISO),FID(NMIXT,NIFISS),FNM(2,NIFISS))
*----
*  Write header
*----
      IPRTL=IPRINT
      IF(IPRTL .GE. 10) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
*----
*  Isotope names identification
*----
      CALL LCMGET(IPMIC,'ISOTOPERNAME',NISON)
*----
*  Fissile isotopes identifier
*----
      CALL LCMSIX(IPMIC,'MACROLIB    ',ILCMUP)
      CALL LCMGET(IPMIC,'FISSIONINDEX',FID)
      CALL LCMGET(IPMIC,'FISSIONNAMES',FNM)
      CALL LCMSIX(IPMIC,'MACROLIB    ',ILCMDN)
      CALL XDISET(IDVF,2*NIS,0)
      CALL XDISET(IDMF,2*NBISO,0)
      DO ISO=1,NIS
*----
*  Test if isotope used in Microlib
*----
        NBIU=0
        DO JSO=1,NBISO
          IF( (NISON(1,JSO) .EQ. NAMISO(1,ISO)) .AND.
     >        (NISON(2,JSO) .EQ. NAMISO(2,ISO)) .AND.
     >        (NISON(3,JSO) .EQ. NAMISO(3,ISO)) ) THEN
            IDMF(1,JSO)=ISO
            NBIU=NBIU+1
          ENDIF
        ENDDO
        IF(NBIU .GT. 0) IDVF(1,ISO)=1
      ENDDO
*----
*  Find fissile isotope id
*----
      DO JSO=1,NBISO
        ISO=IDMF(1,JSO)
        IF(ISO .GT. 0) THEN
          DO IFI=1,NIFISS
            IF( (FNM(1,IFI) .EQ. NISOU(1,JSO)) .AND.
     >          (FNM(2,IFI) .EQ. NISOU(2,JSO))  ) THEN
              IDMF(2,JSO)=IFI
              IDVF(2,ISO)=IFI
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      IF(IPRTL .GE. 2) THEN
        WRITE(IOUT,6001) NAMSBR
      ENDIF
*----
*  Scratch storage deallocation
*----
      DEALLOCATE(FNM,FID,NISON)
      RETURN
*----
*  Formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
      END
