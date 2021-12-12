*DECK EPCRMS
      SUBROUTINE EPCRMS(IPMIC,IPRINT,NGR,NXS,NMIXT,NIFISS,
     >                  NAMDXS,XSMAC)
*
*-----------------------------------------------------------------------
*
*Purpose:
*  Add perturbation to base cross section and save
*  on MACROLIB.
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
* IPMIC   pointer to MICROLIB.
* IPRINT  print level.
* NGR     number of groups.
* NXS     number of cross section types on EPC.
* NMIXT   number of mixtures on MICROLIB.
* NIFISS  number of fissile isotopes on MICROLIB.
* NAMDXS  names of XS.
* XSMAC   macroscopic vector XS contribution.
*
*-----------------------------------------------------------------------
*
      USE              GANLIB
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      TYPE(C_PTR)      IPMIC
      INTEGER          IPRINT,NGR,NXS,NMIXT,NIFISS
      CHARACTER*6      NAMDXS(NXS)
      REAL             XSMAC(NGR,NXS,NMIXT,NIFISS)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='EPCRMS')
*----
*  Local variables
*----
      TYPE(C_PTR)      JPMIC,KPMIC
      INTEGER          IPRTL,IXS,IGR,IMIX,ILCMLN,ILCMTY,ILOC,IFI
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,IPO
      REAL, ALLOCATABLE, DIMENSION(:) :: XSREC,XSSCMP
      REAL, ALLOCATABLE, DIMENSION(:,:) :: XSFIS
*----
*  Scratch storage allocation
*   XSREC   Macroscopic vector XS
*   XSFIS   Macroscopic vector fission XS
*   XSSCMP  Compressed macroscopic scattering matrix
*   IJJ     IJJ scattering index
*   IPO     IPO scattering index
*----
      ALLOCATE(IJJ(NMIXT),IPO(NMIXT))
      ALLOCATE(XSREC(NMIXT),XSFIS(NMIXT,NIFISS),XSSCMP(NGR*NMIXT))
*----
*  Write header
*----
      IPRTL=IPRINT
      IF(IPRTL .GE. 10) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      JPMIC=LCMGID(IPMIC,'GROUP')
      DO IGR=1,NGR
        KPMIC=LCMGIL(JPMIC,IGR)
        DO IXS=1,NXS
          IF(NAMDXS(IXS) .EQ. 'SIGS00') THEN
*  Total and within group scattering
            CALL LCMGET(KPMIC,'SIGS00',XSREC)
            IF(IPRTL .GE. 10) THEN
              write(6,*) NAMDXS(IXS),IGR
              WRITE(IOUT,'(1P,3E12.4)')
     >        (XSREC(IMIX),XSMAC(IGR,IXS,IMIX,1),
     >         XSREC(IMIX)+XSMAC(IGR,IXS,IMIX,1),
     >         IMIX=1,NMIXT)
            ENDIF
            DO IMIX=1,NMIXT
              IF(ABS(XSREC(IMIX)) .LT. ABS(XSMAC(IGR,IXS,IMIX,1))) THEN
               write(6,*) 'Erreur ',NAMDXS(IXS),' ',imix,igr
              ENDIF
              XSREC(IMIX)=XSREC(IMIX)+XSMAC(IGR,IXS,IMIX,1)
            ENDDO
            CALL LCMPUT(KPMIC,'SIGS00',NMIXT,2,XSREC)
            CALL LCMGET(KPMIC,'SIGW00',XSREC)
            DO IMIX=1,NMIXT
              IF(ABS(XSREC(IMIX)) .LT. ABS(XSMAC(IGR,IXS,IMIX,1))) THEN
               write(6,*) 'Erreur ',NAMDXS(IXS),' ',imix,igr
              ENDIF
              XSREC(IMIX)=XSREC(IMIX)+XSMAC(IGR,IXS,IMIX,1)
            ENDDO
            CALL LCMPUT(KPMIC,'SIGW00',NMIXT,2,XSREC)
* Scattering matrix
            CALL LCMLEN(KPMIC,'SCAT00',ILCMLN,ILCMTY)
            CALL LCMGET(KPMIC,'SCAT00',XSSCMP)
            CALL LCMGET(KPMIC,'IJJS00',IJJ)
            CALL LCMGET(KPMIC,'IPOS00',IPO)
            DO IMIX=1,NMIXT
              ILOC=IPO(IMIX)+IJJ(IMIX)-IGR
              XSSCMP(ILOC)=XSSCMP(ILOC)+XSMAC(IGR,IXS,IMIX,1)
            ENDDO
            CALL LCMPUT(KPMIC,'SCAT00',ILCMLN,ILCMTY,XSSCMP)
          ELSE IF(NAMDXS(IXS) .EQ. 'NUSIGF') THEN
*  Fission
            CALL LCMGET(KPMIC,NAMDXS(IXS),XSFIS)
            IF(IPRTL .GE. 10) THEN
              write(6,*) NAMDXS(IXS),IGR
              WRITE(IOUT,'(1P,3E12.4)')
     >       ((XSFIS(IMIX,IFI),XSMAC(IGR,IXS,IMIX,IFI),
     >        XSFIS(IMIX,IFI)+XSMAC(IGR,IXS,IMIX,IFI),
     >        IMIX=1,NMIXT),IFI=1,NIFISS)
            ENDIF
            DO IFI=1,NIFISS
              DO IMIX=1,NMIXT
              IF(ABS(XSFIS(IMIX,IFI)) .LT.
     >           ABS(XSMAC(IGR,IXS,IMIX,IFI))) THEN
               write(6,*) 'Erreur ',NAMDXS(IXS),' ',imix,igr
              ENDIF
              XSFIS(IMIX,IFI)=XSFIS(IMIX,IFI)+XSMAC(IGR,IXS,IMIX,IFI)
              ENDDO
            ENDDO
            CALL LCMPUT(KPMIC,NAMDXS(IXS),NIFISS*NMIXT,2,XSFIS)
          ELSE
            CALL LCMGET(KPMIC,NAMDXS(IXS),XSREC)
            IF(IPRTL .GE. 10) THEN
              WRITE(6,*) NAMDXS(IXS),IGR
              WRITE(IOUT,'(1P,3E12.4)')
     >        (XSREC(IMIX),XSMAC(IGR,IXS,IMIX,1),
     >         XSREC(IMIX)+XSMAC(IGR,IXS,IMIX,1),
     >         IMIX=1,NMIXT)
            ENDIF
            DO IMIX=1,NMIXT
              IF(ABS(XSREC(IMIX)) .LT. ABS(XSMAC(IGR,IXS,IMIX,1))) THEN
               write(6,*) 'Erreur ',NAMDXS(IXS),' ',imix,igr
              ENDIF
              XSREC(IMIX)=XSREC(IMIX)+XSMAC(IGR,IXS,IMIX,1)
            ENDDO
            CALL LCMPUT(KPMIC,NAMDXS(IXS),NMIXT,2,XSREC)
          ENDIF
        ENDDO
      ENDDO
*----
*  Write header and return
*----
      IF(IPRTL .GE. 10) THEN
        WRITE(IOUT,6001) NAMSBR
      ENDIF
*----
*  Scratch storage deallocation
*----
      DEALLOCATE(XSSCMP,XSFIS,XSREC)
      DEALLOCATE(IPO,IJJ)
      RETURN
*----
*  Formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
      END
