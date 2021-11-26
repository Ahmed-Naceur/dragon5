*DECK FMTBRN
      SUBROUTINE FMTBRN(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY,SENTRY,
     >                  IPRINT,NOPT,IOPT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Create Matlab ASCII file for burnup.
*
*Copyright:
* Copyright (C) 2017 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):
* G. Marleau
*
*Parameters: input
* NENTRY  number of data structures transfered to this module.
* HENTRY  name of the data structures.
* IENTRY  data structure type where:
*         =1 for LCM memory object;
*         =2 for XSM file;
*         =3 for sequential binary file;
*         =4 for sequential ASCII file.
* JENTRY  access permission for the data structure where:
*         =0 for a data structure in creation mode;
*         =1 for a data structure in modifications mode;
*         =2 for a data structure in read-only mode.
* KENTRY  data structure pointer.
* SENTRY  data structure signature.
* IPRINT  print level.
* NOPT    number of options.
* IOPT    processing option.
*
*-----------------------------------------------------------------------
*
      USE              GANLIB
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          NENTRY
      CHARACTER        HENTRY(NENTRY)*12
      TYPE(C_PTR)      KENTRY(NENTRY)
      INTEGER          IENTRY(NENTRY),JENTRY(NENTRY)
      CHARACTER        SENTRY(NENTRY)*12
      INTEGER          IPRINT,NOPT,IOPT(NOPT)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='FMTBRN')
      INTEGER          ILCMUP,ILCMDN
      PARAMETER       (ILCMUP=1,ILCMDN=2)
      INTEGER          NSTATE
      PARAMETER       (NSTATE=40)
*----
*  Local variables
*----
      INTEGER          IEN,IKMAT,IKEDI,IKEVO
      CHARACTER        HNEWB*12,HNEWE*12
      INTEGER          ISTATE(NSTATE)
      INTEGER          NSTEP,NBISO,NBMIXB,NGROUP,NBMIXE,MISPRT,ISTEP,
     >                 ISOT,ISOR,IKI
      REAL             TIMEFS
      CHARACTER        FMT1*35,FMT2*35
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:)   :: ISOPRT,ISOMIX
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NAMISO,NAMRD
      REAL, ALLOCATABLE, DIMENSION(:)      :: FUELDN,TIME,KEFF
      REAL, ALLOCATABLE, DIMENSION(:,:)    :: ISODEN
      IF(IPRINT .GE. 1) THEN
        WRITE(IOUT,6000) NAMSBR
        IF(IPRINT .GE. 10) THEN
          IF(NOPT.GT. 1) THEN
            WRITE(IOUT,6002) (IOPT(IKI),IKI=1,NOPT)
          ENDIF
        ENDIF
      ENDIF

*----
*  Validate entry parameters
*----
      IF(NENTRY .LT. 3) CALL XABORT(NAMSBR//
     >  ': At least three data structures required for this option.')
*----
*  Output structure (matlab ascii file)
*----
      IEN=1
      IF(IENTRY(IEN) .EQ. 4) THEN
        IF(JENTRY(IEN) .EQ. 2) CALL XABORT(NAMSBR//
     >  ': Data structure not in update or creation mode.')
        IF(JENTRY(IEN) .EQ. 1) REWIND(FILUNIT(KENTRY(IEN)))
      ELSE
        CALL XABORT(NAMSBR//
     >  ': Invalid data structure format for '//HENTRY(IEN)//'.')
      ENDIF
      IKMAT=FILUNIT(KENTRY(IEN))
*----
*  Input structure
*  Edition and burnup data structures 
*----
      IKEDI=0
      IKEVO=0
      DO IEN=2,NENTRY
        IF(IENTRY(IEN) .EQ. 1 .OR. IENTRY(IEN) .EQ. 2) THEN
          IF(JENTRY(IEN) .NE. 2) CALL XABORT(NAMSBR//
     >    ': Data structure not in read-only mode.')
          IF(SENTRY(IEN) .EQ. 'L_EDIT') THEN
            IKEDI=IEN
          ELSE IF(SENTRY(IEN) .EQ. 'L_BURNUP') THEN
            IKEVO=IEN
          ELSE
            CALL XABORT(NAMSBR//
     >      ': Invalid data structure type '//HENTRY(IEN)//'.')
          ENDIF
        ELSE
          CALL XABORT(NAMSBR//
     >    ': Invalid data structure format for '//HENTRY(IEN)//'.')
        ENDIF
      ENDDO
*----
*  Get STATE-VECTOR from BURNUP data structure
*----
      IF(IKEDI .EQ. 0) CALL XABORT(NAMSBR//
     >      ': No EDITION structure available.')
      IF(IKEVO .EQ. 0) CALL XABORT(NAMSBR//
     >      ': No BURNUP structure available.')
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(KENTRY(IKEVO),'STATE-VECTOR',ISTATE)
      NSTEP=ISTATE(3)
      NBISO=ISTATE(4)
      NBMIXB=ISTATE(8)
*----
*  Get STATE-VECTOR from EDITION data structure
*----
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(KENTRY(IKEDI),'STATE-VECTOR',ISTATE)
      NGROUP=ISTATE(1)
      NBMIXE=ISTATE(2)
*----
*  Allocate memory
*----
      ALLOCATE(FUELDN(NBMIXB),TIME(NSTEP),NAMISO(3,NBISO),
     >         ISOMIX(NBISO),NAMRD(2,NBISO),ISOPRT(NBISO))
      CALL LCMGET(KENTRY(IKEVO),'FUELDEN-MIX ',FUELDN)
      CALL LCMGET(KENTRY(IKEVO),'DEPL-TIMES  ',TIME)
      CALL LCMGET(KENTRY(IKEVO),'ISOTOPESUSED',NAMISO)
      CALL LCMGET(KENTRY(IKEVO),'ISOTOPESMIX ',ISOMIX)
*----
*  Get isotopes to to print
*----
      CALL FMTGIS(IPRINT,NBISO,NAMISO,MISPRT,NAMRD,ISOPRT)
*----
*  Allocate memory for isotopic densities and keff
*----
      ALLOCATE(KEFF(NSTEP),ISODEN(NBISO,NSTEP))
*----
*  Loop over burnup steps and extract isotopic concentration
*  keffective
*----
      DO ISTEP=1,NSTEP
        WRITE(HNEWB,'(A8,I4.4)') 'DEPL-DAT',ISTEP
        WRITE(HNEWE,'(A8,I4.4)') 'REF-CASE',ISTEP
        CALL LCMSIX(KENTRY(IKEVO),HNEWB,ILCMUP)
        CALL LCMGET(KENTRY(IKEVO),'ISOTOPESDENS',ISODEN(1,ISTEP))
        CALL LCMSIX(KENTRY(IKEVO),HNEWB,ILCMDN)
        CALL LCMSIX(KENTRY(IKEDI),HNEWE,ILCMUP)
        CALL LCMSIX(KENTRY(IKEDI),'MACROLIB    ',ILCMUP)
        CALL LCMGET(KENTRY(IKEDI),'K-EFFECTIVE ',KEFF(ISTEP))
        CALL LCMSIX(KENTRY(IKEDI),'MACROLIB    ',ILCMDN)
        CALL LCMSIX(KENTRY(IKEDI),HNEWE,ILCMDN)
      ENDDO
*----
*  Store time and keff on ascii file
*----
      TIMEFS=1E8/(24*60*60)
      FMT1=' '
      FMT2=' '
      WRITE(FMT1,'(A7,I3,A9)') '(A7,1P,',NSTEP,'E15.8,A2)'
      WRITE(IKMAT,FMT1) 'Time =[',
     >  (TIME(ISTEP)*TIMEFS,ISTEP=1,NSTEP),'];'
      WRITE(IKMAT,FMT1) 'Keff =[',
     >  (KEFF(ISTEP),ISTEP=1,NSTEP),'];'
      WRITE(FMT2,'(A17,I3,A9)') '(A3,2A4,A3,I4,A3,1P,',NSTEP,'E15.8,A2)'
      DO ISOR=1,MISPRT
        WRITE(IKMAT,'(A3,2A4,A7,I4,A1,I4,A2)') 'Den',
     >  NAMRD(1,ISOR),NAMRD(2,ISOR),
     >  '=zeros(',NSTEP,',',NBMIXB,');'
        DO ISOT=1,NBISO
          IF(ISOPRT(ISOT).EQ.ISOR) THEN
            WRITE(IKMAT,FMT2) 'Den',
     >      NAMRD(1,ISOR),NAMRD(2,ISOR),
     >      '(:,',ISOMIX(ISOT),')=[',
     >      (ISODEN(ISOT,ISTEP),ISTEP=1,NSTEP),'];'
          ENDIF
        ENDDO
      ENDDO
*----
*  Release memory and return
*----
      DEALLOCATE(ISODEN,KEFF)
      DEALLOCATE(ISOPRT,NAMRD,ISOMIX,NAMISO,TIME,FUELDN)
      IF(IPRINT .GE. 1) THEN
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Warning formats
*----
*----
*  FORMATS
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6002 FORMAT(' Options',10I10)
      END
