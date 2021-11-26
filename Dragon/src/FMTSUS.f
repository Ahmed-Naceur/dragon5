*DECK FMTSUS
      SUBROUTINE FMTSUS(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY,SENTRY,
     >                  IPRINT,NOPT,IOPT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To prepare information for the SUS3D code.
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
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
      PARAMETER       (IOUT=6,NAMSBR='FMTSUS')
      INTEGER          ILCMUP,ILCMDN
      PARAMETER       (ILCMUP=1,ILCMDN=2)
      INTEGER          NSTATE
      PARAMETER       (NSTATE=40)
*----
*  Local variables
*----
      INTEGER          IEN,IKFLU,IKTRK
      CHARACTER        HSIGN*12
      INTEGER          ISTATE(NSTATE)
      INTEGER          NTREG,NTUNK,NTANI,ITROP,LTRK,
     >                 AZMOAQ,AZMQUA,POLQUA,POLOAQ
      INTEGER          NDIM
      INTEGER          ISADJ,NGROUP,NREG,NBAQU
      INTEGER          LCMLN,LCMTY
*----
*  Allocatable arrays
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: VOLUME,XPOL,WPOL,XAZI,WAZI,
     > FLUX,WGHT,MU,ETA
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: XDAZI,WDAZI,AFLUX,
     > TFLUX
*----
*  Validate entry parameters
*----
      IF(NENTRY .LT. 5) CALL XABORT(NAMSBR//
     >  ': At least five data structures required for this option.')
*----
*  Output structure
*----
      IKTRK=0
      IEN=1
      IF(IENTRY(IEN) .EQ. 4) THEN
        IF(JENTRY(IEN) .EQ. 2) CALL XABORT(NAMSBR//
     >  ': Data structure not in update or creation mode.')
        IF(JENTRY(IEN) .EQ. 1) REWIND(FILUNIT(KENTRY(IEN)))
      ELSE
        CALL XABORT(NAMSBR//
     >  ': Invalid data structure format for '//HENTRY(IEN)//'.')
      ENDIF
      DO IEN=2,3
        IF(IENTRY(IEN) .EQ. 3 .OR. IENTRY(IEN) .EQ. 4) THEN
          IF(JENTRY(IEN) .EQ. 2) CALL XABORT(NAMSBR//
     >    ': Data structure not in update or creation mode.')
          IF(JENTRY(IEN) .EQ. 1) REWIND(FILUNIT(KENTRY(IEN)))
        ELSE
          CALL XABORT(NAMSBR//
     >    ': Invalid data structure format for '//HENTRY(IEN)//'.')
        ENDIF
      ENDDO
*----
*  Input structure
*----
      DO IEN=4,NENTRY
        IF(IENTRY(IEN) .EQ. 1 .OR. IENTRY(IEN) .EQ. 2) THEN
          IF(JENTRY(IEN) .NE. 2) CALL XABORT(NAMSBR//
     >    ': Data structure not in read-only mode.')
          IF(SENTRY(IEN) .EQ. 'L_FLUX') THEN
            IKFLU=IEN
          ELSE IF(SENTRY(IEN) .EQ. 'L_TRACK') THEN
            IKTRK=IEN
            CALL LCMGTC(KENTRY(IEN),'TRACK-TYPE',12,1,HSIGN)
            IF(HSIGN .NE. 'EXCELL') CALL XABORT(NAMSBR//
     >      ': Only EXCELL type tracking valid.')
          ELSE
            CALL XABORT(NAMSBR//
     >      ': Invalid signature for '//HENTRY(IEN)//'.')
          ENDIF
        ELSE
          CALL XABORT(NAMSBR//
     >    ': Invalid data structure format for '//HENTRY(IEN)//'.')
        ENDIF
      ENDDO
*----
*  Get STATE-VECTOR from FLUX data structure
*----
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(KENTRY(IKFLU),'STATE-VECTOR',ISTATE)
      ISADJ=ISTATE(8)
      NGROUP=ISTATE(9)
      NREG=ISTATE(10)
*----
*  Get STATE-VECTOR from VOLTRK data structure
*----
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(KENTRY(IKTRK),'STATE-VECTOR',ISTATE)
      NTREG=ISTATE(1)
      NTUNK=ISTATE(2)
      NTANI=ISTATE(6)
      ITROP=ISTATE(7)
      LTRK=ISTATE(9)
      AZMOAQ=ISTATE(11)
      POLQUA=ISTATE(13)
      POLOAQ=ISTATE(14)
      AZMQUA=ISTATE(15)
      IF(NTREG .NE. NTUNK) CALL XABORT(NAMSBR//
     >': Inconsistent number of regions and unknowns.')
      IF(NTREG .NE. NREG) CALL XABORT(NAMSBR//
     >': Number of unknowns in VOLTRK and FLUX inconsistent.')
      IF(ITROP .NE. 4) CALL XABORT(NAMSBR//
     >': Only NXT: tracking permitted.')
      IF(POLQUA .LE. -1) CALL XABORT(NAMSBR//
     >': No polar quadrature provided.')
      IF(AZMQUA .LE.  0) CALL XABORT(NAMSBR//
     >': No azimuthal quadrature provided.')
      ALLOCATE(VOLUME(NTREG))
      CALL LCMGET(KENTRY(IKTRK),'VOLUME      ',VOLUME)
      CALL LCMSIX(KENTRY(IKTRK),'NXTRecords  ',ILCMUP)
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(KENTRY(IKTRK),'G00000001DIM',ISTATE)
      NDIM=ISTATE(1)
*----
*  Allocate memory for direction and weights and read them
*----
      CALL LCMLEN(KENTRY(IKTRK),'TrackingTrkW',LCMLN,LCMTY)
      AZMOAQ=LCMLN/2
      ALLOCATE(XPOL(2*POLOAQ),WPOL(POLOAQ))
      ALLOCATE(XDAZI(NDIM*AZMOAQ),WDAZI(AZMOAQ))
      CALL LCMGET(KENTRY(IKTRK),'TrackingDirc',XDAZI)
      CALL LCMGET(KENTRY(IKTRK),'TrackingTrkW',WDAZI)
      CALL LCMGET(KENTRY(IKTRK),'POLAR MU    ',XPOL)
      CALL LCMGET(KENTRY(IKTRK),'POLAR WEIGHT',WPOL)
      CALL LCMSIX(KENTRY(IKTRK),'NXTRecords  ',ILCMDN)
      ALLOCATE(XAZI(NDIM*AZMOAQ),WAZI(AZMOAQ))
      CALL XDRSDB(NDIM*AZMOAQ,XAZI,XDAZI,1)
      CALL XDRSDB(AZMOAQ,WAZI,WDAZI,1)
      DEALLOCATE(WDAZI,XDAZI)
*----
*  Allocate memory for integrated and angular flux
*---
      NBAQU=POLOAQ*AZMOAQ*2
      ALLOCATE(FLUX(NREG*2*NGROUP*2),WGHT(POLOAQ*AZMOAQ*2),
     >         MU(POLOAQ*AZMOAQ*2),ETA(POLOAQ*AZMOAQ*2))
      ALLOCATE(AFLUX(NREG*NBAQU*2*NGROUP),TFLUX((NREG+1)*NBAQU))
*----
*  Process information
*----
      CALL FMTSUD(NENTRY,IENTRY,KENTRY,IPRINT,NOPT,IOPT,IKFLU,
     >            NREG  ,NGROUP,NDIM  ,POLOAQ,AZMOAQ,
     >            VOLUME,XPOL  ,WPOL , XAZI  ,WAZI  ,
     >            FLUX  ,AFLUX ,TFLUX, WGHT  ,MU    ,ETA  )
*----
*  Release memory and return
*----
      DEALLOCATE(ETA,MU,WGHT,TFLUX,AFLUX,FLUX,WAZI,XAZI,WPOL,XPOL,
     > VOLUME)
      RETURN
*----
*  Warning formats
*----
      END
