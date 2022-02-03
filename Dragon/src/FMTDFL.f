*DECK FMTDFL
      SUBROUTINE FMTDFL(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY,SENTRY,
     >                  IPRINT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To prepare information for the directional flux.
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
      INTEGER          IPRINT
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='FMTDFL')
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
      INTEGER          NTREG,NFUNL,NLIN,NTUNK,NTANI,ITROP,LTRK
      INTEGER          NDIM
      INTEGER          NGROUP,NREG
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KEYFLX
      REAL, ALLOCATABLE, DIMENSION(:) :: VOLUME
*----
*  Validate entry parameters
*----
      IF(NENTRY .LT. 3) CALL XABORT(NAMSBR//
     >  ': At least three data structures required for this option.')
*----
*  Output structure
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
*----
*  Input structure
*----
      IKTRK=0
      DO IEN=2,NENTRY
        IF(IENTRY(IEN) .EQ. 1 .OR. IENTRY(IEN) .EQ. 2) THEN
          IF(JENTRY(IEN) .NE. 2) CALL XABORT(NAMSBR//
     >    ': Data structure not in read-only mode.')
          IF(SENTRY(IEN) .EQ. 'L_FLUX') THEN
            IKFLU=IEN
          ELSE IF(SENTRY(IEN) .EQ. 'L_TRACK') THEN
            IKTRK=IEN
            CALL LCMGTC(KENTRY(IEN),'TRACK-TYPE',12,1,HSIGN)
            IF((HSIGN .NE. 'EXCELL').AND.(HSIGN .NE. 'MCCG')) THEN
              CALL XABORT(NAMSBR//': Only EXCELL type tracking valid.')
            ENDIF
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
      NGROUP=ISTATE(1)
      NREG=ISTATE(2)
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
      CALL LCMGET(KENTRY(IKTRK),'MCCG-STATE',ISTATE)
      NFUNL=ISTATE(19)
      NLIN=ISTATE(20)
      IF(NTUNK .NE. NREG) CALL XABORT(NAMSBR//
     >': Number of unknowns in VOLTRK and FLUX inconsistent.')
      IF(ITROP .NE. 4) CALL XABORT(NAMSBR//
     >': Only NXT: tracking permitted.')
      ALLOCATE(VOLUME(NTREG),KEYFLX(NTREG*NFUNL*NLIN))
      CALL LCMGET(KENTRY(IKTRK),'VOLUME      ',VOLUME)
      CALL LCMGET(KENTRY(IKTRK),'KEYFLX      ',KEYFLX)
      CALL LCMSIX(KENTRY(IKTRK),'NXTRecords  ',ILCMUP)
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(KENTRY(IKTRK),'G00000001DIM',ISTATE)
      NDIM=ISTATE(1)
      CALL LCMSIX(KENTRY(IKTRK),' ',ILCMDN)
*----
*  Process information
*----
      CALL FMTDFD(NENTRY,KENTRY,IPRINT,IKFLU ,NTREG ,
     >            NREG  ,NGROUP,NDIM  ,VOLUME,KEYFLX)
*----
*  Release memory and return
*----
      DEALLOCATE(VOLUME,KEYFLX)
      RETURN
      END
