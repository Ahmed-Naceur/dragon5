*DECK PSP
      SUBROUTINE PSP(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* POSTSCRIPT plot utility module.
*
*Copyright:
* Copyright (C) 1999 Ecole Polytechnique de Montreal
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
*
*Comments:
*  Input requirements
*    NENTRY  >= 2
*      IEN = 1 : structure is a sequential ascii file
*                containing the output POSTSCRIPT. this can be A
*                new file or a file to update.
*                IENTRY(ien) = 4, JENTRY(ien) <= 1
*    FOR PSP BY Mixture or region
*      IEN > 1 : structure is a valid dragon geometry
*                for excelt stored in a linked list or XSM file.
*                this structure must be in read-only mode
*                IENTRY(ien)<= 2, JENTRY(ien) = 2
*    FOR PSP BY Flux
*      IEN = 2 : structure is a valid dragon geometry
*                for excelt stored in a linked list or XSM file.
*                this structure must be in read-only mode
*                IENTRY(ien)<= 2, JENTRY(ien) = 2
*      IEN = 3 : structure is a valid flux structure
*
*-----------------------------------------------------------------------
*
      USE          GANLIB
      IMPLICIT     NONE
      INTEGER      IOUT,NSTATE,ILCMUP,ILCMDN
      CHARACTER    NAMSBR*6
      PARAMETER   (IOUT=6,NSTATE=40,ILCMUP=1,ILCMDN=2,
     >             NAMSBR='PSP   ')
*----
*  ROUTINE PARAMTERS
*----
      INTEGER      NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      TYPE(C_PTR)  KENTRY(NENTRY)
      CHARACTER    HENTRY(NENTRY)*12
*----
*  LOCAL VARIABLES
*----
      INTEGER      ISTATE(NSTATE),IPRINT,IEN,
     >             IENGT,IENFL,ISPSP,ITYPE,ICOLR,
     >             NPAGE,ITRK,NGROUP,NUNKNO,IGR,NGT
      CHARACTER    HSIGN*12,NAMGT*12,NAMLEG*24
      REAL         XYPOS(2)
      CHARACTER    NAMTR2*12,NAMGEO*12
      TYPE(C_PTR)  IPTRK2,IPFL,IPGT
      INTEGER      IMODT2,IMEDT2,ICLST2,IPRIN2
      LOGICAL      LASS,LDRASS
      INTEGER      IMODE,NMODE
      TYPE(C_PTR)  JPFL,KPFL
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ICOND
      REAL, ALLOCATABLE, DIMENSION(:) :: FLUX,TFLX
*----
      NAMTR2='PSPGEOIPTRK2'
      IMODT2=0
      IMEDT2=1
      IPRIN2=0
      ICLST2=2
*----
*  INPUT PARAMETER VALIDATION
*----
      XYPOS(1)=0.5
      XYPOS(2)=0.5
      IF(NENTRY    .LT. 2 ) CALL XABORT(NAMSBR//
     >  ': AT LEAST TWO DATA STRUCTURES REQUIRED')
      ISPSP=FILUNIT(KENTRY(1))
      IF(IENTRY(1) .NE. 4 ) CALL XABORT(NAMSBR//
     >  ': POSTSCRIPT DATA STRUCTURE NOT AN ASCII FILE')
      IF(JENTRY(1) .NE. 0   .AND.
     >   JENTRY(1) .NE. 1 ) CALL XABORT(NAMSBR//
     >  ': POSTSCRIPT DATA STRUCTURE NOT IN CREATE OR MODIFY MODE')
*----
*  Find if one of the structures is a flux
*  and get number of groups if required
*----
      NGROUP=1
      IENFL=0
      DO IEN=1,NENTRY
        IF(IENTRY(IEN) .EQ. 1 .OR. IENTRY(IEN) .EQ. 2 ) THEN
          IF(JENTRY(IEN) .GT. 0 ) THEN
            IPFL=KENTRY(IEN)
            CALL LCMGTC(IPFL,'SIGNATURE',12,1,HSIGN)
            IF(HSIGN .EQ. 'L_FLUX      ') THEN
              IENFL=IEN
              CALL XDISET(ISTATE,NSTATE,0)
              CALL LCMGET(IPFL,'STATE-VECTOR',ISTATE)
              NGROUP=ISTATE(1)
              NUNKNO=ISTATE(2)
            ENDIF
          ENDIF
        ENDIF
      ENDDO
*----
*  READ PSP OPTIONS
*  IPRINT : EDIT LEVEL
*           = 0 NO EDIT
*           = 1 NORMAL EDIT  (DEFAULT)
*           > 1 EDIT FOR DEBUG
*----
      ALLOCATE(ICOND(NGROUP))
      CALL PSPGET(IPRINT,ITYPE,ICOLR,NGROUP,NGT,ICOND)
*----
*  OPEN POSTSCRIPT OUTPUT FILE
*    1) IF THE FIRST DATA STRUCTURE IS IN UPDATE
*       TEST IF IT IS A POSTSCRIPT FILE CREATED BY DRAGON
*       AND PREPARE FILE FOR OUTPUT
*----
      CALL PSPFIL(ISPSP,JENTRY(1),HENTRY(1),NPAGE)
*----
*  SCAN OVER DATA STRUCTURES AND PROCESS STRUCTURE ONE AFTER THE OTHER
*----
      IF(ITYPE .EQ. 0 .OR. ITYPE .EQ. 1 .OR. ITYPE .EQ. 4) THEN
        NUNKNO=1
        IF(ITYPE .EQ. 0) THEN
          NAMLEG='Region'
        ELSE IF(ITYPE .EQ. 1) THEN
          NAMLEG='Mixture'
        ELSE IF(ITYPE .EQ. 4) THEN
          NAMLEG='HMIX'
        ENDIF
        ALLOCATE(FLUX(NUNKNO))
        FLUX=0.0
        DO 100 IENGT=2,NENTRY
*----
*  READ SIGNATURE OF NEXT DATA STRUCTURE AND TEST IF
*  PSP CAN BE USED TO PROCESS THIS DATA STRUCTURE
*----
          IF(IENTRY(IENGT) .NE. 1   .AND.
     >       IENTRY(IENGT) .NE. 2 ) CALL XABORT(NAMSBR//
     >      ': NEXT DATA STRUCTURE NOT A LINKED LIST OR XSM FILE')
          IF(JENTRY(IENGT) .NE. 2 ) CALL XABORT(NAMSBR//
     >      ': NEXT DATA STRUCTURE NOT IN READ-ONLY MODE')
          IPGT=KENTRY(IENGT)
          NAMGT=HENTRY(IENGT)
          CALL LCMGTC(IPGT,'SIGNATURE',12,1,HSIGN)
          ITRK=1
*----
*  TEST IF GEOMETRY OR EXCELL TRACK DATA STRUCTURE
*----
          IF(HSIGN .EQ. 'L_GEOM      ') THEN
            ITRK=0
            NAMGEO=HENTRY(IENGT)
          ELSE IF(HSIGN .EQ. 'L_TRACK     ') THEN
            CALL LCMGTC(IPGT,'TRACK-TYPE',12,1,HSIGN)
            IF(HSIGN .NE. 'EXCELL') ITRK=-1
          ELSE
            GO TO 115
          ENDIF
*----
*  FOR GEOMETRY OPTION CALL AXGGEO
*  TO GENERATE TEMPORARY TRACKING STRUCTURE
*----
          IF(ITRK .EQ. 0) THEN
            LASS=LDRASS(IPGT,IPRINT)
            CALL  LCMOP(IPTRK2,NAMTR2,IMODT2,IMEDT2,IPRIN2)
            CALL AXGGEO(IPGT  ,IPTRK2,IPRINT,NAMGEO)
            IPGT=IPTRK2
          ENDIF
*----
*  CALL PSPTRK TO GENERATE POSTSCRIPT
*----
          CALL PSPTRK(IPRINT,ISPSP ,ITYPE ,ICOLR ,IPGT  ,NAMGT ,
     >                NAMLEG,NUNKNO,FLUX )
          IF(ITRK .EQ. 0) THEN
            CALL LCMCL(IPTRK2,ICLST2)
          ENDIF
          CALL PSCUTP(ISPSP)
          IF(IENGT .NE. NENTRY) THEN
            NPAGE=NPAGE+1
            CALL PSPAGE(ISPSP,NPAGE,XYPOS)
          ENDIF
 115      CONTINUE
 100    CONTINUE
        DEALLOCATE(FLUX)
      ELSE IF(ITYPE .EQ. 2 .OR. ITYPE .EQ. 3 .OR.
     >        ITYPE .EQ. 5 .OR. ITYPE .EQ. 6 ) THEN
*----
*  TEST SECOND DATA STRUCTURE
*----
        IENGT=2
        IF(IENTRY(IENGT) .NE. 1   .AND.
     >     IENTRY(IENGT) .NE. 2 ) CALL XABORT(NAMSBR//
     >    ': SECOND DATA STRUCTURE NOT A LINKED LIST OR XSM FILE')
        IF(JENTRY(IENGT) .NE. 2 ) CALL XABORT(NAMSBR//
     >    ': SECOND DATA STRUCTURE NOT IN READ-ONLY MODE')
        IPGT=KENTRY(IENGT)
        NAMGT=HENTRY(IENGT)
        CALL LCMGTC(IPGT,'SIGNATURE',12,1,HSIGN)
        ITRK=1
*----
*  TEST IF GEOMETRY OR EXCELL TRACK DATA STRUCTURE
*----
        IF(HSIGN .EQ. 'L_GEOM      ') THEN
          ITRK=0
          NAMGEO=HENTRY(IENGT)
        ELSE IF(HSIGN .EQ. 'L_TRACK     ') THEN
          CALL LCMGTC(IPGT,'TRACK-TYPE',12,1,HSIGN)
          IF(HSIGN .NE. 'EXCELL') ITRK=-1
        ENDIF
*----
*  TEST IF FLUX DATA STRUCTURE EXISTS
*----
        IF(IENFL .EQ. 0) CALL XABORT(NAMSBR//
     >    ': No flux data structure available')
        ALLOCATE(FLUX(NUNKNO),TFLX(NUNKNO))
        IF(ITYPE .EQ. 5 .OR. ITYPE .EQ. 6) THEN
          NMODE=ISTATE(4)
          JPFL=LCMGID(IPFL,'MODE')
        ELSE
          NMODE=1
          JPFL=IPFL
        ENDIF
        DO IMODE=1,NMODE
          IF(ITYPE .EQ. 5 .OR. ITYPE .EQ. 6) THEN
            KPFL=LCMGIL(JPFL,IMODE)
          ELSE
            KPFL=JPFL
          ENDIF
        DO IGR=1,NGT
*----
*  Compute condensed flux
*----
          CALL PSPFCD(KPFL,NGROUP,NUNKNO,IGR,ICOND,FLUX,TFLX)
*----
*  FOR GEOMETRY OPTION CALL AXGGEO
*  TO GENERATE TEMPORARY TRACKING STRUCTURE
*----
          IF(NGT .EQ. 1) THEN
            WRITE(NAMLEG,'(A21)') 'Flux: fully condensed'
          ELSE IF(NGT .EQ. NGROUP) THEN
            WRITE(NAMLEG,'(A18,I5)') 'Flux: tran. group ',IGR
          ELSE
            WRITE(NAMLEG,'(A18,I5)') 'Flux: cond. group ',IGR
          ENDIF
          IF(ITRK .EQ. 0) THEN
            CALL LCMOP(IPTRK2,NAMTR2,IMODT2,IMEDT2,IPRIN2)
            CALL AXGGEO(IPGT  ,IPTRK2,IPRINT,NAMGEO)
            IPGT=IPTRK2
          ENDIF
*----
*  CALL PSPTRK TO GENERATE POSTSCRIPT
*----
          CALL PSPTRK(IPRINT,ISPSP ,ITYPE ,ICOLR ,IPGT  ,NAMGT ,
     >                NAMLEG,NUNKNO,FLUX )
          IF(ITRK .EQ. 0) THEN
            CALL LCMCL(IPTRK2,ICLST2)
          ENDIF
          CALL PSCUTP(ISPSP)
          IF(IGR .NE. NGT) THEN
            NPAGE=NPAGE+1
            CALL PSPAGE(ISPSP,NPAGE,XYPOS)
          ENDIF
        ENDDO
        ENDDO
        DEALLOCATE(TFLX,FLUX)
      ENDIF
      WRITE(ISPSP,'(1X)')
      DEALLOCATE(ICOND)
      RETURN
      END
