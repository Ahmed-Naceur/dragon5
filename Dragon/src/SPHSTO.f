*DECK SPHSTO
      SUBROUTINE SPHSTO(IPSAP,ICAL,IMPX,LNEW,HEQUI,HEQNAM,NMIL,NGROUP,
     1 SPH)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Store a new set of SPH factors for an elementary calculation in a
* Saphyb.
*
*Copyright:
* Copyright (C) 2011 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPSAP   pointer to the Saphyb (L_SAPHYB signature).
* ICAL    index of the elementary calculation being considered.
* IMPX    print parameter (equal to zero for no print).
* LNEW    flag set to .TRUE. to allow the overwriting of the existing
*         set of SPH factors named HEQUI.
* HEQUI   LOCKEY name of SPH-factor set to be stored.
* HEQNAM  LOCNAM name of SPH-factor set to be stored.
* NMIL    number of mixtures in the elementary calculation.
* NGROUP  number of energy groups in the elementary calculation.
* SPH     SPH-factor set to be stored the Saphyb.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPSAP
      INTEGER ICAL,IMPX,NMIL,NGROUP
      REAL SPH(NMIL,NGROUP)
      LOGICAL LNEW
      CHARACTER HEQUI*4,HEQNAM*80
*----
*  LOCAL VARIABLES
*----
      PARAMETER(MAXLOC=10)
      INTEGER DIMSAP(50)
      CHARACTER TEXT12*12,HSMG*131,LOCTYP(MAXLOC)*4,LOCNAM(MAXLOC)*80,
     1 LOCKEY(MAXLOC)*4
      INTEGER, ALLOCATABLE, DIMENSION(:) :: LOCAD
      REAL, ALLOCATABLE, DIMENSION(:) :: RVALO
*----
*  RECOVER SAPHYB CHARACTERISTICS
*----
      IF(HEQUI.EQ.' ') CALL XABORT('SPHSTO: HEQUI NOT DEFINED')
      CALL LCMLEN(IPSAP,'DIMSAP',ILENG,ITYLCM)
      IF(ILENG.EQ.0) CALL XABORT('SPHSTO: DIMSAP NOT DEFINED')
      CALL LCMGET(IPSAP,'DIMSAP',DIMSAP)
      NMIL=DIMSAP(7)   ! number of mixtures
      NCALS=DIMSAP(19) ! number of elementary calculations in the SAPHYB
      NGROUP=DIMSAP(20)! number of energy groups
      IF(IMPX.GT.0) THEN
        WRITE(6,'(29H SPHSTO: number of mixtures =,I5)') NMIL
        WRITE(6,'(33H SPHSTO: number of calculations =,I5)') NCALS
        WRITE(6,'(34H SPHSTO: number of energy groups =,I4)') NGROUP
      ENDIF
      IF(ICAL.GT.NCALS) CALL XABORT('SPHSTO: ICAL INDEX OVERFLOW')
*----
*  RECOVER INFORMATION FROM caldir DIRECTORY.
*----
      WRITE(TEXT12,'(4Hcalc,I8)') ICAL
      CALL LCMLEN(IPSAP,TEXT12,ILENG,ITYLCM)
      IF(ILENG.EQ.0) THEN
         WRITE(HSMG,'(29HSPHSTO: MISSING CALCULATION '',A12,2H''.)')
     1   TEXT12
         CALL XABORT(HSMG)
      ENDIF
      CALL LCMSIX(IPSAP,TEXT12,1)
      CALL LCMSIX(IPSAP,'info',1)
      CALL LCMGET(IPSAP,'NLOC',NLOC)
      IF(NLOC+1.GT.MAXLOC) CALL XABORT('SPHSTO: MAXLOC OVERFLOW')
      CALL LCMGTC(IPSAP,'LOCTYP',4,NLOC,LOCTYP)
      CALL LCMGTC(IPSAP,'LOCNAM',80,NLOC,LOCNAM)
      CALL LCMGTC(IPSAP,'LOCKEY',4,NLOC,LOCKEY)
      ALLOCATE(LOCAD(NLOC+2))
      CALL LCMGET(IPSAP,'LOCADR',LOCAD)
      DO ILOC=1,NLOC
        IF ((LOCTYP(ILOC).EQ.'EQUI').AND.(LOCKEY(ILOC).EQ.HEQUI)) THEN
*         SET HEQUI EXISTS.
          IF(LNEW) THEN
            IF(IMPX.GT.0) WRITE(6,'(31H SPHSTO: OVERWRITE SPH-FACTOR S,
     1      9HET NAMED ,A)') HEQUI
            JLOC=ILOC
            GO TO 10
          ELSE
            CALL XABORT('SPHSTO: THIS SPH FACTOR SET EXISTS: '//HEQUI)
          ENDIF
        ENDIF
      ENDDO
*     A NEW SET OF SPH FACTORS IS DEFINED IN THE SAPHYB
      JLOC=NLOC+1
      NLOC=NLOC+1
      LOCTYP(NLOC)='EQUI'
      LOCKEY(NLOC)=HEQUI
      IF(HEQNAM.NE.' ') THEN
        LOCNAM(NLOC)=HEQNAM
      ELSE
        LOCNAM(NLOC)=HEQUI
      ENDIF
      LOCAD(NLOC+1)=LOCAD(NLOC)+NGROUP
      CALL LCMPUT(IPSAP,'NLOC',1,1,NLOC)
      CALL LCMPTC(IPSAP,'LOCTYP',4,NLOC,LOCTYP)
      CALL LCMPTC(IPSAP,'LOCNAM',80,NLOC,LOCNAM)
      CALL LCMPTC(IPSAP,'LOCKEY',4,NLOC,LOCKEY)
      CALL LCMPUT(IPSAP,'LOCADR',NLOC+1,1,LOCAD)
   10 CALL LCMSIX(IPSAP,' ',2)
*----
*  LOOP OVER MIXTURES.
*----
      DO IBM=1,NMIL
        WRITE(TEXT12,'(4Hmili,I8)') IBM
        CALL LCMLEN(IPSAP,TEXT12,ILENG,ITYLCM)
        IF(ILENG.EQ.0) THEN
          WRITE(HSMG,'(29HSPHSTO: MISSING MIXTURE '',A12,2H''.)')
     1    TEXT12
          CALL XABORT(HSMG)
        ENDIF
        CALL LCMSIX(IPSAP,TEXT12,1)
        ALLOCATE(RVALO(LOCAD(JLOC+1)))
        CALL LCMGET(IPSAP,'RVALOC',RVALO)
        DO IGR=1,NGROUP
          RVALO(LOCAD(JLOC)+IGR-1)=SPH(IBM,IGR)
        ENDDO
        CALL LCMPUT(IPSAP,'RVALOC',LOCAD(JLOC+1)-1,2,RVALO)
        DEALLOCATE(RVALO)
        CALL LCMSIX(IPSAP,' ',2)
      ENDDO
      DEALLOCATE(LOCAD)
      CALL LCMSIX(IPSAP,' ',2)
      RETURN
      END
