*DECK SAPFWC
      SUBROUTINE SAPFWC(IPSAP,IPLB1,MAXISO,NBISO,NGA,NISOTA)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To fill the 'constphysiq' directory.
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
* IPSAP   pointer to the Saphyb.
* IPLB1   pointer to the reference microlib.
* MAXISO  allocated length of arrays ISOTA and ISOTYP.
* NBISO   number of isotopic definitions in the reference microlib.
* NGA     number of energy groups in the reference microlib.
*
*Parameters: output
* NISOTA  number of distinct isotopes.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPSAP,IPLB1
      INTEGER MAXISO,NBISO,NGA,NISOTA
*----
*  LOCAL VARIABLES
*----
      CHARACTER TEXT12*12
*----
*  ALLOCATABLE ARRAYS
*----
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: ISOTA
      CHARACTER(LEN=4), ALLOCATABLE, DIMENSION(:) :: ISOTYP
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ITYP
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISONAM
      REAL, ALLOCATABLE, DIMENSION(:) :: ENRGA
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(ISONAM(3,NBISO),ITYP(NBISO))
      ALLOCATE(ENRGA(NGA+1))
      ALLOCATE(ISOTA(MAXISO),ISOTYP(MAXISO))
*
      CALL LCMSIX(IPSAP,'constphysiq',1)
*
      CALL LCMGET(IPLB1,'ENERGY',ENRGA)
      DO 10 I=1,NGA+1
      ENRGA(I)=ENRGA(I)*1.0E-6
   10 CONTINUE
      CALL LCMPUT(IPSAP,'ENRGA',NGA+1,2,ENRGA)
*
      CALL LCMGET(IPLB1,'ISOTOPESUSED',ISONAM)
      CALL LCMGET(IPLB1,'ISOTOPESTYPE',ITYP)
      NISOTA=0
      DO 30 IBISO=1,NBISO
      WRITE(TEXT12,'(3A4)') (ISONAM(I0,IBISO),I0=1,3)
      DO 20 ISO=1,NISOTA
      IF(TEXT12(:8).EQ.ISOTA(ISO)) GO TO 30
   20 CONTINUE
      NISOTA=NISOTA+1
      IF(NISOTA.GT.MAXISO) CALL XABORT('SAPFWC: ISOTA OVERFLOW.')
      ISOTA(NISOTA)=TEXT12(:8)
      IF(ITYP(IBISO).EQ.2) THEN
         ISOTYP(NISOTA)='FISS'
      ELSE IF(ITYP(IBISO).EQ.3) THEN
         ISOTYP(NISOTA)='F.P.'
      ELSE
         ISOTYP(NISOTA)=' '
      ENDIF
   30 CONTINUE
      CALL LCMPTC(IPSAP,'ISOTA',8,NISOTA,ISOTA)
      CALL LCMPTC(IPSAP,'ISOTYP',4,NISOTA,ISOTYP)
*
      CALL LCMSIX(IPSAP,' ',2)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(ISOTYP,ISOTA)
      DEALLOCATE(ENRGA)
      DEALLOCATE(ITYP,ISONAM)
      RETURN
      END
