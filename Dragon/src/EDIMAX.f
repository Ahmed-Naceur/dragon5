*DECK EDIMAX
      SUBROUTINE EDIMAX(NBISO,ISONAM,MIX,IPRINT,NREGIO,NMERGE,MATCOD,
     1 IMERGE,LSISO,LISO,MAXISO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the number of homogenized/condensed isotopes in the output
* microlib.
*
*Copyright:
* Copyright (C) 2017 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NBISO   number of isotopes in the input microlib.
* ISONAM  local names of NBISO isotopes:
*         chars 1 to 8  is the local isotope name;
*         chars 9 to 12 is a suffix function of the mix number.
* MIX     mixture number associated with each isotope.
* IPRINT  print index.
* NREGIO  number of volumes.
* NMERGE  number of merged regions.
* MATCOD  mixture index per volume.
* IMERGE  index of merged regions.
* LSISO   flag for isotopes saved.
* LISO    =.TRUE. if we want to keep all the isotopes after 
*         homogeneization.
*
*Parameters: output
* MAXISO  number of homogenized/condensed isotopes in the output
*         microlib.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER   NBISO,ISONAM(3,NBISO),MIX(NBISO),IPRINT,NREGIO,NMERGE,
     1          MATCOD(NREGIO),IMERGE(NREGIO),LSISO(NBISO),MAXISO
      LOGICAL   LISO
*----
*  LOCAL VARIABLES
*----
      LOGICAL    LOGIC
*----
*  ALLOCATABLE ARRAYS
*----
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: MASK
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(MASK(NBISO))
*
      JJISO=0
      IF(IPRINT.GT.5) THEN
        WRITE(6,'(/34H EDIMAX: MERGING FLAG PER ISOTOPE:/(1X,20I6))')
     1  (LSISO(I),I=1,NBISO)
      ENDIF
      DO 70 INM=1,NMERGE
      CALL XDLSET(MASK,NBISO,.FALSE.)
      DO 60 ISO=1,NBISO
      IF(MASK(ISO).OR.(LSISO(ISO).EQ.0)) GO TO 60
      DO 10 IREGIO=1,NREGIO
      IF((IMERGE(IREGIO).EQ.INM).AND.(MATCOD(IREGIO).EQ.MIX(ISO)))
     1 GO TO 20
   10 CONTINUE
      GO TO 60
   20 LOGIC=.FALSE.
      DO 50 IREGIO=1,NREGIO
      MATNUM=MATCOD(IREGIO)
      IF(IMERGE(IREGIO).EQ.INM) THEN
         DO 40 JSO=ISO,NBISO
         IF((ISONAM(1,ISO).EQ.ISONAM(1,JSO)).AND.
     1      (ISONAM(2,ISO).EQ.ISONAM(2,JSO)).AND.
     2      (MATNUM.EQ.MIX(JSO)).AND.(LSISO(JSO).NE.0)) THEN
            IF(LISO) THEN
              IF(ISONAM(3,ISO).EQ.ISONAM(3,JSO)) GOTO 30
              GOTO 40
            ENDIF
   30       LOGIC=.TRUE.
            MASK(JSO)=.TRUE.
            GO TO 40
         ENDIF
   40    CONTINUE
      ENDIF
   50 CONTINUE
      IF(LOGIC) THEN
        JJISO=JJISO+1
        IF(JJISO.GT.NBISO*NMERGE) CALL XABORT('EDIMAX: JJISO OVERFLOW.')
      ENDIF
   60 CONTINUE
   70 CONTINUE
      MAXISO=JJISO
      IF(IPRINT.GT.1) WRITE(6,100) MAXISO
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(MASK)
      RETURN
  100 FORMAT(/53H EDIMAX: NUMBER OF HOMOGENIZED/CONDENSED ISOTOPES IN ,
     1 20HTHE OUTPUT MICROLIB=,I8)
      END
