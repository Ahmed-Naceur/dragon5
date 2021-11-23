*DECK INFAPL
      SUBROUTINE INFAPL(CFILNA,IPRINT,NBISO,HNAMIS,AWR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To recover mass for isotopes of APOLIB libraries.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): A. Hebert
*
*Parameters: input
* CFILNA  APOLIB1 file name.
* IPRINT  print flag.
* NBISO   number of isotopes.
* HNAMIS  isotope names.
*
*Parameters: output
* AWR     isotope weights.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
* PARAMETERS
*----
      INTEGER    IOUT,MAXIT
      PARAMETER (IOUT=6,MAXIT=1000)
*----
* FUNCTIONS
*----
      INTEGER    KDROPN,KDRCLS
      DOUBLE PRECISION XDRCST
*----
* LOCAL VARIABLES
*----
      INTEGER    NBISO,IPRINT,IUNIT,IISO,INDLOR,NR,NIT,I,K,IMX,NISB,
     1           NRST,ICC,NS1,IC,NRSTR,NN,IER
      INTEGER    IT(MAXIT)
*
      REAL       AA
      REAL       AWR(NBISO)
      CHARACTER  CFILNA*64,HNAMIS(NBISO)*8,HNISOR*8,FORM*4
      EQUIVALENCE(AA,NN)
      REAL       CONVM
*----
* OPEN APOLIB
*----
      CONVM=REAL(XDRCST('Neutron mass','amu'))
      IF( IPRINT.GT.0 ) THEN
        WRITE(IOUT,6000) CFILNA
      ENDIF
      IUNIT=KDROPN(CFILNA,2,2,0)
      IF( IUNIT.LE.0 )THEN
        WRITE(IOUT,9000) CFILNA
        CALL XABORT('INFAPL: APOL LIBRARY CANNOT BE OPENED')
      ENDIF
      IISO= 0
      REWIND(IUNIT)
   50 READ(IUNIT) INDLOR,NR,NIT,(IT(I),I=1,NIT)
      IF( NIT.GT.MAXIT ) CALL XABORT('INFAPL: MAXIT IS TOO SMALL')
      IF(INDLOR.EQ.9999) GO TO 700
      DO 70 IMX=1,NBISO
         HNISOR= HNAMIS(IMX)
         I=INDEX(HNISOR,' ')
         IF(I.EQ.0) THEN
            READ(HNISOR,'(I8)') NISB
         ELSE
            WRITE(FORM,'(2H(I,I1,1H))') I-1
            READ(HNISOR,FORM) NISB
         ENDIF
         IF( NISB.EQ.INDLOR )THEN
            IF( IPRINT.GT.0 ) WRITE(IOUT,6001) HNISOR
            IISO= IISO + 1
            NRST= IT(4)
            NS1= 0
            IF( IT(5).LT.0 ) NS1= -IT(5)
            IC=5+NS1+NRST
            NRSTR=IT(IC)
            ICC=IC+6*NRSTR+1
            NN=IT(ICC)
            AWR(IMX)=AA*CONVM
         ENDIF
   70 CONTINUE
      DO 80 K=1,NR
      READ(IUNIT)
   80 CONTINUE
      GO TO 50
*
*     CHECK IF ALL NBISO ISOTOPES HAVE BEEN PROCESSED.
  700 IF( IISO.NE.NBISO )THEN
         CALL XABORT('INFAPL: SOME ISOTOPES WERE NOT RECOVERED')
      ENDIF
*
*     CLOSE APOLIB FILE.
      IER=KDRCLS(IUNIT,1)
      IF(IER.LT.0) CALL XABORT(
     > 'INFAPL: Impossible to close library '//CFILNA)
      RETURN
*
 9000 FORMAT(/' ERROR IN PROCESSING APOL LIBRARY:',A8)
 6000 FORMAT(/' PROCESSING APOL LIBRARY NAME ',A8)
 6001 FORMAT(/'    PROCESSING ISOTOPE/MATERIAL = ',A12)
      END
