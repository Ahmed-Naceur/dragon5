*DECK INFTR2
      SUBROUTINE INFTR2(CFILNA,IPRINT,NBISO,HNAMIS,AWRISO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To recover mass for isotopes of MATXS type libraries
* use MATXS format from NJOY-91.
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
* CFILNA  file name.
* IPRINT  print flag.
* NBISO   number of isotopes.
* HNAMIS  isotope names.
*
*Parameters: output
* AWRISO  isotope weights
*
*Reference:
* R. E. MACFARLANE, TRANSX-CTR: A code for interfacing
* MATXS cross-section libraries to nuclear transport codes for
* fusion systems analysis, Los Alamos National Laboratory,
* Report LA-9863-MS, New Mexico, February 1984.
*
*-----------------------------------------------------------------------
*
      USE XDRMOD
      IMPLICIT         NONE
      INTEGER          IPRINT,NBISO
      CHARACTER        CFILNA*8,HNAMIS(NBISO)*64
      REAL             AWRISO(NBISO)
C----
C  LOCAL VARIABLES
C----
      INTEGER          IOUT,MULT,MAXA
      CHARACTER        FORM*4
      PARAMETER       (IOUT=6,MULT=2,MAXA=1000,FORM='(A6)')
C----
C FUNCTIONS
C----
      INTEGER          KDROPN,KDRCLS
      DOUBLE PRECISION XDRCST
      INTEGER          NIN,IREC,NWDS,NPART,NTYPE,NMAT,L2,L2H,IRZM,IM,
     >                 ISO,LOC,IER,IA(MAXA)
      CHARACTER        HSMG*131,HMAT*6
      REAL             RA(MAXA)
      DOUBLE PRECISION DA(MAXA/2)
      REAL             CONVM
      EQUIVALENCE     (RA(1),IA(1),DA(1))
C----
C  OPEN MATXS FILE AND INITIALIZE LIBRARY
C----
      CONVM=REAL(XDRCST('Neutron mass','amu'))
      NIN=KDROPN(CFILNA,2,2,0)
      IF(NIN.LE.0) THEN
        WRITE(HSMG,9000) CFILNA
        CALL XABORT(HSMG)
      ENDIF
      IREC=2
      NWDS=6
C-------FILE CONTROL---------------
      CALL XDREED(NIN,IREC,RA,NWDS)
C----------------------------------
      NPART=IA(1)
      NTYPE=IA(2)
      NMAT=IA(4)
      IREC=4
      NWDS=(NPART+NTYPE+NMAT)*MULT+2*NTYPE+NPART+2*NMAT
      IF(NWDS.GT.MAXA) CALL XABORT
     >  ('INFTR2: LENGTH OF RECORD 4 > MAXA ')
C-------FILE DATA------------------
      CALL XDREED(NIN,IREC,RA,NWDS)
C----------------------------------
      IF((NWDS/2)*2.NE.NWDS) NWDS=NWDS+1
      L2=1+NWDS
      L2H=NWDS/MULT+1
      IRZM=5+NPART
C----
C  READ THROUGH MATXS FILE AND GET AWR FOR ISOTOPES
C----
      DO 100 IM=1,NMAT
        WRITE(HMAT,FORM) DA(L2H-1+IM)
        DO 110 ISO=1,NBISO
          IF(HMAT.EQ.HNAMIS(ISO)(:6)) THEN
            LOC=(NPART+NTYPE+NMAT)*MULT+NPART+2*NTYPE+IM
            IREC=IA(LOC+NMAT)+IRZM
            NWDS=MULT+1+6*IA(LOC)
            IF(L2+NWDS-1.GT.MAXA) CALL XABORT
     >        ('INFTR2: LENGTH OF CURRENT RECORD > MAXA ')
C-------------------------------------------
              CALL XDREED(NIN,IREC,RA(L2),NWDS)
C-------------------------------------------
            AWRISO(ISO)=RA(L2+MULT)*CONVM
            IF(IPRINT.GE.100) THEN
              WRITE(IOUT,6000) HNAMIS(ISO),AWRISO(ISO)
            ENDIF
          ENDIF
 110    CONTINUE
 100  CONTINUE
C----
C  CLOSE MATXS FILE.
C----
      CALL XDRCLS(NIN)
      IER=KDRCLS(NIN,1)
      IF(IER.LT.0) THEN
        WRITE(HSMG,9001) CFILNA
        CALL XABORT(HSMG)
      ENDIF
      RETURN
C----
C  PRINT FORMATS
C----
 6000 FORMAT(' MATXS ISOTOPE =',A8,
     >       ' HAS ATOMIC WEIGHT RATIO = ',F10.3)
C----
C  ABORT FORMATS
C----
 9000 FORMAT('INFTR2: UNABLE TO OPEN MATXS LIBRARY FILE ',A64)
 9001 FORMAT('INFTR2: UNABLE TO CLOSE MATXS LIBRARY FILE ',A64)
      END
