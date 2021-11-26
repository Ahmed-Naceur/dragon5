*DECK INFTR1
      SUBROUTINE INFTR1(CFILNA,IPRINT,NBISO,HNAMIS,AWRISO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To recover mass for isotopes of MATXS type libraries
* use MATXS format from NJOY-II or NJOY89.
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
* AWRISO  isotope weights.
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
      CHARACTER        CFILNA*64,HNAMIS(NBISO)*8
      REAL             AWRISO(NBISO)
*----
*  LOCAL VARIABLES
*----
      INTEGER          IOUT,MULT,MAXA
      CHARACTER        FORM*4
      PARAMETER       (IOUT=6,MULT=2,MAXA=1000,FORM='(A6)')
*----
* FUNCTIONS
*----
      INTEGER          KDROPN,KDRCLS
      DOUBLE PRECISION XDRCST
      INTEGER          NIN,IREC,NWDS,NPART,NTYPE,L2,L2H,IRZT,IT,
     >                 NDEX,NMAT,NINP,NING,NOUTP,NOUTG,LOCT,LMC,
     >                 IRZM,IM,ISO,LOC,IER,IA(MAXA)
      CHARACTER        HSMG*131,HTYPE*6,HMAT*6
      REAL             RA(MAXA)
      DOUBLE PRECISION DA(MAXA/2)
      REAL             CONVM
      EQUIVALENCE     (RA(1),IA(1),DA(1))
*----
*  OPEN MATXS FILE AND INITIALIZE LIBRARY
*----
      CONVM=REAL(XDRCST('Neutron mass','amu'))
      NIN=KDROPN(CFILNA,2,2,0)
      IF(NIN.LE.0) THEN
        WRITE(HSMG,9000) CFILNA
        CALL XABORT(HSMG)
      ENDIF
      IREC=2
      NWDS=3
*-------FILE CONTROL---------------
      CALL XDREED(NIN,IREC,RA,NWDS)
*----------------------------------
      NPART=IA(1)
      NTYPE=IA(2)
      IREC=4
      NWDS=(NPART+NTYPE)*MULT+6*NTYPE+NPART
      IF(NWDS.GT.MAXA) CALL XABORT
     >  ('INFTR1: LENGTH OF RECORD 4 > MAXA ')
*-------FILE DATA------------------
      CALL XDREED(NIN,IREC,RA,NWDS)
*----------------------------------
      IF((NWDS/2)*2.NE.NWDS) NWDS=NWDS+1
      L2=1+NWDS
      L2H=(L2-1)/MULT+1
      IRZT=5+NPART
*----
*  DATA TYPE LOOP
*----
      DO 100 IT=1,NTYPE
        WRITE(HTYPE,FORM) DA(NPART+IT)
        CALL XDRCAS('LOWTOUP',HTYPE)
        IF(HTYPE.NE.'NSCAT'.AND.HTYPE.NE.'NTHERM') GO TO 105
        NDEX=(NPART+NTYPE)*MULT+IT
        NMAT=IA(NDEX)
        NDEX=NDEX+NTYPE
        NINP=IA(NDEX)
        NDEX=NDEX+NTYPE
        NING=IA(NDEX)
        NDEX=NDEX+NTYPE
        NOUTP=IA(NDEX)
        NDEX=NDEX+NTYPE
        NOUTG=IA(NDEX)
        NDEX=NDEX+NTYPE
        LOCT=IA(NDEX)
*----
*  DATA TYPE CONTROL
*----
        IREC=LOCT+IRZT
        NWDS=(2+MULT)*NMAT+NINP+NOUTP+1
        IF(L2+NWDS-1.GT.MAXA)  CALL XABORT
     >    ('INFTR1: LENGTH OF CURRENT RECORD > MAXA ')
*----------------------------------------
        CALL XDREED(NIN,IREC,RA(L2),NWDS)
*----------------------------------------
        IF((NWDS/2)*2.NE.NWDS) NWDS=NWDS+1
        LMC=L2+NWDS
        IRZM=IREC+1
*----
*  READ THROUGH MATXS FILE AND GET AWR FOR ISOTOPES
*----
        DO 110 IM=1,NMAT
          WRITE(HMAT,FORM) DA(L2H-1+IM)
          DO 120 ISO=1,NBISO
            IF(HMAT.EQ.HNAMIS(ISO)(:6)) THEN
              LOC=L2-1+MULT*NMAT+IM
              IREC=IA(LOC+NMAT)+IRZM
              NWDS=MULT+1+6*IA(LOC)
              IF(LMC+NWDS-1.GT.MAXA) CALL XABORT
     >          ('INFTR1: LENGTH OF CURRENT RECORD > MAXA ')
*-------------------------------------------
              CALL XDREED(NIN,IREC,RA(LMC),NWDS)
*-------------------------------------------
              AWRISO(ISO)=RA(LMC+MULT)*CONVM
              IF(IPRINT.GE.100) THEN
                WRITE(IOUT,6000) HNAMIS(ISO),AWRISO(ISO)
              ENDIF
            ENDIF
 120      CONTINUE
 110    CONTINUE
 105    CONTINUE
 100  CONTINUE
*----
*  CLOSE MATXS FILE.
*----
      CALL XDRCLS(NIN)
      IER=KDRCLS(NIN,1)
      IF(IER.LT.0) THEN
        WRITE(HSMG,9001) CFILNA
        CALL XABORT(HSMG)
      ENDIF
      RETURN
*----
*  PRINT FORMATS
*----
 6000 FORMAT(' MATXS ISOTOPE =',A8,
     >       ' HAS ATOMIC WEIGHT RATIO = ',F10.3)
*----
*  ABORT FORMATS
*----
 9000 FORMAT('INFTR1: UNABLE TO OPEN MATXS LIBRARY FILE ',A64)
 9001 FORMAT('INFTR1: UNABLE TO CLOSE MATXS LIBRARY FILE ',A64)
      END
