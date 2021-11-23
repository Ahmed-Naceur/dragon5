*DECK LIBDI3
      SUBROUTINE LIBDI3 (MAXDIL,NAMFIL,HNISOR,NDIL,DILUT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Find the dilutions corresponding to a resonant isotope within a
* library in matxs (njoy-91) format.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* MAXDIL  maximum number of dilutions.
* NAMFIL  name of the MATXS file.
* HNISOR  library name of the resonant isotope.
*
*Parameters: output
* NDIL    number of finite dilutions.
* DILUT   dilutions.
*
*-----------------------------------------------------------------------
*
      USE XDRMOD
      USE LIBEEDR
      IMPLICIT CHARACTER*6 (H)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER MAXDIL,NDIL
      CHARACTER NAMFIL*(*),HNISOR*12
      REAL DILUT(MAXDIL)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (MULT=2,MAXA=10000)
      REAL A(MAXA)
      INTEGER IA(MAXA)
      CHARACTER HSMG*131
      DOUBLE PRECISION XHA(MAXA/2)
      EQUIVALENCE (A(1),IA(1),XHA(1))
*
      ILIBIN=2
      IF(NAMFIL(:1).EQ.'_') ILIBIN=3
      NIN=KDROPN(NAMFIL,2,ILIBIN,0)
      IF(NIN.LE.0) THEN
         WRITE (HSMG,'(35HLIBDI3: UNABLE TO OPEN LIBRARY FILE,1X,A16,
     1   6H. NIN=,I4,1H.)') NAMFIL,NIN
         CALL XABORT(HSMG)
      ENDIF
      NDIL=0
      NWDS=1+3*MULT
      IREC=1
*     --FILE IDENTIFICATION--------------
      IF(ILIBIN.EQ.2) THEN
         CALL XDREED (NIN,IREC,A(1),NWDS)
      ELSE IF(ILIBIN.EQ.3) THEN
         CALL LIBEED (NIN,IREC,A(1),NWDS)
      ENDIF
*     -----------------------------------
*
      NWDS=6
      IREC=2
*     --FILE CONTROL---------------------
      IF(ILIBIN.EQ.2) THEN
         CALL XDREED (NIN,IREC,A(1),NWDS)
      ELSE IF(ILIBIN.EQ.3) THEN
         CALL LIBEED (NIN,IREC,A(1),NWDS)
      ENDIF
*     -----------------------------------
      NPART=IA(1)
      NTYPE=IA(2)
      NHOLL=IA(3)
      NMAT=IA(4)
*
      NWDS=NHOLL*MULT
      IF(NWDS.GT.MAXA)
     1 CALL XABORT('LIBDI3: INSUFFICIENT VALUE OF MAXA(1).')
      IREC=3
*     --HOLLERITH IDENTIFICATION---------
      IF(ILIBIN.EQ.2) THEN
         CALL XDREED (NIN,IREC,A(1),NWDS)
      ELSE IF(ILIBIN.EQ.3) THEN
         CALL LIBEED (NIN,IREC,A(1),NWDS)
      ENDIF
*     -----------------------------------
      NWDS=(NPART+NTYPE+NMAT)*MULT+2*NTYPE+NPART+2*NMAT
      IF(NWDS.GT.MAXA)
     1 CALL XABORT('LIBDI3: INSUFFICIENT VALUE OF MAXA(2).')
      IREC=4
*     --FILE DATA------------------------
      IF(ILIBIN.EQ.2) THEN
         CALL XDREED (NIN,IREC,A(1),NWDS)
      ELSE IF(ILIBIN.EQ.3) THEN
         CALL LIBEED (NIN,IREC,A(1),NWDS)
      ENDIF
*     -----------------------------------
      IF((NWDS/2)*2.NE.NWDS) NWDS=NWDS+1
      L2=1+NWDS
      IREC=IREC+NPART
      IRZM=IREC+1
*
      DO 50 IM=1,NMAT
      WRITE (HMAT,'(A6)') XHA(NPART+NTYPE+IM)
      IF(HMAT.NE.HNISOR(:6)) GO TO 50
*
      LOC=(NPART+NTYPE+NMAT)*MULT+NPART+2*NTYPE+IM
      NSUB=IA(LOC)
      LOCM=IA(LOC+NMAT)
      IREC=LOCM+IRZM
      NWDS=MULT+1+6*NSUB
      IF(L2+NWDS-1.GT.MAXA)
     1 CALL XABORT('LIBDI3: INSUFFICIENT VALUE OF MAXA(3).')
*     --MATERIAL CONTROL------------------
      IF(ILIBIN.EQ.2) THEN
         CALL XDREED (NIN,IREC,A(L2),NWDS)
      ELSE IF(ILIBIN.EQ.3) THEN
         CALL LIBEED (NIN,IREC,A(L2),NWDS)
      ENDIF
*     ------------------------------------
      IF((NWDS/2)*2.NE.NWDS) NWDS=NWDS+1
*
      DO 40 ISUBM=1,NSUB
      DILI=A(L2+MULT+6*(ISUBM-1)+2)
      DO 10 I=1,NDIL
      IF(ABS(DILI-DILUT(I)).LT.1.0E-5*ABS(DILI)) GO TO 40
   10 CONTINUE
      DO 30 I=1,NDIL
      IF(DILI.LT.DILUT(I)) THEN
         DO 20 J=NDIL,I,-1
         DILUT(J+1)=DILUT(J)
   20    CONTINUE
         DILUT(I)=DILI
         NDIL=NDIL+1
         IF(NDIL.GT.MAXDIL) CALL XABORT('LIBDI3: MAXDIL IS TOO SMALL.')
         GO TO 40
      ENDIF
   30 CONTINUE
      NDIL=NDIL+1
      IF(NDIL.GT.MAXDIL) CALL XABORT('LIBDI3: MAXDIL IS TOO SMALL.')
      DILUT(NDIL)=DILI
   40 CONTINUE
   50 CONTINUE
      NDIL=NDIL-1
      IF(NDIL.LT.0) CALL XABORT('LIBDI3: UNABLE TO FIND THE TABULATED'
     1 //' DILUTIONS.')
*     --CLOSE CCCC FILE--
      IF(ILIBIN.EQ.2) THEN
         CALL XDRCLS(NIN)
      ELSE IF(ILIBIN.EQ.3) THEN
         CALL LIBCLS()
      ENDIF
*     -------------------
      IER=KDRCLS(NIN,1)
      IF(IER.LT.0) THEN
         WRITE (HSMG,'(36HLIBDI3: UNABLE TO CLOSE LIBRARY FILE,1X,A16,
     1   1H.)') NAMFIL
         CALL XABORT(HSMG)
      ENDIF
      RETURN
      END
