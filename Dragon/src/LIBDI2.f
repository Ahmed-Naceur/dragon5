*DECK LIBDI2
      SUBROUTINE LIBDI2 (MAXDIL,NAMFIL,HNISOR,NDIL,DILUT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Find the dilutions corresponding to a resonant isotope within a
* library in matxs (njoy-89) format.
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
      PARAMETER (MULT=2,MAXA=1000)
      REAL A(MAXA)
      INTEGER IA(MAXA)
      CHARACTER HSMG*131
      DOUBLE PRECISION HA(MAXA/2)
      EQUIVALENCE (A(1),IA(1),HA(1))
*
      NIN=KDROPN(NAMFIL,2,2,0)
      IF(NIN.LE.0) THEN
         WRITE (HSMG,'(35HLIBDI2: UNABLE TO OPEN LIBRARY FILE,1X,A16,
     1   6H. NIN=,I4,1H.)') NAMFIL,NIN
         CALL XABORT(HSMG)
      ENDIF
      NDIL=0
      NWDS=1+3*MULT
      IREC=1
*     --------------------------------
      CALL XDREED (NIN,IREC,A(1),NWDS)
*     --------------------------------
*
      NWDS=3
      IREC=2
*     --------------------------------
      CALL XDREED (NIN,IREC,A(1),NWDS)
*     --------------------------------
      NPART=IA(1)
      NTYPE=IA(2)
      NHOLL=IA(3)
      NWDS=NHOLL*MULT
      IF(NWDS.GT.MAXA)
     1 CALL XABORT('LIBDI2: INSUFFICIENT VALUE OF MAXA(1).')
      IREC=3
*     --------------------------------
      CALL XDREED (NIN,IREC,A(1),NWDS)
*     --------------------------------
      NWDS=(NPART+NTYPE)*MULT+6*NTYPE+NPART
      IF(NWDS.GT.MAXA)
     1 CALL XABORT('LIBDI2: INSUFFICIENT VALUE OF MAXA(2).')
      IREC=4
*     --------------------------------
      CALL XDREED (NIN,IREC,A(1),NWDS)
*     --------------------------------
      NWC=NPART+NTYPE
      IF((NWDS/2)*2.NE.NWDS) NWDS=NWDS+1
      L2=1+NWDS
      L2H=(L2-1)/MULT+1
      NEX1=(NPART+NTYPE)*MULT+6*NTYPE
      IREC=IREC+NPART
      IRZT=5+NPART
      DO 680 IT=1,NTYPE
      WRITE(HTYPE,'(A6)') HA(NPART+IT)
      IF(HTYPE.NE.'NSCAT') GO TO 680
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
      NWDS=(2+MULT)*NMAT+NINP+NOUTP+1
      IF(L2+NWDS-1.GT.MAXA)
     1 CALL XABORT('LIBDI2: INSUFFICIENT VALUE OF MAXA(3).')
      IREC=LOCT+IRZT
*     ---------------------------------
      CALL XDREED (NIN,IREC,A(L2),NWDS)
*     ---------------------------------
      IF((NWDS/2)*2.NE.NWDS) NWDS=NWDS+1
      LMC=L2+NWDS
      LMCH=L2H+NWDS/MULT
      NSBLK=IA(L2+NMAT*(MULT+2)+NINP+NOUTP)
      IRZM=IREC+1
*----
*  MATERIAL/ISOTOPE LOOP
*----
      DO 670 IM=1,NMAT
      WRITE (HMAT,'(A6)') HA(L2H-1+IM)
      IF(HMAT.NE.HNISOR(:6)) GO TO 670
*
      LOC=L2-1+MULT*NMAT+IM
      NSUBM=IA(LOC)
      LOCM=IA(LOC+NMAT)
      IREC=LOCM+IRZM
      NWDS=MULT+1+6*NSUBM
      IF(LMC+NWDS-1.GT.MAXA)
     1 CALL XABORT('LIBDI2: INSUFFICIENT VALUE OF MAXA(4).')
*     ----------------------------------
      CALL XDREED (NIN,IREC,A(LMC),NWDS)
*     ----------------------------------
      NWDS=NWDS+MULT-1
      DO 307 ISUBM=1,NSUBM
      DILI=A(LMC+MULT+6*(ISUBM-1)+2)
      DO 555 I=1,NDIL
      IF(ABS(DILI-DILUT(I)).LT.1.0E-5*ABS(DILI)) GO TO 307
  555 CONTINUE
      DO 556 I=1,NDIL
      IF(DILI.LT.DILUT(I)) THEN
         DO 557 J=NDIL,I,-1
         DILUT(J+1)=DILUT(J)
  557    CONTINUE
         DILUT(I)=DILI
         NDIL=NDIL+1
         IF(NDIL.GT.MAXDIL) CALL XABORT('LIBDI2: MAXDIL IS TOO SMALL.')
         GO TO 307
      ENDIF
  556 CONTINUE
      NDIL=NDIL+1
      IF(NDIL.GT.MAXDIL) CALL XABORT('LIBDI2: MAXDIL IS TOO SMALL.')
      DILUT(NDIL)=DILI
  307 CONTINUE
  670 CONTINUE
  680 CONTINUE
      NDIL=NDIL-1
      IF(NDIL.LT.0) CALL XABORT('LIBDI2: UNABLE TO FIND THE TABULATED'
     1 //' DILUTIONS.')
      CALL XDRCLS(NIN)
      IER=KDRCLS(NIN,1)
      IF(IER.LT.0) THEN
         WRITE (HSMG,'(36HLIBDI2: UNABLE TO CLOSE LIBRARY FILE,1X,A16,
     1   1H.)') NAMFIL
         CALL XABORT(HSMG)
      ENDIF
      RETURN
      END
