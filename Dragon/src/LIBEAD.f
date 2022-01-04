*DECK LIBEAD
      SUBROUTINE LIBEAD (IPLIB,MAXISO,MAXMIX,IMPX,NDEPL,NFISS,NSUPS,
     1 NREAC,NPAR,NBISO,ISONAM,ISONRF,IHLIB,ILLIB,MIX,TN,IEVOL,ITYP,
     2 NCOMB)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Add the missing isotopes from the depletion chain.
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
*Parameters: input/output
* IPLIB   pointer to the lattice microscopic cross section library
*         (L_LIBRARY signature).
* MAXISO  maximum value of nbiso.
* MAXMIX  maximum number of mixtures.
* IMPX    print flag. Equal to zero for no print.
* NDEPL   number of depleting isotopes.
* NFISS   number of fissiles isotopes producing fission products.
* NSUPS   number of non-depleting isotopes producing energy.
* NREAC   maximum number of depletion reactions.
* NPAR    maximum number of parent nuclides in the depletion chain.
* NBISO   old/new number of isotopes present in the calculation
*         domain.
* ISONAM  alias name of isotopes.
* ISONRF  library name of isotopes.
* IHLIB   isotope options.
* ILLIB   xs library index for each isotope.
* MIX     mix number of each isotope (can be zero).
* TN      temperature of each isotope.
* IEVOL   non-depletion mask (=1/2 to suppress/force depletion of an
*         isotope).
* ITYP    isotope type:
*         =1: the isotope is not fissile and not a fission product;
*         =2: the isotope is fissile; =3: is a fission product.
* NCOMB   number of depleting mixtures.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB
      INTEGER MAXISO,MAXMIX,IMPX,NDEPL,NFISS,NSUPS,NREAC,NPAR,
     1 NBISO,ISONAM(3,MAXISO),ISONRF(3,MAXISO),IHLIB(2,MAXISO,4),
     2 ILLIB(MAXISO),MIX(MAXISO),IEVOL(MAXISO),ITYP(MAXISO),NCOMB
      REAL TN(MAXISO)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (IOUT=6)
      CHARACTER TEXT1*12,TEXT2*12,TEXT3*8
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MILVO,IIPAR,KFISS
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IDR,KPAR,HGAR
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(MILVO(MAXMIX),IIPAR(NDEPL),IDR(NREAC,NDEPL),
     1 KPAR(NPAR,NDEPL),KFISS(NFISS),HGAR(3,NDEPL))
*----
*  FIND THE NUMBER OF DEPLETING MIXTURES
*----
      CALL LCMGET(IPLIB,'ISOTOPESDEPL',HGAR)
      IF(NDEPL.GT.MAXISO) CALL XABORT('LIBEAD: TOO MANY DEPLETING ISOT'
     1 //'OPES.')
      CALL LCMGET(IPLIB,'DEPLETE-REAC',IDR)
      CALL LCMGET(IPLIB,'PRODUCE-REAC',KPAR)
      NCOMB=0
      DO 30 ISOT=1,NBISO
      IBM=MIX(ISOT)
      IF(IBM.EQ.0) GO TO 30
      IF((IEVOL(ISOT).NE.1).AND.(ITYP(ISOT).GT.1)) THEN
         DO 10 J=1,NCOMB
         IF(IBM.EQ.MILVO(J)) GO TO 30
   10    CONTINUE
         NCOMB=NCOMB+1
         MILVO(NCOMB)=IBM
         GO TO 30
      ENDIF
      IF((IEVOL(ISOT).EQ.1).OR.(ILLIB(ISOT).EQ.0)) GO TO 30
      DO 20 I=1,NDEPL-NSUPS
      IF((ISONRF(1,ISOT).EQ.HGAR(1,I)).AND.(ISONRF(2,ISOT).EQ.
     1 HGAR(2,I)).AND.(ISONRF(3,ISOT).EQ.HGAR(3,I))) THEN
         ITYP(ISOT)=1
         IF(IEVOL(ISOT).EQ.2) ITYP(ISOT)=3
         IF(MOD(IDR(2,I),100).EQ.3) ITYP(ISOT)=2
         IF(MOD(IDR(2,I),100).EQ.4) ITYP(ISOT)=2
         IF(MOD(IDR(2,I),100).EQ.5) ITYP(ISOT)=3
         DO 15 J=1,NCOMB
         IF(IBM.EQ.MILVO(J)) GO TO 30
   15    CONTINUE
         NCOMB=NCOMB+1
         MILVO(NCOMB)=IBM
         GO TO 30
      ENDIF
   20 CONTINUE
      IEVOL(ISOT)=1
   30 CONTINUE
*----
*  ADD THE MISSING ISOTOPES FROM THE DEPLETION CHAIN
*----
      CALL XDISET(KFISS,NFISS,0)
      DO 35 INUCL=1,NDEPL-NSUPS
      IF(MOD(IDR(2,INUCL),100).EQ.4) THEN
         KDRI=IDR(2,INUCL)/100
         IF(KDRI.GT.NFISS) CALL XABORT('LIBEAD: INVALID NFISS.')
         IF(KDRI.GT.0) KFISS(KDRI)=INUCL
      ENDIF
   35 CONTINUE
      NBOLD=NBISO
      DO 130 ICOMB=1,NCOMB
      IBM=MILVO(ICOMB)
      ITER=0
      IFIRST=0
      DO 36 I=1,NBISO
      IF(MIX(I).EQ.IBM) THEN
         IFIRST=I
         GO TO 40
      ENDIF
   36 CONTINUE
      CALL XABORT('LIBEAD: UNABLE TO FIND A DEPLETING MIXTURE.')
   40 ITER=ITER+1
      IF(ITER.GT.100) CALL XABORT('LIBEAD: UNABLE TO COMPLETE THE BURN'
     1 //'UP CHAINS.')
      NADD=0
      DO 120 INUCL=1,NDEPL-NSUPS
      DO 50 I=1,NBISO
      IF((ISONRF(1,I).EQ.HGAR(1,INUCL)).AND.(ISONRF(2,I).EQ.
     1 HGAR(2,INUCL)).AND.(ISONRF(3,I).EQ.HGAR(3,INUCL)).AND.
     2 (MIX(I).EQ.IBM)) GO TO 120
   50 CONTINUE
      WRITE(TEXT1,'(3A4)') (HGAR(I0,INUCL),I0=1,3)
      I1=INDEX(TEXT1,'_')
      IF(I1.EQ.0) THEN
         TEXT2=TEXT1
       ELSE
         TEXT2=TEXT1(:I1-1)
      ENDIF
      TEXT2(9:12)='    '
      DO 60 I=1,NBISO
      IF(MIX(I).NE.IBM) GO TO 60
      WRITE(TEXT1,'(3A4)') (ISONRF(I0,I),I0=1,3)
      I1=INDEX(TEXT1,'_')
      IF(I1.EQ.0) THEN
         TEXT3=TEXT1(:8)
       ELSE
         TEXT3=TEXT1(:I1-1)
      ENDIF
      IF(TEXT3.EQ.TEXT2(:8)) GO TO 120
   60 CONTINUE
      CALL XDISET(IIPAR,NDEPL-NSUPS,0)
      IF(MOD(IDR(2,INUCL),100).EQ.5) THEN
         DO 70 IFIS=1,NFISS
         IF(KFISS(IFIS).GT.0) IIPAR(KFISS(IFIS))=1
   70    CONTINUE
      ENDIF
      DO 80 IPAR=1,NPAR
      KGAR=KPAR(IPAR,INUCL)
      IF(KGAR.EQ.0) THEN
         GO TO 90
      ELSE
         IIPAR(KGAR/100)=1
      ENDIF
   80 CONTINUE
   90 DO 110 JNUCL=1,NDEPL-NSUPS
      IF(IIPAR(JNUCL).EQ.1) THEN
         NBISOL=NBISO
         DO 100 I=1,NBISOL
         IF((ISONRF(1,I).EQ.HGAR(1,JNUCL)).AND.(ISONRF(2,I).EQ.
     1   HGAR(2,JNUCL)).AND.(ISONRF(3,I).EQ.HGAR(3,JNUCL)).AND.
     2   (MIX(I).EQ.IBM)) THEN
*           A PARENT EXISTS. ADD ONE ISOTOPE IN THE ISOTOPE LIST AND
*           SET ISOTOPE PARAMETERS TO STANDARD VALUES.
            NBISO=NBISO+1
            IF(NBISO.GT.MAXISO) CALL XABORT('LIBEAD: MAXISO TOO SMALL.')
            NADD=NADD+1
            IF(IMPX.GT.8) WRITE(IOUT,'(25H LIBEAD: ADDING ISOTOPE '',
     1      3A4,20H'' TO CHILD ISOTOPE '',3A4,12H'' IN MIXTURE,I5)')
     2      (HGAR(I0,INUCL),I0=1,3),(HGAR(I0,JNUCL),I0=1,3),IBM
*           TEXT2 IS THE NEW ALIAS NAME FOR NBISO-TH ISOTOPE.
            READ(TEXT2,'(3A4)') (ISONAM(I0,NBISO),I0=1,3)
            DO 95 I0=1,3
            ISONRF(I0,NBISO)=HGAR(I0,INUCL)
   95       CONTINUE
            DO 96 I0=1,2
            IHLIB(I0,NBISO,1)=IHLIB(I0,IFIRST,1)
   96       CONTINUE
            ILLIB(NBISO)=ILLIB(IFIRST)
            MIX(NBISO)=IBM
            TN(NBISO)=TN(IFIRST)
            IEVOL(NBISO)=0
            ITYP(NBISO)=1
            IF(MOD(IDR(2,INUCL),100).EQ.3) ITYP(NBISO)=2
            IF(MOD(IDR(2,INUCL),100).EQ.4) ITYP(NBISO)=2
            IF(MOD(IDR(2,INUCL),100).EQ.5) ITYP(NBISO)=3
            GO TO 120
         ENDIF
  100    CONTINUE
      ENDIF
  110 CONTINUE
  120 CONTINUE
      IF(NADD.GT.0) GO TO 40
      IF((IMPX.GT.0).AND.(NBISO-NBOLD.GT.0)) THEN
         WRITE(IOUT,150) NBISO-NBOLD,IBM
      ENDIF
      NBOLD=NBISO
  130 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(HGAR,KFISS,KPAR,IDR,IIPAR,MILVO)
      RETURN
*
  150 FORMAT(/8H LIBEAD:,I5,39H DEPLETING ISOTOPES HAVE BEEN ADDED IN ,
     1 7HMIXTURE,I5,1H.)
      END
