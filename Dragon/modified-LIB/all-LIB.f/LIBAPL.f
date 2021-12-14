*DECK LIBAPL
      SUBROUTINE LIBAPL (IPLIB,NAMFIL,MAXTRA,NGRO,NBISO,NL,ISONAM,
     1 ISONRF,IPISO,ISHINA,MASKI,TN,SN,SB,IMPX,NGF,NGFR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Transcription of the useful interpolated microscopic cross section
* data from APOLIB-1 to LCM data structures.
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
* IPLIB   pointer to the lattice microscopic cross section library
*         (L_LIBRARY signature).
* NAMFIL  name of the apolib file.
* MAXTRA  available storage for apollo compacted 
*         transfer cross sections.
* NGRO    number of energy groups.
* NBISO   number of isotopes present in the calculation domain.
* NL      number of Legendre orders required in the calculation
*         NL=1 or higher.
* ISONAM  alias name of isotopes.
* ISONRF  library reference name of isotopes.
* IPISO   pointer array towards microlib isotopes.
* ISHINA  self-shielding name.
* MASKI   isotopic mask. Isotope with index I is processed if
*         MASKI(I)=.true.
* TN      temperature of each isotope.
* SN      dilution cross section in each energy group of each.
*         isotope. a value of 1.0E10 is used for infinite dilution.
* SB      dilution cross section as used by Livolant and Jeanpierre
*         normalization.
* IMPX    print flag.
*
*Parameters: output
* NGF     number of fast groups without self-shielding.
* NGFR    number of fast and resonance groups.
*
*Reference:
*  A. Hoffmann, F. Jeanpierre, A. Kavenoky, M. Livolant AND H. Lorain,
* 'APOLLO - Code multigroupe de resolution de l'equation du transport
*  pour les neutrons thermiques et rapides', Rapport SERMA 'T' No.
*  1 193, Commissariat a l'Energie Atomique, Saclay (1973).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      PARAMETER(MAXQUA=11,MAXDIL=60)
      CHARACTER*(*) NAMFIL
      TYPE(C_PTR) IPLIB,IPISO(NBISO)
      INTEGER MAXTRA,NGRO,NBISO,NL,ISONAM(3,NBISO),ISONRF(3,NBISO),
     1 ISHINA(3,NBISO),IMPX,NGF,NGFR
      REAL TN(NBISO),SN(NGRO,NBISO),SB(NGRO,NBISO)
      LOGICAL MASKI(NBISO)
*----
*  LOCAL VARIABLES
*----
      CHARACTER FORM*4,HVEC(5)*6,HSMG*131,HNISOR*12,HSHI*12,HNAMIS*12
      PARAMETER (NSYSO=6,MAXIT=1000,MAXVEC=11,MAXTMP=40)
      TYPE(C_PTR) KPLIB
      LOGICAL NOTG,LEXC,LALL,LALL2,LALBIS
      DOUBLE PRECISION X1,X2,DDE,ENER,TMP
      INTEGER IANIS(80),ITY(80),NEXT(80),NEXU(80),NEXV(80),NEXW(80),
     1 III(80),IT(MAXIT),ITYPE(MAXVEC),ITYSEC(MAXVEC),TIT(18),NTETA(4),
     2 NSE(4)
      REAL TETAB(MAXTMP),SIGE(MAXDIL,4),SEAUX(MAXDIL,150),XE(MAXDIL),
     1 GE(MAXDIL)
      EQUIVALENCE(AA,NN)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NISB,NISBEF,ITYPRO
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IPR
      REAL, ALLOCATABLE, DIMENSION(:) :: AWR,VECT,SIG1,SIGA,SIGF,
     1 PRI,VTHER,SSS,SSS1,SS1,SS11,UUU,DELTA
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SIGS,PHI,PP,PP1
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SCAT,SEFF
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: LINF
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IPR(4,NBISO),NISB(NBISO),NISBEF(NBISO),ITYPRO(NL))
      ALLOCATE(AWR(NBISO),VECT(NGRO+1),SIGS(NGRO,NL),SIG1(NGRO),
     1 SIGA(NGRO),SIGF(NGRO),PRI(MAXTRA),VTHER(NGRO),PHI(NGRO,MAXQUA),
     2 PP(NGRO,NGRO+1),PP1(NGRO,NGRO+1),SSS(NGRO),SSS1(NGRO),SS1(NGRO),
     3 SS11(NGRO),UUU(NGRO),DELTA(NGRO),SCAT(NGRO,NGRO,NL),
     4 SEFF(MAXDIL,NGRO,4))
      ALLOCATE(LINF(NGRO))
*
      IQUAN=0
      X1=0.0D0
      X2=0.0D0
      NGF=NGRO+1
      NGFR=0
      DO 10 IMX=1,NBISO
      IPR(1,IMX)=0
      HSHI=' '
      IF(MASKI(IMX)) THEN
         WRITE(HSHI,'(3A4)') (ISHINA(I0,IMX),I0=1,3)
         WRITE(HNISOR,'(3A4)') (ISONRF(I0,IMX),I0=1,3)
         I=INDEX(HNISOR,' ')
         IF(I.EQ.0) THEN
            READ(HNISOR,'(I8)') NISB(IMX)
         ELSE
            WRITE(FORM,'(2H(I,I1,1H))') I-1
            READ(HNISOR,FORM) NISB(IMX)
         ENDIF
         I=INDEX(HSHI,' ')
         IF(HSHI.EQ.' ') THEN
            NISBEF(IMX)=0
         ELSE IF(I.EQ.0) THEN
            READ(HSHI,'(I8)') NISBEF(IMX)
         ELSE
            WRITE(FORM,'(2H(I,I1,1H))') I-1
            READ(HSHI,FORM) NISBEF(IMX)
         ENDIF
      ENDIF
      IF(HSHI.EQ.' ') THEN
         IPR(2,IMX)=1
      ELSE
         IPR(2,IMX)=0
      ENDIF
      IPR(3,IMX)=0
      IPR(4,IMX)=0
   10 CONTINUE
      IF(IMPX.GT.0) WRITE(NSYSO,890) NAMFIL
      NIN=KDROPN(NAMFIL,2,2,0)
      IF(NIN.LE.0) THEN
         WRITE(HSMG,'(36HLIBAPL: UNABLE TO OPEN LIBRARY FILE ,A16,
     1   6H. NIN=,I4,1H.)') NAMFIL,NIN
         CALL XABORT(HSMG)
      ENDIF
*----
*  RECOVER THE GROUP STRUCTURE
*----
   20 READ(NIN) INDLOR,NR,NIT,(IT(I),I=1,NIT)
      IF(INDLOR.EQ.9999) THEN
         WRITE(NSYSO,940)
         CALL LCMGET(IPLIB,'DELTAU',DELTA)
         CALL LCMGET(IPLIB,'ENERGY',VECT)
         E0=1.0E-6*VECT(1)
         DO 25 I=1,NGRO
         UUU(I)=LOG(VECT(1)/VECT(I+1))
   25    CONTINUE
      ELSE IF(IT(3).EQ.0) THEN
         DO 30 K=1,NR
         READ(NIN)
   30    CONTINUE
         GO TO 20
      ELSE
         READ(NIN) E0,DEL,(UUU(I),I=1,NGRO),(DELTA(I),I=1,NGRO)
         NR1=NR-1
         VECT(1)=1.0E6*E0
         DO 40 I=1,NGRO
         VECT(I+1)=1.0E6*E0*EXP(-UUU(I))
   40    CONTINUE
         CALL LCMPUT(IPLIB,'ENERGY',NGRO+1,2,VECT)
         CALL LCMPUT(IPLIB,'DELTAU',NGRO,2,DELTA)
      ENDIF
*----
*  ***MATERIAL/ISOTOPE LOOP***
*----
      NED=0
      LALBIS=.TRUE.
   45 NOTG=.TRUE.
      REWIND(NIN)
      NTITLE=18
   50 READ(NIN) INDLOR,NR,NIT,(IT(I),I=1,NIT),(TIT(I),I=1,NTITLE)
      IF(NIT.GT.MAXIT) THEN
         WRITE(HSMG,960) 'MAXIT'
         CALL XABORT(HSMG)
      ENDIF
      IF(INDLOR.EQ.9999) GO TO 740
      IMAIL=IT(1)
      IF(IMAIL.EQ.99) THEN
         IX=47
      ELSE IF(IMAIL.EQ.142) THEN
         IX=60
      ELSE IF(IMAIL.EQ.172) THEN
         IX=80
      ELSE
         WRITE(HSMG,'(45HLIBAPL: INCONSISTENT GROUP STRUCTURES. IT(1)=,
     1   I5)') IMAIL
         CALL XABORT(HSMG)
      ENDIF
      LALL=.TRUE.
      LALL2=LALBIS
      DO 70 IMX=1,NBISO
      IF(MASKI(IMX)) THEN
         IMT=IMX
         LALL=LALL.AND.(IPR(1,IMX).EQ.1).AND.(IPR(2,IMX).EQ.1)
         LALL2=LALL2.AND.(IPR(1,IMX).EQ.1)
         IF((INDLOR.EQ.NISB(IMX)).AND.(IPR(1,IMX).EQ.0)) GO TO 90
         IF((INDLOR.EQ.NISBEF(IMX)).AND.(IPR(1,IMX).EQ.1).AND.
     1   (IPR(2,IMX).EQ.0)) GO TO 500
      ENDIF
   70 CONTINUE
      IF(LALL) THEN
         GO TO 740
      ELSE IF(LALL2) THEN
         LALBIS=.FALSE.
         GO TO 45
      ELSE
         DO 80 K=1,NR
         READ(NIN)
   80    CONTINUE
         GO TO 50
      ENDIF
*----
*  MATERIAL CONTROL
*----
   90 IPR(1,IMT)=1
      NOTG=.FALSE.
      NR1=NR
      IF(IT(2).NE.NGRO) CALL XABORT('LIBAPL: INCONSISTENT GROUP STRUC'
     1 //'TURES.')
      IF(IT(3).NE.0) THEN
         READ(NIN)
         NR1=NR-1
      ENDIF
      NRST=IT(4)
      KPLIB=IPISO(IMT) ! set IMT-th isotope
      DO 106 J=1,NGRO
      SIG1(J)=0.0
      SIGA(J)=0.0
      SIGF(J)=0.0
      DO 105 IL=1,NL
      SIGS(J,IL)=0.0
  105 CONTINUE
  106 CONTINUE
      NTYPE=0
      NS1=0
      DO 205 IRST=1,NRST
      IF(IRST.GT.1) THEN
         IQUAN=1
         ITYSEC(1)=IT(4+NS1+IRST)
      ELSE IF(IT(5).GT.0) THEN
         IQUAN=1
         ITYSEC(1)=IT(5)
      ELSE IF(IT(5).EQ.0) THEN
         IQUAN=4
         ITYSEC(1)=1
         ITYSEC(2)=2
         ITYSEC(3)=3
         ITYSEC(4)=4
      ELSE IF(IT(5).LT.0) THEN
         NS1=-IT(5)
         IQUAN=NS1
         DO 110 I=1,IQUAN
         ITYSEC(I)=IT(5+I)
  110    CONTINUE
      ENDIF
      IF(IQUAN.GT.MAXQUA) CALL XABORT('LIBAPL: MAXQUA TOO SMALL.')
      READ(NIN)((PHI(J,ISEC),J=1,NGRO),ISEC=1,IQUAN)
      NR1=NR1-1
      DO 200 ISEC=1,IQUAN
      IMMOND=ITYSEC(ISEC)
      DO 120 I=1,NTYPE
      IF(IMMOND.EQ.ITYPE(I)) GO TO 200
  120 CONTINUE
      IF(IMPX.GT.7) THEN
         WRITE(NSYSO,920) NISB(IMT),IRST,IMMOND
         WRITE(NSYSO,930) (PHI(J,ISEC),J=1,NGRO)
      ENDIF
      NTYPE=NTYPE+1
      ITYPE(NTYPE)=IMMOND
      IF(IMMOND.EQ.1) THEN
         DO 140 J=1,NGRO
         SIGS(J,1)=PHI(J,ISEC)
  140    CONTINUE
      ELSE IF(IMMOND.EQ.2) THEN
         DO 150 J=1,NGRO
         SIGA(J)=PHI(J,ISEC)
  150    CONTINUE
      ELSE IF(IMMOND.EQ.3) THEN
         CALL LCMPUT(KPLIB,'NUSIGF',NGRO,2,PHI(1,ISEC))
      ELSE IF(IMMOND.EQ.4) THEN
         DO 155 J=1,NGRO
         PHI(J,ISEC)=PHI(J,ISEC)*DELTA(J)
  155    CONTINUE
         CALL LCMPUT(KPLIB,'CHI',NGRO,2,PHI(1,ISEC))
      ELSE IF(IMMOND.EQ.5) THEN
         CALL LCMPUT(KPLIB,'NG',NGRO,2,PHI(1,ISEC))
         IPR(3,IMT)=1
         DO 160 I=1,NED
         IF(HVEC(I).EQ.'NG') GO TO 200
  160    CONTINUE
         NED=NED+1
         HVEC(NED)='NG'
      ELSE IF(IMMOND.EQ.6) THEN
         CALL LCMPUT(KPLIB,'NFTOT',NGRO,2,PHI(1,ISEC))
         IPR(4,IMT)=1
         DO 170 I=1,NED
         IF(HVEC(I).EQ.'NFTOT') GO TO 200
  170    CONTINUE
         NED=NED+1
         HVEC(NED)='NFTOT'
      ELSE IF(IMMOND.EQ.10) THEN
         DO 180 J=1,NGRO
         SIG1(J)=PHI(J,ISEC)
  180    CONTINUE
      ELSE IF(IMMOND.EQ.11) THEN
         DO 185 J=1,NGRO
         VECT(J)=1.0/(3.0*PHI(J,ISEC))
  185    CONTINUE
         CALL LCMPUT(KPLIB,'STRD',NGRO,2,VECT)
         DO 190 I=1,NED
         IF(HVEC(I).EQ.'STRD') GO TO 200
  190    CONTINUE
         NED=NED+1
         HVEC(NED)='STRD'
      ELSE
         WRITE(NSYSO,920) NISB(IMT),IRST,IMMOND
         CALL XABORT('LIBAPL: UNKNOWN REACTION TYPE.')
      ENDIF
  200 CONTINUE
  205 CONTINUE
*----
*  SCATTERING MATRIX CONTROL
*----
      ITH=0
      IMAT=0
      IC=5+NS1+NRST
      NRSTR=IT(IC)
      ICC=IC+6*NRSTR+1
      NN=IT(ICC)
      AWR(IMT)=AA
      NN=IT(ICC+1)
      AT=AA
      NKDEB=1
      IIIC=1
*
      IF(NRSTR.EQ.0)GO TO 380
      IC=IC-5
      ITH=0
      IMAT=0
      IMAT1=0
      ITROUV=0
      IMAT0=0
      DO 290 IS=1,NRSTR
      IC=IC+6
      IF(IT(IC).GT.NL-1) GO TO 280
      IF(IT(IC+1).EQ.7) THEN
*        DO TEMPERATURE INTERPOLATION FOR THE THERMAL TRANSFER MATRICES.
         IF(IX.NE.IT(IC+3)) THEN
            WRITE(NSYSO,950) IX,IT(IC+3)
            IX=IT(IC+3)
         ENDIF
         IF(IT(IC).EQ.1) IMAT1=IMAT1+1
         IF(IMAT1.EQ.1) THEN
            IMAT=0
            ITROUV=0
         ENDIF
         ITH=1
         TEMPI=REAL(IT(IC+5))+0.16
         TEMPA=TN(IMT)
         IMAT=IMAT+1
         IF(ITROUV.NE.0) GO TO 280
         ITEST=(IMAT/2)*2-IMAT
         IF(ITEST.NE.0) THEN
            READ(NIN)((PP1(K,J),K=1,IX),J=1,IX),(SSS1(K),K=1,IX),
     1      (SS11(K),K=1,IX)
            NR1=NR1-1
            X1=TEMPI
         ELSE
            IF(IT(IC).EQ.0) THEN
               READ(NIN)((PP(K,J),K=1,IX),J=1,IX),(SSS(K),K=1,IX),
     1         (SS1(K),K=1,IX)
               NR1=NR1-1
               X2=TEMPI
            ELSE
               READ(NIN)((PP(K,J),K=1,IX),J=1,IX),(SSS1(K),K=1,IX),
     1         (SS11(K),K=1,IX)
               NR1=NR1-1
               X2=TEMPI
            ENDIF
         ENDIF
         IF(IMAT.EQ.1)GO TO 290
         XX=REAL((TEMPA-X1)*(TEMPA-X2))
         IF(XX.LE.0.)ITROUV=1
         IF((TEMPA.LE.TEMPI).AND.(IMAT.EQ.2))ITROUV=1
         IF(IT(IC+6).EQ.1.AND.IT(IC).NE.1)IMAT0=IMAT
         IF(IMAT.EQ.IMAT0)ITROUV=1
         IF(ITROUV.NE.1)GO TO 290
         XX=REAL((TEMPA-X1)/(X1-X2))
         IF(IMAT.EQ.1)XX=0.
         I2=IIIC+IX*IX-1
         IF(I2.GT.MAXTRA) THEN
            WRITE(HSMG,960) 'MAXTRA'
            CALL XABORT(HSMG)
         ENDIF
         IF(IT(IC).EQ.0) THEN
            DO 215 K=1,IX
            SSS(K)=SSS1(K)+(SSS1(K)-SSS(K))*XX
            KI=NGRO-K+1
            SS1(K)=SS11(K)+(SS11(K)-SS1(K))*XX
            DO 210 J=1,IX
            KJ=NGRO-J+1
            PP(J,K)=PP1(J,K)+(PP1(J,K)-PP(J,K))*XX
  210       CONTINUE
  215       CONTINUE
            DO 225 J=1,IX
            DO 220 K=1,IX
            PRI(IIIC+(J-1)*IX+K-1)=PP(J,K)
  220       CONTINUE
  225       CONTINUE
         ELSE
*           ANISOTROPES
            DO 245 K=1,IX
            KI=NGRO-K+1
            DO 240 J=1,IX
            KJ=NGRO-J+1
            PP1(J,K)=PP1(J,K)+(PP1(J,K)-PP(J,K))*XX
  240       CONTINUE
  245       CONTINUE
            DO 255 J=1,IX
            DO 250 K=1,IX
            PRI(IIIC+(J-1)*IX+K-1)=PP1(J,K)
  250       CONTINUE
  255       CONTINUE
         ENDIF
         IANIS(NKDEB)=IT(IC)
         ITY(NKDEB)=7
         NEXT(NKDEB)=IX*IX
         NEXU(NKDEB)=IX
         NEXV(NKDEB)=IX
         NEXW(NKDEB)=INT(TEMPA)
      ELSE
         IANIS(NKDEB)=IT(IC)
         ITY(NKDEB)=IT(IC+1)
         NEXT(NKDEB)=IT(IC+2)
         NEXU(NKDEB)=IT(IC+3)
         NEXV(NKDEB)=IT(IC+4)
         NEXW(NKDEB)=IT(IC+5)
         I2=IIIC+NEXT(NKDEB)-1
         IF(I2.GE.IIIC) THEN
            IF(I2.GT.MAXTRA) THEN
               WRITE(HSMG,960) 'MAXTRA'
               CALL XABORT(HSMG)
            ENDIF
            READ(NIN)(PRI(J),J=IIIC,I2)
            NR1=NR1-1
         ENDIF
      ENDIF
      III(NKDEB)=IIIC
      IIIC=I2+1
      NKDEB=NKDEB+1
      GOTO 290
  280 READ(NIN)
      NR1=NR1-1
  290 CONTINUE
*----
*  FREE GAS THERMAL DIFFUSION MATRICES.
*----
      IF(IX.EQ.0)GO TO 380
      IF(ITH.NE.0)GO TO 360
      T=TN(IMT)/293.16
      AMT=AWR(IMT)
      IF(AMT.LT.1.0)AMT=1.0
      DO 300 K=NGRO-IX+1,NGRO
      SIG1(K)=SIG1(K)/SIGS(K,1)
  300 CONTINUE
      X1=0.0253D-06
      DDE=-UUU(NGRO-IX)
      ENER=E0*EXP(DDE)
      X2=SQRT(ENER/X1)
      DO 305 J=1,IX
      K=IX+1-J
      IE=NGRO-IX+J
      DDE=-UUU(IE)
      ENER=E0*DEXP(DDE)
      DDE=SQRT(ENER/X1)
      VECT(K)=REAL(X2-DDE)
      X2=DDE
      VTHER(K)=2.0*VECT(K)/DELTA(IE)
  305 CONTINUE
      CALL LIBBAS(1,AT,0.0,AMT,T,IX,VTHER,VECT,NGRO,PP,SSS,SSS1,SS11)
      IF(AMT.GT.100.) THEN
         DO 310 J=1,IX
         K=NGRO-J+1
         SSS(J)=SIGS(K,1)
  310    CONTINUE
      ENDIF
      DO 335 I=1,IX
      K=NGRO-I+1
      SIG1(K)=SIG1(K)*SSS(I)
      SS1(I)=SIG1(K)
      RENORM=0.0
      DO 320 J=1,IX
      RENORM=RENORM+VTHER(J)*VECT(J)*PP(J,I)
  320 CONTINUE
      RENORM=RENORM/(VTHER(I)*VECT(I))
      RENORM=1.0/RENORM
      DO 330 J=1,IX
      PP(J,I)=PP(J,I)*SSS(I)*RENORM
  330 CONTINUE
  335 CONTINUE
      DO 345 J=1,IX
      AUX=VTHER(J)*VTHER(J)
      DO 340 I=1,IX
      PP(I,J)=PP(I,J)/AUX*VTHER(I)*VTHER(I)
  340 CONTINUE
  345 CONTINUE
      I2=IIIC+IX*IX
      IF(I2.GT.MAXTRA) THEN
         WRITE(HSMG,960) 'MAXTRA'
         CALL XABORT(HSMG)
      ENDIF
      DO 355 J=1,IX
      DO 350 K=1,IX
      PRI(IIIC+(J-1)*IX+K-1)=PP(J,K)
  350 CONTINUE
  355 CONTINUE
      IANIS(NKDEB)=0
      ITY(NKDEB)=7
      NEXT(NKDEB)=IX*IX
      NEXU(NKDEB)=IX
      NEXV(NKDEB)=IX
      NEXW(NKDEB)=INT(TN(IMT))
      III(NKDEB)=IIIC
      IIIC=I2+1
      NKDEB=NKDEB+1
*
  360 DO 370 J=1,IX
      K=NGRO-J+1
      SIG1(K)=SS1(J)
      SIGS(K,1)=SSS(J)
  370 CONTINUE
*
  380 IF(NR1.GT.0) THEN
         DO 390 IR=1,NR1
         READ(NIN)
  390    CONTINUE
      ENDIF
      NKDEB=NKDEB-1
      IIIC=IIIC-1
      IF(IMPX.GT.0) THEN
         WRITE(NSYSO,860) (ISONAM(I0,IMT),I0=1,3),(TIT(J),J=1,9),
     1   NISBEF(IMT),IIIC,(ITYPE(L),L=1,NTYPE)
         WRITE(NSYSO,870) (ITY(L),L=1,NKDEB)
         WRITE(NSYSO,880) (TIT(J),J=10,18)
      ENDIF
      IF(IMPX.GT.7) THEN
         DO 395 K1=1,NKDEB
         I1=III(K1)
         I2=I1+NEXT(K1)-1
         WRITE(NSYSO,910) NISB(IMT),K1,ITY(K1),NEXU(K1),NEXV(K1),
     1   NEXW(K1),IANIS(K1),(PRI(K),K=I1,I2)
  395    CONTINUE
      ENDIF
*----
*  SAVE SCATTERING MATRICES ON LCM
*----
      INGRO=0
      DO 396 IG1=1,NGRO
      IF(SIG1(IG1).NE.0.0) INGRO=NL-1
  396 CONTINUE
      DO 480 IL=0,INGRO
      ZL=2.0*REAL(IL)+1.0
      DO 420 IG2=1,NGRO
      CALL LIBSEC(MAXTRA,IG2,IL,NGRO,IX,UUU,DELTA,SIGS(1,1),SIG1,PRI,
     1 NLET,VECT,DEL,NKDEB,IANIS,ITY,NEXT,NEXU,NEXV,NEXW,III)
      DO 400 IG1=1,IG2
      SCAT(IG2,IG1,IL+1)=VECT(IG2-IG1+1)*DELTA(IG2)/(ZL*DELTA(IG1))
  400 CONTINUE
      DO 410 IG1=IG2+1,NGRO
      SCAT(IG2,IG1,IL+1)=VECT(IG2+NGRO-IG1+1)*DELTA(IG2)/(ZL*DELTA(IG1))
  410 CONTINUE
  420 CONTINUE
*
      IF(IL.EQ.0) THEN
*        PROCESS NEXCESS INFORMATION.
         LEXC=.FALSE.
         DO 430 IG1=1,NGRO-IX
         SSS(IG1)=-SIGS(IG1,1)
         DO 425 IG2=1,NGRO
         SSS(IG1)=SSS(IG1)+SCAT(IG2,IG1,1)
  425    CONTINUE
         IF(SSS(IG1)/SIGS(IG1,1).GT.1.0E-5) THEN
            LEXC=.TRUE.
            SIGS(IG1,1)=SIGS(IG1,1)+SSS(IG1)
         ELSE
            SSS(IG1)=0.0
         ENDIF
  430    CONTINUE
         DO 440 IG1=NGRO-IX+1,NGRO
         SSS(IG1)=0.0
  440    CONTINUE
         IF(LEXC) CALL LCMPUT(KPLIB,'N2N',NGRO,2,SSS)
      ENDIF
*
      IF(IL.GT.0) THEN
         DO 455 IG1=1,NGRO
         SIGS(IG1,IL+1)=0.0
         DO 450 IG2=1,NGRO
         SIGS(IG1,IL+1)=SIGS(IG1,IL+1)+SCAT(IG2,IG1,IL+1)
  450    CONTINUE
  455    CONTINUE
      ENDIF
  480 CONTINUE
      DO 490 IG1=1,NGRO
      VECT(IG1)=SIGA(IG1)+SIGS(IG1,1)-SSS(IG1)
  490 CONTINUE
*----
*  SAVE INFINITE-DILUTION X-S INFORMATION.
*----
      WRITE(HNAMIS,'(3A4)') (ISONAM(I0,IMT),I0=1,3)
      CALL LCMPUT(KPLIB,'NTOT0',NGRO,2,VECT)
      CALL LCMPUT(KPLIB,'README',18,3,TIT)
      CALL LCMPTC(KPLIB,'ALIAS',12,1,HNAMIS)
      CALL LCMPUT(KPLIB,'AWR',1,2,AWR(IMT))
      CALL XDRLGS(KPLIB,1,0,0,INGRO,1,NGRO,SIGS,SCAT,ITYPRO)
      GO TO 50
*----
*  SELF-SHIELDING CONTROL.
*----
  500 IPR(2,IMT)=1
      IF((IT(1).NE.IMAIL).OR.(IT(2).NE.NGRO).OR.(IT(3).NE.0))
     1 CALL XABORT('LIBAPL: SELF-SHIELDING FAILURE (1).')
      NS=IT(4)
*----
*  RECOVER INFINITE-DILUTION X-S INFORMATION.
*----
      KPLIB=IPISO(IMT) ! set IMT-th isotope
      CALL LCMGET(KPLIB,'NTOT0',SIGA)
      CALL XDRLGS(KPLIB,-1,0,0,NL-1,1,NGRO,SIGS,SCAT,ITYPRO)
*----
*  COMPUTE P0 TRANSFER PROBABILITIES.
*----
      DO 515 IG2=1,NGRO
      SIGA(IG2)=SIGA(IG2)-SIGS(IG2,1)
      DO 510 IG1=1,NGRO
      SCAT(IG2,IG1,1)=SCAT(IG2,IG1,1)/SIGS(IG1,1)
  510 CONTINUE
  515 CONTINUE
*
      ISS=0
      IAS=0
      IFS=0
      I104=0
      JTYSEC=0
      NTYPE=0
      DO 520 IK=1,NS
      IF(IT(IK+4).NE.JTYSEC) THEN
         NTYPE=NTYPE+1
         IF(NTYPE.GT.4) CALL XABORT('LIBAPL: TOO MANY TYPES.')
         NTETA(NTYPE)=1
         IF(IT(IK+4).EQ.101) ISS=NTYPE
         IF(IT(IK+4).EQ.102) IAS=NTYPE
         IF(IT(IK+4).EQ.103) IFS=NTYPE
         IF(IT(IK+4).EQ.104) I104=NTYPE
         JTYSEC=IT(IK+4)
         ITYPE(NTYPE)=IT(IK+4)
      ELSE
         NTETA(NTYPE)=NTETA(NTYPE)+1
      ENDIF
  520 CONTINUE
      IF(IFS.GT.0) CALL LCMGET(KPLIB,'NUSIGF',SIGF)
      IF(IAS.EQ.0) CALL XABORT('LIBAPL: SELF-SHIELDING FAILURE (2).')
      IF(IMPX.GT.0) THEN
         WRITE(NSYSO,990) NISBEF(IMT),(TIT(I),I=1,9),
     1   (ITYPE(I),I=1,NTYPE)
         WRITE(NSYSO,880) (TIT(I),I=10,18)
      ENDIF
*----
*  TEMPERATURE INTERPOLATION OF EFFECTIVE REACTION RATES.
*----
      DO 590 I=1,NTYPE
      IF(NTETA(I).EQ.1) THEN
         READ (NIN) TEMP,NSEI,(SIGE(K,I),K=1,NSEI),N2,N6,((SEFF(K,J
     1   ,I),J=1,N6),K=1,NSEI)
      ELSE
         IF(NTETA(I).GT.MAXTMP) THEN
            WRITE(HSMG,960) 'MAXTMP'
            CALL XABORT(HSMG)
         ENDIF
         DO 532 ITET=1,NTETA(I)
         READ(NIN) TETAB(ITET),NSEI,(SIGE(K,I),K=1,NSEI),N2,N6,
     1   ((SEAUX(K,J),J=1,N6),K=1,NSEI)
         IF(ITET.NE.1) THEN
            IF(TN(IMT).LT.TETAB(ITET)) GO TO 540
            IF(ITET.EQ.NTETA(I)) GO TO 560
         ENDIF
         DO 531 K=1,NSEI
         DO 530 J=1,N6
         SEFF(K,J,I)=SEAUX(K,J)
  530    CONTINUE
  531    CONTINUE
  532    CONTINUE
  540    ITE=ITET+1
         DO 550 ITT=ITE,NTETA(I)
         READ(NIN)
  550    CONTINUE
  560    DT=SQRT(TN(IMT))-SQRT(TETAB(ITET))
         DT=DT/(SQRT(TETAB(ITET))-SQRT(TETAB(ITET-1)))
         DO 575 K=1,NSEI
         DO 570 J=1,N6
         SEFF(K,J,I)=(SEAUX(K,J)-SEFF(K,J,I))*DT+SEAUX(K,J)
  570    CONTINUE
  575    CONTINUE
      ENDIF
      IF(NSEI.GT.MAXDIL) THEN
         WRITE(HSMG,'(37HLIBAPL: MAXDIL SHOULD BE INCREASED TO,I4)')
     1   NSEI
         CALL XABORT(HSMG)
      ELSE IF(NSEI.GT.1) THEN
         IF(SIGE(1,I).GT.SIGE(2,I)) CALL XABORT('LIBAPL: INVALID ORDER'
     1   //'ING OF THE DILUTIONS.')
      ENDIF
      NGF=MIN(NGF,N2-1)
      NGFR=MAX(NGFR,N2+N6-1)
      IF(I.EQ.I104) THEN
         DO 585 J=1,N6
         DO 580 K=1,NSEI
         IF((SIGE(K,I).LT.1.0E10).OR.(K.EQ.1)) THEN
            SEFF(K,J,I)=(1.0-SEFF(K,J,I))*SIGE(K,I)
         ELSE
            SEFF(K,J,I)=SEFF(K-1,J,I)
         ENDIF
  580    CONTINUE
  585    CONTINUE
      ENDIF
      NSE(I)=NSEI
  590 CONTINUE
*----
*  DILUTION INTERPOLATION OF EFFECTIVE REACTION RATES.
*----
      DO 600 L=1,NGRO
      IF(ISS.NE.0) PHI(L,ISS)=SIGS(L,1)
      IF(IAS.NE.0) PHI(L,IAS)=SIGA(L)
      IF(IFS.NE.0) PHI(L,IFS)=SIGF(L)
      IF(I104.NE.0) PHI(L,I104)=SIGA(L)
      LINF(L)=.FALSE.
      VECT(L)=SIGS(L,1)
  600 CONTINUE
*
      DO 625 LE=1,N6
      L=LE+N2-1
      SEIM=MAX(0.0,SN(L,IMT))
      DO 620 I=1,NTYPE
      IF(NSE(I).EQ.1) THEN
         PHI(L,I)=SEFF(1,LE,I)
      ELSE
         NSEI=NSE(I)
         IF(SIGE(NSE(I),I).GE.1.0E10) NSEI=NSE(I)-1
         IF(SEIM.LT.SIGE(NSEI,I)) THEN
            DO 610 K=1,NSEI
            XE(K)=SQRT(SIGE(K,I))
            GE(K)=SEFF(K,LE,I)
  610       CONTINUE
            CALL LIBLAG(NSEI,XE,GE,SQRT(SEIM),PHI(L,I))
         ELSE IF(NSE(I).GT.NSEI) THEN
            IF(I.EQ.I104) LINF(L)=.TRUE.
            FAC=SIGE(NSEI,I)/SEIM
            PHI(L,I)=FAC*SEFF(NSEI,LE,I)+(1.0-FAC)*SEFF(NSE(I),LE,I)
         ENDIF
      ENDIF
  620 CONTINUE
  625 CONTINUE
*----
*  RECOVER THE EFFECTIVE FLUX.
*----
      IF(IMPX.GT.4) WRITE(NSYSO,1020)
      DO 630 L=1,NGRO
      SS1(L)=1.0
  630 CONTINUE
      DO 660 L=N2,N2+N6-1
      SEIM=SN(L,IMT)
      IF(SEIM.EQ.0.) CALL XABORT('LIBAPL: SELF-SHIELDING FAILURE (3).')
      IF((IAS.NE.0).AND.(ISS.NE.0)) THEN
*        COMPUTE THE EFFECTIVE FLUX.
         TMP1=0.0D0
         DO 640 IG2=1,N2-1
         TMP1=TMP1+SCAT(L,IG2,1)*PHI(IG2,ISS)*DELTA(IG2)/DELTA(L)
  640    CONTINUE
         IF(TMP1.GT.5.0E-3*PHI(L,ISS)) THEN
*           USE A SIMPLIFIED MODEL.
            AUX=PHI(L,IAS)
         ELSE
*           USE A SLOWING-DOWN BALANCE EQUATION.
            TMP=TMP1
            DO 650 IG2=N2,N2+N6-1
            TMP=TMP+SCAT(L,IG2,1)*PHI(IG2,ISS)*DELTA(IG2)/DELTA(L)
  650       CONTINUE
            AUX=REAL(PHI(L,IAS)+PHI(L,ISS)-TMP)
         ENDIF
      ELSE IF(IAS.NE.0) THEN
*        COMPUTE THE EFFECTIVE FLUX USING A SIMPLIFIED MODEL.
         AUX=PHI(L,IAS)
      ELSE
         AUX=0.0
      ENDIF
*
      IF(SB(L,IMT).GE.1.0E10) THEN
*        USE AN INFINITE DILUTION VALUE.
         ZNPHI=0.0
      ELSE IF((I104.NE.0).AND.LINF(L)) THEN
*        USE AN INTERPOLATED VALUE NEAR INFINITE DILUTION.
         NSEI=NSE(I104)
         IF(SIGE(NSE(I104),I104).GE.1.0E10) NSEI=NSE(I104)-1
         FAC=(SIGE(NSEI,I104)/SEIM)**2
         ZNPHI=FAC*PHI(L,I104)+(1.0-FAC)*AUX
      ELSE IF(I104.NE.0) THEN
*        USE AN INTERPOLATED VALUE.
         ZNPHI=PHI(L,I104)
      ELSE
*        USE A CALCULATED VALUE.
         ZNPHI=AUX
      ENDIF
      PHI0=1.0-ZNPHI/SB(L,IMT)
      IF((PHI0.LE.0.0).OR.(PHI0.GT.1.2)) THEN
         WRITE(HSMG,980) PHI0,L,ZNPHI,SEIM,(ISONAM(I0,IMT),I0=1,3)
         WRITE(NSYSO,'(/1X,A131)') HSMG
      ENDIF
      SS1(L)=PHI0
      IF(IFS.GT.0) SIGF(L)=PHI(L,IFS)/PHI0
      IF(IAS.GT.0) SIGA(L)=PHI(L,IAS)/PHI0
      IF(ISS.GT.0) SIGS(L,1)=PHI(L,ISS)/PHI0
      IF(IMPX.GT.4) WRITE(NSYSO,1010) L,PHI0,SIGF(L),SIGA(L),SIGS(L,1),
     1 SEIM,SB(L,IMT),ZNPHI
  660 CONTINUE
      IF(IMPX.GT.4) WRITE(NSYSO,'(/)')
*
      CALL LCMPUT(KPLIB,'NWT0',NGRO,2,SS1)
*----
*  SELF-SHIELDING OF THE TRANSFERT CROSS SECTIONS.
*----
      IF(ISS.NE.0) THEN
         DO 675 IG1=1,NGRO
         DO 670 IG2=1,NGRO
         SCAT(IG2,IG1,1)=SCAT(IG2,IG1,1)*SIGS(IG1,1)
  670    CONTINUE
  675    CONTINUE
         INGRO=NL-1
         DO 680 IL=NL-1,0,-1
         IF(ITYPRO(IL+1).EQ.0) THEN
            INGRO=INGRO-1
         ELSE
            GO TO 685
         ENDIF
  680    CONTINUE
  685    DO 695 IL=1,NL-1
         IF(ITYPRO(IL+1).GT.0) THEN
            DO 691 IG2=1,NGRO
            SIGS(IG2,IL+1)=SIGS(IG2,IL+1)*SIGS(IG2,1)/VECT(IG2)
            DO 690 IG1=1,NGRO
            SCAT(IG2,IG1,IL+1)=SCAT(IG2,IG1,IL+1)*SIGS(IG1,1)/VECT(IG1)
  690       CONTINUE
  691       CONTINUE
         ENDIF
  695    CONTINUE
*
*        SAVE SELF-SHIELDED X-S INFORMATION.
         CALL XDRLGS(KPLIB,1,0,0,INGRO,1,NGRO,SIGS,SCAT,ITYPRO)
      ENDIF
*----
*  SELF-SHIELDING OF THE RADIATIVE CAPTURE CROSS SECTIONS.
*----
      IF(IPR(3,IMT).EQ.1) THEN
         CALL LCMGET(KPLIB,'NTOT0',SS1)
         DO 700 I=1,NGRO
         SS1(I)=SS1(I)-VECT(I)
  700    CONTINUE
         CALL LCMGET(KPLIB,'NG',VECT)
         DO 710 I=1,NGRO
         IF(SS1(I).EQ.0.0) GO TO 710
         VECT(I)=VECT(I)*SIGA(I)/SS1(I)
  710    CONTINUE
         CALL LCMPUT(KPLIB,'NG',NGRO,2,VECT)
      ENDIF
*----
*  SELF-SHIELDING OF THE FISSION CROSS SECTIONS.
*----
      IF(IFS.NE.0) THEN
         IF(IPR(4,IMT).EQ.1) THEN
            CALL LCMGET(KPLIB,'NUSIGF',SS1)
            CALL LCMGET(KPLIB,'NFTOT',VECT)
            DO 720 I=1,NGRO
            IF(SS1(I).EQ.0.0) GO TO 720
            VECT(I)=VECT(I)*SIGF(I)/SS1(I)
  720       CONTINUE
            CALL LCMPUT(KPLIB,'NFTOT',NGRO,2,VECT)
         ENDIF
         CALL LCMPUT(KPLIB,'NUSIGF',NGRO,2,SIGF)
      ENDIF
*
      DO 730 I=1,NGRO
      SIGA(I)=SIGA(I)+SIGS(I,1)
  730 CONTINUE
      CALL LCMPUT(KPLIB,'NTOT0',NGRO,2,SIGA)
      GO TO 50
*----
*  CHECK IF ALL NBISO ISOTOPES HAVE BEEN PROCESSED.
*----
  740 NISOT=0
      DO 750 IMT=1,NBISO
      IF(MASKI(IMT)) THEN
         IF((IPR(1,IMT).EQ.0).AND.(.NOT.NOTG)) THEN
            GO TO 45
         ELSE IF((IPR(1,IMT).EQ.0).AND.NOTG) THEN
            WRITE(NSYSO,900) (ISONAM(I0,IMT),I0=1,3),NAMFIL
            NISOT=NISOT+1
         ELSE IF((IPR(2,IMT).EQ.0).AND.(.NOT.NOTG)) THEN
            GO TO 45
         ELSE IF((IPR(2,IMT).EQ.0).AND.NOTG) THEN
            WRITE(NSYSO,900) (ISHINA(I0,IMT),I0=1,3),NAMFIL
            NISOT=NISOT+1
         ENDIF
      ENDIF
  750 CONTINUE
*----
*  ADD NG CROSS SECTIONS.
*----
      DO 790 IMT=1,NBISO
      IF(MASKI(IMT).AND.(IPR(3,IMT).EQ.0)) THEN
         KPLIB=IPISO(IMT) ! set IMT-th isotope
         CALL LCMGET(KPLIB,'NTOT0',VECT)
         CALL LCMLEN(KPLIB,'SIGS00',LENGT,ITYLCM)
         IF(LENGT.EQ.NGRO) THEN
            CALL LCMGET(KPLIB,'SIGS00',SSS)
            DO 760 IU=1,NGRO
            VECT(IU)=VECT(IU)-SSS(IU)
  760       CONTINUE
         ENDIF
         IF(IPR(4,IMT).EQ.1) THEN
            CALL LCMGET(KPLIB,'NFTOT',SSS)
            DO 770 IU=1,NGRO
            VECT(IU)=VECT(IU)-SSS(IU)
  770       CONTINUE
         ENDIF
         CALL LCMLEN(KPLIB,'N2N',LENGT,ITYLCM)
         IF(LENGT.EQ.NGRO) THEN
            CALL LCMGET(KPLIB,'N2N',SSS)
            DO 780 IU=1,NGRO
            VECT(IU)=VECT(IU)+SSS(IU)
  780       CONTINUE
         ENDIF
         CALL LCMPUT(KPLIB,'NG',NGRO,2,VECT)
      ENDIF
  790 CONTINUE
*----
*  CLOSE THE APOLIB FILE.
*----
      IER=KDRCLS(NIN,1)
      IF(IER.LT.0) THEN
         WRITE(HSMG,'(37HLIBAPL: UNABLE TO CLOSE LIBRARY FILE ,A16,1H.
     1   )') NAMFIL
         CALL XABORT(HSMG)
      ENDIF
      IF((IMPX.GT.0).AND.(NED.GT.0)) WRITE(NSYSO,1030) (HVEC(I),
     1 I=1,NED)
      IF(NISOT.GT.0) CALL XABORT('LIBAPL: MISSING ISOTOPES')
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(LINF)
      DEALLOCATE(SEFF,SCAT,DELTA,UUU,SS11,SS1,SSS1,SSS,PP1,PP,PHI,
     1 VTHER,PRI,SIGF,SIGA,SIG1,SIGS,VECT,AWR)
      DEALLOCATE(ITYPRO,NISBEF,NISB,IPR)
      RETURN
*
  860 FORMAT(1X,3A4,3H * ,9A4,2H *,I10,I12,2X,8I3)
  870 FORMAT(1H+,101X,10I3/(102X,10I3))
  880 FORMAT(14X,2H* ,9A4,2H *)
  890 FORMAT(/35H PROCESSING APOLLO-1 LIBRARY NAMED ,A16,1H.//
     1 55X,14HSELF-SHIELDING,1X,8HTRANSFER/22H ISOTOPE..... LIBRARY ,
     2 7HCONTENT,25(1H.),6X,4HDATA,6X,7HFILL-IN,2X,7HVECTOR ,
     3 5HTYPES,11(1H.),1X,12HMATRIX TYPES,11(1H.)/1X,12(1H-),1X,
     4 40(1H-),1X,14(1H-),1X,8(1H-),1X,23(1H-),1X,23(1H-))
  900 FORMAT(/27H LIBAPL: MATERIAL/ISOTOPE ',3A4,16H' IS MISSING ON ,
     1 17HAPOLIB FILE NAME ,A8,1H.)
  910 FORMAT(//8H ISOTOPE,I12,5X,20HDIFFUSION MATRIX NB.,I3,5X,
     1 6HTYPE =,I3,5X,6HNEXU =,I3,5X,6HNEXV =,I3,5X,6HNEXW =,I4,5X,
     1 12HANISOTROPY =,I3/(1P,10E13.5))
  920 FORMAT(//8H ISOTOPE,I12,5X,8HRECORD =,I3,5X,10HREACTION =,I3)
  930 FORMAT(1X,1P,10E13.5)
  940 FORMAT(/47H LIBAPL: UNABLE TO RECOVER THE GROUP STRUCTURE.)
  950 FORMAT(/53H LIBAPL: *** WARNING *** THE NUMBER OF THERMAL GROUPS,
     1 17H WAS CHANGED FROM,I4,3H TO,I4,1H.)
  960 FORMAT(30HLIBAPL: INSUFFICIENT VALUE OF ,A6,1H.)
  980 FORMAT(47HLIBAPL: *** WARNING *** INVALID VALUE OF PHI0 (,1P,
     1 E11.3,0P,10H) IN GROUP,I4,8H. ZNPHI=,1P E12.3,2X,5HSEIM=,E12.3,
     2 2X,5HISO=',3A4,1H')
  990 FORMAT(1X,I12,3H * ,9A4,2H *,21H SELF-SHIELDING DATA.,4X,8I4)
 1010 FORMAT(5X,I5,1P,8E15.5)
 1020 FORMAT(/5X,'GROUP',11X,'PHI0',10X,'SIGF0',10X,'SIGA0',10X,
     1 'SIGS0',10X,'DILUT',13X,'SB',12X,'ZNPHI')
 1030 FORMAT(/39H EXTRA REACTION EDITS FOUND ON APOLIB: ,5A7)
      END
