*DECK LIBDRB
      SUBROUTINE LIBDRB (IPDRL,NGRO,NL,NDEL,NBESP,SN,SB,NED,HVECT,DELTA,
     1 LBIN,NFS,BENER,IMPX,NGF,NGFR,LSCAT,LSIGF,LADD,LGOLD,SIGS,SCAT,
     2 TOTAL,ZNPHI,SIGF,CHI,CHI4G,SADD,GOLD,BIN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read an interpolate in dilution one isotope in draglib format at a
* selected temperature.
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
* IPDRL   pointer to the draglib (L_DRAGLIB signature).
* NGRO    number of energy groups.
* NL      number of Legendre orders required in the calculation.
*         NL=1 or higher.
* NDEL    number of delayed precursor groups.
* NBESP   number of energy-dependent fission spectra.
* SN      dilution cross section in each energy group. A value of
*         1.0E10 is used for infinite dilution.
* SB      dilution cross section as used in Livolant and Jeanpierre
*         normalization.
* NED     number of extra vector edits.
* HVECT   names of the extra vector edits.
* DELTA   lethargy widths.
* LBIN    number of fine groups.
* NFS     number of fine groups per coarse group.
* BENER   energy limits of the fine groups.
* IMPX    print flag.
*
*Parameters: input/output
* NGF     number of fast groups without self-shielding.
* NGFR    number of fast and resonance groups.
*
*Parameters: output
* LSCAT   scattering mask (=.true. if a given Legendre order of the
*         scattering cross section exists).
* LSIGF   fission mask (=.true. if the isotope can fission).
* LADD    additional cross section mask (=.true. if a given additional
*         cross section exists).
* LGOLD   Goldstein-Cohen mask (=.true. if Goldstein-Cohen parameters
*         exists).
* SIGS    scattering cross sections.
* SCAT    scattering transfer matrices.
* TOTAL   total cross sections.
* ZNPHI   fluxes.
* SIGF    nu*fission cross sections.
* CHI     fission spectrum.
* CHI4G   energy-dependent fission spectra.
* SADD    additional cross sections.
* GOLD    Goldstein-Cohen parameters.
* BIN     BIN(IGR,1): total fine group cross sections;
*         BIN(IGR,2): isotropic scattering fine group cross sections;
*         BIN(IGR,3): nu*fission fine group cross sections.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      PARAMETER(MAXDIL=50)
      CHARACTER*(*) HVECT(NED)
      TYPE(C_PTR) IPDRL
      INTEGER NGRO,NL,NDEL,NBESP,NED,LBIN,NFS(NGRO),IMPX,NGF,NGFR
      REAL SN(NGRO),SB(NGRO),DELTA(NGRO),BENER(LBIN+1),SIGS(NGRO,NL),
     1 SCAT(NGRO,NGRO,NL),TOTAL(NGRO),ZNPHI(NGRO),SIGF(NGRO,0:NDEL),
     2 CHI(NGRO,0:NDEL),CHI4G(NGRO,NBESP),SADD(NGRO,NED),GOLD(NGRO),
     3 BIN(LBIN,3)
      LOGICAL LSCAT(NL),LSIGF,LADD(NED),LGOLD
*----
*  LOCAL VARIABLES
*----
      CHARACTER CM*2,CD*4,HSMG*131,HNUSIG*12,HCHI*12,HTOTAL*5
      PARAMETER (IOUT=6)
      INTEGER KTOTLR,KSIGFR,KCHIR,KPHIR
      LOGICAL LPCAT
      DOUBLE PRECISION TMP,ZNGAR,SQD,SQ0,SQ1,SQ2,SQ3,FACT1,FACT2
      REAL DILUT(MAXDIL)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NJJ,IJJ,KADDR
      REAL, ALLOCATABLE, DIMENSION(:) :: GAR
      REAL, ALLOCATABLE, DIMENSION(:,:) :: TERP,SIGT
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: LSDIL,LPDIL,LINF
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(NJJ(NGRO),IJJ(NGRO),KADDR(NED))
      ALLOCATE(GAR(NGRO**2),TERP(MAXDIL,NGRO),SIGT(NGRO,MAXDIL))
      ALLOCATE(LSDIL(NL),LPDIL(NL),LINF(NGRO))
*
      CALL XDRSET(TOTAL,NGRO,0.0)
      HTOTAL='NTOT0'
      CALL LCMLEN(IPDRL,'NTOT0',LENGT,ITYLCM)
      IF(LENGT.EQ.0) CALL XABORT('LIBDRB: MISSING TOTAL XS INFO.')
      CALL LCMGET(IPDRL,HTOTAL,TOTAL)
      CALL LCMLEN(IPDRL,'NUSIGF',LENGT,ITYLCM)
      LSIGF=(LENGT.GT.0)
      IF(LSIGF) THEN
         DO 10 IDEL=0,NDEL
         IF(IDEL.EQ.0) THEN
            HNUSIG='NUSIGF'
            HCHI='CHI'
         ELSE
            WRITE(HNUSIG,'(6HNUSIGF,I2.2)') IDEL
            WRITE(HCHI,'(3HCHI,I2.2)') IDEL
            CALL LCMLEN(IPDRL,HNUSIG,ILONG,ITYLCM)
            IF(ILONG.EQ.0) CALL XABORT('LIBDRB: MISSING '//HNUSIG//
     1      ' INFO.')
         ENDIF
         CALL XDRSET(SIGF(1,IDEL),NGRO,0.0)
         CALL LCMGET(IPDRL,HNUSIG,SIGF(1,IDEL))
         IF((NBESP.EQ.0).OR.(IDEL.GT.0)) THEN
            CALL XDRSET(CHI(1,IDEL),NGRO,0.0)
            CALL LCMLEN(IPDRL,HCHI,LENGT,ITYLCM)
            IF(LENGT.GT.0) CALL LCMGET(IPDRL,HCHI,CHI(1,IDEL))
         ENDIF
   10    CONTINUE
         DO 11 ISP=1,NBESP
         WRITE(HCHI,'(5HCHI--,I2.2)') ISP
         CALL XDRSET(CHI4G(1,ISP),NGRO,0.0)
         CALL LCMLEN(IPDRL,HCHI,LENGT,ITYLCM)
         IF(LENGT.GT.0) CALL LCMGET(IPDRL,HCHI,CHI4G(1,ISP))
   11    CONTINUE
      ELSE
         CALL XDRSET(SIGF,NGRO*(1+NDEL),0.0)
      ENDIF
      DO 80 IL=0,NL-1
      DO 16 IG2=1,NGRO
      DO 15 IG1=1,NGRO
      SCAT(IG1,IG2,IL+1)=0.0
   15 CONTINUE
   16 CONTINUE
      CALL XDRSET(SIGS(1,IL+1),NGRO,0.0)
      WRITE (CM,'(I2.2)') IL
      CALL LCMLEN(IPDRL,'SCAT'//CM,LENGT,ITYLCM)
      LPCAT=(LENGT.GT.0)
      IF(LPCAT) THEN
         CALL LCMGET(IPDRL,'NJJS'//CM,NJJ)
         CALL LCMGET(IPDRL,'IJJS'//CM,IJJ)
         LENGT=0
         DO 20 I=1,NGRO
         LENGT=LENGT+NJJ(I)
   20    CONTINUE
         CALL XDRSET(GAR,LENGT,0.0)
         CALL LCMGET(IPDRL,'SCAT'//CM,GAR)
         IGAR=0
*        IG2 IS THE SECONDARY GROUP.
         DO 35 IG2=1,NGRO
         DO 30 IG1=IJJ(IG2),IJJ(IG2)-NJJ(IG2)+1,-1
         IGAR=IGAR+1
         SCAT(IG2,IG1,IL+1)=GAR(IGAR)
   30    CONTINUE
   35    CONTINUE
      ENDIF
      CALL LCMLEN(IPDRL,'SIGS'//CM,LENGT,ITYLCM)
      LSCAT(IL+1)=(LENGT.GT.0)
      IF(LSCAT(IL+1)) THEN
         CALL LCMGET(IPDRL,'SIGS'//CM,SIGS(1,IL+1))
         CALL LCMLEN(IPDRL,'PCAT'//CM,LENGT,ITYLCM)
         IF(.NOT.LPCAT.AND.(LENGT.GT.0)) THEN
            CALL LCMGET(IPDRL,'NJJS'//CM,NJJ)
            CALL LCMGET(IPDRL,'IJJS'//CM,IJJ)
            LENGT=0
            DO 40 I=1,NGRO
            LENGT=LENGT+NJJ(I)
   40       CONTINUE
            CALL XDRSET(GAR,LENGT,0.0)
            CALL LCMGET(IPDRL,'PCAT'//CM,GAR)
            IGAR=0
            DO 46 IG2=1,NGRO
            DO 45 IG1=IJJ(IG2),IJJ(IG2)-NJJ(IG2)+1,-1
            IGAR=IGAR+1
            SCAT(IG2,IG1,IL+1)=GAR(IGAR)*SIGS(IG1,IL+1)
   45       CONTINUE
   46       CONTINUE
         ELSE IF(.NOT.LPCAT) THEN
            DO 50 IG1=1,NGRO
            SCAT(IG1,IG1,IL+1)=SIGS(IG1,IL+1)
   50       CONTINUE
         ENDIF
      ELSE IF(LPCAT) THEN
         DO 70 IG1=1,NGRO
         TMP=0.0D0
         DO 60 IG2=1,NGRO
         TMP=TMP+SCAT(IG2,IG1,IL+1)
   60    CONTINUE
         SIGS(IG1,IL+1)=REAL(TMP)
   70    CONTINUE
         LSCAT(IL+1)=.TRUE.
      ENDIF
   80 CONTINUE
      LSCAT(1)=.TRUE.
      DO 90 IED=1,NED
      CALL XDRSET(SADD(1,IED),NGRO,0.0)
      CALL LCMLEN(IPDRL,HVECT(IED),LENGT,ITYLCM)
      LADD(IED)=(LENGT.GT.0)
      IF(LADD(IED)) CALL LCMGET(IPDRL,HVECT(IED),SADD(1,IED))
   90 CONTINUE
      CALL LCMLEN(IPDRL,'NGOLD',LENGT,ITYLCM)
      LGOLD=(LENGT.GT.0)
      IF(LGOLD) THEN
         CALL XDRSET(GOLD,NGRO,0.0)
         CALL LCMGET(IPDRL,'NGOLD',GOLD)
      ELSE
         CALL XDRSET(GOLD,NGRO,1.0)
      ENDIF
      IF(LBIN.GT.0) THEN
         CALL LCMGET(IPDRL,'BIN-'//HTOTAL,BIN(1,1))
         CALL LCMGET(IPDRL,'BIN-SIGS00',BIN(1,2))
         CALL LCMLEN(IPDRL,'BIN-NUSIGF',LENGF,ITYLCM)
         IF(LENGF.GT.0) THEN
            CALL LCMGET(IPDRL,'BIN-NUSIGF',BIN(1,3))
         ELSE
            CALL XDRSET(BIN(1,3),LBIN,0.0)
         ENDIF
         IGF0=0
         DO 190 IG=1,NGRO
         IF(NFS(IG).GT.0) THEN
*           BIN CROSS SECTION NORMALIZATION.
            SQ0=0.0D0
            SQ1=0.0D0
            SQ2=0.0D0
            SQ3=0.0D0
            DO 170 IGF=IGF0+1,IGF0+NFS(IG)
            DELTAU=LOG(BENER(IGF)/BENER(IGF+1))
            SQ0=SQ0+DELTAU
            SQ1=SQ1+BIN(IGF,1)*DELTAU
            SQ2=SQ2+(BIN(IGF,1)-BIN(IGF,2))*DELTAU
            SQ3=SQ3+BIN(IGF,3)*DELTAU
  170       CONTINUE
            FACT1=TOTAL(IG)*(SQ0/SQ1)
            FACT2=(TOTAL(IG)-SIGS(IG,1))*(SQ0/SQ2)
            DO 180 IGF=IGF0+1,IGF0+NFS(IG)
            BIN(IGF,2)=REAL(BIN(IGF,2)*FACT2+BIN(IGF,1)*(FACT1-FACT2))
            BIN(IGF,1)=REAL(BIN(IGF,1)*FACT1)
            IF((LENGF.GT.0).AND.(SQ3.NE.0.0)) THEN
               BIN(IGF,3)=REAL(BIN(IGF,3)*(SIGF(IG,0)*(SQ0/SQ3)))
            ENDIF
  180       CONTINUE
            IGF0=IGF0+NFS(IG)
         ENDIF
  190    CONTINUE
      ENDIF
      KTOTLR=0
      KSIGFR=0
      KPHIR=0
      KCHIR=0
      KADR=0
      CALL XDISET(KADDR,NED,0)
*----
*  PERFORM DILUTION INTERPOLATION.
*----
      CALL LCMLEN(IPDRL,'DILUTION',NDIL,ITYLCM)
      IF(NDIL.GT.0) THEN
         IF(NDIL+1.GT.MAXDIL) CALL XABORT('LIBDRB: INVALID MAXDIL.')
         CALL LCMGET(IPDRL,'DILUTION',DILUT)
         IF(DILUT(NDIL).GE.1.0E10) CALL XABORT('LIBDRB: INVALID DILUTI'
     1   //'ON VALUE.')
*----
*  FIND MAX LENGTH OF VECTORS ON SUBMAT
*  KTOTLR,KSIGFR,KCHIR,KPHIR AND KADDR
*  GIVES LENGTH OF SELF SHIELDING VECTOR
*  FOR TOTAL, SIGF, CHI, PHI AND ADD XS
*----
         CALL XDLSET(LSDIL,NL,.FALSE.)
         CALL XDLSET(LPDIL,NL,.FALSE.)
         DO 240 IDIL=1,NDIL
         WRITE (CD,'(I4.4)') IDIL
         CALL LCMSIX(IPDRL,'SUBMAT'//CD,1)
         CALL XDRSET(SIGT(1,IDIL),NGRO,0.0)
         CALL LCMGET(IPDRL,HTOTAL,SIGT(1,IDIL))
         DO 220 IL=0,NL-1
            WRITE (CM,'(I2.2)') IL
            CALL LCMLEN(IPDRL,'SCAT'//CM,LENGT,ITYLCM)
            IF(.NOT.LSDIL(IL+1))
     >        LSDIL(IL+1)=(LENGT.GT.0).AND.LSCAT(IL+1)
            CALL LCMLEN(IPDRL,'SIGS'//CM,LENGT,ITYLCM)
            IF(.NOT.LPDIL(IL+1))
     >        LPDIL(IL+1)=(LENGT.GT.0).AND.LSCAT(IL+1)
  220    CONTINUE
         CALL LCMLEN(IPDRL,HTOTAL,LENGT,ITYLCM)
         KTOTLR=MAX(KTOTLR,LENGT)
         CALL LCMLEN(IPDRL,'NUSIGF',LENGT,ITYLCM)
         KSIGFR=MAX(KSIGFR,LENGT)
         IF(NBESP.EQ.0) THEN
           CALL LCMLEN(IPDRL,'CHI',LENGT,ITYLCM)
           KCHIR=MAX(KCHIR,LENGT)
         ELSE
           DO 225 ISP=1,NBESP
           WRITE(HCHI,'(5HCHI--,I2.2)') ISP
           CALL LCMLEN(IPDRL,HCHI,LENGT,ITYLCM)
           KCHIR=MAX(KCHIR,LENGT)
  225      CONTINUE
         ENDIF
         CALL LCMLEN(IPDRL,'NWT0',LENGT,ITYLCM)
         KPHIR=MAX(KPHIR,LENGT)
         DO 230 IED=1,NED
            CALL LCMLEN(IPDRL,HVECT(IED),LENGT,ITYLCM)
            IF((LENGT.GT.0).AND.LADD(IED)) THEN
               KADDR(IED)=MAX(KADDR(IED),LENGT)
               KADR=MAX(KADDR(IED),KADR)
            ENDIF
  230    CONTINUE
         CALL LCMSIX(IPDRL,' ',2)
  240    CONTINUE
         NGRRE=MAX(KTOTLR,KSIGFR,KCHIR,KPHIR,KADR)
         IF(NGRRE.GT.NGRO) CALL XABORT('LIBDRB: TOO MANY GROUPS.')
*
         CALL XDRSET(TERP,MAXDIL*NGRO,0.0)
         DILUT(NDIL+1)=1.0E10
         DO 280 IG1=1,NGRRE
         LINF(IG1)=.FALSE.
         ZNPHI(IG1)=0.0
         DILX=MIN(SN(IG1),1.0E10)
         IF(DILX.LE.0.0) THEN
            WRITE (HSMG,930) IG1
            CALL XABORT(HSMG)
         ENDIF
         IFIRST=0
         DO 260 I=1,NDIL+1
         IF(ABS(DILX-DILUT(I)).LE.1.0E-5*ABS(DILX)) THEN
            TERP(I,IG1)=1.0
            GO TO 280
         ELSE IF(DILX.LT.DILUT(I)) THEN
            IFIRST=I-1
            GO TO 270
         ENDIF
  260    CONTINUE
*
  270    SQD=SQRT(DILX)
         IF((IFIRST-1.GE.1).AND.(IFIRST+2.LE.NDIL)) THEN
            SQ0=SQRT(DILUT(IFIRST-1))
            SQ1=SQRT(DILUT(IFIRST))
            SQ2=SQRT(DILUT(IFIRST+1))
            SQ3=SQRT(DILUT(IFIRST+2))
            TERP(IFIRST-1,IG1)=REAL((SQ1-SQD)*(SQ2-SQD)*(SQ3-SQD)/
     1      (SQ1-SQ0)/(SQ2-SQ0)/(SQ3-SQ0))
            TERP(IFIRST,IG1)=REAL((SQ0-SQD)*(SQ2-SQD)*(SQ3-SQD)/
     1      (SQ0-SQ1)/(SQ2-SQ1)/(SQ3-SQ1))
            TERP(IFIRST+1,IG1)=REAL((SQ0-SQD)*(SQ1-SQD)*(SQ3-SQD)/
     1      (SQ0-SQ2)/(SQ1-SQ2)/(SQ3-SQ2))
            TERP(IFIRST+2,IG1)=REAL((SQ0-SQD)*(SQ1-SQD)*(SQ2-SQD)/
     1      (SQ0-SQ3)/(SQ1-SQ3)/(SQ2-SQ3))
            TT=TERP(IFIRST-1,IG1)*SIGT(IG1,IFIRST-1)
     1        +TERP(IFIRST,IG1)*SIGT(IG1,IFIRST)
     2        +TERP(IFIRST+1,IG1)*SIGT(IG1,IFIRST+1)
     3        +TERP(IFIRST+2,IG1)*SIGT(IG1,IFIRST+2)
            YMIN=MIN(SIGT(IG1,IFIRST),SIGT(IG1,IFIRST+1))
            YMAX=MAX(SIGT(IG1,IFIRST),SIGT(IG1,IFIRST+1))
            IF((TT.GT.YMAX).OR.(TT.LT.YMIN)) THEN
               TERP(IFIRST-1,IG1)=0.0
               TERP(IFIRST,IG1)=REAL((SQ2-SQD)/(SQ2-SQ1))
               TERP(IFIRST+1,IG1)=REAL((SQ1-SQD)/(SQ1-SQ2))
               TERP(IFIRST+2,IG1)=0.0
            ENDIF
         ELSE IF((IFIRST.EQ.1).AND.(IFIRST+2.LE.NDIL)) THEN
            SQ1=SQRT(DILUT(1))
            SQ2=SQRT(DILUT(2))
            SQ3=SQRT(DILUT(3))
            TERP(1,IG1)=REAL((SQ2-SQD)*(SQ3-SQD)/(SQ2-SQ1)/(SQ3-SQ1))
            TERP(2,IG1)=REAL((SQ1-SQD)*(SQ3-SQD)/(SQ1-SQ2)/(SQ3-SQ2))
            TERP(3,IG1)=REAL((SQ1-SQD)*(SQ2-SQD)/(SQ1-SQ3)/(SQ2-SQ3))
            TT=TERP(1,IG1)*SIGT(IG1,1)+TERP(2,IG1)*SIGT(IG1,2)
     1        +TERP(3,IG1)*SIGT(IG1,3)
            YMIN=MIN(SIGT(IG1,1),SIGT(IG1,2))
            YMAX=MAX(SIGT(IG1,1),SIGT(IG1,2))
            IF((TT.GT.YMAX).OR.(TT.LT.YMIN)) THEN
               TERP(1,IG1)=REAL((SQ2-SQD)/(SQ2-SQ1))
               TERP(2,IG1)=REAL((SQ1-SQD)/(SQ1-SQ2))
               TERP(3,IG1)=0.0
            ENDIF
         ELSE IF((IFIRST-1.GE.1).AND.(IFIRST+1.EQ.NDIL)) THEN
            SQ0=SQRT(DILUT(NDIL-2))
            SQ1=SQRT(DILUT(NDIL-1))
            SQ2=SQRT(DILUT(NDIL))
          TERP(NDIL-2,IG1)=REAL((SQ1-SQD)*(SQ2-SQD)/(SQ1-SQ0)/(SQ2-SQ0))
          TERP(NDIL-1,IG1)=REAL((SQ0-SQD)*(SQ2-SQD)/(SQ0-SQ1)/(SQ2-SQ1))
            TERP(NDIL,IG1)=REAL((SQ0-SQD)*(SQ1-SQD)/(SQ0-SQ2)/(SQ1-SQ2))
            TT=TERP(NDIL-2,IG1)*SIGT(IG1,NDIL-2)
     1        +TERP(NDIL-1,IG1)*SIGT(IG1,NDIL-1)
     2        +TERP(NDIL,IG1)*SIGT(IG1,NDIL)
            YMIN=MIN(SIGT(IG1,NDIL-1),SIGT(IG1,NDIL))
            YMAX=MAX(SIGT(IG1,NDIL-1),SIGT(IG1,NDIL))
            IF((TT.GT.YMAX).OR.(TT.LT.YMIN)) THEN
               TERP(NDIL-2,IG1)=0.0
               TERP(NDIL-1,IG1)=REAL((SQ2-SQD)/(SQ2-SQ1))
               TERP(NDIL,IG1)=REAL((SQ1-SQD)/(SQ1-SQ2))
            ENDIF
         ELSE IF((IFIRST.EQ.0).OR.((IFIRST.EQ.1).AND.(NDIL.EQ.2))) THEN
            SQ0=SQRT(DILUT(1))
            SQ1=SQRT(DILUT(2))
            TERP(1,IG1)=REAL((SQ1-SQD)/(SQ1-SQ0))
            TERP(2,IG1)=REAL((SQ0-SQD)/(SQ0-SQ1))
         ELSE IF(IFIRST.EQ.NDIL) THEN
            LINF(IG1)=.TRUE.
            TERP(NDIL,IG1)=DILUT(NDIL)/DILX
         ELSE
            CALL XABORT('LIBDRB: FAILURE OF DILUTION INTERPOLATION.')
         ENDIF
  280    CONTINUE
*
         NGRODP=NGRO+1
         NGROIN=0
         DO 330 IDIL=1,NDIL
           NCORF=0
           DO 290 IG1=NGRO,1,-1
             IF(TERP(IDIL,IG1).NE.0.0) THEN
               NCORF=IG1
               GO TO 300
             ENDIF
  290      CONTINUE
  300      NGROIN=MAX(NCORF,NGROIN)
           NCORD=NGRO+1
           DO 310 IG1=1,NGROIN
             IF(TERP(IDIL,IG1).NE.0.0) THEN
               NCORD=IG1
               GO TO 320
             ENDIF
  310      CONTINUE
  320      NGRODP=MIN(NCORD,NGRODP)
  330    CONTINUE
         DO 345 IDIL=1,NDIL
         DO 340 IG1=1,NGRO
         IF(SIGT(IG1,IDIL).NE.0.0) THEN
            NGF=MIN(NGF,IG1-1)
            NGFR=MAX(NGFR,IG1)
         ENDIF
  340    CONTINUE
  345    CONTINUE
         IF(NGROIN.EQ.0.OR.NGRRE.EQ.0) THEN
            CALL XDRSET(ZNPHI,NGRO,1.0)
            GO TO 680
         ENDIF
         KTOTLR=MIN(KTOTLR,NGROIN)
         KSIGFR=MIN(KSIGFR,NGROIN)
         KCHIR=MIN(KCHIR,NGROIN)
         KPHIR=MIN(KPHIR,NGROIN)
         DO 360 IED=1,NED
           KADDR(IED)=MIN(KADDR(IED),NGROIN)
  360    CONTINUE
*----
*  VARIOUS DIMENSION OF VECTORS ARE SET
*  LOOP OVER DILUTION AND SELF-SHIELD XS
*  FROM NGRODP TO NGROIN (THESE CORRESPOND
*  TO CASES WHERE DIL<1.0E10 FOR AT LEAST ONE GROUP
*  HERE ONE ASSUMES THAT TOTAL XS ALWAYS SELF-SHIELDED
*----
         DO 490 IDIL=1,NDIL
           DO 370 IG1=1,NGRO
             IF(TERP(IDIL,IG1).NE.0.0) GO TO 380
  370      CONTINUE
           GO TO 490
  380    WRITE (CD,'(I4.4)') IDIL
         CALL LCMSIX(IPDRL,'SUBMAT'//CD,1)
         DO 390 IG1=NGRODP,NGROIN
         TOTAL(IG1)=TOTAL(IG1)+TERP(IDIL,IG1)*SIGT(IG1,IDIL)
  390    CONTINUE
         IF(KSIGFR.GT.0) THEN
            DO 401 IDEL=0,NDEL
            IF(IDEL.EQ.0) THEN
               HNUSIG='NUSIGF'
            ELSE
               WRITE(HNUSIG,'(6HNUSIGF,I2.2)') IDEL
            ENDIF
            CALL XDRSET(GAR,KSIGFR,0.0)
            CALL LCMGET(IPDRL,HNUSIG,GAR)
            DO 400 IG1=NGRODP,KSIGFR
            SIGF(IG1,IDEL)=SIGF(IG1,IDEL)+TERP(IDIL,IG1)*GAR(IG1)
  400       CONTINUE
  401       CONTINUE
         ENDIF
         IF(KCHIR.GT.0) THEN
            DO 410 IDEL=0,NDEL
            IF(IDEL.EQ.0) THEN
               IF(NBESP.GT.0) GO TO 410
               HCHI='CHI'
            ELSE
               WRITE(HCHI,'(3HCHI,I2.2)') IDEL
            ENDIF
            CALL XDRSET(GAR,KCHIR,0.0)
            CALL LCMGET(IPDRL,HCHI,GAR)
            DO 405 IG1=NGRODP,KCHIR
            CHI(IG1,IDEL)=CHI(IG1,IDEL)+TERP(IDIL,IG1)*GAR(IG1)
  405       CONTINUE
  410       CONTINUE
            DO 416 ISP=1,NBESP
            WRITE(HCHI,'(5HCHI--,I2.2)') ISP
            CALL XDRSET(GAR,KCHIR,0.0)
            CALL LCMLEN(IPDRL,HCHI,ILONG,ITYLCM)
            IF(ILONG.GT.0) THEN
              CALL LCMGET(IPDRL,HCHI,GAR)
              DO 415 IG1=NGRODP,KCHIR
              CHI4G(IG1,ISP)=CHI4G(IG1,ISP)+TERP(IDIL,IG1)*GAR(IG1)
  415         CONTINUE
            ENDIF
  416       CONTINUE
         ENDIF
         DO 450 IL=0,NL-1
         WRITE (CM,'(I2.2)') IL
         IF(LSDIL(IL+1)) THEN
            CALL LCMGET(IPDRL,'NJJS'//CM,NJJ)
            CALL LCMGET(IPDRL,'IJJS'//CM,IJJ)
            LENGT=0
            DO 420 I=1,NGRO
            LENGT=LENGT+NJJ(I)
  420       CONTINUE
            CALL XDRSET(GAR,LENGT,0.0)
            CALL LCMGET(IPDRL,'SCAT'//CM,GAR)
            IGAR=0
            DO 435 IG2=1,NGRO
            DO 430 IG1=IJJ(IG2),IJJ(IG2)-NJJ(IG2)+1,-1
            IGAR=IGAR+1
            SCAT(IG2,IG1,IL+1)=SCAT(IG2,IG1,IL+1)+TERP(IDIL,IG1)
     1      *GAR(IGAR)
  430       CONTINUE
  435       CONTINUE
         ENDIF
         IF(LPDIL(IL+1)) THEN
            CALL XDRSET(GAR,NGRO,0.0)
            CALL LCMGET(IPDRL,'SIGS'//CM,GAR)
            DO 440 IG1=NGRODP,NGROIN
            SIGS(IG1,IL+1)=SIGS(IG1,IL+1)+TERP(IDIL,IG1)*GAR(IG1)
  440       CONTINUE
         ENDIF
  450    CONTINUE
         IF(KPHIR.GT.0) THEN
            CALL XDRSET(GAR,KPHIR,0.0)
            CALL LCMGET(IPDRL,'NWT0',GAR)
            DO 460 IG1=NGRODP,KPHIR
            IF(.NOT.LINF(IG1)) THEN
               ZNPHI(IG1)=ZNPHI(IG1)+TERP(IDIL,IG1)*GAR(IG1)*
     1         DILUT(IDIL)
            ELSE
               ZNPHI(IG1)=GAR(IG1)*DILUT(IDIL)
            ENDIF
  460       CONTINUE
         ENDIF
         DO 480 IED=1,NED
         IF(KADDR(IED).GT.0) THEN
            CALL XDRSET(GAR,KADDR(IED),0.0)
            CALL LCMGET(IPDRL,HVECT(IED),GAR)
            DO 470 IG1=NGRODP,KADDR(IED)
            SADD(IG1,IED)=SADD(IG1,IED)+TERP(IDIL,IG1)*GAR(IG1)
  470       CONTINUE
         ENDIF
  480    CONTINUE
         CALL LCMSIX(IPDRL,' ',2)
  490    CONTINUE
*----
*  COMPUTE MISSING SCATTERING INFORMATION.
*----
         DO 550 IL=0,NL-1
         IF(LPDIL(IL+1).AND.(.NOT.LSDIL(IL+1))) THEN
            WRITE (CM,'(I2.2)') IL
            CALL LCMGET(IPDRL,'NJJS'//CM,NJJ)
            CALL LCMGET(IPDRL,'IJJS'//CM,IJJ)
            LENGT=0
            DO 500 I=1,NGRO
            LENGT=LENGT+NJJ(I)
  500       CONTINUE
            CALL XDRSET(GAR,LENGT,0.0)
            CALL LCMGET(IPDRL,'PCAT'//CM,GAR)
            IGAR=0
            DO 525 IG2=1,NGRO
            DO 510 IG1=1,NGRO
            SCAT(IG2,IG1,IL+1)=0.0
  510       CONTINUE
            DO 520 IG1=IJJ(IG2),IJJ(IG2)-NJJ(IG2)+1,-1
            IGAR=IGAR+1
            SCAT(IG2,IG1,IL+1)=GAR(IGAR)*SIGS(IG1,IL+1)
  520       CONTINUE
  525       CONTINUE
         ELSE IF((.NOT.LPDIL(IL+1)).AND.LSDIL(IL+1)) THEN
            DO 540 IG1=1,NGRO
            TMP=0.0D0
            DO 530 IG2=1,NGRO
            TMP=TMP+SCAT(IG2,IG1,IL+1)
  530       CONTINUE
            IF(IL.EQ.0) THEN
               SIGS(IG1,1)=MIN(REAL(TMP),TOTAL(IG1))
            ELSE
               SIGS(IG1,IL+1)=REAL(TMP)
            ENDIF
  540       CONTINUE
         ENDIF
  550    CONTINUE
*----
*  COMPUTE CONDENSED FINE STRUCTURE FUNCTION.
*----
         DO 580 IG1=1,NGROIN
         IF((.NOT.LSDIL(1)).AND.(.NOT.LPDIL(1))) THEN
*           SCATTERING CROSS SECTIONS ARE NOT SELF-SHIELDED.
            TMP=-TOTAL(IG1)
            DO 560 IG2=1,IG1-1
            TMP=TMP+SCAT(IG1,IG2,1)*ZNPHI(IG2)*DELTA(IG2)/DELTA(IG1)
  560       CONTINUE
            ZNGAR=(TMP+SCAT(IG1,IG1,1))*SB(IG1)/
     1                               (SB(IG1)-SCAT(IG1,IG1,1))
         ELSE
*           SCATTERING CROSS SECTIONS ARE SELF-SHIELDED.
            ZNGAR=-TOTAL(IG1)
            DO 570 IG2=1,IG1
            ZNGAR=ZNGAR+SCAT(IG1,IG2,1)*DELTA(IG2)/DELTA(IG1)
  570       CONTINUE
         ENDIF
         IF(IG1.LT.NGRODP) ZNGAR=0.0
         IF(KPHIR.EQ.0) THEN
*           USE A CALCULATED VALUE.
            ZNPHI(IG1)=REAL(1.0+ZNGAR/SB(IG1))
         ELSE IF(LINF(IG1)) THEN
*           USE AN INTERPOLATED VALUE NEAR INFINITE DILUTION.
            AUX=(DILUT(NDIL)/SB(IG1))**2
            ZNPHI(IG1)=REAL(AUX*ZNPHI(IG1)+(1.0-AUX)*ZNGAR)
            ZNPHI(IG1)=1.0+ZNPHI(IG1)/SB(IG1)
         ELSE
*           USE AN INTERPOLATED VALUE.
            ZNPHI(IG1)=1.0+ZNPHI(IG1)/SB(IG1)
         ENDIF
         IF((ZNPHI(IG1).LE.0.0).OR.(ZNPHI(IG1).GT.3.0)) THEN
            WRITE (HSMG,960) ZNPHI(IG1),IG1,SB(IG1),SN(IG1),KPHIR
            CALL XABORT(HSMG)
         ELSE IF(ZNPHI(IG1).GT.1.2) THEN
            WRITE (HSMG,960) ZNPHI(IG1),IG1,SB(IG1),SN(IG1),KPHIR
            WRITE(6,'(1X,A)') HSMG
         ENDIF
  580    CONTINUE
         DO 590 IG1=NGROIN+1,NGRO
         ZNPHI(IG1)=1.0
  590    CONTINUE
*----
*  DIVIDE EFFECTIVE REACTION RATES BY ZNPHI FOR SELF-SHIELDED
*  REACTION RATES
*----
         DO 620 IL=0,NL-1
         IF(LSCAT(IL+1).AND.(LSDIL(IL+1).OR.LPDIL(IL+1))) THEN
            DO 610 IG1=NGRODP,NGROIN
              SIGS(IG1,IL+1)=SIGS(IG1,IL+1)/ZNPHI(IG1)
              DO 600 IG2=1,NGRO
                SCAT(IG2,IG1,IL+1)=SCAT(IG2,IG1,IL+1)/ZNPHI(IG1)
  600         CONTINUE
  610       CONTINUE
         ENDIF
  620    CONTINUE
         DO 630 IG1=NGRODP,NGROIN
           TOTAL(IG1)=TOTAL(IG1)/ZNPHI(IG1)
  630    CONTINUE
         IF(KSIGFR.GT.0) THEN
           DO 645 IDEL=0,NDEL
           DO 640 IG1=NGRODP,NGROIN
             SIGF(IG1,IDEL)=SIGF(IG1,IDEL)/ZNPHI(IG1)
  640      CONTINUE
  645      CONTINUE
         ENDIF
         DO 660 IED=1,NED
           IF(KADDR(IED).GT.0) THEN
             DO 650 IG1=NGRODP,NGROIN
               SADD(IG1,IED)=SADD(IG1,IED)/ZNPHI(IG1)
  650        CONTINUE
           ENDIF
  660    CONTINUE
         IF(IMPX.GT.4) THEN
            WRITE(IOUT,940)
            DO 670 IG1=1,NGRO
            WRITE (IOUT,950) IG1,SN(IG1),SB(IG1),ZNPHI(IG1),TOTAL(IG1),
     1      SIGS(IG1,1),SIGF(IG1,0),GOLD(IG1)
  670       CONTINUE
            WRITE (IOUT,'(/)')
         ENDIF
      ELSE
         CALL XDRSET(ZNPHI,NGRO,1.0)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
  680 DEALLOCATE(LINF,LPDIL,LSDIL)
      DEALLOCATE(SIGT,TERP,GAR)
      DEALLOCATE(KADDR,IJJ,NJJ)
      RETURN
*
  930 FORMAT(42HLIBDRB: NEGATIVE OR ZERO DILUTION IN GROUP,I4,1H.)
  940 FORMAT(/5X,'GROUP',10X,'DILUT',13X,'SB',11X,'NWT0',10X,'NTOT0',
     1 11X,'SIGS',9X,'NUSIGF',10X,'NGOLD')
  950 FORMAT(5X,I5,1P,8E15.5)
  960 FORMAT(32HLIBDRB: INVALID VALUE OF ZNPHI (,1P,E11.3,
     1 10H) IN GROUP,I4,11H. DILUTION=,E11.3,2H (,E11.3,
     2 9H). KPHIR=,I4,1H.)
      END
