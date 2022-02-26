*DECK B1HXS1
      SUBROUTINE B1HXS1(IPMACR,NGRO,NBM,IAN,NFISSI,IJJ0,IJJ1,NJJ0,NJJ1,
     1 IDEL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Homogenization of the lattice cell nuclear properties before a B-n
* calculation.
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
* IPMACR  pointer to the macrolib LCM object (L_MACROLIB signature).
* NGRO    number of groups.
* NBM     number of mixtures.
* IAN     type of homogenization:
*         =-1: transport corrected P0; =0: P0; =1: P1.
* NFISSI  maximum number of fission spectrum assigned to a mixture.
* NUNKNO  number of flux/current unknowns.
* IPAS    number of volumes.
* VOL     volumes.
* MAT     mixture number of each volume.
* KEYFLX  position of each flux in the unknown vector.
* FLUX    direct unknown vector.
* INORM   type of leakage model:
*         =1: Diffon; =2: Ecco; =3: Tibere.
*
*Parameters: output
* IJJ0    most thermal group in band for P0 scattering.
* NJJ0    number of groups in band for P0 scattering.
* IJJ1    most thermal group in band for P1 scattering.
* NJJ1    number of groups in band for P1 scattering.
* IDEL    dimension of matrices SCAT0 and SCAT1.
* FLXIN   integrated fluxes.
* SA      absorption macroscopic cross sections.
* ST      total macroscopic cross sections.
* SFNU    nu * macroscopic fission cross-sections.
* XHI     fission spectrum.
* SCAT0   packed diffusion P0 macroscopic cross sections.
* SCAT1   packed diffusion P1 macroscopic cross sections.
* NGROIN  number of groups without up-scattering.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMACR
      INTEGER NGRO,NBM,IAN,NFISSI,IJJ0(NGRO),IJJ1(NGRO),NJJ0(NGRO),
     1 NJJ1(NGRO),IDEL(2)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      TYPE(C_PTR) JPMACR,KPMACR
      LOGICAL LOGIC
      CHARACTER CM*2
      INTEGER IDATA(NSTATE)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,NJJ
*----
*  SCRATCH STORAGE ALLOCATION
*   IJJ     last scattering group (IJJ(0) = 0).
*   NJJ     number of scattering group (NJJ(0)=-NGROUP).
*----
      ALLOCATE(IJJ(0:NBM),NJJ(0:NBM))
*
      CALL LCMGET(IPMACR,'STATE-VECTOR',IDATA)
      LOGIC=(NGRO.EQ.IDATA(1)).AND.(NBM.EQ.IDATA(2)).AND.(NFISSI.EQ.
     1 IDATA(4)).AND.(IDATA(3).GE.1)
      IF(.NOT.LOGIC) CALL XABORT('B1HXS1: INCONSISTENT LCM FILE.')
      IANN=IAN
      IF(IAN.LT.0) IANN=-(IAN+1)
      IDEL(1)=0
      IDEL(2)=0
      JPMACR=LCMGID(IPMACR,'GROUP')
      DO 30 LLL=1,NGRO
      KPMACR=LCMGIL(JPMACR,LLL)
      DO 20 M=0,IANN
      WRITE (CM,'(I2.2)') M
      CALL LCMGET(KPMACR,'NJJS'//CM,NJJ(1))
      CALL LCMGET(KPMACR,'IJJS'//CM,IJJ(1))
      IMAX=1
      IMIN=NGRO
      DO 10 I=1,NBM
      IMAX=MAX(IJJ(I),IMAX)
      IMIN=MIN(IJJ(I)-NJJ(I)+1,IMIN)
   10 CONTINUE
      IF(M.EQ.0) THEN
         IJJ0(LLL)=IMAX
         NJJ0(LLL)=IMAX-IMIN+1
      ELSE IF(M.EQ.1) THEN
         IJJ1(LLL)=IMAX
         NJJ1(LLL)=IMAX-IMIN+1
      ENDIF
      IDEL(M+1)=IDEL(M+1)+IMAX-IMIN+1
   20 CONTINUE
   30 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(NJJ,IJJ)
      RETURN
      END
*
      SUBROUTINE B1HXS2(NUNKNO,IPMACR,IPAS,NGRO,NBM,IAN,NFISSI,VOL,MAT,
     1 KEYFLX,FLUX,IJJ0,IJJ1,NJJ0,NJJ1,IDEL,FLXIN,SA,ST,SFNU,XHI,SCAT0,
     2 SCAT1,NGROIN,INORM)
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMACR
      INTEGER NUNKNO,IPAS,NGRO,NBM,IAN,NFISSI,MAT(IPAS),KEYFLX(IPAS),
     1 IJJ0(NGRO),IJJ1(NGRO),NJJ0(NGRO),NJJ1(NGRO),IDEL(2),NGROIN,INORM
      REAL VOL(IPAS),FLUX(NUNKNO,NGRO),SA(NGRO),ST(NGRO),SFNU(NGRO),
     1 XHI(NGRO),SCAT0(IDEL(1)),SCAT1(IDEL(2))
      DOUBLE PRECISION FLXIN(NGRO)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPMACR,KPMACR
      LOGICAL LOGIC
      CHARACTER CM*2
      DOUBLE PRECISION SUM,A11,A13
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,NJJ,IPOS
      REAL, ALLOCATABLE, DIMENSION(:) :: XSCAT,GAR,GARFI
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: GAF
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: A14
*----
*  SCRATCH STORAGE ALLOCATION
*   XSCAT   scattering vector (XSCAT(0)=0.0).
*   IJJ     last scattering group (IJJ(0) = 0).
*   NJJ     number of scattering group (NJJ(0)=-NGROUP).
*   IPOS    position self scattering in XSCAT (IPOS(0)=NGROUP+1).
*----
      ALLOCATE(IJJ(0:NBM),NJJ(0:NBM),IPOS(0:NBM))
      ALLOCATE(XSCAT(0:NBM*NGRO),GAR(0:NBM),GARFI(0:NBM*NFISSI))
      ALLOCATE(A14(NFISSI,0:NBM),GAF(NGRO))
*
      IANN=IAN
      IF(IAN.LT.0) IANN=-(IAN+1)
      NGROIN=0
      SUM=0.0D0
      A13=0.0D0
      DO 45 NF=1,NFISSI
      DO 40 IBM=1,NBM
      A14(NF,IBM)=0.0D0
   40 CONTINUE
   45 CONTINUE
      JPMACR=LCMGID(IPMACR,'GROUP')
      IF(NFISSI.GT.0) THEN
         DO 62 LLL=1,NGRO
         KPMACR=LCMGIL(JPMACR,LLL)
         A13=0.0D0
         CALL LCMGET(KPMACR,'NUSIGF',GARFI(1))
         DO 61 NF=1,NFISSI
         DO 60 I=1,IPAS
         IBM=MAT(I)
         IF(IBM.GT.0) A14(NF,IBM)=A14(NF,IBM)+FLUX(KEYFLX(I),LLL)*
     1   VOL(I)*GARFI((NF-1)*NBM+IBM)
   60    CONTINUE
   61    CONTINUE
   62    CONTINUE
         DO 75 NF=1,NFISSI
         DO 70 IBM=1,NBM
         A13=A13+A14(NF,IBM)
   70    CONTINUE
   75    CONTINUE
      ENDIF
*
      IF(INORM.EQ.1) THEN
        DO 85 LLL=1,NGRO
        A11=0.0D0
        DO 80 I=1,IPAS
        A11=A11+FLUX(KEYFLX(I),LLL)*VOL(I)
   80   CONTINUE
        FLXIN(LLL)=A11
   85   CONTINUE
        IDEL(1)=0
        IDEL(2)=0
      ENDIF
*
      DO 200 LLL=1,NGRO
      KPMACR=LCMGIL(JPMACR,LLL)
      IF(NFISSI.GT.0) THEN
         A11=0.0D0
         CALL LCMGET(KPMACR,'NUSIGF',GARFI(1))
         DO 95 NF=1,NFISSI
         DO 90 I=1,IPAS
         IBM=MAT(I)
         IF(IBM.GT.0) A11=A11+FLUX(KEYFLX(I),LLL)*VOL(I)*
     1   GARFI((NF-1)*NBM+IBM)
   90    CONTINUE
   95    CONTINUE
         SFNU(LLL)=REAL(A11/FLXIN(LLL))
      ELSE
         SFNU(LLL)=0.0
      ENDIF
*
      GAR(0)=0.0
      IF(INORM.EQ.1) THEN
        CALL LCMGET(KPMACR,'NTOT0',GAR(1))
        A11=0.0D0
        DO 100 I=1,IPAS
        A11=A11+FLUX(KEYFLX(I),LLL)*VOL(I)*GAR(MAT(I))
  100   CONTINUE
        ST(LLL)=REAL(A11/FLXIN(LLL))
      ELSE
        A11=ST(LLL)*FLXIN(LLL)
      ENDIF
*
      CALL LCMGET(KPMACR,'SIGS00',GAR(1))
      DO 110 I=1,IPAS
      A11=A11-FLUX(KEYFLX(I),LLL)*VOL(I)*GAR(MAT(I))
  110 CONTINUE
      SA(LLL)=REAL(A11/FLXIN(LLL))
*
      IF(NFISSI.GT.0) THEN
         A11=0.0D0
         CALL LCMGET(KPMACR,'CHI',GARFI(1))
         DO 125 NF=1,NFISSI
         DO 120 IBM=1,NBM
         A11=A11+A14(NF,IBM)*GARFI((NF-1)*NBM+IBM)
  120    CONTINUE
  125    CONTINUE
         XHI(LLL)=REAL(A11/A13)
         SUM=SUM+XHI(LLL)
      ELSE
         XHI(LLL)=0.0
      ENDIF
      IF(INORM.EQ.1) THEN
*----
*  TRANSPORT CORRECTION
*----
        A11=0.0D0
        IF(IAN.EQ.-1) THEN
           GAR(0)=0.0
           CALL LCMGET(KPMACR,'SIGS01',GAR(1))
           DO 130 I=1,IPAS
           A11=A11+FLUX(KEYFLX(I),LLL)*VOL(I)*GAR(MAT(I))
  130      CONTINUE
           ST(LLL)=ST(LLL)-REAL(A11/FLXIN(LLL))
        ENDIF
*
        DO 190 M=0,IANN
        WRITE (CM,'(I2.2)') M
        CALL LCMGET(KPMACR,'NJJS'//CM,NJJ(1))
        CALL LCMGET(KPMACR,'IJJS'//CM,IJJ(1))
        CALL LCMGET(KPMACR,'IPOS'//CM,IPOS(1))
        CALL LCMGET(KPMACR,'SCAT'//CM,XSCAT(1))
        DO 140 IG=1,NGRO
        GAF(IG)=0.0D0
  140   CONTINUE
        DO 160 I=1,IPAS
        IBM=MAT(I)
        IF(IBM.EQ.0) GO TO 160
        DO 150 IG=IJJ(IBM)-NJJ(IBM)+1,IJJ(IBM)
        IGAR=IPOS(IBM)+IJJ(IBM)-IG
        GAF(IG)=GAF(IG)+FLUX(KEYFLX(I),IG)*VOL(I)*XSCAT(IGAR)
  150   CONTINUE
        IF(IAN.EQ.-1) THEN
           IGAR=IPOS(IBM)+IJJ(IBM)-LLL
           GAF(LLL)=GAF(LLL)-FLUX(KEYFLX(I),LLL)*VOL(I)*GAR(IBM)
        ENDIF
  160   CONTINUE
        IF(M.EQ.0) THEN
           DO 170 IG=IJJ0(LLL)-NJJ0(LLL)+1,IJJ0(LLL)
           IGAR=IDEL(1)+1+IJJ0(LLL)-IG
           SCAT0(IGAR)=REAL(GAF(IG)/FLXIN(IG))
  170      CONTINUE
           IDEL(1)=IDEL(1)+NJJ0(LLL)
        ELSE IF(M.EQ.1) THEN
           DO 180 IG=IJJ1(LLL)-NJJ1(LLL)+1,IJJ1(LLL)
           IGAR=IDEL(2)+1+IJJ1(LLL)-IG
           SCAT1(IGAR)=REAL(GAF(IG)/FLXIN(IG))
  180      CONTINUE
           IDEL(2)=IDEL(2)+NJJ1(LLL)
        ENDIF
  190   CONTINUE
        LOGIC=(IJJ0(LLL).LE.LLL).AND.(NGROIN.EQ.LLL-1)
        IF(IANN.GE.1) LOGIC=LOGIC.AND.(IJJ1(LLL).LE.LLL)
        IF(LOGIC) NGROIN=LLL
      ENDIF
  200 CONTINUE
      IF((ABS(1.0D0-SUM).GT.1.0D-3).AND.(NFISSI.GT.0)) THEN
         CALL XABORT('B1HXS2: INCONSISTENT FISSION SPECTRUM.')
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(GAF,A14)
      DEALLOCATE(GARFI,GAR,XSCAT)
      DEALLOCATE(IPOS,NJJ,IJJ)
      RETURN
      END
*
      SUBROUTINE B1HXS3(NUNKNO,IPMACR,IPAS,NGRO,NBM,IAN,VOL,MAT,
     1 KEYFLX,FLUX,IJJ0,IJJ1,NJJ0,NJJ1,IDEL,FLXIN,ST,SCAT0,SCAT1,
     2 NGROIN)
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMACR
      INTEGER NUNKNO,IPAS,NGRO,NBM,IAN,MAT(IPAS),KEYFLX(IPAS),
     1 IJJ0(NGRO),IJJ1(NGRO),NJJ0(NGRO),NJJ1(NGRO),IDEL(2),NGROIN
      REAL VOL(IPAS),FLUX(NUNKNO,NGRO),ST(NGRO),SCAT0(IDEL(1)),
     1 SCAT1(IDEL(2))
      DOUBLE PRECISION FLXIN(NGRO)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPMACR,KPMACR
      LOGICAL LOGIC
      CHARACTER CM*2
      DOUBLE PRECISION A11,A13
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,NJJ,IPOS
      REAL, ALLOCATABLE, DIMENSION(:) :: XSCAT,GAR
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: GAF,CUR
*----
*  SCRATCH STORAGE ALLOCATION
*   XSCAT   scattering vector (XSCAT(0)=0.0).
*   IJJ     last scattering group (IJJ(0) = 0).
*   NJJ     number of scattering group (NJJ(0)=-NGROUP).
*   IPOS    position self scattering in XSCAT (IPOS(0)=NGROUP+1).
*----
      ALLOCATE(IJJ(0:NBM),NJJ(0:NBM),IPOS(0:NBM))
      ALLOCATE(XSCAT(0:NBM*NGRO),GAR(0:NBM))
      ALLOCATE(GAF(NGRO),CUR(NGRO))
*
      IANN=IAN
      IF(IAN.LT.0) IANN=-(IAN+1)
      NGROIN=0
*----
*  FIND HOMOGENISED FLUX AND CURRENTS
*----
      DO 305 LLL=1,NGRO
      FLXIN(LLL)=0.0D0
      DO 300 I=1,IPAS
      FLXIN(LLL)=FLXIN(LLL)+FLUX(KEYFLX(I),LLL)*VOL(I)
  300 CONTINUE
  305 CONTINUE
      DO 320 LLL=1,NGRO
      CUR(LLL)=0.0D0
      A13=0.0D0
      DO 310 I=1,IPAS
      A13=A13+FLUX(NUNKNO/2+KEYFLX(I),LLL)*VOL(I)
  310 CONTINUE
      CUR(LLL)=CUR(LLL)+A13
  320 CONTINUE
*
      IDEL(1)=0
      IDEL(2)=0
      JPMACR=LCMGID(IPMACR,'GROUP')
      DO 410 LLL=1,NGRO
      KPMACR=LCMGIL(JPMACR,LLL)
      GAR(0)=0.0
      CALL LCMGET(KPMACR,'NTOT0',GAR(1))
      A11=0.0D0
      DO 330 I=1,IPAS
      A11=A11+FLUX(KEYFLX(I),LLL)*VOL(I)*GAR(MAT(I))
  330 CONTINUE
      ST(LLL)=REAL(A11/FLXIN(LLL))
      A11=ST(LLL)*CUR(LLL)
      DO 340 I=1,IPAS
      A11=A11-VOL(I)*GAR(MAT(I))*FLUX(NUNKNO/2+KEYFLX(I),LLL)
  340 CONTINUE
*
      DO 400 M=0,IANN
      WRITE (CM,'(I2.2)') M
      CALL LCMGET(KPMACR,'NJJS'//CM,NJJ(1))
      CALL LCMGET(KPMACR,'IJJS'//CM,IJJ(1))
      CALL LCMGET(KPMACR,'IPOS'//CM,IPOS(1))
      CALL LCMGET(KPMACR,'SCAT'//CM,XSCAT(1))
      DO 350 IG=1,NGRO
      GAF(IG)=0.0D0
  350 CONTINUE
      IF(M.EQ.0) THEN
         DO 365 I=1,IPAS
         IBM=MAT(I)
         IF(IBM.EQ.0) GO TO 365
         DO 360 IG=IJJ(IBM)-NJJ(IBM)+1,IJJ(IBM)
         IGAR=IPOS(IBM)+IJJ(IBM)-IG
         GAF(IG)=GAF(IG)+FLUX(KEYFLX(I),IG)*VOL(I)*XSCAT(IGAR)
  360    CONTINUE
  365    CONTINUE
         DO 370 IG=IJJ0(LLL)-NJJ0(LLL)+1,IJJ0(LLL)
         IGAR=IDEL(1)+1+IJJ0(LLL)-IG
         SCAT0(IGAR)=REAL(GAF(IG)/FLXIN(IG))
  370    CONTINUE
         IDEL(1)=IDEL(1)+NJJ0(LLL)
      ELSE IF(M.EQ.1) THEN
         DO 385 I=1,IPAS
         IBM=MAT(I)
         IF(IBM.EQ.0) GO TO 385
         DO 380 IG=IJJ(IBM)-NJJ(IBM)+1,IJJ(IBM)
         IGAR=IPOS(IBM)+IJJ(IBM)-IG
         GAF(IG)=GAF(IG)+VOL(I)*XSCAT(IGAR)*FLUX(NUNKNO/2+KEYFLX(I),IG)
  380    CONTINUE
  385    CONTINUE
         GAF(LLL)=GAF(LLL)+A11
         DO 390 IG=IJJ1(LLL)-NJJ1(LLL)+1,IJJ1(LLL)
         IGAR=IDEL(2)+1+IJJ1(LLL)-IG
         SCAT1(IGAR)=REAL(GAF(IG)/CUR(IG))
  390    CONTINUE
         IDEL(2)=IDEL(2)+NJJ1(LLL)
      ENDIF
  400 CONTINUE
      LOGIC=(IJJ0(LLL).LE.LLL).AND.(NGROIN.EQ.LLL-1)
      IF(IANN.GE.1) LOGIC=LOGIC.AND.(IJJ1(LLL).LE.LLL)
      IF(LOGIC) NGROIN=LLL
  410 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(CUR,GAF)
      DEALLOCATE(GAR,XSCAT)
      DEALLOCATE(IPOS,NJJ,IJJ)
      RETURN
      END
