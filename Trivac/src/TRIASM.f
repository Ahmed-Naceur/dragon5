*DECK TRIASM
      SUBROUTINE TRIASM(HNAMT,IPTRK,IPSYS,IMPX,MAXMIX,NEL,NALBP,IPR,
     1 MAT,VOL,GAMMA,SGD,XSGD)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembly of a single-group system matrix with leakage and removal
* cross sections.
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
* HNAMT   name of the matrix.
* IPTRK   L_TRACK pointer to the TRIVAC tracking information.
* IPSYS   L_SYSTEM pointer to system matrices.
* IMPX    print parameter (equal to zero for no print).
* MAXMIX  first dimension for matrices SGD and XSGD.
* NEL     total number of finite elements.
* NALBP   number of physical albedos.
* IPR     type of assembly:
*         =0: calculation of the system matrices;
*         =1: calculation of the derivative of these matrices;
*         =2: calculation of the first variation of these matrices;
*         =3: identical to IPR=2, but these variation are added to
*         unperturbed system matrices.
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* GAMMA   physical albedo functions.
* SGD     nuclear properties per material mixture.
* XSGD    first variations or derivatives of nuclear properties:
*         if IPR.ge.1, XSGD contain first variations or derivatives
*         of nuclear properties in each material mixture;
*         if IPR=0, XSGD should be equivalenced with SGD. This is
*         obtained using 'CALL TRIASM(...,SGD,SGD)'.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPSYS
      CHARACTER HNAMT*10
      INTEGER IMPX,MAXMIX,NEL,IPR,MAT(NEL)
      REAL VOL(NEL),GAMMA(NALBP),SGD(MAXMIX,4),XSGD(MAXMIX,4)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      LOGICAL CYLIND,CHEX,DIAG,LSGD,LOGY,LOGZ
      CHARACTER TEXT10*10
      INTEGER NCODE(6),ICODE(6),ISTATE(NSTATE)
      REAL ZCODE(6),ZALB(6)
      INTEGER, DIMENSION(:), ALLOCATABLE :: KN,IQFR,MUW,MUZ,MATN,IPERT
      INTEGER, DIMENSION(:), ALLOCATABLE, TARGET :: MUY
      INTEGER, DIMENSION(:), POINTER :: MUX
      REAL, DIMENSION(:), ALLOCATABLE :: VOL2,QFR,XX,YY,ZZ,DD,T,TS,FRZ,
     1 DIF
      REAL, DIMENSION(:,:), ALLOCATABLE :: R,RS,Q,QS,V,RH,QH,RT,QT,DSGD
      REAL, DIMENSION(:), ALLOCATABLE :: RR0,XR0,ANG
      INTEGER, DIMENSION(:), POINTER :: IPW,IPX,IPY,IPZ
      INTEGER, DIMENSION(:), POINTER :: IPBW,IPBX,IPBY,IPBZ
      REAL, DIMENSION(:), POINTER :: TF,WA,AW,XA,AX,YA,AY,ZA,AZ
      REAL, DIMENSION(:), POINTER :: BW,BX,BY,BZ
      TYPE(C_PTR) IPW_PTR,IPX_PTR,IPY_PTR,IPZ_PTR
      TYPE(C_PTR) IPBW_PTR,IPBX_PTR,IPBY_PTR,IPBZ_PTR
      TYPE(C_PTR) TF_PTR,WA_PTR,AW_PTR,XA_PTR,AX_PTR,YA_PTR,AY_PTR,
     1 ZA_PTR,AZ_PTR
      TYPE(C_PTR) BW_PTR,BX_PTR,BY_PTR,BZ_PTR
*----
*  RECOVER TRIVAC SPECIFIC TRACKING INFORMATION
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      ITYPE=ISTATE(6)
      CYLIND=(ITYPE.EQ.3).OR.(ITYPE.EQ.6)
      CHEX=(ITYPE.EQ.8).OR.(ITYPE.EQ.9)
      IDIM=1
      IF((ITYPE.EQ.5).OR.(ITYPE.EQ.6).OR.(ITYPE.EQ.8)) IDIM=2
      IF((ITYPE.EQ.7).OR.(ITYPE.EQ.9)) IDIM=3
      IHEX=ISTATE(7)
      DIAG=(ISTATE(8).EQ.1)
      IELEM=ISTATE(9)
      ICOL=ISTATE(10)
      LL4=ISTATE(11)
      ICHX=ISTATE(12)
      ISPLH=ISTATE(13)
      LX=ISTATE(14)
      LY=ISTATE(15)
      LZ=ISTATE(16)
      ISEG=ISTATE(17)
      IMPV=ISTATE(18)
      NR0=ISTATE(24)
      LL4F=ISTATE(25)
      IF(ICHX.EQ.2) THEN
         ITY=3
         LL4W=ISTATE(26)
         LL4X=ISTATE(27)
         LL4Y=ISTATE(28)
         LL4Z=ISTATE(29)
         LOGY=LL4Y.GT.0
         LOGZ=LL4Z.GT.0
      ELSE
         ITY=2
         LL4W=LL4
         LL4X=LL4
         LL4Y=LL4
         LL4Z=LL4
         LOGY=IDIM.GT.1
         LOGZ=IDIM.GT.2
      ENDIF
      CALL LCMLEN(IPTRK,'KN',MAXKN,ITYLCM)
      CALL LCMLEN(IPTRK,'QFR',MAXQF,ITYLCM)
      ALLOCATE(ZZ(LX*LY*LZ),KN(MAXKN),QFR(MAXQF),IQFR(MAXQF))
      CALL LCMGET(IPTRK,'ZZ',ZZ)
      CALL LCMGET(IPTRK,'KN',KN)
      CALL LCMGET(IPTRK,'QFR',QFR)
      CALL LCMGET(IPTRK,'IQFR',IQFR)
      IF(CHEX) THEN
         CALL LCMGET(IPTRK,'SIDE',SIDE)
         ALLOCATE(MUW(LL4W))
         CALL LCMGET(IPTRK,'MUW',MUW)
      ELSE
         ALLOCATE(XX(LX*LY*LZ),YY(LX*LY*LZ),DD(LX*LY*LZ))
         CALL LCMGET(IPTRK,'XX',XX)
         CALL LCMGET(IPTRK,'YY',YY)
         CALL LCMGET(IPTRK,'DD',DD)
      ENDIF
      IF(LOGY) THEN
         ALLOCATE(MUY(LL4Y))
         CALL LCMGET(IPTRK,'MUY',MUY)
      ENDIF
      IF(.NOT.DIAG) THEN
         ALLOCATE(MUX(LL4X))
         CALL LCMGET(IPTRK,'MUX',MUX)
      ELSE
         MUX=>MUY
      ENDIF
      IF(LOGZ) THEN
         ALLOCATE(MUZ(LL4Z))
         CALL LCMGET(IPTRK,'MUZ',MUZ)
      ENDIF
*----
*  RECOVER UNIT MATRICES
*----
      IF((ICHX.EQ.1).OR.(ICHX.EQ.2)) THEN
         CALL LCMSIX(IPTRK,'BIVCOL',1)
         CALL LCMLEN(IPTRK,'T',LC,ITYLCM)
         ALLOCATE(T(LC),TS(LC),R(LC,LC),RS(LC,LC),Q(LC,LC),QS(LC,LC),
     1   V(LC,LC-1),RH(6,6),QH(6,6),RT(3,3),QT(3,3))
         CALL LCMGET(IPTRK,'T',T)
         CALL LCMGET(IPTRK,'TS',TS)
         CALL LCMGET(IPTRK,'R',R)
         CALL LCMGET(IPTRK,'RS',RS)
         CALL LCMGET(IPTRK,'Q',Q)
         CALL LCMGET(IPTRK,'QS',QS)
         CALL LCMGET(IPTRK,'V',V)
         IF((IELEM.EQ.1).AND.(ICOL.LE.2)) THEN
            CALL LCMGET(IPTRK,'RH',RH)
            CALL LCMGET(IPTRK,'QH',QH)
            CALL LCMGET(IPTRK,'RT',RT)
            CALL LCMGET(IPTRK,'QT',QT)
         ENDIF
         CALL LCMSIX(IPTRK,' ',2)
      ENDIF
*
      TEXT10=HNAMT(:10)
      IF(IMPX.GT.0) WRITE(6,'(/36H TRIASM: ASSEMBLY OF SYMMETRIC MATRI,
     1 3HX '',A10,38H'' IN COMPRESSED DIAGONAL STORAGE MODE.)') TEXT10
      CALL KDRCPU(TK1)
*----
*  COMPUTE THE INVERSE CROSS SECTIONS FOR DUAL FINITE ELEMENT CASES
*----
      IF(ICHX.EQ.2) THEN
         ALLOCATE(DSGD(MAXMIX,4))
         IF(IPR.EQ.0) THEN
            DO 15 J=1,4
            DO 10 I=1,MAXMIX
            IF(SGD(I,J).NE.0.) DSGD(I,J)=1.0/SGD(I,J)
   10       CONTINUE
   15       CONTINUE
         ELSE IF(IPR.EQ.1) THEN
            DO 25 J=1,4
            DO 20 I=1,MAXMIX
            IF(SGD(I,J).NE.0.0) THEN
               DSGD(I,J)=-XSGD(I,J)/(SGD(I,J)**2)
            ENDIF
   20       CONTINUE
   25       CONTINUE
         ELSE
            DO 35 J=1,4
            DO 30 I=1,MAXMIX
            SIGMA=SGD(I,J)+XSGD(I,J)
            IF((SGD(I,J).NE.0.0).AND.(SIGMA.NE.0.0)) THEN
               DSGD(I,J)=1.0/SIGMA-1.0/SGD(I,J)
            ENDIF
   30       CONTINUE
   35       CONTINUE
         ENDIF
      ENDIF
*----
*  DETERMINATION OF THE PERTURBED ELEMENTS AND INCLUSION OF ELEMENTS
*  NEIGHBOUR TO PERTURBED ZONES IN MCFD CASES. NON-PERTURBED ELEMENTS
*  WILL HAVE VOL2(K)=0.0
*----
      ALLOCATE(VOL2(NEL))
      IF((IPR.EQ.0).OR.(NALBP.GT.0)) THEN
         DO 40 K=1,NEL
         VOL2(K)=VOL(K)
   40    CONTINUE
      ELSE
         CALL XDRSET(VOL2,NEL,0.0)
         IF(ICHX.EQ.3) THEN
*           MCFD CASE.
            NUM1=0
            DO 70 L=1,NEL
            IF(MAT(L).EQ.0) GO TO 70
            LSGD=.FALSE.
            DO 50 I=1,4
            LSGD=LSGD.OR.(XSGD(MAT(L),I).NE.0.0)
   50       CONTINUE
            IF(LSGD) THEN
               VOL2(L)=VOL(L)
               DO 60 I=1,6
               K=KN(NUM1+I)
               IF(K.GT.0) THEN
                  IF(K.GT.NEL) CALL XABORT('TRIASM: INVALID BOUNDARY E'
     1            //'LEMENT INDEX.')
                  VOL2(K)=VOL(K)
               ENDIF
   60          CONTINUE
            ENDIF
            NUM1=NUM1+6
   70       CONTINUE
         ELSE
            DO 90 L=1,NEL
            IF(MAT(L).EQ.0) GO TO 90
            LSGD=.FALSE.
            DO 80 I=1,4
            LSGD=LSGD.OR.(XSGD(MAT(L),I).NE.0.0)
   80       CONTINUE
            IF(LSGD) VOL2(L)=VOL(L)
   90       CONTINUE
         ENDIF
      ENDIF
*----
*  APPLY PHYSICAL ALBEDOS AND INTRODUCE THE CYLINDER BOUNDARY
*  APPROXIMATION IN CARTESIAN GEOMETRY
*----
      IF(NR0.GT.0) THEN
         IF(IPR.GT.0) CALL XABORT('TRIASM: PERTURBATION CALCULATION NO'
     1   //'T AVAILABLE WITH CYLINDRICAL CORRECTION.')
         ALLOCATE(RR0(NR0),XR0(NR0),ANG(NR0))
         CALL LCMGET(IPTRK,'RR0',RR0)
         CALL LCMGET(IPTRK,'XR0',XR0)
         CALL LCMGET(IPTRK,'ANG',ANG)
         CALL LCMGET(IPTRK,'NCODE',NCODE)
         CALL LCMGET(IPTRK,'ICODE',ICODE)
         CALL LCMGET(IPTRK,'ZCODE',ZCODE)
         DO IC=1,6
            IF(ICHX.NE.2) THEN
               ZALB(IC)=0.5*(1.0-ZCODE(IC))/(1.0+ZCODE(IC))
            ELSE IF((ICHX.EQ.2).AND.(ZCODE(IC).NE.1.0)) THEN
               ZALB(IC)=2.0*(1.0+ZCODE(IC))/(1.0-ZCODE(IC))
            ELSE IF((ICHX.EQ.2).AND.(ZCODE(IC).EQ.1.0)) THEN
               ZALB(IC)=1.0E20
            ENDIF
         ENDDO
         IF(NALBP.GT.0) THEN
            DO IC=1,6
               IALB=ICODE(IC)
               IF(IALB.NE.0) ZALB(IC)=GAMMA(IALB)
            ENDDO
         ENDIF
         CALL TRICYL(MAXMIX,IMPX,ICHX,IDIM,LX,LY,LZ,XX,YY,ZZ,VOL,MAT,
     1   NCODE,ZALB,NR0,RR0,XR0,ANG,SGD,QFR)
         DEALLOCATE(ANG,XR0,RR0)
      ELSE IF(NALBP.GT.0) THEN
         IF((IPR.GT.0).AND.(ICHX.NE.2)) CALL XABORT('TRIASM: PERTURBAT'
     1   //'ION CALCULATION NOT AVAILABLE WITH PHYSICAL ALBEDOS.')
         DO IQW=1,MAXQF
            IALB=IQFR(IQW)
            IF(IALB.NE.0) QFR(IQW)=QFR(IQW)*GAMMA(IALB)
         ENDDO
      ELSE IF(IPR.GT.0) THEN
         CALL XDRSET(QFR,MAXQF,0.0)
      ENDIF
*----
*  ASSEMBLY OF THE ADI SPLITTED SYSTEM MATRICES
*----
*
* DIMENSION W
      IF(CHEX) THEN
         IF((ICHX.EQ.3).AND.(ISPLH.GT.1)) THEN
            ALLOCATE(MATN(LL4))
            NUM1=0
            DO 110 I=1,LX*LZ
            IF(MAT(I).EQ.0) GO TO 110
            DO 100 J=1,6*(ISPLH-1)**2
            KEL=KN(NUM1+J)
            MATN(KEL)=MAT(I)
  100       CONTINUE
            NUM1=NUM1+18*(ISPLH-1)**2+8
  110       CONTINUE
         ENDIF
         IIMAW=MUW(LL4W)
         IF(IPR.NE.3) THEN
            IF((IPR.EQ.0).OR.(ICHX.NE.2)) THEN
               WA_PTR=LCMARA(IIMAW)
               CALL C_F_POINTER(WA_PTR,WA,(/ IIMAW /))
            ELSE
               ALLOCATE(WA(IIMAW))
            ENDIF
            CALL XDRSET(WA,IIMAW,0.0)
         ELSE
            IF(ISEG.GT.0) CALL MTBLD('W_'//TEXT10,IPTRK,IPSYS,1)
            CALL LCMGPD(IPSYS,'W_'//TEXT10,WA_PTR)
            CALL C_F_POINTER(WA_PTR,WA,(/ IIMAW /))
         ENDIF
         IF(ICHX.EQ.1) THEN
*           MESH CORNER FINITE DIFFERENCES IN HEXAGONAL GEOMETRY.
            CALL LCMGPD(IPTRK,'IPW',IPW_PTR)
            CALL C_F_POINTER(IPW_PTR,IPW,(/ LL4 /))
            CALL TRIRWW(MAXMIX,NEL,LL4,VOL,MAT,XSGD,SIDE,ZZ,KN,QFR,MUW,
     1      WA,ISPLH,R,Q,RH,QH,RT,QT)
         ELSE IF(ICHX.EQ.2) THEN
*           THOMAS-RAVIART-SCHNEIDER FINITE ELEMENTS IN HEXAGONAL
*           GEOMETRY.
            IF(IPR.NE.3) THEN
               TF_PTR=LCMARA(LL4F)
               AW_PTR=LCMARA(IIMAW)
               CALL C_F_POINTER(TF_PTR,TF,(/ LL4F /))
               CALL C_F_POINTER(AW_PTR,AW,(/ IIMAW /))
               CALL XDRSET(TF,LL4F,0.0)
               CALL XDRSET(AW,IIMAW,0.0)
            ELSE
               IF(ISEG.GT.0) CALL MTBLD('WA'//TEXT10,IPTRK,IPSYS,1)
               CALL LCMGPD(IPSYS,'TF'//TEXT10,TF_PTR)
               CALL LCMGPD(IPSYS,'WA'//TEXT10,AW_PTR)
               CALL C_F_POINTER(TF_PTR,TF,(/ LL4F /))
               CALL C_F_POINTER(AW_PTR,AW,(/ IIMAW /))
            ENDIF
            NBLOS=LX*LZ/3
            ALLOCATE(IPERT(NBLOS),FRZ(NBLOS),DIF(NBLOS))
            CALL LCMGPD(IPTRK,'IPBBW',IPBW_PTR)
            CALL LCMGPD(IPTRK,'WB',BW_PTR)
            CALL C_F_POINTER(IPBW_PTR,IPBW,(/ 2*IELEM*LL4W /))
            CALL C_F_POINTER(BW_PTR,BW,(/ 2*IELEM*LL4W /))
            CALL LCMGET(IPTRK,'IPERT',IPERT)
            CALL LCMGET(IPTRK,'FRZ',FRZ)
            DO 120 KEL=1,NBLOS
            DIF(KEL)=0.0
            IF(IPERT(KEL).GT.0) THEN
               IBM=MAT((IPERT(KEL)-1)*3+1)
               DZ=ZZ((IPERT(KEL)-1)*3+1)*FRZ(KEL)
               IF(IBM.GT.0) DIF(KEL)=DZ/SGD(IBM,1)
            ENDIF
  120       CONTINUE
            CALL LCMPUT(IPSYS,'DIFF'//TEXT10,NBLOS,2,DIF)
            CALL TRIHWW(MAXMIX,NBLOS,IELEM,LL4F,LL4W,MAT,SIDE,ZZ,FRZ,
     1      QFR,IPERT,KN,XSGD,DSGD,MUW,IPBW,LC,R,V,BW,TF,AW,WA)
            DEALLOCATE(DIF,FRZ,IPERT)
         ELSE IF(ICHX.EQ.3) THEN
*           MESH CENTERED FINITE DIFFERENCES IN HEXAGONAL GEOMETRY.
            CALL LCMGPD(IPTRK,'IPW',IPW_PTR)
            CALL C_F_POINTER(IPW_PTR,IPW,(/ LL4 /))
            IF(ISPLH.EQ.1) THEN
               CALL TRIMWW(MAXMIX,NEL,LL4,VOL,MAT,SGD,XSGD,SIDE,ZZ,KN,
     1         QFR,MUW,IPW,IPR,WA)
            ELSE
               CALL TRIMTW(ISPLH,MAXMIX,NEL,LL4,VOL,MAT,MATN,SGD,XSGD,
     1         SIDE,ZZ,KN,QFR,MUW,IPW,IPR,WA)
            ENDIF
         ENDIF
         IF((IPR.EQ.0).OR.(IPR.EQ.3).OR.(ICHX.NE.2)) THEN
            CALL LCMPPD(IPSYS,'W_'//TEXT10,IIMAW,2,WA_PTR)
         ELSE
            DEALLOCATE(WA)
         ENDIF
         IF(ICHX.EQ.2) THEN
            CALL LCMPPD(IPSYS,'WA'//TEXT10,IIMAW,2,AW_PTR)
            CALL LCMPPD(IPSYS,'TF'//TEXT10,LL4F,2,TF_PTR)
         ENDIF
      ENDIF
*
* DIMENSION X
      IIMAX=MUX(LL4X)
      IF(CHEX.AND.(ICHX.EQ.2)) THEN
*        THOMAS-RAVIART-SCHNEIDER FINITE ELEMENTS IN HEXAGONAL GEOMETRY.
         IF(IPR.NE.3) THEN
            AX_PTR=LCMARA(IIMAX)
            CALL C_F_POINTER(AX_PTR,AX,(/ IIMAX /))
            CALL XDRSET(AX,IIMAX,0.0)
         ELSE
            IF(ISEG.GT.0) CALL MTBLD('XA'//TEXT10,IPTRK,IPSYS,1)
            CALL LCMGPD(IPSYS,'XA'//TEXT10,AX_PTR)
            CALL C_F_POINTER(AX_PTR,AX,(/ IIMAX /))
         ENDIF
         NBLOS=LX*LZ/3
         ALLOCATE(IPERT(NBLOS),FRZ(NBLOS))
         CALL LCMGPD(IPSYS,'TF'//TEXT10,TF_PTR)
         CALL C_F_POINTER(TF_PTR,TF,(/ LL4F /))
         CALL LCMGPD(IPTRK,'IPBBX',IPBX_PTR)
         CALL LCMGPD(IPTRK,'XB',BX_PTR)
         CALL C_F_POINTER(IPBX_PTR,IPBX,(/ 2*IELEM*LL4X /))
         CALL C_F_POINTER(BX_PTR,BX,(/ 2*IELEM*LL4X /))
         CALL LCMGET(IPTRK,'IPERT',IPERT)
         CALL LCMGET(IPTRK,'FRZ',FRZ)
         IF((IPR.EQ.0).OR.(IPR.EQ.3)) THEN
            XA_PTR=LCMARA(IIMAX)
            CALL C_F_POINTER(XA_PTR,XA,(/ IIMAX /))
         ELSE
            ALLOCATE(XA(IIMAX))
         ENDIF
         CALL TRIHWX(MAXMIX,NBLOS,IELEM,LL4F,LL4W,LL4X,MAT,SIDE,ZZ,FRZ,
     1   QFR,IPERT,KN,DSGD,MUX,IPBX,LC,R,BX,TF,AX,XA)
         DEALLOCATE(FRZ,IPERT)
      ELSE IF(ICHX.EQ.2) THEN
*        THOMAS-RAVIART ADI ITERATIVE METHOD.
         IF(DIAG) THEN
            ALLOCATE(AX(IIMAX))
            IF(IPR.NE.3) THEN
               TF_PTR=LCMARA(LL4F)
               CALL C_F_POINTER(TF_PTR,TF,(/ LL4F /))
               CALL XDRSET(TF,LL4F,0.0)
               CALL XDRSET(AX,IIMAX,0.0)
            ELSE
               IF(ISEG.GT.0) CALL MTBLD('XA'//TEXT10,IPTRK,IPSYS,1)
               CALL LCMGPD(IPSYS,'TF'//TEXT10,TF_PTR)
               CALL C_F_POINTER(TF_PTR,TF,(/ LL4F /))
               CALL LCMGET(IPSYS,'XA'//TEXT10,AX)
            ENDIF
            ALLOCATE(XA(IIMAX))
         ELSE
            IF(IPR.NE.3) THEN
               TF_PTR=LCMARA(LL4F)
               AX_PTR=LCMARA(IIMAX)
               CALL C_F_POINTER(TF_PTR,TF,(/ LL4F /))
               CALL C_F_POINTER(AX_PTR,AX,(/ IIMAX /))
               CALL XDRSET(TF,LL4F,0.0)
               CALL XDRSET(AX,IIMAX,0.0)
            ELSE
               IF(ISEG.GT.0) CALL MTBLD('XA'//TEXT10,IPTRK,IPSYS,1)
               CALL LCMGPD(IPSYS,'TF'//TEXT10,TF_PTR)
               CALL LCMGPD(IPSYS,'XA'//TEXT10,AX_PTR)
               CALL C_F_POINTER(TF_PTR,TF,(/ LL4F /))
               CALL C_F_POINTER(AX_PTR,AX,(/ IIMAX /))
            ENDIF
            IF((IPR.EQ.0).OR.(IPR.EQ.3)) THEN
               XA_PTR=LCMARA(IIMAX)
               CALL C_F_POINTER(XA_PTR,XA,(/ IIMAX /))
            ELSE
               ALLOCATE(XA(IIMAX))
            ENDIF
         ENDIF
         CALL LCMGPD(IPTRK,'IPBBX',IPBX_PTR)
         CALL LCMGPD(IPTRK,'XB',BX_PTR)
         CALL C_F_POINTER(IPBX_PTR,IPBX,(/ 2*IELEM*LL4X /))
         CALL C_F_POINTER(BX_PTR,BX,(/ 2*IELEM*LL4X /))
         CALL TRIDXX(MAXMIX,CYLIND,IELEM,ICOL,NEL,LL4F,LL4X,MAT,VOL2,
     1   XX,YY,ZZ,DD,KN,QFR,XSGD,DSGD,MUX,IPBX,LC,R,V,BX,TF,AX,XA)
      ELSE
*        GENERIC ADI ITERATIVE METHOD.
         CALL LCMGPD(IPTRK,'IPX',IPX_PTR)
         CALL C_F_POINTER(IPX_PTR,IPX,(/ LL4 /))
         IF(DIAG) THEN
            ALLOCATE(XA(IIMAX))
            CALL XDRSET(XA,IIMAX,0.0)
         ELSE IF(IPR.NE.3) THEN
            XA_PTR=LCMARA(IIMAX)
            CALL C_F_POINTER(XA_PTR,XA,(/ IIMAX /))
            CALL XDRSET(XA,IIMAX,0.0)
         ELSE
            IF(ISEG.GT.0) CALL MTBLD('X_'//TEXT10,IPTRK,IPSYS,1)
            CALL LCMGPD(IPSYS,'X_'//TEXT10,XA_PTR)
            CALL C_F_POINTER(XA_PTR,XA,(/ IIMAX /))
         ENDIF
         IF((ICHX.EQ.1).AND.(.NOT.CHEX)) THEN
            CALL TRIPXX(MAXMIX,MAXKN,NEL,LL4,VOL2,MAT,XSGD,XX,YY,ZZ,DD,
     1      KN,QFR,MUX,IPX,CYLIND,LC,T,TS,Q,QS,XA)
         ELSE IF((ICHX.EQ.3).AND.(.NOT.CHEX)) THEN
            CALL TRIMXX(MAXMIX,CYLIND,IELEM,IDIM,NEL,LL4,VOL2,MAT,SGD,
     1      XSGD,XX,YY,ZZ,DD,KN,QFR,MUX,IPX,IPR,XA)
         ELSE IF((ICHX.EQ.1).AND.CHEX) THEN
*           MESH CORNER FINITE DIFFERENCES IN HEXAGONAL GEOMETRY.
            CALL TRIRWX(MAXMIX,NEL,LL4,VOL,MAT,XSGD,SIDE,ZZ,KN,QFR,
     1      MUX,IPX,XA,ISPLH,R,Q,RH,QH,RT,QT)
         ELSE IF((ICHX.EQ.3).AND.CHEX) THEN
*           MESH CENTERED FINITE DIFFERENCES IN HEXAGONAL GEOMETRY.
            IF(ISPLH.EQ.1) THEN
               CALL TRIMWX(MAXMIX,NEL,LL4,VOL,MAT,SGD,XSGD,SIDE,ZZ,KN,
     1         QFR,MUX,IPX,IPR,XA)
            ELSE
               CALL TRIMTX(ISPLH,MAXMIX,NEL,LL4,VOL,MAT,MATN,SGD,XSGD,
     1         SIDE,ZZ,KN,QFR,MUX,IPX,IPR,XA)
            ENDIF
         ENDIF
      ENDIF
      IF(.NOT.DIAG) THEN
         IF((IPR.EQ.0).OR.(IPR.EQ.3).OR.(ICHX.NE.2)) THEN
            CALL LCMPPD(IPSYS,'X_'//TEXT10,IIMAX,2,XA_PTR)
         ELSE
            DEALLOCATE(XA)
         ENDIF
         IF(ICHX.EQ.2) CALL LCMPPD(IPSYS,'XA'//TEXT10,IIMAX,2,AX_PTR)
      ELSE
*        IN DIAGONAL SYMMETRY CASE, DO NOT SAVE THE X-DIRECTED ADI
*        MATRIX COMPONENT SINCE IT IS EQUAL TO THE Y-DIRECTED COMPONENT
         DEALLOCATE(XA)
         IF(ICHX.EQ.2) DEALLOCATE(AX)
      ENDIF
      IF(.NOT.CHEX.AND.(ICHX.EQ.2)) CALL LCMPPD(IPSYS,'TF'//TEXT10,LL4F,
     1 2,TF_PTR)
*
* DIMENSION Y
      IF(LOGY) THEN
         IIMAY=MUY(LL4Y)
         IF(CHEX.AND.(ICHX.EQ.2)) THEN
*           THOMAS-RAVIART-SCHNEIDER FINITE ELEMENTS IN HEXAGONAL
*           GEOMETRY.
            IF(IPR.NE.3) THEN
               AY_PTR=LCMARA(IIMAY)
               CALL C_F_POINTER(AY_PTR,AY,(/ IIMAY /))
               CALL XDRSET(AY,IIMAY,0.0)
            ELSE
               IF(ISEG.GT.0) CALL MTBLD('YA'//TEXT10,IPTRK,IPSYS,1)
               CALL LCMGPD(IPSYS,'YA'//TEXT10,AY_PTR)
               CALL C_F_POINTER(AY_PTR,AY,(/ IIMAY /))
            ENDIF
            NBLOS=LX*LZ/3
            ALLOCATE(IPERT(NBLOS),FRZ(NBLOS))
            CALL LCMGPD(IPSYS,'TF'//TEXT10,TF_PTR)
            CALL C_F_POINTER(TF_PTR,TF,(/ LL4F /))
            CALL LCMGPD(IPTRK,'IPBBY',IPBY_PTR)
            CALL LCMGPD(IPTRK,'YB',BY_PTR)
            CALL C_F_POINTER(IPBY_PTR,IPBY,(/ 2*IELEM*LL4Y /))
            CALL C_F_POINTER(BY_PTR,BY,(/ 2*IELEM*LL4Y /))
            CALL LCMGET(IPTRK,'IPERT',IPERT)
            CALL LCMGET(IPTRK,'FRZ',FRZ)
            IF((IPR.EQ.0).OR.(IPR.EQ.3)) THEN
               YA_PTR=LCMARA(IIMAY)
               CALL C_F_POINTER(YA_PTR,YA,(/ IIMAY /))
            ELSE
               ALLOCATE(YA(IIMAY))
            ENDIF
            CALL TRIHWY(MAXMIX,NBLOS,IELEM,LL4F,LL4W,LL4X,LL4Y,MAT,
     1      SIDE,ZZ,FRZ,QFR,IPERT,KN,DSGD,MUY,IPBY,LC,R,BY,TF,AY,YA)
            DEALLOCATE(FRZ,IPERT)
         ELSE IF(ICHX.EQ.2) THEN
*           THOMAS-RAVIART ADI ITERATIVE METHOD.
            IF(IPR.NE.3) THEN
               AY_PTR=LCMARA(IIMAY)
               CALL C_F_POINTER(AY_PTR,AY,(/ IIMAY /))
               CALL XDRSET(AY,IIMAY,0.0)
            ELSE
               IF(ISEG.GT.0) CALL MTBLD('YA'//TEXT10,IPTRK,IPSYS,1)
               CALL LCMGPD(IPSYS,'YA'//TEXT10,AY_PTR)
               CALL C_F_POINTER(AY_PTR,AY,(/ IIMAY /))
            ENDIF
            CALL LCMGPD(IPSYS,'TF'//TEXT10,TF_PTR)
            CALL LCMGPD(IPTRK,'IPBBY',IPBY_PTR)
            CALL LCMGPD(IPTRK,'YB',BY_PTR)
            CALL C_F_POINTER(TF_PTR,TF,(/ LL4F /))
            CALL C_F_POINTER(IPBY_PTR,IPBY,(/ 2*IELEM*LL4Y /))
            CALL C_F_POINTER(BY_PTR,BY,(/ 2*IELEM*LL4Y /))
            IF((IPR.EQ.0).OR.(IPR.EQ.3)) THEN
               YA_PTR=LCMARA(IIMAY)
               CALL C_F_POINTER(YA_PTR,YA,(/ IIMAY /))
            ELSE
               ALLOCATE(YA(IIMAY))
            ENDIF
            CALL TRIDXY(MAXMIX,IELEM,ICOL,NEL,LL4F,LL4X,LL4Y,MAT,VOL2,
     1      YY,KN,QFR,DSGD,MUY,IPBY,LC,R,BY,TF,AY,YA)
         ELSE
*           GENERIC ADI ITERATIVE METHOD.
            CALL LCMGPD(IPTRK,'IPY',IPY_PTR)
            CALL C_F_POINTER(IPY_PTR,IPY,(/ LL4 /))
            IF(IPR.NE.3) THEN
               YA_PTR=LCMARA(IIMAY)
               CALL C_F_POINTER(YA_PTR,YA,(/ IIMAY /))
               CALL XDRSET(YA,IIMAY,0.0)
             ELSE
               IF(ISEG.GT.0) CALL MTBLD('Y_'//TEXT10,IPTRK,IPSYS,1)
               CALL LCMGPD(IPSYS,'Y_'//TEXT10,YA_PTR)
               CALL C_F_POINTER(YA_PTR,YA,(/ IIMAY /))
            ENDIF
            IF((ICHX.EQ.1).AND.(.NOT.CHEX)) THEN
               CALL TRIPXY(MAXMIX,MAXKN,NEL,LL4,VOL2,MAT,XSGD,XX,YY,ZZ,
     1         DD,KN,QFR,MUY,IPY,CYLIND,LC,T,TS,Q,QS,YA)
            ELSE IF((ICHX.EQ.3).AND.(.NOT.CHEX)) THEN
               CALL TRIMXY(MAXMIX,CYLIND,IELEM,IDIM,NEL,LL4,VOL2,MAT,
     1         SGD,XSGD,XX,YY,ZZ,DD,KN,QFR,MUY,IPY,IPR,YA)
            ELSE IF((ICHX.EQ.1).AND.CHEX) THEN
*              MESH CORNER FINITE DIFFERENCES IN HEXAGONAL GEOMETRY.
               CALL TRIRWY(MAXMIX,NEL,LL4,VOL,MAT,XSGD,SIDE,ZZ,
     1         KN,QFR,MUY,IPY,YA,ISPLH,R,Q,RH,QH,RT,QT)
            ELSE IF((ICHX.EQ.3).AND.CHEX) THEN
*              MESH CENTERED FINITE DIFFERENCES IN HEXAGONAL GEOMETRY.
               IF(ISPLH.EQ.1) THEN
                  CALL TRIMWY(MAXMIX,NEL,LL4,VOL,MAT,SGD,XSGD,SIDE,ZZ,
     1            KN,QFR,MUY,IPY,IPR,YA)
               ELSE
                  CALL TRIMTY(ISPLH,MAXMIX,NEL,LL4,VOL,MAT,MATN,SGD,
     1            XSGD,SIDE,ZZ,KN,QFR,MUY,IPY,IPR,YA)
               ENDIF
            ENDIF
         ENDIF
         IF((IPR.EQ.0).OR.(IPR.EQ.3).OR.(ICHX.NE.2)) THEN
            CALL LCMPPD(IPSYS,'Y_'//TEXT10,IIMAY,2,YA_PTR)
         ELSE
            DEALLOCATE(YA)
         ENDIF
         IF(ICHX.EQ.2) CALL LCMPPD(IPSYS,'YA'//TEXT10,IIMAY,2,AY_PTR)
      ENDIF
*
* DIMENSION Z
      IF(LOGZ) THEN
         IIMAZ=MUZ(LL4Z)
         IF(CHEX.AND.(ICHX.EQ.2)) THEN
*           THOMAS-RAVIART-SCHNEIDER FINITE ELEMENTS IN HEXAGONAL
*           GEOMETRY.
            IF(IPR.NE.3) THEN
               AZ_PTR=LCMARA(IIMAZ)
               CALL C_F_POINTER(AZ_PTR,AZ,(/ IIMAZ /))
               CALL XDRSET(AZ,IIMAZ,0.0)
            ELSE
               IF(ISEG.GT.0) CALL MTBLD('ZA'//TEXT10,IPTRK,IPSYS,1)
               CALL LCMGPD(IPSYS,'ZA'//TEXT10,AZ_PTR)
               CALL C_F_POINTER(AZ_PTR,AZ,(/ IIMAZ /))
            ENDIF
            NBLOS=LX*LZ/3
            ALLOCATE(IPERT(NBLOS),FRZ(NBLOS))
            CALL LCMGPD(IPSYS,'TF'//TEXT10,TF_PTR)
            CALL C_F_POINTER(TF_PTR,TF,(/ LL4F /))
            CALL LCMGPD(IPTRK,'IPBBZ',IPBZ_PTR)
            CALL LCMGPD(IPTRK,'ZB',BZ_PTR)
            CALL C_F_POINTER(IPBZ_PTR,IPBZ,(/ 2*IELEM*LL4Z /))
            CALL C_F_POINTER(BZ_PTR,BZ,(/ 2*IELEM*LL4Z /))
            CALL LCMGET(IPTRK,'IPERT',IPERT)
            CALL LCMGET(IPTRK,'FRZ',FRZ)
            IF((IPR.EQ.0).OR.(IPR.EQ.3)) THEN
               ZA_PTR=LCMARA(IIMAZ)
               CALL C_F_POINTER(ZA_PTR,ZA,(/ IIMAZ /))
            ELSE
               ALLOCATE(ZA(IIMAZ))
            ENDIF
            CALL TRIHWZ(MAXMIX,NBLOS,IELEM,ICOL,LL4F,LL4W,LL4X,LL4Y,
     1      LL4Z,MAT,SIDE,ZZ,FRZ,QFR,IPERT,KN,DSGD,MUZ,IPBZ,LC,R,BZ,
     2      TF,AZ,ZA)
            DEALLOCATE(FRZ,IPERT)
         ELSE IF(ICHX.EQ.2) THEN
*           THOMAS-RAVIART ADI ITERATIVE METHOD.
            IF(IPR.NE.3) THEN
               AZ_PTR=LCMARA(IIMAZ)
               CALL C_F_POINTER(AZ_PTR,AZ,(/ IIMAZ /))
               CALL XDRSET(AZ,IIMAZ,0.0)
            ELSE
               IF(ISEG.GT.0) CALL MTBLD('ZA'//TEXT10,IPTRK,IPSYS,1)
               CALL LCMGPD(IPSYS,'ZA'//TEXT10,AZ_PTR)
               CALL C_F_POINTER(AZ_PTR,AZ,(/ IIMAZ /))
            ENDIF
            CALL LCMGPD(IPSYS,'TF'//TEXT10,TF_PTR)
            CALL LCMGPD(IPTRK,'IPBBZ',IPBZ_PTR)
            CALL LCMGPD(IPTRK,'ZB',BZ_PTR)
            CALL C_F_POINTER(TF_PTR,TF,(/ LL4F /))
            CALL C_F_POINTER(IPBZ_PTR,IPBZ,(/ 2*IELEM*LL4Z /))
            CALL C_F_POINTER(BZ_PTR,BZ,(/ 2*IELEM*LL4Z /))
            IF((IPR.EQ.0).OR.(IPR.EQ.3)) THEN
               ZA_PTR=LCMARA(IIMAZ)
               CALL C_F_POINTER(ZA_PTR,ZA,(/ IIMAZ /))
            ELSE
               ALLOCATE(ZA(IIMAZ))
            ENDIF
            CALL TRIDXZ(MAXMIX,IELEM,ICOL,NEL,LL4F,LL4X,LL4Y,LL4Z,MAT,
     1      VOL2,ZZ,KN,QFR,DSGD,MUZ,IPBZ,LC,R,BZ,TF,AZ,ZA)
         ELSE
            CALL LCMGPD(IPTRK,'IPZ',IPZ_PTR)
            CALL C_F_POINTER(IPZ_PTR,IPZ,(/ LL4 /))
            IF(IPR.NE.3) THEN
               ZA_PTR=LCMARA(IIMAZ)
               CALL C_F_POINTER(ZA_PTR,ZA,(/ IIMAZ /))
               CALL XDRSET(ZA,IIMAZ,0.0)
            ELSE
               IF(ISEG.GT.0) CALL MTBLD('Z_'//TEXT10,IPTRK,IPSYS,1)
               CALL LCMGPD(IPSYS,'Z_'//TEXT10,ZA_PTR)
               CALL C_F_POINTER(ZA_PTR,ZA,(/ IIMAZ /))
            ENDIF
            IF((ICHX.EQ.1).AND.(.NOT.CHEX)) THEN
               CALL TRIPXZ(MAXMIX,MAXKN,NEL,LL4,VOL2,MAT,XSGD,XX,YY,ZZ,
     1         DD,KN,QFR,MUZ,IPZ,CYLIND,LC,T,TS,Q,QS,ZA)
            ELSE IF((ICHX.EQ.3).AND.(.NOT.CHEX)) THEN
               CALL TRIMXZ(MAXMIX,CYLIND,IELEM,NEL,LL4,VOL2,MAT,SGD,
     1         XSGD,XX,YY,ZZ,DD,KN,QFR,MUZ,IPZ,IPR,ZA)
            ELSE IF((ICHX.EQ.1).AND.CHEX) THEN
*              MESH CORNER FINITE DIFFERENCES IN HEXAGONAL GEOMETRY.
               CALL TRIRWZ(MAXMIX,NEL,LL4,VOL,MAT,XSGD,SIDE,ZZ,KN,QFR,
     1         MUZ,IPZ,ZA,ISPLH,R,Q,RH,QH,RT,QT)
            ELSE IF((ICHX.EQ.3).AND.CHEX) THEN
*              MESH CENTERED FINITE DIFFERENCES IN HEXAGONAL GEOMETRY.
               IF(ISPLH.EQ.1) THEN
                  CALL TRIMWZ(MAXMIX,NEL,LL4,VOL,MAT,SGD,XSGD,SIDE,ZZ,
     1            KN,QFR,MUZ,IPZ,IPR,ZA)
               ELSE
                  CALL TRIMTZ(ISPLH,MAXMIX,NEL,LL4,VOL,MAT,MATN,SGD,
     1            XSGD,SIDE,ZZ,KN,QFR,MUZ,IPZ,IPR,ZA)
               ENDIF
            ENDIF
         ENDIF
         IF((IPR.EQ.0).OR.(IPR.EQ.3).OR.(ICHX.NE.2)) THEN
            CALL LCMPPD(IPSYS,'Z_'//TEXT10,IIMAZ,2,ZA_PTR)
         ELSE
            DEALLOCATE(ZA)
         ENDIF
         IF(ICHX.EQ.2) CALL LCMPPD(IPSYS,'ZA'//TEXT10,IIMAZ,2,AZ_PTR)
      ENDIF
      DEALLOCATE(VOL2)
      IF(ICHX.EQ.2) DEALLOCATE(DSGD)
      IF((ICHX.EQ.3).AND.(ISPLH.GT.1).AND.CHEX) DEALLOCATE(MATN)
*----
*  CHECK FOR MATRIX CONSISTENCY
*----
      IF(ICHX.NE.2) CALL TRICHK (TEXT10,IPTRK,IPSYS,IDIM,DIAG,CHEX,
     1 IPR,LL4)
      CALL KDRCPU(TK2)
      IF(IMPX.GT.1) WRITE(6,'(/35H TRIASM: CPU TIME FOR SYSTEM MATRIX,
     1 11H ASSEMBLY =,F9.2,3H S.)') TK2-TK1
*----
*  PERFORM SUPERVECTORIZATION REBUILD OF THE COEFFICIENT MATRICES
*----
      IF(ISEG.GT.0) THEN
         IF((IPR.EQ.0).OR.(IPR.EQ.3).OR.(ICHX.NE.2)) THEN
            IF(CHEX) CALL MTBLD('W_'//TEXT10,IPTRK,IPSYS,3)
            IF(.NOT.DIAG) CALL MTBLD('X_'//TEXT10,IPTRK,IPSYS,3)
            IF(LOGY) CALL MTBLD('Y_'//TEXT10,IPTRK,IPSYS,3)
            IF(LOGZ) CALL MTBLD('Z_'//TEXT10,IPTRK,IPSYS,3)
         ENDIF
         IF(ICHX.EQ.2) THEN
            IF(CHEX) CALL MTBLD('WA'//TEXT10,IPTRK,IPSYS,3)
            IF(.NOT.DIAG) CALL MTBLD('XA'//TEXT10,IPTRK,IPSYS,3)
            IF(LOGY) CALL MTBLD('YA'//TEXT10,IPTRK,IPSYS,3)
            IF(LOGZ) CALL MTBLD('ZA'//TEXT10,IPTRK,IPSYS,3)
         ENDIF
      ENDIF
*----
*  MATRIX FACTORIZATIONS
*----
      IF((IPR.EQ.0).OR.(IPR.EQ.3)) THEN
         CALL KDRCPU(TK1)
         CALL MTLDLF(TEXT10,IPTRK,IPSYS,ITY,IMPX)
         CALL KDRCPU(TK2)
         IF(IMPX.GT.1) WRITE(6,'(/34H TRIASM: CPU TIME FOR LDLT FACTORI,
     1   18HZATION OF MATRIX '',A10,2H''=,F9.2,3H S.)') TEXT10,TK2-TK1
      ENDIF
*----
*  RELEASE UNIT MATRICES
*----
      IF((ICHX.EQ.1).OR.(ICHX.EQ.2)) THEN
         DEALLOCATE(T,TS,R,RS,Q,QS,V,RH,QH,RT,QT)
      ENDIF
*----
*  RELEASE TRIVAC SPECIFIC TRACKING INFORMATION
*----
      DEALLOCATE(IQFR,QFR,KN,ZZ)
      IF(CHEX) THEN
         DEALLOCATE(MUW)
      ELSE
         DEALLOCATE(DD,YY,XX)
      ENDIF
      IF(LOGY) DEALLOCATE(MUY)
      IF(.NOT.DIAG) DEALLOCATE(MUX)
      IF(LOGZ) DEALLOCATE(MUZ)
      RETURN
      END
