*DECK TRIASN
      SUBROUTINE TRIASN(HNAMT,IPTRK,IPSYS,IMPX,NBMIX,NEL,NAN,NALBP,IPR,
     1 MAT,VOL,GAMMA,SIGT,SIGTI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembly of a single-group system matrix with leakage and removal
* cross sections for the simplified PN method.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
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
* NBMIX   number of mixtures.
* NEL     total number of finite elements.
* NAN     number of Legendre orders for the cross sections.
* NALBP   number of physical albedos.
* IPR     type of assembly:
*         =0: calculation of the system matrices;
*         =1: calculation of the derivative of these matrices;
*         =2: calculation of the first variation of these matrices;
*         =3: identical to IPR=2, but these variation are added to
*         unperturbed system matrices.
* MAT     index-number of the mixture type assigned to each volume.
* GAMMA   physical albedo functions.
* VOL     volumes.
* SIGT    total minus self-scattering macroscopic cross sections.
*         SIGT(:,NAN) generally contains the total cross section only.
*         If IPR.gt.0, SIGT contains perturbed or derivative values.
* SIGTI   inverse macroscopic cross sections ordered by mixture.
*         SIGTI(:,NAN) generally contains the inverse total cross
*         section only. If IPR.gt.0, SIGTI contains perturbed or
*         derivative values.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPSYS
      CHARACTER HNAMT*10
      INTEGER IMPX,NBMIX,NEL,NAN,IPR,MAT(NEL)
      REAL VOL(NEL),GAMMA(NALBP),SIGT(NBMIX,NAN),SIGTI(NBMIX,NAN)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      LOGICAL CYLIND,CHEX,DIAG,LSGD
      CHARACTER TEXT10*10
      INTEGER ISTATE(NSTATE)
      INTEGER, DIMENSION(:), ALLOCATABLE :: KN,IQFR,MUW,MUZ,IPERT
      INTEGER, DIMENSION(:), POINTER :: MUX
      INTEGER, DIMENSION(:), ALLOCATABLE, TARGET :: MUY
      REAL, DIMENSION(:), ALLOCATABLE :: VOL2,XX,YY,ZZ,QFR,FRZ,DIF
      REAL, DIMENSION(:,:), ALLOCATABLE :: R,V
      INTEGER, DIMENSION(:), POINTER :: IPBW,IPBX,IPBY,IPBZ
      REAL, DIMENSION(:), POINTER :: TF,AW,AX,AY,AZ,WA,XA,YA,ZA,BW,BX,
     1 BY,BZ
      TYPE(C_PTR) IPBW_PTR,IPBX_PTR,IPBY_PTR,IPBZ_PTR
      TYPE(C_PTR) TF_PTR,AW_PTR,AX_PTR,AY_PTR,AZ_PTR,WA_PTR,XA_PTR,
     1 YA_PTR,ZA_PTR,BW_PTR,BX_PTR,BY_PTR,BZ_PTR
*----
*  RECOVER TRIVAC SPECIFIC TRACKING INFORMATION
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      ITYPE=ISTATE(6)
      CYLIND=(ITYPE.EQ.3).OR.(ITYPE.EQ.6)
      CHEX=(ITYPE.EQ.8).OR.(ITYPE.EQ.9)
      IF(CYLIND) CALL XABORT('TRIASN: GEOMETRY NOT AVAILABLE.')
      IHEX=ISTATE(7)
      DIAG=(ISTATE(8).EQ.1)
      IELEM=ISTATE(9)
      ICOL=ISTATE(10)
      LL4=ISTATE(11)
      ICHX=ISTATE(12)
      IF(ICHX.NE.2) CALL XABORT('TRIASN: DISCRETIZATION NOT AVAILABLE.')
      ISPLH=ISTATE(13)
      LX=ISTATE(14)
      LY=ISTATE(15)
      LZ=ISTATE(16)
      ISEG=ISTATE(17)
      IMPV=ISTATE(18)
      NR0=ISTATE(24)
      LL4F=ISTATE(25)
      ITY=3
      LL4W=ISTATE(26)
      LL4X=ISTATE(27)
      LL4Y=ISTATE(28)
      LL4Z=ISTATE(29)
      NLF=ISTATE(30)
      NVD=ISTATE(34)
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
         ALLOCATE(XX(LX*LY*LZ),YY(LX*LY*LZ))
         CALL LCMGET(IPTRK,'XX',XX)
         CALL LCMGET(IPTRK,'YY',YY)
      ENDIF
      IF(LL4Y.GT.0) THEN
         ALLOCATE(MUY(LL4Y))
         CALL LCMGET(IPTRK,'MUY',MUY)
      ENDIF
      IF(.NOT.DIAG) THEN
         ALLOCATE(MUX(LL4X))
         CALL LCMGET(IPTRK,'MUX',MUX)
      ELSE
         MUX=>MUY
      ENDIF
      IF(LL4Z.GT.0) THEN
         ALLOCATE(MUZ(LL4Z))
         CALL LCMGET(IPTRK,'MUZ',MUZ)
      ENDIF
*----
*  RECOVER UNIT MATRICES
*----
      CALL LCMSIX(IPTRK,'BIVCOL',1)
      CALL LCMLEN(IPTRK,'T',LC,ITYLCM)
      ALLOCATE(R(LC,LC),V(LC,LC-1))
      CALL LCMGET(IPTRK,'R',R)
      CALL LCMGET(IPTRK,'V',V)
      CALL LCMSIX(IPTRK,' ',2)
*
      TEXT10=HNAMT(:10)
      IF(IMPX.GT.0) WRITE(6,'(/36H TRIASN: ASSEMBLY OF SYMMETRIC MATRI,
     1 3HX '',A10,38H'' IN COMPRESSED DIAGONAL STORAGE MODE.)') TEXT10
      CALL KDRCPU(TK1)
*----
*  DETERMINATION OF THE PERTURBED ELEMENTS. NON-PERTURBED ELEMENTS WILL
*  HAVE VOL(K)=0.0
*----
      ALLOCATE(VOL2(NEL))
      IF((IPR.EQ.0).OR.(NALBP.GT.0)) THEN
         DO 25 K=1,NEL
         VOL2(K)=VOL(K)
   25    CONTINUE
      ELSE
         CALL XDRSET(VOL2,NEL,0.0)
         DO 50 L=1,NEL
         IBM=MAT(L)
         IF(IBM.EQ.0) GO TO 50
         LSGD=.FALSE.
         DO 45 I=1,NAN
         LSGD=LSGD.OR.(SIGT(IBM,I).NE.0.0).OR.(SIGTI(IBM,I).NE.0.0)
   45    CONTINUE
         IF(LSGD) VOL2(L)=VOL(L)
   50    CONTINUE
      ENDIF
*----
*  APPLY PHYSICAL ALBEDOS AND INTRODUCE THE CYLINDER BOUNDARY
*  APPROXIMATION IN CARTESIAN GEOMETRY
*----
      IF(NR0.GT.0) THEN
         CALL XABORT('TRIASN: CYLINDRICAL CORRECTION NOT IMPLEMENTED.')
      ELSE IF(NALBP.GT.0) THEN
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
         IIMAW=MUW(LL4W)*NLF/2
         IF(DIAG.OR.(IPR.NE.3)) THEN
            TF_PTR=LCMARA(LL4F*NLF/2)
            AW_PTR=LCMARA(IIMAW)
            CALL C_F_POINTER(TF_PTR,TF,(/ LL4F*NLF/2 /))
            CALL C_F_POINTER(AW_PTR,AW,(/ IIMAW /))
            CALL XDRSET(TF,LL4F*NLF/2,0.0)
            CALL XDRSET(AW,IIMAW,0.0)
         ELSE
            IF(ISEG.GT.0) CALL MTBLD('WA'//TEXT10,IPTRK,IPSYS,1)
            CALL LCMGPD(IPSYS,'TF'//TEXT10,TF_PTR)
            CALL LCMGPD(IPSYS,'WA'//TEXT10,AW_PTR)
            CALL C_F_POINTER(TF_PTR,TF,(/ LL4F*NLF/2 /))
            CALL C_F_POINTER(AW_PTR,AW,(/ IIMAW /))
         ENDIF
         CALL LCMGPD(IPTRK,'IPBBW',IPBW_PTR)
         CALL LCMLEN(IPSYS,'WB',LENWB,ITYL)
         IF(LENWB.EQ.0) THEN
           CALL LCMGPD(IPTRK,'WB',BW_PTR)
         ELSE
           CALL LCMGPD(IPSYS,'WB',BW_PTR)
         ENDIF
         CALL C_F_POINTER(IPBW_PTR,IPBW,(/ 2*IELEM*LL4W /))
         CALL C_F_POINTER(BW_PTR,BW,(/ 2*IELEM*LL4W /))
         NBLOS=LX*LZ/3
         ALLOCATE(IPERT(NBLOS),FRZ(NBLOS),DIF(NBLOS))
         CALL LCMGET(IPTRK,'IPERT',IPERT)
         CALL LCMGET(IPTRK,'FRZ',FRZ)
         IF((IPR.EQ.0).OR.(IPR.EQ.3)) THEN
            WA_PTR=LCMARA(IIMAW)
            CALL C_F_POINTER(WA_PTR,WA,(/ IIMAW /))
         ELSE
            ALLOCATE(WA(IIMAW))
         ENDIF
         DO 60 KEL=1,NBLOS
         DIF(KEL)=0.0
         IF(IPERT(KEL).GT.0) THEN
            IBM=MAT((IPERT(KEL)-1)*3+1)
            DZ=ZZ((IPERT(KEL)-1)*3+1)*FRZ(KEL)
            IF(IBM.GT.0) DIF(KEL)=DZ*SIGT(IBM,1)
         ENDIF
   60    CONTINUE
         CALL LCMPUT(IPSYS,'SIGT'//TEXT10,NBLOS,2,DIF)
         CALL PN3HWW(NBMIX,NBLOS,IELEM,ICOL,NLF,NVD,NAN,LL4F,LL4W,MAT,
     1   SIGT,SIGTI,SIDE,ZZ,FRZ,QFR,IPERT,KN,MUW,IPBW,LC,R,V,BW,TF,AW,
     2   WA)
         IF((IPR.EQ.0).OR.(IPR.EQ.3)) THEN
            CALL LCMPPD(IPSYS,'W_'//TEXT10,IIMAW,2,WA_PTR)
         ELSE
            DEALLOCATE(WA)
         ENDIF
         CALL LCMPPD(IPSYS,'WA'//TEXT10,IIMAW,2,AW_PTR)
         CALL LCMPPD(IPSYS,'TF'//TEXT10,LL4F*NLF/2,2,TF_PTR)
         DEALLOCATE(DIF,FRZ,IPERT)
      ENDIF
*
* DIMENSION X
      IIMAX=MUX(LL4X)*NLF/2
      CALL LCMGPD(IPTRK,'IPBBX',IPBX_PTR)
      CALL LCMLEN(IPSYS,'XB',LENXB,ITYL)
      IF(LENXB.EQ.0) THEN
        CALL LCMGPD(IPTRK,'XB',BX_PTR)
      ELSE
        CALL LCMGPD(IPSYS,'XB',BX_PTR)
      ENDIF
      CALL C_F_POINTER(IPBX_PTR,IPBX,(/ 2*IELEM*LL4X /))
      CALL C_F_POINTER(BX_PTR,BX,(/ 2*IELEM*LL4X /))
      IF(CHEX) THEN
         IF(IPR.NE.3) THEN
            AX_PTR=LCMARA(IIMAX)
            CALL C_F_POINTER(AX_PTR,AX,(/ IIMAX /))
            CALL XDRSET(AX,IIMAX,0.0)
         ELSE         
            IF(ISEG.GT.0) CALL MTBLD('XA'//TEXT10,IPTRK,IPSYS,1)
            CALL LCMGPD(IPSYS,'XA'//TEXT10,AX_PTR)
            CALL C_F_POINTER(AX_PTR,AX,(/ IIMAX /))
         ENDIF
         CALL LCMGPD(IPSYS,'TF'//TEXT10,TF_PTR)
         CALL C_F_POINTER(TF_PTR,TF,(/ LL4F*NLF/2 /))
         NBLOS=LX*LZ/3
         ALLOCATE(IPERT(NBLOS),FRZ(NBLOS))
         CALL LCMGET(IPTRK,'IPERT',IPERT)
         CALL LCMGET(IPTRK,'FRZ',FRZ)
         IF((IPR.EQ.0).OR.(IPR.EQ.3)) THEN
            XA_PTR=LCMARA(IIMAX)
            CALL C_F_POINTER(XA_PTR,XA,(/ IIMAX /))
         ELSE
            ALLOCATE(XA(IIMAX))
         ENDIF
         CALL PN3HWX(NBMIX,NBLOS,IELEM,ICOL,NLF,NVD,NAN,LL4F,LL4W,LL4X,
     1   MAT,SIGT,SIDE,ZZ,FRZ,QFR,IPERT,KN,MUX,IPBX,LC,R,BX,TF,AX,XA)
         DEALLOCATE(FRZ,IPERT)
         IF((IPR.EQ.0).OR.(IPR.EQ.3)) THEN
            CALL LCMPPD(IPSYS,'X_'//TEXT10,IIMAX,2,XA_PTR)
         ELSE
            DEALLOCATE(XA)
         ENDIF
         CALL LCMPPD(IPSYS,'XA'//TEXT10,IIMAX,2,AX_PTR)
      ELSE
         IF(DIAG) THEN
            ALLOCATE(AX(IIMAX))
            IF(IPR.NE.3) THEN
               TF_PTR=LCMARA(LL4F*NLF/2)
               CALL C_F_POINTER(TF_PTR,TF,(/ LL4F*NLF/2 /))
               CALL XDRSET(TF,LL4F*NLF/2,0.0)
               CALL XDRSET(AX,IIMAX,0.0)
            ELSE
               IF(ISEG.GT.0) CALL MTBLD('XA'//TEXT10,IPTRK,IPSYS,1)
               CALL LCMGPD(IPSYS,'TF'//TEXT10,TF_PTR)
               CALL C_F_POINTER(TF_PTR,TF,(/ LL4F*NLF/2 /))
               CALL LCMGET(IPSYS,'XA'//TEXT10,AX)
            ENDIF
            ALLOCATE(XA(IIMAX))
         ELSE
            IF(IPR.NE.3) THEN
               TF_PTR=LCMARA(LL4F*NLF/2)
               AX_PTR=LCMARA(IIMAX)
               CALL C_F_POINTER(TF_PTR,TF,(/ LL4F*NLF/2 /))
               CALL C_F_POINTER(AX_PTR,AX,(/ IIMAX /))
               CALL XDRSET(TF,LL4F*NLF/2,0.0)
               CALL XDRSET(AX,IIMAX,0.0)
            ELSE
               IF(ISEG.GT.0) CALL MTBLD('XA'//TEXT10,IPTRK,IPSYS,1)
               CALL LCMGPD(IPSYS,'TF'//TEXT10,TF_PTR)
               CALL LCMGPD(IPSYS,'XA'//TEXT10,AX_PTR)
               CALL C_F_POINTER(TF_PTR,TF,(/ LL4F*NLF/2 /))
               CALL C_F_POINTER(AX_PTR,AX,(/ IIMAX /))
            ENDIF
            IF((IPR.EQ.0).OR.(IPR.EQ.3)) THEN
               XA_PTR=LCMARA(IIMAX)
               CALL C_F_POINTER(XA_PTR,XA,(/ IIMAX /))
            ELSE
               ALLOCATE(XA(IIMAX))
            ENDIF
         ENDIF
         CALL PN3DXX(NBMIX,IELEM,ICOL,NEL,NLF,NVD,NAN,LL4F,LL4X,SIGT,
     1   SIGTI,MAT,VOL2,XX,YY,ZZ,KN,QFR,MUX,IPBX,LC,R,V,BX,TF,AX,XA)
         IF(.NOT.DIAG) THEN
            IF((IPR.EQ.0).OR.(IPR.EQ.3)) THEN
               CALL LCMPPD(IPSYS,'X_'//TEXT10,IIMAX,2,XA_PTR)
            ELSE
               DEALLOCATE(XA)
            ENDIF
            CALL LCMPPD(IPSYS,'XA'//TEXT10,IIMAX,2,AX_PTR)
         ELSE
*           IN DIAGONAL SYMMETRY CASE, DO NOT SAVE THE X-DIRECTED
*           ADI MATRIX COMPONENT SINCE IT IS EQUAL TO THE Y-DIRECTED
*           COMPONENT
            DEALLOCATE(XA,AX)
         ENDIF
         CALL LCMPPD(IPSYS,'TF'//TEXT10,LL4F*NLF/2,2,TF_PTR)
      ENDIF
*
* DIMENSION Y
      IF(CHEX) THEN
         IIMAY=MUY(LL4Y)*NLF/2
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
         CALL C_F_POINTER(TF_PTR,TF,(/ LL4F*NLF/2 /))
         CALL LCMGPD(IPTRK,'IPBBY',IPBY_PTR)
         CALL LCMLEN(IPSYS,'YB',LENYB,ITYL)
         IF(LENYB.EQ.0) THEN
            CALL LCMGPD(IPTRK,'YB',BY_PTR)
         ELSE
            CALL LCMGPD(IPSYS,'YB',BY_PTR)
         ENDIF
         CALL C_F_POINTER(IPBY_PTR,IPBY,(/ 2*IELEM*LL4Y /))
         CALL C_F_POINTER(BY_PTR,BY,(/ 2*IELEM*LL4Y /))
         NBLOS=LX*LZ/3
         ALLOCATE(IPERT(NBLOS),FRZ(NBLOS))
         CALL LCMGET(IPTRK,'IPERT',IPERT)
         CALL LCMGET(IPTRK,'FRZ',FRZ)
         IF((IPR.EQ.0).OR.(IPR.EQ.3)) THEN
            YA_PTR=LCMARA(IIMAY)
            CALL C_F_POINTER(YA_PTR,YA,(/ IIMAY /))
         ELSE
            ALLOCATE(YA(IIMAY))
         ENDIF
         CALL PN3HWY(NBMIX,NBLOS,IELEM,ICOL,NLF,NVD,NAN,LL4F,LL4W,LL4X,
     1   LL4Y,MAT,SIGT,SIDE,ZZ,FRZ,QFR,IPERT,KN,MUY,IPBY,LC,R,BY,TF,AY,
     2   YA)
         DEALLOCATE(FRZ,IPERT)
         IF((IPR.EQ.0).OR.(IPR.EQ.3)) THEN
            CALL LCMPPD(IPSYS,'Y_'//TEXT10,IIMAY,2,YA_PTR)
         ELSE
            DEALLOCATE(YA)
         ENDIF
         CALL LCMPPD(IPSYS,'YA'//TEXT10,IIMAY,2,AY_PTR)
      ELSE IF(LL4Y.GT.0) THEN
         IIMAY=MUY(LL4Y)*NLF/2
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
         CALL C_F_POINTER(TF_PTR,TF,(/ LL4F*NLF/2 /))
         CALL LCMGPD(IPTRK,'IPBBY',IPBY_PTR)
         CALL LCMGPD(IPTRK,'YB',BY_PTR)
         CALL C_F_POINTER(IPBY_PTR,IPBY,(/ 2*IELEM*LL4Y /))
         CALL C_F_POINTER(BY_PTR,BY,(/ 2*IELEM*LL4Y /))
         IF((IPR.EQ.0).OR.(IPR.EQ.3)) THEN
            YA_PTR=LCMARA(IIMAY)
            CALL C_F_POINTER(YA_PTR,YA,(/ IIMAY /))
         ELSE
            ALLOCATE(YA(IIMAY))
         ENDIF
         CALL PN3DXY(NBMIX,IELEM,ICOL,NEL,NLF,NVD,NAN,LL4F,LL4X,LL4Y,
     1   SIGT,MAT,VOL2,YY,KN,QFR,MUY,IPBY,LC,R,BY,TF,AY,YA)
         IF((IPR.EQ.0).OR.(IPR.EQ.3)) THEN
            CALL LCMPPD(IPSYS,'Y_'//TEXT10,IIMAY,2,YA_PTR)
         ELSE
            DEALLOCATE(YA)
         ENDIF
         CALL LCMPPD(IPSYS,'YA'//TEXT10,IIMAY,2,AY_PTR)
      ENDIF
*
* DIMENSION Z
      IF(LL4Z.GT.0) THEN
         IIMAZ=MUZ(LL4Z)*NLF/2
         IF(CHEX) THEN
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
            CALL C_F_POINTER(TF_PTR,TF,(/ LL4F*NLF/2 /))
            CALL LCMGPD(IPTRK,'IPBBZ',IPBZ_PTR)
            CALL LCMLEN(IPSYS,'ZB',LENZB,ITYL)
            IF(LENZB.EQ.0) THEN
               CALL LCMGPD(IPTRK,'ZB',BZ_PTR)
            ELSE
               CALL LCMGPD(IPSYS,'ZB',BZ_PTR)
            ENDIF
            CALL C_F_POINTER(IPBZ_PTR,IPBZ,(/ 2*IELEM*LL4Z /))
            CALL C_F_POINTER(BZ_PTR,BZ,(/ 2*IELEM*LL4Z /))
            NBLOS=LX*LZ/3
            ALLOCATE(IPERT(NBLOS),FRZ(NBLOS))
            CALL LCMGET(IPTRK,'IPERT',IPERT)
            CALL LCMGET(IPTRK,'FRZ',FRZ)
            IF((IPR.EQ.0).OR.(IPR.EQ.3)) THEN
               ZA_PTR=LCMARA(IIMAZ)
               CALL C_F_POINTER(ZA_PTR,ZA,(/ IIMAZ /))
            ELSE
               ALLOCATE(ZA(IIMAZ))
            ENDIF
            CALL PN3HWZ(NBMIX,NBLOS,IELEM,ICOL,NLF,NVD,NAN,LL4F,LL4W,
     1      LL4X,LL4Y,LL4Z,MAT,SIGT,SIDE,ZZ,FRZ,QFR,IPERT,KN,MUZ,IPBZ,
     2      LC,R,BZ,TF,AZ,ZA)
            DEALLOCATE(FRZ,IPERT)
         ELSE
            IF(IPR.NE.3) THEN
               IF(IPR.EQ.0) THEN
                  AZ_PTR=LCMARA(IIMAZ)
                  CALL C_F_POINTER(AZ_PTR,AZ,(/ IIMAZ /))
               ELSE
                  ALLOCATE(ZA(IIMAZ))
               ENDIF
               CALL XDRSET(AZ,IIMAZ,0.0)
            ELSE            
               IF(ISEG.GT.0) CALL MTBLD('ZA'//TEXT10,IPTRK,IPSYS,1)
               CALL LCMGPD(IPSYS,'ZA'//TEXT10,AZ_PTR)
               CALL C_F_POINTER(AZ_PTR,AZ,(/ IIMAZ /))
            ENDIF
            CALL LCMGPD(IPSYS,'TF'//TEXT10,TF_PTR)
            CALL C_F_POINTER(TF_PTR,TF,(/ LL4F*NLF/2 /))
            CALL LCMGPD(IPTRK,'IPBBZ',IPBZ_PTR)
            CALL LCMGPD(IPTRK,'ZB',BZ_PTR)
            CALL C_F_POINTER(IPBZ_PTR,IPBZ,(/ 2*IELEM*LL4Z /))
            CALL C_F_POINTER(BZ_PTR,BZ,(/ 2*IELEM*LL4Z /))
            IF((IPR.EQ.0).OR.(IPR.EQ.3)) THEN
               ZA_PTR=LCMARA(IIMAZ)
               CALL C_F_POINTER(ZA_PTR,ZA,(/ IIMAZ /))
            ELSE
               ALLOCATE(ZA(IIMAZ))
            ENDIF
            CALL PN3DXZ(NBMIX,IELEM,ICOL,NEL,NLF,NVD,NAN,LL4F,LL4X,
     1      LL4Y,LL4Z,SIGT,MAT,VOL2,ZZ,KN,QFR,MUZ,IPBZ,LC,R,BZ,TF,
     2      AZ,ZA)
         ENDIF
         IF((IPR.EQ.0).OR.(IPR.EQ.3)) THEN
            CALL LCMPPD(IPSYS,'Z_'//TEXT10,IIMAZ,2,ZA_PTR)
         ELSE
            DEALLOCATE(ZA)
         ENDIF
         CALL LCMPPD(IPSYS,'ZA'//TEXT10,IIMAZ,2,AZ_PTR)
      ENDIF
      DEALLOCATE(VOL2)
      CALL KDRCPU(TK2)
      IF(IMPX.GT.1) WRITE(6,'(/35H TRIASN: CPU TIME FOR SYSTEM MATRIX,
     1 11H ASSEMBLY =,F9.2,3H S.)') TK2-TK1
*----
*  PERFORM SUPERVECTORIZATION REBUILD OF THE COEFFICIENT MATRICES
*----
      IF(ISEG.GT.0) THEN
         IF((IPR.EQ.0).OR.(IPR.EQ.3)) THEN
            IF(CHEX) CALL MTBLD('W_'//TEXT10,IPTRK,IPSYS,3)
            IF(.NOT.DIAG) CALL MTBLD('X_'//TEXT10,IPTRK,IPSYS,3)
            IF(LL4Y.GT.0) CALL MTBLD('Y_'//TEXT10,IPTRK,IPSYS,3)
            IF(LL4Z.GT.0) CALL MTBLD('Z_'//TEXT10,IPTRK,IPSYS,3)
         ENDIF
         IF(CHEX) CALL MTBLD('WA'//TEXT10,IPTRK,IPSYS,3)
         IF(.NOT.DIAG) CALL MTBLD('XA'//TEXT10,IPTRK,IPSYS,3)
         IF(LL4Y.GT.0) CALL MTBLD('YA'//TEXT10,IPTRK,IPSYS,3)
         IF(LL4Z.GT.0) CALL MTBLD('ZA'//TEXT10,IPTRK,IPSYS,3)
      ENDIF
*----
*  MATRIX FACTORIZATIONS
*----
      IF((IPR.EQ.0).OR.(IPR.EQ.3)) THEN
         CALL KDRCPU(TK1)
         CALL MTLDLF(TEXT10,IPTRK,IPSYS,ITY,IMPX)
         CALL KDRCPU(TK2)
         IF(IMPX.GT.1) WRITE(6,'(/34H TRIASN: CPU TIME FOR LDLT FACTORI,
     1   18HZATION OF MATRIX '',A10,2H''=,F9.2,3H S.)') TEXT10,TK2-TK1
      ENDIF
*----
*  RELEASE UNIT MATRICES
*----
      DEALLOCATE(V,R)
*----
*  RELEASE TRIVAC SPECIFIC TRACKING INFORMATION
*----
      DEALLOCATE(IQFR,QFR,KN,ZZ)
      IF(CHEX) THEN
         DEALLOCATE(MUW)
      ELSE
         DEALLOCATE(YY,XX)
      ENDIF
      IF(LL4Z.GT.0) DEALLOCATE(MUZ)
      IF(LL4Y.GT.0) DEALLOCATE(MUY)
      IF(.NOT.DIAG) DEALLOCATE(MUX)
      RETURN
      END
