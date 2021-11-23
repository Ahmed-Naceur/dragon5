*DECK BIVASM
      SUBROUTINE BIVASM(HNAMT,ITY,IPTRK,IPSYS,IMPX,NBMIX,NEL,NLF,NDIM,
     1 NALBP,MAT,VOL,GAMMA,SGD)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembling of a single-group system matrix for BIVAC.
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal.
*
*Author(s): A. Hebert
*
*Parameters: input/output
* HNAMT   name of the matrix.
* ITY     type of assembly: =0: leakage-removal matrix assembly;
*         =1: cross section matrix assembly.
* IPTRK   L_TRACK pointer to the BIVAC tracking information.
* IPSYS   L_SYSTEM pointer to system matrices.
* IMPX    print parameter (equal to zero for no print).
* NBMIX   total number of material mixtures.
* NEL     total number of finite elements.
* NLF     number of Legendre orders for the flux (even number). Equal
*         to zero for diffusion theory.
* NDIM    second dimension of matrix SGD.
* NALBP   number of physical albedos.
* MAT     mixture index assigned to each volume.
* VOL     volume of each element.
* GAMMA   physical albedo functions.
* SGD     nuclear properties per material mixture.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER HNAMT*(*)
      TYPE(C_PTR) IPTRK,IPSYS
      INTEGER ITY,IMPX,NBMIX,NEL,NLF,NDIM,NALBP,MAT(NEL)
      REAL VOL(NEL),GAMMA(NALBP),SGD(NBMIX,NDIM)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      LOGICAL CYLIND
      CHARACTER TEXT11*11
      INTEGER ITP(NSTATE)
      INTEGER, DIMENSION(:), ALLOCATABLE :: KN,IQFR,MU,IPERT
      REAL, DIMENSION(:), ALLOCATABLE :: XX,YY,DD,QFR
      REAL, DIMENSION(:,:), ALLOCATABLE :: R,RS,Q,QS,V,H,RH,QH,RT,QT
      REAL, DIMENSION(:), POINTER :: SYS,ASS
      TYPE(C_PTR) SYS_PTR,ASS_PTR
*----
*  RECOVER BIVAC SPECIFIC TRACKING INFORMATION
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',ITP)
      ITYPE=ITP(6)
      CYLIND=(ITYPE.EQ.3).OR.(ITYPE.EQ.6)
      IELEM=ITP(8)
      ICOL=ITP(9)
      ISPLH=ITP(10)
      LL4=ITP(11)
      LX=ITP(12)
      LY=ITP(13)
      NVD=ITP(17)
      CALL LCMLEN(IPTRK,'KN',MAXKN,ITYLCM)
      CALL LCMLEN(IPTRK,'QFR',MAXQF,ITYLCM)
      ALLOCATE(XX(LX*LY),YY(LX*LY),DD(LX*LY),KN(MAXKN),QFR(MAXQF),
     1 IQFR(MAXQF),MU(LL4))
      IF(ITYPE.EQ.8) THEN
         CALL LCMGET(IPTRK,'SIDE',SIDE)
      ELSE
         CALL LCMGET(IPTRK,'XX',XX)
         CALL LCMGET(IPTRK,'YY',YY)
         CALL LCMGET(IPTRK,'DD',DD)
      ENDIF
      CALL LCMGET(IPTRK,'KN',KN)
      CALL LCMGET(IPTRK,'QFR',QFR)
      CALL LCMGET(IPTRK,'IQFR',IQFR)
      CALL LCMGET(IPTRK,'MU',MU)
*----
*  APPLY PHYSICAL ALBEDO FUNCTIONS
*----
      IF(NALBP.GT.0) THEN
         DO IQW=1,MAXQF
            IALB=IQFR(IQW)
            IF(IALB.NE.0) QFR(IQW)=QFR(IQW)*GAMMA(IALB)
         ENDDO
      ENDIF
*
      TEXT11=HNAMT
      IF(IMPX.GT.0) WRITE(6,'(/36H BIVASM: ASSEMBLY OF SYMMETRIC MATRI,
     1 3HX '',A11,38H'' IN COMPRESSED DIAGONAL STORAGE MODE.)') TEXT11
*----
*  ASSEMBLY OF THE SYSTEM MATRICES
*----
      CALL KDRCPU(TK1)
      IIMAX=MU(LL4)
      IF(NLF.NE.0) IIMAX=IIMAX*NLF/2
      SYS_PTR=LCMARA(IIMAX)
      CALL C_F_POINTER(SYS_PTR,SYS,(/ IIMAX /))
      CALL XDRSET(SYS,IIMAX,0.0)
*
      IF((IELEM.LT.0).AND.(ITYPE.NE.8)) THEN
*        MESH CORNER FINITE DIFFERENCES OR LAGRANGIAN FINITE ELEMENTS
         CALL LCMSIX(IPTRK,'BIVCOL',1)
         CALL LCMLEN(IPTRK,'T',LC,ITYLCM)
         ALLOCATE(R(LC,LC),RS(LC,LC),Q(LC,LC),QS(LC,LC))
         CALL LCMGET(IPTRK,'R',R)
         CALL LCMGET(IPTRK,'RS',RS)
         CALL LCMGET(IPTRK,'Q',Q)
         CALL LCMGET(IPTRK,'QS',QS)
         CALL LCMSIX(IPTRK,' ',2)
         CALL BIVA01(ITY,MAXKN,SGD,CYLIND,NEL,LL4,NBMIX,IIMAX,
     1   XX,YY,DD,MAT,KN,QFR,VOL,MU,LC,R,RS,Q,QS,SYS)
         DEALLOCATE(R,RS,Q,QS)
      ELSE IF((IELEM.GT.0).AND.(ITYPE.NE.8).AND.(NLF.GT.0)) THEN
*        MIXED-DUAL FINITE ELEMENTS (SIMPLIFIED PN THEORY)
         CALL LCMSIX(IPTRK,'BIVCOL',1)
         CALL LCMLEN(IPTRK,'T',LC,ITYLCM)
         ALLOCATE(R(LC,LC),V(LC,LC-1))
         CALL LCMGET(IPTRK,'R',R)
         CALL LCMGET(IPTRK,'V',V)
         CALL LCMSIX(IPTRK,' ',2)
         CALL PNDM2E(ITY,NEL,LL4,IELEM,ICOL,MAT,VOL,NBMIX,NLF,NVD,
     1   NDIM/2,SGD(1,1),SGD(1,1+NDIM/2),XX,YY,KN,QFR,MU,IIMAX,LC,
     2   R,V,SYS)
         DEALLOCATE(R,V)
      ELSE IF((IELEM.GT.0).AND.(ITYPE.NE.8)) THEN
*        MIXED-DUAL FINITE ELEMENTS (DIFFUSION THEORY).
         CALL LCMSIX(IPTRK,'BIVCOL',1)
         CALL LCMLEN(IPTRK,'T',LC,ITYLCM)
         ALLOCATE(R(LC,LC),V(LC,LC-1))
         CALL LCMGET(IPTRK,'R',R)
         CALL LCMGET(IPTRK,'V',V)
         CALL LCMSIX(IPTRK,' ',2)
         CALL BIVA02(ITY,SGD,CYLIND,IELEM,ICOL,NEL,LL4,NBMIX,IIMAX,XX,
     1   YY,DD,MAT,KN,QFR,VOL,MU,LC,R,V,SYS)
         DEALLOCATE(R,V)
      ELSE IF((IELEM.LT.0).AND.(ITYPE.EQ.8)) THEN
*        MESH CORNER FINITE DIFFERENCES FOR HEXAGONS
         CALL LCMSIX(IPTRK,'BIVCOL',1)
         ALLOCATE(R(2,2),RH(6,6),QH(6,6),RT(3,3),QT(3,3))
         CALL LCMGET(IPTRK,'R',R)
         CALL LCMGET(IPTRK,'RH',RH)
         CALL LCMGET(IPTRK,'QH',QH)
         CALL LCMGET(IPTRK,'RT',RT)
         CALL LCMGET(IPTRK,'QT',QT)
         CALL LCMSIX(IPTRK,' ',2)
         IF(ISPLH.EQ.1) THEN
            NELEM=MAXKN/7
         ELSE
            NELEM=MAXKN/4
         ENDIF
         CALL BIVA03(ITY,MAXKN,MAXQF,SGD,NEL,LL4,ISPLH,NELEM,NBMIX,
     1   IIMAX,SIDE,MAT,KN,QFR,VOL,MU,R,RH,QH,RT,QT,SYS)
         DEALLOCATE(R,RH,QH,RT,QT)
      ELSE IF((IELEM.GT.0).AND.(ITYPE.EQ.8).AND.(ICOL.EQ.4)) THEN
*        MESH CENTERED FINITE DIFFERENCES FOR HEXAGONS
         CALL BIVA04(ITY,MAXKN,MAXQF,SGD,NEL,LL4,ISPLH,NBMIX,IIMAX,
     1   SIDE,MAT,KN,QFR,VOL,MU,SYS)
      ELSE IF((IELEM.GT.0).AND.(ITYPE.EQ.8).AND.(NLF.GT.0)) THEN
*        THOMAS-RAVIART-SCHNEIDER METHOD FOR HEXAGONS (SIMPLIFIED PN
*        THEORY)
         LXH=LX/(3*ISPLH**2)
         NBLOS=LXH*ISPLH**2
         ALLOCATE(IPERT(NBLOS))
         CALL LCMGET(IPTRK,'IPERT',IPERT)
         CALL LCMSIX(IPTRK,'BIVCOL',1)
         CALL LCMLEN(IPTRK,'T',LC,ITYLCM)
         ALLOCATE(R(LC,LC),V(LC,LC-1),H(LC,LC-1))
         CALL LCMGET(IPTRK,'R',R)
         CALL LCMGET(IPTRK,'V',V)
         CALL LCMGET(IPTRK,'H',H)
         CALL LCMSIX(IPTRK,' ',2)
         CALL PNDH2E(ITY,IELEM,ICOL,NBLOS,LL4,NBMIX,IIMAX,SIDE,MAT,
     1   IPERT,SGD(1,1),KN,QFR,NLF,NVD,NDIM/2,MU,LC,R,V,H,SYS)
         DEALLOCATE(R,V,H,IPERT)
      ELSE IF((IELEM.GT.0).AND.(ITYPE.EQ.8)) THEN
*        THOMAS-RAVIART-SCHNEIDER METHOD FOR HEXAGONS (DIFFUSION THEORY)
         LXH=LX/(3*ISPLH**2)
         NBLOS=LXH*ISPLH**2
         ALLOCATE(IPERT(NBLOS))
         CALL LCMGET(IPTRK,'IPERT',IPERT)
         CALL LCMSIX(IPTRK,'BIVCOL',1)
         CALL LCMLEN(IPTRK,'T',LC,ITYLCM)
         ALLOCATE(R(LC,LC),V(LC,LC-1),H(LC,LC-1))
         CALL LCMGET(IPTRK,'R',R)
         CALL LCMGET(IPTRK,'V',V)
         CALL LCMGET(IPTRK,'H',H)
         CALL LCMSIX(IPTRK,' ',2)
         CALL BIVA05(ITY,SGD,IELEM,NBLOS,LL4,NBMIX,IIMAX,SIDE,MAT,IPERT,
     1   KN,QFR,MU,LC,R,V,H,SYS)
         DEALLOCATE(R,V,H,IPERT)
      ENDIF
      CALL LCMPPD(IPSYS,TEXT11,IIMAX,2,SYS_PTR)
      CALL KDRCPU(TK2)
      IF(IMPX.GT.0) WRITE(6,'(/35H BIVASM: CPU TIME FOR SYSTEM MATRIX,
     1 11H ASSEMBLY =,F9.2,3H S.)') TK2-TK1
*----
*  MATRIX FACTORIZATIONS
*----
      IF((ITY.EQ.0).OR.(TEXT11.EQ.'RM')) THEN
         CALL KDRCPU(TK1)
         ASS_PTR=LCMARA(IIMAX)
         CALL C_F_POINTER(ASS_PTR,ASS,(/ IIMAX /))
         CALL LCMGET(IPSYS,TEXT11,ASS)
         IF(NLF.EQ.0) THEN
            CALL ALLDLF(LL4,ASS,MU)
         ELSE
            IOF=1
            DO 50 IL=0,NLF-2,2
            CALL ALLDLF(LL4,ASS(IOF:),MU)
            IOF=IOF+MU(LL4)
   50       CONTINUE
         ENDIF
         CALL LCMPPD(IPSYS,'I'//TEXT11,IIMAX,2,ASS_PTR)
         CALL KDRCPU(TK2)
         IF(IMPX.GT.1) WRITE(6,'(/34H BIVASM: CPU TIME FOR LDLT FACTORI,
     1   18HZATION OF MATRIX '',A11,2H''=,F9.2,3H S.)') TEXT11,TK2-TK1
      ENDIF
*----
*  RELEASE BIVAC SPECIFIC TRACKING INFORMATION
*----
      DEALLOCATE(MU,IQFR,QFR,KN,DD,XX,YY)
      RETURN
      END
