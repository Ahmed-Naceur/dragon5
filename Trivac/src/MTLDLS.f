*DECK MTLDLS
      SUBROUTINE MTLDLS(NAMP,IPTRK,IPSYS,LL4,ITY,F1)
*
*-----------------------------------------------------------------------
*
*Purpose:
* LCM driver for the solution of a linear system after LDL(t)
* factorization.
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
* NAMP    name of the coefficient matrix.
* IPTRK   L_TRACK pointer to the tracking information.
* IPSYS   L_SYSTEM pointer to system matrices.
* LL4     order of the matrix.
* ITY     type of coefficient matrix (1: Bivac; 2: classical Trivac;
*         3: Raviart-Thomas; 11: SPN/Bivac; 13: SPN/Raviart-Thomas).
* F1      right-hand side of the linear system.
*
*Parameters: output
* F1      solution of the linear system.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPSYS
      CHARACTER NAMP*(*)
      INTEGER LL4,ITY
      REAL F1(LL4)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      CHARACTER NAMT*12,TEXT12*12
      INTEGER ITP(NSTATE),ASS_LEN
      LOGICAL LMU,LMUW,LMUX,LMUY,LMUZ,DIAG
      REAL, DIMENSION(:), ALLOCATABLE :: GAR
      TYPE(C_PTR) MU_PTR,IP_PTR,IPV_PTR,NBL_PTR,LBL_PTR
      TYPE(C_PTR) ASS_PTR,DGV_PTR
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPB
      INTEGER, DIMENSION(:), POINTER :: MU,IP,IPV,NBL,LBL
      REAL, DIMENSION(:), POINTER :: ASS,DGV
*
*-----------------------------------------------------------------------
*
* INFORMATION RECOVERED FROM XSM OR LCM (NON-SPLITTED MATRIX):
*  NAMP       : COEFFICIENT MATRIX.
*  'I'//NAMP  : FACTORIZED COEFFICIENT MATRIX.
*  'MU'       : POSITION OF DIAGONAL ELEMENT IN COEFFICIENT MATRIX.
*
* INFORMATION RECOVERED FROM XSM OR LCM (SPLITTED MATRIX):
*  'W_'//NAMP 'X_'//NAMP 'Y_'//NAMP 'Z_'//NAMP : W-, X-, Y- AND Z-
*                            ORIENTED MATRIX COMPONENTS.
*  'WI'//NAMP 'XI'//NAMP 'YI'//NAMP 'ZI'//NAMP : W-, X-, Y- AND Z-
*                            ORIENTED FACTORIZED MATRIX COMPONENTS.
*
* SCALAR INFORMATION RECOVERED FROM LCM.
*  'MUW' 'MUX' 'MUY' 'MUZ' : POSITION OF DIAGONAL ELEMENT IN W-, X-, Y-
*                            OR Z-ORIENTED MATRIX COMPONENTS.
*  'IPW' 'IPX' 'IPY' 'IPZ' : PERMUTATION INFORMATION FOR W-, X-, Y- OR
*                            Z-ORIENTED MATRIX COMPONENTS.
*
* SUPERVECTORIZATION INFORMATION RECOVERED FROM LCM.
*  'WD'//NAMP 'XD'//NAMP 'YD'//NAMP 'ZD'//NAMP : DIAGONAL ELEMENTS FOR
*                           W-, X-, Y- AND Z-ORIENTED MATRIX COMPONENTS.
*  'LL4VW' 'LL4VX' 'LL4VY' 'LL4VZ' : ORDER OF THE REORDERED W-, X-, Y-
*                            AND Z-ORIENTED MATRIX COMPONENTS.
*  'MUVW' 'MUVX' 'MUVY' 'MUVZ' : POSITION OF DIAGONAL ELEMENT IN W-, X-,
*                            Y-OR Z-ORIENTED MATRIX COMPONENTS.
*  'IPVW' 'IPVX' 'IPVY' 'IPVZ' : PERMUTATION INFORMATION FOR W-, X-, Y-
*                            OR Z-ORIENTED MATRIX COMPONENTS.
*  'NBLW' 'NBLX' 'NBLY' 'NBLZ' : NUMBER OF LINEAR SYSTEMS IN EACH SUPER-
*                            VECTORIAL UNKNOWN GROUP.
*  'LBLW' 'LBLX' 'LBLY' 'LBLZ' : ORDER OF LINEAR SYSTEMS IN EACH SUPER-
*                            VECTORIAL UNKNOWN GROUP.
*
*-----------------------------------------------------------------------
*
      IF(ITY.EQ.1) THEN
*        BIVAC TRACKING.
         ISEG=0
      ELSE IF(ITY.EQ.2) THEN
*        CLASSICAL TRIVAC TRACKING.
         CALL LCMGET(IPTRK,'STATE-VECTOR',ITP)
         ISEG=ITP(17)
         LTSW=ITP(19)
      ELSE IF(ITY.EQ.3) THEN
*        RAVIART-THOMAS/DIFFUSION TRIVAC TRACKING.
         ALLOCATE(GAR(LL4))
         GAR(:LL4)=F1(:LL4)
         F1(:LL4)=0.0
         CALL FLDTRS(NAMP,IPTRK,IPSYS,LL4,GAR,F1,1)
         DEALLOCATE(GAR)
         RETURN
      ELSE IF(ITY.EQ.11) THEN
*        SIMPLIFIED PN BIVAC TRACKING.
         CALL LCMGET(IPSYS,'STATE-VECTOR',ITP)
         NBMIX=ITP(7)
         NAN=ITP(8)
         IF(NAN.EQ.0) CALL XABORT('MTLDLS: SPN-ONLY ALGORITHM(1).')
         CALL FLDBSS(NAMP,IPTRK,IPSYS,LL4,NBMIX,NAN,F1,1)
         RETURN
      ELSE IF(ITY.EQ.13) THEN
*        RAVIART-THOMAS/SIMPLIFIED PN TRIVAC TRACKING.
         CALL LCMGET(IPSYS,'STATE-VECTOR',ITP)
         NBMIX=ITP(7)
         NAN=ITP(8)
         IF(NAN.EQ.0) CALL XABORT('MTLDLS: SPN-ONLY ALGORITHM(2).')
         ALLOCATE(GAR(LL4))
         GAR(:LL4)=F1(:LL4)
         F1(:LL4)=0.0
         CALL FLDSPN(NAMP,IPTRK,IPSYS,LL4,NBMIX,NAN,GAR,F1,1)
         DEALLOCATE(GAR)
         RETURN
      ENDIF
*
      CALL LCMLEN(IPTRK,'MU',IDUM,ITYLCM)
      LMU=(IDUM.NE.0).AND.(ITYLCM.EQ.1)
      CALL LCMLEN(IPTRK,'MUW',IDUM,ITYLCM)
      LMUW=(IDUM.NE.0).AND.(ITYLCM.EQ.1)
      CALL LCMLEN(IPTRK,'MUX',IDUM,ITYLCM)
      LMUX=(IDUM.NE.0).AND.(ITYLCM.EQ.1)
      CALL LCMLEN(IPTRK,'MUY',IDUM,ITYLCM)
      LMUY=(IDUM.NE.0).AND.(ITYLCM.EQ.1)
      CALL LCMLEN(IPTRK,'MUZ',IDUM,ITYLCM)
      LMUZ=(IDUM.NE.0).AND.(ITYLCM.EQ.1)
      DIAG=LMUY.AND.(.NOT.LMUX)
*
      NAMT=NAMP
      IF(LMU) THEN
         CALL LCMLEN(IPTRK,'MU',LL4TS,ITYLCM)
         IF(LL4.NE.LL4TS) CALL XABORT('MTLDLS: INVALID LL4(1).')
         CALL LCMGPD(IPTRK,'MU',MU_PTR)
         CALL LCMGPD(IPSYS,'I'//NAMT,ASS_PTR)
         CALL C_F_POINTER(MU_PTR,MU,(/ LL4 /))
         CALL C_F_POINTER(ASS_PTR,ASS,(/ MU(LL4) /))
         CALL ALLDLS(LL4,MU,ASS,F1)
      ELSE IF(ISEG.EQ.0) THEN
*        SCALAR SOLUTION FOR A W- OR X-ORIENTED LINEAR SYSTEM.
         TEXT12=' '
         IF(LMUW) THEN
            TEXT12='WI'//NAMT(:10)
            CALL LCMGPD(IPTRK,'MUW',MU_PTR)
            CALL LCMGPD(IPTRK,'IPW',IP_PTR)
            CALL LCMLEN(IPTRK,'IPW',LL4TS,ITYLCM)
         ELSE IF(DIAG) THEN
            TEXT12='YI'//NAMT(:10)
            CALL LCMGPD(IPTRK,'MUY',MU_PTR)
            CALL LCMGPD(IPTRK,'IPX',IP_PTR)
            CALL LCMLEN(IPTRK,'IPX',LL4TS,ITYLCM)
         ELSE
            TEXT12='XI'//NAMT(:10)
            CALL LCMGPD(IPTRK,'MUX',MU_PTR)
            CALL LCMGPD(IPTRK,'IPX',IP_PTR)
            CALL LCMLEN(IPTRK,'IPX',LL4TS,ITYLCM)
         ENDIF
         IF(LL4.NE.LL4TS) CALL XABORT('MTLDLS: INVALID LL4(2).')
         CALL C_F_POINTER(MU_PTR,MU,(/ LL4 /))
         CALL C_F_POINTER(IP_PTR,IP,(/ LL4 /))
         ALLOCATE(GAR(LL4))
         DO 10 I=1,LL4
         GAR(IP(I))=F1(I)
   10    CONTINUE
         CALL LCMGPD(IPSYS,TEXT12,ASS_PTR)
         CALL C_F_POINTER(ASS_PTR,ASS,(/ MU(LL4) /))
         CALL ALLDLS(LL4,MU,ASS,GAR)
         DO 20 I=1,LL4
         F1(I)=GAR(IP(I))
   20    CONTINUE
         IF(LMUW) THEN
*           SCALAR SOLUTION FOR A X-ORIENTED LINEAR SYSTEM.
            CALL LCMGPD(IPTRK,'MUX',MU_PTR)
            CALL LCMGPD(IPTRK,'IPX',IP_PTR)
            CALL C_F_POINTER(MU_PTR,MU,(/ LL4 /))
            CALL C_F_POINTER(IP_PTR,IP,(/ LL4 /))
            CALL LCMGPD(IPSYS,'X_'//NAMT,ASS_PTR)
            CALL C_F_POINTER(ASS_PTR,ASS,(/ MU(LL4) /))
            DO 30 I=1,LL4
            II=IP(I)
            GAR(II)=F1(I)*ASS(MU(II))
   30       CONTINUE
            CALL LCMGPD(IPSYS,'XI'//NAMT,ASS_PTR)
            CALL C_F_POINTER(ASS_PTR,ASS,(/ MU(LL4) /))
            CALL ALLDLS(LL4,MU,ASS,GAR)
            DO 50 I=1,LL4
            F1(I)=GAR(IP(I))
   50       CONTINUE
         ENDIF
         IF(LMUY) THEN
*           SCALAR SOLUTION FOR A Y-ORIENTED LINEAR SYSTEM.
            CALL LCMGPD(IPTRK,'MUY',MU_PTR)
            CALL LCMGPD(IPTRK,'IPY',IP_PTR)
            CALL C_F_POINTER(MU_PTR,MU,(/ LL4 /))
            CALL C_F_POINTER(IP_PTR,IP,(/ LL4 /))
            CALL LCMGPD(IPSYS,'Y_'//NAMT,ASS_PTR)
            CALL C_F_POINTER(ASS_PTR,ASS,(/ MU(LL4) /))
            DO 60 I=1,LL4
            II=IP(I)
            GAR(II)=F1(I)*ASS(MU(II))
   60       CONTINUE
            CALL LCMGPD(IPSYS,'YI'//NAMT,ASS_PTR)
            CALL C_F_POINTER(ASS_PTR,ASS,(/ MU(LL4) /))
            CALL ALLDLS(LL4,MU,ASS,GAR)
            DO 80 I=1,LL4
            F1(I)=GAR(IP(I))
   80       CONTINUE
         ENDIF
         IF(LMUZ) THEN
*           SCALAR SOLUTION FOR A Z-ORIENTED LINEAR SYSTEM.
            CALL LCMGPD(IPTRK,'MUZ',MU_PTR)
            CALL LCMGPD(IPTRK,'IPZ',IP_PTR)
            CALL C_F_POINTER(MU_PTR,MU,(/ LL4 /))
            CALL C_F_POINTER(IP_PTR,IP,(/ LL4 /))
            CALL LCMGPD(IPSYS,'Z_'//NAMT,ASS_PTR)
            CALL C_F_POINTER(ASS_PTR,ASS,(/ MU(LL4) /))
            DO 90 I=1,LL4
            II=IP(I)
            GAR(II)=F1(I)*ASS(MU(II))
   90       CONTINUE
            CALL LCMGPD(IPSYS,'ZI'//NAMT,ASS_PTR)
            CALL C_F_POINTER(ASS_PTR,ASS,(/ MU(LL4) /))
            CALL ALLDLS(LL4,MU,ASS,GAR)
            DO 110 I=1,LL4
            F1(I)=GAR(IP(I))
  110       CONTINUE
         ENDIF
         DEALLOCATE(GAR)
      ELSE IF(ISEG.GT.0) THEN
*        SUPERVECTORIAL SOLUTION FOR A W- OR X-ORIENTED LINEAR SYSTEM.
         IF(LMUW) THEN
            CALL LCMGET(IPTRK,'LL4VW',LL4V)
            CALL LCMGPD(IPTRK,'MUVW',MU_PTR)
            TEXT12='WI'//NAMT(:10)
            CALL LCMLEN(IPTRK,'IPW',LL4TS,ITYLCM)
            IF(LL4.NE.LL4TS) CALL XABORT('MTLDLS: INVALID LL4(5).')
            CALL LCMGPD(IPTRK,'IPW',IP_PTR)
            CALL LCMLEN(IPTRK,'IPVW',IPV_LEN,ITYLCM)
            CALL LCMGPD(IPTRK,'IPVW',IPV_PTR)
            CALL LCMLEN(IPTRK,'NBLW',NBL_LEN,ITYLCM)
            CALL LCMGPD(IPTRK,'NBLW',NBL_PTR)
            CALL LCMGPD(IPTRK,'LBLW',LBL_PTR)
         ELSE IF(DIAG) THEN
            CALL LCMGET(IPTRK,'LL4VY',LL4V)
            CALL LCMGPD(IPTRK,'MUVY',MU_PTR)
            TEXT12='YI'//NAMT(:10)
            CALL LCMLEN(IPTRK,'IPX',LL4TS,ITYLCM)
            IF(LL4.NE.LL4TS) CALL XABORT('MTLDLS: INVALID LL4(6).')
            CALL LCMGPD(IPTRK,'IPX',IP_PTR)
            CALL LCMLEN(IPTRK,'IPVY',IPV_LEN,ITYLCM)
            CALL LCMGPD(IPTRK,'IPVY',IPV_PTR)
            CALL LCMLEN(IPTRK,'NBLY',NBL_LEN,ITYLCM)
            CALL LCMGPD(IPTRK,'NBLY',NBL_PTR)
            CALL LCMGPD(IPTRK,'LBLY',LBL_PTR)
         ELSE
            CALL LCMGET(IPTRK,'LL4VX',LL4V)
            CALL LCMGPD(IPTRK,'MUVX',MU_PTR)
            TEXT12='XI'//NAMT(:10)
            CALL LCMLEN(IPTRK,'IPX',LL4TS,ITYLCM)
            IF(LL4.NE.LL4TS) CALL XABORT('MTLDLS: INVALID LL4(7).')
            CALL LCMGPD(IPTRK,'IPX',IP_PTR)
            CALL LCMLEN(IPTRK,'IPVX',IPV_LEN,ITYLCM)
            CALL LCMGPD(IPTRK,'IPVX',IPV_PTR)
            CALL LCMLEN(IPTRK,'NBLX',NBL_LEN,ITYLCM)
            CALL LCMGPD(IPTRK,'NBLX',NBL_PTR)
            CALL LCMGPD(IPTRK,'LBLX',LBL_PTR)
         ENDIF
         CALL C_F_POINTER(MU_PTR,MU,(/ LL4V/ISEG /))
         CALL C_F_POINTER(IP_PTR,IP,(/ LL4 /))
         CALL C_F_POINTER(IPV_PTR,IPV,(/ IPV_LEN /))
         CALL C_F_POINTER(NBL_PTR,NBL,(/ NBL_LEN /))
         CALL C_F_POINTER(LBL_PTR,LBL,(/ NBL_LEN /))
         ALLOCATE(IPB(LL4),GAR(LL4V))
         DO 120 I=1,LL4
         IPB(I)=IPV(IP(I))
  120    CONTINUE
         CALL XDRSET(GAR,LL4V,0.0)
         DO 130 I=1,LL4
         GAR(IPB(I))=F1(I)
  130    CONTINUE
         CALL LCMLEN(IPSYS,TEXT12,ASS_LEN,ITYLCM)
         CALL LCMGPD(IPSYS,TEXT12,ASS_PTR)
         CALL C_F_POINTER(ASS_PTR,ASS,(/ ASS_LEN /))
         CALL ALVDLS(LTSW,MU,ASS,GAR,ISEG,NBL_LEN,NBL,LBL)
         DO 140 I=1,LL4
         F1(I)=GAR(IPB(I))
  140    CONTINUE
         DEALLOCATE(GAR,IPB)
         IF(LMUW) THEN
*           SUPERVECTORIAL SOLUTION FOR A X-ORIENTED LINEAR SYSTEM.
            CALL LCMGET(IPTRK,'LL4VX',LL4V)
            CALL LCMGPD(IPTRK,'MUVX',MU_PTR)
            CALL C_F_POINTER(MU_PTR,MU,(/ LL4V/ISEG /))
            CALL LCMLEN(IPTRK,'IPX',LL4,ITYLCM)
            CALL LCMGPD(IPTRK,'IPX',IP_PTR)
            CALL LCMLEN(IPTRK,'IPVX',IPV_LEN,ITYLCM)
            CALL LCMGPD(IPTRK,'IPVX',IPV_PTR)
            CALL C_F_POINTER(IP_PTR,IP,(/ LL4 /))
            CALL C_F_POINTER(IPV_PTR,IPV,(/ IPV_LEN /))
            CALL LCMLEN(IPTRK,'NBLX',NBL_LEN,ITYLCM)
            CALL LCMGPD(IPTRK,'NBLX',NBL_PTR)
            CALL LCMGPD(IPTRK,'LBLX',LBL_PTR)
            CALL C_F_POINTER(NBL_PTR,NBL,(/ NBL_LEN /))
            CALL C_F_POINTER(LBL_PTR,LBL,(/ NBL_LEN /))
            ALLOCATE(IPB(LL4),GAR(LL4V))
            DO 150 I=1,LL4
            IPB(I)=IPV(IP(I))
  150       CONTINUE
            CALL XDRSET(GAR,LL4V,0.0)
            DO 160 I=1,LL4
            GAR(IPB(I))=F1(I)
  160       CONTINUE
            CALL LCMLEN(IPSYS,'XI'//NAMT,ASS_LEN,ITYLCM)
            CALL LCMGPD(IPSYS,'XI'//NAMT,ASS_PTR)
            CALL LCMGPD(IPSYS,'XD'//NAMT,DGV_PTR)
            CALL C_F_POINTER(ASS_PTR,ASS,(/ ASS_LEN /))
            CALL C_F_POINTER(DGV_PTR,DGV,(/ LL4V /))
CDIR$ IVDEP
            DO 170 I=1,LL4V
            GAR(I)=GAR(I)*DGV(I)
  170       CONTINUE
            CALL ALVDLS(LTSW,MU,ASS,GAR,ISEG,NBL_LEN,NBL,LBL)
            DO 190 I=1,LL4
            F1(I)=GAR(IPB(I))
  190       CONTINUE
            DEALLOCATE(GAR,IPB)
         ENDIF
         IF(LMUY) THEN
*           SUPERVECTORIAL SOLUTION FOR A Y-ORIENTED LINEAR SYSTEM.
            CALL LCMGET(IPTRK,'LL4VY',LL4V)
            CALL LCMGPD(IPTRK,'MUVY',MU_PTR)
            CALL C_F_POINTER(MU_PTR,MU,(/ LL4V/ISEG /))
            CALL LCMLEN(IPTRK,'IPY',LL4,ITYLCM)
            CALL LCMGPD(IPTRK,'IPY',IP_PTR)
            CALL LCMLEN(IPTRK,'IPVY',IPV_LEN,ITYLCM)
            CALL LCMGPD(IPTRK,'IPVY',IPV_PTR)
            CALL C_F_POINTER(IP_PTR,IP,(/ LL4 /))
            CALL C_F_POINTER(IPV_PTR,IPV,(/ IPV_LEN /))
            CALL LCMLEN(IPTRK,'NBLY',NBL_LEN,ITYLCM)
            CALL LCMGPD(IPTRK,'NBLY',NBL_PTR)
            CALL LCMGPD(IPTRK,'LBLY',LBL_PTR)
            CALL C_F_POINTER(NBL_PTR,NBL,(/ NBL_LEN /))
            CALL C_F_POINTER(LBL_PTR,LBL,(/ NBL_LEN /))
            ALLOCATE(IPB(LL4),GAR(LL4V))
            DO 200 I=1,LL4
            IPB(I)=IPV(IP(I))
  200       CONTINUE
            CALL XDRSET(GAR,LL4V,0.0)
            DO 210 I=1,LL4
            GAR(IPB(I))=F1(I)
  210       CONTINUE
            CALL LCMLEN(IPSYS,'YI'//NAMT,ASS_LEN,ITYLCM)
            CALL LCMGPD(IPSYS,'YI'//NAMT,ASS_PTR)
            CALL LCMGPD(IPSYS,'YD'//NAMT,DGV_PTR)
            CALL C_F_POINTER(ASS_PTR,ASS,(/ ASS_LEN /))
            CALL C_F_POINTER(DGV_PTR,DGV,(/ LL4V /))
CDIR$ IVDEP
            DO 220 I=1,LL4V
            GAR(I)=GAR(I)*DGV(I)
  220       CONTINUE
            CALL ALVDLS(LTSW,MU,ASS,GAR,ISEG,NBL_LEN,NBL,LBL)
            DO 240 I=1,LL4
            F1(I)=GAR(IPB(I))
  240       CONTINUE
            DEALLOCATE(GAR,IPB)
         ENDIF
         IF(LMUZ) THEN
*           SUPERVECTORIAL SOLUTION FOR A Z-ORIENTED LINEAR SYSTEM.
            CALL LCMGET(IPTRK,'LL4VZ',LL4V)
            CALL LCMGPD(IPTRK,'MUVZ',MU_PTR)
            CALL C_F_POINTER(MU_PTR,MU,(/ LL4V/ISEG /))
            CALL LCMLEN(IPTRK,'IPZ',LL4,ITYLCM)
            CALL LCMGPD(IPTRK,'IPZ',IP_PTR)
            CALL LCMLEN(IPTRK,'IPVZ',IPV_LEN,ITYLCM)
            CALL LCMGPD(IPTRK,'IPVZ',IPV_PTR)
            CALL C_F_POINTER(IP_PTR,IP,(/ LL4 /))
            CALL C_F_POINTER(IPV_PTR,IPV,(/ IPV_LEN /))
            CALL LCMLEN(IPTRK,'NBLZ',NBL_LEN,ITYLCM)
            CALL LCMGPD(IPTRK,'NBLZ',NBL_PTR)
            CALL LCMGPD(IPTRK,'LBLZ',LBL_PTR)
            CALL C_F_POINTER(NBL_PTR,NBL,(/ NBL_LEN /))
            CALL C_F_POINTER(LBL_PTR,LBL,(/ NBL_LEN /))
            ALLOCATE(IPB(LL4),GAR(LL4V))
            DO 250 I=1,LL4
            IPB(I)=IPV(IP(I))
  250       CONTINUE
            CALL XDRSET(GAR,LL4V,0.0)
            DO 260 I=1,LL4
            GAR(IPB(I))=F1(I)
  260       CONTINUE
            CALL LCMLEN(IPSYS,'ZI'//NAMT,ASS_LEN,ITYLCM)
            CALL LCMGPD(IPSYS,'ZI'//NAMT,ASS_PTR)
            CALL LCMGPD(IPSYS,'ZD'//NAMT,DGV_PTR)
            CALL C_F_POINTER(ASS_PTR,ASS,(/ ASS_LEN /))
            CALL C_F_POINTER(DGV_PTR,DGV,(/ LL4V /))
CDIR$ IVDEP
            DO 270 I=1,LL4V
            GAR(I)=GAR(I)*DGV(I)
  270       CONTINUE
            CALL ALVDLS(LTSW,MU,ASS,GAR,ISEG,NBL_LEN,NBL,LBL)
            DO 290 I=1,LL4
            F1(I)=GAR(IPB(I))
  290       CONTINUE
            DEALLOCATE(GAR,IPB)
         ENDIF
      ENDIF
      RETURN
      END
