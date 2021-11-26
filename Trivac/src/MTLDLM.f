*DECK MTLDLM
      SUBROUTINE MTLDLM(NAMP,IPTRK,IPSYS,LL4,ITY,F2,F3)
*
*-----------------------------------------------------------------------
*
*Purpose:
* LCM driver for the multiplication of a matrix by a vector.
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
* F2      vector to multiply.
*
*Parameters: output
* F3      result of the multiplication.
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
      REAL F2(LL4),F3(LL4)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      CHARACTER NAMT*12,TEXT12*12
      INTEGER ITP(NSTATE),ASS_LEN
      LOGICAL LMU,LMUW,LMUX,LMUY,LMUZ,DIAG
      REAL, DIMENSION(:), ALLOCATABLE :: GAR,GAF
      TYPE(C_PTR) MU_PTR,IP_PTR,IPV_PTR,NBL_PTR,LBL_PTR
      TYPE(C_PTR) ASS_PTR
      INTEGER, DIMENSION(:), ALLOCATABLE :: IPB
      INTEGER, DIMENSION(:), POINTER :: MU,IP,IPV,NBL,LBL
      REAL, DIMENSION(:), POINTER :: ASS
*
*-----------------------------------------------------------------------
*
* INFORMATION RECOVERED FROM XSM OR LCM (SPLITTED MATRIX):
*  'W_'//NAMP 'X_'//NAMP 'Y_'//NAMP 'Z_'//NAMP : W-, X-, Y- AND Z-
*                            ORIENTED MATRIX COMPONENTS.
*
* SCALAR INFORMATION RECOVERED FROM LCM.
*  'MUW' 'MUX' 'MUY' 'MUZ' : POSITION OF DIAGONAL ELEMENT IN W-, X-, Y-
*                            OR Z-ORIENTED MATRIX COMPONENTS.
*  'IPW' 'IPX' 'IPY' 'IPZ' : PERMUTATION INFORMATION FOR W-, X-, Y- OR
*                            Z-ORIENTED MATRIX COMPONENTS.
*
* SUPERVECTORIZATION INFORMATION RECOVERED FROM LCM.
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
*        DIFFUSION BIVAC TRACKING.
         ISEG=0
      ELSE IF(ITY.EQ.2) THEN
*        CLASSICAL TRIVAC TRACKING.
         CALL LCMGET(IPTRK,'STATE-VECTOR',ITP)
         ISEG=ITP(17)
         LTSW=ITP(19)
      ELSE IF(ITY.EQ.3) THEN
*        RAVIART-THOMAS TRIVAC TRACKING.
         CALL FLDTRM(NAMP,IPTRK,IPSYS,LL4,F2,F3)
         RETURN
      ELSE IF(ITY.EQ.11) THEN
*        SIMPLIFIED PN BIVAC TRACKING.
         CALL LCMGET(IPSYS,'STATE-VECTOR',ITP)
         NBMIX=ITP(7)
         NAN=ITP(8)
         IF(NAN.EQ.0) CALL XABORT('MTLDLM: SPN-ONLY ALGORITHM(1).')
         CALL XDRSET(F3,LL4,0.0)
         CALL FLDBSM(NAMP,IPTRK,IPSYS,LL4,NBMIX,NAN,F2,F3)
         RETURN
      ELSE IF(ITY.EQ.13) THEN
*        SIMPLIFIED PN TRIVAC TRACKING.
         CALL LCMGET(IPSYS,'STATE-VECTOR',ITP)
         NBMIX=ITP(7)
         NAN=ITP(8)
         IF(NAN.EQ.0) CALL XABORT('MTLDLM: SPN-ONLY ALGORITHM(2).')
         CALL XDRSET(F3,LL4,0.0)
         CALL FLDTSM(NAMP,IPTRK,IPSYS,LL4,NBMIX,NAN,F2,F3)
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
         IF(LL4.NE.LL4TS) CALL XABORT('MTLDLM: INVALID LL4(1).')
         CALL LCMGPD(IPTRK,'MU',MU_PTR)
         CALL LCMGPD(IPSYS,NAMT,ASS_PTR)
         CALL C_F_POINTER(MU_PTR,MU,(/ LL4 /))
         CALL C_F_POINTER(ASS_PTR,ASS,(/ MU(LL4) /))
         CALL ALLDLM(LL4,ASS,F2,F3,MU,1)
      ELSE IF(ISEG.EQ.0) THEN
*        SCALAR MULTIPLICATION FOR A W- OR X-ORIENTED MATRIX.
         IF(LMUW) THEN
            TEXT12='W_'//NAMT(:10)
            CALL LCMGPD(IPTRK,'MUW',MU_PTR)
            CALL LCMGPD(IPTRK,'IPW',IP_PTR)
            CALL LCMLEN(IPTRK,'IPW',LL4TS,ITYLCM)
         ELSE IF(DIAG) THEN
            TEXT12='Y_'//NAMT(:10)
            CALL LCMGPD(IPTRK,'MUY',MU_PTR)
            CALL LCMGPD(IPTRK,'IPX',IP_PTR)
            CALL LCMLEN(IPTRK,'IPX',LL4TS,ITYLCM)
         ELSE
            TEXT12='X_'//NAMT(:10)
            CALL LCMGPD(IPTRK,'MUX',MU_PTR)
            CALL LCMGPD(IPTRK,'IPX',IP_PTR)
            CALL LCMLEN(IPTRK,'IPX',LL4TS,ITYLCM)
         ENDIF
         IF(LL4.NE.LL4TS) CALL XABORT('MTLDLM: INVALID LL4(2).')
         CALL C_F_POINTER(MU_PTR,MU,(/ LL4 /))
         CALL C_F_POINTER(IP_PTR,IP,(/ LL4 /))
         ALLOCATE(GAR(LL4))
         DO 10 I=1,LL4
         GAR(IP(I))=F2(I)
   10    CONTINUE
         DO 20 I=1,LL4
         F2(I)=GAR(I)
   20    CONTINUE
         CALL LCMGPD(IPSYS,TEXT12,ASS_PTR)
         CALL C_F_POINTER(ASS_PTR,ASS,(/ MU(LL4) /))
         CALL ALLDLM(LL4,ASS,F2,GAR,MU,1)
         DO 30 I=1,LL4
         II=IP(I)
         F3(I)=GAR(II)
         GAR(II)=F2(II)
   30    CONTINUE
         DO 40 I=1,LL4
         F2(I)=GAR(IP(I))
   40    CONTINUE
         IF(LMUW) THEN
*           SCALAR MULTIPLICATION FOR A X-ORIENTED MATRIX.
            CALL LCMGPD(IPTRK,'MUX',MU_PTR)
            CALL LCMLEN(IPTRK,'IPX',LL4TS,ITYLCM)
            IF(LL4.NE.LL4TS) CALL XABORT('MTLDLM: INVALID LL4(5).')
            CALL LCMGPD(IPTRK,'IPX',IP_PTR)
            CALL C_F_POINTER(MU_PTR,MU,(/ LL4 /))
            CALL C_F_POINTER(IP_PTR,IP,(/ LL4 /))
            DO 50 I=1,LL4
            GAR(IP(I))=F2(I)
   50       CONTINUE
            DO 60 I=1,LL4
            II=IP(I)
            F2(II)=GAR(II)
            GAR(II)=F3(I)
   60       CONTINUE
            CALL LCMGPD(IPSYS,'X_'//NAMT,ASS_PTR)
            CALL C_F_POINTER(ASS_PTR,ASS,(/ MU(LL4) /))
            CALL ALLDLM(LL4,ASS,F2,GAR,MU,2)
            DO 70 I=1,LL4
            II=IP(I)
            F3(I)=GAR(II)
            GAR(II)=F2(II)
   70       CONTINUE
            DO 80 I=1,LL4
            F2(I)=GAR(IP(I))
   80       CONTINUE
         ENDIF
         IF(LMUY) THEN
*           SCALAR MULTIPLICATION FOR A Y-ORIENTED MATRIX.
            CALL LCMGPD(IPTRK,'MUY',MU_PTR)
            CALL LCMLEN(IPTRK,'IPY',LL4TS,ITYLCM)
            IF(LL4.NE.LL4TS) CALL XABORT('MTLDLM: INVALID LL4(6).')
            CALL LCMGPD(IPTRK,'IPY',IP_PTR)
            CALL C_F_POINTER(MU_PTR,MU,(/ LL4 /))
            CALL C_F_POINTER(IP_PTR,IP,(/ LL4 /))
            DO 90 I=1,LL4
            GAR(IP(I))=F2(I)
   90       CONTINUE
            DO 100 I=1,LL4
            II=IP(I)
            F2(II)=GAR(II)
            GAR(II)=F3(I)
  100       CONTINUE
            CALL LCMGPD(IPSYS,'Y_'//NAMT,ASS_PTR)
            CALL C_F_POINTER(ASS_PTR,ASS,(/ MU(LL4) /))
            CALL ALLDLM(LL4,ASS,F2,GAR,MU,2)
            DO 110 I=1,LL4
            II=IP(I)
            F3(I)=GAR(II)
            GAR(II)=F2(II)
  110       CONTINUE
            DO 120 I=1,LL4
            F2(I)=GAR(IP(I))
  120       CONTINUE
         ENDIF
         IF(LMUZ) THEN
*           SCALAR MULTIPLICATION FOR A Z-ORIENTED MATRIX.
            CALL LCMGPD(IPTRK,'MUZ',MU_PTR)
            CALL LCMLEN(IPTRK,'IPZ',LL4TS,ITYLCM)
            IF(LL4.NE.LL4TS) CALL XABORT('MTLDLM: INVALID LL4(7).')
            CALL LCMGPD(IPTRK,'IPZ',IP_PTR)
            CALL C_F_POINTER(MU_PTR,MU,(/ LL4 /))
            CALL C_F_POINTER(IP_PTR,IP,(/ LL4 /))
            DO 130 I=1,LL4
            GAR(IP(I))=F2(I)
  130       CONTINUE
            DO 140 I=1,LL4
            II=IP(I)
            F2(II)=GAR(II)
            GAR(II)=F3(I)
  140       CONTINUE
            CALL LCMGPD(IPSYS,'Z_'//NAMT,ASS_PTR)
            CALL C_F_POINTER(ASS_PTR,ASS,(/ MU(LL4) /))
            CALL ALLDLM(LL4,ASS,F2,GAR,MU,2)
            DO 150 I=1,LL4
            II=IP(I)
            F3(I)=GAR(II)
            GAR(II)=F2(II)
  150       CONTINUE
            DO 160 I=1,LL4
            F2(I)=GAR(IP(I))
  160       CONTINUE
         ENDIF
         DEALLOCATE(GAR)
      ELSE IF(ISEG.GT.0) THEN
*        SUPERVECTORIAL MULTIPLICATION FOR A W- OR X-ORIENTED MATRIX.
         IF(LMUW) THEN
            CALL LCMGET(IPTRK,'LL4VW',LL4V)
            CALL LCMGPD(IPTRK,'MUVW',MU_PTR)
            TEXT12='W_'//NAMT(:10)
            CALL LCMLEN(IPTRK,'IPW',LL4TS,ITYLCM)
            IF(LL4.NE.LL4TS) CALL XABORT('MTLDLM: INVALID LL4(8).')
            CALL LCMGPD(IPTRK,'IPW',IP_PTR)
            CALL LCMGPD(IPTRK,'IPVW',IPV_PTR)
            CALL LCMLEN(IPTRK,'NBLW',NBL_LEN,ITYLCM)
            CALL LCMGPD(IPTRK,'NBLW',NBL_PTR)
            CALL LCMGPD(IPTRK,'LBLW',LBL_PTR)
         ELSE IF(DIAG) THEN
            CALL LCMGET(IPTRK,'LL4VY',LL4V)
            CALL LCMGPD(IPTRK,'MUVY',MU_PTR)
            TEXT12='Y_'//NAMT(:10)
            CALL LCMLEN(IPTRK,'IPX',LL4TS,ITYLCM)
            IF(LL4.NE.LL4TS) CALL XABORT('MTLDLM: INVALID LL4(9).')
            CALL LCMGPD(IPTRK,'IPX',IP_PTR)
            CALL LCMGPD(IPTRK,'IPVY',IPV_PTR)
            CALL LCMLEN(IPTRK,'NBLY',NBL_LEN,ITYLCM)
            CALL LCMGPD(IPTRK,'NBLY',NBL_PTR)
            CALL LCMGPD(IPTRK,'LBLY',LBL_PTR)
         ELSE
            CALL LCMGET(IPTRK,'LL4VX',LL4V)
            CALL LCMGPD(IPTRK,'MUVX',MU_PTR)
            TEXT12='X_'//NAMT(:10)
            CALL LCMLEN(IPTRK,'IPX',LL4TS,ITYLCM)
            IF(LL4.NE.LL4TS) CALL XABORT('MTLDLM: INVALID LL4(10).')
            CALL LCMGPD(IPTRK,'IPX',IP_PTR)
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
         ALLOCATE(IPB(LL4),GAR(LL4V),GAF(LL4V))
         DO 165 I=1,LL4
         IPB(I)=IPV(IP(I))
  165    CONTINUE
         CALL XDRSET(GAR,LL4V,0.0)
         DO 180 I=1,LL4
         GAR(IPB(I))=F2(I)
  180    CONTINUE
         CALL LCMLEN(IPSYS,TEXT12,ASS_LEN,ITYLCM)
         CALL LCMGPD(IPSYS,TEXT12,ASS_PTR)
         CALL C_F_POINTER(ASS_PTR,ASS,(/ ASS_LEN /))
         CALL ALVDLM(LTSW,ASS,GAR,GAF,MU,1,ISEG,NBL_LEN,NBL,LBL)
         DO 190 I=1,LL4
         II=IPB(I)
         F2(I)=GAR(II)
         F3(I)=GAF(II)
  190    CONTINUE
         DEALLOCATE(GAF,GAR,IPB)
         IF(LMUW) THEN
*           SUPERVECTORIAL MULTIPLICATION FOR A X-ORIENTED MATRIX.
            CALL LCMGET(IPTRK,'LL4VX',LL4V)
            CALL LCMGPD(IPTRK,'MUVX',MU_PTR)
            CALL C_F_POINTER(MU_PTR,MU,(/ LL4V/ISEG /))
            CALL LCMLEN(IPTRK,'IPX',LL4TS,ITYLCM)
            IF(LL4.NE.LL4TS) CALL XABORT('MTLDLM: INVALID LL4(11).')
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
            ALLOCATE(IPB(LL4),GAR(LL4V),GAF(LL4V))
            DO 200 I=1,LL4
            IPB(I)=IPV(IP(I))
  200       CONTINUE
            CALL XDRSET(GAR,LL4V,0.0)
            CALL XDRSET(GAF,LL4V,0.0)
            DO 220 I=1,LL4
            II=IPB(I)
            GAR(II)=F2(I)
            GAF(II)=F3(I)
  220       CONTINUE
            CALL LCMLEN(IPSYS,'X_'//NAMT,ASS_LEN,ITYLCM)
            CALL LCMGPD(IPSYS,'X_'//NAMT,ASS_PTR)
            CALL C_F_POINTER(ASS_PTR,ASS,(/ ASS_LEN /))
            CALL ALVDLM(LTSW,ASS,GAR,GAF,MU,2,ISEG,NBL_LEN,NBL,LBL)
            DO 230 I=1,LL4
            II=IPB(I)
            F2(I)=GAR(II)
            F3(I)=GAF(II)
  230       CONTINUE
            DEALLOCATE(GAF,GAR,IPB)
         ENDIF
         IF(LMUY) THEN
*           SUPERVECTORIAL MULTIPLICATION FOR A Y-ORIENTED MATRIX.
            CALL LCMGET(IPTRK,'LL4VY',LL4V)
            CALL LCMGPD(IPTRK,'MUVY',MU_PTR)
            CALL C_F_POINTER(MU_PTR,MU,(/ LL4V/ISEG /))
            CALL LCMLEN(IPTRK,'IPY',LL4TS,ITYLCM)
            IF(LL4.NE.LL4TS) CALL XABORT('MTLDLM: INVALID LL4(12).')
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
            ALLOCATE(IPB(LL4),GAR(LL4V),GAF(LL4V))
            DO 235 I=1,LL4
            IPB(I)=IPV(IP(I))
  235       CONTINUE
            CALL XDRSET(GAR,LL4V,0.0)
            CALL XDRSET(GAF,LL4V,0.0)
            DO 260 I=1,LL4
            II=IPB(I)
            GAR(II)=F2(I)
            GAF(II)=F3(I)
  260       CONTINUE
            CALL LCMLEN(IPSYS,'Y_'//NAMT,ASS_LEN,ITYLCM)
            CALL LCMGPD(IPSYS,'Y_'//NAMT,ASS_PTR)
            CALL C_F_POINTER(ASS_PTR,ASS,(/ ASS_LEN /))
            CALL ALVDLM(LTSW,ASS,GAR,GAF,MU,2,ISEG,NBL_LEN,NBL,LBL)
            DO 270 I=1,LL4
            II=IPB(I)
            F2(I)=GAR(II)
            F3(I)=GAF(II)
  270       CONTINUE
            DEALLOCATE(GAF,GAR,IPB)
         ENDIF
         IF(LMUZ) THEN
*           SUPERVECTORIAL MULTIPLICATION FOR A Z-ORIENTED MATRIX.
            CALL LCMGET(IPTRK,'LL4VZ',LL4V)
            CALL LCMGPD(IPTRK,'MUVZ',MU_PTR)
            CALL C_F_POINTER(MU_PTR,MU,(/ LL4V/ISEG /))
            CALL LCMLEN(IPTRK,'IPZ',LL4TS,ITYLCM)
            IF(LL4.NE.LL4TS) CALL XABORT('MTLDLM: INVALID LL4(13).')
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
            ALLOCATE(IPB(LL4),GAR(LL4V),GAF(LL4V))
            DO 275 I=1,LL4
            IPB(I)=IPV(IP(I))
  275       CONTINUE
            CALL XDRSET(GAR,LL4V,0.0)
            CALL XDRSET(GAF,LL4V,0.0)
            DO 300 I=1,LL4
            II=IPB(I)
            GAR(II)=F2(I)
            GAF(II)=F3(I)
  300       CONTINUE
            CALL LCMLEN(IPSYS,'Z_'//NAMT,ASS_LEN,ITYLCM)
            CALL LCMGPD(IPSYS,'Z_'//NAMT,ASS_PTR)
            CALL C_F_POINTER(ASS_PTR,ASS,(/ ASS_LEN /))
            CALL ALVDLM(LTSW,ASS,GAR,GAF,MU,2,ISEG,NBL_LEN,NBL,LBL)
            DO 310 I=1,LL4
            II=IPB(I)
            F2(I)=GAR(II)
            F3(I)=GAF(II)
  310       CONTINUE
            DEALLOCATE(GAF,GAR,IPB)
         ENDIF
      ENDIF
      RETURN
      END
