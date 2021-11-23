*DECK MTLDLF
      SUBROUTINE MTLDLF(NAMP,IPTRK,IPSYS,ITY,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* LCM driver for the L-D-L(t) factorization of a symmetric matrix.
* The factorized matrix is stored on LCM under name 'I'//NAMP.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): Alain Hebert
*
*Parameters: input
* NAMP    name of the coefficient matrix.
* IPTRK   L_TRACK pointer to the tracking information.
* IPSYS   L_SYSTEM pointer to system matrices.
* ITY     type of coefficient matrix (1: Bivac; 2: classical Trivac;
*         3: Thomas-Raviart).
* IMPX    print flag (equal to zero for no print).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPSYS
      CHARACTER NAMP*12
      INTEGER ITY,IMPX
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40,NPREF=5)
      CHARACTER HIN*12,HOUT*12,PREFIX(5)*2,NAMLCM*12,NAMMY*12
      LOGICAL EMPTY,LCM
      INTEGER ITP(NSTATE)
      INTEGER, DIMENSION(:), ALLOCATABLE :: MU,NBL,LBL
      REAL, DIMENSION(:), ALLOCATABLE :: T
      REAL, DIMENSION(:), POINTER :: ASM
      TYPE(C_PTR) ASM_PTR
      DATA (PREFIX(I),I=1,NPREF)/'  ','W_','X_','Y_','Z_'/
*
      IF(ITY.EQ.1) THEN
*        BIVAC TRACKING.
         CALL LCMGET(IPTRK,'STATE-VECTOR',ITP)
         ISEG=0
         NLF=ITP(14)
      ELSE
*        CLASSICAL TRIVAC TRACKING.
         CALL LCMGET(IPTRK,'STATE-VECTOR',ITP)
         ISEG=ITP(17)
         NLF=ITP(30)
      ENDIF
*
      DO 30 IS=1,NPREF
      IF(PREFIX(IS).EQ.'  ') THEN
         HIN=NAMP
         HOUT='I'//NAMP(:11)
      ELSE
         HIN=PREFIX(IS)//NAMP(:11)
         HOUT=PREFIX(IS)(:1)//'I'//NAMP(:10)
      ENDIF
*----
*  PERFORM FACTORIZATION OF MATRICES
*----
      CALL LCMLEN(IPSYS,HIN,ILENG,ITYLCM)
      IF(ILENG.GT.0) THEN
         IF(ISEG.EQ.0) THEN
            CALL LCMLEN(IPTRK,'MU'//PREFIX(IS)(:1),LMU,ITYLCM)
            ALLOCATE(MU(LMU))
            CALL LCMGET(IPTRK,'MU'//PREFIX(IS)(:1),MU)
         ELSE
            CALL LCMLEN(IPTRK,'MUV'//PREFIX(IS)(:1),LMU,ITYLCM)
            ALLOCATE(MU(LMU))
            CALL LCMGET(IPTRK,'MUV'//PREFIX(IS)(:1),MU)
         ENDIF
         ILEN=MU(LMU)
         IF(NLF.GT.0) ILEN=ILEN*NLF/2
         IF(IMPX.GT.0) THEN
            CALL LCMINF(IPSYS,NAMLCM,NAMMY,EMPTY,ILONG,LCM)
            WRITE(6,'(/30H MTLDLF: FACTORIZATION OF LCM ,
     1      8HMATRIX '',A12,23H''. CREATION OF MATRIX '',A12,
     2      14H'' LOCATED IN '',A12,2H''.)') HIN,HOUT,NAMLCM
         ENDIF
         ASM_PTR=LCMARA(ILENG)
         CALL C_F_POINTER(ASM_PTR,ASM,(/ ILENG /))
         CALL LCMGET(IPSYS,HIN,ASM)
         IF(ISEG.EQ.0) THEN
            IF(ILEN.NE.ILENG) CALL XABORT('MTLDLF: INCONSISTENT INF'
     1      //'ORMATION ON LCM (1).')
            IF(NLF.EQ.0) THEN
               CALL ALLDLF(LMU,ASM,MU)
            ELSE
               IOF=1
               DO 10 IL=0,NLF-2,2
               CALL ALLDLF(LMU,ASM(IOF:),MU)
               IOF=IOF+MU(LMU)
   10          CONTINUE
            ENDIF
         ELSE
            IF(ISEG*ILEN.NE.ILENG) CALL XABORT('MTLDLF: INCONSISTEN'
     1      //'T INFORMATION ON LCM (2).')
            CALL LCMLEN(IPTRK,'NBL'//PREFIX(IS)(:1),LON,ITYLCM)
            ALLOCATE(NBL(LON),LBL(LON))
            CALL LCMGET(IPTRK,'NBL'//PREFIX(IS)(:1),NBL)
            CALL LCMGET(IPTRK,'LBL'//PREFIX(IS)(:1),LBL)
            ALLOCATE(T(ISEG))
            IF(NLF.EQ.0) THEN
               CALL ALVDLF(ASM,MU,ISEG,LON,NBL,LBL,T)
            ELSE
               IOF=1
               DO 20 IL=0,NLF-2,2
               CALL ALVDLF(ASM(IOF:),MU,ISEG,LON,NBL,LBL,T)
               IOF=IOF+MU(LMU)
   20          CONTINUE
            ENDIF
            DEALLOCATE(T,LBL,NBL)
         ENDIF
         DEALLOCATE(MU)
         CALL LCMPPD(IPSYS,HOUT,ILENG,2,ASM_PTR)
      ENDIF
   30 CONTINUE
      RETURN
      END
