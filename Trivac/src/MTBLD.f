*DECK MTBLD
      SUBROUTINE MTBLD(HNAME,IPTRK,IPSYS,ITY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* LCM driver for VECBLD.
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
* HNAME   name of the matrix. HNAME(:1) is 'W ', 'X ', 'Y', or 'Z'.
*         In case of a Thomas-Raviart basis, can also be equal to 'WA',
*        'XA', 'YA' or 'ZA'.
* IPTRK   L_TRACK pointer to the tracking information.
* IPSYS   L_SYSTEM pointer to system matrices.
* ITY     type of processing:
*         =1 gather back; =2 scatter forth;
*         =3 scatter forth and store the diagonal elements.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPSYS
      CHARACTER HNAME*(*)
      INTEGER ITY
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      CHARACTER TEXT12*12,HSMG*131,HCHAR*1
      INTEGER ITP(NSTATE)
      REAL DUMMY(1)
      TYPE(C_PTR) ASS_PTR,ASSV_PTR
      INTEGER, DIMENSION(:), ALLOCATABLE :: MU,MUV,IPV,LBL
      REAL, DIMENSION(:), POINTER :: ASS,ASSV
      REAL, DIMENSION(:), ALLOCATABLE :: DGV
*----
*  RECOVER TRACKING INFORMATION FROM LCM
*----
      CALL KDRCPU(TK1)
      HCHAR=HNAME(:1)
      CALL LCMGET(IPTRK,'STATE-VECTOR',ITP)
      ISEG=ITP(17)
      IMPV=ITP(18)
      CALL LCMLEN(IPTRK,'MU'//HCHAR,LL4,ITYLCM)
      CALL LCMLEN(IPTRK,'LBL'//HCHAR,LON,ITYLCM)
      ALLOCATE(MU(LL4),LBL(LON),MUV(LL4),IPV(LL4))
      CALL LCMGET(IPTRK,'MU'//HCHAR,MU)
      CALL LCMGET(IPTRK,'LBL'//HCHAR,LBL)
      CALL LCMGET(IPTRK,'MUV'//HCHAR,MUV)
      CALL LCMGET(IPTRK,'IPV'//HCHAR,IPV)
*
      TEXT12=HNAME
      IIMAX=MU(LL4)
      LBL0=0
      DO 10 I=1,LON
      LBL0=LBL0+LBL(I)
   10 CONTINUE
      IIMAXV=MUV(LBL0)*ISEG
      IF(ITY.EQ.1) THEN
*        SUPERVECTORIAL TO SCALAR REBUILD.
         ASS_PTR=LCMARA(IIMAX)
         CALL LCMLEN(IPSYS,TEXT12,ILEN,ITYLCM)
         IF(ILEN.NE.IIMAXV) THEN
            WRITE(HSMG,'(38HMTBLD: REBUILD FAILURE 1 IN PROCESSING,
     1      9H MATRIX '',A12,2H''.)') TEXT12
            CALL XABORT(HSMG)
         ENDIF
         CALL LCMGPD(IPSYS,TEXT12,ASSV_PTR)
         CALL C_F_POINTER(ASSV_PTR,ASSV,(/ IIMAXV /))
         CALL C_F_POINTER(ASS_PTR,ASS,(/ IIMAX /))
         CALL VECBLD(ISEG,LL4,MU,LON,LBL,MUV,IPV,1,ASS,ASSV,DUMMY)
         CALL LCMPPD(IPSYS,TEXT12,IIMAX,2,ASS_PTR)
      ELSE IF(ITY.GE.2) THEN
*        SCALAR TO SUPERVECTORIAL REBUILD.
         ALLOCATE(DGV(LBL0*ISEG))
         ASSV_PTR=LCMARA(IIMAXV)
         CALL LCMLEN(IPSYS,TEXT12,ILEN,ITYLCM)
         IF(ILEN.NE.IIMAX) THEN
            WRITE(HSMG,'(38HMTBLD: REBUILD FAILURE 2 IN PROCESSING,
     1      9H MATRIX '',A12,2H''.)') TEXT12
            CALL XABORT(HSMG)
         ENDIF
         CALL LCMGPD(IPSYS,TEXT12,ASS_PTR)
         CALL C_F_POINTER(ASSV_PTR,ASSV,(/ IIMAXV /))
         CALL C_F_POINTER(ASS_PTR,ASS,(/ IIMAX /))
         CALL VECBLD(ISEG,LL4,MU,LON,LBL,MUV,IPV,2,ASS,ASSV,DGV)
         IF(ITY.EQ.3) THEN
            CALL LCMPUT(IPSYS,HCHAR//'D'//TEXT12(3:),LBL0*ISEG,2,DGV)
         ENDIF
         DEALLOCATE(DGV)
         CALL LCMPPD(IPSYS,TEXT12,IIMAXV,2,ASSV_PTR)
      ENDIF
      DEALLOCATE(IPV,MUV,LBL,MU)
      CALL KDRCPU(TK2)
      IF(IMPV.GE.3) WRITE (6,'(/18H MTBLD: CPU TIME =,F7.2,3H S.)')
     1 TK2-TK1
      RETURN
      END
