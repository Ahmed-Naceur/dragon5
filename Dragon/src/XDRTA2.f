*DECK XDRTA2
      SUBROUTINE XDRTA2(IPTRK)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover the tabulated functions required by the flux solution.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy and A. Hebert
*
*Parameters: input
* IPTRK   pointer to the tracking information (L_TRACK signature).
*
*-----------------------------------------------------------------------
*
***** OUTPUT:   THE FIVE COMMONS OF BICKLEY QUADRATIC TABLES ARE FILLED
*               AND SAVED WITH RESPECTIVE NAMES:
*                  /BICKL1/,/BICKL2/,/BICKL3/,/BICKL4/,/BICKL5/
*
*               A COMMON FOR LINEAR EXPONENTIAL TABLES IS FILLED
*               AND SAVED WITH NAME: /EXP1/
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK
*----
*  LOCAL VARIABLES
*----
      REAL PASV(5),XLIMV(5),PARAM(3)
      INTEGER MLOG(5)
*
      PARAMETER (MKI1=600,MKI2=600,MKI3=600,MKI4=600,MKI5=600)
      COMMON /BICKL1/BI1(0:MKI1),BI11(0:MKI1),BI12(0:MKI1)
     >              ,PAS1,XLIM1,L1
      COMMON /BICKL2/BI2(0:MKI2),BI21(0:MKI2),BI22(0:MKI2)
     >              ,PAS2,XLIM2,L2
      COMMON /BICKL3/BI3(0:MKI3),BI31(0:MKI3),BI32(0:MKI3)
     >              ,PAS3,XLIM3,L3
      COMMON /BICKL4/BI4(0:MKI4),BI41(0:MKI4),BI42(0:MKI4)
     >              ,PAS4,XLIM4,L4
      COMMON /BICKL5/BI5(0:MKI5),BI51(0:MKI5),BI52(0:MKI5)
     >              ,PAS5,XLIM5,L5
      SAVE  /BICKL1/,/BICKL2/,/BICKL3/,/BICKL4/,/BICKL5/
*
      PARAMETER (MEX1=7936)
      COMMON /EXP1/ E0(0:MEX1),E1(0:MEX1),PASE1,DXE1,XLIME1
      COMMON /EXP0/ E00(0:MEX1),E01(0:MEX1),PASE0,DXE0,XLIME0
      SAVE   /EXP1/,/EXP0/
*----
*  CHARGE BICKLEY TABLES INTO COMMON
*----
      CALL LCMLEN(IPTRK,'FUNC-TABLES',ILENG,ITYLCM)
      IF(ILENG.EQ.0) RETURN
      CALL LCMSIX(IPTRK,'FUNC-TABLES',1)
      CALL LCMLEN(IPTRK,'PAS',ILENG,ITYLCM)
      IF(ILENG.GT.0) THEN
         CALL LCMGET(IPTRK,'PAS',PASV)
         CALL LCMGET(IPTRK,'XLIM',XLIMV)
         CALL LCMGET(IPTRK,'MLOG',MLOG)
         CALL LCMLEN(IPTRK,'BI1',ILENG2,ITYLCM)
         IF(ILENG2.NE.MKI1+1) CALL XABORT('XDRTA2: INVALID BICKLEY.')
         CALL LCMGET(IPTRK,'BI1',BI1(0))
         CALL LCMGET(IPTRK,'BI11',BI11(0))
         CALL LCMGET(IPTRK,'BI12',BI12(0))
         CALL LCMGET(IPTRK,'BI2',BI2(0))
         CALL LCMGET(IPTRK,'BI21',BI21(0))
         CALL LCMGET(IPTRK,'BI22',BI22(0))
         CALL LCMGET(IPTRK,'BI3',BI3(0))
         CALL LCMGET(IPTRK,'BI31',BI31(0))
         CALL LCMGET(IPTRK,'BI32',BI32(0))
         CALL LCMGET(IPTRK,'BI4',BI4(0))
         CALL LCMGET(IPTRK,'BI41',BI41(0))
         CALL LCMGET(IPTRK,'BI42',BI42(0))
         CALL LCMGET(IPTRK,'BI5',BI5(0))
         CALL LCMGET(IPTRK,'BI51',BI51(0))
         CALL LCMGET(IPTRK,'BI52',BI52(0))
         PAS1=PASV(1)
         PAS2=PASV(2)
         PAS3=PASV(3)
         PAS4=PASV(4)
         PAS5=PASV(5)
         XLIM1=XLIMV(1)
         XLIM2=XLIMV(2)
         XLIM3=XLIMV(3)
         XLIM4=XLIMV(4)
         XLIM5=XLIMV(5)
         L1=MLOG(1)
         L2=MLOG(2)
         L3=MLOG(3)
         L4=MLOG(4)
         L5=MLOG(5)
      ENDIF
*----
*  CHARGE EXPONENTIAL TABLES INTO COMMON
*----
      CALL LCMLEN(IPTRK,'PARAM',ILENG,ITYLCM)
      IF(ILENG.GT.0) THEN
         CALL LCMGET(IPTRK,'PARAM',PARAM)
         CALL LCMLEN(IPTRK,'E10',ILENG2,ITYLCM)
         IF(ILENG2.NE.MEX1+1) CALL XABORT('XDRTA2: INVALID EXP.')
         CALL LCMGET(IPTRK,'E10',E0(0))
         CALL LCMGET(IPTRK,'E11',E1(0))
         CALL LCMGET(IPTRK,'E00',E00(0))
         CALL LCMGET(IPTRK,'E01',E01(0))
         PASE1=PARAM(1)
         DXE1=PARAM(2)
         XLIME1=PARAM(3)
         PASE0=PARAM(1)
         DXE0=PARAM(2)
         XLIME0=PARAM(3)
      ENDIF
      CALL LCMSIX(IPTRK,' ',2)
      RETURN
      END
