*DECK MCGREC
      SUBROUTINE MCGREC(NFI,KM,MCUW,MCUI,MCU,LMCU,LMXMCU,IPRINT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Reform connection matrices.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): I. Suslov and R. Le Tellier
*
*Parameters: input
* NFI     total number of volumes and surfaces for which specific values
*         of the neutron flux and reactions rates are required.
* MCUW    undefined.
* MCUI    undefined.
* LMCU    dimension (used) of MCUW.
* LMXMCU  real dimension of MCUW MCUI.
* IPRINT  print level flag.
*
*Parameters: output
* KM      connection matrices.
* MCU     undefined.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NFI,KM(NFI),MCUW(LMXMCU),MCUI(LMXMCU),MCU(LMCU),LMCU,
     1 LMXMCU,IPRINT
*
      IM=0
      DO 10 I=1,NFI
      IP=IM+1
      II=0
      IC=I
      IF( MCUW(I).EQ.0 ) GOTO  9
    3 II=II+1
      IM=IM+1
      IF(IM.GT.LMCU) CALL XABORT('MCGREC: OVERFLOW.')
      MCU(IM)=MCUW(IC)
      IC=MCUI(IC)
      IF(IC.NE.0) GOTO 3
    9 CONTINUE
      KM(I)=II
      IF(II.EQ.0) GOTO 10
      IPP=IP+II-1
      IF(IPRINT.GE.10) WRITE (6,13) I,(MCU(JP),JP=IP,IPP)
      CALL SORTIN(II,MCU(IP))
   10 CONTINUE
*
   13 FORMAT(1X,'I=',I3,'  MCU=',30I4)
      RETURN
      END
