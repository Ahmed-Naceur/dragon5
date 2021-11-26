*DECK FMAC02
      SUBROUTINE FMAC02(IPMACR,NK,NGROUP,ARRAY,HNAME)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Save a 1D cross section in the GROUP list of a MACROLIB.
*
*Copyright:
* Copyright (C) 2020 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPMACR  LCM object address of the MACROLIB.
* NK      number of mixtures.
* NGROUP  number of energy groups.
* ARRAY   array to save.
* HNAME   MACROLIB name of the cross section.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR)  IPMACR
      INTEGER      NK,NGROUP
      REAL         ARRAY(NK,NGROUP)
      CHARACTER(LEN=*) HNAME
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR)  JPMACR,KPMACR
      CHARACTER HSMG*131
*
      CALL LCMLEN(IPMACR,'GROUP',ILENG,ITYLCM)
      ITY=0
      IF(ILENG.EQ.NGROUP) THEN
        ITY=1
      ELSE IF(ILENG+1.EQ.NGROUP) THEN
        ITY=2
      ELSE
        CALL LCMLIB(IPMACR)
        WRITE(HSMG,'(33HFMAC02: INVALID VALUE OF NGROUP (,I6,5H) XS=,
     1  A8,1H.)') NGROUP,HNAME
        CALL XABORT(HSMG)
      ENDIF
      JPMACR=LCMGID(IPMACR,'GROUP')
      DO IG=1,ILENG
        KPMACR=LCMDIL(JPMACR,IG)
        IF(ITY.EQ.1) THEN
          CALL LCMPUT(KPMACR,HNAME,NK,2,ARRAY(:NK,IG))
        ELSE IF(ITY.EQ.2) THEN
          CALL LCMPUT(KPMACR,HNAME,NK*2,2,ARRAY(:NK,IG:IG+1))
        ENDIF
      ENDDO
      RETURN
      END
