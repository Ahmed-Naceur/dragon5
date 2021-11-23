*DECK EDIHMX
      SUBROUTINE EDIHMX(IPTRK,NREGIO,NMERGE,IMERGE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Find merge vector for merge by HMIX.
*
*Copyright:
* Copyright (C) 2001 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau
*
*Parameters: input
* IPTRK   calculation tracking data structure.
* NREGIO  number of regions.
*
*Parameters: output
* NMERGE  final number of merged regions.
* IMERGE  merged region positions.
*
*-----------------------------------------------------------------------
*
      USE         GANLIB
      IMPLICIT    NONE
      INTEGER     IOUT,NSTATE
      CHARACTER   NAMSBR*6
      PARAMETER  (IOUT=6,NSTATE=40,NAMSBR='EDIHMX')
*----
*  ROUTINE PARAMETERS
*----
      TYPE(C_PTR) IPTRK
      INTEGER     NREGIO
      INTEGER     NMERGE
      INTEGER     IMERGE(NREGIO)
*----
*  LOCAL PARAMETERS
*----
      INTEGER     IMRGLN,IMRGTY,IREG
      INTEGER     ISTATE(NSTATE)
*----
*  IMERGE is HOMMATCOD or MATCOD
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      CALL LCMLEN(IPTRK,'HOMMATCOD   ',IMRGLN,IMRGTY)
      IF(IMRGLN .EQ. 0) THEN
        WRITE(IOUT,8000) NAMSBR
        CALL LCMGET(IPTRK,'MATCOD      ',IMERGE)
      ELSE
        CALL LCMGET(IPTRK,'HOMMATCOD   ',IMERGE)
      ENDIF
*----
*  Check for double heterogeneity
*----
      IF(ISTATE(40).EQ.1) THEN
        CALL EDIBHX (NREGIO,IPTRK,IMRGLN,IMERGE)
      ENDIF
      IF(IMRGLN.NE.NREGIO) CALL XABORT('EDIHMX: bad nb of regions')
*----
*  Compute number of merged regions
*----
      NMERGE=0
      DO IREG=1,NREGIO
        NMERGE=MAX(NMERGE,IMERGE(IREG))
      ENDDO
      RETURN
*----
*  WARNING FORMAT
*----
 8000 FORMAT('***** Warning in routine - ',A6,' - *****'/
     >'No HMIX data in GEO: for NXT tracking file'/
     >'Homogenize using MIX instead of HMIX')
      END
