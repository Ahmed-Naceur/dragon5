*DECK COMCAL
      SUBROUTINE COMCAL(IMPX,IPCPO,IPDEPL,IPEDIT,IPEDI2,LMACRO,LISO,
     1 ITRES)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Store the results of an elementary calculation in the multicompo.
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
* IMPX    print parameter.
* IPCPO   pointer to the multicompo.
* IPDEPL  pointer to the burnup object (L_BURNUP signature).
* IPEDIT  pointer to the edition object (L_EDIT signature).
* IPEDI2  pointer to the edition object containing group form factor
*         information (L_EDIT signature).
* LMACRO  flag set to .TRUE. to recover cross sections from the
*         macrolib.
* LISO    =.true. if we want to register the region number of the
*         isotopes.
*
*Parameters: output
* ITRES   creation index for the macroscopic residual (=0: not created;
*         =1: not a FP precursor; =2: is a FP precursor).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IMPX,ITRES
      TYPE(C_PTR) IPCPO,IPDEPL,IPEDIT,IPEDI2
      LOGICAL LMACRO,LISO
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40,MAXISO=100)
      TYPE(C_PTR) JPCPO,KPCPO,LPCPO
      INTEGER ISTATE(NSTATE),IPAR(NSTATE)
      REAL BIRRAD(2)
      CHARACTER CDIRO*12,HSMG*131,NOMISP(MAXISO)*8
*
      CALL LCMGTC(IPEDIT,'LAST-EDIT',12,1,CDIRO)
      CALL LCMSIX(IPEDIT,CDIRO,1)
      IF(LMACRO) THEN
         CALL LCMSIX(IPEDIT,'MACROLIB',1)
         CALL LCMGET(IPEDIT,'STATE-VECTOR',IPAR)
         NMIL=IPAR(2)
         NISOTS=1
         NG=IPAR(1)
         NED=IPAR(5)
         NW=IPAR(10)
         CALL LCMSIX(IPEDIT,' ',2)
      ELSE
         CALL LCMGET(IPEDIT,'STATE-VECTOR',IPAR)
         NMIL=IPAR(1)
         NISOTS=IPAR(2)
         NG=IPAR(3)
         NED=IPAR(13)
         NW=IPAR(25)
      ENDIF
      CALL LCMSIX(IPEDIT,' ',2)
*
      CALL LCMGET(IPCPO,'STATE-VECTOR',ISTATE)
      IF(ISTATE(3).EQ.0) THEN
*        COMPLETE STATE-VECTOR.
         IF(ISTATE(1).EQ.0) THEN
            ISTATE(1)=NMIL
         ELSE IF(NMIL.NE.ISTATE(1)) THEN
            WRITE(HSMG,'(42HCOMCAL: ELEMENTARY CALCULATION WITH AN INV,
     1      22HALIB NB. OF MIXTURES =,I7,3H NE,I7,1H.)') NMIL,ISTATE(1)
            CALL XABORT(HSMG)
         ENDIF
         ISTATE(2)=NG
      ELSE
         IF(NMIL.NE.ISTATE(1)) THEN
            WRITE(HSMG,'(42HCOMCAL: ELEMENTARY CALCULATION WITH AN INV,
     1      22HALIB NB. OF MIXTURES =,I7,3H NE,I7,1H.)') NMIL,ISTATE(1)
            CALL XABORT(HSMG)
         ELSE IF(NG.NE.ISTATE(2)) THEN
            WRITE(HSMG,'(42HCOMCAL: ELEMENTARY CALCULATION WITH AN INV,
     1      20HALIB NB. OF GROUPS =,I7,3H NE,I7,1H.)') NG,ISTATE(2)
            CALL XABORT(HSMG)
         ENDIF
      ENDIF
      ISTATE(3)=ISTATE(3)+1
      IF(ISTATE(3).GT.ISTATE(4)) THEN
         ISTATE(4)=ISTATE(4)+10
         JPCPO=LCMLID(IPCPO,'MIXTURES',NMIL)
         DO 10 IMIL=1,NMIL
         KPCPO=LCMDIL(JPCPO,IMIL)
         LPCPO=LCMLID(KPCPO,'CALCULATIONS',ISTATE(4))
   10    CONTINUE
      ENDIF
      ICAL=ISTATE(3)
      MAXCAL=ISTATE(4)
      NISOP=ISTATE(13)
      NGFF=ISTATE(14)
      NALBP=ISTATE(15)
*----
*  RECOVER THE USER-REQUESTED PARTICULARIZED ISOTOPES
*----
      IF(NISOP.GT.MAXISO) CALL XABORT('COMCAL: MAXISO OVERFLOW.')
      IF(NISOP.GT.0) CALL LCMGTC(IPCPO,'NOMISP',8,NISOP,NOMISP)
*----
*  RECOVER THE MACRO-GEOMETRY
*----
      CALL LCMLEN(IPEDIT,'MACRO-GEOM',ILONG,ITYLCM)
      IF(ILONG.NE.0) THEN
         JPCPO=LCMLID(IPCPO,'GEOMETRIES',MAXCAL)
         KPCPO=LCMDIL(JPCPO,ICAL)
         CALL LCMSIX(IPEDIT,'MACRO-GEOM',1)
         CALL LCMEQU(IPEDIT,KPCPO)
         CALL LCMSIX(IPEDIT,' ',2)
         ISTATE(11)=1
      ENDIF
*----
*  RECOVER THE FLUX NORMALIZATION FACTOR
*----
      IF(C_ASSOCIATED(IPDEPL)) THEN
         CALL LCMGET(IPDEPL,'BURNUP-IRRAD',BIRRAD)
         BURN=BIRRAD(1)
         CALL LCMLEN(IPDEPL,'FLUX-NORM',ILONG,ITYLCM)
         IF(ILONG.EQ.0) THEN
            WRITE(HSMG,'(40HCOMCAL: THE ''FLUX-NORM'' RECORD IS NOT SE,
     1      20HT FOR BURNUP STEP AT,E12.5,14H MW-DAY/TONNE.)') BURN
            CALL XABORT(HSMG)
         ENDIF
         CALL LCMGET(IPDEPL,'FLUX-NORM',FNORM)
         IF(IMPX.GT.0) WRITE(6,100) FNORM,BURN
      ELSE
         FNORM=1.0
         IF(IMPX.GT.0) WRITE(6,110)
      ENDIF
*----
*  RECOVER THE CROSS SECTIONS AND NORMALIZE THE FLUX
*----
      CALL LCMSIX(IPEDIT,CDIRO,1)
      CALL COMMIC(IMPX,IPCPO,IPEDIT,IPEDI2,LMACRO,ICAL,MAXCAL,NMIL,
     1 NISOTS,NG,NED,NW,FNORM,LISO,NISOP,NOMISP,NGFF,NALBP,IDF,ITRES)
      ISTATE(14)=NGFF
      ISTATE(15)=NALBP
      ISTATE(16)=IDF
      CALL LCMSIX(IPEDIT,' ',2)
*----
*  UPDATE THE STATE-VECTOR
*----
      CALL LCMPUT(IPCPO,'STATE-VECTOR',NSTATE,1,ISTATE)
      RETURN
*
  100 FORMAT(45H COMCAL: NORMALIZE THE FLUX WITH THE FACTOR =,1P,E12.5,
     1 26H TAKEN FROM BURNUP STEP AT,E12.5,14H MW-DAY/TONNE.)
  110 FORMAT(36H COMCAL: THE FLUX IS NOT NORMALIZED.)
      END
