*DECK OUT
      SUBROUTINE OUT(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Simple edition module for TRIVAC-3.
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
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1): create or modification type(L_MACROLIB);
*         HENTRY(2): read-only type(L_FLUX);
*         HENTRY(3): read-only type(L_TRACK);
*         HENTRY(4): read-only type(L_MACROLIB);
*         HENTRY(5): read-only type(L_GEOM).
* IENTRY  type of each LCM object or file:
*         =1 LCM memory object; =2 XSM file; =3 sequential binary file;
*         =4 sequential ascii file.
* JENTRY  access of each LCM object or file:
*         =0 the LCM object or file is created;
*         =1 the LCM object or file is open for modifications;
*         =2 the LCM object or file is open in read-only mode.
* KENTRY  LCM object address or file unit number.
*
*Comments:
* The OUT: calling specifications are:
* MACRO2 := OUT: FLUX TRACK MACRO GEOM :: (out\_data) ;
* where
*   MACRO2 : name of the \emph{lcm} object (type L\_MACROLIB) containing the 
*     extended \emph{macrolib}.
*   FLUX   : name of the \emph{lcm} object (type L\_FLUX) containing a solution
*   TRACK  : name of the \emph{lcm} object (type L\_TRACK) containing a 
*     \emph{tracking}.
*   MACRO  : name of the \emph{lcm} object (type L\_MACROLIB) containing the 
*     reference \emph{macrolib}.
*   GEOM   : name of the \emph{lcm} object (type L\_GEOM) containing the 
*     reference \emph{geometry}.
*   out\_data}] : structure containing the data to module OUT:
* 
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER      NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      TYPE(C_PTR)  KENTRY(NENTRY)
      CHARACTER    HENTRY(NENTRY)*12
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      CHARACTER TEXT12*12,TITLE*72,HTRACK*12,HSIGN*12
      INTEGER IGP(NSTATE)
      TYPE(C_PTR) IPMAC1,IPMAC2,IPFLUX,IPTRK,IPGEOM
      INTEGER, DIMENSION(:),ALLOCATABLE :: MAT,IDL
      REAL, DIMENSION(:),ALLOCATABLE :: VOL
*----
*  PARAMETER VALIDATION.
*----
      IF(NENTRY.LE.1) CALL XABORT('OUT: TWO PARAMETERS EXPECTED.')
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2)) CALL XABORT('OUT: LCM '
     1 //'OBJECT EXPECTED AT LHS.')
      IF(JENTRY(1).NE.0) CALL XABORT('OUT: ENTRY IN CREATE MODE EXPECT'
     1 //'ED.')
      IF((JENTRY(2).NE.2).OR.((IENTRY(2).NE.1).AND.(IENTRY(2).NE.2)))
     1 CALL XABORT('OUT: LCM OBJECT IN READ-ONLY MODE EXPECTED AT RHS.')
      IPMAC2=KENTRY(1)
      IPFLUX=KENTRY(2)
      CALL LCMGTC(IPFLUX,'SIGNATURE',12,1,HSIGN)
      TEXT12=HENTRY(2)
      IF(HSIGN.NE.'L_FLUX') THEN
         CALL XABORT('OUT: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1   '. L_FLUX EXPECTED.')
      ENDIF
      HSIGN='L_MACROLIB'
      CALL LCMPTC(IPMAC2,'SIGNATURE',12,1,HSIGN)
      CALL LCMPTC(IPMAC2,'LINK.FLUX',12,1,TEXT12)
*----
*  RECOVER IPGEOM, IPMAC1 AND IPTRK POINTERS.
*----
      CALL LCMGTC(IPFLUX,'LINK.TRACK',12,1,TEXT12)
      DO 10 I=1,NENTRY
      IF(HENTRY(I).EQ.TEXT12) THEN
         IPTRK=KENTRY(I)
         CALL LCMGTC(IPTRK,'SIGNATURE',12,1,HSIGN)
         IF(HSIGN.NE.'L_TRACK') THEN
            TEXT12=HENTRY(I)
            CALL XABORT('OUT: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1      '. L_TRACK EXPECTED.')
         ENDIF
         GO TO 20
      ENDIF
   10 CONTINUE
      CALL XABORT('OUT: UNABLE TO FIND A POINTER TO L_TRACK.')
   20 CALL LCMGTC(IPFLUX,'LINK.MACRO',12,1,TEXT12)
      DO 50 I=1,NENTRY
      IF(HENTRY(I).EQ.TEXT12) THEN
         IPMAC1=KENTRY(I)
         CALL LCMGTC(IPMAC1,'SIGNATURE',12,1,HSIGN)
         IF(HSIGN.NE.'L_MACROLIB') THEN
            TEXT12=HENTRY(I)
            CALL XABORT('OUT: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1      '. L_MACROLIB EXPECTED.')
         ENDIF
         GO TO 60
      ENDIF
   50 CONTINUE
      CALL XABORT('OUT: UNABLE TO FIND A POINTER TO L_MACROLIB.')
   60 DO 70 I=1,NENTRY
      CALL LCMLEN(KENTRY(I),'SIGNATURE',ILONG,ITYLCM)
      IF(ILONG.NE.0) THEN
         CALL LCMGTC(KENTRY(I),'SIGNATURE',12,1,HSIGN)
         IF(HSIGN.EQ.'L_GEOM') THEN
            IPGEOM=KENTRY(I)
            GO TO 80
         ENDIF
      ENDIF
   70 CONTINUE
      CALL XABORT('OUT: UNABLE TO FIND A POINTER TO L_GEOM.')
   80 CALL LCMGET(IPMAC1,'STATE-VECTOR',IGP)
      NGRP=IGP(1)
      NBMIX=IGP(2)
      NL=IGP(3)
      NBFIS=IGP(4)
      NALBP=IGP(8)
*----
*  RECOVER GENERAL TRACKING INFORMATION.
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',IGP)
      NEL=IGP(1)
      NUN=IGP(2)
      IELEM=IGP(8)
      ICOL=IGP(9)
      MAXNEL=NEL
      CALL LCMLEN(IPTRK,'KEYFLX',LKFL,ITYLCM)
      ALLOCATE(MAT(NEL),VOL(NEL),IDL(LKFL))
      CALL LCMGET(IPTRK,'MATCOD',MAT)
      CALL LCMGET(IPTRK,'VOLUME',VOL)
      CALL LCMGET(IPTRK,'KEYFLX',IDL)
      CALL LCMLEN(IPTRK,'TITLE',LENGT,ITYLCM)
      IF(LENGT.GT.0) THEN
         CALL LCMGTC(IPTRK,'TITLE',72,1,TITLE)
         CALL LCMPTC(IPMAC2,'TITLE',72,1,TITLE)
      ELSE
         TITLE='*** NO TITLE PROVIDED ***'
      ENDIF
*----
*  FIND TYPE OF TRACKING.
*----
      CALL LCMGTC(IPTRK,'TRACK-TYPE',12,1,HTRACK)
*----
*  EDITION.
*----
      CALL OUTDRV(IPGEOM,IPMAC1,IPFLUX,IPMAC2,MAXNEL,NBMIX,NL,
     1 NBFIS,NGRP,NEL,NUN,NALBP,HTRACK,IELEM,ICOL,MAT,VOL,IDL,
     2 TITLE)
*----
*  RELEASE GENERAL TRACKING INFORMATION.
*----
      DEALLOCATE(IDL,VOL,MAT)
      RETURN
      END
