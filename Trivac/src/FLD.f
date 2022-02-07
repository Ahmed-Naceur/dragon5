*DECK FLD
      SUBROUTINE FLD(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Multigroup flux solution operator for BIVAC and TRIVAC.
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
*         HENTRY(1): create or modification type(L_FLUX);
*         HENTRY(2): read-only type(L_SYSTEM);
*         HENTRY(3): read-only type(L_TRACK);
*         HENTRY(4): optional read-only type(L_MACROLIB).
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
* The FLUD:  calling specifications are:
* FLUX := FLUD: [ FLUX ] SYST TRACK [ MACRO ] :: (flud\_data) ;
* where
*   FLUX  : name of the \emph{lcm} object (type L\_FLUX) containing the 
*     solution. If FLUX appears on the RHS, the solution previously stored in 
*     FLUX is used to initialize the new iterative process; otherwise, a uniform
*     unknown vector is used.
*   SYST  : name of the \emph{lcm} object (type L\_SYSTEM) containing the 
*     system matrices.
*   TRACK : name of the \emph{lcm} object (type L\_TRACK) containing the 
*     \emph{tracking}.
*   MACRO : name of the optional \emph{lcm} object (type L\_MACROLIB) 
*     containing the cross sections. This object is only used to set a link to 
*     the \emph{macrolib} name inside the \emph{flux} object. By default, the 
*     name of the \emph{macrolib} is recovered   from the link in the 
*     \emph{system} object.
*   flud\_data}] : structure containing the data to module FLUD:}
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
      CHARACTER TEXT12*12,TITLE*72,CMODUL*12,HSIGN*12
      LOGICAL REC,LREL
      INTEGER IGP(NSTATE),ITR(NSTATE)
      TYPE(C_PTR) IPTRK,IPSYS,IPFLUX
      INTEGER, DIMENSION(:), ALLOCATABLE :: MAT,IDL
      REAL, DIMENSION(:), ALLOCATABLE :: VOL
*----
*  PARAMETER VALIDATION
*----
      LREL=(JENTRY(1).EQ.1)
      IF(NENTRY.LE.1) CALL XABORT('FLD: TWO PARAMETERS EXPECTED.')
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2)) CALL XABORT('FLD: LCM '
     1 //'OBJECT EXPECTED AT LHS.')
      IF((JENTRY(1).NE.0).AND.(JENTRY(1).NE.1)) CALL XABORT('FLD: ENTR'
     1 //'Y IN CREATE OR MODIFICATION MODE EXPECTED.')
      IF((JENTRY(2).NE.2).OR.((IENTRY(2).NE.1).AND.(IENTRY(2).NE.2)))
     1 CALL XABORT('FLD: LCM OBJECT IN READ-ONLY MODE EXPECTED AT RHS.')
      IF(JENTRY(1).EQ.1) THEN
         CALL LCMGTC(KENTRY(1),'SIGNATURE',12,1,HSIGN)
         IF(HSIGN.NE.'L_FLUX') THEN
            TEXT12=HENTRY(1)
            CALL XABORT('FLD: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1      '. L_FLUX EXPECTED.')
         ENDIF
      ENDIF
      CALL LCMGTC(KENTRY(2),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_SYSTEM') THEN
         TEXT12=HENTRY(2)
         CALL XABORT('FLD: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1   '. L_SYSTEM EXPECTED.')
      ENDIF
      HSIGN='L_FLUX'
      CALL LCMPTC(KENTRY(1),'SIGNATURE',12,1,HSIGN)
      TEXT12=HENTRY(2)
      CALL LCMPTC(KENTRY(1),'LINK.SYSTEM',12,1,TEXT12)
      IPFLUX=KENTRY(1)
      IPSYS=KENTRY(2)
      REC=(JENTRY(1).EQ.1)
*----
*  RECOVER IPTRK POINTER AND VALIDATE IT
*----
      IF(NENTRY.EQ.4) THEN
         CALL LCMGTC(KENTRY(4),'SIGNATURE',12,1,HSIGN)
         IF(HSIGN.NE.'L_MACROLIB') THEN
            TEXT12=HENTRY(4)
            CALL XABORT('FLD: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1      '. L_MACROLIB EXPECTED.')
         ENDIF
         TEXT12=HENTRY(4)
      ELSE
         CALL LCMGTC(IPSYS,'LINK.MACRO',12,1,TEXT12)
      ENDIF
      CALL LCMPTC(KENTRY(1),'LINK.MACRO',12,1,TEXT12)
      CALL LCMGTC(IPSYS,'LINK.TRACK',12,1,TEXT12)
      CALL LCMPTC(KENTRY(1),'LINK.TRACK',12,1,TEXT12)
      DO 10 I=1,NENTRY
      IF(HENTRY(I).EQ.TEXT12) THEN
         IPTRK=KENTRY(I)
         CALL LCMGTC(IPTRK,'SIGNATURE',12,1,HSIGN)
         IF(HSIGN.NE.'L_TRACK') THEN
            TEXT12=HENTRY(I)
            CALL XABORT('FLD: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1      '. L_TRACK EXPECTED.')
         ENDIF
         CALL LCMGTC(IPTRK,'TRACK-TYPE',12,1,CMODUL)
         IF((IENTRY(I).NE.1).AND.(IENTRY(I).NE.2)) CALL XABORT('FLD: L'
     1   //'CM OBJECT EXPECTED TO CONTAIN THE TRACKING.')
         GO TO 20
      ENDIF
   10 CONTINUE
      CALL XABORT('FLD: UNABLE TO FIND A POINTER TO TRACKING.')
*----
*  RECOVER GENERAL TRACKING INFORMATION
*----
   20 CALL LCMGET(IPTRK,'STATE-VECTOR',IGP)
      NEL=IGP(1)
      NUN=IGP(2)
      NLF=0
      IF(CMODUL.EQ.'BIVAC') THEN
         NLF=IGP(14)
      ELSE IF(CMODUL.EQ.'TRIVAC') THEN
         NLF=IGP(30)
      ENDIF
      ALLOCATE(MAT(NEL),VOL(NEL),IDL(NEL))
      CALL LCMGET(IPTRK,'MATCOD',MAT)
      CALL LCMGET(IPTRK,'VOLUME',VOL)
      CALL LCMGET(IPTRK,'KEYFLX',IDL)
      CALL LCMLEN(IPTRK,'TITLE',LENGT,ITYLCM)
      IF(LENGT.GT.0) THEN
         CALL LCMGTC(IPTRK,'TITLE',72,1,TITLE)
      ELSE
         TITLE='*** NO TITLE PROVIDED ***'
      ENDIF
*----
*  RECOVER GENERAL L_SYSTEM INFORMATION
*----
      CALL LCMGET(IPSYS,'STATE-VECTOR',ITR)
      NGRP=ITR(1)
      LL4=ITR(2)
      ITY=ITR(4)
      NBMIX=ITR(7)
      IF((ITY.EQ.11).OR.(ITY.EQ.13)) LL4=LL4*NLF/2
*----
*  COMPUTE THE FLUX
*----
      CALL FLDDRV(CMODUL,IPTRK,IPSYS,REC,NEL,LL4,ITY,NUN,NBMIX,MAT,VOL,
     1 IDL,NGRP,TITLE,LREL,IPFLUX)
*----
*  RELEASE GENERAL TRACKING INFORMATION
*----
      DEALLOCATE(IDL,VOL,MAT)
      RETURN
      END
