*DECK BIVACA
      SUBROUTINE BIVACA(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* BIVAC type (2-D) system matrix assembly operator.
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
*         HENTRY(1): create or modification type(L_SYSTEM);
*         HENTRY(2): read-only type(L_MACROLIB);
*         HENTRY(3): read-only type(L_TRACK).
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
* The BIVACA: calling specifications are:
* SYST := BIVACA: [ SYST} ] MACRO  TRACK :: (bivaca\_data) ;
* where
*   SYST  : name of the \emph{lcm} object (type L\_SYSTEM) containing the 
*     system matrices. If SYST appears on the RHS, the system matrices 
*     previously stored in SYST} are kept.
*   MACRO : name of the \emph{lcm} object (type L\_MACROLIB) containing the 
*     macroscopic cross sections and diffusion coefficients.
*   TRACK : name of the \emph{lcm} object (type L\_BIVAC) containing the BIVAC
*     \emph{tracking}.
*   bivaca\_data : structure containing the data to module BIVACA:.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      CHARACTER HENTRY(NENTRY)*12
      TYPE(C_PTR) KENTRY(NENTRY)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40)
      CHARACTER TEXT4*4,HSIGN*12,TEXT12*12,HSMG*131,CNAM*12
      DOUBLE PRECISION DFLOTT
      INTEGER IGP(NSTATE),IPAR(NSTATE),IBR(NSTATE)
      LOGICAL LDIFF
      TYPE(C_PTR) IPSYS,IPMACR,JPMACR,KPMACR,IPTRK
      INTEGER, DIMENSION(:), ALLOCATABLE :: MAT
      REAL, DIMENSION(:), ALLOCATABLE :: VOL,UN,VII,GAMMA
*----
*  PARAMETER VALIDATION.
*----
      IF(NENTRY.LE.2) CALL XABORT('BIVACA: THREE PARAMETERS EXPECTED.')
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2)) CALL XABORT('BIVACA: L'
     1 //'CM OBJECT EXPECTED AT LHS.')
      IF((JENTRY(1).NE.0).AND.(JENTRY(1).NE.1)) CALL XABORT('BIVACA: E'
     1 //'NTRY IN CREATE OR MODIFICATION MODE EXPECTED.')
      IF((JENTRY(2).NE.2).OR.((IENTRY(2).NE.1).AND.(IENTRY(2).NE.2)))
     1 CALL XABORT('BIVACA: LCM OBJECT IN READ-ONLY MODE EXPECTED AT S'
     2 //'ECOND RHS.')
      IF((JENTRY(3).NE.2).OR.((IENTRY(3).NE.1).AND.(IENTRY(3).NE.2)))
     1 CALL XABORT('BIVACA: LCM OBJECT IN READ-ONLY MODE EXPECTED AT F'
     2 //'IRST RHS.')
      CALL LCMGTC(KENTRY(3),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_TRACK') THEN
         TEXT12=HENTRY(3)
         CALL XABORT('BIVACA: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1   '. L_TRACK EXPECTED.')
      ENDIF
      CALL LCMGTC(KENTRY(3),'TRACK-TYPE',12,1,HSIGN)
      IF(HSIGN.NE.'BIVAC') THEN
         TEXT12=HENTRY(3)
         CALL XABORT('BIVACA: TRACK-TYPE OF '//TEXT12//' IS '//HSIGN//
     1   '. BIVAC EXPECTED.')
      ENDIF
      HSIGN='L_SYSTEM'
      IPSYS=KENTRY(1)
      CALL LCMPTC(IPSYS,'SIGNATURE',12,1,HSIGN)
      IPMACR=KENTRY(2)
      IPTRK=KENTRY(3)
      TEXT12=HENTRY(2)
      CALL LCMPTC(IPSYS,'LINK.MACRO',12,1,TEXT12)
      TEXT12=HENTRY(3)
      CALL LCMPTC(IPSYS,'LINK.TRACK',12,1,TEXT12)
*----
*  RECOVER GENERAL TRACKING INFORMATION.
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',IGP)
      NEL=IGP(1)
      NLF=IGP(14)
      ISCAT=IGP(16)
      LDIFF=(ISCAT.LT.0)
      ISCAT=ABS(ISCAT)
      IF((NLF.NE.0).AND.(IGP(15).NE.1)) CALL XABORT('BIVACA: ONLY SPN '
     1 //'DISCRETIZATIONS ARE ALLOWED.')
      ALLOCATE(MAT(NEL),VOL(NEL))
      CALL LCMGET(IPTRK,'MATCOD',MAT)
      CALL LCMGET(IPTRK,'VOLUME',VOL)
*----
*  RECOVER MACROLIB PARAMETERS.
*----
      CALL LCMGTC(IPMACR,'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_MACROLIB') THEN
         TEXT12=HENTRY(2)
         CALL XABORT('BIVACA: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1   '. L_MACROLIB EXPECTED.')
      ENDIF
      CALL LCMGET(IPMACR,'STATE-VECTOR',IPAR)
      NGRP=IPAR(1)
      NBMIX=IPAR(2)
      NANI=IPAR(3)
      NBFIS=IPAR(4)
      NALBP=IPAR(8)
      IF(IGP(4).GT.NBMIX) THEN
         WRITE(HSMG,'(46HBIVACA: THE NUMBER OF MIXTURES IN THE TRACKING,
     1   2H (,I5,51H) IS GREATER THAN THE NUMBER OF MIXTURES IN THE MAC,
     2   7HROLIB (,I5,2H).)') IGP(4),NBMIX
         CALL XABORT(HSMG)
      ENDIF
*
      IMPX=1
      IUNIT=0
      IOVEL=0
   10 CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
      IF(INDIC.EQ.10) GO TO 40
      IF(INDIC.NE.3) CALL XABORT('BIVACA: CHARACTER DATA EXPECTED.')
      IF(TEXT4.EQ.'EDIT') THEN
         CALL REDGET(INDIC,IMPX,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('BIVACA: INTEGER DATA EXPECTED(1).')
      ELSE IF(TEXT4.EQ.'UNIT') THEN
*       COMPUTE THE UNITARY WEIGHTING MATRIX.
        IUNIT=1
        ALLOCATE(UN(NBMIX),GAMMA(NALBP))
        CALL XDRSET(UN,NBMIX,1.0)
        CALL BIVASM('RM',1,IPTRK,IPSYS,IMPX,NBMIX,NEL,0,1,0,MAT,VOL,
     1  GAMMA,UN)
        DEALLOCATE(GAMMA,UN)
      ELSE IF(TEXT4.EQ.'OVEL') THEN
*        COMPUTE THE RECIPROCAL NEUTRON VELOCITIES MATRIX.
         IOVEL=1
         JPMACR=LCMGID(IPMACR,'GROUP')
         ALLOCATE(VII(NBMIX),GAMMA(NALBP))
         DO 30 IGR=1,NGRP
         KPMACR=LCMGIL(JPMACR,IGR)
         CALL LCMLEN(KPMACR,'OVERV',LENGT,ITYLCM)
         IF(LENGT.EQ.0) THEN
            CALL XABORT('BIVACA: NO ''VELOCITY'' INFORMATION.')
         ELSE IF(LENGT.GT.NBMIX) THEN
            CALL XABORT('BIVACA: INVALID LENGTH FOR ''VELOCITY'' IN'
     1      //'FORMATION.')
         ENDIF
         CALL LCMGET(KPMACR,'OVERV',VII)
         WRITE(CNAM,'(1HV,2I3.3)') IGR,IGR
         CALL BIVASM(CNAM,1,IPTRK,IPSYS,IMPX,NBMIX,NEL,0,1,0,MAT,VOL,
     1   GAMMA,VII)
   30    CONTINUE
         DEALLOCATE(GAMMA,VII)
      ELSE IF(TEXT4.EQ.';') THEN
         GO TO 40
      ELSE
         CALL XABORT('BIVACA: '//TEXT4//' IS AN INVALID KEY WORD.')
      ENDIF
      GO TO 10
*----
*  SET THE STATE VECTOR FOR THE L_SYSTEM OBJECT
*----
   40 CALL XDISET(IBR,NSTATE,0)
      IBR(1)=NGRP
      IBR(2)=IGP(11)
      IBR(4)=1
      IF(NLF.GT.0) IBR(4)=11
      IF(IUNIT.EQ.1) IBR(5)=1
      IBR(7)=NBMIX
      NAN=MIN(ISCAT,NANI)
      IBR(8)=NLF
      CALL LCMPUT(IPSYS,'STATE-VECTOR',NSTATE,1,IBR)
*----
*  BIVAC SYSTEM MATRIX ASSEMBLY.
*----
      IF(NLF.EQ.0) THEN
*        DIFFUSION THEORY.
         CALL BIVSYS(IPTRK,IPMACR,IPSYS,IMPX,NGRP,NEL,NBFIS,NALBP,MAT,
     1   VOL,NBMIX)
      ELSE
*        SIMPLIFIED PN THEORY.
         CALL BIVSPS(IPTRK,IPMACR,IPSYS,IMPX,NGRP,NEL,NLF,NAN,NBFIS,
     1   NALBP,LDIFF,MAT,VOL,NBMIX)
      ENDIF
*
      IF(IMPX.GE.3) CALL LCMLIB(IPSYS)
*----
*  RELEASE GENERAL TRACKING INFORMATION.
*----
      DEALLOCATE(VOL,MAT)
      RETURN
      END
