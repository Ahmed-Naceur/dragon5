*DECK TRIVAA
      SUBROUTINE TRIVAA(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* TRIVAC type (3-D and ADI) system matrix assembly operator.
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
*         HENTRY(2): read-only type(L_MACROLIB) (unperturbed);
*         HENTRY(3): read-only type(L_TRACK);
*         HENTRY(4): optional read-only type(L_MACROLIB) (perturbed).
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
* The TRIVAA: calling specifications are:
* SYST := TRIVAA: [ SYST ] MACRO  TRACK [ DMACRO ] :: (trivaa\_data) ;
* where
*   SYST   : name of the \emph{lcm} object (type L\_SYSTEM) containing the 
*     system matrices. If SYST appears on the RHS, the system matrices 
*     previously stored in SYST are kept.
*   MACRO  : name of the \emph{lcm} object (type L\_MACROLIB) containing the 
*     macroscopic cross sections and diffusion coefficients.
*   TRACK  : name of the \emph{lcm} object (type L\_TRIVAC) containing the 
*     TRIVAC \emph{tracking}.
*   DMACRO : name of the \emph{lcm} object  (type L\_MACROLIB) containing 
*     derivatives or perturbations of the macroscopic cross sections and 
*     diffusion coefficients. If DMACRO is given, only the derivatives or 
*     perturbations of the system matrices are computed.
*   trivaa\_data : structure containing the data to module TRIVAA:
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
      CHARACTER TEXT4*4,TEXT11*12,TEXT12*12,HSMG*131,TITLE*72,CNAM*12
      DOUBLE PRECISION DFLOTT
      INTEGER IGP(NSTATE),IPAR(NSTATE),ITR(NSTATE)
      LOGICAL LDIFF
      TYPE(C_PTR) IPSYS,JPSYS,KPSYS,IPMACR,JPMACR,KPMACR,IPTRK,IPMACP
      INTEGER, DIMENSION(:), ALLOCATABLE :: MAT
      REAL, DIMENSION(:), ALLOCATABLE :: VOL,UN,VII
*----
*  PARAMETER VALIDATION.
*----
      IF(NENTRY.LE.2) CALL XABORT('TRIVAA: THREE PARAMETERS EXPECTED.')
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2)) CALL XABORT('TRIVAA: L'
     1 //'CM OBJECT EXPECTED AT LHS.')
      IF((JENTRY(1).NE.0).AND.(JENTRY(1).NE.1)) CALL XABORT('TRIVAA: E'
     1 //'NTRY IN CREATE OR MODIFICATION MODE EXPECTED.')
      IF((JENTRY(3).NE.2).OR.((IENTRY(3).NE.1).AND.(IENTRY(3).NE.2)))
     1 CALL XABORT('TRIVAA: LCM OBJECT IN READ-ONLY MODE EXPECTED AT F'
     2 //'IRST RHS.')
      IF((JENTRY(2).NE.2).OR.((IENTRY(2).NE.1).AND.(IENTRY(2).NE.2)))
     1 CALL XABORT('TRIVAA: LCM OBJECT IN READ-ONLY MODE EXPECTED AT S'
     2 //'ECOND RHS.')
      CALL LCMGTC(KENTRY(3),'SIGNATURE',12,1,TEXT11)
      IF(TEXT11.NE.'L_TRACK') THEN
         TEXT12=HENTRY(3)
         CALL XABORT('TRIVAA: SIGNATURE OF '//TEXT12//' IS '//TEXT11//
     1   '. L_TRACK EXPECTED.')
      ENDIF
      CALL LCMGTC(KENTRY(3),'TRACK-TYPE',12,1,TEXT11)
      IF(TEXT11.NE.'TRIVAC') THEN
         TEXT12=HENTRY(3)
         CALL XABORT('TRIVAA: TRACK-TYPE OF '//TEXT12//' IS '//TEXT11
     1   //'. TRIVAC EXPECTED.')
      ENDIF
      TEXT11='L_SYSTEM'
      IPSYS=KENTRY(1)
      CALL LCMPTC(IPSYS,'SIGNATURE',12,1,TEXT11)
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
      NLF=IGP(30)
      ISCAT=IGP(32)
      LDIFF=(ISCAT.LT.0)
      ISCAT=ABS(ISCAT)
      IF((NLF.NE.0).AND.(IGP(31).NE.1)) CALL XABORT('TRIVAA: ONLY SPN '
     1 //'DISCRETIZATIONS ARE ALLOWED.')
      ITY=2
      IF(IGP(12).EQ.2) ITY=3
      ALLOCATE(MAT(NEL),VOL(NEL))
      CALL LCMGET(IPTRK,'MATCOD',MAT)
      CALL LCMGET(IPTRK,'VOLUME',VOL)
      CALL LCMLEN(IPTRK,'TITLE',LENGT,ITYLCM)
      IF(LENGT.GT.0) THEN
         CALL LCMGTC(IPTRK,'TITLE',72,1,TITLE)
      ELSE
         TITLE='*** NO TITLE PROVIDED ***'
      ENDIF
*----
*  RECOVER MACROLIB PARAMETERS.
*----
      CALL LCMGTC(IPMACR,'SIGNATURE',12,1,TEXT11)
      IF(TEXT11.NE.'L_MACROLIB') THEN
         TEXT12=HENTRY(2)
         CALL XABORT('TRIVAA: SIGNATURE OF '//TEXT12//' IS '//TEXT11//
     1   '. L_MACROLIB EXPECTED.')
      ENDIF
      CALL LCMGET(IPMACR,'STATE-VECTOR',IPAR)
      NGRP=IPAR(1)
      NBMIX=IPAR(2)
      NANI=IPAR(3)
      NBFIS=IPAR(4)
      NALBP=IPAR(8)
      IF(IGP(4).GT.NBMIX) THEN
         WRITE(HSMG,'(46HTRIVAA: THE NUMBER OF MIXTURES IN THE TRACKING,
     1   2H (,I5,51H) IS GREATER THAN THE NUMBER OF MIXTURES IN THE MAC,
     2   7HROLIB (,I5,2H).)') IGP(4),NBMIX
         CALL XABORT(HSMG)
      ENDIF
*
      IMPX=1
      IASM=0
      IPR=0
      IUNIT=0
      IOVEL=0
      NSTEP=0
   10 CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
      IF(INDIC.EQ.10) GO TO 30
      IF(INDIC.NE.3) CALL XABORT('TRIVAA: CHARACTER DATA EXPECTED.')
      IF(TEXT4.EQ.'EDIT') THEN
         CALL REDGET(INDIC,IMPX,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('TRIVAA: INTEGER DATA EXPECTED(1).')
      ELSE IF(TEXT4.EQ.'SKIP') THEN
*        OPTION TO SKIP THE SYSTEM MATRIX ASSEMBLY (DO NOT SKIP THE
*        LDLT FACTORIZATION).
         IASM=1
      ELSE IF(TEXT4.EQ.'DERI') THEN
         IPR=1
         WRITE(6,'(/43H TRIVAA: USE DERIVATIVE OF SYSTEM MATRICES.)')
      ELSE IF(TEXT4.EQ.'PERT') THEN
         IPR=2
         WRITE(6,'(/41H TRIVAA: PERTURBATION OF SYSTEM MATRICES.)')
      ELSE IF(TEXT4.EQ.'UNIT') THEN
*       COMPUTE THE UNITARY WEIGHTING MATRIX.
        IUNIT=1
        ALLOCATE(UN(NBMIX))
        CALL XDRSET(UN,NBMIX,1.0)
        CALL TRIDIG('RM',IPTRK,IPSYS,IMPX,NBMIX,NEL,0,MAT,VOL,UN)
        DEALLOCATE(UN)
      ELSE IF(TEXT4.EQ.'OVEL') THEN
*        COMPUTE THE RECIPROCAL NEUTRON VELOCITIES MATRIX.
         IOVEL=1
         JPMACR=LCMGID(IPMACR,'GROUP')
         ALLOCATE(VII(NBMIX))
         DO 25 IGR=1,NGRP
         KPMACR=LCMGIL(JPMACR,IGR)
         CALL LCMLEN(KPMACR,'OVERV',LENGT,ITYLCM)
         IF(LENGT.EQ.0) THEN
            CALL XABORT('TRIVAA: NO ''VELOCITY'' INFORMATION.')
         ELSE IF(LENGT.GT.NBMIX) THEN
            CALL XABORT('TRIVAA: INVALID LENGTH FOR ''VELOCITY'' IN'
     1      //'FORMATION.')
         ENDIF
         CALL LCMGET(KPMACR,'OVERV',VII)
         WRITE (CNAM,'(1HV,2I3.3)') IGR,IGR
         CALL TRIDIG(CNAM,IPTRK,IPSYS,IMPX,NBMIX,NEL,0,MAT,VOL,VII)
   25    CONTINUE
         DEALLOCATE(VII)
      ELSE IF(TEXT4.EQ.';') THEN
         GO TO 30
      ELSE
         CALL XABORT('TRIVAA: '//TEXT4//' IS AN INVALID KEY WORD.')
      ENDIF
      GO TO 10
*----
*  2-MACROLIBS PERTURBATION CALCULATION.
*----
   30 IF(IPR.GT.0) THEN
         IF(NENTRY.LE.3) CALL XABORT('TRIVAA: 4 PARAMETERS EXPECTED WIT'
     1   //'H DERI OR PERT OPTIONS.')
         IF((JENTRY(4).NE.2).OR.((IENTRY(4).NE.1).AND.(IENTRY(4).NE.2)))
     1   CALL XABORT('TRIVAA: LINKED LIST OR XSM FILE IN READ-ONLY MODE'
     2   //' EXPECTED AT THIRD RHS.')
         IPMACP=KENTRY(4)
         CALL LCMGTC(IPMACP,'SIGNATURE',12,1,TEXT11)
         IF(TEXT11.NE.'L_MACROLIB') THEN
            TEXT12=HENTRY(4)
            CALL XABORT('TRIVAA: SIGNATURE OF '//TEXT12//' IS '
     1      //TEXT11//'. L_MACROLIB EXPECTED.')
         ENDIF
         CALL LCMGET(IPMACP,'STATE-VECTOR',IPAR)
         NSTEP=IPAR(11)
         IF((IPAR(1).NE.NGRP).OR.(IPAR(2).GT.NBMIX)) THEN
            WRITE(HSMG,'(43HTRIVAA: INCONSISTENT PERTURBATION MACROLIB ,
     1      1H'',A12,8H''. NGRP=,2I5,7H NBMIX=,2I9)') HENTRY(4),IPAR(1),
     2      NGRP,IPAR(2),NBMIX
            CALL XABORT(HSMG)
         ENDIF
      ENDIF
*----
*  SET THE STATE VECTOR FOR THE L_SYSTEM OBJECT
*----
      CALL XDISET(ITR,NSTATE,0)
      ITR(1)=NGRP
      ITR(2)=IGP(11)
      ITR(3)=0
      ITR(4)=ITY
      IF((NLF.GT.0).AND.(ITY.GE.3)) ITR(4)=10+ITR(4)
      IF(IUNIT.EQ.1) ITR(5)=1
      ITR(6)=NSTEP
      ITR(7)=NBMIX
      NAN=MIN(ISCAT,NANI)
      ITR(8)=NLF
      ITR(9)=IPR
      CALL LCMPUT(IPSYS,'STATE-VECTOR',NSTATE,1,ITR)
*----
*  SYSTEM MATRIX ASSEMBLY.
*----
      IF((IASM.EQ.0).AND.(IPR.EQ.0)) THEN
         IF(NLF.EQ.0) THEN
*           DIFFUSION THEORY.
            CALL TRISYS(IPTRK,IPMACR,IPMACR,IPSYS,IMPX,NGRP,NEL,NBFIS,
     1      NALBP,IPR,MAT,VOL,NBMIX)
         ELSE
*           SIMPLIFIED PN THEORY.
            CALL TRISPS(IPTRK,IPMACR,IPMACR,IPSYS,IMPX,NGRP,NEL,NLF,
     1      NAN,NBFIS,NALBP,LDIFF,IPR,MAT,VOL,NBMIX)
         ENDIF
      ELSE IF((IASM.EQ.1).AND.(IPR.EQ.0)) THEN
*        PERFORM FACTORIZATION WITHOUT ASSEMBLY.
         DO 40 I=1,NGRP
         WRITE(TEXT11,'(1HA,2I3.3)') I,I
         CALL MTLDLF(TEXT11,IPTRK,IPSYS,ITY,IMPX)
   40    CONTINUE
      ELSE IF((IPR.GT.0).AND.(NSTEP.EQ.0)) THEN
*        ASSEMBLY OF PERTURBED SYSTEM MATRICES (NO STEP DIRECTORIES).
         IF(NLF.EQ.0) THEN
*           DIFFUSION THEORY.
            CALL TRISYS(IPTRK,IPMACR,IPMACP,IPSYS,IMPX,NGRP,NEL,NBFIS,
     1      NALBP,IPR,MAT,VOL,NBMIX)
         ELSE
*           SIMPLIFIED PN THEORY.
            CALL TRISPS(IPTRK,IPMACR,IPMACP,IPSYS,IMPX,NGRP,NEL,NLF,
     1      NAN,NBFIS,NALBP,LDIFF,IPR,MAT,VOL,NBMIX)
         ENDIF
      ELSE IF(NSTEP.GT.0) THEN
*        ASSEMBLY OF PERTURBED SYSTEM MATRICES (WITH STEP DIRECTORIES).
         JPMACR=LCMGID(IPMACP,'STEP')
         JPSYS=LCMLID(IPSYS,'STEP',NSTEP)
         DO 50 ISTEP=1,NSTEP
         KPMACR=LCMGIL(JPMACR,ISTEP)
         KPSYS=LCMDIL(JPSYS,ISTEP)
         CALL LCMPUT(KPSYS,'STATE-VECTOR',NSTATE,1,ITR)
         IF(NLF.EQ.0) THEN
*           DIFFUSION THEORY.
            CALL TRISYS(IPTRK,IPMACR,KPMACR,KPSYS,IMPX,NGRP,NEL,NBFIS,
     1      NALBP,IPR,MAT,VOL,NBMIX)
         ELSE
*           SIMPLIFIED PN THEORY.
            CALL TRISPS(IPTRK,IPMACR,KPMACR,KPSYS,IMPX,NGRP,NEL,NLF,
     1      NAN,NBFIS,NALBP,LDIFF,IPR,MAT,VOL,NBMIX)
         ENDIF
   50    CONTINUE
      ELSE
         CALL XABORT('TRIVAA: INVALID REQUEST.')
      ENDIF
*
      IF(IMPX.GE.3) CALL LCMLIB(IPSYS)
*----
*  RELEASE GENERAL TRACKING INFORMATION.
*----
      DEALLOCATE(VOL,MAT)
      RETURN
      END
