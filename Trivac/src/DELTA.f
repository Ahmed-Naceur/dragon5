*DECK DELTA
      SUBROUTINE DELTA(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* calculation of direct or adjoint source components for a fixed source
* eigenvalue problem.
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
*         HENTRY(1): create or modification type(L_SOURCE) (GPT source)
*         HENTRY(2): read-only type(L_FLUX) => unperturbed solution
*         HENTRY(3): read-only type(L_SYSTEM) => unperturbed matrices
*         HENTRY(4): read-only type(L_SYSTEM) => perturbed matrices
*         HENTRY(5): read-only type(L_TRACK) => tracking.
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
* The DELTA: calling specifications are:
* GPT := DELTA: [ GPT ] FLUX0 SYST0 DSYST TRACK :: (delta\_data) ;
* where
*   GPT   : name of the \emph{lcm} object (type L\_GPT) containing the fixed 
*     source. If GPT appears on the RHS, this information is used to initialize 
*     the state vector.
*   FLUX0 : name of the \emph{lcm} object (type L\_FLUX) containing the 
*     unperturbed flux.
*   SYST0 : name of the \emph{lcm} object (type L\_SYSTEM) containing the 
*     unperturbed system matrices.
*   DSYST : name of the \emph{lcm} object (type L\_SYSTEM) containing a 
*     perturbation to the system matrices.
*   TRACK : name of the \emph{lcm} object (type L\_TRACK) containing the 
*     \emph{tracking}.
*   delta\_data}] : structure containing the data to module DELTA:}
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
      PARAMETER (NSTATE=40)
      CHARACTER TEXT12*12,HSIGN*12,CMODUL*12
      LOGICAL REC
      INTEGER ISTATE(NSTATE)
      TYPE(C_PTR) IPGPT,IPFLU0,IPSYS0,IPSYSP,IPTRK
*----
*  PARAMETER VALIDATION.
*----
      IF(NENTRY.LE.4) CALL XABORT('DELTA: FIVE PARAMETERS EXPECTED.')
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2)) CALL XABORT('DELTA: LC'
     1 //'M OBJECT EXPECTED AT LHS.')
      IF((JENTRY(1).NE.0).AND.(JENTRY(1).NE.1)) CALL XABORT('DELTA: EN'
     1 //'TRY IN CREATE OR MODIFICATION MODE EXPECTED.')
      IF((JENTRY(2).NE.2).OR.((IENTRY(2).NE.1).AND.(IENTRY(2).NE.2)))
     1 CALL XABORT('DELTA: LCM OBJECT IN READ-ONLY MODE EXPECTED AT FI'
     2 //'RST RHS.')
      IF((JENTRY(3).NE.2).OR.((IENTRY(3).NE.1).AND.(IENTRY(3).NE.2)))
     1 CALL XABORT('DELTA: LCM OBJECT IN READ-ONLY MODE EXPECTED AT SE'
     2 //'COND RHS.')
      IF((JENTRY(4).NE.2).OR.((IENTRY(4).NE.1).AND.(IENTRY(4).NE.2)))
     1 CALL XABORT('DELTA: LCM OBJECT IN READ-ONLY MODE EXPECTED AT TH'
     2 //'IRD RHS.')
      IF((JENTRY(5).NE.2).OR.((IENTRY(5).NE.1).AND.(IENTRY(5).NE.2)))
     1 CALL XABORT('DELTA: LCM OBJECT IN READ-ONLY MODE EXPECTED AT FO'
     2 //'URTH RHS.')
      REC=(JENTRY(1).EQ.1)
      IF(REC) THEN
         CALL LCMGTC(KENTRY(1),'SIGNATURE',12,1,HSIGN)
         IF(HSIGN.NE.'L_SOURCE') THEN
            TEXT12=HENTRY(1)
            CALL XABORT('DELTA: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1      '. L_SOURCE EXPECTED.')
         ENDIF
      ELSE
         HSIGN='L_SOURCE'
         CALL LCMPTC(KENTRY(1),'SIGNATURE',12,1,HSIGN)
      ENDIF
      CALL LCMGTC(KENTRY(2),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_FLUX') THEN
         TEXT12=HENTRY(2)
         CALL XABORT('DELTA: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1   '. L_FLUX EXPECTED.')
      ENDIF
      CALL LCMGTC(KENTRY(3),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_SYSTEM') THEN
         TEXT12=HENTRY(3)
         CALL XABORT('DELTA: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1   '. L_SYSTEM EXPECTED.')
      ENDIF
      CALL LCMGTC(KENTRY(4),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_SYSTEM') THEN
         TEXT12=HENTRY(4)
         CALL XABORT('DELTA: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1   '. L_SYSTEM EXPECTED.')
      ENDIF
      CALL LCMGTC(KENTRY(5),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_TRACK') THEN
         TEXT12=HENTRY(5)
         CALL XABORT('DELTA: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1   '. L_TRACK EXPECTED.')
      ENDIF
      CALL LCMGTC(KENTRY(2),'LINK.SYSTEM',12,1,TEXT12)
      IF(TEXT12.NE.HENTRY(3)) CALL XABORT('DELTA: OBJECT '//HENTRY(3)//
     1 ' IS NOT AN UNPERTURBED SYSTEM OBJECT.')
      CALL LCMGTC(KENTRY(2),'LINK.TRACK',12,1,TEXT12)
      IF(TEXT12.NE.HENTRY(5)) CALL XABORT('DELTA: OBJECT '//HENTRY(3)//
     1 ' IS NOT A TRACKING OBJECT.')
      TEXT12=HENTRY(2)
      CALL LCMPTC(KENTRY(1),'LINK.FLUX',12,1,TEXT12)
      TEXT12=HENTRY(3)
      CALL LCMPTC(KENTRY(1),'LINK.SYSTEM',12,1,TEXT12)
      TEXT12=HENTRY(4)
      CALL LCMPTC(KENTRY(1),'LINK.TRACK',12,1,TEXT12)
      IPGPT=KENTRY(1)
      IPFLU0=KENTRY(2)
      IPSYS0=KENTRY(3)
      IPSYSP=KENTRY(4)
      IPTRK=KENTRY(5)
*----
*  RECOVER GENERAL TRACKING INFORMATION.
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      NUN=ISTATE(2)
      CALL LCMGTC(IPTRK,'TRACK-TYPE',12,1,CMODUL)
      IF(CMODUL.NE.'TRIVAC') CALL XABORT('DELTA: TRIVAC TRACKING EXPEC'
     1 //'TED.')
      CALL LCMGET(IPSYS0,'STATE-VECTOR',ISTATE)
      NGRP=ISTATE(1)
      LL4=ISTATE(2)
      CALL LCMGET(IPSYSP,'STATE-VECTOR',ISTATE)
      IF(ISTATE(1).NE.NGRP) CALL XABORT('DELTA: INVALID NGRP.')
      IF(ISTATE(2).NE.LL4) CALL XABORT('DELTA: INVALID LL4.')
      NSTEP=ISTATE(6)
*----
*  COMPUTE THE GPT SOLUTION.
*----
      CALL DELDRV(IPTRK,IPSYS0,IPSYSP,IPFLU0,IPGPT,NUN,NGRP,NSTEP)
*----
*  RELEASE GENERAL TRACKING INFORMATION.
*----
      IF(JENTRY(1).EQ.0) THEN
         CALL XDISET(ISTATE,NSTATE,0)
         ISTATE(1)=NGRP
         ISTATE(2)=NUN
         CALL LCMLEN(IPGPT,'DSOUR',ILENG,ITYLCM)
         IF(ILENG.NE.0) ISTATE(3)=ILENG
         CALL LCMLEN(IPGPT,'ASOUR',ILENG,ITYLCM)
         IF(ILENG.NE.0) ISTATE(4)=ILENG
         CALL LCMPUT(IPGPT,'STATE-VECTOR',NSTATE,1,ISTATE)
      ENDIF
      RETURN
      END
