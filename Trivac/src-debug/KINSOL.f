*DECK KINSOL
      SUBROUTINE KINSOL(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* solve the space-time neutron kinetics equations.
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal.
*
*Author(s): D. Sekki
*
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1): modification type(L_KINET);
*         HENTRY(2): read-only type(L_MACROLIB);
*         HENTRY(3): read-only type(L_TRACK);
*         HENTRY(4): read-only type(L_SYSTEM) made with HENTRY(2);
*         HENTRY(5): optional read-only type(L_MACROLIB);
*         HENTRY(6): optional read-only type(L_SYSTEM) made with
*                    HENTRY(5).
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
* The KINSOL: calling specifications are:
* KINET := KINSOL: KINET MACRO TRACK SYST [ MACRO\_0 SYST\_0 ] :: 
*                  (kinsol\_data) ; 
* where
*   KINET : name of the \emph{lcm} object (type L\_KINET) in modification mode.
*   MACRO : name of the \emph{lcm} object (type L\_MACROLIB) containing the 
*     \emph{macrolib} information corresponding to the current time step of a 
*     transient.
*   TRACK : name of the \emph{lcm} object (type L\_TRACK) containing the 
*     \emph{tracking} information.
*   SYST  : name of the \emph{lcm} object (type L\_SYSTEM) corresponding to 
*     \emph{macrolib} MACRO and \emph{tracking} TRACK.
*   MACRO\_0 : name of the \emph{lcm} object (type L\_MACROLIB) containing the 
*     \emph{macrolib} information corresponding to the beginning of step 
*     conditions in case a ramp variation of the cross sections in set. 
*     Beginning of step conditions should not be confused with beginning of
*     transient or initial conditions.} By default, a step variation is set 
*     where cross sections are assumed constant and given by MACRO.
*   SYST\_0  : name of the \emph{lcm} object (type L\_SYSTEM) corresponding to 
*     \emph{macrolib} MACRO\_0   and \emph{tracking} TRACK.
*   kinsol\_data : structure containing the data to module KINSOL: 
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
      CHARACTER TEXT12*12,HSIGN*12,CMODUL*12,HSMG*131
*----
*  PARAMETER VALIDATION
*----
      IF((NENTRY.NE.4).AND.(NENTRY.NE.6))CALL XABORT('@KINSOL:'
     1 //' INVALID NUMBER OF MODULE PARAMETERS.')
      DO 10 IEN=1,NENTRY
        IF((IENTRY(IEN).NE.1).AND.(IENTRY(IEN).NE.2))
     1  CALL XABORT('@KINSOL: LCM OBJECTS EXPECTED.')
   10 CONTINUE
*     L_KINET
      CALL LCMGTC(KENTRY(1),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_KINET')THEN
        TEXT12=HENTRY(1)
        CALL XABORT('@KINSOL: SIGNATURE OF '//TEXT12//' IS '
     1  //HSIGN//'. L_KINET EXPECTED.')
      ENDIF
      IF(JENTRY(1).NE.1)CALL XABORT('@KINSOL: L_KINET IN MODI'
     1 //'FICATION MODE EXPECTED.')
      CALL LCMGTC(KENTRY(1),'TRACK-TYPE',12,1,CMODUL)
*     L_MACROLIB(1)
      CALL LCMGTC(KENTRY(2),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_MACROLIB')THEN
        TEXT12=HENTRY(2)
        CALL XABORT('@KINSOL: SIGNATURE OF '//TEXT12//' IS '
     1  //HSIGN//'. L_MACROLIB EXPECTED(1).')
      ENDIF
      IF(JENTRY(2).NE.2)CALL XABORT('@KINSOL: L_MACROLIB IN R'
     1 //'EAD-ONLY MODE EXPECTED AT RHS(1).')
*     L_TRACK
      CALL LCMGTC(KENTRY(3),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_TRACK')THEN
        TEXT12=HENTRY(3)
        CALL XABORT('@KINSOL: SIGNATURE OF '//TEXT12//' IS '
     1  //HSIGN//'. L_TRACK EXPECTED.')
      ENDIF
      IF(JENTRY(3).NE.2)CALL XABORT('@KINSOL: L_TRACK IN READ'
     1 //'-ONLY MODE EXPECTED AT RHS.')
      CALL LCMGTC(KENTRY(3),'TRACK-TYPE',12,1,HSIGN)
      IF(HSIGN.NE.CMODUL)CALL XABORT('@KINSOL: INVALID TRACKI'
     1 //'NG TYPE IN L_TRACK.')
*     L_SYSTEM(1)
      CALL LCMGTC(KENTRY(4),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_SYSTEM')THEN
        TEXT12=HENTRY(4)
        CALL XABORT('@KINSOL: SIGNATURE OF '//TEXT12//' IS '
     1  //HSIGN//'. L_SYSTEM EXPECTED.')
      ENDIF
      IF(JENTRY(4).NE.2)CALL XABORT('@KINSOL: L_SYSTEM IN READ'
     1 //'-ONLY MODE EXPECTED AT RHS.')
      CALL LCMGTC(KENTRY(4),'LINK.MACRO',12,1,TEXT12)
      IF(HENTRY(2).NE.TEXT12) THEN
        WRITE(HSMG,'(40H@KINSOL: INVALID MACROLIB OBJECT NAME ='',
     1  A12,18H'', EXPECTED NAME='',A12,2H''.)') HENTRY(2),TEXT12
        CALL XABORT(HSMG)
      ENDIF
      CALL LCMGTC(KENTRY(4),'LINK.TRACK',12,1,TEXT12)
      IF(HENTRY(3).NE.TEXT12) THEN
        WRITE(HSMG,'(40H@KINSOL: INVALID TRACKING OBJECT NAME ='',A12,
     1  18H'', EXPECTED NAME='',A12,2H''.)') HENTRY(3),TEXT12
        CALL XABORT(HSMG)
      ENDIF
      CALL LCMPTC(KENTRY(1),'LINK.TRACK',12,1,TEXT12)
*     L_MACROLIB(2)
      IF(NENTRY.EQ.6)THEN
        CALL LCMGTC(KENTRY(5),'SIGNATURE',12,1,HSIGN)
        IF(HSIGN.NE.'L_MACROLIB')THEN
          TEXT12=HENTRY(5)
          CALL XABORT('@KINSOL: SIGNATURE OF '//TEXT12//' IS '
     1    //HSIGN//'. L_MACROLIB EXPECTED(2).')
        ENDIF
        IF(JENTRY(5).NE.2)CALL XABORT('@KINSOL: L_MACROLIB IN'
     1  //' READ-ONLY MODE EXPECTED AT RHS(2).')
*       L_SYSTEM(2)
        CALL LCMGTC(KENTRY(6),'SIGNATURE',12,1,HSIGN)
        IF(HSIGN.NE.'L_SYSTEM')THEN
          TEXT12=HENTRY(6)
          CALL XABORT('@KINSOL: SIGNATURE OF '//TEXT12//' IS '
     1    //HSIGN//'. L_SYSTEM EXPECTED.')
        ENDIF
        IF(JENTRY(6).NE.2)CALL XABORT('@KINSOL: L_SYSTEM IN READ'
     1   //'-ONLY MODE EXPECTED AT RHS.')
        CALL LCMGTC(KENTRY(6),'LINK.TRACK',12,1,TEXT12)
        IF(HENTRY(3).NE.TEXT12) THEN
          WRITE(HSMG,'(40H@KINSOL: INVALID TRACKING OBJECT NAME ='',
     1    A12,18H'', EXPECTED NAME='',A12,2H''.)') HENTRY(3),TEXT12
          CALL XABORT(HSMG)
        ENDIF
      ENDIF
*----
*  READ THE INPUT DATA
*----
      CALL KINRD2(NENTRY,KENTRY,CMODUL)
      RETURN
      END
