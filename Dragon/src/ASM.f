*DECK ASM
      SUBROUTINE ASM(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Multigroup assembly operator for system matrices.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): A. Hebert
*
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1): create or modification type(L_PIJ);
*         HENTRY(2): read-only type(L_MACROLIB or L_LIBRARY);
*         HENTRY(3): read-only type(L_TRACK);
*         HENTRY(4): optional read-only sequential binary tracking file.
* IENTRY  type of each LCM object or file:
*         =1 LCM memory object; =2 XSM file; =3 sequential binary file;
*         =4 sequential ascii file.
* JENTRY  access of each LCM object or file:
*         =0 the LCM object or file is created;
*         =1 the LCM object or file is open for modifications;
*         =2 the LCM object or file is open in read-only mode.
* KENTRY  LCM object address or file unit number.
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
      PARAMETER   (NSTATE=40,IOUT=6)
      CHARACTER    TEXT4*4,HSIGN*12,TEXT12*12,HSMG*131,CDOOR*12,
     1             TITRE*72
      DOUBLE PRECISION DFLOTT
      LOGICAL      LEAKSW,LNORM,LALBS,LDIFF,LADJ
      INTEGER      IGP(NSTATE),IPAR(NSTATE),IPP(NSTATE),NALBP
      TYPE(C_PTR)  IPSYS,IPTRK,IPMACR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MAT,IDL
      REAL, ALLOCATABLE, DIMENSION(:) :: VOL
*----
*  PARAMETER VALIDATION
*----
      IF(NENTRY.LE.2) CALL XABORT('ASM: THREE PARAMETERS EXPECTED.')
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2)) CALL XABORT('ASM: LC'
     1 //'M OBJECT EXPECTED AT LHS.')
      IF((JENTRY(1).NE.0).AND.(JENTRY(1).NE.1)) CALL XABORT('ASM: EN'
     1 //'TRY IN CREATE OR MODIFICATION MODE EXPECTED.')
      IF((JENTRY(2).NE.2).OR.((IENTRY(2).NE.1).AND.(IENTRY(2).NE.2)))
     1 CALL XABORT('ASM: LCM OBJECT IN READ-ONLY MODE EXPECTED AT FI'
     2 //'RST RHS.')
      IF((JENTRY(3).NE.2).OR.((IENTRY(3).NE.1).AND.(IENTRY(3).NE.2)))
     1 CALL XABORT('ASM: LCM OBJECT IN READ-ONLY MODE EXPECTED AT SE'
     2 //'COND RHS.')
      CALL LCMGTC(KENTRY(3),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_TRACK') THEN
         TEXT12=HENTRY(3)
         CALL XABORT('ASM: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1   '. L_TRACK EXPECTED.')
      ENDIF
      CALL LCMGTC(KENTRY(3),'TRACK-TYPE',12,1,CDOOR)
      HSIGN='L_PIJ'
      IPSYS=KENTRY(1)
      IPMACR=KENTRY(2)
      IPTRK=KENTRY(3)
      CALL LCMPTC(IPSYS,'SIGNATURE',12,1,HSIGN)
      TEXT12=HENTRY(2)
      CALL LCMPTC(IPSYS,'LINK.MACRO',12,1,TEXT12)
      TEXT12=HENTRY(3)
      CALL LCMPTC(IPSYS,'LINK.TRACK',12,1,TEXT12)
*----
*  RECOVER TABULATED FUNCTIONS
*----
      CALL XDRTA2(IPTRK)
*----
*  RECOVER TRACKING FILE INFORMATION
*----
      IF(NENTRY.LT.4) THEN
        IFTRAK=0
      ELSE
        TEXT12=HENTRY(4)
        IF(IENTRY(4).EQ.3) THEN
          IF(JENTRY(4).NE.2) CALL XABORT('ASM: BINARY TRACKING FILE NA'
     1    //'MED '//TEXT12//' IS NOT IN REAL-ONLY MODE.')
          IFTRAK=FILUNIT(KENTRY(4))
        ENDIF
      ENDIF
*----
*  RECOVER GENERAL TRACKING INFORMATION
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',IGP)
      NREG=IGP(1)
      NUN=IGP(2)
      LEAKSW=IGP(3).EQ.0
      IBFP=0
      IF(CDOOR.EQ.'MCCG') THEN
*        SET ANISOTROPY LEVEL FOR WITHIN-GROUP SCATTERING XS.
         NANI=IGP(6)
      ELSE IF((CDOOR.EQ.'BIVAC').OR.(CDOOR.EQ.'SN')) THEN
*        SET ANISOTROPY LEVEL FOR TOTAL AND WITHIN-GROUP SCATTERING XS.
         NANI=MAX(1,IGP(16))
         IF(CDOOR.EQ.'SN') IBFP=IGP(31)
      ELSE IF(CDOOR.EQ.'TRIVAC') THEN
*        SET ANISOTROPY LEVEL FOR TOTAL AND WITHIN-GROUP SCATTERING XS.
         NANI=MAX(1,IGP(32))
      ELSE
         NANI=1
      ENDIF
      IF(CDOOR.EQ.'MCCG') THEN
         CALL LCMLEN(IPTRK,'KEYFLX',LKFL,ITYLCM)
         NFUNL=LKFL/NREG
      ELSE
         NFUNL=1
      ENDIF
      ALLOCATE(MAT(NREG),VOL(NREG),IDL(NREG*NFUNL))
      CALL LCMLEN(IPTRK,'MATCOD',ILNLCM,ITYLCM)
      IF(ILNLCM.NE.NREG) THEN
         CALL XABORT( 'ASM: INCOMPATIBLE NUMBER OF REGIONS')
      ENDIF
      CALL LCMGET(IPTRK,'MATCOD',MAT)
      CALL LCMGET(IPTRK,'VOLUME',VOL)
      CALL LCMGET(IPTRK,'KEYFLX',IDL)
      CALL LCMLEN(IPTRK,'TITLE',LENGT,ITYLCM)
      IF(LENGT.GT.0) THEN
         CALL LCMGTC(IPTRK,'TITLE',72,1,TITRE)
      ELSE
         TITRE='*** NO TITLE PROVIDED ***'
      ENDIF
*----
*  RECOVER MACROLIB PARAMETERS
*----
      CALL LCMGTC(IPMACR,'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.EQ.'L_LIBRARY') THEN
         CALL LCMSIX(IPMACR,'MACROLIB',1)
      ELSE IF(HSIGN.NE.'L_MACROLIB') THEN
         TEXT12=HENTRY(2)
         CALL XABORT('ASM: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1   '. L_MACROLIB EXPECTED.')
      ENDIF
      CALL LCMGET(IPMACR,'STATE-VECTOR',IPAR)
      NGRP=IPAR(1)
      MAXMIX=IPAR(2)
      ITRANC=IPAR(6)
      NALBP=IPAR(8)
      LDIFF=IPAR(9).EQ.1
      NW=IPAR(10)
      LADJ=IPAR(13).EQ.1
      IF(IGP(4).GT.MAXMIX) THEN
         WRITE(HSMG,'(45HASM: THE NUMBER OF MIXTURES IN THE TRACKING (,
     1   I5,55H) IS GREATER THAN THE NUMBER OF MIXTURES IN THE MACROLI,
     2   3HB (,I5,2H).)') IGP(4),MAXMIX
         CALL XABORT(HSMG)
      ENDIF
      IF(CDOOR.EQ.'MCCG') THEN
         IF((IGP(38).EQ.1).AND.(.NOT.LADJ)) THEN
            CALL XABORT('ASM: ADJOINT MCCG CALCULATION IS EXPECTING AN'
     1      //' ADJOINT MACROLIB.')
         ENDIF
      ENDIF
*
      ITPIJ=1
      LNORM=.FALSE.
      LALBS=.FALSE.
      IPHASE=2
      ISTRM=1
      KNORM=4
      IF(JENTRY(1).EQ.1) THEN
         CALL LCMGTC(KENTRY(1),'SIGNATURE',12,1,HSIGN)
         IF(HSIGN.NE.'L_PIJ') THEN
            TEXT12=HENTRY(1)
            CALL XABORT('ASM: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1      '. L_PIJ EXPECTED.')
         ENDIF
         CALL LCMGET(KENTRY(1),'STATE-VECTOR',IPP)
         IF((IPP(8).NE.NGRP).OR.(IPP(9).NE.NUN)) THEN
            WRITE(HSMG,'(36HASM: INCONSISTENT NUMBER OF GROUPS (,I3,
     1      3H VS,I4,15H) OR UNKNOWNS (,I5,3H VS,I8,2H).)') IPP(8),
     2      NGRP,IPP(9),NUN
            CALL XABORT(HSMG)
         ENDIF
         ITPIJ=IPP(1)
         LNORM=IPP(2).EQ.0
         LALBS=IPP(3).EQ.0
         IPHASE=IPP(5)
         ISTRM=IPP(6)
         KNORM=IPP(7)
      ELSE IF(JENTRY(1).NE.0) THEN
         CALL XABORT('ASM: NO LHS OBJECT.')
      ENDIF
      IMPX=1
      NANIST=NANI
   15 CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
      IF(INDIC.EQ.10) GO TO 40
   20 IF(INDIC.NE.3) CALL XABORT('ASM: CHARACTER DATA EXPECTED(1).')
      IF(TEXT4.EQ.'EDIT') THEN
         CALL REDGET(INDIC,IMPX,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('ASM: INTEGER DATA EXPECTED(1).')
      ELSE IF(TEXT4.EQ.'PNOR') THEN
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.3) CALL XABORT('ASM: CHARACTER DATA EXPECTED(2).')
         IF(TEXT4.EQ.'NONE') THEN
           KNORM=0
         ELSE IF(TEXT4.EQ.'GELB') THEN
           KNORM=1
         ELSE IF(TEXT4.EQ.'DIAG') THEN
           KNORM=2
         ELSE IF(TEXT4.EQ.'NONL') THEN
           KNORM=3
         ELSE IF(TEXT4.EQ.'HELI') THEN
           KNORM=4
         ELSE
           GO TO 20
         ENDIF
      ELSE IF(TEXT4.EQ.'ARM') THEN
         IPHASE=1
      ELSE IF(TEXT4(1:3).EQ.'PIJ') THEN
         IPHASE=2
         IF(TEXT4(4:4).EQ.'K') THEN
           IF(CDOOR.EQ.'EXCELL') THEN
             ISTRM=3
             ITPIJ=ITPIJ+2
             NANI=MAX(2,NANI)
             NANIST=NANI
           ELSE
             WRITE(IOUT,6300) CDOOR
           ENDIF
         ENDIF
  30     CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.3) CALL XABORT('ASM: CHARACTER DATA EXPECTED(3).')
         IF(TEXT4.EQ.'NORM') THEN
           LNORM=.TRUE.
         ELSE IF(TEXT4.EQ.'ALBS') THEN
           LALBS=.TRUE.
           IF(.NOT.LEAKSW) THEN
            CALL XABORT('ASM: INVALID BOUNDARY CONDITIONS. THE ALSB '
     1      //'OPTION REQUIRES SOME BOUNDARY LEAKAGE.')
           ENDIF
         ELSE
           GO TO 20
         ENDIF
         GO TO 30
      ELSE IF(TEXT4.EQ.'SKIP') THEN
         ITPIJ=ITPIJ+1
      ELSE IF(TEXT4.EQ.'ECCO') THEN
         ISTRM=2
         NANIST=MAX(2,NANI)
      ELSE IF(TEXT4.EQ.'HETE') THEN
         ISTRM=3
         NANIST=MAX(2,NANI)   
      ELSE IF(TEXT4.EQ.';') THEN
         GO TO 40
      ELSE
         CALL XABORT('ASM: '//TEXT4//' IS AN INVALID KEY WORD.')
      ENDIF
      GO TO 15
*----
*  CHECK FOR THE ANISOTROPY SETTINGS COHERENCE
*----
   40 IF((ITRANC.NE.0).AND.(NANI.GT.1).AND.(ISTRM.EQ.1)) THEN
         WRITE(IOUT,6400) CDOOR,NANI
         ITRANC=0
      ENDIF
      IF((IMPX.GE.1).AND.(NW.GT.0)) THEN
         WRITE (IOUT,'(/44H ASM: A LEAKAGE CORRECTION IS PERFORMED (NW=,
     1   I2,2H).)') NW
      ENDIF
*----
*  STORE PIJ PARAMETERS
*----
      CALL XDISET(IPP,NSTATE,0)
      IPP(1)=ITPIJ
      IPP(2)=1
      IF(LNORM) IPP(2)=0
      IPP(3)=1
      IF(LALBS) IPP(3)=0
      IPP(5)=IPHASE
      IPP(6)=ISTRM
      IPP(7)=KNORM
      IPP(8)=NGRP
      IPP(9)=NUN
      IPP(10)=MAXMIX
      IPP(11)=NANI
      IF(LDIFF) IPP(12)=1
      IPP(13)=IBFP
      CALL LCMPUT(IPSYS,'STATE-VECTOR',NSTATE,1,IPP)
*----
*  BUILD COLLISION PROBABILITIES
*----
      CALL ASMDRV(IPSYS,IPTRK,IPMACR,IFTRAK,CDOOR,IMPX,NGRP,MAXMIX,
     1 NREG,NANI,NANIST,NW,MAT,VOL,LEAKSW,ITRANC,LDIFF,IBFP,TITRE,ITPIJ,
     2 LNORM.OR.LALBS,IPHASE,ISTRM,KNORM,NALBP)
*
      IF(IMPX.GE.5) CALL LCMLIB(IPSYS)
*----
*  RELEASE GENERAL TRACKING INFORMATION
*----
      DEALLOCATE(IDL,VOL,MAT)
      CALL LCMSIX(IPMACR,' ',0)
      IF(IMPX.GE.1) THEN
         WRITE (IOUT,6040)IMPX,(IPP(I),I=1,3),(IPP(I),I=5,13)
         WRITE (IOUT,'(5H DOOR,13X,1H(,A,1H))') CDOOR
      ENDIF
      RETURN
*
 6300 FORMAT(//' *** WARNING:  OPTION PIJK IS INVALID FOR DOOR = ',
     >          A12/'               OPTION PIJ USED INSTEAD')
 6400 FORMAT(//' *** WARNING:  DOOR ',A12,'IS USED WITH AN ANISOTROPY',
     > ' LEVEL FROM L_TRACK =',I2,' AND WITH A TRANSPORT CORRECTION S',
     > 'ET IN LIB:.'/15X,'--> THE TRANSPORT CORRECTION IS DISABLED.'/)
 6040 FORMAT(/8H OPTIONS/8H -------/
     1 7H IMPX  ,I8,30H   (0=NO PRINT/1=SHORT/2=MORE)/
     2 7H ITPIJ ,I8,30H   (1=WIJ/2=PIJ/3=WIJK/4=PIJK)/
     3 7H LNORM ,I8,34H   (0=NORMALIZE PIJ TO 1/1=DO NOT)/
     4 7H LALBS ,I8,36H   (0=RECOVER AND SAVE WIS/1=DO NOT)/
     5 7H IPHASE,I8,43H   (1=GENERAL FLUX SOLUTION/2=PIJ APPROACH)/
     6 7H ISTRM ,I8,44H   (1=HOMO BN OR NO LEAKAGE/2=ECCO/3=TIBERE)/
     7 7H KNORM ,I8,46H   (0=NO/1=GELBARD/2=DIAGONAL/3=NON-LINEAR/4=H,
     8 6HELIOS)/
     9 7H NGRP  ,I8,21H   (NUMBER OF GROUPS)/
     1 7H NUN   ,I8,23H   (NUMBER OF UNKNOWNS)/
     2 7H NBMIX ,I8,23H   (NUMBER OF MIXTURES)/
     3 7H NANI  ,I8,44H   (NUMBER OF LEGENDRE ORDERS SCATTERING XS)/
     4 7H IDIFF ,I8,47H   (0/1: DIFFUSION COEFFICIENTS ABSENT/PRESENT)/
     5 7H IBFP  ,I8,44H   (0/1/2: FOKKER-PLANCK SOLUTION OFF/ON/ON))
      END
