*DECK BIVACT
      SUBROUTINE BIVACT(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* BIVAC type (2-D) tracking operator.
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
*         HENTRY(1): create or modification type(L_TRACK);
*         HENTRY(2): read-only type(L_GEOM).
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
* The BIVACT: calling specifications are:
* TRACK := BIVACT: [ TRACK ] GEOM  :: (bivact\_data) ;
* where
*   TRACK : name of the \emph{lcm} object (type L\_BIVAC) containing the 
*     \emph{tracking} information. If TRACK appears on the RHS, the previous 
*     settings will be applied by default.
*   GEOM  : name of the \emph{lcm} object (type L\_GEOM) containing the 
*     geometry.
*   bivact\_data : structure containing the data to module BIVACT:}
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
      PARAMETER(NSTATE=40,IOUT=6)
      CHARACTER TEXT4*4,TEXT12*12,TITLE*72,HSIGN*12
      DOUBLE PRECISION DFLOTT
      LOGICAL LOG,LDIFF
      INTEGER IGP(NSTATE),ISTATE(NSTATE),NCODE(6)
*----
*  PARAMETER VALIDATION
*----
      IF(NENTRY.LE.1) CALL XABORT('BIVACT: TWO PARAMETERS EXPECTED.')
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2)) CALL XABORT('BIVACT: L'
     1 //'CM OBJECT EXPECTED AT LHS.')
      IF((JENTRY(1).NE.0).AND.(JENTRY(1).NE.1)) CALL XABORT('BIVACT: E'
     1 //'NTRY IN CREATE OR MODIFICATION MODE EXPECTED.')
      IF((JENTRY(2).NE.2).OR.((IENTRY(2).NE.1).AND.(IENTRY(2).NE.2)))
     1 CALL XABORT('BIVACT: LCM OBJECT IN READ-ONLY MODE EXPECTED AT R'
     2 //'HS.')
      CALL LCMGTC(KENTRY(2),'SIGNATURE',12,1,HSIGN)
      IF(HSIGN.NE.'L_GEOM') THEN
         TEXT12=HENTRY(2)
         CALL XABORT('BIVACT: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1   '. L_GEOM EXPECTED.')
      ENDIF
      HSIGN='L_TRACK'
      CALL LCMPTC(KENTRY(1),'SIGNATURE',12,1,HSIGN)
      HSIGN='BIVAC'
      CALL LCMPTC(KENTRY(1),'TRACK-TYPE',12,1,HSIGN)
      CALL LCMGET(KENTRY(2),'STATE-VECTOR',ISTATE)
      ITYPE=ISTATE(1)
      CALL LCMLEN(KENTRY(2),'BIHET',ILONG,ITYLCM)
      IF(ILONG.NE.0) CALL XABORT('BIVACT: DOUBLE-HETEROGENEITY NOT SUP'
     1 //'PORTED.')
*
      IMPX=1
      TITLE=' '
      IF(JENTRY(1).EQ.0) THEN
         MAXPTS=ISTATE(6)
         IELEM=1
         ICOL=2
         NLF=0
         ISPN=0
         ISCAT=0
         NVD=0
         CALL LCMGET(KENTRY(2),'NCODE',NCODE)
         LOG=.FALSE.
         DO 10 I=1,4
         LOG=LOG.OR.(NCODE(I).EQ.3)
   10    CONTINUE
         IF(LOG) MAXPTS=2*MAXPTS
      ELSE IF(JENTRY(1).EQ.1) THEN
         CALL LCMGTC(KENTRY(1),'SIGNATURE',12,1,HSIGN)
         IF(HSIGN.NE.'L_TRACK') THEN
            TEXT12=HENTRY(1)
            CALL XABORT('BIVACT: SIGNATURE OF '//TEXT12//' IS '//HSIGN
     1      //'. L_TRACK EXPECTED.')
         ENDIF
         CALL LCMGTC(KENTRY(1),'TRACK-TYPE',12,1,HSIGN)
         IF(HSIGN.NE.'BIVAC') THEN
            TEXT12=HENTRY(1)
            CALL XABORT('BIVACT: TRACK-TYPE OF '//TEXT12//' IS '//HSIGN
     1      //'. BIVAC EXPECTED.')
         ENDIF
         CALL LCMGET(KENTRY(1),'STATE-VECTOR',IGP)
         MAXPTS=IGP(1)
         IELEM=IGP(8)
         ICOL=IGP(9)
         NLF=IGP(14)
         ISPN=IGP(15)
         ISCAT=IGP(16)
         NVD=IGP(17)
         CALL LCMLEN(KENTRY(1),'TITLE',LENGT,ITYLCM)
         IF(LENGT.GT.0) CALL LCMGTC(KENTRY(1),'TITLE',72,1,TITLE)
      ENDIF
   15 CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
      IF(INDIC.EQ.10) GO TO 30
   20 IF(INDIC.NE.3) CALL XABORT('BIVACT: CHARACTER DATA EXPECTED.')
      IF(TEXT4.EQ.'EDIT') THEN
         CALL REDGET(INDIC,IMPX,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('BIVACT: INTEGER DATA EXPECTED(1).')
      ELSE IF(TEXT4.EQ.'TITL') THEN
         CALL REDGET(INDIC,NITMA,FLOTT,TITLE,DFLOTT)
         IF(INDIC.NE.3) CALL XABORT('BIVACT: TITLE EXPECTED.')
      ELSE IF(TEXT4.EQ.'MAXR') THEN
         CALL REDGET(INDIC,MAXPTS,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('BIVACT: INTEGER DATA EXPECTED(2).')
      ELSE IF(TEXT4.EQ.'PRIM') THEN
*        MESH CORNER FINITE DIFFERENCES OR PRIMAL FINITE ELEMENTS.
         IELEM=-1
         ICOL=2
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.EQ.1) THEN
            IELEM=-NITMA
            CALL REDGET(INDIC,ICOL,FLOTT,TEXT4,DFLOTT)
            IF(INDIC.NE.1) CALL XABORT('BIVACT: INTEGER DATA EXPECTED('
     1      //'3).')
         ELSE
            GO TO 20
         ENDIF
      ELSE IF(TEXT4.EQ.'MCFD') THEN
*        MESH CENTERED FINITE DIFFERENCES IN HEXAGONAL GEOMETRY.
         IF(ITYPE.NE.8) CALL XABORT('BIVACT: MCFD OPTION LIMITED TO HE'
     1   //'XAGONAL GEOMETRY.')
         ICOL=4
      ELSE IF(TEXT4.EQ.'DUAL') THEN
*        MESH CENTERED FINITE DIFFERENCES OR MIXED-DUAL FINITE ELEMENTS.
         IELEM=1
         ICOL=2
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.EQ.1) THEN
            IELEM=NITMA
            CALL REDGET(INDIC,ICOL,FLOTT,TEXT4,DFLOTT)
            IF(INDIC.NE.1) CALL XABORT('BIVACT: INTEGER DATA EXPECTED('
     1      //'6).')
         ELSE
            GO TO 20
         ENDIF
      ELSE IF(TEXT4.EQ.'VOID') THEN
         IF(NLF.EQ.0) CALL XABORT('BIVACT: SPN-RELATED OPTION.')
         CALL REDGET(INDIC,NVD,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('BIVACT: INTEGER DATA EXPECTED(8).')
         IF((NVD.LT.0).OR.(NVD.GT.2)) CALL XABORT('BIVACT: INVALID VAL'
     1   //'UE OF NVD (0, 1 OR 2 EXPECTED).')
      ELSE IF(TEXT4.EQ.'PN') THEN
         CALL REDGET(INDIC,NLF,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('BIVACT: INTEGER DATA EXPECTED(9).')
         IF(MOD(NLF,2).EQ.0) CALL XABORT('BIVACT: ODD PN ORDER EXPECT'
     1   //'ED.')
         IF(NLF.GT.0) NLF=NLF+1
         ISCAT=NLF
         ISPN=0
      ELSE IF(TEXT4.EQ.'SPN') THEN
         CALL REDGET(INDIC,NLF,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('BIVACT: INTEGER DATA EXPECTED(10)'
     1   //'.')
         IF(NLF.EQ.0) THEN
*           DIFFUSION THEORY.
            ISCAT=0
            ISPN=0
         ELSE
            IF(MOD(NLF,2).EQ.0) CALL XABORT('BIVACT: ODD SPN ORDER EXP'
     1      //'ECTED.')
            NLF=NLF+1
            ISCAT=NLF
            ISPN=1
         ENDIF
      ELSE IF(TEXT4.EQ.'SCAT') THEN
         IF(NLF.EQ.0) CALL XABORT('BIVACT: DEFINE PN OR SPN FIRST.')
         CALL REDGET(INDIC,ISCAT,FLOTT,TEXT4,DFLOTT)
         IF(ISCAT.LE.0) CALL XABORT('BIVACT: POSITIVE ISCAT EXPECTED.')
         LDIFF=.FALSE.
         IF((INDIC.EQ.3).AND.(TEXT4.EQ.'DIFF')) THEN
            LDIFF=.TRUE.
            CALL REDGET(INDIC,ISCAT,FLOTT,TEXT4,DFLOTT)
         ENDIF
         IF(INDIC.NE.1) CALL XABORT('BIVACT: INTEGER DATA EXPECTED(11)'
     1   //'.')
         IF(LDIFF) ISCAT=-ISCAT
      ELSE IF(TEXT4.EQ.';') THEN
         GO TO 30
      ELSE
         CALL XABORT('BIVACT: '//TEXT4//' IS AN INVALID KEY WORD.')
      ENDIF
      GO TO 15
*
   30 IF(TITLE.NE.' ') THEN
         CALL LCMPTC(KENTRY(1),'TITLE',72,1,TITLE)
      ENDIF
      IF((NLF.GT.0).AND.(IELEM.LT.0)) CALL XABORT('BIVACT: SPN APPROXI'
     1 //'MATIONS LIMITED TO DUAL DISCRETIZATIONS.')
      TEXT12=HENTRY(2)
      CALL LCMPTC(KENTRY(1),'LINK.GEOM',12,1,TEXT12)
      IF(IMPX.GT.1) WRITE(IOUT,100) TITLE
*
      IF(MAXPTS.EQ.0) CALL XABORT('BIVACT: MAXPTS NOT DEFINED.')
      CALL BIVTRK (MAXPTS,KENTRY(1),KENTRY(2),IMPX,IELEM,ICOL,NLF,NVD,
     1 ISPN,ISCAT)
*
      IF(IMPX.GT.1) THEN
         CALL LCMGET(KENTRY(1),'STATE-VECTOR',IGP)
         WRITE(IOUT,110) (IGP(I),I=1,17)
      ENDIF
      RETURN
*
  100 FORMAT(1H1,36HBBBBBB  IIIIII VV  VV   AA    CCCCC ,95(1H*)/
     1 38H BBBBBBB IIIIII VV  VV  AAAA  CCCCCCC ,56(1H*),
     2 38H MULTIGROUP VERSION.  A. HEBERT (1993)/
     3 37H BB   BB   II   VV  VV  AAAA  CC   CC/
     4 37H BBBBB     II   VV  VV AA  AA CC     /
     5 37H BBBBB     II   VV  VV AAAAAA CC     /
     6 37H BB   BB   II   VV  VV AAAAAA CC   CC/
     7 37H BBBBBBB IIIIII  VVVV  AA  AA CCCCCCC/
     8 37H BBBBBB  IIIIII   VV   AA  AA  CCCCC //1X,A72//)
  110 FORMAT(/14H STATE VECTOR:/
     1 7H NREG  ,I6,22H   (NUMBER OF REGIONS)/
     2 7H NUN   ,I6,23H   (NUMBER OF UNKNOWNS)/
     3 7H ILK   ,I6,39H   (0=LEAKAGE PRESENT/1=LEAKAGE ABSENT)/
     4 7H NBMIX ,I6,36H   (MAXIMUM NUMBER OF MIXTURES USED)/
     5 7H NSURF ,I6,29H   (NUMBER OF OUTER SURFACES)/
     6 7H ITYPE ,I6,21H   (TYPE OF GEOMETRY)/
     7 7H IHEX  ,I6,31H   (TYPE OF HEXAGONAL SYMMETRY)/
     8 7H IELEM ,I6,28H   (TYPE OF FINITE ELEMENTS)/
     9 7H ICOL  ,I6,47H   (TYPE OF QUADRATURE USED TO INTEGRATE THE MA,
     1 10HSS MATRIX)/
     2 7H ISPLH ,I6,37H   (TYPE OF HEXAGONAL MESH-SPLITTING)/
     3 7H LL4   ,I6,45H   (ORDER OF THE MATRICES PER GROUP IN BIVAC)/
     4 7H LX    ,I6,40H   (NUMBER OF ELEMENTS ALONG THE X AXIS)/
     5 7H LY    ,I6,40H   (NUMBER OF ELEMENTS ALONG THE Y AXIS)/
     6 7H NLF   ,I6,45H   (0=DIFFUSION/NB OF PN ORDERS FOR THE FLUX)/
     7 7H ISPN  ,I6,34H   (0=COMPLETE PN/1=SIMPLIFIED PN)/
     8 7H ISCAT ,I6,47H   (1=ISOTROPIC SOURCE/2=LINEARLY ANISOTROPIC S,
     9 6HOURCE)/
     1 7H NVD   ,I6,47H   (0=PN-TYPE VOID/1=SN-TYPE VOID/2=DIFFUSION-T,
     2 9HYPE VOID))
      END
