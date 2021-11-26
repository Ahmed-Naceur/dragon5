*DECK GPTFLU
      SUBROUTINE GPTFLU(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Computes generalized adjoints.
* GPTFLU = Generalized Perturbation Theory FLUx calculation
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert, E. Varin and R. Chambon
*
*Parameters: input/ouput
* NENTRY  number of linked lists or files used by the module
* HENTRY  character*12 name of each linked list or file
*         HENTRY(1): create or modification type(L_FLUX) (GPT solution)
*         HENTRY(2): read-only type(L_SOURCE) => GPT fixed source
*         HENTRY(3): read-only type(L_FLUX) => unperturbed solution
*         HENTRY(4): read-only type(L_SYSTEM) => reference matrices
*         HENTRY(5): read-only type(L_TRACK) => TRIVAC tracking.
* IENTRY  =1 linked list; =2 xsm file; =3 sequential binary file;
*         =4 sequential ascii file
* JENTRY  =0 the linked list or file is created
*         =1 the linked list or file is open for modifications;
*         =2 the linked list or file is open in read-only mode
* KENTRY  =file unit number; =linked list address otherwise.
*
*Comments:
* The GPTFLU: calling specifications are:
* FLUX\_GPT := GPTFLU: [ FLUX\_GPT ] GPT FLUX0 SYST TRACK :: (gptflu\_data) ;
* where
*   FLUX\_GPT : name of the \emph{lcm} object (type L\_FLUX) containing the GPT 
*     solution. If FLUX\_GPT} appears on the RHS, the solution previously stored
*     in FLUX\_GPT} is used to initialize the new iterative process; otherwise, 
*     a uniform unknown vector is used.
*   GPT       : name of the \emph{lcm} object (type L\_GPT) containing the 
*     fixed sources.
*   FLUX0     : name of the \emph{lcm} object (type L\_FLUX) containing the 
*     unperturbed flux used to decontaminate the GPT solution.
*   SYST      : name of the \emph{lcm} object (type L\_SYSTEM) containing the 
*     unperturbed system matrices.
*   TRACK     : name of the \emph{lcm} object (type L\_TRACK) containing the 
*     \emph{tracking}.
*   gptflu\_data}] : structure containing the data to module GPTFLU:}
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT  NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      CHARACTER HENTRY(NENTRY)*12
      TYPE(C_PTR) KENTRY(NENTRY)
*----
*  LOCAL VARIABLES
*----
      INTEGER   IPRINT,NITMA,I,J,IGR,ITYP,LENGT,SIGNA(3),ITY,NSTART
      DOUBLE PRECISION DFLOTT
      REAL      FLOTT
      CHARACTER TEXT12*12,CMODUL*12
      LOGICAL   LFLU
      INTEGER   NEL,NUN,NGRP,LL4,ISRC
      CHARACTER TITLE*72
*----
*  STATE-VECTOR VARIABLES
*----
      INTEGER   NSTATE
      PARAMETER (NSTATE=40)
      INTEGER   FLUPRM(NSTATE),SYSPRM(NSTATE),TRKPRM(NSTATE),
     1          GPTPRM(NSTATE)
*----
*  Generalized Adjoint calculation
*----
      INTEGER   SRCFRM,SRCTO,MAXOUT,ICL1,ICL2,NADI,IMPH,NLF,MAXINR
      REAL      FKEFF,FKEFF2,EPSOUT,EPSINR,EPSCON(5)
      LOGICAL   ADJ,REC,RECP
      INTEGER, DIMENSION(:), ALLOCATABLE :: MAT,IDL
      REAL, DIMENSION(:), ALLOCATABLE :: VOL
      REAL, DIMENSION(:,:), ALLOCATABLE :: EVECT,ADECT,EASS,SOUR
      TYPE(C_PTR) IPFLUP,IPFLU,IPGPT,IPTRK,IPSYS,JPFLU1,JPFLU2,JPGPT,
     1          KPGPT,JPFLUP,KPFLUP
*----
*  VALIDITY OF OBJECTS
*----
      IF(NENTRY.LT.4) CALL XABORT('GPTFLU: 5 OBJECTS EXPECTED.')
      IPFLUP=C_NULL_PTR
      IPFLU =C_NULL_PTR
      IPGPT =C_NULL_PTR
      IPTRK =C_NULL_PTR
      IPSYS =C_NULL_PTR
      RECP=(JENTRY(1).EQ.1)
      LFLU=.FALSE.
      DO 2 I=1,NENTRY
         TEXT12=HENTRY(I)
         IF((IENTRY(I).NE.1).AND.(IENTRY(I).NE.2))CALL XABORT('GPTFLU:'
     1      //' LINKED LIST OR XSM FILE EXPECTED ('//TEXT12//')')
         IF((JENTRY(I).EQ.0).AND.(I.EQ.1)) THEN
            TEXT12='L_FLUX'
            READ(TEXT12,'(3A4)') (SIGNA(J),J=1,3)
            CALL LCMPUT(KENTRY(I),'SIGNATURE',3,3,SIGNA)
         ELSE
            CALL LCMGET(KENTRY(I),'SIGNATURE',SIGNA)
            WRITE(TEXT12,'(3A4)') (SIGNA(J),J=1,3)
            IF(JENTRY(I).EQ.0)  CALL XABORT('GPTFLU:'//TEXT12//' IS '
     1      //'NOT ON RHS')
         ENDIF
         IF(TEXT12.EQ.'L_FLUX') THEN 
             IF(LFLU) THEN
                IPFLU=KENTRY(I)
             ELSE
                IPFLUP=KENTRY(I)
                LFLU=.TRUE.
             ENDIF
         ELSEIF(TEXT12.EQ.'L_SOURCE') THEN 
             IPGPT=KENTRY(I)
         ELSEIF(TEXT12.EQ.'L_TRACK') THEN 
             IPTRK=KENTRY(I)
         ELSEIF(TEXT12.EQ.'L_SYSTEM') THEN 
             IPSYS=KENTRY(I)
         ELSE
             CALL XABORT('GPTFLU: NOT GOOD TYPE OF OBJECT')
         ENDIF
    2 CONTINUE
      IF(.NOT.C_ASSOCIATED(IPGPT))
     1 CALL XABORT('GPTFLU: MISSING GPT SOURCE OBJECT.')
      IF(.NOT.C_ASSOCIATED(IPFLU))
     1 CALL XABORT('GPTFLU: MISSING FLUX OBJECT.')
      IF(.NOT.C_ASSOCIATED(IPSYS))
     1 CALL XABORT('GPTFLU: MISSING SYSTEM OBJECT.')
      IF(.NOT.C_ASSOCIATED(IPTRK))
     1 CALL XABORT('GPTFLU: MISSING TRACK OBJECT.')
*----
*  VARIABLE INITIALISATION	
*----
      CALL LCMGET(IPFLU,'STATE-VECTOR',FLUPRM)
      NGRP   = FLUPRM(1)
      NUN    = FLUPRM(2)
      MAXINR = FLUPRM(11)
      MAXOUT = FLUPRM(12)
      CALL LCMGET(IPFLU,'EPS-CONVERGE',EPSCON)
      EPSINR=EPSCON(1)
      EPSOUT=EPSCON(2)
      CALL LCMGET(IPTRK,'STATE-VECTOR',TRKPRM)
      NEL    = TRKPRM(1)
      IF(NUN.NE.TRKPRM(2)) CALL XABORT('GPTFLU: TRACKING AND UNPERTURB'
     + //'ED FLUX HAVE DIFFERENT NUMBER OF UNKNOWS')
      CALL LCMLEN(IPTRK,'TITLE',LENGT,ITYP)
      IF(LENGT.GT.0) THEN
         CALL LCMGTC(IPTRK,'TITLE',72,1,TITLE)
      ELSE
         TITLE='*** NO TITLE PROVIDED ***'
      ENDIF
      CALL LCMGTC(IPTRK,'TRACK-TYPE',12,1,CMODUL)
      IF(CMODUL.NE.'TRIVAC') CALL XABORT('GPTFLU: TRIVAC TRACKING EXPE'
     + //'CTED.')
      LL4 = TRKPRM(11)
      NLF = TRKPRM(30)
      CALL LCMGET(IPSYS,'STATE-VECTOR',SYSPRM)
      IF( SYSPRM(1).NE.NGRP )CALL XABORT('GPTFLU: L_SYSTEM AND L_FLUX'
     + //'FOR UNPERTURBED STATE HAVE DIFFERENT NUMBER OF GROUPS')
      IF( SYSPRM(2).NE.LL4 )CALL XABORT('GPTFLU: UNPERTURBED SYSTEM A'
     + //'ND TRACKING OBJECTS HAVE DIFFERENT NUMBER OF LINEAR ORDER')
      ITY = SYSPRM(4)
      IF(ITY.EQ.13) LL4=LL4*NLF/2
      CALL LCMGET(IPGPT,'STATE-VECTOR',GPTPRM)
      IF( GPTPRM(1).NE.NGRP )CALL XABORT('GPTFLU: L_SOURCE AND L_FLUX '
     + //'HAVE DIFFERENT NUMBER OF GROUPS')
      IF( GPTPRM(2).NE.NUN )CALL XABORT('GPTFLU: L_SOURCE AND L_FLUX H'
     + //'AVE DIFFERENT NUMBER OF UNKNOWS')
*----
*  READ USER INPUT:
*----
      IPRINT=0
      IMPH=0
      IF(RECP) THEN
*        RECOVER EXISTING OPTIONS.
         CALL LCMGET(IPFLU,'STATE-VECTOR',FLUPRM)
         ICL1=FLUPRM(8)
         ICL2=FLUPRM(9)
         MAXINR=FLUPRM(11)
         MAXOUT=FLUPRM(12)
         NADI=FLUPRM(13)
         NSTART=FLUPRM(16)
         CALL LCMGET(IPFLU,'EPS-CONVERGE',EPSCON)
         EPSINR=EPSCON(1)
         EPSOUT=EPSCON(2)
      ELSE
*        DEFAULT OPTIONS.
         ICL1=3
         ICL2=3
         NADI=TRKPRM(33)
         MAXINR=0
         MAXOUT=200
         NSTART=0
         EPSINR=1.0E-2
         EPSOUT=1.0E-4
      ENDIF
      IF(GPTPRM(3).LE.GPTPRM(4)) THEN
         ADJ=.FALSE.
      ELSE
         ADJ=.TRUE.
      ENDIF
      REC=.FALSE.
  505 CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
  506 IF(ITYP.NE.3) CALL XABORT('GPTFLU: CHARACTER DATA EXPECTED.')
      IF(TEXT12.EQ.'EDIT') THEN
         CALL REDGET(ITYP,IPRINT,FLOTT,TEXT12,DFLOTT)
         IF(ITYP.NE.1) CALL XABORT('GPTFLU: *IPRINT* MUST BE INTEGER')
         GO TO 505
      ELSEIF((TEXT12.EQ.'VAR1').OR.(TEXT12.EQ.'ACCE')) THEN
         CALL REDGET(ITYP,ICL1,FLOTT,TEXT12,DFLOTT)
         IF(ITYP.NE.1) CALL XABORT('GPTFLU: INTEGER DATA EXPECTED.')
         CALL REDGET(ITYP,ICL2,FLOTT,TEXT12,DFLOTT)
         IF(ITYP.NE.1) CALL XABORT('GPTFLU: INTEGER DATA EXPECTED.')
         GO TO 505
      ELSEIF(TEXT12.EQ.'GMRES') THEN
         CALL REDGET(ITYP,NSTART,FLOTT,TEXT12,DFLOTT)
         IF(ITYP.NE.1) CALL XABORT('GPTFLU: INTEGER DATA EXPECTED.')
         IF(NSTART.LT.0) CALL XABORT('GPTFLU: POSITIVE VALUE EXPECTED.')
         GO TO 505
      ELSEIF(TEXT12.EQ.'IMPLICIT') THEN
         ADJ=.TRUE.
         GO TO 505
      ELSEIF(TEXT12.EQ.'EXPLICIT') THEN
         ADJ=.FALSE.
         GO TO 505
      ELSEIF(TEXT12.EQ.'RCVR-LAST') THEN
         REC=.TRUE.
         GO TO 505
      ELSEIF(TEXT12.EQ.'EXTE') THEN
  507    CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
         IF(ITYP.EQ.1) THEN
            MAXOUT=NITMA
         ELSE IF(ITYP.EQ.2) THEN
            EPSOUT=FLOTT
         ELSE
            GO TO 506
         ENDIF
         GO TO 507
      ELSEIF(TEXT12.EQ.'ADI') THEN
         CALL REDGET(ITYP,NADI,FLOTT,TEXT12,DFLOTT)
         IF(ITYP.NE.1) CALL XABORT('GPTFLU: INTEGER DATA EXPECTED.')
         GO TO 505
      ELSEIF(TEXT12.EQ.'THER') THEN
         MAXINR = NGRP*2
  508    CALL REDGET(ITYP,NITMA,FLOTT,TEXT12,DFLOTT)
         IF(ITYP.EQ.1) THEN
            MAXINR=NITMA
         ELSEIF(ITYP.EQ.2) THEN
            EPSINR=FLOTT
         ELSE
            GO TO 506
         ENDIF 
         GO TO 508
      ELSEIF(TEXT12.EQ.'HIST') THEN
         CALL REDGET(ITYP,IMPH,FLOTT,TEXT12,DFLOTT)
         IF(ITYP.NE.1) CALL XABORT('GPTFLU: INTEGER DATA EXPECTED.')
         GO TO 505
      ENDIF
      IF(TEXT12.EQ.'FROM-TO') THEN
         CALL REDGET(ITYP,SRCFRM,FLOTT,TEXT12,DFLOTT)
         IF((ITYP.EQ.3).AND.(TEXT12.EQ.'ALL')) THEN
            SRCFRM=1
            IF(ADJ) THEN
               SRCTO=GPTPRM(4)
            ELSE
               SRCTO=GPTPRM(3)
            ENDIF
         ELSE
            IF(ITYP.NE.1) CALL XABORT('GPTFLU: INTEGER DATA EXPECTED.')
            CALL REDGET(ITYP,SRCTO,FLOTT,TEXT12,DFLOTT)
            IF(ITYP.NE.1) CALL XABORT('GPTFLU: INTEGER DATA EXPECTED.')
            IF(ADJ) THEN
               IF(SRCTO.GT.GPTPRM(4)) WRITE(6,*) 'THE NUMBER OF THE '//
     1         'SOURCE ',SRCTO,' IS GREATER THAN THE NUMBER OF CONST'//
     2         'RAINTS +1',GPTPRM(4)
            ELSE
               IF(SRCTO.GT.GPTPRM(3)) WRITE(6,*) 'THE NUMBER OF THE '//
     1         'SOURCE ',SRCTO,' IS GREATER THAN THE NUMBER OF VARIA'//
     2         'BLES',GPTPRM(3)
            ENDIF
            IF(SRCFRM.GT.SRCTO) CALL XABORT('GPTFLU:SCRFRM .GT. SCRTO')
         ENDIF
      ELSEIF(TEXT12.EQ.';') THEN
         GO TO 1000 
      ELSE
         WRITE(6,*) 'Your keyword is : ',TEXT12
         CALL XABORT('GPTFLU:"FROM-TO" or ";" EXPECTED')
      ENDIF
*----
*  RECOVER TRACKING INFORMATION
*----
      ALLOCATE(MAT(NEL),VOL(NEL),IDL(NEL))
      CALL LCMGET(IPTRK,'MATCOD',MAT)
      CALL LCMGET(IPTRK,'VOLUME',VOL)
      CALL LCMGET(IPTRK,'KEYFLX',IDL)
*----
*  RECOVER UNPERTURBED K-EFFECTIVE AND FLUXES.
*----
      ALLOCATE(EVECT(NUN,NGRP),ADECT(NUN,NGRP),EASS(NUN,NGRP))
      CALL LCMGET(IPFLU,'K-EFFECTIVE',FKEFF)
      JPFLU1=LCMGID(IPFLU,'FLUX')
      JPFLU2=LCMGID(IPFLU,'AFLUX')
      DO 510 IGR=1,NGRP
      CALL LCMGDL(JPFLU1,IGR,EVECT(1,IGR))
      CALL LCMGDL(JPFLU2,IGR,ADECT(1,IGR))
  510 CONTINUE
      ALLOCATE(SOUR(NUN,NGRP))
*----
*  RECOVER FIXED SOURCE AND SET INITIAL VALUE OF GPFLUX 
*----
      IF(ADJ) THEN
         JPFLUP=LCMLID(IPFLUP,'ADFLUX',SRCTO)
      ELSE
         JPFLUP=LCMLID(IPFLUP,'DFLUX',SRCTO)
      ENDIF
      DO 590 ISRC=SRCFRM,SRCTO
      IF(ADJ) THEN
         JPGPT=LCMGID(IPGPT,'ASOUR')
      ELSE
         JPGPT=LCMGID(IPGPT,'DSOUR')
      ENDIF
      CALL LCMLEL(JPGPT,ISRC,LENGT,ITYP)
      IF(LENGT.EQ.0) GO TO 590
      KPGPT=LCMGIL(JPGPT,ISRC)
      DO 520 IGR=1,NGRP
      CALL LCMGDL(KPGPT,IGR,SOUR(1,IGR))
  520 CONTINUE
      IF(REC.AND.(IMPH.EQ.0)) THEN
         CALL LCMLEL(JPFLUP,ISRC,LENGT,ITYP)
         IF(LENGT.EQ.0) THEN
            WRITE(TEXT12,'(I4,3H-TH)') ISRC
            CALL XABORT('GPTFLU: '//TEXT12//' GENERALIZED ADJOINT CANN'
     1      //'OT BE RECOVERED.')
         ENDIF
         KPFLUP=LCMGIL(JPFLUP,ISRC)
         DO 530 IGR=1,NGRP
         CALL LCMGDL(KPFLUP,IGR,EASS(1,IGR))
  530    CONTINUE
      ELSE
         DO 540 IGR=1,NGRP
         CALL XDRSET(EASS(1,IGR),NUN,1.0)
  540    CONTINUE
      ENDIF
*----
*  ADJOINT NEUTRON FLUX CALCULATION
*----
      IF(IPRINT.GE.1) WRITE(6,*) 'GPTFLU: ISRC=',ISRC
      IF(ADJ) THEN
         IF(IPRINT.GE.2) WRITE(6,*) 'implicit'
         CALL GPTAFL(IPTRK,IPSYS,IPFLUP,LL4,ITY,NUN,NGRP,ICL1,ICL2,
     1   NSTART,IPRINT,IMPH,TITLE,EPSOUT,MAXINR,EPSINR,NADI,MAXOUT,
     2   FKEFF,EVECT,ADECT,FKEFF2,EASS,SOUR)
      ELSE
         IF(IPRINT.GE.2) WRITE(6,*) 'explicit'
         CALL GPTDFL(IPTRK,IPSYS,IPFLUP,LL4,ITY,NUN,NGRP,ICL1,ICL2,
     1   NSTART,IPRINT,IMPH,TITLE,EPSOUT,MAXINR,EPSINR,NADI,MAXOUT,
     2   FKEFF,EVECT,ADECT,FKEFF2,EASS,SOUR)
      ENDIF
      CALL LCMPUT(IPFLUP,'K-EFFECTIVE',1,2,FKEFF2)
      KPFLUP=LCMLIL(JPFLUP,ISRC,NGRP)
      DO 550 IGR=1,NGRP
      CALL FLDTRI(IPTRK,NEL,NUN,EASS(1,IGR),MAT,VOL,IDL)
      CALL LCMPDL(KPFLUP,IGR,NUN,2,EASS(1,IGR))
  550 CONTINUE
  590 CONTINUE
      DEALLOCATE(EASS,ADECT,EVECT,SOUR,IDL,VOL,MAT)
      GO TO 505
*----
*  END
*----
 1000 IF(ADJ) THEN
         FLUPRM(3)=1000
      ELSE
         FLUPRM(3)=100
      ENDIF
      FLUPRM(5)=SRCTO
      FLUPRM(6)=1
      FLUPRM(16)=NSTART
      IF(IPRINT.GT.0) WRITE(6,2020) (FLUPRM(I),I=1,5),FLUPRM(16)
      CALL LCMPUT(IPFLUP,'STATE-VECTOR',NSTATE,1,FLUPRM)
      RETURN
*
 2020 FORMAT(/8H OPTIONS/8H -------/
     1 7H NGRP  ,I8,28H   (NUMBER OF ENERGY GROUPS)/
     2 7H NUN   ,I8,40H   (NUMBER OF UNKNOWNS PER ENERGY GROUP)/
     3 7H IADJ  ,I8,30H   (=100/1000: DIRECT/ADJOINT)/
     4 7H NMOD  ,I8,13H   (NOT USED)/
     5 7H SRCTO ,I8,48H   (NUMBER OF FIXED-SOURCE EIGENVALUE EQUATIONS)/
     6 7H NSTART,I8,46H   (NUMBER OF GMRES ITERATIONS BEFORE RESTART))
      END
