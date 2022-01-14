*DECK LIB
      SUBROUTINE LIB(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Interpolation of nuclear properties in an internal library.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): A. Hebert and G. Marleau
*
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1): create or modification type(L_LIBRARY)
*         HENTRY(2): optional read-only type(L_MACROLIB) used to
*                    initialize a new lattice code library.
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
      USE          GANLIB
      IMPLICIT     NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER      NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      TYPE(C_PTR)  KENTRY(NENTRY)
      CHARACTER    HENTRY(NENTRY)*12
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR)  IPLIB,IPLIBX
      INTEGER      IOUT,NSTATE,ILCMUP,ILCMDN,MAXED,MAXISD
      CHARACTER    NAMSBR*6
      PARAMETER   (IOUT=6,NSTATE=40,ILCMUP=1,ILCMDN=2,MAXED=50,
     >             MAXISD=300,NAMSBR='LIB   ')
*----
*  INPUT
*----
      INTEGER      ITYPLU,INTLIR
      CHARACTER    CARLIR*12
      REAL         REALIR
      DOUBLE PRECISION DBLLIR
*----
*  LOCAL PARAMETERS
*----
      CHARACTER    TEXT12*12,HSIGN*12,HVECT(MAXED)*8,HADD*8,HLIBX*12
      INTEGER      ISTATE(NSTATE),IPRINT,NBISOX,NBMIXX,MAXMIX,INDREC,
     >             NBISO,NGRO,NGT,NGF,NGFR,NL,ITRANC,ITIME,NLIB,NIDEPL,
     >             NCOMB,NEDMAC,NBMIX,NRES,MAXISM,IEN,ILCMLN,ILCMTY,
     >             IED,JED,KED,IDP,IBSTEP,MAXISO,NDEPL,NEDMA0,ITPROC,
     >             ISOADD,NADDXS,IPROB,IPROC,IMAC,NDEL,NFISS,IPRECI
      REAL         TMPDAY(3),DELT,TIMBRN 
      INTEGER      IKSTEP
      LOGICAL      LEXIST
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IADNAM
      REAL, ALLOCATABLE, DIMENSION(:) :: ENER,BSTD
*----
*  PARAMETER VALIDATION.
*----
      IF(NENTRY .EQ. 0) CALL XABORT(NAMSBR//': PARAMETER EXPECTED.')
      IF(IENTRY(1) .NE. 1 .AND.
     >   IENTRY(1) .NE. 2) CALL XABORT(NAMSBR//
     >': LCM OBJECT OR XSM FILE EXPECTED AT LHS.')
      IF(JENTRY(1) .NE. 0 .AND.
     >   JENTRY(1) .NE. 1) CALL XABORT(NAMSBR//': ENTRY'
     1 //' IN CREATE OR MODIFICATION MODE EXPECTED.')
      IPLIB=KENTRY(1)
*----
*  READ THE INPUT DATA.
*  DEFAULT OPTIONS:
*----
      IPRINT=1
      NBISOX=0
      NBMIXX=0
      IPLIBX=C_NULL_PTR
      IBSTEP=0
      HLIBX=' '
      LEXIST=(JENTRY(1) .EQ. 1)
      NDEPL=0
      IF(JENTRY(1) .EQ. 0) THEN
        MAXMIX=0
        INDREC=1
        NBISO=0
        NGRO=0
        NGT=0
        NGF=9999999
        NGFR=0
        NL=2
        ITRANC=0
        IPROB=0
        ITIME=1
        NLIB=0
        NIDEPL=0
        NCOMB=0
        NEDMAC=0
        NBMIX=0
        NRES=0
        IPROC=0
        IMAC=1
        NDEL=0
        NFISS=0
        ISOADD=0
        MAXISM=MAXISD
        IPRECI=4
*----
*  TRY TO FIND A READ-ONLY MICROLIB OR MACROLIB TO COPY IN THE LIBRARY
*----
        DO 130 IEN=2,NENTRY
          HSIGN=' '
          IF(IENTRY(IEN) .LE. 2 .AND. JENTRY(IEN) .EQ. 2) THEN
            CALL LCMLEN(KENTRY(IEN),'SIGNATURE',ILCMLN,ILCMTY)
            IF(ILCMLN .GT. 0) THEN
              CALL LCMGTC(KENTRY(IEN),'SIGNATURE',12,1,HSIGN)
            ELSE
              CALL LCMLEN(KENTRY(IEN),'MACROLIB',ILCMLN,ILCMTY)
              IF(ILCMLN .NE. 0) THEN
                CALL LCMSIX(KENTRY(IEN),'MACROLIB',1)
                HSIGN='L_MACROLIB'
              ENDIF
            ENDIF
            IF(HSIGN.EQ.'L_LIBRARY') THEN
              CALL LCMEQU(KENTRY(IEN),IPLIB)
              LEXIST=.TRUE.
              GO TO 140
            ELSE IF(HSIGN.EQ.'L_MACROLIB') THEN
              CALL LCMSIX(IPLIB,'MACROLIB',1)
              CALL LCMEQU(KENTRY(IEN),IPLIB)
              CALL LCMSIX(IPLIB,' ',2)
              INDREC=3
              CALL LCMGET(KENTRY(IEN),'STATE-VECTOR',ISTATE)
              NGRO=ISTATE(1)
              NGT=NGRO
              MAXMIX=ISTATE(2)
              NL=ISTATE(3)
              NADDXS=ISTATE(5)
              ITRANC=ISTATE(6)
              NDEL=ISTATE(7)
              IF(NGT .GT. 0) THEN
                ALLOCATE(ENER(2*NGT+1))
                CALL LCMGET(KENTRY(IEN),'ENERGY',ENER)
                CALL LCMGET(KENTRY(IEN),'DELTAU',ENER(NGT+2))
              ENDIF
              CALL LCMSIX(IPLIB,'MACROLIB',ILCMUP)
              CALL LCMEQU(KENTRY(IEN),IPLIB)
              IF(NADDXS .NE. 0) THEN
                IF(NADDXS .GT. MAXED-NEDMAC) CALL XABORT(NAMSBR//
     >          ': TOO MANY EXTRA EDITS REQUESTED')
                ALLOCATE(IADNAM(2*NADDXS))
                CALL LCMGET(IPLIB,'ADDXSNAME-P0',IADNAM)
                JED=0
                DO 120 IED=1,NADDXS
                  WRITE(HADD,'(2A4)') IADNAM(JED+1),IADNAM(JED+2)
                  DO 100 KED=1,NEDMAC
                    IF(HADD.EQ.HVECT(KED)) GO TO 110
 100              CONTINUE
                  NEDMAC=NEDMAC+1
                  HVECT(NEDMAC)=HADD
 110              CONTINUE
                  JED=JED+2
 120            CONTINUE
                DEALLOCATE(IADNAM)
              ENDIF
*----
*  WRITE ENERGY AND DELTAU ON MACROLIB
*----
              IF(NGT .GT. 0) THEN
                CALL LCMPUT(IPLIB,'ENERGY',NGT+1,2,ENER)
                CALL LCMPUT(IPLIB,'DELTAU',NGT,2,ENER(NGT+2))
              ENDIF
              CALL LCMSIX(IPLIB,'MACROLIB',ILCMDN)
              IF(NGT .GT. 0) THEN
                CALL LCMPUT(IPLIB,'ENERGY',NGT+1,2,ENER)
                CALL LCMPUT(IPLIB,'DELTAU',NGT,2,ENER(NGT+2))
                DEALLOCATE(ENER)
              ENDIF
              CALL LCMSIX(KENTRY(IEN),' ',0)
              GO TO 140
            ENDIF
          ENDIF
 130    CONTINUE
      ENDIF
*----
*  RECOVER STATE-VECTOR FROM EXISTING MICROLIB
*----
 140  IF(LEXIST) THEN
        CALL LCMGTC(IPLIB,'SIGNATURE',12,1,HSIGN)
        IF(HSIGN.NE.'L_LIBRARY') THEN
          TEXT12=HENTRY(1)
          CALL XABORT(NAMSBR//
     >    ': SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     >    '. L_LIBRARY EXPECTED.')
        ENDIF
        INDREC=2
        CALL LCMGET(IPLIB,'STATE-VECTOR',ISTATE)
        MAXMIX=ISTATE(1)
        NBISO=ISTATE(2)
        NGRO=ISTATE(3)
        NGT=NGRO
        NL=ISTATE(4)
        ITRANC=ISTATE(5)
        IPROB=ISTATE(6)
        ITIME=ISTATE(7)
        NLIB=ISTATE(8)
        NGF=ISTATE(9)
        NGFR=ISTATE(10)
        NIDEPL=ISTATE(11)
        NCOMB=ISTATE(12)
        NEDMAC=ISTATE(13)
        NBMIX=ISTATE(14)
        NRES=ISTATE(15)
        IPROC=ISTATE(17)
        IMAC=ISTATE(18)
        NDEL=ISTATE(19)
        NFISS=ISTATE(20)
        ISOADD=ISTATE(21)
        MAXISM=ISTATE(22)
        IPRECI=ISTATE(23)
        IF(NEDMAC.GT.0) THEN
          IF(NEDMAC .GT. MAXED) CALL XABORT(NAMSBR//': MAXED OVERFLOW')
          CALL LCMGTC(IPLIB,'ADDXSNAME-P0',8,NEDMAC,HVECT)
        ENDIF
      ENDIF
*----
*  READ LIBRARY DATA
*----
 155  CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
      IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//': KEYWORD EXPECTED')
      IF(CARLIR(1:4) .EQ. 'EDIT') THEN
*---
*  READ THE PRINT INDEX
*----
        CALL REDGET(ITYPLU,IPRINT,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.1) CALL XABORT(NAMSBR//
     >  ': VALUE FOR IPRINT EXPECTED')
      ELSE IF(CARLIR(1:4) .EQ. 'NGRO') THEN
*----
*  READ THE NUMBER OF ENERGY GROUPS.
*----
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >  ': VALUE FOR NGRO EXPECTED')
        IF(INDREC .EQ. 2) THEN
          IF(NGRO .NE. INTLIR) CALL XABORT(NAMSBR//
     >    ': INCOMPATIBLE VALUE OF NGRO')
        ELSE
          NGRO=INTLIR
        ENDIF
      ELSE IF(CARLIR(1:4) .EQ. 'MXIS') THEN
*----
*  CHANGE MAXIMUM NUMBER OF ISOTOPES
*----
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >  ': VALUE FOR MXIS EXPECTED')
        MAXISM=MAX(MAXISM,INTLIR)
      ELSE IF(CARLIR(1:4) .EQ. 'NMIX') THEN
*----
*  READ THE MAXIMUM NUMBER OF MATERIAL MIXTURES
*----
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >  ': VALUE FOR NMIX EXPECTED')
        MAXMIX=MAX(MAXMIX,INTLIR)
      ELSE IF(CARLIR(1:4) .EQ. 'CTRA') THEN
*----
*  READ TRANSPORT CORRECTION TYPE
*----
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >  ': CHARACTER CTRA TYPE EXPECTED')
        IF(CARLIR(1:4) .EQ. 'NONE') THEN
           ITRANC=0
        ELSE IF(CARLIR(1:4) .EQ. 'APOL') THEN
           ITRANC=1
        ELSE IF(CARLIR(1:4) .EQ. 'WIMS') THEN
           ITRANC=2
        ELSE IF(CARLIR(1:4) .EQ. 'OLDW') THEN
           ITRANC=3
        ELSE IF(CARLIR(1:4) .EQ. 'LEAK') THEN
           ITRANC=4
        ELSE
           CALL XABORT(NAMSBR//
     >     ': CTRA TYPE NONE, APOL, WIMS, OLDW OR LEAK EXPECTED.')
        ENDIF
      ELSE IF(CARLIR(1:4) .EQ. 'ANIS') THEN
*----
*  READ THE SCATTERING ANISOTROPY FOR TRANSPORT THEORY CASES
*----
        CALL REDGET(ITYPLU,NL,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >  ': VALUE FOR ANIS EXPECTED')
      ELSE IF(CARLIR(1:3) .EQ. 'ADJ') THEN
        IPROB=1
      ELSE IF(CARLIR(1:4) .EQ. 'PROM') THEN
        ITIME=2
      ELSE IF(CARLIR(1:7) .EQ. 'RDEPCHN') THEN
        ISOADD=1
      ELSE IF(CARLIR(1:7) .EQ. 'CDEPCHN') THEN
        ISOADD=0
      ELSE IF(CARLIR(1:4) .EQ. 'SKIP') THEN
        IPROC=-1
        IMAC=0
      ELSE IF(CARLIR(1:4) .EQ. 'INTR') THEN
        IPROC=0
        IMAC=0
      ELSE IF(CARLIR(1:4) .EQ. 'SUBG') THEN
        IPROC=1
        IMAC=0
      ELSE IF(CARLIR(1:4) .EQ. 'NEWL') THEN
        IPROC=2
        IMAC=0
      ELSE IF(CARLIR(1:4) .EQ. 'PTSL') THEN
        IPROC=4
        IMAC=0
      ELSE IF(CARLIR(1:4) .EQ. 'PTMC') THEN
        IPROC=5
        IMAC=0
      ELSE IF(CARLIR(1:2) .EQ. 'PT') THEN
        IPROC=3
        IMAC=0
      ELSE IF(CARLIR(1:4) .EQ. 'MACR') THEN
        IMAC=1
      ELSE IF(CARLIR(1:4) .EQ. 'CALE') THEN
        CALL REDGET(ITYPLU,IPRECI,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >  ': INTEGER VALUE EXPECTED FOR CALENDF ACCURACY')
      ELSE IF(CARLIR(1:4) .EQ. 'DEPL') THEN
        CALL LIBDEP(IPLIB,IPRINT,NIDEPL)
      ELSE IF(CARLIR.EQ.'ADED') THEN
        CALL REDGET(ITYPLU,NEDMA0,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .NE. 1) CALL XABORT(NAMSBR//
     >  ': VALUE FOR ADED EXPECTED')
        DO 170 IED=1,NEDMA0
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU .NE. 3) CALL XABORT(NAMSBR//
     >    ': TYPE FOR ADED EXPECTED')
          DO 160 JED=1,NEDMAC
            IF(CARLIR(:8) .EQ. HVECT(JED)) GO TO 170
 160      CONTINUE
          NEDMAC=NEDMAC+1
          IF(NEDMAC .GT. MAXED) CALL XABORT(NAMSBR//
     >    ': TOO MANY EXTRA EDITS REQUESTED')
          HVECT(NEDMAC)=CARLIR(:8)
 170    CONTINUE
      ELSE IF(CARLIR(1:4) .EQ. 'MIXS') THEN
         ITPROC=1
         GO TO 240
      ELSE IF(CARLIR(1:4) .EQ. 'MAXS') THEN
        ITPROC=2
        IF(INDREC .NE. 2)  CALL XABORT(NAMSBR//
     >  ': MAXS CAN ONLY BE USE TO UPDATE '//
     >  'A LIBRARY - IT CANNOT CREATE A NEW LIBRARY')
*----
*  TRY TO FIND A SECOND READ-ONLY LIBRARY TO MODIFY
*  ORIGINAL LIBRARY
*----
        DO 180 IEN=2,NENTRY
          IF(IENTRY(IEN) .LE. 2 .AND.
     >       JENTRY(IEN) .EQ.2 ) THEN
            CALL LCMGTC(KENTRY(IEN),'SIGNATURE',12,1,HSIGN)
            IF(HSIGN.EQ.'L_LIBRARY') THEN
              IPLIBX=KENTRY(IEN)
              HLIBX=HENTRY(IEN)
              CALL LCMGET(IPLIBX,'STATE-VECTOR',ISTATE)
              NBMIXX=ISTATE(1)
              NBISOX=ISTATE(2)
            ENDIF
          ENDIF
 180    CONTINUE
        IF(NBMIXX .EQ. 0) THEN
          NBMIXX=MAXMIX
          NBISOX=NBISO
          IPLIBX=IPLIB
        ENDIF
        TMPDAY(1)=0.0
        TMPDAY(2)=0.0
        TMPDAY(3)=0.0
        CALL LCMLEN(IPLIB,'MACROLIB',ILCMLN,ILCMTY)
        IF(ILCMLN .EQ. -1) THEN
          CALL LCMSIX(IPLIB,'MACROLIB',ILCMUP)
          CALL LCMLEN(IPLIB,'TIMESTAMP',ILCMLN,ILCMTY)
          IF(ILCMLN .GT. 0 .AND. ILCMLN .LE. 3) THEN
            CALL LCMGET(IPLIB,'TIMESTAMP',TMPDAY)
          ENDIF
          CALL LCMSIX(IPLIB,'MACROLIB',ILCMDN)
        ENDIF
        GO TO 240
      ELSE IF(CARLIR(1:4) .EQ. 'BURN') THEN
        ITPROC=2
        IF(INDREC .NE. 2) CALL XABORT(NAMSBR//
     >  ': MAXS CAN ONLY BE USE TO UPDATE '//
     >  'A LIBRARY - IT CANNOT CREATE A NEW LIBRARY')
        DO 190 IEN=2,NENTRY
          IF(IENTRY(IEN) .LE. 2 .AND.
     >       JENTRY(IEN) .EQ. 2) THEN
            CALL LCMGTC(KENTRY(IEN),'SIGNATURE',12,1,HSIGN)
            IF(HSIGN.EQ.'L_BURNUP') THEN
              IPLIBX=KENTRY(IEN)
              HLIBX=HENTRY(IEN)
              CALL LCMGET(IPLIBX,'STATE-VECTOR',ISTATE)
              NDEPL=ISTATE(3)
              NBISOX=ISTATE(4)
              NBMIXX=ISTATE(8)
              GO TO 200
            ENDIF
          ENDIF
 190    CONTINUE
        CALL XABORT(NAMSBR//': BURNUP FILE MISSING')
 200    CONTINUE
        ALLOCATE(BSTD(NDEPL))
        CALL LCMGET(IPLIBX,'DEPL-TIMES  ',BSTD)
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU .EQ. 3) CALL XABORT(NAMSBR//': INVALID BURNUP STEP')
        IF(ITYPLU .EQ. 2) THEN
          TMPDAY(1)=REALIR
          TIMBRN=0.000864*TMPDAY(1) 
          IF(TIMBRN .LE. 0.0) THEN  
            IBSTEP=1
            TMPDAY(1)=0.0
          ELSE
            IBSTEP=1
            IKSTEP=0
            DO 210 IDP=1,NDEPL
              DELT=ABS(TIMBRN-BSTD(IDP))
              IF(DELT .LT. 1.0E-6) THEN
                IBSTEP=IDP
                GO TO 220
              ELSE IF(TIMBRN .GT. BSTD(IDP)) THEN
                IKSTEP=IDP  
              ENDIF
 210        CONTINUE
            WRITE(IOUT,9000) TMPDAY
            WRITE(IOUT,9001) (BSTD(IDP)/0.000864,IDP=1,NDEPL)
            IBSTEP=MIN(IKSTEP+1,NDEPL)
            WRITE(IOUT,9002) BSTD(IBSTEP)/0.000864
 220        CONTINUE
          ENDIF
        ELSE IF(ITYPLU .EQ. 1) THEN
          IBSTEP=INTLIR
          IF(IBSTEP .LE. 0 ) THEN 
            WRITE(IOUT,9010) 
            IBSTEP=1
            WRITE(IOUT,9010) BSTD(IBSTEP)/0.000864
          ELSE IF(IBSTEP .GT. NDEPL) THEN 
            IBSTEP=NDEPL 
            WRITE(IOUT,9011) BSTD(IBSTEP)/0.000864
          ENDIF
          TMPDAY(1)=BSTD(IBSTEP)/0.000864
        ENDIF
        DEALLOCATE(BSTD) 
        TMPDAY(2)=0.0
        TMPDAY(3)=0.0
        IF(IPRINT .GE. 1) WRITE(IOUT,6000) IBSTEP,TMPDAY(1)
        GO TO 240
      ELSE IF(CARLIR(1:1).EQ.';') THEN
*       SAVE THE LIBRARY SPECIFIC INFORMATION.
        TEXT12='L_LIBRARY'
        CALL LCMPTC(IPLIB,'SIGNATURE',12,1,TEXT12)
        CALL XDISET(ISTATE,NSTATE,0)
        ISTATE(1)=MAXMIX
        ISTATE(2)=NBISO
        ISTATE(3)=NGRO
        ISTATE(4)=NL
        ISTATE(5)=ITRANC
        ISTATE(6)=IPROB
        ISTATE(7)=ITIME
        ISTATE(8)=NLIB
        ISTATE(9)=NGF
        ISTATE(10)=NGFR
        ISTATE(11)=NIDEPL
        ISTATE(12)=NCOMB
        ISTATE(13)=NEDMAC
        ISTATE(14)=NBMIX
        ISTATE(15)=NRES
        ISTATE(17)=IPROC
        ISTATE(18)=IMAC
        ISTATE(19)=NDEL
        ISTATE(20)=NFISS
        ISTATE(21)=ISOADD
        ISTATE(22)=MAXISM
        ISTATE(23)=IPRECI
        CALL LCMPUT(IPLIB,'STATE-VECTOR',NSTATE,1,ISTATE)
        GO TO 250
      ELSE
        CALL XABORT(NAMSBR//': '//CARLIR//' IS AN INVALID KEY-WORD.')
      ENDIF
      GO TO 155
*----
*  PROCESS THE LIB: MODULE INPUT DATA.
*----
 240  CONTINUE
      IF(MAXMIX.EQ.0) CALL XABORT(NAMSBR//': MAXMIX NOT YET DEFINED.')
      MAXISO=MAX(NIDEPL,MAXISM)*MAXMIX
      IF(ITPROC .EQ. 1) THEN
        CALL LIBINP(MAXMIX,MAXED ,MAXISO,IPLIB ,INDREC,IPRINT,
     >              NBISO ,NGRO  ,NGT   ,NL    ,ITRANC,IPROB ,
     >              ITIME ,NLIB  ,NGF   ,NGFR  ,NIDEPL,NCOMB ,
     >              NEDMAC,NBMIX ,NRES  ,IPROC ,IMAC  ,NDEL  ,
     >              ISOADD,MAXISM,HVECT ,IPRECI)
      ELSE IF(ITPROC .GE. 2) THEN
*----
*  ALLOCATE
*----
        IF(NGRO .EQ. 0) CALL XABORT(NAMSBR//
     >  ': NUMBER OF GROUP REQUIRED FOR MAXS OF BURN')
        CALL LIBMAC(IPLIB ,IPLIBX,IPRINT,MAXISO,NBISO ,NBISOX,
     >              IBSTEP,NBMIX ,NBMIXX,NGRO  ,TMPDAY)
      ENDIF
  250 IF(IPRINT .GE. 5) CALL LCMLIB(IPLIB)
      RETURN
*----
*  FORMATS
*----
 6000 FORMAT(' LIBRARY UPDATE AT BURNUP STEP : ',I5,
     >       ' BURNUP TIME = ',F20.7,' DAYS') 
 9000 FORMAT(' **** WARNING  *****'/
     >       ' INVALID BURNUP TIME =',F20.7,' DAYS'/
     >       ' BURNUP TABULATION (DAYS) ')
 9001 FORMAT(6F20.7)
 9002 FORMAT(' BURNUP STEP SELECTED =',F20.7,' DAYS')
 9010 FORMAT(' **** WARNING  *****'/
     >       ' BURNUP STEP NEGATIVE '/
     >       ' USE FIRST BURNUP STEP AT ',F20.7,' DAYS')
 9011 FORMAT(' **** WARNING  *****'/
     >       ' BURNUP STEP TOO LARGE '/
     >       ' USE LAST BURNUP STEP AT ',F20.7,' DAYS')
      END
