*DECK MAC
      SUBROUTINE MAC(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Macroscopic cross sections and diffusion coefficients input module.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert and G. Marleau
*
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1) create or modification type(L_MACROLIB);
*         HENTRY(2) optional read-only type(L_MACROLIB) or
*                     type(L_OPTIMIZE).
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
      TYPE(C_PTR) IPMACR,JPLIST,KPLIST
      PARAMETER(NSTATE=40,IOUT=6)
      CHARACTER TEXT12*12,HSIGN*12,CARLIR*12
      INTEGER ISTATE(NSTATE)
      INTEGER NALBP
      DOUBLE PRECISION DBLLIR
*----
*  PARAMETER VALIDATION
*----
      IF(NENTRY.EQ.0) CALL XABORT('MAC: PARAMETER EXPECTED.')
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2)) CALL XABORT('MAC: LCM '
     1 //'OBJECT EXPECTED AT LHS.')
      IF((JENTRY(1).NE.0).AND.(JENTRY(1).NE.1)) CALL XABORT('MAC: ENTR'
     1 //'Y IN CREATE OR MODIFICATION MODE EXPECTED.')
      ITYPE=JENTRY(1)
      IPMACR=KENTRY(1)
      NGO=0
*----
*  LOOK FOR OTHER MACROLIB IN SET OF DATA STRUCTURES
*----
      NMACSR=1
      NOLDMX=0
      NGO=0
      NLO=0
      NFO=0
      NEO=0
      ITO=0
      IPMAC2=0
      IF(NENTRY.GT.2) CALL XABORT('MAC: ONLY TWO OBJECTS PERMITTED.')
      IF(NENTRY.EQ.2) THEN
        IPMAC2=NENTRY
        IF((IENTRY(IPMAC2).NE.1).AND.(IENTRY(IPMAC2).NE.2)) THEN
           CALL XABORT('MAC: INVALID STRUCTURE TYPE FOR SECOND OBJECT.')
        ELSE IF(JENTRY(IPMAC2).NE.2) THEN
           CALL XABORT('MAC: DATA STRUCTURE '//HENTRY(IPMAC2)//' NOT '
     1     //'IN READ-ONLY MODE')
        ENDIF
        CALL LCMGTC(KENTRY(IPMAC2),'SIGNATURE',12,1,HSIGN)
        IF(HSIGN.EQ.'L_MACROLIB') THEN
           NMACSR=2
           CALL XDISET(ISTATE,NSTATE,0)
           CALL LCMGET(KENTRY(IPMAC2),'STATE-VECTOR',ISTATE)
           IF(ISTATE(2).GT.0) THEN
              NOLDMX=ISTATE(2)
              NGO=ISTATE(1)
              NLO=ISTATE(3)
              NFO=ISTATE(4)
              NEO=ISTATE(5)
              ITO=ISTATE(6)
           ELSE
              NMACSR=1
           ENDIF
        ELSE IF(HSIGN.EQ.'L_LIBRARY') THEN
           NMACSR=-2
           CALL LCMSIX(KENTRY(IPMAC2),'MACROLIB',1)
           CALL XDISET(ISTATE,NSTATE,0)
           CALL LCMGET(KENTRY(IPMAC2),'STATE-VECTOR',ISTATE)
           IF(ISTATE(2).GT.0) THEN
              NOLDMX=ISTATE(2)
              NGO=ISTATE(1)
              NLO=ISTATE(3)
              NFO=ISTATE(4)
              NEO=ISTATE(5)
              ITO=ISTATE(6)
           ELSE
              CALL LCMSIX(KENTRY(IPMAC2),' ',2)
              NMACSR=1
           ENDIF
        ELSE IF(HSIGN.EQ.'L_OPTIMIZE') THEN
           CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
           IF(ITYPLU.NE.10) CALL XABORT('MAC: NO INPUT DATA EXPECTED.')
           CALL MACOPT(IPMACR,KENTRY(IPMAC2))
           RETURN
        ELSE
           CALL XABORT('MAC: SECOND DATA STRUCTURE HAS INVALID SIGNATU'
     1     //'RE SET TO '//HSIGN//'.')
        ENDIF
      ENDIF
*----
*  READ THE INPUT DATA
*----
*     DEFAULT OPTIONS:
      IPRINT=1
      IF(ITYPE.EQ.0) THEN
         INDREC=1
         NANISO=1
         NGROUP=0
         NBMIX=0
         NIFISS=0
         NEDMAC=0
         ITRANC=0
         NDELG=0
         NALBP=0
         NSTEP=0
         IDF=0
         NPART0=0
      ELSE
         INDREC=2
         CALL LCMGTC(IPMACR,'SIGNATURE',12,1,HSIGN)
         IF(HSIGN.EQ.'L_LIBRARY') THEN
            CALL LCMSIX(IPMACR,'MACROLIB',1)
         ELSE IF(HSIGN.NE.'L_MACROLIB') THEN
            TEXT12=HENTRY(1)
            CALL XABORT('MAC: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1      '. L_MACROLIB OR L_LIBRARY EXPECTED.')
         ENDIF
         CALL XDISET(ISTATE,NSTATE,0)
         CALL LCMGET(IPMACR,'STATE-VECTOR',ISTATE)
         NGROUP=ISTATE(1)
         NBMIX=ISTATE(2)
         NANISO=ISTATE(3)
         NIFISS=ISTATE(4)
         NEDMAC=ISTATE(5)
         ITRANC=ISTATE(6)
         NDELG=ISTATE(7)
         NALBP =ISTATE(8)
         NSTEP=ISTATE(11)
         IDF=ISTATE(12)
         NPART0=ISTATE(17)
      ENDIF
*----
*  PROCESS THE MAC: INPUT DATA
*----
      IF(NMACSR.EQ.1) THEN
         CALL MACDRV(IPMACR,INDREC,IPRINT,IDF,NBMIX,NGROUP,NANISO,
     1   NIFISS,NEDMAC,ITRANC,NDELG,NSTEP,NALBP)
      ELSE
         NNEWMX=0
         NANISO=MAX(NLO,NANISO)
         NIFISS=NFO+NIFISS
         NEDMAC=MAX(NEDMAC,NEO)
         ITRANC=MAX(ITRANC,ITO)
*----
*  TAKE MACROSCOPIC XS FROM OLD MACROLIB
*  READ MAIN INPUT PARAMETERS UNTIL KEYWORD MIX FOUND
*----
 1000    CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
         IF(ITYPLU.NE.3) CALL XABORT('MAC: CHARACTER KEYWORD EXPECTED.')
         IF(CARLIR.EQ.'EDIT') THEN
            CALL REDGET(ITYPLU,IPRINT,REALIR,CARLIR,DBLLIR)
            IF(ITYPLU.NE.1) CALL XABORT('MAC: EDIT LEVEL SHOULD BE AN '
     1      //'INTEGER.')
         ELSE IF(CARLIR.EQ.'NMIX') THEN
            CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
            IF(ITYPLU.NE.1) CALL XABORT('MAC: READ ERROR - NUMBER OF M'
     1      //'IXTURES EXPECTED.')
            NNEWMX=MAX(INTLIR,NNEWMX)
         ELSE IF(CARLIR.EQ.'ANIS') THEN
            CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
            IF(ITYPLU.NE.1) CALL XABORT('MAC: READ ERROR - ANIS LEVEL '
     1      //'EXPECTED.')
            NANISO=MAX(NANISO,INTLIR)
         ELSE IF(CARLIR.EQ.'NIFI') THEN
            CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
            IF(ITYPLU.NE.1) CALL XABORT('MAC: READ ERROR - NUMBER FISS'
     1      //'ILE ISOTOPES EXPECTED.')
            NIFISS=MAX(INTLIR,NIFISS)
         ELSE IF(CARLIR(1:4).EQ.'CTRA') THEN
            CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
            IF(ITYPLU.NE.3) CALL XABORT('MAC: CTRA MUST BE FOLLOWED BY'
     1      //' CHARACTER.')
            IF(CARLIR.EQ.'OFF') THEN
               ITRANC=0
            ELSE IF(CARLIR.EQ.'ON') THEN
               ITRANC=2
            ELSE
               CALL XABORT('MAC: INVALID CTRA OPTION.')
            ENDIF
         ELSE IF(CARLIR(1:3).EQ.'MIX') THEN
            GO TO 1005
         ELSE
            CALL XABORT('MAC: KEYWORD '//CARLIR//' NOT PERMITTED.')
         ENDIF
         GO TO 1000
 1005    CONTINUE
         NTOTMX=NOLDMX+NBMIX+NNEWMX
         IF(NGROUP.EQ.0) THEN
            IF(NGO.EQ.0) CALL XABORT('MAC: MACROLIBS HAVE 0 GROUP.')
            NGROUP=NGO
         ELSE IF(NGROUP.NE.NGO) THEN
            CALL XABORT('MAC: MACROLIBS HAVE DIFFERENT GROUP STRUCTURE'
     1      //'S.')
         ENDIF
         CALL MACUPD(NENTRY,KENTRY,IPRINT,NTOTMX,NBMIX,NGROUP,
     1               NANISO,NIFISS,NEDMAC,ITRANC)
         IF(NMACSR.EQ.-2) CALL LCMSIX(KENTRY(IPMAC2),' ',2)
      ENDIF
*
      IF(ITYPE.EQ.0) THEN
         HSIGN='L_MACROLIB'
         CALL LCMPTC(IPMACR,'SIGNATURE',12,1,HSIGN)
      ENDIF
      IF(ITYPE.NE.2) THEN
         CALL XDISET(ISTATE,NSTATE,0)
         ISTATE(1)=NGROUP
         ISTATE(2)=NBMIX
         ISTATE(3)=NANISO
         ISTATE(4)=NIFISS
         ISTATE(5)=NEDMAC
         ISTATE(6)=ITRANC
         ISTATE(7)=NDELG
         ISTATE(8)=NALBP
         ISTATE(11)=NSTEP
         ISTATE(12)=IDF
         ISTATE(17)=NPART0
         IF(ITRANC.NE.0) ISTATE(6)=2
         JPLIST=LCMGID(IPMACR,'GROUP')
         KPLIST=LCMGIL(JPLIST,1)
         CALL LCMLEN(KPLIST,'DIFF',ILONG,ITYLCM)
         IF(ILONG.GT.0) ISTATE(9)=1
         CALL LCMLEN(KPLIST,'DIFFX',ILONG,ITYLCM)
         IF(ILONG.GT.0) ISTATE(9)=2
         CALL LCMPUT(IPMACR,'STATE-VECTOR',NSTATE,1,ISTATE)
      ENDIF
      IF(IPRINT.GT.1) CALL LCMLIB(IPMACR)
      IF(IPRINT.GT.0) WRITE(IOUT,100) IPRINT,(ISTATE(I),I=1,9),
     1 ISTATE(11),ISTATE(12),ISTATE(17)
      CALL LCMSIX(IPMACR,' ',0)
      RETURN
*
  100 FORMAT(/8H OPTIONS/8H -------/
     1 7H IPRINT,I6,30H   (0=NO PRINT/1=SHORT/2=MORE)/
     2 7H NGROUP,I6,28H   (NUMBER OF ENERGY GROUPS)/
     3 7H NBMIX ,I6,39H   (NUMBER OF MIXTURES IN THE MACROLIB)/
     4 7H NANISO,I6,34H   (MAXIMUM SCATTERING ANISOTROPY)/
     5 7H NIFISS,I6,45H   (MAXIMUM NUMBER OF FISSILE ISOTOPES IN A M,
     6 7HIXTURE)/
     7 7H NEDMAC,I6,34H   (NUMBER OF CROSS SECTION EDITS)/
     8 7H ITRANC,I6,45H   (0=NO TRANSPORT CORRECTION/1=APOLLO TYPE/2,
     9 43H=RECOVER FROM LIBRARY/4=LEAKAGE CORRECTION)/
     1 7H NLG   ,I6,39H   (NUMBER OF DELAYED PRECURSOR GROUPS)/
     2 7H NALB  ,I6,31H   (NUMBER OF PHYSICAL ALBEDOS)/
     3 7H ILEAK ,I6,40H   (1=DIFF AVAILABLE; 2=DIFFX AVAILABLE)/
     4 7H NSTEP ,I6,39H   (NUMBER OF PERTURBATION DIRECTORIES)/
     5 7H IDF   ,I6,48H   (=0/2 BOUNDARY FLUXES FOR ADF ABSENT/PRESENT)/
     6 7H NPART0,I6,34H   (NUMBER OF COMPANION PARTICLES))
      END
