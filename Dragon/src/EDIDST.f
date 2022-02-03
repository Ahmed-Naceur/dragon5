*DECK EDIDST
      SUBROUTINE EDIDST(IPEDIT,IPRINT,NL,NGCOND,NMERGE,NSTATS,ILEAKS,
     >                  EIGENK,B2,VOLMER,WLETYC,RATECM,FLUXCM,SCATTS,
     >                  OLDNAM,NW,NTAUXT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Reaction rates and fluxes statistics.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau
*
*Parameters: input
* IPEDIT  pointer to the edition LCM object.
* IPRINT  print level;
*         = 0 no print;
*         = 1 print fluxes;
*         = 2 1+print reaction rates;
*         = 3 2+print homogenized cross sections.
* NL      number of legendre orders.
* NGCOND  number of groups.
* NMERGE  number of regions.
* NSTATS  statistics options:
*         = 1  flux stats;
*         = 2  reaction rates stats;
*         = 3  flux+reaction rates stats;
*         =-1  delta sigma calculations.
* ILEAKS  type of leakage calculation:
*         = 0 no leakage;
*         = 1 homogeneous leakage (Diffon);
*         = 2 isotropic streaming (Ecco);
*         = 3 anisotropic streaming (Tibere).
* EIGENK  New eigenvalue.
* B2      New buckling.
* VOLMER  volume of merged regions.
* WLETYC  lethargy width.
* RATECM  averaged region/group cross sections:
*         = RATECM(*,1) = total P0;
*         = RATECM(*,2) = total P1;
*         = RATECM(*,NW+2) = absorption;
*         = RATECM(*,NW+3) = fission;
*         = RATECM(*,NW+4) = fixed sources / productions;
*         = RATECM(*,NW+5) = leakage;
*         = RATECM(*,NW+6) = total out of group scattering;
*         = RATECM(*,NW+7) = diagonal scattering x-s;
*         = RATECM(*,NW+8) = chi;
*         = RATECM(*,NW+9) = wims type transport correction;
*         = RATECM(*,NW+10) = x-directed leakage;
*         = RATECM(*,NW+11) = y-directed leakage;
*         = RATECM(*,NW+12) = z-directed leakage;
*         = RATECM(*,NW+13) = nu-sigf for delayed neutrons;
*         = RATECM(*,NW+13+NDEL) = fission spectra for delayed neutrons.
* FLUXCM  integrated region/group fluxes:
*         = FLUXCM(*,1) = fluxes P0;
*         = FLUXCM(*,2) = fluxes P1.
* SCATTS  new scattering matrix.
* OLDNAM  name of reference calculation directory.
* NW      type of weighting for PN cross section info (=0 P0; =1 P1).
* NTAUXT  number of reaction rate edits.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPEDIT
      INTEGER     IPRINT,NL,NGCOND,NMERGE,NSTATS,ILEAKS,NW,NTAUXT
      REAL        EIGENK,B2,VOLMER(NMERGE),WLETYC(NGCOND),
     >            RATECM(NMERGE,NGCOND,NTAUXT),
     >            FLUXCM(NMERGE,NGCOND,NW+1),
     >            SCATTS(NMERGE,NGCOND,NGCOND,NL)
      CHARACTER   OLDNAM*12
*----
*  LOCAL VARIABLES
*----
      PARAMETER  (IUNOUT=6,ILCMUP=1,ILCMDN=2,TRONCE=0.001,NSTATE=40)
      TYPE(C_PTR) JPEDIT,KPEDIT
      CHARACTER   CXSNAM*12,CM*2
      INTEGER     IDATA(NSTATE)
      DOUBLE PRECISION SOUOLD,SOUNEW
      REAL        EIGOLD,B2OLD(4)
      INTEGER     IEGC,IB2C
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INGSCT,IFGSCT,IPOSCT
      REAL, ALLOCATABLE, DIMENSION(:) :: FLXNEW,FLXOLD,OLDRAT,XSCAT,
     > DELSC
*----
*  SCRATCH STORAGE ALLOCATION
*   FLXNEW  new fluxes.
*   FLXOLD  old fluxes.
*   OLDRAT  old rates.
*----
      ALLOCATE(FLXNEW(NMERGE),FLXOLD(NMERGE),OLDRAT(NMERGE+NGCOND))
*----
      CALL LCMLEN(IPEDIT,OLDNAM,ILCMLN,ILCMTY)
      IF(ILCMLN.GE.0) THEN
        WRITE(IUNOUT,7000) OLDNAM
        CALL LCMLIB(IPEDIT)
        RETURN
      ENDIF
      CALL LCMSIX(IPEDIT,OLDNAM,ILCMUP)
      CALL LCMLEN(IPEDIT,'DELTAU',ILCMLN,ILCMTY)
      IF(ILCMLN.GT.0) THEN
        CALL LCMGET(IPEDIT,'DELTAU',OLDRAT)
        ERRREL=0.0
        DO 80 IGR=1,NGCOND
          ERRREL=ERRREL+ABS(WLETYC(IGR)-OLDRAT(IGR))
 80     CONTINUE
        ERRREL=ERRREL/NGCOND
        IF(ERRREL.GT.TRONCE) THEN
          WRITE(IUNOUT,7004) ERRREL,TRONCE
          CALL LCMSIX(IPEDIT,' ',ILCMDN)
          RETURN
        ENDIF
      ENDIF
      CALL LCMSIX(IPEDIT,'MACROLIB',ILCMUP)
      CALL LCMGET(IPEDIT,'STATE-VECTOR',IDATA)
      NGOLD=IDATA(1)
      NROLD=IDATA(2)
      CALL LCMLEN(IPEDIT,'K-EFFECTIVE',ILCMLN,ILCMTY)
      IEGC=0
      IF(ILCMLN .EQ. 1) THEN
        CALL LCMGET(IPEDIT,'K-EFFECTIVE',EIGOLD)
        IEGC=1
      ELSE
        EIGOLD=1.0
      ENDIF
      IB2C=0
      CALL XDRSET(B2OLD,4,0.0)
      CALL LCMLEN(IPEDIT,'B2  HETE',ILCMLN,ILCMTY)
      IF(ILCMLN .EQ. 1) THEN
        CALL LCMGET(IPEDIT,'B2  HETE',B2OLD)
        IB2C=2
      ELSE
        CALL LCMLEN(IPEDIT,'B2  B1HOM',ILCMLN,ILCMTY)
        IF(ILCMLN .EQ. 1) THEN
          CALL LCMGET(IPEDIT,'B2  B1HOM',B2(4))
          IB2C=1
        ENDIF
      ENDIF
      IF( (NROLD.NE.NMERGE).OR.(NGOLD.NE.NGCOND)) THEN
        WRITE(IUNOUT,7001) NMERGE,NROLD,NGCOND,NGOLD
        CALL LCMSIX(IPEDIT,' ',ILCMDN)
        CALL LCMSIX(IPEDIT,' ',ILCMDN)
        RETURN
      ENDIF
*----
*  COMPUTE TOTAL SOURCES FOR RELATIVE FLUX NORMALIZATION
*----
      SOUOLD=0.0D0
      SOUNEW=0.0D0
      JPEDIT=LCMGID(IPEDIT,'GROUP')
      DO 300 IGR=1,NGCOND
        KPEDIT=LCMGIL(JPEDIT,IGR)
        CALL LCMGET(KPEDIT,'FLUX-INTG',FLXOLD)
        DO 310 IREG=1,NMERGE
          FLXNEW(IREG)=FLUXCM(IREG,IGR,1)
 310    CONTINUE
        CALL LCMGET(KPEDIT,'PRODUCTION',OLDRAT)
        DO 320 IREG=1,NMERGE
          SOUOLD=SOUOLD+DBLE(FLXOLD(IREG))*DBLE(OLDRAT(IREG))
          SOUNEW=SOUNEW+DBLE(FLXNEW(IREG))*DBLE(RATECM(IREG,IGR,NW+4))
 320    CONTINUE
 300  CONTINUE
*----
*  CHECK FOR VOLUME CONSISTENCE
*----
      VOLTOT=0.0
      VOLT2=0.0
      CALL LCMGET(IPEDIT,'VOLUME',OLDRAT)
      DO 100 IREG=1,NMERGE
        VOLTOT=VOLTOT+VOLMER(IREG)
        VOLT2=VOLT2+OLDRAT(IREG)
 100  CONTINUE
      VOLREL=VOLT2 /VOLTOT
      DO 101 IREG=1,NMERGE
        VREL1=VOLMER(IREG)/VOLTOT
        VREL2=OLDRAT(IREG)/VOLT2
        ERRREL=ABS(VREL1-VREL2)/VREL2
        IF(ERRREL.GT.TRONCE) THEN
          WRITE(IUNOUT,7002) VOLREL,IREG,ERRREL,TRONCE
          CALL LCMSIX(IPEDIT,' ',ILCMDN)
          CALL LCMSIX(IPEDIT,' ',ILCMDN)
          RETURN
        ENDIF
 101  CONTINUE
      IF((SOUOLD.EQ.0.0D0).OR.(SOUNEW.EQ.0.0D0)) THEN
        WRITE(IUNOUT,7005)
        SOUREL=1.0
      ELSE
        SOUREL=REAL(SOUOLD/SOUNEW)
      ENDIF
      WRITE(IUNOUT,6000) OLDNAM,NGCOND,NMERGE,VOLREL,SOUREL
      IF(IEGC .EQ. 1) THEN
        WRITE(IUNOUT,6010)  EIGOLD,EIGENK,1000.*(EIGENK-EIGOLD)
      ENDIF
      IF(IB2C .GE. 1) THEN
        WRITE(IUNOUT,6011)  B2OLD(4),B2(4),B2(4)-B2OLD(4)
      ENDIF
      IF(NSTATS.EQ.-1) THEN
        ALLOCATE(INGSCT(NMERGE),IFGSCT(NMERGE),IPOSCT(NMERGE))
        ALLOCATE(XSCAT(NMERGE*NGCOND),DELSC(NGCOND))
      ENDIF
      DO 210 IGR=1,NGCOND
        KPEDIT=LCMGIL(JPEDIT,IGR)
        WRITE(IUNOUT,6001) IGR
        CXSNAM='FLUX-INTG'
        CALL LCMGET(KPEDIT,CXSNAM,FLXOLD)
        IF(NSTATS.GE.1) THEN
          IF(NSTATS.NE.2) THEN
            WRITE(IUNOUT,6002) CXSNAM
            ITYPE=1
            CALL EDISTA(IPRINT,NMERGE,ITYPE,VOLMER,SOUREL,VOLT2,
     >                  FLUXCM(1,IGR,1),FLXOLD,RATECM(1,IGR,1),OLDRAT)
          ENDIF
          IF(NSTATS.GE.2) THEN
            DO 102 IREG=1,NMERGE
              FLXNEW(IREG)=FLUXCM(IREG,IGR,1)
 102        CONTINUE
            ITYPE=2
            CXSNAM='NTOT0'
            CALL LCMGET(KPEDIT,CXSNAM,OLDRAT)
            WRITE(IUNOUT,6002) CXSNAM
            CALL EDISTA(IPRINT,NMERGE,ITYPE,VOLMER,SOUREL,VOLT2,
     >                  FLXNEW,FLXOLD,RATECM(1,IGR,1),OLDRAT)
            CXSNAM='ABS'
            CALL LCMGET(KPEDIT,CXSNAM,OLDRAT)
            WRITE(IUNOUT,6002) CXSNAM
            CALL EDISTA(IPRINT,NMERGE,ITYPE,VOLMER,SOUREL,VOLT2,
     >                  FLXNEW,FLXOLD,RATECM(1,IGR,NW+2),OLDRAT)
            CXSNAM='PRODUCTION'
            CALL LCMGET(KPEDIT,CXSNAM,OLDRAT)
            WRITE(IUNOUT,6002) CXSNAM
            CALL EDISTA(IPRINT,NMERGE,ITYPE,VOLMER,SOUREL,VOLT2,
     >                  FLXNEW,FLXOLD,RATECM(1,IGR,NW+4),OLDRAT)
            IF(IDATA(4).EQ.1) THEN
              CXSNAM='NUSIGF'
              CALL LCMGET(KPEDIT,CXSNAM,OLDRAT)
              WRITE(IUNOUT,6002) CXSNAM
              CALL EDISTA(IPRINT,NMERGE,ITYPE,VOLMER,SOUREL,VOLT2,
     >                  FLXNEW,FLXOLD,RATECM(1,IGR,NW+3),OLDRAT)
            ENDIF
            DO 134 IL=1,NL
              WRITE (CM,'(I2.2)') IL-1
              CXSNAM='SIGW'//CM
              CALL LCMLEN(KPEDIT,CXSNAM,ILCMLN,ILCMTY)
              IF(ILCMLN.EQ.NMERGE) THEN
                CALL LCMGET(KPEDIT,CXSNAM,OLDRAT)
                WRITE(IUNOUT,6002) CXSNAM
                CALL EDISTA(IPRINT,NMERGE,ITYPE,VOLMER,SOUREL,VOLT2,
     >                      FLXNEW,FLXOLD,SCATTS(1,IGR,IGR,IL),OLDRAT)
              ENDIF
 134        CONTINUE
            IF(ILEAKS.EQ.3) THEN
              CXSNAM='DIFFX'
              CALL LCMLEN(KPEDIT,CXSNAM,ILCMLN,ILCMTY)
              IF(ILCMLN.GT.0) THEN
                CALL LCMGET(KPEDIT,CXSNAM,OLDRAT)
                WRITE(IUNOUT,6003) CXSNAM
                CALL EDISTA(IPRINT,NMERGE,ITYPE,VOLMER,SOUREL,VOLT2,
     >                      FLXNEW,FLXOLD,RATECM(1,IGR,NW+10),OLDRAT)
              ENDIF
              CXSNAM='DIFFY'
              CALL LCMLEN(KPEDIT,CXSNAM,ILCMLN,ILCMTY)
              IF(ILCMLN.GT.0) THEN
                CALL LCMGET(KPEDIT,CXSNAM,OLDRAT)
                WRITE(IUNOUT,6003) CXSNAM
                CALL EDISTA(IPRINT,NMERGE,ITYPE,VOLMER,SOUREL,VOLT2,
     >                      FLXNEW,FLXOLD,RATECM(1,IGR,NW+11),OLDRAT)
              ENDIF
              CXSNAM='DIFFZ'
              CALL LCMLEN(KPEDIT,CXSNAM,ILCMLN,ILCMTY)
              IF(ILCMLN.GT.0) THEN
                CALL LCMGET(KPEDIT,CXSNAM,OLDRAT)
                WRITE(IUNOUT,6003) CXSNAM
                CALL EDISTA(IPRINT,NMERGE,ITYPE,VOLMER,SOUREL,VOLT2,
     >                      FLXNEW,FLXOLD,RATECM(1,IGR,NW+12),OLDRAT)
              ENDIF
            ENDIF
            CXSNAM='DIFF'
            CALL LCMLEN(KPEDIT,CXSNAM,ILCMLN,ILCMTY)
            IF(ILCMLN.GT.0) THEN
              CALL LCMGET(KPEDIT,CXSNAM,OLDRAT)
              WRITE(IUNOUT,6003) CXSNAM
              CALL EDISTA(IPRINT,NMERGE,ITYPE,VOLMER,SOUREL,VOLT2,
     >                    FLXNEW,FLXOLD,RATECM(1,IGR,NW+5),OLDRAT)
            ENDIF
          ENDIF
        ELSE IF(NSTATS.EQ.-1) THEN
          ITYPE=3
          CXSNAM='NTOT0'
          CALL LCMGET(KPEDIT,CXSNAM,OLDRAT)
          WRITE(IUNOUT,6003) CXSNAM
          CALL EDISTA(IPRINT,NMERGE,ITYPE,VOLMER,SOUREL,VOLT2,
     >                FLXNEW,FLXOLD,RATECM(1,IGR,1),OLDRAT)
          CXSNAM='ABS'
          CALL LCMGET(KPEDIT,CXSNAM,OLDRAT)
          WRITE(IUNOUT,6003) CXSNAM
          CALL EDISTA(IPRINT,NMERGE,ITYPE,VOLMER,SOUREL,VOLT2,
     >                FLXNEW,FLXOLD,RATECM(1,IGR,NW+2),OLDRAT)
          IF(IDATA(4).EQ.1) THEN
            CXSNAM='NUSIGF'
            CALL LCMGET(KPEDIT,CXSNAM,OLDRAT)
            WRITE(IUNOUT,6003) CXSNAM
            CALL EDISTA(IPRINT,NMERGE,ITYPE,VOLMER,SOUREL,VOLT2,
     >                FLXNEW,FLXOLD,RATECM(1,IGR,NW+3),OLDRAT)
          ENDIF
          IF(ILEAKS.EQ.3) THEN
            CXSNAM='DIFFX'
            CALL LCMLEN(KPEDIT,CXSNAM,ILCMLN,ILCMTY)
            IF(ILCMLN.GT.0) THEN
              CALL LCMGET(KPEDIT,CXSNAM,OLDRAT)
              WRITE(IUNOUT,6003) CXSNAM
              CALL EDISTA(IPRINT,NMERGE,ITYPE,VOLMER,SOUREL,VOLT2,
     >                    FLXNEW,FLXOLD,RATECM(1,IGR,NW+10),OLDRAT)
            ENDIF
            CXSNAM='DIFFY'
            CALL LCMLEN(KPEDIT,CXSNAM,ILCMLN,ILCMTY)
            IF(ILCMLN.GT.0) THEN
              CALL LCMGET(KPEDIT,CXSNAM,OLDRAT)
              WRITE(IUNOUT,6003) CXSNAM
              CALL EDISTA(IPRINT,NMERGE,ITYPE,VOLMER,SOUREL,VOLT2,
     >                    FLXNEW,FLXOLD,RATECM(1,IGR,14),OLDRAT)
            ENDIF
            CXSNAM='DIFFZ'
            CALL LCMLEN(KPEDIT,CXSNAM,ILCMLN,ILCMTY)
            IF(ILCMLN.GT.0) THEN
              CALL LCMGET(KPEDIT,CXSNAM,OLDRAT)
              WRITE(IUNOUT,6003) CXSNAM
              CALL EDISTA(IPRINT,NMERGE,ITYPE,VOLMER,SOUREL,VOLT2,
     >                    FLXNEW,FLXOLD,RATECM(1,IGR,15),OLDRAT)
            ENDIF
          ENDIF
          CXSNAM='DIFF'
          CALL LCMLEN(KPEDIT,CXSNAM,ILCMLN,ILCMTY)
          IF(ILCMLN.GT.0) THEN
            CALL LCMGET(KPEDIT,CXSNAM,OLDRAT)
            WRITE(IUNOUT,6003) CXSNAM
            CALL EDISTA(IPRINT,NMERGE,ITYPE,VOLMER,SOUREL,VOLT2,
     >                  FLXNEW,FLXOLD,RATECM(1,IGR,NW+5),OLDRAT)
          ENDIF
          CXSNAM='TRANC'
          CALL LCMLEN(KPEDIT,CXSNAM,ILCMLN,ILCMTY)
          IF(ILCMLN.GT.0) THEN
            CALL LCMGET(KPEDIT,CXSNAM,OLDRAT)
            WRITE(IUNOUT,6003) CXSNAM
            CALL EDISTA(IPRINT,NMERGE,ITYPE,VOLMER,SOUREL,VOLT2,
     >                  FLXNEW,FLXOLD,RATECM(1,IGR,NW+9),OLDRAT)
          ENDIF
          DO 135 IL=1,NL
            WRITE (CM,'(I2.2)') IL-1
            CALL LCMLEN(KPEDIT,'NJJS'//CM,ILCMLN,ILCMTY)
            CXSNAM='SCATTERING'//CM
            IF(ILCMLN.EQ.NMERGE) THEN
              WRITE(IUNOUT,6003) CXSNAM
              CALL LCMGET(KPEDIT,'NJJS'//CM,INGSCT)
              CALL LCMGET(KPEDIT,'IJJS'//CM,IFGSCT)
              CALL LCMGET(KPEDIT,'IPOS'//CM,IPOSCT)
              CALL LCMGET(KPEDIT,'SCAT'//CM,XSCAT)
              CALL EDIDEL(IPRINT,NGCOND,NMERGE,IGR,SCATTS(1,1,1,IL),
     >                 INGSCT,IFGSCT,IPOSCT,XSCAT,DELSC)
            ENDIF
 135      CONTINUE
        ENDIF
 210  CONTINUE
      CALL LCMSIX(IPEDIT,' ',ILCMDN)
      CALL LCMSIX(IPEDIT,' ',ILCMDN)
      IF(NSTATS.EQ.-1) THEN
        DEALLOCATE(DELSC,XSCAT)
        DEALLOCATE(IPOSCT,IFGSCT,INGSCT)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(OLDRAT,FLXOLD,FLXNEW)
      RETURN
*----
*  FORMAT
*----
 7000 FORMAT(' ******  EDIDST  WARNING ',
     >'ROUTINE  ******'/' ******  ',A12,' REFERENCE EXECUTION',
     >' DIRECTORY NOT ON LCM  ******'/'************************',
     >'**********************************************')
 7001 FORMAT(' ******  EDIDST  WARNING ',
     >'ROUTINE  ******'/' ******  NUMBER OF REGION OR NUMBER OF ',
     >'GROUP INCONSISTENT        ******'/
     >6X,I10,' CURRENT REGIONS'/6X,I10,' OLD REGIONS'/6X,I10,' CURRENT',
     >' GROUPS'/6X,I10,' OLD GROUPS'/'************************',
     >'**********************************************')
 7002 FORMAT(' ******  EDIDST  WARNING ',
     >'ROUTINE  ******'/' ******  RELATIVE ERROR IN VOLUME: CURRENT ',
     >'TO REFERENCE TOO LARGE ******'/
     >6X,' CURRENT VOLUME/REFERENCE VOLUME  =',1P,E12.4/
     >6X,' REGION NUMERO   =',I10/6X,' RELATIVE ERROR  =',E12.4/
     >6X,' ERROR ALLOWED   =',E12.4/'************************',
     >'**********************************************')
 7004 FORMAT(' ******  EDIDST  WARNING ',
     >'ROUTINE  ******'/' ******  AVERAGE ABSOLUTE ERROR ON LETHARGY',
     >' WIDTH TOO LARGE      ******'/6X,' RELATIVE ERROR  =',E12.4/
     >6X,' ERROR ALLOWED   =',E12.4/'************************',
     >'**********************************************')
 7005 FORMAT(' **************  EDIDST  WARNING  **************',/
     >       ' TOTAL NEW AND/OR OLD SOURCE IS 0.0 ',/
     >       ' RELATIVE SOURCE NORMALIZATION FACTOR SET TO 1.0',/
     >       '************************************************')
 6000 FORMAT(///20X,'D R A G O N   S T A T I S T I C S'/
     >10X,'LCM REFERENCE CASE NAME : ',A12/
     >10X,'NUMBER OF GROUPS        : ',I10/
     >10X,'NUMBER OF REGIONS       : ',I10/
     >10X,'RELATIVE VOLUMES        : ',1P,E12.4/
     >10X,'RELATIVE SOURCES        : ',1P,E12.4)
 6001 FORMAT(/' ANALYSIS OF GROUP : ',I5)
 6002 FORMAT(/' STATISTICS FOR : ',A12)
 6003 FORMAT(/' DELTA SIGMA FOR : ',A12)
 6010 FORMAT(/
     >10X,'REFERENCE Keff          : ',F15.5/
     >10X,'CURRENT Keff            : ',F15.5/
     >10X,'CHANGE IN Keff          : ',3X,F12.2,' mk')
 6011 FORMAT(1P/
     >10X,'REFERENCE B2            : ',E15.3/
     >10X,'CURRENT B2              : ',E15.3/
     >10X,'CHANGE IN B2            : ',E15.3)
      END
