*DECK EDIACT
      SUBROUTINE EDIACT(IPEDIT,IPRINT,NGROUP,NGCOND,NREGIO,NMERGE,NL,
     >                  NBISO,NED,VOLUME,MIX,IGCOND,IMERGE,FLUXES,
     >                  ITRANC,ISONAM,IPISO,HVECT,CURNAM,NACTI,IACTI,
     >                  EMEVF2,EMEVG2)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Homogenization and condensation of activation cross sections.
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
* IPRINT  print index.
* NGROUP  number of energy groups.
* NGCOND  number of condensed groups.
* NREGIO  number of volumes.
* NMERGE  number of merged regions.
* NL      number of legendre orders required in the calculation
*         (NL=1 or higher).
* NBISO   number of isotopes.
* NED     number of extra vector edits.
* VOLUME  volumes.
* MIX     mixture number associated with each isotope.
* IGCOND  limits of condensed groups.
* IMERGE  index of merged regions.
* FLUXES  fluxes.
* ITRANC  transport correction type (0 -> no transport correction).
* ISONAM  names of the isotopes to be treated.
* IPISO   pointer array towards microlib isotopes.
* HVECT   names of the extra vector edits.
* CURNAM  name of the lcm directory where the microscopic cross
*         sections are stored (blank name implies no save).
* NACTI   number of mixture with WIMS activation edit.
* IACTI   mixtures with activation edits.
* EMEVF2  fission production energy.
* EMEVG2  capture production energy.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPEDIT,IPISO(NBISO)
      CHARACTER   HVECT(*)*8,CURNAM*(*)
      INTEGER     IPRINT,NGROUP,NGCOND,NREGIO,NMERGE,NL,NBISO,NED,
     >            MIX(NBISO),IGCOND(NGCOND),IMERGE(NREGIO),ITRANC,
     >            ISONAM(3,NBISO),NACTI,IACTI(NACTI)
      REAL        VOLUME(NREGIO),FLUXES(NREGIO,NGROUP),
     >            EMEVF2(NBISO),EMEVG2(NBISO)
*----
*  LOCAL VARIABLES
*----
      PARAMETER  (IOUT=6,NSTATE=40)
      TYPE(C_PTR) KPLIB
      INTEGER     IPAR(NSTATE)
      CHARACTER   CACTI*12,CM*2,HMAKE(100)*8,HNEW*12,TEXT12*12,HSMG*131
      LOGICAL     LMEVF,LMEVG,LLCM
      DOUBLE PRECISION DVOL,DFLI,DTMP,QEN,ERR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ISOMIX
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: KCJJ,HNISO
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: KFJJ
      REAL, ALLOCATABLE, DIMENSION(:) :: CXSV,CSCAT
      REAL, ALLOCATABLE, DIMENSION(:,:) :: RXSV,RSCAT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: DFLX,DXSV
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DSCAT
*----
*  SCRATCH STORAGE ALLOCATION
*   RXSV    real microscopic cross section/flux (vector- full group
*           structure):
*           RXSV(ig,1->NL) = total scattering order 0 to NL-1;
*           RXSV(ig,1+NL)  = total xs;
*           RXSV(ig,2+NL)  = nusigf;
*           RXSV(ig,3+NL->2+NL+NED) = additional xs;
*           RXSV(ig,3+NL+NED)= tranc;
*           RXSV(ig,4+NL+NED)= chi.
*   KFJJ    scattering vector index (vector- full group structure).
*   RSCAT   real microscopic scattering x-s (vector- full group
*           structure):
*   DFLX    double flux.
*   DXSV    double microscopic reaction rates (vector- condensed group
*           structure).
*   DSCAT   microscopic scattering rate (vector- condensed group
*           structure).
*   CXSV    real microscopic cross section/flux (vector- condensed
*           group structure).
*   KCJJ    scattering vector index (vector- condensed group structure).
*   CSCAT   real microscopic scattering rate (vector- condensed group
*           structure).
*   HNISO   isotope name vector.
*   ISOMIX  mixture number associated with new isotope.
*----
      ALLOCATE(ISOMIX(NACTI*NMERGE*NBISO),KFJJ(NGROUP,3,NL),
     > KCJJ(NGCOND,2),HNISO(3,NACTI*NMERGE*NBISO))
      ALLOCATE(RXSV(NGROUP,NL+NED+4),RSCAT(NGROUP*NGROUP,NL),
     > CXSV(NGCOND),CSCAT(NGCOND*NGCOND))
      ALLOCATE(DFLX(NGCOND,NMERGE),DXSV(NGCOND,NL+NED+4),
     > DSCAT(NGCOND,NGCOND,NL))
*----
*  EVALUATE INTEGRATED FLUX
*----
      DO 10 INM=1,NMERGE
        DO 11 IGRCND=1,NGCOND
          DFLX(IGRCND,INM)=0.0D0
 11     CONTINUE
        DO 20 IREGIO=1,NREGIO
          IF(IMERGE(IREGIO).EQ.INM) THEN
            IGRFIN=0
            DO 21 IGRCND=1,NGCOND
              IGRDEB=IGRFIN+1
              IGRFIN=IGCOND(IGRCND)
              DTMP=0.0D0
              DO 22 IGR=IGRDEB,IGRFIN
                DTMP=DTMP+DBLE(FLUXES(IREGIO,IGR))
 22           CONTINUE
              DFLX(IGRCND,INM)=DFLX(IGRCND,INM)+
     >          DTMP*DBLE(VOLUME(IREGIO))
 21         CONTINUE
          ENDIF
 20     CONTINUE
 10   CONTINUE
*----
*  LOOP OVER EACH MIXTURE WITH ACTIVATION EDIT
*  FIND ISOTOPES ASSOCIATED WITH THIS MIXTURE
*----
      LLCM=CURNAM.NE.' '
      MAXH=4+NL+NED
      DO 100 IRE=1,NACTI
        IMIXR=IACTI(IRE)
        WRITE(CACTI,'(8HACTIVITY,I4)') IRE
        IF(IPRINT.GT.0) WRITE(IOUT,300) IMIXR,CACTI
        JJISO=0
        DO 110 ISO=1,NBISO
          IF(MIX(ISO).EQ.IMIXR) THEN
            IF(IPRINT.GT.0) WRITE(IOUT,310) (ISONAM(I0,ISO),I0=1,2)
*----
*  THIS ISOTOPE IS ASSOCIATED WITH AN ACTIVATION MIXTURE
*  READ MICROSCOPIC CROSS SECTIONS 'SIGS'//CM, 'SCAT'//CM, 'NTOT0',
* 'NUSIGF', 'CHI', HVECT.
*----
            DO 114 INAM=1,MAXH+NL
              HMAKE(INAM)=' '
 114        CONTINUE
        KPLIB=IPISO(ISO) ! set ISO-th isotope
        IF(.NOT.C_ASSOCIATED(KPLIB)) THEN
          WRITE(HSMG,'(17HEDIACT: ISOTOPE '',3A4,16H'' IS NOT AVAILAB,
     >    19HLE IN THE MICROLIB.)') (ISONAM(I0,ISO),I0=1,3)
          CALL XABORT(HSMG)
        ENDIF
            CALL LCMGET(KPLIB,'AWR',AWR)
            IF(EMEVF2(ISO).GT.0.0) EVF=EMEVF2(ISO)
            CALL LCMLEN(KPLIB,'MEVF',LENGTH,ITYLCM)
            IF(LENGTH.EQ.1) CALL LCMGET(KPLIB,'MEVF',EVF)
            LMEVF=(LENGTH.EQ.1).OR.(EMEVF2(ISO).GT.0.0)
            IF(EMEVG2(ISO).GT.0.0) EVG=EMEVG2(ISO)
            CALL LCMLEN(KPLIB,'MEVG',LENGTH,ITYLCM)
            IF(LENGTH.EQ.1) CALL LCMGET(KPLIB,'MEVG',EVG)
            LMEVG=(LENGTH.EQ.1).OR.(EMEVG2(ISO).GT.0.0)
            DO 111 IL=1,NL
              WRITE (CM,'(I2.2)') IL-1
              CALL LCMLEN(KPLIB,'SIGS'//CM,LENGTH,ITYLCM)
              IF(LENGTH.EQ.NGROUP) THEN
                CALL LCMGET(KPLIB,'SIGS'//CM,RXSV(1,IL))
                HMAKE(IL)='SIGS'//CM
              ELSE
                HMAKE(IL)=' '
              ENDIF
              CALL LCMLEN(KPLIB,'NJJS'//CM,LENGTH,ITYLCM)
              IF(LENGTH.EQ.NGROUP) THEN
                CALL LCMGET(KPLIB,'NJJS'//CM,KFJJ(1,1,IL))
                CALL LCMGET(KPLIB,'IJJS'//CM,KFJJ(1,2,IL))
                CALL LCMGET(KPLIB,'SCAT'//CM,RSCAT(1,IL))
                HMAKE(MAXH+IL)=CM
                IPO=0
                DO 112 IGR=1,NGROUP
                  KFJJ(IGR,3,IL)=IPO+1
                  IPO=IPO+KFJJ(IGR,1,IL)
 112            CONTINUE
              ELSE
                HMAKE(MAXH+IL)=' '
              ENDIF
 111        CONTINUE
            CALL LCMGET(KPLIB,'NTOT0',RXSV(1,1+NL))
            HMAKE(1+NL)='NTOT0'
            CALL LCMLEN(KPLIB,'NUSIGF',LENGTH,ITYLCM)
            IF(LENGTH.EQ.NGROUP) THEN
              CALL LCMGET(KPLIB,'NUSIGF',RXSV(1,2+NL))
              HMAKE(2+NL)='NUSIGF'
            ELSE
              HMAKE(2+NL)=' '
            ENDIF
            CALL LCMLEN(KPLIB,'CHI',LENGTH,ITYLCM)
            IF(LENGTH.EQ.NGROUP) THEN
              CALL LCMGET(KPLIB,'CHI',RXSV(1,MAXH))
              HMAKE(MAXH)='CHI'
            ELSE
              HMAKE(MAXH)=' '
            ENDIF
            DO 113 IED=1,NED
              CALL LCMLEN(KPLIB,HVECT(IED),LENGTH,ITYLCM)
              IF(LENGTH.GT.0) THEN
                CALL LCMGET(KPLIB,HVECT(IED),RXSV(1,2+NL+IED))
                HMAKE(2+NL+IED)=HVECT(IED)
              ELSE
                HMAKE(2+NL+IED)=' '
              ENDIF
 113        CONTINUE
            IF(LLCM) THEN
              CALL LCMSIX(IPEDIT,CURNAM,1)
              CALL LCMSIX(IPEDIT,CACTI,1)
            ENDIF
            DO 120 INM=1,NMERGE
              DVOL=0.0D0
              JJISO=JJISO+1
              DO 121 J=1,MAXH
                DO 122 I=1,NGCOND
                  DXSV(I,J)=0.0D0
 122            CONTINUE
 121          CONTINUE
              DO 123 K=1,NL
                DO 124 J=1,NGCOND
                  DO 125 I=1,NGCOND
                    DSCAT(I,J,K)=0.0D0
 125              CONTINUE
 124            CONTINUE
 123          CONTINUE
*----
*  MERGE/CONDENSE REACTIONS 'SIGS'//CM, 'SCAT'//CM, 'NTOT0',
*  'NUSIGF', 'CHI', AND HVECT.
*----
              DO 130 IREGIO=1,NREGIO
                VOL=VOLUME(IREGIO)
                IF(IMERGE(IREGIO).EQ.INM) THEN
                  DVOL=DVOL+DBLE(VOL)
                  IGRFIN=0
                  DO 150 IGRCND=1,NGCOND
                    IGRDEB=IGRFIN+1
                    IGRFIN=IGCOND(IGRCND)
                    DO 151 IGR=IGRDEB,IGRFIN
                      DFLI=DBLE(FLUXES(IREGIO,IGR)*VOL)
                      DO 160 J=1,MAXH-2
                        IF(HMAKE(J).NE.' ') THEN
                          DXSV(IGRCND,J)=DXSV(IGRCND,J)
     >                     +DBLE(RXSV(IGR,J))*DFLI
                        ENDIF
 160                  CONTINUE
                      DO 152 IL=1,NL
                        IF(HMAKE(MAXH+IL).NE.' ') THEN
*----
*  IGRCND IS THE SECONDARY GROUP.
*----
                          NGSCAT=KFJJ(IGR,1,IL)
                          IGSCAT=KFJJ(IGR,2,IL)
                          JGRFIN=0
                          DO 170 JGRCND=1,NGCOND
*----
*  JGRCND IS THE PRIMARY GROUP.
*----
                            JGRDEB=JGRFIN+1
                            JGRFIN=IGCOND(JGRCND)
                            J2=MIN(JGRFIN,IGSCAT)
                            J1=MAX(JGRDEB,IGSCAT-NGSCAT+1)
                            DTMP=0.0D0
                            IPO=KFJJ(IGR,3,IL)+IGSCAT-J2
                            DO 171 JGR=J2,J1,-1
                              DTMP=DTMP+DBLE(RSCAT(IPO,IL)*
     >                          FLUXES(IREGIO,JGR)*VOL)
                              IPO=IPO+1
 171                        CONTINUE
                            DSCAT(JGRCND,IGRCND,IL)=
     >                        DSCAT(JGRCND,IGRCND,IL)+DTMP
 170                      CONTINUE
                          IF((ITRANC.NE.0).AND.(IL.EQ.2)) THEN
*----
*  INFO USED BY WIMS TYPE TRANSPORT CORRECTION.
*----
                            HMAKE(MAXH-1)='TRANC'
                            DTMP=DBLE(RXSV(IGR,IL))
                            DXSV(IGRCND,MAXH-1)=DXSV(IGRCND,MAXH-1)
     >                        +DTMP*DFLI
                          ENDIF
                        ENDIF
 152                  CONTINUE
 151                CONTINUE
 150              CONTINUE
                ENDIF
 130          CONTINUE
              WRITE(HNEW,'(2A4,I4.4)') (ISONAM(I0,ISO),I0=1,2),INM
              READ(HNEW,'(3A4)') (HNISO(I1,JJISO),I1=1,3)
              ISOMIX(JJISO)=INM
              IF(IPRINT.GT.0) WRITE(IOUT,320) INM,HNEW
*----
*  EVALUATE FEWGROUPS CHI
*----
              IF(HMAKE(MAXH).NE.' ') THEN
                IGRFIN=0
                DO 191 IGRCND=1,NGCOND
                  IGRDEB=IGRFIN+1
                  IGRFIN=IGCOND(IGRCND)
                  DO 192 IGR=IGRDEB,IGRFIN
                    DXSV(IGRCND,MAXH)=DXSV(IGRCND,MAXH)
     >                +DBLE(RXSV(IGR,MAXH))
 192              CONTINUE
 191            CONTINUE
              ENDIF
*----
*  EVALUATE FEWGROUPS MICROSCOPIC XS FOR ACTIVATION ISOTOPES
*----
              DO 200 IGRCND=1,NGCOND
                DO 210 IL=1,NL
                  DTMP=DXSV(IGRCND,IL)
                  DO 211 JGRCND=1,NGCOND
                    IF(JGRCND.NE.IGRCND) THEN
                      DTMP=DTMP-DSCAT(IGRCND,JGRCND,IL)
                    ENDIF
 211              CONTINUE
                  QEN=MAX(ABS(DTMP),ABS(DSCAT(IGRCND,IGRCND,IL)))
                  IF(QEN.GT.0.0D0) THEN
                    ERR=ABS(DTMP-DSCAT(IGRCND,IGRCND,IL))/QEN
                    IF(ERR.GT.1.0D-3) THEN
                      WRITE(IOUT,340) IGRCND,IL-1,100.0*ERR
                    ENDIF
                    DSCAT(IGRCND,IGRCND,IL)=DTMP
                  ENDIF
                  DO 212 JGRCND=1,NGCOND
                    IF(DFLX(IGRCND,INM).GT.0.0D0) THEN
                      DSCAT(IGRCND,JGRCND,IL)=
     >                  DSCAT(IGRCND,JGRCND,IL)/DFLX(IGRCND,INM)
                    ELSE
                      DSCAT(IGRCND,JGRCND,IL)=0.0D0
                    ENDIF
 212              CONTINUE
 210            CONTINUE
                DO 213 J=1,MAXH-1
                  IF(DFLX(IGRCND,INM).GT.0.0D0) THEN
                    DXSV(IGRCND,J)=DXSV(IGRCND,J)/DFLX(IGRCND,INM)
                  ELSE
                    DXSV(IGRCND,J)=0.0D0
                  ENDIF
 213            CONTINUE
 200          CONTINUE
              IF(LLCM) THEN
                CALL LCMSIX(IPEDIT,HNEW,1)
                CALL LCMPUT(IPEDIT,'AWR',1,2,AWR)
                IF(LMEVF) CALL LCMPUT(IPEDIT,'MEVF',1,2,EVF)
                IF(LMEVG) CALL LCMPUT(IPEDIT,'MEVG',1,2,EVG)
                DO 220 J=1,MAXH
                  IF(HMAKE(J).NE.' ') THEN
                    DO 221 IGCD=1,NGCOND
                      CXSV(IGCD)=REAL(DXSV(IGCD,J))
 221                CONTINUE
                    CALL LCMPUT(IPEDIT,HMAKE(J),NGCOND,2,CXSV)
                  ENDIF
 220            CONTINUE
                DO 230 IL=1,NL
                  IF(HMAKE(MAXH+IL).NE.' ') THEN
                    KGAR=0
                    DO 231 IG2=1,NGCOND
                      IGMIN=IG2
                      IGMAX=IG2
                      DO 232 IG1=NGCOND,1,-1
                        IF(DSCAT(IG1,IG2,IL).NE.0.0D0) THEN
                          IGMIN=MIN(IGMIN,IG1)
                          IGMAX=MAX(IGMAX,IG1)
                        ENDIF
 232                  CONTINUE
                      KCJJ(IG2,1)=IGMAX-IGMIN+1
                      KCJJ(IG2,2)=IGMAX
                      DO 233 IG1=IGMAX,IGMIN,-1
                        KGAR=KGAR+1
                        CSCAT(KGAR)=REAL(DSCAT(IG1,IG2,IL))
 233                  CONTINUE
 231                CONTINUE
                    CM=HMAKE(MAXH+IL)(:2)
                    CALL LCMPUT(IPEDIT,'NJJS'//CM,NGCOND,1,KCJJ(1,1))
                    CALL LCMPUT(IPEDIT,'IJJS'//CM,NGCOND,1,KCJJ(1,2))
                    CALL LCMPUT(IPEDIT,'SCAT'//CM,KGAR,2,CSCAT)
                  ENDIF
 230            CONTINUE
                CALL LCMSIX(IPEDIT,' ',2)
              ENDIF
              IF(IPRINT.GT.3) THEN
                WRITE(IOUT,330) 'FLXAVG',(DFLX(I,INM)/DVOL,I=1,NGCOND)
                DO 240 J=1,MAXH
                  IF(HMAKE(J).NE.' ') THEN
                     WRITE(IOUT,330) HMAKE(J),(DXSV(I,J),I=1,NGCOND)
                  ENDIF
 240            CONTINUE
                WRITE(IOUT,330) 'SIGA',
     >            (DXSV(I,1+NL)-DXSV(I,1),I=1,NGCOND)
                WRITE(IOUT,330) 'SIGW00',(DSCAT(I,I,1),I=1,NGCOND)
                IF(NL.GT.1) THEN
                  IF(HMAKE(MAXH+2).NE.' ')
     >            WRITE (6,330) 'SIGW01',(DSCAT(I,I,2),I=1,NGCOND)
                ENDIF
              ENDIF
 120        CONTINUE
            IF(LLCM) THEN
              CALL LCMSIX(IPEDIT,' ',2)
              CALL LCMSIX(IPEDIT,' ',2)
            ENDIF
          ENDIF
 110    CONTINUE
        IF(JJISO.GT.0.AND.LLCM) THEN
          CALL LCMSIX(IPEDIT,CURNAM,1)
          CALL LCMSIX(IPEDIT,CACTI,1)
          TEXT12='L_LIBRARY'
          CALL LCMPTC(IPEDIT,'SIGNATURE',12,1,TEXT12)
          DO 105 I=1,NSTATE
            IPAR(I)=0
 105      CONTINUE
          IPAR(1)=NMERGE
          IPAR(2)=JJISO
          IPAR(3)=NGCOND
          IPAR(4)=NL
          IPAR(5)=ITRANC
          IF(ITRANC.NE.0) IPAR(5)=2
          IPAR(7)=1
          IPAR(13)=NED
          IPAR(14)=NACTI
          CALL LCMPUT(IPEDIT,'STATE-VECTOR',NSTATE,1,IPAR)
          IF(NED.GT.0) CALL LCMPTC(IPEDIT,'ADDXSNAME-P0',8,NED,HVECT)
          CALL LCMPUT(IPEDIT,'ISOTOPESUSED',3*JJISO,3,HNISO)
          CALL LCMPUT(IPEDIT,'ISOTOPESMIX',JJISO,1,ISOMIX)
          CALL LCMSIX(IPEDIT,' ',2)
          CALL LCMSIX(IPEDIT,' ',2)
        ENDIF
 100  CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(DSCAT,DXSV,DFLX)
      DEALLOCATE(CSCAT,CXSV,RSCAT,RXSV)
      DEALLOCATE(HNISO,KCJJ,KFJJ,ISOMIX)
      RETURN
*
 300  FORMAT(//' MICROSCOPIC ACTIVITY XS FOR MATERIAL NUMBER : ',I5/
     >         '                     STORED ON SUB-DIRECTORY : ',A12)
 310  FORMAT(/24X,'ISOTOPE NAME PREFIX  : ',2A4)
 320  FORMAT(31X,'REGION NUMBER : ',I5,5X,'FINAL ISOTOPE NAME : ',A12)
 330  FORMAT(' XS TYPE  ',A8/(1X,1P,10E12.4))
 340  FORMAT(' EDIACT: *** WARNING *** NORMALIZATION OF THE WITHIN-',
     > 'GROUP SCATTERING TRANSFER IN GROUP',I4,' AND ORDER',I3,' BY',
     > F6.2,' %.')
      END
