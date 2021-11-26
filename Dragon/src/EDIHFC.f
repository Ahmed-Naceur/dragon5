*DECK EDIHFC
      SUBROUTINE EDIHFC(IPEDIT,NGROUP,NGCOND,NREGIO,NMERGE,NBISO,
     >                  NDEPL,NREAC,MATCOD,VOLUME,INADPL,ISONAM,ISONRF,
     >                  IPISO,MIX,FLUXES,DEN,IGCOND,IMERGE,RER,EMEVF2,
     >                  EMEVG2,VOLME,IPRINT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Evaluate H-factors using information recovered from the reference
* internal library and store them in the edition macrolib.
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
* NGROUP  number of groups.
* NGCOND  number of condensed groups.
* NREGIO  number of regions.
* NMERGE  number of merged regions.
* NBISO   number of isotopes.
* NDEPL   number of depleting isotopes.
* NREAC   number of depletion reactions.
* MATCOD  material per region.
* VOLUME  volume of region.
* INADPL  name of depleting isotopes.
* ISONAM  isotopes names.
* ISONRF  library name of isotopes.
* IPISO   pointer array towards microlib isotopes.
* MIX     mixture associated with isotopes.
* FLUXES  multigroup fluxes.
* DEN     isotope density.
* IGCOND  limits of condensed groups.
* IMERGE  index of merged region.
* RER     fission and capture production energy (MeV/reaction).
* VOLME   merged volume.
* IPRINT  print level.
*
*Parameters: output
* EMEVF2  fission production energy by isotope.
* EMEVG2  capture production energy by isotope.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPEDIT,IPISO(NBISO)
      INTEGER     IUNOUT
      INTEGER     NGROUP,NGCOND,NREGIO,NMERGE,NBISO,NDEPL,NREAC,
     >            MATCOD(NREGIO),INADPL(3,NDEPL),ISONAM(3,NBISO),
     >            ISONRF(3,NBISO),MIX(NBISO),IGCOND(NGCOND),
     >            IMERGE(NREGIO)
      REAL        VOLUME(NREGIO),FLUXES(NREGIO,NGROUP),DEN(NBISO),
     >            RER(NREAC,NDEPL),EMEVF2(NBISO),EMEVG2(NBISO)
      REAL        VOLME(NMERGE)
      INTEGER     IPRINT
      DOUBLE PRECISION TOTPOW,POWF,POWC,POWT
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INDX
      REAL, ALLOCATABLE, DIMENSION(:) :: SIG,HFACT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: FLXMER
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: WORK
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPEDIT,KPEDIT,KPLIB
      PARAMETER  (IUNOUT=6)
      INTEGER     IGAR(3)
      CHARACTER   HNISOR*12,TEXT12*12,HSMG*131
      LOGICAL     L1,L2
      DOUBLE PRECISION GAR,CONV,XDRCST
*----
*  SCRATCH STORAGE ALLOCATION
*   SIG     fission/capture cross sections.
*   HFACT   H-factor in a macrogroup.
*   FLXMER  merged and condensed flux.
*   WORK    H-factors.
*   INDX    depleting isotope index.
*----
      ALLOCATE(INDX(NBISO))
      ALLOCATE(SIG(NGROUP),HFACT(NMERGE))
      ALLOCATE(FLXMER(NMERGE,NGCOND),WORK(NMERGE,NGCOND,3))
*----
*  COMPUTE THE DEPLETING ISOTOPE INDEX
*----
      DO 20 ISO=1,NBISO
        WRITE(HNISOR,'(3A4)') (ISONRF(I0,ISO),I0=1,3)
        I1=INDEX(HNISOR,'_')
        IF(I1.EQ.0) THEN
          TEXT12=HNISOR
        ELSE
          TEXT12=HNISOR(:I1-1)
        ENDIF
        READ(TEXT12,'(3A4)') (IGAR(I0),I0=1,3)
        DO 10 IDP=1,NDEPL
          L1=((ISONRF(1,ISO).EQ.INADPL(1,IDP)).AND.
     1        (ISONRF(2,ISO).EQ.INADPL(2,IDP)).AND.
     2        (ISONRF(3,ISO).EQ.INADPL(3,IDP)))
          L2=((IGAR(1).EQ.INADPL(1,IDP)).AND.
     1        (IGAR(2).EQ.INADPL(2,IDP)).AND.
     2        (IGAR(3).EQ.INADPL(3,IDP)))
          IF(L1.OR.L2) THEN
            INDX(ISO)=IDP
            GO TO 20
          ENDIF
  10    CONTINUE
        INDX(ISO)=0
  20  CONTINUE
*----
*  COMPUTE H-FACTOR
*----
      CONV=1.0D6*XDRCST('eV','J')
      IZFISS=0
      CALL XDDSET(WORK,NMERGE*NGCOND*3,0.0D0)
      CALL XDDSET(FLXMER,NMERGE*NGCOND,0.0D0)
      DO 160 ISO=1,NBISO
        IDPL=INDX(ISO)
        IF(IDPL.EQ.0) GO TO 160
        KPLIB=IPISO(ISO) ! set ISO-th isotope
        IF(.NOT.C_ASSOCIATED(KPLIB)) THEN
          WRITE(HSMG,'(17HEDIHFC: ISOTOPE '',3A4,16H'' IS NOT AVAILAB,
     >    19HLE IN THE MICROLIB.)') (ISONAM(I0,ISO),I0=1,3)
          CALL XABORT(HSMG)
        ENDIF
*----
*  COMPUTE FISSION ENERGY
*----
        CALL LCMLEN(KPLIB,'NFTOT',ILLCM,ITLCM)
        IF(ILLCM.EQ.NGROUP) THEN
          IZFISS=IZFISS+1
          EMEVF2(ISO)=RER(2,IDPL)
          CALL LCMGET(KPLIB,'NFTOT',SIG)
          DO 120 IREG=1,NREGIO
            IMR=IMERGE(IREG)
            IF((IMR.GT.0).AND.(MATCOD(IREG).EQ.MIX(ISO))) THEN
              IGRFIN=0
              DO 110 IGC=1,NGCOND
                IGRDEB=IGRFIN+1
                IGRFIN=IGCOND(IGC)
                GAR=0.0D0
                DO 100 IGR=IGRDEB,IGRFIN
                  GAR=GAR+FLUXES(IREG,IGR)*DEN(ISO)*VOLUME(IREG)*
     >            SIG(IGR)
 100            CONTINUE
                WORK(IMR,IGC,1)=WORK(IMR,IGC,1)+GAR*RER(2,IDPL)*CONV
 110          CONTINUE
            ENDIF
 120      CONTINUE
        ENDIF
*----
*  COMPUTE CAPTURE ENERGY
*----
        CALL LCMLEN(KPLIB,'NG',ILLCM,ITLCM)
        IF(ILLCM.EQ.NGROUP) THEN
          IZFISS=IZFISS+1
          EMEVG2(ISO)=RER(3,IDPL)
          CALL LCMGET(KPLIB,'NG',SIG)
          DO 150 IREG=1,NREGIO
            IMR=IMERGE(IREG)
            IF((IMR.GT.0).AND.(MATCOD(IREG).EQ.MIX(ISO))) THEN
              IGRFIN=0
              DO 140 IGC=1,NGCOND
                IGRDEB=IGRFIN+1
                IGRFIN=IGCOND(IGC)
                GAR=0.0D0
                DO 130 IGR=IGRDEB,IGRFIN
                  GAR=GAR+FLUXES(IREG,IGR)*DEN(ISO)*VOLUME(IREG)*
     >            SIG(IGR)
 130            CONTINUE
                WORK(IMR,IGC,2)=WORK(IMR,IGC,2)+GAR*RER(3,IDPL)*CONV
 140          CONTINUE
            ENDIF
 150      CONTINUE
        ENDIF
 160  CONTINUE
*----
*  Normalize total power to 1 W
*  Print fission, capture and total power density
*----
      TOTPOW=0.0D0
      DO IGC=1,NGCOND
        DO IMR=1,NMERGE
          WORK(IMR,IGC,3)=WORK(IMR,IGC,1)+WORK(IMR,IGC,2)
          TOTPOW=TOTPOW+WORK(IMR,IGC,3)
        ENDDO
      ENDDO
      IF(TOTPOW.GT.0.0D0) THEN
        IF(ABS(IPRINT).GE.1) THEN
          WRITE(IUNOUT,6000)
          DO IMR=1,NMERGE
            POWF=0.0D0
            POWC=0.0D0
            POWT=0.0D0
            DO IGC=1,NGCOND
              POWF=POWF+WORK(IMR,IGC,1)
              POWC=POWC+WORK(IMR,IGC,2)
              POWT=POWT+WORK(IMR,IGC,3)
            ENDDO
            IF(VOLME(IMR).NE.0.0) THEN
              POWF=POWF/(TOTPOW*VOLME(IMR))
              POWC=POWC/(TOTPOW*VOLME(IMR))
              POWT=POWT/(TOTPOW*VOLME(IMR))
              WRITE(IUNOUT,6001) IMR,VOLME(IMR),POWF,POWC,POWT
            ENDIF
          ENDDO
        ENDIF
      ENDIF
*----
*  COMPUTE THE HOMOGENIZED/CONDENSED FLUX
*----
      IF(IZFISS.NE.0) THEN
        DO 190 IREG=1,NREGIO
          IMR=IMERGE(IREG)
          IF(IMR.GT.0) THEN
            IGRFIN=0
            DO 180 IGC=1,NGCOND
              IGRDEB=IGRFIN+1
              IGRFIN=IGCOND(IGC)
              GAR=0.0D0
              DO 170 IGR=IGRDEB,IGRFIN
                GAR=GAR+FLUXES(IREG,IGR)*VOLUME(IREG)
 170          CONTINUE
              FLXMER(IMR,IGC)=FLXMER(IMR,IGC)+GAR
 180        CONTINUE
          ENDIF
 190    CONTINUE
        DO 210 IGC=1,NGCOND
          DO 200 IMR=1,NMERGE
            IF(FLXMER(IMR,IGC).GT.0.0) THEN
              WORK(IMR,IGC,3)=WORK(IMR,IGC,3)/FLXMER(IMR,IGC)
            ENDIF
 200      CONTINUE
 210    CONTINUE
      ENDIF
*----
*  SAVE ON LCM
*----
      CALL LCMSIX(IPEDIT,'MACROLIB',1)
      JPEDIT=LCMLID(IPEDIT,'GROUP',NGCOND)
      DO 230 IGC=1,NGCOND
        DO 220 IMR=1,NMERGE
        HFACT(IMR)=REAL(WORK(IMR,IGC,3))
 220    CONTINUE
        KPEDIT=LCMDIL(JPEDIT,IGC)
        CALL LCMPUT(KPEDIT,'H-FACTOR',NMERGE,2,HFACT)
 230  CONTINUE
      CALL LCMSIX(IPEDIT,' ',2)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(WORK,FLXMER)
      DEALLOCATE(HFACT,SIG)
      DEALLOCATE(INDX)
*----
*  FORMAT
*----
 6000 FORMAT(/' POWER DENSITY (W/cc) NORMALIZED TO 1 W TOTAL POWER '/
     >1X,'REGION',6X,'VOLUME',7X,'FISSION',7X,'CAPTURE',9X,'TOTAL')
 6001 FORMAT(1X,I4,4F14.8)
      RETURN
      END
