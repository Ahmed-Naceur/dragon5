*DECK EDIPXS
      SUBROUTINE EDIPXS(IPEDIT,IADJ,IPRINT,NL,NDEL,NALBP,ITRANC,NSAVES,
     >                  NGCOND,NMERGE,ILEAKS,NW,NTAUXT,EIGENK,B2,
     >                  CUREIN,NIFISS,CURNAM,NEDMAC,VOLMER,WLETYC,
     >                  WENERG,SCATTD,RATECM,FLUXCM,FADJCM,SIGS,SCATTS,
     >                  DISFCT,ALBP,TAUXE,HVECT,OVERV,HFACT,NENER,TIMEF,
     >                  LH)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Save homogenized/condensed macroscopic cross sections.
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
* IADJ    type of flux weighting:
*         = 0 direct flux weighting;
*         = 1 direct-adjoint flux weighting.
* IPRINT  print level;
*         = 0 no print;
*         = 1 print fluxes;
*         = 2 1+print reaction rates;
*         = 3 2+print homogenized cross sections.
* NL      number of Legendre orders.
* NDEL    number of delayed precursor groups.
* NALBP   number of physical albedos.
* ITRANC  type of transport correction.
* NSAVES  homogenized cross section compute/save flag:
*         = 0  no compute, no save;
*         = 1  compute, no save;
*         = 2  compute and save.
* NGCOND  number of groups condensed.
* NMERGE  number of regions merged.
* ILEAKS  type of leakage calculation:
*         = 0 no leakage;
*         = 1 homogeneous leakage (Diffon);
*         = 2 isotropic streaming (Ecco);
*         = 3 anisotropic streaming (Tibere);
*         = 10 isotropic diffusion coefficients recovered from input
*           macrolib;
*         = 11 anisotropic diffusion coefficients recovered from input
*           macrolib.
* NW      type of weighting for PN cross section info (=0 P0; =1 P1).
* NTAUXT  number of reaction rate edits (=15+2*NDEL).
* EIGENK  eigenvalue for problem.
* B2      square buckling:
*         for ILEAKS=1,2: B2(4) is homogeneous;
*         for ILEAKS=3: B2(1),B2(2),B2(3) are directional heterogeneous
*         and B2(4) is homogeneous.
* CUREIN  infinite multiplication factor.
* NIFISS  number of fissile isotopes.
* CURNAM  name of LCM directory where the merged/condensed cross
*         sections are stored.
* NEDMAC  number of extra edit vectors.
* VOLMER  volume of region merged.
* WLETYC  lethargy width condensed.
* WENERG  energy group limits.
* SCATTD  double precision scattering rates.
* NENER   number of energy groups limits.
* TIMEF   time stamp in day/burnup/irradiation.
* LH      flag set to true if H-factors are set.
*
*Parameters: output
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
* FADJCM  averaged region/group afjoint fluxes:
*         = FADJCM(*,1) = adjoint fluxes P0;
*         = FADJCM(*,2) = adjoint fluxes P1.
* SIGS    Legendre dependent scattering cross sections.
* SCATTS  homogenized scattering cross sections.
* DISFCT  disadvantage factor.
* ALBP    physical albedos.
* TAUXE   extra edit rates.
* HVECT   extra edit names.
* OVERV   1/v merge condensed.
* HFACT   H-factors condensed.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPEDIT
      INTEGER     IADJ,IPRINT,NL,NDEL,NALBP,ITRANC,NSAVES,NGCOND,NMERGE,
     >            ILEAKS,NW,NTAUXT,NIFISS,NEDMAC,NENER
      REAL        EIGENK,B2(4),CUREIN,VOLMER(NMERGE),WLETYC(NGCOND),
     >            WENERG(NGCOND+1),RATECM(NMERGE,NGCOND,NTAUXT),
     >            FLUXCM(NMERGE,NGCOND,NW+1),FADJCM(NMERGE,NGCOND,NW+1),
     >            SIGS(NMERGE,NGCOND,NL),
     >            SCATTS(NMERGE,NGCOND,NGCOND,NL),DISFCT(NGCOND),
     >            ALBP(NALBP,NGCOND,NGCOND),TAUXE(NMERGE,NGCOND,NEDMAC),
     >            OVERV(NMERGE,NGCOND),HFACT(NMERGE,NGCOND),TIMEF(3)
      LOGICAL     LH
      CHARACTER   CURNAM*12,HVECT(NEDMAC)*8
      DOUBLE PRECISION SCATTD(NMERGE,NGCOND,NGCOND,NL)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPEDIT,KPEDIT
      CHARACTER   APG*3
      PARAMETER  (IUNOUT=6,APG=' > ',ILCMUP=1,ILCMDN=2,NSTATE=40)
      CHARACTER   CEDNAM*12,HSIGN*12,CM*2
      INTEGER     IDATA(NSTATE)
      DOUBLE PRECISION SCATWG,SCATTN
      LOGICAL     LAL1D
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) ::IJJ,NJJ,IPOS
      REAL, ALLOCATABLE, DIMENSION(:) ::SCATC,DIFF
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FACT,ALB1
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IJJ(NMERGE),NJJ(NMERGE),IPOS(NMERGE))
      ALLOCATE(SCATC(NMERGE*NGCOND),FACT(NMERGE,NW+1),DIFF(NGCOND))
*----
*  COMPUTE MERGED/CONDENSED X-S
*----
      IF(NSAVES.GE.1) THEN
        IDATA(4)=0
        DO 200 IGR=1,NGCOND
          DO 40 IKK=1,NMERGE
            DO 5 IL=1,NW+1
              IF(FLUXCM(IKK,IGR,1).EQ.0.0) THEN
                FACT(IKK,IL)=0.0
              ELSE
                FACT(IKK,IL)=1.0/FLUXCM(IKK,IGR,IL)
              ENDIF
    5       CONTINUE
            RATECM(IKK,IGR,NW+3)=RATECM(IKK,IGR,NW+3)*FACT(IKK,1)
            IF((RATECM(IKK,IGR,NW+3).NE.0.0).OR.
     >         (RATECM(IKK,IGR,NW+8).NE.0.0)) IDATA(4)=1
            IF(IADJ.EQ.0) THEN
              DO IW=1,NW+1
                RATECM(IKK,IGR,IW)=RATECM(IKK,IGR,IW)*FACT(IKK,IW)
              ENDDO
              RATECM(IKK,IGR,NW+2)=RATECM(IKK,IGR,NW+2)*FACT(IKK,1)
              RATECM(IKK,IGR,NW+4)=RATECM(IKK,IGR,NW+4)*FACT(IKK,1)
              IF(NENER.GT.0) OVERV(IKK,IGR)=OVERV(IKK,IGR)*FACT(IKK,1)
              IF(LH) HFACT(IKK,IGR)=HFACT(IKK,IGR)*FACT(IKK,1)
              IF(ITRANC.NE.0) RATECM(IKK,IGR,NW+9)=RATECM(IKK,IGR,NW+9)
     >        *FACT(IKK,1)
              DO 10 IL=1,NL
                IW=MIN(IL,NW+1,2)
                SIGS(IKK,IGR,IL)=SIGS(IKK,IGR,IL)*FACT(IKK,IW)
  10          CONTINUE
            ELSE IF(IADJ.EQ.1) THEN
              DO IL=1,NW+1
                FAD1=FADJCM(IKK,IGR,IL)
                RATECM(IKK,IGR,IL)=RATECM(IKK,IGR,IL)*FACT(IKK,IL)/FAD1
              ENDDO
              FAD1=FADJCM(IKK,IGR,1)
              RATECM(IKK,IGR,NW+2)=RATECM(IKK,IGR,NW+2)*FACT(IKK,1)/FAD1
              RATECM(IKK,IGR,NW+4)=RATECM(IKK,IGR,NW+4)*FACT(IKK,1)/FAD1
              IF(NENER.GT.0) OVERV(IKK,IGR)=OVERV(IKK,IGR)*FACT(IKK,1)
     >        /FAD1
              IF(LH) HFACT(IKK,IGR)=HFACT(IKK,IGR)*FACT(IKK,1)/FAD1
              IF(ITRANC.NE.0) RATECM(IKK,IGR,NW+9)=RATECM(IKK,IGR,NW+9)
     >        *FACT(IKK,1)/FAD1
              DO 20 IL=1,NL
                IW=MIN(IL,NW+1,2)
                SIGS(IKK,IGR,IL)=SIGS(IKK,IGR,IL)*FACT(IKK,IW)/
     >          FADJCM(IKK,IGR,IW)
  20          CONTINUE
            ENDIF
            DO 30 IDEL=1,NDEL
              K=NW+12+IDEL
              RATECM(IKK,IGR,K)=RATECM(IKK,IGR,K)*FACT(IKK,1)
  30        CONTINUE
  40      CONTINUE
          IF((ILEAKS.EQ.1).OR.(ILEAKS.EQ.2)) THEN
            ZNU=0.0
            DEN=0.0
            IF(IADJ.EQ.0) THEN
              DO 50 IKK=1,NMERGE
                ZNU=ZNU+RATECM(IKK,IGR,NW+5)
                DEN=DEN+FLUXCM(IKK,IGR,1)
                RATECM(IKK,IGR,NW+5)=RATECM(IKK,IGR,NW+5)*FACT(IKK,1)
  50          CONTINUE
              DIFF(IGR)=ZNU/DEN
            ELSE IF(IADJ.EQ.1) THEN
              DEN2=0.0
              DO 60 IKK=1,NMERGE
                ZNU=ZNU+RATECM(IKK,IGR,NW+5)
                DEN=DEN+FLUXCM(IKK,IGR,1)
                DEN2=DEN2+FADJCM(IKK,IGR,1)
                RATECM(IKK,IGR,NW+5)=RATECM(IKK,IGR,NW+5)*FACT(IKK,1)/
     >          FADJCM(IKK,IGR,1)
  60          CONTINUE
              DIFF(IGR)=ZNU/(DEN*DEN2)
            ENDIF
          ELSE IF(ILEAKS.GT.0) THEN
            DO 70 IKK=1,NMERGE
              RATECM(IKK,IGR,NW+5)=RATECM(IKK,IGR,NW+5)*FACT(IKK,1)
              RATECM(IKK,IGR,NW+10)=RATECM(IKK,IGR,NW+10)*FACT(IKK,1)
              RATECM(IKK,IGR,NW+11)=RATECM(IKK,IGR,NW+11)*FACT(IKK,1)
              RATECM(IKK,IGR,NW+12)=RATECM(IKK,IGR,NW+12)*FACT(IKK,1)
  70        CONTINUE
          ENDIF
          DO 100 JGR=1,NGCOND
            DO 90 IKK=1,NMERGE
              DO 80 IL=1,NL
                IW=MIN(IL,NW+1)
                IF(IADJ.EQ.0) THEN
                  SCATTS(IKK,JGR,IGR,IL)=REAL(SCATTD(IKK,JGR,IGR,IL)
     >            *FACT(IKK,IW))
                ELSE IF(IADJ.EQ.1) THEN
                  SCATTS(IKK,JGR,IGR,IL)=REAL(SCATTD(IKK,JGR,IGR,IL)
     >            *FACT(IKK,IW)/FADJCM(IKK,JGR,IW))
                ENDIF
  80          CONTINUE
  90        CONTINUE
 100      CONTINUE
          DO 110 IKK=1,NMERGE
            RATECM(IKK,IGR,NW+7)=SCATTS(IKK,IGR,IGR,1)
 110      CONTINUE
          DO 130 IED=1,NEDMAC
            DO 120 IKK=1,NMERGE
              IF(IADJ.EQ.0) THEN
                TAUXE(IKK,IGR,IED)=TAUXE(IKK,IGR,IED)*FACT(IKK,1)
              ELSE IF(IADJ.EQ.1) THEN
                TAUXE(IKK,IGR,IED)=TAUXE(IKK,IGR,IED)*FACT(IKK,1)/
     >          FADJCM(IKK,IGR,1)
              ENDIF
 120        CONTINUE
 130      CONTINUE
 200    CONTINUE
        IF(NSAVES.EQ.2) THEN
*----
*  SAVE MERGED/CONDENSED X-S ON LCM
*----
          CALL LCMSIX(IPEDIT,CURNAM,ILCMUP)
          CALL LCMSIX(IPEDIT,'MACROLIB',ILCMUP)
          CALL LCMPUT(IPEDIT,'TIMESTAMP',3,2,TIMEF)
          IDATA(1)=NGCOND
          IDATA(2)=NMERGE
          IDATA(3)=NL
          IDATA(5)=NEDMAC
          IDATA(6)=ITRANC
          IDATA(7)=NDEL
          IDATA(15)=IADJ
          IF(NEDMAC.GT.0) THEN
            CALL LCMPTC(IPEDIT,'ADDXSNAME-P0',8,NEDMAC,HVECT)
          ENDIF
          JPEDIT=LCMLID(IPEDIT,'GROUP',NGCOND)
          DO 210 IGR=1,NGCOND
            KPEDIT=LCMDIL(JPEDIT,IGR)
            IF(NEDMAC.GT.0) THEN
              DO 211 IED=1,NEDMAC
                CEDNAM=HVECT(IED)
                IF(CEDNAM(:2).EQ.'NW') GO TO 211
                CALL LCMPUT(KPEDIT,CEDNAM,NMERGE,2,TAUXE(1,IGR,IED))
 211          CONTINUE
            ENDIF
            IF(NENER.GT.0) CALL LCMPUT(KPEDIT,'OVERV',NMERGE,2,
     >      OVERV(1,IGR))
            IF(LH) CALL LCMPUT(KPEDIT,'H-FACTOR',NMERGE,2,HFACT(1,IGR))
            DO IW=1,MIN(NW+1,10)
              WRITE(CEDNAM,'(4HNTOT,I1)') IW-1
              CALL LCMPUT(KPEDIT,CEDNAM,NMERGE,2,RATECM(1,IGR,IW))
            ENDDO
            CALL LCMPUT(KPEDIT,'ABS',NMERGE,2,RATECM(1,IGR,NW+2))
            CALL LCMPUT(KPEDIT,'PRODUCTION',NMERGE,2,RATECM(1,IGR,NW+4))
            DO 212 IKK=1,NMERGE
             RATECM(IKK,IGR,NW+6)=RATECM(IKK,IGR,1)-RATECM(IKK,IGR,NW+2)
 212        CONTINUE
            IF(IDATA(4).EQ.1) THEN
              CALL LCMPUT(KPEDIT,'NUSIGF',NMERGE,2,RATECM(1,IGR,NW+3))
              CALL LCMPUT(KPEDIT,'CHI',NMERGE,2,RATECM(1,IGR,NW+8))
              DO 901 IDEL=1,NDEL
                K=NW+12+IDEL
                WRITE(CEDNAM,'(6HNUSIGF,I2.2)') IDEL
                CALL LCMPUT(KPEDIT,CEDNAM,NMERGE,2,RATECM(1,IGR,K))
                WRITE(CEDNAM,'(3HCHI,I2.2)') IDEL
                CALL LCMPUT(KPEDIT,CEDNAM,NMERGE,2,RATECM(1,IGR,NDEL+K))
 901          CONTINUE
            ENDIF
            IF(ITRANC.NE.0) THEN
              CALL LCMPUT(KPEDIT,'TRANC',NMERGE,2,RATECM(1,IGR,NW+9))
            ENDIF
            IF((ILEAKS.EQ.1).OR.(ILEAKS.EQ.2).OR.(ILEAKS.EQ.10)) THEN
              CALL LCMPUT(KPEDIT,'DIFF',NMERGE,2,RATECM(1,IGR,NW+5))
            ELSE IF(ILEAKS.EQ.3) THEN
              CALL LCMPUT(KPEDIT,'DIFF',NMERGE,2,RATECM(1,IGR,NW+5))
              CALL LCMPUT(KPEDIT,'DIFFX',NMERGE,2,RATECM(1,IGR,NW+10))
              CALL LCMPUT(KPEDIT,'DIFFY',NMERGE,2,RATECM(1,IGR,NW+11))
              CALL LCMPUT(KPEDIT,'DIFFZ',NMERGE,2,RATECM(1,IGR,NW+12))
            ELSE IF(ILEAKS.EQ.11) THEN
              CALL LCMPUT(KPEDIT,'DIFFX',NMERGE,2,RATECM(1,IGR,NW+10))
              CALL LCMPUT(KPEDIT,'DIFFY',NMERGE,2,RATECM(1,IGR,NW+11))
              CALL LCMPUT(KPEDIT,'DIFFZ',NMERGE,2,RATECM(1,IGR,NW+12))
            ENDIF
            CALL LCMPUT(KPEDIT,'FLUX-INTG',NMERGE,2,FLUXCM(1,IGR,1))
            DO IL=2,MIN(NW+1,10)
              WRITE(CEDNAM,'(11HFLUX-INTG-P,I1)') IL-1
              CALL LCMPUT(KPEDIT,CEDNAM,NMERGE,2,FLUXCM(1,IGR,IL))
            ENDDO
            IF(IADJ.EQ.1) THEN
              DO IL=1,MIN(NW+1,10)
                WRITE(CEDNAM,'(4HNWAT,I1)') IL-1
                CALL LCMPUT(KPEDIT,CEDNAM,NMERGE,2,FADJCM(1,IGR,IL))
              ENDDO
            ENDIF
            DO 350 IL=1,NL
            WRITE (CM,'(I2.2)') IL-1
            IPOSIT=0
            DO 214 IKK=1,NMERGE
              J2=IGR
              J1=IGR
              DO 215 JGR=1,NGCOND
                IF(SCATTS(IKK,IGR,JGR,IL).NE.0.0) THEN
                  J2=MAX(J2,JGR)
                  J1=MIN(J1,JGR)
                ENDIF
 215          CONTINUE
              NJJ(IKK)=J2-J1+1
              IJJ(IKK)=J2
              IPOS(IKK)=IPOSIT+1
              DO 216 JGR=J2,J1,-1
                IPOSIT=IPOSIT+1
                SCATC(IPOSIT)=SCATTS(IKK,IGR,JGR,IL)
 216          CONTINUE
 214        CONTINUE
            CALL LCMPUT(KPEDIT,'SIGS'//CM,NMERGE,2,SIGS(1,IGR,IL))
            CALL LCMPUT(KPEDIT,'SIGW'//CM,NMERGE,2,SCATTS(1,IGR,IGR,IL))
            CALL LCMPUT(KPEDIT,'SCAT'//CM,IPOSIT,2,SCATC)
            CALL LCMPUT(KPEDIT,'NJJS'//CM,NMERGE,1,NJJ)
            CALL LCMPUT(KPEDIT,'IJJS'//CM,NMERGE,1,IJJ)
            CALL LCMPUT(KPEDIT,'IPOS'//CM,NMERGE,1,IPOS)
 350        CONTINUE
            IF(IPRINT.GE.4) THEN
              WRITE(IUNOUT,'(/14H G R O U P   :,I4)') IGR
              CALL LCMLIB(KPEDIT)
            ENDIF
 210      CONTINUE
          IF(ILEAKS.EQ.1) THEN
            CALL LCMPUT(IPEDIT,'DIFFB1HOM',NGCOND,2,DIFF)
            CALL LCMPUT(IPEDIT,'B2  B1HOM',1,2,B2(4))
          ELSE IF(ILEAKS.EQ.2) THEN
            CALL LCMPUT(IPEDIT,'B2  B1HOM',1,2,B2(4))
          ELSE IF(ILEAKS.EQ.3) THEN
            CALL LCMPUT(IPEDIT,'B2  B1HOM',1,2,B2(4))
            CALL LCMPUT(IPEDIT,'B2  HETE',3,2,B2)
          ENDIF
          IDATA(8)=NALBP
          DO 217 I=9,NSTATE
          IDATA(I)=0
 217      CONTINUE
          IF((ILEAKS.EQ.1).OR.(ILEAKS.EQ.2).OR.(ILEAKS.EQ.10)) THEN
             IDATA(9)=1
          ELSE IF((ILEAKS.EQ.3).OR.(ILEAKS.EQ.11)) THEN
             IDATA(9)=2
          ENDIF
          IDATA(10)=NW
          CALL LCMPUT(IPEDIT,'STATE-VECTOR',NSTATE,1,IDATA)
          HSIGN='L_MACROLIB'
          CALL LCMPTC(IPEDIT,'SIGNATURE',12,1,HSIGN)
          IF(NENER.GT.0) THEN
            CALL LCMPUT(IPEDIT,'ENERGY',NGCOND+1,2,WENERG)
            CALL LCMPUT(IPEDIT,'DELTAU',NGCOND,2,WLETYC)
          ENDIF
          CALL LCMPUT(IPEDIT,'VOLUME',NMERGE,2,VOLMER)
          IF((EIGENK.NE.0.0).AND.(NIFISS.GT.0)) THEN
            CALL LCMPUT(IPEDIT,'K-EFFECTIVE',1,2,EIGENK)
            CALL LCMPUT(IPEDIT,'K-INFINITY',1,2,CUREIN)
          ENDIF
          CALL LCMPUT(IPEDIT,'FLUXDISAFACT',NGCOND,2,DISFCT)
          IF(NALBP.GT.0) THEN
            LAL1D=.TRUE.
            DO IAL=1,NALBP
              DO IGR=1,NGCOND
                DO JGR=1,NGCOND
                  IF((IGR.NE.JGR).AND.(ALBP(IAL,IGR,JGR).NE.0.0)) THEN
                    LAL1D=.FALSE.
                    GO TO 218
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
 218        IF(LAL1D) THEN
*             diagonal physical albedos
              ALLOCATE(ALB1(NALBP,NGCOND))
              DO IAL=1,NALBP
                DO IGR=1,NGCOND
                  ALB1(IAL,IGR)=ALBP(IAL,IGR,IGR)
                ENDDO
              ENDDO
              CALL LCMPUT(IPEDIT,'ALBEDO',NALBP*NGCOND,2,ALB1)
              DEALLOCATE(ALB1)
            ELSE
*             matrix physical albedos
              CALL LCMPUT(IPEDIT,'ALBEDO',NALBP*NGCOND*NGCOND,2,ALBP)
            ENDIF
          ENDIF
          CALL LCMSIX(IPEDIT,' ',ILCMDN)
          CALL LCMSIX(IPEDIT,' ',ILCMDN)
          IF(IPRINT.GT.0) WRITE(IUNOUT,6031) CURNAM
        ENDIF
      ENDIF
*----
*  PRINT X-S
*----
      IF(IPRINT.GE.3) THEN
        WRITE(IUNOUT,6000)
        DO 170 IGR=1,NGCOND
          IF((ILEAKS.EQ.1).OR.(ILEAKS.EQ.2).OR.(ILEAKS.EQ.10)) THEN
            WRITE(IUNOUT,6020) IGR
          ELSE
            WRITE(IUNOUT,6021) IGR
          ENDIF
          DO 171 IKK=1,NMERGE
*----
*  UNCOMMENT THE 4 LINES TO PERFORM TRANSPORT CORRECTION
*----
            TOTAL=RATECM(IKK,IGR,1)
            SCATWG=SCATTS(IKK,IGR,IGR,1)
*           IF(ITRANC.NE.0) THEN
*             TOTAL=TOTAL-RATECM(IKK,IGR,NW+9)
*             SCATWG=SCATWG-RATECM(IKK,IGR,NW+9)
*           ENDIF
*
            IF (FLUXCM(IKK,IGR,1).NE.0.0) THEN
              FLXAVG=FLUXCM(IKK,IGR,1)/VOLMER(IKK)
              SCATTN=0.0D0
              DO 172 JGR=1,NGCOND
                 IF(JGR.NE.IGR) SCATTN=SCATTN+SCATTS(IKK,JGR,IGR,1)
 172          CONTINUE
              IF((ILEAKS.EQ.1).OR.(ILEAKS.EQ.2).OR.(ILEAKS.EQ.10)) THEN
                WRITE(IUNOUT,6022) IKK,FLXAVG,TOTAL,
     >          RATECM(IKK,IGR,NW+5),RATECM(IKK,IGR,NW+2),
     >          RATECM(IKK,IGR,NW+3),RATECM(IKK,IGR,NW+8),SCATWG,SCATTN
              ELSE
                WRITE(IUNOUT,6022) IKK,FLXAVG,TOTAL,
     >          RATECM(IKK,IGR,NW+2),RATECM(IKK,IGR,NW+3),
     >          RATECM(IKK,IGR,NW+8),SCATWG,SCATTN
              ENDIF
            ENDIF
 171      CONTINUE
          IF((ILEAKS.EQ.3).OR.(ILEAKS.EQ.11)) THEN
            WRITE(IUNOUT,6024)
            DO 173 IKK=1,NMERGE
              WRITE(IUNOUT,6025) IKK,RATECM(IKK,IGR,NW+10),
     >        RATECM(IKK,IGR,NW+11),RATECM(IKK,IGR,NW+12),
     >        RATECM(IKK,IGR,NW+5)
 173        CONTINUE
          ENDIF
          WRITE(IUNOUT,6026) DISFCT(IGR)
 170    CONTINUE
      ENDIF
      IF(IPRINT.GE.4) THEN
        DO 190 IKK=1,NMERGE
          WRITE(IUNOUT,6027) IKK,(JGR,JGR=1,NGCOND)
          DO 180 IGR=1,NGCOND
*----
*  UNCOMMENT THE FOLLOWING LINE TO PERFORM TRANSPORT CORRECTION
*----
            SCATWG=SCATTS(IKK,IGR,IGR,1)
*           IF(ITRANC.NE.0) SCATWG=SCATWG-RATECM(IKK,IGR,NW+9)
*
            WRITE(IUNOUT,6028) IGR,(SCATTS(IKK,JGR,IGR,1),JGR=1,IGR-1),
     >      SCATWG,(SCATTS(IKK,JGR,IGR,1),JGR=IGR+1,NGCOND)
 180      CONTINUE
          WRITE (IUNOUT,'(//)')
 190    CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(DIFF,FACT,SCATC)
      DEALLOCATE(IPOS,NJJ,IJJ)
      RETURN
*----
*  FORMAT
*----
 6000 FORMAT(/' F L U X E S   A N D   H O M O G E N I Z E D   X - S'/
     > 1X,51(1H-))
 6020 FORMAT(/' G R O U P   :',I4/
     >1X,'REGION',3X,'AVERAGE',9X,'NTOT0',7X,'DIFFUSION',5X,
     >'ABSORPTION',5X,'NUSIGF',8X,'FISSION',10X,'SCATTERING X-S'/11X,
     >'FLUX',12X,'X-S',7X,'COEFFICIENT',7X,'X-S',10X,'X-S',10X,
     >'SPECTRUM',2X,'WITHIN GROUP',2X,'OUT OF GROUP')
 6021 FORMAT(/' G R O U P   :',I4/
     >1X,'REGION',3X,'AVERAGE',9X,'NTOT0',7X,
     >'ABSORPTION',5X,'NUSIGF',8X,'FISSION',10X,'SCATTERING X-S'/11X,
     >'FLUX',12X,'X-S',11X,'X-S',10X,'X-S',10X,'SPECTRUM',2X,
     >'WITHIN GROUP',2X,'OUT OF GROUP')
 6022 FORMAT(1X,I4,1P,8E14.5)
 6024 FORMAT(/' REGION     X-LEAKAGE     Y-LEAKAGE     Z-LEAKAGE',
     >'    HOM-LEAKAGE'/'           COEFFICIENT   COEFFICIENT   ',
     >'COEFFICIENT   COEFFICIENT')
 6025 FORMAT(1X,I6,1X,1P,5E14.5)
 6026 FORMAT(/' FLUX DISADVANTAGE FACTOR =',1P,E14.5)
 6027 FORMAT(/47H SCATTERING TRANSFER X-S (I TOWARD J) IN REGION,I5,1H:
     > //(11X,2HJ=,I4,:,6X,2HJ=,I4,:,6X,2HJ=,I4,:,6X,2HJ=,I4,:,6X,2HJ=,
     > I4,:,6X,2HJ=,I4,:,6X,2HJ=,I4,:,6X,2HJ=,I4,:,6X,2HJ=,I4,:,6X,
     > 2HJ=,I4))
 6028 FORMAT(3H I=,I4,2H: ,1P,10E12.4/(9X,10E12.4))
 6031 FORMAT(/53H MERGED/CONDENSED SET OF X-S SAVED IN LCM DIRECTORY ',
     > A12,2H'./)
      END
