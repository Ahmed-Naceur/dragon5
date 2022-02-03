*DECK EDIBAL
      SUBROUTINE EDIBAL(IPEDIT,IPFLUX,IPRINT,NL,IFFAC,NGCOND,NMERGE,
     >                  EIGENK,RATECM,FLUXCM,SCATTS,ILEAKS,B2,NW,
     >                  NTAUXT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Four factor calculation.
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
* IPFLUX  pointer to the flux LCM object.
* IPRINT  print level:
*         =1 neutron balance or four factor;
*         =2 material analysis (if available).
* NL      number of Legendre orders.
* IFFAC   number of neutrons for neutron balance.
* NGCOND  number of condensed groups.
* NMERGE  number of merge regions.
* EIGENK  problem eigenvalue.
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
*         = RATECM(*,NW+12) = z-directed leakage.
* FLUXCM  integrated region/group fluxes:
*         = FLUXCM(*,1) = fluxes P0;
*         = FLUXCM(*,2) = fluxes P1.
* SCATTS  scattering matrix.
* ILEAKS  leakage calculation flag:
*         = 0 no leakage;
*         = 1 homogeneous leakage (Diffon);
*         = 2 isotropic streaming (Ecco);
*         = 3 anisotropic streaming (Tibere).
* B2      square buckling:
*         for ILEAKS=1,2: B2(4) is homogeneous;
*         for ILEAKS=3: B2(1),B2(2),B2(3) are directional heterogeneous
*         and B2(4) is homogeneous.
* NW      type of weighting for PN cross section info (=0 P0; =1 P1).
* NTAUXT  number of reaction rate edits (=12+NW+2*NDEL).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR)    IPEDIT,IPFLUX
      INTEGER        IPRINT,NL,IFFAC,NGCOND,NMERGE,ILEAKS,NW,NTAUXT
      REAL           EIGENK,RATECM(NMERGE,NGCOND,NTAUXT),
     >               FLUXCM(NMERGE,NGCOND,NW+1),
     >               SCATTS(NMERGE,NGCOND,NGCOND,NL),B2(4)
*----
*  LOCAL VARIABLES
*----
      SAVE           CNAMAT
      PARAMETER     (IUNOUT=6,INAMAT=2)
      CHARACTER      CNAMAT(INAMAT)*15
      REAL           XN(8)
      DOUBLE PRECISION DACCK,BUCKL2,XNF(3),XKINF,XN1,XN2,XN3,XN4,XN5,
     > XN6,XKEFF,XA,XLAMF,XLAMTH,XNORMF,XAUX
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ITYPER
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FLXINT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: DLEAK
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DRATE
      DATA          (CNAMAT(JJ),JJ=1,INAMAT)
     >           /'  FUEL         ','NON-FUEL       '/
*----
*  SCRATCH STORAGE ALLOCATION
*   ITYPER  region type.
*   FLXINT  integrated flux.
*   DRATE   group reaction rates:
*           DRATE(1,*,*)  for fuel;
*           DRATE(2,*,*)  for non-fuel;
*           DRATE(*,1,*)  production;
*           DRATE(*,2,*)  absorption;
*           DRATE(*,*,ngcond+1)  group total.
*   DLEAK   group leak rates:
*           DLEAK(1,*)  leakage;
*           DLEAK(2,*)  leakage + absorption;
*           DLEAK(*,ngcond+1)  group total.
*----
      ALLOCATE(ITYPER(NMERGE))
      ALLOCATE(FLXINT(NMERGE,NGCOND))
      ALLOCATE(DRATE(2,2,NGCOND+1),DLEAK(2,NGCOND+1))
*----
*  LOCALIZE FUEL AND NON-FUEL REGION
*----
      NGC1=NGCOND+1
      DO 100 IGR=1,NGC1
        DO 110 IRAT=1,2
          DRATE(1,IRAT,IGR)=0.0D0
          DRATE(2,IRAT,IGR)=0.0D0
 110    CONTINUE
        DLEAK(1,IGR)=0.0D0
        DLEAK(2,IGR)=0.0D0
 100  CONTINUE
      DO 120 IREG=1,NMERGE
        ITYPER(IREG)=2
 120  CONTINUE
      DO 130 IGR=1,NGCOND
        DO 140 IREG=1,NMERGE
          IF(RATECM(IREG,IGR,NW+3).GT.0.0) ITYPER(IREG)=1
 140    CONTINUE
 130  CONTINUE
*----
*  GET BUCKL2
*----
      CALL LCMLEN(IPFLUX,'B2  B1HOM',ILCMLN,ILCMTY)
      IF(ILCMLN.EQ.1) THEN
        CALL LCMGET(IPFLUX,'B2  B1HOM',BL2)
        BUCKL2=DBLE(BL2)
      ELSE
        BUCKL2=0.0D0
      ENDIF
      IF(EIGENK.EQ.0.0) THEN
        WRITE(IUNOUT,7000)
        FLXRGE=1.0
        FLXREN=1.0
        DACCK=1.0D0
        WRITE(IUNOUT,6000)
      ELSE
        FLXRGE=REAL(IFFAC)
        FLXREN=REAL(IFFAC)*EIGENK
        DACCK=DBLE(EIGENK)
        IF(NGCOND.NE.3) THEN
          WRITE(IUNOUT,7001) NGCOND
          WRITE(IUNOUT,6001)
        ELSE
          WRITE(IUNOUT,6002)
        ENDIF
      ENDIF
      IF(IPRINT.EQ.1) THEN
        WRITE(IUNOUT,6100)
        WRITE(IUNOUT,6101)(IREG,CNAMAT(ITYPER(IREG)),IREG=1,NMERGE)
      ENDIF
*----
*  FIND INTEGRATED FLUX DIVIDED BY SPH FACTOR
*----
      DO 150 IGR=1,NGCOND
        DO 160 IREG=1,NMERGE
          FLXINT(IREG,IGR)=FLUXCM(IREG,IGR,1)
 160    CONTINUE
 150  CONTINUE
      IF(IPRINT.GE.1) THEN
        WRITE(IUNOUT,6200)
      ENDIF
      DO 170 IGR=1,NGCOND
*----
*  REACTION RATES PER GROUP AND MATERIAL TYPE
*----
        IF(IPRINT.GE.1) THEN
          WRITE(IUNOUT,6201) IGR,NGCOND
        ENDIF
        DO 180 IREG=1,NMERGE
*----
*  ADD PRODUCTION AND ABSORPTION IN MATERIAL TYPE
*----
          IF(ITYPER(IREG).EQ.1) THEN
            DRATE(1,1,IGR)=DRATE(1,1,IGR)+DBLE(RATECM(IREG,IGR,NW+4))
     >                   *DBLE(FLXINT(IREG,IGR))
            DRATE(1,2,IGR)=DRATE(1,2,IGR)+DBLE(RATECM(IREG,IGR,NW+2))
     >                   *DBLE(FLXINT(IREG,IGR))
          ELSE
            DRATE(2,1,IGR)=DRATE(2,1,IGR)+DBLE(RATECM(IREG,IGR,NW+4))
     >                   *DBLE(FLXINT(IREG,IGR))
            DRATE(2,2,IGR)=DRATE(2,2,IGR)+DBLE(RATECM(IREG,IGR,NW+2))
     >                   *DBLE(FLXINT(IREG,IGR))
          ENDIF
*----
*  PRINT PRODUCTION AND ABSORPTION PER REGION
*----
          IF(IPRINT.GE.2) THEN
            IF(IFFAC.EQ.1000) THEN
              WRITE(IUNOUT,6300) IREG,CNAMAT(ITYPER(IREG)),
     >                  FLXRGE*RATECM(IREG,IGR,NW+4)*FLXINT(IREG,IGR),
     >                  FLXREN*RATECM(IREG,IGR,NW+2)*FLXINT(IREG,IGR)
            ELSE
              WRITE(IUNOUT,6301) IREG,CNAMAT(ITYPER(IREG)),
     >                  FLXRGE*RATECM(IREG,IGR,NW+4)*FLXINT(IREG,IGR),
     >                  FLXREN*RATECM(IREG,IGR,NW+2)*FLXINT(IREG,IGR)
            ENDIF
          ENDIF
*----
*  ADD GROUP LEAKAGE
*----
          IF((ILEAKS.EQ.1).OR.(ILEAKS.EQ.2)) THEN
            DLEAK(1,IGR)=DLEAK(1,IGR)+DBLE(RATECM(IREG,IGR,NW+5))
     >                        *B2(4)*DBLE(FLXINT(IREG,IGR))
          ELSE IF(ILEAKS.EQ.3) THEN
            DLEAK(1,IGR)=DLEAK(1,IGR)+DBLE(FLXINT(IREG,IGR))
     >       *(DBLE(RATECM(IREG,IGR,NW+10))*B2(1)
     >        +DBLE(RATECM(IREG,IGR,NW+11))*B2(2)
     >        +DBLE(RATECM(IREG,IGR,NW+12))*B2(3))
          ELSE
            DO 190 JGR=1,NGCOND
              IF(JGR.NE.IGR)
     >          DLEAK(2,IGR)=DLEAK(2,IGR)+DBLE(SCATTS(IREG,IGR,JGR,1))
     >               *DBLE(FLXINT(IREG,JGR))
 190        CONTINUE
            DLEAK(2,IGR)=DLEAK(2,IGR)+DBLE(FLXINT(IREG,IGR))*
     >                     ( DBLE(RATECM(IREG,IGR,NW+4))/DACCK
     >                      -DBLE(RATECM(IREG,IGR,1))
     >                      +DBLE(RATECM(IREG,IGR,NW+2))
     >                      +DBLE(SCATTS(IREG,IGR,IGR,1)) )
          ENDIF
 180    CONTINUE
        IF((ILEAKS.EQ.1).OR.(ILEAKS.EQ.2)) THEN
          DLEAK(2,IGR)=DLEAK(1,IGR)+DRATE(1,2,IGR)+DRATE(2,2,IGR)
        ELSE
          DLEAK(1,IGR)=DLEAK(2,IGR)-DRATE(1,2,IGR)-DRATE(2,2,IGR)
        ENDIF
*----
*  PRINT PRODUCTION AND ABSORPTION PER MATERIAL TYPE AND TOTAL
*  PRODUCTION, ABSORPTION, LEAKAGE AND LEAKAGE+ABSORPTION
*----
        IF(IPRINT.GE.1) THEN
          IF(IFFAC.EQ.1000) THEN
            WRITE(IUNOUT,6302) 'TOTAL FUEL     ',
     >                  FLXRGE*DRATE(1,1,IGR),
     >                  FLXREN*DRATE(1,2,IGR)
            WRITE(IUNOUT,6302) 'TOTAL NON-FUEL ',
     >                  FLXRGE*DRATE(2,1,IGR),
     >                  FLXREN*DRATE(2,2,IGR)
            WRITE(IUNOUT,6302) 'FUEL + NON-FUEL',
     >                  FLXRGE*(DRATE(1,1,IGR)+DRATE(2,1,IGR)),
     >                  FLXREN*(DRATE(1,2,IGR)+DRATE(2,2,IGR)),
     >                  FLXREN*DLEAK(1,IGR),FLXREN*DLEAK(2,IGR)
          ELSE
            WRITE(IUNOUT,6303) 'TOTAL FUEL     ',
     >                  FLXRGE*DRATE(1,1,IGR),
     >                  FLXREN*DRATE(1,2,IGR)
            WRITE(IUNOUT,6303) 'TOTAL NON-FUEL ',
     >                  FLXRGE*DRATE(2,1,IGR),
     >                  FLXREN*DRATE(2,2,IGR)
            WRITE(IUNOUT,6303) 'FUEL + NON-FUEL',
     >                  FLXRGE*(DRATE(1,1,IGR)+DRATE(2,1,IGR)),
     >                  FLXREN*(DRATE(1,2,IGR)+DRATE(2,2,IGR)),
     >                  FLXREN*DLEAK(1,IGR),FLXREN*DLEAK(2,IGR)
          ENDIF
        ENDIF
*----
*  GROUP SUM
*----
        DRATE(1,1,NGC1)=DRATE(1,1,NGC1)+DRATE(1,1,IGR)
        DRATE(2,1,NGC1)=DRATE(2,1,NGC1)+DRATE(2,1,IGR)
        DRATE(1,2,NGC1)=DRATE(1,2,NGC1)+DRATE(1,2,IGR)
        DRATE(2,2,NGC1)=DRATE(2,2,NGC1)+DRATE(2,2,IGR)
        DLEAK(1,NGC1)=DLEAK(1,NGC1)+DLEAK(1,IGR)
        DLEAK(2,NGC1)=DLEAK(2,NGC1)+DLEAK(2,IGR)
 170  CONTINUE
      IF(NGCOND.GT.1) THEN
        WRITE(IUNOUT,6202)
        IF(IFFAC.EQ.1000) THEN
          WRITE(IUNOUT,6302) 'TOTAL FUEL     ',
     >                FLXRGE*DRATE(1,1,NGC1),
     >                FLXREN*DRATE(1,2,NGC1)
          WRITE(IUNOUT,6302) 'TOTAL NON-FUEL ',
     >                FLXRGE*DRATE(2,1,NGC1),
     >                FLXREN*DRATE(2,2,NGC1)
          WRITE(IUNOUT,6302) 'FUEL + NON-FUEL',
     >                FLXRGE*(DRATE(1,1,NGC1)+DRATE(2,1,NGC1)),
     >                FLXREN*(DRATE(1,2,NGC1)+DRATE(2,2,NGC1)),
     >                FLXREN*DLEAK(1,NGC1),FLXREN*DLEAK(2,NGC1)
        ELSE
          WRITE(IUNOUT,6303) 'TOTAL FUEL     ',
     >                FLXRGE*DRATE(1,1,NGC1),
     >                FLXREN*DRATE(1,2,NGC1)
          WRITE(IUNOUT,6303) 'TOTAL NON-FUEL ',
     >                FLXRGE*DRATE(2,1,NGC1),
     >                FLXREN*DRATE(2,2,NGC1)
          WRITE(IUNOUT,6303) 'FUEL + NON-FUEL',
     >                FLXRGE*(DRATE(1,1,NGC1)+DRATE(2,1,NGC1)),
     >                FLXREN*(DRATE(1,2,NGC1)+DRATE(2,2,NGC1)),
     >                FLXREN*DLEAK(1,NGC1),FLXREN*DLEAK(2,NGC1)
        ENDIF
      ENDIF
      IF( (EIGENK.GT.0.0) .AND. (NGCOND.EQ.3) ) THEN
*----
*  FOUR FACTOR CALCULATION
*----
        XNF(1)=0.0D0
        XNF(2)=0.0D0
        XNF(3)=0.0D0
        DO 200 IREG=1,NMERGE
*----
*  ADD NUSIGF AND TOTAL RATES IN FUEL
*----
          IF(ITYPER(IREG).EQ.1) THEN
            XNF(3)=XNF(3)+DBLE(RATECM(IREG,3,NW+3))*DBLE(FLXINT(IREG,3))
            XNF(2)=XNF(2)+DBLE(RATECM(IREG,2,NW+3))*DBLE(FLXINT(IREG,2))
            XNF(1)=XNF(1)+DBLE(RATECM(IREG,1,NW+3))*DBLE(FLXINT(IREG,1))
          ENDIF
 200    CONTINUE
        XKINF=DRATE(1,1,4)/(DRATE(1,2,4)+DRATE(2,2,4))
        XKEFF=DRATE(1,1,4)/DLEAK(2,4)
        XN1=(DRATE(1,1,4)/XKEFF)-DLEAK(1,4)
        XN2=XN1-(DRATE(1,2,1)+DRATE(2,2,1))
        XN3=XN2-(DRATE(1,2,2)+DRATE(2,2,2))
        XN4=XN3-DRATE(2,2,3)
        XN5=XNF(2)+XNF(3)
        XN6=XN5+XNF(1)
*----
*  COMPUTATION OF EPSILON, F, ETA, P, LAMBDAF, LAMBDATH, KINF, KEFF
*----
        XN(1)=REAL((XN6-DRATE(1,2,1)-DRATE(2,2,1))/XN5)
        XN(5)=REAL(XN4/XN3)
        XN(6)=REAL(XN5/XN4)
        XN(3)=REAL(XKINF)/(XN(1)*XN(5)*XN(6))
        XA=DRATE(1,2,4)+DRATE(2,2,4)
        XLAMF=1/(1+((DLEAK(1,1)+DLEAK(1,2))/XA))
        XLAMTH=1/(1+(DLEAK(1,3)/XA))
        XAUX=XLAMF*XLAMTH
        XNORMF=XKEFF/(XKINF*XAUX)
        XN(2)=REAL(XLAMF*SQRT(XNORMF))
        XN(4)=REAL(XLAMTH*SQRT(XNORMF))
        XN(7)=XN(1)*XN(3)*XN(5)*XN(6)
        XN(8)=XN(7)*XN(2)*XN(4)
        IF(IPRINT.GE.1) THEN
          WRITE(IUNOUT,6400) (XN(JJ),JJ=1,8)
        ENDIF
        CALL LCMPUT(IPEDIT,'FOUR-FACTOR ',8,2,XN)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(DLEAK,DRATE)
      DEALLOCATE(FLXINT)
      DEALLOCATE(ITYPER)
      RETURN
*----
*  PRINT FORMATS
*----
 6000 FORMAT(/20X,'F I X E D    S O U R C E',
     >'    N E U T R O N    B A L A N C E')
 6001 FORMAT(/20X,'S I M P L E    N E U T R O N    B A L A N C E')
 6002 FORMAT(/20X,'F O U R    F A C T O R    C A L C U L A T I O N')
 6100 FORMAT(4(4X,'REGION',5X,'MATERIAL TYPE  '))
 6101 FORMAT(4(5X,I5,5X,A15))
 6200 FORMAT(15X,' REGION',3X,'MATERIAL TYPE  ',3X,'NEUTRON PRODUCTION',
     >       3X,'NEUTRON ABSORPTION',3X,'NEUTRON LEAKAGE   ',
     >       3X,'ABSORPTION+LEAKAGE')
 6201 FORMAT(' GROUP :',I3,'/',I3)
 6202 FORMAT(' SUM OVER GROUPS')
 6300 FORMAT(17X,I5,3X,A15,3X,F12.1,8X,F12.1)
 6301 FORMAT(17X,I5,3X,A15,3X,1P,E15.7,5X,E15.7)
 6302 FORMAT(25X,A15,3X,F12.1,8X,F12.1,8X,F12.1,8X,F12.1)
 6303 FORMAT(25X,A15,3X,1P,E15.7,5X,E15.7,5X,E15.7,5X,E15.7)
 6400 FORMAT(/' FOUR FACTORS'/1P,
     >        ' EPSILON (FAST FISSION FACTOR)  =',E15.7/
     >        ' LAMBDAF (FAST NON-LEAKAGE)     =',E15.7/
     >        ' P (ANTITRAP FACTOR)            =',E15.7/
     >        ' LAMBDAT (THERMAL NON-LEAKAGE)  =',E15.7/
     >        ' F (THERMAL UTILIZATION FACTOR) =',E15.7/
     >        ' ETA (THERMAL REPRODUCTION)     =',E15.7/
     >        ' INFINITE MULTIPLICATION (4F)   =',E15.7/
     >        ' EFFECTIVE MULTIPLICATION       =',E15.7)
*----
*  WARNING FORMATS
*----
 7000 FORMAT(' * * *  W A R N I N G  * * * ',/
     >       ' NO FOUR FACTOR CALCULATION PERMITTED FOR FIXED SOURCE',
     >       ' PROBLEM',/' A SIMPLE GROUP BY GROUP NEUTRON',
     >       ' BALANCE PERFORMED HERE')
 7001 FORMAT(' * * *  W A R N I N G  * * * ',/
     >       ' FOUR FACTOR CALCULATION REQUIRES 3 GROUPS, NUMBER OF',
     >       ' GROUPS HERE IS =',I10/' A SIMPLE GROUP BY GROUP NEUTRON',
     >       ' BALANCE PERFORMED HERE')
      END
