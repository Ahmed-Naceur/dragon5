*DECK USSCOR
      SUBROUTINE USSCOR(MAXNOR,IGRP,IPSYS,IASM,IRES,NBNRS,NIRES,NOR,
     1 CONR,IPPT1,IPPT2,WEIGH,TOTPT,SIGX,VOLMER)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the collision probability matrix taking into account the
* correlation effects between pairs of resonant isotopes in the same
* energy group.
*
*Copyright:
* Copyright (C) 2003 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* MAXNOR  maximum order of the probability tables (PT).
* IGRP    energy group index.
* IPSYS   pointer to the internal microscopic cross section library.
*         builded by the self-shielding module.
* IASM    offset in IPSYS.
* IRES    index of the resonant isotope been processed.
* NBNRS   number of correlated fuel regions.
* NIRES   exact number of correlated resonant isotopes.
* NOR     exact order of the probability table.
* CONR    number density of the resonant isotopes.
* IPPT1   pointer to LCM directory of each resonant isotope.
* IPPT2   information related to each resonant isotope:
*         IPPT2(:,1) index of a resonant region (used with infinite
*         dilution case);
*         IPPT2(:,2:4) alias name of resonant isotope.
* WEIGH   multiband weights.
* TOTPT   base points in total xs.
* SIGX    macroscopic total xs of the non-resonant isotopes in each fuel
*         region.
* VOLMER  volumes of the resonant and non-resonant regions.
*
* Reference: 
*  A. Hebert, "A Mutual Resonance Self-Shielding Model Consistent with
*  Ribon Subgroup Equations", Int. Mtg. on the Physics of Fuel Cycles
*  and Advanced Nuclear Systems: Global Developments. PHYSOR-2004,
*  Chicago, Illinois, April 25 - 29, 2004.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPSYS,IPPT1(NIRES)
      INTEGER MAXNOR,IGRP,IASM,IRES,NBNRS,NIRES,NOR(NIRES),
     1 IPPT2(NIRES,4)
      REAL CONR(NBNRS,NIRES),WEIGH(MAXNOR,NIRES),TOTPT(MAXNOR,NIRES),
     1 SIGX(NBNRS,NIRES),VOLMER(0:NBNRS)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) KPSYS,KPLIB1
      CHARACTER TEXT12*12
      LOGICAL LMOD
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: WCOR,SIGR
      REAL, ALLOCATABLE, DIMENSION(:,:) :: PIJ2,PIJ3,DILW
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: DIL,PIJ4
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(WCOR(MAXNOR**2),SIGR(NBNRS),DIL(0:NBNRS,0:NBNRS,MAXNOR),
     1 PIJ2(0:NBNRS,0:NBNRS),PIJ3(0:NBNRS,0:NBNRS),
     2 PIJ4(0:NBNRS,0:NBNRS,MAXNOR),DILW(0:NBNRS,0:NBNRS))
*----
*  COMPUTE THE MULTIBAND DILUTION MATRICES.
*----
      LMOD=(VOLMER(0).EQ.0.0)
      NQT1=NOR(IRES)
      DO 35 I1=1,NQT1
      DO 20 IND=1,NBNRS
      SIGR(IND)=CONR(IND,IRES)*TOTPT(I1,IRES)
      DO 10 JRES=1,NIRES
      IF(JRES.NE.IRES) SIGR(IND)=SIGR(IND)+SIGX(IND,JRES)
   10 CONTINUE
   20 CONTINUE
      KPSYS=LCMGIL(IPSYS,IASM+I1)
      CALL LCMGET(KPSYS,'DRAGON-PAV',DIL(0,0,I1))
      IF(LMOD) THEN
         CALL ALINV(NBNRS,DIL(1,1,I1),NBNRS+1,IER)
      ELSE
         CALL ALINV(NBNRS+1,DIL(0,0,I1),NBNRS+1,IER)
      ENDIF
      IF(IER.NE.0) CALL XABORT('USSCOR: SINGULAR MATRIX(1).')
      DO 30 IND=1,NBNRS
      DIL(IND,IND,I1)=DIL(IND,IND,I1)-SIGR(IND)
   30 CONTINUE
   35 CONTINUE
*----
*  COMPUTE THE AVERAGED COLLISION PROBABILITY MATRICES.
*----
      IF(NIRES.EQ.2) THEN
         JRES=MOD(IRES,NIRES)+1
         CALL LCMLEL(IPPT1(IRES),IGRP,ILONG1,ITYLCM)
         CALL LCMLEL(IPPT1(JRES),IGRP,ILONG2,ITYLCM)
         IF((ILONG1.NE.0).AND.(ILONG2.NE.0)) THEN
            KPLIB1=LCMGIL(IPPT1(IRES),IGRP)
*
*           COMPUTE THE FULLY CORRELATED PIJ MATRIX.
            WRITE(TEXT12,'(3A4)') (IPPT2(JRES,I0),I0=2,4)
            CALL LCMGET(KPLIB1,TEXT12,WCOR)
            NQT2=NOR(JRES)
            DO 100 I1=1,NQT1
            DO 40 IND=1,NBNRS
            SIGR(IND)=CONR(IND,IRES)*TOTPT(I1,IRES)
   40       CONTINUE
            KPSYS=LCMGIL(IPSYS,IASM+I1)
            CALL LCMGET(KPSYS,'DRAGON-PAV',PIJ3(0,0))
            CALL USSSEK(NBNRS,NQT2,LMOD,SIGR,CONR(1,JRES),WEIGH(1,JRES),
     1      TOTPT(1,JRES),PIJ3(0,0),DIL(0,0,I1))
            CALL XDRSET(PIJ3(0,0),(NBNRS+1)**2,0.0)
            DO 95 I2=1,NQT2
            WWW=WCOR((I2-1)*NQT1+I1)/WEIGH(I1,IRES)
            DO 60 I=0,NBNRS
            DO 50 J=0,NBNRS
            PIJ2(I,J)=DIL(I,J,I1)
   50       CONTINUE
   60       CONTINUE
            DO 70 I=1,NBNRS
            PIJ2(I,I)=PIJ2(I,I)+SIGR(I)+CONR(I,JRES)*TOTPT(I2,JRES)
   70       CONTINUE
            IF(LMOD) THEN
               CALL ALINV(NBNRS,PIJ2(1,1),NBNRS+1,IER)
            ELSE
               CALL ALINV(NBNRS+1,PIJ2(0,0),NBNRS+1,IER)
            ENDIF
            IF(IER.NE.0) CALL XABORT('USSCOR: SINGULAR MATRIX(2).')
            DO 90 I=0,NBNRS
            DO 80 J=0,NBNRS
            PIJ3(I,J)=PIJ3(I,J)+WWW*PIJ2(I,J)
   80       CONTINUE
   90       CONTINUE
   95       CONTINUE
*
*           STORE CORRECTED PIJ MATRIX.
            CALL LCMPUT(KPSYS,'DRAGON-PAV',(NBNRS+1)**2,2,PIJ3(0,0))
  100       CONTINUE
         ENDIF
      ELSE IF(NIRES.GT.1) THEN
         DO 110 I1=1,NQT1
         KPSYS=LCMGIL(IPSYS,IASM+I1)
         CALL LCMGET(KPSYS,'DRAGON-PAV',PIJ4(0,0,I1))
  110    CONTINUE
         DO 200 JRES=1,NIRES
         CALL LCMLEL(IPPT1(IRES),IGRP,ILONG1,ITYLCM)
         CALL LCMLEL(IPPT1(JRES),IGRP,ILONG2,ITYLCM)
         IF((JRES.NE.IRES).AND.(ILONG1.NE.0).AND.(ILONG2.NE.0)) THEN
            KPLIB1=LCMGIL(IPPT1(IRES),IGRP)
*
*           COMPUTE THE FULLY CORRELATED PIJ MATRIX.
            WRITE(TEXT12,'(3A4)') (IPPT2(JRES,I0),I0=2,4)
            CALL LCMGET(KPLIB1,TEXT12,WCOR)
            NQT2=NOR(JRES)
            DO 190 I1=1,NQT1
            DO 130 I=0,NBNRS
            DO 120 J=0,NBNRS
            DILW(I,J)=DIL(I,J,I1)
  120       CONTINUE
  130       CONTINUE
            DO 145 IND=1,NBNRS
            SIGR(IND)=CONR(IND,IRES)*TOTPT(I1,IRES)
            DO 140 KRES=1,NIRES
            IF((KRES.NE.IRES).AND.(KRES.NE.JRES)) THEN
               SIGR(IND)=SIGR(IND)+SIGX(IND,KRES)
            ENDIF
  140       CONTINUE
  145       CONTINUE
            CALL USSSEK(NBNRS,NQT2,LMOD,SIGR,CONR(1,JRES),WEIGH(1,JRES),
     1      TOTPT(1,JRES),PIJ4(0,0,I1),DILW(0,0))
*
            CALL XDRSET(PIJ3,(NBNRS+1)**2,0.0)
            DO 172 I2=1,NQT2
            WWW=WCOR((I2-1)*NQT1+I1)/WEIGH(I1,IRES)
            DO 155 I=0,NBNRS
            DO 150 J=0,NBNRS
            PIJ2(I,J)=DILW(I,J)
  150       CONTINUE
  155       CONTINUE
            DO 160 I=1,NBNRS
            PIJ2(I,I)=PIJ2(I,I)+SIGR(I)+CONR(I,JRES)*TOTPT(I2,JRES)
  160       CONTINUE
            IF(LMOD) THEN
               CALL ALINV(NBNRS,PIJ2(1,1),NBNRS+1,IER)
            ELSE
               CALL ALINV(NBNRS+1,PIJ2(0,0),NBNRS+1,IER)
            ENDIF
            IF(IER.NE.0) CALL XABORT('USSCOR: SINGULAR MATRIX(3).')
            DO 171 I=0,NBNRS
            DO 170 J=0,NBNRS
            PIJ3(I,J)=PIJ3(I,J)+WWW*PIJ2(I,J)
  170       CONTINUE
  171       CONTINUE
  172       CONTINUE
            KPSYS=LCMGIL(IPSYS,IASM+I1)
            CALL LCMGET(KPSYS,'DRAGON-PAV',PIJ2(0,0))
            DO 185 I=0,NBNRS
            DO 180 J=0,NBNRS
            IF(PIJ4(I,J,I1).NE.0.0) THEN
               PIJ2(I,J)=PIJ2(I,J)*PIJ3(I,J)/PIJ4(I,J,I1)
            ENDIF
  180       CONTINUE
  185       CONTINUE
*
*           STORE CORRECTED PIJ MATRIX.
            CALL LCMPUT(KPSYS,'DRAGON-PAV',(NBNRS+1)**2,2,PIJ2(0,0))
  190       CONTINUE
         ENDIF
  200    CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(DILW,PIJ4,PIJ3,PIJ2,DIL,SIGR,WCOR)
      RETURN
      END
