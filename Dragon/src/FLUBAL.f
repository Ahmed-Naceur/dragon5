*DECK FLUBAL
      SUBROUTINE FLUBAL(IPMACR,NGRP,ILEAK,NMAT,NREG,ICREB,NUNKNO,
     1 NANIS,MATCOD,VOL,KEYFLX,XSTRC,XSDIA,XCSOU,IGDEB,B2,DDD,
     2 KEYCUR,MATALB,ALBEDO,SURFAC,FUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Flux rebalancing for non converged groups with up-scattering.
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
* IPMACR  pointer to the macrolib LCM object.
* NGRP    number of energy groups.
* ILEAK   method used to include DB2 effect:
*         <5 uniform DB2 model;
*         =5 Ecco-type isotropic streaming model;
*         >5 Tibere anisotropic streaming model.
* NMAT    number of mixtures.
* NREG    number of regions.
* ICREB   number of outer surfaces where outgoing leakage occurs. If
*         ICREB=0, perfect particle balance is assumed.
* NUNKNO  number of unknowns per energy group including spherical
*         harmonic terms, interface currents and fundamental currents.
* NANIS   maximum cross section Legendre order.
* MATCOD  mixture indices.
* VOL     volumes.
* KEYFLX  index of flux components in unknown vector.
* XSTRC   transport-corrected macroscopic total cross sections.
* XSDIA   transport-corrected macroscopic within-group scattering cross
*         sections.
* XCSOU   source for system of unknown.
* IGDEB   first non-converged group.
* B2      directional buckling.
* DDD     leakage coefficients.
* KEYCUR  index for currents position in FUNKNO. Used if ICREB.GT.0.
* MATALB  albedo indices. Used if ICREB.GT.0.
* ALBEDO  albedo array. Used if ICREB.GT.0.
* SURFAC  numerical surfaces. Used if ICREB.GT.0.
*
*Parameters: input/output
* FUNKNO  neutron flux.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMACR
      INTEGER    NGRP,ILEAK,NMAT,NREG,ICREB,NUNKNO,NANIS,MATCOD(NREG),
     >           KEYFLX(NREG),IGDEB,KEYCUR(ICREB),MATALB(ICREB)
      REAL       VOL(NREG),FUNKNO(NUNKNO,NGRP),XSTRC(0:NMAT,NGRP),
     >           XSDIA(0:NMAT,0:NANIS,NGRP),B2(4),DDD(NGRP),
     >           ALBEDO(6),SURFAC(ICREB)
      DOUBLE PRECISION XCSOU(NGRP)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPMACR,KPMACR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,NJJ,IPOS
      REAL, ALLOCATABLE, DIMENSION(:) :: XSCAT
      REAL, ALLOCATABLE, DIMENSION(:,:) :: REBAL
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IJJ(0:NMAT),NJJ(0:NMAT),IPOS(0:NMAT))
      ALLOCATE(REBAL(NGRP,NGRP+1),XSCAT(0:NMAT*NGRP))
*----
*  INITIALIZE REBALANCE MATRIX
*----
      NGREB=NGRP-IGDEB+1
      DO 30 I=1,NGREB
        DO 31 J=1,NGREB+1
          REBAL(I,J)=0.0
 31     CONTINUE
 30   CONTINUE
*----
*  CREATE REBALANCE MATRIX
*----
      JPMACR=LCMGID(IPMACR,'GROUP')
      DO 70 IGR=IGDEB,NGRP
        IOFF=IGR-IGDEB+1
        KPMACR=LCMGIL(JPMACR,IGR)
*----
*  READ SCATT X-SECTIONS.
*----
        CALL LCMLEN(KPMACR,'NJJS00',ILCMLN,ITYLCM)
        IF(ILCMLN.NE.NMAT) THEN
          CALL LCMLIB(KPMACR)
          CALL XABORT('FLUBAL: READ ERROR ON LCM RECORD = NJJS00')
        ELSE
          CALL LCMGET(KPMACR,'NJJS00',NJJ(1))
          CALL LCMGET(KPMACR,'IJJS00',IJJ(1))
          CALL LCMGET(KPMACR,'IPOS00',IPOS(1))
          CALL LCMGET(KPMACR,'SCAT00',XSCAT(1))
        ENDIF
*----
*  FIXE + FISSION NEUTRON SOURCES
*----
        REBAL(IOFF,NGREB+1)=REAL(XCSOU(IGR))
*----
* SUM OVER SURFACES
*----
        DO 35 ISUR=1,ICREB
           IND=KEYCUR(ISUR)
           REBAL(IOFF,IOFF)=REBAL(IOFF,IOFF)+
     >     (1.0-ALBEDO(-MATALB(ISUR)))*FUNKNO(IND,IGR)*SURFAC(ISUR)
 35     CONTINUE
*----
*  SUM OVER REGIONS
*----
        DO 60 IREG=1,NREG
          IBM=MATCOD(IREG)
          IF(IBM.EQ.0) GO TO 60
          IND=KEYFLX(IREG)
*----
*  INCLUDE SCATTERING SOURCES FROM CONVERGED FLUX IN REBALANCE SOURCE
*----
            NGSCAT=NJJ(IBM)
            IFSCAT=IJJ(IBM)-NGSCAT+1
            ISCATP=IPOS(IBM)-1+NGSCAT
            DO 40 JGR=IFSCAT,IGDEB-1
               REBAL(IOFF,NGREB+1)=REBAL(IOFF,NGREB+1)+
     >           FUNKNO(IND,JGR)*XSCAT(ISCATP)*VOL(IREG)
              ISCATP=ISCATP-1
 40         CONTINUE
*----
*  INCLUDE SCATTERING SOURCES FROM NON CONVERGED FLUX IN REBALANCE
*  MATRIX
*----
            IF(IFSCAT.LT.IGDEB) THEN
              NGSCAT=NGSCAT+IFSCAT-IGDEB
              IFSCAT=IGDEB
            ENDIF
            ISCATP=IPOS(IBM)-1+NGSCAT
            DO 50 JGR=IFSCAT,IJJ(IBM)
              IF(JGR.EQ.IGR) THEN
                REBAL(IOFF,IOFF)=REBAL(IOFF,IOFF)+FUNKNO(IND,IGR)
     >          *(XSTRC(IBM,IGR)-XSDIA(IBM,0,IGR))*VOL(IREG)
              ELSE
                REBAL(IOFF,JGR-IGDEB+1)=REBAL(IOFF,JGR-IGDEB+1)
     >          -FUNKNO(IND,JGR)*XSCAT(ISCATP)*VOL(IREG)
              ENDIF
              ISCATP=ISCATP-1
 50         CONTINUE
 60     CONTINUE
*----
*  FOR ALL REGIONS ADD CONTRIBUTION DUE TO DB2 TERM
*----
        IF(ILEAK.LT.5.AND.ILEAK.GT.0) THEN
          DO 61 IREG=1,NREG
            IND=KEYFLX(IREG)
            REBAL(IOFF,IOFF)=REBAL(IOFF,IOFF)
     >         +FUNKNO(IND,IGR)*DDD(IGR)*B2(4)*VOL(IREG)
 61       CONTINUE
        ELSE IF(ILEAK.EQ.5) THEN
          DO 65 IREG=1,NREG
            IND=KEYFLX(IREG)
            REBAL(IOFF,IOFF)=REBAL(IOFF,IOFF)
     >         +FUNKNO(NUNKNO/2+IND,IGR)*B2(4)*VOL(IREG)
 65       CONTINUE
        ELSE IF(ILEAK.GT.6) THEN
          DO 66 IREG=1,NREG
            IND=KEYFLX(IREG)
            REBAL(IOFF,IOFF)=REBAL(IOFF,IOFF)
     >         +(FUNKNO(NUNKNO/4+IND,IGR)*B2(1)
     >         +FUNKNO(NUNKNO/2+IND,IGR)*B2(2)
     >         +FUNKNO(3*NUNKNO/4+IND,IGR)*B2(3))*VOL(IREG)
 66       CONTINUE
        ENDIF
 70   CONTINUE
*----
*  SOLVE REBALANCE EQUATIONS
*----
      CALL ALSB(NGREB,1,REBAL,IER,NGRP)
      IF(IER.NE.0) THEN
        WRITE(6,'(/36H FLUBAL SINGULAR REBALANCING MATRIX.)')
        GO TO 100
      ENDIF
*----
*  REBALANCE FLUXES
*----
      DO 90 IGR=IGDEB,NGRP
        IOFF=IGR-IGDEB+1
        DO 80 IND=1,NUNKNO
           FUNKNO(IND,IGR)=FUNKNO(IND,IGR)*REBAL(IOFF,NGREB+1)
 80     CONTINUE
 90   CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
100   DEALLOCATE(XSCAT,REBAL)
      DEALLOCATE(IPOS,NJJ,IJJ)
      RETURN
      END
