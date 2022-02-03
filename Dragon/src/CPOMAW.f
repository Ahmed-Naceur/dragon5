*DECK CPOMAW
      SUBROUTINE CPOMAW(IPCPO ,IPRINT,NGROUP,NL    ,NPROC ,INDPRO,
     >                  ITYPRO,DXSMAC,DSCMAC,DXSREM,DSCREM,DISFC ,
     >                  DMJMAC,IFCDIS,DISFAC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Put macroscopic cross section on Compo.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): G. Marleau
*
*Parameters: input
* IPCPO   pointer to the Compo.
* IPRINT  print parameter. Equal to zero for no print.
* NGROUP  number of groups condensed.
* NL      number of Legendre orders.
* NPROC   number of microscopic xs to process.
* INDPRO  identifier for xs processing.
* ITYPRO  identifier for xs processed .
* DXSMAC  macroscopic averaged region/group x-s.
* DSCMAC  macroscopic scattering.
* DXSREM  removed averaged region/group x-s.
* DSCREM  removed scattering rates.
* DISFC   disadvantage factor.
* DMJMAC  energy.
* IFCDIS  discontinuity factor present (1) or absent.
* DISFAC  discontinuity factors.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR)      IPCPO
      INTEGER          IPRINT,NGROUP,NL,NPROC,INDPRO(NPROC),
     >                 ITYPRO(NPROC)
      REAL             DISFC(NGROUP)
      DOUBLE PRECISION DXSMAC(NGROUP,NPROC),
     >                 DSCMAC(NGROUP,NGROUP,NL),
     >                 DXSREM(NGROUP,NPROC),
     >                 DSCREM(NGROUP,NGROUP,NL),
     >                 DMJMAC
      INTEGER          IFCDIS
      DOUBLE PRECISION DISFAC(2,NGROUP,3)
*----
*  LOCAL PARAMETERS
*----
      INTEGER          NDPROC
      REAL             CUTOFF
      PARAMETER       (NDPROC=20,CUTOFF=1.0E-7)
      INTEGER          IXSR,JXSR,KXSR,IL,IGR,JGR,IORD
      REAL             CUTLIM
      DOUBLE PRECISION DNUFI,DNUFT
*----
*  ALLOCATABLE ARRAYS
*   XSREC   micro vector xs
*   XSCAT   compress scattering data
*   DISTMP  temporary storage for discontinuity factors              
*----
      REAL, ALLOCATABLE, DIMENSION(:,:) :: XSREC,DISTMP
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: XSCAT
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(XSREC(NGROUP,NPROC),XSCAT(NGROUP,NGROUP,NL),
     > DISTMP(2,NGROUP))
*----
*   SAVE AVERAGE XS
*----
      CALL XDRSET(ITYPRO,NPROC,1)
      CALL LCMSIX(IPCPO,'MACR',1)
      XSREC(:NGROUP,:NPROC)=0.0
      XSCAT(:NGROUP,:NGROUP,:NL)=0.0
      DO 100 IXSR=1,4
        IF(INDPRO(IXSR).GT.0.AND.MOD(ITYPRO(IXSR),2).EQ.1) THEN
          DO 101 IGR=1,NGROUP
            XSREC(IGR,IXSR)=REAL(DXSMAC(IGR,IXSR)-DXSREM(IGR,IXSR))
            CUTLIM=ABS(REAL(DXSMAC(IGR,IXSR))*CUTOFF)
            IF(XSREC(IGR,IXSR).LT.CUTLIM)  XSREC(IGR,IXSR)=0.0
 101      CONTINUE
        ENDIF
 100  CONTINUE
      IXSR=5
      IF(INDPRO(IXSR).GT.0.AND.MOD(ITYPRO(IXSR),2).EQ.1) THEN
        JXSR=3
        KXSR=16
        DNUFI=0.0D0
        DNUFT=0.0D0
        DO 120 IGR=1,NGROUP
          DNUFT=DNUFT+DXSMAC(IGR,KXSR)*DXSMAC(IGR,JXSR)
          DNUFI=DNUFI+(DXSMAC(IGR,KXSR)-DXSREM(IGR,KXSR))
     >      *(DXSMAC(IGR,JXSR)-DXSREM(IGR,JXSR))
 120    CONTINUE
        CUTLIM=ABS(REAL(DNUFT)*CUTOFF)
        IF(REAL(DNUFI).GT.CUTLIM) THEN
          DNUFI=1.0D0/DNUFI
          DNUFT=1.0D0/DNUFT
          DO 130 IGR=1,NGROUP
            XSREC(IGR,IXSR)=REAL(DNUFI*(DXSMAC(IGR,IXSR)
     >        -DXSREM(IGR,IXSR)))
 130      CONTINUE
        ENDIF
      ENDIF
      IXSR=6
      IF(INDPRO(IXSR).GT.0.AND.MOD(ITYPRO(IXSR),2).EQ.1) THEN
        JXSR=4
        KXSR=3
        DO 140 IGR=1,NGROUP
          DNUFI=DXSMAC(IGR,JXSR)-DXSREM(IGR,JXSR)
          CUTLIM=ABS(REAL(DXSMAC(IGR,JXSR))*CUTOFF)
          IF(REAL(DNUFI).GT.CUTLIM) THEN
            XSREC(IGR,IXSR)=REAL((DXSMAC(IGR,KXSR)
     >        -DXSREM(IGR,KXSR))/DNUFI)
          ENDIF
 140    CONTINUE
      ENDIF
      DO 150 IXSR=7,NDPROC
        IF(INDPRO(IXSR).GT.0.AND.MOD(ITYPRO(IXSR),2).EQ.1) THEN
          DO 160 IGR=1,NGROUP
            XSREC(IGR,IXSR)=REAL(DXSMAC(IGR,IXSR)-DXSREM(IGR,IXSR))
            CUTLIM=ABS(REAL(DXSMAC(IGR,IXSR))*CUTOFF)
            IF(XSREC(IGR,IXSR).LT.CUTLIM)  XSREC(IGR,IXSR)=0.0
 160      CONTINUE
        ENDIF
 150  CONTINUE
      IL=0
      DO 170 IXSR=NDPROC+1,NDPROC+NL
        IL=IL+1
        IF(INDPRO(IXSR).GT.0.AND.MOD(ITYPRO(IXSR),2).EQ.1) THEN
          DO 180 IGR=1,NGROUP
            XSREC(IGR,IXSR)=REAL(DXSMAC(IGR,IXSR)-DXSREM(IGR,IXSR))
            CUTLIM=ABS(REAL(DXSMAC(IGR,IXSR))*CUTOFF)
            IF(ABS(XSREC(IGR,IXSR)).LT.CUTLIM)
     >         XSREC(IGR,IXSR)=0.0
            DO 190 JGR=1,NGROUP
              XSCAT(IGR,JGR,IL)=REAL(DSCMAC(IGR,JGR,IL)
     >          -DSCREM(IGR,JGR,IL))
              CUTLIM=ABS(REAL(DSCMAC(IGR,JGR,IL))*CUTOFF)
              IF(ABS(XSCAT(IGR,JGR,IL)).LT.CUTLIM)
     >          XSCAT(IGR,JGR,IL)=0.0
 190        CONTINUE
 180      CONTINUE
        ENDIF
 170  CONTINUE
*----
*  COMPUTE AVERAGED ENERGY PER FISSION
*----
      JXSR=4
      KXSR=16
      DNUFI=0.0D0
      DNUFT=0.0D0
      DO 200 IGR=1,NGROUP
        DNUFT=DNUFT+DXSMAC(IGR,KXSR)*DXSMAC(IGR,JXSR)
        DNUFI=DNUFI+(DXSMAC(IGR,KXSR)-DXSREM(IGR,KXSR))
     >    *(DXSMAC(IGR,JXSR)-DXSREM(IGR,JXSR))
 200  CONTINUE
      CUTLIM=ABS(REAL(DNUFT)*CUTOFF)
      IF(REAL(DNUFI).GT.CUTLIM) THEN
        DMJMAC=DMJMAC/DNUFI
      ELSE
        DMJMAC=0.0D0
      ENDIF
*----
*  SAVE CPO MICRO
*----
      IORD=1
      CALL XDRLGS(IPCPO,1,IPRINT,0,NL-1,IORD,NGROUP,XSREC(1,NDPROC+1),
     >            XSCAT,ITYPRO(NDPROC+1))
      CALL CPOLGX(IPCPO,1,IPRINT,IORD,NGROUP,INDPRO,XSREC(1,1),ITYPRO)
      CALL LCMSIX(IPCPO,'MACR',2)
      IXSR=NDPROC+NL+1
      DO 210 IGR=1,NGROUP
        XSREC(IGR,1)=REAL(DXSMAC(IGR,IXSR))
 210  CONTINUE
      CALL LCMPUT(IPCPO,'OVERV',NGROUP,2,XSREC)
      IXSR=16
      DO 220 IGR=1,NGROUP
        XSREC(IGR,1)=REAL(DXSMAC(IGR,IXSR))
 220  CONTINUE
      CALL LCMPUT(IPCPO,'FLUX-INTG',NGROUP,2,XSREC)
      CALL LCMPUT(IPCPO,'FLUXDISAFACT',NGROUP,2,DISFC)
      IF(IFCDIS .EQ. 1) THEN
        CALL XDRSDB(2*NGROUP,DISTMP,DISFAC(1,1,1),1)
        CALL LCMPUT(IPCPO,'DISFACX',2*NGROUP,2,DISTMP)
        CALL XDRSDB(2*NGROUP,DISTMP,DISFAC(1,1,2),1)
        CALL LCMPUT(IPCPO,'DISFACY',2*NGROUP,2,DISTMP)
        CALL XDRSDB(2*NGROUP,DISTMP,DISFAC(1,1,3),1)
        CALL LCMPUT(IPCPO,'DISFACZ',2*NGROUP,2,DISTMP)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(DISTMP,XSCAT,XSREC)
      RETURN
      END
