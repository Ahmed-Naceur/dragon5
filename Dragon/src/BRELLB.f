*DECK BRELLB
      SUBROUTINE BRELLB(IPMAC1,NC,NG,NMIX1,ENER,JXM,FHETXM,IPRINT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Implement the 1D Lefebvre-Lebigot reflector model.
*
*Copyright:
* Copyright (C) 2021 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert and V. Salino
*
*Parameters: input
* IPMAC1  nodal macrolib.
* NC      number of sn macrolibs (=2: Lefebvre-Lebigot method).
* NG      number of energy groups.
* NMIX1   number of mixtures in the nodal calculation.
* ENER    energy limits.
* JXM     left boundary currents.
* FHETXM  left boundary fluxes.
* IPRINT  edition flag.
*
*Reference:
*  Edwige Richebois, 'Calculs de coeur REP en transport 3D', PhD,
*  Universite Aix-Marseille, 1999 (p.193).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAC1
      INTEGER NC,NG,NMIX1,IPRINT
      REAL ENER(NG+1),JXM(NMIX1,NG,NC),FHETXM(NMIX1,NG,NC)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      INTEGER ISTATE(NSTATE)
      REAL BF1BF2(2),CUR1BF1(2),CUR2BF1(2),SIGT(2),DIF(2)
      TYPE(C_PTR) JPMAC1,KPMAC1
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,NJJ,IPOS
      REAL, ALLOCATABLE, DIMENSION(:) :: WORK
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SCAT
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IJJ(NMIX1),NJJ(NMIX1),IPOS(NMIX1),WORK(NG*NMIX1),
     1 SCAT(NG,NG))
*----
*  RECOVER FLUX, MACROSCOPIC CROSS SECTIONS AND DIFFUSION COEFFICIENTS
*----
      IF(NC.NE.2) CALL XABORT('BRELLB: NC=2 EXPECTED.')
      IF(NG.NE.2) CALL XABORT('BRELLB: NG=2 EXPECTED.')
      IF(NMIX1.NE.1) CALL XABORT('BRELLB: NMIX1=1 EXPECTED.')
*----
*  COMPUTE EQUIVALENT REFLECTOR
*----
      IBM=1
      DO IC=1,NC
        BF1BF2(IC)=FHETXM(IBM,2,IC)/FHETXM(IBM,1,IC)
        CUR1BF1(IC)=JXM(IBM,1,IC)/FHETXM(IBM,1,IC)
        CUR2BF1(IC)=JXM(IBM,2,IC)/FHETXM(IBM,1,IC)
      ENDDO
      DIF(1)=1.3
      DIF(2)=0.4
      PENTE=(CUR2BF1(1)-CUR2BF1(2))/(BF1BF2(1)-BF1BF2(2))
      ORDORI=CUR2BF1(1)-PENTE*BF1BF2(1)
      R1=CUR1BF1(1)
      R2=PENTE
      R3=-ORDORI
      IF(IPRINT.GT.0) WRITE(6,10) R1,R2,R3
      SIGT(2)=R2*R2/DIF(2)
      SIGSLW=(SQRT(DIF(2)/DIF(1))*R1+SQRT(DIF(1)*SIGT(2)))*R3/
     1 SQRT(DIF(1)*DIF(2))
      SIGT(1)=R1*R1/DIF(1)
      IF(SIGSLW.LT.0.0) THEN
        CALL XABORT('BRELLB: Negative fast SIGS00 XS.')
      ELSE IF(SIGT(1)-SIGSLW.LT.0.0) THEN
        CALL XABORT('BRELLB: Negative fast absorption XS.')
      ENDIF ;
      SCAT(:,:)=0.0
      SCAT(2,1)=SIGSLW
      IF(IPRINT.GT.0) THEN
        WRITE(6,'(/40H BRELLB: LEFEBVRE-LEBIGOT CROSS SECTIONS)')
        WRITE(6,'(/12H MIXTURE=  1)')
        WRITE(6,20) 'DIFF',DIF(:NG)
        WRITE(6,20) 'SIGT',SIGT(:NG)
        WRITE(6,20) 'SCAT',SCAT(:NG,:NG)
      ENDIF
*----
*  SAVE THE OUTPUT NODAL MACROLIB
*----
      IBM=1
      ISTATE(:)=0
      ISTATE(1)=NG
      ISTATE(2)=NMIX1
      ISTATE(3)=1
      ISTATE(9)=1  ! diffusion coefficient information
      CALL LCMPUT(IPMAC1,'STATE-VECTOR',NSTATE,1,ISTATE)
      CALL LCMPUT(IPMAC1,'ENERGY',NG+1,2,ENER)
      WORK(1)=1.0
      CALL LCMPUT(IPMAC1,'VOLUME',NMIX1,2,WORK)
      JPMAC1=LCMLID(IPMAC1,'GROUP',NG)
      DO IGR=1,NG
        KPMAC1=LCMDIL(JPMAC1,IGR)
        WORK(1)=1.0
        CALL LCMPUT(KPMAC1,'FLUX-INTG',NMIX1,2,WORK)
        CALL LCMPUT(KPMAC1,'NTOT0',NMIX1,2,SIGT(IGR))
        WORK(1)=0.0
        CALL LCMPUT(KPMAC1,'SIGW00',NMIX1,2,WORK)
        CALL LCMPUT(KPMAC1,'DIFF',NMIX1,2,DIF(IGR))
        IPOSDE=0
        DO J=1,NMIX1
          IBM=1
          J2=IGR
          J1=IGR
          DO JGR=1,NG
            IF(SCAT(IGR,JGR).NE.0.0) THEN
              J2=MAX(J2,JGR)
              J1=MIN(J1,JGR)
            ENDIF
          ENDDO
          NJJ(J)=J2-J1+1
          IJJ(J)=J2
          IPOS(J)=IPOSDE+1
          DO JGR=J2,J1,-1
            IPOSDE=IPOSDE+1
            IF(IPOSDE.GT.NG*NMIX1) CALL XABORT('BRELLB: SCAT OVERFLOW.')
            WORK(IPOSDE)=SCAT(IGR,JGR)
          ENDDO
        ENDDO
        CALL LCMPUT(KPMAC1,'SCAT00',IPOSDE,2,WORK)
        CALL LCMPUT(KPMAC1,'NJJS00',NMIX1,1,NJJ)
        CALL LCMPUT(KPMAC1,'IJJS00',NMIX1,1,IJJ)
        CALL LCMPUT(KPMAC1,'IPOS00',NMIX1,1,IPOS)
      ENDDO
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(SCAT,WORK,IPOS,NJJ,IJJ)
      RETURN
*
   10 FORMAT(12H BRELLB: R1=,1P,E12.4,4H R2=,E12.4,4H R3=,E12.4)
   20 FORMAT(1X,A9,1P,10E12.4,/(10X,10E12.4))
      END
