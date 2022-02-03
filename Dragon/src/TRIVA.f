*DECK TRIVA
      SUBROUTINE TRIVA(IPSYS,IPTRK,IMPX,NREG,NBMIX,NANI,NW,MAT,VOL,
     1 SIGT0,SIGW0,DIFF)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembly of one-speed finite-difference or finite-element matrices
* for a discretization of the 3D diffusion or SPN equation.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPSYS   pointer to the system matrices.
* IPTRK   pointer to the tracking (L_TRACK signature).
* IMPX    print flag (equal to zero for no print).
* NREG    total number of merged regions for which specific values
*         of the neutron flux and reactions rates are required.
* NBMIX   number of mixtures.
* NANI    number of Legendre orders for the scattering cross sections.
* NW      type of weighting for P1 cross section info (=0 P0 ; =1 P1).
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* SIGT0   P0 and P1 total macroscopic cross sections ordered by mixture.
* SIGW0   within-group scattering macroscopic cross section ordered
*         by mixture.
* DIFF    diffusion coefficients ordered by mixture.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPSYS,IPTRK
      INTEGER IMPX,NREG,NBMIX,NANI,NW,MAT(NREG)
      REAL VOL(NREG),SIGT0(0:NBMIX,NW+1),SIGW0(0:NBMIX,NANI),
     1 DIFF(0:NBMIX)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40)
      INTEGER ISTATE(NSTATE),IGB(8)
      LOGICAL LBIHET
      CHARACTER NAMP*12,TEXT10*10
      REAL, ALLOCATABLE, DIMENSION(:) :: GAMMA
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SGD,SGDI
      PARAMETER(TEXT10='A001001')
*----
*  RECOVER TRIVAC SPECIFIC TRACKING INFORMATION
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      LBIHET=ISTATE(40).NE.0
      IF(LBIHET) THEN
         CALL LCMSIX(IPTRK,'BIHET',1)
         CALL LCMGET(IPTRK,'PARAM',IGB)
         IF(NREG.NE.IGB(3)) CALL XABORT('TRIVA: INVALID VALUE OF NREG('
     1   //'1).')
         CALL LCMSIX(IPTRK,' ',2)
      ELSE
         IF(NREG.NE.ISTATE(1)) CALL XABORT('TRIVA: INVALID VALUE OF NR'
     1   //'EG(2).')
      ENDIF
      ICHX=ISTATE(12)
      NLF=ISTATE(30)
      ISCAT=ABS(ISTATE(32))
*----
*  RECOVER PHYSICAL ALBEDO FUNCTIONS.
*----
      CALL LCMLEN(IPSYS,'ALBEDO-FU',NALBP,ITYLCM)
      IF(NALBP.GT.0) THEN
         ALLOCATE(GAMMA(NALBP))
         CALL LCMGET(IPSYS,'ALBEDO-FU',GAMMA)
      ENDIF
*----
*  COMPUTE THE WITHIN-GROUP SYSTEM MATRICES (LEAKAGE AND REMOVAL).
*  ASSEMBLY OF THE ADI SPLITTED SYSTEM MATRICES
*----
      IF(NLF.EQ.0) THEN
*----
*  ++++ DIFFUSION THEORY ++++
*----
         IF(NANI.GT.1) THEN
            CALL XABORT('TRIVA: SPN MACRO-CALCULATION EXPECTED(1).')
         ENDIF
         ALLOCATE(SGD(NBMIX,4))
         DO 10 IBM=1,NBMIX
         SGD(IBM,1)=DIFF(IBM)
         SGD(IBM,2)=DIFF(IBM)
         SGD(IBM,3)=DIFF(IBM)
         SGD(IBM,4)=SIGT0(IBM,1)-SIGW0(IBM,1)
   10    CONTINUE
*----
*  ASSEMBLY OF A SINGLE-GROUP SYSTEM MATRIX WITH LEAKAGE AND REMOVAL
*  CROSS SECTIONS.
*----
         CALL TRIASM(TEXT10,IPTRK,IPSYS,IMPX,NBMIX,NREG,NALBP,0,MAT,
     1   VOL,GAMMA,SGD,SGD)
         DEALLOCATE(SGD)
      ELSE
*----
*  ++++ PN OR SPN THEORY ++++
*----
         IF(NLF.LT.2) THEN
            CALL XABORT('TRIVA: PN OR SPN KEYWORD EXPECTED.')
         ELSE IF(ICHX.NE.2) THEN
            CALL XABORT('TRIVA: DISCRETIZATION NOT AVAILABLE.')
         ENDIF
         NAN=MIN(ISCAT,NANI)+1
         ALLOCATE(SGD(NBMIX,NAN),SGDI(NBMIX,NAN))
         DO 30 IL=0,NAN-1
         DO 20 IBM=1,NBMIX
         IF(IL.LE.NW) THEN
            GARS=SIGT0(IBM,IL+1)
         ELSE IF((NW.GE.1).AND.(MOD(IL,2).EQ.1)) THEN
            GARS=SIGT0(IBM,2)
         ELSE
            GARS=SIGT0(IBM,1)
         ENDIF
         IF(IL.LE.NAN-2) GARS=GARS-SIGW0(IBM,IL+1)
         SGD(IBM,IL+1)=GARS
         IF(GARS.NE.0.0) THEN
            SGDI(IBM,IL+1)=1.0/GARS
         ELSE
            SGDI(IBM,IL+1)=1.0E10
         ENDIF
   20    CONTINUE
         WRITE(NAMP,'(4HSCAR,I2.2,6H001001)') IL
         CALL LCMPUT(IPSYS,NAMP,NBMIX,2,SGD(1,IL+1))
         WRITE(NAMP,'(4HSCAI,I2.2,6H001001)') IL
         CALL LCMPUT(IPSYS,NAMP,NBMIX,2,SGDI(1,IL+1))
   30    CONTINUE
         CALL XDISET(ISTATE,NSTATE,0)
         ISTATE(7)=NBMIX
         ISTATE(8)=NAN
         CALL LCMPUT(IPSYS,'STATE-VECTOR',NSTATE,1,ISTATE)
*----
*  ASSEMBLY OF A SINGLE-GROUP SYSTEM MATRIX WITH LEAKAGE AND REMOVAL
*  CROSS SECTIONS FOR THE SIMPLIFIED PN METHOD.
*----
         CALL TRIASN(TEXT10,IPTRK,IPSYS,IMPX,NBMIX,NREG,NAN,NALBP,0,
     1   MAT,VOL,GAMMA,SGD,SGDI)
         DEALLOCATE(SGDI,SGD)
      ENDIF
      IF(NALBP.GT.0) DEALLOCATE(GAMMA)
      IF(IMPX.GT.2) CALL LCMLIB(IPSYS)
      IF(IMPX.GT.10) CALL LCMVAL(IPSYS,' ')
      RETURN
      END
