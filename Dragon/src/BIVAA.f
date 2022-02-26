*DECK BIVAA
      SUBROUTINE BIVAA(IPSYS,IPTRK,IMPX,NREG,NBMIX,NANI,NW,MAT,VOL,
     1 SIGT0,SIGW0,DIFF)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembly of one-speed finite-difference or finite-element matrices
* for a discretization of the 2D diffusion equation.
*
*Copyright:
* Copyright (C) 2004 Ecole Polytechnique de Montreal
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
*         of the neutron flux and reaction rates are required.
* NBMIX   number of mixtures.
* NANI    number of Legendre orders for the scattering cross sections.
* NW      type of weighting for P1 cross section info (=0: P0 ; =1: P1).
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
      INTEGER JPAR(NSTATE),IGB(8)
      LOGICAL LBIHET
      CHARACTER NAMP*12,TEXT10*10
      REAL, ALLOCATABLE, DIMENSION(:) :: GAMMA
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SGD
      PARAMETER(TEXT10='A001001')
*----
*  RECOVER BIVAC SPECIFIC TRACKING PARAMETERS.
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',JPAR)
      LBIHET=JPAR(40).NE.0
      IF(LBIHET) THEN
         CALL LCMSIX(IPTRK,'BIHET',1)
         CALL LCMGET(IPTRK,'PARAM',IGB)
         IF(NREG.NE.IGB(3)) CALL XABORT('BIVAA: INVALID VALUE OF NREG('
     1   //'1).')
         CALL LCMSIX(IPTRK,' ',2)
      ELSE
         IF(NREG.NE.JPAR(1)) CALL XABORT('BIVAA: INVALID VALUE OF NREG'
     1   //'(2).')
      ENDIF
      NLF=JPAR(14)
      ISCAT=ABS(JPAR(16))
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
*----
      IF(NLF.EQ.0) THEN
*----
*  ++++ DIFFUSION THEORY ++++
*----
         IF(NANI.GT.1) THEN
            CALL XABORT('BIVAA: SPN MACRO-CALCULATION EXPECTED(1).')
         ENDIF
         ALLOCATE(SGD(NBMIX,3))
         DO 10 IBM=1,NBMIX
         SGD(IBM,1)=DIFF(IBM)
         SGD(IBM,2)=DIFF(IBM)
         SGD(IBM,3)=SIGT0(IBM,1)-SIGW0(IBM,1)
   10    CONTINUE
*----
*  ASSEMBLING OF A SINGLE-GROUP SYSTEM MATRIX FOR BIVAC.
*----
         CALL BIVASM(TEXT10,0,IPTRK,IPSYS,IMPX,NBMIX,NREG,NLF,3,NALBP,
     1   MAT,VOL,GAMMA,SGD)
         DEALLOCATE(SGD)
      ELSE
*----
*  ++++ PN OR SPN THEORY ++++
*----
         IF(NLF.LT.2) THEN
            CALL XABORT('BIVAA: PN OR SPN KEYWORD EXPECTED.')
         ENDIF
         NAN=MIN(ISCAT,NANI)+1
         ALLOCATE(SGD(NBMIX,2*NAN))
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
            SGD(IBM,NAN+IL+1)=1.0/GARS
         ELSE
            SGD(IBM,NAN+IL+1)=1.0E10
         ENDIF
   20    CONTINUE
         WRITE(NAMP,'(4HSCAR,I2.2,6H001001)') IL
         CALL LCMPUT(IPSYS,NAMP,NBMIX,2,SGD(1,NAN+IL+1))
         WRITE(NAMP,'(4HSCAI,I2.2,6H001001)') IL
         CALL LCMPUT(IPSYS,NAMP,NBMIX,2,SGD(1,NAN+IL+1))
   30    CONTINUE
         CALL XDISET(JPAR,NSTATE,0)
         JPAR(7)=NBMIX
         JPAR(8)=NAN
         CALL LCMPUT(IPSYS,'STATE-VECTOR',NSTATE,1,JPAR)
*----
*  ASSEMBLING OF A SINGLE-GROUP SYSTEM MATRIX FOR BIVAC.
*----
         CALL BIVASM(TEXT10,0,IPTRK,IPSYS,IMPX,NBMIX,NREG,NLF,2*NAN,
     1   NALBP,MAT,VOL,GAMMA,SGD)
         DEALLOCATE(SGD)
      ENDIF
      IF(NALBP.GT.0) DEALLOCATE(GAMMA)
      IF(IMPX.GT.2) CALL LCMLIB(IPSYS)
      IF(IMPX.GT.10) CALL LCMVAL(IPSYS,' ')
      RETURN
      END
