*DECK SYBILP
      SUBROUTINE SYBILP (IPTRK,IMPX,NREG,NBMIX,MAT,VOL,SIGT0,SIGW0,
     1 NELPIJ,PIJ,ILK)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of the collision probabilities for Sybil.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPTRK   pointer to the tracking (L_TRACK signature).
* IMPX    print flag (equal to zero for no print).
* NREG    total number of merged blocks for which specific values
*         of the neutron flux and reactions rates are required.
* NBMIX   number of mixtures (NBMIX=max(MAT(i))).
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* SIGT0   total macroscopic cross sections ordered by mixture.
* SIGW0   P0 within-group scattering macroscopic cross sections
*         ordered by mixture.
* NELPIJ  number of elements in pij matrix.
* ILK     leakage flag (=.true. if neutron leakage through external
*         boundary is present).
*
*Parameters: output
* PIJ     reduced and symmetrized collision probabilities.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      LOGICAL ILK
      TYPE(C_PTR) IPTRK
      INTEGER IMPX,NREG,NBMIX,MAT(NREG),NELPIJ
      REAL VOL(NREG),SIGT0(0:NBMIX),SIGW0(0:NBMIX),PIJ(NELPIJ)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (EPS1=1.0E-4,NSTATE=40)
      INTEGER JPAR(NSTATE)
      REAL, ALLOCATABLE, DIMENSION(:) :: SIGT,SIGW
      REAL, ALLOCATABLE, DIMENSION(:,:) :: PP
*----
*  RECOVER SYBIL SPECIFIC PARAMETERS
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',JPAR)
      ITG=JPAR(6)
*
      ALLOCATE(SIGT(NREG),SIGW(NREG),PP(NREG,NREG))
      DO 10 I=1,NREG
      SIGT(I)=SIGT0(MAT(I))
      SIGW(I)=SIGW0(MAT(I))
   10 CONTINUE
      CALL SYBCP1(IPTRK,ITG,IMPX,NREG,SIGT,SIGW,PP)
*
      IF((IMPX.GE.10).OR.(IMPX.LT.0)) THEN
*        CHECK THE RECIPROCITY CONDITIONS.
         VOLTOT=0.0
         DO 20 I=1,NREG
         VOLTOT=VOLTOT+VOL(I)
   20    CONTINUE
         VOLTOT=VOLTOT/REAL(NREG)
         WRK=0.0
         DO 40 I=1,NREG
         DO 30 J=1,NREG
         AAA=PP(I,J)*VOL(I)
         BBB=PP(J,I)*VOL(J)
         WRK=MAX(WRK,ABS(AAA-BBB)/VOLTOT)
   30    CONTINUE
   40    CONTINUE
         IF(WRK.GE.EPS1) WRITE (6,150) WRK
*        CHECK THE CONSERVATION CONDITIONS.
         IF(.NOT.ILK) THEN
            WRK=0.0
            DO 60 I=1,NREG
            F1=1.0
            DO 50 J=1,NREG
            AAA=PP(I,J)
            F1=F1-AAA*(SIGT(J)-SIGW(J))
   50       CONTINUE
            WRK=AMAX1(WRK,ABS(F1))
   60       CONTINUE
            IF(WRK.GE.EPS1) WRITE (6,160) WRK
         ENDIF
      ENDIF
*
      IC=0
      DO 80 IKK=1,NREG
      IOF=(IKK-1)*NREG
      DO 70 JKK=1,IKK
      IC=IC+1
      PIJ(IC)=PP(JKK,IKK)*VOL(JKK)
   70 CONTINUE
   80 CONTINUE
      DEALLOCATE(PP,SIGT,SIGW)
      RETURN
*
  150 FORMAT (/50H THE SCATTERING-REDUCED PIJ DO NOT MEET THE RECIPR,
     1 25HOCITY CONDITIONS. RECIP =,1P,E10.3/)
  160 FORMAT (/50H THE SCATTERING-REDUCED PIJ DO NOT MEET THE CONSER,
     1 25HVATION CONDITIONS. LEAK =,1P,E10.3/)
      END
