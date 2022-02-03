*DECK MCGDYA
      SUBROUTINE MCGDYA(NMU,ZMU,WZMU,NANGL,CAZ1,CAZ2,IL,IM,ILP,IMP,MM)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Returns the dyadic matrix of Eq. (43).
*
*Copyright:
* Copyright (C) 2015 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): A. Hebert
*
*Parameters: input
* NMU     number of polar angles.
* ZMU     polar quadrature set in 2D.
* WZMU    polar quadrature set in 2D.
* NANGL   number of azimuthal angles.
* CAZ1    first azimuthal cosines.
* CAZ2    second azimuthal cosines.
* IL      spherical harmonics index.
* IM      spherical harmonics index.
* ILP     spherical harmonics index.
* IMP     spherical harmonics index.
*
*Parameters: output
* MM      dyadic matrix.
*
*-----------------------------------------------------------------------
*
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NMU, IL, IM, ILP, IMP
      REAL ZMU(NMU),WZMU(NMU)
      DOUBLE PRECISION CAZ1(NANGL),CAZ2(NANGL),MM(2,2)
*
      REAL PNSH
      DOUBLE PRECISION DYAD(2,2)
      REAL, DIMENSION(:), ALLOCATABLE :: WW, MU, ETA, XI
      DOUBLE PRECISION, PARAMETER :: PI=3.14159265358979323846264338328
*
      ALLOCATE(WW(2*NMU*NANGL), MU(2*NMU*NANGL), ETA(2*NMU*NANGL),
     >         XI(2*NMU*NANGL))
      DO IANGLE=1,NANGL
        DO IMU=1,NMU
          MU((IANGLE-1)*NMU+IMU)=SQRT(1.0-1.0/ZMU(IMU)**2)
          ETA((IANGLE-1)*NMU+IMU)=REAL(CAZ1(IANGLE)/ZMU(IMU))
          XI((IANGLE-1)*NMU+IMU)=REAL(CAZ2(IANGLE)/ZMU(IMU))
          WW((IANGLE-1)*NMU+IMU)=WZMU(IMU)*MU((IANGLE-1)*NMU+IMU)
          MU((NANGL+IANGLE-1)*NMU+IMU)=-MU((IANGLE-1)*NMU+IMU)
          ETA((NANGL+IANGLE-1)*NMU+IMU)=-ETA((IANGLE-1)*NMU+IMU)
          XI((NANGL+IANGLE-1)*NMU+IMU)=-XI((IANGLE-1)*NMU+IMU)
          WW((NANGL+IANGLE-1)*NMU+IMU)=WW((IANGLE-1)*NMU+IMU)
        ENDDO
      ENDDO
      MM(:,:)=0.0D0 ; WSUM=0.0D0 ;
      DO IANGLE=1,2*NMU*NANGL
        DYAD=MATMUL(RESHAPE((/ ETA(IANGLE), XI(IANGLE) /),(/ 2, 1 /)),
     >              RESHAPE((/ ETA(IANGLE), XI(IANGLE) /),(/ 1, 2 /)))
        RLM=PNSH(IL,IM,MU(IANGLE),ETA(IANGLE),XI(IANGLE))
        RLMP=PNSH(ILP,IMP,MU(IANGLE),ETA(IANGLE),XI(IANGLE))
        WSUM=WSUM+WW(IANGLE)
        MM=MM+WW(IANGLE)*RLM*RLMP*DYAD
      ENDDO
      MM=(2.0D0*ILP+1.0D0)/WSUM*MM
      DEALLOCATE(WW, MU, ETA, XI)
      RETURN
      END
