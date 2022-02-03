*DECK MCGDSD
      SUBROUTINE MCGDSD(N,NMU,LPS,NFUNL,NMOD,NGEFF,WEI2D,TRHAR,
     1                  H2D,ZMU,WZMU,NOMCEL,NZON,NFI,NREG,NDIM,M,IS,JS,
     2                  PJJ,PSJ,LPJJAN,NPJJM,PJJIND,SIGAL,MUST,MODST,
     3                  PHI1,PHI2,PJJX,PJJY,PJJZ,PJJXI,PJJYI,PJJZI,
     4                  CAZ0,PSJX,PSJY,PSJZ)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of the PJJ and PSJ as well as directional values for
* TIBERE.
*
*Copyright:
* Copyright (C) 2019 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): S. Musongela
*
*Parameters: input
* N       number of elements in the current track.
* NMU     order of the polar quadrature set.
* LPS     first dimension of PSJ.
* NFUNL   number of moments of the flux (in 2D : NFUNL=NANI*(NANI+1)/2).
* NMOD    first dimension of ISGNR.
* NGEFF   number of energy groups to process.
* NFI     total number of volumes for which specific values
*         of the neutron flux and reactions rates are required.
* NREG    number of volumes for which specific values
*         of the neutron flux and reactions rates are required.
* NDIM    number of dimensions for the geometry.
* M       number of material mixtures.
* IS      arrays for surfaces neighbors.
* JS      JS(IS(ISOUT)+1:IS(ISOUT+1)) give the neighboring regions to
*         surface ISOUT.
* NZON    index-number of the mixture type assigned to each volume.
* TRHAR   spherical harmonics components for this angle in the plane.
* WEI2D   track weight.
* NOMCEL  integer tracking elements.
* H2D     real tracking elements.
* ZMU     polar quadrature set.
* WZMU    polar quadrature set.
* LPJJAN  flag for the calculation of anisotropic moments of the pjj.
* NPJJM   number of pjj modes to store for LPJJAN option.
* PJJIND  index of the modes for LPJJAN option.
* SIGAL   albedos and total cross sections array.
* MUST    polar index in TRHAR for 3D geometry.
* MODST   starting angular mode index.
* CAZ0    cosines of the tracking polar angles in 3D.
* PHI1    first cosine of the tracking azimuthal angle.
* PHI2    second cosine of the tracking azimuthal angle.
*
*Parameters: input/output
* PJJ     collision probabilities.
* PJJX    collision probabilities for TIBERE.
* PJJY    collision probabilities for TIBERE.
* PJJZ    collision probabilities for TIBERE.
* PJJXI   collision probabilities for TIBERE.
* PJJYI   collision probabilities for TIBERE.
* PJJZI   collision probabilities for TIBERE.
* PSJ     escape probabilities.
* PSJX    escape probabilities for TIBERE.
* PSJY    escape probabilities for TIBERE.
* PSJZ    escape probabilities for TIBERE.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NGEFF,N,NMU,LPS,NFUNL,NMOD,NOMCEL(N),NZON(NFI),NFI,M,
     1 NREG,NDIM,IS(NFI-NREG+1),JS(LPS),NPJJM,PJJIND(NPJJM,2),MUST,
     2 MODST
      DOUBLE PRECISION WEI2D,ZZZ,H2D(N)
      REAL TRHAR(NMU,NFUNL,NMOD),ZMU(NMU),WZMU(NMU),PSJ(LPS,NGEFF),
     1 SIGAL(-6:M,NGEFF),PSJX(LPS,NGEFF),PSJY(LPS,NGEFF),
     2 PSJZ(LPS,NGEFF)
      DOUBLE PRECISION PJJ(NREG,NPJJM,NGEFF),PHI1,PHI2,CAZ0
      DOUBLE PRECISION PJJX(NREG,NPJJM,NGEFF),PJJY(NREG,NPJJM,NGEFF),
     > PJJZ(NREG,NPJJM,NGEFF),PJJXI(NREG,NPJJM,NGEFF),
     > PJJYI(NREG,NPJJM,NGEFF),PJJZI(NREG,NPJJM,NGEFF)
      LOGICAL LPJJAN
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION W,OMEGA2(3)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: HG
*----
*  CALCULATION OF COEFFICIENTS
*----
      IF (NDIM.EQ.3) THEN
*     3D calculation -> no loop over a polar angle   
         DO II=1,NGEFF
*           MCGDSCB: Step-Characteristics Scheme with Tabulated
*                    exponentials
            OMEGA2(3)=CAZ0*CAZ0
            ZZZ=1.0D0/SQRT(1.0D0-OMEGA2(3))
            OMEGA2(1)=(PHI1/ZZZ)**2
            OMEGA2(2)=(PHI2/ZZZ)**2
            W=WEI2D
            CALL MCGDSCB(M,N,LPS,IS,JS,H2D,NOMCEL,NZON,SIGAL(0,II),W,
     1           NFI,NREG,PJJ(1,1,II),PSJ(1,II),MUST,NMU,NFUNL,NMOD,
     2           NPJJM,TRHAR,LPJJAN,PJJIND,MODST,OMEGA2,PJJX(1,1,II),
     3           PJJY(1,1,II),PJJZ(1,1,II),PJJXI(1,1,II),
     4           PJJYI(1,1,II),PJJZI(1,1,II),PSJX(1,II),PSJY(1,II),
     5           PSJZ(1,II))
         ENDDO
      ELSE
*     2D calculation -> loop over the polar angle
         ALLOCATE(HG(N))
         DO IMU=1,NMU
            OMEGA2(1)=(PHI1/ZMU(IMU))**2
            OMEGA2(2)=(PHI2/ZMU(IMU))**2
            OMEGA2(3)=(1.0-1.0/ZMU(IMU)**2)
            ZMUI=ZMU(IMU)
            W=WEI2D*WZMU(IMU)
            DO I=1,N
               IF(NZON(NOMCEL(I)).GE.0) THEN
                  HG(I)=H2D(I)*ZMUI
               ENDIF
            ENDDO            
            DO II=1,NGEFF
               CALL MCGDSCB(M,N,LPS,IS,JS,HG,NOMCEL,NZON,SIGAL(0,II),W,
     1              NFI,NREG,PJJ(1,1,II),PSJ(1,II),IMU,NMU,NFUNL,NMOD,
     2              NPJJM,TRHAR,LPJJAN,PJJIND,MODST,OMEGA2,PJJX(1,1,II),
     3              PJJY(1,1,II),PJJZ(1,1,II),PJJXI(1,1,II),
     4              PJJYI(1,1,II),PJJZI(1,1,II),PSJX(1,II),PSJY(1,II),
     5              PSJZ(1,II))
            ENDDO
         ENDDO
         DEALLOCATE(HG)
      ENDIF
*
      RETURN
      END
