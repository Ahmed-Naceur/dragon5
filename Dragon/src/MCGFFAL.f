*DECK MCGFFAL
      SUBROUTINE MCGFFAL(SUBSCH,K,KPN,M,N,H,NOM,NZON,WEIGHT,XST,SP,SM,
     1           DSP,DSM,NREG,NMU,NANI,NFUNLX,NMOD,TRHAR,KEYCUR,IMU,B,
     2           MODP,MODM,F,PHIV,DPHIV)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solution of transport equation on a finite track.
* Linear-discontinuous-characteristics approximation.
* Ray-tracing (anisotropic scattering case,'source term isolation' off).
*
*Copyright:
* Copyright (C) 2015 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* SUBSCH  track coefficients calculation subroutine.
* K       total number of volumes for which specific values
*         of the neutron flux and reactions rates are required.
* KPN     total number of unknowns in vectors F.
* M       number of material mixtures.
* N       number of elements in the current track.
* H       vector containing the lenght of the different segments of this
*         track.
* NOM     vector containing the region number of the different segments
*         of this track.
* NZON    index-number of the mixture type assigned to each volume.
* WEIGHT  track weight.
* XST     macroscopic total cross section.
* SP      total source vector for + direction.
* SM      total source vector for - direction.
* DSP     linear component of the total source vector for + direction.
* DSM     linear component of the total source vector for - direction.
* NREG    number of volumes.
* NMU     order of the polar quadrature set.
* NANI    scattering anisotropy (=1 for isotropic scattering).
* NFUNLX  number of moments of the spherical harmonics.
* NMOD    third dimension of TRHAR.
* TRHAR   spherical harmonics components for this angle in the plane.
* KEYCUR  position of current elements in PHI vector.
* IMU     azimuthal angle corresponding to this track.
* MODP    + direction angular mode index.
* MODM    - direction angular mode index.
*
*Parameters: input/output
* F       vector containing the zonal scalar flux (surface components).
* PHIV    vector containing the zonal scalar flux (component 1).
* DPHIV   vector containing the zonal scalar flux (components 2 and 3).
*
*Parameters: scratch
* B       undefined.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER K,KPN,M,N,NOM(N),NZON(K),NMU,NFUNLX,NMOD,NREG,
     1 KEYCUR(K-NREG),IMU,NANI,MODP,MODM
      REAL XST(0:M),TRHAR(NMU,NFUNLX,NMOD)
      DOUBLE PRECISION WEIGHT,H(N),SP(N),SM(N),DSP(N),DSM(N),B(0:5,N),
     1 F(KPN),PHIV(NFUNLX,NREG),DPHIV(2*NFUNLX,NREG)
      EXTERNAL SUBSCH
*----
* LOCAL VARIABLES
*---
      REAL ETA,XI
      DOUBLE PRECISION F0,DF0,SI,SJ,DSI,DSJ,RM,RP,DSIG,DH
      INTEGER I,NOMI,IND1,INDN,JF,J,NOMJ
*----
*     Calculation of the coefficients for this track.
*----
*       MCGSCAL: Linear discontinuous-Characteristics Scheme with
*                Tabulated Exponentials
*       MCGDDFL: Diamond-Differencing DD1 Scheme
*       MCGSCEL: Linear discontinuous-Characteristics Scheme with
*                Exact Exponentials
      IF(NANI.LE.0) CALL XABORT('MCGFFAL: INVALID VALUE OF NANI.')
      CALL SUBSCH(N,K,M,NOM,NZON,H,XST,B)
*----
*     Summation along the track in both directions
*----
*     incoming flux in + direction
      IND1=KEYCUR(-NOM(1))
      RP=SP(1)
*     incoming flux in - direction
      INDN=KEYCUR(-NOM(N))
      RM=SM(N)
*     track angles in 3D
      ETA=TRHAR(IMU,3,MODP)
      XI=TRHAR(IMU,2,MODP)
      DO I=2,N-1
*        + direction
         NOMI=NOM(I)
         DSIG=XST(NZON(NOMI))
         DH=H(I)
         SI=SP(I)
         DSI=DSP(I)
         F0=B(1,I)*RP+B(2,I)*SI+B(3,I)*DSI
         DF0=B(4,I)*RP-B(3,I)*SI/(DH*DH)+B(5,I)*DSIG*DSI
         RP=B(0,I)*RP+B(1,I)*SI-B(3,I)*DSIG*DSI
         DO JF=1,NFUNLX
            PHIV(JF,NOMI)=PHIV(JF,NOMI)+WEIGHT*F0*TRHAR(IMU,JF,MODP)
            DPHIV(JF,NOMI)=DPHIV(JF,NOMI)+WEIGHT*DF0*ETA*
     >                     TRHAR(IMU,JF,MODP)
            DPHIV(NFUNLX+JF,NOMI)=DPHIV(NFUNLX+JF,NOMI)+WEIGHT*DF0*XI*
     >                            TRHAR(IMU,JF,MODP)
         ENDDO
*        - direction
         J=N+1-I
         NOMJ=NOM(J)
         DSIG=XST(NZON(NOMJ))
         DH=H(J)
         SJ=SM(J)
         DSJ=DSM(J)
         F0=B(1,J)*RM+B(2,J)*SJ+B(3,J)*DSJ
         DF0=B(4,J)*RM-B(3,J)*SJ/(DH*DH)+B(5,J)*DSIG*DSJ
         RM=B(0,J)*RM+B(1,J)*SJ-B(3,J)*DSIG*DSJ
         DO JF=1,NFUNLX
            PHIV(JF,NOMJ)=PHIV(JF,NOMJ)+WEIGHT*F0*TRHAR(IMU,JF,MODM)
            DPHIV(JF,NOMJ)=DPHIV(JF,NOMJ)-WEIGHT*DF0*ETA*
     >                     TRHAR(IMU,JF,MODM)
            DPHIV(NFUNLX+JF,NOMJ)=DPHIV(NFUNLX+JF,NOMJ)-WEIGHT*DF0*XI*
     >                            TRHAR(IMU,JF,MODM)
         ENDDO
      ENDDO
*     outgoing flux in + direction
      F(INDN)=F(INDN)+WEIGHT*RP
*     outgoing flux in - direction
      F(IND1)=F(IND1)+WEIGHT*RM
*
      RETURN
      END
