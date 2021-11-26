*DECK MCGFFIS
      SUBROUTINE MCGFFIS(SUBSCH,K,KPN,M,N,H,NOM,NZON,XST,S,NREG,KEYFLX,
     1           KEYCUR,F,B,W,OMEGA2,IDIR,NSOUT,XSI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solution of transport equation on a track
* ray-tracing (isotropic scattering case,'source term isolation' off).
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): I. Suslov and R. Le Tellier
*
*Parameters: input
* SUBSCH  Track coefficients calculation subroutine.
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
* XST     macroscopic total cross section.
* S       total source vector.
* NREG    number of volumes.
* KEYFLX  position of flux elements in PHI vector.
* KEYCUR  position of current elements in PHI vector.
* W       weight associated with this track.
* OMEGA2  square x, y and z-component of the direction
*         Omega for 2D geometry.
* IDIR    direction of fundamental current for TIBERE with MoC
*         (=0,1,2,3).
* NSOUT   number of outer surfaces.
* XSI     x,y and z component of the shape parameter for TIBERE. 
*
*Parameters: input/output
* F       vector containing the zonal scalar flux.
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
      INTEGER K,KPN,M,N,NOM(N),NZON(K),NREG,KEYFLX(NREG,1),
     1 KEYCUR(K-NREG),IDIR
      REAL XST(0:M)
      DOUBLE PRECISION W,H(N),S(KPN),F(KPN),B(2,N),OMEGA2(3)
      INTEGER NSOUT
      DOUBLE PRECISION XSI(NSOUT)
      EXTERNAL SUBSCH
*---
* LOCAL VARIABLES
*---
      DOUBLE PRECISION F0,SI,SJ,RM,RP,WW
      INTEGER I,J,NOMI,IND,IND1,INDN,NOMJ,INDC
*
      WW=DBLE(W)
*----
*     Calculation coefficients for this track.
*----
*          MCGSCAS: Step-Characteristics Scheme with Tabulated Exponentials
*          MCGDDFS: Diamond-Differencing Scheme
*          MCGSCES: Step-Characteristics Scheme with Exact Exponentials
      CALL SUBSCH(N,K,M,NOM,NZON,H,XST,B)
*----
*     Summation along the track in both directions
*----
*     Calculation of the incoming current on surface with the correction XSI
*     incoming flux in + direction
      IND1=KEYCUR(-NOM(1))
      RP=S(IND1)
*     incoming flux in - direction
      INDN=KEYCUR(-NOM(N))
      RM=S(INDN)
      IF(IDIR.GT.0) THEN
*     incoming flux in + direction
        RP=RP*OMEGA2(IDIR)/XSI(IND1-NREG)
*     incoming flux in - direction
        RM=RM*OMEGA2(IDIR)/XSI(INDN-NREG)
      ENDIF
*                        
      DO I=2,N-1
         NOMI=NOM(I)
         SI=S(NOMI)
         J=N+1-I
         NOMJ=NOM(J)
         SJ=S(NOMJ)
         IF(IDIR.EQ.0) THEN
*        + direction
           F0=B(1,I)*RP+B(2,I)*SI
           RP=RP+B(1,I)*(SI-XST(NZON(NOMI))*RP)
           IND=KEYFLX(NOMI,1)
           F(IND)=F(IND)+F0*WW
*        - direction
           F0=B(1,J)*RM+B(2,J)*SJ
           RM=RM+B(1,J)*(SJ-XST(NZON(NOMJ))*RM)
           IND=KEYFLX(NOMJ,1)
           F(IND)=F(IND)+F0*WW
         ELSE
*        Solve current equation for TIBERE
*        and calculate Xi, Yi and Zi 
*        + direction
           F0=B(1,I)*RP+B(2,I)*SI*OMEGA2(IDIR)
           RP=RP+B(1,I)*(SI*OMEGA2(IDIR)-XST(NZON(NOMI))*RP)
           IND=KEYFLX(NOMI,1)
           F(IND)=F(IND)+F0*WW
           INDC=KPN/2+IND
           F(INDC)=F(INDC)+F0*WW*OMEGA2(IDIR)
*        - direction
           F0=B(1,J)*RM+B(2,J)*SJ*OMEGA2(IDIR)
           RM=RM+B(1,J)*(SJ*OMEGA2(IDIR)-XST(NZON(NOMJ))*RM)
           IND=KEYFLX(NOMJ,1)
           F(IND)=F(IND)+F0*WW
           INDC=KPN/2+IND
           F(INDC)=F(INDC)+F0*WW*OMEGA2(IDIR)
         ENDIF
      ENDDO
*     outgoing flux in + direction
       F(INDN)=F(INDN)+RP*WW
*     outgoing flux in - direction
       F(IND1)=F(IND1)+RM*WW    
      RETURN
      END
