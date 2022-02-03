*DECK MOCDS2
      SUBROUTINE MOCDS2(SUBDSC,LC,M,N,H,NOM,NZON,TR,SC,W,NFI,DIAGF,
     1                  DIAGQ,CA,CQ,PREV,NEXT,DINV2,A2,B2)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of contribution in second-order ACA coefficients on one
* track.
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
* SUBDSC  ACA coefficients calculation subroutine.
* LC      dimension of vector MCU.
* M       number of material mixtures.
* N       number of elements for this track.
* H       tracking widths.
* NOM     integer tracking elements.
* NZON    index-number of the mixture type assigned to each volume.
* TR      macroscopic total cross section.
* SC      macroscopic P0 scattering cross section.
* W       weight associated with this track.
* NFI     total number of volumes for which specific values
*         of the neutron flux and reactions rates are required.
*
*Parameters: input/output
* CA      undefined.
* CQ      undefined.
* DIAGQ   undefined.
* DIAGF   undefined.
*
*Parameters: scratch
* PREV    undefined.
* NEXT    undefined.
* DINV2   undefined.
* A2      undefined.
* B2      undefined.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*
*---
* SUBROUTINE ARGUMENTS
*---
      INTEGER LC,M,N,NFI,NZON(NFI),NOM(N),PREV(N),NEXT(N)
      DOUBLE PRECISION W,H(N),CA(LC),DIAGF(NFI),DINV2(N),A2(N),B2(N)
      REAL TR(0:M),SC(0:M),DIAGQ(NFI),CQ(LC)
      EXTERNAL SUBDSC
*---
* LOCAL VARIABLES
*---
      INTEGER IBCV
      DOUBLE PRECISION DMINV,DMAX
      PARAMETER(DMINV=2.D-2,IBCV=-7)
      DOUBLE PRECISION WW,AAW,CAW,AQW,CQW,CAWP,CQWP,B,A,DINV,DH,DHP,
     1 BN,AN,DINVN,BP,AP,DINVP
      INTEGER I,I1,NOMI,NZI,NOMIN,NZIN,NOMIP,NZIP,ICN,ICP
*
      DMAX=1.D0/DMINV
      WW=W
*---
*  CALCULATE COEFFICIENTS OF THIS TRACK
*---
*          MCGDS2A: Tabulated Exponentials
*          MCGDS2E: Exact Exponentials
      CALL SUBDSC(N,M,NFI,NOM,NZON,H,TR,SC,DINV2,B2,A2)
*----
*  CONSTRUCTION OF ACA MATRICES
*---
      NOMIN=0
      NZIN=0
      NOMI=0
      NZI=0
      DHP=0.0D0
      DH=0.0D0
      AN=0.0D0
      BN=0.0D0
      AP=0.0D0
      BP=0.0D0
      DO I=1,N
         ICN=NEXT(I)
         ICP=PREV(I)
         IF (I.GT.1) THEN
            NOMIP=NOMI
            NZIP=NZI
            NOMI=NOMIN
            NZI=NZIN
         ELSE
            NOMIP=NOM(N)
            NZIP=NZON(NOMIP)
            NOMI=NOM(1)
            NZI=NZON(NOMI)
         ENDIF
         IF (I.LT.N) THEN
            NOMIN=NOM(I+1)
            NZIN=NZON(NOMIN)
         ELSE
            NOMIN=NOM(1)
            NZIN=NZON(NOMIN)
         ENDIF
         IF (NZI.GE.0) THEN
*        -----------
*        Volume Cell
*        -----------
            DINV=DINV2(I)
            B=B2(I)
            A=A2(I)
            IF (NZIN.GE.0) THEN
*           next cell is a volume
               IF (I.EQ.N) THEN
                  I1=1
               ELSE
                  I1=I+1
               ENDIF
               DINVN=DINV2(I1)
               BN=B2(I1)
               AN=A2(I1)
               IF (ABS(DINV+DINVN).LT.DMINV) THEN
                  DH=DMAX
               ELSE
                  DH=1.D0/(DINV+DINVN)
               ENDIF
            ELSEIF (NZIN.EQ.IBCV) THEN
*           next cell is a fixed boundary condition
               DH=1.D0/(1.D0+DINV)
            ENDIF
            IF (NZIP.GE.0) THEN
*           previous cell is a volume
               IF (I.EQ.1) THEN
                  I1=N
               ELSE
                  I1=I-1
               ENDIF
               DINVP=DINV2(I1)
               BP=B2(I1)
               AP=A2(I1)               
               IF (ABS(DINV+DINVP).LT.DMINV) THEN
                  DHP=DMAX
               ELSE
                  DHP=1.D0/(DINV+DINVP)
               ENDIF
            ELSEIF (NZIP.EQ.IBCV) THEN
*           previous cell is a fixed boundary condition
               DHP=1.D0/(1.D0+DINV)
            ENDIF
*
*           assembling coefficients
            AAW=DH*A
            AQW=DH*B
            IF(NZIN.GE.0) THEN
*           next cell is a volume
               CAW=DH*AN
               CQW=DH*BN
            ELSE
*           next cell is a voided boundary condition
               CAW=0.D0
               CQW=0.D0
            ENDIF
*
            AAW=AAW+DHP*A
            AQW=AQW+DHP*B
            IF(NZIP.GE.0) THEN
*           previous cell is a volume
               CAWP=DHP*AP
               CQWP=DHP*BP
            ELSE
*           previous cell is a voided boundary condition
               CAWP=0.D0
               CQWP=0.D0
            ENDIF
*
*           assembling matrices
            DIAGF(NOMI)=DIAGF(NOMI)+AAW*WW
            DIAGQ(NOMI)=DIAGQ(NOMI)-REAL(W*AQW)
            IF(ICN.GT.0) THEN
*           next cell is a volume different from this one
               CA(ICN)=CA(ICN)-CAW*WW
               CQ(ICN)=CQ(ICN)+REAL(W*CQW)
            ELSE
*           next cell is a voided boundary or a volume identical to this one
                DIAGF(NOMI)=DIAGF(NOMI)-CAW*WW
                DIAGQ(NOMI)=DIAGQ(NOMI)+REAL(W*CQW)
            ENDIF 
            IF(ICP.GT.0) THEN
*           previous cell is a volume different from this one
               CA(ICP)=CA(ICP)-CAWP*WW
               CQ(ICP)=CQ(ICP)+REAL(W*CQWP)
            ELSE
*           previous cell is a voided boundary or a volume identical to this one
                DIAGF(NOMI)=DIAGF(NOMI)-CAWP*WW
                DIAGQ(NOMI)=DIAGQ(NOMI)+REAL(W*CQWP)
            ENDIF 
         ENDIF
*        -----------
      ENDDO
*
      RETURN
      END
