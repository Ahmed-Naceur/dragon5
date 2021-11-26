*DECK MOCFFIT
      SUBROUTINE MOCFFIT(SUBSCH,NR,NS,NUN,MT,LINE,SEGLEN,NRSEG,NE,
     1                   MATALB,SIGANG,KEYFLX,YG,FLUX,EXPT,EXP2,FLM,
     2                   FLP,CYM,CYP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solution of transport equation on a cyclic track
* (isotropic scattering case, 'MOCC/MCI' integration strategy).
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy and R. Le Tellier
*
*Parameters: input
* SUBSCH  track coefficients calculation subroutine.
* NR      number of volumes.
* NS      number of surfaces.
* NUN     total number of unknowns in vectors FLUX.
* MT      number of material mixtures.
* LINE    number of segments on this tracking line.
* SEGLEN  vector containing the lenght of the different segments of this
*         track.
* NRSEG   vector containing the region number of the different segments
*         of this track.
* NE      order of the polar quadrature set.
* MATALB  index-number of the mixture assigned to each volume
*         and the albedo to each surface.
* SIGANG  total cross-sections and albedos.
* KEYFLX  position of flux elements in FLUX vector.
* YG      inverse of polar quadrature cosines.
* FLM     total source vector for + direction.
* FLP     total source vector for - direction.
*
*Parameters: output
* FLUX    vector containing the zonal flux moments.
*
*Parameters: scratch
* EXPT    track coefficient.
* EXP2    quadratic expansion of (1-exp(-a*L))/L with small argument.
* CYM     undefined.
* CYP     undefined.
*
*-----------------------------------------------------------------------
*
      IMPLICIT  NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NR,NS,NUN,MT,LINE,NRSEG(LINE),NE,MATALB(-NS:NR),KEYFLX(NR)
      REAL SIGANG(-6:MT),YG(NE)
      DOUBLE PRECISION SEGLEN(LINE),FLUX(NUN),EXPT(NE,LINE),
     1 EXP2(NE,LINE),FLM(NE,LINE),FLP(NE,LINE),CYM(NE,LINE),
     2 CYP(NE,LINE)
      EXTERNAL SUBSCH
*----
*  LOCAL VARIABLES
*----
      INTEGER MXE
      PARAMETER (MXE=64)
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.D0,ZERO=0.D0)
      INTEGER IE,NOIL,IL,JL,IND
      DOUBLE PRECISION FL0(MXE),CY0(MXE),FL1(MXE),CY1(MXE)
*
      DO 4210 IE=1,NE
         FL0(IE)= ZERO
         FL1(IE)= ZERO
         CY0(IE)= ONE
         CY1(IE)= ONE
 4210 CONTINUE
*----
*     Calculation of the coefficients for this track.
*----
*       MOCSCAT: Step-Characteristics Scheme with Tabulated Exponentials
*       MOCDDFT: Diamond-Differencing Scheme
*       MOCSCET: Step-Characteristics Scheme with Exact Exponentials
      CALL SUBSCH(LINE,NR,NS,MT,NRSEG,MATALB,SEGLEN,SIGANG(-6),EXPT,
     1     EXP2,NE,YG)
*----
*     Summation along the track in both directions
*----
      DO 4226 IL= 1, LINE
         JL= LINE + 1 - IL
         DO 4225 IE=1,NE
            FLM(IE,IL) = (FL0(IE) - FLM(IE,IL)) * EXP2(IE,IL)
            CYM(IE,IL) =               CY0(IE)  * EXP2(IE,IL)
            FLP(IE,JL) = (FL1(IE) - FLP(IE,JL)) * EXP2(IE,JL)
            CYP(IE,JL) =               CY1(IE)  * EXP2(IE,JL)
            FL0(IE)    = FL0(IE) - FLM(IE,IL)
            CY0(IE)    = CY0(IE) - CYM(IE,IL)
            FL1(IE)    = FL1(IE) - FLP(IE,JL)
            CY1(IE)    = CY1(IE) - CYP(IE,JL)
 4225    CONTINUE
 4226 CONTINUE
      DO 4230 IE=1,NE
         FL0(IE)=  FL0(IE)/(ONE-CY0(IE))
         FL1(IE)=  FL1(IE)/(ONE-CY1(IE))
 4230 CONTINUE
      DO 4240 IL= 1, LINE
         NOIL  = NRSEG(IL)
         IF( NOIL.GT.0 )THEN
            IND=KEYFLX(NOIL)
            DO 4241 IE=1,NE
               FLUX(IND)= FLUX(IND)
     >              + ((FL0(IE)*CYM(IE,IL)+FLM(IE,IL))
     >              +(FL1(IE)*CYP(IE,IL)+FLP(IE,IL)))
 4241       CONTINUE
         ENDIF
 4240 CONTINUE
*
      RETURN
      END
