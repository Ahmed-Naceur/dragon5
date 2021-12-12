*DECK MOCFFAT
      SUBROUTINE MOCFFAT(SUBSCH,NR,NS,NUN,MT,LINE,SEGLEN,NRSEG,NE,NF,
     1                   MATALB,SIGANG,KEYFLX,YG,FLUX,EXPT,EXP2,FLM,FLP,
     2                   CYM,CYP,NPHI,NSUB,KANGL,TRHAR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solution of transport equation on a cyclic track.
* Ray-tracing (anisotropic scattering case, 'MOCC/MCI' integration
* strategy).
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
* SEGLEN  vector containing the lenght of the different segments of
*         this track.
* NRSEG   vector containing the region number of the different segments
*         of this track.
* NE      order of the polar quadrature set.
* NF      number of moments of the flux (in 2D NFUNL=NLF*(NLF+1)/2).
* MATALB  index-number of the mixture assigned to each volume
*         and the albedo to each surface.
* SIGANG  total cross-sections and albedos.
* KEYFLX  position of flux elements in FLUX vector.
* YG      inverse of polar quadrature cosines.
* FLM     total source vector for + direction.
* FLP     total source vector for - direction.
* NPHI    number of angles in the plane.
* NSUB    number of subtracks.
* KANGL   angle indices per subtrack.
* TRHAR   spherical harmonics components for this angle in the plan.
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
      INTEGER NR,NS,NUN,MT,LINE,NRSEG(LINE),NE,NF,MATALB(-NS:NR),
     1 KEYFLX(NR,NF),NPHI,NSUB,KANGL(NSUB)
      REAL SIGANG(-6:MT),YG(NE),TRHAR(NE,NF,NPHI,2)
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
      INTEGER IE,NOIL,IL,JL,JF,IND,ISUB,IANG
      DOUBLE PRECISION FL0(MXE),CY0(MXE),FL1(MXE),CY1(MXE)
      LOGICAL LNEW
*
      DO 2210 IE=1,NE
         FL0(IE)= ZERO
         FL1(IE)= ZERO
         CY0(IE)= ONE
         CY1(IE)= ONE
 2210 CONTINUE
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
      DO 2226 IL= 1, LINE
         JL= LINE + 1 - IL
         DO 2225 IE=1,NE
*           + direction
            FLP(IE,IL) = (FL1(IE) - FLP(IE,IL)) * EXP2(IE,IL)
            CYP(IE,IL) =               CY1(IE)  * EXP2(IE,IL)
            FL1(IE)    = FL1(IE) - FLP(IE,IL)
            CY1(IE)    = CY1(IE) - CYP(IE,IL)
*
*           - direction
            FLM(IE,JL) = (FL0(IE) - FLM(IE,JL)) * EXP2(IE,JL)
            CYM(IE,JL) =               CY0(IE)  * EXP2(IE,JL)
            FL0(IE)    = FL0(IE) - FLM(IE,JL)
            CY0(IE)    = CY0(IE) - CYM(IE,JL)
 2225    CONTINUE
 2226 CONTINUE
      DO 2230 IE=1,NE
         FL0(IE)=  FL0(IE)/(ONE-CY0(IE))
         FL1(IE)=  FL1(IE)/(ONE-CY1(IE))
 2230 CONTINUE
      ISUB=0
      IANG=0
      LNEW=.TRUE.
      DO 2240 IL= 1, LINE
         NOIL  = NRSEG(IL)
         IF(NOIL.LE.0) THEN
            LNEW=.TRUE.
         ELSE
            IF(LNEW) THEN
               ISUB=ISUB+1
               IANG=KANGL(ISUB)
               LNEW=.FALSE.
            ENDIF
            DO 2242 IE=1,NE
            DO 2241 JF=1,NF
               IND=KEYFLX(NOIL,JF)
               FLUX(IND)= FLUX(IND)
     >              +(FLM(IE,IL)+FL0(IE)*CYM(IE,IL))*TRHAR(IE,JF,IANG,2)
     >              +(FLP(IE,IL)+FL1(IE)*CYP(IE,IL))*TRHAR(IE,JF,IANG,1)
 2241       CONTINUE
 2242       CONTINUE
         ENDIF
 2240 CONTINUE
*
      RETURN
      END
