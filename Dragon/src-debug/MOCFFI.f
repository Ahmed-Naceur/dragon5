*DECK MOCFFI
      SUBROUTINE MOCFFI(SUBSCH,NR,NS,NUN,MT,LINE,SEGLEN,NRSEG,NE,MATALB,
     1                  SIGANG,KEYFLX,YG,FLUX,EXPT,EXP2,FLM,FLP,CYM,CYP,
     2                  IDIR,OMG2)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solution of transport equation on a cyclic track
* ray-tracing (isotropic scattering case and 'source term isolation' 
* on).
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
* IDIR    direction of fundamental current for TIBERE with MoC 
*         =0,1,2,3.
* OMG2    x, y and z components of the $3\\Omega^2$. 
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
      INTEGER IDIR
      REAL SIGANG(-6:MT),YG(NE)
      DOUBLE PRECISION SEGLEN(LINE),FLUX(NUN),EXPT(NE,LINE),
     1 EXP2(NE,LINE),FLM(NE,LINE),FLP(NE,LINE),CYM(NE,LINE),
     2 CYP(NE,LINE),OMG2(NE,3)
      EXTERNAL SUBSCH
*----
*  LOCAL VARIABLES
*----
      INTEGER MXE
      PARAMETER (MXE=64)
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.D0,ZERO=0.D0)
      INTEGER IE,NOIL,IL,JL,IND,INDC
      DOUBLE PRECISION TEMP,FL0(MXE),CY0(MXE),FL1(MXE),CY1(MXE)
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
*       MOCSCA: Step-Characteristics Scheme with Tabulated Exponentials
*       MOCDDF: Diamond-Differencing Scheme
*       MOCSCE: Step-Characteristics Scheme with Exact Exponentials
      CALL SUBSCH(LINE,NR,NS,MT,NRSEG,MATALB,SEGLEN,SIGANG(-6),EXPT,
     1     EXP2,NE,YG)  
*----
*     Summation along the track in both directions
*----        
      DO 4226 IL= 1, LINE
         JL= LINE + 1 - IL
         DO 4225 IE=1,NE
*           phi_k
            TEMP       = FL0(IE)
*           phi_{k+1}
            FL0(IE)    = TEMP * EXPT(IE,IL) 
     1                 + FLM(IE,IL) * EXP2(IE,IL)
*           phi_k * ((1 - exp(-tau_k)) / tau_k)
            FLM(IE,IL) = TEMP * EXP2(IE,IL)
*           ((1 - exp(-tau_k)) / tau_k) * exp(-tau_1^{k-1}) 
            CYM(IE,IL) = CY0(IE) * EXP2(IE,IL)
*           exp(-tau_1^{k})
            CY0(IE)    = CY0(IE) * EXPT(IE,IL)
*
            TEMP       = FL1(IE)
            FL1(IE)    = TEMP * EXPT(IE,JL) 
     1                 + FLP(IE,JL) * EXP2(IE,JL)
            FLP(IE,JL) = TEMP * EXP2(IE,JL)
            CYP(IE,JL) = CY1(IE) * EXP2(IE,JL)
            CY1(IE)    = CY1(IE) * EXPT(IE,JL)
 4225    CONTINUE
 4226 CONTINUE
      DO 4230 IE=1,NE
         TEMP=ONE-CY0(IE)
         IF (TEMP.GT.ZERO) THEN
            FL0(IE)= FL0(IE)/TEMP
         ELSE
            FL0(IE)= ZERO
         ENDIF
         TEMP=ONE-CY1(IE)
         IF (TEMP.GT.ZERO) THEN
            FL1(IE)= FL1(IE)/TEMP
         ELSE
            FL1(IE)= ZERO
         ENDIF
 4230 CONTINUE
      DO 4240 IL= 1, LINE
         NOIL  = NRSEG(IL)
         IF( NOIL.GT.0 )THEN
            IND=KEYFLX(NOIL)
            INDC=NUN/2+IND
            DO 4241 IE=1,NE
               FLUX(IND)= FLUX(IND)
     >              + ((FL0(IE)*CYM(IE,IL)+FLM(IE,IL))
     >              +(FL1(IE)*CYP(IE,IL)+FLP(IE,IL)))
*     CALCULATE XI, YI OR ZI FOR TIBERE     
               IF(IDIR.GT.0) THEN
                 FLUX(INDC)= FLUX(INDC)
     >              + ((FL0(IE)*CYM(IE,IL)+FLM(IE,IL))
     >              +(FL1(IE)*CYP(IE,IL)+FLP(IE,IL)))
     >              *OMG2(IE,IDIR)
               ENDIF
 4241       CONTINUE
         ENDIF
 4240 CONTINUE
*
      RETURN
      END
