*DECK MOCFFAL
      SUBROUTINE MOCFFAL(SUBSCH,NR,NS,MT,LINE,SEGLEN,NRSEG,NE,NFX,
     1 MATALB,DWEIG,SIGANG,YG,FLM,FLP,DFLM,DFLP,NPHI,NSUB,KANGL,TRHAR,
     2 PHIV,DPHIV,DSIG,EXPT,EXP2,CYM1,CYP1,CYM2,CYP2)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solution of transport equation on a cyclic track.
* Linear-discontinuous-characteristics approximation.
* Ray-tracing (anisotropic scattering case, 'source term isolation' 
* off).
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
* NR      number of volumes.
* NS      number of surfaces.
* MT      number of material mixtures.
* LINE    number of segments on this tracking line.
* SEGLEN  vector containing the lenght of the different segments of this
*         track.
* NRSEG   vector containing the region number of the different segments
*         of this track.
* NE      order of the polar quadrature set.
* NFX     number of moments of the spherical harmonics.
* MATALB  index-number of the mixture assigned to each volume
*         and the albedo to each surface.
* DWEIG   track weight per polar angle.
* SIGANG  total cross-sections and albedos.
* YG      inverse of polar quadrature cosines.
* FLM     total source vector for + direction.
* FLP     total source vector for - direction.
* DFLM    linear component of the total source vector for + direction.
* DFLP    linear component of the total source vector for - direction.
* NPHI    number of angles in the plane.
* NSUB    number of subtracks.
* KANGL   angle indices per subtrack.
* TRHAR   spherical harmonics components for this angle in the plan.
*
*Parameters: output
* PHIV    vector containing the zonal scalar flux (component 1).
* DPHIV   vector containing the zonal scalar flux (components 2 and 3).
*
*Parameters: scratch
* DSIG    undefined.
* EXPT    track coefficient.
* EXP2    quadratic expansion of (1-exp(-a*L))/L with small argument.
* CYM1    undefined.
* CYP1    undefined.
* CYM2    undefined.
* CYP2    undefined.
*
*-----------------------------------------------------------------------
*
      IMPLICIT  NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NR,NS,MT,LINE,NRSEG(LINE),NE,NFX,MATALB(-NS:NR),NPHI,
     1 NSUB,KANGL(NSUB)
      REAL SIGANG(-6:MT),YG(NE),TRHAR(NE,NFX,NPHI,2)
      DOUBLE PRECISION SEGLEN(LINE),DWEIG(NE),FLM(NE,LINE),FLP(NE,LINE),
     1 DFLM(NE,LINE),DFLP(NE,LINE),PHIV(NFX,NR),DPHIV(2*NFX,NR),
     2 DSIG(LINE),EXPT(NE,LINE),EXP2(5,NE,LINE),CYM1(NE,LINE),
     3 CYP1(NE,LINE),CYM2(NE,LINE),CYP2(NE,LINE)
      EXTERNAL SUBSCH
*----
*  LOCAL VARIABLES
*----
      INTEGER MXE
      PARAMETER (MXE=64)
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.D0,ZERO=0.D0)
      INTEGER IE,NOIL,IL,JL,JF,ISUB,IANG
      REAL ETA,XI
      DOUBLE PRECISION TEMP1,TEMP2,TEMP3,FL0(MXE),CY0(MXE),FL1(MXE),
     1 CY1(MXE),DH
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
*       MOCSCAL: Linear discontinuous-Characteristics Scheme with
*                Tabulated Exponentials
*       MOCDDFL: Diamond-Differencing DD1 Scheme
*       MOCSCEL: Linear discontinuous-Characteristics Scheme with
*                Exact Exponentials
      CALL SUBSCH(LINE,NR,NS,MT,NRSEG,MATALB,SEGLEN,SIGANG(-6),DSIG,
     1     EXPT,EXP2,NE,YG)
*----
*     Summation along the track in both directions
*----
      DO 2226 IL= 1, LINE
         JL= LINE + 1 - IL
         DO 2225 IE=1,NE
*           + direction
            NOIL       = NRSEG(IL)
            TEMP1      = FL1(IE)
            TEMP2      = FLP(IE,IL)
            TEMP3      = DFLP(IE,IL)
            DH         = SEGLEN(IL)*YG(IE)
            FL1(IE)    = TEMP1*EXPT(IE,IL) + TEMP2*EXP2(1,IE,IL)
     1                   - TEMP3*EXP2(3,IE,IL)*DSIG(IL)
            FLP(IE,IL) = TEMP1*EXP2(1,IE,IL) + TEMP2*EXP2(2,IE,IL)
     1                   + TEMP3*EXP2(3,IE,IL)
            DFLP(IE,IL)= TEMP1*EXP2(4,IE,IL)-TEMP2*EXP2(3,IE,IL)/(DH*DH)
     1                   + TEMP3*EXP2(5,IE,IL)*DSIG(IL)
            CYP1(IE,IL)= CY1(IE) * EXP2(1,IE,IL)
            CYP2(IE,IL)= CY1(IE) * EXP2(4,IE,IL)
            CY1(IE)    = CY1(IE) * EXPT(IE,IL)
*
*           - direction
            TEMP1      = FL0(IE)
            TEMP2      = FLM(IE,JL)
            TEMP3      = DFLM(IE,JL)
            DH         = SEGLEN(JL)*YG(IE)
            FL0(IE)    = TEMP1*EXPT(IE,JL) + TEMP2*EXP2(1,IE,JL)
     1                   - TEMP3*EXP2(3,IE,JL)*DSIG(JL)
            FLM(IE,JL) = TEMP1*EXP2(1,IE,JL) + TEMP2*EXP2(2,IE,JL)
     1                   + TEMP3*EXP2(3,IE,JL)
            DFLM(IE,JL)= TEMP1*EXP2(4,IE,JL)-TEMP2*EXP2(3,IE,JL)/(DH*DH)
     1                   + TEMP3*EXP2(5,IE,JL)*DSIG(JL)
            CYM1(IE,JL)= CY0(IE) * EXP2(1,IE,JL)
            CYM2(IE,JL)= CY0(IE) * EXP2(4,IE,JL)
            CY0(IE)    = CY0(IE) * EXPT(IE,JL)
 2225    CONTINUE
 2226 CONTINUE
      DO 2230 IE=1,NE
         TEMP1=ONE-CY0(IE)
         IF (TEMP1.GT.ZERO) THEN
            FL0(IE)= FL0(IE)/TEMP1
         ELSE
            FL0(IE)= ZERO
         ENDIF
         TEMP1=ONE-CY1(IE)
         IF (TEMP1.GT.ZERO) THEN
            FL1(IE)= FL1(IE)/TEMP1
         ELSE
            FL1(IE)= ZERO
         ENDIF
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
            ETA=TRHAR(IE,3,IANG,1)
            XI=TRHAR(IE,2,IANG,1)
            DO 2241 JF=1,NFX
               TEMP1=TRHAR(IE,JF,IANG,1)
               TEMP2=TRHAR(IE,JF,IANG,2)
               PHIV(JF,NOIL)=PHIV(JF,NOIL)+DWEIG(IE)*
     >         ((FLM(IE,IL)+FL0(IE)*CYM1(IE,IL))*TEMP2+
     >          (FLP(IE,IL)+FL1(IE)*CYP1(IE,IL))*TEMP1)
               DPHIV(JF,NOIL)=DPHIV(JF,NOIL)+DWEIG(IE)*ETA*
     >         (-(DFLM(IE,IL)+FL0(IE)*CYM2(IE,IL))*TEMP2+
     >          (DFLP(IE,IL)+FL1(IE)*CYP2(IE,IL))*TEMP1)
               DPHIV(NFX+JF,NOIL)=DPHIV(NFX+JF,NOIL)+DWEIG(IE)*XI*
     >         (-(DFLM(IE,IL)+FL0(IE)*CYM2(IE,IL))*TEMP2+
     >          (DFLP(IE,IL)+FL1(IE)*CYP2(IE,IL))*TEMP1)
 2241       CONTINUE
 2242       CONTINUE
         ENDIF
 2240 CONTINUE
      IF(ISUB.NE.NSUB) CALL XABORT('MOCFFAL: NSUB INCONSISTENCY')
*
      RETURN
      END
