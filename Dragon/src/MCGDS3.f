*DECK MCGDS3
      SUBROUTINE MCGDS3(NLONG,PACA,M,TR,SC,V,NZON,LC,MCU,IM,JU,LC0,IM0,
     1 MCU0,DIAGF,CF,DIAGQ,DIAGFR,CFR,LUDF,LUCF)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Adding of capture and when PACA>=2 calculation of ILU0 preconditioner 
* for BICGSTAB iterations to solve the precontioning system.
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
* NLONG   corrective system dimension.
* PACA    type of preconditioner to solve the ACA corrective system.
* M       number of material mixtures.
* TR      macroscopic total cross section.
* SC      macroscopic P0 scattering cross section.
* V       volumes.
* NZON    index-number of the mixture type assigned to each volume.
* LC      dimension of CF and MCU.
* MCU     used in CDD acceleration.
* IM      used in CDD acceleration.
* JU      used in ilu0 preconditioner.
* LC0     used in ILU0-ACA acceleration.
* IM0     used in ILU0-ACA acceleration.
* MCU0    used in ILU0-ACA acceleration.
* DIAGF   diagonal contribution for D (Double Precision).
* CF      non diagonal contribution for left hand-side matrix D.
*
*Parameters: input/output
* DIAGQ   diagonal contribution for right hand-side matrix
*         used in CDD acceleration.
*
*Parameters: output
* DIAGFR  diagonal contribution for D.
* CFR     non diagonal contribution for D.
* LUDF    diagonal contribution for ilu0 decomposition of D.
* LUCF    non diagonal contribution for ilu0 decomposition of D.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NLONG,PACA,M,NZON(NLONG),LC,MCU(LC),IM(NLONG+1),JU(NLONG),
     1 LC0,IM0(*),MCU0(*)
      DOUBLE PRECISION DIAGF(NLONG),CF(LC)
      REAL TR(0:M),SC(0:M),V(NLONG),DIAGQ(NLONG),DIAGFR(NLONG),CFR(LC),
     1 LUDF(NLONG),LUCF(LC0)
*----
*  LOCAL VARIABLES
*----
      INTEGER IBCV
      PARAMETER(IBCV=-7)
*----
*  ADDING OF CAPTURE
*----
      DO I=1,NLONG
         NZN=NZON(I)
         IF(NZN.GE.0) THEN
            DIAGF(I)=DIAGF(I)+DBLE((TR(NZN)-SC(NZN))*V(I)/2.0)
            DIAGQ(I)=DIAGQ(I)+V(I)/2.0
         ELSEIF(NZN.EQ.IBCV) THEN
            DIAGF(I)=1.D0
         ENDIF
      ENDDO
*----
*  TYPE CONVERSION
*----
      DO I=1,NLONG
         DIAGFR(I)=REAL(DIAGF(I))
      ENDDO
      DO I=1,LC
         CFR(I)=REAL(CF(I))
      ENDDO
      IF(PACA.GE.2) THEN
*----
*  ILU0 DECOMPOSITION : BICGSTAB WILL BE USED TO SOLVE THE SYSTEM
*----
         CALL MSRILU(NLONG,LC,IM,MCU,JU,DIAGF,CF)
         DO I=1,NLONG
            LUDF(I)=REAL(DIAGF(I))
         ENDDO
         IF(PACA.EQ.2) THEN
            DO I=1,LC
               LUCF(I)=REAL(CF(I))
            ENDDO
         ELSEIF(PACA.EQ.3) THEN
            DO I=1,NLONG
               DO IJ=IM0(I)+1,IM0(I+1)
                  J=MCU0(IJ)
                  DO IK=IM(I)+1,IM(I+1)
                     IF(MCU(IK).EQ.J) GOTO 10
                  ENDDO
                  CALL XABORT('MCGDS3: ILU0-ACA PROBLEM')
 10               LUCF(IJ)=REAL(CF(IK))
               ENDDO
            ENDDO
         ENDIF
      ENDIF
*
      RETURN
      END
