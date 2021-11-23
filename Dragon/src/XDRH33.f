*DECK XDRH33
      SUBROUTINE XDRH33 (IBIHET,NUN,IR1,NMILG,NREG,NREG2,NG,NSMAX,
     1 KEYFLX,NS,IDIL,MIXGR,IBI,FRACT,VOLK,SIGMA,P1I,P1DI,P1KI,SIGA1,
     2 FUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of the composite flux for the She-Liu-Shi double
* heterogeneity model (part 3).
*
*Copyright:
* Copyright (C) 2019 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Chambon
*
*Parameters: input
* IBIHET  type of double-heterogeneity method (=1 Sanchez-Pomraning
*         model; =2 Hebert model).
* NUN     number of unknown in the system.
* IR1     number of mixtures in the domain.
* NMILG   number of composite mixtures in the domain.
* NREG    number of volumes in the composite geometry.
* NREG2   number of volumes in the macro geometry.
* NG      number of different kind of micro structures. A kind of
*         micro structure is characterized by the radius of its
*         micro volumes. All the micro volumes of the same kind
*         should own the same nuclear properties in a given macro
*         volume.
* NSMAX   maximum number of volumes (tubes or shells) in each kind of
*         micro structure.
* KEYFLX  flux elements in unknown system.
* NS      number of volumes in each kind of micro structure.
* IDIL    elementary mixture indices in the diluent of the composite
*         mixtures.
* MIXGR   elementary mixture indices in the micro structures.
* IBI     type of mixture in each volume of the macro geometry.
*         If IBI(IKK) is greater than IR1, the volume IKK contains a
*         micro structure.
* FRACT   volumic fractions of the micro volumes.
* VOLK    volumic fractions of the tubes or shells in the micro volumes.
* SIGMA   equivalent total macroscopic cross section in each mixture.
* P1I     non collision probability in subvolume with 1 grain type.    
* P1KI    escape probability from layer k in subvolume
*         with 1 grain type.
* P1DI    escape probability from matrix in subvolume
*         with 1 grain type.
* SIGA1   corrected cross section.
*
*Parameters: input/output
* FUNKNO  macro-flux on input (solution of the transport equation 
*         defined over the macro-geometry) and composite flux on output.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IBIHET,NUN,IR1,NMILG,NREG,NREG2,NG,NSMAX,KEYFLX(NREG),
     1 NS(NG),IDIL(NMILG),MIXGR(NSMAX,NG,NMILG),IBI(NREG2)
      REAL FRACT(NG,IR1+NMILG),VOLK(NG,NSMAX),SIGMA(0:IR1+NMILG),
     1 P1I(NG,NMILG),P1DI(NG,NMILG),P1KI(NSMAX,NG,NMILG),
     2 SIGA1(NG,NMILG),FUNKNO(NUN)
*----
*  LOCAL VARIABLES
*----
      INTEGER IBM,I,J,K,IUNK,IND1,MIL
      REAL FAVG,SIGTT,DILF,FRT
*----
*  SCRATCH STORAGE ALLOCATION
*----
      IF(IBIHET.NE.3) THEN
         CALL XABORT('XDRH33: INVALID DOUBLE HETEROGENEITY METHOD.')
      ENDIF
*----
*  COMPUTE THE COMPOSITE FLUX.
*----
      IND1=NREG2
      DO 60 I=1,NREG2
      MIL=IBI(I)
      IF(MIL.GT.IR1) THEN
         IBM=MIL-IR1
         FAVG=FUNKNO(KEYFLX(I))
         SIGTT=SIGMA(IR1+IBM)
         DILF=1.0
         DO 30 J=1,NG
           FRT=FRACT(J,MIL)
           DILF=DILF-FRT
           IF(FRT.LE.0.00001) GO TO 30
   30      CONTINUE
* flux in composite matrix         
         FUNKNO(KEYFLX(I))=0.0
         IUNK=KEYFLX(I)
         DO 40 J=1,NG
           FRT=FRACT(J,MIL)
           IF(FRT.LE.0.00001) GO TO 40
           FUNKNO(KEYFLX(I))=FUNKNO(KEYFLX(I))+FAVG*SIGA1(J,IBM)
     >      /SIGMA(IDIL(IBM))/(1.0-P1I(J,IBM))*P1DI(J,IBM)/DILF
           DO K=1,NS(J)
             IUNK=KEYFLX(IND1+K)
             FUNKNO(IUNK)=FAVG*SIGA1(J,IBM)/SIGMA(MIXGR(K,J,IBM))
     >        /(1.0-P1I(J,IBM))*P1KI(K,J,IBM)/(1.0-DILF)/VOLK(J,K)
           ENDDO
           IND1=IND1+NS(J)
   40      CONTINUE
      ENDIF
   60 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      RETURN
      END
