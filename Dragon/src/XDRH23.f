*DECK XDRH23
      SUBROUTINE XDRH23 (IBIHET,NUN,IR1,NMILG,NREG,NREG2,NG,NSMAX,
     1 KEYFLX,NS,IDIL,MIXGR,IBI,FRACT,SIGMA,P1I,P1DI,P1KI,SIGA1,SUNKNO,
     2 FLUAS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of the macro-source for the She-Liu-Shi double
* heterogeneity model (part 2).
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
* SIGMA   equivalent total macroscopic cross section in each mixture.
* P1I     non collision probability in subvolume with 1 grain type.    
* P1KI    escape probability from layer k in subvolume
*         with 1 grain type.
* P1DI    escape probability from matrix in subvolume
*         with 1 grain type.
* SIGA1   equiv. total macro. XS with 1 grain type.
* SUNKNO  sources defined in the composite geometry.
*
*Parameters: output
* FLUAS   equivalent macro-source.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IBIHET,NUN,IR1,NMILG,NREG,NREG2,NG,NSMAX,KEYFLX(NREG),
     1 NS(NG),IDIL(NMILG),MIXGR(NSMAX,NG,NMILG),IBI(NREG2)
      REAL FRACT(NG,IR1+NMILG),SIGMA(0:IR1+NMILG),P1I(NG,NMILG),
     1 P1DI(NG,NMILG),P1KI(NSMAX,NG,NMILG),SIGA1(NG,NMILG),
     2 SUNKNO(NUN),FLUAS(NREG2)
*----
*  LOCAL VARIABLES
*----
      INTEGER IND1,I,MIL,K,J,IBM,IUNK
      REAL Q1I,DILF,FRT
*----
*  SCRATCH STORAGE ALLOCATION
*----
      IF(IBIHET.NE.3) THEN
         CALL XABORT('XDRH23: INVALID DOUBLE HETEROGENEITY METHOD.')
      ENDIF
*----
*  COMPUTE THE EQUIVALENT MACRO-SOURCE (SLS METHOD).
*----
      IND1=NREG2
      DO 60 I=1,NREG2
      MIL=IBI(I)
      IF(MIL.GT.IR1) THEN
        FLUAS(I)=0.0
        IBM=MIL-IR1
        DILF=1.0
        Q1I=0.0
        DO 30 J=1,NG
          FRT=FRACT(J,MIL)
          DILF=DILF-FRT
   30     CONTINUE
        DO 50 J=1,NG
          FRT=FRACT(J,MIL)
          IF(FRT.LE.0.00001) GO TO 40
          Q1I=P1DI(J,IBM)/(1.0-P1I(J,IBM))*SUNKNO(KEYFLX(I))
     >     *SIGA1(J,IBM)/SIGMA(IDIL(IBM))
          DO K=1,NS(J)
            IUNK=KEYFLX(IND1+K)
            Q1I=Q1I+P1KI(K,J,IBM)/(1.0-P1I(J,IBM))*SUNKNO(IUNK)
     >       *SIGA1(J,IBM)/SIGMA(MIXGR(K,J,IBM))
          ENDDO
          IND1=IND1+NS(J)
   40     FLUAS(I)=FLUAS(I)+Q1I*FRT
   50   CONTINUE
        FLUAS(I)=FLUAS(I)/(1.0-DILF)
      ELSE
        FLUAS(I)=SUNKNO(KEYFLX(I))
      ENDIF
   60 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      RETURN
      END
