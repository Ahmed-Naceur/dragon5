*DECK XDRH20
      SUBROUTINE XDRH20 (IBIHET,NUN,IR1,NMILG,NREG,NREG2,NG,NSMAX,
     1 KEYFLX,NS,IDIL,MIXGR,IBI,FRACT,VOLK,SIGMA,NCO,RRRR,QKOLD,QKDEL,
     2 PKL,COEF,SUNKNO,FLUAS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of the macro-source for the Hebert or Sanchez-Pomraning
* double heterogeneity model (part 2).
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
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
* NCO     number of volumes in each composite mixture.
* QKOLD   information computed by XDRH11.
* QKDEL   information computed by XDRH11 or XDRH12.
* PKL     information computed by XDRH11 or XDRH12.
* RRRR    information computed by XDRH11 or XDRH12.
* COEF    information computed by XDRH11 or XDRH12.
* SUNKNO  sources defined in the composite geometry.
*
*Parameters: output
* FLUAS   equivalent macro-source.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IBIHET,NUN,IR1,NMILG,NREG,NREG2,NG,NSMAX,KEYFLX(NREG),
     1 NS(NG),IDIL(NMILG),MIXGR(NSMAX,NG,NMILG),IBI(NREG2),NCO(NMILG)
      REAL FRACT(NG,IR1+NMILG),VOLK(NG,NSMAX),SIGMA(0:IR1+NMILG),
     1 RRRR(NMILG),QKOLD(NG,NSMAX,NMILG),QKDEL(NG,NSMAX,NMILG),
     2 PKL(NG,NSMAX,NSMAX,NMILG),SUNKNO(NUN),FLUAS(NREG2)
      DOUBLE PRECISION COEF(1+NG*NSMAX,1+NG*NSMAX,NMILG)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION DP0,DDOT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: RHS
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(RHS(1+NG*NSMAX))
*
      IF(IBIHET.EQ.1) THEN
         GO TO 10
      ELSE IF(IBIHET.EQ.2) THEN
         GO TO 60
      ELSE
         CALL XABORT('XDRH20: INVALID DOUBLE HETEROGENEITY METHOD.')
      ENDIF
*----
*  COMPUTE THE EQUIVALENT MACRO-SOURCE (SANCHEZ-POMRANING METHOD).
*----
   10 IND1=NREG2
      DO 50 I=1,NREG2
      MIL=IBI(I)
      IF(MIL.GT.IR1) THEN
         IBM=MIL-IR1
         INDGAR=IND1
         DILF=1.0
         DP0=0.0D0
         DO 30 J=1,NG
         FRT=FRACT(J,MIL)
         DILF=DILF-FRT
         IF(FRT.LE.0.00001) GO TO 30
         DO 20 K=1,NS(J)
         DP0=DP0+FRT*VOLK(J,K)*QKOLD(J,K,IBM)*SIGMA(MIXGR(K,J,IBM))
   20    CONTINUE
         IND1=IND1+NS(J)
   30    CONTINUE
         DP0=DP0+DILF*SIGMA(IDIL(IBM))
         RHS(1)=DILF*SUNKNO(KEYFLX(I))/DP0
         IND1=INDGAR
         IND2=1
         DO 40 J=1,NG
         FRT=FRACT(J,MIL)
         IF(FRT.LE.0.00001) GO TO 40
         DO K=1,NS(J)
           IUNK=KEYFLX(IND1+K)
           RHS(1)=RHS(1)+FRT*VOLK(J,K)*QKOLD(J,K,IBM)*SUNKNO(IUNK)/DP0
           RHS(IND2+K)=0.0D0
           DO N=1,NS(J)
             IUNK=KEYFLX(IND1+N)
             RHS(IND2+K)=RHS(IND2+K)+PKL(J,K,N,IBM)*SUNKNO(IUNK)
           ENDDO
         ENDDO
         IND1=IND1+NS(J)
         IND2=IND2+NS(J)
   40    CONTINUE
         DP0=DDOT(NCO(IBM),COEF(1,1,IBM),1+NG*NSMAX,RHS,1)
         FLUAS(I)=REAL(DP0)*SIGMA(IBI(I))*RRRR(IBM)
      ELSE
         FLUAS(I)=SUNKNO(KEYFLX(I))
      ENDIF
   50 CONTINUE
      RETURN
*----
*  COMPUTE THE EQUIVALENT MACRO-SOURCE (HEBERT METHOD).
*----
   60 IND1=NREG2
      DO 100 I=1,NREG2
      MIL=IBI(I)
      IF(MIL.GT.IR1) THEN
         IBM=MIL-IR1
         INDGAR=IND1
         DILF=1.0
         DP0=0.0D0
         DO 80 J=1,NG
         FRT=FRACT(J,MIL)
         DILF=DILF-FRT
         IF(FRT.LE.0.00001) GO TO 80
         DO 70 K=1,NS(J)
         DP0=DP0+FRT*VOLK(J,K)*QKDEL(J,K,IBM)*SIGMA(MIXGR(K,J,IBM))
   70    CONTINUE
         IND1=IND1+NS(J)
   80    CONTINUE
         DP0=DP0+DILF*SIGMA(IDIL(IBM))
         RHS(1)=DILF*SUNKNO(KEYFLX(I))/DP0
         IND1=INDGAR
         IND2=1
         DO 90 J=1,NG
         FRT=FRACT(J,MIL)
         IF(FRT.LE.0.00001) GO TO 90
         DO K=1,NS(J)
           IUNK=KEYFLX(IND1+K)
           RHS(1)=RHS(1)+FRT*VOLK(J,K)*QKDEL(J,K,IBM)*SUNKNO(IUNK)/DP0
           RHS(IND2+K)=0.0D0
           DO N=1,NS(J)
             IUNK=KEYFLX(IND1+N)
             RHS(IND2+K)=RHS(IND2+K)+PKL(J,K,N,IBM)*SUNKNO(IUNK)
           ENDDO
         ENDDO
         IND1=IND1+NS(J)
         IND2=IND2+NS(J)
   90    CONTINUE
         DP0=DDOT(NCO(IBM),COEF(1,1,IBM),1+NG*NSMAX,RHS,1)
         FLUAS(I)=REAL(DP0)*SIGMA(IBI(I))
      ELSE
         FLUAS(I)=SUNKNO(KEYFLX(I))
      ENDIF
  100 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(RHS)
      RETURN
      END
