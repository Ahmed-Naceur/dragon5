*DECK PIJRHL
      SUBROUTINE PIJRHL(IPRT,NREG,NSOUT,SIGTAL,PROB)
*
*-----------------------------------------------------------------------
*
*Purpose:
* HELIOS type normalization of collision probs (CP).
*
*Copyright:
* Copyright (C) 1994 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy, E. Varin
*
*Parameters: input
* IPRT    print level.
* NREG    number of zones for geometry.
* NSOUT   number of surfaces for geometry.
* SIGTAL  albedo-sigt vector.
*
*Parameters: input/output
* PROB    CP matrix for all types.
*
*References:
* R. Roy and G. Marleau,
* Normalization Techniques for CP Matrices,
* CONF/PHYSOR-90, Marseille/France, V 2, P IX-40 (1990).
* \\\\
* E.A. Villarino, R.J.J. Stamm'ler, A.A. Ferri and J.J. Casal
* HELIOS: Angularly Dependent Collision Probabilities.
* Nucl.Sci.Eng. 112,16-31, 1992.
*
*-----------------------------------------------------------------------
*
      IMPLICIT   NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER   IPRT,NREG,NSOUT 
      REAL      SIGTAL(-NSOUT:NREG)
      DOUBLE PRECISION PROB(*)
*----
*  LOCAL VARIABLES
*----
      INTEGER    IUNOUT,NITMAX,NIT,IPRINT,IR,JR,IP,IPRB,IND,I,J,CPTLB,
     >           CPTAC,CTOT,NSURC,NSURM,NVOLC,NVOLM
      LOGICAL    NOTCON
      DOUBLE PRECISION NOM,DENOM,DMU,WFSPAD,WFSP,EPSCON,R1,R2,TOTCON,
     >                 TMPCON
      CHARACTER  HSMG*131
      PARAMETER (IUNOUT=6, IPRINT=10, EPSCON=1.0E-6, NITMAX=20)
*----
*  ALLOCATABLE ARRAYS
*----
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: CHI
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: WEIG
*
*----- INTRINSIC FUNCTION FOR POSITION IN CONDENSE PIJ MATRIX
*
      IND(I,J)=(MAX(I+NSOUT+1,J+NSOUT+1)*
     >         (MAX(I+NSOUT+1,J+NSOUT+1)-1))/2
     >         +MIN(I+NSOUT+1,J+NSOUT+1)
*----
*  SCRATCH STORAGE ALLOCATION
*   WEIG : ADDITIVE WEIGHT
*----
      ALLOCATE(WEIG(-NSOUT:NREG,3),CHI(-NSOUT:NREG))
*
      NOTCON= .FALSE.
      CPTLB = 3
      CPTAC = 3
      CTOT = CPTAC+CPTLB
*
*     INITIALISATION OF WEIGHTS
      DO 60 IR=-NSOUT, NREG
         WEIG(IR,1)=0.0D0
         WEIG(IR,2)=0.5D0
         WEIG(IR,3)=0.5D0
   60 CONTINUE
      DO 50 IR=-NSOUT, NREG
         CHI(IR)= 1.0D0
         IF( IR.GE.0.AND.SIGTAL(IR).EQ.0.0D0 )THEN
            CHI(IR)= 0.0D0
         ENDIF
   50 CONTINUE
*
*
*     MAIN ITERATION LOOP
      IF(IPRT.GT.2) WRITE(IUNOUT,'(A24)')
     >       'ITER.     MU      ERROR '
      DO 110 NIT=1,NITMAX
*
         DO 220 IR= -NSOUT, NREG
            WFSPAD = PROB(IND(IR,0))
     >                   + CHI(IR)*PROB(IND(IR,IR))*WEIG(IR,3)
            WFSP = CHI(IR)*PROB(IND(IR,IR))
            DO 200 JR=-NSOUT, NREG
               WFSPAD = WFSPAD - CHI(JR)*WEIG(JR,3)*PROB(IND(IR,JR))
               WFSP = WFSP + CHI(JR)*PROB(IND(IR,JR))
  200       CONTINUE
            WEIG(IR,3) = WFSPAD / WFSP
  220    CONTINUE
*
*        ACCELERATION TECHNIQUE
         IF(  MOD(NIT-1,CTOT).GE.CPTAC )THEN
            NOM   = 0.0D0
            DENOM = 0.0D0
            DO 10 IR=-NSOUT, NREG
               R1= WEIG(IR,2) - WEIG(IR,1)
               R2= WEIG(IR,3) - WEIG(IR,2)
               NOM = NOM + R1*(R2-R1)
               DENOM = DENOM + (R2-R1)*(R2-R1)
   10       CONTINUE
            IF(DENOM.EQ.0.0D0) THEN
              DMU = 1.0D0
            ELSE
              DMU = - NOM / DENOM
            ENDIF
            IF( DMU.GT.50.0D0 .OR. DMU.LT.0.0D0 ) THEN
              WRITE(HSMG,'(37HPIJRHL: PROBLEM OF ACCELERATION (DMU=,1P,
     >        E11.4,2H).)') DMU
              CALL XABORT(HSMG)
            ENDIF
            DO 20 IR=-NSOUT, NREG
               WEIG(IR,3) = WEIG(IR,2) + DMU *
     >                           (WEIG(IR,3) - WEIG(IR,2))
               WEIG(IR,2) = WEIG(IR,1) + DMU *
     >                           (WEIG(IR,2) - WEIG(IR,1))
   20       CONTINUE
         ELSE
            DMU = 1.0D0
         ENDIF
*
*        CALCULATIONS OF SQUARE DISTANCE BETWEEN 2 ITERATIONS
*        AND UPDATING THE SOLUTION
         TOTCON = 0.0D0
         DO 100 IR=-NSOUT, NREG
            TMPCON=ABS(WEIG(IR,3)-WEIG(IR,2))/WEIG(IR,3)
            TOTCON=MAX(TMPCON,TOTCON)
            WEIG(IR,1)= WEIG(IR,2)
            WEIG(IR,2)= WEIG(IR,3)
  100    CONTINUE
         IF( IPRT.GT.2 ) WRITE(IUNOUT,'(I3,F9.5,E15.7)') NIT,DMU,TOTCON
*
*        CONVERGENCE TEST
         IF( TOTCON.LT.EPSCON )GO TO 120
*
  110 CONTINUE
      NOTCON=.TRUE.
      WRITE(IUNOUT,'(35H PIJRHL: WEIGHTS NOT CONVERGED          )')
  120 CONTINUE
*
*     RENORMALIZE "PIJ" SYMMETRIC MATRIX
      IPRB = 0
      DO 240 IR   = -NSOUT, NREG
         DO 230 JR= -NSOUT, IR
            IPRB= IPRB+1
            IF( IR.NE.0.AND.JR.NE.0 )THEN
                PROB(IPRB)=PROB(IPRB)*(WEIG(IR,1)+WEIG(JR,1))
            ENDIF
  230    CONTINUE
  240 CONTINUE
*
*     PRINT WEIGHT FACTORS IF THERE IS A PROBLEM...
      IF( NOTCON .OR. IPRT.GE.IPRINT )THEN
         WRITE(IUNOUT,'(30H0 SURFACE WEIGHTS FACTORS                /)')
         NSURC = -1
         DO 300 IP  = 1, (9 +NSOUT) / 10
            NSURM= MAX( -NSOUT, NSURC-9 )
            WRITE(IUNOUT,'(10X,10( A5,    I6)/)')
     >                     (' SUR ',-IR,IR= NSURC, NSURM, -1)
            WRITE(IUNOUT,'(10H WEIGHT   ,10F11.5)')
     >                   (WEIG(IR,1),IR=NSURC,NSURM,-1)
            NSURC = NSURC - 10
 300     CONTINUE
         WRITE(IUNOUT,'(30H0  VOLUME WEIGHTS FACTORS                /)')
         NVOLC =  1
         DO 310 IP  = 1, (9 + NREG) / 10
            NVOLM= MIN( NREG, NVOLC+9 )
            WRITE(IUNOUT,'(10X,10( A5 ,  I6)/)')
     >                  (' VOL ',IR,IR=NVOLC,NVOLM, 1)
            WRITE(IUNOUT,'(10H WEIGHT   ,10F11.5)')
     >                   (WEIG(IR,1),IR=NVOLC,NVOLM, 1)
            NVOLC = NVOLC + 10
 310     CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(CHI,WEIG)
      RETURN
      END
