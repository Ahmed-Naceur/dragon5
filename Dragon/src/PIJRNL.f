*DECK PIJRNL
      SUBROUTINE PIJRNL(IPRT,NREG,NSOUT,SIGTAL,PROB)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Non-linear type normalization of collision probs (CP).
*
*Copyright:
* Copyright (C) 1994 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy, G. Marleau
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
*
*-----------------------------------------------------------------------
*
      IMPLICIT   NONE
      INTEGER    IPRT,NREG,NSOUT,IUNOUT,NITMAX,NIT,
     >           NUNKNO,IPRB,IPRF,IVOL,IDIA,IUNK,JUNK,IR,JR,
     >           NSURC,NSURM,NVOLC,NVOLM,IP,NPR
      REAL       SIGTAL(-NSOUT:NREG)
      DOUBLE PRECISION PROB(*),EPSCON,TOTCON,WFSPAD
      PARAMETER (IUNOUT=6, EPSCON=1.0E-8, NITMAX=10)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDL
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: CIJ,WSPACE,WFSP,
     > WEIG
C----
C  SCRATCH STORAGE ALLOCATION
C    CIJ   : MODIFIED CP PROB MATRIX
C    WSPACE: NON-LINEAR SYSTEM MATRIX
C    WFSP  : NON LINEAR SYSTEM SOLUTION
C    IDL   : WSPACE DIAGONAL LOCATION
C    WEIG  : NON-LINEAR WEIGHT
C----
      NPR=(NSOUT+NREG+1)*(NSOUT+NREG+2)/2
      ALLOCATE(IDL(NSOUT+NREG+1))
      ALLOCATE(CIJ(NPR),WSPACE(NPR),WFSP(-NSOUT:NREG),WEIG(-NSOUT:NREG))
C
C     CHARGE MATRIX "CIJ"
      NUNKNO=NREG+NSOUT+1
      IPRB= 0
      IUNK= 0
      IVOL= NSOUT*(NSOUT+1)/2
      DO 20 IR   = -NSOUT, NREG
         IUNK= IUNK+1
         IF( IR.LT.0.OR.SIGTAL(IR).GT.0.0 )THEN
            DO 10 JR= -NSOUT, IR-1
               IPRB= IPRB+1
               IF( JR.LT.0.OR.SIGTAL(JR).GT.0.0 )THEN
                  CIJ(IPRB)= PROB(IPRB)
               ELSE
                  CIJ(IPRB)= 0.0D0
               ENDIF
   10       CONTINUE
         ELSE
            DO 15 JR= -NSOUT, IR-1
               IPRB= IPRB+1
               CIJ(IPRB)= 0.0D0
   15       CONTINUE
         ENDIF
         IPRB= IPRB+1
         IDL(IUNK)= IPRB
         IF( IR.LT.0 )THEN
            IVOL= IVOL+1
            CIJ(IPRB)= PROB(IPRB)
         ELSEIF( IR.GT.0 )THEN
            IVOL= IVOL+IUNK-1
            IF( SIGTAL(IR).GT.0.0 )THEN
               CIJ(IPRB)= PROB(IPRB)
            ELSE
               CIJ(IPRB)= PROB(IVOL)
            ENDIF
         ELSE
            IVOL= IVOL+1
            CIJ(IPRB)= 1.0D0
         ENDIF
   20 CONTINUE
C
C     COPY MATRIX "CIJ" IN THE "WSPACE" ARRAY FOR INVERSION
C          AND ADD TO THE DIAGONAL ALL TERMS OF A LINE
      IPRB= 0
      IUNK= 0
      IDIA= 0
      DO 50 IR   = -NSOUT, NREG
         IUNK= IUNK+1
         IDIA= IDIA+IUNK
         WSPACE(IDIA)= CIJ(IDIA) + CIJ(IDIA)
         DO 30 JR= -NSOUT, IR-1
            IPRB= IPRB+1
            WSPACE(IPRB)= CIJ(IPRB)
            WSPACE(IDIA)= WSPACE(IDIA) + CIJ(IPRB)
   30    CONTINUE
         IPRB= IPRB+1
         IPRF= IPRB
         JUNK= IUNK
         DO 40 JR=  IR+1 , NREG
            IPRF= IPRF+JUNK
            JUNK= JUNK+1
            WSPACE(IDIA)= WSPACE(IDIA) + CIJ(IPRF)
   40    CONTINUE
   50 CONTINUE
      IF( IPRT.GT.100 )THEN
         WRITE(IUNOUT,8002)
         IPRB= 0
         DO 55 IR= -NSOUT, NREG
         DO 52 JR= -NSOUT, IR
            IPRB= IPRB+1
            WRITE(IUNOUT,8003) IR, JR, CIJ(IPRB),
     >                            WSPACE(IPRB),PROB(IPRB)
   52    CONTINUE
   55    CONTINUE
      ENDIF
C
C     INVERSION OF THE INITIAL SYSTEM JACOBIAN MATRIX
      CALL ALDDLF(NUNKNO,WSPACE,IDL)
C
C     INITIALISATION OF WEIGHTS
      DO 60 IR=-NSOUT, NREG
         WEIG(IR)=1.0D0
   60 CONTINUE
      WEIG(0)= 0.0D0
C
C     THE NON-LINEAR SYSTEM FOR WEIGHTS IS:
C      F1(W1, W2, ... WN)= W1*(W1*C11+W2*C12+ ... +WN*C1N) - TRUE1
C      F2(W1, W2, ... WN)= W2*(W1*C21+W2*C22+ ... +WN*C2N) - TRUE2
C      ...
C      FN(W1, W2, ... WN)= WN*(W1*CN1+W2*CN2+ ... +WN*CNN) - TRUEN
C     FORMING THE SYSTEM USING WEIGHTS "WEIG" & CONTRIBUTIONS "CIJ"
C
C     MAIN ITERATION LOOP
      DO 110 NIT=1,NITMAX
C
         IPRB= 0
         IUNK= 0
         IVOL= NSOUT*(NSOUT+1)/2
         DO 90 IR=-NSOUT, NREG
            IF( IR.LE.0 )THEN
               IVOL= IVOL+1
            ELSE
               IVOL= IVOL+IUNK
            ENDIF
            IUNK= IUNK+1
            WFSPAD= 0.0D0
            DO 70 JR=-NSOUT, IR
               IPRB= IPRB+1
               WFSPAD=WFSPAD+WEIG(JR)*CIJ(IPRB)
   70       CONTINUE
            IPRF= IPRB
            JUNK= IUNK
            DO 80 JR= IR+1 , NREG
               IPRF= IPRF+JUNK
               JUNK= JUNK+1
               WFSPAD=WFSPAD+WEIG(JR)*CIJ(IPRF)
   80       CONTINUE
            WFSP(IR)=WEIG(IR)*WFSPAD-PROB(IVOL)
   90    CONTINUE
         IF( IPRT.GT.100 )THEN
            WRITE(IUNOUT,9000)
            DO 92 IR= -NSOUT, NREG
               WRITE(IUNOUT,9001) IR, WFSP(IR)
   92       CONTINUE
         ENDIF
         CALL ALDDLS(NUNKNO,IDL,WSPACE,WFSP)
C
C        CALCULATIONS OF SQUARE DISTANCE BETWEEN 2 ITERATIONS
C        AND UPDATING THE SOLUTION
         TOTCON = 0.0D0
         DO 100 IR=-NSOUT, NREG
            TOTCON= TOTCON + WFSP(IR)**2
            WEIG(IR)= WEIG(IR) - WFSP(IR)
  100    CONTINUE
         IF( IPRT.GT.100 )THEN
            WRITE(IUNOUT,9004)
            DO 102 IR= -NSOUT, NREG
               WRITE(IUNOUT,9005) IR, WEIG(IR)
  102       CONTINUE
            WRITE(IUNOUT,'( 8H TOTCON: ,E15.7)') TOTCON
         ENDIF
C
C        CONVERGENCE TEST
         IF( TOTCON.LT.EPSCON )GO TO 120
C
  110 CONTINUE
      WRITE(IUNOUT,'(35H PIJRNL: WEIGHTS NOT CONVERGED          )')
  120 CONTINUE
C
C     RECOMPUTE WEIGHTS FOR VOID REGIONS
      IPRB= (NSOUT+1)*(NSOUT+2)/2
      IVOL= IPRB
      IUNK= NSOUT+1
      DO 220 IR= 1, NREG
         IVOL= IVOL+IUNK
         IUNK= IUNK+1
         IF( SIGTAL(IR).EQ.0.0 )THEN
            WFSPAD= 0.0D0
            DO 200 JR=-NSOUT, IR
               IPRB= IPRB+1
               IF( JR.LT.0.OR.SIGTAL(JR).GT.0.0 )THEN
                  WFSPAD=WFSPAD+WEIG(JR)*PROB(IPRB)
               ENDIF
  200       CONTINUE
            IPRF= IPRB
            JUNK= IUNK
            DO 210 JR= IR+1 , NREG
               IPRF= IPRF+JUNK
               JUNK= JUNK+1
               IF( JR.LT.0.OR.SIGTAL(JR).GT.0.0 )THEN
                  WFSPAD=WFSPAD+WEIG(JR)*PROB(IPRF)
               ENDIF
  210       CONTINUE
            WEIG(IR)=PROB(IVOL)/WFSPAD
         ELSE
            IPRB= IPRB+IUNK
         ENDIF
  220 CONTINUE
C
C     RENORMALIZE "PIJ" SYMMETRIC MATRIX
      IPRB = 0
      DO 240 IR   = -NSOUT, NREG
         DO 230 JR= -NSOUT, IR
            IPRB= IPRB+1
            IF( IR.NE.0.AND.JR.NE.0 )THEN
               PROB(IPRB)=PROB(IPRB)*WEIG(IR)*WEIG(JR)
            ENDIF
  230    CONTINUE
  240 CONTINUE
C
C     PRINT WEIGHT FACTORS IF REQUESTED
      IF( IPRT .GE. 100 )THEN
         WRITE(IUNOUT,'(30H0 SURFACE WEIGHTS FACTORS                /)')
         NSURC = -1
         DO 300 IP  = 1, (9 +NSOUT) / 10
            NSURM= MAX( -NSOUT, NSURC-9 )
            WRITE(IUNOUT,'(10X,10( A5,    I6)/)')
     >                     (' SUR ',-IR,IR= NSURC, NSURM, -1)
            WRITE(IUNOUT,'(10H WEIGHT   ,10F11.5)')
     >                   (WEIG(IR),IR=NSURC,NSURM,-1)
            NSURC = NSURC - 10
 300     CONTINUE
         WRITE(IUNOUT,'(30H0  VOLUME WEIGHTS FACTORS                /)')
         NVOLC =  1
         DO 310 IP  = 1, (9 + NREG) / 10
            NVOLM= MIN( NREG, NVOLC+9 )
            WRITE(IUNOUT,'(10X,10( A5 ,  I6)/)')
     >                  (' VOL ',IR,IR=NVOLC,NVOLM, 1)
            WRITE(IUNOUT,'(10H WEIGHT   ,10F11.5)')
     >                   (WEIG(IR),IR=NVOLC,NVOLM, 1)
            NVOLC = NVOLC + 10
 310     CONTINUE
      ENDIF
C----
C  SCRATCH STORAGE DEALLOCATION
C----
      DEALLOCATE(WEIG,WFSP,WSPACE,CIJ)
      DEALLOCATE(IDL)
      RETURN
C
 8002 FORMAT(//' S U R F / S U R F    C O N T R I B U T I O N S'//
     >9X,'BEGIN S',6X,'END S ',11X,'CIJ.    ',11X, 'WSPACE  ',
     >                         11X,'PROBS.  ')
 8003 FORMAT(6X,I10,5X,I10,5X,1P,E15.7,5X,E15.7,5X,E15.7 )
 9000 FORMAT(//' F U N C T I O N     V A L U E S'//
     >9X,'VOL/SUR',6X,'VALUE')
 9001 FORMAT(6X,I10,5X,F10.4)
 9004 FORMAT(//' W E I G H T E D     V A L U E S'//
     >9X,'VOL/SUR',6X,'VALUE')
 9005 FORMAT(6X,I10,5X,F10.4)
      END
