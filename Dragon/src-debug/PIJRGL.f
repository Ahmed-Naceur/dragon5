*DECK PIJRGL
      SUBROUTINE PIJRGL(IPRT,NREG,NSOUT,SIGTAL,PROB)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Gelbard normalization of collision probs (CP).
*
*Copyright:
* Copyright (C) 1991 Ecole Polytechnique de Montreal
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
*-----------------------------------------------------------------------
*
      IMPLICIT   NONE
      INTEGER    IPRT,NREG,NSOUT,IUNOUT,IPRINT,
     >           IPRB,IPRF,IUNK,JUNK,IVOL,JVOL,IR,JR,
     >           NSURM,NSURC,NVOLM,NVOLC,IP
      PARAMETER (IUNOUT=6, IPRINT=4)
      REAL       SIGTAL(-NSOUT:NREG)
      DOUBLE PRECISION PROB(*),RBARRE,GBARRE
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: RI
C----
C  SCRATCH STORAGE ALLOCATION
C----
      ALLOCATE(RI(-NSOUT:NREG))
C
      RBARRE=0.0
      GBARRE=0.0
      IPRB= 0
      IUNK= 0
      IVOL= NSOUT*(NSOUT+1)/2
C
C     COMPUTE R-SUB(I) FACTORS AND: GBARRE, RBARRE
      DO 30 IR=-NSOUT, NREG
         IUNK= IUNK+1
         RI(IR)=0.0
         DO 10 JR=-NSOUT, IR
            IPRB= IPRB+1
            IF( JR.LT.0.OR.SIGTAL(JR).GT.0.0 )THEN
               RI(IR)=RI(IR)+PROB(IPRB)
            ENDIF
   10    CONTINUE
         IPRF= IPRB
         JUNK= IUNK
         DO 20 JR= IR+1 , NREG
            IPRF= IPRF+JUNK
            JUNK= JUNK+1
            IF( JR.LT.0.OR.SIGTAL(JR).GT.0.0 )THEN
               RI(IR)=RI(IR)+PROB(IPRF)
            ENDIF
   20    CONTINUE
         IF( IR.LT.0 )THEN
            IVOL= IVOL+1
            RI(IR)= PROB(IVOL)-RI(IR)
            GBARRE= GBARRE+PROB(IVOL)
            RBARRE= RBARRE+RI(IR)
         ELSEIF( IR.GT.0 )THEN
            IVOL= IVOL+IUNK-1
            RI(IR)= PROB(IVOL)-RI(IR)
            IF( SIGTAL(IR).GT.0.0 )THEN
               GBARRE= GBARRE+PROB(IVOL)
               RBARRE= RBARRE+RI(IR)
            ENDIF
         ELSE
            IVOL= IVOL+1
            RI(IR)=0.0
         ENDIF
   30 CONTINUE
      GBARRE=1.0/GBARRE
      RBARRE=RBARRE*GBARRE
C
C     RENORMALIZE PROB MATRIX
      IVOL= NSOUT*(NSOUT+1)/2
      IPRB= 0
      IUNK= 0
      DO 210 IR   = -NSOUT, NREG
         IF( IR.LE.0 )THEN
            IVOL= IVOL+1
         ELSE
            IVOL= IVOL+IUNK
         ENDIF
         IUNK= IUNK+1
         JVOL= NSOUT*(NSOUT+1)/2
         JUNK= 0
         DO 200 JR= -NSOUT, IR
            IF( JR.LE.0 )THEN
               JVOL= JVOL+1
            ELSE
               JVOL= JVOL+JUNK
            ENDIF
            JUNK= JUNK+1
            IPRB= IPRB+1
            IF( IR.NE.0.AND.JR.NE.0 )THEN
               PROB(IPRB)= PROB(IPRB)+(PROB(JVOL)*RI(IR)
     >         +PROB(IVOL)*RI(JR)-PROB(IVOL)*PROB(JVOL)*RBARRE)*GBARRE
            ENDIF
 200     CONTINUE
 210  CONTINUE
C
C     PRINT IF REQUESTED
      IF( IPRT.GE.IPRINT )THEN
         WRITE(IUNOUT,'(19H0  GLOBAL FACTORS: ,
     >         8H RBARRE=,1P,F11.5,5X,7HGBARRE=,F11.5)')
     >            RBARRE,               GBARRE
         WRITE(IUNOUT,'(30H0 SURFACE ADJUSTMENT FACTORS             /)')
         NSURC = -1
         DO 300 IP  = 1, (9 +NSOUT) / 10
            NSURM= MAX( -NSOUT, NSURC-9 )
            WRITE(IUNOUT,'(10X,10( A5,    I6)/)')
     >                     (' SUR ',-IR,IR= NSURC, NSURM, -1)
            WRITE(IUNOUT,'(10H R-SUB(I) ,10F11.5)')
     >                   (RI(IR),IR=NSURC,NSURM,-1)
            NSURC = NSURC - 10
 300     CONTINUE
         WRITE(IUNOUT,'(30H0  VOLUME ADJUSTMENT FACTORS             /)')
         NVOLC =  1
         DO 310 IP  = 1, (9 + NREG) / 10
            NVOLM= MIN( NREG, NVOLC+9 )
            WRITE(IUNOUT,'(10X,10( A5 ,  I6)/)')
     >                  (' VOL ',IR,IR=NVOLC,NVOLM, 1)
            WRITE(IUNOUT,'(10H R-SUB(I) ,10F11.5)')
     >                   (RI(IR),IR=NVOLC,NVOLM, 1)
            NVOLC = NVOLC + 10
 310     CONTINUE
      ENDIF
C----
C  SCRATCH STORAGE DEALLOCATION
C----
      DEALLOCATE(RI)
C
      RETURN
      END
