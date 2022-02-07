*DECK NSSDRV
      SUBROUTINE NSSDRV(IPTRK,IPMAC,IPFLX,NG,LX1,NMIX,ITRIAL,EPSOUT,
     1 MAXOUT,LNODF,BB2,IPRINT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Driver for the flux calculation with the nodal expansion method.
*
*Copyright:
* Copyright (C) 2021 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPTRK   nodal tracking.
* IPMAC   nodal macrolib.
* IPFLX   nodal flux.
* NG      number of energy groups.
* LX1     number of nodes in the nodal calculation.
* NMIX    number of mixtures in the nodal calculation.
* ITRIAL  type of expansion functions in the nodal calculation
*         (=1: polynomial; =2: hyperbolic).
* EPSOUT  convergence epsilon for the power method.
* MAXOUT  maximum number of iterations for the power method.
* LNODF   flag set to .true. to force discontinuity factors to one.
* BB2     imposed leakege used in non-regression tests.
* IPRINT  edition flag.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPMAC,IPFLX
      INTEGER NG,LX1,NMIX,ITRIAL(NMIX,NG),MAXOUT,IPRINT
      REAL EPSOUT,BB2
      LOGICAL LNODF
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      INTEGER ISTATE(NSTATE),NCODE(6),ICODE(6)
      REAL ZCODE(6)
      TYPE(C_PTR) JPMAC,KPMAC,JPFLX
      CHARACTER(LEN=8) HADF(2)
      CHARACTER(LEN=72) TITLE
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MAT,IJJ,NJJ,IPOS
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: KN,IQFR
      REAL, ALLOCATABLE, DIMENSION(:) :: XX,VOL,WORK,EVECT
      REAL, ALLOCATABLE, DIMENSION(:,:) :: BETA,DIFF,SIGR,CHI,SIGF,A,B,
     1 AI,A11,QFR,QFR2,COUR
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SCAT,FDXM,FDXP
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: FD
*
      ALB(X)=0.5*(1.0-X)/(1.0+X)
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(DIFF(NMIX,NG),SIGR(NMIX,NG),CHI(NMIX,NG),SIGF(NMIX,NG),
     1 SCAT(NMIX,NG,NG),BETA(NG,NG),FD(NMIX,2,NG,NG),COUR(LX1+1,NG))
      ALLOCATE(MAT(LX1),VOL(LX1),XX(LX1),KN(6,LX1),QFR(6,LX1),
     1 IQFR(6,LX1))
*----
*  RECOVER TRACKING INFORMATION
*----
      TITLE=' '
      CALL LCMLEN(IPTRK,'TITLE',LENGT,ITYLCM)
      IF(LENGT.GT.0) THEN
        CALL LCMGTC(IPTRK,'TITLE',72,1,TITLE)
        IF(IPRINT.GT.0) WRITE(6,'(/9H NSSDRV: ,A72)') TITLE
      ENDIF
      CALL LCMGET(IPTRK,'NCODE',NCODE)
      CALL LCMGET(IPTRK,'ICODE',ICODE)
      CALL LCMGET(IPTRK,'ZCODE',ZCODE)
      CALL LCMGET(IPTRK,'MATCOD',MAT)
      CALL LCMGET(IPTRK,'VOLUME',VOL)
      CALL LCMGET(IPTRK,'XX',XX)
      CALL LCMGET(IPTRK,'KN',KN)
      CALL LCMGET(IPTRK,'QFR',QFR)
      CALL LCMGET(IPTRK,'IQFR',IQFR)
*----
*  RECOVER MACROLIB INFORMATION
*----
      IF(BB2.NE.0.0) THEN
        IF(IPRINT.GT.0) WRITE(6,'(/32H NSSDRV: INCLUDE LEAKAGE IN THE ,
     >  13HMACROLIB (B2=,1P,E12.5,2H).)') BB2
      ENDIF
      CALL LCMGET(IPMAC,'STATE-VECTOR',ISTATE)
      JPMAC=LCMGID(IPMAC,'GROUP')
      ALLOCATE(WORK(NMIX*NG),IJJ(NMIX),NJJ(NMIX),IPOS(NMIX),
     1 FDXM(NMIX,NG,NG),FDXP(NMIX,NG,NG))
      DO IGR=1,NG
        KPMAC=LCMGIL(JPMAC,IGR)
        CALL LCMGET(KPMAC,'NTOT0',SIGR(1,IGR))
        CALL LCMGET(KPMAC,'DIFF',DIFF(1,IGR))
        CALL LCMGET(KPMAC,'CHI',CHI(1,IGR))
        CALL LCMGET(KPMAC,'NUSIGF',SIGF(1,IGR))
        CALL LCMGET(KPMAC,'IJJS00',IJJ)
        CALL LCMGET(KPMAC,'NJJS00',NJJ)
        CALL LCMGET(KPMAC,'IPOS00',IPOS)
        CALL LCMGET(KPMAC,'SCAT00',WORK)
        DO IBM=1,NMIX
          SCAT(IBM,IGR,:)=0.0
          IPOSDE=IPOS(IBM)-1
          DO JGR=IJJ(IBM),IJJ(IBM)-NJJ(IBM)+1,-1
            IPOSDE=IPOSDE+1
            IF(IPOSDE.GT.NMIX*NG) CALL XABORT('NSSDRV: SCAT OVERFLOW.')
            SCAT(IBM,IGR,JGR)=WORK(IPOSDE) ! IGR <-- JGR
          ENDDO
          SIGR(IBM,IGR)=SIGR(IBM,IGR)-SCAT(IBM,IGR,IGR)
        ENDDO
        IF(BB2.NE.0.0) THEN
          DO IBM=1,NMIX
            SIGR(IBM,IGR)=SIGR(IBM,IGR)+BB2*DIFF(IBM,IGR)
          ENDDO
        ENDIF
      ENDDO
      IF(ICODE(2).NE.0) THEN
        CALL LCMLEN(IPMAC,'ALBEDO',ILONG,ITYLCM)
        IF(ILONG.EQ.NG) THEN
          CALL LCMGET(IPMAC,'ALBEDO',WORK)
          BETA(:,:)=1.0
          DO IGR=1,NG
            BETA(IGR,IGR)=WORK(IGR)
          ENDDO
        ELSE IF(ILONG.EQ.NG*NG) THEN
          CALL LCMGET(IPMAC,'ALBEDO',BETA)
        ELSE
          CALL XABORT('NSSDRV: INVALID ALBEDO LENGTH.')
        ENDIF
      ENDIF
      FD(:,:,:,:)=0.0
      IF(LNODF.OR.ISTATE(12).EQ.0) THEN
        DO IBM=1,NMIX
          DO IGR=1,NG
            FD(IBM,1,IGR,IGR)=1.0
            FD(IBM,2,IGR,IGR)=1.0
          ENDDO
        ENDDO
      ELSE IF(ISTATE(12).EQ.1) THEN
        CALL XABORT('NSSDRV: CURRENT INFORMATION NOT SUPPORTED.')
      ELSE IF(ISTATE(12).EQ.2) THEN
        CALL LCMSIX(IPMAC,'ADF',1)
          CALL LCMGTC(IPMAC,'HADF',8,1,HADF(1))
          CALL LCMGET(IPMAC,HADF(1),WORK)
        CALL LCMSIX(IPMAC,' ',2)
        DO IBM=1,NMIX
          DO IGR=1,NG
            FD(IBM,1,IGR,IGR)=WORK((IGR-1)*NMIX+IBM)
            FD(IBM,2,IGR,IGR)=WORK((IGR-1)*NMIX+IBM)
          ENDDO
        ENDDO
      ELSE IF(ISTATE(12).EQ.3) THEN
        CALL LCMSIX(IPMAC,'ADF',1)
          CALL LCMGTC(IPMAC,'HADF',8,2,HADF)
          CALL LCMGET(IPMAC,HADF(1),FDXM)
          CALL LCMGET(IPMAC,HADF(2),FDXP)
        CALL LCMSIX(IPMAC,' ',2)
        DO JGR=1,NG
          DO IGR=1,NG
            FD(:NMIX,1,IGR,JGR)=FDXM(:NMIX,IGR,JGR)
            FD(:NMIX,2,IGR,JGR)=FDXP(:NMIX,IGR,JGR)
          ENDDO
        ENDDO
      ENDIF
      DEALLOCATE(FDXP,FDXM,IPOS,NJJ,IJJ,WORK)
*----
* INITIALIZE IPFLX OBJECT
*----
      ISTATE(:)=0
      ISTATE(1)=NG
      ISTATE(2)=5*LX1
      ISTATE(6)=2
      CALL LCMPUT(IPFLX,'STATE-VECTOR',NSTATE,1,ISTATE)
*----
*  COMPUTE NODAL SOLUTION
*----
      ALLOCATE(A(5*LX1*NG,5*LX1*NG),B(5*LX1*NG,5*LX1*NG))
      A(:5*LX1*NG,:5*LX1*NG)=0.0
      B(:5*LX1*NG,:5*LX1*NG)=0.0
      ALLOCATE(WORK(NMIX),A11(5*LX1,5*LX1),QFR2(6,LX1))
      DO J=1,NG
        IOF1=(J-1)*5*LX1
        DO I=1,NG
          IF(ICODE(2).NE.0) THEN
            DO IQW=1,2
              DO IEL=1,LX1
                IALB=IQFR(IQW,IEL)
                IF(IALB.NE.0) QFR2(IQW,IEL)=QFR(IQW,IEL)*ALB(BETA(I,J))
              ENDDO
            ENDDO
          ELSE
            IF(I == J) THEN
              QFR2(:6,:LX1)=QFR(:6,:LX1)
            ELSE
              QFR2(:6,:LX1)=0.0
            ENDIF
          ENDIF
          DO IBM=1,NMIX
            WORK(IBM)=CHI(IBM,I)*SIGF(IBM,J)
          ENDDO
          IOF2=(I-1)*5*LX1
          IF(I == J) THEN
            CALL NSS1TR(ITRIAL(1,J),LX1,NMIX,MAT,XX,KN,QFR2,DIFF(:,I),
     1      SIGR(:,I),FD(:,:,I,J),A11)
            A(IOF1+1:IOF1+LX1*5,IOF2+1:IOF2+LX1*5)=A11(:,:)
          ELSE
            CALL NSS2TR(ITRIAL(1,J),LX1,NMIX,MAT,XX,KN,QFR2,DIFF(:,J),
     1      SIGR(:,J),SCAT(:,I,J),FD(:,:,I,J),A11)
            A(IOF2+1:IOF2+LX1*5,IOF1+1:IOF1+LX1*5)=-A11(:,:)
          ENDIF
          CALL NSS3TR(ITRIAL(1,J),LX1,NMIX,MAT,XX,DIFF(:,J),SIGR(:,J),
     1    WORK(:),A11)
          B(IOF2+1:IOF2+LX1*5,IOF1+1:IOF1+LX1*5)=A11(:,:)
        ENDDO
      ENDDO
      DEALLOCATE(QFR2,A11,WORK)
*----
*  SOLVE EIGENVALUE MATRIX SYSTEM
*----
      EPS=1.0E-7
      CALL ALINV(5*LX1*NG,A,5*LX1*NG,IER)
      IF(IER.NE.0) CALL XABORT('NSSDRV: SINGULAR MATRIX')
      ALLOCATE(AI(5*LX1*NG,5*LX1*NG),EVECT(5*LX1*NG))
      AI(:5*LX1*NG,:5*LX1*NG)=MATMUL(A(:5*LX1*NG,:5*LX1*NG),
     1 B(:5*LX1*NG,:5*LX1*NG))
      EVECT(:)=1.0
      CALL AL1EIG(5*LX1*NG,AI,EPSOUT,MAXOUT,ITER,EVECT,EVAL)
      IF(IPRINT.GT.0) WRITE(6,10) EVAL,ITER
      DEALLOCATE(AI)
*----
*  NORMALIZE THE FLUX
*----
      FLXMAX=0.0
      DO IUN=1,5*LX1*NG
        IF(ABS(EVECT(IUN)).GT.ABS(FLXMAX)) FLXMAX=EVECT(IUN)
      ENDDO
      EVECT(:)=EVECT(:)/FLXMAX
*----
*  STORE FLUX IN IPFLX
*----
      JPFLX=LCMLID(IPFLX,'FLUX',NG)
      NUM0=0
      DO IG=1,NG
        DO KEL=1,LX1
          IBM=MAT(KEL)
          IOF=NUM0+(KEL-1)*5
          IF(ITRIAL(IBM,IG).EQ.1) THEN
            COUR(KEL,IG)=-(DIFF(IBM,IG)/XX(KEL))*(EVECT(IOF+2)-
     1      EVECT(IOF+3)+EVECT(IOF+4)/2.0-EVECT(IOF+5)/5.0)
          ELSE
            ETA=XX(KEL)*SQRT(SIGR(IBM,IG)/DIFF(IBM,IG))
            COUR(KEL,IG)=-(DIFF(IBM,IG)/XX(KEL))*(EVECT(IOF+2)-
     1      EVECT(IOF+3)+EVECT(IOF+4)*ETA*COSH(ETA/2)-
     2      EVECT(IOF+5)*ETA*SINH(ETA/2))
          ENDIF
        ENDDO
        IBM=MAT(LX1)
        IOF=NUM0+(LX1-1)*5
        IF(ITRIAL(IBM,IG).EQ.1) THEN
          COUR(LX1+1,IG)=-(DIFF(IBM,IG)/XX(LX1))*(EVECT(IOF+2)+
     1    EVECT(IOF+3)+EVECT(IOF+4)/2.0+EVECT(IOF+5)/5.0)
        ELSE
          ETA=XX(LX1)*SQRT(SIGR(IBM,IG)/DIFF(IBM,IG))
          COUR(LX1+1,IG)=-(DIFF(IBM,IG)/XX(LX1))*(EVECT(IOF+2)+
     1    EVECT(IOF+3)+EVECT(IOF+4)*ETA*COSH(ETA/2)+
     2    EVECT(IOF+5)*ETA*SINH(ETA/2))
        ENDIF
        IF(IPRINT.GT.0) THEN
          WRITE(6,'(/33H NSSDRV: AVERAGED FLUXES IN GROUP,I5)') IG
          WRITE(6,'(1P,10e12.4)') (EVECT(NUM0+(I-1)*5+1),I=1,LX1)
          WRITE(6,'(/39H NSSDRV: SURFACIC NET CURRENTS IN GROUP,I5)') IG
          WRITE(6,'(1P,10e12.4)') (COUR(I,IG),I=1,LX1+1)
        ENDIF
*----
* STORE NODAL SOLUTION IN IPFLX
*----
        CALL LCMPDL(JPFLX,IG,5*LX1,2,EVECT(NUM0+1))
        NUM0=NUM0+5*LX1
      ENDDO
      CALL LCMPUT(IPFLX,'K-EFFECTIVE',1,2,EVAL)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(EVECT,B,A)
      DEALLOCATE(COUR,FD,BETA,SCAT,SIGF,CHI,SIGR,DIFF,MAT,XX,VOL,KN,
     1 QFR,IQFR)
      RETURN
   10 FORMAT(14H NSSDRV: KEFF=,F11.8,12H OBTAINED IN,I5,11H ITERATIONS)
      END
