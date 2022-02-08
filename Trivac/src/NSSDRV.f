*DECK NSSDRV
      SUBROUTINE NSSDRV(IPTRK,IPMAC,IPFLX,LNEM,NUN,NG,LX1,NMIX,ITRIAL,
     1 EPSOUT,MAXOUT,LNODF,BB2,IPRINT)
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
* LNEM    solution flag (=.true.: NEM; =.false.:CMFD).
* NUN     number of unknowns.
* NG      number of energy groups.
* LX1     number of nodes in the nodal calculation.
* NMIX    number of mixtures in the nodal calculation.
* ITRIAL  type of expansion functions in the nodal calculation
*         (=1: polynomial; =2: hyperbolic).
* EPSOUT  convergence epsilon for the power method.
* MAXOUT  maximum number of iterations for the power method.
* LNODF   flag set to .true. to force discontinuity factors to one.
* BB2     imposed leakage used in non-regression tests.
* IPRINT  edition flag.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPMAC,IPFLX
      INTEGER NUN,NG,LX1,NMIX,ITRIAL(NMIX,NG),MAXOUT,IPRINT
      REAL EPSOUT,BB2
      LOGICAL LNEM,LNODF
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      INTEGER ISTATE(NSTATE),ICODE(6)
      TYPE(C_PTR) JPMAC,KPMAC,JPFLX
      CHARACTER(LEN=8) HADF(2)
      CHARACTER(LEN=72) TITLE
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MAT,IJJ,NJJ,IPOS
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IQFR
      REAL, ALLOCATABLE, DIMENSION(:) :: XX,WORK
      REAL, ALLOCATABLE, DIMENSION(:,:) :: DIFF,SIGR,CHI,SIGF,QFR,ALBP,
     1 EVECT 
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: BETA,SCAT,FDXM,FDXP
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: FD
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(DIFF(NMIX,NG),SIGR(NMIX,NG),CHI(NMIX,NG),SIGF(NMIX,NG),
     1 SCAT(NMIX,NG,NG),FD(NMIX,2,NG,NG))
      ALLOCATE(MAT(LX1),XX(LX1),QFR(6,LX1),IQFR(6,LX1))
*----
*  RECOVER TRACKING INFORMATION
*----
      TITLE=' '
      CALL LCMLEN(IPTRK,'TITLE',LENGT,ITYLCM)
      IF(LENGT.GT.0) THEN
        CALL LCMGTC(IPTRK,'TITLE',72,1,TITLE)
        IF(IPRINT.GT.0) WRITE(6,'(/9H NSSDRV: ,A72)') TITLE
      ENDIF
      CALL LCMGET(IPTRK,'ICODE',ICODE)
      CALL LCMGET(IPTRK,'MATCOD',MAT)
      CALL LCMGET(IPTRK,'XX',XX)
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
      NALB=ISTATE(8) ! number of physical albedos
      JPMAC=LCMGID(IPMAC,'GROUP')
      ALLOCATE(WORK(NMIX*NG),IJJ(NMIX),NJJ(NMIX),IPOS(NMIX),
     1 FDXM(NMIX,NG,NG),FDXP(NMIX,NG,NG),BETA(NALB,NG,NG))
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
        DO IBM=1,NMIX
          IF(SIGR(IBM,IGR).LE.0.0) CALL XABORT('NSSDRV: SIGR<=0.')
        ENDDO
      ENDDO
      IF(NALB.GT.0) THEN
        CALL LCMLEN(IPMAC,'ALBEDO',ILONG,ITYLCM)
        IF(ILONG.EQ.NALB*NG) THEN
          ALLOCATE(ALBP(NALB,NG))
          CALL LCMGET(IPMAC,'ALBEDO',ALBP)
          BETA(:,:,:)=1.0
          DO IGR=1,NG
            BETA(:NALB,IGR,IGR)=ALBP(:NALB,IGR)
          ENDDO
          DEALLOCATE(ALBP)
        ELSE IF(ILONG.EQ.NALB*NG*NG) THEN
          CALL LCMGET(IPMAC,'ALBEDO',BETA)
        ELSE
          CALL XABORT('NSSDRV: INVALID ALBEDO LENGTH.')
        ENDIF
        IF(IPRINT.GT.1) THEN
          DO IALB=1,NALB
            WRITE(6,'(/35H NSSDRV: PHYSICAL ALBEDO MATRIX ID=,I4)') IALB
            DO IGR=1,NG
              WRITE(6,'(5X,1P,10E12.4)') BETA(IALB,IGR,:)
            ENDDO
          ENDDO
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
      ELSE IF(ISTATE(12).EQ.3) THEN
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
      ELSE IF(ISTATE(12).EQ.4) THEN
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
      ELSE
        CALL XABORT('NSSDRV: FLUX/CURRENT INFORMATION NOT SUPPORTED.')
      ENDIF
      DEALLOCATE(FDXP,FDXM,IPOS,NJJ,IJJ,WORK)
*----
*  COMPUTE THE FLUX AND STORE NODAL SOLUTION IN IPFLX
*----
      ALLOCATE(EVECT(NUN,NG))
      IF(LNEM) THEN
        CALL NSSFL1(NUN,NG,LX1,NMIX,NALB,ITRIAL,EPSOUT,MAXOUT,MAT,XX,
     1  IQFR,QFR,DIFF,SIGR,CHI,SIGF,SCAT,BETA,FD,IPRINT,EVAL,EVECT)
      ELSE
        CALL NSSFL2(NUN,NG,LX1,NMIX,NALB,EPSOUT,MAXOUT,MAT,XX,IQFR,
     1  QFR,DIFF,SIGR,CHI,SIGF,SCAT,BETA,FD,IPRINT,EVAL,EVECT)
      ENDIF
      ISTATE(:)=0
      ISTATE(1)=NG
      ISTATE(2)=NUN
      ISTATE(6)=2
      CALL LCMPUT(IPFLX,'STATE-VECTOR',NSTATE,1,ISTATE)
      JPFLX=LCMLID(IPFLX,'FLUX',NG)
      DO IGR=1,NG
        CALL LCMPDL(JPFLX,IGR,NUN,2,EVECT(1,IGR))
      ENDDO
      CALL LCMPUT(IPFLX,'K-EFFECTIVE',1,2,EVAL)
      DEALLOCATE(EVECT)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(FD,BETA,SCAT,SIGF,CHI,SIGR,DIFF,MAT,XX,QFR,IQFR)
      RETURN
      END
