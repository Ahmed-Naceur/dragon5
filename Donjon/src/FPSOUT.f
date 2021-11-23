*DECK FPSOUT
      SUBROUTINE FPSOUT(IPMAC,IPRINT,NG,NMIL,NFIS,ILEAKS,TEXT9,OUTG)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the leakage rate in each energy group
*
*Copyright:
* Copyright (C) 2019 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* IPMAC   pointer to the macrolib structure.
* IPRINT  print parameter
* NG      number of energy groups.
* NMIL    number of material mixtures.
* NFIS    number of fissile isotopes.
* ILEAKS  type of leakage calculation =0: no leakage; =1: homogeneous
*         leakage (Diffon).
* TEXT9   type of calculation ('REFERENCE' or 'MACRO').
*
*Parameters: output
* OUTG    leakage rates.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAC
      INTEGER IPRINT,NG,NMIL,NFIS,ILEAKS
      CHARACTER TEXT9*9
      REAL OUTG(NG)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPMAC,KPMAC
      CHARACTER HSMG*131
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,NJJ,IPOS
      REAL, ALLOCATABLE, DIMENSION(:) :: GAR,WORK,DIFHOM,DIFF
      REAL, ALLOCATABLE, DIMENSION(:,:) :: PHI,NUF
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: CHI,RHS,LHS
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(PHI(NMIL,NG),RHS(NMIL,NG,NG),LHS(NMIL,NG,NG))
      ALLOCATE(IJJ(NMIL),NJJ(NMIL),IPOS(NMIL),GAR(NMIL),WORK(NMIL*NG),
     > CHI(NMIL,NFIS,NG),NUF(NMIL,NFIS),DIFHOM(NG),DIFF(NMIL))
*----
*  COMPUTE THE ACTUAL AND REFERENCE REACTION RATE MATRICES
*----
      CALL LCMGET(IPMAC,'K-EFFECTIVE',ZKEFF)
      IF(IPRINT.GT.1) WRITE(6,120) TEXT9,ZKEFF
      CALL LCMLEN(IPMAC,'B2  B1HOM',ILCMLN,ITYLCM)
      IF(ILCMLN.EQ.1) THEN
        CALL LCMGET(IPMAC,'B2  B1HOM',B2)
      ELSE
        B2=0.0
      ENDIF
      IF((ILEAKS.EQ.1).AND.(IPRINT.GT.1)) THEN
        WRITE(6,'(/9H FPSOUT: ,A,4H B2=,1P,E12.4)') TEXT9,B2
      ENDIF
      CALL XDRSET(RHS,NMIL*NG*NG,0.0)
      CALL XDRSET(LHS,NMIL*NG*NG,0.0)
      JPMAC=LCMGID(IPMAC,'GROUP')
      DO IG=1,NG
        KPMAC=LCMGIL(JPMAC,IG)
        CALL LCMGET(KPMAC,'CHI',CHI(1,1,IG))
        CALL LCMLEN(KPMAC,'FLUX-INTG',ILG,ITYLCM)
        IF(ILG.NE.NMIL) CALL XABORT('FPSOUT: MISSING REFERENCE FLUX.')
        CALL LCMGET(KPMAC,'FLUX-INTG',PHI(1,IG))
      ENDDO
      DO IG=1,NG
        KPMAC=LCMGIL(JPMAC,IG)
        IF(ILEAKS.EQ.1) THEN
          CALL LCMLEN(KPMAC,'DIFF',ILCMLN,ITYLCM)
          IF(ILCMLN.GT.0) THEN
            CALL LCMGET(KPMAC,'DIFF',DIFF)
          ELSE
            CALL LCMGET(IPMAC,'DIFHOMB1HOM',DIFHOM)
            DO IBM=1,NMIL
              DIFF(IBM)=DIFHOM(IG)
            ENDDO
          ENDIF
        ELSE
          CALL XDRSET(DIFF,NMIL,0.0)
        ENDIF
        CALL LCMGET(KPMAC,'NTOT0',GAR)
        CALL LCMGET(KPMAC,'SCAT00',WORK)
        CALL LCMGET(KPMAC,'NJJS00',NJJ)
        CALL LCMGET(KPMAC,'IJJS00',IJJ)
        CALL LCMGET(KPMAC,'IPOS00',IPOS)
        DO IBM=1,NMIL
          IPOSDE=IPOS(IBM)
          DO JG=IJJ(IBM),IJJ(IBM)-NJJ(IBM)+1,-1
*           IG <-- JG
            RHS(IBM,IG,JG)=RHS(IBM,IG,JG)-WORK(IPOSDE)*PHI(IBM,JG)
            IPOSDE=IPOSDE+1
          ENDDO
          RHS(IBM,IG,IG)=RHS(IBM,IG,IG)+(GAR(IBM)+B2*DIFF(IBM))*
     >    PHI(IBM,IG)
        ENDDO
        CALL LCMGET(KPMAC,'NUSIGF',NUF)
        DO IBM=1,NMIL
          DO IFIS=1,NFIS
            DO JG=1,NG
              LHS(IBM,JG,IG)=LHS(IBM,JG,IG)+CHI(IBM,IFIS,JG)*
     >        NUF(IBM,IFIS)*PHI(IBM,IG)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
*----
*  COMPUTE THE ACTUAL AND REFERENCE ABSORPTION AND FISSION RATES
*----
      DO IG=1,NG
        OUTG(IG)=0.0
        DO IBM=1,NMIL
          OUTG(IG)=OUTG(IG)+SUM(LHS(IBM,IG,:NG))/ZKEFF-
     1    SUM(RHS(IBM,IG,:NG))
        ENDDO
        IF(OUTG(IG).LT.-1.0E-6) THEN
          WRITE(HSMG,'(21HFPSOUT: INCONSISTENT ,A,17H LEAKAGE IN GROUP,
     1    I4,7H. LEAK=,1P,E13.4)') TEXT9,IG,OUTG(IG)
          CALL XABORT(HSMG)
        ENDIF
        IF(IPRINT.GT.1) WRITE(6,130) IG,TEXT9,OUTG(IG)
      ENDDO
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(DIFF,DIFHOM,NUF,CHI,WORK,GAR,IPOS,NJJ,IJJ)
      DEALLOCATE(LHS,RHS,PHI)
      RETURN
*
  120 FORMAT(/9H FPSOUT: ,A,33H EFFECTIVE MULTIPLICATION FACTOR=,1P,
     1 E12.4)
  130 FORMAT(/8H FPSOUT:,5X,6HGROUP=,I4,1X,A,9H LEAKAGE=,1P,E12.4)
      END
