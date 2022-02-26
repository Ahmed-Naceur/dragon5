*DECK BREKOE
      SUBROUTINE BREKOE(IPMAC1,NC,NG,NMIX1,ISPH,B2,ENER,DC1,TOT1,SCAT1,
     1 JXM,FHETXM,IPRINT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Implement the 1D Koebke reflector model.
*
*Copyright:
* Copyright (C) 2021 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert and V. Salino
*
*Parameters: input
* IPMAC1  nodal macrolib.
* NC      number of sn macrolibs (=2: Koebke method).
* NG      number of energy groups.
* NMIX1   number of mixtures in the nodal calculation.
* ISPH    SPH flag (=0: use discontinuity factors; =1: use SPH factors).
* B2      buckling.
* ENER    energy limits.
* TOT1    total cross sections.
* SCAT1   scattering P0 cross sections.
* JXM     left boundary currents.
* FHETXM  left boundary fluxes.
* IPRINT  edition flag.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAC1
      INTEGER NC,NG,NMIX1,ISPH,IPRINT
      REAL B2(NC),ENER(NG+1),DC1(NMIX1,NG,NC),TOT1(NMIX1,NG,NC),
     1 SCAT1(NMIX1,NG,NG,NC),JXM(NMIX1,NG,NC),FHETXM(NMIX1,NG,NC)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      INTEGER ISTATE(NSTATE)
      CHARACTER HADF*8
      DOUBLE PRECISION R11,R21,R22,SIGR1,SIGR2,SIG21,D1,D2,A,B,C,F2
      TYPE(C_PTR) JPMAC1,KPMAC1
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,NJJ,IPOS
      REAL, ALLOCATABLE, DIMENSION(:) :: WORK,FDX,DIF
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IJJ(NMIX1),NJJ(NMIX1),IPOS(NMIX1),WORK(NG*NMIX1))
*----
*  RECOVER FLUX, MACROSCOPIC CROSS SECTIONS AND DIFFUSION COEFFICIENTS
*----
      IF(NC.NE.2) CALL XABORT('BREKOE: NC=2 EXPECTED.')
      IF(NG.NE.2) CALL XABORT('BREKOE: NG=2 EXPECTED.')
      IF(NMIX1.NE.1) CALL XABORT('BREKOE: NMIX1=1 EXPECTED.')
*----
*  COMPUTE EQUIVALENT REFLECTOR
*----
      ALLOCATE(FDX(NG),DIF(NG))
      IBM=1
      R11=.5*(FHETXM(IBM,1,1)/JXM(IBM,1,1)+FHETXM(IBM,1,2)/JXM(IBM,1,2))
      R21=(FHETXM(IBM,2,1)*JXM(IBM,2,2)-FHETXM(IBM,2,2)*JXM(IBM,2,1))/
     1 (JXM(IBM,1,1)*JXM(IBM,2,2)-JXM(IBM,1,2)*JXM(IBM,2,1))
      R22=(FHETXM(IBM,2,2)*JXM(IBM,1,1)-FHETXM(IBM,2,1)*JXM(IBM,1,2))/
     1 (JXM(IBM,1,1)*JXM(IBM,2,2)-JXM(IBM,1,2)*JXM(IBM,2,1))
      IF(IPRINT.GT.0) WRITE(6,10) R11,R21,R22
      SIGR1=.5*(TOT1(IBM,1,1)+TOT1(IBM,1,2)-SCAT1(IBM,1,1,1)-
     1          SCAT1(IBM,1,1,2)+B2(1)*DC1(IBM,1,1)+B2(2)*DC1(IBM,1,2))
      SIGR2=.5*(TOT1(IBM,2,1)+TOT1(IBM,2,2)-SCAT1(IBM,2,2,1)-
     1          SCAT1(IBM,2,2,2)+B2(1)*DC1(IBM,2,1)+B2(2)*DC1(IBM,2,2))
      SIG21=.5*(SCAT1(IBM,2,1,1)+SCAT1(IBM,2,1,2))
      IF(IPRINT.GT.0) WRITE(6,20) SIGR1,SIGR2,SIG21
      D1=1.0/(R11*R11*SIGR1)
      A=(R21*SIGR1-R22*SIG21)*SQRT(SIGR1/SIGR2)/(R22*R22)
      B=SIG21*SQRT(D1*SIGR2)
      C=-R21*D1*SIGR2*SQRT(SIGR1*SIGR2)
      F2=(-B+SQRT(B*B-4.0*A*C))/(2.0*A)
      D2=F2*F2/(R22*R22*SIGR2)
      IF(IPRINT.GT.0) WRITE(6,30) D1,D2,F2
      FDX(1)=1.0
      FDX(2)=REAL(F2)
      DIF(1)=REAL(D1)
      DIF(2)=REAL(D2)
      IF(IPRINT.GT.0) THEN
        WRITE(6,'(/37H BREKOE: KOEBKE DISCONTINUITY FACTORS)')
        WRITE(6,'(/12H MIXTURE=  1)')
        WRITE(6,'(6H  FDX=,1P,10E13.5,/(6X,10E13.5))') FDX(:NG)
      ENDIF
*----
*  APPLY SPH FACTORS
*----
      IF(ISPH.EQ.1) THEN
        DO IGR=1,NG
          TOT1(IBM,IGR,:2)=TOT1(IBM,IGR,:2)/FDX(IGR)
          DIF(IGR)=DIF(IGR)/FDX(IGR)
          DO JGR=1,NG
            SCAT1(IBM,IGR,JGR,:2)=SCAT1(IBM,IGR,JGR,:2)/FDX(JGR)
          ENDDO
        ENDDO
      ENDIF
      IF(IPRINT.GT.0) THEN
        WRITE(6,'(/38H BREKOE: KOEBKE DIFFUSION COEFFICIENTS)')
        WRITE(6,'(/12H MIXTURE=  1)')
        WRITE(6,'(6H DIFF=,1P,10E13.5,/(6X,10E13.5))') DIF(:NG)
      ENDIF
*----
*  SAVE THE OUTPUT NODAL MACROLIB
*----
      IBM=1
      ISTATE(:)=0
      ISTATE(1)=NG
      ISTATE(2)=NMIX1
      ISTATE(3)=1
      ISTATE(9)=1  ! diffusion coefficient information
      IF(ISPH.EQ.0) ISTATE(12)=3 ! discontinuity factor information
      CALL LCMPUT(IPMAC1,'STATE-VECTOR',NSTATE,1,ISTATE)
      CALL LCMPUT(IPMAC1,'ENERGY',NG+1,2,ENER)
      WORK(1)=1.0
      CALL LCMPUT(IPMAC1,'VOLUME',NMIX1,2,WORK)
      IF(ISPH.EQ.0) THEN
        CALL LCMSIX(IPMAC1,'ADF',1)
          NTYPE=1
          HADF='FD_B'
          CALL LCMPUT(IPMAC1,'NTYPE',1,1,NTYPE)
          CALL LCMPTC(IPMAC1,'HADF',8,1,HADF)
          CALL LCMPUT(IPMAC1,HADF,NG,2,FDX)
        CALL LCMSIX(IPMAC1,' ',2)
      ENDIF
      JPMAC1=LCMLID(IPMAC1,'GROUP',NG)
      DO IGR=1,NG
        KPMAC1=LCMDIL(JPMAC1,IGR)
        WORK(1)=1.0
        CALL LCMPUT(KPMAC1,'FLUX-INTG',NMIX1,2,WORK)
        WORK(1)=0.5*(TOT1(IBM,IGR,1)+TOT1(IBM,IGR,2))
        CALL LCMPUT(KPMAC1,'NTOT0',NMIX1,2,WORK)
        WORK(1)=0.0
        CALL LCMPUT(KPMAC1,'SIGW00',NMIX1,2,WORK)
        CALL LCMPUT(KPMAC1,'DIFF',NMIX1,2,DIF(IGR))
        IF(ISPH.EQ.1) THEN
          WORK(1)=1.0/FDX(IGR)
          CALL LCMPUT(KPMAC1,'NSPH',NMIX1,2,WORK)
        ENDIF
        IPOSDE=0
        DO J=1,NMIX1
          IBM=1
          J2=IGR
          J1=IGR
          DO JGR=1,NG
            IF(SCAT1(IBM,IGR,JGR,1)+SCAT1(IBM,IGR,JGR,2).NE.0.0) THEN
              J2=MAX(J2,JGR)
              J1=MIN(J1,JGR)
            ENDIF
          ENDDO
          NJJ(J)=J2-J1+1
          IJJ(J)=J2
          IPOS(J)=IPOSDE+1
          DO JGR=J2,J1,-1
            IPOSDE=IPOSDE+1
            IF(IPOSDE.GT.NG*NMIX1) CALL XABORT('BREKOE: SCAT OVERFLOW.')
            WORK(IPOSDE)=0.5*(SCAT1(IBM,IGR,JGR,1)+SCAT1(IBM,IGR,JGR,2))
          ENDDO
        ENDDO
        CALL LCMPUT(KPMAC1,'SCAT00',IPOSDE,2,WORK)
        CALL LCMPUT(KPMAC1,'NJJS00',NMIX1,1,NJJ)
        CALL LCMPUT(KPMAC1,'IJJS00',NMIX1,1,IJJ)
        CALL LCMPUT(KPMAC1,'IPOS00',NMIX1,1,IPOS)
      ENDDO
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(DIF,FDX,WORK,IPOS,NJJ,IJJ)
      RETURN
*
   10 FORMAT(13H BREKOE: R11=,1P,E12.4,5H R21=,E12.4,5H R22=,E12.4)
   20 FORMAT(15H BREKOE: SIGR1=,1P,E12.4,7H SIGR2=,E12.4,7H SIG21=,
     1 E12.4)
   30 FORMAT(12H BREKOE: D1=,1P,E12.4,4H D2=,E12.4,4H F2=,E12.4)
      END
