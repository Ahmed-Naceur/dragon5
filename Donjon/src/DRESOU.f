*DECK DRESOU
      SUBROUTINE DRESOU(IPRINT,IPGPT,IPMAC1,IPMAC2,IPFLX,IPGRAD,NG,NREG,
     1 NMIL,NUN,MATCOD,KEYFLX,VOL,LNO,RMSD)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the GPT sources corresponding to the gradient of the RMS power
* distribution. Case with no direct effect.
*
*Copyright:
* Copyright (C) 2012 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* IPRINT  print parameter
* IPGPT   pointer to the L_GPT data structure.
* IPMAC1  pointer to the actual macrolib structure.
* IPMAC2  pointer to the reference macrolib structure.
* IPFLX   pointer to the multigroup flux.
* IPGRAD  pointer to the L_OPTIMIZE object.
* NG      number of energy groups.
* NREG    number of regions.
* NMIL    number of material mixtures.
* NUN     number of unknowns per energy group.
* MATCOD  material mixture indices per region.
* KEYFLX  position of averaged fluxes in unknown vector.
* VOL     volumes.
* LNO     flag set to .true. to exit after calculation of RMS.
*
*Parameters: output
* RMSD    RMS error on power distribution.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPGPT,IPMAC1,IPMAC2,IPFLX,IPGRAD
      INTEGER IPRINT,NG,NREG,NMIL,NUN,MATCOD(NREG),KEYFLX(NREG)
      REAL VOL(NREG)
      DOUBLE PRECISION RMSD
      LOGICAL LNO
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPMAC1,JPMAC2,KPMAC1,KPMAC2,JPFLX,JPGPT,KPGPT
      DOUBLE PRECISION SOUT2,SOUTOT,PW1TOT,PW2TOT,DSUM,AIL,BIL
      REAL, ALLOCATABLE, DIMENSION(:) :: POW1,H1,F1,POW2,H2,F2,SUNK,
     1 FLUX
      CHARACTER HSMG*131
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(POW1(NMIL),H1(NMIL),F1(NMIL),POW2(NMIL),H2(NMIL),
     1 F2(NMIL))
*----
*  COMPUTE THE ACTUAL AND REFERENCE POWER DISTRIBUTION
*----
      CALL XDRSET(POW1,NMIL,0.0)
      CALL XDRSET(POW2,NMIL,0.0)
      JPMAC1=LCMGID(IPMAC1,'GROUP')
      JPMAC2=LCMGID(IPMAC2,'GROUP')
      DO IG=1,NG
         KPMAC1=LCMGIL(JPMAC1,IG)
         KPMAC2=LCMGIL(JPMAC2,IG)
         CALL LCMLEN(KPMAC1,'FLUX-INTG',ILG,ITYLCM)
         IF(ILG.EQ.0) CALL XABORT('DRESOU: MISSING ACTUAL FLUX.')
         CALL LCMLEN(KPMAC2,'FLUX-INTG',ILG,ITYLCM)
         IF(ILG.EQ.0) CALL XABORT('DRESOU: MISSING REFERENCE FLUX.')
         CALL LCMLEN(KPMAC1,'H-FACTOR',ILG,ITYLCM)
         IF(ILG.EQ.0) CALL XABORT('DRESOU: MISSING ACTUAL H-FACTOR.')
         CALL LCMLEN(KPMAC2,'H-FACTOR',ILG,ITYLCM)
         IF(ILG.EQ.0) CALL XABORT('DRESOU: MISSING REFERENCE H-FACTOR.')
         CALL LCMGET(KPMAC1,'FLUX-INTG',F1)
         CALL LCMGET(KPMAC2,'FLUX-INTG',F2)
         CALL LCMGET(KPMAC1,'H-FACTOR',H1)
         CALL LCMGET(KPMAC2,'H-FACTOR',H2)
         DO IBM=1,NMIL
            POW1(IBM)=POW1(IBM)+F1(IBM)*H1(IBM)
            POW2(IBM)=POW2(IBM)+F2(IBM)*H2(IBM)
         ENDDO
      ENDDO
*----
*  COMPUTE THE RMS FUNCTIONAL
*----
      PW1TOT=0.0D0
      PW2TOT=0.0D0
      DO IBM=1,NMIL
         PW1TOT=PW1TOT+POW1(IBM)
         PW2TOT=PW2TOT+POW2(IBM)
      ENDDO
      RMSD=0.0D0
      DO IBM=1,NMIL
         RMSD=RMSD+(POW1(IBM)/PW1TOT-POW2(IBM)/PW2TOT)**2
      ENDDO
      CALL LCMPUT(IPGRAD,'FOBJ-CST-VAL',1,4,RMSD)
      IF(IPRINT.GT.0) WRITE(6,100) RMS
      IF(LNO) GO TO 10
*----
*  COMPUTE THE GRADIENT OF THE RMS FUNCTIONAL
*----
      ALLOCATE(SUNK(NUN))
      JPFLX=LCMGID(IPFLX,'FLUX')
      JPGPT=LCMLID(IPGPT,'ASOUR',1)
      KPGPT=LCMLIL(JPGPT,1,NG)
      DO IG=1,NG
         CALL XDRSET(SUNK,NUN,0.0)
         KPMAC1=LCMGIL(JPMAC1,IG)
         CALL LCMGET(KPMAC1,'H-FACTOR',H1)
         DO IR=1,NREG
            IUNK=KEYFLX(IR)
            IF(IUNK.EQ.0) CYCLE
            IBM=MATCOD(IR)
            IF(IBM.EQ.0) CYCLE
            SOUT2=0.0D0
            DO JBM=1,NMIL
              SOUTOT=0.0D0
              IF(IBM.EQ.JBM) SOUTOT=1.0D0
              SOUTOT=SOUTOT-POW1(JBM)/PW1TOT
              SOUT2=SOUT2+SOUTOT*(POW1(JBM)/PW1TOT-POW2(JBM)/PW2TOT)
            ENDDO
            SUNK(IUNK)=2.0*VOL(IR)*H1(IBM)*REAL(SOUT2/PW1TOT)
         ENDDO
         CALL LCMPDL(KPGPT,IG,NUN,2,SUNK)
      ENDDO
*----
*  CHECK SOURCE ORTHOGONALITY
*----
      ALLOCATE(FLUX(NUN))
      AIL=0.0D0
      BIL=0.0D0
      DO IG=1,NG
        CALL LCMGDL(KPGPT,IG,SUNK)
        CALL LCMGDL(JPFLX,IG,FLUX)
        DO IUNK=1,NUN
          GAZ=FLUX(IUNK)*SUNK(IUNK)
          DAZ=FLUX(IUNK)**2
          AIL=AIL+GAZ
          BIL=BIL+DAZ
        ENDDO
      ENDDO
      DSUM=ABS(AIL)/ABS(BIL)/REAL(NUN)
      IF(IPRINT.GT.0) THEN
        WRITE(6,'(/21H DRESOU: DOT PRODUCT=,1P,E11.4)') DSUM
      ENDIF
      IF(ABS(DSUM).GT.1.0E-5) THEN
        WRITE(HSMG,'(36HDRESOU: NON ORTHOGONAL SOURCE (DSUM=,1P,E11.3,
     1  2H).)') DSUM
        CALL XABORT(HSMG)
      ENDIF
      DEALLOCATE(FLUX,SUNK)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
   10 DEALLOCATE(F2,H2,POW2,F1,H1,POW1)
      RETURN
*
  100 FORMAT(/41H DRESOU: RMS ERROR ON POWER DISTRIBUTION=,1P,E11.4)
      END
