*DECK PIJS3D
      SUBROUTINE PIJS3D(NREG,NSOUT,NSLINE,NSBG,WEIGHT,
     >                  RCUTOF,SIGTAL,NPSYS,
     >                  SEGLEN,NRSEG,
     >                  STAYIN,GOSOUT,DPR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Integration for general 3D specular tracking.
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
* NREG    total number of regions.
* NSOUT   number of outer surface.
* NSLINE  number of segemnts on line.
* NSBG    number of subgroup.
* WEIGHT  line weight.
* RCUTOF  MFP cut-off factor (truncate lines).
* SIGTAL  albedo-cross section vector.
* NPSYS   non-converged energy group indices.
* SEGLEN  length of track.
* NRSEG   region crossed by track.
*
*Parameters: output
* DPR     collision probabilities.
*
*Parameters: scratch
* STAYIN  stay-in zone probability.
* GOSOUT  goes-out zone probability.
*
*References:
* R. Roy et al.,
* A Cyclic Tracking Procedure for CP Calculations in 2-D Lattices
* Conf/Advances in Math, Comp & Reactor Physics,
* Pittsburgh, V 1, P 2.2 4-1 (1991).
*
*-----------------------------------------------------------------------
*
      IMPLICIT         NONE
*----
* VARIABLES
*----
      INTEGER          NREG,NSOUT,NSLINE,NSBG
      INTEGER          NRSEG(NSLINE),NPSYS(NSBG)
      REAL             RCUTOF,SIGTAL(-NSOUT:NREG,NSBG)
      DOUBLE PRECISION WEIGHT,SEGLEN(NSLINE),STAYIN(NSLINE),
     >                 GOSOUT(NSLINE)
      DOUBLE PRECISION DPR(-NSOUT:NREG,-NSOUT:NREG,NSBG)
*----
*  Local variables
*----
      INTEGER          ISBG,IL,JL,NOIL,NOJL,ISODD,JSODD,IJDEL
      DOUBLE PRECISION TTOT,XSIL,OPATH,FINV,CUTOF
      REAL             ZERO,ONE,HALF
      PARAMETER       (ZERO=0.0E0, ONE=1.0E0, HALF=0.5E0 )
      REAL             SIXT,CUTEXP
      PARAMETER       (SIXT=HALF/3.0,CUTEXP=0.02)
      DOUBLE PRECISION EXSIL,XSIL2
      TTOT= ONE
      DO 2001 ISBG=1,NSBG
        IF(NPSYS(ISBG).EQ.0) GO TO 2001
*
*1.1)    CHANGE PATHS => GOSOUT AND STAYIN PATHS, INCLUDING ALBEDOS
*        ADD *PII* LOCAL NON-CYCLIC CONTRIBUTIONS
        ISODD=0
        DO 30 IL= 1, NSLINE
          NOIL  = NRSEG(IL)
          IF( NOIL.LT.0 )THEN
            IF(ISODD .EQ. 1) THEN
              ISODD=0
*----
*  FOR SURFACES:
*    OLD VERSION BEFORE SURFACE DOUBLING
*      GOSOUT= ALBEDO * SURFACE WEIGHT
*      WHERE ALL SURFACE WEIGHTS WERE 1.0
*    NEW VERSION WITH SURFACE DOUBLING
*      GOSOUT= ALBEDO
*    STAYIN = 1- ALBEDO * SURFACE WEIGHT
*    TTOT   = PRODUCT OF GOSOUT
*----
              GOSOUT(IL)= SIGTAL(NOIL,ISBG)
              STAYIN(IL)= ONE - GOSOUT(IL)
              TTOT= TTOT * GOSOUT(IL)
            ELSE
              ISODD=1
*----
*  FOR SURFACES:
*    OLD VERSION BEFORE SURFACE DOUBLING
*      GOSOUT= ALBEDO * SURFACE WEIGHT
*      WHERE ALL SURFACE WEIGHTS WERE 1.0
*    NEW VERSION WITH SURFACE DOUBLING
*      GOSOUT= ALBEDO
*    STAYIN = 1- ALBEDO * SURFACE WEIGHT
*    TTOT   = PRODUCT OF GOSOUT
*----
              GOSOUT(IL)= SIGTAL(NOIL,ISBG)
              STAYIN(IL)= ONE
            ENDIF
          ELSE
*----
*  FOR REGIONS
*  STAYIN = 1 -  EXP[ -CROSS SECTION * LENGTH OF NSLINE]
*  GOSOUT = 1 -  STAYIN
*  TTOT   = PRODUCT OF GOSOUT
*----
            XSIL  = SIGTAL(NOIL,ISBG)
            IF( XSIL .EQ. ZERO) THEN
              GOSOUT(IL)= ONE
              STAYIN(IL)= SEGLEN(IL)
              DPR(NOIL,NOIL,ISBG)= DPR(NOIL,NOIL,ISBG)
     >                           + HALF*WEIGHT*STAYIN(IL)*STAYIN(IL)
            ELSE IF( XSIL .LT. CUTEXP) THEN
              OPATH= SIGTAL(NOIL,ISBG)*SEGLEN(IL)
              XSIL2=OPATH*OPATH
              EXSIL=XSIL2*(HALF-SIXT*OPATH+XSIL2/24.0)
              STAYIN(IL)=OPATH-EXSIL
              GOSOUT(IL)= ONE - STAYIN(IL)
              TTOT= TTOT * GOSOUT(IL)
              DPR(NOIL,NOIL,ISBG)= DPR(NOIL,NOIL,ISBG) + WEIGHT*EXSIL
            ELSE
              OPATH= SIGTAL(NOIL,ISBG)*SEGLEN(IL)
              EXSIL= EXP(-OPATH)
              STAYIN(IL)= ONE - EXSIL
              GOSOUT(IL)= EXSIL
              TTOT= TTOT * GOSOUT(IL)
              DPR(NOIL,NOIL,ISBG)= DPR(NOIL,NOIL,ISBG)
     >                           + WEIGHT*(OPATH-STAYIN(IL))
            ENDIF
          ENDIF
 30     CONTINUE
*
*1.2)    COMPUTE CYCLIC FACTORS BY ANGLE
*        USING GLOBAL TRACK ATTENUATION: BETA(TOT)*EXP(-MFP(TOT))
         IF(TTOT .GE. ONE )THEN
           CALL XABORT( 'PIJS3D: ALBEDOS ARE NOT COMPATIBLE')
         ENDIF
         FINV= WEIGHT / (ONE-TTOT)
*
*1.3)    ADD *PIJ* CONTRIBUTIONS FOR FORWARD SOURCES
        ISODD=0
        DO 50 IL= 1, NSLINE
          NOIL  = NRSEG(IL)
          TTOT= FINV * STAYIN(IL)
          CUTOF= RCUTOF*TTOT
          IF( NOIL .LT. 0) THEN
            ISODD=MOD(ISODD+1,2)
            JSODD=ISODD
            DO 70 IJDEL= 1, NSLINE
              JL= MOD(IL+IJDEL-1,NSLINE) + 1
              NOJL=NRSEG(JL)
              IF( NOJL .LT. 0 ) THEN
                JSODD=MOD(JSODD+1,2)
                IF( ISODD.EQ.1 .AND. JSODD .EQ.0) THEN
                  DPR(NOJL,NOIL,ISBG)= DPR(NOJL,NOIL,ISBG)
     >                               + TTOT * STAYIN(JL)
                  TTOT= TTOT * GOSOUT(JL)
                  IF( TTOT.LE.CUTOF ) GO TO 55
                ENDIF
              ELSE IF(ISODD.EQ.1) THEN
                DPR(NOJL,NOIL,ISBG)= DPR(NOJL,NOIL,ISBG)
     >                             + TTOT * STAYIN(JL)
                TTOT= TTOT * GOSOUT(JL)
                IF( TTOT.LE.CUTOF ) GO TO 55
              ENDIF
 70         CONTINUE
          ELSE
            JSODD=ISODD
            DO 80 IJDEL= 1, NSLINE
              JL= MOD(IL+IJDEL-1,NSLINE) + 1
              NOJL=NRSEG(JL)
              IF( NOJL .LT. 0 ) THEN
                JSODD=MOD(JSODD+1,2)
                IF( JSODD .EQ.0) THEN
                  DPR(NOJL,NOIL,ISBG)= DPR(NOJL,NOIL,ISBG)
     >                               + TTOT * STAYIN(JL)
                  TTOT= TTOT * GOSOUT(JL)
                  IF( TTOT.LE.CUTOF ) GO TO 55
                ENDIF
              ELSE
                DPR(NOJL,NOIL,ISBG)= DPR(NOJL,NOIL,ISBG)
     >                             + TTOT* STAYIN(JL)
                TTOT= TTOT * GOSOUT(JL)
                IF( TTOT.LE.CUTOF ) GO TO 55
              ENDIF
 80         CONTINUE
          ENDIF
 55       CONTINUE
 50     CONTINUE
 2001 CONTINUE
      RETURN
      END
