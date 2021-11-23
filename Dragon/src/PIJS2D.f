*DECK PIJS2D
      SUBROUTINE PIJS2D(NREG,NSOUT,NSLINE,NSBG,WEIGHT,
     >                  RCUTOF,NGSS,SIGANG,XGSS,WGSS,NPSYS,
     >                  SEGLEN,NRSEG,
     >                  STAYIN,GOSOUT,DPR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Integration for general 2D specular tracking.
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
* NGSS    number of Gauss points.
* SIGANG  3D albedo-cross section vector.
* XGSS    2D->3D conversion for integration.
* WGSS    2D->3D conversion for integration.
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
      INTEGER          NREG,NSOUT,NSLINE,NSBG,NGSS
      DOUBLE PRECISION WEIGHT,SEGLEN(NSLINE),STAYIN(NGSS,NSLINE),
     >                 GOSOUT(NGSS,NSLINE)
      REAL             RCUTOF
      INTEGER          NRSEG(NSLINE),NPSYS(NSBG)
      REAL             SIGANG(NGSS,-NSOUT:NREG,NSBG)
      DOUBLE PRECISION DPR(-NSOUT:NREG,-NSOUT:NREG,NSBG)
*----
*  Local variables
*----
      INTEGER          ISBG,IL,JL,NOIL,NOJL,ISODD,JSODD,IJDEL
      INTEGER          MXGAUS,IGSS
      PARAMETER       (MXGAUS= 64 )
      REAL             WGSS(MXGAUS),FINV(MXGAUS),XGSS(MXGAUS),
     >                 ZERO,ONE,HALF,XSIL,CUTOF
      DOUBLE PRECISION TTOT(MXGAUS),OPATH
      PARAMETER       (ZERO=0.0E0, ONE=1.0E0, HALF=0.5E0 )
      REAL             SIXT,CUTEXP
      PARAMETER       (SIXT=HALF/3.0,CUTEXP=0.02)
      DOUBLE PRECISION EXSIL,XSIL2
      DO 2001 ISBG=1,NSBG
        IF(NPSYS(ISBG).EQ.0) GO TO 2001
        DO 20 IGSS= 1,NGSS
           TTOT(IGSS)= ONE
 20     CONTINUE
*
*1.1)    CHANGE PATHS => GOSOUT AND STAYIN PATHS, INCLUDING ALBEDOS
*        ADD *PII* LOCAL NON-CYCLIC CONTRIBUTIONS
        ISODD=0
        DO 30 IL= 1, NSLINE
          NOIL  = NRSEG(IL)
          IF( NOIL .LT. 0 )THEN
            IF(ISODD .EQ. 1) THEN
              ISODD=0
              DO 31 IGSS= 1,NGSS
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
                GOSOUT(IGSS,IL)= SIGANG(1,NOIL,ISBG)
                STAYIN(IGSS,IL)= ONE - GOSOUT(IGSS,IL)
                TTOT(IGSS)= TTOT(IGSS) * GOSOUT(IGSS,IL)
 31           CONTINUE
            ELSE
              ISODD=1
              DO 32 IGSS= 1,NGSS
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
                GOSOUT(IGSS,IL)= SIGANG(1,NOIL,ISBG)
                STAYIN(IGSS,IL)= ONE
*                TTOT(IGSS)= TTOT(IGSS) * GOSOUT(IGSS,IL)
 32           CONTINUE
            ENDIF
          ELSE IF(NOIL .GT. 0) THEN
*----
*  FOR REGIONS
*  STAYIN = 1 -  EXP[ -CROSS SECTION * LENGTH OF NSLINE]
*  GOSOUT = 1 -  STAYIN
*  TTOT   = PRODUCT OF GOSOUT
*----
            XSIL  = SIGANG(1,NOIL,ISBG)
            IF( XSIL .EQ. ZERO) THEN
              DO 33 IGSS= 1,NGSS
                GOSOUT(IGSS,IL)= ONE
                STAYIN(IGSS,IL)= SEGLEN(IL)*XGSS(IGSS)
                DPR(NOIL,NOIL,ISBG)= DPR(NOIL,NOIL,ISBG)+HALF*WEIGHT*
     >                 WGSS(IGSS)*STAYIN(IGSS,IL)*STAYIN(IGSS,IL)
 33           CONTINUE
            ELSE IF( XSIL .LT. CUTEXP) THEN
              DO 333 IGSS= 1,NGSS
                OPATH= SIGANG(IGSS,NOIL,ISBG)*SEGLEN(IL)
                XSIL2=OPATH*OPATH
                EXSIL=XSIL2*(HALF-SIXT*OPATH+XSIL2/24.0)
                STAYIN(IGSS,IL)=OPATH-EXSIL
                GOSOUT(IGSS,IL)= ONE - STAYIN(IGSS,IL)
                TTOT(IGSS)= TTOT(IGSS) * GOSOUT(IGSS,IL)
                DPR(NOIL,NOIL,ISBG)= DPR(NOIL,NOIL,ISBG)
     >                + WEIGHT*WGSS(IGSS)*EXSIL
 333          CONTINUE
            ELSE
              DO 34 IGSS=1,NGSS
                OPATH= SIGANG(IGSS,NOIL,ISBG)*SEGLEN(IL)
                EXSIL= EXP(-OPATH)
                STAYIN(IGSS,IL)= ONE - EXSIL
                GOSOUT(IGSS,IL)= EXSIL
                TTOT(IGSS)= TTOT(IGSS) * GOSOUT(IGSS,IL)
                DPR(NOIL,NOIL,ISBG)= DPR(NOIL,NOIL,ISBG)
     >                + WEIGHT*WGSS(IGSS)*(OPATH-STAYIN(IGSS,IL))
 34           CONTINUE
            ENDIF
          ENDIF
 30     CONTINUE
*
*1.2)    COMPUTE CYCLIC FACTORS BY ANGLE
*        USING GLOBAL TRACK ATTENUATION: BETA(TOT)*EXP(-MFP(TOT))
        DO 40 IGSS= 1,NGSS
          IF( TTOT(IGSS).GE.ONE )THEN
            CALL XABORT( 'PIJS2D: ALBEDOS ARE NOT COMPATIBLE')
          ENDIF
          FINV(IGSS)= REAL(WEIGHT * WGSS(IGSS) / (ONE-TTOT(IGSS)))
 40     CONTINUE
*
*1.3)    ADD *PIJ* CONTRIBUTIONS FOR FORWARD SOURCES
        ISODD=0
        DO 50 IL= 1, NSLINE
          NOIL  = NRSEG(IL)
          DO 60 IGSS= 1, NGSS
            TTOT(IGSS)= FINV(IGSS) * STAYIN(IGSS,IL)
 60       CONTINUE
          CUTOF= REAL(RCUTOF*TTOT(1))
          IF( NOIL .LT. 0) THEN
            ISODD=MOD(ISODD+1,2)
            JSODD=ISODD
            DO 70 IJDEL= 1, NSLINE
              JL= MOD(IL+IJDEL-1,NSLINE) + 1
              NOJL=NRSEG(JL)
              IF( NOJL .LT. 0 ) THEN
                JSODD=MOD(JSODD+1,2)
                IF( ISODD.EQ.1 .AND. JSODD .EQ.0) THEN
                  DO 71 IGSS= 1, NGSS
                    DPR(NOJL,NOIL,ISBG)= DPR(NOJL,NOIL,ISBG)
     >                            + TTOT(IGSS) * STAYIN(IGSS,JL)
                    TTOT(IGSS)= TTOT(IGSS) * GOSOUT(IGSS,JL)
 71               CONTINUE
                  IF( TTOT(1).LE.CUTOF ) GO TO 55
                ENDIF
              ELSE IF((ISODD.EQ.1).AND.( NOJL .GT. 0 )) THEN
                DO 72 IGSS= 1, NGSS
                  DPR(NOJL,NOIL,ISBG)= DPR(NOJL,NOIL,ISBG)
     >                          + TTOT(IGSS) * STAYIN(IGSS,JL)
                  TTOT(IGSS)= TTOT(IGSS) * GOSOUT(IGSS,JL)
 72             CONTINUE
                IF( TTOT(1).LE.CUTOF ) GO TO 55
              ENDIF
 70         CONTINUE
          ELSE IF( NOIL .GT. 0) THEN
            JSODD=ISODD
            DO 80 IJDEL= 1, NSLINE
              JL= MOD(IL+IJDEL-1,NSLINE) + 1
              NOJL=NRSEG(JL)
              IF( NOJL .LT. 0 ) THEN
                JSODD=MOD(JSODD+1,2)
                IF( JSODD .EQ.0) THEN
                  DO 81 IGSS= 1, NGSS
                    DPR(NOJL,NOIL,ISBG)= DPR(NOJL,NOIL,ISBG)
     >                            + TTOT(IGSS) * STAYIN(IGSS,JL)
                    TTOT(IGSS)= TTOT(IGSS) * GOSOUT(IGSS,JL)
 81               CONTINUE
                  IF( TTOT(1).LE.CUTOF ) GO TO 55
                ENDIF
              ELSE IF( NOJL .GT. 0 ) THEN
                DO 82 IGSS= 1, NGSS
                  DPR(NOJL,NOIL,ISBG)= DPR(NOJL,NOIL,ISBG)
     >                          + TTOT(IGSS) * STAYIN(IGSS,JL)
                  TTOT(IGSS)= TTOT(IGSS) * GOSOUT(IGSS,JL)
 82             CONTINUE
                IF( TTOT(1).LE.CUTOF ) GO TO 55
              ENDIF
 80         CONTINUE
          ENDIF
 55       CONTINUE
 50     CONTINUE
 2001 CONTINUE
      RETURN
      END
