*DECK PIJI2D
      SUBROUTINE PIJI2D(NREG,NSOUT,NSLINE,NCOR,NSBG,
     >                  SWVOID,SIGTAL,NPSYS,WEIGHT,
     >                  SEGLEN,NRSEG,SEGPAT,DPR,
     >                  MKI0,BIN0,PAS0,L0,
     >                  MKI1,BIN1,PAS1,XLM1,L1,
     >                  MKI2,BIN2,PAS2,XLM2)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Integration for general 2D isotropic tracking.
*
*Copyright:
* Copyright (C) 1991 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* NREG    total number of regions.
* NSOUT   number of outer surface.
* NSLINE  number of segemnts on line.
* NCOR    maximum number of corners.
* NSBG    number of subgroup.
* SWVOID  flag to indicate if there are voids.
* SIGTAL  albedo-cross section vector.
* NPSYS   non-converged energy group indices.
* WEIGHT  line weight.
* SEGLEN  length of track.
* NRSEG   region crossed by track.
* MKI0    nb element quadratic BICKLEY table order N.
* BIN0    elements quadratic BICKLEY table order N.
* PAS0    step quadratic BICKLEY table order N.
* L0      log divergence quadratic BICKLEY table order N.
* MKI1    nb element quadratic BICKLEY table order N+1.
* BIN1    elements quadratic BICKLEY table order N+1.
* PAS1    step quadratic BICKLEY table order N+1.
* XLM1    upper limit quadratic BICKLEY table order N+1.
* L1      log divergence quadratic BICKLEY table order N+1. 
* MKI2    nb element quadratic BICKLEY table order N+2.
* BIN2    elements quadratic BICKLEY table order N+2.
* PAS2    step quadratic BICKLEY table order N+2.
* XLM2    upper limit quadratic BICKLEY table order N+2.
*
*Parameters: output
* DPR     CP matrix. 
*
*Parameters: scratch
* SEGPAT  optical path.
*
*Comments:
*  PIJ  => WITH BICKLEY FUNCTIONS OF ORDER 1,2,3
*  PIJK => WITH BICKLEY FUNCTIONS OF ORDER 3,4,5 
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
* PARAMETERS
*----
      REAL             PI
      PARAMETER       (PI=3.1415926535897932)
*----
* INTERNAL FUNCTIONS
*----
      DOUBLE PRECISION CBIN0,CBIN1,CBIN2
*----
* VARIABLES
*----
      INTEGER          NREG,NSOUT,NSLINE,NCOR,NSBG
      INTEGER          NRSEG(NSLINE),NPSYS(NSBG)
      LOGICAL          SWVOID
      REAL             SIGTAL(-NSOUT:NREG,NSBG),SEGPAT(NSLINE)
      DOUBLE PRECISION WEIGHT,SEGLEN(NSLINE)
      DOUBLE PRECISION DPR(-NSOUT:NREG,-NSOUT:NREG,NSBG)
      INTEGER          MKI0,L0,MKI1,L1,MKI2
      REAL             BIN0(0:MKI0,3),PAS0,
     >                 BIN1(0:MKI1,3),PAS1,XLM1,
     >                 BIN2(0:MKI2,3),PAS2,XLM2
*----
*  Local variables
*----
      INTEGER          IL,ISBG,ISD,ISF,NSEG,ISEG,IREG,
     >                 LPOS,ICSEG,JCSEG,JSEG,IR1
      DOUBLE PRECISION FLOG0,FLOG1,ZREF1,ZREF2,ZREF3,XPOS,
     >                 ZINTP,X3,ZI,X
*----
*  QUADRATIC INTERPOLATION IN BICKLEY TABLES
*----
      CBIN0(IL,X)=       BIN0(IL,1)+X*(BIN0(IL,2)+X*BIN0(IL,3))
      CBIN1(IL,X)=       BIN1(IL,1)+X*(BIN1(IL,2)+X*BIN1(IL,3))
      CBIN2(IL,X)=       BIN2(IL,1)+X*(BIN2(IL,2)+X*BIN2(IL,3))
*----
*  LOGARITHMIC DIVERGENCE FOR KI1 AND KI2
*  FOR KI1 -> X*LOG(X)        FOR K<L0 AND L0>0
*  FOR KI2 -> -X**2*LOG(X)/2. FOR K<L1 AND L1>0
*----
      FLOG0=FLOAT(L0/MAX(1,L0))
      FLOG1=-0.5*FLOAT(L1/MAX(1,L1))
*----
*  Process track required
*----
      NSEG=NSLINE-2*NCOR
      ZI=0.0D0
      DO 2001 ISBG=1,NSBG
        IF(NPSYS(ISBG).EQ.0) GO TO 2001
        ZREF3=CBIN2(0,0.0D0)*WEIGHT
        ZREF2=CBIN1(0,0.0D0)*WEIGHT
        ZREF1=CBIN0(0,0.0D0)*WEIGHT
*----
*  CONVERT SEGMENT LENGHT TO PATH LENGTH
*----
        ICSEG=NCOR
        DO 210 ISEG=1,NSEG
          ICSEG=ICSEG+1
          IREG=NRSEG(ICSEG)
          SEGPAT(ISEG)=REAL(SEGLEN(ICSEG)*SIGTAL(IREG,ISBG))
          DPR(IREG,IREG,ISBG)=DPR(IREG,IREG,ISBG)+ZREF2*SEGPAT(ISEG)
 210    CONTINUE
*----
*  INTEGRATION
*----
        XPOS=0.0D0
        ZINTP=ZREF3
*----
*  INTEGRATE OVER FIRST REGION AND COMPUTE PVS, PSS
*----
        ICSEG=NCOR+1
        DO 220 ISEG=1,NSEG
          IREG=NRSEG(ICSEG)
          IR1=NRSEG(NCOR+1)
          XPOS=XPOS+SEGPAT(ISEG)
          LPOS=MIN(NINT(XPOS*PAS2),MKI2)
          ZI=WEIGHT*CBIN2(LPOS,XPOS)
          X3=ZI-ZINTP
          DPR(IR1,IREG,ISBG)=DPR(IR1,IREG,ISBG)+X3
          DO 225 ISD=1,NCOR
            DPR(IREG,NRSEG(ISD),ISBG)=DPR(IREG,NRSEG(ISD),ISBG)-X3
 225      CONTINUE
          IF(XPOS.GT.XLM2) GO TO 221
          ZINTP=ZI
          ICSEG=ICSEG+1
 220    CONTINUE
 221    CONTINUE
        IR1=NRSEG(NCOR+1)
        DO 226 ISF=NSLINE-NCOR+1,NSLINE
          DPR(IR1,NRSEG(ISF),ISBG)=DPR(IR1,NRSEG(ISF),ISBG)-ZI
          DO 227 ISD=1,NCOR
            DPR(NRSEG(ISD),NRSEG(ISF),ISBG)=
     >      DPR(NRSEG(ISD),NRSEG(ISF),ISBG)+ZI
 227      CONTINUE
 226    CONTINUE
        ICSEG=NCOR+2
        DO 230 ISEG=2,NSEG
          XPOS=0.0D0
          ZINTP=ZREF3
          JCSEG=ICSEG
          DO 240 JSEG=ISEG,NSEG
            XPOS=XPOS+SEGPAT(JSEG)
            LPOS=MIN(NINT(XPOS*PAS2),MKI2)
            ZI=WEIGHT*CBIN2(LPOS,XPOS)
            X3=ZI-ZINTP
            DPR(NRSEG(ICSEG-1),NRSEG(JCSEG),ISBG)=
     >          DPR(NRSEG(ICSEG-1),NRSEG(JCSEG),ISBG)-X3
            DPR(NRSEG(ICSEG),NRSEG(JCSEG),ISBG)=
     >          DPR(NRSEG(ICSEG),NRSEG(JCSEG),ISBG)+X3
            IF(XPOS.GT.XLM2) GO TO 241
            ZINTP=ZI
            JCSEG=JCSEG+1
 240      CONTINUE
 241      CONTINUE
          DO 235 ISF=NSLINE-NCOR+1,NSLINE
            DPR(NRSEG(ICSEG),NRSEG(ISF),ISBG)=
     >            DPR(NRSEG(ICSEG),NRSEG(ISF),ISBG)-ZI
            DPR(NRSEG(ICSEG-1),NRSEG(ISF),ISBG)=
     >            DPR(NRSEG(ICSEG-1),NRSEG(ISF),ISBG)+ZI
 235      CONTINUE
          ICSEG=ICSEG+1
 230    CONTINUE
        ICSEG=NCOR+NSEG
        DO 236 ISF=NSLINE-NCOR+1,NSLINE
          DPR(NRSEG(ICSEG),NRSEG(ISF),ISBG)=
     >          DPR(NRSEG(ICSEG),NRSEG(ISF),ISBG)+ZREF3
 236    CONTINUE
*----
*  FOR VOID REGIONS RESET PROBABILITIES
*----
        IF(SWVOID) THEN
          ICSEG=NCOR+1
          DO 300 ISEG=1,NSEG
            IF(SIGTAL(NRSEG(ICSEG),ISBG).NE.0.0) GO TO 301
            XPOS=0.0D0
            DPR(NRSEG(ICSEG),NRSEG(ICSEG),ISBG)=
     >          DPR(NRSEG(ICSEG),NRSEG(ICSEG),ISBG)
     >         + 0.5*ZREF1*SEGLEN(ICSEG)*SEGLEN(ICSEG)
            ZINTP=ZREF2*SEGLEN(ICSEG)
            JCSEG=ICSEG+1
            DO 310 JSEG=ISEG+1,NSEG
              IF(SIGTAL(NRSEG(JCSEG),ISBG).EQ.0.0) THEN
                LPOS=MIN(NINT(XPOS*PAS0),MKI0)
                IF(LPOS.LT.L0.AND.XPOS.NE.0.0) THEN
                  ZI=WEIGHT*SEGLEN(ICSEG)*SEGLEN(JCSEG)*
     >              (CBIN0(LPOS,XPOS)+FLOG0*XPOS*LOG(XPOS))
                ELSE
                  ZI=WEIGHT*SEGLEN(ICSEG)*SEGLEN(JCSEG)*CBIN0(LPOS,XPOS)
                ENDIF
                X3=0.5*ZI
              ELSE
                XPOS=XPOS+SEGPAT(JSEG)
                LPOS=MIN(NINT(XPOS*PAS1),MKI1)
                IF(LPOS.LT.L1.AND.XPOS.NE.0.0) THEN
                  ZI=WEIGHT*SEGLEN(ICSEG)*
     >              (CBIN1(LPOS,XPOS)+FLOG1*XPOS*XPOS*LOG(XPOS))
                ELSE
                  ZI=WEIGHT*SEGLEN(ICSEG)*CBIN1(LPOS,XPOS)
                ENDIF
                X3=ZINTP-ZI
                ZINTP=ZI
              ENDIF
              DPR(NRSEG(ICSEG),NRSEG(JCSEG),ISBG)=
     >            DPR(NRSEG(ICSEG),NRSEG(JCSEG),ISBG)+X3
              IF(XPOS.GT.XLM1) GO TO 311
              JCSEG=JCSEG+1
 310        CONTINUE
 311        CONTINUE
            DO 320 ISF=NSLINE-NCOR+1,NSLINE
              DPR(NRSEG(ICSEG),NRSEG(ISF),ISBG)=
     >              DPR(NRSEG(ICSEG),NRSEG(ISF),ISBG)+ZINTP
 320        CONTINUE
            XPOS=0.0D0
            ZINTP=ZREF2*SEGLEN(ICSEG)
            JCSEG=ICSEG-1
            DO 330 JSEG=ISEG-1,1,-1
              IF(SIGTAL(NRSEG(JCSEG),ISBG).EQ.0.0) THEN
                LPOS=MIN(NINT(XPOS*PAS0),MKI0)
                IF(LPOS.LT.L0.AND.XPOS.NE.0.0) THEN
                  ZI=WEIGHT*SEGLEN(ICSEG)*SEGLEN(JCSEG)*
     >              (CBIN0(LPOS,XPOS)+FLOG0*XPOS*LOG(XPOS))
                ELSE
                  ZI=WEIGHT*SEGLEN(ICSEG)*SEGLEN(JCSEG)*CBIN0(LPOS,XPOS)
                ENDIF
                X3=0.5*ZI
              ELSE
                XPOS=XPOS+SEGPAT(JSEG)
                LPOS=MIN(NINT(XPOS*PAS1),MKI1)
                IF(LPOS.LT.L1.AND.XPOS.NE.0.0) THEN
                  ZI=WEIGHT*SEGLEN(ICSEG)*
     >              (CBIN1(LPOS,XPOS)+FLOG1*XPOS*XPOS*LOG(XPOS))
                ELSE
                  ZI=WEIGHT*SEGLEN(ICSEG)*CBIN1(LPOS,XPOS)
                ENDIF
                X3=ZINTP-ZI
                ZINTP=ZI
              ENDIF
              DPR(NRSEG(ICSEG),NRSEG(JCSEG),ISBG)=
     >            DPR(NRSEG(ICSEG),NRSEG(JCSEG),ISBG)+X3
              IF(XPOS.GT.XLM1) GO TO 331
              JCSEG=JCSEG-1
 330        CONTINUE
 331        CONTINUE
            DO 340 ISD=1,NCOR
              DPR(NRSEG(ICSEG),NRSEG(ISD),ISBG)=
     >              DPR(NRSEG(ICSEG),NRSEG(ISD),ISBG)+ZINTP
 340        CONTINUE
 301        CONTINUE
            ICSEG=ICSEG+1
 300      CONTINUE
        ENDIF
 2001 CONTINUE
*----
*  RETURN
*----
      RETURN
      END
