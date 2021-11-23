*DECK DREJ02
      SUBROUTINE DREJ02(ITYPE,IELEM,NREG,NUN,MAXKN,MAXQF,MAT,KN,QFR,
     1 IQFR,VOL,FUNKNO,OUT,GAMMA)
*
*-----------------------------------------------------------------------
*
*Purpose:
* One-speed net surfacic current calculation in 3D mixed-dual finite
* element approximation.
*
*Copyright:
* Copyright (C) 2018 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* ITYPE   type of geometry: =5/7 (Cartesian); =8/9 (hexagonal).
* IELEM   degree of the Lagrangian finite elements: =1 (linear);
*          =2 (parabolic); =3 (cubic); =4 (quartic).
* NREG    number of elements.
* NUN     dimension of array FUNKNO.
* MAXKN   dimension of array KN.
* MAXQF   dimension of array QFR.
* MAT     mixture index per region.
* KN      element-ordered unknown list.
* QFR     element-ordered surfaces.
* IQFR    element-ordered physical albedo indices.
* VOL     volume of regions.
* FUNKNO  neutron fluxes.
*
*Parameters: output
* OUT     net surfacic current.
* GAMMA   gamma function.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER ITYPE,IELEM,NREG,NUN,MAXKN,MAXQF,MAT(NREG),KN(MAXKN),
     1 IQFR(MAXQF)
      REAL QFR(MAXQF),VOL(NREG),FUNKNO(NUN),OUT,GAMMA(NUN)
*
      IF(IELEM.LT.0) CALL XABORT('DREJ02: TYPE OF DISCRETIZATION NOT I'
     1 //'MPLEMENTED.')
      CALL XDRSET(GAMMA,NUN,0.0)
      OUT=0.0
      NUM1=0
      NUM2=0
      IF((ITYPE.NE.5).OR.(ITYPE.NE.7)) THEN
        DO 20 K=1,NREG
        IF(MAT(K).EQ.0) GO TO 20
        IF(VOL(K).EQ.0.0) GO TO 10
        IF(NUM2+4.GT.MAXQF) call XABORT('overflow')
        IALB=IQFR(NUM2+1)
        IF((IALB.NE.0).AND.(QFR(NUM2+1).GT.0.0)) THEN
          IND1=KN(NUM1+2)
          IF(IND1.LE.0) CALL XABORT('DREJ02: bug1')
          GAMMA(IND1)=-QFR(NUM2+1)
          OUT=OUT-FUNKNO(IND1)*QFR(NUM2+1)
        ENDIF
        IALB=IQFR(NUM2+2)
        IF((IALB.NE.0).AND.(QFR(NUM2+2).GT.0.0)) THEN
          IND1=KN(NUM1+2+IELEM**2)
          IF(IND1.LE.0) CALL XABORT('DREJ02: bug2')
          GAMMA(IND1)=QFR(NUM2+2)
          OUT=OUT+FUNKNO(IND1)*QFR(NUM2+2)
        ENDIF
        IALB=IQFR(NUM2+3)
        IF((IALB.NE.0).AND.(QFR(NUM2+3).GT.0.0)) THEN
          IND1=KN(NUM1+2+2*IELEM**2)
          IF(IND1.LE.0) CALL XABORT('DREJ02: bug3')
          GAMMA(IND1)=-QFR(NUM2+3)
          OUT=OUT-FUNKNO(IND1)*QFR(NUM2+3)
        ENDIF
        IALB=IQFR(NUM2+4)
        IF((IALB.NE.0).AND.(QFR(NUM2+4).GT.0.0)) THEN
          IND1=KN(NUM1+2+3*IELEM**2)
          IF(IND1.LE.0) CALL XABORT('DREJ02: bug4')
          GAMMA(IND1)=QFR(NUM2+4)
          OUT=OUT+FUNKNO(IND1)*QFR(NUM2+4)
        ENDIF
        IALB=IQFR(NUM2+5)
        IF((IALB.NE.0).AND.(QFR(NUM2+5).GT.0.0)) THEN
          IND1=KN(NUM1+2+4*IELEM**2)
          IF(IND1.LE.0) CALL XABORT('DREJ02: bug5')
          GAMMA(IND1)=-QFR(NUM2+5)
          OUT=OUT-FUNKNO(IND1)*QFR(NUM2+5)
        ENDIF
        IALB=IQFR(NUM2+6)
        IF((IALB.NE.0).AND.(QFR(NUM2+6).GT.0.0)) THEN
          IND1=KN(NUM1+2+5*IELEM**2)
          IF(IND1.LE.0) CALL XABORT('DREJ02: bug6')
          GAMMA(IND1)=QFR(NUM2+6)
          OUT=OUT+FUNKNO(IND1)*QFR(NUM2+6)
        ENDIF
   10   NUM1=NUM1+1+6*IELEM**2
        NUM2=NUM2+6
   20   CONTINUE
      ELSE IF((ITYPE.NE.8).OR.(ITYPE.NE.9)) THEN
        DO 40 K=1,NREG
        IF(MAT(K).EQ.0) GO TO 40
        IF(VOL(K).EQ.0.0) GO TO 30
        IF((IALB.NE.0).AND.(QFR(NUM2+1).GT.0.0)) THEN
          IND1=KN(NUM1+2)
          IF(IND1.LE.0) CALL XABORT('DREJ02: bug7')
          GAMMA(IND1)=-QFR(NUM2+1)
          OUT=OUT-FUNKNO(IND1)*QFR(NUM2+1)
        ENDIF
        IALB=IQFR(NUM2+2)
        IF((IALB.NE.0).AND.(QFR(NUM2+2).GT.0.0)) THEN
          IND1=KN(NUM1+2+IELEM**2)
          IF(IND1.LE.0) CALL XABORT('DREJ02: bug8')
          GAMMA(IND1)=QFR(NUM2+2)
          OUT=OUT+FUNKNO(IND1)*QFR(NUM2+2)
        ENDIF
        IALB=IQFR(NUM2+3)
        IF((IALB.NE.0).AND.(QFR(NUM2+3).GT.0.0)) THEN
          IND1=KN(NUM1+2+2*IELEM**2)
          IF(IND1.LE.0) CALL XABORT('DREJ02: bug9')
          GAMMA(IND1)=-QFR(NUM2+3)
          OUT=OUT-FUNKNO(IND1)*QFR(NUM2+3)
        ENDIF
        IALB=IQFR(NUM2+4)
        IF((IALB.NE.0).AND.(QFR(NUM2+4).GT.0.0)) THEN
          IND1=KN(NUM1+2+3*IELEM**2)
          IF(IND1.LE.0) CALL XABORT('DREJ02: bug10')
          GAMMA(IND1)=QFR(NUM2+4)
          OUT=OUT+FUNKNO(IND1)*QFR(NUM2+4)
        ENDIF
        IALB=IQFR(NUM2+5)
        IF((IALB.NE.0).AND.(QFR(NUM2+5).GT.0.0)) THEN
          IND1=KN(NUM1+2+4*IELEM**2)
          IF(IND1.LE.0) CALL XABORT('DREJ02: bug11')
          GAMMA(IND1)=-QFR(NUM2+5)
          OUT=OUT-FUNKNO(IND1)*QFR(NUM2+5)
        ENDIF
        IALB=IQFR(NUM2+6)
        IF((IALB.NE.0).AND.(QFR(NUM2+6).GT.0.0)) THEN
          IND1=KN(NUM1+2+5*IELEM**2)
          IF(IND1.LE.0) CALL XABORT('DREJ02: bug12')
          GAMMA(IND1)=QFR(NUM2+6)
          OUT=OUT+FUNKNO(IND1)*QFR(NUM2+6)
        ENDIF
        IALB=IQFR(NUM2+7)
        IF((IALB.NE.0).AND.(QFR(NUM2+7).GT.0.0)) THEN
          IND1=KN(NUM1+2+6*IELEM**2)
          IF(IND1.LE.0) CALL XABORT('DREJ02: bug13')
          GAMMA(IND1)=-QFR(NUM2+7)
          OUT=OUT-FUNKNO(IND1)*QFR(NUM2+7)
        ENDIF
        IALB=IQFR(NUM2+8)
        IF((IALB.NE.0).AND.(QFR(NUM2+8).GT.0.0)) THEN
          IND1=KN(NUM1+2+7*IELEM**2)
          IF(IND1.LE.0) CALL XABORT('DREJ02: bug14')
          GAMMA(IND1)=QFR(NUM2+8)
          OUT=OUT+FUNKNO(IND1)*QFR(NUM2+8)
        ENDIF
   30   NUM1=NUM1+1+8*IELEM**2
        NUM2=NUM2+8
   40   CONTINUE
      ELSE
        CALL XABORT('DREJ02: TYPE OF GEOMETRY NOT IMPLEMENTED.')
      ENDIF
      RETURN
      END
