*DECK DETFLU
      SUBROUTINE DETFLU(LHEX,NX,NY,NZ,NEL,NUN,MESHX,MESHY,MESHZ,KEYF,
     > FLUX,NGRP,SPEC,DEVPOS,NHEX,IHEX,RESP,IPRT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute flux at detector site
*
*Copyright:
* Copyright (C) 2010 Ecole Polytechnique de Montreal.
*
*Author(s): 
* E. Varin, M. Guyot
*
*Parameters:
* LHEX   =.TRUE. if hexagonal detectors are present
* NX     number of x mesh-splitted elements 
* NY     number of y mesh-splitted elements 
* NZ     number of z mesh-splitted elements
* NEL    number of finite elements
* NUN    number of unknowns
* MESHX  regions coordinates according to x
* MESHY  regions coordinates according to y
* MESHZ  regions coordinates according to z
* KEYF   keyflux recover from L_TRACk object
* FLUX   flux for each mesh-splitted elements
* NGRP   number of energy groups
* SPEC   spectral information
* DEVPOS detector coordinates
* NHEX   number of hexagons in the detector
* IHEX   index number of hexagons
* COR    center detector coordinates
* RESP   flux reads by the detector
* IPRT   printing index
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NX,NY,NZ,NEL,NUN,NGRP,IPRT,NHEX,NT,KEYF(NEL),IHEX(NHEX)
      REAL MESHX(NX+1),MESHY(NY+1),MESHZ(NZ+1),FLUX(NUN,NGRP),RESP,
     1     DEVPOS(6),SPEC(NGRP)
      LOGICAL LHEX
*----
*  LOCAL VARIABLES
*----
      INTEGER NXP1,NYP1,NZP1,I,J,K,I1,I2,J1,J2,K1,K2,IAM,IGR
      REAL    X1,X2,Y1,Y2,Z1,Z2

      NXP1 = NX+1
      NYP1 = NY+1
      NZP1 = NZ+1

      X1=DEVPOS(1)
      X2=DEVPOS(2)
      Y1=DEVPOS(3)
      Y2=DEVPOS(4)
      Z1=DEVPOS(5)
      Z2=DEVPOS(6)

      IF(.NOT.LHEX) THEN
        IF(X1.LT.MESHX(1)) X1=MESHX(1)
        IF(X2.LT.MESHX(1)) X2=MESHX(1)
        IF(X2.GT.MESHX(NXP1)) X2=MESHX(NXP1)
        IF(X1.GT.MESHX(NXP1)) X1=MESHX(NXP1)

        IF(Y1.LT.MESHY(1)) Y1=MESHY(1)
        IF(Y2.LT.MESHY(1)) Y2=MESHY(1)
        IF(Y2.GT.MESHY(NYP1)) Y2=MESHY(NYP1)
        IF(Y1.GT.MESHY(NYP1)) Y1=MESHY(NYP1)
      ENDIF

      IF(Z1.LT.MESHZ(1)) Z1=MESHZ(1)
      IF(Z2.LT.MESHZ(1)) Z2=MESHZ(1)
      IF(Z2.GT.MESHZ(NZP1)) Z2=MESHZ(NZP1)
      IF(Z1.GT.MESHZ(NZP1)) Z1=MESHZ(NZP1)

      IF(.NOT.LHEX) THEN
        I1=0
        DO 20 I=1,NXP1
          IF(X1.GE.MESHX(I) .AND. X1.LE.MESHX(I+1)) THEN
            I1=I
          ENDIF
          IF(X2.GE.MESHX(I) .AND. X2.LE.MESHX(I+1)) THEN
            I2=I
            GOTO 10
          ENDIF
  20    CONTINUE

  10    DO 30 J=1,NYP1
          IF(Y1.GE.MESHY(J) .AND. Y1.LE.MESHY(J+1)) THEN
              J1=J
          ENDIF
          IF(Y2.GE.MESHY(J) .AND. Y2.LE.MESHY(J+1)) THEN
            J2=J
            GOTO 40
          ENDIF
  30    CONTINUE
  40    CONTINUE
      ELSE
        J1 = 1
        J2 = 1
        I1 = 1
        I2 = NHEX
      ENDIF

      DO 50 K=1,NZP1
          IF(Z1.GE.MESHZ(K) .AND. Z1.LE.MESHZ(K+1)) THEN
             K1=K
          ENDIF
          IF(Z2.GE.MESHZ(K) .AND. Z2.LE.MESHZ(K+1)) THEN
             K2=K
             GOTO 60
          ENDIF
  50  CONTINUE

  60  RESP = 0.0
      NT = 0

      IF(IPRT.GT.4) WRITE(6,*) 'POS GEOM ',I1,I2,J1,J2,K1,K2
      DO 70 K=K1,K2
       DO 71 J=J1,J2
        DO 72 I=I1,I2
          NT = NT+1
          IF(LHEX) THEN
            IAM = (K-1)*NX+IHEX(I)
          ELSE
            IAM=(K-1)*NX*NY+(J-1)*NX+I
          ENDIF
          DO 73 IGR=1,NGRP
            RESP = RESP + SPEC(IGR)*FLUX(KEYF(IAM),IGR)
  73      CONTINUE
          IF(IPRT.GT.4) WRITE(6,*) 'DETFLU: FINITE ELEMENT NUMBER ',
     +                             IAM
  72  CONTINUE
  71  CONTINUE
  70  CONTINUE

      RESP = RESP / FLOAT(NT)

      RETURN
      END
