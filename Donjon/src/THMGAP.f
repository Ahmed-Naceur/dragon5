*DECK THMGAP
      SUBROUTINE THMGAP(POWLIN,BURN,HGAP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the heat exchange coefficient of the gap.
*
*Copyright:
* Copyright (C) 2012 Ecole Polytechnique de Montreal.
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* POWLIN  linear power in W/m
* BURN    fuel burnup in MWday/tonne
*
*Parameters: output
* HGAP    heat exchange coefficient of the gap in W/m^2/K. Values with
*         POWLIN greater than 400 W/cm or BURN greater than 50000
*         MWday/ton and up to 90000 MWday/ton are extrapolated.
*         After 90000 MWday/ton, the setting of a constant HGAP value
*         is required and the thermal mechanic model below is by-passed. 
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      REAL POWLIN,BURN,HGAP
*----
*  LOCAL VARIABLES
*----
      CHARACTER HSMG*300
      REAL TAB1(19),TAB2(11),C(19,11),TERP1(19),TERP2(11),WK1(3,19),
     > WK2(3,11)
      INTEGER I1,I2
*
      DATA TAB1/0.,5000.,10000.,15000.,20000.,25000.,30000.,35000.,
     >          40000.,45000.,50000.,55000.,60000.,65000.,70000.,
     >          75000.,80000.,85000.,90000./
      DATA TAB2/0.,30.,100.,170.,240.,310.,380.,400.,420.,440.,460./
      DATA C/0.657,0.702,0.814,0.987,1.311,2.114,2.445,2.415,2.324,2.229
     >      ,2.137,2.048,1.962,1.880,1.786,1.697,1.608,1.519,1.430
     >      ,0.678,0.726,0.848,1.043,1.444,2.465,2.810,2.790,2.718,2.640
     >      ,2.562,2.484,2.406,2.328,2.250,2.172,2.094,2.016,1.938
     >      ,0.727,0.783,0.927,1.173,1.755,3.283,3.661,3.666,3.637,3.598
     >      ,3.554,3.505,3.453,3.397,3.356,3.307,3.259,3.211,3.163
     >      ,0.787,0.854,1.032,1.373,2.322,3.800,3.790,3.780,3.769,3.756
     >      ,3.741,3.724,3.706,3.687,3.673,3.656,3.640,3.623,3.607
     >      ,0.861,0.949,1.185,1.725,3.385,3.873,3.863,3.854,3.842,3.829
     >      ,3.814,3.797,3.779,3.760,3.746,3.729,3.713,3.696,3.680
     >      ,0.949,1.068,1.415,2.385,3.925,3.910,3.900,3.891,3.879,3.865
     >      ,3.850,3.834,3.817,3.800,3.785,3.769,3.754,3.738,3.722
     >      ,1.071,1.248,1.843,3.686,3.957,3.941,3.929,3.915,3.898,3.875
     >      ,3.847,3.814,3.779,3.742,3.711,3.678,3.644,3.611,3.578
     >      ,1.114,1.317,2.033,3.981,3.964,3.946,3.931,3.911,3.885,3.851
     >      ,3.807,3.754,3.697,3.638,3.589,3.535,3.481,3.428,3.374
     >      ,1.161,1.396,2.264,4.153,4.002,3.950,3.926,3.897,3.856,3.804
     >      ,3.735,3.651,3.560,3.469,3.390,3.306,3.221,3.137,3.052
     >      ,1.212,1.485,2.542,4.155,4.090,3.953,3.913,3.869,3.806,3.729
     >      ,3.624,3.495,3.356,3.219,3.098,2.969,2.841,2.712,2.583
     >      ,1.268,1.586,2.873,3.938,4.243,3.956,3.889,3.826,3.731,3.620
     >      ,3.465,3.273,3.067,2.867,2.687,2.497,2.306,2.116,1.926/
*
      IF(BURN.GT.90000.) THEN
        WRITE(HSMG,'(22HTHMGAP: BURNUP VALUE (,1P,E11.4,
     >      35H) TOO HIGH FOR THE THERMAL MECHANIC,
     >      41H MODEL COMPUTING THE HEAT EXCHANGE OF THE,
     >      38H FUEL-CLADDING GAP (LIMIT 90000MWd/t).,
     >      45H ALTERNATIVELY, YOU CAN SET THE HGAP CONSTANT,
     >      19H IN THE THM MODULE.)') BURN
        CALL XABORT(HSMG)
      ENDIF

      CALL ALTERP(.TRUE.,19,TAB1,BURN,.FALSE.,TERP1,WK1)
      HGAP=0.0
      IF(POWLIN.LE.460.E2) THEN
        CALL ALTERP(.TRUE.,11,TAB2,POWLIN/1.0E2,.FALSE.,TERP2,WK2)
        DO I1=1,19
          DO I2=1,11
            HGAP=HGAP+TERP1(I1)*TERP2(I2)*C(I1,I2)
          ENDDO
        ENDDO
      ELSE
        DO I1=1,19
          HGAP=HGAP+TERP1(I1)*C(I1,11)
        ENDDO
      ENDIF
      HGAP=HGAP*1.0E4
      RETURN
      END
