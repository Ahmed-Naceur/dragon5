*DECK TABKI
      FUNCTION TABKI(L,X)
*
*-----------------------------------------------------------------------
*
* COMPUTES BICKLEY FUNCTION FROM QUADRATIC TABLES.
*
*    L  : ORDER OF THE BICKLEY FUNCTION.
*    X  : ARGUMENT.
*
*-----------------------------------------------------------------------
*
      IMPLICIT  NONE
      INTEGER   IOUT,MKI1,MKI2,MKI3,MKI4,MKI5
      CHARACTER NAMSBR*6
      PARAMETER(IOUT=6,MKI1=600,MKI2=600,MKI3=600,MKI4=600,MKI5=600,
     >          NAMSBR='TABKI ')
C----
C  ROUTINE PARAMETERS
C-----
      INTEGER   L
      REAL      X
C----
C  FUNCTION TYPE
C----
      REAL      TABKI
C----
C  COMMON PARAMETERS
C----
      INTEGER   L1,L2,L3,L4,L5
      REAL      BI1,PAS1,XLIM1,BI2,PAS2,XLIM2,
     >          BI3,PAS3,XLIM3,BI4,PAS4,XLIM4,
     >          BI5,PAS5,XLIM5
      COMMON /BICKL1/BI1(0:MKI1,3),PAS1,XLIM1,L1
      COMMON /BICKL2/BI2(0:MKI2,3),PAS2,XLIM2,L2
      COMMON /BICKL3/BI3(0:MKI3,3),PAS3,XLIM3,L3
      COMMON /BICKL4/BI4(0:MKI4,3),PAS4,XLIM4,L4
      COMMON /BICKL5/BI5(0:MKI5,3),PAS5,XLIM5,L5
C----
C  LOCAL PARAMETERS
C----
      INTEGER   K
      REAL      Y
      IF(X .LT. 0.0) THEN
        WRITE(IOUT,9000) NAMSBR,L,X,NAMSBR,L,0.0
      ENDIF
      Y=MAX(X,0.0)
      TABKI=0.0
      IF(Y.GT.1.E5) RETURN
      IF(L .EQ. 1) THEN
        K=MIN(NINT(Y*PAS1),MKI1)
        TABKI=BI1(K,1)+Y*(BI1(K,2)+Y*BI1(K,3))
        IF(K .LT. L1 ) THEN
          IF(Y .NE. 0.) THEN
            TABKI= TABKI + Y*LOG(Y)
          ENDIF
        ENDIF
      ELSE IF(L .EQ. 2) THEN
        K=MIN(NINT(Y*PAS2),MKI2)
        TABKI=BI2(K,1)+Y*(BI2(K,2)+Y*BI2(K,3))
        IF(K .LT. L2)THEN
          IF(Y .NE. 0. ) THEN
            TABKI= TABKI - 0.5*Y*Y*LOG(Y)
          ENDIF
        ENDIF
      ELSE IF(L .EQ. 3) THEN
        K=MIN(NINT(Y*PAS3),MKI3)
        TABKI=BI3(K,1)+Y*(BI3(K,2)+Y*BI3(K,3))
      ELSE IF(L .EQ. 4) THEN
        K=MIN(NINT(Y*PAS4),MKI4)
        TABKI=BI4(K,1)+Y*(BI4(K,2)+Y*BI4(K,3))
      ELSE IF(L .EQ. 5) THEN
        K=MIN(NINT(Y*PAS5),MKI5)
        TABKI=BI5(K,1)+Y*(BI5(K,2)+Y*BI5(K,2))
      ELSE
        CALL XABORT(NAMSBR//': L > 5 AND L < 1 ARE INVALID')
      ENDIF
C----
C  FORMATS
C----
 9000 FORMAT(1X,' INVALID X IN : ',A6,'(',I1,',',E15.6,')',
     >       5X,' REPLACED BY  : ',A6,'(',I1,',',E15.6,')')
      RETURN
      END
