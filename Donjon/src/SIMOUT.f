*DECK SIMOUT
      SUBROUTINE SIMOUT(IPMAP,IMPX,BURNINS,IZONE,NCH,NB,LX,LY,HHX,IHY,
     > STATE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Print burnup distribution (3D), radial averages or axial averages
*
*Copyright:
* Copyright (C) 2013 Ecole Polytechnique de Montreal
*
*Author(s): 
* V. Salino
*
*Parameters: input
* IPMAP   fuel map object.
* IMPX    print parameter.
* BURNINS instantaneous burnups.
* IZONE   default assembly or quart-of-assembly names as defined in
*         the fuel map.
* NCH     number of assemblies or number of quart-of-assemblies.
* NB      number of axial burnup subdivisions in an assembly.
* LX      number of assemblies along the X axis.
* LY      number of assemblies along the Y axis.
* LXMIN   coordinates on X axis of the first assembly.
* LYMIN   coordinates on Y axis of the first assembly.
* HHX     naval battle indices along X axis.
* IHY     naval battle indices along Y axis.
* STATE   flag indicating whether it is a beginning-of-stage print
*         or a end-of-stage print.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAP
      INTEGER IMPX,IHY(LY),NCH,NB,LX,LY
      CHARACTER HHX(LX)*1,IZONE(NCH)*4,STATE*5
      REAL BURNINS(NCH,NB)
*----
*  LOCAL VARIABLES
*----
      INTEGER INTG2,INTG2B
      REAL MEANR
      CHARACTER TEXT4*4,TEXT1*1,TEXT1B*1
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: MEAN
*
      IF(STATE.EQ.'BEGIN')THEN
        CALL LCMGET(IPMAP,'BURN-INST',BURNINS)
      ENDIF
*----
*  RADIALLY-AVERAGED BURNUP MAP
*----
      IF((STATE.EQ.'BEGIN'.AND.IMPX.GE.8).OR.
     >   (STATE.EQ.'END  '.AND.IMPX.GE.3)) THEN
        IF(STATE.EQ.'BEGIN')THEN
          WRITE(6,100)
        ELSE
          WRITE(6,105)
        ENDIF
        WRITE(6,110) (HHX(I),I=1,LX)
        ICH=1
        DO I=1,LY
          TEXT4=IZONE(ICH)
          READ(TEXT4,'(A1,I2)') TEXT1,INTG2
          DO J=1,LX+1
            NFULL=J
            IF(ICH.EQ.(NCH+1))GO TO 10
            TEXT4=IZONE(ICH)
            READ(TEXT4,'(A1,I2)') TEXT1B,INTG2B
            IF(INTG2.NE.INTG2B)GO TO 10
            ICH=ICH+1
          ENDDO
          CALL XABORT('@SIMOUT: INCOHERENCE IN BASIC ASSEMBLY '
     >    //'LAYOUT GIVEN IN RESINI:.')
  10      CONTINUE
          NFULL=NFULL-1
          NEMPTY=(LX-NFULL)/2
          ALLOCATE(MEAN(NFULL))
          CALL XDRSET(MEAN,NFULL,0.0)
          DO K=1,NFULL
            DO IB=1,NB
              MEAN(K)=MEAN(K)+BURNINS(ICH-1-NFULL+K,IB)/NB
            ENDDO
          ENDDO
          WRITE(6,115,ADVANCE='NO') IHY(I)
          DO K=1,NEMPTY
            WRITE(6,120,ADVANCE='NO')
          ENDDO
          WRITE(6,125) (NINT(MEAN(K)),K=1,NFULL)
          DEALLOCATE(MEAN)
        ENDDO
      ENDIF
*----
*  AXIALLY-AVERAGED BURNUP MAP
*----
      IF((STATE.EQ.'BEGIN'.AND.IMPX.GE.9).OR.
     >   (STATE.EQ.'END  '.AND.IMPX.GE.4))THEN
        IF(STATE.EQ.'BEGIN')THEN
          WRITE(6,130)
        ELSE
          WRITE(6,135)
        ENDIF
        DO IB=1,NB
          MEANR=0.0
          DO ICH=1,NCH
            MEANR=MEANR+BURNINS(ICH,IB)/NCH
          ENDDO
          WRITE(6,140) NINT(MEANR)
        ENDDO
      ENDIF
*----
*  PER-ASSEMBLY 3D BURNUP MAP
*----
      IF((STATE.EQ.'BEGIN'.AND.IMPX.GE.10).OR.
     >   (STATE.EQ.'END  '.AND.IMPX.GE.5))THEN
        IF(STATE.EQ.'BEGIN')THEN
          WRITE(6,150)
        ELSE
          WRITE(6,155)
        ENDIF
        DO ICH=1,NCH
          WRITE(6,160) IZONE(ICH)
          WRITE(6,170) (BURNINS(ICH,IB),IB=1,NB)
        ENDDO
      ENDIF
*
      IF(STATE.EQ.'BEGIN') CALL XDRSET(BURNINS,NCH*NB,0.0)
      RETURN
*
  100 FORMAT(' SIM: BEGINNING-OF-STAGE BURNUP MAP (MW*D/TONNE), ',
     > 'RADIAL VIEW :')
  105 FORMAT(' SIM: END-OF-STAGE BURNUP MAP (MW*D/TONNE), ',
     > 'RADIAL VIEW :')
  110 FORMAT(1X,20(5X,1A1))
  115 FORMAT(1X,I2)
  120 FORMAT(6X)
  125 FORMAT(21I6)
  130 FORMAT(/,' SIM: BEGINNING-OF-STAGE BURNUP MAP (MW*D/TONNE), ',
     > 'AXIAL VIEW :')
  135 FORMAT(/,' SIM: END-OF-STAGE BURNUP MAP (MW*D/TONNE), ',
     > 'AXIAL VIEW :')
  140 FORMAT(1X,I5.1)
  150 FORMAT(/,' SIM: BEGINNING-OF-STAGE 3D BURNUP MAP (MW*D/TONNE) :')
  155 FORMAT(/,' SIM: END-OF-STAGE 3D BURNUP MAP (MW*D/TONNE) :')
  160 FORMAT('   Assembly ',A)
  170 FORMAT(3X,16(1X,F7.1))
      END
