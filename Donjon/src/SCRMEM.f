*DECK SCRMEM
      SUBROUTINE SCRMEM(IPSAP,IPMEM,NCAL,NMIL,NMIX,TERP,MIXC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Copy a Saphyb into memory taking care to keep only required
* calculations and mixtures.
*
*Copyright:
* Copyright (C) 2013 Ecole Polytechnique de Montreal
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* IPSAP   address of the Saphyb object.
* IPMEM   address of the simplified Saphyb in memory created by SCRMEM.
* NCAL    number of elementary calculations in the Saphyb.
* NMIL    number of material mixtures in the Saphyb
* NMIX    maximum number of material mixtures in the microlib.
* TERP    interpolation factors.
* MIXC    mixture index in the Saphyb corresponding to each microlib
*         mixture.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPSAP,IPMEM
      INTEGER NCAL,NMIL,NMIX,MIXC(NMIX)
      REAL TERP(NCAL,NMIX)
*----
*  LOCAL VARIABLES
*----
      INTEGER DIMSAP(50)
      INTEGER IBM, IBMOLD, ICAL, ILONG, ITYLCM
      CHARACTER SIGN*12,TEXT12*12
      TYPE(C_PTR) JPSAP,KPSAP,JPMEM1,JPMEM2,KPMEM1,KPMEM2
*
      CALL LCMOP(IPMEM,'*tempSaphyb*',0,1,0)
      CALL LCMGTC(IPSAP,'SIGNATURE',12,1,SIGN)
      CALL LCMPTC(IPMEM,'SIGNATURE',12,1,SIGN)
      CALL LCMGET(IPSAP,'DIMSAP',DIMSAP)
      CALL LCMPUT(IPMEM,'DIMSAP',50,1,DIMSAP)
      JPSAP=LCMGID(IPSAP,'constphysiq')
      JPMEM1=LCMDID(IPMEM,'constphysiq')
      CALL LCMEQU(JPSAP,JPMEM1)
      JPSAP=LCMGID(IPSAP,'contenu')
      JPMEM1=LCMDID(IPMEM,'contenu')
      CALL LCMEQU(JPSAP,JPMEM1)
      JPSAP=LCMGID(IPSAP,'adresses')
      JPMEM1=LCMDID(IPMEM,'adresses')
      CALL LCMEQU(JPSAP,JPMEM1)
      JPSAP=LCMGID(IPSAP,'geom')
      JPMEM1=LCMDID(IPMEM,'geom')
      CALL LCMEQU(JPSAP,JPMEM1)
      JPMEM1=LCMLID(IPMEM,'calc',NCAL)
      DO 30 ICAL=1,NCAL
        DO IBM=1,NMIX
          IF(TERP(ICAL,IBM).NE.0.0) GO TO 10
        ENDDO
        GO TO 30
   10   WRITE(TEXT12,'(4Hcalc,I8)') ICAL
        JPSAP=LCMGID(IPSAP,TEXT12)
        JPMEM2=LCMDIL(JPMEM1,ICAL)
        KPSAP=LCMGID(JPSAP,'info')
        KPMEM1=LCMDID(JPMEM2,'info')
        CALL LCMEQU(KPSAP,KPMEM1)
        KPSAP=LCMGID(JPSAP,'divers')
        KPMEM1=LCMDID(JPMEM2,'divers')
        CALL LCMEQU(KPSAP,KPMEM1)
        CALL LCMLEN(JPSAP,'outflx',ILONG,ITYLCM)
        IF(ILONG.NE.0) THEN
          KPSAP=LCMGID(JPSAP,'outflx')
          KPMEM1=LCMDID(JPMEM2,'outflx')
          CALL LCMEQU(KPSAP,KPMEM1)
        ENDIF
        KPMEM1=LCMLID(JPMEM2,'mili',NMIL)
        DO IBMOLD=1,NMIL
          DO IBM=1,NMIX
           IF((TERP(ICAL,IBM).NE.0.).AND.(MIXC(IBM).EQ.IBMOLD)) GO TO 20
          ENDDO
          CYCLE
   20     WRITE(TEXT12,'(4Hmili,I8)') IBMOLD
          KPSAP=LCMGID(JPSAP,TEXT12)
          KPMEM2=LCMDIL(KPMEM1,IBMOLD)
          CALL LCMEQU(KPSAP,KPMEM2)
        ENDDO
   30 CONTINUE
      RETURN
      END
