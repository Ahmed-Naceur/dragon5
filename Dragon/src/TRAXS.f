*DECK TRAXS
      SUBROUTINE TRAXS(IPMAC1,IPMAC2,NG,NMIL,NL,NF,NDEL,ISTEP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Macrolib transposition.
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPMAC1  pointer to the transposed macrolib (L_MACROLIB signature).
* IPMAC2  pointer to the original macrolib (L_MACROLIB signature).
* NG      number of energy groups.
* NMIL    number of homogenized mixtures.
* NL      number of Legendre orders.
* NF      number of fissile isotopes.
* NDEL    number of precursor groups.
* ISTEP   number of components in STEP directory.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAC1,IPMAC2
      INTEGER NG,NMIL,NL,NF,NDEL,ISTEP
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NCOPY1=4,NCOPY2=11)
      TYPE(C_PTR) JPMAC1,KPMAC1,JPMAC2,KPMAC2
      CHARACTER TEXT12*12,TCOPY1(NCOPY1)*12,TCOPY2(NCOPY2)*12
      REAL, ALLOCATABLE, DIMENSION(:) :: GAR1,XIOF
      DATA TCOPY1/'ENERGY','DELTAU','FLUXDISAFACT','DIFFB1HOM'/
      DATA TCOPY2/'ADDXSNAME-P0','FISSIONINDEX','ALBEDO','VOLUME',
     1 'LAMBDA-D','BETA-D','K-EFFECTIVE','K-INFINITY','B2  B1HOM',
     2 'B2  HETE','TIMESTAMP'/
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(GAR1(NG+1))
*----
*  COPY THE INFORMATION ON ROOT
*----
      DO ICOPY=1,NCOPY1
        TEXT12=TCOPY1(ICOPY)
        CALL LCMLEN(IPMAC2,TEXT12,ILONG,ITYLCM)
        IF(ILONG.GT.0) THEN
          CALL LCMGET(IPMAC2,TEXT12,GAR1)
          ALLOCATE(XIOF(ILONG))
          DO I=1,ILONG
            XIOF(I)=GAR1(ILONG+1-I)
          ENDDO
          CALL LCMPUT(IPMAC1,TEXT12,ILONG,2,XIOF)
          DEALLOCATE(XIOF)
        ENDIF
      ENDDO
      DO ICOPY=1,NCOPY2
        TEXT12=TCOPY2(ICOPY)
        CALL LCMLEN(IPMAC2,TEXT12,ILONG,ITYLCM)
        IF(ILONG.GT.0) THEN
          ALLOCATE(XIOF(ILONG))
          CALL LCMGET(IPMAC2,TEXT12,XIOF)
          CALL LCMPUT(IPMAC1,TEXT12,ILONG,2,XIOF)
          DEALLOCATE(XIOF)
        ENDIF
      ENDDO
*----
*  COPY THE INFORMATION ON DIRECTORY GROUP
*----
      CALL TRAGRO(IPMAC1,IPMAC2,NG,NMIL,NL,NF,NDEL)
      IF(ISTEP.GT.0) THEN
        JPMAC2=LCMGID(IPMAC2,'STEP')
        JPMAC1=LCMLID(IPMAC1,'STEP',ISTEP)
        DO IS=1,ISTEP
          KPMAC2=LCMGIL(JPMAC2,IS)
          KPMAC1=LCMDIL(JPMAC1,IS)
          CALL TRAGRO(KPMAC1,KPMAC2,NG,NMIL,NL,NF,NDEL)
        ENDDO
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(GAR1)
      RETURN
      END
