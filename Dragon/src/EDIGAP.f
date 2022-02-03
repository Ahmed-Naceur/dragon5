*DECK EDIGAP
      SUBROUTINE EDIGAP(IPADF,TEXT8,NGROUP,NGCOND,NREGIO,VOLUME,IGCOND,
     1 FLUXES,FLUHOM,IPRINT,COURIN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of the averaged boundary fluxes.
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
* IPADF   pointer to the LCM directory containing ADF-related options.
* TEXT8   name of the boundary fluxes used to compute ADF ('FD_B',
*         'FD_C' or 'FD_H').
* NGROUP  number of energy groups in the reference calculation.
* NGCOND  number of condensed groups.
* NREGIO  number of regions in the reference calculation.
* VOLUME  volume of regions.
* IGCOND  limit of condensed groups.
* FLUXES  heterogeneous gap fluxes.
* FLUHOM  homogeneous fluxes.
* IPRINT  print flag.
*
*Parameters: output
* COURIN  averaged boundary fluxes.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPADF
      INTEGER NGROUP,NGCOND,NREGIO,IGCOND(NGCOND),IPRINT
      CHARACTER TEXT8*(*)
      REAL VOLUME(NREGIO),FLUXES(NREGIO,NGROUP),FLUHOM(NGROUP),
     1 COURIN(NGCOND)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION SUM,SUD
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IFGAP
*
      CALL LCMLEN(IPADF,TEXT8,NGAP,ITYLCM)
      ALLOCATE(IFGAP(NGAP))
      CALL LCMGET(IPADF,TEXT8,IFGAP)
      IGRFIN=0
      DO 25 IGRCD=1,NGCOND
      COURIN(IGRCD)=0.0
      IGRDEB=IGRFIN+1
      IGRFIN=IGCOND(IGRCD)
      DO 20 IGR=IGRDEB,IGRFIN
      SUM=0.0D0
      SUD=0.0D0
      DO 10 IGAP=1,NGAP
      IREG=IFGAP(IGAP)
      SUM=SUM+FLUXES(IREG,IGR)*VOLUME(IREG)
      SUD=SUD+VOLUME(IREG)
   10 CONTINUE
      COURIN(IGRCD)=COURIN(IGRCD)+REAL(SUM/SUD)
   20 CONTINUE
      COURIN(IGRCD)=COURIN(IGRCD)/FLUHOM(IGRCD)
   25 CONTINUE
      DEALLOCATE(IFGAP)
      IF(IPRINT.GT.3) THEN
         WRITE(6,'(/19H EDIGAP: VALUES OF ,A,22H FLUXES OR ADF PER MAC,
     1   13HRO-GROUPS ARE)') TEXT8
         WRITE(6,'(1X,1P,10E13.5)') (COURIN(IGRCD),IGRCD=1,NGCOND)
      ENDIF
      RETURN
      END
