*DECK EDISTA
      SUBROUTINE EDISTA(IPRINT,NMERGE,ITYPE,VOLMER,VOLREL,VOLTOT,
     >                  FLXNEW,FLXOLD,RATNEW,RATOLD)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Print homogenized/condensed macroscopic cross sections statistics.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau
*
*Parameters: input
* IPRINT  print level;
*         = 0 no print;
*         = 1 print fluxes;
*         = 2 1+print reaction rates;
*         = 3 2+print homogenized cross sections.
* NMERGE  number of regions.
* ITYPE   type of statistics:
*         = 1 flux relative errors;
*         = 2 reaction rates relative errors;
*         = 3 delta sigma.
* VOLMER  current region merged volumes.
* VOLREL  old volume/new volume.
* VOLTOT  total old volume.
* FLXNEW  new integrated flux.
* FLXOLD  old integrated flux.
* RATNEW  new homogenized cross sections.
* RATOLD  old homogenized cross sections.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER     IPRINT,NMERGE,ITYPE
      REAL        VOLMER(NMERGE),VOLREL,VOLTOT,FLXNEW(NMERGE),
     >            FLXOLD(NMERGE),RATNEW(NMERGE),RATOLD(NMERGE)
*----
*  LOCAL VARIABLES
*----
      PARAMETER  (IUNOUT=6)
      REAL, ALLOCATABLE, DIMENSION(:) :: VALERR
*----
*  SCRATCH STORAGE ALLOCATION
*   VALERR  relative error or delta sigma.
*----
      ALLOCATE(VALERR(NMERGE))
*
      IF(IPRINT.GT.2) THEN
        IF(ITYPE.NE.3) THEN
          WRITE(IUNOUT,6000)
        ELSE
          WRITE(IUNOUT,6001)
        ENDIF
      ENDIF
      EPSMAX=0.0
      EPSAVG=0.0
      RPSAVG=0.0
      DO 100 IREG=1,NMERGE
        IF(ITYPE.EQ.3) THEN
          CURVAL=RATNEW(IREG)
          OLDVAL=RATOLD(IREG)
          VALERR(IREG)=CURVAL-OLDVAL
        ELSE
          IF(ITYPE.EQ.1) THEN
            CURVAL=VOLREL*FLXNEW(IREG)
            OLDVAL=FLXOLD(IREG)
          ELSE IF(ITYPE.EQ.2) THEN
            CURVAL=VOLREL*RATNEW(IREG)*FLXNEW(IREG)
            OLDVAL=RATOLD(IREG)*FLXOLD(IREG)
          ENDIF
          IF(OLDVAL.NE.0.0) THEN
            VALERR(IREG)=100.0*(CURVAL-OLDVAL)/OLDVAL
          ELSE IF(CURVAL.NE.0.0) THEN
            VALERR(IREG)=100.0*(CURVAL-OLDVAL)/CURVAL
          ELSE
            VALERR(IREG)=0.0
          ENDIF
        ENDIF
        IF(IPRINT.GT.2) THEN
          WRITE(IUNOUT,6002) IREG,CURVAL,OLDVAL,VALERR(IREG)
        ENDIF
        IF(ITYPE.NE.3) THEN
          EPSMAX=MAX(EPSMAX,ABS(VALERR(IREG)))
          EPSAVG=EPSAVG+ABS(VALERR(IREG))*VOLMER(IREG)*VOLREL
          RPSAVG=RPSAVG+VALERR(IREG)*VALERR(IREG)
        ENDIF
 100  CONTINUE
      IF(ITYPE.NE.3) THEN
        IF(IPRINT.EQ.2) THEN
          WRITE(IUNOUT,6003)
          WRITE(IUNOUT,6006) (VALERR(IREG),IREG=1,NMERGE)
        ENDIF
        EPSAVG=EPSAVG/VOLTOT
        WRITE(IUNOUT,6005) EPSMAX,EPSAVG,SQRT(RPSAVG/NMERGE)
      ELSE IF(IPRINT.LE.2) THEN
        WRITE(IUNOUT,6004)
        WRITE(IUNOUT,6006) (VALERR(IREG),IREG=1,NMERGE)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(VALERR)
      RETURN
*----
*  FORMAT
*----
 6000 FORMAT(
     >4X,'REGION',13X,'CURRENT VALUE',10X,'REFERENCE',16X,' ERROR (%) ')
 6001 FORMAT(
     >4X,'REGION',13X,'CURRENT VALUE',10X,'REFERENCE',16X,'DELTA SIGMA')
 6002 FORMAT(4X,I5,10X,1P,E14.4,8X,E14.4,8X,E14.4)
 6003 FORMAT(' RELATIVE ERROR (NEW-OLD) ON FLUXES (%)')
 6004 FORMAT(' DELTA SIGMA (NEW-OLD)')
 6005 FORMAT(/4X,'                MAXIMUM ERROR=',F8.2,' %'/
     >        4X,'VOLUME WEIGHTED AVERAGE ERROR=',F8.2,' %'/
     >        4X,'                   RMS  ERROR=',F8.2,' %')
 6006 FORMAT(1P,7(3X,E15.7))
      END
