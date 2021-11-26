*DECK FMTDFD
      SUBROUTINE FMTDFD(NENTRY,KENTRY,IPRINT,IKFLU ,NTREG ,
     >                  NREG  ,NGROUP,NDIM  ,VOLUME,KEYFLX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To process the angular fluxes and generate the directional
* flux file.
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):
* G. Marleau
*
*Parameters: input
* NENTRY  number of data structures transfered to this module.
* KENTRY  data structure pointer.
* IPRINT  print level.
* IKFLU   pointer to the FLUX data structure.
* NTREG   number of regions for problem.
* NREG    number of unknowns for problem.
* NGROUP  number of groups for problem.
* NDIM    number of dimensions of problem.
* VOLUME  regional volumes.
* KEYFLX  index for regional fluxes in unknown vector.
*
*----------
*
      USE              GANLIB
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          NENTRY
      TYPE(C_PTR)      KENTRY(NENTRY)
      INTEGER          IPRINT,IKFLU
      INTEGER          NREG,NTREG,NGROUP,NDIM,KEYFLX(NTREG)
      REAL             VOLUME(NTREG)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='FMTDFD')
      INTEGER          ILCMUP,ILCMDN
      PARAMETER       (ILCMUP=1,ILCMDN=2)
*----
*  Local variables
*----
      INTEGER          ILONG,ITYLCM
      TYPE(C_PTR)      IPU,JPU
      INTEGER          IFPU,IGROUP,IR,NFLUX,IFTT
      CHARACTER*12     NAMFLX(2)
*----
*  Allocatable arrays
*----
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: AFLUX
*----
*  Work storage allocation
*----
      ALLOCATE(AFLUX(NREG,2,NGROUP))
*----
*  Initialize FLUX vectors
*----
      NFLUX=1
      NAMFLX(1)='FLUX        '
      IPU=KENTRY(IKFLU)
      CALL LCMLEN(IPU,'AFLUX',ILONG,ITYLCM)
      write(6,*) 'FLUXADJOINT ',ILONG,ITYLCM
      IF(ILONG .EQ. -1) THEN
        NFLUX=2
        NAMFLX(2)='ADJOINT     '
      ENDIF
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6000) NAMSBR,NGROUP,NDIM,NTREG,NREG
        WRITE(IOUT,6008)
        WRITE(IOUT,6012) (NAMFLX(IFTT),IFTT=1,NFLUX)
      ENDIF
*----
*  Get information from FLUX data structure.
*  1. Angular flux
*  2. Angular adjoint
*----
      JPU=LCMGID(IPU,'FLUX')
      DO IGROUP=1,NGROUP
        CALL LCMGDL(JPU,IGROUP,AFLUX(1,1,IGROUP))
      ENDDO
      IF(NFLUX .GT. 1) THEN
        JPU=LCMGID(IPU,'AFLUX')
        DO IGROUP=1,NGROUP
          CALL LCMGDL(JPU,IGROUP,AFLUX(1,1,IGROUP))
        ENDDO
      ENDIF
*----
*  Create output file
*----
      IFPU=FILUNIT(KENTRY(1))
      WRITE(IFPU,1000) NGROUP,NDIM,NREG,NFLUX
      WRITE(IFPU,1001) (NAMFLX(IFTT),IFTT=1,NFLUX)
*----
*  Print volumes
*----
      WRITE(IFPU,1002) (VOLUME(IR),IR=1,NREG)
*----
*  Print angular flux
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6003) 
      ENDIF
      DO IGROUP=1,NGROUP
        WRITE(IOUT,6002) IGROUP
        WRITE(IOUT,1002) (AFLUX(IR,1,IGROUP),IR=1,NREG)
        WRITE(IFPU,1002) (AFLUX(IR,1,IGROUP),IR=1,NREG)
      ENDDO
*----
*  Print scalar flux
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6004)
        DO IGROUP=1,NGROUP
          WRITE(IOUT,6002) IGROUP
          WRITE(IOUT,1002) (AFLUX(KEYFLX(IR),1,IGROUP),IR=1,NTREG)
          WRITE(IFPU,1002) (AFLUX(KEYFLX(IR),1,IGROUP),IR=1,NTREG)
        ENDDO
      ENDIF
*----
*  Print angular adjoint
*----
      IF(NFLUX .GT. 1) THEN
        IF(IPRINT .GE. 100) THEN
          WRITE(IOUT,6005)
        ENDIF
        DO IGROUP=1,NGROUP
          WRITE(IOUT,6002) IGROUP
          WRITE(IOUT,1002) (AFLUX(IR,2,IGROUP),IR=1,NREG)
          WRITE(IFPU,1002) (AFLUX(IR,2,IGROUP),IR=1,NREG)
        ENDDO
      ENDIF
*----
*  Work storage deallocation
*----
      DEALLOCATE(AFLUX)
*----
*  Processing finished, return
*----
      RETURN
*----
*  Formats
*----
 1000 FORMAT(5I10)
 1001 FORMAT(5(A12,2X))
 1002 FORMAT(1P,5E20.10)
 6000 FORMAT('Output from routine ',A6/
     >       'Number of groups  =',I5/
     >       'Number of dimens  =',I5/
     >       'Number of regions =',I5/
     >       'Number of unknowns=',I5)
 6002 FORMAT('Group = ',I5)
 6003 FORMAT('Direct angular flux per region ')
 6004 FORMAT('Scalar flux per region integrated from angular flux')
 6005 FORMAT('Adjoint angular flux per region')
 6008 FORMAT('Flux record types')
 6012 FORMAT(5(A12,2X))
      END
