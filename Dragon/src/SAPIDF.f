*DECK SAPIDF
      SUBROUTINE SAPIDF(IPSAP,IPEDIT,NG,NMIL,ICAL,FNORM,REGFLX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To store discontinuity factor information in the Saphyb.
*
*Copyright:
* Copyright (C) 2015 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPSAP   pointer to the Saphyb.
* IPEDIT  pointer to the edition object (L_EDIT signature).
* NG      number of condensed energy groups.
* NMIL    number of mixtures.
* ICAL    index of the current elementary calculation.
* FNORM   flux normalization factor.
* REGFLX  averaged flux in the complete geometry.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPSAP,IPEDIT
      INTEGER NG,NMIL,ICAL
      REAL FNORM,REGFLX(NG)
*----
*  LOCAL VARIABLES
*----
      CHARACTER DIRNAM*12,HSMG*131
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: SURF
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SURFLX
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: HADF
*----
*  RECOVER DISCONTINUITY FACTOR INFORMATION FROM MACROLIB
*----
      IF(NMIL.NE.1) CALL XABORT('SAPIDF: NMIL=1 EXPECTED.')
      CALL LCMSIX(IPEDIT,'MACROLIB',1)
      CALL LCMLEN(IPEDIT,'ADF',ILONG,ITYLCM)
      IF(ILONG.EQ.0) CALL XABORT('SAPIDF: MISSING ADF DIRECTORY IN EDI'
     1 //'TION OBJECT.')
      CALL LCMSIX(IPEDIT,'ADF',1)
      CALL LCMGET(IPEDIT,'NTYPE',NSURFD)
      ALLOCATE(SURFLX(NSURFD,NG),SURF(NG),HADF(NSURFD))
      CALL LCMGTC(IPEDIT,'HADF',8,NSURFD,HADF)
      DO I=1,NSURFD
        CALL LCMLEN(IPEDIT,HADF(I),ILONG,ITYLCM)
        IF(ILONG.NE.NG) THEN
          WRITE(HSMG,'(12HSAPIDF: BAD ,A,8H LENGTH=,I5,10H EXPECTED=,
     1    I5,1H.)') HADF(I),ILONG,NG
          CALL XABORT(HSMG)
        ENDIF
        CALL LCMGET(IPEDIT,HADF(I),SURF)
        DO IGR=1,NG
          SURFLX(I,IGR)=SURF(IGR)*FNORM*1.0E13
        ENDDO
      ENDDO
      DEALLOCATE(HADF,SURF)
      CALL LCMSIX(IPEDIT,' ',2)
      CALL LCMSIX(IPEDIT,' ',2)
*----
*  MOVE TO THE 'calc' DIRECTORY.
*----
      WRITE(DIRNAM,'(''calc'',I8)') ICAL
      CALL LCMSIX(IPSAP,DIRNAM,1)
      CALL LCMSIX(IPSAP,'outflx',1)
      CALL LCMPUT(IPSAP,'REGFLX',NG,2,REGFLX)
      CALL LCMPUT(IPSAP,'SURFLX',NSURFD*NG,2,SURFLX)
      CALL LCMSIX(IPSAP,' ',2)
      CALL LCMSIX(IPSAP,' ',2)
      DEALLOCATE(SURFLX)
*----
*  CREATE dummy 'outgeom' DIRECTORY.
*----
      CALL LCMSIX(IPSAP,'geom',1)
      CALL LCMSIX(IPSAP,'outgeom',1)
      ALLOCATE(SURF(NSURFD))
      SURF(:)=1.0
      CALL LCMPUT(IPSAP,'SURF',NSURFD,2,SURF)
      DEALLOCATE(SURF)
      CALL LCMSIX(IPSAP,' ',2)
      CALL LCMSIX(IPSAP,' ',2)
      RETURN
      END
