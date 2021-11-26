!
!---------------------------------------------------------------------
!
!Purpose:
! To release allocated memory in SALT: module.
!
!Copyright:
! Copyright (C) 2014 Ecole Polytechnique de Montreal.
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version
!
!Author(s):
! A. Hebert
!
!---------------------------------------------------------------------
!
SUBROUTINE SALEND()
  USE SAL_GEOMETRY_MOD,   ONLY : GG
  USE SAL_TRACKING_TYPES, ONLY : ITRAC2,RTRAC2,IPART,RPART
  INTEGER :: OK,ELEM
  !----
  !  Release allocated memory for tracking data
  !----
  DEALLOCATE(RPART,RTRAC2,IPART,ITRAC2,STAT =OK)
  IF(OK /= 0) CALL XABORT('SALEND: failure to deallocate tracking storage')
  !----
  !  Release geometry allocated memory
  !----
  DEALLOCATE(GG%IPAR,GG%RPAR,GG%IBC2_ELEM,GG%ISURF2_ELEM,GG%VOL_NODE,GG%PPERIM_NODE, &
       GG%TYPE_BC2,GG%IDATA_BC2,GG%PERIM_MAC2,GG%PPERIM_MAC2,GG%MED,  &
       GG%BCDATA, GG%PERIM_NODE,STAT =OK)
  IF(OK /= 0) CALL XABORT('SALEND: failure to deallocate GG members')
  
  IF(GG%NB_SURF2>0)  THEN
     DEALLOCATE(GG%IBC2_SURF2,GG%IELEM_SURF2,GG%SURF2,STAT =OK)
     IF(OK /= 0) CALL XABORT('SALEND: failure to deallocate GG surf members')
  ENDIF
  DEALLOCATE(GG%NUM_MERGE,STAT =OK)
  IF(OK /= 0) CALL XABORT('SALEND: failure to deallocate GG%NUM_MERGE')
  IF(GG%NBBCDA>0)  THEN
    DO ELEM=1,GG%NBBCDA
      DEALLOCATE(GG%BCDATAREAD(ELEM)%ELEMNB,STAT =OK)
      IF(OK /= 0) CALL XABORT('SALEND: failure to deallocate GG%BCDATAREAD(ELEM)%ELEMNB')
    ENDDO
    DEALLOCATE(GG%BCDATAREAD,STAT =OK)
    IF(OK /= 0) CALL XABORT('SALEND: FAILURE TO DEALLOCATE GG%BCDATAREAD')
  ENDIF
  DEALLOCATE(GG, STAT= OK)
  IF(OK /= 0) CALL XABORT('SALEND: failure to deallocate GG')
END SUBROUTINE SALEND
