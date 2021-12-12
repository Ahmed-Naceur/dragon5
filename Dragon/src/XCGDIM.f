*DECK XCGDIM
      SUBROUTINE XCGDIM(IPGEOM,MREGIO,NSOUT,IROT,IAPP,MAXJ,NVOL,
     >                  NBAN,MNAN,NRT,MSROD,MAROD,NSURF)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Initialize dimension for 2-D cluster geometry.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau
*
*Parameters: input
* IPGEOM  pointer to the geometry.
* MREGIO  maximum number of regions.
* NSOUT   number of surface for outer region.
* IROT    type of PIJ reconstruction:
*         =0  CP calculations;
*         =1 direct JPM reconstruction;
*         =2 rot2 type reconstruction.
* IAPP    type of surface conditions
*         level of dp approximation for jpm:
*         =1 DP0 all;
*         =2 DP1 all (default);
*         =3 DP1 inside DP0 outside symmetry conditions for CP.
* MAXJ    maximum number of currents.
*         Unused for CP calculations.
*
*Parameters: output
* NVOL    number of regions.
* NBAN    number of concentric regions.
* MNAN    maximum nunber of radius to read.
* NRT     number of rod types.
* MSROD   maximum number of subrods per rods.
* MAROD   maximum number of rods an annulus.
* NSURF   maximum number real surfaces.
*         Unused for CP calculation.
*
*----------------------------------------------------------------------
*
      USE         GANLIB
      PARAMETER  (NSTATE=40)
      TYPE(C_PTR) IPGEOM
      INTEGER     MREGIO,NSOUT,IROT,IAPP,MAXJ,NVOL,
     >            NBAN,MNAN,NRT,MSROD,MAROD,NSURF,ISTATE(NSTATE)
      CHARACTER   CMSG*131,TEXT12*12
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: JSPLIT,JGEOM
*----
*  CHECK FOR VALID IROT AND GEOMETRY
*----
      IF(IROT.GT.2.OR.IROT.LT.0)
     >  CALL XABORT('XCGDIM: UNABLE TO PROCESS THE GEOMETRY.')
      CALL XDRSET(ISTATE,NSTATE,0)
      CALL LCMGET(IPGEOM,'STATE-VECTOR',ISTATE)
*----
*  CHECK FOR INVALID GEOMETRY OPTIONS
*  ISTATE( 8) -> CELL IS INVALID
*  ISTATE(10) -> MERGE IS INVALID
*  ISTATE(11) -> SPLIT IS INVALID FOR CLUSTER ANNULUS
*----
      IF ( (ISTATE(8).NE.0).OR.(ISTATE(10).NE.0) )
     >     CALL XABORT('XCGDIM: UNABLE TO PROCESS THE GEOMETRY.')
      IF(ISTATE(11).EQ.0) THEN
        NVOL=ISTATE(6)
      ELSE
        CALL LCMLEN(IPGEOM,'SPLITR',NSPLIT,ITYPE)
        IF(ITYPE.NE.1)
     >    CALL XABORT('XCGDIM: SPLIT RECORD ON LCM IS NOT INTEGER')
        ALLOCATE(JSPLIT(NSPLIT))
        CALL LCMGET(IPGEOM,'SPLITR',JSPLIT)
        IF(NSOUT.GT.1) THEN
          NVOL=1
        ELSE
          NVOL=0
        ENDIF
        DO 135 ISPLIT=1,NSPLIT
          NVOL=NVOL+ABS(JSPLIT(ISPLIT))
 135    CONTINUE
        DEALLOCATE(JSPLIT)
      ENDIF
      NBAN=NVOL
      MNAN=NBAN+1
      IF(NSOUT.EQ.4) THEN
        MNAN=MNAN+3
      ENDIF
      IF(NSOUT.EQ.4) THEN
        NSURF=2*NVOL+2
      ELSE IF(NSOUT.EQ.6) THEN
        NSURF=2*NVOL+4
      ELSE
        NSURF=2*NVOL-1
      ENDIF
*----
*  COUNT NUMBER OF ROD TYPES IN CLUSTER
*----
      CALL LCMLEN(IPGEOM,'CLUSTER',ILONG,ITYPE)
      IF(ITYPE.NE.3)
     >   CALL XABORT('XCGDIM: CLUSTER RECORD ON LCM IS NOT CHARACTER')
      NRT=ILONG/3
      IF(ISTATE(9).LT.NRT) THEN
        WRITE(CMSG,9001) ISTATE(9),NRT
        CALL XABORT(CMSG)
      ENDIF
      ALLOCATE(JGEOM(ILONG))
      IPOS=1
      MSROD=1
      MAROD=1
      CALL LCMGET(IPGEOM,'CLUSTER',JGEOM)
*----
*  FOR EACH ROD TYPE FIND NUMBER OF SUBRODS AND NUMBER OF PINS
*----
      DO 120 IRT=1,NRT
        WRITE(TEXT12(1:4),'(A4)')  JGEOM(IPOS)
        WRITE(TEXT12(5:8),'(A4)')  JGEOM(IPOS+1)
        WRITE(TEXT12(9:12),'(A4)') JGEOM(IPOS+2)
        IPOS=IPOS+3
        CALL LCMSIX(IPGEOM,TEXT12,1)
        CALL XDRSET(ISTATE,NSTATE,0)
        CALL LCMGET(IPGEOM,'STATE-VECTOR',ISTATE)
        CALL LCMGET(IPGEOM,'RPIN',RPIN)
        CALL LCMGET(IPGEOM,'NPIN',NPIN)
        MAROD=MAX(MAROD,NPIN)
        IF(RPIN.EQ.0.0) THEN
          IROTS=1
        ELSE
          IROTS=IROT
        ENDIF
        IF(ISTATE(1).NE.3) THEN
          WRITE(CMSG,9002) ISTATE(1)
          CALL XABORT(CMSG)
        ENDIF
        IF(ISTATE(11).EQ.0) THEN
          NVOL=NVOL+ISTATE(6)
          NMSROD=ISTATE(6)
          IF(IROT.GT.0) NSURF=NSURF+2*IROTS*ISTATE(6)
        ELSE
          CALL LCMLEN(IPGEOM,'SPLITR',NSPLIT,ITYPE)
          IF(ITYPE.NE.1)
     >      CALL XABORT('XCGDIM: SPLIT RECORD ON LCM IS NOT INTEGER')
          ALLOCATE(JSPLIT(NSPLIT))
          CALL LCMGET(IPGEOM,'SPLITR',JSPLIT)
          NMSROD=0
          DO 130 ISPLIT=1,NSPLIT
            NMSROD=NMSROD+ABS(JSPLIT(ISPLIT))
            NVOL=NVOL+ABS(JSPLIT(ISPLIT))
 130      CONTINUE
          IF(IROT.GT.0) NSURF=NSURF+2*IROTS*NMSROD
          DEALLOCATE(JSPLIT)
        ENDIF
        MSROD=MAX(MSROD,NMSROD)
        CALL LCMSIX(IPGEOM,' ',2)
 120  CONTINUE
      MNAN=MAX(MNAN,MSROD+1)
      DEALLOCATE(JGEOM)
*----
*  CHECK IF NUMBER OF REGIONS IS ADEQUATE
*----
      IF (NVOL.GT.MREGIO) THEN
        WRITE(CMSG,9003) MREGIO,NVOL
        CALL XABORT(CMSG)
      ENDIF
      IF(IROT.GT.0) THEN
        IF(IAPP.EQ.3) THEN
          IAPPR=2
        ELSE
          IAPPR=IAPP
        ENDIF
        IF(NSOUT.EQ.4) THEN
          NSURF=NSURF*IAPPR+4
        ELSE IF(NSOUT.EQ.6) THEN
          NSURF=NSURF*IAPPR+6
        ELSE
        NSURF=NSURF*IAPPR
        ENDIF
        IF(MAXJ.LT.NSURF) THEN
          WRITE(CMSG,9004) NSURF,MAXJ
          CALL XABORT(CMSG)
        ENDIF
      ELSE
        NSURF=1
        IF(NSOUT.EQ.6) THEN
          CALL LCMGET(IPGEOM,'IHEX',IHEX)
          IF(IHEX.EQ.1) THEN
            IAPP=12
          ELSE IF(IHEX.EQ.3) THEN
            IAPP=6
          ENDIF
        ENDIF
      ENDIF
      RETURN
*----
*  ERROR MESSAGES FORMATS
*----
 9001 FORMAT('XCGDIM: ONLY ',I10,5X,'SUB GEOMETRIES ON LCM WHILE ',5X,
     >       I10,5X,'SUB GEOMETRIES ARE REQUIRED BY CLUSTER')
 9002 FORMAT('XCGDIM: ',I10,5X,'IS AN ILLEGAL GEOMETRY INSIDE CLUSTER')
 9003 FORMAT('XCGDIM: MAXIMUM NUMBER OF REGION ALLOCATED =',I10,
     >       5X,'NUMBER OF REGION REQUIRED =',I10)
 9004 FORMAT('XCGDIM: NUMBER OF CURRENT=',I10,5X,'IS LARGER THAN ',
     >'ALLOWED MAXIMUM VALUE MAXJ=',I10)
      END
