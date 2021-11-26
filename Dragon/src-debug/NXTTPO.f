*DECK NXTTPO
      SUBROUTINE NXTTPO(IPGEO ,IPRINT,ITYPBC,NBGCLS,NTPIN ,MAXMSH,
     >                  NCDIM ,IDIRR ,DRW   ,OFFCEN,NAGCLS)
*
*----------
*
*Purpose:
* To test that cluster pins do not overlapp.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s):
* G. Marleau.
*
*Parameters: input
* IPGEO   pointer to the GEOMETRY data structure.
* IPRINT  print level.
* ITYPBC  type of boundary conditions where:
*         =0 for geometry with Cartesianb oundaries;
*         =1 for geometry with annular boundary;
*         =2 for geometry with hexagonal boundary.
* NBGCLS  number of cluster sets.
* NTPIN   total number of pins.
* MAXMSH  maximum mesh dimension.
* NCDIM   number of dimensions.
* IDIRR   mesh direction.
* DRW     cell dimensions.
* OFFCEN  off centering of pin and annular regions in cell.
* NAGCLS  the cluster names in an integer format.
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*  \\\\
*  This routine is based on the XELDCL routine written by
*  R. Roy for the EXCELT: module. It contains an additional
*  level for cluster subgeometry analysis.
*
*----------
*
      USE              GANLIB
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      TYPE(C_PTR)      IPGEO
      INTEGER          IPRINT,ITYPBC,NBGCLS,NTPIN ,MAXMSH,IDIRR,NCDIM
      DOUBLE PRECISION DRW(3)
      REAL             OFFCEN(3)
      INTEGER          NAGCLS(3,NBGCLS)
*----
*  Functions
*----
      DOUBLE PRECISION XDRCST,PI
      INTEGER          NXTIRA,NXTIAA,NXTIHA,INTTYP
      DOUBLE PRECISION VOLINT
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTTPO')
      DOUBLE PRECISION DZERO,DONE,DTWO
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0,DTWO=2.0D0)
*----
*  Local variables
*----
      INTEGER          ICLS,NPIN,IX,IY,IZ,IP,JP,ITPIN,ICN,ITC
      INTEGER          ILCMLN,ILCMTY
      CHARACTER        NAMCL*12,NAMREC*12
      REAL             DELTA
      DOUBLE PRECISION ZPIN(2),XYZCAR(6),POSAH(0:2)
*----
*  Allocatable arrays
*----
      REAL, ALLOCATABLE, DIMENSION(:,:) :: PINDIM
      REAL, ALLOCATABLE, DIMENSION(:) :: RPIN,APIN
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: PINPOS
*----
*  Data
*----
      CHARACTER        CDIR(4)*1
      SAVE             CDIR
      DATA             CDIR /'X','Y','Z','R'/
*----
*  Scratch storage allocation
*   PINDIM  temporary storage for pin radius.
*   RPIN    pin radius.
*   APIN    pin angles.
*   PINPOS  pin position and outer radius.
*----
      ALLOCATE(PINDIM(0:MAXMSH,2),RPIN(NTPIN),APIN(NTPIN))
      ALLOCATE(PINPOS(0:4,NTPIN))
*----
*  Processing starts:
*  print routine openning output header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. -500) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      PI=XDRCST('Pi',' ')
*----
*  Define cell limits
*----
      IX=MOD(IDIRR-1,3)+1
      IY=MOD(IDIRR,3)+1
      IZ=MOD(IDIRR+1,3)+1
      IF(ITYPBC .EQ. 0) THEN
        XYZCAR(1)=-DRW(IX)/DTWO-DBLE(OFFCEN(IX))
        XYZCAR(2)=DRW(IX)/DTWO-DBLE(OFFCEN(IX))
        XYZCAR(3)=-DRW(IY)/DTWO-DBLE(OFFCEN(IY))
        XYZCAR(4)=DRW(IY)/DTWO-DBLE(OFFCEN(IY))
      ELSE IF(ITYPBC .EQ. 1) THEN
        POSAH(0)=DRW(IX)
        POSAH(1)=-DBLE(OFFCEN(IX))
        POSAH(2)=-DBLE(OFFCEN(IY))
      ELSE IF(ITYPBC .EQ. 2) THEN
        POSAH(0)=DRW(IX)
        POSAH(1)=-DBLE(OFFCEN(IX))
        POSAH(2)=-DBLE(OFFCEN(IY))
      ENDIF
      IF(NCDIM .EQ. 3) THEN
        XYZCAR(5)=-DRW(IZ)/DTWO-DBLE(OFFCEN(IZ))
        XYZCAR(6)=DRW(IZ)/DTWO-DBLE(OFFCEN(IZ))
        NAMREC='MESH'//CDIR(IZ)//'       '
      ENDIF
*----
*  Find pin locations
*----
      ITPIN=0
      DO ICLS=1,NBGCLS
        ICN=3*(ICLS-1)
        WRITE(NAMCL,'(3A4)') (NAGCLS(ITC,ICLS),ITC=1,3)
        CALL LCMSIX(IPGEO,NAMCL,1)
        CALL LCMGET(IPGEO,'NPIN',NPIN)
        CALL LCMLEN(IPGEO,'RPIN',ILCMLN,ILCMTY)
        IF(ILCMLN .EQ. 0) THEN
           CALL XABORT(NAMSBR//
     >    ': RPIN is absent for pin cluster')
        ELSE IF(ILCMLN .EQ. 1) THEN
          CALL LCMGET(IPGEO,'RPIN',RPIN)
          DO IP=2,NPIN
            RPIN(IP)=RPIN(1)
          ENDDO
        ELSE IF(ILCMLN .EQ. NPIN) THEN
          CALL LCMGET(IPGEO,'RPIN',RPIN)
        ELSE
          CALL XABORT(NAMSBR//
     >    ': Length of RPIN vector is invalid')
        ENDIF
        CALL LCMLEN(IPGEO,'APIN',ILCMLN,ILCMTY)
        IF(ILCMLN .EQ. 0) THEN
          APIN(1)=0.0
          DELTA=REAL((DTWO*PI)/DBLE(NPIN))
          DO IP=2,NPIN
            APIN(IP)=APIN(IP-1)+DELTA
          ENDDO
        ELSE IF(ILCMLN .EQ. 1) THEN
          CALL LCMGET(IPGEO,'APIN',APIN)
          DELTA=REAL((DTWO*PI)/DBLE(NPIN))
          DO IP=2,NPIN
            APIN(IP)=APIN(IP-1)+DELTA
          ENDDO
        ELSE IF(ILCMLN .EQ. NPIN) THEN
          CALL LCMGET(IPGEO,'APIN',APIN)
        ELSE
          CALL XABORT(NAMSBR//
     >    ': Length of APIN vector is invalid')
        ENDIF
        CALL LCMLEN(IPGEO,'RADIUS',ILCMLN,ILCMTY)
        IF(ILCMLN .EQ. 0) CALL XABORT(NAMSBR//
     >    ': RADIUS vector for pin is missing')
        CALL LCMGET(IPGEO,'RADIUS',PINDIM(0,1))
        IF(NCDIM .EQ. 3) THEN
          CALL LCMLEN(IPGEO,NAMREC,ILCMLN,ILCMTY)
          IF(ILCMLN .EQ. 0) CALL XABORT(NAMSBR//
     >      ': '//NAMREC//' vector for pin is missing')
          CALL LCMGET(IPGEO,NAMREC,PINDIM(0,2))
          ZPIN(2)=DBLE(PINDIM(ILCMLN-1,2)-PINDIM(0,2))/DTWO
          ZPIN(1)=-ZPIN(2)
          IF(ZPIN(1) .LT. XYZCAR(5) .OR.
     >       ZPIN(1) .GT. XYZCAR(6) .OR.
     >       ZPIN(2) .LT. XYZCAR(5) .OR.
     >       ZPIN(2) .GT. XYZCAR(6) ) CALL XABORT(NAMSBR//
     >      ': '//NAMREC//' pin extend outside cell')
        ELSE
          ZPIN(2)=DZERO
          ZPIN(1)=-ZPIN(2)
        ENDIF
*----
*  Store information in PINPOS
*----
        DO IP=1,NPIN
          ITPIN=ITPIN+1
          PINPOS(0,ITPIN)=DBLE(PINDIM(ILCMLN-1,1))
          PINPOS(1,ITPIN)=DBLE(RPIN(IP))*COS(APIN(IP))
          PINPOS(2,ITPIN)=DBLE(RPIN(IP))*SIN(APIN(IP))
          PINPOS(3,ITPIN)=ZPIN(1)
          PINPOS(4,ITPIN)=ZPIN(2)
        ENDDO
        CALL LCMSIX(IPGEO,NAMCL,2)
      ENDDO
*----
*  All pin localized, test for overlapp
*----
      NPIN=ITPIN
      IF(IPRINT .GE. 500) THEN
        IF(ITYPBC .EQ. 0) THEN
          WRITE(IOUT,6010) (XYZCAR(IX),IX=1,4)
        ELSE IF(ITYPBC .EQ. 1) THEN
          WRITE(IOUT,6011) (POSAH(IX),IX=0,2)
        ELSE IF(ITYPBC .EQ. 2) THEN
          WRITE(IOUT,6012) (POSAH(IX),IX=0,2)
        ENDIF
        IF(NCDIM .EQ. 3) THEN
          DO IP=1,NPIN
            WRITE(IOUT,6014) IP,(PINPOS(IX,IP),IX=0,4)
          ENDDO
        ELSE
          DO IP=1,NPIN
            WRITE(IOUT,6013) IP,(PINPOS(IX,IP),IX=0,2)
          ENDDO
        ENDIF
      ENDIF
      DO IP=1,NPIN
*----
*  Test if pin inside cell
*----
        INTTYP=-1
        IF(ITYPBC .EQ. 0) THEN
*----
*  Cell is a rectangle
*----
          INTTYP=NXTIRA(XYZCAR,PINPOS(0,IP),VOLINT)
          IF(INTTYP .NE. 2) CALL XABORT(NAMSBR//
     >    ': Rectangular cell does not contain completely the pin')
        ELSE IF(ITYPBC .EQ. 1) THEN
*----
*  Cell is a circle
*----
           INTTYP=NXTIAA(POSAH,PINPOS(0,IP),VOLINT)
          IF(INTTYP .NE. 2) CALL XABORT(NAMSBR//
     >    ': Annular cell does not contain completely the pin')
       ELSE IF(ITYPBC .EQ. 2) THEN
*----
*  Cell is an hexagon
*----
          INTTYP=NXTIHA(POSAH,PINPOS(0,IP),VOLINT)
          IF(INTTYP .NE. 2) CALL XABORT(NAMSBR//
     >    ': Hexagonal cell does not contain completely the pin')
        ENDIF
        IF(INTTYP .NE. 2) CALL XABORT(NAMSBR//
     >  ': Pin outside rectangular cell')
        DO JP=IP+1,NPIN
          IF(NCDIM .EQ. 3) THEN
*----
*  check for z-overlapp if required
*----
            IF(PINPOS(4,JP) .LT. PINPOS(3,IP) .OR.
     >         PINPOS(4,IP) .LT. PINPOS(3,JP) ) GO TO 100
          ENDIF
          INTTYP=NXTIAA(PINPOS(0,IP),PINPOS(0,JP),VOLINT)
          IF(INTTYP .NE. 0) CALL XABORT(NAMSBR//
     >  ': two pins overlapp')
 100      CONTINUE
        ENDDO
      ENDDO
*----
*  Processing finished:
*  print routine output header if required
*  and return
*----
      IF(IPRINT .GE. 500) THEN
        WRITE(IOUT,6001) NAMSBR
      ENDIF
*----
*  Scratch storage deallocation
*----
      DEALLOCATE(PINPOS)
      DEALLOCATE(APIN,RPIN,PINDIM)
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6010 FORMAT(' XYCAR= ',4F20.10)
 6011 FORMAT(' Annular cell   -> radius and x-y center= ',3F20.10)
 6012 FORMAT(' Hexagonal cell -> side and x-y center  = ',3F20.10)
 6013 FORMAT(' PIN ',I10,' -> radius and x-y center= ',3F20.10)
 6014 FORMAT(' PIN ',I10,' -> radius, x-y center and z location= ',
     >       5F20.10)
      END
