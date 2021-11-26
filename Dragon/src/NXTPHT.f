*DECK NXTPHT
      SUBROUTINE NXTPHT(IPRINT,NDIM  ,IDIRC ,MXMESH,MAXSUR,MAXREG,
     >                  MESH  ,DAMESH,NPIN  ,ITPIN ,DPIN  ,
     >                  NBSUR ,NBREG ,INDXSR,SURVOL,POSTRI)
*
*----------
*
*Purpose:
* Remove from the volumes or surfaces
* associated with a hexagonal 2-D or 3-D geometry
* the volumes or surfaces of the overlapping pins.
*
*Copyright:
* Copyright (C) 2010 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s):
* G. Marleau.
*
*Parameters: input
* IPRINT  print level.
* NDIM    dimension of problem.
* IDIRC   the direction of the first axis of a Cartesian geometry
*         assuming the axis are in a cyclic rotation.
* MXMESH  maximum number of spatial subdivision in
*         $X$, $Y$ and $Z$.
* MAXSUR  maximum number of surfaces in the geometry.
* MAXREG  maximum number of regions in the geometry.
* MESH    effective number of spatial subdivision in
*         each direction ($X$, $Y$ and $Z$).
* DAMESH  spatial description of the Cartesian geometry.
* NPIN    number of pins to superimpose on geometry.
* ITPIN   type of pin.
* DPIN    pin location and dimensions.
* NBSUR   number of surfaces in the geometry.
* NBREG   final number of non void regions in the geometry.
* POSTRI  triangle position:
*         POSTRI(1,*,*,*) is X position;
*         POSTRI(2,*,*,*) is Y position;
*         POSTRI(*,1,*,*) is location of first corner;
*         POSTRI(*,2,*,*) is location of second corner;
*         POSTRI(*,3,*,*) is location of third corner;
*         POSTRI(*,*,i,j) is location of triangle i in cector j.
*
*Parameters: input/output
* INDXSR  local indexing of regions.
* SURVOL  volume of regions.
*
*Comments:
*  1- Contents of IDIRC:
*     IDIRC      axes in 1-D   axes in 2-D   axes in 3-D
*     1          x             (x,y)         (x,y,z)
*     2          y             (y,z)         (y,z,x)
*     3          z             (z,x)         (z,x,y)
*  2- Contents of the DAMESH array:
*     hexagonal mesh is DAMESH(i,1) for i=0,MESH(1);
*     mesh in $Z$  is z(k)=DAMESH(k,3) for k=0,MESH(3);
*  3- Contents of the DPIN array for pin IPIN:
*     -> annular pin
*       ->annular regions in the $X-Y$ plane
*         centre (x,y,z)=(DPIN(0,IPIN)*COS(DPIN(-1,IPIN))
*                         DPIN(0,IPIN)*SIN(DPIN(-1,IPIN)),0.0D0)
*         outer pin radius      r=DPIN(4,IPIN)
*         pin height            dz(iz)=DPIN(3,IPIN)
*  4- Contents of the INDXSR array:
*     For i>0
*       INDXSR(1,i)= iu is the $U$ location of region i
*       INDXSR(2,i)= iv is the $V$ location of region i
*       INDXSR(3,i)= iz is the $Z$ location of region i
*       INDXSR(4,i)= ir is the $R$ location of region i
*       INDXSR(5,i)= iw is the $W$ location of region i
*     For i<0
*       INDXSR(1,i)= iu is the $U$ location of surface i
*       INDXSR(2,i)= iv is the $V$ location of surface i
*       INDXSR(3,i)= iz is the $Z$ location of surface i
*       INDXSR(4,i)= ir is the $R$ location of surface i
*       INDXSR(5,i)= iw is the $W$ location of surface i
*       with INDXSR(n,i)=-1 for surface associated with
*                           location 0 in direction n.
*       with INDXSR(n,i)=-2 for surface associated with
*                           location MESH(n) in direction n.
*       Note that for radial regions INDXSR(n,i)=-1 does not
*       exists.
*----------
*
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          IPRINT,NDIM,IDIRC,MXMESH,MAXSUR,MAXREG
      INTEGER          MESH(4)
      DOUBLE PRECISION DAMESH(-1:MXMESH,4)
      INTEGER          NPIN,ITPIN(3,NPIN)
      DOUBLE PRECISION DPIN(-1:4,NPIN)
      INTEGER          NBSUR,NBREG,INDXSR(5,-MAXSUR:MAXREG)
      DOUBLE PRECISION SURVOL(-MAXSUR:MAXREG)
      DOUBLE PRECISION POSTRI(2,3,MXMESH*MXMESH,6)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTPHT')
      INTEGER          MAXDIM
      PARAMETER       (MAXDIM=4)
      DOUBLE PRECISION DCUTOF
      PARAMETER       (DCUTOF=1.0D-8)
      DOUBLE PRECISION DZERO,DTWO
      PARAMETER       (DZERO=0.0D0,DTWO=2.0D0)
*----
*  Functions
*----
      DOUBLE PRECISION XDRCST,PI
      INTEGER          NXTITA
      INTEGER          ITYITA
      DOUBLE PRECISION VOLINT
*----
*  Local variables
*----
      INTEGER          NX,NZ,NRTP,NRP,NSTP,NSP,ISBOT,ISTOP,NRTPP,NRTPS,
     >                 NRP1,IZ,IR,ISECT,IPIN,ISUR,IVOL
      DOUBLE PRECISION ZB,ZT,VOLPIN,PPRMIN,PPRMAX,PPPMIN,PPPMAX,DPP
      DOUBLE PRECISION POSPIN(0:2)
*----
*  Prepare loops over spatial directions as a function
*  of IDIRC and NDIM.
*----
      PI=XDRCST('Pi',' ')
      NX=MESH(1)
      NZ=MESH(3)
      NRTP=NX**2
      NRP=6*NRTP
      NSTP=2*NX-1
      NSP=6*NSTP
      NRP1=1
      NRTPP=NRP*NRP1
      NRTPS=NRTP*NRP1
      ZB=0.0D0
      ZT=0.0D0
      IF(NZ .EQ. 0) THEN
        ISBOT=0
        ISTOP=0
      ELSE
        ISBOT=-NSP*NZ
        ISTOP=ISBOT-NRP
      ENDIF
*      write(6,*) 'MXMESH,MAXSUR,MAXREG',MXMESH,MAXSUR,MAXREG
*----
*  Print mesh if required
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6000) NAMSBR
        WRITE(IOUT,6010) NX,NZ
        WRITE(IOUT,6011) 'MESHH ='
        WRITE(IOUT,6012) (DAMESH(IR,1),IR=-1,2*NX)
        IF(NZ .GT. 0) THEN
          WRITE(IOUT,6011) 'MESHZ ='
          WRITE(IOUT,6012) (DAMESH(IZ,3),IZ=-1,NZ)
        ENDIF
*----
* Pin description
*----
        DO 610 IPIN=1,NPIN
          WRITE(IOUT,6017) 'PinTyp',IPIN,ITPIN(3,IPIN)
          IF(ITPIN(3,IPIN) .EQ. 3) THEN
            WRITE(IOUT,6019) 'PinCXY',IPIN,
     >                     DPIN(0,IPIN)*COS(DPIN(-1,IPIN)),
     >                     DPIN(0,IPIN)*SIN(DPIN(-1,IPIN))
            WRITE(IOUT,6018) 'PinRad',IPIN,DPIN(4,IPIN)
          ENDIF
          IF(NDIM .EQ. 3) THEN
            WRITE(IOUT,6019) 'PinPoZ',IPIN,
     >                    DAMESH(-1,3)-DPIN(3,IPIN)/DTWO,
     >                    DAMESH(-1,3)+DPIN(3,IPIN)/DTWO
          ENDIF
 610    CONTINUE
      ENDIF
      IF(IDIRC .NE. 1) CALL XABORT(NAMSBR//
     >': Only z-directed pins permitted in hexagons')
*----
*  Loop over pins
*----
      DO IPIN=1,NPIN
*----
*  For 3-D problem,
*  Find pin bottom (ZB) and top (ZT) z location.
*----
        IF(NDIM .EQ. 3) THEN
          ZB=DAMESH(-1,3)-DPIN(3,IPIN)/DTWO
          ZT=ZB+DPIN(3,IPIN)
        ENDIF
*----
*  Annular pin properties
*----
        POSPIN(0)=DPIN(4,IPIN)
        POSPIN(1)=DPIN(0,IPIN)*COS(DPIN(-1,IPIN))
        POSPIN(2)=DPIN(0,IPIN)*SIN(DPIN(-1,IPIN))
        VOLPIN=PI*POSPIN(0)*POSPIN(0)
*----
*  1- Loop over sectors
*----
        DO ISECT=1,6
*----
*  Loop over region in sector
*----
          DO IR=1,NRTP
*            write(6,*) 'ANALYSE IPIN,ISECT,IR', IPIN,ISECT,IR
            ITYITA=NXTITA(POSTRI(1,1,IR,ISECT),POSPIN,VOLINT)
            IF(ITYITA .NE. 0) THEN
              VOLPIN=VOLPIN-VOLINT
              IF(IPRINT .GE. 100) THEN
                WRITE(IOUT,6022) IPIN,ISECT,IR,ITYITA,VOLINT,VOLPIN
              ENDIF
*----
*  There is an intersection possible between the pin and
*  the triangle.
*  1- Look for botton and top surface in 3-D
*----
              IF(NDIM .EQ. 3) THEN
                IF(ZB .LE. DAMESH(0,3) .AND.
     >             ZT .GE. DAMESH(0,3)) THEN
*----
*  Remove area contribution from bottom surface
*----
                  ISUR=ISBOT-(ISECT-1)*NRTPS-IR
                  SURVOL(ISUR)=SURVOL(ISUR)-VOLINT
                ENDIF
                IF(ZB .LE. DAMESH(MESH(3),3) .AND.
     >             ZT .GE. DAMESH(MESH(3),3)) THEN
*----
*  Remove area contribution from top surface
*----
                  ISUR=ISTOP-(ISECT-1)*NRTPS-IR
                  SURVOL(ISUR)=SURVOL(ISUR)-VOLINT
                ENDIF
                PPPMIN=ZB
                PPPMAX=PPPMIN+DPIN(3,IPIN)
                DO IZ=1,MESH(3)
                  PPRMIN=MAX(DAMESH(IZ-1,3),PPPMIN)
                  PPRMAX=MIN(DAMESH(IZ,3),PPPMAX)
                  IF(PPRMIN .LT. PPRMAX) THEN
                    DPP=VOLINT*(PPRMAX-PPRMIN)
                    IVOL=(IZ-1)*NRTPP+(ISECT-1)*NRTPS+IR*NRP1
                    SURVOL(IVOL)=SURVOL(IVOL)-DPP
                  ENDIF
                ENDDO
              ELSE
                IVOL=(ISECT-1)*NRTPS+IR*NRP1
*                WRITE(IOUT,'(A12,3I10,3F20.10)') 'Volume id',
*     >ISECT,IR,(ISECT-1)*NRTPS+IR*NRP1,SURVOL(IVOL),VOLINT,
*     >SURVOL(IVOL)-VOLINT
                SURVOL(IVOL)=SURVOL(IVOL)-VOLINT
              ENDIF
            ENDIF
*----
*  If pin all extracted, go to next pin
*----
            IF(VOLPIN .LE. DZERO) GO TO 115
          ENDDO
        ENDDO
 115    CONTINUE
      ENDDO
*----
*  Test for negative surface area and volumes
*----
      DO ISUR=1,NBSUR
        IF(SURVOL(-ISUR) .LT. -DCUTOF) THEN
          WRITE(IOUT,9000) NAMSBR,-ISUR
          WRITE(IOUT,9002) (INDXSR(IR,-ISUR),IR=1,5),SURVOL(-ISUR)
          CALL XABORT(NAMSBR//
     >    ': Region with negative surface area detected')
        ELSE IF(SURVOL(-ISUR) .LT. DCUTOF) THEN
          SURVOL(-ISUR)=DZERO
        ENDIF
      ENDDO
      DO  IVOL=1,NBREG
        IF(SURVOL(IVOL) .LT. -DCUTOF) THEN
          WRITE(IOUT,9001) NAMSBR,IVOL
          WRITE(IOUT,9002) (INDXSR(IR,IVOL),IR=1,5),SURVOL(IVOL)
          CALL XABORT(NAMSBR//
     >    ': Region with negative volume detected')
        ELSE IF(SURVOL(IVOL) .LT. DCUTOF) THEN
          SURVOL(IVOL)=DZERO
        ENDIF
      ENDDO
*----
*  Print volumes if required
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6002) 'SurVol'
        WRITE(IOUT,6005) (IVOL,(INDXSR(IR,IVOL),IR=1,5),SURVOL(IVOL),
     >                     IVOL=-NBSUR,NBREG)
        WRITE(IOUT,6003)
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6002 FORMAT(A12,'={')
 6003 FORMAT('};')
 6005 FORMAT((6(I10,','),F20.10,:,','))
 6010 FORMAT(1X,'MESH DIMENSIONS =',2I10)
 6011 FORMAT(1X,A7)
 6012 FORMAT(5F20.10)
 6017 FORMAT(A6,I4.4,'=',I10,';')
 6018 FORMAT(A6,I4.4,'=',F20.10,';')
 6019 FORMAT(A6,I4.4,'={',F20.10,',',F20.10,'};')
 6022 FORMAT('Pin =',I10,' Sector = ',I10,' Region = ',I10,
     >       ' Intersection =',I10,
     >       ' Volume of intersection= ',F20.10,
     >       ' Volume remaining= ',F20.10)
*----
*  Error and Warning formats
*----
 9000 FORMAT('**** ERROR in -- ',A6,'-- found'/
     >       '     Area of region ',I5,' is negative')
 9001 FORMAT('**** ERROR in -- ',A6,'-- found'/
     >       '     Volume of region ',I5,' is negative')
 9002 FORMAT(5I10,F20.10)
      END
