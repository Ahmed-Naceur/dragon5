*DECK NXTPHC
      SUBROUTINE NXTPHC(IPRINT,NDIM  ,MXMESH,MAXSUR,MAXREG,
     >                  MESH  ,DAMESH,NPIN  ,ITPIN ,DPIN  ,
     >                  NBSUR ,NBREG ,INDXSR,SURVOL,POSTRI)
*
*----------
*
*Purpose:
* Remove from the volumes or surfaces
* associated with an annular/hexagonal 2-D or 3-D geometry
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
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*
*Comments:
*  1- Contents of the DAMESH array:
*     hexagonal mesh is DAMESH(i,1) for i=0,MESH(1);
*     mesh in $Z$  is z(k)=DAMESH(k,3) for k=0,MESH(3);
*  2- Contents of the DPIN array for pin IPIN:
*     -> annular pin
*       ->annular regions in the $X-Y$ plane
*         centre (x,y,z)=(DPIN(0,IPIN)*COS(DPIN(-1,IPIN))
*                         DPIN(0,IPIN)*SIN(DPIN(-1,IPIN)),0.0D0)
*         outer pin radius      r=DPIN(4,IPIN)
*         pin height            dz(iz)=DPIN(3,IPIN)
*  3- Contents of the INDXSR array:
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
*       Note that for radial regions INDXSR(n,i)=-1 does not exists.
*----------
*
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          IPRINT,NDIM,MXMESH,MAXSUR,MAXREG
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
      PARAMETER       (IOUT=6,NAMSBR='NXTPHC')
      INTEGER          MAXDIM
      PARAMETER       (MAXDIM=4)
      DOUBLE PRECISION DCUTOF
      PARAMETER       (DCUTOF=1.0D-8)
      DOUBLE PRECISION DZERO,DONE,DTWO,DHALF,DSQ3O2
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0,DTWO=2.0D0,
     >                 DHALF=0.5D0,DSQ3O2=0.86602540378444D0)
*----
*  Functions
*----
      DOUBLE PRECISION XDRCST,PI
      INTEGER          NXTIAA,NXTITA,NXTPRA
      INTEGER          ITYITP,ITYITA,ITYTAP
      DOUBLE PRECISION VOLITP,VOLITA,VOLTAP
*----
*  Local variables
*----
      INTEGER          NX,NZ,NANN,NRTP,NRP,NSTP,NSP,NRP1,NRTPP,NRTPS,
     >                 ISBOT,ISTOP,IX,IZ,IR,ISECT,IPIN,ISUR,IVOL,ILOCT
      DOUBLE PRECISION ZB,ZT,PPRMIN,PPRMAX,PPPMIN,PPPMAX,DPP
      DOUBLE PRECISION VOLPIN,VOLIAO,VOLIAI
      DOUBLE PRECISION POSPIN(0:2),POSANN(0:2)
      INTEGER          NFACES
      INTEGER          IA
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ITYIAP
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: VOLIAP
*----
*  Get dimensioning information
*----
      ALLOCATE(ITYIAP(MXMESH),VOLIAP(MXMESH))
      NFACES=3
      PI=XDRCST('Pi',' ')
      NX=MESH(1)
      NZ=MESH(3)
      NANN=MESH(4)
      NRTP=NX**2
      NRP=6*NRTP
      NSTP=2*NX-1
      NSP=6*NSTP
      NRP1=NANN+1
      NRTPP=NRP*NRP1
      NRTPS=NRTP*NRP1
      ZB=0.0D0
      ZT=0.0D0
      IF(NDIM .EQ. 3) THEN
        ISBOT=-NSP*NZ
        ISTOP=ISBOT-NRTPP
      ELSE
        ISBOT=0
        ISTOP=0
      ENDIF
      POSANN(1)=0.0D0
      POSANN(2)=0.0D0
*----
*  Print mesh if required
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6000) NAMSBR
        WRITE(IOUT,6010) 'H',NX
        WRITE(IOUT,6011) 'MESHH ='
        WRITE(IOUT,6012) (DAMESH(IX,1),IX=-1,2*NX)
        IF(NZ .GT. 0) THEN
          WRITE(IOUT,6010) 'Z',NZ
          WRITE(IOUT,6011) 'MESHZ ='
          WRITE(IOUT,6012) (DAMESH(IZ,3),IZ=-1,NZ)
        ENDIF
          WRITE(IOUT,6010) 'R',NANN
        WRITE(IOUT,6011) 'MESHR ='
        WRITE(IOUT,6012) (DAMESH(IR,4),IR=-1,NANN)
*----
* Pin description
*----
        DO IPIN=1,NPIN
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
        ENDDO
      ENDIF
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
*  Determine pin-annular regions intersections
*----
        DO IA=1,NANN
          POSANN(0)=DAMESH(IA,4)
          ITYIAP(IA)=0
          VOLIAP(IA)=0.0D0
*----
*  Find 2-D annular region/annular pin intersection
*----
          ITYIAP(IA)=NXTIAA(POSANN,POSPIN,VOLIAP(IA))
        ENDDO
*        IF(NANN .GE. 0) THEN
*          WRITE(IOUT,*) 'Pin-Annular intersection for pin =',IPIN
*          WRITE(IOUT,'(2I10,F20.10)')
*     >     (IA,ITYIAP(IA),VOLIAP(IA),IA=1,NANN)
*        ENDIF
*----
*  1- Loop over sectors
*----
        DO ISECT=1,6
*----
*  Loop over region in sector
*----
          DO IR=1,NRTP
            ITYITP=NXTITA(POSTRI(1,1,IR,ISECT),POSPIN,VOLITP)
*            WRITE(IOUT,'(A25,3I10,F20.10)') 'Pin-Triangle intersection',
*     >      ISECT,IR,ITYITP,VOLITP
            ILOCT=(ISECT-1)*NRTPS+(IR-1)*NRP1
            IF(ITYITP .NE. 0) THEN
              VOLIAO=DZERO
              VOLIAI=DZERO
              DO IA=1,NANN
*                IF(IPRINT .GE. 100) THEN
*                  WRITE(IOUT,6022) IPIN,ISECT,IR,IA,ITYIAP(IA)
*                ENDIF
                POSANN(0)=DAMESH(IA,4)
                VOLIAI=VOLIAO
                IF(ITYIAP(IA) .NE. 0) THEN
                  ITYITA=NXTITA(POSTRI(1,1,IR,ISECT),POSANN,VOLITA)
*                  WRITE(IOUT,*) 'NXTITA ',ITYITA,VOLITA
                  IF(ITYITA .EQ. -1) THEN
*----
*  Partial Cartesian/annular region intersection
*  Examine Rectangle/pin intersection
*  Note: ITYITA=0 already considered above
*----
                    IF(ITYITP .EQ. -1) THEN
*----
*  Partial Cartesian/pin intersection
*  Examine Annular/pin intersection
*  Note: ITYIAP=0 already considered above
*----
                      IF(ITYIAP(IA) .EQ. -1) THEN
*----
*  Partial Annular/pin intersection
*  Find intersection volume of three regions
*----
                        ITYTAP=NXTPRA(NFACES,POSTRI(1,1,IR,ISECT),
     >                                POSANN,POSPIN,VOLTAP)
                        VOLIAO=VOLTAP
                        VOLIAI=VOLIAO-VOLIAI
                      ELSE IF(ITYIAP(IA) .EQ. 1) THEN
*----
*  Annular region in pin
*  Volume is given by Rectangle/annular intersection
*----
                        VOLIAO=VOLITA
                        VOLIAI=VOLIAO-VOLIAI
                      ELSE IF(ITYIAP(IA) .EQ. 2) THEN
*----
*  Annular region contains pin
*  Volume is given by Annular/pin intersection
*----
                        VOLIAO=VOLITP
                        VOLIAI=VOLIAO-VOLIAI
                      ENDIF
                    ELSE IF(ITYITP .EQ. 1) THEN
*----
*  Cartesian region in pin
*  Volume is given by Rectangle/annular intersection
*----
                      VOLIAO=VOLITA
                      VOLIAI=VOLIAO-VOLIAI
                    ELSE IF(ITYITP .EQ. 2) THEN
*----
*  Cartesian region contains pin
*  Volume is given by Annular/pin intersection
*----
                      VOLIAO=VOLIAP(IA)
                      VOLIAI=VOLIAO-VOLIAI
                    ENDIF
                  ELSE IF(ITYITA .EQ. 0) THEN
*----
*  No Cartesian/annular region intersection
*  go to next annular region.
*----
                    GO TO 125
                  ELSE IF(ITYITA .EQ. 1) THEN
*----
*  Cartesian region in annular region
*  Volume is given by Rectangle/pin intersection
*----
                    VOLIAO=VOLITP
                    VOLIAI=VOLIAO-VOLIAI
                  ELSE IF(ITYITA .EQ. 2) THEN
*----
*  Cartesian region contains annular region
*  Volume is given by Annular/pin intersection
*----
                    VOLIAO=VOLIAP(IA)
                    VOLIAI=VOLIAO-VOLIAI
                  ENDIF
                ENDIF
 125            CONTINUE
*----
*  There is an intersection possible between the pin and
*  the triangle.
*  1- Look for botton and top surface in 3-D
*----
*                IF(IPRINT .GE. 100) THEN
*                  WRITE(IOUT,6023) VOLIAI,VOLIAO
*                ENDIF
                IF(NDIM .EQ. 3) THEN
                  IF(ZB .LE. DAMESH(0,3) .AND.
     >               ZT .GE. DAMESH(0,3)) THEN
*----
*  Remove area contribution from bottom surface
*----
                    ISUR=ISBOT-ILOCT-IA
                    SURVOL(ISUR)=SURVOL(ISUR)-VOLIAI
                  ENDIF
                  IF(ZB .LE. DAMESH(MESH(3),3) .AND.
     >               ZT .GE. DAMESH(MESH(3),3)) THEN
*----
*  Remove area contribution from top surface
*----
                    ISUR=ISTOP-ILOCT-IA
                    SURVOL(ISUR)=SURVOL(ISUR)-VOLIAI
                  ENDIF
                  PPPMIN=ZB
                  PPPMAX=PPPMIN+DPIN(3,IPIN)
                  DO IZ=1,MESH(3)
                    PPRMIN=MAX(DAMESH(IZ-1,3),PPPMIN)
                    PPRMAX=MIN(DAMESH(IZ,3),PPPMAX)
                    IF(PPRMIN .LT. PPRMAX) THEN
                      DPP=VOLIAI*(PPRMAX-PPRMIN)
                      IVOL=(IZ-1)*NRTPP+ILOCT+IA
                      SURVOL(IVOL)=SURVOL(IVOL)-DPP
                    ENDIF
                  ENDDO
                ELSE
*                  WRITE(IOUT,'(A12,4I10,3F20.10)') 'Volume id',
*     >ISECT,IR,IA,ILOCT+IA,SURVOL(IVOL),VOLIAI,
*     >SURVOL(IVOL)-VOLIAI
                  IVOL=ILOCT+IA
                  SURVOL(IVOL)=SURVOL(IVOL)-VOLIAI
                ENDIF
*----
*  If pin all extracted, go to next pin
*----
                IF(VOLPIN .EQ. VOLIAO) GO TO 115
*                IF(VOLPIN .LE. DZERO) GO TO 115
              ENDDO
*----
*  Use DELV to correct volumes and surface area for
*  regions inside a rectangle but outside annular ring
*----
*              write(6,*) VOLIAI,VOLIAO,VOLITP,VOLITP-VOLIAI
              VOLIAI=VOLIAO
              VOLIAO=VOLITP
              VOLIAI=VOLIAO-VOLIAI
              IA=NRP1
              IF(NDIM .EQ. 3) THEN
                IF(ZB .LE. DAMESH(0,3) .AND.
     >             ZT .GE. DAMESH(0,3)) THEN
*----
*  Remove area contribution from bottom surface
*----
                  ISUR=ISBOT-ILOCT-IA
                  SURVOL(ISUR)=SURVOL(ISUR)-VOLIAI
                ENDIF
                IF(ZB .LE. DAMESH(NZ,3) .AND.
     >             ZT .GE. DAMESH(NZ,3)) THEN
*----
*  Remove area contribution from top surface
*----
                  ISUR=ISTOP-ILOCT-IA
                  SURVOL(ISUR)=SURVOL(ISUR)-VOLIAI
                ENDIF
                PPPMIN=ZB
                PPPMAX=PPPMIN+DPIN(3,IPIN)
                DO IZ=1,NZ
                  PPRMIN=MAX(DAMESH(IZ-1,3),PPPMIN)
                  PPRMAX=MIN(DAMESH(IZ,3),PPPMAX)
                  IF(PPRMIN .LT. PPRMAX) THEN
                    DPP=VOLIAI*(PPRMAX-PPRMIN)
                    IVOL=(IZ-1)*NRTPP+ILOCT+IA
                    SURVOL(IVOL)=SURVOL(IVOL)-DPP
                  ENDIF
                ENDDO
              ELSE
                IVOL=ILOCT+IA
                SURVOL(IVOL)=SURVOL(IVOL)-VOLIAI
              ENDIF
            ENDIF  
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
 6010 FORMAT(1X,'MESH DIMENSIONS IN ',A1,' =',I10)
 6011 FORMAT(1X,A7)
 6012 FORMAT(5F20.10)
 6017 FORMAT(A6,I4.4,'=',I10,';')
 6018 FORMAT(A6,I4.4,'=',F20.10,';')
 6019 FORMAT(A6,I4.4,'={',F20.10,',',F20.10,'};')
* 6022 FORMAT('Pin =',I10,' Sector = ',I10,' Region = ',I10,
*     >       ' Annulus =',I10,' Intersection type',I10)
* 6023 FORMAT(20X,' Volume of intersection= ',F20.10,
*     >       ' Volume remaining= ',F20.10)
*----
*  Error and Warning formats
*----
 9000 FORMAT('**** ERROR in -- ',A6,'-- found'/
     >       '     Area of region ',I5,' is negative')
 9001 FORMAT('**** ERROR in -- ',A6,'-- found'/
     >       '     Volume of region ',I5,' is negative')
 9002 FORMAT(5I10,F20.10)
      END
