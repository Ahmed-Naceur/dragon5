*DECK NXTIRA
      FUNCTION NXTIRA(XYCAR ,POSPIN,VOLINT)
*
*----------
*
*Purpose:
* Compute the volume of intersection between
* a rectangular region and an annular pin.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): G. Marleau.
*
*Parameters: input
* XYCAR   spatial description of the Cartesian region with:
*         XYCAR(1) for left face; XYCAR(2) for right face;
*         XYCAR(3) for bottom face; XYCAR(4) for top face
*         positions.
* POSPIN  spatial description of the annular pin region with
*         POSPIN(0) the radius; POSPIN(1) the $X$ position
*         of center; POSPIN(2) the $Y$ position
*         of center.
*
*Parameters: output
* NXTIRA  type of intersection between Cartesian region and
*         annular pin or annular region and Cartesian pin, where:
*         = 0 means that there is no intersection
*         between the two regions;
*         = 1 means that the Cartesian region
*         is all located inside the annular pin;
*         = 2 means that the annular pin
*         is all located inside the Cartesian region;
*         =-1 means that the intersection between
*         the annular pin and the Cartesian region is partial.
* VOLINT  2-D volume of intersection (area) between Cartesian region
*         and annular pin.
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*
*----
*
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          NXTIRA
      DOUBLE PRECISION XYCAR(4),POSPIN(0:2)
      DOUBLE PRECISION VOLINT
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTIRA')
      INTEGER          IPRINT
      PARAMETER       (IPRINT=100)
      DOUBLE PRECISION DCUTOF
      PARAMETER       (DCUTOF=1.0D-8)
      DOUBLE PRECISION DZERO,DONE
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0)
*----
*  Functions
*----
      DOUBLE PRECISION XDRCST,PI
*----
*  Local variables
*----
      INTEGER          IFACE,ILOC(4),IDX,IDY
      DOUBLE PRECISION XYFACE(4),RP2,VOLCAR,VOLPIN,FACDIR,VSUB(4),
     >                 TRIANG,ALPHA,FACTX,FACTY,
     >                 DIST,QUARTA,CARTV
      DOUBLE PRECISION DT1,DT2,DT3
*----
*  Initialize NXTIRA and VOLINT and PI
*----
      IF(IPRINT .GE. 200) THEN
        WRITE(IOUT,6000) NAMSBR
        WRITE(IOUT,6010) (XYCAR(IFACE),IFACE=1,4)
        WRITE(IOUT,6011) (POSPIN(IFACE),IFACE=0,2)
      ENDIF
      NXTIRA=0
      VOLINT=DZERO
      PI=XDRCST('Pi',' ')
*----
*  Locate pin center at origin
*----
      XYFACE(1)=XYCAR(1)-POSPIN(1)
      XYFACE(2)=XYCAR(2)-POSPIN(1)
      XYFACE(3)=XYCAR(3)-POSPIN(2)
      XYFACE(4)=XYCAR(4)-POSPIN(2)
*----
*  Find location of each face with respect to annular region
*----
      DO 100 IFACE=1,4
        IF(XYFACE(IFACE) .LE. -POSPIN(0)) THEN
*----
*  Plane to the left or under
*----
          ILOC(IFACE)=(-1)**IFACE
        ELSE IF(XYFACE(IFACE) .GE. POSPIN(0)) THEN
*----
*  Plane to the right or above
*----
          ILOC(IFACE)=(-1)**(IFACE+1)
        ELSE
*----
*  Plane croses annular region
*----
          ILOC(IFACE)=0
        ENDIF
 100  CONTINUE
      IF(ILOC(1) .NE. 1 .AND. ILOC(2) .NE. 1 .AND.
     >   ILOC(3) .NE. 1 .AND. ILOC(4) .NE. 1 ) THEN
        RP2=POSPIN(0)*POSPIN(0)
        VOLPIN=PI*RP2
        VOLCAR=(XYFACE(2)-XYFACE(1))*(XYFACE(4)-XYFACE(3))
        VOLINT=VOLPIN
*----
*  Find annular surface
*  1- to the left of X-
*  2- to the right of X+
*  3- below Y-
*  4- above Y+
*----
        FACDIR=-DONE
        DO 110 IFACE=1,4
          IF(ILOC(IFACE) .EQ. -1) THEN
            VSUB(IFACE)=DZERO
          ELSE
            TRIANG=SQRT(RP2-XYFACE(IFACE)*XYFACE(IFACE))
     >            *FACDIR*XYFACE(IFACE)
            ALPHA=ACOS(FACDIR*XYFACE(IFACE)/POSPIN(0))
            VSUB(IFACE)=RP2*ALPHA-TRIANG
            VOLINT=VOLINT-VSUB(IFACE)
          ENDIF
          FACDIR=-FACDIR
 110    CONTINUE
*----
*  For the case where two faces intersect inside annular region
*  compute intersections between the two surfaces VSUB
*  associated with each of these faces.
*----
        FACTX=DONE
        DO 120 IDX=1,2
          FACTY=DONE
          DO 130 IDY=3,4
            IF(ILOC(IDX) .EQ. 0 .AND. ILOC(IDY) .EQ. 0) THEN
              DIST=XYFACE(IDX)*XYFACE(IDX)+XYFACE(IDY)*XYFACE(IDY)
              IF(DIST .LT. RP2) THEN
                QUARTA=0.25D0*PI*RP2
                CARTV=FACTX*XYFACE(IDY)*FACTY*XYFACE(IDX)
                CARTV=CARTV+0.5D0*(VSUB(IDX)+VSUB(IDY))-QUARTA
              ELSE
                IF(FACTX*XYFACE(IDX) .LT. DZERO) THEN
                  IF(FACTY*XYFACE(IDY) .LT. DZERO) THEN
                    CARTV=0.0
                  ELSE
                    CARTV=VSUB(IDX)
                  ENDIF
                ELSE
                  IF(FACTY*XYFACE(IDY) .LT. DZERO) THEN
                    CARTV=VSUB(IDY)
                  ELSE
                    CARTV=-(VOLPIN-VSUB(IDX)-VSUB(IDY))
                  ENDIF
                ENDIF
              ENDIF
              VOLINT=VOLINT+CARTV
            ENDIF
            FACTY=-FACTY
 130      CONTINUE
          FACTX=-FACTX
 120    CONTINUE
        DT1=ABS(VOLINT-VOLPIN)
        DT2=ABS(VOLINT-VOLCAR)
        DT3=ABS(VOLINT)
        IF(DT1 .LT. DCUTOF) THEN
          VOLINT=VOLPIN
          NXTIRA=2
        ELSE IF(DT2 .LT. DCUTOF) THEN
          VOLINT=VOLCAR
          NXTIRA=1
        ELSE IF(DT3 .LT. DCUTOF) THEN
          VOLINT=DZERO
          NXTIRA=0
        ELSE
          NXTIRA=-1
        ENDIF
      ENDIF
*----
      IF(IPRINT .GE. 200) THEN
        WRITE(IOUT,6012) NAMSBR,NXTIRA,VOLINT
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6010 FORMAT('XYCAR ={',3(F20.10,','),F20.10,'};')
 6011 FORMAT('POSPIN={',2(F20.10,','),F20.10,'};')
 6012 FORMAT(A6,'={',I5,',',F20.10,'};')
      END
