*DECK NXTIHA
      FUNCTION NXTIHA(POSHEX ,PINPOS,VOLINT)
*
*----------
*
*Purpose:
* Compute the volume of intersection between
* a 2--D hexagon and an annular pin.
*
*Copyright:
* Copyright (C) 2010 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): G. Marleau.
*
*Parameters: input
* POSHEX  spatial description of the hexagon with:
*         POSHEX(0) the dimension of one of its sides;
*         POSHEX(1) the $X$ position of hexagon center;
*         POSHEX(2) the $Y$ position of hexagon center.
* PINPOS  spatial description of the annular pin region with:
*         PINPOS(0) the radius of the annular pin;
*         PINPOS(1) the $X$ position of the annular pin center;
*         PINPOS(2) the $Y$ position of the annular pin center.
*
*Parameters: output
* NXTIHA  type of intersection between haxagon and annular pin, where:
*         = 0 means that there is no intersection
*         between the two regions;
*         = 1 means that the hexagon
*         is all located inside the annular pin;
*         = 2 means that the annular pin
*         is all located inside the hexagon;
*         =-1 means that the intersection between
*         the hexagon and the annular pin is partial.
* VOLINT  2-D volume of intersection (area) between hexagon and
*         annular pin.
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
      INTEGER          NXTIHA
      DOUBLE PRECISION POSHEX(0:2),PINPOS(0:2)
      DOUBLE PRECISION VOLINT
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTIHA')
      INTEGER          IPRLOC
      PARAMETER       (IPRLOC=10)
      DOUBLE PRECISION DCUTOF
      PARAMETER       (DCUTOF=1.0D-8)
      DOUBLE PRECISION DZERO,DONE,DHALF,DSQ3O2
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0,
     >                 DHALF=0.5D0,DSQ3O2=0.86602540378444D0)
*----
*  Functions
*----
      INTEGER          NXTITA
      DOUBLE PRECISION XDRCST,PI
*----
*  Local variables
*----
      INTEGER          IDIR,ICORN,IFACE,NFPINS,NCIN,ITRI,INTTRI
      DOUBLE PRECISION POSTRI(2,3),RADC,DISTF,DISTC,RADP2,
     >                 VTPINT
*----
*  Data
*----
      DOUBLE PRECISION CORNRH(2,6),DIRFAC(2,6)
      SAVE             CORNRH,DIRFAC
      DATA   CORNRH   / 0.86602540378444D0,-0.5D0,
     >                  0.86602540378444D0, 0.5D0,
     >                  0.0D0             , 1.0D0,
     >                 -0.86602540378444D0, 0.5D0,
     >                 -0.86602540378444D0,-0.5D0,
     >                  0.0D0             ,-1.0D0/
      DATA   DIRFAC   /-1.0D0             , 0.0D0,
     >                 -0.5D0             ,-0.86602540378444D0,
     >                  0.5D0             ,-0.86602540378444D0,
     >                  1.0D0             , 0.0D0,
     >                  0.5D0             , 0.86602540378444D0,
     >                 -0.5D0             , 0.86602540378444D0/
*----
*  Print header if required
*----
      IF(IPRLOC .GE. 200) THEN
        WRITE(IOUT,6000) NAMSBR
        WRITE(IOUT,6010) (POSHEX(IFACE),IFACE=0,2)
        WRITE(IOUT,6011) (PINPOS(IFACE),IFACE=0,2)
      ENDIF
*----
*  Initialize PI, NXTIHA and VOLINT
*----
      PI=XDRCST('Pi',' ')
      NXTIHA=0
      VOLINT=DZERO
*----
*  Evaluate distance from FACES to pin center
*----
      RADP2=PINPOS(0)**2
      NFPINS=0
      NCIN=0
      ICORN=1
      POSTRI(1,ICORN)=POSHEX(1)
      POSTRI(2,ICORN)=POSHEX(2)
      DO IFACE=1,6
        DISTF=DZERO
        RADC=DZERO
*        write(6,*) 'CORNRH',(CORNRH(IDIR,IFACE),IDIR=1,2)
*        write(6,*) 'DIRFAC',(DIRFAC(IDIR,IFACE),IDIR=1,2)
        DO IDIR=1,2
          DISTC=PINPOS(IDIR)-(POSHEX(IDIR)+POSHEX(0)*CORNRH(IDIR,IFACE))
          DISTF=DISTF+DISTC*DIRFAC(IDIR,IFACE)
          RADC=RADC+DISTC**2
        ENDDO
*        write(6,*) DISTC,DISTF,RADC,PINPOS(0)
        IF(DISTF .LT. -PINPOS(0)) THEN
*----
*  Pin outside hexagon
*  Return
*----
          NXTIHA=0
          VOLINT=DZERO
          RETURN
        ELSE IF(DISTF .GE. PINPOS(0)) THEN
          NFPINS=NFPINS+1
        ENDIF
        IF(RADC .LT. RADP2) THEN
          NCIN=NCIN+1
        ENDIF
      ENDDO
      IF(NFPINS .EQ. 6) THEN
        NXTIHA=2
        VOLINT=PI*RADP2
      ELSE IF(NCIN .EQ. 6) THEN
        NXTIHA=1
        VOLINT=3.0D0*DSQ3O2*POSHEX(0)**2
      ELSE
        NXTIHA=-1
*----
*  First five triangles
*----
        DO ITRI=1,5
          DO ICORN=2,3
            DO IDIR=1,2
              POSTRI(IDIR,ICORN)=POSTRI(IDIR,ICORN)
     >                        +POSHEX(0)*CORNRH(IDIR,IFACE+ICORN-2)
            ENDDO
          ENDDO
          INTTRI=NXTITA(POSTRI,PINPOS,VTPINT)
          VOLINT=VOLINT+VTPINT
        ENDDO
*----
*  Last five triangles
*----
        DO IDIR=1,2
          POSTRI(IDIR,ICORN)=POSTRI(IDIR,ICORN)
     >                    +POSHEX(0)*CORNRH(IDIR,6)
          POSTRI(IDIR,ICORN)=POSTRI(IDIR,ICORN)
     >                    +POSHEX(0)*CORNRH(IDIR,1)
        ENDDO
        ITRI=NXTITA(POSTRI,PINPOS,VTPINT)
        VOLINT=VOLINT+VTPINT
      ENDIF
      IF(IPRLOC .GE. 200) THEN
        WRITE(IOUT,6012) NAMSBR,NXTIHA,VOLINT
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6010 FORMAT('POSHEX={',2(F20.10,','),F20.10,'};')
 6011 FORMAT('PINPOS={',2(F20.10,','),F20.10,'};')
 6012 FORMAT(A6,'={',I5,',',F20.10,'};')
      END
