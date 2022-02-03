*DECK NXTTLO
      SUBROUTINE NXTTLO(IPRINT,MXMESH,
     >                  MESH  ,DMESH ,POSTRI)
*
*----------
*
*Purpose:
* Locate triangles positions for HEXT, HEXTZ, HEXTCEL and HEXTCELZ.
*
*Copyright:
* Copyright (C) 2013 Ecole Polytechnique de Montreal
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
* MXMESH  maximum number of spatial subdivision in
*         $X$, $Y$ and $Z$.
* MESH    effective number of spatial subdivision in
*         each direction ($X$, $Y$ and $Z$).
* DMESH   spatial description of the Cartesian geometry.
* POSTRI  triangle position:
*         POSTRI(1,*,*,*) is X position;
*         POSTRI(2,*,*,*) is Y position;
*         POSTRI(*,1,*,*) is location of first corner;
*         POSTRI(*,2,*,*) is location of second corner;
*         POSTRI(*,3,*,*) is location of third corner;
*         POSTRI(*,*,i,j) is location of triangle i in cector j.
*
*----------
*
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          IPRINT,MXMESH
      INTEGER          MESH(4)
      DOUBLE PRECISION DMESH(-1:MXMESH,4)
      DOUBLE PRECISION POSTRI(2,3,MXMESH*MXMESH,6)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTTLO')
      DOUBLE PRECISION DZERO,DHALF,DSQ3O2
      PARAMETER       (DZERO=0.0D0,DHALF=0.5D0,
     >                 DSQ3O2=0.86602540378444D0)
*----
*  Local variables
*----
      INTEGER          NX,IREG,IX,IR,ISECT,ITRI
      DOUBLE PRECISION SIDET,SIDEL,SIDEH,
     >                 POSCX,POSCXD,POSCY,POSCYD
*----
*  Initialize POSTRI and get triangle side
*----
      IF(IPRINT .GE. 500) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      CALL XDDSET(POSTRI,36*MXMESH*MXMESH,DZERO)
      NX=MESH(1)
      IF(NX .EQ. 1) THEN
        SIDET=DMESH(1,1)-DMESH(0,1)
      ELSE
        SIDET=DMESH(2,1)-DMESH(1,1)
      ENDIF
      SIDEL=SIDET/DSQ3O2
      SIDEH=DHALF*SIDEL
*----
*  1- First sector
*----
      IREG=0
      POSCXD=DZERO
      POSCYD=SIDEH
      DO IX=1,NX
*----
*  Loop over right triangles on the line
*----
        POSCX=POSCXD
        POSCY=POSCYD
        DO IR=1,IX-1
          IREG=IREG+1
          POSTRI(1,1,IREG,1)=POSCX+SIDET
          POSTRI(2,1,IREG,1)=POSCY
          POSTRI(1,2,IREG,1)=POSCXD
          POSTRI(2,2,IREG,1)=POSCY+SIDEH
          POSTRI(1,3,IREG,1)=POSCXD
          POSTRI(2,3,IREG,1)=POSCY-SIDEH
          POSCY=POSCY+SIDEL
        ENDDO
*----
*  Loop over left triangles on the line
*----
        POSCYD=POSCYD-SIDEH
        POSCX=POSCXD
        POSCY=POSCYD
        DO IR=1,IX
          IREG=IREG+1
          POSTRI(1,1,IREG,1)=POSCX
          POSTRI(2,1,IREG,1)=POSCY
          POSTRI(1,2,IREG,1)=POSCX+SIDET
          POSTRI(2,2,IREG,1)=POSCY-SIDEH
          POSTRI(1,3,IREG,1)=POSCX+SIDET
          POSTRI(2,3,IREG,1)=POSCY+SIDEH
          POSCY=POSCY+SIDEL
        ENDDO
        POSCXD=POSCXD+SIDET
      ENDDO
*----
*  Five other sectors are just Pi/3 rotations of first sectors
*  triangles
*----
      DO ISECT=2,6
        DO IR=1,NX*NX
          DO ITRI=1,3
            POSTRI(1,ITRI,IR,ISECT)=DHALF*POSTRI(1,ITRI,IR,ISECT-1)
     >                             -DSQ3O2*POSTRI(2,ITRI,IR,ISECT-1)
            POSTRI(2,ITRI,IR,ISECT)=DSQ3O2*POSTRI(1,ITRI,IR,ISECT-1)
     >                             +DHALF*POSTRI(2,ITRI,IR,ISECT-1)
          ENDDO
        ENDDO
      ENDDO
      IF(IPRINT .GE. 500) THEN
        WRITE(IOUT,6020) 
        DO ISECT=1,6
          DO IR=1,NX*NX
            WRITE(IOUT,6021) ISECT,IR,
     >      (POSTRI(1,ITRI,IR,ISECT),POSTRI(2,ITRI,IR,ISECT),ITRI=1,3)
          ENDDO
        ENDDO
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6020 FORMAT('X-Y location of triangle corners')
 6021 FORMAT(2I5,6F20.10)
      END
