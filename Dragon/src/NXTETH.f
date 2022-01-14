*DECK NXTETH
      SUBROUTINE NXTETH(IPRINT,IDT   ,ILEV  ,NFSUR ,MAXMSP,
     >                  NS1   ,NR1   ,NS2   ,NR2   ,IEDIMG,DAMESH,
     >                  IX1   ,ID1   ,SV1   ,IX2   ,ID2   ,SV2   ,
     >                  MATRT ,NBSD  ,NBST  )
*
*----------
*
*Purpose:
* To built equivalent surface array for translational symmetry.
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
* IPRINT  print level.
* IDT     translation direction.
* ILEV    geometry level.
* NFSUR   final number of surfaces.
* MAXMSP  maximum number of elements in MESH array.
* NS1     maximum number of surfaces in splitted geometry 1.
* NR1     maximum number of regions in splitted geometry 1.
* NS2     maximum number of surfaces in splitted geometry 2.
* NR2     maximum number of regions in splitted geometry 2.
* IEDIMG  geometries state vector.
* DAMESH  final mesh description for geometry.
* IX1     local indexing of surfaces/regions for geometry 1.
* ID1     surface identifier after symmetry for geometry 1.
* SV1     area/volume of regions for geometry 1.
* IX2     local indexing of surfaces/regions for geometry 2.
* ID2     surface identifier after symmetry for geometry 2.
* SV2     area/volume of regions for geometry 2.
*
*Parameters: output
* MATRT   reflection/transmission surface coupling array.
* NBSD    number of direct surfaces considered.
* NBST    number of translated surfaces found.
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*
*----------
*
      IMPLICIT         NONE
      INTEGER          NSTATE
      PARAMETER       (NSTATE=40)
*----
*  Subroutine arguments
*----
      INTEGER          IPRINT,IDT,ILEV,NFSUR,MAXMSP,NS1,NR1,NS2,NR2,
     >                 IEDIMG(NSTATE,2)
      DOUBLE PRECISION DAMESH(-1:MAXMSP,4,2)
      INTEGER          IX1(5,-NS1:NR1),ID1(NS1)
      DOUBLE PRECISION SV1(-NS1:NR1)
      INTEGER          IX2(5,-NS2:NR2),ID2(NS2)
      DOUBLE PRECISION SV2(-NS2:NR2)
      INTEGER          MATRT(NFSUR),NBSD,NBST
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTETH')
*----
*  Local variables
*----
      INTEGER          ISVD,ISVT,
     >                 ISD,IST,IPU,IPV,IPW,IDIR
      DOUBLE PRECISION D1B,D1T,D2B,D2T
*----
*  Processing starts:
*  print routine openning output header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6000) NAMSBR
        WRITE(IOUT,6010) ILEV
      ENDIF
*----
*  Assume mesh in U, V and W identical and stored in element 1
*  of DAMESH
*----
      IDIR=1
      IF(IDT .EQ. 3) THEN
*----
*  Scan for bottom surfaces
*----
        DO ISVD=1,IEDIMG(9,1)
          IF(IX1(IDT,-ISVD) .EQ. -1) THEN
            NBSD=NBSD+1
*----
*  This surface in the good direction
*----
            ISD=ID1(ISVD)
            IF(ISD .GT. 0) THEN
*----
*  This is an external surface/locate mesh position
*----
              IPU=IX1(1,-ISVD)
              IPV=IX1(2,-ISVD)
              IPW=IX1(4,-ISVD)
*----
*  Scan IX2 for top surfaces
*----
              DO ISVT=1,IEDIMG(9,2)
                IF(IX2(IDT,-ISVT) .EQ. -2) THEN
                  IST=ID2(ISVT)
                  IF(IST .GT. 0) THEN
*----
*  This is an external surface
*  test if mesh position is compatible with
*  direct geometry
*----
                    IF(IPU .EQ. IX2(1,-ISVT) .AND.
     >                 IPV .EQ. IX2(2,-ISVT) .AND.
     >                 IPW .EQ. IX2(4,-ISVT) ) THEN
*----
*  This should be the translated surface we are seeking
*  Test if area and dimensions are compatible
*----
                      IF(SV1(-ISVD) .NE. SV2(-ISVT))
     >CALL XABORT(NAMSBR//': Translated surfaces are invalid')
                      IF(IPU .GT. 0) THEN
                        D1B=DAMESH(IPU-1,IDIR,1)-DAMESH(0,IDIR,1)
                        D1T=DAMESH(IPU,IDIR,1)-DAMESH(0,IDIR,1)
                        D2B=DAMESH(IPU-1,IDIR,2)-DAMESH(0,IDIR,2)
                        D2T=DAMESH(IPU,IDIR,2)-DAMESH(0,IDIR,2)
                        IF(D1B .NE. D2B .AND. D1T .NE. D2T)
     >CALL XABORT(NAMSBR//': U mesh for translation is invalid')
                      ENDIF
                      IF(IPV .GT. 0) THEN
                        D1B=DAMESH(IPV-1,IDIR,1)-DAMESH(0,IDIR,1)
                        D1T=DAMESH(IPV,IDIR,1)-DAMESH(0,IDIR,1)
                        D2B=DAMESH(IPV-1,IDIR,2)-DAMESH(0,IDIR,2)
                        D2T=DAMESH(IPV,IDIR,2)-DAMESH(0,IDIR,2)
                        IF(D1B .NE. D2B .AND. D1T .NE. D2T)
     >CALL XABORT(NAMSBR//': V mesh for translation is invalid')
                      ENDIF
                      IF(IPW .GT. 0) THEN
                        D1B=DAMESH(IPW-1,IDIR,1)-DAMESH(0,IDIR,1)
                        D1T=DAMESH(IPW,IDIR,1)-DAMESH(0,IDIR,1)
                        D2B=DAMESH(IPW-1,IDIR,2)-DAMESH(0,IDIR,2)
                        D2T=DAMESH(IPW,IDIR,2)-DAMESH(0,IDIR,2)
                        IF(D1B .NE. D2B .AND. D1T .NE. D2T)
     >CALL XABORT(NAMSBR//': W mesh for translation is invalid')
                      ENDIF
*----
*  Everything seems all right
*  couple surfaces
*----
                      NBST=NBST+1
                      MATRT(ISD)=IST
                      MATRT(IST)=ISD
                      GO TO 105
                    ENDIF
                  ENDIF
                ENDIF
              ENDDO
*----
*  Could not find translated surface
*----
              WRITE(IOUT,9000) ISD
 105          CONTINUE
            ENDIF
          ENDIF
        ENDDO
      ELSE
        CALL XABORT(NAMSBR//': Translation BC for hexagonal faces '//
     >'not programmed yet ')
      ENDIF
*----
*  Processing finished:
*  print routine closing output header if required
*  and return
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6010 FORMAT('Geometry level = ',I5)
 9000 FORMAT(' ***** Warning ***** '/
     >       '       Translated surface for ',I5,1X,'is absent')
      END
