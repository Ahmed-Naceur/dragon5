*DECK NXTETS
      SUBROUTINE NXTETS(IPRINT,IDT   ,ILEV  ,NFSUR ,MAXMSP,
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
      PARAMETER       (IOUT=6,NAMSBR='NXTETS')
*----
*  Local variables
*----
      INTEGER          IDP1,IDP2,IDPR,IGT,ISVD,ISVT,
     >                 ISD,IST,IP1,IP2,IPR,KLEV
      DOUBLE PRECISION D1B,D1T,D2B,D2T
*----
*  Processing starts:
*  print routine openning output header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
*----
*  Scan IX1 for bottom surfaces
*----
      IDP1=MOD(IDT,3)+1
      IDP2=MOD(IDT+1,3)+1
      IDPR=4
      IGT=2
      NBSD=0
      NBST=0
      KLEV=ILEV
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
            IP1=IX1(IDP1,-ISVD)
            IP2=IX1(IDP2,-ISVD)
            IPR=IX1(IDPR,-ISVD)
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
                  IF(IP1 .EQ. IX2(IDP1,-ISVT) .AND.
     >               IP2 .EQ. IX2(IDP2,-ISVT) .AND.
     >               IPR .EQ. IX2(IDPR,-ISVT) ) THEN
*----
*  This should be the translated surface we are seeking
*  Test if area and dimensions are compatible
*----
                    IF(SV1(-ISVD) .NE. SV2(-ISVT))
     >CALL XABORT(NAMSBR//': Translated surfaces are invalid')
                    IF(IP1 .GT. 0) THEN
                      D1B=DAMESH(IP1-1,IDP1,1)-DAMESH(0,IDP1,1)
                      D1T=DAMESH(IP1,IDP1,1)-DAMESH(0,IDP1,1)
                      D2B=DAMESH(IP1-1,IDP1,2)-DAMESH(0,IDP1,2)
                      D2T=DAMESH(IP1,IDP1,2)-DAMESH(0,IDP1,2)
                      IF(D1B .NE. D2B .AND. D1T .NE. D2T) THEN
                        WRITE(IOUT,9001) D1B,D1T,D2B,D2T,
     >                   50.D0*(D1B-D2B)/(D1B+D2B),
     >                   50.D0*(D1T-D2T)/(D1T+D2T)
                         CALL XABORT(NAMSBR//
     >                   ': First mesh for translation is invalid')
                      ENDIF
                    ENDIF
                    IF(IP2 .GT. 0) THEN
                      D1B=DAMESH(IP2-1,IDP2,1)-DAMESH(0,IDP2,1)
                      D1T=DAMESH(IP2,IDP2,1)-DAMESH(0,IDP2,1)
                      D2B=DAMESH(IP2-1,IDP2,2)-DAMESH(0,IDP2,2)
                      D2T=DAMESH(IP2,IDP2,2)-DAMESH(0,IDP2,2)
                      IF(D1B .NE. D2B .AND. D1T .NE. D2T) THEN
                        WRITE(IOUT,9001) D1B,D1T,D2B,D2T,
     >                   50.D0*(D1B-D2B)/(D1B+D2B),
     >                   50.D0*(D1T-D2T)/(D1T+D2T)
                        CALL XABORT(NAMSBR//
     >                  ': Second mesh for translation is invalid')
                      ENDIF
                    ENDIF
                    IF(IPR .GT. 0) THEN
                      D1B=DAMESH(IPR-1,IDPR,1)
                      D1T=DAMESH(IPR,IDPR,1)
                      D2B=DAMESH(IPR-1,IDPR,2)
                      D2T=DAMESH(IPR,IDPR,2)
                      IF(D1B .NE. D2B .AND. D1T .NE. D2T) THEN
                        WRITE(IOUT,9001) D1B,D1T,D2B,D2T,
     >                  50.D0*(D1B-D2B)/(D1B+D2B),
     >                  50.D0*(D1T-D2T)/(D1T+D2T)
                        CALL XABORT(NAMSBR//
     >                  ': Radial mesh for translation is invalid')
                      ENDIF
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
 105        CONTINUE
          ENDIF
        ENDIF
      ENDDO
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
 9000 FORMAT(' ***** Warning ***** '/
     >       '       Translated surface for ',I5,1X,'is absent')
 9001 FORMAT(6F20.10)
      END
