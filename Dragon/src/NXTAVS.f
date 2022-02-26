*DECK NXTAVS
      SUBROUTINE NXTAVS(IPRINT,NDIM  ,ITYPBC,NFSUR ,NFREG ,NSUR  ,
     >                  NREG  ,MIX   ,MIXH  ,INDXSR,IDSUR ,IDREG ,
     >                  SVSGEO,DFACC ,MATALB,SURVOL)
*
*----------
*
*Purpose:
* To add current cell information to global
* surfaces and volumes for geometry.
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
* NDIM    problem dimensions.
* ITYPBC  type of boundary conditions where
*         =0 for geometry with Cartesianb oundaries;
*         =1 for geometry with annular boundary;
*         =2 for geometry with hexagonal boundary.
* NFSUR   final number of surfaces.
* NFREG   final number of regions.
* NSUR    maximum number of surfaces in splitted geometry.
* NREG    maximum number of regions in splitted geometry.
* MIX     geometry mixtures .
* MIXH    homogenization mixtures.
* INDXSR  local indexing of surfaces/regions.
* IDSUR   local surface identifier .
* IDREG   local region identifier.
* SVSGEO  area/volume of regions.
* DFACC   multiplication factor for surface and volume.
*
*Parameters: input/output
* MATALB  global mixture/albedo identification vector (including HMIX).
* SURVOL  global surface volume vector.
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
*----
*  Subroutine arguments
*----
      INTEGER          IPRINT,NDIM,ITYPBC,NFSUR,NFREG,NSUR,NREG
      INTEGER          MIX(NREG),MIXH(NREG),INDXSR(5,-NSUR:NREG),
     >                 IDREG(NREG),IDSUR(NSUR)
      DOUBLE PRECISION SVSGEO(-NSUR:NREG),DFACC
      INTEGER          MATALB(-NFSUR:NFREG,2)
      DOUBLE PRECISION SURVOL(-NFSUR:NFREG)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTAVS')
*----
*  Local variables
*----
      INTEGER          NDSCAN,ISV,IDSV,IFSV,ID,IDSA
      INTEGER          IDALB(10)
*----
*  Processing starts:
*  print routine openning output header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      NDSCAN=NDIM
      IF(ITYPBC .EQ. 2) THEN
        NDSCAN=5
        DO ISV=1,4
          IDALB(ISV)=-1
        ENDDO
        DO ISV=5,6
          IDALB(ISV)=-ISV
        ENDDO
        IDALB(7)= 0
        IDALB(8)=-2
        IDALB(9)=-1
        IDALB(10)=-1
      ELSE
        DO ISV=1,6
          IDALB(ISV)=-ISV
        ENDDO
        IDALB(7)=-1
        IDALB(8)=-2
        IDALB(9)=-1
        IDALB(10)=-1
      ENDIF
*----
*  Add surface contributions
*----
      DO ISV=1,NSUR
        IDSV=IDSUR(ISV)
        IF(IDSV .NE. 0) THEN
          IFSV=-ABS(IDSV)
          SURVOL(IFSV)=SURVOL(IFSV)+DFACC*SVSGEO(-ISV)
          IF(IDSV .GT. 0) THEN
            DO ID=1,NDSCAN
              IF(INDXSR(ID,-ISV) .LT. 0) THEN
                IDSA=2*(ID-1)-INDXSR(ID,-ISV)
                MATALB(IFSV,1)=IDALB(IDSA)
                MATALB(IFSV,2)=IDALB(IDSA)
                GO TO 105
              ENDIF
            ENDDO
*----
*  Albedo type not found
*----
            CALL XABORT(NAMSBR//': Albedo type not found')
 105        CONTINUE
          ENDIF
        ENDIF
      ENDDO
*----
*  Add volume contribution
*----
      DO ISV=1,NREG
        IDSV=IDREG(ISV)
        IF(IDSV .NE. 0) THEN
          IFSV=ABS(IDSV)
          SURVOL(IFSV)=SURVOL(IFSV)+DFACC*SVSGEO(ISV)
          IF(IDSV .GT. 0) THEN
            MATALB(IFSV,1)=MIX(ISV)
            MATALB(IFSV,2)=MIXH(ISV)
          ENDIF
        ENDIF
      ENDDO
*----
*  Processing finished:
*  print routine closing output header if required
*  and return
*----
      IF(IPRINT .GE. 100) THEN
        IF(IPRINT .GE. 500) THEN
          WRITE(IOUT,*) 'MATALB'
          WRITE(IOUT,'(10I10)') (MATALB(ISV,1),ISV=-NFSUR,NFREG)
          WRITE(IOUT,*) 'HOMMATALB'
          WRITE(IOUT,'(10I10)') (MATALB(ISV,2),ISV=-NFSUR,NFREG)
          WRITE(IOUT,*) 'SURVOL'
          WRITE(IOUT,'(1P,5E20.10)') (SURVOL(ISV),ISV=-NFSUR,NFREG)
        ENDIF
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
      END
