*DECK NXTBCG
      SUBROUTINE NXTBCG(IPGEO ,IPTRK ,IPRINT,NDIM  ,ITYPBC,IDIRG ,
     >                  IDIAG ,ISAXIS,IHSYM ,ILEAK ,IPRISM )
*
*----------
*
*Purpose:
* To read boundary conditions and symmetries and verify
* if information is consistent. Also prepare axis for
* symmetry identification.
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
* IPTRK   pointer to the TRACKING data structure in
*         update or creation mode.
* IPRINT  print level.
* NDIM    problem dimensions.
* ITYPBC  type of boundary conditions where:
*         =0 for geometry with Cartesian
*         boundaries;
*         =1 for geometry with annular
*         boundary;
*         =2 for geometry with hexagonal
*         boundary;
* IDIRG   geometry main direction:
*         =1 for $X-Y-Z$ geometry;
*         =2 for $Y-Z-X$ geometry;
*         =3 for $Z-X-Y$ geometry.
* IPRISM  projection axis for prismatic tracking.
*
*Parameters: output
* IDIAG   the diagonal symmetry flag where:
*         =-1 indicates X- DIAG symmetry;
*         = 1 indicates X+ Y- DIAG symmetry;
*         = 0 indicates no DIAG symmetry.
* ISAXIS  symmetry vector for each direction.
* IHSYM   hexagonal symmetry option where:
*         = 0 geometry is not hexagonal;
*         = 1 for S30;
*         = 2 for SA60;
*         = 3 for SB60;
*         = 4 for S90;
*         = 5 for R120;
*         = 6 for R180;
*         = 7 for SA180;
*         = 8 for SB180;
*         = 9 for COMPLETE;
*         =10 for R60.
* ILEAK   leakage option where:
*         =1 indicates that there is no out of cell leakage;
*         =0 indicates that there  no out of cell leakage.
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*
*----------
*
      USE              GANLIB
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      TYPE(C_PTR)      IPGEO,IPTRK
      INTEGER          IPRINT,NDIM,ITYPBC,IDIRG
      INTEGER          IDIAG,ISAXIS(3),IHSYM,ILEAK,IPRISM
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTBCG')
      INTEGER          MAXDIM,NGBC
      PARAMETER       (MAXDIM=3,NGBC=6)
*----
*  Local variables
*----
      INTEGER          IDIRP,IT3,ISUR,IDIR,ISCOMP
      INTEGER          ICODE(NGBC),NCODE(NGBC),JCODE(NGBC)
      REAL             ZCODE(NGBC)
      INTEGER          IVBC(NGBC),MRGSUR(NGBC)
*----
*  Data
*----
      CHARACTER        CDIR(MAXDIM)*1
      SAVE             CDIR
      DATA             CDIR /'X','Y','Z'/
*----
*  Processing starts:
*  print routine opening header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      IHSYM=0
*----
*  Read boundary conditions from geometry
*----
      DO ISUR=1,NGBC
        ICODE(ISUR)=-ISUR
      ENDDO
      CALL LCMGET(IPGEO ,'ICODE',JCODE)
      CALL LCMGET(IPGEO ,'NCODE',NCODE)
      IF (IPRISM.NE.0) THEN
*     special case of a prismatic tracking:
*     the geometry is not unfolded for the symmetries along the projection axis
         CALL LCMPUT(IPTRK,'NCODE',6,1,NCODE)
         DO IDIRG=1,2
            IF ((NCODE(2*(IPRISM-1)+IDIRG).EQ.5).OR.
     1          (NCODE(2*(IPRISM-1)+IDIRG).EQ.10)) THEN
               NCODE(2*(IPRISM-1)+IDIRG)=2
            ENDIF
         ENDDO
      ENDIF
      CALL LCMGET(IPGEO ,'ZCODE',ZCODE)
*----
*  Validate boundary conditions to process
*----
      IDIRP=MOD(IDIRG+1,3)+1
      IT3=NDIM/3
      IF(ITYPBC .EQ. 0) THEN
        DO  ISUR=1,4
          IVBC(ISUR)=1
          MRGSUR( ISUR)=ISUR
        ENDDO
        DO ISUR=5,NGBC
          IVBC(ISUR)=IT3
          MRGSUR( ISUR)=ISUR
        ENDDO
      ELSE IF(ITYPBC .EQ. 1) THEN
        DO ISUR=1,NGBC
          IVBC(ISUR)=0
        ENDDO
        IF(IDIRG .EQ. 1) THEN
          IVBC(2)=1
          IF(NDIM .EQ. 3) THEN
            IVBC(5)=1
            IVBC(6)=1
          ENDIF
        ELSE IF(IDIRG .EQ. 2) THEN
          IVBC(4)=1
          IF(NDIM .EQ. 3) THEN
            IVBC(1)=1
            IVBC(2)=1
          ENDIF
        ELSE IF(IDIRG .EQ. 3) THEN
          IVBC(6)=1
          IF(NDIM .EQ. 3) THEN
            IVBC(3)=1
            IVBC(4)=1
          ENDIF
        ENDIF
      ELSE IF(ITYPBC .EQ. 2) THEN
        CALL LCMGET(IPGEO ,'IHEX',IHSYM)
        IF(IHSYM .NE. 9) CALL XABORT(NAMSBR//': only COMPLETE '//
     >  ' hexagonal symmetry option programed in NXT:.')
        DO ISUR=1,NGBC
          IVBC(ISUR)=0
        ENDDO
        IVBC(1)=1
        IF(NDIM .EQ. 3) THEN
          IVBC(5)=1
          IVBC(6)=1
        ENDIF
      ENDIF
      DO ISUR=1,NGBC
        MRGSUR(ISUR)=ISUR
      ENDDO
*----
*  Find pairs of diagonal and translation B.C.
*----
      IDIAG=0
      ISCOMP=0
      DO ISUR=1,NGBC
        IF(IVBC(ISUR) .EQ. 1) THEN
          IF(JCODE(ISUR) .NE. 0 ) THEN
            ICODE(ISUR)=JCODE(ISUR)
            ZCODE(ISUR)= 1.0
          ELSE IF(NCODE(ISUR) .EQ. 0) THEN
            CALL XABORT(NAMSBR//
     >      ': A boundary condition is missing.')
          ENDIF
          IF(NCODE(ISUR) .EQ. 2) THEN
            ZCODE(ISUR)= 1.0
          ELSE IF(NCODE(ISUR) .EQ. 3) THEN
            IDIAG=IDIAG+1
          ELSE IF(NCODE(ISUR) .EQ. 4)THEN
            ISCOMP=ISCOMP+1
            ZCODE(ISUR)=1.0
          ELSE IF(NCODE(ISUR) .EQ. 6 ) THEN
            NCODE(ISUR)= 1
          ELSE IF(NCODE(ISUR) .EQ. 7 .OR.
     >            NCODE(ISUR) .EQ. 8 .OR.
     >            NCODE(ISUR) .EQ. 9 .OR.
     >            NCODE(ISUR) .GE. 11 ) THEN
            CALL XABORT(NAMSBR//
     >      ': An invalid boundary condition detected.')
          ENDIF
        ENDIF
      ENDDO
*----
*  Analyse DIAG boundary conditions
*  Only X+ DIAG Y- DIAG or X- DIAG Y+ DIAG permitted
*----
      IF(IDIAG .GT. 0) THEN
        IF(ITYPBC .NE. 0) CALL XABORT(NAMSBR//
     >  ': DIAG BC permitted only for Cartesian geometries')
        IF(IDIAG .NE. 2) CALL XABORT(NAMSBR//
     >  ': Only one pair of DIAG boundary conditions permitted')
        IF((NCODE(2) .EQ. 3) .AND. (NCODE(3) .EQ. 3)) THEN
          MRGSUR(2)= 4
          MRGSUR(3)= 1
          NCODE(2)=  NCODE(4)
          NCODE(3)=  NCODE(1)
          ICODE(2)=  ICODE(4)
          ICODE(3)=  ICODE(1)
          ZCODE(2)=  ZCODE(4)
          ZCODE(3)=  ZCODE(1)
          IDIAG=1
        ELSE IF((NCODE(1) .EQ. 3) .AND. (NCODE(4) .EQ. 3)) THEN
          MRGSUR(1)= 3
          MRGSUR(4)= 2
          NCODE(1)=  NCODE(3)
          NCODE(4)=  NCODE(2)
          ICODE(1)=  ICODE(3)
          ICODE(4)=  ICODE(2)
          ZCODE(1)=  ZCODE(3)
          ZCODE(4)=  ZCODE(2)
          IDIAG=-1
        ELSE
          CALL XABORT(NAMSBR//
     >    ': Only (X+ DIAG Y- DIAG) or (X- DIAG Y+ DIAG) permitted')
        ENDIF
      ENDIF
*----
*  Analyse TRAN boundary conditions
*  Only X- TRAN X+ TRAN,  Y- TRAN Y+ TRAN and
*  Z- TRAN Z+ TRAN permitted
*----
      DO IDIR=1,MAXDIM
        ISAXIS(IDIR)=0
      ENDDO
      IF(ISCOMP .GT. 0) THEN
        IF(ITYPBC .NE. 0) CALL XABORT(NAMSBR//
     >  ': TRAN BC permitted only for Cartesian geometries')
        IF(MOD(ISCOMP,2) .EQ. 1) CALL XABORT(NAMSBR//
     >  ': TRAN boundary conditions must come in pairs')
        DO IDIR=1,MAXDIM
          ISUR=2*IDIR
          IF(IVBC(ISUR) .EQ. 1) THEN
            IF(NCODE(ISUR) .EQ. 4 .AND. NCODE(ISUR-1) .EQ. 4) THEN
              MRGSUR(ISUR  )=ISUR-1
              MRGSUR(ISUR-1)=ISUR
              ISCOMP=ISCOMP-2
              ISAXIS(IDIR)=3
            ENDIF
          ENDIF
        ENDDO
        IF(ISCOMP .NE. 0) CALL XABORT(NAMSBR//
     >  ': Illegal pairing of TRAN boundary conditions')
      ENDIF
*----
*  Analyse SYME and SSYM boundary conditions
*----
      DO ISUR=1,NGBC
        IF(IVBC(ISUR) .EQ. 1) THEN
          IDIR=(ISUR+1)/2
          IF(NCODE(ISUR) .EQ. 5)THEN
            IF(MOD(ISUR,2) .EQ. 0) THEN
              ISCOMP=ISUR-1
              ISAXIS(IDIR)=1
            ELSE
              ISCOMP=ISUR+1
              ISAXIS(IDIR)=-1
            ENDIF
            IF(NCODE(ISCOMP) .NE. 1 .AND. NCODE(ISCOMP) .NE. 2 .AND.
     >         NCODE(ISCOMP) .NE. 6) CALL XABORT(NAMSBR//
     >      ': Invalid combination for SYME or SSYM symmetry')
            MRGSUR(ISUR)=ISCOMP
            ZCODE(ISUR)=ZCODE(ISCOMP)
            ICODE(ISUR)=ICODE(ISCOMP)
            NCODE(ISUR)=NCODE(ISCOMP)
          ELSE IF(NCODE(ISUR) .EQ. 10)THEN
            IF(MOD(ISUR,2) .EQ. 0) THEN
              ISCOMP=ISUR-1
              ISAXIS(IDIR)=2
            ELSE
              ISCOMP=ISUR+1
              ISAXIS(IDIR)=-2
            ENDIF
            IF(NCODE(ISCOMP) .NE. 1 .AND. NCODE(ISCOMP) .NE. 2 .AND.
     >         NCODE(ISCOMP) .NE. 6) CALL XABORT(NAMSBR//
     >      ': Invalid combination for SYME or SSYM symmetry')
            MRGSUR(ISUR)=ISCOMP
            ZCODE(ISUR)=ZCODE(ISCOMP)
            ICODE(ISUR)=ICODE(ISCOMP)
            NCODE(ISUR)=NCODE(ISCOMP)
          ENDIF
        ENDIF
      ENDDO
      ILEAK=1
      DO ISUR=1,NGBC
        IF(IVBC(ISUR) .EQ. 1) THEN
          IF(ICODE(ISUR) .GT. 0) THEN
            ILEAK=0
          ELSE IF(ZCODE(ISUR) .NE. 1.0) THEN
            ILEAK=0
          ENDIF
        ENDIF
      ENDDO
*----
*  For combined DIAG/SYME symmetry
*  complete set of symmetry available
*  X/Y SYME -> Y/X SYME
*----
*      IF(IDIAG .EQ. -1) THEN
*        IF(ISAXIS(2) .EQ. -1) THEN
*          ISAXIS(1)=1
*        ELSE IF(ISAXIS(1) .EQ. 1) THEN
*          ISAXIS(2)=-1
*        ENDIF
*      ELSE IF(IDIAG .EQ. 1) THEN
*        IF(ISAXIS(2) .EQ. 1) THEN
*          ISAXIS(1)=-1
*        ELSE IF(ISAXIS(1) .EQ. -1) THEN
*          ISAXIS(2)=1
*        ENDIF
*      ENDIF
*----
*  Save boundary conditions on tracking
*----
      CALL LCMPUT(IPTRK ,'ALBEDO      ',NGBC,2,ZCODE)
      CALL LCMPUT(IPTRK ,'ICODE       ',NGBC,1,ICODE)
*----
*  Processing finished:
*  print routine closing header if required
*  and return
*----
      IF(IPRINT .GE. 10) THEN
        IF(IDIAG .EQ. 1) THEN
          WRITE(IOUT,6010) '(X+, Y-)'
        ELSE IF(IDIAG .EQ. -1) THEN
          WRITE(IOUT,6010) '(X+, Y-)'
        ENDIF
        DO ISUR=1,3
          IF(ISAXIS(ISUR) .EQ. -2) THEN
            WRITE(IOUT,6011) 'SSYM',CDIR(ISUR),'-'
          ELSE IF(ISAXIS(ISUR) .EQ. -1) THEN
            WRITE(IOUT,6011) 'SYME',CDIR(ISUR),'-'
          ELSE IF(ISAXIS(ISUR) .EQ.  1) THEN
            WRITE(IOUT,6011) 'SYME',CDIR(ISUR),'+'
          ELSE IF(ISAXIS(ISUR) .EQ.  2) THEN
            WRITE(IOUT,6011) 'SSYM',CDIR(ISUR),'+'
          ELSE IF(ISAXIS(ISUR) .EQ.  3) THEN
            WRITE(IOUT,6011) 'TRAN',CDIR(ISUR),' '
          ENDIF
        ENDDO
        DO ISUR=1,NGBC
          WRITE(IOUT,6012) ISUR,ICODE(ISUR),NCODE(ISUR),IVBC(ISUR),
     >                  MRGSUR(ISUR),ZCODE(ISUR)
        ENDDO
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6010 FORMAT('Diagonal symmetry : ',A8)
 6011 FORMAT('Symmetry : ',A4,1X,'on surface ',2A1)
 6012 FORMAT('BC[[',I10,']]={',4(I5,','),F20.10,'}')
      END
