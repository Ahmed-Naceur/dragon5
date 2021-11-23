*DECK NXTRCS
      SUBROUTINE NXTRCS(IPTRK ,IPRINT,IGEO  ,ILEV  ,
     >                  NREG  ,NSUR  ,NSURN ,IDFEX ,
     >                  INDXSR,NASUR ,IDSUR )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Renumber cell surfaces.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s):
* G. Marleau
*
*Parameters: input
* IPTRK   pointer to the TRACKING data structure.
* IPRINT  intermediate printing level for output.
* IGEO    geometry number.
* ILEV    geometry level.
* NREG    maximum number of regions in splitted geometry.
* NSUR    maximum number of surfaces in splitted geometry.
* NSURN   number of surfaces in splitted geometry after symmetry.
* IDFEX   flag to identify surface to consider
*         (see NXTCUA for Cartesion geometry
*         and NXTHUA for hexagonal geometry).
* INDXSR  local indexing of surfaces/regions.
*
*Parameters: input/output
* NASUR   last surcace number considered.
* IDSUR   surface identifier after symmetry.
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*
*-----------------------------------------------------------------------
*
      USE              GANLIB
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      TYPE(C_PTR)      IPTRK
      INTEGER          IPRINT,IGEO,ILEV,
     >                 NREG,NSUR,NSURN
      INTEGER          IDFEX(0:10),INDXSR(5,-NSUR:NREG)
      INTEGER          NASUR,IDSUR(NSUR)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTRCS')
      DOUBLE PRECISION DZERO,DONE,DTWO
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0,DTWO=2.0D0)
*----
*  Local variables
*----
      INTEGER          IDGPP,ISUR,ID,IND,INV,LSTSUR,IDS,INS
      CHARACTER        NAMREC*12
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INREN
*----
*  Data
*----
      CHARACTER        CLEV(2)*1
      SAVE             CLEV
      DATA             CLEV /'C','P'/
*----
*  Processing starts:
*  print routine openning output header if required
*  and initialize various parameters.
*----
      ALLOCATE(INREN(NSURN))
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6000) NAMSBR
        WRITE(IOUT,6012)
        WRITE(IOUT,6014)  (IDSUR(IDS),IDS=1,NSUR)
      ENDIF
*----
*  Get rid of surfaces not used
*----
      DO ID=1,10
        IF(IDFEX(ID) .EQ. 0) THEN
          IDGPP=(ID+1)/2
          IND=-(MOD(ID-1,2)+1)
          DO ISUR=1,NSUR
            IF(INDXSR(IDGPP,-ISUR) .EQ. IND) THEN
              IDSUR(ISUR)=0
            ENDIF
          ENDDO
        ENDIF
      ENDDO
*----
*  Renumber surfaces
*----
      CALL XDISET(INREN,NSURN,0)
      INV=0
      DO ISUR=1,NSUR
        IND=IDSUR(ISUR)
        IF(IND .GT. 0) THEN
          INV=INV+1
          INREN(IND)=INV
        ENDIF
      ENDDO
      LSTSUR=INV+NASUR
      DO ISUR=1,NSUR
        IDS=IDSUR(ISUR)
        IF(IDS .NE. 0) THEN
          INS=INREN(ABS(IDS))
          IF(INS .NE. 0) INS=INS+NASUR
          IF(IDS .LT. 0) THEN
            IDSUR(ISUR)=-INS
          ELSE
            IDSUR(ISUR)=INS
          ENDIF
        ENDIF
      ENDDO
      NASUR=LSTSUR
      WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEO,'SID'
      CALL LCMPUT(IPTRK,NAMREC,NSUR,1,IDSUR)
*----
*  Processing finished:
*  print routine closing output header if required
*  and return
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6013)
        WRITE(IOUT,6014) (IDSUR(IDS),IDS=1,NSUR)
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      DEALLOCATE(INREN)
      RETURN
*----
*  FORMATS
*----
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6012 FORMAT(' Original surfaces ID')
 6013 FORMAT(' Final surfaces ID')
 6014 FORMAT(5I15)
      END
