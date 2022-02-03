*DECK NXTRPS
      SUBROUTINE NXTRPS(IPTRK ,IPRINT,ITYPG ,IGEO  ,ILEV  ,
     >                  NREG  ,NSUR  ,NSURN ,IDFEX ,
     >                  INDXSR,DHPIN ,DCMESH,NASUR ,IDSUR )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Renumber pin cluster surfaces.
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
* ITYPG   type of geometry.
* IGEO    geometry number.
* ILEV    geometry level.
* NREG    maximum number of regions in split geometry.
* NSUR    maximum number of surfaces in split geometry.
* NSURN   number of surfaces in splitted geometry after symmetry.
* IDFEX   flag to identify surface to consider.
* INDXSR  local indexing of surfaces/regions.
* DHPIN   pins height.
* DCMESH  cell dimensions.
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
      INTEGER          IPRINT,ITYPG,IGEO,ILEV,
     >                 NREG,NSUR,NSURN
      INTEGER          IDFEX(0:8),INDXSR(5,-NSUR:NREG)
      DOUBLE PRECISION DHPIN,DCMESH(3,2)
      INTEGER          NASUR,IDSUR(NSUR)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTRPS')
      DOUBLE PRECISION DZERO,DONE,DTWO
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0,DTWO=2.0D0)
*----
*  Local variables
*----
      INTEGER          IDGPP,IAS(2),IDIRS(2),ISUR,ID,IND,
     >                 IDIRC,INV,LSTSUR,IDS,INS
      CHARACTER        NAMREC*12
      DOUBLE PRECISION DCTOP,DCBOT,DOFF,DPBOT,DPTOP
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
*  Select main pin direction
*----
      IDGPP=3
      IF(ITYPG .EQ. 10  .OR. ITYPG .EQ. 21) THEN
        IDGPP=1
      ELSE IF(ITYPG .EQ. 11 .OR. ITYPG .EQ. 22) THEN
        IDGPP=2
      ENDIF
*----
*  Find if pin reaches bottom or top surface
*----
      DCTOP=DCMESH(IDGPP,1)/DTWO
      DCBOT=-DCTOP
      DOFF=DCMESH(IDGPP,2)
      DPBOT=DCMESH(IDGPP,2)-DHPIN/DTWO
      DPTOP=DPBOT+DHPIN
      IAS(1)=0
      IF(DPBOT .LE. DCBOT .AND. DPTOP .GE. DCBOT) IAS(1)=1
      IAS(2)=0
      IF(DPBOT .LE. DCTOP .AND. DPTOP .GE. DCTOP) IAS(2)=1
      IDIRS(2)=2*IDGPP
      IDIRS(1)=IDIRS(2)-1
*----
*  get rid of radial surfaces
*----
      DO ISUR=1,NSUR
        IF(INDXSR(4,-ISUR) .EQ. -2) THEN
          IDSUR(ISUR)=0
        ENDIF
      ENDDO
*----
*  Get rid of botton and top surfaces if not used
*----
      DO ID=1,2
        IND=-ID
        IDIRC=IDIRS(ID)
        IF(IAS(ID) .EQ. 0 .OR. IDFEX(IDIRC) .EQ. 0) THEN
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
        WRITE(IOUT,6014)  (IDSUR(IDS),IDS=1,NSUR)
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
