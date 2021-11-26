*DECK NXTCVS
      SUBROUTINE NXTCVS(IPTRK ,IPRINT,NDIM  ,ITYPBC,NBOCEL,
     >                  NFSUR ,NFREG ,MXGSUR,MXGREG,
     >                  KEYMRG,MATALB,SURVOL)
*
*----------
*
*Purpose:
* To compute final surfaces and volumes for geometry
* and to create the EXCELL type MATALB and KEYMRG vector.
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
* IPTRK   pointer to the TRACKING data structure.
* IPRINT  print level.
* NDIM    problem dimensions.
* ITYPBC  type of boundary conditions where:
*         =0 for geometry with Cartesian boundaries;
*         =1 for geometry with annular boundary;
*         =2 for geometry with hexagonal boundary.
* NBOCEL  number of cells in original geometry.
* NFSUR   final number of surfaces.
* NFREG   final number of regions.
* MXGSUR  maximum number of surfaces for any geometry.
* MXGREG  maximum number of region for any geometry.
*
*Parameters: output
* KEYMRG  global merging vector.
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
      USE              GANLIB
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      TYPE(C_PTR)      IPTRK
*      INTEGER          IPTRK
      INTEGER          IPRINT,NDIM,ITYPBC,
     >                 NBOCEL,NFSUR,NFREG,MXGSUR,MXGREG
      INTEGER          KEYMRG(-NFSUR:NFREG),MATALB(-NFSUR:NFREG,2)
      DOUBLE PRECISION SURVOL(-NFSUR:NFREG)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTCVS')
      INTEGER          NSTATE
      PARAMETER       (NSTATE=40)
      DOUBLE PRECISION DZERO,DONE,DTWO
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0,DTWO=2.0D0)
*----
*  Local variables
*----
      INTEGER          ICEL,ICLS,ILEV,ISV,IGEO
      INTEGER          NREG,NSUR,NBGCLS,IGCLS,NUNK,MXRUNK
      INTEGER          IEDIMX(NSTATE),IEDIMP(NSTATE)
      CHARACTER        NAMREC*12
      DOUBLE PRECISION DFACC,DFACP
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDREG,IDSUR,MIX,MIXH
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: INDXSR
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: SVSGEO
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
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      NUNK=NFSUR+NFREG+1
      MXRUNK=MXGSUR+MXGREG+1
      CALL XDDSET(SURVOL,NUNK,DZERO)
      CALL XDISET(MATALB,NUNK*2,0)
*----
*  Here there are no merge
*----
      DO ISV=-NFSUR,NFREG
        KEYMRG(ISV)=ISV
      ENDDO
      DO ICEL=1,NBOCEL
        ILEV=1
        IGEO=ICEL
        WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEO,'DIM'
        CALL XDISET(IEDIMX,NSTATE,0)
        CALL LCMGET(IPTRK,NAMREC,IEDIMX)
        NREG=IEDIMX(8)
        NSUR=IEDIMX(9)
        NBGCLS=IEDIMX(16)
        IGCLS=IEDIMX(17)-1
*----
*   Get MIXTURE
*----
        ALLOCATE(MIX(NREG),MIXH(NREG),INDXSR(5,-NSUR:NREG))
        ALLOCATE(IDREG(NREG),IDSUR(NSUR))
        ALLOCATE(SVSGEO(2*(NSUR+NREG+1)))
      
        WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEO,'MIX'
        CALL LCMGET(IPTRK,NAMREC,MIX)
        WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEO,'HOM'
        CALL LCMGET(IPTRK,NAMREC,MIXH)
*----
*  Get INDEX and SURVOL for pin
*----
        WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEO,'VSI'
        CALL LCMGET(IPTRK,NAMREC,INDXSR)
        WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEO,'VSE'
        CALL LCMGET(IPTRK,NAMREC,SVSGEO)
*----
*   Get IDREG and IDSUR
*----
        WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEO,'RID'
        CALL LCMGET(IPTRK,NAMREC,IDREG)
        WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEO,'SID'
        CALL LCMGET(IPTRK,NAMREC,IDSUR)
        DFACC=DBLE(IEDIMX(19))
        CALL NXTAVS(IPRINT,NDIM  ,ITYPBC,NFSUR ,NFREG ,NSUR  ,
     >              NREG  ,MIX   ,MIXH  ,INDXSR,IDSUR ,IDREG ,
     >              SVSGEO,DFACC ,MATALB,SURVOL)
        DEALLOCATE(SVSGEO,IDSUR,IDREG,INDXSR,MIXH,MIX)
        IF(NBGCLS .NE. 0) THEN
          ILEV=2
          DO ICLS=1,NBGCLS
            IGEO=IGCLS+ICLS
            WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEO,'DIM'
            CALL XDISET(IEDIMP,NSTATE,0)
            CALL LCMGET(IPTRK,NAMREC,IEDIMP)
            NREG=IEDIMP(8)
            NSUR=IEDIMP(9)
            ALLOCATE(MIX(NREG),MIXH(NREG),INDXSR(5,-NSUR:NREG))
            ALLOCATE(IDREG(NREG),IDSUR(NSUR))
            ALLOCATE(SVSGEO(2*(NSUR+NREG+1)))
*----
*   Get MIXTURE
*----
            WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEO,'MIX'
            CALL LCMGET(IPTRK,NAMREC,MIX)
            WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEO,'HOM'
            CALL LCMGET(IPTRK,NAMREC,MIXH)
*----
*   Get INDEX and SURVOL for cell
*----
            WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEO,'VSI'
            CALL LCMGET(IPTRK,NAMREC,INDXSR)
            WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEO,'VSE'
            CALL LCMGET(IPTRK,NAMREC,SVSGEO)
*----
*   Get IDREG and IDSUR
*----
            WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEO,'RID'
            CALL LCMGET(IPTRK,NAMREC,IDREG)
            WRITE(NAMREC,'(A1,I8.8,A3)') CLEV(ILEV),IGEO,'SID'
            CALL LCMGET(IPTRK,NAMREC,IDSUR)
            DFACP=DFACC*DBLE(IEDIMP(17))
            CALL NXTAVS(IPRINT,NDIM  ,ITYPBC,NFSUR ,NFREG ,NSUR  ,
     >                  NREG  ,MIX   ,MIXH  ,INDXSR,IDSUR ,IDREG ,
     >                  SVSGEO,DFACP ,MATALB,SURVOL)
            DEALLOCATE(SVSGEO,IDSUR,IDREG,INDXSR,MIXH,MIX)
          ENDDO
        ENDIF
      ENDDO
*----
*  Save records on IPTRK
*----
      CALL LCMPUT(IPTRK ,'KEYMRG      ',NUNK  ,1,KEYMRG)
      CALL LCMPUT(IPTRK ,'MATALB      ',NUNK  ,1,MATALB(-NFSUR,1))
      CALL LCMPUT(IPTRK ,'HOMMATALB   ',NUNK  ,1,MATALB(-NFSUR,2))
      CALL LCMPUT(IPTRK ,'SAreaRvolume',NUNK,4,SURVOL)
*----
*  Processing finished:
*  print routine closing output header if required
*  and return
*----
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
      END
