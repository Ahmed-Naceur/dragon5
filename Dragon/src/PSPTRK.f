*DECK PSPTRK
      SUBROUTINE PSPTRK(IPRINT,ISPSP ,ITYPE ,ICOLR ,IPTRKT,NAMFIL,
     >                  NAMLEG,NUNKNO,FLUX  )
*
*----------
*
*Purpose:
* To generate a POSTSCRIPT file containing a graphical description
* of a 2-D geometry from an EXCELL generated
* tracking data structure.
*
*Copyright:
* Copyright (C) 1999 Ecole Polytechnique de Montreal.
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):
* G. Marleau.
*
*Parameters: input
* IPRINT  print level.
* ISPSP   POSTSCRIPT file index.
* ITYPE   identifier for the type of graphics where:
*         =0 when the geometry is colored by region;
*         =1 when the geometry is colored by mixture;
*         =2 when the geometry is colored by flux
*            (one group);
*         =3 when the geometry is colored by flux
*            (multigroup);
*         =4 when the geometry is colored by mixture for
*            homogenization.
* ICOLR   color set used where:
*         =-4 HSB filling with no contour;
*         =-3 CYMK filling with no contour;
*         =-2 RGB filling with no contour;
*         =-1 BW filling with no contour;
*         = 0 no filling with contour;
*         = 1 BW filling with contour;
*         = 2 RGB filling with contour;
*         = 3 CMYK filling with contour;
*         = 4 HSB filling with  contour.
* IPTRKT  pointer to the TRACKING data structure.
* NAMFIL  geometry file name.
* NAMLEG  legend name.
* NUNKNO  number of flux unknowns.
*
*Parameters: temporary storage
* FLUX    flux storage array.
*
*----------
*
      USE          GANLIB
      IMPLICIT     NONE
*----
*  Subroutine arguments
*----
      TYPE(C_PTR)  IPTRKT
      INTEGER      IPRINT,ISPSP,ITYPE,ICOLR,NUNKNO
      CHARACTER    NAMFIL*12,NAMLEG*24
      REAL         FLUX(NUNKNO)
*----
*  Local parameters
*----
      INTEGER      IOUT
      CHARACTER    NAMSBR*6
      PARAMETER   (IOUT=6,NAMSBR='PSPTRK')
      INTEGER      NSTATE
      PARAMETER   (NSTATE=40)
      INTEGER      ILCMUP,ILCMDN
      PARAMETER   (ILCMUP=1,ILCMDN=2)
*----
*  Local variables
*----
      INTEGER      ISTATE(NSTATE),IPARAM(NSTATE),ITROP,
     >             NDIM,NVOL,NSUR,NSURX,NBAN,NUNK,NRT,MSROD,
     >             MAROD,NTOTCL,MAXR,NUNKT,NREGT,NNSUR
      REAL         COTE
      INTEGER      IEDIMG(NSTATE),ITYPBC,NBUCEL,NUCELL(3),
     >             MAXMSH,MAXMDH,MAXREG,NBTCLS,MAXPIN,MAXMSP,
     >             MAXRSP,NFSUR,NFREG,MXGREG
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KEYFLX,KEYMRG,MATALB,IUNFLD,
     > NRODS,NRODR,NRINFO,NXRI,MINDIM,MAXDIM,INDEX
      REAL, ALLOCATABLE, DIMENSION(:) :: COLRG,RAN,RODS,RODR,REMSH
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DGMESH
*----
*  Processing starts:
*  print routine openning output header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6000) NAMSBR
        WRITE(IOUT,6002) NAMFIL
      ENDIF
*----
*  Get state vector from tracking
*  and check if a graphical description
*  of the geometry is possible.
*----
      CALL XDISET(ISTATE,NSTATE,0)
      CALL LCMGET(IPTRKT,'STATE-VECTOR',ISTATE)
      IF(ITYPE .EQ. 2 .OR. ITYPE .EQ. 3 .OR.
     >   ITYPE .EQ. 5 .OR. ITYPE .EQ. 6) THEN
        NREGT=ISTATE(1)
        NUNKT=ISTATE(2)
        IF(NUNKNO .NE. NUNKT) CALL XABORT(NAMSBR//
     >  ': Tracking is not consistent with fluxes')
        ALLOCATE(KEYFLX(NREGT))
        CALL LCMGET(IPTRKT,'KEYFLX      ',KEYFLX)
      ELSE
        NREGT=1
        ALLOCATE(KEYFLX(NREGT))
        KEYFLX=0
      ENDIF
      ITROP=ISTATE(7)
      IF(ITROP .EQ. 4) THEN
*----
*  NXT processed geometry
*----
        CALL LCMSIX(IPTRKT,'NXTRecords  ',ILCMUP)
        CALL LCMGET(IPTRKT,'G00000001DIM',IEDIMG)
        NDIM=IEDIMG( 1)
        ITYPBC   =IEDIMG( 2)
        NBUCEL   =IEDIMG( 5)
        NUCELL(1)=IEDIMG(13)
        NUCELL(2)=IEDIMG(14)
        NUCELL(3)=IEDIMG(15)
        MAXMSH   =IEDIMG(16)
        MAXREG   =IEDIMG(17)
        NBTCLS   =IEDIMG(18)
        MAXPIN   =IEDIMG(19)
        MAXMSP   =IEDIMG(20)
        MAXRSP   =IEDIMG(21)
        NFSUR    =IEDIMG(22)
        NFREG    =IEDIMG(23)
        MXGREG   =IEDIMG(25)
        NSUR=NFSUR
        NVOL=NFREG
        NNSUR=-NSUR
        NUNK=NSUR+NVOL+1
        ALLOCATE(COLRG(4*NVOL))
        ALLOCATE(KEYMRG(NUNK),MATALB(NUNK))
        CALL XDISET(IEDIMG,NSTATE,0)
        IF(NDIM .EQ. 2) THEN
          IF(ITYPE .EQ. 4) THEN
            CALL LCMGET(IPTRKT,'HOMMATALB   ',MATALB)
          ELSE
            CALL LCMGET(IPTRKT,'MATALB      ',MATALB)
          ENDIF
          CALL LCMGET(IPTRKT,'KEYMRG      ',KEYMRG)
*----
*  Produce legend
*----
          CALL PSPLEG(IPRINT,ISPSP ,ITYPE ,ICOLR ,NNSUR ,NVOL  ,
     >                NAMLEG,NUNKNO,FLUX  ,NREGT ,
     >                MATALB,KEYMRG,KEYFLX,
     >                COLRG)
*----
*  Produce graphical description of geometry
*----
          NUNK=NFSUR+NFREG+1
          MAXMDH=MAX(MAXMSH,MAXMSP,MAXREG)
          ALLOCATE(IUNFLD(2*NBUCEL),DGMESH((MAXMDH+2)*4))
          CALL LCMGET(IPTRKT,'G00000001CUF',IUNFLD)
          CALL PSPNXT(IPRINT,ISPSP ,ICOLR ,IPTRKT,ITYPBC,MAXMDH,
     >                NDIM  ,NFSUR ,NFREG ,NUCELL,NBUCEL,
     >                MXGREG,MAXPIN,COLRG, IUNFLD,MATALB,DGMESH)
          DEALLOCATE(DGMESH,IUNFLD)
        ELSE
          WRITE(IOUT,9000)
        ENDIF
        CALL LCMSIX(IPTRKT,'NXTRecords  ',ILCMDN)
      ELSE
        CALL LCMSIX(IPTRKT,'EXCELL      ',1)
        CALL XDISET(IPARAM,NSTATE,0)
        CALL LCMGET(IPTRKT,'STATE-VECTOR',IPARAM)
        NDIM=IPARAM(1)
        NSUR=-IPARAM(2)
        NVOL=IPARAM(3)
        NSURX=IPARAM(4)
        NBAN=IPARAM(5)
        NUNK=IPARAM(6)
        ALLOCATE(COLRG(4*NVOL))
        ALLOCATE(KEYMRG(NUNK),MATALB(NUNK))
        CALL LCMGET(IPTRKT,'MATALB      ',MATALB)
        CALL LCMGET(IPTRKT,'KEYMRG      ',KEYMRG)
        IF(ITROP .EQ. 3) THEN
*----
*  EXCELL based CLUSTER geometries
*----
          NRT=IPARAM(7)
          MSROD=IPARAM(8)
          MAROD=IPARAM(9)
          ALLOCATE(NRODS(3*NRT),NRODR(NRT),NRINFO(2*NBAN),
     >    NXRI(NRT*NBAN))
          ALLOCATE(RAN(NBAN),RODS(2*NRT),RODR(MSROD*NRT))
          CALL LCMGET(IPTRKT,'RAN         ',RAN)
          IF(NSURX .EQ. 4)
     >    CALL LCMGET(IPTRKT,'COTE        ',COTE)
          CALL LCMGET(IPTRKT,'NRODS       ',NRODS)
          CALL LCMGET(IPTRKT,'RODS        ',RODS)
          CALL LCMGET(IPTRKT,'NRODR       ',NRODR)
          CALL LCMGET(IPTRKT,'RODR        ',RODR)
          CALL LCMGET(IPTRKT,'NRINFO      ',NRINFO)
          CALL LCMGET(IPTRKT,'NXRI        ',NXRI)
*----
*  Produce legend
*----
          CALL PSPLEG(IPRINT,ISPSP ,ITYPE ,ICOLR ,NSUR  ,NVOL  ,
     >                NAMLEG,NUNKNO,FLUX  ,NREGT ,
     >                MATALB,KEYMRG,KEYFLX,COLRG)
*----
*  Produce graphical description of geometry
*----
          CALL PSPXCG(IPRINT,ISPSP ,ICOLR ,NBAN  ,NRT   ,MSROD ,
     >                NSURX ,NSUR  ,NVOL  ,COTE  ,
     >                RAN   ,NRODS ,RODS  ,RODR  ,NRINFO,NRODR ,
     >                NXRI  ,KEYMRG,COLRG)
          DEALLOCATE(RODR,RODS,RAN)
          DEALLOCATE(NXRI,NRINFO,NRODR,NRODS)
        ELSE IF(ITROP .EQ. 2 ) THEN
*----
*  EXCELL based hexagonal geometries
*  Not available yet
*----
*          CALL PSPXHX(IPRINT,IPTRKT,TITREC)
          WRITE(IOUT,6001)
        ELSE IF(ITROP .EQ. 1 ) THEN
*----
*  EXCELL based Cartesian geometries
*----
          NTOTCL=NSURX
          MAXR=NBAN
          ALLOCATE(MINDIM(NTOTCL),MAXDIM(NTOTCL),INDEX(4*NUNK))
          ALLOCATE(REMSH(MAXR))
          CALL LCMGET(IPTRKT,'MINDIM      ',MINDIM)
          CALL LCMGET(IPTRKT,'MAXDIM      ',MAXDIM)
          CALL LCMGET(IPTRKT,'INDEX       ',INDEX)
          CALL LCMGET(IPTRKT,'REMESH      ',REMSH)
          IF(NDIM .EQ. 2) THEN
*----
*  Produce legend
*----
            CALL PSPLEG(IPRINT,ISPSP ,ITYPE ,ICOLR ,NSUR  ,NVOL  ,
     >                  NAMLEG,NUNKNO,FLUX  ,NREGT ,
     >                  MATALB,KEYMRG,KEYFLX,COLRG)
*----
*  Produce graphical description of geometry
*----
            CALL PSPXEL(IPRINT,ISPSP ,ICOLR ,NDIM  ,NSUR  ,NVOL  ,
     >                  NTOTCL,MAXR  ,MINDIM,MAXDIM,KEYMRG,
     >                  INDEX ,REMSH,COLRG)
          ELSE
            WRITE(IOUT,9000)
          ENDIF
          DEALLOCATE(REMSH)
          DEALLOCATE(INDEX,MAXDIM,MINDIM)
        ENDIF
        CALL LCMSIX(IPTRKT,'EXCELL      ',2)
      ENDIF
      DEALLOCATE(MATALB,KEYMRG,COLRG,KEYFLX)
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
 6002 FORMAT('   Processing geometry ',A12)
 9000 FORMAT(' PSP: does not work yet for 3-D',
     >       ' Cartesian geometries')
      END
