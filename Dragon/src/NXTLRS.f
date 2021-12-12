*DECK NXTLRS
      SUBROUTINE NXTLRS(IPRINT,IOVER ,ITYPG ,LINMAX,MXSLIN,
     >                  NREG  ,NSUR  ,MESH  ,IDSUR ,IDREG ,
     >                  NBCOR ,NBSINT,ISINT ,TRKLSI,
     >                  IBLIN ,IELIN ,NUMERO,DLENGT)
*
*----------
*
*Purpose:
* To store line segments in tracking
* vector with global region and surface identification.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): G. Marleau.
*
*Parameters: input
* IPRINT  print level.
* IOVER   flag to everlap tracks on current lines. In the case where
*         IOVER=0, the tracks are only stored after the current
*         tracks. When IOVER=1, the new tracks segments are
*         overlapped over the existing tracks.
* ITYPG   type of geometry.
* LINMAX  maximum number of segments in a complete track.
* MXSLIN  maximum number of segments in a subgeometry track.
* NREG    maximum number of regions in geometry.
* NSUR    maximum number of surfaces in geometry.
* MESH    effective number of spatial subdivision in $X$
*         $Y$, $Z$ and $R$.
* IDSUR   local surface identifier.
* IDREG   local region identifier.
* NBCOR   number of corner found for each external faces.
* NBSINT  number of surface crossed by track.
* ISINT   direction of plane intersected and
*         the surfaces crossed by the track:
*         ISINT(0,I) not used here
*         IX=ISINT(1,I)
*         IY=ISINT(2,I)
*         IZ=ISINT(3,I)
*         IR=ISINT(4,I)
*         IW=ISINT(5,I) not used.
*         For global Cartesian geometry
*           a) Internal regions
*           (IX,IY,IZ,IR) is location of region.
*           b) External surfaces.
*           One and only one of IX, IY or IZ is negative
*           with:
*              IX=-1 for X- face with position (IY,IZ,IR);
*              IX=-2 for X+ face with position (IY,IZ.IR);
*              IY=-1 for Y- face with position (IX,IZ,IR);
*              IY=-2 for Y+ face with position (IX,IZ,IR);
*              IZ=-1 for Z- face with position (IX,IY,IR);
*              IZ=-2 for Z+ face with position (IX,IY,IR).
*         For global annular geometry
*           a) Internal regions
*           (IX,IY,IZ,IR) is location of region.
*           b) External surfaces.
*           One and only one of IX, IY, IZ or IR is negative
*           with:
*              IX=-1 for X- face with position (IY,IZ,IR);
*              IX=-2 for X+ face with position (IY,IZ,IR);
*              IY=-1 for Y- face with position (IX,IZ,IR);
*              IY=-2 for Y+ face with position (IX,IZ,IR);
*              IZ=-1 for Z- face with position (IX,IY,IR);
*              IZ=-2 for Z+ face with position (IX,IY,IR);
*              IR=-2 for R+ face with position (IX,IY,IZ).
*         For hexagonal geometry
*           a) Internal regions
*           (IH,IZ,IR) is location of region
*           IH=IX=IY=1 permitted only.
*           b) External surfaces.
*           Pair IH=IX=IY < 0 or IZ< 10
*           with
*              IH=-1 for hexagonal face 1 in plane IZ;
*              IH=-2 for hexagonal face 2 in plane IZ;
*              IH=-3 for hexagonal face 3 in plane IZ;
*              IH=-4 for hexagonal face 4 in plane IZ;
*              IH=-5 for hexagonal face 5 in plane IZ;
*              IH=-6 for hexagonal face 6 in plane IZ.
* TRKLSI  the surface intersection distance.
* IBLIN   start position for tracking line.
*
*Parameters: input/output
* IELIN   end position for tracking line.
* NUMERO  region/surface identification number
*         for segment.
* DLENGT  spatial location of each line segment.
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
      INTEGER          IPRINT,IOVER,ITYPG,LINMAX,MXSLIN,NREG,NSUR
      INTEGER          MESH(4),IDSUR(NSUR),IDREG(NREG),
     >                 NBCOR(2),NBSINT,ISINT(0:5,MXSLIN)
      DOUBLE PRECISION TRKLSI(MXSLIN)
      INTEGER          IBLIN,IELIN,NUMERO(LINMAX)
      DOUBLE PRECISION DLENGT(LINMAX)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTLRS')
      DOUBLE PRECISION DCUTOF,DZERO,DONE,DTWO
      PARAMETER       (DCUTOF=1.0D-8,DZERO=0.0D0,DONE=1.0D0,DTWO=2.0D0)
*----
*  Local variables
*----
      INTEGER          NXO,NYO,NZO,NRO,IKLIN,IKWLIN,INELIN,
     >                 IOFXM,IOFXP,IOFYM,IOFYP,IOFZM,IOFZP,IOFRP,ISEG,
     >                 IXR,IYR,IZR,IRR,IDIR,
     >                 NROX,NROY,NROZ,NDIM
      INTEGER          NSBADD,NSEADD,NRGADD,NDISP,IKB,IKE,IPB,IPE
      INTEGER          JSUR,ISUR,KSUR,JREG,IREG,KREG
      DOUBLE PRECISION BEGO1,BEGO2,ENDO1,ENDO2,BEGNEW,ENDNEW,BEGERR
      DOUBLE PRECISION DELLOC
      INTEGER          ILINI
      INTEGER          IHEX,IOFHP
*----
*  Print header if required
*----
      INELIN=0
      IOFXP=0
      IOFYM=0
      IOFYP=0
      ISUR=0
      NSBADD=NBCOR(1)
      NRGADD=NBSINT+1-NBCOR(1)-NBCOR(2)
      NSEADD=NBCOR(2)
      IF(IPRINT .GT. 1000) THEN
        WRITE(IOUT,6000) NAMSBR
        WRITE(IOUT,6020) 'Initial face'
        DO JSUR=1,NBCOR(1)
          WRITE(IOUT,6021) TRKLSI(JSUR),
     >    (ISINT(IDIR,JSUR),IDIR=1,5)
        ENDDO
        WRITE(IOUT,6020) 'Regions '
        DO JSUR=NBCOR(1)+1,NBSINT-NBCOR(2)+1
          WRITE(IOUT,6021) TRKLSI(JSUR),
     >    (ISINT(IDIR,JSUR),IDIR=1,5)
        ENDDO
        WRITE(IOUT,6020) 'Final face'
        DO JSUR=NBSINT-NBCOR(2)+2,NBSINT+1
          WRITE(IOUT,6021) TRKLSI(JSUR),
     >    (ISINT(IDIR,JSUR),IDIR=1,5)
        ENDDO
        WRITE(IOUT,6010)
        WRITE(IOUT,6011) (NUMERO(IKWLIN),IKWLIN=IBLIN-1,IELIN)
        WRITE(IOUT,'(A16)') 'LinePosBefore={ '
        WRITE(IOUT,6012) (DLENGT(IKWLIN),IKWLIN=IBLIN-1,IELIN)
        WRITE(IOUT,'(11X,A2)')'};'
      ENDIF
      IKB=0
      IPB=-1
      IKE=0
      IPE=0
      IF(IOVER .EQ. 1) THEN
*----
*  Find initial and final track locations where line
*  starts to overlapp
*----
        INELIN=IELIN
        BEGNEW=TRKLSI(1)
        BEGO1=DLENGT(IBLIN-1)
        DELLOC=ABS(BEGNEW-BEGO1)
        IF(IPRINT .GT. 2000)
     >  WRITE(IOUT,*) ' Begin ->',BEGNEW,BEGO1,DELLOC
        IKB=0
*        IF(BEGNEW .GT. BEGO1) THEN
        IF(DELLOC .GT. DCUTOF) THEN
          DO ISEG=IBLIN,IELIN-1
            IKB=IKB+1
*        IF(IPRINT .GT. 2000)
*     >  WRITE(IOUT,*) ' + Segment ',ISEG,' ->',IKB,
*     >           NUMERO(ISEG),DLENGT(ISEG)
            IF(NUMERO(ISEG) .GT. 0) THEN
              BEGO2=BEGO1+DLENGT(ISEG)
              DELLOC=ABS(BEGNEW-BEGO2)
*        IF(IPRINT .GT. 2000)
*     >  WRITE(IOUT,*) BEGO1,BEGO2,DELLOC
              IF(DELLOC .LT. DCUTOF) THEN
                BEGO1=BEGNEW-BEGO1
                BEGO2=DZERO
                IPB=ISEG
                GO TO 105
              ELSE IF(BEGNEW .LE. BEGO2) THEN
                BEGO1=BEGNEW-BEGO1
                BEGO2=BEGO2-BEGNEW
                IPB=ISEG
                GO TO 105
              ENDIF
              BEGO1=BEGO2
            ENDIF
          ENDDO
*----
*  Abort: Overlapping for first points is impossible
*----
          BEGO1=DLENGT(IBLIN-1)
          WRITE(IOUT,*) ' Begin ->',BEGNEW,BEGO1
          IKB=0
          DO ISEG=IBLIN,IELIN-1
            IKB=IKB+1
            WRITE(IOUT,*) ' + Segment ',ISEG,' ->',IKB,
     >           NUMERO(ISEG),DLENGT(ISEG)
            IF(NUMERO(ISEG) .GT. 0) THEN
              BEGO2=BEGO1+DLENGT(ISEG)
              WRITE(IOUT,*) BEGO1,BEGO2
              BEGO1=BEGO2
            ENDIF
          ENDDO
          CALL XABORT(NAMSBR//': Overlapping for first point '//
     >    ' is impossible')
 105      CONTINUE
          NSBADD=0
        ELSE
          DO ISEG=IBLIN,IELIN
            IF(NUMERO(ISEG) .GT. 0) GO TO 115
            NSBADD=NSBADD-1
          ENDDO
 115      CONTINUE
          NSBADD=MAX(NSBADD,0)
        ENDIF
        ENDNEW=TRKLSI(NBSINT+1)
        ENDO2=DLENGT(IELIN)
        DELLOC=ABS(ENDNEW-ENDO2)
        IF(IPRINT .GT. 2000)
     >  WRITE(IOUT,*) ' Ends ->',ENDNEW,ENDO2,DELLOC
        IKE=0
*        IF(ENDNEW .LT. ENDO2) THEN
        IF(DELLOC .GT. DCUTOF) THEN
          DO ISEG=IELIN,IBLIN,-1
            IKE=IKE+1
*        IF(IPRINT .GT. 2000)
*     >  WRITE(IOUT,*) ' - Segment ',ISEG,' ->',IKE
            IF(NUMERO(ISEG) .GT. 0) THEN
              ENDO1=ENDO2-DLENGT(ISEG)
              DELLOC=ABS(ENDNEW-ENDO1)
*        IF(IPRINT .GT. 2000)
*     >  WRITE(IOUT,*) DLENGT(ISEG),ENDO1,ENDO2,DELLOC
              IF(DELLOC .LT. DCUTOF) THEN
                ENDO1=DZERO
                ENDO2=ENDO2-ENDNEW
                IPB=ISEG
                GO TO 125
              ELSE IF(ENDNEW .GE. ENDO1) THEN
                ENDO1=ENDNEW-ENDO1
                ENDO2=ENDO2-ENDNEW
                IPB=ISEG
                GO TO 125
              ENDIF
              ENDO2=ENDO1
            ENDIF
          ENDDO
*----
*  Abort: Overlapping for last points is impossible
*----
          ENDO2=DLENGT(IELIN)
          WRITE(IOUT,*) ' Ends ->',ENDNEW,ENDO2
          IKE=0
          DO ISEG=IELIN,IBLIN,-1
            IKE=IKE+1
            WRITE(IOUT,*) ' - Segment ',ISEG,' ->',IKE
            IF(NUMERO(ISEG) .GT. 0) THEN
              ENDO1=ENDO2-DLENGT(ISEG)
              WRITE(IOUT,*) DLENGT(ISEG),ENDO1,ENDO2
              ENDO2=ENDO1
            ENDIF
          ENDDO
          CALL XABORT(NAMSBR//': Overlapping for last point '//
     >    ' is impossible')
 125      CONTINUE
        ENDIF
*        IF(IPRINT .GT. 2000)
*     >  WRITE(IOUT,*) 'IKB,IKE =',IKB,IKE
        IF(IKE .EQ. 0) THEN
*----
*  Last surface is overlapped
*  Next segment position avalable is IBLIN+IKB
*----
          IF(IKB .GT. 0) THEN
*----
*  Cut line segment at IBLIN+IKB-1 if region intersection
*----
            ISEG=IBLIN+IKB-1
            DLENGT(ISEG)=BEGO1
          ENDIF
        ELSE
*----
*  Displace segments starting at IELIN-IKE+1
*  towards end, the displacement is of NRGADD+NSBADD
*----
          NDISP=NRGADD+NSBADD+1
          IF(IPB .EQ. IPE) NDISP=NDISP+1
          DO ISEG=IELIN,IELIN-IKE+1,-1
*        IF(IPRINT .GT. 2000)
*     >  WRITE(IOUT,*) 'Displace =',ISEG,ISEG+NDISP
            DLENGT(ISEG+NDISP)=DLENGT(ISEG)
            NUMERO(ISEG+NDISP)=NUMERO(ISEG)
            DLENGT(ISEG)=DZERO
          ENDDO
          ILINI=IELIN+1
          IF(IKB .GT. 0) THEN
            ILINI=IBLIN+IKB
          ELSE
            ILINI=IBLIN+1
          ENDIF
          DO ISEG=IELIN-IKE+NDISP,ILINI,-1
*        IF(IPRINT .GT. 2000)
*     >  WRITE(IOUT,*) 'Re-initialize =',ISEG
            DLENGT(ISEG)=DZERO
          ENDDO
          INELIN=INELIN+NDISP
*----
*  Cut line segment at IELIN-IKE+1 for region intersection
*----
          ISEG=IELIN-IKE+1+NDISP
          DLENGT(ISEG)=ENDO2
*        IF(IPRINT .GT. 2000)
*     >  WRITE(IOUT,*) 'End cut at =',ISEG,ENDO2
          IF(IKB .GT. 0) THEN
*----
*  Cut line segment at IBLIN+IKB-1 if region intersection
*----
            ISEG=IBLIN+IKB-1
            DLENGT(ISEG)=BEGO1
*        IF(IPRINT .GT. 2000)
*     >  WRITE(IOUT,*) 'start cut at =',ISEG,BEGO1
          ENDIF
        ENDIF
      ENDIF
      NDIM=3
      IF(ITYPG .EQ.  3 .OR. ITYPG .EQ.  5 .OR.
     >   ITYPG .EQ.  8 .OR. ITYPG .EQ. 20 .OR.
     >   ITYPG .EQ. 24) NDIM=2
      NXO=MESH(1)
      NYO=MESH(2)
      NZO=MESH(3)
      IF(NDIM .EQ. 2) NZO=1
      NRO=MESH(4)
      NROX=1
      NROY=1
      NROZ=1
      IHEX=0
      IF(ITYPG .GE. 20) NRO=NRO+1
      IF(ITYPG .EQ. 6) THEN
        NROX=0
        NROY=0
        NROZ=NRO
      ELSE IF(ITYPG .EQ. 8) THEN
        NROY=0
        NROZ=NRO
        IHEX=1
      ELSE IF(ITYPG .EQ. 9) THEN
        NROY=0
        IHEX=1
      ELSE IF(ITYPG .EQ. 10) THEN
        NROX=NRO
        NROY=0
        NROZ=0
      ELSE IF(ITYPG .EQ. 11) THEN
        NROX=0
        NROY=NRO
        NROZ=0
      ELSE IF(ITYPG .EQ. 21) THEN
        NROX=NRO
      ELSE IF(ITYPG .EQ. 22) THEN
        NROY=NRO
      ELSE IF(ITYPG .EQ. 23) THEN
        NROZ=NRO
      ELSE IF(ITYPG .EQ. 24) THEN
        NROY=0
        NROZ=NRO
        IHEX=1
      ELSE IF(ITYPG .EQ. 25) THEN
        NROY=0
        IHEX=1
        NROZ=NRO
      ENDIF
      IF(IHEX .EQ. 0) THEN
        IOFXM=0
        IOFXP=IOFXM+NYO*NZO*NROX
        IOFYM=IOFXP+NYO*NZO*NROX
        IOFYP=IOFYM+NXO*NZO*NROY
        IOFZM=IOFYP+NXO*NZO*NROY
        IOFZP=IOFZM+NXO*NYO*NROZ
        IOFRP=IOFZP+NXO*NYO*NROZ
      ELSE
        IOFHP=0
        IOFZM=IOFHP+6*NZO
        IOFZP=IOFZM+NROZ
        IOFRP=IOFZP+NROZ
      ENDIF
*      IF(IPRINT .GT. 2000)
*     >WRITE(IOUT,*) 'IOF=',NXO,NYO,NZO,NRO,NROX,NROY,NROZ,
*     >IOFXM,IOFXP,IOFYM,IOFYP,IOFZM,IOFZP,IOFRP
*----
*  Initial faces
*----
      IF(IKB .EQ. 0) THEN
        IKLIN=IBLIN+IKB-1
        NSBADD=NBCOR(1)
        BEGERR=ABS(TRKLSI(1)-DLENGT(IKLIN))
        IF(BEGERR .GT. DCUTOF) THEN
          WRITE(IOUT,9000) NAMSBR,IKLIN,TRKLSI(1),DLENGT(IKLIN)
          CALL XABORT(NAMSBR//
     >    ': Initial tracking position is not valid')
        ENDIF
      ELSE
        IKLIN=IBLIN+IKB-1
        NSBADD=0
      ENDIF
*----
*  Initial position
*----
      JSUR=0
*      IF(IPRINT .GT. 2000)
*     >WRITE(IOUT,*) IOFXM,IOFXP,IOFYM,IOFYP,IOFZM,IOFZP
*      IF(IPRINT .GT. 2000)
*     >WRITE(IOUT,*) NRO,NXO,NYO,NZO
      DO KSUR=1,NSBADD
        JSUR=JSUR+1
*        IF(IPRINT .GT. 2000)
*     >  WRITE(IOUT,*) 'Initial face id',JSUR,ISINT(1,JSUR),
*     >             ISINT(2,JSUR),ISINT(3,JSUR),ISINT(4,JSUR)
        IXR=MAX(1,ISINT(1,JSUR))
        IYR=MAX(1,ISINT(2,JSUR))
        IZR=MAX(1,ISINT(3,JSUR))
        IF(ISINT(1,JSUR) .LT. 0 .AND. ISINT(2,JSUR) .LT. 0) THEN
*----
*  H+ faces
*----
*          IF(IPRINT .GT. 2000)
*     >    WRITE(IOUT,*) 'H- faces'
          IF(ISINT(1,JSUR) .NE. ISINT(2,JSUR)) CALL XABORT(NAMSBR//
     >': X and Y face id must be identical for hexagons')
          ISUR=6*(IZR-1)-ISINT(1,JSUR)+IOFHP
        ELSE IF(ISINT(1,JSUR) .EQ. -1) THEN
*----
*  X- faces
*----
*          IF(IPRINT .GT. 2000)
*     >    WRITE(IOUT,*) 'X- faces'
          IRR=NROX
          IF(ISINT(4,JSUR) .GT. 0) IRR=ISINT(4,JSUR)
          ISUR=IRR+((IZR-1)*NYO+IYR-1)*NROX+IOFXM
        ELSE IF(ISINT(1,JSUR) .EQ. -2) THEN
*----
*  X+ faces
*----
*          IF(IPRINT .GT. 2000)
*     >    WRITE(IOUT,*) 'X+ faces'
          IRR=NROX
          IF(ISINT(4,JSUR) .GT. 0) IRR=ISINT(4,JSUR)
          ISUR=IRR+((IZR-1)*NYO+IYR-1)*NROX+IOFXP
        ELSE IF(ISINT(2,JSUR) .EQ. -1) THEN
*----
*  Y- faces
*----
*          IF(IPRINT .GT. 2000)
*     >    WRITE(IOUT,*) 'Y- faces',IZR,IXR,((IXR-1)*NZO+IZR-1)
          IRR=NROY
          IF(ISINT(4,JSUR) .GT. 0) IRR=ISINT(4,JSUR)
          ISUR=IRR+((IXR-1)*NZO+IZR-1)*NROY+IOFYM
        ELSE IF(ISINT(2,JSUR) .EQ. -2) THEN
*----
*  Y+ faces
*----
*          IF(IPRINT .GT. 2000)
*     >    WRITE(IOUT,*) 'Y+ faces',IZR,IXR,((IXR-1)*NZO+IZR-1)
          IRR=NROY
          IF(ISINT(4,JSUR) .GT. 0) IRR=ISINT(4,JSUR)
          ISUR=IRR+((IXR-1)*NZO+IZR-1)*NROY+IOFYP
        ELSE IF(ISINT(3,JSUR) .EQ. -1) THEN
*----
*  Z- faces
*----
*          IF(IPRINT .GT. 2000)
*     >    WRITE(IOUT,*) 'Z- faces'
          IRR=NROZ
          IF(ISINT(4,JSUR) .GT. 0) IRR=ISINT(4,JSUR)
          ISUR=IRR+((IYR-1)*NXO+IXR-1)*NROZ+IOFZM
        ELSE IF(ISINT(3,JSUR) .EQ. -2) THEN
*----
*  Z+ faces
*----
*          IF(IPRINT .GT. 2000)
*     >    WRITE(IOUT,*) 'Z+ faces'
          IRR=NROZ
          IF(ISINT(4,JSUR) .GT. 0) IRR=ISINT(4,JSUR)
          ISUR=IRR+((IYR-1)*NXO+IXR-1)*NROZ+IOFZP
        ELSE IF(ISINT(4,JSUR) .EQ. -2) THEN
*----
*  R+ faces
*----
*          IF(IPRINT .GT. 2000)
*     >    WRITE(IOUT,*) 'R+ faces'
          IRR=NRO
          IF(ISINT(4,JSUR) .GT. 0) IRR=ISINT(4,JSUR)
          ISUR=((IZR-1)*NYO+IYR-1)*NXO+IXR+IOFRP
        ENDIF
*        IF(IPRINT .GT. 2000)
*     >  WRITE(IOUT,*) 'Surface  ',ISUR,' =',
*     >  IXR,IYR,IZR,IRR,IDSUR(ISUR),IKLIN+1
        IF(IDSUR(ISUR) .NE. 0) THEN
          IKLIN=IKLIN+1
          NUMERO(IKLIN)=-ABS(IDSUR(ISUR))
          DLENGT(IKLIN)=DONE/DBLE(NSBADD)
        ENDIF
      ENDDO
*----
*  Regions
*----
      JREG=NBCOR(1)
      DELLOC=DZERO
      DO KREG=1,NRGADD
        JREG=JREG+1
        IXR=MAX(1,ISINT(1,JREG))
        IYR=MAX(1,ISINT(2,JREG))
        IZR=MAX(1,ISINT(3,JREG))
        IRR=NRO
        IF(ISINT(4,JREG) .GT. 0) IRR=ISINT(4,JREG)
        IF(NRO .GT. 0) THEN
          IREG=IRR+((IXR-1)+((IYR-1)+(IZR-1)*NYO)*NXO)*NRO
        ELSE
          IREG=IXR+((IYR-1)+(IZR-1)*NYO)*NXO
        ENDIF
        IF(IPRINT .GT. 2000)
     >  WRITE(IOUT,*) 'Region  ',IREG,' =',
     >  IXR,IYR,IZR,IRR,IREG,JREG,IDREG(IREG),TRKLSI(JREG),IKLIN+1
        IF(IDREG(IREG) .NE. 0) THEN
          IKLIN=IKLIN+1
          NUMERO(IKLIN)=ABS(IDREG(IREG))
          DLENGT(IKLIN)=TRKLSI(JREG)+DELLOC
          DELLOC=DZERO
        ELSE
*----
*  Add distance to DELLOC. This will be added to the next
*  track with non zero region.
*----
          DELLOC=DELLOC+TRKLSI(JREG)
        ENDIF
      ENDDO
*----
*  Final faces
*----
      IF(IKE .EQ. 0) THEN
        NSEADD=NBCOR(2)
        INELIN=IKLIN+NSEADD+1
        DLENGT(INELIN)=TRKLSI(NBSINT+1)
        NUMERO(INELIN)=0
      ELSE
        NSEADD=0
      ENDIF
      JSUR=NBSINT+1-NBCOR(2)
      DO KSUR=1,NSEADD
        JSUR=JSUR+1
*        IF(IPRINT .GT. 2000)
*     >  WRITE(IOUT,*) 'Final face id',JSUR,ISINT(1,JSUR),
*     >             ISINT(2,JSUR),ISINT(3,JSUR),ISINT(4,JSUR)
        IXR=MAX(1,ISINT(1,JSUR))
        IYR=MAX(1,ISINT(2,JSUR))
        IZR=MAX(1,ISINT(3,JSUR))
        IF(ISINT(1,JSUR) .LT. 0 .AND. ISINT(2,JSUR) .LT. 0) THEN
*----
*  H+ faces
*----
*          IF(IPRINT .GT. 2000)
*     >    WRITE(IOUT,*) 'H- faces'
          IF(ISINT(1,JSUR) .NE. ISINT(2,JSUR)) CALL XABORT(NAMSBR//
     >': X and Y face id must be identical for hexagons')
          ISUR=6*(IZR-1)-ISINT(1,JSUR)+IOFHP
        ELSE IF(ISINT(1,JSUR) .EQ. -1) THEN
*----
*  X- faces
*----
*          IF(IPRINT .GT. 2000)
*     >    WRITE(IOUT,*) 'X- faces'
          IRR=NROX
          IF(ISINT(4,JSUR) .GT. 0) IRR=ISINT(4,JSUR)
          ISUR=IRR+((IZR-1)*NYO+IYR-1)*NROX+IOFXM
        ELSE IF(ISINT(1,JSUR) .EQ. -2) THEN
*----
*  X+ faces
*----
*          IF(IPRINT .GT. 2000)
*     >    WRITE(IOUT,*) 'X+ faces'
          IRR=NROX
          IF(ISINT(4,JSUR) .GT. 0) IRR=ISINT(4,JSUR)
          ISUR=IRR+((IZR-1)*NYO+IYR-1)*NROX+IOFXP
        ELSE IF(ISINT(2,JSUR) .EQ. -1) THEN
*----
*  Y- faces
*----
*          IF(IPRINT .GT. 2000)
*     >    WRITE(IOUT,*) 'Y- faces'
          IRR=NROY
          IF(ISINT(4,JSUR) .GT. 0) IRR=ISINT(4,JSUR)
          ISUR=IRR+((IXR-1)*NZO+IZR-1)*NROY+IOFYM
        ELSE IF(ISINT(2,JSUR) .EQ. -2) THEN
*----
*  Y+ faces
*----
*          IF(IPRINT .GT. 2000)
*     >    WRITE(IOUT,*) 'Y+ faces'
          IRR=NROY
          IF(ISINT(4,JSUR) .GT. 0) IRR=ISINT(4,JSUR)
          ISUR=IRR+((IXR-1)*NZO+IZR-1)*NROY+IOFYP
        ELSE IF(ISINT(3,JSUR) .EQ. -1) THEN
*----
*  Z- faces
*----
*          IF(IPRINT .GT. 2000)
*     >    WRITE(IOUT,*) 'Z- faces'
          IRR=NROZ
          IF(ISINT(4,JSUR) .GT. 0) IRR=ISINT(4,JSUR)
          ISUR=IRR+((IYR-1)*NXO+IXR-1)*NROZ+IOFZM
        ELSE IF(ISINT(3,JSUR) .EQ. -2) THEN
*----
*  Z+ faces
*----
*          IF(IPRINT .GT. 2000)
*     >    WRITE(IOUT,*) 'Z+ faces'
          IRR=NROZ
          IF(ISINT(4,JSUR) .GT. 0) IRR=ISINT(4,JSUR)
          ISUR=IRR+((IYR-1)*NXO+IXR-1)*NROZ+IOFZP
        ELSE IF(ISINT(4,JSUR) .EQ. -2) THEN
*----
*  R+ faces
*----
*          IF(IPRINT .GT. 2000)
*     >    WRITE(IOUT,*) 'R+ faces'
          IRR=NRO
          IF(ISINT(4,JSUR) .GT. 0) IRR=ISINT(4,JSUR)
          ISUR=((IZR-1)*NYO+IYR-1)*NXO+IXR+IOFRP
        ENDIF
*        IF(IPRINT .GT. 2000)
*     >  WRITE(IOUT,*) 'Surface  ',ISUR,' =',
*     >  IXR,IYR,IZR,IRR,IDSUR(ISUR),IKLIN+1
        IKLIN=IKLIN+1
        IF(IDSUR(ISUR) .NE. 0) THEN
          NUMERO(IKLIN)=-ABS(IDSUR(ISUR))
          DLENGT(IKLIN)=DONE/DBLE(NSEADD)
        ELSE
          NUMERO(IKLIN)=0
          DLENGT(IKLIN)=DZERO
        ENDIF
      ENDDO
      IF(IOVER .EQ. 0) THEN
        IKLIN=IKLIN+1
        DLENGT(IKLIN)=TRKLSI(NBSINT+1)
        NUMERO(IKLIN)=0
        IELIN=IKLIN
      ELSE
        IELIN=INELIN
      ENDIF
      IF(IPRINT .GT. 1000) THEN
        WRITE(IOUT,6010)
        WRITE(IOUT,6011) (NUMERO(IKWLIN),IKWLIN=IBLIN-1,IELIN)
        WRITE(IOUT,'(A10)') 'LinePos={ '
        WRITE(IOUT,6012) (DLENGT(IKWLIN),IKWLIN=IBLIN-1,IELIN)
        WRITE(IOUT,'(11X,A2)')'};'
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6010 FORMAT(' Integration line :')
 6011 FORMAT(6(I25,:,','))
 6012 FORMAT(6(F25.16,:,','))
 6020 FORMAT(A20)
 6021 FORMAT(1X,F25.16,6I10)
 9000 FORMAT(' **** Error in ',A6,' ****'/
     >       '      Initial face location is invalid :',
     >              I10,1P,2(2X,D20.10))
      END
