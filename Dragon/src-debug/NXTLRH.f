*DECK NXTLRH
      SUBROUTINE NXTLRH(IPRINT,IOVER ,ITYPG ,LINMAX,MXSLIN,
     >                  NREG  ,NSUR  ,MESH  ,IDSUR ,IDREG ,
     >                  NBCOR ,NBSINT,ISINT ,TRKLSI,
     >                  IBLIN ,IELIN ,NUMERO,DLENGT,VSI)
*
*----------
*
*Purpose:
* To store line segments in tracking
* vector for HEXT or HEXTZ geometry
* with global region and surface identification.
*
*Copyright:
* Copyright (C) 2012 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): G. Marleau.
*
*Parameters: input
* IPRINT  print level.
* IOVER   flag to everlap tracks on current lines. Here only the case 
*         IOVER=0 is permitted and the tracks are only stored 
*         after the current tracks.
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
*         ISINT(0,I) not used here;
*         IU=ISINT(1,I);
*         IV=ISINT(2,I);
*         IZ=ISINT(3,I);
*         IR=ISINT(4,I);
*         IW=ISINT(5,I).
*         For global HEXT or HEXTZ geometry
*           a) Internal regions
*           (IU,IV,IZ,IW) is location of region.
*           b) External surfaces.
*           One and only one of IX, IY, IZ or IW is negative
*           with:
*              IU=-1 for U- face with position (IV,IZ,IW);
*              IU=-2 for U+ face with position (IV,IZ.IW);
*              IV=-1 for V- face with position (IU,IZ,IW);
*              IV=-2 for V+ face with position (IU,IZ,IW);
*              IZ=-1 for Z- face with position (IX,IY,IW);
*              IZ=-2 for Z+ face with position (IX,IY,IW);
*              IW=-1 for W- face with position (IX,IY,IZ);
*              IW=-2 for W+ face with position (IX,IY,IZ).
* TRKLSI  the surface intersection distance.
* IBLIN   start position for tracking line.
*
*Parameters: input/output
* IELIN   end position for tracking line.
* NUMERO  region/surface identification number
*         for segment.
* DLENGT  spatial location of each line segment.
* VSI     region-surfaces identifier.
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
      INTEGER          VSI(5,-NSUR:NREG)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTLRH')
      DOUBLE PRECISION DCUTOF,DZERO,DONE,DTWO
      PARAMETER       (DCUTOF=1.0D-8,DZERO=0.0D0,DONE=1.0D0,DTWO=2.0D0)
*----
*  Local variables
*----
      INTEGER          IKLIN,IKWLIN,INELIN,IDIR,NXYZ,
     >                 NSBADD,NSEADD,NRGADD
      INTEGER          JSUR,ISUR,KSUR,JREG,IREG,KREG
      DOUBLE PRECISION BEGERR
      DOUBLE PRECISION DELLOC
*----
*  Print header if required
*----
      NXYZ=MESH(1)
      INELIN=0
      ISUR=0
      NSBADD=NBCOR(1)
      NSEADD=NBCOR(2)
      NRGADD=NBSINT+1-NBCOR(2)
*      write(6,*) 'NBSINT,NSBADD,NSEADD,NRGADD',
*     >NBSINT,NSBADD,NSEADD,NRGADD
      IF(IOVER .NE. 0) CALL XABORT(NAMSBR//
     >    ': Invalid option for overlapping geometry')
      IF(IPRINT .GT. 1000) THEN
        WRITE(IOUT,6000) NAMSBR
        WRITE(IOUT,6022) ITYPG
        WRITE(IOUT,6020) 'Initial face'
        DO JSUR=1,NBCOR(1)
          WRITE(IOUT,6021) TRKLSI(JSUR),
     >    (ISINT(IDIR,JSUR),IDIR=1,5)
        ENDDO
        WRITE(IOUT,6020) 'Regions '
        DO JREG=NBCOR(1)+1,NRGADD
          WRITE(IOUT,6021) TRKLSI(JREG),
     >    (ISINT(IDIR,JREG),IDIR=1,5)
        ENDDO
        WRITE(IOUT,6020) 'Final face'
        DO JSUR=NRGADD+1,NBSINT+1
          WRITE(IOUT,6021) TRKLSI(JSUR),
     >    (ISINT(IDIR,JSUR),IDIR=1,5)
        ENDDO
        IF(IELIN .GT. 0) THEN
          WRITE(IOUT,6010)
          WRITE(IOUT,6011) (NUMERO(IKWLIN),IKWLIN=IBLIN-1,IELIN)
          WRITE(IOUT,'(A16)') 'LinePosBefore={ '
          WRITE(IOUT,6012) (DLENGT(IKWLIN),IKWLIN=IBLIN-1,IELIN)
          WRITE(IOUT,'(11X,A2)')'};'
        ENDIF
      ENDIF
*----
*  Initial faces
*----
      IKLIN=IBLIN-1
      NSBADD=NBCOR(1)
      BEGERR=ABS(TRKLSI(1)-DLENGT(IKLIN))
      IF(BEGERR .GT. DCUTOF) THEN
        WRITE(IOUT,9000) NAMSBR,IKLIN,TRKLSI(1),DLENGT(IKLIN)
        CALL XABORT(NAMSBR//
     >  ': Initial tracking position is not nalid')
      ENDIF
*----
*  Initial position
*----
*      JSUR=0
      DO JSUR=1,NSBADD
*        JSUR=JSUR+1
        IF(IPRINT .GT. 200)
     >  WRITE(IOUT,*) 'Initial face id',JSUR,ISINT(1,JSUR),
     >  ISINT(2,JSUR),ISINT(3,JSUR),ISINT(4,JSUR),ISINT(5,JSUR)
        ISUR=0
        DO KSUR=1,NSUR
          IF(ISINT(1,JSUR) .EQ. VSI(1,-KSUR) .AND.
     >       ISINT(2,JSUR) .EQ. VSI(2,-KSUR) .AND.
     >       ISINT(3,JSUR) .EQ. VSI(3,-KSUR) .AND.
     >       ISINT(4,JSUR) .EQ. VSI(4,-KSUR) .AND.
     >       ISINT(5,JSUR) .EQ. VSI(5,-KSUR) ) THEN
            ISUR=KSUR
            GO TO 100 
          ENDIF
        ENDDO
        WRITE(IOUT,*) 'Surface ',JSUR,' not found'
        WRITE(IOUT,'(4I10)') ISINT(1,JSUR),ISINT(2,JSUR),
     >  ISINT(3,JSUR),ISINT(4,JSUR),ISINT(5,JSUR)
        CALL XABORT(NAMSBR//': unable to identify surface')
 100    CONTINUE
*        WRITE(IOUT,*) 'Surface associated with ',JSUR,' is ',ISUR
        IF(IPRINT .GT. 200) THEN
          IF(ISINT(1,JSUR) .EQ. -1) THEN
*----
*  U- faces
*----
            WRITE(IOUT,*) 'U- faces'
          ELSE IF(ISINT(1,JSUR) .EQ. -2) THEN
*----
*  U+ faces
*----
            WRITE(IOUT,*) 'U+ faces'
          ELSE IF(ISINT(2,JSUR) .EQ. -1) THEN
*----
*  V- faces
*----
            WRITE(IOUT,*) 'V- faces'
          ELSE IF(ISINT(2,JSUR) .EQ. -2) THEN
*----
*  V+ faces
*----
            WRITE(IOUT,*) 'V+ faces'
          ELSE IF(ISINT(3,JSUR) .EQ. -1) THEN
*----
*  Z- faces
*----
            WRITE(IOUT,*) 'Z- faces'
          ELSE IF(ISINT(3,JSUR) .EQ. -2) THEN
*----
*  Z+ faces
*----
            WRITE(IOUT,*) 'Z+ faces'
          ELSE IF(ISINT(5,JSUR) .EQ. -1) THEN
*----
*  W- faces
*----
            WRITE(IOUT,*) 'W- faces'
          ELSE IF(ISINT(5,JSUR) .EQ. -2) THEN
*----
*  W+ faces
*----
            WRITE(IOUT,*) 'W+ faces'
          ENDIF
        ENDIF
        IF(IDSUR(ISUR) .NE. 0) THEN
          IKLIN=IKLIN+1
          NUMERO(IKLIN)=-ABS(IDSUR(ISUR))
          DLENGT(IKLIN)=DONE/DBLE(NSBADD)
        ENDIF
      ENDDO
*----
*  Regions
*----
*      JREG=NBCOR(1)
      DELLOC=DZERO
      DO JREG=NBCOR(1)+1,NRGADD
*        JREG=JREG+1
        IREG=0
        DO KREG=1,NREG
          IF(ISINT(1,JREG) .EQ. VSI(1,KREG) .AND.
     >       ISINT(2,JREG) .EQ. VSI(2,KREG) .AND.
     >       ISINT(3,JREG) .EQ. VSI(3,KREG) .AND.
     >       ISINT(4,JREG) .EQ. VSI(4,KREG) .AND.
     >       ISINT(5,JREG) .EQ. VSI(5,KREG) ) THEN
            IREG=KREG
            GO TO 110 
          ENDIF
        ENDDO
        WRITE(IOUT,*) 'Region ',JREG,' not found'
        WRITE(IOUT,'(4I10)') ISINT(1,JREG),ISINT(2,JREG),
     >  ISINT(3,JREG),ISINT(4,JREG),ISINT(5,JREG)
        CALL XABORT(NAMSBR//': unable to identify region')
 110    CONTINUE
*        WRITE(IOUT,*) 'Region associated with ',JREG,' is ',IREG
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
*  Final face
*----
      NSEADD=NBCOR(2)
      INELIN=IKLIN+NSEADD+1
      DLENGT(INELIN)=TRKLSI(NBSINT+1)
      NUMERO(INELIN)=0
*----
*  Final position
*----
*      JSUR=NBSINT+1-NBCOR(2)
      DO JSUR=NRGADD+1,NBSINT+1
*        JSUR=JSUR+1
        IF(IPRINT .GT. 200)
     >  WRITE(IOUT,*) 'Final face id',JSUR,ISINT(1,JSUR),
     >  ISINT(2,JSUR),ISINT(3,JSUR),ISINT(4,JSUR),ISINT(5,JSUR)
        ISUR=0
        DO KSUR=1,NSUR
          IF(ISINT(1,JSUR) .EQ. VSI(1,-KSUR) .AND.
     >       ISINT(2,JSUR) .EQ. VSI(2,-KSUR) .AND.
     >       ISINT(3,JSUR) .EQ. VSI(3,-KSUR) .AND.
     >       ISINT(4,JSUR) .EQ. VSI(4,-KSUR) .AND.
     >       ISINT(5,JSUR) .EQ. VSI(5,-KSUR) ) THEN
            ISUR=KSUR
            GO TO 120 
          ENDIF
        ENDDO
        WRITE(IOUT,*) 'Surface ',JSUR,' not found'
        WRITE(IOUT,'(4I10)') ISINT(1,JSUR),ISINT(2,JSUR),
     >  ISINT(3,JSUR),ISINT(4,JSUR),ISINT(5,JSUR)
        CALL XABORT(NAMSBR//': unable to identify surface')
 120    CONTINUE
*        WRITE(IOUT,*) 'Surface associated with ',JSUR,' is ',ISUR
        IF(IPRINT .GT. 200) THEN
          IF(ISINT(1,JSUR) .EQ. -1) THEN
*----
*  U- faces
*----
            WRITE(IOUT,*) 'U- faces'
          ELSE IF(ISINT(1,JSUR) .EQ. -2) THEN
*----
*  U+ faces
*----
            WRITE(IOUT,*) 'U+ faces'
          ELSE IF(ISINT(2,JSUR) .EQ. -1) THEN
*----
*  V- faces
*----
            WRITE(IOUT,*) 'V- faces'
          ELSE IF(ISINT(2,JSUR) .EQ. -2) THEN
*----
*  V+ faces
*----
            WRITE(IOUT,*) 'V+ faces'
          ELSE IF(ISINT(3,JSUR) .EQ. -1) THEN
*----
*  Z- faces
*----
            WRITE(IOUT,*) 'Z- faces'
          ELSE IF(ISINT(3,JSUR) .EQ. -2) THEN
*----
*  Z+ faces
*----
            WRITE(IOUT,*) 'Z+ faces'
          ELSE IF(ISINT(5,JSUR) .EQ. -1) THEN
*----
*  W- faces
*----
            WRITE(IOUT,*) 'W- faces'
          ELSE IF(ISINT(5,JSUR) .EQ. -2) THEN
*----
*  W+ faces
*----
            WRITE(IOUT,*) 'W+ faces'
          ENDIF
        ENDIF
        IF(IDSUR(ISUR) .NE. 0) THEN
          IKLIN=IKLIN+1
          NUMERO(IKLIN)=-ABS(IDSUR(ISUR))
          DLENGT(IKLIN)=DONE/DBLE(NSEADD)
        ENDIF
      ENDDO
      IKLIN=IKLIN+1
      DLENGT(IKLIN)=TRKLSI(NBSINT+1)
      NUMERO(IKLIN)=0
      IELIN=IKLIN
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
 6022 FORMAT('Geometry type ',I5)
 9000 FORMAT(' **** Error in ',A6,' ****'/
     >       '      Initial face location is invalid :',
     >              I10,1P,2(2X,D20.10))
      END
