*DECK NXTLCU
      SUBROUTINE NXTLCU(IPRINT,MXSLIN,ICOMB ,
     >                  NBCOR ,NBSINT,ISINT ,TRKLSI)
*
*----------
*
*Purpose:
* To merge two sets of tracks and store the result in tracking
* vector with adequate region identification.
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
* MXSLIN  maximum number of segments in a subgeometry track.
* ICOMB   flag for combination if two sets of tracks with:
*         ICOMB=-1  one considers only the first
*         track set (Cartesian); ICOMB=-2 one considers
*         only the second track set (annular); ICOMB=1, 
*         the Cartesian and annular regions are super imposed while
*         preserving the outer Cartesian boundary;
*         ICOMB=2 the Cartesian and annular regions
*         are super imposed while
*         preserving the outer annular boundary.
*
*Parameters: input/output
* NBCOR   number of corner found for each external faces.
* NBSINT  number of surface crossed by track.
* ISINT   direction of plane intersected and
*         the surfaces crossed by the track.
* TRKLSI  the surface intersection distance.
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
      INTEGER          IPRINT,MXSLIN,ICOMB,NBCOR(2,3),NBSINT(3)
      INTEGER          ISINT(0:5,MXSLIN,3)
      DOUBLE PRECISION TRKLSI(MXSLIN,3)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTLCU')
      DOUBLE PRECISION DCUTOF,DZERO,DONE,DTWO
      PARAMETER       (DCUTOF=1.0D-7,DZERO=0.0D0,DONE=1.0D0,DTWO=2.0D0)
*----
*  Local variables
*----
      INTEGER          IT,IC,IA,IG,ISURC,JSURC,ISURA,ISURG,JSURG,
     >                 IFIC,ILIC,ITIC,IFIA,ILIA,ITIA,ISUR,JSUR,ITB,ITE
      INTEGER          IDIR,KDIR,JDIR,LDIR,IRADG,IFSG
      DOUBLE PRECISION REFLOC,REFLOA,REFLOG,CURLOC,CURLOA,
     >                 THICK,DELCUR
*----
*  Data
*----
      CHARACTER        NAMTYP(3)*12
      SAVE             NAMTYP
      DATA             NAMTYP
     >               /'Cartesian   ','Annular     ','Combined    '/
*----
*  Print header if required
*----
      IF(IPRINT .GT. 1000) THEN
        WRITE(IOUT,6000) NAMSBR
        ITB=1
        ITE=2
        IF(ICOMB .EQ. -1) ITE=1
        IF(ICOMB .EQ. -2) ITB=2
        DO IT=ITB,ITE
          WRITE(IOUT,6016) NAMTYP(IT)
          WRITE(IOUT,6011) 'Initial face'
          DO JSUR=1,NBCOR(1,IT)
            WRITE(IOUT,6010) TRKLSI(JSUR,IT),
     >      (ISINT(IDIR,JSUR,IT),IDIR=1,5)
          ENDDO
          WRITE(IOUT,6011) 'Regions '
          DO JSUR=NBCOR(1,IT)+1,NBSINT(IT)-NBCOR(2,IT)+1
            WRITE(IOUT,6010) TRKLSI(JSUR,IT),
     >      (ISINT(IDIR,JSUR,IT),IDIR=1,5)
          ENDDO
          WRITE(IOUT,6011) 'Final face'
          DO JSUR=NBSINT(IT)-NBCOR(2,IT)+2,NBSINT(IT)+1
            WRITE(IOUT,6010) TRKLSI(JSUR,IT),
     >      (ISINT(IDIR,JSUR,IT),IDIR=1,5)
          ENDDO
        ENDDO
      ENDIF
      IFIA=1
      ILIA=1
      LDIR=1
      JSURC=0
      ISURG=0
*----
*  Create combined tracking
*----
      IC=1
      IA=2
      IG=3
      NBCOR(1,IG)=0
      NBCOR(2,IG)=0
      IF(ABS(ICOMB) .EQ. 1) THEN
*----
*  Save startup Cartesian face
*----
        ISURG=0
        IFIC=1
        ILIC=NBCOR(1,IC)
        DO ISURC=IFIC,ILIC
          ISURG=ISURG+1
          TRKLSI(ISURG,IG)=TRKLSI(ISURC,IC)
          DO IDIR=1,5
            ISINT(IDIR,ISURG,IG)=ISINT(IDIR,ISURC,IC)
          ENDDO
        ENDDO
        NBCOR(1,IG)=ISURG
        ISURC=NBCOR(1,IC)
        REFLOC=TRKLSI(ISURC,IC)
        ISURA=NBCOR(1,2)
        REFLOA=TRKLSI(ISURA,IA)
        KDIR=0
        IRADG=0
        IFIC=NBCOR(1,IC)+1
        ILIC=NBSINT(IC)-NBCOR(2,IC)+1
        REFLOG=REFLOC
        DELCUR=ABS(REFLOC-REFLOA)
        IF(ICOMB .EQ. 1) THEN
*----
*  Find radial surface direction
*----
          DO IDIR=1,5
            IF(ISINT(IDIR,ISURA,IA) .LT. 0) THEN
              KDIR=IDIR
            ENDIF
          ENDDO
          IF(KDIR .EQ. 4) THEN
*----
*  It is an inner face,
*  Scan Cartesian regions
*----
            ITIC=0
            DO ISURC=IFIC,ILIC
              CURLOC=REFLOC+TRKLSI(ISURC,IC)
              ISURG=ISURG+1
              DO JDIR=1,5
                ISINT(JDIR,ISURG,IG)=ISINT(JDIR,ISURC,IC)
              ENDDO
              TRKLSI(ISURG,IG)=CURLOC-REFLOG
              ISINT(4,ISURG,IG)=IRADG
              IF(ICOMB .EQ. 1) THEN
                IF(REFLOA .LE. CURLOC) THEN
                  IRADG=ISINT(4,ISURA+1,IA)
                  THICK=REFLOA-REFLOG
                  TRKLSI(ISURG,IG)=THICK
                  REFLOG=REFLOA
                  GO TO 205
                ENDIF
              ENDIF
              ITIC=ITIC+1
              REFLOC=CURLOC
              REFLOG=REFLOC
            ENDDO
 205        CONTINUE
            IFIC=IFIC+ITIC
          ELSE IF(DELCUR .LT. DCUTOF) THEN
*----
*  It is an outer face,
*  reset radial id for Cartesian face
*----
            REFLOG=REFLOC
            IRADG=ISINT(4,1,IA)
            DO ISURC=1,NBCOR(1,IC)
              IF(ISINT(KDIR,ISURC,IC) .LT. 0) GO TO 105
            ENDDO
            CALL XABORT (NAMSBR//
     >      ': Cartesian and annular face not coherent')
 105        CONTINUE
            DO JSURG=1,NBCOR(1,IG)
              ISINT(4,JSURG,IG)=IRADG
            ENDDO
          ELSE
            CALL XABORT(NAMSBR//
     >      ': Initial Cartesian and annular faces are incoherent')
          ENDIF
          IFIA=NBCOR(1,IA)+1
          ILIA=NBSINT(IA)-NBCOR(2,IA)+1
        ELSE
          REFLOG=REFLOC
        ENDIF
*----
*  Scan Cartesian regions
*----
        DO ISURC=IFIC,ILIC
          CURLOC=REFLOC+TRKLSI(ISURC,IC)
*          IF(IPRINT .GT. 1000)
*     >    write(6,7000) 'ISURC=',ISURC,CURLOC,
*     >    (ISINT(JDIR,ISURC,IC),JDIR=1,5)
          ISURG=ISURG+1
          DO JDIR=1,5
            ISINT(JDIR,ISURG,IG)=ISINT(JDIR,ISURC,IC)
          ENDDO
          TRKLSI(ISURG,IG)=CURLOC-REFLOG
          ISINT(4,ISURG,IG)=IRADG
          IF(ICOMB .EQ. 1) THEN
            ITIA=0
*----
*  Scan annular regions for intermediate mesh
*----
            DO ISURA=IFIA,ILIA
              CURLOA=REFLOA+TRKLSI(ISURA,IA)
              DELCUR=ABS(CURLOC-CURLOA)
*              IF(IPRINT .GT. 1000)
*     >        write(6,7000) 'ISURA=',ISURA,CURLOA,
*     >        (ISINT(JDIR,ISURA,IA),JDIR=0,5)
              IF(DELCUR .LT. DCUTOF) THEN
                ITIA=ITIA+1
                IF(ISURA .EQ. ILIA) THEN
                  IRADG=0
                ELSE
                  IRADG=ISINT(4,ISURA+1,IA)
                ENDIF
*                IRADG=ISINT(4,ISURA,IA)
                DO JDIR=1,3
                  IF(ISINT(JDIR,ISURA,IA) .NE. 0) THEN
                    ISINT(JDIR,ISURG,IG)=ISINT(JDIR,ISURA,IA)
                  ENDIF
                ENDDO
                IF(ISINT(5,ISURA,IA) .NE. 0) THEN
                  ISINT(5,ISURG,IG)=ISINT(5,ISURA,IA)
                ENDIF
*                IF(IPRINT .GT. 1000)
*     >          write(6,7001) 'CURLOA .EQ. CURLOC',ISURG,
*     >          (ISINT(JDIR,ISURG,IG),JDIR=1,5)
                THICK=CURLOA-REFLOG
                REFLOA=CURLOA
                REFLOG=REFLOA
                GO TO 115
              ELSE IF(CURLOA .GT. CURLOC) THEN
*                IF(IPRINT .GT. 1000)
*     >          write(6,7001) 'CURLOA .GT. CURLOC=',ISURG,
*     >          (ISINT(JDIR,ISURG,IG),JDIR=1,5)
                GO TO 115
              ELSE
                ITIA=ITIA+1
                IF(ISURA .EQ. ILIA) THEN
                  IRADG=0
                ELSE
                  IRADG=ISINT(4,ISURA+1,IA)
                ENDIF
                DO JDIR=1,3
                  ISINT(JDIR,ISURG+1,IG)=ISINT(JDIR,ISURG,IG)
                ENDDO
                ISINT(5,ISURG+1,IG)=ISINT(5,ISURG,IG)
              ENDIF
              THICK=CURLOA-REFLOG
              ISINT(4,ISURG+1,IG)=IRADG
              TRKLSI(ISURG+1,IG)=TRKLSI(ISURG,IG)-THICK
              TRKLSI(ISURG,IG)=THICK
              ISURG=ISURG+1
              REFLOA=CURLOA
              REFLOG=REFLOA
*          IF(IPRINT .GT. 1000)
*     >    write(6,7001) 'CURLOA .LT. CURLOC',ISURG,
*     >    (ISINT(JDIR,ISURG,IG),JDIR=1,5)
            ENDDO
 115        CONTINUE
            IFIA=IFIA+ITIA
          ENDIF
          REFLOC=CURLOC
          REFLOG=REFLOC
        ENDDO
*----
*  Save final Cartesian faces
*----
        IFIC=NBSINT(IC)-NBCOR(2,IC)+2
        ILIC=NBSINT(IC)+1
        IFSG=ISURG
        DO ISURC=IFIC,ILIC
          ISURG=ISURG+1
          TRKLSI(ISURG,IG)=TRKLSI(ISURC,IC)
          DO IDIR=1,5
            ISINT(IDIR,ISURG,IG)=ISINT(IDIR,ISURC,IC)
          ENDDO
          ISINT(4,ISURG,IG)=IRADG
        ENDDO
        NBCOR(2,IG)=NBCOR(2,IC)
*----
*  Find radial zone associated with final Cartesian face
*  when cylinder face coincides with final face
*----
        IF(ICOMB .EQ. 1) THEN
*----
*  Find radial surface direction
*----
          ISURC=NBSINT(IC)+1
          REFLOC=TRKLSI(ISURC,IC)
          ISURA=NBSINT(IA)+1
          REFLOA=TRKLSI(ISURA,IA)
          DELCUR=ABS(REFLOC-REFLOA)
          DO IDIR=1,5
            IF(ISINT(IDIR,ISURA,IA) .LT. 0) THEN
              KDIR=IDIR
            ENDIF
          ENDDO
          IF(KDIR .NE. 4) THEN
*----
*  It is an inner face,
*  Scan Cartesian regions
*----
            IF(DELCUR .LT. DCUTOF) THEN
*----
*  It is an outer face,
*  reset radial id for Cartesian face
*----
              REFLOG=REFLOC
              IRADG=ISINT(4,ISURA,IA)
              JSURG=IFSG
              DO ISURC=IFIC,ILIC
                JSURG=JSURG+1
                ISINT(4,JSURG,IG)=IRADG
              ENDDO
            ELSE
              WRITE(IOUT,9002) ISURA,KDIR,REFLOC,REFLOA
              CALL XABORT(NAMSBR//
     >        ': Final Cartesian and annular faces are incoherent')
            ENDIF
          ENDIF
          IFIA=NBCOR(1,IA)+1
          ILIA=NBSINT(IA)-NBCOR(2,IA)+1
        ELSE
          REFLOG=REFLOC
        ENDIF
      ELSE IF(ABS(ICOMB) .EQ. 2) THEN
*----
*  Find radial surface direction
*----
        ISURG=0
        ISURA=NBCOR(1,IA)
        REFLOA=TRKLSI(ISURA,IA)
        REFLOG=REFLOA
        IF(ICOMB .EQ. 2) THEN
          DO IDIR=1,5
            IF(ISINT(IDIR,ISURA,IA) .LT. 0) THEN
              KDIR=IDIR
            ELSE IF(ISINT(IDIR,ISURA,IA) .GT. 0) THEN
              LDIR=IDIR
            ENDIF
          ENDDO
          IF(KDIR .EQ. 4) THEN
*----
*  It is an outer radial face,
*  Determine Cartesian region location for
*  this point
*----
            ISURC=NBCOR(1,IC)
            REFLOC=TRKLSI(ISURC,IC)
            IF(REFLOC .GT. REFLOA) THEN
              WRITE(IOUT,9000) ISURA,ISURC,REFLOA,REFLOC
              CALL XABORT(NAMSBR//
     >': No Cartesian region found for this annulus')
            ENDIF
            IFIC=NBCOR(1,IC)+1
            ILIC=NBSINT(IC)-NBCOR(2,IC)+1
            DO JSURC=IFIC,ILIC
              CURLOC=REFLOC+TRKLSI(JSURC,IC)
              IF(REFLOA .LT. CURLOC) THEN
                ISURC=JSURC
                GO TO 125
              ENDIF
              REFLOC=CURLOC
            ENDDO
            CALL XABORT(NAMSBR//
     >': Impossible to find Cartesian region containing annulus')
 125        CONTINUE
            IFIC=ISURC
          ELSE
*----
*  It is an outer Cartesian face,
*  Determine Cartesian region location for
*  this face
*----
            ISURC=NBCOR(1,IC)
            REFLOC=TRKLSI(ISURC,IC)
            IF(REFLOC .NE. REFLOA) CALL XABORT(NAMSBR//
     >': No compatible Cartesian face found for this annulus')
            IFIC=ISURC+1
          ENDIF
          IFIA=1
          ILIA=NBCOR(1,IA)
          DO ISURA=IFIA,ILIA
            ISURG=ISURG+1
            TRKLSI(ISURG,IG)=TRKLSI(ISURA,IA)
            DO IDIR=1,5
              IF(IDIR .EQ. KDIR .OR. IDIR .EQ. LDIR ) THEN
                ISINT(IDIR,ISURG,IG)=ISINT(IDIR,ISURA,IA)
              ELSE
                ISINT(IDIR,ISURG,IG)=ISINT(IDIR,ISURC,IC)
              ENDIF
            ENDDO
            ISINT(4,ISURG,IG)=ISINT(4,ISURA,IA)
          ENDDO
          NBCOR(1,IG)=ISURG
          IFIA=NBCOR(1,IA)+1
          IRADG=ISINT(4,IFIA,IA)
          ILIA=NBSINT(IA)-NBCOR(2,IA)+1
*----
*  Scan Cartesian regions
*----
          ILIC=NBSINT(IC)-NBCOR(2,IC)+1
          DO ISURC=IFIC,ILIC
            CURLOC=REFLOC+TRKLSI(ISURC,IC)
            ISURG=ISURG+1
            DO JDIR=1,5
              ISINT(JDIR,ISURG,IG)=ISINT(JDIR,ISURC,IC)
            ENDDO
            ISINT(4,ISURG,IG)=IRADG
            THICK=CURLOC-REFLOG
            TRKLSI(ISURG,IG)=THICK
            ITIA=0
            JSURC=ISURC
*----
*  Scan annular regions for intermediate mesh
*----
            DO ISURA=IFIA,ILIA
              CURLOA=REFLOA+TRKLSI(ISURA,IA)
              IF(CURLOA .GT. CURLOC) THEN
                GO TO 145
              ELSE IF(CURLOA .EQ. CURLOC) THEN
                ITIA=ITIA+1
                IRADG=ISINT(4,ISURA,IA)
                DO JDIR=1,5
                  IF(ISINT(JDIR,ISURA,IA) .NE. 0) THEN
                    ISINT(JDIR,ISURG,IG)=ISINT(JDIR,ISURA,IA)
                  ENDIF
                ENDDO
                ISINT(4,ISURG,IG)=IRADG
*                IF(IPRINT .GT. 1000)
*     >          write(6,7001) 'CURLOA .EQ. CURLOC',ISURG,
*     >          (ISINT(JDIR,ISURG,IG),JDIR=1,5)
                THICK=CURLOA-REFLOG
                REFLOA=CURLOA
                REFLOG=REFLOA
                GO TO 145
              ELSE
                IF(ISURA .EQ. ILIA) THEN
                  THICK=CURLOA-REFLOG
                  TRKLSI(ISURG,IG)=THICK
                  IRADG=0
                  GO TO 155
                ELSE
                  ITIA=ITIA+1
                  IRADG=ISINT(4,ISURA+1,IA)
                  THICK=CURLOA-REFLOG
                  DO JDIR=1,5
                    ISINT(JDIR,ISURG+1,IG)=ISINT(JDIR,ISURG,IG)
                  ENDDO
                  ISINT(4,ISURG+1,IG)=IRADG
                  TRKLSI(ISURG+1,IG)=TRKLSI(ISURG,IG)-THICK
                  TRKLSI(ISURG,IG)=THICK
                  ISURG=ISURG+1
                ENDIF
                REFLOA=CURLOA
                REFLOG=REFLOA
              ENDIF
            ENDDO
 145        CONTINUE
            IFIA=IFIA+ITIA
            REFLOC=CURLOC
            REFLOG=REFLOC
          ENDDO
 155      CONTINUE
*----
*  Save final annular faces
*----
          IFIA=NBSINT(IA)-NBCOR(2,IA)+2
          ILIA=NBSINT(IA)+1
          ISURA=IFIA
          DO IDIR=1,5
            IF(ISINT(IDIR,ISURA,IA) .LT. 0) THEN
              KDIR=IDIR
            ELSE IF(ISINT(IDIR,ISURA,IA) .GT. 0) THEN
              LDIR=IDIR
            ENDIF
          ENDDO
          ISURC=JSURC
          IFIA=NBSINT(IA)-NBCOR(2,IA)+2
          ILIA=NBSINT(IA)+1
          DO ISURA=IFIA,ILIA
            ISURG=ISURG+1
            TRKLSI(ISURG,IG)=TRKLSI(ISURA,IA)
             DO IDIR=1,5
              IF(IDIR .EQ. KDIR .OR. IDIR .EQ. LDIR ) THEN
                ISINT(IDIR,ISURG,IG)=ISINT(IDIR,ISURA,IA)
              ELSE
                ISINT(IDIR,ISURG,IG)=ISINT(IDIR,ISURC,IC)
              ENDIF
            ENDDO
            ISINT(4,ISURG,IG)=ISINT(4,ISURA,IA)
          ENDDO
          NBCOR(2,IG)=NBCOR(2,IA)
        ELSE
          IFIA=1
          ILIA=NBCOR(1,IA)
          DO ISURA=IFIA,ILIA
            ISURG=ISURG+1
            TRKLSI(ISURG,IG)=TRKLSI(ISURA,IA)
            DO IDIR=1,5
              ISINT(IDIR,ISURG,IG)=ISINT(IDIR,ISURA,IA)
            ENDDO
          ENDDO
          NBCOR(1,IG)=NBCOR(1,IA)
*----
*  Scan annular regions for intermediate mesh
*----
          IFIA=NBCOR(1,IA)+1
          ILIA=NBSINT(IA)-NBCOR(2,IA)+1
          DO ISURA=IFIA,ILIA
            ISURG=ISURG+1
            CURLOA=REFLOA+TRKLSI(ISURA,IA)
            DO JDIR=1,5
              ISINT(JDIR,ISURG,IG)=ISINT(JDIR,ISURA,IA)
            ENDDO
            TRKLSI(ISURG,IG)=TRKLSI(ISURA,IA)
            REFLOA=CURLOA
          ENDDO
          IFIA=NBSINT(IA)-NBCOR(2,IA)+2
          ILIA=NBSINT(IA)+1
          DO ISURA=IFIA,ILIA
            ISURG=ISURG+1
            TRKLSI(ISURG,IG)=TRKLSI(ISURA,IA)
            DO IDIR=1,5
              ISINT(IDIR,ISURG,IG)=ISINT(IDIR,ISURA,IA)
            ENDDO
          ENDDO
          NBCOR(2,IG)=NBCOR(2,IA)
          NBSINT(IG)=ISURG-1
        ENDIF
      ENDIF
*----
*  Remove track segment with vanishing distance
*----
      IT=3
      NBSINT(IT)=ISURG-1
      ISUR=NBCOR(1,IT)
      DO JSUR=NBCOR(1,IT)+1,NBSINT(IT)-NBCOR(2,IT)+1
        IF(TRKLSI(JSUR,IT) .GT. DCUTOF) THEN
          ISUR=ISUR+1
          TRKLSI(ISUR,IT)=TRKLSI(JSUR,IT)
          DO IDIR=1,5
            ISINT(IDIR,ISUR,IT)=ISINT(IDIR,JSUR,IT)
          ENDDO
        ELSE
          IF(IPRINT .GT. 1000)
     >    WRITE(IOUT,9001) JSUR,IT,ISUR
        ENDIF
      ENDDO
      IF(ISUR .LT. NBSINT(IT)-NBCOR(2,IT)+1) THEN
        DO JSUR=NBSINT(IT)-NBCOR(2,IT)+2,NBSINT(IT)+1
          ISUR=ISUR+1
          TRKLSI(ISUR,IT)=TRKLSI(JSUR,IT)
          DO IDIR=1,5
            ISINT(IDIR,ISUR,IT)=ISINT(IDIR,JSUR,IT)
          ENDDO
        ENDDO
        NBSINT(IT)=ISUR-1
      ENDIF
      IF(IPRINT .GT. 1000) THEN
        WRITE(IOUT,6016) NAMTYP(IT)
        WRITE(IOUT,6011) 'Initial face'
        DO JSUR=1,NBCOR(1,IT)
          WRITE(IOUT,6010) TRKLSI(JSUR,IT),
     >    (ISINT(IDIR,JSUR,IT),IDIR=1,5)
        ENDDO
        WRITE(IOUT,6011) 'Regions '
        DO JSUR=NBCOR(1,IT)+1,NBSINT(IT)-NBCOR(2,IT)+1
          WRITE(IOUT,6010) TRKLSI(JSUR,IT),
     >    (ISINT(IDIR,JSUR,IT),IDIR=1,5)
        ENDDO
        WRITE(IOUT,6011) 'Final face'
        DO JSUR=NBSINT(IT)-NBCOR(2,IT)+2,NBSINT(IT)+1
          WRITE(IOUT,6010) TRKLSI(JSUR,IT),
     >    (ISINT(IDIR,JSUR,IT),IDIR=1,5)
        ENDDO
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6010 FORMAT(1X,F25.16,5I10)
 6011 FORMAT(A20)
 6016 FORMAT(1X,A12)
* 7000 FORMAT('Verify ',A6,I10,F25.16,6I10)
* 7001 FORMAT('Verify ',A18,I10,6I10)
 9000 FORMAT(' Problem with surface ',2I10,2F20.10)
 9001 FORMAT(' Warning : region with vanishing distance found ',
     >       3I10)
 9002 FORMAT('Final Cartesian and annular faces are incoherent',
     >       2(2X,I10),2(2X,F20.10))
      END
