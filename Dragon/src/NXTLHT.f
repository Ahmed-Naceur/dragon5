*DECK NXTLHT
      FUNCTION NXTLHT(IPRINT,ITST  ,NDIM  ,MXMESH,LINMAX,
     >                MESH  ,ORITRK,DIRTRK,DCMESH,TRKLIM,
     >                NBCOR ,NBSINT,ISINT ,TRKLSI)
*
*----------
*
*Purpose:
* To track a triangular hexagon in 2-D or 3-D geometry
* using the NXT tracking procedure.
*
*Copyright:
* Copyright (C) 2010 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): G. Marleau.
*
*Parameters: input
* IPRINT  print level.
* ITST    type of tracking, where:
*         =-1   only the exact geometry
*               is considered taking into account the
*               submesh in each direction;
*         = 0   only the global geometry
*               is considered without taking into account the
*               submesh in each direction;
*         = 1   both the global
*               geometry (as a first step) and the exact geometry
*               are considered taking into account the
*               submesh in each direction.
* NDIM    dimension of problem.
* MXMESH  maximum number of spatial subdivision in
*         $U$, $V$, $Z$ and $W$.
* LINMAX  maximum number of segments in a track.
* MESH    effective number of spatial subdivision in
*         each direction ($U$, $W$, $Z$ and $W$).
* ORITRK  a point on the track (origin). The triangular hexagon is
*         assumed centered at the origin.
* DIRTRK  the track direction (director cosines).
* DCMESH  spatial description of the triangular hexagon.
* TRKLIM  beginning and end of track in this cell.
*
*Parameters: output
* NXTLHT  number of side intersections.
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
      INTEGER          IPRINT,ITST,NDIM,MXMESH,LINMAX
      INTEGER          MESH(NDIM)
      DOUBLE PRECISION ORITRK(NDIM),DIRTRK(NDIM),
     >                 DCMESH(-1:MXMESH,5),TRKLIM(2)
      INTEGER          NBCOR(2),NBSINT
      INTEGER          ISINT(0:5,LINMAX)
      DOUBLE PRECISION TRKLSI(LINMAX)
      INTEGER          NXTLHT
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTLHT')
      DOUBLE PRECISION DCUTOF,DZERO,DONE
      PARAMETER       (DCUTOF=1.0D-9,DZERO=0.0D0,DONE=1.0D0)
*----
*  Local variables
*----
      INTEGER          ITC,IEL,KDIR,IDIR,IFRST,ILAST,INEXT(5),IFACE,
     >                 ISF(5),ISL(5),KEL,JEL,IDIRS,MESHL(5)
      DOUBLE PRECISION DEP,SLP,TRKDIS,DISF,DISL
      DOUBLE PRECISION XXX,YYY
*----
*  Data
*----
      CHARACTER        CDIR(5)*1
      DOUBLE PRECISION VNORMD(3,5)
      SAVE             CDIR,VNORMD
      DATA             CDIR /'U','V','Z','R','W'/
      DATA             VNORMD /
     >           1.0D0,0.0D0               ,0.0D0,
     >           0.5D0,0.866025403784439D0 ,0.0D0,
     >           0.0D0,0.0D0               ,1.0D0,
     >           1.0D0,0.0D0               ,0.0D0,
     >           0.5D0,-0.866025403784439D0,0.0D0/
      
*----
*  Verify ITST option and reset to default value if invalid
*----
      IF(ITST .LT. -1 .OR. ITST .GT. 1) THEN
*----
*  Reset ITST=1 (complete analysis) if the value of ITST is invalid
*----
        ITST=1
      ENDIF
*----
*  Initialise output vectors
*----
      NBCOR(1)=0
      NBCOR(2)=0
      MESHL(1)=MESH(1)
      MESHL(2)=MESH(1)
      MESHL(3)=MESH(3)
      MESHL(5)=MESH(1)
      DISF=0.0D0
      DISL=0.0D0
      CALL XDISET(ISINT ,6*LINMAX,0)
      CALL XDDSET(TRKLSI,LINMAX,DZERO)
*----
*  U,V and W MESH identical
*----
      DO ITC=0,MESH(1)
        DCMESH(ITC,2)=DCMESH(ITC,1)
        DCMESH(ITC,5)=DCMESH(ITC,1)
      ENDDO
*----
*  Print header if required
*----
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6000) NAMSBR
        WRITE(IOUT,6011) 'mesh'//CDIR(1)//'={             '
        WRITE(IOUT,6012) (DCMESH(ITC,1),ITC=0,MESHL(1))
        WRITE(IOUT,6013)
        WRITE(IOUT,6011) 'mesh'//CDIR(2)//'={             '
        WRITE(IOUT,6012) (DCMESH(ITC,2),ITC=0,MESHL(2))
        WRITE(IOUT,6013)
        WRITE(IOUT,6011) 'mesh'//CDIR(5)//'={             '
        WRITE(IOUT,6012) (DCMESH(ITC,5),ITC=0,MESHL(5))
        WRITE(IOUT,6013)
        IF(NDIM .EQ. 3) THEN
          WRITE(IOUT,6011) 'mesh'//CDIR(3)//'={             '
          WRITE(IOUT,6012) (DCMESH(ITC,3),ITC=0,MESHL(3))
          WRITE(IOUT,6013)
        ENDIF
        WRITE(IOUT,6016) 'Htrackorigin={      '
        WRITE(IOUT,6012) ORITRK
        WRITE(IOUT,6013)
        WRITE(IOUT,6016) 'Htrackdirection={   '
        WRITE(IOUT,6012) DIRTRK
        WRITE(IOUT,6013)
      ENDIF
*      write(6,*) 'TRKLIM = ',TRKLIM(1),TRKLIM(2)
*----
*  Scan over directions
*----
      IEL=0
      DO IDIR=1,5
        ISF(IDIR)=0
        ISL(IDIR)=0
        IF(IDIR .EQ. 4) GO TO 100
        DEP=DZERO
        SLP=DZERO
        DO ITC=1,3
          DEP=DEP+ORITRK(ITC)*VNORMD(ITC,IDIR)
          SLP=SLP+DIRTRK(ITC)*VNORMD(ITC,IDIR)
        ENDDO
        IF(IDIR .EQ. 3 .AND. NDIM .NE. 3) GO TO 100
*----
*  Select planes order in direction IDIR (forward or backward)
*----
        IF(SLP .LT. 0) THEN
          INEXT(IDIR)=-1
          IFRST=MESHL(IDIR)+1
          ILAST=1
        ELSE
          INEXT(IDIR)=1
          IFRST=1
          ILAST=MESHL(IDIR)+1
        ENDIF
*----
*  Scan over planes in direction IDIR
*----
*        write(6,'(A10,5X,4I10,2F20.10)')' Direction=',IDIR,IFRST,ILAST,
*     >  INEXT(IDIR),DEP,SLP
        DO IFACE=IFRST,ILAST,INEXT(IDIR) 
*----
*  Compute track length required to reach a face
*----
          IF(ABS(SLP) .LT. DCUTOF) CALL XABORT(NAMSBR//
     >    ': line parallel to face not yet programmed')
          TRKDIS=(DCMESH(IFACE-1,IDIR)-DEP)/SLP
          XXX=ABS(TRKLIM(1)-TRKDIS)          
          YYY=ABS(TRKLIM(2)-TRKDIS)          
*          write(6,'(A10,I10,4F20.10)')
*     >    'Face      ',IFACE,DCMESH(IFACE-1,IDIR),TRKDIS,XXX,YYY
*----
*  Store point only if it is inside cell
*----
          IF(XXX .LT. DCUTOF .OR. YYY .LT. DCUTOF) THEN
*            write(6,*) 'Intersection at ',TRKDIS,' for outer face'
            KEL=IEL
            DO JEL=1,IEL
              IF(TRKDIS .LT. TRKLSI(JEL)) THEN
                DO KEL=IEL,JEL,-1
                  TRKLSI(KEL+1)=TRKLSI(KEL)
                  DO ITC=0,5
                    ISINT(ITC,KEL+1)=ISINT(ITC,KEL)
                    ISINT(ITC,KEL)=0
                  ENDDO
                ENDDO
                KEL=JEL-1
                GO TO 111
              ENDIF
            ENDDO
 111        CONTINUE
            TRKLSI(KEL+1)=TRKDIS
            ISINT(0,KEL+1)=INEXT(IDIR)*IDIR
            ISINT(IDIR,KEL+1)=IFACE
*            write(6,*) KEL+1,ISINT(0,KEL+1),ISINT(IDIR,KEL+1)
            IEL=IEL+1
          ELSE IF(TRKDIS .LT. TRKLIM(1) ) THEN
*----
*GMC  Corriger pour direction de la ligne
*----
            IF(ISF(IDIR) .EQ. 0) THEN
              ISF(IDIR)=IFACE+MIN(INEXT(IDIR),0)
              DISF=TRKDIS
            ELSE IF(TRKDIS .GT. DISF) THEN
              ISF(IDIR)=IFACE+MIN(INEXT(IDIR),0)
              DISF=TRKDIS
            ENDIF
*            write(6,*) 'Before 1 ',ISF(IDIR),DISF
          ELSE IF(TRKDIS .GT. TRKLIM(2) ) THEN
*----
*GMC  Corriger pour direction de la ligne
*----
            IF(ISL(IDIR) .EQ. 0) THEN
              ISL(IDIR)=IFACE-MAX(INEXT(IDIR),0)
              DISL=TRKDIS
            ELSE IF(TRKDIS .LT. DISL) THEN
              ISL(IDIR)=IFACE-MAX(INEXT(IDIR),0)
              DISL=TRKDIS
            ENDIF
*            write(6,*) 'After 2 ',ISL(IDIR),DISF
          ELSE
*----
*  Store point in adequate order in TRKLIS
*----
*            write(6,*) 'Intersection at ',TRKDIS,' inside cell'
            KEL=IEL
            DO JEL=1,IEL
              IF(TRKDIS .LT. TRKLSI(JEL)) THEN
                DO KEL=IEL,JEL,-1
                  TRKLSI(KEL+1)=TRKLSI(KEL)
                  DO ITC=0,5
                    ISINT(ITC,KEL+1)=ISINT(ITC,KEL)
                    ISINT(ITC,KEL)=0
                  ENDDO
                ENDDO
                KEL=JEL-1
                GO TO 110
              ENDIF
            ENDDO
 110        CONTINUE
            TRKLSI(KEL+1)=TRKDIS
            ISINT(0,KEL+1)=INEXT(IDIR)*IDIR
            ISINT(IDIR,KEL+1)=IFACE
*            write(6,*) KEL+1,ISINT(0,KEL+1),ISINT(IDIR,KEL+1)
            IEL=IEL+1
          ENDIF
        ENDDO
 100    CONTINUE
      ENDDO
      NBSINT=IEL
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6016) 'HtrackIntA={        '
        DO JEL=1,NBSINT-1
          WRITE(IOUT,6014) TRKLSI(JEL),
     >    (ISINT(IDIR,JEL),IDIR=0,5)
        ENDDO
        JEL=NBSINT
        WRITE(IOUT,6015) TRKLSI(JEL),
     >  (ISINT(IDIR,JEL),IDIR=0,5)
      ENDIF
*----
*  Process corners inside cell
*----
      KEL=1
      DO JEL=2,NBSINT
        XXX=ABS(TRKLSI(JEL)-TRKLSI(KEL))
        IF(XXX .LT. DCUTOF) THEN
          ISINT(0,KEL)=5
          IDIR=ABS(ISINT(0,JEL))
          ISINT(IDIR,KEL)=ISINT(IDIR,JEL)
        ELSE
          KEL=KEL+1
          TRKLSI(KEL)=TRKLSI(JEL)
          DO IDIR=0,5
            ISINT(IDIR,KEL)=ISINT(IDIR,JEL)
          ENDDO
        ENDIF  
      ENDDO
      NBSINT=KEL
      IF(IPRINT .GE. 100) THEN
        WRITE(IOUT,6016) 'HtrackIntB={        '
        DO JEL=1,NBSINT-1
          WRITE(IOUT,6014) TRKLSI(JEL),
     >    (ISINT(IDIR,JEL),IDIR=0,5)
        ENDDO
        JEL=NBSINT
        WRITE(IOUT,6015) TRKLSI(JEL),
     >  (ISINT(IDIR,JEL),IDIR=0,5)
      ENDIF
*----
*  Identify final faces
*----
      JEL=NBSINT
      IDIR=ABS(ISINT(0,JEL))
      ISINT(0,JEL+1)=ISINT(0,JEL)
      TRKLSI(JEL+1)=TRKLSI(JEL)
      NBCOR(2)=1
      DO IDIRS=1,5
        IF(IDIRS .EQ. IDIR) THEN
          IF(ISINT(IDIRS,JEL) .EQ. 1) THEN
            ISINT(IDIRS,JEL+1)=-1
          ELSE
            ISINT(IDIRS,JEL+1)=-2
          ENDIF
        ELSE
          ISINT(IDIRS,JEL+1)=ISL(IDIRS)
        ENDIF
      ENDDO
*----
*  Regions
*----
      DO JEL=NBSINT,2,-1
        TRKLSI(JEL)=TRKLSI(JEL)-TRKLSI(JEL-1)
        IDIR=ABS(ISINT(0,JEL))
        IF(IDIR .EQ. 5) THEN
          DO IDIRS=1,5
            KDIR=ISINT(IDIRS,JEL)
            IF(KDIR .NE. 0) THEN
              ISINT(IDIRS,JEL)=KDIR-MAX(INEXT(IDIRS),0)
            ELSE
              ISINT(IDIRS,JEL)=ISINT(IDIRS,JEL+1)
            ENDIF
          ENDDO
        ELSE
          DO IDIRS=1,5
            IF(IDIRS .EQ. IDIR) THEN
              ISINT(IDIRS,JEL)=ISINT(IDIRS,JEL)-MAX(INEXT(IDIR),0)
            ELSE
              ISINT(IDIRS,JEL)=ISINT(IDIRS,JEL+1)
            ENDIF
          ENDDO
        ENDIF
      ENDDO
*----
*  Identify initial faces
*----
      JEL=1
      IDIR=ABS(ISINT(0,JEL))
      NBCOR(1)=1
      DO IDIRS=1,5
        IF(IDIRS .EQ. IDIR) THEN
          IF(ISINT(IDIRS,JEL) .EQ. 1) THEN
            ISINT(IDIRS,JEL)=-1
          ELSE
            ISINT(IDIRS,JEL)=-2
          ENDIF
        ELSE
          ISINT(IDIRS,JEL)=ISF(IDIRS)
        ENDIF
      ENDDO
      NBSINT=NBSINT
*----
*  Print final track information
*----
      IF(IPRINT .GT. 100) THEN
        WRITE(IOUT,6011) 'Initial face        '
        JEL=1
        WRITE(IOUT,6010) TRKLSI(JEL),
     >  (ISINT(IDIR,JEL),IDIR=1,5)
        WRITE(IOUT,6011) 'Regions             '
        DO JEL=2,NBSINT
          WRITE(IOUT,6010) TRKLSI(JEL),
     >    (ISINT(IDIR,JEL),IDIR=1,5)
        ENDDO
        WRITE(IOUT,6011) 'Final face          '
        JEL=NBSINT+1
        WRITE(IOUT,6010) TRKLSI(JEL),
     >  (ISINT(IDIR,JEL),IDIR=1,5)
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      NXTLHT=NBSINT
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6010 FORMAT(1X,F25.16,6I10)
 6011 FORMAT(A20)
 6012 FORMAT(6(1X,F25.16,:,','))
 6013 FORMAT('};')
 6014 FORMAT(('{',F25.16,6(',',I5),'},'))
 6015 FORMAT(('{',F25.16,9(',',I5),'}};'))
 6016 FORMAT(A20)
      END
