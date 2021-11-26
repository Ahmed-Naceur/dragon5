*DECK XELVOL
      SUBROUTINE XELVOL(   IPRT,   NDIM, NEXTGE,   NCPC,  MINDO,  MAXDO,
     >                   ICORDO,  NSURO,  NVOLO, IDLGEO, INDEXO,
     >                     MAXC, REMESH, MATGEO, VOLSO )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute volumes and surfaces.
*
*Copyright:
* Copyright (C) 1987 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* IPRT    intermediate printing level for output.      
* NDIM    number of dimensions (2 or 3).                    
* NEXTGE  rectangular(0)/circular(1) boundary.         
* NCPC    dimension for MINDO.                     
* MINDO   min index values for all axes (rect/cyl).    
* MAXDO   max index values for all axes (rect/cyl).    
* ICORDO  principal axis directions (X/Y/Z) for meshes.
* NSURO   number of surface of the geometry.                
* NVOLO   max. number of track segments in a single track.  
* IDLGEO  relative index of geometry in VOLSO.         
* INDEXO  to retrieve zones in geometry.               
* MAXC    dimension of REMESH.                       
* REMESH  real meshes values (rect/cyl).               
* MATGEO  material numbers used in the geometry.            
*
*Parameters: output
* VOLSO   volumes and surfaces for each geometry.
*
*-----------------------------------------------------------------------
*
      IMPLICIT           NONE
*
      INTEGER            IPRT,NDIM,NEXTGE,NCPC,NSURO,NVOLO,IDLGEO,MAXC
      REAL               VOLSO(*),REMESH(MAXC)
      INTEGER            MINDO(NCPC), MAXDO(NCPC), ICORDO(NCPC),
     >                   INDEXO(4,*), MATGEO(*)
*
      DOUBLE PRECISION   PI, PIO2
      PARAMETER        ( PI=3.14159265358979323846D0,PIO2= 0.5D0*PI)
      INTEGER            IOUT
      PARAMETER        ( IOUT=6 )
      INTEGER            ICUR(4), IND, I, KSUR, JX, JY, JZ, JCY, JRAY,
     >                   IVO, IXYZ, ICT, ICX, ICY, IVO1, IVO2, IP, IR,
     >                   NESUR, NSUNN, NSURC, NSURM,
     >                   NEVOL, NVONN, NVOLC, NVOLM
      INTEGER            NVOL0
      DOUBLE PRECISION   CENTEC(2),RAYONC(2),DXX(2),DYY(2),DZZ,
     >                   XYPOS(2,2),DELTA(3),PRODUC,VOLARE,SURF2,AREAI 
      LOGICAL            LELCRN, SWZCYL
      CHARACTER*4        CORIEN(-6:0)
      SAVE               CORIEN
*
      DATA    CORIEN  
     >       / ' Z+ ',' Z- ',' Y+ ',' Y- ',' X+ ',' X- ','    ' /
*
      IND(I)= IDLGEO + I
*
*     VOL & SURF CALCULATION (CARTESIAN MESHES)
      RAYONC(1)= 0.0D0
      KSUR= MOD(NDIM+1,3)
      DO 50 JX  = MINDO(1)-1, MAXDO(1)
         ICUR(1)= JX
         IF( JX.EQ.MINDO(1)-1 .OR. JX.EQ.MAXDO(1) )THEN
            DELTA(1) = 0.25D0
         ELSE
            DELTA(1)= DBLE(REMESH(JX+1))- DBLE(REMESH(JX))
         ENDIF
      DO 40 JY  = MINDO(2)-1, MAXDO(2)
         ICUR(2)= JY
         IF( JY.EQ.MINDO(2)-1 .OR. JY.EQ.MAXDO(2) )THEN
            DELTA(2) = 0.25D0
         ELSE
            DELTA(2)= DBLE(REMESH(JY+1))-DBLE(REMESH(JY))
         ENDIF
      DO 30 JZ  = MINDO(3)-KSUR, MAXDO(3)+KSUR-1
         ICUR(3)= JZ
         IF( JZ.EQ.MINDO(3)-1 .OR. JZ.EQ.MAXDO(3) )THEN
            DELTA(3) = 0.25D0
         ELSE
            DELTA(3)= DBLE(REMESH(JZ+1))-DBLE(REMESH(JZ))
         ENDIF
         PRODUC= DELTA(1)*DELTA(2)*DELTA(3)
         DO 20 IVO = NSURO, NVOLO
            DO 10 IXYZ= 1, 3
               IF(INDEXO(IXYZ,IND(IVO)).NE.ICUR(IXYZ))GO TO 20
   10       CONTINUE
            VOLSO(IND(IVO))= REAL(PRODUC)
   20    CONTINUE
   30 CONTINUE
   40 CONTINUE
   50 CONTINUE
      NVOL0=0
      DO 130 JCY= 4, NCPC
         ICT= ICORDO(JCY)
         ICX= MOD(ICT  , 3) + 1
         ICY= MOD(ICT+1, 3) + 1
         CENTEC(1)= DBLE(REMESH(MINDO(JCY)-2))
         CENTEC(2)= DBLE(REMESH(MINDO(JCY)-1))
         DO 120 JRAY= MAXDO(JCY), MINDO(JCY), -1
            RAYONC(2)= DBLE(REMESH(JRAY))
         DO 110 JX  = MINDO(ICX), MAXDO(ICX)-1
            ICUR(ICX)= JX
            DXX(1)   = DBLE(REMESH(JX))
            DXX(2)   = DBLE(REMESH(JX+1))
            XYPOS(1,1)   = DXX(1)-CENTEC(1)
            XYPOS(2,1)   = DXX(2)-CENTEC(1)
         DO 100 JY  = MINDO(ICY), MAXDO(ICY)-1
            ICUR(ICY)= JY
            DYY(1)   = DBLE(REMESH(JY))
            DYY(2)   = DBLE(REMESH(JY+1))
            XYPOS(1,2)   = DYY(1)-CENTEC(2)
            XYPOS(2,2)   = DYY(2)-CENTEC(2)
            IF( .NOT.LELCRN(CENTEC,RAYONC,DXX,DYY ))
     >         GO TO 100
            CALL XELCRN(IPRT,RAYONC(2),1,1,XYPOS,AREAI)
            DO 90 JZ = MINDO(ICT)-KSUR, MAXDO(ICT)+KSUR-1
               ICUR(ICT)= JZ
               IF( JZ.EQ.MINDO(ICT)-1 .OR. JZ.EQ.MAXDO(ICT) )THEN
                  DZZ    = 0.25D0
               ELSE
                  DZZ    = DBLE(REMESH(JZ+1)) - DBLE(REMESH(JZ))
               ENDIF
               VOLARE= AREAI * DZZ
               DO 85 IVO1= NSURO, NVOLO
                  ICUR(4)= JRAY-1
                  DO 70 IXYZ= 1, 4
                     IF(INDEXO(IXYZ,IND(IVO1)).NE.ICUR(IXYZ)) GO TO 85
   70             CONTINUE
                  ICUR(4)= JRAY
                  DO 80 IVO2= NSURO, NVOLO
                     DO 75 IXYZ= 1, 4
                      IF(INDEXO(IXYZ,IND(IVO2)).NE.ICUR(IXYZ)) GO TO 80
   75                CONTINUE
                     IF(VOLARE .LE. 0.0D0) THEN
                       IF(NVOL0 .EQ. 0) WRITE(IOUT,8000)
                       NVOL0=NVOL0+1 
                       VOLARE=0.0D0
                     ELSE IF(VOLSO(IND(IVO2))-VOLARE .LE. 0.0) THEN
                       IF(NVOL0 .EQ. 0) WRITE(IOUT,8000)
                       NVOL0=NVOL0+1 
                       VOLARE=DBLE(VOLSO(IND(IVO2)))
                     ENDIF
                     VOLSO(IND(IVO1))= REAL(VOLARE)
                     VOLSO(IND(IVO2))= VOLSO(IND(IVO2)) - REAL(VOLARE)
                     GO TO 85
   80             CONTINUE
   85          CONTINUE
   90       CONTINUE
  100    CONTINUE
  110    CONTINUE
  120    CONTINUE
  130 CONTINUE
      NESUR= 0
      NEVOL= 0
      SWZCYL=.TRUE.
      IF( NEXTGE.EQ.1 )THEN
         DO 700 IVO1= NSURO, NVOLO
           IF( IVO1.LT.0.AND.MATGEO(IND(IVO1)).EQ.0 ) NESUR= NESUR-1
           IF( IVO1.GT.0.AND.MATGEO(IND(IVO1)).LT.0 ) NEVOL= NEVOL+1
 700     CONTINUE
         JRAY = MAXDO(4)
         SURF2= PIO2 * SQRT(DBLE(REMESH(JRAY)))
         ICT= ICORDO(4)
         DO 720 JZ = MINDO(ICT), MAXDO(ICT)-1
            DZZ    = DBLE(REMESH(JZ+1)) - DBLE(REMESH(JZ))
            DO 710 IVO1= 1, NVOLO
              IF( INDEXO(ICT,IND(IVO1)).NE.JZ     ) GO TO 710
              IF( INDEXO( 4,IND(IVO1)).NE.JRAY ) GO TO 710
              SWZCYL= SWZCYL.AND.(MATGEO(IND(IVO1)).LT.0)
              VOLSO(IND(MATGEO(IND(IVO1))))= REAL(DZZ*SURF2)
  710       CONTINUE
  720    CONTINUE
      ENDIF
*
      IF( IPRT.GT.1 )THEN
         NSUNN = NSURO-NESUR-NEVOL
         NSURC = -1
         WRITE(IOUT,'(1H )')
         WRITE(IOUT,'(/35H  CELL SURFACES BEFORE ASSEMBLING         )')
         DO 180 IP  = 1, (9 - NSUNN) / 10
            NSURM= MAX( NSUNN, NSURC-9 )
            WRITE(IOUT,'(10X,10(A5,I7))')
     >           (' SUR ',-IR,IR= NSURC, NSURM, -1)
            WRITE(IOUT,'(8H VALUE  ,2X,1P,10E12.4)')
     >           (4.*VOLSO(IND(IR)),IR=NSURC,NSURM,-1)
            WRITE(IOUT,'(8H ORIENT ,2X,10A12)')
     >           (CORIEN(MATGEO(IND(IR))),IR=NSURC,NSURM,-1)
            WRITE(IOUT,'(1H )')
            NSURC = NSURC - 10
  180    CONTINUE
         NVONN= NVOLO-NEVOL
         NVOLC= 1
         WRITE(IOUT,'(1H )')
         WRITE(IOUT,'( 35H  CELL VOLUMES  BEFORE ASSEMBLING         )')
         DO 190 IP  = 1, (9 + NVONN) / 10
            NVOLM= MIN( NVONN, NVOLC+9 )
            WRITE(IOUT,'(10X,10(A5,I7))')
     >                (' VOL ',IR,IR=NVOLC,NVOLM, 1)
            WRITE(IOUT,'(8H VALUE  ,2X,1P,10E12.4)')
     >               (VOLSO(IND(IR)),IR=NVOLC,NVOLM, 1)
            WRITE(IOUT,'(8H MERGE  ,2X,10I12)')
     >               (MATGEO(IND(IR)),IR=NVOLC,NVOLM, 1)
            WRITE(IOUT,'(1H )')
            NVOLC = NVOLC + 10
  190    CONTINUE
          IF( .NOT.SWZCYL )
     >           CALL XABORT( 'XELVOL: '//
     >                        'CIRCULAR ZONES INCORRECTLY DEFINED' )
      ENDIF
*
      RETURN
*----
*  Formats
*----
 8000 FORMAT(1X,'****** WARNING IN XELVOL ******'/
     >       1X,'AT LEAST ONE REGION WITH NEGATIVE VOLUME'/
     >       1X,'THIS VOLUME IS RESET TO 0.0')
      END
