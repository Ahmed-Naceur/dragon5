*DECK XELTI3
      SUBROUTINE XELTI3( IPRT,IFTEMP,NANGLE,DENUSR,ISYMM,ANGLES,DENSTY,
     >                   NTOTCL,NEXTGE,MAXR,REMESH,LINMAX,RCUTOF,
     >                   NSUR,NVOL,INDEL,MINDIM,
     >                   MAXDIM,ICOORD,INCR,ICUR,TRKBEG,CONV,TRKDIR,
     >                   LENGHT,NUMERO,DDENWT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Construct the sequential tape that will contain tracks for  
* isotropic BC in 3-D calculation.
*
*Copyright:
* Copyright (C) 1989 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* IPRT    intermediate printing level for output.      
* IFTEMP  tracking file number.                             
* NANGLE  number of angles used in the tracking process.    
* DENUSR  density of tracks in the plane perpendicular 
*              to the tracking angles.                 
* ISYMM   flag for symetry (1/0 for on/off):         
*          2 reflection plane normal to X axis;          
*          4 reflection plane normal to Y axis;          
*          8 reflection plane normal to X and Y axis;   
*         16 reflection plane normal to Z axis;          
*         18 reflection plane normal to X and Z axis;    
*         20 reflection plane normal to Y and Z axis;    
*         24 reflection plane normal to X, Y and Z axis. 
* ANGLES  3d angle values.                             
* DENSTY  density of tracks angle by angle.            
* NTOTCL  number of cylindres of a type + 3.                
* NEXTGE  for tubez, nextge=1                          
* MAXR    max number of real mesh values in REMESH.  
* REMESH  real mesh values (rect/cyl).                 
* LINMAX  max. number of track segments in a single track.  
* RCUTOF  cutof for corner tracking( 0.25 suggested )  
* NSUR    number of surfaces.                               
* NVOL    number of zones.                                  
* INDEL   numbering of surfaces & zones.                    
* MINDIM  min index values for all axes (rect/cyl).    
* MAXDIM  max index values for all axes (rect/cyl).    
* ICOORD  principal axes direction (X/Y/Z) for meshes. 
*   ICUR  current zonal location for a track segment.  
*   INCR  increment direction for next track segment.  
* TRKBEG  position where a track begins.               
*   CONV  segments of tracks.                          
* TRKDIR  direction of a track in all axes.            
* LENGHT  relative lenght of each segment in a track.  
* NUMERO  material identification of each track segment.
* DDENWT  density of tracks angle by angle.            
*
*-----------------------------------------------------------------------
*
      IMPLICIT      NONE
*
*     DECLARE       DUMMY ARGUMENTS
      INTEGER       IPRT,IFTEMP,NANGLE,NTOTCL,NEXTGE,MAXR,LINMAX,
     >              NSUR,NVOL
      INTEGER       MINDIM(NTOTCL), MAXDIM(NTOTCL), ICUR(NTOTCL),
     >              ICOORD(NTOTCL), INCR(NTOTCL), NUMERO(LINMAX),
     >              INDEL(4,*)
      REAL          DENUSR,REMESH(MAXR),TRKBEG(NTOTCL),TRKDIR(NTOTCL),
     >              CONV(NTOTCL),RCUTOF
      DOUBLE PRECISION ANGLES(3,*),DENSTY(*),LENGHT(LINMAX),
     >                 DDENWT(NANGLE)
      INTEGER       ISYMM
      INTEGER       IQUART(4)
*
*     DECLARE       LOCAL VARIABLES
      REAL          TRKEND(3), TRKORI(3), TRKOR2(3), PROJC2(3),
     >              ANGEQN(3,3), ANGLE2(3), ANGLE3(3), BARY(3),
     >              TRKCUT(3,2), TCUTOF(3,4), TORIC(3)
      INTEGER       NSBEG(4), NSEND(4), NSCUT(2)
      LOGICAL       LANGLE
      EQUIVALENCE ( ANGEQN(1,2), ANGLE2 ), ( ANGEQN(1,3), ANGLE3 )
      CHARACTER     TEDATA*13
      INTEGER       NPAN, IOUT
      PARAMETER   ( NPAN=3, IOUT=6 )
*
      INTEGER       NDIM,         I, J, NPOINT, NPO2, NCUTOF,
     >              NOTRAK, NSOLAN,      IANG,                   ISB,
     >              NANGLS, IPAN, NESTIM, JANG, IX, IY, IREF1, I2, I3,
     >              NSOLMX, NTTRK, IANG0, NDEBS, N, IZZ, IANGL, K, K3,
     >              NCROS, LINE, JC,     ISE, IDIM
      REAL          ANORM2, A, B, DENLIN, RCIRC, R2CIRC, DP,
     >              SURTOT, VOLTOT, DEPART, X, Y, ANN, ODDNXT,
     >              TOTLEN, TOTXXX, DDENST
      REAL          DZ,ZMAX
      DOUBLE PRECISION WEIGHT, WZ
*
      ANORM2( A, B ) = A*A + B*B
      NDIM= 3
*
*     ONE WEIGHT FOR ALL TRACKS
*     DENLIN= # OF TRACKS / CM
      DENLIN= SQRT(DENUSR)
*
*     COMPUTE THE CIRCUMSCRIBED RADIUS AND
*             THE  COORDINATES FOR THE TRUE CENTER OF THE CELL
      R2CIRC= 0.0
      DO 10 I = 1, 3
         BARY(I)= 0.5 * (REMESH(MAXDIM(I)) + REMESH(MINDIM(I)))
         IF( NEXTGE.EQ.1 )THEN
            CALL XABORT('XELTI3: TUBEZ NOT SUPPORTED')
         ELSE
            R2CIRC =  R2CIRC
     >             + (REMESH(MAXDIM(I)) - REMESH(MINDIM(I)))**2
         ENDIF
   10 CONTINUE
      ZMAX=MAX(ABS(REMESH(MAXDIM(3))),ABS(REMESH(MINDIM(3))))
      R2CIRC= 0.25 * R2CIRC
      RCIRC = SQRT(R2CIRC)
*
*     NPOINT= # OF TRACKS ALONG ONE PERPENDICULAR AXIS
      NPOINT= INT( 2. * RCIRC * DENLIN )
***** BEWARE ***** BEWARE ***** BEWARE ***** BEWARE ***** BEWARE *****
***** CHANGE THIS "NPOINT" PARAMETER HAS TREMENDOUS EFFECTS ON TRACKING
***** BEWARE ***** BEWARE ***** BEWARE ***** BEWARE ***** BEWARE *****
*
*     OTHER POSSIBLE CHOICES (EXPLORED WITHOUT SUCCESS) ARE ==>
*1-)  NPOINT=  INT( 2. * RCIRC * DENLIN ) + 1
*2-)  NPOINT= NINT( 2. * RCIRC * DENLIN )
*3-)  NPOINT= NINT( 2. * RCIRC * DENLIN ) + 1
*
*     KEEP "NPOINT" ODD & CORRECT DENSITY
      NPO2  = NPOINT / 2
      NPOINT= 2 * NPO2 + 1
      DP    = 2. * RCIRC / NPOINT
      DENLIN= 1. / DP
      DENUSR= DENLIN**2
      ODDNXT= (2*NPO2+3) / (2.*RCIRC)
      IF( RCUTOF.EQ.0.0 )THEN
         NCUTOF= 1
      ELSE
         NCUTOF= 4
      ENDIF
      NOTRAK= 0
      SURTOT= 0.0
      VOLTOT= 0.0
      IQUART(1)=1
      IQUART(2)=1
      IQUART(3)=1
      IQUART(4)=1
      IF(NEXTGE .EQ. 1) THEN
        IQUART(2)=0
        IQUART(3)=0
        IQUART(4)=0
        NSOLAN= (NANGLE * (NANGLE + 2)) / 2
      ELSE
        IF( ISYMM .EQ. 8 .OR. ISYMM .EQ. 24 )THEN
          NSOLAN= (NANGLE * (NANGLE + 2)) / 8
          IQUART(2)=0
          IQUART(3)=0
          IQUART(4)=0
        ELSE IF( ISYMM .EQ. 4 .OR. ISYMM .EQ. 20 )THEN
          NSOLAN= (NANGLE * (NANGLE + 2)) / 4
          IQUART(2)=0
          IQUART(4)=0
        ELSE IF( ISYMM .EQ. 2 .OR. ISYMM .EQ. 18 )THEN
          NSOLAN= (NANGLE * (NANGLE + 2)) / 4
          IQUART(3)=0
          IQUART(4)=0
        ELSE
          NSOLAN= (NANGLE * (NANGLE + 2)) / 2
        ENDIF
      ENDIF
      DDENST= 1.0/(NPAN*DENUSR)
      NANGLS=  (NANGLE * (NANGLE + 2)) / 2
      CALL XELEQN( 3, 0, ANGEQN )
      IANG= 0
      DO 15 JANG= 1, NANGLS
      DO 16 IPAN= 1, NPAN
         CALL    XELEQN( 3, NANGLE, ANGEQN )
   16 CONTINUE
        IF(IQUART(MOD(JANG-1,4)+1).NE.1 ) GO TO 15
        IANG= IANG+1
        DENSTY(IANG)= REAL(2*NSOLAN)
        ANGLES(1,IANG)= ANGEQN(1,1)
        ANGLES(2,IANG)= ANGEQN(2,1)
        ANGLES(3,IANG)= ANGEQN(3,1)
   15 CONTINUE
*
*     COPY ANGLES AND DENSITIES ON TEMPORARY TRACKING FILE
      WRITE(IFTEMP) ((ANGLES(IDIM,IANG),IDIM=1,NDIM),IANG=1,NSOLAN)
      WRITE(IFTEMP) (DENSTY(IANG)                   ,IANG=1,NSOLAN)
*
*     TO REINITIATE THE EQN ANGLES
      CALL XELEQN( 3, 0, ANGEQN )
      IF( NEXTGE.EQ.1 )THEN
         DDENST= 12.0*DDENST
      ENDIF
      WEIGHT= 0.5*DDENST/DBLE(NSOLAN)
      DO IANG= 1, NSOLAN
        DDENWT(IANG)=WEIGHT
      ENDDO
      NSOLMX= 0
      NDEBS= 0
      IF( IPRT.GT.1 )THEN
*
*        PREPARE & PRINT THE ESTIMATED NUMBER OF TRACKS
         NESTIM= 0
         DEPART= - (NPO2+1) * DP
         X     = DEPART
         DO 25 IX = 1, NPOINT
            X     = X + DP
            Y     = DEPART
            DO 20 IY = 1, NPOINT
               Y = Y + DP
               IF( ANORM2( X, Y ) .LE. R2CIRC ) NESTIM= NESTIM + 1
   20       CONTINUE
   25    CONTINUE
         WRITE(IOUT,'(1H )')
         WRITE(IOUT,'( 8H0ECHO =  ,I8,20H TRACKS/AXIS/ANGLE      )')
     >                          NPOINT
      WRITE(IOUT,'( 8H ECHO =  ,I8,25H TRACKS/CIRCLE/ANGLE          )')
     >                          NESTIM
         NESTIM= NESTIM * NPAN * NSOLAN
      WRITE(IOUT,'(30H ECHO = NEXT ODD    DENSITY >  ,F15.7,4H/CM2)')
     >                                                ODDNXT**2
      WRITE(IOUT,'( 8H ECHO = ,28H ESTIMATED NUMBER OF TRACKS=  ,I8 )')
     >                                                           NESTIM
*
*        PREPARE PRINTING WITH VARIABLE FORMAT
         WRITE(IOUT,'(1H )')
         WRITE(IOUT,'( 8H0ECHO = ,I3,27H SOLID ANGLES TO BE TRACKED )')
     >                   NANGLS
         NSOLMX= MIN(9, NANGLS/10)
         IREF1 = 0
         WRITE(IOUT,'( 1H0,10(I1,9X))') (IREF1, IZZ=0,NSOLMX)
         WRITE(IOUT,'( 1H ,10(I1,9X))') (MOD(IZZ,10), IZZ=0,NSOLMX)
         WRITE(IOUT,'( 2H 0)')
         TEDATA= '(1H+,TXXX,I1)'
      ENDIF
      IANG  = 0
      IANG0 = 0
      NTTRK = 0
      DO 290 IANGL= 1, NANGLS
         IANG=IANG+1
         IF(IQUART(MOD(IANG-1,4)+1).NE.1)THEN
*----
*  Do not track this angle because of the problem symmetry
*----
           LANGLE= .FALSE. 
         ELSE
*----
*  Track this angle
*----
           IANG0= IANG0+1
           LANGLE=.TRUE.
         ENDIF
         IF( IPRT.GT.1) THEN
           IF( MOD(IANGL,100) .EQ. 0 )THEN
             IREF1=IREF1+1
             NDEBS= NSOLMX+1
             NSOLMX=MIN(NDEBS+9, NANGLS/10)
             WRITE(IOUT,'( 1H0,10(I1,9X))')(IREF1,IZZ=NDEBS,NSOLMX)
             WRITE(IOUT,'( 1H ,10(I1,9X))')
     >            (MOD(IZZ,10),IZZ=NDEBS,NSOLMX)
             WRITE(IOUT,'( 2H 0)')
           ELSE
             IF( IPRT.GT.10000.AND.MOD(IANGL,100).NE.0 )THEN
               WRITE(IOUT,'( 1H ,10(I1,9X))')(IREF1,IZZ=NDEBS,NSOLMX)
               WRITE(IOUT,'( 1H ,10(I1,9X))')
     >               (MOD(IZZ,10),IZZ=NDEBS,NSOLMX)
               WRITE(IOUT,'( 2H  )')
             ENDIF
             WRITE(TEDATA(7:9),'(I3.3)') MOD(IANGL,100) + 2
             WRITE(IOUT,TEDATA) MOD(IANGL,10)
           ENDIF
         ENDIF
*
*     NPAN AXES DESCRIPTION (X=0.0, Y=0.0 & Z=0.0)
      DO 250 IPAN= 1, NPAN
*----
* Start tesp print
*      WRITE(IOUT,7001) IANGL,IPAN 
* 7001 FORMAT(' ANGLE = ',I8,5X,'PLAN =',I3)
* Finish test print
*----
         CALL    XELEQN( 3, NANGLE, ANGEQN )
         IF(.NOT.LANGLE) GO TO 250
         IF( NEXTGE.EQ.1 )THEN
            IF( IPAN.NE.2 ) GO TO 250
         ENDIF
         DO 30 I   = 1, 3
            N     = ICOORD(I)
            TRKDIR(N)= ANGEQN(N,1)
            INCR(I)= +1
            IF( TRKDIR(N) .LT. 0.0 ) INCR(I)= -1
*
*           MODIFY ANGLES TO TAKE INTO ACCOUNT DP
            ANGLE2(I)= DP * ANGLE2(I)
            ANGLE3(I)= DP * ANGLE3(I)
            IF( NCUTOF.NE.1 )THEN
               TCUTOF(I,1)= RCUTOF*( ANGLE2(I)+ ANGLE3(I) )
               TCUTOF(I,2)= RCUTOF*( ANGLE2(I)- ANGLE3(I) )
               TCUTOF(I,3)= -TCUTOF(I,2)
               TCUTOF(I,4)= -TCUTOF(I,1)
            ENDIF
*
*           DETERMINE THE ORIGIN OF ALL TRACKS
            TRKOR2(I)= BARY(I) - (NPO2+1)*(ANGLE2(I)+ANGLE3(I))
   30    CONTINUE
         DO 45 I   = 1, 3
            PROJC2(I)= 0.0
            DO 40 J   = 1, 3
               IF( I.EQ.J ) GO TO 40
               PROJC2(I)= PROJC2(I) + TRKDIR(J) * TRKDIR(J)
   40       CONTINUE
   45    CONTINUE
*
*        SCAN ALL TRACKS IN THE PERPENDICULAR PLANE
         DO 180 I2  = 1, NPOINT
            DO 50 J   = 1, 3
               TRKOR2(J)= TRKOR2(J) + ANGLE2(J)
               TRKORI(J)= TRKOR2(J)
   50       CONTINUE
            DO 170 I3  = 1, NPOINT
               ANN   = 0.0
               DO 60 J    = 1, 3
                  TRKORI(J)= TRKORI(J) + ANGLE3(J)
                  ANN= ANN + (TRKORI(J)-BARY(J))**2
   60          CONTINUE
*----
* Start tesp print
*               WRITE(IOUT,7002) I2,I3,(TRKORI(K)-BARY(K),K=1,NDIM)
*7002           FORMAT(' ORIGINE MESH:',I10,5X,I10,5X,3(F11.5))
* Finish test print
*----
               WZ=1.0D0
               DZ=(TRKORI(3)-BARY(3))/ZMAX
*----
* Start Z reflection symmetry
*
*               IF(ISYMM .GE. 16) THEN 
*                 IF (ABS(DZ) .LT. 1.0E-6) THEN
*                   WRITE(IOUT,'(A10)') 'ZERO Z    '
*                 ELSE IF(DZ .LT. 0.0) THEN
*                   WRITE(IOUT,'(A10)') 'NEGATIVE Z'
*                   GO TO 170
*                 ELSE
*                   WRITE(IOUT,'(A10)') 'POSITIVE Z'
*                   WZ=2.0
*                 ENDIF
*               ENDIF
* Finish Z reflection symmetry
*----
*
*              ELIMINATE TRACKS OUTSIDE CIRCUMSCRIBED CIRCLE
               IF( ANN.GT.R2CIRC ) GO TO 170
*
*              WHICH EXTERNAL SURFACES DO THIS TRACK CROSS ?
               NTTRK=NTTRK+1
               CALL XELLSR(  NDIM, NTOTCL, NSUR, MAXR, REMESH,
     >                      INDEL, MINDIM, MAXDIM, ICOORD, ICUR, INCR,
     >                      TRKORI, TRKDIR, TRKCUT, NSCUT, NCROS,
     >                      TOTLEN)
*
*              WHEN NOT SURFACES ARE CROSSED, ELIMINATE THE TRACK
               IF(NCROS.LT.2) GO TO 170
               DO 70 K= 1, NDIM
                  TRKBEG(K)= TRKCUT(K,1)
                  TRKEND(K)= TRKCUT(K,2)
   70          CONTINUE
               DO 75 K= 1, 4
                  NSBEG(K)= NSCUT(1)
                  NSEND(K)= NSCUT(2)
   75          CONTINUE
               IF( NCUTOF.NE.1 )THEN
                  DO 77 K= 1, 4
                     DO 76 K3= 1, 3
                        TORIC(K3)= TRKORI(K3)+TCUTOF(K3,K)
   76                CONTINUE
                     CALL XELLSR( NDIM, NTOTCL, NSUR, MAXR, REMESH,
     >                       INDEL, MINDIM, MAXDIM, ICOORD, ICUR, INCR,
     >                           TORIC, TRKDIR, TRKCUT, NSCUT, NCROS,
     >                           TOTXXX)
                     IF(NSCUT(1).NE.0) NSBEG(K)= NSCUT(1)
                     IF(NSCUT(2).NE.0) NSEND(K)= NSCUT(2)
   77             CONTINUE
               ENDIF
               CALL XELLIN(   NDIM, NTOTCL, MAXR, REMESH,
     >                NSUR,   NVOL,  INDEL, MINDIM, MAXDIM,
     >              ICOORD,   ICUR,   INCR, TRKBEG, TRKEND, TRKDIR,
     >              PROJC2, TOTLEN,
     >                CONV, LINMAX, LENGHT, NUMERO, LINE)
               NOTRAK= NOTRAK+1
*
               WRITE(IFTEMP) 1,LINE+2*NCUTOF,WEIGHT*WZ,IANG0,
     >                     (NSBEG(ISB),ISB=1,NCUTOF),
     >                     (NUMERO(I),I=1,LINE),
     >                     (NSEND(ISE),ISE=1,NCUTOF),
     >                     ( DBLE(1.0/NCUTOF),ISB=1,NCUTOF),
     >                     (LENGHT(I),I=1,LINE),
     >                     ( DBLE(1.0/NCUTOF),ISE=1,NCUTOF)
               IF( IPRT.GT.10000)THEN
                  WRITE(IOUT,6001) NOTRAK,
     >                         NCUTOF,(TRKBEG(JC),JC=1,3),
     >                         NCUTOF,(TRKEND(JC),JC=1,3),
     >                                (TRKDIR(JC),JC=1,3)
                 WRITE(IOUT,6002) (1.0/NCUTOF,NSBEG(ISB),ISB=1,NCUTOF),
     >                         (LENGHT(I), NUMERO(I),I=1,LINE),
     >                         (1.0/NCUTOF,NSEND(ISE),ISE=1,NCUTOF)
               ENDIF
  170       CONTINUE
  180    CONTINUE
  250 CONTINUE
  290 CONTINUE
      IF( IPRT.GT.1 )THEN
         WRITE(IOUT,'(1H )')
         WRITE(IOUT,'(27H0ECHO = TRACKING PROPERTIES )')
         WRITE(IOUT,'( 8H0ECHO = ,I3,20H ANGLES AND DENSITY:,
     >               F9.3,4H/CM2)')
     >                         NANGLE, DENUSR
         WRITE(IOUT,'( 8H0ECHO = ,I10,3H / ,I10,
     >                             23H TRACKS STORED ON TAPE ,I2/)')
     >                         NOTRAK,NTTRK,IFTEMP
      ENDIF
*
      RETURN
 6001 FORMAT(' #',I8,1P,' B',I1,'(',2(E10.2,','),E10.2,')',
     >                  ' E',I1,'(',2(E10.2,','),E10.2,')',
     >                  ' D(',2(E10.2,','),E10.2,')' )
 6002 FORMAT(1P,5(1X,E15.7,1X,I6))
      END
