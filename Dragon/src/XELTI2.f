*DECK XELTI2
      SUBROUTINE XELTI2( IPRT,IFTEMP,NANGLE,DENUSR,ISYMM,ANGLES,DENSTY,
     >                   NTOTCL,MAXR,REMESH,LINMAX,RCUTOF,
     >                   NSUR,NVOL,INDEL,MINDIM,
     >                   MAXDIM,ICOORD,INCR,ICUR,TRKBEG,CONV,TRKDIR,
     >                   LENGHT,NUMERO,DDENWT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Construct the sequential tape that will contain tracks for  
* isotropic BC for 2-D calculation. 
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
* ISYMM   flag for symetry:        
*         = 2 reflection plane normal to X axis; 
*         = 4 reflection plane normal to Y axis;
*         = 8 reflection plane normal to X and Y axis.
* ANGLES  3D angle values.                             
* DENSTY  density of tracks angle by angle.            
* NTOTCL  number of cylindres of a type + 2.                
* MAXR    max number of real mesh values in REMESH.  
* REMESH  real mesh values (rect/cyl).                 
* LINMAX  max. number of track segments in a single track.  
* RCUTOF  cutof for corner tracking.                    
* NSUR    number of surfaces.                               
* NVOL    number of zones.                                  
* INDEL   numbering of surfaces & zones.                    
* MINDIM  min index values for all axes (rect/cyl).    
* MAXDIM  max index values for all axes (rect/cyl).    
* ICOORD  principal axes direction (X/Y/Z) for meshes.
* ICUR    current zonal location for a track segment.  
* INCR    increment direction for next track segment.  
* TRKBEG  position where a track begins.               
* CONV    segments of tracks.                          
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
      INTEGER       IPRT,IFTEMP,NANGLE,NTOTCL,MAXR,LINMAX
      REAL          TRKBEG(NTOTCL), TRKDIR(NTOTCL), CONV(NTOTCL),
     >              REMESH(MAXR),DENUSR
      DOUBLE PRECISION DENSTY(NANGLE),ANGLES(2,NANGLE),LENGHT(LINMAX),
     >                 DDENWT(NANGLE)
      INTEGER       MINDIM(NTOTCL), MAXDIM(NTOTCL), ICUR(NTOTCL),
     >              ICOORD(NTOTCL), INCR(NTOTCL), NUMERO(LINMAX),
     >              INDEL(4,*)
*
*     DECLARE       LOCAL VARIABLES
      REAL          TRKEND(3), TRKORI(3), BARY(2), PROJC2(3), RCUTOF,
     >              ANGEQN(2,2), ANGLE2(2), TRKCUT(3,2),
     >              TCUTOF(2,2), TORIC(3), DENLIN, R2CIRC, RCIRC, DP,
     >              ODDNXT, SURTOT, VOLTOT, DDENST, ANN, TOTLEN, TOTXXX
      DOUBLE PRECISION WEIGHT,WTA,WTD
      INTEGER       NSBEG(4), NSEND(4), NSCUT(2), NDIM, I, NPOINT,
     >              NPO2, NCUTOF, NOTRAK, IANGL, IDIM, NANGLS, NESTIM,
     >              NSOLMX, IREF1, IZZ, IANG, NTTRK, NDEBS, I2, J,
     >              NSUR, NVOL, NCROS, K, K3, LINE, ISB, ISE, JC
      INTEGER       ISYMM
      REAL          WGA,WLA,WGD,WLD
      EQUIVALENCE ( ANGEQN(1,2), ANGLE2 )
      INTEGER       IOUT
      PARAMETER   ( IOUT=6 )
      CHARACTER     TEDATA*13
*
*     ONE WEIGHT FOR ALL TRACKS
*     DENLIN= # OF TRACKS / CM
      DENLIN= DENUSR
      NDIM= 2
      PROJC2(1)= 0.0
      PROJC2(2)= 0.0
      PROJC2(3)= 1.0
      TRKBEG(3)= 0.0
      TRKDIR(3)= 0.0
      TRKEND(3)= 0.0
*
*     COMPUTE THE CIRCUMSCRIBED RADIUS
*             THE COORDINATE FOR THE TRUE CENTER OF THE CELL
      R2CIRC= 0.0
      DO 10 I = 1, 2
         BARY(I)= 0.5 * (REMESH(MAXDIM(I)) + REMESH(MINDIM(I)))
         R2CIRC= R2CIRC
     >         + (REMESH(MAXDIM(I)) - REMESH(MINDIM(I)))**2
   10 CONTINUE
      R2CIRC= 0.25 * R2CIRC
      RCIRC = SQRT(R2CIRC)
*
*     NPOINT= # OF TRACKS ALONG THE PERPENDICULAR AXIS
      NPOINT= INT( 2. * RCIRC * DENLIN )
****** BEWARE ***** BEWARE ***** BEWARE ***** BEWARE ***** BEWARE *****
****** CHANGE THIS "NPOINT" PARAMETER HAS TREMENDOUS EFFECTS ON TRACKING
****** BEWARE ***** BEWARE ***** BEWARE ***** BEWARE ***** BEWARE *****
*
*     POSSIBLE OTHER CHOICES (EXPLORED WITHOUT SUCCESS) ARE ==>
*1-)  NPOINT=  INT( 2. * RCIRC * DENLIN ) + 1
*2-)  NPOINT= NINT( 2. * RCIRC * DENLIN )
*3-)  NPOINT= NINT( 2. * RCIRC * DENLIN ) + 1
*
*     KEEP "NPOINT" ODD & CORRECT DENSITY
      NPO2  = NPOINT / 2
      NPOINT= 2 * NPO2 + 1
      DP    = 2. * RCIRC / NPOINT
      DENLIN= 1. / DP
      DENUSR= DENLIN
      ODDNXT= (2*NPO2+3) / (2.*RCIRC)
      IF( RCUTOF.EQ.0.0 )THEN
         NCUTOF= 1
      ELSE
         NCUTOF= 2
      ENDIF
*
      NOTRAK= 0
      SURTOT= 0.0
      VOLTOT= 0.0
      DDENST= 1.0 / DENUSR
      WEIGHT= 0.5*DDENST/DBLE(NANGLE)
      CALL XELEQN( 2, 0, ANGEQN )
      DO 15 IANGL= 1, NANGLE
         CALL    XELEQN( 2, NANGLE, ANGEQN )
         DENSTY(IANGL)= REAL(2*NANGLE)
         ANGLES(1,IANGL)= ANGEQN(1,1)
         ANGLES(2,IANGL)= ANGEQN(2,1)
         DDENWT(IANGL)=WEIGHT
   15 CONTINUE
*----
*  Optimize tracking by taking into account 
*  symmetry of geometry
*  1) Do not track symmetric lines
*  2) Nodify weight as required for tracks droped
*----
      WGA=1.0
      WLA=1.0
      WGD=1.0
      WLD=1.0
      NANGLS= NANGLE
      IF(ISYMM .EQ. 2 .OR. 
     >   ISYMM .EQ. 4 .OR.
     >   ISYMM .EQ. 8 ) THEN
        NANGLS= (NANGLE+1)/2
        WGA=2.0
        DO 200 IANGL=1,NANGLS-1
          DENSTY(IANGL)=0.5*DENSTY(IANGL)
 200    CONTINUE
        IF(2*NANGLS .EQ. NANGLE) THEN
          WLA=WGA
          DENSTY(NANGLS)=0.5*DENSTY(NANGLS)
        ENDIF
      ENDIF
      IF(ISYMM .EQ. 8 ) THEN
        WGD=2.0
        NPOINT=(NPOINT+1)/2 
      ENDIF
*
*     COPY ANGLES AND DENSITIES ON TEMPORARY TRACKING FILE
      WRITE(IFTEMP) ((ANGLES(IDIM,IANGL),IDIM=1,NDIM),IANGL=1,NANGLS)
      WRITE(IFTEMP)  (DENSTY(IANGL)                  ,IANGL=1,NANGLS)
*
*     TO REINITIATE THE EQN ANGLES ROUTINE
      CALL XELEQN( 2, 0, ANGEQN )
*
      NSOLMX= 0
      NDEBS= 0
      IF( IPRT.GT.1 )THEN
*
*        PREPARE & PRINT THE ESTIMATED NUMBER OF TRACKS
         NESTIM= NPOINT * NANGLS
         WRITE(IOUT,'(1H )')
         WRITE(IOUT,'( 8H ECHO =  ,I8,15H TRACKS/ANGLE      )')
     >                          NPOINT
      WRITE(IOUT,'(30H ECHO = NEXT ODD    DENSITY >  ,F15.7,3H/CM )')
     >                                                ODDNXT
      WRITE(IOUT,'( 8H ECHO = ,28H ESTIMATED NUMBER OF TRACKS=   ,I8)')
     >                                                           NESTIM
*
*        PREPARE PRINTING WITH VARIABLE FORMAT
         WRITE(IOUT,'(1H )')
         WRITE(IOUT,6003) NANGLE
         IF(NANGLS .NE. NANGLE) THEN
           WRITE(IOUT,6004) NANGLS
         ENDIF
         NSOLMX= MIN(9, NANGLS/10)
         IREF1=0
         WRITE(IOUT,'( 1H0,10(I1,9X))') (IREF1, IZZ=0,NSOLMX)
         WRITE(IOUT,'( 1H ,10(I1,9X))') (MOD(IZZ,10), IZZ=0,NSOLMX)
         WRITE(IOUT,'( 2H 0)')
         TEDATA= '(1H+,TXXX,I1)'
      ENDIF
      IANG=0
      NTTRK=0
      DO 290 IANGL= 1, NANGLS
          IANG=IANG+1
          WTA=WGA
          IF(IANGL .EQ. NANGLS) THEN
            WTA=WLA
          ENDIF 
          IF( IPRT.GT.1 )THEN
            IF( MOD(IANG,100) .EQ. 0 )THEN
              IREF1=IREF1+1
              NDEBS= NSOLMX+1
              NSOLMX=MIN(NDEBS+9, NANGLS/10)
              WRITE(IOUT,'( 1H0,10(I1,9X))')(IREF1,IZZ=NDEBS,NSOLMX)
              WRITE(IOUT,'( 1H ,10(I1,9X))')
     >             (MOD(IZZ,10),IZZ=NDEBS,NSOLMX)
              WRITE(IOUT,'( 2H 0)')
            ELSE
              IF( IPRT.GT.10000.AND.MOD(IANG,100).NE.0 )THEN
                 WRITE(IOUT,'( 1H ,10(I1,9X))')(IREF1,IZZ=NDEBS,NSOLMX)
                 WRITE(IOUT,'( 1H ,10(I1,9X))')
     >                (MOD(IZZ,10),IZZ=NDEBS,NSOLMX)
                 WRITE(IOUT,'( 2H  )')
              ENDIF
              WRITE(TEDATA(7:9),'(I3.3)') MOD(IANG,100) + 2
              WRITE(IOUT,TEDATA) MOD(IANG,10)
            ENDIF
          ENDIF
          CALL    XELEQN( 2, NANGLE, ANGEQN )
*
          DO 40 I   = 1, 2
              TRKDIR(I)= ANGEQN(I,1)
*
*             DETERMINE THE DIRECTION OF INCREASE FOR TRKDIR
              INCR(I)= +1
              IF( TRKDIR(I) .LT. 0.0 ) INCR(I)= -1
*
*             MODIFY PERPENDICULAR ANGLES TO TAKE INTO ACCOUNT DP
              ANGLE2(I)= DP * ANGLE2(I)
              IF( NCUTOF.NE.1 )THEN
                 TCUTOF(I,1)= RCUTOF*ANGLE2(I)
                 TCUTOF(I,2)= -TCUTOF(I,1)
              ENDIF
*
*             DETERMINE THE ORIGINE OF ALL TRACKS
              TRKORI(I)= BARY(I) -  (NPO2+1) * ANGLE2(I)
   40     CONTINUE
          DO 180 I2  = 1, NPOINT
              WTD=WGD
              IF(I2 .EQ. NPOINT) THEN
                WTD=WLD
              ENDIF 
              ANN   = 0.0
              DO 50 J   = 1, 2
                  TRKORI(J)= TRKORI(J) + ANGLE2(J)
                  ANN= ANN + (TRKORI(J)-BARY(J))**2
   50         CONTINUE
*
*             ELIMINATE TRACKS OUTSIDE CIRCUMSCRIBED CIRCLE
              IF( ANN.GT.R2CIRC ) GO TO 180
*----
* Start test print
*
*             WRITE(IOUT,7002) I2,I3,(TRKORI(JJ),JJ=1,3)
*7002         FORMAT(' ORIGINE MESH:',I10,5X,I10,5X,3(F11.5)) 
*  Finish test print
*----
*
*             WHICH EXTERNAL SURFACES DO THIS TRACK CROSS ?
              NTTRK=NTTRK+1
              CALL XELLSR( NDIM, NTOTCL, NSUR, MAXR, REMESH,
     >                    INDEL, MINDIM, MAXDIM, ICOORD, ICUR, INCR,
     >                    TRKORI, TRKDIR, TRKCUT, NSCUT, NCROS,
     >                    TOTLEN)
*
*             WHEN NOT SURFACES ARE CROSSED, ELIMINATE THE TRACK
              IF( NCROS.LT.2 ) GO TO 180
              DO 70 K= 1, NDIM
                 TRKBEG(K)= TRKCUT(K,1)
                 TRKEND(K)= TRKCUT(K,2)
   70         CONTINUE
              DO 75 K= 1, 4
                 NSBEG(K)= NSCUT(1)
                 NSEND(K)= NSCUT(2)
   75         CONTINUE
              IF( NCUTOF.NE.1 )THEN
                 DO 77 K= 1, 2
                    DO 76 K3= 1, 2
                       TORIC(K3)= TRKORI(K3)+TCUTOF(K3,K)
   76               CONTINUE
                    CALL XELLSR( NDIM, NTOTCL, NSUR, MAXR, REMESH,
     >                      INDEL, MINDIM, MAXDIM, ICOORD, ICUR, INCR,
     >                          TORIC, TRKDIR, TRKCUT, NSCUT, NCROS,
     >                          TOTXXX)
                    IF(NSCUT(1).NE.0) NSBEG(K)= NSCUT(1)
                    IF(NSCUT(2).NE.0) NSEND(K)= NSCUT(2)
   77            CONTINUE
              ENDIF
              CALL XELLIN(  NDIM, NTOTCL, MAXR, REMESH,
     >                NSUR,   NVOL,  INDEL, MINDIM, MAXDIM,
     >              ICOORD,   ICUR,   INCR, TRKBEG, TRKEND, TRKDIR,
     >              PROJC2, TOTLEN,
     >                CONV, LINMAX, LENGHT, NUMERO, LINE)
              NOTRAK= NOTRAK+1
*
              WRITE(IFTEMP) 1,LINE+2*NCUTOF,WEIGHT*WTA*WTD,IANG,
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
  180   CONTINUE
  290 CONTINUE
      IF( IPRT.GT.1 )THEN
         WRITE(IOUT,'(1H )')
         WRITE(IOUT,'(27H0ECHO = TRACKING PROPERTIES )')
         WRITE(IOUT,'( 8H0ECHO = ,I3,20H ANGLES AND DENSITY:,
     >                F9.3,4H/CM )')
     >                          NANGLS, DENUSR
         WRITE(IOUT,'( 8H0ECHO = ,I10,3H / ,I10,
     >                              23H TRACKS STORED ON TAPE ,I2/)')
     >                          NOTRAK,NTTRK,IFTEMP
      ENDIF
      RETURN
*
 6001 FORMAT(' #',I8,1P,' B',I1,'(',2(E10.2,','),E10.2,')',
     >                  ' E',I1,'(',2(E10.2,','),E10.2,')',
     >                  ' D(',2(E10.2,','),E10.2,')' )
 6002 FORMAT(1P,5(1X,E15.7,1X,I6))
 6003 FORMAT(' '/
     >       ' ECHO = ',I3,' ANGLES TO BE TRACKED IN RANGE 0 TO PI')
 6004 FORMAT('        ',I3,' ANGLES IN RANGE 0 TO PI/2 AFTER SYMMETRY')
      END
