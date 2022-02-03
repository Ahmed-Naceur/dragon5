*DECK XL3TI3
      SUBROUTINE XL3TI3( IPRT,NANGLE,DENUSR,ISYMM,ANGLES,DENSTY,
     >                   NTOTCL,NEXTGE,MAXR,REMESH,LINMAX,RCUTOF,
     >                   NSUR,NVOL,INDEL,MINDIM,
     >                   MAXDIM,ICOORD,INCR,ICUR,TRKBEG,CONV,TRKDIR,
     >                   LENGHT,NUMERO,NPIJ,NGRP,SIGTAL,SWVOID,NORE,
     >                   NRMV,VOLTRK,KEYMRG,NSOUT,NREG,NPSYS,DBLPIJ )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Construct the sequential tape that will contain tracks for 
* isotropic BC 3-D calculations.
*
*Copyright:
* Copyright (C) 1996 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* IPRT    intermediate printing level for output.      
* NANGLE  number of angles used in the tracking process.    
* DENUSR  density of tracks in the plane perpendicular 
*         to the tracking angles.                 
* ISYMM   flag for symetry (1/0 for on/off):         
*         = 2 reflection plane normal to X axis;
*         = 4 reflection plane normal to Y axis;
*         = 8 reflection plane normal to X and Y axis;  
*         =16 reflection plane normal to Z axis;
*         =18 reflection plane normal to X and Z axis;
*         =20 reflection plane normal to Y and Z axis;
*         =24 reflection plane normal to X, Y and Z axis.  
* ANGLES  3D angle values.                             
* DENSTY  density of tracks angle by angle.            
* NTOTCL  number of cylindres of a type + 3.                
* NEXTGE  for tubez, nextge=1                          
* MAXR    max number of real mesh values in REMESH.  
* REMESH  real mesh values (rect/cyl).                 
* LINMAX  max. number of track segments in a single track.  
* RCUTOF  cutof for corner tracking(0.25 suggested).  
* NSUR    number of surfaces.                               
* NVOL    number of zones.                                  
* INDEL   numbering of surfaces and zones.                    
* MINDIM  min index values for all axes (rect/cyl).    
* MAXDIM  max index values for all axes (rect/cyl).    
* ICOORD  principal axes direction (X/Y/Z) for meshes. 
* ICUR    current zonal location for a track segment.  
* INCR    increment direction for next track segment.  
* TRKBEG  position where a track begins.               
* CONV  s  egments of tracks.                          
* TRKDIR  direction of a track in all axes.            
* LENGHT  relative lenght of each segment in a track.  
* NUMERO  material identification of each track segment.
* NPIJ    number of probabilities in one group         
* NORE    track normalization (-1 yes; 1 no).            
* NRMV    volume factors only (0 no; 1 yes).            
* NGRP    number of groups.
* SIGTAL  total XS.
* SWVOID  flag for void regions.
* KEYMRG  merge keys.
* NSOUT   number of outer surfaces.
* NREG    number of regions.
* NPSYS   undefined.
* DBLPIJ  collision probabilities.
*
*Parameters: output
* VOLTRK  volume factors.
*
*Reference:
* R.Roy, A. Hebert and G. Marleau
* A transport method for treating 3-d lattices of heterogeneous cells,
* Nuclear Science and Engineering, 101, 217 (1989). 
*
*-----------------------------------------------------------------------
*
      IMPLICIT      NONE
*
      INTEGER       IPRT,NANGLE,NTOTCL,NEXTGE,MAXR,LINMAX,
     >              NSUR,NVOL,NPIJ,NGRP,NORE,NRMV,NSOUT,NREG
      REAL          DENUSR
*
      REAL          TRKBEG(NTOTCL), TRKDIR(NTOTCL), CONV(NTOTCL),
     >              REMESH(MAXR), DENSTY(*), RCUTOF, TRKEND(3),
     >              TRKORI(3), TRKOR2(3), PROJC2(3), ANGEQN(3,3),
     >              ANGLE2(3), ANGLE3(3), BARY(3),   TRKCUT(3,2),
     >              TCUTOF(3,4), TORIC(3), ANGLES(3,*),
     >              SIGTAL(NSOUT:NREG,NGRP)
      DOUBLE PRECISION LENGHT(LINMAX)
      DOUBLE PRECISION VOLTRK(NSUR:NVOL,0:*),DBLPIJ(NPIJ,NGRP)
      INTEGER       MINDIM(NTOTCL), MAXDIM(NTOTCL), ICUR(NTOTCL),
     >              ICOORD(NTOTCL), INCR(NTOTCL), NUMERO(LINMAX),
     >              INDEL(4,*), NSBEG(4), NSEND(4), NSCUT(2),
     >              KEYMRG(NSUR:NVOL), NCSEG, NOLDS, NNEWS, NPSYS(NGRP)
      INTEGER       ISYMM
      INTEGER       IQUART(4)
      LOGICAL       LANGLE, SWVOID
      EQUIVALENCE ( ANGEQN(1,2), ANGLE2 ), ( ANGEQN(1,3), ANGLE3 )
      CHARACTER     TEDATA*13
      INTEGER       NPAN, IOUT
      DOUBLE PRECISION DZERO
      PARAMETER   ( NPAN=3, IOUT=6, DZERO=0.0D0 )
*
      INTEGER       NDIM, NSOUTM, I, J, NPOINT, NPO2, NCUTOF, LINM2,
     >              NOTRAK, NSOLAN, IVS, IANG, IGRP, ISB,
     >              NANGLS, IPAN, NESTIM, JANG, IX, IY, IREF1, I2, I3,
     >              NSOLMX, NTTRK, IANG0, NDEBS, N, IZZ, IANGL, K, K3,
     >              NCROS, LINE, JC, IL
      REAL          ANORM2, A, B, DENLIN, RCIRC, R2CIRC, DP, SURTOT,
     >              VOLTOT, WEIGHT, DEPART, X, Y, ANN, ODDNXT, TOTLEN,
     >              TOTXXX, DDENST
      DOUBLE PRECISION FCUTOF
*
*
      ANORM2( A, B ) = A*A + B*B
      NDIM= 3
      NSOUTM= -NSOUT
      NSOLMX= 0
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
            CALL XABORT('XL3TI3: TUBEZ NOT SUPPORTED')
         ELSE
            R2CIRC =  R2CIRC
     >             + (REMESH(MAXDIM(I)) - REMESH(MINDIM(I)))**2
         ENDIF
   10 CONTINUE
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
      LINM2= LINMAX-2*NCUTOF
      FCUTOF= 1.0/DBLE(NCUTOF)
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
*
*     INITIALIZE NORMALIZED FACTORS
      IF( NRMV.EQ.1 )THEN
         DO 12 IVS= NSUR, NVOL
           DO 11 IANG= 0, NSOLAN
              VOLTRK(IVS,IANG)= DZERO
   11      CONTINUE
   12    CONTINUE
      ELSE
         DO 23 IGRP= 1, NGRP
           DO 21 IVS= 1,NPIJ
             DBLPIJ(IVS,IGRP)= DZERO
   21      CONTINUE
   23    CONTINUE
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
*     TO REINITIATE THE EQN ANGLES
      CALL XELEQN( 3, 0, ANGEQN )
      IF( NEXTGE.EQ.1 )THEN
         DDENST= 12.0*DDENST
      ENDIF
      WEIGHT= 0.5*DDENST/REAL(NSOLAN)
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
     >                   NSOLAN
         NSOLMX= MIN(9, NSOLAN/10)
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
         LANGLE= .FALSE.
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
         IF(( IPRT.GT.1).AND.( MOD(IANGL,100) .EQ. 0 ))THEN
           IREF1=IREF1+1
           NDEBS= NSOLMX+1
           NSOLMX=MIN(NDEBS+9, NANGLS/10)
           WRITE(IOUT,'( 1H0,10(I1,9X))')(IREF1,IZZ=NDEBS,NSOLMX)
           WRITE(IOUT,'( 1H ,10(I1,9X))')
     >          (MOD(IZZ,10),IZZ=NDEBS,NSOLMX)
           WRITE(IOUT,'( 2H 0)')
         ENDIF
*
*     NPAN AXES DESCRIPTION (X=0.0, Y=0.0 & Z=0.0)
      DO 250 IPAN= 1, NPAN
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
*
*              ELIMINATE TRACKS OUTSIDE CIRCUMSCRIBED CIRCLE
               IF( ANN.GT.R2CIRC ) GO TO 170
*
*              WRITE(IOUT,7002) I2,I3,(TRKORI(JJ),JJ=1,3)
*7002          FORMAT(' ORIGINE MESH:',I10,5X,I10,5X,3(F11.5))
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
     >              PROJC2, TOTLEN,   CONV,  LINM2,
     >              LENGHT(NCUTOF+1), NUMERO(NCUTOF+1), LINE)
               NOTRAK= NOTRAK+1
*
               DO 78 ISB=1,NCUTOF
                  NUMERO(ISB)= NSBEG(ISB)
                  LENGHT(ISB)= FCUTOF
                  NUMERO(NCUTOF+LINE+ISB)= NSEND(ISB)
                  LENGHT(NCUTOF+LINE+ISB)= FCUTOF
   78          CONTINUE
               LINE= LINE+2*NCUTOF
               IF( IPRT.GT.10000)THEN
                  WRITE(IOUT,6001) NOTRAK,
     >                         NCUTOF,(TRKBEG(JC),JC=1,3),
     >                         NCUTOF,(TRKEND(JC),JC=1,3),
     >                                (TRKDIR(JC),JC=1,3)
                  WRITE(IOUT,6002) (LENGHT(I),NUMERO(I),I=1,LINE)
               ENDIF
               IF( NRMV.EQ.1 )THEN
                  DO 301 I= 1, LINE
                     VOLTRK(NUMERO(I),IANG0)= VOLTRK(NUMERO(I),IANG0)
     >                                  + DBLE(WEIGHT) * LENGHT(I)
  301             CONTINUE
               ELSE
                  IF( NORE.EQ.-1 )THEN
                     DO 302 I= 1, LINE
                        IF( NUMERO(I).GT.0 )THEN
                           LENGHT(I) = LENGHT(I) *
     >                                 SNGL(VOLTRK(NUMERO(I),IANG0))
                        ENDIF
  302                CONTINUE
                  ENDIF
                  DO 303 I= 1, LINE
                     NUMERO(I)= KEYMRG(NUMERO(I))
  303             CONTINUE
*
*                 START MODIFICATIONS 98/06 (G.M.)
*                 COMPRESS TRACKING FILE FOR
*                          SUCCESSIVE IDENTICAL REGIONS
                  NOLDS=NUMERO(1)
                  NCSEG=1
                  DO 304 IL = 2, LINE
                     NNEWS=NUMERO(IL)
                     IF( NNEWS.LT.0 .OR. NNEWS.NE.NOLDS )THEN
                        NOLDS=NNEWS
                        NCSEG=NCSEG+1
                        NUMERO(NCSEG)=NUMERO(IL)
                        LENGHT(NCSEG)=LENGHT(IL)
                     ELSEIF( NNEWS.EQ.NOLDS )THEN
                        LENGHT(NCSEG)=LENGHT(NCSEG)+LENGHT(IL)
                     ENDIF
  304             CONTINUE
                  CALL QIJI3D(NREG,NSOUTM,NPIJ,NGRP,LINMAX,NCUTOF,
     >                  SWVOID,NCSEG,WEIGHT,NUMERO,LENGHT,
     >                  SIGTAL,NPSYS,DBLPIJ)
*                 END MODIFICATIONS 98/06 (G.M.)
*
               ENDIF
  170       CONTINUE
  180    CONTINUE
  250 CONTINUE
  290 CONTINUE
*
      IF( IPRT.GT.1 )THEN
         WRITE(IOUT,'(1H )')
         WRITE(IOUT,'(27H0ECHO = TRACKING PROPERTIES )')
         WRITE(IOUT,'( 8H0ECHO = ,I3,20H ANGLES AND DENSITY:,
     >               F9.6,4H/CM2)')
     >                         NANGLE, DENUSR
         WRITE(IOUT,'( 8H0ECHO = ,I10,3H / ,I10,
     >                             27H TRACKS NOT STORED ON TAPE /)')
     >                         NOTRAK,NTTRK
      ENDIF
*
      RETURN
*
 6001 FORMAT(' #',I8,1P,' B',I1,'(',2(E10.2,','),E10.2,')',
     >                  ' E',I1,'(',2(E10.2,','),E10.2,')',
     >                  ' D(',2(E10.2,','),E10.2,')' )
 6002 FORMAT(1P,5(1X,E15.7,1X,I6))
      END
