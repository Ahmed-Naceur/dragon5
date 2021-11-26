*DECK XELLSR
      SUBROUTINE XELLSR(  NDIM,    NCP,   NSUR, MAXREM, REMESH,
     >                   INDEL, MINDIM, MAXDIM, ICOORD,   ICUR,   INCR,
     >                  TRKORI, TRKDIR, TRKCUT,  NSCUT,  NCROS,
     >                  TOTLEN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Find the beginning and ending surfaces crossed by a track.
*
*Copyright:
* Copyright (C) 1990 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* NDIM    number of dimension (2 or 3).                     
* NCP     number of cylindres of a type + 3.
* NSUR    number of surfaces.                               
* MAXREM  max number of real mesh values in REMESH.  
* REMESH  real mesh values (rect/cyl).                 
* INDEL   numbering of surfaces and zones.                    
* MINDIM  min index values for all axes (rect/cyl).    
* MAXDIM  max index values for all axes (rect/cyl).    
* ICOORD  principal axes direction (X/Y/Z) for meshes. 
* ICUR    current zonal location for a track segment.  
* INCR    increment direction for next track segment.  
* TRKORI  origin of a track.                           
* TRKDIR  direction of a track in all axes.            
*
*Parameters: output
* TRKCUT  points where track cut the domain.           
* NSCUT   surface where the track begins/ends.        
* NCROS   number of surface crossing.                        
* TOTLEN  total length of the track.                    
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*
      INTEGER NDIM, NCP, NSUR, MAXREM, NCROS
      REAL    TRKCUT(3,2), REMESH(MAXREM), TRKDIR(3), TRKORI(3), TOTLEN,
     >        TKBEG1, TKBEG2, R2BEG
      INTEGER MINDIM(NCP), MAXDIM(NCP), ICUR(NCP), INCR(NCP),
     >        ICOORD(NCP), IFACUT(2), ISFCUT(2),
     >        IORD(4), NSCUT(2), INDEL(4,*)
      INTEGER     IOUT
      PARAMETER ( IOUT=6 )
      REAL        ANORM2, CENTRE, A, B, XYZP2, CONBEG, CONEND, CON,
     >            XYZP1
      INTEGER     I, J, NUM, NCRBEG, NCREND, NP1, NUMP1, NP2, NUMP2,
     >            N, NUMP0, K, IBEGIN, KELSUR, KWW, IDM
*
      ANORM2(A,B)= A*A + B*B
      CENTRE(I,J)= REMESH( MAXDIM(I-1) + J )
      NUM(J)= J + 1 - NSUR
      NUMP2= 0
      CALL XDISET(IFACUT,2,0)
      CALL XDISET(ISFCUT,2,0)
*
      IF( NDIM.EQ.2 )THEN
         NP2= 3
         NUMP2 = ICOORD(NP2)
         XYZP2= 0.0
      ENDIF
*
*     IF THERE ARE NO CYLINDER AT ALL
      NSCUT(1)= 0
      NSCUT(2)= 0
      NCRBEG= 0
      NCREND= 0
      CONBEG=+1.0E+36
      CONEND=-1.0E+36
*
*     FING BEGINNING AND ENDING POINTS OF THE TRACK
      DO 75 N   = 1, NDIM
         NUMP0 = ICOORD(N  )
         IF( INCR(NUMP0).EQ.0 ) GO TO 75
         NP1   = MOD(N   ,NDIM)+1
         NUMP1 = ICOORD(NP1)
         IF( NDIM.EQ.3 )THEN
            NP2   = MOD(N+1 ,NDIM)+1
            NUMP2 = ICOORD(NP2)
         ENDIF
         DO 70 IDM = MINDIM(N), MAXDIM(N), MAXDIM(N)-MINDIM(N)
            CON   = (REMESH(IDM) - TRKORI(NUMP0)) / TRKDIR(NUMP0)
            XYZP1 = TRKORI(NUMP1) + CON * TRKDIR(NUMP1)
            IF( XYZP1.LT.REMESH(MINDIM(NP1)).OR.
     >          XYZP1.GT.REMESH(MAXDIM(NP1))) GO TO 70
            IF( NDIM.EQ.3 )THEN
               XYZP2 = TRKORI(NUMP2) + CON * TRKDIR(NUMP2)
               IF( XYZP2.LT.REMESH(MINDIM(NP2)).OR.
     >             XYZP2.GT.REMESH(MAXDIM(NP2))) GO TO 70
            ENDIF
            IF( CON.LT.CONBEG )THEN
               NCRBEG=1
               NCREND=MAX(1,NCREND)
               IFACUT(1)= NUMP0
               ISFCUT(1)= IDM
               IF( IDM.EQ.MINDIM(N) ) ISFCUT(1)= ISFCUT(1)-1
               CONBEG=CON
               TRKCUT(NUMP0,1)= REMESH(IDM)
               TRKCUT(NUMP1,1)= XYZP1
               TRKCUT(NUMP2,1)= XYZP2
            ENDIF
            IF( CON.GT.CONEND )THEN
               NCREND=2
               NCRBEG=MIN(2,NCRBEG)
               IFACUT(2)= NUMP0
               ISFCUT(2)= IDM
               IF( IDM.EQ.MINDIM(N) ) ISFCUT(2)= ISFCUT(2)-1
               CONEND=CON
               TRKCUT(NUMP0,2)= REMESH(IDM)
               TRKCUT(NUMP1,2)= XYZP1
               TRKCUT(NUMP2,2)= XYZP2
            ENDIF
   70    CONTINUE
   75 CONTINUE
      NCROS = NCREND + NCRBEG
      TOTLEN= CONEND - CONBEG
      IF( NCROS.EQ.0 ) GO TO 1000
      NCROS = NCREND + 1 - NCRBEG
*
*     FIND BEGINNING AND ENDING SURFACE NUMBERS
      DO 900 K= NCRBEG, NCREND
         DO 90 I   = 1, NDIM
            N      = ICOORD(I)
            ICUR(I)= MINDIM(I)
            DO 80 J = MINDIM(I), MAXDIM(I)-1
               IF(TRKCUT(N,K).GE.REMESH(J)) ICUR(I)= J
   80       CONTINUE
   90    CONTINUE
         ICUR(IFACUT(K))= ISFCUT(K)
         IBEGIN= MAXDIM(3) + 3
         DO 110 I  = 4, NCP
            N     = ICOORD(I)
            NP1   = MOD(N  ,3) + 1
            NP2   = MOD(N+1,3) + 1
            TKBEG1= CENTRE(I,1) - TRKCUT(NP1,K)
            TKBEG2= CENTRE(I,2) - TRKCUT(NP2,K)
            R2BEG = ANORM2(TKBEG1,TKBEG2)
            ICUR(I)  = IBEGIN - 1
            DO 100 J  = IBEGIN, MAXDIM(I)
               IF( R2BEG    .GE. REMESH(J) )ICUR(I)= J
  100       CONTINUE
            IBEGIN= MAXDIM(I) + 3
  110    CONTINUE
*
*        FIND IORD(4) FOR LOCATION IN THE INDEX VECTOR
         DO 115 I= 1,NCP
            IORD(MIN(4,I))= ICUR(I)
            IF( I.GT.3.AND.ICUR(I).LT.MAXDIM(I)) GOTO 116
  115    CONTINUE
         IORD(4)= 0
  116    CONTINUE
*
*        FIND NSCUT=BEGINNING/ENDING SURFACE #S
         KELSUR= NSUR
         INDEL(1,NUM(0))= IORD(1)
         INDEL(2,NUM(0))= IORD(2)
         INDEL(3,NUM(0))= IORD(3)
         INDEL(4,NUM(0))= IORD(4)
  880    CONTINUE
            IF( IORD(1).EQ.INDEL(1,NUM(KELSUR)).AND.
     >          IORD(2).EQ.INDEL(2,NUM(KELSUR)).AND.
     >          IORD(3).EQ.INDEL(3,NUM(KELSUR)).AND.
     >          IORD(4).EQ.INDEL(4,NUM(KELSUR)) ) GO TO 890
            KELSUR= KELSUR + 1
         GO TO 880
  890    NSCUT(K)= KELSUR
         IF( KELSUR.EQ.0 )THEN
            WRITE(IOUT,*) '         BAD SURFACE IDENTIFICATION'
            WRITE(IOUT,*) ' NSCUT=', NSCUT(K)
            WRITE(IOUT,*) 'TRKCUT=', (TRKCUT(KWW,K),KWW=1,3)
            WRITE(IOUT,*) '  IORD=', IORD
            WRITE(IOUT,*) '  ICUR=', (ICUR(KWW),KWW=1,NCP)
            CALL XABORT('XELLSR: BAD SURFACE IDENTIFICATION')
         ENDIF
  900 CONTINUE
*
 1000 RETURN
      END
