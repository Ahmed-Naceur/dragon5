*DECK TLMVPL
      FUNCTION TLMVPL(NDIM,NANGL,NPLOTS,IPLOT,IPLP,DPLPR,DANGLT,XYZL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To verify that the plane selected crosses the geometry.
*
*Copyright:
* Copyright (C) 2006 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):
* C. Plamondon, G. Marleau
*
*Parameters: input
* NDIM    number of dimensions for problem.
* NANGL   number of direction for tracking.
* NPLOTS  number of plots.
* IPLOT   plot number being processed.
* IPLP    integer plot parameters.
* DPLPR   real plot parameters.
* DANGLT  track directions.
* XYZL    mesh limits.
*
*Parameters: output.
* TLMVPL  flag to indicate intersection with
*         no intersection if TLPVPL<0.
*
*----------
*
*
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          NDIM,NANGL,NPLOTS,IPLOT
      INTEGER          IPLP(6,NPLOTS)
      DOUBLE PRECISION DPLPR(4,NPLOTS),DANGLT(NDIM,NANGL),XYZL(2,3)
*----
*  Function type
*----
      INTEGER          TLMVPL
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='TLMVPL')
      DOUBLE PRECISION DZERO,DONE
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0)
*----
*  Local variables
*----
      INTEGER           IDIR
      DOUBLE PRECISION  DPLP(4),A,B,C,D,DROITE,G(4),
     >                  XMIN,XMAX,YMIN,YMAX,ZMIN,ZMAX
*----
*  Processing starts:
*  Initialize TLMVPL
*  Verification que le plan croise bien la region 3D
*----
      IF(ABS(IPLP(1,IPLOT)) .EQ. 3) THEN
*----
*  PLANA
        DO IDIR=1,4
          DPLP(IDIR)=DPLPR(IDIR,IPLOT)
        ENDDO
*----
      ELSE IF(ABS(IPLP(1,IPLOT)) .EQ. 4) THEN
*----
*  PLANP lines
*----
        DO IDIR=1,NDIM
          DPLP(IDIR)=DANGLT(IDIR,IPLP(2,IPLOT))
        ENDDO
        DO IDIR=NDIM+1,3
          DPLP(IDIR)=DZERO
        ENDDO
      ENDIF
      A=DPLP(1)
      B=DPLP(2)
      C=DPLP(3)
      D=DPLP(4)
      XMIN=XYZL(1,1)
      XMAX=XYZL(2,1)
      YMIN=XYZL(1,2)
      YMAX=XYZL(2,2)
      ZMIN=XYZL(1,3)
      ZMAX=XYZL(2,3)
      TLMVPL=0
      IF(NDIM .EQ. 3) THEN
*----
*  Verification que le plan
*  A*B+B*Y+C*Z=D
*  croise bien la region 3D
*----
        IF(A .EQ. DZERO) THEN
          IF(B .EQ. DZERO) THEN
            IF(C .EQ. DZERO) THEN
              TLMVPL=-8
              WRITE(IOUT,9000) NAMSBR,A,B,C,D
            ELSE
              DROITE=D/C
              IF(DROITE .LE. ZMAX .AND. DROITE .GE. ZMIN) THEN
                TLMVPL=1
              ELSE
                TLMVPL=-1
                WRITE(IOUT,9001) NAMSBR,'ZZ',A,B,C,D
              ENDIF
            ENDIF
          ELSE IF(C .EQ. DZERO) THEN
            DROITE=D/B
            IF(DROITE .LE. YMAX .AND. DROITE .GE. YMIN) THEN
              TLMVPL=1
            ELSE
              TLMVPL=-2
              WRITE(IOUT,9001) NAMSBR,'YY',A,B,C,D
            ENDIF
          ELSE
            G(1)=(D-B*YMIN)/C
            G(2)=(D-B*YMAX)/C
            G(3)=(D-C*ZMIN)/B
            IF(G(1) .LT. ZMAX .AND. G(1) .GT. ZMIN) THEN
              TLMVPL=1
            ELSE IF(G(2) .LT. ZMAX .AND. G(2) .GT. ZMIN) THEN
              TLMVPL=1
            ELSE IF(G(3) .LT. YMAX .AND. G(3) .GT. YMIN) THEN
              TLMVPL=1
            ELSE
              G(1)=(D-(B*(YMIN+YMAX)/2))/C
              G(2)=(D-(C*(ZMIN+ZMAX)/2))/B
              IF(G(1) .LE. ZMAX .AND. G(1) .GE. ZMIN) THEN
                TLMVPL=1
              ELSE IF(G(2) .LE. YMAX .AND. G(2) .GE. YMIN) THEN
                TLMVPL=1
              ELSE
                TLMVPL=-3
                WRITE(IOUT,9001) NAMSBR,'YZ',A,B,C,D
              ENDIF
            ENDIF
          ENDIF
        ELSE IF(B .EQ. 0) THEN
          IF(C .EQ. 0) THEN
            DROITE=D/A
            IF(DROITE .LE. XMAX .AND. DROITE .GE. XMIN) THEN
              TLMVPL=1
            ELSE
              TLMVPL=-4
              WRITE(IOUT,9001) NAMSBR,'XX',A,B,C,D
            ENDIF
          ELSE
            G(1)=(D-A*XMIN)/C
            G(2)=(D-A*XMAX)/C
            G(3)=(D-C*ZMIN)/A
            IF(G(1) .LT. ZMAX .AND. G(1) .GT. ZMIN) THEN
              TLMVPL=1
            ELSE IF(G(2) .LT. ZMAX .AND. G(2) .GT. ZMIN) THEN
              TLMVPL=1
            ELSE IF(G(3) .LT. XMAX .AND. G(3) .GT. XMIN) THEN
              TLMVPL=1
            ELSE
              G(1)=(D-(A*(XMIN+XMAX)/2))/C
              G(2)=(D-(C*(ZMIN+ZMAX)/2))/A
              IF(G(1) .LE. ZMAX .AND. G(1) .GE. ZMIN) THEN
                TLMVPL=1
              ELSE IF(G(2) .LE. XMAX .AND. G(2) .GE. XMIN) THEN
                TLMVPL=1
              ELSE
                TLMVPL=-5
                WRITE(IOUT,9001) NAMSBR,'XY',A,B,C,D
              ENDIF
            ENDIF
          ENDIF
        ELSE IF(C .EQ. 0) THEN
          G(1)=(D-A*XMIN)/B
          G(2)=(D-A*XMAX)/B
          G(3)=(D-B*YMIN)/A
          IF(G(1) .LT. YMAX .AND. G(1) .GT. YMIN) THEN
            TLMVPL=1
          ELSE IF(G(2) .LT. YMAX .AND. G(2) .GT. YMIN) THEN
            TLMVPL=1
          ELSE IF(G(3) .LT. XMAX .AND. G(3) .GT. XMIN) THEN
            TLMVPL=1
          ELSE
            G(1)=(D-(A*(XMIN+XMAX)/2))/B
            G(2)=(D-(B*(YMIN+YMAX)/2))/A
            IF(G(1) .LE. YMAX .AND. G(1) .GE. YMIN) THEN
              TLMVPL=1
            ELSE IF(G(2) .LE. XMAX .AND. G(2) .GE. XMIN) THEN
              TLMVPL=1
            ELSE
              TLMVPL=-6
              WRITE(IOUT,9001) NAMSBR,'XZ',A,B,C,D
            ENDIF
          ENDIF
        ELSE
          G(1)=(D-A*XMIN-B*YMIN)/C
          G(2)=(D-A*XMIN-B*YMAX)/C
          G(3)=(D-A*XMAX-B*YMIN)/C
          G(4)=(D-A*XMAX-B*YMAX)/C
          IF(G(1) .LT. ZMAX .AND. G(1) .GT. ZMIN) THEN
            TLMVPL=1
          ELSE IF(G(2) .LT. ZMAX .AND. G(2) .GT. ZMIN) THEN
            TLMVPL=1
          ELSE IF(G(3) .LT. ZMAX .AND. G(3) .GT. ZMIN) THEN
            TLMVPL=1
          ELSE IF(G(4) .LT. YMAX .AND. G(4) .GT. YMIN) THEN
            TLMVPL=1
          ELSE
            G(1)=(D-C*ZMIN-B*YMIN)/A
            G(2)=(D-C*ZMIN-B*YMAX)/A
            G(3)=(D-C*ZMAX-B*YMIN)/A
            G(4)=(D-C*ZMAX-B*YMAX)/A
            IF(G(1) .LT. XMAX .AND. G(1) .GT. XMIN) THEN
              TLMVPL=1
            ELSE IF(G(2) .LT. XMAX .AND. G(2) .GT. XMIN) THEN
              TLMVPL=1
            ELSE IF(G(3) .LT. XMAX .AND. G(3) .GT. XMIN) THEN
              TLMVPL=1
            ELSE IF(G(4) .LT. XMAX .AND. G(4) .GT. XMIN) THEN
              TLMVPL=1
            ELSE
              G(1)=(D-C*ZMIN-A*XMIN)/B
              G(2)=(D-C*ZMIN-A*XMAX)/B
              G(3)=(D-C*ZMAX-A*XMIN)/B
              G(4)=(D-C*ZMAX-A*XMAX)/B
               IF(G(1) .LT. YMAX .AND. G(1) .GT. YMIN) THEN
                 TLMVPL=1
               ELSE IF(G(2) .LT. YMAX .AND. G(2) .GT. YMIN) THEN
                 TLMVPL=1
               ELSE IF(G(3) .LT. YMAX .AND. G(3) .GT. YMIN) THEN
                 TLMVPL=1
               ELSE IF(G(4) .LT. YMAX .AND. G(4) .GT. YMIN) THEN
                 TLMVPL=1
               ELSE
                 G(1)=(D-((C*(ZMIN+ZMAX)/2)+(B*(YMIN+YMAX)/2)))/A
                 G(2)=(D-((A*(XMIN+XMAX)/2)+(C*(ZMIN+ZMAX)/2)))/B
                 G(3)=(D-((A*(XMIN+XMAX)/2)+(B*(YMIN+YMAX)/2)))/C
                 IF(G(1) .LE. XMAX .AND. G(1) .GE. XMIN) THEN
                   TLMVPL=1
                 ELSE IF(G(2) .LE. YMAX .AND. G(2) .GE. YMIN) THEN
                   TLMVPL=1
                 ELSE IF(G(2) .LE. ZMAX .AND. G(2) .GE. ZMIN) THEN
                   TLMVPL=1
                 ELSE
                   TLMVPL=-7
                   WRITE(IOUT,9001) NAMSBR,'ZZ',A,B,C,D
                 ENDIF
               ENDIF
             ENDIF
           ENDIF
         ENDIF
       ELSE IF(NDIM .EQ. 2) THEN
*----
*  Verification que le plan
*  A*B+B*Y=D
*  croise bien la region 2D (2 PLANS)
*----
         IF(A .EQ. 0) THEN
*----
*  LIGNE PARALLELE A Y
*----
           IF(B .EQ. 0) THEN
             TLMVPL=-9
             WRITE(IOUT,9010) NAMSBR,A,B,D
           ELSE
             DROITE=D/B
             IF(DROITE .LE. YMAX .AND. DROITE .GE. YMIN) THEN
               TLMVPL=1
             ELSE
               TLMVPL=-10
               WRITE(IOUT,9011) NAMSBR,'YY',A,B,D
             ENDIF
           ENDIF
         ELSE IF(B .EQ. 0) THEN
*----
*  LIGNE PARALLELE A X
*----
           DROITE=D/A
           IF(DROITE .LE. XMAX .AND. DROITE .GE. XMIN) THEN
             TLMVPL=1
           ELSE
             TLMVPL=-11
             WRITE(IOUT,9011) NAMSBR,'XX',A,B,D
           ENDIF
         ELSE
*----
*  LIGNE DIAGONALE
*----
           G(1)=(D-A*XMIN)/B
           G(2)=(D-A*XMAX)/B
           G(3)=(D-B*YMIN)/A
           IF(G(1) .LT. YMAX .AND. G(1) .GT. YMIN) THEN
*----
*  PLAN XMIN + UN AUTRE
*----
             TLMVPL=1
           ELSE IF(G(2) .LT. YMAX .AND. G(2) .GT. YMIN) THEN
*----
*  PLAN XMAX + 1 AUTRE
*----
             TLMVPL=1
           ELSE IF(G(3) .LT. XMAX .AND. G(3) .GT. XMIN) THEN
*----
*  PLAN YMIN + 1 AUTRE
*----
             TLMVPL=1
           ELSE
*----
*  0 OU 1 PLAN
*  VERIFIER POUR COINS
*  EN TROUVANT L'INTERSECTION AVEC LE PLAN CENTRAL
*----
             G(1)=(D-A*(XMIN+XMAX)/2)/B
             G(2)=(D-(B*(YMIN+YMAX)/2))/A
             IF(G(1) .LT. YMAX .AND. G(1) .GT. YMIN) THEN
*----
*  1 COIN + PLAN CENTRAL EN X
*----
               TLMVPL=1
             ELSE IF(G(2) .LT. XMAX .AND. G(2) .GT. XMIN) THEN
*----
*  1 COIN + PLAN CENTRAL EN Y
*----
               TLMVPL=1
             ELSE
*----
*  PAS D'INTERSECTION
*----
               TLMVPL=-12
               WRITE(IOUT,9011) NAMSBR,'XY',A,B,D
             ENDIF
           ENDIF
        ENDIF
      ENDIF
*----
*  Processing finished, return
*----
      RETURN
*----
*  Output formats
*----
 9000 FORMAT(1X,'***** Warning in ',A6,' *****'/
     >       1X,'Invalid equation for 3-D plane : '/
     >       1X,F20.10,'*X + ',F20.10,'*Y + ',F20.10,'*Z = ',F20.10)
 9001 FORMAT(1X,'***** Warning in ',A6,' *****'/
     >       1X,'No intersection between region and plane in ',A2/
     >       1X,F20.10,'*X + ',F20.10,'*Y + ',F20.10,'*Z = ',F20.10)
 9010 FORMAT(1X,'***** Warning in ',A6,' *****'/
     >       1X,'Invalid equation for 2-D plane : '/
     >       1X,F20.10,'*X + ',F20.10,'*Y +  = ',F20.10)
 9011 FORMAT(1X,'***** Warning in ',A6,' *****'/
     >       1X,'No intersection between region and LINE in ',A2/
     >       1X,F20.10,'*X + ',F20.10,'*Y +  = ',F20.10)
      END
