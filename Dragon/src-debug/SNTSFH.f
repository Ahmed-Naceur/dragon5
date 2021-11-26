*DECK SNTSFH
      SUBROUTINE SNTSFH (IMPX,LXH,ISPLH,MAT,IZGLOB,ILOZSWP,CONNEC,
     1 CONFROM)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Output matrices/arrays/graphs needed for resolution of the discrete
* ordinates transport equation in hexagonal geometry, including the 
* hexagon and lozenge sweep orders which both depend on the direction
* of sweep.
*
*Copyright:
* Copyright (C) 2019 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. A. Calloo
*
*Parameters: input
* IMPX    print parameter.
* LXH     number of hexagons.
* ISPLH   type of hexagonal mesh-splitting:
*         =1 no mesh splitting (complete hexagons);
*         =K 6*(K-1)*(K-1) triangles per hexagon.
* MAT     mixture index assigned to each element.
*
*Parameters: output
* IZGLOB  hexagon sweep order depending on direction.
* ILOZSWP lozenge sweep order depending on direction.
* CONNEC  connectivity matrix for flux swapping -- which lozenges is the
*         lozenge under consideration connected to; in order to pass the
*         flux along. This is dependent on direction.
* CONFROM matrix for incoming flux -- which lozenges are feeding into
*         the lozenge under consideration. This is dependent on
*         direction.
*
*Comments:
* The lozenge under consideration is given by the position within the
* the matrix. See user manual and/or data manual and/or thesis
*                                 _____
*                                /   / \
*                               / B /   \
*                         ,----(----  A  )----.
*                        /      \ C \   /      \
*                       /        \___\_/        \
*                       \   4    /     \   2    /
*                        \      /       \      /
*                         )----(   1     )----(
*                        /      \       /      \
*                       /        \_____/        \
*                       \   5    /     \   7    /
*                        \      /       \      /
*                         `----(   6     )----'
*                               \       /
*                                \_____/
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IMPX,LXH,ISPLH,MAT(3,ISPLH**2,LXH),IZGLOB(LXH,6)
      INTEGER CONNEC(3,(LXH*3)*2,6)
      INTEGER, DIMENSION(2,3,6) :: CONFROM
*----
*  LOCAL VARIABLES
*----
      INTEGER, DIMENSION(:),ALLOCATABLE :: I1,I3,I4
      INTEGER, DIMENSION(3,6) :: ILOZSWP
      INTEGER, DIMENSION(2,3) :: IFROM
*----
*  DISCRETE ORDINATES SWEEP SPECIFIC NUMEROTATION
*----
*
*--------
*  Numbering of order of sweep of lozenges in each hexagon depending
*  on direction
*--------
      ILOZSWP = RESHAPE((/ 3, 2, 1, 3, 1, 2, 1, 3, 2, 1, 2, 3, 2, 1,
     1   3, 2, 3, 1 /), SHAPE(ILOZSWP))
      DO I=1,3
         WRITE(*,*) ILOZSWP(I,:)
      ENDDO
*--------
*  Numbering of order of hexagons depending on direction
*--------
      NBC=INT((SQRT(  REAL((4*LXH-1)/3)  )+1.)/2.)
      NCOL=2*NBC -1
      NTOTHEX = (((NBC-1)*(NBC)/2)*6)+1
      IF(NBC.EQ.1) CALL XABORT('does not work for single hexagon yet')
      IF(LXH.NE.1+3*NBC*(NBC-1)) CALL XABORT('SNTSFH: INVALID VALUE OF'
     1 //' LXH(1).')
*
      ALLOCATE(I1(6*LXH),I3(LXH),I4(LXH))
*
      LT3=0
      LT4=0
      DO 30 I=1,LXH
         I3(I)=I
         I4(I)=0
         LT3=LT3+1
         IF(MAT(1,1,I).GT.0) THEN
            LT4=LT4+1
            I4(I)=I
         ENDIF
   30 CONTINUE

      CALL XDISET(IZGLOB,6*LXH,0)

      DO IND=1,6
         J1 = LXH - ((NBC-2)*3 +2) + (NBC-1)*(IND-1)
         IF(J1.GT.LXH) J1=J1 - ((NBC-1)*6)
         CALL SNTPER(J1,IND,LXH,LT3,I1((IND-1)*LXH +1),I3)

         DO I=1,LXH
            IOFIDIR = I1((IND-1)*LXH +I)
            IZGLOB(IOFIDIR,IND) = I3(I)
         ENDDO
      ENDDO
*--------
*  Numbering of connectivity of lozenges depending on directions
*--------
! I2    - hexagon currently under consideration
! J2    - lozenge currently under consideration
! IHEX1 - hexagon connected to side X of J2 lozenge
! IHEX2 - hexagon connected to side Y of J2 lozenge
! ILOZ1 - lozenge connected to side X of J2 lozenge
! ILOZ2 - lozenge connected to side Y of J2 lozenge
      IFACE1 = 0
      IFACE2 = 0
      IFACE3 = 0
      DO IND=1,6
         IF(IND.EQ.1) THEN
            IFACE1 = 6
            IFACE2 = 2
            IFACE3 = 1
         ELSE IF(IND.EQ.2) THEN
            IFACE1 = 1
            IFACE2 = 3
            IFACE3 = 2
         ELSE IF(IND.EQ.3) THEN
            IFACE1 = 2
            IFACE2 = 4
            IFACE3 = 3
         ELSE IF(IND.EQ.4) THEN
            IFACE1 = 3
            IFACE2 = 5
            IFACE3 = 4
         ELSE IF(IND.EQ.5) THEN
            IFACE1 = 4
            IFACE2 = 6
            IFACE3 = 5
         ELSE IF(IND.EQ.6) THEN
            IFACE1 = 5
            IFACE2 = 1
            IFACE3 = 6
         ELSE
            CALL XABORT('SNTSFH: INVALID DATA')
         ENDIF

         IHEX1 = 0
         IHEX2 = 0
         ISIDE1 = 0
         ISIDE2 = 0
         DO I=1,LXH
            I2=IZGLOB(I,IND)
            DO J=1,3
               J2=ILOZSWP(J,IND)
               JINDEX = ((I2-1)*3*2) +((J2-1)*2) +1
*
               IF(J.EQ.1)THEN
                  IHEX1 = I2
                  IHEX2 = I2
               ELSEIF(J.EQ.2)THEN
                  IF((IND.EQ.1).OR.(IND.EQ.3).OR.(IND.EQ.5))THEN
                     IFACE=IFACE2
                  ELSE
                     IFACE=IFACE1
                  ENDIF
                  IF((IND.EQ.1).OR.(IND.EQ.2).OR.(IND.EQ.6))THEN
                     IHEX1 = I2
                     IHEX2 = NEIGHB(I2,IFACE,9,LXH,POIDS)
                  ELSE
                     IHEX1 = NEIGHB(I2,IFACE,9,LXH,POIDS)
                     IHEX2 = I2
                  ENDIF
               ELSEIF(J.EQ.3)THEN
                  IF((IND.EQ.1).OR.(IND.EQ.5))THEN
                     IFACEi=IFACE1
                     IFACEj=IFACE3
                  ELSEIF((IND.EQ.4).OR.(IND.EQ.6))THEN
                     IFACEi=IFACE3
                     IFACEj=IFACE2
                  ELSEIF((IND.EQ.2))THEN
                     IFACEi=IFACE2
                     IFACEj=IFACE3
                  ELSEIF((IND.EQ.3))THEN
                     IFACEi=IFACE3
                     IFACEj=IFACE1
                  ENDIF
                  IHEX1 = NEIGHB(I2,IFACEi,9,LXH,POIDS)
                  IHEX2 = NEIGHB(I2,IFACEj,9,LXH,POIDS)
               ENDIF
*
               ILOZ1=J2-2
               ILOZ2=J2-1
               IF(ILOZ1.LT.1) ILOZ1=ILOZ1+3
               IF(ILOZ2.LT.1) ILOZ2=ILOZ2+3
               IF(J2.EQ.2)THEN
                  ILOZTEMP = ILOZ1
                  ILOZ1 = ILOZ2
                  ILOZ2 = ILOZTEMP
               ENDIF
*
               IF((J2.EQ.3))THEN
                  ISIDE1=2
                  ISIDE2=2
               ELSEIF(J2.EQ.1)THEN
                  ISIDE1=1
                  ISIDE2=1
               ELSEIF(J2.EQ.2)THEN
                  ISIDE1=1
                  ISIDE2=2
               ENDIF
*
               CONNEC(1,JINDEX,IND)   = IHEX1
               CONNEC(2,JINDEX,IND)   = ILOZ1
               CONNEC(3,JINDEX,IND)   = ISIDE1
               CONNEC(1,JINDEX+1,IND) = IHEX2
               CONNEC(2,JINDEX+1,IND) = ILOZ2
               CONNEC(3,JINDEX+1,IND) = ISIDE2
*
            ENDDO
         ENDDO
      ENDDO

      IFROM(1,:)=0
      IFROM(2,:)=1

      DO IND=1,6

         IFROM(1,:)=IFROM(1,:) +1
         IFROM(2,:)=IFROM(2,:) +1

         DO J=1,3
            IF(J.EQ.1)THEN
               IF((IND.EQ.4))THEN
                  IFROM(1,J) = IFROM(1,J) -1
                  IFROM(2,J) = IFROM(2,J) -1
               ENDIF
            ELSEIF(J.EQ.2)THEN
               IF((IND.EQ.3).OR.(IND.EQ.6))THEN
                  IFROM(1,J) = IFROM(1,J) -1
                  IFROM(2,J) = IFROM(2,J) -1
               ENDIF
            ELSEIF(J.EQ.3)THEN
               IF((IND.EQ.2).OR.(IND.EQ.5))THEN
                  IFROM(1,J) = IFROM(1,J) -1
                  IFROM(2,J) = IFROM(2,J) -1
               ENDIF
            ENDIF
            IF(IFROM(1,J).GT.4) IFROM(1,J)=IFROM(1,J)-4
            IF(IFROM(2,J).GT.4) IFROM(2,J)=IFROM(2,J)-4

            IF(MOD(IFROM(1,J),2).EQ.1)THEN
               CONFROM(1,J,IND) = IFROM(1,J)
               CONFROM(2,J,IND) = IFROM(2,J)
            ELSEIF(MOD(IFROM(1,J),2).EQ.0)THEN
               CONFROM(1,J,IND) = IFROM(2,J)
               CONFROM(2,J,IND) = IFROM(1,J)
            ENDIF
         ENDDO
      ENDDO
*----
*  PRINT A FEW GEOMETRY CHARACTERISTICS
*----
      IF(IMPX.GT.0) THEN
         write(6,*) ' '
         write(6,*) 'NBC   =',NBC
      ENDIF
      RETURN
      END
