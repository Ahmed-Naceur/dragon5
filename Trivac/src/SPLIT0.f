*DECK SPLIT0
      SUBROUTINE SPLIT0 (MAXPTS,ITYPE,NCODE,LXOLD,LYOLD,LZOLD,ISPLTX,
     1 ISPLTY,ISPLTZ,ISPLTH,ISPLTL,NMBLK,LX,LY,LZ,SIDE,XXX,YYY,ZZZ,
     2 MAT,ITYP,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Generalized mesh-splitting algorithm.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):A. Hebert
*
*Parameters: input/output
* MAXPTS  dimension of vector MAT.
* ITYPE   type of geometry.
* NCODE   boundary condition relative to each side of the domain.
* LXOLD   number of parallelepipeds along the X-axis as given in the
*         input data.
* LYOLD   number of parallelepipeds along the Y-axis.
* LZOLD   number of parallelepipeds along the Z-axis.
* ISPLTX  mesh-splitting data for parallelepipeds along the X-axis
*         negative value is used for equal-volume splitting of tubes.
* ISPLTY  mesh-splitting data for parallelepipeds along the Y-axis.
* ISPLTZ  mesh-splitting data for parallelepipeds along the Z-axis.
* ISPLTH  mesh-splitting index for hexagons into triangles.
* ISPLTL  mesh-splitting index for hexagons into lozenges.
* NMBLK   number of parallelepipeds in the domain.
* LX      number of parallelepipeds along the X-axis after mesh-
*         splitting.
* LY      number of parallelepipeds along the Y-axis.
* LZ      number of parallelepipeds along the Z-axis.
* XXX     Cartesian coordinates of the domain along the X-axis.
* YYY     Cartesian coordinates of the domain along the Y-axis.
* ZZZ     Cartesian coordinates of the domain along the Z-axis.
* MAT     index-number of the mixture type assigned to each volume
*         before and after mesh-splitting.
* ITYP    modification flag:
*         =.true. modification of XXX, YYY, ZZZ and MAT;
*         =.false. modification of MAT only.
* IMPX    print flag. Minimum printing if IMPX=0.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER MAXPTS,ITYPE,NCODE(6),LXOLD,LYOLD,LZOLD,ISPLTX(LX),
     1 ISPLTY(LY),ISPLTZ(LZ),ISPLTL,NMBLK,LX,LY,LZ,MAT(MAXPTS),IMPX
      REAL XXX(LX+1),YYY(LY+1),ZZZ(LZ+1)
      LOGICAL ITYP
*----
*  LOCAL VARIABLES
*----
      CHARACTER HSMG*130
      LOGICAL LL1,LL2,NEWCOD(6),LTRI,LLOZ
      DOUBLE PRECISION DEL,GAR
*----
*  SETTING LOGICAL PARAMETERS
*----
      LL1=(NCODE(2).EQ.3).AND.(NCODE(3).EQ.3)
      LL2=(NCODE(1).EQ.3).AND.(NCODE(4).EQ.3)
      LTRI=(ISPLTH.NE.0).AND.((ITYPE.EQ.8).OR.(ITYPE.EQ.9))
      LLOZ=(ISPLTL.NE.0).AND.((ITYPE.EQ.8).OR.(ITYPE.EQ.9))
*
      IF(ITYP) THEN
         IF(LL1.OR.LL2) THEN
*           DIAGONAL SYMMETRY: CHECK IF ISPLTY(I)=ISPLTX(I)
            DO 10 I=1,LX
            IF(ISPLTX(I).NE.ISPLTY(I)) CALL XABORT('SPLIT0: INCONSTEN'
     1      //'T MESH-SPLITTING INPUT DATA.')
   10       CONTINUE
         ENDIF
*        DETERMINATION OF THE NEW BOUNDARY CONDITIONS.
         DO 20 I=1,6
         NEWCOD(I)=.FALSE.
   20    CONTINUE
         IF((NCODE(1).EQ.5).OR.(LL2.AND.(NCODE(3).EQ.5))) THEN
            DEL=XXX(2)-XXX(1)
            IF(MOD(ISPLTX(1),2).EQ.0) THEN
               ISPLTX(1)=ISPLTX(1)/2
               NEWCOD(1)=.TRUE.
               XXX(1)=XXX(2)-REAL(0.5*DEL)
            ELSE
               IGAR=ISPLTX(1)
               ISPLTX(1)=(ISPLTX(1)+1)/2
               XXX(1)=XXX(2)-REAL(DEL*(DBLE(ISPLTX(1))/DBLE(IGAR)))
            ENDIF
         ENDIF
         IF((NCODE(2).EQ.5).OR.(LL1.AND.(NCODE(4).EQ.5))) THEN
            DEL=XXX(LX+1)-XXX(LX)
            IF(MOD(ISPLTX(LX),2).EQ.0) THEN
               ISPLTX(LX)=ISPLTX(LX)/2
               NEWCOD(2)=.TRUE.
               XXX(LX+1)=XXX(LX)+REAL(0.5*DEL)
            ELSE
               IGAR=ISPLTX(LX)
               ISPLTX(LX)=(ISPLTX(LX)+1)/2
               XXX(LX+1)=XXX(LX)+REAL(DEL*(DBLE(ISPLTX(LX))/DBLE(IGAR)))
            ENDIF
         ENDIF
         IF((NCODE(3).EQ.5).OR.(LL1.AND.(NCODE(1).EQ.5))) THEN
            DEL=YYY(2)-YYY(1)
            IF(MOD(ISPLTY(1),2).EQ.0) THEN
               ISPLTY(1)=ISPLTY(1)/2
               NEWCOD(3)=.TRUE.
               YYY(1)=YYY(2)-REAL(0.5*DEL)
            ELSE
               IGAR=ISPLTY(1)
               ISPLTY(1)=(ISPLTY(1)+1)/2
               YYY(1)=YYY(2)-REAL(DEL*(DBLE(ISPLTY(1))/DBLE(IGAR)))
            ENDIF
         ENDIF
         IF((NCODE(4).EQ.5).OR.(LL2.AND.(NCODE(2).EQ.5))) THEN
            DEL=YYY(LY+1)-YYY(LY)
            IF(MOD(ISPLTY(LY),2).EQ.0) THEN
               ISPLTY(LY)=ISPLTY(LY)/2
               NEWCOD(4)=.TRUE.
               YYY(LY+1)=YYY(LY)+REAL(0.5*DEL)
            ELSE
               IGAR=ISPLTY(LY)
               ISPLTY(LY)=(ISPLTY(LY)+1)/2
               YYY(LY+1)=YYY(LY)+REAL(DEL*(DBLE(ISPLTY(LY))/DBLE(IGAR)))
            ENDIF
         ENDIF
         IF(NCODE(5).EQ.5) THEN
            DEL=ZZZ(2)-ZZZ(1)
            IF(MOD(ISPLTZ(1),2).EQ.0) THEN
               ISPLTZ(1)=ISPLTZ(1)/2
               NEWCOD(5)=.TRUE.
               ZZZ(1)=ZZZ(2)-REAL(0.5*DEL)
            ELSE
               IGAR=ISPLTZ(1)
               ISPLTZ(1)=(ISPLTZ(1)+1)/2
               ZZZ(1)=ZZZ(2)-REAL(DEL*(DBLE(ISPLTZ(1))/DBLE(IGAR)))
            ENDIF
         ENDIF
         IF(NCODE(6).EQ.5) THEN
            DEL=ZZZ(LZ+1)-ZZZ(LZ)
            IF(MOD(ISPLTZ(LZ),2).EQ.0) THEN
               ISPLTZ(LZ)=ISPLTZ(LZ)/2
               NEWCOD(6)=.TRUE.
               ZZZ(LZ+1)=ZZZ(LZ)+REAL(0.5*DEL)
            ELSE
               IGAR=ISPLTZ(LZ)
               ISPLTZ(LZ)=(ISPLTZ(LZ)+1)/2
               ZZZ(LZ+1)=ZZZ(LZ)+REAL(DEL*(DBLE(ISPLTZ(LZ))/DBLE(IGAR)))
            ENDIF
         ENDIF
         IF((.NOT.LL2).AND.NEWCOD(1)) NCODE(1)=2
         IF((.NOT.LL1).AND.NEWCOD(2)) NCODE(2)=2
         IF((.NOT.LL1).AND.NEWCOD(3)) NCODE(3)=2
         IF((.NOT.LL2).AND.NEWCOD(4)) NCODE(4)=2
         IF(NEWCOD(5)) NCODE(5)=2
         IF(NEWCOD(6)) NCODE(6)=2
*
*        COMPUTE THE NEW VALUES OF LX, LY AND LZ.
         LXOLD=LX
         LYOLD=LY
         LZOLD=LZ
         IF(ITYPE.LT.8) THEN
            LX=0
            DO 40 IOLD=1,LXOLD
            LX=LX+ABS(ISPLTX(IOLD))
   40       CONTINUE
            LY=0
            DO 50 IOLD=1,LYOLD
            LY=LY+ISPLTY(IOLD)
   50       CONTINUE
         ELSEIF(LTRI) THEN
            LX=LXOLD*6*(ISPLTH**2)
         ELSEIF(LLOZ) THEN
            LX=LXOLD*3*(ISPLTL**2)
         ENDIF
         LZ=0
         DO 55 IOLD=1,LZOLD
         LZ=LZ+ISPLTZ(IOLD)
   55    CONTINUE
*
*        COMPUTE THE NEW VALUES OF XXX, YYY AND ZZZ.
         IF(ITYPE.LT.8) THEN
            K=LX+1
            GAR=XXX(LXOLD+1)
            DO 61 IOLD=LXOLD,1,-1
            ISP=ISPLTX(IOLD)
            DEL=(GAR-XXX(IOLD))/DBLE(ABS(ISP))
            IF(ISP.LT.0) THEN
              IF((ITYPE.EQ.3).OR.(ITYPE.EQ.6)) DEL=DEL*(GAR+XXX(IOLD))
              IF(ITYPE.EQ.4) DEL=DEL*(GAR**2+GAR*XXX(IOLD)+XXX(IOLD)**2)
            ENDIF
            GAR=XXX(IOLD)
            DO 60 I=ABS(ISP),1,-1
            IF(ISP.GT.0) THEN
               XXX(K)=REAL(GAR+DEL*DBLE(I))
            ELSE IF((ITYPE.EQ.3).OR.(ITYPE.EQ.6)) THEN
               XXX(K)=REAL(SQRT(GAR*GAR+DEL*DBLE(I)))
            ELSE IF(ITYPE.EQ.4) THEN
               XXX(K)=REAL((GAR**3+DEL*DBLE(I))**(1.0D0/3.0D0))
            ELSE
            CALL XABORT('SPLIT0: INVALID MESH-SPLITTING INDEX.')
            ENDIF
            K=K-1
   60       CONTINUE
   61       CONTINUE
            K=LY+1
            GAR=YYY(LYOLD+1)
            DO 71 IOLD=LYOLD,1,-1
            ISP=ISPLTY(IOLD)
            DEL=(GAR-YYY(IOLD))/DBLE(ISP)
            GAR=YYY(IOLD)
            DO 70 I=ISP,1,-1
            YYY(K)=REAL(GAR+DEL*DBLE(I))
            K=K-1
   70       CONTINUE
   71       CONTINUE
         ELSEIF(LTRI) THEN
            SIDE=SIDE/REAL(ISPLTH)
         ELSEIF(LLOZ) THEN
            SIDE=SIDE/REAL(ISPLTL)
         ENDIF
         K=LZ+1
         GAR=ZZZ(LZOLD+1)
         DO 76 IOLD=LZOLD,1,-1
         ISP=ISPLTZ(IOLD)
         DEL=(GAR-ZZZ(IOLD))/DBLE(ISP)
         GAR=ZZZ(IOLD)
         DO 75 I=ISP,1,-1
         ZZZ(K)=REAL(GAR+DEL*DBLE(I))
         K=K-1
   75    CONTINUE
   76    CONTINUE
*
*        COMPUTE THE NUMBER OF PARALLEPIPEDS AFTER MESH-SPLITTING.
         IF(LL1.OR.LL2) THEN
            IF(LX.EQ.LY) THEN
               NMBLK=LZ*((LX+1)*LX)/2
            ELSE
               CALL XABORT('SPLIT0: LX ET LY SHOULD BE EQUAL.')
            ENDIF
         ELSE IF(ITYPE.LT.8) THEN
            NMBLK=LX*LY*LZ
         ELSE
            NMBLK=LX*LZ
         ENDIF
         IF(IMPX.GE.3) THEN
            WRITE (6,200) LX,LY,LZ,NMBLK,(NCODE(I),I=1,6)
            IF(ITYPE.LT.8) THEN
               WRITE (6,210) 'XXX',(XXX(I),I=1,LX+1)
               WRITE (6,210) 'YYY',(YYY(I),I=1,LY+1)
            ELSE
               WRITE (6,210) 'SIDE',SIDE
            ENDIF
            WRITE (6,210) 'ZZZ',(ZZZ(I),I=1,LZ+1)
         ENDIF
      ENDIF
*----
*  COMPUTE THE NEW MIXTURE NUMBERS MAT(I).
*----
      IF(ITYPE.LT.8) THEN
         IF(LL1.OR.LL2) THEN
            KOLD=LZOLD*((LXOLD+1)*LXOLD)/2
            KNEW=LZ*((LX+1)*LX)/2
         ELSE
            KOLD=LXOLD*LYOLD*LZOLD
            KNEW=LX*LY*LZ
         ENDIF
         NMBLK=KNEW
         IF(KNEW.GT.MAXPTS) THEN
            WRITE (HSMG,230) 'NMBLK',KNEW,'MAXPTS',MAXPTS
            CALL XABORT(HSMG)
         ENDIF
         DO 103 K0=LZOLD,1,-1
         KIOFZ=KOLD
         DO 102 K=ISPLTZ(K0),1,-1
         KOLD=KIOFZ
         DO 101 K1=LYOLD,1,-1
         KIOFY=KOLD
         DO 100 J=ISPLTY(K1),1,-1
         KOLD=KIOFY
         DO 90 K2=LXOLD,1,-1
         IF(LL1.AND.(K1.LT.K2)) GO TO 90
         IF(LL2.AND.(K1.GT.K2)) GO TO 90
         IGAR=MAT(KOLD)
         DO 80 I=ABS(ISPLTX(K2)),1,-1
         IF(LL1.AND.(J.LT.I).AND.(K1.EQ.K2)) GO TO 80
         IF(LL2.AND.(J.GT.I).AND.(K1.EQ.K2)) GO TO 80
         MAT(KNEW)=IGAR
         MAT(KNEW)=IGAR
         KNEW=KNEW-1
   80    CONTINUE
         KOLD=KOLD-1
   90    CONTINUE
  100    CONTINUE
  101    CONTINUE
  102    CONTINUE
  103    CONTINUE
      ELSEIF(LTRI) THEN
*        HEXAGONAL GEOMETRY WITH TRIANGULAR SUBMESH.
         KOLD=LXOLD*LZOLD
         KNEW=LXOLD*6*(ISPLTH**2)*LZ
         NMBLK=KNEW
         IF(KNEW.GT.MAXPTS) THEN
            WRITE (HSMG,230) 'NMBLK',KNEW,'MAXPTS',MAXPTS
            CALL XABORT(HSMG)
         ENDIF
         DO 135 K0=LZOLD,1,-1
         KIOFZ=KOLD
         DO 130 K=ISPLTZ(K0),1,-1
         KOLD=KIOFZ
         DO 120 K2=LXOLD,1,-1
         IGAR=MAT(KOLD)
         DO 110 I=(6*ISPLTH**2),1,-1
         MAT(KNEW)=IGAR
         KNEW=KNEW-1
  110    CONTINUE
         KOLD=KOLD-1
  120    CONTINUE
  130    CONTINUE
  135    CONTINUE
      ELSEIF(LLOZ) THEN
*        HEXAGONAL GEOMETRY WITH LOZENGE SUBMESH.
         KOLD=LXOLD*LZOLD
         KNEW=LXOLD*3*(ISPLTL**2)*LZ
         NMBLK=KNEW
         IF(KNEW.GT.MAXPTS) THEN
            WRITE (HSMG,230) 'NMBLK',KNEW,'MAXPTS',MAXPTS
            CALL XABORT(HSMG)
         ENDIF
         DO 165 K0=LZOLD,1,-1
         KIOFZ=KOLD
         DO 160 K=ISPLTZ(K0),1,-1
         KOLD=KIOFZ
         DO 150 K2=LXOLD,1,-1
         IGAR=MAT(KOLD)
         DO 140 I=(3*ISPLTL**2),1,-1
         MAT(KNEW)=IGAR
         KNEW=KNEW-1
  140    CONTINUE
         KOLD=KOLD-1
  150    CONTINUE
  160    CONTINUE
  165    CONTINUE
      ELSE
*        HEXAGONAL GEOMETRY.
         KOLD=LXOLD*LZOLD
         KNEW=LXOLD*LZ
         NMBLK=KNEW
         IF(KNEW.GT.MAXPTS) THEN
            WRITE (HSMG,230) 'NMBLK',KNEW,'MAXPTS',MAXPTS
            CALL XABORT(HSMG)
         ENDIF
         DO 185 K0=LZOLD,1,-1
         KIOFZ=KOLD
         DO 180 K=ISPLTZ(K0),1,-1
         KOLD=KIOFZ
         DO 170 K2=LXOLD,1,-1
         MAT(KNEW)=MAT(KOLD)
         KNEW=KNEW-1
         KOLD=KOLD-1
  170    CONTINUE
  180    CONTINUE
  185    CONTINUE
      ENDIF
      IF(IMPX.GE.3) WRITE (6,220) (MAT(I),I=1,NMBLK)
      RETURN
*
  200 FORMAT (//4H LX=,I4,4X,3HLY=,I4,4X,3HLZ=,I4,4X,6HNMBLK=,I5,
     1 4X,9HNCODE(1)=,I2,3X,9HNCODE(2)=,I2,3X,9HNCODE(3)=,I2,3X,
     2 9HNCODE(4)=,I2,3X,9HNCODE(5)=,I2,3X,9HNCODE(6)=,I2/)
  210 FORMAT (//1X,A4/(1X,1P,10E12.4))
  220 FORMAT (//4H MAT/(1X,20I6))
  230 FORMAT (29HSPLIT0: INSUFFICIENT STORAGE.,5X,A6,1H=,I9,8H ; AVAIL,
     1 13HABLE STORAGE ,A6,1H=,I9)
      END
