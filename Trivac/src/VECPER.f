*DECK VECPER
      SUBROUTINE VECPER(HNAME,IMPV,ISEG,L4,MUIN,LON,LTSW,NBL,LBL,MUV,
     1 IPV)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Ordering of matrix elements for supervectorial operations on a matrix
* in compressed diagonal storage mode.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* HNAME   name of the matrix (for edition purpose only).
* IMPV    print parameter for statistics (equal to zero for no print).
* ISEG    number of elements in a vector register.
* L4      matrix order.
* MUIN    position of each diagonal element in non-ordered matrix.
*
*Parameters: output
* LON    number of groups of linear systems.
* LTSW   maximum bandwidth (=2 for tridiagonal systems).
* NBL    number of linear systems in each group.
* LBL    number of unknowns in each group.
* MUV    position of each diagonal element in ordered matrix.
* IPV    permutation vector for the ordered unknowns.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER HNAME*4
      INTEGER IMPV,ISEG,L4,MUIN(L4),LON,LTSW,NBL(LON),LBL(LON),MUV(L4),
     1 IPV(L4)
*----
*  LOCAL VARIABLES
*----
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISET,ISORT,IOFSET,IORD
*----
*  DETERMINE THE TOTAL NUMBER OF LINEAR SYSTEMS AND COMPUTE THE ORDER
*  OF EACH OF THEM
*----
      ALLOCATE(ISET(L4))
      ISET(1)=0
      K1=MUIN(1)+1
      DO 10 I=2,L4
      ISET(I)=0
      K2=MUIN(I)
      DO 5 J=I-K2+K1,I-1
      ISET(J)=1
    5 CONTINUE
      K1=K2+1
   10 CONTINUE
      NSYS=0
      DO 15 I=1,L4
      IPV(I)=0
      MUV(I)=0
      IF(ISET(I).EQ.0) NSYS=NSYS+1
   15 CONTINUE
      LON=1+(NSYS-1)/ISEG
      IF(IMPV.GE.2) WRITE (6,'(/35H VECPER: NUMBER OF INDEPENDANT LINE,
     1 22HAR SYSTEMS IN MATRIX '',A4,3H'' =,I7/9X,17HNUMBER OF GROUPS ,
     1 19HOF LINEAR SYSTEMS =,I6)') HNAME,NSYS,LON
      ALLOCATE(IORD(NSYS),IOFSET(NSYS))
      ISYS=0
      IORD0=0
      IOFSET(1)=1
      DO 20 I=1,L4
      IF(ISET(I).EQ.0) THEN
         ISYS=ISYS+1
         IORD(ISYS)=I-IORD0
         IF(I.NE.L4) IOFSET(ISYS+1)=I+1
         IORD0=I
      ENDIF
   20 CONTINUE
      DEALLOCATE(ISET)
*----
*  SORT THE LINEAR SYSTEMS BY DECREASING ORDER
*----
      ALLOCATE(ISORT(NSYS))
      JNEW=NSYS
      DO 25 ISYS=NSYS,1,-1
      IF(IORD(ISYS).EQ.1) THEN
         ISORT(JNEW)=ISYS
         JNEW=JNEW-1
      ENDIF
   25 CONTINUE
      INEW=0
   30 IBIG=0
      DO 50 ISYS=1,NSYS
      IF(IORD(ISYS).EQ.1) GO TO 50
      DO 40 KSYS=1,INEW
      IF(ISORT(KSYS).EQ.ISYS) GO TO 50
   40 CONTINUE
      IBIG=MAX(IBIG,IORD(ISYS))
   50 CONTINUE
      IF(IBIG.LE.1) GO TO 70
      DO 60 ISYS=1,NSYS
      IF(IORD(ISYS).EQ.IBIG) THEN
         INEW=INEW+1
         ISORT(INEW)=ISYS
      ENDIF
   60 CONTINUE
      GO TO 30
   70 IF(INEW.NE.JNEW) CALL XABORT('VECPER: ALGORITHM FAILURE 1')
      DO 80 I=1,LON
      ISYS=ISORT((I-1)*ISEG+1)
      NBL(I)=ISEG
      LBL(I)=IORD(ISYS)
   80 CONTINUE
      NBL(LON)=NSYS-(LON-1)*ISEG
      IF(IMPV.GE.2) WRITE (6,'(9X,33HMAXIMUM ORDER OF AN INDEPENDANT L,
     1 14HINEAR SYSTEM =,I9)') LBL(1)
      IF(IMPV.GE.3) THEN
         I1=1
         DO 90 I=1,(LON-1)/8+1
         I2=I1+7
         IF(I2.GT.LON) I2=LON
         WRITE (6,200) (J,J=I1,I2)
         WRITE (6,210) (NBL(J),J=I1,I2)
         WRITE (6,220) (LBL(J),J=I1,I2)
         I1=I1+8
   90    CONTINUE
      ENDIF
*----
*  COMPUTE THE PERMUTATION MATRIX
*----
      LBL0=0
      KSYS=0
      DO 105 J=1,LON
      DO 101 K=1,NBL(J)
      KSYS=KSYS+1
      ISYS=ISORT(KSYS)
      IOF0=IOFSET(ISYS)
      IOF1=IOF0+IORD(ISYS)-1
      IF(IOF1.GT.L4) CALL XABORT('VECPER: ALGORITHM FAILURE 2')
      DO 100 I=IOF0,IOF1
      IPV(I)=(LBL0+I-IOF0)*ISEG+K
  100 CONTINUE
  101 CONTINUE
      LBL0=LBL0+LBL(J)
  105 CONTINUE
      DO 110 I=1,L4
      IF(IPV(I).LE.0) CALL XABORT('VECPER: ALGORITHM FAILURE 3')
      IF(IPV(I).GT.LBL0*ISEG) CALL XABORT('VECPER: ALGORITHM FAILURE 4')
  110 CONTINUE
      L4NEW=0
      DO 115 J=1,LON
      L4NEW=L4NEW+LBL(J)*NBL(J)
  115 CONTINUE
      IF(IMPV.GE.2) WRITE (6,'(/35H VECPER: INCREASING NUMBER OF UNKNO,
     1 8HWNS FROM,I7,3H TO,I7,11H. FILL-IN =,F7.2,3H %.)') L4,L4NEW,
     2 100.0*(REAL(L4NEW)/REAL(L4)-1.0)
*----
*  COMPUTE THE VECTORIAL BANDWIDTH
*----
      LBL0=0
      KSYS=0
      IIMAX=0
      LTSW=0
      MAXNEW=0
      MUVOLD=0
      DO 150 J=1,LON
      DO 120 I=1,LBL(J)
      MUV(LBL0+I)=1
  120 CONTINUE
      DO 131 K=1,NBL(J)
      KSYS=KSYS+1
      ISYS=ISORT(KSYS)
      IOF0=IOFSET(ISYS)-1
      DO 130 I=2,IORD(ISYS)
      IBIG=MUIN(IOF0+I)-MUIN(IOF0+I-1)
      IF(IBIG.GT.MUV(LBL0+I)) MUV(LBL0+I)=IBIG
  130 CONTINUE
  131 CONTINUE
      DO 140 I=1,LBL(J)
      LTSW=MAX(LTSW,MUV(LBL0+I))
      IIMAX=IIMAX+MUV(LBL0+I)
      MUV(LBL0+I)=IIMAX
  140 CONTINUE
      LBL0=LBL0+LBL(J)
      MAXNEW=MAXNEW+(MUV(LBL0)-MUVOLD)*NBL(J)
      MUVOLD=MUV(LBL0)
  150 CONTINUE
      IF(IMPV.GE.2) WRITE (6,'(/35H VECPER: INCREASING NUMBER OF TERMS,
     1 17H IN MATRICES FROM,I9,3H TO,I9,11H. FILL-IN =,F7.2,3H %./9X,
     2 19HMAXIMUM BANDWIDTH =,I4)') MUIN(L4),MAXNEW,
     3 100.0*(REAL(MAXNEW)/REAL(MUIN(L4))-1.0),LTSW
*
      DEALLOCATE(ISORT,IOFSET,IORD)
      RETURN
*
  200 FORMAT (//13H GROUP       ,8(I8,5X,1HI))
  210 FORMAT (  13H NB. SYSTEMS ,8(I8,5X,1HI))
  220 FORMAT (  13H NB. UNKNOWNS,8(I8,5X,1HI))
      END
