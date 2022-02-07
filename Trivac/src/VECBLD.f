*DECK VECBLD
      SUBROUTINE VECBLD(ISEG,L4,MUIN,LON,LBL,MUV,IPV,ITY,ASSIN,ASSV,
     1 DGV)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Rebuild a matrix stored in compressed diagonal storage mode in a form
* compatible with supervectorial calculations.
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
* ISEG    number of elements in a vector register.
* L4      ASSIN matrix order.
* MUIN    position of each diagonal element in matrix ASSIN.
* LON     number of groups of linear systems.
* LBL     number of unknowns in each group.
* MUV     position of each diagonal element in matrix ASSV.
* IPV     permutation vector for the ordered unknowns.
* ITY     type of operation: =1 gather back; =2 scatter forth.
*
*Parameters: input/output
* ASSIN   input (ITY=2) or output (ITY=1) matrix in scalar compressed
*         diagonal storage mode. Dimensionned to MUIN(L4).
* ASSV    input (ITY=1) or output (ITY=2) matrix in supervectorial
*         compressed diagonal storage mode. The second dimension is
*         equal to MUV(SUM(LBL(I))).
*
*Parameters: output
* DGV     diagonal of ASSV. This information is produced only if ITY=2.
*         The second dimension is equal to SUM(LBL(I)).
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER ISEG,L4,MUIN(L4),LON,LBL(LON),MUV(L4),IPV(L4),ITY
      REAL ASSIN(*),ASSV(ISEG,*),DGV(ISEG,*)
*----
*  REBUILD THE MATRIX
*----
      IF(ITY.EQ.1) THEN
         IOF0=0
         DO 20 IND1=1,L4
         IPOS=1+(IPV(IND1)-1)/ISEG
         IBANC=1+MOD(IPV(IND1)-1,ISEG)
         IOF1=MUIN(IND1)
         DO 10 JND1=1-IOF1+IOF0,0
         ASSIN(IOF1+JND1)=ASSV(IBANC,MUV(IPOS)+JND1)
   10    CONTINUE
         IOF0=IOF1
   20    CONTINUE
      ELSE IF(ITY.EQ.2) THEN
         LBL0=0
         DO 30 I=1,LON
         LBL0=LBL0+LBL(I)
   30    CONTINUE
         DO 45 J=1,MUV(LBL0)
         DO 40 I=1,ISEG
         ASSV(I,J)=0.0
   40    CONTINUE
   45    CONTINUE
         LBL0=0
         DO 60 J=1,LON
         DO 55 K=1,LBL(J)
         DO 50 I=1,ISEG
         ASSV(I,MUV(LBL0+K))=1.0
         DGV(I,LBL0+K)=1.0
   50    CONTINUE
   55    CONTINUE
         LBL0=LBL0+LBL(J)
   60    CONTINUE
         IOF0=0
         DO 80 IND1=1,L4
         IPOS=1+(IPV(IND1)-1)/ISEG
         IBANC=1+MOD(IPV(IND1)-1,ISEG)
         IOF1=MUIN(IND1)
         DO 70 JND1=1-IOF1+IOF0,0
         ASSV(IBANC,MUV(IPOS)+JND1)=ASSIN(IOF1+JND1)
   70    CONTINUE
         DGV(IBANC,IPOS)=ASSIN(IOF1)
         IOF0=IOF1
   80    CONTINUE
      ENDIF
      RETURN
      END
