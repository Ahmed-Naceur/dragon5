*DECK NUMER3
      SUBROUTINE NUMER3 (NCOUR,MULTC,NCODE,ZCODE,LX,LY,LZ,IORI,ISM,
     1 POURCE,IMPX,NMBLK,IFR,ALB,SUR,NMERGE,INUM,MIX,DVX,NGEN,IGEN,
     2 XX,YY,ZZ)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Surface renumbering algorithm for Cartesian geometry.
* The 3-D DP-1 approximation is not implemented.
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
*Parameters: input/output
* NCOUR   number of surfaces per block (input); number of out-currents
*         per block (output).
* MULTC   type of multicell approximation:
*         =1 Roth;    =2 Roth X ncour;    =3 DP-0;    =4 DP-1.
* NCODE   type of boundary condition on each side of the domain:
*         =0 not used;   =1 VOID;     =2 REFL;
*         =3 DIAG;       =4 TRAN;     =5 SYME.
* ZCODE   value of the albedo on each side of the domain.
* LX      number of blocks along the X-axis.
* LY      number of blocks along the Y-axis.
* LZ      number of blocks along the Z-axis.
* IORI    orientation of the blocks.
* ISM     permutation index corresponding to each orientation
*         (ISM(I,N)=I is the natural orientation).
* POURCE  weight associated with each merged block.
* IMPX    print flag (equal to 0 for no print).
* NMBLK   total number of blocks in the domain.
* IFR     index-number of in-currents.
* ALB     transmission/albedo associated with each in-current.
* SUR     surface associated with each in-current.
* NMERGE  total number of merged cells for which specific values
*         of the neutron flux and reactions rates are required.
*         Many cells with different position in the domain can
*         be merged before the neutron flux calculation if they
*         own the same generating cell (NMERGE.le.NMBLK).
* INUM    index-number of the merged cell associated to each cell.
* MIX     index-number of out-currents.
* DVX     weight associated with each out-current.
*         Note: IFR, ALB, MIX and DVX contains information to rebuild
*         the geometrical 'A' matrix.
* NGEN    total number of generating blocks in the cartesian domain.
* IGEN    index-number of the generating block associated with each
*         merged block.
* XX      X-thickness of the generating blocks.
* YY      Y-thickness of the generating blocks.
* ZZ      Z-thickness of the generating blocks.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NCOUR,MULTC,NCODE(6),LX,LY,LZ,IORI(NMBLK),ISM(6,8),
     1 IMPX,NMBLK,IFR(12*NMBLK),NMERGE,INUM(NMBLK),MIX(12*NMERGE),
     2 NGEN,IGEN(NMERGE)
      REAL ZCODE(6),POURCE(NMERGE),ALB(12*NMBLK),SUR(12*NMBLK),
     1 DVX(12*NMERGE),XX(NGEN),YY(NGEN),ZZ(NGEN)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (EPS=1.0E-5)
      LOGICAL LL1,LL2,LOG1,LOG2,LOG3
      CHARACTER DIRR(6)*2,DIRZ(12)*2,HSMG*131
      INTEGER IDDD(6),ISMZ(12)
      REAL DDD(6)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: JF2
      REAL, ALLOCATABLE, DIMENSION(:) :: GG3
      SAVE DIRR
      DATA DIRR/'X-','X+','Y-','Y+','Z-','Z+'/
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(JF2(6*NMBLK),GG3(NMBLK))
*
      DO 100 I=1,NCOUR*NMERGE
      MIX(I)=I
      DVX(I)=1.0
100   CONTINUE
      IS1=0
      IS2=0
      LXY=LX*LY
      LL1=((NCODE(2).EQ.3).AND.(NCODE(3).EQ.3))
      LL2=((NCODE(1).EQ.3).AND.(NCODE(4).EQ.3))
      IF (LL1) THEN
         IS1=1
         LXY=LX*(LX+1)/2
      ELSE IF (LL2) THEN
         IS2=1
         LXY=LX*(LX+1)/2
      ENDIF
      IBLK=0
      DO 280 K0=1,LZ
      DO 275 K1=1,LY
      LXM=1
      LXP=LX
      IF (LL1) LXP=K1
      IF (LL2) LXM=K1
      DO 270 K2=LXM,LXP
      IBLK=IBLK+1
      IKK=INUM(IBLK)
      FRX=1.0
      FRY=1.0
      FRZ=1.0
      IF (IKK.EQ.0) GO TO 265
      IS=NCOUR*(IBLK-1)
      IT=NCOUR*(IKK-1)
      II=IORI(IBLK)
      DO 110 IC=1,6
      IDDD(IC)=-1
110   CONTINUE
      IF (K2.GT.1)  IDDD(1)=IBLK-1
      IF (K2.LT.LX) IDDD(2)=IBLK+1
      IF (K1.GT.1)  IDDD(3)=IBLK-(LXP-LXM+1)+IS1
      IF (K1.LT.LY) IDDD(4)=IBLK+(LXP-LXM+1)-IS2
      IF (K0.GT.1)  IDDD(5)=IBLK-LXY
      IF (K0.LT.LZ) IDDD(6)=IBLK+LXY
*
      DO 120 IC=1,NCOUR
      ALB(IS+IC)=1.0
      SUR(IS+IC)=0.0
      JBLK=IDDD(IC)
      IF (JBLK.GT.0) THEN
         JKK=INUM(JBLK)
         JT=NCOUR*(JKK-1)
         IF ((MOD(IC,2).EQ.1).AND.(JKK.GT.0)) THEN
            IFR(IS+ISM(IC,II))=JT+ISM(IC+1,IORI(JBLK))
         ELSE IF ((MOD(IC,2).EQ.0).AND.(JKK.GT.0)) THEN
            IFR(IS+ISM(IC,II))=JT+ISM(IC-1,IORI(JBLK))
         ELSE
            IFR(IS+ISM(IC,II))=0
         ENDIF
         IDDD(IC)=JKK
      ELSE
         IFR(IS+ISM(IC,II))=0
      ENDIF
120   CONTINUE
*----
*  VOID OR REFL BOUNDARY CONDITIONS
*----
      IKG=IGEN(IKK)
      LOG1=(K2.EQ.1).OR.(IDDD(1).EQ.0)
      LOG2=(NCODE(1).EQ.1).OR.(LL2.AND.(NCODE(3).EQ.1))
      LOG3=(NCODE(1).EQ.2).OR.(LL2.AND.(NCODE(3).EQ.2))
      IF (LOG1.AND.LOG2) THEN
         ALB(IS+ISM(1,II))=-ZCODE(1)
         IFR(IS+ISM(1,II))=IT+ISM(1,II)
      ELSE IF (LOG1.AND.LOG3) THEN
         ALB(IS+ISM(1,II))=-1.0
         IFR(IS+ISM(1,II))=IT+ISM(1,II)
      ENDIF
      IF(LOG1.AND.(NCODE(1).EQ.1)) SUR(IS+ISM(1,II))=YY(IKG)*ZZ(IKG)
      IF(LOG1.AND.(NCODE(1).EQ.2)) SUR(IS+ISM(1,II))=YY(IKG)*ZZ(IKG)
      LOG1=(K2.EQ.LX).OR.(IDDD(2).EQ.0)
      LOG2=(NCODE(2).EQ.1).OR.(LL1.AND.(NCODE(4).EQ.1))
      LOG3=(NCODE(2).EQ.2).OR.(LL1.AND.(NCODE(4).EQ.2))
      IF (LOG1.AND.LOG2) THEN
         ALB(IS+ISM(2,II))=-ZCODE(2)
         IFR(IS+ISM(2,II))=IT+ISM(2,II)
      ELSE IF (LOG1.AND.LOG3) THEN
         ALB(IS+ISM(2,II))=-1.0
         IFR(IS+ISM(2,II))=IT+ISM(2,II)
      ENDIF
      IF(LOG1.AND.(NCODE(2).EQ.1)) SUR(IS+ISM(2,II))=YY(IKG)*ZZ(IKG)
      IF(LOG1.AND.(NCODE(2).EQ.2)) SUR(IS+ISM(2,II))=YY(IKG)*ZZ(IKG)
      LOG1=(K1.EQ.1).OR.(IDDD(3).EQ.0)
      LOG2=(NCODE(3).EQ.1).OR.(LL1.AND.(NCODE(1).EQ.1))
      LOG3=(NCODE(3).EQ.2).OR.(LL1.AND.(NCODE(1).EQ.2))
      IF (LOG1.AND.LOG2) THEN
         ALB(IS+ISM(3,II))=-ZCODE(3)
         IFR(IS+ISM(3,II))=IT+ISM(3,II)
      ELSE IF (LOG1.AND.LOG3) THEN
         ALB(IS+ISM(3,II))=-1.0
         IFR(IS+ISM(3,II))=IT+ISM(3,II)
      ENDIF
      IF(LOG1.AND.(NCODE(3).EQ.1)) SUR(IS+ISM(3,II))=XX(IKG)*ZZ(IKG)
      IF(LOG1.AND.(NCODE(3).EQ.2)) SUR(IS+ISM(3,II))=XX(IKG)*ZZ(IKG)
      LOG1=(K1.EQ.LY).OR.(IDDD(4).EQ.0)
      LOG2=(NCODE(4).EQ.1).OR.(LL2.AND.(NCODE(2).EQ.1))
      LOG3=(NCODE(4).EQ.2).OR.(LL2.AND.(NCODE(2).EQ.2))
      IF (LOG1.AND.LOG2) THEN
         ALB(IS+ISM(4,II))=-ZCODE(4)
         IFR(IS+ISM(4,II))=IT+ISM(4,II)
      ELSE IF (LOG1.AND.LOG3) THEN
         ALB(IS+ISM(4,II))=-1.0
         IFR(IS+ISM(4,II))=IT+ISM(4,II)
      ENDIF
      IF(LOG1.AND.(NCODE(4).EQ.1)) SUR(IS+ISM(4,II))=XX(IKG)*ZZ(IKG)
      IF(LOG1.AND.(NCODE(4).EQ.2)) SUR(IS+ISM(4,II))=XX(IKG)*ZZ(IKG)
      LOG1=(K0.EQ.1).OR.(IDDD(5).EQ.0)
      IF (LOG1.AND.(NCODE(5).EQ.1)) THEN
         ALB(IS+ISM(5,II))=-ZCODE(5)
         IFR(IS+ISM(5,II))=IT+ISM(5,II)
      ELSE IF (LOG1.AND.(NCODE(5).EQ.2)) THEN
         ALB(IS+ISM(5,II))=-1.0
         IFR(IS+ISM(5,II))=IT+ISM(5,II)
      ENDIF
      IF(LOG1.AND.(NCODE(5).EQ.1)) SUR(IS+ISM(5,II))=XX(IKG)*YY(IKG)
      IF(LOG1.AND.(NCODE(5).EQ.2)) SUR(IS+ISM(5,II))=XX(IKG)*YY(IKG)
      LOG1=(K0.EQ.LZ).OR.(IDDD(6).EQ.0)
      IF (LOG1.AND.(NCODE(6).EQ.1)) THEN
         ALB(IS+ISM(6,II))=-ZCODE(6)
         IFR(IS+ISM(6,II))=IT+ISM(6,II)
      ELSE IF (LOG1.AND.(NCODE(6).EQ.2)) THEN
         ALB(IS+ISM(6,II))=-1.0
         IFR(IS+ISM(6,II))=IT+ISM(6,II)
      ENDIF
      IF(LOG1.AND.(NCODE(6).EQ.1)) SUR(IS+ISM(6,II))=XX(IKG)*YY(IKG)
      IF(LOG1.AND.(NCODE(6).EQ.2)) SUR(IS+ISM(6,II))=XX(IKG)*YY(IKG)
*----
*  CORRECT THE PARITY OF THE INTERFACE CURRENTS FOR DP-1 CASES WITH
*  'MIRROR' ORIENTATION
*----
      DO 125 IC=1,NCOUR
      IF(II.GE.5) ALB(IS+IC)=-ALB(IS+IC)
      JBLK=IBLK
      IF((K2.GT.1).AND.(IC.EQ.1))  JBLK=IBLK-1
      IF((K2.LT.LX).AND.(IC.EQ.2)) JBLK=IBLK+1
      IF((K1.GT.1).AND.(IC.EQ.3))  JBLK=IBLK-(LXP-LXM+1)+IS1
      IF((K1.LT.LY).AND.(IC.EQ.4)) JBLK=IBLK+(LXP-LXM+1)-IS2
      IF((K0.GT.1).AND.(IC.EQ.5))  JBLK=IBLK-LXY
      IF((K0.LT.LZ).AND.(IC.EQ.6)) JBLK=IBLK+LXY
      IF(IORI(JBLK).GE.5) ALB(IS+ISM(IC,II))=-ALB(IS+ISM(IC,II))
125   CONTINUE
*----
*  DIAG BOUNDARY CONDITION
*----
      IF (K1.EQ.K2) THEN
         IF(LL1.OR.LL2) THEN
            IKG=IGEN(IKK)
            IF(XX(IKG).NE.YY(IKG)) CALL XABORT('NUMER3: A CELL ON THE '
     1      //'DIAGONAL SYMMETRY AXIS IS NOT SQUARE.')
         ENDIF
         IF ((K1.EQ.1).AND.(NCODE(1).EQ.3).AND.(NCODE(3).EQ.5)) THEN
            FRX=0.25
            ALB(IS+ISM(1,II))=-ALB(IS+ISM(2,II))
            IFR(IS+ISM(1,II))=IFR(IS+ISM(2,II))
            DVX(IT+ISM(1,II))=-DVX(IT+ISM(2,II))
            MIXNEW=MIX(IT+ISM(2,II))
            MIXOLD=MIX(IT+ISM(1,II))
            DO 130 KC=1,NCOUR
            IF(MIX(IT+KC).EQ.MIXOLD) MIX(IT+KC)=MIXNEW
130         CONTINUE
         ELSE IF (NCODE(1).EQ.3) THEN
            FRX=0.5
            ALB(IS+ISM(1,II))=-ALB(IS+ISM(3,II))
            IFR(IS+ISM(1,II))=IFR(IS+ISM(3,II))
            DVX(IT+ISM(1,II))=-DVX(IT+ISM(3,II))
            MIXNEW=MIX(IT+ISM(3,II))
            MIXOLD=MIX(IT+ISM(1,II))
            DO 140 KC=1,NCOUR
            IF(MIX(IT+KC).EQ.MIXOLD) MIX(IT+KC)=MIXNEW
140         CONTINUE
         ENDIF
         IF ((K1.EQ.LY).AND.(NCODE(2).EQ.3).AND.(NCODE(4).EQ.5)) THEN
            FRX=0.25
            ALB(IS+ISM(2,II))=-ALB(IS+ISM(1,II))
            IFR(IS+ISM(2,II))=IFR(IS+ISM(1,II))
            DVX(IT+ISM(2,II))=-DVX(IT+ISM(1,II))
            MIXNEW=MIX(IT+ISM(1,II))
            MIXOLD=MIX(IT+ISM(2,II))
            DO 150 KC=1,NCOUR
            IF(MIX(IT+KC).EQ.MIXOLD) MIX(IT+KC)=MIXNEW
150         CONTINUE
         ELSE IF (NCODE(2).EQ.3) THEN
            ALB(IS+ISM(2,II))=-ALB(IS+ISM(4,II))
            IFR(IS+ISM(2,II))=IFR(IS+ISM(4,II))
            DVX(IT+ISM(2,II))=-DVX(IT+ISM(4,II))
            MIXNEW=MIX(IT+ISM(4,II))
            MIXOLD=MIX(IT+ISM(2,II))
            DO 160 KC=1,NCOUR
            IF(MIX(IT+KC).EQ.MIXOLD) MIX(IT+KC)=MIXNEW
160         CONTINUE
         ENDIF
         IF ((K1.EQ.1).AND.(NCODE(3).EQ.3).AND.(NCODE(1).EQ.5)) THEN
            FRY=0.25
            ALB(IS+ISM(3,II))=-ALB(IS+ISM(4,II))
            IFR(IS+ISM(3,II))=IFR(IS+ISM(4,II))
            DVX(IT+ISM(3,II))=-DVX(IT+ISM(4,II))
            MIXNEW=MIX(IT+ISM(4,II))
            MIXOLD=MIX(IT+ISM(3,II))
            DO 170 KC=1,NCOUR
            IF(MIX(IT+KC).EQ.MIXOLD) MIX(IT+KC)=MIXNEW
170         CONTINUE
         ELSE IF (NCODE(3).EQ.3) THEN
            FRY=0.5
            ALB(IS+ISM(3,II))=-ALB(IS+ISM(1,II))
            IFR(IS+ISM(3,II))=IFR(IS+ISM(1,II))
            DVX(IT+ISM(3,II))=-DVX(IT+ISM(1,II))
            MIXNEW=MIX(IT+ISM(1,II))
            MIXOLD=MIX(IT+ISM(3,II))
            DO 180 KC=1,NCOUR
            IF(MIX(IT+KC).EQ.MIXOLD) MIX(IT+KC)=MIXNEW
180         CONTINUE
         ENDIF
         IF ((K1.EQ.LY).AND.(NCODE(4).EQ.3).AND.(NCODE(2).EQ.5)) THEN
            FRY=0.25
            ALB(IS+ISM(4,II))=-ALB(IS+ISM(3,II))
            IFR(IS+ISM(4,II))=IFR(IS+ISM(3,II))
            DVX(IT+ISM(4,II))=-DVX(IT+ISM(3,II))
            MIXNEW=MIX(IT+ISM(3,II))
            MIXOLD=MIX(IT+ISM(4,II))
            DO 190 KC=1,NCOUR
            IF(MIX(IT+KC).EQ.MIXOLD) MIX(IT+KC)=MIXNEW
190         CONTINUE
         ELSE IF (NCODE(4).EQ.3) THEN
            ALB(IS+ISM(4,II))=-ALB(IS+ISM(2,II))
            IFR(IS+ISM(4,II))=IFR(IS+ISM(2,II))
            DVX(IT+ISM(4,II))=-DVX(IT+ISM(2,II))
            MIXNEW=MIX(IT+ISM(2,II))
            MIXOLD=MIX(IT+ISM(4,II))
            DO 200 KC=1,NCOUR
            IF(MIX(IT+KC).EQ.MIXOLD) MIX(IT+KC)=MIXNEW
200         CONTINUE
         ENDIF
      ENDIF
*----
*  TRAN BOUNDARY CONDITION
*----
      IF ((K2.EQ.1).AND.(NCODE(1).EQ.4)) THEN
         JBLK=IBLK+LXP-LXM
         IFR(IS+ISM(1,II))=NCOUR*(INUM(JBLK)-1)+ISM(2,IORI(JBLK))
      ENDIF
      IF ((K2.EQ.LX).AND.(NCODE(2).EQ.4)) THEN
         JBLK=IBLK+LXM-LXP
         IFR(IS+ISM(2,II))=NCOUR*(INUM(JBLK)-1)+ISM(1,IORI(JBLK))
      ENDIF
      IF ((K1.EQ.1).AND.(NCODE(3).EQ.4)) THEN
         JBLK=IBLK+(LY-1)*LX
         IFR(IS+ISM(3,II))=NCOUR*(INUM(JBLK)-1)+ISM(4,IORI(JBLK))
      ENDIF
      IF ((K1.EQ.LY).AND.(NCODE(4).EQ.4)) THEN
         JBLK=IBLK-(LY-1)*LX
         IFR(IS+ISM(4,II))=NCOUR*(INUM(JBLK)-1)+ISM(3,IORI(JBLK))
      ENDIF
      IF ((K0.EQ.1).AND.(NCODE(5).EQ.4)) THEN
         JBLK=IBLK+(LZ-1)*LXY
         IFR(IS+ISM(5,II))=NCOUR*(INUM(JBLK)-1)+ISM(6,IORI(JBLK))
      ENDIF
      IF ((K0.EQ.LZ).AND.(NCODE(6).EQ.4)) THEN
         JBLK=IBLK-(LZ-1)*LXY
         IFR(IS+ISM(6,II))=NCOUR*(INUM(JBLK)-1)+ISM(5,IORI(JBLK))
      ENDIF
*----
*  SYME BOUNDARY CONDITION
*----
      IF ((K2.EQ.1).AND.(NCODE(1).EQ.5)) THEN
         FRX=0.5
         ALB(IS+ISM(1,II))=-ALB(IS+ISM(2,II))
         IFR(IS+ISM(1,II))=IFR(IS+ISM(2,II))
         SUR(IS+ISM(3,II))=0.5*SUR(IS+ISM(3,II))
         SUR(IS+ISM(4,II))=0.5*SUR(IS+ISM(4,II))
         IF(ISM(5,II).NE.0) SUR(IS+ISM(5,II))=0.5*SUR(IS+ISM(5,II))
         IF(ISM(6,II).NE.0) SUR(IS+ISM(6,II))=0.5*SUR(IS+ISM(6,II))
         DVX(IT+ISM(1,II))=-DVX(IT+ISM(2,II))
         MIXNEW=MIX(IT+ISM(2,II))
         MIXOLD=MIX(IT+ISM(1,II))
         DO 210 KC=1,NCOUR
         IF(MIX(IT+KC).EQ.MIXOLD) MIX(IT+KC)=MIXNEW
210      CONTINUE
      ELSE IF ((K2.EQ.LX).AND.(NCODE(2).EQ.5)) THEN
         FRX=0.5
         ALB(IS+ISM(2,II))=-ALB(IS+ISM(1,II))
         IFR(IS+ISM(2,II))=IFR(IS+ISM(1,II))
         SUR(IS+ISM(3,II))=0.5*SUR(IS+ISM(3,II))
         SUR(IS+ISM(4,II))=0.5*SUR(IS+ISM(4,II))
         IF(ISM(5,II).NE.0) SUR(IS+ISM(5,II))=0.5*SUR(IS+ISM(5,II))
         IF(ISM(6,II).NE.0) SUR(IS+ISM(6,II))=0.5*SUR(IS+ISM(6,II))
         DVX(IT+ISM(2,II))=-DVX(IT+ISM(1,II))
         MIXNEW=MIX(IT+ISM(1,II))
         MIXOLD=MIX(IT+ISM(2,II))
         DO 220 KC=1,NCOUR
         IF(MIX(IT+KC).EQ.MIXOLD) MIX(IT+KC)=MIXNEW
220      CONTINUE
      ENDIF
      IF ((K1.EQ.1).AND.(NCODE(3).EQ.5)) THEN
         FRY=0.5
         ALB(IS+ISM(3,II))=-ALB(IS+ISM(4,II))
         IFR(IS+ISM(3,II))=IFR(IS+ISM(4,II))
         SUR(IS+ISM(1,II))=0.5*SUR(IS+ISM(1,II))
         SUR(IS+ISM(2,II))=0.5*SUR(IS+ISM(2,II))
         IF(ISM(5,II).NE.0) SUR(IS+ISM(5,II))=0.5*SUR(IS+ISM(5,II))
         IF(ISM(6,II).NE.0) SUR(IS+ISM(6,II))=0.5*SUR(IS+ISM(6,II))
         DVX(IT+ISM(3,II))=-DVX(IT+ISM(4,II))
         MIXNEW=MIX(IT+ISM(4,II))
         MIXOLD=MIX(IT+ISM(3,II))
         DO 230 KC=1,NCOUR
         IF(MIX(IT+KC).EQ.MIXOLD) MIX(IT+KC)=MIXNEW
230      CONTINUE
      ELSE IF ((K1.EQ.LY).AND.(NCODE(4).EQ.5)) THEN
         FRY=0.5
         ALB(IS+ISM(4,II))=-ALB(IS+ISM(3,II))
         IFR(IS+ISM(4,II))=IFR(IS+ISM(3,II))
         SUR(IS+ISM(1,II))=0.5*SUR(IS+ISM(1,II))
         SUR(IS+ISM(2,II))=0.5*SUR(IS+ISM(2,II))
         IF(ISM(5,II).NE.0) SUR(IS+ISM(5,II))=0.5*SUR(IS+ISM(5,II))
         IF(ISM(6,II).NE.0) SUR(IS+ISM(6,II))=0.5*SUR(IS+ISM(6,II))
         DVX(IT+ISM(4,II))=-DVX(IT+ISM(3,II))
         MIXNEW=MIX(IT+ISM(3,II))
         MIXOLD=MIX(IT+ISM(4,II))
         DO 240 KC=1,NCOUR
         IF(MIX(IT+KC).EQ.MIXOLD) MIX(IT+KC)=MIXNEW
240      CONTINUE
      ENDIF
      IF ((K0.EQ.1).AND.(NCODE(5).EQ.5)) THEN
         FRZ=0.5
         ALB(IS+ISM(5,II))=-ALB(IS+ISM(6,II))
         IFR(IS+ISM(5,II))=IFR(IS+ISM(6,II))
         SUR(IS+ISM(1,II))=0.5*SUR(IS+ISM(1,II))
         SUR(IS+ISM(2,II))=0.5*SUR(IS+ISM(2,II))
         SUR(IS+ISM(3,II))=0.5*SUR(IS+ISM(3,II))
         SUR(IS+ISM(4,II))=0.5*SUR(IS+ISM(4,II))
         DVX(IT+ISM(5,II))=-DVX(IT+ISM(6,II))
         MIXNEW=MIX(IT+ISM(6,II))
         MIXOLD=MIX(IT+ISM(5,II))
         DO 250 KC=1,NCOUR
         IF(MIX(IT+KC).EQ.MIXOLD) MIX(IT+KC)=MIXNEW
250      CONTINUE
      ELSE IF ((K0.EQ.LZ).AND.(NCODE(6).EQ.5)) THEN
         FRZ=0.5
         ALB(IS+ISM(6,II))=-ALB(IS+ISM(5,II))
         IFR(IS+ISM(6,II))=IFR(IS+ISM(5,II))
         SUR(IS+ISM(1,II))=0.5*SUR(IS+ISM(1,II))
         SUR(IS+ISM(2,II))=0.5*SUR(IS+ISM(2,II))
         SUR(IS+ISM(3,II))=0.5*SUR(IS+ISM(3,II))
         SUR(IS+ISM(4,II))=0.5*SUR(IS+ISM(4,II))
         DVX(IT+ISM(6,II))=-DVX(IT+ISM(5,II))
         MIXNEW=MIX(IT+ISM(5,II))
         MIXOLD=MIX(IT+ISM(6,II))
         DO 260 KC=1,NCOUR
         IF(MIX(IT+KC).EQ.MIXOLD) MIX(IT+KC)=MIXNEW
260      CONTINUE
      ENDIF
*
265   GG3(IBLK)=FRX*FRY*FRZ
270   CONTINUE
275   CONTINUE
280   CONTINUE
      DO 285 I=1,NCOUR*NMBLK
      IFR(I)=MIX(IFR(I))
285   CONTINUE
*----
*  ELIMINATION OF THE BLOCKS OUTSIDE THE DOMAIN
*----
      JBLK=0
      DO 300 IBLK=1,NMBLK
      IKK=INUM(IBLK)
      IF (IKK.GT.0) THEN
         JBLK=JBLK+1
         INUM(JBLK)=IKK
         IORI(JBLK)=IORI(IBLK)
         GG3(JBLK)=GG3(IBLK)
         IS=NCOUR*(IBLK-1)
         JS=NCOUR*(JBLK-1)
         DO 290 IC=1,NCOUR
         IFR(JS+IC)=IFR(IS+IC)
         ALB(JS+IC)=ALB(IS+IC)
         SUR(JS+IC)=SUR(IS+IC)
290      CONTINUE
      ENDIF
300   CONTINUE
      NMBLK=JBLK
*
      DO 310 IKK=1,NMERGE
      POURCE(IKK)=0.0
310   CONTINUE
      DO 330 IBLK=1,NMBLK
      IKK=INUM(IBLK)
      POURCE(IKK)=POURCE(IKK)+GG3(IBLK)
330   CONTINUE
*----
*  VALIDATION OF VECTOR IFR
*----
      DO 345 IBLK=1,NMBLK
      IS=NCOUR*(IBLK-1)
      DO 340 IC=1,NCOUR
      ISURF=IFR(IS+IC)
      IF (ISURF.EQ.0) THEN
         WRITE (HSMG,'(44HNUMER3: FAILURE OF THE SURFACE RENUMBERING A,
     1   12HLGORITHM(1).)')
         GO TO 570
      ENDIF
      JC=1+MOD(ISURF-1,NCOUR)
      JT=NCOUR*((ISURF-1)/NCOUR)
      IF (MIX(JT+JC).NE.ISURF) THEN
         WRITE (HSMG,'(44HNUMER3: FAILURE OF THE SURFACE RENUMBERING A,
     1   12HLGORITHM(2).)')
         GO TO 570
      ENDIF
340   CONTINUE
345   CONTINUE
*----
*  VALIDATION OF THE GEOMETRICAL RECIPROCITY AND MODIFICATION OF ALBEDOS
*----
      DO 355 IBLK=1,NMBLK
      IKK=INUM(IBLK)
      IKG=IGEN(IKK)
      IS=NCOUR*(IBLK-1)
      IT=NCOUR*(IKK-1)
      DO 350 IC=1,NCOUR
      ALB(IS+IC)=ALB(IS+IC)*GG3(IBLK)/POURCE(IKK)
      IP=1+MOD(MIX(IT+IC)-1,NCOUR)
      FR1=YY(IKG)*ZZ(IKG)
      IF ((IP.EQ.3).OR.(IP.EQ.4)) FR1=XX(IKG)*ZZ(IKG)
      IF ((IP.EQ.5).OR.(IP.EQ.6)) FR1=XX(IKG)*YY(IKG)
      JP=1+MOD(IFR(IS+IC)-1,NCOUR)
      JKG=IGEN(1+(IFR(IS+IC)-1)/NCOUR)
      FR2=YY(JKG)*ZZ(JKG)
      IF ((JP.EQ.3).OR.(JP.EQ.4)) FR2=XX(JKG)*ZZ(JKG)
      IF ((JP.EQ.5).OR.(JP.EQ.6)) FR2=XX(JKG)*YY(JKG)
      DELTA=ABS(FR1-FR2)
      IF (ABS(FR1-FR2).GT.EPS) THEN
         WRITE (HSMG,680) DIRR(IP),IKG,DIRR(JP),JKG
         GO TO 570
      ENDIF
350   CONTINUE
355   CONTINUE
*----
*  COMPUTE VECTOR DVX
*----
      DO 395 IKK=1,NMERGE
      IKG=IGEN(IKK)
      IF (NCOUR.EQ.2) THEN
         DDD(1)=0.5
         DDD(2)=0.5
      ELSE IF (NCOUR.EQ.4) THEN
         SURFA=2.0*(XX(IKG)+YY(IKG))
         DO 360 IC=1,NCOUR
         FR1=YY(IKG)
         IF ((IC.EQ.3).OR.(IC.EQ.4)) FR1=XX(IKG)
         DDD(IC)=FR1/SURFA
360      CONTINUE
      ELSE IF (NCOUR.EQ.6) THEN
         SURFA=2.0*(XX(IKG)*ZZ(IKG)+YY(IKG)*ZZ(IKG)+XX(IKG)*YY(IKG))
         DO 370 IC=1,NCOUR
         FR1=YY(IKG)*ZZ(IKG)
         IF ((IC.EQ.3).OR.(IC.EQ.4)) FR1=XX(IKG)*ZZ(IKG)
         IF ((IC.EQ.5).OR.(IC.EQ.6)) FR1=XX(IKG)*YY(IKG)
         DDD(IC)=FR1/SURFA
370      CONTINUE
      ENDIF
      IT=NCOUR*(IKK-1)
      DO 390 IC=1,NCOUR
      IF (MULTC.EQ.1) THEN
*        ROTH APPROXIMATION.
         DVX(IT+IC)=DDD(IC)
         MIX(IT+IC)=IKK
      ELSE
         DELTA=0.0
         I1=MIX(IT+IC)
         DO 380 JC=1,NCOUR
         IF (MIX(IT+JC).EQ.I1) DELTA=DELTA+DDD(JC)
380      CONTINUE
         ZSIGN=SIGN(1.0,DVX(IT+IC))
         DVX(IT+IC)=ZSIGN*DDD(IC)/DELTA
      ENDIF
390   CONTINUE
395   CONTINUE
      IJAS=NCOUR*NMBLK
      IJAR=NCOUR*NMERGE
*----
*  RECOMPUTE VECTOR IFR FOR ROTH APPROXIMATION
*----
      IF (MULTC.EQ.1) THEN
         DO 400 I=1,IJAS
         IFR(I)=1+(IFR(I)-1)/NCOUR
400      CONTINUE
      ENDIF
*----
*  REMOVE THE UNUSED SURFACE NUMBERS
*----
      DO 410 I=1,IJAS
      JF2(I)=0
410   CONTINUE
      IJAT=0
      DO 420 I=1,IJAR
      J=MIX(I)
      IF (J.GT.IJAS) THEN
         WRITE (HSMG,'(44HNUMER3: FAILURE OF THE SURFACE RENUMBERING A,
     1   12HLGORITHM(3).)')
         GO TO 570
      ENDIF
      IF (JF2(J).EQ.0) THEN
         IJAT=IJAT+1
         JF2(J)=IJAT
      ENDIF
420   CONTINUE
      DO 430 I=1,IJAR
      MIX(I)=JF2(MIX(I))
430   CONTINUE
      DO 440 I=1,IJAS
      IFR(I)=JF2(IFR(I))
440   CONTINUE
*----
*  INCLUDE THE DP-1 APPROXIMATION
*----
      IF ((MULTC.EQ.4).AND.(NCOUR.EQ.2)) THEN
*        DP-1 APPROXIMATION IN 1-D.
         DO 455 I1=IJAR,1,-1
         FR1=ABS(DVX(I1))
         JND=(MIX(I1)-1)*2
         DO 450 JCOUR=1,2
         JSURF=(I1-1)*2+JCOUR
         DVX(JSURF)=FR1
         MIX(JSURF)=JND+JCOUR
450      CONTINUE
455      CONTINUE
         DO 465 I1=IJAS,1,-1
         FR1=ABS(ALB(I1))
         FR2=SUR(I1)
         JND=(IFR(I1)-1)*2
         DO 460 JCOUR=1,2
         JSURF=(I1-1)*2+JCOUR
         ALB(JSURF)=FR1
         SUR(JSURF)=FR2
         IFR(JSURF)=JND+JCOUR
460      CONTINUE
465      CONTINUE
         NCOUR=4
      ELSE IF ((MULTC.EQ.4).AND.(NCOUR.EQ.4)) THEN
*        DP-1 APPROXIMATION IN 2-D.
         DO 480 I1=IJAR,1,-1
         ZSIGN=SIGN(1.0,DVX(I1))
         FR1=ABS(DVX(I1))
         FR2=SUR(I1)
         JND=(MIX(I1)-1)*3
         DO 470 JCOUR=1,3
         JSURF=(I1-1)*3+JCOUR
         DVX(JSURF)=FR1
         MIX(JSURF)=JND+JCOUR
470      CONTINUE
         DVX(JSURF)=ZSIGN*FR1
480      CONTINUE
         DO 500 I1=IJAS,1,-1
         ZSIGN=SIGN(1.0,ALB(I1))
         FR1=ABS(ALB(I1))
         FR2=SUR(I1)
         JND=(IFR(I1)-1)*3
         DO 490 JCOUR=1,3
         JSURF=(I1-1)*3+JCOUR
         ALB(JSURF)=FR1
         SUR(JSURF)=FR2
         IFR(JSURF)=JND+JCOUR
490      CONTINUE
         ALB(JSURF)=ZSIGN*FR1
500      CONTINUE
         NCOUR=12
      ELSE IF ((MULTC.EQ.4).AND.(NCOUR.EQ.6)) THEN
         CALL XABORT('NUMER3: INVALID OPTION.')
      ELSE
         DO 510 I=1,IJAS
         ALB(I)=ABS(ALB(I))
510      CONTINUE
         DO 520 I=1,IJAR
         DVX(I)=ABS(DVX(I))
520      CONTINUE
      ENDIF
      IJAS=NCOUR*NMBLK
      IJAR=NCOUR*NMERGE
*----
*  PRINT THE SURFACE NUMBERS AFTER MERGING
*----
      IF (IMPX.GT.2) THEN
         WRITE (6,620)
         MIN6=MIN(6,NCOUR)
         WRITE (6,650) ('----------------',I=1,MIN6)
         DO 560 IBLK=1,NMBLK
         IKK=INUM(IBLK)
         WRITE (6,630) IBLK,IKK,IGEN(IKK)
         I1=IORI(IBLK)
         IF ((MULTC.EQ.4).AND.(NCOUR.EQ.12)) THEN
            DO 530 I=1,12
            ISMZ(I)=3*ISM(1+(I-1)/3,I1)+MOD(I-1,3)-2
            DIRZ(I)=DIRR(1+(I-1)/3)
530         CONTINUE
         ELSE IF ((MULTC.EQ.4).AND.(NCOUR.EQ.4)) THEN
            DO 540 I=1,4
            ISMZ(I)=2*ISM(1+(I-1)/2,I1)+MOD(I-1,2)-1
            DIRZ(I)=DIRR(1+(I-1)/2)
540         CONTINUE
         ELSE
            DO 550 I=1,NCOUR
            ISMZ(I)=ISM(I,I1)
            DIRZ(I)=DIRR(I)
550         CONTINUE
         ENDIF
         IT0=NCOUR*(IBLK-1)
         IT1=NCOUR*(IKK-1)
         WRITE (6,660) (DIRZ(I),I=1,MIN6)
         WRITE (6,635) (MIX(IT1+ISMZ(I)),IFR(IT0+ISMZ(I)),I=1,MIN6)
         WRITE (6,640) (ALB(IT0+ISMZ(I)),I=1,MIN6)
         WRITE (6,645) (DVX(IT1+ISMZ(I)),I=1,MIN6)
         IF (NCOUR.EQ.12) THEN
            WRITE (6,660) (DIRZ(I),I=7,12)
            WRITE (6,635) (MIX(IT1+ISMZ(I)),IFR(IT0+ISMZ(I)),I=7,12)
            WRITE (6,640) (ALB(IT0+ISMZ(I)),I=7,12)
            WRITE (6,645) (DVX(IT1+ISMZ(I)),I=7,12)
         ENDIF
         WRITE (6,650) ('----------------',I=1,MIN6)
560      CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(GG3,JF2)
      RETURN
*
570   WRITE (6,620)
      WRITE (6,650) ('----------------',I=1,NCOUR)
      DO 580 IBLK=1,NMBLK
      IKK=INUM(IBLK)
      WRITE (6,630) IBLK,IKK,IGEN(IKK)
      I1=IORI(IBLK)
      IT0=NCOUR*(IBLK-1)
      IT1=NCOUR*(IKK-1)
      WRITE (6,660) (DIRR(I),I=1,NCOUR)
      WRITE (6,635) (MIX(IT1+ISM(I,I1)),IFR(IT0+ISM(I,I1)),I=1,NCOUR)
      WRITE (6,640) (ALB(IT0+ISM(I,I1)),I=1,NCOUR)
      WRITE (6,645) (DVX(IT1+ISM(I,I1)),I=1,NCOUR)
      WRITE (6,650) ('----------------',I=1,NCOUR)
580   CONTINUE
      CALL XABORT(HSMG)
*
620   FORMAT (///31H SURFACE NUMBERS AFTER MERGING./)
630   FORMAT (7H BLOCK=,I5,5X,13HMERGED BLOCK=,I5,5X,12HGENERATING B,
     1 5HLOCK=,I5)
635   FORMAT (8H IN/OUT:,6(I6,2H /,I5,3H  I))
640   FORMAT (8H ALBEDO:,1P,6(E13.5,3H  I))
645   FORMAT (8H    DVX:,1P,6(E13.5,3H  I))
650   FORMAT (8H -------,6(A16))
660   FORMAT (/8H   SIDE:,6(A9,6X,1HI))
680   FORMAT (49HNUMER3: GEOMETRICAL RECIPROCITY CONDITION IS VIOL,
     1 10HATED (SIDE,A3,20H OF GENERATING BLOCK,I5,8H VS SIDE,A3,
     2 20H OF GENERATING BLOCK,I5,3H ).)
      END
