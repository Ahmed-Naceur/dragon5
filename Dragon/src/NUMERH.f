*DECK NUMERH
      SUBROUTINE NUMERH (NCOUR,MULTC,NCODE,ZCODE,IHEX,LX,LZ,IORI,ISM,
     1 POURCE,IMPX,NMBLK,IFR,ALB,SUR,NMERGE,INUM,MIX,DVX,NGEN,IGEN,XX,
     2 ZZ)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Surface renumbering algorithm for hexagonal geometry.
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
*         =3 not used;   =4 TRAN;     =5 SYME.
* ZCODE   value of the albedo on each side of the domain (only
*         components 1, 5 and 6 are used).
* IHEX    type of symmetry:
*         =1 S30;   =2 SA60;   =3 SB60;   =4 S90;   =5 R120;
*         =6 R180;  =7 SA180;  =8 SB180;  =9 COMPLETE.
* LX      number of blocks in the X-Y (hexagonal) plane.
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
* XX      side of the generating hexagons.
* ZZ      Z-thickness of the generating blocks (hexagons).
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NCOUR,MULTC,NCODE(6),IHEX,LX,LZ,IORI(NMBLK),ISM(8,12),
     1 IMPX,NMBLK,IFR(18*NMBLK),NMERGE,INUM(NMBLK),MIX(18*NMERGE),
     2 NGEN,IGEN(NMERGE)
      REAL ZCODE(6),POURCE(NMERGE),ALB(18*NMBLK),SUR(18*NMBLK),
     1 DVX(18*NMERGE),XX(NGEN),ZZ(NGEN)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (EPS=1.0E-5)
      LOGICAL LOG1,LTRAN
      CHARACTER DIRR(8)*2,DIRZ(18)*2,HSMG*131
      INTEGER ICCC(6),IDDD(8),ISMZ(18)
      REAL DDD(8)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: JF2
      REAL, ALLOCATABLE, DIMENSION(:) :: GG3
      SAVE DIRR
      DATA DIRR/'H1','H2','H3','H4','H5','H6','Z-','Z+'/
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(JF2(8*NMBLK),GG3(NMBLK))
*
      DO 5 I=1,NCOUR*NMERGE
      MIX(I)=I
      DVX(I)=1.0
5     CONTINUE
      LTRAN=(IHEX.EQ.5).OR.(IHEX.EQ.6)
      ZALB=1.0
      IF(NCODE(1).EQ.1) ZALB=ZCODE(1)
      IBLK=0
      DO 52 K0=1,LZ
      DO 51 K1=1,LX
      IBLK=IBLK+1
      IKK=INUM(IBLK)
      FRZ=1.0
      IF(IKK.EQ.0) GO TO 50
      IS=NCOUR*(IBLK-1)
      IT=NCOUR*(IKK-1)
      II=IORI(IBLK)
      DO 10 IC=1,NCOUR
      IDDD(IC)=-1
      ALB(IS+IC)=1.0
      SUR(IS+IC)=0.0
      IF(IC.LE.6) THEN
         ICCC(IC)=NEIGHB(K1,IC,IHEX,LX,FRX)
         IF((ICCC(IC).GT.0).AND.(ICCC(IC).NE.K1).AND.(ICCC(IC).LE.LX))
     1   IDDD(IC)=ICCC(IC)+(K0-1)*LX
      ELSE IF((IC.EQ.7).AND.(K0.GT.1)) THEN
         IDDD(7)=IBLK-LX
      ELSE IF((IC.EQ.8).AND.(K0.LT.LZ)) THEN
         IDDD(8)=IBLK+LX
      ENDIF
      JBLK=IDDD(IC)
      IF(JBLK.GT.0) THEN
         JKK=INUM(JBLK)
         JT=NCOUR*(JKK-1)
         IF((IC.LE.6).AND.(JKK.GT.0)) THEN
            IFR(IS+ISM(IC,II))=JT+ISM(1+MOD(IC+2,6),IORI(JBLK))
         ELSE IF((IC.EQ.7).AND.(JKK.GT.0)) THEN
            IFR(IS+ISM(7,II))=JT+ISM(8,IORI(JBLK))
         ELSE IF((IC.EQ.8).AND.(JKK.GT.0)) THEN
            IFR(IS+ISM(8,II))=JT+ISM(7,IORI(JBLK))
         ELSE
            IFR(IS+ISM(IC,II))=0
         ENDIF
         IDDD(IC)=JKK
      ELSE
         IFR(IS+ISM(IC,II))=0
      ENDIF
10    CONTINUE
      IF(IMPX.GT.10) WRITE(6,610) K1,(ICCC(IC),IC=1,6),FRX
*----
*  VOID OR REFL BOUNDARY CONDITIONS
*----
      IKG=IGEN(IKK)
      DO 15 IC=1,6
      IF(ICCC(IC).GT.LX) THEN
         ALB(IS+ISM(IC,II))=-ZALB
         IFR(IS+ISM(IC,II))=IT+ISM(IC,II)
         SUR(IS+ISM(IC,II))=XX(IKG)*ZZ(IKG)
      ELSE IF(ICCC(IC).EQ.K1) THEN
         ALB(IS+ISM(IC,II))=-1.0
         IFR(IS+ISM(IC,II))=IT+ISM(IC,II)
      ENDIF
15    CONTINUE
      SURFA=1.5*SQRT(3.0)*XX(IKG)*XX(IKG)
      LOG1=(K0.EQ.1).OR.(IDDD(7).EQ.0)
      IF(LOG1.AND.(NCODE(5).EQ.1)) THEN
         ALB(IS+ISM(7,II))=-ZCODE(5)
         IFR(IS+ISM(7,II))=IT+ISM(7,II)
      ELSE IF(LOG1.AND.(NCODE(5).EQ.2)) THEN
         ALB(IS+ISM(7,II))=-1.0
         IFR(IS+ISM(7,II))=IT+ISM(7,II)
      ENDIF
      IF(LOG1.AND.(NCODE(5).EQ.1)) SUR(IS+ISM(7,II))=SURFA
      IF(LOG1.AND.(NCODE(5).EQ.2)) SUR(IS+ISM(7,II))=SURFA
      LOG1=(K0.EQ.LZ).OR.(IDDD(8).EQ.0)
      IF(LOG1.AND.(NCODE(6).EQ.1)) THEN
         ALB(IS+ISM(8,II))=-ZCODE(6)
         IFR(IS+ISM(8,II))=IT+ISM(8,II)
      ELSE IF(LOG1.AND.(NCODE(6).EQ.2)) THEN
         ALB(IS+ISM(8,II))=-1.0
         IFR(IS+ISM(8,II))=IT+ISM(8,II)
      ENDIF
      IF(LOG1.AND.(NCODE(6).EQ.1)) SUR(IS+ISM(8,II))=SURFA
      IF(LOG1.AND.(NCODE(6).EQ.2)) SUR(IS+ISM(8,II))=SURFA
*----
*  CORRECT THE PARITY OF THE INTERFACE CURRENTS FOR DP-1 CASES WITH
*  'MIRROR' ORIENTATION
*----
      DO 20 IC=1,NCOUR
      IF(II.GE.7) ALB(IS+IC)=-ALB(IS+IC)
      JBLK=IBLK
      IF(IC.LE.6) THEN
         IF((ICCC(IC).GT.0).AND.(ICCC(IC).NE.K1).AND.(ICCC(IC).LE.LX))
     1   JBLK=ICCC(IC)+(K0-1)*LX
      ELSE IF((IC.EQ.7).AND.(K0.GT.1)) THEN
         JBLK=IBLK-LX
      ELSE IF((IC.EQ.8).AND.(K0.LT.LZ)) THEN
         JBLK=IBLK+LX
      ENDIF
      IF(IORI(JBLK).GE.7) ALB(IS+ISM(IC,II))=-ALB(IS+ISM(IC,II))
20    CONTINUE
*----
*  TRAN BOUNDARY CONDITION
*----
      DO 30 IC=1,6
      IF(LTRAN.AND.(ICCC(IC).LT.0)) THEN
         JBLK=-ICCC(IC)+(K0-1)*LX
         JT=NCOUR*(INUM(JBLK)-1)
         K1M=-K1
         IF(K1.EQ.1) K1M=1
         JC=0
         DO 25 KC=1,6
         IF(NEIGHB(-ICCC(IC),KC,IHEX,LX,FRW).EQ.K1M) JC=KC
25       CONTINUE
         IF(JC.EQ.0) CALL XABORT('NUMERH: FAILURE NB. 1.')
         IFR(IS+ISM(IC,II))=JT+ISM(JC,IORI(JBLK))
      ENDIF
30    CONTINUE
      IF((K0.EQ.1).AND.(NCODE(5).EQ.4)) THEN
         JBLK=IBLK+(LZ-1)*LX
         IFR(IS+ISM(7,II))=NCOUR*(INUM(JBLK)-1)+ISM(8,IORI(JBLK))
      ENDIF
      IF((K0.EQ.LZ).AND.(NCODE(6).EQ.4)) THEN
         JBLK=IBLK-(LZ-1)*LX
         IFR(IS+ISM(8,II))=NCOUR*(INUM(JBLK)-1)+ISM(7,IORI(JBLK))
      ENDIF
*----
*  SYME BOUNDARY CONDITION
*----
      DO 40 IC=1,6
      IF((.NOT.LTRAN).AND.(ICCC(IC).LT.0)) THEN
         JC=0
         DO 35 KC=1,6
         IF(-ICCC(IC).EQ.ICCC(KC)) JC=KC
35       CONTINUE
         IF(JC.EQ.0) CALL XABORT('NUMERH: FAILURE NB. 2.')
         ZSIGN=-1.0
         IF(K1.EQ.1) THEN
            IF((IHEX.EQ.3).AND.(IC.EQ.3)) ZSIGN=1.0
            IF((IHEX.EQ.3).AND.(IC.EQ.5)) ZSIGN=1.0
            IF((IHEX.EQ.4).AND.(IC.EQ.4)) ZSIGN=1.0
         ENDIF
         IF((IHEX.EQ.10).AND.(FRX.EQ.0.25).AND.(IC.EQ.6)) ZSIGN=1.0
         ALB(IS+ISM(IC,II))=ALB(IS+ISM(JC,II))*ZSIGN
         IFR(IS+ISM(IC,II))=IFR(IS+ISM(JC,II))
         DVX(IT+ISM(IC,II))=DVX(IT+ISM(JC,II))*ZSIGN
         MIXOLD=MIX(IT+ISM(IC,II))
         MIXNEW=MIX(IT+ISM(JC,II))
         DO 36 KC=1,NCOUR
         IF(MIX(IT+KC).EQ.MIXOLD) MIX(IT+KC)=MIXNEW
36       CONTINUE
      ENDIF
40    CONTINUE
      IF((ICCC(2).EQ.-ICCC(6)).AND.(ICCC(3).EQ.-ICCC(5))) THEN
         SUR(IS+ISM(1,II))=0.5*SUR(IS+ISM(1,II))
         SUR(IS+ISM(4,II))=0.5*SUR(IS+ISM(4,II))
      ELSE IF((ICCC(1).EQ.-ICCC(3)).AND.(ICCC(4).EQ.-ICCC(6))) THEN
         SUR(IS+ISM(2,II))=0.5*SUR(IS+ISM(2,II))
         SUR(IS+ISM(5,II))=0.5*SUR(IS+ISM(5,II))
      ELSE IF((ICCC(1).EQ.-ICCC(5)).AND.(ICCC(2).EQ.-ICCC(4))) THEN
         SUR(IS+ISM(3,II))=0.5*SUR(IS+ISM(3,II))
         SUR(IS+ISM(6,II))=0.5*SUR(IS+ISM(6,II))
      ENDIF
      IF((K0.EQ.1).AND.(NCODE(5).EQ.5)) THEN
         FRZ=0.5
         ALB(IS+ISM(7,II))=-ALB(IS+ISM(8,II))
         IFR(IS+ISM(7,II))=IFR(IS+ISM(8,II))
         SUR(IS+ISM(1:6,II))=0.5*SUR(IS+ISM(1:6,II))
         DVX(IT+ISM(7,II))=-DVX(IT+ISM(8,II))
         MIXOLD=MIX(IT+ISM(7,II))
         MIXNEW=MIX(IT+ISM(8,II))
      ELSE IF((K0.EQ.LZ).AND.(NCODE(6).EQ.5)) THEN
         FRZ=0.5
         ALB(IS+ISM(8,II))=-ALB(IS+ISM(7,II))
         IFR(IS+ISM(8,II))=IFR(IS+ISM(7,II))
         SUR(IS+ISM(1:6,II))=0.5*SUR(IS+ISM(1:6,II))
         DVX(IT+ISM(8,II))=-DVX(IT+ISM(7,II))
         MIXOLD=MIX(IT+ISM(8,II))
         MIXNEW=MIX(IT+ISM(7,II))
      ENDIF
      IF(FRZ.EQ.0.5) THEN
         DO 45 KC=1,NCOUR
         IF(MIX(IT+KC).EQ.MIXOLD) MIX(IT+KC)=MIXNEW
45       CONTINUE
      ENDIF
50    GG3(IBLK)=FRX*FRZ
51    CONTINUE
52    CONTINUE
      DO 60 I=1,NCOUR*NMBLK
      IFR(I)=MIX(IFR(I))
60    CONTINUE
*----
*  ELIMINATION OF THE BLOCKS OUTSIDE THE DOMAIN
*----
      JBLK=0
      DO 80 IBLK=1,NMBLK
      IKK=INUM(IBLK)
      IF(IKK.GT.0) THEN
         JBLK=JBLK+1
         INUM(JBLK)=IKK
         IORI(JBLK)=IORI(IBLK)
         GG3(JBLK)=GG3(IBLK)
         IS=NCOUR*(IBLK-1)
         JS=NCOUR*(JBLK-1)
         DO 70 IC=1,NCOUR
         IFR(JS+IC)=IFR(IS+IC)
         ALB(JS+IC)=ALB(IS+IC)
         SUR(JS+IC)=SUR(IS+IC)
70       CONTINUE
      ENDIF
80    CONTINUE
      NMBLK=JBLK
*
      DO 90 IKK=1,NMERGE
      POURCE(IKK)=0.0
90    CONTINUE
      DO 110 IBLK=1,NMBLK
      IKK=INUM(IBLK)
      POURCE(IKK)=POURCE(IKK)+GG3(IBLK)
110   CONTINUE
*----
*  VALIDATION OF VECTOR IFR
*----
      DO 125 IBLK=1,NMBLK
      IS=NCOUR*(IBLK-1)
      DO 120 IC=1,NCOUR
      ISURF=IFR(IS+IC)
      IF(ISURF.EQ.0) THEN
         GO TO 570
      ENDIF
      JC=1+MOD(ISURF-1,NCOUR)
      JT=NCOUR*((ISURF-1)/NCOUR)
      IF(MIX(JT+JC).NE.ISURF) THEN
         GO TO 570
      ENDIF
120   CONTINUE
125   CONTINUE
*----
*  VALIDATION OF THE GEOMETRICAL RECIPROCITY AND MODIFICATION OF ALBEDOS
*----
      DO 135 IBLK=1,NMBLK
      IKK=INUM(IBLK)
      IKG=IGEN(IKK)
      IS=NCOUR*(IBLK-1)
      IT=NCOUR*(IKK-1)
      DO 130 IC=1,NCOUR
      ALB(IS+IC)=ALB(IS+IC)*GG3(IBLK)/POURCE(IKK)
      IP=1+MOD(MIX(IT+IC)-1,NCOUR)
      FR1=XX(IKG)*ZZ(IKG)
      IF((IP.EQ.7).OR.(IP.EQ.8)) FR1=1.5*SQRT(3.0)*XX(IKG)*XX(IKG)
      JP=1+MOD(IFR(IS+IC)-1,NCOUR)
      JKG=IGEN(1+(IFR(IS+IC)-1)/NCOUR)
      FR2=XX(IKG)*ZZ(JKG)
      IF((JP.EQ.7).OR.(JP.EQ.8)) FR2=1.5*SQRT(3.0)*XX(IKG)*XX(IKG)
      DELTA=ABS(FR1-FR2)
      IF(ABS(FR1-FR2).GT.EPS) THEN
         WRITE (HSMG,680) DIRR(IP),IKG,DIRR(JP),JKG
         CALL XABORT(HSMG)
      ENDIF
130   CONTINUE
135   CONTINUE
*----
*  COMPUTE VECTOR DVX
*----
      DO 205 IKK=1,NMERGE
      IKG=IGEN(IKK)
      IF(NCOUR.EQ.6) THEN
         DO 140 IC=1,NCOUR
         DDD(IC)=1.0/6.0
140      CONTINUE
      ELSE IF(NCOUR.EQ.8) THEN
         SURFA=1.5*SQRT(3.0)*XX(IKG)*XX(IKG)
         DO 150 IC=1,6
         DDD(IC)=XX(IKG)*ZZ(IKG)/(2.0*SURFA+6.0*XX(IKG)*ZZ(IKG))
150      CONTINUE
         DDD(7)=SURFA/(2.0*SURFA+6.0*XX(IKG)*ZZ(IKG))
         DDD(8)=DDD(7)
      ENDIF
      IT=NCOUR*(IKK-1)
      DO 200 IC=1,NCOUR
      IF(MULTC.EQ.1) THEN
*        ROTH APPROXIMATION.
         DVX(IT+IC)=DDD(IC)
         MIX(IT+IC)=IKK
      ELSE
         DELTA=0.0
         I1=MIX(IT+IC)
         DO 180 JC=1,NCOUR
         IF(MIX(IT+JC).EQ.I1) DELTA=DELTA+DDD(JC)
180      CONTINUE
         ZSIGN=SIGN(1.0,DVX(IT+IC))
         DVX(IT+IC)=ZSIGN*DDD(IC)/DELTA
      ENDIF
200   CONTINUE
205   CONTINUE
      IJAS=NCOUR*NMBLK
      IJAR=NCOUR*NMERGE
*----
*  RECOMPUTE VECTOR IFR FOR ROTH APPROXIMATION
*----
      IF(MULTC.EQ.1) THEN
         DO 210 I=1,IJAS
         IFR(I)=1+(IFR(I)-1)/NCOUR
210      CONTINUE
      ENDIF
*----
*  REMOVE THE UNUSED SURFACE NUMBERS
*----
      DO 240 I=1,IJAS
      JF2(I)=0
240   CONTINUE
      IJAT=0
      DO 250 I=1,IJAR
      J=MIX(I)
      IF(J.GT.IJAS) THEN
        GO TO 570
      ENDIF
      IF(JF2(J).EQ.0) THEN
         IJAT=IJAT+1
         JF2(J)=IJAT
      ENDIF
250   CONTINUE
      DO 260 I=1,IJAR
      MIX(I)=JF2(MIX(I))
260   CONTINUE
      DO 270 I=1,IJAS
      IFR(I)=JF2(IFR(I))
270   CONTINUE
*----
*  INCLUDE THE DP-1 APPROXIMATION
*----
      IF((MULTC.EQ.4).AND.(NCOUR.EQ.6)) THEN
*        DP-1 APPROXIMATION IN 2-D.
         DO 330 I1=IJAR,1,-1
         ZSIGN=SIGN(1.0,DVX(I1))
         FR1=ABS(DVX(I1))
         FR2=SUR(I1)
         JND=(MIX(I1)-1)*3
         DO 320 JCOUR=1,3
         JSURF=(I1-1)*3+JCOUR
         DVX(JSURF)=FR1
         MIX(JSURF)=JND+JCOUR
320      CONTINUE
         DVX(JSURF)=ZSIGN*FR1
330      CONTINUE
         DO 350 I1=IJAS,1,-1
         ZSIGN=SIGN(1.0,ALB(I1))
         FR1=ABS(ALB(I1))
         FR2=SUR(I1)
         JND=(IFR(I1)-1)*3
         DO 340 JCOUR=1,3
         JSURF=(I1-1)*3+JCOUR
         ALB(JSURF)=FR1
         SUR(JSURF)=FR2
         IFR(JSURF)=JND+JCOUR
340      CONTINUE
         ALB(JSURF)=ZSIGN*FR1
350      CONTINUE
         NCOUR=18
      ELSE IF((MULTC.EQ.4).AND.(NCOUR.EQ.8)) THEN
         CALL XABORT('NUMERH: INVALID OPTION.')
      ELSE
         DO 360 I=1,IJAS
         ALB(I)=ABS(ALB(I))
360      CONTINUE
         DO 370 I=1,IJAR
         DVX(I)=ABS(DVX(I))
370      CONTINUE
      ENDIF
      IJAS=NCOUR*NMBLK
      IJAR=NCOUR*NMERGE
*----
*  PRINT THE SURFACE NUMBERS AFTER MERGING
*----
      IF(IMPX.GT.2) THEN
         WRITE (6,620)
         MIN8=MIN(8,NCOUR)
         WRITE (6,650) ('---------------',I=1,MIN8)
         DO 410 IBLK=1,NMBLK
         IKK=INUM(IBLK)
         WRITE (6,630) IBLK,IKK,IGEN(IKK),POURCE(IKK)
         I1=IORI(IBLK)
         IF((MULTC.EQ.4).AND.(NCOUR.EQ.18)) THEN
            DO 380 I=1,18
            ISMZ(I)=3*ISM(1+(I-1)/3,I1)+MOD(I-1,3)-2
            DIRZ(I)=DIRR(1+(I-1)/3)
380      CONTINUE
         ELSE
            DO 400 I=1,NCOUR
            ISMZ(I)=ISM(I,I1)
            DIRZ(I)=DIRR(I)
400      CONTINUE
         ENDIF
         IT0=NCOUR*(IBLK-1)
         IT1=NCOUR*(IKK-1)
         WRITE (6,660) (DIRZ(I),I=1,MIN8)
         WRITE (6,635) (MIX(IT1+ISMZ(I)),IFR(IT0+ISMZ(I)),I=1,MIN8)
         WRITE (6,640) (ALB(IT0+ISMZ(I)),I=1,MIN8)
         WRITE (6,645) (DVX(IT1+ISMZ(I)),I=1,MIN8)
         IF(NCOUR.EQ.18) THEN
            WRITE (6,660) (DIRZ(I),I=9,16)
            WRITE (6,635) (MIX(IT1+ISMZ(I)),IFR(IT0+ISMZ(I)),I=9,16)
            WRITE (6,640) (ALB(IT0+ISMZ(I)),I=9,16)
            WRITE (6,645) (DVX(IT1+ISMZ(I)),I=9,16)
            WRITE (6,660) (DIRZ(I),I=17,18)
            WRITE (6,635) (MIX(IT1+ISMZ(I)),IFR(IT0+ISMZ(I)),I=17,18)
            WRITE (6,640) (ALB(IT0+ISMZ(I)),I=17,18)
            WRITE (6,645) (DVX(IT1+ISMZ(I)),I=17,18)
         ENDIF
         WRITE (6,650) ('---------------',I=1,MIN8)
410      CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(GG3,JF2)
      RETURN
*
570   WRITE (6,620)
      WRITE (6,650) ('---------------',I=1,NCOUR)
      DO 580 IBLK=1,NMBLK
      IKK=INUM(IBLK)
      WRITE (6,630) IBLK,IKK,IGEN(IKK),POURCE(IKK)
      I1=IORI(IBLK)
      IT0=NCOUR*(IBLK-1)
      IT1=NCOUR*(IKK-1)
      WRITE (6,660) (DIRR(I),I=1,NCOUR)
      WRITE (6,635) (MIX(IT1+ISM(I,I1)),IFR(IT0+ISM(I,I1)),I=1,NCOUR)
      WRITE (6,640) (ALB(IT0+ISM(I,I1)),I=1,NCOUR)
      WRITE (6,645) (DVX(IT1+ISM(I,I1)),I=1,NCOUR)
      WRITE (6,650) ('---------------',I=1,NCOUR)
580   CONTINUE
      CALL XABORT('NUMERH: FAILURE OF THE SURFACE RENUMBERING ALGORITH'
     1 //'M.')
*
610   FORMAT(9H HEXAGON=,I5,5X,11HNEIGHBOURS=,6I5,5X,7HWEIGHT=,F7.4)
620   FORMAT (///31H SURFACE NUMBERS AFTER MERGING./)
630   FORMAT (7H BLOCK=,I5,5X,13HMERGED BLOCK=,I5,5X,12HGENERATING B,
     1 5HLOCK=,I5,5X,7HPOURCE=,1P,E12.4)
635   FORMAT (8H IN/OUT:,8(I5,2H /,I5,3H  I))
640   FORMAT (8H ALBEDO:,1P,8(E12.4,3H  I))
645   FORMAT (8H    DVX:,1P,8(E12.4,3H  I))
650   FORMAT (8H -------,8(A15))
660   FORMAT (/8H   SIDE:,8(A9,5X,1HI))
680   FORMAT (49HNUMERH: GEOMETRICAL RECIPROCITY CONDITION IS VIOL,
     1 10HATED (SIDE,A3,20H OF GENERATING BLOCK,I5,8H VS SIDE,A3,
     2 20H OF GENERATING BLOCK,I5,3H ).)
      END
