*DECK READEU
      SUBROUTINE READEU (MAXPTS,MAXCEL,IPGEOM,IR,MAT,ILK,NMCEL,NMERGE,
     1 NGEN,INUM,IGEN,NMBLK,LX,LY,XX,YY,LSECT,RAYRE,NMC,NMCR,IORI,NCODE,
     2 ZCODE,IHEX,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover input data for the description of a 2-D assembly (Eurydice-2).
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
* MAXPTS  allocated storage for arrays of dimension NREG.
* MAXCEL  allocated storage for arrays of dimension NMCEL, NMERGE or
*         NGEN.
* IPGEOM  pointer to the geometry LCM object (L_GEOM signature).
* IR      number of mixtures.
* MAT     index-number of the mixture type assigned to each volume.
* ILK     leakage flag (ILK=.true. if neutron leakage through external
*         boundary is present).
* NMCEL   total number of cells in the domain.
* NMERGE  total number of merged cells for which specific values
*         of the neutron flux and reactions rates are required.
*         Many cells with different position in the domain can
*         be merged before the neutron flux calculation if they
*         own the same generating cell (NMERGE.le.NMCEL).
* NGEN    total number of generating cells. A generating cell is
*         defined by its material and dimensions, irrespective of
*         its position in the domain (NGEN.le.NMERGE).
* INUM    index-number of the merged cell associated to each cell.
* IGEN    index-number of the generating cell associated with each
*         merged cell.
* NMBLK   total number of volumes in all the merged cells.
* LX      number of cells along the X-axis.
* LY      number of cells along the Y-axis.
* XX      X-thickness of the generating cells or side of the hexagons.
* YY      Y-thickness of the generating cells.
* LSECT   type of sectorization:
*         =0     no sectorization / specialized treatment;
*         =-999  no sectorization / processed as a sectorized cell;
*         =-101  X-type sectorization of the coolant;
*         =-1    X-type sectorization of the cell;
*         =101   +-type sectorization of the coolant;
*         =1     +-type sectorization of the cell;
*         =102   + and X-type sectorization of the coolant;
*         =2     + and X-type sectorization of the cell.
* RAYRE   radius of the tubes in the generating cells.
* NMC     offsets of the first zone index in each generating cell.
* NMCR    offsets of the first radius index in each generating cell.
* IORI    orientation of the cells.
* NCODE   boundary condition relative to each side of the domain:
*         =0 not used;   =1 VOID;     =2 REFL;
*         =3 DIAG;       =4: TRAN     =5: SYME.
* ZCODE   albedo relative to each side of the domain.
* IHEX    type of symmetry for hexagonal geometry:
*         =0 Cartesian geometry;
*         =1 S30;   =2 SA60;   =3 SB60;   =4 S90;   =5 R120;
*         =6 R180;  =7 SA180;  =8 SB180;  =9 COMPLETE.
* IMPX    print flag (equal to 0 for no print).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPGEOM
      INTEGER MAXPTS,MAXCEL,IR,MAT(MAXPTS),NMCEL,NMERGE,NGEN,
     1 INUM(MAXCEL),IGEN(MAXCEL),NMBLK,LX,LY,LSECT(MAXCEL),
     2 NMC(MAXCEL+1),NMCR(MAXCEL+1),IORI(MAXCEL),NCODE(6),IHEX,IMPX
      REAL XX(MAXCEL),YY(MAXCEL),RAYRE(MAXPTS),ZCODE(6)
      LOGICAL ILK
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40)
      LOGICAL LL1,LL2,EMPTY,LCM
      CHARACTER GEONAM*12,TEXT12*12,HSMG*131
      INTEGER ISTATE(NSTATE),ISTAT2(NSTATE),ICODE(6)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MILIE,NBREG,MILIEU,CELL
      REAL, ALLOCATABLE, DIMENSION(:) :: RAYON,XXX,YYY
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(MILIE(MAXPTS),NBREG(MAXPTS),MILIEU(MAXPTS),
     1 CELL(3*MAXPTS))
      ALLOCATE(RAYON(MAXPTS),XXX(MAXPTS+1),YYY(MAXPTS+1))
*
      CALL LCMINF(IPGEOM,GEONAM,TEXT12,EMPTY,ILONG,LCM)
      IHEX=0
      CALL LCMGET(IPGEOM,'STATE-VECTOR',ISTATE)
      IF((ISTATE(1).EQ.8).OR.(ISTATE(1).EQ.24)) THEN
         CALL LCMGET(IPGEOM,'IHEX',IHEX)
      ENDIF
      IF((ISTATE(5).NE.0).OR.((ISTATE(11).NE.0).AND.(ISTATE(1).EQ.5))
     1 .OR.((ISTATE(11).NE.0).AND.(ISTATE(1).EQ.8)).OR.(ISTATE(13).NE.0)
     2 ) CALL XABORT('READEU: UNABLE TO PROCESS THE GEOMETRY.')
      LX=ISTATE(3)
      LY=ISTATE(4)
      IF(LX.GT.MAXPTS) CALL XABORT('READEU: INSUFFICIENT STORAGE(1).')
      IF(LY.GT.MAXPTS) CALL XABORT('READEU: INSUFFICIENT STORAGE(2).')
      IF(ISTATE(6).GT.MAXPTS) CALL XABORT('READEU: INSUFFICIENT STORA'
     1 //'GE(3).')
*----
*  RECOVER THE BOUNDARY CONDITIONS
*----
      CALL LCMGET(IPGEOM,'NCODE',NCODE)
      CALL LCMGET(IPGEOM,'ZCODE',ZCODE)
      CALL LCMGET(IPGEOM,'ICODE',ICODE)
      I2=0
      DO 10 IC=1,4
      IF(ICODE(IC).NE.0) THEN
         CALL XABORT('READEU: MACROLIB DEFINED ALBEDOS ARE NOT IMPLEMEN'
     1   //'TED.')
      ENDIF
      IF(NCODE(IC).EQ.10) NCODE(IC)=2
      IF(NCODE(IC).EQ.2) ZCODE(IC)=1.0
      IF(NCODE(IC).EQ.6) NCODE(IC)=1
      IF(NCODE(IC).GE.7) CALL XABORT('READEU: INVALID TYPE OF B.C.')
      IF(NCODE(IC).EQ.3) I2=I2+1
   10 CONTINUE
      IF(NCODE(1).EQ.0) GO TO 550
      LL1=.FALSE.
      LL2=.FALSE.
      IF(IHEX.EQ.0) THEN
         IF((NCODE(2).EQ.0).OR.(NCODE(3).EQ.0).OR.(NCODE(4).EQ.0))
     1   GO TO 550
         NSUPCE=LX*LY
         IF(I2.GT.0) THEN
            IF(I2.NE.2) GO TO 560
            IF(LX.NE.LY) CALL XABORT('READEU: LX=LY WITH A DIAGONAL S'
     1      //'YMMETRY.')
            NSUPCE=(LX+1)*LX/2
            LL1=((NCODE(2).EQ.3).AND.(NCODE(3).EQ.3))
            LL2=((NCODE(1).EQ.3).AND.(NCODE(4).EQ.3))
            IF((.NOT.LL1).AND.(.NOT.LL2)) GO TO 560
         ENDIF
         DO 20 IC=1,4,2
         IF((NCODE(IC).EQ.4).AND.(NCODE(IC+1).NE.4))
     1      CALL XABORT('READEU: THE TRANSLATION CONDITIONS X- TRAN X+'
     2      //' TRAN AND Y- TRAN Y+ TRAN ARE THE ONLY PERMITTED.')
   20    CONTINUE
      ELSE
         NSUPCE=LX
         IF((NCODE(2).NE.0).OR.(NCODE(3).NE.0).OR.(NCODE(4).NE.0))
     1   CALL XABORT('READEU: INVALID TYPE OF HEXAGONAL B.C.')
         IF(NCODE(1).EQ.5) THEN
            IF(IHEX.EQ.1) THEN
               IHEX=10
            ELSE IF(IHEX.EQ.2) THEN
               IHEX=11
            ELSE
               CALL XABORT('READEU: BOUNDARY CONDITION HBC WITH OPTION'
     1         //' SYME CAN ONLY BE USED WITH OPTION S30 OR SA60.')
            ENDIF
         ELSE IF(NCODE(1).GT.2) THEN
            CALL XABORT('READEU: BOUNDARY CONDITION HBC CAN ONLY BE US'
     1      //'ED WITH OPTIONS VOID, REFL, SYME OR ALBE.')
         ENDIF
      ENDIF
*
      NMC(1)=0
      NMCR(1)=0
      IG=0
      IGR=0
      IR=0
      NMERGE=1
      NGEN=1
      IF((ISTATE(1).EQ.5).OR.(ISTATE(1).EQ.8)) THEN
         NMCEL=ISTATE(6)
         IF(NMCEL.GT.MAXPTS) THEN
            WRITE(HSMG,'(36HREADEU: INSUFFICIENT STORAGE. NMCEL=,I8,
     1      8H MAXPTS=,I8)') NMCEL,MAXPTS
            CALL XABORT(HSMG)
         ENDIF
         DO 30 IKK=1,NMCEL
         IGEN(IKK)=0
         INUM(IKK)=IKK
         IORI(IKK)=1
   30    CONTINUE
         IF(ISTATE(8).EQ.1) THEN
*           MIXED GEOMETRY.
            CALL LCMLEN(IPGEOM,'MERGE',ILEN,ITYLCM)
            IF(ILEN.GT.0) CALL LCMGET(IPGEOM,'MERGE',INUM)
            CALL LCMLEN(IPGEOM,'TURN',ILEN,ITYLCM)
            IF(ILEN.GT.0) CALL LCMGET(IPGEOM,'TURN',IORI)
            CALL LCMGET(IPGEOM,'CELL',CELL)
            CALL LCMGET(IPGEOM,'MIX',NBREG)
            DO 36 IKK=1,NMCEL
            IF(-NBREG(IKK).LT.0) CALL XABORT('READEU: GENERATING CELL '
     1      //'EXPECTED.')
            IGEN(INUM(IKK))=-NBREG(IKK)
            NGEN=MAX(NGEN,-NBREG(IKK))
            NMERGE=MAX(NMERGE,INUM(IKK))
            DO 35 JKK=1,NMCEL
            IF(INUM(IKK).EQ.INUM(JKK)) THEN
              IF(NBREG(IKK).NE.NBREG(JKK)) THEN
                WRITE(HSMG,'(38HREADEU: TWO CELLS WITH THE SAME MERGED,
     1          46H NUMBER DO NOT HAVE THE SAME GENERATING CELL (,2I7,
     2          2H).)') IKK,JKK
                CALL XABORT(HSMG)
              ENDIF
            ENDIF
   35       CONTINUE
   36       CONTINUE
            IF(NGEN.GT.ISTATE(9)) CALL XABORT('READEU: INVALID NUMBER'
     1      //' OF SUB GEOMETRIES.')
            DO 70 IKG=1,NGEN
            WRITE (TEXT12(:4),'(A4)') CELL(3*(IKG-1)+1)
            WRITE (TEXT12(5:8),'(A4)') CELL(3*(IKG-1)+2)
            WRITE (TEXT12(9:),'(A4)') CELL(3*(IKG-1)+3)
            CALL LCMLEN(IPGEOM,TEXT12,ILEN,ITYLCM)
            IF((ILEN.EQ.0).OR.(ITYLCM.NE.0)) CALL XABORT('READEU: SUB'
     1      //' GEOMETRY '//TEXT12//' IS MISSING FROM L_GEOM.')
            CALL LCMSIX(IPGEOM,TEXT12,1)
            CALL LCMGET(IPGEOM,'STATE-VECTOR',ISTAT2)
            IF(ISTAT2(6).GT.MAXPTS) CALL XABORT('READEU: INSUFFICI'
     1      //'ENT STORAGE(4).')
            ISECTO=ISTAT2(14)
            JSECTO=ISTAT2(15)
            IF((IHEX.EQ.0).OR.(IHEX.NE.0)) THEN
               NZONE=ISTAT2(2)+1
               CALL LCMLEN(IPGEOM,'SPLITR',ILENN,ITYLCM)
               IF(ILENN.GT.0) THEN
                  CALL LCMGET(IPGEOM,'SPLITR',NBREG)
               ELSE
                  DO 40 I=1,NZONE
                  NBREG(I)=1
   40             CONTINUE
               ENDIF
               IF(NZONE.GT.1) CALL LCMGET(IPGEOM,'RADIUS',RAYON)
            ELSE IF(((ISTAT2(1).EQ.5).OR.(ISTAT2(1).EQ.8)).AND.
     1      (ISTAT2(6).EQ.1)) THEN
               NZONE=1
               NBREG(1)=1
            ELSE
               CALL XABORT('READEU: INVALID SUB GEOMETRY.')
            ENDIF
            LS1=0
            LS2=0
            NZONES=0
            LSECT(IKG)=ISECTO
            IF((ISECTO.EQ.0).OR.(ISECTO.EQ.-999)) THEN
*              NO SECTORIZATION.
               LS1=1
               LS2=1
               NZONES=NZONE
            ELSE IF((IHEX.EQ.0).AND.(ISECTO.EQ.-1).AND.
     1              (JSECTO.EQ.0)) THEN
*              X-TYPE SECTORIZATION IN CARTESIAN CELL.
               LS1=4
               LS2=4
               NZONES=4*NZONE
            ELSE IF((IHEX.EQ.0).AND.(ISECTO.EQ.-1).AND.
     1              (JSECTO.EQ.NZONE-1)) THEN
*              X-TYPE COOLANT SECTORIZATION IN CARTESIAN CELL.
               LSECT(IKG)=-101
               LS1=1
               LS2=4
               NZONES=NZONE+3
            ELSE IF((IHEX.EQ.0).AND.(ISECTO.EQ.1).AND.
     1              (JSECTO.EQ.0)) THEN
*              +-TYPE SECTORIZATION IN CARTESIAN CELL.
               LS1=4
               LS2=4
               NZONES=4*NZONE
            ELSE IF((IHEX.EQ.0).AND.(ISECTO.EQ.1).AND.
     1              (JSECTO.EQ.NZONE-1)) THEN
*              +-TYPE COOLANT SECTORIZATION IN CARTESIAN CELL.
               LSECT(IKG)=101
               LS1=1
               LS2=4
               NZONES=NZONE+3
            ELSE IF((IHEX.EQ.0).AND.(ISECTO.EQ.2).AND.
     1              (JSECTO.EQ.0)) THEN
*              +-TYPE SECTORIZATION IN CARTESIAN CELL.
               LS1=8
               LS2=8
               NZONES=8*NZONE
            ELSE IF((IHEX.EQ.0).AND.(ISECTO.EQ.2).AND.
     1              (JSECTO.EQ.NZONE-1)) THEN
*              +-TYPE COOLANT SECTORIZATION IN CARTESIAN CELL.
               LSECT(IKG)=102
               LS1=1
               LS2=8
               NZONES=NZONE+7
            ELSE IF((IHEX.GT.0).AND.(ISECTO.EQ.-1).AND.
     1              (JSECTO.EQ.0)) THEN
*              X-TYPE SECTORIZATION IN HEXAGONAL CELL.
               LS1=6
               LS2=6
               NZONES=6*NZONE
            ELSE IF((IHEX.GT.0).AND.(ISECTO.EQ.-1).AND.
     1              (JSECTO.EQ.NZONE-1)) THEN
*              X-TYPE COOLANT SECTORIZATION IN HEXAGONAL CELL.
               LSECT(IKG)=-101
               LS1=1
               LS2=6
               NZONES=NZONE+5
            ELSE
               CALL XABORT('READEU: INVALID TYPE OF SECTORIZATION.')
            ENDIF
            IF(NZONES.GT.MAXPTS) CALL XABORT('READEU: MAXPTS OVERFLOW.')
            CALL LCMLEN(IPGEOM,'MIX',ILENG,ITYLCM)
            IF(ILENG.NE.NZONES) CALL XABORT('READEU: BAD MIX LENGTH.')
            CALL LCMGET(IPGEOM,'MIX',MILIEU)
            IF(IHEX.EQ.0) THEN
               CALL LCMGET(IPGEOM,'MESHX',XXX)
               CALL LCMGET(IPGEOM,'MESHY',YYY)
            ELSE
               CALL LCMGET(IPGEOM,'SIDE',SIDE)
            ENDIF
            CALL LCMSIX(IPGEOM,' ',2)
*
            RJ=0.0
            RAYRE(IGR+1)=0.0
            DO 60 I=1,NZONE-1
            IF(RAYON(I+1).LE.RJ) GO TO 520
            PAS=(RAYON(I+1)-RJ)/REAL(ABS(NBREG(I)))
            IF(NBREG(I).LT.0) PAS=PAS*(RAYON(I+1)+RJ)
            DO 50 J=1,ABS(NBREG(I))
            IGR=IGR+1
            DO 45 ISEC=1,LS1
            IG=IG+1
            IF(IG.GT.MAXPTS) THEN
               WRITE(HSMG,'(28HREADEU: INCREASE MAXPTS FROM,I7,3H TO,
     1         I7,5H (1).)') MAXPTS,IG
               CALL XABORT(HSMG)
            ENDIF
            MILIE(IG)=MILIEU((I-1)*LS1+ISEC)
            IR=MAX(IR,MILIE(IG))
   45       CONTINUE
            IF(NBREG(I).GT.0) THEN
               RJ=RJ+PAS
            ELSE
               RJ=SQRT(RJ*RJ+PAS)
            ENDIF
            RAYRE(IGR+1)=RJ
   50       CONTINUE
            RJ=RAYON(I+1)
   60       CONTINUE
            IGR=IGR+1
            DO 65 ISEC=1,LS2
            IG=IG+1
            IF(IG.GT.MAXPTS) THEN
               WRITE(HSMG,'(28HREADEU: INCREASE MAXPTS FROM,I7,3H TO,
     1         I7,5H (2).)') MAXPTS,IG
               CALL XABORT(HSMG)
            ENDIF
            MILIE(IG)=MILIEU((NZONE-1)*LS1+ISEC)
            IR=MAX(IR,MILIE(IG))
   65       CONTINUE
            IF(IHEX.EQ.0) THEN
               XX(IKG)=XXX(2)-XXX(1)
               YY(IKG)=YYY(2)-YYY(1)
               TEST=2.0*RAYRE(IGR)
               IF(SQRT(XX(IKG)**2+YY(IKG)**2).LE.TEST) GO TO 520
            ELSE
               XX(IKG)=SIDE
               YY(IKG)=0.0
               IF(SIDE.LE.RAYRE(IGR)) GO TO 520
            ENDIF
            NMC(IKG+1)=IG
            NMCR(IKG+1)=IGR
   70       CONTINUE
*           COMPUTE THE MIXTURE NUMBERS IN THE MERGED CELLS.
            NMBLK=0
            DO 90 IKK=1,NMERGE
            IKG=IGEN(IKK)
            IF(IKG.EQ.0) THEN
               WRITE(HSMG,'(14HREADEU: VOLUME,I5,16H NOT DEFINED(1).)')
     1         IKK
               CALL XABORT(HSMG)
            ENDIF
            I1=NMC(IKG)
            I2=NMC(IKG+1)-I1
            IF(NMBLK+I2.GT.MAXPTS) THEN
               WRITE(HSMG,'(28HREADEU: INCREASE MAXPTS FROM,I7,3H TO,
     1         I7,5H (3).)') MAXPTS,NMBLK+I2
               CALL XABORT(HSMG)
            ENDIF
            DO 80 I=1,I2
            MAT(NMBLK+I)=MILIE(I1+I)
   80       CONTINUE
            NMBLK=NMBLK+I2
   90       CONTINUE
         ELSE IF(IHEX.EQ.0) THEN
*           PURE CARTESIAN GEOMETRY.
            CALL LCMGET(IPGEOM,'MIX',MILIEU)
            DO 100 I=1,NMCEL
            MAT(INUM(I))=MILIEU(I)
  100       CONTINUE
            CALL LCMGET(IPGEOM,'MESHX',XXX)
            IF(LL1.OR.LL2) THEN
               CALL LCMGET(IPGEOM,'MESHX',YYY)
            ELSE
               CALL LCMGET(IPGEOM,'MESHY',YYY)
            ENDIF
            NGEN=0
            IBLK=0
            DO 125 K1=1,LY
            LXM=1
            LXP=LX
            IF(LL1) LXP=K1
            IF(LL2) LXM=K1
            DO 120 K2=LXM,LXP
            IBLK=IBLK+1
            IF(MAT(INUM(IBLK)).EQ.0) GO TO 120
            IKK=INUM(IBLK)
            NMERGE=MAX(NMERGE,IKK)
            A=XXX(K2+1)-XXX(K2)
            B=YYY(K1+1)-YYY(K1)
            DO 110 JBLK=1,IBLK-1
            JKG=IGEN(INUM(JBLK))
            IF(MAT(INUM(IBLK)).EQ.MAT(INUM(JBLK))) THEN
               IF((A.EQ.XX(JKG)).AND.(B.EQ.YY(JKG))) THEN
                  IGEN(IKK)=JKG
                  GO TO 120
               ELSE IF((B.EQ.XX(JKG)).AND.(A.EQ.YY(JKG))) THEN
                  IGEN(IKK)=JKG
                  IORI(IBLK)=2
                  GO TO 120
               ENDIF
            ENDIF
  110       CONTINUE
            NGEN=NGEN+1
            IGEN(IKK)=NGEN
            XX(NGEN)=A
            YY(NGEN)=B
            LSECT(NGEN)=0
            NMC(NGEN+1)=NMC(NGEN)+1
            NMCR(NGEN+1)=NMCR(NGEN)+1
            RAYRE(NGEN)=0.0
  120       CONTINUE
  125       CONTINUE
            NMBLK=NMERGE
         ELSE IF(IHEX.GT.0) THEN
*           PURE HEXAGONAL GEOMETRY.
            CALL LCMGET(IPGEOM,'MIX',MILIEU)
            DO 130 I=1,NMCEL
            MAT(INUM(I))=MILIEU(I)
  130       CONTINUE
            CALL LCMGET(IPGEOM,'SIDE',SIDE)
            NGEN=0
            DO 140 IBLK=1,LX
            IF(MAT(INUM(IBLK)).EQ.0) GO TO 140
            IKK=INUM(IBLK)
            NMERGE=MAX(NMERGE,IKK)
            DO 135 JBLK=1,IBLK-1
            IF(MAT(INUM(IBLK)).EQ.MAT(INUM(JBLK))) THEN
               IGEN(IKK)=IGEN(INUM(JBLK))
               GO TO 140
            ENDIF
  135       CONTINUE
            NGEN=NGEN+1
            IGEN(IKK)=NGEN
            XX(NGEN)=SIDE
            YY(NGEN)=0.0
            LSECT(NGEN)=0
            NMC(NGEN+1)=NMC(NGEN)+1
            NMCR(NGEN+1)=NMCR(NGEN)+1
            RAYRE(NGEN)=0.0
  140       CONTINUE
            NMBLK=NMERGE
         ENDIF
      ELSE IF((ISTATE(1).EQ.20).OR.(ISTATE(1).EQ.24)) THEN
         NZONE=ISTATE(2)+1
         ISECTO=ISTATE(14)
         JSECTO=ISTATE(15)
         NMCEL=1
         IGEN(1)=1
         INUM(1)=1
         IORI(1)=1
         CALL LCMLEN(IPGEOM,'SPLITR',ILENN,ITYLCM)
         IF(ILENN.GT.0) THEN
            CALL LCMGET(IPGEOM,'SPLITR',NBREG)
         ELSE
            DO 150 I=1,NZONE
            NBREG(I)=1
  150       CONTINUE
         ENDIF
         IF(NZONE.GT.1) CALL LCMGET(IPGEOM,'RADIUS',RAYON)
         LS1=0
         LS2=0
         NZONES=0
         LSECT(1)=ISECTO
         IF((ISECTO.EQ.0).OR.(ISECTO.EQ.-999)) THEN
*           NO SECTORIZATION.
            LS1=1
            LS2=1
            NZONES=NZONE
         ELSE IF((IHEX.EQ.0).AND.(ISECTO.EQ.-1).AND.(JSECTO.EQ.0)) THEN
*           X-TYPE SECTORIZATION IN CARTESIAN CELL.
            LS1=4
            LS2=4
            NZONES=4*NZONE
         ELSE IF((IHEX.EQ.0).AND.(ISECTO.EQ.-1).AND.
     1           (JSECTO.EQ.NZONE-1)) THEN
*           X-TYPE COOLANT SECTORIZATION IN CARTESIAN CELL.
            LSECT(1)=-101
            LS1=1
            LS2=4
            NZONES=NZONE+3
         ELSE IF((IHEX.EQ.0).AND.(ISECTO.EQ.1).AND.(JSECTO.EQ.0)) THEN
*           +-TYPE SECTORIZATION IN CARTESIAN CELL.
            LS1=4
            LS2=4
            NZONES=4*NZONE
         ELSE IF((IHEX.EQ.0).AND.(ISECTO.EQ.1).AND.
     1           (JSECTO.EQ.NZONE-1)) THEN
*           +-TYPE COOLANT SECTORIZATION IN CARTESIAN CELL.
            LSECT(1)=101
            LS1=1
            LS2=4
            NZONES=NZONE+3
         ELSE IF((IHEX.EQ.0).AND.(ISECTO.EQ.2).AND.(JSECTO.EQ.0)) THEN
*           +-TYPE SECTORIZATION IN CARTESIAN CELL.
            LS1=8
            LS2=8
            NZONES=8*NZONE
         ELSE IF((IHEX.EQ.0).AND.(ISECTO.EQ.2).AND.
     1           (JSECTO.EQ.NZONE-1)) THEN
*           +-TYPE COOLANT SECTORIZATION IN CARTESIAN CELL.
            LSECT(1)=102
            LS1=1
            LS2=8
            NZONES=NZONE+7
         ELSE IF((IHEX.GT.0).AND.(ISECTO.EQ.-1).AND.(JSECTO.EQ.0)) THEN
*           X-TYPE SECTORIZATION IN HEXAGONAL CELL.
            LS1=6
            LS2=6
            NZONES=6*NZONE
         ELSE IF((IHEX.GT.0).AND.(ISECTO.EQ.-1).AND.
     1           (JSECTO.EQ.NZONE-1)) THEN
*           X-TYPE COOLANT SECTORIZATION IN HEXAGONAL CELL.
            LSECT(1)=-101
            LS1=1
            LS2=6
            NZONES=NZONE+5
         ELSE
            CALL XABORT('READEU: INVALID TYPE OF SECTORIZATION(2).')
         ENDIF
         IF(NZONES.GT.MAXPTS) CALL XABORT('READEU: MAXPTS OVERFLOW(2).')
         CALL LCMLEN(IPGEOM,'MIX',ILENG,ITYLCM)
         IF(ILENG.NE.NZONES) CALL XABORT('READEU: BAD MIX LENGTH(2).')
         CALL LCMGET(IPGEOM,'MIX',MILIEU)
         IF(IHEX.EQ.0) THEN
            CALL LCMGET(IPGEOM,'MESHX',XXX)
            CALL LCMGET(IPGEOM,'MESHY',YYY)
         ELSE
            CALL LCMGET(IPGEOM,'SIDE',SIDE)
         ENDIF
*
         RJ=0.0
         RAYRE(1)=0.0
         DO 170 I=1,NZONE-1
         IF(RAYON(I+1).LE.RJ) GO TO 520
         PAS=(RAYON(I+1)-RJ)/REAL(ABS(NBREG(I)))
         IF(NBREG(I).LT.0) PAS=PAS*(RAYON(I+1)+RJ)
         DO 160 J=1,ABS(NBREG(I))
         IGR=IGR+1
         DO 155 ISEC=1,LS1
         IG=IG+1
         IF(IG.GT.MAXPTS) THEN
            WRITE(HSMG,'(28HREADEU: INCREASE MAXPTS FROM,I7,3H TO,
     1      I7,5H (4).)') MAXPTS,IG
            CALL XABORT(HSMG)
         ENDIF
         MAT(IG)=MILIEU((I-1)*LS1+ISEC)
         IR=MAX(IR,MAT(IG))
  155    CONTINUE
         IF(NBREG(I).GT.0) THEN
            RJ=RJ+PAS
         ELSE
            RJ=SQRT(RJ*RJ+PAS)
         ENDIF
         RAYRE(IGR+1)=RJ
  160    CONTINUE
         RJ=RAYON(I+1)
  170    CONTINUE
         IGR=IGR+1
         DO 180 ISEC=1,LS2
         IG=IG+1
         IF(IG.GT.MAXPTS) THEN
            WRITE(HSMG,'(28HREADEU: INCREASE MAXPTS FROM,I7,3H TO,
     1      I7,5H (5).)') MAXPTS,IG
            CALL XABORT(HSMG)
         ENDIF
         MAT(IG)=MILIEU((NZONE-1)*LS1+ISEC)
         IR=MAX(IR,MAT(IG))
  180    CONTINUE
         IF(IHEX.EQ.0) THEN
            XX(1)=XXX(2)-XXX(1)
            YY(1)=YYY(2)-YYY(1)
            TEST=2.0*RAYRE(IGR)
            IF(SQRT(XX(1)**2+YY(1)**2).LE.TEST) GO TO 520
         ELSE
            XX(1)=SIDE
            YY(1)=0.0
            IF(SIDE.LE.RAYRE(IGR)) GO TO 520
         ENDIF
         NMC(2)=IG
         NMCR(2)=IGR
         NMBLK=IG
      ELSE
         CALL XABORT('READEU: INVALID PRIMARY GEOMETRY.')
      ENDIF
      IF(NSUPCE.NE.NMCEL) CALL XABORT('READEU: THE CALCULATED NUMBER O'
     1 //'F CELLS IS INCONSISTENT.')
*
      ILK=((NCODE(1).EQ.1).AND.(ZCODE(1).NE.1.0)).OR.
     1    ((NCODE(2).EQ.1).AND.(ZCODE(2).NE.1.0)).OR.
     2    ((NCODE(3).EQ.1).AND.(ZCODE(3).NE.1.0)).OR.
     3    ((NCODE(4).EQ.1).AND.(ZCODE(4).NE.1.0))
      IF(IMPX.GT.0) THEN
         IF(IHEX.EQ.0) THEN
            WRITE (6,'(/43H CARTESIAN MULTICELL OPTION (EURYDICE-2) BA,
     1      28HSED ON GEOMETRY LOCATED IN '',A12,2H''./)') GEONAM
         ELSE
            WRITE (6,'(/43H HEXAGONAL MULTICELL OPTION (EURYDICE-2) BA,
     1      28HSED ON GEOMETRY LOCATED IN '',A12,2H''./)') GEONAM
         ENDIF
         IF(IHEX.EQ.0) THEN
            WRITE (6,670) LX,LY,MAXPTS,NMBLK,IR
         ELSE
            WRITE (6,680) LX,MAXPTS,NMBLK,IR
         ENDIF
         WRITE (6,630)
         DO 190 IKG=1,NGEN
         IF((ISTATE(1).EQ.20).OR.(ISTATE(1).EQ.24)) THEN
            TEXT12=GEONAM
         ELSE IF(ISTATE(8).EQ.1) THEN
            WRITE (TEXT12(:4),'(A4)') CELL(3*(IKG-1)+1)
            WRITE (TEXT12(5:8),'(A4)') CELL(3*(IKG-1)+2)
            WRITE (TEXT12(9:),'(A4)') CELL(3*(IKG-1)+3)
         ELSE
            WRITE (TEXT12,'(4HCELL,I5)') IKG
         ENDIF
         I1=NMCR(IKG)+1
         I2=NMCR(IKG+1)
         IF(I1.EQ.I2) THEN
            IF(IHEX.EQ.0) THEN
               WRITE (6,660) IKG,TEXT12,XX(IKG),YY(IKG)
            ELSE
               WRITE (6,665) IKG,TEXT12,XX(IKG)
            ENDIF
         ELSE
            WRITE (6,640) IKG,TEXT12,(RAYRE(I),I=I1,I2)
            IF(IHEX.EQ.0) THEN
               WRITE (6,650) XX(IKG),YY(IKG)
            ELSE
               WRITE (6,655) XX(IKG)
            ENDIF
         ENDIF
  190    CONTINUE
         WRITE (6,'(/)')
         IF(.NOT.ILK) WRITE (6,'(17H INFINITE DOMAIN./)')
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(YYY,XXX,RAYON)
      DEALLOCATE(CELL,MILIEU,NBREG,MILIE)
      RETURN
*
  520 CALL XABORT('READEU: RADIUS ARE INCONSISTENTS.')
  550 CALL XABORT('READEU: A BOUNDARY CONDITION IS MISSING.')
  560 CALL XABORT('READEU: THE DIAGONAL CONDITIONS X+ DIAG Y- DIAG AND'
     1 //' X- DIAG Y+ DIAG ARE THE ONLY PERMITTED.')
*
  630 FORMAT (/5X,24HGENERATING CELL   RADIUS)
  640 FORMAT (/1X,I4,2H ',A12,3H'  ,1P,9E12.5/(22X,9E12.5))
  650 FORMAT (23X,3HA =,1P,E12.5,6H   B =,E12.5)
  655 FORMAT (23X,6HSIDE =,1P,E12.5)
  660 FORMAT (/1X,I4,2H ',A12,1H',3X,3HA =,1P,E12.5,6H   B =,E12.5)
  665 FORMAT (/1X,I4,2H ',A12,1H',3X,6HSIDE =,1P,E12.5)
  670 FORMAT (/35H NUMBER OF CELLS ALONG THE X-AXIS =,I4/17X,
     1 18HALONG THE Y-AXIS =,I4,5X,26HAVAILABLE STORAGE MAXPTS =,I4/
     2 27H NUMBER OF MERGED VOLUMES =,I5/
     3 39H NUMBER OF DISTINCT PHYSICAL MIXTURES =,I5/)
  680 FORMAT (/34H NUMBER OF HEXAGONS IN ONE PLANE =,I4,5X,9HAVAILABLE,
     1 17H STORAGE MAXPTS =,I4/27H NUMBER OF MERGED VOLUMES =,I5/
     2 39H NUMBER OF DISTINCT PHYSICAL MIXTURES =,I5/)
      END
