*DECK EDIGEO
      SUBROUTINE EDIGEO(MAXPTS,MAXMER,IPGEO1,IPGEO2,IPRINT,NREGIO,IEUR,
     1 NMERGE,IMERGE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover the geometry information from LCM (stored by operator GEO) and
* compute a corresponding macro-geometry (i.e. a pure geometry that can
* be treated by Bivac or Trivac or a two-level Eurydice geometry)
* together with NMERGE and array IMERGE.
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
* MAXPTS  allocated storage for arrays of dimension NMCEL.
* MAXMER  allocated storage for arrays of dimension NMERGE.
* IPGEO1  pointer to the original geometry (L_GEOM signature).
* IPGEO2  pointer to the macro-geometry (L_GEOM signature).
* IPRINT  print flag (equal to 0 for no print).
* NREGIO  number of regions.
* IEUR    type of tracking operator for the macro-geometry:
*         =1: 2-level macro-geometry with MERGE;
*         =2: 2-level NXT-type macro-geometry without MERGE;
*         =3: 1-level TRIVAC-type macro-geometry without MERGE;
*         =4: 2-level macro-geometry with MERGE and unfolding.
*
*Parameters: input/output
* NMERGE  macro-calculation merging flag (=1 to indicate a homogeneous
*         macro-calculation) on input and.
*         number of merged regions (equal to the number of physical
*         cells in Eurydice) at output.
*
*Parameters: output
* IMERGE  index of merged regions.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR), TARGET :: IPGEO1
      TYPE(C_PTR)  IPGEO2
      INTEGER      MAXPTS,MAXMER,IPRINT,NREGIO,IEUR,NMERGE,
     >             IMERGE(NREGIO)
*----
*  LOCAL VARIABLES
*----
      PARAMETER   (NSTATE=40)
      CHARACTER    HSMG*131,TYPE(0:9)*16,HSIGN*12,TEXT12*12
      LOGICAL      ILK,LL1,LL2,LHOM,LBIHET
      INTEGER      NCODE(6),ICODE(6),ISTATE(NSTATE),JSTATE(NSTATE)
      REAL         ZCODE(6),ALBEDO(6),XXX(2),YYY(2)
      LOGICAL      LDIAG,LSYMX,LSYMY
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INUM,IGEN,IORI,IMAT,NMC,IBI,
     1 LSECT,NMCR,NRINFO,NRODS,KEYMRG,MATALB,NRODR,NXRS,NXRI,NS,IMILIE,
     2 IHCEL
      REAL, ALLOCATABLE, DIMENSION(:) :: XX,YY,RAYRE,RAN,VS,RODS,RODR,
     1 MESHX,MESHY
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FRACT
      TYPE(C_PTR), POINTER :: IPGEOD,IPGEOX,IPGEOY
*----
*  DATA STATEMENTS
*----
      DATA TYPE/'VIRTUAL','HOMOGENEOUS','CARTESIAN 1-D','TUBE 1-D',
     1 'SPHERE 1-D','CARTESIAN 2-D','TUBE 2-D','CARTESIAN 3-D',
     2 'HEXAGONAL 2-D','HEXAGONAL 3-D'/
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(INUM(2*MAXPTS),IGEN(MAXPTS),IORI(2*MAXPTS),
     1 IMAT(2*MAXPTS),NMC(MAXPTS+1),IBI(NREGIO))
      ALLOCATE(XX(MAXPTS),YY(MAXPTS))
*----
*  RECOVER STATE VECTOR
*----
      ICLS=-99
      LR=-99
      LHOM=(NMERGE.EQ.1)
      CALL LCMGET(IPGEO1,'STATE-VECTOR',ISTATE)
      LBIHET=(ISTATE(12).EQ.1)
*----
*  UNFOLD GEOMETRY
*----
      IF(IEUR.EQ.4) THEN
         IF((ISTATE(1).NE.5).OR.(ISTATE(8).NE.1)) THEN
            CALL XABORT('EDIGEO: UNFOLDING NOT IMPLEMENTED FOR THIS TY'
     1      //'PE OF GEOMETRY.')
         ENDIF
         CALL LCMGET(IPGEO1,'NCODE',NCODE)
         LDIAG=(NCODE(1).EQ.3).AND.(NCODE(4).EQ.3)
         IF(LDIAG) THEN
            ALLOCATE(IPGEOD)
            CALL LCMOP(IPGEOD,'GEOM_DIAG',0,1,0)
            CALL EDIUNF(IPGEO1,IPGEOD,'DIAG')
         ELSE
            IPGEOD=>IPGEO1
         ENDIF
         CALL LCMGET(IPGEOD,'NCODE',NCODE)
         LSYMX=(NCODE(3).EQ.5).AND.(NCODE(1).NE.3)
         IF(LSYMX) THEN
            ALLOCATE(IPGEOX)
            CALL LCMOP(IPGEOX,'GEOM_SYMX',0,1,0)
            CALL EDIUNF(IPGEOD,IPGEOX,'SYMX')
         ELSE
            IPGEOX=>IPGEOD
         ENDIF
         CALL LCMGET(IPGEOX,'NCODE',NCODE)
         LSYMY=(NCODE(1).EQ.5).AND.(NCODE(4).NE.3)
         IF(LSYMY) THEN
            ALLOCATE(IPGEOY)
            CALL LCMOP(IPGEOY,'GEOM_SYMY',0,1,0)
            CALL EDIUNF(IPGEOX,IPGEOY,'SYMY')
         ELSE
            IPGEOY=>IPGEOX
         ENDIF
         ICLS=0
         LR=0
         ALLOCATE(RAYRE(MAXPTS),LSECT(MAXPTS),NMCR(MAXPTS))
         CALL READEU(MAXPTS,MAXPTS,IPGEOY,IR2,IBI,ILK,NMCEL,NMERGE,
     1   NGEN,INUM,IGEN,IPAS,LX,LY,XX,YY,LSECT,RAYRE,NMC,NMCR,IORI,
     2   NCODE,ZCODE,IHEX,IPRINT-5)
         DEALLOCATE(RAYRE,LSECT,NMCR)
         IF(LDIAG) THEN
            CALL LCMCL(IPGEOD,2)
            DEALLOCATE(IPGEOD)
         ENDIF
         IF(LSYMX) THEN
            CALL LCMCL(IPGEOX,2)
            DEALLOCATE(IPGEOX)
         ENDIF
         IF(LSYMY) THEN
            CALL LCMCL(IPGEOY,2)
            DEALLOCATE(IPGEOY)
         ENDIF
         IEUR=3
         GO TO 25
      ENDIF
*----
*  RECOVER ASSEMBLY INFORMATION FROM THE ORIGINAL GEOMETRY
*----
      IF( ((ISTATE(1).EQ.5).AND.(ISTATE(8).EQ.1)).OR.
     1     ((ISTATE(1).EQ.8).AND.(ISTATE(8).EQ.1)).OR.
     2     ((ISTATE(1).EQ.20).AND.(ISTATE(13).EQ.0)).OR.
     3     ((ISTATE(1).EQ.24).AND.(ISTATE(13).EQ.0)) ) THEN
         ICLS=0
         LR=0
         ALLOCATE(RAYRE(MAXPTS),LSECT(MAXPTS),NMCR(MAXPTS))
         CALL READEU(MAXPTS,MAXPTS,IPGEO1,IR2,IBI,ILK,NMCEL,NMERGE,
     1   NGEN,INUM,IGEN,IPAS,LX,LY,XX,YY,LSECT,RAYRE,NMC,NMCR,IORI,
     2   NCODE,ZCODE,IHEX,IPRINT-5)
         DEALLOCATE(RAYRE,LSECT,NMCR)
      ELSE IF( ((ISTATE(1).EQ.3).OR.(ISTATE(1).EQ.20).OR.
     1          (ISTATE(1).EQ.24)).AND.(ISTATE(13).GE.1)) THEN
        IROT=3
        LX=0
        LY=0
        IF(ISTATE(1).EQ.3) THEN
          NSURW=1
          ICLS=1
        ELSE IF(ISTATE(1).EQ.20) THEN
          NSURW=4
          ICLS=2
          LX=1
          LY=1
        ELSE IF(ISTATE(1).EQ.24) THEN
          NSURW=6
          ICLS=3
        ENDIF
        IAPP=1
        MAXJ=1
        CALL XCGDIM(IPGEO1,MAXPTS,NSURW,IROT,IAPP,MAXJ,NVOLW,
     1              NANW,MNANW,NRTW,MSRODW,MARODW,NSURFW)
*----
*  ALLOCATE BLOCK FOR READING GEOMETRY INFORMATION
*----
        ALLOCATE(NRINFO(2*NANW),NRODS(3*NRTW),KEYMRG(NSURW+NVOLW+1),
     1  MATALB(NSURW+NVOLW+1),NRODR(NRTW),NXRS(NRTW),NXRI(NRTW*NANW))
        ALLOCATE(RAN(NANW),VS(NSURW+NVOLW+1),RODS(2*NRTW),
     1  RODR(MSRODW*NRTW))
*
        CALL XCGGEO(IPGEO1,IROT,NSURW,NVOLW,NANW,MNANW,NRTW,MSRODW,
     1  IPRINT,ILK,NMATW,RAN,NRODS,RODS,NRODR,RODR,NRINFO,MATALB,VS,
     2  COTE,RADMIN,NCODE,ICODE,ZCODE,ALBEDO,KEYMRG,NXRS,NXRI)
*
        IR2=0
        DO 10 I=1,NVOLW
        IBI(I)=MATALB((NSURW+1)+I)
        IR2=MAX(IR2,IBI(I))
   10   CONTINUE
        DEALLOCATE(RODR,RODS,VS)
        DEALLOCATE(NXRI,NXRS,NRODR,MATALB,KEYMRG,NRODS)
*----
*  STORE INFORMATION FOR NEW GEOMETRY
*----
        IPAS=NVOLW
        NMCEL=NANW
        NMERGE=NANW
        IHEX=0
        XX(1)=0.0
        NMC(1)=0
        DO 20 II=1,NANW
          XX(II+1)=RAN(II)
          INUM(II)=II
          IGEN(II)=II
          NMC(II+1)=NRINFO(2*(II-1)+1)
   20   CONTINUE
        IF(NSURW.EQ.1) THEN
          LR=NANW
        ELSE IF(NSURW.EQ.4) THEN
          LR=NANW-1
          XXX(1)=0.0
          XXX(2)=XX(NANW+1)
          YYY(1)=0.0
          YYY(2)=COTE
        ELSE
          LR=NANW-1
          SIDE=XX(NANW+1)
        ENDIF
*
        DEALLOCATE(RAN,NRINFO)
      ELSE
         CALL XABORT('EDIGEO: INVALID TYPE OF INPUT GEOMETRY.')
      ENDIF
*----
*  COMPUTE IMERGE ARRAY
*----
   25 IF(LHOM) THEN
*        HOMOGENEOUS MACRO-CALCULATION.
         DO 30 I=1,NREGIO
         IMERGE(I)=1
   30    CONTINUE
         NMERGE=1
         DO 40 I=1,NMCEL
         IMAT(I)=1
   40    CONTINUE
      ELSE
*        CELL-WISE HETEROGENEOUS MACRO-CALCULATION.
         IPAS2=0
         DO 60 IKK=1,NMERGE
         IKG=IGEN(IKK)
         I2=NMC(IKG+1)-NMC(IKG)
         DO 50 I=1,I2
         IMERGE(IPAS2+I)=IKK
   50    CONTINUE
         IPAS2=IPAS2+I2
   60    CONTINUE
         IF(IPAS2.NE.IPAS) THEN
            CALL XABORT('EDIGEO: INCONSISTENT ARRAY NMC.')
         ENDIF
         DO 70 I=1,NMCEL
         IMAT(I)=INUM(I)
   70    CONTINUE
      ENDIF
*----
*  IF DOUBLE HETEROGENEITY IS PRESENT IN THE REFERENCE GEOMETRY
*----
      IF(LBIHET) THEN
         CALL LCMSIX(IPGEO1,'BIHET',1)
         CALL LCMGET(IPGEO1,'STATE-VECTOR',ISTATE)
         NG=ISTATE(1)
         NMILG=ISTATE(3)
         ALLOCATE(NS(NG),IMILIE(NMILG),FRACT(NG,NMILG))
         CALL LCMGET(IPGEO1,'NS',NS)
         CALL LCMGET(IPGEO1,'MILIE',IMILIE)
         CALL LCMGET(IPGEO1,'FRACT',FRACT)
         CALL LCMSIX(IPGEO1,' ',2)
         IND1=IR2
         DO 75 I=1,NMILG
         IF(IMILIE(I).GT.IR2) THEN
            WRITE (HSMG,500) IMILIE(I),IR2
            CALL XABORT(HSMG)
         ENDIF
         IND1=MIN(IND1,IMILIE(I))
   75    CONTINUE
         IPAS2=IPAS
         DO 110 IKK=1,IPAS2
         IF(IBI(IKK).GE.IND1) THEN
            IND=0
            DO 80 I=1,NMILG
            IF(IMILIE(I).EQ.IBI(IKK)) IND=I
   80       CONTINUE
            IF(IND.EQ.0) THEN
               WRITE(HSMG,'(29HEDIGEO: A COMPOSITE MIXTURE (,I5,4H) IS,
     1         13H NOT DEFINED.)') IBI(IKK)
               CALL XABORT(HSMG)
            ENDIF
            DO 100 J=1,NG
            IF(FRACT(J,IND).GT.0.00001) THEN
               DO 90 K=1,NS(J)
               IPAS=IPAS+1
               IMERGE(IPAS)=IMERGE(IKK)
   90          CONTINUE
             ENDIF
  100        CONTINUE
         ENDIF
  110    CONTINUE
         DEALLOCATE(FRACT,IMILIE,NS)
      ENDIF
      IF(IPAS.NE.NREGIO) THEN
         WRITE(HSMG,'(45HEDIGEO: INCONSISTENT NUMBER OF REGIONS. IPAS=,
     1   I6,8H NREGIO=,I6)') IPAS,NREGIO
         CALL XABORT(HSMG)
      ELSE IF(NMERGE.GT.MAXMER) THEN
         WRITE(HSMG,'(37HEDIGEO: INSUFFICIENT STORAGE (MAXMER=,I5,
     1   47H). THE MERGE OPTION SHOULD NOT BE USED IN EDI:.)') MAXMER
         CALL XABORT(HSMG)
      ENDIF
      IF(IPRINT.GT.0) THEN
         WRITE (6,'(/23H EDIGEO: MERGING INDEX:/(1X,14I5))') (IMERGE(I),
     1   I=1,NREGIO)
      ENDIF
*----
*  COMPUTE THE MACRO-GEOMETRY
*----
      CALL XDISET(ISTATE,NSTATE,0)
      CALL XDISET(ICODE,6,0)
      CALL LCMPUT(IPGEO2,'NCODE',6,1,NCODE)
      CALL LCMPUT(IPGEO2,'ZCODE',6,2,ZCODE)
      CALL LCMPUT(IPGEO2,'ICODE',6,1,ICODE)
      IF((ICLS.EQ.0).AND.(IEUR.LE.2)) THEN
*       TWO-LEVEL MACRO-GEOMETRY.
        ALLOCATE(IHCEL(3*NMERGE))
        CALL XDISET(JSTATE,NSTATE,0)
        JSTATE(6)=1
        HSIGN='L_GEOM'
        CALL LCMPTC(IPGEO2,'SIGNATURE',12,1,HSIGN)
        IF(IHEX.EQ.0) THEN
          ISTATE(1)=5
          JSTATE(1)=20
          JSTATE(3)=1
          JSTATE(4)=1
          DO 120 IKK=1,NMERGE
          WRITE(TEXT12,'(4HCELL,I8.8)') IKK
          READ(TEXT12,'(3A4)') (IHCEL((IKK-1)*3+I0),I0=1,3)
          CALL LCMSIX(IPGEO2,TEXT12,1)
          CALL LCMPUT(IPGEO2,'MIX',1,1,IKK)
          XXX(1)=0.0
          XXX(2)=XX(IGEN(IKK))
          YYY(1)=0.0
          YYY(2)=YY(IGEN(IKK))
          CALL LCMPUT(IPGEO2,'MESHX',2,2,XXX)
          CALL LCMPUT(IPGEO2,'MESHY',2,2,YYY)
          JSTATE(7)=IKK
          CALL LCMPUT(IPGEO2,'STATE-VECTOR',NSTATE,1,JSTATE)
          CALL LCMSIX(IPGEO2,' ',2)
  120     CONTINUE
        ELSE
          ISTATE(1)=8
          JSTATE(1)=24
          JSTATE(3)=1
          DO 130 IKK=1,NMERGE
          WRITE(TEXT12,'(4HCELL,I8.8)') IKK
          READ(TEXT12,'(3A4)') (IHCEL((IKK-1)*3+I0),I0=1,3)
          CALL LCMSIX(IPGEO2,TEXT12,1)
          CALL LCMPUT(IPGEO2,'MIX',1,1,IKK)
          CALL LCMPUT(IPGEO2,'SIDE',1,2,XX(1))
          JSTATE(7)=IKK
          CALL LCMPUT(IPGEO2,'STATE-VECTOR',NSTATE,1,JSTATE)
          CALL LCMSIX(IPGEO2,' ',2)
  130     CONTINUE
          CALL LCMPUT(IPGEO2,'IHEX',1,1,IHEX)
        ENDIF
        CALL LCMPUT(IPGEO2,'CELL',3*NMERGE,3,IHCEL)
        IF(IEUR.EQ.1) CALL LCMPUT(IPGEO2,'MERGE',NMCEL,1,INUM)
        CALL LCMPUT(IPGEO2,'TURN',NMCEL,1,IORI)
        DO 140 IBLK=1,NMCEL
        INUM(IBLK)=-INUM(IBLK)
  140   CONTINUE
        CALL LCMPUT(IPGEO2,'MIX',NMCEL,1,INUM)
        DEALLOCATE(IHCEL)
        ISTATE(8)=1
        ISTATE(9)=NMERGE
        ISTATE(10)=1
      ELSE IF(ICLS.EQ.0) THEN
*       ONE-LEVEL (TRIVAC-TYPE) MACRO-GEOMETRY.
        CALL LCMPUT(IPGEO2,'MIX',NMCEL,1,IMAT)
        IF(IHEX.EQ.0) THEN
          ISTATE(1)=5
          ALLOCATE(MESHX(LX+1),MESHY(LY+1))
          MESHX(1)=0.0
          MESHY(1)=0.0
          LL1=((NCODE(2).EQ.3).AND.(NCODE(3).EQ.3))
          LL2=((NCODE(1).EQ.3).AND.(NCODE(4).EQ.3))
          IBLK=0
          DO 155 K1=1,LY
            LXM=1
            LXP=LX
            IF(LL1) LXP=K1
            IF(LL2) LXM=K1
            DO 150 K2=LXM,LXP
              IBLK=IBLK+1
              IKK=INUM(IBLK)
              IF(IKK.EQ.0) GO TO 150
              IF(MOD(IORI(IBLK)-1,2).EQ.0) THEN
                MESHX(K2+1)=MESHX(K2)+XX(IGEN(IKK))
                MESHY(K1+1)=MESHY(K1)+YY(IGEN(IKK))
              ELSE
                MESHX(K2+1)=MESHX(K2)+YY(IGEN(IKK))
                MESHY(K1+1)=MESHY(K1)+XX(IGEN(IKK))
              ENDIF
  150       CONTINUE
  155     CONTINUE
          CALL LCMPUT(IPGEO2,'MESHX',LX+1,2,MESHX)
          CALL LCMPUT(IPGEO2,'MESHY',LY+1,2,MESHY)
          DEALLOCATE(MESHY,MESHX)
        ELSE
          ISTATE(1)=8
          IF(IHEX.EQ.10) THEN
            IHEX=1
          ELSE IF(IHEX.EQ.11) THEN
            IHEX=2
          ENDIF
          CALL LCMPUT(IPGEO2,'IHEX',1,1,IHEX)
          CALL LCMPUT(IPGEO2,'SIDE',1,2,XX(1))
        ENDIF
      ELSE
        CALL LCMPUT(IPGEO2,'MIX',NMCEL,1,IMAT)
        CALL LCMPUT(IPGEO2,'RADIUS',LR+1,2,XX)
        IF(ICLS.EQ.1) THEN
          ISTATE(1)=3
        ELSE IF(ICLS.EQ.2) THEN
          ISTATE(1)=20
          CALL LCMPUT(IPGEO2,'MESHX',LX+1,2,XXX)
          CALL LCMPUT(IPGEO2,'MESHY',LY+1,2,YYY)
        ELSE
          ISTATE(1)=24
          IHEX=9
          CALL LCMPUT(IPGEO2,'IHEX',1,1,IHEX)
          CALL LCMPUT(IPGEO2,'SIDE',1,2,SIDE)
        ENDIF
      ENDIF
      ISTATE(2)=LR
      ISTATE(3)=LX
      ISTATE(4)=LY
      ISTATE(6)=NMCEL
      ISTATE(7)=NMERGE
      CALL LCMPUT(IPGEO2,'STATE-VECTOR',NSTATE,1,ISTATE)
      HSIGN='L_GEOM'
      CALL LCMPTC(IPGEO2,'SIGNATURE',12,1,HSIGN)
      IF(IPRINT.GT.0) THEN
         HSIGN='MACRO$GEO'
         WRITE (6,510) HSIGN,1,TYPE(ISTATE(1))
         IF(IPRINT.GT.1) THEN
            WRITE (6,520) ISTATE(1),TYPE(ISTATE(1)),(ISTATE(I),I=2,7)
         ENDIF
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(YY,XX)
      DEALLOCATE(IBI,NMC,IMAT,IORI,IGEN,INUM)
      RETURN
*
  500 FORMAT (34HEDIGEO: THE INPUT MIXTURE NUMBER (,I4,12H) IS GREATER,
     1 10H THAN IR (,I4,2H ))
  510 FORMAT(/39H EDIGEO: CREATION OF A GEOMETRY NAMED ',A12,7H' ON LE,
     1 3HVEL,I3,11H WITH TYPE ,A16,1H.)
  520 FORMAT(/31H STATE VECTOR (MACRO GEOMETRY):/
     1 7H ITYPE ,I6, 4H   (,A16,1H)/
     2 7H LR    ,I6,20H   (NUMBER OF TUBES)/
     3 7H LX    ,I6,22H   (X-DIMENSION INDEX)/
     4 7H LY    ,I6,22H   (Y-DIMENSION INDEX)/
     5 7H LZ    ,I6,22H   (Z-DIMENSION INDEX)/
     6 7H LREG  ,I6,22H   (NUMBER OF REGIONS)/
     7 7H MAXMIX,I6,25H   (MAX. NB. OF MIXTURES)/1X,60(1H-)/)
      END
