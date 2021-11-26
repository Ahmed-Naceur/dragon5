*DECK SYBEUR
      SUBROUTINE SYBEUR (MAXPTS,MAXCEL,MAXJ,MAXZ,IPGEOM,NREG,IR,MAT,
     1 VOL,ILK,IMPX,IHEX,NCOUR,LMAILI,LMAILR,NMCEL,NMERGE,NGEN,IJAT,
     2 MULTC,IWIGN,IHALT,ILIGN,INORM,IRECT,IQW,IQUAD,XX,YY,LSECT,NMC,
     3 NMCR,RAYRE,MAIL,IZMAIL,RZMAIL,IFR,ALB,SUR,INUM,MIX,DVX,IGEN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover and analysis of the geometry for the description of a 2-D
* assembly (Eurydice-2) for the interface current method.
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
* MAXJ    allocated storage for arrays of dimension IJAS, IJAZ or IJAT.
* MAXZ    allocated storage for arrays of dimension LMAILI or LMAILR.
* IPGEOM  pointer to the geometry LCM object (L_GEOM signature).
* NREG    total number of volumes.
* IR      number of mixtures.
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* ILK     leakage flag (ILK=.TRUE. if neutron leakage through external
*         boundary is present).
* IMPX    print flag (equal to 0 for no print).
* IHEX    type of lattice: =0 Cartesian lattice; .ne.0 hexagonal
*         lattice.
* NCOUR   number of out-currents per cell.
* LMAILI  space required to store the integer tracking information.
* LMAILR  space required to store the real tracking information.
* MULTC   type of multicell approximation.
* IWIGN   type of cylinderization (=1/2/3: Askew/Wigner/Sanchez
*         cylinderization).
* IHALT   halt switch (=1 to stop after calculation of the maximum
*         tracking storage).
* ILIGN   tracking print flag (=1 to print the tracking).
* INORM   track normalization flag (=1 to avoid track normalization).
* IRECT   rectangular flag (=1 to avoid considering symmetries of
*         square cells).
* IQW     equal weight quadrature flag (=1 to use equal weight
*         quadratures in angle and space).
* IQUAD   quadrature parameters.
* NMCEL   total number of cells in the domain.
* IFR     index-number of in-currents.
* ALB     transmission/albedo associated with each in-current.
* SUR     surface associated with each in-current.
* NMERGE  total number of merged cells for which specific values
*         of the neutron flux and reactions rates are required.
*         Many cells with different position in the domain can
*         be merged before the neutron flux calculation if they
*         own the same generating cell (NMERGE.le.NMCEL).
* IJAT    total number of distinct out-currents.
* INUM    index-number of the merged cell associated to each cell.
* MIX     index-number of out-currents.
* DVX     weight associated with each out-current.
*         Note: IFR, ALB, MIX and DVX contains information to rebuild
*         the geometrical 'A' matrix.
* NGEN    total number of generating cells. A generating cell is
*         defined by its material and dimensions, irrespective of
*         its position in the domain (NGEN.le.NMERGE).
* XX      X-thickness of the generating cells.
* YY      Y-thickness of the generating cells.
* LSECT   type of sectorization.
* NMC     offset of the first volume in each generating cell.
* NMCR    offset of the first radius in each generating cell
*         (identical to NMC in no-sectorization cases).
* RAYRE   radius of the tubes in each generating cell.
* MAIL    offset of the first tracking information in each generating
*         cell.
* IZMAIL  integer tracking information.
* RZMAIL  real tracking information.
* IGEN    index-number of the generating cell associated with each
*         merged cell.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPGEOM
      INTEGER MAXPTS,MAXCEL,MAXJ,MAXZ,NREG,IR,MAT(MAXPTS),IMPX,IHEX,
     1 NCOUR,LMAILI,LMAILR,NMCEL,NMERGE,NGEN,IJAT,MULTC,IWIGN,IHALT,
     2 ILIGN,INORM,IRECT,IQW,IQUAD(4),LSECT(MAXCEL),NMC(MAXCEL+1),
     3 NMCR(MAXCEL+1),MAIL(2,MAXCEL),IZMAIL(MAXZ),IFR(MAXJ),
     4 INUM(MAXCEL),MIX(MAXJ),IGEN(MAXCEL)
      REAL VOL(MAXPTS),XX(MAXCEL),YY(MAXCEL),RAYRE(MAXPTS),RZMAIL(MAXZ),
     1 ALB(MAXJ),SUR(MAXJ),DVX(MAXJ)
      LOGICAL ILK
*----
*  LOCAL VARIABLES
*----
      PARAMETER (PI=3.141592654)
      LOGICAL LHEX
      CHARACTER HSMG*131,HDOM*8
      INTEGER NCODE(6),IIT(24),ISMR(6,8),ISMH(8,12),ZZI
      REAL ZCODE(6)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IORI
      REAL, ALLOCATABLE, DIMENSION(:) :: POURCE,ZZ,VOL2,RAD
*----
*  DATA STATEMENTS
*----
      SAVE ISMR,ISMH,IIT
      DATA ISMR/1,2,3,4,0,0,3,4,2,1,0,0,2,1,4,3,0,0,4,3,1,2,0,0,
     1          2,1,3,4,0,0,3,4,1,2,0,0,1,2,4,3,0,0,4,3,2,1,0,0/
      DATA ISMH/1,2,3,4,5,6,0,0,6,1,2,3,4,5,0,0,5,6,1,2,3,4,0,0,
     1          4,5,6,1,2,3,0,0,3,4,5,6,1,2,0,0,2,3,4,5,6,1,0,0,
     2          1,6,5,4,3,2,0,0,2,1,6,5,4,3,0,0,3,2,1,6,5,4,0,0,
     3          4,3,2,1,6,5,0,0,5,4,3,2,1,6,0,0,6,5,4,3,2,1,0,0/
      DATA IIT/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,
     1          24,28,32,64/
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IORI(MAXCEL),POURCE(MAXCEL))
*
      CALL READEU(MAXPTS,MAXCEL,IPGEOM,IR,MAT,ILK,NMCEL,NMERGE,
     1 NGEN,INUM,IGEN,NREG,LX,LY,XX,YY,LSECT,RAYRE,NMC,NMCR,IORI,
     2 NCODE,ZCODE,IHEX,IMPX)
      LHEX=IHEX.NE.0
      IF(MULTC.LE.2) THEN
*        MODIFICATION OF THE RADIUS IN CASES WHERE THEY INTERCEPT
*        THE CELL BOUNDARY. CYLINDERIZED TUBES ARE COMPUTED SO AS
*        TO PRESERVE THE EXACT VOLUMES.
         DO 25 IKG=1,NGEN
         I1=NMC(IKG)+1
         I2=NMC(IKG+1)
         DO 20 I=I1+1,I2
         RJ=RAYRE(I)
         VTOT=PI*RJ*RJ
         IF(LHEX) THEN
            A=XX(IKG)*SQRT(3.0)
            IF(2.0*RJ.GT.A) THEN
               THETA=ACOS(0.5*A/RJ)
               VT=0.5*A*RJ*SIN(THETA)
               VC=RJ*RJ*THETA
               VTOT=VTOT-6.0*(VC-VT)
               RJ=SQRT(VTOT/PI)
            ENDIF
         ELSE
            IF(2.0*RJ.GT.XX(IKG)) THEN
               THETA=ACOS(0.5*XX(IKG)/RJ)
               VT=0.5*XX(IKG)*RJ*SIN(THETA)
               VC=RJ*RJ*THETA
               VTOT=VTOT-2.0*(VC-VT)
            ENDIF
            IF(2.0*RJ.GT.YY(IKG)) THEN
               THETA=ACOS(0.5*YY(IKG)/RJ)
               VT=0.5*YY(IKG)*RJ*SIN(THETA)
               VC=RJ*RJ*THETA
               VTOT=VTOT-2.0*(VC-VT)
            ENDIF
            IF(2.0*RJ.GT.MIN(XX(IKG),YY(IKG))) RJ=SQRT(VTOT/PI)
         ENDIF
         RAYRE(I)=RJ
20       CONTINUE
25       CONTINUE
*
         DO 30 IKG=1,NGEN
            IF (LSECT(IKG).NE.0) CALL XABORT('SYBEUR: SECTORIZATION FO'
     1      //'RBIDDEN.')
30       CONTINUE
*
         IF(.NOT.LHEX) THEN
            A=0.0
            B=0.0
            DO 40 IKG=1,NGEN
            IF((A.NE.0.0).AND.(A.NE.XX(IKG)).AND.(IWIGN.EQ.2)) THEN
               CALL XABORT('SYBEUR: OPTION ASKE OR SANC REQUIRED(1).')
            ENDIF
            A=XX(IKG)
            IF((B.NE.0.0).AND.(B.NE.YY(IKG)).AND.(IWIGN.EQ.2)) THEN
               CALL XABORT('SYBEUR: OPTION ASKE OR SANC REQUIRED(2).')
            ENDIF
            B=YY(IKG)
40          CONTINUE
         ENDIF
      ELSE
*        SWITCH TO LSECT(IKG)=-999 IN CASES WHERE THE RADIUS INTERCEPT
*        THE CELL BOUNDARY.
         DO 55 IKG=1,NGEN
         I1=NMC(IKG)+1
         I2=NMC(IKG+1)
         DO 50 I=I1+1,I2
         RJ=RAYRE(I)
         IF(LHEX.AND.(LSECT(IKG).EQ.0)) THEN
            IF(2.0*RJ.GT.XX(IKG)*SQRT(3.0)) LSECT(IKG)=-999
         ELSE IF(LSECT(IKG).EQ.0) THEN
            IF(2.0*RJ.GT.MIN(XX(IKG),YY(IKG))) LSECT(IKG)=-999
         ENDIF
50       CONTINUE
55       CONTINUE
      ENDIF
*
      IF(MULTC.EQ.1) THEN
         IF(IMPX.GE.1) WRITE (6,'(/20H ROTH APPROXIMATION.)')
         IJAS=NMCEL
      ELSE IF(LHEX.AND.(MULTC.EQ.2)) THEN
         IF(IMPX.GE.1) WRITE (6,'(/24H ROTH X 6 APPROXIMATION.)')
         IJAS=6*NMCEL
      ELSE IF(MULTC.EQ.2) THEN
         IF(IMPX.GE.1) WRITE (6,'(/24H ROTH X 4 APPROXIMATION.)')
         IJAS=4*NMCEL
      ELSE IF(LHEX.AND.(MULTC.EQ.3)) THEN
         IF(IMPX.GE.1) WRITE (6,'(/27H DP-0 APPROXIMATION AROUND ,
     1   13HEACH HEXAGON.)')
         IJAS=6*NMCEL
      ELSE IF(MULTC.EQ.3) THEN
         IF(IMPX.GE.1) WRITE (6,'(/27H DP-0 APPROXIMATION AROUND ,
     1   20HEACH CARTESIAN CELL.)')
         IJAS=4*NMCEL
      ELSE IF(LHEX.AND.(MULTC.EQ.4)) THEN
         IF(IMPX.GE.1) WRITE (6,'(/27H DP-1 APPROXIMATION AROUND ,
     1   13HEACH HEXAGON.)')
         IJAS=18*NMCEL
      ELSE IF(MULTC.EQ.4) THEN
         IF(IMPX.GE.1) WRITE (6,'(/27H DP-1 APPROXIMATION AROUND ,
     1   20HEACH CARTESIAN CELL.)')
         IJAS=12*NMCEL
      ENDIF
      IF(IJAS.GT.MAXJ) THEN
         WRITE (HSMG,660) IJAS
         CALL XABORT(HSMG)
      ENDIF
      IF(IMPX.GE.1) THEN
         IF((MULTC.LE.2).AND.(IWIGN.EQ.1)) THEN
            WRITE (6,'(23H ASKEW CYLINDERIZATION./)')
         ELSE IF((MULTC.LE.2).AND.(IWIGN.EQ.2)) THEN
            WRITE (6,'(24H WIGNER CYLINDERIZATION./)')
         ELSE IF((MULTC.LE.2).AND.(IWIGN.EQ.3)) THEN
            WRITE (6,'(25H SANCHEZ CYLINDERIZATION./)')
         ENDIF
         IF(INORM.EQ.1) WRITE (6,670)
         IF(IRECT.EQ.1) WRITE (6,'(/30H SYMMETRIES OF SQUARE CELLS AR,
     1   17HE NOT CONSIDERED.)')
         IF(IQW.EQ.1) WRITE (6,'(/32H USE EQUAL WEIGHT QUADRATURES IN,
     1   17H ANGLE AND SPACE.)')
         HDOM='(0,PI/4)'
         IF(LHEX) HDOM='(0,PI/6)'
         WRITE (6,690) IQUAD(1),HDOM,IQUAD(2),IQUAD(3),IQUAD(4)
      ENDIF
      DO 100 IQQ=1,4
      DO 60 I=1,24
      IF(IQUAD(IQQ).EQ.IIT(I)) GO TO 100
60    CONTINUE
      CALL XABORT('SYBEUR: INVALID QUADRATURE PARAMETER.')
100   CONTINUE
*
      NCODE(5)=0
      NCODE(6)=0
      ALLOCATE(ZZ(NGEN))
      CALL XDRSET(ZZ,NGEN,1.0)
      IF(LHEX) THEN
         NCOUR=6
         CALL NUMERH(NCOUR,MULTC,NCODE,ZCODE,IHEX,LX,1,IORI,ISMH,POURCE,
     1   IMPX,NMCEL,IFR,ALB,SUR,NMERGE,INUM,MIX,DVX,NGEN,IGEN,XX,ZZ)
      ELSE
         NCOUR=4
         CALL NUMER3(NCOUR,MULTC,NCODE,ZCODE,LX,LY,1,IORI,ISMR,POURCE,
     1   IMPX,NMCEL,IFR,ALB,SUR,NMERGE,INUM,MIX,DVX,NGEN,IGEN,XX,YY,ZZ)
      ENDIF
      DEALLOCATE(ZZ)
*----
*  COMPUTE THE NUMBER OF DISTINCT OUT CURRENTS BETWEEN THE CELLS
*----
      IJAT=0
      DO 130 I=1,NCOUR*NMERGE
      IJAT=MAX(IJAT,MIX(I))
130   CONTINUE
      IF(IMPX.GE.1) WRITE (6,760) NMCEL,NMERGE,NGEN,MAXCEL
      IF(NMCEL.GT.MAXCEL) CALL XABORT('SYBEUR: INSUFFICIENT MAXCEL.')
      IF(IMPX.GE.1) WRITE (6,770) IJAT,IJAS,MAXJ
      IF(IJAT.GT.MAXJ) CALL XABORT('SYBEUR: INSUFFICIENT MAXJ.')
*----
*  ESTIMATION OF THE MEMORY REQUIRED TO STORE THE TRACKING INFORMATION
*----
      NA=IQUAD(1)
      NX=IQUAD(2)
      MR=IQUAD(4)
      LMAILI=0
      LMAILR=0
      DO 150 IKG=1,NGEN
      J2=NMC(IKG+1)-NMC(IKG)
      J1R=NMCR(IKG)
      J2R=NMCR(IKG+1)-J1R
      IF((J2.EQ.1).AND.(LSECT(IKG).NE.-999)) THEN
         GO TO 150
      ELSE IF(LHEX.AND.LSECT(IKG).NE.0) THEN
*        TRACKING FOR AN HEXAGONAL SECTORIZED CELL.
         IF(LSECT(IKG).EQ.-999) THEN
            NSECT=6
         ELSE
            NSECT=6*MOD(ABS(LSECT(IKG)),100)
         ENDIF
         DELR=XX(IKG)/REAL(NX*J2R)
         ALLOCATE(VOL2(J2))
         CALL SYB7TS(NA,J2R,NSECT,LSECT(IKG),J2,XX(IKG),RAYRE(J1R+2),
     1   ILIGN,IQW,DELR,.FALSE.,VOL2,LR,ZZR,LI,ZZI)
         DEALLOCATE(VOL2)
         LMAILI=LMAILI+LI
         LMAILR=LMAILR+LR
      ELSE IF(LSECT(IKG).NE.0) THEN
*        TRACKING FOR A CARTESIAN SECTORIZED CELL.
         IF(LSECT(IKG).EQ.-999) THEN
            NSECT=4
         ELSE IF((LSECT(IKG).EQ.-1).OR.(LSECT(IKG).EQ.-101)) THEN
            NSECT=8
         ELSE
            NSECT=4*MOD(ABS(LSECT(IKG)),100)
         ENDIF
         DELR=SQRT(XX(IKG)**2+YY(IKG)**2)/REAL(2*NX*J2R)
         ALLOCATE(VOL2(J2))
         CALL SYB4TS(NA,J2R,NSECT,LSECT(IKG),J2,XX(IKG),YY(IKG),
     1   RAYRE(J1R+2),ILIGN,IQW,DELR,.FALSE.,VOL2,LR,ZZR,LI,ZZI)
         DEALLOCATE(VOL2)
         LMAILI=LMAILI+LI
         LMAILR=LMAILR+LR
      ELSE IF(MULTC.LE.2) THEN
         NPIJ=J2R
         LMAILR=LMAILR+1+MR*((NPIJ*(5+NPIJ))/2)
      ELSE IF((MULTC.GE.3).AND.(J2R.GT.1)) THEN
         NPIJ=J2R-1
         LMAILR=LMAILR+1+MR*((NPIJ*(5+NPIJ))/2)
         IF(LHEX) THEN
            LMAILI=LMAILI+4+3*NA*(2+(J2R+1)*(3+2*J2R))
            LMAILR=LMAILR+J2R+4+3*NA*(13+2*(J2R+1)*NX*J2R)
         ELSE IF((XX(IKG).EQ.YY(IKG)).AND.(IRECT.NE.1)) THEN
            LMAILI=LMAILI+4+2*NA*(2+(J2R+1)*(3+2*J2R))
            LMAILR=LMAILR+J2R+4+2*NA*(9+2*(J2R+1)*NX*J2R)
         ELSE
            LMAILI=LMAILI+4+2*NA*(2+(2*J2R+1)*(3+2*J2R))
            LMAILR=LMAILR+J2R+4+2*NA*(9+2*(2*J2R+1)*NX*J2R)
         ENDIF
      ENDIF
150   CONTINUE
      IF((IMPX.GE.1).OR.(LMAILR.GT.MAXZ).OR.(LMAILI.GT.MAXZ)) THEN
         WRITE (6,680) LMAILI,LMAILR,MAXZ
      ENDIF
      IF((LMAILR.GT.MAXZ).OR.(LMAILR.GT.MAXZ)) WRITE (6,'(/9H SYBEUR: ,
     1 60H*** WARNING *** MAXZ IS LESS THAN THE ESTIMATED TRACKING LEN,
     2 4HGTH.)')
      CALL XDISET(IZMAIL,LMAILI,0)
      CALL XDRSET(RZMAIL,LMAILR,0.0)
*----
*  COMPUTE THE VOLUMES AND THE TRACKING INFORMATION
*----
      IF(IHALT.EQ.1) CALL XABORT('SYBEUR: STOP REQUESTED.')
      LMAILI=0
      LMAILR=0
      PRECC=0.0
      ALLOCATE(VOL2(NMC(NGEN+1)))
      DO 210 IKG=1,NGEN
      MAIL(1,IKG)=0
      MAIL(2,IKG)=0
      J1=NMC(IKG)
      J2=NMC(IKG+1)-J1
      J1R=NMCR(IKG)
      J2R=NMCR(IKG+1)-J1R
      IF(MULTC.LE.2) THEN
         ALLOCATE(RAD(J2R+1))
         DO 180 I=1,J2R
         RAD(I)=RAYRE(J1R+I)
180      CONTINUE
         IF(LHEX.AND.(IWIGN.EQ.1)) THEN
*           ASKEW CYLINDERIZATION HEXAGONAL.
            RAD(J2R+1)=3.0*XX(IKG)/PI
         ELSE IF(LHEX.AND.(IWIGN.GE.2)) THEN
*           WIGNER OR SANCHEZ CYLINDERIZATION HEXAGONAL.
            RAD(J2R+1)=SQRT(1.5*SQRT(3.0)/PI)*XX(IKG)
         ELSE IF(IWIGN.EQ.1) THEN
*           ASKEW CYLINDERIZATION CARTESIAN.
            RAD(J2R+1)=(XX(IKG)+YY(IKG))/PI
         ELSE IF(IWIGN.GE.2) THEN
*           WIGNER OR SANCHEZ CYLINDERIZATION CARTESIAN.
            RAD(J2R+1)=SQRT(XX(IKG)*YY(IKG)/PI)
         ENDIF
         CALL SYBT1D(J2R,RAD,.FALSE.,MR,RZMAIL(LMAILR+1))
         DEALLOCATE(RAD)
         LI=0
         LR=1+MR*((J2R*(5+J2R))/2)
      ELSE
         IF((J2.EQ.1).AND.(LSECT(IKG).NE.-999)) THEN
            GO TO 190
         ELSE IF((LSECT(IKG).NE.0).AND.LHEX) THEN
*           TRACKING FOR AN HEXAGONAL SECTORIZED CELL.
            IF(LSECT(IKG).EQ.-999) THEN
               NSECT=6
            ELSE
               NSECT=6*MOD(ABS(LSECT(IKG)),100)
            ENDIF
            DELR=XX(IKG)/REAL(NX*J2R)
            CALL SYB7TS(NA,J2R,NSECT,LSECT(IKG),J2,XX(IKG),RAYRE(J1R+2),
     1      ILIGN,IQW,DELR,.TRUE.,VOL2(J1+1),LR,RZMAIL(LMAILR+1),
     2      LI,IZMAIL(LMAILI+1))
         ELSE IF(LSECT(IKG).NE.0) THEN
*           TRACKING FOR A CARTESIAN SECTORIZED CELL.
            IF(LSECT(IKG).EQ.-999) THEN
               NSECT=4
            ELSE IF((LSECT(IKG).EQ.-1).OR.(LSECT(IKG).EQ.-101)) THEN
               NSECT=8
            ELSE
               NSECT=4*MOD(ABS(LSECT(IKG)),100)
            ENDIF
            DELR=SQRT(XX(IKG)**2+YY(IKG)**2)/REAL(2*NX*J2R)
            CALL SYB4TS(NA,J2R,NSECT,LSECT(IKG),J2,XX(IKG),YY(IKG),
     1      RAYRE(J1R+2),ILIGN,IQW,DELR,.TRUE.,VOL2(J1+1),LR,
     2      RZMAIL(LMAILR+1),LI,IZMAIL(LMAILI+1))
         ELSE IF(LHEX) THEN
            CALL SYBHTK(NA,NX,J2R,XX(IKG),RAYRE(J1R+1),ILIGN,INORM,
     1      IQW,LR,RZMAIL(LMAILR+1),LI,IZMAIL(LMAILI+1),PREC)
         ELSE
            CALL SYBRTK(NA,NX,J2R,XX(IKG),YY(IKG),RAYRE(J1R+1),ILIGN,
     1      INORM,IRECT,IQW,LR,RZMAIL(LMAILR+1),LI,IZMAIL(LMAILI+1),
     2      PREC)
         ENDIF
         PRECC=AMAX1(PRECC,PREC)
      ENDIF
      MAIL(1,IKG)=LMAILI+1
      MAIL(2,IKG)=LMAILR+1
      LMAILI=LMAILI+LI
      LMAILR=LMAILR+LR
      IF(LMAILI.GT.MAXZ) CALL XABORT('SYBEUR: INSUFFICIENT MAXZ(1).')
      IF(LMAILR.GT.MAXZ) CALL XABORT('SYBEUR: INSUFFICIENT MAXZ(2).')
*
190   IF(LSECT(IKG).EQ.0) THEN
*        COMPUTE THE VOLUMES IN A NON-SECTORIZED GENERATING CELL.
         RJ=0.0
         DO 200 I=1,J2R-1
         RJ1=RAYRE(J1R+I+1)**2
         VOL2(J1+I)=PI*(RJ1-RJ)
         RJ=RJ1
200      CONTINUE
         IF(LHEX) THEN
            VOL2(J1+J2)=1.5*SQRT(3.0)*XX(IKG)**2-PI*RJ
         ELSE
            VOL2(J1+J2)=XX(IKG)*YY(IKG)-PI*RJ
         ENDIF
      ENDIF
210   CONTINUE
      IF(IMPX.GE.1) WRITE (6,710) LMAILI,LMAILR,PRECC
*----
*  WEIGHTING VOLUMES BY POURCE
*----
      I1=0
      DO 230 IKK=1,NMERGE
      IKG=IGEN(IKK)
      J1=NMC(IKG)
      I2=NMC(IKG+1)-J1
      IF(I1+I2.GT.MAXPTS) THEN
         WRITE(HSMG,790) I1+I2
         CALL XABORT(HSMG)
      ENDIF
      DO 220 I=1,I2
      I1=I1+1
      IF(VOL2(J1+I).LE.0.0) THEN
         WRITE(HSMG,'(41HSYBEUR: NEGATIVE OR ZERO VOLUME IN REGION,
     1   I5,1H.)') J1+I
         CALL XABORT(HSMG)
      ENDIF
      VOL(I1)=VOL2(J1+I)*POURCE(IKK)
220   CONTINUE
230   CONTINUE
      DEALLOCATE(VOL2)
      IF(I1.NE.NREG) CALL XABORT('SYBEUR: WRONG NUMBER OF REGIONS.')
      IF(IMPX.GE.1) WRITE (6,780) NREG,MAXPTS
      IF(NREG.GT.MAXPTS) CALL XABORT('SYBEUR: INSUFFICIENT MAXPTS.')
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(POURCE,IORI)
      RETURN
*
660   FORMAT (42HSYBEUR: MAXJ SHOULD BE GREATER OR EQUAL TO,I7,1H.)
670   FORMAT (52H THE TRACKS ARE NOT NORMALIZED AND THE GEOMETRIC ACC,
     1 43HURACY OF THE TRACKING MESH IS NOT COMPUTED./)
680   FORMAT (/49H SYBEUR: MAXIMUM STORAGE REQUIRED BY THE TRACKS =,I8,
     1 2H +,I7,8X,24HAVAILABLE STORAGE MAXZ =,I8/)
690   FORMAT (/23H QUADRATURE PARAMETERS:/24H HETEROGENEOUS CELLS   -,
     1 17H    ANGLES   NA =,I3,4H IN ,A8/24X,17H  SEGMENTS   NX =,
     2 I3//41H HOMOGENEOUS CELLS     -    ANGLES   NH =,I3//
     3 41H CYLINDRICAL CELLS     -  SEGMENTS   MR =,I3//)
710   FORMAT (/47H SYBEUR: EXACT STORAGE REQUIRED BY THE TRACKS =,I7,
     1 6H INT +,I7,5H REAL,8X,10HACCURACY =,1P,E11.3/)
760   FORMAT (/38H TOTAL NUMBER OF CELLS IN THE DOMAIN =,I7/
     1         31H TOTAL NUMBER OF MERGED CELLS =,I7/
     2         35H TOTAL NUMBER OF GENERATING CELLS =,I7,5X,
     3         26HAVAILABLE STORAGE MAXCEL =,I7/)
770   FORMAT (/40H TOTAL NUMBER OF DISTINCT OUT-CURRENTS =,I7/
     1 57H TOTAL NUMBER OF OUT-CURRENTS SURROUNDING ALL THE CELLS =,I7,
     2 5X,24HAVAILABLE STORAGE MAXJ =,I7/)
780   FORMAT (/20H NUMBER OF VOLUMES =,I7,5X,22HAVAILABLE STORAGE MAXB,
     1 4HLK =,I7/)
790   FORMAT (44HSYBEUR: MAXPTS SHOULD BE GREATER OR EQUAL TO,I7,1H.)
      END
