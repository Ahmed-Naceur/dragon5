*DECK SYBTRK
      SUBROUTINE SYBTRK (IPTRK,IPGEOM,IMPX,MAXPTS,MAXJ,MAXZ,MULTC,IWIGN,
     1 IHALT,ILIGN,INORM,IRECT,IQW,JQUA1,JQUA2,IQUA10,IBIHET,FRTM)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover of the geometry and tracking for Sybil modules.
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
* IPTRK   pointer to the Sybil tracking (L_TRACK signature).
* IPGEOM  pointer to the geometry (L_GEOM signature).
* IMPX    print flag.
* MAXPTS  allocated storage for arrays of dimension NREG.
* MAXJ    allocated storage for interface current arrays.
* MAXZ    allocated storage for tracking arrays.
* MULTC   type of multicell approximation in Eurydice.
* IWIGN   type of cylinderization.
* IHALT   stop flag at the end of tracking (set with IHALT=1).
* ILIGN   on/off switch for track printout.
* INORM   on/off switch for track normalization.
* IRECT   on/off switch for using symmetries in square cells.
* IQW     type of quadrature.
* JQUA1   1-D quadrature parameter.
* JQUA2   2-D quadrature parameters.
* IQUA10  quadrature parameter for micro-structures in Bihet.
* IBIHET  type of double-heterogeneity method (=1 Sanchez-Pomraning
*         model; =2 Hebert model; =3 She-Liu-Shi model (no shadow);
*         =4 She-Liu-Shi model (with shadow)).
* FRTM    minimum volume fraction of the grain in the representative 
*         volume for She-Liu-Shi model.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPGEOM
      INTEGER IMPX,MAXPTS,MAXJ,MAXZ,MULTC,IWIGN,IHALT,ILIGN,INORM,
     1 IRECT,IQW,JQUA1,JQUA2(2),IQUA10,IBIHET
      REAL FRTM
*----
*  LOCAL VARIABLES
*----
      PARAMETER (PI=3.141592654,NSTATE=40)
      LOGICAL ILK,LBIHET
      INTEGER ISTATE(NSTATE),IQUAD(4),IGP(NSTATE),IPARAM(16)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MAT,IDL,NCODE,ICODE,ISPLX,
     1 ISPLY,ISPLZ,NMC3,LSECT,NMC4,NMCR4,MAIL,IZMAI,IFR,INUM,MIX,IGEN
      REAL, ALLOCATABLE, DIMENSION(:) :: VOL,XX2,YY2,ZZ2,ZCODE,RAYR3,
     1 PROCE,POURC,SURFA,XX4,YY4,RAYR4,RZMAI,ALB,SUR,DVX
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(MAT(MAXPTS),IDL(MAXPTS))
      ALLOCATE(VOL(MAXPTS))
*
      ITG=0
      CALL LCMGET(IPGEOM,'STATE-VECTOR',ISTATE)
      LBIHET=(ISTATE(12).EQ.1)
      IF(ISTATE(1).EQ.1) THEN
         ITG=1
         IF(ISTATE(6).NE.1) CALL XABORT('SYBTRK: INVALID NUMBER OF REGI'
     1   //'ONS.')
         CALL LCMLEN(IPGEOM,'MIX',NMILG,ITYLCM)
         IF(NMILG.NE.1) CALL XABORT('SYBTRK: INVALID MIX VECTOR.')
         CALL LCMGET(IPGEOM,'MIX',MAT(1))
         IR=MAT(1)
         VOL(1)=1.0
         ILK=.FALSE.
         NREG=1
      ELSE IF((ISTATE(1).GT.1).AND.(ISTATE(1).LT.5).AND.(ISTATE(9).EQ.
     1 0)) THEN
         ITG=2
         ALLOCATE(NCODE(6),ICODE(6))
         ALLOCATE(XX2(MAXPTS+1),YY2(MAXPTS+1),ZZ2(MAXPTS+1),ZCODE(6))
*
         ALLOCATE(ISPLX(MAXPTS),ISPLY(MAXPTS),ISPLZ(MAXPTS))
         CALL READ3D(MAXPTS,MAXPTS,MAXPTS,MAXPTS,IPGEOM,IHEX,IR,ILK,
     1   SIDE,XX2,YY2,ZZ2,IMPX,LX,LY,LZ,MAT,NREG,NCODE,ICODE,ZCODE,
     2   ISPLX,ISPLY,ISPLZ,ISPLH,ISPLL)
         DEALLOCATE(ISPLZ,ISPLY,ISPLX)
         DO 10 IC=1,6
         IF(NCODE(IC).EQ.7) CALL XABORT('SYBTRK: ZERO FLUX BOUNDARY CO'
     1   //'NDITION NOT PERMITTED.')
   10    CONTINUE
*
*        COMPUTATION OF THE VOLUMES.
         DO 20 IKK=1,NREG
         A=XX2(IKK)
         B=XX2(IKK+1)
         IF(ISTATE(1).EQ.2) THEN
            VOL(IKK)=(B-A)
         ELSE IF(ISTATE(1).EQ.3) THEN
            VOL(IKK)=PI*(B-A)*(B+A)
         ELSE IF(ISTATE(1).EQ.4) THEN
            VOL(IKK)=4.0*PI*(B-A)*(A*A+A*B+B*B)/3.0
         ENDIF
   20    CONTINUE
         IF(IMPX.GE.1) WRITE (6,'(/29H QUADRATURE PARAMETER JQUA1 =,
     1   I2/)') JQUA1
*
         IPARAM(1)=ISTATE(1)
         IPARAM(2)=IHEX
         IPARAM(3)=JQUA1
         IPARAM(4)=LX
         IPARAM(5)=LY
         IPARAM(6)=LZ
         CALL LCMSIX(IPTRK,'PURE-GEOM',1)
         CALL LCMPUT(IPTRK,'PARAM',6,1,IPARAM)
         CALL LCMPUT(IPTRK,'XXX',LX+1,2,XX2)
         CALL LCMPUT(IPTRK,'YYY',LY+1,2,YY2)
         CALL LCMPUT(IPTRK,'ZZZ',LZ+1,2,ZZ2)
         CALL LCMPUT(IPTRK,'NCODE',6,1,NCODE)
         CALL LCMPUT(IPTRK,'ICODE',6,1,ICODE)
         CALL LCMPUT(IPTRK,'ZCODE',6,2,ZCODE)
         DEALLOCATE(ZCODE,ICODE,NCODE)
         IF(ISTATE(1).GE.8) CALL LCMPUT(IPTRK,'SIDE',1,2,SIDE)
         CALL LCMSIX(IPTRK,' ',2)
         DEALLOCATE(ZZ2,YY2,XX2)
      ELSE IF(ISTATE(1).EQ.30) THEN
         ITG=3
         ALLOCATE(NMC3(1+MAXPTS))
         ALLOCATE(RAYR3(2*MAXPTS),PROCE(MAXPTS**2),POURC(MAXPTS),
     1   SURFA(MAXPTS))
*
         CALL READMT (MAXPTS,IPGEOM,IR,MAT,VOL,ILK,ISTAT,NSUPCE,
     1   NREG,NMC3,RAYR3,PROCE,POURC,
     2   SURFA,IMPX)
         IF(IMPX.GE.1) WRITE (6,'(/29H QUADRATURE PARAMETER JQUA1 =,
     1   I2/)') JQUA1
*
         IPARAM(1)=NSUPCE
         IPARAM(2)=JQUA1
         IPARAM(3)=ISTAT
         CALL LCMSIX(IPTRK,'DOITYOURSELF',1)
         CALL LCMPUT(IPTRK,'PARAM',3,1,IPARAM)
         CALL LCMPUT(IPTRK,'NMC',1+NSUPCE,1,NMC3)
         CALL LCMPUT(IPTRK,'RAYRE',NREG+NSUPCE,2,RAYR3)
         CALL LCMPUT(IPTRK,'PROCEL',NSUPCE**2,2,PROCE)
         CALL LCMPUT(IPTRK,'POURCE',NSUPCE,2,POURC)
         CALL LCMPUT(IPTRK,'SURFA',NSUPCE,2,SURFA)
         CALL LCMSIX(IPTRK,' ',2)
         DEALLOCATE(SURFA,POURC,PROCE,RAYR3,NMC3)
      ELSE IF( (ISTATE(1).EQ.5).OR.(ISTATE(1).EQ.8) .OR.
     1         ((ISTATE(1).EQ.20).AND.(ISTATE(13).EQ.0)) .OR.
     2         ((ISTATE(1).EQ.24).AND.(ISTATE(13).EQ.0)) ) THEN
         ITG=4
         MAXCEL=MAXPTS
         ALLOCATE(LSECT(MAXCEL),NMC4(MAXCEL+1),NMCR4(MAXCEL+1),
     1   MAIL(2*MAXCEL),IZMAI(MAXZ),IFR(MAXJ),INUM(MAXCEL),MIX(MAXJ),
     2   IGEN(MAXCEL))
         ALLOCATE(XX4(MAXCEL),YY4(MAXCEL),RAYR4(MAXPTS),RZMAI(MAXZ),
     1   ALB(MAXJ),SUR(MAXJ),DVX(MAXJ))
         IQUAD(1)=JQUA2(1)
         IQUAD(2)=JQUA2(2)
         IQUAD(3)=JQUA1
         IQUAD(4)=JQUA1
*
         CALL SYBEUR(MAXPTS,MAXCEL,MAXJ,MAXZ,IPGEOM,NREG,IR,MAT,VOL,
     1   ILK,IMPX,IHEX,NCOUR,LMAILI,LMAILR,NMCEL,NMERGE,NGEN,IJAT,MULTC,
     2   IWIGN,IHALT,ILIGN,INORM,IRECT,IQW,IQUAD,XX4,YY4,LSECT,NMC4,
     3   NMCR4,RAYR4,MAIL,IZMAI,RZMAI,IFR,ALB,SUR,INUM,MIX,DVX,IGEN)
*
         IPARAM(1)=IHEX
         IPARAM(2)=MULTC
         IPARAM(3)=IWIGN
         IPARAM(4)=NMCEL
         IPARAM(5)=NMERGE
         IPARAM(6)=NGEN
         IPARAM(7)=IJAT
         IPARAM(8)=IQUAD(1)
         IPARAM(9)=IQUAD(2)
         IPARAM(10)=IQUAD(3)
         IPARAM(11)=IQUAD(4)
         IPARAM(12)=INORM
         IPARAM(13)=IQW
         IPARAM(14)=NCOUR
         IPARAM(15)=LMAILI
         IPARAM(16)=LMAILR
         IRDIM=NMCR4(NGEN+1)
         CALL LCMSIX(IPTRK,'EURYDICE',1)
         CALL LCMPUT(IPTRK,'PARAM',16,1,IPARAM)
         CALL LCMPUT(IPTRK,'XX',NGEN,2,XX4)
         CALL LCMPUT(IPTRK,'YY',NGEN,2,YY4)
         CALL LCMPUT(IPTRK,'LSECT',NGEN,1,LSECT)
         CALL LCMPUT(IPTRK,'NMC',1+NGEN,1,NMC4)
         CALL LCMPUT(IPTRK,'NMCR',1+NGEN,1,NMCR4)
         CALL LCMPUT(IPTRK,'RAYRE',IRDIM,2,RAYR4)
         CALL LCMPUT(IPTRK,'MAIL',2*NGEN,1,MAIL)
         IF(LMAILI.GT.0) THEN
            CALL LCMPUT(IPTRK,'ZMAILI',LMAILI,1,IZMAI)
         ENDIF
         IF(LMAILR.GT.0) THEN
            CALL LCMPUT(IPTRK,'ZMAILR',LMAILR,2,RZMAI)
         ENDIF
         CALL LCMPUT(IPTRK,'IFR',NCOUR*NMCEL,1,IFR)
         CALL LCMPUT(IPTRK,'ALB',NCOUR*NMCEL,2,ALB)
         CALL LCMPUT(IPTRK,'SUR',NCOUR*NMCEL,2,SUR)
         CALL LCMPUT(IPTRK,'INUM',NMCEL,1,INUM)
         CALL LCMPUT(IPTRK,'MIX',NCOUR*NMERGE,1,MIX)
         CALL LCMPUT(IPTRK,'DVX',NCOUR*NMERGE,2,DVX)
         CALL LCMPUT(IPTRK,'IGEN',NMERGE,1,IGEN)
         CALL LCMSIX(IPTRK,' ',2)
         DEALLOCATE(DVX,SUR,ALB,RZMAI,RAYR4,YY4,XX4)
         DEALLOCATE(IGEN,MIX,INUM,IFR,IZMAI,MAIL,NMCR4,NMC4,LSECT)
      ELSE
         CALL XABORT('SYBTRK: INVALID GEOMETRY MODULE.')
      ENDIF
*----
*  SAVE GENERAL AND SYBIL-SPECIFIC TRACKING INFORMATION
*----
      DO 30 I=1,NREG
      IDL(I)=I
   30 CONTINUE
      IF(ITG.EQ.3) THEN
         NUNCUR=NSUPCE
      ELSE IF(ITG.EQ.4) THEN
         NUNCUR=IJAT
      ELSE
         NUNCUR=0
      ENDIF
      CALL XDISET(IGP,NSTATE,0)
      IGP(1)=NREG
      IGP(2)=NREG
      IF(ILK) THEN
         IGP(3)=0
      ELSE
         IGP(3)=1
      ENDIF
      IGP(4)=IR
      IGP(5)=0
      IGP(6)=ITG
      IGP(7)=MAXZ
      IGP(8)=MAXJ
      IGP(9)=NUNCUR
      IF(LBIHET) IGP(40)=1
      CALL LCMPUT(IPTRK,'STATE-VECTOR',NSTATE,1,IGP)
      CALL LCMPUT(IPTRK,'MATCOD',NREG,1,MAT)
      CALL LCMPUT(IPTRK,'VOLUME',NREG,2,VOL)
      CALL LCMPUT(IPTRK,'KEYFLX',NREG,1,IDL)
*----
*  DOUBLE HETEROGENEITY OPTION
*----
      IF(LBIHET) CALL XDRTBH(IPGEOM,IPTRK,IQUA10,IBIHET,IMPX,FRTM)
*----
*  PRINT TRACKING ARRAYS
*----
      IF(IMPX.GT.5) THEN
         CALL LCMGET(IPTRK,'STATE-VECTOR',IGP)
         NREG=IGP(1)
         CALL LCMGET(IPTRK,'MATCOD',MAT)
         CALL LCMGET(IPTRK,'VOLUME',VOL)
         CALL LCMGET(IPTRK,'KEYFLX',IDL)
         I1=1
         DO 60 I=1,(NREG-1)/8+1
         I2=I1+7
         IF(I2.GT.NREG) I2=NREG
         WRITE (6,620) (J,J=I1,I2)
         WRITE (6,630) (MAT(J),J=I1,I2)
         WRITE (6,640) (IDL(J),J=I1,I2)
         WRITE (6,650) (VOL(J),J=I1,I2)
         I1=I1+8
   60    CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(VOL)
      DEALLOCATE(IDL,MAT)
      RETURN
*
  620 FORMAT (///11H REGION    ,8(I8,6X,1HI))
  630 FORMAT (   11H MIXTURE   ,8(I8,6X,1HI))
  640 FORMAT (   11H POINTER   ,8(I8,6X,1HI))
  650 FORMAT (   11H VOLUME    ,8(1P,E13.6,2H I))
      END
