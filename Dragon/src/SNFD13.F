*DECK SNFD13
      SUBROUTINE SNFD13(IOMP,NUN,NGEFF,IMPX,INCONV,NGIND,LX,LY,LZ,
     1 IELEM,NMAT,NPQ,NSCT,MAT,VOL,TOTAL,NCODE,ZCODE,QEXT,LFIXUP,DU,
     2 DE,DZ,W,MRMX,MRMY,MRMZ,DC,DB,DA,PL,FUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform one inner iteration for solving SN equations in 3D Cartesian
* geometry for the HODD method. Domino-type multithreading. Albedo
* boundary conditions.
*
*Copyright:
* Copyright (C) 2020 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): N. Martin and A. Hebert
*
*Parameters: input
* IOMP    number of cells per axis in a Domino macrocell.
* NUN     total number of unknowns in vector FUNKNO.
* NGEFF   number of energy groups processed in parallel.
* IMPX    print flag (equal to zero for no print).
* INCONV  energy group convergence flag (set to .FALSE. if converged).
* NGIND   energy group indices assign to the NGEFF set.
* LX      number of meshes along X axis.
* LY      number of meshes along Y axis.
* LZ      number of meshes along Z axis.
* IELEM   degree of the spatial approximation (=1 diamond scheme).
* NMAT    number of material mixtures.
* NPQ     number of SN directions in height octants.
* NSCT    maximum number of spherical harmonics moments of the flux.
* MAT     material mixture index in each region.
* VOL     volumes of each region.
* TOTAL   macroscopic total cross sections.
* NCODE   boundary condition indices.
* ZCODE   albedos.
* QEXT    Legendre components of the fixed source.
* LFIXUP  flag to enable negative flux fixup.
* DU      first direction cosines ($\\mu$).
* DE      second direction cosines ($\\eta$).
* DZ      third direction cosines ($\\xi$).
* W       weights.
* MRMX    quadrature index.
* MRMY    quadrature index.
* MRMZ    quadrature index.
* DC      diamond-scheme parameter.
* DB      diamond-scheme parameter.
* DA      diamond-scheme parameter.
* PL      discrete values of the spherical harmonics corresponding
*         to the 3D SN quadrature.
*
*Parameters: input/output
* FUNKNO  Legendre components of the flux and boundary fluxes.
*
*-----------------------------------------------------------------------
*
#if defined(_OPENMP)
      USE omp_lib
#endif
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IOMP,NUN,NGEFF,IMPX,NGIND(NGEFF),LX,LY,LZ,IELEM,NMAT,
     1 NPQ,NSCT,MAT(LX,LY,LZ),NCODE(6),MRMX(NPQ),MRMY(NPQ),MRMZ(NPQ)
      LOGICAL INCONV(NGEFF)
      REAL VOL(LX,LY,LZ),TOTAL(0:NMAT,NGEFF),ZCODE(6),QEXT(NUN,NGEFF),
     1 DU(NPQ),DE(NPQ),DZ(NPQ),W(NPQ),DC(LX,LY,NPQ),DB(LX,LZ,NPQ),
     2 DA(LY,LZ,NPQ),PL(NSCT,NPQ),FUNKNO(NUN,NGEFF)
      LOGICAL LFIXUP
*----
*  LOCAL VARIABLES
*----
      INTEGER NPQD(8),IIND(8)
      PARAMETER(IUNOUT=6,RLOG=1.0E-8,PI=3.141592654)
      DOUBLE PRECISION VT,CONST0,CONST1,CONST2
      DOUBLE PRECISION Q(IELEM**3),Q2(IELEM**3,IELEM**3+1)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: III,JJJ
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: INDANG
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: FLUX
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: FLUX_G,
     1 XNI,XNJ,XNK
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(INDANG(NPQ,8))
      ALLOCATE(XNI(IELEM,IELEM,LY,LZ,NPQ,NGEFF),
     1 XNJ(IELEM,IELEM,LX,LZ,NPQ,NGEFF),
     2 XNK(IELEM,IELEM,LX,LY,NPQ,NGEFF),
     3 FLUX(IELEM**3,NSCT),FLUX_G(IELEM**3,NSCT,LX,LY,LZ,NGEFF))
*----
*  DEFINITION OF CONSTANTS.
*----
      NMX=1+(LX-1)/IOMP
      NMY=1+(LY-1)/IOMP
      NMZ=1+(LZ-1)/IOMP
      L4=IELEM*IELEM*IELEM*LX*LY*LZ*NSCT
      IOFY=IELEM*IELEM*LY*LZ*NPQ
      IOFZ=IELEM*IELEM*(LX+LY)*LZ*NPQ
      CONST0=2.0D0*DSQRT(3.0D0)
      CONST1=2.0D0*DSQRT(5.0D0)
      CONST2=2.0D0*DSQRT(15.0D0)
*----
*  PARAMETER VALIDATION.
*----
      IF(IELEM.GT.3) CALL XABORT('SNFD13: INVALID IELEM (DIAM) VALUE. '
     1 //'CHECK INPUT DATA FILE.')
*----
*  SET OCTANT SWAPPING ORDER.
*----
      NPQD(:8)=0
      INDANG(:NPQ,:8)=0
      DO 10 M=1,NPQ
        VU=DU(M)
        VE=DE(M)
        VZ=DZ(M)
        IF((VU.GE.0.0).AND.(VE.GE.0.0).AND.(VZ.GE.0.0)) THEN
          IND=1
          JND=8
        ELSE IF((VU.LE.0.0).AND.(VE.GE.0.0).AND.(VZ.GE.0.0)) THEN
          IND=2
          JND=7
        ELSE IF((VU.LE.0.0).AND.(VE.LE.0.0).AND.(VZ.GE.0.0)) THEN
          IND=3
          JND=5
        ELSE IF((VU.GE.0.0).AND.(VE.LE.0.0).AND.(VZ.GE.0.0)) THEN
          IND=4
          JND=6
        ELSE IF((VU.GE.0.0).AND.(VE.GE.0.0).AND.(VZ.LE.0.0)) THEN
          IND=5
          JND=4
        ELSE IF((VU.LE.0.0).AND.(VE.GE.0.0).AND.(VZ.LE.0.0)) THEN
          IND=6
          JND=3
        ELSE IF((VU.LE.0.0).AND.(VE.LE.0.0).AND.(VZ.LE.0.0)) THEN
          IND=7
          JND=1
        ELSE
          IND=8
          JND=2
        ENDIF
        IIND(JND)=IND
        NPQD(IND)=NPQD(IND)+1
        INDANG(NPQD(IND),IND)=M
   10 CONTINUE
*----
*  MAIN LOOP OVER OCTANTS.
*----
      FLUX_G(:IELEM**3,:NSCT,:LX,:LY,:LZ,:NGEFF)=0.0D0
      DO 520 JND=1,8
      IND=IIND(JND)
*----
*  PRELIMINARY LOOPS FOR SETTING BOUNDARY CONDITIONS.
*----
*$OMP  PARALLEL DO
*$OMP1 PRIVATE(M,IG,VU,VE,VZ,M1,E1,IOF,JOF,IEL,JEL,I,J,K)
*$OMP2 SHARED(FUNKNO) COLLAPSE(2)
      DO 150 IG=1,NGEFF
      DO 140 IPQD=1,NPQD(IND)
      IF(.NOT.INCONV(IG)) GO TO 140
      M=INDANG(IPQD,IND)
      VU=DU(M)
      VE=DE(M)
      VZ=DZ(M)
*
      IF(VU.GT.0.0)THEN
         M1=MRMX(M)
         IF(NCODE(1).NE.4)THEN
         DO JEL=1,IELEM
         DO IEL=1,IELEM
            DO J=1,LY
            DO K=1,LZ
             IOF=((((M-1)*LZ+(K-1))*LY+(J-1))*IELEM+(JEL-1))*IELEM+IEL
             JOF=((((M1-1)*LZ+(K-1))*LY+(J-1))*IELEM+(JEL-1))*IELEM+IEL
             FUNKNO(L4+IOF,IG)=FUNKNO(L4+JOF,IG)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDIF
      ELSEIF(VU.LT.0.0)THEN
         M1=MRMX(M)
         IF(NCODE(2).NE.4)THEN
         DO JEL=1,IELEM
         DO IEL=1,IELEM
            DO J=1,LY
            DO K=1,LZ
             IOF=((((M-1)*LZ+(K-1))*LY+(J-1))*IELEM+(JEL-1))*IELEM+IEL
             JOF=((((M1-1)*LZ+(K-1))*LY+(J-1))*IELEM+(JEL-1))*IELEM+IEL
             FUNKNO(L4+IOF,IG)=FUNKNO(L4+JOF,IG)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDIF
      ENDIF
      IF(VE.GT.0.0)THEN
         M1=MRMY(M)
         IF(NCODE(3).NE.4)THEN
         DO IEL=1,IELEM
         DO JEL=1,IELEM
            DO I=1,LX
            DO K=1,LZ
             IOF=((((M-1)*LZ+(K-1))*LX+(I-1))*IELEM+(JEL-1))*IELEM+IEL
             JOF=((((M1-1)*LZ+(K-1))*LX+(I-1))*IELEM+(JEL-1))*IELEM+IEL
             FUNKNO(L4+IOFY+IOF,IG)=FUNKNO(L4+IOFY+JOF,IG)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDIF
      ELSEIF(VE.LT.0.0)THEN
         M1=MRMY(M)
         IF(NCODE(4).NE.4)THEN
         DO IEL=1,IELEM
         DO JEL=1,IELEM
            DO I=1,LX
            DO K=1,LZ
             IOF=((((M-1)*LZ+(K-1))*LX+(I-1))*IELEM+(JEL-1))*IELEM+IEL
             JOF=((((M1-1)*LZ+(K-1))*LX+(I-1))*IELEM+(JEL-1))*IELEM+IEL
             FUNKNO(L4+IOFY+IOF,IG)=FUNKNO(L4+IOFY+JOF,IG)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDIF
      ENDIF
      IF(VZ.GT.0.0)THEN
         M1=MRMZ(M)
         IF(NCODE(5).NE.4)THEN
         DO IEL=1,IELEM
         DO JEL=1,IELEM
            DO I=1,LX
            DO J=1,LY
             IOF=((((M-1)*LY+(J-1))*LX+(I-1))*IELEM+(JEL-1))*IELEM+IEL
             JOF=((((M1-1)*LY+(J-1))*LX+(I-1))*IELEM+(JEL-1))*IELEM+IEL
             E1=FUNKNO(L4+IOFZ+IOF,IG)
             FUNKNO(L4+IOFZ+IOF,IG)=FUNKNO(L4+IOFZ+JOF,IG)
             FUNKNO(L4+IOFZ+JOF,IG)=E1
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDIF
      ELSEIF(VZ.LT.0.0)THEN
         M1=MRMZ(M)
         IF(NCODE(6).NE.4)THEN
         DO IEL=1,IELEM
         DO JEL=1,IELEM
            DO I=1,LX
            DO J=1,LY
             IOF=((((M-1)*LY+(J-1))*LX+(I-1))*IELEM+(JEL-1))*IELEM+IEL
             JOF=((((M1-1)*LY+(J-1))*LX+(I-1))*IELEM+(JEL-1))*IELEM+IEL
             E1=FUNKNO(L4+IOFZ+IOF,IG)
             FUNKNO(L4+IOFZ+IOF,IG)=FUNKNO(L4+IOFZ+JOF,IG)
             FUNKNO(L4+IOFZ+JOF,IG)=E1
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDIF
      ENDIF
*
  140 CONTINUE
  150 CONTINUE
*$OMP END PARALLEL DO
*----
*  DOMINO DIAGONAL LOOP
*----
      CALL XDDSET(XNI,IELEM*IELEM*LY*LZ*NPQ*NGEFF,0.0D0)
      CALL XDDSET(XNJ,IELEM*IELEM*LX*LZ*NPQ*NGEFF,0.0D0)
      CALL XDDSET(XNK,IELEM*IELEM*LX*LY*NPQ*NGEFF,0.0D0)
      DO 510 IDI=1,NMX+NMY+NMZ-2
*----
*  SET DOMINO SWAPPING INDICES
*----
      NMAX=MIN(NMX,IDI)*MIN(NMY,IDI)
      ALLOCATE(III(NMAX),JJJ(NMAX))
      NCEL=0
      DO 170 II=MAX(1,IDI-NMY-NMZ+2),MIN(NMX,IDI)
      DO 160 JJ=MAX(1,IDI-II-NMZ+2),MIN(NMY,IDI-II+1)
      NCEL=NCEL+1
      IF(NCEL.GT.NMAX) CALL XABORT('SNFT3D: NMAX OVERFLOW.')
      III(NCEL)=II
      JJJ(NCEL)=JJ
  160 CONTINUE
  170 CONTINUE
*----
*  MAIN SWAPPING LOOPS FOR SN FLUX CALCULATION
*----
*$OMP  PARALLEL DO
*$OMP1 PRIVATE(ITID,FLUX,M,IG,IMX,IMY,IMZ,Q,Q2,IOF,IER,II,JJ,KK,IEL,JEL)
*$OMP2 PRIVATE(I,J,K,L,VT) SHARED(FUNKNO,XNI,XNJ,XNK)
*$OMP3 REDUCTION(+:FLUX_G) COLLAPSE(3)
      DO 500 IG=1,NGEFF
      DO 490 IPQD=1,NPQD(IND)
      DO 480 ICEL=1,NCEL
#if defined(_OPENMP)
        ITID=omp_get_thread_num()
#else
        ITID=0
#endif
      II=III(ICEL)
      JJ=JJJ(ICEL)
      KK=IDI-II-JJ+2
*
      IF(IMPX.GT.5) WRITE(IUNOUT,600) ITID,NGIND(IG),IPQD,II,JJ,KK
      IF(.NOT.INCONV(IG)) GO TO 480
      M=INDANG(IPQD,IND)
      IF(W(M).EQ.0.0) GO TO 480
*
      IF(KK.EQ.1) THEN
        DO 210 IMX=1,MIN(IOMP,LX-(II-1)*IOMP)
        I=(II-1)*IOMP+IMX
        IF((IND.EQ.2).OR.(IND.EQ.3).OR.(IND.EQ.6).OR.(IND.EQ.7))I=LX+1-I
        DO 200 IMY=1,MIN(IOMP,LY-(JJ-1)*IOMP)
        J=(JJ-1)*IOMP+IMY
        IF((IND.EQ.3).OR.(IND.EQ.4).OR.(IND.EQ.7).OR.(IND.EQ.8))J=LY+1-J
        DO 190 JEL=1,IELEM
        DO 180 IEL=1,IELEM
        IOF=((((M-1)*LY+(J-1))*LX+(I-1))*IELEM+(JEL-1))*IELEM+IEL
        IF((IND.EQ.1).OR.(IND.EQ.2).OR.(IND.EQ.3).OR.(IND.EQ.4)) THEN
          XNK(IEL,JEL,I,J,IPQD,IG)=FUNKNO(L4+IOFZ+IOF,IG)*ZCODE(5)
        ELSE
          XNK(IEL,JEL,I,J,IPQD,IG)=FUNKNO(L4+IOFZ+IOF,IG)*ZCODE(6)
        ENDIF
  180   CONTINUE
  190   CONTINUE
  200   CONTINUE
  210   CONTINUE
      ENDIF
      IF(JJ.EQ.1) THEN
        DO 250 IMX=1,MIN(IOMP,LX-(II-1)*IOMP)
        I=(II-1)*IOMP+IMX
        IF((IND.EQ.2).OR.(IND.EQ.3).OR.(IND.EQ.6).OR.(IND.EQ.7))I=LX+1-I
        DO 240 IMZ=1,MIN(IOMP,LZ-(KK-1)*IOMP)
        K=(KK-1)*IOMP+IMZ
        IF((IND.EQ.5).OR.(IND.EQ.6).OR.(IND.EQ.7).OR.(IND.EQ.8))K=LZ+1-K
        DO 230 JEL=1,IELEM
        DO 220 IEL=1,IELEM
        IOF=((((M-1)*LZ+(K-1))*LX+(I-1))*IELEM+(JEL-1))*IELEM+IEL
        IF((IND.EQ.1).OR.(IND.EQ.2).OR.(IND.EQ.5).OR.(IND.EQ.6)) THEN
          XNJ(IEL,JEL,I,K,IPQD,IG)=FUNKNO(L4+IOFY+IOF,IG)*ZCODE(3)
        ELSE
          XNJ(IEL,JEL,I,K,IPQD,IG)=FUNKNO(L4+IOFY+IOF,IG)*ZCODE(4)
        ENDIF
  220   CONTINUE
  230   CONTINUE
  240   CONTINUE
  250   CONTINUE
      ENDIF
      IF(II.EQ.1) THEN
        DO 290 IMY=1,MIN(IOMP,LY-(JJ-1)*IOMP)
        J=(JJ-1)*IOMP+IMY
        IF((IND.EQ.3).OR.(IND.EQ.4).OR.(IND.EQ.7).OR.(IND.EQ.8))J=LY+1-J
        DO 280 IMZ=1,MIN(IOMP,LZ-(KK-1)*IOMP)
        K=(KK-1)*IOMP+IMZ
        IF((IND.EQ.5).OR.(IND.EQ.6).OR.(IND.EQ.7).OR.(IND.EQ.8))K=LZ+1-K
        DO 270 JEL=1,IELEM
        DO 260 IEL=1,IELEM
        IOF=((((M-1)*LZ+(K-1))*LY+(J-1))*IELEM+(JEL-1))*IELEM+IEL
        IF((IND.EQ.1).OR.(IND.EQ.4).OR.(IND.EQ.5).OR.(IND.EQ.8)) THEN
          XNI(IEL,JEL,J,K,IPQD,IG)=FUNKNO(L4+IOF,IG)*ZCODE(1)
        ELSE
          XNI(IEL,JEL,J,K,IPQD,IG)=FUNKNO(L4+IOF,IG)*ZCODE(2)
        ENDIF
  260   CONTINUE
  270   CONTINUE
  280   CONTINUE
  290   CONTINUE
      ENDIF
*---
      DO 470 IMX=1,MIN(IOMP,LX-(II-1)*IOMP)
      I=(II-1)*IOMP+IMX
      IF((IND.EQ.2).OR.(IND.EQ.3).OR.(IND.EQ.6).OR.(IND.EQ.7)) I=LX+1-I
      DO 460 IMY=1,MIN(IOMP,LY-(JJ-1)*IOMP)
      J=(JJ-1)*IOMP+IMY
      IF((IND.EQ.3).OR.(IND.EQ.4).OR.(IND.EQ.7).OR.(IND.EQ.8)) J=LY+1-J
      DO 450 IMZ=1,MIN(IOMP,LZ-(KK-1)*IOMP)
      K=(KK-1)*IOMP+IMZ
      IF((IND.EQ.5).OR.(IND.EQ.6).OR.(IND.EQ.7).OR.(IND.EQ.8)) K=LZ+1-K
      FLUX(:IELEM**3,:NSCT)=0.0D0
      IF(MAT(I,J,K).EQ.0) GO TO 480
      DO 340 IEL=1,IELEM**3
      Q(IEL)=0.0D0
      DO 330 L=1,NSCT
      IOF=((((K-1)*LY+(J-1))*LX+(I-1))*NSCT+(L-1))*IELEM**3+IEL
      Q(IEL)=Q(IEL)+QEXT(IOF,IG)*PL(L,M)/(4.0D0*PI)
  330 CONTINUE
  340 CONTINUE
      VT=VOL(I,J,K)*TOTAL(MAT(I,J,K),IG)
      CALL XDDSET(Q2,(IELEM**3)*((IELEM**3)+1),0.0D0)
*---------------------------------------
      IF(IELEM.EQ.1) THEN
      Q2(1,1)=2.0D0*ABS(DA(J,K,M))+2.0D0*ABS(DB(I,K,M))+
     1      2.0D0*ABS(DC(I,J,M))+VT
      Q2(1,2)=2.0D0*ABS(DA(J,K,M))*XNI(1,1,J,K,IPQD,IG)+
     1      2.0D0*ABS(DB(I,K,M))*XNJ(1,1,I,K,IPQD,IG)+
     2      2.0D0*ABS(DC(I,J,M))*XNK(1,1,I,J,IPQD,IG)+VOL(I,J,K)*Q(1)
      ELSE IF(IELEM.EQ.2) THEN
      Q2(1,1)=VT
      Q2(1,2)=CONST0*DA(J,K,M)
      Q2(1,3)=CONST0*DB(I,K,M)
      Q2(1,4)=CONST0*DC(I,J,M)
      Q2(2,2)=-VT-6.0D0*ABS(DA(J,K,M))
      Q2(2,5)=-CONST0*DB(I,K,M)
      Q2(2,6)=-CONST0*DC(I,J,M)
      Q2(3,3)=-VT-6.0D0*ABS(DB(I,K,M))
      Q2(3,5)=-CONST0*DA(J,K,M)
      Q2(3,7)=-CONST0*DC(I,J,M)
      Q2(4,4)=-VT-6.0D0*ABS(DC(I,J,M))
      Q2(4,6)=-CONST0*DA(J,K,M)
      Q2(4,7)=-CONST0*DB(I,K,M)
      Q2(5,5)=VT+6.0D0*ABS(DA(J,K,M))+6.0D0*ABS(DB(I,K,M))
      Q2(5,8)=CONST0*DC(I,J,M)
      Q2(6,6)=VT+6.0D0*ABS(DA(J,K,M))+6.0D0*ABS(DC(I,J,M))
      Q2(6,8)=CONST0*DB(I,K,M)
      Q2(7,7)=VT+6.0D0*ABS(DB(I,K,M))+6.0D0*ABS(DC(I,J,M))
      Q2(7,8)=CONST0*DA(J,K,M)
      Q2(8,8)=-VT-6.0D0*ABS(DA(J,K,M))-6.0D0*ABS(DB(I,K,M))
     1        -6.0D0*ABS(DC(I,J,M))
*---------------
      Q2(1,9)=VOL(I,J,K)*Q(1)
      Q2(2,9)=-VOL(I,J,K)*Q(2)+CONST0*DA(J,K,M)*XNI(1,1,J,K,IPQD,IG)
      Q2(3,9)=-VOL(I,J,K)*Q(3)+CONST0*DB(I,K,M)*XNJ(1,1,I,K,IPQD,IG)
      Q2(4,9)=-VOL(I,J,K)*Q(4)+CONST0*DC(I,J,M)*XNK(1,1,I,J,IPQD,IG)
      Q2(5,9)=VOL(I,J,K)*Q(5)-CONST0*DA(J,K,M)*XNI(2,1,J,K,IPQD,IG)
     1       -CONST0*DB(I,K,M)*XNJ(2,1,I,K,IPQD,IG)
      Q2(6,9)=VOL(I,J,K)*Q(6)-CONST0*DA(J,K,M)*XNI(1,2,J,K,IPQD,IG)
     1       -CONST0*DC(I,J,M)*XNK(2,1,I,J,IPQD,IG)
      Q2(7,9)=VOL(I,J,K)*Q(7)-CONST0*DB(I,K,M)*XNJ(1,2,I,K,IPQD,IG)
     1       -CONST0*DC(I,J,M)*XNK(1,2,I,J,IPQD,IG)
      Q2(8,9)=-VOL(I,J,K)*Q(8)+CONST0*DA(J,K,M)*XNI(2,2,J,K,IPQD,IG)
     1       +CONST0*DB(I,K,M)*XNJ(2,2,I,K,IPQD,IG)
     2       +CONST0*DC(I,J,M)*XNK(2,2,I,J,IPQD,IG)

      ELSE IF(IELEM.EQ.3) THEN
* UPPER DIAGONAL TERMS
      Q2(1,1)=-VT-2.0D0*ABS(DA(J,K,M))-2.0D0*ABS(DB(I,K,M))-
     1        2.0D0*ABS(DC(I,J,M))
      Q2(1,9)=-CONST1*ABS(DA(J,K,M))
      Q2(1,10)=-CONST1*ABS(DB(I,K,M))
      Q2(1,11)=-CONST1*ABS(DC(I,J,M))
      Q2(2,2)=VT+2.0D0*ABS(DB(I,K,M))+2.0D0*ABS(DC(I,J,M))
      Q2(2,9)=CONST2*DA(J,K,M)
      Q2(2,16)=CONST1*ABS(DC(I,J,M))
      Q2(2,17)=CONST1*ABS(DB(I,K,M))
      Q2(3,3)=VT+2.0D0*ABS(DA(J,K,M))+2.0D0*ABS(DC(I,J,M))
      Q2(3,10)=CONST2*DB(I,K,M)
      Q2(3,12)=CONST1*ABS(DA(J,K,M))
      Q2(3,15)=CONST1*ABS(DC(I,J,M))
      Q2(4,4)=VT+2.0D0*ABS(DA(J,K,M))+2.0D0*ABS(DB(I,K,M))
      Q2(4,11)=CONST2*DC(I,J,M)
      Q2(4,13)=CONST1*ABS(DA(J,K,M))
      Q2(4,19)=CONST1*ABS(DB(I,K,M))
      Q2(5,5)=-VT-2.0D0*ABS(DC(I,J,M))
      Q2(5,12)=-CONST2*DA(J,K,M)
      Q2(5,17)=-CONST2*DB(I,K,M)
      Q2(5,20)=-CONST1*ABS(DC(I,J,M))
      Q2(6,6)=-VT-2.0D0*ABS(DA(J,K,M))
      Q2(6,14)=-CONST1*ABS(DA(J,K,M))
      Q2(6,15)=-CONST2*DC(I,J,M)
      Q2(6,19)=-CONST2*DB(I,K,M)
      Q2(7,7)=-VT-2.0D0*ABS(DB(I,K,M))
      Q2(7,13)=-CONST2*DA(J,K,M)
      Q2(7,16)=-CONST2*DC(I,J,M)
      Q2(7,18)=-CONST1*ABS(DB(I,K,M))
      Q2(8,8)=VT
      Q2(8,14)=CONST2*DA(J,K,M)
      Q2(8,18)=CONST2*DB(I,K,M)
      Q2(8,20)=CONST2*DC(I,J,M)
      Q2(9,9)=-VT-1.0D1*ABS(DA(J,K,M))-2.0D0*ABS(DB(I,K,M))-
     1         2.0D0*ABS(DC(I,J,M))
      Q2(9,21)=-CONST1*ABS(DB(I,K,M))
      Q2(9,23)=-CONST1*ABS(DC(I,J,M))
      Q2(10,10)=-VT-2.0D0*ABS(DA(J,K,M))-1.0D1*ABS(DB(I,K,M))-
     1        2.0D0*ABS(DC(I,J,M))
      Q2(10,21)=-CONST1*ABS(DA(J,K,M))
      Q2(10,22)=-CONST1*ABS(DC(I,J,M))
      Q2(11,11)=-VT-2.0D0*ABS(DA(J,K,M))-2.0D0*ABS(DB(I,K,M))-
     1        1.0D1*ABS(DC(I,J,M))
      Q2(11,22)=-CONST1*ABS(DB(I,K,M))
      Q2(11,23)=-CONST1*ABS(DA(J,K,M))
      Q2(12,12)=VT+1.0D1*ABS(DA(J,K,M))+2.0D0*ABS(DC(I,J,M))
      Q2(12,21)=CONST2*DB(I,K,M)
      Q2(12,25)=CONST1*ABS(DC(I,J,M))
      Q2(13,13)=VT+1.0D1*ABS(DA(J,K,M))+2.0D0*ABS(DB(I,K,M))
      Q2(13,23)=CONST2*DC(I,J,M)
      Q2(13,26)=CONST1*ABS(DB(I,K,M))
      Q2(14,14)=-VT-1.0D1*ABS(DA(J,K,M))
      Q2(14,25)=-CONST2*DC(I,J,M)
      Q2(14,26)=-CONST2*DB(I,K,M)
      Q2(15,15)=VT+2.0D0*ABS(DA(J,K,M))+1.0D1*ABS(DC(I,J,M))
      Q2(15,22)=CONST2*DB(I,K,M)
      Q2(15,25)=CONST1*ABS(DA(J,K,M))
      Q2(16,16)=VT+2.0D0*ABS(DB(I,K,M))+1.0D1*ABS(DC(I,J,M))
      Q2(16,23)=CONST2*DA(J,K,M)
      Q2(16,24)=CONST1*ABS(DB(I,K,M))
      Q2(17,17)=VT+1.0D1*ABS(DB(I,K,M))+2.0D0*ABS(DC(I,J,M))
      Q2(17,21)=CONST2*DA(J,K,M)
      Q2(17,24)=CONST1*ABS(DC(I,J,M))
      Q2(18,18)=-VT-1.0D1*ABS(DB(I,K,M))
      Q2(18,24)=-CONST2*DC(I,J,M)
      Q2(18,26)=-CONST2*DA(J,K,M)
      Q2(19,19)=VT+2.0D0*ABS(DA(J,K,M))+1.0D1*ABS(DB(I,K,M))
      Q2(19,22)=CONST2*DC(I,J,M)
      Q2(19,26)=CONST1*ABS(DA(J,K,M))
      Q2(20,20)=-VT-1.0D1*ABS(DC(I,J,M))
      Q2(20,24)=-CONST2*DB(I,K,M)
      Q2(20,25)=-CONST2*DA(J,K,M)
      Q2(21,21)=-VT-1.0D1*ABS(DA(J,K,M))-1.0D1*ABS(DB(I,K,M))-
     1        2.0D0*ABS(DC(I,J,M))
      Q2(21,27)=-CONST1*ABS(DC(I,J,M))
      Q2(22,22)=-VT-1.0D1*ABS(DB(I,K,M))-1.0D1*ABS(DC(I,J,M))-
     1        2.0D0*ABS(DA(J,K,M))
      Q2(22,27)=-CONST1*ABS(DA(J,K,M))
      Q2(23,23)=-VT-1.0D1*ABS(DA(J,K,M))-2.0D0*ABS(DB(I,K,M))-
     1        1.0D1*ABS(DC(I,J,M))
      Q2(23,27)=-CONST1*ABS(DB(I,K,M))
      Q2(24,24)=VT+1.0D1*ABS(DB(I,K,M))+1.0D1*ABS(DC(I,J,M))
      Q2(24,27)=CONST2*DA(J,K,M)
      Q2(25,25)=VT+1.0D1*ABS(DA(J,K,M))+1.0D1*ABS(DC(I,J,M))
      Q2(25,27)=CONST2*DB(I,K,M)
      Q2(26,26)=VT+1.0D1*ABS(DA(J,K,M))+1.0D1*ABS(DB(I,K,M))
      Q2(26,27)=CONST2*DC(I,J,M)
      Q2(27,27)=-VT-1.0D1*ABS(DA(J,K,M))-1.0D1*ABS(DB(I,K,M))-
     1        1.0D1*ABS(DC(I,J,M))
*---------------------
      Q2(1,28)=-VOL(I,J,K)*Q(1)-
     1        2.0D0*ABS(DA(J,K,M))*XNI(1,1,J,K,IPQD,IG)-
     2        2.0D0*ABS(DB(I,K,M))*XNJ(1,1,I,K,IPQD,IG)-
     3        2.0D0*ABS(DC(I,J,M))*XNK(1,1,I,J,IPQD,IG)
      Q2(2,28)=VOL(I,J,K)*Q(2)+
     1        2.0D0*ABS(DB(I,K,M))*XNJ(2,1,I,K,IPQD,IG)+
     2        2.0D0*ABS(DC(I,J,M))*XNK(2,1,I,J,IPQD,IG)
      Q2(3,28)=VOL(I,J,K)*Q(3)+
     1        2.0D0*ABS(DA(J,K,M))*XNI(2,1,J,K,IPQD,IG)+
     2        2.0D0*ABS(DC(I,J,M))*XNK(1,2,I,J,IPQD,IG)
      Q2(4,28)=VOL(I,J,K)*Q(4)+
     1        2.0D0*ABS(DA(J,K,M))*XNI(1,2,J,K,IPQD,IG)+
     2        2.0D0*ABS(DB(I,K,M))*XNJ(1,2,I,K,IPQD,IG)
      Q2(5,28)=-VOL(I,J,K)*Q(5)-
     1        2.0D0*ABS(DC(I,J,M))*XNK(2,2,I,J,IPQD,IG)
      Q2(6,28)=-VOL(I,J,K)*Q(6)-
     1        2.0D0*ABS(DA(J,K,M))*XNI(2,2,J,K,IPQD,IG)
      Q2(7,28)=-VOL(I,J,K)*Q(7)-
     1        2.0D0*ABS(DB(I,K,M))*XNJ(2,2,I,K,IPQD,IG)
      Q2(8,28)=VOL(I,J,K)*Q(8)
      Q2(9,28)=-VOL(I,J,K)*Q(9)-
     1        CONST1*ABS(DA(J,K,M))*XNI(1,1,J,K,IPQD,IG)-
     2        2.0D0*ABS(DB(I,K,M))*XNJ(3,1,I,K,IPQD,IG)-
     3        2.0D0*ABS(DC(I,J,M))*XNK(3,1,I,J,IPQD,IG)
      Q2(10,28)=-VOL(I,J,K)*Q(10)-
     1        CONST1*ABS(DB(I,K,M))*XNJ(1,1,I,K,IPQD,IG)-
     2        2.0D0*ABS(DA(J,K,M))*XNI(3,1,J,K,IPQD,IG)-
     3        2.0D0*ABS(DC(I,J,M))*XNK(1,3,I,J,IPQD,IG)
      Q2(11,28)=-VOL(I,J,K)*Q(11)-
     1        2.0D0*ABS(DA(J,K,M))*XNI(1,3,J,K,IPQD,IG)-
     2        2.0D0*ABS(DB(I,K,M))*XNJ(1,3,I,K,IPQD,IG)-
     3        CONST1*ABS(DC(I,J,M))*XNK(1,1,I,J,IPQD,IG)
      Q2(12,28)=VOL(I,J,K)*Q(12)+
     1        CONST1*ABS(DA(J,K,M))*XNI(2,1,J,K,IPQD,IG)+
     2        2.0D0*ABS(DC(I,J,M))*XNK(3,2,I,J,IPQD,IG)
      Q2(13,28)=VOL(I,J,K)*Q(13)+
     1        2.0D0*ABS(DB(I,K,M))*XNJ(3,2,I,K,IPQD,IG)+
     2        CONST1*ABS(DA(J,K,M))*XNI(1,2,J,K,IPQD,IG)
      Q2(14,28)=-VOL(I,J,K)*Q(14)-
     1        CONST1*ABS(DA(J,K,M))*XNI(2,2,J,K,IPQD,IG)
      Q2(15,28)=VOL(I,J,K)*Q(15)+
     1        2.0D0*ABS(DA(J,K,M))*XNI(2,3,J,K,IPQD,IG)+
     2        CONST1*ABS(DC(I,J,M))*XNK(1,2,I,J,IPQD,IG)
      Q2(16,28)=VOL(I,J,K)*Q(16)+
     1        2.0D0*ABS(DB(I,K,M))*XNJ(2,3,I,K,IPQD,IG)+
     2        CONST1*ABS(DC(I,J,M))*XNK(2,1,I,J,IPQD,IG)
      Q2(17,28)=VOL(I,J,K)*Q(17)+
     1        CONST1*ABS(DB(I,K,M))*XNJ(2,1,I,K,IPQD,IG)+
     2        CONST1*ABS(DC(I,J,M))*XNK(2,3,I,J,IPQD,IG)
      Q2(18,28)=-VOL(I,J,K)*Q(18)-
     1        CONST1*ABS(DB(I,K,M))*XNJ(2,2,I,K,IPQD,IG)
      Q2(19,28)=VOL(I,J,K)*Q(19)+
     1        2.0D0*ABS(DA(J,K,M))*XNI(3,2,J,K,IPQD,IG)+
     2        CONST1*ABS(DB(I,K,M))*XNJ(1,2,I,K,IPQD,IG)
      Q2(20,28)=-VOL(I,J,K)*Q(20)-
     1        CONST1*ABS(DC(I,J,M))*XNK(2,2,I,J,IPQD,IG)
      Q2(21,28)=-VOL(I,J,K)*Q(21)-
     1        CONST1*ABS(DA(J,K,M))*XNI(3,1,J,K,IPQD,IG)-
     2        CONST1*ABS(DB(I,K,M))*XNJ(3,1,I,K,IPQD,IG)-
     3        2.0D0*CONST1*ABS(DC(I,J,M))*XNK(3,3,I,J,IPQD,IG)
      Q2(22,28)=-VOL(I,J,K)*Q(22)-
     1        2.0D0*ABS(DA(J,K,M))*XNI(3,3,J,K,IPQD,IG)
     2        -CONST1*ABS(DB(I,K,M))*XNJ(1,3,I,K,IPQD,IG)
     3        -CONST1*ABS(DC(I,J,M))*XNK(1,3,I,J,IPQD,IG)
      Q2(23,28)=-VOL(I,J,K)*Q(23)-
     1        CONST1*ABS(DA(J,K,M))*XNI(1,3,J,K,IPQD,IG)-
     2        2.0D0*ABS(DB(I,K,M))*XNJ(3,3,I,K,IPQD,IG)-
     3        CONST1*ABS(DC(I,J,M))*XNK(3,1,I,J,IPQD,IG)
      Q2(24,28)=VOL(I,J,K)*Q(24)+
     1        CONST1*ABS(DB(I,K,M))*XNJ(2,3,I,K,IPQD,IG)+
     2        CONST1*ABS(DC(I,J,M))*XNK(2,3,I,J,IPQD,IG)
      Q2(25,28)=VOL(I,J,K)*Q(25)+
     1        CONST1*ABS(DA(J,K,M))*XNI(2,3,J,K,IPQD,IG)+
     2        CONST1*ABS(DC(I,J,M))*XNK(3,2,I,J,IPQD,IG)
      Q2(26,28)=VOL(I,J,K)*Q(26)+
     1        CONST1*ABS(DA(J,K,M))*XNI(3,2,J,K,IPQD,IG)+
     2        CONST1*ABS(DB(I,K,M))*XNJ(3,2,I,K,IPQD,IG)
      Q2(27,28)=-VOL(I,J,K)*Q(27)-
     1        CONST1*ABS(DA(J,K,M))*XNI(3,3,J,K,IPQD,IG)-
     2        CONST1*ABS(DB(I,K,M))*XNJ(3,3,I,K,IPQD,IG)-
     3        CONST1*ABS(DC(I,J,M))*XNK(3,3,I,J,IPQD,IG)
      ENDIF
      DO 360 IEL=1,IELEM**3
      DO 350 JEL=IEL+1,IELEM**3
      Q2(JEL,IEL)=Q2(IEL,JEL)
 350  CONTINUE
 360  CONTINUE
      CALL ALSBD(IELEM**3,1,Q2,IER,IELEM**3)
      IF(IER.NE.0) CALL XABORT('SNFD13: SINGULAR MATRIX.')
      IF(IELEM.EQ.1) THEN
        IF(LFIXUP.AND.(Q2(1,2).LE.RLOG)) Q2(1,2)=0.0
        XNI(1,1,J,K,IPQD,IG)=2.0D0*Q2(1,2)-XNI(1,1,J,K,IPQD,IG)
        XNJ(1,1,I,K,IPQD,IG)=2.0D0*Q2(1,2)-XNJ(1,1,I,K,IPQD,IG)
        XNK(1,1,I,J,IPQD,IG)=2.0D0*Q2(1,2)-XNK(1,1,I,J,IPQD,IG)
        IF(LFIXUP.AND.(XNI(1,1,J,K,IPQD,IG).LE.RLOG))
     1    XNI(1,1,J,K,IPQD,IG)=0.0D0
        IF(LFIXUP.AND.(XNJ(1,1,I,K,IPQD,IG).LE.RLOG))
     1    XNJ(1,1,I,K,IPQD,IG)=0.0D0
        IF(LFIXUP.AND.(XNK(1,1,I,J,IPQD,IG).LE.RLOG))
     1    XNK(1,1,I,J,IPQD,IG)=0.0D0
      ELSE IF(IELEM.EQ.2) THEN
        XNI(1,1,J,K,IPQD,IG)=XNI(1,1,J,K,IPQD,IG)+SIGN(1.0,DU(M))*
     1    CONST0*Q2(2,9)
        XNI(1,2,J,K,IPQD,IG)=XNI(1,2,J,K,IPQD,IG)+SIGN(1.0,DU(M))*
     1    CONST0*Q2(6,9)
        XNI(2,1,J,K,IPQD,IG)=XNI(2,1,J,K,IPQD,IG)+SIGN(1.0,DU(M))*
     1    CONST0*Q2(5,9)
        XNI(2,2,J,K,IPQD,IG)=XNI(2,2,J,K,IPQD,IG)+SIGN(1.0,DU(M))*
     1    CONST0*Q2(8,9)
        XNJ(1,1,I,K,IPQD,IG)=XNJ(1,1,I,K,IPQD,IG)+SIGN(1.0,DE(M))*
     1    CONST0*Q2(3,9)
        XNJ(1,2,I,K,IPQD,IG)=XNJ(1,2,I,K,IPQD,IG)+SIGN(1.0,DE(M))*
     1    CONST0*Q2(7,9)
        XNJ(2,1,I,K,IPQD,IG)=XNJ(2,1,I,K,IPQD,IG)+SIGN(1.0,DE(M))*
     1    CONST0*Q2(5,9)
        XNJ(2,2,I,K,IPQD,IG)=XNJ(2,2,I,K,IPQD,IG)+SIGN(1.0,DE(M))*
     1    CONST0*Q2(8,9)
        XNK(1,1,I,J,IPQD,IG)=XNK(1,1,I,J,IPQD,IG)+SIGN(1.0,DZ(M))*
     1    CONST0*Q2(4,9)
        XNK(1,2,I,J,IPQD,IG)=XNK(1,2,I,J,IPQD,IG)+SIGN(1.0,DZ(M))*
     1    CONST0*Q2(7,9)
        XNK(2,1,I,J,IPQD,IG)=XNK(2,1,I,J,IPQD,IG)+SIGN(1.0,DZ(M))*
     1    CONST0*Q2(6,9)
        XNK(2,2,I,J,IPQD,IG)=XNK(2,2,I,J,IPQD,IG)+SIGN(1.0,DZ(M))*
     1    CONST0*Q2(8,9)
      ELSE IF(IELEM.EQ.3) THEN
        XNI(1,1,J,K,IPQD,IG)=2.0D0*Q2(1,28)+CONST1*Q2(9,28)-
     1    XNI(1,1,J,K,IPQD,IG)
        XNI(2,1,J,K,IPQD,IG)=2.0D0*Q2(3,28)+CONST1*Q2(12,28)-
     1    XNI(2,1,J,K,IPQD,IG)
        XNI(1,2,J,K,IPQD,IG)=2.0D0*Q2(4,28)+CONST1*Q2(13,28)-
     1    XNI(1,2,J,K,IPQD,IG)
        XNI(2,2,J,K,IPQD,IG)=2.0D0*Q2(6,28)+CONST1*Q2(14,28)-
     1    XNI(2,2,J,K,IPQD,IG)
        XNI(1,3,J,K,IPQD,IG)=2.0D0*Q2(11,28)+CONST1*Q2(23,28)-
     1    XNI(1,3,J,K,IPQD,IG)
        XNI(3,1,J,K,IPQD,IG)=2.0D0*Q2(10,28)+CONST1*Q2(21,28)-
     1    XNI(3,1,J,K,IPQD,IG)
        XNI(3,2,J,K,IPQD,IG)=2.0D0*Q2(19,28)+CONST1*Q2(26,28)-
     1    XNI(3,2,J,K,IPQD,IG)
        XNI(2,3,J,K,IPQD,IG)=2.0D0*Q2(15,28)+CONST1*Q2(25,28)-
     1    XNI(2,3,J,K,IPQD,IG)
        XNI(3,3,J,K,IPQD,IG)=2.0D0*Q2(22,28)+CONST1*Q2(27,28)-
     1    XNI(3,3,J,K,IPQD,IG)
        XNJ(1,1,I,K,IPQD,IG)=2.0D0*Q2(1,28)+CONST1*Q2(10,28)-
     1    XNJ(1,1,I,K,IPQD,IG)
        XNJ(2,1,I,K,IPQD,IG)=2.0D0*Q2(2,28)+CONST1*Q2(17,28)-
     1    XNJ(2,1,I,K,IPQD,IG)
        XNJ(1,2,I,K,IPQD,IG)=2.0D0*Q2(4,28)+CONST1*Q2(19,28)-
     1    XNJ(1,2,I,K,IPQD,IG)
        XNJ(2,2,I,K,IPQD,IG)=2.0D0*Q2(7,28)+CONST1*Q2(18,28)-
     1    XNJ(2,2,I,K,IPQD,IG)
        XNJ(3,1,I,K,IPQD,IG)=2.0D0*Q2(9,28)+CONST1*Q2(21,28)-
     1    XNJ(3,1,I,K,IPQD,IG)
        XNJ(1,3,I,K,IPQD,IG)=2.0D0*Q2(11,28)+CONST1*Q2(22,28)-
     1    XNJ(1,3,I,K,IPQD,IG)
        XNJ(3,2,I,K,IPQD,IG)=2.0D0*Q2(13,28)+CONST1*Q2(26,28)-
     1    XNJ(3,2,I,K,IPQD,IG)
        XNJ(2,3,I,K,IPQD,IG)=2.0D0*Q2(16,28)+CONST1*Q2(24,28)-
     1    XNJ(2,3,I,K,IPQD,IG)
        XNJ(3,3,I,K,IPQD,IG)=2.0D0*Q2(23,28)+CONST1*Q2(27,28)-
     1    XNJ(3,3,I,K,IPQD,IG)
        XNK(1,1,I,J,IPQD,IG)=2.0D0*Q2(1,28)+CONST1*Q2(11,28)-
     1    XNK(1,1,I,J,IPQD,IG)
        XNK(2,1,I,J,IPQD,IG)=2.0D0*Q2(2,28)+CONST1*Q2(16,28)-
     1    XNK(2,1,I,J,IPQD,IG)
        XNK(1,2,I,J,IPQD,IG)=2.0D0*Q2(3,28)+CONST1*Q2(15,28)-
     1    XNK(1,2,I,J,IPQD,IG)
        XNK(2,2,I,J,IPQD,IG)=2.0D0*Q2(5,28)+CONST1*Q2(20,28)-
     1    XNK(2,2,I,J,IPQD,IG)
        XNK(3,1,I,J,IPQD,IG)=2.0D0*Q2(9,28)+CONST1*Q2(23,28)-
     1    XNK(3,1,I,J,IPQD,IG)
        XNK(1,3,I,J,IPQD,IG)=2.0D0*Q2(10,28)+CONST1*Q2(22,28)-
     1    XNK(1,3,I,J,IPQD,IG)
        XNK(3,2,I,J,IPQD,IG)=2.0D0*Q2(12,28)+CONST1*Q2(25,28)-
     1    XNK(3,2,I,J,IPQD,IG)
        XNK(2,3,I,J,IPQD,IG)=2.0D0*Q2(17,28)+CONST1*Q2(24,28)-
     1    XNK(2,3,I,J,IPQD,IG)
        XNK(3,3,I,J,IPQD,IG)=2.0D0*Q2(21,28)+CONST1*Q2(27,28)-
     1    XNK(3,3,I,J,IPQD,IG)
      ENDIF
      DO 380 L=1,NSCT
      DO 370 IEL=1,IELEM**3
      FLUX(IEL,L)=FLUX(IEL,L)+W(M)*Q2(IEL,IELEM**3+1)*PL(L,M)
  370 CONTINUE
  380 CONTINUE
      IF((II-1)*IOMP+IMX.EQ.LX) THEN
        DO 400 JEL=1,IELEM
        DO 390 IEL=1,IELEM
        IOF=((((M-1)*LZ+(K-1))*LY+(J-1))*IELEM+(JEL-1))*IELEM+IEL
        FUNKNO(L4+IOF,IG)=REAL(XNI(IEL,JEL,J,K,IPQD,IG))
  390   CONTINUE
  400   CONTINUE
      ENDIF
      IF((JJ-1)*IOMP+IMY.EQ.LY) THEN
        DO 420 JEL=1,IELEM
        DO 410 IEL=1,IELEM
        IOF=((((M-1)*LZ+(K-1))*LX+(I-1))*IELEM+(JEL-1))*IELEM+IEL
        FUNKNO(L4+IOFY+IOF,IG)=REAL(XNJ(IEL,JEL,I,K,IPQD,IG))
  410   CONTINUE
  420   CONTINUE
      ENDIF
      IF((KK-1)*IOMP+IMZ.EQ.LZ) THEN
        DO 440 JEL=1,IELEM
        DO 430 IEL=1,IELEM
        IOF=((((M-1)*LY+(J-1))*LX+(I-1))*IELEM+(JEL-1))*IELEM+IEL
        FUNKNO(L4+IOFZ+IOF,IG)=REAL(XNK(IEL,JEL,I,J,IPQD,IG))
  430   CONTINUE
  440   CONTINUE
      ENDIF
      FLUX_G(:,:,I,J,K,IG)=FLUX_G(:,:,I,J,K,IG)+FLUX(:,:)
  450 CONTINUE
  460 CONTINUE
  470 CONTINUE
*---
  480 CONTINUE
  490 CONTINUE
  500 CONTINUE
*$OMP END PARALLEL DO
      DEALLOCATE(JJJ,III)
  510 CONTINUE
  520 CONTINUE
      DO 530 IG=1,NGEFF
        IF(.NOT.INCONV(IG)) GO TO 530
        FUNKNO(:L4,IG)=
     1  RESHAPE(REAL(FLUX_G(:IELEM**3,:NSCT,:LX,:LY,:LZ,IG)), (/L4/))
  530 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(FLUX_G,FLUX,XNK,XNJ,XNI)
      DEALLOCATE(INDANG)
      RETURN
*
  600 FORMAT(16H SNFD13: thread=,I8,12H --->(group=,I4,7H angle=,I4,
     1 11H macrocell=,3I5,1H))
      END
