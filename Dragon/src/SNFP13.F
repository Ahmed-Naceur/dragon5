*DECK SNFP13
      SUBROUTINE SNFP13(NUN,NGRP,NGEFF,IMPX,INCONV,INGIND,LX,LY,LZ,
     1 IELEM,IBFP,NMAT,NPQ,NSCT,MAT,VOL,TOTAL,ESTOPW,NCODE,ZCODE,
     2 DELTAE,QEXT,LFIXUP,DU,DE,DZ,W,MRMX,MRMY,MRMZ,DC,DB,DA,PL,FUNKNO,
     3 DCUTOFF)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform one inner iteration for solving SN equations in 3D Cartesian
* geometry for the HODD method. Energy-angle multithreading. Albedo
* boundary conditions. Boltzmann-Fokker-Planck (BFP) discretization.
*
*Copyright:
* Copyright (C) 2021 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NUN     total number of unknowns in vector FUNKNO.
* NGRP    number of energy groups.
* NGEFF   dimension of arrays INCONV and INGIND.
* IMPX    print flag (equal to zero for no print).
* INCONV  energy group convergence flag (set to .FALSE. if converged).
* INGIND  energy group index assign to 1:NGEFF arrays.
* LX      number of meshes along X axis.
* LY      number of meshes along Y axis.
* LZ      number of meshes along Z axis.
* IELEM   degree of the spatial approximation (=1 diamond scheme).
* IBFP    type of energy proparation relation:
*         =1 Galerkin type;
*         =2 heuristic Przybylski and Ligou type.
* NMAT    number of material mixtures.
* NPQ     number of SN directions in height octants.
* NSCT    maximum number of spherical harmonics moments of the flux.
* MAT     material mixture index in each region.
* VOL     volumes of each region.
* TOTAL   macroscopic total cross sections.
* ESTOPW  stopping power.
* NCODE   boundary condition indices.
* ZCODE   albedos.
* DELTAE  energy group width in MeV.
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
* DCUTOFF energy deposition under the energy cutoff.
*
*-----------------------------------------------------------------------
*
#if defined(_OPENMP)
      USE omp_lib
#endif
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NUN,NGRP,NGEFF,IMPX,INGIND(NGEFF),LX,LY,LZ,IELEM,IBFP,
     1 NMAT,NPQ,NSCT,MAT(LX,LY,LZ),NCODE(6),MRMX(NPQ),MRMY(NPQ),
     2 MRMZ(NPQ)
      LOGICAL INCONV(NGEFF)
      REAL VOL(LX,LY,LZ),TOTAL(0:NMAT,NGRP),ESTOPW(0:NMAT,2,NGRP),
     1 ZCODE(6),DELTAE(NGRP),QEXT(NUN,NGRP),DU(NPQ),DE(NPQ),DZ(NPQ),
     2 W(NPQ),DC(LX,LY,NPQ),DB(LX,LZ,NPQ),DA(LY,LZ,NPQ),PL(NSCT,NPQ),
     3 FUNKNO(NUN,NGRP),DCUTOFF(LX,LY,LZ)
      LOGICAL LFIXUP
*----
*  LOCAL VARIABLES
*----
      INTEGER NPQD(8),IIND(8)
      PARAMETER(IUNOUT=6,RLOG=1.0E-8,PI=3.141592654)
      DOUBLE PRECISION VT,CONST0,CONST1,CONST2,FN,FD,G1,G2,Q(IELEM**3),
     1 Q2(IELEM**3,IELEM**3+1),FEP(IELEM**3),XNK(IELEM,IELEM)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: INDANG
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: FLUX,FLUX0
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: FLUX_G,
     1 FLUX0_G
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: XNI
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: XNJ
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(INDANG(NPQ,8))
      ALLOCATE(XNI(IELEM,IELEM,LY,LZ),XNJ(IELEM,IELEM,LZ))
      ALLOCATE(FLUX(IELEM**3,NSCT,LX,LY,LZ),
     1 FLUX0(IELEM**3,NPQ,LX,LY,LZ))
      ALLOCATE(FLUX_G(IELEM**3,NSCT,LX,LY,LZ,NGEFF),
     1 FLUX0_G(IELEM**3,NPQ,LX,LY,LZ,NGEFF))
*----
*  DEFINITION OF CONSTANTS.
*----
      L4=IELEM*IELEM*IELEM*LX*LY*LZ*NSCT
      L6=IELEM*IELEM*IELEM*LX*LY*LZ*NPQ
      L5=NUN-L6
      IOFY=IELEM*IELEM*LY*LZ*NPQ
      IOFZ=IELEM*IELEM*(LX+LY)*LZ*NPQ
      CONST0=2.0D0*DSQRT(3.0D0)
      CONST1=2.0D0*DSQRT(5.0D0)
      CONST2=2.0D0*DSQRT(15.0D0)
*----
*  PARAMETER VALIDATION.
*----
      IF(IELEM.GT.3) CALL XABORT('SNFP13: INVALID IELEM (DIAM) VALUE. '
     1 //'CHECK INPUT DATA FILE.')
      FLUX_G(:IELEM**3,:NSCT,:LX,:LY,:LZ,:NGEFF)=0.0D0
      FLUX0_G(:IELEM**3,:NPQ,:LX,:LY,:LZ,:NGEFF)=0.0D0
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
      DO 420 JND=1,8
      IND=IIND(JND)
*----
*  PRELIMINARY LOOPS FOR SETTING BOUNDARY CONDITIONS.
*----
*$OMP  PARALLEL DO
*$OMP1 PRIVATE(M,IG,VU,VE,VZ,M1,E1,IOF,JOF,IEL,JEL,I,J,K)
*$OMP2 SHARED(FUNKNO) COLLAPSE(2)
      DO 150 IGEFF=1,NGEFF
      DO 140 IPQD=1,NPQD(IND)
      IF(.NOT.INCONV(IGEFF)) GO TO 140
      M=INDANG(IPQD,IND)
      IG=INGIND(IGEFF)
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
*  MAIN SWAPPING LOOPS FOR SN FLUX CALCULATION
*----
*$OMP  PARALLEL DO
*$OMP1 PRIVATE(ITID,FLUX,M,IG,XNI,XNJ,XNK,Q,Q2,IOF,IER,IX,IY,IZ,IEL,JEL)
*$OMP1 PRIVATE(FEP,FLUX0,FN,FD,G1,G2)
*$OMP2 PRIVATE(I,J,K,L,VT) SHARED(FUNKNO) REDUCTION(+:FLUX_G,FLUX0_G)
*$OMP3 COLLAPSE(2)
      DO 410 IGEFF=1,NGEFF
      DO 400 IPQD=1,NPQD(IND)
#if defined(_OPENMP)
        ITID=omp_get_thread_num()
#else
        ITID=0
#endif
      IF(IMPX.GT.5) WRITE(IUNOUT,500) ITID,IGEFF,IPQD
      IF(.NOT.INCONV(IGEFF)) GO TO 400
      M=INDANG(IPQD,IND)
      IG=INGIND(IGEFF)
      FLUX(:IELEM**3,:NSCT,:LX,:LY,:LZ)=0.0D0
      FLUX0(:IELEM**3,:NPQ,:LX,:LY,:LZ)=0.0D0
      IF(W(M).EQ.0.0) GO TO 400
*----
*  LOOP OVER X-, Y- AND Z-DIRECTED AXES.
*----
      DO 350 IX=1,LX
      I=IX
      IF((IND.EQ.2).OR.(IND.EQ.3).OR.(IND.EQ.6).OR.(IND.EQ.7)) I=LX+1-I
      DO 310 IY=1,LY
      J=IY
      IF((IND.EQ.3).OR.(IND.EQ.4).OR.(IND.EQ.7).OR.(IND.EQ.8)) J=LY+1-J
      DO 170 JEL=1,IELEM
      DO 160 IEL=1,IELEM
      IOF=((((M-1)*LY+(J-1))*LX+(I-1))*IELEM+(JEL-1))*IELEM+IEL
      IF((IND.EQ.1).OR.(IND.EQ.2).OR.(IND.EQ.3).OR.(IND.EQ.4)) THEN
        XNK(IEL,JEL)=FUNKNO(L4+IOFZ+IOF,IG)*ZCODE(5)
      ELSE
        XNK(IEL,JEL)=FUNKNO(L4+IOFZ+IOF,IG)*ZCODE(6)
      ENDIF
  160 CONTINUE
  170 CONTINUE
      DO 280 IZ=1,LZ
      K=IZ
      IF((IND.EQ.5).OR.(IND.EQ.6).OR.(IND.EQ.7).OR.(IND.EQ.8)) K=LZ+1-K
      IF(IY.EQ.1) THEN
        DO 190 JEL=1,IELEM
        DO 180 IEL=1,IELEM
        IOF=((((M-1)*LZ+(K-1))*LX+(I-1))*IELEM+(JEL-1))*IELEM+IEL
        IF((IND.EQ.1).OR.(IND.EQ.2).OR.(IND.EQ.5).OR.(IND.EQ.6)) THEN
          XNJ(IEL,JEL,K)=FUNKNO(L4+IOFY+IOF,IG)*ZCODE(3)
        ELSE
          XNJ(IEL,JEL,K)=FUNKNO(L4+IOFY+IOF,IG)*ZCODE(4)
        ENDIF
  180   CONTINUE
  190   CONTINUE
      ENDIF
      IF(IX.EQ.1) THEN
        DO 210 JEL=1,IELEM
        DO 200 IEL=1,IELEM
        IOF=((((M-1)*LZ+(K-1))*LY+(J-1))*IELEM+(JEL-1))*IELEM+IEL
        IF((IND.EQ.1).OR.(IND.EQ.4).OR.(IND.EQ.5).OR.(IND.EQ.8)) THEN
          XNI(IEL,JEL,J,K)=FUNKNO(L4+IOF,IG)*ZCODE(1)
        ELSE
          XNI(IEL,JEL,J,K)=FUNKNO(L4+IOF,IG)*ZCODE(2)
        ENDIF
  200   CONTINUE
  210   CONTINUE
      ENDIF
      IBM=MAT(I,J,K)
      IF(IBM.EQ.0) GO TO 280
      DO 230 IEL=1,IELEM**3
      Q(IEL)=0.0D0
      DO 220 L=1,NSCT
      IOF=((((K-1)*LY+(J-1))*LX+(I-1))*NSCT+(L-1))*IELEM**3+IEL
      Q(IEL)=Q(IEL)+QEXT(IOF,IG)*PL(L,M)/(4.0D0*PI)
  220 CONTINUE
      IOF=((((K-1)*LY+(J-1))*LX+(I-1))*NPQ+(M-1))*IELEM**3+IEL 
      FEP(IEL)=QEXT(L5+IOF,IG)
  230 CONTINUE
      IF(IBFP.EQ.1) THEN
*       Galerkin type
        FN=2.0D0*ESTOPW(IBM,1,IG)/DELTAE(IG)
        FD=(ESTOPW(IBM,1,IG)+ESTOPW(IBM,2,IG))/DELTAE(IG)
        G1=ESTOPW(IBM,1,IG)/ESTOPW(IBM,2,IG)
        G2=(ESTOPW(IBM,1,IG)+ESTOPW(IBM,2,IG))/ESTOPW(IBM,2,IG)
      ELSE
*       Przybylski and Ligou type
        FN=(ESTOPW(IBM,1,IG)+ESTOPW(IBM,2,IG))/DELTAE(IG)
        FD=2.0D0*ESTOPW(IBM,2,IG)/DELTAE(IG)
        G1=1.0D0
        G2=2.0D0
      ENDIF
      VT=VOL(I,J,K)*(FD+TOTAL(IBM,IG))
      CALL XDDSET(Q2,(IELEM**3)*((IELEM**3)+1),0.0D0)
*---------------------------------------
      IF(IELEM.EQ.1) THEN
      Q2(1,1)=2.0D0*ABS(DA(J,K,M))+2.0D0*ABS(DB(I,K,M))+
     1        2.0D0*ABS(DC(I,J,M))+VT
      Q2(1,2)=2.0D0*ABS(DA(J,K,M))*XNI(1,1,J,K)+
     1        2.0D0*ABS(DB(I,K,M))*XNJ(1,1,K)+
     2        2.0D0*ABS(DC(I,J,M))*XNK(1,1)+FN*FEP(1)*VOL(I,J,K)+
     3        VOL(I,J,K)*Q(1)
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
      Q2(1,9)=FN*FEP(1)*VOL(I,J,K)+VOL(I,J,K)*Q(1)
      Q2(2,9)=-FN*FEP(2)*VOL(I,J,K)-VOL(I,J,K)*Q(2)+CONST0*DA(J,K,M)*
     1        XNI(1,1,J,K)
      Q2(3,9)=-FN*FEP(3)*VOL(I,J,K)-VOL(I,J,K)*Q(3)+CONST0*DB(I,K,M)*
     1        XNJ(1,1,K)
      Q2(4,9)=-FN*FEP(4)*VOL(I,J,K)-VOL(I,J,K)*Q(4)+CONST0*DC(I,J,M)*
     1        XNK(1,1)
      Q2(5,9)=FN*FEP(5)*VOL(I,J,K)+VOL(I,J,K)*Q(5)-CONST0*DA(J,K,M)*
     1        XNI(2,1,J,K)-CONST0*DB(I,K,M)*XNJ(2,1,K)
      Q2(6,9)=FN*FEP(6)*VOL(I,J,K)+VOL(I,J,K)*Q(6)-CONST0*DA(J,K,M)*
     1        XNI(1,2,J,K)-CONST0*DC(I,J,M)*XNK(2,1)
      Q2(7,9)=FN*FEP(7)*VOL(I,J,K)+VOL(I,J,K)*Q(7)-CONST0*DB(I,K,M)*
     1        XNJ(1,2,K)-CONST0*DC(I,J,M)*XNK(1,2)
      Q2(8,9)=-FN*FEP(8)*VOL(I,J,K)-VOL(I,J,K)*Q(8)+CONST0*DA(J,K,M)*
     1        XNI(2,2,J,K)+CONST0*DB(I,K,M)*XNJ(2,2,K)
     2        +CONST0*DC(I,J,M)*XNK(2,2)

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
      Q2(1,28)=-FN*FEP(1)*VOL(I,J,K)-VOL(I,J,K)*Q(1)-2.0D0*
     1        ABS(DA(J,K,M))*XNI(1,1,J,K)-2.0D0*ABS(DB(I,K,M))*
     2        XNJ(1,1,K)-2.0D0*ABS(DC(I,J,M))*XNK(1,1)
      Q2(2,28)=FN*FEP(2)*VOL(I,J,K)+VOL(I,J,K)*Q(2)+2.0D0*
     1        ABS(DB(I,K,M))*XNJ(2,1,K)+2.0D0*ABS(DC(I,J,M))*XNK(2,1)
      Q2(3,28)=FN*FEP(3)*VOL(I,J,K)+VOL(I,J,K)*Q(3)+2.0D0*
     1        ABS(DA(J,K,M))*XNI(2,1,J,K)+2.0D0*ABS(DC(I,J,M))*XNK(1,2)
      Q2(4,28)=FN*FEP(4)*VOL(I,J,K)+VOL(I,J,K)*Q(4)+2.0D0*
     1        ABS(DA(J,K,M))*XNI(1,2,J,K)+2.0D0*ABS(DB(I,K,M))*
     2        XNJ(1,2,K)
      Q2(5,28)=-FN*FEP(5)*VOL(I,J,K)-VOL(I,J,K)*Q(5)-2.0D0*
     1        ABS(DC(I,J,M))*XNK(2,2)
      Q2(6,28)=-FN*FEP(6)*VOL(I,J,K)-VOL(I,J,K)*Q(6)-2.0D0*
     1        ABS(DA(J,K,M))*XNI(2,2,J,K)
      Q2(7,28)=-FN*FEP(7)*VOL(I,J,K)-VOL(I,J,K)*Q(7)-2.0D0*
     1        ABS(DB(I,K,M))*XNJ(2,2,K)
      Q2(8,28)=FN*FEP(8)*VOL(I,J,K)+VOL(I,J,K)*Q(8)
      Q2(9,28)=-FN*FEP(9)*VOL(I,J,K)-VOL(I,J,K)*Q(9)-CONST1*
     1        ABS(DA(J,K,M))*XNI(1,1,J,K)-2.0D0*ABS(DB(I,K,M))*
     2        XNJ(3,1,K)-2.0D0*ABS(DC(I,J,M))*XNK(3,1)
      Q2(10,28)=-FN*FEP(10)*VOL(I,J,K)-VOL(I,J,K)*Q(10)-CONST1*
     1        ABS(DB(I,K,M))*XNJ(1,1,K)-2.0D0*ABS(DA(J,K,M))*
     2        XNI(3,1,J,K)-2.0D0*ABS(DC(I,J,M))*XNK(1,3)
      Q2(11,28)=-FN*FEP(11)*VOL(I,J,K)-VOL(I,J,K)*Q(11)-2.0D0*
     1        ABS(DA(J,K,M))*XNI(1,3,J,K)-2.0D0*ABS(DB(I,K,M))*
     2        XNJ(1,3,K)-CONST1*ABS(DC(I,J,M))*XNK(1,1)
      Q2(12,28)=FN*FEP(12)*VOL(I,J,K)+VOL(I,J,K)*Q(12)+CONST1*
     1        ABS(DA(J,K,M))*XNI(2,1,J,K)+2.0D0*ABS(DC(I,J,M))*XNK(3,2)
      Q2(13,28)=FN*FEP(13)*VOL(I,J,K)+VOL(I,J,K)*Q(13)+2.0D0*
     1        ABS(DB(I,K,M))*XNJ(3,2,K)+CONST1*ABS(DA(J,K,M))*
     2        XNI(1,2,J,K)
      Q2(14,28)=-FN*FEP(14)*VOL(I,J,K)-VOL(I,J,K)*Q(14)-CONST1*
     1        ABS(DA(J,K,M))*XNI(2,2,J,K)
      Q2(15,28)=FN*FEP(15)*VOL(I,J,K)+VOL(I,J,K)*Q(15)+2.0D0*
     1        ABS(DA(J,K,M))*XNI(2,3,J,K)+CONST1*ABS(DC(I,J,M))*XNK(1,2)
      Q2(16,28)=FN*FEP(16)*VOL(I,J,K)+VOL(I,J,K)*Q(16)+2.0D0*
     1        ABS(DB(I,K,M))*XNJ(2,3,K)+CONST1*ABS(DC(I,J,M))*XNK(2,1)
      Q2(17,28)=FN*FEP(17)*VOL(I,J,K)+VOL(I,J,K)*Q(17)+CONST1*
     1        ABS(DB(I,K,M))*XNJ(2,1,K)+CONST1*ABS(DC(I,J,M))*XNK(2,3)
      Q2(18,28)=-FN*FEP(18)*VOL(I,J,K)-VOL(I,J,K)*Q(18)-CONST1*
     1        ABS(DB(I,K,M))*XNJ(2,2,K)
      Q2(19,28)=FN*FEP(19)*VOL(I,J,K)+VOL(I,J,K)*Q(19)+2.0D0*
     1        ABS(DA(J,K,M))*XNI(3,2,J,K)+CONST1*ABS(DB(I,K,M))*
     2        XNJ(1,2,K)
      Q2(20,28)=-FN*FEP(20)*VOL(I,J,K)-VOL(I,J,K)*Q(20)-CONST1*
     1        ABS(DC(I,J,M))*XNK(2,2)
      Q2(21,28)=-FN*FEP(21)*VOL(I,J,K)-VOL(I,J,K)*Q(21)-
     1        CONST1*ABS(DA(J,K,M))*XNI(3,1,J,K)-
     2        CONST1*ABS(DB(I,K,M))*XNJ(3,1,K)-
     3        2.0D0*CONST1*ABS(DC(I,J,M))*XNK(3,3)
      Q2(22,28)=-FN*FEP(22)*VOL(I,J,K)-VOL(I,J,K)*Q(22)-2.0D0*
     1        ABS(DA(J,K,M))*XNI(3,3,J,K)-CONST1*ABS(DB(I,K,M))*
     2        XNJ(1,3,K)-CONST1*ABS(DC(I,J,M))*XNK(1,3)
      Q2(23,28)=-FN*FEP(23)*VOL(I,J,K)-VOL(I,J,K)*Q(23)-
     1        CONST1*ABS(DA(J,K,M))*XNI(1,3,J,K)-
     2        2.0D0*ABS(DB(I,K,M))*XNJ(3,3,K)-
     3        CONST1*ABS(DC(I,J,M))*XNK(3,1)
      Q2(24,28)=FN*FEP(24)*VOL(I,J,K)+VOL(I,J,K)*Q(24)+CONST1*
     1        ABS(DB(I,K,M))*XNJ(2,3,K)+CONST1*ABS(DC(I,J,M))*XNK(2,3)
      Q2(25,28)=FN*FEP(25)*VOL(I,J,K)+VOL(I,J,K)*Q(25)+CONST1*
     1        ABS(DA(J,K,M))*XNI(2,3,J,K)+CONST1*ABS(DC(I,J,M))*XNK(3,2)
      Q2(26,28)=FN*FEP(26)*VOL(I,J,K)+VOL(I,J,K)*Q(26)+CONST1*
     1        ABS(DA(J,K,M))*XNI(3,2,J,K)+CONST1*ABS(DB(I,K,M))*
     2        XNJ(3,2,K)
      Q2(27,28)=-FN*FEP(27)*VOL(I,J,K)-VOL(I,J,K)*Q(27)-
     1        CONST1*ABS(DA(J,K,M))*XNI(3,3,J,K)-
     2        CONST1*ABS(DB(I,K,M))*XNJ(3,3,K)-
     3        CONST1*ABS(DC(I,J,M))*XNK(3,3)
      ENDIF
      DO 250 IEL=1,IELEM**3
      DO 240 JEL=IEL+1,IELEM**3
      Q2(JEL,IEL)=Q2(IEL,JEL)
 240  CONTINUE
 250  CONTINUE
      CALL ALSBD(IELEM**3,1,Q2,IER,IELEM**3)
      IF(IER.NE.0) CALL XABORT('SNFP13: SINGULAR MATRIX.')
      IF(IELEM.EQ.1) THEN
        IF(LFIXUP.AND.(Q2(1,2).LE.RLOG)) Q2(1,2)=0.0
        XNI(1,1,J,K)=2.0D0*Q2(1,2)-XNI(1,1,J,K)
        XNJ(1,1,K)=2.0D0*Q2(1,2)-XNJ(1,1,K)
        XNK(1,1)=2.0D0*Q2(1,2)-XNK(1,1)
        IF(LFIXUP.AND.(XNI(1,1,J,K).LE.RLOG)) XNI(1,1,J,K)=0.0
        IF(LFIXUP.AND.(XNJ(1,1,K).LE.RLOG)) XNJ(1,1,K)=0.0
        IF(LFIXUP.AND.(XNK(1,1).LE.RLOG)) XNK(1,1)=0.0
      ELSE IF(IELEM.EQ.2) THEN
        XNI(1,1,J,K)=XNI(1,1,J,K)+SIGN(1.0,DU(M))*CONST0*Q2(2,9)
        XNI(1,2,J,K)=XNI(1,2,J,K)+SIGN(1.0,DU(M))*CONST0*Q2(6,9)
        XNI(2,1,J,K)=XNI(2,1,J,K)+SIGN(1.0,DU(M))*CONST0*Q2(5,9)
        XNI(2,2,J,K)=XNI(2,2,J,K)+SIGN(1.0,DU(M))*CONST0*Q2(8,9)
        XNJ(1,1,K)=XNJ(1,1,K)+SIGN(1.0,DE(M))*CONST0*Q2(3,9)
        XNJ(1,2,K)=XNJ(1,2,K)+SIGN(1.0,DE(M))*CONST0*Q2(7,9)
        XNJ(2,1,K)=XNJ(2,1,K)+SIGN(1.0,DE(M))*CONST0*Q2(5,9)
        XNJ(2,2,K)=XNJ(2,2,K)+SIGN(1.0,DE(M))*CONST0*Q2(8,9)
        XNK(1,1)=XNK(1,1)+SIGN(1.0,DZ(M))*CONST0*Q2(4,9)
        XNK(1,2)=XNK(1,2)+SIGN(1.0,DZ(M))*CONST0*Q2(7,9)
        XNK(2,1)=XNK(2,1)+SIGN(1.0,DZ(M))*CONST0*Q2(6,9)
        XNK(2,2)=XNK(2,2)+SIGN(1.0,DZ(M))*CONST0*Q2(8,9)
      ELSE IF(IELEM.EQ.3) THEN
        XNI(1,1,J,K)=2.0D0*Q2(1,28)+CONST1*Q2(9,28)-XNI(1,1,J,K)
        XNI(2,1,J,K)=2.0D0*Q2(3,28)+CONST1*Q2(12,28)-XNI(2,1,J,K)
        XNI(1,2,J,K)=2.0D0*Q2(4,28)+CONST1*Q2(13,28)-XNI(1,2,J,K)
        XNI(2,2,J,K)=2.0D0*Q2(6,28)+CONST1*Q2(14,28)-XNI(2,2,J,K)
        XNI(1,3,J,K)=2.0D0*Q2(11,28)+CONST1*Q2(23,28)-XNI(1,3,J,K)
        XNI(3,1,J,K)=2.0D0*Q2(10,28)+CONST1*Q2(21,28)-XNI(3,1,J,K)
        XNI(3,2,J,K)=2.0D0*Q2(19,28)+CONST1*Q2(26,28)-XNI(3,2,J,K)
        XNI(2,3,J,K)=2.0D0*Q2(15,28)+CONST1*Q2(25,28)-XNI(2,3,J,K)
        XNI(3,3,J,K)=2.0D0*Q2(22,28)+CONST1*Q2(27,28)-XNI(3,3,J,K)
        XNJ(1,1,K)=2.0D0*Q2(1,28)+CONST1*Q2(10,28)-XNJ(1,1,K)
        XNJ(2,1,K)=2.0D0*Q2(2,28)+CONST1*Q2(17,28)-XNJ(2,1,K)
        XNJ(1,2,K)=2.0D0*Q2(4,28)+CONST1*Q2(19,28)-XNJ(1,2,K)
        XNJ(2,2,K)=2.0D0*Q2(7,28)+CONST1*Q2(18,28)-XNJ(2,2,K)
        XNJ(3,1,K)=2.0D0*Q2(9,28)+CONST1*Q2(21,28)-XNJ(3,1,K)
        XNJ(1,3,K)=2.0D0*Q2(11,28)+CONST1*Q2(22,28)-XNJ(1,3,K)
        XNJ(3,2,K)=2.0D0*Q2(13,28)+CONST1*Q2(26,28)-XNJ(3,2,K)
        XNJ(2,3,K)=2.0D0*Q2(16,28)+CONST1*Q2(24,28)-XNJ(2,3,K)
        XNJ(3,3,K)=2.0D0*Q2(23,28)+CONST1*Q2(27,28)-XNJ(3,3,K)
        XNK(1,1)=2.0D0*Q2(1,28)+CONST1*Q2(11,28)-XNK(1,1)
        XNK(2,1)=2.0D0*Q2(2,28)+CONST1*Q2(16,28)-XNK(2,1)
        XNK(1,2)=2.0D0*Q2(3,28)+CONST1*Q2(15,28)-XNK(1,2)
        XNK(2,2)=2.0D0*Q2(5,28)+CONST1*Q2(20,28)-XNK(2,2)
        XNK(3,1)=2.0D0*Q2(9,28)+CONST1*Q2(23,28)-XNK(3,1)
        XNK(1,3)=2.0D0*Q2(10,28)+CONST1*Q2(22,28)-XNK(1,3)
        XNK(3,2)=2.0D0*Q2(12,28)+CONST1*Q2(25,28)-XNK(3,2)
        XNK(2,3)=2.0D0*Q2(17,28)+CONST1*Q2(24,28)-XNK(2,3)
        XNK(3,3)=2.0D0*Q2(21,28)+CONST1*Q2(27,28)-XNK(3,3)
      ENDIF
      FEP(:IELEM**3)=G2*Q2(:IELEM**3,IELEM**3+1)-G1*FEP(:IELEM**3)
      FLUX0(:IELEM**3,M,I,J,K)=FEP(:IELEM**3)/DELTAE(IG)
      IF(IG.EQ.NGRP) THEN
      DCUTOFF(I,J,K)=DCUTOFF(I,J,K)+W(M)*REAL(FLUX0(1,M,I,J,K))*PL(1,M)
      ENDIF
      DO 270 L=1,NSCT
      DO 260 IEL=1,IELEM**3
      FLUX(IEL,L,I,J,K)=FLUX(IEL,L,I,J,K)+W(M)*Q2(IEL,IELEM**3+1)*
     1 PL(L,M)
  260 CONTINUE
  270 CONTINUE
*--
  280 CONTINUE
      DO 300 JEL=1,IELEM
      DO 290 IEL=1,IELEM
      IOF=((((M-1)*LY+(J-1))*LX+(I-1))*IELEM+(JEL-1))*IELEM+IEL
      FUNKNO(L4+IOFZ+IOF,IG)=REAL(XNK(IEL,JEL))
  290 CONTINUE
  300 CONTINUE
*--
  310 CONTINUE
      DO 340 K=1,LZ
      DO 330 JEL=1,IELEM
      DO 320 IEL=1,IELEM
      IOF=((((M-1)*LZ+(K-1))*LX+(I-1))*IELEM+(JEL-1))*IELEM+IEL
      FUNKNO(L4+IOFY+IOF,IG)=REAL(XNJ(IEL,JEL,K))
  320 CONTINUE
  330 CONTINUE
  340 CONTINUE
*--
  350 CONTINUE
      DO 390 K=1,LZ
      DO 380 J=1,LY
      DO 370 JEL=1,IELEM
      DO 360 IEL=1,IELEM
      IOF=((((M-1)*LZ+(K-1))*LY+(J-1))*IELEM+(JEL-1))*IELEM+IEL
      FUNKNO(L4+IOF,IG)=REAL(XNI(IEL,JEL,J,K))
  360 CONTINUE
  370 CONTINUE
  380 CONTINUE
  390 CONTINUE
      FLUX_G(:,:,:,:,:,IGEFF)=FLUX_G(:,:,:,:,:,IGEFF)+FLUX(:,:,:,:,:)
      FLUX0_G(:,:,:,:,:,IGEFF)=FLUX0_G(:,:,:,:,:,IGEFF)+FLUX0(:,:,:,:,:)
  400 CONTINUE
  410 CONTINUE
*$OMP END PARALLEL DO
  420 CONTINUE
      DO 430 IGEFF=1,NGEFF
        IF(.NOT.INCONV(IGEFF)) GO TO 430
        IG=INGIND(IGEFF)
        FUNKNO(:L4,IG)=
     1  RESHAPE(REAL(FLUX_G(:IELEM**3,:NSCT,:LX,:LY,:LZ,IGEFF)), (/L4/))
        FUNKNO(L5+1:L5+L6,IG)=
     1  RESHAPE(REAL(FLUX0_G(:IELEM**3,:NPQ,:LX,:LY,:LZ,IGEFF)),
     2  (/ L6 /) )
  430 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(FLUX0_G,FLUX_G,FLUX0,FLUX,XNJ,XNI,INDANG)
      RETURN
  500 FORMAT(16H SNFP13: thread=,I8,12H --->(group=,I4,7H angle=,I4,1H))
      END
