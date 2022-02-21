*DECK SNFT22
      SUBROUTINE SNFT22(NUN,NGEFF,IMPX,INCONV,NGIND,LX,LY,IELEM,NMAT,
     1 NPQ,NSCT,MAT,VOL,TOTAL,NCODE,ZCODE,QEXT,DU,DE,W,MRM,MRMY,DB,DA,
     2 PL,FUNKNO,ISBS,NBS,ISBSM,BS,MAXL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform one inner iteration for solving SN equations in 2D Cartesian
* for the DISCONTINUOUS GALERKIN method. Energy-angle multithreading.
* Albedo boundary conditions.
*
*Copyright:
* Copyright (C) 2019 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. A. Calloo
*
*Parameters: input
* NUN     total number of unknowns in vector FUNKNO.
* NGEFF   number of energy groups processed in parallel.
* IMPX    print flag (equal to zero for no print).
* INCONV  energy group convergence flag (set to .FALSE. if converged).
* NGIND   energy group indices assign to the NGEFF set.
* LX      number of meshes along X axis.
* LY      number of meshes along Y axis.
* IELEM   measure of order of the spatial approximation polynomial:
*         =2 linear - default for DG;
*         =3 parabolic;
*         =4 cubic - only for DG.
* NMAT    number of material mixtures.
* NPQ     number of SN directions in four octants (including zero-weight
*         directions).
* NSCT    maximum number of spherical harmonics moments of the flux.
* MAT     material mixture index in each region.
* VOL     volumes of each region.
* TOTAL   macroscopic total cross sections.
* NCODE   boundary condition indices.
* ZCODE   albedos.
* QEXT    Legendre components of the fixed source.
* DU      first direction cosines ($\\mu$).
* DE      second direction cosines ($\\eta$).
* W       weights.
* MRM     quadrature index.
* MRMY    quadrature index.
* DB      diamond-scheme parameter.
* DA      diamond-scheme parameter.
* PL      discrete values of the spherical harmonics corresponding
*         to the 2D SN quadrature.
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
      INTEGER NUN,NGEFF,IMPX,NGIND(NGEFF),LX,LY,IELEM,NMAT,NPQ,NSCT,
     1 MAT(LX,LY),NCODE(4),MRM(NPQ),MRMY(NPQ),ISBS,NBS,
     2 ISBSM(4*ISBS,NPQ*ISBS,NGEFF*ISBS),MAXL
      LOGICAL INCONV(NGEFF)
      REAL VOL(LX,LY),TOTAL(0:NMAT,NGEFF),ZCODE(4),QEXT(NUN,NGEFF),
     1 DU(NPQ),DE(NPQ),W(NPQ),DB(LX,NPQ),DA(LX,LY,NPQ),PL(NSCT,NPQ),
     2 FUNKNO(NUN,NGEFF),BS(MAXL*ISBS,NBS*ISBS)
*----
*  LOCAL VARIABLES
*----
      INTEGER NPQD(4),IIND(4)
      DOUBLE PRECISION Q(IELEM**2),Q2(IELEM**2,(IELEM**2)+1),XNJ(IELEM),
     1 VT
      PARAMETER(IUNOUT=6,PI=3.141592654)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: INDANG
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: FLUX
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: FLUX_G
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: XNI
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(INDANG(NPQ,4))
      ALLOCATE(XNI(IELEM,LY),FLUX(IELEM**2,NSCT,LX,LY))
      ALLOCATE(FLUX_G(IELEM**2,NSCT,LX,LY,NGEFF))
*----
*  DEFINITION OF CONSTANTS.
*----
      L4=IELEM*IELEM*LX*LY*NSCT
*----
*  PARAMETER VALIDATION.
*----
      IF((IELEM.LT.1).OR.(IELEM.GT.4)) CALL XABORT('SNFT22: INVALID IE'
     1 //'LEM (DIAM) VALUE. CHECK INPUT DATA FILE.')
      FLUX_G(:IELEM**2,:NSCT,:LX,:LY,:NGEFF)=0.0D0
*----
*  SET OCTANT SWAPPING ORDER.
*----
      NPQD(:4)=0
      INDANG(:NPQ,:4)=0
      DO M=1,NPQ
        VU=DU(M)
        VE=DE(M)
        IF((VU.GE.0.0).AND.(VE.GE.0.0)) THEN
          IND=1
          JND=4
        ELSE IF((VU.LE.0.0).AND.(VE.GE.0.0)) THEN
          IND=2
          JND=3
        ELSE IF((VU.LE.0.0).AND.(VE.LE.0.0)) THEN
          IND=3
          JND=1
        ELSE
          IND=4
          JND=2
        ENDIF
        IIND(JND)=IND
        NPQD(IND)=NPQD(IND)+1
        INDANG(NPQD(IND),IND)=M
      ENDDO
*----
*  MAIN LOOP OVER OCTANTS.
*----
      DO 190 JND=1,4
      IND=IIND(JND)
*----
*  PRELIMINARY LOOPS FOR SETTING BOUNDARY CONDITIONS.
*----
*$OMP  PARALLEL DO
*$OMP1 PRIVATE(M,IG,WEIGHT,VU,VE,M1,E1,IOF,JOF,IEL,I,J)
*$OMP2 SHARED(FUNKNO) COLLAPSE(2)
      DO 70 IG=1,NGEFF
      DO 60 IPQD=1,NPQD(IND)
      IF(.NOT.INCONV(IG)) GO TO 60
      M=INDANG(IPQD,IND)
      WEIGHT=W(M)
      VU=DU(M)
      VE=DE(M)
      IF(VU.GT.0.0)THEN
         M1=MRM(M)
         IF((NCODE(1).NE.4))THEN
            DO IEL=1,IELEM
               DO J=1,LY
                  IOF=((M-1)*LY+(J-1))*IELEM+IEL
                  JOF=((M1-1)*LY+(J-1))*IELEM+IEL
                  FUNKNO(L4+IOF,IG)=FUNKNO(L4+JOF,IG)
               ENDDO
            ENDDO
         ENDIF
      ELSEIF(VU.LT.0.0)THEN
         M1=MRM(M)
         IF((NCODE(2).NE.4))THEN
            DO IEL=1,IELEM
               DO J=1,LY
                  IOF=((M-1)*LY+(J-1))*IELEM+IEL
                  JOF=((M1-1)*LY+(J-1))*IELEM+IEL
                  FUNKNO(L4+IOF,IG)=FUNKNO(L4+JOF,IG)
               ENDDO
            ENDDO
         ENDIF
      ENDIF
      IF(VE.GT.0.0)THEN
         M1=MRMY(M)
         IF((NCODE(3).NE.4))THEN
            DO IEL=1,IELEM
               DO I=1,LX
                  IOF=((M-1)*LX+(I-1))*IELEM+IEL
                  JOF=((M1-1)*LX+(I-1))*IELEM+IEL
                  FUNKNO(L4+IELEM*LY*NPQ+IOF,IG)=
     >               FUNKNO(L4+IELEM*LY*NPQ+JOF,IG)
               ENDDO
            ENDDO
         ENDIF
      ELSEIF(VE.LT.0.0)THEN
         M1=MRMY(M)
         IF((NCODE(4).NE.4))THEN
            DO IEL=1,IELEM
               DO I=1,LX
                  IOF=((M-1)*LX+(I-1))*IELEM+IEL
                  JOF=((M1-1)*LX+(I-1))*IELEM+IEL
                  FUNKNO(L4+IELEM*LY*NPQ+IOF,IG)=
     >               FUNKNO(L4+IELEM*LY*NPQ+JOF,IG)
               ENDDO
            ENDDO
         ENDIF
      ENDIF
   60 CONTINUE
   70 CONTINUE
*$OMP END PARALLEL DO
*----
*  MAIN SWAPPING LOOPS FOR SN FLUX CALCULATION
*----
*$OMP  PARALLEL DO
*$OMP1 PRIVATE(ITID,FLUX,M,IG,XNI,XNJ,Q,Q2,IOF,IER,II,JJ,IEL)
*$OMP2 PRIVATE(JEL,I,J,K,VT) SHARED(FUNKNO) REDUCTION(+:FLUX_G)
*$OMP3 COLLAPSE(2)
      DO 180 IG=1,NGEFF
      DO 170 IPQD=1,NPQD(IND)
#if defined(_OPENMP)
        ITID=omp_get_thread_num()
#else
        ITID=0
#endif
      IF(IMPX.GT.5) WRITE(IUNOUT,400) ITID,NGIND(IG),IPQD
      IF(.NOT.INCONV(IG)) GO TO 170
      M=INDANG(IPQD,IND)
      FLUX(:IELEM**2,:NSCT,:LX,:LY)=0.0D0
      IF(W(M).EQ.0.0) GO TO 170
*----
*  LOOP OVER X- AND Y-DIRECTED AXES.
*----
      DO 155 II=1,LX
      I=II
      IF((IND.EQ.2).OR.(IND.EQ.3)) I=LX+1-I
      DO 100 IEL=1,IELEM
      IOF=(M-1)*IELEM*LX+(I-1)*IELEM+IEL
      IF((IND.EQ.1).OR.(IND.EQ.2)) THEN
         XNJ(IEL)=FUNKNO(L4+IELEM*LY*NPQ+IOF,IG)*ZCODE(3)
      ELSE
         XNJ(IEL)=FUNKNO(L4+IELEM*LY*NPQ+IOF,IG)*ZCODE(4)
      ENDIF
  100 CONTINUE
      IF(ISBS.EQ.1) THEN
        IF((IND.EQ.3.OR.IND.EQ.4).AND.ISBSM(4,M,IG).NE.0) THEN
          XNJ(1)=XNJ(1)+BS(I,ISBSM(4,M,IG))
        ELSE IF((IND.EQ.1.OR.IND.EQ.2).AND.ISBSM(3,M,IG).NE.0) THEN
          XNJ(1)=XNJ(1)+BS(I,ISBSM(3,M,IG))
        ENDIF
      ENDIF
      DO 140 JJ=1,LY
      J=JJ
      IF((IND.EQ.3).OR.(IND.EQ.4)) J=LY+1-J
      DO 105 IEL=1,IELEM
      IF(II.EQ.1) THEN
         IOF=(M-1)*IELEM*LY+(J-1)*IELEM+IEL
         IF((IND.EQ.1).OR.(IND.EQ.4)) THEN
            XNI(IEL,J)=FUNKNO(L4+IOF,IG)*ZCODE(1)
         ELSE
            XNI(IEL,J)=FUNKNO(L4+IOF,IG)*ZCODE(2)
         ENDIF
      ENDIF
  105 CONTINUE
      IF(ISBS.EQ.1.AND.II.EQ.1) THEN
        IF((IND.EQ.2.OR.IND.EQ.3).AND.ISBSM(2,M,IG).NE.0) THEN
          XNI(1,J)=XNI(1,J)+BS(J,ISBSM(2,M,IG))
        ELSE IF((IND.EQ.1.OR.IND.EQ.4).AND.ISBSM(1,M,IG).NE.0) THEN
          XNI(1,J)=XNI(1,J)+BS(J,ISBSM(1,M,IG))
        ENDIF
      ENDIF 
      IF(MAT(I,J).EQ.0) GO TO 140
      DO 115 IEL=1,IELEM**2
      Q(IEL)=0.0D0
      DO 110 K=1,NSCT
      IOF=((J-1)*LX*NSCT+(I-1)*NSCT+(K-1))*IELEM*IELEM+IEL
      Q(IEL)=Q(IEL)+QEXT(IOF,IG)*PL(K,M)/(4.0D0*PI)
  110 CONTINUE
  115 CONTINUE

      VT=VOL(I,J)*TOTAL(MAT(I,J),IG)
      CALL XDDSET(Q2,(IELEM**2)*((IELEM**2)+1),0.0D0)
*----------
* BUILD MATRIX FOR THE RELEVANT ORDER
* Note: the matrix does seem to be skew-symmetric (or antisymmetric), 
* i.e. A^T = -A. However, from a numerical perspective, it might be 
* more efficient to assign values than to sweep over the half matrix.
      IF(IELEM.EQ.1) THEN
      Q2(1,1) = VT + ABS(DB(I,M)) + ABS(DA(I,J,M)) 
      Q2(1,2) = Q(1)*VOL(I,J) 

      Q2(1,2) = Q2(1,2) + DB(I,M)*SIGN(1.0,DE(M))*XNJ(1) + 
     >   DA(I,J,M)*SIGN(1.0,DU(M))*XNI(1,J) 
      ELSEIF(IELEM.EQ.2) THEN
*     alpha=1 , beta=1
      Q2(1,1) = VT + ABS(DB(I,M)) + ABS(DA(I,J,M)) 
      Q2(1,2) = 3**(0.5D0)*DA(I,J,M) 
      Q2(1,3) = 3**(0.5D0)*DB(I,M) 
*     alpha=2 , beta=1
      Q2(2,1) = -3**(0.5D0)*DA(I,J,M) 
      Q2(2,2) = VT + ABS(DB(I,M)) + 3*ABS(DA(I,J,M)) 
      Q2(2,4) = 3**(0.5D0)*DB(I,M) 
*     alpha=1 , beta=2
      Q2(3,1) = -3**(0.5D0)*DB(I,M) 
      Q2(3,3) = VT + 3*ABS(DB(I,M)) + ABS(DA(I,J,M)) 
      Q2(3,4) = 3**(0.5D0)*DA(I,J,M) 
*     alpha=2 , beta=2
      Q2(4,2) = -3**(0.5D0)*DB(I,M) 
      Q2(4,3) = -3**(0.5D0)*DA(I,J,M) 
      Q2(4,4) = VT + 3*ABS(DB(I,M)) + 3*ABS(DA(I,J,M)) 
*     source
      Q2(1,5) = Q(1)*VOL(I,J)
      Q2(2,5) = Q(2)*VOL(I,J)
      Q2(3,5) = Q(3)*VOL(I,J)
      Q2(4,5) = Q(4)*VOL(I,J)
*     jump terms (dependent on direction)
      Q2(1,5) = Q2(1,5) + DB(I,M)*SIGN(1.0,DE(M))*XNJ(1) + DA(I,J
     >   ,M)*SIGN(1.0,DU(M))*XNI(1,J) 
      Q2(2,5) = Q2(2,5) + DB(I,M)*SIGN(1.0,DE(M))*XNJ(2) - 3**(0.
     >   5D0)*DA(I,J,M)*XNI(1,J) 
      Q2(3,5) = Q2(3,5) + DA(I,J,M)*SIGN(1.0,DU(M))*XNI(2,J) - 3**(
     >   0.5D0)*DB(I,M)*XNJ(1) 
      Q2(4,5) = Q2(4,5) - 3**(0.5D0)*DB(I,M)*
     >   XNJ(2) - 3**(0.5D0)*DA(I,J,M)*XNI(2,J) 
*** ---------------------------------------------------------------- ***
      ELSEIF(IELEM.EQ.3) THEN
      Q2(1,1) = VT + ABS(DB(I,M)) + ABS(DA(I,J,M)) 
      Q2(1,2) = 3**(0.5D0)*DA(I,J,M) 
      Q2(1,3) = 5**(0.5D0)*ABS(DA(I,J,M)) 
      Q2(1,4) = 3**(0.5D0)*DB(I,M) 
      Q2(1,7) = 5**(0.5D0)*ABS(DB(I,M)) 
      Q2(2,1) = -3**(0.5D0)*DA(I,J,M) 
      Q2(2,2) = VT + ABS(DB(I,M)) + 3*ABS(DA(I,J,M)) 
      Q2(2,3) = 15**(0.5D0)*DA(I,J,M) 
      Q2(2,5) = 3**(0.5D0)*DB(I,M) 
      Q2(2,8) = 5**(0.5D0)*ABS(DB(I,M)) 
      Q2(3,1) = 5**(0.5D0)*ABS(DA(I,J,M)) 
      Q2(3,2) = -15**(0.5D0)*DA(I,J,M) 
      Q2(3,3) = VT + ABS(DB(I,M)) + 5*ABS(DA(I,J,M)) 
      Q2(3,6) = 3**(0.5D0)*DB(I,M) 
      Q2(3,9) = 5**(0.5D0)*ABS(DB(I,M)) 
      Q2(4,1) = -3**(0.5D0)*DB(I,M) 
      Q2(4,4) = VT + 3*ABS(DB(I,M)) + ABS(DA(I,J,M)) 
      Q2(4,5) = 3**(0.5D0)*DA(I,J,M) 
      Q2(4,6) = 5**(0.5D0)*ABS(DA(I,J,M)) 
      Q2(4,7) = 15**(0.5D0)*DB(I,M) 
      Q2(5,2) = -3**(0.5D0)*DB(I,M) 
      Q2(5,4) = -3**(0.5D0)*DA(I,J,M) 
      Q2(5,5) = VT + 3*ABS(DB(I,M)) + 3*ABS(DA(I,J,M)) 
      Q2(5,6) = 15**(0.5D0)*DA(I,J,M) 
      Q2(5,8) = 15**(0.5D0)*DB(I,M) 
      Q2(6,3) = -3**(0.5D0)*DB(I,M) 
      Q2(6,4) = 5**(0.5D0)*ABS(DA(I,J,M)) 
      Q2(6,5) = -3**(0.5D0)*5**(0.5D0)*DA(I,J,M) 
      Q2(6,6) = VT + 3*ABS(DB(I,M)) + 5*ABS(DA(I,J,M)) 
      Q2(6,9) = 2*15**(0.5D0)*DB(I,M) - 3**(0.5D0)*5**(0.5D0)*DB(I,
     >   M) 
      Q2(7,1) = 5**(0.5D0)*ABS(DB(I,M)) 
      Q2(7,4) = -15**(0.5D0)*DB(I,M) 
      Q2(7,7) = VT + 5*ABS(DB(I,M)) + ABS(DA(I,J,M)) 
      Q2(7,8) = 3**(0.5D0)*DA(I,J,M) 
      Q2(7,9) = 5**(0.5D0)*ABS(DA(I,J,M)) 
      Q2(8,2) = 5**(0.5D0)*ABS(DB(I,M)) 
      Q2(8,5) = -3**(0.5D0)*5**(0.5D0)*DB(I,M) 
      Q2(8,7) = -3**(0.5D0)*DA(I,J,M) 
      Q2(8,8) = VT + 5*ABS(DB(I,M)) + 3*ABS(DA(I,J,M)) 
      Q2(8,9) = 2*15**(0.5D0)*DA(I,J,M) - 3**(0.5D0)*5**(0.5D0)*DA(
     >   I,J,M) 
      Q2(9,3) = 5**(0.5D0)*ABS(DB(I,M)) 
      Q2(9,6) = -15**(0.5D0)*DB(I,M) 
      Q2(9,7) = 5**(0.5D0)*ABS(DA(I,J,M)) 
      Q2(9,8) = -15**(0.5D0)*DA(I,J,M) 
      Q2(9,9) = VT + 5*ABS(DB(I,M)) + 5*ABS(DA(I,J,M)) 
      Q2(1,10) = Q(1)*VOL(I,J) 
      Q2(2,10) = Q(2)*VOL(I,J) 
      Q2(3,10) = Q(3)*VOL(I,J) 
      Q2(4,10) = Q(4)*VOL(I,J) 
      Q2(5,10) = Q(5)*VOL(I,J) 
      Q2(6,10) = Q(6)*VOL(I,J) 
      Q2(7,10) = Q(7)*VOL(I,J) 
      Q2(8,10) = Q(8)*VOL(I,J) 
      Q2(9,10) = Q(9)*VOL(I,J) 
      Q2(1,10) = Q2(1,10) + DB(I,M)*XNJ(1)*SIGN(1.0,DE(M)) + DA(I
     >   ,J,M)*XNI(1,J)*SIGN(1.0,DU(M)) 
      Q2(2,10) = Q2(2,10) + DB(I,M)*XNJ(2)*SIGN(1.0,DE(M)) - 3**(
     >   0.5D0)*DA(I,J,M)*XNI(1,J) 
      Q2(3,10) = Q2(3,10) + DB(I,M)*XNJ(3)*SIGN(1.0,DE(M)) + 5**(
     >   0.5D0)*DA(I,J,M)*XNI(1,J)*SIGN(1.0,DU(M)) 
      Q2(4,10) = Q2(4,10) + DA(I,J,M)*XNI(2,J)*SIGN(1.0,DU(M)) - 3*
     >   *(0.5D0)*DB(I,M)*XNJ(1) 
      Q2(5,10) = Q2(5,10) - 3**(0.5D0)*DB(I,M
     >   )*XNJ(2) - 3**(0.5D0)*DA(I,J,M)*XNI(2,J) 
      Q2(6,10) = Q2(6,10) + 5**(0.5D0)*DA(I,J,M
     >   )*XNI(2,J)*SIGN(1.0,DU(M)) - 3**(0.5D0)*DB(I,M)*XNJ(3) 
      Q2(7,10) = Q2(7,10) + DA(I,J,M)*XNI(3,J)*SIGN(1.0,DU(M)) + 5*
     >   *(0.5D0)*DB(I,M)*XNJ(1)*SIGN(1.0,DE(M)) 
      Q2(8,10) = Q2(8,10) + 5**(0.5D0)*DB(I,M)*SIGN(1.0,DE(M))*
     >   XNJ(2) - 3**(0.5D0)*DA(I,J,M)*XNI(3,J) 
      Q2(9,10) = Q2(9,10) + 5**(0.5D0)*DB(I,M)*SIGN(1.0,DE(M))*
     >   XNJ(3) + 5**(0.5D0)*DA(I,J,M)*XNI(3,J)*SIGN(1.0,DU(M)) 
*** ---------------------------------------------------------------- ***
      ELSEIF(IELEM.EQ.4) THEN
      Q2(1,1) = VT + ABS(DB(I,M)) + ABS(DA(I,J,M)) 
      Q2(1,2) = 3**(0.5D0)*DA(I,J,M) 
      Q2(1,3) = 5**(0.5D0)*ABS(DA(I,J,M)) 
      Q2(1,4) = 7**(0.5D0)*DA(I,J,M) 
      Q2(1,5) = 3**(0.5D0)*DB(I,M) 
      Q2(1,9) = 5**(0.5D0)*ABS(DB(I,M)) 
      Q2(1,13) = 7**(0.5D0)*DB(I,M) 

      Q2(2,1) = -3**(0.5D0)*DA(I,J,M) 
      Q2(2,2) = VT + ABS(DB(I,M)) + 3*ABS(DA(I,J,M)) 
      Q2(2,3) = 15**(0.5D0)*DA(I,J,M) 
      Q2(2,4) = 21**(0.5D0)*ABS(DA(I,J,M)) 
      Q2(2,6) = 3**(0.5D0)*DB(I,M) 
      Q2(2,10) = 5**(0.5D0)*ABS(DB(I,M)) 
      Q2(2,14) = 7**(0.5D0)*DB(I,M) 

      Q2(3,1) = 5**(0.5D0)*ABS(DA(I,J,M)) 
      Q2(3,2) = -15**(0.5D0)*DA(I,J,M) 
      Q2(3,3) = VT + ABS(DB(I,M)) + 5*ABS(DA(I,J,M)) 
      Q2(3,4) = 35**(0.5D0)*DA(I,J,M) 
      Q2(3,7) = 3**(0.5D0)*DB(I,M) 
      Q2(3,11) = 5**(0.5D0)*ABS(DB(I,M)) 
      Q2(3,15) = 7**(0.5D0)*DB(I,M) 

      Q2(4,1) = -7**(0.5D0)*DA(I,J,M) 
      Q2(4,2) = 21**(0.5D0)*ABS(DA(I,J,M)) 
      Q2(4,3) = -35**(0.5D0)*DA(I,J,M) 
      Q2(4,4) = VT + ABS(DB(I,M)) + 7*ABS(DA(I,J,M)) 
      Q2(4,8) = 3**(0.5D0)*DB(I,M) 
      Q2(4,12) = 5**(0.5D0)*ABS(DB(I,M)) 
      Q2(4,16) = 7**(0.5D0)*DB(I,M) 

      Q2(5,1) = -3**(0.5D0)*DB(I,M) 
      Q2(5,5) = VT + 3*ABS(DB(I,M)) + ABS(DA(I,J,M)) 
      Q2(5,6) = 3**(0.5D0)*DA(I,J,M) 
      Q2(5,7) = 5**(0.5D0)*ABS(DA(I,J,M)) 
      Q2(5,8) = 7**(0.5D0)*DA(I,J,M) 
      Q2(5,9) = 15**(0.5D0)*DB(I,M) 
      Q2(5,13) = 21**(0.5D0)*ABS(DB(I,M)) 

      Q2(6,2) = -3**(0.5D0)*DB(I,M) 
      Q2(6,5) = -3**(0.5D0)*DA(I,J,M) 
      Q2(6,6) = VT + 3*ABS(DB(I,M)) + 3*ABS(DA(I,J,M)) 
      Q2(6,7) = 15**(0.5D0)*DA(I,J,M) 
      Q2(6,8) = 21**(0.5D0)*ABS(DA(I,J,M)) 
      Q2(6,10) = 15**(0.5D0)*DB(I,M) 
      Q2(6,14) = 21**(0.5D0)*ABS(DB(I,M)) 

      Q2(7,3) = -3**(0.5D0)*DB(I,M) 
      Q2(7,5) = 5**(0.5D0)*ABS(DA(I,J,M)) 
      Q2(7,6) = -15**(0.5D0)*DA(I,J,M) 
      Q2(7,7) = VT + 3*ABS(DB(I,M)) + 5*ABS(DA(I,J,M)) 
      Q2(7,8) = 35**(0.5D0)*DA(I,J,M) 
      Q2(7,11) = 2*15**(0.5D0)*DB(I,M) - 3**(0.5D0)*5**(0.5D0)*DB(I
     >   ,M) 
      Q2(7,15) = 3**(0.5D0)*7**(0.5D0)*ABS(DB(I,M)) 

      Q2(8,4) = -3**(0.5D0)*DB(I,M) 
      Q2(8,5) = -7**(0.5D0)*DA(I,J,M) 
      Q2(8,6) = 21**(0.5D0)*ABS(DA(I,J,M)) 
      Q2(8,7) = -35**(0.5D0)*DA(I,J,M) 
      Q2(8,8) = VT + 3*ABS(DB(I,M)) + 7*ABS(DA(I,J,M)) 
      Q2(8,12) = 2*15**(0.5D0)*DB(I,M) - 3**(0.5D0)*5**(0.5D0)*DB(I
     >   ,M) 
      Q2(8,16) = 3**(0.5D0)*7**(0.5D0)*ABS(DB(I,M)) 

      Q2(9,1) = 5**(0.5D0)*ABS(DB(I,M)) 
      Q2(9,5) = -15**(0.5D0)*DB(I,M) 
      Q2(9,9) = VT + 5*ABS(DB(I,M)) + ABS(DA(I,J,M)) 
      Q2(9,10) = 3**(0.5D0)*DA(I,J,M) 
      Q2(9,11) = 5**(0.5D0)*ABS(DA(I,J,M)) 
      Q2(9,12) = 7**(0.5D0)*DA(I,J,M) 
      Q2(9,13) = 35**(0.5D0)*DB(I,M) 

      Q2(10,2) = 5**(0.5D0)*ABS(DB(I,M)) 
      Q2(10,6) = -15**(0.5D0)*DB(I,M) 
      Q2(10,9) = -3**(0.5D0)*DA(I,J,M) 
      Q2(10,10) = VT + 5*ABS(DB(I,M)) + 3*ABS(DA(I,J,M)) 
      Q2(10,11) = 2*15**(0.5D0)*DA(I,J,M) - 3**(0.5D0)*5**(0.5D0)*D
     >   A(I,J,M) 
      Q2(10,12) = 3**(0.5D0)*7**(0.5D0)*ABS(DA(I,J,M)) 
      Q2(10,14) = 35**(0.5D0)*DB(I,M) 

      Q2(11,3) = 5**(0.5D0)*ABS(DB(I,M)) 
      Q2(11,7) = -15**(0.5D0)*DB(I,M) 
      Q2(11,9) = 5**(0.5D0)*ABS(DA(I,J,M)) 
      Q2(11,10) = -15**(0.5D0)*DA(I,J,M) 
      Q2(11,11) = VT + 5*ABS(DB(I,M)) + 5*ABS(DA(I,J,M)) 
      Q2(11,12) = 35**(0.5D0)*DA(I,J,M) 
      Q2(11,15) = 35**(0.5D0)*DB(I,M) 

      Q2(12,4) = 5**(0.5D0)*ABS(DB(I,M)) 
      Q2(12,8) = -3**(0.5D0)*5**(0.5D0)*DB(I,M) 
      Q2(12,9) = -7**(0.5D0)*DA(I,J,M) 
      Q2(12,10) = 3**(0.5D0)*7**(0.5D0)*ABS(DA(I,J,M)) 
      Q2(12,11) = -5**(0.5D0)*7**(0.5D0)*DA(I,J,M) 
      Q2(12,12) = VT + 5*ABS(DB(I,M)) + 7*ABS(DA(I,J,M)) 
      Q2(12,16) = 2*35**(0.5D0)*DB(I,M) - 5**(0.5D0)*7**(0.5D0)*DB(
     >   I,M) 

      Q2(13,1) = -7**(0.5D0)*DB(I,M) 
      Q2(13,5) = 21**(0.5D0)*ABS(DB(I,M)) 
      Q2(13,9) = -35**(0.5D0)*DB(I,M) 
      Q2(13,13) = VT + 7*ABS(DB(I,M)) + ABS(DA(I,J,M)) 
      Q2(13,14) = 3**(0.5D0)*DA(I,J,M) 
      Q2(13,15) = 5**(0.5D0)*ABS(DA(I,J,M)) 
      Q2(13,16) = 7**(0.5D0)*DA(I,J,M) 

      Q2(14,2) = -7**(0.5D0)*DB(I,M) 
      Q2(14,6) = 21**(0.5D0)*ABS(DB(I,M)) 
      Q2(14,10) = -35**(0.5D0)*DB(I,M) 
      Q2(14,13) = -3**(0.5D0)*DA(I,J,M) 
      Q2(14,14) = VT + 7*ABS(DB(I,M)) + 3*ABS(DA(I,J,M)) 
      Q2(14,15) = 2*15**(0.5D0)*DA(I,J,M) - 3**(0.5D0)*5**(0.5D0)*D
     >   A(I,J,M) 
      Q2(14,16) = 3**(0.5D0)*7**(0.5D0)*ABS(DA(I,J,M)) 

      Q2(15,3) = -7**(0.5D0)*DB(I,M) 
      Q2(15,7) = 3**(0.5D0)*7**(0.5D0)*ABS(DB(I,M)) 
      Q2(15,11) = -5**(0.5D0)*7**(0.5D0)*DB(I,M) 
      Q2(15,13) = 5**(0.5D0)*ABS(DA(I,J,M)) 
      Q2(15,14) = -3**(0.5D0)*5**(0.5D0)*DA(I,J,M) 
      Q2(15,15) = VT + 7*ABS(DB(I,M)) + 5*ABS(DA(I,J,M)) 
      Q2(15,16) = 2*35**(0.5D0)*DA(I,J,M) - 5**(0.5D0)*7**(0.5D0)*D
     >   A(I,J,M) 

      Q2(16,4) = -7**(0.5D0)*DB(I,M) 
      Q2(16,8) = 21**(0.5D0)*ABS(DB(I,M)) 
      Q2(16,12) = -35**(0.5D0)*DB(I,M) 
      Q2(16,13) = -7**(0.5D0)*DA(I,J,M) 
      Q2(16,14) = 21**(0.5D0)*ABS(DA(I,J,M)) 
      Q2(16,15) = -35**(0.5D0)*DA(I,J,M) 
      Q2(16,16) = VT + 7*ABS(DB(I,M)) + 7*ABS(DA(I,J,M)) 

      Q2(1,17)  =  Q(1)*VOL(I,J) 
      Q2(2,17)  =  Q(2)*VOL(I,J) 
      Q2(3,17)  =  Q(3)*VOL(I,J) 
      Q2(4,17)  =  Q(4)*VOL(I,J) 
      Q2(5,17)  =  Q(5)*VOL(I,J) 
      Q2(6,17)  =  Q(6)*VOL(I,J) 
      Q2(7,17)  =  Q(7)*VOL(I,J) 
      Q2(8,17)  =  Q(8)*VOL(I,J) 
      Q2(9,17)  =  Q(9)*VOL(I,J) 
      Q2(10,17) = Q(10)*VOL(I,J) 
      Q2(11,17) = Q(11)*VOL(I,J) 
      Q2(12,17) = Q(12)*VOL(I,J) 
      Q2(13,17) = Q(13)*VOL(I,J) 
      Q2(14,17) = Q(14)*VOL(I,J) 
      Q2(15,17) = Q(15)*VOL(I,J) 
      Q2(16,17) = Q(16)*VOL(I,J) 

      Q2(1,17) = Q2(1,17) + DB(I,M)*XNJ(1)*SIGN(1.0,DE(M)) + 
     >   DA(I,J,M)*XNI(1,J)*SIGN(1.0,DU(M))  
      Q2(2,17) = Q2(2,17) + DB(I,M)*XNJ(2)*SIGN(1.0,DE(M)) - 
     >   3**(0.5D0)*DA(I,J,M)*XNI(1,J) 
      Q2(3,17) = Q2(3,17) + DB(I,M)*XNJ(3)*SIGN(1.0,DE(M)) + 
     >   5**(0.5D0)*DA(I,J,M)*XNI(1,J)*SIGN(1.0,DU(M))  
      Q2(4,17) = Q2(4,17) + DB(I,M)*XNJ(4)*SIGN(1.0,DE(M)) - 
     >   7**(0.5D0)*DA(I,J,M)*XNI(1,J) 
      Q2(5,17) = Q2(5,17) + DA(I,J,M)*XNI(2,J)*SIGN(1.0,DU(M))  - 
     >   3**(0.5D0)*DB(I,M)*XNJ(1) 
      Q2(6,17) = Q2(6,17) - 3**(0.5D0)*DB(I,M)*XNJ(2) - 3**(0.5D0)*
     >   DA(I,J,M)*XNI(2,J) 
      Q2(7,17) = Q2(7,17) + 5**(0.5D0)*DA(I,J,M)*XNI(2,J)*
     >   SIGN(1.0,DU(M)) - 3**(0.5D0)*DB(I,M)*XNJ(3) 
      Q2(8,17) = Q2(8,17) - 3**(0.5D0)*DB(I,M)*XNJ(4) - 7**(0.5D0)*
     >   DA(I,J,M)*XNI(2,J) 
      Q2(9,17) = Q2(9,17) + DA(I,J,M)*XNI(3,J)*SIGN(1.0,DU(M))  + 
     >   5**(0.5D0)*DB(I,M)*XNJ(1)*SIGN(1.0,DE(M)) 
      Q2(10,17) = Q2(10,17) + 5**(0.5D0)*DB(I,M)*XNJ(2)*SIGN(1.0,DE(M))
     >    - 3**(0.5D0)*DA(I,J,M)*XNI(3,J) 
      Q2(11,17) = Q2(11,17) + 5**(0.5D0)*DB(I,M)*XNJ(3)*SIGN(1.0,DE(M))
     >    + 5**(0.5D0)*DA(I,J,M)*XNI(3,J)*SIGN(1.0,DU(M))  
      Q2(12,17) = Q2(12,17) + 5**(0.5D0)*DB(I,M)*XNJ(4)*SIGN(1.0,DE(M))
     >    - 7**(0.5D0)*DA(I,J,M)*XNI(3,J) 
      Q2(13,17) = Q2(13,17) + DA(I,J,M)*XNI(4,J)*SIGN(1.0,DU(M)) - 
     >   7**(0.5D0)*DB(I,M)*XNJ(1) 
      Q2(14,17) = Q2(14,17) - 7**(0.5D0)*DB(I,M)*XNJ(2) - 3**(0.5D0
     >   )*DA(I,J,M)*XNI(4,J) 
      Q2(15,17) = Q2(15,17) + 5**(0.5D0)*DA(I,J,M)*XNI(4,J)*
     >   SIGN(1.0,DU(M)) - 7**(0.5D0)*DB(I,M)*XNJ(3) 
      Q2(16,17) = Q2(16,17) - 7**(0.5D0)*DB(I,M)*XNJ(4) - 7**(0.5D0
     >   )*DA(I,J,M)*XNI(4,J) 
      ENDIF
*
      CALL ALSBD(IELEM**2,1,Q2,IER,IELEM**2)
      IF(IER.NE.0) CALL XABORT('SNFT22: SINGULAR MATRIX.')
*
      IF(IELEM.EQ.1)THEN
         XNI(1,J)=Q2(1,2)
         XNJ(1)  =Q2(1,2)
      ELSEIF(IELEM.EQ.2)THEN
         XNI(1,J)=Q2(1,5) + SIGN(1.0,DU(M))*Q2(2,5)*(3**(0.5D0))
         XNI(2,J)=Q2(3,5) + SIGN(1.0,DU(M))*Q2(4,5)*(3**(0.5D0))
         XNJ(1)  =Q2(1,5) + SIGN(1.0,DE(M))*Q2(3,5)*(3**(0.5D0))
         XNJ(2)  =Q2(2,5) + SIGN(1.0,DE(M))*Q2(4,5)*(3**(0.5D0))
      ELSEIF(IELEM.EQ.3)THEN
         XNI(1,J)=Q2(1,10)+SIGN(1.0,DU(M))*(3**(0.5D0))*Q2(2,10)+
     >      (5**(0.5D0))*Q2(3,10)
         XNI(2,J)=Q2(4,10)+SIGN(1.0,DU(M))*(3**(0.5D0))*Q2(5,10)+
     >      (5**(0.5D0))*Q2(6,10)
         XNI(3,J)=Q2(7,10)+SIGN(1.0,DU(M))*(3**(0.5D0))*Q2(8,10)+
     >      (5**(0.5D0))*Q2(9,10)
         XNJ(1)  =Q2(1,10)+SIGN(1.0,DE(M))*(3**(0.5D0))*Q2(4,10)+
     >      (5**(0.5D0))*Q2(7,10)
         XNJ(2)  =Q2(2,10)+SIGN(1.0,DE(M))*(3**(0.5D0))*Q2(5,10)+
     >      (5**(0.5D0))*Q2(8,10)
         XNJ(3)  =Q2(3,10)+SIGN(1.0,DE(M))*(3**(0.5D0))*Q2(6,10)+
     >      (5**(0.5D0))*Q2(9,10)
      ELSEIF(IELEM.EQ.4)THEN
         XNI(1,J)=Q2(1,17)+SIGN(1.0,DU(M))*(3**(0.5D0))*Q2(2,17)+
     >    (5**(0.5D0))*Q2(3,17)+SIGN(1.0,DU(M))*(7**(0.5D0))*Q2(4,17)
         XNI(2,J)=Q2(5,17)+SIGN(1.0,DU(M))*(3**(0.5D0))*Q2(6,17)+
     >    (5**(0.5D0))*Q2(7,17)+SIGN(1.0,DU(M))*(7**(0.5D0))*Q2(8,17)
         XNI(3,J)=Q2(9,17)+SIGN(1.0,DU(M))*(3**(0.5D0))*Q2(10,17)+
     >    (5**(0.5D0))*Q2(11,17)+SIGN(1.0,DU(M))*(7**(0.5D0))*Q2(12,17)
         XNI(4,J)=Q2(13,17)+SIGN(1.0,DU(M))*(3**(0.5D0))*Q2(14,17)+
     >    (5**(0.5D0))*Q2(15,17)+SIGN(1.0,DU(M))*(7**(0.5D0))*Q2(16,17)
         XNJ(1)=Q2(1,17)+SIGN(1.0,DE(M))*(3**(0.5D0))*Q2(5,17)+
     >    (5**(0.5D0))*Q2(9,17)+SIGN(1.0,DE(M))*(7**(0.5D0))*Q2(13,17)
         XNJ(2)=Q2(2,17)+SIGN(1.0,DE(M))*(3**(0.5D0))*Q2(6,17)+
     >    (5**(0.5D0))*Q2(10,17)+SIGN(1.0,DE(M))*(7**(0.5D0))*Q2(14,17)
         XNJ(3)=Q2(3,17)+SIGN(1.0,DE(M))*(3**(0.5D0))*Q2(7,17)+
     >    (5**(0.5D0))*Q2(11,17)+SIGN(1.0,DE(M))*(7**(0.5D0))*Q2(15,17)
         XNJ(4)=Q2(4,17)+SIGN(1.0,DE(M))*(3**(0.5D0))*Q2(8,17)+
     >    (5**(0.5D0))*Q2(12,17)+SIGN(1.0,DE(M))*(7**(0.5D0))*Q2(16,17)
      ENDIF
*
      DO 135 K=1,NSCT
      DO 130 IEL=1,IELEM**2
      FLUX(IEL,K,I,J)=FLUX(IEL,K,I,J)+2.0D0*W(M)*Q2(IEL,IELEM**2+1)*
     1 PL(K,M)
  130 CONTINUE
  135 CONTINUE
*--
  140 CONTINUE
      DO 150 IEL=1,IELEM
      IOF=(M-1)*IELEM*LX+(I-1)*IELEM+IEL
      FUNKNO(L4+IELEM*LY*NPQ+IOF,IG)=REAL(XNJ(IEL))
  150 CONTINUE
*--
  155 CONTINUE
      DO 165 J=1,LY
      DO 160 IEL=1,IELEM
      IOF=(M-1)*IELEM*LY+(J-1)*IELEM+IEL
      FUNKNO(L4+IOF,IG)=REAL(XNI(IEL,J))
  160 CONTINUE
  165 CONTINUE
      FLUX_G(:,:,:,:,IG)=FLUX_G(:,:,:,:,IG)+FLUX(:,:,:,:)
  170 CONTINUE
  180 CONTINUE
*$OMP END PARALLEL DO
  190 CONTINUE
      DO 200 IG=1,NGEFF
        IF(.NOT.INCONV(IG)) GO TO 200
        FUNKNO(:L4,IG)=
     1  RESHAPE(REAL(FLUX_G(:IELEM**2,:NSCT,:LX,:LY,IG)), (/ L4 /) )
  200 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(XNI,FLUX_G,FLUX,INDANG)
      RETURN
  400 FORMAT(16H SNFT22: thread=,I8,12H --->(group=,I4,7H angle=,I4,1H))
      END
