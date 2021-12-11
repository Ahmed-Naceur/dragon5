*DECK SNFTH4
      SUBROUTINE SNFTH4(NUN,NGEFF,IMPX,INCONV,NGIND,NHEX,LZ,IELEM,
     1 ISPLH,SIDE,NMAT,NPQ,NSCT,MAT,VOL,NCODE,ZCODE,TOTAL,QEXT,DU,DE,
     2 DZ,W,MRMZ,DC,DB,DA,PL,CONNEC,IZGLOB,CONFROM,FUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform one inner iteration for solving SN equations in 3D hexagonal
* geometry for the DISCONTINUOUS GALERKIN method. VOID boundary
* conditions on sides, and, VOID or REFL boundary conditions on top and
* bottom. Energy-angle multithreading. Albedo boundary conditions.
*
*Copyright:
* Copyright (C) 2020 Ecole Polytechnique de Montreal
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
* NHEX    number of hexagons in X-Y plane.
* LZ      number of meshes along Z axis.
* IELEM   measure of order of the spatial approximation polynomial:
*         =1 constant - classical diamond scheme, default for HODD;
*         =2 linear - default for DG;
*         =3 parabolic;
*         =4 cubic - only for DG.
* SIDE    side of an hexagon.
* NMAT    number of material mixtures.
* ISPLH   splitting option for hexagons.
* NPQ     number of SN directions in 12 dodecants (including zero-weight
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
* DZ      third direction cosines ($\\xi$).
* W       weights.
* MRMZ    quadrature index.
* DC      diamond-scheme parameter.
* DB      diamond-scheme parameter.
* DA      diamond-scheme parameter.
* PL      discrete values of the spherical harmonics corresponding
*         to the 2D SN quadrature.
* CONNEC  connectivity matrix for flux swapping -- which lozenges is the
*         lozenge under consideration connected to; in order to pass the
*         flux along. This is dependent on direction.
* IZGLOB  hexagon sweep order depending on direction.
* CONFROM matrix for incoming flux -- which lozenges are feeding into
*         the lozenge under consideration. This is dependent on
*         direction.
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
      INTEGER NUN,NGEFF,IMPX,NGIND(NGEFF),NHEX,LZ,IELEM,ISPLH,
     1 NMAT,NPQ,NSCT,IZGLOB(NHEX,6),MAT(ISPLH,ISPLH,3,NHEX,LZ),
     2 NCODE(6),MRMZ(NPQ),CONNEC(3,(NHEX*3)*2,6),CONFROM(2,3,6)
      LOGICAL INCONV(NGEFF)
      REAL SIDE,VOL(ISPLH,ISPLH,3,NHEX,LZ),ZCODE(6),TOTAL(0:NMAT,NGEFF),
     1 QEXT(NUN,NGEFF),DU(NPQ),DE(NPQ),DZ(NPQ),W(NPQ),
     2 DC(ISPLH*ISPLH*3*NHEX,1,NPQ),DB(ISPLH*ISPLH*3*NHEX,LZ,NPQ),
     3 DA(1,LZ,NPQ),PL(NSCT,NPQ),FUNKNO(NUN,NGEFF)
*----
*  LOCAL VARIABLES
*----
      INTEGER NPQD(12),IIND(12)
      DOUBLE PRECISION Q(IELEM**3), Q2(IELEM**3,(IELEM**3)+1)
      PARAMETER(IUNOUT=6,PI=3.141592654)
      INTEGER :: ILOZSWP(3,6), IFROMI, IFROMJ
      REAL :: JAC(2,2,3), MUH, ETAH, XIH, AAA, BBB, CCC, DDD, MUHTEMP,
     1 ETAHTEMP
      DOUBLE PRECISION :: VT, THETA, XNI(IELEM,IELEM,ISPLH,LZ),
     1 XNJ(IELEM,IELEM,LZ),XNK(IELEM,IELEM)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: INDANG
      DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:,:,:,:,:) :: FLUX
      DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: FLUX_G
      DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: BFLUX
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(INDANG(NPQ,12))
      ALLOCATE(BFLUX(2,NHEX*3,ISPLH,IELEM,IELEM,LZ),
     1 FLUX(IELEM**3,NSCT,3*ISPLH**2,NHEX,LZ))
      ALLOCATE(FLUX_G(IELEM**3,NSCT,3*ISPLH**2,NHEX,LZ,NGEFF))
*----
*  CONSTRUCT JACOBIAN MATRIX FOR EACH LOZENGE
*----
      ILOZSWP = RESHAPE((/ 3, 2, 1, 3, 1, 2, 1, 3, 2, 1, 2, 3, 2, 1,
     1   3, 2, 3, 1 /), SHAPE(ILOZSWP))
      JAC = RESHAPE((/ 1., -SQRT(3.), 1., SQRT(3.), 2., 0., 1.,
     1    SQRT(3.), 2., 0., -1., SQRT(3.) /), SHAPE(JAC))
      JAC = (SIDE/2.)*JAC
*----
*  DEFINITION OF CONSTANTS.
*----
      L4=3*(IELEM**3)*(ISPLH**2)*NHEX*LZ*NSCT
*----
*  PARAMETER VALIDATION.
*----
      IF((IELEM.LT.1).OR.(IELEM.GT.4))
     1   CALL XABORT('SNFTH4: INVALID IELEM (DIAM) VALUE. '
     2   //'CHECK INPUT DATA FILE.')
      FLUX_G(:IELEM**3,:NSCT,:3*ISPLH**2,:NHEX,:LZ,:NGEFF)=0.0D0
*----
*  SET DODECANT SWAPPING ORDER.
*----
      NPQD(:12)=0
      INDANG(:NPQ,:12)=0
      DO 10 M=1,NPQ
        VU=DU(M)
        VE=DE(M)
        VZ=DZ(M)
        THETA=0.0D0
        IF(VE.GT.0.0)THEN
          IF(VU.EQ.0.0)THEN
            THETA = PI/2
          ELSEIF(VU.GT.0.0)THEN
            THETA = ATAN(ABS(VE/VU))
          ELSEIF(VU.LT.0.0)THEN
            THETA = PI - ATAN(ABS(VE/VU))
          ENDIF
        ELSEIF(VE.LT.0.0)THEN
          IF(VU.EQ.0.0)THEN
            THETA = 3*PI/2
          ELSEIF(VU.LT.0.0)THEN
            THETA = PI + ATAN(ABS(VE/VU))
          ELSEIF(VU.GT.0.0)THEN
            THETA = 2.*PI - ATAN(ABS(VE/VU))
          ENDIF
        ENDIF
*---- UNFOLD DODECANTS
        IND=0
        IF(VZ.GE.0.0)THEN
          IF((THETA.GT.0.0).AND.(THETA.LT.(PI/3.)))THEN
            IND=1
          ELSEIF((THETA.GT.(PI/3.)).AND.(THETA.LT.(2.*PI/3.)))THEN
            IND=2
          ELSEIF((THETA.GT.(2.*PI/3.)).AND.(THETA.LT.(PI)))THEN
            IND=3
          ELSEIF((THETA.GT.(PI)).AND.(THETA.LT.(4.*PI/3.)))THEN
            IND=4
          ELSEIF((THETA.GT.(4.*PI/3.)).AND.(THETA.LT.(5.*PI/3.)))THEN
            IND=5
          ELSEIF((THETA.GT.(5.*PI/3.)).AND.(THETA.LT.(2.*PI)))THEN
            IND=6
          ENDIF
        ELSEIF(VZ.LT.0.0)THEN
          IF((THETA.GT.0.0).AND.(THETA.LT.(PI/3.)))THEN
            IND=7
          ELSEIF((THETA.GT.(PI/3.)).AND.(THETA.LT.(2.*PI/3.)))THEN
            IND=8
          ELSEIF((THETA.GT.(2.*PI/3.)).AND.(THETA.LT.(PI)))THEN
            IND=9
          ELSEIF((THETA.GT.(PI)).AND.(THETA.LT.(4.*PI/3.)))THEN
            IND=10
          ELSEIF((THETA.GT.(4.*PI/3.)).AND.(THETA.LT.(5.*PI/3.)))THEN
            IND=11
          ELSEIF((THETA.GT.(5.*PI/3.)).AND.(THETA.LT.(2.*PI)))THEN
            IND=12
          ENDIF
        ENDIF
        IIND(IND)=IND ! assume IIND(I)=I in hexagonal geometry
        NPQD(IND)=NPQD(IND)+1
        INDANG(NPQD(IND),IND)=M
   10 CONTINUE
*----
*  MAIN LOOP OVER DODECANTS.
*----
      DO 310 JND=1,12
      IND=IIND(JND)
      IND2=MOD(IND-1,6)+1
      IF(IND.EQ.0) GO TO 310 ! needed because of S2 LS (8 dir. for 12 dodecants)
*----
*  PRELIMINARY LOOPS FOR SETTING BOUNDARY CONDITIONS.
*----
*$OMP  PARALLEL DO
*$OMP1 PRIVATE(M,IG,VZ,M1,E1,IOF,JOF,IEL,JEL,I,J,IL,JL)
*$OMP2 SHARED(FUNKNO) COLLAPSE(2)
      DO 150 IG=1,NGEFF
      DO 140 IPQD=1,NPQD(IND)
      IF(.NOT.INCONV(IG)) GO TO 140
      M=INDANG(IPQD,IND)
      VZ=DZ(M)
*
      IF(VZ.GT.0.0)THEN
         M1=MRMZ(M)
         IF(NCODE(5).NE.4)THEN
         DO IEL=1,IELEM
         DO JEL=1,IELEM
           DO I=1,NHEX
           DO J=1,3
           DO IL=1,ISPLH
           DO JL=1,ISPLH
           IOF=((((((M-1)*NHEX+(I-1))*3+(J-1))*ISPLH+(JL-1))*ISPLH+
     >      (IL-1))*IELEM+(JEL-1))*IELEM+IEL
           JOF=((((((M1-1)*NHEX+(I-1))*3+(J-1))*ISPLH+(JL-1))*ISPLH+
     >      (IL-1))*IELEM+(JEL-1))*IELEM+IEL
           FUNKNO(L4+IOF,IG)=FUNKNO(L4+JOF,IG)
           ENDDO
           ENDDO
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
           DO I=1,NHEX
           DO J=1,3
           DO IL=1,ISPLH
           DO JL=1,ISPLH
           IOF=((((((M-1)*NHEX+(I-1))*3+(J-1))*ISPLH+(JL-1))*ISPLH+
     >      (IL-1))*IELEM+(JEL-1))*IELEM+IEL
           JOF=((((((M1-1)*NHEX+(I-1))*3+(J-1))*ISPLH+(JL-1))*ISPLH+
     >      (IL-1))*IELEM+(JEL-1))*IELEM+IEL
           FUNKNO(L4+IOF,IG)=FUNKNO(L4+JOF,IG)
           ENDDO
           ENDDO
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
*$OMP1 PRIVATE(ITID,FLUX,BFLUX,M,IG,XNI,XNJ,XNK,Q,Q2,IOF,IER,II,JJ,IZ)
*$OMP2 PRIVATE(IEL,JEL,I,J,K,IL,JL,VT,AAA,BBB,CCC,DDD,IHEXI,ILOZI,ILOZ)
*$OMP3 PRIVATE(ISIDEI,IHEXJ,ILOZJ,ISIDEJ,IFROMI,IFROMJ,INDEXI,INDEXJ)
*$OMP4 PRIVATE(MUHTEMP,ETAHTEMP,MUH,ETAH,XIH,I2,I3,J2,J3,K2,C1,CORNERQ)
*$OMP5 SHARED(FUNKNO) REDUCTION(+:FLUX_G) COLLAPSE(2)
      DO 300 IG=1,NGEFF
      DO 290 IPQD=1,NPQD(IND)
#if defined(_OPENMP)
        ITID=omp_get_thread_num()
#else
        ITID=0
#endif
      IF(IMPX.GT.5) WRITE(IUNOUT,500) ITID,NGIND(IG),IPQD
      IF(.NOT.INCONV(IG)) GO TO 290
      M=INDANG(IPQD,IND)
      FLUX(:IELEM**3,:NSCT,:3*ISPLH**2,:NHEX,:LZ)=0.0D0
      BFLUX(:2,:NHEX*3,:ISPLH,:IELEM,:IELEM,:LZ)=0.0D0
      IF(W(M).EQ.0.0) GO TO 290
*----
*  LOOP OVER X-Y PLANE AND Z-DIRECTED AXIS.
*----
* SECOND loop over hexagons
      DO 280 II=1,NHEX
      I=IZGLOB(II,IND2)

* THIRD loop over lozenges
      DO 270 JJ=1,3
      J=ILOZSWP(JJ,IND2)
*
      AAA = JAC(1,1,J)
      BBB = JAC(1,2,J)
      CCC = JAC(2,1,J)
      DDD = JAC(2,2,J)
*
      IHEXI  = CONNEC(1,((I-1)*3*2) + ((J-1)*2) +1,IND2)
      ILOZI  = CONNEC(2,((I-1)*3*2) + ((J-1)*2) +1,IND2)
      ISIDEI = CONNEC(3,((I-1)*3*2) + ((J-1)*2) +1,IND2)
      IHEXJ  = CONNEC(1,((I-1)*3*2) + ((J-1)*2) +2,IND2)
      ILOZJ  = CONNEC(2,((I-1)*3*2) + ((J-1)*2) +2,IND2)
      ISIDEJ = CONNEC(3,((I-1)*3*2) + ((J-1)*2) +2,IND2)
      IFROMI = CONFROM(1,J,IND2)
      IFROMJ = CONFROM(2,J,IND2)
      INDEXI = ((IHEXI-1)*3)+ILOZI
      INDEXJ = ((IHEXJ-1)*3)+ILOZJ
*
      DO 260 IL=1,ISPLH
      I2=IL
      IF(CONFROM(1,J,IND2).EQ.3) I2=ISPLH+1-I2
      DO 250 JL=1,ISPLH
      J2=JL
      IF(CONFROM(2,J,IND2).EQ.4) J2=ISPLH+1-J2
      DO IEL=1,IELEM
      DO JEL=1,IELEM
        IOF=((((((M-1)*NHEX+(I-1))*3+(J-1))*ISPLH+(J2-1))*ISPLH+
     1      (I2-1))*IELEM+(JEL-1))*IELEM+IEL
        IF(IND.LT.7) THEN
          XNK(IEL,JEL)=FUNKNO(L4+IOF,IG)*ZCODE(5)
        ELSE
          XNK(IEL,JEL)=FUNKNO(L4+IOF,IG)*ZCODE(6)
        ENDIF
      ENDDO
      ENDDO
      DO 220 IZ=1,LZ
      K2= IZ
      IF(IND.GE.7) K2=LZ+1-IZ
      IF(JL.EQ.1)THEN
         DO IEL=1,IELEM
            DO JEL=1,IELEM
               XNJ(IEL,JEL,K2) = BFLUX(2,((I-1)*3)+J,I2,IEL,JEL,K2)
            ENDDO
         ENDDO
      ENDIF
      IF(IL.EQ.1)THEN
         DO IEL=1,IELEM
            DO JEL=1,IELEM
               XNI(IEL,JEL,J2,K2) = BFLUX(1,((I-1)*3)+J,J2,IEL,JEL,K2)
            ENDDO
         ENDDO
      ENDIF
*
      MUHTEMP  =  DA(1,K2,M)
      ETAHTEMP =  DB(1,K2,M)
      MUH  = (MUHTEMP*DDD) - (ETAHTEMP*BBB)
      ETAH = (-MUHTEMP*CCC) + (ETAHTEMP*AAA)
      XIH  = DC(1,1,M)
*
      IF(MAT(I2,J2,J,I,K2).EQ.0) GO TO 220

      DO 170 IEL=1,IELEM**3
      Q(IEL)=0.0D0
      DO 160 K=1,NSCT
      IOF=(((((((K2-1)*NHEX+(I-1))*3+(J-1))*ISPLH+(J2-1))*ISPLH+
     1 (I2-1))*NSCT+(K-1))*IELEM**3)+IEL
      Q(IEL)=Q(IEL)+QEXT(IOF,IG)*PL(K,M)/(4.0D0*PI)
  160 CONTINUE
  170 CONTINUE
*
      VT=VOL(I2,J2,J,I,K2)*TOTAL(MAT(I2,J2,J,I,K2),IG)
      CALL XDDSET(Q2,(IELEM**3)*((IELEM**3)+1),0.0D0)
*
*----------
* BUILD MATRIX FOR THE RELEVANT ORDER
* Note: the matrix does seem to be skew-symmetric (or antisymmetric), 
* i.e. A^T = -A. However, from a numerical perspective, it might be 
* more efficient to assign values than to sweep over the half matrix.
      IF(IELEM.EQ.1) THEN  
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
      Q2(1,1) = VT + ABS(ETAH) + ABS(MUH) + ABS(XIH) 
      Q2(1,2) = Q(1)*VOL(I2,J2,J,I,K2) 
      Q2(1,2) = Q2(1,2) + XNK(1,1)*XIH*SIGN(1.0,XIH) + 
     >   XNJ(1,1,K2)*ETAH*SIGN(1.0,ETAH) + 
     >   MUH*XNI(1,1,J2,K2)*SIGN(1.0,MUH)
      ELSEIF(IELEM.EQ.2) THEN  
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
      Q2(1,1) = VT + ABS(ETAH) + ABS(MUH) + ABS(XIH) 
      Q2(1,2) = 3**(0.5D0)*MUH 
      Q2(1,3) = 3**(0.5D0)*ETAH 
      Q2(1,5) = 3**(0.5D0)*XIH 
      Q2(2,1) = -3**(0.5D0)*MUH 
      Q2(2,2) = VT + ABS(ETAH) + 3*ABS(MUH) + ABS(XIH) 
      Q2(2,4) = 3**(0.5D0)*ETAH 
      Q2(2,6) = 3**(0.5D0)*XIH 
      Q2(3,1) = -3**(0.5D0)*ETAH 
      Q2(3,3) = VT + 3*ABS(ETAH) + ABS(MUH) + ABS(XIH) 
      Q2(3,4) = 3**(0.5D0)*MUH 
      Q2(3,7) = 3**(0.5D0)*XIH 
      Q2(4,2) = -3**(0.5D0)*ETAH 
      Q2(4,3) = -3**(0.5D0)*MUH 
      Q2(4,4) = VT + 3*ABS(ETAH) + 3*ABS(MUH) + ABS(XIH) 
      Q2(4,8) = 3**(0.5D0)*XIH 
      Q2(5,1) = -3**(0.5D0)*XIH 
      Q2(5,5) = VT + ABS(ETAH) + ABS(MUH) + 3*ABS(XIH) 
      Q2(5,6) = 3**(0.5D0)*MUH 
      Q2(5,7) = 3**(0.5D0)*ETAH 
      Q2(6,2) = -3**(0.5D0)*XIH 
      Q2(6,5) = -3**(0.5D0)*MUH 
      Q2(6,6) = VT + ABS(ETAH) + 3*ABS(MUH) + 3*ABS(XIH) 
      Q2(6,8) = 3**(0.5D0)*ETAH 
      Q2(7,3) = -3**(0.5D0)*XIH 
      Q2(7,5) = -3**(0.5D0)*ETAH 
      Q2(7,7) = VT + 3*ABS(ETAH) + ABS(MUH) + 3*ABS(XIH) 
      Q2(7,8) = 3**(0.5D0)*MUH 
      Q2(8,4) = -3**(0.5D0)*XIH 
      Q2(8,6) = -3**(0.5D0)*ETAH 
      Q2(8,7) = -3**(0.5D0)*MUH 
      Q2(8,8) = VT + 3*ABS(ETAH) + 3*ABS(MUH) + 
     >   3*ABS(XIH) 

      Q2(1,9) = Q(1)*VOL(I2,J2,J,I,K2) 
      Q2(2,9) = Q(2)*VOL(I2,J2,J,I,K2) 
      Q2(3,9) = Q(3)*VOL(I2,J2,J,I,K2) 
      Q2(4,9) = Q(4)*VOL(I2,J2,J,I,K2) 
      Q2(5,9) = Q(5)*VOL(I2,J2,J,I,K2) 
      Q2(6,9) = Q(6)*VOL(I2,J2,J,I,K2) 
      Q2(7,9) = Q(7)*VOL(I2,J2,J,I,K2) 
      Q2(8,9) = Q(8)*VOL(I2,J2,J,I,K2) 

      Q2(1,9) = Q2(1,9) + XNK(1,1)*XIH*SIGN(1.0,XIH) + 
     >   ETAH*XNJ(1,1,K2)*SIGN(1.0,ETAH) + MUH*XNI(1,1,J2,K2)*
     >   SIGN(1.0,MUH) 
      Q2(2,9) = Q2(2,9) + XNK(2,1)*XIH*SIGN(1.0,XIH) + 
     >   ETAH*XNJ(2,1,K2)*SIGN(1.0,ETAH) - 3**(0.5D0)*MUH*
     >   XNI(1,1,J2,K2) 
      Q2(3,9) = Q2(3,9) + XNK(1,2)*XIH*SIGN(1.0,XIH) + 
     >   MUH*XNI(2,1,J2,K2)*SIGN(1.0,MUH) - 3**(0.5D0)*ETAH*
     >   XNJ(1,1,K2) 
      Q2(4,9) = Q2(4,9) + XNK(2,2)*XIH*SIGN(1.0,XIH) - 3**(
     >   0.5D0)*ETAH*XNJ(2,1,K2) - 3**(0.5D0)*MUH*
     >   XNI(2,1,J2,K2) 
      Q2(5,9) = Q2(5,9) + ETAH*XNJ(1,2,K2)*SIGN(1.0,ETAH) + 
     >   MUH*XNI(1,2,J2,K2)*SIGN(1.0,MUH) - 3**(0.5D0)*XNK(1,1)*
     >   XIH 
      Q2(6,9) = Q2(6,9) + ETAH*XNJ(2,2,K2)*SIGN(1.0,ETAH) - 3*
     >   *(0.5D0)*MUH*XNI(1,2,J2,K2) - 3**(0.5D0)*XNK(2,1)*XIH 
      Q2(7,9) = Q2(7,9) + MUH*XNI(2,2,J2,K2)*SIGN(1.0,MUH) - 
     >   3**(0.5D0)*ETAH*XNJ(1,2,K2) - 3**(0.5D0)*XNK(1,2)*XIH 
      Q2(8,9) = Q2(8,9) - 3**(0.5D0)*ETAH*XNJ(2,2,K2) - 3**(0
     >   .5D0)*MUH*XNI(2,2,J2,K2) - 3**(0.5D0)*XNK(2,2)*XIH 
      ELSE IF(IELEM.EQ.3) THEN
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
      Q2(1,1) = VT + ABS(ETAH) + ABS(MUH) + ABS(XIH) 
      Q2(1,2) = 3**(0.5D0)*MUH 
      Q2(1,3) = 5**(0.5D0)*ABS(MUH) 
      Q2(1,4) = 3**(0.5D0)*ETAH 
      Q2(1,7) = 5**(0.5D0)*ABS(ETAH) 
      Q2(1,10) = 3**(0.5D0)*XIH 
      Q2(1,19) = 5**(0.5D0)*ABS(XIH) 
      Q2(2,1) = -3**(0.5D0)*MUH 
      Q2(2,2) = VT + ABS(ETAH) + 3*ABS(MUH) + ABS(XIH) 
      Q2(2,3) = 15**(0.5D0)*MUH 
      Q2(2,5) = 3**(0.5D0)*ETAH 
      Q2(2,8) = 5**(0.5D0)*ABS(ETAH) 
      Q2(2,11) = 3**(0.5D0)*XIH 
      Q2(2,20) = 5**(0.5D0)*ABS(XIH) 
      Q2(3,1) = 5**(0.5D0)*ABS(MUH) 
      Q2(3,2) = -15**(0.5D0)*MUH 
      Q2(3,3) = VT + ABS(ETAH) + 5*ABS(MUH) + ABS(XIH) 
      Q2(3,6) = 3**(0.5D0)*ETAH 
      Q2(3,9) = 5**(0.5D0)*ABS(ETAH) 
      Q2(3,12) = 3**(0.5D0)*XIH 
      Q2(3,21) = 5**(0.5D0)*ABS(XIH) 
      Q2(4,1) = -3**(0.5D0)*ETAH 
      Q2(4,4) = VT + 3*ABS(ETAH) + ABS(MUH) + ABS(XIH) 
      Q2(4,5) = 3**(0.5D0)*MUH 
      Q2(4,6) = 5**(0.5D0)*ABS(MUH) 
      Q2(4,7) = 15**(0.5D0)*ETAH 
      Q2(4,13) = 3**(0.5D0)*XIH 
      Q2(4,22) = 5**(0.5D0)*ABS(XIH) 
      Q2(5,2) = -3**(0.5D0)*ETAH 
      Q2(5,4) = -3**(0.5D0)*MUH 
      Q2(5,5) = VT + 3*ABS(ETAH) + 3*ABS(MUH) + ABS(XIH) 
      Q2(5,6) = 15**(0.5D0)*MUH 
      Q2(5,8) = 15**(0.5D0)*ETAH 
      Q2(5,14) = 3**(0.5D0)*XIH 
      Q2(5,23) = 5**(0.5D0)*ABS(XIH) 
      Q2(6,3) = -3**(0.5D0)*ETAH 
      Q2(6,4) = 5**(0.5D0)*ABS(MUH) 
      Q2(6,5) = -15**(0.5D0)*MUH 
      Q2(6,6) = VT + 3*ABS(ETAH) + 5*ABS(MUH) + ABS(XIH) 
      Q2(6,9) = 15**(0.5D0)*ETAH 
      Q2(6,15) = 3**(0.5D0)*XIH 
      Q2(6,24) = 5**(0.5D0)*ABS(XIH) 
      Q2(7,1) = 5**(0.5D0)*ABS(ETAH) 
      Q2(7,4) = -15**(0.5D0)*ETAH 
      Q2(7,7) = VT + 5*ABS(ETAH) + ABS(MUH) + ABS(XIH) 
      Q2(7,8) = 3**(0.5D0)*MUH 
      Q2(7,9) = 5**(0.5D0)*ABS(MUH) 
      Q2(7,16) = 3**(0.5D0)*XIH 
      Q2(7,25) = 5**(0.5D0)*ABS(XIH) 
      Q2(8,2) = 5**(0.5D0)*ABS(ETAH) 
      Q2(8,5) = -15**(0.5D0)*ETAH 
      Q2(8,7) = -3**(0.5D0)*MUH 
      Q2(8,8) = VT + 5*ABS(ETAH) + 3*ABS(MUH) + (3**(0.5
     >   D0)*5**(0.5D0)*15**(0.5D0)*ABS(XIH))/15 
      Q2(8,9) = 15**(0.5D0)*MUH 
      Q2(8,17) = 3**(0.5D0)*XIH 
      Q2(8,26) = 5**(0.5D0)*ABS(XIH) 
      Q2(9,3) = 5**(0.5D0)*ABS(ETAH) 
      Q2(9,6) = -15**(0.5D0)*ETAH 
      Q2(9,7) = 5**(0.5D0)*ABS(MUH) 
      Q2(9,8) = -15**(0.5D0)*MUH 
      Q2(9,9) = VT + 5*ABS(ETAH) + 5*ABS(MUH) + ABS(XIH) 
      Q2(9,18) = 3**(0.5D0)*XIH 
      Q2(9,27) = 5**(0.5D0)*ABS(XIH) 
      Q2(10,1) = -3**(0.5D0)*XIH 
      Q2(10,10) = VT + ABS(ETAH) + ABS(MUH) + 3*ABS(XIH) 
      Q2(10,11) = 3**(0.5D0)*MUH 
      Q2(10,12) = 5**(0.5D0)*ABS(MUH) 
      Q2(10,13) = 3**(0.5D0)*ETAH 
      Q2(10,16) = 5**(0.5D0)*ABS(ETAH) 
      Q2(10,19) = 15**(0.5D0)*XIH 
      Q2(11,2) = -3**(0.5D0)*XIH 
      Q2(11,10) = -3**(0.5D0)*MUH 
      Q2(11,11) = VT + ABS(ETAH) + 3*ABS(MUH) + 3*ABS(
     >   XIH) 
      Q2(11,12) = 2*15**(0.5D0)*MUH - 3**(0.5D0)*5**(0.5D0)*
     >   MUH 
      Q2(11,14) = 3**(0.5D0)*ETAH 
      Q2(11,17) = 5**(0.5D0)*ABS(ETAH) 
      Q2(11,20) = 15**(0.5D0)*XIH 
      Q2(12,3) = -3**(0.5D0)*XIH 
      Q2(12,10) = 5**(0.5D0)*ABS(MUH) 
      Q2(12,11) = -3**(0.5D0)*5**(0.5D0)*MUH 
      Q2(12,12) = VT + ABS(ETAH) + 5*ABS(MUH) + 3*ABS(
     >   XIH) 
      Q2(12,15) = 3**(0.5D0)*ETAH 
      Q2(12,18) = 5**(0.5D0)*ABS(ETAH) 
      Q2(12,21) = 15**(0.5D0)*XIH 
      Q2(13,4) = -3**(0.5D0)*XIH 
      Q2(13,10) = -3**(0.5D0)*ETAH 
      Q2(13,13) = VT + 3*ABS(ETAH) + ABS(MUH) + 3*ABS(
     >   XIH) 
      Q2(13,14) = 3**(0.5D0)*MUH 
      Q2(13,15) = 5**(0.5D0)*ABS(MUH) 
      Q2(13,16) = 2*15**(0.5D0)*ETAH - 3**(0.5D0)*5**(0.5D0)*
     >   ETAH 
      Q2(13,22) = 2*15**(0.5D0)*XIH - 3**(0.5D0)*5**(0.5D0)*
     >   XIH 
      Q2(14,5) = -3**(0.5D0)*XIH 
      Q2(14,11) = -3**(0.5D0)*ETAH 
      Q2(14,13) = -3**(0.5D0)*MUH 
      Q2(14,14) = VT + 3*ABS(ETAH) + 3*ABS(MUH) + 3*ABS(
     >   XIH) 
      Q2(14,15) = 2*15**(0.5D0)*MUH - 3**(0.5D0)*5**(0.5D0)*
     >   MUH 
      Q2(14,17) = 2*15**(0.5D0)*ETAH - 3**(0.5D0)*5**(0.5D0)*
     >   ETAH 
      Q2(14,23) = 2*15**(0.5D0)*XIH - 3**(0.5D0)*5**(0.5D0)*
     >   XIH 
      Q2(15,6) = -3**(0.5D0)*XIH 
      Q2(15,12) = -3**(0.5D0)*ETAH 
      Q2(15,13) = 5**(0.5D0)*ABS(MUH) 
      Q2(15,14) = -3**(0.5D0)*5**(0.5D0)*MUH 
      Q2(15,15) = VT + 3*ABS(ETAH) + 5*ABS(MUH) + 3*ABS(
     >   XIH) 
      Q2(15,18) = 2*15**(0.5D0)*ETAH - 3**(0.5D0)*5**(0.5D0)*
     >   ETAH 
      Q2(15,24) = 2*15**(0.5D0)*XIH - 3**(0.5D0)*5**(0.5D0)*
     >   XIH 
      Q2(16,7) = -3**(0.5D0)*XIH 
      Q2(16,10) = 5**(0.5D0)*ABS(ETAH) 
      Q2(16,13) = -3**(0.5D0)*5**(0.5D0)*ETAH 
      Q2(16,16) = VT + 5*ABS(ETAH) + ABS(MUH) + 3*ABS(
     >   XIH) 
      Q2(16,17) = 3**(0.5D0)*MUH 
      Q2(16,18) = 5**(0.5D0)*ABS(MUH) 
      Q2(16,25) = 2*15**(0.5D0)*XIH - 3**(0.5D0)*5**(0.5D0)*
     >   XIH 
      Q2(17,8) = -3**(0.5D0)*XIH 
      Q2(17,11) = 5**(0.5D0)*ABS(ETAH) 
      Q2(17,14) = -3**(0.5D0)*5**(0.5D0)*ETAH 
      Q2(17,16) = -3**(0.5D0)*MUH 
      Q2(17,17) = VT + 5*ABS(ETAH) + 3*ABS(MUH) + 3*ABS(
     >   XIH) 
      Q2(17,18) = 2*15**(0.5D0)*MUH - 3**(0.5D0)*5**(0.5D0)*
     >   MUH 
      Q2(17,26) = 2*15**(0.5D0)*XIH - 3**(0.5D0)*5**(0.5D0)*
     >   XIH 
      Q2(18,9) = -3**(0.5D0)*XIH 
      Q2(18,12) = 5**(0.5D0)*ABS(ETAH) 
      Q2(18,15) = -3**(0.5D0)*5**(0.5D0)*ETAH 
      Q2(18,16) = 5**(0.5D0)*ABS(MUH) 
      Q2(18,17) = -3**(0.5D0)*5**(0.5D0)*MUH 
      Q2(18,18) = VT + 5*ABS(ETAH) + 5*ABS(MUH) + 3*ABS(
     >   XIH) 
      Q2(18,27) = 2*15**(0.5D0)*XIH - 3**(0.5D0)*5**(0.5D0)*
     >   XIH 
      Q2(19,1) = 5**(0.5D0)*ABS(XIH) 
      Q2(19,10) = -15**(0.5D0)*XIH 
      Q2(19,19) = VT + ABS(ETAH) + ABS(MUH) + 5*ABS(XIH) 
      Q2(19,20) = 3**(0.5D0)*MUH 
      Q2(19,21) = 5**(0.5D0)*ABS(MUH) 
      Q2(19,22) = 3**(0.5D0)*ETAH 
      Q2(19,25) = 5**(0.5D0)*ABS(ETAH) 
      Q2(20,2) = 5**(0.5D0)*ABS(XIH) 
      Q2(20,11) = -15**(0.5D0)*XIH 
      Q2(20,19) = -3**(0.5D0)*MUH 
      Q2(20,20) = VT + ABS(ETAH) + 3*ABS(MUH) + 5*ABS(
     >   XIH) 
      Q2(20,21) = 2*15**(0.5D0)*MUH - 3**(0.5D0)*5**(0.5D0)*
     >   MUH 
      Q2(20,23) = 3**(0.5D0)*ETAH 
      Q2(20,26) = 5**(0.5D0)*ABS(ETAH) 
      Q2(21,3) = 5**(0.5D0)*ABS(XIH) 
      Q2(21,12) = -15**(0.5D0)*XIH 
      Q2(21,19) = 5**(0.5D0)*ABS(MUH) 
      Q2(21,20) = -3**(0.5D0)*5**(0.5D0)*MUH 
      Q2(21,21) = VT + ABS(ETAH) + 5*ABS(MUH) + 5*ABS(
     >   XIH) 
      Q2(21,24) = 3**(0.5D0)*ETAH 
      Q2(21,27) = 5**(0.5D0)*ABS(ETAH) 
      Q2(22,4) = 5**(0.5D0)*ABS(XIH) 
      Q2(22,13) = -3**(0.5D0)*5**(0.5D0)*XIH 
      Q2(22,19) = -3**(0.5D0)*ETAH 
      Q2(22,22) = VT + 3*ABS(ETAH) + ABS(MUH) + 5*ABS(
     >   XIH) 
      Q2(22,23) = 3**(0.5D0)*MUH 
      Q2(22,24) = 5**(0.5D0)*ABS(MUH) 
      Q2(22,25) = 2*15**(0.5D0)*ETAH - 3**(0.5D0)*5**(0.5D0)*
     >   ETAH 
      Q2(23,5) = (3**(0.5D0)*15**(0.5D0)*ABS(XIH))/3 
      Q2(23,14) = -15**(0.5D0)*XIH 
      Q2(23,20) = -3**(0.5D0)*ETAH 
      Q2(23,22) = -3**(0.5D0)*MUH 
      Q2(23,23) = VT + 3*ABS(ETAH) + 3*ABS(MUH) + (3**(0
     >   .5D0)*5**(0.5D0)*15**(0.5D0)*ABS(XIH))/3 
      Q2(23,24) = 2*15**(0.5D0)*MUH - 3**(0.5D0)*5**(0.5D0)*
     >   MUH 
      Q2(23,26) = 2*15**(0.5D0)*ETAH - 3**(0.5D0)*5**(0.5D0)*
     >   ETAH 
      Q2(24,6) = 5**(0.5D0)*ABS(XIH) 
      Q2(24,15) = -3**(0.5D0)*5**(0.5D0)*XIH 
      Q2(24,21) = -3**(0.5D0)*ETAH 
      Q2(24,22) = 5**(0.5D0)*ABS(MUH) 
      Q2(24,23) = -3**(0.5D0)*5**(0.5D0)*MUH 
      Q2(24,24) = VT + 3*ABS(ETAH) + 5*ABS(MUH) + 5*ABS(
     >   XIH) 
      Q2(24,27) = 2*15**(0.5D0)*ETAH - 3**(0.5D0)*5**(0.5D0)*
     >   ETAH 
      Q2(25,7) = 5**(0.5D0)*ABS(XIH) 
      Q2(25,16) = -3**(0.5D0)*5**(0.5D0)*XIH 
      Q2(25,19) = 5**(0.5D0)*ABS(ETAH) 
      Q2(25,22) = -3**(0.5D0)*5**(0.5D0)*ETAH 
      Q2(25,25) = VT + 5*ABS(ETAH) + ABS(MUH) + 5*ABS(
     >   XIH) 
      Q2(25,26) = 3**(0.5D0)*MUH 
      Q2(25,27) = 5**(0.5D0)*ABS(MUH) 
      Q2(26,8) = 5**(0.5D0)*ABS(XIH) 
      Q2(26,17) = -15**(0.5D0)*XIH 
      Q2(26,20) = 5**(0.5D0)*ABS(ETAH) 
      Q2(26,23) = -3**(0.5D0)*5**(0.5D0)*ETAH 
      Q2(26,25) = -3**(0.5D0)*MUH 
      Q2(26,26) = VT + 5*ABS(ETAH) + 3*ABS(MUH) + (3**(0
     >   .5D0)*5**(0.5D0)*15**(0.5D0)*ABS(XIH))/3 
      Q2(26,27) = 2*15**(0.5D0)*MUH - 3**(0.5D0)*5**(0.5D0)*
     >   MUH 
      Q2(27,9) = 5**(0.5D0)*ABS(XIH) 
      Q2(27,18) = -3**(0.5D0)*5**(0.5D0)*XIH 
      Q2(27,21) = 5**(0.5D0)*ABS(ETAH) 
      Q2(27,24) = -3**(0.5D0)*5**(0.5D0)*ETAH 
      Q2(27,25) = 5**(0.5D0)*ABS(MUH) 
      Q2(27,26) = -3**(0.5D0)*5**(0.5D0)*MUH 
      Q2(27,27) = VT + 5*ABS(ETAH) + 5*ABS(MUH) + 5*ABS(
     >   XIH) 

      Q2(01,28) = Q(01)*VOL(I2,J2,J,I,K2) 
      Q2(02,28) = Q(02)*VOL(I2,J2,J,I,K2) 
      Q2(03,28) = Q(03)*VOL(I2,J2,J,I,K2) 
      Q2(04,28) = Q(04)*VOL(I2,J2,J,I,K2) 
      Q2(05,28) = Q(05)*VOL(I2,J2,J,I,K2) 
      Q2(06,28) = Q(06)*VOL(I2,J2,J,I,K2) 
      Q2(07,28) = Q(07)*VOL(I2,J2,J,I,K2) 
      Q2(08,28) = Q(08)*VOL(I2,J2,J,I,K2) 
      Q2(09,28) = Q(09)*VOL(I2,J2,J,I,K2) 
      Q2(10,28) = Q(10)*VOL(I2,J2,J,I,K2) 
      Q2(11,28) = Q(11)*VOL(I2,J2,J,I,K2) 
      Q2(12,28) = Q(12)*VOL(I2,J2,J,I,K2) 
      Q2(13,28) = Q(13)*VOL(I2,J2,J,I,K2) 
      Q2(14,28) = Q(14)*VOL(I2,J2,J,I,K2) 
      Q2(15,28) = Q(15)*VOL(I2,J2,J,I,K2) 
      Q2(16,28) = Q(16)*VOL(I2,J2,J,I,K2) 
      Q2(17,28) = Q(17)*VOL(I2,J2,J,I,K2) 
      Q2(18,28) = Q(18)*VOL(I2,J2,J,I,K2) 
      Q2(19,28) = Q(19)*VOL(I2,J2,J,I,K2) 
      Q2(20,28) = Q(20)*VOL(I2,J2,J,I,K2) 
      Q2(21,28) = Q(21)*VOL(I2,J2,J,I,K2) 
      Q2(22,28) = Q(22)*VOL(I2,J2,J,I,K2) 
      Q2(23,28) = Q(23)*VOL(I2,J2,J,I,K2) 
      Q2(24,28) = Q(24)*VOL(I2,J2,J,I,K2) 
      Q2(25,28) = Q(25)*VOL(I2,J2,J,I,K2) 
      Q2(26,28) = Q(26)*VOL(I2,J2,J,I,K2) 
      Q2(27,28) = Q(27)*VOL(I2,J2,J,I,K2) 

      Q2(1,28) = Q2(1,28) + XNK(1,1)*XIH*SIGN(1.0,XIH) + 
     >   ETAH*XNJ(1,1,K2)*SIGN(1.0,ETAH) + MUH*XNI(1,1,J2,K2)*
     >   SIGN(1.0,MUH) 
      Q2(2,28) = Q2(2,28) + XNK(2,1)*XIH*SIGN(1.0,XIH) + 
     >   ETAH*XNJ(2,1,K2)*SIGN(1.0,ETAH) - 3**(0.5D0)*MUH*
     >   XNI(1,1,J2,K2) 
      Q2(3,28) = Q2(3,28) + XNK(3,1)*XIH*SIGN(1.0,XIH) + 
     >   ETAH*XNJ(3,1,K2)*SIGN(1.0,ETAH) + 5**(0.5D0)*MUH*
     >   XNI(1,1,J2,K2)*SIGN(1.0,MUH) 
      Q2(4,28) = Q2(4,28) + XNK(1,2)*XIH*SIGN(1.0,XIH) + 
     >   MUH*XNI(2,1,J2,K2)*SIGN(1.0,MUH) - 3**(0.5D0)*ETAH
     >   *XNJ(1,1,K2) 
      Q2(5,28) = Q2(5,28) + XNK(2,2)*XIH*SIGN(1.0,XIH) - 3*
     >   *(0.5D0)*ETAH*XNJ(2,1,K2) - 3**(0.5D0)*MUH*
     >   XNI(2,1,J2,K2) 
      Q2(6,28) = Q2(6,28) + XNK(3,2)*XIH*SIGN(1.0,XIH) - 3*
     >   *(0.5D0)*ETAH*XNJ(3,1,K2) + 5**(0.5D0)*MUH*
     >   XNI(2,1,J2,K2)*SIGN(1.0,MUH) 
      Q2(7,28) = Q2(7,28) + XNK(1,3)*XIH*SIGN(1.0,XIH) + 
     >   MUH*XNI(3,1,J2,K2)*SIGN(1.0,MUH) + 5**(0.5D0)*ETAH
     >   *XNJ(1,1,K2)*SIGN(1.0,ETAH) 
      Q2(8,28) = Q2(8,28) + XNK(2,3)*XIH*SIGN(1.0,XIH) + 5*
     >   *(0.5D0)*ETAH*XNJ(2,1,K2)*SIGN(1.0,ETAH) - 3**(0.5D0)*
     >   MUH*XNI(3,1,J2,K2) 
      Q2(9,28) = Q2(9,28) + XNK(3,3)*XIH*SIGN(1.0,XIH) + 5*
     >   *(0.5D0)*ETAH*XNJ(3,1,K2)*SIGN(1.0,ETAH) + 5**(0.5D0)*
     >   MUH*XNI(3,1,J2,K2)*SIGN(1.0,MUH) 
      Q2(10,28) = Q2(10,28) + ETAH*XNJ(1,2,K2)*SIGN(1.0,ETAH) 
     >   + MUH*XNI(1,2,J2,K2)*SIGN(1.0,MUH) - 3**(0.5D0)*
     >   XNK(1,1)*XIH 
      Q2(11,28) = Q2(11,28) + ETAH*XNJ(2,2,K2)*SIGN(1.0,ETAH) 
     >   - 3**(0.5D0)*MUH*XNI(1,2,J2,K2) - 3**(0.5D0)*XNK(2,1)*
     >   XIH 
      Q2(12,28) = Q2(12,28) + ETAH*XNJ(3,2,K2)*SIGN(1.0,ETAH) 
     >   + 5**(0.5D0)*MUH*XNI(1,2,J2,K2)*SIGN(1.0,MUH) - 3**(0.5
     >   D0)*XNK(3,1)*XIH 
      Q2(13,28) = Q2(13,28) + MUH*XNI(2,2,J2,K2)*SIGN(1.0,MUH
     >   ) - 3**(0.5D0)*ETAH*XNJ(1,2,K2) - 3**(0.5D0)*XNK(1,2)*
     >   XIH 
      Q2(14,28) = Q2(14,28) - 3**(0.5D0)*ETAH*XNJ(2,2,K2) - 3
     >   **(0.5D0)*MUH*XNI(2,2,J2,K2) - 3**(0.5D0)*XNK(2,2)*
     >   XIH 
      Q2(15,28) = Q2(15,28) + 5**(0.5D0)*MUH*XNI(2,2,J2,K2)*SIG
     >   N(1.0,MUH) - 3**(0.5D0)*ETAH*XNJ(3,2,K2) - 3**(0.5D0)*
     >   XNK(3,2)*XIH 
      Q2(16,28) = Q2(16,28) + MUH*XNI(3,2,J2,K2)*SIGN(1.0,MUH
     >   ) + 5**(0.5D0)*ETAH*XNJ(1,2,K2)*SIGN(1.0,ETAH) - 3**(0.5
     >   D0)*XNK(1,3)*XIH 
      Q2(17,28) = Q2(17,28) + 5**(0.5D0)*ETAH*XNJ(2,2,K2)*SIGN(
     >   1.0,ETAH) - 3**(0.5D0)*MUH*XNI(3,2,J2,K2) - 3**(0.5D0)*
     >   XNK(2,3)*XIH 
      Q2(18,28) = Q2(18,28) + 5**(0.5D0)*ETAH*XNJ(3,2,K2)*SIGN(
     >   1.0,ETAH) + 5**(0.5D0)*MUH*XNI(3,2,J2,K2)*SIGN(1.0,MUH
     >   ) - 3**(0.5D0)*XNK(3,3)*XIH 
      Q2(19,28) = Q2(19,28) + ETAH*XNJ(1,3,K2)*SIGN(1.0,ETAH) 
     >   + MUH*XNI(1,3,J2,K2)*SIGN(1.0,MUH) + 5**(0.5D0)*XNK(1,1)
     >   *XIH*SIGN(1.0,XIH) 
      Q2(20,28) = Q2(20,28) + ETAH*XNJ(2,3,K2)*SIGN(1.0,ETAH) 
     >   - 3**(0.5D0)*MUH*XNI(1,3,J2,K2) + 5**(0.5D0)*XNK(2,1)*
     >   XIH*SIGN(1.0,XIH) 
      Q2(21,28) = Q2(21,28) + ETAH*XNJ(3,3,K2)*SIGN(1.0,ETAH) 
     >   + 5**(0.5D0)*MUH*XNI(1,3,J2,K2)*SIGN(1.0,MUH) + 5**(0.5
     >   D0)*XNK(3,1)*XIH*SIGN(1.0,XIH) 
      Q2(22,28) = Q2(22,28) + MUH*XNI(2,3,J2,K2)*SIGN(1.0,MUH
     >   ) - 3**(0.5D0)*ETAH*XNJ(1,3,K2) + 5**(0.5D0)*XNK(1,2)*
     >   XIH*SIGN(1.0,XIH) 
      Q2(23,28) = Q2(23,28) + SIGN(1.0,XIH)*(3**(0.5D0)*15**(0.5D
     >   0)*XNK(2,2)*XIH)/3 - 3**(0.5D0)*MUH*XNI(2,3,J2,K2) 
     >   - 3**(0.5D0)*ETAH*XNJ(2,3,K2) 
      Q2(24,28) = Q2(24,28) + 5**(0.5D0)*MUH*XNI(2,3,J2,K2)*SIG
     >   N(1.0,MUH) - 3**(0.5D0)*ETAH*XNJ(3,3,K2) + 5**(0.5D0)*X
     >   NK(3,2)*XIH*SIGN(1.0,XIH) 
      Q2(25,28) = Q2(25,28) + MUH*XNI(3,3,J2,K2)*SIGN(1.0,MUH
     >   ) + 5**(0.5D0)*ETAH*XNJ(1,3,K2)*SIGN(1.0,ETAH) + 5**(0.5
     >   D0)*XNK(1,3)*XIH*SIGN(1.0,XIH) 
      Q2(26,28) = Q2(26,28) + 5**(0.5D0)*ETAH*XNJ(2,3,K2)*SIGN(
     >   1.0,ETAH) - 3**(0.5D0)*MUH*XNI(3,3,J2,K2) + 5**(0.5D0)*X
     >   NK(2,3)*XIH*SIGN(1.0,XIH) 
      Q2(27,28) = Q2(27,28) + 5**(0.5D0)*ETAH*XNJ(3,3,K2)*SIGN(
     >   1.0,ETAH) + 5**(0.5D0)*MUH*XNI(3,3,J2,K2)*SIGN(1.0,MUH
     >   ) + 5**(0.5D0)*XNK(3,3)*XIH*SIGN(1.0,XIH) 
      ELSE IF (IELEM.EQ.4) THEN
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
      Q2(1,1) = VT + ABS(ETAH) + ABS(MUH) + ABS(XIH)
      Q2(1,2) = 3**(0.5D0)*MUH 
      Q2(1,3) = 5**(0.5D0)*ABS(MUH) 
      Q2(1,4) = 7**(0.5D0)*MUH 
      Q2(1,5) = 3**(0.5D0)*ETAH 
      Q2(1,9) = 5**(0.5D0)*ABS(ETAH) 
      Q2(1,13) = 7**(0.5D0)*ETAH 
      Q2(1,17) = 3**(0.5D0)*XIH 
      Q2(1,33) = 5**(0.5D0)*ABS(XIH) 
      Q2(1,49) = 7**(0.5D0)*XIH 
      Q2(2,1) = -3**(0.5D0)*MUH 
      Q2(2,2) = VT + ABS(ETAH) + 3*ABS(MUH) + ABS(XIH)
      Q2(2,3) = 15**(0.5D0)*MUH 
      Q2(2,4) = 21**(0.5D0)*ABS(MUH) 
      Q2(2,6) = 3**(0.5D0)*ETAH 
      Q2(2,10) = 5**(0.5D0)*ABS(ETAH) 
      Q2(2,14) = 7**(0.5D0)*ETAH 
      Q2(2,18) = 3**(0.5D0)*XIH 
      Q2(2,34) = 5**(0.5D0)*ABS(XIH) 
      Q2(2,50) = 7**(0.5D0)*XIH 
      Q2(3,1) = 5**(0.5D0)*ABS(MUH) 
      Q2(3,2) = -15**(0.5D0)*MUH 
      Q2(3,3) = VT + ABS(ETAH) + 5*ABS(MUH) + ABS(XIH)
      Q2(3,4) = 35**(0.5D0)*MUH 
      Q2(3,7) = 3**(0.5D0)*ETAH 
      Q2(3,11) = 5**(0.5D0)*ABS(ETAH) 
      Q2(3,15) = 7**(0.5D0)*ETAH 
      Q2(3,19) = 3**(0.5D0)*XIH 
      Q2(3,35) = 5**(0.5D0)*ABS(XIH) 
      Q2(3,51) = 7**(0.5D0)*XIH 
      Q2(4,1) = -7**(0.5D0)*MUH 
      Q2(4,2) = 21**(0.5D0)*ABS(MUH) 
      Q2(4,3) = -35**(0.5D0)*MUH 
      Q2(4,4) = VT + ABS(ETAH) + 7*ABS(MUH) + ABS(XIH)
      Q2(4,8) = 3**(0.5D0)*ETAH 
      Q2(4,12) = 5**(0.5D0)*ABS(ETAH) 
      Q2(4,16) = 7**(0.5D0)*ETAH 
      Q2(4,20) = 3**(0.5D0)*XIH 
      Q2(4,36) = 5**(0.5D0)*ABS(XIH) 
      Q2(4,52) = 7**(0.5D0)*XIH 
      Q2(5,1) = -3**(0.5D0)*ETAH 
      Q2(5,5) = VT + 3*ABS(ETAH) + ABS(MUH) + ABS(XIH)
      Q2(5,6) = 3**(0.5D0)*MUH 
      Q2(5,7) = 5**(0.5D0)*ABS(MUH) 
      Q2(5,8) = 7**(0.5D0)*MUH 
      Q2(5,9) = 15**(0.5D0)*ETAH 
      Q2(5,13) = 21**(0.5D0)*ABS(ETAH) 
      Q2(5,21) = 3**(0.5D0)*XIH 
      Q2(5,37) = 5**(0.5D0)*ABS(XIH) 
      Q2(5,53) = 7**(0.5D0)*XIH 
      Q2(6,2) = -3**(0.5D0)*ETAH 
      Q2(6,5) = -3**(0.5D0)*MUH 
      Q2(6,6) = VT + 3*ABS(ETAH) + 3*ABS(MUH) + ABS(XIH) 
      Q2(6,7) = 15**(0.5D0)*MUH 
      Q2(6,8) = 21**(0.5D0)*ABS(MUH) 
      Q2(6,10) = 15**(0.5D0)*ETAH 
      Q2(6,14) = 21**(0.5D0)*ABS(ETAH) 
      Q2(6,22) = 3**(0.5D0)*XIH 
      Q2(6,38) = 5**(0.5D0)*ABS(XIH) 
      Q2(6,54) = 7**(0.5D0)*XIH 
      Q2(7,3) = -3**(0.5D0)*ETAH 
      Q2(7,5) = 5**(0.5D0)*ABS(MUH) 
      Q2(7,6) = -15**(0.5D0)*MUH 
      Q2(7,7) = VT + 3*ABS(ETAH) + 5*ABS(MUH) + ABS(XIH) 
      Q2(7,8) = 35**(0.5D0)*MUH 
      Q2(7,11) = 15**(0.5D0)*ETAH 
      Q2(7,15) = 21**(0.5D0)*ABS(ETAH) 
      Q2(7,23) = 3**(0.5D0)*XIH 
      Q2(7,39) = 5**(0.5D0)*ABS(XIH) 
      Q2(7,55) = 7**(0.5D0)*XIH 
      Q2(8,4) = -3**(0.5D0)*ETAH 
      Q2(8,5) = -7**(0.5D0)*MUH 
      Q2(8,6) = 21**(0.5D0)*ABS(MUH) 
      Q2(8,7) = -35**(0.5D0)*MUH 
      Q2(8,8) = VT + 3*ABS(ETAH) + 7*ABS(MUH) + ABS(XIH) 
      Q2(8,12) = 15**(0.5D0)*ETAH 
      Q2(8,16) = 21**(0.5D0)*ABS(ETAH) 
      Q2(8,24) = 3**(0.5D0)*XIH 
      Q2(8,40) = 5**(0.5D0)*ABS(XIH) 
      Q2(8,56) = 7**(0.5D0)*XIH 
      Q2(9,1) = 5**(0.5D0)*ABS(ETAH) 
      Q2(9,5) = -15**(0.5D0)*ETAH 
      Q2(9,9) = VT + 5*ABS(ETAH) + ABS(MUH) + ABS(XIH)
      Q2(9,10) = 3**(0.5D0)*MUH 
      Q2(9,11) = 5**(0.5D0)*ABS(MUH) 
      Q2(9,12) = 7**(0.5D0)*MUH 
      Q2(9,13) = 35**(0.5D0)*ETAH 
      Q2(9,25) = 3**(0.5D0)*XIH 
      Q2(9,41) = 5**(0.5D0)*ABS(XIH) 
      Q2(9,57) = 7**(0.5D0)*XIH 
      Q2(10,2) = 5**(0.5D0)*ABS(ETAH) 
      Q2(10,6) = -15**(0.5D0)*ETAH 
      Q2(10,9) = -3**(0.5D0)*MUH 
      Q2(10,10) = VT + 5*ABS(ETAH) + 3*ABS(MUH) + ABS(
     >   XIH) 
      Q2(10,11) = 15**(0.5D0)*MUH 
      Q2(10,12) = 21**(0.5D0)*ABS(MUH) 
      Q2(10,14) = 35**(0.5D0)*ETAH 
      Q2(10,26) = 3**(0.5D0)*XIH 
      Q2(10,42) = 5**(0.5D0)*ABS(XIH) 
      Q2(10,58) = 7**(0.5D0)*XIH 
      Q2(11,3) = 5**(0.5D0)*ABS(ETAH) 
      Q2(11,7) = -15**(0.5D0)*ETAH 
      Q2(11,9) = 5**(0.5D0)*ABS(MUH) 
      Q2(11,10) = -15**(0.5D0)*MUH 
      Q2(11,11) = VT + 5*ABS(ETAH) + 5*ABS(MUH) + ABS(
     >   XIH) 
      Q2(11,12) = 35**(0.5D0)*MUH 
      Q2(11,15) = 35**(0.5D0)*ETAH 
      Q2(11,27) = 3**(0.5D0)*XIH 
      Q2(11,43) = 5**(0.5D0)*ABS(XIH) 
      Q2(11,59) = 7**(0.5D0)*XIH 
      Q2(12,4) = 5**(0.5D0)*ABS(ETAH) 
      Q2(12,8) = -15**(0.5D0)*ETAH 
      Q2(12,9) = -7**(0.5D0)*MUH 
      Q2(12,10) = 21**(0.5D0)*ABS(MUH) 
      Q2(12,11) = -35**(0.5D0)*MUH 
      Q2(12,12) = VT + 5*ABS(ETAH) + 7*ABS(MUH) + ABS(
     >   XIH) 
      Q2(12,16) = 35**(0.5D0)*ETAH 
      Q2(12,28) = 3**(0.5D0)*XIH 
      Q2(12,44) = 5**(0.5D0)*ABS(XIH) 
      Q2(12,60) = 7**(0.5D0)*XIH 
      Q2(13,1) = -7**(0.5D0)*ETAH 
      Q2(13,5) = 21**(0.5D0)*ABS(ETAH) 
      Q2(13,9) = -35**(0.5D0)*ETAH 
      Q2(13,13) = VT + 7*ABS(ETAH) + ABS(MUH) + ABS(XIH) 
      Q2(13,14) = 3**(0.5D0)*MUH 
      Q2(13,15) = 5**(0.5D0)*ABS(MUH) 
      Q2(13,16) = 7**(0.5D0)*MUH 
      Q2(13,29) = 3**(0.5D0)*XIH 
      Q2(13,45) = 5**(0.5D0)*ABS(XIH) 
      Q2(13,61) = 7**(0.5D0)*XIH 
      Q2(14,2) = -7**(0.5D0)*ETAH 
      Q2(14,6) = 21**(0.5D0)*ABS(ETAH) 
      Q2(14,10) = -35**(0.5D0)*ETAH 
      Q2(14,13) = -3**(0.5D0)*MUH 
      Q2(14,14) = VT + 7*ABS(ETAH) + 3*ABS(MUH) + ABS(
     >   XIH) 
      Q2(14,15) = 15**(0.5D0)*MUH 
      Q2(14,16) = 21**(0.5D0)*ABS(MUH) 
      Q2(14,30) = 3**(0.5D0)*XIH 
      Q2(14,46) = 5**(0.5D0)*ABS(XIH) 
      Q2(14,62) = 7**(0.5D0)*XIH 
      Q2(15,3) = -7**(0.5D0)*ETAH 
      Q2(15,7) = 21**(0.5D0)*ABS(ETAH) 
      Q2(15,11) = -35**(0.5D0)*ETAH 
      Q2(15,13) = 5**(0.5D0)*ABS(MUH) 
      Q2(15,14) = -15**(0.5D0)*MUH 
      Q2(15,15) = VT + 7*ABS(ETAH) + 5*ABS(MUH) + ABS(
     >   XIH) 
      Q2(15,16) = 35**(0.5D0)*MUH 
      Q2(15,31) = 3**(0.5D0)*XIH 
      Q2(15,47) = 5**(0.5D0)*ABS(XIH) 
      Q2(15,63) = 7**(0.5D0)*XIH 
      Q2(16,4) = -7**(0.5D0)*ETAH 
      Q2(16,8) = 21**(0.5D0)*ABS(ETAH) 
      Q2(16,12) = -35**(0.5D0)*ETAH 
      Q2(16,13) = -7**(0.5D0)*MUH 
      Q2(16,14) = 21**(0.5D0)*ABS(MUH) 
      Q2(16,15) = -35**(0.5D0)*MUH 
      Q2(16,16) = VT + 7*ABS(ETAH) + 7*ABS(MUH) + ABS(
     >   XIH) 
      Q2(16,32) = 3**(0.5D0)*XIH 
      Q2(16,48) = 5**(0.5D0)*ABS(XIH) 
      Q2(16,64) = 7**(0.5D0)*XIH 
      Q2(17,1) = -3**(0.5D0)*XIH 
      Q2(17,17) = VT + ABS(ETAH) + ABS(MUH) + 3*ABS(XIH) 
      Q2(17,18) = 3**(0.5D0)*MUH 
      Q2(17,19) = 5**(0.5D0)*ABS(MUH) 
      Q2(17,20) = 7**(0.5D0)*MUH 
      Q2(17,21) = 3**(0.5D0)*ETAH 
      Q2(17,25) = 5**(0.5D0)*ABS(ETAH) 
      Q2(17,29) = 7**(0.5D0)*ETAH 
      Q2(17,33) = 15**(0.5D0)*XIH 
      Q2(17,49) = 21**(0.5D0)*ABS(XIH) 
      Q2(18,2) = -3**(0.5D0)*XIH 
      Q2(18,17) = -3**(0.5D0)*MUH 
      Q2(18,18) = VT + ABS(ETAH) + 3*ABS(MUH) + 3*ABS(
     >   XIH) 
      Q2(18,19) = 15**(0.5D0)*MUH 
      Q2(18,20) = 21**(0.5D0)*ABS(MUH) 
      Q2(18,22) = 3**(0.5D0)*ETAH 
      Q2(18,26) = 5**(0.5D0)*ABS(ETAH) 
      Q2(18,30) = 7**(0.5D0)*ETAH 
      Q2(18,34) = 15**(0.5D0)*XIH 
      Q2(18,50) = 21**(0.5D0)*ABS(XIH) 
      Q2(19,3) = -3**(0.5D0)*XIH 
      Q2(19,17) = 5**(0.5D0)*ABS(MUH) 
      Q2(19,18) = -15**(0.5D0)*MUH 
      Q2(19,19) = VT + ABS(ETAH) + 5*ABS(MUH) + 3*ABS(
     >   XIH) 
      Q2(19,20) = 35**(0.5D0)*MUH 
      Q2(19,23) = 3**(0.5D0)*ETAH 
      Q2(19,27) = 5**(0.5D0)*ABS(ETAH) 
      Q2(19,31) = 7**(0.5D0)*ETAH 
      Q2(19,35) = 15**(0.5D0)*XIH 
      Q2(19,51) = 21**(0.5D0)*ABS(XIH) 
      Q2(20,4) = -3**(0.5D0)*XIH 
      Q2(20,17) = -7**(0.5D0)*MUH 
      Q2(20,18) = 21**(0.5D0)*ABS(MUH) 
      Q2(20,19) = -35**(0.5D0)*MUH 
      Q2(20,20) = VT + ABS(ETAH) + 7*ABS(MUH) + 3*ABS(
     >   XIH) 
      Q2(20,24) = 3**(0.5D0)*ETAH 
      Q2(20,28) = 5**(0.5D0)*ABS(ETAH) 
      Q2(20,32) = 7**(0.5D0)*ETAH 
      Q2(20,36) = 15**(0.5D0)*XIH 
      Q2(20,52) = 21**(0.5D0)*ABS(XIH) 
      Q2(21,5) = -3**(0.5D0)*XIH 
      Q2(21,17) = -3**(0.5D0)*ETAH 
      Q2(21,21) = VT + 3*ABS(ETAH) + ABS(MUH) + 3*ABS(
     >   XIH) 
      Q2(21,22) = 3**(0.5D0)*MUH 
      Q2(21,23) = 5**(0.5D0)*ABS(MUH) 
      Q2(21,24) = 7**(0.5D0)*MUH 
      Q2(21,25) = 15**(0.5D0)*ETAH 
      Q2(21,29) = 21**(0.5D0)*ABS(ETAH) 
      Q2(21,37) = 15**(0.5D0)*XIH 
      Q2(21,53) = 21**(0.5D0)*ABS(XIH) 
      Q2(22,6) = -3**(0.5D0)*XIH 
      Q2(22,18) = -3**(0.5D0)*ETAH 
      Q2(22,21) = -3**(0.5D0)*MUH 
      Q2(22,22) = VT + 3*ABS(ETAH) + 3*ABS(MUH) + 3*ABS(
     >   XIH) 
      Q2(22,23) = 15**(0.5D0)*MUH 
      Q2(22,24) = 21**(0.5D0)*ABS(MUH) 
      Q2(22,26) = 15**(0.5D0)*ETAH 
      Q2(22,30) = 21**(0.5D0)*ABS(ETAH) 
      Q2(22,38) = 15**(0.5D0)*XIH 
      Q2(22,54) = 21**(0.5D0)*ABS(XIH) 
      Q2(23,7) = -3**(0.5D0)*XIH 
      Q2(23,19) = -3**(0.5D0)*ETAH 
      Q2(23,21) = 5**(0.5D0)*ABS(MUH) 
      Q2(23,22) = -15**(0.5D0)*MUH 
      Q2(23,23) = VT + 3*ABS(ETAH) + 5*ABS(MUH) + 3*ABS(
     >   XIH) 
      Q2(23,24) = 35**(0.5D0)*MUH 
      Q2(23,27) = 15**(0.5D0)*ETAH 
      Q2(23,31) = 21**(0.5D0)*ABS(ETAH) 
      Q2(23,39) = 15**(0.5D0)*XIH 
      Q2(23,55) = 21**(0.5D0)*ABS(XIH) 
      Q2(24,8) = -3**(0.5D0)*XIH 
      Q2(24,20) = -3**(0.5D0)*ETAH 
      Q2(24,21) = -7**(0.5D0)*MUH 
      Q2(24,22) = 21**(0.5D0)*ABS(MUH) 
      Q2(24,23) = -35**(0.5D0)*MUH 
      Q2(24,24) = VT + 3*ABS(ETAH) + 7*ABS(MUH) + 3*ABS(
     >   XIH) 
      Q2(24,28) = 15**(0.5D0)*ETAH 
      Q2(24,32) = 21**(0.5D0)*ABS(ETAH) 
      Q2(24,40) = 15**(0.5D0)*XIH 
      Q2(24,56) = 21**(0.5D0)*ABS(XIH) 
      Q2(25,9) = -3**(0.5D0)*XIH 
      Q2(25,17) = 5**(0.5D0)*ABS(ETAH) 
      Q2(25,21) = -15**(0.5D0)*ETAH 
      Q2(25,25) = VT + 5*ABS(ETAH) + ABS(MUH) + 3*ABS(
     >   XIH) 
      Q2(25,26) = 3**(0.5D0)*MUH 
      Q2(25,27) = 5**(0.5D0)*ABS(MUH) 
      Q2(25,28) = 7**(0.5D0)*MUH 
      Q2(25,29) = 35**(0.5D0)*ETAH 
      Q2(25,41) = 2*15**(0.5D0)*XIH - 3**(0.5D0)*5**(0.5D0)*
     >   XIH 
      Q2(25,57) = 3**(0.5D0)*7**(0.5D0)*ABS(XIH) 
      Q2(26,10) = -3**(0.5D0)*XIH 
      Q2(26,18) = 5**(0.5D0)*ABS(ETAH) 
      Q2(26,22) = -15**(0.5D0)*ETAH 
      Q2(26,25) = -3**(0.5D0)*MUH 
      Q2(26,26) = VT + 5*ABS(ETAH) + 3*ABS(MUH) + 3*ABS(
     >   XIH) 
      Q2(26,27) = 15**(0.5D0)*MUH 
      Q2(26,28) = 21**(0.5D0)*ABS(MUH) 
      Q2(26,30) = 35**(0.5D0)*ETAH 
      Q2(26,42) = 2*15**(0.5D0)*XIH - 3**(0.5D0)*5**(0.5D0)*
     >   XIH 
      Q2(26,58) = 3**(0.5D0)*7**(0.5D0)*ABS(XIH) 
      Q2(27,11) = -3**(0.5D0)*XIH 
      Q2(27,19) = 5**(0.5D0)*ABS(ETAH) 
      Q2(27,23) = -15**(0.5D0)*ETAH 
      Q2(27,25) = 5**(0.5D0)*ABS(MUH) 
      Q2(27,26) = -15**(0.5D0)*MUH 
      Q2(27,27) = VT + 5*ABS(ETAH) + 5*ABS(MUH) + 3*ABS(
     >   XIH) 
      Q2(27,28) = 35**(0.5D0)*MUH 
      Q2(27,31) = 35**(0.5D0)*ETAH 
      Q2(27,43) = 2*15**(0.5D0)*XIH - 3**(0.5D0)*5**(0.5D0)*
     >   XIH 
      Q2(27,59) = 3**(0.5D0)*7**(0.5D0)*ABS(XIH) 
      Q2(28,12) = -3**(0.5D0)*XIH 
      Q2(28,20) = 5**(0.5D0)*ABS(ETAH) 
      Q2(28,24) = -15**(0.5D0)*ETAH 
      Q2(28,25) = -7**(0.5D0)*MUH 
      Q2(28,26) = 21**(0.5D0)*ABS(MUH) 
      Q2(28,27) = -35**(0.5D0)*MUH 
      Q2(28,28) = VT + 5*ABS(ETAH) + 7*ABS(MUH) + 3*ABS(
     >   XIH) 
      Q2(28,32) = 35**(0.5D0)*ETAH 
      Q2(28,44) = 2*15**(0.5D0)*XIH - 3**(0.5D0)*5**(0.5D0)*
     >   XIH 
      Q2(28,60) = 3**(0.5D0)*7**(0.5D0)*ABS(XIH) 
      Q2(29,13) = -3**(0.5D0)*XIH 
      Q2(29,17) = -7**(0.5D0)*ETAH 
      Q2(29,21) = 21**(0.5D0)*ABS(ETAH) 
      Q2(29,25) = -35**(0.5D0)*ETAH 
      Q2(29,29) = VT + 7*ABS(ETAH) + ABS(MUH) + 3*ABS(
     >   XIH) 
      Q2(29,30) = 3**(0.5D0)*MUH 
      Q2(29,31) = 5**(0.5D0)*ABS(MUH) 
      Q2(29,32) = 7**(0.5D0)*MUH 
      Q2(29,45) = 2*15**(0.5D0)*XIH - 3**(0.5D0)*5**(0.5D0)*
     >   XIH 
      Q2(29,61) = 3**(0.5D0)*7**(0.5D0)*ABS(XIH) 
      Q2(30,14) = -3**(0.5D0)*XIH 
      Q2(30,18) = -7**(0.5D0)*ETAH 
      Q2(30,22) = 21**(0.5D0)*ABS(ETAH) 
      Q2(30,26) = -35**(0.5D0)*ETAH 
      Q2(30,29) = -3**(0.5D0)*MUH 
      Q2(30,30) = VT + 7*ABS(ETAH) + 3*ABS(MUH) + 3*ABS(
     >   XIH) 
      Q2(30,31) = 15**(0.5D0)*MUH 
      Q2(30,32) = 21**(0.5D0)*ABS(MUH) 
      Q2(30,46) = 2*15**(0.5D0)*XIH - 3**(0.5D0)*5**(0.5D0)*
     >   XIH 
      Q2(30,62) = 3**(0.5D0)*7**(0.5D0)*ABS(XIH) 
      Q2(31,15) = -3**(0.5D0)*XIH 
      Q2(31,19) = -7**(0.5D0)*ETAH 
      Q2(31,23) = 21**(0.5D0)*ABS(ETAH) 
      Q2(31,27) = -35**(0.5D0)*ETAH 
      Q2(31,29) = 5**(0.5D0)*ABS(MUH) 
      Q2(31,30) = -15**(0.5D0)*MUH 
      Q2(31,31) = VT + 7*ABS(ETAH) + 5*ABS(MUH) + 3*ABS(
     >   XIH) 
      Q2(31,32) = 35**(0.5D0)*MUH 
      Q2(31,47) = 2*15**(0.5D0)*XIH - 3**(0.5D0)*5**(0.5D0)*
     >   XIH 
      Q2(31,63) = 3**(0.5D0)*7**(0.5D0)*ABS(XIH) 
      Q2(32,16) = -3**(0.5D0)*XIH 
      Q2(32,20) = -7**(0.5D0)*ETAH 
      Q2(32,24) = 21**(0.5D0)*ABS(ETAH) 
      Q2(32,28) = -35**(0.5D0)*ETAH 
      Q2(32,29) = -7**(0.5D0)*MUH 
      Q2(32,30) = 21**(0.5D0)*ABS(MUH) 
      Q2(32,31) = -35**(0.5D0)*MUH 
      Q2(32,32) = VT + 7*ABS(ETAH) + 7*ABS(MUH) + 3*ABS(
     >   XIH) 
      Q2(32,48) = 2*15**(0.5D0)*XIH - 3**(0.5D0)*5**(0.5D0)*
     >   XIH 
      Q2(32,64) = 3**(0.5D0)*7**(0.5D0)*ABS(XIH) 
      Q2(33,1) = 5**(0.5D0)*ABS(XIH) 
      Q2(33,17) = -15**(0.5D0)*XIH 
      Q2(33,33) = VT + ABS(ETAH) + ABS(MUH) + 5*ABS(XIH) 
      Q2(33,34) = 3**(0.5D0)*MUH 
      Q2(33,35) = 5**(0.5D0)*ABS(MUH) 
      Q2(33,36) = 7**(0.5D0)*MUH 
      Q2(33,37) = 3**(0.5D0)*ETAH 
      Q2(33,41) = 5**(0.5D0)*ABS(ETAH) 
      Q2(33,45) = 7**(0.5D0)*ETAH 
      Q2(33,49) = 35**(0.5D0)*XIH 
      Q2(34,2) = 5**(0.5D0)*ABS(XIH) 
      Q2(34,18) = -15**(0.5D0)*XIH 
      Q2(34,33) = -3**(0.5D0)*MUH 
      Q2(34,34) = VT + ABS(ETAH) + 3*ABS(MUH) + 5*ABS(
     >   XIH) 
      Q2(34,35) = 2*15**(0.5D0)*MUH - 3**(0.5D0)*5**(0.5D0)*
     >   MUH 
      Q2(34,36) = 3**(0.5D0)*7**(0.5D0)*ABS(MUH) 
      Q2(34,38) = 3**(0.5D0)*ETAH 
      Q2(34,42) = 5**(0.5D0)*ABS(ETAH) 
      Q2(34,46) = 7**(0.5D0)*ETAH 
      Q2(34,50) = 35**(0.5D0)*XIH 
      Q2(35,3) = 5**(0.5D0)*ABS(XIH) 
      Q2(35,19) = -15**(0.5D0)*XIH 
      Q2(35,33) = 5**(0.5D0)*ABS(MUH) 
      Q2(35,34) = -3**(0.5D0)*5**(0.5D0)*MUH 
      Q2(35,35) = VT + ABS(ETAH) + 5*ABS(MUH) + 5*ABS(
     >   XIH) 
      Q2(35,36) = 2*35**(0.5D0)*MUH - 5**(0.5D0)*7**(0.5D0)*
     >   MUH 
      Q2(35,39) = 3**(0.5D0)*ETAH 
      Q2(35,43) = 5**(0.5D0)*ABS(ETAH) 
      Q2(35,47) = 7**(0.5D0)*ETAH 
      Q2(35,51) = 35**(0.5D0)*XIH 
      Q2(36,4) = 5**(0.5D0)*ABS(XIH) 
      Q2(36,20) = -15**(0.5D0)*XIH 
      Q2(36,33) = -7**(0.5D0)*MUH 
      Q2(36,34) = 3**(0.5D0)*7**(0.5D0)*ABS(MUH) 
      Q2(36,35) = -5**(0.5D0)*7**(0.5D0)*MUH 
      Q2(36,36) = VT + ABS(ETAH) + 7*ABS(MUH) + 5*ABS(
     >   XIH) 
      Q2(36,40) = 3**(0.5D0)*ETAH 
      Q2(36,44) = 5**(0.5D0)*ABS(ETAH) 
      Q2(36,48) = 7**(0.5D0)*ETAH 
      Q2(36,52) = 35**(0.5D0)*XIH 
      Q2(37,5) = 5**(0.5D0)*ABS(XIH) 
      Q2(37,21) = -15**(0.5D0)*XIH 
      Q2(37,33) = -3**(0.5D0)*ETAH 
      Q2(37,37) = VT + 3*ABS(ETAH) + ABS(MUH) + 5*ABS(
     >   XIH) 
      Q2(37,38) = 3**(0.5D0)*MUH 
      Q2(37,39) = 5**(0.5D0)*ABS(MUH) 
      Q2(37,40) = 7**(0.5D0)*MUH 
      Q2(37,41) = 2*15**(0.5D0)*ETAH - 3**(0.5D0)*5**(0.5D0)*
     >   ETAH 
      Q2(37,45) = 3**(0.5D0)*7**(0.5D0)*ABS(ETAH) 
      Q2(37,53) = 35**(0.5D0)*XIH 
      Q2(38,6) = 5**(0.5D0)*ABS(XIH) 
      Q2(38,22) = -15**(0.5D0)*XIH 
      Q2(38,34) = -3**(0.5D0)*ETAH 
      Q2(38,37) = -3**(0.5D0)*MUH 
      Q2(38,38) = VT + 3*ABS(ETAH) + 3*ABS(MUH) + 5*ABS(
     >   XIH) 
      Q2(38,39) = 2*15**(0.5D0)*MUH - 3**(0.5D0)*5**(0.5D0)*
     >   MUH 
      Q2(38,40) = 3**(0.5D0)*7**(0.5D0)*ABS(MUH) 
      Q2(38,42) = 2*15**(0.5D0)*ETAH - 3**(0.5D0)*5**(0.5D0)*
     >   ETAH 
      Q2(38,46) = 3**(0.5D0)*7**(0.5D0)*ABS(ETAH) 
      Q2(38,54) = 35**(0.5D0)*XIH 
      Q2(39,7) = 5**(0.5D0)*ABS(XIH) 
      Q2(39,23) = -15**(0.5D0)*XIH 
      Q2(39,35) = -3**(0.5D0)*ETAH 
      Q2(39,37) = 5**(0.5D0)*ABS(MUH) 
      Q2(39,38) = -3**(0.5D0)*5**(0.5D0)*MUH 
      Q2(39,39) = VT + 3*ABS(ETAH) + 5*ABS(MUH) + 5*ABS(
     >   XIH) 
      Q2(39,40) = 2*35**(0.5D0)*MUH - 5**(0.5D0)*7**(0.5D0)*
     >   MUH 
      Q2(39,43) = 2*15**(0.5D0)*ETAH - 3**(0.5D0)*5**(0.5D0)*
     >   ETAH 
      Q2(39,47) = 3**(0.5D0)*7**(0.5D0)*ABS(ETAH) 
      Q2(39,55) = 35**(0.5D0)*XIH 
      Q2(40,8) = 5**(0.5D0)*ABS(XIH) 
      Q2(40,24) = -15**(0.5D0)*XIH 
      Q2(40,36) = -3**(0.5D0)*ETAH 
      Q2(40,37) = -7**(0.5D0)*MUH 
      Q2(40,38) = 3**(0.5D0)*7**(0.5D0)*ABS(MUH) 
      Q2(40,39) = -5**(0.5D0)*7**(0.5D0)*MUH 
      Q2(40,40) = VT + 3*ABS(ETAH) + 7*ABS(MUH) + 5*ABS(
     >   XIH) 
      Q2(40,44) = 2*15**(0.5D0)*ETAH - 3**(0.5D0)*5**(0.5D0)*
     >   ETAH 
      Q2(40,48) = 3**(0.5D0)*7**(0.5D0)*ABS(ETAH) 
      Q2(40,56) = 35**(0.5D0)*XIH 
      Q2(41,9) = 5**(0.5D0)*ABS(XIH) 
      Q2(41,25) = -3**(0.5D0)*5**(0.5D0)*XIH 
      Q2(41,33) = 5**(0.5D0)*ABS(ETAH) 
      Q2(41,37) = -3**(0.5D0)*5**(0.5D0)*ETAH 
      Q2(41,41) = VT + 5*ABS(ETAH) + ABS(MUH) + 5*ABS(
     >   XIH) 
      Q2(41,42) = 3**(0.5D0)*MUH 
      Q2(41,43) = 5**(0.5D0)*ABS(MUH) 
      Q2(41,44) = 7**(0.5D0)*MUH 
      Q2(41,45) = 2*35**(0.5D0)*ETAH - 5**(0.5D0)*7**(0.5D0)*
     >   ETAH 
      Q2(41,57) = 2*35**(0.5D0)*XIH - 5**(0.5D0)*7**(0.5D0)*
     >   XIH 
      Q2(42,10) = 5**(0.5D0)*ABS(XIH) 
      Q2(42,26) = -3**(0.5D0)*5**(0.5D0)*XIH 
      Q2(42,34) = 5**(0.5D0)*ABS(ETAH) 
      Q2(42,38) = -3**(0.5D0)*5**(0.5D0)*ETAH 
      Q2(42,41) = -3**(0.5D0)*MUH 
      Q2(42,42) = VT + 5*ABS(ETAH) + 3*ABS(MUH) + 5*ABS(
     >   XIH) 
      Q2(42,43) = 2*15**(0.5D0)*MUH - 3**(0.5D0)*5**(0.5D0)*
     >   MUH 
      Q2(42,44) = 3**(0.5D0)*7**(0.5D0)*ABS(MUH) 
      Q2(42,46) = 2*35**(0.5D0)*ETAH - 5**(0.5D0)*7**(0.5D0)*
     >   ETAH 
      Q2(42,58) = 2*35**(0.5D0)*XIH - 5**(0.5D0)*7**(0.5D0)*
     >   XIH 
      Q2(43,11) = 5**(0.5D0)*ABS(XIH) 
      Q2(43,27) = -3**(0.5D0)*5**(0.5D0)*XIH 
      Q2(43,35) = 5**(0.5D0)*ABS(ETAH) 
      Q2(43,39) = -3**(0.5D0)*5**(0.5D0)*ETAH 
      Q2(43,41) = 5**(0.5D0)*ABS(MUH) 
      Q2(43,42) = -3**(0.5D0)*5**(0.5D0)*MUH 
      Q2(43,43) = VT + 5*ABS(ETAH) + 5*ABS(MUH) + 5*ABS(
     >   XIH) 
      Q2(43,44) = 2*35**(0.5D0)*MUH - 5**(0.5D0)*7**(0.5D0)*
     >   MUH 
      Q2(43,47) = 2*35**(0.5D0)*ETAH - 5**(0.5D0)*7**(0.5D0)*
     >   ETAH 
      Q2(43,59) = 2*35**(0.5D0)*XIH - 5**(0.5D0)*7**(0.5D0)*
     >   XIH 
      Q2(44,12) = 5**(0.5D0)*ABS(XIH) 
      Q2(44,28) = -3**(0.5D0)*5**(0.5D0)*XIH 
      Q2(44,36) = 5**(0.5D0)*ABS(ETAH) 
      Q2(44,40) = -3**(0.5D0)*5**(0.5D0)*ETAH 
      Q2(44,41) = -7**(0.5D0)*MUH 
      Q2(44,42) = 3**(0.5D0)*7**(0.5D0)*ABS(MUH) 
      Q2(44,43) = -5**(0.5D0)*7**(0.5D0)*MUH 
      Q2(44,44) = VT + 5*ABS(ETAH) + 7*ABS(MUH) + 5*ABS(
     >   XIH) 
      Q2(44,48) = 2*35**(0.5D0)*ETAH - 5**(0.5D0)*7**(0.5D0)*
     >   ETAH 
      Q2(44,60) = 2*35**(0.5D0)*XIH - 5**(0.5D0)*7**(0.5D0)*
     >   XIH 
      Q2(45,13) = 5**(0.5D0)*ABS(XIH) 
      Q2(45,29) = -3**(0.5D0)*5**(0.5D0)*XIH 
      Q2(45,33) = -7**(0.5D0)*ETAH 
      Q2(45,37) = 3**(0.5D0)*7**(0.5D0)*ABS(ETAH) 
      Q2(45,41) = -5**(0.5D0)*7**(0.5D0)*ETAH 
      Q2(45,45) = VT + 7*ABS(ETAH) + ABS(MUH) + 5*ABS(
     >   XIH) 
      Q2(45,46) = 3**(0.5D0)*MUH 
      Q2(45,47) = 5**(0.5D0)*ABS(MUH) 
      Q2(45,48) = 7**(0.5D0)*MUH 
      Q2(45,61) = 2*35**(0.5D0)*XIH - 5**(0.5D0)*7**(0.5D0)*
     >   XIH 
      Q2(46,14) = 5**(0.5D0)*ABS(XIH) 
      Q2(46,30) = -3**(0.5D0)*5**(0.5D0)*XIH 
      Q2(46,34) = -7**(0.5D0)*ETAH 
      Q2(46,38) = 3**(0.5D0)*7**(0.5D0)*ABS(ETAH) 
      Q2(46,42) = -5**(0.5D0)*7**(0.5D0)*ETAH 
      Q2(46,45) = -3**(0.5D0)*MUH 
      Q2(46,46) = VT + 7*ABS(ETAH) + 3*ABS(MUH) + 5*ABS(
     >   XIH) 
      Q2(46,47) = 2*15**(0.5D0)*MUH - 3**(0.5D0)*5**(0.5D0)*
     >   MUH 
      Q2(46,48) = 3**(0.5D0)*7**(0.5D0)*ABS(MUH) 
      Q2(46,62) = 2*35**(0.5D0)*XIH - 5**(0.5D0)*7**(0.5D0)*
     >   XIH 
      Q2(47,15) = 5**(0.5D0)*ABS(XIH) 
      Q2(47,31) = -3**(0.5D0)*5**(0.5D0)*XIH 
      Q2(47,35) = -7**(0.5D0)*ETAH 
      Q2(47,39) = 3**(0.5D0)*7**(0.5D0)*ABS(ETAH) 
      Q2(47,43) = -5**(0.5D0)*7**(0.5D0)*ETAH 
      Q2(47,45) = 5**(0.5D0)*ABS(MUH) 
      Q2(47,46) = -3**(0.5D0)*5**(0.5D0)*MUH 
      Q2(47,47) = VT + 7*ABS(ETAH) + 5*ABS(MUH) + 5*ABS(
     >   XIH) 
      Q2(47,48) = 2*35**(0.5D0)*MUH - 5**(0.5D0)*7**(0.5D0)*MUH 
      Q2(47,63) = 2*35**(0.5D0)*XIH - 5**(0.5D0)*7**(0.5D0)*XIH 
      Q2(48,16) = 5**(0.5D0)*ABS(XIH) 
      Q2(48,32) = -3**(0.5D0)*5**(0.5D0)*XIH 
      Q2(48,36) = -7**(0.5D0)*ETAH 
      Q2(48,40) = 3**(0.5D0)*7**(0.5D0)*ABS(ETAH) 
      Q2(48,44) = -5**(0.5D0)*7**(0.5D0)*ETAH 
      Q2(48,45) = -7**(0.5D0)*MUH 
      Q2(48,46) = 3**(0.5D0)*7**(0.5D0)*ABS(MUH) 
      Q2(48,47) = -5**(0.5D0)*7**(0.5D0)*MUH 
      Q2(48,48) = VT + 7*ABS(ETAH) + 7*ABS(MUH) + 5*ABS(XIH) 
      Q2(48,64) = 2*35**(0.5D0)*XIH - 5**(0.5D0)*7**(0.5D0)*XIH 
      Q2(49,1) = -7**(0.5D0)*XIH 
      Q2(49,17) = 21**(0.5D0)*ABS(XIH) 
      Q2(49,33) = -35**(0.5D0)*XIH 
      Q2(49,49) = VT + ABS(ETAH) + ABS(MUH) + 7*ABS(XIH) 
      Q2(49,50) = 3**(0.5D0)*MUH 
      Q2(49,51) = 5**(0.5D0)*ABS(MUH) 
      Q2(49,52) = 7**(0.5D0)*MUH 
      Q2(49,53) = 3**(0.5D0)*ETAH 
      Q2(49,57) = 5**(0.5D0)*ABS(ETAH) 
      Q2(49,61) = 7**(0.5D0)*ETAH 
      Q2(50,2) = -7**(0.5D0)*XIH 
      Q2(50,18) = 21**(0.5D0)*ABS(XIH) 
      Q2(50,34) = -35**(0.5D0)*XIH 
      Q2(50,49) = -3**(0.5D0)*MUH 
      Q2(50,50) = VT + ABS(ETAH) + 3*ABS(MUH) + 7*ABS(XIH) 
      Q2(50,51) = 2*15**(0.5D0)*MUH - 3**(0.5D0)*5**(0.5D0)*MUH 
      Q2(50,52) = 3**(0.5D0)*7**(0.5D0)*ABS(MUH) 
      Q2(50,54) = 3**(0.5D0)*ETAH 
      Q2(50,58) = 5**(0.5D0)*ABS(ETAH) 
      Q2(50,62) = 7**(0.5D0)*ETAH 
      Q2(51,3) = -7**(0.5D0)*XIH 
      Q2(51,19) = 21**(0.5D0)*ABS(XIH) 
      Q2(51,35) = -35**(0.5D0)*XIH 
      Q2(51,49) = 5**(0.5D0)*ABS(MUH) 
      Q2(51,50) = -3**(0.5D0)*5**(0.5D0)*MUH 
      Q2(51,51) = VT + ABS(ETAH) + 5*ABS(MUH) + 7*ABS(XIH) 
      Q2(51,52) = 2*35**(0.5D0)*MUH - 5**(0.5D0)*7**(0.5D0)*MUH 
      Q2(51,55) = 3**(0.5D0)*ETAH 
      Q2(51,59) = 5**(0.5D0)*ABS(ETAH) 
      Q2(51,63) = 7**(0.5D0)*ETAH 
      Q2(52,4) = -7**(0.5D0)*XIH 
      Q2(52,20) = 21**(0.5D0)*ABS(XIH) 
      Q2(52,36) = -35**(0.5D0)*XIH 
      Q2(52,49) = -7**(0.5D0)*MUH 
      Q2(52,50) = 3**(0.5D0)*7**(0.5D0)*ABS(MUH) 
      Q2(52,51) = -5**(0.5D0)*7**(0.5D0)*MUH 
      Q2(52,52) = VT + ABS(ETAH) + 7*ABS(MUH) + 7*ABS(XIH) 
      Q2(52,56) = 3**(0.5D0)*ETAH 
      Q2(52,60) = 5**(0.5D0)*ABS(ETAH) 
      Q2(52,64) = 7**(0.5D0)*ETAH 
      Q2(53,5) = -7**(0.5D0)*XIH 
      Q2(53,21) = 21**(0.5D0)*ABS(XIH) 
      Q2(53,37) = -35**(0.5D0)*XIH 
      Q2(53,49) = -3**(0.5D0)*ETAH 
      Q2(53,53) = VT + 3*ABS(ETAH) + ABS(MUH) + 7*ABS(XIH) 
      Q2(53,54) = 3**(0.5D0)*MUH 
      Q2(53,55) = 5**(0.5D0)*ABS(MUH) 
      Q2(53,56) = 7**(0.5D0)*MUH 
      Q2(53,57) = 2*15**(0.5D0)*ETAH - 3**(0.5D0)*5**(0.5D0)*ETAH 
      Q2(53,61) = 3**(0.5D0)*7**(0.5D0)*ABS(ETAH) 
      Q2(54,6) = -7**(0.5D0)*XIH 
      Q2(54,22) = 21**(0.5D0)*ABS(XIH) 
      Q2(54,38) = -35**(0.5D0)*XIH 
      Q2(54,50) = -3**(0.5D0)*ETAH 
      Q2(54,53) = -3**(0.5D0)*MUH 
      Q2(54,54) = VT + 3*ABS(ETAH) + 3*ABS(MUH) + 7*ABS(XIH) 
      Q2(54,55) = 2*15**(0.5D0)*MUH - 3**(0.5D0)*5**(0.5D0)*MUH 
      Q2(54,56) = 3**(0.5D0)*7**(0.5D0)*ABS(MUH) 
      Q2(54,58) = 2*15**(0.5D0)*ETAH - 3**(0.5D0)*5**(0.5D0)*ETAH 
      Q2(54,62) = 3**(0.5D0)*7**(0.5D0)*ABS(ETAH) 
      Q2(55,7) = -7**(0.5D0)*XIH 
      Q2(55,23) = 21**(0.5D0)*ABS(XIH) 
      Q2(55,39) = -35**(0.5D0)*XIH 
      Q2(55,51) = -3**(0.5D0)*ETAH 
      Q2(55,53) = 5**(0.5D0)*ABS(MUH) 
      Q2(55,54) = -3**(0.5D0)*5**(0.5D0)*MUH 
      Q2(55,55) = VT + 3*ABS(ETAH) + 5*ABS(MUH) + 7*ABS(XIH) 
      Q2(55,56) = 2*35**(0.5D0)*MUH - 5**(0.5D0)*7**(0.5D0)*MUH 
      Q2(55,59) = 2*15**(0.5D0)*ETAH - 3**(0.5D0)*5**(0.5D0)*ETAH 
      Q2(55,63) = 3**(0.5D0)*7**(0.5D0)*ABS(ETAH) 
      Q2(56,8) = -7**(0.5D0)*XIH 
      Q2(56,24) = 21**(0.5D0)*ABS(XIH) 
      Q2(56,40) = -35**(0.5D0)*XIH 
      Q2(56,52) = -3**(0.5D0)*ETAH 
      Q2(56,53) = -7**(0.5D0)*MUH 
      Q2(56,54) = 3**(0.5D0)*7**(0.5D0)*ABS(MUH) 
      Q2(56,55) = -5**(0.5D0)*7**(0.5D0)*MUH 
      Q2(56,56) = VT + 3*ABS(ETAH) + 7*ABS(MUH) + 7*ABS(XIH) 
      Q2(56,60) = 2*15**(0.5D0)*ETAH - 3**(0.5D0)*5**(0.5D0)*ETAH 
      Q2(56,64) = 3**(0.5D0)*7**(0.5D0)*ABS(ETAH) 
      Q2(57,9) = -7**(0.5D0)*XIH 
      Q2(57,25) = 3**(0.5D0)*7**(0.5D0)*ABS(XIH) 
      Q2(57,41) = -5**(0.5D0)*7**(0.5D0)*XIH 
      Q2(57,49) = 5**(0.5D0)*ABS(ETAH) 
      Q2(57,53) = -3**(0.5D0)*5**(0.5D0)*ETAH 
      Q2(57,57) = VT + 5*ABS(ETAH) + ABS(MUH) + 7*ABS(XIH) 
      Q2(57,58) = 3**(0.5D0)*MUH 
      Q2(57,59) = 5**(0.5D0)*ABS(MUH) 
      Q2(57,60) = 7**(0.5D0)*MUH 
      Q2(57,61) = 2*35**(0.5D0)*ETAH - 5**(0.5D0)*7**(0.5D0)*ETAH 
      Q2(58,10) = -7**(0.5D0)*XIH 
      Q2(58,26) = 3**(0.5D0)*7**(0.5D0)*ABS(XIH) 
      Q2(58,42) = -5**(0.5D0)*7**(0.5D0)*XIH 
      Q2(58,50) = 5**(0.5D0)*ABS(ETAH) 
      Q2(58,54) = -3**(0.5D0)*5**(0.5D0)*ETAH 
      Q2(58,57) = -3**(0.5D0)*MUH 
      Q2(58,58) = VT + 5*ABS(ETAH) + 3*ABS(MUH) + 7*ABS(XIH) 
      Q2(58,59) = 2*15**(0.5D0)*MUH - 3**(0.5D0)*5**(0.5D0)*MUH 
      Q2(58,60) = 3**(0.5D0)*7**(0.5D0)*ABS(MUH) 
      Q2(58,62) = 2*35**(0.5D0)*ETAH - 5**(0.5D0)*7**(0.5D0)*ETAH 
      Q2(59,11) = -7**(0.5D0)*XIH 
      Q2(59,27) = 3**(0.5D0)*7**(0.5D0)*ABS(XIH) 
      Q2(59,43) = -5**(0.5D0)*7**(0.5D0)*XIH 
      Q2(59,51) = 5**(0.5D0)*ABS(ETAH) 
      Q2(59,55) = -3**(0.5D0)*5**(0.5D0)*ETAH 
      Q2(59,57) = 5**(0.5D0)*ABS(MUH) 
      Q2(59,58) = -3**(0.5D0)*5**(0.5D0)*MUH 
      Q2(59,59) = VT + 5*ABS(ETAH) + 5*ABS(MUH) + 7*ABS(XIH) 
      Q2(59,60) = 2*35**(0.5D0)*MUH - 5**(0.5D0)*7**(0.5D0)*MUH 
      Q2(59,63) = 2*35**(0.5D0)*ETAH - 5**(0.5D0)*7**(0.5D0)*ETAH 
      Q2(60,12) = -7**(0.5D0)*XIH 
      Q2(60,28) = 3**(0.5D0)*7**(0.5D0)*ABS(XIH) 
      Q2(60,44) = -5**(0.5D0)*7**(0.5D0)*XIH 
      Q2(60,52) = 5**(0.5D0)*ABS(ETAH) 
      Q2(60,56) = -3**(0.5D0)*5**(0.5D0)*ETAH 
      Q2(60,57) = -7**(0.5D0)*MUH 
      Q2(60,58) = 3**(0.5D0)*7**(0.5D0)*ABS(MUH) 
      Q2(60,59) = -5**(0.5D0)*7**(0.5D0)*MUH 
      Q2(60,60) = VT + 5*ABS(ETAH) + 7*ABS(MUH) + 7*ABS(XIH) 
      Q2(60,64) = 2*35**(0.5D0)*ETAH - 5**(0.5D0)*7**(0.5D0)*ETAH 
      Q2(61,13) = -7**(0.5D0)*XIH 
      Q2(61,29) = 3**(0.5D0)*7**(0.5D0)*ABS(XIH) 
      Q2(61,45) = -5**(0.5D0)*7**(0.5D0)*XIH 
      Q2(61,49) = -7**(0.5D0)*ETAH 
      Q2(61,53) = 3**(0.5D0)*7**(0.5D0)*ABS(ETAH) 
      Q2(61,57) = -5**(0.5D0)*7**(0.5D0)*ETAH 
      Q2(61,61) = VT + 7*ABS(ETAH) + ABS(MUH) + 7*ABS(XIH) 
      Q2(61,62) = 3**(0.5D0)*MUH 
      Q2(61,63) = 5**(0.5D0)*ABS(MUH) 
      Q2(61,64) = 7**(0.5D0)*MUH 
      Q2(62,14) = -7**(0.5D0)*XIH 
      Q2(62,30) = 3**(0.5D0)*7**(0.5D0)*ABS(XIH) 
      Q2(62,46) = -5**(0.5D0)*7**(0.5D0)*XIH 
      Q2(62,50) = -7**(0.5D0)*ETAH 
      Q2(62,54) = 3**(0.5D0)*7**(0.5D0)*ABS(ETAH) 
      Q2(62,58) = -5**(0.5D0)*7**(0.5D0)*ETAH 
      Q2(62,61) = -3**(0.5D0)*MUH 
      Q2(62,62) = VT + 7*ABS(ETAH) + 3*ABS(MUH) + 7*ABS(XIH) 
      Q2(62,63) = 2*15**(0.5D0)*MUH - 3**(0.5D0)*5**(0.5D0)*MUH 
      Q2(62,64) = 3**(0.5D0)*7**(0.5D0)*ABS(MUH) 
      Q2(63,15) = -7**(0.5D0)*XIH 
      Q2(63,31) = 3**(0.5D0)*7**(0.5D0)*ABS(XIH) 
      Q2(63,47) = -5**(0.5D0)*7**(0.5D0)*XIH 
      Q2(63,51) = -7**(0.5D0)*ETAH 
      Q2(63,55) = 3**(0.5D0)*7**(0.5D0)*ABS(ETAH) 
      Q2(63,59) = -5**(0.5D0)*7**(0.5D0)*ETAH 
      Q2(63,61) = 5**(0.5D0)*ABS(MUH) 
      Q2(63,62) = -3**(0.5D0)*5**(0.5D0)*MUH 
      Q2(63,63) = VT + 7*ABS(ETAH) + 5*ABS(MUH) + 7*ABS(XIH) 
      Q2(63,64) = 2*35**(0.5D0)*MUH - 5**(0.5D0)*7**(0.5D0)*MUH 
      Q2(64,16) = -7**(0.5D0)*XIH 
      Q2(64,32) = 3**(0.5D0)*7**(0.5D0)*ABS(XIH) 
      Q2(64,48) = -5**(0.5D0)*7**(0.5D0)*XIH 
      Q2(64,52) = -7**(0.5D0)*ETAH 
      Q2(64,56) = 3**(0.5D0)*7**(0.5D0)*ABS(ETAH) 
      Q2(64,60) = -5**(0.5D0)*7**(0.5D0)*ETAH 
      Q2(64,61) = -7**(0.5D0)*MUH 
      Q2(64,62) = 3**(0.5D0)*7**(0.5D0)*ABS(MUH) 
      Q2(64,63) = -5**(0.5D0)*7**(0.5D0)*MUH 
      Q2(64,64) = VT + 7*ABS(ETAH) + 7*ABS(MUH) + 7*ABS(XIH) 
      Q2(01,65) = Q(01)*VOL(I2,J2,J,I,K2) 
      Q2(02,65) = Q(02)*VOL(I2,J2,J,I,K2) 
      Q2(03,65) = Q(03)*VOL(I2,J2,J,I,K2) 
      Q2(04,65) = Q(04)*VOL(I2,J2,J,I,K2) 
      Q2(05,65) = Q(05)*VOL(I2,J2,J,I,K2) 
      Q2(06,65) = Q(06)*VOL(I2,J2,J,I,K2) 
      Q2(07,65) = Q(07)*VOL(I2,J2,J,I,K2) 
      Q2(08,65) = Q(08)*VOL(I2,J2,J,I,K2) 
      Q2(09,65) = Q(09)*VOL(I2,J2,J,I,K2) 
      Q2(10,65) = Q(10)*VOL(I2,J2,J,I,K2) 
      Q2(11,65) = Q(11)*VOL(I2,J2,J,I,K2) 
      Q2(12,65) = Q(12)*VOL(I2,J2,J,I,K2) 
      Q2(13,65) = Q(13)*VOL(I2,J2,J,I,K2) 
      Q2(14,65) = Q(14)*VOL(I2,J2,J,I,K2) 
      Q2(15,65) = Q(15)*VOL(I2,J2,J,I,K2) 
      Q2(16,65) = Q(16)*VOL(I2,J2,J,I,K2) 
      Q2(17,65) = Q(17)*VOL(I2,J2,J,I,K2) 
      Q2(18,65) = Q(18)*VOL(I2,J2,J,I,K2) 
      Q2(19,65) = Q(19)*VOL(I2,J2,J,I,K2) 
      Q2(20,65) = Q(20)*VOL(I2,J2,J,I,K2) 
      Q2(21,65) = Q(21)*VOL(I2,J2,J,I,K2) 
      Q2(22,65) = Q(22)*VOL(I2,J2,J,I,K2) 
      Q2(23,65) = Q(23)*VOL(I2,J2,J,I,K2) 
      Q2(24,65) = Q(24)*VOL(I2,J2,J,I,K2) 
      Q2(25,65) = Q(25)*VOL(I2,J2,J,I,K2) 
      Q2(26,65) = Q(26)*VOL(I2,J2,J,I,K2) 
      Q2(27,65) = Q(27)*VOL(I2,J2,J,I,K2) 
      Q2(28,65) = Q(28)*VOL(I2,J2,J,I,K2) 
      Q2(29,65) = Q(29)*VOL(I2,J2,J,I,K2) 
      Q2(30,65) = Q(30)*VOL(I2,J2,J,I,K2) 
      Q2(31,65) = Q(31)*VOL(I2,J2,J,I,K2) 
      Q2(32,65) = Q(32)*VOL(I2,J2,J,I,K2) 
      Q2(33,65) = Q(33)*VOL(I2,J2,J,I,K2) 
      Q2(34,65) = Q(34)*VOL(I2,J2,J,I,K2) 
      Q2(35,65) = Q(35)*VOL(I2,J2,J,I,K2) 
      Q2(36,65) = Q(36)*VOL(I2,J2,J,I,K2) 
      Q2(37,65) = Q(37)*VOL(I2,J2,J,I,K2) 
      Q2(38,65) = Q(38)*VOL(I2,J2,J,I,K2) 
      Q2(39,65) = Q(39)*VOL(I2,J2,J,I,K2) 
      Q2(40,65) = Q(40)*VOL(I2,J2,J,I,K2) 
      Q2(41,65) = Q(41)*VOL(I2,J2,J,I,K2) 
      Q2(42,65) = Q(42)*VOL(I2,J2,J,I,K2) 
      Q2(43,65) = Q(43)*VOL(I2,J2,J,I,K2) 
      Q2(44,65) = Q(44)*VOL(I2,J2,J,I,K2) 
      Q2(45,65) = Q(45)*VOL(I2,J2,J,I,K2) 
      Q2(46,65) = Q(46)*VOL(I2,J2,J,I,K2) 
      Q2(47,65) = Q(47)*VOL(I2,J2,J,I,K2) 
      Q2(48,65) = Q(48)*VOL(I2,J2,J,I,K2) 
      Q2(49,65) = Q(49)*VOL(I2,J2,J,I,K2) 
      Q2(50,65) = Q(50)*VOL(I2,J2,J,I,K2) 
      Q2(51,65) = Q(51)*VOL(I2,J2,J,I,K2) 
      Q2(52,65) = Q(52)*VOL(I2,J2,J,I,K2) 
      Q2(53,65) = Q(53)*VOL(I2,J2,J,I,K2) 
      Q2(54,65) = Q(54)*VOL(I2,J2,J,I,K2) 
      Q2(55,65) = Q(55)*VOL(I2,J2,J,I,K2) 
      Q2(56,65) = Q(56)*VOL(I2,J2,J,I,K2) 
      Q2(57,65) = Q(57)*VOL(I2,J2,J,I,K2) 
      Q2(58,65) = Q(58)*VOL(I2,J2,J,I,K2) 
      Q2(59,65) = Q(59)*VOL(I2,J2,J,I,K2) 
      Q2(60,65) = Q(60)*VOL(I2,J2,J,I,K2) 
      Q2(61,65) = Q(61)*VOL(I2,J2,J,I,K2) 
      Q2(62,65) = Q(62)*VOL(I2,J2,J,I,K2) 
      Q2(63,65) = Q(63)*VOL(I2,J2,J,I,K2) 
      Q2(64,65) = Q(64)*VOL(I2,J2,J,I,K2) 

      Q2(01,65) = Q2(01,65) + MUH*XNI(1,1,J2,K2)*SIGN(1.0,MUH
     >   ) + XIH*XNK(1,1)*SIGN(1.0,XIH) + ETAH*XNJ(1,1,K2)
     >   *SIGN(1.0,ETAH) 
      Q2(02,65) = Q2(02,65) + XIH*XNK(2,1)*SIGN(1.0,XIH) + 
     >   ETAH*XNJ(2,1,K2)*SIGN(1.0,ETAH) - 3**(0.5D0)*MUH*X
     >   NI(1,1,J2,K2) 
      Q2(03,65) = Q2(03,65) + XIH*XNK(3,1)*SIGN(1.0,XIH) + 
     >   ETAH*XNJ(3,1,K2)*SIGN(1.0,ETAH) + 5**(0.5D0)*MUH*X
     >   NI(1,1,J2,K2)*SIGN(1.0,MUH) 
      Q2(04,65) = Q2(04,65) + XIH*XNK(4,1)*SIGN(1.0,XIH) + 
     >   ETAH*XNJ(4,1,K2)*SIGN(1.0,ETAH) - 7**(0.5D0)*MUH*X
     >   NI(1,1,J2,K2) 
      Q2(05,65) = Q2(05,65) + MUH*XNI(2,1,J2,K2)*SIGN(1.0,MUH
     >   ) + XIH*XNK(1,2)*SIGN(1.0,XIH) - 3**(0.5D0)*ETAH
     >   *XNJ(1,1,K2) 
      Q2(06,65) = Q2(06,65) + XIH*XNK(2,2)*SIGN(1.0,XIH) - 
     >   3**(0.5D0)*ETAH*XNJ(2,1,K2) - 3**(0.5D0)*MUH*XNI(2,
     >   1,J2,K2) 
      Q2(07,65) = Q2(07,65) + XIH*XNK(3,2)*SIGN(1.0,XIH) - 
     >   3**(0.5D0)*ETAH*XNJ(3,1,K2) + 5**(0.5D0)*MUH*XNI(2,
     >   1,J2,K2)*SIGN(1.0,MUH) 
      Q2(08,65) = Q2(08,65) + XIH*XNK(4,2)*SIGN(1.0,XIH) - 
     >   3**(0.5D0)*ETAH*XNJ(4,1,K2) - 7**(0.5D0)*MUH*XNI(2,
     >   1,J2,K2) 
      Q2(09,65) = Q2(09,65) + MUH*XNI(3,1,J2,K2)*SIGN(1.0,MUH
     >   ) + XIH*XNK(1,3)*SIGN(1.0,XIH) + 5**(0.5D0)*ETAH
     >   *XNJ(1,1,K2)*SIGN(1.0,ETAH) 
      Q2(10,65) = Q2(10,65) + XIH*XNK(2,3)*SIGN(1.0,XIH) + 
     >   5**(0.5D0)*ETAH*XNJ(2,1,K2)*SIGN(1.0,ETAH) - 3**(0.5D0)*
     >   MUH*XNI(3,1,J2,K2) 
      Q2(11,65) = Q2(11,65) + XIH*XNK(3,3)*SIGN(1.0,XIH) + 
     >   5**(0.5D0)*ETAH*XNJ(3,1,K2)*SIGN(1.0,ETAH) + 5**(0.5D0)*
     >   MUH*XNI(3,1,J2,K2)*SIGN(1.0,MUH) 
      Q2(12,65) = Q2(12,65) + XIH*XNK(4,3)*SIGN(1.0,XIH) + 
     >   5**(0.5D0)*ETAH*XNJ(4,1,K2)*SIGN(1.0,ETAH) - 7**(0.5D0)*
     >   MUH*XNI(3,1,J2,K2) 
      Q2(13,65) = Q2(13,65) + MUH*XNI(4,1,J2,K2)*SIGN(1.0,MUH
     >   ) + XIH*XNK(1,4)*SIGN(1.0,XIH) - 7**(0.5D0)*ETAH
     >   *XNJ(1,1,K2) 
      Q2(14,65) = Q2(14,65) + XIH*XNK(2,4)*SIGN(1.0,XIH) - 
     >   7**(0.5D0)*ETAH*XNJ(2,1,K2) - 3**(0.5D0)*MUH*XNI(4,
     >   1,J2,K2) 
      Q2(15,65) = Q2(15,65) + XIH*XNK(3,4)*SIGN(1.0,XIH) - 
     >   7**(0.5D0)*ETAH*XNJ(3,1,K2) + 5**(0.5D0)*MUH*XNI(4,
     >   1,J2,K2)*SIGN(1.0,MUH) 
      Q2(16,65) = Q2(16,65) + XIH*XNK(4,4)*SIGN(1.0,XIH) - 
     >   7**(0.5D0)*ETAH*XNJ(4,1,K2) - 7**(0.5D0)*MUH*XNI(4,
     >   1,J2,K2) 
      Q2(17,65) = Q2(17,65) + MUH*XNI(1,2,J2,K2)*SIGN(1.0,MUH
     >   ) + ETAH*XNJ(1,2,K2)*SIGN(1.0,ETAH) - 3**(0.5D0)*
     >   XIH*XNK(1,1) 
      Q2(18,65) = Q2(18,65) + ETAH*XNJ(2,2,K2)*SIGN(1.0,ETAH) 
     >   - 3**(0.5D0)*MUH*XNI(1,2,J2,K2) - 3**(0.5D0)*XIH*XN
     >   K(2,1) 
      Q2(19,65) = Q2(19,65) + ETAH*XNJ(3,2,K2)*SIGN(1.0,ETAH) 
     >   + 5**(0.5D0)*MUH*XNI(1,2,J2,K2)*SIGN(1.0,MUH) - 3**(0.5
     >   D0)*XIH*XNK(3,1) 
      Q2(20,65) = Q2(20,65) + ETAH*XNJ(4,2,K2)*SIGN(1.0,ETAH) 
     >   - 7**(0.5D0)*MUH*XNI(1,2,J2,K2) - 3**(0.5D0)*XIH*XN
     >   K(4,1) 
      Q2(21,65) = Q2(21,65) + MUH*XNI(2,2,J2,K2)*SIGN(1.0,MUH
     >   ) - 3**(0.5D0)*ETAH*XNJ(1,2,K2) - 3**(0.5D0)*XIH*XN
     >   K(1,2) 
      Q2(22,65) = Q2(22,65) - 3**(0.5D0)*ETAH*XNJ(2,2,K2) - 3**
     >   (0.5D0)*MUH*XNI(2,2,J2,K2) - 3**(0.5D0)*XIH*XNK(2,2
     >   ) 
      Q2(23,65) = Q2(23,65) + 5**(0.5D0)*MUH*XNI(2,2,J2,K2)*SIG
     >   N(1.0,MUH) - 3**(0.5D0)*ETAH*XNJ(3,2,K2) - 3**(0.5D0)*
     >   XIH*XNK(3,2) 
      Q2(24,65) = Q2(24,65) - 3**(0.5D0)*ETAH*XNJ(4,2,K2) - 7**
     >   (0.5D0)*MUH*XNI(2,2,J2,K2) - 3**(0.5D0)*XIH*XNK(4,2
     >   ) 
      Q2(25,65) = Q2(25,65) + MUH*XNI(3,2,J2,K2)*SIGN(1.0,MUH
     >   ) + 5**(0.5D0)*ETAH*XNJ(1,2,K2)*SIGN(1.0,ETAH) - 3**(0.5
     >   D0)*XIH*XNK(1,3) 
      Q2(26,65) = Q2(26,65) + 5**(0.5D0)*ETAH*XNJ(2,2,K2)*SIGN(
     >   1.0,ETAH) - 3**(0.5D0)*MUH*XNI(3,2,J2,K2) - 3**(0.5D0)*
     >   XIH*XNK(2,3) 
      Q2(27,65) = Q2(27,65) + 5**(0.5D0)*ETAH*XNJ(3,2,K2)*SIGN(
     >   1.0,ETAH) + 5**(0.5D0)*MUH*XNI(3,2,J2,K2)*SIGN(1.0,MUH
     >   ) - 3**(0.5D0)*XIH*XNK(3,3) 
      Q2(28,65) = Q2(28,65) + 5**(0.5D0)*ETAH*XNJ(4,2,K2)*SIGN(
     >   1.0,ETAH) - 7**(0.5D0)*MUH*XNI(3,2,J2,K2) - 3**(0.5D0)*
     >   XIH*XNK(4,3) 
      Q2(29,65) = Q2(29,65) + MUH*XNI(4,2,J2,K2)*SIGN(1.0,MUH
     >   ) - 7**(0.5D0)*ETAH*XNJ(1,2,K2) - 3**(0.5D0)*XIH*XN
     >   K(1,4) 
      Q2(30,65) = Q2(30,65) - 7**(0.5D0)*ETAH*XNJ(2,2,K2) - 3**
     >   (0.5D0)*MUH*XNI(4,2,J2,K2) - 3**(0.5D0)*XIH*XNK(2,4
     >   ) 
      Q2(31,65) = Q2(31,65) + 5**(0.5D0)*MUH*XNI(4,2,J2,K2)*SIG
     >   N(1.0,MUH) - 7**(0.5D0)*ETAH*XNJ(3,2,K2) - 3**(0.5D0)*
     >   XIH*XNK(3,4) 
      Q2(32,65) = Q2(32,65) - 7**(0.5D0)*ETAH*XNJ(4,2,K2) - 7**
     >   (0.5D0)*MUH*XNI(4,2,J2,K2) - 3**(0.5D0)*XIH*XNK(4,4
     >   ) 
      Q2(33,65) = Q2(33,65) + MUH*XNI(1,3,J2,K2)*SIGN(1.0,MUH
     >   ) + ETAH*XNJ(1,3,K2)*SIGN(1.0,ETAH) + 5**(0.5D0)*
     >   XIH*XNK(1,1)*SIGN(1.0,XIH) 
      Q2(34,65) = Q2(34,65) + ETAH*XNJ(2,3,K2)*SIGN(1.0,ETAH) 
     >   - 3**(0.5D0)*MUH*XNI(1,3,J2,K2) + 5**(0.5D0)*XIH*XN
     >   K(2,1)*SIGN(1.0,XIH) 
      Q2(35,65) = Q2(35,65) + ETAH*XNJ(3,3,K2)*SIGN(1.0,ETAH) 
     >   + 5**(0.5D0)*MUH*XNI(1,3,J2,K2)*SIGN(1.0,MUH) + 5**(0.5
     >   D0)*XIH*XNK(3,1)*SIGN(1.0,XIH) 
      Q2(36,65) = Q2(36,65) + ETAH*XNJ(4,3,K2)*SIGN(1.0,ETAH) 
     >   - 7**(0.5D0)*MUH*XNI(1,3,J2,K2) + 5**(0.5D0)*XIH*XN
     >   K(4,1)*SIGN(1.0,XIH) 
      Q2(37,65) = Q2(37,65) + MUH*XNI(2,3,J2,K2)*SIGN(1.0,MUH
     >   ) - 3**(0.5D0)*ETAH*XNJ(1,3,K2) + 5**(0.5D0)*XIH*XN
     >   K(1,2)*SIGN(1.0,XIH) 
      Q2(38,65) = Q2(38,65) + 5**(0.5D0)*XIH*XNK(2,2)*SIGN(1.
     >   0,XIH) - 3**(0.5D0)*MUH*XNI(2,3,J2,K2) - 3**(0.5D0)*
     >   ETAH*XNJ(2,3,K2) 
      Q2(39,65) = Q2(39,65) + 5**(0.5D0)*MUH*XNI(2,3,J2,K2)*SIG
     >   N(1.0,MUH) - 3**(0.5D0)*ETAH*XNJ(3,3,K2) + 5**(0.5D0)*
     >   XIH*XNK(3,2)*SIGN(1.0,XIH) 
      Q2(40,65) = Q2(40,65) + 5**(0.5D0)*XIH*XNK(4,2)*SIGN(1.
     >   0,XIH) - 7**(0.5D0)*MUH*XNI(2,3,J2,K2) - 3**(0.5D0)*
     >   ETAH*XNJ(4,3,K2) 
      Q2(41,65) = Q2(41,65) + MUH*XNI(3,3,J2,K2)*SIGN(1.0,MUH
     >   ) + 5**(0.5D0)*ETAH*XNJ(1,3,K2)*SIGN(1.0,ETAH) + 5**(0.5
     >   D0)*XIH*XNK(1,3)*SIGN(1.0,XIH) 
      Q2(42,65) = Q2(42,65) + 5**(0.5D0)*ETAH*XNJ(2,3,K2)*SIGN(
     >   1.0,ETAH) - 3**(0.5D0)*MUH*XNI(3,3,J2,K2) + 5**(0.5D0)*
     >   XIH*XNK(2,3)*SIGN(1.0,XIH) 
      Q2(43,65) = Q2(43,65) + 5**(0.5D0)*ETAH*XNJ(3,3,K2)*SIGN(
     >   1.0,ETAH) + 5**(0.5D0)*MUH*XNI(3,3,J2,K2)*SIGN(1.0,MUH
     >   ) + 5**(0.5D0)*XIH*XNK(3,3)*SIGN(1.0,XIH) 
      Q2(44,65) = Q2(44,65) + 5**(0.5D0)*ETAH*XNJ(4,3,K2)*SIGN(
     >   1.0,ETAH) - 7**(0.5D0)*MUH*XNI(3,3,J2,K2) + 5**(0.5D0)*
     >   XIH*XNK(4,3)*SIGN(1.0,XIH) 
      Q2(45,65) = Q2(45,65) + MUH*XNI(4,3,J2,K2)*SIGN(1.0,MUH
     >   ) - 7**(0.5D0)*ETAH*XNJ(1,3,K2) + 5**(0.5D0)*XIH*XN
     >   K(1,4)*SIGN(1.0,XIH) 
      Q2(46,65) = Q2(46,65) + 5**(0.5D0)*XIH*XNK(2,4)*SIGN(1.
     >   0,XIH) - 3**(0.5D0)*MUH*XNI(4,3,J2,K2) - 7**(0.5D0)*
     >   ETAH*XNJ(2,3,K2) 
      Q2(47,65) = Q2(47,65) + 5**(0.5D0)*MUH*XNI(4,3,J2,K2)*SIG
     >   N(1.0,MUH) - 7**(0.5D0)*ETAH*XNJ(3,3,K2) + 5**(0.5D0)*
     >   XIH*XNK(3,4)*SIGN(1.0,XIH) 
      Q2(48,65) = Q2(48,65) + 5**(0.5D0)*XIH*XNK(4,4)*SIGN(1.
     >   0,XIH) - 7**(0.5D0)*MUH*XNI(4,3,J2,K2) - 7**(0.5D0)*
     >   ETAH*XNJ(4,3,K2) 
      Q2(49,65) = Q2(49,65) + MUH*XNI(1,4,J2,K2)*SIGN(1.0,MUH
     >   ) + ETAH*XNJ(1,4,K2)*SIGN(1.0,ETAH) - 7**(0.5D0)*
     >   XIH*XNK(1,1) 
      Q2(50,65) = Q2(50,65) + ETAH*XNJ(2,4,K2)*SIGN(1.0,ETAH) 
     >   - 3**(0.5D0)*MUH*XNI(1,4,J2,K2) - 7**(0.5D0)*XIH*XN
     >   K(2,1) 
      Q2(51,65) = Q2(51,65) + ETAH*XNJ(3,4,K2)*SIGN(1.0,ETAH) 
     >   + 5**(0.5D0)*MUH*XNI(1,4,J2,K2)*SIGN(1.0,MUH) - 7**(0.5
     >   D0)*XIH*XNK(3,1) 
      Q2(52,65) = Q2(52,65) + ETAH*XNJ(4,4,K2)*SIGN(1.0,ETAH) 
     >   - 7**(0.5D0)*MUH*XNI(1,4,J2,K2) - 7**(0.5D0)*XIH*XN
     >   K(4,1) 
      Q2(53,65) = Q2(53,65) + MUH*XNI(2,4,J2,K2)*SIGN(1.0,MUH
     >   ) - 3**(0.5D0)*ETAH*XNJ(1,4,K2) - 7**(0.5D0)*XIH*XN
     >   K(1,2) 
      Q2(54,65) = Q2(54,65) - 3**(0.5D0)*ETAH*XNJ(2,4,K2) - 3**
     >   (0.5D0)*MUH*XNI(2,4,J2,K2) - 7**(0.5D0)*XIH*XNK(2,2
     >   ) 
      Q2(55,65) = Q2(55,65) + 5**(0.5D0)*MUH*XNI(2,4,J2,K2)*SIG
     >   N(1.0,MUH) - 3**(0.5D0)*ETAH*XNJ(3,4,K2) - 7**(0.5D0)*
     >   XIH*XNK(3,2) 
      Q2(56,65) = Q2(56,65) - 3**(0.5D0)*ETAH*XNJ(4,4,K2) - 7**
     >   (0.5D0)*MUH*XNI(2,4,J2,K2) - 7**(0.5D0)*XIH*XNK(4,2
     >   ) 
      Q2(57,65) = Q2(57,65) + MUH*XNI(3,4,J2,K2)*SIGN(1.0,MUH
     >   ) + 5**(0.5D0)*ETAH*XNJ(1,4,K2)*SIGN(1.0,ETAH) - 7**(0.5
     >   D0)*XIH*XNK(1,3) 
      Q2(58,65) = Q2(58,65) + 5**(0.5D0)*ETAH*XNJ(2,4,K2)*SIGN(
     >   1.0,ETAH) - 3**(0.5D0)*MUH*XNI(3,4,J2,K2) - 7**(0.5D0)*
     >   XIH*XNK(2,3) 
      Q2(59,65) = Q2(59,65) + 5**(0.5D0)*ETAH*XNJ(3,4,K2)*SIGN(
     >   1.0,ETAH) + 5**(0.5D0)*MUH*XNI(3,4,J2,K2)*SIGN(1.0,MUH
     >   ) - 7**(0.5D0)*XIH*XNK(3,3) 
      Q2(60,65) = Q2(60,65) + 5**(0.5D0)*ETAH*XNJ(4,4,K2)*SIGN(
     >   1.0,ETAH) - 7**(0.5D0)*MUH*XNI(3,4,J2,K2) - 7**(0.5D0)*
     >   XIH*XNK(4,3) 
      Q2(61,65) = Q2(61,65) + MUH*XNI(4,4,J2,K2)*SIGN(1.0,MUH
     >   ) - 7**(0.5D0)*ETAH*XNJ(1,4,K2) - 7**(0.5D0)*XIH*XN
     >   K(1,4) 
      Q2(62,65) = Q2(62,65) - 7**(0.5D0)*ETAH*XNJ(2,4,K2) - 3**
     >   (0.5D0)*MUH*XNI(4,4,J2,K2) - 7**(0.5D0)*XIH*XNK(2,4
     >   ) 
      Q2(63,65) = Q2(63,65) + 5**(0.5D0)*MUH*XNI(4,4,J2,K2)*SIG
     >   N(1.0,MUH) - 7**(0.5D0)*ETAH*XNJ(3,4,K2) - 7**(0.5D0)*
     >   XIH*XNK(3,4) 
      Q2(64,65) = Q2(64,65) - 7**(0.5D0)*ETAH*XNJ(4,4,K2) - 7**
     >   (0.5D0)*MUH*XNI(4,4,J2,K2) - 7**(0.5D0)*XIH*XNK(4,4
     >   ) 
      ENDIF
*
      CALL ALSBD(IELEM**3,1,Q2,IER,IELEM**3)
      IF(IER.NE.0) CALL XABORT('SNFTH4: SINGULAR MATRIX.')
*
      IF(IELEM.EQ.1)THEN
      XNK(1,1)     = Q2(1,2) 
      ELSEIF(IELEM.EQ.2)THEN
      XNK(1,1)     = Q2(1,9) + SIGN(1.0,XIH)*3**(0.5D0)*Q2(5,9) 
      XNK(2,1)     = Q2(2,9) + SIGN(1.0,XIH)*3**(0.5D0)*Q2(6,9) 
      XNK(1,2)     = Q2(3,9) + SIGN(1.0,XIH)*3**(0.5D0)*Q2(7,9) 
      XNK(2,2)     = Q2(4,9) + SIGN(1.0,XIH)*3**(0.5D0)*Q2(8,9) 
      ELSEIF(IELEM.EQ.3)THEN
      XNK(1,1)     = Q2(01,28) + SIGN(1.0,XIH)*3**(0.5D0)*Q2(10,2
     >   8) + 5**(0.5D0)*Q2(19,28)
      XNK(2,1)     = Q2(02,28) + SIGN(1.0,XIH)*3**(0.5D0)*Q2(11,2
     >   8) + 5**(0.5D0)*Q2(20,28)
      XNK(3,1)     = Q2(03,28) + SIGN(1.0,XIH)*3**(0.5D0)*Q2(12,2
     >   8) + 5**(0.5D0)*Q2(21,28)
      XNK(1,2)     = Q2(04,28) + SIGN(1.0,XIH)*3**(0.5D0)*Q2(13,2
     >   8) + 5**(0.5D0)*Q2(22,28)
      XNK(2,2)     = Q2(05,28) + SIGN(1.0,XIH)*3**(0.5D0)*Q2(14,2
     >   8) + 5**(0.5D0)*Q2(23,28)
      XNK(3,2)     = Q2(06,28) + SIGN(1.0,XIH)*3**(0.5D0)*Q2(15,2
     >   8) + 5**(0.5D0)*Q2(24,28)
      XNK(1,3)     = Q2(07,28) + SIGN(1.0,XIH)*3**(0.5D0)*Q2(16,2
     >   8) + 5**(0.5D0)*Q2(25,28)
      XNK(2,3)     = Q2(08,28) + SIGN(1.0,XIH)*3**(0.5D0)*Q2(17,2
     >   8) + 5**(0.5D0)*Q2(26,28)
      XNK(3,3)     = Q2(09,28) + SIGN(1.0,XIH)*3**(0.5D0)*Q2(18,2
     >   8) + 5**(0.5D0)*Q2(27,28)
      ELSEIF(IELEM.EQ.4)THEN
      XNK(1,1)     = Q2(01,65) + SIGN(1.0,XIH)*3**(0.5D0)*Q2(17,6
     >   5) + 5**(0.5D0)*Q2(33,65) + SIGN(1.0,XIH)*7**(0.5D0)*Q2(49,
     >   65)
      XNK(2,1)     = Q2(02,65) + SIGN(1.0,XIH)*3**(0.5D0)*Q2(18,6
     >   5) + 5**(0.5D0)*Q2(34,65) + SIGN(1.0,XIH)*7**(0.5D0)*Q2(50,
     >   65)
      XNK(3,1)     = Q2(03,65) + SIGN(1.0,XIH)*3**(0.5D0)*Q2(19,6
     >   5) + 5**(0.5D0)*Q2(35,65) + SIGN(1.0,XIH)*7**(0.5D0)*Q2(51,
     >   65)
      XNK(4,1)     = Q2(04,65) + SIGN(1.0,XIH)*3**(0.5D0)*Q2(20,6
     >   5) + 5**(0.5D0)*Q2(36,65) + SIGN(1.0,XIH)*7**(0.5D0)*Q2(52,
     >   65)
      XNK(1,2)     = Q2(05,65) + SIGN(1.0,XIH)*3**(0.5D0)*Q2(21,6
     >   5) + 5**(0.5D0)*Q2(37,65) + SIGN(1.0,XIH)*7**(0.5D0)*Q2(53,
     >   65)
      XNK(2,2)     = Q2(06,65) + SIGN(1.0,XIH)*3**(0.5D0)*Q2(22,6
     >   5) + 5**(0.5D0)*Q2(38,65) + SIGN(1.0,XIH)*7**(0.5D0)*Q2(54,
     >   65)
      XNK(3,2)     = Q2(07,65) + SIGN(1.0,XIH)*3**(0.5D0)*Q2(23,6
     >   5) + 5**(0.5D0)*Q2(39,65) + SIGN(1.0,XIH)*7**(0.5D0)*Q2(55,
     >   65)
      XNK(4,2)     = Q2(08,65) + SIGN(1.0,XIH)*3**(0.5D0)*Q2(24,6
     >   5) + 5**(0.5D0)*Q2(40,65) + SIGN(1.0,XIH)*7**(0.5D0)*Q2(56,
     >   65)
      XNK(1,3)     = Q2(09,65) + SIGN(1.0,XIH)*3**(0.5D0)*Q2(25,6
     >   5) + 5**(0.5D0)*Q2(41,65) + SIGN(1.0,XIH)*7**(0.5D0)*Q2(57,
     >   65)
      XNK(2,3)     = Q2(10,65) + SIGN(1.0,XIH)*3**(0.5D0)*Q2(26,6
     >   5) + 5**(0.5D0)*Q2(42,65) + SIGN(1.0,XIH)*7**(0.5D0)*Q2(58,
     >   65)
      XNK(3,3)     = Q2(11,65) + SIGN(1.0,XIH)*3**(0.5D0)*Q2(27,6
     >   5) + 5**(0.5D0)*Q2(43,65) + SIGN(1.0,XIH)*7**(0.5D0)*Q2(59,
     >   65)
      XNK(4,3)     = Q2(12,65) + SIGN(1.0,XIH)*3**(0.5D0)*Q2(28,6
     >   5) + 5**(0.5D0)*Q2(44,65) + SIGN(1.0,XIH)*7**(0.5D0)*Q2(60,
     >   65)
      XNK(1,4)     = Q2(13,65) + SIGN(1.0,XIH)*3**(0.5D0)*Q2(29,6
     >   5) + 5**(0.5D0)*Q2(45,65) + SIGN(1.0,XIH)*7**(0.5D0)*Q2(61,
     >   65)
      XNK(2,4)     = Q2(14,65) + SIGN(1.0,XIH)*3**(0.5D0)*Q2(30,6
     >   5) + 5**(0.5D0)*Q2(46,65) + SIGN(1.0,XIH)*7**(0.5D0)*Q2(62,
     >   65)
      XNK(3,4)     = Q2(15,65) + SIGN(1.0,XIH)*3**(0.5D0)*Q2(31,6
     >   5) + 5**(0.5D0)*Q2(47,65) + SIGN(1.0,XIH)*7**(0.5D0)*Q2(63,
     >   65)
      XNK(4,4)     = Q2(16,65) + SIGN(1.0,XIH)*3**(0.5D0)*Q2(32,6
     >   5) + 5**(0.5D0)*Q2(48,65) + SIGN(1.0,XIH)*7**(0.5D0)*Q2(64,
     >   65)
      ENDIF
*
      IF(JL.LT.ISPLH)THEN
      IF(IELEM.EQ.1)THEN
      XNJ(1,1,K2)   = Q2(1,2) 
      ELSEIF(IELEM.EQ.2)THEN
      XNJ(1,1,K2)   = Q2(1,9) + SIGN(1.0,ETAH)*3**(0.5D0)*Q2(3,9) 
      XNJ(2,1,K2)   = Q2(2,9) + SIGN(1.0,ETAH)*3**(0.5D0)*Q2(4,9) 
      XNJ(1,2,K2)   = Q2(5,9) + SIGN(1.0,ETAH)*3**(0.5D0)*Q2(7,9) 
      XNJ(2,2,K2)   = Q2(6,9) + SIGN(1.0,ETAH)*3**(0.5D0)*Q2(8,9) 
      ELSEIF(IELEM.EQ.3)THEN
      XNJ(1,1,K2)   = Q2(01,28) + SIGN(1.0,ETAH)*3**(0.5D0)*Q2(04,2
     >   8) + 5**(0.5D0)*Q2(07,28)
      XNJ(2,1,K2)   = Q2(02,28) + SIGN(1.0,ETAH)*3**(0.5D0)*Q2(05,2
     >   8) + 5**(0.5D0)*Q2(08,28)
      XNJ(3,1,K2)   = Q2(03,28) + SIGN(1.0,ETAH)*3**(0.5D0)*Q2(06,2
     >   8) + 5**(0.5D0)*Q2(09,28)
      XNJ(1,2,K2)   = Q2(10,28) + SIGN(1.0,ETAH)*3**(0.5D0)*Q2(13,2
     >   8) + 5**(0.5D0)*Q2(16,28)
      XNJ(2,2,K2)   = Q2(11,28) + SIGN(1.0,ETAH)*3**(0.5D0)*Q2(14,2
     >   8) + 5**(0.5D0)*Q2(17,28)
      XNJ(3,2,K2)   = Q2(12,28) + SIGN(1.0,ETAH)*3**(0.5D0)*Q2(15,2
     >   8) + 5**(0.5D0)*Q2(18,28)
      XNJ(1,3,K2)   = Q2(19,28) + SIGN(1.0,ETAH)*3**(0.5D0)*Q2(22,2
     >   8) + 5**(0.5D0)*Q2(25,28)
      XNJ(2,3,K2)   = Q2(20,28) + SIGN(1.0,ETAH)*3**(0.5D0)*Q2(23,2
     >   8) + 5**(0.5D0)*Q2(26,28)
      XNJ(3,3,K2)   = Q2(21,28) + SIGN(1.0,ETAH)*3**(0.5D0)*Q2(24,2
     >   8) + 5**(0.5D0)*Q2(27,28)
      ELSEIF(IELEM.EQ.4)THEN
      XNJ(1,1,K2)   = Q2(01,65) + SIGN(1.0,ETAH)*3**(0.5D0)*Q2(05,6
     >   5) + 5**(0.5D0)*Q2(09,65) + SIGN(1.0,ETAH)*7**(0.5D0)*Q2(13,
     >   65)
      XNJ(2,1,K2)   = Q2(02,65) + SIGN(1.0,ETAH)*3**(0.5D0)*Q2(06,6
     >   5) + 5**(0.5D0)*Q2(10,65) + SIGN(1.0,ETAH)*7**(0.5D0)*Q2(14,
     >   65)
      XNJ(3,1,K2)   = Q2(03,65) + SIGN(1.0,ETAH)*3**(0.5D0)*Q2(07,6
     >   5) + 5**(0.5D0)*Q2(11,65) + SIGN(1.0,ETAH)*7**(0.5D0)*Q2(15,
     >   65)
      XNJ(4,1,K2)   = Q2(04,65) + SIGN(1.0,ETAH)*3**(0.5D0)*Q2(08,6
     >   5) + 5**(0.5D0)*Q2(12,65) + SIGN(1.0,ETAH)*7**(0.5D0)*Q2(16,
     >   65)
      XNJ(1,2,K2)   = Q2(17,65) + SIGN(1.0,ETAH)*3**(0.5D0)*Q2(21,6
     >   5) + 5**(0.5D0)*Q2(25,65) + SIGN(1.0,ETAH)*7**(0.5D0)*Q2(29,
     >   65)
      XNJ(2,2,K2)   = Q2(18,65) + SIGN(1.0,ETAH)*3**(0.5D0)*Q2(22,6
     >   5) + 5**(0.5D0)*Q2(26,65) + SIGN(1.0,ETAH)*7**(0.5D0)*Q2(30,
     >   65)
      XNJ(3,2,K2)   = Q2(19,65) + SIGN(1.0,ETAH)*3**(0.5D0)*Q2(23,6
     >   5) + 5**(0.5D0)*Q2(27,65) + SIGN(1.0,ETAH)*7**(0.5D0)*Q2(31,
     >   65)
      XNJ(4,2,K2)   = Q2(20,65) + SIGN(1.0,ETAH)*3**(0.5D0)*Q2(24,6
     >   5) + 5**(0.5D0)*Q2(28,65) + SIGN(1.0,ETAH)*7**(0.5D0)*Q2(32,
     >   65)
      XNJ(1,3,K2)   = Q2(33,65) + SIGN(1.0,ETAH)*3**(0.5D0)*Q2(37,6
     >   5) + 5**(0.5D0)*Q2(41,65) + SIGN(1.0,ETAH)*7**(0.5D0)*Q2(45,
     >   65)
      XNJ(2,3,K2)   = Q2(34,65) + SIGN(1.0,ETAH)*3**(0.5D0)*Q2(38,6
     >   5) + 5**(0.5D0)*Q2(42,65) + SIGN(1.0,ETAH)*7**(0.5D0)*Q2(46,
     >   65)
      XNJ(3,3,K2)   = Q2(35,65) + SIGN(1.0,ETAH)*3**(0.5D0)*Q2(39,6
     >   5) + 5**(0.5D0)*Q2(43,65) + SIGN(1.0,ETAH)*7**(0.5D0)*Q2(47,
     >   65)
      XNJ(4,3,K2)   = Q2(36,65) + SIGN(1.0,ETAH)*3**(0.5D0)*Q2(40,6
     >   5) + 5**(0.5D0)*Q2(44,65) + SIGN(1.0,ETAH)*7**(0.5D0)*Q2(48,
     >   65)
      XNJ(1,4,K2)   = Q2(49,65) + SIGN(1.0,ETAH)*3**(0.5D0)*Q2(53,6
     >   5) + 5**(0.5D0)*Q2(57,65) + SIGN(1.0,ETAH)*7**(0.5D0)*Q2(61,
     >   65)
      XNJ(2,4,K2)   = Q2(50,65) + SIGN(1.0,ETAH)*3**(0.5D0)*Q2(54,6
     >   5) + 5**(0.5D0)*Q2(58,65) + SIGN(1.0,ETAH)*7**(0.5D0)*Q2(62,
     >   65)
      XNJ(3,4,K2)   = Q2(51,65) + SIGN(1.0,ETAH)*3**(0.5D0)*Q2(55,6
     >   5) + 5**(0.5D0)*Q2(59,65) + SIGN(1.0,ETAH)*7**(0.5D0)*Q2(63,
     >   65)
      XNJ(4,4,K2)   = Q2(52,65) + SIGN(1.0,ETAH)*3**(0.5D0)*Q2(56,6
     >   5) + 5**(0.5D0)*Q2(60,65) + SIGN(1.0,ETAH)*7**(0.5D0)*Q2(64,
     >   65)
      ENDIF
      ELSEIF((JL.EQ.ISPLH).AND.(IHEXJ.LE.NHEX))THEN
         I3=I2
         C1=1.0D0
         IF((J.EQ.1).AND.(ILOZJ.EQ.3)) THEN
            I3=ISPLH+1 -I2
            C1=-1.0D0
         ENDIF
      IF(IELEM.EQ.1)THEN
      BFLUX(ISIDEJ,INDEXJ,I3,1,1,K2) = Q2(1,2) 
      ELSEIF(IELEM.EQ.2)THEN
      BFLUX(ISIDEJ,INDEXJ,I3,1,1,K2) = Q2(1,9) + SIGN(1.0,ETAH)
     >   *3**(0.5D0)*Q2(3,9) 
      BFLUX(ISIDEJ,INDEXJ,I3,2,1,K2) = (Q2(2,9) + SIGN(1.0,ETAH)
     >   *3**(0.5D0)*Q2(4,9))*C1 
      BFLUX(ISIDEJ,INDEXJ,I3,1,2,K2) = (Q2(5,9) + SIGN(1.0,ETAH)
     >   *3**(0.5D0)*Q2(7,9))*C1 
      BFLUX(ISIDEJ,INDEXJ,I3,2,2,K2) = Q2(6,9) + SIGN(1.0,ETAH)
     >   *3**(0.5D0)*Q2(8,9) 
      ELSEIF(IELEM.EQ.3)THEN
      BFLUX(ISIDEJ,INDEXJ,I3,1,1,K2) = Q2(01,28) + SIGN(1.0,ETAH)*
     >   3**(0.5D0)*Q2(04,28) + 5**(0.5D0)*Q2(07,28)
      BFLUX(ISIDEJ,INDEXJ,I3,2,1,K2) = (Q2(02,28) + SIGN(1.0,ETAH)*
     >   3**(0.5D0)*Q2(05,28) + 5**(0.5D0)*Q2(08,28))*C1
      BFLUX(ISIDEJ,INDEXJ,I3,3,1,K2) = Q2(03,28) + SIGN(1.0,ETAH)*
     >   3**(0.5D0)*Q2(06,28) + 5**(0.5D0)*Q2(09,28)
      BFLUX(ISIDEJ,INDEXJ,I3,1,2,K2) = (Q2(10,28) + SIGN(1.0,ETAH)*
     >   3**(0.5D0)*Q2(13,28) + 5**(0.5D0)*Q2(16,28))*C1
      BFLUX(ISIDEJ,INDEXJ,I3,2,2,K2) = Q2(11,28) + SIGN(1.0,ETAH)*
     >   3**(0.5D0)*Q2(14,28) + 5**(0.5D0)*Q2(17,28)
      BFLUX(ISIDEJ,INDEXJ,I3,3,2,K2) = (Q2(12,28) + SIGN(1.0,ETAH)*
     >   3**(0.5D0)*Q2(15,28) + 5**(0.5D0)*Q2(18,28))*C1
      BFLUX(ISIDEJ,INDEXJ,I3,1,3,K2) = Q2(19,28) + SIGN(1.0,ETAH)*
     >   3**(0.5D0)*Q2(22,28) + 5**(0.5D0)*Q2(25,28)
      BFLUX(ISIDEJ,INDEXJ,I3,2,3,K2) = (Q2(20,28) + SIGN(1.0,ETAH)*
     >   3**(0.5D0)*Q2(23,28) + 5**(0.5D0)*Q2(26,28))*C1
      BFLUX(ISIDEJ,INDEXJ,I3,3,3,K2) = Q2(21,28) + SIGN(1.0,ETAH)*
     >   3**(0.5D0)*Q2(24,28) + 5**(0.5D0)*Q2(27,28)
      ELSEIF(IELEM.EQ.4)THEN
      BFLUX(ISIDEJ,INDEXJ,I3,1,1,K2) = Q2(01,65) + SIGN(1.0,ETAH)*
     >   3**(0.5D0)*Q2(05,65) + 5**(0.5D0)*Q2(09,65) + 
     >   SIGN(1.0,ETAH)*7**(0.5D0)*Q2(13,65)
      BFLUX(ISIDEJ,INDEXJ,I3,2,1,K2) = (Q2(02,65) + SIGN(1.0,ETAH)*
     >   3**(0.5D0)*Q2(06,65) + 5**(0.5D0)*Q2(10,65) + 
     >   SIGN(1.0,ETAH)*7**(0.5D0)*Q2(14,65))*C1
      BFLUX(ISIDEJ,INDEXJ,I3,3,1,K2) = Q2(03,65) + SIGN(1.0,ETAH)*
     >   3**(0.5D0)*Q2(07,65) + 5**(0.5D0)*Q2(11,65) + 
     >   SIGN(1.0,ETAH)*7**(0.5D0)*Q2(15,65)
      BFLUX(ISIDEJ,INDEXJ,I3,4,1,K2) = (Q2(04,65) + SIGN(1.0,ETAH)*
     >   3**(0.5D0)*Q2(08,65) + 5**(0.5D0)*Q2(12,65) + 
     >   SIGN(1.0,ETAH)*7**(0.5D0)*Q2(16,65))*C1
      BFLUX(ISIDEJ,INDEXJ,I3,1,2,K2) = (Q2(17,65) + SIGN(1.0,ETAH)*
     >   3**(0.5D0)*Q2(21,65) + 5**(0.5D0)*Q2(25,65) + 
     >   SIGN(1.0,ETAH)*7**(0.5D0)*Q2(29,65))*C1
      BFLUX(ISIDEJ,INDEXJ,I3,2,2,K2) = Q2(18,65) + SIGN(1.0,ETAH)*
     >   3**(0.5D0)*Q2(22,65) + 5**(0.5D0)*Q2(26,65) + 
     >   SIGN(1.0,ETAH)*7**(0.5D0)*Q2(30,65)
      BFLUX(ISIDEJ,INDEXJ,I3,3,2,K2) = (Q2(19,65) + SIGN(1.0,ETAH)*
     >   3**(0.5D0)*Q2(23,65) + 5**(0.5D0)*Q2(27,65) + 
     >   SIGN(1.0,ETAH)*7**(0.5D0)*Q2(31,65))*C1
      BFLUX(ISIDEJ,INDEXJ,I3,4,2,K2) = Q2(20,65) + SIGN(1.0,ETAH)*
     >   3**(0.5D0)*Q2(24,65) + 5**(0.5D0)*Q2(28,65) + 
     >   SIGN(1.0,ETAH)*7**(0.5D0)*Q2(32,65)
      BFLUX(ISIDEJ,INDEXJ,I3,1,3,K2) = Q2(33,65) + SIGN(1.0,ETAH)*
     >   3**(0.5D0)*Q2(37,65) + 5**(0.5D0)*Q2(41,65) + 
     >   SIGN(1.0,ETAH)*7**(0.5D0)*Q2(45,65)
      BFLUX(ISIDEJ,INDEXJ,I3,2,3,K2) = (Q2(34,65) + SIGN(1.0,ETAH)*
     >   3**(0.5D0)*Q2(38,65) + 5**(0.5D0)*Q2(42,65) + 
     >   SIGN(1.0,ETAH)*7**(0.5D0)*Q2(46,65))*C1
      BFLUX(ISIDEJ,INDEXJ,I3,3,3,K2) = Q2(35,65) + SIGN(1.0,ETAH)*
     >   3**(0.5D0)*Q2(39,65) + 5**(0.5D0)*Q2(43,65) + 
     >   SIGN(1.0,ETAH)*7**(0.5D0)*Q2(47,65)
      BFLUX(ISIDEJ,INDEXJ,I3,4,3,K2) = (Q2(36,65) + SIGN(1.0,ETAH)*
     >   3**(0.5D0)*Q2(40,65) + 5**(0.5D0)*Q2(44,65) + 
     >   SIGN(1.0,ETAH)*7**(0.5D0)*Q2(48,65))*C1
      BFLUX(ISIDEJ,INDEXJ,I3,1,4,K2) = (Q2(49,65) + SIGN(1.0,ETAH)*
     >   3**(0.5D0)*Q2(53,65) + 5**(0.5D0)*Q2(57,65) + 
     >   SIGN(1.0,ETAH)*7**(0.5D0)*Q2(61,65))*C1
      BFLUX(ISIDEJ,INDEXJ,I3,2,4,K2) = Q2(50,65) + SIGN(1.0,ETAH)*
     >   3**(0.5D0)*Q2(54,65) + 5**(0.5D0)*Q2(58,65) + 
     >   SIGN(1.0,ETAH)*7**(0.5D0)*Q2(62,65)
      BFLUX(ISIDEJ,INDEXJ,I3,3,4,K2) = (Q2(51,65) + SIGN(1.0,ETAH)*
     >   3**(0.5D0)*Q2(55,65) + 5**(0.5D0)*Q2(59,65) + 
     >   SIGN(1.0,ETAH)*7**(0.5D0)*Q2(63,65))*C1
      BFLUX(ISIDEJ,INDEXJ,I3,4,4,K2) = Q2(52,65) + SIGN(1.0,ETAH)*
     >   3**(0.5D0)*Q2(56,65) + 5**(0.5D0)*Q2(60,65) + 
     >   SIGN(1.0,ETAH)*7**(0.5D0)*Q2(64,65)
      ENDIF
      ENDIF
*
      IF(IL.LT.ISPLH)THEN
      IF(IELEM.EQ.1)THEN
      XNI(1,1,J2,K2) = Q2(1,2) 
      ELSEIF(IELEM.EQ.2)THEN
      XNI(1,1,J2,K2) = Q2(1,9) + SIGN(1.0,MUH)*3**(0.5D0)*Q2(2,9) 
      XNI(2,1,J2,K2) = Q2(3,9) + SIGN(1.0,MUH)*3**(0.5D0)*Q2(4,9) 
      XNI(1,2,J2,K2) = Q2(5,9) + SIGN(1.0,MUH)*3**(0.5D0)*Q2(6,9) 
      XNI(2,2,J2,K2) = Q2(7,9) + SIGN(1.0,MUH)*3**(0.5D0)*Q2(8,9) 
      ELSEIF(IELEM.EQ.3)THEN
      XNI(1,1,J2,K2) = Q2(01,28) + SIGN(1.0,MUH)*3**(0.5D0)*Q2(02,2
     >   8) + 5**(0.5D0)*Q2(03,28)
      XNI(2,1,J2,K2) = Q2(04,28) + SIGN(1.0,MUH)*3**(0.5D0)*Q2(05,2
     >   8) + 5**(0.5D0)*Q2(06,28)
      XNI(3,1,J2,K2) = Q2(07,28) + SIGN(1.0,MUH)*3**(0.5D0)*Q2(08,2
     >   8) + 5**(0.5D0)*Q2(09,28)
      XNI(1,2,J2,K2) = Q2(10,28) + SIGN(1.0,MUH)*3**(0.5D0)*Q2(11,2
     >   8) + 5**(0.5D0)*Q2(12,28)
      XNI(2,2,J2,K2) = Q2(13,28) + SIGN(1.0,MUH)*3**(0.5D0)*Q2(14,2
     >   8) + 5**(0.5D0)*Q2(15,28)
      XNI(3,2,J2,K2) = Q2(16,28) + SIGN(1.0,MUH)*3**(0.5D0)*Q2(17,2
     >   8) + 5**(0.5D0)*Q2(18,28)
      XNI(1,3,J2,K2) = Q2(19,28) + SIGN(1.0,MUH)*3**(0.5D0)*Q2(20,2
     >   8) + 5**(0.5D0)*Q2(21,28)
      XNI(2,3,J2,K2) = Q2(22,28) + SIGN(1.0,MUH)*3**(0.5D0)*Q2(23,2
     >   8) + 5**(0.5D0)*Q2(24,28)
      XNI(3,3,J2,K2) = Q2(25,28) + SIGN(1.0,MUH)*3**(0.5D0)*Q2(26,2
     >   8) + 5**(0.5D0)*Q2(27,28)
      ELSEIF(IELEM.EQ.4)THEN
      XNI(1,1,J2,K2) = Q2(01,65) + SIGN(1.0,MUH)*3**(0.5D0)*Q2(02,6
     >   5) + 5**(0.5D0)*Q2(03,65) + SIGN(1.0,MUH)*7**(0.5D0)*Q2(04,
     >   65)
      XNI(2,1,J2,K2) = Q2(05,65) + SIGN(1.0,MUH)*3**(0.5D0)*Q2(06,6
     >   5) + 5**(0.5D0)*Q2(07,65) + SIGN(1.0,MUH)*7**(0.5D0)*Q2(08,
     >   65)
      XNI(3,1,J2,K2) = Q2(09,65) + SIGN(1.0,MUH)*3**(0.5D0)*Q2(10,6
     >   5) + 5**(0.5D0)*Q2(11,65) + SIGN(1.0,MUH)*7**(0.5D0)*Q2(12,
     >   65)
      XNI(4,1,J2,K2) = Q2(13,65) + SIGN(1.0,MUH)*3**(0.5D0)*Q2(14,6
     >   5) + 5**(0.5D0)*Q2(15,65) + SIGN(1.0,MUH)*7**(0.5D0)*Q2(16,
     >   65)
      XNI(1,2,J2,K2) = Q2(17,65) + SIGN(1.0,MUH)*3**(0.5D0)*Q2(18,6
     >   5) + 5**(0.5D0)*Q2(19,65) + SIGN(1.0,MUH)*7**(0.5D0)*Q2(20,
     >   65)
      XNI(2,2,J2,K2) = Q2(21,65) + SIGN(1.0,MUH)*3**(0.5D0)*Q2(22,6
     >   5) + 5**(0.5D0)*Q2(23,65) + SIGN(1.0,MUH)*7**(0.5D0)*Q2(24,
     >   65)
      XNI(3,2,J2,K2) = Q2(25,65) + SIGN(1.0,MUH)*3**(0.5D0)*Q2(26,6
     >   5) + 5**(0.5D0)*Q2(27,65) + SIGN(1.0,MUH)*7**(0.5D0)*Q2(28,
     >   65)
      XNI(4,2,J2,K2) = Q2(29,65) + SIGN(1.0,MUH)*3**(0.5D0)*Q2(30,6
     >   5) + 5**(0.5D0)*Q2(31,65) + SIGN(1.0,MUH)*7**(0.5D0)*Q2(32,
     >   65)
      XNI(1,3,J2,K2) = Q2(33,65) + SIGN(1.0,MUH)*3**(0.5D0)*Q2(34,6
     >   5) + 5**(0.5D0)*Q2(35,65) + SIGN(1.0,MUH)*7**(0.5D0)*Q2(36,
     >   65)
      XNI(2,3,J2,K2) = Q2(37,65) + SIGN(1.0,MUH)*3**(0.5D0)*Q2(38,6
     >   5) + 5**(0.5D0)*Q2(39,65) + SIGN(1.0,MUH)*7**(0.5D0)*Q2(40,
     >   65)
      XNI(3,3,J2,K2) = Q2(41,65) + SIGN(1.0,MUH)*3**(0.5D0)*Q2(42,6
     >   5) + 5**(0.5D0)*Q2(43,65) + SIGN(1.0,MUH)*7**(0.5D0)*Q2(44,
     >   65)
      XNI(4,3,J2,K2) = Q2(45,65) + SIGN(1.0,MUH)*3**(0.5D0)*Q2(46,6
     >   5) + 5**(0.5D0)*Q2(47,65) + SIGN(1.0,MUH)*7**(0.5D0)*Q2(48,
     >   65)
      XNI(1,4,J2,K2) = Q2(49,65) + SIGN(1.0,MUH)*3**(0.5D0)*Q2(50,6
     >   5) + 5**(0.5D0)*Q2(51,65) + SIGN(1.0,MUH)*7**(0.5D0)*Q2(52,
     >   65)
      XNI(2,4,J2,K2) = Q2(53,65) + SIGN(1.0,MUH)*3**(0.5D0)*Q2(54,6
     >   5) + 5**(0.5D0)*Q2(55,65) + SIGN(1.0,MUH)*7**(0.5D0)*Q2(56,
     >   65)
      XNI(3,4,J2,K2) = Q2(57,65) + SIGN(1.0,MUH)*3**(0.5D0)*Q2(58,6
     >   5) + 5**(0.5D0)*Q2(59,65) + SIGN(1.0,MUH)*7**(0.5D0)*Q2(60,
     >   65)
      XNI(4,4,J2,K2) = Q2(61,65) + SIGN(1.0,MUH)*3**(0.5D0)*Q2(62,6
     >   5) + 5**(0.5D0)*Q2(63,65) + SIGN(1.0,MUH)*7**(0.5D0)*Q2(64,
     >   65)
      ENDIF
      ELSEIF((IL.EQ.ISPLH).AND.(IHEXI.LE.NHEX))THEN
         J3=J2
         C1=1.0D0
         IF((J.EQ.3).AND.(ILOZI.EQ.1)) THEN
            J3=ISPLH+1-J2
            C1=-1.0D0
         ENDIF
      IF(IELEM.EQ.1)THEN
      BFLUX(ISIDEI,INDEXI,J3,1,1,K2) = Q2(1,2) 
      ELSEIF(IELEM.EQ.2)THEN
      BFLUX(ISIDEI,INDEXI,J3,1,1,K2) = Q2(1,9) + SIGN(1.0,MUH)
     >   *3**(0.5D0)*Q2(2,9) 
      BFLUX(ISIDEI,INDEXI,J3,2,1,K2) = (Q2(3,9) + SIGN(1.0,MUH)
     >   *3**(0.5D0)*Q2(4,9))*C1
      BFLUX(ISIDEI,INDEXI,J3,1,2,K2) = (Q2(5,9) + SIGN(1.0,MUH)
     >   *3**(0.5D0)*Q2(6,9))*C1
      BFLUX(ISIDEI,INDEXI,J3,2,2,K2) = Q2(7,9) + SIGN(1.0,MUH)
     >   *3**(0.5D0)*Q2(8,9) 
      ELSEIF(IELEM.EQ.3)THEN
      BFLUX(ISIDEI,INDEXI,J3,1,1,K2) = Q2(01,28) + SIGN(1.0,MUH)
     >   *3**(0.5D0)*Q2(02,28) + 5**(0.5D0)*Q2(03,28)
      BFLUX(ISIDEI,INDEXI,J3,2,1,K2) = (Q2(04,28) + SIGN(1.0,MUH)
     >   *3**(0.5D0)*Q2(05,28) + 5**(0.5D0)*Q2(06,28))*C1
      BFLUX(ISIDEI,INDEXI,J3,3,1,K2) = Q2(07,28) + SIGN(1.0,MUH)
     >   *3**(0.5D0)*Q2(08,28) + 5**(0.5D0)*Q2(09,28)
      BFLUX(ISIDEI,INDEXI,J3,1,2,K2) = (Q2(10,28) + SIGN(1.0,MUH)
     >   *3**(0.5D0)*Q2(11,28) + 5**(0.5D0)*Q2(12,28))*C1
      BFLUX(ISIDEI,INDEXI,J3,2,2,K2) = Q2(13,28) + SIGN(1.0,MUH)
     >   *3**(0.5D0)*Q2(14,28) + 5**(0.5D0)*Q2(15,28)
      BFLUX(ISIDEI,INDEXI,J3,3,2,K2) = (Q2(16,28) + SIGN(1.0,MUH)
     >   *3**(0.5D0)*Q2(17,28) + 5**(0.5D0)*Q2(18,28))*C1
      BFLUX(ISIDEI,INDEXI,J3,1,3,K2) = Q2(19,28) + SIGN(1.0,MUH)
     >   *3**(0.5D0)*Q2(20,28) + 5**(0.5D0)*Q2(21,28)
      BFLUX(ISIDEI,INDEXI,J3,2,3,K2) = (Q2(22,28) + SIGN(1.0,MUH)
     >   *3**(0.5D0)*Q2(23,28) + 5**(0.5D0)*Q2(24,28))*C1
      BFLUX(ISIDEI,INDEXI,J3,3,3,K2) = Q2(25,28) + SIGN(1.0,MUH)
     >   *3**(0.5D0)*Q2(26,28) + 5**(0.5D0)*Q2(27,28)
      ELSEIF(IELEM.EQ.4)THEN
      BFLUX(ISIDEI,INDEXI,J3,1,1,K2) = Q2(01,65) + SIGN(1.0,MUH)
     >   *3**(0.5D0)*Q2(02,65) + 5**(0.5D0)*Q2(03,65) + 
     >   SIGN(1.0,MUH)*7**(0.5D0)*Q2(04,65)
      BFLUX(ISIDEI,INDEXI,J3,2,1,K2) = (Q2(05,65) + SIGN(1.0,MUH)
     >   *3**(0.5D0)*Q2(06,65) + 5**(0.5D0)*Q2(07,65) + 
     >   SIGN(1.0,MUH)*7**(0.5D0)*Q2(08,65))*C1
      BFLUX(ISIDEI,INDEXI,J3,3,1,K2) = Q2(09,65) + SIGN(1.0,MUH)
     >   *3**(0.5D0)*Q2(10,65) + 5**(0.5D0)*Q2(11,65) + 
     >   SIGN(1.0,MUH)*7**(0.5D0)*Q2(12,65)
      BFLUX(ISIDEI,INDEXI,J3,4,1,K2) = (Q2(13,65) + SIGN(1.0,MUH)
     >   *3**(0.5D0)*Q2(14,65) + 5**(0.5D0)*Q2(15,65) + 
     >   SIGN(1.0,MUH)*7**(0.5D0)*Q2(16,65))*C1
      BFLUX(ISIDEI,INDEXI,J3,1,2,K2) = (Q2(17,65) + SIGN(1.0,MUH)
     >   *3**(0.5D0)*Q2(18,65) + 5**(0.5D0)*Q2(19,65) + 
     >   SIGN(1.0,MUH)*7**(0.5D0)*Q2(20,65))*C1
      BFLUX(ISIDEI,INDEXI,J3,2,2,K2) = Q2(21,65) + SIGN(1.0,MUH)
     >   *3**(0.5D0)*Q2(22,65) + 5**(0.5D0)*Q2(23,65) + 
     >   SIGN(1.0,MUH)*7**(0.5D0)*Q2(24,65)
      BFLUX(ISIDEI,INDEXI,J3,3,2,K2) = (Q2(25,65) + SIGN(1.0,MUH)
     >   *3**(0.5D0)*Q2(26,65) + 5**(0.5D0)*Q2(27,65) + 
     >   SIGN(1.0,MUH)*7**(0.5D0)*Q2(28,65))*C1
      BFLUX(ISIDEI,INDEXI,J3,4,2,K2) = Q2(29,65) + SIGN(1.0,MUH)
     >   *3**(0.5D0)*Q2(30,65) + 5**(0.5D0)*Q2(31,65) + 
     >   SIGN(1.0,MUH)*7**(0.5D0)*Q2(32,65)
      BFLUX(ISIDEI,INDEXI,J3,1,3,K2) = Q2(33,65) + SIGN(1.0,MUH)
     >   *3**(0.5D0)*Q2(34,65) + 5**(0.5D0)*Q2(35,65) + 
     >   SIGN(1.0,MUH)*7**(0.5D0)*Q2(36,65)
      BFLUX(ISIDEI,INDEXI,J3,2,3,K2) = (Q2(37,65) + SIGN(1.0,MUH)
     >   *3**(0.5D0)*Q2(38,65) + 5**(0.5D0)*Q2(39,65) + 
     >   SIGN(1.0,MUH)*7**(0.5D0)*Q2(40,65))*C1
      BFLUX(ISIDEI,INDEXI,J3,3,3,K2) = Q2(41,65) + SIGN(1.0,MUH)
     >   *3**(0.5D0)*Q2(42,65) + 5**(0.5D0)*Q2(43,65) + 
     >   SIGN(1.0,MUH)*7**(0.5D0)*Q2(44,65)
      BFLUX(ISIDEI,INDEXI,J3,4,3,K2) = (Q2(45,65) + SIGN(1.0,MUH)
     >   *3**(0.5D0)*Q2(46,65) + 5**(0.5D0)*Q2(47,65) + 
     >   SIGN(1.0,MUH)*7**(0.5D0)*Q2(48,65))*C1
      BFLUX(ISIDEI,INDEXI,J3,1,4,K2) = (Q2(49,65) + SIGN(1.0,MUH)
     >   *3**(0.5D0)*Q2(50,65) + 5**(0.5D0)*Q2(51,65) + 
     >   SIGN(1.0,MUH)*7**(0.5D0)*Q2(52,65))*C1
      BFLUX(ISIDEI,INDEXI,J3,2,4,K2) = Q2(53,65) + SIGN(1.0,MUH)
     >   *3**(0.5D0)*Q2(54,65) + 5**(0.5D0)*Q2(55,65) + 
     >   SIGN(1.0,MUH)*7**(0.5D0)*Q2(56,65)
      BFLUX(ISIDEI,INDEXI,J3,3,4,K2) = (Q2(57,65) + SIGN(1.0,MUH)
     >   *3**(0.5D0)*Q2(58,65) + 5**(0.5D0)*Q2(59,65) + 
     >   SIGN(1.0,MUH)*7**(0.5D0)*Q2(60,65))*C1
      BFLUX(ISIDEI,INDEXI,J3,4,4,K2) = Q2(61,65) + SIGN(1.0,MUH)
     >   *3**(0.5D0)*Q2(62,65) + 5**(0.5D0)*Q2(63,65) + 
     >   SIGN(1.0,MUH)*7**(0.5D0)*Q2(64,65)
      ENDIF
      ENDIF
*
      ILOZ=((J-1)*ISPLH+(J2-1))*ISPLH+I2
      DO 210 K=1,NSCT
      DO 200 IEL=1,IELEM**3
      FLUX(IEL,K,ILOZ,I,K2) = FLUX(IEL,K,ILOZ,I,K2) +
     1 W(M)*Q2(IEL,IELEM**3+1)*PL(K,M)
  200 CONTINUE
  210 CONTINUE
*
  220 CONTINUE
      DO 240 IEL=1,IELEM
      DO 230 JEL=1,IELEM
      IOF=((((((M-1)*NHEX+(I-1))*3+(J-1))*ISPLH+(J2-1))*ISPLH+
     1    (I2-1))*IELEM+(JEL-1))*IELEM+IEL
      FUNKNO(L4+IOF,IG)=REAL(XNK(IEL,JEL))
  230 CONTINUE
  240 CONTINUE
*
  250 CONTINUE
  260 CONTINUE
  270 CONTINUE
  280 CONTINUE
      FLUX_G(:,:,:,:,:,IG)=FLUX_G(:,:,:,:,:,IG)+FLUX(:,:,:,:,:)
  290 CONTINUE
  300 CONTINUE
*$OMP END PARALLEL DO
  310 CONTINUE
      DO 320 IG=1,NGEFF
        IF(.NOT.INCONV(IG)) GO TO 320
        FUNKNO(:L4,IG)=
     1  RESHAPE(REAL(FLUX_G(:IELEM**3,:NSCT,:3*ISPLH**2,:NHEX,:LZ,IG)),
     2  (/L4/))
  320 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(FLUX_G,FLUX,BFLUX,INDANG)
      RETURN
  500 FORMAT(16H SNFTH4: thread=,I8,12H --->(group=,I4,7H angle=,I4,1H))
      END
