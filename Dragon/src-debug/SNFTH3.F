*DECK SNFTH3
      SUBROUTINE SNFTH3(NUN,NGRP,NGEFF,IMPX,INCONV,INGIND,NHEX,LZ,IELEM,
     1 ISPLH,SIDE,NMAT,NPQ,NSCT,MAT,VOL,NCODE,ZCODE,TOTAL,QEXT,DU,DE,DZ,
     2 W,MRMZ,DC,DB,DA,PL,CONNEC,IZGLOB,CONFROM,FUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform one inner iteration for solving SN equations in 3D hexagonal
* geometry for the HODD method. VOID boundary conditions on sides, and,
* VOID or REFL boundary conditions on top and bottom. Energy-angle
* multithreading. Albedo boundary conditions.
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
* NGRP    number of energy groups.
* NGEFF   dimension of arrays INCONV and INGIND.
* IMPX    print flag (equal to zero for no print).
* INCONV  energy group convergence flag (set to .FALSE. if converged).
* INGIND  energy group index assign to 1:NGEFF arrays.
* NHEX    number of hexagons in X-Y plane.
* LZ      number of meshes along Z axis.
* IELEM   measure of order of the spatial approximation polynomial:
*         =1 constant (classical diamond scheme - default for HODD);
*         =2 linear (default for DG);
*         =3 parabolic;
*         =4 cubic (only for DG).
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
      INTEGER NUN,NGRP,NGEFF,IMPX,INGIND(NGEFF),NHEX,LZ,IELEM,ISPLH,
     1 NMAT,NPQ,NSCT,IZGLOB(NHEX,6),MAT(ISPLH,ISPLH,3,NHEX,LZ),
     2 NCODE(6),MRMZ(NPQ),CONNEC(3,(NHEX*3)*2,6),CONFROM(2,3,6)
      LOGICAL INCONV(NGEFF)
      REAL SIDE,VOL(ISPLH,ISPLH,3,NHEX,LZ),ZCODE(6),TOTAL(0:NMAT,NGRP),
     1 QEXT(NUN,NGRP),DU(NPQ),DE(NPQ),DZ(NPQ),W(NPQ),
     2 DC(ISPLH*ISPLH*3*NHEX,1,NPQ),DB(ISPLH*ISPLH*3*NHEX,LZ,NPQ),
     3 DA(1,LZ,NPQ),PL(NSCT,NPQ),FUNKNO(NUN,NGRP)
*----
*  LOCAL VARIABLES
*----
      INTEGER NPQD(12),IIND(12)
      DOUBLE PRECISION Q(IELEM**3), Q2(IELEM**3,(IELEM**3)+1),
     1 VT, CONST0, CONST1, CONST2
      PARAMETER(IUNOUT=6,PI=3.141592654)
      INTEGER :: ILOZSWP(3,6), IFROMI, IFROMJ
      REAL :: JAC(2,2,3), MUH, ETAH, XIH, AAA, BBB, CCC, DDD, MUHTEMP,
     1 ETAHTEMP
      DOUBLE PRECISION :: THETA, XNI(IELEM,IELEM,ISPLH,LZ),
     1 XNJ(IELEM,IELEM,LZ),XNK(IELEM,IELEM),C1
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
      CONST0=2.0D0*DSQRT(3.0D0)
      CONST1=2.0D0*DSQRT(5.0D0)
      CONST2=2.0D0*DSQRT(15.0D0)
*----
*  PARAMETER VALIDATION.
*----
      IF((IELEM.LT.1).OR.(IELEM.GT.3))
     1   CALL XABORT('SNFTH3: INVALID IELEM (DIAM) VALUE. '
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
      DO 150 IGEFF=1,NGEFF
      DO 140 IPQD=1,NPQD(IND)
      IF(.NOT.INCONV(IGEFF)) GO TO 140
      M=INDANG(IPQD,IND)
      IG=INGIND(IGEFF)
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
*$OMP4 PRIVATE(MUHTEMP,ETAHTEMP,MUH,ETAH,XIH,I2,I3,J2,J3,K2,C1)
*$OMP5 SHARED(FUNKNO) REDUCTION(+:FLUX_G) COLLAPSE(2)
      DO 300 IGEFF=1,NGEFF
      DO 290 IPQD=1,NPQD(IND)
#if defined(_OPENMP)
        ITID=omp_get_thread_num()
#else
        ITID=0
#endif
      IF(IMPX.GT.5) WRITE(IUNOUT,500) ITID,IGEFF,IPQD
      IF(.NOT.INCONV(IGEFF)) GO TO 290
      M=INDANG(IPQD,IND)
      IG=INGIND(IGEFF)
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

*     -----------------------------------------------------
      DO 170 IEL=1,IELEM**3
      Q(IEL)=0.0D0
      DO 160 K=1,NSCT
      IOF=(((((((K2-1)*NHEX+(I-1))*3+(J-1))*ISPLH+(J2-1))*ISPLH+
     1 (I2-1))*NSCT+(K-1))*IELEM**3)+IEL
      Q(IEL)=Q(IEL)+QEXT(IOF,IG)*PL(K,M)/(4.0D0*PI)
  160 CONTINUE
  170 CONTINUE
*     -----------------------------------------------------
*
      VT=VOL(I2,J2,J,I,K2)*TOTAL(MAT(I2,J2,J,I,K2),IG)
      CALL XDDSET(Q2,(IELEM**3)*((IELEM**3)+1),0.0D0)
*
*     -----------------------------------------------------
      IF(IELEM.EQ.1) THEN
      Q2(1,1)=2.0D0*ABS(MUH)+2.0D0*ABS(ETAH)+
     1        2.0D0*ABS(XIH)+VT
      Q2(1,2)=2.0D0*ABS(MUH)*XNI(1,1,J2,K2)+
     1        2.0D0*ABS(ETAH)*XNJ(1,1,K2)+
     2        2.0D0*ABS(XIH)*XNK(1,1)+VOL(I2,J2,J,I,K2)*Q(1)
      ELSE IF(IELEM.EQ.2) THEN
      Q2(1,1)=VT
      Q2(1,2)=CONST0*MUH
      Q2(1,3)=CONST0*ETAH
      Q2(1,4)=CONST0*XIH
      Q2(2,2)=-VT-6.0D0*ABS(MUH)
      Q2(2,5)=-CONST0*ETAH
      Q2(2,6)=-CONST0*XIH
      Q2(3,3)=-VT-6.0D0*ABS(ETAH)
      Q2(3,5)=-CONST0*MUH
      Q2(3,7)=-CONST0*XIH
      Q2(4,4)=-VT-6.0D0*ABS(XIH)
      Q2(4,6)=-CONST0*MUH
      Q2(4,7)=-CONST0*ETAH
      Q2(5,5)=VT+6.0D0*ABS(MUH)+6.0D0*ABS(ETAH)
      Q2(5,8)=CONST0*XIH
      Q2(6,6)=VT+6.0D0*ABS(MUH)+6.0D0*ABS(XIH)
      Q2(6,8)=CONST0*ETAH
      Q2(7,7)=VT+6.0D0*ABS(ETAH)+6.0D0*ABS(XIH)
      Q2(7,8)=CONST0*MUH
      Q2(8,8)=-VT-6.0D0*ABS(MUH)-6.0D0*ABS(ETAH)
     1        -6.0D0*ABS(XIH)
*---------------
      Q2(1,9)=VOL(I2,J2,J,I,K2)*Q(1)
      Q2(2,9)=-VOL(I2,J2,J,I,K2)*Q(2)+CONST0*MUH*XNI(1,1,J2,K2)
      Q2(3,9)=-VOL(I2,J2,J,I,K2)*Q(3)+CONST0*ETAH*XNJ(1,1,K2)
      Q2(4,9)=-VOL(I2,J2,J,I,K2)*Q(4)+CONST0*XIH*XNK(1,1)
      Q2(5,9)=VOL(I2,J2,J,I,K2)*Q(5)-CONST0*MUH*XNI(2,1,J2,K2)
     1       -CONST0*ETAH*XNJ(2,1,K2)
      Q2(6,9)=VOL(I2,J2,J,I,K2)*Q(6)-CONST0*MUH*XNI(1,2,J2,K2)
     1       -CONST0*XIH*XNK(2,1)
      Q2(7,9)=VOL(I2,J2,J,I,K2)*Q(7)-CONST0*ETAH*XNJ(1,2,K2)
     1       -CONST0*XIH*XNK(1,2)
      Q2(8,9)=-VOL(I2,J2,J,I,K2)*Q(8)+CONST0*MUH*XNI(2,2,J2,K2)
     1       +CONST0*ETAH*XNJ(2,2,K2)
     2       +CONST0*XIH*XNK(2,2)
      ELSE IF(IELEM.EQ.3) THEN
*--------------- UPPER DIAGONAL TERMS
      Q2(1,1)=-VT-2.0D0*ABS(MUH)-2.0D0*ABS(ETAH)-
     1        2.0D0*ABS(XIH)
      Q2(1,9)=-CONST1*ABS(MUH)
      Q2(1,10)=-CONST1*ABS(ETAH)
      Q2(1,11)=-CONST1*ABS(XIH)
      Q2(2,2)=VT+2.0D0*ABS(ETAH)+2.0D0*ABS(XIH)
      Q2(2,9)=CONST2*MUH
      Q2(2,16)=CONST1*ABS(XIH)
      Q2(2,17)=CONST1*ABS(ETAH)
      Q2(3,3)=VT+2.0D0*ABS(MUH)+2.0D0*ABS(XIH)
      Q2(3,10)=CONST2*ETAH
      Q2(3,12)=CONST1*ABS(MUH)
      Q2(3,15)=CONST1*ABS(XIH)
      Q2(4,4)=VT+2.0D0*ABS(MUH)+2.0D0*ABS(ETAH)
      Q2(4,11)=CONST2*XIH
      Q2(4,13)=CONST1*ABS(MUH)
      Q2(4,19)=CONST1*ABS(ETAH)
      Q2(5,5)=-VT-2.0D0*ABS(XIH)
      Q2(5,12)=-CONST2*MUH
      Q2(5,17)=-CONST2*ETAH
      Q2(5,20)=-CONST1*ABS(XIH)
      Q2(6,6)=-VT-2.0D0*ABS(MUH)
      Q2(6,14)=-CONST1*ABS(MUH)
      Q2(6,15)=-CONST2*XIH
      Q2(6,19)=-CONST2*ETAH
      Q2(7,7)=-VT-2.0D0*ABS(ETAH)
      Q2(7,13)=-CONST2*MUH
      Q2(7,16)=-CONST2*XIH
      Q2(7,18)=-CONST1*ABS(ETAH)
      Q2(8,8)=VT
      Q2(8,14)=CONST2*MUH
      Q2(8,18)=CONST2*ETAH
      Q2(8,20)=CONST2*XIH
      Q2(9,9)=-VT-1.0D1*ABS(MUH)-2.0D0*ABS(ETAH)-
     1         2.0D0*ABS(XIH)
      Q2(9,21)=-CONST1*ABS(ETAH)
      Q2(9,23)=-CONST1*ABS(XIH)
      Q2(10,10)=-VT-2.0D0*ABS(MUH)-1.0D1*ABS(ETAH)-
     1        2.0D0*ABS(XIH)
      Q2(10,21)=-CONST1*ABS(MUH)
      Q2(10,22)=-CONST1*ABS(XIH)
      Q2(11,11)=-VT-2.0D0*ABS(MUH)-2.0D0*ABS(ETAH)-
     1        1.0D1*ABS(XIH)
      Q2(11,22)=-CONST1*ABS(ETAH)
      Q2(11,23)=-CONST1*ABS(MUH)
      Q2(12,12)=VT+1.0D1*ABS(MUH)+2.0D0*ABS(XIH)
      Q2(12,21)=CONST2*ETAH
      Q2(12,25)=CONST1*ABS(XIH)
      Q2(13,13)=VT+1.0D1*ABS(MUH)+2.0D0*ABS(ETAH)
      Q2(13,23)=CONST2*XIH
      Q2(13,26)=CONST1*ABS(ETAH)
      Q2(14,14)=-VT-1.0D1*ABS(MUH)
      Q2(14,25)=-CONST2*XIH
      Q2(14,26)=-CONST2*ETAH
      Q2(15,15)=VT+2.0D0*ABS(MUH)+1.0D1*ABS(XIH)
      Q2(15,22)=CONST2*ETAH
      Q2(15,25)=CONST1*ABS(MUH)
      Q2(16,16)=VT+2.0D0*ABS(ETAH)+1.0D1*ABS(XIH)
      Q2(16,23)=CONST2*MUH
      Q2(16,24)=CONST1*ABS(ETAH)
      Q2(17,17)=VT+1.0D1*ABS(ETAH)+2.0D0*ABS(XIH)
      Q2(17,21)=CONST2*MUH
      Q2(17,24)=CONST1*ABS(XIH)
      Q2(18,18)=-VT-1.0D1*ABS(ETAH)
      Q2(18,24)=-CONST2*XIH
      Q2(18,26)=-CONST2*MUH
      Q2(19,19)=VT+2.0D0*ABS(MUH)+1.0D1*ABS(ETAH)
      Q2(19,22)=CONST2*XIH
      Q2(19,26)=CONST1*ABS(MUH)
      Q2(20,20)=-VT-1.0D1*ABS(XIH)
      Q2(20,24)=-CONST2*ETAH
      Q2(20,25)=-CONST2*MUH
      Q2(21,21)=-VT-1.0D1*ABS(MUH)-1.0D1*ABS(ETAH)-
     1        2.0D0*ABS(XIH)
      Q2(21,27)=-CONST1*ABS(XIH)
      Q2(22,22)=-VT-1.0D1*ABS(ETAH)-1.0D1*ABS(XIH)-
     1        2.0D0*ABS(MUH)
      Q2(22,27)=-CONST1*ABS(MUH)
      Q2(23,23)=-VT-1.0D1*ABS(MUH)-2.0D0*ABS(ETAH)-
     1        1.0D1*ABS(XIH)
      Q2(23,27)=-CONST1*ABS(ETAH)
      Q2(24,24)=VT+1.0D1*ABS(ETAH)+1.0D1*ABS(XIH)
      Q2(24,27)=CONST2*MUH
      Q2(25,25)=VT+1.0D1*ABS(MUH)+1.0D1*ABS(XIH)
      Q2(25,27)=CONST2*ETAH
      Q2(26,26)=VT+1.0D1*ABS(MUH)+1.0D1*ABS(ETAH)
      Q2(26,27)=CONST2*XIH
      Q2(27,27)=-VT-1.0D1*ABS(MUH)-1.0D1*ABS(ETAH)-
     1        1.0D1*ABS(XIH)

*---------------------
      Q2(1,28)=-VOL(I2,J2,J,I,K2)*Q(1)-2.0D0*ABS(MUH)*XNI(1,1,J2,K2)-
     1        2.0D0*ABS(ETAH)*XNJ(1,1,K2)-
     2        2.0D0*ABS(XIH)*XNK(1,1)
      Q2(2,28)=VOL(I2,J2,J,I,K2)*Q(2)+2.0D0*ABS(ETAH)*XNJ(2,1,K2)+
     1        2.0D0*ABS(XIH)*XNK(2,1)
      Q2(3,28)=VOL(I2,J2,J,I,K2)*Q(3)+2.0D0*ABS(MUH)*XNI(2,1,J2,K2)+
     1        2.0D0*ABS(XIH)*XNK(1,2)
      Q2(4,28)=VOL(I2,J2,J,I,K2)*Q(4)+2.0D0*ABS(MUH)*XNI(1,2,J2,K2)+
     1        2.0D0*ABS(ETAH)*XNJ(1,2,K2)
      Q2(5,28)=-VOL(I2,J2,J,I,K2)*Q(5)-2.0D0*ABS(XIH)*XNK(2,2)
      Q2(6,28)=-VOL(I2,J2,J,I,K2)*Q(6)-2.0D0*ABS(MUH)*XNI(2,2,J2,K2)
      Q2(7,28)=-VOL(I2,J2,J,I,K2)*Q(7)-2.0D0*ABS(ETAH)*XNJ(2,2,K2)
      Q2(8,28)=VOL(I2,J2,J,I,K2)*Q(8)
      Q2(9,28)=-VOL(I2,J2,J,I,K2)*Q(9)-CONST1*ABS(MUH)*XNI(1,1,J2,K2)-
     1        2.0D0*ABS(ETAH)*XNJ(3,1,K2)-
     2        2.0D0*ABS(XIH)*XNK(3,1)
      Q2(10,28)=-VOL(I2,J2,J,I,K2)*Q(10)-CONST1*ABS(ETAH)*XNJ(1,1,K2)-
     1        2.0D0*ABS(MUH)*XNI(3,1,J2,K2)-
     2        2.0D0*ABS(XIH)*XNK(1,3)
      Q2(11,28)=-VOL(I2,J2,J,I,K2)*Q(11)-2.0D0*ABS(MUH)*XNI(1,3,J2,K2)-
     1        2.0D0*ABS(ETAH)*XNJ(1,3,K2)-
     2        CONST1*ABS(XIH)*XNK(1,1)
      Q2(12,28)=VOL(I2,J2,J,I,K2)*Q(12)+CONST1*ABS(MUH)*XNI(2,1,J2,K2)+
     1        2.0D0*ABS(XIH)*XNK(3,2)
      Q2(13,28)=VOL(I2,J2,J,I,K2)*Q(13)+2.0D0*ABS(ETAH)*XNJ(3,2,K2)+
     1        CONST1*ABS(MUH)*XNI(1,2,J2,K2)
      Q2(14,28)=-VOL(I2,J2,J,I,K2)*Q(14)-CONST1*ABS(MUH)*XNI(2,2,J2,K2)
      Q2(15,28)=VOL(I2,J2,J,I,K2)*Q(15)+2.0D0*ABS(MUH)*XNI(2,3,J2,K2)+
     1        CONST1*ABS(XIH)*XNK(1,2)
      Q2(16,28)=VOL(I2,J2,J,I,K2)*Q(16)+2.0D0*ABS(ETAH)*XNJ(2,3,K2)+
     1        CONST1*ABS(XIH)*XNK(2,1)
      Q2(17,28)=VOL(I2,J2,J,I,K2)*Q(17)+CONST1*ABS(ETAH)*XNJ(2,1,K2)+
     1        CONST1*ABS(XIH)*XNK(2,3)
      Q2(18,28)=-VOL(I2,J2,J,I,K2)*Q(18)-CONST1*ABS(ETAH)*XNJ(2,2,K2)
      Q2(19,28)=VOL(I2,J2,J,I,K2)*Q(19)+2.0D0*ABS(MUH)*XNI(3,2,J2,K2)+
     1        CONST1*ABS(ETAH)*XNJ(1,2,K2)
      Q2(20,28)=-VOL(I2,J2,J,I,K2)*Q(20)-CONST1*ABS(XIH)*XNK(2,2)
      Q2(21,28)=-VOL(I2,J2,J,I,K2)*Q(21)-
     1        CONST1*ABS(MUH)*XNI(3,1,J2,K2)-
     2        CONST1*ABS(ETAH)*XNJ(3,1,K2)-
     3        2.0D0*CONST1*ABS(XIH)*XNK(3,3)
      Q2(22,28)=-VOL(I2,J2,J,I,K2)*Q(22)-2.0D0*ABS(MUH)*XNI(3,3,J2,K2)
     1        -CONST1*ABS(ETAH)*XNJ(1,3,K2)
     2        -CONST1*ABS(XIH)*XNK(1,3)
      Q2(23,28)=-VOL(I2,J2,J,I,K2)*Q(23)-
     1        CONST1*ABS(MUH)*XNI(1,3,J2,K2)-
     2        2.0D0*ABS(ETAH)*XNJ(3,3,K2)-
     3        CONST1*ABS(XIH)*XNK(3,1)
      Q2(24,28)=VOL(I2,J2,J,I,K2)*Q(24)+CONST1*ABS(ETAH)*XNJ(2,3,K2)+
     1        CONST1*ABS(XIH)*XNK(2,3)
      Q2(25,28)=VOL(I2,J2,J,I,K2)*Q(25)+CONST1*ABS(MUH)*XNI(2,3,J2,K2)+
     1        CONST1*ABS(XIH)*XNK(3,2)
      Q2(26,28)=VOL(I2,J2,J,I,K2)*Q(26)+CONST1*ABS(MUH)*XNI(3,2,J2,K2)+
     1        CONST1*ABS(ETAH)*XNJ(3,2,K2)
      Q2(27,28)=-VOL(I2,J2,J,I,K2)*Q(27)-
     1        CONST1*ABS(MUH)*XNI(3,3,J2,K2)-
     2        CONST1*ABS(ETAH)*XNJ(3,3,K2)-
     3        CONST1*ABS(XIH)*XNK(3,3)
      ENDIF
      DO 190 IEL=1,IELEM**3
      DO 180 JEL=IEL+1,IELEM**3
      Q2(JEL,IEL)=Q2(IEL,JEL)
  180 CONTINUE
  190 CONTINUE
*
      CALL ALSBD(IELEM**3,1,Q2,IER,IELEM**3)
      IF(IER.NE.0) CALL XABORT('SNFTH3: SINGULAR MATRIX.')
*
      IF(IELEM.EQ.1) THEN
         XNK(1,1)=2.0D0*Q2(1,2)-XNK(1,1)
      ELSE IF(IELEM.EQ.2) THEN
         XNK(1,1)=XNK(1,1)+SIGN(1.0,XIH)*CONST0*Q2(4,9)
         XNK(1,2)=XNK(1,2)+SIGN(1.0,XIH)*CONST0*Q2(7,9)
         XNK(2,1)=XNK(2,1)+SIGN(1.0,XIH)*CONST0*Q2(6,9)
         XNK(2,2)=XNK(2,2)+SIGN(1.0,XIH)*CONST0*Q2(8,9)
      ELSE IF(IELEM.EQ.3) THEN
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
*
      IF(JL.LT.ISPLH)THEN
         IF(IELEM.EQ.1) THEN
            XNJ(1,1,K2)=2.0D0*Q2(1,2)-XNJ(1,1,K2)
         ELSEIF(IELEM.EQ.2) THEN
            XNJ(1,1,K2)=XNJ(1,1,K2)+SIGN(1.0,ETAH)*CONST0*Q2(3,9)
            XNJ(1,2,K2)=XNJ(1,2,K2)+SIGN(1.0,ETAH)*CONST0*Q2(7,9)
            XNJ(2,1,K2)=XNJ(2,1,K2)+SIGN(1.0,ETAH)*CONST0*Q2(5,9)
            XNJ(2,2,K2)=XNJ(2,2,K2)+SIGN(1.0,ETAH)*CONST0*Q2(8,9)
         ELSEIF(IELEM.EQ.3) THEN
            XNJ(1,1,K2)=2.0D0*Q2(1,28)+CONST1*Q2(10,28)-XNJ(1,1,K2)
            XNJ(2,1,K2)=2.0D0*Q2(2,28)+CONST1*Q2(17,28)-XNJ(2,1,K2)
            XNJ(1,2,K2)=2.0D0*Q2(4,28)+CONST1*Q2(19,28)-XNJ(1,2,K2)
            XNJ(2,2,K2)=2.0D0*Q2(7,28)+CONST1*Q2(18,28)-XNJ(2,2,K2)
            XNJ(3,1,K2)=2.0D0*Q2(9,28)+CONST1*Q2(21,28)-XNJ(3,1,K2)
            XNJ(1,3,K2)=2.0D0*Q2(11,28)+CONST1*Q2(22,28)-XNJ(1,3,K2)
            XNJ(3,2,K2)=2.0D0*Q2(13,28)+CONST1*Q2(26,28)-XNJ(3,2,K2)
            XNJ(2,3,K2)=2.0D0*Q2(16,28)+CONST1*Q2(24,28)-XNJ(2,3,K2)
            XNJ(3,3,K2)=2.0D0*Q2(23,28)+CONST1*Q2(27,28)-XNJ(3,3,K2)
         ENDIF
      ELSEIF((JL.EQ.ISPLH).AND.(IHEXJ.LE.NHEX))THEN
         I3=I2
         C1=1.0D0
         IF((J.EQ.1).AND.(ILOZJ.EQ.3)) THEN
            I3=ISPLH+1 -I2
            C1=-1.0D0
         ENDIF
         IF(IELEM.EQ.1) THEN
            BFLUX(ISIDEJ,INDEXJ,I3,1,1,K2) = 2.0D0*Q2(1,2)-XNJ(1,1,K2)
         ELSEIF(IELEM.EQ.2) THEN
            BFLUX(ISIDEJ,INDEXJ,I3,1,1,K2)= XNJ(1,1,K2)+SIGN(1.0,ETAH)*
     1         CONST0*Q2(3,9)
            BFLUX(ISIDEJ,INDEXJ,I3,1,2,K2)=(XNJ(1,2,K2)+SIGN(1.0,ETAH)*
     1         CONST0*Q2(7,9))*C1
            BFLUX(ISIDEJ,INDEXJ,I3,2,1,K2)=(XNJ(2,1,K2)+SIGN(1.0,ETAH)*
     1         CONST0*Q2(5,9))*C1
            BFLUX(ISIDEJ,INDEXJ,I3,2,2,K2)= XNJ(2,2,K2)+SIGN(1.0,ETAH)*
     1         CONST0*Q2(8,9)
         ELSEIF(IELEM.EQ.3) THEN
            BFLUX(ISIDEJ,INDEXJ,I3,1,1,K2)= 2.0D0*Q2(1,28)+CONST1*
     1         Q2(10,28)-XNJ(1,1,K2)
            BFLUX(ISIDEJ,INDEXJ,I3,2,1,K2)=(2.0D0*Q2(2,28)+CONST1*
     1         Q2(17,28)-XNJ(2,1,K2))*C1
            BFLUX(ISIDEJ,INDEXJ,I3,1,2,K2)=(2.0D0*Q2(4,28)+CONST1*
     1         Q2(19,28)-XNJ(1,2,K2))*C1
            BFLUX(ISIDEJ,INDEXJ,I3,2,2,K2)= 2.0D0*Q2(7,28)+CONST1*
     1         Q2(18,28)-XNJ(2,2,K2)
            BFLUX(ISIDEJ,INDEXJ,I3,3,1,K2)= 2.0D0*Q2(9,28)+CONST1*
     1         Q2(21,28)-XNJ(3,1,K2)
            BFLUX(ISIDEJ,INDEXJ,I3,1,3,K2)= 2.0D0*Q2(11,28)+CONST1*
     1         Q2(22,28)-XNJ(1,3,K2)
            BFLUX(ISIDEJ,INDEXJ,I3,3,2,K2)=(2.0D0*Q2(13,28)+CONST1*
     1         Q2(26,28)-XNJ(3,2,K2))*C1
            BFLUX(ISIDEJ,INDEXJ,I3,2,3,K2)=(2.0D0*Q2(16,28)+CONST1*
     1         Q2(24,28)-XNJ(2,3,K2))*C1
            BFLUX(ISIDEJ,INDEXJ,I3,3,3,K2)= 2.0D0*Q2(23,28)+CONST1*
     1         Q2(27,28)-XNJ(3,3,K2)
         ENDIF
      ENDIF
*
      IF(IL.LT.ISPLH)THEN
         IF(IELEM.EQ.1) THEN
            XNI(1,1,J2,K2) = 2.0D0*Q2(1,2)-XNI(1,1,J2,K2)
         ELSEIF(IELEM.EQ.2) THEN
            XNI(1,1,J2,K2)=XNI(1,1,J2,K2)+SIGN(1.0,MUH)*CONST0*Q2(2,9)
            XNI(1,2,J2,K2)=XNI(1,2,J2,K2)+SIGN(1.0,MUH)*CONST0*Q2(6,9)
            XNI(2,1,J2,K2)=XNI(2,1,J2,K2)+SIGN(1.0,MUH)*CONST0*Q2(5,9)
            XNI(2,2,J2,K2)=XNI(2,2,J2,K2)+SIGN(1.0,MUH)*CONST0*Q2(8,9)
         ELSEIF(IELEM.EQ.3) THEN
            XNI(1,1,J2,K2)=2.0D0*Q2(1,28)+CONST1*Q2(9,28)-
     1         XNI(1,1,J2,K2)
            XNI(2,1,J2,K2)=2.0D0*Q2(3,28)+CONST1*Q2(12,28)-
     1         XNI(2,1,J2,K2)
            XNI(1,2,J2,K2)=2.0D0*Q2(4,28)+CONST1*Q2(13,28)-
     1         XNI(1,2,J2,K2)
            XNI(2,2,J2,K2)=2.0D0*Q2(6,28)+CONST1*Q2(14,28)-
     1         XNI(2,2,J2,K2)
            XNI(1,3,J2,K2)=2.0D0*Q2(11,28)+CONST1*Q2(23,28)-
     1         XNI(1,3,J2,K2)
            XNI(3,1,J2,K2)=2.0D0*Q2(10,28)+CONST1*Q2(21,28)-
     1         XNI(3,1,J2,K2)
            XNI(3,2,J2,K2)=2.0D0*Q2(19,28)+CONST1*Q2(26,28)-
     1         XNI(3,2,J2,K2)
            XNI(2,3,J2,K2)=2.0D0*Q2(15,28)+CONST1*Q2(25,28)-
     1         XNI(2,3,J2,K2)
            XNI(3,3,J2,K2)=2.0D0*Q2(22,28)+CONST1*Q2(27,28)-
     1         XNI(3,3,J2,K2)
         ENDIF
      ELSEIF((IL.EQ.ISPLH).AND.(IHEXI.LE.NHEX))THEN
         J3=J2
         C1=1.0D0
         IF((J.EQ.3).AND.(ILOZI.EQ.1)) THEN
            J3=ISPLH+1-J2
            C1=-1.0D0
         ENDIF
         IF(IELEM.EQ.1) THEN
            BFLUX(ISIDEI,INDEXI,J3,1,1,K2)=2.0D0*Q2(1,2)-XNI(1,1,J2,K2)
         ELSEIF(IELEM.EQ.2) THEN
            BFLUX(ISIDEI,INDEXI,J3,1,1,K2)= XNI(1,1,J2,K2)+
     1         SIGN(1.0,MUH)*CONST0*Q2(2,9)
            BFLUX(ISIDEI,INDEXI,J3,1,2,K2)=(XNI(1,2,J2,K2)+
     1         SIGN(1.0,MUH)*CONST0*Q2(6,9))*C1
            BFLUX(ISIDEI,INDEXI,J3,2,1,K2)=(XNI(2,1,J2,K2)+
     1         SIGN(1.0,MUH)*CONST0*Q2(5,9))*C1
            BFLUX(ISIDEI,INDEXI,J3,2,2,K2)= XNI(2,2,J2,K2)+
     1         SIGN(1.0,MUH)*CONST0*Q2(8,9)
         ELSEIF(IELEM.EQ.3) THEN
            BFLUX(ISIDEI,INDEXI,J3,1,1,K2)= 2.0D0*Q2(1,28)+CONST1*
     1         Q2(9,28)-XNI(1,1,J2,K2)
            BFLUX(ISIDEI,INDEXI,J3,2,1,K2)=(2.0D0*Q2(3,28)+CONST1*
     1         Q2(12,28)-XNI(2,1,J2,K2))*C1
            BFLUX(ISIDEI,INDEXI,J3,1,2,K2)=(2.0D0*Q2(4,28)+CONST1*
     1         Q2(13,28)-XNI(1,2,J2,K2))*C1
            BFLUX(ISIDEI,INDEXI,J3,2,2,K2)= 2.0D0*Q2(6,28)+CONST1*
     1         Q2(14,28)-XNI(2,2,J2,K2)
            BFLUX(ISIDEI,INDEXI,J3,1,3,K2)= 2.0D0*Q2(11,28)+CONST1*
     1         Q2(23,28)-XNI(1,3,J2,K2)
            BFLUX(ISIDEI,INDEXI,J3,3,1,K2)= 2.0D0*Q2(10,28)+CONST1*
     1         Q2(21,28)-XNI(3,1,J2,K2)
            BFLUX(ISIDEI,INDEXI,J3,3,2,K2)=(2.0D0*Q2(19,28)+CONST1*
     1         Q2(26,28)-XNI(3,2,J2,K2))*C1
            BFLUX(ISIDEI,INDEXI,J3,2,3,K2)=(2.0D0*Q2(15,28)+CONST1*
     1         Q2(25,28)-XNI(2,3,J2,K2))*C1
            BFLUX(ISIDEI,INDEXI,J3,3,3,K2)= 2.0D0*Q2(22,28)+CONST1*
     1         Q2(27,28)-XNI(3,3,J2,K2)
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
      FLUX_G(:,:,:,:,:,IGEFF)=FLUX_G(:,:,:,:,:,IGEFF)+FLUX(:,:,:,:,:)
  290 CONTINUE
  300 CONTINUE
*$OMP END PARALLEL DO
  310 CONTINUE
      DO 320 IGEFF=1,NGEFF
        IF(.NOT.INCONV(IGEFF)) GO TO 320
        IG=INGIND(IGEFF)
        FUNKNO(:L4,IG)=
     1  RESHAPE(REAL(FLUX_G(:IELEM**3,:NSCT,:3*ISPLH**2,:NHEX,:LZ,
     2  IGEFF)), (/L4/))
  320 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(FLUX_G,FLUX,BFLUX,INDANG)
      RETURN
  500 FORMAT(16H SNFTH3: thread=,I8,12H --->(group=,I4,7H angle=,I4,1H))
      END
