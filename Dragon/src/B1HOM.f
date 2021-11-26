*DECK B1HOM
      SUBROUTINE B1HOM (IPMACR,LEAKSW,NUNKNO,OPTION,TYPE,NGRO,IPAS,NBM,
     1                  NFISSI,VOL,MAT,KEYFLX,FLUX,REFKEF,IMPX,D,GAMMA,
     2                  ALAM1,INORM,B2,OLDBIL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Homogenization of the unit cell and solution of the B-n equations.
* The cross section information is found on LCM.
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
* IPMACR  pointer to the macrolib LCM object (L_MACROLIB signature).
* LEAKSW  leakage flag (=.TRUE. if leakage is present on the outer
*         surface).
* NUNKNO  number of flux/current unknowns.
* OPTION  type of leakage coefficients; can be 'LKRD' (recover leakage
*         coefficients in Macrolib), 'RHS' (recover leakage coefficients
*         in RHS flux object), 'B0' (B-0), 'P0' (P-0), 'B1' (B-1),
*         'P1' (P-1), 'B0TR' (B-0 with transport correction) or 'P0TR'
*         (P-0 with transport correction).
* TYPE    type of buckling iteration.
*         Can be 'DIFF' (do a B-0 calculation of D(NGRO) and exit);
*                'K' (do a B-n calculation with keff search);
*                'B' (do a B-n calculation with buckling search);
*                'L' (do a B-n calculation with buckling search
*                     for a problem with few or no fission).
* NGRO    number of groups.
* IPAS    number of volumes.
* NBM     number of mixtures.
* NFISSI  maximum number of fission spectrum assigned to a mixture.
* VOL     volumes.
* MAT     mixture number of each volume.
* KEYFLX  position of each flux in the unknown vector.
* FLUX    direct unknown vector.
* REFKEF  target K-effective for type B or type L calculations
* IMPX    print flag.
* INORM   type of leakage model:
*         =1: Diffon; =2: Ecco; =3: Tibere.
* B2      original direction dependant buckling.
* OLDBIL  previous norm of the flux.
*
*Parameters: output
* D       diffusion coefficients.
* GAMMA   gamma factors.
* ALAM1   effective multiplication factor.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER*4 OPTION,TYPE
      TYPE(C_PTR) IPMACR
      LOGICAL LEAKSW
      INTEGER NUNKNO,NGRO,IPAS,NBM,NFISSI,MAT(IPAS),KEYFLX(IPAS),IMPX,
     1 INORM
      REAL VOL(IPAS),FLUX(NUNKNO,NGRO),D(NGRO),GAMMA(NGRO),B2(4)
      DOUBLE PRECISION REFKEF,ALAM1,OLDBIL
*----
*  LOCAL VARIABLES
*----
      INTEGER IDEL(2)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ0,IJJ1,NJJ0,NJJ1
      REAL, ALLOCATABLE, DIMENSION(:) :: ST,SA,SFNU,XHI,SCAT0,SCAT1,FL2
      DOUBLE PRECISION B2HOM,CAET,A2,CURN,B2T(3)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: PHI
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(PHI(NGRO))
*
      IAN=0
      IF ((OPTION.EQ.'B0').OR.(OPTION.EQ.'P0')) THEN
         IAN=0
      ELSE IF ((OPTION.EQ.'B1').OR.(OPTION.EQ.'P1')) THEN
         IAN=1
      ELSE IF ((OPTION.EQ.'B0TR').OR.(OPTION.EQ.'P0TR')) THEN
         IAN=-1
      ENDIF
      ALLOCATE(IJJ0(NGRO),IJJ1(NGRO),NJJ0(NGRO),NJJ1(NGRO))
      CALL B1HXS1(IPMACR,NGRO,NBM,IAN,NFISSI,IJJ0,IJJ1,NJJ0,NJJ1,IDEL)
*
      ALLOCATE(ST(NGRO),SA(NGRO),SFNU(NGRO),XHI(NGRO),SCAT0(IDEL(1)),
     1 SCAT1(IDEL(2)))
      IF(INORM.EQ.2) THEN
*        ECCO-TYPE ISOTROPIC STREAMING.
         CALL B1HXS3(NUNKNO,IPMACR,IPAS,NGRO,NBM,IAN,VOL,MAT,KEYFLX,
     1   FLUX,IJJ0,IJJ1,NJJ0,NJJ1,IDEL,PHI,ST,SCAT0,SCAT1,NGROIN)
      ELSE IF(INORM.EQ.3) THEN
*        TIBERE-TYPE ANISOTROPIC STREAMING.
         IF(B2(4).EQ.0.0) THEN
            B2T(1)=0.33333333333333D0
            B2T(2)=B2T(1)
            B2T(3)=B2T(1)
         ELSE
            B2T(1)=DBLE(B2(1))/DBLE(B2(4))
            B2T(2)=DBLE(B2(2))/DBLE(B2(4))
            B2T(3)=DBLE(B2(3))/DBLE(B2(4))
         ENDIF
         ALLOCATE(FL2(2*NUNKNO*NGRO))
         IOF=0
         DO 30 IGRO=1,NGRO
           DO 10 IUNK=1,NUNKNO/4
           IOF=IOF+1
           FL2(IOF)=FLUX(IUNK,IGRO)
   10      CONTINUE
           DO 20 IUNK=1,NUNKNO/4
           IOF=IOF+1
           CURN=0.0D0
           DO 15 IDIR=1,3
           CURN=CURN+B2T(IDIR)*FLUX(NUNKNO/4*IDIR+IUNK,IGRO)
   15      CONTINUE
           FL2(IOF)=REAL(CURN)
   20      CONTINUE
   30    CONTINUE
         CALL B1HXS3(NUNKNO/2,IPMACR,IPAS,NGRO,NBM,IAN,VOL,MAT,KEYFLX,
     1   FL2,IJJ0,IJJ1,NJJ0,NJJ1,IDEL,PHI,ST,SCAT0,SCAT1,NGROIN)
         DEALLOCATE(FL2)
      ENDIF
      CALL B1HXS2(NUNKNO,IPMACR,IPAS,NGRO,NBM,IAN,NFISSI,VOL,MAT,
     1 KEYFLX,FLUX,IJJ0,IJJ1,NJJ0,NJJ1,IDEL,PHI,SA,ST,SFNU,XHI,SCAT0,
     2 SCAT1,NGROIN,INORM)
*
      IF(LEAKSW) THEN
*        Obtain leakage coefficients using a Todorova approximation.
         CALL B1TODO(OPTION,TYPE,IMPX,NGRO,IJJ1,NJJ1,IDEL,PHI,ST,SFNU,
     1   SCAT1,OLDBIL,B2(4),D,ALAM1)
         GAMMA(:NGRO)=1.0
         GO TO 130
      ENDIF
      B2HOM=DBLE(B2(4))
      CALL B1DIF(OPTION,TYPE,NGRO,ST,SFNU,XHI,IJJ0,IJJ1,NJJ0,NJJ1,SCAT0,
     1 SCAT1,REFKEF,IMPX,D,GAMMA,B2HOM,ALAM1,CAET,A2,PHI)
      B2(4)=REAL(B2HOM)
*
      IF (TYPE.EQ.'DIFF') GO TO 130
*----
*  NORMALIZE THE DRAGON FLUX USING THE FUNDAMENTAL B1 SOLUTION
*----
      IF(INORM.EQ.1) THEN
        DO 60 I=1,NGRO
          CAET=0.0D0
          DO 40 L=1,IPAS
            CAET=CAET+VOL(L)*FLUX(KEYFLX(L),I)
   40     CONTINUE
          CAET=PHI(I)/CAET
          DO 50 L=1,NUNKNO
            FLUX(L,I)=FLUX(L,I)*REAL(CAET)
   50     CONTINUE
   60   CONTINUE
      ELSE IF(INORM.EQ.2) THEN
        DO 90 I=1,NGRO
          CAET=0.0D0
          CURN=0.0D0
          DO 70 L=1,IPAS
            CAET=CAET+VOL(L)*FLUX(KEYFLX(L),I)
            CURN=CURN+VOL(L)*FLUX(KEYFLX(L)+NUNKNO/2,I)
   70     CONTINUE
          CAET=PHI(I)/CAET
          CURN=PHI(I)*D(I)/CURN
          DO 80 L=1,NUNKNO/2
            FLUX(L,I)=FLUX(L,I)*REAL(CAET)
            FLUX(L+NUNKNO/2,I)=FLUX(L+NUNKNO/2,I)*REAL(CURN)
   80     CONTINUE
   90   CONTINUE
      ELSE IF(INORM.EQ.3) THEN
        IF(B2(4).EQ.0.0.OR.
     >    (B2(1).EQ.0.0.AND.B2(2).EQ.0.0.AND.B2(3).EQ.0.0)) THEN
          B2T(1)=0.33333333333333D0
          B2T(2)=B2T(1)
          B2T(3)=B2T(1)
        ELSE
          B2HOM=1.0D0/(DBLE(B2(1))+DBLE(B2(2))+DBLE(B2(3)))
          B2T(1)=B2HOM*DBLE(B2(1))
          B2T(2)=B2HOM*DBLE(B2(2))
          B2T(3)=B2HOM*DBLE(B2(3))
        ENDIF
        DO 120 I=1,NGRO
          CAET=0.0D0
          CURN=0.0D0
          DO 100 L=1,IPAS
          CAET=CAET+VOL(L)*FLUX(KEYFLX(L),I)
          CURN=CURN+B2T(1)*FLUX(KEYFLX(L)+NUNKNO/4,I)*VOL(L)
     >             +B2T(2)*FLUX(KEYFLX(L)+NUNKNO/2,I)*VOL(L)
     >             +B2T(3)*FLUX(KEYFLX(L)+3*NUNKNO/4,I)*VOL(L)
  100     CONTINUE
          CAET=PHI(I)/CAET
          CURN=PHI(I)*D(I)/CURN
          DO 110 L=1,IPAS
          FLUX(KEYFLX(L),I)=FLUX(KEYFLX(L),I)*REAL(CAET)
          FLUX(KEYFLX(L)+NUNKNO/4,I)=FLUX(KEYFLX(L)+NUNKNO/4,I)*
     1    REAL(CURN)
          FLUX(KEYFLX(L)+NUNKNO/2,I)=FLUX(KEYFLX(L)+NUNKNO/2,I)*
     1    REAL(CURN)
          FLUX(KEYFLX(L)+3*NUNKNO/4,I)=FLUX(KEYFLX(L)+3*NUNKNO/4,I)*
     1    REAL(CURN)
  110     CONTINUE
  120   CONTINUE
      ENDIF
*
  130 DEALLOCATE(SCAT1,SCAT0,XHI,SFNU,SA,ST)
      DEALLOCATE(NJJ1,NJJ0,IJJ1,IJJ0)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(PHI)
      RETURN
      END
