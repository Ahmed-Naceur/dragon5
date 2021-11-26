*DECK MOCCHR
      SUBROUTINE MOCCHR(NDIM,NANI,NFUNL,NMU,XMUANG,PHI1,PHI2,R)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Generates all spherical harmonics R(L,M) at point (XMUANG,PHI1).
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* NDIM    number of dimensions for the geometry.
* NANI    scattering anisotropy (=0 for isotropic scattering).
* NFUNL   number of spherical harmonics per polar angle.
* NMU     order of the polar quadrature in 2D and 1 in 3D.
* XMUANG  cosines of the different tracking polar angles
*         (polar quadrature in 2D and tracking angles in 3D).
* PHI1    first cosine of the tracking azimuthal angle.
* PHI2    second cosine of the tracking azimuthal angle.
*
*Parameters: output
* R       spherical harmonics.
*
*Comments:
*  for 0 <= L <= NANI (and for -L <= M <= L)
*  and for all angles (for 1 <= IMU <= NMU).
*    Definition of spherical harmonics R(L,M):
*    R(L, M)= FACT(L,M)*P(L,M)*COS(M*PHI1)  for 0 < M <= L
*    R(L, 0)=           P(L,0)
*    R(L,-M)= FACT(L,M)*P(L,M)*SIN(M*PHI1)  for 0 < M <= L
*    where FACT(L,M)= SQRT( 2*(L-M)!/(L+M)! )
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER          NDIM,NANI,NFUNL,NMU
      REAL             XMUANG(NMU),R(NMU,NFUNL)
      DOUBLE PRECISION PHI1,PHI2
*----
*  LOCAL VARIABLES
*----
      INTEGER          IMU,L,M,LPM,LMMP1,IND,NSELEC
      LOGICAL          LROK
      DOUBLE PRECISION DPHI,DCOP2,DSIP2,DL00,DL10,DL01,DL11,DM0,DM1,
     >                 DM2,DCOS,DSIN,DCOT,DFAC,DZERO,DONE
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: RWORK
      PARAMETER ( DZERO= 0.0D0, DONE= 1.0D0 )
*
*     INDEX FOR MATRIX 'RWORK'
      IND(L,M)= L*(L+1) + M + 1
*
***** FIRST, COMPUTES ALL SPHERICAL HARMONICS
*
*     GENERATES ALL LEGENDRE POLYNOMIALS P(L,M) AT POINT XMU
*               FOR 0 <= L <= NANI (AND FOR 0 <= M <= L)
*     USING THE RECURENCE RELATIONS:
*     P(L+1,M)= ((2*L+1)*XMU*P(L,M)-(L+M)*P(L-1,M))/(L-M+1)
*     P(L,M+1)= 2*M*XMU*P(L,M)/SQRT(1-XMU**2)-(L+M)*(L-M+1)*P(L,M-1)
*
*     ESTABLISH SIMPLE CASES
      ALLOCATE(RWORK(NMU,(NANI+1)*(NANI+1)))
      IF( NANI.GE.0 )THEN
         DO 10 IMU= 1, NMU
            RWORK(IMU,IND(0,0))= DONE
   10    CONTINUE
      ELSE
         CALL XABORT('MOCCHR: THE FIRST ARGUMENT MUST NON NEGATIVE')
      ENDIF
*
      IF( NANI.GE.1 )THEN
         DPHI  = SIGN(ACOS(PHI1),PHI2)
         DO 40 IMU= 1, NMU
            DCOS=   DBLE(XMUANG(IMU))
            DSIN=   SQRT( DONE - DCOS  *DCOS   )
            RWORK(IMU,IND(1,-1))= DSIN*PHI2
            RWORK(IMU,IND(1, 0))= DCOS
            RWORK(IMU,IND(1, 1))= DSIN*PHI1
*
            IF( NANI.GE.2 )THEN
*              RECURENCE PLG(IND(L,0)),PLG(IND(L,1) FOR 2 <= L <= NANI
               DCOT= DCOS/DSIN
               DL00= DONE
               DL10= DCOS
               DL01= DZERO
               DL11= DSIN
               DO 30 L= 1, NANI-1
*                 IF M=1, L=L+1 THEN L+M=L+2, L-M+1=L+1
                  LPM=   L + 2
                  LMMP1= L + 1
                  DFAC= (DONE+DONE)/DBLE(LPM*LMMP1)
                  DM0=(DBLE(2*L+1)*DCOS*DL10 - DBLE(L)*DL00)/DBLE(L+1)
                  DM1=(DBLE(2*L+1)*DCOS*DL11 - DBLE(L+1)*DL01)/DBLE(L)
*                 ESTABLISH RELATIONS FOR L=L+1, ABS(M)<=1
                  RWORK(IMU,IND(L+1,-1))= SQRT(DFAC)*DM1*PHI2
                  RWORK(IMU,IND(L+1, 0))= DM0
                  RWORK(IMU,IND(L+1, 1))= SQRT(DFAC)*DM1*PHI1
                  DL00= DL10
                  DL01= DL11
                  DL10= DM0
                  DL11= DM1
*                 RECURENCE PLG(IND(L,M)) FOR 2 <= M <= L
                  DO 20 M= 1, L
*                    HERE DM0=PLG(L+1,0), DM1=PLG(L+1,1)
*                    ESTABLISH RELATIONS FOR L=L+1, ABS(M)>1
                     DFAC= DFAC/DBLE((LMMP1-1)*(LPM+1))
                     DCOP2= COS(DBLE(M+1)*DPHI)
                     DSIP2= SIN(DBLE(M+1)*DPHI)
                     DM2= DBLE(2*M)*DCOT*DM1 - DBLE(LPM*LMMP1)*DM0
                     RWORK(IMU,IND(L+1,-(M+1)))= SQRT(DFAC)*DM2*DSIP2
                     RWORK(IMU,IND(L+1, M+1 ))= SQRT(DFAC)*DM2*DCOP2
                     DM0= DM1
                     DM1= DM2
                     LPM=   LPM   + 1
                     LMMP1= LMMP1 - 1
   20             CONTINUE
   30          CONTINUE
            ENDIF
   40    CONTINUE
      ENDIF
*
***** SELECTS THE GOOD SPHERICAL HARMONICS RWORK(L,M) FUNCTIONS
*             FOR NDIM=1(SLAB),2(TWO-D RECT),3(THREE-D).
*     COMPRESSES RWORK INTO R.
*
      NSELEC= 0
      LROK= .FALSE.
      DO 80 L= 0, NANI
         DO 70 M= -L, L
            IF( NDIM.EQ.1 )THEN
               LROK= M.EQ.0
            ELSEIF( NDIM.EQ.2 )THEN
               LROK= MOD(L+M,2).EQ.0
            ELSEIF( NDIM.EQ.3 )THEN
               LROK= .TRUE.
            ENDIF
            IF( LROK )THEN
               NSELEC= NSELEC+1
               DO 50 IMU= 1, NMU
                  R(IMU,NSELEC)= REAL(RWORK(IMU,IND(L,M)))
   50          CONTINUE
            ENDIF
   70   CONTINUE
   80 CONTINUE
      IF(NSELEC.NE.NFUNL) CALL XABORT('MOCCHR: INVALID NSELEC')
      DEALLOCATE(RWORK)
*
      RETURN
      END
