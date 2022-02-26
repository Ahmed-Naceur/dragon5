*DECK NXTQAS
      SUBROUTINE NXTQAS(IPRINT,NDIM  ,AZMQUA,NANGL ,NQUAD ,NBANGL,
     >                  DQUAD ,DANGLT,DDENWT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To define quadrature angles for a given tracking option.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s):
*  G. Marleau, R. Roy, M. Hampartzounian
*
*Parameters: input
* IPRINT  print level.
* NDIM    number of dimensions for geometry.
* AZMQUA  quadrature type.
* NANGL   quadrature order.
* NQUAD   number of quadrant (in 3-D) and quarter (in 2-D).
* NBANGL  number of angles.
* DQUAD   relative density of each quadrant.
*
*Parameters: output
* DANGLT  angles.
* DDENWT  angular density for each angle.
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*  Extracted from the subroutine XELTS2 of EXCELL.
*
*----------
*
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          IPRINT,NDIM,AZMQUA,NANGL,NQUAD,NBANGL
      DOUBLE PRECISION DQUAD(NQUAD)
      DOUBLE PRECISION DANGLT(NDIM,NQUAD,NBANGL),DDENWT(NQUAD,NBANGL)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTQAS')
      DOUBLE PRECISION DZERO,DONE,DTWO
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0,DTWO=2.0D0)
*----
*  Functions
*----
      DOUBLE PRECISION XDRCST,PI
*----
*  Local variables
*----
      INTEGER          IANG,IQUAD,IDIR
      DOUBLE PRECISION DTHETA,THETA,DDA,COST,SINT
*----
*  Processing starts:
*  print routine openning output header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6000) NAMSBR
        WRITE(IOUT,6010) NDIM,AZMQUA,NANGL,NQUAD
      ENDIF
      PI=XDRCST('Pi',' ')
      IF(NDIM .EQ. 2) THEN
        IF(AZMQUA .EQ. 1) THEN
*----
*  Trapezoidal quadrature
*  NBANGL point quadrature for both:
*  (1/2*Pi)*Integral(0,Pi/2) and (1/2*Pi)*Integral(Pi/2,Pi)
*  Quadrature weight = Pi/(2*NBANGL)
*  DENSITY=(2*Pi)/weight=4*NBANGL
*----
          THETA =PI/DBLE(4*NANGL)
          DTHETA=DTWO*THETA
          DDA=DBLE(4*NBANGL)
          DO IANG=1,NBANGL
            COST=COS(THETA)
            SINT=SIN(THETA)
            DANGLT(1,1,IANG)=COS(THETA)
            DANGLT(2,1,IANG)=SIN(THETA)
            DDENWT(1,IANG)=DQUAD(1)*DDA
            DANGLT(1,2,IANG)=-SIN(THETA)
            DANGLT(2,2,IANG)=COS(THETA)
            DDENWT(2,IANG)=DQUAD(2)*DDA
            THETA=THETA+DTHETA
          ENDDO
        ELSE
          WRITE(IOUT,9000) NAMSBR,AZMQUA
          CALL XABORT(NAMSBR//': INVALID QUADRATURE OPTION IN 2D')
        ENDIF
      ELSE IF(NDIM .EQ. 3) THEN
        IF(AZMQUA .EQ. 1) THEN
          CALL NXTQEW(NDIM  ,NANGL ,NQUAD ,NBANGL,DQUAD ,
     >                DANGLT,DDENWT)
        ELSE IF(AZMQUA .EQ. 4) THEN
          CALL NXTQLC(NDIM  ,NANGL ,NQUAD ,NBANGL,DQUAD ,
     >                DANGLT,DDENWT)
        ELSE IF(AZMQUA .EQ. 5) THEN
          CALL NXTQLT(NDIM  ,NANGL ,NQUAD ,NBANGL,DQUAD ,
     >                DANGLT,DDENWT)
        ELSE IF(AZMQUA .EQ. 6) THEN
          CALL NXTLSN(NDIM  ,NANGL ,NQUAD ,NBANGL,DQUAD ,
     >                DANGLT,DDENWT)
        ELSE IF(AZMQUA .EQ. 7) THEN
          CALL NXTQRN(NDIM  ,NANGL ,NQUAD ,NBANGL,DQUAD ,
     >                DANGLT,DDENWT)
        ENDIF
      ENDIF
*----
*  Processing finished: return
*----
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6011)
        DO IANG=1,NBANGL
          DO IQUAD=1,NQUAD
            IF(DDENWT(IQUAD,IANG) .GT. DZERO) THEN
              WRITE(IOUT,6012) IANG,IQUAD,
     >                         (DANGLT(IDIR,IQUAD,IANG),IDIR=1,NDIM),
     >                          DDENWT(IQUAD,IANG)
            ENDIF
          ENDDO
        ENDDO
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6010 FORMAT(1X,'NDIM  =',I8,1X,'AZMQUA=',I8,
     >       1X,'NANGL =',I8,1X,'NQUAD =',I8)
 6011 FORMAT(' Tracking directions and weights '/
     >       1X,'     Angle',1X,'  Quadrant',
     >       1X,'     Directions and weight')
 6012 FORMAT(2(1X,I10),4(2X,F24.14))
 9000 FORMAT(A6,': AZMQUA=',I5,' is invalid in 2D') 
      END
