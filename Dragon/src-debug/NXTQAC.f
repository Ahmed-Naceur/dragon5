*DECK NXTQAC
      SUBROUTINE NXTQAC(IPRINT,NDIM  ,NANGL ,NBANGL,DENUSR,ABSC  ,
     >                  RCIRC ,AZMQUA,IPER  ,
     >                  DANGLT,DDENWT,DNSANG,NBSANG,DDANG )
*
*-----------------------------------------------------------------------
*
*Purpose:
* To define quadrature angles for cyclic tracking.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s):
*  G. Marleau, R. Roy
*
*Parameters: input
* IPRINT  print level.
* NDIM    number of dimensions for geometry.
* NANGL   quadrature order.
* NBANGL  number of angles.
* DENUSR  requested density for spatial tracking.
* ABSC    multidimensional width of the cell.
* RCIRC   radius of circle surrounding geometry.
* AZMQUA  tracking type.
* IPER    cell periodicity factor in each direction.
*
*Parameters: output
* DANGLT  director cosines of angles.
* DDENWT  angular density for each angle.
* DNSANG  spatial density required.
* NBSANG  number of segments for each angles.
* DDANG   angles.
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
      INTEGER          IPRINT,NDIM,NANGL,NBANGL
      DOUBLE PRECISION DENUSR,ABSC(NDIM),RCIRC
      INTEGER          AZMQUA,IPER(3)
      DOUBLE PRECISION DANGLT(NDIM,NBANGL,4),DDENWT(NBANGL,4),
     >                 DNSANG(NBANGL)
      INTEGER          NBSANG(5,NBANGL)
      DOUBLE PRECISION DDANG(NBANGL)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTQAC')
      DOUBLE PRECISION DZERO,DONE,DTWO
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0,DTWO=2.0D0)
*----
*  Local variables
*----
      INTEGER          ISUM,IDEB,ISTRID,IANG,ITX,ITY,INDC(2),NTRAC,IOF,
     >                 IA,IB,IC,IPERG,IDIR
      DOUBLE PRECISION DDAX(9),DENLIN
*----
*  Start processing
*----
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6000) NAMSBR
        WRITE(IOUT,6010) NDIM,AZMQUA,NANGL,NBANGL
      ENDIF
      IF(NDIM .NE. 2) CALL XABORT(NAMSBR//
     >': Cyclic tracking works only in 2-D')
*----
*  1. Define angles
*     IWT option is:
*       0<= theta <= Pi/2
*       ITX=0,NBANGL-1
*       ITY=NBANGL-ITX
*     MEDI option is
*       0< theta < Pi/2
*       ITX=1,2*NBANGL,2
*       ITY=2*NBANGL-ITX
*     EQW2 option is
*       0< theta < Pi/2
*       ITX=1,NBANGL
*       ITY=NBANGL-ITX+1
*----
      ISUM=0
      IDEB=0
      ISTRID=0
      IOF=0
      IF(AZMQUA .EQ. 1) THEN
        ISUM=NBANGL-1
        IDEB=0
        ISTRID=1
      ELSE IF(AZMQUA .EQ. 3) THEN
        ISUM=2*NBANGL-1
        IDEB=1
        ISTRID=2
        IOF=1
      ELSE IF(AZMQUA .EQ. 8) THEN
        ISUM=NBANGL
        IDEB=1
        ISTRID=1
        IOF=1
      ELSE
        CALL XABORT(NAMSBR//': Invalid quadrature')
      ENDIF
      IPERG=MIN(IPER(1)*IPER(2),2)
      IANG=0
      DO ITX=IDEB,ISUM,ISTRID
        INDC(1)=ITX
        ITY=ISUM-ITX+IOF
        INDC(2)=ITY
*----
*  Read angle
*----
        IANG=IANG+1
        CALL XELTSA(NDIM  ,ABSC  ,INDC  ,DNSANG(IANG) ,DDAX)
        DO IDIR=1,NDIM
          DANGLT(IDIR,IANG,1)=DDAX(IDIR)
        ENDDO
        DENLIN=DNSANG(IANG)/RCIRC
        IF(ITX .EQ. 0 .OR. ITY .EQ. 0) THEN
          DNSANG(IANG)= DENUSR
        ELSE
          NTRAC=MAX(1,INT(DENUSR/DENLIN+0.5D0))
          DNSANG(IANG)=DBLE(NTRAC)*DENLIN
        ENDIF
        DDANG(IANG)=DANGLT(1,IANG,1)
        NBSANG(1,IANG)=ITX
        NBSANG(2,IANG)=ITY
        NBSANG(3,IANG)=IPERG
        NBSANG(4,IANG)=1
        NBSANG(5,IANG)=0
        IF((AZMQUA .EQ. 1).AND.(ITX .EQ. ISUM)) THEN
          NBSANG(3,IANG)=IPER(1)
        ELSE IF((AZMQUA .EQ. 1).AND.(ITY .EQ. ISUM)) THEN
          NBSANG(3,IANG)=IPER(2)
        ELSE
*         Find the least common multiple of ITX and ITY
          IA=ITX
          IB=ITY
          DO WHILE (IB.NE.0)
            IC = MOD(IA,IB)
            IA = IB
            IB = IC       
          ENDDO
          IC=ABS(IA)
          NBSANG(1,IANG)=ITX/IC
          NBSANG(2,IANG)=ITY/IC
          NBSANG(4,IANG)=(ITX+ITY)/IC
        ENDIF
        NBSANG(4,IANG)=NBSANG(4,IANG)*NBSANG(3,IANG)
      ENDDO
*----
*  2. Get weights
*----
      CALL XELTCW(NBANGL,DDANG ,DDENWT)
*----
*  3. Compute density and store in adequate location
*----
      DO IANG=1,NBANGL
        DDENWT(IANG,1)=DTWO/DDENWT(IANG,1)
      ENDDO
*----
*  Processing finished: return
*----
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6011)
        DO IANG=1,NBANGL
          IF(DDENWT(IANG,1) .GT. DZERO) THEN
            WRITE(IOUT,6012) IANG,(NBSANG(IDIR,IANG),IDIR=1,4),
     >                      (DANGLT(IDIR,IANG,1),IDIR=1,NDIM),
     >                       DDENWT(IANG,1)
          ENDIF
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
     >       1X,'NANGL =',I8,1X,'NBANGL=',I8)
 6011 FORMAT(' Tracking directions and weights '/
     >       1X,'     Angle',1X,'  Segments',33X,
     >       1X,'     Directions and weight')
 6012 FORMAT(5(1X,I10),4(2X,F24.14))
      END
