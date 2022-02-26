*DECK NXTITA
      FUNCTION NXTITA(POSTRI ,PINPOS,VOLINT)
*
*----------
*
*Purpose:
* Compute the volume of intersection between
* a 2--D triangle and an annular pin.
*
*Copyright:
* Copyright (C) 2010 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): G. Marleau.
*
*Parameters: input
* POSTRI  spatial description of the triangle with:
*         POSTRI(1,J) $X$ location of corner $j$ of triangle;
*         POSTRI(2,J) $Y$ location of corner $j$ of triangle.
* PINPOS  spatial description of the annular pin region with:
*         PINPOS(0) the radius of the annular pin;
*         PINPOS(1) the $X$ position of the annular pin center;
*         PINPOS(2) the $Y$ position of the annular pin center.
*
*Parameters: output
* NXTITA  type of intersection between haxagon and annular pin, where:
*         = 0 means that there is no intersection
*         between the two regions;
*         = 1 means that the hexagon
*         is all located inside the annular pin;
*         = 2 means that the annular pin
*         is all located inside the hexagon;
*         =-1 means that the intersection between
*         the hexagon and the annular pin is partial.
* VOLINT  2-D volume of intersection (area) between hexagon and
*         annular pin.
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*
*----
*
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          NXTITA
      DOUBLE PRECISION POSTRI(2,3),PINPOS(0:2)
      DOUBLE PRECISION VOLINT
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTITA')
      INTEGER          IPRINT
      PARAMETER       (IPRINT=100)
      DOUBLE PRECISION DCUTOF
      PARAMETER       (DCUTOF=1.0D-8)
      DOUBLE PRECISION DZERO,DONE,DHALF,DSQ3O2
      PARAMETER       (DZERO=0.0D0,DONE=1.0D0,
     >                 DHALF=0.5D0,DSQ3O2=0.86602540378444D0)
*----
*  Functions
*----
      DOUBLE PRECISION XDRCST,PI
*----
*  Local variables
*----
      INTEGER          IDIR,IFACE,NFPINS,NCIN,NKPINS,ICOR,ICORIN(3),IC
      DOUBLE PRECISION VOLOUT,RADP2,X1,Y1,X2,Y2,VOLTRI,ALPHA
      DOUBLE PRECISION CPPP(2),CORPOS(2),DIRFAC(2,3),DIRLIN(2,3),RADC,
     >                 DIRN,DISTM,DISTC,DNOR(3),DTAN(3),
     >                 VT,AFPINT(3),XFPINT(3),
     >                 YFPINT(3),TC1(2),TC2(2),DISTS,DISTF
*----
*  Print header if required
*----
      IF(IPRINT .GE. 200) THEN
        WRITE(IOUT,6000) NAMSBR
        WRITE(IOUT,6010) (POSTRI(1,ICOR),POSTRI(2,ICOR),ICOR=1,3)
        WRITE(IOUT,6011) (PINPOS(IFACE),IFACE=0,2)
      ENDIF
*----
*  Initialize PI, NXTITA and VOLINT
*----
      PI=XDRCST('Pi',' ')
      NXTITA=0
      VOLINT=DZERO
      VOLOUT=DZERO
*----
*  Evaluate distance from FACES to pin center
*----
      RADP2=PINPOS(0)**2
      NFPINS=0
      NCIN=0
      NKPINS=0
*----
*  Compute volume of triangle
*----
      X1=POSTRI(1,2)-POSTRI(1,1)
      X2=POSTRI(1,3)-POSTRI(1,1)
      Y1=POSTRI(2,2)-POSTRI(2,1)
      Y2=POSTRI(2,3)-POSTRI(2,1)
      VOLTRI=ABS(X1*Y2-X2*Y1)/2.0D0
*----
*  Analyze each face
*----
      CPPP(1)=POSTRI(1,2)
      CPPP(2)=POSTRI(2,2)
      CORPOS(1)=POSTRI(1,3)
      CORPOS(2)=POSTRI(2,3)
      DO IFACE=1,3        
        ICOR=IFACE
        ICORIN(ICOR)=0
*        write(6,'(A6,5x,2i5)') 'FACE  ',IFACE,ICOR
        AFPINT(IFACE)=DZERO
        XFPINT(IFACE)=DZERO
        YFPINT(IFACE)=DZERO
        DNOR(IFACE)=DZERO
        DTAN(IFACE)=DZERO
        RADC=DZERO
        DISTM=DZERO
*----
*  Find direction of face and its normal directed inward
*  the triangle
*----
        DIRLIN(1,IFACE)=(POSTRI(1,IFACE)-CORPOS(1))
        DIRLIN(2,IFACE)=(POSTRI(2,IFACE)-CORPOS(2))
        DIRN=SQRT(DIRLIN(1,IFACE)**2+DIRLIN(2,IFACE)**2)
        DIRLIN(1,IFACE)=DIRLIN(1,IFACE)/DIRN
        DIRLIN(2,IFACE)=DIRLIN(2,IFACE)/DIRN
        DIRFAC(1,IFACE)=-DIRLIN(2,IFACE)
        DIRFAC(2,IFACE)=DIRLIN(1,IFACE)
*        write(6,*) 'Corner and pin position with respect to corner'
*        write(6,'(4F20.10)') 
*     > CORPOS(1),CORPOS(2),PINPOS(1)-CORPOS(1),
*     > PINPOS(2)-CORPOS(2)
*        write(6,*) 'Face tangent and normal'
*        write(6,'(4F20.10)') 
*     > DIRLIN(1,IFACE),DIRLIN(2,IFACE),DIRFAC(1,IFACE),DIRFAC(2,IFACE)
        DO IDIR=1,2
          DISTC=PINPOS(IDIR)-CORPOS(IDIR)
          DISTM=DISTM+(CPPP(IDIR)-CORPOS(IDIR))*DIRFAC(IDIR,IFACE)
          DNOR(IFACE)=DNOR(IFACE)+DISTC*DIRFAC(IDIR,IFACE)
          DTAN(IFACE)=DTAN(IFACE)+DISTC*DIRLIN(IDIR,IFACE)
          RADC=RADC+DISTC**2
        ENDDO
*        write(6,*) 'Distance of center to face',DNOR(IFACE)
        IF(DNOR(IFACE) .GT. DISTM+PINPOS(0)) THEN
          NXTITA=0
          VOLINT=DZERO
          IF(IPRINT .GE. 200) THEN
            WRITE(IOUT,6012) NAMSBR,NXTITA,VOLINT
            WRITE(IOUT,6001) NAMSBR
          ENDIF
          RETURN
        ELSE IF(DNOR(IFACE) .LT. -PINPOS(0)) THEN
          NXTITA=0
          VOLINT=DZERO
          IF(IPRINT .GE. 200) THEN
            WRITE(IOUT,6012) NAMSBR,NXTITA,VOLINT
            WRITE(IOUT,6001) NAMSBR
          ENDIF
          RETURN
        ELSE IF(DNOR(IFACE) .GE. PINPOS(0)) THEN
          NFPINS=NFPINS+1
          NKPINS=NKPINS+1
        ELSE
*----
*  Point of intersection of current face with pin 
*----
          XFPINT(IFACE)=DNOR(IFACE)
          YFPINT(IFACE)=SQRT(RADP2-XFPINT(IFACE)**2)
          AFPINT(IFACE)=ACOS(XFPINT(IFACE)/PINPOS(0))
          VT=XFPINT(IFACE)*YFPINT(IFACE)
          DISTS=DTAN(IFACE)-YFPINT(IFACE)
          DISTF=DTAN(IFACE)+YFPINT(IFACE)
*          write(6,*) 'DISTS/DISTF/DIRN=',DISTS,DISTF,DIRN
          IF(DISTS .LE. DIRN .AND. DISTF .GT. DZERO) THEN
*            write(6,'(A9,2X,7F20.10)') 
*     >      'Volout 1=',XFPINT(IFACE),YFPINT(IFACE),
*     >      AFPINT(IFACE),VT,VOLOUT,
*     >      RADP2*AFPINT(IFACE)-VT,VOLOUT+RADP2*AFPINT(IFACE)-VT
            VOLOUT=VOLOUT+RADP2*AFPINT(IFACE)-VT
          ELSE
            NKPINS=NKPINS+1
*            write(6,'(A9,2X,7F20.10)') 
*     >      'Volout 3=',XFPINT(IFACE),YFPINT(IFACE),
*     >      AFPINT(IFACE),VT,VOLOUT,
*     >      RADP2*AFPINT(IFACE)-VT,VOLOUT
            VOLOUT=VOLOUT
          ENDIF
        ENDIF
        IF(RADC .LT. RADP2) THEN
*----
*  Identify corners in pin
*----
          ICORIN(ICOR)=1
          NCIN=NCIN+1
        ENDIF
        CPPP(1)=CORPOS(1)
        CPPP(2)=CORPOS(2)
        CORPOS(1)=POSTRI(1,IFACE)
        CORPOS(2)=POSTRI(2,IFACE)
      ENDDO
      IF(NFPINS .EQ. 3) THEN
*----
*  Pin completely inside triangle
*----
        NXTITA=2
        VOLINT=PI*RADP2
      ELSE IF(NCIN .EQ. 3) THEN
*----
*  Triangle completely inside pin
*----
        NXTITA=1
        VOLINT=VOLTRI
      ELSE IF(NKPINS .EQ. 3) THEN
*----
*  No intersection between triangle and pin
*----
        NXTITA=0
        VOLINT=DZERO
      ELSE
*----
*  For corners inside pin, find intersection of outside surfaces
*  and remove from VOLOUT
*----
        NXTITA=-1
        DO ICOR=1,3
*          write(6,*) 'ICORIN=',ICOR,ICORIN(ICOR)
          IF(ICORIN(ICOR) .EQ. 1) THEN
*----
*  Point of intersection of previous face with pin in the positive
*  direction
*----
            IFACE=ICOR-1
            IF(IFACE .LE. 0) IFACE=3+IFACE
            IC=ICOR-2
            IF(IC .LE. 0) IC=3+IC
            DO IDIR=1,2
              TC2(IDIR)=POSTRI(IDIR,IC)
     >                 +(DTAN(IFACE)+YFPINT(IFACE))*DIRLIN(IDIR,IFACE)
            ENDDO
*            write(6,*) 'TC2=',TC2(1),TC2(2)
*----
*  Point of intersection of current face with pin in the negative
*  direction
*----
            IFACE=ICOR
            IC=ICOR-1
            IF(IC .LE. 0) IC=3+IC
            DO IDIR=1,2
              TC1(IDIR)=POSTRI(IDIR,IC)
     >                 +(DTAN(IFACE)-YFPINT(IFACE))*DIRLIN(IDIR,IFACE)
            ENDDO
*            write(6,*) 'TC1=',TC1(1),TC1(2)
*----
*  Triangle outside is identified by CORPOS(*,IC),TC1,TC2
*  Compute its volume
*----
            X1=TC1(1)-POSTRI(1,IC)
            X2=TC2(1)-POSTRI(1,IC)
            Y1=TC1(2)-POSTRI(2,IC)
            Y2=TC2(2)-POSTRI(2,IC)
            VOLTRI=ABS(X1*Y2-X2*Y1)/2.0D0
*            write(6,*) 'Triangle 1=',X1,X2,Y1,Y2,VOLTRI
*----
*  Add contribution from annular sector between T1 and T2
*  Compute distance of T1 to T2
*----
            Y2=DZERO
            DO IDIR=1,2
              Y2=Y2+(TC1(IDIR)-TC2(IDIR))**2
            ENDDO
            Y2=Y2/4.0D0
            X2=SQRT(RADP2-Y2)
            Y2=SQRT(Y2)
            ALPHA=ACOS(X2/PINPOS(0))
*            write(6,*) 'Triangle 2=',X2,Y2,ALPHA,
*     >      VOLOUT,RADP2*ALPHA,X2*Y2,RADP2*ALPHA-X2*Y2+VOLTRI
            VOLOUT=VOLOUT-RADP2*ALPHA+X2*Y2-VOLTRI
          ENDIF
        ENDDO
        VOLINT=PI*RADP2-VOLOUT
*        write(6,*) 'Triangle 3=',VOLINT,PI*RADP2,VOLOUT
      ENDIF
      IF(IPRINT .GE. 200) THEN
        WRITE(IOUT,6012) NAMSBR,NXTITA,VOLINT
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6010 FORMAT('POSTRI={',5(F20.10,','),F20.10,'};')
 6011 FORMAT('PINPOS={',2(F20.10,','),F20.10,'};')
 6012 FORMAT(A6,'={',I5,',',F20.10,'};')
      END
