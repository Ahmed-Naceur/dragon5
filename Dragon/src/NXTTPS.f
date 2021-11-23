*DECK NXTTPS
      SUBROUTINE NXTTPS(IPRINT,NPIN  ,IDGPP ,ITSYM ,DRAPIN)
*
*----------
*
*Purpose:
* To test if pins satisfy required symmetry.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s):
* G. Marleau.
*
*Parameters: input
* IPRINT  print level.
* NPIN    number of pins.
* IDGPP   pin direction.
* ITSYM   flag for symmetries to test.
* DRAPIN  pin position/angle/height/radius.
*
*Reference:
*  G. Marleau,
*  New Geometries Processing in DRAGON: The NXT: Module,
*  Report IGE-260, Polytechnique Montreal,
*  Montreal, 2005.
*
*----------
*
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      INTEGER          IPRINT,NPIN,IDGPP,ITSYM(4)
      DOUBLE PRECISION DRAPIN(-1:4,NPIN)
*----
*  Functions
*----
      DOUBLE PRECISION XDRCST,PI
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='NXTTPS')
      DOUBLE PRECISION DCUTOF,DZERO,DONE,DTWO
      PARAMETER       (DCUTOF=1.0D-6,DZERO=0.0D0,DONE=1.0D0,DTWO=2.0D0)
*----
*  Local variables
*----
      INTEGER          IPLOC
      INTEGER          IS,IP,JP,NPIR
      DOUBLE PRECISION DNAP,PIO2,TWOPI
      INTEGER          ITSYR(4)
*----
*  Processing starts:
*  print routine openning output header if required
*  and initialize various parameters.
*----
      IPLOC=IPRINT
      IF(IPLOC .GE. 100) THEN
        WRITE(IOUT,6000) NAMSBR
        IF(IPLOC .GE. 1000) THEN
          WRITE(IOUT,*) 'Symmetry =',ITSYM
          WRITE(IOUT,*) 'Npin =',NPIN
          DO IP=1,NPIN
            WRITE(IOUT,'(6F20.12)') (DRAPIN(JP,IP),JP=-1,4)
          ENDDO
        ENDIF
      ENDIF
*----
*  Rotate symmetry factors for cell directions
*----
      DO IS=1,4
        ITSYR(IS)=ITSYM(IS)
      ENDDO
      IF(IDGPP .EQ. 1) THEN
        ITSYR(1)=ITSYM(2)
        ITSYR(2)=ITSYM(3)
        ITSYR(3)=ITSYM(1)
        ITSYR(4)=ITSYM(4)
        IF(ITSYR(4) .EQ. 1) CALL XABORT(NAMSBR//
     >  ': X=Y symmetry invalid for pin in direction X')
      ELSE IF(IDGPP .EQ. 2) THEN
        ITSYR(1)=ITSYM(3)
        ITSYR(2)=ITSYM(1)
        ITSYR(3)=ITSYM(2)
        ITSYR(4)=ITSYM(4)
        IF(ITSYR(4) .EQ. 1) CALL XABORT(NAMSBR//
     >  ': X=Y symmetry invalid for pin in direction Y')
      ENDIF
      PI=XDRCST('Pi',' ')
      PIO2=PI/DTWO
      TWOPI=DTWO*PI
*----
*  Scan over symmetrization options (in plane only)
*----
      DO 100 IS=1,3
        IF(ITSYR(IS) .EQ. 1) THEN
*----
*  Scan over symmetrized pin
*----
          IF(IPLOC .GE. 1000) THEN
            WRITE(IOUT,'(A3,1X,I8,I8)') 'IS=',IS,ITSYR(IS)
          ENDIF
          DO 110 IP=1,NPIN
*----
*  Find location of pin after symmetrisation
*----
            IF(IS .EQ. 1) THEN
*----
*  X symmetry
*  Symmetric pin should be at \pi-\varphi
*----
              DNAP=PI-DRAPIN(-1,IP)
            ELSE IF(IS .EQ. 2) THEN
*----
*  Y symmetry
*  Symmetric pin should be at -\varphi
*----
              DNAP=-DRAPIN(-1,IP)
            ELSE IF(IS .EQ. 3) THEN
*----
*  Z symmetry
*  Symmetric pin should be at \varphi
*----
              DNAP=DRAPIN(-1,IP)
            ELSE IF(IS .EQ. 4) THEN
*----
*  X=Y symmetry
*  Symmetric pin should be at \pi/2-\varphi
*----
              DNAP=PIO2-DRAPIN(-1,IP)
            ENDIF
*----
*  Position angle in range 0 to 2*Pi
*----
            IF(ABS(DNAP) .LE. DCUTOF) THEN
              DNAP=DZERO
            ELSE IF(DNAP .GT. DCUTOF) THEN
              NPIR=INT((DNAP+DCUTOF)/TWOPI)
              DNAP=DNAP-DBLE(NPIR)*TWOPI
            ELSE
              NPIR=INT((DNAP-DCUTOF)/TWOPI)
              DNAP=DNAP-DBLE(NPIR-1)*TWOPI
            ENDIF
            IF(IPLOC .GE. 1000) THEN
              WRITE(IOUT,'(A3,1X,I8,2F20.12)')
     >        'IP=',IP,DRAPIN(0,IP),DNAP
            ENDIF
            IF(DRAPIN(0,IP) .LT. DCUTOF) THEN
*----
*  For centered pin, test for radial position only
*----
              DO JP=1,NPIN
*----
*  Verify if pin coincide
*  with symmetrized pin
*----
                IF(IPLOC .GE. 1000) THEN
                  WRITE(IOUT,'(A3,1X,I8,3F20.12)')
     >           'JP=',JP,DRAPIN(0,JP),DRAPIN(-1,JP),
     >           ABS(DRAPIN(0,IP)-DRAPIN(0,JP))
                ENDIF
                IF(ABS(DRAPIN(0,IP)-DRAPIN(0,JP)) .LT. DCUTOF)
     >          GO TO 115
              ENDDO
*----
*  no pin coincide, symmetry not satisfied, abort
*----
              CALL XABORT(NAMSBR//': Symmetric pin not found (C)')
            ELSE
*----
*  For pin not centered, test for angular and radial position
*----
              DO JP=1,NPIN
*----
*  Verify if pin coincide
*  with symmetrized pin
*----
                IF(IPLOC .GE. 1000) THEN
                  WRITE(IOUT,'(A3,1X,I8,4F20.8)')
     >           'JP=',JP,DRAPIN(0,JP),DRAPIN(-1,JP),
     >            ABS(DRAPIN(0,IP)-DRAPIN(0,JP)),
     >            ABS(DNAP-DRAPIN(-1,JP))
                ENDIF
                IF(ABS(DRAPIN(0,IP)-DRAPIN(0,JP)) .LT. DCUTOF) THEN
                  IF(ABS(DNAP-DRAPIN(-1,JP)) .LT. DCUTOF) GO TO 115
                ENDIF
              ENDDO
*----
*  no pin coincide, symmetry not satisfied, abort
*----
              CALL XABORT(NAMSBR//': Symmetric pin not found (O-C)')
            ENDIF
 115        CONTINUE
 110      CONTINUE
        ENDIF
 100  CONTINUE
*----
*  Processing finished:
*  print routine closing output header if required
*  and return
*----
      IF(IPLOC .GE. 100) THEN
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
      END
