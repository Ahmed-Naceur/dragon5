!
!-----------------------------------------------------------------------
!
!Purpose:
! Return the real spherical harmonics corresponding to a set of
! direction cosines.
!
!Copyright:
! Copyright (C) 2004 Ecole Polytechnique de Montreal
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version
!
!Author(s): A. Hebert
!
!Parameters: input
! L       Legendre order.
! M       azimuthal order.
! ZMU     X-directed direction cosine.
! ETA     Y-directed direction cosine.
! XI      Z-directed direction cosine.
!
!Parameters: output
! PNOUT   value of the spherical harmonics.
!
!-----------------------------------------------------------------------
!
REAL FUNCTION PNSH(L,M,ZMU,ETA,XI) RESULT(PNOUT)
  !----
  !  FUNCTION ARGUMENTS
  !----
  INTEGER L,M
  REAL ZMU,ETA,XI
  !----
  !  LOCAL VARIABLES
  !----
  PARAMETER (PI=3.141592653589793)
  DOUBLE PRECISION FACTOR,PHI,PNL,DZMU,COEF,ALPLGN
  !
  TEST=ZMU*ZMU+ETA*ETA+XI*XI
  IF(ABS(TEST-1.0).GT.1.0E-5) THEN
    CALL XABORT('PNSH: INVALID DIRECTION COSINES.')
  ENDIF
  PNOUT=0.0
  IF((L.EQ.0).AND.(M.EQ.0)) THEN
    PNOUT=1.0
  ELSE IF((L.EQ.1).AND.(M.EQ.-1)) THEN
    PNOUT=XI
  ELSE IF((L.EQ.1).AND.(M.EQ.0)) THEN
    PNOUT=ZMU
  ELSE IF((L.EQ.1).AND.(M.EQ.1)) THEN
    PNOUT=ETA
  ELSE IF((L.EQ.2).AND.(M.EQ.-2)) THEN
    PNOUT=SQRT(3.0)*ETA*XI
  ELSE IF((L.EQ.2).AND.(M.EQ.-1)) THEN
    PNOUT=SQRT(3.0)*ZMU*XI
  ELSE IF((L.EQ.2).AND.(M.EQ.0)) THEN
    PNOUT=0.5*(3.0*ZMU*ZMU-1.0)
  ELSE IF((L.EQ.2).AND.(M.EQ.1)) THEN
    PNOUT=SQRT(3.0)*ZMU*ETA
  ELSE IF((L.EQ.2).AND.(M.EQ.2)) THEN
    PNOUT=0.5*SQRT(3.0)*(ETA*ETA-XI*XI)
  ELSE IF((L.EQ.3).AND.(M.EQ.-3)) THEN
    PNOUT=SQRT(5./8.)*XI*(3.0*ETA*ETA-XI*XI)
  ELSE IF((L.EQ.3).AND.(M.EQ.-2)) THEN
    PNOUT=SQRT(15.0)*ETA*XI*ZMU
  ELSE IF((L.EQ.3).AND.(M.EQ.-1)) THEN
    PNOUT=SQRT(3./8.)*XI*(5.0*ZMU*ZMU-1.0)
  ELSE IF((L.EQ.3).AND.(M.EQ.0)) THEN
    PNOUT=0.5*ZMU*(5.0*ZMU*ZMU-3.0)
  ELSE IF((L.EQ.3).AND.(M.EQ.1)) THEN
    PNOUT=SQRT(3./8.)*ETA*(5.0*ZMU*ZMU-1.0)
  ELSE IF((L.EQ.3).AND.(M.EQ.2)) THEN
    PNOUT=SQRT(15.0/4.0)*ZMU*(ETA*ETA-XI*XI)
  ELSE IF((L.EQ.3).AND.(M.EQ.3)) THEN
    PNOUT=SQRT(5./8.)*ETA*(ETA*ETA-3.0*XI*XI)
  ELSE IF((L.EQ.4).AND.(M.EQ.-4)) THEN
    PNOUT=0.5*SQRT(35.)*ETA*XI*(ETA*ETA-XI*XI)
  ELSE IF((L.EQ.4).AND.(M.EQ.-3)) THEN
    PNOUT=0.5*SQRT(0.5*35.)*ZMU*XI*(3.*ETA*ETA-XI*XI)
  ELSE IF((L.EQ.4).AND.(M.EQ.-2)) THEN
    PNOUT=SQRT(5.)*(21.*ZMU*ZMU-3.)*ETA*XI/6.
  ELSE IF((L.EQ.4).AND.(M.EQ.-1)) THEN
    PNOUT=0.5*SQRT(2.5)*ZMU*XI*(7.*ZMU*ZMU-3.)
  ELSE IF((L.EQ.4).AND.(M.EQ.0)) THEN
    PNOUT=(35.*ZMU**4-30.*ZMU*ZMU+3.)/8.
  ELSE IF((L.EQ.4).AND.(M.EQ.1)) THEN
    PNOUT=0.5*SQRT(2.5)*ZMU*ETA*(7.*ZMU*ZMU-3.)
  ELSE IF((L.EQ.4).AND.(M.EQ.2)) THEN
    PNOUT=SQRT(5.)*(21.*ZMU*ZMU-3.)*(ETA*ETA-XI*XI)/12.
  ELSE IF((L.EQ.4).AND.(M.EQ.3)) THEN
    PNOUT=0.5*SQRT(0.5*35.)*ZMU*ETA*(ETA*ETA-3.*XI*XI)
  ELSE IF((L.EQ.4).AND.(M.EQ.4)) THEN
    PNOUT=SQRT(35.)*(ETA**4-6.*(ETA*XI)**2+XI**4)/8.
  ELSE IF((L.EQ.5).AND.(M.EQ.-5)) THEN
    PNOUT=21.*XI*(5.*ETA**4-10.*(ETA*XI)**2+XI**4)/(8.*SQRT(14.))
  ELSE IF((L.EQ.5).AND.(M.EQ.-4)) THEN
    PNOUT=0.5*105.*ZMU*ETA*XI*(ETA*ETA-XI*XI)/SQRT(35.)
  ELSE IF((L.EQ.5).AND.(M.EQ.-3)) THEN
    PNOUT=35.*(9*ZMU*ZMU-1.)*XI*(3.*ETA*ETA-XI*XI)/(8.*SQRT(70.))
  ELSE IF((L.EQ.5).AND.(M.EQ.-2)) THEN
    PNOUT=0.5*SQRT(105.)*ZMU*(3.*ZMU*ZMU-1.)*ETA*XI
  ELSE IF((L.EQ.5).AND.(M.EQ.-1)) THEN
    PNOUT=SQRT(15.)*XI*(21.*ZMU**4-14.*ZMU*ZMU+1.)/8.
  ELSE IF((L.EQ.5).AND.(M.EQ.0)) THEN
    PNOUT=ZMU*(63.*ZMU**4-70.*ZMU*ZMU+15.)/8.
  ELSE IF((L.EQ.5).AND.(M.EQ.1)) THEN
    PNOUT=SQRT(15.)*ETA*(21.*ZMU**4-14.*ZMU*ZMU+1.)/8.
  ELSE IF((L.EQ.5).AND.(M.EQ.2)) THEN
    PNOUT=0.25*SQRT(105.)*ZMU*(3.*ZMU*ZMU-1.)*(ETA*ETA-XI*XI)
  ELSE IF((L.EQ.5).AND.(M.EQ.3)) THEN
    PNOUT=35.*(9*ZMU*ZMU-1.)*ETA*(ETA*ETA-3.*XI*XI)/(8.*SQRT(70.))
  ELSE IF((L.EQ.5).AND.(M.EQ.4)) THEN
    PNOUT=105.*ZMU*(ETA**4-6.*(ETA*XI)**2+XI**4)/(8.*SQRT(35.))
  ELSE IF((L.EQ.5).AND.(M.EQ.5)) THEN
    PNOUT=21.*ETA*(ETA**4-10.*(ETA*XI)**2+5.*XI**4)/(8.*SQRT(14.))
  ELSE
    FACTOR=SQRT(1.0D0-ZMU*ZMU)
    PHI=0.0D0
    IF(XI.GE.0) THEN
      PHI=ACOS(ETA/FACTOR)
    ELSE IF(XI.LT.0) THEN
      PHI=2.0D0*PI-ACOS(ETA/FACTOR)
    ENDIF
    COEF=SQRT(2.0D0*ALFACT(L-ABS(M))/ALFACT(L+ABS(M)))
    DZMU=ZMU
    PNL=ALPLGN(L,ABS(M),DZMU)
    IF(M.GT.0) THEN
      PNOUT=REAL(COEF*PNL*COS(M*PHI))
    ELSE IF(M.EQ.0) THEN
      PNOUT=REAL(PNL)
    ELSE IF(M.LT.0) THEN
      PNOUT=REAL(COEF*PNL*SIN(-M*PHI))
    ENDIF
  ENDIF
  RETURN
  !
  CONTAINS
  RECURSIVE DOUBLE PRECISION FUNCTION ALFACT(N) RESULT(OUT)
  ! return the factorial of N
  INTEGER N
  IF(N.LE.1) THEN
    OUT=1.0D0
  ELSE
    OUT=N*ALFACT(N-1)
  ENDIF
  END FUNCTION ALFACT
END FUNCTION PNSH
