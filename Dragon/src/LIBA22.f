*DECK LIBA22
      SUBROUTINE LIBA22(NG,TT,NT0,NSECT0,FGTD,TEMP,SECT0,SECT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Temperature interpolation of a cross section array stored in the
* APOLIB-2 format.
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
* NG      number of energy groups.
* TT      temperature of isotope.
* NT0     number of tabulated temperatures.
* NSECT0  size of vector SECT0.
* FGTD    first temperature-dependent energy group.
* TEMP    tabulated temperatures.
* SECT0   input cross section data in APOLIB-2 compressed format.
*
*Parameters: output
* SECT    interpolated transfer matrix.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NG,NT0,NSECT0
      REAL TT,TEMP(NT0),SECT0(NSECT0),SECT(NG)
*----
*  LOCAL VARIABLES
*----
      CHARACTER HSMG*131
      PARAMETER (NINT=2,DTMIN=1.0)
      INTEGER FGTD
      DOUBLE PRECISION S
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DTEMP,WEIJHT
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(DTEMP(NT0),WEIJHT(NT0))
*
      IF(NSECT0.EQ.NG) THEN
        DO 10 I=1,NG
        SECT(I)=SECT0(I)
   10   CONTINUE
        RETURN
      ENDIF
*
      DO 15 I=1,NT0
      DTEMP(I)=TEMP(I)
   15 CONTINUE
      IF(NT0.EQ.1) THEN
        IPROX=1
        IGTFIX=1
      ELSE
        CALL LIBA28(TT,DTEMP,NT0,NINT,WEIJHT,IORD,IPROX,I0)
        IF(ABS(TT-TEMP(IPROX)).LE.DTMIN) THEN
          IGTFIX=1
        ELSE IF((TT.LT.TEMP(1)).OR.(TT.GT.TEMP(NT0))) THEN
          WRITE(HSMG,'(A,F8.2,A,F8.2,A,F8.2)')
     1    'LIBA22: A TEMPERATURE', TT,'K IS NOT INCLUDED BETWEEN ',
     2    TEMP(1),' AND ',TEMP(NT0)
          WRITE(6,'(/1X,A)') HSMG
          IGTFIX=2
        ELSE
          IGTFIX=0
        ENDIF
      ENDIF
*
      IDIS=NG+1-FGTD
      IPID=(IPROX-1)*IDIS
      IF(FGTD.GT.1) THEN
        DO 20 I=1,FGTD-1
        SECT(I)=SECT0(I)
   20   CONTINUE
      ENDIF
      IF(IGTFIX.EQ.1) THEN
       ISECT0=FGTD+IPID
       DO 30 I=1,IDIS
       SECT(FGTD+I-1)=SECT0(ISECT0+I-1)
   30  CONTINUE
      ELSE
        DO 50 I=FGTD,NG
        S=0.D0
        ID=I+I0*IDIS
        IDP=I+IPID
        DO 40 J=1,IORD
        S=S+WEIJHT(J)*SECT0(ID)
        ID=ID+IDIS
   40   CONTINUE
        IF(IGTFIX.EQ.2) THEN
          IF(SECT0(IDP).GE.0.) THEN
             S=MAX(0.D0,S)
          ELSE
             S=MIN(S,0.D0)
          ENDIF
        ENDIF
        SECT(I)=REAL(S)
   50   CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(WEIJHT,DTEMP)
      RETURN
      END
