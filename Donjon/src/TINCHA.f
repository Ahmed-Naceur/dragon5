*DECK TINCHA
      SUBROUTINE TINCHA(IPMAP,NCH,IMPX,NAMCHA,TTIME,RFCHAN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute 'REF-CHAN' record in L_MAP object for history-based cases.
*
*Copyright:
* Copyright (C) 2013 Ecole Polytechnique de Montreal
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* IPMAP   pointer to fuel-map information.
* NCH     number of channels
* IMPX    print flag
* NAMCHA  channel name
* TTIME   refuelling time
*
*Parameters: output
* RFCHAN  time values at which channels are refueled inside a refueling
*         time period
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAP
      INTEGER NCH,IMPX
      CHARACTER*(*) NAMCHA
      REAL RFCHAN(NCH)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40)
      INTEGER ISTATE(NSTATE)
      CHARACTER XNAM*4,YNAM*4,TEXT4*4
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MIX,IXN,IYN
*
      CALL LCMSIX(IPMAP,'GEOMAP',1)
      CALL LCMGET(IPMAP,'STATE-VECTOR',ISTATE)
      IF(ISTATE(1).NE.7) CALL XABORT('TINCHA: 3-D CARTESIAN GEOMETRY'
     +    //' REQUIRED')
      NX = ISTATE(3)
      NY = ISTATE(4)
      NREG = ISTATE(6)
      ALLOCATE(MIX(NREG),IXN(NX),IYN(NY))
      CALL LCMGET(IPMAP,'MIX',MIX)
      CALL LCMSIX(IPMAP,' ',2)
      CALL LCMGET(IPMAP,'XNAME',IXN)
      CALL LCMGET(IPMAP,'YNAME',IYN)
      TEXT4 = NAMCHA(2:3)
      IX = 1
      IY = 1
      DO 10 I=1,NX
        WRITE(XNAM,'(A4)') IXN(I)
        IF (XNAM.EQ.TEXT4) THEN
           IX = I
           GOTO 20
        ENDIF
  10  CONTINUE
  20  TEXT4 = NAMCHA(1:1)
      DO 30 I=1,NY
        WRITE(YNAM,'(A4)') IYN(I)
        IF (YNAM.EQ.TEXT4) THEN
           IY = I
           GOTO 40
        ENDIF
  30  CONTINUE
*
  40  I = (IY-1)*NX + IX
      ICHANAM = MIX(I)
      IF(ICHANAM.EQ.0) CALL XABORT('TINCHA: WRONG CHANNEL NAME')
      DEALLOCATE(IYN,IXN,MIX)
      RFCHAN(ICHANAM) = TTIME
      IF(IMPX.GT.0) THEN
        WRITE(6,*) 'TINCHA: REFUEL ',NAMCHA,' NUMBER ',I,' AT TIME ',
     1  TTIME
      ENDIF
      RETURN
      END
