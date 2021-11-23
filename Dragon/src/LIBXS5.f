*DECK LIBXS5
      SUBROUTINE LIBXS5(IG,NGBIN,IPAP,NSIGF,TT,NTEMPS,TEMPS,DELTF,
     1 SIGTF,SIGAF,SIGFF,DELINF,SGTINF,SGAINF,SGFINF)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Temperature interpolation of autolib (bin cross sections) information.
*
*Copyright:
* Copyright (C) 2014 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IG      coarse energy group under consideration.
* NGBIN   number of coarse energy groups.
* IPAP    APOLIB-XSM pointer.
* NSIGF   number of fine energy groups.
* TT      temperature of isotope.
* NTEMPS  number of tabulated temperatures.
* TEMPS   tabulated temperatures.
*
*Parameters: output
* DELTF   fine group lethargy widths.
* SIGTF   fine group total x-s.
* SIGAF   fine group absorption x-s.
* SIGFF   fine group fission x-s.
* DELINF  calculated lethargy width for group IG.
* SGTINF  calculated infinite-dilution total x-s for group IG.
* SGAINF  calculated infinite-dilution absorption x-s for group IG.
* SGFINF  calculated infinite-dilution fission x-s for group IG.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPAP
      INTEGER IG,NGBIN,NSIGF,NTEMPS
      REAL TT,TEMPS(NTEMPS),DELTF(NSIGF),SIGTF(NSIGF),SIGAF(NSIGF),
     1 SIGFF(NSIGF),DELINF,SGTINF,SGAINF
*----
*  LOCAL VARIABLES
*----
      CHARACTER HSMG*131,TEXT12*12
      PARAMETER (NINT=2,DTMIN=1.0)
      DOUBLE PRECISION D1,D2,D3,D4
      LOGICAL LOK
      REAL, ALLOCATABLE, DIMENSION(:) :: DT,DA,DF,DD
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: SQRTEM,WEIJHT
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(WEIJHT(NTEMPS),SQRTEM(NTEMPS))
*----
*  COMPUTE THE WEIGHTS.
*----
      DO 10 I=1,NTEMPS
      SQRTEM(I)=SQRT(TEMPS(I))
   10 CONTINUE
      IF(NTEMPS.EQ.1) THEN
        IPROX=1
        IGTFIX=1
      ELSE
        STT=SQRT(TT)
        CALL LIBA28(STT,SQRTEM,NTEMPS,NINT,WEIJHT,IORD,IPROX,I0)
        IF(ABS(TT-TEMPS(IPROX)).LE.DTMIN) THEN
          IGTFIX=1
        ELSEIF((STT.LT.SQRTEM(1)).OR.(STT.GT.SQRTEM(NTEMPS))) THEN
          WRITE(HSMG,'(A,F8.2,A,F8.2,A,F8.2)')
     1    'LIBXS5: A TEMPERATURE', TT,'K IS NOT INCLUDED BETWEEN ',
     2    TEMPS(1),' AND ',TEMPS(NTEMPS)
          WRITE(6,'(/1X,A)') HSMG
          IGTFIX=2
        ELSE
          IGTFIX=0
        ENDIF
      ENDIF
*----
*  LOOP OVER TABULATED TEMPERATURES.
*----
      ALLOCATE(DT(NSIGF),DA(NSIGF),DF(NSIGF),DD(NSIGF))
      D1=0.0D0
      CALL XDRSET(SIGTF,NSIGF,0.0)
      CALL XDRSET(SIGAF,NSIGF,0.0)
      CALL XDRSET(SIGFF,NSIGF,0.0)
      DO 50 J=1,IORD
      IT=I0+J
      WRITE(TEXT12,'(6HNTEMPS,I6.6)') IT
      CALL LCMSIX(IPAP,TEXT12,1)
      CALL LCMLEN(IPAP,'DELTF',NV,ITYLCM)
      CALL LCMLEN(IPAP,'SIGFF',NF,ITYLCM)
      IF(NV.NE.NSIGF) CALL XABORT('LIBXS5: INVALID NSIGF.')
      CALL LCMGET(IPAP,'SIGTF',DT)
      CALL LCMGET(IPAP,'SIGAF',DA)
      IF(NF.EQ.NSIGF) CALL LCMGET(IPAP,'SIGFF',DF)
      CALL LCMGET(IPAP,'DELTF',DD)
      CALL LCMSIX(IPAP,' ',2)
      IS=(IT-1)*NGBIN+IG
      IF(IT.EQ.I0+1) THEN
         D1=0.0D0
         DO 20 I=1,NSIGF
         DELTF(I)=DD(I)
         D1=D1+DELTF(I)
   20    CONTINUE
      ELSE
         LOK=.TRUE.
         DO 30 I=1,NSIGF
         LOK=LOK.AND.(DELTF(I).EQ.DD(I))
   30    CONTINUE
         IF(.NOT.LOK) CALL XABORT('LIBXS5: INVALID AUTOLIB MESH.')
      ENDIF
      DO 40 I=1,NSIGF
      SIGTF(I)=SIGTF(I)+REAL(WEIJHT(J)*DT(I))
      SIGAF(I)=SIGAF(I)+REAL(WEIJHT(J)*DA(I))
      IF(NF.EQ.NSIGF) SIGFF(I)=SIGFF(I)+REAL(WEIJHT(J)*DF(I))
   40 CONTINUE
   50 CONTINUE
      D2=0.0D0
      D3=0.0D0
      D4=0.0D0
      DO 60 I=1,NSIGF
      SIGTF(I)=MAX(SIGTF(I),0.0)
      SIGAF(I)=MAX(SIGAF(I),0.0)
      SIGFF(I)=MAX(SIGFF(I),0.0)
      D2=D2+SIGTF(I)*DELTF(I)
      D3=D3+SIGAF(I)*DELTF(I)
      D4=D4+SIGFF(I)*DELTF(I)
   60 CONTINUE
      DELINF=REAL(D1)
      SGTINF=REAL(D2/D1)
      SGAINF=REAL(D3/D1)
      SGFINF=REAL(D4/D1)
      DEALLOCATE(DT,DA,DF,DD)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(SQRTEM,WEIJHT)
      RETURN
      END
