*DECK LIBA23
      SUBROUTINE LIBA23(NG,NANI,TT,NT0,NGTD,NPSN0,TEMP,FGTD,ID2,FAGG,
     1 LAGG,FDGG,WGAL,FAG,LAG,FDG,IAD,DEPL,PSN0,SCAT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembly and temperature interpolation of a transfer matrix stored
* in the APOLIB-2 format.
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
* NANI    anisotropy level. NANI=1 for isotropic scattering.
* TT      temperature of isotope.
* NT0     number of tabulated temperatures.
* NGTD    temperature dependence flag: =0 if no dependence;
*         =NG+1 otherwise.
* NPSN0   size of vector PSN0.
* TEMP    tabulated temperatures.
* FGTD    first temperature-dependent group.
* ID2     number of temperature-dependent terms in the matrix.
* FAGG    first incoming group for the galoche.
* LAGG    last incoming group for the galoche.
* FDGG    first outgoing group for the galoche.
* WGAL    galoche width. The last outgoing group is FDGG+WGAL-1.
* FAG     first incoming group for the rest of the matrix.
* LAG     last incoming group for the rest of the matrix.
* FDG     first outgoing group per incoming group for the rest of
*         the matrix.
* IAD     offset in vector PSN of the data related to each incoming
*         group.
* DEPL    displacement of the IAD offset for the first two
*         temperatures.
* PSN0    input cross section data in APOLIB-2 compressed format.
*
*Parameters: output
* SCAT    interpolated transfer matrix.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NG,NANI,NT0,NGTD,NPSN0,FGTD,ID2,FAGG,LAGG,FDGG,WGAL,FAG,
     1 LAG,FDG(NG),IAD(NG+1),DEPL(NGTD)
      REAL TT,TEMP(NT0),PSN0(NPSN0),SCAT(NG,NG)
*----
*  LOCAL VARIABLES
*----
      CHARACTER HSMG*131
      PARAMETER (NINT=2,DTMIN=1.0)
      LOGICAL LGTP,LGAUX
      DOUBLE PRECISION S
      REAL, ALLOCATABLE, DIMENSION(:) :: PSN
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DTEMP,WEIJHT
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(DTEMP(NT0),WEIJHT(NT0))
*
      NPSN=IAD(NG+1)-1
      IF(NT0.EQ.1) THEN
        IPROX=1
        IGTFIX=1
      ELSE
        DO 10 I=1,NT0
        DTEMP(I)=TEMP(I)
   10   CONTINUE
        CALL LIBA28(TT,DTEMP,NT0,NINT,WEIJHT,IORD,IPROX,I0)
        IF(ABS(TT-TEMP(IPROX)).LE.DTMIN) THEN
          IGTFIX=1
        ELSE IF((TT.LT.TEMP(1)).OR.(TT.GT.TEMP(NT0))) THEN
          WRITE(HSMG,'(A,F8.2,A,F8.2,A,F8.2)')
     1    'LIBA23: A TEMPERATURE', TT,'K IS NOT INCLUDED BETWEEN ',
     2    TEMP(1),' AND ',TEMP(NT0)
          WRITE(6,'(/1X,A)') HSMG
          IGTFIX=2
        ELSE
          IGTFIX=0
        ENDIF
      ENDIF
      ALLOCATE(PSN(NPSN))
      LGTP=I0.GT.0
*----
*  GALOCHE
*----
      IF(WGAL.NE.0) THEN
        DO 15 I=1,WGAL*(LAGG+1-FAGG)
        PSN(I)=PSN0(I)
   15   CONTINUE
      ENDIF
      DO 50 IGA=FAG,LAG
       IPGD=FDG(IGA)
       IDGD=IPGD+IAD(IGA+1)-IAD(IGA)-1
*----
*  PART INDEPENDENT OF TEMPERATURE OF LENGTH LONG FROM IPGD TO IGD
*----
       IF(IPGD.LT.FGTD)THEN
        IGD=MIN0(IDGD,FGTD-1)
        LONG=IGD+1-IPGD
        DO 20 I=1,LONG
        PSN(IAD(IGA)+I-1)=PSN0(IAD(IGA)+I-1)
   20   CONTINUE
       ELSE
        IGD=IPGD-1
        LONG=0
       ENDIF
       IF(IGD.LT.IDGD)THEN
        LONT=IDGD-IGD
*----
*  PART DEPENDENT OF TEMPERATURE
*----
        DO 40 IG=1,LONT
         ID=IAD(IGA)+LONG+IG-1
         ID0=ID
         IDP=ID
         IF(IPROX.GT.1)IDP=IDP+DEPL(IGA)+ID2*(IPROX-2)
         IF(IGTFIX .EQ. 1) THEN
          PSN(ID0)=PSN0(IDP)
         ELSE
          S=0.0D0
          IF(LGTP)ID=ID+DEPL(IGA)+ID2*(I0-1)
          SP=PSN0(IDP)
          LGAUX=.NOT.LGTP
          DO 30 J=1,IORD
           S=S+PSN0(ID)*WEIJHT(J)
           IF(LGAUX)THEN
            ID=ID+DEPL(IGA)
            LGAUX=.FALSE.
           ELSE
            ID=ID+ID2
           ENDIF
   30     CONTINUE
          IF(IGTFIX.EQ.2) THEN
            IF(SP.GE.0.) THEN
               S=MAX(0.D0,S)
            ELSE
               S=MIN(S,0.D0)
            ENDIF
          ENDIF
          PSN(ID0)=REAL(S)
         ENDIF
   40   CONTINUE
       ENDIF
   50 CONTINUE
*----
*  BUILD THE COMPLETE TRANSFER MATRIX SCAT(IG->JG).
*----
      DO 70 IG=1,NG
      DO 60 JG=1,NG
      RAUX=0.
      IF((JG.GE.FAGG).AND.(JG.LE.LAGG).AND.
     1    (IG.GE.FDGG).AND.(IG.LE.(FDGG+WGAL-1))) THEN
         RAUX=PSN((JG-FAGG)*WGAL+IG-FDGG+1)
      ELSE
         IF((IG.GE.FDG(JG)) .AND.
     1       (IG.LE.(IAD(JG+1)-IAD(JG)+FDG(JG)-1))
     2       .AND.(JG.GE.FAG).AND.(JG.LE.LAG))
     3       RAUX=PSN(IAD(JG)+IG-FDG(JG))
      ENDIF
      SCAT(JG,IG)=RAUX/REAL(2*NANI-1)
  60  CONTINUE
  70  CONTINUE
      DEALLOCATE(PSN)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(WEIJHT,DTEMP)
      RETURN
      END
