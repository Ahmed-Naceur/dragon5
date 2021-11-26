*DECK SYBUQ0
      SUBROUTINE SYBUQ0(ZZR,ZZI,NSURF,NREG,SIGT,MNA,LGFULL,PSI,QSS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the one-group DP-1 leakage and transmission
* probabilities in a sectorized Cartesian or hexagonal cell.
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
* ZZR     real tracking elements.
* ZZI     integer tracking elements.
* NSURF   number of surfaces.
* NREG    number of regions.
* SIGT    total macroscopic cross section.
* MNA     number of angles.
* LGFULL  voided region flag (=.TRUE. in voided regions).
*
*Parameters: output
* PSI     escape probabilities.
* QSS     transmission probabilities.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER ZZI(*),NSURF,NREG,MNA
      REAL ZZR(*),SIGT(NREG),PSI(0:NSURF-1,3,NREG),
     1 QSS(9*NSURF*(NSURF-1)/2)
      LOGICAL LGFULL(NREG)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (MKI3=600,MKI4=600,MKI5=600)
      PARAMETER (ZI30=0.785398164,ZI40=0.666666667)
      LOGICAL LGEMPT
      REAL, ALLOCATABLE, DIMENSION(:) :: POPTX
      REAL, ALLOCATABLE, DIMENSION(:,:) :: COSINU
      COMMON /BICKL3/BI3(0:MKI3),BI31(0:MKI3),BI32(0:MKI3),PAS3,XLIM3,L3
      COMMON /BICKL4/BI4(0:MKI4),BI41(0:MKI4),BI42(0:MKI4),PAS4,XLIM4,L4
      COMMON /BICKL5/BI5(0:MKI5),BI51(0:MKI5),BI52(0:MKI5),PAS5,XLIM5,L5
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(POPTX(2*NREG+4),COSINU(2,0:NSURF-1))
*----
*  CHECK FOR VOIDED REGIONS
*----
      LGEMPT=.FALSE.
      DO 10 IR=1,NREG
      LGEMPT=LGEMPT.OR.(.NOT.LGFULL(IR))
   10 CONTINUE
*----
*  INITIALIZATION
*----
      IZ0=0
      IZR=4
      DO 22 I=1,NREG
      DO 21 J=1,3
      DO 20 K=0,NSURF-1
      PSI(K,J,I)=0.0
   20 CONTINUE
   21 CONTINUE
   22 CONTINUE
      DO 30 I=1,9*NSURF*(NSURF-1)/2
      QSS(I)=0.0
   30 CONTINUE
*
      IF(LGEMPT) THEN
*        VOIDED REGION DETECTED.
*        LOOP OVER THE ANGLE (FROM 0 TO PI/3).
         DO 100 IA=1,MNA
         MNT=ZZI(IZ0+1)
         IZ0=IZ0+2
         IZFINI=IZ0
         IOF=IZR
         DO 40 I=0,NSURF-1
         IOF=IOF+1
         COSINU(1,I)=ZZR(IOF)
         IOF=IOF+1
         COSINU(2,I)=ZZR(IOF)
   40    CONTINUE
         IZR=IZR+NSURF*2+1
*
*        LOOP OVER COMPLETE TRACKS.
         DO 90 ITT=1,MNT
         IZ0=IZ0+1
         NH=ZZI(IZ0)
         IZ0=IZ0+1
         NX=ZZI(IZ0)
*
*        WEST SIDE.
         IZ0=IZ0+1
         ISW=ZZI(IZ0)
         IZDEBU=IZ0
         IZFINI=IZ0+NH+1
*
         DO 80 ITX=1,NX
         IZR=IZR+1
         WITRAJ=ZZR(IZR)
         IZ0=IZDEBU
*
*        COMPUTE OPTICAL LENGTHS.
         DO 50 IX=1,NH
         IZ0=IZ0+1
         IR=ZZI(IZ0)+1-NSURF
         IZR=IZR+1
         SEG=ZZR(IZR)
         POPTX(IX)=SEG*SIGT(IR)
   50    CONTINUE
*
         IZ0=IZ0+1
         ISE=ZZI(IZ0)
*
*        EXTERNAL IS IXI=0
         IZ0=IZDEBU
         PKI3=ZI30*WITRAJ
         PKI4=ZI40*WITRAJ
         TAUX=0.
         DO 60 IXI=1,NH
         QIJ3=PKI3
         QIJ4=PKI4
         IZ0=IZ0+1
         IRS=ZZI(IZ0)+1-NSURF
         TAUP=TAUX
         TAUX=TAUX+POPTX(IXI)
         IF(TAUX.GE.XLIM3) THEN
            PKI3=0.0
            PKI4=0.0
         ELSE
            K=NINT(TAUX*PAS3)
            PKI3=(BI3(K)+TAUX*(BI31(K)+TAUX*BI32(K)))*WITRAJ
            PKI4=(BI4(K)+TAUX*(BI41(K)+TAUX*BI42(K)))*WITRAJ
         ENDIF
         SIJ3=QIJ3-PKI3
         SIJ4=QIJ4-PKI4
*
*        COMPUTE LEAKAGE PROBABILITIES.
         IF(.NOT. LGFULL(IRS)) THEN
            CALL SYB32C(SIJ3,TAUP,POPTX(IXI),2)
            SIJ3=SIJ3*WITRAJ
            CALL SYB32C(SIJ4,TAUP,POPTX(IXI),3)
            SIJ4=SIJ4*WITRAJ
         ENDIF
         PSI(ISW,1,IRS)=PSI(ISW,1,IRS)+SIJ3
         PSI(ISW,2,IRS)=PSI(ISW,2,IRS)+SIJ4*COSINU(2,ISW)
         PSI(ISW,3,IRS)=PSI(ISW,3,IRS)+SIJ4*COSINU(1,ISW)
   60    CONTINUE
*
*        EXTERNAL IS IXI=NH
         IZ0=IZFINI
         PKI3=ZI30*WITRAJ
         PKI4=ZI40*WITRAJ
         TAUX=0.
         DO 70 IXI=NH,1,-1
         QIJ3=PKI3
         QIJ4=PKI4
         IZ0=IZ0-1
         IRS=ZZI(IZ0)+1-NSURF
         TAUP=TAUX
         TAUX=TAUX+POPTX(IXI)
         IF(TAUX.GE.XLIM3) THEN
            PKI3=0.0
            PKI4=0.0
         ELSE
            K=NINT(TAUX*PAS3)
            PKI3=(BI3(K)+TAUX*(BI31(K)+TAUX*BI32(K)))*WITRAJ
            PKI4=(BI4(K)+TAUX*(BI41(K)+TAUX*BI42(K)))*WITRAJ
         ENDIF
         SIJ3=QIJ3-PKI3
         SIJ4=QIJ4-PKI4
         IF(.NOT. LGFULL(IRS)) THEN
            CALL SYB32C(SIJ3,TAUP,POPTX(IXI),2)
            SIJ3=SIJ3*WITRAJ
            CALL SYB32C(SIJ4,TAUP,POPTX(IXI),3)
            SIJ4=SIJ4*WITRAJ
         ENDIF
         PSI(ISE,1,IRS)=PSI(ISE,1,IRS)+SIJ3
         PSI(ISE,2,IRS)=PSI(ISE,2,IRS)+SIJ4*COSINU(2,ISE)
         PSI(ISE,3,IRS)=PSI(ISE,3,IRS)+SIJ4*COSINU(1,ISE)
   70    CONTINUE
*
         IF(TAUX.GE.XLIM5) THEN
            PKI5=0.0
         ELSE
            K=NINT(TAUX*PAS5)
            PKI5=(BI5(K)+TAUX*(BI51(K)+TAUX*BI52(K)))*WITRAJ
         ENDIF
*
*        COMPUTE TRANSMISSION PROBABILITIES.
         IIJ=0
         IS1=0
         IS2=0
         IF(ISW .LT. ISE) THEN
            IIJ=((ISE-1)*ISE)/2+ISW
            IS1=ISW
            IS2=ISE
         ELSE IF(ISE .LT. ISW) THEN
            IIJ=((ISW-1)*ISW)/2+ISE
            IS2=ISW
            IS1=ISE
         ELSE
            CALL XABORT('SYBUQ0: IDENTICAL INCOMING AND OUTCOMING SUR'
     1      //'FACES(1).')
         ENDIF
         IIJ=IIJ*9
         QSS(IIJ+1)=QSS(IIJ+1)+PKI3
         QSS(IIJ+2)=QSS(IIJ+2)+PKI4*COSINU(2,IS1)
         QSS(IIJ+3)=QSS(IIJ+3)+PKI4*COSINU(1,IS1)
         QSS(IIJ+4)=QSS(IIJ+4)+PKI4*COSINU(2,IS2)
         QSS(IIJ+7)=QSS(IIJ+7)+PKI4*COSINU(1,IS2)
         QSS(IIJ+5)=QSS(IIJ+5)+PKI5*COSINU(2,IS1)*COSINU(2,IS2)
         QSS(IIJ+6)=QSS(IIJ+6)+PKI5*COSINU(1,IS1)*COSINU(2,IS2)
         QSS(IIJ+8)=QSS(IIJ+8)+PKI5*COSINU(2,IS1)*COSINU(1,IS2)
         QSS(IIJ+9)=QSS(IIJ+9)+PKI5*COSINU(1,IS1)*COSINU(1,IS2)
*        END OF TRACK.
   80    CONTINUE
         IZ0=IZFINI
   90    CONTINUE
*        END OF ANGLE
  100    CONTINUE
      ELSE
*        NO VOIDED REGION DETECTED. FAST INTEGRATION METHOD.
*        LOOP OVER THE ANGLE (FROM 0 TO PI/3).
         DO 170 IA=1,MNA
         MNT=ZZI(IZ0+1)
         IZ0=IZ0+2
         IZFINI=IZ0
         IOF=IZR
         DO 110 I=0,NSURF-1
         IOF=IOF+1
         COSINU(1,I)=ZZR(IOF)
         IOF=IOF+1
         COSINU(2,I)=ZZR(IOF)
  110    CONTINUE
         IZR=IZR+NSURF*2+1
*
*        LOOP OVER COMPLETE TRACKS.
         DO 160 ITT=1,MNT
         IZ0=IZ0+1
         NH=ZZI(IZ0)
         IZ0=IZ0+1
         NX=ZZI(IZ0)
*
*        WEST SIDE.
         IZ0=IZ0+1
         ISW=ZZI(IZ0)
         IZDEBU=IZ0
         IZFINI=IZ0+NH+1
         DO 150 ITX=1,NX
         IZR=IZR+1
         WITRAJ=ZZR(IZR)
         IZ0=IZDEBU
*
*        COMPUTE OPTICAL LENGTHS.
         DO 120 IX=1,NH
         IZ0=IZ0+1
         IR=ZZI(IZ0)+1-NSURF
         IZR=IZR+1
         SEG=ZZR(IZR)
         POPTX(IX)=SEG*SIGT(IR)
  120    CONTINUE
*
         IZ0=IZ0+1
         ISE=ZZI(IZ0)
*
*        EXTERNAL IS IXI=0
         IZ0=IZDEBU
         PKI3=ZI30*WITRAJ
         PKI4=ZI40*WITRAJ
         TAUX=0.
         DO 130 IXI=1,NH
         QIJ3=PKI3
         QIJ4=PKI4
         IZ0=IZ0+1
         IRS=ZZI(IZ0)+1-NSURF
         TAUX=TAUX+POPTX(IXI)
         IF(TAUX.GE.XLIM3) THEN
            PKI3=0.0
            PKI4=0.0
         ELSE
            K=NINT(TAUX*PAS3)
            PKI3=(BI3(K)+TAUX*(BI31(K)+TAUX*BI32(K)))*WITRAJ
            PKI4=(BI4(K)+TAUX*(BI41(K)+TAUX*BI42(K)))*WITRAJ
         ENDIF
         SIJ3=QIJ3-PKI3
         SIJ4=QIJ4-PKI4
*
*        COMPUTE LEAKAGE PROBABILITIES.
         PSI(ISW,1,IRS)=PSI(ISW,1,IRS)+SIJ3
         PSI(ISW,2,IRS)=PSI(ISW,2,IRS)+SIJ4*COSINU(2,ISW)
         PSI(ISW,3,IRS)=PSI(ISW,3,IRS)+SIJ4*COSINU(1,ISW)
  130    CONTINUE
*
*        EXTERNAL IS IXI=NH
         IZ0=IZFINI
         PKI3=ZI30*WITRAJ
         PKI4=ZI40*WITRAJ
         TAUX=0.
         DO 140 IXI=NH,1,-1
         QIJ3=PKI3
         QIJ4=PKI4
         IZ0=IZ0-1
         IRS=ZZI(IZ0)+1-NSURF
         TAUX=TAUX+POPTX(IXI)
         IF(TAUX.GE.XLIM3) THEN
            PKI3=0.0
            PKI4=0.0
         ELSE
            K=NINT(TAUX*PAS3)
            PKI3=(BI3(K)+TAUX*(BI31(K)+TAUX*BI32(K)))*WITRAJ
            PKI4=(BI4(K)+TAUX*(BI41(K)+TAUX*BI42(K)))*WITRAJ
         ENDIF
         SIJ3=QIJ3-PKI3
         SIJ4=QIJ4-PKI4
         PSI(ISE,1,IRS)=PSI(ISE,1,IRS)+SIJ3
         PSI(ISE,2,IRS)=PSI(ISE,2,IRS)+SIJ4*COSINU(2,ISE)
         PSI(ISE,3,IRS)=PSI(ISE,3,IRS)+SIJ4*COSINU(1,ISE)
  140    CONTINUE
*
         IF(TAUX.GE.XLIM5) THEN
            PKI5=0.0
         ELSE
            K=NINT(TAUX*PAS5)
            PKI5=(BI5(K)+TAUX*(BI51(K)+TAUX*BI52(K)))*WITRAJ
         ENDIF
*
*        COMPUTE TRANSMISSION PROBABILITIES.
         IIJ=0
         IS1=0
         IS2=0
         IF(ISW .LT. ISE) THEN
           IIJ=((ISE-1)*ISE)/2+ISW
           IS1=ISW
           IS2=ISE
         ELSE IF(ISE .LT. ISW) THEN
           IIJ=((ISW-1)*ISW)/2+ISE
           IS2=ISW
           IS1=ISE
         ELSE
            CALL XABORT('SYBUQ0: IDENTICAL INCOMING AND OUTCOMING SUR'
     1      //'FACES(2).')
         ENDIF
         IIJ=IIJ*9
         QSS(IIJ+1)=QSS(IIJ+1)+PKI3
         QSS(IIJ+2)=QSS(IIJ+2)+PKI4*COSINU(2,IS1)
         QSS(IIJ+3)=QSS(IIJ+3)+PKI4*COSINU(1,IS1)
         QSS(IIJ+4)=QSS(IIJ+4)+PKI4*COSINU(2,IS2)
         QSS(IIJ+7)=QSS(IIJ+7)+PKI4*COSINU(1,IS2)
         QSS(IIJ+5)=QSS(IIJ+5)+PKI5*COSINU(2,IS1)*COSINU(2,IS2)
         QSS(IIJ+6)=QSS(IIJ+6)+PKI5*COSINU(1,IS1)*COSINU(2,IS2)
         QSS(IIJ+8)=QSS(IIJ+8)+PKI5*COSINU(2,IS1)*COSINU(1,IS2)
         QSS(IIJ+9)=QSS(IIJ+9)+PKI5*COSINU(1,IS1)*COSINU(1,IS2)
*        END OF TRACK.
  150    CONTINUE
         IZ0=IZFINI
  160    CONTINUE
*        END OF ANGLE.
  170    CONTINUE
      ENDIF
*----
*  NUMERICAL ORTHONORMALIZATION
*----
      Z1=ZZR(1)
      Z2=ZZR(2)
      Z3=ZZR(3)
      Z4=ZZR(4)
      DO 185 IRS=1,NREG
      DO 180 ISE=0,NSURF-1
      DEN0=PSI(ISE,1,IRS)
      PSI(ISE,1,IRS)=0.25*Z1*Z1*DEN0
      PSI(ISE,2,IRS)=0.25*Z1*Z4*PSI(ISE,2,IRS)
      PSI(ISE,3,IRS)=0.25*Z1*(Z2*PSI(ISE,3,IRS)-Z3*DEN0)
  180 CONTINUE
  185 CONTINUE
      IIJ=0
      DO 190 IS=1,NSURF*(NSURF-1)/2
      DEN0=QSS(IIJ+1)
      DEN1=QSS(IIJ+3)
      DEN2=QSS(IIJ+7)
      QSS(IIJ+1)=Z1*Z1*DEN0
      QSS(IIJ+3)=Z1*(Z2*DEN1-Z3*DEN0)
      QSS(IIJ+7)=Z1*(Z2*DEN2-Z3*DEN0)
      QSS(IIJ+9)=Z2*Z2*QSS(IIJ+9)-Z2*Z3*(DEN1+DEN2)+Z3*Z3*DEN0
      DEN1=QSS(IIJ+2)
      DEN2=QSS(IIJ+4)
      QSS(IIJ+2)=Z1*Z4*DEN1
      QSS(IIJ+4)=Z1*Z4*DEN2
      QSS(IIJ+6)=(Z2*QSS(IIJ+6)-Z3*DEN2)*Z4
      QSS(IIJ+8)=(Z2*QSS(IIJ+8)-Z3*DEN1)*Z4
      QSS(IIJ+5)=Z4*Z4*QSS(IIJ+5)
      IIJ=IIJ+9
  190 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(COSINU,POPTX)
      RETURN
      END
