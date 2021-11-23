*DECK SYBUQV
      SUBROUTINE SYBUQV(ZZR,ZZI,NSURF,NREG,SIGT,MNA,LGFULL,PIJ)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the one-group collision, DP-0 leakage and DP-0 transmission
* probabilities in a Cartesian or hexagonal cell.
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
* NSURF   number of surfaces (4 or 6).
* NREG    number of regions.
* SIGT    total macroscopic cross section.
* MNA     number of angles.
* LGFULL  void flad (=.TRUE. in voided regions).
*
*Parameters: output
* PIJ     collision, DP-0 leakage and DP-0 transmission probabilities
*         in lower triangular form.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER   ZZI(*),NSURF,NREG,MNA
      REAL      ZZR(*),SIGT(NSURF:NSURF+NREG-1),
     >          PIJ(0:(NREG+NSURF)*(NREG+NSURF+1)/2-1)
      LOGICAL   LGFULL(NSURF:NSURF+NREG-1)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (MKI3=600,ZI30=0.785398164)
      LOGICAL   LGEMPT
      REAL, ALLOCATABLE, DIMENSION(:) :: POPTX
      COMMON /BICKL3/BI3(0:MKI3),BI31(0:MKI3),BI32(0:MKI3),PAS3,XLIM3,L3
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(POPTX(2*NREG+4))
*----
*  CHECK FOR VOIDED REGIONS
*----
      LGEMPT=.FALSE.
      DO 10 IR=1,NREG
      LGEMPT=LGEMPT.OR.(.NOT.LGFULL(IR+NSURF-1))
   10 CONTINUE
*----
*  INITIALIZATION
*----
      IZ0=0
      IZR=4
      CALL XDRSET(PIJ(0),(NREG+NSURF)*(NREG+NSURF+1)/2,0.0)
      IF(LGEMPT) THEN
*        VOIDED REGION DETECTED.
         DO 52 IA=1,MNA
         MNT=ZZI(IZ0+1)
         IZ0=IZ0+2
         IZR=IZR+NSURF*2+1
*
*        LOOP OVER COMPLETE TRACKS.
         DO 51 ITT=1,MNT
         IZ0=IZ0+1
         NH=ZZI(IZ0)
         IZ0=IZ0+1
         NX=ZZI(IZ0)
*
         IZ0=IZ0+1
         ISW=ZZI(IZ0) ! WEST SIDE
         JZI=IZ0
         DO 50 ITX=1,NX
         IZR=IZR+1
         WITRAJ=ZZR(IZR)/4.0
         IZ0=JZI
*
*        COMPUTE OPTICAL LENGTHS.
         DO 20 IX=1,NH
         IR=ZZI(IZ0+IX)
         POPTX(IX)=ZZR(IZR+IX)*SIGT(IR)
   20    CONTINUE
         ISE=ZZI(IZ0+NH+1) ! EAST SIDE
         IZR=IZR+NH
*
*        EXTERNAL IS IR=0
         IRS=ISW
         POP1=0.
         DO 40 IXI=1,NH
         IRT=IRS
         IRS=ZZI(IZ0+IXI)
*
*        BEGINNING WITH IRS == JRS, P(i,i)=Taup  -(B3(0)-P3)
         IIJ=IRS*(IRS+3)/2
         CALL SYB33C(PPLUS,POPTX(IXI))
         PIJ(IIJ)=PIJ(IIJ)+2.0*PPLUS*WITRAJ
         TAUX=0.
         DO 30 IXJ=IXI,NH
         JRS=ZZI(IZ0+IXJ)
         IIJ=MAX(IRT*(IRT+1)/2+JRS,JRS*(JRS+1)/2+IRT)
         IF(IRT.LT.NSURF) THEN
            CALL SYB32C(PPLUS,TAUX,POPTX(IXJ),2)
            PIJ(IIJ)=PIJ(IIJ)+PPLUS*WITRAJ
         ELSE
            CALL SYB31C(PPLUS,TAUX,POP1,POPTX(IXJ))
            IF(JRS.EQ.IRT) PPLUS=2.0*PPLUS
            PIJ(IIJ)=PIJ(IIJ)+PPLUS*WITRAJ
         ENDIF
         TAUX=TAUX+POPTX(IXJ)
   30    CONTINUE
         IIJ=MAX(IRT*(IRT+1)/2+ISE,ISE*(ISE+1)/2+IRT)
         IF(IRT.LT.NSURF) THEN
            PIJ(IIJ)=PIJ(IIJ)+TABKI(3,TAUX)*WITRAJ
         ELSE
            CALL SYB32C(PPLUS,TAUX,POP1,2)
            PIJ(IIJ)=PIJ(IIJ)+PPLUS*WITRAJ
         ENDIF
         POP1=POPTX(IXI)
   40    CONTINUE
         IZ0=IZ0+NH+1
*
*        COMPUTE REMAINING PSI FROM LAST REGION TO EAST SIDE.
         IIJ=MAX(IRS*(IRS+1)/2+ISE,ISE*(ISE+1)/2+IRS)
         CALL SYB32C(PPLUS,0.0,POP1,2)
         PIJ(IIJ)=PIJ(IIJ)+PPLUS*WITRAJ
   50    CONTINUE
   51    CONTINUE
   52    CONTINUE
      ELSE
*        NO VOIDED REGION DETECTED. FAST INTEGRATION METHOD.
*        LOOP OVER THE ANGLE(FROM 0 TO PI/3).
         DO 92 IA=1,MNA
         MNT=ZZI(IZ0+1)
         IZ0=IZ0+2
         IZR=IZR+NSURF*2+1
*
*        LOOP OVER COMPLETE TRACKS.
         DO 91 ITT=1,MNT
         IZ0=IZ0+1
         NH=ZZI(IZ0)
         IZ0=IZ0+1
         NX=ZZI(IZ0)
*
*        WEST SIDE
         IZ0=IZ0+1
         ISW=ZZI(IZ0)
         JZI=IZ0
         DO 90 ITX=1,NX
         IZR=IZR+1
         WITRAJ=ZZR(IZR)/4.0
         JZR=IZR
         IZ0=JZI
         IF(WITRAJ.LT.0.0) CALL XABORT('SYBUQV: FAILURE 3.')
         IF(NH.GT.2*NREG+4) CALL XABORT('SYBUQV: FAILURE 4.')
*
*        COMPUTE OPTICAL LENGTHS.
         DO 60 IX=1,NH
         IZ0=IZ0+1
         IR=ZZI(IZ0)
         IZR=IZR+1
         SEG=ZZR(IZR)
         POPTX(IX)=SEG*SIGT(IR)
   60    CONTINUE
*
*        EXTERNAL IS IR=0
         IZ0=JZI
         IRS=ISW
         JZ0=0
         DO 80 IXI=1,NH
         IRT=IRS
         IZ00=IZ0
         IZ0=IZ0+1
         IRS=ZZI(IZ0)
*
*        PREVIOUS REGION J
         JZ0=IZ00
         PKI3=ZI30*WITRAJ
         TAUX=0.
*
*        BEGINNING WITH IRS == JRS,P(i,i)=Taup  -(B3(0)-P3)
         IIJ=IRS*(IRS+3)/2
         PIJ(IIJ)=PIJ(IIJ)+POPTX(IXI)*WITRAJ
         DO 70 IXJ=IXI,NH
         JZ0=JZ0+1
         JRS=ZZI(JZ0)
         QIJ3=PKI3
         TAUX=TAUX+POPTX(IXJ)
         IF(TAUX.GE.XLIM3) THEN
            PKI3=0.0
         ELSE
            K=NINT(TAUX*PAS3)
            PKI3=(BI3(K)+TAUX*(BI31(K)+TAUX*BI32(K)))*WITRAJ
         ENDIF
         SIJ3=QIJ3-PKI3
*
*        BEGINNING OF PIJ : REGION IRS,JRS
         IIJ=MAX(IRS*(IRS+1)/2+JRS,JRS*(JRS+1)/2+IRS)
         PIJ(IIJ)=PIJ(IIJ)-SIJ3
*
*        REMAINING OF PIJ : Region IRS,JRS
*        OR ... PSI IFF IXI=1, IRT IS THE WEST SIDE
         IIJ=MAX(IRT*(IRT+1)/2+JRS,JRS*(JRS+1)/2+IRT)
         PIJ(IIJ)=PIJ(IIJ)+SIJ3
   70    CONTINUE
*
*        COMPUTE LEAKAGE AND TRANSMISSION PROBABILITIES PSI, PSS
         JZ0=JZ0+1
         ISE=ZZI(JZ0)
         IIJ=IRS*(IRS+1)/2+ISE
         PIJ(IIJ)=PIJ(IIJ)-PKI3
         IIJ=MAX(IRT*(IRT+1)/2+ISE,ISE*(ISE+1)/2+IRT)
         PIJ(IIJ)=PIJ(IIJ)+PKI3
   80    CONTINUE
*
*        COMPUTE REMAINING PSI FROM LAST REGION TO EAST SIDE.
         IZ0=IZ0+1
         ISE=ZZI(JZ0)
         IIJ=MAX(IRS*(IRS+1)/2+ISE,ISE*(ISE+1)/2+IRS)
         PIJ(IIJ)=PIJ(IIJ)+ZI30*WITRAJ
   90    CONTINUE
   91    CONTINUE
   92    CONTINUE
         IIJ=NSURF*(NSURF+1)/2-1
         DO 100 JR=1,NREG
         IIJ=IIJ+NSURF+JR
         PIJ(IIJ)=PIJ(IIJ)*2.0
  100    CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(POPTX)
      RETURN
      END
