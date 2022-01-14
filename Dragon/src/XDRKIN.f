*DECK XDRKIN
      SUBROUTINE XDRKIN(IPTRK, DX, NBX, MLOG)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Construct Bickley tables for KI1(X), KI2(X), KI3(X), KI4(X), KI5(X),
* taking into account logarithmic singularities for KI1(X) and KI2(X).
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* IPTRK   pointer to the LCM table.
* DX      step for tables (here, DX=0.02d0).
* NBX     order for tables (here, NBX=600).
* MLOG    interval for logarithmic singularities (suggested values:
*         MLOG(1)=30, MLOG(2)=15, MLOG(3)= 0, MLOG(4)= 0, MLOG(5)= 0).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
***** NOTE:     YOU MAY UNCOMMENT THE FOLLOWING FORTRAN INSTRUCTION:
*     IMPLICIT NONE
*
***** OUTPUT:   THE FIVE COMMONS OF BICKLEY QUADRATIC TABLES ARE STORED
*               IN THE LCM TABLE AT ADDRESS IPTRK
*
***** CALLS:    *AKIN10* ROUTINE FOR ACCURATE KIN(X) BICKLEY VALUES
*               *AK0BES* ROUTINE FOR ACCURATE K0(X)   BESSEL VALUES
*               *AK1BES* ROUTINE FOR ACCURATE K1(X)   BESSEL VALUES
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK
      INTEGER  NBX, MLOG(5)
      DOUBLE PRECISION DX
*----
*  LOCAL VARIABLES
*----
      INTEGER  I, KI
      DOUBLE PRECISION X, AKIN(-1:10), AK0BES, AK1BES
      DOUBLE PRECISION GAMMA, PIO2, PAS, XLIM
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: C0, C1, C2
      PARAMETER ( GAMMA=0.57721566490153D0,
     >            PIO2= 1.57079632679490D0 )
*----
*  SPACE AND VALUES FOR BICKLEY KI1 TO KI5 TABLES
*----
      INTEGER MKI1,MKI2,MKI3,MKI4,MKI5
      PARAMETER (MKI1=600,MKI2=600,MKI3=600,MKI4=600,MKI5=600)
      REAL    BI1(0:MKI1),BI11(0:MKI1),BI12(0:MKI1),BI2(0:MKI2),
     >        BI21(0:MKI2),BI22(0:MKI2),BI3(0:MKI3),BI31(0:MKI3),
     >        BI32(0:MKI3), BI4(0:MKI4),BI41(0:MKI4),BI42(0:MKI4),
     >        BI5(0:MKI5),BI51(0:MKI5),BI52(0:MKI5)
      REAL    XLIMV(5),PASV(5)
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(C0(5,0:NBX), C1(5,0:NBX), C2(5,0:NBX))
*
      CALL LCMSIX(IPTRK,'FUNC-TABLES',1)
      CALL LCMLEN(IPTRK,'PAS',ILENG,ITYLCM)
      IF(ILENG.GT.0) GO TO 80
*
      IF( MLOG(3).GT.0.OR.MLOG(4).GT.0.OR.MLOG(5).GT.0 )GOTO 99
*
      IF( NBX.NE.   600 )GOTO 97
      IF( DX .NE.0.02D0 )GOTO 98
      PAS=   1.0D0/DX
      XLIM=  DBLE(NBX)*DX
      DO 10 I=1,5
      XLIMV(I)=REAL(XLIM)
      PASV(I)=REAL(PAS)
   10 CONTINUE
*----
*  FIRST, WE CONSTRUCT THE TABLES USING ACCURATE *AKIN10* VALUES
*----
      X= 0.D0
      CALL AKIN10(X,AKIN(1))
      AKIN( 0)= 0.D0
      AKIN(-1)= 0.D0
      DO 30 I= 0, NBX-1
         DO 20 KI= 1, 5
            C2(KI,I)=   0.5D0 * AKIN(KI-2)
            C1(KI,I)= -(AKIN(KI-1)+X*AKIN(KI-2))
            C0(KI,I)=   AKIN(KI)+X*(AKIN(KI-1)+X*C2(KI,I))
   20    CONTINUE
         X= X + DX
         CALL AKIN10(X,AKIN(1))
         AKIN( 0)= AK0BES(X)
         AKIN(-1)= AK1BES(X)
   30 CONTINUE
      DO 40 KI= 1, 5
         C0(KI,NBX)= 0.D0
         C1(KI,NBX)= 0.D0
         C2(KI,NBX)= 0.D0
   40 CONTINUE
*----
*  KI1(X) ADJUSTMENTS
*----
      X= 0.D0
      DO 50 I= 1, MLOG(1)-1
         X= X + DX
         C0(1,I)= C0(1,I) + 0.5D0*X
         C1(1,I)= C1(1,I) - LOG(X)
         C2(1,I)= C2(1,I) - 0.5D0/X
   50 CONTINUE
*
      C1(1,0)= GAMMA-(LOG(2.D0)+1.D0)
*----
*  KI2(X) ADJUSTMENTS
*----
      X= 0.D0
      DO 60 I= 1, MLOG(2)-1
         X= X + DX
         C0(2,I)= C0(2,I) + 0.25D0*X*X
         C1(2,I)= C1(2,I) - X
         C2(2,I)= C2(2,I) + 0.5D0*LOG(X) +0.75D0
   60 CONTINUE
*
      C1(2,0)=  -PIO2
      C2(2,0)=   0.5D0*(LOG(2.D0)+1.5D0-GAMMA)
*----
*  CHARGE LCM
*----
      DO 70 I= 0, NBX
         BI1(I)=  REAL(C0(1,I))
         BI11(I)= REAL(C1(1,I))
         BI12(I)= REAL(C2(1,I))
         BI2(I)=  REAL(C0(2,I))
         BI21(I)= REAL(C1(2,I))
         BI22(I)= REAL(C2(2,I))
         BI3(I)=  REAL(C0(3,I))
         BI31(I)= REAL(C1(3,I))
         BI32(I)= REAL(C2(3,I))
         BI4(I)=  REAL(C0(4,I))
         BI41(I)= REAL(C1(4,I))
         BI42(I)= REAL(C2(4,I))
         BI5(I)=  REAL(C0(5,I))
         BI51(I)= REAL(C1(5,I))
         BI52(I)= REAL(C2(5,I))
   70 CONTINUE
      CALL LCMPUT(IPTRK,'PAS',5,2,PASV)
      CALL LCMPUT(IPTRK,'XLIM',5,2,XLIMV)
      CALL LCMPUT(IPTRK,'MLOG',5,1,MLOG)
      CALL LCMPUT(IPTRK,'BI1',NBX+1,2,BI1(0))
      CALL LCMPUT(IPTRK,'BI11',NBX+1,2,BI11(0))
      CALL LCMPUT(IPTRK,'BI12',NBX+1,2,BI12(0))
      CALL LCMPUT(IPTRK,'BI2',NBX+1,2,BI2(0))
      CALL LCMPUT(IPTRK,'BI21',NBX+1,2,BI21(0))
      CALL LCMPUT(IPTRK,'BI22',NBX+1,2,BI22(0))
      CALL LCMPUT(IPTRK,'BI3',NBX+1,2,BI3(0))
      CALL LCMPUT(IPTRK,'BI31',NBX+1,2,BI31(0))
      CALL LCMPUT(IPTRK,'BI32',NBX+1,2,BI32(0))
      CALL LCMPUT(IPTRK,'BI4',NBX+1,2,BI4(0))
      CALL LCMPUT(IPTRK,'BI41',NBX+1,2,BI41(0))
      CALL LCMPUT(IPTRK,'BI42',NBX+1,2,BI42(0))
      CALL LCMPUT(IPTRK,'BI5',NBX+1,2,BI5(0))
      CALL LCMPUT(IPTRK,'BI51',NBX+1,2,BI51(0))
      CALL LCMPUT(IPTRK,'BI52',NBX+1,2,BI52(0))
   80 CALL LCMSIX(IPTRK,' ',2)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(C2, C1, C0)
      RETURN
*----
*  ERROR SECTION
*----
   97 CALL XABORT('XDRKIN: KIN TABLES HAVE 601 ELEMENTS')
   98 CALL XABORT('XDRKIN: KIN TABLES HAVE A STEP OF 0.02')
   99 CALL XABORT('XDRKIN: NO LOG SINGULARITY TAKEN FOR KI345')
      END
