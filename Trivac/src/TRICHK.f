*DECK TRICHK
      SUBROUTINE TRICHK (HNAMT,IPTRK,IPSYS,IDIM,DIAG,CHEX,IPR,LL4)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Partial consistency check for an ADI-splitted system matrix.
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
* HNAMT   name of the matrix to check.
* IPTRK   L_TRACK pointer to the tracking information.
* IPSYS   L_SYSTEM pointer to system matrices.
* IDIM    number of dimensions.
* DIAG    diagonal symmetry flag for cartesian geometries.
* CHEX    hexagonal geometry flag.
* IPR     perturbation flag (if IPR.ne.0, matrix may contain
*         perturbation values).
* LL4     order of system matrices.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPSYS
      CHARACTER HNAMT*10
      INTEGER IDIM,IPR,LL4
      LOGICAL DIAG,CHEX
*----
*  LOCAL VARIABLES
*----
      PARAMETER (EPSMAX=5.0E-5)
      CHARACTER TEXT10*10,HSMG*60,TEXT8*8
      INTEGER, DIMENSION(:), ALLOCATABLE :: MU,IP
      REAL, DIMENSION(:), ALLOCATABLE :: XTT1
      REAL, DIMENSION(:), POINTER :: A11
      TYPE(C_PTR) A11_PTR
*
      TEXT10=HNAMT(:10)
*----
*  DIMENSION X
*----
      ALLOCATE(XTT1(LL4),MU(LL4),IP(LL4))
      CALL LCMGET(IPTRK,'IPX',IP)
      IF(.NOT.DIAG) THEN
         CALL LCMGET(IPTRK,'MUX',MU)
         CALL LCMGPD(IPSYS,'X_'//TEXT10,A11_PTR)
      ELSE
*        DIAGONAL SYMMETRY
         CALL LCMGET(IPTRK,'MUY',MU)
         CALL LCMGPD(IPSYS,'Y_'//TEXT10,A11_PTR)
      ENDIF
      CALL C_F_POINTER(A11_PTR,A11,(/ MU(LL4) /))
      DO 10 I=1,LL4
      IGAR=MU(IP(I))
      XTT1(I)=A11(IGAR)
      IF((IPR.EQ.0).AND.(XTT1(I).EQ.0.0)) THEN
         WRITE (TEXT8,'(I8)') I
         CALL XABORT('TRICHK: ZERO ELEMENT ON DIAGONAL ELEMENT'//
     1   TEXT8//' OF MATRIX '//TEXT10//'.')
      ENDIF
   10 CONTINUE
      DEALLOCATE(IP,MU)
      IF(IDIM.EQ.1) GO TO 50
*----
*  DIMENSION W
*----
      IF(CHEX) THEN
         ALLOCATE(MU(LL4),IP(LL4))
         CALL LCMGET(IPTRK,'MUW',MU)
         CALL LCMGET(IPTRK,'IPW',IP)
         CALL LCMGPD(IPSYS,'W_'//TEXT10,A11_PTR)
         CALL C_F_POINTER(A11_PTR,A11,(/ MU(LL4) /))
         DO 20 I=1,LL4
         RR=XTT1(I)
         IGAR=MU(IP(I))
         IF(ABS(RR-A11(IGAR)).GT.ABS(RR)*EPSMAX) THEN
            WRITE(HSMG,'(8H: DIAGX(,I6,3H )=,1P,E12.5,7H DIAGW(,I6,
     1      3H )=,E12.5)') I,RR,I,A11(IGAR)
            CALL XABORT('TRICHK: W-AXIS INCONSISTENT ASSEMBLY(1)'//HSMG)
         ENDIF
   20    CONTINUE
         DEALLOCATE(IP,MU)
      ENDIF
*----
*  DIMENSION Y
*----
      ALLOCATE(MU(LL4),IP(LL4))
      CALL LCMGET(IPTRK,'MUY',MU)
      CALL LCMGET(IPTRK,'IPY',IP)
      CALL LCMGPD(IPSYS,'Y_'//TEXT10,A11_PTR)
      CALL C_F_POINTER(A11_PTR,A11,(/ MU(LL4) /))
      DO 30 I=1,LL4
      RR=XTT1(I)
      IGAR=MU(IP(I))
      IF(ABS(RR-A11(IGAR)).GT.ABS(RR)*EPSMAX) THEN
         WRITE(HSMG,'(8H: DIAGX(,I6,3H )=,1P,E12.5,7H DIAGY(,I6,3H )=,
     1   E12.5)') I,RR,I,A11(IGAR)
         CALL XABORT('TRICHK: Y-AXIS INCONSISTENT ASSEMBLY(1)'//HSMG)
      ENDIF
   30 CONTINUE
      DEALLOCATE(IP,MU)
*----
*  DIMENSION Z
*----
      IF(IDIM.GT.2) THEN
         ALLOCATE(MU(LL4),IP(LL4))
         CALL LCMGET(IPTRK,'MUZ',MU)
         CALL LCMGET(IPTRK,'IPZ',IP)
         CALL LCMGPD(IPSYS,'Z_'//TEXT10,A11_PTR)
         CALL C_F_POINTER(A11_PTR,A11,(/ MU(LL4) /))
         DO 40 I=1,LL4
         RR=XTT1(I)
         IGAR=MU(IP(I))
         IF(ABS(RR-A11(IGAR)).GT.ABS(RR)*EPSMAX) THEN
            WRITE(HSMG,'(8H: DIAGX(,I6,3H )=,1P,E12.5,7H DIAGZ(,I6,
     1      3H )=,E12.5)') I,RR,I,A11(IGAR)
            CALL XABORT('TRICHK: Z-AXIS INCONSISTENT ASSEMBLY(1)'//HSMG)
         ENDIF
   40    CONTINUE
         DEALLOCATE(IP,MU)
      ENDIF
   50 DEALLOCATE(XTT1)
      RETURN
      END
