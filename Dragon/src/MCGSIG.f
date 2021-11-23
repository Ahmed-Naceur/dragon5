*DECK MCGSIG
      SUBROUTINE MCGSIG(IPTRK,NMAT,NGEFF,NALBP,KPSYS,SIGAL,LVOID)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Construct total cross sections and albedos array and check for void.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Le Tellier
*
*Parameters: input
* IPTRK   pointer to the tracking (L_TRACK signature).
* NMAT    number of mixtures.
* NGEFF   effective number of energy groups.
* NALBP   number of physical albedos.
* KPSYS   pointer array for each group properties.
*
*Parameters: output
* SIGAL   total cross sections and albedos array.
* LVOID   void flag.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*---
* SUBROUTINES ARGUMENTS
*---
      TYPE(C_PTR) IPTRK,KPSYS(NGEFF)
      INTEGER NMAT,NGEFF,NALBP
      REAL SIGAL(-6:NMAT,NGEFF)
      LOGICAL LVOID
*---
* LOCAL VARIABLES
*---
      TYPE(C_PTR) JPSYS
      INTEGER I,II,ISA,ICODE(6)
      REAL ALBG(6),ALBEDO(6)
      REAL, ALLOCATABLE, DIMENSION(:) :: ALBP
*---
* RECOVER ALBEDO INFORMATION FROM TRACKING
*---
      CALL LCMGET(IPTRK,'ICODE',ICODE)
      CALL LCMGET(IPTRK,'ALBEDO',ALBG)
*
      LVOID=.FALSE.
      ALLOCATE(ALBP(NALBP))
      DO II=1,NGEFF
         JPSYS=KPSYS(II)
         DO ISA=1,6
           ALBEDO(ISA)=ALBG(ISA)
         ENDDO
         IF(NALBP .GT. 0) THEN
           CALL LCMGET(JPSYS,'ALBEDO',ALBP)
           DO ISA=1,6
             IF(ICODE(ISA).GT.0) ALBEDO(ISA)=ALBP(ICODE(ISA))
           ENDDO
         ENDIF
         CALL LCMGET(JPSYS,'DRAGON-TXSC',SIGAL(0,II))
         DO I=1,NMAT
            IF (SIGAL(I,II).EQ.0.0) LVOID=.TRUE.
         ENDDO
         DO ISA=-6,-1
            SIGAL(ISA,II)=ALBEDO(-ISA)
         ENDDO
      ENDDO
      DEALLOCATE(ALBP)
*
      RETURN
      END
