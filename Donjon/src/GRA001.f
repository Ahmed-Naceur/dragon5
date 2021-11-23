*DECK GRA001
      SUBROUTINE GRA001(IPFLX,IPGPT,NVAR,NCST,DERIV)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute new gradients of system characteristics (part 2).
*
*Copyright:
* Copyright (C) 2012 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* IPFLX   pointer of a LCM object containing a set solutions of
*         fixed-source eigenvalue problems.
* IPGPT   pointer of a LCM object containing a set of fixed sources.
* NVAR    number of control variables.
* NCST    number of constraints with indirect effects (can be zero).
*
*Parameters: output
* DERIV   gradient matrix.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPFLX,IPGPT
      INTEGER NVAR,NCST
      DOUBLE PRECISION DERIV(NVAR,NCST+1)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPFLX,KPFLX,JPGPT,KPGPT
      PARAMETER (NSTATE=40)
      INTEGER ISTATE(NSTATE)
      DOUBLE PRECISION SUM
      REAL, ALLOCATABLE, DIMENSION(:) :: DFLUX,SOUR
*
      CALL LCMGET(IPFLX,'STATE-VECTOR',ISTATE)
      NGRP=ISTATE(1)
      NUN=ISTATE(2)
      ITYPE=ISTATE(3)
      NGPT=ISTATE(5)
      CALL LCMGET(IPGPT,'STATE-VECTOR',ISTATE)
      IF(ISTATE(1).NE.NGRP) CALL XABORT('GRA001: INVALID NGRP')
      IF(ISTATE(2).NE.NUN) CALL XABORT('GRA001: INVALID NUN')
      ND=ISTATE(3)
      NA=ISTATE(4)
*
      ALLOCATE(SOUR(NUN),DFLUX(NUN))
      CALL XDDSET(DERIV,NVAR*(NCST+1),0.0D0)
      IF(ITYPE.EQ.100) THEN
*       EXPLICIT APPROACH
        IF(NVAR.NE.NGPT) CALL XABORT('GRA001: INVALID NGPT(1)')
        IF(NCST+1.NE.NA) CALL XABORT('GRA001: INVALID NA(1)')
        JPFLX=LCMGID(IPFLX,'DFLUX')
        JPGPT=LCMGID(IPGPT,'ASOUR')
        DO 25 IVAR=1,NVAR
        KPFLX=LCMGIL(JPFLX,IVAR)
        DO 20 ICST=1,NCST+1
        KPGPT=LCMGIL(JPGPT,ICST)
        SUM=0.0D0
        DO 15 IGR=1,NGRP
        CALL LCMGDL(KPGPT,IGR,SOUR)
        CALL LCMGDL(KPFLX,IGR,DFLUX)
        DO 10 IUN=1,NUN
        SUM=SUM+SOUR(IUN)*DFLUX(IUN)
   10   CONTINUE
   15   CONTINUE
        DERIV(IVAR,ICST)=SUM
   20   CONTINUE
   25   CONTINUE
      ELSE IF(ITYPE.EQ.1000) THEN
*       IMPLICIT APPROACH
        IF(NVAR.NE.ND) CALL XABORT('GRA001: INVALID ND(2)')
        IF(NCST+1.NE.NGPT) CALL XABORT('GRA001: INVALID NGPT(2)')
        JPFLX=LCMGID(IPFLX,'ADFLUX')
        JPGPT=LCMGID(IPGPT,'DSOUR')
        DO 45 ICST=1,NCST+1
        KPFLX=LCMGIL(JPFLX,ICST)
        DO 40 IVAR=1,NVAR
        KPGPT=LCMGIL(JPGPT,IVAR)
        SUM=0.0D0
        DO 35 IGR=1,NGRP
        CALL LCMGDL(KPGPT,IGR,SOUR)
        CALL LCMGDL(KPFLX,IGR,DFLUX)
        DO 30 IUN=1,NUN
        SUM=SUM+SOUR(IUN)*DFLUX(IUN)
   30   CONTINUE
   35   CONTINUE
        DERIV(IVAR,ICST)=SUM
   40   CONTINUE
   45   CONTINUE
      ELSE
        CALL XABORT('GRA001: INVALID FLUX OBJECT')
      ENDIF
      DEALLOCATE(DFLUX,SOUR)
      RETURN
      END
