*DECK LIBCOM
      SUBROUTINE LIBCOM(NFS,DELTA,SIGAF,SIGTF,NORA,NOR,COMOM)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute a set of comoments ((SIGA**P)*(SIGT**Q)).
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
* NFS    number of fine energy groups.
* DELTA  lethargy widths of the fine groups.
* SIGAF  microscopic absorption x-sections in the fine groups.
* SIGTF  microscopic total x-sections in the fine groups.
* NORA   related to the number of absorption moments to preserve.
* NOR    related to the number of total moments to preserve:
*        (2-NORA)/2 <= P <= (NORA+1)/2 and (2-NOR)/2 <= Q <= (NOR+1)/2.
*
*Parameters: output
* COMOM  comoments.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NFS,NORA,NOR
      REAL DELTA(NFS),SIGAF(NFS),SIGTF(NFS)
      DOUBLE PRECISION COMOM(NORA,NOR)
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION DEL,T,T0,SIGT,SIGA
      INTEGER PNOR,QNOR
*
      DEL=0.0D0
      DO 15 PNOR=1,NORA
      DO 10 QNOR=1,NOR
      COMOM(PNOR,QNOR)=0.0D0
   10 CONTINUE
   15 CONTINUE
*
      DO 80 IGF=1,NFS
        DELF=DELTA(IGF)
        SIGT=MAX(0.001,SIGTF(IGF))
        SIGA=MAX(0.001,SIGAF(IGF))
        DEL=DEL+DELF
        T0=DELF
        DO 40 PNOR=MAX(1,NORA/2),NORA
          T=T0
          DO 20 QNOR=MAX(1,NOR/2),NOR
            COMOM(PNOR,QNOR)=COMOM(PNOR,QNOR)+T
            T=T*SIGT
   20     CONTINUE
          T=T0/SIGT
          DO 30 QNOR=NOR/2-1,1,-1
            COMOM(PNOR,QNOR)=COMOM(PNOR,QNOR)+T
            T=T/SIGT
   30       CONTINUE
          T0=T0*SIGA
   40   CONTINUE
        T0=DELF/SIGA
        DO 70 PNOR=NORA/2-1,1,-1
          T=T0
          DO 50 QNOR=MAX(1,NOR/2),NOR
            COMOM(PNOR,QNOR)=COMOM(PNOR,QNOR)+T
            T=T*SIGT
   50     CONTINUE
          T=T0/SIGT
          DO 60 QNOR=NOR/2-1,1,-1
            COMOM(PNOR,QNOR)=COMOM(PNOR,QNOR)+T
            T=T/SIGT
   60     CONTINUE
          T0=T0/SIGA
   70   CONTINUE
   80 CONTINUE
*
      IF(DEL.EQ.0.0) CALL XABORT('LIBCOM: ALGORITHM FAILURE.')
      DO 100 PNOR=1,NORA
      DO 90 QNOR=1,NOR
        COMOM(PNOR,QNOR)=COMOM(PNOR,QNOR)/DEL
   90 CONTINUE
  100 CONTINUE
      RETURN
      END
