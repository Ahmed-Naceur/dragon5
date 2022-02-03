*DECK OUTFLX
      SUBROUTINE OUTFLX(IPFLUX,ITYP,NGRP,NUN,LMOD,IMPX,EVECT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover the direct or adjoint flux.
*
*Copyright:
* Copyright (C) 2012 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A Hebert
*
*Parameters: input
* IPFLUX  L_FLUX pointer to the solution.
* ITYP    type of flux (=0: direct; =1: adjoint).
* NGRP    total number of energy groups.
* NUN     total number of unknowns per group.
* LMOD    index of mode.
* IMPX    print flag.
*
*Parameters: output
* EVECT   flux.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPFLUX
      INTEGER ITYP,NGRP,NUN,LMOD,IMPX
      REAL EVECT(NUN,NGRP)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPFLUX,KPFLUX,MPFLUX
*
      IF(ITYP.EQ.0) THEN
*        RECOVER THE DIRECT FLUX.
         IF(IMPX.GT.0) WRITE(6,20) 'DIRECT'
         CALL LCMLEN(IPFLUX,'K-EFFECTIVE',ILEN,ITYLCM)
         IF(ILEN.GT.0) CALL LCMGET(IPFLUX,'K-EFFECTIVE',FKEFF)
         CALL LCMLEN(IPFLUX,'FLUX',LENGT,ITYLCM)
         IF(LENGT.GT.0) THEN
            MPFLUX=LCMGID(IPFLUX,'FLUX')
         ELSE
            CALL LCMLEN(IPFLUX,'MODE',LENGT,ITYLCM)
            IF(LENGT.GT.0) THEN
               IF(LMOD.LE.0) CALL XABORT('OUTFLX: INVALID MODE INDEX.')
               JPFLUX=LCMGID(IPFLUX,'MODE')
               KPFLUX=LCMGIL(JPFLUX,LMOD)
               MPFLUX=LCMGID(KPFLUX,'FLUX')
            ELSE
               CALL LCMLIB(IPFLUX)
               CALL XABORT('OUTFLX: UNABLE TO RECOVER A DIRECT FLUX.')
            ENDIF
         ENDIF
      ELSE IF(ITYP.EQ.1) THEN
*        RECOVER THE ADJOINT FLUX.
         IF(IMPX.GT.0) WRITE(6,20) 'ADJOINT'
         CALL LCMLEN(IPFLUX,'AK-EFFECTIVE',ILEN,ITYLCM)
         IF(ILEN.GT.0) CALL LCMGET(IPFLUX,'AK-EFFECTIVE',FKEFF)
         CALL LCMLEN(IPFLUX,'AFLUX',LENGT,ITYLCM)
         IF(LENGT.GT.0) THEN
            MPFLUX=LCMGID(IPFLUX,'AFLUX')
         ELSE
            CALL LCMLEN(IPFLUX,'MODE',LENGT,ITYLCM)
            IF(LENGT.GT.0) THEN
               IF(LMOD.LE.0) CALL XABORT('OUTFLX: INVALID MODE INDEX.')
               JPFLUX=LCMGID(IPFLUX,'MODE')
               KPFLUX=LCMGIL(JPFLUX,LMOD)
               MPFLUX=LCMGID(KPFLUX,'AFLUX')
            ELSE
               CALL LCMLIB(IPFLUX)
               CALL XABORT('OUTFLX: UNABLE TO RECOVER AN ADJOINT FLUX.')
            ENDIF
         ENDIF
      ENDIF
      DO 10 IGR=1,NGRP
      CALL LCMGDL(MPFLUX,IGR,EVECT(1,IGR))
   10 CONTINUE
      RETURN
   20 FORMAT(/21H OUTFLX: RECOVER THE ,A,6H FLUX.)
      END
