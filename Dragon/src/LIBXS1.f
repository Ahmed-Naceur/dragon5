*DECK LIBXS1
      SUBROUTINE LIBXS1(CFILNA,NEL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Initialize dimensions for depletion data with APOLIB-XSM.
*
*Copyright:
* Copyright (C) 2014 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): A. Hebert
*
*Parameters: input
* CFILNA  APOLIB-XSM file name.
*
*Parameters: output
* NEL     number of isotopes on library.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER CFILNA*(*)
      INTEGER NEL
*----
*  Local variables
*----
      TYPE(C_PTR) IPAP
*
      CALL LCMOP(IPAP,CFILNA,2,2,0)
      CALL LCMSIX(IPAP,'PHEAD',1)
      CALL LCMLEN(IPAP,'NOM',NV,ITYLCM)
      NEL=NV/5
      CALL LCMSIX(IPAP,' ',2)
      CALL LCMCL(IPAP,1)
      RETURN
      END
