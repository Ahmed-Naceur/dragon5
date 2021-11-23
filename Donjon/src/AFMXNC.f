*DECK AFMXNC
          SUBROUTINE AFMXNC (NGRP,SIGX,SIGF,FLUX,XXE,XNP,FLUR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Computation of Xenon and Neptunium concentrations.
*
*Copyright:
* Copyright (C) 1996 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): 
* M.T. Sissaoui
*
*Parameters: input
* NGRP
* SIGX  Xenon absorption micro-x-section dimension (ngrp).
* SIGF  fission macro-x-section dimension (ngrp).
* FLUX  flux dimension (ngrp)
*
*Parameters: output
* XXE   Xenon concentration
* XNP   Neptunium concentration
* FLUR
*
*-----------------------------------------------------------------------
*
      DIMENSION FLUX(NGRP),SIGF(NGRP),SIGX(NGRP),FLUR(NGRP)
      REAL      CF
* SET THE YIELD AND THE DECAY CONSTANTE FOR XENON AND NEPTUNIUM
      XLAMBDAX = 2.09E-5
      XLAMBDAI = 2.85E-5
      GAMMAI = 0.0631
      GAMMAX = 0.0045
* CF=1.E-24(barn)
      CF=1.0E-24
      CINTG=1.0E+13
*  CALCUL DES TAUX DE FISSION
      TAUF=0.0
      TAUAX=0.0
      FLR=0.0
      FLX=0.0
      DO 10 IGR = 1,NGRP
        TAUF = TAUF+FLUX(IGR)*SIGF(IGR)
        TAUAX = TAUAX+FLUX(IGR)*SIGX(IGR)
        FLR=FLR+FLUR(IGR)*CINTG
        FLX=FLX+FLUX(IGR)
 10   CONTINUE
*  COMPUTE THE XENON CONCENTRATION
      XXE=CF*(GAMMAX+GAMMAI)*TAUF/(XLAMBDAX+TAUAX*CF)
*  COMPUTE THE NEPTUNIUM CONCENTRATION
      XNP=XNP*FLX/FLR
*
      RETURN
      END
