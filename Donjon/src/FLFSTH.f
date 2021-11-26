*DECK FLFSTH
      SUBROUTINE FLFSTH(PTOT,POWER,POWC,POWB,FLUX,NGRP,NCH,
     +            NB,NEL,FSTH,FLUB)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Update the fuel average fluxes and the channel and bundle powers
* over the fuel lattice using FTSH
*
*Copyright:
* Copyright (C) 2010 Ecole Polytechnique de Montreal.
*
*Author(s): 
* M. Guyot
*
*Parameters: input
* PTOT   total power in MW
* POWER  total power computed with H-factors in MW
* POWC   channel powers in kW
* POWB   bundle powers in kW
* FLUX   average fluxes per regions
* NGRP   number of energy groups
* NCH    number of reactor channels.
* NB     number of fuel bundles per channel.
* NEL    total number of finite elements.
* FSTH   thermal to fission ratio power
* FLUB   average fluxers per bundles
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NCH,NB,NGRP,NEL
      REAL FLUX(NEL,NGRP),FSTH,FLUB(NCH,NB,NGRP),
     1     POWB(NCH,NB),POWC(NCH)
      DOUBLE PRECISION POWER,FACT,PTOT
*----
*  LOCAL VARIABLES
*----
      INTEGER I,J,K
*
      FACT=PTOT/POWER
      FACT= FACT/FSTH
      DO 10 I=1,NCH
        POWC(I)=POWC(I)*REAL(FACT)
        DO 20 J=1,NB
          POWB(I,J)=POWB(I,J)*REAL(FACT)
          DO 30 K=1,NGRP
            FLUB(I,J,K)=FLUB(I,J,K)*REAL(FACT)  
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE
      DO 40 I=1,NEL
        DO 50 J=1,NGRP
          FLUX(I,J)=FLUX(I,J)*REAL(FACT)
   50   CONTINUE
   40 CONTINUE
      RETURN
      END
