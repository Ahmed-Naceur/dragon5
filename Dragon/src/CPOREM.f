*DECK CPOREM
      SUBROUTINE CPOREM(NGROUP,NL    ,NPROC ,INDPRO,DENCPO,
     >                  DXSMIC,DSCMIC,DXSREM,DSCREM)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Remove Compo isotope xs from macroscopic xs.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): G. Marleau
*
*Parameters: input
* NGROUP  number of groups condensed.
* NL      number of Legendre orders.
* NPROC   number of microscopic xs to process.
* INDPRO  identifier for xs processing.
* DENCPO  Compo isotopes concentration.
* DXSMIC  microscopic vector xs.
* DSCMIC  microscopic scat matrix xs.
*
*Parameters: input/output
* DXSREM  averaged region/group x-s.
* DSCREM  scattering rates.
*
*-----------------------------------------------------------------------
*
      IMPLICIT         NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER          NGROUP,NL,NPROC,INDPRO(NPROC)
      DOUBLE PRECISION DENCPO,
     >                 DXSMIC(NGROUP,NPROC),
     >                 DSCMIC(NGROUP,NGROUP,NL),
     >                 DXSREM(NGROUP,NPROC),
     >                 DSCREM(NGROUP,NGROUP,NL)
*----
*  LOCAL PARAMETERS
*----
      INTEGER          NDPROC
      PARAMETER       (NDPROC=20)
      INTEGER          IGR,JGR,IXSR,IL
*----
*  REMOVE STANDARD XS
*----
      DO 100 IXSR=1,NDPROC
        IF(IXSR.NE.16.AND.INDPRO(IXSR).GT.0) THEN
          DO 110 IGR=1,NGROUP
            DXSREM(IGR,IXSR)=DXSREM(IGR,IXSR)
     >                      +DENCPO*DXSMIC(IGR,IXSR)
 110      CONTINUE
        ENDIF
 100  CONTINUE
*----
*  REMOVE SCATTERING XS
*----
      IL=0
      DO 120 IXSR=NDPROC+1,NDPROC+NL
        IL=IL+1
        IF(INDPRO(IXSR).GT.0) THEN
          DO 130 IGR=1,NGROUP
            DXSREM(IGR,IXSR)=DXSREM(IGR,IXSR)
     >        +DENCPO*DXSMIC(IGR,IXSR)
            DO 131 JGR=1,NGROUP
              DSCREM(IGR,JGR,IL)=DSCREM(IGR,JGR,IL)
     >          +DENCPO*DSCMIC(IGR,JGR,IL)
 131        CONTINUE
 130      CONTINUE
        ENDIF
 120  CONTINUE
      RETURN
      END
