*DECK EDIRAT
      SUBROUTINE EDIRAT(IOPERA,NREGIO,NBMIX,MATCOD,FLXINT,AFLUX,RATES,
     >                  SIGMAX,IMERGE,NMERGE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Evaluate reaction rates from cross sections.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau
*
*Parameters: input
* IOPERA  type of action taken:
*         = 2  add cross section (no flux);
*         = 1  add reaction rates;
*         = 0  evaluate integrated flux;
*         =-1  subtract reaction rates.
* NREGIO  number of regions.
* NBMIX   number of mixtures.
* MATCOD  material per region.
* FLXINT  integrated fluxes.
* AFLUX   adjoint fluxes.
* SIGMAX  cross section array.
* IMERGE  region merging matrix.
* NMERGE  number of merged regions.
*
*Parameters: input/output
* RATES   initial and final reaction rates.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER     IOPERA,NREGIO,NBMIX,MATCOD(NREGIO),IMERGE(NREGIO),
     >            NMERGE
      REAL        FLXINT(NREGIO),AFLUX(NREGIO),RATES(NMERGE),
     >            SIGMAX(0:NBMIX)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) ::  DBRAT
*----
*  SCRATCH STORAGE ALLOCATION
*   DBRAT   double-precision reaction rates.
*----
      ALLOCATE(DBRAT(0:NMERGE))
*----
*  INITIALIZE DOUBLE PRECISION REACTION RATES
*----
      DO 90 IREG=0,NMERGE
        DBRAT(IREG)=0.0D0
 90   CONTINUE
      IF(IOPERA.EQ.0) THEN
*----
*  INTEGRATED FLUXES
*----
        DO 100 IREG=1,NREGIO
          IRATME=IMERGE(IREG)
          DBRAT(IRATME)=DBRAT(IRATME)+DBLE(FLXINT(IREG))
 100    CONTINUE
      ELSE IF(IOPERA.EQ.1) THEN
*----
*  SUM REACTION RATES
*----
        DO 110 IREG=1,NREGIO
          IRATME=IMERGE(IREG)
          MATNUM=MATCOD(IREG)
          DBRAT(IRATME)=DBRAT(IRATME)+DBLE(FLXINT(IREG))
     >                 *DBLE(SIGMAX(MATNUM))
 110    CONTINUE
      ELSE IF(IOPERA.EQ.-1) THEN
*----
*  SUBSTRACT REACTION RATES
*----
        DO 120 IREG=1,NREGIO
          IRATME=IMERGE(IREG)
          MATNUM=MATCOD(IREG)
          DBRAT(IRATME)=DBRAT(IRATME)-DBLE(FLXINT(IREG))
     >                 *DBLE(SIGMAX(MATNUM))
 120    CONTINUE
      ELSE IF(IOPERA.EQ.2) THEN
*----
*  ADD CROSS SECTIONS
*----
        DO 130 IREG=1,NREGIO
          IRATME=IMERGE(IREG)
          MATNUM=MATCOD(IREG)
          DBRAT(IRATME)=DBRAT(IRATME)+DBLE(SIGMAX(MATNUM))
 130    CONTINUE
      ELSE IF(IOPERA.EQ.10) THEN
*----
*  INTEGRATED ADJOINT FLUXES
*----
        DO 140 IREG=1,NREGIO
          IRATME=IMERGE(IREG)
          DBRAT(IRATME)=DBRAT(IRATME)+DBLE(FLXINT(IREG))
     >                 *DBLE(AFLUX(IREG))
 140    CONTINUE
      ELSE IF(IOPERA.EQ.11) THEN
*----
*  SUM ADJOINT-WEIGHTED REACTION RATES
*----
        DO 150 IREG=1,NREGIO
          IRATME=IMERGE(IREG)
          MATNUM=MATCOD(IREG)
          DBRAT(IRATME)=DBRAT(IRATME)+DBLE(FLXINT(IREG))
     >                 *DBLE(AFLUX(IREG))*DBLE(SIGMAX(MATNUM))
 150    CONTINUE
      ELSE IF(IOPERA.EQ.-11) THEN
*----
*  SUBSTRACT ADJOINT-WEIGHTED REACTION RATES
*----
        DO 160 IREG=1,NREGIO
          IRATME=IMERGE(IREG)
          MATNUM=MATCOD(IREG)
          DBRAT(IRATME)=DBRAT(IRATME)-DBLE(FLXINT(IREG))
     >                 *DBLE(AFLUX(IREG))*DBLE(SIGMAX(MATNUM))
 160    CONTINUE
      ENDIF
*----
*  STORE DOUBLE PRECISION REACTION RATES IN SINGLE PRECISION VECTOR
*----
      DO 170 IREG=1,NMERGE
        RATES(IREG)=RATES(IREG)+REAL(DBRAT(IREG))
 170  CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(DBRAT)
      RETURN
      END
