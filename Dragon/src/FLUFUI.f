*DECK FLUFUI
      FUNCTION FLUFUI(IPSYS,NREGIO,NUNKNO,MATCOD,VOLUME,KEYFLX,FUNKNO,
     >                SUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the DB2 value corresponding to the actual leakage.
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
* IPSYS   pointer to the pij LCM object.
* NREGIO  number of regions considered.
* NUNKNO  number of unknown in the system.
* MATCOD  material code in region.
* VOLUME  volume of region.
* KEYFLX  flux elements in unknown system.
* FUNKNO  unknown vector.
* SUNKNO  source for system of unknown.
*
*Parameters: output
* FLUFUI  leakage DB2 parameter.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPSYS
      INTEGER     NREGIO,NUNKNO,MATCOD(NREGIO),KEYFLX(NREGIO)
      REAL        VOLUME(NREGIO),FUNKNO(NUNKNO),SUNKNO(NUNKNO),FLUFUI
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: SIGT0,SIGS0
*----
*  READ THE TRANSPORT CORRECTED CROSS SECTIONS.
*----
      CALL LCMLEN(IPSYS,'DRAGON-TXSC',ILCTXS,ITYLCM)
      ALLOCATE(SIGT0(0:ILCTXS-1))
      CALL LCMGET(IPSYS,'DRAGON-TXSC',SIGT0(0))
      CALL LCMLEN(IPSYS,'DRAGON-S0XSC',ILCS0X,ITYLCM)
      ALLOCATE(SIGS0(0:ILCS0X-1))
      CALL LCMGET(IPSYS,'DRAGON-S0XSC',SIGS0(0))
*----
*  COMPUTE DB2
*----
      ZNUM=0.0
      ZDEN=0.0
      DO 10 I=1,NREGIO
      SSS=SIGT0(MATCOD(I))-SIGS0(MATCOD(I))
      ZNUM=ZNUM+VOLUME(I)*(SUNKNO(KEYFLX(I))-SSS*FUNKNO(KEYFLX(I)))
      ZDEN=ZDEN+VOLUME(I)*FUNKNO(KEYFLX(I))
   10 CONTINUE
      FLUFUI=0.0
      IF(ZDEN.GT.0.0) FLUFUI=ZNUM/ZDEN
*
      DEALLOCATE(SIGS0,SIGT0)
      RETURN
      END
