*DECK LIBIPS
      SUBROUTINE LIBIPS(IPLIB,NBISO,IPISO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Construct the pointer array towards microlib isotopes.
*
*Copyright:
* Copyright (C) 2017 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPLIB   pointer to the lattice microscopic cross section library
*         (L_LIBRARY signature).
* NBISO   number of isotopes present in the calculation domain.
*
*Parameters: output
* IPISO   pointer array towards microlib isotopes.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NBISO
      TYPE(C_PTR) IPLIB,IPISO(NBISO)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPLIB
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISONAM
*
      ALLOCATE(ISONAM(3,NBISO))
      CALL LCMGET(IPLIB,'ISOTOPESUSED',ISONAM)
      JPLIB=LCMGID(IPLIB,'ISOTOPESLIST')
      DO 20 ISOT=1,NBISO
      CALL LCMLEL(JPLIB,ISOT,LENGTH,ITYLCM)
      IF(LENGTH.NE.0) THEN
        IPISO(ISOT)=LCMGIL(JPLIB,ISOT) ! set ISOT-th isotope
      ELSE
        DO 10 JSOT=1,ISOT-1
        IF((ISONAM(1,ISOT).EQ.ISONAM(1,JSOT)).AND.(ISONAM(2,ISOT).EQ.
     >  ISONAM(2,JSOT)).AND.(ISONAM(3,ISOT).EQ.ISONAM(3,JSOT))) THEN
          IPISO(ISOT)=IPISO(JSOT) ! set JSOT-th isotope
          GO TO 20
        ENDIF
   10   CONTINUE
        IPISO(ISOT)=C_NULL_PTR
      ENDIF
   20 CONTINUE
      RETURN
      END
