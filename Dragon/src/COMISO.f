*DECK COMISO
      SUBROUTINE COMISO(ITYP,MAXISO,IPLIB,NISO,NOMISO,NOMEVO,TYPISO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover the names of the isotopes stored in a microlib.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* ITYP    type of operation:
*         =0:  check the values of the isotope names and types;
*         =-1: recover all isotopes;
*         =-2: recover fissiles isotopes;
*         =-3: recover fission products;
*         >0:  recover all isotopes in mixture ITYP.
* MAXISO  dimension of arrays NOMISO and TYPISO.
* IPLIB   pointer to the microlib (L_LIBRARY signature).
*
*Parameters: input/output
* NISO    number of particularized isotopes.
* NOMISO  alias names of the particularized isotopes.
*
*Parameters: output
* NOMEVO  library names of the particularized isotopes.
* TYPISO  type of each isotope:
*         =1: the isotope is not fissile and not a fission product;
*         =2: the isotope is fissile;
*         =3: the isotope is a fission product.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB
      INTEGER ITYP,MAXISO,NISO,TYPISO(MAXISO)
      CHARACTER NOMISO(MAXISO)*8,NOMEVO(MAXISO)*12
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      CHARACTER HNAME*8
      INTEGER ISTATE(NSTATE)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ISUSED,ISNEVO,ISMIX,ISTYP
*
      IF(.NOT.C_ASSOCIATED(IPLIB)) RETURN
      CALL LCMGET(IPLIB,'STATE-VECTOR',ISTATE)
      NBISOT=ISTATE(2)
      ALLOCATE(ISUSED(3*NBISOT),ISNEVO(3*NBISOT),ISMIX(NBISOT),
     1 ISTYP(NBISOT))
      CALL LCMGET(IPLIB,'ISOTOPESUSED',ISUSED)
      CALL LCMGET(IPLIB,'ISOTOPERNAME',ISNEVO)
      CALL LCMGET(IPLIB,'ISOTOPESMIX',ISMIX)
      CALL LCMGET(IPLIB,'ISOTOPESTYPE',ISTYP)
      IF(ITYP.EQ.0) THEN
         DO 15 ISOT=1,NBISOT
         WRITE(HNAME,'(2A4)') (ISUSED((ISOT-1)*3+I0),I0=1,2)
         DO 10 I=1,NISO
         IF(NOMISO(I).EQ.HNAME) THEN
            TYPISO(I)=MAX(TYPISO(I),ISTYP(ISOT))
            WRITE(NOMEVO(I),'(3A4)') (ISNEVO((ISOT-1)*3+I0),I0=1,3)
         ENDIF
   10    CONTINUE
   15    CONTINUE
         DO 20 I=1,NISO
         IF(TYPISO(I).EQ.0) THEN
            HNAME=NOMISO(I)
            CALL XABORT('COMISO: UNABLE TO FIND ISOTOPE '//HNAME//
     1      ' IN THE MICROLIB.')
         ENDIF
   20    CONTINUE
      ELSE
         DO 40 ISOT=1,NBISOT
         WRITE(HNAME,'(2A4)') (ISUSED((ISOT-1)*3+I0),I0=1,2)
         DO 30 I=1,NISO
         IF(NOMISO(I).EQ.HNAME) GO TO 40
   30    CONTINUE
         IMIX=ISMIX(ISOT)
         JTYP=ISTYP(ISOT)
         IF((ITYP.EQ.-1).OR.(ITYP.EQ.-JTYP).OR.(ITYP.EQ.IMIX)) THEN
            NISO=NISO+1
            NOMISO(NISO)=HNAME
            WRITE(NOMEVO(NISO),'(3A4)') (ISNEVO((ISOT-1)*3+I0),I0=1,3)
            TYPISO(NISO)=0
         ENDIF
   40    CONTINUE
      ENDIF
      DEALLOCATE(ISTYP,ISMIX,ISNEVO,ISUSED)
      RETURN
      END
