*DECK ACRNDF
      SUBROUTINE ACRNDF(IMPX,NBISO1,ISO,IBM,HNOMIS,IPAPX,IPLIB,MY1,MY2,
     1 YLDS,ISTYP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Store records PYNAM, PYMIX and PYIELD into a Microlib.
*
*Copyright:
* Copyright (C) 2021 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* IMPX    print parameter (equal to zero for no print).
* NBISO1  number of particularized isotopes.
* ISO     particularized isotope index.
* IBM     material mixture.
* HNOMIS  array containing the names of the particularized isotopes.
* IPAPX   address of the Apex file.
* IPLIB   address of the output microlib LCM object.
* MY1     number of fissile isotopes including macroscopic sets.
* MY2     number of fission fragment.
* YLDS    fission yields.
*
*Parameters: output
* ISTYP   type of isotope ISO (=1: stable;=2: fissile; =3: fission
*         product).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE hdf5_wrap
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPAPX,IPLIB
      INTEGER IMPX,NBISO1,ISO,IBM,MY1,MY2,ISTYP
      DOUBLE PRECISION YLDS(MY1,MY2)
      CHARACTER(LEN=8) HNOMIS(NBISO1)
*----
*  LOCAL VARIABLES
*----
      INTEGER I, IOF, NBISO
*----
*  ALLOCATABLE AYYAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IPYMIX
      REAL, ALLOCATABLE, DIMENSION(:) :: PYIELD
      CHARACTER(LEN=4), ALLOCATABLE, DIMENSION(:) :: TYPISO
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: HPYNAM
*
      IF(hdf5_group_exists(IPAPX,"/physconst/")) THEN
        CALL hdf5_read_data(IPAPX,"/physconst/ISOTYP",TYPISO)
      ELSE IF(hdf5_group_exists(IPAPX,"/physco001/")) THEN
        CALL hdf5_read_data(IPAPX,"/physco001/ISOTYP",TYPISO)
      ELSE
        CALL XABORT('ACRNDF: GROUP physconst NOT FOUND IN HDF5 FILE.')
      ENDIF
      NBISO=SIZE(TYPISO,1)
      IF(ISO.LE.NBISO) THEN
        IF(TYPISO(ISO).EQ.'OTHE') ISTYP=1
        IF(TYPISO(ISO).EQ.'FISS') ISTYP=2
        IF(TYPISO(ISO).EQ.'F.P.') ISTYP=3
      ELSE
        ISTYP=1
      ENDIF
      IF(ISTYP.EQ.3) THEN
        ALLOCATE(HPYNAM(MY1),PYIELD(MY1),IPYMIX(MY1))
        IOF=0
        DO I=1,NBISO
          IF(TYPISO(I).EQ.'FISS') THEN
            IOF=IOF+1
            IF(IOF.GT.MY1) CALL XABORT('ACRNDF: MY1 OVERFLOW.')
            HPYNAM(IOF)=HNOMIS(I)
            IPYMIX(IOF)=IBM
            PYIELD(IOF)=REAL(YLDS(IOF,ISO))
          ENDIF
        ENDDO
        DO I=NBISO+1,NBISO1
            IOF=IOF+1
            IF(IOF.GT.MY1) CALL XABORT('ACRNDF: MY1 OVERFLOW.')
            HPYNAM(IOF)=HNOMIS(I)
            IPYMIX(IOF)=IBM
            PYIELD(IOF)=0.0
        ENDDO
        IF(IOF.NE.MY1) CALL XABORT('ACRNDF: MY1 COUNT ERROR.')
        CALL LCMPTC(IPLIB,'PYNAM',8,MY1,HPYNAM)
        CALL LCMPUT(IPLIB,'PYMIX',MY1,1,IPYMIX)
        CALL LCMPUT(IPLIB,'PYIELD',MY1,2,PYIELD)
        IF(IMPX.GT.2) THEN
          WRITE(6,'(3X,7HPYIELD=,1P,8E12.4/(8X,10E12.4))') (PYIELD(I),
     1    I=1,MY1)
        ENDIF
        DEALLOCATE(IPYMIX,PYIELD,HPYNAM)
      ENDIF
      DEALLOCATE(TYPISO)
      RETURN
      END
