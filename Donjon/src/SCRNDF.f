*DECK SCRNDF
      SUBROUTINE SCRNDF(IMPX,NBISO1,ISO,IBM,INOMIS,IPMEM,IPLIB,NCAL,
     1 TERP,MY1,MY2,YLDS,ISTYP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Store records PYNAM, PYMIX and PYIELD into a Microlib.
*
*Copyright:
* Copyright (C) 2015 Ecole Polytechnique de Montreal
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
* INOMIS  array containing the names of the particularized isotopes.
* IPMEM   pointer to the memory-resident Saphyb object.
* IPLIB   address of the output microlib LCM object.
* NCAL    number of elementary calculations in the Saphyb.
* TERP    interpolation factors.
* MY1     number of fissile isotopes including macroscopic sets.
* MY2     number of fission fragment.
* YLDS    fission yields.
*
*Parameters: output
* ISTYP   type of isotope ISO (=2: fissile; =3: fission product).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMEM,IPLIB
      INTEGER IMPX,NBISO1,ISO,IBM,INOMIS(2,NBISO1),NCAL,MY1,MY2,ISTYP
      REAL TERP(NCAL)
      DOUBLE PRECISION YLDS(MY1,MY2)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPMEM,KPMEM
      INTEGER I, I0, ICAL, IY1, IY2, JSO, NISY
*----
*  ALLOCATABLE AYYAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ADRY,IPYMIX
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IPYNAM
      REAL, ALLOCATABLE, DIMENSION(:) :: PYIELD
*
      JPMEM=LCMGID(IPMEM,'calc')
      ISTYP=0
      DO 10 ICAL=NCAL,1,-1
      IF(TERP(ICAL).EQ.0.0) GO TO 10
      KPMEM=LCMGIL(JPMEM,ICAL)
      CALL LCMSIX(KPMEM,'info',1)
      CALL LCMGET(KPMEM,'NISY',NISY)
      IF(ISO.GT.NISY) CALL XABORT('SCRNDF: NISY OVERFLOW.')
      ALLOCATE(ADRY(NISY))
      CALL LCMGET(KPMEM,'ADRY',ADRY)
      CALL LCMSIX(KPMEM,' ',2)
      IF(ADRY(ISO).GT.0) THEN
*       ISO is a fissile isotope
        ISTYP=2
      ELSE IF(ADRY(ISO).LT.0) THEN
*       ISO is a fission product
        ISTYP=3
        IY2=-ADRY(ISO)
        IF(IY2.GT.MY2) CALL XABORT('SCRNDF: MY2 OVERFLOW.')
        ALLOCATE(IPYNAM(2,MY1),IPYMIX(MY1),PYIELD(MY1))
        IPYNAM(:2,:MY1)=0
        IPYMIX(:MY1)=0
        PYIELD(:MY1)=0.0
        IF(IMPX.GT.2) THEN
          WRITE(6,'(25H SCRNDF: fission product=,2A4,9H mixture=,I8)')
     1    (INOMIS(I0,ISO),I0=1,2),IBM
        ENDIF
        DO JSO=1,NISY
          IF(ADRY(JSO).GT.0) THEN
            IY1=ADRY(JSO)
            IF(IY1.GT.MY1) CALL XABORT('SCRNDF: MY1 OVERFLOW.')
            IPYNAM(1,IY1)=INOMIS(1,JSO)
            IPYNAM(2,IY1)=INOMIS(2,JSO)
            IPYMIX(IY1)=IBM
            PYIELD(IY1)=REAL(YLDS(IY1,IY2))
            IF(IMPX.GT.2) THEN
              WRITE(6,'(9X,16Hfissile isotope(,I4,2H)=,2A4,9H mixture=,
     1        I8)') IY1,(IPYNAM(I0,IY1),I0=1,2),IPYMIX(IY1)
            ENDIF
          ENDIF
        ENDDO
        CALL LCMPUT(IPLIB,'PYNAM',2*MY1,3,IPYNAM)
        CALL LCMPUT(IPLIB,'PYMIX',MY1,1,IPYMIX)
        CALL LCMPUT(IPLIB,'PYIELD',MY1,2,PYIELD)
        IF(IMPX.GT.2) THEN
          WRITE(6,'(3X,7HPYIELD=,1P,8E12.4/(8X,10E12.4))') (PYIELD(I),
     1    I=1,MY1)
        ENDIF
        DEALLOCATE(PYIELD,IPYMIX,IPYNAM)
      ENDIF
      DEALLOCATE(ADRY)
      RETURN
   10 CONTINUE
      CALL XABORT('SCRNDF: UNABLE TO FIND A CALCULATION DIRECTORY.')
      RETURN
      END
