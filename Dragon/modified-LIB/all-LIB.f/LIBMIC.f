*DECK LIBMIC
      SUBROUTINE LIBMIC (IPLIB,IPMIC,NAMFIL,NGRO,NBISO,ISONAM,ISONRF,
     1 IPISO,MASKI,IMPX,NGF,NGFR,NDEL,NBESP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Transcription of the useful microscopic cross section data from a
* microlib to LCM data structures.
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
* IPLIB   pointer to the lattice microscopic cross section library
*         (L_LIBRARY signature).
* IPMIC   pointer to the draglib (L_DRAGLIB signature).
* NAMFIL  name of the Dragon library file.
* NGRO    number of energy groups.
* NBISO   number of isotopes present in the calculation domain.
* ISONAM  alias name of isotopes.
* ISONRF  library name of isotopes.
* IPISO   pointer array towards microlib isotopes.
* MASKI   isotopic mask. Isotope with index I is processed if
*         MASKI(I)=.true.
* IMPX    print flag.
*
*Parameters: output
* NGF     number of fast groups without self-shielding.
* NGFR    number of fast and resonance groups.
* NDEL    number of precursor groups for delayed neutrons.
* NBESP   number of energy-dependent fission spectra.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      PARAMETER(MAXESP=4)
      CHARACTER*(*) NAMFIL
      TYPE(C_PTR) IPLIB,IPMIC,IPISO(NBISO)
      INTEGER NGRO,NBISO,ISONAM(3,NBISO),ISONRF(3,NBISO),IMPX,
     1 NGF,NGFR,NDEL,NBESP
      LOGICAL MASKI(NBISO)
*----
*  LOCAL VARIABLES
*----
      CHARACTER HSMG*131,HTITLE*80,HNISOR*12,HNAMIS*12
      PARAMETER (IOUT=6,NOTX=3,NSTATE=40)
      TYPE(C_PTR) KPLIB,JPMIC,KPMIC
      INTEGER IESP(MAXESP+1),ISTATE(NSTATE)
      REAL EESP(MAXESP+1)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ITITLE
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: JSONAM
      REAL, ALLOCATABLE, DIMENSION(:) :: DELTA,ENER
*----
*  RECOVER THE GROUP STRUCTURE.
*----
      NGF=NGRO+1
      NGFR=0
      NDEL=0
      IF(IMPX.GT.0) WRITE (IOUT,900) NAMFIL
      CALL LCMLEN(IPMIC,'README',LENGT,ITYLCM)
      IF((IMPX.GT.0).AND.(LENGT.GT.0)) THEN
         ALLOCATE(ITITLE(LENGT))
         CALL LCMGET(IPMIC,'README',ITITLE)
         WRITE (IOUT,940)
         I2=0
         DO 10 J=0,LENGT/20
         I1=I2+1
         I2=MIN(I1+19,LENGT)
         WRITE (HTITLE,'(20A4)') (ITITLE(I),I=I1,I2)
         WRITE (IOUT,'(1X,A80)') HTITLE
   10    CONTINUE
         DEALLOCATE(ITITLE)
      ENDIF
      ALLOCATE(DELTA(NGRO),ENER(NGRO+1))
      CALL LCMLEN(IPMIC,'ENERGY',LENGT,ITYLCM)
      LENGT=LENGT-1
      IF(LENGT.NE.NGRO) CALL XABORT('LIBMIC: INVALID GROUP STRUCTURE.')
      CALL LCMGET(IPMIC,'ENERGY',ENER)
      CALL LCMLEN(IPMIC,'DELTAU',LENGT,ITYLCM)
      IF(LENGT.EQ.NGRO) THEN
         CALL LCMGET(IPMIC,'DELTAU',DELTA)
      ELSE IF(LENGT.EQ.0) THEN
         IF(ENER(NGRO+1).EQ.0.0) ENER(NGRO+1)=1.0E-5
         DO 20 J=1,NGRO
         DELTA(J)=LOG(ENER(J)/ENER(J+1))
   20    CONTINUE
      ENDIF
      CALL LCMPUT(IPLIB,'ENERGY',NGRO+1,2,ENER)
      CALL LCMPUT(IPLIB,'DELTAU',NGRO,2,DELTA)
      DEALLOCATE(ENER,DELTA)
      CALL LCMLEN(IPMIC,'CHI-LIMITS',NBESP,ITYLCM)
      IF(NBESP.GT.0) THEN
         NBESP=NBESP-1
         IF(NBESP.GT.MAXESP) CALL XABORT('LIBMIC: MAXESP OVERFLOW.')
         CALL LCMGET(IPMIC,'CHI-LIMITS',IESP)
         CALL LCMPUT(IPLIB,'CHI-LIMITS',NBESP+1,1,IESP)
         CALL LCMGET(IPMIC,'CHI-ENERGY',EESP)
         CALL LCMPUT(IPLIB,'CHI-ENERGY',NBESP+1,2,EESP)
      ENDIF
*----
*  READ THROUGH MICROLIB AND ACCUMULATE CROSS SECTIONS.
*----
      CALL LCMGET(IPMIC,'STATE-VECTOR',ISTATE)
      NBML=ISTATE(2)
      ALLOCATE(JSONAM(3,NBML))
      CALL LCMGET(IPMIC,'ISOTOPESUSED',JSONAM)
      JPMIC=LCMGID(IPMIC,'ISOTOPESLIST')
      DO 40 IMX=1,NBISO
      IF(MASKI(IMX)) THEN
         WRITE(HNAMIS,'(3A4)') (ISONAM(I0,IMX),I0=1,3)
         WRITE(HNISOR,'(3A4)') (ISONRF(I0,IMX),I0=1,3)
         KML=0
         DO IML=1,NBML
           IF((ISONAM(1,IMX).EQ.JSONAM(1,IML)).AND.
     1        (ISONAM(2,IMX).EQ.JSONAM(2,IML)).AND.
     2        (ISONAM(3,IMX).EQ.JSONAM(3,IML))) THEN
             KML=IML
             GO TO 30
           ENDIF
         ENDDO
         WRITE (HSMG,910) HNAMIS,HNISOR,NAMFIL
         CALL XABORT(HSMG)
   30    KPMIC=LCMGIL(JPMIC,KML) ! set KML-th isotope
         KPLIB=IPISO(IMX) ! set IMX-th isotope
         IF(.NOT.C_ASSOCIATED(KPLIB)) THEN
           WRITE(HSMG,'(17HLIBMIC: ISOTOPE '',3A4,7H'' (ISO=,I8,
     1     35H) IS NOT AVAILABLE IN THE MICROLIB.)') (ISONAM(I0,IMX),
     2     I0=1,3),IMX
           CALL XABORT(HSMG)
         ENDIF
         CALL LCMEQU(KPMIC,KPLIB)
      ENDIF
   40 CONTINUE
      DEALLOCATE(JSONAM)
      RETURN
*
  900 FORMAT(/27H PROCESSING MICROLIB NAMED ,A12,1H.)
  910 FORMAT(26HLIBMIC: MATERIAL/ISOTOPE ',A12,5H' = ',A12,9H' IS MISS,
     1 22HING ON MICROLIB NAMED ,A12,1H.)
  940 FORMAT(/24H X-SECTION LIBRARY INFO:)
      END
