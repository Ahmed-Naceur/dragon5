*DECK LIBADJ
      SUBROUTINE LIBADJ (IPLIB,NGRO,NBISO,NL,NDEL,NBESP,IPISO,NED,
     1 NAMEAD)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Transposition of the usefull interpolated microscopic cross section
* for producing an adjoint problem.
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
* NGRO    number of energy groups.
* NBISO   number of isotopes present in the calculation domain.
* NL      number of Legendre orders required in the calculation
*         NL=1 or higher.
* NDEL    number of delayed precursor groups.
* NBESP   number of energy-dependent fission spectra.
* IPISO   pointer array towards microlib isotopes.
* NED     number of extra vector edits from matxs.
* NAMEAD  matxs names of the extra vector edits.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB,IPISO(NBISO)
      INTEGER NGRO,NBISO,NL,NDEL,NBESP,NED,NAMEAD(2,NED)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPLIB
      INTEGER I,J,I0,IED,IDEL,IL,IMPX,IMT,INGRO,LENGT,ITYLCM
      REAL SUM
      CHARACTER TEXT8*8,HNUSIG*12,HCHI*12
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ITYPRO
      REAL, ALLOCATABLE, DIMENSION(:,:) :: GA1,GA2,SIGS
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SCAT
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(ITYPRO(NL),GA1(NGRO,2),GA2(NGRO,NGRO),SIGS(NGRO,NL),
     1 SCAT(NGRO,NGRO,NL))
*----
*  ***MATERIAL/ISOTOPE LOOP***
*----
      IF(NBESP.NE.0) CALL XABORT('LIBADJ: MULTIPLE FISSION SPECTRA NOT'
     1 //' IMPLEMENTED.')
      IMPX=0
      DO 200 IMT=1,NBISO
      JPLIB=IPISO(IMT)
      IF(.NOT.C_ASSOCIATED(JPLIB)) GO TO 200
      CALL XDRLGS(JPLIB,-1,IMPX,0,NL-1,1,NGRO,SIGS,SCAT,ITYPRO)
      INGRO=NL-1
      DO 10 IL=NL-1,0,-1
      IF(ITYPRO(IL+1).EQ.0) THEN
         INGRO=INGRO-1
      ELSE
         GO TO 20
      ENDIF
   10 CONTINUE
   20 DO 50 IL=0,INGRO
      IF(ITYPRO(IL+1).GT.0) THEN
         DO 35 I=1,NGRO
         GA1(I,1)=SIGS(NGRO-I+1,IL+1)
         DO 30 J=1,NGRO
         GA2(I,J)=SCAT(NGRO-J+1,NGRO-I+1,IL+1)
   30    CONTINUE
   35    CONTINUE
         DO 45 I=1,NGRO
         SIGS(I,IL+1)=GA1(I,1)
         DO 40 J=1,NGRO
         SCAT(NGRO-J+1,NGRO-I+1,IL+1)=GA2(J,I)
   40    CONTINUE
   45    CONTINUE
      ENDIF
   50 CONTINUE
      CALL XDRLGS(JPLIB,1,IMPX,0,INGRO,1,NGRO,SIGS,SCAT,ITYPRO)
*
      CALL LCMLEN(JPLIB,'TRANC',LENGT,ITYLCM)
      IF (LENGT.GT.0) THEN
         CALL LCMGET(JPLIB,'TRANC',GA1(1,1))
         DO 130 I=1,NGRO
         GA1(I,2)=GA1(NGRO-I+1,1)
  130    CONTINUE
         CALL LCMPUT(JPLIB,'TRANC',NGRO,2,GA1(1,2))
      ENDIF
*
      CALL LCMGET(JPLIB,'NTOT0',GA1(1,1))
      DO 140 I=1,NGRO
      GA1(I,2)=GA1(NGRO-I+1,1)
  140 CONTINUE
      CALL LCMPUT(JPLIB,'NTOT0',NGRO,2,GA1(1,2))
*
      DO 175 IDEL=0,NDEL
      IF(IDEL.EQ.0) THEN
         HNUSIG='NUSIGF'
         HCHI='CHI'
      ELSE
         WRITE(HNUSIG,'(6HNUSIGF,I2.2)') IDEL
         WRITE(HCHI,'(3HCHI,I2.2)') IDEL
      ENDIF
      CALL LCMLEN(JPLIB,HNUSIG,LENGT,ITYLCM)
      IF (LENGT.GT.0) THEN
         CALL LCMGET(JPLIB,HNUSIG,GA1(1,1))
         SUM=0.0
         DO 150 I=1,NGRO
         SUM=SUM+GA1(I,1)
  150    CONTINUE
         DO 160 I=1,NGRO
         GA1(I,2)=GA1(NGRO-I+1,1)/SUM
  160    CONTINUE
         CALL LCMGET(JPLIB,HCHI,GA1(1,1))
         CALL LCMPUT(JPLIB,HCHI,NGRO,2,GA1(1,2))
         DO 170 I=1,NGRO
         GA1(I,2)=GA1(NGRO-I+1,1)*SUM
  170    CONTINUE
         CALL LCMPUT(JPLIB,HNUSIG,NGRO,2,GA1(1,2))
      ENDIF
  175 CONTINUE
*
      DO 190 IED=1,NED
      WRITE(TEXT8,'(2A4)') (NAMEAD(I0,IED),I0=1,2)
      IF((TEXT8.EQ.'TRANC').OR.(TEXT8.EQ.'NTOT0').OR.
     1 (TEXT8(:6).EQ.'NUSIGF').OR.(TEXT8(:3).EQ.'CHI'))
     2 GO TO 190
      CALL LCMLEN(JPLIB,TEXT8,LENGT,ITYLCM)
      IF (LENGT.GT.0) THEN
         CALL LCMGET(JPLIB,TEXT8,GA1(1,1))
         DO 180 I=1,NGRO
         GA1(I,2)=GA1(NGRO-I+1,1)
  180    CONTINUE
         CALL LCMPUT(JPLIB,TEXT8,NGRO,2,GA1(1,2))
      ENDIF
  190 CONTINUE
  200 CONTINUE
*
      CALL LCMGET(IPLIB,'DELTAU',GA1(1,1))
      DO 210 I=1,NGRO
      GA1(I,2)=GA1(NGRO-I+1,1)
  210 CONTINUE
      CALL LCMPUT(IPLIB,'DELTAU',NGRO,2,GA1(1,2))
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(SCAT,SIGS,GA2,GA1,ITYPRO)
      RETURN
      END
