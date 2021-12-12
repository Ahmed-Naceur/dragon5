*DECK LIBNOT
      SUBROUTINE LIBNOT (IPLIB,NGRO,NL,NDIL,NED,NDEL,IMPX,LSCAT,LSIGF,
     1 LADD,DILUT,FLUX,TOTAL,SIGF,SIGS,SCAT,SADD,ZDEL,HVECT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Write the incremental x-s data on a temperature-independant Draglib.
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
* IPLIB   pointer to the internal library (L_LIBRARY signature).
* NGRO    number of energy groups.
* NL      number of Legendre orders required in the calculation
*         (NL=1 or higher).
* NDIL    number of finite dilutions.
* NED     number of extra vector edits.
* NDEL    number of delayed neutron precursor groups.
* IMPX    print flag.
* LSCAT   Legendre flag (=.true. if a given Legendre order of the
*         scattering cross section exists).
* LSIGF   fission flag (=.true. if the isotope can fission).
* LADD    additional xs flag (=.true. if a given additional cross
*         section exists).
* DILUT   dilutions.
* FLUX    weighting flux.
* TOTAL   total cross sections.
* SIGF    nu*fission cross sections.
* SIGS    diffusion cross sections.
* SCAT    scattering transfer matrices (sec,prim,Legendre,dilution).
* SADD    additional cross sections.
* ZDEL    delayed nu-sigf cross sections.
* HVECT   names of the extra vector edits.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB
      INTEGER NGRO,NL,NDIL,NED,NDEL,IMPX
      REAL DILUT(NDIL+1),FLUX(NGRO,NDIL+1),TOTAL(NGRO,NDIL+1),
     1 SIGF(NGRO,NDIL+1),SIGS(NGRO,NL,NDIL+1),SCAT(NGRO,NGRO,NL,NDIL+1),
     2 SADD(NGRO,NED,NDIL+1),ZDEL(NGRO,NDEL,NDIL+1)
      LOGICAL LSIGF,LSCAT(NL),LADD(NED)
      CHARACTER HVECT(NED)*8
*----
*  LOCAL VARIABLES
*----
      CHARACTER TEXT12*12,CD*4
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ITYPRO
      REAL, ALLOCATABLE, DIMENSION(:) :: GAS
      REAL, ALLOCATABLE, DIMENSION(:,:) :: GA1
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: GA2
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(ITYPRO(NL))
      ALLOCATE(GAS(NGRO),GA1(NGRO,NL),GA2(NGRO,NGRO,NL))
*
      CALL LCMPUT(IPLIB,'DILUTION',NDIL,2,DILUT)
      DO 130 IDIL=1,NDIL
      WRITE (CD,'(I4.4)') IDIL
      CALL LCMSIX(IPLIB,'SUBMAT'//CD,1)
      DO 10 IG1=1,NGRO
      GAS(IG1)=FLUX(IG1,IDIL)-1.0
   10 CONTINUE
      CALL LCMPUT(IPLIB,'NWT0',NGRO,2,GAS)
      DO 20 IG1=1,NGRO
      GAS(IG1)=TOTAL(IG1,IDIL)*FLUX(IG1,IDIL)-TOTAL(IG1,NDIL+1)
   20 CONTINUE
      CALL LCMPUT(IPLIB,'NTOT0',NGRO,2,GAS)
      IF(LSIGF) THEN
         DO 30 IG1=1,NGRO
         GAS(IG1)=SIGF(IG1,IDIL)*FLUX(IG1,IDIL)-SIGF(IG1,NDIL+1)
   30    CONTINUE
         CALL LCMPUT(IPLIB,'NUSIGF',NGRO,2,GAS)
      ENDIF
      INGRO=NL-1
      DO 40 IL=NL-1,0,-1
      IF(.NOT.LSCAT(IL+1)) THEN
         INGRO=INGRO-1
      ELSE
         GO TO 50
      ENDIF
   40 CONTINUE
   50 DO 80 IL=1,INGRO+1
      IF(LSCAT(IL)) THEN
         DO 65 IG1=1,NGRO
         GA1(IG1,IL)=SIGS(IG1,IL,IDIL)*FLUX(IG1,IDIL)-
     1   SIGS(IG1,IL,NDIL+1)
         DO 60 IG2=1,NGRO
         GA2(IG1,IG2,IL)=SCAT(IG1,IG2,IL,IDIL)*FLUX(IG2,IDIL)-
     1   SCAT(IG1,IG2,IL,NDIL+1)
   60    CONTINUE
   65    CONTINUE
      ELSE
         DO 75 IG1=1,NGRO
         GA1(IG1,IL)=0.0
         DO 70 IG2=1,NGRO
         GA2(IG1,IG2,IL)=0.0
   70    CONTINUE
   75    CONTINUE
      ENDIF
   80 CONTINUE
      CALL XDRLGS(IPLIB,1,IMPX,0,INGRO,1,NGRO,GA1,GA2,ITYPRO)
      DO 100 IED=1,NED
      IF(LADD(IED)) THEN
         DO 90 IG1=1,NGRO
         GAS(IG1)=SADD(IG1,IED,IDIL)*FLUX(IG1,IDIL)-SADD(IG1,IED,NDIL+1)
   90    CONTINUE
         CALL LCMPUT(IPLIB,HVECT(IED),NGRO,2,GAS)
      ENDIF
  100 CONTINUE
      DO 120 IDEL=1,NDEL
      WRITE(TEXT12,'(6HNUSIGF,I2.2)') IDEL
      DO 110 IG1=1,NGRO
      GAS(IG1)=ZDEL(IG1,IDEL,IDIL)*FLUX(IG1,IDIL)-ZDEL(IG1,IDEL,NDIL+1)
  110 CONTINUE
      CALL LCMPUT(IPLIB,TEXT12,NGRO,2,GAS)
  120 CONTINUE
      CALL LCMSIX(IPLIB,' ',2)
  130 CONTINUE
      IF(IMPX.GT.3) CALL LCMLIB(IPLIB)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(GA2,GA1,GAS)
      DEALLOCATE(ITYPRO)
      RETURN
      END
