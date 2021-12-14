*DECK LIBEXT
      SUBROUTINE LIBEXT (IPDRL,NGRO,NL,NDIL,NED,HVECT,NDEL,IMPX,DILUT,
     1 MDIL,LSCAT,LSIGF,LADD,LGOLD,FLUX,TOTAL,SIGF,SIGS,SCAT,SADD,ZDEL,
     2 DELTG,GOLD,ISMIN,ISMAX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read dilution-dependent information of one isotope in multi-dilution
* internal library format.
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
* IPDRL   pointer to the multi-dilution internal library.
* NGRO    number of energy groups.
* NL      number of Legendre orders required in the calculation
*         (NL=1 or higher).
* NDIL    number of finite dilutions.
* NED     number of extra vector edits.
* HVECT   names of the extra vector edits.
* NDEL    number of delayed neutron precursor groups.
* IMPX    print flag.
*
*Parameters: input/output
* DILUT   dilutions.
*
*Parameters: output
* MDIL    number of finite dilutions used.
* LSCAT   Legendre flag (=.true. if a given Legendre order of the
*         scattering cross section exists).
* LSIGF   fission flag (=.true. if the isotope can fission).
* LADD    additional xs flag (=.true. if a given additional cross
*         section exists).
* LGOLD   Goldstein-Cohen flag (=.true. if Goldstein-Cohen parameters
*         exists).
* FLUX    weighting flux.
* TOTAL   total cross sections.
* SIGF    nu*fission cross sections.
* SIGS    scattering cross sections.
* SCAT    scattering transfer matrices (sec,prim,Legendre,dilution).
* SADD    additional cross sections.
* ZDEL    delayed nu-sigf cross sections.
* DELTG   lethargy widths.
* GOLD    Goldstein-Cohen parameters.
* ISMIN   minimum secondary group corresponding to each primary group.
* ISMAX   maximum secondary group corresponding to each primary group.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPDRL
      INTEGER NGRO,NL,NDIL,NED,NDEL,IMPX,MDIL,ISMIN(NL,NGRO),
     1 ISMAX(NL,NGRO)
      REAL DILUT(NDIL+1),FLUX(NGRO,NDIL+1),TOTAL(NGRO,NDIL+1),
     1 SIGF(NGRO,NDIL+1),SIGS(NGRO,NL,NDIL+1),SCAT(NGRO,NGRO,NL,NDIL+1),
     2 SADD(NGRO,NED,NDIL+1),ZDEL(NGRO,NDEL,NDIL+1),DELTG(NGRO),
     3 GOLD(NGRO)
      CHARACTER HVECT(NED)*8
      LOGICAL LSIGF,LSCAT(NL),LADD(NED),LGOLD
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ITYPRO,IPDIL
*----
*  LOCAL VARIABLES
*----
      PARAMETER(MAXTIT=10)
      TYPE(C_PTR) JPDRL,KPDRL
      CHARACTER TEXNUD*12
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(ITYPRO(NL),IPDIL(NDIL+1))
*
      DO 10 IL=1,NL
      LSCAT(IL)=.FALSE.
   10 CONTINUE
      LSIGF=.FALSE.
      DO 20 IED=1,NED
      LADD(IED)=.FALSE.
   20 CONTINUE
      CALL LCMGET(IPDRL,'DELTAU',DELTG)
*----
*  RECOVER DILUTION-DEPENDENT VALUES.
*----
      JPDRL=LCMGID(IPDRL,'ISOTOPESLIST')
      DO 80 IDIL=1,NDIL+1
      KPDRL=LCMGIL(JPDRL,IDIL) ! set IDIL-th isotope
      CALL LCMGET(KPDRL,'NWT0',FLUX(1,IDIL))
      CALL LCMGET(KPDRL,'NTOT0',TOTAL(1,IDIL))
      CALL LCMLEN(KPDRL,'NUSIGF',LENGT,ITYLCM)
      LSIGF=LSIGF.OR.(LENGT.GT.0)
      IF(LENGT.GT.0) THEN
         CALL LCMGET(KPDRL,'NUSIGF',SIGF(1,IDIL))
      ELSE
         CALL XDRSET(SIGF(1,IDIL),NGRO,0.0)
      ENDIF
      CALL XDRLGS(KPDRL,-1,IMPX,0,NL-1,1,NGRO,SIGS(1,1,IDIL),
     1 SCAT(1,1,1,IDIL),ITYPRO)
      DO 30 IL=0,NL-1
      LSCAT(IL+1)=LSCAT(IL+1).OR.(ITYPRO(IL+1).GT.0)
   30 CONTINUE
      DO 50 IED=1,NED
      DO 40 IG1=1,NGRO
      SADD(IG1,IED,IDIL)=0.0
   40 CONTINUE
      CALL LCMLEN(KPDRL,HVECT(IED),LENGT,ITYLCM)
      LADD(IED)=LADD(IED).OR.(LENGT.GT.0)
      IF(LENGT.GT.0) CALL LCMGET(KPDRL,HVECT(IED),SADD(1,IED,IDIL))
   50 CONTINUE
      DO 70 IDEL=1,NDEL
      WRITE(TEXNUD,'(6HNUSIGF,I2.2)') IDEL
      DO 60 IG1=1,NGRO
      ZDEL(IG1,IDEL,IDIL)=0.0
   60 CONTINUE
      CALL LCMLEN(KPDRL,TEXNUD,LENGT,ITYLCM)
      IF(LENGT.GT.0) CALL LCMGET(KPDRL,TEXNUD,ZDEL(1,IDEL,IDIL))
   70 CONTINUE
      IF(IDIL.EQ.NDIL+1) THEN
         CALL LCMLEN(KPDRL,'NGOLD',LENGT,ITYLCM)
         LGOLD=LENGT.GT.0
         IF(LGOLD) THEN
            CALL LCMGET(KPDRL,'NGOLD',GOLD)
         ELSE
            CALL XDRSET(GOLD,NGRO,1.0)
         ENDIF
      ENDIF
   80 CONTINUE
*----
*  SET THE SIGNIFICANT DILUTIONS.
*----
      MDIL=0
      DO 90 IDIL=1,NDIL
      IF(DILUT(IDIL).LT.1.5) THEN
         CONTINUE
      ELSE IF((DILUT(IDIL).GT.1.0E5).AND.(DILUT(IDIL).LT.1.0E10)) THEN
         CONTINUE
      ELSE
         MDIL=MDIL+1
         IPDIL(MDIL)=IDIL
      ENDIF
   90 CONTINUE
      IPDIL(MDIL+1)=NDIL+1
      DO 122 IDIL=1,MDIL+1
      DILUT(IDIL)=DILUT(IPDIL(IDIL))
      DO 121 IG1=1,NGRO
      FLUX(IG1,IDIL)=FLUX(IG1,IPDIL(IDIL))
      TOTAL(IG1,IDIL)=TOTAL(IG1,IPDIL(IDIL))
      SIGF(IG1,IDIL)=SIGF(IG1,IPDIL(IDIL))
      DO 105 IL=1,NL
      SIGS(IG1,IL,IDIL)=SIGS(IG1,IL,IPDIL(IDIL))
      DO 100 IG2=1,NGRO
      SCAT(IG2,IG1,IL,IDIL)=SCAT(IG2,IG1,IL,IPDIL(IDIL))
  100 CONTINUE
  105 CONTINUE
      DO 110 IED=1,NED
      SADD(IG1,IED,IDIL)=SADD(IG1,IED,IPDIL(IDIL))
  110 CONTINUE
      DO 120 IDEL=1,NDEL
      ZDEL(IG1,IDEL,IDIL)=ZDEL(IG1,IDEL,IPDIL(IDIL))
  120 CONTINUE
  121 CONTINUE
  122 CONTINUE
*----
*  COMPUTE THE SCATTERING BANDWIDTH AND MOST THERMAL GROUPS.
*----
      DO 160 IL=1,NL
      IF(LSCAT(IL)) THEN
         DO 130 IG1=1,NGRO
         ISMIN(IL,IG1)=NGRO
         ISMAX(IL,IG1)=1
  130    CONTINUE
         DO 142 IG2=1,NGRO
         DO 141 IDIL=1,MDIL+1
         DO 140 IG1=NGRO,1,-1
         IF(SCAT(IG2,IG1,IL,IDIL).NE.0.0) THEN
            ISMIN(IL,IG1)=MIN(ISMIN(IL,IG1),IG2)
            ISMAX(IL,IG1)=MAX(ISMAX(IL,IG1),IG2)
         ENDIF
  140    CONTINUE
  141    CONTINUE
  142    CONTINUE
      ELSE
         DO 150 IG1=1,NGRO
         ISMIN(IL,IG1)=NGRO+1
         ISMAX(IL,IG1)=0
  150    CONTINUE
      ENDIF
  160 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IPDIL,ITYPRO)
      RETURN
      END
