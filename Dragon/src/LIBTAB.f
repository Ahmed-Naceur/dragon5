*DECK LIBTAB
      SUBROUTINE LIBTAB (IGRP,NGRO,NL,NDIL,NPART,NED,NDEL,HNAMIS,IMPX,
     1 LSCAT,LSIGF,LADD,DILUT,TOTAL,SIGF,SIGS,SCAT,SADD,ZDEL,GOLD,ISMIN,
     2 ISMAX,NOR,SIGP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Transform dilution dependent information into probability tables.
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
* IGRP    energy group index where the transformation occurs.
* NGRO    number of energy groups.
* NL      number of Legendre orders required in the calculation
*         (NL=1 or higher).
* NDIL    number of finite dilutions.
* NPART   2 + number of partial cross sections.
* NED     number of extra vector edits.
* NDEL    number of delayed neutron precursor groups.
* HNAMIS  local name of the isotope:
*         HNAMIS(1:8)  the local isotope name;
*         HNAMIS(9:12) suffix function of the mix number.
* IMPX    print flag.
* LSCAT   Legendre flag (=.true. if a given Legendre order of the
*         scattering cross section exists).
* LSIGF   fission flag (=.true. if the isotope can fission).
* LADD    additional xs flag (=.true. if a given additional cross
*         section exists).
* DILUT   dilutions.
* TOTAL   total cross sections.
* SIGF    nu*fission cross sections.
* SIGS    scattering cross sections.
* SCAT    scattering transfer matrices (sec,prim,Legendre,dilution).
* SADD    additional cross sections.
* ZDEL    delayed nu-sigf cross sections.
* GOLD    Goldstein-Cohen parameter.
* ISMIN   minimum secondary group corresponding to each primary group.
* ISMAX   maximum secondary group corresponding to each primary group.
*
*Parameters: output
* NOR     order of the probability table.
* SIGP    partial cross sections.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      PARAMETER (MAXNOR=12)
      INTEGER IGRP,NGRO,NL,NDIL,NPART,NED,NDEL,IMPX,ISMIN(NL,NGRO),
     1 ISMAX(NL,NGRO),NOR
      REAL DILUT(NDIL+1),TOTAL(NGRO,NDIL+1),SIGF(NGRO,NDIL+1),
     1 SIGS(NGRO,NL,NDIL+1),SCAT(NGRO,NGRO,NL,NDIL+1),
     2 SADD(NGRO,NED,NDIL+1),ZDEL(NGRO,NDEL,NDIL+1),GOLD,
     3 SIGP(MAXNOR,NPART)
      LOGICAL LSIGF,LSCAT(NL),LADD(NED)
      CHARACTER HNAMIS*12
*----
*  LOCAL VARIABLES
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: DILUT2,XSDIL
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: LDIL
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(LDIL(NDIL+1),DILUT2(NDIL+1))
*----
*  REMOVE BADLY BEHAVED COLLOCATIONS POINTS.
*----
      MDIL=NDIL
      DO 10 IDIL=1,NDIL+1
      LDIL(IDIL)=.TRUE.
   10 CONTINUE
      TEST=TOTAL(IGRP,NDIL+1)
      DO 20 IDIL=NDIL,1,-1
      IF(ABS(TOTAL(IGRP,IDIL)-TEST).LE.2.0E-4*ABS(TEST)) THEN
         MDIL=MDIL-1
         LDIL(IDIL)=.FALSE.
      ELSE IF(DILUT(IDIL).LT.1.0) THEN
         MDIL=MDIL-1
         LDIL(IDIL)=.FALSE.
      ELSE IF((DILUT(IDIL).GT.1.0E5).AND.(DILUT(IDIL).LT.1.0E10)) THEN
         MDIL=MDIL-1
         LDIL(IDIL)=.FALSE.
      ELSE IF(TOTAL(IGRP,IDIL).LE.0.0) THEN
         MDIL=MDIL-1
         LDIL(IDIL)=.FALSE.
      ELSE IF(TOTAL(IGRP,IDIL)-(1.0-GOLD)*SIGS(IGRP,1,IDIL).LE.0.0) THEN
         MDIL=MDIL-1
         LDIL(IDIL)=.FALSE.
      ELSE
         TEST=TOTAL(IGRP,IDIL)
      ENDIF
   20 CONTINUE
*
      ALLOCATE(XSDIL((MDIL+1)*(NPART+1)))
      IOFSET=-1
      IDD=0
      DO 30 IDIL=1,NDIL+1
      IF(LDIL(IDIL)) THEN
         IDD=IDD+1
         DILUT2(IDD)=DILUT(IDIL)
         IOFSET=IOFSET+1
         XSDIL(IOFSET+1)=TOTAL(IGRP,IDIL)
      ENDIF
   30 CONTINUE
      IF(IDD.NE.MDIL+1) CALL XABORT('LIBTAB: INTERNAL ERROR.')
      DO 40 IDIL=1,NDIL+1
      IF(LDIL(IDIL)) THEN
         IOFSET=IOFSET+1
         IF(LSIGF) THEN
            XSDIL(IOFSET+1)=SIGF(IGRP,IDIL)
         ELSE
            XSDIL(IOFSET+1)=0.0
         ENDIF
      ENDIF
   40 CONTINUE
      DO 55 IL=1,NL
      DO 50 IDIL=1,NDIL+1
      IF(LDIL(IDIL)) THEN
         IOFSET=IOFSET+1
         IF(LSCAT(IL)) THEN
            XSDIL(IOFSET+1)=SIGS(IGRP,IL,IDIL)
         ELSE
            XSDIL(IOFSET+1)=0.0
         ENDIF
      ENDIF
   50 CONTINUE
   55 CONTINUE
      IF(NPART.EQ.3+NL) GO TO 100
      DO 70 IL=1,NL
      IF(LSCAT(IL)) THEN
         DO 65 IG2=ISMIN(IL,IGRP),ISMAX(IL,IGRP)
         DO 60 IDIL=1,NDIL+1
         IF(LDIL(IDIL)) THEN
            IOFSET=IOFSET+1
            XSDIL(IOFSET+1)=SCAT(IG2,IGRP,IL,IDIL)
         ENDIF
   60    CONTINUE
   65    CONTINUE
      ENDIF
   70 CONTINUE
      DO 85 IED=1,NED
      DO 80 IDIL=1,NDIL+1
      IF(LDIL(IDIL)) THEN
         IOFSET=IOFSET+1
         IF(LADD(IED)) THEN
            XSDIL(IOFSET+1)=SADD(IGRP,IED,IDIL)
         ELSE
            XSDIL(IOFSET+1)=0.0
         ENDIF
      ENDIF
   80 CONTINUE
   85 CONTINUE
      DO 95 IDEL=1,NDEL
      DO 90 IDIL=1,NDIL+1
      IF(LDIL(IDIL)) THEN
         IOFSET=IOFSET+1
         XSDIL(IOFSET+1)=ZDEL(IGRP,IDEL,IDIL)
      ENDIF
   90 CONTINUE
   95 CONTINUE
*
  100 DO 115 IPART=1,NPART
      DO 110 INOR=1,MAXNOR
      SIGP(INOR,IPART)=0.0
  110 CONTINUE
  115 CONTINUE
      CALL LIBPTT(IGRP,MDIL,NPART-2,DILUT2,XSDIL,GOLD,HNAMIS,
     1 IMPX,NOR,SIGP(1,1),SIGP(1,2),SIGP(1,3))
*
      DEALLOCATE(XSDIL)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(DILUT2,LDIL)
      RETURN
      END
