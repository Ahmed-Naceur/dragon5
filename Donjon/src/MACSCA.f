*DECK MACSCA
      SUBROUTINE MACSCA(KPOLD,KPNEW,SCAT,SCAT2,CM,JGR,IL,MIX,NMXNEW,
     1 NTOT,NMXOLD,NL,NGRP,LMAP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover scattering matrices and store them in a new macrolib for
* a given anistropic level and energy group.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* J. Koclas, E. Varin, D. Sekki
*
*Parameters: input
* KPOLD   pointer to group directory in the initial macrolib.
* NL      number of legendre orders (=1 for isotropic scattering).
* NGRP    number of energy groups.
* NMXOLD  number of material mixtures in the initial macrolib.
* NMXNEW  number of material mixtures in the final macrolib.
* MIX     index of all (material and virtual) mixtures per region.
* NTOT    total number of all (material and virtual) mixtures.
* SCAT    scattering matrices in the initial macrolib.
* SCAT2   scattering matrices in the final macrolib.
* IL      anisotropic level to be treated.
* JGR     energy group to be treated.
* CM      anisotropic level in I2.2 format.
* LMAP    flag for the initial macrolib:
*          =.true. if the fuel map macrolib.
*
*Parameters: output
* KPNEW   pointer to group directory in the final macrolib.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) KPOLD,KPNEW
      REAL SCAT(NMXOLD,NL,NGRP,NGRP),SCAT2(NMXNEW,NL,NGRP,NGRP)
      INTEGER MIX(NTOT)
      CHARACTER CM*2
      LOGICAL LMAP
*----
*  LOCAL VARIABLES
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IPOS,IPOS2
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: IJJ,IJJ2,NJJ,NJJ2
      REAL, ALLOCATABLE, DIMENSION(:) :: WORK,WORK2
      CHARACTER HSMG*131
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IPOS(NMXOLD),IPOS2(NMXNEW),IJJ(NMXOLD,NL,NGRP),
     1 IJJ2(NMXNEW,NL,NGRP),NJJ(NMXOLD,NL,NGRP),NJJ2(NMXNEW,NL,NGRP))
      ALLOCATE(WORK(NMXOLD*NGRP),WORK2(NMXNEW*NGRP))
      CALL XDRSET(WORK,NMXOLD*NGRP,0.)
      CALL XDRSET(WORK2,NMXNEW*NGRP,0.)
*----
*  RECOVER EXISTING DATA
*----
      CALL LCMLEN(KPNEW,'NJJS'//CM,ILENG,ITYP)
      IF(LMAP.AND.(ILENG.GT.0))THEN
        IF(ILENG.NE.NMXNEW)CALL XABORT('@MACSCA: INVALID MACROLIB(1).')
        CALL LCMGET(KPNEW,'SCAT'//CM,WORK2(1))
        CALL LCMGET(KPNEW,'NJJS'//CM,NJJ2(1,IL,JGR))
        CALL LCMGET(KPNEW,'IJJS'//CM,IJJ2(1,IL,JGR))
        CALL LCMGET(KPNEW,'IPOS'//CM,IPOS2(1))
        DO 15 IBM=1,NMXNEW
          IJJ0=IJJ2(IBM,IL,JGR)
          IPOSDE=IPOS2(IBM)
          DO 10 IGR=IJJ0,IJJ0-NJJ2(IBM,IL,JGR)+1,-1
            SCAT2(IBM,IL,IGR,JGR)=WORK2(IPOSDE)
            IPOSDE=IPOSDE+1
   10     CONTINUE
   15   CONTINUE
      ENDIF
*----
*  RECOVER SCAT,IJJ,NJJ,IPOS
*----
      CALL LCMLEN(KPOLD,'NJJS'//CM,ILENG,ITYP)
      IF(ILENG.EQ.0)CALL XABORT('@MACSCA: INVALID MACROLIB(2).')
      CALL LCMGET(KPOLD,'SCAT'//CM,WORK(1))
      CALL LCMGET(KPOLD,'NJJS'//CM,NJJ(1,IL,JGR))
      CALL LCMGET(KPOLD,'IJJS'//CM,IJJ(1,IL,JGR))
      CALL LCMGET(KPOLD,'IPOS'//CM,IPOS(1))
      DO 25 IBM=1,NMXOLD
        IJJ0=IJJ(IBM,IL,JGR)
        IPOSDE=IPOS(IBM)
        DO 20 IGR=IJJ0,IJJ0-NJJ(IBM,IL,JGR)+1,-1
          SCAT(IBM,IL,IGR,JGR)=WORK(IPOSDE)
          IPOSDE=IPOSDE+1
   20   CONTINUE
   25 CONTINUE
*----
*  NEW SCAT2
*----
      ITOT=0
      DO 50 IBM=1,NTOT
      IF(MIX(IBM).EQ.0)GOTO 50
      ITOT=ITOT+1
      IF(LMAP)THEN
*     ONLY FUEL DATA WILL BE COPIED
        IF(MIX(IBM).GT.0)GOTO 50
        J=-MIX(IBM)
        IF(J.GT.NMXOLD) THEN
          WRITE(HSMG,'(25HMACSCA: A MIXTURE INDEX (,I6,12H) IS GREATER,
     >    36H THAN THE TOTAL NUMBER OF MIXTURES (,I6,14H) IN 2ND RHS M,
     >    8HACROLIB.)') J,NMXOLD
          CALL XABORT(HSMG)
        ENDIF
      ELSE
*     FUEL DATA WILL NOT BE COPIED
        IF(MIX(IBM).LT.0)GOTO 50
        J=MIX(IBM)
        IF(J.GT.NMXOLD) THEN
          WRITE(HSMG,'(25HMACSCA: A MIXTURE INDEX (,I6,12H) IS GREATER,
     >    36H THAN THE TOTAL NUMBER OF MIXTURES (,I6,14H) IN 1ST RHS M,
     >    8HACROLIB.)') J,NMXOLD
          CALL XABORT(HSMG)
        ENDIF
      ENDIF
*     COPY DATA
      IJJ0=IJJ(J,IL,JGR)
      DO 40 IGR=IJJ0,IJJ0-NJJ(J,IL,JGR)+1,-1
        SCAT2(ITOT,IL,IGR,JGR)=SCAT(J,IL,IGR,JGR)
   40 CONTINUE
   50 CONTINUE
*----
*  NEW IJJ2 AND NJJ2
*----
      DO 70 IBM=1,NMXNEW
        IGMIN=JGR
        IGMAX=JGR
        DO 60 IGR=NGRP,1,-1
        IF(SCAT2(IBM,IL,IGR,JGR).NE.0.)THEN
          IGMIN=MIN(IGMIN,IGR)
          IGMAX=MAX(IGMAX,IGR)
        ENDIF
   60   CONTINUE
      IJJ2(IBM,IL,JGR)=IGMAX
      NJJ2(IBM,IL,JGR)=IGMAX-IGMIN+1
   70 CONTINUE
*----
*  STORE SCAT2,IJJ2,NJJ2,IPOS2
*----
      IPOSDE=0
      DO 85 IBM=1,NMXNEW
        IPOS2(IBM)=IPOSDE+1
        DO 80 IGR=IJJ2(IBM,IL,JGR),IJJ2(IBM,IL,JGR)-
     1            NJJ2(IBM,IL,JGR)+1,-1
          IPOSDE=IPOSDE+1
          WORK2(IPOSDE)=SCAT2(IBM,IL,IGR,JGR)
   80   CONTINUE
   85 CONTINUE
      CALL LCMPUT(KPNEW,'SCAT'//CM,IPOSDE,2,WORK2)
      CALL LCMPUT(KPNEW,'IPOS'//CM,NMXNEW,1,IPOS2)
      CALL LCMPUT(KPNEW,'NJJS'//CM,NMXNEW,1,NJJ2(1,IL,JGR))
      CALL LCMPUT(KPNEW,'IJJS'//CM,NMXNEW,1,IJJ2(1,IL,JGR))
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(WORK2,WORK)
      DEALLOCATE(NJJ2,NJJ,IJJ2,IJJ,IPOS2,IPOS)
      RETURN
      END
