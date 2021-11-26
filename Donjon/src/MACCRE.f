*DECK MACCRE
      SUBROUTINE MACCRE(IPOLD,IPNEW,NL,NW,NF,NGRP,NMXOLD,NMXNEW,NTOT,
     1 MIX,LMAP,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover nuclear properties from an initial macrolib and store them
* in a new one containing one mixture per region.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* J. Koclas, E. Varin, D. Sekki
*
*Parameters: input
* IPOLD   pointer to the initial macrolib.
* NL      number of legendre orders (=1 for isotropic scattering).
* NW      legendre order of NWT information (=0: NTOT0; =1: NTOT1).
* NF      number of fissile isotopes.
* NGRP    number of energy groups.
* NMXOLD  number of material mixtures in the initial macrolib.
* NMXNEW  number of material mixtures in the final macrolib.
* NTOT    total number of all (material and virtual) mixtures.
* MIX     index of all (material and virtual) mixtures.
* LMAP    flag for the initial macrolib:
*          =.true. if the fuel map macrolib.
* IMPX    printing index (=0 for no print).
*
*Parameters: output
* IPNEW   pointer to the final macrolib.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPOLD,IPNEW
      INTEGER NL,NW,NF,NGRP,NMXOLD,NMXNEW,NTOT,MIX(NTOT)
      LOGICAL LMAP
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IOUT=6)
      CHARACTER CM*2,NAME*12,FIRST*12
      TYPE(C_PTR) JPOLD,JPNEW,KPOLD,KPNEW
      REAL, ALLOCATABLE, DIMENSION(:) ::SCAT,SCAT2,DATA,DATA2
*
      ALLOCATE(SCAT(NMXOLD*NL*NGRP*NGRP))
      ALLOCATE(SCAT2(NMXNEW*NL*NGRP*NGRP))
      CALL XDRSET(SCAT,NMXOLD*NL*NGRP*NGRP,0.)
      CALL XDRSET(SCAT2,NMXNEW*NL*NGRP*NGRP,0.)
*----
*  RECOVER MACROLIB DATA
*----
      JPOLD=LCMGID(IPOLD,'GROUP')
      JPNEW=LCMLID(IPNEW,'GROUP',NGRP)
      DO 100 JGR=1,NGRP
      KPOLD=LCMGIL(JPOLD,JGR)
      KPNEW=LCMDIL(JPNEW,JGR)
      IF(IMPX.GT.3)CALL LCMLIB(KPOLD)
      IF(IMPX.GT.2)WRITE(IOUT,*)'** TREATING ENERGY GROUP #',JGR
      NAME=' '
      CALL LCMNXT(KPOLD,NAME)
      FIRST=NAME
   10 CALL LCMLEN(KPOLD,NAME,LENGT,ITYP)
      IF((INDEX(NAME,'NTOT0').EQ.1).OR.(INDEX(NAME,'DIF').EQ.1).OR.
     1   (INDEX(NAME,'NFT').EQ.1).OR.(INDEX(NAME,'OVE').EQ.1).OR.
     2   (INDEX(NAME,'H-F').EQ.1).OR.(INDEX(NAME,'SIG').EQ.1))THEN
*     RECOVER THESE PROPERTIES
        IF(IMPX.GT.2)WRITE(IOUT,*)'PROPERTY NAME : ',NAME
        IF(LENGT.EQ.NMXOLD)THEN
          ALLOCATE(DATA(NMXOLD),DATA2(NMXNEW))
          CALL XDRSET(DATA,NMXOLD,0.)
          CALL LCMGET(KPOLD,NAME,DATA)
          CALL XDRSET(DATA2,NMXNEW,0.)
          IF(LMAP)THEN
*         RECOVER EXISTING DATA
            CALL LCMLEN(KPNEW,NAME,LENGT2,ITYP2)
            IF(LENGT2.NE.0)CALL LCMGET(KPNEW,NAME,DATA2)
          ENDIF
          ITOT=0
          DO 20 IBM=1,NTOT
          IF(MIX(IBM).EQ.0)GOTO 20
          ITOT=ITOT+1
          IF(LMAP)THEN
*         ONLY FUEL DATA WILL BE COPIED
            IF(MIX(IBM).GT.0)GOTO 20
            J=-MIX(IBM)
          ELSE
*         FUEL DATA WILL NOT BE COPIED
            IF(MIX(IBM).LT.0)GOTO 20
            J=MIX(IBM)
          ENDIF
*         COPY DATA
          DATA2(ITOT)=DATA(J)
   20     CONTINUE
*         STORE DATA
          CALL LCMPUT(KPNEW,NAME,NMXNEW,ITYP,DATA2)
          DEALLOCATE(DATA,DATA2)
        ELSEIF(LENGT.EQ.-1)THEN
          CALL XABORT('@MACCRE: '//NAME//' IS A DIRECTORY.')
        ELSEIF(LENGT.NE.0)THEN
          CALL XABORT('@MACCRE: INVALID INPUT MACROLIB(1).')
        ENDIF
      ELSE IF((INDEX(NAME,'NUS').EQ.1).OR.(INDEX(NAME,'CHI').EQ.1))THEN
*       RECOVER FISSION-RELATED PROPERTIES
        IF(IMPX.GT.2)WRITE(IOUT,*)'PROPERTY NAME : ',NAME
        IF(LENGT.EQ.NMXOLD*NF)THEN
          ALLOCATE(DATA(NMXOLD*NF),DATA2(NMXNEW*NF))
          CALL XDRSET(DATA,NMXOLD*NF,0.)
          CALL LCMGET(KPOLD,NAME,DATA)
          CALL XDRSET(DATA2,NMXNEW*NF,0.)
          IF(LMAP)THEN
*         RECOVER EXISTING DATA
            CALL LCMLEN(KPNEW,NAME,LENGT2,ITYP2)
            IF(LENGT2.NE.0)CALL LCMGET(KPNEW,NAME,DATA2)
          ENDIF
          ITOT=0
          DO 35 INF=1,NF
          DO 30 IBM=1,NTOT
          IF(MIX(IBM).EQ.0)GOTO 30
          ITOT=ITOT+1
          IF(LMAP)THEN
*         ONLY FUEL DATA WILL BE COPIED
            IF(MIX(IBM).GT.0)GOTO 30
            J=-MIX(IBM)
          ELSE
*         FUEL DATA WILL NOT BE COPIED
            IF(MIX(IBM).LT.0)GOTO 30
            J=MIX(IBM)
          ENDIF
*         COPY DATA
          J1=(INF-1)*NMXOLD+J
          DATA2(ITOT)=DATA(J1)
   30     CONTINUE
   35     CONTINUE
*         STORE DATA
          CALL LCMPUT(KPNEW,NAME,NMXNEW*NF,ITYP,DATA2)
          DEALLOCATE(DATA,DATA2)
        ELSEIF(LENGT.EQ.-1)THEN
          CALL XABORT('@MACCRE: '//NAME//' IS A DIRECTORY.')
        ELSEIF(LENGT.NE.0)THEN
          CALL XABORT('@MACCRE: INVALID INPUT MACROLIB(2).')
        ENDIF
      ENDIF
      CALL LCMNXT(KPOLD,NAME)
      IF(FIRST.EQ.NAME)GOTO 40
      GOTO 10
*     RECOVER SCAT,IJJ,NJJ,IPOS
   40 IF(IMPX.GT.2)WRITE(IOUT,*)'RECOVERING OF SCAT,IJJ,NJJ,IPOS'
      DO IL=1,NL
        WRITE (CM,'(I2.2)') IL-1
        CALL LCMLEN(KPOLD,'SCAT'//CM,LENGT,ITYP)
        IF(LENGT.EQ.0)THEN
          EXIT
        ELSEIF(LENGT.GT.NMXOLD*NL*NGRP*NGRP)THEN
          CALL XABORT('@MACCRE: INVALID INPUT MACROLIB(3).')
        ELSEIF(LENGT.GT.0)THEN
          CALL MACSCA(KPOLD,KPNEW,SCAT,SCAT2,CM,JGR,IL,MIX,NMXNEW,NTOT,
     1    NMXOLD,NL,NGRP,LMAP)
        ENDIF
      ENDDO
*     RECOVER NTOT1 information
      IF(NW.GT.0) THEN
        ALLOCATE(DATA(NMXOLD),DATA2(NMXNEW))
        CALL XDRSET(DATA,NMXOLD,0.)
        CALL LCMGET(KPOLD,'NTOT1',DATA)
        CALL XDRSET(DATA2,NMXNEW,0.)
        IF(LMAP)THEN
*       RECOVER EXISTING DATA
          CALL LCMLEN(KPNEW,'NTOT0',LENGT1,ITYP1)
          CALL LCMLEN(KPNEW,'NTOT1',LENGT2,ITYP2)
          IF(LENGT2.NE.0) THEN
            CALL LCMGET(KPNEW,'NTOT1',DATA2)
          ELSE IF(LENGT1.NE.0) THEN
            CALL LCMGET(KPNEW,'NTOT0',DATA2)
          ENDIF
        ENDIF
        ITOT=0
        DO 50 IBM=1,NTOT
        IF(MIX(IBM).EQ.0)GOTO 50
        ITOT=ITOT+1
        IF(LMAP)THEN
*       ONLY FUEL DATA WILL BE COPIED
          IF(MIX(IBM).GT.0)GOTO 50
          J=-MIX(IBM)
        ELSE
*       FUEL DATA WILL NOT BE COPIED
          IF(MIX(IBM).LT.0)GOTO 50
          J=MIX(IBM)
        ENDIF
*       COPY DATA
        DATA2(ITOT)=DATA(J)
   50   CONTINUE
*       STORE DATA
        CALL LCMPUT(KPNEW,'NTOT1',NMXNEW,ITYP,DATA2)
        DEALLOCATE(DATA,DATA2)
      ENDIF
      IF(IMPX.GT.3)CALL LCMLIB(KPNEW)
  100 CONTINUE
      DEALLOCATE(SCAT,SCAT2)
      RETURN
      END
