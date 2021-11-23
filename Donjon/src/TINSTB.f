*DECK TINSTB
      SUBROUTINE TINSTB(IPMAP,TIME,BURNSTP,NCH,NB,NF,BUNDPOW,BURNAVG,
     1 BURNINS,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute new burnup for each bundle given either an average burnup step
* or a burning time
*
*Copyright:
* Copyright (C) 2009 Ecole Polytechnique de Montreal
*
*Author(s): 
* B. Toueg
*
*Parameters: input/output
* IPMAP   pointer to fuel-map information.
* TIME    time to burn
* BURNSTP average burnup step
* NCH     number of reactor channels.
* NB      number of fuel bundles.
* NF      number of fuel types.
* BUNDPOW bundle powers.
* BURNAVG average burnup.
* BURNINS instantaneous burnups.
* IMPX    printing index (=0 for no print).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAP
      INTEGER NCH,NB,NF,IMPX
      REAL BUNDPOW(NCH,NB), BURNINS(NCH,NB)
      REAL TIME,BURNSTP, BURNAVG, PTOT, MASSTOT, WEIGHT
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40,IOUT=6)
      TYPE(C_PTR) JPMAP,KPMAP
      CHARACTER HSMG*131
      INTEGER, ALLOCATABLE, DIMENSION(:) :: FLMIX,IFLRANK
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: BUNDMIX
      REAL, ALLOCATABLE, DIMENSION(:) :: FLWEIGHT
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(BUNDMIX(NCH,NB),FLMIX(NF),FLWEIGHT(NF))
*----
*  RECOVER INFORMATION
*----
      CALL LCMLIB(IPMAP)
*     FUEL MIX
      CALL XDISET(BUNDMIX,NCH*NB,0)
      CALL LCMGET(IPMAP,'FLMIX',BUNDMIX)
*     BURN-INST
      CALL XDRSET(BURNINS,NCH*NB,0.)
      CALL LCMGET(IPMAP,'BURN-INST',BURNINS)
*     FUEL INFORMATION (WEIGHT & MIX)
      JPMAP=LCMGID(IPMAP,'FUEL')
      MAXFL=0 ! maximum fuel mix number
      DO IFL=1,NF
        KPMAP=LCMGIL(JPMAP,IFL)
        CALL LCMGET(KPMAP,'MIX',FLMIX(IFL))
        MAXFL=MAX(MAXFL,FLMIX(IFL))
        CALL LCMGET(KPMAP,'WEIGHT',FLWEIGHT(IFL))
      ENDDO
      IF(MAXFL.LT.NF)THEN
        WRITE(HSMG,'(38H@TINSTB: FOUND MAX FUEL MIX NUMBER : (,I6,
     1  8H) THOUGH,I7,23H FUEL MIXES ARE DEFINED)')
        CALL XABORT(HSMG)
      ENDIF
* the mix stored in FLMIX field of /FMAP/
* is not the rank of the fuel in FUEL Dir list of /FMAP/
      ALLOCATE(IFLRANK(MAXFL))
      CALL XDISET(IFLRANK,MAXFL,0)
      DO IFL=1,NF
        IFLRANK(FLMIX(IFL))=IFL
      ENDDO
*----
*  COMPUTE BURNAVG, PTOT, MASSTOT, ( TIME if BURNSTEP is specified)
*----
      BURNAVG=0.
      PTOT=0.
      MASSTOT=0.
      NTOT=0
      DO ICH=1,NCH
        DO IB=1,NB
          IBD=BUNDMIX(ICH,IB)
          IF(IBD.EQ.0) CYCLE
          IFL=IFLRANK(IBD)
          IF(IFL.EQ.0) CYCLE
          NTOT=NTOT+1
          WEIGHT = FLWEIGHT(IFL)
          BURNAVG=BURNAVG+BURNINS(ICH,IB)
          PTOT=PTOT+BUNDPOW(ICH,IB)
          MASSTOT=MASSTOT+WEIGHT
        ENDDO
      ENDDO
      BURNAVG=BURNAVG/REAL(NTOT)
      IF(TIME.EQ.0.)THEN
        TIME = BURNSTP*MASSTOT/PTOT
      ENDIF
      IF(IMPX.GT.0)THEN
        WRITE(IOUT,*)'@TINSTB: TOTAL POWER = ',PTOT,' kW'
        WRITE(IOUT,*)'@TINSTB: TOTAL FUEL MASS = ',MASSTOT,' kg'
        WRITE(IOUT,*)'@TINSTB: AVERAGE BURN UP BEFORE = ',
     1  BURNAVG,'MWd/t'
      ENDIF
*----
*  COMPUTE NEW BURN-INST GIVEN TIME
*----
      BURNAVG=0.
      NTOT=0
      DO ICH=1,NCH
        DO IB=1,NB
          IBD=BUNDMIX(ICH,IB)
          IF(IBD.EQ.0) CYCLE
          IFL=IFLRANK(IBD)
          IF(IFL.EQ.0) CYCLE
          NTOT=NTOT+1
          WEIGHT = FLWEIGHT(IFL)
          IF(WEIGHT.GT.0.)THEN
            BURNINS(ICH,IB)=BURNINS(ICH,IB)
     1      +(BUNDPOW(ICH,IB)/WEIGHT)*TIME
            BURNAVG=BURNAVG+BURNINS(ICH,IB)
          ELSE
            IF(IMPX.GT.0)THEN
              WRITE(IOUT,*)'@TINSTB: WARNING MIX ',
     1        BUNDMIX(ICH,IB),' WEIGHS ',WEIGHT,'kg'
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      BURNAVG=BURNAVG/REAL(NTOT)
      IF(IMPX.GT.0)THEN
        WRITE(IOUT,*)'@TINSTB: AVERAGE BURN UP AFTER = ',BURNAVG,'MWd/t'
      ENDIF
      CALL LCMPUT(IPMAP,'BURN-INST',NCH*NB,2,BURNINS)
*----
*  RELEASE MEMORY AND RETURN
*----
      DEALLOCATE(IFLRANK)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(FLWEIGHT,FLMIX,BUNDMIX)
      RETURN
      END
