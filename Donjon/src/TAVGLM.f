*DECK TAVGLM
      SUBROUTINE TAVGLM(NB,SHIFT,BCHAN,PSI,BURN0,BURN1,IVECT,NSCH)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the burnup integration limits for a given channel.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
*
*Author(s): 
* D.Rozon, M.Beaudet, D.Sekki
*
*Parameters: input
* NB     number of fuel bundles.
* SHIFT  number of bundles to refuel (bundle-shift).
* PSI    axial shape over each bundle.
* NSCH   refuelling scheme of a given channel.
* BCHAN  average exit burnup for a given channel.
* IVECT  refuelling pattern vector for a given channel.
*
*Parameters: output
* BURN0  lower burnup integration limit.
* BURN1  upper burnup integration limit.
*
*Parameters: scratch
* DELT   incremental burnup over each fuel bundle.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NB,SHIFT,NSCH,IVECT(NB)
      REAL BURN0(NB),BURN1(NB),PSI(NB),BCHAN
*----
*  LOCAL VARIABLES
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: DELT
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(DELT(NB))
*----
*  COMPUTE BURNUP LIMITS
*----
      CALL XDRSET(BURN0,NB,0.)
      CALL XDRSET(BURN1,NB,0.)
      CALL XDRSET(DELT,NB,0.)
      DO 10 IB=1,NB
      DELT(IB)=SHIFT*BCHAN*PSI(IB)
   10 CONTINUE
*     NEGATIVE DIRECTION
      IF(NSCH.LT.0)THEN
        DO 20 IB=1,NB
        KK=NB-IB+1
        KA=NB-IVECT(IB)+1
        IF(IVECT(IB).LE.0)THEN
          BURN0(KK)=0.
        ELSE
          BURN0(KK)=BURN1(KA)
        ENDIF
        BURN1(KK)=BURN0(KK)+DELT(KK)
   20   CONTINUE
*     POSITIVE DIRECTION
      ELSE
        DO 30 IB=1,NB
        IF(IVECT(IB).LE.0)THEN
          BURN0(IB)=0.
        ELSE
          BURN0(IB)=BURN1(IVECT(IB))
        ENDIF
        BURN1(IB)=BURN0(IB)+DELT(IB)
   30   CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(DELT)
      RETURN
      END
