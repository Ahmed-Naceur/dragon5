*DECK KINPOW
      SUBROUTINE KINPOW(IPMAC,NGR,NBM,NUN,NEL,MAT,VOL,IDL,EVECT,POWTOT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute reactor power.
*
*Copyright:
* Copyright (C) 2011 Ecole Polytechnique de Montreal.
*
*Author(s): R. Chambon
*
*Parameters: input
* IPMAC   addresses of L_MACROLIB object.
* NGR     number of energy groups.
* NBM     number of material mixtures.
* NUN     total number of unknowns per energy group.
* NEL     total number of finite elements.
* MAT     mixture index assigned to each finite element.
* VOL     volume of each element.
* IDL     position of the average flux component associated with each
*         finite element.
* EVECT   neutron flux.
* IMPX    print parameter (equal to zero for no print).
*
*Parameters: output
* POWTOT  power.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NGR,NBM,NUN,NEL,MAT(NEL),IDL(NEL)
      TYPE(C_PTR) IPMAC
      REAL EVECT(NUN,NGR),VOL(NEL),POWTOT
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION POWD
      INTEGER IGR,IEL,ITYLCM,LENGT
      TYPE(C_PTR) JPMAC,KPMAC
      REAL, DIMENSION(:), ALLOCATABLE :: HF
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(HF(NBM))
*
      POWTOT=0.0
      JPMAC=LCMGID(IPMAC,'GROUP')
      KPMAC=LCMGIL(JPMAC,1)
      CALL LCMLEN(KPMAC,'H-FACTOR',LENGT,ITYLCM)
      IF(LENGT.EQ.0) RETURN
      IF(LENGT.NE.NBM) CALL XABORT('@KINPOW: INVALID LENGTH FO'
     1  //'R H-FACTOR INFORMATION.')
*----
*  Compute power as H*Phi*Vol.
*----
      CALL XDRSET(HF,NBM,0.0)
      POWD=0.0D0
      DO 20 IGR=1,NGR
        KPMAC=LCMGIL(JPMAC,IGR)
        CALL LCMGET(KPMAC,'H-FACTOR',HF)
        DO 10 IEL=1,NEL
          IF(MAT(IEL).GT.0) THEN
            POWD=POWD+VOL(IEL)*HF(MAT(IEL))*EVECT(IDL(IEL),IGR)
          ENDIF
   10   CONTINUE
   20 CONTINUE
      POWTOT=REAL(POWD)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(HF)
      RETURN
      END
