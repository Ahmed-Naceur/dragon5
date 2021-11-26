*DECK RANDDN
      SUBROUTINE RANDDN(ISEED,NRAND,DRANDN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* This subroutine generates pseudo-random numbers
* from a normal distribution of width 1 centered at 0.0.
* 
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal
*
*Author(s): G. Marleau
*
*Parameters: input
* ISEED   the seed for the generation of random numbers. 
* NRAND   number of random number requested.
*                 
*Parameters: ouput
* DRANDN  random numbers between picked from 
*         a normal distribution of width 1 centered at 0 .
*
*Reference:                     
* Box-Muller method.
*                          
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER          ISEED,NRAND
      DOUBLE PRECISION DRANDN(NRAND)
*----
*  Parameters
*----
      DOUBLE PRECISION TWOPI
      PARAMETER       (TWOPI=6.283185307179586D0)
*----
*  LOCAL VARIABLES
*----
      INTEGER          IRAND,NSTEP,ISTEP
      DOUBLE PRECISION DXRAND(2),DSRAND
*----
*  Saved variables
*----
      INTEGER          INEXT
      DOUBLE PRECISION DLAST
      SAVE             INEXT,DLAST
      DATA             INEXT/0/
      DATA             DLAST/0.0D0/
*----
*  Pick 2 random numbers
*----
      NSTEP=(NRAND-INEXT)/2
      IF(INEXT .EQ. 1) THEN
        DRANDN(INEXT)=DLAST
      ENDIF
      IRAND=INEXT+1
      DO ISTEP=1,NSTEP
        CALL RANDD(ISEED,2,DXRAND)
        DSRAND=SQRT(-2*LOG(DXRAND(1)))
        DRANDN(IRAND)=DSRAND*COS(TWOPI*DXRAND(2))
        DRANDN(IRAND+1)=DSRAND*SIN(TWOPI*DXRAND(2))
        IRAND=IRAND+2
      ENDDO
      IF(MOD(NRAND-INEXT,2) .EQ. 1) THEN
        CALL RANDD(ISEED,2,DXRAND)
        DSRAND=SQRT(-2*LOG(DXRAND(1)))
        DRANDN(IRAND)=DSRAND*COS(TWOPI*DXRAND(2))
        DLAST=DSRAND*SIN(TWOPI*DXRAND(2))
        INEXT=1
      ELSE
        INEXT=0
      ENDIF
      RETURN
      END
