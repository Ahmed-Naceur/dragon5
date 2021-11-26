*DECK USPMIX
      SUBROUTINE USPMIX(IPMTX,NEL,NREFL,NFUEL,MAT,RMIX,FMIX,INDX,NMIX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover and check the material mixtures.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* D. Sekki
*
*Parameters: input
* IPMTX  pointer to matex information.
* NEL    total number of volumes in reactor geometry.
* NREFL  total number of reflector types.
* NFUEL  total number of fuel types.
* MAT    material index from geometry.
* RMIX   reflector-type mixtures indices.
* FMIX   fuel-type mixtures indices.
*
*Parameters: output
* INDX   renumbered material index.
* NMIX   total number of non-virtual volumes.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMTX
      INTEGER NEL,NREFL,NFUEL,MAT(NEL),RMIX(NREFL),FMIX(NFUEL),
     1 INDX(NEL),NMIX
*----
*  LOCAL VARIABLES
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: RTOT,FTOT
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(RTOT(NREFL),FTOT(NFUEL))
*----
*  MATERIAL INDEX
*----
      CALL XDISET(RTOT,NREFL,0)
      CALL XDISET(FTOT,NFUEL,0)
      CALL XDISET(INDX,NEL,0)
      NMIX=0
      DO IEL=1,NEL
        IF(MAT(IEL).NE.0)THEN
          NMIX=NMIX+1
          INDX(IEL)=NMIX
        ENDIF
      ENDDO
      IF((NFUEL.EQ.0).AND.(NREFL.EQ.0)) GOTO 20
*     CHECK MIXTURES
      DO 10 IEL=1,NEL
      IMIX=MAT(IEL)
      IF(IMIX.EQ.0)GOTO 10
      IF(NREFL.EQ.0)GOTO 5
      DO IREFL=1,NREFL
        IF(IMIX.EQ.RMIX(IREFL))THEN
          RTOT(IREFL)=RTOT(IREFL)+1
          GOTO 10
        ENDIF
      ENDDO
    5 IF(NFUEL.EQ.0)GOTO 10
      DO IFUEL=1,NFUEL
        IF(IMIX.EQ.FMIX(IFUEL))THEN
          FTOT(IFUEL)=FTOT(IFUEL)+1
          GOTO 10
        ENDIF
      ENDDO
   10 CONTINUE
*     STORAGE
20    CALL LCMPUT(IPMTX,'MAT',NEL,1,MAT)
      CALL LCMPUT(IPMTX,'INDEX',NEL,1,INDX)
      IF(NREFL.NE.0) THEN
      	CALL LCMPUT(IPMTX,'RMIX',NREFL,1,RMIX)
      	CALL LCMPUT(IPMTX,'RTOT',NREFL,1,RTOT)
      ENDIF
      IF(NFUEL.NE.0) THEN
      	CALL LCMPUT(IPMTX,'FMIX',NFUEL,1,FMIX)
      	CALL LCMPUT(IPMTX,'FTOT',NFUEL,1,FTOT)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(FTOT,RTOT)
      RETURN
      END
