*DECK RESROD
      SUBROUTINE RESROD(NB,NZ,ZZ,IND,ZLEVEL,ITOP,VB)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Move a control rod over a fuel channel.
*
*Copyright:
* Copyright (C) 2017 Ecole Polytechnique de Montreal.
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* NB      number of fuel bundles per channel.
* NZ      number of axial meshes.
* ZZ      axial meshes.
* IND     bundle index of each axial mesh.
* ZLEVEL  insertion parameter of the control rod in the channel (set
*         between 0.0 and 1.0).
* ITOP    direction flag for the rod (=1: from top; =-1: from bottom).
*
*Parameters: output
* VB      insertion parameter corresponding to each bundle.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NB,NZ,IND(NZ),ITOP
      REAL ZZ(NZ+1),ZLEVEL,VB(NB)
*
      ZMIN=ZZ(NZ+1)
      ZMAX=ZZ(1)
      DO IZ=1,NZ
        IF(IND(IZ).EQ.0) CYCLE
        ZMIN=MIN(ZZ(IZ),ZMIN)
        ZMAX=MAX(ZZ(IZ+1),ZMAX)
      ENDDO
      IF(ITOP.EQ.1) THEN
        CALL XDRSET(VB,NB,0.0)
        ZPOS=ZMAX-ZLEVEL*(ZMAX-ZMIN)
        DO IB=1,NB
          ZBMIN=ZZ(NZ+1)
          ZBMAX=ZZ(1)
          DO IZ=1,NZ
            IF(IND(IZ).EQ.IB) THEN
              ZBMIN=MIN(ZZ(IZ),ZBMIN)
              ZBMAX=MAX(ZZ(IZ+1),ZBMAX)
            ENDIF
          ENDDO
          IF((ZPOS.GE.ZBMIN).AND.(ZPOS.LE.ZBMAX)) THEN
            VB(IB)=1.0-(ZPOS-ZBMIN)/(ZBMAX-ZBMIN)
            CALL XDRSET(VB(IB+1),NB-IB,1.0)
            EXIT
          ENDIF
        ENDDO
      ELSEIF(ITOP.EQ.-1) THEN
        CALL XDRSET(VB,NB,1.0)
        ZPOS=ZMIN+ZLEVEL*(ZMAX-ZMIN)
        DO IB=1,NB
          ZBMIN=ZZ(NZ+1)
          ZBMAX=ZZ(1)
          DO IZ=1,NZ
            IF(IND(IZ).EQ.IB) THEN
              ZBMIN=MIN(ZZ(IZ),ZBMIN)
              ZBMAX=MAX(ZZ(IZ+1),ZBMAX)
            ENDIF
          ENDDO
          IF((ZPOS.GE.ZBMIN).AND.(ZPOS.LE.ZBMAX)) THEN
            VB(IB)=(ZPOS-ZBMIN)/(ZBMAX-ZBMIN)
            CALL XDRSET(VB(IB+1),NB-IB,0.0)
            EXIT
          ENDIF
        ENDDO
      ENDIF
      RETURN
      END
