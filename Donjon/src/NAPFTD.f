*DECK NAPFTD
      SUBROUTINE NAPFTD(NXP,MXP,NXD,MXD,FXTD)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform a projection of second geometry on first one to compute 
*     fraction of region of the first geometry occupied by the second 
*     geometry regions
*
*Copyright:
* Copyright (C) 2014 Ecole Polytechnique de Montreal.
*
*Author(s): 
* R. Chambon
*
*Parameters: input/output
*    for core with heterogeneous mixture
* NXP    number of region along X direction for first geometry
* MXP    mesh of region along X direction for first geometry
* NXD    number of region along X direction for second geometry
* MXD    mesh of region along X direction for second geometry
* FXTD   fraction of region along X direction
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER      NXP,NXD
      REAL         MXP(NXP),MXD(NXD),FXTD(NXP,NXD)
*----
*  LOCAL VARIABLES
*----
      INTEGER IP,ID
      REAL DXP

      DO IP=1,NXP
      DXP=MXP(IP+1)-MXP(IP)
      DO ID=1,NXD
        IF((MXD(ID).LE.MXP(IP)).AND.(MXD(ID+1).GE.MXP(IP+1))) THEN
          FXTD(IP,ID)=1.0
        ELSEIF ((MXD(ID).LE.MXP(IP)).AND.(MXD(ID+1).GT.MXP(IP))) THEN
          FXTD(IP,ID)=(MXD(ID+1)-MXP(IP))/DXP
        ELSEIF ((MXD(ID).GE.MXP(IP)).AND.
     1          (MXD(ID+1).LE.MXP(IP+1))) THEN
          FXTD(IP,ID)=(MXD(ID+1)-MXD(ID))/DXP
        ELSEIF ((MXD(ID).LT.MXP(IP+1)).AND.
     1          (MXD(ID+1).GE.MXP(IP+1))) THEN
          FXTD(IP,ID)=(MXP(IP+1)-MXD(ID))/DXP
        ENDIF
      ENDDO
      ENDDO

      RETURN
      END
