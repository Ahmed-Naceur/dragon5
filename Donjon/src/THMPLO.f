*DECK THMPLO
      SUBROUTINE THMPLO(P,X,PHIL0)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the value of the corrective factor for two phase calculation
* of frictional pressure loss based on an homogeneous flow correlation
*
*Copyright:
* Copyright (C) 2014 Ecole Polytechnique de Montreal.
*
*Author(s): 
* P. Gallet
*
*Parameters: input
* P       pressure (Pa)
* X       steam quality   
*
*Parameters: output
* PHIL0   corrective factor for two phase pressure loss calculation
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      REAL P,X,PHIL0
*----
*  LOCAL VARIABLES
*----
      REAL TSAT,RHOL,RHOG,MUL,MUG,TG,TL,R1,R2,R3
*----
*  COMPUTE VALUE OF THE CORRECTIVE FACTOR USING DENSITIES AND
*  VISCOSITIES OF BOTH SATURATED WATER AND DRY SATURATED STEAM
*----
*     compute the values of the thermodynamic parameters of steam and
*     liquid phases using freesteam steam tables
      CALL THMSAT(P,TSAT)
      TG=TSAT+0.1
      TL=TSAT-0.1
      CALL THMPT(P,TL,RHOL,R1,R2,MUL,R3)
      CALL THMPT(P,TG,RHOG,R1,R2,MUG,R3)
      PHIL0=(1+X*(RHOL/RHOG-1))/((1+X*(MUL/MUG-1))**0.25)
      RETURN
      END