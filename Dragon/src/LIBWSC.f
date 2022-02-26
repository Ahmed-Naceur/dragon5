*DECK LIBWSC
      SUBROUTINE LIBWSC(NGROUP,NGD,NGF,NSCT,CSCAT,XSCAT,SIGS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Expand WIMS format scattering cross sections.
*
*Copyright:
* Copyright (C) 1997 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): 
* G. Marleau
*
*Parameters: input
* NGROUP  number of groups.                         
* NGD     starting group number.                    
* NGF     finishing group number.                   
* NSCT    number of elements in CSCAT.              
* CSCAT   WIMS condense scattering at input.        
*
*Parameters: output
* XSCAT   DRAGON format expanded scattering.        
*         SCAT(JG,IG) is from IG to JG.              
* SIGS    total scattering out of group.            
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
* INTERFACE VARIABLES
*----
      INTEGER  NGROUP,NGD,NGF,NSCT
      REAL     CSCAT(NSCT),XSCAT(NGROUP,NGROUP),SIGS(NGROUP)
*----
* LOCAL VARIABLES
*----
      INTEGER          LG,IG1,IG2,N2,IGG
      DOUBLE PRECISION SUMSCT
*
*----
*
      LG=0
      DO 100 IG1=NGD,NGF
        CALL XDRSET(XSCAT(1,IG1),NGROUP,0.0)
        LG=LG+1
        IG2=IG1-INT(CSCAT(LG)+0.1)
        LG=LG+1
        N2=INT(CSCAT(LG)+0.1)
        SUMSCT=0.0D0
        DO 110 IGG=1,N2
          LG=LG+1
          IG2=IG2+1
          IF(IG2.LT.1) THEN
            CALL XABORT('LIBWSC: IG2 < 1')
          ELSE IF(IG2.GT.NGROUP) THEN
            CALL XABORT('LIBWSC: IG2 > NGROUP')
          ENDIF
          XSCAT(IG2,IG1)=CSCAT(LG)
          SUMSCT=SUMSCT+DBLE(CSCAT(LG))
 110    CONTINUE
        SIGS(IG1)=REAL(SUMSCT)
 100  CONTINUE
      IF(LG.NE.NSCT) CALL XABORT('LIBWSC: INVALID COUNT')
*----
*  RETURN LIBWSC
*----
      RETURN
      END
