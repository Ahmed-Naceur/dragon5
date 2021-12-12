*DECK XL3SIG
      SUBROUTINE XL3SIG(  NGRT,  NBMIX,  XSSIGT,ALBEDO,  NPSYS,
     >                    NGRP,     NS,     NR, MATALB,    VOL,
     >                  SIGTAL, SIGVOL, SWVOID, SWNZBC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Unfold cross-section data which becomes available by subset of groups.
*
*Copyright:
* Copyright (C) 1996 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* NGRT    total number of groups.                       
* NBMIX   number of mixtures in the MACROLIB.                
* XSSIGT  total XS for mixtures in the MACROLIB.        
* ALBEDO  geometric albedos.                           
* NPSYS   group masks.                                  
* NGRP    number of groups.                             
* NS      number of surfaces in the assembly.               
* NR      number of zones in the assembly.                  
* MATALB  material numbers for zones in the supercell.       
* VOL     volumes.                                      
*
*Parameters: output
* SIGTAL  total XS and albedos by region & surface.       
* SIGVOL  volume times total XS by region.              
* SWVOID  logical switch (.TRUE. if void regions).     
* SWNZBC  logical switch (.TRUE. if non-zero BC).     
*
*-----------------------------------------------------------------------
*
      IMPLICIT   NONE
*
      INTEGER    NGRT,NBMIX,NGRP,NS,NR,NPSYS(NGRP)
      REAL       XSSIGT(0:NBMIX,NGRT),VOL(NR),SIGTAL(NS:NR,NGRP),
     >           SIGVOL(NR,NGRP),ALBEDO(6)
      INTEGER    MATALB(NS:NR)
      INTEGER    IUN,JG,I
      LOGICAL    SWVOID,SWNZBC
      REAL       ZERO
      INTEGER    IOUT
      PARAMETER (ZERO=0.0,IOUT=6)
*
      SWVOID= .FALSE.
      SWNZBC= .FALSE.
*
      IF( NS.GT.0 ) CALL XABORT('XL3SIG: # OF SURFACES IS > 0')
      IF( NR.LT.0 ) CALL XABORT('XL3SIG: # OF REGIONS  IS < 0')
*
      DO 10 IUN= NS, NR
      DO 20 JG= 1, NGRP
         IF(NPSYS(JG).EQ.0) GO TO 20
         IF( IUN.LT.0 )THEN
            SIGTAL(IUN,JG)= ALBEDO(-MATALB(IUN))
            SWNZBC=SWNZBC.OR.(ALBEDO(-MATALB(IUN)).NE.ZERO)
         ELSEIF( IUN.EQ.0 )THEN
            SIGTAL(IUN,JG)= ZERO
         ELSE
            IF( MATALB(IUN).LT.0.OR.MATALB(IUN).GT.NBMIX)THEN
               WRITE(IOUT,*) 'NBMIX=',NBMIX
               WRITE(IOUT,*) 'IG/NGRT=',JG,NGRT
               WRITE(IOUT,*) 'XSSIGT=',(XSSIGT(I,JG),I=0,NBMIX)
               WRITE(IOUT,*) 'MATALB<=0 =',(MATALB(I),I=NS,0)
               WRITE(IOUT,*) 'MATALB >0 =',(MATALB(I),I=1,NR)
               CALL XABORT('XL3SIG: INVALID NUMBER OF MIXTURES')
            ENDIF
            SIGTAL(IUN,JG)= XSSIGT(MATALB(IUN),JG)
            SIGVOL(IUN,JG)= SIGTAL(IUN,JG)*VOL(IUN)
         ENDIF
   20 CONTINUE
   10 CONTINUE
*
      RETURN
      END
