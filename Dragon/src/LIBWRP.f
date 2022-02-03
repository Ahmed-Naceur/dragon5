*DECK LIBWRP
      SUBROUTINE LIBWRP(IPRINT,NTYP,NGR,NRTOT,MAXTEM,MAXDIL,IGR,IRES,
     >                  ITYP,DSIGPL,NTM,NDI,RTMP,RDIL,RESI,NTMPR,NDILR,
     >                  TMPT,DILT,REST)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Prepare WIMS-D4 resonance data.
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
* IPRINT  print flag.           
* NTYP    number of resonance tables per isotopes. 
* NGR     number of resonance groups.           
* NRTOT   maminum number of resonant isotopes.  
* MAXTEM  maminum number of temperature.        
* MAXDIL  maminum number of dilutions.          
* IGR     resonance group number.               
* IRES    resonance isotope set.                
* ITYP    XS type.                              
* DSIGPL  background XS.                        
* NTM     number of temperatures.               
* NDI     number of dilutions.                  
* RTMP    resonance temperature.                
* RDIL    resonance dilution.                   
* RESI    resonance integrals.                  
*
*Parameters: output
* NTMPR   number of local temperatures.         
* NDILR   number of local dilutions.            
* TMPT    work temperature.                     
* DILT    work dilution.                       
* REST    work resonance integrals.             
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
* PARAMETERS
*----
      INTEGER   IOUT
      CHARACTER NAMSBR*6
      PARAMETER (IOUT=6,NAMSBR='LIBWRP')
*----
* INTERFACE VARIABLES
*----
      INTEGER   IPRINT,NTYP,NGR,NRTOT,MAXTEM,MAXDIL,IGR,IRES,ITYP,
     1          NTMPR,NDILR
      INTEGER   NTM(NTYP,NRTOT,NGR),NDI(NTYP,NRTOT,NGR)
      REAL      DSIGPL
      REAL      RTMP(MAXTEM,NTYP,NRTOT,NGR),RDIL(MAXDIL,NTYP,NRTOT,NGR),
     1          RESI(MAXDIL,MAXTEM,NTYP,NRTOT,NGR),TMPT(MAXTEM),
     2          DILT(MAXDIL),REST(MAXDIL*MAXTEM)
*----
* LOCAL VARIABLES
*----
      INTEGER   ITT,IT,IPOS,ID
      REAL      XDIL
*
*----
*
      NTMPR=NTM(ITYP,IRES,IGR)
      NDILR=NDI(ITYP,IRES,IGR)
      IF(ABS(IPRINT) .GE. 100) THEN
        WRITE(IOUT,6010) NAMSBR
        WRITE(IOUT,6000)
        WRITE(IOUT,6002) (RTMP(ITT,ITYP,IRES,IGR),ITT=1,NTMPR)
        WRITE(IOUT,6001)
        WRITE(IOUT,6002) (RDIL(ITT,ITYP,IRES,IGR),ITT=1,NDILR)
      ENDIF
      DO 100 IT=1,NTMPR
        TMPT(IT)=SQRT(RTMP(IT,ITYP,IRES,IGR))
 100  CONTINUE
      DO 110 ID=1,NDILR
        XDIL=RDIL(ID,ITYP,IRES,IGR)-DSIGPL
        IF(XDIL.GT.0.0) THEN
          DILT(ID)=SQRT(XDIL)
        ELSE
          DILT(ID)=0.0
        ENDIF
 110  CONTINUE
      IPOS=0
      DO 120 IT=1,NTMPR
        DO 121 ID=1,NDILR
          IPOS=IPOS+1
          REST(IPOS)=RESI(ID,IT,ITYP,IRES,IGR)
 121    CONTINUE
 120  CONTINUE
      IF(ABS(IPRINT) .GE. 100) THEN
        WRITE(IOUT,6011) NAMSBR
      ENDIF
      RETURN
*----
*  FORMAT
*----
 6000 FORMAT('   RESONANCE TEMPERATURE TABULATION = ')
 6001 FORMAT('   RESONANCE DILUTIONS TABULATION   = ')
 6002 FORMAT(1P,5E15.7)
 6010 FORMAT('(* Output from --',A6,'-- follows ')
 6011 FORMAT('   Output from --',A6,'-- completed *)')
      END
