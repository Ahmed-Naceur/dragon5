*DECK LIBWTE
      SUBROUTINE LIBWTE(IACT,ITXS,NGROUP,NGTHER,NTMP,NF,TERP,SCAT,
     >                  SIGS,XSNG,SIGF,XSFI,TRAN,TMPXS,TMPSC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform temperature interpolation for WIMS-AECL or WIMS-D4 XS.
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
* IACT    Action:
*         = 1 initialize before adding;
*         = 2 only add.
* ITXS    type:
*         = 1 all cross sections;
*         = 2 only scattering.
* NGROUP  number of groups.
* NGTHER  number of thermal groups.
* NTMP    number of temperature.
* NF      flag for fissile.
* TERP    temperature coefficients.
*
*Parameters: input/output
* SCAT    complete scattering matrix           
*         SCAT(JG,IG) is from IG to JG.           
* SIGS    total scattering out of group.        
* XSNG    (n,g) XS.                               
* SIGF    nu*fission XS.
* XSFI    fission XS.                           
* TRAN    transport XS.                         
*
*Parameters: scratch
* TMPXS   temperature dependent vect XS.        
* TMPSC   temperature dependent scat XS.        
*
*Comments:
*   WIMS-AECL library parameters
*   MAXISO : max. nb. of iso = 246                
*   MLDEP  : maximum number of reaction per       
*            isotope = MAXISO +4
*   LPZ    : length of parameter array = 9   
*   LMASTB : length of mst tab = MAXISO+9         
*   LMASIN : length of mst idx = LMASTB-4         
*   LGENTB : length of gen tab = 6                
*   LGENIN : length of gen idx = LGENTB
*   MASTER : master index array                   
*   GENINX : general index array
*   NPZ    : list of main parameters              
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
* INTERFACE VARIABLES
*----
      INTEGER          IACT,ITXS,NGROUP,NGTHER,NTMP,NF
      DOUBLE PRECISION TERP(NTMP)
      REAL             SCAT(NGROUP,NGROUP),SIGS(NGROUP),
     1                 XSNG(NGROUP),SIGF(NGROUP),XSFI(NGROUP),
     2                 TRAN(NGROUP),TMPXS(NGROUP,5,NTMP),
     3                 TMPSC(NGROUP,NGROUP,NTMP)
*----
* LOCAL VARIABLES
*----
      INTEGER          IGF,ITM,IGD,NGD
      REAL             RTERP
*----
*  INITIALIZED IF REQUIRED
*----
      NGD=NGROUP-NGTHER+1
      IF(IACT.EQ.1) THEN
        IF(ITXS.EQ.1) THEN
          CALL XDRSET(XSNG(NGD),NGTHER,0.0)
          CALL XDRSET(TRAN(NGD),NGTHER,0.0)
          IF(NF.GT.1) THEN
            CALL XDRSET(SIGF(NGD),NGTHER,0.0)
            CALL XDRSET(XSFI(NGD),NGTHER,0.0)
          ENDIF
        ENDIF
        IF(ITXS.GE.1) THEN
          CALL XDRSET(SIGS(NGD),NGTHER,0.0)
          DO 110 IGD=NGD,NGROUP
            CALL XDRSET(SCAT(1,IGD),NGROUP,0.0)
 110      CONTINUE
        ENDIF
      ENDIF
*----
*  INTERPOLATE STANDARD CROSS SECTIONS IN TEMPERATURE
*----
      IF(ITXS.EQ.1) THEN
        DO 120 ITM=1,NTMP
          RTERP=REAL(TERP(ITM))
          IF(RTERP.NE.0.0) THEN
            DO 121 IGD=NGD,NGROUP
              TRAN(IGD)=TRAN(IGD)+RTERP*TMPXS(IGD,1,ITM)
              XSNG(IGD)=XSNG(IGD)+RTERP*TMPXS(IGD,2,ITM)
              IF(NF.GT.1) THEN
                SIGF(IGD)=SIGF(IGD)+RTERP*TMPXS(IGD,3,ITM)
                XSFI(IGD)=XSFI(IGD)+RTERP*TMPXS(IGD,4,ITM)
              ENDIF
 121        CONTINUE
          ENDIF
 120    CONTINUE
      ENDIF
*----
*  INTERPOLATE SCATTERING CROSS SECTIONS IN TEMPERATURE
*----
      IF(ITXS.GE.1) THEN
        DO 130 ITM=1,NTMP
          RTERP=REAL(TERP(ITM))
          IF(RTERP.NE.0.0D0) THEN
            DO 131 IGD=NGD,NGROUP
              SIGS(IGD)=SIGS(IGD)+RTERP*TMPXS(IGD,5,ITM)
              DO 132 IGF=1,NGROUP
                SCAT(IGF,IGD)=SCAT(IGF,IGD)+RTERP*TMPSC(IGF,IGD,ITM)
 132          CONTINUE
 131        CONTINUE
          ENDIF
 130    CONTINUE
      ENDIF
      RETURN
      END
