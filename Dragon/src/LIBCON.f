*DECK LIBCON
      SUBROUTINE LIBCON(IPLIB,IMX,NBISO,ISOMIX,DENISO,DENMIX,IN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Convert weight percent to atomic density.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau
*
*Parameters: input
* IPLIB   pointer to the internal library.
* IMX     mixture index to process.
* NBISO   number of isotopes present in the calculation domain.
* ISOMIX  mix number of each isotope.
* IN      type of conversion:
*         =1 conversion of wgt% to nb atoms with denmix;
*         =2 conversion of nb atoms to wgt% and denmix.
*
*Parameters: input/output
* DENISO  number density (if IN=1) or weight percent (if IN=2) for
*         isotopes present in mixture IMX on input. On optput, 
*         number density.
* DENMIX  mixture density g*cm**(-3) (if IN=2).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB
      INTEGER IMX,NBISO,ISOMIX(NBISO),IN
      REAL DENISO(NBISO),DENMIX
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) KPLIB
      DOUBLE PRECISION XDRCST,AVCON
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: IPISO
*----
*  INTERNAL PARAMETERS
*----
      PARAMETER       (IOUT=6)
*------
*  COMPUTE NUMBER DENSITIES FOR ISOTOPES
*------
      AVCON=1.0D-24*XDRCST('Avogadro','N/moles')
     >             /XDRCST('Neutron mass','amu')
      IF(IN.EQ.1) THEN
        IF(DENMIX.EQ.-1.0) CALL XABORT('LIBCON: DENMIX NOT DEFINED')
         TWPC=0.0
         DO 120 ISO=1,NBISO
           IF(ISOMIX(ISO).EQ.IMX) TWPC=DENISO(ISO)+TWPC
 120     CONTINUE
         IF(TWPC.EQ.0.0) THEN
           IF(DENMIX.EQ.0.0) THEN
             RETURN
           ELSE
             CALL XABORT('LIBCON: A MIXTURE OF DENSITY > 0.0 '//
     >         'HAS ALL ITS ISOTOPIC WEIGHT PERCENT = 0.0')
           ENDIF
         ENDIF
         WMIX=DENMIX*REAL(AVCON)/TWPC
         IF(NBISO.GT.0) THEN
           ALLOCATE(IPISO(NBISO))
           CALL LIBIPS(IPLIB,NBISO,IPISO)
           DO 130 ISO=1,NBISO
             IF(ISOMIX(ISO).EQ.IMX) THEN
               KPLIB=IPISO(ISO) ! set ISO-th isotope
               CALL LCMGET(KPLIB,'AWR',AWRISO)
               IF(AWRISO.GT.0.0) THEN
                 DENISO(ISO)=DENISO(ISO)*WMIX/AWRISO
               ELSE
                 DENISO(ISO)=0.0
               ENDIF
             ENDIF
 130       CONTINUE
           DEALLOCATE(IPISO)
         ENDIF
      ELSE IF(IN.EQ.2) THEN
*------
*  COMPUTE MIXTURE DENSITIES AND ISOTOPIC WEIGHT PERCENTS
*  (NORMALIZED TO 100.)
*------
         DENMIX=0.0
         IF(NBISO.GT.0) THEN
           ALLOCATE(IPISO(NBISO))
           CALL LIBIPS(IPLIB,NBISO,IPISO)
           DO 220 ISO=1,NBISO
             IF(ISOMIX(ISO).EQ.IMX) THEN
               KPLIB=IPISO(ISO) ! set ISO-th isotope
               CALL LCMGET(KPLIB,'AWR',AWRISO)
               DENISO(ISO)=DENISO(ISO)*AWRISO/REAL(AVCON)
               DENMIX=DENMIX+DENISO(ISO)
             ENDIF
 220       CONTINUE
           DEALLOCATE(IPISO)
         ENDIF
         IF(DENMIX.NE.0.0) THEN
           DO 230 ISO=1,NBISO
             IF(ISOMIX(ISO).EQ.IMX)
     >         DENISO(ISO)=100.*DENISO(ISO)/DENMIX
 230       CONTINUE
         ENDIF
      ELSE
         CALL XABORT('LIBCON: INVALID *IN* VALUE')
      ENDIF
      RETURN
      END
