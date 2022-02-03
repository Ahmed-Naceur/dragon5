*DECK DOORFB3
      SUBROUTINE DOORFB3(IPSYS,IPTRK,IMPX,NBMIX,NREG,NUN,KEYFLX,SUNKNO,
     1 FUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Double heterogeneity treatment (part 3). One-group version.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPSYS   pointer to the assembly LCM object (L_PIJ signature). IPSYS is
*         a list of directories.
* IPTRK   pointer to the tracking (L_TRACK signature).
* IMPX    print flag (equal to zero for no print).
* NBMIX   number of composite mixtures in the domain. Equal to the
*         number of mixtures in the internal library.
* NREG    number of volumes in the composite geometry.
* NUN     total number of unknowns in vector SUNKNO.
* KEYFLX  index of flux components in unknown vector.
* SUNKNO  equivalent macro-sources.
*
*Parameters: input/output
* FUNKNO  equivalent macro-fluxes on input and 
*         composite fluxes on output.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPSYS,IPTRK
      INTEGER IMPX,NBMIX,NREG,NUN,KEYFLX(NREG)
      REAL SUNKNO(NUN)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IOUT=6)
      INTEGER IPAR(8)
*----
*  ALLOCATABLE ARRAYS
*----
      TYPE(C_PTR) NS_PTR,IBI_PTR,FRACT_PTR,VOLK_PTR,IDIL_PTR,MIXGR_PTR,
     > NCO_PTR,RRRR_PTR,QKOLD_PTR,QKDEL_PTR,PKL_PTR,COEF_PTR,P1I_PTR,
     > P1DI_PTR,P1KI_PTR,SIGA1_PTR
      INTEGER, DIMENSION(:), POINTER :: NS,IBI,IDIL,MIXGR,NCO
      REAL, DIMENSION(:), POINTER :: FRACT,VOLK,RRRR,QKOLD,QKDEL,PKL,
     > P1I,P1DI,P1KI,SIGA1
      REAL, DIMENSION(:), ALLOCATABLE :: SGAR,SGAS
      DOUBLE PRECISION, DIMENSION(:), POINTER :: COEF
*----
*  RECOVER DOUBLE HETEROGENEITY DATA
*----
      IF(IMPX.GT.50) THEN
         WRITE(6,'(/38H DOORFB3: DOUBLE HETEROGENEITY OPTION.)')
      ENDIF
      CALL LCMSIX(IPTRK,'BIHET',1)
      CALL LCMGET(IPTRK,'PARAM',IPAR)
      IR1=IPAR(1)
      IR2=IPAR(2)
      NREG2=IPAR(3)
      NG=IPAR(4)
      NSMAX=IPAR(5)
      IBIHET=IPAR(6)
      IF(IR1.NE.NBMIX) CALL XABORT('DOORFB3: INVALID DATA IN TRACKING.')
      IF(NREG2.GT.NUN) CALL XABORT('DOORFB3: NUN OVERFLOW.')
      NMILG=IR2-IR1
      CALL LCMGPD(IPTRK,'NS',NS_PTR)
      CALL LCMGPD(IPTRK,'IBI',IBI_PTR)
      CALL LCMGPD(IPTRK,'FRACT',FRACT_PTR)
      CALL LCMGPD(IPTRK,'VOLK',VOLK_PTR)
      CALL LCMGPD(IPTRK,'IDIL',IDIL_PTR)
      CALL LCMGPD(IPTRK,'MIXGR',MIXGR_PTR)
      CALL LCMSIX(IPTRK,' ',2)
*
      CALL C_F_POINTER(NS_PTR,NS,(/ NG /))
      CALL C_F_POINTER(IBI_PTR,IBI,(/ NREG2 /))
      CALL C_F_POINTER(FRACT_PTR,FRACT,(/ NG*(NBMIX+NMILG) /))
      CALL C_F_POINTER(VOLK_PTR,VOLK,(/ NG*NSMAX /))
      CALL C_F_POINTER(IDIL_PTR,IDIL,(/ NMILG /))
      CALL C_F_POINTER(MIXGR_PTR,MIXGR,(/ NSMAX*NG*NMILG /))
*----
*   RECOVER GROUP-DEPENDENT BIHET INFORMATION
*----
      IF((IBIHET.EQ.1) .OR. (IBIHET.EQ.2)) THEN
        CALL LCMGPD(IPSYS,'NCO',NCO_PTR)
        CALL LCMGPD(IPSYS,'RRRR',RRRR_PTR)
        CALL LCMGPD(IPSYS,'QKOLD',QKOLD_PTR)
        CALL LCMGPD(IPSYS,'QKDEL',QKDEL_PTR)
        CALL LCMGPD(IPSYS,'PKL',PKL_PTR)
        CALL LCMGPD(IPSYS,'COEF',COEF_PTR)
*  
        CALL C_F_POINTER(NCO_PTR,NCO,(/ NMILG /))
        CALL C_F_POINTER(RRRR_PTR,RRRR,(/ NMILG /))
        CALL C_F_POINTER(QKOLD_PTR,QKOLD,(/ NG*NSMAX*NMILG /))
        CALL C_F_POINTER(QKDEL_PTR,QKDEL,(/ NG*NSMAX*NMILG /))
        CALL C_F_POINTER(PKL_PTR,PKL,(/ NG*NSMAX*NSMAX*NMILG /))
        CALL C_F_POINTER(COEF_PTR,COEF,(/ NMILG*(1+NG*NSMAX)**2 /))
      ELSE IF((IBIHET.EQ.3) .OR. (IBIHET.EQ.4)) THEN
        CALL LCMGPD(IPSYS,'P1I',P1I_PTR)
        CALL LCMGPD(IPSYS,'P1DI',P1DI_PTR)
        CALL LCMGPD(IPSYS,'P1KI',P1KI_PTR)
        CALL LCMGPD(IPSYS,'SIGA1',SIGA1_PTR)
*  
        CALL C_F_POINTER(P1I_PTR,P1I,(/ NG*NMILG /))
        CALL C_F_POINTER(P1DI_PTR,P1DI,(/ NG*NMILG /))
        CALL C_F_POINTER(P1KI_PTR,P1KI,(/ NSMAX*NG*NMILG /))
        CALL C_F_POINTER(SIGA1_PTR,SIGA1,(/ NG*NMILG /))
      ENDIF
*----
*  COMPUTE THE EQUIVALENT CROSS SECTIONS IN COMPOSITE REGIONS
*----
      CALL LCMSIX(IPSYS,'BIHET',1)
      NB1=NBMIX+1
      CALL LCMLEN(IPSYS,'DRAGON-S0XSC',ILONG,ITYLCM)
      NANI=ILONG/(NB1+NMILG)
      ALLOCATE(SGAR(NB1+NMILG),SGAS((NB1+NMILG)*NANI))
      CALL XDRSET(SGAS,(NB1+NMILG)*NANI,0.0)
      CALL LCMGET(IPSYS,'DRAGON-TXSC',SGAR)
      CALL LCMGET(IPSYS,'DRAGON-S0XSC',SGAS)
      CALL LCMSIX(IPSYS,' ',2)
*----
*   DOUBLE HETEROGENEITY TREATMENT -- PART 2
*----
      IF((IBIHET.EQ.1) .OR. (IBIHET.EQ.2)) THEN
      CALL XDRH30(IBIHET,NUN,NBMIX,NMILG,NREG,NREG2,NG,NSMAX,KEYFLX,
     > NS,IDIL,MIXGR,IBI,FRACT,VOLK,SGAR,SGAS,NCO,RRRR,QKOLD,QKDEL,
     > PKL,COEF,SUNKNO,FUNKNO)
      ELSE IF((IBIHET.EQ.3) .OR. (IBIHET.EQ.4)) THEN
      CALL XDRH33(IBIHET,NUN,NBMIX,NMILG,NREG,NREG2,NG,NSMAX,KEYFLX,
     > NS,IDIL,MIXGR,IBI,FRACT,VOLK,SGAR,P1I,P1DI,P1KI,SIGA1,FUNKNO)
      ENDIF
*----
*   MEMORY RELEASE
*----
      DEALLOCATE(SGAS,SGAR)
      RETURN
      END
