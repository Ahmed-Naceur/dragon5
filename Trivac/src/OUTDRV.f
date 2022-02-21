*DECK OUTDRV
      SUBROUTINE OUTDRV (IPGEOM,IPMAC1,IPFLUX,IPMAC2,MAXNEL,NBMIX,NL,
     1 NBFIS,NGRP,NEL,NUN,NALBP,HTRACK,IELEM,ICOL,MAT,VOL,IDL,TITR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Driver for the post-treatment of reactor calculation results.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A Hebert
*
*Parameters: input
* IPGEOM  L_GEOM pointer to the geometry.
* IPMAC1  L_MACROLIB pointer to the nuclear properties.
* IPFLUX  L_FLUX pointer to the solution.
* IPMAC2  L_MACROLIB pointer to the edition information.
* MAXNEL  maximum number of finite elements.
* NBMIX   number of material mixtures.
* NL      scattering anisotropy.
* NBFIS   number of fissionable isotopes.
* NGRP    total number of energy groups.
* NEL     total number of finite elements.
* NUN     total number of unknowns per group.
* NALBP   number of physical albedos.
* HTRACK  type of tracking (equal to 'BIVAC' or 'TRIVAC').
* IELEM   degree of the Lagrangian finite elements:
* ICOL    type of quadrature used to integrate the mass matrix
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* IDL     position of the average flux component associated with
*         each volume.
* TITR    title.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPGEOM,IPMAC1,IPMAC2,IPFLUX
      CHARACTER TITR*72,HTRACK*12
      INTEGER MAXNEL,NBMIX,NL,NBFIS,NGRP,NEL,NUN,NALBP,IELEM,ICOL,
     1 MAT(NEL),IDL(NEL)
      REAL VOL(NEL)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPMAC1,KPMAC1
      CHARACTER TEXT4*4
      REAL NORM
      DOUBLE PRECISION DFLOTT,ZNORM
      INTEGER, DIMENSION(:), ALLOCATABLE :: IHOM,IGCOND,MATCOD
      REAL, DIMENSION(:), ALLOCATABLE :: SGD,FLUXC
      REAL, DIMENSION(:,:), ALLOCATABLE :: EVECT,ADECT,ZUFIS
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IHOM(NEL),IGCOND(NGRP),EVECT(NUN,NGRP),SGD(NBMIX),
     1 FLUXC(NEL),MATCOD(NEL))
*
      TKR=0.0
      IMPX=1
      IADJ=0
      NGCOND=NGRP
      DO IGR=1,NGRP
        IGCOND(IGR)=IGR
      ENDDO
      LMOD=0
      CALL KDRCPU(TK1)
*----
*  RECOVER THE K-EFFECTIVE AND THE DIRECT FLUX.
*----
      CALL LCMLEN(IPFLUX,'K-EFFECTIVE',ILEN,ITYLCM)
      IF(ILEN.GT.0) THEN
         CALL LCMGET(IPFLUX,'K-EFFECTIVE',FKEFF)
         CALL LCMPUT(IPMAC2,'K-EFFECTIVE',1,2,FKEFF)
      ENDIF
      CALL LCMLEN(IPFLUX,'NORM-FS',ILEN,ITYLCM)
      IF(ILEN.GT.0) THEN
         CALL LCMGET(IPFLUX,'NORM-FS',NORM)
         CALL LCMPUT(IPMAC2,'NORM-FS',1,2,NORM)
         CALL LCMGET(IPFLUX,'MATCOD',MATCOD)
         CALL LCMPUT(IPMAC2,'MATCOD',NEL,1,MATCOD)
      ENDIF
      CALL LCMLEN(IPFLUX,'FLUXC',ILEN,ITYLCM)
      IF(ILEN.GT.0) THEN
         CALL LCMGET(IPFLUX,'FLUXC',FLUXC)
         CALL LCMPUT(IPMAC2,'FLUXC',NEL,2,FLUXC)
         CALL LCMGET(IPFLUX,'ECUTOFF',ECUTOFF)
         CALL LCMPUT(IPMAC2,'ECUTOFF',1,2,ECUTOFF)
      ENDIF
*
   20 CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
      IF(INDIC.NE.3) CALL XABORT('OUTDRV: CHARACTER DATA EXPECTED.')
*
      IF(TEXT4.EQ.'EDIT') THEN
         CALL REDGET(INDIC,IMPX,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('OUT: INTEGER DATA EXPECTED.')
      ELSE IF(TEXT4.EQ.'MODE') THEN
         CALL REDGET(INDIC,LMOD,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('OUT: INTEGER DATA EXPECTED.')
      ELSE IF(TEXT4.EQ.'DIRE') THEN
         IADJ=0
      ELSE IF(TEXT4.EQ.'PROD') THEN
         IADJ=1
      ELSE
         CALL OUTFLX(IPFLUX,0,NGRP,NUN,LMOD,IMPX,EVECT)
         GO TO 40
      ENDIF
      GO TO 20
*
   30 CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
      IF(INDIC.NE.3) CALL XABORT('OUTDRV: CHARACTER DATA EXPECTED.')
*
   40 IF(TEXT4.EQ.'POWR') THEN
*        NORMALIZATION TO A GIVEN FISSION POWER.
         CALL REDGET (INDIC,NITMA,POWER,TEXT4,DFLOTT)
         IF(INDIC.NE.2) CALL XABORT('OUTDRV: REAL DATA EXPECTED.')
*        NORMALIZATION FACTOR FOR THE DIRECT FLUX.
         ZNORM=0.0D0
         JPMAC1=LCMGID(IPMAC1,'GROUP')
         DO 60 IGR=1,NGRP
         KPMAC1=LCMGIL(JPMAC1,IGR)
         CALL LCMLEN(KPMAC1,'H-FACTOR',LENGT,ITYLCM)
         IF(LENGT.GT.0) THEN
            CALL LCMGET(KPMAC1,'H-FACTOR',SGD)
         ELSE
            WRITE(6,'(/43H OUTDRV: *** WARNING *** NO H-FACTOR FOUND ,
     1      28HON LCM. USE NU*SIGF INSTEAD.)')
            ALLOCATE(ZUFIS(NBMIX,NBFIS))
            CALL XDRSET(SGD,NBMIX,0.0)
            CALL LCMGET(KPMAC1,'NUSIGF',ZUFIS)
            DO IBM=1,NBMIX
              DO IFISS=1,NBFIS
                SGD(IBM)=SGD(IBM)+ZUFIS(IBM,IFISS)
              ENDDO
            ENDDO
            DEALLOCATE(ZUFIS)
         ENDIF
         DO 50 K=1,NEL
         L=MAT(K)
         IF((L.EQ.0).OR.(IDL(K).EQ.0)) GO TO 50
         ZNORM=ZNORM+EVECT(IDL(K),IGR)*VOL(K)*SGD(L)
   50    CONTINUE
   60    CONTINUE
         ZNORM=POWER/ZNORM
         WRITE(6,300) ' DIRECT',ZNORM
         DO 80 IGR=1,NGRP
         DO 70 I=1,NUN
         EVECT(I,IGR)=EVECT(I,IGR)*REAL(ZNORM)
   70    CONTINUE
   80    CONTINUE
      ELSE IF(TEXT4.EQ.'SOUR') THEN
*        NORMALIZATION TO A GIVEN SOURCE INTENSITY.
         CALL REDGET (INDIC,NITMA,SNUMB,TEXT4,DFLOTT)
         IF(INDIC.NE.2) CALL XABORT('OUTDRV: REAL DATA EXPECTED.')
*        NORMALIZATION FACTOR FOR THE DIRECT FLUX.
         ZNORM=0.0D0
         JPMAC1=LCMGID(IPMAC1,'GROUP')
         DO 100 IGR=1,NGRP
         KPMAC1=LCMGIL(JPMAC1,IGR)
         CALL LCMLEN(KPMAC1,'FIXE',LENGT,ITYLCM)
         IF(LENGT.EQ.0) THEN
            CALL LCMLIB(KPMAC1)
            CALL XABORT('OUTDRV: SOURCE RECORD MISSING IN MACROLIB.')
         ENDIF
         CALL LCMGET(KPMAC1,'FIXE',SGD)
         DO 90 K=1,NEL
         L=MAT(K)
         IF(L.GT.0) ZNORM=ZNORM+VOL(K)*SGD(L)
   90    CONTINUE
  100    CONTINUE
         ZNORM=SNUMB/ZNORM
         WRITE(6,305) ' DIRECT',ZNORM
         DO 120 IGR=1,NGRP
         DO 110 I=1,NUN
         EVECT(I,IGR)=EVECT(I,IGR)*REAL(ZNORM)
  110    CONTINUE
  120    CONTINUE
      ELSE IF(TEXT4.EQ.'COND') THEN
         NGCOND=0
         CALL REDGET (INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.EQ.3) THEN
           IF(TEXT4.EQ.'NONE') THEN
             NGCOND=NGRP
             DO IGR=1,NGRP
               IGCOND(IGR)=IGR
             ENDDO
             GO TO 30
           ENDIF
           NGCOND=1
           IGCOND(NGCOND)=NGRP
           GO TO 40
         ELSE IF(INDIC.EQ.1) THEN
  130      IF(NITMA.GT.NGRP) NITMA=NGRP
           NGCOND=NGCOND+1
           IGCOND(NGCOND)=NITMA
           CALL REDGET (INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
           IF(INDIC.EQ.1) THEN
             GO TO 130
           ELSE IF(INDIC.EQ.3) THEN
             GO TO 30
           ELSE
             CALL XABORT('OUTDRV: INTEGER OR CHARACTER DATA EXPECTED.')
           ENDIF
         ELSE
           CALL XABORT('OUTDRV: INTEGER OR CHARACTER DATA EXPECTED.')
         ENDIF
      ELSE IF(TEXT4.EQ.'INTG') THEN
*        COMPUTE AND DISPLAY THE MACRO-ZONE REACTION RATES.
*        READ THE MACRO-ZONES DEFINITION.
         IF(IMPX.GT.0) WRITE(6,330) (IGCOND(IG),IG=1,NGCOND)
         CALL OUTHOM (MAXNEL,IPGEOM,IMPX,NEL,IELEM,ICOL,HTRACK,MAT,NZS,
     1   IHOM)
         IF(NZS.GT.NEL) CALL XABORT('OUTDRV: INVALID VALUE OF NZS.')
         IF(IMPX.GT.0) WRITE(6,320) TITR
         IF(IADJ.EQ.0) THEN
           CALL OUTAUX(IPMAC1,IPMAC2,NBMIX,NL,NBFIS,NGRP,NEL,NUN,NALBP,
     1     NZS,NGCOND,MAT,VOL,IDL,EVECT,IHOM,IGCOND,IMPX)
         ELSE IF(IADJ.EQ.1) THEN
           ALLOCATE(ADECT(NUN,NGRP))
           CALL OUTFLX(IPFLUX,1,NGRP,NUN,LMOD,IMPX,ADECT)
           CALL OUTPRO(IPMAC1,IPMAC2,NBMIX,NL,NBFIS,NGRP,NEL,NUN,NALBP,
     1     NZS,NGCOND,MAT,VOL,IDL,EVECT,ADECT,IHOM,IGCOND,IMPX)
           DEALLOCATE(ADECT)
         ENDIF
      ELSE IF(TEXT4.EQ.';') THEN
         CALL KDRCPU(TK2)
         TKR=TK2-TK1
         WRITE(6,310) TKR
         GO TO 140
      ELSE
         CALL XABORT('OUTDRV: '//TEXT4//' IS AN INVALID KEY WORD.')
      ENDIF
      GO TO 30
*----
*  SCRATCH STORAGE DEALLOCATION
*----
  140 DEALLOCATE(FLUXC,SGD,EVECT,IGCOND,IHOM)
      RETURN
*
  300 FORMAT(/9H OUTDRV: ,A7,28H FLUX NORMALIZATION FACTOR =,1P,E13.5)
  305 FORMAT(/9H OUTDRV: ,A7,30H SOURCE NORMALIZATION FACTOR =,1P,E13.5)
  310 FORMAT(/49H OUTDRV: CPU TIME FOR REACTION RATE CALCULATION =,F7.3)
  320 FORMAT(/12H OUTDRV: ***,A72,3H***)
  330 FORMAT(/20H CONDENSATION INDEX:/(1X,14I5))
      END
