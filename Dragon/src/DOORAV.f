*DECK DOORAV
      SUBROUTINE DOORAV (CDOOR,JPSYS,NPSYS,IPTRK,IFTRAK,IMPX,NGRP,NREG,
     1 NBMIX,NANI,NW,MAT,VOL,KNORM,LEAKSW,TITR,NALBP,ISTRM)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Assembly of system matrices. Vectorial version.
*
*Copyright:
* Copyright (C) 2004 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* CDOOR   name of the geometry/solution operator.
* JPSYS   pointer to the PIJ LCM object (L_PIJ signature). JPSYS is
*         a list of directories.
* NPSYS   index array pointing to the JPSYS list component corresponding
*         to each energy group. Set to zero if a group is not to be
*         processed. Usually, NPSYS(I)=I.
* IPTRK   pointer to the tracking (L_TRACK signature).
* IFTRAK  unit of the sequential binary tracking file.
* IMPX    print flag (equal to zero for no print).
* NGRP    number of energy groups.
* NREG    total number of merged blocks for which specific values
*         of the neutron flux and reactions rates are required.
* NBMIX   number of mixtures (NBMIX=max(MAT(i))).
* NANI    number of Legendre orders.
* NW      type of weighting for P1 cross section info (=0: P0 ; =1: P1).
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* KNORM   normalization scheme.
* LEAKSW  leakage flag (=.true. if neutron leakage through external
*         boundary is present).
* TITR    title.
* NALBP   number of physical albedos.
* ISTRM   type of streaming effect:
*         =1 no streaming effect;
*         =2 isotropic streaming effect;
*         =3 anisotropic streaming effect.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER CDOOR*12,TITR*72
      LOGICAL LEAKSW
      TYPE(C_PTR) JPSYS,IPTRK
      INTEGER NPSYS(NGRP),IFTRAK,IMPX,NGRP,NREG,NBMIX,NANI,NW,MAT(NREG),
     > KNORM,NALBP,ISTRM
      REAL VOL(NREG)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40)
      INTEGER ISTATE(NSTATE)
      LOGICAL LBIHET
      CHARACTER TEXT12*12
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MAT2
      REAL, ALLOCATABLE, DIMENSION(:) :: VOL2,SGAR,SGAS,SGAD,ALBP,GAMMA
      TYPE(C_PTR) KPSYS
*
      ALB(X)=0.5*(1.0-X)/(1.0+X)
*
      IF(IMPX.GT.5) THEN
         WRITE(6,'(/36H DOORAV: ASSEMBLY OF SYSTEM MATRICES//9X,A72)')
     1   TITR
         WRITE(6,'(/30H DOORAV: NORMALIZATION SCHEME=,I2,9H LEAKAGE ,
     1   7HSWITCH=,L2)') KNORM,LEAKSW
      ENDIF
*----
*  DOUBLE HETEROGENEITY TREATMENT
*----
      NREGAR=0
      NBMIXG=0
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      LBIHET=ISTATE(40).NE.0
      IF(LBIHET) THEN
         ALLOCATE(MAT2(NREG),VOL2(NREG))
         DO I=1,NREG
            MAT2(I)=MAT(I)
            VOL2(I)=VOL(I)
         ENDDO
         NREGAR=NREG
         NBMIXG=NBMIX
         CALL DOORAB(CDOOR,JPSYS,NPSYS,IPTRK,IMPX,NGRP,NREG,NBMIX,NANI,
     1   MAT,VOL)
      ENDIF
*
      IF (CDOOR.EQ.'SYBIL') THEN
         ALLOCATE(SGAR(NBMIX+1),SGAS((NBMIX+1)*NANI))
         DO 100 IGR=1,NGRP
            IOFSET=NPSYS(IGR)
            IF(IOFSET.NE.0) THEN
               KPSYS=LCMGIL(JPSYS,IOFSET)
               IF(LBIHET) CALL LCMSIX(KPSYS,'BIHET',1)
               CALL LCMLEN(KPSYS,'DRAGON-TXSC',ILONG,ITYLCM)
               IF(ILONG.NE.NBMIX+1) CALL XABORT('DOORAV: INVALID TXSC '
     1         //'LENGTH(1).')
               CALL LCMGET(KPSYS,'DRAGON-TXSC',SGAR)
               CALL LCMGET(KPSYS,'DRAGON-S0XSC',SGAS)
               CALL SYBILA(KPSYS,IPTRK,IMPX,NREG,NBMIX,MAT,SGAR,SGAS)
               IF(LBIHET) CALL LCMSIX(KPSYS,' ',2)
            ENDIF
 100     CONTINUE
         DEALLOCATE(SGAS,SGAR)
      ELSE IF (CDOOR.EQ.'SN') THEN
         CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
         ITYPE=ISTATE(6)
         IF(ISTATE(19).EQ.1) THEN
*           SYNTHETIC ACCELERATION.
            IF(IMPX.GT.0) WRITE (6,'(/29H DOORAV: SYNTHETIC ACCELERATI,
     1      20HON ASSEMBLY FOLLOWS:)')
            CALL LCMSIX(IPTRK,'DSA',1)
            ALLOCATE(SGAR((NBMIX+1)*(NW+1)),SGAS((NBMIX+1)*NANI),
     1      SGAD(NBMIX+1))
            DO 150 IGR=1,NGRP
               IOFSET=NPSYS(IGR)
               IF(IOFSET.NE.0) THEN
                  KPSYS=LCMGIL(JPSYS,IOFSET)
                  CALL XDISET(ISTATE,NSTATE,0)
                  CALL LCMPUT(KPSYS,'STATE-VECTOR',NSTATE,1,ISTATE)
                  IF(LBIHET) CALL LCMSIX(KPSYS,'BIHET',1)
                  CALL LCMLEN(KPSYS,'DRAGON-TXSC',ILONG,ITYLCM)
                  IF(ILONG.NE.NBMIX+1) CALL XABORT('DOORAV: INVALID TX'
     1            //'SC LENGTH(2).')
                  CALL LCMGET(KPSYS,'DRAGON-TXSC',SGAR)
                  DO 110 IW=2,MIN(NW+1,10)
                  IOF=(NBMIX+1)*(IW-1)+1
                  WRITE(TEXT12,'(8HDRAGON-T,I1,3HXSC)') IW-1
                  CALL LCMLEN(KPSYS,TEXT12,ILONG,ITYLCM)
                  IF(ILONG.EQ.NBMIX+1) THEN
                     CALL LCMGET(KPSYS,TEXT12,SGAR(IOF))
                  ELSE
                     CALL LCMGET(KPSYS,'DRAGON-TXSC',SGAR(IOF))
                  ENDIF
  110             CONTINUE
                  CALL LCMGET(KPSYS,'DRAGON-S0XSC',SGAS)
                  IF((ITYPE.EQ.7).OR.(ITYPE.EQ.9)) THEN
                     CALL TRIVA(KPSYS,IPTRK,IMPX,NREG,NBMIX,NANI,NW,MAT,
     1                          VOL,SGAR,SGAS,SGAD)
                  ELSE
                     CALL BIVAA(KPSYS,IPTRK,IMPX,NREG,NBMIX,NANI,NW,MAT,
     1                          VOL,SGAR,SGAS,SGAD)
                  ENDIF
                  IF(LBIHET) CALL LCMSIX(KPSYS,' ',2)
               ENDIF
 150        CONTINUE
            DEALLOCATE(SGAD,SGAS,SGAR)
            CALL LCMSIX(IPTRK,' ',2)
         ENDIF
      ELSE IF (CDOOR.EQ.'BIVAC') THEN
         CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
         IELEM=ISTATE(8)
         ICOL=ISTATE(9)
         NLF=ISTATE(14)
         ISCAT=ISTATE(16)
         ALLOCATE(SGAR((NBMIX+1)*(NW+1)),SGAS((NBMIX+1)*NANI),
     1   SGAD(NBMIX+1))
         DO 190 IGR=1,NGRP
            IOFSET=NPSYS(IGR)
            IF(IOFSET.NE.0) THEN
               KPSYS=LCMGIL(JPSYS,IOFSET)
*              SAVE ALBEDO FUNCTIONS ON KPSYS
               IF(NALBP.GT.0) THEN
                 ALLOCATE(ALBP(NALBP),GAMMA(NALBP))
                 CALL LCMGET(KPSYS,'ALBEDO',ALBP)
                 DO IALB=1,NALBP
                   IF((IELEM.LT.0).OR.(ICOL.EQ.4)) THEN
                     GAMMA(IALB)=ALB(ALBP(IALB))
                   ELSE IF(ALBP(IALB).NE.1.0) THEN
                     GAMMA(IALB)=1.0/ALB(ALBP(IALB))
                   ELSE
                     GAMMA(IALB)=1.0E20
                   ENDIF
                 ENDDO
                 CALL LCMPUT(KPSYS,'ALBEDO-FU',NALBP,2,GAMMA)
                 DEALLOCATE(GAMMA,ALBP)
               ENDIF
*
               IF(LBIHET) CALL LCMSIX(KPSYS,'BIHET',1)
               CALL LCMLEN(KPSYS,'DRAGON-TXSC',ILONG,ITYLCM)
               IF(ILONG.NE.NBMIX+1) CALL XABORT('DOORAV: INVALID TXSC '
     1         //'LENGTH(3).')
               CALL LCMGET(KPSYS,'DRAGON-TXSC',SGAR)
               IF(NLF.EQ.0) THEN
                  CALL LCMGET(KPSYS,'DRAGON-DIFF',SGAD)
               ELSE IF(ISCAT.LT.0) THEN
                  CALL LCMGET(KPSYS,'DRAGON-DIFF',SGAD)
                  SGAR(NBMIX+2)=1.0E10
                  DO 180 IMIX=1,NBMIX
                  SGAR(NBMIX+2+IMIX)=1.0/(3.0*SGAD(IMIX+1))
 180              CONTINUE
               ELSE IF(ISCAT.GT.0) THEN
                  DO 185 IW=2,MIN(NW+1,10)
                  IOF=(NBMIX+1)*(IW-1)+1
                  WRITE(TEXT12,'(8HDRAGON-T,I1,3HXSC)') IW-1
                  CALL LCMLEN(KPSYS,TEXT12,ILONG,ITYLCM)
                  IF(ILONG.NE.0) THEN
                     CALL LCMGET(KPSYS,'DRAGON-T1XSC',SGAR(IOF))
                  ELSE
                     CALL LCMGET(KPSYS,'DRAGON-TXSC',SGAR(IOF))
                  ENDIF
 185              CONTINUE
               ENDIF
               CALL LCMGET(KPSYS,'DRAGON-S0XSC',SGAS)
               CALL BIVAA(KPSYS,IPTRK,IMPX,NREG,NBMIX,NANI,NW,MAT,VOL,
     1                    SGAR,SGAS,SGAD)
               IF(LBIHET) CALL LCMSIX(KPSYS,' ',2)
            ENDIF
 190     CONTINUE
         DEALLOCATE(SGAD,SGAS,SGAR)
      ELSE IF (CDOOR.EQ.'TRIVAC') THEN
         CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
         ICHX=ISTATE(12)
         NLF=ISTATE(30)
         ISCAT=ISTATE(32)
         ALLOCATE(SGAR((NBMIX+1)*2),SGAS((NBMIX+1)*NANI),SGAD(NBMIX+1))
         DO 210 IGR=1,NGRP
            IOFSET=NPSYS(IGR)
            IF(IOFSET.NE.0) THEN
               KPSYS=LCMGIL(JPSYS,IOFSET)
*              SAVE ALBEDO FUNCTIONS ON KPSYS
               IF(NALBP.GT.0) THEN
                 ALLOCATE(ALBP(NALBP),GAMMA(NALBP))
                 CALL LCMGET(KPSYS,'ALBEDO',ALBP)
                 DO IALB=1,NALBP
                   IF(ICHX.NE.2) THEN
                     GAMMA(IALB)=ALB(ALBP(IALB))
                   ELSE IF(ALBP(IALB).NE.1.0) THEN
                     GAMMA(IALB)=1.0/ALB(ALBP(IALB))
                   ELSE
                     GAMMA(IALB)=1.0E20
                   ENDIF
                 ENDDO
                 CALL LCMPUT(KPSYS,'ALBEDO-FU',NALBP,2,GAMMA)
                 DEALLOCATE(GAMMA,ALBP)
               ENDIF
*
               CALL XDISET(ISTATE,NSTATE,0)
               CALL LCMPUT(KPSYS,'STATE-VECTOR',NSTATE,1,ISTATE)
               IF(LBIHET) CALL LCMSIX(KPSYS,'BIHET',1)
               CALL LCMLEN(KPSYS,'DRAGON-TXSC',ILONG,ITYLCM)
               IF(ILONG.NE.NBMIX+1) CALL XABORT('DOORAV: INVALID TXSC '
     1         //'LENGTH(4).')
               CALL LCMGET(KPSYS,'DRAGON-TXSC',SGAR)
               IF(NLF.EQ.0) THEN
                  CALL LCMGET(KPSYS,'DRAGON-DIFF',SGAD)
               ELSE IF(ISCAT.LT.0) THEN
                  CALL LCMGET(KPSYS,'DRAGON-DIFF',SGAD)
                  SGAR(NBMIX+2)=1.0E10
                  DO 200 IMIX=1,NBMIX
                  SGAR(NBMIX+2+IMIX)=1.0/(3.0*SGAD(IMIX+1))
 200              CONTINUE
               ELSE IF(ISCAT.GT.0) THEN
                  DO 205 IW=2,MIN(NW+1,10)
                  IOF=(NBMIX+1)*(IW-1)+1
                  WRITE(TEXT12,'(8HDRAGON-T,I1,3HXSC)') IW-1
                  CALL LCMLEN(KPSYS,TEXT12,ILONG,ITYLCM)
                  IF(ILONG.NE.0) THEN
                     CALL LCMGET(KPSYS,'DRAGON-T1XSC',SGAR(IOF))
                  ELSE
                     CALL LCMGET(KPSYS,'DRAGON-TXSC',SGAR(IOF))
                  ENDIF
 205              CONTINUE
               ENDIF
               CALL LCMGET(KPSYS,'DRAGON-S0XSC',SGAS)
               CALL TRIVA(KPSYS,IPTRK,IMPX,NREG,NBMIX,NANI,NW,MAT,VOL,
     1                    SGAR,SGAS,SGAD)
               IF(LBIHET) CALL LCMSIX(KPSYS,' ',2)
            ENDIF
 210     CONTINUE
         DEALLOCATE(SGAD,SGAS,SGAR)
      ELSE IF ((CDOOR.EQ.'MCCG').OR.(CDOOR.EQ.'EXCELL')) THEN
         CALL MCCGA(JPSYS,NPSYS,IPTRK,IFTRAK,IMPX,NGRP,NBMIX,NANI,
     1   NALBP,ISTRM)
      ELSE
         CALL XABORT('DOORAV: UNKNOWN DOOR:'//CDOOR//'.')
      ENDIF
*----
*  DOUBLE HETEROGENEITY TREATMENT
*----
      IF(LBIHET) THEN
         NREG=NREGAR
         NBMIX=NBMIXG
         DO I=1,NREG
            MAT(I)=MAT2(I)
            VOL(I)=VOL2(I)
         ENDDO
         DEALLOCATE(MAT2,VOL2)
      ENDIF
      RETURN
      END
