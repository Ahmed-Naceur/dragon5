*DECK DOORPV
      SUBROUTINE DOORPV (CDOOR,JPSYS,NPSYS,IPTRK,IFTRAK,IMPX,NGRP,NREG,
     1 NBMIX,NANI,MAT,VOL,KNORM,IPIJK,LEAKSW,LNORM,TITR,NALBP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of the collision probabilities. Vectorial version.
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
* NANI    number of Legendre orders (usually equal to one).
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* KNORM   normalization scheme for PIJ matrices.
* IPIJK   pij option (=1 pij, =4 pijk).
* LEAKSW  leakage flag (=.true. if neutron leakage through external
*         boundary is present).
* LNORM   logical switch for removing leakage from collision
*         probabilities and keeping the PIS information.
* TITR    title.
* NALBP   number of physical albedos.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER CDOOR*12,TITR*72
      LOGICAL LEAKSW,LNORM
      TYPE(C_PTR) JPSYS,IPTRK
      INTEGER NPSYS(NGRP),IFTRAK,IMPX,NGRP,NREG,NBMIX,NANI,MAT(NREG),
     > KNORM,IPIJK,NALBP
      REAL VOL(NREG)
      INTEGER NNPSYS(1)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40,IOUT=6)
      CHARACTER CMS(3)*1
      INTEGER ISTATE(NSTATE)
      LOGICAL LBIHET,LNXT
      TYPE(C_PTR) KPSYS
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MAT2
      REAL, ALLOCATABLE, DIMENSION(:) :: VOL2,PGAR,SGAR,SGAS,PROBKS,
     > PIS
      REAL, ALLOCATABLE, DIMENSION(:,:) :: PGARG,ALBP
      REAL, POINTER, DIMENSION(:,:) :: PREG
      TYPE(C_PTR) :: PREG_PTR
*----
*  DATA STATEMENT AND INLINE FUNCTION
*----
      SAVE CMS
      DATA CMS/'1','2','3'/
      INDPOS(I,J)=MAX(I,J)*(MAX(I,J)-1)/2+MIN(I,J)
*----
*  DOUBLE HETEROGENEITY TREATMENT
*----
      NNPSYS(1)=1
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      LBIHET=ISTATE(40).NE.0
      NREGAR=0
      NBMIXG=0
      IF(LBIHET) THEN
        ALLOCATE(MAT2(NREG),VOL2(NREG))
        DO I=1,NREG
          MAT2(I)=MAT(I)
          VOL2(I)=VOL(I)
        ENDDO
        NREGAR=NREG
        NBMIXG=NBMIX
        CALL DOORAB(CDOOR,JPSYS,NPSYS,IPTRK,IMPX,NGRP,NREG,NBMIX,NANI,
     1  MAT,VOL)
      ENDIF
*
      NELPIJ=NREG*(NREG+1)/2
      NB1=NBMIX+1
      IF(CDOOR.EQ.'EXCELL') THEN
        CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
        INSB=ISTATE(22)
        LNXT=(ISTATE(7).EQ.4)
        IF(INSB.GE.1) GO TO 110
      ENDIF
*----
*  COMPUTE THE REDUCED PIJ MATRIX -- NON-VECTORIAL ALGORITHM.
*----
      ALLOCATE(PGAR(IPIJK*NELPIJ),SGAR(NB1),SGAS(NB1*NANI))
      ALLOCATE(ALBP(NALBP,1))
      DO 100 IGR=1,NGRP
      IOFSET=NPSYS(IGR)
      IF(IOFSET.NE.0) THEN
        IF(IMPX.GT.10) WRITE(IOUT,'(/25H DOORPV: PROCESSING GROUP,I5,
     >  6H WITH ,A,1H.)') IGR,CDOOR
        KPSYS=LCMGIL(JPSYS,IOFSET)
        IF(LBIHET) CALL LCMSIX(KPSYS,'BIHET',1)
        CALL LCMGET(KPSYS,'DRAGON-TXSC',SGAR)
        CALL LCMGET(KPSYS,'DRAGON-S0XSC',SGAS)
        IF(NALBP.GT.0) CALL LCMGET(KPSYS,'ALBEDO',ALBP)
        IF(CDOOR.EQ.'EXCELL') THEN
          IF(IPIJK.EQ.4) THEN
            NNREG=NREG*NREG
            NPST=3*NNREG
            ALLOCATE(PROBKS(NPST))
          ELSE
            NPST=1
          ENDIF
          CALL EXCELP(IPTRK,IFTRAK,IMPX,NREG,NBMIX,NANI,MAT,VOL,KNORM,
     >    SGAR,SGAS,NELPIJ,IPIJK,PGAR,LEAKSW,1,NNPSYS,NPST,PROBKS,
     >    TITR,NALBP,ALBP)
          IF(IPIJK.EQ.4) THEN
            CALL LCMPUT(KPSYS,'DRAGON1P*SCT',NNREG,2,PROBKS)
            CALL LCMPUT(KPSYS,'DRAGON2P*SCT',NNREG,2,PROBKS(NNREG+1))
            CALL LCMPUT(KPSYS,'DRAGON3P*SCT',NNREG,2,PROBKS(2*NNREG+1))
            DEALLOCATE(PROBKS)
          ENDIF
        ELSE IF(CDOOR(:5).EQ.'SYBIL') THEN
          CALL SYBILP(IPTRK,IMPX,NREG,NBMIX,MAT,VOL,SGAR,SGAS,NELPIJ,
     >    PGAR,LEAKSW)
        ELSE
          CALL XABORT('DOORPV: UNKNOWN PIJ DOOR NAMED '//CDOOR//'.')
        ENDIF
*----
*  REMOVE LEAKAGE FROM THE SCATTERING-REDUCED CP MATRIX.
*----
        IF(LNORM) THEN
          ALLOCATE(PIS(NREG))
          CALL XDRNRM(NREG,NBMIX,MAT,VOL,SGAR(2),SGAS(2),PGAR,PIS)
          IF(LEAKSW) CALL LCMPUT(KPSYS,'DRAGON-WIS',NREG,2,PIS)
          DEALLOCATE(PIS)
        ENDIF
*----
*  FORMAT THE REDUCED PIJ MATRIX AS A SQUARE MATRIX.
*----
        PREG_PTR=LCMARA(NREG*NREG)
        CALL C_F_POINTER(PREG_PTR,PREG,(/ NREG,NREG /))
        DO 15 I=1,NREG
        FACT=1.0/VOL(I)
        DO 10 J=1,NREG
        PREG(I,J)=PGAR(INDPOS(I,J))*FACT
   10   CONTINUE
   15   CONTINUE
        CALL LCMPPD(KPSYS,'DRAGON-PCSCT',NREG*NREG,2,PREG_PTR)
        IF(IPIJK.EQ.4) THEN
          DO 30 IJKS=1,3
          PREG_PTR=LCMARA(NREG*NREG)
          CALL C_F_POINTER(PREG_PTR,PREG,(/ NREG,NREG /))
          DO 25 I=1,NREG
          FACT=1.0/VOL(I)
          DO 20 J=1,NREG
          KS=NELPIJ*IJKS+INDPOS(I,J)
          PREG(I,J)=PGAR(KS)*FACT
   20     CONTINUE
   25     CONTINUE
          CALL LCMPPD(KPSYS,'DRAGON'//CMS(IJKS)//'PCSCT',NREG*NREG,2,
     >    PREG_PTR)
   30     CONTINUE
        ENDIF
        IF(LBIHET) CALL LCMSIX(KPSYS,' ',2)
      ENDIF
  100 CONTINUE
      DEALLOCATE(ALBP,SGAS,SGAR,PGAR)
      GO TO 210
*----
*  COMPUTE THE REDUCED PIJ MATRIX -- VECTORIAL ALGORITHM FOR EXCELP.
*----
  110 ALLOCATE(SGAR(NB1*NGRP),SGAS(NB1*NANI*NGRP))
      ALLOCATE(ALBP(NALBP,NGRP))
      DO 120 IGR=1,NGRP
      IOFSET=NPSYS(IGR)
      IF(IOFSET.NE.0) THEN
        KPSYS=LCMGIL(JPSYS,IOFSET)
        IF(LBIHET) CALL LCMSIX(KPSYS,'BIHET',1)
        CALL LCMGET(KPSYS,'DRAGON-TXSC',SGAR((IGR-1)*NB1+1))
        CALL LCMGET(KPSYS,'DRAGON-S0XSC',SGAS((IGR-1)*NB1*NANI+1))
        IF(NALBP.GT.0) CALL LCMGET(KPSYS,'ALBEDO',ALBP(1,IGR))
        IF(LBIHET) CALL LCMSIX(KPSYS,' ',2)
      ENDIF
  120 CONTINUE
      ALLOCATE(PGARG(IPIJK*NELPIJ,NGRP))
      IF(LNXT.OR.(INSB.EQ.1)) THEN
*       --- ALLG VECTORIZATION
        IF(IPIJK.EQ.4) THEN
          NNREG=NREG*NREG
          NPST=3*NNREG
          ALLOCATE(PROBKS(NPST*NGRP))
        ELSE
          NPST=1
        ENDIF
        CALL EXCELP(IPTRK,IFTRAK,IMPX,NREG,NBMIX,NANI,MAT,VOL,KNORM,
     >  SGAR,SGAS,NELPIJ,IPIJK,PGARG,LEAKSW,NGRP,NPSYS,NPST,PROBKS,TITR,
     >  NALBP,ALBP)
      ELSE IF(INSB.EQ.2) THEN
*       --- XCLL VECTORIZATION
        IF(IPIJK.NE.1) CALL XABORT('DOORPV: INVALID VALUE OF IPIJK')
        CALL PIJXL3(IPTRK,IMPX,NGRP,NANI,NBMIX,NPSYS,KNORM,LEAKSW,
     >  SGAR,SGAS,NELPIJ,PGARG)
      ELSE
         CALL XABORT('DOORPV: INVALID VALUE OF INSB')
      ENDIF
      KPST=0
      DO 200 IGR=1,NGRP
      IOFSET=NPSYS(IGR)
      IF(IOFSET.NE.0) THEN
        KPSYS=LCMGIL(JPSYS,IOFSET)
        IF(LBIHET) CALL LCMSIX(KPSYS,'BIHET',1)
        IF(IPIJK.EQ.4) THEN
          CALL LCMPUT(KPSYS,'DRAGON1P*SCT',NNREG,2,PROBKS(KPST+1))
          CALL LCMPUT(KPSYS,'DRAGON2P*SCT',NNREG,2,PROBKS(KPST+NNREG+1))
          CALL LCMPUT(KPSYS,'DRAGON3P*SCT',NNREG,2,
     >                PROBKS(KPST+2*NNREG+1))
        ENDIF
*----
*  REMOVE LEAKAGE FROM THE SCATTERING-REDUCED CP MATRIX.
*----
        IF(LNORM) THEN
          ALLOCATE(PIS(NREG))
          CALL XDRNRM(NREG,NBMIX,MAT,VOL,SGAR((IGR-1)*NB1+2),
     >    SGAS((IGR-1)*NB1*NANI+2),PGARG(1,IGR),PIS)
          IF(LEAKSW) CALL LCMPUT(KPSYS,'DRAGON-WIS',NREG,2,PIS)
          DEALLOCATE(PIS)
        ENDIF
*----
*  FORMAT THE REDUCED PIJ MATRIX AS A SQUARE MATRIX.
*----
        PREG_PTR=LCMARA(NREG*NREG)
        CALL C_F_POINTER(PREG_PTR,PREG,(/ NREG,NREG /))
        DO 135 I=1,NREG
        FACT=1.0/VOL(I)
        DO 130 J=1,NREG
        PREG(I,J)=PGARG(INDPOS(I,J),IGR)*FACT
  130   CONTINUE
  135   CONTINUE
        CALL LCMPPD(KPSYS,'DRAGON-PCSCT',NREG*NREG,2,PREG_PTR)
        IF(IPIJK.EQ.4) THEN
          DO 150 IJKS=1,3
          PREG_PTR=LCMARA(NREG*NREG)
          CALL C_F_POINTER(PREG_PTR,PREG,(/ NREG,NREG /))
          DO 145 I=1,NREG
          FACT=1.0/VOL(I)
          DO 140 J=1,NREG
          KS=NELPIJ*IJKS+INDPOS(I,J)
          PREG(I,J)=PGARG(KS,IGR)*FACT
  140     CONTINUE
  145     CONTINUE
          CALL LCMPPD(KPSYS,'DRAGON'//CMS(IJKS)//'PCSCT',NREG*NREG,2,
     >    PREG_PTR)
  150     CONTINUE
        ENDIF
        IF(LBIHET) CALL LCMSIX(KPSYS,' ',2)
      ENDIF
      KPST=KPST+NPST
  200 CONTINUE
      IF(IPIJK.EQ.4) DEALLOCATE(PROBKS)
      DEALLOCATE(ALBP,SGAS,SGAR,PGARG)
*----
*  DOUBLE HETEROGENEITY TREATMENT
*----
  210 IF(LBIHET) THEN
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
