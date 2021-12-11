*DECK DOORFV
      SUBROUTINE DOORFV(CDOOR,IPSYS,NPSYS,IPTRK,IFTRAK,IMPX,NGRP,NMAT,
     1 IDIR,NREG,NUN,IPHASE,LEXAC,MAT,VOL,KEYFLX,TITR,SUNKNO,FUNKNO,
     2 IPMACR,REBFLG)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of the flux. Vectorial version. Multigroup rebalancing
* option.
*
*Copyright:
* Copyright (C) 2004 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Le Tellier
*
*Parameters: input
* CDOOR   name of the geometry/solution operator.
* IPSYS   pointer to the assembly LCM object (L_PIJ signature). IPSYS is
*         a list of directories.
* NPSYS   index array pointing to the IPSYS list component corresponding
*         to each energy group. Set to zero if a group is not to be
*         processed. Usually, NPSYS(I)=I.
* IPTRK   pointer to the tracking (L_TRACK signature).
* IFTRAK  unit of the sequential binary tracking file.
* IMPX    print flag (equal to zero for no print).
* NGRP    number of energy groups.
* NMAT    number of mixtures in the internal library.
* IDIR    directional collision probability flag:
*         =0 for pij or wij;
*         =k for pijk or wijk k=1,2,3.
*         direction of fundamental current for TIBERE with MoC 
*         (=0,1,2,3).  
* NREG    total number of merged blocks for which specific values
*         of the neutron flux and reactions rates are required.
* NUN     total number of unknowns in vectors SUNKNO and FUNKNO.
* IPHASE  type of flux solution (=1: use a native flux solution door;
*         =2: use collision probabilities).
* LEXAC   type of exponential function calculation (=.false. to compute
*         exponential functions using tables).
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* KEYFLX  index of L-th order flux components in unknown vector.
* TITR    title.
* SUNKNO  input source vector.
* FUNKNO  unknown vector.
* IPMACR  pointer to the macrolib LCM object.
* REBFLG  ACA or SCR rebalancing flag.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPSYS,IPTRK,IPMACR
      CHARACTER CDOOR*12,TITR*72
      INTEGER NPSYS(NGRP),IFTRAK,IMPX,NGRP,NMAT,IDIR,NREG,NUN,IPHASE,
     > MAT(NREG),KEYFLX(NREG)
      REAL VOL(NREG)
      REAL, TARGET :: SUNKNO(NUN,NGRP),FUNKNO(NUN,NGRP)
      LOGICAL LEXAC,REBFLG
*----
*  LOCAL VARIABLES
*----
      PARAMETER (IUNOUT=6,NSTATE=40)
      INTEGER IPAR(NSTATE)
      LOGICAL LBIHET
      CHARACTER CNFDIR(0:3)*9
      SAVE CNFDIR
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NGIND
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: INCONV
      REAL, ALLOCATABLE, DIMENSION(:) :: FGAR
      REAL, POINTER, DIMENSION(:,:) :: SUNKN,SUNKNO2,FUNKNO2
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: JPSYS,KPSYS
*----
*  DATA STATEMENTS
*----
      DATA CNFDIR /'F L U X  ','C U R - X','C U R - Y','C U R - Z'/
*----
*  GATHER INITIAL FLUXES AND SOURCES FOR PARALLEL PROCESSING
*----
      NGEFF=0
      JJ=0
      IIG=0
      DO IG=1,NGRP
         IF(NPSYS(IG).NE.0) THEN
            NGEFF=NGEFF+1
            IIG=IG
         ENDIF
      ENDDO
      IF(NGEFF.EQ.0) RETURN
      ALLOCATE(NGIND(NGEFF),JPSYS(NGEFF))
      IF(NGEFF.EQ.1) THEN
         NGIND(1)=IIG
         SUNKNO2=>SUNKNO(1:NUN,IIG:IIG)
         FUNKNO2=>FUNKNO(1:NUN,IIG:IIG)
         JPSYS(1)=LCMGIL(IPSYS,NPSYS(IIG))
      ELSE IF(NGEFF.EQ.NGRP) THEN
         SUNKNO2=>SUNKNO
         FUNKNO2=>FUNKNO
         DO IG=1,NGRP
            NGIND(IG)=IG
            JPSYS(IG)=LCMGIL(IPSYS,NPSYS(IG))
         ENDDO
      ELSE
         ALLOCATE(SUNKNO2(NUN,NGEFF),FUNKNO2(NUN,NGEFF))
         JJ=0
         DO IG=1,NGRP
            IF(NPSYS(IG).NE.0) THEN
               JJ=JJ+1
               NGIND(JJ)=IG
               SUNKNO2(:NUN,JJ)=SUNKNO(:NUN,IG)
               FUNKNO2(:NUN,JJ)=FUNKNO(:NUN,IG)
               JPSYS(JJ)=LCMGIL(IPSYS,NPSYS(IG))
            ENDIF
         ENDDO
      ENDIF
*
      IF(IMPX.GT.3) THEN
        WRITE(IUNOUT,'(//11H DOORFV: **,A,3H** ,A72)') CDOOR,TITR
        ALLOCATE(FGAR(NREG))
        DO II=1,NGEFF
          CALL XDRSET(FGAR,NREG,0.0)
          DO I=1,NREG
            IF(KEYFLX(I).NE.0) FGAR(I)=SUNKNO(KEYFLX(I),II)
          ENDDO
          WRITE(IUNOUT,'(/33H N E U T R O N    S O U R C E S (,I5,
     1    3H ):,4X,A9)') NGIND(II),CNFDIR(IDIR)
          WRITE(IUNOUT,'(1P,6(5X,E15.7))') (FGAR(I),I=1,NREG)
        ENDDO
        DEALLOCATE(FGAR)
      ENDIF
*---
* RECOVER FLUXES FROM A PREVIOUS SELF-SHIELDING CALCULATION IF AVAILABLE
*---
      DO JJ=1,NGEFF
         CALL LCMLEN(JPSYS(JJ),'FUNKNO$USS',ILENG,ITYLCM)
         IF(ILENG.EQ.NUN) THEN
            CALL LCMGET(JPSYS(JJ),'FUNKNO$USS',FUNKNO2(1,JJ))
         ENDIF
      ENDDO
*----
*  RECOVER STATE VECTOR INFORMATION
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',IPAR)
      NLF=IPAR(14)
      LBIHET=IPAR(40).NE.0
      IF(NREG.NE.IPAR(1)) CALL XABORT('DOORFV: INVALID NREG ON LCM.')
*----
*  DOUBLE HETEROGENEITY TREATMENT. REDEFINE THE SOURCE.
*----
      NMATG=0
      NREGG=0
      IF(LBIHET) THEN
         ALLOCATE(SUNKN(NUN,NGEFF))
         NMATG=NMAT
         NREGG=NREG
         DO II=1,NGEFF
           SUNKN(:NUN,II)=SUNKNO2(:NUN,II)
           IF(CDOOR.EQ.'MCCG') THEN
             CALL DOORFB2(JPSYS(II),IPTRK,IMPX,NMAT,NREG,NUN,KEYFLX,
     1       NMAT2,NREG2,SUNKNO(1,NGIND(II)))
           ELSE
             CALL DOORFB2(JPSYS(II),IPTRK,IMPX,NMAT,NREG,NUN,KEYFLX,
     1       NMAT2,NREG2,SUNKNO2(1,II))
           ENDIF
         ENDDO
         NMAT=NMAT2
         NREG=NREG2
      ELSE
         ALLOCATE(SUNKN(1,1))
      ENDIF
*---
*  RECOVER POINTERS TO EACH GROUP PROPERTIES
*---
      ALLOCATE(KPSYS(NGEFF))
      DO II=1,NGEFF
         IF(LBIHET) THEN
            KPSYS(II)=LCMGID(JPSYS(II),'BIHET')
         ELSE
            KPSYS(II)=JPSYS(II)
         ENDIF
      ENDDO
*----
*  COMPUTE NEW FLUXES
*----
      IF(IPHASE.EQ.1) THEN
*     USE A NATIVE DOOR
         IF(CDOOR.EQ.'SYBIL') THEN
            CALL SYBILF(KPSYS,IPTRK,IFTRAK,IMPX,NGEFF,NGIND,IDIR,NREG,
     >                  NUN,MAT,VOL,FUNKNO2,SUNKNO2,TITR)
         ELSE IF(CDOOR.EQ.'BIVAC') THEN
            IF(NLF.EQ.0) THEN
               CALL BIVAF(KPSYS,IPTRK,IFTRAK,IMPX,NGEFF,NGIND,IDIR,NREG,
     >                    NUN,MAT,VOL,KEYFLX,FUNKNO2,SUNKNO2,TITR)
            ELSE
               IF(IDIR.NE.0) CALL XABORT('DOORFV: EXPECTING IDIR=0')
               IF(IFTRAK.NE.0) CALL XABORT('DOORFV: EXPECTING IFTRAK=0')
               ALLOCATE(INCONV(NGEFF))
               INCONV(:NGEFF)=.TRUE.
               CALL PNFLV(KPSYS,INCONV,NGIND,IPTRK,IMPX,NGEFF,NREG,NMAT,
     >                    NUN,MAT,VOL,KEYFLX,FUNKNO2,SUNKNO2)
               DEALLOCATE(INCONV)
            ENDIF
         ELSE IF(CDOOR.EQ.'TRIVAC') THEN
            IF(IDIR.NE.0) CALL XABORT('DOORFV: EXPECTING IDIR=0')
            IF(IFTRAK.NE.0) CALL XABORT('DOORFV: EXPECTING IFTRAK=0')
            ALLOCATE(INCONV(NGEFF))
            INCONV(:NGEFF)=.TRUE.
            CALL TRIFLV(KPSYS,INCONV,NGIND,IPTRK,IMPX,NGEFF,NREG,NUN,
     >                  MAT,VOL,KEYFLX,FUNKNO2,SUNKNO2)
            DEALLOCATE(INCONV)
         ELSE IF(CDOOR.EQ.'SN') THEN
            CALL SNF(KPSYS,IPTRK,IFTRAK,IMPX,NGEFF,NGIND,IDIR,NREG,NMAT,
     >               NUN,MAT,VOL,KEYFLX,FUNKNO2,SUNKNO2,TITR)
         ELSE IF(CDOOR.EQ.'MCCG') THEN
            CALL MCCGF(KPSYS,IPTRK,IFTRAK,IPMACR,IMPX,NGRP,NGEFF,NGIND,
     >                 IDIR,NREG,NMAT,NUN,LEXAC,MAT,VOL,KEYFLX,FUNKNO2,
     >                 SUNKNO,TITR,REBFLG)
         ENDIF
      ELSE IF(IPHASE.EQ.2) THEN
         CALL TRFICF(KPSYS,IFTRAK,IMPX,NGEFF,NGIND,IDIR,NREG,NUN,MAT,
     >               VOL,KEYFLX,FUNKNO2,SUNKNO2,TITR)
      ENDIF
*----
* DOUBLE HETEROGENEITY TREATMENT
*----
      IF(LBIHET) THEN
         NMAT=NMATG
         NREG=NREGG
         DO II=1,NGEFF
           SUNKNO2(:NUN,II)=SUNKN(:NUN,II)
           CALL DOORFB3(JPSYS(II),IPTRK,IMPX,NMAT,NREG,NUN,KEYFLX,
     1     SUNKNO2(1,II),FUNKNO2(1,II))
         ENDDO
      ENDIF
      DEALLOCATE(SUNKN)
*
      IF(IMPX.GT.3) THEN
        ALLOCATE(FGAR(NREG))
        DO II=1,NGEFF
          CALL XDRSET(FGAR,NREG,0.0)
          DO I=1,NREG
            IF(KEYFLX(I).NE.0) FGAR(I)=FUNKNO(KEYFLX(I),II)
          ENDDO
          IF(IMPX.GT.4) THEN
            WRITE(IUNOUT,'(/31H U N K N O W N    F L U X E S (,I5,
     1      3H ):,6X,A9)') NGIND(II),CNFDIR(IDIR)
            WRITE(IUNOUT,300) (FUNKNO2(I,II),I=1,NUN)
          ELSE
            WRITE(IUNOUT,'(/31H N E U T R O N    F L U X E S (,I5,
     1      3H ):,6X,A9)') NGIND(II),CNFDIR(IDIR)
            WRITE(IUNOUT,300) (FGAR(I),I=1,NREG)
          ENDIF
        ENDDO
        DEALLOCATE(FGAR)
      ENDIF
*----
*  SCATTER NEW FLUXES
*----
      IF(CDOOR.EQ.'MCCG') THEN
         DO JJ=1,NGEFF
            FUNKNO(:NUN,NGIND(JJ))=FUNKNO2(:NUN,JJ)
         ENDDO
      ELSE
         DO JJ=1,NGEFF
            IG=NGIND(JJ)
            SUNKNO(:NUN,IG)=SUNKNO2(:NUN,JJ)
            FUNKNO(:NUN,IG)=FUNKNO2(:NUN,JJ)
         ENDDO
      ENDIF
      IF((NGEFF.GT.1).AND.(NGEFF.LT.NGRP)) DEALLOCATE(FUNKNO2,SUNKNO2)
      DEALLOCATE(KPSYS,JPSYS,NGIND)
      RETURN
 300  FORMAT(1P,6(5X,E15.7))
      END
