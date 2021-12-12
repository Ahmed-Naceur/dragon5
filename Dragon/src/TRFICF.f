*DECK TRFICF
      SUBROUTINE TRFICF(IPSYS,NPSYS,IPTRK,IFTRAK,IPRNTF,NGRP,IDIR,
     >                  NREGIO,NBMIX,NUNKNO,MATCOD,VOLUME,KEYFLX,
     >                  FUNKNO,SUNKNO,TITRE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solve N-group transport equation for fluxes using the scattering
* modified collision probability matrix.
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
* IPSYS   pointer to the pij matrices (L_PIJ signature). IPSYS is
*         a list of directories.
* NPSYS   index array pointing to the IPSYS list component corresponding
*         to each energy group. Set to zero if a group is not to be
*         processed. Usually, NPSYS(I)=I.
* IPTRK   pointer to the tracking LCM object.
* IFTRAK  not used.
* IPRNTF  print selection for flux modules.
* NGRP    number of energy groups.
* IDIR    directional collision probability flag:
*         =0 for pij or wij;
*         =k for pijk or wijk k=1,2,3.
* NREGIO  number of regions considered.
* NBMIX   number of mixtures.
* NUNKNO  number of unknown in the system.
* MATCOD  mixture code in region.
* VOLUME  volume of region.
* KEYFLX  flux elements in unknown system.
* SUNKNO  source for system of unknown.
* TITRE   title.
*
*Parameters: input/output
* FUNKNO  unknown vector solved for.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPSYS,IPTRK
      CHARACTER   TITRE*72
      INTEGER     NGRP,NPSYS(NGRP),IFTRAK,IPRNTF,IDIR,NREGIO,NBMIX,
     >            NUNKNO,MATCOD(NREGIO),KEYFLX(NREGIO)
      REAL        VOLUME(NREGIO),FUNKNO(NUNKNO,NGRP),SUNKNO(NUNKNO,NGRP)
*----
*  LOCAL VARIABLES
*----
      PARAMETER  (IUNOUT=6,NSTATE=40)
      TYPE(C_PTR) JPSYS
      CHARACTER   CNS(0:3)*1,CNFDIR(0:3)*9,NAMLCM*12,NAMMY*12
      INTEGER     ISTATE(NSTATE),ILCMLN
      LOGICAL     EMPTY,LCM,LBIHET
      SAVE        CNS,CNFDIR
*----
*  ALLOCATABLE ARRAYS
*----
      TYPE(C_PTR) CPMAT_PTR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INGIND
      REAL, ALLOCATABLE, DIMENSION(:) :: FGAR
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SUNKN
      REAL, POINTER, DIMENSION(:) :: CPMAT
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: KPSYS
*----
*  DATA STATEMENTS
*----
      DATA   CNS  /'-','1','2','3'/
      DATA CNFDIR /'F L U X  ','C U R - X','C U R - Y','C U R - Z'/
*----
*  RECOVER TRAFIC SPECIFIC PARAMETERS
*----
      CALL LCMINF(IPSYS,NAMLCM,NAMMY,EMPTY,ILONG,LCM)
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      LBIHET=ISTATE(40).NE.0
      IF(IFTRAK.LT.0) CALL XABORT('TRFICF: EXPECTING IFTRAK>=0')
      IF(MATCOD(1).LT.0) CALL XABORT('TRFICF: EXPECTING MATCOD(1)>=0')
      IF(VOLUME(1).LT.0.0) CALL XABORT('TRFICF: EXPECTING VOLUME(1)>=0')
*
      IF(IPRNTF.GT.3) THEN
        WRITE(IUNOUT,'(//9H TRFICF: ,A72)') TITRE
        ALLOCATE(FGAR(NREGIO))
        DO IG=1,NGRP
          IF(NPSYS(IG).EQ.0) CYCLE
          CALL XDRSET(FGAR,NREGIO,0.0)
          DO I=1,NREGIO
            IF(KEYFLX(I).NE.0) FGAR(I)=SUNKNO(KEYFLX(I),IG)
          ENDDO
          WRITE(IUNOUT,'(/33H N E U T R O N    S O U R C E S (,I5,
     1    3H ):,4X,A9)') IG,CNFDIR(IDIR)
          WRITE(IUNOUT,'(1P,6(5X,E15.7))') (FGAR(I),I=1,NREGIO)
        ENDDO
        DEALLOCATE(FGAR)
      ENDIF
*----
* DOUBLE HETEROGENEITY TREATMENT
*----
      NBMIXG=0
      NREGAR=0
      IF(LBIHET) THEN
         ALLOCATE(SUNKN(NUNKNO,NGRP))
         NBMIXG=NBMIX
         NREGAR=NREGIO
         DO IG=1,NGRP
           IOFSET=NPSYS(IG)
           IF(IOFSET.NE.0) THEN
             JPSYS=LCMGIL(IPSYS,IOFSET)
             SUNKN(:NUNKNO,IG)=SUNKNO(:NUNKNO,IG)
             CALL DOORFB2(JPSYS,IPTRK,IPRNTF,NBMIX,NREGIO,NUNKNO,KEYFLX,
     1       NBMIX2,NREG2,SUNKNO(1,IG))
           ENDIF
         ENDDO
         NBMIX=NBMIX2
         NREGIO=NREG2
      ENDIF
*---
* DETERMINE THE NUMBER OF GROUPS TO BE PROCESSED
*---
      NGEFF=0
      DO IG=1,NGRP
         IOFSET=NPSYS(IG)
         IF(IOFSET.NE.0) NGEFF=NGEFF+1
      ENDDO
      IF(NGEFF.EQ.0) GO TO 100
*---
* RECOVER POINTERS TO EACH GROUP PROPERTIES
* CREATE AN INDEX FOR THE GROUPS TO BE PROCESSED
*---
      ALLOCATE(INGIND(NGEFF),KPSYS(NGEFF))
      II=1
      DO IG=1,NGRP
         IOFSET=NPSYS(IG)
         IF(IOFSET.NE.0) THEN
            INGIND(II)=IG
            IF(LBIHET) THEN
               JPSYS=LCMGIL(IPSYS,IOFSET)
               KPSYS(II)=LCMGID(JPSYS,'BIHET')
            ELSE
               KPSYS(II)=LCMGIL(IPSYS,IOFSET)
            ENDIF
            II=II+1
         ENDIF
      ENDDO
*----
*  MAIN LOOP OVER ENERGY GROUPS.
*----
      IF(.NOT.LCM) ALLOCATE(CPMAT(NREGIO*NREGIO))
      DO 90 II=1,NGEFF
      JPSYS=KPSYS(II)
      IG=INGIND(II)
      IOFSET=NPSYS(IG)
      IF(IPRNTF.GT.1) WRITE(IUNOUT,'(/25H TRFICF: PROCESSING GROUP,I5,
     1 6H WITH ,A,1H.)') IG,'TRAFIC'
*----
*  READ SCATTERING MODIFIED COLLISION PROBABILITIES
*----
      CALL LCMLEN(JPSYS,'DRAGON'//CNS(IDIR)//'PCSCT',ILCMLN,ITYLCM)
      IF((ILCMLN.GT.0).AND.LCM) THEN
        CALL LCMGPD(JPSYS,'DRAGON'//CNS(IDIR)//'PCSCT',CPMAT_PTR)
        CALL C_F_POINTER(CPMAT_PTR,CPMAT,(/ NREGIO*NREGIO /))
      ELSE IF(ILCMLN.GT.0) THEN
        CALL LCMGET(JPSYS,'DRAGON'//CNS(IDIR)//'PCSCT',CPMAT)
      ELSE
        CALL XABORT('TRFICF: RECORD DRAGON'//CNS(IDIR)//
     >              'PCSCT ABSENT FROM LCM')
      ENDIF
*----
*  SOLVE TRANSPORT EQUATION
*----
      JCPMAT=0
      DO 30 I=1,NREGIO
        FUNKNO(KEYFLX(I),IG)=0.0
   30 CONTINUE
      DO 50 I=1,NREGIO
        IPOS=KEYFLX(I)
        DO 40 J=1,NREGIO
          JPOS=KEYFLX(J)
          JCPMAT=JCPMAT+1
          FUNKNO(JPOS,IG)=FUNKNO(JPOS,IG)+SUNKNO(IPOS,IG)*CPMAT(JCPMAT)
   40   CONTINUE
   50 CONTINUE
*----
* END OF LOOP OVER ENERGY GROUPS
*----
   90 CONTINUE
      DEALLOCATE(KPSYS,INGIND)
*----
* DOUBLE HETEROGENEITY TREATMENT
*----
  100 IF(LBIHET) THEN
         NBMIX=NBMIXG
         NREGIO=NREGAR
         DO IG=1,NGRP
           IOFSET=NPSYS(IG)
           IF(IOFSET.NE.0) THEN
             SUNKNO(:NUNKNO,IG)=SUNKN(:NUNKNO,IG)
             JPSYS=LCMGIL(IPSYS,IOFSET)
             CALL DOORFB3(JPSYS,IPTRK,IPRNTF,NBMIX,NREGIO,NUNKNO,KEYFLX,
     1       SUNKNO(1,IG),FUNKNO(1,IG))
           ENDIF
         ENDDO
         DEALLOCATE(SUNKN)
      ENDIF
*
      IF(IPRNTF.GT.3) THEN
        ALLOCATE(FGAR(NREGIO))
        DO IG=1,NGRP
          IF(NPSYS(IG).EQ.0) CYCLE
          CALL XDRSET(FGAR,NREGIO,0.0)
          DO I=1,NREGIO
            IF(KEYFLX(I).NE.0) FGAR(I)=FUNKNO(KEYFLX(I),IG)
          ENDDO
          WRITE(IUNOUT,'(//31H N E U T R O N    F L U X E S (,I5,
     1    3H ):,6X,A9)') IG,CNFDIR(IDIR)
          WRITE(IUNOUT,'(1P,6(5X,E15.7))') (FGAR(I),I=1,NREGIO)
        ENDDO
        DEALLOCATE(FGAR)
      ENDIF
      RETURN
      END
