*DECK TRFICF
      SUBROUTINE TRFICF(KPSYS,IFTRAK,IPRNTF,NGEFF,NGIND,IDIR,NREGIO,
     >                  NUNKNO,MATCOD,VOLUME,KEYFLX,FUNKNO,SUNKNO,
     >                  TITRE)
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
* KPSYS   pointer to the pij matrices (L_PIJ signature). KPSYS is
*         an array of directories.
* IFTRAK  not used.
* IPRNTF  print selection for flux modules.
* NGEFF   number of energy groups processed in parallel.
* NGIND   energy group indices assign to the NGEFF set.
* IDIR    directional collision probability flag:
*         =0 for pij or wij;
*         =k for pijk or wijk k=1,2,3.
* NREGIO  number of regions considered.
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
      TYPE(C_PTR) KPSYS(NGEFF)
      CHARACTER   TITRE*72
      INTEGER     NGEFF,NGIND(NGEFF),IFTRAK,IPRNTF,IDIR,NREGIO,NUNKNO,
     >            MATCOD(NREGIO),KEYFLX(NREGIO)
      REAL        VOLUME(NREGIO),FUNKNO(NUNKNO,NGEFF),
     >            SUNKNO(NUNKNO,NGEFF)
*----
*  LOCAL VARIABLES
*----
      PARAMETER  (IUNOUT=6)
      CHARACTER   CNS(0:3)*1,NAMLCM*12,NAMMY*12
      INTEGER     ILCMLN
      LOGICAL     EMPTY,LCM
      SAVE        CNS
*----
*  ALLOCATABLE ARRAYS
*----
      TYPE(C_PTR) CPMAT_PTR
      REAL, POINTER, DIMENSION(:) :: CPMAT
*----
*  DATA STATEMENTS
*----
      DATA   CNS  /'-','1','2','3'/
*----
*  RECOVER TRAFIC SPECIFIC PARAMETERS
*----
      IF(IPRNTF.GT.2) WRITE(IUNOUT,'(//9H TRFICF: ,A72)') TITRE
      CALL LCMINF(KPSYS(1),NAMLCM,NAMMY,EMPTY,ILONG,LCM)
      IF(IFTRAK.LT.0) CALL XABORT('TRFICF: EXPECTING IFTRAK>=0')
      IF(MATCOD(1).LT.0) CALL XABORT('TRFICF: EXPECTING MATCOD(1)>=0')
      IF(VOLUME(1).LT.0.0) CALL XABORT('TRFICF: EXPECTING VOLUME(1)>=0')
*----
*  MAIN LOOP OVER ENERGY GROUPS.
*----
      IF(.NOT.LCM) ALLOCATE(CPMAT(NREGIO*NREGIO))
      DO 60 II=1,NGEFF
      IF(IPRNTF.GT.1) WRITE(IUNOUT,'(/25H TRFICF: PROCESSING GROUP,I5,
     1 6H WITH ,A,1H.)') NGIND(II),'TRAFIC'
*----
*  READ SCATTERING MODIFIED COLLISION PROBABILITIES
*----
      CALL LCMLEN(KPSYS(II),'DRAGON'//CNS(IDIR)//'PCSCT',ILCMLN,ITYLCM)
      IF((ILCMLN.GT.0).AND.LCM) THEN
        CALL LCMGPD(KPSYS(II),'DRAGON'//CNS(IDIR)//'PCSCT',CPMAT_PTR)
        CALL C_F_POINTER(CPMAT_PTR,CPMAT,(/ NREGIO*NREGIO /))
      ELSE IF(ILCMLN.GT.0) THEN
        CALL LCMGET(KPSYS(II),'DRAGON'//CNS(IDIR)//'PCSCT',CPMAT)
      ELSE
        CALL XABORT('TRFICF: RECORD DRAGON'//CNS(IDIR)//
     >              'PCSCT ABSENT FROM LCM')
      ENDIF
*----
*  SOLVE TRANSPORT EQUATION
*----
      JCPMAT=0
      DO 30 I=1,NREGIO
        FUNKNO(KEYFLX(I),II)=0.0
   30 CONTINUE
      DO 50 I=1,NREGIO
        IPOS=KEYFLX(I)
        DO 40 J=1,NREGIO
          JPOS=KEYFLX(J)
          JCPMAT=JCPMAT+1
          FUNKNO(JPOS,II)=FUNKNO(JPOS,II)+SUNKNO(IPOS,II)*CPMAT(JCPMAT)
   40   CONTINUE
   50 CONTINUE
*----
* END OF LOOP OVER ENERGY GROUPS
*----
   60 CONTINUE
      RETURN
      END
