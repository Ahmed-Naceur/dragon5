*DECK LIBADD
      SUBROUTINE LIBADD (IPLIB,NBISO,MASKI,IMPX,NGRO,NL,ITRANC,ISONAM,
     1 IPISO,NIR,GIR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Add transport correction and Goldstein-Cohen data to a /microlib/
* directory.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPLIB   pointer to the lattice microscopic cross section library
*         (L_LIBRARY signature).
* NBISO   number of isotopes present in the calculation domain.
* MASKI   isotopic mask. Isotope with index I is processed if
*         MASKI(I)=.true.
* IMPX    print flag.
* NGRO    number of energy groups.
* NL      number of Legendre orders required in the calculation
*         NL=1 (for isotropic scattering) or higher.
* ITRANC  transport correction option (=0: no correction; =1: Apollo-
*         type; =2: recover TRANC record; =3: Wims-type; =4: leakage
*         correction alone).
* ISONAM  alias name of each isotope.
* IPISO   pointer array towards microlib isotopes.
* NIR     group index with an imposed IR slowing-down model (=0 for no
*         IR model).
* GIR     value of the imposed Goldstein-Cohen parameter for groups
*         with an IR model.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB,IPISO(NBISO)
      INTEGER NBISO,IMPX,NGRO,NL,ITRANC,ISONAM(3,NBISO),NIR(NBISO)
      LOGICAL MASKI(NBISO)
      REAL GIR(NBISO)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (IOUT=6)
      TYPE(C_PTR) KPLIB
      CHARACTER TEXT12*12
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: WORK,WR2,DELTA
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SCAT
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(WORK(NGRO),WR2(NGRO),SCAT(NGRO,NGRO),DELTA(NGRO+1))
*----
*  RECOVER THE ENERGY GRID.
*----
      CALL LCMLEN(IPLIB,'ENERGY',LENGT,ITYLCM)
      IF(LENGT.EQ.0) CALL XABORT('LIBADD: NO GROUP STRUCTURE AVAILABLE')
      CALL LCMGET(IPLIB,'ENERGY',DELTA)
      NGX=0
      DO 10 IGR=1,NGRO
      IF((NGX.EQ.0).AND.(DELTA(IGR+1).LT.4.0)) NGX=IGR-1
  10  CONTINUE
      DO 15 IGR=1,NGRO
      DELTA(IGR)=LOG(DELTA(IGR)/DELTA(IGR+1))
  15  CONTINUE
*
      DO 110 ISO=1,NBISO
      IF(MASKI(ISO)) THEN
        WRITE(TEXT12,'(3A4)') ISONAM(1,ISO),ISONAM(2,ISO),ISONAM(3,ISO)
        KPLIB=IPISO(ISO) ! set ISO-th isotope
        IF(.NOT.C_ASSOCIATED(KPLIB)) GO TO 110
*
*       REDIFINE THE GOLDSTEIN-COHEN PARAMETERS.
        IF(NIR(ISO).GT.0) THEN
           DO 20 IGR=1,MIN(NGRO,NIR(ISO)-1)
           WORK(IGR)=1.0
   20      CONTINUE
           DO 30 IGR=NIR(ISO),NGRO
           WORK(IGR)=GIR(ISO)
   30      CONTINUE
           CALL LCMPUT(KPLIB,'NGOLD',NGRO,2,WORK)
           IF(IMPX.GT.1) THEN
              IF(GIR(ISO).EQ.-998.0) THEN
                 WRITE(IOUT,210) TEXT12,'PT',NIR(ISO)
              ELSE IF(GIR(ISO).EQ.-999.0) THEN
                 WRITE(IOUT,210) TEXT12,'PTSL',NIR(ISO)
              ELSE IF(GIR(ISO).EQ.-1000.0) THEN
                 WRITE(IOUT,210) TEXT12,'PTMC',NIR(ISO)
              ELSE
                 WRITE(IOUT,200) TEXT12,GIR(ISO),NIR(ISO)
              ENDIF
           ENDIF
        ENDIF
*
*       COMPUTE OR RECOVER THE TRANSPORT CORRECTION.
        IF(ITRANC.EQ.2) THEN
*          RECOVER THE TRANSPORT CORRECTION FROM THE LIBRARY.
           CALL LCMLEN(KPLIB,'TRANC',ILENG,ITYLCM)
           IF(ILENG.EQ.0) THEN
              CALL XDRSET(WORK,NGRO,0.0)
              CALL LCMPUT(KPLIB,'TRANC',NGRO,2,WORK)
           ENDIF
        ELSE IF(ITRANC.NE.0) THEN
           CALL XDRSET(WORK,NGRO,0.0)
           CALL LCMLEN(KPLIB,'NTOT1',ILENG,ITYLCM)
           IF(ILENG.NE.0) THEN
*             LEAKAGE CORRECTION.
              CALL LCMGET(KPLIB,'NTOT1',WORK)
              CALL LCMGET(KPLIB,'NTOT0',WR2)
              DO 40 IG1=1,NGRO
              WORK(IG1)=WR2(IG1)-WORK(IG1)
   40         CONTINUE
           ENDIF
           IF((NL.GE.2).AND.(ITRANC.NE.4)) THEN
              CALL LCMLEN(KPLIB,'SCAT-SAVED',ILENG,ITYLCM)
              IF(ILENG.EQ.0) THEN
                 CALL XABORT('LIBADD: NO SCAT-SAVED RECORD FOR ISOTO'//
     1           'PE '//TEXT12//'.')
              ENDIF
              CALL XDRLGS(KPLIB,-1,0,1,1,1,NGRO,WR2,SCAT,ITY)
              IF(ITRANC.EQ.1) THEN
*                APOLLO-TYPE TRANSPORT CORRECTION. USE THE MICRO-
*                REVERSIBILITY PRINCIPLE AT ALL ENERGIES.
                 DO 50 IG1=1,NGRO
                 WORK(IG1)=WORK(IG1)+WR2(IG1)
   50            CONTINUE
              ELSE IF(ITRANC.EQ.3) THEN
*                WIMS-TYPE TRANSPORT CORRECTION. USE THE MICRO-
*                REVERSIBILITY PRINCIPLE BELOW 4 EV AND A 1/E SPECTRUM
*                ABOVE.
                 DO 65 IG1=1,MIN(NGRO,NGX)
                 DO 60 IG2=1,NGRO
                 WORK(IG1)=WORK(IG1)+SCAT(IG1,IG2)*DELTA(IG2)/DELTA(IG1)
   60            CONTINUE
   65            CONTINUE
                 DO 70 IG1=NGX+1,NGRO
                 WORK(IG1)=WORK(IG1)+WR2(IG1)
   70            CONTINUE
              ELSE
                 CALL XABORT('LIBADD: UNKNOWN TYPE OF CORRECTION.')
              ENDIF
           ENDIF
*          ***CAUTION*** 'TRANC' CONTAINS BOTH TRANSPORT AND LEAKAGE
*          CORRECTIONS.
           CALL LCMPUT(KPLIB,'TRANC',NGRO,2,WORK)
        ENDIF
      ENDIF
  110 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(DELTA,SCAT,WR2,WORK)
      RETURN
*
  200 FORMAT(/51H LIBADD: THE GOLDSTEIN-COHEN PARAMETER OF ISOTOPE ',
     1 A12,12H' WAS SET TO,F5.2,33H FOR GROUPS WITH INDEX GREATER OR,
     2 9H EQUAL TO,I4,1H.)
  210 FORMAT(/18H LIBADD: ISOTOPE ',A12,20H' IS PROCESSED WITH ,A,
     1 48H METHOD IN GROUPS WITH INDEX GREATER OR EQUAL TO,I4,1H.)
      END
