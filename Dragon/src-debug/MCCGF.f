*DECK MCCGF
      SUBROUTINE MCCGF(IPSYS,NPSYS,IPTRK,IFTRAK,IPMACR,IMPX,NGRP,IDIR,
     1                 NBREG,NBMIX,NUNKNO,LEXAC,MAT,VOL,KEYFLX,FUNKNO,
     2                 SUNKNO,TITR,REBFLG)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solve N-group transport equation for fluxes using the method of
* characteristics (vectorial version).
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert and R. Le Tellier
*
*Parameters: input/output
* IPSYS   pointer to the assembly LCM object (L_PIJ signature). IPSYS is
*         a list of directories.
* NPSYS   index array pointing to the IPSYS list component corresponding
*         to each energy group. Set to zero if a group is not to be
*         processed. Usually, NPSYS(I)=I.
* IPTRK   pointer to the tracking (L_TRACK signature).
* IPMACR  pointer to the macrolib LCM object.
* IFTRAK  tracking file unit number.
* IMPX    print flag (equal to zero for no print).
* NGRP    number of energy groups.
* IDIR    direction of fundamental current for TIBERE with MoC 
*         (=0,1,2,3). 
* NBREG   total number of volumes for which specific values of the
*         neutron flux and reactions rates are required.
* NBMIX   number of mixtures (NBMIX=max(MAT(i))).
* NUNKNO  total number of unknowns in vectors SUNKNO and FUNKNO.
* LEXAC   type of exponential function calculation (=.false. to compute
*         exponential functions using tables).
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* KEYFLX  position of flux elements in FUNKNO vector.
* FUNKNO  unknown vector.
* SUNKNO  input source vector.
* TITR    title.
* REBFLG  ACA or SCR rebalancing flag.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPSYS,IPTRK,IPMACR
      INTEGER NGRP,NPSYS(NGRP),IFTRAK,IMPX,IDIR,NBREG,NBMIX,NUNKNO,
     1 MAT(NBREG),KEYFLX(NBREG)
      REAL VOL(NBREG),FUNKNO(NUNKNO,NGRP),SUNKNO(NUNKNO,NGRP)
      CHARACTER TITR*72
      LOGICAL LEXAC,REBFLG
*----
*  LOCAL VARIABLES
*----
      PARAMETER  (IUNOUT=6,NSTATE=40,MXNMU=64)
      TYPE(C_PTR) JPSYS
      CHARACTER   TEXT4*4,TOPT*72
      INTEGER     JPAR(NSTATE),PACA,STIS,TRTY,IGB(8)
      REAL        ZREAL(4),HDD,DELU,FACSYM
      LOGICAL     CYCLIC,LVOID,LEXF,LFORW,LPRISM,LBIHET
      EXTERNAL    MOCFFI,MOCFFA,MOCFFIS,MOCFFAS,MOCFFIT,MOCFFAT,
     1            MCGFFI,MCGFFA,MCGFFIS,MCGFFAS,MCGFFIT,MCGFFAT
      EXTERNAL    MOCSCA,MOCDDF,MOCSCE,MOCSCAS,MOCDDFS,MOCSCES,MOCSCAT,
     1            MOCDDFT,MOCSCET,MOCFFAL,MOCSCEL,MOCSCAL,MOCDDFL,
     2            MCGSCA,MCGDDF,MCGSCE,MCGSCAS,MCGDDFS,MCGSCES,MCGSCAT,
     3            MCGDDFT,MCGSCET,MCGFFAL,MCGSCEL,MCGSCAL,MCGDDFL
      INTEGER, TARGET, SAVE, DIMENSION(1) :: IDUMMY
*----
*  ALLOCATABLE ARRAYS
*----
      TYPE(C_PTR) WZMU_PTR,ZMU_PTR,V_PTR,NZON_PTR,KEY_PTR,KEYCUR_PTR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: INGIND,ITST,MATALB
      REAL, ALLOCATABLE, DIMENSION(:) :: SUNKN,SIGAL,REPS,EPS,CPO
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: CAZ0,CAZ1,CAZ2
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: INCONV
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: KPSYS
      INTEGER, POINTER, DIMENSION(:) :: NZON,KEY,KEYCUR
      REAL, POINTER, DIMENSION(:) :: WZMU,ZMU,V
*
      IF(MAT(1).LT.0) CALL XABORT('MCCGF: EXPECTING MAT(1)>=0')
      IF(VOL(1).LT.0.0) CALL XABORT('MCCGF: EXPECTING VOL(1)>=0')
      IF(IMPX.GT.3) WRITE(IUNOUT,'(//8H MCCGF: ,A72/)') TITR
*----
* DOUBLE HETEROGENEITY TREATMENT
*----
      NBMIXG=0
      NREGAR=0
      CALL LCMGET(IPTRK,'STATE-VECTOR',JPAR)
      LBIHET=JPAR(40).NE.0
      IF(LBIHET) THEN
         ALLOCATE(SUNKN(NUNKNO*NGRP))
         NBMIXG=NBMIX
         NREGAR=NBREG
         DO IG=1,NGRP
           IOFSET=NPSYS(IG)
           IF(IOFSET.NE.0) THEN
             JPSYS=LCMGIL(IPSYS,IOFSET)
             DO I=1,NUNKNO
               SUNKN((IG-1)*NUNKNO+I)=SUNKNO(I,IG)
             ENDDO
             CALL DOORFB2(JPSYS,IPTRK,IMPX,NBMIX,NBREG,NUNKNO,KEYFLX,
     1       NBMIX2,NBREG2,SUNKNO(1,IG))
           ENDIF
         ENDDO
         NBMIX=NBMIX2
         NBREG=NBREG2
      ENDIF
*---
* DETERMINE THE NUMBER OF GROUPS TO BE PROCESSED
* RECOVER FLUXES FROM A PREVIOUS SELF-SHIELDING CALCULATION IF AVAILABLE
*---
      NGEFF=0
      DO IG=1,NGRP
         IOFSET=NPSYS(IG)
         IF(IOFSET.NE.0) THEN
            NGEFF=NGEFF+1
            JPSYS=LCMGIL(IPSYS,IOFSET)
            CALL LCMLEN(JPSYS,'FUNKNO$USS',ILENG,ITYLCM)
            IF(ILENG.EQ.NUNKNO) THEN
               CALL LCMGET(JPSYS,'FUNKNO$USS',FUNKNO(1,IG))
            ENDIF
         ENDIF
      ENDDO
      IF(NGEFF.EQ.0) GO TO 40
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
*  RECOVER MCCG SPECIFIC PARAMETERS
*----
*     check for cross-sections in SYS object
      JPSYS=KPSYS(1)
      CALL LCMLEN(JPSYS,'DRAGON-TXSC',ILENG,ITYLCM)
      IF(ILENG.NE.NBMIX+1) CALL XABORT('MCCGF: INVALID VALUE OF NBMIX.')
      IF(JPAR(4).GT.NBMIX) CALL XABORT('MCCGF: MIXTURE OVERFLOW.')
*     check for a tracking binary file
      IF(IFTRAK.LE.0) CALL XABORT('MCCGF: INVALID TRACKING FILE.')
*     recover state-vector information
      IF(JPAR(40).EQ.1) THEN
         CALL LCMSIX(IPTRK,'BIHET',1)
         CALL LCMGET(IPTRK,'PARAM',IGB)
         NREG=IGB(3)
         CALL LCMSIX(IPTRK,' ',2)
      ELSE
         NREG=JPAR(1)
      ENDIF
      NSOU=JPAR(5)
      NFI=NREG+NSOU
      IF(JPAR(2).GT.NUNKNO) 
     1 CALL XABORT('MCCGF: UNKNOWN VECTOR OVERFLOW.')
      NANI=JPAR(6)
      TRTY=JPAR(9)
      IF(TRTY.EQ.1) THEN
         CYCLIC=.TRUE.
         NLONG=NREG
      ELSE
         CYCLIC=.FALSE.
         NLONG=NFI
      ENDIF
      NZP=JPAR(39)
      LPRISM=(NZP.NE.0)
      CALL LCMGET(IPTRK,'MCCG-STATE',JPAR)
      NMU=JPAR(2)
      IF(NMU.GT.MXNMU)
     1 CALL XABORT('MCCGF: POLAR ANGLE QUADRATURE OVERFLOW') 
      NMAX=JPAR(5)
      MAXI=JPAR(13)
      STIS=JPAR(15)
      LC=JPAR(6)
      IAAC=JPAR(7)
      KRYL=JPAR(3)
      IDIFC=JPAR(4)
      ISCR=JPAR(8)
      LPS=JPAR(9)
      PACA=JPAR(10)
      LEXF=(JPAR(12).EQ.1)
      LFORW=(JPAR(18).EQ.0)
      NFUNL=JPAR(19)
      NLIN=JPAR(20)
*     to be coherent with the exponential function used for the Pjj calculation
      IF((LEXAC).AND.(.NOT.LEXF).AND.(STIS.EQ.1)) STIS=0 
      NPJJM=JPAR(16)
*     recover real parameters
      CALL LCMGET(IPTRK,'REAL-PARAM',ZREAL)
      EPSI=ZREAL(1)
      DELU=ZREAL(3)
      FACSYM=ZREAL(4)
!!! temporary
      HDD=ZREAL(2)
      IF(HDD.GT.0.0) THEN
         ISCH=0
      ELSEIF(LEXF) THEN
         ISCH=-1
      ELSE
         ISCH=1
      ENDIF
!!!
*----
* RECOVER TRACKING FILE INFORMATION
*----
      REWIND IFTRAK
      READ(IFTRAK) TEXT4,NCOMNT,NBTR,IFMT
      DO ICOM=1,NCOMNT
         READ(IFTRAK)
      ENDDO
      READ(IFTRAK) NDIM,ISPEC,N2REG,N2SOU,NALBG,NCOR,NANGL,MXSUB,MXSEG
      IF(NCOR.NE.1) 
     1 CALL XABORT('MCCGF: INVALID TRACKING FILE: NCOR.NE.1')
      ALLOCATE(MATALB(N2REG+N2SOU+1))
      READ(IFTRAK)
      READ(IFTRAK) (MATALB(JJ),JJ=1,N2REG+N2SOU+1)
      READ(IFTRAK)
      READ(IFTRAK)
      ALLOCATE(CAZ0(NANGL),CAZ1(NANGL),CAZ2(NANGL),CPO(NMU))
      IF(NDIM.EQ.2) THEN
         CALL LCMGET(IPTRK,'XMU$MCCG',CPO)
         READ(IFTRAK) (CAZ1(JJ),CAZ2(JJ),JJ=1,NANGL)
      ELSE ! NDIM.EQ.3
**        correction Sylvie Musongela, december 2019
         READ(IFTRAK) (CAZ1(JJ),CAZ2(JJ),CAZ0(JJ),JJ=1,NANGL)
         DO JJ=1,NANGL
            CAZ1(JJ)=CAZ1(JJ)/SQRT(1.0D0-CAZ0(JJ)*CAZ0(JJ))
            CAZ2(JJ)=CAZ2(JJ)/SQRT(1.0D0-CAZ0(JJ)*CAZ0(JJ))
         ENDDO
      ENDIF
*----
* RECOVER TRACKING TABLE INFORMATION
*----
*     recover polar quadrature
      CALL LCMGPD(IPTRK,'WZMU$MCCG',WZMU_PTR)
      CALL LCMGPD(IPTRK,'ZMU$MCCG',ZMU_PTR)
*     recover modified MATALB, VOLSUR and KEYFLX
      CALL LCMGPD(IPTRK,'V$MCCG',V_PTR)
      CALL LCMGPD(IPTRK,'NZON$MCCG',NZON_PTR)
      CALL LCMGPD(IPTRK,'KEYFLX$ANIS',KEY_PTR)
*     recover index for the currents in FUNKNO (non-cyclic case)
      IF(.NOT.CYCLIC) CALL LCMGPD(IPTRK,'KEYCUR$MCCG',KEYCUR_PTR)
*
      CALL C_F_POINTER(WZMU_PTR,WZMU,(/ NMU /))
      CALL C_F_POINTER(ZMU_PTR,ZMU,(/ NMU /))
      CALL C_F_POINTER(V_PTR,V,(/ NLONG /))
      CALL C_F_POINTER(NZON_PTR,NZON,(/ NLONG /))
      CALL C_F_POINTER(KEY_PTR,KEY,(/ NREG*NLIN*NFUNL /))
      IF(.NOT.CYCLIC) THEN
         CALL C_F_POINTER(KEYCUR_PTR,KEYCUR,(/ NLONG-NBREG /))
      ELSE
         KEYCUR=>IDUMMY
      ENDIF
*----
*  CONSTRUCT TOTAL CROSS SECTIONS ARRAY AND CHECK FOR ZERO CROSS SECTION
*----
      CALL LCMLEN(KPSYS(1),'ALBEDO',NALBP,ITYLCM)
      ALLOCATE(SIGAL((NBMIX+7)*NGEFF))
      CALL MCGSIG(IPTRK,NBMIX,NGEFF,NALBP,KPSYS,SIGAL,LVOID)
      IF((LVOID).AND.(STIS.EQ.-1)) THEN
         IF(IMPX.GT.0) 
     1     WRITE(IUNOUT,*) 'VOID EXISTS -> STIS SET TO 1 INSTEAD OF -1'
         STIS=1
      ENDIF
      ISCH=ISCH+10*STIS+100*(NLIN-1)
*----
*  PERFORM INNER ITERATIONS TO COMPUTE THE NEUTRON FLUX IN THE DIFFERENT
*  GROUPS
*----
      ALLOCATE(REPS(MAXI*NGEFF),EPS(NGEFF),ITST(NGEFF),INCONV(NGEFF))
      CALL XDLSET(INCONV,NGEFF,.TRUE.)
      LNCONV=NGEFF
*
      IF(IDIFC.EQ.1) THEN
*     ------------------------------------
*     ACA-Simplified Transport Calculation
*     ------------------------------------
         TOPT='ACA-SIMPLIFIED TRANSPORT OPERATOR'
         CALL MCGFLS(IMPX,IPTRK,IPMACR,NUNKNO,NFI,NBREG,NLONG,NBMIX,
     1        NGRP,NGEFF,LC,LFORW,PACA,NZON,KEY,KEYCUR,INGIND,KPSYS,
     2        INCONV,EPSI,MAXI,FUNKNO,SUNKNO)
*     ------------------------------------
      ELSE
      IF(CYCLIC) THEN
*     --------------------------------
*     Method of Cyclic Characteristics
*     --------------------------------
*********'Source Term Isolation' Strategy turned off 
         IF(ISCH.EQ.1) THEN
*          Step-Characteristics Scheme with Tabulated Exponentials
           TOPT='CYCLIC - STIS 0 - SC SCHEME - TABULATED EXP'
           CALL MCGFLX(MOCFFIS,MOCFFAS,MOCSCAS,MOCFFAL,CYCLIC,KPSYS,
     1          IMPX,IPTRK,IFTRAK,IPMACR,NDIM,NFI,NUNKNO,NLONG,NBREG,
     2          NSOU,NGRP,NGEFF,INGIND,NZON,MATALB,V,FUNKNO,SUNKNO,
     3          NBMIX,NANI,MAXI,IAAC,KRYL,ISCR,NMU,NANGL,NMAX,LC,EPSI,
     4          CAZ0,CAZ1,CAZ2,CPO,ZMU,WZMU,LFORW,PACA,NLIN,NFUNL,KEY,
     5          KEYCUR,SIGAL,LPS,REPS,EPS,ITST,INCONV,LNCONV,REBFLG,
     6          STIS,NPJJM,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
         ELSEIF(ISCH.EQ.0) THEN
*          Diamond-Differencing Scheme
           TOPT='CYCLIC - STIS 0 - DD0 SCHEME'
           CALL MCGFLX(MOCFFIS,MOCFFAS,MOCDDFS,MOCFFAL,CYCLIC,KPSYS,
     1          IMPX,IPTRK,IFTRAK,IPMACR,NDIM,NFI,NUNKNO,NLONG,NBREG,
     2          NSOU,NGRP,NGEFF,INGIND,NZON,MATALB,V,FUNKNO,SUNKNO,
     3          NBMIX,NANI,MAXI,IAAC,KRYL,ISCR,NMU,NANGL,NMAX,LC,EPSI,
     4          CAZ0,CAZ1,CAZ2,CPO,ZMU,WZMU,LFORW,PACA,NLIN,NFUNL,KEY,
     5          KEYCUR,SIGAL,LPS,REPS,EPS,ITST,INCONV,LNCONV,REBFLG,
     6          STIS,NPJJM,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
         ELSEIF(ISCH.EQ.-1) THEN
*          Step-Characteristics Scheme with Exact Exponentials
           TOPT='CYCLIC - STIS 0 - SC SCHEME - EXACT EXP'
           CALL MCGFLX(MOCFFIS,MOCFFAS,MOCSCES,MOCFFAL,CYCLIC,KPSYS,
     1          IMPX,IPTRK,IFTRAK,IPMACR,NDIM,NFI,NUNKNO,NLONG,NBREG,
     2          NSOU,NGRP,NGEFF,INGIND,NZON,MATALB,V,FUNKNO,SUNKNO,
     3          NBMIX,NANI,MAXI,IAAC,KRYL,ISCR,NMU,NANGL,NMAX,LC,EPSI,
     4          CAZ0,CAZ1,CAZ2,CPO,ZMU,WZMU,LFORW,PACA,NLIN,NFUNL,KEY,
     5          KEYCUR,SIGAL,LPS,REPS,EPS,ITST,INCONV,LNCONV,REBFLG,
     6          STIS,NPJJM,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
*********'Source Term Isolation' Strategy turned on
         ELSEIF(ISCH.EQ.11) THEN
*          Step-Characteristics Scheme with Tabulated Exponentials
           TOPT='CYCLIC - STIS 1 - SC SCHEME - TABULATED EXP'
           CALL MCGFLX(MOCFFI,MOCFFA,MOCSCA,MOCFFAL,CYCLIC,KPSYS,IMPX,
     1          IPTRK,IFTRAK,IPMACR,NDIM,NFI,NUNKNO,NLONG,NBREG,NSOU,
     2          NGRP,NGEFF,INGIND,NZON,MATALB,V,FUNKNO,SUNKNO,NBMIX,
     3          NANI,MAXI,IAAC,KRYL,ISCR,NMU,NANGL,NMAX,LC,EPSI,CAZ0,
     4          CAZ1,CAZ2,CPO,ZMU,WZMU,LFORW,PACA,NLIN,NFUNL,KEY,KEYCUR,
     5          SIGAL,LPS,REPS,EPS,ITST,INCONV,LNCONV,REBFLG,STIS,
     6          NPJJM,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
         ELSEIF(ISCH.EQ.10) THEN
*          Diamond-Differencing Scheme
           TOPT='CYCLIC - STIS 1 - DD0 SCHEME'
           CALL MCGFLX(MOCFFI,MOCFFA,MOCDDF,MOCFFAL,CYCLIC,KPSYS,IMPX,
     1          IPTRK,IFTRAK,IPMACR,NDIM,NFI,NUNKNO,NLONG,NBREG,NSOU,
     2          NGRP,NGEFF,INGIND,NZON,MATALB,V,FUNKNO,SUNKNO,NBMIX,
     3          NANI,MAXI,IAAC,KRYL,ISCR,NMU,NANGL,NMAX,LC,EPSI,CAZ0,
     4          CAZ1,CAZ2,CPO,ZMU,WZMU,LFORW,PACA,NLIN,NFUNL,KEY,KEYCUR,
     5          SIGAL,LPS,REPS,EPS,ITST,INCONV,LNCONV,REBFLG,STIS,
     6          NPJJM,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
         ELSEIF(ISCH.EQ.9) THEN
*          Step-Characteristics Scheme with Exact Exponentials
           TOPT='CYCLIC - STIS 1 - SC SCHEME - EXACT EXP'
           CALL MCGFLX(MOCFFI,MOCFFA,MOCSCE,MOCFFAL,CYCLIC,KPSYS,IMPX,
     1          IPTRK,IFTRAK,IPMACR,NDIM,NFI,NUNKNO,NLONG,NBREG,NSOU,
     2          NGRP,NGEFF,INGIND,NZON,MATALB,V,FUNKNO,SUNKNO,NBMIX,
     3          NANI,MAXI,IAAC,KRYL,ISCR,NMU,NANGL,NMAX,LC,EPSI,CAZ0,
     4          CAZ1,CAZ2,CPO,ZMU,WZMU,LFORW,PACA,NLIN,NFUNL,KEY,KEYCUR,
     5          SIGAL,LPS,REPS,EPS,ITST,INCONV,LNCONV,REBFLG,STIS,
     6          NPJJM,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
*********'MOCC/MCI' Iterative Strategy
         ELSEIF(ISCH.EQ.-9) THEN
*          Step-Characteristics Scheme with Tabulated Exponentials
           TOPT='CYCLIC - STIS -1 - SC SCHEME - TABULATED EXP'
           CALL MCGFLX(MOCFFIT,MOCFFAT,MOCSCAT,MOCFFAL,CYCLIC,KPSYS,
     1          IMPX,IPTRK,IFTRAK,IPMACR,NDIM,NFI,NUNKNO,NLONG,NBREG,
     2          NSOU,NGRP,NGEFF,INGIND,NZON,MATALB,V,FUNKNO,SUNKNO,
     3          NBMIX,NANI,MAXI,IAAC,KRYL,ISCR,NMU,NANGL,NMAX,LC,EPSI,
     4          CAZ0,CAZ1,CAZ2,CPO,ZMU,WZMU,LFORW,PACA,NLIN,NFUNL,KEY,
     5          KEYCUR,SIGAL,LPS,REPS,EPS,ITST,INCONV,LNCONV,REBFLG,
     6          STIS,NPJJM,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
         ELSEIF(ISCH.EQ.-10) THEN
*          Diamond-Differencing Scheme
           TOPT='CYCLIC - STIS -1 - DD0 SCHEME'
           CALL MCGFLX(MOCFFIT,MOCFFAT,MOCDDFT,MOCFFAL,CYCLIC,KPSYS,
     1          IMPX,IPTRK,IFTRAK,IPMACR,NDIM,NFI,NUNKNO,NLONG,NBREG,
     2          NSOU,NGRP,NGEFF,INGIND,NZON,MATALB,V,FUNKNO,SUNKNO,
     3          NBMIX,NANI,MAXI,IAAC,KRYL,ISCR,NMU,NANGL,NMAX,LC,EPSI,
     4          CAZ0,CAZ1,CAZ2,CPO,ZMU,WZMU,LFORW,PACA,NLIN,NFUNL,KEY,
     5          KEYCUR,SIGAL,LPS,REPS,EPS,ITST,INCONV,LNCONV,REBFLG,
     6          STIS,NPJJM,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
         ELSEIF(ISCH.EQ.-11) THEN
*          Step-Characteristics Scheme with Exact Exponentials
           TOPT='CYCLIC - STIS -1 - SC SCHEME - EXACT EXP'
           CALL MCGFLX(MOCFFIT,MOCFFAT,MOCSCET,MOCFFAL,CYCLIC,KPSYS,
     1          IMPX,IPTRK,IFTRAK,IPMACR,NDIM,NFI,NUNKNO,NLONG,NBREG,
     2          NSOU,NGRP,NGEFF,INGIND,NZON,MATALB,V,FUNKNO,SUNKNO,
     3          NBMIX,NANI,MAXI,IAAC,KRYL,ISCR,NMU,NANGL,NMAX,LC,EPSI,
     4          CAZ0,CAZ1,CAZ2,CPO,ZMU,WZMU,LFORW,PACA,NLIN,NFUNL,KEY,
     5          KEYCUR,SIGAL,LPS,REPS,EPS,ITST,INCONV,LNCONV,REBFLG,
     6          STIS,NPJJM,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
         ELSEIF(ISCH.EQ.199) THEN
*          Lin.-Disc.-Characteristics Scheme with Exact Exponentials
           TOPT='CYCLIC - STIS 0 - LDC SCHEME - EXACT EXP'
           CALL MCGFLX(MOCFFIT,MOCFFAT,MOCSCEL,MOCFFAL,CYCLIC,KPSYS,
     1          IMPX,IPTRK,IFTRAK,IPMACR,NDIM,NFI,NUNKNO,NLONG,NBREG,
     2          NSOU,NGRP,NGEFF,INGIND,NZON,MATALB,V,FUNKNO,SUNKNO,
     3          NBMIX,NANI,MAXI,IAAC,KRYL,ISCR,NMU,NANGL,NMAX,LC,EPSI,
     4          CAZ0,CAZ1,CAZ2,CPO,ZMU,WZMU,LFORW,PACA,NLIN,NFUNL,KEY,
     5          KEYCUR,SIGAL,LPS,REPS,EPS,ITST,INCONV,LNCONV,REBFLG,
     6          STIS,NPJJM,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
         ELSEIF(ISCH.EQ.200) THEN
*          Lin.-Disc.-Characteristics Scheme with Exact Exponentials
           TOPT='CYCLIC - STIS 0 - LDC SCHEME - DD1 SCHEME'
           CALL MCGFLX(MOCFFIT,MOCFFAT,MOCDDFL,MOCFFAL,CYCLIC,KPSYS,
     1          IMPX,IPTRK,IFTRAK,IPMACR,NDIM,NFI,NUNKNO,NLONG,NBREG,
     2          NSOU,NGRP,NGEFF,INGIND,NZON,MATALB,V,FUNKNO,SUNKNO,
     3          NBMIX,NANI,MAXI,IAAC,KRYL,ISCR,NMU,NANGL,NMAX,LC,EPSI,
     4          CAZ0,CAZ1,CAZ2,CPO,ZMU,WZMU,LFORW,PACA,NLIN,NFUNL,KEY,
     5          KEYCUR,SIGAL,LPS,REPS,EPS,ITST,INCONV,LNCONV,REBFLG,
     6          STIS,NPJJM,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
         ELSEIF(ISCH.EQ.201) THEN
*          Lin.-Disc.-Characteristics Scheme with Exact Exponentials
           TOPT='CYCLIC - STIS 0 - LDC SCHEME - TABULATED EXP'
           CALL MCGFLX(MOCFFIT,MOCFFAT,MOCSCAL,MOCFFAL,CYCLIC,KPSYS,
     1          IMPX,IPTRK,IFTRAK,IPMACR,NDIM,NFI,NUNKNO,NLONG,NBREG,
     2          NSOU,NGRP,NGEFF,INGIND,NZON,MATALB,V,FUNKNO,SUNKNO,
     3          NBMIX,NANI,MAXI,IAAC,KRYL,ISCR,NMU,NANGL,NMAX,LC,EPSI,
     4          CAZ0,CAZ1,CAZ2,CPO,ZMU,WZMU,LFORW,PACA,NLIN,NFUNL,KEY,
     5          KEYCUR,SIGAL,LPS,REPS,EPS,ITST,INCONV,LNCONV,REBFLG,
     6          STIS,NPJJM,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
         ELSE
           CALL XABORT('MCCGF: CYCLIC SCHEME NOT IMPLEMENTED')
         ENDIF   
*     --------------------------------
      ELSE
*     ------------------------------------
*     Method of Non-Cyclic Characteristics
*     ------------------------------------
*********'Source Term Isolation' Strategy turned off 
         IF(ISCH.EQ.1) THEN
*          Step-Characteristics Scheme with Tabulated Exponentials
           TOPT='NON CYCLIC - STIS 0 - SC SCHEME - TABULATED EXP'
           CALL MCGFLX(MCGFFIS,MCGFFAS,MCGSCAS,MCGFFAL,CYCLIC,KPSYS,
     1          IMPX,IPTRK,IFTRAK,IPMACR,NDIM,NFI,NUNKNO,NLONG,NBREG,
     2          NSOU,NGRP,NGEFF,INGIND,NZON,MATALB,V,FUNKNO,SUNKNO,
     3          NBMIX,NANI,MAXI,IAAC,KRYL,ISCR,NMU,NANGL,NMAX,LC,EPSI,
     4          CAZ0,CAZ1,CAZ2,CPO,ZMU,WZMU,LFORW,PACA,NLIN,NFUNL,KEY,
     5          KEYCUR,SIGAL,LPS,REPS,EPS,ITST,INCONV,LNCONV,REBFLG,
     6          STIS,NPJJM,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
         ELSEIF(ISCH.EQ.0) THEN
*          Diamond-Differencing Scheme
           TOPT='NON CYCLIC - STIS 0 - DD0 SCHEME'
           CALL MCGFLX(MCGFFIS,MCGFFAS,MCGDDFS,MCGFFAL,CYCLIC,KPSYS,
     1          IMPX,IPTRK,IFTRAK,IPMACR,NDIM,NFI,NUNKNO,NLONG,NBREG,
     2          NSOU,NGRP,NGEFF,INGIND,NZON,MATALB,V,FUNKNO,SUNKNO,
     3          NBMIX,NANI,MAXI,IAAC,KRYL,ISCR,NMU,NANGL,NMAX,LC,EPSI,
     4          CAZ0,CAZ1,CAZ2,CPO,ZMU,WZMU,LFORW,PACA,NLIN,NFUNL,KEY,
     5          KEYCUR,SIGAL,LPS,REPS,EPS,ITST,INCONV,LNCONV,REBFLG,
     6          STIS,NPJJM,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
         ELSEIF(ISCH.EQ.-1) THEN
*          Step-Characteristics Scheme with Exact Exponentials
           TOPT='NON CYCLIC - STIS 0 - SC SCHEME - EXACT EXP'
           CALL MCGFLX(MCGFFIS,MCGFFAS,MCGSCES,MCGFFAL,CYCLIC,KPSYS,
     1          IMPX,IPTRK,IFTRAK,IPMACR,NDIM,NFI,NUNKNO,NLONG,NBREG,
     2          NSOU,NGRP,NGEFF,INGIND,NZON,MATALB,V,FUNKNO,SUNKNO,
     3          NBMIX,NANI,MAXI,IAAC,KRYL,ISCR,NMU,NANGL,NMAX,LC,EPSI,
     4          CAZ0,CAZ1,CAZ2,CPO,ZMU,WZMU,LFORW,PACA,NLIN,NFUNL,KEY,
     5          KEYCUR,SIGAL,LPS,REPS,EPS,ITST,INCONV,LNCONV,REBFLG,
     6          STIS,NPJJM,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
*********'Source Term Isolation' Strategy turned on
         ELSEIF(ISCH.EQ.11) THEN
*          Step-Characteristics Scheme with Tabulated Exponentials
           TOPT='NON CYCLIC - STIS 1 - SC SCHEME - TABULATED EXP'
           CALL MCGFLX(MCGFFI,MCGFFA,MCGSCA,MCGFFAL,CYCLIC,KPSYS,IMPX,
     1          IPTRK,IFTRAK,IPMACR,NDIM,NFI,NUNKNO,NLONG,NBREG,NSOU,
     2          NGRP,NGEFF,INGIND,NZON,MATALB,V,FUNKNO,SUNKNO,NBMIX,
     3          NANI,MAXI,IAAC,KRYL,ISCR,NMU,NANGL,NMAX,LC,EPSI,CAZ0,
     4          CAZ1,CAZ2,CPO,ZMU,WZMU,LFORW,PACA,NLIN,NFUNL,KEY,KEYCUR,
     5          SIGAL,LPS,REPS,EPS,ITST,INCONV,LNCONV,REBFLG,STIS,
     6          NPJJM,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
         ELSEIF(ISCH.EQ.10) THEN
*          Diamond-Differencing Scheme
           TOPT='NON CYCLIC - STIS 1 - DD0 SCHEME'
           CALL MCGFLX(MCGFFI,MCGFFA,MCGDDF,MCGFFAL,CYCLIC,KPSYS,IMPX,
     1          IPTRK,IFTRAK,IPMACR,NDIM,NFI,NUNKNO,NLONG,NBREG,NSOU,
     2          NGRP,NGEFF,INGIND,NZON,MATALB,V,FUNKNO,SUNKNO,NBMIX,
     3          NANI,MAXI,IAAC,KRYL,ISCR,NMU,NANGL,NMAX,LC,EPSI,CAZ0,
     4          CAZ1,CAZ2,CPO,ZMU,WZMU,LFORW,PACA,NLIN,NFUNL,KEY,KEYCUR,
     5          SIGAL,LPS,REPS,EPS,ITST,INCONV,LNCONV,REBFLG,STIS,
     6          NPJJM,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
         ELSEIF(ISCH.EQ.9) THEN
*          Step-Characteristics Scheme with Exact Exponentials
           TOPT='NON CYCLIC - STIS 1 - SC SCHEME - EXACT EXP'
           CALL MCGFLX(MCGFFI,MCGFFA,MCGSCE,MCGFFAL,CYCLIC,KPSYS,IMPX,
     1          IPTRK,IFTRAK,IPMACR,NDIM,NFI,NUNKNO,NLONG,NBREG,NSOU,
     2          NGRP,NGEFF,INGIND,NZON,MATALB,V,FUNKNO,SUNKNO,NBMIX,
     3          NANI,MAXI,IAAC,KRYL,ISCR,NMU,NANGL,NMAX,LC,EPSI,CAZ0,
     4          CAZ1,CAZ2,CPO,ZMU,WZMU,LFORW,PACA,NLIN,NFUNL,KEY,KEYCUR,
     5          SIGAL,LPS,REPS,EPS,ITST,INCONV,LNCONV,REBFLG,STIS,
     6          NPJJM,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
*********'MOCC/MCI' Iterative Strategy
         ELSEIF(ISCH.EQ.-9) THEN
*          Step-Characteristics Scheme with Tabulated Exponentials
           TOPT='NON CYCLIC - STIS -1 - SC SCHEME - TABULATED EXP'
           CALL MCGFLX(MCGFFIT,MCGFFAT,MCGSCAT,MCGFFAL,CYCLIC,KPSYS,
     1          IMPX,IPTRK,IFTRAK,IPMACR,NDIM,NFI,NUNKNO,NLONG,NBREG,
     2          NSOU,NGRP,NGEFF,INGIND,NZON,MATALB,V,FUNKNO,SUNKNO,
     3          NBMIX,NANI,MAXI,IAAC,KRYL,ISCR,NMU,NANGL,NMAX,LC,EPSI,
     4          CAZ0,CAZ1,CAZ2,CPO,ZMU,WZMU,LFORW,PACA,NLIN,NFUNL,KEY,
     5          KEYCUR,SIGAL,LPS,REPS,EPS,ITST,INCONV,LNCONV,REBFLG,
     6          STIS,NPJJM,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
         ELSEIF(ISCH.EQ.-10) THEN
*          Diamond-Differencing Scheme
           TOPT='NON CYCLIC - STIS -1 - DD0 SCHEME'
           CALL MCGFLX(MCGFFIT,MCGFFAT,MCGDDFT,MCGFFAL,CYCLIC,KPSYS,
     1          IMPX,IPTRK,IFTRAK,IPMACR,NDIM,NFI,NUNKNO,NLONG,NBREG,
     2          NSOU,NGRP,NGEFF,INGIND,NZON,MATALB,V,FUNKNO,SUNKNO,
     3          NBMIX,NANI,MAXI,IAAC,KRYL,ISCR,NMU,NANGL,NMAX,LC,EPSI,
     4          CAZ0,CAZ1,CAZ2,CPO,ZMU,WZMU,LFORW,PACA,NLIN,NFUNL,KEY,
     5          KEYCUR,SIGAL,LPS,REPS,EPS,ITST,INCONV,LNCONV,REBFLG,
     6          STIS,NPJJM,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
         ELSEIF(ISCH.EQ.-11) THEN
*          Step-Characteristics Scheme with Exact Exponentials
           TOPT='NON CYCLIC - STIS -1 - SC SCHEME - EXACT EXP'
           CALL MCGFLX(MCGFFIT,MCGFFAT,MCGSCET,MCGFFAL,CYCLIC,KPSYS,
     1          IMPX,IPTRK,IFTRAK,IPMACR,NDIM,NFI,NUNKNO,NLONG,NBREG,
     2          NSOU,NGRP,NGEFF,INGIND,NZON,MATALB,V,FUNKNO,SUNKNO,
     3          NBMIX,NANI,MAXI,IAAC,KRYL,ISCR,NMU,NANGL,NMAX,LC,EPSI,
     4          CAZ0,CAZ1,CAZ2,CPO,ZMU,WZMU,LFORW,PACA,NLIN,NFUNL,KEY,
     5          KEYCUR,SIGAL,LPS,REPS,EPS,ITST,INCONV,LNCONV,REBFLG,
     6          STIS,NPJJM,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
         ELSEIF(ISCH.EQ.199) THEN
*          Lin.-Disc.-Characteristics Scheme with Exact Exponentials
           TOPT='NON CYCLIC - STIS 0 - LDC SCHEME - EXACT EXP'
           CALL MCGFLX(MCGFFIT,MCGFFAT,MCGSCEL,MCGFFAL,CYCLIC,KPSYS,
     1          IMPX,IPTRK,IFTRAK,IPMACR,NDIM,NFI,NUNKNO,NLONG,NBREG,
     2          NSOU,NGRP,NGEFF,INGIND,NZON,MATALB,V,FUNKNO,SUNKNO,
     3          NBMIX,NANI,MAXI,IAAC,KRYL,ISCR,NMU,NANGL,NMAX,LC,EPSI,
     4          CAZ0,CAZ1,CAZ2,CPO,ZMU,WZMU,LFORW,PACA,NLIN,NFUNL,KEY,
     5          KEYCUR,SIGAL,LPS,REPS,EPS,ITST,INCONV,LNCONV,REBFLG,
     6          STIS,NPJJM,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
         ELSEIF(ISCH.EQ.200) THEN
*          Diamond-Differencing Scheme
           TOPT='NON CYCLIC - STIS 0 - DD1 SCHEME'
           CALL MCGFLX(MCGFFIT,MCGFFAT,MCGDDFL,MCGFFAL,CYCLIC,KPSYS,
     1          IMPX,IPTRK,IFTRAK,IPMACR,NDIM,NFI,NUNKNO,NLONG,NBREG,
     2          NSOU,NGRP,NGEFF,INGIND,NZON,MATALB,V,FUNKNO,SUNKNO,
     3          NBMIX,NANI,MAXI,IAAC,KRYL,ISCR,NMU,NANGL,NMAX,LC,EPSI,
     4          CAZ0,CAZ1,CAZ2,CPO,ZMU,WZMU,LFORW,PACA,NLIN,NFUNL,KEY,
     5          KEYCUR,SIGAL,LPS,REPS,EPS,ITST,INCONV,LNCONV,REBFLG,
     6          STIS,NPJJM,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
         ELSEIF(ISCH.EQ.201) THEN
*          Lin.-Disc.-Characteristics Scheme with Tabulated Exponentials
           TOPT='NON CYCLIC - STIS 0 - LDC SCHEME - TABULATED EXP'
           CALL MCGFLX(MCGFFIT,MCGFFAT,MCGSCAL,MCGFFAL,CYCLIC,KPSYS,
     1          IMPX,IPTRK,IFTRAK,IPMACR,NDIM,NFI,NUNKNO,NLONG,NBREG,
     2          NSOU,NGRP,NGEFF,INGIND,NZON,MATALB,V,FUNKNO,SUNKNO,
     3          NBMIX,NANI,MAXI,IAAC,KRYL,ISCR,NMU,NANGL,NMAX,LC,EPSI,
     4          CAZ0,CAZ1,CAZ2,CPO,ZMU,WZMU,LFORW,PACA,NLIN,NFUNL,KEY,
     5          KEYCUR,SIGAL,LPS,REPS,EPS,ITST,INCONV,LNCONV,REBFLG,
     6          STIS,NPJJM,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
         ELSE
           CALL XABORT('MCCGF: NON-CYCLIC SCHEME NOT IMPLEMENTED')
         ENDIF
*     ------------------------------------
      ENDIF
      ENDIF
      DEALLOCATE(INCONV,SIGAL,CPO,CAZ2,CAZ1,CAZ0,MATALB,KPSYS)
*---
* PRINT RESULTS
*---
      IF(IMPX.GT.0) WRITE(IUNOUT,50) TOPT
      IF((IDIFC.EQ.0).AND.(MAXI.GT.1).AND.(IMPX.GT.1)) THEN     
      DO II=1,NGEFF
         IG=INGIND(II)
         WRITE(IUNOUT,100) IG
         IF(IMPX.GT.3) 
     1      WRITE(IUNOUT,200) (SUNKNO(KEYFLX(I),IG),I=1,NBREG)
         ITEMP=ITST(II)
         TEMP=EPS(II)
         IF((ITEMP.EQ.MAXI).AND.(TEMP.GT.EPSI)) WRITE(IUNOUT,60) 
         WRITE(IUNOUT,70) ITEMP,TEMP
         IF(IMPX.GT.2) THEN
           WRITE(IUNOUT,101) (REPS(MAXI*(II-1)+I),I=1,MIN(ITEMP,100))
         ENDIF
         IF(IMPX.GT.4) THEN
            WRITE(IUNOUT,400) (FUNKNO(I,IG),I=1,NUNKNO)
         ELSEIF(IMPX.GT.3) THEN
            WRITE(IUNOUT,300) (FUNKNO(KEYFLX(I),IG),I=1,NBREG)
         ENDIF
      ENDDO
      ENDIF
      DEALLOCATE(ITST,EPS,REPS,INGIND)
*----
* DOUBLE HETEROGENEITY TREATMENT
*----
 40   IF(LBIHET) THEN
         NBMIX=NBMIXG
         NBREG=NREGAR
         DO IG=1,NGRP
           IOFSET=NPSYS(IG)
           IF(IOFSET.NE.0) THEN
             DO I=1,NUNKNO
               SUNKNO(I,IG)=SUNKN((IG-1)*NUNKNO+I)
             ENDDO
             JPSYS=LCMGIL(IPSYS,IOFSET)
             CALL DOORFB3(JPSYS,IPTRK,IMPX,NBMIX,NBREG,NUNKNO,KEYFLX,
     1       SUNKNO(1,IG),FUNKNO(1,IG))
           ENDIF
         ENDDO
         DEALLOCATE(SUNKN)
      ENDIF
      RETURN
*
 50   FORMAT(9X,18H M O C PARAMETERS:,2X,A72)
 60   FORMAT(49H *** WARNING *** MAXIMUM NUMBER OF MCCG ITERATION,
     1       10HS REACHED.)
 70   FORMAT(34H MCCGF: NUMBER OF MCCG ITERATIONS=,I4,
     1       11H  ACCURACY=,1P,E11.4,1H.)
 100  FORMAT(9X,8H GROUP (,I4,4H) : )
 101  FORMAT('----  EPS  ----'/(1P,6E16.6))
 200  FORMAT(33H N E U T R O N    S O U R C E S :/(1P,6(5X,E15.7)))
 300  FORMAT(31H N E U T R O N    F L U X E S :/(1P,6(5X,E15.7)))
 400  FORMAT(31H U N K N O W N    V E C T O R :/(1P,6(5X,E15.7)))
      END