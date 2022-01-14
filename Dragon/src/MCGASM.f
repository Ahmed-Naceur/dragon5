*DECK MCGASM
      SUBROUTINE MCGASM(SUBPJJ,SUBDS2,SUBDSP,SUBDSC,IPTRK,KPSYS,IPRINT,
     1                  IFTRAK,NANI,NGEFF,NFI,NREG,NLONG,M,NMU,NANGL,
     2                  N2MAX,LC,NDIM,NGIND,CYCLIC,ISCR,CAZ0,CAZ1,CAZ2,
     3                  CPO,LC0,PACA,LPS,LTMT,NPJJM,LACA,LPJJ,LPJJAN,
     4                  SIGAL,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,ISTRM)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Preconditioning matrices calculation based on the Algebraic Collapsing
* Acceleration developed by I. R. Suslov and R. Le Tellier 
* or Self-Collision Probabilities acceleration developed by G.J. Wu 
* and R. Roy.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): I. Suslov and R. Le Tellier
*
*Parameters: input/output
* SUBPJJ  PJJ calculation subroutine.
* SUBDS2  ACA coefficients summation subroutine.
* SUBDSP  ACA coefficients position subroutine.
* SUBDSC  ACA coefficients calculation subroutine.
* IPTRK   pointer to the tracking (L_TRACK signature).
* KPSYS   pointer array for each group properties.
* IPRINT  print parameter (equal to zero for no print).
* IFTRAK  tracking file unit number if IOFSET=0.
* NANI    number of Legendre orders.
* NGEFF   number of groups to process.
* NFI     total number of volumes and surfaces for which specific values
*         of the neutron flux and reactions rates are required.
* NREG    number of volumes for which specific values
*         of the neutron flux and reactions rates are required.
* NLONG   order of the corrective system.
* M       number of material mixtures.
* NMU     order of the polar quadrature in 2D / 1 in 3D.
* NANGL   number of tracking angles in the plan.
* N2MAX   maximum number of elements in a track.
* LC      dimension of vector MCU.
* NDIM    number of dimensions for the geometry.
* NGIND   index of the groups to process.
* CYCLIC  flag set to .true. for cyclic tracking.
* ISCR    SCR preconditionning flag.
* CAZ0    cosines of the tracking polar angles in 3D.
* CAZ1    first cosines of the different tracking azimuthal angles.
* CAZ2    second cosines of the different tracking azimuthal angles.
* CPO     cosines of the different tracking polar angles in 2D.
* PACA    type of preconditioner to solve the ACA corrective system.
* LC0     used in ILU0-ACA acceleration.
* LPS     dimension of JS.
* LTMT    tracking merging flag.
* NPJJM   number of pjj modes to store for STIS option.
* LACA    ACA flag.
* LPJJ    PJJ flag.
* LPJJAN  anisotropic PJJ flag.
* SIGAL   albedos and total cross sections array.
* LPRISM  3D prismatic extended tracking flag.
* N2REG   number of regions in the 2D tracking if LPRISM.
* N2SOU   number of external surfaces in the 2D tracking if LPRISM.
* NZP     number of z-plans if LPRISM.
* DELU    input track spacing for 3D track reconstruction if LPRISM.
* FACSYM  tracking symmetry factor for maximum track length if LPRISM.
* ISTRM   type of streaming effect:
*         =1 no streaming effect;
*         =2 isotropic streaming effect;
*         =3 anisotropic streaming effect.
*
*Reference:
*  Igor R. Suslov, "An Algebraic Collapsing Acceleration in Long 
*  Characteristics Transport Theory" Proc. of 11-th Symposium of 
*  AER, 178/9-188, Csopak, September 2001.
*  \\\\
*  Igor R. Suslov, "Solution of Transport Equation in 2- and 3-
*  dimensional Irregular Geometry by the Method of Characteristics"
*  Int. Conf. Mathematical. Methods and Supercomputing in Nuclear
*  Applications, Karlsruhe, 1994.
*  \\\\
*  G.J. Wu and R. Roy, "Acceleration Techniques for Trajectory-based
*  Deterministic 3D Transport Solvers",
*  Ann. Nucl. Energy, 30, 567-583 (2003).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,KPSYS(NGEFF)
      INTEGER IPRINT,IFTRAK,NGEFF,NFI,NREG,NLONG,M,NMU,NANGL,N2MAX,LC,
     1 NDIM,NGIND(NGEFF),LC0,PACA,LPS,NPJJM,N2REG,N2SOU,NZP,ISTRM
      REAL CPO(NMU),SIGAL(-6:M,NGEFF),DELU,FACSYM
      DOUBLE PRECISION CAZ0(NANGL),CAZ1(NANGL),CAZ2(NANGL)
      LOGICAL LTMT,LACA,LPJJ,LPJJAN,LPRISM,LFORC,CYCLIC
      EXTERNAL SUBPJJ,SUBDS2,SUBDSP,SUBDSC
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPSYS
      INTEGER NCODE(6),SSYM
      REAL T1,T2,T3,XMUANG(1)
      DOUBLE PRECISION WEIGHT,WEIGHT0
      CHARACTER TEXT4*4
      INTEGER, TARGET, SAVE, DIMENSION(1) :: IDUMMY
      INTEGER, TARGET, SAVE, DIMENSION(1,1) :: I2DUMMY
      REAL, TARGET, SAVE, DIMENSION(1) :: DUMMY
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, POINTER, DIMENSION(:) :: NZON,KM,IM,MCU,NZONA,IPERM,
     1 JU,IM0,MCU0,IS,JS
      INTEGER, POINTER, DIMENSION(:,:) :: PJJIND
      REAL, POINTER, DIMENSION(:) :: ZMU,WZMU,V,VA
      TYPE(C_PTR) :: ZMU_PTR,WZMU_PTR,NZON_PTR,V_PTR,KM_PTR,IM_PTR,
     1 MCU_PTR,NZONA_PTR,VA_PTR,IPERM_PTR,JU_PTR,IM0_PTR,MCU0_PTR,
     2 PJJIND_PTR,IS_PTR,JS_PTR
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NOM,NOM0,PREV,NEXT,NOM3D,
     1 NOM3D0,KANGL,KEYANI
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISGNR
      REAL, ALLOCATABLE, DIMENSION(:) :: ZMUA,WZMUA,PJJ,PSJ,XSW,CQ,
     1 DIAGQ,DIAGFR,CFR,WORK,LUDF,LUCF,PJJX,PJJY,PJJZ,PJJXI,PJJYI,
     2 PJJZI,PSJX,PSJY,PSJZ
      REAL, ALLOCATABLE, DIMENSION(:,:) :: RHARM
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: TRHAR
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: H,H0,HH,H3D,H3D0,
     1 PJJD,PJJDX,PJJDY,PJJDZ,PJJDXI,PJJDYI,PJJDZI
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: DIAGF,CF
      INTEGER, POINTER, DIMENSION(:) :: INDREG
      REAL, POINTER, DIMENSION(:) :: ZZZ
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: VNORF,CMU,CMUI,SMU,
     1 SMUI,TMU,TMUI
      TYPE(C_PTR) :: INDREG_PTR,ZZZ_PTR,VNORF_PTR,CMU_PTR,CMUI_PTR,
     1 SMU_PTR,SMUI_PTR,TMU_PTR,TMUI_PTR
*----
*  INITIALIZE POINTERS
*----
      KM=>IDUMMY
      IM=>IDUMMY
      MCU=>IDUMMY
      NZONA=>IDUMMY
      VA=>DUMMY
      IPERM=>IDUMMY
      JU=>IDUMMY
      IM0=>IDUMMY
      MCU0=>IDUMMY
      IS=>IDUMMY
      JS=>IDUMMY
      PJJIND=>I2DUMMY
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(DIAGF(NLONG,NGEFF),CF(LC,NGEFF))
*----
*  TRACKING INFORMATION PINNING AND MEMORY ALLOCATION
*   NZON    index-number of the mixture type assigned to each volume.
*   NZONA   index-number of the mixture type assigned to each volume for
*           ACA.
*   ISGNR   array of the spherical harmonics signs for the different
*           reflections.
*   V       volumes and surfaces.
*   VA      volumes and surfaces reordered for ACA.
*   ZMU     polar quadrature set in 2D.
*   WZMU    polar quadrature set in 2D.
*   KM      used in ACA acceleration.
*   MCU     used in ACA acceleration.
*   IM      used in ACA acceleration.
*   JU      used in ACA  acceleration for ilu0.
*   IPERM   permutation array for the unknowns of the corrective system
*           for ilu0.
*   LC0     used in ILU0-ACA acceleration.
*   IM0     used in ILU0-ACA acceleration.
*   MCU0    used in ILU0-ACA acceleration.
*   IS      arrays for surfaces neighbors
*   JS      JS(IS(ISOUT)+1:IS(ISOUT+1)) give the neighboring regions to
*           surface ISOUT.
*   PJJIND  index of the modes for STIS option.
*---- 
*---
* GENERATE ALL SIGNS FOR SPHERICAL HARMONICS
*---
      IF(NDIM.EQ.1) THEN
         NFUNL=NANI
         NMOD=2
      ELSE IF((.NOT.LPRISM).AND.(NDIM.EQ.2)) THEN
         NFUNL=NANI*(NANI+1)/2
         NMOD=4
      ELSE ! NDIM.EQ.3
         NFUNL=NANI*NANI
         NMOD=8
      ENDIF
      ALLOCATE(ISGNR(NMOD,NFUNL),KEYANI(NFUNL))
      CALL MOCIK3(NANI-1,NFUNL,NMOD,ISGNR,KEYANI) 
      DEALLOCATE(KEYANI)
*---
* ASSEMBLY OF SCR AND ACA MATRICES
*---
*     recover polar quadrature
      CALL LCMGPD(IPTRK,'ZMU$MCCG',ZMU_PTR)
      CALL LCMGPD(IPTRK,'WZMU$MCCG',WZMU_PTR)
      CALL C_F_POINTER(ZMU_PTR,ZMU,(/ NMU /))
      CALL C_F_POINTER(WZMU_PTR,WZMU,(/ NMU /))
*     recover MATALB and VOLSUR
      CALL LCMGPD(IPTRK,'NZON$MCCG',NZON_PTR)
      CALL LCMGPD(IPTRK,'V$MCCG',V_PTR)
      CALL C_F_POINTER(NZON_PTR,NZON,(/ NFI /))
      CALL C_F_POINTER(V_PTR,V,(/ NFI /))
      IF(LACA) THEN
*        recover connection matrices
         CALL LCMGPD(IPTRK,'KM$MCCG',KM_PTR)
         CALL LCMGPD(IPTRK,'IM$MCCG',IM_PTR)
         CALL LCMGPD(IPTRK,'MCU$MCCG',MCU_PTR)
         CALL C_F_POINTER(KM_PTR,KM,(/ NFI /))
         CALL C_F_POINTER(IM_PTR,IM,(/ NLONG+1 /))
         CALL C_F_POINTER(MCU_PTR,MCU,(/ LC /))
*        recover modified MATALB and VOLSUR
         CALL LCMGPD(IPTRK,'NZONA$MCCG',NZONA_PTR)
         CALL LCMGPD(IPTRK,'VA$MCCG',VA_PTR)
         CALL C_F_POINTER(NZONA_PTR,NZONA,(/ NFI /))
         CALL C_F_POINTER(VA_PTR,VA,(/ NFI /))
*        recover permutation array
         CALL LCMGPD(IPTRK,'INVPI$MCCG',IPERM_PTR)
         CALL C_F_POINTER(IPERM_PTR,IPERM,(/ NFI /))
         IF(PACA.GE.2) THEN
            CALL LCMGPD(IPTRK,'JU$MCCG',JU_PTR)
            CALL C_F_POINTER(JU_PTR,JU,(/ NLONG /))
            IF(PACA.EQ.3) THEN
               CALL LCMLEN(IPTRK,'IM0$MCCG',ILONG1,ITYLCM)
               CALL LCMLEN(IPTRK,'MCU0$MCCG',ILONG2,ITYLCM)
               CALL LCMGPD(IPTRK,'IM0$MCCG',IM0_PTR)
               CALL LCMGPD(IPTRK,'MCU0$MCCG',MCU0_PTR)
               CALL C_F_POINTER(IM0_PTR,IM0,(/ ILONG1 /))
               CALL C_F_POINTER(MCU0_PTR,MCU0,(/ ILONG2 /))
            ENDIF
         ENDIF
      ENDIF
      IF((.NOT.CYCLIC).AND.(ISCR.GT.0)) THEN
*        recover (IS,JS) array for surfaces neighbors identification
         CALL LCMGPD(IPTRK,'IS$MCCG',IS_PTR)
         CALL LCMGPD(IPTRK,'JS$MCCG',JS_PTR)
         CALL C_F_POINTER(IS_PTR,IS,(/ NFI-NREG+1 /))
         CALL C_F_POINTER(JS_PTR,JS,(/ LPS /))
      ENDIF
      IF(LPJJAN) THEN
         CALL LCMGPD(IPTRK,'PJJIND$MCCG',PJJIND_PTR)
         CALL C_F_POINTER(PJJIND_PTR,PJJIND,(/ NPJJM,2 /))
      ENDIF
      IF(LPRISM) THEN
         CALL LCMGET(IPTRK,'NCODE',NCODE)
         IF(NCODE(6).EQ.30) THEN
            IF(NCODE(5).EQ.30) THEN
*           Z- and Z+ surfaces symmetry
               SSYM=2
            ELSE
*           Z+ symmetry
               SSYM=1
            ENDIF
         ELSE
            SSYM=0
         ENDIF
         NMAX=(INT(FACSYM)+1)*N2MAX*(NZP+2)
      ELSE
         NMAX=N2MAX
      ENDIF
      ALLOCATE(NOM(N2MAX))
      ALLOCATE(H(N2MAX),HH(N2MAX),ZMUA(NMU),WZMUA(NMU))
      ALLOCATE(TRHAR(NMU,NFUNL,NMOD),RHARM(NMU,NFUNL))
      IF(LPJJ) THEN
*     Self-Collision Probabilities
         ALLOCATE(PJJ(NREG*NPJJM*NGEFF),PJJX(NREG*NPJJM*NGEFF),
     >            PJJXI(NREG*NPJJM*NGEFF),PJJY(NREG*NPJJM*NGEFF),
     >            PJJYI(NREG*NPJJM*NGEFF),PJJZ(NREG*NPJJM*NGEFF),
     >            PJJZI(NREG*NPJJM*NGEFF))
         IF(LPS.GT.0) THEN
            ALLOCATE(PSJ(LPS*NGEFF),PSJX(LPS*NGEFF),PSJY(LPS*NGEFF),
     >               PSJZ(LPS*NGEFF))
            CALL XDRSET(PSJ,LPS*NGEFF,0.0)
            CALL XDRSET(PSJX,LPS*NGEFF,0.0)
            CALL XDRSET(PSJY,LPS*NGEFF,0.0)
            CALL XDRSET(PSJZ,LPS*NGEFF,0.0)
         ENDIF
         ALLOCATE(PJJD(NREG*NPJJM*NGEFF),PJJDX(NREG*NPJJM*NGEFF),
     >            PJJDXI(NREG*NPJJM*NGEFF),PJJDY(NREG*NPJJM*NGEFF),
     >            PJJDYI(NREG*NPJJM*NGEFF),PJJDZ(NREG*NPJJM*NGEFF),
     >            PJJDZI(NREG*NPJJM*NGEFF))
         CALL XDDSET(PJJD,NREG*NPJJM*NGEFF,0.D0)
         CALL XDDSET(PJJDX,NREG*NPJJM*NGEFF,0.D0)
         CALL XDDSET(PJJDXI,NREG*NPJJM*NGEFF,0.D0)
         CALL XDDSET(PJJDY,NREG*NPJJM*NGEFF,0.D0)
         CALL XDDSET(PJJDYI,NREG*NPJJM*NGEFF,0.D0)
         CALL XDDSET(PJJDZ,NREG*NPJJM*NGEFF,0.D0)
         CALL XDDSET(PJJDZI,NREG*NPJJM*NGEFF,0.D0)
      ENDIF
      IF(LACA) THEN
*     Algebraic Collapsing Acceleration
         ALLOCATE(XSW((M+1)*NANI*NGEFF))
         IF(LTMT) ALLOCATE(NOM0(N2MAX),H0(N2MAX))
         ALLOCATE(PREV(NMAX),NEXT(NMAX))
         ALLOCATE(CQ(LC*NGEFF),DIAGQ(NLONG*NGEFF),DIAGFR(NLONG),CFR(LC),
     1   WORK(6*NMAX))
         IF(PACA.GE.2) THEN
            ALLOCATE(LUDF(NLONG))
            IF(LC0.GT.0) ALLOCATE(LUCF(LC0))
         ENDIF
         DO II=1,NGEFF
            JPSYS=KPSYS(II)
            CALL LCMGET(JPSYS,'DRAGON-S0XSC',XSW((II-1)*(M+1)+1))
         ENDDO
         CALL XDRSET(DIAGQ,NLONG*NGEFF,0.0)
         CALL XDRSET(CQ,LC*NGEFF,0.0)
         CALL XDDSET(DIAGF,NLONG*NGEFF,0.D0)
         CALL XDDSET(CF,LC*NGEFF,0.D0)
      ENDIF 
*
      NMUA=NMU
      DO IE=1,NMUA
         ZMUA(IE)=ZMU(IE)
         WZMUA(IE)=WZMU(IE)       
      ENDDO
      IF((.NOT.LPRISM).AND.(NDIM.EQ.2).AND.LTMT) THEN
         NMUA=1
         ZMUI=ZMUA(1)
         W=WZMUA(1)
         TEMP=ZMUI*W
         DO IE=2,NMU
            ZMUI=ZMUA(IE)
            WZMUI=WZMUA(IE)
            W=W+WZMUI
            TEMP=TEMP+ZMUI*WZMUI
         ENDDO               
         ZMUA=TEMP/W
         WZMUA=W
      ENDIF
*---
* SUMMATION UPON THE TRACKING
*---
      REWIND IFTRAK
      READ(IFTRAK) TEXT4,NCOMNT,NBTR,IFMT
      DO ICOM=1,NCOMNT
         READ(IFTRAK)
      ENDDO
      READ(IFTRAK) (NITMA,II=1,7),MXSUB,NITMA
      DO ICOM=1,6
         READ(IFTRAK)
      ENDDO
      CALL KDRCPU(T1)
      IANGL0=0
      IPANG=1
      NMERG=0
      NTR=0
      NTRTMT=0
      NSE=0
      NSETMT=0
      LFORC=.FALSE.
      NTPROC=1
*
      ALLOCATE(KANGL(MXSUB))
      IF(LPRISM) THEN
*     3D PRISMATIC GEOMETRY CONSTRUCTED FROM A 2D TRACKING
      N3TR=0
      N3TRTMT=0
      N3SE=0
      N3SETMT=0
      ALLOCATE(NOM3D(NMAX),H3D(NMAX))
      IF(LTMT) ALLOCATE(NOM3D0(2*NMAX),H3D0(2*NMAX))
      CALL LCMSIX(IPTRK,'PROJECTION',1)
      CALL LCMGPD(IPTRK,'IND2T3',INDREG_PTR)
      CALL LCMGPD(IPTRK,'ZCOORD',ZZZ_PTR)
      CALL LCMGPD(IPTRK,'VNORF',VNORF_PTR)
      CALL LCMGPD(IPTRK,'CMU',CMU_PTR)
      CALL LCMGPD(IPTRK,'CMUI',CMUI_PTR)
      CALL LCMGPD(IPTRK,'SMU',SMU_PTR)
      CALL LCMGPD(IPTRK,'SMUI',SMUI_PTR)
      CALL LCMGPD(IPTRK,'TMU',TMU_PTR)
      CALL LCMGPD(IPTRK,'TMUI',TMUI_PTR)
      CALL C_F_POINTER(INDREG_PTR,INDREG,(/ (N2SOU+N2REG+1)*(NZP+1) /))
      CALL C_F_POINTER(ZZZ_PTR,ZZZ,(/ NZP+1 /))
      CALL C_F_POINTER(VNORF_PTR,VNORF,(/ NREG*NANGL*NMU*2 /))
      CALL C_F_POINTER(CMU_PTR,CMU,(/ NMU /))
      CALL C_F_POINTER(CMUI_PTR,CMUI,(/ NMU /))
      CALL C_F_POINTER(SMU_PTR,SMU,(/ NMU /))
      CALL C_F_POINTER(SMUI_PTR,SMUI,(/ NMU /))
      CALL C_F_POINTER(TMU_PTR,TMU,(/ NMU /))
      CALL C_F_POINTER(TMUI_PTR,TMUI,(/ NMU /))
      CALL LCMSIX(IPTRK,'PROJECTION',2)
      DO 10 ILINE=1,NBTR
         READ(IFTRAK) NSUB,NSEG,WEIGHT,(KANGL(II),II=1,NSUB),
     1        (NOM(II),II=1,NSEG),(H(II),II=1,NSEG)
         IF(NSUB.GT.MXSUB) CALL XABORT('MCGASM: MXSUB OVERFLOW.')
         IANGL=KANGL(1)
         IF(LPJJ) THEN
*        ----------------------------
*        Self-Collision Probabilities
*        ----------------------------
            NR2SE=NSEG-2
            HH=0.0
            DO II=0,NR2SE
               HH(II+2)=HH(II+1)+H(II+2)
            ENDDO
            CALL MCGPTS(SUBPJJ,NFI,NREG,M,NANI,NFUNL,NANGL,NMU,NMOD,
     1           LPS,NPJJM,NGEFF,IANGL,IANGL0,NSEG,ISGNR,NZON,NOM,
     2           IS,JS,PJJIND,WEIGHT,CPO,CAZ1,CAZ2,RHARM,ZMU,WZMU,TRHAR,
     3           SIGAL,HH,PSJ,PJJD,LPJJAN,NR2SE,NMAX,NZP,N2REG,N2SOU,
     4           DELU,INDREG,NOM3D,H3D,ZZZ,VNORF,CMU,CMUI,SMU,SMUI,TMU,
     5           TMUI,SSYM)
         ENDIF
         IF(LTMT) THEN
*        ----------------
*        Tracking Merging (angle by angle) for ACA
*        ----------------
            NTR=NTR+1
            NSE=NSE+NSEG
            LFORC=(IPANG.NE.IANGL)
            IF(LFORC) THEN
               ITEMP=IANGL
               IANGL=IPANG
               IPANG=ITEMP
            ENDIF
            IF(ILINE.EQ.NBTR) LFORC=.TRUE.
            CALL MCGTMT(NMERG,NTRTMT,NSETMT,NSEG,NSEG0,NOM,NOM0,WEIGHT,
     1           WEIGHT0,H,H0,LFORC,NTPROC)
            IF(NTPROC.EQ.0) GOTO 10
         ENDIF
         IF(LACA) THEN
*        ---------------------------------
*        Algebraic Collapsing Acceleration
*        ---------------------------------
            NR2SE=NSEG-2
            HH=0.0
            DO II=0,NR2SE
               HH(II+2)=HH(II+1)+H(II+2)
            ENDDO
            CALL MCGPTA(NFI,NREG,NLONG,M,NANGL,NMU,LC,NGEFF,
     1           IANGL,NSEG,NOM,NZONA,IPERM,KM,IM,MCU,PREV,NEXT,
     2           WEIGHT,ZMU,WZMU,SIGAL,XSW,HH,DIAGQ,CQ,DIAGF,
     3           CF,WORK,LTMT,SUBDS2,SUBDSP,SUBDSC,NR2SE,NMAX,
     4           NZP,N2REG,N2SOU,DELU,INDREG,NOM3D,NOM3D0,H3D,
     5           H3D0,ZZZ,VNORF,CMU,CMUI,SMU,SMUI,TMU,TMUI,N3TR,
     6           N3TRTMT,N3SE,N3SETMT,NTPROC,SSYM)
         ENDIF
 10   CONTINUE
      IF(LTMT) THEN
*     process last integration line for ACA TMT
         CALL MCGTMT(NMERG,NTRTMT,NSETMT,NSEG,NSEG0,NOM,NOM0,WEIGHT,
     1        WEIGHT0,H,H0,LFORC,NTPROC) 
         NR2SE=NSEG-2
         HH=0.0
         DO II=0,NR2SE
            HH(II+2)=HH(II+1)+H(II+2)
         ENDDO
         CALL MCGPTA(NFI,NREG,NLONG,M,NANGL,NMU,LC,NGEFF,IANGL,NSEG,
     1        NOM,NZONA,IPERM,KM,IM,MCU,PREV,NEXT,WEIGHT,ZMU,WZMU,
     2        SIGAL,XSW,HH,DIAGQ,CQ,DIAGF,CF,WORK,LTMT,SUBDS2,SUBDSP,
     3        SUBDSC,NR2SE,NMAX,NZP,N2REG,N2SOU,DELU,INDREG,NOM3D,
     4        NOM3D0,H3D,H3D0,ZZZ,VNORF,CMU,CMUI,SMU,SMUI,TMU,TMUI,
     5        N3TR,N3TRTMT,N3SE,N3SETMT,NTPROC,SSYM)
         DEALLOCATE(H3D0,NOM3D0)
      ENDIF
      DEALLOCATE(H3D,NOM3D)
      ELSE
*     REGULAR 2D OR 3D GEOMETRY
      DO 20 ILINE=1,NBTR
         READ(IFTRAK) NSUB,NSEG,WEIGHT,(KANGL(II),II=1,NSUB),
     1        (NOM(II),II=1,NSEG),(H(II),II=1,NSEG)
         IF(NSUB.GT.MXSUB) CALL XABORT('MCGASM: MXSUB OVERFLOW.')
         IANGL=KANGL(1)
         DO II=1,NSEG
            IF(NOM(II).LT.0) THEN
               NOM(II)=NREG-NOM(II)
            ELSE IF(NOM(II).EQ.0) THEN
               NOM(II)=NREG+1
            ENDIF
         ENDDO
         IF(LPJJ) THEN
*        ----------------------------
*        Self-Collision Probabilities
*        ----------------------------
            IF(LPJJAN) THEN
            IF(IANGL.NE.IANGL0) THEN
               IANGL0=IANGL
               IF(NDIM.EQ.2) THEN
                 CALL MOCCHR(NDIM,NANI-1,NFUNL,NMU,CPO,CAZ1(IANGL),
     1                       CAZ2(IANGL),RHARM)
               ELSE
                 XMUANG(1)=REAL(CAZ0(IANGL))
                 CALL MOCCHR(NDIM,NANI-1,NFUNL,1,XMUANG,CAZ1(IANGL),
     1                       CAZ2(IANGL),RHARM)
               ENDIF
               DO 23 JM=1,NMOD
               DO 22 JF=1,NFUNL
               DO 21 IE=1,NMU
                  IND2=IE+(JF-1)*NMU
                  IND1=IND2+(JM-1)*NFUNL*NMU
                  TRHAR(IE,JF,JM)=ISGNR(JM,JF)*RHARM(IE,JF)
 21            CONTINUE
 22            CONTINUE
 23            CONTINUE
            ENDIF
            ENDIF
            IF(ISTRM.LE.2) THEN
              CALL MCGDS4(SUBPJJ,NSEG,NMU,LPS,NFUNL,NMOD,NGEFF,WEIGHT,
     1           TRHAR,H,ZMU,WZMU,NOM,NZON,NFI,NREG,NDIM,M,IS,JS,PJJD,
     2           PSJ,LPJJAN,NPJJM,PJJIND,SIGAL,1,1)
            ELSE IF(ISTRM.EQ.3) THEN
*             TIBERE model
              CALL MCGDSD(NSEG,NMU,LPS,NFUNL,NMOD,NGEFF,WEIGHT,
     1           TRHAR,H,ZMU,WZMU,NOM,NZON,NFI,NREG,NDIM,M,IS,JS,PJJD,
     2           PSJ,LPJJAN,NPJJM,PJJIND,SIGAL,1,1,CAZ1(IANGL),
     3           CAZ2(IANGL),PJJDX,PJJDY,PJJDZ,PJJDXI,PJJDYI,PJJDZI,
     4           CAZ0(IANGL),PSJX,PSJY,PSJZ)
            ENDIF
         ENDIF
         IF(LTMT) THEN
*        ----------------
*        Tracking Merging (angle by angle) for ACA
*        ----------------
            NTR=NTR+1
            NSE=NSE+NSEG
            IF(ILINE.EQ.NBTR) LFORC=.TRUE.
            CALL MCGTMT(NMERG,NTRTMT,NSETMT,NSEG,NSEG0,NOM,NOM0,WEIGHT,
     1           WEIGHT0,H,H0,LFORC,NTPROC)
            IF(NTPROC.EQ.0) GOTO 20
         ENDIF
         IF(LACA) THEN
*        ---------------------------------
*        Algebraic Collapsing Acceleration
*        ---------------------------------
            DO II=1,NSEG
               NOM(II)=IPERM(NOM(II))
            ENDDO
            CALL MCGDS1(SUBDS2,SUBDSP,SUBDSC,NSEG,NMUA,NGEFF,WEIGHT,H,
     1           ZMUA,WZMUA,NOM,NZONA,NLONG,NFI,NDIM,LC,M,KM,IM,MCU,
     2           DIAGF,DIAGQ,CF,CQ,PREV,NEXT,SIGAL,XSW,WORK)
         ENDIF
 20   CONTINUE
      IF(LTMT) THEN
*     process last integration line
         CALL MCGTMT(NMERG,NTRTMT,NSETMT,NSEG,NSEG0,NOM,NOM0,WEIGHT,
     1        WEIGHT0,H,H0,LFORC,NTPROC)
         DO II=1,NSEG
            NOM(II)=IPERM(NOM(II))
         ENDDO
         CALL MCGDS1(SUBDS2,SUBDSP,SUBDSC,NSEG,NMUA,NGEFF,WEIGHT,H,
     1        ZMUA,WZMUA,NOM,NZONA,NLONG,NFI,NDIM,LC,M,KM,IM,MCU,
     2        DIAGF,DIAGQ,CF,CQ,PREV,NEXT,SIGAL,XSW,WORK)
      ENDIF
      ENDIF
*
      IF((LTMT).AND.(IPRINT.GT.1)) THEN
         WRITE(6,*) 'TRACKING MERGING: FROM TRACKING FILE'
         WRITE(6,*) '                 ',
     1              NTR,' TRACKS ->',
     2              NTRTMT,' EQUIVALENT TRACKS'
         WRITE(6,*) '                 ',
     1              NSE,' SEGMENTS ->',
     2              NSETMT, ' SEGMENTS'
         IF((.NOT.LPRISM).AND.(NDIM.EQ.2))
     1    WRITE(6,*) '                 ',
     2              NMU,' POLAR ANGLES ->',
     3              ' 1 EQUIVALENT POLAR ANGLE'
         IF(LPRISM) THEN
            WRITE(6,*) 'TRACKING MERGING: 3D PRISMATIC EXTENSION'
            WRITE(6,*) '                 ',
     1                 N3TR,' TRACKS ->',
     2                 N3TRTMT,' EQUIVALENT TRACKS'
            WRITE(6,*) '                 ',
     1                 N3SE,' SEGMENTS ->',
     2                 N3SETMT, ' SEGMENTS'
         ENDIF
      ENDIF
      CALL KDRCPU(T2)
      IF(IPRINT.GT.0) THEN
         WRITE(6,100) 'SUMMATION UPON THE TRACKING FOR ACA/SCR ',(T2-T1)
      ENDIF
*---
*  FOR ACA: (IF PACA.GE.2)
*  CALCULATION OF ILU0 PRECONDITIONER FOR BICGSTAB ITERATIONS
*---
      IF(LACA) THEN
         DO II=1,NGEFF
           IG=NGIND(II)
           JPSYS=KPSYS(II)
           IF(IPRINT.GT.2) WRITE(6,200) 'GROUP',IG
           CALL MCGDS3(NLONG,PACA,M,SIGAL(0,II),
     1          XSW((II-1)*(M+1)+1),VA,NZONA,LC,MCU,IM,JU,LC0,IM0,
     2          MCU0,DIAGF(1,II),CF(1,II),DIAGQ((II-1)*NLONG+1),DIAGFR,
     3          CFR,LUDF,LUCF)
*     
           IF(IPRINT.GT.3) THEN
              CALL PRINAM('DIAGF ',DIAGFR,NLONG)
              CALL PRINAM('CF    ',CFR,LC)
           ENDIF
           CALL LCMPUT(JPSYS,'DIAGF$MCCG',NLONG,2,DIAGFR)
           CALL LCMPUT(JPSYS,'CF$MCCG',LC,2,CFR)
           IF(PACA.GE.2) THEN
              CALL LCMPUT(JPSYS,'ILUDF$MCCG',NLONG,2,LUDF)
              IF(LC0.GT.0) CALL LCMPUT(JPSYS,'ILUCF$MCCG',LC0,2,LUCF)
           ENDIF
           IF(IPRINT.GT.3) THEN
              CALL PRINAM('DIAGQ ',DIAGQ((II-1)*NLONG+1),NLONG)
              CALL PRINAM('CQ    ',CQ((II-1)*LC+1),LC)
           ENDIF
           CALL LCMPUT(JPSYS,'CQ$MCCG',LC,2,CQ((II-1)*LC+1))
           CALL LCMPUT(JPSYS,'DIAGQ$MCCG',NLONG,2,DIAGQ((II-1)*NLONG+1))
         ENDDO
         CALL KDRCPU(T3)
         IF(IPRINT.GT.0) THEN
            WRITE(6,100) 'CALCULATION OF ACA PRECONDITIONER       ',
     1                   (T3-T2)
         ENDIF
*
         IF(PACA.GE.2) THEN
            IF(LC0.GT.0) DEALLOCATE(LUCF)
            DEALLOCATE(LUDF)
         ENDIF
         DEALLOCATE(WORK,CFR,DIAGFR,DIAGQ,CQ)
         DEALLOCATE(NEXT,PREV)
         IF(LTMT) DEALLOCATE(H0,NOM0)
         DEALLOCATE(XSW)
      ENDIF
*----
*  FOR SCR/STIS:
*  VECTORS NORMALIZATION
*----
      IF(LPJJ) THEN
         CALL MCGDS6(NGEFF,NPJJM,NREG,PJJD,V,PJJ)
         CALL MCGDS6(NGEFF,NPJJM,NREG,PJJDX,V,PJJX)
         CALL MCGDS6(NGEFF,NPJJM,NREG,PJJDY,V,PJJY)
         CALL MCGDS6(NGEFF,NPJJM,NREG,PJJDZ,V,PJJZ)
         CALL MCGDS6(NGEFF,NPJJM,NREG,PJJDXI,V,PJJXI)
         CALL MCGDS6(NGEFF,NPJJM,NREG,PJJDYI,V,PJJYI)
         CALL MCGDS6(NGEFF,NPJJM,NREG,PJJDZI,V,PJJZI)
         DEALLOCATE(PJJD,PJJDX,PJJDY,PJJDZ,PJJDXI,PJJDYI,PJJDZI)
*
         DO II=1,NGEFF
            JPSYS=KPSYS(II)
            IF(IPRINT.GT.3) THEN
               CALL PRINAM('PJJ   ',PJJ((II-1)*NREG*NPJJM+1),NREG*NPJJM)
               IF(LPS.GT.0) CALL PRINAM('PSJ   ',PSJ((II-1)*LPS+1),LPS)
            ENDIF
            CALL LCMPUT(JPSYS,'PJJ$MCCG',NREG*NPJJM,2,
     1                  PJJ((II-1)*NREG*NPJJM+1))
            CALL LCMPUT(JPSYS,'PJJX$MCCG',NREG*NPJJM,2,
     1                  PJJX((II-1)*NREG*NPJJM+1))
            CALL LCMPUT(JPSYS,'PJJY$MCCG',NREG*NPJJM,2,
     1                  PJJY((II-1)*NREG*NPJJM+1))
            CALL LCMPUT(JPSYS,'PJJZ$MCCG',NREG*NPJJM,2,
     1                  PJJZ((II-1)*NREG*NPJJM+1))
            CALL LCMPUT(JPSYS,'PJJXI$MCCG',NREG*NPJJM,2,
     1                  PJJXI((II-1)*NREG*NPJJM+1))
            CALL LCMPUT(JPSYS,'PJJYI$MCCG',NREG*NPJJM,2,
     1                  PJJYI((II-1)*NREG*NPJJM+1))
            CALL LCMPUT(JPSYS,'PJJZI$MCCG',NREG*NPJJM,2,
     1                  PJJZI((II-1)*NREG*NPJJM+1))
            IF(LPS.GT.0) THEN
               CALL LCMPUT(JPSYS,'PSJ$MCCG',LPS,2,PSJ((II-1)*LPS+1))
               CALL LCMPUT(JPSYS,'PSJX$MCCG',LPS,2,PSJX((II-1)*LPS+1))
               CALL LCMPUT(JPSYS,'PSJY$MCCG',LPS,2,PSJY((II-1)*LPS+1))
               CALL LCMPUT(JPSYS,'PSJZ$MCCG',LPS,2,PSJZ((II-1)*LPS+1))
            ENDIF
         ENDDO
*
         IF(LPS.GT.0) DEALLOCATE(PSJ,PSJX,PSJY,PSJZ)
         DEALLOCATE(PJJX,PJJXI,PJJY,PJJYI,PJJZ,PJJZI)
         DEALLOCATE(PJJ)
      ENDIF
*
      DEALLOCATE(RHARM,TRHAR)
      DEALLOCATE(KANGL,WZMUA,ZMUA,HH,H,NOM)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(CF,DIAGF)
      RETURN
*
 100  FORMAT('   -->>TIME SPENT IN: ',A40,':',F13.3,' s.')
 200  FORMAT(1X,A6,1X,I4)
      END
