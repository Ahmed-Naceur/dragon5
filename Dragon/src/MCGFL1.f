*DECK MCGFL1
      SUBROUTINE MCGFL1(SUBFFI,SUBFFA,SUBLDC,SUBSCH,CYCLIC,KPSYS,IPRINT,
     1           IPTRK,IFTRAK,IPMACR,NDIM,K,KPN,NLONG,PHIOUT,NZON,
     2           MATALB,M,NANI,NMU,N2MAX,NANGL,NREG,NSOUT,NG,NGEFF,
     3           NGIND,S,IAAC,ISCR,LC,LFORW,PACA,EPSACC,MAXACC,NLIN,
     4           NFUNL,KEYFLX,KEYCUR,QFR,PHIIN,CAZ0,CAZ1,CAZ2,CPO,ZMU,
     5           WZMU,V,SIGAL,LPS,NCONV,FORM,LAST,STIS,NPJJM,REBFLG,
     6           LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,IDIR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform one inner iteration of the characteristics method.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Le Tellier
*
*Parameters: input/output
* SUBFFI  flux integration subroutine with isotropic source.
* SUBFFA  flux integration subroutine with anisotropic source.
* SUBLDC  flux integration subroutine with linear-discontinuous source.
* SUBSCH  track coefficients calculation subroutine.
* CYCLIC  cyclic tracking flag.
* KPSYS   pointer array for each group properties.
* IPRINT  print parameter (equal to zero for no print).
* IPTRK   pointer to the tracking (L_TRACK signature).
* IFTRAK  tracking file unit number.
* IPMACR  pointer to the macrolib LCM object.
* NDIM    number of dimensions for the geometry.
* K       total number of volumes for which specific values
*         of the neutron flux and reactions rates are required.
* KPN     total number of unknowns per group in vector PHIIN.
* NLONG   number of spatial unknowns.
* PHIOUT  output flux vector.
* NZON    mixture-albedo index array in MCCG format.
* MATALB  albedo-mixture index array in MOCC format.
* M       number of material mixtures.
* NANI    scattering anisotropy (=1 for isotropic scattering).
* NMU     order of the polar quadrature in 2D / 1 in 3D.
* N2MAX   maximum number of elements in a track.
* NANGL   number of tracking angles in the plane.
* NREG    number of volumes.
* NSOUT   number of outer surfaces.
* NG      number of groups.
* NGEFF   number of groups to process.
* NGIND   index of the groups to process.
* S       scratch.
* IAAC    no acceleration / CDD acceleration of inner iterations (0/1).
* ISCR    no acceleration / SCR acceleration of inner iterations (0/1).
* LC      dimension of profiled matrices MCU and CQ.
* LFORW   flag set to .false. to transpose the coefficient matrix.
* PACA    type of preconditioner to solve the ACA corrective system.
* EPSACC  stopping criterion for BICGSTAB in ACA resolution.
* MAXACC  maximum number of iterations allowed for BICGSTAB in ACA
*         resolution.
* NLIN    number of polynomial components in flux spatial expansion.
* NFUNL   number of moments of the flux (in 2D: NFUNL=NANI*(NANI+1)/2).
* KEYFLX  position of flux elements in PHIIN vector.
* KEYCUR  position of current elements in PHIIN vector.
* QFR     input source vector.
* PHIIN   input flux vector.
* CAZ0    cosines of the tracking polar angles in 3D.
* CAZ1    first cosines of the different tracking azimuthal angles.
* CAZ2    second cosines of the different tracking azimuthal angles.
* CPO     cosines of the different tracking polar angles in 2D.
* ZMU     polar quadrature set in 2D.
* WZMU    polar quadrature set in 2D.
* V       volumes.
* SIGAL   total cross-section and albedo array.
* LPS     used in scr acceleration.
* NCONV   logical array of convergence status for each group (.TRUE.
*         not converged).
* FORM    input flux format flag (.TRUE. same format as output flux ;
*         .FALSE. same format as input source).
* LAST    flag for SCR and ACA rebalancing.
* STIS    source term isolation' option for flux integration.
* NPJJM   number of pjj modes to store for STIS option.
* REBFLG  ACA or SCR rebalancing flag.
* LPRISM  3D prismatic extended tracking flag.
* N2REG   number of regions in the 2D tracking if LPRISM.
* N2SOU   number of external surfaces in the 2D tracking if LPRISM.
* NZP     number of z-plans if LPRISM.
* DELU    input track spacing for 3D track reconstruction if LPRISM.
* FACSYM  tracking symmetry factor for maximum track length if LPRISM.
* IDIR    direction of fundamental current for TIBERE with MoC 
*         (=0,1,2,3). 
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*---
* SUBROUTINES ARGUMENTS
*---
      TYPE(C_PTR) KPSYS(NGEFF),IPTRK,IPMACR
      INTEGER NGEFF,IPRINT,IFTRAK,NDIM,K,KPN,NLONG,NG,NGIND(NGEFF),
     1 NZON(NLONG),M,NANI,NMU,N2MAX,NANGL,IAAC,LC,PACA,NREG,NSOUT,
     2 ISCR,LPS,NLIN,NFUNL,KEYFLX(NREG,NLIN,NFUNL),KEYCUR(NLONG-NREG),
     3 MAXACC,STIS,NPJJM,MATALB(-NSOUT:NREG),N2REG,N2SOU,NZP,IDIR
      REAL QFR(KPN,NG),PHIIN(KPN,*),CPO(NMU),ZMU(NMU),WZMU(NMU),
     1 V(NLONG),SIGAL(-6:M,NGEFF),EPSACC,DELU,FACSYM
      DOUBLE PRECISION CAZ0(NANGL),CAZ1(NANGL),CAZ2(NANGL),
     1 PHIOUT(KPN,NGEFF),S(KPN,NGEFF)
      LOGICAL LFORW,CYCLIC,NCONV(NGEFF),FORM,LAST,REBFLG,LPRISM
      EXTERNAL SUBFFI,SUBFFA,SUBLDC,SUBSCH
*---
* LOCAL VARIABLES
*---
      TYPE(C_PTR) JPSYS
      REAL T1,T2,T3
      INTEGER NCODE(6),SSYM
      INTEGER, DIMENSION(1), TARGET :: IDUMMY
      CHARACTER TEXT4*4
      LOGICAL MACFLG,COMBFLG
      INTEGER ICREB,ITYLCM
*----
*  ALLOCATABLE ARRAYS
*----
      TYPE(C_PTR) IBC_PTR,XSSC_PTR,PJJM_PTR
      TYPE(C_PTR) INDREG_PTR,Z_PTR,VNORF_PTR,CMU_PTR,CMUI_PTR,SMU_PTR,
     1 SMUI_PTR,TMU_PTR,TMUI_PTR
      INTEGER, ALLOCATABLE, DIMENSION(:) ::KEYANI
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISGNR
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: XSIXYZ
      INTEGER, POINTER, DIMENSION(:) :: IBC,INDREG,PJJM
      REAL, POINTER, DIMENSION(:) :: XSSC,Z
      DOUBLE PRECISION, POINTER, DIMENSION(:) :: VNORF,CMU,CMUI,SMU,
     1 SMUI,TMU,TMUI
*
      CALL XDDSET(PHIOUT,KPN*NGEFF,0.D0)
*
      REWIND IFTRAK
      READ(IFTRAK) TEXT4,NCOMNT,NBTR,IFMT
      DO ICOM=1,NCOMNT
         READ(IFTRAK)
      ENDDO
      READ(IFTRAK) (NITMA, II=1,7),MXSUB,NITMA
      DO ICOM=1,6
         READ(IFTRAK)
      ENDDO
*
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
         NDIMB=3
         NMAX=(INT(FACSYM)+1)*N2MAX*(NZP+2)
      ELSE
         NDIMB=NDIM
         NMAX=N2MAX
      ENDIF
*----
*  SOURCE ELEMENTS CALCULATION
*----
      IF(NLONG-NREG.GT.0) THEN
        CALL LCMGPD(IPTRK,'BC-REFL+TRAN',IBC_PTR)
        CALL C_F_POINTER(IBC_PTR,IBC,(/ NLONG-NREG /))
      ELSE
        IBC => IDUMMY
      ENDIF
      DO II=1,NGEFF
         IF(NCONV(II)) THEN
           IG=NGIND(II)
           IGG=IG
           IF(FORM) IGG=II
           JPSYS=KPSYS(II)
           CALL LCMGPD(JPSYS,'DRAGON-S0XSC',XSSC_PTR)
           CALL C_F_POINTER(XSSC_PTR,XSSC,(/ (M+1)*NANI /))
           CALL MCGFCS(NLONG,NDIMB,NZON,QFR(1,IG),PHIIN(1,IGG),M,NANI,
     1        NLIN,NFUNL,XSSC,S(1,II),KPN,NREG,IPRINT,KEYFLX,KEYCUR,
     2        IBC,SIGAL(-6,II),STIS)
         ENDIF
      ENDDO
*---
* GENERATE ALL SIGNS FOR SPHERICAL HARMONICS
*---
      NANIX=NANI
      IF(NLIN.EQ.3) NANIX=MAX(NANI,3)
      IF(NDIM.EQ.1) THEN
         NFUNLX=NANIX
         NMOD=2
      ELSE IF((.NOT.LPRISM).AND.(NDIM.EQ.2)) THEN
         NFUNLX=NANIX*(NANIX+1)/2
         NMOD=4
      ELSE ! NDIM.EQ.3
         NFUNLX=NANIX*NANIX
         NMOD=8
      ENDIF
      ALLOCATE(ISGNR(NMOD,NFUNLX),KEYANI(NFUNLX))
      CALL MOCIK3(NANIX-1,NFUNLX,NMOD,ISGNR,KEYANI)
*----
*  FLUX INTEGRATION UPON THE TRACKING FILE
*----
      CALL KDRCPU(T1)
      IF(CYCLIC) THEN
*     --------------------------------
*     Method of Cyclic Characteristics
*     --------------------------------
         CALL MOCFCF(SUBFFI,SUBFFA,SUBLDC,SUBSCH,IFTRAK,NBTR,MXSUB,
     1        N2MAX,KPN,NREG,NSOUT,M,6,NGEFF,NANGL,NMU,NANI,NFUNL,
     2        NMOD,NANIX,NLIN,NFUNLX,KEYFLX,MATALB,NCONV,SIGAL,CAZ1,
     3        CAZ2,CPO,ZMU,WZMU,PHIOUT,S,ISGNR,IDIR)
      ELSE
*     ------------------------------------
*     Method of Non-Cyclic Characteristics
*     ------------------------------------
         IF(LPRISM) THEN
*        3D PRISMATIC GEOMETRY CONSTRUCTED FROM A 2D TRACKING
            CALL LCMSIX(IPTRK,'PROJECTION',1)
            CALL LCMGPD(IPTRK,'IND2T3',INDREG_PTR)
            CALL LCMGPD(IPTRK,'ZCOORD',Z_PTR)
            CALL LCMGPD(IPTRK,'VNORF',VNORF_PTR)
            CALL LCMGPD(IPTRK,'CMU',CMU_PTR)
            CALL LCMGPD(IPTRK,'CMUI',CMUI_PTR)
            CALL LCMGPD(IPTRK,'SMU',SMU_PTR)
            CALL LCMGPD(IPTRK,'SMUI',SMUI_PTR)
            CALL LCMGPD(IPTRK,'TMU',TMU_PTR)
            CALL LCMGPD(IPTRK,'TMUI',TMUI_PTR)
            CALL LCMSIX(IPTRK,'PROJECTION',2)
*
            CALL C_F_POINTER(INDREG_PTR,INDREG,
     1                       (/ (N2REG+N2SOU+1)*(NZP+2) /))
            CALL C_F_POINTER(Z_PTR,Z,(/ NZP+1 /))
            CALL C_F_POINTER(VNORF_PTR,VNORF,(/ NREG*NANGL*NMU*2 /))
            CALL C_F_POINTER(CMU_PTR,CMU,(/ NMU /))
            CALL C_F_POINTER(CMUI_PTR,CMUI,(/ NMU /))
            CALL C_F_POINTER(SMU_PTR,SMU,(/ NMU /))
            CALL C_F_POINTER(SMUI_PTR,SMUI,(/ NMU /))
            CALL C_F_POINTER(TMU_PTR,TMU,(/ NMU /))
            CALL C_F_POINTER(TMUI_PTR,TMUI,(/ NMU /))
            CALL MCGPTF(SUBFFI,SUBFFA,SUBSCH,IFTRAK,NBTR,N2MAX,KPN,
     1           K,NREG,M,NGEFF,NANGL,NMU,NANI,NFUNL,NMOD,KEYFLX,
     2           KEYCUR,NZON,NCONV,CAZ1,CAZ2,CPO,WZMU,PHIOUT,S,SIGAL,
     3           ISGNR,NMAX,NZP,N2REG,N2SOU,DELU,INDREG,Z,VNORF,CMU,
     4           CMUI,SMU,SMUI,TMU,TMUI,SSYM,IDIR)
         ELSE
*        REGULAR 2D OR 3D GEOMETRY
            ALLOCATE(XSIXYZ(NSOUT,3))
            CALL XDDSET(XSIXYZ,3*NSOUT,0.D0)
            CALL LCMLEN(IPTRK,'XSI$MCCG',ICREB,ITYLCM)
            IF(ICREB.EQ.3*NSOUT) CALL LCMGET(IPTRK,'XSI$MCCG',XSIXYZ)
            CALL MCGFCF(SUBFFI,SUBFFA,SUBLDC,SUBSCH,IFTRAK,NBTR,N2MAX,
     1           NDIM,KPN,K,NREG,M,NGEFF,NANGL,NMU,NANI,NFUNL,NMOD,
     2           NANIX,NLIN,NFUNLX,KEYFLX,KEYCUR,NZON,NCONV,CAZ0,CAZ1,
     3           CAZ2,CPO,ZMU,WZMU,PHIOUT,S,SIGAL,ISGNR,IDIR,NSOUT,
     4           XSIXYZ(1,IDIR))
            DEALLOCATE(XSIXYZ)
         ENDIF
      ENDIF
*
      IF(STIS.EQ.1) THEN
         CALL LCMGPD(IPTRK,'PJJIND$MCCG',PJJM_PTR)
         CALL C_F_POINTER(PJJM_PTR,PJJM,(/ NPJJM*2 /))
         CALL MCGFST(NGEFF,KPSYS,NCONV,KPN,NLONG,NREG,NANI,NFUNL,NPJJM,
     1        KEYFLX,KEYCUR,PJJM,NZON,V,S,PHIOUT,IDIR)
      ELSEIF(STIS.EQ.-1) THEN
         DO II=1,NGEFF
         IF(NCONV(II)) THEN
            CALL MCGFMC(KPN,NLONG,NREG,M,NANI,NFUNL,NZON,KEYFLX,KEYCUR,
     1           PHIOUT(1,II),V,S(1,II),SIGAL(0,II),KEYANI)
         ENDIF
         ENDDO
      ELSE
         DO II=1,NGEFF
         IF(NCONV(II)) THEN
            DO I=1,NLONG
            IF(V(I).GT.0.) THEN 
               IF(NZON(I).LT.0) THEN
                  IND=KEYCUR(I-NREG)
                  PHIOUT(IND,II)=PHIOUT(IND,II)/V(I)
               ELSE
                  DO IL=1,NFUNL
                     DO IU=1,NLIN
                        IND=KEYFLX(I,IU,IL)
                        PHIOUT(IND,II)=PHIOUT(IND,II)/V(I)
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF
            ENDDO
         ENDIF
         ENDDO
      ENDIF
      CALL KDRCPU(T2)
*----
*  PRECONDITIONING TECHNIQUES
*----
      MACFLG=((LAST).AND.(C_ASSOCIATED(IPMACR)).AND.((IAAC.GT.1).OR.
     1       (ISCR.GT.1)))
      REBFLG=(REBFLG.AND.MACFLG)
      COMBFLG=((IAAC.GT.0).AND.(ISCR.GT.0))
      IF(IAAC.GT.0) THEN
*     ---------------------------------
*     Algebraic Collapsing Acceleration
*     ---------------------------------
         IF(REBFLG) MAXACC=IAAC
         CALL MCGFCA(IPTRK,KPSYS,IPMACR,IPRINT,NLONG,NG,NGEFF,KPN,
     1        NREG,NANI,NFUNL,M,LC,LFORW,PACA,KEYFLX,KEYCUR,NZON,NGIND,
     2        NCONV,FORM,MAXACC,EPSACC,MACFLG,REBFLG,PHIOUT,PHIIN,
     3        COMBFLG,NPJJM,KEYANI,IDIR)
         CALL KDRCPU(T3)
*     ---------------------------------
      ENDIF
      NLON2=NLONG
      IF(COMBFLG) NLON2=NREG
      IF(ISCR.GT.0) THEN
*     ---------------------------------
*     Self-Collision Rebalancing Method
*     ---------------------------------
         IF(REBFLG) THEN
*        rebalancing            
            MAXSCR=ISCR
         ELSE
            MAXSCR=1
         ENDIF
         CALL MCGSCR(IPTRK,KPSYS,IPMACR,IPRINT,NLON2,NG,NGEFF,KPN,K,
     1        NREG,NANI,NFUNL,M,LPS,KEYFLX,KEYCUR,NZON,NGIND,NCONV,
     2        FORM,MAXSCR,EPSACC,MACFLG,PHIOUT,PHIIN,V,NPJJM,KEYANI,
     3        IDIR)
         CALL KDRCPU(T3)
*     ---------------------------------
      ENDIF
      DEALLOCATE(KEYANI,ISGNR)
      IF((IPRINT.GT.1).AND.(LAST)) THEN
         WRITE(6,100) ' FLUX INTEGRATION       ',(T2-T1)
         IF((IAAC.GT.0).OR.(ISCR.GT.0)) THEN
            WRITE(6,100) ' ACCELERATION           ',(T3-T2)
         ENDIF
      ENDIF
      RETURN
*
 100  FORMAT('   -->>TIME SPENT IN ',A24,':',F13.3)
      END
