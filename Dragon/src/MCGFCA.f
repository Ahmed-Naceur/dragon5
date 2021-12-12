*DECK MCGFCA
      SUBROUTINE MCGFCA(IPTRK,KPSYS,IPMACR,IPRINT,N1,NG,NGEFF,KPN,NREG,
     1                  NANI,NFUNL,M,LC,LFORW,PACA,KEYFLX,KEYCUR,NZON,
     2                  NGIND,NCONV,FORM,MXACA,EPSACA,MACFLG,REBFLG,
     3                  PHIOUT,PHIIN,COMBFLG,NPJJM,KEYANI,IDIR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Acceleration of iterations (ACA method).
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
*Parameters: input
* IPTRK   pointer to the tracking LCM object.
* KPSYS   pointer array for each group properties.
* IPMACR  pointer to the macrolib LCM object.
* IPRINT  print parameter (equal to zero for no print).
* N1      number of unknowns per group of the corrective system.
* NG      number of groups.
* NGEFF   number of groups to process.
* KPN     total number of unknowns in vectors SUNKNO and FUNKNO.
* NREG    number of volumes.
* NANI    scattering anisotropy (=1 for isotropic scattering).
* NFUNL   number of moments of the flux (in 2D: NFUNL=NANI*(NANI+1)/2).
* M       number of material mixtures.
* LC      dimension of profiled matrices MCU and CQ.
* LFORW   flag set to .false. to transpose the coefficient matrix.
* PACA    type of preconditioner to solve the ACA corrective system.
* KEYFLX  position of flux elements in FI vector.
* KEYCUR  position of current elements in FI vector.
* NZON    index-number of the mixture type assigned to each volume.
* NGIND   index of the groups to process.
* NCONV   logical array of convergence status for each group (.TRUE.
*         not converged).
* FORM    input flux format flag (.TRUE.  same format as output flux;
*         .FALSE. same format as input source).
* MXACA   maximum number of iterations.
* EPSACA  convergence criterion.
* MACFLG  multigroup cross section flag.
* REBFLG  rebalancing form flag for ACA.
* PHIIN   initial guess (for this iteration) of zonal scalar flux.
* COMBFLG flag for three-step scheme in combination wih SCR.
* NPJJM   second dimension of PJJ.
* KEYANI 'mode to l' index: l=KEYANI(nu).
* IDIR    direction of fundamental current for TIBERE with MoC 
*         (=0,1,2,3). 
*
*Parameters: input/output
* PHIOUT  zonal scalar flux.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,KPSYS(NGEFF),IPMACR
      INTEGER N1,N2,NGEFF,NG,IPRINT,KPN,NREG,NANI,NFUNL,M,LC,PACA,
     1 KEYFLX(NREG,NFUNL),KEYCUR(*),NZON(N1),NGIND(NGEFF),MXACA,NPJJM,
     2 KEYANI(NFUNL),IDIR
      REAL EPSACA,PHIIN(KPN,NG)
      DOUBLE PRECISION PHIOUT(KPN,NGEFF)
      LOGICAL LFORW,FORM,NCONV(NGEFF),MACFLG,REBFLG,COMBFLG
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPMACR,KPMACR,JPSYS
      INTEGER LC0
      REAL FLXN
      CHARACTER*12 NGTYP
      INTEGER, TARGET, SAVE, DIMENSION(1) :: IDUMMY
      REAL, TARGET, SAVE, DIMENSION(1) :: DUMMY
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NGINDV,NJJ,IJJ,IPOS
      REAL, ALLOCATABLE, DIMENSION(:) :: XSCAT
      REAL, ALLOCATABLE, DIMENSION(:,:) :: XSW,PJJ
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: AR,PSI,ARSCR
*
      TYPE(C_PTR) PJJIND_PTR,IM_PTR,MCU_PTR,IPERM_PTR,JU_PTR,IM0_PTR,
     1 MCU0_PTR
      TYPE(C_PTR) DIAGQ_PTR,CQ_PTR,LUDF_PTR,LUCF_PTR,CF_PTR,DIAGF_PTR
      INTEGER, POINTER, DIMENSION(:) :: IM,MCU,IPERM,JU,IM0,MCU0
      INTEGER, POINTER, DIMENSION(:,:) :: PJJIND
      REAL, POINTER, DIMENSION(:) :: DIAGQ,CQ,LUDF,LUCF,CF,DIAGF
*----
*  INITIALIZE POINTERS
*----
      JU=>IDUMMY
      IM0=>IDUMMY
      MCU0=>IDUMMY
      LUDF=>DUMMY
      LUCF=>DUMMY
      CF=>DUMMY
      DIAGF=>DUMMY
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(NGINDV(NG),AR(N1,NGEFF),PSI(N1,NGEFF),XSW(0:M,NANI),
     1 ARSCR(KPN,NGEFF))
      CALL XDRSET(XSW,(M+1)*NANI,0.0)
      CALL XDDSET(AR,N1*NGEFF,0.0D0)
      CALL XDDSET(ARSCR,KPN*NGEFF,0.0D0)
*     recover connection matrices
      CALL LCMGPD(IPTRK,'IM$MCCG',IM_PTR)
      CALL LCMGPD(IPTRK,'MCU$MCCG',MCU_PTR)
      CALL C_F_POINTER(IM_PTR,IM,(/ N1+1 /))
      CALL C_F_POINTER(MCU_PTR,MCU,(/ LC /))
*     recover permutation array
      CALL LCMGPD(IPTRK,'PI$MCCG',IPERM_PTR)
      CALL C_F_POINTER(IPERM_PTR,IPERM,(/ N1 /))
      IF(PACA.GE.2) THEN
         CALL LCMGPD(IPTRK,'JU$MCCG',JU_PTR)
         CALL C_F_POINTER(JU_PTR,JU,(/ N1 /))
      ENDIF
      IF(PACA.EQ.3) THEN
         CALL LCMLEN(IPTRK,'IM0$MCCG',LIM0,ITYLCM)
         CALL LCMLEN(IPTRK,'MCU0$MCCG',LC0,ITYLCM)
         CALL LCMGPD(IPTRK,'IM0$MCCG',IM0_PTR)
         CALL LCMGPD(IPTRK,'MCU0$MCCG',MCU0_PTR)
         CALL C_F_POINTER(IM0_PTR,IM0,(/ LIM0 /))
         CALL C_F_POINTER(MCU0_PTR,MCU0,(/ LC0 /))
      ELSE
         LIM0=0
         LC0=0
      ENDIF
      IF(MACFLG) THEN
         JPMACR=LCMGID(IPMACR,'GROUP')
         ALLOCATE(NJJ(0:M),IJJ(0:M),IPOS(0:M),XSCAT(0:M*NG))
      ENDIF
      IF(REBFLG) THEN
*     N2: number of groups to treat at the same time.
*     rebalancing
         N2=NGEFF      
      ELSE
*     inner iterations acceleration
         N2=1
      ENDIF
*----
*  CONSTRUCT NGINDV (index to pass from "NGEFF format" to "NG format").
*----
      CALL XDISET(NGINDV,NG,0)
      DO II=1,NGEFF
         IF(NCONV(II)) THEN
            IG=NGIND(II)
            NGINDV(IG)=II
         ENDIF
      ENDDO
*----
*  COMPUTE RESIDUAL OF THE PREVIOUS FREE ITERATION FOR RHS WITHOUT
*  COMBFLG OPTION
*----
      IF(IPRINT.GT.10) WRITE(6,*) 'Direction',IDIR
      IF(.NOT.COMBFLG) THEN
         DO II=1,NGEFF
            IF(NCONV(II)) THEN
               IG=NGIND(II)
               JPSYS=KPSYS(II)
               CALL LCMGET(JPSYS,'DRAGON-S0XSC',XSW)
               IF(MACFLG) THEN
                  KPMACR=LCMGIL(JPMACR,IG)
                  CALL LCMGET(KPMACR,'NJJS00',NJJ(1))
                  CALL LCMGET(KPMACR,'IJJS00',IJJ(1))
                  CALL LCMGET(KPMACR,'IPOS00',IPOS(1))
                  CALL LCMGET(KPMACR,'SCAT00',XSCAT(1))
               ENDIF
*              residual for ACA system
               CALL MCGFCR(IPRINT,IG,II,NG,NGEFF,KPN,N1,NREG,NANI,NFUNL,
     1           M,.TRUE.,KEYFLX,KEYCUR,NZON,NGINDV,FORM,MACFLG,PHIOUT,
     2           PHIIN,XSW,IPERM,NJJ,IJJ,IPOS,XSCAT,AR(1,II))
            ENDIF
         ENDDO
*----
*  COMPUTE RESIDUAL OF THE PREVIOUS FREE ITERATION FOR RHS WITH COMBFLG
*  OPTION
*----
      ELSE
         ALLOCATE(PJJ(NREG,NPJJM))
         CALL LCMGPD(IPTRK,'PJJIND$MCCG',PJJIND_PTR)
         CALL C_F_POINTER(PJJIND_PTR,PJJIND,(/ NPJJM,2 /))
         DO II=1,NGEFF
            IF(NCONV(II)) THEN
               IG=NGIND(II)
               JPSYS=KPSYS(II)
               CALL LCMGET(JPSYS,'DRAGON-S0XSC',XSW)
               IF(MACFLG) THEN
                  KPMACR=LCMGIL(JPMACR,IG)
                  CALL LCMGET(KPMACR,'NJJS00',NJJ(1))
                  CALL LCMGET(KPMACR,'IJJS00',IJJ(1))
                  CALL LCMGET(KPMACR,'IPOS00',IPOS(1))
                  CALL LCMGET(KPMACR,'SCAT00',XSCAT(1))
               ENDIF
*              residual for ACA system
               CALL MCGFCR(IPRINT,IG,II,NG,NGEFF,KPN,N1,NREG,NANI,NFUNL,
     1           M,.TRUE.,KEYFLX,KEYCUR,NZON,NGINDV,FORM,MACFLG,PHIOUT,
     2           PHIIN,XSW,IPERM,NJJ,IJJ,IPOS,XSCAT,AR(1,II))
*              residual for SCR-combined scheme
               CALL MCGFCR(IPRINT,IG,II,NG,NGEFF,KPN,N1,NREG,NANI,NFUNL,
     1           M,.FALSE.,KEYFLX,KEYCUR,NZON,NGINDV,FORM,MACFLG,
     2           PHIOUT,PHIIN,XSW,KEYANI,NJJ,IJJ,IPOS,XSCAT,ARSCR(1,II))
               IF(NANI.GT.1) THEN
               IF(IDIR.EQ.0) THEN
                 CALL LCMGET(JPSYS,'PJJ$MCCG',PJJ)
               ELSEIF(IDIR.EQ.1) THEN
                 CALL LCMGET(JPSYS,'PJJX$MCCG',PJJ)
               ELSEIF(IDIR.EQ.2) THEN
                 CALL LCMGET(JPSYS,'PJJY$MCCG',PJJ)
               ELSEIF(IDIR.EQ.3) THEN
                 CALL LCMGET(JPSYS,'PJJZ$MCCG',PJJ)
               ENDIF
               DO I=1,N1
                  J=IPERM(I)
                  IBM=NZON(J)
                  IF(IBM.GE.0) THEN
                     DO IMOD=1,NPJJM
                        INU1=PJJIND(IMOD,1)
                        INU2=PJJIND(IMOD,2)
                        IF((INU1.EQ.1).AND.(INU2.NE.1)) THEN
                           IND2=KEYFLX(J,INU2)
                           AR(I,II)=AR(I,II)+PJJ(J,IMOD)*ARSCR(IND2,II)
                        ELSEIF((INU2.EQ.1).AND.(INU1.NE.1)) THEN
                           IND1=KEYFLX(J,INU1)
                           AR(I,II)=AR(I,II)+PJJ(J,IMOD)*ARSCR(IND1,II)
                        ENDIF
                     ENDDO
                  ENDIF
               ENDDO
               ENDIF
            ENDIF
         ENDDO
         DEALLOCATE(PJJ)
      ENDIF
*---
*  ITERATIVE APPROACH TO SOLVE THE PRECONDITIONING SYSTEM
*---
*     ---
*     GROUP PER GROUP PROCEDURE
*     ---
      IF(MACFLG) THEN
*     MULTIGROUP REBALANCING (GAUSS-SEIDEL SCHEME)
         NGTYP='GAUSS-SEIDEL'
         CALL XDDSET(PSI,N1*NGEFF,0.D0)
         NGFAST=NGEFF
         IF(REBFLG) THEN
*        ONLY FOR FAST GROUPS (thermal group will be treated iteratively)
         DO II=1,NGEFF
            IF(NCONV(II)) THEN
               IG=NGIND(II)
               KPMACR=LCMGIL(JPMACR,IG)
               CALL LCMGET(KPMACR,'IJJS00',IJJ(1))
               DO IBM=1,M
                  IF(IJJ(IBM).GT.IG) THEN
                     NGFAST=II-1 ! last fast group index in NGEFF format
                     GOTO 5
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
         ENDIF
      ELSE
*     INNER ITERATION ACCELERATION
         NGTYP='   ONE-GROUP'
         NGFAST=NGEFF
      ENDIF
 5    CONTINUE
      DO II=1,NGFAST
         IF(NCONV(II)) THEN
*        infinite norm of group scalar flux
            FLXN=0.0
            DO I=1,NREG
               IND=KEYFLX(I,1)
               TEMP=REAL(ABS(PHIOUT(IND,II)))
               FLXN=MAX(TEMP,FLXN)
            ENDDO
            IF(MACFLG) THEN
*           contribution from other groups (Gauss-Seidel multigroup
*           scheme without iterations)
               IG=NGIND(II)
               KPMACR=LCMGIL(JPMACR,IG)
               CALL LCMGET(KPMACR,'NJJS00',NJJ(1))
               CALL LCMGET(KPMACR,'IJJS00',IJJ(1))
               CALL LCMGET(KPMACR,'IPOS00',IPOS(1))
               CALL LCMGET(KPMACR,'SCAT00',XSCAT(1))
               DO I=1,N1
                  J=IPERM(I)
                  IBM=NZON(J)
                  IF(IBM.GT.0) THEN
                     JG=IJJ(IBM)
                     DO 10 JND=1,NJJ(IBM)
                     IF(JG.NE.IG) THEN
                        JJ=NGINDV(JG)
                        IF(JJ.GT.0) THEN
                           AR(I,II)=AR(I,II)+XSCAT(IPOS(IBM)+JND-1)*
     1                              PSI(I,JJ)
                        ENDIF
                     ENDIF
                     JG=JG-1
 10                  CONTINUE
                  ENDIF
               ENDDO
            ENDIF
*           apply preconditioner to RHS            
            IG=NGIND(II)
            JPSYS=KPSYS(II)
            CALL LCMGPD(JPSYS,'DIAGQ$MCCG',DIAGQ_PTR)
            CALL LCMGPD(JPSYS,'CQ$MCCG',CQ_PTR)
            CALL C_F_POINTER(DIAGQ_PTR,DIAGQ,(/ N1 /))
            CALL C_F_POINTER(CQ_PTR,CQ,(/ LC /))
            IF(PACA.GE.2) THEN
               CALL LCMGPD(JPSYS,'ILUDF$MCCG',LUDF_PTR)
               CALL C_F_POINTER(LUDF_PTR,LUDF,(/ N1 /))
               IF(PACA.LT.4) THEN
                  CALL LCMGPD(JPSYS,'ILUCF$MCCG',LUCF_PTR)
                  CALL C_F_POINTER(LUCF_PTR,LUCF,(/ LC /))
               ENDIF
               IF(PACA.GE.3) THEN
                  CALL LCMGPD(JPSYS,'CF$MCCG',CF_PTR)
                  CALL C_F_POINTER(CF_PTR,CF,(/ N1 /))
               ENDIF
            ELSE IF(PACA.EQ.1) THEN
               CALL LCMGPD(JPSYS,'DIAGF$MCCG',DIAGF_PTR)
               CALL C_F_POINTER(DIAGF_PTR,DIAGF,(/ LC /))
            ENDIF
            CALL MCGPRA(LFORW,3,PACA,.TRUE.,N1,LC,IM,MCU,JU,DIAGQ,CQ,
     1           LUDF,LUCF,DIAGF,AR(1,II),PSI(1,II),LC0,IM0,MCU0,CF)
*           group per group BICGSTAB
            JPSYS=KPSYS(II)
            CALL LCMGPD(JPSYS,'DIAGF$MCCG',DIAGF_PTR)
            CALL LCMGPD(JPSYS,'CF$MCCG',CF_PTR)
            CALL C_F_POINTER(DIAGF_PTR,DIAGF,(/ N1 /))
            CALL C_F_POINTER(CF_PTR,CF,(/ LC /))
            IF(PACA.GE.2) THEN
               CALL LCMGPD(JPSYS,'ILUDF$MCCG',LUDF_PTR)
               CALL C_F_POINTER(LUDF_PTR,LUDF,(/ N1 /))
               IF(PACA.LT.4) THEN
                  CALL LCMGPD(JPSYS,'ILUCF$MCCG',LUCF_PTR)
                  CALL C_F_POINTER(LUCF_PTR,LUCF,(/ LC /))
               ENDIF
            ENDIF        
            CALL MCGABG(IPRINT,LFORW,PACA,N1,LC,EPSACA,MXACA,IM,MCU,
     1           JU,DIAGF,CF,LUDF,LUCF,AR(1,II),PSI(1,II),FLXN,LC0,
     2           IM0,MCU0)
         ENDIF
      ENDDO
*
      IF((REBFLG).AND.(IPRINT.GT.0)) THEN
         IF(NGFAST.GT.0) WRITE(6,100) NGIND(1),NGIND(NGFAST),NGTYP
      ELSE
         IF(IPRINT.GT.1) WRITE(6,100) NGIND(1),NGIND(NGFAST),NGTYP
      ENDIF
*
      IF((REBFLG).AND.(NGFAST.LT.NGEFF)) THEN
*     ---
*     MULTIGROUP PROCEDURE
*     ---
*     THERMAL GROUPS REBALANCING
         FLXN=0.0
         NFIRST=NGFAST+1
         DO II=NFIRST,NGEFF
            IF(NCONV(II)) THEN
*              infinite norm of multigroup (thermal groups) scalar flux
               DO I=1,NREG
                  IND=KEYFLX(I,1)
                  TEMP=REAL(ABS(PHIOUT(IND,II)))
                  FLXN=MAX(TEMP,FLXN)
               ENDDO
*              contribution from fast groups to rhs
               IG=NGIND(II)
               KPMACR=LCMGIL(JPMACR,IG)
               CALL LCMGET(KPMACR,'NJJS00',NJJ(1))
               CALL LCMGET(KPMACR,'IJJS00',IJJ(1))
               CALL LCMGET(KPMACR,'IPOS00',IPOS(1))
               CALL LCMGET(KPMACR,'SCAT00',XSCAT(1))
               DO I=1,N1
                  J=IPERM(I)
                  IBM=NZON(J)
                  IF(IBM.GT.0) THEN
                     JG=IJJ(IBM)
                     DO 20 JND=1,NJJ(IBM)
                        IF(JG.NE.IG) THEN
                        JJ=NGINDV(JG)
                        IF((JJ.GT.0).AND.(JJ.LE.NGFAST)) THEN
                           AR(I,II)=AR(I,II)+XSCAT(IPOS(IBM)+JND-1)*
     1                          PSI(I,JJ)
                        ENDIF
                        ENDIF
                        JG=JG-1
 20                  CONTINUE
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
*        apply preconditioner to RHS
         DO II=NFIRST,NGEFF
         IF(NCONV(II)) THEN
            JPSYS=KPSYS(II)
            CALL LCMGPD(JPSYS,'DIAGQ$MCCG',DIAGQ_PTR)
            CALL LCMGPD(JPSYS,'CQ$MCCG',CQ_PTR)
            CALL C_F_POINTER(DIAGQ_PTR,DIAGQ,(/ 1 /))
            CALL C_F_POINTER(CQ_PTR,CQ,(/ 1 /))
            IF(PACA.GE.2) THEN
               CALL LCMGPD(JPSYS,'ILUDF$MCCG',LUDF_PTR)
               CALL C_F_POINTER(LUDF_PTR,LUDF,(/ 1 /))
               IF(PACA.LT.4) THEN
                  CALL LCMGPD(JPSYS,'ILUCF$MCCG',LUCF_PTR)
                  CALL C_F_POINTER(LUCF_PTR,LUCF,(/ 1 /))
               ENDIF
               IF(PACA.GE.3) THEN
                  CALL LCMGPD(JPSYS,'CF$MCCG',CF_PTR)
                  CALL C_F_POINTER(CF_PTR,CF,(/ 1 /))
               ENDIF
            ELSEIF(PACA.EQ.1) THEN
               CALL LCMGPD(JPSYS,'DIAGF$MCCG',DIAGF_PTR)
               CALL C_F_POINTER(DIAGF_PTR,DIAGF,(/ 1 /))
            ENDIF
            CALL MCGPRA(LFORW,3,PACA,.TRUE.,N1,LC,IM,MCU,JU,DIAGQ,CQ,
     1           LUDF,LUCF,DIAGF,AR(1,II),PSI(1,II),LC0,IM0,MCU0,CF)
         ENDIF
         ENDDO
*        multigroup BICGSTAB
         CALL MCGABGR(IPRINT,LFORW,PACA,N1,NG,NFIRST,NGEFF,M,LC,NGIND,
     1        NGINDV,NCONV,KPSYS,JPMACR,NZON,IPERM,IM,MCU,JU,EPSACA,
     2        MXACA,AR,PSI,FLXN,LC0,IM0,MCU0)
      ENDIF
*----
* PERFORM THE CORRECTION
*----
      IF(COMBFLG) THEN
*     -----------------------------------------------
*     ACA is combined in a three-step scheme with SCR
*     -----------------------------------------------
         ALLOCATE(PJJ(NREG,NPJJM))
         CALL LCMGPD(IPTRK,'PJJIND$MCCG',PJJIND_PTR)
         CALL C_F_POINTER(PJJIND_PTR,PJJIND,(/ NPJJM,2 /))
         DO II=1,NGEFF
         IF(NCONV(II)) THEN
            IG=NGIND(II)
            JPSYS=KPSYS(II)
            CALL LCMGET(JPSYS,'DRAGON-S0XSC',XSW)
            IF(MACFLG) THEN
               KPMACR=LCMGIL(JPMACR,IG)
               CALL LCMGET(KPMACR,'NJJS00',NJJ(1))
               CALL LCMGET(KPMACR,'IJJS00',IJJ(1))
               CALL LCMGET(KPMACR,'IPOS00',IPOS(1))
               CALL LCMGET(KPMACR,'SCAT00',XSCAT(1))
            ENDIF
            IF(IDIR.EQ.0) THEN
              CALL LCMGET(JPSYS,'PJJ$MCCG',PJJ)
            ELSEIF(IDIR.EQ.1) THEN
              CALL LCMGET(JPSYS,'PJJX$MCCG',PJJ)
            ELSEIF(IDIR.EQ.2) THEN
              CALL LCMGET(JPSYS,'PJJY$MCCG',PJJ)
            ELSEIF(IDIR.EQ.3) THEN
              CALL LCMGET(JPSYS,'PJJZ$MCCG',PJJ)
            ENDIF
            DO I=1,N1
            J=IPERM(I)
            IBM=NZON(J)
            IF(IBM.GE.0) THEN
*           Flux Correction
               IND=KEYFLX(J,1)
               PHIOUT(IND,II)=PHIOUT(IND,II)
     1                       +(1.0-PJJ(J,1)*XSW(IBM,1))*PSI(I,II)
               DO IMOD=1,NPJJM
                  INU1=PJJIND(IMOD,1)
                  INU2=PJJIND(IMOD,2)
                  IF(INU1.EQ.1) THEN
                     IND2=KEYFLX(J,INU2)
                     PHIOUT(IND,II)=PHIOUT(IND,II)
     1                    -PJJ(J,IMOD)*ARSCR(IND2,II)
                  ELSEIF(INU2.EQ.1) THEN
                     IND1=KEYFLX(J,INU1)
                     PHIOUT(IND,II)=PHIOUT(IND,II)
     1                    -PJJ(J,IMOD)*ARSCR(IND1,II)
                  ENDIF
               ENDDO
               IF(MACFLG) THEN
                  JG=IJJ(IBM)
                  DO 30 JND=1,NJJ(IBM)
                  IF(JG.NE.IG) THEN
                     JJ=NGINDV(JG)
                     IF(JJ.GT.0) THEN
                        PHIOUT(IND,II)=PHIOUT(IND,II)-PJJ(J,1)*
     1                       XSCAT(IPOS(IBM)+JND-1)*PSI(I,JJ)
                     ENDIF
                  ENDIF
                  JG=JG-1
 30               CONTINUE
               ENDIF
            ELSE
*           Current Correction
               IND=KEYCUR(J-NREG)
               PHIOUT(IND,II)=PHIOUT(IND,II)+PSI(I,II)
            ENDIF
            ENDDO
         ENDIF
         ENDDO
         DEALLOCATE(PJJ)
      ELSE
*     -----------------
*     ACA is used alone
*     -----------------
         DO II=1,NGEFF
         IF(NCONV(II)) THEN
            DO I=1,N1
               J=IPERM(I)
               IF(NZON(J).GE.0) THEN
*              Flux Correction
                  IND=KEYFLX(J,1)
                  PHIOUT(IND,II)=PHIOUT(IND,II)+PSI(I,II)
               ELSE
*              Current Correction
                  IND=KEYCUR(J-NREG)
                  PHIOUT(IND,II)=PHIOUT(IND,II)+PSI(I,II)
               ENDIF
            ENDDO
         ENDIF
         ENDDO
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      IF(MACFLG) DEALLOCATE(XSCAT,IPOS,IJJ,NJJ)
      DEALLOCATE(ARSCR,XSW,PSI,AR,NGINDV)
      RETURN
*
 100  FORMAT(10X,11HACA: GROUPS,I4,3H TO,I4,2H: ,A12,7H SCHEME)
      END
