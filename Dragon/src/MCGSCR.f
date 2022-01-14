*DECK MCGSCR
      SUBROUTINE MCGSCR(IPTRK,KPSYS,IPMACR,IPRINT,N1,NG,NGEFF,KPN,K,
     1                  NREG,NANI,NFUNL,M,LPS,KEYFLX,KEYCUR,NZON,NGIND,
     2                  NCONV,FORM,MXSCR,EPSSCR,REBAL,PHIOUT,PHIIN,V,
     3                  NPJJM,KEYANI,IDIR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Acceleration of inner iteration (SCR method).
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
* K       total number of volumes for which specific values
*         of the neutron flux and reactions rates are required.
* NREG    number of volumes.
* NANI    scattering anisotropy (=1 for isotropic scattering).
* NFUNL   number of moments of the flux (in 2D NFUNL=NANI*(NANI+1)/2).
* M       number of material mixtures.
* LPS     dimension of PSJ.
* KEYFLX  position of flux elements in FI vector.
* KEYCUR  position of current elements in FI vector.
* NZON    index-number of the mixture type assigned to each volume.
* NGIND   index of the groups to process.
* NCONV   logical array of convergence status for each group (.TRUE.
*         not converged).
* FORM    input flux format flag (.TRUE. same format as output flux ;
*        .FALSE. same format as input source).
* MXSCR   maximum number of iterations for rebalancing system.
* EPSSCR  convergence criterion for rebalancing system.
* REBAL   type of acceleration (.TRUE. rebalancing ; .FALSE. inner
*         iterations acceleration).
* PHIIN   initial guess (for this iteration) of zonal scalar flux.
* V       volumes.
* NPJJM   second dimension of PJJ.
* KEYANI  'mode to l' index l=KEYANI(nu).
* IDIR    direction of fundamental current for TIBERE with MoC 
*         =0,1,2,3. 
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
      INTEGER N1,NGEFF,NG,IPRINT,KPN,K,NREG,NANI,NFUNL,M,LPS,
     1 KEYFLX(NREG,NFUNL),KEYCUR(*),NZON(K),NGIND(NGEFF),MXSCR,NPJJM,
     2 KEYANI(NFUNL),IDIR
      REAL EPSSCR,PHIIN(KPN,NG),V(N1)
      DOUBLE PRECISION PHIOUT(KPN,NGEFF)
      LOGICAL FORM,NCONV(NGEFF),REBAL
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPMACR,KPMACR,JPSYS
      DOUBLE PRECISION TEMP
      CHARACTER*12 NGTYP
      CHARACTER*12 NAMPJJ,NAMPSJ
      INTEGER, TARGET, SAVE, DIMENSION(1) :: IDUMMY
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NGINDV,NJJ,IJJ,IPOS
      REAL, ALLOCATABLE, DIMENSION(:) ::   XSCAT,MATR
      REAL, ALLOCATABLE, DIMENSION(:) :: PJJ,PSJ
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SC
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: AR,PSI
*
      TYPE(C_PTR) PJJIND_PTR,IS_PTR,JS_PTR
      INTEGER, POINTER, DIMENSION(:) :: IS,JS
      INTEGER, POINTER, DIMENSION(:,:) :: PJJIND
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(NGINDV(NG),SC(0:M,NANI),AR(KPN,NGEFF,2),PSI(KPN,NGEFF,2),
     1 PJJ(NREG*NPJJM),PSJ(LPS))
      CALL XDDSET(PSI,(2*KPN*NGEFF),0.D0)
      CALL LCMGPD(IPTRK,'PJJIND$MCCG',PJJIND_PTR)
      CALL C_F_POINTER(PJJIND_PTR,PJJIND,(/ NPJJM,2 /))
      IF(N1.GT.NREG) THEN
*     recover (IS,JS) arrays
*     IS:  arrays for surfaces neighbors
*     JS:  JS(IS(ISOUT)+1:IS(ISOUT+1)) give the neighboring regions to
*          surface ISOUT.
         CALL LCMGPD(IPTRK,'IS$MCCG',IS_PTR)
         CALL LCMGPD(IPTRK,'JS$MCCG',JS_PTR)
         CALL C_F_POINTER(IS_PTR,IS,(/ N1-NREG+1 /))
         CALL C_F_POINTER(JS_PTR,JS,(/ LPS /))
      ELSE
         IS=>IDUMMY
         JS=>IDUMMY
      ENDIF
      IF(REBAL) THEN
         JPMACR=LCMGID(IPMACR,'GROUP')
         ALLOCATE(NJJ(0:M),IJJ(0:M),IPOS(0:M),XSCAT(0:M*NG))
      ENDIF
      IF(IDIR .EQ.0) THEN
        NAMPJJ='PJJ$MCCG'
        NAMPSJ='PSJ$MCCG'
      ELSEIF(IDIR .EQ. 1) THEN
        NAMPJJ='PJJX$MCCG'
        NAMPSJ='PSJX$MCCG'
      ELSEIF(IDIR .EQ. 2) THEN
        NAMPJJ='PJJY$MCCG'
        NAMPSJ='PSJY$MCCG'
      ELSE
        NAMPJJ='PJJZ$MCCG'
        NAMPSJ='PSJZ$MCCG'
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
*---
*  COMPUTE RESIDUAL OF THE PREVIOUS FREE ITERATION FOR RHS
*---
      DO II=1,NGEFF
         IF(NCONV(II)) THEN
            IG=NGIND(II)
            JPSYS=KPSYS(II)
            CALL LCMGET(JPSYS,'DRAGON-S0XSC',SC(0,1))
            IF(REBAL) THEN
               KPMACR=LCMGIL(JPMACR,IG)
               CALL LCMGET(KPMACR,'NJJS00',NJJ(1))
               CALL LCMGET(KPMACR,'IJJS00',IJJ(1))
               CALL LCMGET(KPMACR,'IPOS00',IPOS(1))
               CALL LCMGET(KPMACR,'SCAT00',XSCAT(1))
            ENDIF
            CALL MCGFCR(IPRINT,IG,II,NG,NGEFF,KPN,N1,NREG,NANI,NFUNL,
     1           M,.FALSE.,KEYFLX,KEYCUR,NZON,NGINDV,FORM,REBAL,PHIOUT,
     2           PHIIN,SC,KEYANI,NJJ,IJJ,IPOS,XSCAT,AR(1,II,1))
         ENDIF
      ENDDO
*---
*  GAUSS SEIDEL ITERATIVE APPROACH TO SOLVE THE REBALANCING SYSTEM
*---
      IF(REBAL) THEN
         NGTYP='GAUSS-SEIDEL' 
         NFIRST=NGEFF+1
         DO II=1,NGEFF
            IF(NCONV(II)) THEN
               IG=NGIND(II)
               KPMACR=LCMGIL(JPMACR,IG)
               CALL LCMGET(KPMACR,'IJJS00',IJJ(1))
               DO IBM=1,M
                  IF(IJJ(IBM).GT.IG) THEN
                     NFIRST=II ! first thermal group index in NGEFF format
                     GOTO 5
                  ENDIF
               ENDDO
            ENDIF
         ENDDO
      ELSE
         NGTYP='   ONE-GROUP'
         NFIRST=1
      ENDIF
 5    CONTINUE
*
      IF(NANI.GT.1) ALLOCATE(MATR(NFUNL*(NFUNL+1)*NREG))
      DO ITSCR=1,MXSCR
         DO 20 II=1,NGEFF
         IF(NCONV(II)) THEN
            IF((II.LT.NFIRST).AND.(ITSCR.GT.1)) GOTO 20
            IG=NGIND(II)
            IGG=IG
            IF(FORM) IGG=II
            JPSYS=KPSYS(II)
            CALL LCMGET(JPSYS,'DRAGON-S0XSC',SC(0,1))
            CALL LCMGET(JPSYS,NAMPJJ,PJJ)
            IF(REBAL) THEN
               KPMACR=LCMGIL(JPMACR,IG)
               CALL LCMGET(KPMACR,'NJJS00',NJJ(1))
               CALL LCMGET(KPMACR,'IJJS00',IJJ(1))
               CALL LCMGET(KPMACR,'IPOS00',IPOS(1))
               CALL LCMGET(KPMACR,'SCAT00',XSCAT(1))
            ENDIF
            DO I=1,NREG
               IBM=NZON(I)
               DO INU=1,NFUNL
                  IND=KEYFLX(I,INU)
                  AR(IND,II,2)=AR(IND,II,1)
               ENDDO
               IF(REBAL) THEN
*              rebalancing option on : contribution from others groups.
                  IF(IBM.GT.0) THEN
                     IND=KEYFLX(I,1)
                     JG=IJJ(IBM)
                     DO 10 JND=1,NJJ(IBM)
                        IF(JG.NE.IG) THEN
                        JJ=NGINDV(JG)
                        IF(JJ.GT.0) THEN
                           AR(IND,II,2)=AR(IND,II,2)+
     1                              XSCAT(IPOS(IBM)+JND-1)*PSI(I,JJ,1)
                        ENDIF
                        ENDIF
                        JG=JG-1
 10                  CONTINUE
                  ENDIF
               ENDIF
            ENDDO
            IF(NANI.EQ.1) THEN
              DO I=1,NREG
                IBM=NZON(I)
                IND=KEYFLX(I,1)
                PSI(IND,II,1)=AR(IND,II,2)
     1                       *PJJ(I)/(1.0-SC(IBM,1)*PJJ(I))
              ENDDO
            ELSE
              CALL MCGSCS(KPN,K,NREG,M,NANI,NFUNL,NPJJM,KEYFLX,KEYANI,
     1             PJJIND,NZON,SC(0,1),PJJ,AR(1,II,2),PSI(1,II,1),MATR)
            ENDIF
         ENDIF
 20      CONTINUE
         IF(REBAL) THEN
            ERRSCR=0.0
            DO II=NFIRST,NGEFF
            IF(NCONV(II)) THEN
               ERR1=0.0
               ERR2=0.0
               DO I=1,NREG
                  DO INU=1,NFUNL
                     IND=KEYFLX(I,INU)
                     TEMP1=REAL(ABS(PSI(IND,II,1)-PSI(IND,II,2)))
                     TEMP2=REAL(ABS(PSI(IND,II,1)))
                     ERR1=MAX(ERR1,TEMP1)
                     ERR2=MAX(ERR2,TEMP2)
                     PSI(IND,II,2)=PSI(IND,II,1)
                  ENDDO
               ENDDO
               IF(ERR2.GT.0.0) ERRSCR=MAX(ERRSCR,ERR1/ERR2)
            ENDIF
            ENDDO
            IF(ERRSCR.LT.EPSSCR) GO TO 30
         ENDIF
      ENDDO
 30   CONTINUE
      IF(NANI.GT.1) DEALLOCATE(MATR)
*
      IF((REBAL).AND.(IPRINT.GT.0)) THEN
         IF(NFIRST.GT.1) WRITE(6,100) NGIND(1),NGIND(NFIRST-1),NGTYP
         IF((MXSCR.GT.1).AND.(NFIRST.LE.NGEFF)) THEN 
            WRITE(6,200) NGTYP,ERRSCR,(ITSCR-1)
         ENDIF
      ELSE
         IF(IPRINT.GT.1) WRITE(6,100) NGIND(1),NGIND(NGEFF),NGTYP
      ENDIF
*----
* PERFORM THE CORRECTION
*----
*     Flux Correction
      DO II=1,NGEFF
      IF(NCONV(II)) THEN
         IG=NGIND(II)
         DO I=1,NREG
            DO INU=1,NFUNL
               IND=KEYFLX(I,INU)
               PHIOUT(IND,II)=PHIOUT(IND,II)+PSI(IND,II,1)
            ENDDO
         ENDDO
      ENDIF
      ENDDO
*     Current Correction
      IF(N1.GT.NREG) THEN
      DO II=1,NGEFF
      IF(NCONV(II)) THEN
         IG=NGIND(II)
         IGG=IG
         IF(FORM) IGG=II
         JPSYS=KPSYS(II)
         CALL LCMGET(JPSYS,NAMPSJ,PSJ)
         CALL LCMGET(JPSYS,'DRAGON-S0XSC',SC(0,1))
         IF(REBAL) THEN
            KPMACR=LCMGIL(JPMACR,IG)
            CALL LCMGET(KPMACR,'NJJS00',NJJ(1))
            CALL LCMGET(KPMACR,'IJJS00',IJJ(1))
            CALL LCMGET(KPMACR,'IPOS00',IPOS(1))
            CALL LCMGET(KPMACR,'SCAT00',XSCAT(1))
         ENDIF
         DO I=NREG+1,N1
         IIS=I-NREG
         INDC=KEYCUR(IIS)
         DO J=IS(IIS)+1,IS(IIS+1)
            MCUI=JS(J)
            IND=KEYFLX(MCUI,1)
            IBM=NZON(MCUI)
            IF(IBM.GT.0) THEN
               TEMP=SC(IBM,1)*(PHIOUT(IND,II)-PHIIN(IND,IGG))
               IF(REBAL) THEN
                  JG=IJJ(IBM)
                  DO 40 JND=1,NJJ(IBM)
                     IF(JG.NE.IG) THEN
                     JJ=NGINDV(JG)
                     IF(JJ.GT.0) THEN
                        JGG=JG
                        IF(FORM) JGG=JJ
                        TEMP=TEMP+XSCAT(IPOS(IBM)+JND-1)
     1                       *(PHIOUT(IND,JJ)-PHIIN(IND,JGG))
                     ENDIF
                     ENDIF
                     JG=JG-1
 40               CONTINUE
               ENDIF
               PHIOUT(INDC,II)=PHIOUT(INDC,II)+PSJ(J)*TEMP/V(I)
            ENDIF
         ENDDO
         ENDDO
      ENDIF
      ENDDO
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      IF(REBAL) DEALLOCATE(XSCAT,IPOS,IJJ,NJJ)
      DEALLOCATE(PSJ,PJJ,PSI,AR,SC,NGINDV)
      RETURN
*
 100  FORMAT(10X,11HSCR: GROUPS,I4,3H TO,I4,2H: ,A12,7H SCHEME)
 200  FORMAT(10X,24HSCR: UP-SCATTE. GROUPS: ,A12,17H ITERATIONS: PRC:,
     1       E9.2,2H (,I4,12H ITERATIONS))
      END
