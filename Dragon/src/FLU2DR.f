*DECK FLU2DR
      SUBROUTINE FLU2DR(IPRT,IPMACR,IPFLUX,IPSYS,IPTRK,IPFLUP,IPSOU,
     1 IGPT,IFTRAK,CXDOOR,TITLE,NUNKNO,NREG,NSOUT,NANIS,NLF,NLIN,NFUNL,
     2 NGRP,NMAT,NIFIS,LFORW,LEAKSW,MAXINR,EPSINR,MAXOUT,EPSUNK,EPSOUT,
     3 NCPTL,NCPTA,ITYPEC,IPHASE,ITPIJ,ILEAK,OPTION,REFKEF,MATCOD,
     4 KEYFLX,VOL,XSTOT,XSTRC,XSDIA,XSNUF,XSCHI,LREBAL,INITFL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Fixed source problem or inverse power method for K-effective or
* buckling iteration. Perform thermal iterations.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* IPRT    print flag.
* IPMACR  pointer to the macrolib LCM object.
* IPFLUX  pointer to the flux LCM object.
* IPSYS   pointer to the system LCM object.
* IPTRK   pointer to the tracking LCM object.
* IPFLUP  pointer to the unperturbed flux LCM object (if ITYPEC=1).
* IPSOU   pointer to the fixed source LCM object (if ITYPEC=0 or 1).
* IGPT    index of the fixed source eigenvalue problem to solve.
* IFTRAK  tracking file unit number.
* CXDOOR  character name of the flux solution door.
* TITLE   title.
* NUNKNO  number of unknowns per energy group including spherical
*         harmonic terms, interface currents and fundamental
*         currents.
* NREG    number of regions.
* NSOUT   number of outer surfaces.
* NANIS   maximum cross section Legendre order.
* NLF     number of Legendre orders for the flux.
* NLIN    number of polynomial components in flux spatial expansion.
* NFUNL   number of spherical harmonics components.
* NGRP    number of energy groups.
* NMAT    number of mixtures.
* NIFIS   number of fissile isotopes.
* LFORW   flag set to .false. to solve an adjoint problem.
* LEAKSW  leakage flag (=.true. if leakage is present on the outer
*         surface).
* MAXINR  maximum number of thermal iterations.
* EPSINR  thermal iterations epsilon.
* MAXOUT  maximum number of outer iterations.
* EPSUNK  outer iterations eigenvector epsilon.
* EPSOUT  outer iterations eigenvalue epsilon.
* NCPTL   number of free iterations in an acceleration cycle.
* NCPTA   number of accelerated iterations in an acceleration cycle.
* ITYPEC  type of flux evaluation:
*         =-1 skip the flux calculation;
*         =0 fixed sources;
*         =1 fixed source eigenvalue problem (GPT type);
*         =2 fission sources/k effective convergence;
*         =3 fission sources/k effective convergence/
*             db2 buckling evaluation;
*         =4 fission sources/db2 buckling convergence;
*         =5 b2 sources/db2 buckling convergence.
* IPHASE  type of flux solution door (1 for asm 2 for pij).
* ITPIJ   type of cp available:
*         =1 scatt mod pij (wij);
*         =2 stand. pij;
*         =3 scatt mod pij+pijk (wij,wijk);
*         =4 standard pij+pijk.
* ILEAK   method used to include DB2 effect:
*         =1 the scattering modified cp matrix is multiplied by PNLR;
*         =2 the reduced cp matrix is multiplied by PNL;
*         =3 sigs0-db2 approximation;
*         =4 albedo approximation;
*         =5 Ecco-type isotropic streaming model;
*         >5 Tibere type anisotropic streaming model.
* OPTION  type of leakage coefficients:
*         'LKRD' (recover leakage coefficients in Macrolib);
*         'RHS' (recover leakage coefficients in RHS flux object);
*         'B0' (B-0), 'P0' (P-0), 'B1' (B-1),
*         'P1' (P-1), 'B0TR' (B-0 with transport correction) or 'P0TR'
*         (P-0 with transport correction).
* REFKEF  target effective multiplication factor.
* MATCOD  mixture indices.
* KEYFLX  index of L-th order flux components in unknown vector.
* VOL     volumes.
* XSTOT   non transport-corrected macroscopic total cross sections.
* XSTRC   transport-corrected macroscopic total cross sections.
* XSDIA   transport-corrected macroscopic within-group scattering cross
*         sections.
* XSNUF   nu*macroscopic fission cross sections.
* XSCHI   fission spectrum.
* LREBAL  thermal iteration rebalancing flag (=.true. if thermal
*         rebalancing required).
* INITFL  flux initialization flag (=0/1/2: uniform flux/LCM/DSA).
*
*-----------------------------------------------------------------------
*
*----
*  INTERNAL PARAMETERS:
*   SYBILF : SYBIL FLUX SOLUTION DOOR                 EXT ROUTINE
*   TRFICF : DEFAULT CP FLUX SOLUTION DOOR            EXT ROUTINE
*   BIVAF  : DEFAULT 2D DIFFUSION FLUX SOLUTION DOOR  EXT ROUTINE
*   TRIVAF : DEFAULT 3D DIFFUSION FLUX SOLUTION DOOR  EXT ROUTINE
*   PNF    : DEFAULT PN/SPN FLUX SOLUTION DOOR        EXT ROUTINE
*   SNF    : DEFAULT SN FLUX SOLUTION DOOR            EXT ROUTINE
*----
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMACR,IPFLUX,IPSYS,IPTRK,IPFLUP,IPSOU
      INTEGER IPRT,IGPT,IFTRAK,NUNKNO,NREG,NSOUT,NANIS,NLF,NLIN,NFUNL,
     1 NGRP,NMAT,NIFIS,MAXINR,MAXOUT,NCPTL,NCPTA,ITYPEC,IPHASE,ITPIJ,
     2 ILEAK,MATCOD(NREG),KEYFLX(NREG,NLIN,NFUNL),INITFL
      REAL EPSINR,EPSUNK,EPSOUT,VOL(NREG),XSTOT(0:NMAT,NGRP),
     1 XSTRC(0:NMAT,NGRP),XSDIA(0:NMAT,0:NANIS,NGRP),
     2 XSNUF(0:NMAT,NIFIS,NGRP),XSCHI(0:NMAT,NIFIS,NGRP),B2(4)
      CHARACTER CXDOOR*12,TITLE*72,OPTION*4
      LOGICAL LFORW,LEAKSW,LREBAL,CFLI,CEXE
      DOUBLE PRECISION REFKEF,XCSOU(NGRP)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40)
      TYPE(C_PTR) IPREB,J1,JPSOU,JPFLUX,JPMACR,KPMACR,JPSYS,KPSYS,IPSTR,
     1 JPSTR,KPSTR,JPFLUP1,JPFLUP2
      INTEGER JPAR(NSTATE),KEYSPN(NREG)
      CHARACTER CAN(0:19)*2,MESSIN*5,MESSOU*5,HTYPE(0:5)*4
      INTEGER INDD(3)
      DOUBLE PRECISION AKEEP(8),FISOUR,OLDBIL,AKEFF,AKEFFO,AFLNOR,
     1 BFLNOR,DDELN1,DDELD1
      LOGICAL LSCAL,LEXAC
      REAL ALBEDO(6), FLUXC(NREG)
*
************************************************************************
*                                                                      *
*   ICHAR       : COUNTER FOR NUM. OF OUTER ITERATIONS                 *
*   ICTOT       : TOTAL NUMBER OF FLUX CALCULATIONS                    *
*                                                                      *
************************************************************************
*
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,NJJ,IPOS,NPSYS,KEYCUR,
     1 MATALB
      REAL, ALLOCATABLE, DIMENSION(:) :: FXSOR,DIFF,XSCAT,GAMMA,V,FL
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: FLUX
      LOGICAL REBFLG
*----
*  DATA STATEMENTS
*----
      SAVE CAN,HTYPE
      DATA CAN /'00','01','02','03','04','05','06','07','08','09',
     >          '10','11','12','13','14','15','16','17','18','19'/
      DATA (HTYPE(JJ),JJ=0,5)/'S   ','P   ',2*'K   ','B   ','L   '/
*----
*  SCRATCH STORAGE ALLOCATION
*   DIFF    homogeneous leakage coefficients.
*   FLUX    iteration flux:
*           FLUX(:,:,1) <=old     outer;
*           FLUX(:,:,2) <=present outer;
*           FLUX(:,:,3) <=new     outer;
*           FLUX(:,:,4) <=source  outer;
*           FLUX(:,:,5) <=old     inner;
*           FLUX(:,:,6) <=present inner;
*           FLUX(:,:,7) <=new     inner;
*           FLUX(:,:,8) <=source  inner.
*----
      ALLOCATE(IJJ(0:NMAT),NJJ(0:NMAT),IPOS(0:NMAT),NPSYS(NGRP))
      ALLOCATE(DIFF(NGRP),FLUX(NUNKNO,NGRP,8),XSCAT(0:NMAT*NGRP),
     1 GAMMA(NGRP))
*
      REBFLG=.TRUE.
      IPREB=IPMACR
*
      ICHAR=0
      ICTOT=0
*----
*  RECOVER INDEX FOR THE CURRENTS IN FLUX, NUMERICAL SURFACES,
*  ALBEDO IF NEEDED BY THE REBALANCING.
*  RECOVER THE SPATIAL APPROXIMATION ORDER (IELEM).
*----
      IELEM=1
      ICREB=0
      NNN=0
      INSB=0
      IBFP=0
      IF(CXDOOR.EQ.'MCCG') THEN
         CALL LCMLEN(IPTRK,'KEYCUR$MCCG',ICREB,ITYLCM)
         IF(ICREB.GT.0) THEN
            CALL LCMLEN(IPTRK,'NZON$MCCG',ILONG,ITYLCM)
            NNN=ILONG-ICREB
            ALLOCATE(KEYCUR(NSOUT),V(ILONG),MATALB(ILONG))
            CALL LCMGET(IPTRK,'KEYCUR$MCCG',KEYCUR)
            CALL LCMGET(IPTRK,'V$MCCG',V)
            CALL LCMGET(IPTRK,'ALBEDO',ALBEDO)
            CALL LCMGET(IPTRK,'NZON$MCCG',MATALB)
         ENDIF
      ELSE IF(CXDOOR.EQ.'SN') THEN
         CALL LCMGET(IPTRK,'STATE-VECTOR',JPAR)
         ITYPE=JPAR(6)
         IELEM=JPAR(8)
         IF(ITYPE.EQ.5) IELEM=IELEM**2
         IF(ITYPE.EQ.7) IELEM=IELEM**3
         INSB=JPAR(27)
         IBFP=JPAR(31)
      ELSE IF((CXDOOR.EQ.'BIVAC').OR.(CXDOOR.EQ.'TRIVAC')) THEN
         IF(IELEM.GT.1) CALL XABORT('FLU2DR: ONLY IELEM=1 AVAILABLE.')
      ENDIF
*----
*  SELECT THE CALCULATION DOORS FOR WHICH A GROUP-BY-GROUP SCALAR
*  PROCEDURE WILL BE USED. A VECTORIAL APPROACH WILL BE USED WITH
*  OTHER DOORS.
*----
      LSCAL=(CXDOOR.EQ.'TRAFIC').OR.(CXDOOR.EQ.'SYBIL').OR.
     > (CXDOOR.EQ.'BIVAC').OR.(CXDOOR.EQ.'TRIVAC').OR.(CXDOOR.EQ.'SN')
      IF(INSB.EQ.1) LSCAL=.FALSE.
*
      CALL KDRCPU(CPU0)
      IF(ILEAK.LT.5) THEN
        INORM=1
      ELSE IF(ILEAK.EQ.5) THEN
        INORM=2
      ELSE IF(ILEAK.GT.5) THEN
        INORM=3
      ENDIF
      LEXAC=.FALSE.
      AKEEP(5)=1.0D0
      AKEEP(6)=1.0D0
      AKEEP(7)=1.0D0
      CALL XDRSET(DIFF,NGRP,0.0)
      CALL XDRSET(GAMMA,NGRP,1.0)
*----
*  EXTERNAL FLUX(:,:,2) INITIALISATION AND FIXED-EXTERNAL SOURCE IN
*  FLUX(:,:,4)
*----
      IF(ITYPEC.GE.3) THEN
         CALL LCMGET(IPFLUX,'B2  B1HOM',B2(4))
         IF(ILEAK.GE.6) CALL LCMGET(IPFLUX,'B2  HETE',B2)
      ELSE
         CALL XDRSET(B2,4,0.0)
      ENDIF
      AKEFFO=0.0D0
      IF(LFORW) THEN
         IF(ITYPEC.EQ.1) THEN
           J1=LCMGID(IPSOU,'DSOUR')
           JPSOU=LCMGIL(J1,IGPT)
           J1=LCMGID(IPFLUX,'DFLUX')
           JPFLUX=LCMGIL(J1,IGPT)
         ELSE
           IF(C_ASSOCIATED(IPSOU)) THEN
             J1=LCMGID(IPSOU,'DSOUR')
             JPSOU=LCMGIL(J1,1)
           ENDIF
           JPFLUX=LCMGID(IPFLUX,'FLUX')
         ENDIF
      ELSE
         IF(ITYPEC.EQ.1) THEN
           J1=LCMGID(IPSOU,'ASOUR')
           JPSOU=LCMGIL(J1,IGPT)
           J1=LCMGID(IPFLUX,'ADFLUX')
           JPFLUX=LCMGIL(J1,IGPT)
         ELSE
           IF(C_ASSOCIATED(IPSOU)) THEN
             J1=LCMGID(IPSOU,'ASOUR')
             JPSOU=LCMGIL(J1,1)
           ENDIF
           JPFLUX=LCMGID(IPFLUX,'AFLUX')
         ENDIF
      ENDIF
      ALLOCATE(FXSOR(0:NMAT))
      DO 20 IG=1,NGRP
      JPMACR=LCMGID(IPMACR,'GROUP')
      CALL XDRSET(FLUX(1,IG,2),NUNKNO,0.0)
      CALL XDRSET(FLUX(1,IG,4),NUNKNO,0.0)
      CALL LCMLEL(JPFLUX,1,ILINIT,ITYLCM)
      IF(LFORW) THEN
         CALL LCMGDL(JPFLUX,IG,FLUX(1,IG,2))
      ELSE
         CALL LCMGDL(JPFLUX,NGRP-IG+1,FLUX(1,IG,2))
      ENDIF
      IF((ITYPEC.EQ.0).AND.(.NOT.C_ASSOCIATED(IPSOU))) THEN
         KPMACR=LCMGIL(JPMACR,IG)
         FXSOR(0)=0.0
         CALL LCMGET(KPMACR,'FIXE',FXSOR(1))
         DO 15 IE=1,NLIN
         DO 10 IR=1,NREG
         IND=KEYFLX(IR,IE,1)
         IF(IND.NE.0) FLUX(IND,IG,4)=FXSOR(MATCOD(IR))
   10    CONTINUE
   15    CONTINUE
      ELSE IF(((ITYPEC.EQ.0).OR.(ITYPEC.EQ.1)).AND.C_ASSOCIATED(IPSOU))
     1 THEN
         IF(LFORW) THEN
            CALL LCMGDL(JPSOU,IG,FLUX(1,IG,4))
         ELSE
            CALL LCMGDL(JPSOU,NGRP-IG+1,FLUX(1,IG,4))
         ENDIF
      ENDIF
   20 CONTINUE
      DEALLOCATE(FXSOR)
*-------
*  IF IMPORTED FLUX PRESENT FOR SN, REORDER FLUX.
*-------
      IF((CXDOOR.EQ.'SN').AND.(INITFL.EQ.2)) THEN
         CALL LCMLEN(IPFLUX,'KEYFLX',ILINIT,ITYLCM)
         IF(ILINIT.NE.NREG) THEN
            WRITE(*,*) NREG, ILINIT
            CALL XABORT('FLU2DR: NUMBER OF REGIONS FROM SPN CALCULATION'
     1      //' (OBTAINED FROM LENGTH OF KEYFLX) DOES NOT MATCH NUMBER '
     2      //'OF REGIONS IN SN CALCULATION. CHECK INPUT FILE FOR '
     3      //'POTENTIAL ERRORS.')
         ENDIF
         KEYSPN(:) = 0
         CALL LCMGET(IPFLUX,'KEYFLX',KEYSPN)
         DO 25 IG=1,NGRP
            CALL SNEST(IPTRK,IPRT,NREG,NUNKNO,MATCOD,IG,KEYFLX,KEYSPN,
     1         FLUX(:,IG,2))
   25    CONTINUE
      ENDIF
*----
*  COMPUTE FIRST K-EFFECTIVE
*----
      IF((ITYPEC.EQ.0).OR.(ITYPEC.EQ.5)) THEN
         AKEFFO=1.0D0
         AKEFF=1.0D0
         AFLNOR=1.0D0
      ELSE IF(ITYPEC.EQ.1) THEN
         CALL LCMGET(IPFLUP,'STATE-VECTOR',JPAR)
         IF(JPAR(6).GE.2) THEN
            CALL LCMGET(IPFLUP,'K-EFFECTIVE',RKEFF)
            CALL LCMGET(IPFLUP,'K-INFINITY',CUREIN)
         ENDIF
         AKEFF=RKEFF
         IF(JPAR(6).GE.3) THEN
            CALL XDRSET(B2,4,0.0)
            CALL LCMGET(IPFLUP,'B2  B1HOM',B2(4))
            CALL LCMGET(IPFLUP,'DIFFB1HOM',DIFF)
         ENDIF
         IF((JPAR(6).GT.2).AND.(JPAR(7).GE.6)) THEN
            CALL LCMGET(IPFLUP,'B2  HETE',B2)
         ENDIF
         IF((JPAR(6).GT.2).AND.(JPAR(7).GE.5)) THEN
            CALL LCMGET(IPFLUP,'GAMMA',GAMMA)
         ENDIF
         AKEFFO=AKEFF
         AKEEP(2)=AKEFF
         AFLNOR=1.0D0/RKEFF
      ELSE
         OLDBIL=0.0D0
         CALL FLUKEF(IPRT,IPMACR,NGRP,NREG,NUNKNO,NMAT,NIFIS,MATCOD(1),
     1   VOL,KEYFLX(1,1,1),XSTOT,XSNUF,XSCHI,DIFF,FLUX(1,1,2),B2,ILEAK,
     2   LEAKSW,OLDBIL,AKEFF,AFLNOR)
         AKEFFO=AKEFF
         AKEEP(2)=AKEFF
      ENDIF
      B2VALO=B2(4)
*
      NCTOT=NCPTA+NCPTL
      IF(NCPTA.EQ.0) THEN
         NCPTM=NCTOT+1
      ELSE
         NCPTM=NCPTL
      ENDIF
      MESSOU='     '
      IF(IPRT.GT.0) WRITE(6,1090) 0,1.0,EPSOUT,AKEFFO,B2(4)
*----
*  CALCULATION OF THE INITIAL LEAKAGE COEFFICIENTS
*----
      IF(ITYPEC.GT.2) THEN
         CALL XDRSET(DIFF,NGRP,0.0)
         IF(OPTION.EQ.'LKRD') THEN
            CALL LCMLEN(IPMACR,'DIFFB1HOM',ILONG,ITYLCM)
            IF(ILONG.EQ.NGRP) THEN
               CALL LCMGET(IPMACR,'DIFFB1HOM',DIFF)
            ELSE
               CALL XABORT('FLU2DR: UNABE TO RECOVER THE DIFFB1HOM REC'
     >         //'ORD IN THE MACROLIB OBJECT.')
            ENDIF
            CALL XDRSET(GAMMA,NGRP,1.0)
         ELSE IF(OPTION.EQ.'RHS') THEN
            CALL LCMLEN(IPFLUX,'DIFFB1HOM',ILONG,ITYLCM)
            IF(ILONG.EQ.NGRP) THEN
              IF(LFORW) THEN
                 CALL LCMGET(IPFLUX,'DIFFB1HOM',DIFF)
              ELSE
*                Permute the diffusion coefficients if the LKRD option
*                is set for an adjoint calculation
                 CALL LCMGET(IPFLUX,'DIFFB1HOM',GAMMA)
                 DO IG=1,NGRP
                   DIFF(IG)=GAMMA(NGRP-IG+1)
                 ENDDO
               ENDIF
            ELSE
               CALL XABORT('FLU2DR: UNABE TO RECOVER THE DIFFB1HOM REC'
     >         //'ORD IN THE FLUX OBJECT.')
            ENDIF
         ELSE
            CALL B1HOM(IPMACR,LEAKSW,NUNKNO,OPTION,'DIFF',NGRP,NREG,
     1      NMAT,NIFIS,VOL,MATCOD,KEYFLX(1,1,1),FLUX(1,1,2),REFKEF,
     2      IPRT,DIFF,GAMMA,AKEFF,INORM,B2,OLDBIL)
            CALL LCMPUT(IPFLUX,'B2  B1HOM',1,2,B2(4))
         ENDIF
         CALL LCMPUT(IPFLUX,'DIFFB1HOM',NGRP,2,DIFF)
      ENDIF
*
****  OUTER LOOP  ******************************************
      IGDEB=1
      CFLI=.FALSE.
      CEXE=.FALSE.
      DO 400 IT=1,MAXOUT
      CALL KDRCPU(CPU1)
      MESSIN='     '
*----
*  FISSION SOURCE CALCULATION IN FLUX(:,:,4)
*----
      IF((ITYPEC.EQ.0).AND.(.NOT.C_ASSOCIATED(IPSOU))) THEN
        ALLOCATE(FXSOR(0:NMAT))
        JPMACR=LCMGID(IPMACR,'GROUP')
        CALL XDRSET(FLUX(1,1,4),NUNKNO*NGRP,0.0)
        DO IG=1,NGRP
          KPMACR=LCMGIL(JPMACR,IG)
          FXSOR(0)=0.0
          CALL LCMGET(KPMACR,'FIXE',FXSOR(1))
          DO IR=1,NREG
            IND=KEYFLX(IR,1,1)
            IF(IND.GT.0) FLUX(IND,IG,4)=FXSOR(MATCOD(IR))
          ENDDO
        ENDDO
        DEALLOCATE(FXSOR)
      ELSE IF(ITYPEC.EQ.0) THEN
        DO IG=1,NGRP
          IF(LFORW) THEN
            CALL LCMGDL(JPSOU,IG,FLUX(1,IG,4))
          ELSE
            CALL LCMGDL(JPSOU,NGRP-IG+1,FLUX(1,IG,4))
          ENDIF
        ENDDO
      ELSE IF(ITYPEC.EQ.1) THEN
        DO IG=1,NGRP
          IF(LFORW) THEN
            CALL LCMGDL(JPSOU,IG,FLUX(1,IG,4))
          ELSE
            CALL LCMGDL(JPSOU,NGRP-IG+1,FLUX(1,IG,4))
          ENDIF
          DO IUN=1,NUNKNO
            FLUX(IUN,IG,4)=-FLUX(IUN,IG,4)
          ENDDO
          DO IA=1,NFUNL
            DO IE=1,NLIN
              DO IR=1,NREG
                IUN=KEYFLX(IR,IE,IA)
                IF(IUN.GT.0) FLUX(IUN,IG,4)=FLUX(IUN,IG,4)/VOL(IR)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ELSE
        CALL XDRSET(FLUX(1,1,4),NUNKNO*NGRP,0.0)
      ENDIF
      IF(NIFIS.GT.0) THEN
        DO 30 IR=1,NREG
        IBM=MATCOD(IR)
        IF(IBM.GT.0) THEN
          DO IE=1,NLIN
            IND=KEYFLX(IR,IE,1)
            DO IS=1,NIFIS
              FISOUR=0.0D0
              DO IG=1,NGRP
                FISOUR=FISOUR+XSNUF(IBM,IS,IG)*FLUX(IND,IG,2)
              ENDDO
              FISOUR=FISOUR*AFLNOR
              DO  IG=1,NGRP
                FLUX(IND,IG,4)=FLUX(IND,IG,4)+XSCHI(IBM,IS,IG)*
     >          REAL(FISOUR)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
   30   CONTINUE
      ENDIF
      DO IG=1,NGRP
        XCSOU(IG)=0.0D0
        DO IR=1,NREG
          IND=KEYFLX(IR,1,1)
          IF(IND.GT.0) XCSOU(IG)=XCSOU(IG)+FLUX(IND,IG,4)*VOL(IR)
        ENDDO
      ENDDO
      ISBS=-1
      IF(C_ASSOCIATED(IPSOU)) THEN
        CALL LCMLEN(IPSOU,'NBS',ISBS,ITYLCM)
      ENDIF
*----
*  SET THE STARTING ENERGY GROUP
*----
      DO 40 IG=1,NGRP
      IGDEB=IG
      IF(XCSOU(IG).NE.0.0.OR.ISBS.NE.0) GO TO 45
   40 CONTINUE
*----
*  DOWNLOAD FROM EXTERNAL FLUX(:,:,2) TO PRESENT INTERNAL FLUX(:,:,6)
*----
   45 DO 55 IG=1,NGRP
      DO 50 IND=1,NUNKNO
      FLUX(IND,IG,6)=FLUX(IND,IG,2)
   50 CONTINUE
   55 CONTINUE
*
****  INNER LOOP  ******************************************
      DO 270 JT=1,MAXINR
      JPMACR=LCMGID(IPMACR,'GROUP')
      DO 140 IG=IGDEB,NGRP
      DO IND=1,NUNKNO
        FLUX(IND,IG,8)=FLUX(IND,IG,4)
      ENDDO
*----
*  PROCESS SELF-SCATTERING REDUCTION IN INNER SOURCES.
*----
      IF((ITPIJ.EQ.2).OR.(ITPIJ.EQ.4)) THEN
         IF((CXDOOR.EQ.'BIVAC').OR.(CXDOOR.EQ.'TRIVAC')) THEN
            CALL XABORT('FLU2DR: SCATTERING REDUCTION IS MANDATORY.')
         ENDIF
         DO 63 IR=1,NREG
         IBM=MATCOD(IR)
         DO 62 IE=1,NLIN
         DO 61 IAL=0,MIN(NLF-1,NANIS)
         XXS=XSDIA(IBM,IAL,IG)*REAL(2*IAL+1)
         DO 60 IAM=0,IAL
         IND=KEYFLX(IR,IE,1+IAL*(IAL+1)/2+IAM)
         IF(IND.GT.0) FLUX(IND,IG,8)=FLUX(IND,IG,8)+XXS*FLUX(IND,IG,6)
   60    CONTINUE
   61    CONTINUE
   62    CONTINUE
   63    CONTINUE
      ENDIF
      IF(ILEAK.EQ.5) THEN
*        ECCO ISOTROPIC STREAMING MODEL.
         CCLBD=0.0
         IF((ITPIJ.EQ.1).OR.(ITPIJ.EQ.3).AND.(OPTION.EQ.'B1')) 
     >    CCLBD=1.0-GAMMA(IG)    
         DO 75 IE=1,NLIN
         DO 70 IR=1,NREG
         IBM=MATCOD(IR)
         IND=NUNKNO/2+KEYFLX(IR,IE,1)
         IF(IND.EQ.NUNKNO/2) GO TO 70
         IF(OPTION(2:2).EQ.'1') THEN
*           B1 OR P1 CASE.
            IF(ITPIJ.EQ.2) THEN
               FLUX(IND,IG,8)=FLUX(IND,IG,8)+XSDIA(IBM,1,IG)*
     >         FLUX(IND,IG,6)
            ENDIF
         ELSE IF(ITPIJ.EQ.1) THEN
*           B0, P0, B0TR OR P0TR CASE.
            FLUX(IND,IG,8)=FLUX(IND,IG,8)-XSDIA(IBM,1,IG)*
     >      FLUX(IND,IG,6)*GAMMA(IG)
         ENDIF
         FLUX(IND,IG,8)=FLUX(IND,IG,8)+CCLBD*XSDIA(IBM,1,IG)*
     >   FLUX(IND,IG,6)
  70     CONTINUE
  75     CONTINUE
      ELSE IF(ILEAK.GE.6) THEN
*        TIBERE ANISOTROPIC STREAMING MODEL.
         CCLBD=0.0
         IF((ITPIJ.EQ.3).AND.(OPTION.EQ.'B1')) CCLBD=1.0-GAMMA(IG)
         DO 86 IE=1,NLIN
         DO 85 IR=1,NREG
         IND0=KEYFLX(IR,IE,1)
         IF(IND0.EQ.0) GO TO 85
         IBM=MATCOD(IR)
         INDD(1)=NUNKNO/4+IND0
         INDD(2)=NUNKNO/2+IND0
         INDD(3)=3*NUNKNO/4+IND0
         DO 80 IDIR=1,3
         IND=INDD(IDIR)
         IF(OPTION(2:2).EQ.'1') THEN
*           B1 OR P1 CASE.
            IF(ITPIJ.EQ.4) THEN
               FLUX(IND,IG,8)=FLUX(IND,IG,8)+XSDIA(IBM,1,IG)*
     >         FLUX(IND,IG,6)
            ENDIF
         ELSE IF(ITPIJ.EQ.3) THEN
*           B0, P0, B0TR OR P0TR CASE.
            FLUX(IND,IG,8)=FLUX(IND,IG,8)-XSDIA(IBM,1,IG)*
     >      FLUX(IND,IG,6)*GAMMA(IG)
         ENDIF
         FLUX(IND,IG,8)=FLUX(IND,IG,8)+CCLBD*XSDIA(IBM,1,IG)*
     >   FLUX(IND,IG,6)
  80     CONTINUE
  85     CONTINUE
  86     CONTINUE
      ENDIF
*----
*  COMPUTE INNER SOURCES ASSUMING SELF-SCATTERING REDUCTION.
*  LSCAL is true for TRAFIC, SYBIL, BIVAC, TRIVAC, and SN
*----
      IF(.NOT.LSCAL) THEN 
         KPMACR=LCMGIL(JPMACR,IG)
         IF((CXDOOR.EQ.'SN').AND.(IBFP.EQ.0)) THEN
            NUNK2=NUNKNO
            IF(ILEAK.EQ.5) NUNK2=NUNKNO/2
            CALL SNSOUR(NUNKNO,IG,IPTRK,KPMACR,NANIS,NREG,NMAT,NUNK2,
     1      NGRP,MATCOD,FLUX(1,1,6),FLUX(1,1,8))
         ELSE IF(CXDOOR.EQ.'SN') THEN
            NUNK2=NUNKNO
            IF(ILEAK.EQ.5) NUNK2=NUNKNO/2
            IPSTR=LCMGID(IPSYS,'STREAMING')
            JPSTR=LCMGID(IPSTR,'GROUP')
            KPSYS=LCMGIL(JPSTR,IG)
            CALL SNSBFP(IG,IPTRK,KPMACR,KPSYS,NANIS,IELEM,NLF,NREG,
     1      NMAT,NUNK2,NGRP,MATCOD,FLUX(1,1,6),FLUX(1,1,8))
         ELSE
            DO 105 IAL=0,MIN(NLF-1,NANIS)
            CALL LCMGET(KPMACR,'NJJS'//CAN(IAL),NJJ(1))
            CALL LCMGET(KPMACR,'IJJS'//CAN(IAL),IJJ(1))
            CALL LCMGET(KPMACR,'IPOS'//CAN(IAL),IPOS(1))
            CALL LCMGET(KPMACR,'SCAT'//CAN(IAL),XSCAT(1))
            DO 100 IR=1,NREG
            IBM=MATCOD(IR)
            IF(IBM.GT.0) THEN
               DO 92 IE=1,NLIN
               DO 91 IAM=0,IAL
               IND=KEYFLX(IR,IE,1+IAL*(IAL+1)/2+IAM)
               JG=IJJ(IBM)
               DO 90 JND=1,NJJ(IBM)
               IF(JG.NE.IG) THEN
                  FLUX(IND,IG,8)=FLUX(IND,IG,8)+REAL(2*IAL+1)*
     >            XSCAT(IPOS(IBM)+JND-1)*FLUX(IND,JG,6)
               ENDIF
               JG=JG-1
  90           CONTINUE
  91           CONTINUE
  92           CONTINUE
            ENDIF
 100        CONTINUE
 105        CONTINUE
         ENDIF
         IF((ILEAK.EQ.5).AND.(OPTION(2:2).EQ.'1')) THEN
*           ECCO ISOTROPIC STREAMING MODEL.
            CALL LCMGET(KPMACR,'NJJS01',NJJ(1))
            CALL LCMGET(KPMACR,'IJJS01',IJJ(1))
            CALL LCMGET(KPMACR,'IPOS01',IPOS(1))
            CALL LCMGET(KPMACR,'SCAT01',XSCAT(1))
            DO 125 IE=1,NLIN
            DO 120 IR=1,NREG
            IBM=MATCOD(IR)
            IF(IBM.GT.0) THEN
               IND=NUNKNO/2+KEYFLX(IR,IE,1)
               JG=IJJ(IBM)
               DO 110 JND=1,NJJ(IBM)
               IF(JG.NE.IG) THEN
                  FLUX(IND,IG,8)=FLUX(IND,IG,8)+XSCAT(IPOS(IBM)+JND-1)*
     >            FLUX(IND,JG,6)
               ENDIF
               JG=JG-1
  110          CONTINUE
            ENDIF
  120       CONTINUE
  125       CONTINUE
         ELSE IF(ILEAK.GE.6) THEN
*           TIBERE ANISOTROPIC STREAMING MODEL.
            CALL LCMGET(KPMACR,'NJJS01',NJJ(1))
            CALL LCMGET(KPMACR,'IJJS01',IJJ(1))
            CALL LCMGET(KPMACR,'IPOS01',IPOS(1))
            CALL LCMGET(KPMACR,'SCAT01',XSCAT(1))
            DO IE=1,NLIN
              DO IR=1,NREG
                IBM=MATCOD(IR)
                IF(IBM.GT.0) THEN
                  INDD(1)=NUNKNO/4+KEYFLX(IR,IE,1)
                  INDD(2)=NUNKNO/2+KEYFLX(IR,IE,1)
                  INDD(3)=3*NUNKNO/4+KEYFLX(IR,IE,1)
                  DO IDIR=1,3
                    IND=INDD(IDIR)
                    JG=IJJ(IBM)
                    DO JND=1,NJJ(IBM)
                      IF(JG.NE.IG) THEN
                        FLUX(IND,IG,8)=FLUX(IND,IG,8)+
     >                  XSCAT(IPOS(IBM)+JND-1)*FLUX(IND,JG,6)
                      ENDIF
                      JG=JG-1
                    ENDDO
                  ENDDO
                ENDIF
              ENDDO
            ENDDO
         ENDIF
         DO IND=1,NUNKNO
           FLUX(IND,IG,7)=FLUX(IND,IG,6)
         ENDDO
      ENDIF
  140 CONTINUE
*----
*  FLUX COMPUTATION
*----
      CALL XDISET(NPSYS,NGRP,0)
      DO 150 IG=IGDEB,NGRP
      NPSYS(IG)=IG
  150 CONTINUE
      JPSTR=C_NULL_PTR
      IF(C_ASSOCIATED(IPSYS)) THEN
         JPSYS=LCMGID(IPSYS,'GROUP')
         IF(ILEAK.EQ.5.OR.((MOD(ILEAK,10).EQ.6).AND.(IPHASE.EQ.1))) THEN
            IPSTR=LCMGID(IPSYS,'STREAMING')
            JPSTR=LCMGID(IPSTR,'GROUP')
         ENDIF
      ENDIF
      IF((.NOT.LSCAL).AND.(ILEAK.EQ.0)) THEN
         IDIR=0
         CALL DOORFV(CXDOOR,JPSYS,NPSYS,IPTRK,IFTRAK,IPRT,NGRP,
     1   NMAT,IDIR,NREG,NUNKNO,IPHASE,LEXAC,MATCOD,VOL,KEYFLX,TITLE,
     2   FLUX(1,1,8),FLUX(1,1,7),IPREB,IPSOU,REBFLG,FLUXC)
      ELSE IF(.NOT.LSCAL) THEN
         CALL FLUDBV(CXDOOR,IPHASE,JPSYS,JPSTR,NPSYS,IPTRK,IFTRAK,IPRT,
     1   NREG,NUNKNO,NFUNL,NGRP,NMAT,LEXAC,MATCOD,VOL,KEYFLX,TITLE,
     2   ILEAK,B2,DIFF,GAMMA,FLUX(1,1,6),IPREB,IPSOU,REBFLG,FLUXC)
      ELSE
*        A GROUP-BY-GROUP SCALAR PROCEDURE IS BEEN USED.
         IF(.NOT.C_ASSOCIATED(IPSYS)) THEN
            CALL XABORT('FLU2DR: MISSING L_PIJ OBJECT.')
         ENDIF
         DO 165 IG=1,NGRP
         DO 160 IND=1,NUNKNO
         FLUX(IND,IG,7)=FLUX(IND,IG,6)
  160    CONTINUE
  165    CONTINUE
         KPSTR=C_NULL_PTR
         DO 230 IG=IGDEB,NGRP
         IF(IPRT.GT.10) WRITE(6,'(/25H FLU2DR: PROCESSING GROUP,I5,
     >   6H WITH ,A,1H.)') IG,CXDOOR
         KPMACR=LCMGIL(JPMACR,IG)
         NUNK2=NUNKNO
         IF(ILEAK.EQ.5) NUNK2=NUNKNO/2
         IF(CXDOOR.EQ.'BIVAC') THEN
            CALL BIVSOU(NUNKNO,IG,IPTRK,KPMACR,NANIS,NREG,NMAT,NUNK2,
     1      NGRP,MATCOD,VOL,FLUX(1,1,7),FLUX(1,1,8))
         ELSE IF(CXDOOR.EQ.'TRIVAC') THEN
            CALL TRIVSO(NUNKNO,IG,IPTRK,KPMACR,NANIS,NREG,NMAT,NUNK2,
     1      NGRP,MATCOD,KEYFLX(1,1,1),VOL,FLUX(1,1,7),FLUX(1,1,8))
         ELSE IF((CXDOOR.EQ.'SN').AND.(IBFP.EQ.0)) THEN
            CALL SNSOUR(NUNKNO,IG,IPTRK,KPMACR,NANIS,NREG,NMAT,NUNK2,
     1      NGRP,MATCOD,FLUX(1,1,7),FLUX(1,1,8))
         ELSE IF(CXDOOR.EQ.'SN') THEN
            JPSYS=LCMGID(IPSYS,'GROUP')
            KPSYS=LCMGIL(JPSYS,IG)
            CALL SNSBFP(IG,IPTRK,KPMACR,KPSYS,NANIS,IELEM,NLF,NREG,
     1      NMAT,NUNK2,NGRP,MATCOD,FLUX(1,1,7),FLUX(1,1,8))
         ELSE
            DO 182 IAL=0,MIN(NLF-1,NANIS)
            CALL LCMGET(KPMACR,'NJJS'//CAN(IAL),NJJ(1))
            CALL LCMGET(KPMACR,'IJJS'//CAN(IAL),IJJ(1))
            CALL LCMGET(KPMACR,'IPOS'//CAN(IAL),IPOS(1))
            CALL LCMGET(KPMACR,'SCAT'//CAN(IAL),XSCAT(1))
            DO 181 IE=1,NLIN
            DO 180 IR=1,NREG
            IBM=MATCOD(IR)
            IF(IBM.GT.0) THEN
               DO 175 IAM=0,IAL
               IND=KEYFLX(IR,IE,1+IAL*(IAL+1)/2+IAM)
               JG=IJJ(IBM)
               DO 170 JND=1,NJJ(IBM)
               IF(JG.NE.IG) THEN
                  FLUX(IND,IG,8)=FLUX(IND,IG,8)+REAL(2*IAL+1)*
     >            XSCAT(IPOS(IBM)+JND-1)*FLUX(IND,JG,7)
               ENDIF
               JG=JG-1
  170          CONTINUE
  175          CONTINUE
            ENDIF
  180       CONTINUE
  181       CONTINUE
  182       CONTINUE
         ENDIF
         IF((ILEAK.GE.5).AND.(OPTION(2:2).EQ.'1')) THEN
            CALL LCMGET(KPMACR,'NJJS01',NJJ(1))
            CALL LCMGET(KPMACR,'IJJS01',IJJ(1))
            CALL LCMGET(KPMACR,'IPOS01',IPOS(1))
            CALL LCMGET(KPMACR,'SCAT01',XSCAT(1))
            IF(ILEAK.EQ.5) THEN
*              ECCO ISOTROPIC STREAMING MODEL.
               KPSTR=LCMGIL(JPSTR,IG)
               DO 205 IE=1,NLIN
               DO 200 IR=1,NREG
               IBM=MATCOD(IR)
               IF(IBM.GT.0) THEN
                  IND=NUNKNO/2+KEYFLX(IR,IE,1)
                  JG=IJJ(IBM)
                  DO 190 JND=1,NJJ(IBM)
                  IF(JG.NE.IG) THEN
                     FLUX(IND,IG,8)=FLUX(IND,IG,8)+
     >               XSCAT(IPOS(IBM)+JND-1)*FLUX(IND,JG,7)
                  ENDIF
                  JG=JG-1
  190             CONTINUE
               ENDIF
  200          CONTINUE
  205          CONTINUE
            ELSE IF(ILEAK.GE.6) THEN
*              TIBERE ANISOTROPIC STREAMING MODEL.
               DO 225 IE=1,NLIN
               DO 220 IR=1,NREG
               IBM=MATCOD(IR)
               IF(IBM.GT.0) THEN
                  INDD(1)=NUNKNO/4+KEYFLX(IR,IE,1)
                  INDD(2)=NUNKNO/2+KEYFLX(IR,IE,1)
                  INDD(3)=3*NUNKNO/4+KEYFLX(IR,IE,1)
                  DO 215 IDIR=1,3
                  IND=INDD(IDIR)
                  JG=IJJ(IBM)
                  DO 210 JND=1,NJJ(IBM)
                  IF(JG.NE.IG) THEN
                     FLUX(IND,IG,8)=FLUX(IND,IG,8)+
     >               XSCAT(IPOS(IBM)+JND-1)*FLUX(IND,JG,7)
                  ENDIF
                  JG=JG-1
  210             CONTINUE
  215             CONTINUE
               ENDIF
  220          CONTINUE
  225          CONTINUE
            ENDIF
         ENDIF
*
         NPSYS(:NGRP)=0
         NPSYS(IG)=IG
         IF(ILEAK.EQ.0) THEN
           IDIR=0
           CALL DOORFV(CXDOOR,JPSYS,NPSYS,IPTRK,IFTRAK,IPRT,NGRP,
     1     NMAT,IDIR,NREG,NUNKNO,IPHASE,LEXAC,MATCOD,VOL,KEYFLX,TITLE,
     2     FLUX(1,1,8),FLUX(1,1,7),IPREB,IPSOU,REBFLG,FLUXC)
         ELSE
           CALL FLUDBV(CXDOOR,IPHASE,JPSYS,JPSTR,NPSYS,IPTRK,IFTRAK,
     1     IPRT,NREG,NUNKNO,NFUNL,NGRP,NMAT,LEXAC,MATCOD,VOL,KEYFLX,
     2     TITLE,ILEAK,B2,DIFF,GAMMA,FLUX(1,1,6),IPREB,IPSOU,REBFLG,
     3     FLUXC)
         ENDIF
  230    CONTINUE
      ENDIF
      IF(LREBAL.AND.(ITYPEC.NE.5)) THEN
         CALL FLUBAL(IPMACR,NGRP,ILEAK,NMAT,NREG,ICREB,NUNKNO,NANIS,
     1   MATCOD,VOL,KEYFLX(1,1,1),XSTRC,XSDIA,XCSOU,IGDEB,B2,DIFF,
     2   KEYCUR,MATALB(NNN+1),ALBEDO,V(NNN+1),FLUX(1,1,7))
      ENDIF
*----
*  ACCELERATING INNER ITERATIONS CYCLICALLY DEPENDING ON PARAM.
*----
      IF(MOD(JT-1,NCTOT).GE.NCPTM) THEN
         CALL FLU2AC(NGRP,NUNKNO,IGDEB,FLUX(1,1,5),AKEEP(5),ZMU)
      ELSE
         ZMU=1.0
      ENDIF
*----
*  CALCULATING ERROR AND PREC BETWEEN PRESENT AND NEW FLUX FOR
*  EACH GROUP. RETAIN LARGEST ERROR BETWEEN ANY GROUP.
*----
      EINN=0.0
      ICHAR=ICHAR+1
      ICTOT=ICTOT+NGRP-IGDEB+1
      IGDEBO=IGDEB
      DO 260 IG=IGDEBO,NGRP
      GINN=0.0
      FINN=0.0
      DO 240 IR=1,NREG
      IND=KEYFLX(IR,1,1)
      IF(IND.EQ.0) GO TO 240
      GINN=MAX(GINN,ABS(FLUX(IND,IG,6)-FLUX(IND,IG,7)))
      FINN=MAX(FINN,ABS(FLUX(IND,IG,7)))
  240 CONTINUE
      DO 250 IND=1,NUNKNO
      FLUX(IND,IG,5)=FLUX(IND,IG,6)
      FLUX(IND,IG,6)=FLUX(IND,IG,7)
  250 CONTINUE
      GINN=GINN/FINN
      IF((GINN.LT.EPSINR).AND.(IGDEB.EQ.IG)) THEN
         IGDEB=IGDEB+1
      ELSEIF((IGDEB.EQ.IG).AND.(IG.LT.NGRP)) THEN
         ERRDEB1=GINN
      ENDIF
      EINN=MAX(EINN,GINN)
  260 CONTINUE
*
      ITERF=JT
      IF(IPRT.GT.0) WRITE(6,1080) JT,EINN,EPSINR,IGDEB,ZMU
      IF((IPRT.GT.0).AND.(IGDEB.GT.1).AND.(IGDEB.LE.NGRP)) THEN
         WRITE(6,1082) ERRDEB1 
      ENDIF
      IF(EINN.LT.EPSINR) THEN
*        thermal convergence is reached
         CFLI=CEXE
         GOTO 280
      ENDIF
*     near convergence (eps < 10.0 criterion) a new outer iteration
*     is started
      IF((IGDEB.GT.1).AND.(EINN.LT.10.*EPSINR)) GOTO 281
  270 CONTINUE
****  END OF INNER LOOP  ******************************************
*
  281 MESSIN='*NOT*'
  280 IF(LREBAL) THEN
         IF(LEAKSW) THEN
            IF(ICREB.EQ.0) THEN
               WRITE(6,*) ' *** INCOMPATIBILITY ON LEAKAGE SWITCH ***'
               CALL XABORT('FLU2DR: ERROR ON LEAKAGE SWITCH')
            ELSE
               IF(IPRT.GT.0) 
     &            WRITE(6,*) 'FLU2DR: LEAKAGE & ICREB -> REBALANCING ON'
            ENDIF
         ELSE
            IF(IPRT.GT.0) 
     &         WRITE(6,*) 'FLU2DR: NO LEAKAGE-> REBALANCING ON'
         ENDIF
      ELSE
         IF(IPRT.GT.0) WRITE(6,*) 'FLU2DR:    LEAKAGE-> REBALANCING OFF'
      ENDIF
      CALL KDRCPU(CPU2)
      IF(IPRT.GT.0) WRITE(6,1040) CPU2-CPU1,'INTERNAL',MESSIN,ITERF
*----
*  PROMOTE FROM NEW INTERNAL FLUX(,,,7) TO NEW EXTERNAL FLUX(,,,3)
*----
      DO 295 IG=1, NGRP
      DO 290 IND=1,NUNKNO
      FLUX(IND,IG,3)=FLUX(IND,IG,7)
  290 CONTINUE
  295 CONTINUE
*----
*  HOTELLING DEFLATION IN GPT CASES.
*----
      IF(ITYPEC.EQ.1) THEN
        JPFLUP1=LCMGID(IPFLUP,'FLUX')
        JPFLUP2=LCMGID(IPFLUP,'AFLUX')
        DDELN1=0.0D0
        DDELD1=0.0D0
        DO 300 IG=1,NGRP
        IF(LFORW) THEN
          CALL LCMGDL(JPFLUP1,IG,FLUX(1,IG,5)) ! EVECT
          CALL LCMGDL(JPFLUP2,IG,FLUX(1,IG,6)) ! ADECT
        ELSE
          CALL LCMGDL(JPFLUP2,NGRP-IG+1,FLUX(1,IG,5)) ! ADECT
          CALL LCMGDL(JPFLUP1,NGRP-IG+1,FLUX(1,IG,6)) ! EVECT
        ENDIF
  300   CONTINUE
        DO 335 IG=1,NGRP
        CALL XDRSET(FLUX(1,IG,7),NUNKNO,0.0)
        DO 325 IE=1,NLIN
        DO 320 IR=1,NREG
        IBM=MATCOD(IR)
        IF(IBM.GT.0) THEN
          IND=KEYFLX(IR,IE,1)
          DO 315 IS=1,NIFIS
          DO 310 JG=1,NGRP
          FLUX(IND,IG,7)=FLUX(IND,IG,7)+VOL(IR)*XSNUF(IBM,IS,IG)*
     1    XSCHI(IBM,IS,JG)*FLUX(IND,JG,6)
  310     CONTINUE
  315     CONTINUE
        ENDIF
  320   CONTINUE
  325   CONTINUE
        DO 330 IND=1,NUNKNO
        DDELN1=DDELN1+FLUX(IND,IG,7)*FLUX(IND,IG,3)
        DDELD1=DDELD1+FLUX(IND,IG,7)*FLUX(IND,IG,5)
  330   CONTINUE
  335   CONTINUE
        DO 345 IG=1,NGRP
        DO 340 IND=1,NUNKNO
        FLUX(IND,IG,3)=FLUX(IND,IG,3)-REAL(DDELN1/DDELD1)*FLUX(IND,IG,5)
  340   CONTINUE
  345   CONTINUE
      ENDIF
*
      IF(ITYPEC.EQ.2) THEN
*        NO B-N LEAKAGE CALCULATION REQUIRED
         IF(B2(4).NE.0.0) CALL XABORT('FLU2DR: NON ZERO BUCKLING.')
         CALL FLUKEF(IPRT,IPMACR,NGRP,NREG,NUNKNO,NMAT,NIFIS,MATCOD(1),
     1   VOL,KEYFLX(1,1,1),XSTOT,XSNUF,XSCHI,DIFF,FLUX(1,1,3),B2,ILEAK,
     2   LEAKSW,OLDBIL,AKEFF,AFLNOR)
      ELSE IF(ITYPEC.GT.2) THEN
*        PERFORM THE B-N LEAKAGE CALCULATION.
         CALL LCMGET(IPFLUX,'DIFFB1HOM',DIFF)
         CALL B1HOM(IPMACR,LEAKSW,NUNKNO,OPTION,HTYPE(ITYPEC),NGRP,NREG,
     1   NMAT,NIFIS,VOL,MATCOD,KEYFLX(1,1,1),FLUX(1,1,3),REFKEF,IPRT,
     2   DIFF,GAMMA,AKEFF,INORM,B2,OLDBIL)
         IF(ILEAK.GE.6) THEN
*           COMPUTE THE DIRECTIONNAL BUCKLING COMPONENTS FOR TIBERE.
            IHETL=ILEAK/10-1
            IF(IHETL.GT.0) THEN
               CALL FLUBLN(IPMACR,IPRT,NGRP,NMAT,NREG,NUNKNO,NIFIS,
     1         MATCOD,VOL,KEYFLX(1,1,1),FLUX(1,1,3),IHETL,REFKEF,B2)
            ENDIF
         ENDIF
         CALL LCMPUT(IPFLUX,'DIFFB1HOM',NGRP,2,DIFF)
         CALL LCMPUT(IPFLUX,'B2  B1HOM',1,2,B2(4))
      ENDIF
      IF(ITYPEC.GE.3) THEN
         IF(B2(4).EQ.0.0) THEN
            BFLNOR=1.0D0
         ELSE
            BFLNOR=1.0D0/ABS(B2(4))
         ENDIF
         EEXT=REAL(ABS(B2(4)-B2VALO)*BFLNOR)
         B2VALO=B2(4)
      ENDIF
      IEXTF=IT
      IF((ITYPEC.GT.1).AND.(ITYPEC.LT.5)) THEN
         AFLNOR=1.0D0/AKEFF
         EEXT=REAL(ABS(AKEFF-AKEFFO)/AKEFF)
      ELSE
         EEXT=0.0
      ENDIF
      AKEEP(3)=AKEFF
*----
*  ACCELERATING INNER ITERATIONS CYCLICALLY DEPENDING ON PARAM.
*----
      IF(MOD(IT-1,NCTOT).GE.NCPTM) THEN
         CALL FLU2AC(NGRP,NUNKNO,1,FLUX(1,1,1),AKEEP(1),ZMU)
      ELSE
         ZMU=1.0
      ENDIF
*
      EINN=0.0
      IF(IPRT.GT.0) WRITE(6,1090) IT,EEXT,EPSOUT,AKEFF,B2(4)
      IF(EEXT.LT.EPSOUT) THEN
*        COMPARE FLUX FOR OUTER ITERATIONS
         DO 370 IG=1,NGRP
         GINN=0.0
         FINN=0.0
         DO 350 IR=1,NREG
         IND=KEYFLX(IR,1,1)
         IF(IND.EQ.0) GO TO 350
         GINN=MAX(GINN,ABS(FLUX(IND,IG,2)-FLUX(IND,IG,3)))
         FINN=MAX(FINN,ABS(FLUX(IND,IG,3)))
  350    CONTINUE
         DO 360 IND=1,NUNKNO
         FLUX(IND,IG,1)=FLUX(IND,IG,2)
         FLUX(IND,IG,2)=FLUX(IND,IG,3)
  360    CONTINUE
         GINN=GINN/FINN
         EINN=MAX(EINN,GINN)
  370    CONTINUE
         IF(IPRT.GT.0) WRITE(6,1100) IT,EINN,EPSUNK,AFLNOR,ZMU
         CEXE=.TRUE.
      ELSE
         DO 385 IG=1,NGRP
         DO 380 IND=1,NUNKNO
         FLUX(IND,IG,1)=FLUX(IND,IG,2)
         FLUX(IND,IG,2)=FLUX(IND,IG,3)
  380    CONTINUE
  385    CONTINUE
         IF(IPRT.GT.0) WRITE(6,1110) IT,AFLNOR,ZMU
      ENDIF
      IF(ITYPEC.GE.2) AFLNOR=(AKEFF/AKEEP(3))*AFLNOR
      AKEEP(1)=AKEEP(2)
      AKEEP(2)=AKEEP(3)
*----
*  UPDATE KEFF
*----
      AKEFFO=AKEFF
      IF((EEXT.LT.EPSOUT).AND.(EINN.LT.EPSUNK)) GO TO 410
  400 CONTINUE
      WRITE(6,*) '*** FLU2DR: CONVERGENCE NOT REACHED ***'
      WRITE(6,*) '*** FLU2DR: CONVERGENCE NOT REACHED ***'
      WRITE(6,*) '*** FLU2DR: CONVERGENCE NOT REACHED ***'
      MESSOU='*NOT*'
*
****  CONVERGENCE REACHED  ******************************************
  410 RKEFF=REAL(AKEFF)
      IF(IPRT.GE.3) THEN
         WRITE(6,1010) (IR,IR=1,NREG)
         ALLOCATE(FL(NREG))
         DO 445 IG=1,NGRP
         WRITE(6,1070) IG
         CALL XDRSET(FL,NREG,0.0)
         DO 420 IR=1,NREG
         IND=KEYFLX(IR,1,1)
         IF(IND.GT.0) FL(IR)=FLUX(IND,IG,3)*REAL(AFLNOR)
  420    CONTINUE
         WRITE(6,1020) (FL(IR),IR=1,NREG)
         DO 440 IA=2,NFUNL
         CALL XDRSET(FL,NREG,0.0)
         DO 430 IR=1,NREG
         IND=KEYFLX(IR,1,IA)
         IF(IND.GT.0) FL(IR)=FLUX(IND,IG,3)*REAL(AFLNOR)
  430    CONTINUE
         WRITE(6,1030) IA,(FL(IR),IR=1,NREG)
  440    CONTINUE
  445    CONTINUE
         DEALLOCATE(FL)
      ENDIF
*----
*  COMPUTE K-INF
*----
      IF(ITYPEC.GE.2) THEN
         FISOUR=0.0D0
         OLDBIL=0.0D0
         DO 490 IG=1,NGRP
         DO 460 IR=1,NREG
         IND=KEYFLX(IR,1,1)
         IF(IND.EQ.0) GO TO 460
         DO 450 IS=1,NIFIS
         FISOUR=FISOUR+XSNUF(MATCOD(IR),IS,IG)*FLUX(IND,IG,3)*VOL(IR)
  450    CONTINUE
         OLDBIL=OLDBIL+XSTOT(MATCOD(IR),IG)*FLUX(IND,IG,3)*VOL(IR)
  460    CONTINUE
         KPMACR=LCMGIL(JPMACR,IG)
         CALL LCMGET(KPMACR,'NJJS00',NJJ(1))
         CALL LCMGET(KPMACR,'IJJS00',IJJ(1))
         CALL LCMGET(KPMACR,'IPOS00',IPOS(1))
         CALL LCMGET(KPMACR,'SCAT00',XSCAT(1))
         DO 480 IR=1,NREG
         IBM=MATCOD(IR)
         IF(IBM.GT.0) THEN
            IND=KEYFLX(IR,1,1)
            JG=IJJ(IBM)
            DO 470 JND=1,NJJ(IBM)
            OLDBIL=OLDBIL-XSCAT(IPOS(IBM)+JND-1)*FLUX(IND,JG,3)*VOL(IR)
            JG=JG-1
  470       CONTINUE
         ENDIF
  480    CONTINUE
  490    CONTINUE
         CUREIN=REAL(FISOUR/OLDBIL)
*
*        FLUX NORMALIZATION TO KEFF.
         IF(ITYPEC.LT.5) THEN
            DO 505 IG=1,NGRP
            DO 500 IND=1,NUNKNO
            FLUX(IND,IG,3)=FLUX(IND,IG,3)*REAL(AKEFF/FISOUR)
  500       CONTINUE
  505       CONTINUE
         ENDIF
      ENDIF
      CALL KDRCPU(CPU1)
      IF(IPRT.GE.1) THEN
         WRITE(6,1040) CPU1-CPU0,'EXTERNAL',MESSOU,IEXTF
         IF(ITYPEC.EQ.0) THEN
            WRITE(6,1050) ICHAR,EEXT
         ELSE
            WRITE(6,1060) ICHAR,CUREIN,AKEFF,B2(4),EEXT
         ENDIF
         WRITE(6,1120) ICTOT
      ENDIF
*----
*  RELEASE ARRAYS
*----
      IF(CXDOOR.EQ.'MCCG') THEN
         IF(ICREB.GT.0) DEALLOCATE(MATALB,V,KEYCUR)
      ENDIF
*----
*  SAVE THE SOLUTION
*----
      DO 510 IG=1,NGRP
      IF(LFORW) THEN
         CALL LCMPDL(JPFLUX,IG,NUNKNO,2,FLUX(1,IG,3))
      ELSE
         CALL LCMPDL(JPFLUX,NGRP-IG+1,NUNKNO,2,FLUX(1,IG,3))
      ENDIF
  510 CONTINUE
      IF(C_ASSOCIATED(IPSOU)) THEN
        CALL LCMLEN(IPSOU,'NORM-FS',ILEN,ITYLCM)
        IF(ILEN.GT.0) THEN
          CALL LCMGET(IPSOU,'NORM-FS',ZNORM)
          CALL LCMPUT(IPFLUX,'NORM-FS',1,2,ZNORM)
          CALL LCMPUT(IPFLUX,'MATCOD',NREG,1,MATCOD)
        ENDIF
      ENDIF
      IF(IBFP.NE.0) THEN
        CALL LCMGET(IPSYS,'ECUTOFF',ECUTOFF)
        CALL LCMPUT(IPFLUX,'ECUTOFF',1,2,ECUTOFF)
        CALL LCMPUT(IPFLUX,'FLUXC',NREG,2,FLUXC)
      ENDIF
      IF(ITYPEC.GE.2) THEN
         CALL LCMPUT(IPFLUX,'K-EFFECTIVE',1,2,RKEFF)
         CALL LCMPUT(IPFLUX,'K-INFINITY',1,2,CUREIN)
      ENDIF
      IF(ITYPEC.GE.3) THEN
         CALL LCMPUT(IPFLUX,'B2  B1HOM',1,2,B2(4))
         CALL LCMPUT(IPFLUX,'DIFFB1HOM',NGRP,2,DIFF)
      ENDIF
      IF((ITYPEC.GT.2).AND.(ILEAK.GE.6)) THEN
         CALL LCMPUT(IPFLUX,'B2  HETE',3,2,B2)
      ENDIF
      IF((ITYPEC.GT.2).AND.(ILEAK.GE.5)) THEN
         CALL LCMPUT(IPFLUX,'GAMMA',NGRP,2,GAMMA)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(GAMMA,XSCAT,FLUX,DIFF)
      DEALLOCATE(NPSYS,IPOS,NJJ,IJJ)
      RETURN
*
 1010 FORMAT (/28H FLUXES AVERAGED           :/
     1 (9X,2H/=,I4,:,6X,2H/=,I4,:,6X,2H/=,I4,:,6X,2H/=,I4,:,6X,2H/=,
     2 I4,:,6X,2H/=,I4,:,6X,2H/=,I4,:,6X,2H/=,I4,:,6X,2H/=,I4,:,6X,
     3 2H/=,I4))
 1020 FORMAT (7H FLUX  ,2H: ,1P,10E12.5/(9X,10E12.5))
 1030 FORMAT (5H CUR ,I2,2H: ,1P,10E12.5/(9X,10E12.5))
 1040 FORMAT (18H FLU2DR: CPU TIME=, F10.0,1X,A8,
     1        13H CONVERGENCE ,A5,14H REACHED AFTER ,I6,
     2        12H ITERATIONS.  )
 1050 FORMAT (/20H ++ TRACKING CALLED=,I4,6H TIMES ,
     1         11H PRECISION=,E9.2)
 1060 FORMAT (/20H ++ TRACKING CALLED=,I4,6H TIMES ,
     1         12H FINAL KINF=,1P,E13.6,
     2         12H FINAL KEFF=,E13.6,4H B2=,E12.5,
     3         11H PRECISION=,E9.2)
 1070 FORMAT (/14H ENERGY GROUP ,I6)
 1080 FORMAT (10X,3HIN(,I3,6H) FLX:,5H PRC=,1P,E9.2,5H TAR=,E9.2,
     1 7H IGDEB=, I13,6H ACCE=,0P,F12.5)
 1082 FORMAT (18X,28HFIRST UNCONVERGED GROUP PRC=,E9.2)
 1090 FORMAT (5H OUT(,I3,6H) EIG:,5H PRC=,1P,E9.2,5H TAR=,E9.2,
     1 6H KEFF=,E13.6,6H BUCK=,E12.5)
 1100 FORMAT (5H OUT(,I3,6H) FLX:,5H PRC=,1P,E9.2,5H TAR=,E9.2,
     1 6H FNOR=,E13.6,6H ACCE=,0P,F12.5)
 1110 FORMAT (5H OUT(,I3,6H) FLX:,28X,6H FNOR=,1P,E13.6,6H ACCE=,
     1 0P,F12.5)
 1120 FORMAT (38H ++ TOTAL NUMBER OF FLUX CALCULATIONS=,I10)
      END
