*DECK FLUDRV
      SUBROUTINE FLUDRV(IPRT,IPFLUX,IPTRK,IPMACR,IPSOU,IFTRAK,IPSYS,
     1 IPHASE,ITPIJ,CXDOOR,ITRANC,TITLE,B2,INITFL,LFORW,LEAKSW,IREBAL,
     2 NGRP,NMAT,NIFIS,NANIS,NLF,NLIN,NFUNL,OPTION,NUN,MAXINR,EPSINR,
     3 MAXOUT,EPSUNK,EPSOUT,IFRITR,IACITR,ITYPEC,ILEAK,NREG,NSOUT,
     4 MATCOD,KEYFLX,VOL,REFKEF)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Driver for Boltzmann equation solvers.
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
* IPFLUX  pointer to the flux LCM object.
* IPTRK   pointer to the tracking LCM object.
* IPMACR  pointer to the macrolib LCM object.
* IPSOU   pointer to the fixed source LCM object.
* IFTRAK  tracking file unit number.
* IPSYS   pointer to the system LCM object (=0 for the method of
*         characteristics).
* IPHASE  1 for asm 2 for pij.
* ITPIJ   type of collision probability information available:
*         =1 scattering modified pij (wij);
*         =2 standard pij;
*         =3 scattering modified pij+pijk (wij,wijk);
*         =4 standard pij+pijk.
* CXDOOR  name of the flux solution door.
* ITRANC  type of transport correction (>0 to perform a correction).
* TITLE   title.
* B2      initial or imposed directional bucklings.
* INITFL  flux initialization flag (=0/1/2: uniform flux/LCM/DSA).
* LFORW   flag set to .false. to solve an adjoint problem.
* LEAKSW  leakage flag (=.true. if leakage is present on the outer
*         surface).
* IREBAL  flux rebalancing flag (=1: perform rebalancing).
* NGRP    number of energy groups.
* NMAT    number of mixtures.
* NIFIS   number of fissile isotopes.
* NANIS   maximum cross section Legendre order.
* NLF     number of Legendre orders for the flux.
* NLIN    number of polynomial components in flux spatial expansion.
* NFUNL   number of spherical harmonics components.
* OPTION  type of leakage coefficients:
*         'LKRD' (recover leakage coefficients in Macrolib);
*         'RHS' (recover leakage coefficients in RHS flux object);
*         'B0' (B-0), 'P0' (P-0), 'B1' (B-1),
*         'P1' (P-1), 'B0TR' (B-0 with transport correction) or 'P0TR'
*         (P-0 with transport correction).
* NUN     number of unknowns per energy group including spherical
*         harmonic terms, interface currents and fundamental
*         currents.
* MAXINR  maximum number of thermal iterations.
* EPSINR  thermal iterations epsilon.
* MAXOUT  maximum number of outer iterations.
* EPSUNK  outer iterations eigenvector epsilon.
* EPSOUT  outer iterations eigenvalue epsilon.
* IFRITR  number of free iterations in an acceleration cycle.
* IACITR  number of accelerated iterations in an acceleration cycle.
* ITYPEC  type of flux evaluation:
*         =-1 skip the flux calculation;
*         =0  fixed sources;
*         =1  fixed source eigenvalue problem (GPT type);
*         =2  fission sources/K-eff convergence;
*         =3  fission sources/K-eff convergence/db2 buckling evaluation;
*         =4  fission sources/db2 buckling convergence;
*         =5  b2 sources/db2 buckling convergence.
* ILEAK   method used to include DB2 effect:
*         =1 the scattering modified cp matrix is multiplied by PNLR;
*         =2 the reduced cp matrix is multiplied by PNL;
*         =3 sigs0-db2 approximation;
*         =4 albedo approximation;
*         =5 Ecco-type isotropic streaming model;
*         >5 Tibere anisotropic streaming model.
* NREG    number of regions.
* NSOUT   number of outer surfaces.
* MATCOD  mixture indices.
* KEYFLX  index of L-th order flux components in unknown vector.
* VOL     volumes.
* REFKEF  target effective multiplication factor (K-eff).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER   CXDOOR*12,TITLE*72,OPTION*4
      TYPE(C_PTR) IPFLUX,IPTRK,IPMACR,IPSOU,IPSYS
      INTEGER     IPRT,IFTRAK,IPHASE,ITPIJ,ITRANC,INITFL,IREBAL,NGRP,
     >            NMAT,NIFIS,NANIS,NLF,NLIN,NFUNL,NUN,MAXINR,MAXOUT,
     >            IFRITR,IACITR,ITYPEC,ILEAK,NREG,NSOUT,MATCOD(NREG),
     >            KEYFLX(NREG,NLIN,NFUNL)
      REAL        EPSUNK,EPSINR,B2(4),VOL(NREG)
      LOGICAL     LFORW,LEAKSW
      DOUBLE PRECISION REFKEF
*----
*  LOCAL VARIABLES
*----
      PARAMETER  (IUNOUT=6,NSTATE=40,NDIMO=2,IGPT=0)
      TYPE(C_PTR) JPMACR,KPMACR,JPFLUX,IPFLUP
      LOGICAL     LREBAL,LTRANC,LTP1
      INTEGER     ISTATE(NSTATE)
      REAL        EPSCON(5)
      CHARACTER   CAN(0:19)*2
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: VSCAT
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FLUXO,XSTOT,XSTRC,XSTK
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: XSDIA,XSCHI,XSNUF
*----
*  DATA STATEMENTS
*----
      SAVE CAN
      DATA CAN /'00','01','02','03','04','05','06','07','08','09',
     >          '10','11','12','13','14','15','16','17','18','19'/
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(FLUXO(NUN,NGRP),XSTOT(0:NMAT,NGRP),XSTRC(0:NMAT,NGRP),
     > XSDIA(0:NMAT,0:NANIS,NGRP),XSCHI(0:NMAT,NIFIS,NGRP),
     > XSNUF(0:NMAT,NIFIS,NGRP),XSTK(NMAT,NIFIS),VSCAT(NMAT))
*
      IF(IPRT.GE.3) THEN
        WRITE(IUNOUT,6000)
        WRITE(IUNOUT,6001) (IREGIO,VOL(IREGIO),MATCOD(IREGIO),
     >                      IREGIO=1,NREG)
      ENDIF
      IPFLUP=C_NULL_PTR
*----
*  RECOVER CROSS SECTIONS.
*----
      JPMACR=LCMGID(IPMACR,'GROUP')
      DO 110 IGR=1,NGRP
        KPMACR=LCMGIL(JPMACR,IGR)
        CALL LCMLEN(KPMACR,'NTOT0',ILNLCM,ITYLCM)
        IF(ILNLCM.GT.0) THEN
          CALL LCMGET(KPMACR,'NTOT0',XSTRC(1,IGR))
          XSTRC(0,IGR)=0.0
          DO 10 IR= 0, NMAT
          XSTOT(IR,IGR)=XSTRC(IR,IGR)
   10     CONTINUE
        ELSE
          CALL LCMLIB(KPMACR)
          CALL XABORT('FLUDRV: READ ERROR ON LCM RECORD= TOTAL')
        ENDIF
        DO 20 IANI=0,NANIS
        CALL LCMLEN(KPMACR,'SIGW'//CAN(IANI),ILNLCM,ITYLCM)
        IF( ILNLCM.GT.0 )THEN
          CALL LCMGET(KPMACR,'SIGW'//CAN(IANI),XSDIA(1,IANI,IGR))
          XSDIA(0,IANI,IGR)=0.0
        ELSE
          CALL XDRSET(XSDIA(0,IANI,IGR),NMAT+1,0.0)
        ENDIF
   20   CONTINUE
        CALL LCMLEN(KPMACR,'TRANC',ILCMLN,ITYLCM)
        LTRANC=ILCMLN.GT.0
        CALL LCMLEN(KPMACR,'NTOT1',ILCMLN,ITYLCM)
        LTP1=ILCMLN.GT.0
        IF( (ITRANC.NE.0).AND.LTRANC )THEN
*         TRANSPORT CORRECTION (INCLUDE THE LEAKAGE CORRECTION).
          CALL LCMGET(KPMACR,'TRANC',VSCAT)
        ELSE IF( LTP1 )THEN
*         PERFORM A LEAKAGE CORRECTION.
          CALL LCMGET(KPMACR,'NTOT1',VSCAT)
          DO 30 IMAT=1,NMAT
          VSCAT(IMAT)=XSTRC(IMAT,IGR)-VSCAT(IMAT)
   30     CONTINUE
        ELSE
          CALL XDRSET(VSCAT,NMAT,0.0)
        ENDIF
        DO 40 IMAT=1,NMAT
          XSTRC(IMAT,IGR)=XSTRC(IMAT,IGR)-VSCAT(IMAT)
          XSDIA(IMAT,0,IGR)=XSDIA(IMAT,0,IGR)-VSCAT(IMAT)
*         Tibere is using transport-corrected XS for the second
*         equation. Scattering reduction must be performed with
*         transport-corrected SIGS1 values.
          IF(ILEAK.GE.6) XSDIA(IMAT,1,IGR)=XSDIA(IMAT,1,IGR)-VSCAT(IMAT)
   40   CONTINUE
        IF(IPRT.GE.3) THEN
          WRITE(IUNOUT,6002)  IGR
          WRITE(IUNOUT,6003) (IMAT,XSTRC(IMAT,IGR),
     >                             XSDIA(IMAT,0,IGR),IMAT=1,NMAT)
        ENDIF
        CALL LCMLEN(KPMACR,'CHI',ILONG,ITYLCM)
        IF( ILONG.EQ.0 )THEN
           IF(( ITYPEC.NE.0 ).AND.( ITYPEC.NE.5 ))THEN
              CALL XABORT('FLUDRV: TYPE S OR L REQUESTED')
           ENDIF
           CALL XDRSET(XSCHI,NGRP*(NMAT+1)*NIFIS,0.0)
           CALL XDRSET(XSNUF,NGRP*(NMAT+1)*NIFIS,0.0)
        ELSE
           CALL LCMGET(KPMACR,'CHI',XSTK)
           DO 60 IFIS= 1, NIFIS
              XSCHI(0,IFIS,IGR)= 0.0
              DO 50 IMAT= 1, NMAT
                 XSCHI(IMAT,IFIS,IGR)= XSTK(IMAT,IFIS)
   50         CONTINUE
   60      CONTINUE
           CALL LCMGET(KPMACR,'NUSIGF',XSTK)
           DO 80 IFIS= 1, NIFIS
              XSNUF(0,IFIS,IGR)= 0.0
              DO 70 IMAT= 1, NMAT
                 XSNUF(IMAT,IFIS,IGR)= XSTK(IMAT,IFIS)
   70         CONTINUE
   80      CONTINUE
           IF( IPRT.GT.3 )THEN
              WRITE(IUNOUT,6004) (IMAT,XSNUF(IMAT,1,IGR),
     >                                 XSCHI(IMAT,1,IGR),IMAT=1,NMAT)
           ENDIF
        ENDIF
        DO 100 IANI=0,NANIS
          CALL LCMLEN(KPMACR,'NJJS'//CAN(IANI),ILNLCM,ITYLCM)
          IF(ILNLCM.NE.NMAT) THEN
            CALL LCMLIB(KPMACR)
            CALL XABORT('FLUDRV: FLUX CALCULATION ERROR, SCATTERING'//
     >      ' MATRIX OF ORDER ANIS ='//CAN(IANI)//' NOT ON LCM')
          ENDIF
  100   CONTINUE
  110 CONTINUE
*----
*  FLUX INITIALIZATION
*----
      IF(LFORW) THEN
        CALL LCMLEN(IPFLUX,'FLUX',ILINIT,ITYLCM)
      ELSE
        CALL LCMLEN(IPFLUX,'AFLUX',ILINIT,ITYLCM)
      ENDIF
      IF((ILINIT.EQ.0).OR.(INITFL.EQ.0)) THEN
        IF(LFORW) THEN
          JPFLUX=LCMLID(IPFLUX,'FLUX',NGRP)
        ELSE
          JPFLUX=LCMLID(IPFLUX,'AFLUX',NGRP)
        ENDIF
        DO 130 IGR=1,NGRP
          CALL XDRSET(FLUXO(1,IGR),NUN,0.0)
          IF(ITYPEC.GT.0) THEN
            DO 120 IREGIO=1,NREG
            IND=KEYFLX(IREGIO,1,1)
            IF(IND.GT.0) FLUXO(IND,IGR)=1.0
  120       CONTINUE
          ENDIF
          IF(LFORW) THEN
            CALL LCMPDL(JPFLUX,IGR,NUN,2,FLUXO(1,IGR))
          ELSE
            CALL LCMPDL(JPFLUX,NGRP-IGR+1,NUN,2,FLUXO(1,IGR))
          ENDIF
  130   CONTINUE
      ENDIF
      IF(ITYPEC.GE.2) THEN
        EIGENK=REAL(REFKEF)
        CALL LCMPUT(IPFLUX,'K-EFFECTIVE',1,2,EIGENK)
      ENDIF
      IF(ITYPEC.GE.3) THEN
        CALL LCMPUT(IPFLUX,'B2  B1HOM',1,2,B2(4))
        IF(ILEAK.GE.6) CALL LCMPUT(IPFLUX,'B2  HETE',3,2,B2)
      ENDIF
      IF(ITYPEC.EQ.-1) GO TO 1001
*
      IF(ILEAK.GE.6.AND.ITYPEC.GE.3) THEN
        IF(ITRANC.EQ.0) THEN
          IF(OPTION.EQ.'B0TR'.OR.OPTION.EQ.'P0TR'.OR.OPTION.EQ.'LKRD'
     >     .OR.OPTION.EQ.'RHS')
     >      CALL XABORT('FLUDRV: ILLEGAL OPTION = '//OPTION//
     >         ' FOR HETEROGENEOUS LEAKAGE CALCULATION'//
     >         ' WITHOUT TRANSPORT CORRECTED CROSS SECTIONS')
        ENDIF
      ENDIF
*
      IF (CXDOOR.EQ.'MCCG') THEN
         CALL LCMLEN(IPTRK,'KEYCUR$MCCG',ICREB,ITYLCM)
         CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
         IF ((ICREB.GT.0).AND.(ISTATE(24).EQ.0)) THEN
            LREBAL=(IREBAL.EQ.1)
         ELSE
            LREBAL=(IREBAL.EQ.1).AND.(.NOT.LEAKSW)
         ENDIF
      ELSE
         LREBAL=(IREBAL.EQ.1).AND.(.NOT.LEAKSW)
      ENDIF
      CALL FLU2DR(IPRT,IPMACR,IPFLUX,IPSYS,IPTRK,IPFLUP,IPSOU,IGPT,
     1 IFTRAK,CXDOOR,TITLE,NUN,NREG,NSOUT,NANIS,NLF,NLIN,NFUNL,NGRP,
     2 NMAT,NIFIS,LFORW,LEAKSW,MAXINR,EPSINR,MAXOUT,EPSUNK,EPSOUT,
     3 IFRITR,IACITR,ITYPEC,IPHASE,ITPIJ,ILEAK,OPTION,REFKEF,MATCOD,
     4 KEYFLX,VOL,XSTOT,XSTRC,XSDIA,XSNUF,XSCHI,LREBAL,INITFL)
*
 1001 CALL LCMLEN(IPFLUX,'FLUX',ILON1,ITYLCM)
      CALL LCMLEN(IPFLUX,'AFLUX',ILON2,ITYLCM)
      CALL XDISET(ISTATE,NSTATE,0)
      ISTATE(1)=NGRP
      ISTATE(2)=NUN
      IF((ILON1.GT.0).AND.(ILON2.GT.0)) THEN
        ISTATE(3)=11
      ELSE IF(ILON1.GT.0) THEN
        ISTATE(3)=1
      ELSE IF(ILON2.GT.0) THEN
        ISTATE(3)=10
      ENDIF
      ISTATE(4)=0
      ISTATE(5)=0
      ISTATE(6)=ITYPEC
      ISTATE(7)=ILEAK
      ISTATE(8)=IFRITR
      ISTATE(9)=IACITR
      ISTATE(10)=IREBAL
      ISTATE(11)=MAXINR
      ISTATE(12)=MAXOUT
      CALL LCMPUT(IPFLUX,'STATE-VECTOR',NSTATE,1,ISTATE)
      EPSCON(1)=EPSINR
      EPSCON(2)=EPSUNK
      EPSCON(3)=EPSOUT
      EPSCON(4:5)=0.0
      CALL LCMPUT(IPFLUX,'EPS-CONVERGE',5,2,EPSCON)
      CALL LCMPUT(IPFLUX,'KEYFLX',NREG,1,KEYFLX)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(VSCAT,XSTK,XSNUF,XSCHI,XSDIA,XSTRC,XSTOT,FLUXO)
      RETURN
*----
*  FORMATS
*----
 6000 FORMAT(//30X,' EDITION REGION/VOL/MIXTURE '//
     >3(5X,'REGION',5X,'VOL  ',5X,'MIXTURE')/)
 6001 FORMAT(1P,3(1X,I8,4X,E12.5,I8,4X))
 6002 FORMAT(//30X,' G R O U P : ',I5//
     >30X,' TOTAL MACROSCOPIC CROSS SECTIONS PER MIXTURE '/)
 6003 FORMAT(3(1X,'MIXTURE',4X,'NTOT0',11X,'SIGW',3X)/
     >1P,3(1X,I4,3X,E12.5,3X,E12.5))
 6004 FORMAT(3(1X,'MIXTURE',4X,'NUSIGF',11X,'CHI ',3X)/
     >1P,3(1X,I4,3X,E12.5,3X,E12.5))
      END
