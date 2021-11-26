*DECK FLUDBV
      SUBROUTINE FLUDBV(CDOOR,IPHASE,JPSYS,JPSTR,NPSYS,IPTRK,IFTRAK,
     1 IPRT,NREG,NUNKNO,NFUNL,NGRP,NMAT,LEXAC,MATCOD,VOL,KEYFLX,TITLE,
     2 ILEAK,B2,DDD,GAMMA,FLUX,IPMACR,REBFLG)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Find a leakage parameter to match the input DB2 value and find the
* corresponding flux. Vectorial version.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy and A. Hebert
*
*Parameters: input
* CDOOR   name of the geometry/solution operator.
* IPHASE  type of flux solution door (1 for asm 2 for pij).
* JPSYS   pointer to the system LCM list object.
* JPSTR   pointer to the system LCM list object containing isotropic
*         streaming information (=0 if not required).
* NPSYS   non-converged energy group indices.
* IPTRK   pointer to the tracking LCM object.
* IFTRAK  tracking file unit number.
* IPRT    print flag.
* NREG    number of regions.
* NFUNL   second dimension of matrix KEYFLX.
* NGRP    number of energy groups.
* NUNKNO  number of flux/sources unknowns per energy group.
* NMAT    number of mixtures in the internal library.
* LEXAC   type of exponential function calculation (=.false. to compute
*         exponential functions using tables).
* MATCOD  mixture indices.
* VOL     volumes.
* KEYFLX  index of L-th order flux components in unknown vector.
* TITLE   title.
* ILEAK   method used to include db2 effect:
*         =1 the scattering modified cp matrix is multiplied by PNLR;
*         =2 the reduced cp matrix is multiplied by PNL;
*         =3 sigs0-db2 approximation;
*         =4 (not available);
*         =5 Ecco-type isotropic streaming model;
*         >5 Tibere anisotropic streaming model.
* B2      buckling.
* DDD     leakage coefficients.
* GAMMA   gamma factors.
* IPMACR  pointer to the macrolib LCM object.
* REBFLG  ACA or SCR rebalancing flag.
*
*Parameters: input/output
* FLUX    neutron flux:
*         FLUX(:,:,1)  present inner;
*         FLUX(:,:,2)  new     inner;
*         FLUX(:,:,3)  source  inner.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER CDOOR*12,TITLE*72
      LOGICAL LEXAC,REBFLG
      TYPE(C_PTR) JPSYS,JPSTR,IPTRK,IPMACR
      INTEGER IPHASE,NPSYS(NGRP),IFTRAK,IPRT,NREG,NUNKNO,NFUNL,NGRP,
     1 NMAT,MATCOD(NREG),KEYFLX(NREG,NFUNL),ILEAK
      REAL VOL(NREG),FLUX(NUNKNO,NGRP,3),B2(4),DDD(NGRP),GAMMA(NGRP)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IOUT=6)
      TYPE(C_PTR) KPSYS,KPSTR
      CHARACTER   TEXT12*12
      INTEGER INDD(3),INDC(3),INDB(3)
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: SIGT0,SIGS0,SUNKNO,FUNKNO,F1,
     1 F2,PP
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SUNKN
*
      ALLOCATE(SUNKN(NUNKNO,NGRP))
      SUNKN(:NUNKNO,:NGRP)=FLUX(:NUNKNO,:NGRP,3)
*
      IF(ILEAK.EQ.1) THEN
         DO 30 IGR=1,NGRP
         IF((NPSYS(IGR).NE.0).AND.(B2(4).NE.0.0)) THEN
            KPSYS=LCMGIL(JPSYS,IGR)
            CALL LCMLEN(KPSYS,'DRAGON-TXSC',ILCTXS,ITYLCM)
            ALLOCATE(SIGT0(0:ILCTXS-1))
            CALL LCMGET(KPSYS,'DRAGON-TXSC',SIGT0(0))
            CALL LCMLEN(KPSYS,'DRAGON-S0XSC',ILCS0X,ITYLCM)
            ALLOCATE(SIGS0(0:ILCS0X-1))
            CALL LCMGET(KPSYS,'DRAGON-S0XSC',SIGS0(0))
            ZNUM=0.0
            ZDEN=0.0
            DO 10 IR=1,NREG
            IBM=MATCOD(IR)
            IND=KEYFLX(IR,1)
            ZNUM=ZNUM+(SIGT0(IBM)-SIGS0(IBM))*FLUX(IND,IGR,1)*VOL(IR)
            ZDEN=ZDEN+FLUX(IND,IGR,1)*VOL(IR)
   10       CONTINUE
            DEALLOCATE(SIGS0)
            DEALLOCATE(SIGT0)
            ALP1=ZNUM/(ZNUM+DDD(IGR)*B2(4)*ZDEN)
            DO 20 IR=1,NREG
            IND=KEYFLX(IR,1)
            FLUX(IND,IGR,3)=ALP1*FLUX(IND,IGR,3)
   20       CONTINUE
         ENDIF
   30    CONTINUE
         IDIR=0
         CALL DOORFVR(CDOOR,JPSYS,NPSYS,IPTRK,IFTRAK,IPRT,NGRP,
     1        NMAT,IDIR,NREG,NUNKNO,IPHASE,LEXAC,MATCOD,VOL,KEYFLX,
     2        TITLE,FLUX(1,1,3),FLUX(1,1,2),IPMACR,REBFLG)
      ELSE IF(ILEAK.EQ.2) THEN
         DO 50 IGR=1,NGRP
         IF((NPSYS(IGR).NE.0).AND.(B2(4).NE.0.0)) THEN
            KPSYS=LCMGIL(JPSYS,IGR)
            CALL LCMLEN(KPSYS,'DRAGON-TXSC',ILCTXS,ITYLCM)
            ALLOCATE(SIGT0(0:ILCTXS-1))
            CALL LCMGET(KPSYS,'DRAGON-TXSC',SIGT0(0))
            CALL LCMLEN(KPSYS,'DRAGON-S0XSC',ILCS0X,ITYLCM)
            ALLOCATE(SIGS0(0:ILCS0X-1))
            CALL LCMGET(KPSYS,'DRAGON-S0XSC',SIGS0(0))
            ZNUM=0.0
            ZDEN=0.0
            DO 40 IR=1,NREG
            IBM=MATCOD(IR)
            IND=KEYFLX(IR,1)
            ZNUM=ZNUM+SIGT0(IBM)*FLUX(IND,IGR,1)*VOL(IR)
            ZDEN=ZDEN+FLUX(IND,IGR,1)*VOL(IR)
   40       CONTINUE
            ALP1=ZNUM/(ZNUM+DDD(IGR)*B2(4)*ZDEN)
            DO 45 IR=1,NREG
            IND=KEYFLX(IR,1)
            FLUX(IND,IGR,3)=ALP1*FLUX(IND,IGR,3)-(1.0-ALP1)
     >                 *SIGS0(MATCOD(IR))*FLUX(IND,IGR,1)
   45       CONTINUE
            DEALLOCATE(SIGS0,SIGT0)
         ENDIF
   50    CONTINUE
         IDIR=0
         CALL DOORFVR(CDOOR,JPSYS,NPSYS,IPTRK,IFTRAK,IPRT,NGRP,
     1        NMAT,IDIR,NREG,NUNKNO,IPHASE,LEXAC,MATCOD,VOL,KEYFLX,
     2        TITLE,FLUX(1,1,3),FLUX(1,1,2),IPMACR,REBFLG)
      ELSE IF(ILEAK.EQ.3) THEN
         DO 70 IGR=1,NGRP
         IF(NPSYS(IGR).NE.0) THEN
            BB=B2(4)
            DO 60 IR=1,NREG
            IND=KEYFLX(IR,1)
            FLUX(IND,IGR,3)=FLUX(IND,IGR,3)-DDD(IGR)*BB*FLUX(IND,IGR,1)
   60       CONTINUE
         ENDIF
   70    CONTINUE
         IDIR=0
         CALL DOORFVR(CDOOR,JPSYS,NPSYS,IPTRK,IFTRAK,IPRT,NGRP,
     1        NMAT,IDIR,NREG,NUNKNO,IPHASE,LEXAC,MATCOD,VOL,KEYFLX,
     2        TITLE,FLUX(1,1,3),FLUX(1,1,2),IPMACR,REBFLG)
      ELSE IF(ILEAK.EQ.4) THEN
         ALLOCATE(F1(NREG),F2(NREG))
         DO 80 IGR=1,NGRP
         IF(NPSYS(IGR).NE.0) THEN
           KPSYS=LCMGIL(JPSYS,IGR)
           CALL LCMLEN(KPSYS,'DRAGON-TXSC',ILCTXS,ITYLCM)
           ALLOCATE(SIGT0(0:ILCTXS-1))
           CALL LCMGET(KPSYS,'DRAGON-TXSC',SIGT0(0))
           CALL LCMLEN(KPSYS,'DRAGON-S0XSC',ILCS0X,ITYLCM)
           ALLOCATE(SIGS0(0:ILCS0X-1))
           CALL LCMGET(KPSYS,'DRAGON-S0XSC',SIGS0(0))
           CALL FLUALB(KPSYS,NREG,NUNKNO,ILCTXS,MATCOD,VOL,KEYFLX,
     >     FLUX(1,IGR,1),FLUX(1,IGR,3),SIGS0(0),SIGT0(0),F1,F2)
           DEALLOCATE(SIGS0,SIGT0)
*
           IF(IPRT.GT.2) THEN
             WRITE(IOUT,'(//33H N E U T R O N    S O U R C E S :)')
             WRITE(IOUT,'(1P,6(5X,E15.7))') (FLUX(KEYFLX(I,1),IGR,3),
     >       I=1,NREG)
           ENDIF
           CALL XDRSET(FLUX(1,IGR,2),NUNKNO,0.0)
           DO 75 I=1,NREG
           FLUX(KEYFLX(I,1),IGR,2)=F1(I)+DDD(IGR)*B2(4)*F2(I)
   75      CONTINUE
           IF(IPRT.GT.2) THEN
             WRITE(IOUT,'(//33H N E U T R O N    F L U X E S   :)')
             WRITE(IOUT,'(1P,6(5X,E15.7))') (FLUX(KEYFLX(I,1),IGR,2),
     >       I=1,NREG)
           ENDIF
         ENDIF
   80    CONTINUE
         DEALLOCATE(F2,F1)
      ELSE IF(ILEAK.EQ.5) THEN
*        ISOTROPIC STREAMING MODEL (ECCO).
         IF(.NOT.C_ASSOCIATED(JPSTR)) THEN
            CALL XABORT('FLUDBV: MISSING STREAMING INFO(1).')
         ENDIF
         DO 95 IGR=1,NGRP
         IF(NPSYS(IGR).NE.0) THEN
            BB=B2(4)
            DO 90 IR=1,NREG
            IND=KEYFLX(IR,1)
            FLUX(IND,IGR,3)=FLUX(IND,IGR,3)-BB*FLUX(NUNKNO/2+IND,IGR,1)
   90       CONTINUE
         ENDIF
   95    CONTINUE
         IF(IPRT.GE.3) WRITE(IOUT,'(28H FLUDBV: FUNDAMENTAL FLUXES.)')
         IDIR=0
         CALL DOORFVR(CDOOR,JPSYS,NPSYS,IPTRK,IFTRAK,IPRT,NGRP,
     1        NMAT,IDIR,NREG,NUNKNO,IPHASE,LEXAC,MATCOD,VOL,KEYFLX,
     2        TITLE,FLUX(1,1,3),FLUX(1,1,2),IPMACR,REBFLG)
         DO 130 IGR=1,NGRP
         IF(NPSYS(IGR).NE.0) THEN
            KPSTR=LCMGIL(JPSTR,IGR)
            CALL LCMLEN(KPSTR,'DRAGON-TXSC',ILCTXS,ITYLCM)
            ALLOCATE(SIGT0(0:ILCTXS-1))
            CALL LCMGET(KPSTR,'DRAGON-TXSC',SIGT0(0))
            ZNUM=0.0
            ZDEN=0.0
            DO 100 IR=1,NREG
               IBM=MATCOD(IR)
               IND=KEYFLX(IR,1)
               ZNUM=ZNUM+SIGT0(IBM)*FLUX(IND,IGR,2)*VOL(IR)
               ZDEN=ZDEN+FLUX(IND,IGR,2)*VOL(IR)
  100       CONTINUE
            DO 110 IR=1,NREG
            IBM=MATCOD(IR)
            IND=KEYFLX(IR,1)
            FLUX(NUNKNO/2+IND,IGR,3)=FLUX(NUNKNO/2+IND,IGR,3)+
     1      (1.0-GAMMA(IGR))*(ZNUM/ZDEN-SIGT0(IBM))*
     2      FLUX(NUNKNO/2+IND,IGR,2)
  110       CONTINUE
            DEALLOCATE(SIGT0)
            DO 120 IR=1,NREG
            IND=KEYFLX(IR,1)
            FLUX(NUNKNO/2+IND,IGR,3)=(FLUX(NUNKNO/2+IND,IGR,3)
     1      +FLUX(IND,IGR,2)/3.0)/GAMMA(IGR)
  120       CONTINUE
         ENDIF
  130    CONTINUE
         IF(IPRT.GE.3) WRITE(IOUT,'(30H FLUDBV: FUNDAMENTAL CURRENTS.)')
         ALLOCATE(SUNKNO((NUNKNO/2)*NGRP),FUNKNO((NUNKNO/2)*NGRP))
         IOF=0
         DO 145 IGR=1,NGRP
         DO 140 IND=1,NUNKNO/2
         IOF=IOF+1
         SUNKNO(IOF)=FLUX(NUNKNO/2+IND,IGR,3)
         FUNKNO(IOF)=FLUX(NUNKNO/2+IND,IGR,2)
  140    CONTINUE
  145    CONTINUE
         IDIR=0
         CALL DOORFVR(CDOOR,JPSTR,NPSYS,IPTRK,IFTRAK,IPRT,NGRP,
     1        NMAT,IDIR,NREG,NUNKNO/2,IPHASE,LEXAC,MATCOD,VOL,KEYFLX,
     2        TITLE,SUNKNO,FUNKNO,IPMACR,REBFLG)
         IOF=0
         DO 155 IGR=1,NGRP
         DO 150 IND=1,NUNKNO/2
         IOF=IOF+1
         FLUX(NUNKNO/2+IND,IGR,3)=SUNKNO(IOF)
         FLUX(NUNKNO/2+IND,IGR,2)=FUNKNO(IOF)
  150    CONTINUE
  155    CONTINUE
         DEALLOCATE(FUNKNO,SUNKNO)
      ELSE IF((MOD(ILEAK,10).EQ.6).AND.(IPHASE.EQ.1)) THEN
*        ----
*        TIBERE ANISOTROPIC STREAMING MODEL FOR MOC.
*        ----
         IF(.NOT.C_ASSOCIATED(JPSTR)) THEN
            CALL XABORT('FLUDBV: MISSING STREAMING INFO(2).')
         ENDIF
* ADD SOURCES FOR FLUX EQUATION
         DO IGR=1,NGRP
           IF(NPSYS(IGR).NE.0) THEN 
             IF((B2(1).NE.0.0).AND.(B2(2).NE.0.0).AND.
     1       (B2(3).NE.0.0)) THEN
               S=0.0
               DO IR=1,NREG 
                 IND=KEYFLX(IR,1)
                 INDC(1)=3*NUNKNO/8+IND
                 INDC(2)=5*NUNKNO/8+IND
                 INDC(3)=7*NUNKNO/8+IND
                 FLUX(IND,IGR,3)=FLUX(IND,IGR,3)-(B2(1)
     1           *FLUX(INDC(1),IGR,1)+B2(2)*FLUX(INDC(2),IGR,1)+
     2           B2(3)*FLUX(INDC(3),IGR,1))
               ENDDO 
             ENDIF  
           ENDIF 
         ENDDO
         IF(IPRT.GE.3) WRITE(IOUT,'(28H FLUDBV: FUNDAMENTAL FLUXES.)')
           IDIR=0
           NUNKNO4=NUNKNO/4
           ALLOCATE(SUNKNO(NUNKNO4*NGRP),FUNKNO(NUNKNO4*NGRP))
           IOF=0
           DO IGR=1,NGRP
             DO IND=1,NUNKNO4
               IOF=IOF+1
               SUNKNO(IOF)=FLUX(IND,IGR,3)
               FUNKNO(IOF)=FLUX(IND,IGR,2)
             ENDDO
           ENDDO
           CALL DOORFVR(CDOOR,JPSYS,NPSYS,IPTRK,IFTRAK,IPRT,NGRP,
     1        NMAT,IDIR,NREG,NUNKNO4,IPHASE,LEXAC,MATCOD,VOL,KEYFLX,
     2        TITLE,SUNKNO,FUNKNO,IPMACR,REBFLG)
           IOF=0
           DO IGR=1,NGRP
             DO IND=1,NUNKNO4
               IOF=IOF+1
               FLUX(IND,IGR,3)=SUNKNO(IOF)
               FLUX(IND,IGR,2)=FUNKNO(IOF) 
             ENDDO
           ENDDO
           DEALLOCATE(FUNKNO,SUNKNO)   
* ADD SOURCES FOR CURRENT EQUATIONS
           DO IGR=1,NGRP
             IF(NPSYS(IGR).NE.0) THEN 
             KPSTR=LCMGIL(JPSTR,IGR)
             CALL LCMLEN(KPSTR,'DRAGON-TXSC',ILCTXS,ITYLCM)
             ALLOCATE(SIGT0(0:ILCTXS-1))
             CALL LCMGET(KPSTR,'DRAGON-TXSC',SIGT0(0))
             ZNUM=0.0
             ZDEN=0.0
             DO IR=1,NREG
               IBM=MATCOD(IR)
               IND=KEYFLX(IR,1)
               ZNUM=ZNUM+SIGT0(IBM)*FLUX(IND,IGR,2)*VOL(IR)
               ZDEN=ZDEN+FLUX(IND,IGR,2)*VOL(IR)
             ENDDO
             DO IR=1,NREG
               IBM=MATCOD(IR)
               INDD(1)=NUNKNO/4+KEYFLX(IR,1)
               INDD(2)=NUNKNO/2+KEYFLX(IR,1)
               INDD(3)=3*NUNKNO/4+KEYFLX(IR,1)
               DO IDIR=1,3
                 FLUX(INDD(IDIR),IGR,3)=FLUX(INDD(IDIR),IGR,3)+
     1           (1.0-GAMMA(IGR))*(ZNUM/ZDEN-SIGT0(IBM))*
     2            FLUX(INDD(IDIR),IGR,2)
               ENDDO
             ENDDO
             DO IR=1,NREG
               INDD(1)=NUNKNO/4+KEYFLX(IR,1)
               INDD(2)=NUNKNO/2+KEYFLX(IR,1)
               INDD(3)=3*NUNKNO/4+KEYFLX(IR,1)
               DO IDIR=1,3
                 FLUX(INDD(IDIR),IGR,3)=(FLUX(INDD(IDIR),IGR,3)
     1           +FLUX(KEYFLX(IR,1),IGR,2)/3.0)/GAMMA(IGR)
               ENDDO
             ENDDO
             DEALLOCATE(SIGT0)
           ENDIF
           ENDDO
         DO IDIR=1,3
           IF(IPRT.GE.3) 
     >     WRITE(IOUT,'(30H FLUDBV: FUNDAMENTAL CURRENTS.)') 
           IF(IDIR.EQ.1) WRITE(6,*)'FUNDAMENTAL CURRENT X '
           IF(IDIR.EQ.2) WRITE(6,*)'FUNDAMENTAL CURRENT Y '
           IF(IDIR.EQ.3) WRITE(6,*)'FUNDAMENTAL CURRENT Z '
           NUNKNO4=NUNKNO/4
           ALLOCATE(SUNKNO(NUNKNO4*NGRP),FUNKNO(NUNKNO4*NGRP))
           IOF=0
           DO IGR=1,NGRP
             DO IND=1,NUNKNO4
               INDB(1)=NUNKNO/4+IND
               INDB(2)=NUNKNO/2+IND
               INDB(3)=3*NUNKNO/4+IND
               IOF=IOF+1
               SUNKNO(IOF)=FLUX(INDB(IDIR),IGR,3)
               FUNKNO(IOF)=FLUX(INDB(IDIR),IGR,2)
             ENDDO
           ENDDO
           CALL DOORFVR(CDOOR,JPSTR,NPSYS,IPTRK,IFTRAK,IPRT,NGRP,
     1       NMAT,IDIR,NREG,NUNKNO4,IPHASE,LEXAC,MATCOD,VOL,KEYFLX,
     2       TITLE,SUNKNO,FUNKNO,IPMACR,REBFLG)
           IOF=0
           DO IGR=1,NGRP
             DO IND=1,NUNKNO4
               INDB(1)=NUNKNO/4+IND
               INDB(2)=NUNKNO/2+IND
               INDB(3)=3*NUNKNO/4+IND
               IOF=IOF+1
               FLUX(INDB(IDIR),IGR,3)=SUNKNO(IOF)
               FLUX(INDB(IDIR),IGR,2)=FUNKNO(IOF)   
             ENDDO
           ENDDO
           DEALLOCATE(FUNKNO,SUNKNO)  
         ENDDO
      ELSE IF((MOD(ILEAK,10).EQ.6).AND.(IPHASE.EQ.2)) THEN
*        ----
*        TIBERE ANISOTROPIC STREAMING MODEL FOR PIJ.
*        ----
         INDD(1)=NUNKNO/4
         INDD(2)=NUNKNO/2
         INDD(3)=3*NUNKNO/4
         NUN4=NUNKNO/4
         ALLOCATE(PP(NREG*NREG))
         DO 210 IGR=1,NGRP
         IF(NPSYS(IGR).NE.0) THEN
           KPSYS=LCMGIL(JPSYS,IGR)
           DO 200 IDIR=1,3
           IF(B2(IDIR).NE.0.0) THEN
             WRITE(TEXT12,'(6HDRAGON,I1,5HP*SCT)') IDIR
             CALL LCMGET(KPSYS,TEXT12,PP)
             DO 190 IREG=1,NREG
             IND=KEYFLX(IREG,1)
             S=0.0
             DO 180 JREG=1,NREG
             JND=KEYFLX(JREG,1)
             S=S+FLUX(INDD(IDIR)+JND,IGR,1)*PP((JREG-1)*NREG+IREG)
 180         CONTINUE
             FLUX(IND,IGR,3)=FLUX(IND,IGR,3)-B2(IDIR)*S
 190         CONTINUE
           ENDIF
 200     CONTINUE
         ENDIF
 210     CONTINUE
         DEALLOCATE(PP)
         IDIR=0
         CALL DOORFVR(CDOOR,JPSYS,NPSYS,IPTRK,IFTRAK,IPRT,NGRP,
     1        NMAT,IDIR,NREG,NUNKNO,IPHASE,LEXAC,MATCOD,VOL,KEYFLX,
     2        TITLE,FLUX(1,1,3),FLUX(1,1,2),IPMACR,REBFLG)
         DO 260 IDIR=1,3
         DO 250 IGR=1,NGRP
         IF(NPSYS(IGR).NE.0) THEN
           KPSYS=LCMGIL(JPSYS,IGR)
           CALL LCMLEN(KPSYS,'DRAGON-TXSC',ILCTXS,ITYLCM)
           ALLOCATE(SIGT0(0:ILCTXS-1))
           CALL LCMGET(KPSYS,'DRAGON-TXSC',SIGT0(0))
           ZNUM=0.0
           ZDEN=0.0
           DO 220 IREG=1,NREG
             IBM=MATCOD(IREG)
             IND=KEYFLX(IREG,1)
             ZNUM=ZNUM+SIGT0(IBM)*FLUX(IND,IGR,1)*VOL(IREG)
             ZDEN=ZDEN+FLUX(IND,IGR,1)*VOL(IREG)
 220       CONTINUE
           DO 230 IREG=1,NREG
           IBM=MATCOD(IREG)
           IND=KEYFLX(IREG,1)
           IND2=INDD(IDIR)+IND
           FLUX(IND2,IGR,3)=FLUX(IND2,IGR,3)+(1.0-GAMMA(IGR))*
     1     (ZNUM/ZDEN-SIGT0(IBM))*FLUX(IND2,IGR,1)
 230       CONTINUE
           DEALLOCATE(SIGT0)
           DO 240 IND=1,NUN4
           IND2=INDD(IDIR)+IND
           FLUX(IND2,IGR,3)=(FLUX(IND2,IGR,3)+FLUX(IND,IGR,2)/3.0)/
     1     GAMMA(IGR)
 240       CONTINUE
         ENDIF
 250     CONTINUE
         CALL DOORFVR(CDOOR,JPSYS,NPSYS,IPTRK,IFTRAK,IPRT,NGRP,
     1        NMAT,IDIR,NREG,NUNKNO,IPHASE,LEXAC,MATCOD,VOL,KEYFLX,
     2        TITLE,FLUX(INDD(IDIR)+1,1,3),FLUX(INDD(IDIR)+1,1,2),
     3        IPMACR,REBFLG)
 260     CONTINUE
      ELSE
         CALL XABORT('FLUDBV: TYPE OF LEAKAGE NOT IMPLEMENTED.')
      ENDIF
*----
*  COMPUTE DB2 PARAMETER CORRESPONDING TO ACTUAL LEAKAGE
*----
      IF(IPRT.GT.10) THEN
        NUN=NUNKNO
        IF(ILEAK.EQ.5) NUN=NUNKNO/2
        IF(ILEAK.GE.6) NUN=NUNKNO/4
        DO 270 IGR=1,NGRP
        IF(NPSYS(IGR).EQ.0) GO TO 270
        KPSYS=LCMGIL(JPSYS,IGR)
        DB2NEW=FLUFUI(KPSYS,NREG,NUN,MATCOD,VOL,KEYFLX,FLUX(1,IGR,2),
     >                SUNKN(1,IGR))
        WRITE(IOUT,'(15H FLUDBV: GROUP=,I5,24H DB2 LEAKAGE PARAMETER F,
     >  12HROM DIFFON =,1P,E13.4/26X,30HACTUAL DB2 LEAKAGE PARAMETER =,
     >  E13.4)') IGR,DDD(IGR)*B2(4),DB2NEW
  270   CONTINUE
      ENDIF
      DEALLOCATE(SUNKN)
      RETURN
      END
