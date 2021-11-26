*DECK NEWMXS
      SUBROUTINE NEWMXS(NTOT0,NTOT1,ZNUS,CHI,ZSIGF,DIFFX,DIFFY,DIFFZ,
     1 HFAC,SCAT,IBM,IBM1,IBM2,NGRP,NMIX,NL,NDEL,LEAK,VF,XFAC,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute nuclear properties perturbed by the device insertion.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* D. Sekki
*
*Parameters: input/output
* NTOT0   flux-weighted total macroscopic x-sections.
* NTOT1   current-weighted total macroscopic x-sections.
* ZNUS    nu*fission macroscopic x-sections.
* CHI     fission spectra.
* ZSIGF   fission macroscopic x-sections.
* DIFFX   x-directed diffusion coefficients.
* DIFFY   y-directed diffusion coefficients.
* DIFFZ   z-directed diffusion coefficients.
* HFAC    h-factors (kappa*fission macroscopic x-sections).
* SCAT    scattering macroscopic x-sections.
* IBM     mixture index for physical region.
* IBM1    device mixture index for inserted device.
* IBM2    device mixture index for extracted device.
* NGRP    number of energy  groups.
* NMIX    maximum number of material mixtures.
* NL      number of legendre orders (=1 for isotropic scattering).
* NDEL    number of precursor groups for delayed neutron.
* LEAK    diffusion coefficient flag (=1: isotropic; =2: anisotropic).
* VF      volume fraction occupied by the device.
* XFAC    corrective factor for delta sigmas.
* IMPX    printing index (=0 for no print).
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IBM,IBM1,IBM2,NGRP,NL,NDEL,LEAK,NMIX,IMPX
      REAL NTOT0(NMIX,NGRP),NTOT1(NMIX,NGRP),ZNUS(NMIX,NGRP,NDEL+1),
     1     CHI(NMIX,NGRP,NDEL+1),ZSIGF(NMIX,NGRP),DIFFX(NMIX,NGRP),
     2     HFAC(NMIX,NGRP),VF,SCAT(NMIX,NL,NGRP,NGRP),DIFFY(NMIX,NGRP),
     3     DIFFZ(NMIX,NGRP),XFAC
      PARAMETER(IOUT=6,EPSI=1.0E-4)
*----
*  UPDATE PROPERTIES
*----
      IF(IMPX.GT.4)WRITE(IOUT,*)' UPDATING PROPERTIES'
      DO 70 JGR=1,NGRP
      IF(IMPX.GT.4)WRITE(IOUT,*)' '
      IF(IMPX.GT.4)WRITE(IOUT,*)' PROCESSING ENERGY GROUP # ',JGR
*----
*  NTOT0
*----
      IF(IMPX.GT.4)WRITE(IOUT,*)' NTOT0 BEFORE  : ',NTOT0(IBM,JGR)
      DELT=NTOT0(IBM1,JGR)-NTOT0(IBM2,JGR)
      NTOT0(IBM,JGR)=NTOT0(IBM,JGR)+XFAC*VF*DELT
      IF(IMPX.GT.4)WRITE(IOUT,*)' DELTA-SIGMA   : ',DELT
      IF(IMPX.GT.4)WRITE(IOUT,*)' NTOT0 AFTER   : ',NTOT0(IBM,JGR)
*----
*  NTOT1
*----
      IF(NTOT1(IBM,JGR).EQ.0.)GOTO 10
      IF(IMPX.GT.4)WRITE(IOUT,*)' NTOT1 BEFORE  : ',NTOT1(IBM,JGR)
      DELT=NTOT1(IBM1,JGR)-NTOT1(IBM2,JGR)
      NTOT1(IBM,JGR)=NTOT1(IBM,JGR)+XFAC*VF*DELT
      IF(IMPX.GT.4)WRITE(IOUT,*)' DELTA-SIGMA   : ',DELT
      IF(IMPX.GT.4)WRITE(IOUT,*)' NTOT1 AFTER   : ',NTOT1(IBM,JGR)
*----
*  NUSIGF
*----
   10 IF(ZNUS(IBM,JGR,1).EQ.0.)GOTO 15
      IF(IMPX.GT.4)WRITE(IOUT,*)' NUSIGF BEFORE : ',ZNUS(IBM,JGR,1)
      DELT=ZNUS(IBM1,JGR,1)-ZNUS(IBM2,JGR,1)
      ZNUS(IBM,JGR,1)=ZNUS(IBM,JGR,1)+XFAC*VF*DELT
      IF(IMPX.GT.4)WRITE(IOUT,*)' DELTA-SIGMA   : ',DELT
      IF(IMPX.GT.4)WRITE(IOUT,*)' NUSIGF AFTER  : ',ZNUS(IBM,JGR,1)
      DO IDEL=1,NDEL
        IF(IMPX.GT.4)WRITE(IOUT,*)' NUSIGF',IDEL,' BEFORE : ',
     >  ZNUS(IBM,JGR,IDEL+1)
        DELT=ZNUS(IBM1,JGR,IDEL+1)-ZNUS(IBM2,JGR,IDEL+1)
        ZNUS(IBM,JGR,IDEL+1)=ZNUS(IBM,JGR,IDEL+1)+XFAC*VF*DELT
        IF(IMPX.GT.4)WRITE(IOUT,*)' DELTA-SIGMA   : ',DELT
        IF(IMPX.GT.4)WRITE(IOUT,*)' NUSIGF',IDEL,' AFTER  : ',
     >  ZNUS(IBM,JGR,IDEL+1)
      ENDDO
*----
*  CHI
*----
   15 IF(CHI(IBM,JGR,1).EQ.0.)GOTO 20
      IF(IMPX.GT.4)WRITE(IOUT,*)' CHI BEFORE : ',CHI(IBM,JGR,1)
      DELT=CHI(IBM1,JGR,1)-CHI(IBM2,JGR,1)
      CHI(IBM,JGR,1)=CHI(IBM,JGR,1)+XFAC*VF*DELT
      IF(IMPX.GT.4)WRITE(IOUT,*)' DELTA-SIGMA   : ',DELT
      IF(IMPX.GT.4)WRITE(IOUT,*)' CHI AFTER  : ',CHI(IBM,JGR,1)
      DO IDEL=1,NDEL
        IF(IMPX.GT.4)WRITE(IOUT,*)' CHI',IDEL,' BEFORE : ',
     >  CHI(IBM,JGR,IDEL+1)
        DELT=CHI(IBM1,JGR,IDEL+1)-CHI(IBM2,JGR,IDEL+1)
        CHI(IBM,JGR,IDEL+1)=CHI(IBM,JGR,IDEL+1)+XFAC*VF*DELT
        IF(IMPX.GT.4)WRITE(IOUT,*)' DELTA-SIGMA   : ',DELT
        IF(IMPX.GT.4)WRITE(IOUT,*)' CHI',IDEL,' AFTER  : ',
     >  CHI(IBM,JGR,IDEL+1)
      ENDDO
*----
*  NFTOT
*----
   20 IF(ZSIGF(IBM,JGR).EQ.0.)GOTO 30
      IF(IMPX.GT.4)WRITE(IOUT,*)' NFTOT BEFORE  : ',ZSIGF(IBM,JGR)
      DELT=ZSIGF(IBM1,JGR)-ZSIGF(IBM2,JGR)
      ZSIGF(IBM,JGR)=ZSIGF(IBM,JGR)+XFAC*VF*DELT
      IF(IMPX.GT.4)WRITE(IOUT,*)' DELTA-SIGMA   : ',DELT
      IF(IMPX.GT.4)WRITE(IOUT,*)' NFTOT AFTER   : ',ZSIGF(IBM,JGR)
*----
*  HFAC
*----
   30 IF(HFAC(IBM,JGR).EQ.0.)GOTO 40
      IF(IMPX.GT.4)WRITE(IOUT,*)' H-FACTOR BEFORE : ',HFAC(IBM,JGR)
      DELT=HFAC(IBM1,JGR)-HFAC(IBM2,JGR)
      HFAC(IBM,JGR)=HFAC(IBM,JGR)+XFAC*VF*DELT
      IF(IMPX.GT.4)WRITE(IOUT,*)' DELTA-SIGMA    : ',DELT
      IF(IMPX.GT.4)WRITE(IOUT,*)' H-FACTOR AFTER : ',HFAC(IBM,JGR)
*----
*  DIFFX
*----
   40 IF(DIFFX(IBM,JGR).LT.EPSI)GOTO 50
      IF(IMPX.GT.4)WRITE(IOUT,*)' DIFFX BEFORE  : ',DIFFX(IBM,JGR)
      IF(DIFFX(IBM1,JGR).LT.EPSI)THEN
        DEL1=0.
      ELSE
        DEL1=1./DIFFX(IBM1,JGR)
      ENDIF
      IF(DIFFX(IBM2,JGR).LT.EPSI)THEN
        DEL2=0.
      ELSE
        DEL2=1./DIFFX(IBM2,JGR)
      ENDIF
      DELT=DEL1-DEL2
      DIFFX(IBM,JGR)=1./DIFFX(IBM,JGR)+XFAC*VF*DELT
      DIFFX(IBM,JGR)=1./DIFFX(IBM,JGR)
      IF(IMPX.GT.4)WRITE(IOUT,*)' DELTA-SIGMA   : ',DELT
      IF(IMPX.GT.4)WRITE(IOUT,*)' DIFFX AFTER   : ',DIFFX(IBM,JGR)
*----
*  DIFFY
*----
   50 IF(LEAK.NE.2)GOTO 70
      IF(DIFFY(IBM,JGR).LT.EPSI)GOTO 60
      IF(IMPX.GT.4)WRITE(IOUT,*)' DIFFY BEFORE  : ',DIFFY(IBM,JGR)
      IF(DIFFY(IBM1,JGR).LT.EPSI)THEN
        DEL1=0.
      ELSE
        DEL1=1./DIFFY(IBM1,JGR)
      ENDIF
      IF(DIFFY(IBM2,JGR).LT.EPSI)THEN
        DEL2=0.
      ELSE
        DEL2=1./DIFFY(IBM2,JGR)
      ENDIF
      DELT=DEL1-DEL2
      DIFFY(IBM,JGR)=1./DIFFY(IBM,JGR)+XFAC*VF*DELT
      DIFFY(IBM,JGR)=1./DIFFY(IBM,JGR)
      IF(IMPX.GT.4)WRITE(IOUT,*)' DELTA-SIGMA   : ',DELT
      IF(IMPX.GT.4)WRITE(IOUT,*)' DIFFY AFTER   : ',DIFFY(IBM,JGR)
*----
*  DIFFZ
*----
   60 IF(DIFFZ(IBM,JGR).LT.EPSI)GOTO 70
      IF(IMPX.GT.4)WRITE(IOUT,*)' DIFFZ BEFORE  : ',DIFFZ(IBM,JGR)
      IF(DIFFZ(IBM1,JGR).LT.EPSI)THEN
        DEL1=0.
      ELSE
        DEL1=1./DIFFZ(IBM1,JGR)
      ENDIF
      IF(DIFFZ(IBM2,JGR).LT.EPSI)THEN
        DEL2=0.
      ELSE
        DEL2=1./DIFFZ(IBM2,JGR)
      ENDIF
      DELT=DEL1-DEL2
      DIFFZ(IBM,JGR)=1./DIFFZ(IBM,JGR)+XFAC*VF*DELT
      DIFFZ(IBM,JGR)=1./DIFFZ(IBM,JGR)
      IF(IMPX.GT.4)WRITE(IOUT,*)' DELTA-SIGMA   : ',DELT
      IF(IMPX.GT.4)WRITE(IOUT,*)' DIFFZ AFTER   : ',DIFFZ(IBM,JGR)
   70 CONTINUE
*----
*  SCAT
*----
      DO 82 IL=1,NL
      DO 81 IGR=1,NGRP
      DO 80 JGR=1,NGRP
      DELT=SCAT(IBM1,IL,IGR,JGR)-SCAT(IBM2,IL,IGR,JGR)
      IF((SCAT(IBM,IL,IGR,JGR).NE.0.0).AND.(DELT.NE.0.0)) THEN
        IF(IMPX.GT.4)WRITE(IOUT,*)' PROCESSING ENERGY GROUP # ',JGR,
     >  '<-',IGR
        IF(IMPX.GT.4)WRITE(IOUT,*)' SCAT',IL,' BEFORE   : ',
     >  SCAT(IBM,IL,IGR,JGR)
      ENDIF
      SCAT(IBM,IL,IGR,JGR)=SCAT(IBM,IL,IGR,JGR)+XFAC*VF*DELT
      IF(DELT.NE.0.0) THEN
        IF(IMPX.GT.4)WRITE(IOUT,*)' DELTA-SIGMA   : ',DELT
        IF(IMPX.GT.4)WRITE(IOUT,*)' SCAT',IL,' AFTER    : ',
     >  SCAT(IBM,IL,IGR,JGR)
      ENDIF
   80 CONTINUE
   81 CONTINUE
   82 CONTINUE
      IF(IMPX.GT.4)WRITE(IOUT,*)' ALL PROPERTIES UPDATED.'
      RETURN
      END
