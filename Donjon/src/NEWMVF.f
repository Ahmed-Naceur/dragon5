*DECK NEWMVF
      SUBROUTINE NEWMVF(INDX,DPOS,DMIX,NGRP,NL,NDEL,LEAK,NEL,NMIX,LX,
     1  LY,LZ,MESHX,MESHY,MESHZ,NTOT0,NTOT1,ZNUS,CHI,ZSIGF,DIFFX,DIFFY,
     2  DIFFZ,HFAC,SCAT,XFAC,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover material regions affected by the device insertion.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal.
*
*Author(s): 
* J. Koclas, D. Sekki
*
*Parameters: input/output
* INDX    index number of each material volume (=0 for virtual regions).
* DPOS    device position in cm in the core.
* DMIX    device mixtures for insertion and extraction.
* NGRP    number of energy  groups.
* NL      number of legendre orders (=1 for isotropic scattering).
* NDEL    number of precursor groups for delayed neutron.
* LEAK    diffusion coefficient flag (=1: isotropic; =2: anisotropic).
* NEL     total number of elements.
* NMIX    maximum number of material mixtures.
* LX      number of elements along x-axis.
* LY      number of elements along y-axis.
* LZ      number of elements along z-axis.
* MESHX   mesh coordinates along x-axis.
* MESHY   mesh coordinates along y-axis.
* MESHZ   mesh coordinates along z-axis.
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
* XFAC    corrective factor for delta sigmas.
* IMPX    printing index (=0 for no print).
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NGRP,NL,NDEL,LEAK,NEL,NMIX,INDX(NEL),DMIX(2),LX,LY,LZ,IMPX
      REAL MESHX(LX+1),MESHY(LY+1),MESHZ(LZ+1),DIFFX(NMIX,NGRP),
     1 ZSIGF(NMIX,NGRP),NTOT1(NMIX,NGRP),ZNUS(NMIX,NGRP,NDEL+1),
     2 CHI(NMIX,NGRP,NDEL+1),DPOS(6),NTOT0(NMIX,NGRP),HFAC(NMIX,NGRP),
     3 SCAT(NMIX,NL,NGRP,NGRP),DIFFY(NMIX,NGRP),DIFFZ(NMIX,NGRP),XFAC
      PARAMETER(IOUT=6,EPSI=1.0E-4)
*----
*  RECOVER REGIONS WHERE DEVICE IS INSERTED
*----
      IF(IMPX.GT.4)WRITE(IOUT,*)'RECOVER REGIONS AFFECTED BY DEVICE'
*     INSERTED COORDINATES
      DX1=DPOS(1)
      DX2=DPOS(2)
      DY1=DPOS(3)
      DY2=DPOS(4)
      DZ1=DPOS(5)
      DZ2=DPOS(6)
      IF(DX1.LT.MESHX(1)) DX1=MESHX(1)
      IF(DX2.LT.MESHX(1)) DX2=MESHX(1)
      IF(DX2.GT.MESHX(LX+1)) DX2=MESHX(LX+1)
      IF(DX1.GT.MESHX(LX+1)) DX1=MESHX(LX+1)
      IF(ABS(DX1-DX2).LT.EPSI) RETURN
      IF(DY1.LT.MESHY(1)) DY1=MESHY(1)
      IF(DY2.LT.MESHY(1)) DY2=MESHY(1)
      IF(DY2.GT.MESHY(LY+1)) DY2=MESHY(LY+1)
      IF(DY1.GT.MESHY(LY+1)) DY1=MESHY(LY+1)
      IF(ABS(DY1-DY2).LT.EPSI) RETURN
      IF(DZ1.LT.MESHZ(1)) DZ1=MESHZ(1)
      IF(DZ2.LT.MESHZ(1)) DZ2=MESHZ(1)
      IF(DZ2.GT.MESHZ(LZ+1)) DZ2=MESHZ(LZ+1)
      IF(DZ1.GT.MESHZ(LZ+1)) DZ1=MESHZ(LZ+1)
      IF(ABS(DZ1-DZ2).LT.EPSI) RETURN
      I1=0
      I2=0
*     CHECK X-AXIS
      DO I=1,LX
      IF(ABS(DX1-MESHX(I)).LT.EPSI) DX1=MESHX(I)
      IF(ABS(DX2-MESHX(I)).LT.EPSI) DX2=MESHX(I)
      IF(ABS(DX1-MESHX(I+1)).LT.EPSI) DX1=MESHX(I+1)
      IF(ABS(DX2-MESHX(I+1)).LT.EPSI) DX2=MESHX(I+1)
      IF((DX1.GE.MESHX(I)).AND.(DX1.LT.MESHX(I+1)))I1=I
      IF((DX2.GT.MESHX(I)).AND.(DX2.LE.MESHX(I+1)))THEN
        I2=I
        GOTO 10
      ENDIF
      ENDDO
   10 IF(IMPX.GT.4)WRITE(IOUT,*)' I1=',I1,', I2=',I2
      IF((I1.EQ.0).OR.(I2.EQ.0))CALL XABORT('@NEWMVF: WR'
     1 //'ONG NUMBER OF AFFECTED REGIONS ALONG X-AXIS.')
      J1=0
      J2=0
*     CHECK Y-AXIS
      DO J=1,LY
      IF(ABS(DY1-MESHY(J)).LT.EPSI) DY1=MESHY(J)
      IF(ABS(DY2-MESHY(J)).LT.EPSI) DY2=MESHY(J)
      IF(ABS(DY1-MESHY(J+1)).LT.EPSI) DY1=MESHY(J+1)
      IF(ABS(DY2-MESHY(J+1)).LT.EPSI) DY2=MESHY(J+1)
      IF((DY1.GE.MESHY(J)).AND.(DY1.LT.MESHY(J+1)))J1=J
      IF((DY2.GT.MESHY(J)).AND.(DY2.LE.MESHY(J+1)))THEN
        J2=J
        GOTO 20
      ENDIF
      ENDDO
   20 IF(IMPX.GT.4)WRITE(IOUT,*)' J1=',J1,', J2=',J2
      IF((J1.EQ.0).OR.(J2.EQ.0))CALL XABORT('@NEWMVF: WR'
     1 //'ONG NUMBER OF AFFECTED REGIONS ALONG Y-AXIS.')
      K1=0
      K2=0
*     CHECK Z-AXIS
      DO K=1,LZ
      IF(ABS(DZ1-MESHZ(K)).LT.EPSI) DZ1=MESHZ(K)
      IF(ABS(DZ2-MESHZ(K)).LT.EPSI) DZ2=MESHZ(K)
      IF(ABS(DZ1-MESHZ(K+1)).LT.EPSI) DZ1=MESHZ(K+1)
      IF(ABS(DZ2-MESHZ(K+1)).LT.EPSI) DZ2=MESHZ(K+1)
      IF((DZ1.GE.MESHZ(K)).AND.(DZ1.LT.MESHZ(K+1)))K1=K
      IF((DZ2.GT.MESHZ(K)).AND.(DZ2.LE.MESHZ(K+1)))THEN
        K2=K
        GOTO 30
      ENDIF
      ENDDO
   30 IF(IMPX.GT.4)WRITE(IOUT,*)' K1=',K1,', K2=',K2
      IF((K1.EQ.0).OR.(K2.EQ.0))CALL XABORT('@NEWMVF: WR'
     1 //'ONG NUMBER OF AFFECTED REGIONS ALONG Z-AXIS.')
*----
*  COMPUTE OCCUPIED VOLUME FRACTION
*----
      DO 42 K=K1,K2
      DO 41 J=J1,J2
      DO 40 I=I1,I2
      IEL=(K-1)*LX*LY+(J-1)*LX+I
      IBM=INDX(IEL)
      IF(IMPX.GT.4)WRITE(IOUT,*)'AFFECTED ELEM #',IEL,' MIX #',IBM
      IF(IBM.NE.0)THEN
        FX=0.
*       FRACTION ALONG X-AXIS
        IF((DX1.GE.MESHX(I)).AND.(DX2.GT.MESHX(I+1)))THEN
          FX=(MESHX(I+1)-DX1)/(MESHX(I+1)-MESHX(I))
        ELSEIF((DX1.GE.MESHX(I)).AND.(DX2.LE.MESHX(I+1)))THEN
          FX=(DX2-DX1)/(MESHX(I+1)-MESHX(I))
        ELSEIF((DX1.LT.MESHX(I)).AND.(DX2.GT.MESHX(I+1)))THEN
          FX=1.
        ELSEIF((DX1.LT.MESHX(I)).AND.(DX2.LE.MESHX(I+1)))THEN
          FX=(DX2-MESHX(I))/(MESHX(I+1)-MESHX(I))
        ENDIF
        FY=0.
*       FRACTION ALONG Y-AXIS
        IF((DY1.GE.MESHY(J)).AND.(DY2.GT.MESHY(J+1)))THEN
          FY=(MESHY(J+1)-DY1)/(MESHY(J+1)-MESHY(J))
        ELSEIF((DY1.GE.MESHY(J)).AND.(DY2.LE.MESHY(J+1)))THEN
          FY=(DY2-DY1)/(MESHY(J+1)-MESHY(J))
        ELSEIF((DY1.LT.MESHY(J)).AND.(DY2.GT.MESHY(J+1)))THEN
          FY=1.
        ELSEIF((DY1.LT.MESHY(J)).AND.(DY2.LE.MESHY(J+1)))THEN
          FY=(DY2-MESHY(J))/(MESHY(J+1)-MESHY(J))
        ENDIF
        FZ=0.
*       FRACTION ALONG Z-AXIS
        IF((DZ1.GE.MESHZ(K)).AND.(DZ2.GT.MESHZ(K+1)))THEN
          FZ=(MESHZ(K+1)-DZ1)/(MESHZ(K+1)-MESHZ(K))
        ELSEIF((DZ1.GE.MESHZ(K)).AND.(DZ2.LE.MESHZ(K+1)))THEN
          FZ=(DZ2-DZ1)/(MESHZ(K+1)-MESHZ(K))
        ELSEIF((DZ1.LT.MESHZ(K)).AND.(DZ2.GT.MESHZ(K+1)))THEN
          FZ=1.
        ELSEIF((DZ1.LT.MESHZ(K)).AND.(DZ2.LE.MESHZ(K+1)))THEN
          FZ=(DZ2-MESHZ(K))/(MESHZ(K+1)-MESHZ(K))
        ENDIF
*       VOLUME FRACTION
        VF=FX*FY*FZ
        IF((IMPX.GT.4).AND.(VF.GT.EPSI))
     1    WRITE(IOUT,*)'INSERTED DEVICE VOLUME FRACTION ',VF
*       UPDATE PROPERTIES
        IF(VF.GT.EPSI)
     1    CALL NEWMXS(NTOT0,NTOT1,ZNUS,CHI,ZSIGF,DIFFX,DIFFY,DIFFZ,HFAC,
     2    SCAT,IBM,DMIX(1),DMIX(2),NGRP,NMIX,NL,NDEL,LEAK,VF,XFAC,IMPX)
      ENDIF
   40 CONTINUE
   41 CONTINUE
   42 CONTINUE
      RETURN
      END
