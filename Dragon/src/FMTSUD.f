*DECK FMTSUD
      SUBROUTINE FMTSUD(NENTRY,IENTRY,KENTRY,IPRINT,NOPT,IOPT,IKFLU,
     >                  NREG  ,NGROUP,NDIM  ,POLOAQ,AZMOAQ,
     >                  VOLUME,XPOL  ,WPOL  ,XAZI  ,WAZI  ,
     >                  FLUX  ,AFLUX ,TFLUX ,WGHT  ,MU    ,ETA   )
*
*-----------------------------------------------------------------------
*
*Purpose:
* To process the angular fluxes and generate the SUS3D file.
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*
*Author(s):
* G. Marleau
*
*Parameters: input
* NENTRY  number of data structures transfered to this module.
* IENTRY  data structure type where:
*         =1 for LCM memory object;
*         =2 for XSM file;
*         =3 for sequential binary file;
*         =4 for sequential ASCII file.
* KENTRY  data structure pointer.
* IPRINT  print level.
* NOPT    number of options.
* IOPT    processing option.
* IKFLU   pointer to the FLUX data structure.
* NREG    number of regions for problem.
* NGROUP  number of groups for problem.
* NDIM    number of dimensions of problem.
* POLOAQ  polar quadrature order.
* AZMOAQ  azimuthal quadrature order.
* VOLUME  regional volumes.
* XPOL    polar quadrature points.
* WPOL    polar quadrature weights.
* XAZI    azimuthal quadrature points.
* WAZI    azimuthal quadrature weights.
*
*Parameters: output
* FLUX    direct and adjoint flux.
* AFLUX   angular components of the direct and adjoint flux.
* TFLUX   temporary flux vector.
* WGHT    temporary weight.
* MU      temporary mu.
* ETA     temporary eta.
*
*----------
*
      USE              GANLIB
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      TYPE(C_PTR)      KENTRY(NENTRY)
      INTEGER          NENTRY,IENTRY(NENTRY)
      INTEGER          IPRINT,NOPT,IOPT(NOPT),IKFLU
      INTEGER          NREG,NGROUP,NDIM,POLOAQ,AZMOAQ
      REAL             VOLUME(NREG),XPOL(POLOAQ,2),WPOL(POLOAQ),
     >                 XAZI(NDIM,AZMOAQ),WAZI(AZMOAQ)
      REAL             FLUX(NREG,2,NGROUP,2)
      DOUBLE PRECISION AFLUX(NREG,POLOAQ,AZMOAQ*2,2,NGROUP)
      DOUBLE PRECISION TFLUX(NREG+1,POLOAQ,AZMOAQ*2)
      REAL             WGHT(POLOAQ*AZMOAQ*2),MU(POLOAQ*AZMOAQ*2),
     >                 ETA(POLOAQ*AZMOAQ*2)
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='FMTSUD')
      INTEGER          ILCMUP,ILCMDN
      PARAMETER       (ILCMUP=1,ILCMDN=2)
*----
*  Local variables
*----
      TYPE(C_PTR)      IPU
      INTEGER          IFPU,IAPU,IFORM,IGROUP,IKU,IAKU,NAZ,IA,IP,IR,IQUA
      CHARACTER*12     NAMREC
      REAL             RELERR,ERRMAX,ETATMP,FNORM
*----
*  Initialize FLUX vectors
*----
      NAZ=AZMOAQ
      DO IP=1,POLOAQ
        XPOL(IP,2)=SQRT(1.-XPOL(IP,1)*XPOL(IP,1))
      ENDDO
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6000) NAMSBR,NGROUP,NDIM,2*AZMOAQ,POLOAQ,NREG
        WRITE(IOUT,6001)
        WRITE(IOUT,6011)
     >  (XAZI(1,IA),XAZI(2,IA),WAZI(IA),IA=1,AZMOAQ)
        WRITE(IOUT,6002)
        WRITE(IOUT,6011)
     >  (XPOL(IP,1),XPOL(IP,2),WPOL(IP),IP=1,POLOAQ)
      ENDIF
      CALL XDRSET(FLUX,NREG*2*NGROUP*2,0.0)
      CALL XDDSET(AFLUX,NREG*POLOAQ*AZMOAQ*2*2*NGROUP,0.0D0)
*----
*  Get information from FLUX data structure.
*  1. Flux
*  2. Angular flux
*  3. Adjoint
*  4. Angular adjoint
*----
      IPU=KENTRY(IKFLU)
      CALL LCMSIX(IPU,'FLUXDIRECT  ',ILCMUP)
      DO IGROUP=1,NGROUP
        WRITE(NAMREC,'(A4,I3,5X)') 'FLUX',IGROUP
        CALL LCMGET(IPU,NAMREC,FLUX(1,1,IGROUP,1))
      ENDDO
      CALL LCMSIX(IPU,'ANGULAR DIR ',ILCMUP)
      DO IGROUP=1,NGROUP
        WRITE(NAMREC,'(A5,I3,4X)') 'AFLUX',IGROUP
        CALL LCMGET(IPU,NAMREC,AFLUX(1,1,1,1,IGROUP))
      ENDDO
      CALL LCMSIX(IPU,'ANGULAR DIR ',ILCMDN)
      CALL LCMSIX(IPU,'FLUXDIRECT  ',ILCMDN)
      CALL LCMSIX(IPU,'FLUXADJOINT ',ILCMUP)
      DO IGROUP=1,NGROUP
        WRITE(NAMREC,'(A4,I3,5X)') 'FLUX',IGROUP
        CALL LCMGET(IPU,NAMREC,FLUX(1,2,IGROUP,1))
      ENDDO
      CALL LCMSIX(IPU,'ANGULAR ADJ ',ILCMUP)
      DO IGROUP=1,NGROUP
        WRITE(NAMREC,'(A5,I3,4X)') 'AFLUX',IGROUP
        CALL LCMGET(IPU,NAMREC,AFLUX(1,1,1,2,IGROUP))
      ENDDO
      CALL LCMSIX(IPU,'ANGULAR ADJ ',ILCMDN)
      CALL LCMSIX(IPU,'FLUXADJOINT ',ILCMDN)
*----
*  Create first SUS file
*  Volume, and tracking directions
*----
      IFPU=FILUNIT(KENTRY(1))
      WRITE(IFPU,1000) NREG
      WRITE(IFPU,1001) 2*AZMOAQ*POLOAQ
*----
*  -\Omega (\varphi+\pi,\pi-\theta)
*----
      IQUA=0
      DO IA=1,AZMOAQ
        ETATMP=XAZI(2,NAZ+1-IA)
        IF(ETATMP .EQ. 0.0) THEN
          ETATMP=1.0E-10
        ENDIF
        DO IP=1,POLOAQ
          IQUA=IQUA+1
          MU(IQUA)=-XPOL(IP,2)*XAZI(1,NAZ+1-IA)
          ETA(IQUA)=-XPOL(IP,2)*ETATMP
          WGHT(IQUA)=WPOL(IP)/WAZI(NAZ+1-IA)
        ENDDO
      ENDDO
*----
*  +\Omega (\varphi,\theta)
*----
      DO IA=1,AZMOAQ
        ETATMP=XAZI(2,IA)
        IF(ETATMP .EQ. 0.0) THEN
          ETATMP=1.0E-10
        ENDIF
        DO IP=1,POLOAQ
          IQUA=IQUA+1
          MU(IQUA)=XPOL(IP,2)*XAZI(1,IA)
          ETA(IQUA)=XPOL(IP,2)*ETATMP
          WGHT(IQUA)=WPOL(IP)/WAZI(IA)
        ENDDO
      ENDDO
      WRITE(IFPU,1010) 'wght',
     >                (WGHT(IQUA),IQUA=1,POLOAQ*AZMOAQ*2)
      WRITE(IFPU,1010) 'mu  ',
     >                (MU(IQUA),IQUA=1,POLOAQ*AZMOAQ*2)
      WRITE(IFPU,1010) 'eta ',
     >                (ETA(IQUA),IQUA=1,POLOAQ*AZMOAQ*2)
      WRITE(IFPU,1010) 'vol ',
     >                (VOLUME(IR),IR=1,NREG)
*----
*  Print angular flux
*----
      IFPU=FILUNIT(KENTRY(2))
      IKU=IENTRY(2)
      IFORM=IOPT(2)
      ERRMAX=0.0
      IF(IKU .EQ. 3) THEN
        DO IGROUP=1,NGROUP
          WRITE(IFPU) 'Group',IGROUP
        ENDDO
        DO IGROUP=NGROUP+1,NGROUP+6
          WRITE(IFPU) 'Comments'
        ENDDO
      ELSE
        DO IGROUP=1,NGROUP
          WRITE(IFPU,'(A8,4X,I8)') 'Groups =',IGROUP
        ENDDO
        DO IGROUP=NGROUP+1,NGROUP+6
          WRITE(IFPU,'(A8)') 'Comments'
        ENDDO
      ENDIF
      DO IGROUP=1,NGROUP
        DO IA=1,AZMOAQ*2
          DO IP=1,POLOAQ
            IF(IFORM .EQ. 1) THEN
              DO IR=1,NREG
                TFLUX(IR,IP,IA)=AFLUX(IR,IP,IA,1,IGROUP)
                IF(IA .LE. AZMOAQ) THEN
                  FLUX(IR,1,IGROUP,2)=FLUX(IR,1,IGROUP,2)
     >           +REAL(AFLUX(IR,IP,IA,1,IGROUP)*WPOL(IP)/WAZI(NAZ+1-IA))
                ELSE
                  FLUX(IR,1,IGROUP,2)=FLUX(IR,1,IGROUP,2)
     >           +REAL(AFLUX(IR,IP,IA,1,IGROUP)*WPOL(IP)/WAZI(IA-NAZ))
                ENDIF
              ENDDO
              IF(IPRINT .GE. 100) THEN
                WRITE(IOUT,6003) IP,IA,IGROUP
                WRITE(IOUT,6010) (TFLUX(IR,IP,IA),IR=1,NREG)
              ENDIF
            ELSE
              IR=1
              TFLUX(IR,IP,IA)=AFLUX(IR,IP,IA,1,IGROUP)
              DO IR=1,NREG
                TFLUX(IR+1,IP,IA)=
     >            2.0D0*AFLUX(IR,IP,IA,1,IGROUP)-TFLUX(IR,IP,IA)
                IF(IA .LE. AZMOAQ) THEN
                  FLUX(IR,1,IGROUP,2)=FLUX(IR,1,IGROUP,2)
     >           +REAL(AFLUX(IR,IP,IA,1,IGROUP)*WPOL(IP)/WAZI(NAZ+1-IA))
                ELSE
                  FLUX(IR,1,IGROUP,2)=FLUX(IR,1,IGROUP,2)
     >           +REAL(AFLUX(IR,IP,IA,1,IGROUP)*WPOL(IP)/WAZI(IA-NAZ))
                ENDIF
              ENDDO
              IF(IPRINT .GE. 100) THEN
                WRITE(IOUT,6003) IP,IA,IGROUP
                WRITE(IOUT,6010)
     >          ((TFLUX(IR,IP,IA)+TFLUX(IR+1,IP,IA))/2.0,IR=1,NREG)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        IF(IFORM .EQ. 1) THEN
          IF(IKU .EQ. 3) THEN
            WRITE(IFPU) (((TFLUX(IR,IP,IA),IR=1,NREG),
     >                   IP=1,POLOAQ),IA=1,AZMOAQ)
            WRITE(IFPU) (((TFLUX(IR,IP,IA),IR=1,NREG),
     >                   IP=1,POLOAQ),IA=AZMOAQ*2,AZMOAQ+1,-1)
          ELSE
            WRITE(IFPU,1002) (((TFLUX(IR,IP,IA),IR=1,NREG),
     >                   IP=1,POLOAQ),IA=1,AZMOAQ)
            WRITE(IFPU,1002) (((TFLUX(IR,IP,IA),IR=1,NREG),
     >                   IP=1,POLOAQ),IA=AZMOAQ*2,AZMOAQ+1,-1)
          ENDIF
        ELSE
          IF(IKU .EQ. 3) THEN
            WRITE(IFPU) (((TFLUX(IR,IP,IA),IR=1,NREG+1),
     >                   IP=1,POLOAQ),IA=1,AZMOAQ)
            WRITE(IFPU) (((TFLUX(IR,IP,IA),IR=1,NREG+1),
     >                   IP=1,POLOAQ),IA=AZMOAQ*2,AZMOAQ+1,-1)
          ELSE
            WRITE(IFPU,1002) (((TFLUX(IR,IP,IA),IR=1,NREG+1),
     >                   IP=1,POLOAQ),IA=1,AZMOAQ)
            WRITE(IFPU,1002) (((TFLUX(IR,IP,IA),IR=1,NREG+1),
     >                   IP=1,POLOAQ),IA=AZMOAQ*2,AZMOAQ+1,-1)
          ENDIF
        ENDIF
      ENDDO
      IF(IPRINT .GE. 10) THEN
        FNORM=FLUX(1,1,1,2)/FLUX(1,1,1,1)
        DO IGROUP=1,NGROUP
          WRITE(6,6030) IGROUP
          DO IR=1,NREG
            RELERR=100.0*(FLUX(IR,1,IGROUP,2)
     >            -FNORM*FLUX(IR,1,IGROUP,1))
     >            /FLUX(IR,1,IGROUP,2)
            ERRMAX=MAX(ERRMAX,ABS(RELERR))
            IF(IPRINT .GE. 20) THEN
              WRITE(6,6031) IR,FNORM*FLUX(IR,1,IGROUP,1),
     >                      FLUX(IR,1,IGROUP,2),RELERR
            ENDIF
          ENDDO
        ENDDO
        WRITE(6,6020) ERRMAX
      ENDIF
*----
*  Print angular adjoint
*----
      IAPU=FILUNIT(KENTRY(3))
      IAKU=IENTRY(3)
      ERRMAX=0.0
      IF(IAKU .EQ. 3) THEN
        DO IGROUP=NGROUP,1,-1
          WRITE(IAPU) 'Group',IGROUP
        ENDDO
        DO IGROUP=NGROUP+1,NGROUP+6
          WRITE(IAPU) 'Comments'
        ENDDO
      ELSE
        DO IGROUP=NGROUP,1,-1
          WRITE(IAPU,'(A8,4X,I8)') 'Groups =',IGROUP
        ENDDO
        DO IGROUP=NGROUP+1,NGROUP+6
          WRITE(IAPU,'(A8)') 'Comments'
        ENDDO
      ENDIF
      DO IGROUP=NGROUP,1,-1
        DO IA=1,AZMOAQ*2
          DO IP=1,POLOAQ
            IF(IFORM .EQ. 1) THEN
              DO IR=1,NREG
                TFLUX(IR,IP,IA)=AFLUX(IR,IP,IA,2,IGROUP)
                IF(IA .LE. AZMOAQ) THEN
                  FLUX(IR,2,IGROUP,2)=FLUX(IR,2,IGROUP,2)
     >           +REAL(AFLUX(IR,IP,IA,2,IGROUP)*WPOL(IP)/WAZI(NAZ+1-IA))
                ELSE
                  FLUX(IR,2,IGROUP,2)=FLUX(IR,2,IGROUP,2)
     >           +REAL(AFLUX(IR,IP,IA,2,IGROUP)*WPOL(IP)/WAZI(IA-NAZ))
                ENDIF
              ENDDO
              IF(IPRINT .GE. 100) THEN
                WRITE(IOUT,6004) IP,IA,IGROUP
                WRITE(IOUT,6010) (TFLUX(IR,IP,IA),IR=1,NREG)
              ENDIF
            ELSE
              IR=1
              TFLUX(IR,IP,IA)=AFLUX(IR,IP,IA,2,IGROUP)
              DO IR=1,NREG
                TFLUX(IR+1,IP,IA)=
     >            2.0D0*AFLUX(IR,IP,IA,2,IGROUP)-TFLUX(IR,IP,IA)
                IF(IA .LE. AZMOAQ) THEN
                  FLUX(IR,2,IGROUP,2)=FLUX(IR,2,IGROUP,2)
     >           +REAL(AFLUX(IR,IP,IA,2,IGROUP)*WPOL(IP)/WAZI(NAZ+1-IA))
                ELSE
                  FLUX(IR,2,IGROUP,2)=FLUX(IR,2,IGROUP,2)
     >           +REAL(AFLUX(IR,IP,IA,2,IGROUP)*WPOL(IP)/WAZI(IA-NAZ))
                ENDIF
              ENDDO
              IF(IPRINT .GE. 100) THEN
                WRITE(IOUT,6004) IP,IA,IGROUP
                WRITE(IOUT,6010)
     >          ((TFLUX(IR,IP,IA)+TFLUX(IR+1,IP,IA))/2.0,IR=1,NREG)
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        IF(IFORM .EQ. 1) THEN
          IF(IAKU .EQ. 3) THEN
            WRITE(IAPU) (((TFLUX(IR,IP,IA),IR=1,NREG),
     >                   IP=1,POLOAQ),IA=1,AZMOAQ)
            WRITE(IAPU) (((TFLUX(IR,IP,IA),IR=1,NREG),
     >                   IP=1,POLOAQ),IA=AZMOAQ*2,AZMOAQ+1,-1)
          ELSE
            WRITE(IAPU,1002) (((TFLUX(IR,IP,IA),IR=1,NREG),
     >                   IP=1,POLOAQ),IA=1,AZMOAQ)
            WRITE(IAPU,1002) (((TFLUX(IR,IP,IA),IR=1,NREG),
     >                   IP=1,POLOAQ),IA=AZMOAQ*2,AZMOAQ+1,-1)
          ENDIF
        ELSE
          IF(IAKU .EQ. 3) THEN
            WRITE(IAPU) (((TFLUX(IR,IP,IA),IR=1,NREG+1),
     >                   IP=1,POLOAQ),IA=1,AZMOAQ)
            WRITE(IAPU) (((TFLUX(IR,IP,IA),IR=1,NREG+1),
     >                   IP=1,POLOAQ),IA=AZMOAQ*2,AZMOAQ+1,-1)
          ELSE
            WRITE(IAPU,1002) (((TFLUX(IR,IP,IA),IR=1,NREG+1),
     >                   IP=1,POLOAQ),IA=1,AZMOAQ)
            WRITE(IAPU,1002) (((TFLUX(IR,IP,IA),IR=1,NREG+1),
     >                   IP=1,POLOAQ),IA=AZMOAQ*2,AZMOAQ+1,-1)
          ENDIF
        ENDIF
      ENDDO
      IF(IPRINT .GE. 10) THEN
        FNORM=FLUX(1,2,1,2)/FLUX(1,2,1,1)
        DO IGROUP=NGROUP,1,-1
          WRITE(6,6032) IGROUP
          DO IR=1,NREG
            RELERR=100.0*(FLUX(IR,2,IGROUP,2)
     >            -FNORM*FLUX(IR,2,IGROUP,1))
     >            /FLUX(IR,1,IGROUP,2)
            ERRMAX=MAX(ERRMAX,ABS(RELERR))
            IF(IPRINT .GE. 20) THEN
              WRITE(6,6031) IR,FNORM*FLUX(IR,2,IGROUP,1),
     >                      FLUX(IR,2,IGROUP,2),RELERR
            ENDIF
          ENDDO
        ENDDO
        WRITE(6,6020) ERRMAX
      ENDIF
*----
*  Processing finished, return
*----
      RETURN
*----
*  Formats
*----
 1000 FORMAT('Number of regions =',I5)
 1001 FORMAT('Nomber of angles  =',I5)
 1002 FORMAT(1P,5E20.10)
 1010 FORMAT(A4,1P/(3E15.7))
 6000 FORMAT('Output from routine ',A6/
     >       'Number of groups  =',I5/
     >       'Number of dimens  =',I5/
     >       'Number of azimuth =',I5/
     >       'Number of polar   =',I5/
     >       'Number of regions =',I5)
 6001 FORMAT('Azimuthal quadrature')
 6002 FORMAT('Polar quadrature')
 6003 FORMAT('Polar =',I5,2X,'Azim  = ',I5,2X,'Group = ',I5/
     >       'Direct angular flux per region ')
 6004 FORMAT('Polar =',I5,2X,'Azim  = ',I5,2X,'Group = ',I5/
     >       'Adjoint angular flux per region')
 6010 FORMAT(1P,5E20.10)
 6011 FORMAT(1P,3E20.10)
 6020 FORMAT('Maximum relative on flux    (%) = ',F15.7)
 6030 FORMAT(' Flux for group    ',I5)
 6031 FORMAT(' Region ',I5,1P,2E20.10,5X,0P,F20.10)
 6032 FORMAT(' Adjoint for group ',I5)
      END
