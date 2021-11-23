*DECK IDET01
      SUBROUTINE IDET01(IPTRK,IPFLU,IPLIB,IPMAP,IMPX,NDETC,MAXNI,NINX,
     > NINY,NINZ,COORD1,COORD2,COORD3,DETNAM,REANAM,DETECT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute detector integrated response on Cartesian geometry
*
*Copyright:
* Copyright (C) 2019 Ecole Polytechnique de Montreal
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* IPTRK   pointer to the tracking.
* IPFLU   pointer to the finite-element flux.
* IPLIB   pointer to the interpolated microlib.
* IPMAP   pointer to the fuelmap.
* IMPX    print parameter.
* NDETC   number of detectors
* MAXNI   first dimension of matrices NIN and COORD.
* NINX    number of interpolation points per detector along x axis.
* NINY    number of interpolation points per detector along y axis.
* NINZ    number of interpolation points per detector along z axis.
* COORD1  interpolation points per detector along x axis.
* COORD2  interpolation points per detector along y axis.
* COORD3  interpolation points per detector along z axis.
* DETNAM  character*12 alias name of the isotope used as detector.
* REANAM  character*12 name of the nuclear reaction used as detector.
*
*Parameters: output
* DETECT  detector response.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,IPFLU,IPLIB,IPMAP
      INTEGER IMPX,NDETC,MAXNI,NINX(MAXNI),NINY(MAXNI),NINZ(MAXNI)
      REAL COORD1(MAXNI,NDETC),COORD2(MAXNI,NDETC),COORD3(MAXNI,NDETC),
     > DETECT(NDETC)
      CHARACTER DETNAM*12,REANAM*12
*----
*  LOCAL VARIABLES
*----
      INTEGER NSTATE
      PARAMETER (NSTATE=40)
      INTEGER ISTATE(NSTATE)
      LOGICAL L3D
      CHARACTER HSMG*131
      TYPE(C_PTR) JPFLU,JPLIB,KPLIB
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, DIMENSION(:), ALLOCATABLE :: MAT,KFLX,KN,IMIX
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: MIXT
      REAL, DIMENSION(:), ALLOCATABLE :: XX,YY,ZZ,XXX,YYY,ZZZ,MXD,MYD,
     > MZD,FLXD,GAR,TERPX,TERPY,TERPZ
      REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: TFLUX,AFLUX,DFLUX,SIGF
      CHARACTER(LEN=12), DIMENSION(:), ALLOCATABLE :: HNAMIS
      TYPE(C_PTR), DIMENSION(:), ALLOCATABLE :: IPISO
*----
*  RECOVER GENERAL TRACKING INFORMATION
*----
      IF(.NOT.C_ASSOCIATED(IPTRK)) CALL XABORT('IDET01: IPTRK not set.')
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      NREG=ISTATE(1)
      NUN=ISTATE(2)
      ITYPE=ISTATE(6)
      NLF=0
      ICHX=0
      IDIM=1
      IF(ITYPE.EQ.5) THEN
        IDIM=2
      ELSE IF(ITYPE.EQ.7) THEN
        IDIM=3
      ELSE
        CALL XABORT('IDET01: Cartesian geometry expected.')
      ENDIF
      IELEM=ISTATE(9)
      L4=ISTATE(11)
      ICHX=ISTATE(12)
      NLF=ISTATE(30)
      NXD=ISTATE(14)
      NYD=ISTATE(15)
      NZD=ISTATE(16)
      L3D=(NZD.GT.0)
      IF(.NOT.L3D) CALL XABORT('IDET01: 3D geometry expected.')
      NZD=MAX(1,NZD)
      ALLOCATE(MAT(NREG),KFLX(NREG))
      CALL LCMGET(IPTRK,'MATCOD',MAT)
      CALL LCMGET(IPTRK,'KEYFLX',KFLX)
      ALLOCATE(MXD(NXD+1),MYD(NYD+1),MZD(NZD+1))
      ALLOCATE(XX(NREG),YY(NREG),ZZ(NREG))
      CALL LCMGET(IPTRK,'XX',XX)
      CALL LCMGET(IPTRK,'YY',YY)
      CALL LCMGET(IPTRK,'ZZ',ZZ)
*----
*  RECOVER FINITE-ELEMENT FLUX INFORMATION
*----
      IF(.NOT.C_ASSOCIATED(IPFLU)) CALL XABORT('IDET01: IPFLU not set.')
      CALL LCMGET(IPFLU,'STATE-VECTOR',ISTATE)
      NG=ISTATE(1)
*----
*  RECOVER RENUMBERED MIXTURE INDICES FROM FUELMAP
*----
      IF(C_ASSOCIATED(IPMAP)) THEN
        CALL LCMGET(IPMAP,'STATE-VECTOR',ISTATE)
        IF(ISTATE(4).NE.NG) CALL XABORT('IDET01: invalid group nb(1).')
        CALL LCMSIX(IPMAP,'GEOMAP',1)
          CALL LCMGET(IPMAP,'STATE-VECTOR',ISTATE)
          NEL=NXD*NYD*NZD
          IF(ISTATE(3).NE.NXD) CALL XABORT('IDET01: invalid NXD.')
          IF(ISTATE(4).NE.NYD) CALL XABORT('IDET01: invalid NYD.')
          IF(ISTATE(5).NE.NZD) CALL XABORT('IDET01: invalid NZD.')
          IF(ISTATE(6).NE.NEL) CALL XABORT('IDET01: invalid NEL.')
          IF(NREG.NE.NEL) CALL XABORT('IDET01: invalid NREG.')
        CALL LCMSIX(IPMAP,' ',2)
        CALL LCMGET(IPMAP,'BMIX',MAT)
      ENDIF
*----
*  RECOVER MICROLIB INFORMATION
*----
      IF(.NOT.C_ASSOCIATED(IPLIB)) CALL XABORT('IDET01: IPLIB not set.')
      CALL LCMGET(IPLIB,'STATE-VECTOR',ISTATE)
      MAXMIX=ISTATE(1)
      NBISO=ISTATE(2)
      IF(ISTATE(3).NE.NG) CALL XABORT('IDET01: invalid group nb(2).')
      NMIX=ISTATE(14)
      ALLOCATE(IMIX(NBISO),HNAMIS(NBISO),IPISO(NBISO),GAR(NG))
      CALL LCMGET(IPLIB,'ISOTOPESMIX',IMIX)
      CALL LCMGTC(IPLIB,'ISOTOPESUSED',12,NBISO,HNAMIS)
      DO ISO=1,NBISO
        IF(HNAMIS(ISO).EQ.DETNAM) GO TO 10
      ENDDO
      WRITE(HSMG,'(48HIDET01: NO DETECTOR ISOTOPE FOUND IN MICROLIB WI,
     > 8HTH NAME=,A12)') DETNAM
      CALL XABORT(HSMG)
   10 CALL LIBIPS(IPLIB,NBISO,IPISO)
      JPLIB=LCMGID(IPLIB,'ISOTOPESLIST')
*----
*  COMPUTE MESH FROM L_TRACK
*----
      ALLOCATE(XXX(NXD),YYY(NYD))
      CALL XDRSET(XXX,NXD,0.0)
      CALL XDRSET(YYY,NYD,0.0)
      IREG=0
      IF(L3D) THEN
        ALLOCATE(ZZZ(NZD))
        CALL XDRSET(ZZZ,NZD,0.0)
        DO K=1,NZD
          DO J=1,NYD
            DO I=1,NXD
              IREG=IREG+1
              IF(XX(IREG).NE.0.0) THEN
                IF(XXX(I).EQ.0.0) THEN
                  XXX(I)=XX(IREG)
                ELSE IF(ABS(XXX(I)-XX(IREG)).GT.1.0E-6) THEN
                  CALL XABORT('IDET01: inconsistent tracking in X')
                ENDIF
              ENDIF
              IF(YY(IREG).NE.0.0) THEN
                IF(YYY(J).EQ.0.0) THEN
                  YYY(J)=YY(IREG)
                ELSE IF(ABS(YYY(J)-YY(IREG)).GT.1.0E-6) THEN
                  CALL XABORT('IDET01: inconsistent tracking in Y')
                ENDIF
              ENDIF
              IF(ZZ(IREG).NE.0.0) THEN
                IF(ZZZ(K).EQ.0.0) THEN
                  ZZZ(K)=ZZ(IREG)
                ELSE IF(ABS(ZZZ(K)-ZZ(IREG)).GT.1.0E-6) THEN
                  CALL XABORT('IDET01: inconsistent tracking in Z')
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO J=1,NYD
          DO I=1,NXD
            IREG=IREG+1
            IF(XX(IREG).NE.0.0) THEN
              IF(XXX(I).EQ.0.0) THEN
                XXX(I)=XX(IREG)
              ELSE IF(ABS(XXX(I)-XX(IREG)).GT.1.0E-6) THEN
                CALL XABORT('IDET01: inconsistent tracking in X')
              ENDIF
            ENDIF
            IF(YY(IREG).NE.0.0) THEN
              IF(YYY(J).EQ.0.0) THEN
                YYY(J)=YY(IREG)
              ELSE IF(ABS(YYY(J)-YY(IREG)).GT.1.0E-6) THEN
                CALL XABORT('IDET01: inconsistent tracking in Y')
              ENDIF
            ENDIF
          ENDDO
        ENDDO
      ENDIF
      IF(IREG.NE.NREG) CALL XABORT('IDET01: invalid tracking')
      MXD(1)=0.0
      MYD(1)=0.0
      MZD(1)=0.0
      DO I=1,NXD
        MXD(I+1)=MXD(I)+XXX(I)
      ENDDO
      MYD(1)=0.0
      DO I=1,NYD
        MYD(I+1)=MYD(I)+YYY(I)
      ENDDO
      MZD(1)=0.0
      IF(L3D) THEN
        DO I=1,NZD
          MZD(I+1)=MZD(I)+ZZZ(I)
        ENDDO
        DEALLOCATE(ZZZ)
      ELSE
        MZD(2)=0.0
      ENDIF
      DEALLOCATE(YYY,XXX)
*----
*  PERFORM FLUX INTERPOLATION OVER DETECTOR LOCATIONS
*----
      IF(IMPX.GT.1) THEN
        WRITE(6,'(/29H IDET01: DETECTOR INFORMATION)')
        WRITE(6,'(5X,12HENERGY GROUP,1X,8HDETECTOR,2X,7HMIXTURE,5X,
     1  13HDETECTOR FLUX,3X,11HDONJON FLUX,5X,9HFLUX RATO,7X,
     2  11HDRAGON FLUX,5X,10HFISSION XS)')
      ENDIF
      ALLOCATE(FLXD(NUN))
      JPFLU=LCMGID(IPFLU,'FLUX')
      DO I=1,NDETC
        ININX=NINX(I)
        ININY=NINY(I)
        ININZ=NINZ(I)
        ALLOCATE(TFLUX(ININX,ININY,ININZ,NG))
        DO IG=1,NG
          CALL LCMGDL(JPFLU,IG,FLXD)
          IF(ICHX.EQ.1) THEN
*           Variational collocation method
            CALL LCMLEN(IPTRK,'KN',MAXKN,ITYLCM)
            MKN=MAXKN/(NXD*NYD*NZD)
            ALLOCATE(KN(MAXKN))
            CALL LCMGET(IPTRK,'KN',KN)
            CALL LCMSIX(IPTRK,'BIVCOL',1)
            CALL LCMLEN(IPTRK,'T',LC,ITYLCM)
            CALL LCMGET(IPTRK,'E',E)
            CALL LCMSIX(IPTRK,' ',2)
            CALL VALUE2(LC,MKN,NXD,NYD,NZD,L4,COORD1(1,I),COORD2(1,I),
     1      COORD3(1,I),MXD,MYD,MZD,FLXD,MAT,KN,ININX,ININY,ININZ,E,
     2      TFLUX(1,1,1,IG))
            DEALLOCATE(KN)
          ELSE IF(ICHX.EQ.2) THEN
*           Raviart-Thomas finite element method
            CALL VALUE4(IELEM,NUN,NXD,NYD,NZD,COORD1(1,I),COORD2(1,I),
     1      COORD3(1,I),MXD,MYD,MZD,FLXD,MAT,KFLX,ININX,ININY,ININZ,
     2      TFLUX(1,1,1,IG))
          ELSE IF(ICHX.EQ.3) THEN
*           Nodal collocation method (MCFD)
            CALL VALUE1(IDIM,NXD,NYD,NZD,L4,COORD1(1,I),COORD2(1,I),
     1      COORD3(1,I),MXD,MYD,MZD,FLXD,MAT,IELEM,ININX,ININY,ININZ,
     2      TFLUX(1,1,1,IG))
          ELSE
            CALL XABORT('IDET01: interpolation not implemented(1).')
          ENDIF
        ENDDO
*----
*  RECOVER AVERAGED FLUX FROM FINITE-ELEMENT CALCULATION
*----
        ALLOCATE(AFLUX(ININX,ININY,ININZ,NG),MIXT(ININX,ININY,ININZ))
        DO INX=1,ININX
          NX=0
          DO IX=1,NXD
            IF(COORD1(INX,I).LE.MXD(IX)) EXIT
            NX=NX+1
          ENDDO
          DO INY=1,ININY
            NY=0
            DO IY=1,NYD
              IF(COORD2(INY,I).LE.MXD(IY)) EXIT
              NY=NY+1
            ENDDO
            DO INZ=1,ININZ
              NZ=0
              DO IZ=1,NZD
                IF(COORD3(INZ,I).LE.MZD(IZ)) EXIT
                NZ=NZ+1
              ENDDO
              IF(NX*NY*NZ.EQ.0) THEN
                WRITE(HSMG,'(38HIDET01: element not found for detector,
     1          I5,7h(1). x=,1p,e12.4,3H y=,e12.4,3H z=,e12.4)') I,
     2          COORD1(INX,I),COORD2(INY,I),COORD3(INZ,I)
                CALL XABORT(HSMG)
              ENDIF
              IEL=(NZ-1)*NXD*NYD+(NY-1)*NXD+NX
              IF(MAT(IEL).EQ.0) THEN
                WRITE(HSMG,'(38HIDET01: element not found for detector,
     1          I5,7h(1). x=,1p,e12.4,3H y=,e12.4,3H z=,e12.4)') I,
     2          COORD1(INX,I),COORD2(INY,I),COORD3(INZ,I)
                CALL XABORT(HSMG)
              ENDIF
              MIXT(INX,INY,INZ)=MAT(IEL)
              IUN=KFLX(IEL)
              IF(IUN.EQ.0) CALL XABORT('IDET01: flux not defined.')
              DO IG=1,NG
                CALL LCMGDL(JPFLU,IG,FLXD)
                AFLUX(INX,INY,INZ,IG)=FLXD(IUN)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
*----
*  RECOVER FLUX AND FISSION CROSS SECTION FROM MICROLIB
*----
        ALLOCATE(DFLUX(ININX,ININY,ININZ,NG),SIGF(ININX,ININY,ININZ,NG))
        DO INX=1,ININX
          DO INY=1,ININY
            DO INZ=1,ININZ
              IBM=MIXT(INX,INY,INZ)
              DFLUX(INX,INY,INZ,:NG)=0.0
              SIGF(INX,INY,INZ,:NG)=0.0
              DO ISO=1,NBISO
                IF((HNAMIS(ISO).EQ.DETNAM).AND.(IMIX(ISO).EQ.IBM)) THEN
                  KPLIB=IPISO(ISO) ! set ISO-th isotope
                  CALL LCMLEN(KPLIB,REANAM,LENGT,ITYLCM)
                  IF(LENGT.NE.NG) THEN
                    CALL LCMLIB(KPLIB)
                    WRITE(HSMG,'(23HIDET01: unable to find ,A,6H for i,
     >              7Hsotope ,A,11H in mixture,I6)') REANAM,DETNAM,IBM
                    CALL XABORT(HSMG)
                  ENDIF
                  CALL LCMGET(KPLIB,'NWT0',GAR)
                  DFLUX(INX,INY,INZ,:NG)=GAR(:NG)
                  CALL LCMGET(KPLIB,REANAM,GAR)
                  SIGF(INX,INY,INZ,:NG)=GAR(:NG)
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
*----
*  PRINT DETECTOR-DEPENDENT VALUES
*----
        IF(IMPX.GT.1) THEN
          DO IG=1,NG
            DO INZ=1,ININZ
              DO INY=1,ININY
                DO INX=1,ININX
                  IBM=MIXT(INX,INY,INZ)
                  TFLUX_I=TFLUX(INX,INY,INZ,IG)
                  AFLUX_I=AFLUX(INX,INY,INZ,IG)
                  DFLUX_I=DFLUX(INX,INY,INZ,IG)
                  SIGF_I=SIGF(INX,INY,INZ,IG)
                  WRITE(6,'(8X,3I9,1P,5E16.5)') IG,I,IBM,TFLUX_I,
     >            AFLUX_I,TFLUX_I/AFLUX_I,DFLUX_I,SIGF_I
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDIF
*----
*  COMPUTE DETECTOR RESPONSE
*----
        ALLOCATE(TERPX(ININX),TERPY(ININY),TERPZ(ININZ))
        IF(ININX.EQ.1) THEN
          TERPX(1)=1.0
        ELSE
          CALL ALTERI(.TRUE.,ININX,COORD1(1,I),COORD1(1,I),
     >    COORD1(ININX,I),TERPX)
        ENDIF
        IF(ININY.EQ.1) THEN
          TERPY(1)=1.0
        ELSE
          CALL ALTERI(.TRUE.,ININY,COORD2(1,I),COORD2(1,I),
     >    COORD2(ININY,I),TERPY)
        ENDIF
        IF(ININZ.EQ.1) THEN
          TERPZ(1)=1.0
        ELSE
          CALL ALTERI(.TRUE.,ININZ,COORD3(1,I),COORD3(1,I),
     >    COORD3(ININZ,I),TERPZ)
        ENDIF
*       integrate along axial direction
        DETECT(I)=0.0
        DO IG=1,NG
          ZNUM=0.0
          ZDEN=0.0
          DO INX=1,ININX
            DO INY=1,ININY
              DO INZ=1,ININZ
                TRP=TERPX(INX)*TERPY(INY)*TERPZ(INZ)
                ZNUM=ZNUM+TRP*TFLUX(INX,INY,INZ,IG)*SIGF(INX,INY,INZ,IG)
                ZDEN=ZDEN+TRP
              ENDDO
            ENDDO
          ENDDO
          DETECT(I)=DETECT(I)+ZNUM/ZDEN
        ENDDO
        DEALLOCATE(TERPZ,TERPY,TERPX)
        DEALLOCATE(SIGF,DFLUX,MIXT,AFLUX,TFLUX)
      ENDDO
*----
*  RELEASE MEMORY
*----
      DEALLOCATE(FLXD)
      DEALLOCATE(GAR,IPISO,HNAMIS,IMIX)
      DEALLOCATE(XX,YY,ZZ,KFLX,MAT)
      RETURN
      END
