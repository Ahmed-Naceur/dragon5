*DECK SENDRV
      SUBROUTINE SENDRV(IPSENS,IPTRACK,IPLIB,IPFLUX,IPAFLUX,IPRINT,
     >                  NR,NU,NI,NG,NANIS,NLTERM)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Creation of sensitivity profiles to
* cross-section on the reactivity using first order perturbation
* method using the adjoint calculation.
*
*Copyright:
* Copyright (C) 2011 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): C. Laville, G. Marleau
*
*Parameters: input/output
* IPSENS  LCM Sensitivity object address.
* IPTRACK LCM Tracking object address.
* IPLIB   LCM Library object address.
* IPFLUX  LCM Flux object address.
* IPAFLUX LCM Adjoin flux object address.
* IPRINT  print level.
* NR      number of region in Tracking object.
* NU      number of unkwnow in Tracking/Flux objects.
* NI      number of isotopes in Library object.
* NG      number of energy group in Library object.
* NANIS   anisotropy order kept.
* NLTERM  total number of term of the flux.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT      NONE
*----
*  Suboutine arguements
*----
      TYPE(C_PTR)  IPTRACK,IPLIB,IPFLUX,IPAFLUX
      INTEGER      IPSENS,NR,NU,NI,NG,IPRINT,NANIS,NLTERM
*----
*  Local variables
*----
      TYPE(C_PTR)  KPISO,JPFLUX,JPAFLUX
      INTEGER      IG,JG,IR,KL,KKEYL,IL,IM,II
      INTEGER      ILENG,ILON,ITYLCM,NSENI,NSENS,NIC
      REAL         KEFF,D
      CHARACTER    ISONAM*12,REV*48,DATE*64,HSMG*131
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MAT,KEY,MELISO,ISOC
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NAMISO,NAMISC
      REAL, ALLOCATABLE, DIMENSION(:) :: VOL,ENERG,DENISO,SOLFLU,NUSIGF,
     > CHI
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: FLUXI,AFLUX
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: P
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: IPISO
*----
*  Parameters
*----
      INTEGER      IOUT
      CHARACTER    NAMSBR*6
      PARAMETER    (IOUT=6,NAMSBR='SENDRV')
*----
*  Scratch storage allocation
*   MAT      Mixture of each region in Tracking object.
*   KEY      Localisation of the flux/aflux for each region in Tracking
*            object.
*   MELISO   Mixture of each isotope/region in Library object.
*   NAMISO   Alias of each isotope/region in Library object.
*   NAMISC   Independent isotopes names.
*   ISOC     Independent isotope number associated with isotope/mixture.
*   VOL      Volume of each region in Tracking object.
*   ENERG    Energy boundaries in Library object.
*   DENISO   Density of each isotope/region in Library object.
*   SOLFLU   Flux solution of each region in Flux object.
*   FLUXI    Flux of each region multiplied by the volume of the region.
*   AFLUX    Adjoint flux of each region.
*   P        Matrix calcution for sensitivity analysis.
*   NUSIGF   Vector for nubar multiplied by the fission cross-section
*   CHI      Vector for fission spectra.
*----
      ALLOCATE(MAT(NR),KEY(NU),MELISO(NI),NAMISO(3,NI),NAMISC(2,NI),
     > ISOC(NI))
      ALLOCATE(VOL(NR),ENERG(NG+1),DENISO(NI),SOLFLU(NU),
     > FLUXI(NR,NLTERM,NG),AFLUX(NR,NLTERM,NG),P(NR,NANIS,NG,NG),
     > NUSIGF(NG),CHI(NG))
*--------
* Step 1:
* Recover Tracking informations
* a) VOLUME
* b) MATCOL
* c) KEYFLX
*----
      CALL XDRSET(VOL,NR,0.0)
      CALL LCMGET(IPTRACK,'VOLUME',VOL)
*----
      CALL XDISET(MAT,NR,0)
      CALL LCMGET(IPTRACK,'MATCOD',MAT)
*----
      CALL XDISET(KEY,NU,0)
      CALL LCMGET(IPTRACK,'KEYFLX',KEY)
*----
*  Print information if required
*----
      IF(IPRINT .GE. 5) THEN
        WRITE(IOUT,'(3(A12,8X))') 'Region #    ','VOLUMES     ',
     >                            'MATCOD      '
        WRITE(IOUT,'(I10,10X,E16.7,4X,I8)')
     >  (IR,VOL(IR),MAT(IR),IR=1,NR)
      ENDIF
*--------
* Step 2:
* Recover library informations
* a) ENERGY GROUP BOUNDARIES
* b) MIXTURE ASSOCIATED WITH ISOTOPE
* c) ISOTOPES ALIAS
* d) ISOTOPES DENSITIES
*----
      CALL XDRSET(ENERG,NG+1,0.0)
      CALL LCMGET(IPLIB,'ENERGY',ENERG)
*----
      CALL XDISET(MELISO,NI,0)
      CALL LCMGET(IPLIB,'ISOTOPESMIX',MELISO)
*----
      CALL LCMGET(IPLIB,'ISOTOPESUSED',NAMISO)
*----
      CALL XDRSET(DENISO,NI,0.0)
      CALL LCMGET(IPLIB,'ISOTOPESDENS',DENISO)
*----
*  Scan library to determine number of sensitivity profiles:
*  NIC    is number of independent isotopes
*  ISOC   is the independent isotope number associated with
*         each isotope/mixture
*  NAMISC is the independent isotope name
*----
      CALL SENCNT(IPLIB,NI,NAMISO,MELISO,
     >            NSENS,NSENI,NAMISC,ISOC,NIC)
      IF(IPRINT.GE.1) THEN
        WRITE(IOUT,6000) NSENS+NSENI,NSENI
      ENDIF
*----
*  Print information if required
*----
      IF(IPRINT .GE. 5) THEN
        WRITE(IOUT,'(A12)') 'ENERGY      '
        WRITE(IOUT,'(5E14.6)') (ENERG(IG),IG=1,NG+1)
        WRITE(IOUT,'(6(A12,8X))') 'Isotope #   ','ISOTOPESMIX ',
     >                            'ISOTOPESUSED','ISOTOPESDENS',
     >                            'Name Int Sen','IsoComb     '
        WRITE(IOUT,'(I8,12X,I8,12X,3A4,8X,E16.7,4X,2A4,12X,I8)')
     >  (II,MELISO(II),NAMISO(1,II),NAMISO(2,II),NAMISO(3,II),
     >   DENISO(II),NAMISC(1,ABS(ISOC(II))),NAMISC(2,ABS(ISOC(II))),
     >   ISOC(II),II=1,NI)
      ENDIF
*--------
* Step 3:
* Recover K-effective, flux and adjoint,
* and compute integrated flux.
*----
      KEFF=0
      CALL LCMGET(IPFLUX,'K-EFFECTIVE',KEFF)
*----
      CALL XDRSET(FLUXI,NR*NLTERM*NG,0.0)
      CALL XDRSET(AFLUX,NR*NLTERM*NG,0.0)
      CALL LCMLEN(IPFLUX,'FLUX',ILON,ITYLCM)
      IF(ILON.EQ.0) CALL XABORT(NAMSBR//': MISSING FLUX INFO.')
      JPFLUX=LCMGID(IPFLUX,'FLUX')
      DO IG=1,NG
        CALL LCMGDL(JPFLUX,IG,SOLFLU)
        KL=0
        KKEYL=0
        DO IL=1,NANIS
          DO IM=1,IL  ! 2*IL-1 in 3D
            KL=KL+1
            DO IR=1,NR
              FLUXI(IR,KL,IG)=SOLFLU(KEY(IR+KKEYL))*VOL(IR)
            ENDDO
            KKEYL=KKEYL+NR
          ENDDO
        ENDDO
      ENDDO
      CALL LCMLEN(IPAFLUX,'AFLUX',ILON,ITYLCM)
      IF(ILON.EQ.0) CALL XABORT(NAMSBR//': MISSING AFLUX INFO.')
      JPAFLUX=LCMGID(IPAFLUX,'AFLUX')
      DO IG=1,NG
        CALL LCMGDL(JPAFLUX,IG,SOLFLU)
        KL=0
        KKEYL=0
        DO IL=1,NANIS
          DO IM=1,IL  ! 2*IL-1 in 3D
            KL=KL+1
            DO IR=1,NR
              AFLUX(IR,KL,IG)=SOLFLU(KEY(IR+KKEYL))
            ENDDO
            KKEYL=KKEYL+NR
          ENDDO
        ENDDO
      ENDDO
*--------
* Step 4:
* Computation of the P matrix used for sensitivity calculation
*----
      CALL XDRSET(P,NR*NANIS*NG*NG,0.0)
      DO IG=1,NG
        DO JG=1,NG
          KL=0
          DO IL=1,NANIS
            DO IM=1,IL  ! 2*IL-1 en 3D
              KL=KL+1
              DO IR=1,NR
                P(IR,IL,IG,JG)=P(IR,IL,IG,JG)+
     >                         FLUXI(IR,KL,IG)*AFLUX(IR,KL,JG)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
*--------
* Step 5:
* Computation of the D factor
*----
      D=0
      ALLOCATE(IPISO(NI))
      CALL LIBIPS(IPLIB,NI,IPISO)
      DO II=1,NI
        KPISO=IPISO(II) ! set II-th isotope
        IF(.NOT.C_ASSOCIATED(KPISO)) THEN
          WRITE(ISONAM,'(3A4)') NAMISO(1,II),NAMISO(2,II),NAMISO(3,II)
          WRITE(HSMG,'(17HSENDRV: ISOTOPE '',A12,7H'' (ISO=,I8,5H) IS ,
     1    30HNOT AVAILABLE IN THE MICROLIB.)') ISONAM,II
          CALL XABORT(HSMG)
        ENDIF
        CALL LCMLEN(KPISO,'NUSIGF',ILENG,ITYLCM)
        IF(ILENG.GT.0) THEN
          CALL XDRSET(NUSIGF,NG,0.0)
          CALL XDRSET(CHI,NG,0.0)
          CALL LCMGET(KPISO,'NUSIGF',NUSIGF)
          CALL LCMGET(KPISO,'CHI',CHI)
          DO IR=1,NR
            IF(MAT(IR).EQ.MELISO(II)) THEN
              DO IG=1,NG
                DO JG=1,NG
                  D=D+FLUXI(IR,1,IG)*DENISO(II)
     >               *NUSIGF(IG)*CHI(JG)*AFLUX(IR,1,JG)
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ENDDO
      DEALLOCATE(IPISO)
*--------
* Step 6:
* Initialize the top of the .sdf file
*----
      WRITE(IPSENS,7000) NG,NSENS+NSENI,NSENI,KEFF
      WRITE(IPSENS,7001)
      WRITE(IPSENS,7002) (ENERG(IG),IG=1,NG+1)
*--------
* Step 7:
* Sensitivity analysis for each isotope/mixture
*----
      CALL SENCAL(IPSENS,IPLIB,IPRINT,NR,NG,NI,NANIS,NAMISO,MELISO,
     >            MAT,DENISO,KEFF,P,D,NAMISC,ISOC,NIC)
*--------
* Step 9:
* Finish the creation of the .sdf file
*----
      CALL KDRVER(REV,DATE)
      WRITE(IPSENS,7010) REV,DATE
*----
*  Scratch storage deallocation
*----
      DEALLOCATE(CHI,NUSIGF,P,AFLUX,FLUXI,SOLFLU,DENISO,ENERG,VOL)
      DEALLOCATE(ISOC,NAMISC,NAMISO,MELISO,KEY,MAT)
      RETURN
*----
*  Formats
*----
 6000 FORMAT('Total number of profils (reactions,isotopes,mixtures) :',
     >5X,I10/'Total number of integrated profiles ',
     >'(reactions,isotopes) :',3X,I10)
 7000 FORMAT(
     >I10,3X,'number of neutron groups'/
     >I10,3X,'number of sensitivity profiles  ',
     >I10,3X,'are region integrated'/
     >F10.6,3X,'k-eff from the forward case')
 7001 FORMAT('energy boundaries:')
 7002 FORMAT(1P,5E14.6)
 7010 FORMAT(//1X,'file verification information'/
     >1X,'code system:    DRAGON'/
     >1X,'Version:        ',A48/
     >1X,'program:        SENS:'/
     >3X,'creation date:  ',A64/
     >3X,'library:        UNKWNOW'/
     >3X,'this is not a scale configuration controlled code')
      END
