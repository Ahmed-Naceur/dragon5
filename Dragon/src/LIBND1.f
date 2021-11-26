*DECK LIBND1
      SUBROUTINE LIBND1 (IPLIB,NAMFIL,NGRO,NBISO,NL,ISONAM,ISONRF,
     1 IPISO,MASKI,TN,SN,SB,IMPX,NGF,NGFR,NDEL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Transcription of the useful interpolated microscopic cross section
* data from NDAS to LCM data structures. Memory allocation interface.
*
*Copyright:
* Copyright (C) 2006 Ecole Polytechnique de Montreal
*
*Author(s): A. Hebert
*
*Parameters: input
* IPLIB   pointer to the lattice microscopic cross section library
*         (L_LIBRARY signature).
* NAMFIL  name of the NDAS file.
* NGRO    number of energy groups.
* NBISO   number of isotopes present in the calculation domain.
* NL      number of Legendre orders required in the calculation
*         NL=1 or higher.
* ISONAM  alias name of isotopes.
* ISONRF  library reference name of isotopes.
* IPISO   pointer array towards microlib isotopes.
* MASKI   isotopic mask. Isotope with index I is processed if
*         MASKI(I)=.true.
* TN      temperature of each isotope.
* SN      dilution cross section in each energy group of each
*         isotope. A value of 1.0E10 is used for infinite dilution.
* SB      dilution cross section as used by Livolant and Jeanpierre
*         normalization.
* IMPX    print flag
*
*Parameters: output
* NGF     number of fast groups without self-shielding.
* NGFR    number of fast and resonance groups.
* NDEL    number of precursor groups for delayed neutrons.
*
*Reference:
* P. J. Laughton, "NJOYPREP and WILMAPREP: UNIX-Based Tools for WIMS-
* AECL Cross-Section Library Production," Atomic Energy of Canada,
* Report COG-92-414 (Rev. 0), June 1993.
* Copyright (C) from NDAS Atomic Energy of Canada Limited utility (2006)
*
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE FSDF
      IMPLICIT NONE
*----
*  Subroutine arguments
*----
      TYPE(C_PTR) IPLIB,IPISO(NBISO)
      INTEGER NGRO,NBISO,NL,ISONAM(3,NBISO),ISONRF(3,NBISO),IMPX,NGF,
     1 NGFR,NDEL
      REAL TN(NBISO),SN(NGRO,NBISO),SB(NGRO,NBISO)
      CHARACTER NAMFIL*(*)
      LOGICAL MASKI(NBISO)
*----
*  Local variables
*----
      INTEGER IOUT
      PARAMETER(IOUT=6)
      TYPE(C_PTR) KPLIB
      INTEGER I,I0,IERR,HEADER(16),NISOLB,NGFIS,NGTHER,MAXTMP,MAXDIL,
     1 MAXTDN,MAXPN,NF,NP1,IND,IHEAD(200),NBTEM,NBDIL,ISOID,IG,IG1,NL2,
     2 IJ,IM,IMX,IOF,J,ITYPRO(2)
      REAL RHEAD(200),WW,SUM
      DOUBLE PRECISION XDRCST,ANEUT
      CHARACTER TEXT8*8,HSMG*131,HNAMIS*12,HNISOR*12
      LOGICAL LCUBIC
      PARAMETER(LCUBIC=.TRUE.)
      EXTERNAL XDRCST
      EQUIVALENCE(RHEAD(1),IHEAD(1))
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: HNAM
      REAL, ALLOCATABLE, DIMENSION(:) :: DELTA,TEMPS,DILUS,TERPT,LOAD,
     1 ENER,CHI,WT0,GC,RESD
      REAL, ALLOCATABLE, DIMENSION(:,:) :: GAR1,GAR2,THERXS,XA,XS,XF,XN
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SCAT
*----
*  Read NDAS library parameters
*----
      CALL XSDOPN(NAMFIL,IERR)
      IF(IERR.NE.0) CALL XABORT('LIBND1: XSDOPN could not open Library'
     >  //' files')
      CALL XSDBLD(6001,HEADER,IERR)
      IF(IERR.NE.0) CALL XABORT('LIBND1: XSDBLD could not read library'
     > //' parameters')
      IF(NGRO.NE.HEADER(2)) CALL XABORT('LIBND1: Invalid number of ene'
     > //'rgy groups')
      NISOLB=HEADER(1)
      NGFIS=HEADER(3)
      NGF=HEADER(4)
      NGFR=HEADER(4)+HEADER(5)
      NGTHER=HEADER(6)
      MAXTMP=HEADER(11)
      MAXDIL=HEADER(12)
      MAXTDN=HEADER(13)
      IF(HEADER(14).NE.2) CALL XABORT('LIBND1: Old NDAS format not sup'
     > //'ported')
      MAXPN=MAX(HEADER(15),HEADER(16))
      NDEL=0
      IF(IMPX.GT.1) WRITE(IOUT,100) (HEADER(I),I=1,16)
*----
*  Scratch storage allocation
*   HNAM    isotope names in NDAS library
*   DELTA   lethargy widths
*   TEMPS   temperature base points
*   DILUS   dilution base points
*   TERPT   interpolation factors in temperature
*   GAR1    vector xs components (1: transport corr. total;
*           2: absorption; 3: fission; 4: nu*fission; 5: P0 scattering;
*           6: P1 scattering, 7: (n,2n)
*   GAR2    self-shielded xs returned by LIBND3
*   SCAT    P0 and P1 differential scattering xs components
*   THERXS  temperature-dependent thermal cross section components
*   LOAD    storage array containing differential scattering components
*   XA      dilution-dependent absorption effective cross sections
*   XS      dilution-dependent scattering effective cross sections
*   XF      dilution-dependent nu*fission effective cross sections
*   XN      dilution-dependent NJOY fluxes
*----
      ALLOCATE(HNAM(2,NISOLB))
      ALLOCATE(DELTA(NGRO),TEMPS(MAXTMP),DILUS(MAXDIL),TERPT(MAXTMP),
     1 GAR1(NGRO,7),GAR2(NGRO,6),SCAT(NGRO,NGRO,2),THERXS(NGTHER,2),
     2 LOAD(MAXPN),XA(NGFR-NGF,MAXDIL),XS(NGFR-NGF,MAXDIL),
     3 XF(NGFR-NGF,MAXDIL),XN(NGFR-NGF,MAXDIL))
*----
*  Recover the group structure
*----
      ANEUT=XDRCST('Neutron mass','amu')
      ALLOCATE(ENER(NGRO+1))
      CALL XSDBLD(5019,ENER,IERR)
      IF(IERR.NE.0) CALL XABORT('LIBND1: xsdbld could not read energy '
     > //'group limits')
      IF(ENER(NGRO+1).EQ.0.0) ENER(NGRO+1)=1.0E-5
      DO I=1,NGRO
        DELTA(I)=LOG(ENER(I)/ENER(I+1))
      ENDDO
      CALL LCMPUT(IPLIB,'ENERGY',NGRO+1,2,ENER)
      CALL LCMPUT(IPLIB,'DELTAU',NGRO,2,DELTA)
      DEALLOCATE(ENER)
*----
*  Recover the isotope names and identifiers from the library
*----
      DO I=1,NISOLB
        CALL XSDNAM(I,ISOID,TEXT8,IERR)
        IF(IERR.NE.0) CALL XABORT('LIBND1: XSDNAM index overflow')
        READ(TEXT8,'(2A4)') HNAM(1,I),HNAM(2,I)
      ENDDO
*----
*  Read through NDAS file and accumulate cross sections for this range
*  of MATS, Legendre orders, and groups
*----
      DO IMX=1,NBISO
        IF(MASKI(IMX)) THEN
          WRITE(HNAMIS,'(3A4)') (ISONAM(I0,IMX),I0=1,3)
          WRITE(HNISOR,'(3A4)') (ISONRF(I0,IMX),I0=1,3)
          IND=0
          DO I=1,NISOLB
            IF((ISONRF(1,IMX).EQ.HNAM(1,I)).AND.
     >         (ISONRF(2,IMX).EQ.HNAM(2,I))) THEN
              IND=I
              GO TO 10
            ENDIF
          ENDDO
          WRITE (HSMG,130) HNAMIS,HNISOR,NAMFIL
          CALL XABORT(HSMG)
   10     IF(IMPX.GT.9) CALL LCMLIB(IPLIB)
          KPLIB=IPISO(IMX) ! set IMX-th isotope
*         Load nuclide header
          CALL XSDISO(7000,6001,IND,RHEAD,IERR)
          IF(IMPX.GT.0) THEN
            WRITE(IOUT,110) HNAMIS,HNISOR
          ENDIF
          IF(IMPX.GT.5) THEN
            WRITE(IOUT,120) IHEAD(1),IHEAD(2),RHEAD(3),(IHEAD(I),I=4,12)
          ENDIF
          CALL LCMPTC(KPLIB,'ALIAS',12,1,HNAMIS)
          CALL LCMPUT(KPLIB,'AWR',1,2,RHEAD(3)/REAL(ANEUT))
          NF=IHEAD(5)
          NBTEM=IHEAD(6)
          NP1=IHEAD(10)
          IF(NBTEM.GT.MAXTMP) CALL XABORT('LIBND1: MAXTMP overflow(1)')
          IF(NBTEM.EQ.1) THEN
            TERPT(1)=1.0
          ELSE
*           Thermal temperatures
            CALL XSDISO(7000,5017,IND,TEMPS,IERR)
            CALL ALTERP(LCUBIC,NBTEM,TEMPS,TN(IMX),.FALSE.,TERPT)
          ENDIF
*
*         transport-corrected total and absorption XS
          CALL XSDISO(7002,5013,IND,GAR1(:,1),IERR)
          CALL XSDISO(7002,5004,IND,GAR1(:,2),IERR)
          IF(NGTHER.GT.0) THEN
            CALL XDRSET(GAR1(NGFR+1,1),NGTHER,0.0)
            CALL XDRSET(GAR1(NGFR+1,2),NGTHER,0.0)
            DO I=1,NBTEM
              WW=TERPT(I)
              IF(ABS(WW).GT.1.0E-6) THEN
                CALL XSDTHE(7004,5013,-1,I,THERXS(:,1),IERR)
                CALL XSDTHE(7004,5004,-1,I,THERXS(:,2),IERR)
                DO I0=1,NGTHER
                  IOF=NGFR+I0
                  GAR1(IOF,1)=GAR1(IOF,1)+WW*THERXS(I0,1)
                  GAR1(IOF,2)=GAR1(IOF,2)+WW*THERXS(I0,2)
                ENDDO
              ENDIF
            ENDDO
          ENDIF
*
*         fission spectrum, fission XS and nu*fission XS
          IF(IHEAD(11).EQ.1) THEN
            ALLOCATE(CHI(NGRO))
            CALL XDRSET(CHI,NGRO,0.0)
            CALL XSDISO(7002,5008,IND,CHI,IERR)
            SUM=0.0
            DO I=1,NGRO
              SUM=SUM+CHI(I)
            ENDDO
            IF(ABS(SUM-1.0).GT.1.0E-5) CALL XABORT('LIBND1: Fission sp'
     >      //'ectrum does not sum to one')
            CALL LCMPUT(KPLIB,'CHI',NGRO,2,CHI)
            DEALLOCATE(CHI)
*
            CALL XSDISO(7002,5005,IND,GAR1(:,3),IERR)
            CALL XSDISO(7002,5006,IND,GAR1(:,4),IERR)
            IF(NGTHER.GT.0) THEN
              CALL XDRSET(GAR1(NGFR+1,3),NGTHER,0.0)
              CALL XDRSET(GAR1(NGFR+1,4),NGTHER,0.0)
              DO I=1,NBTEM
                WW=TERPT(I)
                IF(ABS(WW).GT.1.0E-6) THEN
                  CALL XSDTHE(7004,5005,-1,I,THERXS(:,1),IERR)
                  CALL XSDTHE(7004,5006,-1,I,THERXS(:,2),IERR)
                  DO I0=1,NGTHER
                    IOF=NGFR+I0
                    GAR1(IOF,3)=GAR1(IOF,3)+WW*THERXS(I0,1)
                    GAR1(IOF,4)=GAR1(IOF,4)+WW*THERXS(I0,2)
                  ENDDO
                ENDIF
              ENDDO
            ENDIF
          ELSE
            CALL XDRSET(GAR1(1,3),NGRO,0.0)
            CALL XDRSET(GAR1(1,4),NGRO,0.0)
          ENDIF
*
*         (n,2n) XS
          CALL XSDISO(7001,5007,IND,GAR1(:,7),IERR)
          CALL XDRSET(GAR1(NGF+1,7),NGRO-NGF,0.0)
          CALL LCMPUT(KPLIB,'N2N',NGRO,2,GAR1(1,7))
*
*         P0 differential scattering XS
          CALL XSDISO(7002,5015,IND,LOAD,IERR)
          CALL XDRSET(GAR1(1,5),NGRO,0.0)
          CALL XDRSET(SCAT(1,1,1),NGRO*NGRO,0.0)
          IJ=0
          DO IG=1,NGFR
            IM=NINT(LOAD(IJ+2))
            IG1=-NINT(LOAD(IJ+1))+IG
            IJ=IJ+2
            DO I0=1,IM
*             --- IG is the primary group
              SCAT(IG1+I0,IG,1)=LOAD(IJ+I0)
              GAR1(IG,5)=GAR1(IG,5)+LOAD(IJ+I0)
            ENDDO
            IJ=IJ+IM
          ENDDO
          IF(NGTHER.GT.0) THEN
            CALL XDRSET(GAR1(NGFR+1,5),NGTHER,0.0)
            CALL XDRSET(SCAT(1,NGFR+1,1),NGRO*NGTHER,0.0)
            DO I=1,NBTEM
              WW=TERPT(I)
              IF(ABS(WW).GT.1.0E-6) THEN
                CALL XSDTHE(7004,5015,-1,I,LOAD,IERR)
                IJ=0
                DO IG=1,NGTHER
                  IM=NINT(LOAD(IJ+2))
                  IG1=-NINT(LOAD(IJ+1))+NGFR+IG
                  IJ=IJ+2
                  DO I0=1,IM
*                   --- NGFR+IG is the primary group
                    SCAT(IG1+I0,NGFR+IG,1)=SCAT(IG1+I0,NGFR+IG,1)+WW*
     >              LOAD(IJ+I0)
                    GAR1(NGFR+IG,5)=GAR1(NGFR+IG,5)+WW*LOAD(IJ+I0)
                  ENDDO
                  IJ=IJ+IM
                ENDDO
              ENDIF
            ENDDO
          ENDIF
          IF(NP1.GT.0) THEN
*           P1 differential scattering XS
            CALL XSDISO(7002,5016,IND,LOAD,IERR)
            CALL XDRSET(GAR1(1,6),NGRO,0.0)
            CALL XDRSET(SCAT(1,1,2),NGRO*NGRO,0.0)
            IJ=0
            DO IG=1,NGFR
              IM=NINT(LOAD(IJ+2))
              IG1=-NINT(LOAD(IJ+1))+IG
              IJ=IJ+2
              DO I0=1,IM
*               --- IG is the primary group
                SCAT(IG1+I0,IG,2)=LOAD(IJ+I0)
                GAR1(IG,6)=GAR1(IG,6)+LOAD(IJ+I0)
              ENDDO
              IJ=IJ+IM
            ENDDO
            IF(NGTHER.GT.0) THEN
              CALL XDRSET(GAR1(NGFR+1,6),NGTHER,0.0)
              CALL XDRSET(SCAT(1,NGFR+1,2),NGRO*NGTHER,0.0)
              DO I=1,NBTEM
                WW=TERPT(I)
                IF(ABS(WW).GT.1.0E-6) THEN
                  CALL XSDTHE(7004,5016,-1,I,LOAD,IERR)
                  IJ=0
                  DO IG=1,NGTHER
                    IM=NINT(LOAD(IJ+2))
                    IG1=-NINT(LOAD(IJ+1))+NGFR+IG
                    IJ=IJ+2
                    DO I0=1,IM
*                     --- NGFR+IG is the primary group
                      SCAT(IG1+I0,NGFR+IG,2)=SCAT(IG1+I0,NGFR+IG,2)+WW*
     >                LOAD(IJ+I0)
                      GAR1(NGFR+IG,6)=GAR1(NGFR+IG,6)+WW*LOAD(IJ+I0)
                    ENDDO
                    IJ=IJ+IM
                  ENDDO
                ENDIF
              ENDDO
            ENDIF
          ENDIF
*----
*  Recover self-shielding data
*----
          ALLOCATE(WT0(NGRO))
          CALL XDRSET(WT0,NGRO,1.0)
          IF((NF.GE.1).AND.(NF.LE.3)) THEN
*
*           --- Recover Goldstein-Sehgal parameters
            ALLOCATE(GC(NGRO))
            CALL XDRSET(GC,NGRO,1.0)
            CALL XSDISO(7000,5012,IND,GC(NGF+1:),IERR)
            CALL LCMPUT(KPLIB,'NGOLD',NGRO,2,GC)
            DEALLOCATE(GC)
*
            CALL XSDRES(IND,IHEAD,IERR)
            NBTEM=IHEAD(1)
            NBDIL=IHEAD(2)
            IF(NBTEM.GT.MAXTMP) CALL XABORT('LIBND1: MAXTMP overflow')
            IF(NBDIL.GT.MAXDIL) CALL XABORT('LIBND1: MAXDIL overflow')
*
*           --- Temperature interpolation
            IF(NBTEM.EQ.1) THEN
              TERPT(1)=1.0
            ELSE
*             Resonance temperatures
              DO I=1,NBTEM
                TEMPS(I)=RHEAD(2+I)
              ENDDO
              CALL ALTERP(LCUBIC,NBTEM,TEMPS,TN(IMX),.FALSE.,TERPT)
            ENDIF
            ALLOCATE(RESD(MAXTDN))
            DO I=1,NBDIL
              DILUS(I)=RHEAD(2+NBTEM+I)
            ENDDO
            CALL XDRSET(XA,(NGFR-NGF)*NBDIL,0.0)
            CALL XDRSET(XS,(NGFR-NGF)*NBDIL,0.0)
            CALL XDRSET(XF,(NGFR-NGF)*NBDIL,0.0)
            CALL XDRSET(XN,(NGFR-NGF)*NBDIL,0.0)
            DO IG=1,NGFR-NGF
*             --- Absorption
              CALL XSDTAB(5004,IND,IG,RESD,IERR)
              DO I=1,NBTEM
                WW=TERPT(I)
                IF(ABS(WW).GT.1.0E-6) THEN
                  IOF=(I-1)*NBDIL
                  DO J=1,NBDIL
                    XA(IG,J)=XA(IG,J)+WW*RESD(IOF+J)
                  ENDDO
                ENDIF
              ENDDO
*             --- Scattering
              CALL XSDTAB(5015,IND,IG,RESD,IERR)
              DO I=1,NBTEM
                WW=TERPT(I)
                IF(ABS(WW).GT.1.0E-6) THEN
                  IOF=(I-1)*NBDIL
                  DO J=1,NBDIL
                    XS(IG,J)=XS(IG,J)+WW*RESD(IOF+J)
                  ENDDO
                ENDIF
              ENDDO
              IF(NF.EQ.3) THEN
*               --- Nu*Fission
                CALL XSDTAB(5006,IND,IG,RESD,IERR)
                DO I=1,NBTEM
                  WW=TERPT(I)
                  IF(ABS(WW).GT.1.0E-6) THEN
                    IOF=(I-1)*NBDIL
                    DO J=1,NBDIL
                      XF(IG,J)=XF(IG,J)+WW*RESD(IOF+J)
                    ENDDO
                  ENDIF
                ENDDO
              ENDIF
*             --- NJOY Flux
              CALL XSDTAB(5021,IND,IG,RESD,IERR)
              DO I=1,NBTEM
                WW=TERPT(I)
                IF(ABS(WW).GT.1.0E-6) THEN
                  IOF=(I-1)*NBDIL
                  DO J=1,NBDIL
                    XN(IG,J)=XN(IG,J)+WW*RESD(IOF+J)
                  ENDDO
                ENDIF
              ENDDO
            ENDDO
            DEALLOCATE(RESD)
*
*           --- Dilution interpolation and Livolant-Jeanpierre
*               normalization
            CALL LIBND3(NGF,NGFR,NGRO,NBDIL,SN(1,IMX),SB(1,IMX),DILUS,
     >      DELTA,NF,XA,XS,XF,XN,GAR1,SCAT(1,1,1),GAR2,WT0)
*
*           --- Apply self-shielding on SCAT and GAR1
            DO IG=NGF+1,NGFR
              DO IM=1,2
                WW=GAR2(IG,4+IM)/GAR1(IG,4+IM)
                DO IG1=1,NGRO
                  SCAT(IG1,IG,IM)=SCAT(IG1,IG,IM)*WW
                ENDDO
              ENDDO
              DO I=1,6
                GAR1(IG,I)=GAR2(IG,I)
              ENDDO
            ENDDO
          ENDIF
          CALL LCMPUT(KPLIB,'NWT0',NGRO,2,WT0)
          DEALLOCATE(WT0)
*----
*  Save xs information on the microlib
*----
          DO IG=1,NGRO
*           (n,g) xs
            GAR1(IG,7)=GAR1(IG,2)+GAR1(IG,7)-GAR1(IG,3)
*           Total xs
            GAR1(IG,2)=GAR1(IG,2)+GAR1(IG,5)
*           Transport correction
            GAR1(IG,1)=GAR1(IG,2)-GAR1(IG,1)
          ENDDO
          CALL LCMPUT(KPLIB,'TRANC',NGRO,2,GAR1(1,1))
          CALL LCMPUT(KPLIB,'NTOT0',NGRO,2,GAR1(1,2))
          CALL LCMPUT(KPLIB,'NG',NGRO,2,GAR1(1,7))
          IF(NF.EQ.3) THEN
            CALL LCMPUT(KPLIB,'NFTOT',NGRO,2,GAR1(1,3))
            CALL LCMPUT(KPLIB,'NUSIGF',NGRO,2,GAR1(1,4))
          ENDIF
          NL2=1
          IF((NL.GE.2).AND.(NP1.GT.0)) NL2=2
          CALL XDRLGS(KPLIB,1,IMPX,0,NL2-1,1,NGRO,GAR1(1,5),SCAT,ITYPRO)
          IF(IMPX.GT.5) CALL LCMLIB(KPLIB)
        ENDIF
      ENDDO
      CALL XSDCL()
*----
*  Scratch storage deallocation
*----
      DEALLOCATE(XN,XF,XS,XA,LOAD,THERXS,SCAT,GAR2,GAR1,TERPT,DILUS,
     1 TEMPS,DELTA)
      DEALLOCATE(HNAM)
      RETURN
*
  100 FORMAT(/21H NDAS LIBRARY OPTIONS/21H --------------------/
     1 7H NISOLB,I6,39H   (Number of isotopes in NDAS library)/
     2 7H NGRO  ,I6,28H   (Number of energy groups)/
     3 7H NGFIS ,I6,29H   (Number of fission groups)/
     4 7H NGF   ,I6,26H   (Number of fast groups)/
     5 7H NGRES ,I6,31H   (Number of resonance groups)/
     6 7H NGTHER,I6,29H   (Number of thermal groups)/
     7 7H NBFISS,I6,31H   (Number of fissile nuclides)/
     8 7H NBFP  ,I6,31H   (Number of fission products)/
     9 7H NBP1  ,I6,40H   (Number of nuclides with P1 matrices)/
     1 7H NBRES ,I6,33H   (Number of resonance nuclides)/
     2 7H MAXTMP,I6,40H   (Maximum number of temperature nodes)/
     3 7H MAXDIL,I6,37H   (Maximum number of dilution nodes)/
     4 7H MAXTDN,I6,36H   (Maximum number of product nodes)/
     5 7H IOLD  ,I6,32H   (Library type: old=1, new>=2)/
     6 7H MAXP0 ,I6,34H   (Maximum length of P0 matrices)/
     7 7H MAXP1 ,I6,34H   (Maximum length of P1 matrices))
  110 FORMAT(/30H Processing isotope/material ',A12,11H' (HNISOR=',A12,
     1 3H').)
  120 FORMAT(/16H ISOTOPE OPTIONS/16H ---------------/
     1 7H NBURN ,I6,46H   (Number of daughters in burnup calculation)/
     2 7H ID    ,I6,15H   (Numeric ID)/
     3 5H AW  ,1P E10.2,14H (Atomic mass)/
     4 7H IZ    ,I6,22H   (Number of protons)/
     5 7H NF    ,I6,24H   (Self-shielding flag)/
     6 7H NT    ,I6,38H   (Number of thermal xs temperatures)/
     7 7H NR    ,I6/
     8 7H NDAT2 ,I6/
     9 7H NDAT3 ,I6/
     1 7H NP1   ,I6,23H   (P1 scattering flag)/
     2 7H NS    ,I6,25H   (Fissile isotope flag)/
     3 7H IENDFB,I6,28H   (Type of evaluation file))
  130 FORMAT(26HLIBND1: Material/isotope ',A12,5H' = ',A12,9H' is miss,
     1 23Hing on NDAS file named ,A24,1H.)
      END
