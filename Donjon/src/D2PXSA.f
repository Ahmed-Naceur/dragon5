*DECK D2PXSA
      SUBROUTINE D2PXSA(IPDAT,IPSAP,ICAL,IPRINT,NGP,NREA,NISO,NMAC,
     1 NMIL,NANI,NVAR,NADRX,STAIDX,B2,ADF_T,NSF,LABS,SCAT,LADF)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover FISSION cross sections of an elementary calculation and store
* in INFO/BRANCH_INFO/MACROLIB_XS/SFI.
* WARNING: the GET_SFI_XS subroutine cannot recover FISSION XS in the
* case where cross sections are ineterpolated by the SCR: module
*
*Author(s): 
* J. Taforeau
*
*Parameters: input
* IPDAT   address of the INFO data block
* ICAL    number of the elementary calculation in which fission cross
*         sections is to be recovered
* IPSAP   address of the Saphyb object
* NGP     number of group energies in Saphyb
* NREA    number of reactions in Saphyb
* NISO    number of isotopes in Saphyb
* NMAC    number of macros in Saphyb
* NMIL    number of mixtures in Saphyb
* NANI    number of Legendre orders in Saphyb
* NADRX   concerne cross section vector (ADRX)
* STAIDX  index of current branch state values
* NSF     number of elements of the tranfert matrix
* LABS    content of absorption xs LABS(1) : ABS XS = TOTAL - SIGS00 ;
*         LABS(2) abs xs recovered from sap ; LABS (3)  abs xs recovered
*         from SAP minus excess xs

*Parameters: output
* SFI     fission cross sections of the current BRANCH:
*         INFO/BRANCH_INFO/MACROLIB_XS/SFI
*
*Parameters: 
* IPRINT  
* NVAR    
* B2      
* ADF_T   
* SCAT    
* LADF    
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPDAT,IPSAP
      INTEGER ICAL,IPRINT,NGP,NREA,NISO,NMAC, NMIL,NANI,NVAR,
     1 STAIDX(NVAR),NSF,NADRX
      REAL B2
      CHARACTER*3 ADF_T
      LOGICAL LABS(3),SCAT,LADF
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) IPTH,KPTH
      ! order numbers of current : reaction , isotope, macro
      INTEGER iprf,isot,imil,imac,iabs,iexc,idif,itra
      INTEGER idifc
      INTEGER ani, il, nj, ii, it, i1,i2, j1, j2, iadc,g
      INTEGER ND
      ! location of excess cross sections in RDATAX
      INTEGER ADR_EXC
      ! location of absorption cross sections in RDATAX
      INTEGER ADR_ABS
      ! location of profil cross sections in RDATAX
      INTEGER ADR_PRF
      ! location of TRANSFERT cross sections in RDATAX
      INTEGER ADR_TRA
      INTEGER ADR_DIF
      ! number of group energies in Saphyb
      INTEGER NG
      ! type of data recovered from GANLIB subroutines
      INTEGER ITYLCM
      ! name of : isotopes, macros
      CHARACTER(LEN=8) NOM_MAC(NMAC)
      ! name of reactions
      CHARACTER(LEN=10) NOM_REA(NREA)
      ! residual macro
      INTEGER RESMAC(NMIL)
      ! 3rd index of ADRX
      INTEGER ISADRX(NMIL)
      ! number of elements in RADATAX
      INTEGER LENGDX(NMIL)
      ! name of total macro
      INTEGER TOTMAC(NMIL)
      ! number of elements in IDATAP
      INTEGER LENGDP(NMIL)
      ! contains the adress of the 1st element in RDATAX
      INTEGER ADRX (NREA+2,NISO+NMAC,NADRX)
      REAL ABSORPTION(NGP)
      REAL TRANSFERT (NANI,NGP*NGP)
      REAL DIFC(NGP)
      CHARACTER(LEN=12) CALDIR
      INTEGER fagg, lagg,fdgg,wgal,fag,lag   ! CF SAPHTOOL MANUAL
      INTEGER fdg(NGP),adr(NGP+1)            ! CF SAPHTOOL MANUAL
      INTEGER NSCAT
      REAL CURRN(NSF,NGP,2)
      REAL SRFLX(NSF,NGP)
      REAL ZAFLX(NMIL,NGP)
      DOUBLE PRECISION RPAR (6,NSF)
      INTEGER IPAR (3,NSF)
      REAL ADF(NGP,NSF,10)
      REAL SCAT_MAT(NGP*NGP)
      ! transfert matrix
      INTEGER ,ALLOCATABLE, DIMENSION(:) :: IDATAP
      ! contains values of cross sections of an elementary calculation
      REAL,ALLOCATABLE, DIMENSION(:) :: RDATAX

      TRANSFERT(:,:) = 0

      WRITE(CALDIR,'("calc", I8)') ICAL
      CALL LCMSIX(IPSAP,' ',0)
      CALL LCMSIX(IPSAP,'contenu',1)
      IF(NMIL.NE.1) THEN
        ! the number of mixtures must be equal to one for converting
        ! Saphyb into PMAXS format
        CALL XABORT('@D2P: MORE THAN ONE MIXTRURE IN SAPHYB')
      ENDIF
      CALL LCMGTC(IPSAP,'NOMREA',10,NREA,NOM_REA)
      CALL LCMGTC(IPSAP,'NOMMAC',8,NMAC,NOM_MAC)
      CALL LCMGET(IPSAP,'RESMAC',RESMAC)
      CALL LCMGET(IPSAP,'TOTMAC',TOTMAC)
      CALL LCMSIX(IPSAP,' ',0)

      CALL LCMSIX(IPSAP,'adresses',1)
      CALL LCMGET(IPSAP,'ADRX',ADRX)
      CALL LCMSIX(IPSAP,' ',0)
      CALL LCMSIX(IPSAP,CALDIR,1)
      CALL LCMSIX(IPSAP,'info',1)
      CALL LCMGET(IPSAP,'ISADRX',ISADRX)
      CALL LCMGET(IPSAP,'LENGDX',LENGDX)
      CALL LCMGET(IPSAP,'LENGDP',LENGDP)
      ALLOCATE (RDATAX(LENGDX(1)),IDATAP(LENGDP(1)))
      imac=0
      IF(RESMAC(1).NE.0) THEN
        imac=RESMAC(1)      ! recover name of residual macro
      ELSE IF(TOTMAC(1).NE.0) THEN
        imac=TOTMAC(1)      ! recover name of total macro
      ELSE
        CALL XABORT('@D2P: NO MACRO DEFINED')
      ENDIF
      isot=NISO+imac     ! we interest in macro fission cross sections
      imil=1             ! set the mixtures number to 1
      iprf=0
      iexc=0
      iabs=0
      idif=0
      iadc=0
      itra=0
      idifc=0
      !TEST HFATC
      NSCAT=1

      DO ir=1,NREA
        ! store the order numbers of PROFIL matrix
        IF((SCAT) .and. NOM_REA(ir)=="PROFIL") iprf=ir
        IF((SCAT) .and. NOM_REA(ir)=="DIFFUSION") idif=ir
        IF((SCAT) .and. NOM_REA(ir)=="TRANSFERT") itra=ir
        IF(NOM_REA(ir)=="NU*FISSION") iabs=ir
        ! store the order numbers of EXCESS matrix
        IF(LABS(3).and. NOM_REA(ir)=="EXCESS")  iexc=ir
        IF((LADF) .and. NOM_REA(ir)=="FUITES") idifc=ir
      ENDDO
      IF(iabs==0) CALL XABORT ('@D2P: NO ABSORPTION XS AVAILABLE')
      IF(LABS(3).and.iexc==0) THEN
        CALL XABORT('@D2P: NO EXCESS XS AVAILABLE')
      ENDIF
      IF(SCAT .and. iprf==0) THEN
        CALL XABORT('@D2P: NO PROFIL XS AVAILABLE')
      ENDIF
      IF(SCAT .and. idif==0) THEN
        CALL XABORT('@D2P: NO DIFFUSION XS AVAILABLE')
      ENDIF
      IF(SCAT .and. itra==0) THEN
        CALL XABORT('@D2P: NO TRANSFERT XS AVAILABLE')
      ENDIF
      IF((LADF) .and. idifc==0) THEN
        CALL XABORT('@D2P: NO FUITES XS AVAILABLE')
      ENDIF
      NANI=ADRX(NREA+2,isot,ISADRX(imil))-1
      ND=ADRX(NREA+1,isot,ISADRX(imil))
      IF(MOD(idif,NREA+1).GT.0  .AND. ND.GE.1) THEN
        iadc=ADRX(idif,isot,ISADRX(imil))+NGP
      ENDIF
      ! address in RDATAX of ABSORPTION XS
      ADR_ABS=ADRX(iabs,isot,ISADRX(imil))


      ! address in RDATAX of EXCESS XS
      ADR_EXC=0
      IF(LABS(3)) ADR_EXC=ADRX(iexc,isot,ISADRX(imil))
      ! address in RDATAX of PROFIL XS
      ADR_PRF=0
      IF(SCAT) ADR_PRF=ADRX(iprf,isot,ISADRX(imil))
      ! address in RDATAX of TRANSFERT XS
      ADR_TRA=0
      IF(SCAT) ADR_TRA=ADRX(itra,isot,ISADRX(imil))
      ! address in RDATAX of FUITES XS
      ADR_DIF=0
      IF(LADF) ADR_DIF=ADRX(idifc,isot,ISADRX(imil))

      ! moving in the saphyyb object to recover RDATAX information
      CALL LCMSIX(IPSAP,' ',0)
      CALL LCMSIX(IPSAP,CALDIR,1)
      CALL LCMSIX(IPSAP,'mili       1',1)
      CALL LCMGET(IPSAP,'RDATAX',RDATAX)
      CALL LCMGET(IPSAP,'IDATAP',IDATAP)

      ! LOOP over energy groups
      DO ig=1, NGP
       ABSORPTION(ig)=RDATAX(ADR_ABS+ig-1)
       IF(LADF) DIFC(ig)=RDATAX(ADR_DIF+ig-1)
       IF(LABS(3)) THEN
         ABSORPTION(ig)=ABSORPTION(ig)-RDATAX(ADR_EXC+ig-1)
       ENDIF
      ENDDO


      IF(SCAT)THEN   ! recover the scattering XS from Saphyb
         ii = ADR_PRF
         nj = IDATAP(ii+6+2*NGP)-1

        DO ani=0, NANI
         il = ADR_TRA + (ani) * nj
         fagg =IDATAP(ii)
         lagg =IDATAP(ii+1)
         fdgg =IDATAP(ii+2)
         wgal =IDATAP(ii+3)
         fag  =IDATAP(ii+4)
         lag  =IDATAP(ii+5)
         fdg  =IDATAP(ii+6:ii+5+NGP)
         adr  =IDATAP(ii+6+NGP:ii+6+2*NGP)
         IF(wgal.GT.0)THEN
          it=il
          DO g=fagg,lagg
           i1=(g-1)*NGP+fdgg
           i2=(g-1)*NGP+fdgg+wgal-1
           TRANSFERT(ani+1,i1:i2)=RDATAX(it:it+wgal-1)
           it=it+wgal
          ENDDO
         ENDIF
         DO g=fag,lag
          i1=(g-1)*NGP+fdg(g)
          i2=(g-1)*NGP+fdg(g)+adr(g+1)-adr(g)-1
          j1=il-1+adr(g)
          j2=il-1+adr(g+1)-1
          TRANSFERT(ani+1,i1:i2)=RDATAX(j1:j2)
         ENDDO
        ENDDO
        IF(iadc.NE.0)THEN
         NG=NGP
!         TRANSFERT(1,1:NG*NG:NG+1) =
!     >         TRANSFERT(1,1:NG*NG:NG+1) - RDATAX(iadc:iadc+NGP-1)
        ENDIF
      ENDIF
      DO g=1, NGP
        DO ig=1, NGP
         SCAT_MAT(NSCAT) = TRANSFERT(1,g+(ig-1)*NGP)
         NSCAT=NSCAT+1
        ENDDO
      ENDDO

      ! RECOVER ADF IN SAPHYB (IF AVAILABLE) (adapted from !
      ! saphyb_browser of UPM)

      IF(LADF) THEN
        IF((ADF_T.EQ.'SEL').OR.(ADF_T.EQ.'GET')) THEN
         CALL LCMSIX (IPSAP,' ',0)
         CALL LCMSIX (IPSAP,'geom        ',1)
         CALL LCMSIX (IPSAP,'outgeom     ',1)
         CALL LCMLEN(IPSAP,'SURF',NSURF,ITYLCM)
         IF(NSF.NE.NSURF) THEN
           WRITE(6,*) "@D2P: ERROR IN NUMBER OF ASSEMBLY SURFACES"
           WRITE(6,*) "THE NUMBER OF SURFACES IN SAP (",
     1     NSURF,") IF DIFFERENT FROM DRAG2PARCS INPUT (",NSF,")"
           CALL XABORT('')
         ENDIF

         CALL LCMGET(IPSAP,'IPAR',IPAR)
         CALL LCMGET(IPSAP,'RPAR',RPAR)

         CALL LCMSIX(IPSAP,' ',0)
         CALL LCMSIX(IPSAP,CALDIR,1)
         CALL LCMSIX(IPSAP,'outflx      ',1)

         CALL LCMGET(IPSAP,'CURRM',CURRN(:,:,2))
         CALL LCMGET(IPSAP,'CURRP',CURRN(:,:,1))
         CALL LCMGET(IPSAP,'SURFLX',SRFLX(:,:))
         CALL LCMGET(IPSAP,'REGFLX',ZAFLX(:,:))
         ADF = 0.
         DIFC(:)=DIFC(:)/B2

         ! CALL to GET_SFI_XS to recover ADF
         CALL D2PADF(IPDAT,IPRINT,NGP,NMIL, ADF, NSF, DIFC,CURRN,SRFLX,
     1   ZAFLX,RPAR,IPAR,ADF_T,STAIDX,NVAR)
        ENDIF
      ENDIF
      ! STORE RESULTS IN INFO DATA BLOCK
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'BRANCH_INFO',1)
      IPTH=LCMGID(IPDAT,'CROSS_SECT')
      KPTH=LCMDIL(IPTH,STAIDX(NVAR))
      CALL LCMSIX(KPTH,'MACROLIB_XS',1)
      IF(LABS(2)) CALL LCMPUT(KPTH,'ABSORPTION',NGP,2,ABSORPTION)
      IF(SCAT) CALL LCMPUT(KPTH,'SCAT',NGP*NGP,2,SCAT_MAT)
      IF(LABS(2)) WRITE(6,*) "ABSORPTION EXCESS :", ABSORPTION
      IF(SCAT) WRITE(6,*) "SCATTERING MATRIX  :", SCAT_MAT
      DEALLOCATE (RDATAX,IDATAP)
      END
