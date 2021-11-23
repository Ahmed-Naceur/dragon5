*DECK D2PRFL
      SUBROUTINE D2PRFL(  IPDAT, IPMIC , IPRINT,    NBU,   NGP,   NBMIX,
     >                    NANI,   NVAR,  STAIDX,   LADF,  NADF,   NTYPE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover macroscopic and microscopic cross sections from a microlib
* object and write cross sections for one branch at a fixed burnup point
* in the INFO data block.
* WARNING: 04/2014 The information recovered by this routine is exactly
* the same than GET_MACROLIB_XS but is used for reflector case, in this
* case the following reactions are set to zero :
*  DET(IGR) = 0
*  SFI(IGR) = 0
*  KAPPA_FI(IGR)= 0
*  FLUX(IGR) = 0
*  VELINV(IGR) = 0
*  CHI_SPEC(IGR) = 0
*  X_NU_FI(IGR) = 0
*  KAPPA_FI(IGR) = 0
*  XENG(IGR)=0
*  SMNG(IGR)=0
* NB : for reflector case, the upscattering is fixed to zero
*
*Author(s): 
* J. Taforeau
*
*Parameters: input
* IPDAT   address of info data block
* IPMIC   address of the microlib object
* NBU     number of burnup points
* NBMIX   number of mixturess
* NGP     number of energy groups
* NANI    number of anisotropy
* NVAR    number of state variables
* STAIDX  table of states index order
* NADF    number of ADF to be recovered
* NTYPE   number  of adf type
* LADF    flag for adf
*
*Parameters: 
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPDAT,IPMIC
      INTEGER STAIDX(NVAR)
      INTEGER NBU,NVAR,NBMIX,NGP,NANI,NADF,NTYPE
      LOGICAL LADF
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPMIC,KPMIC,IPTH,KPTH
      INTEGER NSCAT,MIX
      INTEGER IJJ(NBMIX),NJJ(NBMIX),IPOS(NBMIX)
      REAL GAR2(NGP,NGP,NBMIX,NANI),GAR3(NBMIX*NGP)
      REAL XSECT(NGP,NBMIX)          ! TOTAL CROSS SECTIONS
      REAL KAPPA_FI(NGP)             ! KAPPA FISSION CROSS SECTIONS
      REAL X_NU_FI(NGP)              ! NU SIGMA FISSION CROSS SECTIONS
      REAL XTR(NGP)                  ! TRANSPORT CROSS SECTIONS
      REAL DIFF(NGP,NBMIX)           ! DIFFUSION COEFF
      REAL SCAT(NGP,NBMIX)           ! SCATTERING CROSS SECTIONS
      REAL DET(NGP)                  ! DETECTOR CROSS SECTIONS
      REAL SFI(NGP)                  ! FISSION CROSS SECTIONS
      REAL ABSORPTION(NGP)           ! ABSORPTION CROSS SECTIONS
      REAL SCAT_MAT(NGP*NGP)         ! SCATTERING MATRIX
      REAL SCAT_TMP(NGP,NGP,NBMIX,NANI)  ! TEMPORARY SCATTERING MATRIX
      REAL FLUX(NGP)
      REAL VELINV(NGP)
      REAL XENG(NGP)
      REAL CHI_SPEC(NGP),VOLUME(NBMIX)
      REAL SMNG(NGP),FLXHET(NGP*NBMIX),FLXHOM(NGP,NBMIX)
      REAL FLXL(NGP),FLXR(NGP),CURL(NGP),CURR(NGP)
      REAL ADF(NADF,NGP)
      CHARACTER CM*2,ADF_T*3
      CHARACTER*8 ADFD(NADF),HADF(NTYPE),HFLX(2),HCUR(2)
      IF(IPRINT > 0)  THEN
        WRITE(6,*)
        WRITE(6,*) "****** RECOVER REFLECTOR CROSS SECTIONS ******"
      ENDIF
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'SAPHYB_INFO',1)
      IF(LADF) THEN
       ADF_T="   "
       CALL LCMGTC(IPDAT,'ADF_TYPE',3,1,ADF_T)
       IF ((ADF_T.NE.'DRA').AND.(ADF_T.NE.'GEN')) THEN
        WRITE(6,*)'@D2PRFL:',ADF_T,'ADF NOT SUPPORTED ',
     >  'WITH REFL CALCULATION'
        CALL XABORT('')
       ENDIF
       IF ((ADF_T.EQ.'DRA')) THEN
        CALL LCMGTC(IPDAT,'HADF',8,NADF,ADFD)
       ELSE IF ((ADF_T.EQ.'GEN')) THEN
        CALL LCMGTC(IPDAT,'HFLX',8,2,HFLX)
        CALL LCMGTC(IPDAT,'HCUR',8,2,HCUR)
       ENDIF

      ENDIF

      CALL LCMGET(IPDAT,'MIX',MIX)
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPMIC,' ',0)
      CALL LCMSIX(IPMIC,'MACROLIB',1)
      CALL LCMGET(IPMIC,'VOLUME',VOLUME)

      IF (LADF) THEN
        CALL LCMSIX(IPMIC,'ADF',1)
        CALL LCMGTC(IPMIC,'HADF',8,NTYPE,HADF)
        ITYPE=1
       IF ((ADF_T.EQ.'DRA')) THEN
        DO ITYPE=1,NTYPE
         CALL LCMGET(IPMIC,HADF(ITYPE),FLXHET)
         DO I=1,NADF
          IF(HADF(ITYPE).EQ.ADFD(I))THEN
           DO IGR=1, NGP
             ADF(I,IGR)= FLXHET((IGR-1)*NBMIX+MIX)
           ENDDO
          ENDIF
         ENDDO
        ENDDO
       ELSE IF ((ADF_T.EQ.'GEN')) THEN
        DO ITYPE=1,NTYPE
         CALL LCMGET(IPMIC,HADF(ITYPE),FLXHET)
         IF(HADF(ITYPE).EQ.HFLX(1))THEN
          FLXL(:)=FLXHET
         ENDIF
         IF (HADF(ITYPE).EQ.HFLX(2))THEN
          FLXR(:)=FLXHET
         ENDIF
         IF (HADF(ITYPE).EQ.HCUR(1))THEN
          CURL(:)=FLXHET
         ENDIF
         IF (HADF(ITYPE).EQ.HCUR(2))THEN
          CURR(:)=FLXHET
         ENDIF
        ENDDO
       ENDIF
       CALL LCMSIX(IPMIC,'',2)

      ENDIF

      JPMIC=LCMGID(IPMIC,'GROUP')

      !  RECOVER CROSS SECTIONS INFORMATION
      DO IGR=1,NGP
       WRITE(6,'(/28H PROCESS ENERGY GROUP NUMBER,I4)') IGR
       KPMIC=LCMGIL(JPMIC,IGR)
       CALL LCMLEN(KPMIC,'NTOT0',ILONG,ITYLCM)

       IF(ILONG.NE.NBMIX) THEN
         CALL XABORT('D2P: MORE THAN ONE MIXTURE IN SAPHYB')
       ENDIF
       CALL LCMGET(KPMIC,'FLUX-INTG',FLXHOM(IGR,1:NBMIX))
       CALL LCMGET(KPMIC,'NTOT0',XSECT(IGR,1:NBMIX))
       CALL LCMGET(KPMIC,'SIGS00',SCAT(IGR,1:NBMIX))
       CALL LCMGET(KPMIC,'DIFF',DIFF(IGR,1:NBMIX))
       ABSORPTION(IGR)=XSECT(IGR,MIX)-SCAT(IGR,MIX)
       IF (LADF) ADF(:,IGR)= VOLUME * ADF(:,IGR) / FLXHOM(IGR,MIX)
       DET(IGR) = 0
       SFI(IGR) = 0
       KAPPA_FI(IGR)= 0
       FLUX(IGR) = 0
       VELINV(IGR) = 0
       CHI_SPEC(IGR) = 0
       X_NU_FI(IGR) = 0
       KAPPA_FI(IGR) = 0
       XENG(IGR)=0
       SMNG(IGR)=0
       XTR(IGR)=1/(3*DIFF(IGR,MIX))

       CALL XDRSET(GAR2,NGP*NGP*NBMIX*NANI,0.0)
       DO IL=1,NANI
          WRITE(CM,'(I2.2)') IL-1
          LENGTH=1
          IF(IL.GT.1) CALL LCMLEN(KPMIC,'SCAT'//CM,LENGTH,ITYLCM)
          IF(LENGTH.GT.0) THEN
           CALL LCMGET(KPMIC,'SCAT'//CM,GAR3)
           CALL LCMGET(KPMIC,'NJJS'//CM,NJJ)
           CALL LCMGET(KPMIC,'IJJS'//CM,IJJ)
           CALL LCMGET(KPMIC,'IPOS'//CM,IPOS)
           DO IMIL=1,NBMIX
            IPOSDE=IPOS(IMIL)
            DO JGR=IJJ(IMIL),IJJ(IMIL)-NJJ(IMIL)+1,-1
             GAR2(IGR,JGR,IMIL,IL)=GAR3(IPOSDE) ! IGR <-- JGR
             SCAT_TMP(IGR,JGR,IMIL,IL)=GAR2(IGR,JGR,IMIL,IL)
             IPOSDE=IPOSDE+1
            ENDDO
           ENDDO
          ENDIF
         ENDDO
      ENDDO

      NSCAT=1
      DO J=1, NGP
         DO I=1, NGP
          SCAT_MAT(NSCAT)=SCAT_TMP(I,J,MIX,1) ! I <-- J
          IF(NSCAT==3) SCAT_MAT(NSCAT)=0
          NSCAT=NSCAT+1
         ENDDO
      ENDDO

      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'BRANCH_INFO',1)
      IF(STAIDX(NVAR)==1) THEN
        IPTH=LCMLID(IPDAT,'CROSS_SECT',NBU)
      ELSE
        IPTH=LCMGID(IPDAT,'CROSS_SECT')
      ENDIF

      KPTH=LCMDIL(IPTH,STAIDX(NVAR))
      CALL LCMSIX(KPTH,'MICROLIB_XS',1)

      CALL LCMPUT(KPTH,'XENG',NGP,2,XENG)
      CALL LCMPUT(KPTH,'SMNG',NGP,2,SMNG)

      CALL LCMSIX(KPTH,' ',2)
      CALL LCMSIX(KPTH,'MACROLIB_XS',1)

      CALL LCMPUT(KPTH,'XTR',NGP,2,XTR)
      CALL LCMPUT(KPTH,'ABSORPTION',NGP,2,ABSORPTION)
      CALL LCMPUT(KPTH,'X_NU_FI',NGP,2,X_NU_FI)
      CALL LCMPUT(KPTH,'KAPPA_FI',NGP,2,KAPPA_FI)
      CALL LCMPUT(KPTH,'SFI',NGP,2,SFI)
      CALL LCMPUT(KPTH,'DET',NGP,2,DET)
      CALL LCMPUT(KPTH,'SCAT',NGP*NGP,2,SCAT_MAT)
      IF (LADF) THEN
       IF (ADF_T.EQ.'DRA') THEN
        CALL LCMPUT(KPTH,'ADF',NADF*NGP,2,ADF)
       ELSE IF (ADF_T.EQ.'GEN') THEN
        CALL LCMPUT(KPTH,'FLXL',NGP,2,FLXL)
        CALL LCMPUT(KPTH,'FLXR',NGP,2,FLXR)
        CALL LCMPUT(KPTH,'CURL',NGP,2,CURL)
        CALL LCMPUT(KPTH,'CURR',NGP,2,CURR)
       ENDIF
      ENDIF
      IF(IPRINT>1) THEN
         WRITE(6,*)
         WRITE(6,*) "**** MACROSCOPIC cross sections (1:NGP) ****"
         WRITE(6,*) "TOTALE                  :",XSECT(:,MIX)
         WRITE(6,*) "DIFFUSION               :",DIFF(:,MIX)
         WRITE(6,*) "TRANSPORT               :",XTR
         WRITE(6,*) "ABSORPTION              :",ABSORPTION
         WRITE(6,*) "NU FISSION              :",X_NU_FI
         WRITE(6,*) "KAPPA FISSION           :",KAPPA_FI
         WRITE(6,*) "DETECTOR                :",DET
         WRITE(6,*) "SCATTERING (g to g')    :",SCAT_MAT
         IF (LADF) THEN
          IF (ADF_T.EQ.'DRA') THEN
           WRITE(6,*) "ADF([N/E/W/S]||[W/E])   :",ADF
          ELSE IF (ADF_T.EQ.'GEN') THEN
           WRITE(6,*) "WEST FLUX BOUNDARY      :",FLXL
           WRITE(6,*) "EST FLUX BOUNDARY       :",FLXR
           WRITE(6,*) "WEST CURRENT BOUNDARY   :",CURL
           WRITE(6,*) "EST CURRENT BOUNDARY    :",CURR
          ENDIF
         ENDIF
      ENDIF
      END
