*DECK D2PMAC
      SUBROUTINE D2PMAC(  IPDAT, IPMIC , IPRINT,    NBU,   NGP,   NBMIX,
     >                     NADD,   NANI,  NVAR,  STAIDX,  LADF,    NADF,
     >                     NTYPE,  LCDF,  NCDF,    LGFF,  NGFF,    NPIN,
     >                      FLUX                                       )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover macroscopic cross sections from a microlib object and write
* cross sections for one branch at a fixed burnup point in the INFO
* data block.
*
*Author(s): 
* J. Taforeau
*
*Parameters: input
* IPDAT   address of info data block
* IPMIC   address of the microlib object
* IPRINT
* NBU     number of burnup points
* NBMIX   number of mixturess
* NBISO   number of isotopes
* NGP     number of energy groups
* NADD    number of additional cross sections
* NDEL    number of delayed neutron groups
* NANI    number of anisotropy
* NVAR    number of state variables
* STAIDX  table of states index order
* LADF    flag for assembly discontinuity factor
* NADF    number of assembly discontinuity factor per energy groups
* NTYPE   number of type of assembly discontinuity factor
* LCDF    flag for corner discontinuity factor
* NCDF    number of corner discontinuity factor per energy groups
* LGFF    flag for group form factor
* NGFF    number of group form factor per energy groups
* NPIN    number of pin on each side of the assembly
*         (note: if NADF, NCDF, NGFF or NPIN are not defined
*           a fake value of 1 is assigned for allocation memory issue)
*
*Parameters: 
* FLUX
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPDAT,IPMIC
      INTEGER STAIDX(NVAR)
      INTEGER NBU,NADD,NVAR,NBMIX,NGP,NANI,NADF,NCDF,NGFF,NPIN
      LOGICAL LADF,LCDF,LGFF
      REAL FLUX (NGP)
*----
*  LOCAL VARIABLES
*----
      INTEGER   NSTATE
      PARAMETER (NSTATE=40)
      TYPE(C_PTR) JPMIC,KPMIC,IPTH,KPTH
      INTEGER NSCAT,ITYLCM,ILONG,IUPS
      INTEGER IJJ(NBMIX),NJJ(NBMIX),IPOS(NBMIX)
      REAL GAR2(NGP,NGP,NBMIX,NANI),GAR3(NBMIX*NGP)
      REAL XSECT(NGP)                ! TOTAL CROSS SECTIONS
      REAL KAPPA_FI(NGP)             ! KAPPA FISSION CROSS SECTIONS
      REAL X_NU_FI(NGP)              ! NU SIGMA FISSION CROSS SECTIONS
      REAL XTR(NGP)                  ! TRANSPORT CROSS SECTIONS
      REAL DIFF(NGP)                 ! DIFFUSION COEFF
      REAL SCAT(NGP)                 ! SCATTERING CROSS SECTIONS
      !REAL TRANC(NGP)                ! TRANSPORT CORRECTION
      REAL ABSORPTION(NGP)           ! ABSORPTION CROSS SECTIONS
      REAL SCAT_MAT(NGP*NGP)         ! SCATTERING MATRIX
      REAL SCAT_TMP(NGP,NGP,NBMIX,NANI)  ! TEMPORARY SCATTERING MATRIX
      REAL SIGW00(NGP)
      DOUBLE PRECISION SUMSCAT(NGP)
      ! AVERAGE HOMOGENE SURFACIC FLUX (FLUX-INTG/VOLUME) and
      ! HETEROGENE
      REAL FLXHOM(NGP),FLXHET(NGP)
      REAL VOLUME
      CHARACTER(len=8) ADDXSNAM(NADD)
      CHARACTER*8 :: HFLX(8) = 'NUL'
      CHARACTER*8 :: HCUR(8) = 'NUL'
      CHARACTER CM*2,ADF_T*3,CDF_T*3,GFF_T*3
      CHARACTER(LEN=8) ADFD(NADF),CDFD(NCDF)
      CHARACTER(LEN=8) HADF(NTYPE)   ! ADF NAME IN MACROLIB
      REAL ADF(NADF,NGP)             ! ASSEMBLY AND CORNER DF
                                     ! NADF=1 for DRA,  NTYPE=1 for SEL
                                     ! and GET
      REAL CDF(NCDF,NGP)             ! ASSEMBLY AND CORNER DF
      REAL GFFC(NGFF,NGP)            ! GROUP FORM FACTORS GFF by mixture
      REAL KFC(NGFF,NGP)             !   h-factor
!      REAL VOLG(NGFF)                !   volume of GROUP FORM FACTORS
      REAL GFF(NPIN,NPIN,NGP)        ! GFF pin by pin
            ! GFF geometry
      INTEGER MIXG(NPIN,NPIN)        ! mixture

      CALL LCMSIX(IPMIC,' ',0)
      CALL LCMSIX(IPMIC,'MACROLIB',1)
      CALL LCMGET(IPMIC,'VOLUME',VOLUME)
      IF(NADD.GT.0)CALL LCMGTC(IPMIC,'ADDXSNAME-P0',8,NADD,ADDXSNAM)
      IF(NBMIX.NE.1) THEN
         ! SAPHYB MUST CONTAIN ONLY ONE MIXTURES
         CALL XABORT('D2PMAC: MORE THAN ONE MIXTURE IN SAPHYB')
      ENDIF
      JPMIC=LCMGID(IPMIC,'GROUP')
      SUMSCAT=0.0D0
      CALL XDRSET(SCAT_TMP,NGP*NGP*NBMIX*NANI,0.0)
      ! LOOP OVER ENERGY GROUPS
      DO IGR=1,NGP
        KPMIC=LCMGIL(JPMIC,IGR)
        CALL LCMLEN(KPMIC,'NTOT0',ILONG,ITYLCM)
        IF(ILONG.NE.NBMIX) THEN
          CALL XABORT('@D2PMAC: MORE THAN ONE MIXTURE IN SAP/MCO')
        ENDIF
        ! RECOVER CROSS SECTIONS INFORMATION
        CALL LCMGET(KPMIC,'NTOT0',XSECT(IGR))
        CALL LCMGET(KPMIC,'SIGS00',SCAT(IGR))
        CALL LCMGET(KPMIC,'SIGW00',SIGW00(IGR))
       ! CALL LCMGET(KPMIC,'TRANC',TRANC(IGR))
        CALL LCMGET(KPMIC,'NUSIGF',X_NU_FI(IGR))

        CALL LCMGET(KPMIC,'H-FACTOR',KAPPA_FI(IGR))
        CALL LCMLEN(KPMIC,'DIFF',ILONG,ITYLCM)
        IF (ILONG>0) THEN
         PRINT*,'ILONG DIFF ',ILONG
         CALL LCMGET(KPMIC,'DIFF',DIFF(IGR))
         XTR(IGR)=1/(3*DIFF(IGR))
        ELSE
         DIFF(:)=0
         CALL LCMLEN(KPMIC,'NTOT1',ILONG,ITYLCM)
         IF (ILONG.EQ.NGP) THEN
          CALL LCMGET(KPMIC,'NTOT1',XTR(IGR))
          WRITE(6,*) "WARNING : NTOT1 RECOVERED AS TRANSPORT
     > CROSS SECTION (SUITABLE FOR SPn WITH NG>=2)"
         ELSE
          CALL LCMGET(KPMIC,'NTOT0',XTR(IGR))
          WRITE(6,*) "WARNING : NTOT0 RECOVERED AS TRANSPORT
     > CROSS SECTION (SUITABLE FOR SPn WITH NG>2)"
         ENDIF
        ENDIF

        CALL LCMGET(KPMIC,'FLUX-INTG',FLXHOM(IGR))

        ! INITIALIZATION OF GAR2 VECTOR
        CALL XDRSET(GAR2,NGP*NGP*NBMIX*NANI,0.0)

        ! LOOP OVER ANISOTROPY COMPONENT
        DO IL=1,NANI
          WRITE(CM,'(I2.2)') IL-1
          LENGTH=1
          IF(IL.GT.1) CALL LCMLEN(KPMIC,'SCAT'//CM,LENGTH,ITYLCM)
          IF(LENGTH.GT.0) THEN
           CALL LCMGET(KPMIC,'SCAT'//CM,GAR3)
           CALL LCMGET(KPMIC,'NJJS'//CM,NJJ)
           CALL LCMGET(KPMIC,'IJJS'//CM,IJJ)
           CALL LCMGET(KPMIC,'IPOS'//CM,IPOS)
            ! LOOP OVER MIXTRURE
            DO IMIL=1,NBMIX
             IPOSDE=IPOS(IMIL)
             ! LOOP OVER ENERGY GROUPS
             DO JGR=IJJ(IMIL),IJJ(IMIL)-NJJ(IMIL)+1,-1
              GAR2(IGR,JGR,IMIL,IL)=GAR3(IPOSDE) ! IGR <-- JGR

              ! ELEMENTS OF THE SCATTERING MATRIX
              SCAT_TMP(IGR,JGR,IMIL,IL)=GAR2(IGR,JGR,IMIL,IL)
              IPOSDE=IPOSDE+1
             ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDDO


      ! STORE THE SCATTERING MATRIX CORRESPONDING TO L=0 AND MIX=1
      ! IN SCAT_MAT
      NSCAT=1
      DO J=1, NGP
         DO I=1, NGP

          SCAT_MAT(NSCAT)=SCAT_TMP(J,I,1,1) ! I <-- J 1<-1 2<-1

          IF (SCAT_MAT(NSCAT)<0) THEN
           SUMSCAT(J)=SUMSCAT(J)+SCAT_MAT(NSCAT)
           SCAT_MAT(NSCAT)=0
           WRITE(6,*) "WARNING : NEGATIVE VALUES FOR SCATTERING MATRIX
     >     ELEMENT (",J,"->",I,")."
          ENDIF
          NSCAT=NSCAT+1
         ENDDO
         XTR(J)=XTR(J)+REAL(SUMSCAT(J))
         SUMSCAT=0.0D0

      ENDDO

      DO I=1, NGP
       ABSORPTION(I)=XSECT(I)-SCAT(I)
      ENDDO

      ! STORE CROSS SECTIONS IN INFO/CROSS_SECT/MACROLIB_XS
      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'GENPMAXS_INP',1)
      CALL LCMGET(IPDAT,'IUPS',IUPS)
      IF ((IUPS.EQ.2).AND.(NGP.EQ.2)) THEN
       SCAT_MAT(2)=SCAT_MAT(2)-FLXHOM(2)/FLXHOM(1)*SCAT_MAT(3)
       SCAT_MAT(3)=0.
      ENDIF

      CALL LCMSIX(IPDAT,' ',0)
      CALL LCMSIX(IPDAT,'BRANCH_INFO',1)

      IF(STAIDX(NVAR)==1) THEN
        IPTH=LCMLID(IPDAT,'CROSS_SECT',NBU)
      ELSE
        IPTH=LCMGID(IPDAT,'CROSS_SECT')
      ENDIF

      KPTH=LCMDIL(IPTH,STAIDX(NVAR))

      CALL LCMSIX(KPTH,'MACROLIB_XS',1)
      CALL LCMPUT(KPTH,'XTR',NGP,2,XTR)
      CALL LCMPUT(KPTH,'ABSORPTION',NGP,2,ABSORPTION)
      CALL LCMPUT(KPTH,'X_NU_FI',NGP,2,X_NU_FI)
      CALL LCMPUT(KPTH,'KAPPA_FI',NGP,2,KAPPA_FI)
      CALL LCMPUT(KPTH,'SCAT',NGP*NGP,2,SCAT_MAT)

      ! RECOVER THE ASSEMBLY DISCONTINUITY FACTOR IF ADF DRA IS SET
      ! BY THE USER
      IF((LADF).OR.(LCDF)) THEN
        FLXHOM(:)=FLXHOM(:) / VOLUME
        CALL LCMSIX (IPDAT,' ',0)
        CALL LCMSIX (IPDAT,'SAPHYB_INFO',1)
        ADF_T="   "
        IF(LADF) THEN
          CALL LCMGTC(IPDAT,'ADF_TYPE',3,1,ADF_T)
          IF (ADF_T.EQ.'DRA') THEN
           CALL LCMGTC(IPDAT,'HADF',8*NADF,1,ADFD)
          ELSE IF (ADF_T.EQ.'GEN') THEN
           CALL LCMLEN(IPDAT,'HFLX',NFLX,ITYLCM)
           CALL LCMGTC(IPDAT,'HFLX',8*NFLX,1,HFLX(1:NFLX))
           CALL LCMGTC(IPDAT,'HCUR',8*NFLX,1,HCUR(1:NFLX))
          ENDIF
        ENDIF
        CDF_T="   "
        IF(LCDF) THEN
          CALL LCMGTC(IPDAT,'CDF_TYPE',3,1,CDF_T)
          CALL LCMGTC(IPDAT,'HCDF',8*NCDF,1,CDFD)
        ENDIF
        IF((ADF_T(:3) .EQ. 'DRA').OR.(CDF_T(:3) .EQ. 'DRA')
     >         .OR.(ADF_T(:3) .EQ. 'GEN' ) )THEN
           ! NADF = 1 or 4, NCDF = 1 or 4
          CALL LCMSIX(IPMIC,' ',0)
          CALL LCMSIX(IPMIC,'MACROLIB',1)
          CALL LCMSIX(IPMIC,'ADF',1)
          CALL LCMGTC(IPMIC,'HADF',8,NTYPE,HADF)
          DO ITYPE=1,NTYPE
            CALL LCMGET(IPMIC,HADF(ITYPE),FLXHET)
            IF(LADF) THEN
             IF (ADF_T(:3) .EQ. 'DRA') THEN
              DO I=1,NADF
              IF(HADF(ITYPE).EQ.ADFD(I))THEN
                DO IGR=1, NGP
                  ADF(I,IGR)= FLXHET(IGR)/FLXHOM(IGR)
                ENDDO
              ENDIF
              ENDDO
             ELSE IF ((ADF_T(:3) .EQ. 'GEN')) THEN
              IF(HADF(ITYPE).EQ.HFLX(1))THEN
               CALL LCMPUT(KPTH,'FLXL',NGP,2,FLXHET)
              ENDIF
              IF(HADF(ITYPE).EQ.HFLX(2))THEN
               CALL LCMPUT(KPTH,'FLXR',NGP,2,FLXHET)
              ENDIF
              IF (HADF(ITYPE).EQ.HCUR(1))THEN
               CALL LCMPUT(KPTH,'CURL',NGP,2,FLXHET)
              ENDIF
              IF (HADF(ITYPE).EQ.HCUR(2))THEN
               CALL LCMPUT(KPTH,'CURR',NGP,2,FLXHET)
              ENDIF
             ENDIF
            ENDIF
            IF(LCDF) THEN
              DO I=1,NCDF
              IF(HADF(ITYPE).EQ.CDFD(I))THEN
                DO IGR=1, NGP
                  CDF(I,IGR)= FLXHET(IGR)/FLXHOM(IGR)
                ENDDO
              ENDIF
              ENDDO
            ENDIF
          ENDDO
          IF(LADF) CALL LCMPUT(KPTH,'ADF',NADF*NGP,2,ADF)
          IF(LCDF) CALL LCMPUT(KPTH,'CDF',NCDF*NGP,2,CDF)
          IF(IPRINT>1) THEN
            WRITE(6,*)
            IF(LADF) WRITE(6,*)"ADF :",ADF
            IF(LCDF) WRITE(6,*)"CDF :",CDF
          ENDIF
        ENDIF
        FLXHOM(:)=FLXHOM(:) * VOLUME
      ENDIF
      IF(LGFF) THEN
        FLXHOM(:)=FLXHOM(:) / VOLUME
        CALL LCMSIX (IPDAT,' ',0)
        CALL LCMSIX (IPDAT,'SAPHYB_INFO',1)
        CALL LCMGTC(IPDAT,'GFF_TYPE',3,1,GFF_T)


        IF(GFF_T .EQ. 'DRA') THEN
          CALL LCMSIX(IPMIC,' ',0)
          CALL LCMSIX(IPMIC,'MACROLIB',1)
          CALL LCMSIX(IPMIC,'GFF',1)
          CALL LCMSIX(IPMIC,'GFF-GEOM',1)
          CALL LCMGET(IPMIC,'MIX',MIXG)
          CALL LCMSIX(IPMIC,'GFF-GEOM',2)
          CALL LCMLEN(IPMIC,'NWT0',ILONG,ITYLCM)
          IF (ILONG .NE. NGP*NGFF) THEN
           CALL XABORT("@D2PMAC : ERROR IN NUMBER OF GFF IN MCO")
          ENDIF
          CALL LCMGET(IPMIC,'NWT0',GFFC)
!          CALL LCMGET(IPMIC,'VOLUME',VOLG)
          CALL LCMGET(IPMIC,'H-FACTOR',KFC)
          DO J=1,NPIN
            DO I=1,NPIN
              DO IG=1,NGP
                GFF(I,J,IG)=GFFC(MIXG(I,J),IG)*KFC(MIXG(I,J),IG)
     >                      /FLXHOM(IG)/KAPPA_FI(IG)
              ENDDO
            ENDDO
          ENDDO

          IF(IPRINT>1) THEN
            WRITE(6,*)
            WRITE(6,*)"GFF                     :"
            DO IG=1,NGP
              WRITE(6,*)"Group                     :",IG
              DO J=1,NPIN
                WRITE(6,*)GFF(:,J,IG)
              ENDDO
            ENDDO
          ENDIF
          CALL LCMPUT(KPTH,'GFF',NPIN*NPIN*NGP,2,GFF)
        ENDIF
        FLXHOM(:)=FLXHOM(:) * VOLUME
      ENDIF

      FLUX(:)=FLXHOM(:)
      CALL LCMSIX(KPTH,' ',0)
      CALL LCMSIX(IPDAT,' ',0)

      IF(IPRINT>1) THEN
       WRITE(6,'(A)',advance="no")  "Energy group     :"
       DO I=1,NGP
        WRITE(6,'(5X,I12)',advance="no") I
       ENDDO
       WRITE(6,*)
       WRITE(6,'(A,8(5X,ES12.5E2))') "SIGWOO           :",SIGW00
       WRITE(6,'(A,8(5X,ES12.5E2))') "SIGSOO           :",SCAT
       WRITE(6,'(A,8(5X,ES12.5E2))') "TOTALE           :",XSECT
       WRITE(6,'(A,8(5X,ES12.5E2))') "DIFF             :",DIFF
       WRITE(6,'(A,8(5X,ES12.5E2))') "TRANSPORT        :",XTR
       WRITE(6,'(A,8(5X,ES12.5E2))') "ABSORPTION       :",ABSORPTION
       WRITE(6,'(A,8(5X,ES12.5E2))') "NU FISSION       :",X_NU_FI
       WRITE(6,'(A,8(5X,ES12.5E2))') "KAPPA FISSION    :",KAPPA_FI
       WRITE(6,'(A,8(5X,ES12.5E2))') "SCATTERING g->g' :"
       WRITE(6,'(8(5X,ES12.5E2))')SCAT_MAT
      ENDIF
      END
