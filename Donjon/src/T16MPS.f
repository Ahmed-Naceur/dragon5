*DECK T16MPS
      SUBROUTINE T16MPS(IPCPO ,IPRINT,MAXMIX,MNLOCP,MNCPLP,MNPERT,
     >                  NALOCP,IDLCPL,NCMIXS,NGCCPO,TITLE ,SUBTIT ,
     >                  ENECPO,NAMMIX,MIXRCI,PARRCI,PARPER)
*
*----
*
*Purpose:
*  Save mixture processing option on CPO.
*
*Author(s): 
* G. Marleau
*
*Parameters: input
* IPCPO   pointer to CPO data structure.
* IPRINT  print level.
* MAXMIX  maximum number of mixtures.
* MNLOCP  maximum number of local parameters.
* MNCPLP  maximum number of coupled parameters.
* MNPERT  maximum number of perturbations per local parameter.
* NALOCP  local parameter names allowed.
* IDLCPL  local ID for perturbation parameters.
* NCMIXS  number of current mixtures.
* NGCCPO  number of edit groups.
* TITLE   title.
* SUBTIT  subtitle.
* ENECPO  final energy group structure for CPO.
* NAMMIX  names of mixtures.
* MIXRCI  reference information for mixtures where:
*         =0 no information for mixture;
*         >0 information not updated; 
*         <0 information to be updated. 
* PARRCI  reference local parameters for mixtures.
* PARPER  perturbation parameters for mixtures.
*
*----
*
      USE GANLIB
      IMPLICIT         NONE
      TYPE(C_PTR)      IPCPO
      INTEGER          IPRINT,MAXMIX,MNLOCP,MNCPLP,MNPERT,
     >                 NCMIXS,NGCCPO
      CHARACTER        NALOCP(MNLOCP+MNCPLP)*4
      INTEGER          IDLCPL(2,MNLOCP+MNCPLP)
      CHARACTER        TITLE*72,SUBTIT*240
      INTEGER          NAMMIX(2,MAXMIX),
     >                 MIXRCI(2+MNLOCP+MNCPLP,MAXMIX)
      REAL             ENECPO(NGCCPO+1),PARRCI(MNLOCP,MAXMIX),
     >                 PARPER(MNPERT,2,MNLOCP+MNCPLP,MAXMIX)
*----
*  MEMORY ALLOCATION
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDSUF,NBPER
      REAL, ALLOCATABLE, DIMENSION(:) :: LOCALP
*----
*  LOCAL VARIABLES
*----
      INTEGER          IOUT,NPARAM,NTC,ILCMUP,ILCMDN
      CHARACTER        NAMSBR*6,NAMMAC*12
      PARAMETER       (IOUT=6,NPARAM=4,ILCMUP=1,ILCMDN=2,
     >                 NTC=3,NAMSBR='T16MPS',NAMMAC='MACR        ')
      CHARACTER        NAMDIR*12
      INTEGER          IPARAM(NPARAM),KCHAR(NTC),
     >                 NBURN,ILOCP,IMIX,IGC,
     >                 IPER,NBPERP,ILCMLN,ITYLCM,NLPAR,ILPAR,ILOCL,
     >                 IR,NPERTN,NMODRC
*----
*  SAVE MIXTURE NAMES
*----
      NAMDIR=NAMMAC
      READ(NAMDIR,'(3A4)') (KCHAR(IR),IR=1,NTC)
      NPERTN=MNLOCP+MNCPLP
      CALL LCMPUT(IPCPO,'MIXTURE-PREF',2*NCMIXS,3,NAMMIX)
*----
*  SAVE PERTURBATION SUFFIX NAMES
*----
      ALLOCATE(IDSUF(NPERTN))
      DO 100 ILOCP=1,NPERTN
        READ(NALOCP(ILOCP),'(A4)') IDSUF(ILOCP)
 100  CONTINUE
      CALL LCMPUT(IPCPO,'PERTURB-SUFX',NPERTN,3,IDSUF)
      DEALLOCATE(IDSUF)
*----
*  SAVE NUMBER OF PERTURBATION PER LOCAL PARAMETER PER MIXTURE
*----
      ALLOCATE(NBPER(NPERTN*NCMIXS))
      IPER=0
      DO IMIX=1,NCMIXS
        DO ILOCP=1,NPERTN
          IPER=IPER+1
          NBPER(IPER)=ABS(MIXRCI(2+ILOCP,IMIX))
        ENDDO
      ENDDO
      CALL LCMPUT(IPCPO,'PERTURB-NUMB',NPERTN*NCMIXS,
     >            1,NBPER)
      DEALLOCATE(NBPER)
      ALLOCATE(LOCALP(MNLOCP))
*----
*  SCAN OVER MIXTURES
*----
      DO IMIX=1,NCMIXS
*----
*  MIXTURE TO UPDATE
*----
        WRITE(NAMDIR,'(A4,A2,A6)')
     >  NAMMIX(1,IMIX),NAMMIX(2,IMIX),'RC    '
        NBURN=ABS(MIXRCI(2,IMIX))
        IPARAM(1)=NGCCPO
        IPARAM(2)=1
        IPARAM(3)=2
        IPARAM(4)=NBURN
        NMODRC=0
        DO ILOCP=1,NPERTN
          IF(MIXRCI(2+ILOCP,IMIX) .LT. 0) THEN
            NMODRC=NMODRC-1
          ENDIF
        ENDDO
        IF(MIXRCI(2,IMIX) .LT. 0 ) THEN
          CALL LCMSIX(IPCPO,NAMDIR,ILCMUP)
          CALL LCMLEN(IPCPO,'PARAM       ',ILCMLN,ITYLCM)
          IF(ILCMLN .EQ. 0) THEN
            CALL LCMPUT(IPCPO,'PARAM       ',NPARAM,1,IPARAM)
            CALL LCMPUT(IPCPO,'ENERGY      ',NGCCPO+1,2,ENECPO)
          ELSE
            CALL LCMGET(IPCPO,'PARAM       ',IPARAM)
            IF(IPARAM(1) .NE. NGCCPO .OR.
     >         IPARAM(2) .NE. 1      .OR.
     >         IPARAM(3) .NE. 2      .OR.
     >         IPARAM(4) .NE. NBURN      ) THEN
*----
*  ABORT SINCE REFERENCE CASE PARAMETERS HAVE CHANGED
*----
              CALL XABORT(NAMSBR//
     >        ': INCOMPATIBLE PARAMETERS FOR '//NAMDIR)
            ENDIF
          ENDIF
          CALL XDRSET(LOCALP,MNLOCP,0.0)
          DO ILOCP=1,MNLOCP
            LOCALP(ILOCP)=PARRCI(ILOCP,IMIX)
          ENDDO
          CALL LCMPUT(IPCPO,'LOCAL-PARAMS',MNLOCP, 2,LOCALP)
          CALL LCMPTC(IPCPO,'TITLE       ',    72, 1,        TITLE)
          CALL LCMPTC(IPCPO,'SUB-TITLE   ',   240, 1,       SUBTIT)
          CALL LCMPUT(IPCPO,'ISOTOPESNAME',   NTC, 3,        KCHAR)
          CALL LCMSIX(IPCPO,NAMDIR,ILCMDN)
        ELSE IF(NMODRC .LT. 0) THEN
          CALL LCMSIX(IPCPO,NAMDIR,ILCMUP)
          CALL XDRSET(LOCALP,MNLOCP,0.0)
          DO ILOCP=1,MNLOCP
            LOCALP(ILOCP)=PARRCI(ILOCP,IMIX)
          ENDDO
          CALL LCMPUT(IPCPO,'LOCAL-PARAMS',MNLOCP, 2,LOCALP)
          CALL LCMSIX(IPCPO,NAMDIR,ILCMDN)
        ENDIF
*----
*  PERTURBATIONS TO UPDATE
*----
        DO ILOCP=1,NPERTN
          NBPERP=ABS(MIXRCI(2+ILOCP,IMIX))
          IF(MIXRCI(2+ILOCP,IMIX) .LT. 0) THEN
            NLPAR=1
            IF(ILOCP .GT. MNLOCP) NLPAR=2
            DO IPER=1,NBPERP
              WRITE(NAMDIR,'(A4,A2,A4,I2)')
     >        NAMMIX(1,IMIX),NAMMIX(2,IMIX),NALOCP(ILOCP),IPER
              CALL LCMSIX(IPCPO,NAMDIR,ILCMUP)
              CALL LCMLEN(IPCPO,'PARAM       ',ILCMLN,ITYLCM)
              IF(ILCMLN .EQ. 0) THEN
                CALL LCMPUT(IPCPO,'PARAM       ',NPARAM,1,IPARAM)
                CALL LCMPUT(IPCPO,'ENERGY      ',NGCCPO+1,2,ENECPO)
              ELSE
                CALL LCMGET(IPCPO,'PARAM       ',IPARAM)
                IF(IPARAM(1) .NE. NGCCPO .OR.
     >             IPARAM(2) .NE. 1      .OR.
     >             IPARAM(3) .NE. 2      .OR.
     >             IPARAM(4) .NE. NBURN      ) THEN
*----
*  ABORT SINCE PERTURBATION PARAMETERS HAVE CHANGED
*----
                  CALL XABORT(NAMSBR//
     >            ': INCOMPATIBLE PARAMETERS FOR '//NAMDIR)
                ENDIF
              ENDIF
*----
*  TRANSFER REFERENCE PARAMETERS
*----
              DO ILOCL=1,MNLOCP
                LOCALP(ILOCL)=PARRCI(ILOCL,IMIX)
              ENDDO
*----
*  TRANSFER PERTURBED PARAMETERS
*----
              DO ILPAR=1,NLPAR
                ILOCL=IDLCPL(ILPAR,ILOCP)
                LOCALP(ILOCL)=PARPER(IPER,ILPAR,ILOCP,IMIX)
              ENDDO
              CALL LCMPUT(IPCPO,'LOCAL-PARAMS',MNLOCP, 2,       LOCALP)
              CALL LCMPTC(IPCPO,'TITLE       ',    72, 1,        TITLE)
              CALL LCMPTC(IPCPO,'SUB-TITLE   ',   240, 1,       SUBTIT)
              CALL LCMPUT(IPCPO,'ISOTOPESNAME',   NTC, 3,        KCHAR)
              CALL LCMSIX(IPCPO,NAMDIR,ILCMDN)
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      DEALLOCATE(LOCALP)
*----
*  PROCESS PRINT LEVEL
*----
      IF(IPRINT .GE. 1) THEN
        WRITE(IOUT,6000) NAMSBR,NCMIXS,NGCCPO
        WRITE(IOUT,6030) (ENECPO(IGC),IGC=1,NGCCPO+1)
        DO IMIX=1,NCMIXS
          IF(MIXRCI(2,IMIX) .LT. 0) THEN
            WRITE(IOUT,6010) (NAMMIX(IR,IMIX),IR=1,2),
     >                        ABS(MIXRCI(2,IMIX))
            WRITE(IOUT,6020)
     >        (NALOCP(ILOCP),PARRCI(ILOCP,IMIX),ILOCP=1,MNLOCP)
          ELSE IF(MIXRCI(2,IMIX) .GT. 0) THEN
            WRITE(IOUT,6011) (NAMMIX(IR,IMIX),IR=1,2),
     >                        ABS(MIXRCI(2,IMIX))
            WRITE(IOUT,6020)
     >      (NALOCP(ILOCP),PARRCI(ILOCP,IMIX),ILOCP=1,MNLOCP)
          ENDIF
          DO ILOCP=1,NPERTN
            NBPERP=ABS(MIXRCI(2+ILOCP,IMIX))
            NLPAR=1
            IF(ILOCP .GT. MNLOCP) NLPAR=2
            IF(MIXRCI(2+ILOCP,IMIX) .LT. 0) THEN
              WRITE(IOUT,6012) NBPERP,NALOCP(ILOCP)
              DO IPER=1,NBPERP
                WRITE(IOUT,6021) IPER,
     >          (NALOCP(IDLCPL(ILPAR,ILOCP)),
     >           PARPER(IPER,ILPAR,ILOCP,IMIX),
     >           ILPAR=1,NLPAR)
              ENDDO
            ELSE IF(MIXRCI(2+ILOCP,IMIX) .GT. 0) THEN
              WRITE(IOUT,6013) NBPERP,NALOCP(ILOCP)
              DO IPER=1,NBPERP
                WRITE(IOUT,6021) IPER,
     >          (NALOCP(IDLCPL(ILPAR,ILOCP)),
     >           PARPER(IPER,ILPAR,ILOCP,IMIX),
     >           ILPAR=1,NLPAR)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
        WRITE(IOUT,6001)
      ENDIF
*----
*  PRINT FORMAT
*----
 6000 FORMAT(1X,5('*'),'OUTPUT FROM ',A6,5('*')/
     >       6X,'CONTENTS OF CPO AFTER UPDATE'/
     >       6X,'NUMBER OF MIXTURES     = ',I10/
     >       6X,'NUMBER OF GROUPS       = ',I10)
 6001 FORMAT(1X,28('*'))
 6010 FORMAT(6X,'NAME OF MIXTURES    = ',A4,A4,
     >       2X,'NUMBER OF BURNUP    = ',I10,
     >       2X,'UPDATED FROM CURRENT TAPE16')
 6011 FORMAT(6X,'NAME OF MIXTURES    = ',A4,A4,
     >       2X,'NUMBER OF BURNUP    = ',I10,
     >       2X,'TAKEN FROM OLD CPO')
 6012 FORMAT(6X,I10,' PERTURBATIONS FOR ',A4,
     >       2X,'UPDATED FROM CURRENT TAPE16')
 6013 FORMAT(6X,I10,' PERTURBATIONS FOR ',A4,
     >       2X,'UPDATED FROM OLD CPO')
 6020 FORMAT(6X,'LOCAL PARAMETER FOR REFERENCE CASE'/
     >       1P,6(2X,A4,1X,E10.3))
 6021 FORMAT(6X,'LOCAL PARAMETER PERTURBATION = ',I2,
     >       1P,2(2x,A4,1X,E10.3))
 6030 FORMAT(6X,'CPO ENERGY STRUCTURE (EV)'/
     >1P,10(2X,E10.3))
      RETURN
      END
