*DECK T16MPI
      SUBROUTINE T16MPI(IPCPO ,IPRINT,MAXMIX,MNLOCP,MNCPLP,MNPERT,
     >                  NALOCP,IDLCPL,NCMIXS,NGCCPO,ENECPO,NAMMIX,
     >                  MIXRCI,PARRCI,PARPER)
*
*----
*
*Purpose:
*  Initialize mixture processing option.
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
* ENECPO  final energy group structure for CPO.
* NAMMIX  names of mixtures.
*
*Parameters: output
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
      INTEGER          IDLCPL(2,MNLOCP+MNCPLP),NAMMIX(2,MAXMIX),
     >                 MIXRCI(2+MNLOCP+MNCPLP,MAXMIX)
      REAL             ENECPO(NGCCPO+1),PARRCI(MNLOCP,MAXMIX),
     >                 PARPER(MNPERT,2,MNLOCP+MNCPLP,MAXMIX)
*----
*  MEMORY ALLOCATION
*----
      INTEGER ,ALLOCATABLE, DIMENSION(:) :: IDSUF,NBPER
      REAL,ALLOCATABLE, DIMENSION(:) :: ENEMIX,LOCALP
*----
*  LOCAL VARIABLES
*----
      INTEGER          IOUT,NPARAM,NTC,ILCMUP,ILCMDN
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NPARAM=4,ILCMUP=1,ILCMDN=2,
     >                 NTC=3,NAMSBR='T16MPI')
      CHARACTER        NAMDIR*12
      INTEGER          IPARAM(NPARAM),
     >                 ILOCP,IMIX,IPER,NBPERP,IGR,
     >                 NLPAR,ILPAR,ILOCL,NPERTN,IR
*----
*  GET MIXTURE NAMES
*----
      CALL LCMGET(IPCPO,'MIXTURE-PREF',NAMMIX)
*----
*  GET PERTURBATION SUFFIX NAMES AND COMPARE WITH
*  REFERENCE NAMES
*----
      NPERTN=MNLOCP+MNCPLP
      ALLOCATE(ENEMIX(NGCCPO+1),IDSUF(NPERTN))
      CALL LCMGET(IPCPO,'PERTURB-SUFX',IDSUF)
      DO ILOCP=1,NPERTN
        WRITE(NAMDIR,'(A4)') IDSUF(ILOCP)
        IF(NAMDIR .NE. NALOCP(ILOCP)) CALL XABORT(NAMSBR//
     >  ': INCOHERENT PERTURBATION NAMES')
      ENDDO
      DEALLOCATE(IDSUF)
*----
*  GET NUMBER OF PERTURBATION PER LOCAL PARAMETER PER MIXTURE
*----
      ALLOCATE(NBPER(NPERTN*NCMIXS))
      CALL LCMGET(IPCPO,'PERTURB-NUMB',NBPER)
      IPER=0
      DO IMIX=1,NCMIXS
        DO ILOCP=1,NPERTN
          IPER=IPER+1
          MIXRCI(2+ILOCP,IMIX)=NBPER(IPER)
        ENDDO
      ENDDO
      DEALLOCATE(NBPER)
      ALLOCATE(LOCALP(MNLOCP))
*----
*  GET LOCAL PARAMETERS
*----
      DO IMIX=1,NCMIXS
*----
*  GET REFERENCE MIXTURE PARAMETERS
*----
        WRITE(NAMDIR,'(A4,A2,A6)')
     >  NAMMIX(1,IMIX),NAMMIX(2,IMIX),'RC    '
        CALL LCMSIX(IPCPO,NAMDIR,ILCMUP)
        CALL LCMGET(IPCPO,'PARAM       ',IPARAM)
        IF(IPARAM(1) .NE. NGCCPO) THEN
          CALL XABORT(NAMSBR//
     >    ': INVALID NUMBER OF GROUPS FOR REFERENCE')
        ENDIF
*----
*  TEST ENERGY GROUP STRUCTURE
*----
        CALL LCMGET(IPCPO,'ENERGY      ',ENEMIX)
        DO IGR=1,NGCCPO+1
          IF(ENECPO(IGR) .NE. ENEMIX(IGR)) THEN
            CALL XABORT(NAMSBR//
     >      ': INVALID GROUP STRUCTURE FOR REFERENCE')
          ENDIF
        ENDDO
*----
*  READ LOCAL PARAMETERS AND TRANSFER
*  TO LOCAL TABLE
*----
        MIXRCI(2,IMIX)=IPARAM(4)
        CALL LCMGET(IPCPO,'LOCAL-PARAMS',LOCALP)
        DO ILOCP=1,MNLOCP
          PARRCI(ILOCP,IMIX)=LOCALP(ILOCP)
        ENDDO
        CALL LCMSIX(IPCPO,NAMDIR,ILCMDN)
*----
*  GET PERTURBATIONS PARAMETERS
*----
        DO ILOCP=1,NPERTN
          NBPERP=MIXRCI(2+ILOCP,MAXMIX)
          IF(NBPERP .GT. 0) THEN
            NLPAR=1
            IF(ILOCP .GT. MNLOCP) NLPAR=2
            DO IPER=1,NBPERP
              WRITE(NAMDIR,'(A4,A2,A4,I2)')
     >        NAMMIX(1,IMIX),NAMMIX(2,IMIX),NALOCP(ILOCP),IPER
              CALL LCMSIX(IPCPO,NAMDIR,ILCMUP)
              CALL LCMGET(IPCPO,'PARAM       ',IPARAM)
              IF(IPARAM(1) .NE. NGCCPO) THEN
                CALL XABORT(NAMSBR//
     >          ': INVALID NUMBER OF GROUPS FOR PERTURBATION')
              ENDIF
*----
*  TEST ENERGY GROUP STRUCTURE
*----
              CALL LCMGET(IPCPO,'ENERGY      ',ENEMIX)
              DO IGR=1,NGCCPO+1
                IF(ENECPO(IGR) .NE. ENEMIX(IGR)) THEN
                  CALL XABORT(NAMSBR//
     >            ': INVALID GROUP STRUCTURE FOR REFERENCE')
                ENDIF
              ENDDO
              CALL LCMGET(IPCPO,'LOCAL-PARAMS',LOCALP)
*----
*  READ PERTURBED LOCAL PARAMETERS AND TRANSFER
*  TO LOCAL TABLE
*----
              DO ILPAR=1,NLPAR
                ILOCL=IDLCPL(ILPAR,ILOCP)
                PARPER(IPER,ILPAR,ILOCP,IMIX)=LOCALP(ILOCL)
              ENDDO
              CALL LCMSIX(IPCPO,NAMDIR,ILCMDN)
            ENDDO
          ENDIF
        ENDDO
      ENDDO
      DEALLOCATE(LOCALP,ENEMIX)
*----
*  PROCESS PRINT LEVEL
*----
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6000) NAMSBR,NCMIXS
        DO IMIX=1,NCMIXS
          IF(MIXRCI(2,IMIX) .GT. 0) THEN
            WRITE(IOUT,6010) (NAMMIX(IR,IMIX),IR=1,2),
     >                        MIXRCI(2,IMIX)
            WRITE(IOUT,6020)
     >      (NALOCP(ILOCP),PARRCI(ILOCP,IMIX),ILOCP=1,MNLOCP)
          ENDIF
          DO ILOCP=1,NPERTN
            NBPERP=MIXRCI(2+ILOCP,IMIX)
            NLPAR=1
            IF(ILOCP .GT. MNLOCP) NLPAR=2
            IF(NBPERP .GT. 0) THEN
              WRITE(IOUT,6011) NBPERP,NALOCP(ILOCP)
              DO IPER=1,NBPERP
                WRITE(IOUT,6022) IPER,
     >          (NALOCP(IDLCPL(ILPAR,ILOCP)),
     >           PARPER(IPER,ILPAR,ILOCP,IMIX),
     >           ILPAR=1,NLPAR)
              ENDDO
            ENDIF
          ENDDO
        ENDDO
        WRITE(IOUT,6001)
      ENDIF
      RETURN
*----
*  PRINT FORMAT
*----
 6000 FORMAT(1X,5('*'),'OUTPUT FROM ',A6,5('*')/
     >       6X,'CONTENTS OF CPO BEFORE UPDATE'/
     >       6X,'NUMBER OF MIXTURES  = ',I10)
 6001 FORMAT(1X,28('*'))
 6010 FORMAT(6X,'NAME OF MIXTURES    = ',A4,A4,
     >       2X,'NUMBER OF BURNUP    = ',I10,
     >       2X,'ALREADY STORED ON CPO FILE')
 6011 FORMAT(6X,I10,' PERTURBATIONS FOR ',A4,
     >       2X,'ALREADY STORED ON CPO FILE')
 6020 FORMAT(6X,'LOCAL PARAMETER FOR REFERENCE CASE'/
     >       1P,6(2x,A4,1X,E10.3))
 6022 FORMAT(6X,'LOCAL PARAMETER PERTURBATION = ',I2/
     >       1P,2(2x,A4,1X,E10.3))
      END
