*DECK T16WDS
      SUBROUTINE T16WDS(IPCPO ,NGCCPO,NVXSR ,NMXSR ,IBURN ,EFJ   ,
     >                  NAMDXS,ITYXS ,FLXINT,FLXDIS,OVERV ,RECXSV,
     >                  IDRXSM,RECXSM,RECSCA)
*
*----
*
*Purpose:
*  Write properties to CPO data structure.
*
*Author(s):
* G. Marleau
*
*Parameters: input
* IPCPO   pointer to CPO data structure.
* NGCCPO  number of edit groups.
* NVXSR   number of vector cross sections.
* NMXSR   number of matrix cross sections.
* IBURN   burnup step.
* EFJ     energy of fission in joules.
* NAMDXS  name of vector cross sections.
* ITYXS   types of cross sections saved.
* FLXINT  volume integrated fluxes.
* FLXDIS  flux disadvantage factor.
* OVERV   1/V cross sections.
* RECXSV  vector cross sections.
* IDRXSM  compression vector for matrix cross sections.
* RECXSM  matrix cross sections.
* RECSCA  dummy matrix cross sections.
*
*----
*
      USE GANLIB
      IMPLICIT         NONE
      TYPE(C_PTR)      IPCPO
      INTEGER          NGCCPO,NVXSR,NMXSR,IBURN
      CHARACTER        NAMDXS(NVXSR+NMXSR)*12
      INTEGER          IDRXSM(NGCCPO,2),ITYXS(NVXSR+NMXSR)
      REAL             EFJ,FLXINT(NGCCPO),
     >                 FLXDIS(NGCCPO),OVERV(NGCCPO),
     >                 RECXSV(NGCCPO,NVXSR+NMXSR),
     >                 RECXSM(NGCCPO,NGCCPO,NMXSR),
     >                 RECSCA(NGCCPO*NGCCPO)
*----
*  LOCAL VARIABLES
*----
      INTEGER          IOUT,ILCMUP,ILCMDN
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,ILCMUP=1,ILCMDN=2,NAMSBR='T16WDS')
      CHARACTER        NAMBRN*12,NAMMAC*12,NAMLEG*2
      INTEGER          IVXS,IMXS,IGTO,IGFROM,IGMIN,IGMAX,NXSCMP
      REAL             DENMAC
*----
*  SET UP BURUP DIRECTORY
*----
      WRITE(NAMBRN,'(A8,I4)') 'BURN    ',IBURN
      CALL LCMSIX(IPCPO ,NAMBRN,ILCMUP)
*----
*  SAVE ISOTOPES DENSITY, ENERGY, INTEGRATED FLUX,
*  DISADVANTAGE FACTOR AND OVERV ON MAIN DIRECTORY
*----
      DENMAC=1.0
      CALL LCMPUT(IPCPO ,'ISOTOPESDENS',     1,2,DENMAC)
      CALL LCMPUT(IPCPO ,'ISOTOPES-EFJ',     1,2,EFJ)
      CALL LCMPUT(IPCPO ,'FLUX-INTG   ',NGCCPO,2,FLXINT)
      CALL LCMPUT(IPCPO ,'FLUXDISAFACT',NGCCPO,2,FLXDIS)
      CALL LCMPUT(IPCPO ,'OVERV       ',NGCCPO,2,OVERV)
      NAMMAC='MACR        '
      CALL LCMSIX(IPCPO ,NAMMAC,ILCMUP)
*----
*  FIND IF VECTOR XS NOT ALL 0.0
*  AND INITIALIZE ITYXS ACCORDINGLY
*  SAVE XS
*----
      DO IVXS=1,NVXSR
        ITYXS(IVXS)=0
        DO IGFROM=1,NGCCPO
          IF(RECXSV(IGFROM,IVXS) .NE. 0.0) THEN
            ITYXS(IVXS)=1
            CALL LCMPUT(IPCPO ,NAMDXS(IVXS),
     >                  NGCCPO,2,RECXSV(1,IVXS))
          ENDIF
        ENDDO
      ENDDO
*----
*  FIND IF SCATTERING XS NOT ALL 0.0
*  AND INITIALIZE ITYXS ACCORDINGLY
*----
      DO IMXS=1,NMXSR
        ITYXS(IMXS+NVXSR)=0
        DO IGTO=1,NGCCPO
          DO IGFROM=1,NGCCPO
            IF(RECXSM(IGTO,IGFROM,IMXS) .NE. 0.0) THEN
              ITYXS(IMXS+NVXSR)=1
              CALL LCMPUT(IPCPO ,NAMDXS(IMXS+NVXSR),
     >                    NGCCPO,2,RECXSV(1,IMXS+NVXSR))
              GO TO 105
            ENDIF
          ENDDO
        ENDDO
 105    CONTINUE
      ENDDO
*----
*  SAVE ITYXS
*----
      CALL LCMPUT(IPCPO ,'XS-SAVED    ',NVXSR+NMXSR,1,ITYXS)
*----
*  COMPRESS SCATTERING MATRIX
*  RECXSM(IGTO,IGFROM,IMXS) REPRESENT SCATTERING CROSS SECTION
*    FROM GROUP "IGFROM" TO GROUP "IGTO"
*  IDRXSM(IGTO,1) IS MAXIMUM GROUP NUMBER
*    WITH SCATTERING TO "IGTO" GROUP
*  IDRXSM(IGTO,2) IS NUMBER OF GROUPS
*    WITH SCATTERING TO "IGTO" GROUP
*  RECSCA(IX) IS COMPRESSED SCATTERING MATRIX
*  IX CAN BE LOCALIZED IN RECXSM(IGTO,IGFROM) USING
*  IF(IGTO=1) THEN
*    IPOSD=1
*  ELSE
*    IPOSD=1+SUM( IDRXSM(IGF,2) , IGF=1,IGTO-1)
*  ENDIF
*  IF(IGFROM.GT.IDRXSM(IGTO,1)) THEN
*    XSSCMP NOT STORED
*  ELSE IF(IGFROM.LT.IDRXSM(IGTO,1)-IDRXSM(IGTO,2)+1) THEN
*    XSSCMP NOT STORED
*  ELSE
*    IX=IPOSD+IDRXSM(IGTO,1)-IGFROM
*    RECSCA(IX)=RECXSM(IGTO,IGFROM)
*  ENDIF
*----
      DO IMXS=1,NMXSR
        NXSCMP=0
        DO IGTO=1,NGCCPO
          IGMIN=IGTO
          IGMAX=IGTO
          DO IGFROM=1,NGCCPO
            IF(RECXSM(IGTO,IGFROM,IMXS) .NE. 0.0) THEN
              IGMIN=MIN(IGMIN,IGFROM)
              IGMAX=MAX(IGMAX,IGFROM)
            ENDIF
          ENDDO
          IDRXSM(IGTO,1)=IGMAX
          IDRXSM(IGTO,2)=IGMAX-IGMIN+1
          DO IGFROM=IGMAX,IGMIN,-1
            NXSCMP=NXSCMP+1
            RECSCA(NXSCMP)=RECXSM(IGTO,IGFROM,IMXS)
          ENDDO   
        ENDDO
        WRITE(NAMLEG,'(I2)') IMXS-1
        CALL LCMPUT(IPCPO,'NJJ '//NAMLEG//'      ',NGCCPO,1,IDRXSM(1,1))
        CALL LCMPUT(IPCPO,'IJJ '//NAMLEG//'      ',NGCCPO,1,IDRXSM(1,2))
        CALL LCMPUT(IPCPO,'SCAT'//NAMLEG//'      ',NXSCMP,2,RECSCA)
      ENDDO
      CALL LCMSIX(IPCPO ,NAMMAC,ILCMDN)
      CALL LCMSIX(IPCPO ,NAMBRN,ILCMDN)
      RETURN
      END
