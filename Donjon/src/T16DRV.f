*DECK T16DRV
      SUBROUTINE T16DRV(IPCPO ,IFT16 ,IPRINT,MNLOCP,MNCPLP,MNPERT,
     >                  NALOCP,NCMIXS,NGCCPO,MNBURN,NG    ,NGMTR ,
     >                  NMATZ ,MTRMSH,NZONE ,IFGMTR,VELMTR,NAMMIX,
     >                  MIXRCI,MIXPER,MIXREG)
*
*----
*
*Purpose:
*  Main driver for the transfer of cross sections from tape16 to CPO.
*
*Author(s):
* G. Marleau
*
*Parameters: input
* IPCPO   pointer to CPO data structure.
* IFT16   tape16 file unit.
* IPRINT  print level where:
*         =0 for no print; >=  1 print processing option.
* MNLOCP  maximum number of local parameters.
* MNCPLP  maximum number of coupled parameters.
* MNPERT  maximum number of perturbations per local parameter.
* NALOCP  local parameter names allowed.
* NCMIXS  number of current mixtures.
* NGCCPO  number of edit groups.
* MNBURN  maximum number or burnup steps.
* NG      number of groups on X-S library.
* NGMTR   number of main transport groups.
* NMATZ   number of mixtures.
* MTRMSH  number of main transport mesh points.
* NZONE   number of zones.
* IFGMTR  fewgroups for main transport.
* VELMTR  velocity for main transport.
* NAMMIX  names of mixtures.
* MIXRCI  reference information for mixtures where:
*         =0 no information for mixture;
*         >0 information not updated; 
*         <0 information to be updated. 
* MIXPER  perturbation information for mixtures.
*         =0 no information for mixture;
*         >0 information not updated; 
*         <0 information to be updated. 
* MIXREG  mixture update identifier where:
*          =0 do not update;
*          =-1 update using CELLAV information;
*          > 0 update using specified region number.
*
*----
*
      USE GANLIB
      IMPLICIT         NONE
      TYPE(C_PTR)      IPCPO
      INTEGER          IFT16,IPRINT,MNLOCP,MNCPLP,MNPERT,
     >                 NCMIXS,NGCCPO,MNBURN,NG,NGMTR,NMATZ,
     >                 MTRMSH,NZONE
      CHARACTER        NALOCP(MNLOCP+MNCPLP)*4
      INTEGER          IFGMTR(NGCCPO),NAMMIX(2,NCMIXS),
     >                 MIXRCI(2+MNLOCP+MNCPLP,NCMIXS),
     >                 MIXPER(MNPERT,MNLOCP+MNCPLP,NCMIXS),
     >                 MIXREG(NCMIXS)
      REAL             VELMTR(NGMTR)
*----
*  MEMORY ALLOCATION
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: BURNUP,WNKB
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KMSPEC,MATMSH,IDRXSM
      REAL, ALLOCATABLE, DIMENSION(:) :: VQLE,FLXINT,FLXDIS,
     >                                   OVERV,RECXSV,RECXSM,
     >                                   RECTMP,RECSCA,ZONVOL,ZONRAD
*----
*  T16 PARAMETERS
*----
      INTEGER          MAXKEY
      PARAMETER       (MAXKEY=2)
      CHARACTER        TKEY1(MAXKEY)*10,TKEY2(MAXKEY)*10,
     >                 RKEY1*10,RKEY2*10
      INTEGER          NKEY,IOPT,NBE,IR,ZONNUM
*----
*  LOCAL VARIABLES
*----
      INTEGER          IOUT,ILCMUP,ILCMDN,NVXSR,NMXSR
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,ILCMUP=1,ILCMDN=2,NVXSR=20,NMXSR=2,
     >                 NAMSBR='T16DRV')
      CHARACTER        NAMDIR*12
      INTEGER          IMIX,ILOCP,NBPERP,IPER,NBURN,IBURN,NPERTN,
     >                 INEXTR,ITYXS(NVXSR+NMXSR),IMIREG,
     >                 MMXM
      REAL             VOLUME,B2CRI(3),BRNIRR(3),EFJ
*----
*  DATA
*----
      CHARACTER        NAMDXS(NVXSR+NMXSR)*12
      SAVE             NAMDXS
      DATA             NAMDXS
     > /'TOTAL       ','TRANC       ','NUSIGF      ','NFTOT       ',
     >  'CHI         ','NU          ','NG          ','NHEAT       ',
     >  'N2N         ','N3N         ','N4N         ','NP          ',
     >  'NA          ','GOLD        ','ABS         ','NWT0        ',
     >  'STRD        ','STRD X      ','STRD Y      ','STRD Z      ',
     >  'SIGS 0      ','SIGS 1      '/
*----
*  ALLOCATE MEMORY
*----
      MMXM=MAX(NZONE,MTRMSH)
      NPERTN=MNLOCP+MNCPLP
      ALLOCATE(BURNUP(MNBURN),WNKB(MNBURN))
      ALLOCATE(KMSPEC(NMATZ),MATMSH(MMXM),IDRXSM(NGCCPO*2))
      ALLOCATE(VQLE(MMXM),FLXINT(NGCCPO),FLXDIS(NGCCPO),
     >         OVERV(NGCCPO),RECXSV(NGCCPO*(NVXSR+NMXSR)),
     >         RECXSM(NGCCPO*NGCCPO*NMXSR),RECTMP(4*NGMTR),
     >         RECSCA(NGMTR*NGMTR))
      ALLOCATE(ZONVOL(NZONE),ZONRAD(NZONE))
*----
*  FIND MATERIAL SPECTRUM
*  REQUIRED FOR FLUX DISADVANTAGE FACTOR
*----     
      IOPT=1
      TKEY1(1)='MATERIAL  '
      TKEY2(1)='SPECTRUM  '
      NKEY=1
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >            NBE   )
      IF( NBE .NE. NMATZ ) THEN
        WRITE(6,'(128A1)') 'PLEASE RE-RUN WIMS-AECL BECAUSE '//
     > 'T16CPO UTILITY NEEDS A RECORD: '//TKEY1(1)//TKEY1(2)
        CALL XABORT(NAMSBR//
     >  ': CANNOT FIND '//TKEY1(1)//' '//TKEY2(1))
      ELSE
        READ(IFT16) RKEY1,RKEY2,NBE,(KMSPEC(IR),IR=1,NMATZ)
      ENDIF
*----
*  SCAN OVER MIXTURES
*----
      DO IMIX=1,NCMIXS
*----
*  MIXTURE TO UPDATE
*----
        CALL XDRSET(BRNIRR,3,0.0)
        IMIREG=MIXREG(IMIX)
        WRITE(NAMDIR,'(A4,A2,A6)') 
     >  NAMMIX(1,IMIX),NAMMIX(2,IMIX),'RC    '
        NBURN=ABS(MIXRCI(2,IMIX))
        IF(MIXRCI(2,IMIX) .LT. 0) THEN
          CALL LCMSIX(IPCPO,NAMDIR,ILCMUP)
          INEXTR=MIXRCI(1,IMIX)
          DO IBURN=1,NBURN
*----
*  BURNUP STEP TO UPDATE
*----
            CALL T16REC(IFT16 ,IPRINT,INEXTR)
            CALL T16FLX(IFT16 ,IPRINT,NGCCPO,NG    ,NGMTR ,NMATZ ,
     >                  MMXM  ,MTRMSH,NZONE ,IFGMTR,VELMTR,IMIREG,
     >                  VOLUME,B2CRI ,FLXINT,FLXDIS,OVERV ,KMSPEC,
     >                  MATMSH,VQLE  ,ZONNUM, ZONRAD,ZONVOL)
            IF(IMIREG .GT. 0) THEN
              CALL T16REC(IFT16 ,IPRINT,INEXTR)
              CALL T16RRE(IFT16 ,IPRINT,NGCCPO,NGMTR ,IFGMTR,NVXSR ,
     >                    NMXSR ,IMIREG,VELMTR,B2CRI ,BRNIRR,FLXINT,
     >                    OVERV ,RECXSV,RECXSM,RECTMP,RECSCA)
            ELSE
              CALL T16REC(IFT16 ,IPRINT,INEXTR)
              CALL T16RCA(IFT16 ,IPRINT,NGCCPO,NGMTR ,IFGMTR,NVXSR ,
     >                    NMXSR ,B2CRI ,BRNIRR,NZONE ,RECXSV,RECXSM,
     >                    RECTMP,RECSCA,ZONVOL)
            ENDIF
            BURNUP(IBURN)=BRNIRR(1)
            WNKB(IBURN)=BRNIRR(2)
            EFJ=BRNIRR(3)
            CALL T16WDS(IPCPO ,NGCCPO,NVXSR ,NMXSR ,IBURN ,EFJ   ,
     >                  NAMDXS,ITYXS ,FLXINT,FLXDIS,OVERV ,RECXSV,
     >                  IDRXSM,RECXSM,RECSCA)
            INEXTR=INEXTR+1
          ENDDO
          CALL LCMPUT(IPCPO ,'BURNUP      ',NBURN,2,BURNUP)
          CALL LCMPUT(IPCPO ,'N/KB        ',NBURN,2,WNKB)
          CALL LCMSIX(IPCPO,NAMDIR,ILCMDN)
        ENDIF
*----
*  PERTURBATIONS TO UPDATE
*----
        DO ILOCP=1,NPERTN
          NBPERP=ABS(MIXRCI(2+ILOCP,IMIX))
          IF(MIXRCI(2+ILOCP,IMIX) .LT. 0) THEN
            DO IPER=1,NBPERP
              INEXTR=MIXPER(IPER,ILOCP,IMIX)
              WRITE(NAMDIR,'(A4,A2,A4,I2)')
     >        NAMMIX(1,IMIX),NAMMIX(2,IMIX),NALOCP(ILOCP),IPER
              CALL LCMSIX(IPCPO,NAMDIR,ILCMUP)
              DO IBURN=1,NBURN
                 CALL T16REC(IFT16 ,IPRINT,INEXTR)
                 CALL T16FLX(IFT16 ,IPRINT,NGCCPO,NG    ,NGMTR ,NMATZ ,
     >                       MMXM  ,MTRMSH,NZONE ,IFGMTR,VELMTR,IMIREG,
     >                       VOLUME,B2CRI ,FLXINT,FLXDIS,OVERV ,KMSPEC,
     >                       MATMSH,VQLE  ,ZONNUM, ZONRAD,ZONVOL)
                IF(IMIREG .GT. 0) THEN
                  CALL T16REC(IFT16 ,IPRINT,INEXTR)
                  CALL T16RRE(IFT16 ,IPRINT,NGCCPO,NGMTR ,IFGMTR,NVXSR ,
     >                        NMXSR ,IMIREG,VELMTR,B2CRI ,BRNIRR,FLXINT,
     >                        OVERV,RECXSV,RECXSM,RECTMP,RECSCA)
                ELSE
                  CALL T16REC(IFT16 ,IPRINT,INEXTR)
                  CALL T16RCA(IFT16 ,IPRINT,NGCCPO,NGMTR ,IFGMTR,NVXSR ,
     >                        NMXSR ,B2CRI ,BRNIRR,NZONE ,RECXSV,RECXSM,
     >                        RECTMP,RECSCA,ZONVOL)
                ENDIF
                BURNUP(IBURN)=BRNIRR(1)
                WNKB(IBURN)=BRNIRR(2)
                EFJ=BRNIRR(3)
                CALL T16WDS(IPCPO ,NGCCPO,NVXSR ,NMXSR ,IBURN ,EFJ   ,
     >                      NAMDXS,ITYXS ,FLXINT,FLXDIS,OVERV ,RECXSV,
     >                      IDRXSM,RECXSM,RECSCA)
                INEXTR=INEXTR+1
              ENDDO
              CALL LCMPUT(IPCPO ,'BURNUP      ',NBURN,2,BURNUP)
              CALL LCMPUT(IPCPO ,'N/KB        ',NBURN,2,WNKB)
              CALL LCMSIX(IPCPO,NAMDIR,ILCMDN)
            ENDDO
          ENDIF
        ENDDO
      ENDDO
*----
*  RELEASE MEMORY
*----
      DEALLOCATE(ZONRAD,ZONVOL)
      DEALLOCATE(RECSCA,RECTMP,RECXSM,RECXSV,OVERV,FLXDIS,FLXINT,VQLE)
      DEALLOCATE(IDRXSM,MATMSH,KMSPEC)
      DEALLOCATE(WNKB,BURNUP)
      RETURN
      END
