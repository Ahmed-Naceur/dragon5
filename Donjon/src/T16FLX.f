*DECK T16FLX
      SUBROUTINE T16FLX(IFT16 ,IPRINT,NGCCPO,NG    ,NGMTR ,NMATZ ,
     >                  MMXM  ,MTRMSH,NZONE ,IFGMTR,VELMTR,IMIREG,
     >                  VOLUME,B2CRI ,FLXINT,FLXDIS,OVERV ,KMSPEC,
     >                  MATMSH,VQLE  ,ZONNUM, ZONRAD,ZONVOL)
*
*----
*
*Purpose:
*  Read main transport flux and compute integrated flux, 
*  flux disadvantage factor and 1/V cross sections.
*
*Author(s):
* G. Marleau
*
*Parameters: input
* IFT16   tape16 file unit.
* IPRINT  print level where:
*         =0 for no print; >=  1 print processing option.
* NGCCPO  number of edit groups.
* NG      number of groups on X-S library.
* NGMTR   number of main transport groups.
* NMATZ   number of mixtures.
* MMXM    maximum number of zones and main transport meshes.
* MTRMSH  number of main transport mesh points.
* NZONE   number of zones.
* IFGMTR  fewgroups for main transport.
* VELMTR  velocity for main transport.
* IMIREG  mixture update identifier where
*          =0 do not update;
*          =-1 update using CELLAV information;
*          > 0 update using specified region number.
*
*Parameters: output
* VOLUME  total volume.
* B2CRI   critical bucklings.
* FLXINT  volume integrated fluxes.
* FLXDIS  flux disadvantage factor.
* OVERV   1/V cross sections.
* KMSPEC  material types.
* MATMSH  material in each mesh.
* VQLE    volume of each mesh.
* ZONNUM  zone number.
* ZONRAD  zone radius.
* ZONVOL  zone volume.
*
*----
*
      IMPLICIT         NONE
      INTEGER          IFT16,IPRINT,NGCCPO,NG,
     >                 NGMTR,NMATZ,MMXM,MTRMSH,NZONE,IMIREG
      INTEGER          IFGMTR(NGCCPO),
     >                 KMSPEC(NMATZ),MATMSH(MMXM)
      REAL             VELMTR(NGMTR),VOLUME,B2CRI(3),
     >                 FLXINT(NGCCPO),FLXDIS(NGCCPO),
     >                 OVERV(NGCCPO),VQLE(MMXM)
      INTEGER          ZONNUM
      REAL             ZONVOL(NZONE), ZONRAD(NZONE)
*----
*  MEMORY ALLOCATION
*----
      REAL, ALLOCATABLE, DIMENSION(:,:) :: PHI
*----
*  T16 PARAMETERS
*----
      INTEGER          MAXKEY
      PARAMETER       (MAXKEY=3)
      CHARACTER        TKEY1(MAXKEY)*10,TKEY2(MAXKEY)*10,
     >                 RKEY1*10,RKEY2*10
      INTEGER          NKEY,IOPT,NBE,NID,IR
      REAL             RID
*----
*  LOCAL VARIABLES
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='T16FLX')
      INTEGER          IGR,IGC,IGD,IGF,IMIX,ITRFL,IBUCK
      REAL             B2INI(3)
      INTEGER          IZ
      REAL             CELLV
*----
*  SET END RECORDS FOR THIS SEARCH
*----
      ALLOCATE(PHI(NG,MMXM))
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6000) NAMSBR
      ENDIF
      IOPT=0
      CALL XDRSET(FLXINT,NGCCPO              ,0.0)
      CALL XDRSET(FLXDIS,NGCCPO              ,0.0)
      CALL XDRSET(OVERV ,NGCCPO              ,0.0)
*----
*  CELL VOLUME PER UNIT LENGTH 
*----
      REWIND(IFT16) 
      TKEY1(1)='REGION    '
      TKEY2(1)='DESCRIPTON'
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1  ,TKEY2  ,
     >            NBE   )
      IF(NBE .GT. 0) THEN
        READ(IFT16) RKEY1,RKEY2,NBE,ZONNUM, 
     >  (ZONRAD(IZ),IZ=1,NZONE),(ZONVOL(IZ), IZ=1, NZONE)
      ELSE
        CALL XABORT(NAMSBR//': KEYS '//TKEY1(1)//','//
     >              TKEY2(1)//' NOT FOUND ON TAPE16')
      ENDIF
      CELLV=0.0
      DO IZ=1, ZONNUM 
         CELLV=CELLV+ZONVOL(IZ)
      END DO 
*----
*  MTRFLX RECORDS
*----
      REWIND(IFT16)    
      NKEY=2
      TKEY1(2)='REGION    '
      TKEY2(2)='DESCRIPTON'
      TKEY1(1)='MTRFLX    '
      TKEY2(1)='FLUX      '
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >            NBE   )
      ITRFL=0
      IF(NBE .GT. 0 ) THEN
        ITRFL=1
      ELSE IF( NBE .LT. -1 ) THEN
        READ(IFT16) RKEY1,RKEY2,NBE
        IF(IMIREG .GT. 0) THEN
*----
*  Update mixture if IMIREG>0
*----
          TKEY1(2)='CELLAV    '
          TKEY2(2)='NGROUPS   '
          TKEY1(1)='REGION    '
          TKEY2(1)='FLUX      '
          DO IR=1,IMIREG-1
            CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >                  NBE   )
            IF( NBE .LE. 0 ) CALL XABORT(NAMSBR//
     >      ': REGION FLUX NOT AVAILABLE')
            READ(IFT16) RKEY1,RKEY2
          ENDDO
          CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >                NBE   )
          IF( NBE .GT. 0 ) ITRFL=2
        ELSE IF(IMIREG .LT. 0) THEN
*----
*  Update mixture using CELLAV information if IMIREG<0
*----
          TKEY1(2)='CELLAV    '
          TKEY2(2)='K         '
          TKEY1(1)='CELLAV    '
          TKEY2(1)='FLUX      '
          CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >                NBE   )
          IF( NBE .GT. 0 ) ITRFL=3
        ENDIF
      ENDIF
      IF( ITRFL .EQ. 0 ) THEN
         CALL XABORT(NAMSBR//
     >   ': CANNOT FIND '//TKEY1(1)//' '//TKEY2(1)//' OR '//
     >   TKEY1(2)//' '//TKEY2(2))
      ELSE IF(ITRFL .EQ. 1) THEN
*----
*  USE MTRFLX
*  1) CONDENSE AND HOMOGENIZE FLUX
*  2) COMPUTE FLUX DISADVANTAGE FACTOR
*  3) COMPUTE VOLUME
*  4) COMPUTE OVERV
*----
        READ(IFT16) RKEY1,RKEY2,NBE,NID,NID,
     >   (MATMSH(IR),VQLE(IR),
     >   (PHI(IGR,IR),IGR=1,NGMTR),IR=1,MTRMSH)
        VOLUME=0.0
        DO IR=1,MTRMSH
          IGF=0
          VOLUME=VOLUME+VQLE(IR)
          IF(IPRINT .GE. 100) THEN
            WRITE(IOUT,6100) IR,VQLE(IR)
            WRITE(IOUT,6110)(PHI(IGR,IR),IGR=1,NGMTR)
          ENDIF
          DO IGC=1,NGCCPO
            IGD=IGF+1
            IGF=IFGMTR(IGC)
            DO IGR=IGD,IGF
              FLXINT(IGC)=FLXINT(IGC)+PHI(IGR,IR)*VQLE(IR)
              OVERV(IGC)=OVERV(IGC)
     >                  +PHI(IGR,IR)*VQLE(IR)/VELMTR(IGR)
            ENDDO
          ENDDO
          IMIX=MATMSH(IR)
          IF(KMSPEC(IMIX) .EQ. 1) THEN
            IGF=0
            DO IGC=1,NGCCPO
              IGD=IGF+1
              IGF=IFGMTR(IGC)
              DO IGR=IGD,IGF
                FLXDIS(IGC)=FLXDIS(IGC)+PHI(IGR,IR)*VQLE(IR)
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ELSE IF(ITRFL .EQ. 2) THEN
*----
*  USE REGION FLUX
*  1) CONDENSE AND HOMOGENIZE FLUX
*  2) COMPUTE FLUX DISADVANTAGE FACTOR
*  4) COMPUTE OVERV
*----
        READ(IFT16) RKEY1,RKEY2,NBE,NID,NID,VOLUME,
     >    (PHI(IGR,1),IGR=1,NGMTR)
        IR=IMIREG
        IGF=0
        IF(IPRINT .GE. 100) THEN
          WRITE(IOUT,6100) IR,VOLUME
          WRITE(IOUT,6110)(PHI(IGR,IR),IGR=1,NGMTR)
        ENDIF
        DO 120 IGC=1,NGCCPO
          IGD=IGF+1
          IGF=IFGMTR(IGC)
          DO 121 IGR=IGD,IGF
            FLXINT(IGC)=FLXINT(IGC)+PHI(IGR,IR)*VOLUME
            OVERV(IGC)=OVERV(IGC)
     >              +(PHI(IGR,IR)*VOLUME)/VELMTR(IGR)
 121      CONTINUE
          FLXDIS(IGC)=FLXINT(IGC)
 120    CONTINUE
      ELSE
*----
*  USE CELLAV FLUX
*  1) CONDENSE AND HOMOGENIZE FLUX
*  2) COMPUTE FLUX DISADVANTAGE FACTOR
*  3) COMPUTE VOLUME
*  4) COMPUTE OVERV
*----
        IR=1
        VOLUME=1.0
        READ(IFT16) RKEY1,RKEY2,NBE,
     >   (PHI(IGR,IR),IGR=1,NGMTR)
        IF(IPRINT .GE. 100) THEN
          WRITE(IOUT,6101)
          WRITE(IOUT,6110)(PHI(IGR,IR),IGR=1,NGMTR)
        ENDIF
        DO IGR=1, NGMTR 
          PHI(IGR,IR)=PHI(IGR,IR)/CELLV
        ENDDO 
        IGF=0
        DO 130 IGC=1,NGCCPO
          IGD=IGF+1
          IGF=IFGMTR(IGC)
          DO 131 IGR=IGD,IGF
            FLXINT(IGC)=FLXINT(IGC)+PHI(IGR,IR)
            OVERV(IGC)=OVERV(IGC)+PHI(IGR,IR)/VELMTR(IGR)
 131      CONTINUE
          FLXDIS(IGC)=FLXINT(IGC)
 130    CONTINUE
      ENDIF
      DO 140 IGC=1,NGCCPO
        FLXDIS(IGC)=FLXDIS(IGC)/FLXINT(IGC)
        OVERV(IGC)=OVERV(IGC)/FLXINT(IGC)
 140  CONTINUE
*----
*  RADIAL AND AXIAL DIFFUSION COEFFICIENTS
*  AND BUCKLING
*----
      TKEY1(2)='CELLAV    '
      TKEY2(2)='K         '
      TKEY1(1)='CELLAV    '
      TKEY2(1)='DIFFUSION '
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >            NBE   )
      IF( NBE .NE. 5*NGMTR+5 ) CALL XABORT(NAMSBR//
     >': CANNOT FIND '//TKEY1(1)//' '//TKEY2(1))
      READ(IFT16) RKEY1,RKEY2,NBE,(NID,IR=1,3),
     >           (RID,IGR=1,NGMTR),
     >           (RID,IGR=1,NGMTR),
     >           (RID,IGR=1,NGMTR),
     >           (B2INI(IR),IR=1,2)
      TKEY1(1)='CELLAV    '
      TKEY2(1)='CRITICALB '
      CALL T16FND(IFT16 ,IPRINT,IOPT  ,NKEY  ,TKEY1 ,TKEY2 ,
     >            NBE   )
      IF( NBE .EQ. 2*NGMTR+4 ) THEN
        READ(IFT16) RKEY1,RKEY2,NBE,IBUCK,
     >             (B2CRI(IR),IR=1,3)
        IF(IBUCK .EQ. 2) THEN
          B2CRI(3)=B2INI(1)+B2CRI(2)
          B2CRI(1)=B2INI(1)/B2CRI(3)
          B2CRI(2)=B2CRI(2)/B2CRI(3)
        ELSE IF(IBUCK .EQ. 3) THEN
          B2CRI(3)=B2CRI(1)+B2INI(2)
          B2CRI(1)=B2CRI(1)/B2CRI(3)
          B2CRI(2)=B2INI(2)/B2CRI(3)
        ELSE
          B2CRI(1)=B2CRI(1)/B2CRI(3)
          B2CRI(2)=B2CRI(2)/B2CRI(3)
        ENDIF
      ELSE IF(NBE .EQ. -2) THEN
        B2CRI(3)=B2INI(1)+B2INI(2)
        B2CRI(1)=B2INI(1)/B2CRI(3)
        B2CRI(2)=B2INI(2)/B2CRI(3)
      ELSE
        CALL XABORT(NAMSBR//
     >': CANNOT FIND '//TKEY1(2)//' '//TKEY2(2))
      ENDIF
      IF(IPRINT .GE. 10) THEN
        WRITE(IOUT,6010) (FLXINT(IGC),IGC=1,NGCCPO)
        WRITE(IOUT,6011) (FLXDIS(IGC),IGC=1,NGCCPO)
        WRITE(IOUT,6012) (OVERV(IGC),IGC=1,NGCCPO)
        WRITE(IOUT,6013) (B2CRI(IR),IR=1,3)
        WRITE(IOUT,6001)
      ENDIF
      RETURN
*----
*  PRINT FORMAT
*----
 6000 FORMAT(1X,5('*'),' OUTPUT FROM ',A6,1X,5('*'))
 6001 FORMAT(1X,30('*'))
 6010 FORMAT(6X,'INTEGRATED FLUXES'/
     >1P,10(2X,E10.3))
 6011 FORMAT(6X,'FLUX DISADVANTAGE FACTORS'/
     >1P,10(2X,E10.3))
 6012 FORMAT(6X,'1/V '/
     >1P,10(2X,E10.3))
 6013 FORMAT(6X,'CRITICAL BUCKLINGS'/
     >1P,3(2X,E10.3))
 6100 FORMAT(6X,'MAIN TRANSPORT GROUP FLUX IN REGION  = ',I10,
     >       5X,'OF VOLUME = ',1P,E10.3)
 6101 FORMAT(6X,'CELLAV MAIN TRANSPORT GROUP FLUX  ')
 6110 FORMAT(1P,10(2X,E10.3))
      END
